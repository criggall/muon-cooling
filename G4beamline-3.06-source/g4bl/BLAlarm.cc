//	BLAlarm.cc
/*
This source file is part of G4beamline, http://g4beamline.muonsinc.com
Copyright (C) 2002-2013 by Tom Roberts, all rights reserved.

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

http://www.gnu.org/copyleft/gpl.html
*/

#include "BLAlarm.hh"
#include "BLRunManager.hh"
#include "BLTime.hh"

void (*BLAlarm::customHandler)() = 0;
const unsigned long NEVER=0xFFFFFFFFL;
static volatile unsigned long alarm_time=NEVER;

#ifdef WIN32

/**	On Windows we need to implement the alarm generation in a Thread.
 *
 *	SIGALRM does not exist, use SIGABRT.
 **/

#include "Windows.h"
#include "process.h"
#include "signal.h"

#include "globals.hh"

static void (*prev_handler)(int)=0;

static void thread(void *param)
{
	for(;;) {
		unsigned long now=BLTime::time();
		if(alarm_time > now+1) {
			Sleep(1000); // ms
		} else if(alarm_time > now) {
			Sleep(50); // ms
		} else {
			// Cannot call G4Exception here, wrong thread
			raise(SIGABRT);
			Sleep(1000); // ms
		}
	}
}

void BLAlarm::clear()
{
	set(0);
}

void BLAlarm::set(int seconds)
{
	if(seconds <= 0)
		alarm_time = NEVER;
	else
		alarm_time = BLTime::time() + seconds;
}

void BLAlarm::init()
{
	prev_handler = signal(SIGABRT,sighandler); 
	_beginthread(thread,0,NULL);
}

void BLAlarm::sighandler(int sig)
{ 
	if(alarm_time > BLTime::time()) {
		(*prev_handler)(sig);
		return;
	}
	signal(sig,sighandler); 

	if(customHandler != 0) {
		(*customHandler)();
		return;
	}

	BLRunManager::getObject()->abandonCurrentEvent();
	G4Exception("BLAlarm","Alarm Signal",FatalException,
			"SIGALRM fired, cannot be recovered");
}

void (*BLAlarm::setCustomHandler(void (*handler)()))()
{
	void (*t)() = customHandler;
	customHandler = handler;
	return t;
} 

int BLAlarm::timeRemaining()
{
	unsigned long t = BLTime::time();
	if(t > alarm_time || alarm_time == NEVER)
		return -1;
	return alarm_time - t;
}


#else // WIN32

/**	Linux and Mac OS X handle alarm signals natively.
 **/

#include <stdio.h>
#include <unistd.h>
#include <signal.h>
#include <time.h>

#include "globals.hh"

void BLAlarm::clear()
{
	set(0);
}

void BLAlarm::set(int seconds)
{
	alarm(seconds < 0 ? 0 : seconds);
	if(seconds <= 0)
		alarm_time = NEVER;
	else
		alarm_time = BLTime::time() + seconds;
}

void BLAlarm::init()
{
	signal(SIGALRM,sighandler); 
}

void BLAlarm::sighandler(int sig)
{ 
	alarm(0);
	alarm_time = NEVER;
	signal(sig,sighandler); 

	if(customHandler != 0) {
		(*customHandler)();
		return;
	}

	BLRunManager::getObject()->abandonCurrentEvent();
	G4Exception("BLAlarm","Alarm Signal",FatalException,
			"SIGALRM fired, cannot be recovered");
}

void (*BLAlarm::setCustomHandler(void (*handler)()))()
{
	void (*t)() = customHandler;
	customHandler = handler;
	return t;
} 

int BLAlarm::timeRemaining()
{
	unsigned long t = BLTime::time();
	if(t > alarm_time || alarm_time == NEVER)
		return -1;
	return alarm_time - t;

}

#endif // ! WIN32
