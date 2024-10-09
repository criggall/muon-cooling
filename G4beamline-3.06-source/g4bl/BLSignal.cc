//	BLSignal.cc
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

#include <string.h>
#include "BLSignal.hh"

#ifdef WIN32
#include <Windows.h>
#include <signal.h>
#include <stdio.h>
#include "WinStackWalker.hh"
#else
#include <stdlib.h>
#include <unistd.h>
#include <signal.h>
#include <stdio.h>
#include <execinfo.h>
#include <dlfcn.h>
#include <cxxabi.h>
#endif

extern void g4bl_exit(int value);

bool BLSignal::signalReceived = false;
bool BLSignal::usr1 = false;
bool BLSignal::usr2 = false;
int BLSignal::prev = 0;
void *BLSignal::stackTrace[256] = { 0 };
int BLSignal::nStackTrace = 0;
int BLSignal::sigStackTrace = 0;

#ifdef WIN32
static BOOL CtrlHandler( DWORD fdwCtrlType )
{
	BLSignal::setSignalReceived(true);
	return TRUE;
}
#endif

void BLSignal::init() 
{
#ifdef WIN32
	static int sigs[] = {SIGINT, SIGILL, SIGABRT,
		SIGSEGV, SIGTERM, SIGFPE};
#else
	static int sigs[] = {SIGHUP, SIGINT, SIGQUIT, SIGILL, SIGABRT,
		SIGKILL, SIGBUS, SIGSEGV, SIGSYS, SIGPIPE, SIGTERM, SIGXCPU,
		SIGUSR1, SIGUSR2};
#endif

	signalReceived = false;
	prev = 0;

	for(unsigned i=0; i<sizeof(sigs)/sizeof(sigs[0]); ++i)
		signal(sigs[i],sighandler); 

#ifdef WIN32
	SetConsoleCtrlHandler( (PHANDLER_ROUTINE) CtrlHandler, TRUE );
#endif
}

void BLSignal::sighandler(int sig) 
{
	signal(sig,sighandler);

#ifndef WIN32
	if(sig == SIGUSR1) {
		usr1 = true;
		return;
	}
	if(sig == SIGUSR2) {
		usr2 = true;
		return;
	}

	// If this is the first signal, get stackTrace
	if(nStackTrace == 0) {
		sigStackTrace = sig;
		nStackTrace = backtrace(stackTrace,
				sizeof(stackTrace)/sizeof(stackTrace[0]));
	}
#else
	// WinStackWalker with an additional output to the console:
	class MyStackWalker : public WinStackWalker {
	public:
  		MyStackWalker() : WinStackWalker() {}
  		virtual void OnOutput(LPCSTR szText) {
			fprintf(stderr,szText);
			WinStackWalker::OnOutput(szText); // to IDE debug window
		}
	};

	fflush(stdout);
	fprintf(stderr, "\n"
		"********************************************\n"
		"***  Signal received (%d), stack trace:  ***\n"
		"********************************************\n",
	sig);
 
	MyStackWalker sw;
	sw.ShowCallstack();

	fprintf(stderr,"g4beamline is exiting after a signal\n");

	g4bl_exit(99);
#endif

	if(signalReceived) {
		static bool moreThanTwo=false;
		if(moreThanTwo) 
			abort();
		moreThanTwo = true;
		fflush(stdout);
		if(prev == 0 && sig == 2) { // ^C during input -- suppress print
			nStackTrace = 0;
			fprintf(stderr,"\n");
		} else {
		    fprintf(stderr, "\n"
		    "***************************************************\n"
		    "***  Multiple signals received (%d,%d), exiting  ***\n"
		    "***************************************************\n",
		    prev,sig);
		    writeStackTrace(2);
		}
		g4bl_exit(99);
	}

	prev = sig;
	signalReceived = true;
}

void BLSignal::generateStackTrace()
{
#ifndef WIN32
	sigStackTrace = -1;
	nStackTrace = backtrace(stackTrace,
				sizeof(stackTrace)/sizeof(stackTrace[0]));
#else
	sighandler(13579);
#endif
}

void BLSignal::writeStackTrace(int fd)
{
	fflush(stdout);
	fflush(stderr);
#ifndef WIN32
	if(nStackTrace > 0) {
		int n = nStackTrace;
		nStackTrace = 0;

		fprintf(stderr,"\nStack Trace for signal %d:\n",sigStackTrace);
		for(int i=0; i<n; ++i) {
			const char *fname = "???";
			const char *symbol = "???";
			long offset = 0L;
			Dl_info dlinfo;
			if(dladdr(stackTrace[i],&dlinfo)) {
				offset = (char *)stackTrace[i] - 
						(char *)dlinfo.dli_saddr;
				fname = strrchr(dlinfo.dli_fname,'/');
				if(!fname) fname = strrchr(dlinfo.dli_fname,'\\');
				if(fname) ++fname; // omit leading '/' or '\\'
				if(!fname) fname = dlinfo.dli_fname;
				symbol = dlinfo.dli_sname;
				int status;
				char *demangled = 
				    abi::__cxa_demangle(symbol,NULL,0,&status);
				if(status == 0 && demangled != 0)
					symbol = demangled;
				// Deliberate memory leak: don't free demangled
				// (so symbol remains valid).
			}
			fprintf(stderr,
			    (sizeof(void*)==8 ? "%-2d  %18p  %s  %s + %ld\n"
					      : "%-2d  %10p  %s  %s + %ld\n"),
					i,stackTrace[i],fname,symbol,offset);

		}
		fprintf(stderr,"\n");
	}
#endif // !WIN32
}
