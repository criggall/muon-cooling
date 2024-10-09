//	BLAlarm.hh
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

#ifndef BLALARM_HH
#define BLALARM_HH

/**	class BLAlarm implements an alarm clock.
 *	used to prevent an infinite tracking loop from hanging an entire job.
 *
 *	If BLSignal is also used, BLSignal::init() must be called first.
 *
 *	This class is 100% static.
 **/
class BLAlarm {
	static void (*customHandler)();
public:
	/// clear will clear any alarm. 
	static void clear();

	/// set() sets an alarm in the future.
	/// Implicitly clears any previous alarm.
	/// seconds <= 0 does a clear().
	static void set(int seconds);

	/// init() will setup the signal handler.
	/// Must be called before set().
	static void init();

	/// setCustomhandler() sets a custom alarm handler.
	/// returns the previous handler, or 0 if none.
	static void (*setCustomHandler(void (*handler)()))();

	/// timeRemaining() returns the time remaining, in seconds.
	/// returns -1 if no alarm is set.
	int timeRemaining();

private:
	/// sighandler() handles the alarm signal
	static void sighandler(int sig);
};

#endif // BLALARM_HH
