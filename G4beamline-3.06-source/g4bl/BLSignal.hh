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

#ifndef BLSIGNAL_HH
#define BLSIGNAL_HH

/**	class BLSignal implements a signal handler to permit G4beamline
 *	to exit cleanly when a UNIX signal is received. The first signal
 *	received merely sets an intrnal flag which BLManager can check;
 *	a second signal will exit the program immediately.
 *
 *	SIGUSR1 and SIGUSR2 are special, and do not exit the program;
 *	they are merely reported via functions sigusr1() and sigusr2().
 *
 *	Note that a second received signal causes an immediate call to
 *	g4bl_exit() from the signal handler.
 *
 *	If BLAlarm is also used, BLSignal::init() must be called first.
 *	In this case, the SIGALRM signal is not handled by this class.
 *
 *	Note that BLManager checks received() in UserSteppingAction(), and
 *	issues a fatal exception when any signal is received.
 **/
class BLSignal {
public:
	/// init() will setup the signal handler
	static void init();

	/// received() returns true if a signal has been received.
	/// Does not include SIGUSR1 or SIGUSR2. Does not include SIGALRM
	/// if BLAlarm is used.
	static bool received() { return signalReceived; };

	/// value() returns the signal # of the last sginal received.
	/// Does not include SIGUSR1 or SIGUSR2. Does not include SIGALRM
	/// if BLAlarm is used.
	static int value() { return prev; }

	/// clear() will clear the signalReceived flag.
	void clear() { signalReceived = false; prev = 0; }

	/// setSignalReceived() will set the signalReceived flag.
	/// (Used surrounding console input so one ^C will exit immediately.)
	static void setSignalReceived(bool v) { signalReceived = v; }

	/// sigusr1() returns true if one or more SIGUSR1 signals was
	/// received since the last call;
	static bool sigusr1() { bool v=usr1; usr1=false; return v; }

	/// sigusr2() returns true if one or more SIGUSR2 signals was
	/// received since the last call;
	static bool sigusr2() { bool v=usr2; usr2=false; return v; }

	/// generateStackTrace() will generate a stack trace from its caller.
	static void generateStackTrace();

	/// writeStackTrace() writes sighandler()'s stack trace to a fd.
	/// If no signal has been received, nothing is written.
	static void writeStackTrace(int fd);

private:
	static bool signalReceived;
	static bool usr1, usr2;
	static int prev;
	static void *stackTrace[256];
	static int nStackTrace;
	static int sigStackTrace;

	/// sighandler() handles the alarm signal
	static void sighandler(int sig);
};

#endif // BLSIGNAL_HH
