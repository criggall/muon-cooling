//	BLRealTimeMonitor.hh
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

#ifndef BLREALTIMEMONITOR_HH
#define BLREALTIMEMONITOR_HH

#include "BLTime.hh"

/**	class BLRealTimeMonitor creates a versatile multi-state monitor
 *	of real time usage, with resolution of microseconds (or the
 *	system clock resolution).
 *
 *	This entire class is inline.
 **/
class BLRealTimeMonitor {
	long long total[8];
	int state;
public:
	/// constructor and destructor.
	BLRealTimeMonitor() 
		{ for(int i=0; i<8; ++i) total[i]=0LL; state=-1; }
	~BLRealTimeMonitor() { }

	/// setState() stops accumulating in the previous state, and starts
	/// accumulating in the new state.
	/// state can range from 0 to 7; -1 means no accumulation.
	/// initial state is -1.
	void setState(int newState) {
		long long now = BLTime::timeus();
		if(state >= 0 && state < 8) total[state] += now;
		state = newState;
		if(state >= 0 && state < 8) total[state] -= now;
	}

	/// incrState() adds the increment to the current state.
	void incrState(int incr=1) { setState(state+incr); }

	/// getTime() returns the total elapsed time (in seconds) while in the
	/// specified state.
	float getTime(int _state) {
		if(_state < 0 || _state >= 8) return 0.0;
		long long v=total[_state];
		if(state == _state) v += BLTime::timeus();
		return (float)v/1.0E6;
	}

	/// estimateResolution() will estimate the resolution of the
	/// real-time monitor.
	static float estimateResolution() {
		float v=0.0;
		for(int i=0; i<100; ++i) {
			v += estimateResolution1();
			if(v > 0.1)
				return v/(i+1);
		}
		return v/100.0;
	}

	/// estimateResolution1() will make a single estimate of the
	/// resolution of the real-time monitor.
	static float estimateResolution1() {
		long long now=0LL, prev=BLTime::timeus();
		do {
			now = BLTime::timeus();
		} while(now == prev);
		prev = now;
		do {
			now = BLTime::timeus();
		} while(now == prev);
		return (float)(now-prev)/1.0E6;
	}
};

#endif // BLREALTIMEMONITOR_HH
