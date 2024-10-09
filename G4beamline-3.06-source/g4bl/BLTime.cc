//	BLTime.cc

#include "BLTime.hh"

#ifdef __APPLE__
#include <unistd.h>
#include <sys/time.h>
#endif
#ifdef __linux
#include <unistd.h>
#include <time.h>
#include <sys/time.h>
#endif
#ifdef WIN32
#include <windows.h>
#include <time.h>
#undef max
#undef min
#endif

long BLTime::time() { return ::time(0); }

long BLTime::timems() {
#ifndef WIN32
	struct timeval tv;
	gettimeofday(&tv,0);
	return tv.tv_sec*1000+tv.tv_usec/1000;
#else
/***
	struct _timeb timebuffer;
	_ftime( &timebuffer );
	return timebuffer.time*1000 + timebuffer.millitm;
***/
	return ::time(0)*1000;
#endif
}

long long BLTime::timeus() {
#ifndef WIN32
	struct timeval tv;
	gettimeofday(&tv,0);
	return tv.tv_sec*1000000LL+tv.tv_usec;
#else
/***
	struct _timeb timebuffer;
	_ftime( &timebuffer );
	return timebuffer.time*1000000LL + timebuffer.millitm*1000LL;
***/
	return ::time(0)*1000000LL;
#endif
}

	/// sleepms() will sleep for a specified number of milliseconds.
void BLTime::sleepms(int ms) {
#ifdef WIN32
	Sleep(ms);
#else
	usleep(ms*1000);
#endif
}
