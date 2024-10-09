//	BLTime.hh - system-independent time from epoch

#ifndef BLTIME_HH
#define BLTIME_HH

/**	BLTime - get time since epoch.
 *
 *	Entire class could be inline, but is not to prevent bad optimization.
 *
 *	NOTE: On Windows the resolution of timems() and timeus() is only
 *	a second! On Mac OS X and on Linux they have the corresponding
 *	resolutions.
 **/
class BLTime {
public:
	/// time() returns time in seconds since epoch (1970)
	static long time();

	/// timems() returns the time in milliseconds since epoch (1970)
	static long timems();

	/// timeus() returns the time in microseconds since epoch (1970)
	static long long timeus();

	/// sleepms() will sleep for a specified number of milliseconds.
	static void sleepms(int ms);
};

#endif // BLTIME_HH
