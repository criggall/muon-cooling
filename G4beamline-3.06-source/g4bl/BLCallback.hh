//	BLCallback.hh

#ifndef BLCALLBACK_HH
#define BLCALLBACK_HH

/**	class BLCallback is for callbacks from the BLManager.
 *
 **/
class BLCallback {
public:
	BLCallback() { }
	virtual ~BLCallback() { }

	/// callback() implements the user action.
	/// Copied from BLManager.hh:
	/// type=0 for pre-Tune particle,
	/// type=1 for post-Reference (pre-beam tracking),
	/// type=2 for post-beam tracking.
	/// type=3 for replacing the main program loop.
	/// type=4 for visualization.
	/// type=-1 for inside physics definition
	virtual void callback(int type) { }
};

#endif // BLCALLBACK_HH
