//	BLMinimize.cc

#include "BLMinimize.hh"

#ifdef USE_TMINUIT
BLMinimize *BLMinimize::current = 0;
#else
int BLMinimize_dummy = 0;
#endif
