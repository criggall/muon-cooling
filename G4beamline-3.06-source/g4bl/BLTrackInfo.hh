//	BLTrackInfo.hh
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

#ifndef BLTRACKINFO
#define BLTRACKINFO

#include "BLCoordinates.hh"

/**	class BLTrackInfo carries information about the track.
 *	It is derived from BLCoordinates so all the code that assumes it
 *	is a BLCoordinates will still work without change.
 **/
class BLTrackInfo : public BLCoordinates {
	int externalTrackID;
	int externalParentID;
	static G4Allocator<BLTrackInfo> alloc;
public:
	/// Constructor.
	BLTrackInfo() : BLCoordinates() { externalTrackID=externalParentID=0; }

	/// Destructor
	virtual ~BLTrackInfo() { }

	/// operator new - using G4Allocator.
	/// NOTE: any derived class MUST re-implement this.
	void *operator new(size_t) { return (void *)alloc.MallocSingle(); }

	/// operator delete - using G4Allocator.
	/// NOTE: any derived class MUST re-implement this.
	void operator delete(void *p) { alloc.FreeSingle((BLTrackInfo *)p); }

	/// get external trackID
	int getExternalTrackID() const { return externalTrackID; }

	/// set external trackID
	void setExternalTrackID(int v) { externalTrackID = v; }

	/// get external ParentID
	int getExternalParentID() const { return externalParentID; }

	/// set external ParentID
	void setExternalParentID(int v) { externalParentID = v; }
};

#endif // BLTRACKINFO
