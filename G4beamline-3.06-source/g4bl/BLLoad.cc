//	BLLoad.cc
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

#ifdef USE_SHARED_OBJECTS

#include "BLLoad.hh"

std::vector<BLLoad*> BLLoad::list;

BLLoad::BLLoad(G4String name, void *p)
{
	libname=name;
	handle=p;
	list.push_back(this);
}

#else
int BLLoad_dummy;
#endif // USE_SHARED_OBJECTS
