//	BLAssert.hh
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

#ifndef BLASSERT_H
#define BLASSERT_H

/// Macro BLAssert() replaces the assert() macro with a G4Exception.
#define BLAssert(EXPR) \
((EXPR) ? (void)0 : (G4Exception(__FILE__,"Assert Failure",FatalException,#EXPR), (void)0) )

///#define BLAssert(E) assert(E)

#endif // BLASSERT_H
