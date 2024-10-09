//	BLWindowShape.cc
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

#include <stdio.h>
#include <fstream>

#include "BLWindowShape.hh"
#include "BLCommand.hh"

BLWindowShape::BLWindowShape(G4String filename) : r(), z(), t()
{
	flangeInnerRadius = 0.0;
	flangeOuterRadius = 0.0;
	flangeInsideZ = 0.0;
	flangeOutsideZ = 0.0;

	G4String file = filename;

	std::ifstream in;
	if(filename != "")
		in.open(file.c_str());
	if(!in.good()) {
		BLCommand::printError("BLWindowShape Cannot open file '%s'",
				file.c_str());
		return;
	}

	bool flangeLine = true;
	int lineno = 0;
	while(in.good()) {
		char line[1024+1];
		line[0] = '\0';
		in.getline(line,sizeof(line));
		++lineno;
		if(line[0] == '*')
			printf("%s\n",line);
		if(line[0] == '\0' || line[0] == '#' || line[0] == '*')
			continue;
		int n = 0;
		if(flangeLine) {
			flangeLine = false;
			n = sscanf(line,"%lf %lf %lf %lf",&flangeInnerRadius,
			    &flangeOuterRadius,&flangeInsideZ,&flangeOutsideZ);
			n -= 1;	// 4 items here, check is for 3
		} else {
			G4double rr=0.0,zz=0.0,tt=0.0;
			n = sscanf(line,"%lf %lf %lf",&rr,&zz,&tt);
			r.push_back(rr);
			z.push_back(zz);
			t.push_back(tt);
		}
		if(n != 3) {
			BLCommand::printError("BLWindowShape Invalid syntax"
					" '%s' line %d\n",file.c_str(),lineno);
			break;
		}
	}
	in.close();
}
