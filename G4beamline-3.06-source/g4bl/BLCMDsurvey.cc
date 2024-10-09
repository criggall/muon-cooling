//	BLCMDsurvey.cc
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

#include <map>
#include <stdio.h>

#include "BLManager.hh"
#include "BLCommand.hh"
#include "BLElement.hh"
#include "BLCallback.hh"
#include "BLGroup.hh"
#include "BLCoordinates.hh"
#include "BLWriteAsciiFile.hh"

class BLCMDsurvey : public BLCommand, public BLCallback {
	struct Entry {
		G4String command;
		G4String name;
		G4ThreeVector front;
		G4ThreeVector rear;
	};
	G4String coordinates;
	G4String filename;
	std::multimap<G4double,Entry> list;
	FILE *out;
	BLCoordinateType coordType;
public:
	BLCMDsurvey();

	G4String commandName() { return "survey"; }
	
	int command(BLArgumentVector& argv, BLArgumentMap& namedArgs);

	void defineNamedArgs();

	void callback(int type);

	void child(const BLChild &c, const G4RotationMatrix &parentRot,
					const G4ThreeVector &parentPos);
};
BLCMDsurvey defaultSurveyCommand;

BLCMDsurvey::BLCMDsurvey()
{
	registerCommand(BLCMDTYPE_OTHER);
	setSynopsis("survey command.");
	setDescription(
	"The survey command writes the coordinates of the front and rear points "
	"of each element to a file (stdout if '-'). For most elements the points "
	"written are the geometrical centers of their front and rear faces. "
	"Elements are all objects put into the world via the place command.\n\n"
	"Entries are sorted by Zfront. The format of each line is:\n"
	"  command name Xfront Yfront Zfront Xrear Yrear Zrear\n\n"
	"Coordinates can be global, centerline, or reference (if present).");
	coordinates = "global";
	filename = "-";
	out = stdout;
}

int BLCMDsurvey::command(BLArgumentVector& argv, BLArgumentMap& namedArgs)
{
	BLCMDsurvey *s = new BLCMDsurvey(defaultSurveyCommand);
	s->handleNamedArgs(namedArgs);

	BLManager::getObject()->registerCallback(s,1);

	printf("survey\n");

	return 0;
}

void BLCMDsurvey::defineNamedArgs()
{
	argString(coordinates,"coordinates","Coordinate type (global).");
	argString(filename,"filename","Filename to write (-).");
	argString(filename,"file","Synonym for filename.");
	argString(coordinates,"coord","Synonym for coordinates.");
}

void BLCMDsurvey::callback(int type)
{
	if(type != 1) return;

	out = BLWriteAsciiFile::fopen(filename.c_str());
	if(!out)
		G4Exception("survey","Cannot write file",
				FatalException,filename);
	fprintf(out,"#survey using %s coordinates\n",coordinates.c_str());
	fprintf(out,"#command name Xfront Yfront Zfront Xrear Yrear Zrear\n");
	fprintf(out,"#- - mm mm mm mm mm mm\n");

	list.clear();

	coordType = BLCoordinates::getCoordinateType(coordinates);

	const std::vector<BLChild> &worldChildren = 
					BLGroup::getWorld()->getChildren();
	std::vector<BLChild>::const_iterator it;
	for(it=worldChildren.begin(); it!=worldChildren.end(); ++it)
		child(*it,G4RotationMatrix(),G4ThreeVector());

	std::multimap<G4double,Entry>::iterator jt;
	for(jt=list.begin(); jt!=list.end(); ++jt) {
		Entry& e = jt->second;
		const char *q = "";
		if(e.name.find(" ") != e.name.npos)
			q = "\"";
		if(e.name.find("\"") != e.name.npos)
			q = "'";
		fprintf(out,"%s %s%s%s    %.3f %.3f %.3f    %.3f %.3f %.3f\n",
			e.command.c_str(),q,e.name.c_str(),q,
			e.front.x(),e.front.y(),e.front.z(),
			e.rear.x(),e.rear.y(),e.rear.z());
	}
}

void BLCMDsurvey::child(const BLChild &c, const G4RotationMatrix &parentRot,
						const G4ThreeVector &parentPos)
{
	G4RotationMatrix rot = parentRot * c.rot;
	G4ThreeVector pos = parentRot * c.offset + parentPos;
	G4ThreeVector front = c.element->getSurveyPoint(0);
	front = rot * front + pos;
	G4ThreeVector rear = c.element->getSurveyPoint(1);
	rear = rot * rear + pos;

	Entry e;
	BLCoordinates coord;
	coord.setGlobal(front,0.0);
	coord.getCoords(coordType,e.front);
	coord.setGlobal(rear,0.0);
	coord.getCoords(coordType,e.rear);
	e.command = c.element->commandName();
	e.name = c.name;
	list.insert(std::pair<G4double,Entry>(e.front[2],e));

	const std::vector<BLChild> &children = c.element->getChildren();
	std::vector<BLChild>::const_iterator it;
	for(it=children.begin(); it!=children.end(); ++it)
		child(*it,rot,pos);
}
