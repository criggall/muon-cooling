//	BLParam.cc
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

#include "BLParam.hh"
#include "BLEvaluator.hh"

std::map<G4String,G4String> *BLParam::paramMap;
std::map<G4String,G4String> *BLParam::paramHelpText;

BLParam Param;	// the global Param object

void BLParam::init() {
	if (!paramMap) paramMap = new std::map<G4String,G4String>;
	if (!paramHelpText) paramHelpText = new std::map<G4String,G4String>;
}

void BLParam::printParam()
{
	init();
	printf("\nPARAMETERS:\n");
	std::map<G4String,G4String>::iterator i;
	for(i=(*paramMap).begin(); i!=(*paramMap).end(); ++i)
		printf("%15s=%s\n",i->first.c_str(),i->second.c_str());
}

G4String BLParam::getString(G4String name)
{
	init();
	if((*paramMap).count(name) == 0) {
		// define it from the environment, if possible
		char *p = getenv(name.c_str());
		if(p) {
			setParam(name,p);
		} else {
			BLCommand::printError("ERROR: Unknown parameter '%s'",
						name.c_str());
		}
	}
	return (*paramMap)[name];	// "" if not found
}

G4double BLParam::getDouble(G4String name)
{
	G4String v = getString(name);
	if(v.size() == 0) return -HUGE_VAL;
	BLEvaluator e;
	G4double val = e.evaluate(v);
	if(e.status() != HepTool::Evaluator::OK) return -HUGE_VAL;
	return val;
}

G4int BLParam::getInt(G4String name)
{
	G4double v = getDouble(name);
	G4int i = (G4int)(v+0.5);
	if(v < 0) --i;
	if(v == -HUGE_VAL || fabs(v-i) > 0.001)
		BLCommand::printError("ERROR non-integer value %.4g used as int in '%s'",v,name.c_str());
	return i;
}

void BLParam::setParam(G4String name, G4String value)
{
	init();
	(*paramMap)[name] = value;
}

void BLParam::setParam(G4String name, G4double value)
{
	char tmp[32];
	sprintf(tmp,"%g",value);
	setParam(name,tmp);
}

void BLParam::setParam(G4String name, G4int value)
{
	char tmp[32];
	sprintf(tmp,"%d",value);
	setParam(name,tmp);
}

bool BLParam::isDefined(G4String name)
{
	if(paramMap)
		return paramMap->count(name) != 0;
	return false;
}

G4String BLParam::expand(G4String str)
{
	static G4String nameChars("ABCDEFGHIJKLMNOPQRSTUVWXYZ_"
				  "abcdefghijklmnopqrstuvwxyz0123456789");
	G4String out;

	G4String::size_type place=0;
	while(place < str.size()) {
		G4String::size_type i=str.find('$',place);
		if(i == str.npos) {
			out += str.substr(place);
			break;
		}
		out += str.substr(place,i-place);
		place = i + 1;
		// leave $1 - $9 and $# alone (for define command)
		if(isdigit(str[place]) || str[place] == '#') {
			out += "$";
			continue;
		}
		// replace $$ with $ (delayed evaluation of $param)
		if(str[place] == '$') {
			out += "$";
			++place;
			continue;
		}
		G4String::size_type j=str.find_first_not_of(nameChars,place);
		if(j == str.npos) j = str.size();
		G4String name=str.substr(place,j-place);
		if(j == place || isdigit(str[place]))
		    BLCommand::printError("ERROR: Invalid parameter name '%s'",
		    					name.c_str());
		else
		    out += getString(name);
		place = j;
	}

	return out;
}

void BLParam::setHelpText(G4String name, G4String text)
{
	init();
	(*paramHelpText)[name] = text;
}

G4String BLParam::getHelpText()
{
	G4String s;
	std::map<G4String,G4String>::iterator i;
	for(i=(*paramHelpText).begin(); i!=(*paramHelpText).end(); ++i) {
		G4String::size_type n = s.size();
		s += "    ";
		s += i->first;
		s += " ";
		while(s.size()-n < 24) s += " ";
		s += i->second;
		s += "\n";
	}

	s +=
"\nviewer can be just a name, in which case the Geant4 /run/beamOn command \n"
"must be used to generate events to be displayed. viewer can also have two \n"
"or three parts: 'name,1,1' where the first 1 is the number of events per \n"
"image and the second 1 is the number of images to generate. For the 'best'\n"
"viewer, which is Open Inventor, 999999 images are generated, and the user\n"
"does File/Escape to advance from one image to the next.\n\n";

	return s;
}
