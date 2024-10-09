//	BLCMDplace.cc
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

#include "CLHEP/Geometry/Transform3D.h"

#include "BLAssert.hh"
#include "BLCommand.hh"
#include "BLElement.hh"
#include "BLGroup.hh"
#include "BLCoordinates.hh"

const double UNDETERMINED = -1.0e10;

/**	class BLCMDplace implements the place command.
 *
 **/
class BLCMDplace : public BLCommand {
	G4double x;
	G4double y;
	G4double z;
	G4String parent;
	G4String rename;
	G4String rotation;
	G4String coordinates;
	G4int copies;
	G4int front;
	BLElement *element;
	std::map<G4String,G4String> elementArgs;
public:
	BLCMDplace();

	G4String commandName() { return "place"; }
	
	int command(BLArgumentVector& argv, BLArgumentMap& namedArgs);

	void defineNamedArgs();

	/// special version of handleNamedArgs().
	int handleNamedArgs(BLArgumentMap& namedArgs);

	void printElementArgs();
};

BLCMDplace definePlace;

BLCMDplace::BLCMDplace() : BLCommand(), elementArgs()
{
	registerCommand(BLCMDTYPE_PLACE);
	setSynopsis("places an element into the current group (or world).");
	setDescription("Every element can be placed multiple times into the\n"
		"beamline. For most elements the geometrical centerpoint is\n"
		"placed; for polycone the local x=y=z=0 point is placed. "
		"If front is nonzero then the front of the element is "
		"placed (and any rotation is applied at the front).\n"
		"If z is specified, then the element is placed\n"
		"at that z position relative to the center of the enclosing\n"
		"group. If z is not specified, then the element is placed\n"
		"immediately downstream (higher z) of the previous element\n"
		"in the group, or at the upstream edge of the group if this\n"
		"is the first element in the group.\n\n"
		"The 'rename' argument can be used to change the name of the "
		"element (applies to traces and other uses of object names, "
		"such as the NTuple name of a virtualdetector). When placing "
		"into a group or other object, the rename argument should "
		"normally begin with '+' to include the parent's name; "
		"otherwise multiple placements of the parent will generate "
		"multiple objects with identical names -- that should "
		"be avoided for output objects like virtualdetector. "
		"Without a rename argument, the parent's name is included "
		"automatically.\n\n"
		"When multiple copies are placed, z refers to the first, and "
		"the rest are placed sequentially along z.\n"
		"When placing an element into the World group, Centerline\n"
		"coordinates are used unless coordinates=global is present.\n"
		"When centerline coordinates are used, the parameter 'Zcl' is "
		"set to the highest Z value used; this is normally the Z value "
		"for the front of the next element when placed sequentially "
		"(i.e. with no z value given).\n"
		"\nRotations:\n"
		"The rotation parameter can be used to rotate this element \n"
		"relative to "
		"the enclosing group. The object is rotated, not the axes.\n"
		"Rotations are specified as a comma-separated list of axes\n"
		"and angles (in degrees): rotate=Z90,X45 rotates first by 90\n"
		"degrees around Z and then by 45 degrees around X. The axes\n"
		"are the local X,Y,Z coordinate axes of the enclosing group\n"
		"(centerline or global coordinate axes for the World group);\n"
		"groups can be rotated when placed, and are rotated as a\n"
		"rigid unit (including children).\n\n"
		"If parent=name is present, then the name must be an element\n"
		"that accepts children, and it is used as the enclosing group;\n"
		"in this case the size of the group is implied by the size\n"
		"of the parent element, and z must be given (defaults to 0).\n"
		"Note that a given element cannot be the parent of any other\n"
		"element once it has been placed, so you must place children\n"
		"into their parent before placing their parent.\n\n"
		"If the special element 'OFFSET' is given, x, y, and z\n"
		"specify offsets for every following place command into the\n"
		"current group (incl. World), that gives a z position.");
	x = 0.0;
	y = 0.0;
	z = 0.0;
	parent = "";
	rename = "";
	rotation = "";
	coordinates = "Centerline";
	copies = 1;
	front = 0;
	element = 0;
}

int BLCMDplace::command(BLArgumentVector& argv, BLArgumentMap& namedArgs)
{
	if(argv.size() != 1) {
		printError("Invalid place command -- need one element name");
		return -1;
	}

	// handle offsets in x, y, and z
	if(argv[0] == "OFFSET") {
		x = 0.0;
		y = 0.0;
		z = 0.0;
		handleNamedArgs(namedArgs);
		G4ThreeVector offset(x,y,z);
		BLGroup::getCurrent()->setOffset(offset);
		printf("place   %-7s x=%.1f y=%.1f z=%.1f\n",
			"OFFSET",x,y,z);
		return 0;
	}

	element = BLElement::find(argv[0]);
	if(!element) {
		printError("place: cannot find element '%s'",argv[0].c_str());
		return -1;
	}

	// clear previous values
	x = UNDETERMINED;
	y = UNDETERMINED;
	z = UNDETERMINED;
	parent = "";
	rename = NO_RENAME;
	copies = 1;
	rotation = "";
	coordinates = "Centerline";
	front = 0;
	elementArgs.clear();

	handleNamedArgs(namedArgs);

	if(z != UNDETERMINED && front != 0)  {
		double len = element->getLength();
		if(len <= 0.0)
			printError("place: element length is <= 0 -- cannot use front=1");
		z += len/2.0;
	}

	if(rename != NO_RENAME && rename.find_first_of('+') != 0) {
		if(BLGroup::getCurrent() != BLGroup::getWorld())
			G4Exception("place command","No rename",JustWarning,
			  "Multiple placements of the enclosing group "
			  "will create multiple objects with identical names.");
		else if(parent != "")
			G4Exception("place command","No rename",JustWarning,
			  "Multiple placements of the parent "
			  "will create multiple objects with identical names.");
	}

	if(parent != "") {
		BLGroupElement *ge = BLGroupElement::find(parent);
		if(!ge) {
			if(BLElement::find(parent) != 0) {
				printError("place: parent '%s' cannot have children",
					parent.c_str());
			} else {
				printError("place: cannot find parent '%s'",
					parent.c_str());
			}
			return -1;
		}
		if(ge->getPlacedFlag()) {
			printError("place: parent '%s' has already been placed",
					parent.c_str());
			return -1;
		}
		if(z == UNDETERMINED && element->getLength() == 0) {
			printError("place: element length is <= 0 -- cannot be placed sequentially");
		}
		if(x == UNDETERMINED) x = 0.0;
		if(y == UNDETERMINED) y = 0.0;
		if(z == UNDETERMINED) z = 0.0;
		G4ThreeVector offset(x,y,z);
		G4RotationMatrix *rot = 0;
		if(rotation != "")
			rot = stringToRotationMatrix(rotation);
		// Jean-Francois Ostiguy's fix for front=1 to rotate around
		// the placement point (not the center).
		if(front != 0 && rot != 0) {
			G4ThreeVector dz(0.0,0.0,0.5*element->getLength());
			offset += (*rot) * dz - dz; 
		}
		for(int i=0; i<copies; ++i) {
		    ge->placeChild(element,rot,offset,rename);
		    offset[2] += element->getLength();
		}
		printf("place   %-7s parent=%s copies=%d x=%.1f y=%.1f z=%.1f ",
				argv[0].c_str(),parent.c_str(),copies,x,y,z);
		if(rename != NO_RENAME)
			printf("rename='%s'",rename.c_str());
		if(rotation != "")
			printf("rotation='%s'",rotation.c_str());
		printf("\n");
		element->setPlacedFlag();
		printElementArgs();
		return 0;
	}

	BLCoordinateType coordType = 
				BLCoordinates::getCoordinateType(coordinates);
	if(coordType!=BLCOORD_GLOBAL && coordType!=BLCOORD_CENTERLINE)
		printError("place: invalid coordinate type");

	// if placing into a group, indent by 2 spaces
	const char *p;
	if(BLGroup::getCurrent() == BLGroup::getWorld())
		p = "place  ";
	else
		p = "  place";

	if(z == UNDETERMINED) {
		if(rotation != "")
			printError("place without z being given, rotation is"
					" ignored");
		if(x != UNDETERMINED || y != UNDETERMINED)
			printError("place without z being given, x and y are"
						" ignored");
		if(element->getLength() == 0) {
			printError("place: element length is <= 0 -- cannot be placed sequentially");
		}
		for(int i=0; i<copies; ++i)
			BLGroup::getCurrent()->placeElement(element,rename,
						coordType==BLCOORD_GLOBAL);
		printf("%-7s %-7s copies=%d sequentially along z; ",
				p,argv[0].c_str(),copies);
		if(rename != NO_RENAME)
			printf("rename='%s'",rename.c_str());
	} else {
		if(x == UNDETERMINED) x = 0.0;
		if(y == UNDETERMINED) y = 0.0;
		G4ThreeVector offset(x,y,z);
		offset += BLGroup::getCurrent()->getOffset();
		G4RotationMatrix *rot = 0;
		if(rotation != "")
			rot = stringToRotationMatrix(rotation);
		// Jean-Francois Ostiguy's fix for front=1 to rotate around
		// the placement point (not the center).
		if(front != 0 && rot != 0) {
			G4ThreeVector dz(0.0,0.0,0.5*element->getLength());
			offset += (*rot) * dz - dz; 
		}
		for(int i=0; i<copies; ++i) {
		    BLGroup::getCurrent()->placeElement(element,rot,offset,
    					rename, coordType==BLCOORD_GLOBAL);
		    offset[2] += element->getLength();
		}
		printf("%-7s %-7s copies=%d x=%.1f y=%.1f z=%.1f ",
				p,argv[0].c_str(),copies,x,y,z);
		if(rename != NO_RENAME)
			printf("rename='%s'",rename.c_str());
		if(rotation != "")
			printf("rotation='%s'",rotation.c_str());
	}
	printf("\n");

	element->setPlacedFlag();

	printElementArgs();
	return 0;
}

void BLCMDplace::defineNamedArgs()
{
	argDouble(z,"z","Z position of element's center relative to the "
		"center of the enclosing group (mm).");
	argDouble(x,"x","X position of element's center [default=0] (mm)");
	argDouble(y,"y","Y position of element's center [default=0] (mm)");
	argString(parent,"parent","Parent element name (must accept children).");
	argString(rename,"rename","Name to use for this placement; '#' will "
		"be substituted by the number of placements. If the value "
		"begins with '+', it is replaced with the parent's name.");
	argInt(copies,"copies","Number of copies (placed sequentially along z).");
	argInt(front,"front","Nonzero to specify z for the front, not the center.");
	argString(rotation,"rotation","Rotation of this object.");
	argString(coordinates,"coordinates","Coordinates: global or centerline (default=c).");
}

int BLCMDplace::handleNamedArgs(BLArgumentMap& args)
{
	int retval = 0;
	argMode = PROCESS;
	bool isCloned = false;

	BLArgumentMap::iterator i;
	for(i=args.begin(); i!=args.end(); ++i) {
		argName = i->first;
		argValue = i->second;
		argFound = false;
		defineNamedArgs();
		if(!argFound) {
			// handle element args
			if(!isCloned) {
				element = element->clone();
				isCloned = true;
			}
			element->argName = argName;
			element->argValue = argValue;
			element->argFound = false;
			element->argMode = CHANGE;
			element->defineNamedArgs();
			if(!element->argFound) {
				printError("Invalid argument '%s' to place",
					argName.c_str());
				retval = -1;
			} else {
				elementArgs[argName] = argValue;
			}
		}
	}

	if(isCloned) element->argChanged();

	return retval;
}

void BLCMDplace::printElementArgs()
{
	G4String line;
	std::map<G4String,G4String>::iterator it;
	for(it=elementArgs.begin(); it!=elementArgs.end(); ++it) {
		G4String n = it->first;
		G4String v = it->second;
		line += (line.size()==0? "" : " ");
		line += n + "=" + v;
	}
	if(line.size() > 0) {
		G4String indent("                ");
		printf("%s",wrapWords(line,indent,indent).c_str());
	}
}
