//	BLGroup.cc
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
//	NOTE: Separating this into a command class and an infrastructure class
//	was not easy, so class BLGroup does both. In keeping with the other
//	classes, BLCMDgroup.cc contains the command functions, and BLGroup.cc
//	contains the infrastructure functions.

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4PVPlacement.hh"
#include "G4UserLimits.hh"
#include "G4Material.hh"
#include "BLGroup.hh"
#include "BLParam.hh"
#include "BLCoordinates.hh"

BLSetParam Zcl("Zcl","0.0","Last centerline Z position used (updated continuously)");

BLGroup *BLGroup::world = 0;
BLGroup *BLGroup::current = 0;
G4VPhysicalVolume *BLGroup::worldPhysVol = 0;

extern BLGroup defaultGroup;

void BLGroup::placeElement(BLElement* element, G4RotationMatrix *rot,
			G4ThreeVector& offset, G4String rename, G4bool global)
{
	BLGroupElement::Child c;
	c.element = element;
	c.rename = rename;

	// convert centerline coordinates to global coordinates, if world
	if(this == world && !global) {
		BLCoordinates::getCurrentGlobal(offset,c.offset);
		if(rot)
			c.rot = *BLCoordinates::getCurrentRotation() * *rot;
		else
			c.rot = *BLCoordinates::getCurrentRotation();
	} else {
		if(rot) c.rot = *rot;
		c.offset = offset;
	}

	// add this child
	child.push_back(c);

	// handle renamme
	if(rename == NO_RENAME)
		c.rename = element->getName();
	else if(rename.find('+') == 0)
		c.rename = getName()+rename.substr(1);

	// expand group to element's size.
	G4ThreeVector corner(element->getWidth()/2.0,element->getHeight()/2.0,
						element->getLength()/2.0);
	G4ThreeVector parent;
	// loop over all 8 corners of the bounding box of the element.
	for(int i=0; i<2; ++i) {
	    corner[0] = -corner[0];
	    for(int j=0; j<2; ++j) {
		corner[1] = -corner[1];
	    	for(int k=0; k<2; ++k) {
			corner[2] = -corner[2];
			if(rot)
				parent = *rot * corner;
			else
				parent = corner;
			G4double x = offset[0]+parent[0];
			if(!fixedWidth)
				setMinWidth(2.0*fabs(x));
			else if(fabs(x) > width/2.0)
			    printError("Element '%s' extends outside '%s' in x",
			    	c.rename.c_str(),getName().c_str());
			G4double y = offset[1]+parent[1];
			if(!fixedHeight)
				setMinHeight(2.0*fabs(y));
			else if(fabs(y) > height/2.0)
			    printError("Element '%s' extends outside '%s' in y",
			    	c.rename.c_str(),getName().c_str());
			G4double z = offset[2]+parent[2];
			if(!fixedLength)
				setMinLength(this==world ? fabs(z) : z);
			else if(fabs(z) > length/2.0)
			    printError("Element '%s' extends outside '%s' in z",
			    	c.rename.c_str(),getName().c_str());
			if(i+j+k == 0 || z > lastZ)
				lastZ = z;
			if(this == world && !global) {
				G4ThreeVector gp;
				BLCoordinates::getCurrentGlobal(offset+parent,gp);
				setMinWidth(2.0*fabs(gp[0]));
				setMinHeight(2.0*fabs(gp[1]));
				setMinLength(gp[2]);
			}
		}
	    }
	}

	// update z positions used
	if(lastZ > highZ)
		highZ = lastZ;
	if(this == world && !global)
		Param.setParam("Zcl",lastZ);
}

void BLGroup::placeElement(BLElement* element, G4ThreeVector& offset, 
						G4String rename, G4bool global)
{
	placeElement(element,0,offset,rename,global);
}

void BLGroup::placeElement(BLElement* element, G4String rename, G4bool global)
{
	G4double l = element->getLength();
	G4ThreeVector offset(0.0,0.0,lastZ + l/2.0);
	placeElement(element,0,offset,rename,global);
}

void BLGroup::end()
{
	// compute length if not set by argument, and adjust child offsets
	if(length == 0.0) {
		if(this == world) {
			// z=0 is center
			length = highZ * 2.0;
		} else {
			// z=0 is where we started
			length = highZ;
			std::vector<BLGroupElement::Child>::iterator i;
			for(i=child.begin(); i!=child.end(); ++i) {
				i->offset[2] -= length/2.0;
			}
		}
	}

	// expand world by 20.1357 cm
	if(this == world) {
		printf("\nWorld size (before incrementing by 201.357 mm): "
		    "%.1f H  %.1f W  %.1f L\n",height/mm,width/mm,length/mm);
		width += 20.1357*cm;
		height += 20.1357*cm;
		length += 20.1357*cm;
	}

	// pop current group
	if(current != world && prevCurrent != 0)
		current = prevCurrent;
	prevCurrent = 0;
}

BLGroup *BLGroup::getCurrent()
{
	if(!current) current = getWorld();
	return current;
}

BLGroup *BLGroup::getWorld()
{
	if(!world) {
		world = new BLGroup(defaultGroup);
		world->setName("World");
	}
	return world;
}

G4VPhysicalVolume *BLGroup::constructWorld()
{
	G4ThreeVector abs;
	getWorld()->construct(0,abs,0,"",0,abs);
	worldPhysVol->GetLogicalVolume()->
			SetVisAttributes(new G4VisAttributes(false));
	return worldPhysVol;
}

G4bool BLGroup::isOutside(G4ThreeVector &local, G4double tolerance)
{
	if(radius > 0.0)
		return sqrt(local[0]*local[0]+local[1]*local[1]) > radius-tolerance ||
			fabs(local[2]) > length/2.0-tolerance;
	else
		return fabs(local[0]) > width/2.0-tolerance ||
			fabs(local[1]) > height/2.0-tolerance ||
			fabs(local[2]) > length/2.0-tolerance;
}

void BLGroup::generatePoints(int npoints, std::vector<G4ThreeVector> &v)
{
	if(radius > 0.0)
		generateTubs(npoints, 0.0, radius, 0.0, 2.0*pi, length, v);
	else
		generateBox(npoints,width,height,length,v);
}

G4bool BLGroup::isWithin(G4ThreeVector &local, G4double tolerance)
{
	if(radius > 0.0)
		return sqrt(local[0]*local[0]+local[1]*local[1]) < radius+tolerance ||
			fabs(local[2]) < length/2.0+tolerance;
	else
		return fabs(local[0]) < width/2.0+tolerance &&
			fabs(local[1]) < height/2.0+tolerance &&
			fabs(local[2]) < length/2.0+tolerance;
}
