//	BLGroupElement.cc
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

#include "BLGroupElement.hh"
#include "BLMarkers.hh"

std::map<G4String,BLGroupElement*> BLGroupElement::mapGroupElement;

void BLGroupElement::setName(G4String _name)
{
	mapGroupElement[_name] = this;
	BLElement::setName(_name);
}

BLGroupElement *BLGroupElement::find(G4String _name)
{
	return mapGroupElement[_name];
}

void BLGroupElement::constructChildren(G4LogicalVolume *parent, 
			G4String parentName,
			G4RotationMatrix *parentRotation,
			G4ThreeVector parentPosition)
{
	static int level=0;

	if(++level > 64) {
		G4Exception("BLGroupElement","Nesting > 64",FatalException, "");
	}

	int counter = 1;
	std::vector<Child>::iterator i;
	for(i=child.begin(); i!=child.end(); ++i) {
		BLElement *e = i->element;
		G4String savename = e->getName();
		G4String pname = parentName;
		// handle i->rename
		if(i->rename != NO_RENAME) {
			if(i->rename.find('+') == 0) {
				e->setName(i->rename.substr(1));
			} else {
				e->setName(i->rename);
				pname = "";
			}
		}
		// handle '#' in name
		G4String::size_type j = e->getName().find('#');
		if(j != e->getName().npos) {
			char tmp[32];
			sprintf(tmp,"%d",counter++);
			G4String tmpname(e->getName());
			tmpname.replace(j,1,tmp);
			e->setName(tmpname);
		}
		i->name = pname + e->getName();
		e->construct(&i->rot,i->offset,parent,pname,parentRotation,
							parentPosition);
		e->setName(savename);
	}

	--level;
}

void BLGroupElement::placeChild(BLElement* element, G4RotationMatrix *rot,
		G4ThreeVector& offset, G4String rename)
{
	Child c;
	c.element = element;
	c.rename = rename;
	if(rot) c.rot = *rot;
	c.offset = offset;
	child.push_back(c);
}

int BLGroupElement::testGeometry(int npoints, G4double tolerance, bool visual,
			G4RotationMatrix rotation, G4ThreeVector offset)
{
	static BLMarkers *failureMarkers = 0;
	if(failureMarkers == 0) failureMarkers = new BLMarkers("1,0,1",12);

	// no sense testing an empty group
	if(child.size() == 0) return 0;

	printf("Testing geometry for children of %s '%s':\n",
		commandName().c_str(), getName().c_str());

	int err = 0;
	std::vector<G4ThreeVector> points;
	G4RotationMatrix inverse=rotation.inverse();
	int icounter = 1;
	std::vector<Child>::iterator i;
	for(i=child.begin(); i!=child.end(); ++i) {
		BLElement *ie = i->element;
		G4String iname = ie->getName();
		// handle i->rename
		if(i->rename != NO_RENAME)
			iname = i->rename;
		// handle '#' in iname
		G4String::size_type k = iname.find('#');
		if(k != iname.npos) {
			char tmp[32];
			sprintf(tmp,"%d",icounter++);
			iname.replace(k,1,tmp);
		}
		ie->generatePoints(npoints,points);
		// convert local coords of ie to local coords of this
		bool isRotated = !i->rot.isIdentity();
		for(k=0; k<points.size(); ++k) {
			if(isRotated) points[k] = i->rot * points[k];
			points[k] += i->offset;
		}
		if(visual) {
			BLMarkers *m = new BLMarkers("0,.8,0",6);
			for(k=0; k<points.size(); ++k) {
				G4ThreeVector pos = points[k];
				pos = inverse * pos + offset; // local->global
				m->addMarker(pos);
			}
		}
		int errParent=0;
		for(k=0; k<points.size(); ++k) {
		    if(!isWithin(points[k],tolerance)) {
			G4ThreeVector global = inverse * points[k] + offset;
			failureMarkers->addMarker(global);
			++errParent;
		    }
		}
		if(errParent > 0) {
			char tmp[1024];
			sprintf(tmp,"%s '%s' extends outside %s '%s'",
				ie->commandName().c_str(),
				iname.c_str(),
				commandName().c_str(),
				getName().c_str());
			G4Exception("geometry","Geometry Error",
							JustWarning, tmp);
			++err;
		}
		int jcounter = 1;
		std::vector<Child>::iterator j;
		for(j=child.begin(); j!=child.end(); ++j) {
			if(j == i) continue;
			BLElement *je = j->element;
			G4String jname = je->getName();
			// handle j->rename
			if(j->rename != NO_RENAME)
				jname = j->rename;
			// handle '#' in jname
			G4String::size_type k = jname.find('#');
			if(k != jname.npos) {
				char tmp[32];
				sprintf(tmp,"%d",jcounter++);
				jname.replace(k,1,tmp);
			}
			isRotated = !j->rot.isIdentity();
			int errSibling=0;
			for(k=0; k<points.size(); ++k) {
				G4ThreeVector pt = points[k] - j->offset;
				if(isRotated) pt = j->rot.inverse() * pt;
				if(!je->isOutside(pt,tolerance)) {
				    G4ThreeVector global = inverse * points[k] 
				    				+ offset;
				    failureMarkers->addMarker(global);
				    ++errSibling;
				}
			}
			if(errSibling > 0) {
				char tmp[1024];
				sprintf(tmp,"%s '%s' intersects %s '%s'",
					ie->commandName().c_str(),
					iname.c_str(),
					je->commandName().c_str(),
					jname.c_str());
				G4Exception("geometry","Geometry Error",
				    		 JustWarning, tmp);
				++err;
			}
		}
	}
	for(i=child.begin(); i!=child.end(); ++i) {
		BLElement *ie = i->element;
		if(!ie->isGroupElement()) continue;
		BLGroupElement *ge = dynamic_cast<BLGroupElement*>(ie);
		if(ge) {
			G4ThreeVector off = offset + rotation.inverse() * i->offset;
			G4RotationMatrix rot = i->rot.inverse() * rotation;
			err += ge->testGeometry(npoints,tolerance,visual,
								rot,off);
		}
	}

	return err;
}
