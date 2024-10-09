//	BLElement.cc
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

#include "Randomize.hh"

#include "BLElement.hh"

std::map<G4String,BLElement*> BLElement::mapElement;
bool BLElement::surfaceCheck = true;//enables the surface check in G4PVPlacement

void BLElement::setName(G4String _name)
{
	name = _name;
	mapElement[name] = this;
}

BLElement *BLElement::find(G4String _name)
{
	return mapElement[_name];
}

G4bool BLElement::allOK()
{
	std::map<G4String,BLElement*>::iterator i;
	for(i=mapElement.begin(); i!=mapElement.end(); ++i) {
		if(!i->second->isOK()) return false;
	}
	return true;
}

const std::vector<BLChild>& BLElement::getChildren() const
{
	static std::vector<BLChild> v;
	return v;
}

#define rand G4UniformRand

void BLElement::generateBox(unsigned int npoints, G4double width,
		G4double height, G4double length, std::vector<G4ThreeVector> &v)
{
	v.clear();
	G4double w=width/2.0, h=height/2.0, l=length/2.0;
	// the 8 corners
	v.push_back(G4ThreeVector(w,h,l));
	v.push_back(G4ThreeVector(w,h,-l));
	v.push_back(G4ThreeVector(w,-h,l));
	v.push_back(G4ThreeVector(w,-h,-l));
	v.push_back(G4ThreeVector(-w,h,l));
	v.push_back(G4ThreeVector(-w,h,-l));
	v.push_back(G4ThreeVector(-w,-h,l));
	v.push_back(G4ThreeVector(-w,-h,-l));
	// the 12 edge midpoints
	v.push_back(G4ThreeVector(0,h,l));
	v.push_back(G4ThreeVector(0,h,-l));
	v.push_back(G4ThreeVector(0,-h,l));
	v.push_back(G4ThreeVector(0,-h,-l));
	v.push_back(G4ThreeVector(w,0,l));
	v.push_back(G4ThreeVector(w,0,-l));
	v.push_back(G4ThreeVector(-w,0,l));
	v.push_back(G4ThreeVector(-w,0,-l));
	v.push_back(G4ThreeVector(w,h,0));
	v.push_back(G4ThreeVector(w,-h,0));
	v.push_back(G4ThreeVector(-w,h,0));
	v.push_back(G4ThreeVector(-w,-h,0));
	// the 6 face centers
	v.push_back(G4ThreeVector(w,0,0));
	v.push_back(G4ThreeVector(-w,0,0));
	v.push_back(G4ThreeVector(0,h,0));
	v.push_back(G4ThreeVector(0,-h,0));
	// BLCMDgenericbend requires these two be last:
	v.push_back(G4ThreeVector(0,0,l));
	v.push_back(G4ThreeVector(0,0,-l));
	// random points on the surfaces
	while(v.size() < npoints) {
		G4double x=0,y=0,z=0;
		switch(v.size() % 6) {
		case 0: x=w;  y=rand()*height-h; z=rand()*length-l; break;
		case 1: x=-w; y=rand()*height-h; z=rand()*length-l; break;
		case 2: x=rand()*width-w; y=h;  z=rand()*length-l; break;
		case 3: x=rand()*width-w; y=-h; z=rand()*length-l; break;
		case 4: x=rand()*width-w; y=rand()*height-h; z=l;  break;
		case 5: x=rand()*width-w; y=rand()*height-h; z=-l; break;
		}
		v.push_back(G4ThreeVector(x,y,z));
	}
}

void BLElement::generateTubs(unsigned int npoints, G4double innerRadius, 
		G4double outerRadius, G4double initialPhi, G4double finalPhi,
		G4double length, std::vector<G4ThreeVector> &v)
{
	v.clear();
	G4double l=length/2.0;

	// 8 points evenly spaced around outer (+inner) surface
	G4double dPhi = (finalPhi-initialPhi)/8.0;
	for(int i=0; i<9; ++i) {
		// if full cylinder, i=8 duplicates i=0
		if(i == 8 && dPhi >= 3.1415/4.0) break;
		G4double phi = (double)i*dPhi + initialPhi;
		G4double x = outerRadius * cos(phi);
		G4double y = outerRadius * sin(phi);
		v.push_back(G4ThreeVector(x,y,l));
		v.push_back(G4ThreeVector(x,y,-l));
		v.push_back(G4ThreeVector(x,y,0.0));
		if(innerRadius > 0.0) {
			x = innerRadius * cos(phi);
			y = innerRadius * sin(phi);
			v.push_back(G4ThreeVector(x,y,l));
			v.push_back(G4ThreeVector(x,y,-l));
			v.push_back(G4ThreeVector(x,y,0.0));
		}
	}
	// center of +z and =z ends, if innerRadius == 0
	if(innerRadius == 0.0) {
		v.push_back(G4ThreeVector(0.0,0.0,l));
		v.push_back(G4ThreeVector(0.0,0.0,-l));
	}
	// random points on the surface
	dPhi = finalPhi - initialPhi;
	G4double dR = outerRadius - innerRadius;
	while(v.size() < npoints) {
		G4double z=0.0,phi=0.0,r=0.0;
		switch(v.size() % 6) {
		case 0:	// inner surface
inner:			if(innerRadius > 0.0) {
				phi = rand()*dPhi + initialPhi;
				r = innerRadius;
				z = rand()*length - l;
				break;
			}
			// flow into case 1
		case 1:	// outer surface
outer:			phi = rand()*dPhi + initialPhi;
			r = outerRadius;
			z = rand()*length - l;
			break;
		case 2:	// +z end
			phi = rand()*dPhi + initialPhi;
			r = rand()*dR + innerRadius;
			z = l;
			break;
		case 3:	// -z end
			phi = rand()*dPhi + initialPhi;
			r = rand()*dR + innerRadius;
			z = -l;
			break;
		case 4:	// initialPhi face
			if(dPhi >= 3.1415*2.0) goto inner;
			phi = initialPhi;
			r = rand()*dR + innerRadius;
			z = rand()*length - l;
			break;
		case 5:	// finalPhi face
			if(dPhi >= 3.1415*2.0) goto outer;
			phi = finalPhi;
			r = rand()*dR + innerRadius;
			z = rand()*length - l;
			break;
		}
		v.push_back(G4ThreeVector(r*cos(phi),r*sin(phi),z));
	}
}
