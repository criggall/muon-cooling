//	BLGroup.hh
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

#ifndef BLGROUP_HH
#define BLGROUP_HH

#include <vector>
#include "G4RotationMatrix.hh"
#include "G4VPhysicalVolume.hh"
#include "G4ThreeVector.hh"
#include "BLGroupElement.hh"


/**	BLGroup - define one group (container for BLElement-s)
 *
 *	NOTE: the class for the group command is BLGroup, not BLCMDgroup.
 *	This class is too intertwinded between command and infrastructure,
 *	so a single class was used.
 *
 *	Implements the group command, which begins defining the contents
 *	of the group; the endgroup command ends the group's definition, and
 *	computes the size of the group (if it wasn't specified via arguments).
 *	The group expands to hold its contents, unless its size is specified
 *	in arguments to the group command. The World is a expandable group.
 *
 *	If radius is set to 0, uses a Tubs with radius equal to the larger
 *	of the widest or highest element placed into the group; if radius
 *	is set >0 that is the fixed radius of the Tubs. If radius is not
 *	set then uses a Box.
 *
 *	Note that for an expanding group, initially its front is at z=0,
 *	and no object can extend to z<0. The engroup command makes its length
 *	be highZ, and shifts all children forward by highZ/2.
 **/
class BLGroup : public BLGroupElement {
	static BLGroup *current;
	static BLGroup *world;
	static G4VPhysicalVolume *worldPhysVol;
	BLGroup *prevCurrent;
	G4double highZ;
	G4double lastZ;
	G4double length;
	G4double width;
	G4double height;
	G4double radius;
	G4String material;
	G4String color;
	G4double maxStep;
	G4bool fixedLength;
	G4bool fixedWidth;
	G4bool fixedHeight;
	G4ThreeVector offset;	// exclusively for the place command
	friend class BLCMDendgroup;
public:
	/// Default Constructor. registers the command and provides help text.
	BLGroup();

	/// Copy Constructor for a real BLGroup.
	BLGroup(BLGroup& r) : BLGroupElement(r), offset(0,0,0)
		{ prevCurrent=0; highZ=0.0; lastZ=0.0; length=r.length;
		  width=r.width; height=r.height; material=r.material;
		  radius=r.radius; color = r.color; maxStep = r.maxStep; 
		  fixedLength = r.fixedLength;
		  fixedWidth = r.fixedWidth;
		  fixedHeight = r.fixedHeight; }

	/// Destructor.
	~BLGroup();

	/// clone() - invalid.
	BLElement *clone() { return 0; }

	/// commandName() returns "group";
	virtual G4String commandName() { return "group"; }

	/// command() implements the group command.
	int command(BLArgumentVector& argv, BLArgumentMap& namedArgs);

	/// defineNamedArgs() defines the arguments to the command.
	void defineNamedArgs();

	/// construct() will construct an instance of this group, and construct
	/// its contents recursively.
	virtual void construct(G4RotationMatrix *relativeRotation,
			G4ThreeVector relativePosition, 
			G4LogicalVolume *parent, 
			G4String parentName,
			G4RotationMatrix *parentRotation,
			G4ThreeVector parentPosition);

	/// getLength() returns this group's Length along the Z axis.
	G4double getLength() { return length; }

	/// getWidth() returns this group's Width along the X axis.
	G4double getWidth() { return width; }

	G4ThreeVector getSurveyPoint(int index) {
		if(index == 0) return G4ThreeVector(0.0,0.0,-length/2.0);
		if(index == 1) return G4ThreeVector(0.0,0.0,length/2.0);
		throw "UNIMPLEMENTED";
	}

	/// getHeight() returns this group's height along the Y axis.
	G4double getHeight() { return height; }

	/// setMinWidth() will increase the width, if necessary.
	void setMinWidth(G4double minWidth) 
		{ if(width < minWidth) width = minWidth; }

	/// setMinHeight() will increase the height, if necessary.
	void setMinHeight(G4double minHeight) 
		{ if(height < minHeight) height = minHeight; }

	/// setMinLength() will increase the length, if necessary.
	void setMinLength(G4double minLength) 
		{ if(highZ < minLength) highZ = minLength; }

	/// setMaterial() sets the name of the group's material.
	void setMaterial(G4String mat) { material = mat; }

	/// isOK() returns true; individual elements are responsible for
	/// their own status.
	G4bool isOK() { return true; }

	/// isOutside() from BLElement.
	bool isOutside(G4ThreeVector &local, G4double tolerance);

	/// generatePoints from BLElement.
	void generatePoints(int npoints, std::vector<G4ThreeVector> &v);

	/// isWithin() from BLGroupElement.
	bool isWithin(G4ThreeVector &local, G4double tolerance);

	/// placeElement() will place an element within this group.
	/// With a rotation, an offset must be given, and the length of the
	/// group must have been specified via an argument.
	/// Called by the place command.
	/// The object being placed is rotated wrt the fixed axes of the
	/// group's local coordinates.
	/// global=true means use global coordinates for the world group,
	/// not centerline coords.
	virtual void placeElement(BLElement* element, G4RotationMatrix *rot,
		G4ThreeVector& offset, G4String rename, G4bool global=false);

	/// placeElement() will place an element within this group.
	/// With an offset, the length of the group must have been specified via
	/// an argument.
	/// Called by the place command.
	virtual void placeElement(BLElement* element, G4ThreeVector& offset, 
				G4String rename, G4bool global=false);

	/// placeElement() will place an element within this group.
	/// With no offset, this element is placed immediately downstream
	/// of the previous element. If necessary, the endgroup command
	/// will compute the size of the group and readjust offsets so
	/// they are relative to the center of the group.
	/// Called by the place command.
	virtual void placeElement(BLElement* element, G4String rename, 
							G4bool global=false);

	/// end() will compute the size of the group, if necessary, and
	/// pop the current group.
	void end();

	/// getOffset() returns the offset for the place command.
	G4ThreeVector& getOffset() { return offset; }

	/// setOffset() sets the offset for the place command.
	void setOffset(G4ThreeVector& _offset) { this->offset = _offset; }

	/// getCurrent() returns a pointer to the current group (innermost
	/// active group command). Creates current and world if necessary.
	/// Used by the place command.
	static BLGroup *getCurrent();

	/// getWorld() returns a pointer to the World group.
	///  Creates current and world if necessary.
	static BLGroup *getWorld();

	/// constructWorld() constructs the beamline and returns a pointer
	/// to its world PhysicalVolume.
	static G4VPhysicalVolume *constructWorld();
private:
	BLGroup& operator=(const BLGroup&);	// undefined
};

#endif // BLGROUP_HH
