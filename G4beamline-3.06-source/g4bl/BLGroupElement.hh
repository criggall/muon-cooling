//	BLGroupElement.hh
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

#ifndef BLGROUPELEMENT_HH
#define BLGROUPELEMENT_HH

#include <vector>
#include <map>
#include "G4RotationMatrix.hh"
#include "G4VPhysicalVolume.hh"
#include "BLElement.hh"

const G4String NO_RENAME("\x01");	// permits "" to be legal rename

/**	BLGroupElement - Base class for BLElement-s that can contain other 
 *			 BLElement-s
 **/
class BLGroupElement : public BLElement {
	typedef BLChild Child;
	std::vector<Child> child;
	static std::map<G4String,BLGroupElement*> mapGroupElement;
	BLGroupElement& operator=(const BLGroupElement&);	// undefined
	// BLGroup grows automatically, so it needs to handle children itself
	friend class BLGroup;
public:
	/// Default Constructor.
	BLGroupElement() : BLElement(), child() { }

	/// Copy Constructor for a BLGroupElement.
	BLGroupElement(const BLGroupElement& r) : BLElement(r), child(r.child) { }

	/// Destructor.
	~BLGroupElement() { }

	/// setName() adds this BLGroupElement to mapGroupElement; calls
	/// BLElement::setName().
	void setName(G4String _name);

	/// getNChildren() returns the number of children.
	int getNChildren() { return child.size(); }

	/// getChildren() returns the vector of children.
	const std::vector<BLChild>& getChildren() const { return child; }

	/// constructChildren() will construct the children
	/// The enclosing BLGroupElement must already be constructed, so
	/// this is normally called at the end of construct().
	virtual void constructChildren(G4LogicalVolume *parent, 
			G4String parentName,
			G4RotationMatrix *parentRotation,
			G4ThreeVector parentPosition);

	/// placeChild() will place a child within this GroupElement.
	/// The object being placed is rotated wrt the fixed axes of the
	/// group's local coordinates.
	virtual void placeChild(BLElement* element, G4RotationMatrix *rot,
		G4ThreeVector& offset, G4String rename);

	/// testGeometry() will check that all children are inside this
	/// element, and that no children intersect each other. It then
	/// calls testGeometry() for any children that are GroupElement-s.
	/// npoints is the # of points per child to test.
	/// returns the number of errors detected.
	int testGeometry(int npoints, G4double tolerance, bool visual=false,
				G4RotationMatrix rotation=G4RotationMatrix(),
				G4ThreeVector offset=G4ThreeVector());

	/// isGroupElement() returns true if this is a BLGroupElement.
	virtual G4bool isGroupElement() { return true; }

	/// isWithin() tests that a point is within this element, or is
	/// within tolerance of being inside.
	/// local[] is in local coordinates of this element.
	virtual G4bool isWithin(G4ThreeVector &local, G4double tolerance) = 0;

	/// find() will find a GroupElement by name
	static BLGroupElement *find(G4String name);
};

#endif // BLGROUPELEMENT_HH
