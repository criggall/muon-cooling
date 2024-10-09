//	BLElement.hh
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

#ifndef BLELEMENT_HH
#define BLELEMENT_HH

#include <vector>
#include "G4ThreeVector.hh"
#include "G4LogicalVolume.hh"

#include "BLCoordinateTransform.hh"
#include "BLCommand.hh"

class BLElement;

/**	struct BLChild represents a child BLElement.
 **/
struct BLChild {
	G4RotationMatrix rot;
	G4ThreeVector offset;
	BLElement *element;
	G4String rename;
	G4String name;
	BLChild() : rot(), offset(), element(), rename(), name() { }
};

/**	class BLElement - interface class for all g4beamline Elements
 *
 *	BLElement defines elements to be placed into the geometry of g4beamline.
 *
 *	Normally an element implementation will derive a class from this one,
 *	and if it has an EM field, also derives from BLElementField. Usually
 *	all are in a single file named for the element command class.
 *	The default constructor of the derived class is used to implement
 *	the element command (see BLCommand), and the copy constructor
 *	is used by the element command to create instances of the element.
 *
 *	Note the constructor of the derived class must contain enough 
 *	information so the size of the element is determined (getLength(), 
 *	etc.). Or at least these must be known before the element command
 *	returns.
 *
 *	A non-default instance of this class represents a specific type
 *	of the element, complete with all argument values; it is created
 *	by the element command of the derived class. This instance can
 *	be placed multiple times by the place command and by place commands
 *	applied to groups in which the instance has been placed. After the
 *	command file has been read (i.e. all elements have been created and
 *	placed), then the World group is constructed, and that results in a
 *	traversal of the placement tree and the construction of all elements.
 *	Thus a single instance of this class can appear multiple times at
 *	multiple locations within the overall geometry.
 *
 *	When element argments are given on the place command, clone() is used
 *	to copy the element, and then defineNamedArgs() is called for
 *	each element argument. Then the cloned element is placed.
 *
 *	NOTE: Derived classes MUST use the copy constructor, not the
 *	default constructor! ONLY the default object should use the default
 *	constructor.
 **/
class BLElement : public BLCommand {
	static std::map<G4String,BLElement*> mapElement;
	G4String name;
	bool placed;
protected:
	static bool surfaceCheck; // enables the surface check in G4PVPlacement
public:
	/// Default constructor.
	BLElement() : BLCommand(), name() { placed = false; }

	/// Destructor.
	virtual ~BLElement() { }

	/// Copy constructor.
	BLElement(const BLElement& r) : BLCommand(r), name(r.name) 
		{ placed = r.placed; }

	/// clone() will clone a copy of an element. Used when element
	/// arguments are given in a place command.
	virtual BLElement *clone() = 0;

	/// getName() returns the element's name; do not confuse this with
	/// commandName(). Not used by the default instance.
	virtual G4String getName() const { return name; }

	/// setName() sets the element's name;
	/// Adds the BLElement to mapElement. Used by instances of the derived
	/// class representing real elements; not used by the default instance.
	/// virtual so derived base classes can keep a list of their elements.
	virtual void setName(G4String _name);

	/// find() finds the BLElement with a given name.
	static BLElement *find(G4String name);

	/// allOK() returns true iff all elements are ready to process beam.
	static G4bool allOK();

	// General functions used by elements:
	
	// virtual functions to be defined for individual BLElement-s:

	/// getSolid() will return the G4VSolid for this element.
	/// NOTE: only simple elements implement this, all complex elements
	/// return 0. Used by BLCMDboolean.
	virtual G4VSolid *getSolid() { return 0; }

	/// getMaterialName()  returns the name of the material.
	virtual G4String getMaterialName() const { return "Vacuum"; }

	/// construct() will construct a physical implementation of this 
	/// element at a specific location within the overall geometry.
	/// The name of the implementation should be parentName plus the name
	/// of the element; the physical volume should be linked into the 
	/// parent in the usual way. relativeRotation is an active rotation
	/// of this element wrt the parent's coordinate axes, and
	/// relativePosition is in the parent's coords as usual (if the parent
	/// is the world, the centerline coordinate transform is included) --
	/// the rotation for G4VPlacement is relativeRotation->inverse().
	/// parentRotation and parentPosition are the rotation and position
	/// of the parent, expressed in global coordinates -- most elements
	/// can ignore them, as they are normally used to construct a
	/// BLCoordinateTransform from global coordinates to the element's
	/// local coordinates (e.g. for GlobalField).
	virtual void construct(G4RotationMatrix *relativeRotation,
			G4ThreeVector relativePosition, 
			G4LogicalVolume *parent, 
			G4String parentName,
			G4RotationMatrix *parentRotation,
			G4ThreeVector parentPosition) = 0;

	/// getLength() returns this element's Length along the Z axis.
	virtual G4double getLength() = 0;

	/// getWidth() returns this element's Width along the X axis.
	/// For asymmetric elements, the width is twice the element's
	/// futhest point from X=0.
	virtual G4double getWidth() = 0;

	/// getHeight() returns this element's height along the Y axis.
	/// For asymmetric elements, the width is twice the element's
	/// futhest point from Y=0.
	virtual G4double getHeight() = 0;

	/// isOK() returns true iff this element is ready to process beam.
	/// It is queried AFTER the reference particle is tracked, so it
	/// can reflect the status of tuning performed by the reference 
	/// particle.
	virtual G4bool isOK() = 0;

	/// isGroupElement() returns true if this is a BLGroupElement.
	/// That is, if this element can be the parent of other elements.
	virtual G4bool isGroupElement() { return false; }

	/// getChildren() returns a vector of children. If this
	/// element cannot have children, an empty vector is returned.
	virtual const std::vector<BLChild>& getChildren() const;

	/// getSurveyPoint() returns a point in LOCAL coordinates for
	/// reporting by the survey command. index=0 for front, index=1
	/// for rear.
	virtual G4ThreeVector getSurveyPoint(int index) = 0;

	/// getPlacedFlag() returns true if this element has been placed.
	bool getPlacedFlag() { return placed; }

	/// setPlacedFlag() sets the flag for placing this element.
	void setPlacedFlag(bool flag=true) { placed = flag; }

	// Geometry test functions.

	/// generatePoints() will generate points for the testGeometry()
	/// function (of BLGroupElement) for this element. It generates
	/// however many points are necessary to test extremes of its surface,
	/// plus randomly distributed ones on the surface, up to npoints total.
	/// Each point is in the local coordinates of the element.
	virtual void generatePoints(int npoints, std::vector<G4ThreeVector> &v)
		= 0;

	/// isOutside() returns true if the point is outside this element,
	/// or is within tolerance of being outside.
	/// local[] is the point in local coordinates of this element.
	virtual G4bool isOutside(G4ThreeVector &local, G4double tolerance) = 0;

	/// generateBox() generates geometry test points for a Box
	/// (a convenience method for derived classes to use)
	void generateBox(unsigned int npoints, G4double width, G4double height,
			G4double length, std::vector<G4ThreeVector> &v);

	/// generateTubs() generates geometry test points for a Tubs
	/// (a convenience method for derived classes to use)
	void generateTubs(unsigned int npoints, G4double innerRadius, 
		G4double outerRadius, G4double initialPhi, G4double finalPhi,
			G4double length, std::vector<G4ThreeVector> &v);
};

#endif // BLELEMENT_HH
