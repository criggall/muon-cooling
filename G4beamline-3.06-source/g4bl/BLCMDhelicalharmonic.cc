//	BLCMDhelicalharmonic.cc
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

Original by Katsuya Yonehara, used with permission.
*/
//  03/05/10 --- Helical Harmonic of given order by Vasiliy Morozov

#ifdef G4BL_GSL

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif

#include <math.h>
#include <gsl/gsl_sf_bessel.h>

#include "G4VisAttributes.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4Color.hh"
#include "G4UserLimits.hh"
#include "G4Polymarker.hh"
#include "G4VVisManager.hh"

#include "BLElement.hh"
#include "BLElementField.hh"
#include "BLGlobalField.hh"
#include "BLParam.hh"
#include "BLManager.hh"

#ifndef MARKER_SIZE
#define MARKER_SIZE 5	/* pixels */
#endif

//	BLCMDhelicalharmonic implements magnetic helical harmonic of given order

class BLCMDhelicalharmonic : public BLElement, public BLManager::RunAction {
	G4double radius;
	G4double length;
	int n; 
	G4double b;
	G4double lambda;
	G4double phi0;
	BLCoordinateTransform global2local;
	G4Polymarker markers;
	friend class HelicalHarmonicField;
public:
	/// Default constructor. Defines the command, args, etc.
	BLCMDhelicalharmonic();

	/// Destructor.
	virtual ~BLCMDhelicalharmonic() { }

	/// Copy constructor.
	BLCMDhelicalharmonic(const BLCMDhelicalharmonic& r);

	/// clone()
	BLElement *clone() { return new BLCMDhelicalharmonic(*this); }

	/// commandName() returns "helicalharmonic".
	G4String commandName() { return "helicalharmonic"; }

	/// command() implements the helicalharmonic command.
	int command(BLArgumentVector& argv, BLArgumentMap& namedArgs);

	/// defineNamedArgs() defines the named arguments for the command.
	void defineNamedArgs();

	/// construct() - construct the helicalharmonic magnet
	void construct(G4RotationMatrix *relativeRotation,
			G4ThreeVector relativePosition, 
			G4LogicalVolume *parent, 
			G4String parentName,
			G4RotationMatrix *parentRotation,
			G4ThreeVector parentPosition);

	/// getLength() returns the fieldLength of the hh
	G4double getLength() { return length; }

	/// getWidth() returns the outer radius of the hh
	G4double getWidth() { return radius*2.0; }

	/// getHeight() returns the outer radius of the hh
	G4double getHeight() { return radius*2.0; }

	/// getSurveyPoint() returns points in LOCAL coordinates.
	G4ThreeVector getSurveyPoint(int index) {
		if(index == 0) return G4ThreeVector(0.0,0.0,-getLength()/2.0);
		if(index == 1) return G4ThreeVector(0.0,0.0,getLength()/2.0);
		throw "UNIMPLEMENTED";
	}

	/// isOK() returns true.
	G4bool isOK() { return true; }

	/// isOutside() from BLElement.
	bool isOutside(G4ThreeVector &local, G4double tolerance)
		{ return true; }

	/// generatePoints() from BLElement.
	void generatePoints(int npoints, std::vector<G4ThreeVector> &v)
		{ v.clear(); }

	/// BeginOfRunAction() from BLManager::RunAction.
	void BeginOfRunAction(const G4Run *run);

	/// EndOfRunAction() from BLManager::RunAction.
	void EndOfRunAction(const G4Run *run);
};

BLCMDhelicalharmonic defaultHelicalHarmonic;	// default object


/**	HelicalHarmonicField represents one placement of a helicalharmonic magnet.
 *
 **/
class HelicalHarmonicField : public BLElementField {
	G4double radius;
	G4double halflength;
	int n; 
	G4double b;
	G4double lambda;
	G4double phi0;
	BLCoordinateTransform global2local;
	G4RotationMatrix rotation;
public:
	/// constructor. 
	HelicalHarmonicField(BLCoordinateTransform& _global2local, BLCMDhelicalharmonic *hh);

	/// addFieldValue() adds the field for this hh into field[].
	/// point[] is in global coordinates.
	void addFieldValue(const G4double point[4], G4double field[6]) const;
};


// Default constructor - be sure to use the default constructor BLElement()
BLCMDhelicalharmonic::BLCMDhelicalharmonic() : BLElement(), BLManager::RunAction()
{
	// register the commandName(), and its synopsis and description.
	registerCommand(BLCMDTYPE_ELEMENT);
	setSynopsis("construct a helicalharmonic magnet.");
	setDescription(
		"Creates a cylindrical region containing the field of \n"
		"a magnetic helical harmonic of given order [n]. \n"
		"The field is defined by the value of the (n-1) order \n"
		"derivative [b] of the vertical field component (when \n"
		"the initial phase is 0) with respect to the horizontal \n"
		"coordinate at the center of the helix: \n"
		"  b=d^(n-1)B_phi/dr^(n-1) @ [r=0 & phi-k*z+phi0=0], \n" 
		"where k=2*pi/lambda is the helix's wave number, \n" 
		"[lambda] is the length of the helix's period, and \n" 
		"phi0 is the initial phase. \n"
		"The field components in the cylindrical frame are given by: \n"
		"  B_phi=(2/(n*k))^(n-1)*b*(I[n-1](n*k*r)-I[n+1](n*k*r))* \n"
		"        cos(n*(phi-k*z+phi0)), \n"
		"  B_r  =(2/(n*k))^(n-1)*b*(I[n-1](n*k*r)+I[n+1](n*k*r))* \n"
		"        sin(n*(phi-k*z+phi0)), \n"
		"  B_z  =-2*(2/(n*k))^(n-1)*b*I[n](n*k*r)*cos(n*(phi-k*z+phi0)), \n"
		"where I[n](x) is the modified Bessel function of the first kind \n"
		"of order [n]. \n\n"
		"This element must be placed (via the place command).\n\n"
		"Note that this Element generates magnetic field only,\n"
		"and only within the cylinder defined by length and radius.\n"
		"So it has no solid associated with it, and is invisible.\n");

	// provide initial values for fields
	radius = 0.0;
	length = 0.0;
	n=1; 
	b = 0.0;
	lambda = 1.0;
	phi0 = 0.0;
}

// Copy constructor - be sure to use the copy constructor BLElement(r)
BLCMDhelicalharmonic::BLCMDhelicalharmonic(const BLCMDhelicalharmonic& r) : BLElement(r), 
							BLManager::RunAction(r)
{
	// copy fields one at a time (transfers default values from the
	// default object to this new object).
	radius = r.radius;
	length = r.length;
	n = r.n;
	b = r.b;
	lambda = r.lambda;
	phi0 = r.phi0;
}

int BLCMDhelicalharmonic::command(BLArgumentVector& argv, BLArgumentMap& namedArgs)
{
	const char *CoDE="BLCMDhelicalharmonic::command";

	if(argv.size() != 1) {
		printError("helicalharmonic: Invalid command, must have name");
		return -1;
	}

	if(argv[0] == "default") {
		return defaultHelicalHarmonic.handleNamedArgs(namedArgs);
	}

	BLCMDhelicalharmonic *t = new BLCMDhelicalharmonic(defaultHelicalHarmonic);
	t->setName(argv[0]);
	t->handleNamedArgs(namedArgs);	// call it twice to ensure that n is
					// set for units of b.
	int retval = t->handleNamedArgs(namedArgs);

	t->print(argv[0]);

	return retval;
}

void BLCMDhelicalharmonic::defineNamedArgs()
{
	argDouble(radius,"radius","The radius of the field region (mm)",mm);
	argDouble(length,"length","The length of the field region (mm)",mm);
	argInt(n,"n","Order of helical harmonic (i.e. n=1 for dipole)");
	// see comment in command() about setting n
	argDouble(b,"b","(n-1)-order derivative of the field at the center "
		"(T/m^(n-1))",tesla/pow(meter,(n-1)));
	argDouble(lambda,"lambda","Helix period along the Z axis (mm).",mm);
	argDouble(phi0,"phi0","The phase of the XY field at the entrance (rad).");
}

void BLCMDhelicalharmonic::construct(G4RotationMatrix *relativeRotation,
	G4ThreeVector relativePosition, 
	G4LogicalVolume *parent, 
	G4String parentName,
	G4RotationMatrix *parentRotation,
	G4ThreeVector parentPosition)

{
	G4String thisname = parentName+getName();

	// get globalRotation and globalPosition
	G4RotationMatrix *globalRotation = 0;
	if(relativeRotation && parentRotation) {
		globalRotation = 
		new G4RotationMatrix(*parentRotation * *relativeRotation);
	} else if(relativeRotation) {
		globalRotation = relativeRotation;
	} else if(parentRotation) {
		globalRotation = parentRotation;
	}
	G4ThreeVector globalPosition(relativePosition + parentPosition);
	if(parentRotation)
		globalPosition = *parentRotation * relativePosition +
				parentPosition;

	global2local = BLCoordinateTransform(globalRotation,globalPosition);

	G4double zmin = globalPosition[2]-getLength()/2.0;
	G4double zmax = globalPosition[2]+getLength()/2.0;

	HelicalHarmonicField *p = new HelicalHarmonicField(global2local,this);
	BLGlobalField::getObject()->addElementField(p);

	printf("BLCMDhelicalharmonic::Construct %s parent=%s relZ=%.1f globZ=%.1f\n"
			"\tzmin=%.1f zmax=%.1f\n",
		thisname.c_str(),parentName.c_str(),relativePosition[2],
		globalPosition[2], zmin,zmax);

	BLManager::getObject()->registerRunAction(this,false);
}

void BLCMDhelicalharmonic::BeginOfRunAction(const G4Run *run)
{
	markers.clear();
}

void BLCMDhelicalharmonic::EndOfRunAction(const G4Run *run)
{
	G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
	if (!pVVisManager) return;

#ifdef STUB // omit markers
	double dz = lambda/10.0;
	int n = (int)(length/dz) + 1;
	for(int i=0; i<n; ++i) {
		G4double local[4], global[4];
		G4double phi = phi0 + i * dz * 2.0*pi/lambda;
		local[0] = radius/2.0*cos(phi);
		local[1] = radius/2.0*sin(phi);
		local[2] = -length/2.0 + i * dz;
		local[3] = 0.0;
		global2local.getGlobal(local,global);
		G4ThreeVector point(global[0],global[1],global[2]);
		markers.push_back(point);
	}
	markers.SetMarkerType(G4Polymarker::circles);
	markers.SetScreenSize(MARKER_SIZE);
	markers.SetFillStyle(G4VMarker::filled);
	G4VisAttributes va(G4Colour(1.,1.,1.));  // white
	markers.SetVisAttributes(&va);

	pVVisManager->Draw(markers);
#endif // STUB
}

HelicalHarmonicField::HelicalHarmonicField(BLCoordinateTransform& _global2local,
					BLCMDhelicalharmonic *hd) :
					BLElementField(), rotation()
{
	radius = hd->radius;
	halflength = hd->length/2.0;
	n = hd->n;
	b = hd->b;
	lambda = hd->lambda;
	phi0 = hd->phi0;
	global2local = _global2local;
	rotation = global2local.getRotation().inverse();

	// set global bounding box
	G4double local[4], global[4];
	local[3] = 0.0;
	for(int i=0; i<2; ++i) {
		local[0] = (i==0 ? -1.0 : 1.0) * radius;
		for(int j=0; j<2; ++j) {
			local[1] = (j==0 ? -1.0 : 1.0) * radius;
			for(int k=0; k<2; ++k) {
				local[2] = (k==0 ? -1.0 : 1.0) * halflength;
				global2local.getGlobal(local,global);
				setGlobalPoint(global);
			}
		}
	}
}

void HelicalHarmonicField::addFieldValue(const G4double point[4], G4double field[6]) const
{
	const char *CoDE="HelicalHarmonicField::addFieldValue";
	const int X=0,Y=1,Z=2;
	const int PHI=0,RHO=1;                      /* notation */

	G4ThreeVector global(point[X],point[Y],point[Z]);
	G4ThreeVector local;

	global2local.getLocal(local,global);   /*fetch local location in mm*/
	G4double phi = atan2( local[Y], local[X]);
	G4double rho = sqrt( local[X]*local[X] + local[Y]*local[Y] );  /*mm*/
  
	if( rho>radius || fabs(local[Z]) > halflength)  return;

	G4double cosphi = cos(phi); 
	G4double sinphi = sin(phi); 
	G4double k = 2*pi/lambda; 
	G4double kz = k*(local[Z]+halflength); 
	G4double psi = n*(phi-kz+phi0); 
	G4double cospsi = cos(psi); 
	G4double sinpsi = sin(psi); 
	G4double nkrho = n*k*rho;
	G4double Inm1=gsl_sf_bessel_In(n-1,nkrho); 
	G4double In=gsl_sf_bessel_In(n,nkrho); 
	G4double Inp1=gsl_sf_bessel_In(n+1,nkrho); 
	G4double bb=b*pow(2/(n*k),(n-1)); 
	double Bcyl_phi=bb*(Inm1-Inp1)*cospsi; 
	double Bcyl_rho=bb*(Inm1+Inp1)*sinpsi;
	G4ThreeVector B;
	B[Z]=-2*bb*In*cospsi; 
	B[X]=Bcyl_rho*cosphi-Bcyl_phi*sinphi; 
	B[Y]=Bcyl_rho*sinphi+Bcyl_phi*cosphi;

	/* Rotation if applicable */ 
	if(global2local.isRotated())  B = rotation * B;
  
	field[0] += B[X];        /* update the field */
	field[1] += B[Y];
	field[2] += B[Z];
}

#else // G4BL_GSL
int dummyhelicalharmonic=0;
#endif // G4BL_GSL
