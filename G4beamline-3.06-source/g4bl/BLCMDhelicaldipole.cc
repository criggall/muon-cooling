//	BLCMDhelicaldipole.cc
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
//  7/16/04 --- Add ICOOL mode (model=3) ky
//  9/14/04 --- Fix ICOOL FIELD ky
//                Compared field params in ICOOL and HelicalDipole modules. 
//                Both numbers are same.
//                A small discrepancy (0.1%) is in gamma function.  
//  9/15/04 --- Add Base Electric field ez ky
//  1/20/05 --- Add modulated bD and bQ (model 4) ky

#ifndef G4BL_GSL
int BLCMDhelicaldipole_dummy=0;
#else

#define _USE_MATH_DEFINES
#include <math.h>

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
#include "BLHelicalUtils.hh"

#define MARKER_SIZE 5	/* pixels */

/* how cylindrical coordinates stored */
#define cylPHI  0
#define cylRHO  1
#define cylZ    2

/* how rectangular coordinates stored */
#define rectX    0
#define rectY    1
#define rectZ    2

/* small area to avoid singularities */
#define ESSENTIALLY_ZERO   1.E-33
const G4double PI=3.141592653589793238462643;

/* pitch convention */
#define  LEFT_HANDED_THREAD   (-1)
#define RIGHT_HANDED_THREAD   (+1)

#define MuonsInc_PITCH_CONVENTION    RIGHT_HANDED_THREAD

#define MuonsInc_Target_Radius_mm    159.15494

#define KBB_ENABLE             1
#define KBB_DISABLE            0
#define KBB_bugz               KBB_DISABLE

/**	BLCMDhelicaldipole implements helical dipole magnet with a cylindrical
 *	field volume.
 *
 *      model 1
 *	The magnetic field is an ideal dipole that rotates helically along
 *	z. This is a pure field -- there is no iron or physical volume.
 *      model 2
 *      The field is determined from solutions of cylindrical Maxwellian.
 *      RPJ and YSD designed this. (MuNote0284)
 *      model 3
 *      The field is determined from solutions of cylindrical Maxwellian. 
 *      This time you take into account the boundary condition. 
 *      T. Tominaka et al designed this. (NIM A459:398)
 *      model 4
 *      The field strength can be modulated as z. 
 **/
class BLCMDhelicaldipole : public BLElement, public BLManager::RunAction {
	G4double radius;
	G4double length;
	G4double b;
        G4double bprime;
        G4double bpp;
        G4double rr;
	G4double lambda;
        G4double model;
        G4double psi0;
	G4double phi0;
	G4double Bsolenoid;
	G4double ez;
	BLCoordinateTransform global2local;
	G4Polymarker markers;
	friend class HelicalDipoleField;
public:
	/// Default constructor. Defines the command, args, etc.
	BLCMDhelicaldipole();

	/// Destructor.
	virtual ~BLCMDhelicaldipole() { }

	/// Copy constructor.
	BLCMDhelicaldipole(const BLCMDhelicaldipole& r);

	/// clone()
	BLElement *clone() { return new BLCMDhelicaldipole(*this); }

	/// commandName() returns "helicaldipole".
	G4String commandName() { return "helicaldipole"; }

	/// command() implements the helicaldipole command.
	int command(BLArgumentVector& argv, BLArgumentMap& namedArgs);

	/// defineNamedArgs() defines the named arguments for the command.
	void defineNamedArgs();

	/// construct() - construct the helicaldipole magnet
	void construct(G4RotationMatrix *relativeRotation,
			G4ThreeVector relativePosition, 
			G4LogicalVolume *parent, 
			G4String parentName,
			G4RotationMatrix *parentRotation,
			G4ThreeVector parentPosition);

	/// getLength() returns the fieldLength of the hd
	G4double getLength() { return length; }

	/// getWidth() returns the outer radius of the hd
	G4double getWidth() { return radius*2.0; }

	/// getHeight() returns the outer radius of the hd
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

BLCMDhelicaldipole defaultHelicalDipole;	// default object


/**	HelicalDipoleField represents one placement of a helicaldipole magnet.
 *
 **/
class HelicalDipoleField : public BLElementField {
	G4double radius;
	G4double halflength;
	G4double b;
        G4double bprime;
        G4double bpp;
        G4double rr;
        G4double psi0;
        G4double lambda;
	G4double phi0;
	G4double Bsolenoid;
        G4double model;
        G4double ez;
	BLCoordinateTransform global2local;
	G4RotationMatrix rotation;
public:
	/// constructor. 
	HelicalDipoleField(BLCoordinateTransform& _global2local, BLCMDhelicaldipole *hd);

	/// addFieldValue() adds the field for this solenoid into field[].
	/// point[] is in global coordinates.
	void addFieldValue(const G4double point[4], G4double field[6]) const;
};


// Default constructor - be sure to use the default constructor BLElement()
BLCMDhelicaldipole::BLCMDhelicaldipole() : BLElement(), BLManager::RunAction()
{
	// register the commandName(), and its synopsis and description.
	registerCommand(BLCMDTYPE_ELEMENT);
	setSynopsis("construct a helicaldipole magnet.");
	setDescription(
		"The field region is a cylinder with a helical dipole\n"
		"field plus a solenoid field.  The simple model=1 \n"
		"provides just a sine and cosine transverse dependence,\n"
		"while the maxwellian model=2 has both dipole and quadrupole\n"
		"terms.  Both the dipole scale bD [T] and quadrupole scale bQ [T/m]\n"
		"are now at rho=0; the user must determine the correct values externally.\n\n"
		"Note that this Element generates a magnetic field only,\n"
		"and only within the cylinder defined by length and radius.\n"
		"So it has no solid associated with it, and is invisible.\n");

	// provide initial values for fields
	radius = 0.0;
	length = 0.0;
	b = 0.0;
        bprime = 0.0;
	bpp    = 0.0;
	rr     = 159.155;
        model = 3. ;
	lambda = 0.0;
	phi0 = 0.0;
	psi0 = 0.0;
	Bsolenoid = 0.0;
	ez = 0.0;
}

// Copy constructor - be sure to use the copy constructor BLElement(r)
BLCMDhelicaldipole::BLCMDhelicaldipole(const BLCMDhelicaldipole& r) : BLElement(r), 
							BLManager::RunAction(r)
{
	// copy fields one at a time (transfers default values from the
	// default object to this new object).
	radius = r.radius;
	length = r.length;
	b = r.b;
        bprime = r.bprime;
        bpp    = r.bpp;
        rr     = r.rr;
	lambda = r.lambda;
	phi0 = r.phi0;
        model = r.model;
	Bsolenoid = r.Bsolenoid;
	ez = r.ez;
	psi0 = r.psi0;
}

int BLCMDhelicaldipole::command(BLArgumentVector& argv, BLArgumentMap& namedArgs)
{
        const char *CoDE="BLCMDhelicaldipole::command";

	if(argv.size() != 1) {
		printError("helicaldipole: Invalid command, must have name");
		return -1;
	}

	if(argv[0] == "default") {
		return defaultHelicalDipole.handleNamedArgs(namedArgs);
	}

	BLCMDhelicaldipole *t = new BLCMDhelicaldipole(defaultHelicalDipole);
	t->setName(argv[0]);
	int retval = t->handleNamedArgs(namedArgs);

	t->print(argv[0]);

	return retval;
}

void BLCMDhelicaldipole::defineNamedArgs()
{
	argDouble(radius,"radius","The radius of the field region (mm)",mm);
	argDouble(model,"model","The model of field calculated(simple=1, rpj and ysd model=2, ttominaka et al model=3), modulations in bd,bq,bz= 4");
	argDouble(length,"length","The length of the field region (mm)",mm);
	argDouble(b,"bD","The dipole magnitude at rho=0 (Tesla).",tesla);
	argDouble(lambda,"lambda","Helix period along the Z axis (mm).",mm);
	argDouble(phi0,"phi0","The phase of the XY field at the entrance (deg).",deg);
	argDouble(Bsolenoid,"Bsolenoid","The value of Bsolenoid (Tesla).",tesla);
        argDouble(bprime,"bQ","The quadrupole magnitude at rho=0 (Tesla).",tesla/meter);
        argDouble(bpp,"bs","The sextupole magnitude at rho=0 (Tesla).",tesla/meter/meter);
        argDouble(rr,"rr","Reference radius (mm)",mm);
	argDouble(psi0,"psi0","The offset between the dipole term and the quadrupole term (Degrees).",deg);
	argDouble(ez,"ez","The base electric field inside the helix channel (GV/m).",1000*megavolt/meter);
}

void BLCMDhelicaldipole::construct(G4RotationMatrix *relativeRotation,
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

	HelicalDipoleField *p = new HelicalDipoleField(global2local,this);
	BLGlobalField::getObject()->addElementField(p);

	printf("BLCMDhelicaldipole::Construct %s parent=%s relZ=%.1f globZ=%.1f\n"
			"\tzmin=%.1f zmax=%.1f\n",
		thisname.c_str(),parentName.c_str(),relativePosition[2],
		globalPosition[2], zmin,zmax);

	BLManager::getObject()->registerRunAction(this,false);
}

void BLCMDhelicaldipole::BeginOfRunAction(const G4Run *run)
{
	markers.clear();
}

void BLCMDhelicaldipole::EndOfRunAction(const G4Run *run)
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



HelicalDipoleField::HelicalDipoleField(BLCoordinateTransform& _global2local,
					BLCMDhelicaldipole *hd) :
					BLElementField(), rotation()
{
	radius = hd->radius;
	halflength = hd->length/2.0;
	b = hd->b;
        bprime = hd->bprime;
        bpp    = hd->bpp;
        rr     = hd->rr;
        model = hd->model;
  	lambda = hd->lambda;
	phi0 = hd->phi0;
	psi0 = hd->psi0;
	Bsolenoid = hd->Bsolenoid;
	ez = hd->ez;
	global2local = _global2local;
	rotation = global2local.getRotation().inverse();

//printf("HelicalDipoleField: radius=%.1f mm, halflength=%.1f mm, b=%.4f T, bprime=%.4f T/m\n\tmodel=%.1f, lambda=%.1f mm, phi0=%.1f deg, psi0=%.1f deg, Bsolenoid=%.4f T\n",radius/mm,halflength/mm,b/tesla,bprime/(tesla/meter),model,lambda/mm,phi0/deg,psi0/deg,Bsolenoid/tesla);

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


/* The following two routines were added from a public domain repository */
/* which were written in 1996 and allow double precision modified bessel */
/* to be used in the field map for the helical dipole magnet RPJ 24JAN04 */
/* (slightly modified KBB 15APR04) */ 

G4double modified_bessel0( G4double x )
{ 
  /* Returns the modified Bessel function I_0(x) for any real x */
  /* modified bessel function of order 0 */ 
  /* slightly modified DBESI0 for type double, KBB 4/04  */
    double ax;
    double anser;
    double y;
 
    if ((ax=fabs(x)) < 3.75) {
      y = x/3.75;
      y*=y;
      anser = 1+y*(3.5156229+y*(3.0899424+y*(1.2067492
               +y*(0.2659732+y*(0.360768e-1+y*0.45813e-2)))));
    }
    else {
      y=3.75/ax;
      anser = (exp(ax)/sqrt(ax))*(0.39894228+y*(0.1328592e-1
              +y*(0.225319e-2+y*(-.157565e-2+y*(0.916281e-2
              +y*(-0.2057706e-1+y*(0.2635537e-1+y*(-0.1647633e-1
              +y*0.392377e-2)))))))); 
    }
    return anser;
}


G4double modified_bessel1( G4double x )
{   
  /* Returns the modified Bessel function I_1(x) for any real x */
  /* modified bessel function of order 1 */ 
  /* slightly modified DBESI1 for type double, KBB 4/04 */
    double ax;
    double anser;
    double y;
 
    if ((ax=fabs(x)) < 3.75) {
      y = x/3.75;
      y*=y;
      anser = ax*(0.5+y*(0.87890594+y*(0.51498869+y*(0.15084934
	      +y*(0.2658733e-1+y*(0.301532e-2+y*0.32411e-3))))));
    }
    else {
      y=3.75/ax;
      anser = 0.2282967e-1+y*(-0.2895312e-1+y*(0.1787654e-1
	      -y*0.420059e-2));
      anser = 0.39894228+y*(-0.3988024e-1+y*(-0.362018e-2
	      +y*(0.163801e-2+y*(-0.1031555e-1+y*anser))));
      anser *= (exp(ax)/sqrt(ax));
         }
    return x < 0.0 ? -anser : anser;
}

G4double modified_besselN( int o, G4double x )
{
  const char *CoDE="modified_besselN";
  int n;
  G4double modified_besselN_x( int n, G4double x ); /*prototype*/
  G4double Io=0;

  n= o-1;
  if( o>1 ) 
      Io= -2*n * modified_besselN_x(n,x) + modified_besselN(n-1,x);
  else if( o==1 )
      Io= modified_bessel1(x);
  else if( o==0 )
      Io= modified_bessel0(x);
  else
      printf("%s: bad order=%i\n", CoDE, n);

  return( Io );
}

G4double modified_besselN_x( int n, G4double x )
{
  G4double In_t;
  int i;

  if( -ESSENTIALLY_ZERO<x && x<ESSENTIALLY_ZERO )
    { 
      for(In_t=1/2.,i=2; i<=n; ++i)  In_t*= (x/2)/i;
    }
  else
    {
      In_t= modified_besselN(n,x)/x;
    }

  return( In_t );
}


G4double derivative_modified_besselN( int n, G4double x )
{
  const char *CoDE="derivative_modified_besselN";

  if( n<1 )  printf("%s: order#n < 1\n", CoDE);
  return( modified_besselN(n-1,x) - n * modified_besselN_x(n,x) );
}

G4ThreeVector DIPOLEFIELD( G4double bd /* magnitude at zero radius [T] */,
			   G4double k /* +|- pitch of magnet [radians/mm] */,
			   G4double rho /* local angular position [mm] */,
			   G4double psiangle /* net phase angle [radians] */ )
{
	  const int PHI=cylPHI,RHO=cylRHO,Z=cylZ;
	  G4double t,I0,I1_t,psi;
          G4ThreeVector Bcyl;

	  t= k * rho;

	  I0= modified_besselN(0,t);

	  I1_t= modified_besselN_x(1,t);

	  /* Helical Cooling Channel Simnulation Parameters and Fields */
	  /* 6D Cooling Note                                           */
	  /* using cylindrical components  - equation 1.2              */ 

	  Bcyl[PHI]= 2 * bd * I1_t    * cos(psiangle);

	  Bcyl[RHO]= 2 * bd * (I0-I1_t) * sin(psiangle);

	  /* -2 * bd * I1(t) * cos(phiangle) <==> -t * Bcyl[PHI] */

	  Bcyl[Z]=   - k * rho * Bcyl[PHI];
	  //	  printf("1; (%f, %f) %e %e %e \n",rho,psiangle,Bcyl[0],Bcyl[1],Bcyl[2]);
	  return Bcyl;
}


G4ThreeVector QUADRUPOLEFIELD(G4double bprime /* gradient [T/m#] */,
			      G4double k /* +|- pitch [radians/mm] */, 
			      G4double rho /* current radial position [mm] */,
                              G4double psiangle /* local angular position [radians] */ )
{
  /* 
   #note that the NUMERICAL value of bprime is NOT in [T/m];
   one must convert bprime/(tesla/meter) to get units of [T/m],
   but [T=0.001] and [meter=1000]! 

   corrected a mistake in the formula 4/27/04 KBB,RJ
  */

  const char *CoDE="QUADRUPOLEFIELD";
  const int PHI=cylPHI,RHO=cylRHO,Z=cylZ;
  static int messageCount=0;
  G4double t,I1,I2_t,Bprime;
  G4ThreeVector Bcyl;


  if( 0==messageCount++ )  printf("%s: version of 4/27/04\n", CoDE);


  //      <----- T/m--------->   <----- 1/1000 --->   <---- T/mm ------>

  Bprime= bprime/(tesla/meter) * (millimeter/meter) * (tesla/millimeter);

  t= 2 * k * rho;
  
  I2_t= modified_besselN_x(2,t);
    
  I1= modified_besselN(1,t);
  
  /* Helical Cooling Channel Simnulation Parameters and Fields */
  /* 6D Cooling Note                                           */
  /* using cylindrical components  - equation 1.7              */ 
  /* BUT THERE IS A TYPO in equantion 1.7 -  it should say     */
  /* 1/k NOT 1/2 in the 2nd line!                              */
  
  /* I (2*k*rho)/(k*k*rho)=> I (2*k*rho)/(2*k*rho) * 2/k */
  /*  2                       2                          */
 
  Bcyl[PHI] = 2/k * Bprime * I2_t * cos(2*psiangle); 
  
  Bcyl[RHO] = 1/k * Bprime * (I1 - 2*I2_t) * sin(2*psiangle);  // <== fixed! 
  
  Bcyl[Z] = - k * rho * Bcyl[PHI];

  if( KBB_bugz && 0.95*MuonsInc_Target_Radius_mm<rho && rho<1.05*MuonsInc_Target_Radius_mm ) 
    {
      printf("%s: bprime=%f[T/m] %f[native] t=%f I2_t=%f psiangle=%f ", 
	     CoDE, bprime/(tesla/meter), bprime, t, I2_t, psiangle );
      printf("k=%f/mm rho=%fmm Bcyl=(%f,%f,%f)[native]=(%f,%f,%f)[T]\n", 
	     k*mm, rho/mm, Bcyl[PHI], Bcyl[RHO], Bcyl[Z],
	     Bcyl[PHI]/tesla, Bcyl[RHO]/tesla, Bcyl[Z]/tesla);
    }

  //  printf("2; (%f, %f) %e %e %e \n",rho,psiangle,Bcyl[0],Bcyl[1],Bcyl[2]);
  return Bcyl;
}

void HelicalDipoleField::addFieldValue(const G4double point[4], G4double field[6]) const
{
  const char *CoDE="HelicalDipoleField::addFieldValue";
  const int X=rectX,Y=rectY,Z=rectZ;
  const int PHI=cylPHI,RHO=cylRHO;                      /* notation */
  const int node=0;
  int ini;
  G4ThreeVector B,BcylDipole,BcylQuad,BcylSext,BxyzDipole,BxyzQuad,BxyzSext;
  G4ThreeVector global(point[X],point[Y],point[Z]);
  G4ThreeVector local;
  G4double kz,phi,bd,rho,kH,Bprime,Bpp,pZ;
  G4double mul[3];

  global2local.getLocal(local,global);                  /*fetch local location in mm*/
  rho = sqrt( local[X]*local[X] + local[Y]*local[Y] );  /*mm*/
  
  if( rho>radius || fabs(local[Z]) > halflength)  return;

  for(ini=0;ini<=2;ini++){
    BcylDipole[ini]=0.0;
    BxyzDipole[ini]=0.0;
    BcylQuad[ini]=0.0;
    BxyzQuad[ini]=0.0;
    BcylSext[ini]=0.0;
    BxyzSext[ini]=0.0;
    B[ini]=0.0;
  }
  
  /* the Muons,Inc. preferred convention -> lambda>0 == right handed screw */

  kH= MuonsInc_PITCH_CONVENTION *2*PI/lambda;
  //  kH= sqrt(4.0*PI*PI/lambda/lambda);
  //  kz= kH * (local[Z]+halflength)*lambda/sqrt(lambda*lambda);
  kz= kH * (local[Z]+halflength);
  pZ = (local[Z]+halflength)*1e-3;
  phi= atan2( local[Y], local[X]);

  bd = b;   /* the radial correction for bQ, bprime, will now be done externally by the user */ 

  /* Simple SIN/COS Dipole Model */
  if( int(model)==1 )
    {
      B= SIMPLEFIELD(b,kz+phi0,Bsolenoid,model);
    }
  
  /* Dipole Maxwellian Model  */
  else if ( int(model)==2 )
    {   
      BcylDipole= DIPOLEFIELD( bd, kH, rho, phi-kz+phi0 );
      BxyzDipole = CYLTOCARTESIAN(BcylDipole,phi);

      BcylQuad= QUADRUPOLEFIELD( bprime, kH, rho, phi-kz+phi0+psi0 );
      BxyzQuad = CYLTOCARTESIAN(BcylQuad,phi);

      B[X]= BxyzDipole[X] + BxyzQuad[X];
      B[Y]= BxyzDipole[Y] + BxyzQuad[Y];
      B[Z]= BxyzDipole[Z] + BxyzQuad[Z] + Bsolenoid;

      if( KBB_bugz )
	{
	  printf( "%s: BcylDipole=(%f,%f,%f) BcylQuad=(%f,%f,%f)\n", 
		  CoDE, 
		  BcylDipole[PHI]/tesla, BcylDipole[RHO]/tesla, BcylDipole[Z]/tesla,
		  BcylQuad[PHI]/tesla, BcylQuad[RHO]/tesla, BcylQuad[Z]/tesla );
	  
	  printf( "%s: Bxyz=(%f+%f=%f, %f+%f=%f, %f+%f+%f=%f)\n",
		  CoDE, 
		  BxyzDipole[X]/tesla, BxyzQuad[X]/tesla, B[X]/tesla,
		  BxyzDipole[Y]/tesla, BxyzQuad[Y]/tesla, B[Y]/tesla,
		  BxyzDipole[Z]/tesla, BxyzQuad[Z]/tesla, Bsolenoid, B[Z]/tesla );
	}
    }

  // Extended Dipole Maxwellian (ICOOL) Model
  else if ( int(model)==3 )
    {
      // TO DO: Extend to multiple components
      //        Use more realistic formula
      BcylDipole = bd*ICOOLFIELD(1,rho/mm,phi-kz-phi0,kH/mm,rr/mm);
      BxyzDipole = CYLTOCARTESIAN(BcylDipole,phi);
      Bprime= bprime*1e3;
      BcylQuad = Bprime*ICOOLFIELD(2,rho/mm,phi-kz-phi0+psi0,kH/mm,rr/mm);
      BxyzQuad   = CYLTOCARTESIAN(BcylQuad,phi);
      Bpp= bpp*1e6;
      BcylSext = Bpp*ICOOLFIELD(3,rho/mm,phi-kz-phi0+psi0,kH/mm,rr/mm);
      BxyzSext = CYLTOCARTESIAN(BcylSext,phi);

      B[PHI] =BcylDipole[PHI]+BcylQuad[PHI]+BcylSext[PHI];
      B[RHO] =BcylDipole[RHO]+BcylQuad[RHO]+BcylSext[RHO];
      B[Z]   =BcylDipole[Z]+BcylQuad[Z]+BcylSext[Z];

      B[X]= BxyzDipole[X] + BxyzQuad[X] + BxyzSext[X];  
      B[Y]= BxyzDipole[Y] + BxyzQuad[Y] + BxyzSext[Y];
      B[Z]= BxyzDipole[Z] + BxyzQuad[Z] + BxyzSext[Z] + Bsolenoid;

    }

  // Field map calculation for MANX experiment
  else if ( int(model)==4 )
    {
      // TO DO: Extend to multiple components
      //        Use more realistic formula
      // Add new terms to modulate bd,bq,bz (1/20/05 KY)
      // mul[0]=bd, mul[1]=bq, mul[2]=bz

      // Original values
      //mul[0] = 3.24178-0.369821*pZ+1.24318e-3*pZ*pZ-1.91452e-3*pZ*pZ*pZ;
      //mul[1] = -0.616436+6.69245e-2*pZ-5.36372e-4*pZ*pZ+3.23624e-4*pZ*pZ*pZ;
      //mul[2] = -11.5474+1.27729*pZ-1.03380e-2*pZ*pZ+6.48178e-3*pZ*pZ*pZ;

      // Values for 300 MeV/c input
      mul[2] = -8.53621+1.33517*pZ+3.51035e-2*pZ*pZ+6.48178e-3*pZ*pZ*pZ;
      mul[0] = 2.34731-0.409774*pZ-9.50945e-3*pZ*pZ-1.85969e-3*pZ*pZ*pZ;
      mul[1] = -0.442645+8.65446e-2*pZ+6.15396e-4*pZ*pZ+2.78646e-4*pZ*pZ*pZ;

      mul[0]*= bd;
      BcylDipole = mul[0]*ICOOLFIELD(1,rho/mm,phi-kz-phi0,kH/mm,159.155/mm);
      BxyzDipole = CYLTOCARTESIAN(BcylDipole,phi);
      mul[1]*= 1e3*bprime;
      BcylQuad = mul[1]*ICOOLFIELD(2,rho/mm,phi-kz-phi0+psi0,kH/mm,159.155/mm);
      BxyzQuad   = CYLTOCARTESIAN(BcylQuad,phi);
//@@@ mul[3] is OUT OF BOUNDS! It is also UNINITIALIZED!
//@@@ mul[3]*= bpp*1e6;
//@@@ BcylSext = mul[3]*ICOOLFIELD(3,rho/mm,phi-kz-phi0+psi0,kH/mm,rr/mm);
G4Exception("BLCMDhelicaldipole","Invalid model",FatalException, "Internal coding error");
//@@@ BxyzSext = CYLTOCARTESIAN(BcylSext,phi);
      mul[2]*= Bsolenoid;

      //      printf("z=%e m1=%e m2=%e m3=%e\n",pZ,mul[0],mul[1],mul[2]);

      B[PHI] =BcylDipole[PHI]+BcylQuad[PHI]+BcylSext[PHI];
      B[RHO] =BcylDipole[RHO]+BcylQuad[RHO]+BcylSext[RHO];
      B[Z]   =BcylDipole[Z]+BcylQuad[Z]+BcylSext[Z];

      B[X]= BxyzDipole[X] + BxyzQuad[X] + BxyzSext[X];  
      B[Y]= BxyzDipole[Y] + BxyzQuad[Y] + BxyzSext[Y];
      B[Z]= BxyzDipole[Z] + BxyzQuad[Z] + BxyzSext[Z] + mul[2];

    }
  
  else
    {
      fprintf(stderr,"HelicalDipoleField::addFieldValue bad model#%d (1|2)\n", int(model) );
      fflush(stderr);
    }
  
  /* Rotation if applicable */
  
  if(global2local.isRotated())  B = rotation * B;
  
  field[0] += B[X];        /* update the field */
  field[1] += B[Y];
  field[2] += B[Z];

  field[5] += ez;

}

#endif // G4BL_GSL
