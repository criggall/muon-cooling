//	BLCMDrfdevice.cc
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

Mar-Jun 2011- K.B.Beard, beard@muonsinc.com - extensive modifications 
              for better finding timing, phasing, and gradient,
	      support for multicell structures 3jun2011

KBB - It is important to note that when copies of "rfdevice"'s are
      used, multiple rfdeviceField:: share an RFdevice:: - so
      the final result of autoTiming is stored in RFdevicefield::
      not RFdevice::.  
      Multiple "tune" passes require that RFdevicefield:: be re-autoTimed 
      for each pass, as the initial conditions would have changed.
*/

#ifdef G4BL_GSL

#include <vector>
#include <stdio.h>

#include "G4VisAttributes.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4Color.hh"
#include "G4UserLimits.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4SteppingManager.hh"
#include "G4Material.hh"

#include "G4Transportation.hh"
#include "gsl/gsl_sf_bessel.h"

#include "BLElement.hh"
#include "BLElementField.hh"
#include "BLGlobalField.hh"
#include "BLManager.hh"
#include "BLParam.hh"
#include "BLFieldMap.hh"
#include "BLTune.hh"
#include "BLKillTrack.hh"

extern void g4bl_exit(int);

static long long global_2TraceCounter=1;

static G4bool TrackOfInterest=true;

static G4double lastFound_maxGradient=0;


#define TINY_NONZERO_GRADIENT  (pi*pi*1.E-10*megavolt/meter)  //avoid NaN's in G4Transport
                                                              //very unlikely # that anybody
                                                              //might actually use

const G4int ITERATION_LIMIT = 1000;	// limit on iterations/rfdevice autoTiming & adjustments


//bits to represent what is fixed
#define ATFIX_OFFSET       1
#define ATFIX_OFFSETINCR   2
#define ATFIX_GRADIENT     4
#define ATFIX_PHASE        8
#define ATFIX_OUTPUT      16   //only iff 1 output
#define ATFIX_DE          32
#define ATFIX_P           64
#define ATFIX_DT         128
#define ATFIX_XDEFL      256
#define ATFIX_YDEFL      512
#define ATFIX_ZLOCAL    1024
#define ATFIX_ANYOUT     (ATFIX_DE|ATFIX_P|ATFIX_DT|ATFIX_XDEFL|ATFIX_YDEFL)

#define ESSENTIALLY_UNCHANGED    1.E-6               //ignore

enum AutoTimingMethod { AUTOTIMING_OFF, AUTOTIMING_ATZ, 
			AUTOTIMING_MINENERGY, AUTOTIMING_NOENERGY, AUTOTIMING_MAXENERGY, 
			AUTOTIMING_MINTIME, AUTOTIMING_NOMTIME, AUTOTIMING_MAXTIME, 
			AUTOTIMING_MAXDXP, AUTOTIMING_NODXP, AUTOTIMING_MINDXP,
			AUTOTIMING_MAXDYP, AUTOTIMING_NODYP, AUTOTIMING_MINDYP,
			AUTOTIMING_UNKNOWN };

enum AutoTimingWhatToAdjust { AUTOADJUST_UNKNOWN, AUTOADJUST_PHASE, 
			      AUTOADJUST_GRADIENT, AUTOADJUST_OUTPUT, AUTOADJUST_UNNEEDED };


enum AutoTimingState { ATWORKING_UNKNOWN, ATWORKING_RESET, ATWORKING_OFFSET, ATWORKING_BYPASS, 
		       ATWORKING_BASE, ATWORKING_ESTIMATE, ATWORKING_FINETUNE, 
		       ATWORKING_FINAL, ATWORKING_DONE, ATWORKING_INFLUX };


const char *autotimingstate2text(AutoTimingState state)
{
  //returns pointer to text form of state
  switch(state) 
    {
    case ATWORKING_UNKNOWN:
      return("ATWORKING_UNKNOWN");
    case ATWORKING_RESET:
      return("ATWORKING_RESET");
    case ATWORKING_OFFSET:
      return("ATWORKING_OFFSET");
    case ATWORKING_BYPASS:
      return("ATWORKING_BYPASS");
    case ATWORKING_BASE:
      return("ATWORKING_BASE");
    case ATWORKING_ESTIMATE:
      return("ATWORKING_ESTIMATE");
    case ATWORKING_FINETUNE:
      return("ATWORKING_FINETUNE");
    case ATWORKING_FINAL:
      return("ATWORKING_FINAL");
    case ATWORKING_DONE:
      return("ATWORKING_DONE");
    case ATWORKING_INFLUX:
      return("ATWORKING_INFLUX");
    }

  return("?");
}

const char *state2desc(AutoTimingState state)
{
  return     state==ATWORKING_UNKNOWN ? "unknown" :
	     state==ATWORKING_RESET ? "reset" : 
	     state==ATWORKING_OFFSET ? "offset" : 
	     state==ATWORKING_BYPASS ? "bypass" :
	     state==ATWORKING_BASE ? "base" : 
	     state==ATWORKING_ESTIMATE ? "estimate" : 
	     state==ATWORKING_FINETUNE ? "finetune" :
	     state==ATWORKING_FINAL ? "final" :
             state==ATWORKING_DONE ? "done" :
             state==ATWORKING_INFLUX ? "influx" : "?????";
}


typedef struct            //snapshot at entry to timing Volume
{
  G4double q;                   /* particle charge */
  G4double x;                   /* global x location */
  G4double y;                   /* global y location */
  G4double z;                   /* global z location */
  G4double Px;                  /* local X momentum */
  G4double Py;                  /* local Y momentum */
  G4double Pz;                  /* local Z momentum */
  G4double t;                   /* global time */
  G4double KE;                  /* kinetic energy */
  G4double Ptot;                /* |momentum| */
  G4double v;                   /* |velocity| */
  G4double m;                   /* rest mass */
  G4double estTransit;          /* estimated transit time from initial velocity */
  G4double timeOff;             /* current time offset of RF */
  G4bool within;                /* whether within timing volume */
  AutoTimingWhatToAdjust what;  /* what to adjust */
} AutoTimingSnapshot;                     /* KBB save info upon entering timing volume */


typedef struct
{
  G4double dE;       //energy gain
  G4double Pout;     //total momentum output
  G4double dT;       //transit time
  G4double Xdefl;    //Xdeflection
  G4double Ydefl;    //Ydeflection
} AutoTimingChange;


AutoTimingChange autoTimingIn2Out(AutoTimingSnapshot in, AutoTimingSnapshot out)
{
  AutoTimingChange dif;

  dif.dE= out.KE - in.KE;
  dif.Pout= out.Ptot;
  dif.dT= out.t - in.t;
  
  dif.Xdefl= asin((in.Pz*out.Px-in.Px*out.Pz)/in.Ptot/out.Ptot);  
  dif.Ydefl= -asin((in.Pz*out.Py-in.Py*out.Pz)/in.Ptot/out.Ptot);

  return(dif);
}


G4int howmanyfixed(G4int whatisfixed)
{
  // return # of fixed quantities in bit-packed list (excluding ZLOCAL)

  int N=0;
  if(whatisfixed & ATFIX_P) ++N; 
  if(whatisfixed & ATFIX_DE) ++N;
  if(whatisfixed & ATFIX_DT) ++N;
  if(whatisfixed & ATFIX_XDEFL) ++N;
  if(whatisfixed & ATFIX_YDEFL) ++N;
  return(N);
}


G4bool whatFixedOk(char *rfname, G4int whatisfixed)
{
  // given the bit-packed list of fixed quantities,
  // return whether it is sufficient and unambiguous - 
  // actually causes abort upon failure
  G4bool ok,okOff,ok1out,okPh,okV,noOut,okIncr;

  const char *where= rfname!=NULL ? rfname : "whatFixedOk";

  //look at fixed quantities
  int Nfixed= howmanyfixed(whatisfixed);
  
  ok1out= Nfixed==1;
  noOut= Nfixed==0;
  okPh= whatisfixed & ATFIX_PHASE;
  okV= whatisfixed & ATFIX_GRADIENT;

  okIncr= whatisfixed & ATFIX_OFFSETINCR;
  okOff= whatisfixed & ATFIX_OFFSET;

  ok= !(okOff && okIncr);
  if(!ok)
    {
      G4Exception(where,"Cannot set both timeOffset and timeIncrement",FatalException,"");  //abort
      return(ok);
    }

  //know time & voltage or require exactly 2 of 3 
  ok= okOff && okV || okV && okPh && noOut || ok1out && okV && !okPh || ok1out && !okV && okPh;

  if(ok)  return(ok);

  if(!ok1out && !noOut)
    {
      printf("rfdevice(%s):whatFixedOk: ",where);
      if(whatisfixed & ATFIX_P)  printf("fixP=Y ");
      if(whatisfixed & ATFIX_DE)  printf("fixDE=Y ");
      if(whatisfixed & ATFIX_DT)  printf("fixD& ");
      if(whatisfixed & ATFIX_XDEFL)  printf("fixXdeflection=Y ");
      if(whatisfixed & ATFIX_YDEFL)  printf("fixYdeflection=Y ");
      printf("- only 1 allowed!\n");
      G4Exception(where,"Maximum of 1 fixed output setting allowed",FatalException,"");  //abort
      return(ok);
    }
  else if(noOut)
    {
      printf("rfdevice(%s):whatFixedOk: ",where);
      printf("maxGradient=%s ", whatisfixed & ATFIX_GRADIENT? "Y": "N");
      printf("phaseAcc=%s ", whatisfixed & ATFIX_DE? "Y": "N");
      printf("- require both iff no output settings requested.\n");
      G4Exception(where,"Insufficient info to specify timing",FatalException,"");  //abort
      return(ok);
    }
  else
    {
      if(okPh && okV)
	G4Exception(where,"Cannot set both phase and gradient with fixed output setting",FatalException,"");
      else
	G4Exception(where,"Require phase xor gradient with fixed output setting",FatalException,"");
    }

  return(ok);
}


AutoTimingSnapshot step2snapshot(const G4Step *step, G4double toff, G4double L, G4RotationMatrix *rot, G4bool within)
{
  AutoTimingSnapshot snapshot;

  G4Track *trk = step->GetTrack();

  G4StepPoint *pnt = step->GetPostStepPoint();

  const G4ThreeVector& globalxyz= pnt->GetPosition();

  const G4ThreeVector& globalPxyz= pnt->GetMomentum();

  G4ThreeVector localPxyz= globalPxyz;

  if(rot) 
    localPxyz= *rot * localPxyz;

  snapshot.x= globalxyz.getX();
  snapshot.y= globalxyz.getY();
  snapshot.z= globalxyz.getZ();

  snapshot.Px= localPxyz.getX();
  snapshot.Py= localPxyz.getY();
  snapshot.Pz= localPxyz.getZ();
  snapshot.Ptot= sqrt(snapshot.Px*snapshot.Px+snapshot.Py*snapshot.Py+snapshot.Pz*snapshot.Pz);

  snapshot.m= trk->GetDefinition()->GetPDGMass();
  snapshot.t= pnt->GetGlobalTime();
  snapshot.KE= trk->GetKineticEnergy();
  snapshot.v= trk->GetVelocity();
  snapshot.q= trk->GetDefinition()->GetPDGCharge();

  snapshot.timeOff= toff;
  snapshot.estTransit= L/trk->GetVelocity(); /*estimated transit time over timing volume*/
  snapshot.within= within;

  return(snapshot);
}


AutoTimingMethod autotiming2method(G4String autotimingpreference)
{
  //KBB convert string into enumerated type to control the automatic RF tuning

        if(!strcasecmp(autotimingpreference,"atZlocal") || !strcasecmp(autotimingpreference,"atZ"))
	  return(AUTOTIMING_ATZ);                            //drift to center (usually)
	else if(!strcasecmp(autotimingpreference,"maxEnergyGain") || !strcasecmp(autotimingpreference,"maxE"))
	  return(AUTOTIMING_MAXENERGY);                           //maximum dE
	else if(!strcasecmp(autotimingpreference,"noEnergyGain") || !strcasecmp(autotimingpreference,"noE"))
	  return(AUTOTIMING_NOENERGY);                            //dE=0
	else if(!strcasecmp(autotimingpreference,"minEnergyGain") || !strcasecmp(autotimingpreference,"minE")) 
	  return(AUTOTIMING_MINENERGY);                           //minimum dE
	else if(!strcasecmp(autotimingpreference,"minTransitTime") || !strcasecmp(autotimingpreference,"minT"))
	  return(AUTOTIMING_MINTIME);                            //minimum dT
	else if(!strcasecmp(autotimingpreference,"nomTransitTime") || !strcasecmp(autotimingpreference,"noT"))
	  return(AUTOTIMING_NOMTIME);                            //drift dT
	else if(!strcasecmp(autotimingpreference,"maxTransitTime") || !strcasecmp(autotimingpreference,"maxT"))
	  return(AUTOTIMING_MAXTIME);                            //maximum dT
	else if(!strcasecmp(autotimingpreference,"maxXdeflection") || !strcasecmp(autotimingpreference,"maxX"))
	  return(AUTOTIMING_MAXDXP);                             //max local delta-x' change
	else if(!strcasecmp(autotimingpreference,"minXdeflection") || !strcasecmp(autotimingpreference,"minX"))
	  return(AUTOTIMING_MINDXP);                             //min local delta-x' 
	else if(!strcasecmp(autotimingpreference,"noXdeflection") || !strcasecmp(autotimingpreference,"noX"))
	  return(AUTOTIMING_NODXP);                             //no local delta-x'
	else if(!strcasecmp(autotimingpreference,"maxYdeflection") || !strcasecmp(autotimingpreference,"maxY"))
	  return(AUTOTIMING_MAXDYP);                             //max local delta-y' change
	else if(!strcasecmp(autotimingpreference,"minYdeflection") || !strcasecmp(autotimingpreference,"minY"))
	  return(AUTOTIMING_MINDYP);                             //min local delta-y' change
	else if(!strcasecmp(autotimingpreference,"noYdeflection") || !strcasecmp(autotimingpreference,"noY"))
	  return(AUTOTIMING_NODYP);                              //no local delta-y' change

	G4Exception("autotiming2method:","Unsupported method",FatalException,"");  //abort

	return(AUTOTIMING_UNKNOWN);
}


const char *automethod2text(AutoTimingMethod method)
{
  //KBB convert enumerated type to a string

 switch(method) 
   {
   case AUTOTIMING_OFF: 
     return("off");
   case AUTOTIMING_ATZ: 
     return("atZlocal");
   case AUTOTIMING_MINENERGY: 
     return("minEnergyGain");
   case AUTOTIMING_NOENERGY: 
     return("noEnergyGain");
   case AUTOTIMING_MAXENERGY: 
     return("maxEnergyGain");
   case AUTOTIMING_MINTIME: 
     return("minTransitTime");
   case AUTOTIMING_NOMTIME: 
     return("nomTransitTime");
   case AUTOTIMING_MAXTIME: 
     return("maxTransitTime");
   case AUTOTIMING_MAXDXP:
     return("maxXdeflection");
   case AUTOTIMING_NODXP:
     return("noXdeflection");
   case AUTOTIMING_MINDXP:
     return("minXdeflection");
   case AUTOTIMING_MAXDYP:
     return("maxYdeflection");
   case AUTOTIMING_NODYP:
     return("noYdeflection");
   case AUTOTIMING_MINDYP:
     return("minYdeflection");
   case AUTOTIMING_UNKNOWN:
     return("unknown");
   }

  return("unknown?");
}


G4double automethod2timingphase(AutoTimingMethod method, G4double q)
{
  // returns the g4bl-style phase angle corresponding
  // to the tuning method and particle charge
  G4bool pos= q>=0;


  switch(method)
    {
    case AUTOTIMING_OFF: 
      return(0*deg);
    case AUTOTIMING_ATZ:
      return(90*deg);
    case AUTOTIMING_MAXENERGY: 
    case AUTOTIMING_MINTIME: 
    case AUTOTIMING_MAXDXP:
    case AUTOTIMING_MAXDYP:
      return(pos? 90*deg : -90*deg);
    case AUTOTIMING_MINENERGY: 
    case AUTOTIMING_MAXTIME: 
    case AUTOTIMING_MINDXP: 
    case AUTOTIMING_MINDYP: 
      return(pos? -90*deg : 90*deg);
    case AUTOTIMING_NOENERGY: 
    case AUTOTIMING_NOMTIME: 
    case AUTOTIMING_NODXP: 
    case AUTOTIMING_NODYP: 
    case AUTOTIMING_UNKNOWN:
      return(pos? 0*deg : 180*deg);
    default:
      return(0*deg);
    }
  
  return(0*deg);
}

class BLCMDrfdevice; // forward reference

/**	RFdeviceField represents one placement of an rfdevice.
 *
 *	frequency=0.0 is accepted and yields an rfdevice with a constant
 *	Ez = maxGradient (useful to verify units are correct).
 **/
class RFdeviceField : public BLElementField, public BLManager::SteppingAction {
	G4String name;
	BLCMDrfdevice *rfdevice;
	G4VPhysicalVolume *timingPV;        //KBB:* typically either 1/2 or full volume depending on method
	BLCoordinateTransform global2local;
	G4RotationMatrix *rotation;
                                            //KBB: timeOffset now within rfdevice
	G4Track saveTrack;
	G4int timeCount;
	G4bool validSaveTrack;
	BLFieldMap *fieldMap;
        AutoTimingMethod automethod;              //KBB:9mar11: autoTiming method
	friend class BLCMDrfdevice;

private:                              
        AutoTimingState stateOfPlacement;  //KBB:status this particular placement
        G4double maxGradient;
        G4double timeOffset;
        G4int redo;

        char pbname[256];                  //name of rfdevice
        int timeCount000;                  //starting count of revisited timing
        G4bool revisited_timingZero;       //has timingZero been re-evalued yet?
        double timingZero;                 //working time offset for 0 degRF
        double timingPhase;                //current phase for timing
        double timingZeroStep;             //current step for adjusting timingZero
        double timingPhaseStep;            //current step for adjusting phaseAcc
        double timingGrad;                 //current gradient for timing
        double timingGradStep;             //current step for adjusting gradient
        AutoTimingChange previous;         //previous value of a change
        int cycle;                         //counter for adjustment cycles
        double ampl;                       //working output amplitude 
        double drift;                      //drift time
        AutoTimingSnapshot entrance,exit;  //misc. info taken at a point in time & space
        AutoTimingSnapshot base,estimate;  //misc. info FYI
  

public:
	/// constructor. _zmin/max are global coordinates.
	RFdeviceField(G4String& _name,BLCoordinateTransform& _global2local,
			BLCMDrfdevice *_rfdevice, G4VPhysicalVolume *_timingPV);

	/// getName() returns the name of this placement.
	G4String getName() { return name; }

	/// addFieldValue() adds the field for this rfdevice into field[].
	/// point[] is in global coordinates.
	void addFieldValue(const G4double point[4], G4double field[6]) const;

	/// doStep() handles a step of the timing particle. 
	void UserSteppingAction(const G4Step *step);

        /// show1step() prints handy messages about autoTiming
        void show1step(AutoTimingState state,const char *pre,G4double wrkval,const char *suf);

        /// reset() resets various internal used by autoTime to initial settings
        void reset();
};


/**	BLCMDrfdevice implements an rfdevice.
 *
 *	Each placement of a BLCMDrfdevice creates an RFdeviceField that is linked
 *	into BLGlobalField. Each RFdeviceField sets its rfdevice->timeOffset via the
 *	timing particle (re-stepping through its TimingVol if necessary).
 *
 *	If fieldMapFile is non-null, it is read as a BLFieldMap and determines
 *	the peak field; both B and E are multiplied by maxGradient of the 
 *	cavity.
 **/
class BLCMDrfdevice : public BLElement {
	G4double maxGradient;
	G4double phaseAcc;

	G4String color;
	G4double frequency;
	G4double innerLength;
	G4double innerRadius;
	G4double pipeThick;
	G4double wallThick;
	G4String wallMat;
	G4double irisRadius;
	G4double collarRadialThick;
	G4double collarThick;
	G4double win1Thick;
	G4double win1OuterRadius;
	G4double win2Thick;
	G4String winMat;
	G4double skinDepth;
	G4double timingTolerance;
	G4double maxStep;
	G4String cavityMaterial;
	G4double timeIncrement;
	G4String fieldMapFile;

        G4double timeOffset;                //KBB: associated with rfdevice rather than field

        G4String autoTimingMethod;          //KBB: alternative tuning methods
        G4double autoTimingFixP;            //KBB: requested fixed |P| output 
        G4double autoTimingFixDE;           //KBB: requested fixed dKE output
        G4double autoTimingFixDT;           //KBB: requested fixed transit time output
        G4double autoTimingFixXdeflection;  //KBB: requested fixed local XZ angle
        G4double autoTimingFixYdeflection;  //KBB: requested fixed local YZ angle
        G4double autoTimingFailureLimit;    //KBB: tolerable difference betwen goal and acceptable
        G4double autoTimingZlocal;          //KBB: local Z location for method=atZlocal 
        G4int verbose;            //KBB: display timingVolume and timing info

	G4int kill;

	// non-arguments:
        AutoTimingState state;              //current working state
        G4int whatFixed;                    //bits specifing fixed input quantities 
                                            //not to be updated during calculations since
                                            //may have to start over using external tune
        AutoTimingWhatToAdjust what2adjust; //goal of adjustment - do not update

        G4double LtimingVolume;
	G4double overallLength;
	G4double outerRadius;
	G4double omega;
	G4double rFactor;
	G4double E2Bfactor;
	BLFieldMap *fieldMap;
	std::vector<RFdeviceField *> rfdeviceField;
	friend class RFdeviceField;

private:
        int changeCount;                   //counter for argChange

public:
	/// Default constructor. Defines the command, etc.
	BLCMDrfdevice();

	/// Destructor.
	virtual ~BLCMDrfdevice() { }

	/// Copy constructor.
	BLCMDrfdevice(const BLCMDrfdevice& r);

	/// clone()
	BLElement *clone() { return new BLCMDrfdevice(*this); }

	/// commandName() returns "rfdevice".
	G4String commandName() { return "rfdevice"; }

	/// command() implements the rfdevice command.
	int command(BLArgumentVector& argv, BLArgumentMap& namedArgs);

	/// defineNamedArgs() defines the named arguments of the command.
	void defineNamedArgs();

	/// argChanged() does internal computations after some arg changed
	void argChanged();

	/// construct() - construct the rfdevice.
	/// Creates a new RFdeviceField and adds it to BLGlobalField.
	void construct(G4RotationMatrix *relativeRotation,
			G4ThreeVector relativePosition, 
			G4LogicalVolume* parent,
			G4String parentName,
			G4RotationMatrix *parentRotation,
			G4ThreeVector parentPosition);

	/// getLength() returns the overallLength of the rfdevice
	G4double getLength() { return overallLength; }

	/// getWidth() returns the outer radius of the rfdevice
	G4double getWidth() { return outerRadius*2.0; }

	/// getHeight() returns the outer radius of the rfdevice
	G4double getHeight() { return outerRadius*2.0; }

	/// getSurveyPoint() returns points in LOCAL coordinates.
	G4ThreeVector getSurveyPoint(int index) {
		if(index == 0) return G4ThreeVector(0.0,0.0,-getLength()/2.0);
		if(index == 1) return G4ThreeVector(0.0,0.0,getLength()/2.0);
		throw "UNIMPLEMENTED";
	}

	/// isOK() returns true only if all placed rfdevices of this
	/// type have had their timing set, and if the tuning (if
	/// used) converged.
	G4bool isOK();

	/// trackOfInterest() returns true only if rfdevice timing just completed
	G4bool trackOfInterest();

	/// generatePoints() from BLElement
	void generatePoints(int npoints, std::vector<G4ThreeVector> &v);

	/// isOutside() from BLElement
	G4bool isOutside(G4ThreeVector &local, G4double tolerance);
};


BLCMDrfdevice defaultRFdevice;	// default object


// Default constructor - be sure to use the default constructor BLElement()
BLCMDrfdevice::BLCMDrfdevice() : BLElement(), rfdeviceField()
{
	// register the commandName(), and its synopsis and description.
	registerCommand(BLCMDTYPE_ELEMENT);
	setSynopsis("Defines an rfdevice (RF cavity)");
	setDescription("An rfdevice (RF cavity) is the basic RF element used to "
		"construct a linac.  The G4beamline convention is that 0degRF is the positive "
                "going zero crossing of the electric field, so generally phaseAcc=90 (degRF) is on-crest.\n\n"

                "The timeOffset parameter, if set, fixes the overall global absolute timing of the cavity "
                "relative to time=0 of the simulation.  If unspecified, 0degRF is determined via the "
		"timingMethod setting.  The default, timingMethod=atZlocal, defines 0degRF such that the "
                "test particle will arrive at the timingAtZ=# location then.\n\n"

                "For longitudinal cavities, timingMethod=maxEnergyGain emulates how most cavities "
                "in linacs have their overall timing determined; " 
		"while maxX would be appropriate "
                "for a horizontal transverse deflecting cavity.\n\n"

                "Independent of how 0degRF is found, exactly two of the set of maxGradient, phaseAcc, "
                "and one fixed output quantity (fixMomentum, fixEnergyGain, fixTransitTime, fixXdeflection, "
                "or fixYdeflection) must be specified to deterimine the final rfdevice timing.  For example, "
		"with maxGradient and phaseAcc set, the energy gain would be determined, while if maxGradient "
	        "and fixEnergyGain were set, the phaseAcc would be determined.\n\n"
 		
		"The pipe, walls, and collars are made of copper by default. "
		"Pipe, wall, collar, win1, and win2 may be omitted by setting their thickness to 0. "
		"Common usage is to set the collar values such that, by placing multiple rfdevices sequentially, "
		"the collars form a beam pipe between them.\n\n"

		"Note that section 4.4 of the User's Guide has a dimensioned drawing of a pillbox. "
                "Due to the presence of an (usually) invisible timing volume, care must be taken when placing "
                "objects within an rfdevice.\n\n"
		"This element must be placed (via the place command), and "
		"children can be placed inside it.\n\n"
		"See the User's Guide for details on how to use this complex "
		"element.");

	// provide initial values for fields
	timeOffset = undefined();                //[ns] (includes phase)
	maxGradient = undefined();               //[MV/m] (peak field)
	phaseAcc = undefined();                  //[degRF] (90degRF==oncrest)

	color = "1.0,0.0,0.0";
	frequency = undefined();
	innerLength = undefined();
	innerRadius = undefined();
	pipeThick = 3.0*mm;
	wallThick = 3.0*mm;
	irisRadius = 11.0*cm;
	collarRadialThick = 5.0*mm;
	collarThick = 2.5*mm;
	win1Thick = 0.200*mm;
	win1OuterRadius = 5.0*cm;
	win2Thick = 0.500*mm;
	winMat = "Be";
	wallMat= "Cu";                           //KBB-15mar11: allow for other cavity metals
	skinDepth = 0.002*mm;
	timingTolerance = 0.001*ns;
	maxStep = -1.0;
	cavityMaterial = "Vacuum";


	autoTimingMethod= "atZlocal";             //KBB-9mar11: alternative tuning methods
        autoTimingZlocal= undefined();           //KBB  center of rfdevice in local coordinates

        autoTimingFixP = undefined();            //KBB: requested fixed |P| output 
        autoTimingFixDE = undefined();           //KBB: requested fixed dKE output
        autoTimingFixDT = undefined();           //KBB: requested fixed dT output
        autoTimingFixXdeflection = undefined();  //KBB: requested fixed local XZ angle
        autoTimingFixYdeflection = undefined();  //KBB: requested fixed local YZ angle

        autoTimingFailureLimit = 1.E-3;           //KBB: tolerable error between fixed setting and result
        verbose = 1;                    //KBB: display timing volume and print info

        timeIncrement = 0;
	fieldMapFile = "";
	kill = 0;
 	// non-arguments:
	state= ATWORKING_UNKNOWN;               //may need reset later
        whatFixed= 0;                           //nothing yet
        what2adjust = AUTOADJUST_UNKNOWN;       //unknown yet
	overallLength = undefined();
	outerRadius = undefined();
	LtimingVolume= undefined();
	omega = undefined();
	rFactor = undefined();
	E2Bfactor = undefined();
	fieldMap = 0;
	changeCount= 0;
}

// Copy constructor - be sure to use the copy constructor BLElement(r)
BLCMDrfdevice::BLCMDrfdevice(const BLCMDrfdevice& r) : BLElement(r), rfdeviceField()
{
        G4int whatFixedNow= r.whatFixed;

	// copy fields one at a time (transfers default values from the
	// default object to this new object).

        //allow maxGradient to be externally tuned if not autoTiming
        if(whatFixedNow & ATFIX_GRADIENT && !isUndefined(r.maxGradient)) 
	  {
	    BLTune::copyTunableArg(&maxGradient,&r.maxGradient);
	  }
	else
	  {
	    whatFixedNow &= ~ATFIX_GRADIENT;
	    maxGradient= undefined();
	  }

	//allow timeOffset to be externally tuned if not autoTiming
	if(whatFixedNow & ATFIX_OFFSET && !isUndefined(r.timeOffset))
	  {
	    BLTune::copyTunableArg(&timeOffset,&r.timeOffset);
	  }
	else
	  {
	    whatFixedNow &= ~ATFIX_OFFSET;
	    timeOffset= undefined();
	  }

	//do not allow phaseAcc to be tuned, as it's within timeOffset
	if(whatFixedNow & ATFIX_PHASE)
	  phaseAcc = r.phaseAcc;
	else
	  phaseAcc = undefined();

	color = r.color;
	frequency = r.frequency;
	innerLength = r.innerLength;
	innerRadius = r.innerRadius;
	pipeThick = r.pipeThick;
	wallThick = r.wallThick;
	wallMat = r.wallMat;
	irisRadius = r.irisRadius;
	collarRadialThick = r.collarRadialThick;
	collarThick = r.collarThick;
	win1Thick = r.win1Thick;
	win1OuterRadius = r.win1OuterRadius;
	win2Thick = r.win2Thick;
	winMat = r.winMat;
	skinDepth = r.skinDepth;
	timingTolerance = r.timingTolerance;
	maxStep = r.maxStep;
	cavityMaterial = r.cavityMaterial;
	fieldMapFile = r.fieldMapFile;

        timeIncrement = r.timeIncrement;

        autoTimingMethod = r.autoTimingMethod;

	if(whatFixedNow & ATFIX_ZLOCAL)
	  autoTimingZlocal = r.autoTimingZlocal;
	else
	  autoTimingZlocal = undefined();

        autoTimingFixP = r.autoTimingFixP;
        autoTimingFixDE = r.autoTimingFixDE;
        autoTimingFixDT = r.autoTimingFixDT;
        autoTimingFixXdeflection = r.autoTimingFixXdeflection;
        autoTimingFixYdeflection = r.autoTimingFixYdeflection;

        autoTimingFailureLimit = r.autoTimingFailureLimit;
	verbose = r.verbose;

	kill = r.kill;

	// non-arguments:
	state= ATWORKING_UNKNOWN;     //copies need to be updated individually
        whatFixed = whatFixedNow;   //possibly updated
	what2adjust = r.what2adjust;
	LtimingVolume = r.LtimingVolume;
	overallLength = r.overallLength;
	outerRadius = r.outerRadius;
	omega = r.omega;
	rFactor = r.rFactor;
	E2Bfactor = r.E2Bfactor;
	fieldMap = r.fieldMap;
	changeCount= r.changeCount;
}


int BLCMDrfdevice::command(BLArgumentVector& argv, BLArgumentMap& namedArgs)
{
	// If the name of the cavity is "testUnits", perform the test and exit.
	if(argv[0] == "testUnits") {
		// set the test case values, from Al Moretti's SuperFish example
		maxGradient = 1.0 * megavolt/meter;
		frequency = 0.956188/ns;
		double Hcorrect = 1545.0; // A/m
		double BcorrectTesla = 4.0*pi*1E-7 * Hcorrect; // MKS
		argChanged();
		// the first maximum of J1() has a value of 0.5819
		double Bmax = maxGradient*E2Bfactor*0.5819;
		printf("pillobx testUnits:\n           "         
			"Freq(GHz)  Emax(MV/m)  Bmax(Tesla)  Radius(mm)\n");
		printf("Values:  %10.4f %10.4f %12.4e %11.4f\n",
			frequency/(1e9*hertz), maxGradient/(megavolt/meter),
			Bmax/tesla, innerRadius/mm);
		// these values already are in the correct units:
		printf("Correct: %10.4f %10.4f %12.4e %8.1f\n",
			0.956188,1.0,BcorrectTesla, 120.0);
		printf("rfdevice testUnits: exiting\n");
		g4bl_exit(0);
	}

	if(argv.size() != 1) {
		printError("rfdevice: Invalid command, must have name");
		return -1;
	}

	if(argv[0] == "default") {
		return defaultRFdevice.handleNamedArgs(namedArgs);
	}

	BLCMDrfdevice *p = new BLCMDrfdevice(defaultRFdevice);
	p->state= ATWORKING_UNKNOWN;
	p->setName(argv[0]);
	int retval = p->handleNamedArgs(namedArgs);

	if(p->fieldMapFile != "") {
		p->fieldMap = new BLFieldMap();
		p->fieldMap->readFile(p->fieldMapFile);
	}

	// check material exists
	if(p->cavityMaterial.size() > 0) getMaterial(p->cavityMaterial);

	p->print(argv[0]);

	return retval;
}


void BLCMDrfdevice::defineNamedArgs()
{
  	argTunable(maxGradient,"maxGradient","The peak gradient of the cavity (MV/m)",megavolt/meter);
	argString(color,"color","The color of the cavity");
	argDouble(frequency,"frequency","The frequency of the cavity (GHz)",1/ns,"",false);
	argDouble(innerLength,"innerLength","The inside length of the cavity (mm)",mm,"",false);
	argDouble(innerRadius,"innerRadius","The inside radius of the cavity (mm)",mm,"",false);
	argDouble(pipeThick,"pipeThick","The thickness of the pipe wall (mm)",mm,"",false);
	argDouble(wallThick,"wallThick","The thickness of the cavity walls (mm)",mm,"",false);
	argString(wallMat,"wallMat","The material of all the walls [Cu]");
	argDouble(irisRadius,"irisRadius","The radius of the iris (mm)",mm,"",false);
	argDouble(collarRadialThick,"collarRadialThick","The radial thickness of the collar (mm)",mm,"",false);
	argDouble(collarThick,"collarThick","The thickness of the collar along z(mm)",mm,"",false);
	argDouble(win1Thick,"win1Thick","The thickness of the central portion of the\n"
		"windows; zero for no window (mm)",mm,"",false);
	argDouble(win1OuterRadius,"win1OuterRadius","The radius of the central portion of\n"
		"the windows (mm)",mm,"",false);
	argDouble(win2Thick,"win2Thick","The thickness of the outer portion of the\n"
		"windows; zero for no window (mm)",mm,"",false);
	argString(winMat,"winMat","The material of the windows [Be].");
	argDouble(phaseAcc,"phaseAcc","The reference phase of the cavity (degrees).",deg);
	argDouble(skinDepth,"skinDepth","The skin depth (mm).",mm,"",false);
	argDouble(timingTolerance,"timingTolerance","Tolerance for timing tuning (ns)",ns);
	argDouble(maxStep,"maxStep","The maximum stepsize in the element (mm).",mm);
	argString(cavityMaterial,"cavityMaterial","Material of cavity volume [Vacuum].");

  	argTunable(timeOffset,"timeOffset","Time offset for cavity [set via timingMethod] (ns).",ns);
        argDouble(timeIncrement,"timeIncrement","Increment to timeOffset, applied AFTER tuning. (ns).",ns);

 	argString(autoTimingMethod,"timingMethod","Method for determining the nominal timeOffset "
		  "{atZ, maxE, noE, minE, maxT, nomT, minT, maxX, noX, minX, maxY, noY, minY}.");

	argDouble(autoTimingZlocal,"timingAtZ","Local Z location for timing (mm).",mm);

	argDouble(autoTimingFixP,"fixMomentum","Specify total output momentum (MeV/c).",MeV);
	argDouble(autoTimingFixDE,"fixEnergyGain","Specify energy gain (MeV).",MeV);
	argDouble(autoTimingFixDT,"fixTransitTime","Specify transit time (ns).",ns);
	argDouble(autoTimingFixXdeflection,"fixXdeflection","Specify local output XZ angle (deg).",deg);
	argDouble(autoTimingFixYdeflection,"fixYdeflection","Specify local output YZ angle (deg).",deg);
	argDouble(autoTimingFailureLimit,"fixTolerance","Specify allowable error on fixed settings [1.e-3].");

        argInt(verbose,"verbose","Set nonzero to show timing volume and print info messages [1].");

	argString(fieldMapFile,"fieldMapFile","Filename for BLFieldMap (pillbox if null).",false);
	argInt(kill,"kill","Set nonzero to kill tracks that hit the pipe, walls, or collars [0].");
}


void BLCMDrfdevice::argChanged()
{
        /* if called again, reinitialize */

	// compute non-argument variables
	omega = 2.0 * pi * frequency;
	if(isUndefined(innerRadius)) {
		if(frequency > 0.0) {
			// the first zero of J0(r) is at r=2.405
			innerRadius =  2.405 / (omega/c_light);
		} else {
			printError("rfdevice: innerRadius is undefined");
		}
	}
	double wall=collarThick;
	if(wallThick > wall) wall = wallThick;
	if(win1Thick > wall) wall = win1Thick;
	if(win2Thick > wall) wall = win2Thick;
	overallLength = innerLength + 2.0 * wall;
	outerRadius = innerRadius + pipeThick;
	rFactor = omega/c_light;
	// Previous versions had the incorrect value from the BeamTools:
	//     E2Bfactor = 1.0e-9 / c_light;   // MeV/m, MKS units
	// C.f. Eq. 8.77 in J.D.Jackson, _Classical_Electrodynamics_,
	//	 first edition (1962). His value is 1/c in vacuum.
	E2Bfactor = 1.0/c_light;   // Geant4 units
	if(maxStep < 0.0) maxStep = Param.getDouble("maxStep");


        G4bool printTimingSteps=verbose>1;  //print misc info

        // whether quantities were specified or not;
        // set starting guesses if required

	if(state==ATWORKING_RESET)  //reset all
	  {
	    if(!(whatFixed & ATFIX_OFFSET))   timeOffset= undefined();
	    if(!(whatFixed & ATFIX_GRADIENT))  maxGradient= undefined();
	    if(!(whatFixed & ATFIX_PHASE))  phaseAcc= undefined();
	    if(!(whatFixed & ATFIX_P))  autoTimingFixP= undefined();
	    if(!(whatFixed & ATFIX_DE))  autoTimingFixDE= undefined();
	    if(!(whatFixed & ATFIX_DT))  autoTimingFixDT= undefined();
	    if(!(whatFixed & ATFIX_XDEFL))  autoTimingFixXdeflection= undefined();
	    if(!(whatFixed & ATFIX_YDEFL))  autoTimingFixYdeflection= undefined();
            if(!(whatFixed & ATFIX_ZLOCAL))  autoTimingZlocal= undefined();
	  }

	state= ATWORKING_UNKNOWN;     
	what2adjust= AUTOADJUST_UNKNOWN;
	G4int prevFixed= whatFixed;
	whatFixed= 0;  //reset

	G4bool setOffset= !isUndefined(timeOffset);
	G4bool setOffsetIncr= timeIncrement!=0;
      	G4bool setGrad= !isUndefined(maxGradient);
        G4bool setPhase= !isUndefined(phaseAcc);
        G4bool setOff= !isUndefined(timeOffset);
        G4bool setP= !isUndefined(autoTimingFixP);
        G4bool setDE= !isUndefined(autoTimingFixDE);
        G4bool setDT= !isUndefined(autoTimingFixDT);
        G4bool setXd= !isUndefined(autoTimingFixXdeflection);
        G4bool setYd= !isUndefined(autoTimingFixYdeflection);
	G4bool setZ= !isUndefined(autoTimingZlocal);

	if(setOffset)  whatFixed |= ATFIX_OFFSET;
	if(setOffsetIncr)  whatFixed |= ATFIX_OFFSETINCR;
	if(setGrad)  whatFixed |= ATFIX_GRADIENT;
	if(setPhase)  whatFixed |= ATFIX_PHASE;
	if(setDE)  whatFixed |= ATFIX_DE;
	if(setP)  whatFixed |= ATFIX_P;
	if(setDT)  whatFixed |= ATFIX_DT;
	if(setXd)  whatFixed |= ATFIX_XDEFL;
	if(setYd)  whatFixed |= ATFIX_YDEFL;
	if(setZ)  whatFixed |= ATFIX_ZLOCAL;

        G4bool setFix= howmanyfixed(whatFixed) == 1;
	if(setFix)  whatFixed |= ATFIX_OUTPUT;  //ok output flag


	G4int justFixed= whatFixed & ~prevFixed;     //just newly fixed items

        int Njustfixed= howmanyfixed(justFixed);

	if(Njustfixed==1)  //retain only 1 fixed output - override any previous fixed outputs
	  { //          remove old output              insert new output     ok output flag
	    whatFixed= whatFixed & ~ATFIX_ANYOUT | justFixed & ATFIX_ANYOUT | ATFIX_OUTPUT; 
	  }
	else if(Njustfixed>1)
	  {
	    G4Exception("argChanged:","Maximum of 1 fixed output allowed",FatalException,"");  //abort
	  }

	if(setOffset && setGrad)
	  {
	    what2adjust= AUTOADJUST_UNNEEDED;
	  }
        else if(setGrad && setPhase && !setFix)
	  {
	    what2adjust= AUTOADJUST_OUTPUT;
	  }  
	else if(setGrad && !setPhase && setFix)
	  {
	    what2adjust= AUTOADJUST_PHASE;
	  }
	else if(!setGrad && setPhase && setFix)
	  {
	    what2adjust= AUTOADJUST_GRADIENT;
	  }
	else
	  {
	    what2adjust= AUTOADJUST_OUTPUT;
	  }
		
	if(printTimingSteps)
	  {
	    printf("BLCMDrfdevice::argChanged: changeCount=%d prevFixed=%d=%Xx justFixed=%d=%Xx ",
		   changeCount,prevFixed,prevFixed,justFixed,justFixed);

	    printf("whatFixed=%d=%Xx\n",whatFixed,whatFixed);

	    printf("setP=%d setDE=%d setDT=%d setXd=%d setYd=%d setGrad=%d ",
		   setP,setDE,setDT,setXd,setYd,setGrad);

	    printf("setPhase=%d setFix=%d setOff=%d\n",setPhase,setFix,setOff);
	    
	    printf("usa:timeOffset=%g ns. gradient=%g MV/m phase=%g degRF\n", 
		   timeOffset/ns, maxGradient/(megavolt/meter), phaseAcc/deg);
	    
	    fflush(stdout);
	  }

	/* Important to note that copies of rfdevice may be made later;
	   so whatFixed and what2adjust ought not to be updated */

	++changeCount;
}



void BLCMDrfdevice::construct(G4RotationMatrix *relativeRotation,
			G4ThreeVector relativePosition, 
			G4LogicalVolume* parent,
			G4String parentName,
			G4RotationMatrix *parentRotation,
			G4ThreeVector parentPosition)
{
  if(verbose>1)
    {
      printf("BLCMDrfdevice::construct...1\n");
      fflush(stdout);
    }

	G4String thisname = parentName+getName();

	G4Material *cu = getMaterial(wallMat);   //KBB 15mar11- no longer just copper
	G4Material *wm = getMaterial(winMat);
	G4Material *cavityMat = getMaterial(cavityMaterial);
	G4UserLimits *userLimits = new G4UserLimits(maxStep);

	userLimits->SetUserMinRange(0.0001*mm); //KBB 31may11 - minimum step size

	// find enclosing Tubs
	G4double halflen = innerLength/2+wallThick;
	halflen = std::max(halflen,innerLength/2+collarThick);
	halflen = std::max(halflen,innerLength/2+win1Thick);
	halflen = std::max(halflen,innerLength/2+win2Thick);
	G4Tubs *tubs = new G4Tubs(getName(),0.0,innerRadius+pipeThick,
					halflen,0.0,2.0*pi);
	G4LogicalVolume *rfdevice = new G4LogicalVolume(tubs,cavityMat,
					parentName+getName(),0,0,userLimits);

  if(verbose>1)
    {
      printf("BLCMDrfdevice::construct...2\n");
      fflush(stdout);
    }

	G4ThreeVector loc;
	G4VPhysicalVolume *pv;


	// The outer pipe
	G4LogicalVolume *pipe=0;
	if(pipeThick > 0.001*mm) {
		tubs = new G4Tubs(getName()+"pipe",innerRadius,
					innerRadius+pipeThick,
					innerLength/2+wallThick,0.0,2.0*pi);
		pipe = new G4LogicalVolume(tubs,cu,
					getName()+"pipeLog",0,0,userLimits);
		pv = new G4PVPlacement(0,loc,pipe,parentName+getName()+"Pipe",
					rfdevice,false,0,surfaceCheck);
		if(kill) BLManager::getObject()->
			registerSteppingAction(pv,new BLKillTrack(thisname));
	}

  if(verbose>1)
    {
      printf("BLCMDrfdevice::construct...2.1\n");
      fflush(stdout);
    }

	// the wall (us and ds) - radially: collarOuterRadius to innerRadius
	G4LogicalVolume *wall=0;
	if(wallThick > 0.001*mm) {
		tubs = new G4Tubs(getName()+"wall",irisRadius+collarRadialThick,
					innerRadius,wallThick/2,0.0,2.0*pi);
		wall = new G4LogicalVolume(tubs,cu,
					getName()+"wallLog",0,0,userLimits);
		G4double wallCenterZ = innerLength/2.0 + wallThick/2.0;
		loc.setZ(-wallCenterZ);
		pv = new G4PVPlacement(0,loc,wall,parentName+getName()+"UsWall",
					rfdevice,false,0,surfaceCheck);
		if(kill) BLManager::getObject()->
			registerSteppingAction(pv,new BLKillTrack(thisname));
		loc.setZ(+wallCenterZ);
		pv = new G4PVPlacement(0,loc,wall,parentName+getName()+"DsWall",
					rfdevice,false,0,surfaceCheck);
		if(kill) BLManager::getObject()->
			registerSteppingAction(pv,new BLKillTrack(thisname));
	}
	
	// the collar (us and ds) - includes wall thickness (i.e. the inner
	// edge of the collar is the same as the inner edge of the wall)
	G4LogicalVolume *collar=0;
	if(collarThick > 0.001*mm) {
		tubs = new G4Tubs(getName()+"collar",irisRadius,
					irisRadius+collarRadialThick,
					collarThick/2,0.0,2.0*pi);
		collar = new G4LogicalVolume(tubs,cu,
					getName()+"collarLog",0,0,userLimits);
		G4double collarCenterZ = innerLength/2.0 + collarThick/2.0;
		loc.setZ(-collarCenterZ);
		pv = new G4PVPlacement(0,loc,collar,parentName+getName()+"UsCollar",
					rfdevice,false,0,surfaceCheck);
		if(kill) BLManager::getObject()->
			registerSteppingAction(pv,new BLKillTrack(thisname));
		loc.setZ(+collarCenterZ);
		pv = new G4PVPlacement(0,loc,collar,parentName+getName()+"DsCollar",
					rfdevice,false,0,surfaceCheck);
		if(kill) BLManager::getObject()->
			registerSteppingAction(pv,new BLKillTrack(thisname));
	}
	
  if(verbose>1)
    {
      printf("BLCMDrfdevice::construct...2.2\n");
      fflush(stdout);
    }

	// win1 (the inner window)
	G4LogicalVolume *win1=0;
	if(win1Thick > 0.001*mm) {
		tubs = new G4Tubs(getName()+"win1",0.0,win1OuterRadius,
					win1Thick/2,0.0,2.0*pi);
		win1 = new G4LogicalVolume(tubs,wm,
					getName()+"win1Log",0,0,userLimits);
		G4double win1CenterZ = innerLength/2.0 + win1Thick/2.0;
		loc.setZ(-win1CenterZ);
		new G4PVPlacement(0,loc,win1,parentName+getName()+"UsWin1",
					rfdevice,false,0,surfaceCheck);
		loc.setZ(+win1CenterZ);
		new G4PVPlacement(0,loc,win1,parentName+getName()+"DsWin1",
					rfdevice,false,0,surfaceCheck);
	}

  if(verbose>1)
    {
      printf("BLCMDrfdevice::construct...2.3 win2Thick=%g mm\n",win2Thick/mm);
      fflush(stdout);
    }

	// win2 (the outer window)
	G4LogicalVolume *win2=0;
	if(win2Thick > 0.001*mm) {
		tubs = new G4Tubs(getName()+"win2",win1OuterRadius,irisRadius,
					win2Thick/2,0.0,2.0*pi);
		win2 = new G4LogicalVolume(tubs,wm,
					getName()+"win2Log",0,0,userLimits);
		G4double win2CenterZ = innerLength/2.0 + win2Thick/2.0;
		loc.setZ(-win2CenterZ);
		new G4PVPlacement(0,loc,win2,parentName+getName()+"UsWin2",
					rfdevice,false,0,surfaceCheck);
		loc.setZ(+win2CenterZ);
		new G4PVPlacement(0,loc,win2,parentName+getName()+"DsWin2",
					rfdevice,false,0,surfaceCheck);
	}

  if(verbose>1)
    {
      printf("BLCMDrfdevice::construct...3 autoTimingZlocal=");
      fflush(stdout);
      if(isUndefined(autoTimingZlocal))
	printf("undefined()\n");
      else
	printf("%g ns\n",autoTimingZlocal/ns);
    }

	G4double smidgen = 0.001*mm; //clearance

	AutoTimingMethod method= autotiming2method(autoTimingMethod);

	//generally atZlocal method uses Z located at center of cell
	//while other methods use near (a smidgen off) the exit of the rfdevice

	if(isUndefined(autoTimingZlocal))
	  autoTimingZlocal= method==AUTOTIMING_ATZ ? 0.0 : innerLength/2; 

  if(verbose>1)
    {
      printf("BLCMDrfdevice::construct...3.1 method=%s autoTimingZlocal=%g mm\n",
	     automethod2text(method), autoTimingZlocal/mm);

      fflush(stdout);
    }

	if(fabs(autoTimingZlocal) > innerLength/2)
	  {
	    printf(" autoTimingZlocal=%g mm\n",autoTimingZlocal/mm);
	    fflush(stdout);
	    G4Exception("rfdevice","timingAtZ outside of rfdevice!",FatalException,"");
	  }

	if(verbose>1)
	  {
	    printf("BLCMDrfdevice::construct...3.5 autoTimingZlocal=%g mm\n",autoTimingZlocal/mm);
	    fflush(stdout);
	  }

	// nominal length from entrance to exit

	G4double Ltub = innerLength/2 + autoTimingZlocal;

  if(verbose>1)
    {
      printf("BLCMDrfdevice::construct...4 autoTimingZlocal=%g mm, Ltub=%g mm\n",autoTimingZlocal/mm,Ltub/mm);
      fflush(stdout);
    }

	//KBB In general, the "rfdevice" could be placed entirely within an RF structure
        //KBB (such as a 9-cell telsa cavity) - the drawback is that outside that
        //KBB massless rfdevice the fields are zero...
        //KBB If the structure were placed within the rfdevice, but outside the timing
        //KBB volume, the fields would have their full extent.

	G4double tV_radius = innerRadius/2; //only use 1/2 rfdevice radius
	                                    //leave room for structures within rfdevice

	//reduce length from nominal just a little to avoid overlap issues
	LtimingVolume= Ltub-smidgen;

	tubs = new G4Tubs(getName()+"tiVo",0.0,tV_radius,LtimingVolume/2,0.0,2.0*pi);
	
	G4LogicalVolume *tiVo = new G4LogicalVolume(tubs,cavityMat,getName()+"tiVo",0,0,userLimits);

	loc.setZ(0.0);  //KBB - center of rfdevice

	loc.setZ((Ltub-innerLength)/2); //KBB - position upstream: begins just inside rfdevice entrance

	//timing volume used for timing the RF

	G4VPhysicalVolume *timingPV = new G4PVPlacement(0,loc,tiVo,
		parentName+getName()+"tiVo",rfdevice,false,0,surfaceCheck);

	if(verbose>1)
	  {
	    printf("BLCMDrfdevice::construct...5 innerLength=%g mm, LtimingVolume=%g mm\n",innerLength,LtimingVolume/mm);

	    printf("BLCMDrfdevice::construct...6 timeOffset=%g ns, phaseAcc=%g degRF, maxGradient=%g MV/m\n",
		   timeOffset/ns, phaseAcc/deg, maxGradient/(megavolt/meter));

	    fflush(stdout);
	  }

#ifdef G4BL_VISUAL
	G4VisAttributes *c1 = new G4VisAttributes(true,G4Color(0.5,0.0,0.0));
	G4VisAttributes *c2 = new G4VisAttributes(true,G4Color(0.7,0.0,0.0));
	G4VisAttributes *c3 = new G4VisAttributes(true,G4Color(1.0,0.0,0.0));
	const G4VisAttributes *invisible = BLCommand::getVisAttrib("Invisible");

	if(color)   //KBB override default colors
	  {
	    if(pipe) pipe->SetVisAttributes(getVisAttrib(color));
	    if(wall) wall->SetVisAttributes(getVisAttrib(color));
	    if(collar) collar->SetVisAttributes(getVisAttrib(color));
	    if(win1) win1->SetVisAttributes(getVisAttrib(color));
	    if(win2) win2->SetVisAttributes(getVisAttrib(color));
	  }
	else
	  {
	    if(pipe) pipe->SetVisAttributes(c1);
	    if(wall) wall->SetVisAttributes(c1);
	    if(collar) collar->SetVisAttributes(c1);
	    if(win1) win1->SetVisAttributes(c3);
	    if(win2) win2->SetVisAttributes(c2);
	  }

	//KBB - timing volume shown as ghost grey smoke or invisible
	G4VisAttributes *smoke = new G4VisAttributes(true,G4Color(0.6,0.5,0.5,0.45));

	tiVo->SetVisAttributes(verbose>0 ? smoke : invisible); 
	rfdevice->SetVisAttributes(invisible);
#endif

	// geant4 rotation convention is backwards from g4beamline
	G4RotationMatrix *g4rot = 0;
	if(relativeRotation)
		g4rot = new G4RotationMatrix(relativeRotation->inverse());

	// place the rfdevice into the parent
	new G4PVPlacement(g4rot,relativePosition,rfdevice,
			parentName+getName(),parent,false,0,surfaceCheck);

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

	BLCoordinateTransform global2local(globalRotation,globalPosition);

	// add this rfdevice field to the GlobalField, and to BLManager
	G4double dz = innerLength/2.0 + wallThick;
	// if rotated, make an overestimate of the field occupancy along z
	if(global2local.isRotated())
		dz += innerRadius;
	RFdeviceField *pf = 
		new RFdeviceField(thisname,global2local,this,timingPV);
	BLGlobalField::getObject()->addElementField(pf);
	rfdeviceField.push_back(pf);
	BLManager::getObject()->registerTuneParticleStep(timingPV,pf);

	printf("BLCMDrfdevice::construct %s parent=%s relZ=%.1f globZ=%.1f\n"
			"\tzmin=%.1f zmax=%.1f\n",
		thisname.c_str(),parentName.c_str(),relativePosition[2],
		global2local.getPosition()[2], global2local.getPosition()[2]-dz,
		global2local.getPosition()[2]+dz);
}


G4bool BLCMDrfdevice::isOK()
{
	G4bool retval = true;

	// verify all placed rfdeviceField-s have set their timeOffset 
	std::vector<RFdeviceField *>::iterator i;
	for(i=rfdeviceField.begin(); i!=rfdeviceField.end(); ++i) 
	  {
	    if((*i)->stateOfPlacement!=ATWORKING_DONE) 
	      {
		printf("*** BLCMDrfdevice %s {%s->%s} has not determined everything\n",
		       (*i)->getName().c_str(),
		       state2desc((*i)->stateOfPlacement), 
		       state2desc((*i)->rfdevice->state));

		retval = false;
	      }
	  }

	return retval;
}


G4bool BLCMDrfdevice::trackOfInterest()
{
  //a flag to indicate this tuning track of the last pass

       return TrackOfInterest;
}



RFdeviceField::RFdeviceField(G4String& _name, BLCoordinateTransform& _global2local,
			BLCMDrfdevice *_rfdevice,
			G4VPhysicalVolume *_timingPV) :
			BLElementField(), BLManager::SteppingAction()
{
	name = _name;
	stateOfPlacement= ATWORKING_UNKNOWN;
	maxGradient= BLCommand::undefined();
	timeOffset= BLCommand::undefined();
	redo= 0;

	rfdevice = _rfdevice;
	timingPV = _timingPV;
	global2local = _global2local;
	rotation = 0;

	//use rfdevicefield #s if available, then rfdevice #s, then safe values

	timeCount = 0;
	validSaveTrack = false;
	fieldMap = rfdevice->fieldMap;

	if(global2local.isRotated()) {
		rotation = new G4RotationMatrix(global2local.getRotation());
		rotation->invert();
	}

	// set global bounding box
	G4double local[4], global[4];
	local[3] = 0.0;
	G4double dz = rfdevice->innerLength/2.0 + rfdevice->wallThick;
	for(int i=0; i<2; ++i) {
		local[0] = (i==0 ? -1.0 : 1.0) * rfdevice->innerRadius;
		for(int j=0; j<2; ++j) {
			local[1] = (j==0 ? -1.0 : 1.0) * rfdevice->innerRadius;
			for(int k=0; k<2; ++k) {
				local[2] = (k==0 ? -1.0 : 1.0) * dz;
				global2local.getGlobal(local,global);
				setGlobalPoint(global);
			}
		}
	}

        sprintf(pbname,"%s",getName().c_str());
        int timeCount000= 0;         
        G4bool revisited_timingZero= false;   
        double timingZero= 0;           
        double timingPhase= 0;           
        double timingZeroStep= 0;         
        double timingPhaseStep= 0;   
        double timingGrad= 0;       
        double timingGradStep=0;      
        int cycle= 0;              
        double ampl= 0;               
        double drift= 0;                 
}


void RFdeviceField::reset()
{
	stateOfPlacement= ATWORKING_UNKNOWN;
	maxGradient= BLCommand::undefined();
	timeOffset= BLCommand::undefined();

	timeCount = 0;
	validSaveTrack = false;

        int timeCount000= 0;         
        G4bool revisited_timingZero= false;   
        double timingZero= 0;           
        double timingPhase= 0;           
        double timingZeroStep= 0;         
        double timingPhaseStep= 0;   
        double timingGrad= 0;       
        double timingGradStep=0;      
        int cycle= 0;              
        double ampl= 0;               
        double drift= 0;                 
}


void RFdeviceField::addFieldValue(const G4double point[4], G4double field[6]) const
{
	G4double local[4];
	global2local.getLocal(local,point);

	//use locked rfdevicefield values if available, 
	//then working rfdevice value, then a safe value

	G4double tOff = !BLCommand::isUndefined(timeOffset)? timeOffset : 
	  !BLCommand::isUndefined(rfdevice->timeOffset) ? rfdevice->timeOffset : 0;

	G4double Vo = !BLCommand::isUndefined(maxGradient)? maxGradient :
	  !BLCommand::isUndefined(rfdevice->maxGradient) ? rfdevice->maxGradient : TINY_NONZERO_GRADIENT;

	if(fieldMap) {
		G4double myField[6];

		//KBB-8mar11: this getFieldValue applies time & gradient!
		//
		//		fieldMap->getFieldValue(local,myField,rfdevice->maxGradient,
		//				rfdevice->maxGradient);

		fieldMap->getFieldValueNoTimeNoScale(local,myField);  //KBB-8mar11: w/o time, scale applied

		// frequency == 0.0 means a d.c. electric field
		G4double timeE=1.0, timeB=0.0;
		if(rfdevice->frequency != 0) {
		  double tArg = rfdevice->omega * (point[3] - tOff);
		  timeE = cos(tArg);
		  timeB = -sin(tArg);  //KBB 6/17/2011
		}

		timeE *= Vo/(megavolt/meter);
		timeB *= Vo/(megavolt/meter);

		if(rotation) {
			G4ThreeVector B(myField[0],myField[1],myField[2]);
			G4ThreeVector E(myField[3],myField[4],myField[5]);
			B = *rotation * B;
			E = *rotation * E;
			field[0] += B[0] * timeB;
			field[1] += B[1] * timeB;
			field[2] += B[2] * timeB;
			field[3] += E[0] * timeE;
			field[4] += E[1] * timeE;
			field[5] += E[2] * timeE;
		} else {
			field[0] += myField[0] * timeB;
			field[1] += myField[1] * timeB;
			field[2] += myField[2] * timeB;
			field[3] += myField[3] * timeE;
			field[4] += myField[4] * timeE;
			field[5] += myField[5] * timeE;
		}
		return;
	}

	G4double x = local[0];
	G4double y = local[1];
	G4double z = local[2];
	G4double r = sqrt(x*x + y*y);
	G4double dz = rfdevice->innerLength/2.0 + rfdevice->wallThick;
	if(z < -dz || z > dz || r > rfdevice->innerRadius) return;

	// frequency == 0.0 means a d.c. electric field
	if(rfdevice->frequency == 0.0) {
		// zmin and zmax include the walls, so re-test z position
		if(z >= -rfdevice->innerLength/2.0 && 
		   z <=  rfdevice->innerLength/2.0) {
			if(rotation) {
	 		        G4ThreeVector E(0.0, 0.0, Vo);
				E = *rotation * E;
				field[3] += E[0];
				field[4] += E[1];
				field[5] += E[2];
			} else {
			        field[5] += Vo;
			}
		}
		return;
	}

	// compute normal RF rfdevice field
	double arg = rfdevice->rFactor*r;
	if(fabs(arg) < 1.0e-20) arg = 0.0;
	double j0 = gsl_sf_bessel_J0(arg);
	double j1 = gsl_sf_bessel_J1(arg);
	double tArg = rfdevice->omega * (point[3] - tOff);
	double ez = Vo * j0 * cos(tArg);
	double bphi = -Vo * rfdevice->E2Bfactor * j1 * sin(tArg);
	// handle skin depth if in wall or window
	double f = 1.0;
	if(z > rfdevice->innerLength/2.0)
		f=exp(-(z-rfdevice->innerLength/2.0)/rfdevice->skinDepth);
	else if(z < -rfdevice->innerLength/2.0)
		f=exp(-(-z-rfdevice->innerLength/2.0)/rfdevice->skinDepth);
	//@ double phi = atan2(point[1],point[0]);
	//@ TJR 3/29/2011 - this surely should be local, not point
	double phi = atan2(local[1],local[0]);
	if(rotation) {
		G4ThreeVector B(-bphi*sin(phi)*f, bphi*cos(phi)*f, 0.0);
		G4ThreeVector E(0.0, 0.0, ez*f);
		B = *rotation * B;
		E = *rotation * E;
		field[0] += B[0];
		field[1] += B[1];
		field[2] += B[2];
		field[3] += E[0];
		field[4] += E[1];
		field[5] += E[2];
	} else {
		field[0] += -bphi * sin(phi) * f;
		field[1] += bphi * cos(phi) * f;
		field[5] += ez * f;
	}
}

void BLCMDrfdevice::generatePoints(int npoints, std::vector<G4ThreeVector> &v)
{
	generateTubs(npoints, 0.0, outerRadius, 0.0, 360.0*deg,
			overallLength, v);
}

G4bool BLCMDrfdevice::isOutside(G4ThreeVector &local, G4double tolerance)
{
	G4double r = sqrt(local[0]*local[0]+local[1]*local[1]);

	return r > outerRadius-tolerance ||
		fabs(local[2]) > overallLength/2.0-tolerance;
}



void RFdeviceField::show1step(AutoTimingState state,const char *pre,G4double wrkval,const char *suf)
{
  //show info on current autoTiming step

  if(!rfdevice->verbose)  return;

  int g= global_2TraceCounter;  //problems printing I*8

  printf("kbb_{%s}_%d:-2.%d:%s->%s ",pbname,timeCount,g,
	 autotimingstate2text(stateOfPlacement),autotimingstate2text(state));

  //alltrace numbers UserSteppingActions as Event#-2.1,-2.1000,-2.1001,...
  if(global_2TraceCounter<1000)
    global_2TraceCounter= 1000;
  else
    ++global_2TraceCounter; 

  G4int fx= rfdevice->whatFixed & ~ATFIX_OUTPUT;

  if(state==ATWORKING_UNKNOWN && fx!=0)
    {
      printf("fixed:");
      if(fx & ATFIX_OFFSET)  printf("%gnS",rfdevice->timeOffset/ns);
      fx&= ~ATFIX_OFFSET;

      if(fx & ATFIX_OFFSETINCR)
	printf("timeIncrement=%gnS",rfdevice->timeIncrement/ns);

      fx&= ~ATFIX_OFFSETINCR;

      if(fx!=0)  printf(",");
      if(fx & ATFIX_GRADIENT)  printf("%gMV/m",rfdevice->maxGradient/(megavolt/meter));
      fx&= ~ATFIX_GRADIENT;
      if(fx!=0)  printf(",");
      if(fx & ATFIX_PHASE)  printf("%gdegRF",rfdevice->phaseAcc/deg);
      fx&= ~ATFIX_PHASE;
      if(fx!=0)  printf(",");
      if(fx & ATFIX_DE)  printf("%gMeV",rfdevice->autoTimingFixDE/MeV);
      fx&= ~ATFIX_DE;
      if(fx!=0)  printf(",");
      if(fx & ATFIX_P)  printf("%gMeV/c",rfdevice->autoTimingFixP/MeV);
      fx&= ~ATFIX_P;
      if(fx!=0)  printf(",");
      if(fx & ATFIX_DT)  printf("%gns",rfdevice->autoTimingFixDT/ns);
      fx&= ~ATFIX_DT;
      if(fx!=0)  printf(",");
      if(fx & ATFIX_XDEFL)  printf("%gdegX",rfdevice->autoTimingFixXdeflection/deg);
      fx&= ~ATFIX_XDEFL;
      if(fx!=0)  printf(",");
      if(fx & ATFIX_YDEFL)  printf("%gdegY",rfdevice->autoTimingFixYdeflection/deg);
      fx&= ~ATFIX_YDEFL;
      if(fx!=0)  printf(",");
    }

  AutoTimingMethod mthd= autotiming2method( rfdevice->autoTimingMethod);
  printf("%s",automethod2text(mthd));
  if(!BLCommand::isUndefined(rfdevice->autoTimingZlocal))
    printf(",atZlocal=%gmm ",rfdevice->autoTimingZlocal);
  else
    printf(" ");

  printf("%s", 
	 rfdevice->what2adjust==AUTOADJUST_UNKNOWN ? "AUTOADJUST_UNKNOWN" :
	 rfdevice->what2adjust==AUTOADJUST_PHASE ? "AUTOADJUST_PHASE" :
	 rfdevice->what2adjust==AUTOADJUST_GRADIENT ? "AUTOADJUST_GRADIENT" :
	 rfdevice->what2adjust==AUTOADJUST_OUTPUT ? "AUTOADJUST_OUTPUT" : "?");

  if(rfdevice->what2adjust==AUTOADJUST_GRADIENT)
    {
      printf("=%gMV/m ",rfdevice->maxGradient/(megavolt/meter));
    }
  else if(rfdevice->what2adjust==AUTOADJUST_PHASE)
    {
      if(state==ATWORKING_ESTIMATE || state==ATWORKING_FINETUNE)
	{
	  printf("=%gdegRF ",timingPhase/deg);  //FYI only
	}
      else if(state==ATWORKING_FINAL || state==ATWORKING_DONE)
	{
	  printf("=%gdegRF ",rfdevice->phaseAcc/deg);  //FYI only
	}
      else if(state==ATWORKING_BASE)
	{
	  printf("=0degRF ");  //phase temporarially set to zero for base tracking
	}
      else
	{
	  printf(" ");
	}
    }
  else
    {
      printf(" ");
    }

  printf("%s%g%s\n",pre,wrkval,suf);
}
 

// called ONLY in timing particle mode, for this rfdevice's timingVol
void  RFdeviceField::UserSteppingAction(const G4Step *step)
{
  /*
    The traditional pillbox autophasing only worked for a 
    single cell cavity; here it's
    been generallized to any RF structure (via an specified RF map)
    and supports determining the timeOffset in a number of ways.

    In general, it won't matter much which method you specify to
    determine timeOffset - the AtZ is sets the arrival time for
    a local Z location (usually the center z=0), but isn't appropriate
    for all cases.  Usually people adjust the phase for maximum 
    energy gain and call that "on crest", 90degRF in the G4beamline
    notation.
    
    KBB March-April 2011
  */
  static int stepcall=0;

  TrackOfInterest= false;

  if(rfdevice->verbose>4)
    {
      printf("\nRFdeviceField::UserSteppingAction %s call#%d timeCount=%d %d=%s->rfdevice->state=%d=%s ", 
	     getName().c_str(), stepcall, timeCount, stateOfPlacement, state2desc(stateOfPlacement),
	     rfdevice->state, state2desc(rfdevice->state));

      printf("whatFixed=%d what2adjust=%d ",rfdevice->whatFixed,rfdevice->what2adjust);

      printf("frequency=%gGHz maxGradient=%gMV/m timeOffset=%gns\n",
	     rfdevice->frequency*ns, maxGradient/(megavolt/meter), timeOffset/ns);

      fflush(stdout);
    }

  /* 
     need to find out position w.r.t. rfdevice in case "tune" is in use 
     to update TrackOfInterest - otherwise could wait and do this later 
  */

  G4Track *track = step->GetTrack();
  G4StepPoint *prePoint = step->GetPreStepPoint();
  if(!prePoint) return;
  G4VPhysicalVolume *prePV = prePoint->GetPhysicalVolume();
  G4StepPoint *postPoint = step->GetPostStepPoint();
  if(!postPoint) return;
  G4VPhysicalVolume *postPV = postPoint->GetPhysicalVolume();
  G4SteppingManager *steppingMgr = BLManager::getObject()->getSteppingManager();

  enum WhereROI { ROI_OUTSIDE, ROI_ENTERING, ROI_INSIDE, ROI_EXITING };
  WhereROI whereROI;  
  
  if(prePV == timingPV && postPV == timingPV)
    whereROI= ROI_INSIDE;
  else if(prePV != timingPV && postPV == timingPV)
    whereROI= ROI_ENTERING;
  else if(prePV == timingPV && postPV != timingPV)
    whereROI= ROI_EXITING;
  else
    whereROI= ROI_OUTSIDE;

  G4bool Interesting=  prePV==timingPV || postPV==timingPV;

  // DC field specified
  if(rfdevice->frequency<=0)  
    {
      stateOfPlacement= rfdevice->state= ATWORKING_DONE;
      timeOffset= rfdevice->timeOffset= 0;
      maxGradient= rfdevice->maxGradient;
      TrackOfInterest= Interesting;
      return;
    }
 
  // zero field specified
  if(rfdevice->whatFixed & ATFIX_GRADIENT && rfdevice->maxGradient==0)  
    {
      stateOfPlacement= rfdevice->state= ATWORKING_DONE;
      maxGradient= rfdevice->maxGradient;
      timeOffset= rfdevice->timeOffset= 0;
      TrackOfInterest= Interesting;
      return;
    }

  // everything needed specified in advance
  if(rfdevice->what2adjust==AUTOADJUST_UNNEEDED)
    {
      stateOfPlacement= rfdevice->state= ATWORKING_DONE;
      maxGradient= rfdevice->maxGradient;
      timeOffset= rfdevice->timeOffset;
      if(rfdevice->whatFixed & ATFIX_OFFSETINCR)  timeOffset+= rfdevice->timeIncrement;
      TrackOfInterest= Interesting;
      return;
    }

  ++stepcall;
  
  //"tune" 2nd pass requires a re-autoTiming from beginning
  //(can improve efficiency later by comparing to previous entrance snapshot-KBB)
  if(whereROI==ROI_ENTERING)  //restart when reentering timingVolume - w/o "tune" only once
    {
      rfdevice->state= ATWORKING_RESET;
    }

  if(stateOfPlacement==ATWORKING_UNKNOWN || 
      rfdevice->state==ATWORKING_UNKNOWN || rfdevice->state==ATWORKING_RESET)  //clean start
    {
      if(stateOfPlacement==ATWORKING_DONE) ++redo;
      rfdevice->state= ATWORKING_RESET;  //force to initial state
      rfdevice->argChanged();
      rfdevice->state= ATWORKING_UNKNOWN;
      this->reset();
   }

  //no more work to do - rfdevice may be in use elsewhere
  if(stateOfPlacement==ATWORKING_DONE && !BLCommand::isUndefined(maxGradient) && !BLCommand::isUndefined(timeOffset)) 
    return;

  stateOfPlacement= ATWORKING_INFLUX;

  if(!timingPV)  G4Exception("rfdevice","Missing timing volume",FatalException,"");
  
  AutoTimingMethod method= autotiming2method(rfdevice->autoTimingMethod);
  
  if(prePV == postPV) return;     // neither entering nor leaving, do nothing

  if(rfdevice->verbose>2)
    {
      printf("\nRFdeviceField::UserSteppingAction %s call#%d redo#%d:timeCount=%d ", 
	     getName().c_str(), stepcall, redo, timeCount);

      printf("stateOfPlacement=%d=%s->rfdevice->state=%d=%s ", 
	     stateOfPlacement,state2desc(stateOfPlacement),
	     rfdevice->state,state2desc(rfdevice->state));

      printf("what2adjust=%d ",rfdevice->what2adjust);

      printf("frequency=%gGHz pf.maxGradient=%g p.maxGradient=%gMV/m pf.timeOffset=%g p.timeOffset=%gns\n",
	     rfdevice->frequency*ns, maxGradient/(megavolt/meter), rfdevice->maxGradient/(megavolt/meter), 
	     timeOffset/ns, rfdevice->timeOffset/ns);

      fflush(stdout);
    }

  //KBB-mar11 complete rethink of adjusting RF 
 
  const G4bool RevisitTimingAfterGradientAdjust= true;  //once correct gradient found, 
                                                        //refind timing zero using that gradient

  const G4bool TurnTransportationThresholdDownward= false;  //Geant4 thresholds of ~100-250MeV usually OK
                                                            //but in some rare case where particle is
                                                            //supposed to be trapped may be too high -
                                                            //does no harm to turn it on, except may prevent
                                                            //inhibit detection of trapped particles
  G4bool firstVisit=false,AnotherPass=false;
  G4int verbose= rfdevice->verbose;
  static G4bool negativeParticle=false;                
    
  AutoTimingChange chng;                  //change over timing volume

  double RFperiod= 1./rfdevice->frequency; //1 RFperiod in native time units
  double degRF= RFperiod/360.;            //1 degRF in native time units

  double dt=0,phi=0;

  if(rfdevice->state==ATWORKING_UNKNOWN)  //initialize all static info
    {
      if(verbose)
	{
	  printf("RFdeviceField::UserSteppingAction - %s - initialize ",pbname);
	  fflush(stdout);
	
	  printf("whatFixed=0x%X what2adjust=0x%X\n",rfdevice->whatFixed,rfdevice->what2adjust);
	}

      if(rfdevice->what2adjust==AUTOADJUST_UNKNOWN)  //abort
	G4Exception("RFdeviceField::UserSteppingAction","what2adjust==AUTOADJUST_UNKNOWN",FatalException,"");

      if(!whatFixedOk(pbname, rfdevice->whatFixed))
	G4Exception("RFdeviceField::UserSteppingAction","whatFixed not OK",FatalException,"");  //abort

      stepcall= 0;
      cycle= 1;
      timingZero=0; 
      timingPhase= BLCommand::undefined();
      timingZeroStep=0; 
      chng.dE= chng.Pout= chng.dT= chng.Xdefl= chng.Ydefl= 0;
      previous= chng;
      entrance.within= exit.within= false;
            
      firstVisit= true;

      if(rfdevice->state & ATFIX_OFFSET)
	{
	  timingZero= rfdevice->timeOffset;
	  rfdevice->state= ATWORKING_BYPASS;     //timeOffset known already; satisfy other conditions
	}
      else
	{
	  rfdevice->state= ATWORKING_OFFSET;     //start working on determining timeOffset
	  timingZero= rfdevice->timeOffset= 0;   //initial value
	}
    }
  else
    {
      firstVisit= false;
    }
  
  if(verbose)  // show progress for now
    {
      const G4ThreeVector& tmp_xyz= postPoint->GetPosition();
      G4double xyzt[4];
      xyzt[0]= tmp_xyz.getX();
      xyzt[1]= tmp_xyz.getY();
      xyzt[2]= tmp_xyz.getZ();
      xyzt[3]= postPoint->GetGlobalTime();
      
      const G4ThreeVector& tmp_mom= postPoint->GetMomentum();
      G4double Pxyz[3];
      Pxyz[0]= tmp_mom.getX();
      Pxyz[1]= tmp_mom.getY();
      Pxyz[2]= tmp_mom.getZ();
      
      G4double herefield[6];
      BLGlobalField::getObject()->GetFieldValue(xyzt,herefield);
      if(this->fieldMap)
	this->fieldMap->getFieldValue(xyzt,herefield);
      
      if(firstVisit)  //only print once 
	{
	  if(verbose>1)
	    {
	      printf("kbb:ROI:status(name) #x y z Px Py Pz t - - - - - Bx By Bz Ex Ey Ez count KE step\n");
	      
	      printf("kbb:method=%d=%s\n", method, automethod2text(method));

	      fflush(stdout);
	    }
	}
      
      if(verbose>1)
	{
	  printf("kbb:ROI:%s(%s) ", 
		 whereROI==ROI_INSIDE ? "inside" : 
		 whereROI==ROI_ENTERING ? "entering" :
		 whereROI==ROI_EXITING ? "EXITing" : "outside", pbname);
	  
	  printf("%g %g %g ",xyzt[0]/mm,xyzt[1]/mm,xyzt[2]/mm);
	  
	  printf("%g %g %g ",Pxyz[0]/MeV,Pxyz[1]/MeV,Pxyz[2]/MeV);
	  
	  printf("%g 0 0 0 0 0 ",xyzt[3]/ns);
	  
	  printf("%g %g %g ",herefield[0]/tesla,herefield[1]/tesla,herefield[2]/tesla);
	  
	  printf("%g %g %g ",herefield[3]/(megavolt/meter),
		 herefield[4]/(megavolt/meter),herefield[5]/(megavolt/meter));
	  
	  int g= stepcall;
	  printf("timeCount=%d %g stepcall=%d ",timeCount,track->GetKineticEnergy()/MeV,g);
	  
	  printf("%d=%s\n",rfdevice->state, state2desc(rfdevice->state));
	  
	  fflush(stdout);
	}
    }
  
  // Region-Of-Interest //

  switch(whereROI) 
    {
    case ROI_INSIDE:               //entirely inside - keep going
    case ROI_OUTSIDE:              //entirely outside - keep going
      return;

    case ROI_ENTERING:             //entering into timing volume - only once/rfdevice
      {
	G4RotationMatrix *rot= rotation;

	entrance= step2snapshot(step, rfdevice->timeOffset, rfdevice->LtimingVolume, rot, true);
	
	drift= entrance.estTransit;  //drift time
	
	timingPhase= automethod2timingphase(method, entrance.q);
	
	negativeParticle= entrance.q<0; 

	if(verbose>1 && firstVisit) 
	  {
	    printf("kbb------ %s ------- %s ",
		   pbname, automethod2text(method));
	    
	    printf("Gradient=%.3f MV/m",rfdevice->maxGradient/(megavolt/meter));
	    if(rfdevice->whatFixed&ATFIX_PHASE)
	      printf(", phaseAcc=%.1f degRF ",rfdevice->phaseAcc/deg);

	    printf(", timingPhase=%gdeg ",timingPhase/deg);

	    printf("\n");
	    fflush(stdout);
	  } 

	//In some rare future case a nearly trapped particle may be part of the plan...
	//but not in any of the current concepts.  This is primarially just FYI:
	if(firstVisit && TurnTransportationThresholdDownward)
	  {
	    // sometimes the track gets stuck after the rfdevice -
	    // it is unclear why - so tweak Process->Transportation settings

	    G4ProcessManager *procMgr = saveTrack.GetDefinition()->GetProcessManager();

	    G4ProcessVector *procList= procMgr->GetProcessList();
	    
	    for(int i=0; i<procList->size(); ++i) //look through whole list
	      {
		if(verbose>1)
		  printf("kbb::process#%d %s\n", i, (*procList)[i]->GetProcessName().data());
		
		if((*procList)[i]->GetProcessName()=="Transportation")
		  {
		    G4Transportation *trans= (G4Transportation *) (*procList)[i];
		    
		    G4double init_trans_ThresholdWarningEnergy= trans->GetThresholdWarningEnergy();
		    G4double init_trans_ThresholdImportantEnergy= trans->GetThresholdImportantEnergy();
		    G4int init_trans_ThresholdTrials= trans->GetThresholdTrials();
		    
		    // avoid Tune being killed for getting stuck or looping near boundry
		    
		    trans->SetThresholdWarningEnergy(0.1*keV);  //typically 100MeV
		    trans->SetThresholdImportantEnergy(0.250*keV);       //typically 250MeV
		    trans->SetThresholdTrials(3050);            //typically 10
		    
		    if(verbose>1)
		      {
			printf("kbb::trans->ThresholdWarningEnergy=%g=>%gMeV ", 
			       init_trans_ThresholdWarningEnergy/MeV,trans->GetThresholdWarningEnergy()/MeV);
			
			printf("ThresholdImportantEnergy=%g=>%gMeV ", 
			       init_trans_ThresholdImportantEnergy/MeV,trans->GetThresholdImportantEnergy()/MeV);
			
			printf("GetThresholdTrials=%d=>%d\n", 
			       init_trans_ThresholdTrials,trans->GetThresholdTrials());
		      }
		  }
		
		fflush(stdout);
	      }
	  }
	    
	// save track for setting up timeOffset, etc.
	// will return to this track later

	saveTrack.CopyTrackInfo(*step->GetTrack());
	saveTrack.SetUserInformation(0);
	validSaveTrack = true;
	
	if(verbose>1)  // show progress for now
	  {
	    printf("kbb-enter:%d timeOffset=%g\n",timeCount,rfdevice->timeOffset/ns);
	    fflush(stdout);
	  } 
	
	if(firstVisit)
	  {
	    if(rfdevice->state==ATWORKING_OFFSET)          // first visit to this rfdevice
	      {
		timingZero= rfdevice->timeOffset= 0;            // convenient settings
		timingZeroStep= 0*degRF;              // 0 degRF
	      }
	    else
	      {
		timingZero= rfdevice->timeOffset;
	      }
	    
	    if(rfdevice->what2adjust==AUTOADJUST_GRADIENT)  //estimate a reasonable starting value
	      {
		double DISCARD=ESSENTIALLY_UNCHANGED;
		double q= entrance.q;
		double c= c_light;
		double L= rfdevice->LtimingVolume;
		double to= entrance.estTransit;
		double m= entrance.m;
		double Po= entrance.Ptot;
		double Eo= sqrt(Po*Po+m*m);
		double w= !BLCommand::isUndefined(rfdevice->phaseAcc)? rfdevice->phaseAcc : pi/2;
		double dt= rfdevice->whatFixed & ATFIX_DT ? rfdevice->autoTimingFixDT - to : 0;
		double dE= rfdevice->whatFixed & ATFIX_DE ? rfdevice->autoTimingFixDE : 0;
		double Pf= rfdevice->whatFixed & ATFIX_P ? rfdevice->autoTimingFixP : Po;
		double xdf= rfdevice->whatFixed & ATFIX_XDEFL ? rfdevice->autoTimingFixXdeflection : 0;
		double ydf= rfdevice->whatFixed & ATFIX_YDEFL ? rfdevice->autoTimingFixYdeflection : 0;
		
		//avoid divide-by-zero issues; let it work on crest
		double sin_w= fabs(sin(w)) > DISCARD ? sin(w) : 1;
		double aw= fabs(ampl*sin(w)) > DISCARD ? ampl*sin(w) : 1;
		
		if(rfdevice->whatFixed & ATFIX_P)  //estimate dE from fixP & Po
		  {
		    dE= sqrt(Pf*Pf+m*m) - Eo;
		  }
		
		//simplest possible model as a starting point
		
		if(rfdevice->whatFixed & (ATFIX_P | ATFIX_DE)) 
		  {
		    double wavelength= 2*pi*c_light/rfdevice->omega;

		    if(2*L<wavelength)
		      timingGrad= dE*pi/(q*wavelength*sin(pi*L/wavelength)*sin_w);  //single cell?
		    else
		      timingGrad= dE*pi/(q*2*L*sin_w);      //multicell?
		  }
		else if(rfdevice->whatFixed & ATFIX_DT)        //total transit time
		  { 
		    timingGrad= (m/sqrt(1-1/pow(c*dt/L+Eo/Po,2.))-Eo)/(q*L*sin_w);
		  }
		else if(rfdevice->whatFixed & ATFIX_XDEFL)     //deflection angle in radians
		  {
		    timingGrad= (xdf/radian)*(Po/MeV)/q/(to/ns) * 
		      0.3*megavolt/meter/sin_w;
		  }
		else if(rfdevice->whatFixed & ATFIX_YDEFL)
		  {
		    timingGrad= (ydf/radian)*(Po/MeV)/q/(to/ns) *
		      0.3*megavolt/meter/sin_w;
		  }
		else if(lastFound_maxGradient!=0)
		  {
		    timingGrad= lastFound_maxGradient;
		  }
		else
		  {
		    timingGrad= 1*megavolt/meter;     //no idea - ought to never happen
		  }
		
		rfdevice->maxGradient= timingGrad;
		
		if(verbose>1)
		  {
		    printf("kbb:ROI:(%s) guess at working Gradient=%.3f MV/m\n", 
			   pbname, rfdevice->maxGradient/(megavolt/meter));
		  }
	      }
	  }	    

	return;                    //keep going on toward exit
      }
	
    case ROI_EXITING:
      break;                       //stop & process
    }
  
  // ROI_EXITING

  ++timeCount;

  exit= step2snapshot(step, rfdevice->timeOffset, rfdevice->LtimingVolume, rotation, true);
  chng= autoTimingIn2Out(entrance, exit);
  
  if(timeCount > ITERATION_LIMIT)
    G4Exception("rfdevice","Iteration Limit",FatalException,"");  //abort
  
  AnotherPass= true;
  
  switch(rfdevice->state)
    {
    case ATWORKING_OFFSET:
      {
	//determine basic timeOffset; RF ~ sin(w*(t-timeOffset))
	
	if(exit.within!=true)
	  G4Exception("rfdevice","Exit w/o Enter",FatalException,"");  //abort
	
	/*********************************************************************************
	       KBB: Assume that dE, dT, Xdeflection, Ydeflection approximately follow a sine-like curve - 
	       from 2 points 90 degrees apart in phase, estimate time offset for a maximum
	       f(wto)= K sin(w(t-to)+phi) = K sin(phi)
	       f(wto+RFperiod/4)= K sin(w(t-to)+phi+TT/2)=  K cos(w(t-to)+phi)= K cos(phi)
	       f(wto)/f(wto+RFperiod/4) = tan(phi) 
	       f(t=to) => sin(w(t-to)+phi)= 0 -> w(t-to)-phi=0 -> t-to=phi/w
	*********************************************************************************/
	
	if(timeCount<=1 || timeCount000+1==timeCount)  /*KBB always get 1st of 2 points to start */
	  {
	    //timeOffset and timingZero were 0 
	    
	    show1step(rfdevice->state,"(1of3)timingzero=",timingZero/ns,"ns");
	    
	    rfdevice->timeOffset= timingZero + pi/2/rfdevice->omega;  // was 0 -> now shift +90degRF
	    
	    if(verbose>1)
	      printf("kbb.1: was timeZero=timeOffset=0 => 90degRF\n");

	    rfdevice->state= ATWORKING_OFFSET; 
	  }
	else if(timeCount==2 || timeCount000+2==timeCount) /*KBB always from 2nd point, estimate overall phase */
	  {
	    show1step(rfdevice->state,"(2of3)timingzero=",timingZero/ns,"ns");
	    
	    //timeOffset was 90degRF from first
	    
	    if(method==AUTOTIMING_ATZ)      //KBB: based on TJR's original method of arrival time
	      {                               //at a local Z (timingVolume exit); usually middle of cell
	                                      // time of arrival   
		double t = track->GetGlobalTime();                  //time now (usually at the center)
				
		double dt= t - timingZero - timingPhase/rfdevice->omega; 
		
		while(dt>+RFperiod/2)  dt-= RFperiod;

		while(dt<-RFperiod/2)  dt+= RFperiod;

		timingZero+= dt;
		
		phi= atan2(previous.dE, chng.dE);      //kinetic energy change (for subsequent adjustments)

		if(verbose>1)
		  printf("kbb.2: dt=%gns timingZero->%gns\n",dt/ns,timingZero/ns);
	      }
	    else if(method==AUTOTIMING_MAXENERGY || method==AUTOTIMING_NOENERGY || method==AUTOTIMING_MINENERGY)
	      {
		phi= atan2(previous.dE, chng.dE);      //kinetic energy change  
		ampl= sqrt(previous.dE*previous.dE + chng.dE*chng.dE);
	      }
	    else if(method==AUTOTIMING_MINTIME || method==AUTOTIMING_NOMTIME || method==AUTOTIMING_MAXTIME)
	      {
		phi= atan2(previous.dT-drift, chng.dT-drift);        //(transit.time-drift) change
		ampl= sqrt((previous.dT-drift)*(previous.dT-drift) + (chng.dT-drift)*(chng.dT-drift));
	      }
	    else if(method==AUTOTIMING_MAXDXP || method==AUTOTIMING_NODXP || method==AUTOTIMING_MINDXP)
	      {
		phi= atan2(previous.Xdefl, chng.Xdefl);        //X deflection
		ampl= sqrt(previous.Xdefl*previous.Xdefl + chng.Xdefl*chng.Xdefl);
	      }
	    else if(method==AUTOTIMING_MAXDYP || method==AUTOTIMING_NODYP || method==AUTOTIMING_MINDYP)
	      {
		phi= atan2(previous.Ydefl, chng.Ydefl);        //Y deflection
		ampl= sqrt(previous.Ydefl*previous.Ydefl + chng.Ydefl*chng.Ydefl);
	      }
	    else
	      {
		phi= BLCommand::undefined();  //unknown method
		ampl= BLCommand::undefined();
	      }

	    if(rfdevice->whatFixed & (ATFIX_DE | ATFIX_P) || method==AUTOTIMING_ATZ)
	      {
		ampl= sqrt(previous.dE*previous.dE + chng.dE*chng.dE);
	      }
	    else if(rfdevice->whatFixed & ATFIX_DT )
	      {
		ampl= sqrt((previous.dT-drift)*(previous.dT-drift) + (chng.dT-drift)*(chng.dT-drift));
	      }
	    else if(rfdevice->whatFixed & ATFIX_XDEFL)
	      {
		ampl= sqrt(previous.Xdefl*previous.Xdefl + chng.Xdefl*chng.Xdefl);
	      }
	    else if(rfdevice->whatFixed & ATFIX_YDEFL)
	      {
		ampl= sqrt(previous.Ydefl*previous.Ydefl + chng.Ydefl*chng.Ydefl);
	      }
	    
	    if(verbose>1)
	      {
		printf("kbb.2: pass1-2 amplitude=");
		!BLCommand::isUndefined(ampl) ? printf("%g",ampl) : printf("undefined()");
		printf(" delta-phi-zero=");
		!BLCommand::isUndefined(phi) ? printf("%gdegRF",phi/deg) : printf("undefined()");
		printf(" timingZero=%gns rfdevice->timeOffset=%gns",timingZero/ns,rfdevice->timeOffset/ns);
		printf("\n");
	      }

	    if(method!=AUTOTIMING_ATZ && !BLCommand::isUndefined(phi))  //phi of the sine-like signal, not necc. overall phase
	      {
		//redefine timingZero to new value
		timingZero= phi/rfdevice->omega;   //zero crossing time w.r.t. previous 1st pass timeOffset=0

		if(negativeParticle)
		  timingZero+= pi/rfdevice->omega;  //negative particle shifted by 180deg

		rfdevice->timeOffset= timingZero - timingPhase/rfdevice->omega;

		if(verbose>1)
		  {
		    printf("kbb.2a: pass1-2 amplitude=");
		    !BLCommand::isUndefined(ampl) ? printf("%g",ampl) : printf("undefined()");
		    printf(" delta-phi-zero=");
		    !BLCommand::isUndefined(phi) ? printf("%gdegRF",phi/deg) : printf("undefined()");
		    printf(" timingZero=%gns rfdevice->timeOffset=%gns",timingZero/ns,rfdevice->timeOffset/ns);
		    printf("\n");
		  }
	      }

	    rfdevice->state= ATWORKING_OFFSET; 
	  }
	else    //all other passes...
	  {
	    if(method==AUTOTIMING_ATZ)      //KBB: based on TJR's original method of arrival time
	      {                               //at a local Z (timingVolume exit); usually middle of cell
	                                      // time of arrival   
		double t = track->GetGlobalTime();                  //time now (usually at the center)
		
		show1step(rfdevice->state,"timingzero=",timingZero/ns,"ns");
		
		double dt= t - (timingZero - timingPhase/rfdevice->omega); 
		
		while(dt>+RFperiod/2)  dt-= RFperiod;

		while(dt<-RFperiod/2)  dt+= RFperiod;

		timingZero+= dt;
		
		AnotherPass= fabs(dt) > rfdevice->timingTolerance;		
		
		if(verbose>1)
		  printf("kbb.%d: dt=%gns timingZero->%gns Another=%s\n",
			 timeCount,dt/ns,timingZero/ns,AnotherPass?"Y":"N");

	      }
	    else if(timeCount==3 || timeCount000+3==timeCount) //ought to have timingZero to be pretty close
	      {
		//ought to be parked on right point
		
		timingZeroStep= 10*degRF;     //overshoot a bit
		timingZero-= timingZeroStep;  //displace so begin working back toward estimated point
		
		show1step(rfdevice->state,"(3of3)timingZero=",timingZero/ns,"ns");
		  
		rfdevice->timeOffset= timingZero - timingPhase/rfdevice->omega;
		  
		if(verbose>1)
		  printf("kbb.3 timingZeroStep=%gdegRF timingZero=%gns\n",
			 timingZeroStep/degRF,timingZero/ns);

		AnotherPass= true;
	      }
	    else   //KBB - simple search about the method
	      {
		G4bool same=true;  //try#3,4 ought to be close

		switch(method)
		  {
		  case AUTOTIMING_MAXENERGY:
		    same= chng.dE > previous.dE;
		    break;
		  case AUTOTIMING_NOENERGY:
		    same= fabs(chng.dE) < fabs(previous.dE);
		    break;
		  case AUTOTIMING_MINENERGY:
		    same= chng.dE < previous.dE;
		    break;
		  case AUTOTIMING_MAXTIME:
		    same= chng.dT > previous.dT;
		    break;
		  case AUTOTIMING_NOMTIME:
		    same= fabs(chng.dT) < fabs(previous.dT);
		    break;
		  case AUTOTIMING_MINTIME:
		    same= chng.dT < previous.dT;
		    break;
		  case AUTOTIMING_NODXP:
		    same= fabs(chng.Xdefl) < fabs(previous.Xdefl);
		    break;
		  case AUTOTIMING_MAXDXP:
		    same= chng.Xdefl > previous.Xdefl;
		    break;
		  case AUTOTIMING_MINDXP:
		    same= chng.Xdefl < previous.Xdefl;
		    break;
		  case AUTOTIMING_NODYP:
		    same= fabs(chng.Ydefl) < fabs(previous.Ydefl);
		    break;
		  case AUTOTIMING_MAXDYP:
		    same= chng.Ydefl > previous.Ydefl;
		    break;
		  case AUTOTIMING_MINDYP:
		    same= chng.Ydefl < previous.Ydefl;
		    break;
		  default:
		    break;
		  }
		
		show1step(rfdevice->state,"timingzero=",timingZero/ns,"ns");

		AnotherPass= fabs(timingZeroStep) > rfdevice->timingTolerance;		
		
		if(AnotherPass)
		  {
		    if(!same)  timingZeroStep/= -1.5;     //backup and reverse direction, 2/3 step size
		    timingZero+= timingZeroStep;
		  }

		if(verbose>1)
		  printf("kbb.%d timingZeroStep=%gns Another=%s negativeParticle=%s\n",
			 timeCount,timingZeroStep/ns,AnotherPass?"y":"n",
			 negativeParticle?"Y":"N");
	      }

	    if(AnotherPass)
	      {
		rfdevice->timeOffset= timingZero - timingPhase/rfdevice->omega;  //adjustment point
		rfdevice->state= ATWORKING_OFFSET; //continue on
	      }
	    else
	      {
		/* 
		   the definition of 0degRF is the zero crossing on the rising slope,
		   independent of the sign of the particle 
		*/

		//bracket -RFperiod/2<timingZero<+RFperiod/2
		while(timingZero>RFperiod/2) timingZero-= RFperiod;
		while(timingZero<-RFperiod/2) timingZero+= RFperiod;

		rfdevice->timeOffset= timingZero; 

		rfdevice->state= ATWORKING_BASE;  //set next pass at 0degRF
	      }
	  }

	break;
      }

    case ATWORKING_BYPASS:   //base timeOffset set explicitly
      {
	show1step(rfdevice->state,"timingzero=",timingZero/ns,"ns");

	timingZero= rfdevice->timeOffset;  //time explicitly fixed
	timingPhase= 0;

	if(rfdevice->what2adjust==AUTOADJUST_GRADIENT)
	  {
	    rfdevice->state= ATWORKING_BASE;
	  }
	else  //nothing left to adjust
	  {
	    rfdevice->state= ATWORKING_FINAL;

	    if(rfdevice->whatFixed & ATFIX_OFFSETINCR)  
	      rfdevice->timeOffset= rfdevice->timeOffset+= rfdevice->timeIncrement;
	  }

	break;
      }

    case ATWORKING_BASE:  //ought to be final 0degRF baseline
      {
	/* 
	   Because 0degRF is DEFINED to be when voltage is at zero crosing and
	   rising with a postive slope, for both positive and negative reference
	   particles, the BASE ought to show a zero crossing with a rising slope. KBB
	*/

	show1step(rfdevice->state,"timingZero=",timingZero/ns,"ns");

	base= step2snapshot(step, rfdevice->timeOffset, rfdevice->LtimingVolume, rotation, true);
	
	AutoTimingChange diba= autoTimingIn2Out(entrance, base);
	
	//set approximate value to satisfy conditions
	
	switch(rfdevice->what2adjust)
	  {
	  case AUTOADJUST_OUTPUT:   //fixed gradient and phase
	    {
	      timingPhase= rfdevice->phaseAcc;
	      rfdevice->timeOffset= timingZero - rfdevice->phaseAcc/rfdevice->omega;

	      if(verbose>1)
		printf("kbb.%d BASE - adjust output phaseAcc=%gdegRF\n",timeCount,rfdevice->phaseAcc/deg);

	      rfdevice->state= ATWORKING_FINAL;  //done adjusting

	      if(rfdevice->whatFixed & ATFIX_OFFSETINCR)  
		rfdevice->timeOffset+= rfdevice->timeIncrement;

	      break;
	    }
	    
	  case AUTOADJUST_PHASE:    //fixed gradient and output
	    {
	      double signal=0;
	      
	      if(BLCommand::isUndefined(ampl))
		{
		  signal= 0;  //just start from base
		}
	      else if(rfdevice->whatFixed & ATFIX_P)  //estimate P swing from dE swing
		{
		  double Eo= sqrt(base.Ptot*base.Ptot+entrance.m*entrance.m);
		  ampl= sqrt((Eo+ampl)*(Eo+ampl)-entrance.m*entrance.m) - base.Ptot;
		  signal= (rfdevice->autoTimingFixP - base.Ptot)/ampl;
		}
	      else if(rfdevice->whatFixed & ATFIX_DE)  //base.DE ~ 0
		{
		  signal= (rfdevice->autoTimingFixDE - diba.dE)/ampl;
		}
	      else if(rfdevice->whatFixed & ATFIX_DT)  //total transit times
		{
		  signal= (rfdevice->autoTimingFixDT - diba.dT)/ampl;
		}
	      else if(rfdevice->whatFixed & ATFIX_XDEFL)  //base.Xdefl ~ 0
		{
		  signal= (rfdevice->autoTimingFixXdeflection - diba.Xdefl)/ampl;
		}
	      else if(rfdevice->whatFixed & ATFIX_YDEFL)  //base.Ydefl ~ 0
		{
		  signal= (rfdevice->autoTimingFixYdeflection - diba.Ydefl)/ampl;
		}
	      
	      if(signal>1.)
		signal= 1.;
	      else if(signal<-1.)
		signal= -1.;

	      timingPhase= asin(signal);

	      if(negativeParticle)  timingPhase+= pi;

	      if(verbose>1)
		printf("kbb.%d BASE - ampl=%g signal=%g timingPhase=%gdegRF\n",
		       timeCount,ampl,signal,timingPhase/deg);
	      
	      rfdevice->timeOffset= timingZero - timingPhase/rfdevice->omega;
	      
	      rfdevice->state= ATWORKING_ESTIMATE;  //best guess at phase
	      break;
	    }
	    
	  case AUTOADJUST_GRADIENT:    //fixed phase and output
	    {
	      //currently at phase 0, so no energy gain... reuse estimate

	      timingPhase= rfdevice->phaseAcc;

	      if(rfdevice->whatFixed & ATFIX_PHASE)
		{
		  rfdevice->timeOffset= timingZero - rfdevice->phaseAcc/rfdevice->omega;  //fixed phase & timing
		  timingPhase= rfdevice->phaseAcc;
		}
	      else if(rfdevice->whatFixed & ATFIX_OFFSET)
		{
		  timingZero= rfdevice->timeOffset;
		  timingPhase= 0;
		}
	      else   //ought to have phaseAcc or timeOffset set
		{
		  rfdevice->timeOffset= timingZero;
		  timingPhase= 0;
		}

	      timingGrad= rfdevice->maxGradient;
	      rfdevice->state= ATWORKING_ESTIMATE;  //best guess at gradient

	      if(verbose>1)
		{
		  printf("kbb:ROI:(%s) estimated maxGradient=%.6f MV/m\n", 
			 pbname, rfdevice->maxGradient/(megavolt/meter));
		}

	      break;
	    }

	  case AUTOADJUST_UNKNOWN:
	  case AUTOADJUST_UNNEEDED:
	    break;
	  }
	
	break;
      }
      
    case ATWORKING_ESTIMATE:     //extrapolated
      {
	estimate= step2snapshot(step, rfdevice->timeOffset, rfdevice->LtimingVolume, rotation, true);
	
	AutoTimingChange dies= autoTimingIn2Out(entrance, estimate);
	
	switch(rfdevice->what2adjust)
	  {
	  case AUTOADJUST_PHASE:    //fixed gradient and output
	    {
	      show1step(rfdevice->state,"timingZero=",timingZero/ns,"ns");
	      timingPhaseStep= 5*deg;         //back off a little
	      timingPhase-= timingPhaseStep;  //displace so begin working back toward estimated point	
	      rfdevice->timeOffset= timingZero - timingPhase/rfdevice->omega;
	      rfdevice->state= ATWORKING_FINETUNE;
	      rfdevice->phaseAcc= timingPhase; //FYI only - determined by timeOffset

	      if(verbose>1)
		{
		  printf("kbb.%d estimate (shift 10deg) timingPhaseStep=%gdegRF timingPhase=%gdegRF\n", 
			 timeCount, timingPhaseStep/deg, timingPhase/deg);
		}

	      break;
	    }
	    
	  case AUTOADJUST_GRADIENT:    //fixed phase and output
	    {
	      show1step(rfdevice->state,"timingGrad=",timingGrad/(megavolt/meter),"MV/m");
	      timingGradStep= timingGrad/3;  
	      timingGrad-= timingGradStep;  //displace so begin working back toward estimated point
	      rfdevice->maxGradient= timingGrad;
	      rfdevice->state= ATWORKING_FINETUNE;

	      if(verbose>1)
		{
		  printf("kbb.%d estimate (shift 1/3) timingGradStep=%gMV/m timingGrad=%gMV/m\n",
			 timeCount, timingGradStep/(megavolt/meter), timingGrad/(megavolt/meter));
		}

	      break;
	    }

	  case AUTOADJUST_UNKNOWN:  //ought to have been handled already
	  case AUTOADJUST_UNNEEDED:
	  case AUTOADJUST_OUTPUT:   //ought to have been handled already
	    {
	      show1step(rfdevice->state,"timeOffset=",rfdevice->timeOffset/ns,"ns");

	      rfdevice->state= ATWORKING_FINAL;  //done adjusting

	      if(rfdevice->whatFixed & ATFIX_OFFSETINCR)  
		rfdevice->timeOffset+= rfdevice->timeIncrement;

	      if(verbose>1)
		{
		  printf("kbb.%d estimate? what2adjust=0x%X\n", timeCount, rfdevice->what2adjust);
		}

	      break;
	    }	      
	  }

	break;
      }
      
    case ATWORKING_FINETUNE:  //fine tune
      {
	//offset has been determined; working on adjustments from there
	//start parked at predicted phase
	
	G4double doup=0,dout=0;  //diffence between previous,current goal and result
	
	if(rfdevice->whatFixed & ATFIX_P)
	  {
	    dout= chng.Pout - rfdevice->autoTimingFixP;
	    doup= previous.Pout - rfdevice->autoTimingFixP;
	    show1step(rfdevice->state,"Pout=",chng.Pout/MeV,"MeV/c");
	  }
	else if(rfdevice->whatFixed & ATFIX_DE)
	  {
	    dout= chng.dE - rfdevice->autoTimingFixDE;
	    doup= previous.dE - rfdevice->autoTimingFixDE;
	    show1step(rfdevice->state,"dE=",chng.dE/MeV,"MeV");
	  }
	else if(rfdevice->whatFixed & ATFIX_DT)
	  {
	    dout= chng.dT - rfdevice->autoTimingFixDT;
	    doup= previous.dT - rfdevice->autoTimingFixDT;
	    show1step(rfdevice->state,"dT=",chng.dT/ns,"ns");
	  }
	else if(rfdevice->whatFixed & ATFIX_XDEFL)
	  {
	    dout= chng.Xdefl - rfdevice->autoTimingFixXdeflection;
	    doup= previous.Xdefl - rfdevice->autoTimingFixXdeflection;
	    show1step(rfdevice->state,"dX=",chng.Xdefl/deg,"deg");
	  }
	else if(rfdevice->whatFixed & ATFIX_YDEFL)
	  {
	    dout= chng.Ydefl - rfdevice->autoTimingFixYdeflection;
	    doup= previous.Ydefl - rfdevice->autoTimingFixYdeflection;
	    show1step(rfdevice->state,"dY=",chng.Ydefl/deg,"deg");
	  }
	else
	  {
	    show1step(rfdevice->state,"dE=",chng.dE/MeV,"MeV");
	  }

	G4double same= fabs(dout) < fabs(doup);    //KBB - simple search 
	
	switch(rfdevice->what2adjust)
	  {
	  case AUTOADJUST_PHASE:    //fix gradient & output
	    {
	      AnotherPass= fabs(timingPhaseStep/rfdevice->omega) > rfdevice->timingTolerance ||
		fabs(dout) > rfdevice->autoTimingFailureLimit;
	      
	      if(AnotherPass)
		{
		  if(!same)  timingPhaseStep/= -1.5;  //backup and reverse direction, reduce step size
		  timingPhase+= timingPhaseStep;
		  rfdevice->timeOffset= timingZero - timingPhase/rfdevice->omega;

		  if(verbose>1)
		    {
		      printf("kbb.%d FINETUNE timingPhaseStep=%gdeg timingPhase=%.9fdeg Another=Y\n", 
			     timeCount, timingPhaseStep/deg, timingPhase/deg);
		    }
		}
	      else
		{
		  if(verbose>1)
		    {
		      printf("kbb.%d FINETUNE timingPhase=%.9fdeg Another=N\n", 
			     timeCount,timingPhase/deg);
		    }

		  rfdevice->state= ATWORKING_FINAL;  //done adjusting

		  if(rfdevice->whatFixed & ATFIX_OFFSETINCR)  
		    rfdevice->timeOffset= rfdevice->timeOffset+= rfdevice->timeIncrement;
		}
	      
	      break;
	    }
	    
	  case AUTOADJUST_GRADIENT:
	    {
	      AnotherPass= fabs(dout) > rfdevice->autoTimingFailureLimit;
	      
	      if(AnotherPass)
		{
		  if(!same)  timingGradStep/= -1.5;     //backup and reverse direction, decrease step size
		  timingGrad+= timingGradStep;
		  rfdevice->maxGradient= timingGrad;

		  if(verbose>1)
		    {
		      printf("kbb.%d FINETUNE timingGradStep=%gMV/m timingGrad=%.9fMV/m Another=Y\n", 
			     timeCount, timingGradStep/(megavolt/meter), timingGrad/(megavolt/m));
		    }
		}
	      else
		{
		  if(!RevisitTimingAfterGradientAdjust || revisited_timingZero)
		    {
		      if(verbose>1)
			{
			  printf("kbb.%d FINETUNE timingGrad=%.9fMV/m Another=N Revisit=N\n", 
				 timeCount, timingGrad/(megavolt/m));
			}

		      lastFound_maxGradient= rfdevice->maxGradient= timingGrad;

		      rfdevice->state= ATWORKING_FINAL;  //done adjusting

		      if(rfdevice->whatFixed & ATFIX_OFFSETINCR)  
			rfdevice->timeOffset+= rfdevice->timeIncrement;
		    }
		  else  //start all over deterimining timeZero, but using current gradient
		    {	      
		      rfdevice->state= ATWORKING_OFFSET;
		      rfdevice->timeOffset= timingZero= 0;
		      timingPhase= automethod2timingphase(method, entrance.q);  //reset
		      revisited_timingZero= true;
		      timeCount000= timeCount;  //last of 1st of 2 passes

		      if(verbose>1)
			{
			  printf("kbb.%d:STARTOVER:(%s) start over using maxGradient=%.6f MV/m\n", 
				 timeCount, pbname, rfdevice->maxGradient/(megavolt/meter));
			}
		    }
		}

	      break;
	    }
	  case AUTOADJUST_UNKNOWN:  //ought to have been done already
	  case AUTOADJUST_UNNEEDED:
	  case AUTOADJUST_OUTPUT:   //ought to have been done already
	    {
	      if(verbose>1)
		{
		  printf("kbb.%d FINETUNE - don't know what?\n", timeCount);
		}

	      rfdevice->state= ATWORKING_FINAL;  //done adjusting

	      if(rfdevice->whatFixed & ATFIX_OFFSETINCR)  
		rfdevice->timeOffset+= rfdevice->timeIncrement;

	      break;
	    }
	  }

	break;
      }
      
    case ATWORKING_FINAL:   //found solution, but need to do one more pass to update all
      {
	G4bool ok=true;
	if(BLCommand::isUndefined(rfdevice->autoTimingFailureLimit))
	  {
	    ok= true;
	  }
	else if(rfdevice->whatFixed & ATFIX_P)
	  {
	    ok= fabs(rfdevice->autoTimingFixP-chng.Pout)/MeV <  rfdevice->autoTimingFailureLimit;
	    show1step(rfdevice->state,"Pout=",chng.Pout/MeV,"MeV/c");
	  }
	else if(rfdevice->whatFixed & ATFIX_DE)
	  {
	    ok= fabs(rfdevice->autoTimingFixDE-chng.dE)/MeV <  rfdevice->autoTimingFailureLimit;
	    show1step(rfdevice->state,"dE=",chng.dE/MeV,"MeV");
	  }
	else if(rfdevice->whatFixed & ATFIX_DT)
	  {
	    ok= fabs(rfdevice->autoTimingFixDT-chng.dT)/ns <  rfdevice->autoTimingFailureLimit;
	    show1step(rfdevice->state,"dT=",chng.dT/ns,"ns");
	  }
	else if(rfdevice->whatFixed & ATFIX_XDEFL)
	  {
	    ok= fabs(rfdevice->autoTimingFixXdeflection-chng.Xdefl)/deg <  rfdevice->autoTimingFailureLimit;
	    show1step(rfdevice->state,"dX=",chng.Xdefl,"deg");
	  }
	else if(rfdevice->whatFixed & ATFIX_YDEFL)
	  {
	    ok= fabs(rfdevice->autoTimingFixYdeflection-chng.Ydefl)/deg <  rfdevice->autoTimingFailureLimit;
	    show1step(rfdevice->state,"dY=",chng.Ydefl,"deg");
	  }
	else
	  {
	    show1step(rfdevice->state,"dE=",chng.dE/MeV,"MeV");
	  }	
	
	printf("rfdevice::autoTiming:%d",redo);

	printf("[%s,%.1fdeg]ok=%s %s ",
	       automethod2text(method), automethod2timingphase(method,entrance.q)/deg,
	       ok? "Yes":"No", pbname);
	
	printf("timeOffset=%.9fns", rfdevice->timeOffset/ns);

	if(rfdevice->whatFixed & ATFIX_OFFSETINCR)
	  printf("{timeIncrement=%.9fns}",rfdevice->timeIncrement/ns);
	
	printf(",maxGradient=%gMV/m,phaseAcc=%gdegRF -> ",
	       rfdevice->maxGradient/(megavolt/meter),rfdevice->phaseAcc/deg);
	
	printf("%gMeV/c,%gMeV,%gns,%gXdeg,%gYdeg ",
	       chng.Pout/MeV,chng.dE/MeV,chng.dT/ns,chng.Xdefl/deg,chng.Ydefl/deg);
	
	printf("timingZero=%.6fns,exit@%gmm,%gmm,%gmm,%gns\n",
	       timingZero/ns,exit.x/mm,exit.y/mm,exit.z/mm,exit.t/ns);

	if(!ok)
	    G4Exception("rfdevice","Final fixed value out of limits",FatalException,"");  //abort

	TrackOfInterest= true;

	if(verbose>1)
	  {
	    printf("kbb.%d FINAL - all done\n", timeCount);
	    fflush(stdout);
	  }

	/* lock final values into rfdevicefield, as copies of rfdevicefield reuse rfdevice... */
	maxGradient= rfdevice->maxGradient;
	timeOffset= rfdevice->timeOffset;
	stateOfPlacement= rfdevice->state= ATWORKING_DONE;
	break;
      }

    case ATWORKING_DONE:  
      {
	maxGradient= rfdevice->maxGradient;
	timeOffset= rfdevice->timeOffset;
	stateOfPlacement= rfdevice->state= ATWORKING_DONE;
	return;	/* continue onward */
      }	
      
    case ATWORKING_UNKNOWN:  //never ought to happen
    case ATWORKING_INFLUX:
    case ATWORKING_RESET:
    default:
      {
	if(verbose>1)
	  {
	    printf("kbb.%d FINAL - ought not to be here\n", timeCount);
	  }
      }

      return;
    }

    previous= chng;
 
    // restore saved data (i.e. jump back to when the track 
    // first entered TimingVol)

    if(!validSaveTrack)
      G4Exception("rfdevice","Invalid Step",FatalException,"");
 
    track->SetTrackStatus(fStopAndKill);  //kill off current track

    steppingMgr->GetfSecondary()->push_back(new G4Track(saveTrack));  //KBB* g4bl 2.06-2.08 change
 
    if(verbose>0)  fflush(stdout);

    return;
}

#else // G4BL_GSL
int dummyrfdevice=0;
#endif // G4BL_GSL
