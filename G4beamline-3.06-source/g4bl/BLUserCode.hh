//	BLUserCode.hh
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

#ifndef BLUSERCODE_HH
#define BLUSERCODE_HH

#include <vector>

#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4Track.hh"

/**	class BLUserCode is a base class for a user-defined classes.
 *
 *	See the section in the User's Guide titled "User Code for various
 *	commands".
 *
 *	This is used as follows: a G4beamline command that needs
 *	to call user-supplied code will have a name, and it will expect
 *	all user-supplied code to arrange so getType() returns that name.
 *	The G4beamline author will write a class BLUserXXX that is
 *	derived from BLUserCode and implements getType() with the
 *	correct value -- so users should NOT override this function.
 *	Usually this class declaration will be put into this file.
 *	BLUserXXX will also define some additional virtual function(s)
 *	for the user to implement, appropriate to the function to be 
 *	performed. The user must implement these functions and getName().
 *	There's one other requirement: the user class must call 
 *	registerUserCode() at the end of its default constructor, and must
 *	define a static instance of the class using that constructor --
 *	this provides linkage into G4beamline so that the underlying command
 *	can find and call the user code.
 *
 *	The user code should be compiled into a shared object, using the
 *	g4blmake command. The G4beamline input file should use the "load"
 *	command to load the shared object BEFORE any command references
 *	the user class.
 *	NOTE: the "load" command uses the filename built by g4blmake;
 *	the other command(s) use the getName() implemented by
 *	the user class. Multiple instances of user code can be put into a
 *	single .cc file (or multiple .cc files linked into one shared object);
 *	they must each have a unique getName() return value.
 *
 *	NOTE: this class will interface between shared objects, so no
 *	inline function can reference any static data member. This does NOT
 *	apply to user-supplied derived classes that necessarily reside in a 
 *	single shared object.
 **/
class BLUserCode {
	static std::vector<BLUserCode*> *list;
public:
	/// Constructor. 
	BLUserCode() { }

	/// Destructor.
	virtual ~BLUserCode() { }

	/// registerUserCode() -- register this instance with the BLManager so
	/// other commands can find it. Note that every final derived
	/// class must call this at the end of its default constructor,
	/// and must define a static object using that constructor.
	void registerUserCode();

	/// Name of the type -- MUST be overridden in derived classes that
	/// define an interface to user code from a G4beamline command.
	virtual const char *getType() = 0;

	/// Name of the instance -- MUST be overridden in user-supplied
	/// derived classes.
	virtual const char *getName() = 0;
};

/**	class BLUserTrackFilter -- used by BLCMDusertrackfilter.
 *
 *	User code can use any library used by G4beamline simply by putting the
 *	appropriate #include into the .cc file and building with g4blmake --
 *	it automatically adds the appropriate include paths and libraries.
 *	These are: CLHEP, Geant4, Root, Gsl, and the C/C++ libraries.
 *	Coin and Xwindows are used by G4beamline for visualization, but are
 *	not available to user code.
 *
 *	Additional user libraries can be used if they are added to the
 *	CCFLAGS and EXTRALIBS environment variables to g4blmake -- do NOT
 *	do this for libraries used by G4beamline. Do not attempt to use the
 *	Coin or Xwindows libraries.
 *
 *	The user .cc file(s) can contain local class definitions and static 
 *	instances of them; the latter will be constructed when the shared-
 *	object is loaded. NOTE: there is no guarantee that their destructors
 *	will be called (and they usually aren't), so if destruction is
 *	important, use setup() and complete().
 *
 *	filter() will be called by the "usertrackfilter" element
 *	whenever a track enters the physical volume of the element. The
 *	input track will be setup. The function can kill the track, and/or
 *	add pointers to new G4Track-s to the secondaries vector, as necessary.
 *
 *	Note that the Tune particle has eventID -2 and the Reference particle
 *	has eventID -1 -- these events may need special treatment, and
 *	their tracks should not be killed or modified.
 *
 *	The secondaries vector is initially empty; add pointers to G4Track-s
 *	to it to create secondaries.
 *
 *	verbose is the value of the parameter steppingVerbose. If it is 
 *	nonzero you should print debugging information.
 *
 *	Note: for secondaries the trackID will be set to a unique value and
 *	the parentID will be set to the trackID of the input track,
 *	after filter() returns; all secondaries remain part of the
 *	current event. So there's no need to set trackID or parentID
 *	in G4Track-s put into the secondaries vector, and whatever values
 *	they have will be ignored.
 *
 *	The secondaries are created inside the volume, and do not
 *	get filtered.
 **/
class BLUserTrackFilter : public BLUserCode {
public:
	/// Constructor
	BLUserTrackFilter() : BLUserCode() { }

	/// Destructor
	virtual ~BLUserTrackFilter() { }

	/// getType() return "usertrackfilter"
	const char *getType() { return "usertrackfilter"; }

	// Convenience functions for user code to use:

	/// getPDGid() returns the PDGid of a G4Track.
	/// User code can of course use the usual Geant4
	/// functions to do this, this is merely for convenience.
	static int getPDGid(G4Track *track);

	/// setMomentum() will set the momentum of a G4Track.
	/// User code can of course use the usual Geant4
	/// functions to do this, this is merely for convenience.
	static void setMomentum(G4Track *track, G4ThreeVector momentum);

	/// getMass() returns the mass of a particle with known PDGid.
	/// User code can of course use the usual Geant4
	/// functions to do this, this is merely for convenience.
	static G4double getMass(int PDGid);

	/// killTrack() will kill a track.
	static void killTrack(G4Track *track);

	/// constructTrack() creates a new G4Track from position, momentum,
	/// time, and PDGid. Global coordinates are used. User code can of
	/// course use the usual Geant4 functions to do this, this is merely
	/// for convenience. Note that eventID is not contained in a G4Track.
	static G4Track *constructTrack(G4ThreeVector position,
			G4ThreeVector momentum, G4double time, int PDGid,
			G4double weight=1.0);

	// User-supplied virtual functions:

	/// getName() returns the name of the user-supplied instance
	virtual const char *getName() = 0;

	/// user function to filter tracks.
	/// Global coordinates are used.
	/// NOTE: do NOT delete the input track, or those you create and
	/// put into secondaries -- they are deleted by Geant4 internally.
	virtual void filter(G4Track *track, int eventID,
			std::vector<G4Track*> &secondaries, int verbose) = 0;

	/// setup() will be called during initialization; its
	/// argument is the string passed to the "usertrackfilter" command
	/// as its "init" argument. The init string can be used in any way
	/// by the user code (it passes information from the G4beamline input
	/// file to the user code, in a manner defined by the user code; the
	/// input file should be written to match this usage).
	/// If desired, multiple usertrackfilter commands can use a single
	/// user-supplied class instance, using the init string to
	/// differentiate among them. The init string need not be used at all.
	virtual void setup(const char *init) = 0;

	/// complete() will be called after tracking is complete.
	/// It is intended for printing summaries and generally closing down.
	/// If multiple usertrackfilter commands reference this class instance,
	/// complete() will be called for each one, with the corresponding 
	/// init string (do not depend on any particular ordering).
	virtual void complete(const char *init) = 0;
};

//
// Below are definitions for the convenience of the user-defined functions
// above.
//

/** 	Particle IDs for stable particles. Use negative for anti-particle.
 *	ELECTRON and MUON are negative particles; PION, KAON, and PROTON
 *	are positive particles; KAON0, NEUTRON, and NU_* are neutrals.
 *	There are many other PDGid-s, listed in the Geant4 and G4beamline 
 *	documentation; these are just the most common ones.
 **/
const int ELECTRON=11, MUON=13, PION=211, KAON=321, KAON0=311, NEUTRON=2112,
	PROTON=2212, NU_E=12, NU_MU=14, NU_TAU=16, GAMMA=22 ;

#endif // BLUSERCODE_HH 
