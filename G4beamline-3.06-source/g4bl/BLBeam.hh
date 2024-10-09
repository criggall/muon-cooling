//	BLBeam.hh
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

#ifndef BLBEAM_HH
#define BLBEAM_HH

#include "globals.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleDefinition.hh"
#include "Randomize.hh"

#include "BLManager.hh"
#include "BLCommand.hh"
#include "BLTime.hh"

/**	class BLBeam is an interface class for a beam.
 *
 *	NOTE: all functions are either inline or abstract, so no .cc file is 
 *	needed.
 **/
class BLBeam {
protected:
	G4int nEvents;
public:
	/// Constructor.
	BLBeam() { nEvents = 0; }

	/// Utility: setRandomSeedToGenerate() sets the random seed for
	/// the generation of an event.
	static void setRandomSeedToGenerate(int evid) // reference has evid<0
		{ switch(BLManager::getObject()->getPRNGSeedMethod()) {
		  case EVENT_NUMBER:
			CLHEP::HepRandom::setTheSeed(evid>0?evid:0x7FFFFFFE);
			break;
		  case TIME_US:
			CLHEP::HepRandom::setTheSeed(
				(BLTime::timems()&0x7FFFFFFF) | 1);
			// do it only once, to avoid duplication
			BLManager::getObject()->setPRNGSeedMethod(NO_SEED);
			break;
		  case NO_SEED:
		  	return;
		  }
		  CLHEP::RandGauss::setFlag(false);
		  G4UniformRand(); // eat one for luck (don't get seed back)
		}

	/// Utility: setRandomSeedToTrack() sets the random seed for
	/// the tracking of an event.
	static void setRandomSeedToTrack(int evid) // reference has evid<0
		{ setRandomSeedToGenerate(evid);
		  // allow 16 random numbers for generation
		  for(int i=0; i<16; ++i)
		  	G4UniformRand();
		}

	/// getNEvents() returns the # events to process.
	virtual int getNEvents() const { return nEvents; }

	/// init() will initialize internal variables.
	virtual void init() = 0;

	/// generateReferenceParticle() generates the reference particle.
	/// BLManager calls this once for the TuneParticle during event -2,
	/// and once for the ReferenceParticle during event -1.
	/// Should always call setRandomSeedToTrack() before returning.
	/// Returns true if event generated; false if not.
	virtual bool generateReferenceParticle(G4Event *event) = 0;

	/// nextBeamEvent() generates the next beam event.
	/// If generating randomly, should call setRandomSeedToGenerate()
	/// and use no more than 16 random numbers.
	/// If read from a file, it should call event->SetEventID().
	/// Should always call setRandomSeedToTrack() before returning.
	/// Returns true if event generated; false if no more events.
	virtual bool nextBeamEvent(G4Event *event) = 0;

	/// summary() will print a summary, if necessary.
	virtual void summary() { }
};

#endif// BLBEAM_HH
