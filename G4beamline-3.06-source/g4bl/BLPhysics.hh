//	BLPhysics.hh
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

#ifndef BLPHYSICS_HH
#define BLPHYSICS_HH

#include "globals.hh"
#include "G4VUserPhysicsList.hh"
#include "G4ProcessTable.hh"

#include "CLHEP/Units/SystemOfUnits.h"
using namespace CLHEP;

enum BLSetValue { FORCE_ON, FORCE_OFF, NORMAL };

/**	class BLPhysics is the interface class for selecting the physics
 *	processes to be used in the simulation.
 *
 *	Note: all functions are either inline or abstract, so no .cc file 
 *	is needed.
 **/
class BLPhysics {
protected:
	bool stochasticsEnabled;
	G4double minRangeCut;
	G4double maxTime;
public:
	/// default constructor.
	BLPhysics() { 
		stochasticsEnabled = true;
		minRangeCut = 1.0*mm;
		maxTime = 1.0*millisecond;
	}

	/// destructor.
	virtual ~BLPhysics() { }

	/// setDoStochastics() sets whether or not stochastic processes
	/// are to be enabled. Note that Decay is a stochastic process.
	/// The argument can be FORCE_ON, FORCE_OFF, and NORMAL (which
	/// means use the doStochastics parameter of the physics command).
	/// warn controls how warnings are issued: 0=none, 1=G4Exception,
	/// 2=G4Exception.
	virtual void setDoStochastics(BLSetValue value, G4int warn=1) = 0;

	/// getPhysicsList() returns the G4PhysicsList.
	virtual G4VUserPhysicsList *getPhysicsList() = 0;

	/// getStochasticsEnabled() returns true if stochastics are enabled,
	/// false if not.
	bool getStochasticsEnabled() { return stochasticsEnabled; }

	/// getRangeCut() returns the range cut.
	G4double getRangeCut() { return minRangeCut; }

	/// isStochasticProcess() returns true if the process is stochastic.
	bool isStochasticProcess(G4VProcess *process) 
		{ return isStochasticProcess(process->GetProcessName()); }

	/// isStochasticProcess() returns true if the process is stochastic.
	/// Ionization energy loss processes are considered to NOT be 
	/// stochastic, because they will be set to return their mean energy
	/// loss when stochastics are off.
	bool isStochasticProcess(G4String name) {
		if(name.find("Trans") < name.size()) return false;
		if(name.find("Ioni") < name.size()) return false;
		if(name.find("Limiter") < name.size()) return false;
		if(name.find("BLCMD") < name.size()) return false;
		return true;
	}

	/// setProcessEnable() will disable or enable all processes that
	/// contain pattern in their name.
	void setProcessEnable(G4String pattern, bool onoff) {
		G4ProcessTable *pt = G4ProcessTable::GetProcessTable();
		std::vector<G4String> &pnv = *pt->GetNameList();
		for(unsigned i=0; i<pnv.size(); ++i) {
			G4String name = pnv[i];
			if(name.find(pattern) >= name.size()) continue;
			pt->SetProcessActivation(name,onoff);
		}
	}

	/// setMaxTime() / getMaxTime() control the maximum tracking time.
	void setMaxTime(G4double v) { maxTime = v; }
	G4double getMaxTime() { return maxTime; }

	/// isSpinTrackingEnabled() returns true if so.
	virtual bool isSpinTrackingEnabled() { return false; }
};

#endif // BLPHYSICS_HH
