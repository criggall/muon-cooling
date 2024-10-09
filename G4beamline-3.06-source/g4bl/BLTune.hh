//	BLTune.hh
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

#ifndef BLTUNE_HH
#define BLTUNE_HH

#include <map>
#include <set>
#include "globals.hh"
#include "BLEvaluator.hh"


/**	class BLTune contains the general interface to tuning arguments.
 *
 *	The tune command (class BLCMDtune) will tune arguments to commands
 *	with the goal of achieving some specified measure of performance.
 *	Tuning occurs when the Tune Particle is tracked, and only those
 *	aspects of performance that can be achieved during that tracking
 *	can be tuned. The example used here is to tune the current in a
 *	genericbend so the reference particle had dxdz=0 after the magnet.
 *
 *	Tunable arguments to a command must:
 *	 a) be declared with argTunable() in defineArgs()
 *	 b) be referenced indirectly by each placement class of the element;
 *	    that is, the placement object cannot keep a local copy, but
 *	    must always reference the actual parameter variable of the Element.
 *	 c) be copied using BLTune::copyTunableArg() in the element's copy
 *	    constructor
 *
 *	Tunable arguments are set on either the element or place command
 *	using the usual name=expression syntax. The expression can contain
 *	the name of any enclosing tune command as a variable, and the value
 *	of the variable will always be the current value during tuning, and
 *	will be the final result after tuning is complete.
**/
class BLTune {
	static BLEvaluator eval;
	static std::map<G4double*,G4String> tuneExpr;
	static std::map<G4double*,G4double> tuneUnits;
	static std::set<G4String> names;
public:
	/// Constructor
	BLTune() { }

	/// Destructor
	virtual ~BLTune() { }

	/// defineTunableArg() will define a tunable argument.
	static void defineTunableArg(G4double& arg, G4double units, G4String expr);

	/// copyTunableArg() will copy a tunable argument; for use in the copy
	/// constructor of Element classes.
	static void copyTunableArg(G4double *newArg, const G4double *oldArg);

	/// isSet() returns true if the varaible name is set.
	static bool isSet(G4String name);

	/// set() will set the value of a tune variable, and update all
	/// arguments that use it.
	static void set(G4String& name, G4double value);

	/// unset() will remove a tune variable.
	static void unset(G4String& name);

	/// update() will update all tunable args with their current values.
	static void update();

	/// getNames() will return the set of names.
	static const std::set<G4String> &getNames() { return names; }

	/// getValue() will return the value of a name.
	/// Value is effectively in internal Geant4 units (tune variables
	/// don't have units, the argument expressions do).
	static double getValue(G4String& name) 
		{ return eval.evaluate(name.c_str()); }
};

#endif // BLTUNE_HH
