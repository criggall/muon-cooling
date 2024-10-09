//	Param.hh
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

#ifndef BLPARAM_HH
#define BLPARAM_HH

#include <math.h>
#include "globals.hh"
#include "G4Types.hh"
#include "BLCommand.hh"

/**	BLParam - parameter class
 *
 *	The global object Param is used to access parameter values defined
 *	via the param command. The param command and the setParam()
 *	functions create parameter name/value pairs known to the Param object.
 *
 *	The BLSetParam class can be used to set initial values of parameters.
 *
 *	The internal representation of every value is a G4String, but that
 *	can be converted to double or int via get...() functions.
 *
 *	NOTE: the get*() functions are intended to be used during 
 *	initialization, not during tracking. They perform a format conversion
 *	and can be expensive -- a single call to getDouble() in
 *	UserSteppingAction() slowed tracking down by MORE THAN A FACTOR OF TEN!
 *
 *	NOTE: static and global constructors CAN use Param; this is
 *	explicitly allowed because init() guarantees initialization order.
 **/
class BLParam {
	static std::map<G4String,G4String> *paramMap;
	static std::map<G4String,G4String> *paramHelpText;
	static void init();
public:
	/// Constructor. 
	BLParam() { init(); }

	/// printParam() will display all defined parameters on stdout.
	void printParam();

	/// getString() returns the G4String value of a parameter.
	/// If name is not defined as a parameter, the environment is searched.
	/// Returns "" for unknown parameters, and displays the error on stdout.
	/// DO NOT USE DURING TRACKING!
	G4String getString(G4String name);

	/// getDouble() returns the G4double value of a parameter.
	/// returns -HUGE_VAL for unknown or invalid-format parameters,
	/// and displays the error on stdout.
	/// DO NOT USE DURING TRACKING!
	G4double getDouble(G4String name);

	/// getInt() returns the G4int value of a parameter.
	/// returns 0x80000000 for unknown or invalid-format parameters, and
	/// displays the error on stdout.
	/// DO NOT USE DURING TRACKING!
	G4int getInt(G4String name);

	/// setParam() sets the value of a parameter.
	/// DO NOT USE DURING TRACKING!
	void setParam(G4String name, G4String value);

	/// setParam() sets the value of a parameter.
	/// DO NOT USE DURING TRACKING!
	void setParam(G4String name, G4double value);

	/// setParam() sets the value of a parameter.
	/// DO NOT USE DURING TRACKING!
	void setParam(G4String name, G4int value);

	/// isDefined() returns true if the name has a defined value.
	bool isDefined(G4String name);

	/// expand() will expand every $name in str with parameter values.
	/// DO NOT USE DURING TRACKING!
	G4String expand(G4String str);

	/// setHelpText() will set help text for a parameter
	/// DO NOT USE DURING TRACKING!
	void setHelpText(G4String name, G4String text);

	/// getHelpText() returns the (complete) help text string.
	/// DO NOT USE DURING TRACKING!
	G4String getHelpText();
};

extern BLParam Param;	/// the global Param object

/**	class BLSetParam -- Parameter initialization class
 *
 *	Classes that use Parameters should declare a static instance of
 *	this class to provide a default value for each such parameter.
 *	Take care that the value is of the proper format.
 *
 *	These values are set initally, so "param -unset" will not set them.
 **/
class BLSetParam {
public:
	BLSetParam(G4String name, G4String value, G4String helpText="")
		{ Param.setParam(name,value); 
		  if(helpText.size() > 0) Param.setHelpText(name,helpText);
		}
};

#endif // BLPARAM_HH
