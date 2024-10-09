//	BLCommand.hh
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

#ifndef BLCOMMAND_HH
#define BLCOMMAND_HH

#include <vector>
#include <map>
#include <iostream>
#include <fstream>

#include "globals.hh"
#include "G4RotationMatrix.hh"
#include "G4VisAttributes.hh"
#include "G4Material.hh"

#include "BLEvaluator.hh"

typedef std::vector<G4String> BLArgumentVector;
typedef std::map<G4String,G4String> BLArgumentMap;

const G4String IndentDesc("             ");
const G4String IndentArg("                           ");

const G4String DefaultDoubleFmt = "%.4g";

/**	class BLCommandPos - postion of the command input
 *	(used to get the input position for do loops)
 *	see also: BLCommand::setPos().
 **/
class BLCommandPos {
	std::istream *in;
	std::streampos pos;
	friend class BLCommand;
public:
	BLCommandPos(std::istream *_in) {
		in = _in;
		pos = in->tellg();
	}
	bool isValid() { return pos != (std::streampos)-1; }
};

/**	enum BLCmdType -- types of commands (for help index).
 **/
enum BLCmdType {BLCMDTYPE_HELP,BLCMDTYPE_CONTROL,BLCMDTYPE_LAYOUT,
		BLCMDTYPE_BEAM,BLCMDTYPE_AUX,BLCMDTYPE_ELEMENT,
		BLCMDTYPE_PLACE,BLCMDTYPE_CUTS,BLCMDTYPE_DATA,
		BLCMDTYPE_PHYSICS,BLCMDTYPE_OTHER };


/**	class BLCommand - interface to define a command.
 *
 *	The standard usage for a command without a BLElement is to derive
 *	from this class and define a static instance using the default 
 *	constructor. In the derived-class default constructor the command is 
 *	registered, calls to setSynopsis() and setDescription() provide help
 *	text, and initial values are given to all class variables. This
 *	instance is used to execute the command, and normally no other
 *	instances are needed.
 *
 *	The static instance is used to keep the default arguments for the
 *	element.
 *
 *	Standard usage for a command for a BLElement is to derive a class
 *	from BLElement (which is derived from BLCommand), and do the above
 *	in its default constructor. The derived command() function then
 *	creates an instance of the class. Only the static instance is used
 *	to execute the command, but other instances implement the various
 *	elements (to be placed later in the input file).
 *
 *	See BLCMDdemo.cc and BLCMDtubs.cc for examples (.hh files are
 *	not needed for most commands, as their classes are not needed
 *	anywhere but in their implementation file).
 *
 *	NOTE: ONLY the default object should use the default constructor!
 **/
class BLCommand {
	enum TokenType { NONE, ARGNAME, ARGVALUE };
	static std::map<G4String,BLCommand*> *mapCommand;
	static int errors;
	/// nextToken() will return the next token in line; called only by
	/// doCommand(). Side-effect: updates place and sets type.
	static G4String nextToken(const G4String& line, 
				G4String::size_type& place, TokenType& type);
	/// printArgDesc will print the description of an argument, using
	/// wrapWords().
	void printArgDesc(G4String name, G4String _description);
	static std::istream *in;
	friend class BLCMDplace;
protected:
	enum ArgumentMode { PROCESS, COUNT, HELP, PRINT, CHANGE };
	ArgumentMode argMode;
	bool argFound;
	G4String argName;
	G4String argValue;
	G4String synopsis;
	G4String description;
	G4String printArgString;
	int nArgs;
	int nFixed;
	int nTunable;
	BLCmdType type;
	friend class BLCMDhelp;
protected:
	/// Default Constructor.
	BLCommand() : synopsis(), description()
		{ }

	/// Destructor.
	virtual ~BLCommand();

	/// Copy Constructor.
	BLCommand(const BLCommand& r) : synopsis(), description()
		{ }

	// general functions for use by derived classes:

	/// registerCommand() registers the command.
	/// Used only in the default constructor -- ensures mapCommand is
	/// properly initialized regardless of order.
	void registerCommand(BLCmdType _type);

	/// deleteCommand() will delete a command by name
	void deleteCommand(G4String &name)
		{ if(mapCommand) mapCommand->erase(name); }

	/// setSynopsis() gives a 1-line (64 char) description of what the 
	/// command does. Used by the help command with no arguments.
	/// The string should be 64 characters or less, and should not contain
	/// any newline; no word wrapping.
	/// Used only in the default constructor.
	void setSynopsis(const char * _synopsis) 
		{ synopsis = _synopsis; }

	/// setDescription() gives a more complete description of the command,
	/// when help for this specific command is requested.
	/// When printed by the help() function, word wrapping is performed,
	/// but lines beginning with whitespace are kept intact.
	/// Named arguments are automatically described by the help() function,
	/// but positional arguments should be described here.
	/// Words are wrapped using wrapWords().
	/// Used only in the default constructor.
	void setDescription(const char * _description) 
		{ description = _description; }

	/// argString() declares a string-valued argument. name should be a
	/// valid C identifier 15 chars or less, and description should be
	/// a short description (wrapWords() is used).
	/// Used only in the defineNamedArgs() function.
	/// permitChange affects whether this argument can be set when
	/// an element is placed.
	void argString(G4String& var, G4String name, G4String description,
			bool permitChange=true);

	/// argDouble() declares a G4double-valued argument. name should be a
	/// valid C identifier 15 chars or less, and description should be
	/// a short description (wrapWords() is used).
	/// Used only in the defineNamedArgs() function.
	/// units should be the human units for the argument, because
	/// the value of var is always in internal units. fmt and units are
	/// used when argMode=PRINT.
	/// permitChange affects whether this argument can be set when
	/// an element is placed.
	void argDouble(G4double& var, G4String name, G4String description,
		G4double units=1.0, G4String fmt=DefaultDoubleFmt,
		bool permitChange=true);

	/// argTunable() declares a G4double-valued argument that is tunable 
	/// (see the tune command in BLCMDtune.cc). name should be a
	/// valid C identifier 15 chars or less, and description should be
	/// a short description (wrapWords() is used).
	/// Used only in the defineNamedArgs() function.
	/// units should be the human units for the argument, because
	/// the value of var is always in internal units. fmt and units are
	/// used when argMode=PRINT.
	/// Tunable arguments can always be set when an element is placed,
	/// and will be changed during the tracking of the Tune Particle
	/// if their value includes a tune variable. See BLTune.hh for a 
	/// description of how tunable arguments are handled.
	void argTunable(G4double& var, G4String name, G4String description,
		G4double units=1.0, G4String fmt=DefaultDoubleFmt);

	/// argInt() declares a integer-valued argument. name should be a
	/// valid C identifier 15 chars or less, and description should be
	/// a short description (wrapWords() is used).
	/// Used only in the defineNamedArgs() function.
	/// permitChange affects whether this argument can be set when
	/// an element is placed.
	void argInt(G4int& var, G4String name, G4String description,
			bool permitChange=true);

	/// handleNamedArgs() handles all named arguments to the command (i.e.
	/// those defined in defineNamedArgs()). For use in the command() 
	/// function.
	/// returns 0 if OK, -1 on error.
	int handleNamedArgs(BLArgumentMap& namedArgs);

	/// help() prints the help information for this command.
	/// detailed=false prints only the name and synopsis; true prints
	/// the description and arguments also.
	virtual void help(bool detailed);

	/// wrapWords() will wrap words for descriptions in help().
	/// In text, lines beginning with whitespace are preserved as is,
	/// including blank lines. Otherwise, newlines aren't needed in text.
	/// indent1 is the start of the first line, and indent is used
	/// at the start of succeeding lines ("" is OK); indents are counted
	/// in the line width. The returned string always ends with "\n".
	G4String wrapWords(G4String text, G4String indent1, G4String indent,
							G4String::size_type width=79);

	/// printArgs() will print (to stdout) the named arguments of the
	/// command entered, for log purposes.
	/// Calls defineNamedArgs() with argMode=PRINT to print the args.
	/// indent1 is the indent for the first line.
	void printArgs(G4String indent1);
public:
	// Functions implemented by the derived class:

	/// commandName() returns the name of the command.
	virtual G4String commandName() = 0;

	/// command() executes the command associated with this element.
	/// argv and namedArgs come from parsing the command-line:
	/// argv has a vector of any positional arguments (without "=");
	/// namedArgs has a map of any named arguments of the form name=value.
	/// Returns 0 if OK, -1 on error.
	/// NOTE: for classes derived via BLElement.hh, arguments can be given
	/// when an object is defined and when it is placed. For that to work,
	/// argument values must not be used in the command() function unless
	/// they are defined with permitChange=false (in the above argString...
	/// functions). Also the correct units must be given to argDouble().
	virtual int command(BLArgumentVector& argv,
						BLArgumentMap& namedArgs) = 0;

	/// defineNamedArgs() defines the named arguments for the command. 
	/// This function should consist ONLY of a series of calls to
	/// argString(), argDouble(), and argInt().
	/// Used by handleNamedArgs(), help(), and print() (they manipulate
	/// the behavior of the argXXX() functions via argMode).
	virtual void defineNamedArgs();

	/// argChanged() will be called whenever some argument value has
	/// changed. Derived classes should define it and perform any internal
	/// computations that depend on argument values. Many/most derived
	/// classes won't need to use it. In particular, after argChanged()
	/// returns the functions getHeight(), getWidth() and getLength()
	/// must return the correct values for the current arguments.
	/// BEWARE: a single argument change may call this multiple times.
	/// Applies only to named arguments, not positional ones.
	virtual void argChanged() { }

	/// print() will print to stdout the command. The default versions
	/// will do for most commands and elements. name can be "".
	/// prints commandname() and name in the space of IndentArg, followed
	/// by the values of all named arguments, in the order they appear
	/// in defineNamedArgs().
	virtual void print(G4String name);

	/// Print() -- alternate version for commands with special arguments.
	/// does not use defineNamedArgs() to print arguments, but rather
	/// prints namedArgs.
	void print(G4String name, BLArgumentMap& namedArgs);

	/// getName() returns the element's name; do not confuse this with
	/// commandName(). Re-defined by BLElement.
	virtual G4String getName() const { return "??"; }

	/// undefined() returns a quiet Not-A-Number, for initializing
	/// arguments to determine if the user set them or not.
	static double undefined();

	/// isUndefined() determines wheter an argument has been set.
	static bool isUndefined(double v);
public:
	// general functions:

	/// getCmdType() returns the type of this command
	BLCmdType getCmdType() { return type; }

	/// isValidCommand() returns true if cmd is a valid command.
	static bool isValidCommand(G4String& cmd)
		{ return mapCommand != 0 && mapCommand->count(cmd) > 0; }

	/// doCommand() will parse line and perform a single command.
	/// Returns 0 if OK, -1 if error.
	static int doCommand(G4String& line);

	/// parseArgs() will parse the arguments of a command line.
	/// Arguments are appended to argv[] and added to namedArgs.
	/// Returns 0 if OK, -1 if syntax error.
	static int parseArgs(const G4String &line, BLArgumentVector &argv,
						BLArgumentMap &namedArgs);

	/// readFile() reads a file of commands and executes them.
	/// "-" is stdin. Returns the number of errors encountered,
	/// or -1 if filename cannot be read.
	static int readFile(G4String filename);

	/// printError() prints an error message, using sprintf-style args.
	/// It also increments the errorCount.
	/// DEPRECATED -- use G4Exception() instead.
	static void printError(const char *fmt, ...);

	/// getErrorCount() will return the number of errors so far.
	static int getErrorCount() { return errors; }

	/// getNextCommand() returns the next command. backslash-newline has
	/// been handled, and the line is stripped of whitespace at both ends.
	/// returns 0 if error or EOF. The pointer should NOT be deleted.
	static G4String *getNextCommand();

	/// getPos() returns the current position of the Command input stream
	static BLCommandPos getPos();

	/// setPos() will set the Command input stream to a specified position
	void setPos(BLCommandPos &pos);

	// general utility functions for use by other classes

	/// getMaterial() searches the G4Material list and returns the entry
	/// corresponding to name. Searches the NIST database if the material
	/// is not already defined. Prepends "G4_" if necessary.
	/// prints an error and exits if the material is not found, so all
	/// materials MUST be defined before they are used, unless they can
	/// be found in the geant4 NIST database.
	static G4Material *getMaterial(G4String materialName, 
							bool ignoreError=false);

	/// getVisAttrib() returns the appropriate G4VisAttributes.
	static G4VisAttributes *getVisAttrib(G4String color);

	/// stringToRotationMatrix() converts a string "X90,Y45" into a 
	/// G4RotationMatrix.
	/// This is an active rotation, in that the object is first rotated
	/// around the parent's X axis by 90 degrees, then the object is
	/// further rotated around the parent's Y axis by 45 degrees.
	/// The return value points to a G4RotationMatrix on the heap, so
	/// it is persistent. Angles are in degrees, can have decimals,
	/// and can be negative. Axes are X, Y, Z.
	/// An alternate specification is "A1,2,3,4" where [1,2,3] specify
	/// the axis of rotation, and 4 is the angle of rotation around it
	/// (in degrees); the axis need not be normalized.
	/// NOTE: This is expensive and should not be used during tracking.
	static G4RotationMatrix *stringToRotationMatrix(G4String rotation);

	/// dumpRotation() dumps info about a rotation matrix to stdout.
	/// NOTE: This is expensive and should not be used during tracking.
	static void dumpRotation(const G4RotationMatrix *rot, const char *str);

	/// splitString() splits a string at each instance of any character
	/// in delim, returning a vector<G4String>. If trim is true, each
	/// returned string has whitespace removed from its front and rear.
	/// Successive delimiters, or delimiters at the beginning or end,
	/// generate an empty string.
	/// NOTE: This is expensive and should not be used during tracking.
	static std::vector<G4String> splitString(const G4String &s, 
					const G4String &delim, bool trim=false);

	/// getDouble() converts a string expression to a double, returning
	/// nan if the string is not a valid expression using double constants
	/// and standard C operators and functions; whitespace is permitted as
	/// in C, but an empty or all-whitespace string is invalid.
	/// Callers can use std::isnan(double) (from <cmath>) to test
	/// for validity (the isnan() MACRO from <math.h> is unreliable).
	/// NOTE: This is expensive and should not be used during tracking.
	static double getDouble(const G4String &s)
		{ BLEvaluator e; return e.evaluate(s); }

	/// getList() converts a string containing a delimited list
	/// of double expressions to a vector<G4double>. Any character in
	/// delim separates entries. Functions and operators known to
	/// BLEvaluator can be used in each entry of the list.
	/// Empty entries (null or consisting only of whitespace) are invalid.
	/// Returns an empty list if any entry is an invalid expression.
	/// NOTE: This is expensive and should not be used during tracking.
	static std::vector<G4double> getList(const G4String &s,
							const G4String &delim);

	/// matchList() returns true if s matches any pattern in
	/// patternList, whch is a comma-separated list of file-glob patterns.
	/// NOTE: This is expensive and should not be used during tracking.
	static bool matchList(const G4String s, const G4String patternList);
};

#endif // BLCOMMAND_HH
