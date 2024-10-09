//	BLNTuple.hh
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

#ifndef BLNTUPLE_HH
#define BLNTUPLE_HH

#include <map>
#include <vector>
#include <set>
#include "globals.hh"
#include "G4ThreeVector.hh"

class BLNTupleHandler; // forward declarations
class BLNTupleCallback;

/**	class BLNTuple is an interface class to an NTuple, independent of
 *	implementation (HistoScope, BLTrackFile, Root, ...).
 *
 *	Each type of NTuple will derive a class from this one, and implement
 *	the virtual functions.
 **/
class BLNTuple {
protected:
	G4String name;
	G4String fields;
	int nData;
	std::vector<BLNTupleCallback*> callback;
	/// Constructor.
	BLNTuple(G4String _name, G4String _fields) :
		name(_name), fields(_fields), callback() { nData = 0; }
	static bool enableOutput;
	static std::map<G4String,BLNTupleHandler*> *handler;
	static std::vector<BLNTuple *> list;
	static std::set<G4String> nameSet;
	static G4String defaultType;
	static BLNTupleHandler *forceHandler;
	friend class BLNTupleHandler;
public:
	/// Destructor.
	virtual ~BLNTuple() { callback.clear(); }

	/// create() will create a new NTuple for writing.
	/// type must be a known type of NTuple implementation ("" => default).
	/// "category/name" is the name of the NTuple (the "/" may be omitted 
	/// if category is null, depending on type).
	/// fields is a : separated list of the field names.
	/// filename is a hint, and is not used by all types.
	static BLNTuple *create(G4String type, G4String category, G4String name,
					G4String fields, G4String filename);

	/// read() will open an NTuple for reading.
	/// type must be a known type of NTuple implementation ("" is default).
	/// "category/name" is the name of the NTuple (the "/" may be omitted 
	/// if category is null, depending on type).
	/// fields is a : separated list of the field names (not used by all
	/// types, in which case it is up to the user to make sure the file
	/// fields match the required ones).
	/// filename is a hint, and is not used by all types.
	/// The NTuple should be explicitly closed by calling its close().
	static BLNTuple *read(G4String type, G4String category, G4String name,
					G4String fields, G4String filename);

	/// appendRow() adds a row to an NTuple created for writing.
	/// NOTE: data[] must contain n doubles, and n must be the same
	/// as the number of fields in the NTuple.
	virtual void appendRow(double data[], int n) = 0;

	/// readRow() reads a row from an NTuple open for reading.
	/// returns true if valid row, false if EOF or error.
	/// ndata must be equal to the value of getNData(), and data[]
	/// must be at least that long.
	virtual bool readRow(double data[], int ndata) = 0;

	/// getName() returns the name of the NTuple (not including category)
	G4String getName() const { return name; }

	/// getFields() returns a colon-separated list of the fields.
	G4String getFields() const { return fields; }

	/// getNData() returns the # data items required.
	int getNData() const { return nData; }

	/// doCallbacks() will call all registered callbacks
	void doCallbacks(double data[], int ndata);

	/// annotate() will add an annotation line to any ASCII format.
	/// Normally line should begin with a "#", and have no "\r" or "\n".
	/// line == "" will output an empty line to ASCII files.
	/// Ignored for NTuples that don't support annotations.
	virtual void annotate(G4String line) { }

	/// needsReference() returns true if this NTuple needs the Reference
	/// particle.
	virtual bool needsReference() { return false; }

	/// flush() will flush the NTuple to disk, if supported.
	virtual void flush() = 0;

	/// close() will close the NTuple.
	virtual void close() = 0;

	/// do Summary() will print a 1-line summary to stdout, if supported.
	virtual void doSummary() = 0;

	/// registerCallback() registers a callback with this BLNTuple.
	void registerCallback(BLNTupleCallback *cb) { callback.push_back(cb); }

	/// getList() returns the list of all BLNTuples.
	static std::vector<BLNTuple*> getList() { return list; }

	/// isUniqueName() returns true if the name is unique. Also enters
	/// the name into the nameSet.
	static bool isUniqueName(G4String name);

	/// disableOutputFiles() prevents any output files from being written.
	static void disableOutputFiles() { enableOutput = false; }

	/// getEnableOutput() returns enableOutput;
	static bool getEnableOutput() { return enableOutput; }

	/// registerHandler() registers a BLNTupleHandler. This defines
	/// a new type of NTuple that can be created. Can be called in
	/// static initializers without worrying about ordering.
	static void registerHandler(G4String typeName, BLNTupleHandler *h);

	/// getHandler() returns the BLNTupleHandler for a given type.
	/// returns 0 for unknown type.
	static BLNTupleHandler *getHandler(G4String type);

	/// flushAll() writes all NTuples to file(s). This is a "checkpoint
	/// dump".
	/// Loops over all NTuples in the order they were created.
	/// Ignores NTuples open for reading.
	static void flushAll();

	/// closeAll() writes all NTuples to file(s) and closes them.
	/// Loops over all NTuples in the order they were created.
	/// Ignores NTuples open for reading.
	static void closeAll();

	/// summary() will summarize the NTuples to stdout
	/// Loops over all NTuples in the order they were created.
	/// Ignores NTuples open for reading.
	static void summary();

	/// update() will update real-time histogram interface(s)
	static void update();

	/// setDefaultType() sets the default type for NTuples
	static void setDefaultType(G4String name) { defaultType = name; }

	/// getFormatList() returns a list of known formats.
	static G4String getFormatList();
	
	/// getForceHandler() will return the forced handler (NULL if none).
	static BLNTupleHandler *getForceHandler() { return forceHandler; }
	
	/// setForceHandler() will force the specified handler to be used,
	/// regardless of type. Used to force an MPI handler when in MPI mode.
	static void setForceHandler(BLNTupleHandler *f) { forceHandler = f; }
};

/** class BLNTupleHandler defines a handler for a single type of NTuple.
 *  To be useful it must be registered with BLNTuple::registerHandler()
 *  -- this is normally done in the class default constructor, and a
 *  static object is defined.
 **/
class BLNTupleHandler {
protected:
	std::map<G4String,BLNTupleHandler*> *handler()
		{ return BLNTuple::handler; }
public:
	virtual BLNTuple *create(G4String type, G4String category, 
			G4String name, G4String fields, G4String filename) = 0;
	virtual BLNTuple *read(G4String type, G4String category, G4String name,
				G4String fields, G4String filename) = 0;
	virtual bool ignoreMultipleCreates() { return false; }
};

/** class BLNTupleCallback defines a callback for BLNTuple::appendRow().
 **/
class BLNTupleCallback {
public:
	/// callback() is called by BLNTuple::appendRow() for every registered
	/// callback. The first argument is always this, identifying the
	/// BLNTuple that is doing the calling.
	virtual void ntupleCallback(BLNTuple *ntuple, double data[], int ndata)
									= 0;
};

#endif // BLNTUPLE_HH
