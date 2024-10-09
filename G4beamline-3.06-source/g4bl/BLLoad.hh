//	BLLoad.hh
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


#ifndef BLLOAD_HH
#define BLLOAD_HH

#include <vector>

#include "globals.hh"

#if defined(__linux) || defined(__APPLE__)
#include <dlfcn.h>      // for dlopen()
#endif

typedef void (*VoidFunctionPointer)();

/** class BLLoad will load a shared object.
 *
 *  All static object constructors in the shared object are executed
 *  while loading, and they can register the resources implemented.
 *  The corresponding destructors are probably not executed -- use atexit()
 *  or register an appropriate callback with BLManager..
 *
 *  NOTE: the release 2 build architecture makes this unnecessary (user
 *  shared objects are not used). Several .cc files have residual calls
 *  to this class, surrounded by #ifdef USE_SHARED_OBJECTS. That may be
 *  removed when it's clear this is no longer needed.
 **/
class BLLoad {
	G4String libname;
	void *handle;
	static std::vector<BLLoad*> list;
	// Constructor.
	BLLoad(G4String name, void *p);
public:
	/// load() is a static function to load a shared object.
	/// It returns a pointer to a BLLoad object, null if error.
	/// If no '/' or '\' is present in the name, it will be
	/// searched for in the usual way.
	/// Normally no extension is used (.so, .dylib, or .dll);
	/// this function then appends those extensions in order until a
	/// file is found.
	/// Loading an already loaded shared object does nothing (no error),
	/// returning a new BLLoad object with the same handle.
	static BLLoad *load(const G4String &name) {
		G4String n = name;
		void *p = load1(n);
		if(p) return new BLLoad(n,p);
		n = name + ".so";
		p = load1(n);
		if(p) return new BLLoad(n,p);
		n = name + ".dylib";
		p = load1(n);
		if(p) return new BLLoad(n,p);
		n = name + ".dll";
		p = load1(n);
		if(p) return new BLLoad(n,p);

		return 0;
	}

	/// errorString() returns a string describing the previous 
	/// error. Static.
	static const char *errorString() {
#ifdef WIN32
		return "Unsupported OS";
#endif
#ifdef __CYGWIN__
		return "Unsupported OS";
#endif
#if defined(__linux) || defined(__APPLE__)
		return dlerror(); 
#endif
	}

	/// getList() returns the list of BLLoad pointers. Static.
	static std::vector<BLLoad*> getList() { return list; }

	/// getName() returns the library name of this DLLoad object.
	G4String getName() { return libname; }

public:
	/// getDataSymbol() returns a pointer to a DATA symbol in the shared
	/// object. Returns null if not found.
	/// The caller must cast the return value to the appropriate type.
	void *getDataSymbol(G4String sym) 
		{ return dlsym(handle,sym.c_str()); }

	/// getFunctionSymbol() returns a pointer to a FUNCTION symbol in the
	/// shared object. Returns null if not found.
	/// The caller must cast the return value to the appropriate type.
	VoidFunctionPointer getFunctionSymbol(G4String sym) {
		union { void *a; VoidFunctionPointer b; };
		a = dlsym(handle,sym.c_str());
		return b;
	}

	// Destructor.
	~BLLoad() { if(handle != 0) dlclose(handle); handle=0; }
private: 
	/// load1() is used internally by load()
	static void *load1(G4String n) {

#ifdef WIN32
		return 0;
#endif
#ifdef __CYGWIN__
		return 0;
#endif
#if defined(__linux) || defined(__APPLE__)
		return dlopen(n.c_str(),RTLD_NOW);
#endif
	}
};

#endif // BLLOAD_HH
