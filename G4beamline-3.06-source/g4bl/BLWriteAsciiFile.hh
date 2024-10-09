//	BLWriteAsciiFile.hh
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

#ifndef BLASCIIFILE
#define BLASCIIFILE

#include <map>
#include <set>
#include <stdio.h>

#include "globals.hh"

/**	class BLWriteAsciiFile manages file-descriptors for writing ASCII files.
 *
 *	Multiple fopen()-s of a given filename will open it ONCE, and always
 *	return the same FILE *.
 *	Multiple fclose()-s of the FILE* are handled properly -- fhe N-th 
 *	one closes the file, where N is the number of fopen-s.
 *
 *	"-" means stdout (not closed).
 *
 *	The filename-s must match exactly as given, no conversion to absolute
 *	paths (or any other canonicalization) is performed.
 *
 *	Only does writing.
 **/
class BLWriteAsciiFile {
	static std::map<G4String,FILE*> fds;
	static std::map<G4String,int> num_opens;
public:
	/// fopen() will open a file for writing, returning the same FILE *
	/// for multiple opens.
	static FILE *fopen(G4String name) {
		if(name == "-") return stdout;
		if(fds.count(name) == 0) {
			FILE *f = ::fopen(name.c_str(),"w");
			if(f != 0)
				fds[name] = f;
		}
		++num_opens[name];
		return fds[name];
	}

	/// fclose() will close the FILE*, handling multiple calls.
	static void fclose(FILE *f) {
		if(f == stdout) return;
		std::map<G4String,FILE*>::iterator i;
		for(i=fds.begin(); i!=fds.end(); ++i) {
			if(i->second != f) continue;
			if(--num_opens[i->first] > 0) break;
			::fclose(f);
			fds.erase(i);
			break;
		}
	}

	/// closeAll() will close all files.
	static void closeAll() {
		std::map<G4String,FILE*>::iterator i;
		for(i=fds.begin(); i!=fds.end(); ++i) {
			::fclose(i->second);
		}
		fds.clear();
	}
};

#endif // BLASCIIFILE
