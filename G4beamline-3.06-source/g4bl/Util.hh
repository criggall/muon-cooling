//	Util.hh

#ifndef UTIL_HH
#define UTIL_HH

// these routines always use std::string
#include <string>
#include <vector>

/// Returns the full path of the G4beamline install directory.
/// prints a message and exits if it cannot be determined.
/// returns (const char *) because some callers want G4String, some want
/// std::string, and some want QString.
const char *getG4blDir();

/// dirList() returns a list of files and directories in dirPath.
std::vector<std::string> dirList(std::string dirPath);

#endif // UTIL_HH
