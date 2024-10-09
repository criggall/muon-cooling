//	Util.cc

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <time.h>

#ifdef WIN32
#include <windows.h>
#include <tchar.h> 
#include <strsafe.h>
#include <process.h>
#undef min
#undef max
#else
#include <unistd.h>
#include <dirent.h>
#include <string.h>
#endif
#ifdef __APPLE__
#include <mach-o/dyld.h>
#endif

#include "Util.hh"
#include "Geant4Data.hh"

static std::string g4bl_dir; // there is a global in g4bl.cc

const char *getG4blDir()
{
	if(getenv("G4BL_DIR") != 0) {
		g4bl_dir = getenv("G4BL_DIR");
	} else {
#ifdef __linux__
          char arg1[32];
	  char exepath[PATH_MAX+1] = "";
	  sprintf( arg1, "/proc/%d/exe", getpid() );
	  ssize_t n=readlink( arg1, exepath, PATH_MAX );
	  if(n > 0) exepath[n] = '\0';
	  g4bl_dir = exepath;
	  size_t p = g4bl_dir.find_last_of("/");
	  if(p != g4bl_dir.npos) g4bl_dir.erase(p,g4bl_dir.npos);
	  p = g4bl_dir.find_last_of("/");
	  if(p != g4bl_dir.npos) g4bl_dir.erase(p,g4bl_dir.npos);
#endif
#ifdef __APPLE__
	  char path[4096]; uint32_t size = sizeof(path);
	  if(_NSGetExecutablePath(path, &size) == 0)
		g4bl_dir = path;
	  size_t p = g4bl_dir.find_last_of("/");
	  if(p != g4bl_dir.npos) g4bl_dir.erase(p,g4bl_dir.npos);
	  p = g4bl_dir.find_last_of("/");
	  if(p != g4bl_dir.npos) g4bl_dir.erase(p,g4bl_dir.npos);
#endif
#ifdef WIN32
	  char path[4096];
	  GetModuleFileName(NULL,path,sizeof(path));
	  g4bl_dir = path;
	  size_t p = g4bl_dir.find_last_of("/\\");
	  if(p != g4bl_dir.npos) g4bl_dir.erase(p,g4bl_dir.npos);
	  p = g4bl_dir.find_last_of("/\\");
	  if(p != g4bl_dir.npos) g4bl_dir.erase(p,g4bl_dir.npos);
#endif
	}

	// verify g4bl_dir by opening viewer.def
	std::string viewer = g4bl_dir + "/share/g4beamline/viewer.def";
	FILE *f = fopen(viewer.c_str(),"r");
	if(!f) { fprintf(stderr,
	      "************************************************************\n"
	      "***  Cannot determine the G4beamline install directory.  ***\n"
	      "***  g4bl must be run via an absolute path, or you can   ***\n"
	      "***  set G4BL_DIR in the environment.                    ***\n"
	      "************************************************************\n");
	    fprintf(stderr,"best-guess g4bl_dir='%s'\n",g4bl_dir.c_str());
	    // There is no simple sleep() in Windows, so spin for ~ 10 seconds
	    // to give the user time to read the error message.
	    time_t end = time(0)+10;
	    while(time(0) < end) ;
	    exit(1);
	}
	fclose(f);

	return g4bl_dir.c_str();
}

/**	dirList() returns a list of files and directories in a directory
 **/
std::vector<std::string> dirList(std::string dirPath)
{
	std::vector<std::string> ret;
#ifdef WIN32
	WIN32_FIND_DATA ffd;
	TCHAR szDir[MAX_PATH];
	HANDLE hFind = INVALID_HANDLE_VALUE;
	StringCchCopy(szDir, MAX_PATH, dirPath.c_str());
	StringCchCat(szDir, MAX_PATH, "/*");
	hFind = FindFirstFile(szDir, &ffd);
	if(hFind != INVALID_HANDLE_VALUE) {
		do {
			ret.push_back(ffd.cFileName);
		} while(FindNextFile(hFind, &ffd) != 0);
		FindClose(hFind);
	}
#else // !WIN32
	DIR *dir = opendir(dirPath.c_str());
	if(dir != 0) {
		for(;;) {
			struct dirent *d = readdir(dir);
			if(!d) break;
			ret.push_back(d->d_name);
		}
		closedir(dir);
	}
#endif // !WIN32
	return ret;
}

