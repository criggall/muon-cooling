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
	    // There is no simple sleep() in Windows, so spin for ~ 5 seconds
	    // to give the user time to read the error message.
	    time_t end = time(0)+5;
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


/**	findGeant4Data() will find the Geant4 data sets and put them into
 *	the environment.
 **/
bool findGeant4Data()
{
	// check for the always-required datasets in the environment
	if(getenv("G4LEVELGAMMADATA") != 0 &&
	   getenv("G4LEDATA") != 0 && 
	   getenv("G4SAIDXSDATA") != 0 && 
	   getenv("G4RADIOACTIVEDATA") != 0 && 
	   getenv("G4ENSDFSTATEDATA") != 0 && 
	   getenv("G4NEUTRONXSDATA") != 0) return true;

	std::vector<std::string> fileList; // list of files and dirs
	std::string dataDir;		// full path of data directory
	// errorMsg is printed only if search fails; keeps log of attempts
	std::string errorMsg("Search for Geant4Data failed:\n");

	// First, look for the file .data in g4bl_dir. If present
	// it contains the path to the Geant4Data directory.
	if(fileList.size() == 0) {
		std::string dotData = g4bl_dir + "/.data";
		FILE *f = fopen(dotData.c_str(),"r");
		if(f != 0) {
			char line[1024];
			if(fgets(line,sizeof(line),f)) {
			    dataDir = line;
			    std::string::size_type p=dataDir.find_first_of("\r\n");
			    if(p != dataDir.npos) dataDir.erase(p,dataDir.npos);
			}
			fclose(f);
		}
		if(dataDir.size() > 0) fileList = dirList(dataDir);
		if(fileList.size() == 0)
			errorMsg += "No file '" + dotData + "'\n";
	}

	// Second, look for Geant4Data in HOME
	if(fileList.size() == 0 && getenv("HOME") != 0) {
		dataDir = getenv("HOME");
		dataDir += "/Geant4Data";
		fileList = dirList(dataDir);
		if(fileList.size() == 0)
			errorMsg += "No directory '" + dataDir + "'\n";
	}

#ifdef WIN32
	// Third, look in C:\Geant4Data
	if(fileList.size() == 0) {
		dataDir = "C:\\Geant4Data";
		fileList = dirList(dataDir);
		if(fileList.size() == 0)
			errorMsg += "No directory '" + dataDir + "'\n";
	}
#endif // WIN32

	// loop over fileList, ading entries into the environment
	int n_req=0;
#define SET(N,V) { char tmp[1024]; sprintf(tmp,"%s=%s",N,V); putenv(strdup(tmp)); printf("%s=%s\n",N,V); }
	for(unsigned i=0; i<fileList.size(); ++i) {
		std::string fn = fileList[i];
		std::string path = dataDir + "/" + fn;
		if(fn.find("G4ABLA") == 0) {
			SET("G4ABLADATA",path.c_str())
		} else if(fn.find("G4EMLOW") == 0) {
			SET("G4LEDATA",path.c_str())
			++n_req;
		} else if(fn.find("G4ENSDFSTATE") == 0) {
			SET("G4ENSDFSTATEDATA",path.c_str())
			++n_req;
		} else if(fn.find("G4NDL") == 0) {
			SET("G4NEUTRONHPDATA",path.c_str())
		} else if(fn.find("G4NEUTRONXS") == 0) {
			SET("G4NEUTRONXSDATA",path.c_str())
			++n_req;
		} else if(fn.find("G4PII") == 0) {
			SET("G4PIIDATA",path.c_str())
		} else if(fn.find("G4SAIDDATA") == 0) {
			SET("G4SAIDXSDATA",path.c_str())
			++n_req;
		} else if(fn.find("PhotonEvaporation") == 0) {
			SET("G4LEVELGAMMADATA",path.c_str())
			++n_req;
		} else if(fn.find("RadioactiveDecay") == 0) {
			SET("G4RADIOACTIVEDATA",path.c_str())
			++n_req;
		} else if(fn.find("RealSurface") == 0) {
			SET("G4REALSURFACEDATA",path.c_str())
		}
	}
#undef SET
	// print error message, if no data found
	if(fileList.size() == 0)
		fprintf(stderr,"%s\n",errorMsg.c_str());

	// return true only if all always-required datasets are found.
	// (may be multiple versions of some)
	return n_req >= 6;
}

void runG4blData()
{
	std::string path = "\"" + g4bl_dir + "/bin/g4bldata" + "\"";
	printf("Running %s\n",path.c_str());
	system(path.c_str());
	printf("%s returns\n",path.c_str());
	// now that the data are probably there, try to find them again
	findGeant4Data();
}

