//	g4blmpi.cc
/**
 *	g4blmpi runs G4beamline in MPI mode. See USAGE below.
 *
 *	MPI is not supported on Windows, and on Mac OS X and Linux G4beamline
 *	must be built from source at the destination site:
 *		configure the site's MPI installation in MPI.cmake
 *		cmake -DG4BL_MPI=ON ...
 *		make install
 *
 *	Note that with an MPI build, users don't need to use g4blmpi,
 *	they can do "srun -n 1234 g4bl ..." directly. g4blmpi merely
 *	figures out which program (mpirun,aprun,srun) to use, and then
 *	runs it.
 **/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#if defined(WIN32) || ! defined(G4BL_MPI)

int main(int argc, char *argv[])
{
	if(argc > 1 && strcmp(argv[1],"--check") == 0) {
		exit(1);
	}

	fprintf(stderr,"g4blmpi: MPI not available\n");
	exit(1);
}

#else // !WIN32 and G4BL_MPI

#include <unistd.h>
#include <string>
#include <vector>
#include <sys/stat.h>

#include "Util.hh"

std::string g4bl_dir;

// potential MPI run program names
static const char *progs[] = { 
	"mpirun",	// OpenMPI (e.g. on Mac OS X, but also Linux)
	"srun",		// Newer supercomputers at NERSC (Cori, new Edison)
	"aprun",	// Older supercomputers at NERSC (Franklin, old Edison)
};

std::string findInPath(std::string prog)
{
	std::string path(getenv("PATH"));
	while(path.size() > 0) {
		// get next dir
		std::string::size_type i = path.find(":");	// UNIX only
		if(i == path.npos) i = path.size();
		std::string dir = path.substr(0,i);
		if(i == 0) dir = ".";
		path.erase(0,i+1);
		// check if prog exists in dir
		std::string file = dir + "/" + prog;
		struct stat s;
		if(stat(file.c_str(),&s) != 0) continue; // no file
		return file;
	}

	return std::string();
}

int main(int argc, char *argv[])
{
	// determine how to run an MPI program
	std::string mpirun;
	for(int i=0; i<sizeof(progs)/sizeof(progs[0]); ++i) {
		mpirun = findInPath(progs[i]);
		if(mpirun.size() != 0) break;
	}

	// print usage if no arguments
	if(argc <= 1) {
		fprintf(stderr,"g4blmpi USAGE:\n"
		"    g4blmpi 24 input.file [params...]\n"
		"        Runs G4beamline in MPI mode with 24 ranks (23 workers).\n"
		"    g4blmpi --check\n"
		"        Exit code 0 = MPI available; 1 = MPI not available.\n"
		);
		if(mpirun.size() == 0) {
			fprintf(stderr,"NOTE: no mpirun program is found\n");
			fprintf(stderr,"Perhaps you need to add it to PATH\n");
		}
		exit(1);
	}

	// find g4bl executable
	g4bl_dir = getG4blDir();
	std::string g4bl = g4bl_dir + "/bin/g4bl";

	// check for --check
	if(strcmp(argv[1],"--check") == 0) {
		if(mpirun.size() > 0 && g4bl.size() > 0) exit(0);
		exit(1);
	}

	// check we found all the pieces
	if(mpirun.size() == 0 || g4bl.size() == 0) {
		fprintf(stderr,"g4blmpi: MPI not available\n");
		exit(1);
	}

	// check we have sufficient arguments
	if(argc < 3) {
		fprintf(stderr,"g4blmpi: invalid arguments\n");
		exit(1);
	}
	int nRanks = strtol(argv[1],0,0);
	if(nRanks < 2) {
		fprintf(stderr,"g4blmpi: invalid # ranks\n");
		exit(1);
	}

	// exec mpirun 
	std::vector<char *> args;
	args.push_back(strdup(mpirun.c_str()));
	args.push_back((char *)"-n");
	args.push_back(argv[1]);
	args.push_back(strdup(g4bl.c_str()));
	for(int i=2; i<argc; ++i)
		args.push_back(argv[i]);
	args.push_back(0);
	execv(mpirun.c_str(),&args[0]);
	fprintf(stderr,"g4blmpi: cannot exec '%s'\n",mpirun.c_str());
	exit(1);
}

#endif // WIN32
