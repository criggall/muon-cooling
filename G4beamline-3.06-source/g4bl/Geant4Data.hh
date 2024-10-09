//	Geant4Data.hh - list of Geant4 datasets
#ifndef GEANT4DATA_HH
#define GEANT4DATA_HH

#include <vector>
#include <string>

/**
 *	Geant4Data represents a Geant4 dataset.
 *	   list() returns all datasets for the version of Geant4 in use.
 *	   setup() handles only the datasets that are present in directory().
 *	It does NOT download or unpack them (see g4bldata for that).
 **/
struct Geant4Data {
	std::string env;	// name of the environment variable
	std::string name;	// name of the dataset (incl version)
	std::string dir;	// directory name under directory()
	long dl_size;		// size of the download (bytes)
	std::string disk_size;	// size on the disk "xxM"
	bool required;
	static std::string baseURL; // base URL at CERN for downloads
				    // download baseURL+"/"+name+".tar.gz"
				    // into directory(), then un-tar it to 
				    // create directory()/dir.
public:
	/// sets up the environment variables for all datasets that are found.
	/// returns false if some required dataset is not found.
	/// If overwriteEnv is true, environment variables are overwritten;
	/// if false, the values of variables for datasets are preserved.
	static bool setup(bool overwriteEnv=true);

	/// returns the directory containing the Geant4 datasets.
	/// if check is true, returns an empty string if not found.
	/// if check is false, returns the default directory name.
	/// Returns:
	///    the contents of getG4blDir()/.data, if present, or
	///    $HOME/Geant4Data
	/// On Windows, if $HOME/Geant4Data is not present then
	/// C:\Geant4Data is returned.
	static std::string directory(bool check=true);

	/// returns the list of all Geant4 datasets for the version of Geant4
	/// being used. Does not check if datasets are present.
	static std::vector<Geant4Data> list();

	/// prints the environment variables defining datasets.
	static void printEnv();

	/// Runs the g4bldata program
	static void runG4bldata();
};

#endif // GEANT4DATA_HH
