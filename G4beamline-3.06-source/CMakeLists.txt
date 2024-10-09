### RPATH trial
#	CMakeLists.txt for G4beamline
#
#	Options:
#		G4BL_VISUAL	enable visualization, including OGLSQt.
#				Geant4 must be built with visualization
#				and QT; G4BL_GUI must be ON.
#		G4BL_GUI	enables the GUI programs. Requires Qt 5.X,
#				which cmake must be able to find. Can be OFF
#				only on Linux (e.g. for a supercomputer).
#		G4BL_ROOT	enables using Root. ROOT_DIR must be valid
#				(defaults to $ENV{ROOTSYS}).
#		G4BL_GSL	enables commands that use the GSL library.
#				GSL_DIR must be valid (defaults to 
#				$ENV{GSL_DIR}).
#		G4BL_FFTW	enables commands that use the FFTW library.
#				FFTW_DIR must be valid (defaults to 
#				$ENV{FFTW_DIR}).
#		G4BL_MPI	enables MPI. The file MPI.cmake must have a
#				configuration for MPI specific to the system
#				on which the build is running.
#	These options are defined in this file, but their corresponding
#	_DIR variable is set where they are used.
#
#	Geant4 is found in GEANT4_DIR, which defaults to $ENV{GEANT4_DIR}.
#
#	NOTE: the install directory MUST be the same as CMAKE_BINARY_DIR.
#	The programs are run from there. (This is set below.)
#
#	NOTE: "make install" on Mac OS X constructs Help.html, so it MUST
#	be run before any other OS. It also constructs the source distribution.
#
#	NOTE: the finalize directory is used for CMake sequencing -- it comes
#	last. finalize/CMakeLists.txt constructs the appropriate distribution
#	in "make install".
#
#	NOTE: On Windows the GnuWin32 tools must be installed, so g4bldata
#	can include its bsdtar in the install package.
#
cmake_minimum_required(VERSION 2.8)
cmake_policy(SET CMP0043 NEW)
cmake_policy(SET CMP0020 NEW)
project(G4beamline)
###
### This is the only definition of G4BLVERSION
### (Everybody else gets it from here. Some use grep+sed; don't change syntax.)
###
set(G4BLVERSION 3.06)
message(STATUS "G4BLVERSION=${G4BLVERSION}")
 
# install must be the same as $PWD (the binary dir)
set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}" CACHE PATH "Install dir" FORCE)

# only support Release builds
set(CMAKE_CONFIGURATION_TYPES "Release" CACHE STRING "Release only" FORCE)
set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Release only" FORCE)
###set(CMAKE_CONFIGURATION_TYPES "Debug" CACHE STRING "Debug only" FORCE)
###set(CMAKE_BUILD_TYPE "Debug" CACHE STRING "Debug only" FORCE)

# Configuration options
message(STATUS "")
set(G4BL_VISUAL ON CACHE BOOL "Build with visualization (incl. OGLSQt).")
message(STATUS "G4BL_VISUAL=${G4BL_VISUAL}")
set(G4BL_GUI ON CACHE BOOL "Build the GUI programs (g4blgui, g4bldata).")
message(STATUS "G4BL_GUI=${G4BL_GUI}")
set(G4BL_ROOT ON CACHE BOOL "Build with Root.")
message(STATUS "G4BL_ROOT=${G4BL_ROOT}")
set(G4BL_GSL ON CACHE BOOL "Use GSL Library.")
message(STATUS "G4BL_GSL=${G4BL_GSL}")
set(G4BL_FFTW ON CACHE BOOL "Use FFTW Library")
message(STATUS "G4BL_FFTW=${G4BL_FFTW}")
set(G4BL_MPI OFF CACHE BOOL "Build with MPI.")
message(STATUS "G4BL_MPI=${G4BL_MPI}")
set(GEANT4_DIR "$ENV{GEANT4_DIR}" CACHE PATH "Directory containing Geant4")
message(STATUS "GEANT4_DIR=${GEANT4_DIR}")

# check consistency of options
if(NOT G4BL_GUI)
	if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
		message(FATAL_ERROR "Mac OS X needs G4BL_GUI=ON")
	elseif(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
		message(STATUS "G4BL_GUI=OFF, so the Qt library is not used.")
	elseif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
		message(FATAL_ERROR "Windows needs G4BL_GUI=ON")
	endif()
	if(G4BL_VISUAL)
		message(FATAL_ERROR "G4BL_VISUAL needs G4BL_GUI=ON")
	endif()
endif()

if(G4BL_ROOT)
	set(ROOT_DIR "$ENV{ROOTSYS}" CACHE PATH "Directory containing Root")
	if(NOT IS_DIRECTORY ${ROOT_DIR})
		message(FATAL_ERROR "ROOT_DIR=${ROOT_DIR} does not exist")
	endif()
endif()

###
### RPATH
###
if(G4BL_ROOT)
	if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
		SET(CMAKE_INSTALL_RPATH 
			"@executable_path/../root/lib;@executable_path/../lib")
	elseif(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
		SET(CMAKE_INSTALL_RPATH "$ORIGIN/../root/lib;$ORIGIN/../lib")
	endif()
else()
	if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
		SET(CMAKE_INSTALL_RPATH "@executable_path/../lib")
	elseif(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
		SET(CMAKE_INSTALL_RPATH "$ORIGIN/../lib")
	endif()
endif()

###
### add the directory for each executable to be built
###
add_subdirectory(g4bl)
add_subdirectory(g4bltest)
add_subdirectory(g4blmpi)  # always built for "g4blmpi --check"
if(G4BL_GUI)
	add_subdirectory(g4blgui)
	add_subdirectory(g4bldata)
endif()

###
### installation
###
install(DIRECTORY Resources bin doc examples share test validation 
	DESTINATION . USE_SOURCE_PERMISSIONS)
install(FILES doc/README.txt DESTINATION .)

add_subdirectory(finalize)
