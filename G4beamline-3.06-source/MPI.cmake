#	MPI.cmake
#
#	Used ONLY if G4BL_MPI=ON.
#
#	MPI must be configured for every site that uses it, and G4beamline
#	must be built from source on the site.
#
#	Put this just ahead of the else():
#	elseif(${SITE_NAME} MATCHES XXX)   # XXX = the result of `uname -n`
#	    set(LIBS ${LIBS} mpi_cxx mpi)  # the libraries your MPI uses
#	    include_directories(/opt/openmpi-64/include)  # the MPI include dir
#	    link_directories(/opt/openmpi-64/lib)  # the directory for the libs
#
#	If any of these are in default places (e.g. /usr, /usr/local), or
#	known to the build system (e.g. at NERSC), you can omit them.
#	But you must have the elseif(...) to avoid the FATAL_ERROR.

site_name(SITE_NAME)

if(${SITE_NAME} MATCHES tjrob)		# Mac OS X 10.9.5 with OpenMPI
	set(LIBS ${LIBS} mpi)
	# installed in /usr/local by brew
elseif(${SITE_NAME} MATCHES edison)	# NERSC
	set(LIBS ${LIBS} -dynamic)	# compiler defaults to -static
	message(STATUS "Dynamic link -- copy the library .so-s into lib.")
elseif(${SITE_NAME} MATCHES cori)	# NERSC
	set(LIBS ${LIBS} -dynamic)	# compiler defaults to -static
	message(STATUS "Dynamic link -- copy the library .so-s into lib.")
else()
	message(FATAL_ERROR "MPI.cmake: MPI not configured for '${SITE_NAME}'")
endif()
