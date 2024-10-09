#	CMakeLists.txt for finalize - install code that runs after
#	all other directories are installed

###
### Qt5 plugins 
###
if(G4BL_GUI)
	execute_process(COMMAND qmake -query QT_INSTALL_PLUGINS 
		OUTPUT_VARIABLE plugins OUTPUT_STRIP_TRAILING_WHITESPACE)
	install(DIRECTORY ${plugins} 
		DESTINATION . USE_SOURCE_PERMISSIONS)
else()
	install(MAKE_DIRECTORY DESTINATION plugins)
endif()

###
### Darwin (Mac OS X)
###
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    # copy in all shared libraries and modify programs and libs to use them
    install(CODE "execute_process(COMMAND 
    			${CMAKE_SOURCE_DIR}/scripts/macosx_libraries .)")

    if(NOT G4BL_MPI)
	# Create Help.html in CMAKE_SOURCE_DIR, copy it into CMAKE_BINARY_DIR
	# Must be in finalize, as make-help runs g4bl.
	# (this way it is in CMAKE_SOURCE_DIR for other platforms, incl MPI)
	install(CODE "execute_process(COMMAND ${CMAKE_SOURCE_DIR}/scripts/make-help ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR})")
	# need to overwrite what has already been copied in
	install(CODE "execute_process(COMMAND cp ${CMAKE_SOURCE_DIR}/share/g4beamline/Help.html ${CMAKE_BINARY_DIR}/share/g4beamline)")

	# construct the Source Distribution
	install(CODE "execute_process(COMMAND ${CMAKE_SOURCE_DIR}/scripts/source-distrib ${CMAKE_SOURCE_DIR})")

	# construct the G4beamline distribution for Mac OS X
	install(CODE "execute_process(COMMAND ${CMAKE_SOURCE_DIR}/scripts/distrib-darwin ${G4BLVERSION})")
    else(NOT G4BL_MPI)
	message(STATUS "MPI Build, no app or distribution created.")
    endif(NOT G4BL_MPI)

###
### Linux
###
elseif(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
	# copy in all shared libraries 
	install(CODE "execute_process(COMMAND ${CMAKE_SOURCE_DIR}/scripts/linux_libraries --all)")

	# construct the G4beamline distribution for Linux
	install(CODE "execute_process(COMMAND ${CMAKE_SOURCE_DIR}/scripts/distrib-linux ${G4BLVERSION})")

###
### Windows
###
elseif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	message(STATUS "Need to construct Windows distribution")

###
### Unknown OS
###
else()
	message(FATAL_ERROR "unknown OS '${CMAKE_SYSTEM_NAME}'")
endif()
