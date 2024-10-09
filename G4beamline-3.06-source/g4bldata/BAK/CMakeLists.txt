cmake_minimum_required(VERSION 2.8)
cmake_policy(SET CMP0043 NEW)
set(CMAKE_AUTOMOC ON)
set(CMAKE_INCLUDE_CURRENT_DIR ON)
find_package(Qt5Widgets REQUIRED)
find_package(Qt5Network REQUIRED)

add_executable(g4bldata
	g4bldata.cc
	G4blDataWindow.cc
	Reader.cc
	../g4bl/Util.cc
	../g4bl/Geant4Data.cc
	)

if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
	target_link_libraries(g4bldata "-Wl,-rpath,$ORIGIN/../lib")
	set(CMAKE_SKIP_RPATH ON)
elseif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	target_link_libraries(g4bldata "-Wl,-rpath,@executable_path/../lib")
	set(CMAKE_SKIP_RPATH ON)
endif()
target_link_libraries(g4bldata ${LIBS})
qt5_use_modules(g4bldata Widgets Network)
install(TARGETS g4bldata RUNTIME DESTINATION bin)

###
### Windows needs bsdtar from GnuWin32
###
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	if(NOT EXISTS "C:/Program Files (x86)/GnuWin32/bin/bsdtar.exe")
		message( FATAL_ERROR "bsdtar.exe not found")
	endif()
	install(PROGRAMS "C:/Program Files (x86)/GnuWin32/bin/bsdtar.exe"
		"C:/Program Files (x86)/GnuWin32/bin/bzip2.dll"
		"C:/Program Files (x86)/GnuWin32/bin/libarchive2.dll"
		"C:/Program Files (x86)/GnuWin32/bin/zlib1.dll"
		DESTINATION bin)
endif()
