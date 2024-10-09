cmake_minimum_required(VERSION 2.8)
cmake_policy(SET CMP0043 NEW)
cmake_policy(SET CMP0020 NEW)
set(CMAKE_AUTOMOC ON)
set(CMAKE_INCLUDE_CURRENT_DIR ON)
find_package(Qt5Widgets REQUIRED)

if(G4BL_VISUAL)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DG4BL_VISUAL")
endif()

add_executable(g4blgui
	g4blgui.cc
	G4blGuiWindow.cc
	RecentFile.cc
	../g4bl/Util.cc
	../g4bl/Geant4Data.cc
	)

target_link_libraries(g4blgui ${LIBS})
qt5_use_modules(g4blgui Widgets)

install(TARGETS g4blgui RUNTIME DESTINATION bin)
