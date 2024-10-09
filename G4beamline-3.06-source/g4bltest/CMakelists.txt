#	CMakeLists.txt for g4bltest

cmake_minimum_required(VERSION 2.8)
cmake_policy(SET CMP0043 NEW)
cmake_policy(SET CMP0020 NEW)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

if(G4BL_MPI)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DG4BL_MPI")
endif()
if(G4BL_ROOT)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DG4BL_ROOT")
endif()
if(G4BL_GSL)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DG4BL_GSL")
endif()
if(G4BL_FFTW)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DG4BL_FFTW")
endif()

include_directories(${PROJECT_SOURCE_DIR}/g4bl)
add_executable(g4bltest g4bltest.cc ../g4bl/Util.cc ../g4bl/Geant4Data.cc)

install(TARGETS g4bltest RUNTIME DESTINATION bin)
