cmake_minimum_required(VERSION 2.8)
cmake_policy(SET CMP0043 NEW)
cmake_policy(SET CMP0020 NEW)

if(G4BL_MPI)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DG4BL_MPI")
endif()

include_directories(${PROJECT_SOURCE_DIR}/g4bl)
add_executable(g4blmpi g4blmpi.cc ../g4bl/Util.cc ../g4bl/Geant4Data.cc)

install(TARGETS g4blmpi RUNTIME DESTINATION bin)
