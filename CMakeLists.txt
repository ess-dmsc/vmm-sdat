# CMakeLists.txt for vmm-hdf5-to-root program. It creates a library with dictionary and a main program
cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(vmm-hdf5-to-root)

set(CMAKE_BUILD_TYPE RelWithDebInfo)
set(CMAKE_CXX_STANDARD 11)

SET(HDF5_ROOT_DIR $ENV{HDF5_ROOT_DIR})

if(ENABLE_DTRACE)
    message(STATUS "DTrace enabled")
    add_definitions(-DTRC_LEVEL=TRC_L_DEB)   
 endif()  


# You need to tell CMake where to find the ROOT installation. This can be done in a number of ways:
#   - ROOT built with classic configure/make use the provided $ROOTSYS/etc/cmake/FindROOT.cmake
#   - ROOT built with CMake. Add in CMAKE_PREFIX_PATH the installation prefix for ROOT
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})

find_package(Boost REQUIRED COMPONENTS system)

find_package(HDF5 1.10 REQUIRED)

find_package(h5cpp REQUIRED)

#---Locate the ROOT package and defines a number of variables (e.g. ROOT_INCLUDE_DIRS)
find_package(ROOT REQUIRED)

#---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
include(${ROOT_USE_FILE})

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
ROOT_GENERATE_DICTIONARY(vmm-hdf5-to-rootDict Clusterer.h RootFile.h Readout.h DataStructures.h Configuration.h Statistics.h LINKDEF convertLinkDef.h MODULE vmm-hdf5-to-root)

#---Create a shared library with geneated dictionary
add_library(vmm-hdf5-to-root SHARED Clusterer.cpp RootFile.cpp Configuration.cpp Statistics.cpp vmm-hdf5-to-rootDict.cxx)
target_link_libraries(vmm-hdf5-to-root ${ROOT_LIBRARIES} ${Boost_SYSTEM_LIBRARY} ${HDF5_LIBRARIES} h5cpp)

#---Create  a main program using the library
add_executable(convertFile convertFile.cpp)
add_executable(accessTree accessTree.cpp)

target_link_libraries(convertFile vmm-hdf5-to-root)
target_link_libraries(accessTree vmm-hdf5-to-root)


