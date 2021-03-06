# CMakeLists.txt for vmm-sdat program. It creates a library with dictionary and a main program
cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(vmm-sdat)

# Temporary hack allowing conditional compilation of efu code
add_definitions("-DVMM_SDAT")

set(CMAKE_BUILD_TYPE RelWithDebInfo)
set(CMAKE_CXX_STANDARD 14)
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -fsanitize=address -O0 -ggdb")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ggdb")


#=============================================================================
# Augment CMake with our custom scripts
#=============================================================================
set(EXTRA_MODULES_DIR ${CMAKE_CURRENT_LIST_DIR}/cmake)
list(APPEND CMAKE_MODULE_PATH ${EXTRA_MODULES_DIR}/modules)

#=============================================================================
# Versioning
#=============================================================================
include(${EXTRA_MODULES_DIR}/Versioning.cmake)

#set_version()
set(VERSION_INCLUDE_DIR ${CMAKE_BINARY_DIR}/version)
create_version_header(
  ${EXTRA_MODULES_DIR}/templates/version_num.h.in
  ${VERSION_INCLUDE_DIR}/common/version_num.h
)
create_version_header(
  ${EXTRA_MODULES_DIR}/templates/Version.h.in
  ${VERSION_INCLUDE_DIR}/common/Version.h
)

#=============================================================================
# Conan
#=============================================================================

SET(CONAN_PROFILE "default" CACHE STRING "Name of conan profile to use, uses default by default")
SET(CONAN "AUTO" CACHE STRING "conan options AUTO (conan must be in path), MANUAL (expects conanbuildinfo.cmake in build directory) or DISABLE")
if(${CONAN} MATCHES "AUTO")
  include(${EXTRA_MODULES_DIR}/modules/conan.cmake)
  conan_cmake_run(CONANFILE conanfile.txt
      PROFILE ${CONAN_PROFILE}
      BASIC_SETUP NO_OUTPUT_DIRS KEEP_RPATHS NO_IMPORTS
      BUILD_TYPE "None"
      BUILD outdated)
elseif(${CONAN} MATCHES "MANUAL")
  if(EXISTS "${CMAKE_BINARY_DIR}/conanbuildinfo.cmake")
    include("${CMAKE_BINARY_DIR}/conanbuildinfo.cmake")
    conan_basic_setup(NO_OUTPUT_DIRS KEEP_RPATHS)
  else()
    MESSAGE(FATAL_ERROR "CONAN set to MANUAL but no file named conanbuildinfo.cmake found in build directory")
  endif()
elseif(NOT ${CONAN} MATCHES "DISABLE")
  MESSAGE(FATAL_ERROR "Unrecognised option for CONAN (${CONAN}), use AUTO, MANUAL or DISABLE")
endif()



if(ENABLE_DTRACE)
    message(STATUS "DTrace enabled")
    add_definitions(-DTRC_LEVEL=TRC_L_DEB)
 endif()


 #
 # Specify location of EFU repository
 #
 if (NOT DEFINED EFU_DIR)
   set(EFU_DIR "$ENV{HOME}/essdaq/efu/event-formation-unit"  CACHE STRING "EFU directory")
 endif()

 if (NOT EXISTS ${EFU_DIR})
   message(FATAL_ERROR "efu directory ${EFU_DIR} does not exist. Use -DEFU_DIR=efupath")
 endif()



# You need to tell CMake where to find the ROOT installation. This can be done in a number of ways:
#   - ROOT built with classic configure/make use the provided $ROOTSYS/etc/cmake/FindROOT.cmake
#   - ROOT built with CMake. Add in CMAKE_PREFIX_PATH the installation prefix for ROOT
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})

#find_package(Boost REQUIRED COMPONENTS system)
find_package(HDF5 1.10 REQUIRED)
find_package(h5cpp REQUIRED)

#---Locate the ROOT package and defines a number of variables (e.g. ROOT_INCLUDE_DIRS)
find_package(ROOT REQUIRED)

#---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
include(${ROOT_USE_FILE})

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)
include_directories(${EFU_DIR}/src)
include_directories(${EFU_DIR}/src/modules)
include_directories(${EFU_DIR}/src/udpgenpcap)
include_directories(${VERSION_INCLUDE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/analysis)


ROOT_GENERATE_DICTIONARY(vmm-sdatDict ${EFU_DIR}/src/modules/gdgem/nmx/Readout.h ${EFU_DIR}/src/modules/gdgem/srs/ParserVMM3.h ${EFU_DIR}/src/modules/gdgem/srs/SRSTime.h ${EFU_DIR}/src/modules/gdgem/srs/CalibrationFile.h ${EFU_DIR}/src/udpgenpcap/ReaderPcap.h  ${CMAKE_CURRENT_SOURCE_DIR}/src/Clusterer.h ${CMAKE_CURRENT_SOURCE_DIR}/src/RootFile.h ${CMAKE_CURRENT_SOURCE_DIR}/src/DataStructures.h ${CMAKE_CURRENT_SOURCE_DIR}/src/Configuration.h ${CMAKE_CURRENT_SOURCE_DIR}/src/Statistics.h LINKDEF ${CMAKE_CURRENT_SOURCE_DIR}/src/convertLinkDef.h MODULE vmm-sdat)


#---Create a shared library with generated dictionary
add_library(vmm-sdat SHARED ${EFU_DIR}/src/modules/gdgem/srs/ParserVMM3.cpp ${EFU_DIR}/src/modules/gdgem/srs/SRSTime.cpp ${EFU_DIR}/src/modules/gdgem/srs/CalibrationFile.cpp ${EFU_DIR}/src/udpgenpcap/ReaderPcap.cpp ${CMAKE_CURRENT_SOURCE_DIR}/src/Clusterer.cpp ${CMAKE_CURRENT_SOURCE_DIR}/src/RootFile.cpp ${CMAKE_CURRENT_SOURCE_DIR}/src/Configuration.cpp ${CMAKE_CURRENT_SOURCE_DIR}/src/Statistics.cpp vmm-sdatDict.cxx)
target_link_libraries(vmm-sdat ${ROOT_LIBRARIES} ${Boost_SYSTEM_LIBRARY} ${HDF5_LIBRARIES} h5cpp fmt graylog_logger_static pcap)

#---Create  a main program using the library
add_executable(convertFile ${CMAKE_CURRENT_SOURCE_DIR}/src/convertFile.cpp)
add_executable(accessTree ${CMAKE_CURRENT_SOURCE_DIR}/analysis/accessTree.cpp)


target_link_libraries(convertFile vmm-sdat)
target_link_libraries(accessTree vmm-sdat)

target_compile_definitions(vmm-sdat PUBLIC __FAVOR_BSD)

# Final reporting
#
message(STATUS "#######################################################")
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "EFU directory: ${EFU_DIR}")
message(STATUS "#######################################################")
