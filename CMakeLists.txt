cmake_minimum_required(VERSION 3.10 FATAL_ERROR)

project(vmm-sdat LANGUAGES CXX)

# -----------------------------------------------------------------------------
# Build settings
# -----------------------------------------------------------------------------
set(CMAKE_BUILD_TYPE RelWithDebInfo)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# -----------------------------------------------------------------------------
# Extra CMake modules
# -----------------------------------------------------------------------------
set(EXTRA_MODULES_DIR ${CMAKE_CURRENT_LIST_DIR}/cmake)
list(APPEND CMAKE_MODULE_PATH ${EXTRA_MODULES_DIR}/modules)

# -----------------------------------------------------------------------------
# Conan (unchanged)
# -----------------------------------------------------------------------------
set(CONAN_DISABLE_CHECK_COMPILER TRUE)
set(CONAN_PROFILE "default" CACHE STRING "Name of conan profile to use")
set(CONAN "AUTO" CACHE STRING "AUTO | MANUAL | DISABLE")

if(${CONAN} MATCHES "AUTO")
  include(${EXTRA_MODULES_DIR}/modules/conan.cmake)
  conan_cmake_run(
      CONANFILE conanfile.txt
      PROFILE ${CONAN_PROFILE}
      BASIC_SETUP NO_OUTPUT_DIRS KEEP_RPATHS NO_IMPORTS
      BUILD_TYPE "None"
      BUILD outdated
  )
elseif(${CONAN} MATCHES "MANUAL")
  if(EXISTS "${CMAKE_BINARY_DIR}/conanbuildinfo.cmake")
    include("${CMAKE_BINARY_DIR}/conanbuildinfo.cmake")
    conan_basic_setup(NO_OUTPUT_DIRS KEEP_RPATHS)
  else()
    message(FATAL_ERROR "CONAN=MANUAL but conanbuildinfo.cmake not found")
  endif()
elseif(NOT ${CONAN} MATCHES "DISABLE")
  message(FATAL_ERROR "Invalid CONAN option: ${CONAN}")
endif()

# -----------------------------------------------------------------------------
# ROOT (force correct discovery for old installations)
# -----------------------------------------------------------------------------
if(DEFINED ENV{ROOTSYS})
  list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
endif()

find_package(ROOT REQUIRED)

message(STATUS "ROOT version: ${ROOT_VERSION}")
message(STATUS "ROOT include dirs: ${ROOT_INCLUDE_DIRS}")
message(STATUS "ROOT libraries: ${ROOT_LIBRARIES}")

include(${ROOT_USE_FILE})

# -----------------------------------------------------------------------------
# Include paths
# -----------------------------------------------------------------------------
include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}/src
  ${CMAKE_CURRENT_SOURCE_DIR}/analysis
  ${ROOT_INCLUDE_DIRS}
  ${VERSION_INCLUDE_DIR}
)

# -----------------------------------------------------------------------------
# ROOT dictionary
# -----------------------------------------------------------------------------
ROOT_GENERATE_DICTIONARY(
  vmm-sdatDict
  ${CMAKE_CURRENT_SOURCE_DIR}/src/DataStructures.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/RootFile.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/convertLinkDef.h
  MODULE vmm-sdat
)

# -----------------------------------------------------------------------------
# Shared library
# -----------------------------------------------------------------------------
add_library(vmm-sdat SHARED
  src/parser/ParserSRS.cpp
  src/parser/CalibrationFile.cpp
  src/parser/ReaderPcap.cpp
  src/Clusterer.cpp
  src/log.cpp
  src/RootFile.cpp
  src/Configuration.cpp
  src/Statistics.cpp
  vmm-sdatDict.cxx
)

target_link_libraries(vmm-sdat
  PRIVATE
    ${ROOT_LIBRARIES}
    pcap
)

target_compile_definitions(vmm-sdat PUBLIC __FAVOR_BSD)

# -----------------------------------------------------------------------------
# Executables
# -----------------------------------------------------------------------------
add_executable(convertFile src/convertFile.cpp)
add_executable(accessTree analysis/examples_read_data/accessTree.cpp)
add_executable(normalizeHisto analysis/examples_read_data/normalizeHisto.cpp)

target_link_libraries(convertFile PRIVATE vmm-sdat)
target_link_libraries(accessTree PRIVATE vmm-sdat)
target_link_libraries(normalizeHisto PRIVATE vmm-sdat)

# -----------------------------------------------------------------------------
# Final report
# -----------------------------------------------------------------------------
message(STATUS "#######################################################")
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "#######################################################")
