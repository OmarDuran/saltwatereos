# 0. CMake Minimum version
cmake_minimum_required(VERSION 3.3...3.12 FATAL_ERROR)

# Detect native architecture to avoid mismatch with TBB/system libraries
if(NOT CMAKE_OSX_ARCHITECTURES)
  execute_process(COMMAND uname -m OUTPUT_VARIABLE HOST_ARCH OUTPUT_STRIP_TRAILING_WHITESPACE)
  set(CMAKE_OSX_ARCHITECTURES "${HOST_ARCH}" CACHE STRING "Build architecture" FORCE)
  message(STATUS "Setting build architecture to native: ${HOST_ARCH}")
endif()

set(PROGRAM_NAME swEOS)

project(SaltWaterEOS LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(SWEOS_LIB_PATH "../Library/lib")
set(SWEOS_LIB_NAME "eosH2ONaCl")
# options
option(USE_OMP "Enable OpenMP to parallel computing" OFF)
option(USE_TBB "Enable Intel TBB for parallel computing" ON)

set(SWEOS_MAJOR_VERSION 1)
set(SWEOS_MINOR_VERSION 7)
set(SWEOS_PATCH_VERSION 0)
set(SWEOS_EXTRA_VERSION "")
set(SWEOS_EXTRA_VERSION_TEXI "${SWEOS_EXTRA_VERSION}")
if(NOT SWEOS_RELEASE)
  find_package(Git)
  if(GIT_FOUND)
    execute_process(COMMAND ${GIT_EXECUTABLE} log -1 --format=%h
                    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR} ERROR_QUIET
                    OUTPUT_VARIABLE GIT_COMMIT_HASH
                    OUTPUT_STRIP_TRAILING_WHITESPACE)
  endif()
  if(GIT_COMMIT_HASH)
    set(SWEOS_EXTRA_VERSION "${SWEOS_EXTRA_VERSION}-git-${GIT_COMMIT_HASH}")
  endif()
  set(SWEOS_EXTRA_VERSION_TEXI "${SWEOS_EXTRA_VERSION_TEXI} (development version)")
endif()
set(SWEOS_VERSION "${SWEOS_MAJOR_VERSION}.${SWEOS_MINOR_VERSION}")
set(SWEOS_VERSION "${SWEOS_VERSION}.${SWEOS_PATCH_VERSION}${SWEOS_EXTRA_VERSION}")
set(SWEOS_SHORT_LICENSE "GNU General Public License")

if(NOT DATE)
  string(TIMESTAMP DATE "%Y-%m-%d")
endif()
set(SWEOS_DATE "${DATE}")

if(NOT SWEOS_HOST)
  execute_process(COMMAND hostname OUTPUT_VARIABLE HOSTNAME
                  OUTPUT_STRIP_TRAILING_WHITESPACE)
  if(NOT HOSTNAME)
    set(HOSTNAME "unknown")
  endif()
  set(SWEOS_HOST "${HOSTNAME}")
endif()

if(NOT SWEOS_PACKAGER)
  execute_process(COMMAND whoami OUTPUT_VARIABLE PACKAGER
                  OUTPUT_STRIP_TRAILING_WHITESPACE)
  if(NOT PACKAGER)
    set(PACKAGER "unknown")
  endif()
  string(REPLACE "\\" " " PACKAGER ${PACKAGER})
  set(SWEOS_PACKAGER "${PACKAGER}")
endif()

# Using OpenMP
if(USE_OMP)
  include_directories(/usr/local/include /opt/local/include)
  find_package(OpenMP REQUIRED)
  if(APPLE)
      set(CMAKE_C_COMPILER clang)
      set(CMAKE_CXX_COMPILER clang++)

      if(CMAKE_C_COMPILER_ID MATCHES "Clang\$")
          set(OpenMP_C_FLAGS "-Xpreprocessor -fopenmp")
          set(OpenMP_C_LIB_NAMES "omp")
          set(OpenMP_omp_LIBRARY omp)
      endif()

      if(CMAKE_CXX_COMPILER_ID MATCHES "Clang\$")
          set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp")
          set(OpenMP_CXX_LIB_NAMES "omp")
          set(OpenMP_omp_LIBRARY omp)
      endif()
  endif()
  if(OPENMP_FOUND)
      message("OPENMP FOUND")
      set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
      set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
  endif()
  #windows vs2017 doesn't need opm library link
  set(LIBomp "")
  if (APPLE)
    set(LIBomp "omp")
  endif()
  # add MACRO definition
  add_definitions("-DUSE_OMP=1")
  link_directories("/usr/local/lib" "/opt/local/lib")
endif()

# Using TBB
set(LIBtbb "")
if(USE_TBB)
  find_package(TBB REQUIRED)
  if(TBB_FOUND)
    message(STATUS "TBB FOUND: ${TBB_VERSION}")
    add_definitions("-DUSE_TBB=1")
    set(LIBtbb TBB::tbb)
  else()
    message(FATAL_ERROR "TBB not found. Please install TBB (e.g. brew install tbb)")
  endif()
endif()


#2. Source code dir

# add H2ONaCl source and library
include_directories(../Library/include)
# Use absolute path and find the library directly to avoid Xcode appending config subdirectory
find_library(SWEOS_LIB ${SWEOS_LIB_NAME} PATHS "${CMAKE_CURRENT_SOURCE_DIR}/../Library/lib" NO_DEFAULT_PATH)
if(NOT SWEOS_LIB)
  message(FATAL_ERROR "Library ${SWEOS_LIB_NAME} not found in ${CMAKE_CURRENT_SOURCE_DIR}/../Library/lib")
endif()
message(STATUS "Found ${SWEOS_LIB_NAME}: ${SWEOS_LIB}")
aux_source_directory(. SRC_MAIN)
set( Srcs_bash ./SWEOSbash.cpp ./MultiProgressBar.cpp ./getopt.cpp)
set( Hdrs_bash ./SWEOSbash.h ./MultiProgressBar.h ./getopt.h)

# add Executable target
add_executable(${PROGRAM_NAME}  ${Srcs_bash}  main.cpp)

target_link_libraries(${PROGRAM_NAME} ${SWEOS_LIB} ${LIBomp} ${LIBtbb})

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/../swEOSVersion.h.in
               ${CMAKE_CURRENT_SOURCE_DIR}/swEOSVersion.h)

message(STATUS "")
message(STATUS "SWEOS ${SWEOS_VERSION} has been configured for ${SWEOS_OS}")
message(STATUS "")
message(STATUS " * Build type: " ${CMAKE_BUILD_TYPE})
message(STATUS " * C compiler: " ${CMAKE_C_COMPILER})
message(STATUS " * C++ compiler: " ${CMAKE_CXX_COMPILER})
message(STATUS " * Install prefix: " ${CMAKE_INSTALL_PREFIX})
message(STATUS "")
