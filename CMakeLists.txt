#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 3.16...3.21)

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

project(IMBP)


option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  message("with ui_all vis_all")
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  message("Geant4 REQUIRED")
  find_package(Geant4 REQUIRED)
endif()

#----------------------------------------------------------------------------
# Find ROOT (required package)
#
find_package(ROOT REQUIRED)
if(NOT ROOT_FOUND)
  message(STATUS "G4 Examples: ROOT package not found. --> exMPI02 example disabled")
  return()
endif()
#message(":: ROOT_DEFINITIONS = " ${ROOT_DEFINITIONS})
#message(":: ROOT_CC_FLAGS = " ${ROOT_CC_FLAGS})
#message(":: ROOT_FOUND = " ${ROOT_FOUND})
#message(":: ROOT_USE_FILE = " ${ROOT_USE_FILE})
#message(":: ROOT_CXX_FLAGS = " ${ROOT_CXX_FLAGS})
#message(":: ROOT_INCLUDE_DIRS = " ${ROOT_INCLUDE_DIRS})
#message(":: ROOT_DIR = " ${ROOT_DIR})

#message(":: ROOT_LIBRARY_DIR = " ${ROOT_LIBRARY_DIR})
#include("${ROOT_DIR}/modules/RootNewMacros.cmake")


#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
#
include(${Geant4_USE_FILE})
include(${ROOT_USE_FILE})

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
#include_directories(${PROJECT_SOURCE_DIR}/include
#                    ${Geant4_INCLUDE_DIR}
#                    ${ROOT_INCLUDE_DIRS})
include_directories(${PROJECT_SOURCE_DIR}/include
                    ${Geant4_INCLUDE_DIR})
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)


add_executable(RunGun RunGun.cc ${sources} ${headers})
target_link_libraries(RunGun ${Geant4_LIBRARIES} )
#target_link_libraries(RunGun ${Geant4_LIBRARIES} ${ROOT_LIBRARIES} )

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build B1. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(EXAMPLEB1_SCRIPTS
  RunGun.in
  RunGun.out
  init_vis.mac
  run1.mac
  run2.mac
  vis.mac
  )

foreach(_script ${EXAMPLEB1_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# For internal Geant4 use - but has no effect if you build this
# example standalone
#
add_custom_target(B1 DEPENDS exampleB1)

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS RunGun DESTINATION bin)
