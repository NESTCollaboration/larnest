
# Set minimum CMake version required and prevent building in source
cmake_minimum_required(VERSION 3.8...3.16 FATAL_ERROR)

###############        Get version from git        #####################
execute_process(COMMAND git describe --tag
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        OUTPUT_VARIABLE PROJECT_GIT_VERSION
        OUTPUT_STRIP_TRAILING_WHITESPACE)

function(get_versions versionString version)
    if ("${versionString}" STREQUAL "")
        set(version "0.0.0" PARENT_SCOPE)
        return()
    endif ()

    string(REGEX REPLACE "v([0-9]*)([.][0-9]*[.][0-9]*-?.*)$" "\\1" numbers ${versionString})
    set(major ${numbers})
    string(REGEX REPLACE "v([0-9]*[.])([0-9]*)([.][0-9]*-?.*)$" "\\2" numbers ${versionString})
    set(minor ${numbers})
    string(REGEX REPLACE "v([0-9]*[.][0-9]*[.])([0-9]*)(-?.*)$" "\\2" numbers ${versionString})
    set(patch ${numbers})
    set(version "${major}.${minor}.${patch}" PARENT_SCOPE)
endfunction()

get_versions("${PROJECT_GIT_VERSION}" version)
set(PROJECT_GIT_VERSION ${version})

###############    project name and version        #####################
project(LArNEST VERSION ${PROJECT_GIT_VERSION} LANGUAGES CXX)

option(G4 "Build integration with Geant4" OFF)
option(BUILD_ROOT "Build ROOT tools" OFF)
option(BUILD_BENCHMARKS "Build LArNEST benchmarks" ON)
option(BUILD_GARFIELD "Build integration with Garfield++" OFF)

# enable position independent code so shared libraries can link against 
# the static version of LArNEST::Core 
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

include(GNUInstallDirs)
# include(CMakePackageConfigHelpers)

set(LArNEST_CORE_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/src/LArNEST.cxx
  ${CMAKE_CURRENT_SOURCE_DIR}/src/Logger.cxx
  ${CMAKE_CURRENT_SOURCE_DIR}/src/RandomGen.cxx
)
set(LArNEST_CORE_HEADERS
  ${CMAKE_CURRENT_SOURCE_DIR}/include/LArNEST.h
  ${CMAKE_CURRENT_SOURCE_DIR}/include/LArParameters.h
  ${CMAKE_CURRENT_SOURCE_DIR}/include/Logger.h
  ${CMAKE_CURRENT_SOURCE_DIR}/include/RandomGen.h
)

add_library(LArNESTCore ${LArNEST_CORE_SOURCES} ${LArNEST_CORE_HEADERS})
add_library(LArNEST::Core ALIAS LArNESTCore)

list(APPEND LArNEST_TARGETS "LArNESTCore")
list(APPEND LArNEST_LIBRARY_TARGETS "LArNEST::Core")

target_include_directories(
  LArNESTCore PUBLIC
  $<INSTALL_INTERFACE:include>
  $<INSTALL_INTERFACE:include/LArNEST>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/LArNEST>
)

target_compile_features(LArNESTCore PUBLIC cxx_std_17)

# add_executable(execLArNEST
#   ${CMAKE_CURRENT_SOURCE_DIR}/src/execLArNEST.cxx
#   ${CMAKE_CURRENT_SOURCE_DIR}/include/execLArNEST.h
# )
# target_link_libraries(execLArNEST PUBLIC LArNEST::Core)
# list(APPEND LArNEST_TARGETS "execLArNEST")

# ## scripts
# configure_file(examples/loopLArNEST.in loopLArNEST)
# configure_file(examples/220RnCalib.in 220RnCalib)

## Optional targets

if(BUILD_BENCHMARKS)
  add_subdirectory(benchmarks)
endif()

# if(BUILD_GARFIELD)
#   add_subdirectory(GarfieldppIntegration)
# endif()

# if(G4)
#     add_subdirectory(G4integration)
#     list(APPEND LArNEST_TARGETS "LArNESTG4")
#     set(LArNEST_LIBRARY_TARGETS "${LArNEST_LIBRARY_TARGETS} LArNEST::LArNESTG4")
# endif()

# installation, Targets.cmake and Config.cmake
message(STATUS "LArNEST_LIBRARY_TARGETS: ${LArNEST_LIBRARY_TARGETS}")
set(INSTALL_CONFIGDIR ${CMAKE_INSTALL_LIBDIR}/cmake/LArNEST)

# write_basic_package_version_file(
#     LArNESTConfigVersion.cmake
#     VERSION ${PACKAGE_VERSION}
#     COMPATIBILITY ExactVersion
# )

# configure_package_config_file(${CMAKE_CURRENT_LIST_DIR}/cmake/LArNESTConfig.cmake.in
#     ${CMAKE_CURRENT_BINARY_DIR}/LArNESTConfig.cmake
#     INSTALL_DESTINATION ${INSTALL_CONFIGDIR}
# )

install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/LArNESTConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/LArNESTConfigVersion.cmake
    DESTINATION ${INSTALL_CONFIGDIR}
)

message(STATUS "CMAKE_INSTALL_LIBDIR: ${CMAKE_INSTALL_LIBDIR}")
message(STATUS "CMAKE_INSTALL_LIBDIR: ${CMAKE_INSTALL_BINDIR}")
message(STATUS "CMAKE_INSTALL_LIBDIR: ${CMAKE_INSTALL_INCLUDEDIR}")
install(
  TARGETS ${LArNEST_TARGETS}
  EXPORT LArNEST-Targets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(
  EXPORT LArNEST-Targets
  FILE LArNESTTargets.cmake
  NAMESPACE LArNEST::
  DESTINATION ${INSTALL_CONFIGDIR}
)

## install scripts

install(FILES
#   ${CMAKE_CURRENT_BINARY_DIR}/loopLArNEST
#   ${CMAKE_CURRENT_BINARY_DIR}/220RnCalib
  PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
  DESTINATION ${CMAKE_INSTALL_BINDIR}
)

##############################################

export(
  EXPORT LArNEST-Targets
  FILE ${CMAKE_CURRENT_BINARY_DIR}/LArNESTTargets.cmake
  NAMESPACE LArNEST::
)

# Register package in user's package registry
set(CMAKE_EXPORT_PACKAGE_REGISTRY ON)
export(PACKAGE LArNEST)
