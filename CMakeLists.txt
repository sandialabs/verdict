# Copyright 2003,2006,2019 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
# Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software.

cmake_minimum_required(VERSION 3.16)

project(verdict VERSION 1.4.0)

# Includes
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

if (POLICY CMP0063)
  cmake_policy(SET CMP0063 NEW)
endif ()

set(verdict_VERSION_FLAT "${CMAKE_PROJECT_VERSION_MAJOR}${CMAKE_PROJECT_VERSION_MINOR}${CMAKE_PROJECT_VERSION_PATCH}")
set(verdict_VERSION "${CMAKE_PROJECT_VERSION}")

option(VERDICT_BUILD_DOC "Build the 2007 Verdict User Manual" OFF)
option(VERDICT_MANGLE "Mangle verdict names for inclusion in a larger library?" OFF)
if (VERDICT_MANGLE)
  set(VERDICT_MANGLE_PREFIX "verdict" CACHE STRING "The namespace enclosing verdict function names and classes.")
  mark_as_advanced(VERDICT_MANGLE_PREFIX)
endif ()
mark_as_advanced(VERDICT_MANGLE)

option(VERDICT_ENABLE_TESTING "Should tests of the VERDICT library be built?" ON)

set(verdict_HEADERS
  verdict.h
  VerdictVector.hpp
  verdict_defines.hpp)

set(verdict_SOURCES
  V_EdgeMetric.cpp
  V_GaussIntegration.cpp
  V_GaussIntegration.hpp
  V_HexMetric.cpp
  V_KnifeMetric.cpp
  V_PyramidMetric.cpp
  V_QuadMetric.cpp
  V_TetMetric.cpp
  V_TriMetric.cpp
  V_WedgeMetric.cpp
  VerdictVector.cpp)

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/verdict_config.h.in
  ${CMAKE_CURRENT_BINARY_DIR}/verdict_config.h
  @ONLY)

add_library(verdict ${verdict_SOURCES})
target_include_directories(verdict PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>)


# Setting the VERSION and SOVERSION of a library will include
# version information either in the library, or in the library
# name (depending on the platform). You may choose to exclude
# this information.
if (NOT VERDICT_NO_LIBRARY_VERSION)
  set_target_properties(verdict PROPERTIES
    VERSION "${verdict_VERSION}"
    SOVERSION "${CMAKE_PROJECT_VERSION_MAJOR}.${CMAKE_PROJECT_VERSION_MINOR}")
endif ()
if (NOT VERDICT_EXPORT_GROUP)
  set(VERDICT_EXPORT_GROUP VerdictExport)
endif ()

if (VERDICT_ENABLE_TESTING)
  enable_testing()
  if (NOT TARGET GTest::GTest)
    find_package(GTest REQUIRED)
  endif ()

  function(add_verdict_unittests DIR)
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/${DIR} ${CMAKE_CURRENT_BINARY_DIR}/${DIR})

    foreach (TEST ${ARGN})
      add_test(NAME UT-${TEST}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${DIR}
        COMMAND ${TEST})
    endforeach (TEST ${ARGN})
  endfunction()

  add_verdict_unittests(unittests unittests_verdict)
endif ()

if (VERDICT_BUILD_DOC)
  add_subdirectory(docs/VerdictUserManual2007)
endif ()

#########################
# Installation commands #
#########################
if (NOT verdict_INSTALL_DOC_DIR)
  set(verdict_INSTALL_DOC_DIR ${CMAKE_INSTALL_DOCDIR})
endif ()
if (NOT verdict_INSTALL_INCLUDE_DIR)
  set(verdict_INSTALL_INCLUDE_DIR ${CMAKE_INSTALL_INCLUDEDIR})
endif ()
if (NOT verdict_INSTALL_INCLUDE_SUBDIR)
  set(verdict_INSTALL_INCLUDE_SUBDIR verdict)
endif ()
if (NOT verdict_INSTALL_BIN_DIR)
  set(verdict_INSTALL_BIN_DIR ${CMAKE_INSTALL_BINDIR})
endif ()
if (NOT verdict_INSTALL_LIB_DIR)
  set(verdict_INSTALL_LIB_DIR ${CMAKE_INSTALL_LIBDIR})
endif ()
if (NOT verdict_INSTALL_CMAKE_DIR)
  set(verdict_INSTALL_CMAKE_DIR ${CMAKE_INSTALL_LIBDIR}/cmake/verdict)
endif ()

# Install documentation
install(FILES
  README.md
  DESTINATION ${verdict_INSTALL_DOC_DIR}/verdict/${verdict_VERSION}/ COMPONENT VerdictDevelopment)

# Install required header files
install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/verdict_config.h
  verdict.h
  DESTINATION ${verdict_INSTALL_INCLUDE_DIR} COMPONENT VerdictDevelopment)

# Install library
install(TARGETS verdict EXPORT ${VERDICT_EXPORT_GROUP}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Runtime # .exe, .dll
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Runtime # .so, .dll
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT VerdictDevelopment #[[.a, .lib]])

# Export Targets file
export(EXPORT ${VERDICT_EXPORT_GROUP}
  FILE "${CMAKE_CURRENT_BINARY_DIR}/cmake/VerdictTargets.cmake"
  NAMESPACE Verdict::)

# Install Targets file
install(EXPORT ${VERDICT_EXPORT_GROUP}
  FILE "VerdictTargets.cmake"
  NAMESPACE Verdict::
  DESTINATION ${verdict_INSTALL_CMAKE_DIR} COMPONENT VerdictDevelopment)

# Create Config file
configure_package_config_file(${CMAKE_CURRENT_SOURCE_DIR}/VerdictConfig.cmake.in
  "${CMAKE_CURRENT_BINARY_DIR}/VerdictConfig.cmake"
  INSTALL_DESTINATION ${verdict_INSTALL_CMAKE_DIR})

# Generate the version file for the Config file
write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/VerdictConfigVersion.cmake"
  VERSION "${verdict_VERSION}"
  COMPATIBILITY AnyNewerVersion)

# Install Config and ConfigVersion files
install(FILES
  "${CMAKE_CURRENT_BINARY_DIR}/VerdictConfig.cmake"
  "${CMAKE_CURRENT_BINARY_DIR}/VerdictConfigVersion.cmake"
  DESTINATION ${verdict_INSTALL_CMAKE_DIR} COMPONENT VerdictDevelopment)
