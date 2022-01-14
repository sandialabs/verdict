# Copyright 2003,2006,2019 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
# Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights in this software.

cmake_minimum_required( VERSION 3.16 )

project(verdict)
  
if(POLICY CMP0063)
  cmake_policy(SET CMP0063 NEW)
endif()

set( verdict_MAJOR_VERSION "1")
set( verdict_MINOR_VERSION "4")
set( verdict_BUILD_VERSION "0")
set( verdict_VERSION_FLAT "${verdict_MAJOR_VERSION}${verdict_MINOR_VERSION}${verdict_BUILD_VERSION}" )
set( verdict_VERSION "${verdict_MAJOR_VERSION}.${verdict_MINOR_VERSION}.${verdict_BUILD_VERSION}" )

option(VERDICT_BUILD_DOC "Build the 2007 Verdict User Manual" OFF)
option(VERDICT_MANGLE "Mangle verdict names for inclusion in a larger library?" OFF)
if ( VERDICT_MANGLE )
  set( VERDICT_MANGLE_PREFIX "verdict" CACHE STRING "The namespace enclosing verdict function names and classes." )
  mark_as_advanced( VERDICT_MANGLE_PREFIX )
endif ()
mark_as_advanced( VERDICT_MANGLE )

option( VERDICT_ENABLE_TESTING "Should tests of the VERDICT library be built?" ON )

set( verdict_SRCS
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
  verdict.h
  VerdictVector.cpp
  VerdictVector.hpp
  verdict_defines.hpp
  )

configure_file(
  ${verdict_SOURCE_DIR}/verdict_config.h.in
  ${verdict_BINARY_DIR}/verdict_config.h
  @ONLY
)

add_library( verdict ${verdict_SRCS} )
target_include_directories(verdict PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>)


# Setting the VERSION and SOVERSION of a library will include
# version information either in the library, or in the library
# name (depending on the platform). You may choose to exclude
# this information.
if ( NOT VERDICT_NO_LIBRARY_VERSION )
  set_target_properties( verdict PROPERTIES
    VERSION "${verdict_VERSION}"
    SOVERSION "${verdict_MAJOR_VERSION}.${verdict_MINOR_VERSION}"
  )
endif ()

if ( NOT VERDICT_EXPORT_GROUP )
  set(VERDICT_EXPORT_GROUP VerdictExport)
endif ()

if ( VERDICT_ENABLE_TESTING )
  enable_testing()
  if(NOT TARGET GTest::GTest)
    find_package(GTest REQUIRED)
  endif()

  function(ADD_VERDICT_UNITTESTS DIR)
    ADD_SUBDIRECTORY(${verdict_SOURCE_DIR}/${DIR} ${verdict_BINARY_DIR}/${DIR})
  
    FOREACH(TEST ${ARGN})
      ADD_TEST(NAME UT-${TEST}
               WORKING_DIRECTORY ${verdict_BINARY_DIR}/${DIR}
               COMMAND ${TEST})
    ENDFOREACH(TEST ${ARGN})
  endfunction()

  ADD_VERDICT_UNITTESTS(unittests unittests_verdict)
endif ()

if ( NOT verdict_INSTALL_DOC_DIR )
  set (verdict_INSTALL_DOC_DIR doc)
endif ()
if ( NOT verdict_INSTALL_INCLUDE_DIR)
  set (verdict_INSTALL_INCLUDE_DIR include)
endif ()
if ( NOT verdict_INSTALL_INCLUDE_SUBDIR)
  set (verdict_INSTALL_INCLUDE_SUBDIR verdict)
endif ()
if ( NOT verdict_INSTALL_BIN_DIR)
  set (verdict_INSTALL_BIN_DIR bin)
endif ()
if ( NOT verdict_INSTALL_LIB_DIR)
  set (verdict_INSTALL_LIB_DIR lib)
endif ()

if ( VERDICT_BUILD_DOC )
  add_subdirectory( docs/VerdictUserManual2007 )
endif ()

#
# Installation stuff
#
install(FILES 
  README
  DESTINATION ${verdict_INSTALL_DOC_DIR}/verdict/${verdict_VERSION}/
  COMPONENT VerdictDevelopment
)

install(
  FILES ${verdict_BINARY_DIR}/verdict_config.h
  DESTINATION ${verdict_INSTALL_INCLUDE_DIR}
  COMPONENT VerdictDevelopment
)

install(TARGETS verdict EXPORT ${VERDICT_EXPORT_GROUP}
      RUNTIME DESTINATION ${verdict_INSTALL_BIN_DIR} COMPONENT Runtime # .exe, .dll
      LIBRARY DESTINATION ${verdict_INSTALL_LIB_DIR} COMPONENT Runtime # .so, .dll
      ARCHIVE DESTINATION ${verdict_INSTALL_LIB_DIR} COMPONENT VerdictDevelopment      # .a, .lib
)

export(EXPORT ${VERDICT_EXPORT_GROUP} FILE VerdictConfig.cmake)

install(EXPORT ${VERDICT_EXPORT_GROUP} FILE VerdictConfig.cmake
      DESTINATION ${verdict_INSTALL_LIB_DIR} COMPONENT VerdictDevelopment
)
