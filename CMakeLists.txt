# Copyright (c) 2003,2006 Sandia National Laboratories <cubit@sandia.gov>
project(verdict)

cmake_minimum_required( VERSION 2.8 )

set(CMAKE_MACOSX_RPATH 1)

set( verdict_MAJOR_VERSION "1")
set( verdict_MINOR_VERSION "2")
set( verdict_BUILD_VERSION "0")
set( verdict_VERSION_FLAT "${verdict_MAJOR_VERSION}${verdict_MINOR_VERSION}${verdict_BUILD_VERSION}" )
set( verdict_VERSION "${verdict_MAJOR_VERSION}.${verdict_MINOR_VERSION}.${verdict_BUILD_VERSION}" )

option(BUILD_SHARED_LIBS "Build Verdict with shared libraries." OFF)
option(VERDICT_BUILD_DOC "Build the 2007 Verdict User Manual" OFF)

option( VERDICT_ENABLE_TESTING "Should tests of the VERDICT library be built?" ON )

include_directories(
  ${verdict_BINARY_DIR}
  ${verdict_SOURCE_DIR}
)

set( verdict_SRCS
  V_EdgeMetric.cpp
  V_GaussIntegration.cpp
  V_GaussIntegration.hpp
  V_HexMetric.cpp
  V_HexMetric.hpp
  V_KnifeMetric.cpp
  V_PyramidMetric.cpp
  V_QuadMetric.cpp
  V_TetMetric.cpp
  V_TetMetric.hpp
  V_TriMetric.cpp
  v_vector.h
  V_WedgeMetric.cpp
  VerdictVector.cpp
  VerdictVector.hpp
  verdict_defines.hpp
  )

configure_file(
  ${verdict_SOURCE_DIR}/verdict.h.in
  ${verdict_BINARY_DIR}/verdict.h
  @ONLY
)

set( verdict_LIBRARY "verdict" )

# Setting the VERSION and SOVERSION of a library will include
# version information either in the library, or in the library
# name (depending on the platform). You may choose to exclude
# this information.
if ( NOT VERDICT_NO_LIBRARY_VERSION )
  set( VERDICT_LIBRARY_PROPERTIES
    ${VERDICT_LIBRARY_PROPERTIES}
    VERSION "${verdict_VERSION}"
    SOVERSION "${verdict_MAJOR_VERSION}.${verdict_MINOR_VERSION}"
  )
endif ( NOT VERDICT_NO_LIBRARY_VERSION )

add_library( ${verdict_LIBRARY} ${verdict_SRCS} )
set_target_properties( ${verdict_LIBRARY} PROPERTIES DEFINE_SYMBOL verdict_EXPORTS)
set_target_properties( ${verdict_LIBRARY} PROPERTIES INTERFACE_INCLUDE_DIRECTORIES
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR};${CMAKE_CURRENT_BINARY_DIR}>")

# Apply user-defined properties to the library targets.
if ( VERDICT_LIBRARY_PROPERTIES )
  set_target_properties( ${verdict_LIBRARY}
    PROPERTIES ${VERDICT_LIBRARY_PROPERTIES}
  )
endif ( VERDICT_LIBRARY_PROPERTIES )

if ( VERDICT_ENABLE_TESTING )
  enable_testing()

  function(ADD_VERDICT_UNITTESTS DIR)
    ADD_SUBDIRECTORY(${verdict_SOURCE_DIR}/${DIR} ${verdict_BINARY_DIR}/${DIR})
  
    FOREACH(TEST ${ARGN})
      ADD_TEST(NAME UT-${TEST}
               WORKING_DIRECTORY ${verdict_BINARY_DIR}/${DIR}
               COMMAND ${TEST})
      set_tests_properties(UT-${TEST} PROPERTIES ENVIRONMENT "${CUBIT_TEST_ENV}")
    ENDFOREACH(TEST ${ARGN})
  endfunction()

  ADD_VERDICT_UNITTESTS(unittests unittests_verdict)
endif ( VERDICT_ENABLE_TESTING )

if ( NOT verdict_INSTALL_DOC_DIR )
  set (verdict_INSTALL_DOC_DIR doc)
endif( NOT verdict_INSTALL_DOC_DIR )
if ( NOT verdict_INSTALL_INCLUDE_DIR)
  set (verdict_INSTALL_INCLUDE_DIR include)
endif ( NOT verdict_INSTALL_INCLUDE_DIR)
if ( NOT verdict_INSTALL_INCLUDE_SUBDIR)
  set (verdict_INSTALL_INCLUDE_SUBDIR verdict)
endif ( NOT verdict_INSTALL_INCLUDE_SUBDIR)
if ( NOT verdict_INSTALL_BIN_DIR)
  set (verdict_INSTALL_BIN_DIR bin)
endif ( NOT verdict_INSTALL_BIN_DIR)
if ( NOT verdict_INSTALL_LIB_DIR)
  set (verdict_INSTALL_LIB_DIR lib)
endif ( NOT verdict_INSTALL_LIB_DIR)

if ( VERDICT_BUILD_DOC )
  add_subdirectory( docs/VerdictUserManual2007 )
endif ( VERDICT_BUILD_DOC )

#
# Installation stuff
#
if(0)
install(FILES 
  README
  Verdict.htm
  Verdict.doc
  verdict_test.cpp
  DESTINATION ${verdict_INSTALL_DOC_DIR}/verdict/${verdict_VERSION}/
  COMPONENT Development
)

install(
  FILES ${verdict_BINARY_DIR}/verdict.h
  DESTINATION ${verdict_INSTALL_INCLUDE_DIR}
  COMPONENT Development
)
endif()

if(BUILD_SHARED_LIBS)
install(TARGETS ${verdict_LIBRARY} ${VERDICT_EXPORT_GROUP}
      RUNTIME DESTINATION ${verdict_INSTALL_BIN_DIR} COMPONENT Runtime # .exe, .dll
      LIBRARY DESTINATION ${verdict_INSTALL_LIB_DIR} COMPONENT Runtime # .so, mod.dll
      #ARCHIVE DESTINATION ${verdict_INSTALL_LIB_DIR} COMPONENT Development      # .a, .lib
)
endif(BUILD_SHARED_LIBS)

#
# Packing stuff
#
if ("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_CURRENT_SOURCE_DIR}")
  if (EXISTS "${CMAKE_ROOT}/Modules/InstallRequiredSystemLibraries.cmake")
    set (CMAKE_INSTALL_MFC_LIBRARIES 1)
#      include (InstallRequiredSystemLibraries)
  endif (EXISTS "${CMAKE_ROOT}/Modules/InstallRequiredSystemLibraries.cmake")
  set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "VERDICT - a geometric quality functions library")
  set (CPACK_PACKAGE_VENDOR "Sandia National Laboratories")
  set (CPACK_PACKAGE_INSTALL_DIRECTORY "Verdict ${verdict_MAJOR_VERSION}.${verdict_MINOR_VERSION}")
  set (CPACK_SOURCE_PACKAGE_FILE_NAME "verdict-${verdict_VERSION}")
  set (CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/README")
  set (CPACK_PACKAGE_VERSION_MAJOR "${verdict_MAJOR_VERSION}")
  set (CPACK_PACKAGE_VERSION_MINOR "${verdict_MINOR_VERSION}")
  set (CPACK_PACKAGE_VERSION_PATCH "${verdict_BUILD_VERSION}")
  set (CPACK_SYSTEM_NAME ${CMAKE_SYSTEM_NAME})
  if (${CMAKE_SYSTEM_NAME} MATCHES Windows)
    if (CMAKE_CL_64)
      set (CPACK_SYSTEM_NAME win64)
    else (CMAKE_CL_64)
      set (CPACK_SYSTEM_NAME win32)
    endif (CMAKE_CL_64)
  endif (${CMAKE_SYSTEM_NAME} MATCHES Windows)

  set (CPACK_PACKAGE_FILE_NAME "${CPACK_SOURCE_PACKAGE_FILE_NAME}-${CPACK_SYSTEM_NAME}")
  if (WIN32 AND NOT UNIX)
    # There is a bug in NSI that does not handle full unix paths properly.
    # Make sure there is at least one set of four (4) backslashes.
    set (CPACK_PACKAGE_ICON "${verdict_SOURCE_DIR}/Utilities/Release\\\\VerdictIcon.bmp")
    #set (CPACK_PACKAGE_EXECUTABLES "CMakeSetup" "CMake")
    #set (CPACK_NSIS_INSTALLED_ICON_NAME "bin\\\\CMakeSetup.exe")
    set (CPACK_NSIS_DISPLAY_NAME "${CPACK_PACKAGE_INSTALL_DIRECTORY} a geometric quality functions library")
    set (CPACK_NSIS_HELP_LINK "http:\\\\\\\\cubit.sandia.gov/verdict.html")
    set (CPACK_NSIS_URL_INFO_ABOUT "http:\\\\\\\\www.sandia.gov")
    set (CPACK_NSIS_CONTACT "cubit@sandia.gov")
    set (CPACK_NSIS_MODIFY_PATH ON)
  else (WIN32 AND NOT UNIX)
    #set (CPACK_STRIP_FILES "bin/ccmake;bin/cmake;bin/cpack;bin/ctest")
    #set (CPACK_PACKAGE_EXECUTABLES "")
  endif (WIN32 AND NOT UNIX)
  include(CPack)
endif ("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_CURRENT_SOURCE_DIR}")
