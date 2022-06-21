# ==============================================================================
#
#                                 r  i  t  a
#
#            An environment for Modelling and Numerical Simulation
#
# ==============================================================================
#
#   Copyright (C) 2021 - 2022 Rachid Touzani
#
#   This file is part of rita.
#
#   rita is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   rita is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
# ==============================================================================
#
#                        cmake for the rita project
#
# ==============================================================================

cmake_minimum_required (VERSION 3.16.0)
message (STATUS "This is CMake ${CMAKE_VERSION}")
message (STATUS "")

# General information

set (RITA_PACKAGE_VENDOR
     "The rita authors <http://ofeli.org/rita/>"
    )
set (RITA_PACKAGE_DESCRIPTION
     "A console application to solve pde's"
    )

# Set up rita Version

file (STRINGS "${CMAKE_SOURCE_DIR}/VERSION" _version LIMIT_COUNT 2)
set (RITA_VERSION "${_version}")
string (REGEX REPLACE "^([0-9]+)\\..*" "\\1"
        RITA_VERSION_MAJOR "${RITA_VERSION}"
       )
string (REGEX REPLACE "^[0-9]+\\.([0-9]+).*" "\\1"
        RITA_VERSION_MINOR "${RITA_VERSION}"
       )
string (REGEX REPLACE "^[0-9]+\\.[0-9]+\\.([0-9]+).*" "\\1"
        RITA_VERSION_SUBMINOR "${RITA_VERSION}"
       )

string (SUBSTRING "${RITA_VERSION}" 6 -1 RITA_VERSION_DATE)
set (RITA_VERSION ${RITA_VERSION_MAJOR}.${RITA_VERSION_MINOR}.${RITA_VERSION_SUBMINOR})
string (REPLACE "^" "" "${RITA_VERSION_DATE}" RITA_VERSION_DATE)

project (rita)

message (STATUS "Building release ${PROJECT_NAME}-${RITA_VERSION}, ${RITA_VERSION_DATE}")

# Configuration

set (CMAKE_BUILD_TYPE "Release")
set (CMAKE_EXPORT_COMPILE_COMMANDS ON)
set (CMAKE_CXX_STANDARD 11)
enable_testing ()

option (BUILD_TESTS "Build test programs" OFF)

add_definitions (-DUSE_GMSH -Wall -fPIE -Ofast -mtune=native -D_REENTRANT -DNDEBUG -Wunknown-attributes -fstrict-aliasing -fno-strict-overflow)

# Define installation directories
#
include (GNUInstallDirs)
set (OFELI_LIBDIR ${CMAKE_INSTALL_PREFIX}/lib)
set (OFELI_INCDIR ${CMAKE_INSTALL_PREFIX}/include/ofeli)
set (GMSH_LIBDIR ${CMAKE_INSTALL_PREFIX}/lib)
set (GMSH_INCDIR ${CMAKE_INSTALL_PREFIX}/include)
set (INSTALL_BINDIR ${CMAKE_INSTALL_PREFIX}/bin)
set (INSTALL_RITADIR ${CMAKE_INSTALL_PREFIX}/share/rita)
set (INSTALL_DOCDIR ${INSTALL_RITADIR}/doc)
set (INSTALL_TUTORIALDIR ${INSTALL_RITADIR}/tutorial)

configure_file (${CMAKE_CURRENT_SOURCE_DIR}/src/config.h.in ${CMAKE_CURRENT_BINARY_DIR}/config.h)
file (APPEND ${CMAKE_CURRENT_BINARY_DIR}/config.h "\n\n#define RITA_VERSION \"${RITA_VERSION}\"\n")
file (APPEND ${CMAKE_CURRENT_BINARY_DIR}/config.h "#define RITA_DATE \"${RITA_VERSION_DATE}\"")

include_directories (${OFELI_INCDIR} ${GMSH_INCDIR})

add_subdirectory (src)
add_subdirectory (muparserx)
add_subdirectory (tutorial)

# Directories to install
install (DIRECTORY doc/ DESTINATION ${INSTALL_DOCDIR})
install (DIRECTORY tutorial/ DESTINATION ${INSTALL_TUTORIALDIR})

# Preparing Cpack
include (${CMAKE_SOURCE_DIR}/cmake/setup_cpack.cmake)

# uninstall target
if (NOT TARGET uninstall)
   configure_file (
                   "${CMAKE_CURRENT_SOURCE_DIR}/cmake/uninstall.cmake.in"
                   "${CMAKE_CURRENT_BINARY_DIR}/uninstall.cmake"
                   IMMEDIATE @ONLY
                  )

   add_custom_target (uninstall COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake/uninstall.cmake)
endif ()
