cmake_minimum_required(VERSION 3.16.0)

project(CSpectra VERSION 1.0.1
                 DESCRIPTION "A header-only C++ library for large scale eigenvalue problems"
                 LANGUAGES C CXX)

# Make CMake look into the ./cmake/ folder for configuration files
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

# Look for supporting libraries

if (NOT TARGET Eigen3::Eigen)
  find_package(Eigen3 NO_MODULE REQUIRED)
  message(STATUS "Found Eigen3 Version: ${Eigen3_VERSION} Path: ${Eigen3_DIR}")
endif()

if (NOT TARGET Spectra::Spectra)
  find_package(Spectra NO_MODULE REQUIRED)
  message(STATUS "Found Spectra Version: ${Spectra_VERSION} Path: ${Spectra_DIR}")
endif()

# Setup library

add_library(cspectra)

target_sources(cspectra
  PRIVATE
    src/spectra_di_ns.cpp
    src/spectra_di_sg.cpp
    src/spectra_di_ss.cpp
  PUBLIC
    src/common.h
    src/spectra_di_ns.h
    src/spectra_di_sg.h
    src/spectra_di_ss.h)

target_link_libraries(cspectra PUBLIC Eigen3::Eigen)
target_link_libraries(cspectra PUBLIC Spectra::Spectra)
