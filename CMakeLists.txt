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

target_link_libraries(cspectra PUBLIC Eigen3::Eigen)
target_link_libraries(cspectra PUBLIC Spectra::Spectra)

target_include_directories(cspectra
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>/src
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

target_sources(cspectra
  PRIVATE
    src/spectra_di_ns.cpp
    src/spectra_di_sg.cpp
    src/spectra_di_ss.cpp)

# Install the library (relative to the CMAKE_INSTALL_PREFIX)

include(GNUInstallDirs)

set(CSPECTRA_INSTALL_INCLUDEDIR "${CMAKE_INSTALL_INCLUDEDIR}/cspectra")
set(CSPECTRA_INSTALL_CMAKEDIR "${CMAKE_INSTALL_LIBDIR}/cmake/cspectra")

install(
  TARGETS cspectra
  EXPORT CSpectra-targets
  INCLUDES DESTINATION include
)

install(
  DIRECTORY ${PROJECT_SOURCE_DIR}/src/
  DESTINATION ${CSPECTRA_INSTALL_INCLUDEDIR}
  FILES_MATCHING PATTERN "*.h")

install(
  EXPORT CSpectra-targets
  FILE cspectra-targets.cmake
  NAMESPACE CSpectra::
  DESTINATION ${CSPECTRA_INSTALL_CMAKEDIR}
)

# Configure package

include(CMakePackageConfigHelpers)

configure_package_config_file(
  ${PROJECT_SOURCE_DIR}/cmake/cspectra-config.cmake.in
  ${CMAKE_BINARY_DIR}/cmake/cspectra-config.cmake
  INSTALL_DESTINATION ${CSPECTRA_INSTALL_CMAKEDIR})

write_basic_package_version_file(
  ${CMAKE_BINARY_DIR}/cmake/cspectra-config-version.cmake
  VERSION ${Spectra_VERSION}
  COMPATIBILITY AnyNewerVersion)

install(
  FILES
    ${CMAKE_BINARY_DIR}/cmake/cspectra-config.cmake
    ${CMAKE_BINARY_DIR}/cmake/cspectra-config-version.cmake
  DESTINATION ${CSPECTRA_INSTALL_CMAKEDIR})
