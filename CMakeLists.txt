# syntaxis: grammars of geometry and topology — graphs and complexes in space, their grammars learned, sampled and parsed.
#
# Every module is a header-only target syntaxis::<module>. syntaxis uses the targets of calx, morphe and physis: inside
# the infrastructure, the siblings under foundations/ (each added here when no larger project has added it). Built on its
# own, syntaxis builds its checks (SYNTAXIS_CHECKS); added to a larger project, only when asked.

cmake_minimum_required(VERSION 3.20)
project(syntaxis VERSION 0.1.0 LANGUAGES CXX)

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  set(SYNTAXIS_TOP_LEVEL ON)
else()
  set(SYNTAXIS_TOP_LEVEL OFF)
endif()
option(SYNTAXIS_CHECKS "Build syntaxis's checks" ${SYNTAXIS_TOP_LEVEL})

if(NOT DEFINED CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD 20)
  set(CMAKE_CXX_STANDARD_REQUIRED ON)
  set(CMAKE_CXX_EXTENSIONS OFF)
endif()
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release)
endif()

# syntaxis_sibling(<directory> <a target it defines>): the sibling library, added once.
function(syntaxis_sibling directory target)
  if(NOT TARGET ${target})
    if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/../${directory}/CMakeLists.txt)
      add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../${directory} ${CMAKE_BINARY_DIR}/${directory})
    else()
      message(FATAL_ERROR "syntaxis needs ${directory}: build it inside the infrastructure (foundations/${directory}).")
    endif()
  endif()
endfunction()

# The siblings in their order: calx and morphe, then physis.
syntaxis_sibling(calx calx::linalg)
syntaxis_sibling(morphe morphe::topology)
syntaxis_sibling(physis physis::units)
if(SYNTAXIS_CHECKS)
  enable_testing()
endif()

# syntaxis_module(<name> [LIBRARIES ...]): the target syntaxis::<name> and, with the checks, check_syntaxis_<name>, which knows where
# its reference values are (SYNTAXIS_FIXTURES = <name>/tests/fixtures).
function(syntaxis_module name)
  cmake_parse_arguments(ARG "" "" "LIBRARIES" ${ARGN})
  add_library(syntaxis_${name} INTERFACE)
  add_library(syntaxis::${name} ALIAS syntaxis_${name})
  target_include_directories(syntaxis_${name} INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/${name}/include)
  target_compile_features(syntaxis_${name} INTERFACE cxx_std_20)
  target_link_libraries(syntaxis_${name} INTERFACE ${ARG_LIBRARIES})
  if(SYNTAXIS_CHECKS)
    calx_add_check(check_syntaxis_${name} SOURCES ${name}/tests/check_${name}.cpp LIBRARIES syntaxis::${name} calx::random)
    target_compile_definitions(check_syntaxis_${name} PRIVATE SYNTAXIS_FIXTURES="${CMAKE_CURRENT_SOURCE_DIR}/${name}/tests/fixtures")
  endif()
endfunction()

# The modules in their order (decision 0003).
syntaxis_module(shape LIBRARIES physis::computational_geometry)
syntaxis_module(graphs LIBRARIES syntaxis::shape)
syntaxis_module(induction LIBRARIES syntaxis::graphs)
syntaxis_module(derivation LIBRARIES syntaxis::induction)
syntaxis_module(parsing LIBRARIES syntaxis::induction syntaxis::derivation)
syntaxis_module(stochastic LIBRARIES syntaxis::induction syntaxis::derivation calx::random)
syntaxis_module(simplicial LIBRARIES syntaxis::graphs physis::computational_geometry calx::random)
