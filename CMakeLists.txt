cmake_minimum_required (VERSION 2.8.12)
project (pomerol2triqs CXX)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
#include(CommonDefs)

# Disable build in source
if (${CMAKE_BINARY_DIR} STREQUAL ${CMAKE_SOURCE_DIR})
    message(FATAL_ERROR "In source builds are disabled. Please use a separate build directory.")
endif()
set(CMAKE_DISABLE_SOURCE_CHANGES ON)
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)

SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

find_package(pomerol)
message(STATUS "Pomerol includes : ${pomerol_INCLUDE_DIRS}")
message(STATUS "Pomerol libraries : ${pomerol_LIBRARIES}")

include_directories(${pomerol_INCLUDE_DIRS})
include_directories(include)

add_definitions(-std=c++11)

set (sources
    src/ReadWrite.cpp
    src/Params.cpp
    src/OperatorPresetsExtra.cpp
)

add_executable(pomerol2dcore src/pomerol2dcore.cpp ${sources})
target_link_libraries(pomerol2dcore ${pomerol_LIBRARIES})

install(TARGETS pomerol2dcore DESTINATION bin)
