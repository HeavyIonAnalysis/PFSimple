cmake_minimum_required(VERSION 3.14)
project(PFSimple CXX)
set(PROJECT_VERSION 1.0)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

if(NOT DEFINED CMAKE_CXX_STANDARD)
    set(CMAKE_CXX_STANDARD 17)
endif()

include(ExternalProject)

macro(add_target_property _target _prop _value)
   get_target_property(_oldprop "${_target}" ${_prop})
   if(NOT _oldprop)
      set_target_properties("${_target}" PROPERTIES ${_prop} "${_value}")
   else(NOT _oldprop)
      set_target_properties("${_target}" PROPERTIES ${_prop} "${_oldprop} ${_value}")
   endif(NOT _oldprop)
endmacro(add_target_property)

set(CMAKE_CXX_STANDARD_REQUIRED ON)
message(STATUS "Using C++${CMAKE_CXX_STANDARD}")

# by default build optimized code with debug symbols
if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE RELEASE)
endif ()

option(PFSimple_BUNDLED_AT "Build AnalysisTree" OFF)
option(PFSimple_BUNDLED_AT_VERSION "Version of AnalysisTree to build" "v2.2.7")
option(PFSimple_BUILD_TESTS "Build tests for the Centrality" ON)

# in DEBUG mode make verbose Makefile
if (CMAKE_BUILD_TYPE MATCHES DEBUG)
    set(CMAKE_VERBOSE_MAKEFILE ON)
endif ()

# set(CMAKE_CXX_FLAGS_DEBUG "-O0 -ggdb -g -DDEBUG -D__DEBUG -Wall")
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -ggdb -DDEBUG -D__DEBUG -Wall")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELEASE} -ggdb")
# set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -ftree-vectorize -ffast-math -DNODEBUG")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -ftree-vectorize -ffast-math -DNODEBUG")
message(STATUS "Using CXX flags for ${CMAKE_BUILD_TYPE}: ${CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE}}")
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})

set(FIXTARGET TRUE CACHE BOOL "Compile for fix target geometry.")

find_package(ROOT REQUIRED RIO)
find_package(KFParticle QUIET)
find_package(Vc QUIET)

set(EXTERNAL_DIR ${CMAKE_BINARY_DIR}/external)
set(EXTERNAL_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/external)

if(KFParticle_FOUND)
    Message("KFParticle found")
    find_package(Vc REQUIRED)
    include_directories(${KFParticle_INCLUDE_DIR} ${Vc_INCLUDE_DIR})
else()
    Message("KFParticle is not found. It will be installed as external package")
    include(KFParticle)
    include_directories(${CMAKE_SOURCE_DIR} ${ROOT_INCLUDE_DIRS} ${PROJECT_INCLUDE_DIRECTORIES} ${Vc_INCLUDE_DIR})
endIf()
message(STATUS "PROJECT_INCLUDE_DIRECTORIES ${PROJECT_INCLUDE_DIRECTORIES}")

message(STATUS "Using ROOT: ${ROOT_VERSION} <${ROOT_CONFIG}>")
include_directories(${CMAKE_SOURCE_DIR} ${ROOT_INCLUDE_DIRS})
include(${ROOT_USE_FILE})

include(AnalysisTree)

if (AnalysisTree_FOUND OR PFSimple_BUNDLED_AT)
    add_subdirectory(at_interface)
endif ()

add_subdirectory(src)
add_subdirectory(src/interface)
include(TestsSetupPFSimple)

include(CMakePackageConfigHelpers)
write_basic_package_version_file(
        "${CMAKE_CURRENT_BINARY_DIR}/KFParticleSimple/KFParticleSimpleConfigVersion.cmake"
        VERSION ${PROJECT_VERSION}
        COMPATIBILITY AnyNewerVersion
)

export(EXPORT KFParticleSimpleTargets
        FILE "${CMAKE_CURRENT_BINARY_DIR}/KFParticleSimple/KFParticleSimpleTargets.cmake")

message(STATUS "KFParticleSimpleTargets: ${KFParticleSimpleTargets}")
        
set(ConfigPackageLocation lib/cmake/KFParticleSimple)

set(KFParticleSimple_INCLUDE_DIR "include")
set(KFParticleSimple_LIBRARY_DIR "lib")

configure_package_config_file(
        KFParticleSimpleConfig.cmake.in ${CMAKE_CURRENT_BINARY_DIR}/KFParticleSimple/KFParticleSimpleConfig.cmake
        INSTALL_DESTINATION ${ConfigPackageLocation}
        PATH_VARS KFParticleSimple_INCLUDE_DIR KFParticleSimple_LIBRARY_DIR
)

install(EXPORT KFParticleSimpleTargets
        FILE
        KFParticleSimpleTargets.cmake
        DESTINATION
        ${ConfigPackageLocation}
        )

install(FILES
        "${CMAKE_CURRENT_BINARY_DIR}/KFParticleSimple/KFParticleSimpleConfig.cmake"
        "${CMAKE_CURRENT_BINARY_DIR}/KFParticleSimple/KFParticleSimpleConfigVersion.cmake"
        DESTINATION
        ${ConfigPackageLocation}
        COMPONENT
        Devel
        )

