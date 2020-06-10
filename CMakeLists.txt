cmake_minimum_required(VERSION 3.0)
project(KFParticleSimple CXX)
set(CMAKE_CXX_STANDARD 11)
set(PROJECT_VERSION 1.0)

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
    set(CMAKE_BUILD_TYPE RELWITHDEBINFO)
endif ()

# in DEBUG mode make verbose Makefile
if (CMAKE_BUILD_TYPE MATCHES DEBUG)
    set(CMAKE_VERBOSE_MAKEFILE ON)
endif ()

# set(CMAKE_CXX_FLAGS_DEBUG "-O0 -ggdb -g -DDEBUG -D__DEBUG -Wall")
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -ggdb -DDEBUG -D__DEBUG -Wall")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELEASE} -ggdb")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -ftree-vectorize -ffast-math -DNODEBUG")
message(STATUS "Using CXX flags for ${CMAKE_BUILD_TYPE}: ${CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE}}")
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})

set(FIXTARGET TRUE CACHE BOOL "Compile for fix target geometry.")

find_package(ROOT REQUIRED RIO)
find_package(KFParticle QUIET)
find_package(Vc QUIET)

set(EXTERNAL_DIR ${CMAKE_BINARY_DIR}/external)
set(EXTERNAL_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/external)

if(Vc_FOUND)
    Message("Vc found")
Else()
    Message("Vc is not found. It will be installed as external package")
    include(cmake_modules/Vc.cmake)
EndIf()

if(KFParticle_FOUND)
    Message("KFParticle found")
Else()
    Message("KFParticle is not found. It will be installed as external package")
    include(cmake_modules/KFParticle.cmake)
EndIf()

if (ROOT_FOUND)
    message(STATUS "Using ROOT: ${ROOT_VERSION} <${ROOT_CONFIG}>")
    include_directories(${CMAKE_SOURCE_DIR} ${ROOT_INCLUDE_DIR} ${ROOT_INCLUDE_DIRS})
    include(${ROOT_USE_FILE})
endif (ROOT_FOUND)

include_directories(${CMAKE_SOURCE_DIR} ${KFParticle_INCLUDE_DIR} ${Vc_INCLUDE_DIR})

add_subdirectory(KFSimple)
add_subdirectory(Interface)

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
