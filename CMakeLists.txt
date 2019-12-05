cmake_minimum_required(VERSION 2.8.8 FATAL_ERROR)

project(KFParticleSimpleLibrary)

set(PROJECT_VERSION 1.0)

# by default build optimized code with debug symbols
if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE RELWITHDEBINFO)
endif ()

# in DEBUG mode make verbose Makefile
if (CMAKE_BUILD_TYPE MATCHES DEBUG)
    set(CMAKE_VERBOSE_MAKEFILE ON)
endif ()

# set(CMAKE_CXX_FLAGS_DEBUG "-O0 -ggdb -g -DDEBUG -D__DEBUG -Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -ggdb -DDEBUG -D__DEBUG -Wall")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELEASE} -ggdb")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -ftree-vectorize -ffast-math -DNODEBUG")
message(STATUS "Using CXX flags for ${CMAKE_BUILD_TYPE}: ${CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE}}")
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})

macro(add_target_property _target _prop _value)
  get_target_property(_oldprop "${_target}" ${_prop})
  if(NOT _oldprop)
    set_target_properties("${_target}" PROPERTIES ${_prop} "${_value}")
  else(NOT _oldprop)
    set_target_properties("${_target}" PROPERTIES ${_prop} "${_oldprop} ${_value}")
  endif(NOT _oldprop)
endmacro(add_target_property)

if (CMAKE_VERSION VERSION_LESS "3.1")
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  set (CMAKE_CXX_FLAGS "--std=gnu++11 ${CMAKE_CXX_FLAGS}")
endif ()
else ()
  set (CMAKE_CXX_STANDARD 11)
endif ()

list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS} ${ROOTSYS})

# set(FIXTARGET FALSE CACHE BOOL "Compile for collider geometry.")
set(FIXTARGET TRUE CACHE BOOL "Compile for fix target geometry.")

find_package(ROOT REQUIRED COMPONENTS Core Hist Matrix Physics EG Gpad Graf Graf3d RIO MathCore)
If(DEFINED Vc_INCLUDE_DIR AND Vc_LIBRARIES)
  Message("Vc found")
Else()
  find_package(Vc REQUIRED)
EndIf()

if (ROOT_FOUND)
  message(STATUS "Using ROOT: ${ROOT_VERSION} <${ROOT_CONFIG}>")
  include_directories(${CMAKE_SOURCE_DIR} ${ROOT_INCLUDE_DIR} ${ROOT_INCLUDE_DIRS})
  include(${ROOT_USE_FILE})
endif (ROOT_FOUND)

find_package(KFParticle REQUIRED)

include(${ROOT_USE_FILE})
set(CMAKE_INCLUDE_DIRECTORIES_BEFORE ON)
Include_Directories( ${CMAKE_SOURCE_DIR} ${Vc_INCLUDE_DIR} ${INCLUDE_DIRECTORIES} KFSimple Interface)
set(CMAKE_INCLUDE_DIRECTORIES_BEFORE OFF)

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

# install(FILES ${HEADERS} ${NODICT_HEADERS} DESTINATION include)
# install(FILES ${CMAKE_BINARY_DIR}/${CMAKE_SHARED_LIBRARY_PREFIX}KFParticle_rdict.pcm DESTINATION lib OPTIONAL)
# install(FILES ${CMAKE_BINARY_DIR}/${CMAKE_SHARED_LIBRARY_PREFIX}KFParticle.rootmap DESTINATION lib)
# install(TARGETS KFParticle DESTINATION lib)
