set(SOURCES
    InputContainer.cpp
    OutputContainer.cpp
    )

string(REPLACE ".cpp" ".hpp" HEADERS "${SOURCES}")
include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/src)

include_directories(${PROJECT_INCLUDE_DIRECTORIES})
link_directories(${PROJECT_LINK_DIRECTORIES})

add_library(KFParticleInterface SHARED ${SOURCES} G__KFParticleInterface.cxx)

if(PROJECT_DEPENDENCIES)
    add_dependencies(KFParticleInterface ${PROJECT_DEPENDENCIES})
endif()

ROOT_GENERATE_DICTIONARY(G__KFParticleInterface ${HEADERS} LINKDEF KFParticleInterfaceLinkDef.h OPTIONS "-DDO_TPCCATRACKER_EFF_PERFORMANCE" "-DNonhomogeneousField" "-DCBM" "-DUSE_TIMERS")
target_link_libraries(KFParticleInterface ${ROOT_LIBRARIES} KFParticle ${Vc_LIBRARIES})

add_target_property(KFParticleInterface COMPILE_FLAGS "-DDO_TPCCATRACKER_EFF_PERFORMANCE -DNonhomogeneousField -DCBM -DUSE_TIMERS")

if(PFSimple_BUILD_TESTS)
    enable_testing()
    include(GoogleTest)

    set(TEST_SOURCES
                    InputContainer.test.cpp
                    OutputContainer.test.cpp
                    )

    add_executable(PFSimpleInterface_UnitTests ${TEST_SOURCES})
    add_dependencies(PFSimpleInterface_UnitTests KFParticleInterface)
    add_target_property(PFSimpleInterface_UnitTests COMPILE_FLAGS "-DDO_TPCCATRACKER_EFF_PERFORMANCE -DNonhomogeneousField -DCBM -DUSE_TIMERS")
    target_link_libraries(PFSimpleInterface_UnitTests  PRIVATE KFParticleInterface KFParticle Vc gtest_main)
    gtest_discover_tests(PFSimpleInterface_UnitTests TEST_PREFIX PFSimpleInterface_)
endif()

set(PCM_FILE_NAME libKFParticleInterface)
install(FILES ${HEADERS} DESTINATION include)

install(TARGETS KFParticleInterface EXPORT KFParticleSimpleTargets
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
        RUNTIME DESTINATION bin
        INCLUDES DESTINATION include
        )

install(FILES
        "${CMAKE_CURRENT_BINARY_DIR}/${PCM_FILE_NAME}_rdict.pcm"
        DESTINATION
        lib
        OPTIONAL
        )
