set(SOURCES
        SimpleFinder.cpp
        Daughter.cpp
        Mother.cpp
        Decay.cpp
        NonLinearCutBase.cpp
        )

string(REPLACE ".cpp" ".hpp" HEADERS "${SOURCES}")
list(APPEND HEADERS "Constants.hpp")

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/interface ${CMAKE_CURRENT_SOURCE_DIR} ${PROJECT_INCLUDE_DIRECTORIES} )
link_directories(${PROJECT_LINK_DIRECTORIES})

add_library(KFParticleSimple SHARED ${SOURCES} G__KFParticleSimple.cxx)

ROOT_GENERATE_DICTIONARY(G__KFParticleSimple ${HEADERS} LINKDEF KFSimpleLinkDef.h OPTIONS "-DDO_TPCCATRACKER_EFF_PERFORMANCE" "-DNonhomogeneousField" "-DCBM" "-DUSE_TIMERS")
add_dependencies(KFParticleSimple KFParticleInterface)
target_link_libraries(KFParticleSimple ${ROOT_LIBRARIES} KFParticle KFParticleInterface)

add_target_property(KFParticleSimple COMPILE_FLAGS "-DDO_TPCCATRACKER_EFF_PERFORMANCE -DNonhomogeneousField -DCBM -DUSE_TIMERS")

if(PFSimple_BUILD_TESTS)
    enable_testing()
    include(GoogleTest)

    set(TEST_SOURCES
                    Daughter.test.cpp
                    Mother.test.cpp
                    Decay.test.cpp
                    )

    add_executable(PFSimple_UnitTests ${TEST_SOURCES})
    target_link_libraries(PFSimple_UnitTests  PRIVATE KFParticleSimple gtest_main)
    gtest_discover_tests(PFSimple_UnitTests TEST_PREFIX PFSimple_)
endif()

set(PCM_FILE_NAME libKFParticleSimple)
install(FILES ${HEADERS} ${NODICT_HEADERS} DESTINATION include)
#install(FILES ${CMAKE_BINARY_DIR}/${CMAKE_SHARED_LIBRARY_PREFIX}KFParticle.rootmap DESTINATION lib)

install(TARGETS KFParticleSimple EXPORT KFParticleSimpleTargets
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

#set(PCM_FILE_NAME libKFParticle)
#install(FILES ${HEADERS} ${NODICT_HEADERS} DESTINATION include)
##install(FILES ${CMAKE_BINARY_DIR}/${CMAKE_SHARED_LIBRARY_PREFIX}KFParticle.rootmap DESTINATION lib)
#
#install(TARGETS KFParticle EXPORT KFParticleTargets
#        LIBRARY DESTINATION lib
#        ARCHIVE DESTINATION lib
#        RUNTIME DESTINATION bin
#        INCLUDES DESTINATION include
#        )
#
install(FILES
        "${CMAKE_CURRENT_BINARY_DIR}/${PCM_FILE_NAME}_rdict.pcm"
        DESTINATION
        lib
        OPTIONAL
)
