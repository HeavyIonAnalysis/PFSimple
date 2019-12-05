set(SOURCES
    CutsContainer.cxx
    InputContainer.cxx
    OutputContainer.cxx
    )

string(REPLACE ".cxx" ".h" HEADERS "${SOURCES}")
include_directories(${CMAKE_SOURCE_DIR}/KFParticle ${CMAKE_SOURCE_DIR}/KFSimple )

add_library(KFParticleInterface SHARED ${SOURCES} G__KFParticleInterface.cxx)

if(ROOT_FOUND)
    ROOT_GENERATE_DICTIONARY(G__KFParticleInterface ${HEADERS} LINKDEF KFParticleInterfaceLinkDef.h OPTIONS "-DDO_TPCCATRACKER_EFF_PERFORMANCE" "-DNonhomogeneousField" "-DCBM" "-DUSE_TIMERS")
    add_dependencies(KFParticleInterface KFParticle)
    target_link_libraries(KFParticleInterface ${ROOT_LIBRARIES} KFParticle ${Vc_LIBRARIES})
endif(ROOT_FOUND)

add_target_property(KFParticleInterface COMPILE_FLAGS "-DDO_TPCCATRACKER_EFF_PERFORMANCE -DNonhomogeneousField -DCBM -DUSE_TIMERS")

set(PCM_FILE_NAME libKFParticleInterface)
install(FILES ${HEADERS} ${NODICT_HEADERS} DESTINATION include)
#install(FILES ${CMAKE_BINARY_DIR}/${CMAKE_SHARED_LIBRARY_PREFIX}KFParticle.rootmap DESTINATION lib)

install(TARGETS KFParticleInterface EXPORT KFParticleTargets
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
