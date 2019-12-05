set(SOURCES
        SimpleFinder.cxx
        )

string(REPLACE ".cxx" ".h" HEADERS "${SOURCES}")
list(APPEND HEADERS "Constants.h")

include_directories(${CMAKE_SOURCE_DIR}/KFParticle
                    ${CMAKE_SOURCE_DIR}/Interface)

add_library(KFParticleSimple SHARED ${SOURCES} G__KFParticleSimple.cxx)

if(ROOT_FOUND)
    ROOT_GENERATE_DICTIONARY(G__KFParticleSimple ${HEADERS} LINKDEF KFParticleSimpleLinkDef.h OPTIONS "-DDO_TPCCATRACKER_EFF_PERFORMANCE" "-DNonhomogeneousField" "-DCBM" "-DUSE_TIMERS")
    add_dependencies(KFParticleSimple KFParticle)
    target_link_libraries(KFParticleSimple ${ROOT_LIBRARIES} KFParticle ${Vc_LIBRARIES})
endif(ROOT_FOUND)

add_target_property(KFParticleSimple COMPILE_FLAGS "-DDO_TPCCATRACKER_EFF_PERFORMANCE -DNonhomogeneousField -DCBM -DUSE_TIMERS")

set(PCM_FILE_NAME libKFParticleSimple)
install(FILES ${HEADERS} ${NODICT_HEADERS} DESTINATION include)
#install(FILES ${CMAKE_BINARY_DIR}/${CMAKE_SHARED_LIBRARY_PREFIX}KFParticle.rootmap DESTINATION lib)

install(TARGETS KFParticleSimple EXPORT KFParticleTargets
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
