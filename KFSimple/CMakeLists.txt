set(SOURCES
        SimpleFinder.cpp
        SimpleFinderNew.cpp
        DaughterCuts.cpp
        MotherCuts.cpp
        Decay.cpp
        )

string(REPLACE ".cpp" ".hpp" HEADERS "${SOURCES}")
list(APPEND HEADERS "Constants.hpp")

include_directories(${CMAKE_SOURCE_DIR}/Interface ${CMAKE_SOURCE_DIR}/KFSimple ${PROJECT_INCLUDE_DIRECTORIES} )
link_directories(${PROJECT_LINK_DIRECTORIES})

add_library(KFParticleSimple SHARED ${SOURCES} G__KFParticleSimple.cxx)

ROOT_GENERATE_DICTIONARY(G__KFParticleSimple ${HEADERS} LINKDEF KFSimpleLinkDef.h OPTIONS "-DDO_TPCCATRACKER_EFF_PERFORMANCE" "-DNonhomogeneousField" "-DCBM" "-DUSE_TIMERS")
add_dependencies(KFParticleSimple KFParticleInterface)
target_link_libraries(KFParticleSimple ${ROOT_LIBRARIES} KFParticle KFParticleInterface)

add_target_property(KFParticleSimple COMPILE_FLAGS "-DDO_TPCCATRACKER_EFF_PERFORMANCE -DNonhomogeneousField -DCBM -DUSE_TIMERS")

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
