set(SOURCES
    ConverterIn.cpp
    ConverterOut.cpp
    PFSimpleTask.cpp
  )

string(REPLACE ".cpp" ".hpp" HEADERS "${SOURCES}")

include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/src ${CMAKE_SOURCE_DIR}/src/interface ${AnalysisTree_INCLUDE_DIR})

if(PROJECT_LINK_DIRECTORIES)
    include_directories(${PROJECT_INCLUDE_DIRECTORIES})
    link_directories(${PROJECT_LINK_DIRECTORIES} ${Vc_LIB_DIR})
endif()

add_library(KFMan SHARED ${SOURCES} G__KFMan.cxx)

ROOT_GENERATE_DICTIONARY(G__KFMan ${HEADERS} LINKDEF KFManLinkDef.h OPTIONS "-DDO_TPCCATRACKER_EFF_PERFORMANCE" "-DNonhomogeneousField" "-DCBM" "-DUSE_TIMERS")

target_link_libraries(KFMan ${ROOT_LIBRARIES} AnalysisTreeBase AnalysisTreeInfra KFParticleInterface KFParticleSimple KFParticle)
add_dependencies(KFMan KFParticleSimple)
add_target_property(KFMan COMPILE_FLAGS "-DDO_TPCCATRACKER_EFF_PERFORMANCE -DNonhomogeneousField -DCBM -DUSE_TIMERS")

add_executable(main_json main_json.cpp)
add_dependencies(main_json KFMan)
add_target_property(main_json COMPILE_FLAGS "-DDO_TPCCATRACKER_EFF_PERFORMANCE -DNonhomogeneousField -DCBM -DUSE_TIMERS")
target_include_directories(main_json PRIVATE ${ROOT_INCLUDE_DIRS})
target_link_libraries(main_json KFMan KFParticleSimple KFParticleInterface KFParticle Vc)

install(TARGETS KFMan EXPORT KFManTargets
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
        RUNTIME DESTINATION bin
        INCLUDES DESTINATION include
        )
install(
        FILES
        ${HEADERS}
        DESTINATION
        include
        COMPONENT
        Devel
)

set(PCM_FILE_NAME libKFMan)
install(
        FILES
        "${CMAKE_CURRENT_BINARY_DIR}/${PCM_FILE_NAME}_rdict.pcm"
        DESTINATION
        lib
        OPTIONAL
)
install (TARGETS main_json RUNTIME DESTINATION bin)
