if(PFSimple_BUNDLED_AT)
    message("-- Building bundled AnalysisTree")
    include(FetchContent)

    FetchContent_Declare(
            AnalysisTree
            GIT_REPOSITORY ${PFSimple_BUNDLED_AT_URL}
            GIT_TAG ${PFSimple_BUNDLED_AT_VERSION}
            GIT_SHALLOW ${PFSimple_BUNDLED_AT_GIT_SHALLOW}
    )
    FetchContent_MakeAvailable(AnalysisTree)
    list(APPEND PROJECT_INCLUDE_DIRECTORIES ${AnalysisTree_BINARY_DIR}/include)
    list(APPEND PROJECT_LINK_LIBRARIES AnalysisTreeBase AnalysisTreeInfra)
    message(STATUS "AT: ${AnalysisTree_BINARY_DIR} ${AnalysisTree_SOURCE_DIR}")
else()
    find_package(AnalysisTree QUIET)
    if(AnalysisTree_FOUND)
        list(APPEND CMAKE_PREFIX_PATH ${ANALYSISTREE_HOME})
        list(APPEND CMAKE_PREFIX_PATH $ENV{ANALYSISTREE_HOME})
        list(APPEND PROJECT_INCLUDE_DIRECTORIES ${AnalysisTree_INCLUDE_DIR})
        list(APPEND PROJECT_LINK_LIBRARIES AnalysisTreeBase AnalysisTreeInfra)
    else()
        message(FATAL_ERROR "Set AnalysisTree_DIR to a directory containing AnalysisTreeConfig.cmake file\n"
                            "OR\nuse -DPFSimple_BUNDLED_AT=ON option to build AnalysisTree automatically")
    endif()
endif()
