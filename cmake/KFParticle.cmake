set(KFParticle_INSTALL_DIR ${EXTERNAL_INSTALL_DIR})
set(KFParticle_INCLUDE_DIR ${KFParticle_INSTALL_DIR}/include)
set(KFParticle_LIBRARY_DIR ${KFParticle_INSTALL_DIR}/lib)

ExternalProject_Add(KFParticle_Ext
        GIT_REPOSITORY  "https://git.cbm.gsi.de/pwg-c2f/analysis/KFParticle.git"
        GIT_TAG         "cmake"
        UPDATE_DISCONNECTED ${UPDATE_DISCONNECTED}
        SOURCE_DIR      "${EXTERNAL_DIR}/KFParticle_src"
        BINARY_DIR      "${EXTERNAL_DIR}/KFParticle_build"
        INSTALL_DIR     "${KFParticle_INSTALL_DIR}"
        CMAKE_ARGS
        "-DEXTERNAL_INSTALL_DIR=${KFParticle_INSTALL_DIR}"
        "-DCMAKE_INSTALL_PREFIX=${KFParticle_INSTALL_DIR}"
        "-DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}"
        "-DCMAKE_PREFIX_PATH=${CMAKE_PREFIX_PATH}"
        "-DCMAKE_CXX_STANDARD=${CMAKE_CXX_STANDARD}"
        )
list(APPEND PROJECT_DEPENDENCIES KFParticle_Ext)
list(APPEND PROJECT_LINK_DIRECTORIES ${KFParticle_LIBRARY_DIR})
list(APPEND PROJECT_INCLUDE_DIRECTORIES ${KFParticle_INCLUDE_DIR})
