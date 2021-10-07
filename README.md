# PFSimple

## General information

The Paticle Finder Simple package is simplified version of the KF Partice package based on its mathematical apparatus. It is developed for the complete reconstruction of short-lived particles with their momentum, energy, mass, lifetime, decaylength, rapidity, etc.

## Pre-requirements

### Root

ROOT6 is needed for installation:

https://root.cern/install/build_from_source/

Follow instructions
    
### AnalysisTree (optional)

AnalysisTree is needed for usage the interface based on AnalysisTree. If it is not installed, the PFSimple will be installed without this interface.

https://github.com/HeavyIonAnalysis/AnalysisTree

Follow instructions. Use v2.1.3 tag.

## Installation

Clone PFSimple

    git clone git@git.cbm.gsi.de:pwg-c2f/analysis/pf_simple.git
    
Source ROOT

    source /path-to-root/install/bin/thisroot.sh
    
Export AnalysisTree libraries

    export AnalysisTree_DIR=/path-to-analysistree/install/lib/cmake/AnalysisTree
    
Install PFSimple
    
    mkdir build
    cd build
    cmake -DCMAKE_INSTALL_PREFIX=/path-to-install-pfsimple /path-to-source-pfsimple
    make -j install
    
## First run

Use AnalysisTreeInterface/main.cxx to configure the reconstruction: select cuts, number of events to analyze etc.