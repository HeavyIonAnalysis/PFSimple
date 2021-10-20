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

Follow instructions. Use v2.2.0 tag.

## Installation

Clone PFSimple

    git clone git@git.cbm.gsi.de:pwg-c2f/analysis/pf_simple.git
    
Source ROOT

    source /path-to-root/install/bin/thisroot.sh
    
Export AnalysisTree libraries

    export AnalysisTree_DIR=/path-to-analysistree/install/lib/cmake/AnalysisTree
    
You need to source root and export AnalysisTree each time when you are compiling project from 0 (perform cmake command) but have no need to do it when just recompiling project (perform just make).
    
Install PFSimple
    
    mkdir build
    cd build
    cmake -DCMAKE_INSTALL_PREFIX=/path-to-install-pfsimple /path-to-source-pfsimple
    make -j install
    
## First run

Use AnalysisTreeInterface/main.cxx to configure the reconstruction: select cuts, number of events to analyze etc.

Each time before running the prepared executable you should set the environment variables to let your system know where to find libraries:

    source /path-to-root/install/bin/thisroot.sh
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/path-to-pfsimple-installation/lib/
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/path-to-pfsimple-installation/external/lib/
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/path-to-analysistree-installation/lib/
 
Then run the executable:

    ./main2 filelist.txt
    
where filelist.txt must be a text file with names (including paths) to files which you want to analyze with PFSimple. Each file name should be on the next line, and the last symbol should also be a switch to next line. The example of filelist is in the source directory of PFSimple (filelist_example.txt)