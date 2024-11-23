# PFSimple

## General information

The Paticle Finder Simple package is simplified version of the KF Partice package based on its mathematical apparatus. It is developed for the complete reconstruction of short-lived particles with their momentum, energy, mass, lifetime, decaylength, rapidity, etc.

## Requirements

### Root

https://root.cern/install/

Version 6.20 or newer is recommended. Root with c++17 standard is strongly recommended.

### AnalysisTree

Optional. Needed for building AnalysisTree interface.

https://github.com/HeavyIonAnalysis/AnalysisTree

Both pre-installed AnalysisTree and automatically installed together with PFSimple framework can be used (see chapter Installation)

AnalysisTree with c++17 standard is strongly recommended.

## Installation

Clone PFSimple
    
Source ROOT

    source /path-to-root/install/bin/thisroot.sh
    
To use pre-installed AnalysisTree set the environment variable:
    export AnalysisTree_DIR=/path-to-analysistree/install/lib/cmake/AnalysisTree
To build AnalysisTree automatically together with PFSimple use following cmake keys:
    -DPFSimple_BUNDLED_AT=ON
    -DPFSimple_BUNDLED_AT_VERSION=v2.2.7
    
You need to source root and export AnalysisTree each time when you are compiling project from 0 (perform cmake command) but have no need to do it when just recompiling project (perform just make).
    
Install PFSimple
    
    mkdir build
    cd build
    cmake -DCMAKE_INSTALL_PREFIX=/path-to-install-pfsimple /path-to-source-pfsimple
    make -j install
    
## Configuration of decay settings

The reconstruction will be executed with at_interface/main2.cpp for
2-body decays and at_interface/main3.cpp for 3-body decays.

The reconstruction can be configured with the parameter-file
at_interface/parfile2.txt for 2-body decays and at_interface/parfile3.txt for 3-body decays.

The user can select:
1) Decay: mother and daughters

   Option for daughters: alternative pdgs can be used in addition
   for the reconstruction
   - for background: set to "1" / "-1" for pos / neg daughters
   - for tracks without TOF id: set to "2" / "-2" for pos / neg
   daughters (in pid-mode 2, 3, 4)
   - In pid-mode 0 optional daughers "1" / "-1" for pos / neg pdgs must be added.

   If no alternative pdgs should be used, set to "0".

   If more optional pdgs than 2 should be added, change in main2.cpp /
   main3.cpp and recompile.

2) Geometrical cuts to be applied and cut-values

   If a cut should not be used, set to "-1".

   For more information on specific cuts see " Constants.hpp".

3) Pid-method:
   
   - 0: no pid (only charge information is used ! optional daughers must
   be set to "1" / "-1" for pos / neg pdgs !)
   
   - 1: mc-pid
 
   - 2: reconstructed TOF pid - default from Pid-framework
   
   - 3: reconstructed TOF pid - pdg with max. pdg-purity is selected, if pdg-purity >min. purity
   
   - 4: reconstructed TOF pid - pdg is selected, if pdg-purity > min. purity (pdg-specfic purities are possible)

   for 2-4: Pid-framework needs to be applied first to
   include TOF-pid & probabilities into the analysistree input-file

4) Minium purities for pid-mode 3 & 4:
   - for pid-mode 3: minimum purity is set for all pdgs (pdg-spefic purities
   need to be set to "-1")
   - for pid-mode 4: minimum purity can be set for all pdgs or specifically
   for every pdg (general purity will be overwritten, if pdg-specific
   purity is not set to "-1")

5) Input/Output information:
   - treename in input aTree, default names: "rTree" for standard aTree,
   "aTree" after running Pid-framework
   - branchname of reconstructed tracks in input analysistree, default names: "VtxTracks"
     for standard analyistree, "RecParticles" after running Pid-framework
   - number of events to be processed (-1 = all events)
   - Output format: default output is analysistree format
      options:
       - plain tree: contains information about the candidates (is produced in addition to analysistree)
	   - root tree: contains information about candidates, matched mc-particles and events (no analyistree is written)

If more than one decay should be reconstructed in the same run, additional
parfiles at_interface/parfile2_add.txt etc. can be added where 1.) & 2.) can be
selected. 3.) - 5.) can not be changed in the same run.

For 3-body decays, the list of available cuts is extended, amoung others, to cuts
on the secondary mothers (SM) of combinations of 2 daughters. E.g. for
H3L->p + pi + d, cuts on the SM of p-pi can be applied. The cuts
that test the mother against the PV, e.g. chi2topo, are inverted for
the SM to exclude SMs that come from the PV.

Several example parameter files can be found in the folder at_interface.

## First run

Each time before running the prepared executable you should set the environment variables to let your system know where to find libraries:

    source /path-to-root/install/bin/thisroot.sh
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/path-to-pfsimple-installation/lib/
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/path-to-pfsimple-installation/external/lib/
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/path-to-analysistree-installation/lib/
 
Then run the executable:

    ./main2 filelist.txt dataset parfile2.txt (parfile2_add.txt)

or

    ./main3 filelist.txt dataset parfile3.txt (parfile3_add.txt)
	
where filelist.txt must be a text file with names (including paths) to
files which you want to analyze with PFSimple. Each file name should
be on the next line, and the last symbol should also be a switch to
next line. The example of filelist is in the source directory of
PFSimple (filelist_example.txt)

and where [dataset] specifies the names of the outputfiles:
- analyistree: [dataset].PFSimpleOutput.root
- plain tree: [dataset].PFSimplePlainTree.root
- root tree: [dataset].PFSimpleTree.root

