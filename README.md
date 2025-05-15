# PFSimple

## General information

The Paticle Finder Simple package is simplified version of the KF Partice package based on its mathematical apparatus. It is developed for the complete reconstruction of short-lived particles with their momentum, energy, mass, lifetime, decaylength, rapidity, etc. 2- and 3-body decays can get reconstructed. It is also possible to reconstruct multiple decays in one run as well as cascade decays.

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
```
source /path-to-root/install/bin/thisroot.sh
```
To use pre-installed AnalysisTree set the environment variable:
```
export AnalysisTree_DIR=/path-to-analysistree/install/lib/cmake/AnalysisTree
```
To build AnalysisTree automatically together with PFSimple use following cmake keys:
```
-DPFSimple_BUNDLED_AT=ON
-DPFSimple_BUNDLED_AT_VERSION=v2.3.4
```
    
You need to source root and export AnalysisTree each time when you are compiling project from 0 (perform cmake command) but have no need to do it when just recompiling project (perform just make).
    
Install PFSimple
```
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=/path-to-install-pfsimple /path-to-source-pfsimple
make -j install
```
## Configuration of decay settings
The program is configured with a JSON config file which contains various settings related to input and output as well as the decays to reconstruct. In the following, the different blocks are explained.

### io
| Key                            | Description                                                                                                                                               |
| ------------------------------ | --------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `input_treename`               | Name of the input AnalysisTree (for standard analysistree `rTree`, after running PID framework `pTree`)                                                   |
| `rectracks_branchname`         | branchname of reconstructed tracks in input analysistree, default names: `VtxTracks` for standard analyistree, `RecParticles` after running Pid-framework |
| `n_events`                     | Number of events to be processed, set to `-1` to process all events                                                                                       |
| `make_plain_tree` (optional)   | normal root tree containing the candidates is written besides the default output tree in AnalysisTree format (set to `true`)                              |
| `write_detailed_bg` (optional) | optional detailed background information in output variable `generation` (see below)                                                                      |

The format of the detailed background information in the `generation` variable is as follows:
format for bg for 2-body-decay: -m12 0 d2 d1  
format for bg for 3-body-decay: -m13 m23 m12 d3 d2 d1  
with m12: mother of daughter 1 & 2 etc.  
code for daughters (d):  
1 - reco daughter is unmatched to mc  
2 - reco daughter is matched, but primary  
3 - reco daughter is secondary, produced not in decay from mother with not expected pdg  
4 - reco daughter is secondary, produced not in decay from mother with expected pdg  
5 - reco daughter is secondary, produced in decay from mother with not expected pdg  
6 - reco daughter is secondary, produced in decay from mother with expected pdg  
code for mothers (m):  
1 - daughters have the same mother  
2 - daughters have different mothers  
0 - at least one daughter does not have mother (e.g. primary)  

If detailed background is not selected or if mc match for mother is found, value of variable "generation" is:  
0 - no mc match for mother  
1 - mother is primary particle  
2 - mother is second generation etc.

### pid_mode
The method for particle identification is configured in the `pid_mode` variable and can be one of the following ones:
| `pid_mode` | Description                                                                                              |
| ---------- | -------------------------------------------------------------------------------------------------------- |
| `0`        | no pid (only charge information is used! optional daughers must be set to 1 / -1 for pos / neg pdgs!)    |
| `1`        | mc-pid                                                                                                   |
| `2`        | reconstructed TOF pid - default from Pid-framework                                                       |
| `3`        | reconstructed TOF pid - pdg with max. pdg-purity is selected, if pdg-purity >min. purity                 |
| `4`        | reconstructed TOF pid - pdg is selected, if pdg-purity > min. purity (pdg-specfic purities are possible) |

**for 2-4: Pid-framework needs to be applied first to include TOF-pid & probabilities into the analysistree input-file**

### pid_purity
This block is required in pid-mode 3 & 4 and defines the minium purities to be used:
   - for pid-mode 3: minimum purity is set for all pdgs (pdg-spefic purities need to be omitted)
   - for pid-mode 4: minimum purity can be set for all pdgs or specifically for every pdg (general purity will be overwritten, if pdg-specific purity is given)

### decays
The decay block contains the specific decays PFSimple shall reconstruct. Each decay is defined by several sub-blocks:

#### mother
The mother block must contain the fields `pdg_code` and `name` (chosen by user). You can optionally set the keys `mass` (Unit: GeV/c^2) and `mass_sigma` to be used when calculating mass decrepancies or setting nonlinear mass constrains, otherwise values from KFParticleDatabase are taken.
It is possible to define several cuts to be applied to the mother particle in the `cuts` subblock, see `at_interface/config/config_dummy.json` for the available JSON keys. For more information on specific cuts see `Constants.hpp`.

#### daughters
The daughters block is a list consisting of 2 to 3 subblocks, each defining a daugher particle of the mother. It is required to give a JSON list in the `pdg_codes` field which contains the of PDG codes to be considered as possible daughter candidates of the decay in the. Optionally, alternative pdgs can be considered or further options can be selected for the reconstruction:
- for background: set to "1" / "-1" for pos / neg daughters
- for tracks without TOF id: set to "2" / "-2" for pos / neg daughters (in pid-mode 2, 3, 4)
- In pid-mode 0 optional daughers "1" / "-1" for pos / neg pdgs must be added

Analogously to the mother you can also apply cuts for the daughters in the `cuts` subblock (`chi2prim` and `cos`).

#### daughter_combination_cuts (3-body-decays)
For 3-body decays, the list of available cuts is extended, amoung others, to cuts on the secondary mothers (SM) of combinations of 2 daughters. E.g. for H3L->p + pi + d, cuts on the SM of p-pi can be applied. The cuts that test the mother against the PV, e.g. chi2topo, are inverted for the SM to exclude SMs that come from the PV. They are defined in the optional `secondary_mother_cuts` block of the decay, which can contain for each of the 4 possible variables `chi2geo`, `cosopen`, `chi2topo` and `costopo` a list with up to 3 entries. The entries respectively reflect the cuts to be applied on the secondary mothers of daughter combinations 1 & 2, 1 & 3 and 2 & 3.

#### save_options
The optional `save_options` block can contain specific string values which enable features related to the saving of the reconstructed mother:
| Flag                  | Description                                                                                                                                                                               |
| --------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `mass_constraint`     | A nonlinear mass constraint can be applied on the mother particle which modifies the mother's state vector and covariance matrix to assure the 4-momentum is consitent with the pdg mass. |
| `transport_to_pv`     | The mother is transported to the PV before saving                                                                                                                                         |
| `do_not_write_mother` | The writing of the mother into the output-file (e.g. the lower-level mother in cascade decays) can get surpressed. But if this option is selected, it is not possible to peform the MC-matching for the upper-level mother. |

### output_cuts
The optional `output_cuts` block contains AnalysisTree field cuts to be applied to the output candidates before storing them to disk. One can either select a certain range for a cut variable with `from` and `to` or select a specific value with `equal`. It is possible to chain as many cuts as you want. 

### Cascade decays
The reconstruction of a cascade decay works by adding all consecutive decays in one config file. The order of the decays need to be from last generation to
first generation. Cascade decays with multiple stages as well as combinations of
2- and 3-body-decays can get reconstructed.

Several example parameter files can be found in the folder `at_interface/configs`.

## Example Configs
In the `at_interface/configs` directory are some example configs which have the following purpose:
| File                             | Particle(s)        | `pid_mode` | Comment                                                                                                                        |
|----------------------------------|--------------------|------------|--------------------------------------------------------------------------------------------------------------------------------|
| `config_dummy.json`              | Lambda             | `1`        | shows all available JSON fields, cuts are unphysical                                                                           |
| `config_lambda_pidmode1_v1.json` | Lambda             | `1`        | using KFParticle cuts                                                                                                          |
| `config_lambda_pidmode1_v2.json` | Lambda             | `1`        | using optimized cuts                                                                                                           |
| `config_lambda_pidmode0.json`    | Lambda             | `0`        |                                                                                                                                |
| `config_lambda_pidmode2.json`    | Lambda             | `2`        | using PID frameworks default PID mode                                                                                          |
| `config_lambda_pidmode1_v3.json` | Lambda             | `1`        | with very loose cuts needed for machine learning applications, saves only mc-true signal candidates using `save_options` block |
| `config_H3L_pidmode1.json`       | H3L                | `1`        | with MC-PID mode (mode 1), example for 3-body decay and secondary mother cuts                                                  |
| `config_H3L_pidmode4.json`       | H3L                | `4`        | shows application of minimum purities                                                                                          |
| `config_He5L_pidmode1.json`      | He5L               | `1`        | example for 3-body decay and secondary mother cuts                                                                             |
| `config_xi_pidmode1.json`        | Xi, Lambda         | `1`        | basic example for the the cascade decay feature                                                                                |
| `config_omega_pidmode1.json`     | Omega, Lambda      | `1`        | basic example for the the cascade decay feature                                                                                |
| `config_omegastar_pidmode1.json` | Omega*, Xi, Lambda | `1`        | advanced example for the the cascade decay feature including a 3-body decay and secondary mother cuts                          |



## First run
Each time before running the prepared executable you should set the environment variables to let your system know where to find libraries:
```
source /path-to-root/install/bin/thisroot.sh
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/path-to-pfsimple-installation/lib/
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/path-to-pfsimple-installation/external/lib/
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/path-to-analysistree-installation/lib/ #can be ommitted if PFSimple_BUNDLED_AT was set to ON
```

The reconstruction will be executed with the following command:
```
at_interface/main_json filelist.txt config.json [--output <file>] [--plain-output <file>] [--set <key>=<value>]
```
The binary accepts 2 positional arguments: The path of the input filelist.txt file and the config file in JSON format. It is also possible to modify the output file paths with the `--output` and `--plain-output` optional arguments as well as to override keys in the JSON config file by adding keys like `--set decays[0].mother.cuts.LdL=4.0`.

### filelist.txt
filelist.txt must be a text file with names (including paths) to AnalysisTree
files which you want to analyze with PFSimple. Each file name should
be on the next line, and the last symbol should also be a switch to
next line. The example of filelist is in the source directory of
PFSimple (filelist_example.txt)