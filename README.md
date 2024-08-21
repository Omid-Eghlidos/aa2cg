# aa2cg

**A**ll **A**tomic system to **C**oarse-**G**raining system.

This program maps a classical all-atomistic (AA) molecular system to it corresponding coarse-grained (CG) system using the following two methods:

1) Mapping matrix - implemented only for isotactic polypropylene (iPP)
2) Mapping scheme

The mapping matrix will be a separate input text file, while mapping scheme will be defined in the input text file.

### Mapping Matrix

In this method, CG configuration is obtained via multiplying a given mapping matrix defined for a chain inside of one unit-cell of the alpha-modification of the iPP crystal.
The columns of the mapping matrix denotes number of atoms in a single chain, which is 9 carbons (C) and 18 hydrogens (H) for a total of 27 atoms, and its rows denote number of CG sites, or beads, in the corresponding CG system, while its elements are weights of atoms inside of each bead.
In order to apply the mapping matrix to both amorphous and crystaline phases of an AA system, the following strategies were employed:
* In the amorphous phase, similar chain patterns, as defined by Simplified Molecular Input Line Entry System (SMILES), are found and multiplied by the mapping matrix.
* In the crystalline phase, similar chain patterbs are found using the crystal lattice symmetry.

### Mapping Scheme

In this method, CG configuration is obtained by finding molecular patterns defined by a scheme for bonded atoms and their weights.
The mapping scheme determine what bonded atoms to be lumped together and what is their wight inside of each bead.
The position and the weight of each bead will be the weighted average of the position and sum of the weights of its involved atoms, as defined by the scheme.


## Prerequisites

The following are required to deploy the code:

1) libeigen3-dev \n
2) cmake
3) gcc-9 or newer
4) Python 3 (for running scripts)

## Deployment Guideline

1) Get the libraries

    - cd lib
    - ./get_fmt.sh

2) Create a build folder

    - mkdir release
    - cd release

3) Compile

    - cmake -DCMAKE_BUILD_TYPE=Release ..
    - make -j4

4) Deploy (inside the release folder)

    - ./aa2cg [input file] [mapping_matrix]

### Inputs

Depending on the method chosen in the input file needs the followings:

    - **input file** : Path to the input file (.ini) comprising settings.
    - **mapping_matrix**: Path to the input mapping matrix text file, if this method is used.

