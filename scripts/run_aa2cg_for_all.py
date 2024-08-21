#!/usr/bin/env python3
"""
Goal
----
Run aa2cg code to get the CG system for all the atomistic system. It will read
the LAMMPS format data and trajectory files from the /data and /dump folders,
respectively. Then run the aa2cg for each pair of them.

Dependencies
------------
For mapping matrix method, the mapping matrix file is required.

NOTE: The template input script at the bottom of this program needs to be
changed according to the system and the mapping scheme/matrix.

Inputs
------
1) LAMMPS format data files: data/*.lammps
2) LAMMPS format trajectory files: dump/*.lammpstrj

Deployment
----------
./run_aa2cg_for_all [Method]
    * Method:   1) For mapping scheme method enter '-ms'
                2) For mapping matrix method enter '-mm'
    * Options: To clean files from previous run use '--clean'
"""


import os
import sys
from glob import glob
import subprocess


def main():
    """ Run aa2cg code for all the atomistic systems """
    # Remove any files from previous run
    # Choose the method of running aa2cg from the command line
    # Methods: 1) mapping scheme (-ms) <default>, 2) mapping matrix (-mm)
    method = "-ms"
    if len(sys.argv) == 2 and "--clean" not in sys.argv:
        method = sys.argv[1]
    elif len(sys.argv) > 2 and "--clean" in sys.argv:
        old_files = []
        old_files.append(glob("aa2cg_*.ini"))
        old_files.append(glob("*.lammps"))
        old_files.append(glob("*.lammpstrj"))
        for f in old_files:
            for i in f:
                os.remove(i)
        sys.exit(0)
    else:
        print("Wrong number of inputs.")
        print("Try: ./run_aa2cg_for_all.py --clean  -> Remove files from previous run")
        print("Try: ./run_aa2cg_for_all.py -ms -> Mapping scheme method")
        print("Try: ./run_aa2cg_for_all.py -mm -> Mapping matrix method")
        sys.exit(0)

    # Make the input files for aa2cg code and write them into a file
    systems = make_aa2cg_input_files(method)
    # Use the created input files and run the aa2cg code for each system
    run_aa2cg_for_all(systems)


def make_aa2cg_input_files(method):
    data_extension = ".lammps"
    dump_extension = ".lammpstrj"
    # Read data and dump files from their folders
    if len(glob("data/*" + data_extension)) != 0:
        data_files = sorted(glob("data/*" + data_extension))
    else:
        print("Data files do not exist.")
        exit()
    if len(glob("dump/*" + dump_extension)) != 0:
        dump_files = sorted(glob("dump/*" + dump_extension))
    else:
        print("Dump files do not exist.")
        exit()

    # Store each systems data and dump file and its input file
    systems = []
    # Pair data and dump file for each pair of the atomistic system
    for data in data_files:
        system_name = data.split('/')[-1].split('.')[0]
        system_no = int(system_name.split('_')[-1])
        dump = "dump/" + system_name + dump_extension
        if dump in dump_files:
            system_data = "./data/" + system_name + data_extension
            system_dump = "./dump/" + system_name + dump_extension
            system_output = "cg_" + system_name
            input_path = "aa2cg_{}.ini".format(system_no)
            system_input = open(input_path, 'w')
            # Mapping scheme method input
            if method == "-ms":
                input_template = mapping_scheme_template
            # Mapping matrix method input
            elif method == "-mm":
                input_template = mapping_matrix_template
            # Write into a file for each system
            system_input.write(input_template.format(system_data=system_data
                                                   , system_dump=system_dump
                                                   , system_output=system_output))
            systems.append("aa2cg_{}.ini".format(system_no))
        else:
            print("Trajectory for the {} system is missing.\n".format(system_name))
    return systems


def run_aa2cg_for_all(systems):
    """ Run the aa2cg code for all the systems """
    # Bash command to run the aa2cg code
    aa2cg = "/home/oeghlido/development/aa2cg/release/aa2cg"
    for sys in sorted(systems):
        sys_no = sys.split('_')[1].strip(".ini")
        print("\n##################################")
        print("Run aa2cg code for system {} ...".format(sys_no))
        print("##################################\n")
        process = subprocess.run([aa2cg, sys], check=True, text=True
                                , stdout=subprocess.PIPE)
        print(process.stdout)


def make_directory(dir_name):
    """ Make directory in the current working directory if it does not exist """
    dir_path = os.path.join(os.getcwd(), dir_name)
    try:
        os.mkdir(dir_path)
        exist = False
    except OSError:
        exist = True
    return dir_path


mapping_scheme_template = """
################### Inputs ########################
# Path to LAMMPS data file (*.lammps)
lammps_data     {system_data}

# Path to LAMMPS trajectory (dump) file (*.lammpstrj)
lammps_dump     {system_dump}

################ System Settings ###################
# System phase (crystal, amorphous(default))
phase 	amorphous

# Input file atom types definition
atom_type 	1	hc
atom_type	2	c3
atom_type	3	c2

# Bead types definition
bead_type	    A	c3 c2
bead_weights	A	0.5 0.5
bead_type	    B	c2 c2
bead_weights	B	0.5 0.5

# Symmetrical beads in the mapping scheme
symmetrical_beads	1  2

################### Outputs ########################
# Outputs tag ("CG" is default)
tag     {system_output}

# Data file atom_style (molecular, full(default))
atom_style 	molecular

"""


mapping_matrix_template = """
################### Inputs ########################
# Path to LAMMPS data file (*.lammps)
lammps_data     {system_data}

# Path to LAMMPS trajectory (dump) file (*.lammpstrj)
lammps_dump     {system_dump}

################ Mapping matrix #####################
# Path to mapping matrix file (*.txt)
mapping_matrix  mapping_matrix_1.txt

# Symmetrical beads in the mapping matrix
symmetrical_beads	1  2  3

################# System Settings ##################
# System phase (crystal , amorphous(default))
phase 	amorphous

# Simulation box periodicity of the crystal phase (finite/infinite)
#crystal_periodicity infinite

# Atom types (what are atom types in the data file)
# amorphous : C 1 2 4 -&- H 3
type    C   2 3 4
type    H   1

# SMILE block pattern
# amorphous (try CC(C)CC(C)CC or C(C)CC(C)CC(C)C)
smiles_block CC(C)CC(C)CC(C)

################### Outputs ########################
# Outputs tag ("CG" is default)
tag     {system_output}

# Data file atom_style (molecular, full(default))
atom_style 	molecular

"""


if __name__ == "__main__":
    main()

