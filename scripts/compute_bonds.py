#!/usr/bin/env python3
"""
Goal
----
Check the bond length between each two bonded pair and show the ones that are
longer than a specific threshold.
"""


import sys
from glob import glob
import numpy


def compute_bonds():
    """ Go through each bonded pair and calculate the bond length and compare """
    # Default values for upper and lower limit of wrong bond length range
    wrong_bond_length_upper_limit = 2.0
    wrong_bond_length_lower_limit = 1.0
    if len(sys.argv) == 1:
        files = glob("*.lammps")
    elif len(sys.argv) > 1  and len(sys.argv) <= 6:
        files = [sys.argv[1]]
        for i in range(2, len(sys.argv)):
            if sys.argv[i] == "-ul":
                wrong_bond_length_upper_limit = float(sys.argv[i+1])
            elif sys.argv[i] == "-ll":
                wrong_bond_length_lower_limit = float(sys.argv[i+1])
    else:
        print("You need to enter a LAMMPS data file and the wrong bond length.")
        sys.exit(0)
    
    # Masses of C and H for atomistic systems
    C_mass = 12.01120
    H_mass =  1.00797

    for f in sorted(files):
        masses = {}
        fid = open(f, 'r')
        for line in fid:
            args = line.split()
            if len(args) == 0 or "Coeffs" in line:
                continue
            if 'atoms' in line:
                num_atoms = int(args[0])
                atoms = numpy.zeros((num_atoms, 3))
                atom_types = numpy.zeros(num_atoms)
            elif "atom types" in line:
                num_atom_types = int(args[0])
            if 'xlo xhi' in line:
                xlo, xhi = [float(s) for s in args[:2]]
            elif 'ylo yhi' in line:
                ylo, yhi = [float(s) for s in args[:2]]
            elif 'zlo zhi' in line:
                zlo, zhi = [float(s) for s in args[:2]]
            if "Masses" in line:
                # Simulation box dimension
                dx, dy, dz = xhi - xlo, yhi - ylo, zhi - zlo
                next(fid)
                for line in fid:
                    args = line.split()
                    if len(args) == 0:
                        break
                    masses[float(args[1])] = int(args[0])
            # Read and store atom coordinates and types
            elif args[0] == "Atoms":
                next(fid)
                for line in fid:
                    cols = line.split()
                    if len(cols) == 0:
                        break
                    i = int(cols[0]) - 1
                    atom_types[i] = int(cols[2])
                    atoms[i] = [float(s) for s in cols[-3:]]
                    # Map atom inside the PBC
                    while atoms[i,0] < xlo:
                        atoms[i,0] += dx
                    while atoms[i,0] > xhi:
                        atoms[i,0] -= dx

                    while atoms[i,1] < ylo:
                        atoms[i,1] += dy
                    while atoms[i,1] > yhi:
                        atoms[i,1] -= dy

                    while atoms[i,2] < zlo:
                        atoms[i,2] += dz
                    while atoms[i,2] > zhi:
                        atoms[i,2] -= dz
            # Read and store the bonds
            elif args[0] == "Bonds":
                next(fid)
                bonds = []
                for line in fid:
                    cols = line.split()
                    if len(cols) == 0:
                        break
                    i, j = [int(s) for s in cols[-2:]]
                    bonds.append((i-1, j-1))

        # Ensure atoms are inside PBC
        assert all(numpy.min(atoms, axis=0) >= [xlo, ylo, zlo])
        assert all(numpy.max(atoms, axis=0) <= [xhi, yhi, zhi])

        # Compute bond length for each pair of atoms
        count = 0
        bond_lengths = {}
        for i, j in bonds:
            r = atoms[i,:] - atoms[j,:]
            # Find the PBC distance
            while r[0] > 0.5*dx:
                r[0] -= dx
            while r[0] < -0.5*dx:
                r[0] += dx

            while r[1] > 0.5*dy:
                r[1] -= dy
            while r[1] < -0.5*dy:
                r[1] += dy

            while r[2] > 0.5*dz:
                r[2] -= dz
            while r[2] < -0.5*dz:
                r[2] += dz
            # Bond length value
            bond_len = round(numpy.linalg.norm(r), 2)
            if bond_len > wrong_bond_length_upper_limit or bond_len < wrong_bond_length_lower_limit:
                print(f'({i+1}, {j+1}): {bond_len:.2f}')
                count += 1
	    # Type of bond length	
            bond_type = "{:d}-{:d}".format(int(atom_types[i]), int(atom_types[j]))
            if bond_type not in bond_lengths:
                bond_lengths[bond_type] = {}
            if bond_len not in bond_lengths[bond_type]:
                bond_lengths[bond_type][bond_len] = 0
            bond_lengths[bond_type][bond_len] += 1

        # Show the types of bond lengths, their values, and their number of repetition
        for bond_type in bond_lengths:
            for bl in bond_lengths[bond_type]:
                print("Number of bond type {} with length {:4.2f} = {}".format(\
                                   bond_type, bl, bond_lengths[bond_type][bl]))

        # Show number of wrong bond lengths
        print("------------------------------------------------------")
        print("Number of wrong bonds in {} = {}".format(f, count))


if __name__ == "__main__":
    compute_bonds()

