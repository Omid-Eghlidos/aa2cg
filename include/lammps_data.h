#pragma once
#include <Eigen/Dense>
#include <string>
#include <vector>
#include <map>
#include "inputs.h"


using Eigen::MatrixXd;
using Eigen::VectorXd;


// Reads and stores LAMMPS data file.
class LMP_data {
    public:
        LMP_data(const Inputs &inputs);
        // Delete bonds that cross the simulation cell (for periodic crystal).
        void remove_bonds_across_pbc();
        // Connect bonds that cross the simulation cell (for periodic crystal).
        void connect_bonds_across_pbc();
        // Find the mass for each atom type
        double get_atom_mass(int i) const {
            return atom_masses.find(i)->second;
        }
        // Find the defined atom types (e.g., type 2 is C)
        std::string element(int i) const;
        // Find number of hydrogens bonded to a specific atom
        int num_bonded_hydrogens(int atom) const;

        // Store the initial simulation box
        MatrixXd box = MatrixXd::Zero(3, 3);
        // Store the initial coordinates of the atoms
        MatrixXd atom_coords;
        // Store the molecule number of each atom
        std::vector<int> atom_mols;
        // Store the types of each atom
        std::vector<int> atom_types;
        // Store the charge of each atom
        std::vector<double> atom_charges;
        // Store the mass of each atom type
        std::map<int, double> atom_masses;
        // Store the bonds of atoms by using a bond table
        std::vector<std::vector<int>> bond_table;
        // Store atoms of each chains
        std::vector<std::vector<int>> chains;
        // Removed bonds
        std::map<int, int> removed_bonds;
    private:
        // Find the chains and their atoms
        void _find_chains();

        // Input settings
        const Inputs &_inputs;
};

