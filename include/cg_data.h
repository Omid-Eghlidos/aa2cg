#pragma once
#include <bits/stdc++.h>
#include <Eigen/Dense>
#include <string>
#include <vector>
#include "inputs.h"
#include "lammps_data.h"
#include "lammps_dump.h"
#include "mapping_matrix.h"
#include "reorder.h"


using Eigen::MatrixXd;
using Eigen::Vector3d;


typedef std::vector<std::vector<MatrixXd>> BlockAtoms;


// Mapping scheme
namespace MappingScheme {

    // Stores parameters of the beads
    struct Bead {
        // Index of the bead
        int id;
        // Chain number that the bead belongs to.
        int chain;
        // Type of the bead
        int type;
        // Bead's alphabetic name
        char name;
        // Mass of the bead
        double mass = 0.0;
        // Charge of the bead as sum of all the atom charges inside of it
        double charge = 0.0;
        // List of atoms inside the bead
        std::vector<int> atoms;
        // Bead mapping weights
        std::vector<double> atom_weights;
        // Coordinates of the bead (can be unwrapped, scaled, or unscaled)
        Vector3d coords = Vector3d::Zero(3);
    };


    // Class to store specifications of the coarse grained system
    class CG_data {
        public:
            CG_data(const Inputs &inputs, LMP_data &data);
            // Determine the bead type considering any existing symmetry in the mapping
            int get_bead_type(int type) const;
            // Find the bond type for two beads among the existing bond types
            int get_bond_type(int i, int j) const;
            // Find the angle type for three beads among the existing angle types
            int get_angle_type(int i, int j, int k) const;
            // Find the torsion type for four beads among the existing torsion types
            int get_torsion_type(int i, int j, int k, int l) const;
            // Find the improper type for four beads among the existing improper types
            int get_improper_type(int i, int j, int k, int l) const;

            // Beads defined according to the mapping matrix
            std::vector<Bead> beads;
            // Store the unique bead types
            std::vector<int> bead_types;
            // Bead types masses
            std::map<int, double> bead_types_masses;
            // Store the bead bonds
            std::vector<std::vector<int>> bead_bonds;
            // Store the unique bead bond types
            std::vector<std::tuple<int, int>> bead_bond_types;
            // Store the bead angles
            std::vector<std::tuple<int, int, int>> bead_angles;
            // Store the unique bead angle types
            std::vector<std::tuple<int, int, int>> bead_angle_types;
            // Store the bead torsions
            std::vector<std::tuple<int, int, int, int>> bead_torsions;
            // Store the unique bead torsion types
            std::vector<std::tuple<int, int, int, int>> bead_torsion_types;
            // Store the bead impropers
            std::vector<std::tuple<int, int, int, int>> bead_impropers;
            // Store the unique bead improper types
            std::vector<std::tuple<int, int, int, int>> bead_improper_types;
        private:
            // Make beads based on the input definition
            void _make_beads();
            // Compute the beads mass to conserve atomistic momentum
            double _compute_bead_mass(std::vector<int> atoms, std::vector<double> weights);
            // Find beads for each defined type
            void _find_beads_of_type(char bead_name, BeadType bead);
            // Find the bead coordinates for each bead
            Vector3d _find_beads_coords(std::vector<int> atoms, std::vector<double> weights);
            // Find the bonds between the beads
            void _identify_beads_bonds();
            // Determine unique bond types for bead bonds
            void _bond_types(int i, int j);
            // Determine if two beads are bonded together
            bool _are_bonded(std::vector<int> bead1_atoms, std::vector<int> bead2_atoms);
            // Find the angle of the beads
            void _identify_beads_angles();
            // Determine unique angle types for bead angles
            void _angle_types(int i, int j, int k);
            // Find the torsion angle of the beads
            void _identify_beads_torsions();
            // Determine unique torsion types for bead torsions
            void _torsion_types(int i, int j, int k, int l);
            // Find the improper angle of the beads
            void _identify_beads_impropers();
            // Determine unique improper types for bead torsions
            void _improper_types(int i, int j, int k, int l);

            // Input settings
            const Inputs &_inputs;
            LMP_data &_data;
            std::vector<bool> _visited;
    };
}


// Mapping matrix method
namespace MappingMatrix {

    // Stores parameters of the beads
    struct Bead {
        // Index of the bead
        int id;
        // Chain number that the bead belongs to.
        int chain;
        // Type of the bead
        int type;
        // Charge of the bead as sum of all the atom charges inside of it
        double charge;
        // Mass of the bead as the sum of all the atom charges inside of it
        double mass;
        // List of atoms inside the bead
        std::vector<int> atoms;
        // Coordinates of the bead (can be unwrapped, scaled, or unscaled)
        Vector3d coords = Vector3d::Zero(3);
    };


    // Class to store specifications of the coarse grained system
    class CG_data {
        public:
            CG_data(const Inputs &inputs, const LMP_data &data, const Mapping &map
                  , const Reorder blocks);
            // Determine the bead type considering any existing symmetry in the mapping
            int get_bead_type(int type) const;
            // Find the bond type for two beads among the existing bond types
            int get_bond_type(int i, int j) const;
            // Find the angle type for three beads among the existing angle types
            int get_angle_type(int i, int j, int k) const;
            // Find the torsion type for four beads among the existing torsion types
            int get_torsion_type(int i, int j, int k, int l) const;
            // Find the improper type for four beads among the existing improper types
            int get_improper_type(int i, int j, int k, int l) const;

            // Beads defined according to the mapping matrix
            std::vector<Bead> beads;
            // Bead types masses
            std::map<int, double> bead_types_masses;
            // Store the unique bead types
            std::vector<int> bead_types;
            // Coordinates of the beads
            BlockAtoms beads_coords;
            // Store the bead bonds
            std::vector<std::vector<int>> bead_bonds;
            // Store the unique bead bond types
            std::vector<std::tuple<int, int>> bead_bond_types;
            // Store the bead angles
            std::vector<std::tuple<int, int, int>> bead_angles;
            // Store the unique bead angle types
            std::vector<std::tuple<int, int, int>> bead_angle_types;
            // Store the bead torsions
            std::vector<std::tuple<int, int, int, int>> bead_torsions;
            // Store the unique bead torsion types
            std::vector<std::tuple<int, int, int, int>> bead_torsion_types;
            // Store the bead impropers
            std::vector<std::tuple<int, int, int, int>> bead_impropers;
            // Store the unique bead improper types
            std::vector<std::tuple<int, int, int, int>> bead_improper_types;
        private:
            // Make a bead based on the mapping matrix
            void _make_beads();
            // Compute the beads mass to conserve atomistic momentum
            double _compute_bead_mass(int i);
            // Compute bead charges to conserve system's charge balance
            double _compute_bead_charge(int i);
            // Find the coordinates of the beads by multiplying mapping to atomic coords
            void _find_beads_coords();
            // Find the bonds between the beads
            void _identify_beads_bonds();
            // Determine unique bond types for bead bonds
            void _bond_types(int i, int j);
            // Determine if two beads are bonded together
            bool _are_bonded(std::vector<int> bead1_atoms, std::vector<int> bead2_atoms);
            // Find the angle of the beads
            void _identify_beads_angles();
            // Determine unique angle types for bead angles
            void _angle_types(int i, int j, int k);
            // Find the torsion angle of the beads
            void _identify_beads_torsions();
            // Determine unique torsion types for bead torsions
            void _torsion_types(int i, int j, int k, int l);
            // Find the improper angle of the beads
            void _identify_beads_impropers();
            // Determine unique improper types for bead torsions
            void _improper_types(int i, int j, int k, int l);

            // Input settings
            const Inputs &_inputs;
            const LMP_data &_data;
            const Mapping &_mapping;
            Reorder _blocks;
    };
}

