#pragma once
#include <Eigen/Dense>
#include <iostream>
#include <string>
#include <vector>
#include "lammps_data.h"
#include "lammps_dump.h"
#include "inputs.h"

using Eigen::MatrixXd;

typedef std::vector<std::vector<MatrixXd>> BlockAtoms;


// This class is only for the mapping matrix method
namespace MappingMatrix {

    class Blocks {
        public:
            Blocks(const Inputs &i, LMP_data &d);

            // Stores atom ids of blocks of each chain
            BlockAtoms block_atoms;
            // Stores atom coordinates of blocks of each chain
            BlockAtoms block_coords;
            // Number of patterns/blocks found
            int num_blocks = 0;
        private:
            // Find the block atoms and their coordinates
            void _find_blocks();
            // Determine the chain number that the pattern belongs to
            int _block_chain_number(const std::vector<int> &pattern_atoms);
            // Find the C atoms of the SMILES pattern
            std::vector<int> _find_SMILES_atoms(int atom);
            // Make SMILES pattern from the atoms found
            std::string _make_SMILES(std::vector<int> atoms);

            // SMILES pattern type order
            std::vector<int> _pattern_order;
            // Atoms already visited
            std::vector<bool> _visited;
            const Inputs _inputs;
            LMP_data _data;
    };
}

