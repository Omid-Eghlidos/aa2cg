#pragma once
#include <Eigen/Dense>
#include <string>
#include "inputs.h"
#include "lammps_data.h"
#include "mapping_matrix.h"
#include "blocks.h"


using Eigen::Vector3d;


typedef std::vector<std::vector<MatrixXd>> BlockAtoms;
typedef std::vector<std::vector<std::vector<std::string>>> BlockOrder;
typedef std::vector<Vector3d> Coords;


// This class is only for the mapping matrix method
namespace MappingMatrix {

    class Reorder {
        public:
            Reorder(const Inputs &inputs, const LMP_data &data
                   , const Mapping &map, const Blocks &blocks);

            // Store the blocks synced with the order of the mapping matrix
            BlockAtoms block_atoms;
            // Stores atom coordinates of blocks of each chain
            BlockAtoms block_coords;
            // Number of reordered blocks
            int num_blocks = 0;
        private:
            // Match the order of the blocks with the mapping matrix
            void _find_blocks_orders(bool amorphous);
            // Function to comapre and reorder the blocks
            void _reorder_blocks();
            // Find the type of the H atoms in the block
            std::string _get_hydrogen_type(int c, int h);
            // Find the methy H atoms
            Coords _methyl_hydrogens(Coords c, Coords h);
            // Find the gemini H atoms
            Coords _gemini_hydrogens(Coords c, Coords h);
            // Find the hydrogens order
            std::vector<int> _hydrogens_order(Coords c, Coords h);
            // Find the wrap bond vector inside the simulation box
            Vector3d _wrap_vector(Vector3d bv);

            // Input settings
            const Inputs &_inputs;
            const LMP_data &_data;
            const Mapping &_map;
            const Blocks &_blocks;

            // Hydrogen bond length
            double _lh = 1.10;
            // Keep the type order of atoms inside each block
            BlockOrder _blocks_orders;
    };
}

