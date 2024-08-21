#pragma once
#include <Eigen/Dense>
#include <fmt/format.h>
#include <string>
#include <map>
#include "string_tools.h"


using Eigen::MatrixXd;


// This class is only for the mapping matrix method
namespace MappingMatrix {

    // Read and store mapping matrix (each column holds the atoms weights of a bead)
    class Mapping {
        public:
            Mapping(std::string path);

            // Store the mapping matrix read from the file
            MatrixXd map_mat;
            // Mapping matrix atom type order
            std::vector<std::string> map_order;
            // Number of beads (columns of the mapping matrix)
            int num_beads;
            // C-H atom connection for a crsytal chain order
            std::map<int, std::vector<int>> connections;
            // Mapping matrix with no overlaped atoms
            MatrixXd bead_map;
        private:
            // Read the mapping matrix from the file
            void _read_mapping_matrix(std::string path);
            // Find the bead map
            void _find_bead_map();
    };
}
