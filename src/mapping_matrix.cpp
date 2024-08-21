#include "mapping_matrix.h"
#include <iostream>


using Eigen::MatrixXd;


// This class is only for the mapping matrix method
namespace MappingMatrix {

    // Mapping constructor
    Mapping::Mapping(std::string path) {
        // Find the mapping matrix dimension to initialize
        std::fstream fid(path);
        num_beads = 0;
        int num_chain_atoms = 0;
        while (fid) {
            auto row = split(read_line(fid));
            if (row.empty() || row[0] == "Atom") {
                continue;
            }
            else if (row[0] == "id") {
                row = split(read_line(fid));
                num_beads = row.size() - 2;
            }
            else {
                num_chain_atoms++;
            }
        }
        map_mat = MatrixXd::Zero(num_chain_atoms, num_beads);

        // Read the mapping matrix from the file and store it
        _read_mapping_matrix(path);
        // Cretae the C - H connections in a crystal chain (for iPP)
        connections[0] =  {9, 10, 11}; connections[1] = {12}; connections[2] = {13, 14};
        connections[3] = {15, 16, 17}; connections[4] = {18}; connections[5] = {19, 20};
        connections[6] = {21, 22, 23}; connections[7] = {24}; connections[8] = {25, 26};
        connections[9] = {0}; connections[10] = {0}; connections[11] = {0};
        connections[12] = {1}; connections[13] = {2}; connections[14] = {2};
        connections[15] = {3}; connections[16] = {3}; connections[17] = {3};
        connections[18] = {4}; connections[19] = {5}; connections[20] = {5};
        connections[21] = {6}; connections[22] = {6}; connections[23] = {6};
        connections[24] = {7}; connections[25] = {8}; connections[26] = {8};
        // Find a mapping with beads have no overlaped atoms
        _find_bead_map();
    }


    void Mapping::_read_mapping_matrix(std::string path) {
        // Read the mappping matrix input file with its first column as the
        // atom type order of one crystal chain and other columns the weights of
        // atoms in each bead
        std::fstream fid(path);
        if (!fid) {
            std::cerr << "\nMapping matrix file does not exist.\n";
            exit(0);
        }
        int counter = 0;
        while (fid) {
            auto row = split(read_line(fid));
            if (row.empty() || row[0] == "Atom" || startswith(row[0], "-") || row[0] == "id") {
                continue;
            }
            else {
                map_order.push_back(row[1]);
                for (int b=0; b<map_mat.cols(); b++) {
                    map_mat(counter, b) = from_string<double>(row[b+2]);
                }
                counter++;
            }
        }
        fmt::print("Read the mapping matrix from {}.\n", path);
    }


    void Mapping::_find_bead_map() {
        // Overlaped atom will be put in the bead which the atom has the max weight
        bead_map = MatrixXd::Zero(map_mat.rows(), map_mat.cols());
        for (int i=0; i<map_mat.rows(); i++) {
            auto max_wt = map_mat.row(i).maxCoeff();
            for (int j=0; j<map_mat.cols(); j++) {
                if (map_mat(i, j) == max_wt) {
                    bead_map(i, j) = max_wt;
                }
            }
        }
        // Normalize beads weights
        for (int j=0; j<bead_map.cols(); j++) {
            bead_map.col(j) /= bead_map.col(j).sum();
        }
    }
}

