#include "inputs.h"
#include "string_tools.h"
#include <string>
#include <vector>
#include <iostream>


Inputs::Inputs(std::string path) {
    std::fstream fid(path);
    if (!fid) {
        std::cerr << "\nThe input file does not exist.\n";
        exit(0);
    }
    while (fid) {
        auto row = split(read_line(fid));
        if (row.empty()) {
            continue;
        }
        if (row[0] == "lammps_data") {
            lmp_data = row[1];
        }
        else if (row[0] == "lammps_dump") {
            lmp_dump = row[1];
        }
        // System settings
        else if (row[0] == "phase") {
            if (row[1] == "crystal")
                amorphous = false;
        }
        else if (row[0] == "crystal_periodicity" && !amorphous) {
            if (row[1] == "finite") {
                finite_crystal = true;
            }
        }
        else if (row[0] == "symmetrical_beads") {
            std::vector<int> symmetric;
            for (size_t i=1; i<row.size(); i++) {
                symmetric.push_back(from_string<int>(row[i]) - 1);
            }
            symmetrical_beads.push_back(symmetric);
        }
        // Mapping scheme settings (if given)
        else if (row[0] == "atom_type") {
            atom_types[row[2]] = from_string<int>(row[1]);
        }
        else if (row[0] == "bead_type") {
            BeadType bead;
            for (size_t i=2; i<row.size(); i++) {
                bead.atom_types.push_back(atom_types.find(row[i])->second);
            }
            bead.type = bead_types.size();
            bead_types[row[1][0]] = bead;
        }
        else if (row[0] == "bead_weights") {
            for (size_t i=2; i<row.size(); i++) {
                bead_types[row[1][0]].weights.push_back(from_string<double>(row[i]));
            }
        }
        // Mapping matrix settings (if given)
        else if (row[0] == "mapping_matrix") {
            mapping_matrix = row[1];
        }
        else if (row[0] == "type") {
            for (size_t i=2; i<row.size(); i++) {
                elements[from_string<int>(row[i])] = row[1];
            }
        }
        else if (row[0] == "smiles_block") {
            smiles_block = row[1];
        }
        // Output settings
        else if (row[0] == "atom_style") {
            atom_style = row[1];
        }
        else if (row[0] == "tag") {
            tag = row[1];
        }
    }
}

