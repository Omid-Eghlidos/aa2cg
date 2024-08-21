#pragma once
#include <string>
#include <vector>
#include <map>


// Store bead definitions if a mapping scheme is defined and not a mapping matrix
struct BeadType {
    // Bead type
    int type;
    // Atom types that make up the bead
    std::vector<int> atom_types;
    // CG weights for each atom inside the bead
    std::vector<double> weights;
};


// Class to store the input parameters
class Inputs {
    public:
        Inputs(std::string path);

        // Input files
        std::string lmp_data;
        std::string lmp_dump;
        // System settings
        bool amorphous = true;
        bool finite_crystal = false;
        // Mapping scheme settings (type definitions)
        std::map<std::string, int> atom_types;
        std::map<char, BeadType> bead_types;
        // Mapping matrix settings and possible symmetrical beads
        std::string mapping_matrix;
        std::vector<std::vector<int>> symmetrical_beads;
        std::map<int, std::string> elements;
        std::string smiles_block;
        // Output settings
        std::string atom_style = "full";
        std::string tag = "CG";
};

