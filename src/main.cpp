#include "inputs.h"
#include "lammps_data.h"
#include "lammps_dump.h"
#include "mapping_matrix.h"
#include "blocks.h"
#include "reorder.h"
#include "cg_data.h"
#include "cg_dump.h"
#include "outputs.h"


int main(int argc, char **argv) {
    Eigen::initParallel();
    // Take the path to the input file
    std::string input_file;
    if (argc != 2) {
        std::cerr << "Input file is not entered.\n"
                     "Try ./aa2cg <input_file>\n";
        return 1;
    }
    input_file = argv[1];
    // Store the input settings
    Inputs inputs(input_file);
    // Read lammps format atomistic data file
    LMP_data lmp_data(inputs);
    // Read lammps format atomistic dump file
    LMP_dump lmp_dump(inputs.lmp_dump);
    // Coarse graining using the mapping scheme method
    if (inputs.mapping_matrix.empty()) {
        // Generate the lammps format coarse-grained data file
        MappingScheme::CG_data cg_data(inputs, lmp_data);
        // Generate the lammps format coarse-grained dump file
        MappingScheme::CG_dump cg_dump(lmp_data, lmp_dump, cg_data);
        // Writing the outputs
        MappingScheme::Outputs outputs(inputs, lmp_data, lmp_dump, cg_data, cg_dump);
    }
    // Coarse grained using mapping matrix method
    else {
        // Read and store the mapping matrix
        MappingMatrix::Mapping mapping(inputs.mapping_matrix);
        // Find chains and their SMILE blocks
        MappingMatrix::Blocks blocks(inputs, lmp_data);
        // Match the order of the blocks with the mapping matrix
        MappingMatrix::Reorder reordered_blocks(inputs, lmp_data, mapping, blocks);
        // Generate the lammps format coarse-grained data file
        MappingMatrix::CG_data cg_data(inputs, lmp_data, mapping, reordered_blocks);
        // Generate the lammps format coarse-grained dump file
        MappingMatrix::CG_dump cg_dump(lmp_data, lmp_dump, mapping
                                               , reordered_blocks, cg_data);
        // Writing the outputs
        MappingMatrix::Outputs outputs(inputs, lmp_data, lmp_dump, mapping
                                     , reordered_blocks, cg_data, cg_dump);
    }
}

