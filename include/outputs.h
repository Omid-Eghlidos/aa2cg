#pragma once
#include <fmt/format.h>
#include <Eigen/Dense>
#include <tuple>
#include <fstream>
#include "fmt/ostream.h"
#include "inputs.h"
#include "lammps_data.h"
#include "lammps_dump.h"
#include "mapping_matrix.h"
#include "cg_data.h"
#include "cg_dump.h"
#include "reorder.h"


using Eigen::MatrixXd;


namespace MappingScheme {

    // Mapping scheme method class to generate the outputs
    class Outputs {
        public:
            // Constructor
            Outputs(const Inputs &inputs, const LMP_data &data, const LMP_dump &dump,
                    const CG_data &cg_data, const CG_dump &cg_dump);
        private:
            // Write coarse grained data into a file
            void _write_cg_data() const;
            // Write coarse grained timesteps into a file
            void _write_cg_dump() const;

            // Input settings
            const Inputs &_inputs;
            const LMP_data &_data;
            const LMP_dump &_dump;
            const CG_data &_cg_data;
            const CG_dump &_cg_dump;
    };
}


// Mapping matrix method
namespace MappingMatrix {

    // Mapping matrix method class to generate the outputs
    class Outputs {
        public:
            // Constructor
            Outputs(const Inputs &inputs, const LMP_data &data, const LMP_dump &dump,
                    const Mapping &map, const Reorder &blocks, const CG_data &cg_data,
                    const CG_dump &cg_dump);
        private:
            // Write coarse grained data into a file
            void _write_cg_data() const;
            // Write coarse grained timesteps into a file
            void _write_cg_dump() const;

            // Input settings from input file
            const Inputs &_inputs;
            const LMP_data &_data;
            const LMP_dump &_dump;
            const Mapping &_mapping;
            const Reorder &_blocks;
            const CG_data &_cg_data;
            const CG_dump &_cg_dump;
    };
}

