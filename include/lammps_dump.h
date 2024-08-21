#pragma once
#include <Eigen/Dense>
#include <string>
#include <vector>

using Eigen::MatrixXd;
using Eigen::Vector3d;

// Reads and stors the LAMMPS dump file.
class LMP_dump {
    public:
        LMP_dump(std::string path);

        // Number of atoms in the system
        int num_atoms;
        // Store the timestep number
        std::vector<int> timesteps;
        // Store the calculated simulation box dimension in each timestep
        std::vector<MatrixXd> box;
        // Store the box bounds of each timestep from the dump file
        std::vector<MatrixXd> bounds;
        // Determine the type of the coordinates in the dump file
        std::string coords_type;
        // Store the atom coordinates of each timestep
        std::vector<MatrixXd> atom_coords;
        // Return true if the simulation box is triclinic
        bool triclinic = false;
    private:
        // Holds the box for each timestep
        MatrixXd _timestep_box;
};

