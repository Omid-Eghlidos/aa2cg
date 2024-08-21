#include "cg_dump.h"
#include <fmt/format.h>
#include <fmt/ostream.h>
#include <tuple>
#include <iostream>
#include <fstream>


// Mapping scheme method
namespace MappingScheme {

    // Class to store dumped coarse grained timesteps
    CG_dump::CG_dump(const LMP_data &data, const LMP_dump &dump, const CG_data &cgd) :
        _data(data), _dump(dump), _cg_data(cgd) {
        _coarsen_timesteps();
    }


    void CG_dump::_coarsen_timesteps () {
        cg_ts = std::vector<Timestep> (_dump.timesteps.size());
        for (size_t ts=0; ts<_dump.timesteps.size(); ts++) {
            cg_ts[ts].timestep = _dump.timesteps[ts];
            cg_ts[ts].box = _dump.box[ts];
            cg_ts[ts].coords = BeadsCoords (_cg_data.beads.size());
            for (auto bead : _cg_data.beads) {
                // Position of the first atom in the bead
                Vector3d x0 = _dump.atom_coords[ts].row(bead.atoms[0]);
                Vector3d bead_coords = x0 * bead.atom_weights[0];
                for (size_t i=1; i<bead.atoms.size(); i++) {
                    Vector3d x = _dump.atom_coords[ts].row(bead.atoms[i]);
                    Vector3d dx = x - x0;
                    Vector3d ds = (_dump.box[ts].inverse() * dx).array().round().matrix();
                    dx = dx - _dump.box[ts] * ds;
                    bead_coords += (x0 + dx) * bead.atom_weights[i];
                }
                cg_ts[ts].coords[bead.id] = bead_coords;
            }
        }
        fmt::print("Found coarse coordinates for {} timsteps.\n", cg_ts.size());
    }
}


// Mapping matrix method
namespace MappingMatrix {

    // Class to store dumped coarse grained timesteps
    CG_dump::CG_dump(const LMP_data &data, LMP_dump &dump, const Mapping &map,
                     const Reorder &blocks, const CG_data &cgd) :
        _data(data), _dump(dump), _mapping(map), _blocks(blocks), _cg_data(cgd) {
        _coarsen_timesteps();
    }


    void CG_dump::_coarsen_timesteps () {
        cg_ts = std::vector<Timestep> (_dump.timesteps.size());
        for (size_t t=0; t<_dump.timesteps.size(); t++) {
            cg_ts[t].timestep = _dump.timesteps[t];
            cg_ts[t].box = _dump.box[t];
            BlockAtoms chains_coords (_blocks.block_atoms.size());
            for (size_t i = 0; i < _blocks.block_atoms.size(); i++) {
                for (auto block : _blocks.block_atoms[i]) {
                    MatrixXd block_coords(block.rows(), 3);
                    // Map atoms into the periodic cell w.r.t first atom of the chain
                    _map_atoms_to_periodic_cell(block, t);
                    for (int j = 0; j < block.rows(); j++) {
                        block_coords.row(j) = _dump.atom_coords[t].row(block(j));
                    }
                    MatrixXd bead_coords = _mapping.map_mat.transpose() * block_coords;
                    // Adjust the beads coordinates to be inside the periodic cell
                    for (int k=0; k<bead_coords.rows(); k++) {
                        Vector3d ds = (_dump.box[t].inverse()
                                    * bead_coords.row(k)).array().round().matrix();
                        bead_coords.row(k) -= _dump.box[t] * ds;
                    }
                    chains_coords[i].push_back(bead_coords);
                }
            }
            cg_ts[t].coords = chains_coords;
        }
        fmt::print("Found coarse coordinates for {} timsteps.\n", cg_ts.size());
    }


    void CG_dump::_map_atoms_to_periodic_cell(MatrixXd block, int ts) {
        MatrixXd box = _dump.box[ts];
        for (int i=0; i<block.rows(); i++) {
            Vector3d x0 = _dump.atom_coords[ts].row(block(0));
            for (auto j : _data.bond_table[block(i)]) {
                if ((block.array() == j).any()) {
                    Vector3d x = _dump.atom_coords[ts].row(j);
                    Vector3d dx = x - x0;
                    // Find the nearest atom in the box with whole unit cell coordinates
                    Vector3d ds = (box.inverse() * dx).array().round().matrix();
                    // Find the nearest PBC image
                    dx -= box * ds;
                    // New position of the atom
                    x = x0 + dx;
                    _dump.atom_coords[ts].row(j) = x;
                }
            }
        }
    }
}

