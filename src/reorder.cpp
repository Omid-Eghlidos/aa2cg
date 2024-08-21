#include "reorder.h"
#include <fmt/format.h>
#include <unordered_set>
#include <math.h>
#include "string_tools.h"


// This class is only for the mapping matrix method
namespace MappingMatrix {

    // Constructor
    Reorder::Reorder(const Inputs &i, const LMP_data &d, const Mapping &m
        , const Blocks &b) : _inputs(i), _data(d), _map(m), _blocks(b) {
        block_atoms = BlockAtoms(_data.chains.size());
        block_coords = BlockAtoms(_data.chains.size());
        _blocks_orders = BlockOrder(_blocks.block_atoms.size());
        _find_blocks_orders(_inputs.amorphous);
        _reorder_blocks();
    }


    void Reorder::_reorder_blocks() {
        for (size_t bo = 0; bo < _blocks_orders.size(); bo++) {
            for (size_t b = 0; b < _blocks_orders[bo].size(); b++) {
                MatrixXd block = _blocks.block_atoms[bo][b];
                MatrixXd coords = MatrixXd::Zero(27, 3);
                // Carbons
                for (size_t i = 0; i < 9; i += 3) {
                    for (int j = 0; j < 3; j++) {
                        if (_blocks_orders[bo][b][i+j] != _map.map_order[i+j]) {
                            for (int k = 0; k < 3; k++) {
                                if (k != j && _blocks_orders[bo][b][i+k] ==
                                              _map.map_order[i+j]) {
                                    block(i+j) = _blocks.block_atoms[bo][b](i+k);
                                    coords.row(i+j) = _blocks.block_coords[bo][b].row(i+k);
                                }
                            }
                        }
                        else {
                            block(i+j) = _blocks.block_atoms[bo][b](i+j);
                            coords.row(i+j) = _blocks.block_coords[bo][b].row(i+j);
                        }
                    }
                }
                // Hydrogens
                for (size_t i = 9; i < 27; i += 6) {
                    for (int j = 0; j < 6; j++) {
                        if (_blocks_orders[bo][b][i+j] != _map.map_order[i+j]) {
                            for (int k = 0; k < 6; k++) {
                                if (k != j && _blocks_orders[bo][b][i+k] ==
                                              _map.map_order[i+j]) {
                                    block(i+j) = _blocks.block_atoms[bo][b](i+k);
                                    coords.row(i+j) = _blocks.block_coords[bo][b].row(i+k);
                                }
                            }
                        }
                        else {
                            block(i+j) = _blocks.block_atoms[bo][b](i+j);
                            coords.row(i+j) = _blocks.block_coords[bo][b].row(i+j);
                        }
                    }
                }
                block_atoms[bo].push_back(block);
                block_coords[bo].push_back(coords);
                num_blocks++;
            }
        }
        fmt::print("Reordered {} blocks to match the order of the mapping matrix.\n"
                  , num_blocks);
    }


    void Reorder::_find_blocks_orders(bool amorphous) {
        for (size_t i = 0; i < _blocks.block_atoms.size(); i++) {
            for (auto ba : _blocks.block_atoms[i]) {
                std::vector<std::string> block(ba.rows());
                // Hydrogens are queued after the 9 carbons
                if (amorphous) {
                    int nh = 9;
                    for (int j = 0; j < 9; j++) {
                        if (_data.num_bonded_hydrogens(ba(j)) == 1) block[j] = "C1";
                        if (_data.num_bonded_hydrogens(ba(j)) == 2) block[j] = "C2";
                        if (_data.num_bonded_hydrogens(ba(j)) == 3) block[j] = "C3";
                        for (auto h : _data.bond_table[ba(j)]) {
                            if (_data.element(h) == "H") {
                                block[nh] = _get_hydrogen_type(ba(j), h);
                                nh++;
                            }
                        }
                    }
                }
                else {
                    std::vector<std::string> C_type = {"C3", "C1", "C2"};
                    std::vector<std::string> H_type = {"HM1", "HM2", "HM3",
                                                       "HC", "HG1", "HG2"};
                    int nh = 9;
                    for (int i = 0; i < 9; i++) {
                        int ct = i % C_type.size();
                        block[i] = C_type[ct];
                    }
                    for (int j = 0; j < 18; j++) {
                        int ht = j % H_type.size();
                        block[nh] = H_type[ht];
                        nh++;
                    }
                }
                _blocks_orders[i].push_back(block);
            }
        }
    }


    std::string Reorder::_get_hydrogen_type(int C, int H) {
        Vector3d H_coords = _data.atom_coords.row(H);
        Vector3d C_coords = _data.atom_coords.row(C);
        std::vector<Vector3d> c, h;
        c.push_back(C_coords);
        if (_data.num_bonded_hydrogens(C) == 3) {
            std::vector<std::string> HM_type = {"HM1", "HM2", "HM3"};
            std::string HM = "HM";
            for (auto i : _data.bond_table[C]) {
                if (_data.element(i) != "H") {
                    c.push_back(_data.atom_coords.row(i));
                    for (auto j : _data.bond_table[i]) {
                        if (_data.element(j) != "H" && j != C) {
                            // Two C2 atoms bonded to the C1 atom
                            c.push_back(_data.atom_coords.row(j));
                        }
                    }
                }
                else { h.push_back(_data.atom_coords.row(i)); }
            }
            auto order = _hydrogens_order(c, h);
            // Order of input H atom is equal the order of similar h
            for (size_t i = 0; i < h.size(); i++) {
                if (H_coords == h[i]) HM = HM_type[order[i]];
            }
            return HM;
        }
        else if (_data.num_bonded_hydrogens(C) == 2) {
            std::vector<std::string> HG_type = {"HG1", "HG2"};
            std::string HG = "HG";
            for (auto i : _data.bond_table[C]) {
                if (_data.element(i) != "H") {
                    c.push_back(_data.atom_coords.row(i));
                }
                else h.push_back(_data.atom_coords.row(i));
            }
            auto order = _hydrogens_order(c, h);
            for (size_t i = 0; i < h.size(); i++) {
                if (H_coords == h[i]) HG = HG_type[order[i]];
            }
            return HG;
        }
        else { return "HC"; }
    }


    std::vector<int> Reorder::_hydrogens_order(Coords c, Coords h) {
        std::vector<int> order;
        Coords rh;
        if (h.size() == 3) {
            order = {-1, -1, -1};
            rh = _methyl_hydrogens(c, h);
        }
        else {
            order = {-1, -1};
            rh = _gemini_hydrogens(c, h);
        }
        // Order of h atoms are determined w.r.t the order of rh/mapping matrix.
        // If the j-th C-rh bond vector has the maximum overlap with
        // i-th existing C-h bond, then the order of bond i will be j.
        for (size_t i = 0; i < h.size(); i++) {
            double err = 1.0;
            Vector3d bh2 = _wrap_vector(c[0] - h[i]);
            for (size_t j = 0; j < rh.size(); j++) {
                Vector3d bh1 = _wrap_vector(c[0] - rh[j]);
                double e = 1.0 - abs((bh1 / bh1.norm()).dot(bh2 / bh2.norm()));
                if(e < err) {
                    err = e;
                    order[i] = j;
                }
            }
        }
        return order;
    }


    Coords Reorder::_methyl_hydrogens(Coords ccc, Coords hhh) {
        // Bond from C_i-1 to C_i+1 (Fig 1a - Theodorou(1998) MD of aPP melts).
        Vector3d rcc = _wrap_vector(ccc[2] - ccc[3]);
        if (ccc.size() == 3) {
            // For the end C3 carbons with no more available skeleton carbons
            rcc << 0, 0, 1;
        }
        Vector3d bv = -_wrap_vector(ccc[0] - ccc[1]);
        // u,v,w make a local coordinate about the methyl carbon.
        Vector3d u = bv / bv.norm();
        Vector3d v = u.cross(rcc) / (u.cross(rcc)).norm();
        Vector3d w = u.cross(v);
        double phi = ((109.0 - 90) / 180.0) * M_PI;
        std::vector<double> theta = {0.0, (2.0 / 3.0) * M_PI, (4.0 / 3.0) * M_PI};
        Coords rh;
        for (auto t : theta) {
            t += M_PI / 6;
            Vector3d r = sin(phi) * u + cos(phi) * (cos(t) * v + sin(t) * w);
            // Order of rh is similar to the mapping matrix(e.g., rh[0] is MH1)
            rh.push_back(ccc[0] + _lh * r);
        }
        return rh;
    }


    Coords Reorder:: _gemini_hydrogens(Coords cc, Coords hh) {
        // Bond from C_i to C_i-1 (Fig 1b - Theodorou(1998) MD of aPP melts).
        Vector3d bi = _wrap_vector(cc[0] - cc[1]);
        bi /= bi.norm();
        Vector3d bj = -_wrap_vector(cc[0] - cc[2]);
        bj /= bj.norm();
        Vector3d u = (bi - bj) / (sqrt(2.0 * (1.0 - bi.dot(bj))));
        Vector3d v = (bi.cross(bj)) / ((bi.cross(bj)).norm());
        double theta_h = 1.28;
        std::vector<Vector3d> rh(2);
        rh[0] = cc[0] + _lh * (sin(theta_h / 2) * u + cos(theta_h / 2) *v);
        rh[1] = cc[0] + _lh * (sin(theta_h / 2) * u - cos(theta_h / 2) *v);

        return rh;
    }


    Vector3d Reorder::_wrap_vector(Vector3d bv) {
        Eigen::Array3d s = _data.box.transpose().inverse() * bv;
        Vector3d e = s.round().matrix() * _data.box.transpose();
        return bv - e;
    }
}

