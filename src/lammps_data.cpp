#include "lammps_data.h"
#include "string_tools.h"
#include <fmt/format.h>
#include <iostream>


LMP_data::LMP_data(const Inputs &inputs) : _inputs(inputs) {
    fmt::print("\n#################### Initializing #######################\n");
    std::fstream fid;
    fid.open(_inputs.lmp_data, std::fstream::in);
    if (!fid) {
        std::cout << "Data file does not exist.\n";
        exit(0);
    }
    // Triclinic box (covers the cubical box as well)
    while (fid) {
        auto args = split(before(read_line(fid), "#"));
        if (args.empty()) {
            continue;
        }
        else if (args.size() == 2 && args.back() == "atoms") {
            int n = from_string<int>(args[0]);
            atom_coords = MatrixXd::Zero(n, 3);
            atom_mols.assign(n, 0);
            atom_types.assign(n, 0);
            atom_charges.assign(n, 0.0);
            bond_table.assign(n, {});
        }
        else if (args.size() == 4 && args[3] == "xhi") {
            box(0,0) = str2dbl(args[1]) - str2dbl(args[0]);
        }
        else if (args.size() == 4 && args[3] == "yhi") {
            box(1,1) = str2dbl(args[1]) - str2dbl(args[0]);
        }
        else if (args.size() == 4 && args[3] == "zhi") {
            box(2,2) = str2dbl(args[1]) - str2dbl(args[0]);
        }
        else if (args.size() == 6 && args[5] == "yz") {
            box(0,1) = str2dbl(args[0]);
            box(0,2) = str2dbl(args[1]);
            box(1,2) = str2dbl(args[2]);
        }
        else if (args.size() == 1 && args[0] == "Masses") {
            while (fid) {
                auto pos = fid.tellp();
                args = split(before(read_line(fid), "#"));
                if (args.size() == 1) {
                    fid.seekp(pos);
                    break;
                }
                else if (args.size() == 2) {
                    auto t = from_string<int>(args[0]) - 1;
                    atom_masses[t] = from_string<double>(args[1]);
                }
            }
        }
        else if (args.size() == 1 && args[0] == "Atoms") {
            while (fid) {
                auto pos = fid.tellp();
                args = split(before(read_line(fid), "#"));
                if (args.size() == 1) {
                    fid.seekp(pos);
                    break;
                }
                if (args.size() > 6) {
                    auto i = from_string<int>(args[0]) - 1;
                    atom_mols[i] = from_string<int>(args[1]);
                    atom_types[i] = from_string<int>(args[2]);
                    atom_charges[i] = from_string<double>(args[3]);
                    atom_coords(i,0) = from_string<double>(args[4]);
                    atom_coords(i,1) = from_string<double>(args[5]);
                    atom_coords(i,2) = from_string<double>(args[6]);
                }
            }
        }
        else if (args.size() == 1 && args[0] == "Bonds") {
            while (fid) {
                auto pos = fid.tellp();
                args = split(before(read_line(fid), "#"));
                if (args.size() == 1) {
                    fid.seekp(pos);
                    break;
                }
                if (args.size() == 4) {
                    int i = from_string<int>(args[2]) - 1;
                    int j = from_string<int>(args[3]) - 1;
                    bond_table[i].push_back(j);
                    bond_table[j].push_back(i);
                }
            }
        }
    }
    for (auto &b : bond_table) {
        std::sort(b.begin(), b.end());
    }
    fmt::print("Read {} atoms from {}.\n", atom_coords.rows(), _inputs.lmp_data);
    // Find the chains
    _find_chains();
}


void LMP_data::remove_bonds_across_pbc() {
    const auto max_bond_length = 2.0;
    int bonds_removed = 0;
    for (size_t i=0; i<bond_table.size(); ++i) {
        for (size_t j = 0; j < bond_table[i].size(); j++) {
            auto rij = (atom_coords.row(i) - atom_coords.row(bond_table[i][j])).norm();
            if (rij > max_bond_length) {
                removed_bonds[i] = bond_table[i][j];
                bond_table[i].erase(bond_table[i].begin() + j);
                bonds_removed += 1;
            }
        }
    }
    fmt::print("Deleted {} bonds across PBC.\n", bonds_removed);
}


void LMP_data::connect_bonds_across_pbc() {
    int reconnected_bonds = 0;
    for (auto const &[i, j] : removed_bonds) {
        bond_table[i].push_back(j);
        reconnected_bonds++;
    }
    fmt::print("Reconnected {} bonds across PBC.\n", reconnected_bonds);
}


std::string LMP_data::element(int i) const {
    for (auto e : _inputs.elements) {
        if (atom_types[i] == e.first) {
            return e.second;
        }
    }
    return {};
}


int LMP_data::num_bonded_hydrogens(int atom) const {
    int nh = 0;
    for (auto bonded : bond_table[atom]) {
        if (element(bonded) == "H") {
            nh++;
        }
    }
    return nh;
}


void LMP_data::_find_chains() {
    if (!_inputs.amorphous && !_inputs.mapping_matrix.empty()) {
        // Remove the periodic boundary conditions for crystal system
        remove_bonds_across_pbc();
        // In crystal data file each chain is considered as a molecule
        int num_chains = *std::max_element(atom_mols.begin(), atom_mols.end());
        for (int i = 1; i < num_chains + 1; i++) {
            std::vector<int> chain;
            for (size_t j = 0; j < atom_mols.size(); j++) {
                if (atom_mols[j] == i) {
                    chain.push_back(j);
                }
            }
            chains.push_back(chain);
        }
    }
    else {
        std::vector<bool> done(bond_table.size(), false);
        for (int i = 0; i < (int)bond_table.size(); i++) {
            if (done[i]) continue;
            chains.push_back({i});
            done[i] = true;
            std::vector<int> queue = bond_table[i];
            while (!queue.empty()) {
                auto j = queue.back();
                queue.pop_back();
                if (done[j]) continue;
                chains.back().push_back(j);
                done[j] = true;
                for (auto k : bond_table[j]) {
                    queue.push_back(k);
                }
            }
        }
    }
    for (auto &c : chains) {
        std::sort(c.begin(), c.end());
    }
    fmt::print("Found {} chains in this system.\n", chains.size());
}

