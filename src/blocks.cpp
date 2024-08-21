#include "blocks.h"
#include <fmt/format.h>
#include <iostream>
#include <algorithm>
#include "string_tools.h"


using Eigen::MatrixXd;

typedef std::vector<std::vector<MatrixXd>> BlockAtoms;


// This class is only for the mapping matrix method
namespace MappingMatrix {

    // Constructor
    Blocks::Blocks(const Inputs &i, LMP_data &d) : _inputs(i), _data(d) {
        // SMILES pattern type order
        auto SMILES = _inputs.smiles_block;
        for (size_t i = 0; i < SMILES.size(); i++) {
            if (SMILES[i] == '(' || SMILES[i] == ')') continue;
            else if (SMILES[i] == 'C' && SMILES[i+1] == '(') {
                 _pattern_order.push_back(1);
            }
            else if (SMILES[i] == 'C' && SMILES[i+1] == ')') {
                _pattern_order.push_back(3);
            }
            else {
                _pattern_order.push_back(2);
            }
        }
        _visited.assign(_data.atom_types.size(), false);
        _find_blocks();
        // Update the bond table of the system to account for the removed extra H
        d.bond_table = _data.bond_table;
        // Reconnect the removed bonds across PBC for crystal system
        if (!_inputs.amorphous && !_inputs.finite_crystal) {
            d.connect_bonds_across_pbc();
        }
    }


    void Blocks::_find_blocks() {
        block_atoms = BlockAtoms (_data.chains.size());
        block_coords = BlockAtoms (_data.chains.size());
        for (int start = 0; start < (int)_data.bond_table.size(); start++) {
            if (_visited[start] || _data.element(start) != "C") {
                continue;
            }
            auto atoms = _find_SMILES_atoms(start);
            if (atoms.size() != _pattern_order.size()) {
                continue;
            }
            int chain_number = _block_chain_number(atoms);
            // Add the H atoms
            int n = atoms.size();
            for (int i=0; i<n; i++) {
                for (auto j : _data.bond_table[atoms[i]]) {
                    if (_data.element(j) == "H") {
                        atoms.push_back(j);
                        _visited[j] = true;
                    }
                }
            }
            MatrixXd b_atoms(27, 1);
            MatrixXd b_coords(27, 3);
            if (_make_SMILES(atoms) == _inputs.smiles_block) {
                // The following needs to be checked for every system
                if (atoms.size() > 27) {
                    // Tail atoms have an extra H that needs to be removed
                    int count = 0;
                    // The 7th element is the C with the extra H
                    for (auto b : _data.bond_table[atoms[6]]) {
                        if (_data.element(b) == "H" ) {
                            _data.bond_table[atoms[6]].erase(
                                    _data.bond_table[atoms[6]].begin() + count);
                            break;
                        }
                        count++;
                    }
                    atoms.erase(atoms.begin() + 21);
                }
                for (size_t i=0; i<atoms.size(); i++) {
                    b_atoms(i) = atoms[i];
                    b_coords.row(i) = _data.atom_coords.row(atoms[i]);
                }
                block_atoms[chain_number].push_back(b_atoms);
                block_coords[chain_number].push_back(b_coords);
                num_blocks++;
            }
            else {
                if (atoms.size() > 27) {
                    // Head atoms have an extra H that needs to be removed
                    int count = 0;
                    // H iype is 3: The 3rd element has the C with the extra H (atoms[2])
                    // H type is 1: the 1st element has the C with the extra H (atoms[0])
                    for (auto b : _data.bond_table[atoms[0]]) {
                        if (_data.element(b) == "H") {
                            _data.bond_table[atoms[0]].erase(
                                    _data.bond_table[atoms[0]].begin() + count);
                            break;
                        }
                        count++;
                    }
                    // Find the location of the extra H in the atoms
                    atoms.erase(atoms.begin() + 13);
                    if (_make_SMILES(atoms) == _inputs.smiles_block) {
                        for (size_t i=0; i<atoms.size(); i++) {
                            b_atoms(i) = atoms[i];
                            b_coords.row(i) = _data.atom_coords.row(atoms[i]);
                        }
                        block_atoms[chain_number].push_back(b_atoms);
                        block_coords[chain_number].push_back(b_coords);
                        num_blocks++;
                    }
                }
            }
        }
        fmt::print("Found {} blocks in {} chains.\n", num_blocks, _data.chains.size());
    }


    int Blocks::_block_chain_number(const std::vector<int> &block) {
        for (size_t i = 0; i < _data.chains.size(); i++) {
            for (auto j : _data.chains[i]) {
                if (std::find(block.begin(), block.end(), j) != block.end()) {
                    return i;
                }
            }
        }
        return {};
    }


    std::string Blocks::_make_SMILES(std::vector<int> atoms) {
        std::string pattern = "";
        for (size_t i=0; i<_pattern_order.size(); i++) {
            if (_data.num_bonded_hydrogens(atoms[i]) == 3) {
                pattern += "(" + _data.element(atoms[i]) + ")";
                continue;
            }
            pattern += _data.element(atoms[i]);
        }
        return pattern;
    }


    std::vector<int> Blocks::_find_SMILES_atoms(int start) {
        std::vector<int> atoms;
        atoms.push_back(start);
        _visited[start] = true;
        size_t counter = 0;
        while (counter < _pattern_order.size()) {
            for (auto i : _data.bond_table[start]) {
                auto nh = _data.num_bonded_hydrogens(i);
                if (_visited[i] || _data.element(i) == "H") continue;
                if (nh == 3) {
                    if (atoms.size() >= _pattern_order.size()) break;
                    atoms.push_back(i);
                    _visited[i] = true;
                }
                else {
                    if (atoms.size() >= _pattern_order.size()) break;
                    atoms.push_back(i);
                    _visited[i] = true;
                    start = i;
                }
            }
            counter++;
        }
        std::sort(atoms.begin(), atoms.end());
        return atoms;
    }
}

