#include "cg_data.h"
#include <fmt/format.h>
#include <fmt/ostream.h>
#include <tuple>
#include <iostream>
#include <fstream>


// Mapping scheme method
namespace MappingScheme {

    // Class to store specifications of the coarse grained system
    CG_data::CG_data(const Inputs &i, LMP_data &data): _inputs(i), _data(data) {
        fmt::print("\n################## Coarse Graining ######################\n");
        _make_beads();
        _identify_beads_bonds();
        _identify_beads_angles();
        _identify_beads_torsions();
        _identify_beads_impropers();

    }


    void CG_data::_make_beads() {
        _visited = std::vector<bool> (_data.atom_types.size(), false);
        // Find the beads in each system for each defined bead type
        for (auto b : _inputs.bead_types) {
            _find_beads_of_type(b.first, b.second);
        }
        if (_inputs.symmetrical_beads.size() != 0) {
            for (size_t i=0; i<_inputs.symmetrical_beads.size(); i++) {
                fmt::print("Similar Bead Types: \t");
                for (auto t : _inputs.symmetrical_beads[i]) {
                    fmt::print("{}\t", t + 1);
                }
                fmt::print("\n");
            }
        }
    }


    void CG_data::_find_beads_of_type(char bead_name, BeadType bead) {
        fmt::print("Beads of type '{}': ", bead_name);
        // Go through all the atoms and find the ones that matches the bead definition
        int bead_type_count = 0;
        for (size_t i=0; i < _data.chains.size(); i++) {
            for (size_t j=0; j < _data.chains[i].size(); j++) {
                std::vector<int> potential_bead(bead.atom_types.size(), -1);
                std::vector<int> potential_bead_types(bead.atom_types.size(), -1);
                // Continue if already visited this atom or the type of this atom
                // is not equal to the first one in the bead
                int a = _data.chains[i][j];
                if (_visited[a] || (bead.atom_types[0] != _data.atom_types[a])) {
                    continue;
                }
                potential_bead[0] = a;
                potential_bead_types[0] = _data.atom_types[a];
                size_t k = 1;
                for (size_t l=0; l < _data.bond_table[a].size(); l++) {
                    if (k >= bead.atom_types.size()) {
                        break;
                    }
                    int b = _data.bond_table[a][l];
                    if (!_visited[b] && bead.atom_types[k] == _data.atom_types[b]) {
                        potential_bead[k] = b;
                        potential_bead_types[k] = _data.atom_types[b];
                        a = b;
                        _visited[b] = true;
                        k++;
                    }
                }
                if (potential_bead_types == bead.atom_types) {
                    Bead bead_;
                    bead_.id = beads.size();
                    bead_.chain = i;
                    bead_.type = get_bead_type(bead.type);
                    // Check it does not already exist
                    if (!(std::find(bead_types.begin(), bead_types.end()
                                         , bead_.type) != bead_types.end())) {
                        bead_types.push_back(bead_.type);
                    }
                    bead_.name = bead_name;
                    bead_.atoms = potential_bead;
                    for (auto b : potential_bead) {
                        bead_.charge += _data.atom_charges[b];
                    }
                    bead_.atom_weights = bead.weights;
                    bead_.mass = _compute_bead_mass(bead_.atoms, bead_.atom_weights);
                    bead_types_masses[bead_.type] = bead_.mass;
                    bead_.coords = _find_beads_coords(potential_bead, bead.weights);
                    beads.push_back(bead_);
                    bead_type_count++;
                }
            }
        }
        fmt::print("Found {} and their coordinates.\n", bead_type_count);
    }


    int CG_data::get_bead_type(int type) const {
        // Get the bead type by checking for the possible symmetry
        if (_inputs.symmetrical_beads.size() != 0) {
            for (size_t i=0; i<_inputs.symmetrical_beads.size(); i++) {
                for (auto t : _inputs.symmetrical_beads[i])
                    if (type == t) {
                        return _inputs.symmetrical_beads[i][0];
                    }
            }
        }
        return type;
    }


    double CG_data::_compute_bead_mass(std::vector<int> atoms, std::vector<double> weights) {
        // Find the CG bead mass to conserve the momentum of the target system following
        // 2008 - Noid - The multiscale coarse-graining method. I. A rigorous
        // bridge between atomistic and coarse-grained models
        // Go through the definition of elements and find the mass
        double C_mass = 0.0, H_mass = 0.0;
        int H_type = 0;
        for (auto const& [element, type] : _inputs.atom_types) {
            if (element[0] == 'c' || element[0] == 'C') {
                C_mass = _data.get_atom_mass(type-1);
            }
            if (element[0] == 'h' || element[0] == 'H') {
                H_mass = _data.get_atom_mass(type-1);
                H_type = type;
            }
        }
        // Find the mass of each bead
        double M = 0.0;
        for (size_t i=0; i<atoms.size(); i++) {
            double c = weights[i];
            // Find mass of the atom inside the bead
            double m;
            if (_data.atom_types[atoms[i]] != H_type) {
                m = C_mass;
            }
            else {
                m = H_mass;
            }
            M += pow(c, 2) / m;
        }
        return 1 / M;
    }


    Vector3d CG_data::_find_beads_coords(std::vector<int> atoms
                                       , std::vector<double> weights) {
        // Map atoms into the periodic cell w.r.t first atom of the bead
        Vector3d x0 = _data.atom_coords.row(atoms[0]);
        Vector3d coords = x0 * weights[0];
        for (size_t j = 1; j < atoms.size(); j++) {
            Vector3d x = _data.atom_coords.row(atoms[j]);
            Vector3d dx = x - x0;
            Vector3d ds = (_data.box.inverse() * dx).array().round().matrix();
            dx -= _data.box * ds;
            coords += (x0 + dx) * weights[j];
        }
        return coords;
    }


    // Two beads are bonded if any atoms in one is bonded to any atom in the other
    void CG_data::_identify_beads_bonds() {
        bead_bonds.assign(beads.size(), {});
        for (auto b1 : beads) {
            for (auto b2 : beads) {
                if ((b1.chain == b2.chain) && (b1.id != b2.id)) {
                    if (_are_bonded(b1.atoms, b2.atoms)) {
                        bead_bonds[b1.id].push_back(b2.id);
                        _bond_types(b1.type, b2.type);
                    }
                }
            }
        }
        int num_bonds = 0;
        for (size_t i=0; i<bead_bonds.size(); i++) {
            for (size_t j=0; j<bead_bonds[i].size(); j++) {
                if ((int)i < bead_bonds[i][j]) {
                    num_bonds++;
                }
            }
        }
        fmt::print("Found {} bonds and {} bond types.\n"
                  , num_bonds, bead_bond_types.size());
    }


    bool CG_data::_are_bonded(std::vector<int> b1_atoms, std::vector<int> b2_atoms) {
        for (auto atom1 : b1_atoms) {
            for (auto bonded_to_atom1 : _data.bond_table[atom1]) {
                for (auto atom2 : b2_atoms) {
                    if (bonded_to_atom1 == atom2) {
                        return true;
                    }
                }
            }
        }
        return false;
    }


    void CG_data::_bond_types(int i, int j) {
        std::tuple<int, int> type;
        if (i < j) {
            type = std::make_tuple(i, j);
        }
        else {
            type = std::make_tuple(j, i);
        }
        // Check it does not already exist
        if (!(std::find(bead_bond_types.begin(), bead_bond_types.end(), type) !=
                                                          bead_bond_types.end())) {
            bead_bond_types.push_back(type);
        }
    }


    int CG_data::get_bond_type(int i, int j) const {
        size_t bbt;
        std::tuple<int, int> type;
        if (beads[i].type < beads[j].type) {
            type = std::make_tuple(beads[i].type, beads[j].type);
        }
        else {
            type = std::make_tuple(beads[j].type, beads[i].type);
        }
        for (bbt=0; bbt<bead_bond_types.size(); bbt++) {
            if (type == bead_bond_types[bbt]) {
                break;
            }
        }
        return bbt;
    }


    // Three linearly bonded beads (no loop) make an angle
    void CG_data::_identify_beads_angles() {
        for (auto b : beads) {
            for (auto i : bead_bonds[b.id]) {
                for (auto j : bead_bonds[b.id]) {
                    if (i < j) {
                        bead_angles.push_back(std::make_tuple(i, b.id, j));
                        _angle_types(i, b.id, j);
                    }
                }
            }
        }
        fmt::print("Found {} bond angles and {} angle types.\n"
                  , bead_angles.size(), bead_angle_types.size());
    }


    void CG_data::_angle_types(int i, int j, int k) {
        std::tuple<int, int, int> type;
        if (beads[i].type < beads[k].type) {
            type = std::make_tuple(beads[i].type, beads[j].type, beads[k].type);
        }
        else if (beads[i].type > beads[k].type) {
            type = std::make_tuple(beads[k].type, beads[j].type, beads[i].type);
        }
        // Check it does not already exist
        if (!(std::find(bead_angle_types.begin(), bead_angle_types.end(), type) !=
                                                          bead_angle_types.end())) {
            bead_angle_types.push_back(type);
        }
    }


    int CG_data::get_angle_type(int i, int j, int k) const {
        int at = 0;
        for (auto bat : bead_angle_types) {
            at++;
            if (std::make_tuple(beads[i].type, beads[j].type, beads[k].type) == bat) {
                break;
            }
        }
        return at;
    }


    // Four linearly bonded beads (no loop) make a torsion/dihedral angle
    void CG_data::_identify_beads_torsions() {
        for (int j=0; j<(int)bead_bonds.size(); j++) {
            for (auto i : bead_bonds[j]) {
                for (auto k : bead_bonds[j]) {
                    if ((k > j) && (k != i)) {
                        for (auto l : bead_bonds[k]) {
                            if ((l != j) && (l != i)){
                                bead_torsions.push_back(std::make_tuple(i, j, k, l));
                                _torsion_types(i, j, k, l);
                            }
                        }
                    }
                }
            }
        }
        fmt::print("Found {} torsion angles and {} torsion angle types.\n"
                  , bead_torsions.size(), bead_torsion_types.size());
    }


    void CG_data::_torsion_types(int i, int j, int k, int l) {
        std::tuple<int, int, int, int> type;
        if (beads[i].type > beads[l].type) {
            type = std::make_tuple(beads[l].type, beads[k].type,
                                   beads[j].type, beads[i].type);
        }
        else if ((beads[i].type == beads[l].type) && (beads[j].type > beads[k].type)) {
            type = std::make_tuple(beads[l].type, beads[k].type,
                                   beads[j].type, beads[i].type);
        }
        else {
            type = std::make_tuple(beads[i].type, beads[j].type,
                                   beads[k].type, beads[l].type);
        }
        // Check it does not already exist
        if (!(std::find(bead_torsion_types.begin(), bead_torsion_types.end(), type)
                                                    != bead_torsion_types.end())) {
            bead_torsion_types.push_back(type);
        }
    }


    int CG_data::get_torsion_type(int i, int j, int k, int l) const {
        int tt = 0;
        for (auto btt : bead_torsion_types) {
            tt++;
            if (std::make_tuple(beads[i].type, beads[j].type
                              , beads[k].type, beads[l].type) == btt) {
                break;
            }
        }
        return tt;
    }


    // Three beads connected to a bead with three of them in a plane and one in another
    void CG_data::_identify_beads_impropers() {
        for (int j=0; j<(int)bead_bonds.size(); j++) {
            // Only beads with three bonds or more can form improper angle
            if (bead_bonds[j].size() > 2) {
                for (auto i : bead_bonds[j]) {
                    for (auto k : bead_bonds[j]) {
                        for (auto l : bead_bonds[j]) {
                            if (i < j && i < l && l < k && j < k){
                                bead_impropers.push_back(std::make_tuple(i, j, k, l));
                                _improper_types(i, j, k, l);
                            }
                        }
                    }
                }
            }
        }
        fmt::print("Found {} improper angles and {} improper angle types.\n"
                  , bead_impropers.size(), bead_improper_types.size());
    }


    void CG_data::_improper_types(int i, int j, int k, int l) {
        std::tuple<int, int, int, int> type;
        // TODO: Needs to be double checked
        if (beads[i].type == beads[j].type && beads[j].type == beads[l].type &&
            beads[k].type > beads[l].type) {
            type = std::make_tuple(beads[i].type, beads[j].type,
                                   beads[k].type, beads[l].type);
            // Check it does not already exist
            if (!(std::find(bead_improper_types.begin(), bead_improper_types.end(), type)
                                                    != bead_improper_types.end())) {
                bead_improper_types.push_back(type);
            }
        }
    }


    int CG_data::get_improper_type(int i, int j, int k, int l) const {
        int it = 0;
        for (auto bit : bead_improper_types) {
            it++;
            if (std::make_tuple(beads[i].type, beads[j].type
                              , beads[k].type, beads[l].type) == bit) {
                break;
            }
        }
        return it;
    }
}


// Mapping matrix method
namespace MappingMatrix {

    // Class to store specifications of the coarse grained system
    CG_data::CG_data(const Inputs &i, const LMP_data &d, const Mapping &map
                   , const Reorder b) :  _inputs(i), _data(d), _mapping(map)
                   , _blocks(b)  {
        fmt::print("\n################## Coarse Graining ######################\n");
        _make_beads();
        _find_beads_coords();
        _identify_beads_bonds();
        _identify_beads_angles();
        _identify_beads_torsions();
        _identify_beads_impropers();
    }


    void CG_data::_make_beads() {
        // Atoms with nonzero weights in each col. of mapping matrix are in the bead
        int bead_id = 0;
        for (size_t k=0; k < _blocks.block_atoms.size(); k++) {
            for (auto block : _blocks.block_atoms[k]) {
                for (int j=0; j<_mapping.map_mat.cols(); j++) {
                    Bead bead;
                    bead.id = bead_id;
                    bead.type = get_bead_type(j);
                    // Check it does not already exist
                    if (!(std::find(bead_types.begin(), bead_types.end()
                                        , bead.type) != bead_types.end())) {
                            bead_types.push_back(bead.type);
                    }
                    bead.chain = k;
                    for (int i=0; i<_mapping.map_mat.rows(); i++) {
                        if (_mapping.bead_map(i, j) != 0.0) {
                            bead.atoms.push_back(block(i));
                        }
                    }
                    bead.mass = _compute_bead_mass(j);
                    bead.charge = _compute_bead_charge(j);
                    bead_types_masses[bead.type] = bead.mass;
                    beads.push_back(bead);
                    bead_id++;
               }
            }
        }
        fmt::print("Found {} beads and {} bead types in {} blocks.\n"
                   , beads.size(), bead_types.size(), _blocks.num_blocks);
    }


    int CG_data::get_bead_type(int type) const {
        // Get the bead type by checking for the possible symmetry
        if (_inputs.symmetrical_beads.size() != 0) {
            for (size_t i=0; i<_inputs.symmetrical_beads.size(); i++) {
                for (auto t : _inputs.symmetrical_beads[i])
                    if (type == t) {
                        return _inputs.symmetrical_beads[i][0];
                    }
            }
        }
        return type;
    }


    double CG_data::_compute_bead_mass(int j) {
        // Find the CG bead mass to conserve the momentum of the target system following
        // 2008 - Noid - The multiscale coarse-graining method. I. A rigorous b
        // ridge between atomistic and coarse-grained models
        // Go through the definition of elements and find the mass
        double C_mass = 0.0, H_mass = 0.0;
        for (auto const& [type, element] : _inputs.elements) {
            if (element == "C") {
                C_mass = _data.get_atom_mass(type-1);
            }
            if (element == "H") {
                H_mass = _data.get_atom_mass(type-1);
            }
        }
        // Find the mass of each bead
        double M = 0.0;
        for (int i=0; i<_mapping.map_mat.rows(); i++) {
            double c = _mapping.map_mat(i, j);
            // Find mass of the atom inside the bead
            double m;
            if (_mapping.map_order[i] == "C1" ||
                _mapping.map_order[i] == "C2" ||
                _mapping.map_order[i] == "C3") {
                m = C_mass;
            }
            else {
                m = H_mass;
            }
            M += pow(c, 2) / m;
        }
        return 1 / M;
    }


    double CG_data::_compute_bead_charge(int j) {
        // Sum all atom charges in the bead with or without any weight in
        double charge = 0.0;
        for (int i=0; i<_mapping.map_mat.rows(); i++) {
            auto k = _blocks.block_atoms[0][0];
            if (_mapping.map_mat(i, j) != 0.0) {
                charge += _data.atom_charges[k(i)];
                for (auto l : _mapping.connections.at(i)) {
                    if (_mapping.map_mat(l, j) == 0.0) {
                        charge += _data.atom_charges[k(l)];
                    }
                }
            }
        }
        // LAMMPS consider charges with 4 decimals
        if (charge < 1e-4) charge = 0.0;
        return charge;
    }


    void CG_data::_find_beads_coords() {
        beads_coords = BlockAtoms (_blocks.num_blocks);
        if (_inputs.amorphous) {
            for (size_t i=0; i<_blocks.block_atoms.size(); i++) {
                int j = 0;
                for (auto block : _blocks.block_atoms[i]) {
                    // Map atoms into the periodic cell w.r.t first atom of the chain
                    auto x0 = _data.atom_coords.row(block(0));
                    for (int k=1; k<block.rows(); k++) {
                        MatrixXd x = _data.atom_coords.row(block(k));
                        for (int l=0; l<3; l++) {
                            if (x(l) - x0(l) > 0.5 * _data.box(l, l)) {
                                x(l) -= _data.box(l, l);
                            }
                            if (x(l) - x0(l) < -0.5 * _data.box(l, l)) {
                                x(l) += _data.box(l, l);
                            }
                        }
                        _blocks.block_coords[i][j].row(k) = x;
                    }
                    MatrixXd bead_coords = _mapping.map_mat.transpose()
                                         * _blocks.block_coords[i][j];
                    // Adjust the bead coordinates to be inside the periodic cell
                    for (int k=0; k<bead_coords.rows(); k++) {
                        for (int l=0; l<3; l++) {
                            if (bead_coords(k, l) < 0.0) {
                                bead_coords(k, l) += _data.box(l, l);
                            }
                            if (bead_coords(k, l) >= _data.box(l, l)) {
                                bead_coords(k, l) -= _data.box(l, l);
                            }
                        }
                    }
                    beads_coords[i].push_back(bead_coords);
                    j++;
                }
            }
        }
        else {
            // Alpha-ipp crystal phase has a triclinic box and is not orthogonal
            for (size_t i=0; i<_blocks.block_coords.size(); i++) {
                for (auto block : _blocks.block_coords[i]) {
                    Vector3d x0 = block.row(0);
                    for (int j=1; j<block.rows(); j++) {
                        Vector3d x = block.row(j);
                        Vector3d dx = x - x0;
                        Vector3d ds = (_data.box.inverse() * dx).array().round().matrix();
                        dx -= _data.box * ds;
                        x = x0 + dx;
                        block.row(j) = x;
                    }
                    beads_coords[i].push_back(_mapping.map_mat.transpose() * block);
                }
            }
        }
        fmt::print("Found coordinates for {} beads.\n"
                  , beads.size(), _blocks.num_blocks);
    }


    // Two beads are bonded if any atoms in one is bonded to any atom in the other
    void CG_data::_identify_beads_bonds() {
        bead_bonds.assign(beads.size(), {});
        for (auto b1 : beads) {
            for (auto b2 : beads) {
                if ((b1.chain == b2.chain) && (b1.id != b2.id)) {
                    if (_are_bonded(b1.atoms, b2.atoms)) {
                        bead_bonds[b1.id].push_back(b2.id);
                        _bond_types(b1.type, b2.type);
                    }
                }
            }
        }
        int num_bonds = 0;
        for (size_t i=0; i<bead_bonds.size(); i++) {
            for (size_t j=0; j<bead_bonds[i].size(); j++) {
                if ((int)i < bead_bonds[i][j]) {
                    num_bonds++;
                }
            }
        }
        fmt::print("Found {} bonds and {} bond types.\n"
                  , num_bonds, bead_bond_types.size());
    }


    bool CG_data::_are_bonded(std::vector<int> b1_atoms, std::vector<int> b2_atoms) {
        for (auto atom1 : b1_atoms) {
            if (_data.element(atom1) == "H") continue;
            for (auto bonded_to_atom1 : _data.bond_table[atom1]) {
                if (_data.element(bonded_to_atom1) == "H") continue;
                for (auto atom2 : b2_atoms) {
                    if (_data.element(atom2) == "H") continue;
                    if (bonded_to_atom1 == atom2) {
                        return true;
                    }
                }
            }
        }
        return false;
    }


    void CG_data::_bond_types(int i, int j) {
        std::tuple<int, int> type;
        if (i < j) {
            type = std::make_tuple(i, j);
        }
        else {
            type = std::make_tuple(j, i);
        }
        // Check it does not already exist
        if (!(std::find(bead_bond_types.begin(), bead_bond_types.end(), type) !=
                                                          bead_bond_types.end())) {
            bead_bond_types.push_back(type);
        }
    }


    int CG_data::get_bond_type(int i, int j) const {
        size_t bbt;
        std::tuple<int, int> type;
        if (beads[i].type < beads[j].type) {
            type = std::make_tuple(beads[i].type, beads[j].type);
        }
        else {
            type = std::make_tuple(beads[j].type, beads[i].type);
        }
        for (bbt=0; bbt<bead_bond_types.size(); bbt++) {
            if (type == bead_bond_types[bbt]) {
                break;
            }
        }
        return bbt;
    }


    // Three linearly bonded beads (no loop) make an angle
    void CG_data::_identify_beads_angles() {
        for (auto b : beads) {
            for (auto i : bead_bonds[b.id]) {
                for (auto j : bead_bonds[b.id]) {
                    if (i < j) {
                        bead_angles.push_back(std::make_tuple(i, b.id, j));
                        _angle_types(i, b.id, j);
                    }
                }
            }
        }
        fmt::print("Found {} bond angles and {} angle types.\n"
                  , bead_angles.size(), bead_angle_types.size());
    }


    void CG_data::_angle_types(int i, int j, int k) {
        std::tuple<int, int, int> type;
        if (beads[i].type < beads[k].type) {
            type = std::make_tuple(beads[i].type, beads[j].type, beads[k].type);
        }
        else if (beads[i].type > beads[k].type) {
            type = std::make_tuple(beads[k].type, beads[j].type, beads[i].type);
        }
        // Check it does not already exist
        if (!(std::find(bead_angle_types.begin(), bead_angle_types.end(), type) !=
                                                          bead_angle_types.end())) {
            bead_angle_types.push_back(type);
        }
    }


    int CG_data::get_angle_type(int i, int j, int k) const {
        int at = 0;
        for (auto bat : bead_angle_types) {
            at++;
            if (std::make_tuple(beads[i].type, beads[j].type, beads[k].type) == bat) {
                break;
            }
        }
        return at;
    }


    // Four linearly bonded beads (no loop) make a torsion/dihedral angle
    void CG_data::_identify_beads_torsions() {
        for (int j=0; j<(int)bead_bonds.size(); j++) {
            for (auto i : bead_bonds[j]) {
                for (auto k : bead_bonds[j]) {
                    if ((k > j) && (k != i)) {
                        for (auto l : bead_bonds[k]) {
                            if ((l != j) && (l != i)){
                                bead_torsions.push_back(std::make_tuple(i, j, k, l));
                                _torsion_types(i, j, k, l);
                            }
                        }
                    }
                }
            }
        }
        fmt::print("Found {} torsion angles and {} torsion angle types.\n"
                  , bead_torsions.size(), bead_torsion_types.size());
    }


    void CG_data::_torsion_types(int i, int j, int k, int l) {
        std::tuple<int, int, int, int> type;
        if (beads[i].type > beads[l].type) {
            type = std::make_tuple(beads[l].type, beads[k].type,
                                   beads[j].type, beads[i].type);
        }
        else if ((beads[i].type == beads[l].type) && (beads[j].type > beads[k].type)) {
            type = std::make_tuple(beads[l].type, beads[k].type,
                                   beads[j].type, beads[i].type);
        }
        else {
            type = std::make_tuple(beads[i].type, beads[j].type,
                                   beads[k].type, beads[l].type);
        }
        // Check it does not already exist
        if (!(std::find(bead_torsion_types.begin(), bead_torsion_types.end(), type)
                                                    != bead_torsion_types.end())) {
            bead_torsion_types.push_back(type);
        }
    }


    int CG_data::get_torsion_type(int i, int j, int k, int l) const {
        int tt = 0;
        for (auto btt : bead_torsion_types) {
            tt++;
            if (std::make_tuple(beads[i].type, beads[j].type
                              , beads[k].type, beads[l].type) == btt) {
                break;
            }
        }
        return tt;
    }


    // Three beads connected to a bead with three of them in a plane and one in another
    void CG_data::_identify_beads_impropers() {
        for (int j=0; j<(int)bead_bonds.size(); j++) {
            // Only beads with three bonds or more can form improper angle
            if (bead_bonds[j].size() > 2) {
                for (auto i : bead_bonds[j]) {
                    for (auto k : bead_bonds[j]) {
                        for (auto l : bead_bonds[j]) {
                            if (i < j && i < l && l < k && j < k){
                                bead_impropers.push_back(std::make_tuple(i, j, k, l));
                                _improper_types(i, j, k, l);
                            }
                        }
                    }
                }
            }
        }
        fmt::print("Found {} improper angles and {} improper angle types.\n"
                  , bead_impropers.size(), bead_improper_types.size());
    }


    void CG_data::_improper_types(int i, int j, int k, int l) {
        std::tuple<int, int, int, int> type;
        // TODO: Needs to be double checked
        if (beads[i].type == beads[j].type && beads[j].type == beads[k].type &&
            beads[j].type > beads[l].type) {
            type = std::make_tuple(beads[i].type, beads[j].type,
                                   beads[k].type, beads[l].type);
            // Check it does not already exist
            if (!(std::find(bead_improper_types.begin(), bead_improper_types.end(), type)
                                                    != bead_improper_types.end())) {
                bead_improper_types.push_back(type);
            }
        }
    }


    int CG_data::get_improper_type(int i, int j, int k, int l) const {
        int it = 0;
        for (auto bit : bead_improper_types) {
            it++;
            if (std::make_tuple(beads[i].type, beads[j].type
                              , beads[k].type, beads[l].type) == bit) {
                break;
            }
        }
        return it;
    }
}

