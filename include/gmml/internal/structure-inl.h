// Author: Robert Davis

#ifndef GMML_STRUCTURE_INL_H_
#define GMML_STRUCTURE_INL_H_

// TODO: move more stuff in here

namespace gmml {

inline void Structure::set_head(int atom_index) {
    if (atom_index >= size()) {
        throw std::invalid_argument("Invalid atom index.");
    }
    head_ = atom_index;
}

inline void Structure::set_tail(int atom_index) {
    if (atom_index >= size()) {
        throw std::invalid_argument("Invalid atom index.");
    }
    tail_ = atom_index;
}

inline void Structure::set_head(int residue_index,
                                const std::string& atom_name) {
    if (residue_index >= residue_count()) {
        throw std::invalid_argument("Invalid residue index.");
    }
    int atom_index = get_atom_index(residue_index, atom_name);
    if (atom_index == -1) {
        throw std::invalid_argument("Invalid atom name.");
    }
    head_ = atom_index;
}

inline void Structure::set_tail(int residue_index,
                                const std::string& atom_name) {
    if (residue_index >= residue_count()) {
        throw std::invalid_argument("Invalid residue index.");
    }
    int atom_index = get_atom_index(residue_index, atom_name);
    if (atom_index == -1) {
        throw std::invalid_argument("Invalid atom name.");
    }
    tail_ = atom_index;
}

inline std::vector<int> *Structure::get_adjacent_residues_by_atom(
        int residue_index,
        const std::string& atom_name) const {
    return get_adjacent_residues_by_atom(get_atom_index(residue_index,
                                                        atom_name));
}

}  // namespace gmml

#endif  // GMML_STRUCTURE_INL_H
