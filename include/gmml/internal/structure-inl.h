// Author: Robert Davis

#ifndef GMML_STRUCTURE_INL_H_
#define GMML_STRUCTURE_INL_H_

// TODO: move more stuff in here

namespace gmml {

inline std::vector<int> *Structure::get_adjacent_residues_by_atom(
        int residue_index,
        const std::string& atom_name) const {
    return get_adjacent_residues_by_atom(get_atom_index(residue_index,
                                                        atom_name));
}

}  // namespace gmml

#endif  // GMML_STRUCTURE_INL_H
