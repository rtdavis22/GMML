// Author: Robert Davis

#ifndef GMML_STRUCTURE_INL_H_
#define GMML_STRUCTURE_INL_H_

// To do: move more stuff in here

namespace gmml {

inline size_t Structure::get_residue_index(size_t atom_index) const {
    for (size_t i = 1; i < residues_->size(); i++) {
        if (atom_index < residues_->at(i)->start_index)
            return i - 1;
    }
    return residues_->size() - 1;
}

}  // namespace gmml

#endif  // GMML_STRUCTURE_INL_H
