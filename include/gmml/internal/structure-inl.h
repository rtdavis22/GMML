// Copyright (c) 2012 The University of Georgia
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.

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
