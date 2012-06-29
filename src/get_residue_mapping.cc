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

#include "gmml/internal/get_residue_mapping.h"

#include <algorithm>
#include <string>

#include "gmml/internal/atom.h"
#include "gmml/internal/structure.h"

using std::string;
using std::vector;

namespace gmml {
namespace {

class GetResidueMapping {
  public:
    GetResidueMapping(const Structure *structure1, const Structure *structure2)
            : structure1_(structure1), structure2_(structure2) {
    }

    vector<int> *operator()() {
        if (structure1_->size() != structure2_->size()) {
            return NULL;
        }

        mapping_ = new vector<int>(structure1_->residue_count(), -1);

        if (find_mapping()) {
            return mapping_;
        } else {
            delete mapping_;
            return NULL;
        }
    }

  private:
    bool find_mapping() {
        vector<vector<int> > *molecules1 = structure1_->extract_molecules();
        vector<vector<int> > *molecules2 = structure2_->extract_molecules();

        vector<int> unused_molecules;
        for (int i = 0; i < molecules2->size(); i++) {
            unused_molecules.push_back(i);
        }

        for (int i = 0; i < molecules1->size(); i++) {
            bool mapping_found = false;
            for (int j = 0; j < unused_molecules.size(); j++) {
                const vector<int>& molecule1 = (*molecules1)[i];
                const vector<int>& molecule2 =
                        (*molecules2)[unused_molecules[j]];
                if (molecule1.size() != molecule2.size()) {
                    continue;
                }
                if (map_molecules(molecule1, molecule2)) {
                    unused_molecules.erase(unused_molecules.begin() + j);
                    mapping_found = true;
                    break;
                }
            }
            if (!mapping_found) {
                return false;
            }
        }
        return unused_molecules.empty();
    }

    bool map_molecules(const vector<int>& residues1,
                       const vector<int>& residues2) {
        string residue1_name = structure1_->residues(residues1[0])->name();
        for (int i = 0; i < residues2.size(); i++) {
            if (map_residues(residues1[0], residues2[i]))
                return true;
        }
        return false;
    }

    struct ResidueLink {
        ResidueLink(const string& atom_name, int residue_number)
                : atom_name(atom_name), residue_number(residue_number) {}

        std::string atom_name;
        int residue_number;
    };

    struct ResidueLinkComparer {
        bool operator()(const ResidueLink& lhs, const ResidueLink& rhs) {
            return lhs.atom_name < rhs.atom_name;
        }
    };

    bool map_residues(int residue1, int residue2) {
        if ((*mapping_)[residue1] != -1) {
            return (*mapping_)[residue1] == residue2;
        }

        if (structure1_->residues(residue1)->name() !=
                structure2_->residues(residue2)->name()) {
            return false;
        }

        vector<ResidueLink> *links1 = get_sorted_links(structure1_, residue1);
        vector<ResidueLink> *links2 = get_sorted_links(structure2_, residue2);
        if (links1->size() != links2->size()) {
            return false;
        }

        (*mapping_)[residue1] = residue2;

        bool ret_val = true;
        for (int i = 0; i < links1->size(); i++) {
            const ResidueLink& link1 = (*links1)[i];
            const ResidueLink& link2 = (*links2)[i];
            if (link1.atom_name != link2.atom_name) {
                ret_val = false;
                break;
            }
            bool success = map_residues(link1.residue_number,
                                        link2.residue_number);
            if (!success) {
                ret_val = false;
                break;
            }
        }

        delete links1;
        delete links2;

        if (!ret_val) {
            (*mapping_)[residue1] = -1;
        }

        return ret_val;
    }

    static vector<ResidueLink> *get_sorted_links(const Structure *structure,
                                                 int index) {
        vector<ResidueLink> *links = get_links(structure, index);
        std::sort(links->begin(), links->end(), ResidueLinkComparer());
        return links;
    }

    static vector<ResidueLink> *get_links(const Structure *structure,
                                          int index) {
        vector<ResidueLink> *links = new vector<ResidueLink>;
        int atom1 = structure->get_atom_index(index, 0);
        int residue_size = structure->residues(index)->size();
        for (int i = atom1; i < atom1 + residue_size; i++) {
            const Structure::AdjList& adj_atoms = structure->bonds(i);
            for (int j = 0; j < adj_atoms.size(); j++) {
                int other_atom = adj_atoms[j];
                int other_residue_index =
                        structure->get_residue_index(other_atom);
                if (structure->get_residue_index(other_atom) != index) {
                    links->push_back(ResidueLink(structure->atoms(i)->name(),
                                                 other_residue_index));
                                                
                }
            }
        }
        return links;
    }

    vector<int> *mapping_;
    const Structure *structure1_;
    const Structure *structure2_;
};

}  // namespace

vector<int> *get_residue_mapping(const Structure *structure1,
                                 const Structure *structure2) {
    return GetResidueMapping(structure1, structure2)();
}

}  // namespace gmml
