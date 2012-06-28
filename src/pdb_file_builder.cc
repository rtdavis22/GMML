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

#include "gmml/internal/pdb_file_builder.h"

#include <algorithm>
#include <deque>
#include <stack>
#include <string>
#include <utility>
#include <vector>

#include "gmml/internal/atom.h"
#include "gmml/internal/element.h"
#include "gmml/internal/graph.h"
#include "gmml/internal/pdb_file.h"
#include "gmml/internal/proteins.h"
#include "gmml/internal/residue.h"
#include "gmml/internal/structure.h"

namespace gmml {

using std::auto_ptr;
using std::deque;
using std::pair;
using std::stack;
using std::string;
using std::vector;

class BuildPdbFile {
  public:
    BuildPdbFile(const Structure& structure,
                 const PdbFileBuilder& builder)
            : structure_(structure), builder_(builder),
              atom_map_(structure.size(), -1),
              residue_map_(structure.residue_count(), -1) {}

    PdbFile *operator()() {
        file_ = new PdbFile;

        build_atom_section();
        build_connect_section();
        build_link_section();

        file_->insert_card(new PdbEndCard);

        return file_;
    }

  private:
    void build_atom_section() {
        add_atoms();
        add_hetatms();
    }

    void add_atoms() {
        vector<vector<int>*> *proteins = find_proteins(structure_);
        for (int i = 0; i < proteins->size(); i++) {
            visit_amino_acid_chain(*(*proteins)[i]);
            delete (*proteins)[i];
            file_->insert_card(new PdbTerCard);
        }
        delete proteins;
    }

    void visit_amino_acid_chain(const vector<int>& indices) {
        for (int i = 0; i < indices.size(); i++)
            visit_residue(indices[i]);
    }

    void visit_residue(int index) {
        const Residue *residue = structure_.residues(index);
        residue_sequence_.push_back(index);
        residue_map_[index] = residue_sequence_.size();

        int start_atom = structure_.get_atom_index(index, 0);
        for (int i = start_atom; i < start_atom + residue->size(); i++) {
            if (structure_.atoms(i)->element() == Element("H") &&
                    !builder_.hydrogens_included()) {
                continue;
            }
            visit_atom(i, index);
        }
    }

    void visit_atom(int index, int residue_index) {
        atom_sequence_.push_back(index);
        int file_index = atom_sequence_.size();
        atom_map_[index] = file_index;
        const Atom *atom = structure_.atoms(index);
        const Residue *residue = structure_.residues(residue_index);
        PdbAtomCardBuilder builder;
        builder.initialize_from_atom(*atom);
        builder.set_serial(file_index);
        builder.set_res_name(residue->name());
        builder.set_res_seq(residue_sequence_.size());
        if (!is_amino_acid_code(residue->name())) {
            builder.set_hetatm(true);
        }
        file_->insert_card(builder.build());
    }

    bool is_amino_acid_code(const string& name) const {
        static AminoAcidCodeSet amino_acids;
        return amino_acids.lookup(name);
    }

    void add_hetatms() {
        int residue_count = structure_.residue_count();
        Graph *link_graph = structure_.get_link_graph();
        deque<bool> marked(residue_count, false);
        for (int i = 0; i < residue_sequence_.size(); i++)
            marked[residue_sequence_[i]] = true;
        stack<int> st;
        int molecule_start = 0;
        while (true) {
            if (st.empty()) {
                while (molecule_start < residue_count && marked[molecule_start])
                    molecule_start++;
                if (molecule_start == residue_count)
                    break;
                st.push(molecule_start);
                marked[molecule_start] = true;
            }
            int cur_residue = st.top();
            st.pop();
            visit_residue(cur_residue);
            const Graph::AdjList& adj = link_graph->edges(cur_residue);
            for (int i = 0; i < adj.size(); i++) {
                if (!marked[adj[i]]) {
                    marked[adj[i]] = true;
                    st.push(adj[i]);
                }
            }
        }
        delete link_graph;
    }

    void build_connect_section() {
        for (int i = 0; i < atom_sequence_.size(); i++) {
            const Structure::AdjList& adj_atoms =
                    structure_.bonds(atom_sequence_[i]);
            vector<int> bond_list;
            for (int j = 0; j < adj_atoms.size(); j++) {
                int mapped_atom = atom_map_[adj_atoms[j]];
                if (mapped_atom != -1 && test_connect_pair(atom_sequence_[i],
                                                           adj_atoms[j]))
                    bond_list.push_back(mapped_atom);
            }
            std::sort(bond_list.begin(), bond_list.end());

            vector<PdbConnectCard*> *cards =
                    PdbConnectCard::create_cards(i + 1, bond_list);

            for (int j = 0; j < cards->size(); j++)
                file_->insert_card((*cards)[j]);

            delete cards;
        }
    }

    bool test_connect_pair(int atom1_index, int atom2_index) const {
        static AminoAcidCodeSet amino_acids;
        const Atom *atom1 = structure_.atoms(atom1_index);
        const Atom *atom2 = structure_.atoms(atom2_index);
        const Residue *residue1 = structure_.residues(
                structure_.get_residue_index(atom1_index));
        const Residue *residue2 = structure_.residues(
                structure_.get_residue_index(atom2_index));
        if (amino_acids.lookup(residue1->name()) &&
                amino_acids.lookup(residue2->name())) {
            return atom1->element() == Element("S") &&
                   atom2->element() == Element("S");
        }
        return true;
    }

    void build_link_section() {
        AminoAcidCodeSet amino_acids;
        vector<size_t> *residue_index_table =
                structure_.get_residue_index_table();
        vector<PdbLinkCard*> cards;
        for (int i = 0; i < structure_.size(); i++) {
            const Structure::AdjList& adj_atoms = structure_.bonds(i);
            for (int j = 0; j < adj_atoms.size(); j++) {
                if (i > adj_atoms[j])
                    continue;
                int residue1_index = (*residue_index_table)[i];
                int residue2_index = (*residue_index_table)[adj_atoms[j]];
                if (residue1_index == residue2_index)
                    continue;
                string name1 = structure_.residues(residue1_index)->name();
                string name2 = structure_.residues(residue2_index)->name();
                if (amino_acids.lookup(name1) && amino_acids.lookup(name2)) {
                    continue;
                }
                cards.push_back(new PdbLinkCard(structure_.atoms(i)->name(),
                                name1,
                                residue_map_[residue1_index],
                                structure_.atoms(adj_atoms[j])->name(),
                                name2,
                                residue_map_[residue2_index]));
            }
        }
        delete residue_index_table;
        for (int i = cards.size() - 1; i >= 0; i--) {
            file_->insert_at_front(cards[i]);
        }
    }

    const Structure& structure_;
    const PdbFileBuilder& builder_;
    PdbFile *file_;
    // The atom indices in the order they're inserted into the file.
    vector<int> atom_sequence_;
    // The i'th element is the sequence id of the atom i in the structure.
    vector<int> atom_map_;
    vector<int> residue_sequence_;
    vector<int> residue_map_;
};

PdbFile *PdbFileBuilder::build(const Structure& structure) const {
    return BuildPdbFile(structure, *this)();
}

}  // namespace gmml
