// Author: Robert Davis

#include "gmml/internal/pdb_file_builder.h"

#include <algorithm>
#include <deque>
#include <stack>
#include <string>
#include <utility>

#include "gmml/internal/atom.h"
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

using boost::shared_ptr;

struct BuildPdbFile::Impl {
    Impl(const Structure& structure, const PdbFileBuilderBuilder& builder)
            : structure_(structure), builder_(builder),
              atom_map_(structure.size(), -1),
              residue_map_(structure.residue_count(), -1) {}

    void build_atom_section() {
        add_atoms();
        add_hetatms();
    }

    void add_atoms() {
        vector<vector<int>*> *proteins = find_proteins(structure_);
        for (int i = 0; i < proteins->size(); i++) {
            visit_amino_acid_chain(*(*proteins)[i]);
            delete (*proteins)[i];
            file_->insert_card(PdbFile::CardPtr(new PdbTerCard));
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
            if (structure_.atoms(i)->element() == kElementH &&
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
        PdbFile::AtomCardPtr card_ptr(new PdbAtomCard(
                file_index, atom->name(), residue->name(),
                residue_sequence_.size(),
                atom->coordinate().x, atom->coordinate().y,
                atom->coordinate().z,
                get_element_symbol(atom->element())));
        if (!is_amino_acid_code(residue->name())) {
            card_ptr->set_hetatm();
        }
        file_->insert_atom_card(card_ptr);
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
            int complete_rows = bond_list.size()/4;
            for (int j = 0; j < complete_rows; j++) {
                file_->insert_connect_card(PdbFile::ConnectCardPtr(
                        new PdbConnectCard(i + 1,
                                           bond_list[j*4],
                                           bond_list[j*4 + 1], 
                                           bond_list[j*4 + 2],
                                           bond_list[j*4 + 3])
                ));
            }
           
            PdbConnectCard card;
            card.connect1 = i + 1;
            // The fall-through is intentional.
            switch (bond_list.size()%4) {
                case 0:
                    break;
                case 3:
                    card.connect4 = bond_list[complete_rows*4 + 2];
                case 2:
                    card.connect3 = bond_list[complete_rows*4 + 1];
                case 1:
                    card.connect2 = bond_list[complete_rows*4];
                    file_->insert_connect_card(PdbFile::ConnectCardPtr(
                            new PdbConnectCard(card)));
                    break;
            }
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
            return atom1->element() == kElementS &&
                   atom2->element() == kElementS;
        }
        return true;
    }

    void build_link_section() {
        AminoAcidCodeSet amino_acids;
        vector<size_t> *residue_index_table =
                structure_.get_residue_index_table();
        vector<shared_ptr<PdbLinkCard> > cards;
        for (int i = 0; i < structure_.size(); i++) {
            const Structure::AdjList& adj_atoms = structure_.bonds(i);
            for (int j = 0; j < adj_atoms.size(); j++) {
                if (i > adj_atoms[i])
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
                cards.push_back(shared_ptr<PdbLinkCard>(
                        new PdbLinkCard(structure_.atoms(i)->name(),
                                        name1,
                                        residue_map_[residue1_index],
                                        structure_.atoms(adj_atoms[j])->name(),
                                        name2,
                                        residue_map_[residue2_index])));
            }
        }
        delete residue_index_table;
        for (int i = cards.size() - 1; i >= 0; i--) {
            file_->insert_at_front(cards[i]);
        }
    }

    const Structure& structure_;
    const PdbFileBuilderBuilder& builder_;
    PdbFile *file_;
    // The atom indices in the order they're inserted into the file.
    vector<int> atom_sequence_;
    // The i'th element is the sequence id of the atom i in the structure.
    vector<int> atom_map_;
    vector<int> residue_sequence_;
    vector<int> residue_map_;
};

BuildPdbFile::BuildPdbFile(const Structure& structure,
                           const PdbFileBuilderBuilder& builder)
        : impl_(new Impl(structure, builder)) {
}

BuildPdbFile::~BuildPdbFile() {
}

PdbFile *BuildPdbFile::operator()() {
    impl_->file_ = new PdbFile;

    impl_->build_atom_section();
    impl_->build_connect_section();
    impl_->build_link_section();
    impl_->file_->insert_card(PdbFile::CardPtr(new PdbEndCard));

    return impl_->file_;
}








/*

PdbFile *PdbFileBuilder::build(const Structure& structure) {
    PdbFile *file = new PdbFile;

    //build_link_section(structure);
    pair<vector<int>*, vector<int>*> atom_ret =
        build_atom_section(file, structure);
    build_connect_section(file, structure, *atom_ret.first, *atom_ret.second);

    delete atom_ret.first;
    delete atom_ret.second;

    return file;
}


void PdbFileBuilder::build_link_section(PdbFile *file,
                                        const Structure& structure) {
}

pair<vector<int>*, vector<int>*> PdbFileBuilder::build_atom_section(
        PdbFile *file, const Structure& structure) {
    vector<int> *sequence = new vector<int>;
    sequence->reserve(structure.size());
    vector<int> *atom_map = new vector<int>(structure.size(), -1);

    Graph *link_graph = structure.get_link_graph();
    
    int residue_count = structure.residue_count();
    deque<bool> marked(residue_count, false);
    marked[0] = true;
    stack<int> st;
    st.push(0);
    int cur_atom = 1;
    int cur_residue = 1;
    // This is the index of the first residue in the current molecule.
    int molecule_start = 0;
    while (true) {
        // If the stack is empty, that means we're done with the current
        // molecule, so we'll look for the next one.
        if (st.empty()) {
            while (molecule_start < residue_count && marked[molecule_start])
                molecule_start++;
            if (molecule_start == residue_count)
                break;
            st.push(molecule_start);
            marked[molecule_start] = true;
        }
        int top = st.top();
        const Residue *residue = structure.residues(top);
        for (Residue::const_iterator it = residue->begin();
                it != residue->end(); ++it) {
            //int atom_index = std::distance(structure.begin(), it);
            int atom_index = structure.get_atom_index(
                    top, std::distance(residue->begin(), it));
            (*atom_map)[atom_index] = cur_atom;
            sequence->push_back(atom_index);
            file->insert_atom_card(PdbFile::AtomCardPtr(new PdbAtomCard(
                cur_atom++, (*it)->name(), residue->name(), cur_residue,
                (*it)->coordinate().x, (*it)->coordinate().y,
                (*it)->coordinate().z,
                get_element_symbol((*it)->element()))));
        }
        st.pop();
        bool is_terminal = true;
        const Graph::AdjList& adj_list = link_graph->edges(top);
        for (int i = 0; i < adj_list.size(); i++) {
            if (!marked[adj_list[i]]) {
                marked[adj_list[i]] = true;
                st.push(adj_list[i]);
                is_terminal = false;
            }
        }
        if (is_terminal)
            file->insert_card(PdbFile::CardPtr(new PdbTerCard));
        cur_residue++;
    }

    delete link_graph;

    return std::make_pair(sequence, atom_map);
}

void PdbFileBuilder::build_connect_section(PdbFile *file,
                                           const Structure& structure,
                                           const vector<int>& sequence,
                                           const vector<int>& atom_map) {
    for (int i = 0; i < sequence.size(); i++) {
        Structure::AdjList row = structure.bonds(sequence[i]);
        for (int j = 0; j < row.size(); j++)
            row[j] = atom_map[row[j]];
        std::sort(row.begin(), row.end());
        int complete_rows = row.size()/4;
        for (int j = 0; j < complete_rows; j++) {
            file->insert_connect_card(PdbFile::ConnectCardPtr(
                new PdbConnectCard(i + 1, row[j*4], row[j*4 + 1], 
                                   row[j*4 + 2], row[j*4 + 3])
            ));
        }
           
        PdbConnectCard card;
        card.connect1 = i + 1;
        // The fall-through is intentional.
        switch (row.size()%4) {
          case 0:
            break;
          case 3:
            card.connect4 = row[complete_rows*4 + 2];
          case 2:
            card.connect3 = row[complete_rows*4 + 1];
          case 1:
            card.connect2 = row[complete_rows*4];
            file->insert_connect_card(PdbFile::ConnectCardPtr(
                    new PdbConnectCard(card)));
            break;
        }
    }
}
*/


}  // namespace gmml
