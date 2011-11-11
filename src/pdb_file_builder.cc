// Author: Robert Davis

#include "gmml/internal/pdb_file_builder.h"

#include <algorithm>
#include <deque>
#include <memory>
#include <stack>
#include <string>
#include <utility>

#include "gmml/internal/graph.h"
#include "gmml/internal/pdb_file.h"
#include "gmml/internal/structure.h"

namespace gmml {

using std::auto_ptr;
using std::deque;
using std::pair;
using std::stack;
using std::vector;

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
    
    int residue_count = structure.get_residue_count();
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
        auto_ptr<Structure::InternalResidue> residue = structure.residues(top);
        for (Structure::AtomList::const_iterator it = residue->begin();
                it != residue->end(); ++it) {
            int atom_index = std::distance(structure.begin(), it);
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

}  // namespace gmml
