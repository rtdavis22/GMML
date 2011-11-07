#include "gmml/internal/pdb_file_builder.h"

#include <algorithm>
#include <memory>
#include <stack>
#include <string>

#include "gmml/internal/pdb_file.h"
#include "gmml/internal/structure.h"

namespace gmml {

using std::auto_ptr;
using std::stack;
using std::vector;

PdbFile *PdbFileBuilder::build(const Structure& structure) {
    PdbFile *file = new PdbFile;

    //build_link_section(structure);
    build_atom_section(file, structure);
    build_connect_section(file, structure);

    return file;
}


void PdbFileBuilder::build_link_section(PdbFile *file,
                                        const Structure& structure) {
/*
    const vector<Atom*> atoms = structure.atoms();
    const vector<vector<int> >& bond_table = structure.bond_table();
    vector<PdbLinkCard> link_cards;
    for (int i = 0; i < bond_table.size(); i++) {
       vector<int> links = bond_table[i];
       for (int j = 0; j < links.size(); j++) {
           if (i >= links[j]) continue;
           int from = structure.get_residue_index(i);
           int to = structure.get_residue_index(links[j]);
           if (from != to) {
               const Residue& start_residue = structure.get_residue(from);
               const Residue& end_residue = structure.get_residue(to);
               PdbLinkCard new_card(atoms[i]->name, start_residue.name(),
                                    from + 1, atoms[links[j]]->name,
                                    end_residue.name(), to + 1);
               link_cards.push_back(new_card);
           }
       }
    }
    std::sort(link_cards.begin(), link_cards.end());
    for (int i = 0; i < link_cards.size(); i++)
        file->insert_card(PdbFile::CardPtr(new PdbLinkCard(link_cards[i])));
*/
}

//add multiple components
void PdbFileBuilder::build_atom_section(PdbFile *file,
                                        const Structure& structure) {
    int cur_index = 1;
    for (int i = 0; i < structure.get_residue_count(); i++) {
        auto_ptr<Structure::InternalResidue> residue =
            structure.residues(i);
        for (Structure::AtomList::const_iterator it = residue->begin();
                it != residue->end(); ++it) {
            file->insert_atom_card(PdbFile::AtomCardPtr(new PdbAtomCard(
                    cur_index++, (*it)->name(), residue->name(),
                    i, (*it)->coordinate().x, (*it)->coordinate().y,
                    (*it)->coordinate().z,
                    get_element_symbol((*it)->element()))));
        }
    }
/*
    for (int i = 0; i < structure.size(); i++) {
        const Structure::AtomPtr atom = structure.atoms(i);
        PdbFile::AtomCardPtr(new PdbAtomCard(
                i, atom->name(), structure.get_residue_name(i),
        file->insert_atom_card(PdbFile::Card(i, atom->name(),
                 
    }
*/
/*
    vector<vector<int> > *link_table = structure.get_residue_link_table();

    vector<bool> marked(link_table->size());
    marked[0] = true;
    stack<int> st;
    st.push(0);
    int current_atom_number = 1;
    while (!st.empty()) {
        int front = st.top();
        const Residue& residue = structure.get_residue(front);
        vector<Atom*> atoms = residue.atoms();
        for (int i = 0; i < atoms.size(); i++) {
            file->insert_card(PdbFile::CardPtr(
                new PdbAtomCard(current_atom_number++, atoms[i]->name,
                                residue.name(), front + 1, 
                                atoms[i]->coordinate.x, atoms[i]->coordinate.y, 
                                atoms[i]->coordinate.z, atoms[i]->element)
            ));
        }
        st.pop();
        bool is_chain_terminal = true;
        for (int i = (*link_table)[front].size() - 1; i >= 0; i--) {
            if (!marked[(*link_table)[front][i]]) {
                marked[(*link_table)[front][i]] = true;
                st.push((*link_table)[front][i]);
                is_chain_terminal = false;
            }
        }
        if (is_chain_terminal)
            file->insert_card(PdbFile::CardPtr(new PdbTerCard()));
    }

    delete link_table;
*/
}

void PdbFileBuilder::build_connect_section(PdbFile *file,
                                           const Structure& structure) {
    for (int i = 0; i < structure.size(); i++) {
        const Structure::AdjList& row = structure.bonds(i);
        int complete_rows = row.size()/4;
        for (int j = 0; j < complete_rows; j++) {
            file->insert_connect_card(PdbFile::ConnectCardPtr(
                new PdbConnectCard(i + 1, row[j*4] + 1, row[j*4 + 1] + 1, 
                                   row[j*4 + 2] + 1, row[j*4 + 3] + 1)
            ));
        }
           
        PdbConnectCard card;
        card.connect1 = i + 1;
        //note the fall-through
        switch (row.size()%4) {
          case 0:
            break;
          case 3:
            card.connect4 = row[complete_rows*4 + 2] + 1;
          case 2:
            card.connect3 = row[complete_rows*4 + 1] + 1;
          case 1:
            card.connect2 = row[complete_rows*4] + 1;
            file->insert_connect_card(PdbFile::ConnectCardPtr(
                    new PdbConnectCard(card)));
            break;
        }
    }
}

}  // namespace gmml
