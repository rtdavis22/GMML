#include "gmml/internal/pdb_file_builder.h"

#include <algorithm>
#include <stack>
#include <string>

#include "gmml/internal/pdb_file.h"
#include "gmml/internal/structure.h"

namespace gmml
{

using std::stack;
using std::vector;

namespace
{
PdbFile *file;
}

PdbFile *PdbFileBuilder::build(const Structure& structure) {
    file = new PdbFile();

    //build_link_section(structure);
    build_atom_section(structure);
    build_connect_section(structure);

    return file;
}


void PdbFileBuilder::build_link_section(const Structure& structure) {
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
void PdbFileBuilder::build_atom_section(const Structure& structure) {
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

void PdbFileBuilder::build_connect_section(const Structure& structure) {
/*
    const vector<vector<int> >& bond_table = structure.bond_table();
    for (int i = 0; i < bond_table.size(); i++) {
        vector<int> row = bond_table[i];
        std::sort(row.begin(), row.end());
        int complete_rows = row.size()/4;
        for (int j = 0; j < complete_rows; j++) {
            file->insert_card(PdbFile::CardPtr(
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
            file->insert_card(PdbFile::CardPtr(new PdbConnectCard(card)));
            break;
        }
    }
*/
}

} //namespace gmml
