#include "gmml/internal/residue.h"

#include "gmml/internal/atom.h"
#include "gmml/internal/graph.h"

using std::string;
using std::vector;

namespace gmml {

Residue::Residue() : name_(""), bonds_(new Graph), head_(-1), tail_(-1) {}

Residue::Residue(const string& name) : name_(name), bonds_(new Graph),
                                       head_(-1), tail_(-1) {}

Residue::Residue(const string& name, const vector<Atom*> *atoms,
                 const Graph *bonds) : name_(name), head_(-1), tail_(-1) {
    set_atoms(atoms);
    if (bonds != NULL)
        bonds_ = bonds->clone();
    else
        bonds_ = new Graph(size());
}

Residue::~Residue() {
    STLDeleteContainerPointers(atoms_.begin(), atoms_.end());
    if (bonds_ != NULL)
        delete bonds_;
}

void Residue::append(const Atom *atom) {
    atoms_.push_back(atom->clone());
}

int Residue::get_index(const string& atom_name) const {
    for (const_iterator it = begin(); it != end(); ++it) {
        if ((*it)->name() == atom_name) {
            return std::distance(begin(), it);
        }
    }
    return -1;
}

void Residue::set_atoms(const vector<Atom*> *atoms) {
    STLDeleteContainerPointers(atoms_.begin(), atoms_.end());
    atoms_.clear();
    std::vector<Atom*>::const_iterator it = atoms->begin();
    while (it != atoms->end()) {
        atoms_.push_back((*it)->clone());
        ++it;
    }
}

void Residue::set_bonds(const Graph *bonds) {
    if (bonds_ != NULL)
        delete bonds_;
    if (bonds != NULL)
        bonds_ = bonds->clone();
    else
        bonds_ = NULL;
}

void Residue::clone_from(const Residue *residue) {
    name_ = residue->name();
    set_atoms(&residue->atoms_);
    set_bonds(residue->bonds_);
    head_ = residue->head();
    tail_ = residue->tail();
}

void IndexedResidue::set_atom_index(const string& atom_name, int index) {
    for (int i = 0; i < static_cast<int>(size()); i++) {
        if (atoms(i)->name() == atom_name) {
            indices_[i] = index;
            break;
        }
    }
}

int IndexedResidue::get_atom_index(const string& atom_name) const {
    for (int i = 0; i < static_cast<int>(size()); i++) {
        if (atoms(i)->name() == atom_name) {
            return indices_[i];
        }
    }
    return -1;
}

}  // namespace gmml
