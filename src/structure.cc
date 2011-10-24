#include "gmml/internal/structure.h"

#include <algorithm>
#include <vector>

#include "gmml/internal/geometry.h"
#include "gmml/internal/utilities.h"

using std::vector;

namespace gmml {

Structure *Structure::clone() const {
    using detail::StructureResidue;

    Structure *structure = new Structure;
    structure->atoms_.reserve(atoms_.size());
    for (size_t i = 0; i < atoms_.size(); i++)
	structure->atoms_.push_back(AtomPtr(atoms_[i]->clone()));

    structure->residues_->reserve(residues_->size());
    for (size_t i = 0; i < residues_->size(); i++)
	structure->residues_->push_back(
            new StructureResidue(*residues_->at(i))
        );

    structure->bonds_ = bonds_->clone();

    return structure;
}

void Structure::append(const Structure& rhs) {
    size_t old_size = atoms_.size();
    atoms_.reserve(old_size + rhs.atoms_.size());

    //std::copy(rhs.atoms_->begin(), rhs.atoms_->end(), atoms_->end());

    bonds_->append(*rhs.bonds_);

    using detail::StructureResidue;
    residues_->reserve(residues_->size() + rhs.residues_->size());
    for (size_t i = 0; i < rhs.residues_->size(); i++) {
	residues_->push_back(
	    new StructureResidue(rhs.residues_->at(i)->name,
				 rhs.residues_->at(i)->start_index + old_size,
				 rhs.residues_->at(i)->size)
	);
    }
}

void Structure::append(const Residue *residue) {
    size_t old_size = atoms_.size();
    atoms_.reserve(old_size + residue->size());
    atoms_.insert(atoms_.end(), residue->begin(), residue->end());
    bonds_->append(*residue->bonds());
    using detail::StructureResidue;
    residues_->push_back(new StructureResidue(residue->name(), old_size,
                                              residue->size())); 
}

void Structure::set_dihedral(size_t atom1, size_t atom2, size_t atom3,
                             size_t atom4, double degrees) {
    vector<size_t> *atoms = bonds_->edge_bfs(atom2, atom3);
    double current_dihedral = 
        measure(atoms_[atom1]->coordinate(), 
                atoms_[atom2]->coordinate(),
                atoms_[atom3]->coordinate(), 
                atoms_[atom4]->coordinate());
    RotationMatrix matrix(atoms_[atom2]->coordinate(),
                          Vector<3>(atoms_[atom3]->coordinate(),
                                    atoms_[atom2]->coordinate()),
                          current_dihedral - to_radians(degrees));
    for (vector<size_t>::iterator it = atoms->begin();
            it != atoms->end(); ++it)
        matrix.apply(atoms_[*it]->coordinate());

    delete atoms;
}

Graph *Structure::get_residue_link_table() const {
    vector<size_t> *residue_index_table = get_residue_index_table();

    Graph *links = new Graph(residues_->size());
    for (size_t i = 0; i < atoms_.size(); i++) {
        for (size_t j = 0; j < bonds_->edges(i).size(); j++)
            if (residue_index_table->at(bonds_->edges(i)[j]) !=
                    residue_index_table->at(i))
                links->add_edge(residue_index_table->at(i),
                                residue_index_table->at(bonds_->edges(i)[j]));
    }

    delete residue_index_table;

    return links;
}

vector<size_t> *Structure::get_residue_index_table() const {
    vector<size_t> *table = new vector<size_t>(atoms_.size());
    size_t current_index = 0;
    for (size_t i = 0; i < residues_->size(); i++)
        for (size_t j = 0; j < residues_->at(i)->size; j++)
            (*table)[current_index++] = i;
    return table;
}

}
