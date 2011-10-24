#ifndef STRUCTURE_H
#define STRUCTURE_H

#include <algorithm>
#include <list>
#include <map>
#include <string>
#include <vector>

#include "boost/shared_ptr.hpp"

#include "atom.h"
#include "graph.h"
#include "residue.h"
#include "utilities.h"

namespace gmml {
namespace detail {

// put this in structure, probably (private)
struct StructureResidue {
    StructureResidue(const std::string& name, size_t start_index, size_t size)
            : name(name), start_index(start_index), size(size) {}

    std::string name;
    size_t start_index;
    size_t size;
};

}

class Structure {
  public:
    typedef boost::shared_ptr<Atom> AtomPtr;
    typedef std::vector<AtomPtr> AtomList;
    typedef Graph::AdjList AdjList;

    Structure() : atoms_(0),
                  bonds_(new Graph),
                  residues_(new std::vector<detail::StructureResidue*>) {}

    virtual ~Structure() {
        delete bonds_;
        std::for_each(residues_->begin(), residues_->end(), DeletePtr());
        delete residues_;
    }

    virtual Structure *clone() const;

    virtual void append(const Structure& rhs);

    virtual void append(const Residue *residue);

    void remove_residue(int index);

    void set_dihedral(size_t atom1, size_t atom2, size_t atom3, size_t atom4,
                      double degrees);

    Graph *get_residue_link_table() const;

    size_t get_residue_index(size_t atom_index) const;
    size_t get_residue_count() const { return residues_->size(); }
    size_t get_residue_size(size_t i) const { return residues_->at(i)->size; }
    std::string get_residue_name(size_t i) const {
        return residues_->at(i)->name;
    }

    size_t size() const { return atoms_.size(); }
    const AtomList& atoms() const { return atoms_; }
    const AtomPtr atoms(size_t index) const { return atoms_[index]; }
    const Graph *bonds() const { return bonds_; }
    const AdjList& bonds(size_t index) const { return bonds_->edges(index); }

  protected:
    // vector[i] is the residue index of atom i
    std::vector<size_t> *get_residue_index_table() const;

    AtomList atoms_;
    Graph *bonds_;
    std::vector<detail::StructureResidue*> *residues_;

  private:
    DISALLOW_COPY_AND_ASSIGN(Structure);
};

inline size_t Structure::get_residue_index(size_t atom_index) const {
    for (size_t i = 1; i < residues_->size(); i++) {
        if (atom_index < residues_->at(i)->start_index)
            return i - 1;
    }
    return residues_->size() - 1;
}

}  // namespace gmml

#endif  // STRUCTURE_H
