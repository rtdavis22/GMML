// Author: Robert Davis

#ifndef GMML_INTERNAL_RESIDUE_H_
#define GMML_INTERNAL_RESIDUE_H_

#include <string>
#include <vector>

#include "gmml/internal/stubs/common.h"

namespace gmml {

class Atom;
class Graph;

class Residue {
  public:
    typedef std::vector<Atom*>::iterator iterator;
    typedef std::vector<Atom*>::const_iterator const_iterator;

    Residue();

    explicit Residue(const Residue *residue) { clone_from(residue); }

    explicit Residue(const std::string& name);

    Residue(const std::string& name, const std::vector<Atom*> *atoms,
            const Graph *bonds);

    virtual ~Residue();

    virtual Residue *clone() const {
        return new Residue(name_, &atoms_, bonds_);
    }

    void append(const Atom *atom);

    int get_index(const std::string& atom_name) const;

    //
    // Mutators
    //
    void set_name(const std::string& name) { name_ = name; }

    void set_atoms(const std::vector<Atom*> *atoms);

    void set_bonds(const Graph *bonds);

    void set_head(int head) { head_ = head; }

    void set_tail(int tail) { tail_ = tail; }

    //
    // Accessors
    //
    std::string name() const { return name_; }

    size_t size() const { return atoms_.size(); }

    Atom *atoms(int i) { return atoms_[i]; }
    const Atom *atoms(int i) const { return atoms_[i]; }

    iterator begin() { return atoms_.begin(); }
    const_iterator begin() const { return atoms_.begin(); }

    iterator end() { return atoms_.end(); }
    const_iterator end() const { return atoms_.end(); }

    const Graph *bonds() const { return bonds_; }

    int head() const { return head_; }
    int tail() const { return tail_; }

  protected:
    void clone_from(const Residue *residue);

  private:
    std::string name_;
    std::vector<Atom*> atoms_;
    Graph *bonds_;
    int head_;
    int tail_;

    DISALLOW_COPY_AND_ASSIGN(Residue);
};

class IndexedResidue : public Residue {
  public:
    IndexedResidue(const std::string& name, const std::vector<Atom*> *atoms,
                   const Graph *bonds) : Residue(name, atoms, bonds),
                                         indices_(atoms->size(), -1) {}

    explicit IndexedResidue(const Residue *residue)
            : Residue(residue), indices_(residue->size(), -1) {}

    explicit IndexedResidue(const std::string& name) : Residue(name) {}

    virtual ~IndexedResidue() {}

    void append(const Atom *atom, int index) {
        Residue::append(atom);
        indices_.push_back(index);
    }

    void set_atom_index(int atom, int index) { indices_[atom] = index; }
    void set_atom_index(const std::string& atom_name, int index);

    int get_atom_index(int atom) const { return indices_[atom]; }
    int get_atom_index(const std::string& atom_name) const;

  private:
    std::vector<int> indices_;
};

} //namespace gmml

#endif  // GMML_INTERNAL_RESIDUE_H_
