#ifndef STRUCTURE_H
#define STRUCTURE_H

#include <algorithm>
#include <list>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "boost/shared_ptr.hpp"

#include "atom.h"
#include "geometry.h"  // try to remove this
#include "graph.h"  // try to remove this
#include "residue.h"
#include "utilities.h"

namespace gmml {

class AmberTopFile;
class CoordinateFile;
struct MinimizationResults;
class ParameterFileSet;
class PdbFile;

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
    typedef AtomList::iterator iterator;
    typedef AtomList::const_iterator const_iterator;

    class InternalResidue {
      public:
        InternalResidue(const std::string& name,
                        AtomList::const_iterator begin,
                        AtomList::const_iterator end)
                : name_(name), begin_(begin), end_(end) {}

        std::string name() const { return name_; }
        //AtomList::iterator begin() { return begin_; }
        AtomList::const_iterator begin() const { return begin_; }
        //AtomList::iterator end() { return end_; }
        AtomList::const_iterator end() const { return end_; }

      private:
        std::string name_;
        AtomList::const_iterator begin_;
        AtomList::const_iterator end_;
    };

    Structure() : atoms_(0),
                  bonds_(new Graph),
                  residues_(new std::vector<detail::StructureResidue*>) {}

    virtual ~Structure() {
        delete bonds_;
        std::for_each(residues_->begin(), residues_->end(), DeletePtr());
        delete residues_;
    }

    iterator begin() { return atoms_.begin(); }
    const_iterator begin() const { return atoms_.begin(); }

    iterator end() { return atoms_.end(); }
    const_iterator end() const { return atoms_.end(); }

    static Structure *build_from_pdb(const PdbFile& pdb_file);

    virtual Structure *clone() const;

    void add_bond(int atom1, int atom2) { bonds_->add_edge(atom1, atom2); }

    virtual void append(const Structure& rhs);

    virtual int append(const Residue *residue);
    virtual int append(const std::string& prep_code);


    virtual int attach(Residue *new_residue, const std::string& new_atom_name,
                       int residue_index, const std::string& target_atom_name);

    // Returns the index of the residue that's attached
    virtual int attach(const std::string& prep_code,
                       const std::string& new_atom_name,
                       int residue_index,
                       const std::string& target_atom_name);

    MinimizationResults *minimize(const std::string& input_file);

    void translate_residue(int residue_index, double x, double y, double z);

    void remove_residue(int index);

    void set_dihedral(size_t atom1, size_t atom2, size_t atom3, size_t atom4,
                      double degrees);
    void set_dihedral(int residue1_index, const std::string& atom1_name,
                      int residue2_index, const std::string& atom2_name,
                      int residue3_index, const std::string& atom3_name,
                      int residue4_index, const std::string& atom4_name,
                      double degrees);

    // These make assumptions about the names of the atoms. This should be
    // changed.
    bool set_phi(int residue_index, double degrees);
    bool set_psi(int residue_index, double degrees);
    bool set_omega(int residue_index, double degrees);

    Graph *get_residue_link_table() const;

    int get_atom_index(int residue_index, int atom_index) const;
    int get_atom_index(int residue_index, const std::string& atom_name) const;

    // Returns the index of the anomeric atom of the residue. This works by
    // looking for a non-oxygen atom which is bound to an oxygen from another
    // residue. If no such atom is found, -1 is returned.
    int get_anomeric_index(int residue_index) const;

    // Returns the index of the oxygen from the previous function. -1 is
    // returned if there is no such oxygen.
    int get_parent_atom(int residue_index) const;

    // Returns the residue index of the oxygen atom in the previous function.
    // -1 is returned if there is no such atom.
    int get_parent_residue(int residue_index) const;

    // This returns a list of residue indices that are attached flexibly.
    // That is, their parent atom is an oxygen which is adjacent to an
    // exocyclic carbon of the same residue.
    std::vector<int> *get_flexible_linkages() const;

    // Return true if the atom is part of a ring.
    bool is_cyclic(int atom_index) const;

    size_t get_residue_index(size_t atom_index) const;
    size_t get_residue_count() const { return residues_->size(); }
    size_t get_residue_size(size_t i) const { return residues_->at(i)->size; }
    size_t get_residue_start(size_t i) const {
        return residues_->at(i)->start_index;
    }
    std::string get_residue_name(size_t i) const {
        return residues_->at(i)->name;
    }

    size_t size() const { return atoms_.size(); }
    const AtomList& atoms() const { return atoms_; }
    const AtomPtr atoms(size_t index) const { return atoms_[index]; }
    AtomPtr atoms(size_t index) { return atoms_[index]; }
    const Graph *bonds() const { return bonds_; }
    const AdjList& bonds(size_t index) const { return bonds_->edges(index); }

    // File operations
    AmberTopFile *build_amber_top_file(const ParameterFileSet& parm_set) const;
    AmberTopFile *build_amber_top_file() const;
    void print_amber_top_file(const std::string& file_name,
                              const ParameterFileSet& parm_set) const;
    void print_amber_top_file(const std::string& file_name) const;
    void print_amber_top_file() const;

    CoordinateFile *build_coordinate_file() const;
    void print_coordinate_file(const std::string& file_name) const;
    void print_coordinate_file() const;
    void load_coordinates(const CoordinateFile& coordinate_file);

    void set_residue_angle(int atom1, int atom2, int atom3, int residue_index,
                           double radians);
    void set_name(const std::string& name) { name_ = name; }

    std::auto_ptr<InternalResidue> residues(int index) const {
        AtomList::const_iterator first = atoms_.begin() +
                                         residues_->at(index)->start_index;
        return std::auto_ptr<InternalResidue>(
            new InternalResidue(residues_->at(index)->name,
                                first,
                                first + residues_->at(index)->size));
    }

  protected:
    // vector[i] is the residue index of atom i
    std::vector<size_t> *get_residue_index_table() const;

    void clone_from(const Structure& structure);

    AtomList atoms_;
    Graph *bonds_;
    std::vector<detail::StructureResidue*> *residues_;
    std::string name_;

    // bump up
    struct IndexedAtom {
        IndexedAtom(AtomPtr atom, int index) : atom(atom), index(index) {}
        AtomPtr atom;
        int index;
    };
    struct IndexedResidue {
        IndexedResidue() : name("") {}
        std::vector<IndexedAtom*> atoms;
        std::string name;
    };
    //
    void add_indexed_residue(std::map<int, int>& atom_map,
                             const IndexedResidue& residue);

  private:
    DISALLOW_COPY_AND_ASSIGN(Structure);
};

struct StructureAttach {
    int operator()(Structure& structure, Residue *residue,
                    const std::string& new_atom_name, int residue_index,
                    const std::string& atom_name) const;

    Vector<3> get_connection_direction(const Structure& structure,
                                       int source_index,
                                       int target_index) const;

    void set_dihedrals(Structure& structure, int new_residue_index,
                       int target_residue_index, int carbon_number,
                       int oxygen_number) const;
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
