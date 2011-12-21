#ifndef GMML_INTERNAL_STRUCTURE_H_
#define GMML_INTERNAL_STRUCTURE_H_

#include <algorithm>
#include <iosfwd>
#include <list>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "boost/shared_ptr.hpp"

#include "gmml/internal/atom.h"
#include "gmml/internal/graph.h"
#include "gmml/internal/residue.h"
#include "gmml/internal/stubs/common.h"

namespace gmml {

class AmberTopFile;
class CoordinateFile;
struct MinimizationResults;
class ParameterFileSet;
class PdbFile;
struct PdbMappingInfo;

class Structure {
  public:
    typedef boost::shared_ptr<Atom> AtomPtr;
    typedef std::vector<AtomPtr> AtomList;
    typedef Graph::AdjList AdjList;
    typedef AtomList::iterator iterator;
    typedef AtomList::const_iterator const_iterator;

    // This class provides a convenient representation of a residue in the
    // structure.
    class InternalResidue {
      public:
        InternalResidue(const std::string& name,
                        AtomList::const_iterator begin,
                        AtomList::const_iterator end)
                : name_(name), begin_(begin), end_(end) {}

        std::string name() const { return name_; }
        AtomList::const_iterator begin() const { return begin_; }
        AtomList::const_iterator end() const { return end_; }

      private:
        std::string name_;
        AtomList::const_iterator begin_;
        AtomList::const_iterator end_;
    };

    Structure() : atoms_(0),
                  bonds_(new Graph),
                  residues_(new std::vector<StructureResidue*>) {}

    virtual ~Structure() {
        delete bonds_;
        std::for_each(residues_->begin(), residues_->end(), DeletePtr());
        delete residues_;
    }

    //
    // Construction methods
    //
    static Structure *build_from_pdb(const PdbFile& pdb_file,
                                     const PdbMappingInfo& mapping_info);

    // Build the Structure represented by the pdb file. This does not do
    // anything "smart", like infer bonding.
    static Structure *build_from_pdb(const PdbFile& pdb_file);

    // 
    static Structure *build_from_pdb(const PdbFile& pdb_file,
                                     bool infer_protein_bonds);

    virtual Structure *clone() const;

    //
    // Iteration
    //
    // These are for iterating over the atoms.
    iterator begin() { return atoms_.begin(); }
    const_iterator begin() const { return atoms_.begin(); }

    iterator end() { return atoms_.end(); }
    const_iterator end() const { return atoms_.end(); }

    //
    // Bonding-related operations
    //
    // Add or remove the bond between the atoms with the given atom indices.
    void add_bond(int atom1, int atom2) { bonds_->add_edge(atom1, atom2); }
    void remove_bond(int atom1, int atom2);

    // Returns true if the atom with the given index is part of a cycle.
    bool is_cyclic(int atom) const;

    // This returns a graph respresenting the links between residues.
    Graph *get_residue_link_table() const;

    // To me: what's the difference between this and the above?
    Graph *get_link_graph() const;

    // This returns a list of residue indices that are attached flexibly.
    // That is, their parent atom is an oxygen which is adjacent to an
    // exocyclic carbon of the same residue.
    std::vector<int> *get_flexible_linkages() const;

    //
    // Geometric operations
    //
    // Shift the entire structure |x| units in the x direction, |y| units in
    // the y direction, and |z| units in the z direction.
    void shift(double x, double y, double z);

    // Translate a given residue.
    void translate_residue(int residue_index, double x, double y, double z);

    // Set the dihedral of the atoms with the given atom indices. Note that the
    // dihedral angle should be given in degrees.
    void set_dihedral(size_t atom1, size_t atom2, size_t atom3, size_t atom4,
                      double degrees);

    // This is similar to the above function, but atoms are identified by
    // their residue index and name.
    void set_dihedral(int residue1_index, const std::string& atom1_name,
                      int residue2_index, const std::string& atom2_name,
                      int residue3_index, const std::string& atom3_name,
                      int residue4_index, const std::string& atom4_name,
                      double degrees);

    // Set the angle in the structure, but only modify the atoms of the
    // given residue.
    // Note: change this to use degrees.
    void set_residue_angle(int atom1, int atom2, int atom3, int residue_index,
                           double radians);

    // The following three functions set the glycosidic angles between the
    // residue with the given residue index and it's parent (the residue at
    // the reducing end). It should be noted that these make some basic
    // assumptions about the names of some atoms. For example, the n'th carbon
    // should be named Cn.
    //
    // Phi: H1-C1-O-CX'
    //      C1-C2-O-CX'
    bool set_phi(int residue_index, double degrees);

    // Psi: C1-O-CX'-HX'
    //      C1-O-C6'-C5'
    bool set_psi(int residue_index, double degrees);

    // Omega: O-C6'-C5'-O5'
    bool set_omega(int residue_index, double degrees);

    //
    // Augmenting operations
    //
    // The append operations do not modify the positions of any atoms or change
    // bonding.
    //
    // Insert the atoms of the structure at the end of the atoms of this
    // structure. The atoms are not copied. You must explicitly Clone() the
    // structure if you want a copy of the atoms.
    virtual void append(const Structure& rhs);

    // This is similar to the above, but with a residue. The index of the new
    // residue in the structure is returned.
    virtual int append(const Residue *residue);

    // This builds the residue from a prep file in the prep file set of the
    // default environment and appends it. This should probably be changed to
    // look in library files (and more) as well. If the prep file code doesn't
    // exist in the default environment, -1 is returned.
    virtual int append(const std::string& prep_code);

    // The attach operations reposition the atoms that are being attached to
    // the current structure according to StructureAttach below, and a bond
    // is created.
    //
    // Attach |new_atom_name| of residue |new_residue| to |target_atom_name| of
    // residue |residue_index|. The index of the appended residue within the
    // structure is returned.
    virtual int attach(Residue *new_residue, const std::string& new_atom_name,
                       int residue_index, const std::string& target_atom_name);

    // This is similar to the above function, but the residue's prep file code
    // is given.
    virtual int attach(const std::string& prep_code,
                       const std::string& new_atom_name,
                       int residue_index,
                       const std::string& target_atom_name);

    //
    // Removal operations
    //
    // Remove the residue and any bonds involving the residue from the
    // structure.
    void remove_residue(int index);

    // Remove multiple residue from the structure. This is more efficient than
    // calling the above procedure with each residue index.
    void remove_residues(const std::vector<int>& residues);

    //
    // Advanced modification operations
    //
    // Minimize the structure with SANDER using the given SANDER input
    // minimization file and representation of the results of the
    // minimization. NULL is returned if the minimization failed for any
    // reason.
    MinimizationResults *minimize(const std::string& input_file);

    //
    // Atom-related query operations
    //
    // Returns the atom index within the structure, given a residue index and
    // the index of the atom within the residue. -1 is returned if the atom
    // is not found.
    int get_atom_index(int residue_index, int atom_index) const;

    // Returns the atom index within the structure, given a residue index and
    // the atom name. -1 is returned if the atom is not found.
    int get_atom_index(int residue_index, const std::string& atom_name) const;

    // Returns the index of the anomeric atom of the residue. This works by
    // looking for a non-oxygen atom which is bound to an oxygen from another
    // residue. If no such atom is found, -1 is returned.
    int get_anomeric_index(int residue_index) const;

    // Returns the index of the oxygen from the previous function. -1 is
    // returned if there is no such oxygen.
    int get_parent_atom(int residue_index) const;

    // The number at index i of the return vertex is the residue index of
    // atom i.
    std::vector<size_t> *get_residue_index_table() const;

    //
    // Residue-related query operations
    //
    // Return information associated with the residue at index |index|.
    std::auto_ptr<InternalResidue> residues(int index) const;

    // Returns the residue index of the oxygen atom in get_parent_atom().
    // -1 is returned if there is no such atom.
    int get_parent_residue(int residue_index) const;

    // THESE ARE DEPRECATED!
    size_t get_residue_index(size_t atom_index) const;
    size_t get_residue_count() const { return residues_->size(); }
    size_t get_residue_size(size_t i) const { return residues_->at(i)->size; }
    size_t get_residue_start(size_t i) const {
        return residues_->at(i)->start_index;
    }
    std::string get_residue_name(size_t i) const {
        return residues_->at(i)->name;
    }
    // DONT USE THE ABOVE

    //
    // Accessors
    //
    // Returns the number of atoms in the structure.
    size_t size() const { return atoms_.size(); }

    const AtomList& atoms() const { return atoms_; }

    const AtomPtr atoms(size_t index) const { return atoms_[index]; }
    AtomPtr atoms(size_t index) { return atoms_[index]; }

    const Graph *bonds() const { return bonds_; }
    const AdjList& bonds(size_t index) const { return bonds_->edges(index); }

    //
    // File operations
    //
    // Return a pdb file that represents this structure.
    virtual PdbFile *build_pdb_file() const;

    // Write the pdb file representing this structure to a file with the
    // given name.
    void print_pdb_file(const std::string& file_name) const;

    // Write the pdb file to stdout.
    void print_pdb_file() const;

    // Return the AMBER topology file that represents this structure, using
    // the given parameter set.
    AmberTopFile *build_amber_top_file(const ParameterFileSet& parm_set) const;

    // Return the AMBER topology file that represents this structure, using
    // the parameter set in the default environment.
    virtual AmberTopFile *build_amber_top_file() const;

    // Write the AMBER topology file to a file with the given name, using
    // the given parameter set.
    void print_amber_top_file(const std::string& file_name,
                              const ParameterFileSet& parm_set) const;

    // Write the AMBER topology file to a file with the given name, using
    // the parameter set in the default environment.
    void print_amber_top_file(const std::string& file_name) const;

    // Write the AMBER topology file to stdout, using the parameter set in the
    // default environment.
    void print_amber_top_file() const;

    // Return a coordinate file with coordinates of the atoms from this
    // structure in the order they are in the structure.
    virtual CoordinateFile *build_coordinate_file() const;

    // Print the coordinate file to the file with the given name.
    void print_coordinate_file(const std::string& file_name) const;

    // Print the coordinate file to stdout.
    void print_coordinate_file() const;

    // Assign the coordinates from the coordinate file to the atoms of the
    // structure. Other information in the coordinate file (box information,
    // velocites) are ignored.
    void load_coordinates(const CoordinateFile& coordinate_file);

    // Print a representation of the structure.
    void print(std::ostream&) const;

  protected:
    // This is an internal representation of the residues in the structure.
    struct StructureResidue {
        StructureResidue(const std::string& name, size_t start_index,
                         size_t size)
                : name(name), start_index(start_index), size(size) {}

        std::string name;
        size_t start_index;
        size_t size;
    };

    // It is often useful to attach an index to an atom or residue, so these
    // are provided.
    struct IndexedAtom {
        IndexedAtom(AtomPtr atom, int index) : atom(atom), index(index) {}

        AtomPtr atom;
        int index;
    };
    struct IndexedResidue {
        IndexedResidue() : name(""), bonds(NULL) {}

        int get_atom_index_by_name(const std::string& name) const;

        std::vector<IndexedAtom*> atoms;
        std::string name;
        Graph *bonds;
    };

    // This is an alternative to the assignment operator.
    void clone_from(const Structure& structure);

    // Add the atoms from the indexed residue to the structure and update the
    // map in which the keys of the map are the indices of the indexed atoms
    // and the values are the indices of the atoms in the structure.
    // The return value is the residue index of the residue in the structure.
    int add_indexed_residue(std::map<int, int>& atom_map,
                            const IndexedResidue& residue);

    // This function attempts to fill in the missing atoms of the input residue
    // using the structure with the given name. If it is unsuccessful for any
    // reason, NULL is returned. If atoms are added to the input residue, their
    // index is -1.
    static void replace_residue(IndexedResidue *residue,
                                const std::string& name);

    // Find the longest chain (up to length 3) of non-NULL atoms that start
    // with the atom index.
    static std::vector<int> *find_chain(int index,
                                        const std::vector<IndexedAtom*>& atoms,
                                        const Graph *bonds);

    // The random access list of the atoms of the structure.
    AtomList atoms_;

    // A graph representing the bonding information of the structure.
    Graph *bonds_;

    // The internal representation of the residues within a structure.
    std::vector<StructureResidue*> *residues_;

  private:
    DISALLOW_COPY_AND_ASSIGN(Structure);
};


// The is the default attachment functor. It repositions the residue and
// attaches it to the structure. If you need to make other modifications to
// the structure, it is recommended that you make a wrapper around this
// functor.
class StructureAttach {
  public:
    // The return value is -1 if the attachment failed for any reason, such as
    // the atom names don't exist in the residue.
    int operator()(Structure& structure, Residue *residue,
                    const std::string& new_atom_name, int residue_index,
                    const std::string& atom_name) const;

  private:
    Vector<3> get_connection_direction(const Structure& structure,
                                       int source_index,
                                       int target_index) const;

    void set_dihedrals(Structure& structure, int new_residue_index,
                       int target_residue_index, int carbon_number,
                       int oxygen_number) const;
};

}  // namespace gmml

#include "structure-inl.h"

#endif  // GMML_INTERNAL_STRUCTURE_H_
