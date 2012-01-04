#ifndef GMML_INTERNAL_STRUCTURE_H_
#define GMML_INTERNAL_STRUCTURE_H_

#include <iosfwd>
#include <map>
#include <string>
#include <vector>

#include "gmml/internal/stubs/common.h"

namespace gmml {

class Atom;
class AmberTopFile;
class Box;
class CoordinateFile;
class Graph;
class IndexedResidue;
struct MinimizationResults;
class ParameterFileSet;
class PdbFile;
struct PdbMappingInfo;
class Residue;

class Structure {
  public:
    typedef std::vector<Atom*> AtomList;
    typedef AtomList::iterator iterator;
    typedef AtomList::const_iterator const_iterator;
    typedef std::vector<size_t> AdjList;

    Structure();

    virtual ~Structure();

    //
    // Construction methods
    //
    virtual Structure *clone() const;

    static Structure *build_from_pdb(const PdbFile& pdb_file,
                                     const PdbMappingInfo& mapping_info);

    // Build the Structure represented by the pdb file. This does not do
    // anything "smart", like infer bonding.
    static Structure *build_from_pdb(const PdbFile& pdb_file);

    // 
    static Structure *build_from_pdb(const PdbFile& pdb_file,
                                     bool infer_protein_bonds);

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
    void add_bond(int atom1, int atom2);
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
    virtual int append(const Structure *structure);

    // Maybe this should be private or something?
    virtual int append(const Residue *residue, bool load_bonds);

    virtual int append(const Residue *residue) { return append(residue, true); }

    // Add the atoms from the indexed residue to the structure and update the
    // map in which the keys of the map are the indices of the indexed atoms
    // and the values are the indices of the atoms in the structure.
    // The return value is the residue index of the residue in the structure.
    int append(std::map<int, int>& atom_map, const IndexedResidue *residue);

    // This builds the residue from a prep file in the prep file set of the
    // default environment and appends it. This should probably be changed to
    // look in library files (and more) as well. If the prep file code doesn't
    // exist in the default environment, -1 is returned.
    virtual int append(const std::string& name);

    // The attach operations reposition the atoms that are being attached to
    // the current structure according to StructureAttach below, and a bond
    // is created.
    // TODO: Should these really be virtual?
    // TODO: First arg should be a const Residue* for all, I think.
    //
    // Attach |head_name| of residue |new_residue| to |tail_name| of
    // residue |residue|. The index of the appended residue within the
    // structure is returned.

    // Attach residue's head atom to the structure's tail atom.
    int attach(Residue *residue);

    // Attach the atom in the residue with the given head name to the atom in
    // the given target residue of the structure with the specified tail name.
    int attach(Residue *residue, const std::string& head_name,
               int target_residue, const std::string& tail_name);

    // Attach the atom in the residue with the given head index to the
    // atom in the structure with the given tail index.
    int attach(Residue *new_residue, int head, int tail);

    // Attach the head atom in the residue to the atom with the given tail
    // name in the specified target residue of the structure.
    int attach_from_head(Residue *new_residue, int target_residue,
                         const std::string& tail_name);

    // Attach the head atom in the residue to the given tail atom index of
    // the structure.
    int attach_from_head(Residue *new_residue, int tail_atom);

    // Attach the atom in the residue with the given name to the structure's
    // tail atom.
    int attach_to_tail(Residue *new_residue, const std::string& head_name);

    // Attach the atom in the residue with the given head index to the
    // structure's tail atom.
    int attach_to_tail(Residue *new_residue, int head_index);

/*
    int attach(Residue *new_residue, const std::string& head_name,
               int residue, const std::string& tail_name);

    int attach(const std::string& code, int residue);

    // Attach the head atom of the new residue to the specified tail atom in
    // the specified target residue.
    int attach_from_head(Residue *new_residue, int target_residue,
                         const std::string& tail_name);

*/




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
    // Note: do away with this and just use head() and tail().
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
    const Residue *residues(int index) const;

    Residue *residues(int index);

    int residue_count() const { return residues_.size(); }

    size_t get_residue_index(size_t atom_index) const;

    // Returns the residue index of the oxygen atom in get_parent_atom().
    // -1 is returned if there is no such atom.
    int get_parent_residue(int residue_index) const;

    //
    // Accessors
    //
    // Returns the number of atoms in the structure.
    size_t size() const { return atoms_.size(); }

    const AtomList& atoms() const { return atoms_; }

    const Atom *atoms(size_t index) const { return atoms_[index]; }
    Atom *atoms(size_t index) { return atoms_[index]; }

    const Graph *bonds() const { return bonds_; }
    const AdjList& bonds(size_t index) const;

    int head() const { return head_; }
    int tail() const { return tail_; }

    // This might be a good way to check if a solvent that's a Structure
    // has a box.
    virtual const Box *box() const { return NULL; }

    //
    // Mutators
    //
    void set_bonds(const Graph *bonds);

    void set_head(int head) { head_ = head; }
    void set_tail(int tail) { tail_ = tail; }

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

    // Write the coordinate file to the file with the given name.
    void print_coordinate_file(const std::string& file_name) const;

    // Write the coordinate file to stdout.
    void print_coordinate_file() const;

    // Assign the coordinates from the coordinate file to the atoms of the
    // structure. Other information in the coordinate file (box information,
    // velocites) are ignored.
    void load_coordinates(const CoordinateFile& coordinate_file);

    // Print a representation of the structure.
    void print(std::ostream&) const;

  protected:
    // This is an alternative to the assignment operator.
    void clone_from(const Structure& structure);

  private:
    struct InternalResidue;

    // The random access list of atoms of the structure. This is redundant with
    // the residues. It is here for efficiency.
    AtomList atoms_;

    // A graph representing the bonding information of the structure.
    Graph *bonds_;

    // The residues of the structure.
    std::vector<InternalResidue*> residues_;

    // The indices of the head and tail atoms.
    int head_;
    int tail_;

    DISALLOW_COPY_AND_ASSIGN(Structure);
};


// The is the default attachment functor. It repositions the residue and
// attaches it to the structure. If you need to make other modifications to
// the structure, it is recommended that you make a wrapper around this
// functor. EDIT: you should wrap around Structure::attach instead?
// TODO: This functor should be changed to be general enough to only
// depend on the hybridizations of the atoms.
class StructureAttach {
  public:
    // The return value is -1 if the attachment failed for any reason, such as
    // the atom names don't exist in the residue.
    //int operator()(Structure& structure, Residue *new_residue,
    //               const std::string& head_name, int residue,
    //               const std::string& tail_name) const;

    //int operator()(Structure& structure, Residue *new_residue,
    //               int head, int residue, int tail) const;

    //int operator()(Structure& structure, const std::string& code,
    //               int residue) const;

    // Attach the head 
    //int operator()(Structure& structure, const std::string& code,
    //               int head) const;

    //int operator()(Structure& structure, const std::string& code,
    //               const std::string& head_name) const;

    // Attach the head atom of the new residue to the structure's tail atom.
    int operator()(Structure& structure, Residue *new_residue,
                   int head, int tail) const;

  private:
    struct Impl;
};

}  // namespace gmml

#include "structure-inl.h"

#endif  // GMML_INTERNAL_STRUCTURE_H_
