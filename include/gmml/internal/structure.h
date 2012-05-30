// Copyright (c) 2012 The University of Georgia
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.

// Author: Robert Davis

#ifndef GMML_INTERNAL_STRUCTURE_H_
#define GMML_INTERNAL_STRUCTURE_H_

#include <iosfwd>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

#include "gmml/internal/residue.h"
#include "gmml/internal/stubs/common.h"

namespace gmml {

class Atom;
class AmberTopFile;
class Box;
class CoordinateFile;
class Graph;
class IndexedResidue;
struct MinimizationResults;
class ParameterSet;
class PdbFile;
struct PdbMappingInfo;
class Residue;

class DihedralAtoms;

class Structure {
  public:
    typedef std::vector<Atom*> AtomList;
    typedef AtomList::iterator iterator;
    typedef AtomList::const_iterator const_iterator;
    typedef std::vector<size_t> AdjList;

    Structure();

    Structure(const Residue *residue);

    virtual ~Structure();

    //
    // Construction methods
    //
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
    // Returns true if the atoms with the given indices are bonded.
    bool is_bonded(int atom1, int atom2) const;

    // Add or remove the bond between the atoms with the given atom indices.
    bool add_bond(int atom1, int atom2);
    bool remove_bond(int atom1, int atom2);

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

    // Translate all residues after the given residue index.
    void translate_residues_after(int index, double x, double y, double z);

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

    void set_outside_dihedral(size_t atom1, size_t atom2, size_t atom3,
                              size_t atom4, double degrees);
    void set_outside_dihedral(int residue1, const std::string& atom1,
                              int residue2, const std::string& atom2,
                              int residue3, const std::string& atom3,
                              int residue4, const std::string& atom4,
                              double degrees);

    // This syntax experimental. It should be the function that does the work.
    void set_dihedral(const DihedralAtoms& atoms, double value) {
        
    }

    // Set the angle in the structure, but only modify the atoms of the
    // given residue.
    // Note: change this to use degrees.
    void set_residue_angle(int atom1, int atom2, int atom3, int residue_index,
                           double radians);

    void set_angle_in_range(int start_atom, int end_atom, int atom1, int atom2,
                            int atom3, double radians);

    // Modify all residues including and after the given residue index.
    void set_angle_after(int residue, int atom1, int atom2, int atom3,
                         double radians);

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
    //int append(std::map<int, int>& atom_map, const IndexedResidue *residue);

    template<class ATOM_MAP>
    int append(ATOM_MAP& atom_map, const IndexedResidue *residue) {
        int prev_size = size();
        int index = append(residue);
        for (int i = 0; i < residue->size(); i++) {
            typename ATOM_MAP::value_type entry(residue->get_atom_index(i),
                                                prev_size + i);
            atom_map.insert(entry);
        }
        return index;
    }

    // This builds the residue from a prep file in the prep file set of the
    // default environment and appends it. This should probably be changed to
    // look in library files (and more) as well. If the prep file code doesn't
    // exist in the default environment, -1 is returned.
    virtual int append(const std::string& name);

    // The attach operations reposition the atoms that are being attached to
    // the current structure according to StructureAttach below, and a bond
    // is created.
    //
    // Attach |head_name| of residue |new_residue| to |tail_name| of
    // residue |residue|. The index of the appended residue within the
    // structure is returned.

    // Attach residue's head atom to the structure's tail atom.
    int attach(const Structure *structure);
    int attach(const Residue *residue);
    int attach(const std::string& code);

    int attach(const Structure *new_structure,
               int head_residue, const std::string& head_name,
               int tail_residue, const std::string& tail_name);

    // Attach the atom in the residue with the given head name to the atom in
    // the given target residue of the structure with the specified tail name.
    int attach(const Residue *residue, const std::string& head_name,
               int target_residue, const std::string& tail_name);

    // Attach the atom in the residue with the given head index to the
    // atom in the structure with the given tail index.
    int attach(const Structure *residue, int head, int tail);

    // Attach the head atom in the residue to the atom with the given tail
    // name in the specified target residue of the structure.
    int attach_from_head(const Structure *structure, int target_residue,
                         const std::string& tail_name);

    // Attach the head atom in the residue to the given tail atom index of
    // the structure.
    int attach_from_head(const Structure *structure, int tail_atom);

    int attach_from_head(const std::string& code, const std::string& tail_atom);

    int attach_to_tail(const Structure *structure, int head_residue,
                       const std::string& head_name);

    // Attach the atom in the residue with the given name to the structure's
    // tail atom.
    int attach_to_tail(const Residue *residue, const std::string& head_name);

    // Attach the atom in the residue with the given head index to the
    // structure's tail atom.
    int attach_to_tail(const Structure *residue, int head_index);


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
    // minimization file and return a representation of the results of the
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

    // If the specified atom is bonded to another residues, these functions
    // return the index of the residue. Otherwise, they return -1;
    std::vector<int> *get_adjacent_residues_by_atom(int atom_index) const;
    std::vector<int> *get_adjacent_residues_by_atom(
            int residue_index,
            const std::string& atom_name) const;


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

    bool head_not_set() const { return head_ == -1; }
    bool tail_not_set() const { return tail_ == -1; }

    // This may not be the best place for this, but right now it is here and
    // virtual so that box information can be retrieved from structures with
    // static type Structure and dynamic type LibraryFileStructure or
    // SolvatedStructure.
    virtual const Box *box() const { return NULL; }

    //
    // Mutators
    //
    void set_bonds(const Graph *bonds);

    void set_head(int head);
    void set_tail(int tail);

    void set_head(int residue, const std::string& atom);
    void set_tail(int residue, const std::string& atom);

    void unset_head() { head_ = -1; }
    void unset_tail() { tail_ = -1; }

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
    AmberTopFile *build_amber_top_file(const ParameterSet& parm_set) const;

    // Return the AMBER topology file that represents this structure, using
    // the parameter set in the default environment.
    virtual AmberTopFile *build_amber_top_file() const;

    // Write the AMBER topology file to a file with the given name, using
    // the given parameter set.
    void print_amber_top_file(const std::string& file_name,
                              const ParameterSet& parm_set) const;

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


// This functor attaches a copy of a structure to another structure. A bond is
// formed between the specified head atom of the attached structure and the
// specified tail atom of the structure being attached to. The new bond's
// length as well as the angles and dihedrals involved in the connection are
// set appropriately using the parameter set of the default environment, if
// possible.
// TODO: This needs to be smarter about how it sets the angles and dihedrals,
// especially when the parameter set doesn't have this information, but I'm
// not sure what to do about it.
class StructureAttach {
  public:
    int operator()(Structure& structure, const Structure *new_structure,
                   int head, int tail) const;

  private:
    struct Impl;
};



// This design is questionable..
// structure.set_dihedral(DihedralAtoms(structure).atom1(..)..., value);
class DihedralAtoms {
  public:
    explicit DihedralAtoms(const Structure& structure)
            : structure_(structure), atom1_(-1), atom2_(-1), atom3_(-1),
              atom4_(-1) {}

    DihedralAtoms& atom1(int residue_index, const std::string& atom_name) {
        atom1_ = helper(residue_index, atom_name);
        return *this;
    }
    DihedralAtoms& atom2(int residue_index, const std::string& atom_name) {
        atom2_ = helper(residue_index, atom_name);
        return *this;
    }
    DihedralAtoms& atom3(int residue_index, const std::string& atom_name) {
        atom3_ = helper(residue_index, atom_name);
        return *this;
    }
    DihedralAtoms& atom4(int residue_index, const std::string& atom_name) {
        atom4_ = helper(residue_index, atom_name);
        return *this;
    }

    int atom1() const { return atom1_; }
    int atom2() const { return atom2_; }
    int atom3() const { return atom3_; }
    int atom4() const { return atom4_; }

 private:
    int helper(int residue_index, const std::string& atom_name) {
        return structure_.get_atom_index(residue_index, atom_name);
    }

    const Structure& structure_;
    int atom1_;
    int atom2_;
    int atom3_;
    int atom4_;
};

}  // namespace gmml

#include "structure-inl.h"

#endif  // GMML_INTERNAL_STRUCTURE_H_
