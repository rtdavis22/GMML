// Author: Robert Davis

#ifndef GLYCAM_CODE_SET_H_
#define GLYCAM_CODE_SET_H_

#include <bitset>
#include <map>
#include <string>
#include <vector>

#include "array_tree.h"
#include "residue_classification.h"
#include "tree.hh" // see if I can get out of including this
#include "utilities.h"

namespace gmml {

class ParsedResidue;
class Residue;
class Structure;
class TreeResidue;

class GlycamCodeSet {
  public:
    // The default constructor populates the code maps. One maps from a residue
    // name to the corresponding letter and the other is letter to residue name.
    GlycamCodeSet();

    // These parsed residues are created from GlycamParser. It is important
    // to note that the name member of the parsed residue does include the
    // ring type (GalpNAc, for example). Open valences is an unordered list
    // of oxygen positions that are open valences.
    std::string get_code(const ParsedResidue& parsed_residue,
                         const std::vector<int>& open_valences) const;

    // The residue name here is the residue name as it would appear in the
    // glycam sequence without the ring type. So the residue name of
    // GalpNAc is GalNAc.
    std::string get_code(const std::string& residue_name,
                         ResidueClassification::Isomer,
                         ResidueClassification::RingType,
                         ResidueClassification::Configuration,
                         const std::vector<int>& open_valences) const;

    // This residue name does include the ring type. It is the the name as it
    // would appear in the glycam sequence.
    std::string get_code(const std::string& residue_name,
                         ResidueClassification::Isomer,
                         ResidueClassification::Configuration,
                         const std::vector<int>& open_valences) const;

    // This returns the code for the name of an aglycon.
    std::string get_terminal_code(const std::string& terminal_name) const;

    std::string get_derivative_code(char derivative_letter) const;

    std::string get_name_from_code(const std::string& code) const;

    tree<TreeResidue*> *build_residue_tree(
            tree<ParsedResidue*> *parsed_tree) const;

    ArrayTree<TreeResidue*> *build_array_tree(
            ArrayTree<ParsedResidue*> *parsed_tree) const;

  private:
    // Because of the way the code mapping was designed, the the second letter
    // can actually be two letters. Therefore all of these return strings
    // instead of chars to keep it consistent. Also it is possible that the
    // others may be updated to return more than one of no letters.
    //
    // It's easier to work with bitsets internally, so they are used here,
    // but they should not be a part of the public interface. Also, it should
    // be changed to a dynamic bitset when C++0x is more widespread.
    std::string get_first_letter(const std::bitset<10>& open_valences) const;

    // This residue name does not include the ring type.
    std::string get_second_letter(const std::string& residue_name,
                                  ResidueClassification::Isomer) const;

    std::string get_third_letter(ResidueClassification::Configuration,
                                 ResidueClassification::RingType) const;

    // This builds a TreeResidue from a ParsedResidue. A tree iterator is
    // used because the name of the residue depends on the children of
    // ParsedResidue.
    TreeResidue *build_tree_residue(tree<ParsedResidue*>::iterator) const;

    TreeResidue *build_tree_residue(const ParsedResidue& parsed_residue,
                                    const std::vector<int>& valences,
                                    const std::string& parent_name) const;

    TreeResidue *get_derivative_tree_residue(char derivative, int pos) const;

    std::string get_oxygen_name(const std::string& residue_name,
                                int position) const;

    std::map<std::string, std::string> name_to_letter_;
    std::map<std::string, std::string> letter_to_name_;

    DISALLOW_COPY_AND_ASSIGN(GlycamCodeSet);
};

// This functor is a wrapper around the default attachment functor. It accounts
// for changes that need to be made in the charges because of derivatives.
struct GlycamAttach {
    int operator()(Structure& structure, Residue *residue,
                   const std::string& new_atom_name, int residue_index,
                   const std::string& atom_name) const;
};

// The residues in structures typically form a tree. These functions build the
// Structure from the tree. The second function visits and attaches the residues
// in pre-order, so the root residue of the tree is the residue at index 0
// of the structure.
Structure *glycam_build(ArrayTree<TreeResidue*> *residue_tree);
Structure *glycam_build(tree<TreeResidue*> *residue_tree);

// These functions are invoked by the next two functions but may be useful
// by themselves.
Structure *glycam_build(ArrayTree<ParsedResidue*> *parsed_tree);
Structure *glycam_build(tree<ParsedResidue*> *parsed_tree);

// The output of this function is the Structure representation of the
// sequence in GLYCAM condensed notation. The residues in the Structure are
// in the reverse order compared to their position in the sequence (the residue
// at index 0 in the structure is the aglycon). Be sure to load the GLYCAM
// prep file before calling this function.
Structure *glycam_build(const std::string& sequence);

// This builds the sequence using an ArrayTree rather than a tree, as in the
// previous function. The performance difference is probably negligable.
// The resulting structures should be equivalent.
Structure *glycam_build_with_array_tree(const std::string& sequence);

}  // namespace gmml

#endif  // GLYCAM_CODE_SET_H_
