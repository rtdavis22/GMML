// This program shows how to load a pdb file with residue mappings.
// The global residue mappings are similar to what would be found in a LEaP
// configuration file. See include/gmml/internal/pdb_file_structure.h for
// more information.

#include <iostream>

#include "gmml/gmml.h"

using std::cout;
using std::endl;

using namespace gmml;

// The head mappings are used to map the residues in the pdb that are
// immediately after a TER card.
const char *kHeadMap[] = {
  "ALA", "NALA",
  "ARG", "NARG",
  "ASN", "NASN",
  "ASP", "NASP",
  "CYS", "NCYS",
  "CYX", "NCYX",
  "GLN", "NGLN",
  "GLU", "NGLU",
  "GLY", "NGLY",
  "HID", "NHID",
  "HIE", "NHIE",
  "HIP", "NHIP",
  "ILE", "NILE",
  "LEU", "NLEU",
  "LYS", "NLYS",
  "MET", "NMET",
  "PHE", "NPHE",
  "PRO", "NPRO",
  "SER", "NSER",
  "THR", "NTHR",
  "TRP", "NTRP",
  "TYR", "NTYR",
  "VAL", "NVAL",
  "HIS", "NHIS",
  "GUA", "DG5",
  "ADE", "DA5",
  "CYT", "DC5",
  "THY", "DT5",
  "G", "RG5",
  "A", "RA5",
  "C", "RC5",
  "U", "RU5",
  "DG", "DG5",
  "DA", "DA5",
  "DC", "DC5",
  "DT", "DT5"
};

// The tail mappings are used to map residues in the pdb that are immediately
// before a TER card.
const char *kTailMap[] = {
  "ALA", "CALA",
  "ARG", "CARG",
  "ASN", "CASN",
  "ASP", "CASP",
  "CYS", "CCYS",
  "CYX", "CCYX",
  "GLN", "CGLN",
  "GLU", "CGLU",
  "GLY", "CGLY",
  "HID", "CHID",
  "HIE", "CHIE",
  "HIP", "CHIP",
  "ILE", "CILE",
  "LEU", "CLEU",
  "LYS", "CLYS",
  "MET", "CMET",
  "PHE", "CPHE",
  "PRO", "CPRO",
  "SER", "CSER",
  "THR", "CTHR",
  "TRP", "CTRP",
  "TYR", "CTYR",
  "VAL", "CVAL",
  "HIS", "CHIS",
  "GUA", "DG3",
  "ADE", "DA3",
  "CYT", "DC3",
  "THY", "DT3",
  "G", "RG3",
  "A", "RA3",
  "C", "RC3",
  "U", "RU3",
  "DG", "DG3",
  "DA", "DA3",
  "DC", "DC3",
  "DT", "DT3"
};

// These are generic residue mappings.
const char *kResidueMap[] = {
  "HIS", "HIE",
};

#ifndef ARRAY_SIZE
#define ARRAY_SIZE(x) (sizeof(x)/sizeof(x[0]))
#endif

int main() {
    add_path("dat/");
    load_parameter_file("parm99.dat");
    load_parameter_file("Glycam_06h.dat");
    load_prep_file("Glycam_06.prep");
    load_prep_file("HOH.prep");
    load_library_file("all_amino94.lib");
    load_library_file("all_aminoct94.lib");
    load_library_file("all_aminont94.lib");

    PdbMappingInfo mapping_info;

    for (int i = 0; i < ARRAY_SIZE(kHeadMap); i += 2) {
        add_head_mapping(kHeadMap[i], kHeadMap[i + 1]);
    }

    for (int i = 0; i < ARRAY_SIZE(kTailMap); i += 2) {
        add_tail_mapping(kTailMap[i], kTailMap[i + 1]);
    }

    for (int i = 0; i < ARRAY_SIZE(kResidueMap); i += 2) {
        add_residue_mapping(kResidueMap[i], kResidueMap[i + 1]);
    }

    PdbFile pdb(File("1RVZ_New.pdb"));

    // PdbStructureBuilder allows us to specify custom mappings for this
    // particular pdb file.
    PdbStructureBuilder builder(pdb);

    // This mapping overrides the global mapping from HIS to HIE.
    builder.add_mapping("HIS", "HIP");

    // This mapping overrides the previous mapping for the residue with
    // chain_id 'A' and res_num 12 in the pdb file.
    builder.add_mapping(PdbResidueId('A', 12), "HIE");

    // These residues in the pdb file are disulfide bonded so we want to map
    // them both to CYX.
    builder.add_mapping(PdbResidueId('A', 8), "CYX");
    builder.add_mapping(PdbResidueId('B', 637), "CYX");

    PdbFileStructure *structure = builder.build();

    // We can search for a particular residue in the structure using its
    // pdb identifier:
    int res_index = structure->map_residue(PdbResidueId('A', 12));
    const Residue *residue = structure->residues(res_index);
    cout << "Atoms in " << residue->name() << ":" << endl;
    for (int i = 0; i < residue->size(); i++) {
        cout << residue->atoms(i)->name() << endl;
    }

    structure->print_amber_top_file("protein.top");
    structure->print_coordinate_file("protein.rst");

    delete structure;

    return 0;
}
