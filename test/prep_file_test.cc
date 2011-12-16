// Author: Robert Davis

#include "gmml/gmml.h"
#include "gtest/gtest.h"

using gmml::PrepFile;
using gmml::PrepFileAtom;
using gmml::PrepFileResidue;
using gmml::PrepFileSet;
using gmml::Residue;

TEST(PrepFileTest, Constructor) {
    PrepFile file("dat/test.prep");
    PrepFile::const_iterator begin = file.begin();
    EXPECT_EQ(begin->first, "PeA");
    PrepFile::const_iterator end = file.end();
    --end;
    EXPECT_EQ(end->first, "PeA");
    EXPECT_EQ(begin->second, end->second);
}

TEST(PrepFileTest, PrepFileResidue) {
    PrepFile file("dat/test.prep");
    PrepFile::const_iterator it = file.begin();
    PrepFile::ResiduePtr residue = it->second;
    EXPECT_EQ(residue->name, "PeA");
    EXPECT_EQ(residue->coordinate_type, PrepFileResidue::kINT);
    EXPECT_EQ(residue->output_format, PrepFileResidue::kFormatted);
    EXPECT_EQ(residue->geometry_type, PrepFileResidue::kGeometryCorrect);
    EXPECT_EQ(residue->dummy_atom_omission, PrepFileResidue::kOmit);
    EXPECT_EQ(residue->dummy_atom_type, "DU");
    EXPECT_EQ(residue->dummy_atom_position, PrepFileResidue::kPositionBeg);
    EXPECT_EQ(residue->cutoff, -0.5820);
    EXPECT_EQ(residue->atoms.size(), 21);
    EXPECT_EQ(residue->loops.size(), 1);
    EXPECT_EQ(residue->loops[0].from, 17);
    EXPECT_EQ(residue->loops[0].to, 3);
}

TEST(PrepFileTest, FirstNonDummyAtom) {
    PrepFile file("dat/test.prep");
    PrepFile::const_iterator it = file.begin();
    PrepFile::ResiduePtr residue = it->second;
    PrepFileAtom *atom = residue->atoms[3];
    EXPECT_EQ(atom->index, 4);
    EXPECT_EQ(atom->name, "C1");
    EXPECT_EQ(atom->type, "CG");
    EXPECT_EQ(atom->topological_type, PrepFileAtom::kTopTypeM);
    EXPECT_EQ(atom->bond_index, 3);
    EXPECT_EQ(atom->angle_index, 2);
    EXPECT_EQ(atom->dihedral_index, 1);
    EXPECT_EQ(atom->bond_length, 1.4);
    EXPECT_EQ(atom->angle, 113.3);
    EXPECT_EQ(atom->dihedral, -60.0);
    EXPECT_EQ(atom->charge, 0.48);
}

TEST(PrepFileTest, LastAtom) {
    PrepFile file("dat/test.prep");
    PrepFile::const_iterator it = file.begin();
    PrepFile::ResiduePtr residue = it->second;
    PrepFileAtom *atom = residue->atoms[3];
    atom = residue->atoms[20];
    EXPECT_EQ(atom->index, 21);
    EXPECT_EQ(atom->name, "O6");
    EXPECT_EQ(atom->type, "OS");
    EXPECT_EQ(atom->topological_type, PrepFileAtom::kTopTypeM);
    EXPECT_EQ(atom->bond_index, 9);
    EXPECT_EQ(atom->angle_index, 7);
    EXPECT_EQ(atom->dihedral_index, 6);
    EXPECT_EQ(atom->bond_length, 1.413);
    EXPECT_EQ(atom->angle, 112.7);
    EXPECT_EQ(atom->dihedral, -56.7);
    EXPECT_EQ(atom->charge, -0.451);
}

TEST(PrepFileSet, Constructor) {
    PrepFileSet set;
    PrepFile::iterator begin = set.begin();
    PrepFile::iterator end = set.end();
    EXPECT_EQ(begin == end, true);
}

TEST(PrepFileSet, Exists) {
    PrepFileSet set;
    set.load("dat/test.prep");
    EXPECT_EQ(set.exists("PeA"), true);
    EXPECT_EQ(set.exists("OME"), false);
}

TEST(PrepFileSet, Brackets) {
    PrepFileSet set;
    set.load("dat/test.prep");
    const PrepFileResidue& residue = set["PeA"];
    EXPECT_EQ(residue.name, "PeA");
}

TEST(PrepFileSet, Lookup) {
    PrepFileSet set;
    set.load("dat/test.prep");
    PrepFile::ResiduePtr residue = set.lookup("PeA");
    EXPECT_NE(residue, PrepFile::ResiduePtr());
    EXPECT_EQ(residue->name, "PeA");
    EXPECT_EQ(set.lookup("OME"), PrepFile::ResiduePtr());
}

TEST(BuildPrepFileResidue, TypicalResidue) {
    PrepFileSet set;
    set.load("dat/test.prep");
    Residue *residue = gmml::build_prep_file(set["PeA"]);
    EXPECT_EQ(residue->name(), "PeA");
}
