// Author: Robert Davis

#include "gmml/gmml.h"
#include "gtest/gtest.h"

using gmml::Atom;
using gmml::File;
using gmml::PdbFile;
using gmml::PdbFileStructure;
using gmml::PdbResidueId;
using gmml::PdbStructureBuilder;
using gmml::Residue;

class PdbStructureTest : public ::testing::Test {
  protected:
    virtual void SetUp() {
        PdbFile pdb_file(File("dat/1RVZ_New.pdb"));
        PdbStructureBuilder builder(pdb_file);
        structure = builder.build();
    }

    virtual void TearDown() {
        delete structure;
    }

    PdbFileStructure *structure;
};

TEST_F(PdbStructureTest, NotNull) {
    EXPECT_TRUE(structure != NULL);
}

TEST_F(PdbStructureTest, MapResidueValid) {
    int residue_index = structure->map_residue('A', 7);
    ASSERT_GE(residue_index, 0);
    ASSERT_LT(residue_index, structure->residue_count());
    Residue *residue = structure->residues(residue_index);
    EXPECT_EQ("ILE", residue->name());
}

TEST_F(PdbStructureTest, MapResidueEmptyICode) {
    EXPECT_EQ(structure->map_residue('A', 7),
              structure->map_residue('A', 7, ' '));
}

TEST_F(PdbStructureTest, MapResidueInvalid) {
    EXPECT_EQ(-1, structure->map_residue('Z', 3));
}

TEST_F(PdbStructureTest, MapAtomValid) {
    int atom_index = structure->map_atom(53);
    ASSERT_GE(atom_index, 0);
    ASSERT_LT(atom_index, structure->size());
    const Atom *atom = structure->atoms(atom_index);
    EXPECT_EQ("OH", atom->name());
    EXPECT_EQ(gmml::kElementO, atom->element());
    EXPECT_EQ(16.653, atom->coordinate().x);
    EXPECT_EQ(-25.833, atom->coordinate().y);
    EXPECT_EQ(15.868, atom->coordinate().z);
}

TEST(PdbStructureBuilder, IsToBeRemoved1) {
    PdbFile pdb(File("dat/1RVZ_New.pdb"));
    PdbStructureBuilder builder(pdb);
    PdbResidueId id('A', 123, ' ');
    builder.add_residue_to_remove(&id);
    EXPECT_TRUE(builder.is_to_be_removed(&id));
}

TEST(PdbStructureBuilder, IsToBeRemoved2) {
    PdbFile pdb(File("dat/1RVZ_New.pdb"));
    PdbStructureBuilder builder(pdb);
    PdbResidueId id('A', 123, ' ');
    EXPECT_FALSE(builder.is_to_be_removed(&id));
}
