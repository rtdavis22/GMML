// Author: Robert Davis

#include "gmml/gmml.h"
#include "gtest/gtest.h"

using gmml::PdbFileStructure;
using gmml::Residue;

class PdbStructureTest : public ::testing::Test {
  protected:
    virtual void SetUp() {
        structure = PdbFileStructure::build("dat/1RVZ_New.pdb");
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

TEST_F(PdbStructureTest, MapAtomValue) {
    
}
