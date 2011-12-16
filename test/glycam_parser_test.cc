// Author: Robert Davis

#include "gmml/gmml.h"
#include "gtest/gtest.h"

using gmml::ArrayTree;
using gmml::GlycamParser;
using gmml::ParsedResidue;

TEST(GlycamParser, ArrayTreeEmptySequence) {
    GlycamParser parser;
    ArrayTree<ParsedResidue*> *tree = parser.get_array_tree("");
    EXPECT_EQ(tree->size(), 0);
}

TEST(GlycamParser, ParseEmptySequence) {
    GlycamParser parser;
    gmml::tree<ParsedResidue*> *tree = parser.parse("");
    EXPECT_EQ(tree->size(), 0);
}

TEST(GlycamParser, ArrayTreeAglycon) {
    GlycamParser parser;
    ArrayTree<ParsedResidue*> *tree = parser.get_array_tree("OH");
    EXPECT_EQ(tree->size(), 1);
    ParsedResidue *residue = tree->begin()->first;
    EXPECT_EQ(residue->is_terminal, true);
    EXPECT_EQ(residue->name, "OH");
}

TEST(GlycamParser, ParseAglycon) {
    GlycamParser parser;
    gmml::tree<ParsedResidue*> *tree = parser.parse("OH");
    EXPECT_EQ(tree->size(), 1);
    ParsedResidue *residue = *tree->begin();
    EXPECT_EQ(residue->is_terminal, true);
    EXPECT_EQ(residue->name, "OH");
}
