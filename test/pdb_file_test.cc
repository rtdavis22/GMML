// Author: Kyle Forrester

#include "gmml/gmml.h"
#include "gtest/gtest.h"
#include <string>
#include <sstream>

using namespace gmml;

TEST(PdbSeqresCardBuilder, AccessMethods) {
    PdbSeqresCardBuilder *builder = new PdbSeqresCardBuilder('N');
    EXPECT_EQ('N', builder->chain_id());
    builder->add_residue("PHE");
    builder->add_residue("ILE");
    builder->add_residue("CYS");
    builder->add_residue("U");
    builder->add_residue("C");
    builder->add_residue("LEU");
    builder->add_residue("C");
    EXPECT_EQ("PHE", builder->at(0));
    EXPECT_EQ("ILE", builder->at(1));
    EXPECT_EQ("CYS", builder->at(2));
    EXPECT_EQ("U", builder->at(3));
    EXPECT_EQ("C", builder->at(4));
    EXPECT_EQ("LEU", builder->at(5));
    EXPECT_EQ("C", builder->at(6));
    EXPECT_EQ(7, builder->size());
}

TEST(PdbSeqresCard, CreateCards) {
    PdbSeqresCardBuilder *builder = new PdbSeqresCardBuilder('V');
    builder->add_residue("PHE");
    builder->add_residue("VAL");
    builder->add_residue("ASN");
    builder->add_residue("GLN");
    builder->add_residue("HIS");
    builder->add_residue("LEU");
    builder->add_residue("CYS");
    builder->add_residue("GLY");
    builder->add_residue("SER");
    builder->add_residue("HIS");
    builder->add_residue("LEU");
    builder->add_residue("VAL");
    builder->add_residue("GLU");
    builder->add_residue("ALA");
    builder->add_residue("LEU");
    builder->add_residue("TYR");
    builder->add_residue("LEU");
    builder->add_residue("VAL");
    builder->add_residue("CYS");
    builder->add_residue("GLY");
    builder->add_residue("GLU");
    builder->add_residue("ARG");
    builder->add_residue("GLY");
    builder->add_residue("PHE");
    builder->add_residue("PHE");
    builder->add_residue("TYR");
    builder->add_residue("THR");
    builder->add_residue("PRO");
    builder->add_residue("LYS");
    builder->add_residue("ALA");
    std::vector<PdbSeqresCard*> builder_vector = builder->build();
    PdbSeqresCard *current_card = builder_vector[0];
    EXPECT_EQ('V', current_card->chain_id());
    EXPECT_EQ(1, current_card->serial_number());
    EXPECT_EQ(13, current_card->residues_size());
    EXPECT_EQ(30, current_card->number_of_chain_residues());
    EXPECT_EQ("PHE", current_card->at(0));
    EXPECT_EQ("VAL", current_card->at(1));
    EXPECT_EQ("ASN", current_card->at(2));
    EXPECT_EQ("GLN", current_card->at(3));
    EXPECT_EQ("HIS", current_card->at(4));
    EXPECT_EQ("LEU", current_card->at(5));
    EXPECT_EQ("CYS", current_card->at(6));
    EXPECT_EQ("GLY", current_card->at(7));
    EXPECT_EQ("SER", current_card->at(8));
    EXPECT_EQ("HIS", current_card->at(9));
    EXPECT_EQ("LEU", current_card->at(10));
    EXPECT_EQ("VAL", current_card->at(11));
    EXPECT_EQ("GLU", current_card->at(12));
    current_card = builder_vector[1];
    EXPECT_EQ('V', current_card->chain_id());
    EXPECT_EQ(2, current_card->serial_number());
    EXPECT_EQ(13, current_card->residues_size());
    EXPECT_EQ(30, current_card->number_of_chain_residues());
    EXPECT_EQ("ALA", current_card->at(0));
    EXPECT_EQ("LEU", current_card->at(1));
    EXPECT_EQ("TYR", current_card->at(2));
    EXPECT_EQ("LEU", current_card->at(3));
    EXPECT_EQ("VAL", current_card->at(4));
    EXPECT_EQ("CYS", current_card->at(5));
    EXPECT_EQ("GLY", current_card->at(6));
    EXPECT_EQ("GLU", current_card->at(7));
    EXPECT_EQ("ARG", current_card->at(8));
    EXPECT_EQ("GLY", current_card->at(9));
    EXPECT_EQ("PHE", current_card->at(10));
    EXPECT_EQ("PHE", current_card->at(11));
    EXPECT_EQ("TYR", current_card->at(12));
    current_card = builder_vector[2];
    EXPECT_EQ('V', current_card->chain_id());
    EXPECT_EQ(1, current_card->serial_number());
    EXPECT_EQ(4, current_card->residues_size());
    EXPECT_EQ(30, current_card->number_of_chain_residues());
    EXPECT_EQ("THR", current_card->at(0));
    EXPECT_EQ("PRO", current_card->at(1));
    EXPECT_EQ("LYS", current_card->at(2));
    EXPECT_EQ("ALA", current_card->at(3));
}

TEST(PdbSeqresCard, Read) {
    PdbLine *pdb_line = new PdbLine("SEQRES   2 A   21  TYR GLN LEU GLU ASN TYR 
            CYS ASN");
    PdbSeqresCard *current_card = new PdbSeqresCard(*pdb_line);
    EXPECT_EQ('A', current_card->chain_id());
    EXPECT_EQ(2, current_card->serial_number());
    EXPECT_EQ(8, current_card->residues_size());
    EXPECT_EQ(21, current_card->number_of_chain_residues());
    EXPECT_EQ("TYR", current_card->at(0));
    EXPECT_EQ("GLN", current_card->at(1));
    EXPECT_EQ("LEU", current_card->at(2));
    EXPECT_EQ("GLU", current_card->at(3));
    EXPECT_EQ("ASN", current_card->at(4));
    EXPECT_EQ("TYR", current_card->at(5));
    EXPECT_EQ("CYS", current_card->at(6));
    EXPECT_EQ("ASN", current_card->at(7));
}

TEST(PdbSeqresCard, Write) {
    std::string str = "SEQRES   2 X   39    U   A   G   C   G   G   C   G   U   G   G   A   A";
    PdbLine *pdb_line = new PdbLine(str);
    PdbSeqresCard *seqres_card = new PdbSeqresCard(*pdb_line);
    std::ostringstream stream;
    seqres_card->write(stream);
    ASSERT_EQ(str, stream.str());
}

TEST(PdbModresCard, ConstructorWrite) {
    PdbModresCard *modres_card = new PdbModresCard("2R0L", "ASN", 74, ' ', "ASN");
    modres_card->set_comment("GLYCOSYLATION SITE");
    std::string str = "MODRES 2R0L ASN A   74  ASN  GLYCOSYLATION SITE";
    std::ostringstream stream;
    modres_card->write(stream);
    ASSERT_EQ(str, stream.str());
}

TEST(PdbModresCard, Read) {
    std::string str = "MODRES 1IL2 1MG D 1937    G  1N-METHYLGUANOSINE-5'-MONOPHOSPHATE";
    PdbLine *pdb_line = new PdbLine(str);
    PdbModresCard *modres_card = new PdbModresCard(*pdb_line);
    std::ostringstream stream;
    modres_card->write(stream);
    ASSERT_EQ(str, stream.str());
}

TEST(PdbSsbondCard, ConstructorWrite) {
    std::string str = "SSBOND   1 CYS A    6    CYS A  127                          1555   1555  2.03";
    PdbSsbondCard *ssbond_card = new PdbSsbondCard(1, "CYS", 'A', 6, ' ', "CYS", 'A', 127, ' ', 1555, 1555, 2.03);
    std::ostringstream stream;
    ssbond_card->write(stream);
    ASSERT_EQ(str, stream.str());
}

TEST(PdbSsbondCard, Read) {
    std::string str = "SSBOND   3 CYS A   64    CYS A   80                          1555   1555  2.06";
    PdbLine *pdb_line = new PdbLine(str);
    PdbSsbondCard *ssbond_card = new PdbSsbondCard(*pdb_line);
    std::ostringstream stream;
    ssbond_card->write(stream);
    EXPECT_EQ(str, stream.str());
}

TEST(PdbSiteCardBuilder, AccessMethods) {
    PdbSiteCardBuilder *builder = new PdbSiteCardBuilder("AC2");
    NamedPdbResidueId *res0 = new NamedPdbResidueId("ASN", 'A', 62, ' ');
    NamedPdbResidueId *res1 = new NamedPdbResidueId("GLY", 'A', 63, ' ');
    NamedPdbResidueId *res2 = new NamedPdbResidueId("HIS", 'A', 64, ' ');
    NamedPdbResidueId *res3 = new NamedPdbResidueId("HOH", 'A', 328, ' ');
    NamedPdbResidueId *res4 = new NamedPdbResidueId("HOH", 'A', 634, ' ');
    builder->add_residue(*res0);
    builder->add_residue(*res1);
    builder->add_residue(*res2);
    builder->add_residue(*res3);
    builder->add_residue(*res4);
    EXPECT_EQ(5, builder->size());
    EXPECT_EQ("AC2", builder->site_name);
    EXPECT_EQ(*res0, builder.at(0));
    EXPECT_EQ(*res1, builder.at(1));
    EXPECT_EQ(*res2, builder.at(2));
    EXPECT_EQ(*res3, builder.at(3));
    EXPECT_EQ(*res4, builder.at(4));
}

TEST(PdbSiteCardBuilder, CreateCards) {
    PdbSiteCardBuilder *builder = new PdbSiteCardBuilder("AC4");
    NamedPdbResidueId *res0 = new NamedPdbResidueId("HIS", 'A', 64, ' ');
    NamedPdbResidueId *res1 = new NamedPdbResidueId("HIS", 'A', 94, ' ');
    NamedPdbResidueId *res2 = new NamedPdbResidueId("HIS", 'A', 96, ' ');
    NamedPdbResidueId *res3 = new NamedPdbResidueId("HIS", 'A', 119, ' ');
    NamedPdbResidueId *res4 = new NamedPdbResidueId("LEU", 'A', 198, ' ');
    NamedPdbResidueId *res5 = new NamedPdbResidueId("THR", 'A', 199, ' ');
    NamedPdbResidueId *res6 = new NamedPdbResidueId("THR", 'A', 200, ' ');
    NamedPdbResidueId *res7 = new NamedPdbResidueId("TRP", 'A', 209, ' ');
    NamedPdbResidueId *res8 = new NamedPdbResidueId("HOH", 'A', 572, ' ');
    NamedPdbResidueId *res9 = new NamedPdbResidueId("HOH", 'A', 582, ' ');
    NamedPdbResidueId *res10 = new NamedPdbResidueId("HOH", 'A', 635, ' ');
    builder->add_residue(*res0);
    builder->add_residue(*res1);
    builder->add_residue(*res2);
    builder->add_residue(*res3);
    builder->add_residue(*res4);
    builder->add_residue(*res5);
    builder->add_residue(*res6);
    builder->add_residue(*res7);
    builder->add_residue(*res8);
    builder->add_residue(*res9);
    builder->add_residue(*res10);
    std::vector<PdbSiteCard*> site_cards = builder->build();
    EXPECT_EQ(*res0, site_cards[0]);
    EXPECT_EQ(*res1, site_cards[1]);
    EXPECT_EQ(*res2, site_cards[2]);
    EXPECT_EQ(*res3, site_cards[3]);
    EXPECT_EQ(*res4, site_cards[4]);
    EXPECT_EQ(*res5, site_cards[5]);
    EXPECT_EQ(*res6, site_cards[6]);
    EXPECT_EQ(*res7, site_cards[7]);
    EXPECT_EQ(*res8, site_cards[8]);
    EXPECT_EQ(*res9, site_cards[9]);
    EXPECT_EQ(*res10, site_cards[10]);
}

TEST(PdbSiteCard, ReadandWrite) {
    string str = "SITE     3 AC4 11 HOH A 572  HOH A 582  HOH A 635           ";
    trim(str);
    PdbLine *pdb_line = new PdbLine(str);
    PdbSiteCard *site_card = new PdbSiteCard(*pdb_line);
    std::ostringstream stream;
    site_card->write(stream);
    ASSERT_EQ(str, stream.str());
}

TEST(PdbModelCard, ConstructorWrite) {
    string str = "MODEL        1    ";
    trim(str);
    PdbModelCard *model_card = new PdbModelCard(1);
    std::ostringstream stream;
    model_card->write(stream);
    ASSERT_EQ(str, stream.str());
}

TEST(PdbModelCard, Read) {
    string str = "MODEL        5  ";
    trim(str);
    PdbLine *pdb_line = new PdbLine(str);
    PdbModelCard *model_card = new PdbModelCard(*pdb_line);
    std::ostringstream stream;
    model_card->write(stream);
    ASSERT_EQ(str, stream.str());
}

}
