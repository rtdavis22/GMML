#include "gmml/gmml.h"
#include "gtest/gtest.h"

using gmml::FastaSequence;

TEST(FastaSequence, CreateEmptySequence) {
    FastaSequence *sequence = FastaSequence::create_empty_sequence();
    EXPECT_EQ("", sequence->sequence());
}
