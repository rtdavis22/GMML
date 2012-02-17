// Author: Robert Davis

#include "gmml/internal/fasta_sequence.h"

#include <stdexcept>

#include "gmml/internal/proteins.h"
#include "gmml/internal/residue.h"
#include "gmml/internal/structure.h"
#include "gmml/internal/stubs/utils.h"

using std::string;
using std::vector;

namespace gmml {

struct FastaSequence::Impl {
    explicit Impl(const string& sequence) : sequence(sequence) {}

    // This is lazily initialized and not freed.
    static const AminoAcidCodeSet *amino_acid_codes;

    string sequence;
};

const AminoAcidCodeSet *FastaSequence::Impl::amino_acid_codes = NULL;

FastaSequence::~FastaSequence() {
}

FastaSequence *FastaSequence::clone() const {
    return FastaSequence::create(impl_->sequence);
}

FastaSequence *FastaSequence::create(const string& sequence) {
    return new FastaSequence(sequence);
}

FastaSequence *FastaSequence::create(const Structure& structure,
                                     const vector<int>& residue_indices) {
    vector<string> codes;
    for (int i = 0; i < structure.residue_count(); i++) {
        int index = residue_indices[i];
        codes.push_back(structure.residues(index)->name());
    }
    return create(codes);
}

FastaSequence *FastaSequence::create(const vector<string>& codes) {
    FastaSequence *sequence = create_empty_sequence();
    for (int i = 0; i < codes.size(); i++) {
        sequence->add_amino_acid_code(codes[i]);
    }
    return sequence;
}

FastaSequence *FastaSequence::create_empty_sequence() {
    return create("");
}

void FastaSequence::add_amino_acid_code(const string& code) {
    if (impl_->amino_acid_codes == NULL)
        impl_->amino_acid_codes = new AminoAcidCodeSet;
    add_letter(impl_->amino_acid_codes->get_fasta_letter(code));
}

void FastaSequence::add_letter(char letter) {
    impl_->sequence += letter;
}

string FastaSequence::sequence() const {
    return impl_->sequence;
}

FastaSequence::FastaSequence(const string& sequence) :
        impl_(new Impl(sequence)) {
}

}  // namespace gmml
