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
    FastaSequence *sequence = create_empty_sequence();
    int residue_count = structure.residue_count();
    for (int i = 0; i < residue_indices.size(); i++) {
        int index = residue_indices[i];
        if (index < 0 || index >= residue_count) {
            delete sequence;
            throw std::invalid_argument("Invalid index " + to_string(index));
        }
        sequence->add_amino_acid_code(structure.residues(index)->name());
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
