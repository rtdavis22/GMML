// Author: Robert Davis

#ifndef GMML_INTERNAL_FASTA_SEQUENCE_H_
#define GMML_INTERNAL_FASTA_SEQUENCE_H_

#include <memory>
#include <string>
#include <vector>

#include "gmml/internal/stubs/common.h"

namespace gmml {

class Structure;

class FastaSequence {
  public:
    virtual ~FastaSequence();

    virtual FastaSequence *clone() const;

    static FastaSequence *create(const std::string& sequence);
    static FastaSequence *create(const Structure& structure,
                                 const std::vector<int>& residue_indices);
    static FastaSequence *create(const std::vector<std::string>& codes);
    static FastaSequence *create_empty_sequence();

    void add_amino_acid_code(const std::string& code);

    void add_letter(char letter);

    std::string sequence() const;

  private:
    explicit FastaSequence(const std::string& sequence);

    struct Impl;
    std::auto_ptr<Impl> impl_;

    DISALLOW_COPY_AND_ASSIGN(FastaSequence);
};

}  // namespace gmml

#endif  // GMML_INTERNAL_FASTA_SEQUENCE_H_
