// Author: Robert Davis

#ifndef GMML_INTERNAL_PROTEINS_H_
#define GMML_INTERNAL_PROTEINS_H_

#include <memory>
#include <string>
#include <vector>

#include "gmml/internal/stubs/common.h"

namespace gmml {

class Structure;

class AminoAcidCodeSet {
  public:
    AminoAcidCodeSet();
    ~AminoAcidCodeSet();

    bool lookup(const std::string& code) const;
    char get_fasta_letter(const std::string& code) const;

  private:
    struct Impl;
    std::auto_ptr<Impl> impl_;

    DISALLOW_COPY_AND_ASSIGN(AminoAcidCodeSet);
};

std::vector<std::vector<int>*> *find_proteins(const Structure& structure);

std::vector<int> *get_asns_with_sequon(const Structure& structure);

void load_amino_acid_mappings();

}  // namespace gmml

#endif  // GMML_INTERNAL_PROTEINS_H_
