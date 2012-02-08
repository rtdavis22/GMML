#ifndef GMML_INTERNAL_PROTEINS_H_
#define GMML_INTERNAL_PROTEINS_H_

#include <memory>
#include <string>
#include <vector>

namespace gmml {

class Structure;

class AminoAcidCodes {
  public:
    AminoAcidCodes();
    ~AminoAcidCodes();

    bool is_amino_acid(const std::string& code);
    std::string get_fasta_code(const std::string& code);

  private:
    struct Impl;
    std::auto_ptr<Impl> impl_;
};

std::vector<std::vector<int>*> *find_proteins(const Structure& structure);

std::vector<int> *get_asns_with_sequon(const Structure& structure);

}  // namespace gmml

#endif  // GMML_INTERNAL_PROTEINS_H_
