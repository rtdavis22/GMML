// Author: Robert Davis

#ifndef GMML_INTERNAL_GLYCAM_PARSER_H_
#define GMML_INTERNAL_GLYCAM_PARSER_H_

#include <string>

#include "gmml/internal/array_tree.h"
#include "gmml/internal/sequence_parser.h"
#include "gmml/internal/stubs/common.h"

namespace gmml {

// I think both of the public functions should go in SequenceParser, and
// get_parse_info should be purely virtual in SequenceParser.
class GlycamParser : public SequenceParser {
  public:
    GlycamParser() {}

    // This returns a "first-child, next-sibling" n-ary tree.
    virtual tree<ParsedResidue*> *parse(const std::string& sequence) const;

    // This representation is often more convenient because it preserves the
    // order of the residues in the sequence.
    ArrayTree<ParsedResidue*> *get_array_tree(
            const std::string& sequence) const;

  private:
    enum TokenType { kTokenLeft, kTokenRight, kTokenResidue };
    struct ParseInfo {
        ParseInfo(std::vector<ParsedResidue*> *residues,
                  std::vector<TokenType> *tokens) : residues(residues),
                                                    tokens(tokens) {}
        std::vector<ParsedResidue*> *residues;
        std::vector<TokenType> *tokens;
    };

    ParseInfo *get_parse_info(const std::string& sequence) const;
    
    ParsedResidue *parse_residue(const std::string& residue_string) const;

    DISALLOW_COPY_AND_ASSIGN(GlycamParser);
};

}  // namespace gmml

#endif  // GMML_INTERNAL_GLYCAM_PARSER_H
