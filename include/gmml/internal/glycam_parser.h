// Author: Robert Davis

#ifndef GLYCAM_PARSER_H
#define GLYCAM_PARSER_H

#include <string>

#include "array_tree.h"
#include "sequence_parser.h"

namespace gmml {

class GlycamParser : public SequenceParser {
  public:
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
};

}  // namespace gmml

#endif  // GLYCAM_PARSER_H
