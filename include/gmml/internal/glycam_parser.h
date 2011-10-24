#ifndef GLYCAM_PARSER_H
#define GLYCAM_PARSER_H

#include <string>

#include "sequence_parser.h"

namespace gmml {

class GlycamParser : public SequenceParser {
  public:
    virtual tree<ParsedResidue*> *parse(const std::string& sequence) const;

  private:
    enum TokenType { kTokenLeft, kTokenRight, kTokenResidue };

    ParsedResidue *parse_residue(const std::string& residue_string) const;
};

}  // namespace gmml

#endif  // GLYCAM_PARSER_H
