// Copyright (c) 2012 The University of Georgia
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.

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
    GlycamParser() : parse_derivatives_(true) {}

    // This returns a "first-child, next-sibling" n-ary tree.
    virtual tree<ParsedResidue*> *parse(const std::string& sequence) const;

    // This representation is often more convenient because it preserves the
    // order of the residues in the sequence.
    ArrayTree<ParsedResidue*> *get_array_tree(
            const std::string& sequence) const;

    void dont_parse_derivatives() { parse_derivatives_ = false; }
    void parse_derivatives() { parse_derivatives_ = true; }

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

    bool parse_derivatives_;

    DISALLOW_COPY_AND_ASSIGN(GlycamParser);
};

}  // namespace gmml

#endif  // GMML_INTERNAL_GLYCAM_PARSER_H
