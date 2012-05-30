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

#ifndef GMML_INTERNAL_SEQUENCE_PARSER_H_
#define GMML_INTERNAL_SEQUENCE_PARSER_H_

#include <exception>
#include <map>
#include <string>

#include "gmml/internal/residue_classification.h"
#include "gmml/internal/stubs/common.h"
#include "gmml/internal/tree.h"

namespace gmml {

struct ParseException : public std::exception {
  public:
    explicit ParseException(const std::string& what) { what_ = what; }

    virtual const char *what() const throw() { return what_.c_str(); }

    virtual ~ParseException() throw() {}

  private:
    std::string what_;
};

// If is_terminal is true, the only field that is valid is name.
// The ring type is embedded in the name of the residue because it may not be
// easy to extract it from the residue name. It should probably be included.
struct ParsedResidue {
    ParsedResidue() : is_terminal(false) {}

    ResidueClassification::Isomer isomer;
    ResidueClassification::Configuration configuration;
    std::string name;
    int anomeric_carbon;
    int oxygen_position;
    std::map<int, std::string> derivatives;
    bool is_terminal;
};

// This is an abstract base class for parsers that parse a sequence of
// residues.
class SequenceParser {
  public:
    SequenceParser() {}

    virtual tree<ParsedResidue*> *parse(const std::string& sequence) const = 0;
};

}  // namespace gmml

#endif  // GMML_INTERNAL_SEQUENCE_PARSER_H_
