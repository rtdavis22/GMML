// Author: Robert Davis

#ifndef SEQUENCE_PARSER_H
#define SEQUENCE_PARSER_H

#include <exception>
#include <map>
#include <string>

#include "residue_classification.h"
#include "tree.hh"
#include "utilities.h"

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
    std::map<int, char> derivatives;
    bool is_terminal;
};

class SequenceParser {
  public:
    SequenceParser() {}

    virtual tree<ParsedResidue*> *parse(const std::string& sequence) const = 0;
};

}  // namespace gmml

#endif  // SEQUENCE_PARSER_H
