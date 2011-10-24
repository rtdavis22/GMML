#ifndef GENERIC_TYPE_H
#define GENERIC_TYPE_H

#include <iosfwd>
#include <string>

#include "utilities.h"

namespace gmml
{

//This class is only used by AmberTopFile and is probably not a great idea,
//in retrospect
class GenericType {
  public:
    enum Type { kDouble, kInt, kString };

    GenericType() : type_(kDouble), dval_(kNotSet) {}
    GenericType(int n) : type_(kInt), dval_(n) {}
    GenericType(double n) : type_(kDouble), dval_(n) {}
    GenericType(const std::string& s) : type_(kString), sval_(s) {}

    GenericType& operator=(const GenericType&);
    GenericType& operator=(double);
    GenericType& operator=(const std::string&);

    operator std::string() const;
    operator double() const;

    friend std::ostream& operator<<(std::ostream&, GenericType&);

  private:
    Type type_;
    double dval_;
    std::string sval_;
};

} //namespace gmml

#endif
