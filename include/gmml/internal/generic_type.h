#ifndef GMML_INTERNAL_GENERIC_TYPE_H_
#define GMML_INTERNAL_GENERIC_TYPE_H_

#include <iosfwd>
#include <string>

#include "gmml/internal/stubs/common.h"

namespace gmml {
namespace internal {

// This class is only used by AmberTopFile and is probably not a great idea,
// in retrospect. It shouldn't be used by clients.
class GenericType {
  public:
    enum Type { kDouble, kInt, kString };

    GenericType() : type_(kDouble), dval_(kNotSet) {}
    GenericType(int n) : type_(kInt), dval_(n) {}
    GenericType(double n) : type_(kDouble), dval_(n) {}
    GenericType(const std::string& s) : type_(kString), sval_(s) {}

    GenericType& operator=(const GenericType& value);
    GenericType& operator=(double value);
    GenericType& operator=(const std::string& value);

    operator std::string() const;
    operator double() const;

    friend std::ostream& operator<<(std::ostream&, GenericType&);

  private:
    Type type_;
    double dval_;
    std::string sval_;
};

}  // namespace internal

using internal::GenericType;

}  // namespace gmml

#endif  // GMML_INTERNAL_GENERIC_TYPE_H_
