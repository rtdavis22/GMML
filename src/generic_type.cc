#include "gmml/internal/generic_type.h"

#include <iostream>
#include <stdexcept>
#include <string>

namespace gmml
{

using std::string;

GenericType& GenericType::operator=(const GenericType& rhs) {
    if ((type_ != kString && rhs.type_ == kString) || 
            (type_ == kString && rhs.type_ != kString))
        throw std::logic_error("GenericType: error in assignment");
    if (type_ != kString)
        dval_ = rhs.dval_;
    else
        sval_ = rhs.sval_; 
    return *this;
}

GenericType& GenericType::operator=(double val) {
    if (type_ != kString)
        dval_ = val;
    else
        throw std::logic_error("GenericType: error in assignment");
    return *this;
}

GenericType& GenericType::operator=(const string& val) {
    if (type_ == kString)
        sval_ = val;
    else
        throw std::logic_error("GenericType: error in assignment");
    return *this;
}

GenericType::operator string() const {
    if (type_ == kString)
        return sval_;
    throw std::logic_error("GenericType: value is not a string");
}

GenericType::operator double() const {
    if (type_ != kString)
        return dval_;
    throw std::logic_error("GenericType: value is not a double");
}

std::ostream& operator<<(std::ostream &out, GenericType& t) {
    if (t.type_ == GenericType::kInt)
        out << static_cast<int>(t.dval_);
    else if (t.type_ == GenericType::kDouble)
        out << t.dval_;
    else
        out << t.sval_;
    return out;
}

} // namespace gmml 
