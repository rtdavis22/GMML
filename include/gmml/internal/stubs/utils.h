// Author: Robert Davis
//
// This file should mostly only include things that are needed by public
// headers.

#ifndef GMML_INTERNAL_STUBS_UTILS_H_
#define GMML_INTERNAL_STUBS_UTILS_H_

#include <exception>
#include <sstream>
#include <string>

namespace gmml {

//
// String utilities
//
class ConversionException : public std::invalid_argument {
  public:
    explicit ConversionException(const std::string& what_arg)
            : std::invalid_argument(what_arg) {}
};

template<typename T>
inline std::string to_string(const T& val) {
    std::stringstream ss;
    if (ss << val)
        return ss.str();
    throw ConversionException("to_string: invalid conversion");
}

template<typename T>
inline T convert_string(const std::string& str) {
    T val;
    std::stringstream ss(str);
    if (ss >> val)
        return val;

    throw ConversionException("convert_string: invalid conversion of string " +
                              str);
}


template<typename T, typename U = T, typename V = T>
struct Triplet {
    Triplet(const T& first, const U& second, const V& third)
            : first(first), second(second), third(third) {}

    T first;
    U second;
    V third;
};

template<typename T, typename U = T, typename V = T>
struct TripletLess {
    bool operator()(const Triplet<T, U, V>& lhs,
                    const Triplet<T, U, V>& rhs) const {
        if (lhs.first == rhs.first) {
            if (lhs.second == rhs.second) {
                return lhs.third < rhs.third;
            } else {
                return lhs.second < rhs.second;
            }
        } else {
            return lhs.first < rhs.first;
        }
    }
};

template<typename T, typename U = T, typename V = T>
struct TripletPtrLess {
    bool operator()(const Triplet<T, U, V> *lhs,
                    const Triplet<T, U, V> *rhs) const {
        return TripletLess<T, U, V>()(*lhs, *rhs);
    }
};

}  // namespace gmml

#endif  // GMML_INTERNAL_STUBS_UTILS_H_
