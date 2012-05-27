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
//
// This file should mostly only include things that are needed by public
// headers.

#ifndef GMML_INTERNAL_STUBS_UTILS_H_
#define GMML_INTERNAL_STUBS_UTILS_H_

#include <exception>
#include <sstream>
#include <stdexcept>
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
