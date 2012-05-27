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

#ifndef GMML_SRC_UTILITIES_H_
#define GMML_SRC_UTILITIES_H_

#include <exception>
#include <iomanip>
#include <iostream>
#include <istream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "gmml/internal/stubs/common.h"
#include "gmml/internal/stubs/utils.h"

namespace gmml {

// Standard Library Additions

// This makes static map creation easy with a call such as
// std::map<T, U> map = CreateMap<T,U>(t1, u1)(t2, u2), ...;
template<typename T, typename U>
class CreateMap {
  public:
    CreateMap(const T& key, const U& value) { map_[key] = value; }

    CreateMap<T, U>& operator()(const T& key, const U& value) {
        map_[key] = value;
        return *this;
    }

    operator std::map<T, U>() { return map_; }

  private:
    std::map<T, U> map_;
};

// A function to efficiently add or update a map. If a key exists in the map
// already, updating its value with operator[] is faster. Otherwise, insert()
// is faster. This function does the right thing and returns an iterator
// to the added or updated pair. This code is from Effective STL by
// Scott Meyers.
template<typename MapType, typename KeyArgType, typename ValueArgType>
typename MapType::iterator add_or_update_map(MapType& map,
                                             const KeyArgType& key,
                                             const ValueArgType& value) {
    typename MapType::iterator it = map.lower_bound(key);
    if (it != map.end() && !(map.key_comp()(key, it->first))) {
        it->second = value;
        return it;
    } else {
        return map.insert(it, typename MapType::value_type(key, value));
    }
}


// Error utilities

// An exception that is not caught internally
struct ExitException : public std::exception {
  public:
    explicit ExitException(int status);

    int status() const { return status_; }

    virtual const char *what() const throw() { return what_.c_str(); }

    virtual ~ExitException() throw() {}

  private:
    int status_;
    std::string what_;
};

inline void die() {
    throw ExitException(-1);
}

// String utilities

// Remove spaces on both sides of the string.
inline std::string& trim(std::string& str) {
    str.erase(str.find_last_not_of(" ") + 1);
    str.erase(0, str.find_first_not_of(" "));
    return str;
}

// Place val into str between str[index] and str[index + length]. The
// alignment of val in string is determined by the last argument, which should
// either be 'L' (left-aligned) or 'R' (right-aligned).
inline void set_in_string(std::string& str, const std::string& val,
                          size_t index, size_t length, char alignment) {
    if (val.size() > length || index + length - 1 >= str.size()) {
        return;
    }
    if (alignment == 'R') {
        str.replace(index + length - val.size(), val.size(), val);
    } else {
        str.replace(index, val.size(), val);
    }
}

// This overload is for placing an integer in the string.
inline void set_in_string(std::string& str, int number, size_t index,
                          int length, char alignment) {
    set_in_string(str, to_string(number), index, length, alignment);
}

// This overload is for placing a double in a string.
inline void set_in_string(std::string& str, double number, size_t index,
                          size_t length, char alignment, int decimal_places) {
    std::stringstream ss;
    ss << std::setprecision(decimal_places) << std::fixed << number;
    set_in_string(str, ss.str(), index, length, alignment);
}

// Split a string into substring using the delimiter
inline std::vector<std::string>& split(const std::string& str, char delimiter,
                                       std::vector<std::string>& elements) {
    std::stringstream ss(str);
    std::string item;
    while (getline(ss, item, delimiter))
        elements.push_back(item);
    return elements;
}

inline bool is_number(char c) {
    return '0' <= c && c <= '9';
}

inline int char_to_number(char c) {
    if (is_number(c))
        return c - '0';
    return -1;
}

// An exception that is intended not to be caught. It should be used instead
// of something like std::abort so that the call stack is unwound.
inline ExitException::ExitException(int status) : status_(status) {
    what_ = "Exit status: " + to_string(status);
}


template<class charT, class traits>
inline std::basic_istream<charT, traits>&
ignore_line(std::basic_istream<charT, traits>& strm) {
    strm.ignore(std::numeric_limits<int>::max(), strm.widen('\n'));
    return strm;
}

class make_string {
  public:
    template <typename T>
    make_string& operator<<(T const& datum) {
        buffer_ << datum;
        return *this;
    }
    operator std::string () const {
        return buffer_.str();
    }

  private:
    std::ostringstream buffer_;
};

}  // namespace gmml

#endif  // GMML_SRC_UTILITIES_H_
