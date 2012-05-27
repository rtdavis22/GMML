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

#ifndef GMML_INTERNAL_STUBS_COMMON_H_
#define GMML_INTERNAL_STUBS_COMMON_H_

namespace gmml {

#undef DISALLOW_COPY_AND_ASSIGN
#define DISALLOW_COPY_AND_ASSIGN(TypeName) \
    TypeName(const TypeName&);             \
    void operator=(const TypeName&);

// This macro returns the number of elements in a static array. If you
// (wrongfully) pass it a pointer, the last line attempts to generate a
// compiler warning. This was found in the source code of Protocol Buffers.
#undef ARRAY_SIZE
#define ARRAY_SIZE(a) \
    ((sizeof(a) / sizeof(*(a))) / \
     static_cast<size_t>(!(sizeof(a) % sizeof(*(a)))))

// A general-purpose status type.
enum Status { kStatusOK, kStatusError };

// Numeric utilities

// A value for numbers to indicate that they are not set. This obviously isn't
// the safest thing to do, but it is useful in many situations.
const double kNotSet = 123456789.0;

template<typename T>
inline bool is_not_set(T val) {
    return val == kNotSet;
}

template<typename T>
inline bool is_set(T val) {
    return !is_not_set(val);
}

const double kPi = 3.1415926535897932385;

template<typename T>
inline T to_radians(T degrees) {
    return degrees/180.0*kPi;
}

template<typename T>
inline T to_degrees(T radians) {
    return radians*180.0/kPi;
}

// Similar functions are available in <cctype>, but I don't particularly want
// to pull the whole header in.
inline bool is_uppercase(char c) {
    return 'A' <= c && c <= 'Z';
}

inline bool is_lowercase(char c) {
    return 'a' <= c && c <= 'z';
}

inline bool is_letter(char c) {
    return is_uppercase(c) || is_lowercase(c);
}

// Handy functors
struct DeletePtr {
    template<class T>
    void operator()(const T* ptr) const { delete ptr; }
};

struct DereferenceLess {
    template<typename PtrType>
    bool operator()(PtrType lhs, PtrType rhs) const { return *lhs < *rhs; }
};

// This is here so we don't have to include all of <algorithm> just to use
// std::for_each() (with DeletePtr).
template <class ForwardIterator>
void STLDeleteContainerPointers(ForwardIterator begin, ForwardIterator end) {
    while (begin != end) {
        ForwardIterator temp = begin;
        ++begin;
        delete *temp;
    }
}

enum LogLevel {
  LOGLEVEL_INFO,     // Informational.
  LOGLEVEL_WARNING,  // Warns about issues that, although not technically a
                     // problem now, could cause problems in the future.
  LOGLEVEL_ERROR,    // An error occurred which should never happen during
                     // normal use.
  LOGLEVEL_FATAL,    // An error occurred from which the library cannot
                     // recover.

#ifdef NDEBUG
  LOGLEVEL_DFATAL = LOGLEVEL_ERROR
#else
  LOGLEVEL_DFATAL = LOGLEVEL_FATAL
#endif
};

extern LogLevel kDefaultLogLevel;

inline void set_log_level(LogLevel level) {
    kDefaultLogLevel = level;
}

}  // namespace gmml

#endif  // GMML_INTERNAL_STUBS_COMMON_H_
