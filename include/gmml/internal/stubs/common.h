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

inline double to_radians(double degrees) {
    return degrees/180.0*kPi;
}

inline double to_degrees(double radians) {
    return radians*180.0/kPi;
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
  LOGLEVEL_INFO,     // Informational.  This is never actually used by
                     // libprotobuf.
  LOGLEVEL_WARNING,  // Warns about issues that, although not technically a
                     // problem now, could cause problems in the future.  For
                     // example, a // warning will be printed when parsing a
                     // message that is near the message size limit.
  LOGLEVEL_ERROR,    // An error occurred which should never happen during
                     // normal use.
  LOGLEVEL_FATAL,    // An error occurred from which the library cannot
                     // recover.  This usually indicates a programming error
                     // in the code which calls the library, especially when
                     // compiled in debug mode.

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
