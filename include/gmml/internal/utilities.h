// Author: Robert Davis

#ifndef UTILITIES_H
#define UTILITIES_H

#include <exception>
#include <iomanip>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace gmml {

#undef DISALLOW_COPY_AND_ASSIGN
#define DISALLOW_COPY_AND_ASSIGN(TypeName) \
    TypeName(const TypeName&);             \
    void operator=(const TypeName&);

// This macro returns the number of elements in a static array. If you
// (wrongfully) pass it a pointer, the last line attempts to generate a
// compiler warning. I found this in the source code of Protocol Buffers.
#undef GOOGLE_ARRAYSIZE
#define GOOGLE_ARRAYSIZE(a) \
    ((sizeof(a) / sizeof(*(a))) / \
     static_cast<size_t>(!(sizeof(a) % sizeof(*(a)))))

// Standard Library Additions

// This makes static map creation easy with a call such as
// std::map<T,U> map = CreateMap<T,U>(t1, u1)(t2, u2), ...;
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

// Handy functors
struct DeletePtr {
    template<class T>
    void operator()(const T* ptr) const { delete ptr; }
};

struct DereferenceLess {
    template<typename PtrType>
    bool operator()(PtrType lhs, PtrType rhs) const { return *lhs < *rhs; }
};

// File utilities

class FileNotFoundException : public std::invalid_argument {
  public:
    explicit FileNotFoundException(const std::string& file_name)
            : std::invalid_argument("File not found: " + file_name),
              file_name_(file_name) {}

    ~FileNotFoundException() throw() {}

    std::string file_name() const { return file_name_; }

  private:
    std::string file_name_;
};

// Error utilities

// general-purpose status type
enum Status { kStatusOK, kStatusError };

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

extern bool kEnableWarnings;

inline void enable_warnings() {
    kEnableWarnings = true;
}

inline void disable_warnings() {
    kEnableWarnings = false;
}

void warning(const std::string& message);

void error(const std::string& message);

// String utilities

class ConversionException : public std::invalid_argument {
  public:
    explicit ConversionException(const std::string& what_arg)
            : std::invalid_argument(what_arg) {}
};

// remove spaces on either side
inline std::string& trim(std::string& str) {
    str.erase(str.find_last_not_of(" ") + 1);
    str.erase(0, str.find_first_not_of(" "));
    return str;
}

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

// Place val into str between str[index] and str[index + length]. The
// alignment of val in string is determined by the last argument, which should
// either be 'L' (left-aligned) or 'R' (right-aligned).
inline void set_in_string(std::string& str, const std::string& val,
                          size_t index, size_t length, char alignment) {
    if (val.size() > length || index + length - 1 >= str.size())
        return;
    if (alignment == 'R')
        str.replace(index + length - val.size(), val.size(), val);
    else
        str.replace(index, val.size(), val);
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
    return c >= '0' && c <= '9';
}

inline int char_to_number(char c) {
    if (is_number(c))
        return c - '0';
    return -1;
}

// An exception that is intended not to be caught. It should be used instead
// of something like std::abort to ensure that objects are destructed.
inline ExitException::ExitException(int status)
        : status_(status) {
    what_ = "Exit status: " + to_string(status);
}

// Numeric utilities

// A value for numbers to indicate that they are not set. This obviously isn't
// the safest thing to do, but it is useful in many situations.
const double kNotSet = 12345678.0;

template<typename T>
inline bool is_not_set(T val) {
    return val == kNotSet;
}

template<typename T>
inline bool is_set(T val) {
    return !is_not_set(val);
}

const double kPi = 3.14159265359;

inline double to_radians(double degrees) {
    return degrees/180.0*kPi;
}

inline double to_degrees(double radians) {
    return radians*180.0/kPi;
}

}  // namespace gmml

#endif  // UTILITIES_H
