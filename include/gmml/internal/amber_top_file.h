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
// This file contains a representation of the AMBER topology file.
// See http://ambermd.org/formats.html#topology for the specification.
// Note that this file only specifies the format of the file. The contents of
// the file do not necessarily contain anything to do with AMBER or the
// topology of anything. See amber_top_builder.h for that stuff.

#ifndef GMML_INTERNAL_AMBER_TOP_FILE_H_
#define GMML_INTERNAL_AMBER_TOP_FILE_H_

#include <algorithm>
#include <functional>
#include <iosfwd>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "gmml/internal/stubs/common.h"
#include "gmml/internal/stubs/file.h"

namespace gmml {

/**
 This class represents a type of section in the AMBER topology file.

 Subclasses of this class determine how the data in the section is stored,
 read, and written.
 */
class AmberTopSection {
  public:
    virtual ~AmberTopSection() {}

    // Make this append_line_of_data?
    /**
     Appends a line of data from topology file to the section. This will
     typically consist of multiple elements.
     */
    virtual void append(const std::string&) = 0;

    /**
     Removes all the elements of this section.
     */
    virtual void clear() = 0;

    /**
     Prints this section to the given output stream.
     */
    virtual void print(std::ostream& out) = 0;

    /**
     Returns the number of elements in the section.
     */
    virtual size_t size() const = 0;

    /**
     Returns name of the section. It follows "%FLAG" in the topology file.
     */
    std::string name() const { return name_; }

    /**
     Returns the number of elements per line.
     */
    int count_per_line() const { return count_per_line_; }

    /**
     Returns the width of each elements (in columns).
     */
    int width() const { return width_; }

  protected:
    /**
     Creates a section with the given name, number of elements per line,
     and width of each elements (in columns).
     */
    AmberTopSection(const std::string& name, int count_per_line, int width)
            : name_(name), count_per_line_(count_per_line), width_(width) {
    }

  private:
    std::string name_;
    int count_per_line_;
    int width_;

    DISALLOW_COPY_AND_ASSIGN(AmberTopSection);
};

/**
 This is an integer section in the topology file.
 */
class AmberTopIntSection : public AmberTopSection {
  public:
    AmberTopIntSection(const std::string& name, int count_per_line, int width)
            : AmberTopSection(name, count_per_line, width) {
    }

    int sum() const;

    void set(int index, int value) { elements_.at(index) = value; }

    int get(int index) const { return elements_.at(index); }

    virtual void append(const std::string&);

    void insert(int value) { elements_.push_back(value); }

    virtual void print(std::ostream&);

    virtual void clear() { elements_.clear(); }

    virtual size_t size() const { return elements_.size(); }

  private:
    std::vector<int> elements_;

    DISALLOW_COPY_AND_ASSIGN(AmberTopIntSection);
};

/**
 This is a floating point section in the topology file.
 */
class AmberTopDoubleSection : public AmberTopSection {
  public:
    AmberTopDoubleSection(const std::string& name, int count_per_line,
                          int width, int decimal_places)
            : AmberTopSection(name, count_per_line, width),
              decimal_places_(decimal_places) {
    }

    void set(int index, double value) { elements_.at(index) = value; }

    double get(int index) const { return elements_.at(index); }

    virtual void append(const std::string&);

    void insert(double value) { elements_.push_back(value); }

    virtual void print(std::ostream&);

    virtual void clear() { elements_.clear(); }

    virtual size_t size() const { return elements_.size(); }

  private:
    int decimal_places_;
    std::vector<double> elements_;

    DISALLOW_COPY_AND_ASSIGN(AmberTopDoubleSection);
};

/**
 This is a string section in the topology file.
 */
class AmberTopStringSection : public AmberTopSection {
  public:
    AmberTopStringSection(const std::string& name, int count_per_line,
                          int width)
            : AmberTopSection(name, count_per_line, width) {
    }

    void set(int index, const std::string& value) {
        elements_.at(index) = value;
    }

    std::string get(int index) const { return elements_.at(index); }

    virtual void append(const std::string&);

    void insert(const std::string& value) { elements_.push_back(value); }

    virtual void print(std::ostream&);

    virtual void clear() { return elements_.clear(); }

    virtual size_t size() const { return elements_.size(); }

  private:
    std::vector<std::string> elements_;

    DISALLOW_COPY_AND_ASSIGN(AmberTopStringSection);
};

// Put this at top of file~!!!
/**
 This class represents an AMBER topology file.
 */
class AmberTopFile : public Readable, public Writeable {
  public:
    /**
     Creates a topology file with no sections.
     */
    AmberTopFile();

    /**
     Creates a topology file from the given File.
     */
    explicit AmberTopFile(const File& file);

    virtual ~AmberTopFile();

    /**
     Creates an integer section with the given name and FORTRAN format string,
     appends it to the topology file, and returns it.
     */
    AmberTopIntSection *create_int_section(const std::string& name,
                                           const std::string& format);

    /**
     Creates an integer section with the given name, number of elements per
     line, and element width (in columns).
     */
    AmberTopIntSection *create_int_section(const std::string& name,
                                           int count_per_line, int width);

    AmberTopDoubleSection *create_double_section(const std::string& name,
                                                 const std::string& format);
    AmberTopDoubleSection *create_double_section(const std::string& name,
                                                 int count_per_line, int width,
                                                 int decimal_places);

    AmberTopStringSection *create_string_section(const std::string& name,
                                                 const std::string& format);
    AmberTopStringSection *create_string_section(const std::string& name,
                                                 int count_per_line, int width);

    /**
     Removes the section of the topology file with the given name.
     */
    bool remove_section(const std::string& name);

    /**
     Returns the integer section with the given name.
     */
    AmberTopIntSection *get_int_section(const std::string& name);
    AmberTopDoubleSection *get_double_section(const std::string& name);
    AmberTopStringSection *get_string_section(const std::string& name);

    /**
     Sorts the sections of the file, given a comparison function.
     */
    void sort(bool (*comp)(const AmberTopSection *lhs,
                           const AmberTopSection *rhs));

  private:
    struct Impl;
    std::auto_ptr<Impl> impl_;

    //void read(const std::string& file_name);

    virtual void read(std::istream&);
    virtual void write(std::ostream&) const;
/*
    enum SectionType get_section_type(const std::string&);
    enum CardType get_card_type(const std::string&);
    std::string extract_format(const std::string&);

    std::string extract_title(const std::string&);
    std::string extract_version(const std::string&);
    void process_section(std::istream&, AmberTopSection *section);
*/
    DISALLOW_COPY_AND_ASSIGN(AmberTopFile);
};

/*
inline AmberTopIntSection *AmberTopFile::get_int_section(
        const std::string& name) {
    for (int i = 0; i < int_sections_.size(); i++) {
        if (int_sections_[i]->name() == name)
            return int_sections_[i];
    }
    return NULL;
}

inline AmberTopDoubleSection *AmberTopFile::get_double_section(
        const std::string& name) {
    for (int i = 0; i < double_sections_.size(); i++) {
        if (double_sections_[i]->name() == name)
            return double_sections_[i];
    }
    return NULL;
}

inline AmberTopStringSection *AmberTopFile::get_string_section(
        const std::string& name) {
    for (int i = 0; i < string_sections_.size(); i++) {
        if (string_sections_[i]->name() == name)
            return string_sections_[i];
    }
    return NULL;
}
*/

}  // namespace gmml

#endif  // GMML_INTERNAL_AMBER_TOP_FILE_H_
