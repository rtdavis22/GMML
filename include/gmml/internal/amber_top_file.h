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
#include <iosfwd>
#include <map>
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

// Should implement Writeable, but there is some problem with making write
// const and GenericType.
class AmberTopFile : public Readable  {
  public:
    // The three types of sections that are possible in the file.
    enum SectionType { kStringSection, kIntSection, kDoubleSection };

    // These classify the control lines of the file ("%FLAG", "%FORMAT")
    enum CardType { kFlagCard, kFormatCard, kOtherCard };

    // Create a topology file with no sections
    AmberTopFile() {}

    // Create a topology file by reading in a file on disk.
    explicit AmberTopFile(const std::string& file_name) {
        Readable::read(file_name);
    }

    virtual ~AmberTopFile() {
        STLDeleteContainerPointers(section_order_.begin(),
                                   section_order_.end());
    }

    // Create and return a pointer to a section with the given name and
    // FORTRAN format string.
    AmberTopIntSection *create_int_section(const std::string& name,
                                           const std::string& format);
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

    // Remove a section by specifying the section name.
    bool remove_section(const std::string& name);

    AmberTopIntSection *get_int_section(const std::string& name) {
        for (int i = 0; i < int_sections_.size(); i++) {
            if (int_sections_[i]->name() == name)
                return int_sections_[i];
        }
        return NULL;
    }

    AmberTopDoubleSection *get_double_section(const std::string& name) {
        for (int i = 0; i < double_sections_.size(); i++) {
            if (double_sections_[i]->name() == name)
                return double_sections_[i];
        }
        return NULL;
    }

    AmberTopStringSection *get_string_section(const std::string& name) {
        for (int i = 0; i < string_sections_.size(); i++) {
            if (string_sections_[i]->name() == name)
                return string_sections_[i];
        }
        return NULL;
    }

    // Returns true if a section with the given name exists in the file.
    //bool exists(const std::string& name) const {
    //    return sections_.find(name) != sections_.end();
    //}

    // Sort the sections of the file according to a sorting criterion.
    template <class Compare>
    void sort(Compare comp) {
        std::sort(section_order_.begin(), section_order_.end(), comp);
    }

    void print(const std::string& file_name);
    void print();

  private:
    void read(const std::string& file_name);
    virtual void read(std::istream&);
    virtual void write(std::ostream&);
    enum SectionType get_section_type(const std::string&);
    enum CardType get_card_type(const std::string&);
    std::string extract_format(const std::string&);

    void add_int_section(AmberTopIntSection *section) {
        int_sections_.push_back(section);
        section_order_.push_back(section);
    }

    void add_double_section(AmberTopDoubleSection *section) {
        double_sections_.push_back(section);
        section_order_.push_back(section);
    }

    void add_string_section(AmberTopStringSection *section) {
        string_sections_.push_back(section);
        section_order_.push_back(section);
    }

    std::string extract_title(const std::string&);
    std::string extract_version(const std::string&);
    void process_section(std::istream&, AmberTopSection *section);
    std::string get_version_string() const;

    std::vector<AmberTopIntSection*> int_sections_;
    std::vector<AmberTopDoubleSection*> double_sections_;
    std::vector<AmberTopStringSection*> string_sections_;
    std::vector<AmberTopSection*> section_order_;

    DISALLOW_COPY_AND_ASSIGN(AmberTopFile);
};

}  // namespace gmml

#endif  // GMML_INTERNAL_AMBER_TOP_FILE_H_
