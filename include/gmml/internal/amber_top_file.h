// Author: Robert Davis

#ifndef GMML_INTERNAL_AMBER_TOP_FILE_H_
#define GMML_INTERNAL_AMBER_TOP_FILE_H_

#include <algorithm>
#include <iosfwd>
#include <map>
#include <string>
#include <vector>

#include <boost/shared_ptr.hpp>

#include "generic_type.h"
#include "utilities.h"

namespace gmml {

// This class represents a section in the topology file. It holds elements
// of GenericType type, which can be strings, ints, or doubles. It also
// contains the format of the section, which determines how the section is
// read and written.
class AmberTopSection {
  public:
    virtual ~AmberTopSection() {}

    // This is the name of the section. It follows "%FLAG" in the
    // topology file.
    std::string name() const { return name_; }

    // Append a line from topology file, which will typically consist of
    // multiple elements, to the topology file.
    virtual Status append(const std::string&) = 0;

    virtual void insert(GenericType element) { elements_.push_back(element); }
    void clear() { elements_.clear(); }
    virtual void print(std::ostream& out) = 0;

    // Compute sum and minimum and maximum values of the section.
    GenericType sum() const;
    GenericType max() const;
    GenericType min() const;

    size_t size() const { return elements_.size(); }
    GenericType& operator[](size_t i) { return elements_[i]; }
    const GenericType& operator[](size_t i) const { return elements_[i]; }

  protected:
    explicit AmberTopSection(const std::string& name) : name_(name) {}
    AmberTopSection(const std::string& name, size_t size)
            : name_(name), elements_(size) {}

    std::string name_;
    size_t count_per_line_;
    size_t width_;
    std::vector<GenericType> elements_;
};

class AmberTopIntSection : public AmberTopSection {
  public:
    AmberTopIntSection(const std::string& name, const std::string& format,
                       size_t size = 0);
    Status append(const std::string&);
    virtual void print(std::ostream&);
};

class AmberTopDoubleSection : public AmberTopSection {
  public:
    AmberTopDoubleSection(const std::string& name, const std::string& format,
                          size_t size = 0);
    Status append(const std::string&);
    void print(std::ostream&);
  private:
    size_t decimal_places_;
};

class AmberTopStringSection : public AmberTopSection {
  public:
    AmberTopStringSection(const std::string& name, const std::string& format,
                          size_t size = 0);
    Status append(const std::string&);
    void print(std::ostream&);
};

class AmberTopFile {
  public:
    // The 3 types of sections that are possible in the file.
    enum SectionType { kStringSection, kIntSection, kDoubleSection };

    // These classify the control lines of the file ("%FLAG", "%FORMAT")
    enum CardType { kFlagCard, kFormatCard, kOtherCard };

    typedef boost::shared_ptr<AmberTopSection> SectionPtr;
    typedef std::map<std::string, SectionPtr> SectionMap;

    // Create a topology file with no sections
    AmberTopFile() {}

    // Create a topology file by reading in a file on disk.
    explicit AmberTopFile(const std::string& file_name) { read(file_name); }

    virtual ~AmberTopFile() {}

    // Write the topology file to disk or stdout.
    void print(const std::string& file_name);
    void print();

    // Create and return a pointer to a section with the given name and
    // FORTRAN format string.
    SectionPtr create_section(const std::string& name,
                              const std::string& format, size_t size = 0);

    // Remove a section by specifying the section name.
    bool remove_section(const std::string& name);

    // Returns true if a section with the given name exists in the file.
    bool exists(const std::string& name) const {
        return sections_.find(name) != sections_.end();
    }

    // Sort the sections of the file according to a sorting criterion.
    template <class Compare>
    void sort(Compare comp) {
        std::sort(section_list_.begin(), section_list_.end(), comp);
    }

    // Access the sections of the file
    AmberTopSection& operator[](const std::string& s) { return *sections_[s]; }
    const AmberTopSection& operator[](const std::string& s) const {
        return *sections_.find(s)->second;
    }

  private:
    void read(const std::string& file_name);
    void read(std::istream&);
    void write(std::ostream&);
    enum SectionType get_section_type(const std::string&);
    enum CardType get_card_type(const std::string&);
    std::string extract_format(const std::string&);
    std::string extract_title(const std::string&);
    std::string extract_version(const std::string&);
    void process_section(std::istream&, SectionPtr);
    std::string get_version_string() const;

    SectionMap sections_;
    std::vector<std::string> section_list_;
};

}  // namespace gmml

#endif  // GMML_INTERNAL_AMBER_TOP_FILE_H_
