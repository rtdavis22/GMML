#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "gmml/internal/amber_top_file.h"

#include <ctime>

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <string>
#include <vector>

#include "gmml/internal/environment.h"
#include "gmml/internal/generic_type.h"
#include "utilities.h"

namespace gmml {

using std::map;
using std::string;
using std::vector;

GenericType AmberTopSection::sum() const {
    return std::accumulate(elements_.begin(), elements_.end(), 0.0);
}

GenericType AmberTopSection::max() const {
    return *std::max_element(elements_.begin(), elements_.end());
}

GenericType AmberTopSection::min() const {
    return *std::min_element(elements_.begin(), elements_.end());
}

AmberTopIntSection::AmberTopIntSection(const string& name, const string& format,
                                       size_t size)
        : AmberTopSection(name, size) {
    size_t i = format.find("I");
    count_per_line_ = convert_string<size_t>(format.substr(0, i));
    width_ = convert_string<size_t>(format.substr(i + 1));
}

AmberTopDoubleSection::AmberTopDoubleSection(const string& name,
                                             const string& format, size_t size)
        : AmberTopSection(name, size) {
    size_t i = format.find("E");
    size_t j = format.find(".");
    count_per_line_ = convert_string<size_t>(format.substr(0, i));
    width_ = convert_string<size_t>(format.substr(i + 1, j - i - 1));
    decimal_places_ = convert_string<size_t>(format.substr(j + 1));
}

AmberTopStringSection::AmberTopStringSection(const string& name,
                                             const string& format, size_t size)
        : AmberTopSection(name, size) {
    size_t i = format.find("a");
    count_per_line_ = convert_string<size_t>(format.substr(0, i));
    width_ = convert_string<size_t>(format.substr(i + 1));
}


Status AmberTopIntSection::append(const string& line) {
    int count = line.size()/width_;
    for (int j = 0; j < count; j++) {
        string s = line.substr(j*width_, width_);
        int num = convert_string<int>(s);
        elements_.push_back(GenericType(num));
    }
    return kStatusOK;
}

void AmberTopIntSection::print(std::ostream& out) {
    out << "%FLAG " << name_ << std::endl;
    out << "%FORMAT(" << count_per_line_ << "I" << width_ << ")" << std::endl;
    if (size() == 0) {
        out << std::endl;
        return;
    }
    size_t lines = size()/count_per_line_;
    for (size_t i = 0; i < lines; i++) {
        for (size_t j = 0; j < count_per_line_; j++) {
            int val = elements_[i*count_per_line_ + j];
            out << std::right << std::setw(width_) << val;
        }
        out << std::endl;
    }
    if (size()%count_per_line_ == 0)
        return;
    for (size_t i = 0; i < size()%count_per_line_; i++) {
        int val = elements_[lines*count_per_line_ + i];
        out << std::right << std::setw(width_) << val;
    }
    out << std::endl;
}

Status AmberTopDoubleSection::append(const string& line) {
    int count = line.size()/width_;
    for (int j = 0; j < count; j++) {
        string s = line.substr(j*width_, width_);
        double num = convert_string<double>(s);
        elements_.push_back(GenericType(num));
    }
    return kStatusOK;
}

void AmberTopDoubleSection::print(std::ostream& out) {
    out << "%FLAG " << name_ << std::endl;
    out << "%FORMAT(" << count_per_line_ << "E" << width_ << "." <<
           decimal_places_ << ")" << std::endl;
    if (size() == 0) {
        out << std::endl;
        return;
    }
    size_t lines = size()/count_per_line_;
    for (size_t i = 0; i < lines; i++) {
        for (size_t j = 0; j < count_per_line_; j++)
            out << std::right << std::setw(width_) <<
                   std::setprecision(decimal_places_) << std::scientific <<
                   elements_[i*count_per_line_ + j];
        out << std::endl;
    }
    if (size()%count_per_line_ == 0)
        return;
    for (size_t i = 0; i < size()%count_per_line_; i++)
        out << std::right << std::setw(width_) <<
               std::setprecision(decimal_places_) <<
               std::scientific <<
               elements_[(lines)*count_per_line_ + i];
    out << std::endl;
}

Status AmberTopStringSection::append(const string& line) {
    size_t count = line.size()/width_;
    for (size_t j = 0; j < count; j++) {
        string s = line.substr(j*width_, width_);
        trim(s);
        if (s.empty())
            return kStatusOK;
        elements_.push_back(GenericType(s));
    }
    return kStatusOK;
}

void AmberTopStringSection::print(std::ostream& out) {
    out << "%FLAG " << name_ << std::endl;
    out << "%FORMAT(" << count_per_line_ << "a" << width_ << ")" << std::endl;
    if (size() == 0) {
        out << std::endl;
        return;
    }
    size_t lines = size()/count_per_line_;
    for (size_t i = 0; i < lines; i++) {
        for (size_t j = 0; j < count_per_line_; j++)
            out << std::setw(width_) << std::left <<
                   elements_[i*count_per_line_ + j];
        out << std::endl;
    }
    if (size()%count_per_line_ == 0)
        return;
    for (size_t i = 0; i < size()%count_per_line_; i++)
        out << std::setw(width_) << std::left <<
               elements_[(lines)*count_per_line_ + i];
    out << std::endl;
}

void AmberTopFile::read(const string& file_name) {
    std::ifstream stream(find_file(file_name).c_str());
    read(stream);
    stream.close();
}

AmberTopFile::SectionPtr AmberTopFile::create_section(const string& name,
                                                      const string& format,
                                                      size_t size) {
    SectionPtr s;
    switch (get_section_type(format)) {
      case kDoubleSection:
        s.reset(new AmberTopDoubleSection(name, format, size));
        break;
      case kIntSection:
        s.reset(new AmberTopIntSection(name, format, size));
        break;
      case kStringSection:
        s.reset(new AmberTopStringSection(name, format, size));
        break;
    }
    sections_.insert(SectionMap::value_type(s->name(), s));
    section_list_.push_back(s->name());
    return s;
}

bool AmberTopFile::remove_section(const string& name) {
    section_list_.erase(std::remove(section_list_.begin(), section_list_.end(),
                                    name),
                        section_list_.end());
    return sections_.erase(name);
}

AmberTopFile::SectionType AmberTopFile::get_section_type(const string& line) {
    if (line.find('I') != string::npos)
        return kIntSection;
    else if (line.find('E') != string::npos)
        return kDoubleSection;
    else
        return kStringSection;
}

AmberTopFile::CardType AmberTopFile::get_card_type(const string& line) {
    if (line.find("%FLAG") == 0)
        return kFlagCard;
    else if (line.find("%FORMAT") == 0)
        return kFormatCard;
    else
        return kOtherCard;
}

string AmberTopFile::extract_title(const string& line) {
    if (line.size() > 6) {
        string ret = string(line.substr(6));
        trim(ret);
        return ret;
    }
    return string("");
}

string AmberTopFile::extract_format(const string& line) {
    int left = line.find("(");
    int right = line.find(")");
    string ret = line.substr(left + 1, right - left - 1);
    trim(ret);
    return ret;
}

string AmberTopFile::extract_version(const string& line) {
    if (line.size() > 9)
        return line.substr(9);
    return "";
}

void AmberTopFile::read(std::istream& input) {
    string line;
    getline(input, line);
    // Ignore the version
    extract_version(line);
    while (true) {
        while (get_card_type(line) != kFlagCard && !input.eof())
            getline(input, line);
        if (input.eof()) {
            error("Invalid topology file");
            return;
        }
        string title = extract_title(line);
        trim(title);
        if (getline(input, line) && get_card_type(line) != kFormatCard) {
            error("Bad format in section " + title);
            continue;
        }
        string format = extract_format(line);

        SectionPtr section;
        switch (get_section_type(format)) {
            case kIntSection:
                section.reset(new AmberTopIntSection(title, format));
                break;
            case kDoubleSection:
                section.reset(new AmberTopDoubleSection(title, format));
                break;
            case kStringSection:
                section.reset(new AmberTopStringSection(title, format));
                break;
        }
        process_section(input, section);
        if (input.eof())
            return;
    }
}

void AmberTopFile::process_section(std::istream& input, SectionPtr section) {
    string line;
    while (input.peek() != '%' && getline(input, line)) {
        if (section->append(line) != 0) {
            error("Error processing section " + section->name());
            return;
        }
    }
    section_list_.push_back(section->name());
    sections_.insert(SectionMap::value_type(section->name(), section));
}

void AmberTopFile::print() {
    write(std::cout);
}

void AmberTopFile::print(const string& file_name) {
    std::ofstream out;
    out.open(file_name.c_str());
    write(out);
    out.close();
}

void AmberTopFile::write(std::ostream& out) {
    out << "%VERSION " << get_version_string() << std::endl;
    for (size_t i = 0; i < section_list_.size(); i++)
       sections_[section_list_[i]]->print(out);
}

std::string AmberTopFile::get_version_string() const {
    char buffer[80];
    tm time_info;
    time_t ltime;
    time(&ltime);
#ifdef HAVE_LOCALTIME_R
    localtime_r(&ltime, &time_info);
#else
    time_info = *localtime(&ltime);
#endif
    strftime(buffer, 80,
             " VERSION_STAMP = V0001.000  DATE = %m/%d/%y  %H:%M:%S\0",
             &time_info);
    return string(buffer);
}

}  // namespace gmml
