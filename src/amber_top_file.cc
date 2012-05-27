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
#include <stdexcept>
#include <string>
#include <vector>

#include "gmml/internal/environment.h"
#include "utilities.h"

namespace gmml {

using std::map;
using std::string;
using std::vector;

void AmberTopIntSection::append(const string& line) {
    int count = line.size()/width_;
    for (int j = 0; j < count; j++) {
        string s = line.substr(j*width_, width_);
        int num = convert_string<int>(s);
        elements_.push_back(num);
    }
}

void AmberTopIntSection::print(std::ostream& out) {
    out << "%FLAG " << name_ << std::endl;
    out << "%FORMAT(" << count_per_line_ << "I" << width_ << ")" << std::endl;
    if (elements_.empty()) {
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

int AmberTopIntSection::sum() const {
    return std::accumulate(elements_.begin(), elements_.end(), 0);
}

void AmberTopDoubleSection::append(const string& line) {
    int count = line.size()/width_;
    for (int j = 0; j < count; j++) {
        string s = line.substr(j*width_, width_);
        double num = convert_string<double>(s);
        elements_.push_back(num);
    }
}

void AmberTopDoubleSection::print(std::ostream& out) {
    out << "%FLAG " << name_ << std::endl;
    out << "%FORMAT(" << count_per_line_ << "E" << width_ << "." <<
           decimal_places_ << ")" << std::endl;
    if (elements_.empty()) {
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

void AmberTopStringSection::append(const string& line) {
    size_t count = line.size()/width_;
    for (size_t j = 0; j < count; j++) {
        string s = line.substr(j*width_, width_);
        trim(s);
        if (s.empty()) {
            return;
        }
        elements_.push_back(s);
    }
}

void AmberTopStringSection::print(std::ostream& out) {
    out << "%FLAG " << name_ << std::endl;
    out << "%FORMAT(" << count_per_line_ << "a" << width_ << ")" << std::endl;
    if (elements_.empty()) {
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

// Inline these prolly
AmberTopIntSection *AmberTopFile::create_int_section(
        const string& name, int count_per_line, int width) {
    AmberTopIntSection *section = new AmberTopIntSection(name, count_per_line,
                                                         width);
    add_int_section(section);
    return section;
}

AmberTopIntSection *AmberTopFile::create_int_section(const string& name,
                                                     const string& format) {
    size_t i = format.find("I");
    int count_per_line = convert_string<size_t>(format.substr(0, i));
    int width = convert_string<size_t>(format.substr(i + 1));

    return create_int_section(name, count_per_line, width);
}

AmberTopDoubleSection *AmberTopFile::create_double_section(
        const string& name, int count_per_line, int width, int decimal_places) {
    AmberTopDoubleSection *section =
            new AmberTopDoubleSection(name, count_per_line, width,
                                      decimal_places);
    add_double_section(section);
    return section;
}

AmberTopDoubleSection *AmberTopFile::create_double_section(
        const string& name, const string& format) {
    size_t i = format.find("E");
    size_t j = format.find(".");
    int count_per_line = convert_string<size_t>(format.substr(0, i));
    int width = convert_string<size_t>(format.substr(i + 1, j - i - 1));
    int decimal_places = convert_string<size_t>(format.substr(j + 1));

    return create_double_section(name, count_per_line, width, decimal_places);
}

AmberTopStringSection *AmberTopFile::create_string_section(
        const string& name, int count_per_line, int width) {
    AmberTopStringSection *section =
            new AmberTopStringSection(name, count_per_line, width);
    add_string_section(section);
    return section;
}

AmberTopStringSection *AmberTopFile::create_string_section(
        const string& name, const string& format) {
    size_t i = format.find("a");
    int count_per_line = convert_string<size_t>(format.substr(0, i));
    int width = convert_string<size_t>(format.substr(i + 1));

    return create_string_section(name, count_per_line, width);
}

struct IsName {
    explicit IsName(const string& name) : name(name) {}

    bool operator()(const AmberTopSection *section) {
        return section->name() == name;
    }

    string name;
};

template<typename T>
void remove_from_arr(T& vector, const string& name) {
    typename T::iterator it = std::remove_if(vector.begin(), vector.end(), IsName(name));
    vector.erase(it, vector.end());
}

bool AmberTopFile::remove_section(const string& name) {
    remove_from_arr(int_sections_, name);
    remove_from_arr(double_sections_, name);
    remove_from_arr(string_sections_, name);
    vector<AmberTopSection*>::iterator it;
    it = section_order_.erase(std::remove_if(section_order_.begin(),
                                             section_order_.end(),
                                             IsName(name)),
                              section_order_.end());
    // delete these
    return true;
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
    // Ignore the version.
    extract_version(line);
    while (true) {
        while (get_card_type(line) != kFlagCard && !input.eof())
            getline(input, line);
        if (input.eof()) {
            throw std::invalid_argument("Invalid topology file.");
            return;
        }
        string title = extract_title(line);
        trim(title);
        if (getline(input, line) && get_card_type(line) != kFormatCard) {
            throw std::invalid_argument("Bad format in section " + title);
            continue;
        }
        string format = extract_format(line);

        AmberTopSection *section;
        switch (get_section_type(format)) {
            case kIntSection:
                section = create_int_section(title, format);
                break;
            case kDoubleSection:
                section = create_double_section(title, format);
                break;
            case kStringSection:
                section = create_string_section(title, format);
                break;
        }
        process_section(input, section);
        if (input.eof())
            return;
    }
}

void AmberTopFile::process_section(std::istream& input,
                                   AmberTopSection *section) {
    string line;
    while (input.peek() != '%' && getline(input, line)) {
        section->append(line);
    }
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
    for (size_t i = 0; i < section_order_.size(); i++)
        section_order_[i]->print(out);
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
