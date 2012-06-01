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

#include "gmml/internal/environment.h"
#include "utilities.h"

namespace gmml {

using std::map;
using std::string;
using std::vector;

void AmberTopIntSection::append(const string& line) {
    int count = line.size()/width();
    for (int j = 0; j < count; j++) {
        string s = line.substr(j*width(), width());
        int num = convert_string<int>(s);
        elements_.push_back(num);
    }
}

void AmberTopIntSection::print(std::ostream& out) {
    out << "%FLAG " << name() << std::endl;
    out << "%FORMAT(" << count_per_line() << "I" << width() << ")" << std::endl;
    if (elements_.empty()) {
        out << std::endl;
        return;
    }
    int lines = size()/count_per_line();
    for (int i = 0; i < lines; i++) {
        for (int j = 0; j < count_per_line(); j++) {
            int val = elements_[i*count_per_line() + j];
            out << std::right << std::setw(width()) << val;
        }
        out << std::endl;
    }
    if (size()%count_per_line() == 0)
        return;
    for (int i = 0; i < size()%count_per_line(); i++) {
        int val = elements_[lines*count_per_line() + i];
        out << std::right << std::setw(width()) << val;
    }
    out << std::endl;
}

int AmberTopIntSection::sum() const {
    return std::accumulate(elements_.begin(), elements_.end(), 0);
}

void AmberTopDoubleSection::append(const string& line) {
    int count = line.size()/width();
    for (int j = 0; j < count; j++) {
        string s = line.substr(j*width(), width());
        double num = convert_string<double>(s);
        elements_.push_back(num);
    }
}

void AmberTopDoubleSection::print(std::ostream& out) {
    out << "%FLAG " << name() << std::endl;
    out << "%FORMAT(" << count_per_line() << "E" << width() << "." <<
           decimal_places_ << ")" << std::endl;
    if (elements_.empty()) {
        out << std::endl;
        return;
    }
    int lines = size()/count_per_line();
    for (int i = 0; i < lines; i++) {
        for (int j = 0; j < count_per_line(); j++)
            out << std::right << std::setw(width()) <<
                   std::setprecision(decimal_places_) << std::scientific <<
                   elements_[i*count_per_line() + j];
        out << std::endl;
    }
    if (size()%count_per_line() == 0)
        return;
    for (int i = 0; i < size()%count_per_line(); i++)
        out << std::right << std::setw(width()) <<
               std::setprecision(decimal_places_) <<
               std::scientific <<
               elements_[(lines)*count_per_line() + i];
    out << std::endl;
}

void AmberTopStringSection::append(const string& line) {
    int count = line.size()/width();
    for (int j = 0; j < count; j++) {
        string s = line.substr(j*width(), width());
        trim(s);
        if (s.empty()) {
            return;
        }
        elements_.push_back(s);
    }
}

void AmberTopStringSection::print(std::ostream& out) {
    out << "%FLAG " << name() << std::endl;
    out << "%FORMAT(" << count_per_line() << "a" << width() << ")" << std::endl;
    if (elements_.empty()) {
        out << std::endl;
        return;
    }
    int lines = size()/count_per_line();
    for (int i = 0; i < lines; i++) {
        for (int j = 0; j < count_per_line(); j++)
            out << std::setw(width()) << std::left <<
                   elements_[i*count_per_line() + j];
        out << std::endl;
    }
    if (size()%count_per_line() == 0)
        return;
    for (int i = 0; i < size()%count_per_line(); i++)
        out << std::setw(width()) << std::left <<
               elements_[(lines)*count_per_line() + i];
    out << std::endl;
}

struct AmberTopFile::Impl {
    enum SectionType { kStringSection, kIntSection, kDoubleSection };

    // These classify the control lines of the file ("%FLAG", "%FORMAT").
    enum CardType { kFlagCard, kFormatCard, kOtherCard };

    Impl() {}

    ~Impl() {
        std::for_each(section_order.begin(), section_order.end(), DeletePtr());
    }

    AmberTopIntSection *create_int_section(const string& name,
                                           const string& format) {
        size_t i = format.find("I");
        int count_per_line = convert_string<int>(format.substr(0, i));
        int width = convert_string<int>(format.substr(i + 1));

        return create_int_section(name, count_per_line, width);
    }

    AmberTopIntSection *create_int_section(const string& name,
                                           int count_per_line, int width) {
        AmberTopIntSection *section = new AmberTopIntSection(name,
                                                             count_per_line,
                                                             width);
        int_sections.push_back(section);
        section_order.push_back(section);
        return section;
    }

    AmberTopDoubleSection *create_double_section(
            const string& name, int count_per_line, int width,
            int decimal_places) {
        AmberTopDoubleSection *section =
                new AmberTopDoubleSection(name, count_per_line, width,
                                          decimal_places);
        double_sections.push_back(section);
        section_order.push_back(section);
        return section;
    }

    AmberTopDoubleSection *create_double_section(
            const string& name, const string& format) {
        size_t i = format.find("E");
        size_t j = format.find(".");
        int count_per_line = convert_string<int>(format.substr(0, i));
        int width = convert_string<int>(format.substr(i + 1, j - i - 1));
        int decimal_places = convert_string<int>(format.substr(j + 1));

        return create_double_section(name, count_per_line, width,
                                     decimal_places);
    }

    AmberTopStringSection *create_string_section(const string& name,
                                                 int count_per_line,
                                                 int width) {
        AmberTopStringSection *section =
                new AmberTopStringSection(name, count_per_line, width);
        string_sections.push_back(section);
        section_order.push_back(section);
        return section;
    }

    AmberTopStringSection *create_string_section(const string& name,
                                                 const string& format) {
        size_t i = format.find("a");
        int count_per_line = convert_string<int>(format.substr(0, i));
        int width = convert_string<int>(format.substr(i + 1));

        return create_string_section(name, count_per_line, width);
    }

    void read(std::istream& input) {
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

    string extract_version(const string& line) {
        if (line.size() > 9)
            return line.substr(9);
        return "";
    }

    string extract_title(const string& line) {
        if (line.size() > 6) {
            string ret = string(line.substr(6));
            trim(ret);
            return ret;
        }
        return string("");
    }

    string extract_format(const string& line) {
        int left = line.find("(");
        int right = line.find(")");
        string ret = line.substr(left + 1, right - left - 1);
        trim(ret);
        return ret;
    }

    SectionType get_section_type(const string& line) {
        if (line.find('I') != string::npos) {
            return kIntSection;
        } else if (line.find('E') != string::npos) {
            return kDoubleSection;
        } else {
            return kStringSection;
        }
    }

    CardType get_card_type(const string& line) {
        if (line.find("%FLAG") == 0) {
            return kFlagCard;
        } else if (line.find("%FORMAT") == 0) {
            return kFormatCard;
        } else {
            return kOtherCard;
        }
    }

    void process_section(std::istream& input, AmberTopSection *section) {
        string line;
        while (input.peek() != '%' && getline(input, line)) {
            section->append(line);
        }
    }

    void write(std::ostream& out) const {
        out << "%VERSION " << get_version_string() << std::endl;
        for (int i = 0; i < section_order.size(); i++)
            section_order[i]->print(out);
    }

    static string get_version_string() {
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

    void sort(bool (*comp)(const AmberTopSection *lhs,
                           const AmberTopSection *rhs)) {
        std::sort(section_order.begin(),
                  section_order.end(), comp);
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
        typename T::iterator it =
                std::remove_if(vector.begin(), vector.end(), IsName(name));
        vector.erase(it, vector.end());
    }

    bool remove_section(const string& name) {
        remove_from_arr(int_sections, name);
        remove_from_arr(double_sections, name);
        remove_from_arr(string_sections, name);

        vector<AmberTopSection*>::iterator it =
                std::remove_if(section_order.begin(), section_order.end(),
                               IsName(name));
        while (it != section_order.end()) {
            delete *it;
            it = section_order.erase(it);
        }
        return true;
    }

    AmberTopIntSection *get_int_section(const string& name) {
        for (int i = 0; i < int_sections.size(); i++) {
            if (int_sections[i]->name() == name)
                return int_sections[i];
        }
        return NULL;
    }

    AmberTopDoubleSection *get_double_section(const string& name) {
        for (int i = 0; i < double_sections.size(); i++) {
            if (double_sections[i]->name() == name)
                return double_sections[i];
        }
        return NULL;
    }

    AmberTopStringSection *get_string_section(const string& name) {
        for (int i = 0; i < string_sections.size(); i++) {
            if (string_sections[i]->name() == name)
                return string_sections[i];
        }
        return NULL;
    }

    vector<AmberTopIntSection*> int_sections;
    vector<AmberTopDoubleSection*> double_sections;
    vector<AmberTopStringSection*> string_sections;
    vector<AmberTopSection*> section_order;
};

AmberTopFile::AmberTopFile() : impl_(new Impl) {
}

AmberTopFile::AmberTopFile(const File& file) : impl_(new Impl) {
    Readable::read(file);
}

AmberTopFile::~AmberTopFile() {
}

AmberTopIntSection *AmberTopFile::create_int_section(const string& name,
                                                     const string& format) {
    return impl_->create_int_section(name, format);
}

AmberTopIntSection *AmberTopFile::create_int_section(const string& name,
                                                     int count_per_line,
                                                     int width) {
    return impl_->create_int_section(name, count_per_line, width);
}

AmberTopDoubleSection *AmberTopFile::create_double_section(
        const string& name, const string& format) {
    return impl_->create_double_section(name, format);
}

AmberTopDoubleSection *AmberTopFile::create_double_section(
        const string& name, int count_per_line, int width, int decimal_places) {
    return impl_->create_double_section(name, count_per_line, width,
                                        decimal_places);
}

AmberTopStringSection *AmberTopFile::create_string_section(
        const string& name, const string& format) {
    return impl_->create_string_section(name, format);
}

AmberTopStringSection *AmberTopFile::create_string_section(
        const string& name, int count_per_line, int width) {
    return impl_->create_string_section(name, count_per_line, width);
}

bool AmberTopFile::remove_section(const string& name) {
    return impl_->remove_section(name);
}


AmberTopIntSection *AmberTopFile::get_int_section(const string& name) {
    return impl_->get_int_section(name);
}


AmberTopDoubleSection *AmberTopFile::get_double_section(const string& name) {
    return impl_->get_double_section(name);
}

AmberTopStringSection *AmberTopFile::get_string_section(const string& name) {
    return impl_->get_string_section(name);
}

void AmberTopFile::sort(bool (*comp)(const AmberTopSection *lhs,
                                     const AmberTopSection *rhs)) {
    impl_->sort(comp);
}

void AmberTopFile::read(std::istream& in) {
    impl_->read(in);
}

void AmberTopFile::write(std::ostream& out) const {
    impl_->write(out);
}

}  // namespace gmml
