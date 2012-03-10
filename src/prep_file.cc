// Author: Robert Davis

#include "gmml/internal/prep_file.h"

#include <cmath>

#include <algorithm>
#include <deque>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stack>

#include "gmml/internal/atom.h"
#include "gmml/internal/environment.h"
#include "gmml/internal/geometry.h"
#include "gmml/internal/graph.h"
#include "gmml/internal/residue.h"
#include "gmml/internal/stubs/logging.h"
#include "utilities.h"

namespace gmml {

using std::cout; //remove
using std::endl; //remove
using std::deque;
using std::map;
using std::stack;
using std::string;
using std::vector;

void PrepFile::write(std::ostream& out) const {
    out << header1_ << endl;
    out << header2_ << endl;
    for (const_iterator it = begin(); it != end(); ++it) {
        it->second->write(out);
    }
    out << "STOP" << endl;
}

void PrepFile::print(const std::string& file_name) const {
    std::ofstream out;
    out.open(file_name.c_str());
    write(out);
}

void PrepFile::add_residue(ResiduePtr residue) {
    residues_[residue->name()] = residue;
}

void PrepFile::read(const string& file_name) {
    std::ifstream stream(find_file(file_name).c_str());
    read(stream);
}

void PrepFile::read(std::istream& in) {
    getline(in, header1_);
    getline(in, header2_);
    PrepFileResidue *residue = NULL;
    while ((residue = PrepFileResidue::read_from_stream(in)) != NULL)
        residues_[residue->name()] = ResiduePtr(residue);
}


PrepFile::ResiduePtr PrepFileSet::lookup(const string& name) const {
    PrepFile::const_iterator it;
    if ((it = residues_.find(name)) == residues_.end()) {
        return PrepFile::ResiduePtr();
    }
    return it->second;
}

void PrepFileSet::load(const PrepFile& prep_file) {
    PrepFile::const_iterator it;
    for (it = prep_file.begin(); it != prep_file.end(); ++it) {
        if (residues_.find(it->first) == residues_.end()) {
            residues_.insert(*it);
        } else  {
            residues_[it->first] = it->second;
            LOG(WARNING) << "Overwriting prep file " << it->first <<
                            " in prep file set.";
        }
    }
}







struct PrepFileResidue::Impl {
    enum OtherSectionType { kSectionLoop, kSectionImproper, kSectionDone,
                            kSectionOther };

    Impl() : cutoff(0.0), dummy_atom_type("DU") {}

    ~Impl() {
        std::for_each(atoms.begin(), atoms.end(), DeletePtr());
        std::for_each(loops.begin(), loops.end(), DeletePtr());
    }

    int get_atom_index(const std::string& name) const;

    void set_header_info(std::istream&);
    void set_loops(std::istream&);
    void set_improper_dihedrals(std::istream&);

    static CoordinateType extract_coordinate_type(std::istream&);
    static DummyAtomOmission extract_dummy_omission(std::istream&);
    static DummyAtomPosition extract_dummy_position(std::istream&);
    static OutputFormat extract_output_format(std::istream&);
    static GeometryType extract_geometry_type(std::istream&);
    static OtherSectionType get_other_section_type(const string& line);

    std::string header;
    std::string file;
    std::string name;
    CoordinateType coordinate_type;
    OutputFormat output_format;
    GeometryType geometry_type;
    DummyAtomOmission dummy_atom_omission;
    std::string dummy_atom_type;
    DummyAtomPosition dummy_atom_position;
    double cutoff;
    vector<PrepFileAtom*> atoms;
    vector<ImproperDihedral> improper_dihedrals;

    vector<Loop*> loops;
};

class PrepFileResidue::CreatePrepFile {
  public:
    CreatePrepFile(PrepFileResidue *prep_residue, const Residue& residue)
            : prep_residue_(prep_residue), residue_(residue),
              visited_(residue.size(), false),
              index_in_file_(residue.size() + 3, -1),
              parent_list_(residue.size() + 3, -1) {
        Coordinate dummy1(-5.0, -5.0, -5.0);
        Coordinate dummy2(-4.5, -5.0, -5.0);
        Coordinate dummy3(-4.7, -4.2, -4.0);
        coordinate_list_.push_back(dummy1);
        coordinate_list_.push_back(dummy2);
        coordinate_list_.push_back(dummy3);
        for (int i = 0; i < residue.size(); i++)
            coordinate_list_.push_back(residue.atoms(i)->coordinate());
        create();
    }

  private:
    void create();
    void add_atom(int index, int parent1, int parent2, int parent3, int count);
    void add_loop(int index1, int index2);
    void visit_atom(int index, int parent1, int parent2, int parent3);
    int get_index_in_file(int index) const;
    void add_dummy_atom(int index);
    void set_main_chain();

    struct SubatomCompare {
        SubatomCompare(const Residue *residue) : residue_(residue) {}
        bool operator()(int atom1, int atom2) const {
            string name = residue_->name();
            if (residue_->name().size() < 2) {
                return atom1 < atom2;
            }
            bool is_l = 'a' < name[1] && name[1] < 'z';
            Element element1 = residue_->atoms(atom1)->element();
            Element element2 = residue_->atoms(atom2)->element();
            if (is_l)
                return element1 == kElementO;
            else
                return element2 == kElementO;
        };

        const Residue *residue_;
    };

    vector<int> parent_list_;
    vector<Coordinate> coordinate_list_;
    PrepFileResidue *prep_residue_;
    const Residue& residue_;
    deque<bool> visited_;
    vector<int> index_in_file_;
};


//Public impl
PrepFileResidue::PrepFileResidue(const Residue& residue) : impl_(new Impl) {
    CreatePrepFile(this, residue);
}

PrepFileResidue::PrepFileResidue() : impl_(new Impl) {
}

PrepFileResidue::~PrepFileResidue() {
}

PrepFileResidue *PrepFileResidue::read_from_stream(std::istream& in) {
    string line;

    getline(in, line);
    if (trim(line) == "STOP")
        return NULL;
    PrepFileResidue *residue = new PrepFileResidue;

    residue->impl_->header = line;

    residue->impl_->set_header_info(in);

    while (getline(in, line) && !trim(line).empty()) {
        PrepFileAtom *atom = new PrepFileAtom(line);
        residue->impl_->atoms.push_back(atom);
    }

    bool done = false;
    while (!done) {
        while (getline(in, line) && trim(line).empty()) {}
        switch (Impl::get_other_section_type(line)) {
            case Impl::kSectionLoop:
                residue->impl_->set_loops(in);
                break;
            case Impl::kSectionImproper:
                residue->impl_->set_improper_dihedrals(in);
                break;
            case Impl::kSectionDone:
                done = true;
                break;
            case Impl::kSectionOther:
                LOG(WARNING) << "Unrecognized section in prep file";
                break;
        }
    }
    return residue;
}

void PrepFileResidue::write(std::ostream& out) const {
    out << impl_->header << endl;
    out << impl_->file << endl;
    out << impl_->name;
    out << "  ";
    out << ((impl_->coordinate_type == kINT)?"INT":"XYZ");
    out << " ";
    out << ((impl_->output_format == kFormatted)?"0":"1");
    out << endl;
    out << ((impl_->geometry_type == kGeometryCorrect)?"CORRECT":"CHANGE");
    out << " ";
    out << ((impl_->dummy_atom_omission == kOmit)?"OMIT":"NOMIT");
    out << " ";
    out << impl_->dummy_atom_type;
    out << " ";
    out << ((impl_->dummy_atom_position == kPositionAll)?"ALL":"BEG");
    out << endl;  

    double total_charge = 0.0;
    for (int i = 0; i < impl_->atoms.size(); i++)
        total_charge += impl_->atoms[i]->charge;
    out << std::setw(7) << std::setprecision(4) << std::fixed << std::right <<
           total_charge << endl;

    for (int i = 0; i < impl_->atoms.size(); i++) {
        impl_->atoms[i]->write(out);
    }

    out << endl;
    if (impl_->loops.size() > 0) {
        out << "LOOP" << endl;
        for (int i = 0; i < loop_count(); i++) {
            const Loop *loop = loops(i);
            out << impl_->atoms[loop->from()]->name << " " <<
                   impl_->atoms[loop->to()]->name << endl;
        }
        out << endl;
    }
    out << "DONE" << endl;
}

int PrepFileResidue::find(const string& name) const {
    return impl_->get_atom_index(name);
}

int PrepFileResidue::atom_count() const {
    return impl_->atoms.size();
}

const PrepFileAtom *PrepFileResidue::atoms(int index) const {
    return impl_->atoms.at(index);
}

int PrepFileResidue::loop_count() const {
    return impl_->loops.size();
}

const PrepFileResidue::Loop *PrepFileResidue::loops(int index) const {
    return impl_->loops.at(index);
}

void PrepFileResidue::set_header(const std::string& header) {
    impl_->header = header;
}

std::string PrepFileResidue::header() const {
    return impl_->header;
}

std::string PrepFileResidue::file() const {
    return impl_->file;
}

std::string PrepFileResidue::name() const {
    return impl_->name;
}

PrepFileResidue::CoordinateType PrepFileResidue::coordinate_type() const {
    return impl_->coordinate_type;
}

PrepFileResidue::OutputFormat PrepFileResidue::output_format() const {
    return impl_->output_format;
}

PrepFileResidue::GeometryType PrepFileResidue::geometry_type() const {
    return impl_->geometry_type;
}

PrepFileResidue::DummyAtomOmission
PrepFileResidue::dummy_atom_omission() const {
    return impl_->dummy_atom_omission;
}

std::string PrepFileResidue::dummy_atom_type() const {
    return impl_->dummy_atom_type;
}

PrepFileResidue::DummyAtomPosition
PrepFileResidue::dummy_atom_position() const {
    return impl_->dummy_atom_position;
}

double PrepFileResidue::cutoff() const {
    return impl_->cutoff;
}


int PrepFileResidue::Impl::get_atom_index(const string& name) const {
    for (int i = 0; i < atoms.size(); i++)
        if (name == atoms[i]->name)
            return i;
    return -1;
}


void PrepFileResidue::Impl::set_header_info(std::istream& in) {
    getline(in, file);

    string line;
    getline(in, line);
    std::istringstream ss;
    ss.str(line);
    ss >> name;
    coordinate_type = extract_coordinate_type(ss);
    output_format = extract_output_format(ss);

    ss.clear();
    getline(in, line);
    ss.str(line);
    geometry_type = extract_geometry_type(ss);
    dummy_atom_omission = extract_dummy_omission(ss);
    ss >> dummy_atom_type;
    dummy_atom_position = extract_dummy_position(ss);

    getline(in, line);
    cutoff = convert_string<double>(line);
}

void PrepFileResidue::Impl::set_loops(std::istream& in) {
    string line;
    std::stringstream ss;
    while (getline(in, line) && !trim(line).empty()) {
        ss.clear();
        ss.str(line);
        string atom_names[2];
        ss >> atom_names[0] >> atom_names[1];
        int from = get_atom_index(atom_names[0]);
        int to = get_atom_index(atom_names[1]);
        if (from == -1 || to == -1) {
            //throw error here, unknown atom names
        }
        loops.push_back(new Loop(from, to));
    }
}

void PrepFileResidue::Impl::set_improper_dihedrals(std::istream& in) {
    string line;
    std::stringstream ss;
    while (getline(in, line) && !trim(line).empty()) {
        ImproperDihedral dihedral;
        ss.clear();
        ss.str(line);
        ss >> dihedral.atom_names[0] >> dihedral.atom_names[1] >>
              dihedral.atom_names[2] >> dihedral.atom_names[3];
        improper_dihedrals.push_back(dihedral);
    }
}


PrepFileResidue::CoordinateType PrepFileResidue::Impl::extract_coordinate_type(
        std::istream& in) {
    string s;
    in >> s;
    if (s == "XYZ")
        return kXYZ;
    else
        return kINT;
}

PrepFileResidue::DummyAtomOmission
PrepFileResidue::Impl::extract_dummy_omission(
        std::istream& in) {
    string s;
    in >> s;
    if (s == "NOMIT")
        return kNomit;
    else
        return kOmit;
}

PrepFileResidue::DummyAtomPosition
PrepFileResidue::Impl::extract_dummy_position(
        std::istream& in) {
    string s;
    in >> s;
    if (s == "ALL")
        return kPositionAll;
    else
        return kPositionBeg;
}

PrepFileResidue::OutputFormat PrepFileResidue::Impl::extract_output_format(
        std::istream& in) {
    int val;
    in >> val;
    if (val == 1)
        return kBinary;
    else
        return kFormatted;
}

PrepFileResidue::GeometryType PrepFileResidue::Impl::extract_geometry_type(
        std::istream& in) {
    string s;
    in >> s;
    if (s == "CHANGE")
        return kGeometryChange;
    else
        return kGeometryCorrect;
}

PrepFileResidue::Impl::OtherSectionType
PrepFileResidue::Impl::get_other_section_type(const string& line) {
    if (line == "LOOP")
        return kSectionLoop;
    else if (line == "IMPROPER")
        return kSectionImproper;
    else if (line == "DONE")
        return kSectionDone;
    return kSectionOther;
}


void PrepFileResidue::CreatePrepFile::create() {
    add_dummy_atom(0);
    add_dummy_atom(1);
    add_dummy_atom(2);

    prep_residue_->impl_->header = "";
    prep_residue_->impl_->file = "";
    prep_residue_->impl_->name = residue_.name();
    prep_residue_->impl_->coordinate_type = kINT;
    prep_residue_->impl_->output_format = kFormatted;
    prep_residue_->impl_->geometry_type = kGeometryCorrect;
    prep_residue_->impl_->dummy_atom_omission = kOmit;
    prep_residue_->impl_->dummy_atom_position = kPositionBeg;
    prep_residue_->impl_->dummy_atom_type = "DU";

    int start_index = 0;
    if (residue_.head() != -1)
        start_index = residue_.head();
    visit_atom(start_index, -1, -2, -3);

    set_main_chain();
}

void PrepFileResidue::CreatePrepFile::set_main_chain() {
    if (residue_.head() != -1 && residue_.tail() != -1) {
        int tail = residue_.tail();
        while (true) {
            int index_in_file = index_in_file_[tail];
            prep_residue_->impl_->atoms[index_in_file]->topological_type =
                    PrepFileAtom::kTopTypeM;
            if (tail == residue_.head()) {
                break;
            }
            tail = parent_list_[tail];
        }
    } else if (residue_.head() != -1) {
        prep_residue_->impl_->atoms[index_in_file_[residue_.head()]]->
                       topological_type = PrepFileAtom::kTopTypeM;
    }
}

void PrepFileResidue::CreatePrepFile::add_dummy_atom(int index) {
    PrepFileAtom *prep_atom = new PrepFileAtom;
    prep_atom->index = index + 1;
    prep_atom->name = "DUMM";
    prep_atom->type = "DU";
    prep_atom->topological_type = PrepFileAtom::kTopTypeM;
    prep_atom->bond_index = index;
    prep_atom->angle_index = index - 1;
    prep_atom->dihedral_index = index - 2;
    if (index >= 1) {
        prep_atom->bond_length = measure(coordinate_list_[index],
                                         coordinate_list_[index - 1]);
    } else {
        prep_atom->bond_length = 0.0;
    }
    if (index >= 2) {
        prep_atom->angle = to_degrees(measure(coordinate_list_[index],
                                              coordinate_list_[index - 1]));
    } else {
        prep_atom->angle = 0.0;
    }
    prep_atom->dihedral = 0.0;
    prep_atom->charge = 0.0;
    prep_residue_->impl_->atoms.push_back(prep_atom);
}

int PrepFileResidue::CreatePrepFile::get_index_in_file(int index) const {
    if (index < 0)
        return index + 3;
    return index_in_file_[index];
}

void PrepFileResidue::CreatePrepFile::add_atom(int index, int parent1,
                                               int parent2, int parent3,
                                               int count) {
    int index_in_file = prep_residue_->atom_count();
    index_in_file_[index] = index_in_file;
    PrepFileAtom *prep_atom = new PrepFileAtom;
    prep_atom->index = index_in_file + 1;

    const Atom* atom = residue_.atoms(index);
    prep_atom->name = atom->name();
    prep_atom->type = atom->type();
    prep_atom->topological_type = PrepFileAtom::get_topological_type(count);
    prep_atom->bond_index = get_index_in_file(parent1) + 1;
    prep_atom->angle_index = get_index_in_file(parent2) + 1;
    prep_atom->dihedral_index = get_index_in_file(parent3) + 1;
    prep_atom->bond_length = measure(coordinate_list_[index + 3],
                                     coordinate_list_[parent1 + 3]);
    prep_atom->angle = to_degrees(measure(coordinate_list_[index + 3],
                               coordinate_list_[parent1 + 3],
                               coordinate_list_[parent2 + 3]));
    prep_atom->dihedral = to_degrees(measure(coordinate_list_[index + 3],
                                  coordinate_list_[parent1 + 3],
                                  coordinate_list_[parent2 + 3],
                                  coordinate_list_[parent3 + 3]));
    prep_atom->charge = atom->charge();

    prep_residue_->impl_->atoms.push_back(prep_atom);
}

void PrepFileResidue::CreatePrepFile::add_loop(int index1, int index2) {
    prep_residue_->impl_->loops.push_back(new Loop(index1, index2));
}

void PrepFileResidue::CreatePrepFile::visit_atom(int index, int parent1,
                                                 int parent2, int parent3) {
    if (visited_[index])
        return;
    visited_[index] = true;
    int count = 0;
    const vector<size_t>& bonds = residue_.bonds(index);
    vector<int> sub_atoms;
    vector<int> loops_to_add;
    for (int i = 0; i < bonds.size(); i++) {
        if (visited_[bonds[i]] && bonds[i] != parent1) {
            loops_to_add.push_back(bonds[i]);
            //add_loop(index_in_file_[index] + 1, index_in_file_[bonds[i]] + 1);

            //cout << "loop from " << index << " (" <<
            //        residue_.atoms(index)->name() << ") to " <<
            //        bonds[i] << " (" << residue_.atoms(bonds[i])->name() <<
            //        ")" << endl;
        } else if (!visited_[bonds[i]]) {
            count++;
            sub_atoms.push_back(bonds[i]);
        }
    }

    std::sort(sub_atoms.begin(), sub_atoms.end(), SubatomCompare(&residue_));
    for (int i = 0; i < sub_atoms.size(); i++) {
        if (residue_.tail() == sub_atoms[i]) {
            std::swap(sub_atoms[i], sub_atoms[sub_atoms.size() - 1]);
            break;
        }
    }

  
    parent_list_[index] = parent1;
    add_atom(index, parent1, parent2, parent3, count);

    //cout << "Atom " << index << " has subatom count " << count << " " <<
    //        measure(coordinate_list_[index + 3], coordinate_list_[parent1 + 3]) << endl;
    for (int i = 0; i < sub_atoms.size(); i++) {
        visit_atom(sub_atoms[i], index, parent1, parent2);
    }

    for (int i = 0; i < loops_to_add.size(); i++) {
        add_loop(index_in_file_[index],
                 index_in_file_[loops_to_add[i]]);
    }
}

struct PrepFileAtom::Impl {
    TopologicalType extract_topological_type(std::istream& in) const;
};

PrepFileAtom::PrepFileAtom(const string& line) : impl_(new Impl) {
    std::stringstream ss(line);
    ss >> index >> name >> type;
    topological_type = impl_->extract_topological_type(ss);
    ss >> bond_index >> angle_index >> dihedral_index >>
          bond_length >> angle >> dihedral >>
          charge;
}

PrepFileAtom::~PrepFileAtom() {
}

PrepFileAtom::TopologicalType PrepFileAtom::get_topological_type(int count) {
    return static_cast<TopologicalType>(count);
}

string PrepFileAtom::get_top_type_string(TopologicalType type) {
    switch (type) {
        case kTopTypeE:
            return "E";
        case kTopTypeS:
            return "S";
        case kTopTypeB:
            return "B";
        case kTopType3:
            return "3";
        case kTopType4:
            return "4";
        case kTopTypeM:
            return "M";
        default:
            throw std::invalid_argument("Invalid topological type.");
    }
}

void PrepFileAtom::write(std::ostream& out) const {
    out << std::setw(2) << std::right << index;
    out << " ";
    out << std::setw(4) << std::left << name;
    out << " ";
    out << std::setw(2) << type;
    out << std::setw(3) << std::right << get_top_type_string(topological_type);
    out << " ";
    out << std::setw(2) << std::right << bond_index;
    out << " ";
    out << std::setw(2) << std::right << angle_index;
    out << " ";
    out << std::setw(2) << std::right << dihedral_index;
    out << std::setw(7) << std::right << std::setprecision(3) <<
           std::fixed << bond_length;
    out << " ";
    out << std::setw(7) << std::right << std::setprecision(3) <<
           std::fixed << angle;
    out << " ";
    out << std::setw(7) << std::right << std::setprecision(2) <<
           std::fixed << dihedral;
    out << " ";
    out << std::setw(10) << std::right << std::setprecision(4) <<
           std::fixed << charge;
    out << endl;
}

PrepFileAtom::TopologicalType PrepFileAtom::Impl::extract_topological_type(
        std::istream& in) const {
    string s;
    in >> s;
    if (s == "M")
        return kTopTypeM;
    else if (s == "S")
        return kTopTypeS;
    else if (s == "B")
        return kTopTypeB;
    else if (s == "E")
        return kTopTypeE;
    else
        return kTopType3;
}







namespace {

void set_parent_list(const PrepFileResidue& residue,
                     vector<int>& parent_list) {
    stack<int> st;
    for (int i = residue.atom_count() - 1; i >= 0; i--) {
        size_t stack_size = st.size();
        switch (residue.atoms(i)->topological_type) {
            case PrepFileAtom::kTopTypeM:
                while (!st.empty()) {
                    parent_list[st.top()] = i;
                    st.pop();
                }
                break;
            case PrepFileAtom::kTopTypeS:
                if (stack_size > 0) {
                    parent_list[st.top()] = i;
                    st.pop();
                }
                break;
            case PrepFileAtom::kTopTypeB:
                for (size_t j = 0; j < std::min((size_t)2, stack_size); j++) {
                    parent_list[st.top()] = i;
                    st.pop();
                }
                break;
            case PrepFileAtom::kTopTypeE:
                break;
            case PrepFileAtom::kTopType3:
                for (size_t j = 0; j < std::min((size_t)3, stack_size); j++) {
                    parent_list[st.top()] = i;
                    st.pop();
                }
                break;
            default:
                LOG(WARNING) << "Unrecognized topological type";
                break;
        }
        st.push(i);
    }
}

void set_dummy_coordinates(const PrepFileResidue& residue,
                           vector<Coordinate*>& coordinates) {
    coordinates[0] = new Coordinate(0.0, 0.0, 0.0);

    double dist1 = residue.atoms(1)->bond_length;
    coordinates[1] = new Coordinate(dist1, 0.0, 0.0);

    double dist2 = residue.atoms(2)->bond_length;
    double angle = residue.atoms(2)->angle;
    coordinates[2] = new Coordinate(dist1 - cos(to_radians(angle))*dist2,
                                    sin(to_radians(angle))*dist2,
                                    0.0);
}

}  // namespace

Residue *BuildPrepFileResidue::operator()(
        const PrepFileResidue& prep_file_residue) const {
    //const vector<PrepFileAtom*> atom_list = prep_file_residue.atoms;
    // The index of the parent of each atom
    vector<int> parent_list(prep_file_residue.atom_count());
    set_parent_list(prep_file_residue, parent_list);

    vector<Coordinate*> coordinates(prep_file_residue.atom_count());
    set_dummy_coordinates(prep_file_residue, coordinates);

    Graph *bonds = new Graph(prep_file_residue.atom_count() - 3);
    for (size_t i = 3; i < prep_file_residue.atom_count(); i++) {
        int parent1 = parent_list[i];
        // Grandparent of i
        int parent2 = parent_list[parent1];
        // Great-grandparent of i
        int parent3 = parent_list[parent2];
        if (parent1 > 2)
            bonds->add_edge(i - 3, parent1 - 3);

        coordinates[i] = new Coordinate(
                calculate_point(*coordinates[parent3], *coordinates[parent2],
                                *coordinates[parent1],
                                to_radians(prep_file_residue.atoms(i)->angle),
                               to_radians(prep_file_residue.atoms(i)->dihedral),
                                prep_file_residue.atoms(i)->bond_length)
        );
    }
    for (size_t i = 0; i < prep_file_residue.loop_count(); i++) {
        const PrepFileResidue::Loop *loop = prep_file_residue.loops(i);
        bonds->add_edge(loop->from() - 3, loop->to() - 3);
    }

    vector<Atom*> *atoms = new vector<Atom*>;
    atoms->reserve(prep_file_residue.atom_count() - 3);

    int head = -1;
    int tail = -1;
    for (size_t i = 3; i < prep_file_residue.atom_count(); i++) {
        const PrepFileAtom *prep_atom = prep_file_residue.atoms(i);
        Element element = get_element_by_char(prep_atom->name[0]);
        Atom *atom = new Atom(element, *coordinates[i], prep_atom->name,
                              prep_atom->type, prep_atom->charge);
        atoms->push_back(atom);

        if (prep_atom->topological_type == PrepFileAtom::kTopTypeM) {
            if (head == -1)
                head = i - 3;
            if (head != i - 3)
                tail = i - 3;
        }
    }

    std::for_each(coordinates.begin(), coordinates.end(), DeletePtr());

    Residue *result = new Residue(prep_file_residue.name(), atoms, bonds);
    result->set_head(head);
    result->set_tail(tail);

    return result;
}

}  // namespace gmml
