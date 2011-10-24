#include "gmml/gmml.h"

#include <iostream>
#include <iterator>
#include <algorithm>
#include <string>

using namespace gmml;

using namespace std;

struct TestClass {

    ~TestClass() { cout << "dtor" << endl; }

};

/*
void print_node(ParsedResidue *residue, int level);

void print_tree(Tree<ParsedResidue*> *tree) {
    Tree<ParsedResidue*>::iterator it = tree->root();
    int level = 0;
    while (1) {
        
    }
}
*/

int main() {

enable_warnings();

PrepFile prep_file("/home/robert/Documents/gmml/dat/prep_files/Glycam_06.prep");

PrepFileSet prep_files;
prep_files.load(prep_file);

ParameterFile param_file("/home/robert/Documents/gmml/dat/param_files/Glycam_06g.dat");

load_prep_file("/home/robert/Documents/gmml/dat/prep_files/Glycam_06.prep");

ParameterFileSet parm_set;

parm_set.load("/home/robert/Documents/gmml/dat/param_files/parm99.dat.mod");


parm_set.load(param_file);


Residue *residue = build_prep_file("OME");

/*
Structure s;
s.append(residue);
s.append(build_prep_file("0MA"));
AmberTopBuilder b(parm_set);
AmberTopFile *top_file = b.build(s);
top_file->print("temp.top");
*/

//cout << endl;

//s.bonds()->print();

//std::copy(sections.begin(), sections.end(),
//          ostream_iterator<string>(cout, "\n"));

enable_warnings();

SequenceParser *parser = new GlycamParser;
tree<ParsedResidue*> *residue_tree = parser->parse("DGalpb1-4DGlcpNAcb1-6[DGalpb1-4DGlcpNAcb1-2]DManpa1-6[DGalpb1-4DGlcpNAcb1-2DManpa1-3]DManpb1-4DGlcpNAcb1-OH");

tree<ParsedResidue*>::iterator it;
for (it = residue_tree->begin(); it != residue_tree->end(); ++it) {
    int depth = residue_tree->depth(it);
    for (int i = 0; i < depth; i++)
        cout << " ";
    cout << (*it)->anomeric_carbon << "-" << (*it)->oxygen_position <<
            " " << (*it)->name << endl;
}

GlycamCodeSet code_set;

/*
  tree<string> tr;
   tree<string>::iterator top, one, two, loc, banana;
   
   top=tr.begin();
   one=tr.insert(top, "one");
   two=tr.append_child(one, "two");
   tr.append_child(two, "apple");
   banana=tr.append_child(two, "banana");
   tr.append_child(banana,"cherry");
   tr.append_child(two, "peach");
   tr.append_child(one,"three");
*/
/*
PreOrderIterator<ParsedResidue*> it(tree->root());
while (!it.at_end()) {
    //ParsedResidue *res = *it; // *it.ptr;
    TreeNode<ParsedResidue*> *node = 
    cout << (*it)->name_ << endl;;
    ++it;
}
*/

//Residue<TypedAtom> *residue = build_prep_file(prep_files["OME"]);

//Structure<TypedAtom> structure;

//structure.set_dihedral(1, 1, 1, 1, 30.0);

//PdbFileBuilder b;
//PdbFile *pdb_file = b.build(structure);

//error("asd");
/*
PreOrderIterator<int> ioi = t.root();
while(!ioi.at_end()) {
     cout << *ioi << endl;
     ++ioi;
}

cout << endl;
ioi = c->root();
while(!ioi.at_end()) {
     cout << *ioi << endl;
     ++ioi;
}
*/

}
