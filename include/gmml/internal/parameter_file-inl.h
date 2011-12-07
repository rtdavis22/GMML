#ifndef GMML_INTERNAL_PARAMETER_FILE_INL_H_
#define GMML_INTERNAL_PARAMETER_FILE_INL_H

namespace gmml {

inline bool ParameterFileImproperDihedral::matches(
        const std::vector<std::string>& rhs) const {
    std::vector<std::string> rhs_types(types);
    std::sort(rhs_types.begin(), rhs_types.end());
    //no Xs
    if (types.size() == 3)
	return rhs[0] == types[0] && rhs[1] == types[1] &&
	       rhs[2] == types[2];
    //one X
    else if (types.size() == 2)
	return (rhs[0] == types[0] && rhs[1] == types[1]) ||
	       (rhs[1] == types[0] && rhs[2] == types[1]) ||
	       (rhs[0] == types[0] && rhs[2] == types[1]);
    //two Xs
    else if (types.size() == 1)
	return rhs[0] == types[0] || rhs[1] == types[0] || 
	       rhs[2] == types[0];
    //all Xs
    else 
        return true;
}

inline std::pair<ParameterFileDihedralTerm,bool>
ImproperDihedralCollection::lookup(
        const std::string& center, 
        const std::vector<std::string>& types) const {
    ParameterFileDihedralTerm term;
    const_iterator it = dihedrals_.find(center);
    if (it != dihedrals_.end()) {
        for (int i = it->second.size() - 1; i >= 0; i--) {
            if (it->second[i]->matches(types))
                return std::make_pair(it->second[i]->term, true);
        }
    }
    return std::make_pair(term, false);
}

inline void ImproperDihedralCollection::append(
        const ImproperDihedralCollection& coll) {
    std::vector<ParameterFileImproperDihedral*>::const_iterator it;
    const_iterator coll_it;
    for (coll_it = coll.begin(); coll_it != coll.end(); ++coll_it) {
        for (it = coll_it->second.begin(); it != coll_it->second.end(); ++it) {
            dihedrals_[coll_it->first].push_back(
                new ParameterFileImproperDihedral(**it)
            );
        }
    }
}

inline ParameterFile::~ParameterFile() {
    for (AtomTypeMap::iterator it = atom_types_.begin();
            it != atom_types_.end(); ++it)
        delete it->second;
    std::for_each(bonds_.begin(), bonds_.end(), DeletePtr());
    std::for_each(angles_.begin(), angles_.end(), DeletePtr());
    std::for_each(dihedrals_.begin(), dihedrals_.end(), DeletePtr());
    std::for_each(generic_dihedrals_.begin(), generic_dihedrals_.end(),
                  DeletePtr());
}

} // namespace gmml

#endif  // GMML_INTERNAL_PARAMETER_FILE_INL_H_
