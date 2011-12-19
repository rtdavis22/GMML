#ifndef GMML_INTERNAL_BIMAP_H_
#define GMML_INTERNAL_BIMAP_H_

#include <map>
#include <utility>

namespace gmml {
namespace internal {

template <typename T1, typename T2>
class Bimap {
  public:
    template<typename T>
    std::pair<T2, bool> get(const T& key) {
        typename std::map<T2, T1>::iterator it;
        if ((it = forward_map_.find(key)) != forward_map_.end())
            return std::make_pair(*it, true);
        else
            return std::make_pair(T2(), false);
    }

    template<typename T>
    std::pair<T1, bool> get_inverse(const T& value) {
        typename std::map<T1, T2>::iterator it;
        if ((it = inverse_map_.find(value)) != inverse_map_.end())
            return std::make_pair(*it, true);
        else
            return std::make_pair(T1(), false);
    }

    template<typename T, typename U>
    void put(const T& key, const U& value) {
        forward_map_[key] = value;
        inverse_map_[value] = key;
    }

  private:
    std::map<T1, T2> forward_map_;
    std::map<T2, T1> inverse_map_;
};

}  // namespace internal

using internal::Bimap;

}  // namespace gmml

#endif  // GMML_INTERNAL_BIMAP_H_
