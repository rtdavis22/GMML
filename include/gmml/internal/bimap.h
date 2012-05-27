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

#ifndef GMML_INTERNAL_BIMAP_H_
#define GMML_INTERNAL_BIMAP_H_

#include <map>
#include <utility>

namespace gmml {
namespace internal {

template <typename T1, typename T2>
class Bimap {
  public:
    std::pair<T2, bool> get(const T1& key) const {
        typename std::map<T1, T2>::const_iterator it;
        if ((it = forward_map_.find(key)) != forward_map_.end())
            return std::make_pair(it->second, true);
        else
            return std::make_pair(T2(), false);
    }

    std::pair<T1, bool> get_inverse(const T2& value) const {
        typename std::map<T2, T1>::const_iterator it;
        if ((it = inverse_map_.find(value)) != inverse_map_.end())
            return std::make_pair(it->second, true);
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
