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

// The file contains classes and functions for 3D coordinates and
// statically-sized vectors and matrices.

#ifndef GMML_INTERNAL_GEOMETRY_H_
#define GMML_INTERNAL_GEOMETRY_H_

#include <cmath>

#include <algorithm>
#include <functional>
#include <numeric>
#include <stdexcept>
#include <vector>

#include "gmml/internal/stubs/common.h"

namespace gmml {

// A coordinate in 3 dimensions.
struct Coordinate {
    Coordinate() : x(0.0), y(0.0), z(0.0) {}
    Coordinate(double x, double y, double z) : x(x), y(y), z(z) {}
    void translate(double x, double y, double z) {
        this->x += x;
        this->y += y;
        this->z += z;
    }

    double x, y, z;
};

// A base class for statically-sized vectors with numeric elements.
// TODO: There is a design flaw with VectorBase, Vector, and the specialized
// Vector<3, T>. It was intended for Vector<3, T> to include cross-product and
// potentially other things that don't have a useful analog in n dimensions,
// and for VectorBase not to be part of the public interface. There probably
// shouldn't be a VectorBase class at all. We may want to rethink this design in
// the future or use an external library like Boost. At any rate, these classes
// are currently used so little that it isn't worth the time.
template<size_t size_, typename T = double>
class VectorBase {
  public:
    VectorBase() : vector_(size_) {}

    // Dot product.
    template<typename U>
    T dot(const VectorBase<size_, U>& rhs) const {
        return std::inner_product(vector_.begin(), vector_.end(),
                                  rhs.vector_.begin(), T());
    }

    // Euclidean norm.
    double norm() const { return sqrt(dot(*this)); }

    const T& operator[](size_t i) const { return vector_[i]; }
    T& operator[](size_t i) { return vector_[i]; }

    VectorBase& operator*=(double scalar) {
        std::transform(vector_.begin(), vector_.end(), vector_.begin(),
                       std::bind1st(std::multiplies<T>(), scalar));
        return *this;
    }
    VectorBase& operator/=(double scalar) {
        std::transform(vector_.begin(), vector_.end(), vector_.begin(),
                       std::bind2nd(std::divides<T>(), scalar));
        return *this;
    }
    template<typename U>
    VectorBase& operator+=(const VectorBase<size_, U>& rhs) {
        std::transform(vector_.begin(), vector_.end(), rhs.vector_.begin(),
                       vector_.begin(), std::plus<T>());
        return *this;
    }
    VectorBase& operator+=(double rhs) {
        std::transform(vector_.begin(), vector_.end(), vector_.begin(),
                       std::bind1st(std::plus<T>(), rhs));
        return *this;
    }
    template<typename U>
    VectorBase& operator-=(const VectorBase<size_, U>& rhs) {
        std::transform(vector_.begin(), vector_.end(), rhs.vector_.begin(),
                       vector_.begin(), std::minus<T>());
        return *this;
    }
    VectorBase& operator-=(double rhs) {
        std:transform(vector_.begin(), vector_.end(), vector_.begin(),
                      std::bind1st(std::minus<T>(), rhs));
        return *this;
    }

    // Make the vector a unit vector in the same direction.
    VectorBase& normalize() {
        double norm_;
        if ((norm_ = norm()) != 0.0)
            operator/=(norm_);
        return *this;
    }

    template<size_t, typename>
    friend class VectorBase;

  protected:
    // The underlying vector.
    std::vector<T> vector_;
};

template<size_t size, typename T>
inline const VectorBase<size, T> operator*(const VectorBase<size, T>& lhs,
                                           double rhs) {
    VectorBase<size, T> t(lhs);
    t*=rhs;
    return t;
}

template<size_t size, typename T>
inline const VectorBase<size, T> operator/(const VectorBase<size, T>& lhs,
                                           double rhs) {
    VectorBase<size, T> t(lhs);
    t/=rhs;
    return t;
}

template<size_t size, typename T, typename U>
inline const VectorBase<size, T> operator+(const VectorBase<size, T>& lhs,
                                           const VectorBase<size, U>& rhs) {
    VectorBase<size, T> t(lhs);
    lhs += rhs;
    return t;
}

template<size_t size, typename T>
inline const VectorBase<size, T> operator+(const VectorBase<size, T>& lhs,
                                           double rhs) {
    VectorBase<size, T> t(lhs);
    t += rhs;
    return t;
}

template<size_t size, typename T, typename U>
inline const VectorBase<size, T> operator-(const VectorBase<size, T>& lhs,
                                           const VectorBase<size, U>& rhs) {
    VectorBase<size, T> t(lhs);
    t -= rhs;
    return t;
}

template<size_t size, typename T>
inline const VectorBase<size, T> operator-(const VectorBase<size, T>& lhs,
                                           double rhs) {
    VectorBase<size, T> t(lhs);
    t -= rhs;
    return t;
}

template<size_t size_, typename T = double>
class Vector : public VectorBase<size_, T> {};

// This is mostly for the cross product but can also include optimizations
// to the VectorBase operations.
template<typename T>
class Vector<3, T> : public VectorBase<3, T> {
  public:
    Vector() {}
    explicit Vector(const Coordinate& c) {
        this->vector_[0] = c.x;
        this->vector_[1] = c.y;
        this->vector_[2] = c.z;
    }
    Vector(const Coordinate& c1, const Coordinate& c2) {
        this->vector_[0] = c2.x - c1.x;
        this->vector_[1] = c2.y - c1.y;
        this->vector_[2] = c2.z - c1.z;
    }
    Vector(T x, T y, T z) {
        this->vector_[0] = x;
        this->vector_[1] = y;
        this->vector_[2] = z;
    }

    // Cross-product.
    template<typename U>
    Vector<3, T>& cross(const Vector<3, U>& rhs) {
        std::vector<T> t(this->vector_);
        this->vector_[0] = t[1]*rhs[2] - t[2]*rhs[1];
        this->vector_[1] = t[2]*rhs[0] - t[0]*rhs[2];
        this->vector_[2] = t[0]*rhs[1] - t[1]*rhs[0];
        return *this;
    }
};

template<typename T, typename U>
inline Vector<3, T> cross(const Vector<3, T>& vec1, const Vector<3, U>& vec2) {
    Vector<3, T> new_vec(vec1);
    return new_vec.cross(vec2);
}

// A statically sized matrix with numeric elements.
template<size_t dim1, size_t dim2, typename T = double>
class Matrix {
  public:
    Matrix() : matrix_(dim1) {}

    // Extract an entry.
    T& operator()(size_t i, size_t j) { return matrix_[i][j]; }
    const T& operator()(size_t i, size_t j) const { return matrix_[i][j]; }

    // Extract a column.
    Vector<dim1> column(size_t index) const {
        Vector<dim1> column_;
        for (size_t i = 0; i < dim1; i++)
           column_[i] = matrix_[i][index];
        return column_;
    }

    template<typename U>
    void set_column(size_t index, const VectorBase<dim1, U>& column) {
        for (size_t i = 0; i < dim1; i++)
            matrix_[i][index] = column[i];
    }

    // Extract a row.
    Vector<dim2>& operator[](size_t i) { return matrix_[i]; }
    const Vector<dim2>& operator[](size_t i) const { return matrix_[i]; }

  protected:
    std::vector<Vector<dim2, T> > matrix_;
};

template<size_t dim1, size_t dim2, size_t dim3, typename T, typename U>
inline Matrix<dim1, dim2, T> operator*(const Matrix<dim1, dim3, T>& lhs,
                                       const Matrix<dim3, dim2, U>& rhs) {
    Matrix<dim1, dim2, T> matrix;
    for (size_t i = 0; i < dim1; i++)
        for (size_t j = 0; j < dim2; j++)
            matrix[i][j] = lhs[i].dot(rhs.column(j));
    return matrix;
}

template<size_t dim1, size_t dim2, typename T, typename U>
inline Vector<dim1, U> operator*(const Matrix<dim1, dim2, T>& lhs,
                                 const Vector<dim2, U>& rhs) {
    Vector<dim1, U> vec;
    for (size_t i = 0; i < dim1; i++)
        vec[i] = lhs[i].dot(rhs);
    return vec;
}

// A rotation matrix that can be applied to Coordinates.
class RotationMatrix : Matrix<3, 4, double> {
  public:
    // Create a matrix for rotating @radians about a coordinate in a
    // direction
    RotationMatrix(const Coordinate& coordinate, const Vector<3> direction,
                   double radians);

    void apply(Coordinate& coordinate) const {
        Vector<4> vec;
        vec[0] = coordinate.x;
        vec[1] = coordinate.y;
        vec[2] = coordinate.z;
        vec[3] = 1.0;
        Vector<3> result = (*this)*vec;
        coordinate.x = result[0];
        coordinate.y = result[1];
        coordinate.z = result[2];
    }
};

// Measure the angle between two vectors.
template<size_t size_, typename T, typename U>
inline T measure(const Vector<size_, T>& vec1, const Vector<size_, U>& vec2) {
    T norm1 = vec1.norm();
    U norm2 = vec2.norm();
    if (norm1 != 0.0 && norm2 != 0.0)
        return acos(vec1.dot(vec2)/norm1/norm2);
    throw std::invalid_argument("measure(Vector&, Vector&)");
}

// Measure the distance between two coordinates.
inline double measure(const Coordinate& c1, const Coordinate& c2) {
    Vector<3, double> vec1(c1.x, c1.y, c1.z);
    Vector<3, double> vec2(c2.x, c2.y, c2.z);
    return (vec1 - vec2).norm();
}

// Measure the angle formed by the coordinates. The return value is in radians.
inline double measure(const Coordinate& c1, const Coordinate& c2,
                      const Coordinate& c3) {
    Vector<3> b1(c2, c1);
    Vector<3> b2(c2, c3);
    return acos(b1.dot(b2)/b1.norm()/b2.norm());
}

// Measure the dihedral of the coordinates. The return value is in radians.
inline double measure(const Coordinate& c1, const Coordinate& c2,
                      const Coordinate& c3, const Coordinate& c4) {
    Vector<3> b1(c1, c2);
    Vector<3> b2(c2, c3);
    Vector<3> b3(c3, c4);
    Vector<3> b2xb3 = cross(b2, b3);
    return atan2((b1*b2.norm()).dot(b2xb3), cross(b1, b2).dot(b2xb3));
}

// Set the third coordinate so that the angle formed by the coordinates has
// the given measure in radians. This is here mostly for testing purposes.
inline void set_angle(const Coordinate *c1, const Coordinate *c2,
                      Coordinate *c3, double radians) {
    double cur_angle = measure(*c1, *c2, *c3);
    RotationMatrix matrix(*c2, Vector<3>(*c1, *c2).cross(Vector<3>(*c3, *c2)),
                          radians - cur_angle);
    matrix.apply(*c3);
}

// Set the fourth coordinate so that the dihedral has the given measure in
// radians. This is here mostly for testing purposes.
inline void set_dihedral(const Coordinate *c1, const Coordinate *c2,
                         const Coordinate *c3, Coordinate *c4, double radians) {
    double cur_angle = measure(*c1, *c2, *c3, *c4);
    RotationMatrix matrix(*c2, Vector<3>(*c3, *c2), cur_angle - radians);
    matrix.apply(*c4);
}

// Calculate the unique coordinate d so that the following conditions are
// satisfied:
// 1. d is |distance| units from |c|.
// 2. The angle d-c-b is |angle| radians.
// 3. The dihedral d-c-b-a is |dihedral| radians.
// Note that the dihedral angle vanishes (doesn't exist) when the input
// coordinates are collinear. In this case, there are infinitely many
// possibilities for the fourth coordinate (and they lie on a circle).
// I don't see a good, deterministic way of selecting an element from this set
// so, officially, the fourth coordinate can be any of them. Also, the third
// condition cannot be met in this case.
Coordinate calculate_point(const Coordinate& a, const Coordinate& b,
                           const Coordinate& c, double angle,
                           double dihedral, double distance);

}  // namespace gmml

#endif  // GMML_INTERNAL_GEOMETRY_H_
