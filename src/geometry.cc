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

#include "gmml/internal/geometry.h"

namespace gmml {

// Reference: http://inside.mines.edu/~gmurray/ArbitraryAxisRotation/ArbitraryAxisRotation.html
RotationMatrix::RotationMatrix(const Coordinate& point, Vector<3> direction,
                               double radians) {
    direction.normalize();
    double u = direction[0];
    double v = direction[1];
    double w = direction[2];

    double a = point.x;
    double b = point.y;
    double c = point.z;

    double u2 = u*u;
    double v2 = v*v;
    double w2 = w*w;
    double cos_t = cos(radians);
    double sin_t = sin(radians);

    matrix_[0][3] = a*(v2 + w2) - u*(b*v + c*w) +
                 (u*(b*v + c*w) - a*(v2 + w2))*cos_t + (b*w - c*v)*sin_t;

    matrix_[1][3] = b*(u2 + w2) - v*(a*u + c*w) +
                 (v*(a*u + c*w) - b*(u2 + w2))*cos_t + (c*u - a*w)*sin_t;

    matrix_[2][3] = c*(u2 + v2) - w*(a*u + b*v) +
                  (w*(a*u + b*v) - c*(u2 + v2))*cos_t + (a*v - b*u)*sin_t;

    matrix_[0][0] = u2 + (v2 + w2)*cos_t;
    matrix_[0][1] = u*v*(1 - cos_t) - w*sin_t;
    matrix_[0][2] = u*w*(1 - cos_t) + v*sin_t;

    matrix_[1][0] = u*v*(1 - cos_t) + w*sin_t;
    matrix_[1][1] = v2 + (u2 + w2)*cos_t;
    matrix_[1][2] = v*w*(1 - cos_t) - u*sin_t;

    matrix_[2][0] = u*w*(1 - cos_t) - v*sin_t;
    matrix_[2][1] = v*w*(1 - cos_t) + u*sin_t;
    matrix_[2][2] = w2 + (u2 + v2)*cos_t;
}

// See the NeRF algorithm in 
// Jerod Parsons, J. Bradley Holmes, J. Maurice Rojas, Jerry Tsai, Charlie E. M. Strauss (2005). "Practical conversion from torsion space to Cartesian space for in silico protein synthesis". Journal of Computational Chemistry 26 (10): 1063–1068
// TODO: fix this up.
Coordinate calculate_point(const Coordinate& a, const Coordinate& b,
                           const Coordinate& c, double angle,
                           double dihedral, double distance) {
    dihedral = kPi - dihedral;

    Matrix<3, 4> m;
    Vector<3> v1(b, a);
    Vector<3> v2(c, b);

    v1.normalize();
    v2.normalize();

    // if v1 is a multiple of v2, a, b, and c are collinear
    // change to numeric utilities
    if (std::abs(v1[0] + v2[0]) < 0.001 &&
            std::abs(v1[1] + v2[1]) < 0.001 &&
            std::abs(v1[2] + v2[2]) < 0.001) {
        // use an arbitary dummy coordinate
        Coordinate a_p(a.x+10, a.y-1, a.z+3);
        v1 = Vector<3>(b, a_p);
    }

    m.set_column(1, cross(v1, v2).normalize());
    m.set_column(2, v2);
    m.set_column(0, cross(m.column(1), m.column(2)));
    m.set_column(3, Vector<3>(c));
    Vector<4> v;
    v[0] = distance*sin(angle)*cos(dihedral);
    v[1] = distance*sin(angle)*sin(dihedral);
    v[2] = distance*cos(angle);
    v[3] = 1;

    Vector<3> ret = m*v;

    return Coordinate(ret[0], ret[1], ret[2]);
}

}  // namespace gmml
