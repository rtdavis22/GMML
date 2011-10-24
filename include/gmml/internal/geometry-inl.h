#ifndef GEOMETRY_INL_H
#define GEOMETRY_INL_H

namespace gmml
{

inline RotationMatrix::RotationMatrix(const Coordinate& point,
                                      Vector<3> direction,
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

} //namespace gmml

#endif
