/* 
 * Copyright 2017 Jeffrey Davis
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 * 
 * File:   Point_2D.h
 * Author: Jeffrey Davis
 * Version: 1.0
 */

#ifndef POINT_2D_H
#define POINT_2D_H
#include "Vector_2D.h"
#include <cmath>

using namespace std;

namespace VCAD_lib {
    
    class Vector_2D;
    
    /*
     * Represents a 2 Dimensional point.
     * 
     * Exception safety: no-throw guarantee
     */
    class Point_2D {
    public:
        typedef double Measurement;
        typedef double Angle_Meas;
        Point_2D(const Measurement x_val, const Measurement y_val);
        Measurement get_x() const { return x; }
        Measurement get_y() const { return y; }
        Measurement get_r() const;
        Measurement get_r(const Point_2D& origin) const;
        Angle_Meas get_theta() const;
        Angle_Meas get_theta(const Point_2D& origin) const;
        Point_2D& rotate(const Angle_Meas angle);
        Point_2D& rotate(const Angle_Meas angle, const Point_2D& origin);
        Point_2D& scale(const Measurement x_scalar, const Measurement y_scalar);
        Point_2D& scale(const Measurement x_scalar, const Measurement y_scalar, const Point_2D& origin);
        Point_2D& translate(const Measurement x_val, const Measurement y_val);
        Point_2D& translate(const Vector_2D& v);
        Point_2D& move(const Point_2D& new_origin, const Vector_2D& axis,
                const bool is_x_axis, const Point_2D& ref_origin=Point_2D(0,0));
        Point_2D& operator+=(const Vector_2D&);
        Point_2D& operator-=(const Vector_2D&);
        Point_2D& operator*=(const Measurement);
    private:
        static const Angle_Meas two_pi;
        Measurement x;
        Measurement y;
    };
    
    // exception safety: no-throw guarantee
    const Point_2D operator+(const Point_2D&, const Vector_2D&);
    // exception safety: no-throw guarantee
    const Point_2D operator-(const Point_2D&, const Vector_2D&);
    // exception safety: no-throw guarantee
    const Point_2D operator*(const Point_2D&, const Point_2D::Measurement);
    
    /* 
     * determines if two points are within a rounding error of each other
     */
    const bool is_equal(const Point_2D& point1, const Point_2D& point2, 
            const Point_2D::Measurement precision);
    
    // exception safety: strong guarantee.  Throws domain_error if r is negative
    const Point_2D polar_point(const Point_2D::Measurement r, const Point_2D::Angle_Meas theta);

    // exception safety: strong guarantee.  Throws domain_error if r is negative
    const Point_2D polar_point(const Point_2D::Measurement r, const Point_2D::Angle_Meas theta, const Point_2D& origin);
    
    struct Vector_2D_idata
    {
        int num; // number of intersection points found (0,1, or 2))
        Point_2D p1;
        Point_2D p2;
        Vector_2D_idata();
    };
}
#endif /* POINT_2D_H */

