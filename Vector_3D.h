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
 * File:   Vector_3D.h
 * Author: Jeffrey Davis
 */

#ifndef VECTOR_3D_H
#define VECTOR_3D_H
#include <iostream>
#include <cmath>
#include "Point_3D.h"

using namespace std;

namespace VCAD_lib
{
    class Point_3D;
    
    // exception safety: no throw unless otherwise indicated
    class Vector_3D {
    public:
        typedef double Measurement;
        typedef double Angle_Meas;
        Vector_3D(const Measurement x_val, const Measurement y_val, const Measurement z_val);
        Vector_3D(const Point_3D& end_point);
        Vector_3D(const Point_3D& start_point, const Point_3D& end_point);
//        Vector_3D(const Vector_3D& orig) : x(orig.x), y(orig.y), z(orig.z) {}
        Measurement get_x() const { return x; }
        Measurement get_y() const { return y; }
        Measurement get_z() const { return z; }
        Measurement length() const;
        // exception safety: strong guarantee - throws length_error if length() is zero
        Vector_3D& normalize();
        Vector_3D& operator+=(const Vector_3D&);
        Vector_3D& operator-=(const Vector_3D&);
        Vector_3D operator-() { return Vector_3D(-x, -y, -z); }
        Vector_3D& operator*=(const Measurement);
    private:
        Measurement x;
        Measurement y;
        Measurement z;
    };
    
    // exception safety: no throw
    const Vector_3D cross_product(const Vector_3D&, const Vector_3D&);
    // exception safety: no throw
    const Vector_3D::Measurement dot_product(const Vector_3D&, const Vector_3D&);
    // exception safety: no throw
    const Vector_3D::Angle_Meas angle_between(const Vector_3D&, const Vector_3D&);
    // exception safety: no throw
    const Vector_3D operator+(const Vector_3D&, const Vector_3D&); // translate
    // exception safety: no throw
    const Vector_3D operator-(const Vector_3D&, const Vector_3D&); // translate
    // exception safety: no throw
    const Vector_3D operator*(const Vector_3D&, const Vector_3D::Measurement); // scale
    // exception safety: no throw
    const Vector_3D operator*(const Vector_3D::Measurement, const Vector_3D&); // scale
    // orthogonal projection of a onto b
    // exception safety: strong guarantee - throws length error if b.length() is zero
    const Vector_3D orthogonal_projection(const Vector_3D& a, const Vector_3D& b);
    
    /* 
     * determines if two vectors are within a rounding error of each other
     */
//    bool within_round(const Vector_3D& v1, const Vector_3D& v2, 
//            const Vector_3D::Measurement precision);
    
    const bool is_equal(const Vector_3D& vector1, const Vector_3D& vector2, 
            const Vector_3D::Measurement precision);
    
    const bool is_cp_zero(const Vector_3D& v1, const Vector_3D& v2, 
            const Vector_3D::Measurement precision);
    
    const bool is_same_line(const Point_3D& v1_start, const Point_3D& v1_end, 
            const Point_3D& v2_start, const Point_3D& v2_end, bool& same_direction,
            const Vector_3D::Measurement precision);
    
    const bool is_pt_on_vector(const Point_3D& pt, const Point_3D& v_start, 
            const Point_3D& v_end, const Vector_3D::Measurement precision);
    
//    // exception safety: no throw
//    const bool is_pt_btwn_vectors(const Point_3D& p, const Point_3D& v_start, 
//            const Point_3D& v1_end, const Point_3D& v2_end, const Vector_3D::Measurement precision);
    
    // struct defined in Point_3D
    struct Vector_3D_idata;
    
    /*
     * Intersect v1 into v2.  Possible for 0, 1, or 2 intersection points. If
     * true is returned, one or two points will be set in idata. 
     * 
     * exception safety: no throw
     */
    const bool intersect_vectors(
            const Point_3D& v1_start, const Point_3D& v1_end, 
            const Point_3D& v2_start, const Point_3D& v2_end, 
            Vector_3D_idata& idata, const Vector_3D::Measurement precision);
    
    // intersect_vectors helper functions
    
    // intersect_vectors same line
    const bool intersect_vectors_sl(const Point_3D& v1_start, 
            const Point_3D& v1_end, const Point_3D& v2_start, 
            const Point_3D& v2_end, const bool same_direction,
            Vector_3D_idata& idata, const Vector_3D::Measurement precision);
    
    enum axis { X, Y, Z };
    
    // intersect_vectors different lines, solve for t1 first
    const bool intersect_vectors_dl_t1(const Point_3D& o1, 
            const Vector_3D& v1, const Point_3D& o2, const Vector_3D& v2,
            const axis axis_a, const axis axis_b, Point_3D& i_point, 
            const Vector_3D::Measurement precision);
    
    // intersect_vectors different lines, solve for t2 first
    const bool intersect_vectors_dl_t2(const Point_3D& o1, 
            const Vector_3D& v1, const Point_3D& o2, const Vector_3D& v2,
            const axis axis_a, const axis axis_b, Point_3D& i_point, 
            const Vector_3D::Measurement precision);
    
    // intersect_vectors different lines, solve for cross product equal to zero
    const bool intersect_vectors_dl_cpz(const Point_3D& v1_start, const Point_3D& v1_end, 
            const Point_3D& v2_start, const Point_3D& v2_end, Point_3D& i_point, 
            const Vector_3D::Measurement precision);
    
    // intersect_vectors different lines
    const bool intersect_vectors_dl(const Point_3D& v1_start, 
            const Point_3D& v1_end, const Point_3D& v2_start, 
            const Point_3D& v2_end, Vector_3D_idata& idata, 
            const Vector_3D::Measurement precision);
}

#endif /* VECTOR_3D_H */

