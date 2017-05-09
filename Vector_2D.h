/* 
 * File:   Vector_2D.h
 * Author: Jeffrey Davis
 * Version: 1.0
 *
 * Created on February 11, 2016, 5:42 PM
 */

#ifndef VECTOR_2D_H
#define VECTOR_2D_H
#include "Point_2D.h"
#include <cmath>
#include <iostream>

using namespace std;

namespace VCAD_lib {

    class Point_2D;
    
    /*
     * A vector indicates a length and a direction, but does not directly 
     * correspond to a point.
     * exception safety: no throw guarantee unless otherwise indicated
     */
    class Vector_2D {
    public:
        typedef double Measurement;
        typedef double Angle_Meas;
        Vector_2D(const Measurement x_val, const Measurement y_val);
        Vector_2D(const Point_2D& end_point);
        Vector_2D(const Point_2D& start_point, const Point_2D& end_point);
        Measurement get_x() const { return x; }
        Measurement get_y() const { return y; }
        Measurement length() const;
        // exception safety: strong gurantee - throws length_error if length is zero.
        Vector_2D& normalize();
        Vector_2D& operator+=(const Vector_2D&);
        Vector_2D& operator-=(const Vector_2D&);
        Vector_2D operator-() { return Vector_2D(-x, -y); }
        Vector_2D& operator*=(const Measurement);
    private:
        Measurement x;
        Measurement y;
    };
    
    // exception safety: no throw
    const Vector_2D::Measurement cross_product(const Vector_2D&, const Vector_2D&);
    // exception safety: no throw
    const Vector_2D::Measurement dot_product(const Vector_2D&, const Vector_2D&);
    // exception safety: no throw
    const Vector_2D::Angle_Meas angle_between(const Vector_2D&, const Vector_2D&);
    // orthogonal projection of a onto b
    // exception safety: strong guarantee - throws length_error if b.length() is zero
    const Vector_2D orthogonal_projection(const Vector_2D& a, const Vector_2D& b);
    // exception safety: no throw
    const bool operator==(const Vector_2D&, const Vector_2D&);
    // exception safety: no throw
    const bool operator!=(const Vector_2D&, const Vector_2D&);
    // exception safety: no throw
    const Vector_2D operator+(const Vector_2D&, const Vector_2D&); // translate
    // exception safety: no throw
    const Vector_2D operator-(const Vector_2D&, const Vector_2D&); // translate
    // exception safety: no throw
    const Vector_2D operator*(const Vector_2D&, const Vector_2D::Measurement); // scale

    const bool is_equal(const Vector_2D& vector1, const Vector_2D& vector2, 
            const Vector_2D::Measurement precision);
    
    const bool is_cp_zero(const Vector_2D& v1, const Vector_2D& v2, 
            const Vector_2D::Measurement precision);
    
    const bool is_same_line(const Point_2D& v1_start, const Point_2D& v1_end, 
            const Point_2D& v2_start, const Point_2D& v2_end, bool& same_direction,
            const Vector_2D::Measurement precision);
    
    const bool is_pt_on_vector(const Point_2D& pt, const Point_2D& v_start, 
            const Point_2D& v_end, const Vector_2D::Measurement precision);
    
//    // exception safety: no throw
//    const bool point_on_vector(const Point_2D& p, const Vector_2D& v, const Point_2D& o, 
//            const Vector_2D::Measurement precision);
//    
//    // exception safety: no throw
//    const bool point_btwn_vectors(const Point_2D& p, const Vector_2D& v1, const Vector_2D& v2, 
//            const Point_2D& o, const Vector_2D::Measurement precision);
    
    struct Vector_2D_idata;
    
    const bool intersect_vectors_sl(const Point_2D& v1_start, const Point_2D& v1_end, 
            const Point_2D& v2_start, const Point_2D& v2_end, const bool same_direction,
            Vector_2D_idata& result, const Vector_2D::Measurement precision);
    
    const bool intersect_vectors_dl_cpz(const Point_2D& v1_start, const Point_2D& v1_end, 
            const Point_2D& v2_start, const Point_2D& v2_end, Point_2D& i_point, 
            const Vector_2D::Measurement precision);
    
    const bool intersect_vectors_dl(const Point_2D& v1_start, 
            const Point_2D& v1_end, const Point_2D& v2_start, 
            const Point_2D& v2_end, Vector_2D_idata& result, 
            const Vector_2D::Measurement precision);
    
    /*
     * Intersect v1 into v2 and return 0, 1, or 2 intersection points. if the 
     * result is true, then 1 or 2 points will be stored in idata.  If false,
     * no intersection points were found.
     * 
     * exception safety: no throw
     */
    const bool intersect_vectors(const Point_2D& v1_start, const Point_2D& v1_end, 
            const Point_2D& v2_start, const Point_2D& v2_end, Vector_2D_idata& idata, 
            const Vector_2D::Measurement precision);
}

#endif /* VECTOR_2D_H */

