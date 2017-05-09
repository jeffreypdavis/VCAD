/* 
 * File:   Point_3D.h
 * Author: Jeffrey Davis
 * Version: 1.0
 *
 * Created on February 23, 2016, 12:59 AM
 */

#ifndef POINT_3D_H
#define POINT_3D_H

#include "Vector_3D.h"
#include <cmath>

using namespace std;

namespace VCAD_lib
{
    class Vector_3D;
    
    // holds rotation angles for x, y, and z axis
    struct Angle;
    
    // exception safety: no throw guarantee
    class Point_3D {
    public:
        typedef double Measurement;
        typedef double Angle_Meas;
        Point_3D(const Measurement x_val, const Measurement y_val, 
                const Measurement z_val);
        Point_3D(const Vector_3D& v); // 0,0,0 + vector v
        // origin + vector v
        Point_3D(const Point_3D& origin, const Vector_3D& v); 
        Measurement get_x() const { return x; }
        Measurement get_y() const { return y; }
        Measurement get_z() const { return z; }
        Measurement get_r() const; // x-y radius for cylindrical points
        // x-y radius for cylindrical points
        Measurement get_r(const Point_3D& origin) const; 
        Measurement get_p() const; // x-y-z radius for spherical points
        // x-y-z radius for spherical points
        Measurement get_p(const Point_3D& origin) const; 
        Measurement get_theta() const;
        Measurement get_theta(const Point_3D& origin) const;
        Measurement get_phi() const; // for spherical points
        // for spherical points
        Measurement get_phi(const Point_3D& origin) const; 
        Point_3D& rotate(const Angle& angle);
        Point_3D& rotate(const Angle_Meas angle, const Vector_3D& axis);
        Point_3D& rotate(const Angle& angle, const Point_3D& origin);
        Point_3D& rotate(const Angle_Meas angle, const Vector_3D& axis, 
                const Point_3D& origin);
        Point_3D& scale(const Measurement x_scalar, 
                const Measurement y_scalar, const Measurement z_scalar);
        Point_3D& scale(const Measurement x_scalar, 
                const Measurement y_scalar, const Measurement z_scalar, 
                const Point_3D& origin);
        Point_3D& translate(const Measurement x_val, const Measurement y_val, 
                const Measurement z_val);
        Point_3D& translate(const Vector_3D&);
        /*
         * move the point to a different coordinate system and origin.
         * A new coordinate system must be supplied and this is done by 
         * specifying a new origin, a new axis, and one more point that is
         * not on the axis to define a plane.  ref_origin is an option to 
         * move based on a different origin than 0,0,0.
         * 
         * There are six methods to allow for specifying the new coordinate 
         * system. The names of the methods start with move.  Then next 
         * letter is the axis specified in the vector argument - x, y, or z 
         * axis.  The last part of the name is identifying where the third 
         * point is located - xy plane, xz plane, or yz plane.
         */
        Point_3D& move_x_pxy(const Point_3D& new_origin, 
                const Vector_3D& x_axis, const Point_3D& pt_xy_plane, 
                const Point_3D& ref_origin=Point_3D(0,0,0));
        Point_3D& move_x_pxz(const Point_3D& new_origin, 
                const Vector_3D& x_axis, const Point_3D& pt_xz_plane, 
                const Point_3D& ref_origin=Point_3D(0,0,0));
        Point_3D& move_y_pxy(const Point_3D& new_origin, 
                const Vector_3D& y_axis, const Point_3D& pt_xy_plane, 
                const Point_3D& ref_origin=Point_3D(0,0,0));
        Point_3D& move_y_pyz(const Point_3D& new_origin, 
                const Vector_3D& y_axis, const Point_3D& pt_yz_plane, 
                const Point_3D& ref_origin=Point_3D(0,0,0));
        Point_3D& move_z_pxz(const Point_3D& new_origin, 
                const Vector_3D& z_axis, const Point_3D& pt_xz_plane, 
                const Point_3D& ref_origin=Point_3D(0,0,0));
        Point_3D& move_z_pyz(const Point_3D& new_origin, 
                const Vector_3D& z_axis, const Point_3D& pt_yz_plane, 
                const Point_3D& ref_origin=Point_3D(0,0,0));
        Point_3D& operator+=(const Vector_3D&);
        Point_3D& operator-=(const Vector_3D&);
        Point_3D& operator*=(const Measurement);
    private:
        static const Angle_Meas pi;
        Measurement x;
        Measurement y;
        Measurement z;
        // rotate point about an arbitrary axis using right hand rule with 
        // thumb in direction of axis
        void rotate_about_axis(Point_3D& p, const Angle_Meas angle, 
                const Vector_3D& axis, const Point_3D& origin);
        
        void move(const Point_3D& new_origin, const Vector_3D& new_x_axis, 
                const Vector_3D& new_y_axis, const Vector_3D& new_z_axis, 
                const Point_3D& ref_origin=Point_3D(0,0,0));
    };
    
    // exception safety: no throw
    const Point_3D operator+(const Point_3D&, const Vector_3D&);
    // exception safety: no throw
    const Point_3D operator-(const Point_3D&, const Vector_3D&);
    // exception safety: no throw
    const Point_3D operator*(const Point_3D&, const Point_3D::Measurement);

    // holds rotation angles for x, y, and z axis
    struct Angle
    {
        Point_3D::Angle_Meas angle_x;
        Point_3D::Angle_Meas angle_y;
        Point_3D::Angle_Meas angle_z;
        Angle(const Point_3D::Angle_Meas& x, const Point_3D::Angle_Meas& y, 
                const Point_3D::Angle_Meas& z);
    };

    /* 
     * determines if two points are within a rounding error of each other
     */
    const bool is_equal(const Point_3D& point1, const Point_3D& point2, 
            const Point_3D::Measurement precision);
    
    // exception safety: strong guarantee - throws domain_error if r is 
    // negative
    const Point_3D cylindrical_point(const Point_3D::Measurement r, 
            const Point_3D::Angle_Meas theta, const Point_3D::Measurement z);
    // exception safety: strong guarantee - throws domain_error if r is 
    // negative
    const Point_3D cylindrical_point(const Point_3D::Measurement r, 
            const Point_3D::Angle_Meas theta, const Point_3D::Measurement z, 
            const Point_3D& origin);
    // exception safety: strong guarantee - throws domain_error if p is 
    // negative
    const Point_3D spherical_point(Point_3D::Measurement p, 
            Point_3D::Angle_Meas theta, Point_3D::Angle_Meas phi);
    // exception safety: strong guarantee - throws domain_error if p is 
    // negative
    const Point_3D spherical_point(Point_3D::Measurement p, 
            Point_3D::Angle_Meas theta, Point_3D::Angle_Meas phi,
            const Point_3D& origin);

    struct Vector_3D_idata
    {
        int num; // number of intersection points found (0,1, or 2))
        Point_3D p1;
        Point_3D p2;
        Vector_3D_idata();
    };
}
#endif /* POINT_3D_H */

