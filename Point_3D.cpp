/* 
 * File:   Point_3D.cpp
 * Author: Jeffrey Davis
 * Version: 1.0
 * 
 * Created on February 23, 2016, 12:59 AM
 */

#include "Point_3D.h"
#include <iostream>
#include <stdexcept>
#include <sstream>
#include <cmath>
#include <cfloat>

namespace VCAD_lib
{
    Point_3D::Point_3D(const Measurement x_val, const Measurement y_val, const Measurement z_val)
                : x(x_val), y(y_val), z(z_val) {}
    
    const Point_3D::Angle_Meas Point_3D::pi = 3.141592653589793238462643;
    
    Point_3D::Point_3D(const Vector_3D& v) : x(v.get_x()), y(v.get_y()), z(v.get_z()) {}
    
    Point_3D::Point_3D(const Point_3D& origin, const Vector_3D& v) : 
            x(origin.x + v.get_x()),
            y(origin.y + v.get_y()),
            z(origin.z + v.get_z()) {}
    
    Point_3D::Measurement Point_3D::get_r() const
    {
        if (x == 0 && y == 0)
            return 0;
        else
            return sqrt(pow(x,2) + pow(y,2));
    }

    Point_3D::Measurement Point_3D::get_r(const Point_3D& origin) const
    {
        if (x == origin.get_x() && y == origin.get_y())
            return 0;
        else
        {
            Measurement x_val = x - origin.get_x();
            Measurement y_val = y - origin.get_y();
            return sqrt(pow(x_val,2) + pow(y_val,2));
        }
    }
    
    Point_3D::Measurement Point_3D::get_p() const
    {
        if (x == 0 && y == 0 && z == 0)
            return 0;
        else
            return sqrt(pow(x,2) + pow(y,2) + pow(z,2));
    }

    Point_3D::Measurement Point_3D::get_p(const Point_3D& origin) const
    {
        if (x == origin.get_x() && y == origin.get_y() && z == origin.get_z())
            return 0;
        else
        {
            Measurement x_val = x - origin.get_x();
            Measurement y_val = y - origin.get_y();
            Measurement z_val = z - origin.get_z();
            return sqrt(pow(x_val,2) + pow(y_val,2) + pow(z_val,2));
        }
    }

    Point_3D::Angle_Meas Point_3D::get_theta() const
    {
        if (x == 0 && y == 0)
            return 0;
        else
        {
            Angle_Meas theta = acos(x / this->get_r());
            if (y < 0) // third or fourth quadrant
                theta = (pi * 2) - theta;
            return theta;
        }
    }

    Point_3D::Angle_Meas Point_3D::get_theta(const Point_3D& origin) const
    {
        if (x == origin.get_x() && y == origin.get_y())
            return 0;
        else
        {
            Angle_Meas theta = acos((x - origin.get_x()) / this->get_r(origin));
            if (y - origin.get_y() < 0) // third or fourth quadrant
                theta = (pi * 2) - theta;
            return theta;
        }
    }

    Point_3D::Angle_Meas Point_3D::get_phi() const
    {
        if (x == 0 && y == 0)
            return 0;
        else
        {
            Angle_Meas phi = acos(dot_product(Vector_3D(x,y,z), Vector_3D(0,0,1)) / this->get_p());
            if (phi < 0)
                phi = (pi / 2) - phi;
            return phi;
        }
    }

    Point_3D::Angle_Meas Point_3D::get_phi(const Point_3D& origin) const
    {
        if (x == origin.get_x() && y == origin.get_y())
            return 0;
        else
        {
            Angle_Meas phi = acos(dot_product(Vector_3D(x - origin.get_x(), y - origin.get_y(), z - origin.get_z()), Vector_3D(0,0,1)) / this->get_p(origin));
            if (phi < 0)
                phi = (pi / 2) - phi;
            return phi;
        }
    }
    
    void Point_3D::rotate_about_axis(Point_3D& p, const Point_3D::Angle_Meas angle, 
            const Vector_3D& axis, const Point_3D& origin)
    {
        Vector_3D v(origin, p);

        // only process if point is not on the axis
        if (cross_product(v, axis).length() != 0)
        {
            // rotation using right hand rule with thumb pointing in direction of axis
            
            // orthog_v is an orthogonal projection of v along the axis ('z' axis)
            Vector_3D orthog_v = orthogonal_projection(v,axis);
            if (orthog_v.length() == 0)
            {
                // v is the 'x' axis
                Vector_3D other_axis = cross_product(axis,v); // 'y' axis
                // axis is 'z' axis
                
                Measurement radius = v.length();
                
                // 'x' axis
                v.normalize();
                v *= radius * cos(angle);
                
                // 'y' axis
                other_axis.normalize();
                other_axis *= radius * sin(angle);
                
                p = origin + v + other_axis;
            }
            else 
            {
                // radius_v is a perpendicular vector from the line origin + axis to point p ('x' axis)
                // other_axis is the third axis ('y' axis)
                Vector_3D radius_v = v - orthog_v;
                Vector_3D other_axis = cross_product(axis, radius_v);

                Measurement radius = radius_v.length();

                // 'x' vector
                radius_v.normalize();
                radius_v *= radius * cos(angle);

                // 'y' vector
                other_axis.normalize();
                other_axis *= radius * sin(angle);

                p = origin + orthog_v + other_axis + radius_v;
            }
        }
    }
    
    Point_3D& Point_3D::rotate(const Angle& angle)
    {
        return this->rotate(angle, Point_3D(0,0,0));
    }
    
    Point_3D& Point_3D::rotate(const Angle_Meas angle, const Vector_3D& axis)
    {
        return this->rotate(angle, axis, Point_3D(0,0,0));
    }

    Point_3D& Point_3D::rotate(const Angle& angle, const Point_3D& origin)
    {
        if (angle.angle_x != 0)
            this->rotate_about_axis(*this, angle.angle_x, Vector_3D(1,0,0), origin);

        if (angle.angle_y != 0)
            this->rotate_about_axis(*this, angle.angle_y, Vector_3D(0,1,0), origin);

        if (angle.angle_z != 0)
            this->rotate_about_axis(*this, angle.angle_z, Vector_3D(0,0,1), origin);

        return *this;
    }
    
    Point_3D& Point_3D::rotate(const Angle_Meas angle, const Vector_3D& axis, const Point_3D& origin)
    {
        this->rotate_about_axis(*this, angle, axis, origin);
        return *this;
    }

    Point_3D& Point_3D::scale(const Point_3D::Measurement x_scalar, const Point_3D::Measurement y_scalar, const Point_3D::Measurement z_scalar)
    {
        x *= x_scalar;
        y *= y_scalar;
        z *= z_scalar;
        return *this;
    }

    Point_3D& Point_3D::scale(const Point_3D::Measurement x_scalar, const Point_3D::Measurement y_scalar, const Point_3D::Measurement z_scalar,
            const Point_3D& origin)
    {
        x = origin.x + (x - origin.x) * x_scalar;
        y = origin.y + (y - origin.y) * y_scalar;
        z = origin.z + (z - origin.z) * z_scalar;
        return *this;
    }

    Point_3D& Point_3D::translate(const Point_3D::Measurement x_val, const Point_3D::Measurement y_val, const Point_3D::Measurement z_val)
    {
        x += x_val;
        y += y_val;
        z += z_val;
        return *this;
    }

    Point_3D& Point_3D::translate(const Vector_3D& vector)
    {
        x += vector.get_x();
        y += vector.get_y();
        z += vector.get_z();
        return *this;
    }
    
    Point_3D& Point_3D::move_x_pxy(const Point_3D& new_origin, const Vector_3D& x_axis, 
            const Point_3D& pt_xy_plane, const Point_3D& ref_origin)
    {
        // determine y_axis by using pt_xy_plane point
        Vector_3D y_axis(new_origin, pt_xy_plane);
        y_axis -= orthogonal_projection(y_axis, x_axis);
        
        move(new_origin, x_axis, y_axis, cross_product(x_axis, y_axis), ref_origin);
        return *this;
    }
    
    Point_3D& Point_3D::move_x_pxz(const Point_3D& new_origin, const Vector_3D& x_axis, 
            const Point_3D& pt_xz_plane, const Point_3D& ref_origin)
    {
        // determine z_axis by using pt_xz_plane point
        Vector_3D z_axis(new_origin, pt_xz_plane);
        z_axis -= orthogonal_projection(z_axis, x_axis);
        
        move(new_origin, x_axis, cross_product(z_axis, x_axis), z_axis, ref_origin);
        return *this;
    }
    
    Point_3D& Point_3D::move_y_pxy(const Point_3D& new_origin, const Vector_3D& y_axis, 
            const Point_3D& pt_xy_plane, const Point_3D& ref_origin)
    {
        // determine x_axis by using pt_xy_plane point
        Vector_3D x_axis(new_origin, pt_xy_plane);
        x_axis -= orthogonal_projection(x_axis, y_axis);
        
        move(new_origin, x_axis, y_axis, cross_product(x_axis, y_axis), ref_origin);
        return *this;
    }
    
    Point_3D& Point_3D::move_y_pyz(const Point_3D& new_origin, const Vector_3D& y_axis, 
            const Point_3D& pt_yz_plane, const Point_3D& ref_origin)
    {
        // determine z_axis by using pt_yz_plane point
        Vector_3D z_axis(new_origin, pt_yz_plane);
        z_axis -= orthogonal_projection(z_axis, y_axis);
        
        move(new_origin, cross_product(y_axis, z_axis), y_axis, z_axis, ref_origin);
        return *this;
    }
    
    Point_3D& Point_3D::move_z_pxz(const Point_3D& new_origin, const Vector_3D& z_axis, 
            const Point_3D& pt_xz_plane, const Point_3D& ref_origin)
    {
        // determine x_axis by using pt_xz_plane point
        Vector_3D x_axis(new_origin, pt_xz_plane);
        x_axis -= orthogonal_projection(x_axis, z_axis);
        
        move(new_origin, x_axis, cross_product(z_axis, x_axis), z_axis, ref_origin);
        return *this;
    }
    
    Point_3D& Point_3D::move_z_pyz(const Point_3D& new_origin, const Vector_3D& z_axis, 
            const Point_3D& pt_yz_plane, const Point_3D& ref_origin)
    {
        // determine y_axis by using pt_yz_plane point
        Vector_3D y_axis(new_origin, pt_yz_plane);
        y_axis -= orthogonal_projection(y_axis, z_axis);
        
        move(new_origin, cross_product(y_axis, z_axis), y_axis, z_axis, ref_origin);
        return *this;
    }
    
    void Point_3D::move(const Point_3D& new_origin, const Vector_3D& new_x_axis, 
            const Vector_3D& new_y_axis, const Vector_3D& new_z_axis, 
            const Point_3D& ref_origin)
    {
        Vector_3D v(ref_origin, *this);

        Vector_3D x_axis(new_x_axis);
        x_axis.normalize();
        x_axis *= v.get_x();
        
        Vector_3D y_axis(new_y_axis);
        y_axis.normalize();
        y_axis *= v.get_y();
        
        Vector_3D z_axis(new_z_axis);
        z_axis.normalize();
        z_axis *= v.get_z();
        
        Point_3D p = new_origin + x_axis + y_axis + z_axis;
        x = p.get_x();
        y = p.get_y();
        z = p.get_z();
    }
    
    Point_3D& Point_3D::operator+=(const Vector_3D& v)
    {
        x += v.get_x();
        y += v.get_y();
        z += v.get_z();
        return *this;
    }
    
    Point_3D& Point_3D::operator-=(const Vector_3D& v)
    {
        x -= v.get_x();
        y -= v.get_y();
        z -= v.get_z();
        return *this;
    }
    
    Point_3D& Point_3D::operator*=(const Point_3D::Measurement val)
    {
        x *= val;
        y *= val;
        z *= val;
        return *this;
    }
    
    const Point_3D operator+(const Point_3D& p, const Vector_3D& v)
    {
        return Point_3D(p.get_x() + v.get_x(), p.get_y() + v.get_y(), p.get_z() + v.get_z());
    }
    
    const Point_3D operator-(const Point_3D& p, const Vector_3D& v)
    {
        return Point_3D(p.get_x() - v.get_x(), p.get_y() - v.get_y(), p.get_z() - v.get_z());
    }
    
    const Point_3D operator*(const Point_3D& p, const Point_3D::Measurement val)
    {
        return Point_3D(p.get_x() * val, p.get_y() * val, p.get_z() * val);
    }
    
    Angle::Angle(const Point_3D::Angle_Meas& x, const Point_3D::Angle_Meas& y, const Point_3D::Angle_Meas& z) : 
            angle_x(x), angle_y(y), angle_z(z) {}
    
    const bool is_equal(const Point_3D& point1, const Point_3D& point2, 
            const Point_3D::Measurement precision)
    {
#ifdef DEBUG_POINT_3D_IS_EQUAL
        cout << "Point_3D is_equal begin\n";
        cout << "Point_3D is_equal point1 x: " << point1.get_x() << " y: " << point1.get_y() << " z: " << point1.get_z() << "\n";
        cout << "Point_3D is_equal point2 x: " << point2.get_x() << " y: " << point2.get_y() << " z: " << point2.get_z() << "\n";
#endif
        Point_3D::Measurement largest(fabs(point1.get_x()));
        Point_3D::Measurement error_bound(fabs(point2.get_x()));
        if (error_bound > largest)
            largest = error_bound;

        error_bound = largest > 1.0 ? (precision * largest) : precision;

        if (fabs(point1.get_x() - point2.get_x()) <= error_bound)
        {
            largest = fabs(point1.get_y());
            error_bound = fabs(point2.get_y());
            if (error_bound > largest)
                largest = error_bound;

            error_bound = largest > 1.0 ? (precision * largest) : precision;

            if (fabs(point1.get_y() - point2.get_y()) <= error_bound)
            {
                largest = fabs(point1.get_z());
                error_bound = fabs(point2.get_z());
                if (error_bound > largest)
                    largest = error_bound;
                
                error_bound = largest > 1.0 ? (precision * largest) : precision;
#ifdef DEBUG_POINT_3D_IS_EQUAL
                if (fabs(point1.get_z() - point2.get_z()) > error_bound)
                {
                    // determine what multiple of DBL_EPSILON is needed to succeed
                    Vector_3D::Measurement val(fabs(point1.get_z() - point2.get_z()));
                    Vector_3D::Measurement divisor(DBL_EPSILON);
                    if (largest > 1.0)
                        divisor *= largest;
                    cout << "Point_3D is_equal z (|p1z: " << point1.get_z() << " - p2z: " << point2.get_z() << "| = " << val << ") FALSE needs a precision of at least: " << (val / divisor) << " * DBL_EPSILON\n";
                }
#endif                
                return fabs(point1.get_z() - point2.get_z()) <= error_bound;
            }
#ifdef DEBUG_POINT_3D_IS_EQUAL
            else
            {
                // determine what multiple of DBL_EPSILON is needed to succeed
                Vector_3D::Measurement val(fabs(point1.get_y() - point2.get_y()));
                Vector_3D::Measurement divisor(DBL_EPSILON);
                if (largest > 1.0)
                    divisor *= largest;
                cout << "Point_3D is_equal y (|p1y: " << point1.get_y() << " - p2y: " << point2.get_y() << "| = " << val << ") FALSE needs a precision of at least: " << (val / divisor) << " * DBL_EPSILON\n";
            }
#endif
        }
#ifdef DEBUG_POINT_3D_IS_EQUAL
        else
        {
            // determine what multiple of DBL_EPSILON is needed to succeed
            Vector_3D::Measurement val(fabs(point1.get_x() - point2.get_x()));
            Vector_3D::Measurement divisor(DBL_EPSILON);
            if (largest > 1.0)
                divisor *= largest;
            cout << "Point_3D is_equal x (|p1x: " << point1.get_x() << " - p2x: " << point2.get_x() << "| = " << val << ") FALSE needs a precision of at least: " << (val / divisor) << " * DBL_EPSILON\n";
        }
#endif
        
        return false;
    }
    
//
//    const Point_3D c_coord_sys(const Point_3D& point, const Point_3D& new_origin, 
//            const Vector_3D& new_x_axis, const Vector_3D& new_y_axis,
//            const Vector_3D& new_z_axis)
//    {
//        // vector from new origin to point
//        Vector_3D v(new_origin, point);
//        Vector_3D v_x = orthogonal_projection(v, new_x_axis);
//        Vector_3D v_y = orthogonal_projection(v, new_y_axis);
//        Vector_3D v_z = orthogonal_projection(v, new_z_axis);
//        Point_3D::Measurement x = v_x.length();
//        Point_3D::Measurement y = v_y.length();
//        Point_3D::Measurement z = v_z.length();
//        if (dot_product(v_x, new_x_axis) < 0)
//            x = -x;
//        if (dot_product(v_y, new_y_axis) < 0)
//            y = -y;
//        if (dot_product(v_z, new_z_axis) < 0)
//            z = -z;
//        return Point_3D(x, y, z); 
//    }
//    
//    const Point_3D move(const Point_3D& point, const Point_3D& ref, 
//            const Point_3D& origin, const Vector_3D& x_axis, 
//            const Vector_3D& y_axis, const Vector_3D& z_axis)
//    {
//        // alter based off of ref point
//        Vector_3D v(ref,point);
//        // go back to original coordinate system
//        // convert point back to each alternate axis vectors
//        Vector_3D alt_x(x_axis);
//        alt_x.normalize();
//        alt_x *= v.get_x();
//        Vector_3D alt_y(y_axis);
//        alt_y.normalize();
//        alt_y *= v.get_y();
//        Vector_3D alt_z(z_axis);
//        alt_z.normalize();
//        alt_z *= v.get_z();
//        // get new point by taking origin and adding alt axis vectors
//        Point_3D p(origin, alt_x);
//        p += alt_y;
//        p += alt_z;
//        return p;
//    }
    
    const Point_3D cylindrical_point(const Point_3D::Measurement r, const Point_3D::Angle_Meas theta, const Point_3D::Measurement z)
    {
        if (r < 0)
        {
            stringstream ss;
            ss << "r must be greater than or equal to zero. r: " << r;
            throw domain_error(ss.str());
        }
        return Point_3D(r * cos(theta), r * sin(theta), z);
    }
    
    const Point_3D cylindrical_point(const Point_3D::Measurement r, const Point_3D::Angle_Meas theta, const Point_3D::Measurement z,
            const Point_3D& origin)
    {
        if (r < 0)
        {
            stringstream ss;
            ss << "r must be greater than or equal to zero. r: " << r;
            throw domain_error(ss.str());
        }
        return Point_3D(origin.get_x() + (r * cos(theta)), 
                origin.get_y() + (r * sin(theta)), origin.get_z() + z);
    }
    
    const Point_3D spherical_point(Point_3D::Measurement p, Point_3D::Angle_Meas theta, Point_3D::Angle_Meas phi)
    {
        if (p < 0)
        {
            stringstream ss;
            ss << "p must be greater or equal to zero. p: " << p;
            throw domain_error(ss.str());
        }
        const Point_3D::Angle_Meas pi = 3.1415926535897932384;
        const Point_3D::Angle_Meas two_pi = pi * 2;
        // 0 <= theta < 2 pi
        // 0 <= phi <= pi
        while (phi > two_pi)
            phi -= two_pi;
        while (phi < 0)
            phi += two_pi;
        if (phi > pi)
            phi -= pi;
        
        return Point_3D(p * sin(phi) * cos(theta), p * sin(phi) * sin(theta), p * cos(phi));
    }
    
    const Point_3D spherical_point(Point_3D::Measurement p, Point_3D::Angle_Meas theta, Point_3D::Angle_Meas phi,
            const Point_3D& origin)
    {
        if (p < 0)
        {
            stringstream ss;
            ss << "p must be greater than or equal to zero. p: " << p;
            throw domain_error(ss.str());
        }
        const Point_3D::Angle_Meas pi = 3.1415926535897932384;
        const Point_3D::Angle_Meas two_pi = pi * 2;
        // 0 <= theta < 2 pi
        // 0 <= phi <= pi
        while (phi > two_pi)
            phi -= two_pi;
        while (phi < 0)
            phi += two_pi;
        if (phi > pi)
            phi -= pi;
        
        return Point_3D(origin.get_x() + (p * sin(phi) * cos(theta)), 
                origin.get_y() + (p * sin(phi) * sin(theta)), 
                origin.get_z() + (p * cos(phi)));
    }
    
    Vector_3D_idata::Vector_3D_idata() : num(0), p1(0,0,0), p2(0,0,0) {}
}
