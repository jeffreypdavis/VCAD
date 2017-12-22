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
 * File:   Point_2D.cpp
 * Author: Jeffrey Davis
 */

#include "Point_2D.h"
#include <sstream>
#include <stdexcept>

namespace VCAD_lib
{
    Point_2D::Point_2D(const Measurement x_val, const Measurement y_val) : 
            x(x_val), y(y_val) {}
    
    const Point_2D::Angle_Meas Point_2D::two_pi = 6.283185307179586;

    Point_2D::Measurement Point_2D::get_r() const
    {
        if (x == 0 && y == 0)
            return 0;
        else
            return sqrt(pow(x,2) + pow(y,2));
    }

    Point_2D::Measurement Point_2D::get_r(const Point_2D& origin) const
    {
        if (x == origin.get_x() && y == origin.get_y())
            return 0;
        else
        {
            Measurement xVal = x - origin.get_x();
            Measurement yVal = y - origin.get_y();
            return sqrt(pow(xVal,2) + pow(yVal,2));
        }
    }

    Point_2D::Angle_Meas Point_2D::get_theta() const
    {
        if (x == 0 && y == 0)
            return 0;
        else
        {
            Angle_Meas theta = acos(x / this->get_r());
            if (y < 0) // third or fourth quadrant
                theta = two_pi - theta;
            return theta;
        }
    }

    Point_2D::Angle_Meas Point_2D::get_theta(const Point_2D& origin) const
    {
        if (x == origin.get_x() && y == origin.get_y())
            return 0;
        else
        {
            Angle_Meas theta = acos((x - origin.get_x()) / this->get_r(origin));
            if (y - origin.get_y() < 0) // third or fourth quadrant
                theta = two_pi - theta;
            return theta;
        }
    }

    Point_2D& Point_2D::rotate(const Angle_Meas angle)
    {
        if (!(x == 0 && y == 0))
        {
            Angle_Meas theta = this->get_theta() + angle;
            Measurement r = this->get_r();
            x = r * cos(theta);
            y = r * sin(theta);
        }
        return *this;
    }

    Point_2D& Point_2D::rotate(const Angle_Meas angle, const Point_2D& origin)
    {
        if (!(x == origin.get_x() && y == origin.get_y()))
        {
            Angle_Meas theta = this->get_theta(origin) + angle;
            Measurement r = this->get_r(origin);
            x = r * cos(theta) + origin.get_x();
            y = r * sin(theta) + origin.get_y();
        }
        return *this;
    }

    Point_2D& Point_2D::scale(const Measurement x_scalar, const Measurement y_scalar)
    {
        x *= x_scalar;
        y *= y_scalar;
        return *this;
    }

    Point_2D& Point_2D::scale(const Measurement x_scalar, const Measurement y_scalar,
            const Point_2D& origin)
    {
        x = origin.x + (x - origin.x) * x_scalar;
        y = origin.y + (y - origin.y) * y_scalar;
        return *this;
    }

    Point_2D& Point_2D::translate(const Measurement x_val, const Measurement y_val)
    {
        x += x_val;
        y += y_val;
        return *this;
    }

    Point_2D& Point_2D::translate(const Vector_2D& vector)
    {
        return this->translate(vector.get_x(), vector.get_y());
    }
    
    Point_2D& Point_2D::move(const Point_2D& new_origin, const Vector_2D& axis,
            const bool is_x_axis, const Point_2D& ref_origin)
    {
        
        Vector_2D v(ref_origin, *this);

        if (is_x_axis)
        {
            // determine y-axis
            Point_2D y_axis_pt(axis.get_x(), axis.get_y());
            y_axis_pt.rotate(two_pi / 4);
            
            Vector_2D y_axis(y_axis_pt);
//            y_axis -= orthogonal_projection(y_axis, axis);
//            if (cross_product(axis, y_axis) < 0)
//                y_axis = -y_axis;
            
            Vector_2D x_axis(axis);
            x_axis.normalize();
            x_axis *= v.get_x();

            y_axis.normalize();
            y_axis *= v.get_y();

            Point_2D p = new_origin + x_axis + y_axis;
            x = p.get_x();
            y = p.get_y();
        }
        else
        {
            // determine x-axis
            Point_2D x_axis_pt(axis.get_x(), axis.get_y());
            x_axis_pt.rotate(-two_pi / 4);
            
            Vector_2D x_axis(x_axis_pt);
//            x_axis -= orthogonal_projection(x_axis, axis);
//            if (cross_product(x_axis, axis) < 0)
//                x_axis = -x_axis;

            x_axis.normalize();
            x_axis *= v.get_x();
            
            Vector_2D y_axis(axis);
            y_axis.normalize();
            y_axis *= v.get_y();

            Point_2D p = new_origin + x_axis + y_axis;
            x = p.get_x();
            y = p.get_y();
        }
        
        return *this;
    }
    
    Point_2D& Point_2D::operator+=(const Vector_2D& v)
    {
        x += v.get_x();
        y += v.get_y();
        return *this;
    }
    
    Point_2D& Point_2D::operator-=(const Vector_2D& v)
    {
        x -= v.get_x();
        y -= v.get_y();
        return *this;
    }
    
    Point_2D& Point_2D::operator*=(const double val)
    {
        x *= val;
        y *= val;
        return *this;
    }
    
    const Point_2D operator+(const Point_2D& p, const Vector_2D& v)
    {
        return Point_2D(p.get_x() + v.get_x(), p.get_y() + v.get_y());
    }
    
    const Point_2D operator-(const Point_2D& p, const Vector_2D& v)
    {
        return Point_2D(p.get_x() - v.get_x(), p.get_y() - v.get_y());
    }
    
    const Point_2D operator*(const Point_2D& p, const Point_2D::Measurement val)
    {
        return Point_2D(p.get_x() * val, p.get_y() * val);
    }
    
    const bool is_equal(const Point_2D& point1, const Point_2D& point2, 
            const Point_2D::Measurement precision)
    {
#ifdef DEBUG_POINT_2D_IS_EQUAL
        cout << "Point_2D is_equal begin\n";
        cout << "Point_2D is_equal point1 x: " << point1.get_x() << " y: " << point1.get_y() << "\n";
        cout << "Point_2D is_equal point2 x: " << point2.get_x() << " y: " << point2.get_y() << "\n";
#endif
        Point_2D::Measurement largest(fabs(point1.get_x()));
        Point_2D::Measurement error_bound(fabs(point2.get_x()));
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

#ifdef DEBUG_POINT_2D_IS_EQUAL
            if (fabs(point1.get_y() - point2.get_y()) > error_bound)
            {
                // determine what multiple of DBL_EPSILON is needed to succeed
                Vector_2D::Measurement val(fabs(point1.get_y() - point2.get_y()));
                Vector_2D::Measurement divisor(DBL_EPSILON);
                if (largest > 1.0)
                    divisor *= largest;
                cout << "Point_2D is_equal y (|p1y: " << point1.get_y() << " - p2y: " << point2.get_y() << "| = " << val << ") FALSE needs a precision of at least: " << (val / divisor) << " * DBL_EPSILON\n";
            }
#endif                
            return fabs(point1.get_y() - point2.get_y()) <= error_bound;
        }
#ifdef DEBUG_POINT_2D_IS_EQUAL
        else
        {
            // determine what multiple of DBL_EPSILON is needed to succeed
            Vector_2D::Measurement val(fabs(point1.get_x() - point2.get_x()));
            Vector_2D::Measurement divisor(DBL_EPSILON);
            if (largest > 1.0)
                divisor *= largest;
            cout << "Point_2D is_equal x (|p1x: " << point1.get_x() << " - p2x: " << point2.get_x() << "| = " << val << ") FALSE needs a precision of at least: " << (val / divisor) << " * DBL_EPSILON\n";
        }
#endif
        
        return false;
    }
    
    const Point_2D polar_point(const Point_2D::Measurement r, const Point_2D::Angle_Meas theta)
    {
        if (r < 0)
        {
            stringstream ss;
            ss << "r must be greater than or equal to zero. r: " << r;
            throw domain_error(ss.str());
        }
        return Point_2D(r * cos(theta), r * sin(theta));
    }
    
    const Point_2D polar_point(const Point_2D::Measurement r, const Point_2D::Angle_Meas theta, const Point_2D& origin)
    {
        if (r < 0)
        {
            stringstream ss;
            ss << "r must be greater than or equal to zero. r: " << r;
            throw domain_error(ss.str());
        }
        return Point_2D(r * cos(theta) + origin.get_x(), r * sin(theta) + origin.get_y());
    }

    Vector_2D_idata::Vector_2D_idata() : num(0), p1(0,0), p2(0,0) {}
}