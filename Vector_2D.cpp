/* 
 * File:   Vector_2D.cpp
 * Author: Jeffrey Davis
 * Version: 1.0
 * 
 * Created on February 11, 2016, 5:42 PM
 */

#include "Vector_2D.h"
#include <stdexcept>

namespace VCAD_lib
{
    Vector_2D::Vector_2D(const Measurement x_val, const Measurement y_val) : x(x_val), y(y_val) {}
    
    Vector_2D::Vector_2D(const Point_2D& end_point) 
    {
        x = end_point.get_x();
        y = end_point.get_y();
    }

    Vector_2D::Vector_2D(const Point_2D& start_point, const Point_2D& end_point) 
    {
        x = end_point.get_x() - start_point.get_x();
        y = end_point.get_y() - start_point.get_y();
    }

    Vector_2D::Measurement Vector_2D::length() const
    {
        return sqrt(pow(x, 2) + pow(y, 2));
    }

    Vector_2D& Vector_2D::normalize()
    {
        Measurement length = this->length();
        if (length == 0)
            throw length_error("Cannot normalize a zero length vector");
        x /= length;
        y /= length;
        return *this;
    }
    
    Vector_2D& Vector_2D::operator+=(const Vector_2D& v)
    {
        x += v.get_x();
        y += v.get_y();
        return *this;
    }
    
    Vector_2D& Vector_2D::operator-=(const Vector_2D& v)
    {
        x -= v.get_x();
        y -= v.get_y();
        return *this;
    }
    
    Vector_2D& Vector_2D::operator*=(const Vector_2D::Measurement val)
    {
        x *= val;
        y *= val;
        return *this;
    }

    const Vector_2D::Measurement cross_product(const Vector_2D& v1, const Vector_2D& v2)
    {
        return v1.get_x() * v2.get_y() - v1.get_y() * v2.get_x();
    }

    const Vector_2D::Measurement dot_product(const Vector_2D& v1, const Vector_2D& v2)
    {
        return (v1.get_x() * v2.get_x()) + (v1.get_y() * v2.get_y());
    }

    const Vector_2D::Angle_Meas angle_between(const Vector_2D& v1, const Vector_2D& v2)
    {
        return acos(dot_product(v1, v2) / (v1.length() * v2.length()));
    }
    
    const Vector_2D orthogonal_projection(const Vector_2D& a, const Vector_2D& b)
    {
        Vector_2D::Measurement length = b.length();
        if (length == 0)
            throw length_error("Cannot orthogonal project a zero length vector");
        return b * (dot_product(a, b) / pow(length ,2));
    }
    
    const bool operator==(const Vector_2D& v1, const Vector_2D& v2)
    {
        return (v1.get_x() == v2.get_x() && v1.get_y() == v2.get_y());
    }
    
    const bool operator!=(const Vector_2D& v1, const Vector_2D& v2)
    {
        return !(v1 == v2);
    }
    
    const Vector_2D operator+(const Vector_2D& v1, const Vector_2D& v2)
    {
        return Vector_2D(v1.get_x() + v2.get_x(), v1.get_y() + v2.get_y());
    }
    
    const Vector_2D operator-(const Vector_2D& v1, const Vector_2D& v2)
    {
        return Vector_2D(v1.get_x() - v2.get_x(), v1.get_y() - v2.get_y());
    }
    
    const Vector_2D operator*(const Vector_2D& v, const Vector_2D::Measurement scalar)
    {
        return Vector_2D(v.get_x() * scalar, v.get_y() * scalar);
    }
    
    
    const bool is_equal(const Vector_2D& vector1, const Vector_2D& vector2, 
            const Vector_2D::Measurement precision)
    {
#ifdef DEBUG_VECTOR_2D_IS_EQUAL
        cout << "Vector_2D is_equal begin\n";
        cout << "Vector_2D is_equal vector1 x: " << vector1.get_x() << " y: " << vector1.get_y() << "\n";
        cout << "vector_2D is_equal vector2 x: " << vector2.get_x() << " y: " << vector2.get_y() << "\n";
#endif
        Vector_2D::Measurement largest(fabs(vector1.get_x()));
        Vector_2D::Measurement error_bound(fabs(vector2.get_x()));
        if (error_bound > largest)
            largest = error_bound;
        
        error_bound = largest > 1.0 ? (largest * precision) : precision;
        
        if (fabs(vector1.get_x() - vector2.get_x()) <= error_bound)
        {
            largest = fabs(vector1.get_y());
            error_bound = fabs(vector2.get_y());
            if (error_bound > largest)
                largest = error_bound;

            error_bound = largest > 1.0 ? (largest * precision) : precision;
#ifdef DEBUG_VECTOR_2D_IS_EQUAL
            if (fabs(vector1.get_y() - vector2.get_y()) > error_bound)
            {
                // determine what multiple of DBL_EPSILON is needed to succeed
                Vector_2D::Measurement val(fabs(vector1.get_y() - vector2.get_y()));
                Vector_2D::Measurement divisor(DBL_EPSILON);

                if (largest > 1.0)
                    divisor *= largest;
                cout << "Vector_2D is_equal vector y (" << val << ") FALSE needs a precision of at least: " << (val / divisor) << " * DBL_EPSILON\n";
            }
#endif                
            return fabs(vector1.get_y() - vector2.get_y()) <= error_bound;
        }
#ifdef DEBUG_VECTOR_3D_IS_EQUAL
        else
        {
            // determine what multiple of DBL_EPSILON is needed to succeed
            Vector_2D::Measurement val(fabs(vector1.get_x() - vector2.get_x()));
            Vector_2D::Measurement divisor(DBL_EPSILON);
            if (largest > 1.0)
                divisor *= largest;
            cout << "Vector_2D is_equal vector x (" << val << ") FALSE needs a precision of at least: " << (val / divisor) << " * DBL_EPSILON\n";
        }
#endif
        return false;
    }
    
    const bool is_cp_zero(const Vector_2D& v1, const Vector_2D& v2, 
            const Vector_2D::Measurement precision)
    {
        Vector_2D::Measurement one(v1.get_x() * v2.get_y());
        Vector_2D::Measurement two(v1.get_y() * v2.get_x());
        
        Vector_2D::Measurement largest(fabs(one));
        Vector_2D::Measurement error_bound(fabs(two));
        if (error_bound > largest)
            largest = error_bound;
        
        error_bound = largest > 1.0 ? (largest * precision) : precision;
        
        return fabs(one - two) <= error_bound; // zero
    }
    
    const bool is_same_line(const Point_2D& v1_start, const Point_2D& v1_end, 
            const Point_2D& v2_start, const Point_2D& v2_end, bool& same_direction,
            const Vector_2D::Measurement precision)
    {
        if (is_equal(v1_start, v1_end, precision)) // zero length vector1
        {
            if (is_equal(v2_start, v2_end, precision)) // zero length vector2
            {
                if (is_equal(v1_start, v2_start, precision))
                {
                    same_direction = true;
                    return true;
                } // else return false
                same_direction = false;
                return false;
            } // else
            
            if (is_equal(v1_start, v2_start, precision))
            {
                same_direction = true;
                return true;
            } // else return if v1_start is on v2 vector
            
            if (is_cp_zero(Vector_2D(v2_start, v1_start), Vector_2D(v2_start, v2_end), precision))
            {
                same_direction = true;
                return true;
            } // else return false
            same_direction = false;
            return false;
        }
        
        if (is_equal(v2_start, v2_end, precision)) // zero length vector2
        {
            if (is_equal(v1_start, v2_start, precision))
            {
                same_direction = true;
                return true;
            } // else return if v2_start is on v1 vector
            
            if (is_cp_zero(Vector_2D(v1_start, v2_start), Vector_2D(v1_start, v1_end), precision))
            {
                same_direction = true;
                return true;
            } // else return false
            same_direction = false;
            return false;
        }
        
        Vector_2D v1(v1_start, v1_end);
        Vector_2D v2(v2_start, v2_end);
        if (is_cp_zero(v1, v2, precision))
        {
            // same general line check if they start at the same origin
            if (is_equal(v1_start, v2_start, precision))
            {
                // start at same origin
                same_direction = dot_product(v1, v2) > 0;
                return true;
            } // else
            if (is_cp_zero(Vector_2D(v1_start, v2_start), v1, precision))
            {
                same_direction = dot_product(v1, v2) > 0;
                return true;
            } // else return false
        }
        
        same_direction = false;
        return false;
    }

    const bool is_pt_on_vector(const Point_2D& pt, const Point_2D& v_start, 
            const Point_2D& v_end, const Point_2D::Measurement precision)
    {
#ifdef DEBUG_IS_PT_ON_VECTOR
        cout << "    is_pt_on_vector beginning\n";
        cout << "    is_pt_on_vector checking if v has zero length\n";
#endif
        if (is_equal(v_start, v_end, precision)) // v is zero length
            return is_equal(pt, v_start, precision);
            
#ifdef DEBUG_IS_PT_ON_VECTOR
        cout << "    is_pt_on_vector checking if pt is equal to v_start\n";
#endif
        if (is_equal(pt, v_start, precision))
            return true;
        
#ifdef DEBUG_IS_PT_ON_VECTOR
        cout << "    is_pt_on_vector checking if pt is equal to v_end\n";
#endif
        if (is_equal(pt, v_end, precision))
            return true;
        
        Vector_2D v_pt_f_start(v_start, pt);
        Vector_2D v_pt_f_end(v_end, pt);
        
#ifdef DEBUG_IS_PT_ON_VECTOR
        cout << "    is_pt_on_vector Testing from one end\n";
#endif
        Vector_2D v(v_start, v_end);

        Vector_2D::Measurement largest(fabs(v_start.get_x()));
        Vector_2D::Measurement error_bound(fabs(v_end.get_x()));
        if (error_bound > largest)
            largest = error_bound;

        error_bound = largest > 1.0 ? (largest * precision) : precision;
        if (fabs(v.get_x()) > error_bound) // if x is not zero
        {
            // ptX = oX + tX*vX
            // tX = (ptX - oX) / vX
            Vector_2D::Measurement tX((pt.get_x() - v_start.get_x()) / v.get_x());
#ifdef DEBUG_IS_PT_ON_VECTOR
            cout << "    is_pt_on_vector tX: " << tX << " " << (tX > 0 && tX < 1) << " ip x: " << ip.get_x() << " y: " << ip.get_y() << " z: " << ip.get_z() << " pt x: " << pt.get_x() << " y: " << pt.get_y() << " z: " << pt.get_z() << "\n";
#endif
            if (tX > 0 && tX < 1 && is_equal(v_start + v*tX, pt, precision))
                return true;
        }

        largest = fabs(v_start.get_y());
        error_bound = fabs(v_end.get_y());
        if (error_bound > largest)
            largest = error_bound;

        error_bound = largest > 1.0 ? (largest * precision) : precision;
        if (fabs(v.get_y()) > error_bound) // if y is not zero
        {
            // ptY = oY + tY*vY
            // tY = (ptY - oY) / vY
            Vector_2D::Measurement tY((pt.get_y() - v_start.get_y()) / v.get_y());
//            Point_2D ip(v_start + v*tY);
#ifdef DEBUG_IS_PT_ON_VECTOR
            cout << "    is_pt_on_vector tY: " << tY << " " << (tY > 0 && tY < 1) << " ip x: " << ip.get_x() << " y: " << ip.get_y() << " pt x: " << pt.get_x() << " y: " << pt.get_y() << "\n";
#endif
            if (tY > 0 && tY < 1 && is_equal(v_start + v*tY, pt, precision))
                return true;
        }

#ifdef DEBUG_IS_PT_ON_VECTOR
        cout << "    is_pt_on_vector Testing from other end\n";
#endif
        v = Vector_2D(v_end, v_start);

        largest = fabs(v_start.get_x());
        error_bound = fabs(v_end.get_x());
        if (error_bound > largest)
            largest = error_bound;

        error_bound = largest > 1.0 ? (largest * precision) : precision;
        if (fabs(v.get_x()) > error_bound) // if x is not zero
        {
            // ptX = oX + tX*vX
            // tX = (ptX - oX) / vX
            Vector_2D::Measurement tX((pt.get_x() - v_end.get_x()) / v.get_x());
//            Point_3D ip(v_end + v*tX);
#ifdef DEBUG_IS_PT_ON_VECTOR
            cout << "    is_pt_on_vector tX: " << tX << " " << (tX > 0 && tX < 1) << " ip x: " << ip.get_x() << " y: " << ip.get_y() << " pt x: " << pt.get_x() << " y: " << pt.get_y() << "\n";
#endif
            if (tX > 0 && tX < 1 && is_equal(v_end + v*tX, pt, precision))
                return true;
        }

        largest = fabs(v_start.get_y());
        error_bound = fabs(v_end.get_y());
        if (error_bound > largest)
            largest = error_bound;

        error_bound = largest > 1.0 ? (largest * precision) : precision;
        if (fabs(v.get_y()) > error_bound) // if y is not zero
        {
            // ptY = oY + tY*vY
            // tY = (ptY - oY) / vY
            Vector_2D::Measurement tY((pt.get_y() - v_end.get_y()) / v.get_y());
//            Point_3D ip(v_end + v*tY);
#ifdef DEBUG_IS_PT_ON_VECTOR
            cout << "    is_pt_on_vector tY: " << tY << " " << (tY > 0 && tY < 1) << " ip x: " << ip.get_x() << " y: " << ip.get_y() << " pt x: " << pt.get_x() << " y: " << pt.get_y() << "\n";
#endif
            if (tY > 0 && tY < 1 && is_equal(v_end + v*tY, pt, precision))
                return true;
        }
        
#ifdef DEBUG_IS_PT_ON_VECTOR
        cout << "    is_pt_on_vector returning false\n";
#endif
        return false;
    }
    
    const bool intersect_vectors_sl(const Point_2D& v1_start, const Point_2D& v1_end, 
            const Point_2D& v2_start, const Point_2D& v2_end, const bool same_direction,
            Vector_2D_idata& result, const Vector_2D::Measurement precision)
    {
#ifdef DEBUG_INTERSECT_VECTORS
        cout << "intersect_vectors_sl same_line\n";
#endif
        if (same_direction)
        {
            if (is_equal(v1_start, v2_start, precision)) // vector from v1_start to v2_start has zero length
            {
                if (is_pt_on_vector(v2_end, v1_start, v1_end, precision))
                {
                    result.p2 = v2_end;
                    result.num = 2;
                }
                else if (is_pt_on_vector(v1_end, v2_start, v2_end, precision))
                {
                    result.p2 = v1_end;
                    result.num = 2;
                }
                else
                    result.num = 1;

                result.p1 = v1_start;
                return true;
            }

            // if v1 comes before v2
            if (dot_product(Vector_2D(v1_start, v1_end), Vector_2D(v1_start, v2_start)) > 0)
            {
                if (is_pt_on_vector(v2_start, v1_start, v1_end, precision))
                {
                    if (is_pt_on_vector(v2_end, v1_start, v1_end, precision))
                    {
                        result.p2 = v2_end;
                        result.num = 2;
                    }
                    else if (is_pt_on_vector(v1_end, v2_start, v2_end, precision) && !is_equal(v1_end, v2_start, precision))
                    {
                        result.p2 = v1_end;
                        result.num = 2;
                    }
                    else
                        result.num = 1;

                    result.p1 = v2_start;
                    return true;
                }

                return false;
            }
            else // v2 comes before v1
            {
                if (is_pt_on_vector(v1_start, v2_start, v2_end, precision))
                {
                    if (is_pt_on_vector(v1_end, v2_start, v2_end, precision))
                    {
                        result.p2 = v1_end;
                        result.num = 2;
                    }
                    else if (is_pt_on_vector(v2_end, v1_start, v1_end, precision) && !is_equal(v2_end, v1_start, precision))
                    {
                        result.p2 = v2_end;
                        result.num = 2;
                    }
                    else
                        result.num = 1;

                    result.p1 = v1_start;
                    return true;
                }

                return false;
            }
        }
        else
        {
            if (is_equal(v1_start, v2_end, precision)) // vector from v1_start to v2_end has zero length
            {
                if (is_pt_on_vector(v2_start, v1_start, v1_end, precision))
                {
                    result.p2 = v2_start;
                    result.num = 2;
                }
                else if (is_pt_on_vector(v1_end, v2_start, v2_end, precision))
                {
                    result.p2 = v1_end;
                    result.num = 2;
                }
                else
                    result.num = 1;

                result.p1 = v1_start;
                return true;
            }

            // if v1 comes before v2
            if (dot_product(Vector_2D(v1_start, v1_end), Vector_2D(v1_start, v2_end)) > 0)
            {
                if (is_pt_on_vector(v2_end, v1_start, v1_end, precision))
                {
                    if (is_pt_on_vector(v2_start, v1_start, v1_end, precision))
                    {
                        result.p2 = v2_start;
                        result.num = 2;
                    }
                    else if (is_pt_on_vector(v1_end, v2_start, v2_end, precision) && !is_equal(v1_end, v2_end, precision))
                    {
                       result.p2 = v1_end;
                        result.num = 2;
                    }
                    else
                        result.num = 1;

                    result.p1 = v2_end;
                    return true;
                }

                return false;
            }
            else // v2 comes before v1
            {
                if (is_pt_on_vector(v1_start, v2_start, v2_end, precision))
                {
                    if (is_pt_on_vector(v1_end, v2_start, v2_end, precision))
                    {
                        result.p2 = v1_end;
                        result.num = 2;
                    }
                    else if (is_pt_on_vector(v2_start, v1_start, v1_end, precision) && !is_equal(v2_start, v1_start, precision))
                    {
                        result.p2 = v2_start;
                        result.num = 2;
                    }
                    else
                        result.num = 1;

                    result.p1 = v1_start;
                    return true;
                }

                return false;
            }
        }
    }
    
    const bool intersect_vectors_dl_cpz(const Point_2D& v1_start, const Point_2D& v1_end, 
            const Point_2D& v2_start, const Point_2D& v2_end, Point_2D& i_point, 
            const Vector_2D::Measurement precision)
    {
#ifdef DEBUG_INTERSECT_VECTORS
        cout << "intersect_vectors_dl_cpz begin\n";
#endif
        // vectors cross, but may not actually intersect
        // find the point on v2 that is also on v1
        // point on v2 is defined as
        // v1 = o1 + t*v1
        // so
        // pX = o1X + t*v1X
        // pY = o1Y + t*v1Y
        // where t is from 0 to 1
        //
        // cross product will be of V3 and v2 which should equal zero
        // v3 is defined as
        // v3 = p - o2
        // one equation must be zero
        // cp = (v2X * v3Y) - (v2Y * v3X) = 0
        //
        // (v2X * v3Y) - (v2Y * v3X) = 0
        // (v2X * (o1Y + t*v1Y - o2Y)) - (v2Y * (o1X + t*v1X - o2X)) = 0
        // v2X*o1Y + t*v1Y*v2X - o2Y*v2X - v2Y*o1X - t*v1X*v2Y + o2X*v2Y = 0
        // t(v1Y*v2X - v1X*v2Y) = o2Y*v2X + v2Y*o1X - v2X*o1Y - o2X*v2Y
        // t = (o2Y*v2X + v2Y*o1X - v2X*o1Y - o2X*v2Y) / (v1Y*v2X - v1X*v2Y)
        // t = (v2X*(o2Y - o1Y) + v2Y*(o1X - o2X)) / (v1Y*v2X - v1X*v2Y)
        //
        
        Vector_2D v1(v1_start, v1_end);
        Vector_2D v2(v2_start, v2_end);
        // t = (v2X*(o2Y - o1Y) + v2Y*(o1X - o2X)) / (v1Y*v2X - v1X*v2Y)
        // check for zero division
        Vector_2D::Measurement bottom(v1.get_y() * v2.get_x());
        Vector_2D::Measurement largest(fabs(bottom));
        Vector_2D::Measurement next(v1.get_x() * v2.get_y());
        bottom -= next;
        Vector_2D::Measurement error_bound(fabs(next));
        if (error_bound > largest)
            largest = error_bound;
        error_bound = largest > 1.0 ? (largest * precision) : precision;
//        cout << "intersect_vectors_dl_cpz fabs(bottom): " << fabs(bottom) << " error_bound: " << error_bound << "\n";
        if (fabs(bottom) > error_bound)
        {
            Vector_2D::Measurement t(v2.get_x() * (v2_start.get_y() - v1_start.get_y()) + v2.get_y() * (v1_start.get_x() - v2_start.get_x()));
            t /= bottom;
#ifdef DEBUG_INTERSECT_VECTORS
            cout << "intersect_vectors_dl_cpz bottom: " << bottom << " t: " << t << "\n";
#endif
            if (t >= 0 && t <= 1)
            {
                // point is on v1, now check if it is in v2
                i_point = v1_start + v1 * t;
#ifdef DEBUG_INTERSECT_VECTORS
                cout << "intersect_vectors_dl_cpz i_point_cpz_z x: " << i_point.get_x() << " y: " << i_point.get_y() << " z: " << i_point.get_z() << "\n";
                cout << "intersect_vectors_dl_cpz is_pt_on_vector1: " << is_pt_on_vector(i_point, v1_start, v1_end, precision) << "\n";
                cout << "intersect_vectors_dl_cpz is_pt_on_vector2: " << is_pt_on_vector(i_point, v2_start, v2_end, precision) << "\n";
#endif
                if (is_equal(i_point, v2_start, precision))
                {
                    i_point = v2_start;
                    if (is_pt_on_vector(i_point, v1_start, v1_end, precision))
                        return true;
                }
                else if (is_equal(i_point, v2_end, precision))
                {
                    i_point = v2_end;
                    if (is_pt_on_vector(i_point, v1_start, v1_end, precision))
                        return true;
                }
                else if (is_equal(i_point, v1_start, precision))
                {
                    i_point = v1_start;
                    if (is_pt_on_vector(i_point, v2_start, v2_end, precision))
                        return true;
                }
                else if (is_equal(i_point, v1_end, precision))
                {
                    i_point = v1_end;
                    if (is_pt_on_vector(i_point, v2_start, v2_end, precision))
                        return true;
                }
                else if (is_pt_on_vector(i_point, v1_start, v1_end, precision) && is_pt_on_vector(i_point, v2_start, v2_end, precision))
                    return true;
            }
        }
        
#ifdef DEBUG_INTERSECT_VECTORS
        cout << "intersect_vectors_dl_cpz returning false\n";
#endif
        return false;
    }
    
    const bool intersect_vectors_dl(const Point_2D& v1_start, 
            const Point_2D& v1_end, const Point_2D& v2_start, 
            const Point_2D& v2_end, Vector_2D_idata& result, 
            const Vector_2D::Measurement precision)
    {
#ifdef DEBUG_INTERSECT_VECTORS
        cout << "intersect_vectors_dl begin\n";
        cout << "intersect_vectors_dl testing if end points of vectors are on the other vector\n";
#endif
        if (is_pt_on_vector(v1_start, v2_start, v2_end, precision))
        {
            result.p1 = v1_start;
            result.num = 1;
            return true;
        }
        else if (is_pt_on_vector(v1_end, v2_start, v2_end, precision))
        {
            result.p1 = v1_end;
            result.num = 1;
            return true;
        }
        else if (is_pt_on_vector(v2_start, v1_start, v1_end, precision))
        {
            result.p1 = v2_start;
            result.num = 1;
            return true;
        }
        else if (is_pt_on_vector(v2_end, v1_start, v1_end, precision))
        {
            result.p1 = v2_end;
            result.num = 1;
            return true;
        }

        // cross product equal to zero method
        Point_2D i_point(0,0);
        if (intersect_vectors_dl_cpz(v1_start, v1_end, v2_start, v2_end, i_point, precision))
        {
            result.p1 = i_point;
            result.num = 1;
            return true;
        }
#ifdef DEBUG_INTERSECT_VECTORS
        cout << "intersect_vectors_dl Testing cross product zero method v1 negative\n";
#endif
        if (intersect_vectors_dl_cpz(v1_end, v1_start, v2_start, v2_end, i_point, precision))
        {
            result.p1 = i_point;
            result.num = 1;
            return true;
        }
#ifdef DEBUG_INTERSECT_VECTORS
        cout << "intersect_vectors_dl Testing cross product zero method v2 negative\n";
#endif
        if (intersect_vectors_dl_cpz(v1_start, v1_end, v2_end, v2_start, i_point, precision))
        {
            result.p1 = i_point;
            result.num = 1;
            return true;
        }
#ifdef DEBUG_INTERSECT_VECTORS
        cout << "intersect_vectors_dl Testing cross product zero method v1 & v2 negative\n";
#endif
        if (intersect_vectors_dl_cpz(v1_end, v1_start, v2_end, v2_start, i_point, precision))
        {
            result.p1 = i_point;
            result.num = 1;
            return true;
        }
        
#ifdef DEBUG_INTERSECT_VECTORS
        cout << "intersect_vectors_dl returning false\n";
#endif
        return false;
    }
    
    const bool intersect_vectors(const Point_2D& v1_start, const Point_2D& v1_end,
            const Point_2D& v2_start, const Point_2D& v2_end, Vector_2D_idata& result, 
            const Vector_2D::Measurement precision)
    {
#ifdef DEBUG_INTERSECT_VECTORS
        cout << "intersect_vector_vector begin\n";
        cout << "intersect_vector_vector v1_start x: " << v1_start.get_x() << " y: " << v1_start.get_y() << " z: " << v1_start.get_z() << "\n";
        cout << "intersect_vector_vector v1_end x: " << v1_end.get_x() << " y: " << v1_end.get_y() << " z: " << v1_end.get_z() << "\n";
        cout << "intersect_vector_vector v2_start x: " << v2_start.get_x() << " y: " << v2_start.get_y() << " z: " << v2_start.get_z() << "\n";
        cout << "intersect_vector_vector v2_end x: " << v2_end.get_x() << " y: " << v2_end.get_y() << " z: " << v2_end.get_z() << "\n";
#endif

        if (is_equal(v1_start, v1_end, precision)) // v1 has zero length
        {
            if (is_equal(v2_start, v2_end, precision)) // v2 has zero length
            {
                if (is_equal(v1_start, v2_start, precision))
                {
                    result.p1 = v1_start;
                    result.num = 1;
                    return true;
                }
                
                return false;
            } // else v2 has length
            
            if (is_pt_on_vector(v1_start, v2_start, v2_end, precision))
            {
                result.p1 = v1_start;
                result.num = 1;
                return true;
            }
            
            return false;
        }
        
        if (is_equal(v2_start, v2_end, precision)) // v2 has zero length
        {
            if (is_pt_on_vector(v2_start, v1_start, v1_end, precision))
            {
                result.p1 = v2_start;
                result.num = 1;
                return true;
            }
            
            return false;
        }
        
        bool same_direction(false);
        if (is_same_line(v1_start, v1_end, v2_start, v2_end, same_direction, precision))
            return intersect_vectors_sl(v1_start, v1_end, v2_start, v2_end, same_direction, result, precision);
        else // different lines
            return intersect_vectors_dl(v1_start, v1_end, v2_start, v2_end, result, precision);
//        
//        idata.num = 0; // initialize idata
//        
//        // determine end points
//        Point_2D ep1 = o1 + v1;
//        Point_2D ep2 = o2 + v2;
//
////        cout << "o1 x: " << o1.get_x() << " y: " << o1.get_y() << " ep1 x: " << ep1.get_x() << " y: " << ep1.get_y() << "\n";
////        cout << "o2 x: " << o2.get_x() << " y: " << o2.get_y() << " ep2 x: " << ep2.get_x() << " y: " << ep2.get_y() << "\n";
//        
//        // determine if vectors are in the same line
//        Vector_2D vo(o2,o1);
//        Vector_2D vep(o2,ep1);
//        if (cross_product_zero(v2, vo, precision) && 
//                cross_product_zero(v2, vep, precision)) // same line
//        {
//            // maximum of two intersection points to define the overlap - if any
//            if (dot_product(v2, v1) > 0) // v1 in same direction as v2
//            {
//                if (dot_product(v2, vo) >= 0) // o1 in direction of v2
//                {
//                    if (point_on_vector(o1, v2, o2, precision))
//                    {
//                        idata.p1 = o1;
//                        idata.num++;
//                        // see if either end point is common
//                        if (!within_round(ep1, o1, precision) && point_on_vector(ep1, v2, o2, precision)) // v1 ep1 is on v2
//                        {
//                            idata.p2 = ep1;
//                            idata.num++;
//                        }
//                        else if (!within_round(ep2, o1, precision) && point_on_vector(ep2, v1, o1, precision)) // v2 ep2 is on v1
//                        {
//                            idata.p2 = ep2;
//                            idata.num++;
//                        }
//                    } // else no other intersection possible
//                }
//                else // o1 not in direction of v2
//                {
//                    if (dot_product(v2, vep) >= 0) // ep1 past o2 in v2 direction
//                    {
//                        if (point_on_vector(ep1, v2, o2, precision))
//                        {
//                            idata.p1 = ep1;
//                            idata.num++;
//                            if (!within_round(o2, ep1, precision) && point_on_vector(o2, v1, o1, precision))
//                            {
//                                idata.p2 = o2;
//                                idata.num++;
//                            }
//                        }
//                        else // v1 completely covers v2
//                        {
//                            if (point_on_vector(o2, v1, o1, precision))
//                            {
//                                idata.p1 = o2;
//                                idata.num++;
//                            }
//                            if (!within_round(ep2, o2, precision) && point_on_vector(ep2, v1, o1, precision))
//                            {
//                                if (idata.num == 0)
//                                    idata.p1 = ep2;
//                                else
//                                    idata.p2 = ep2;
//                                idata.num++;
//                            }
//                        }
//                    }
//                }
//            }
//            else // v1 in opposite direction as v2
//            {
//                if (dot_product(v2, vep) >= 0) // ep1 in v2 direction
//                {
//                    if (point_on_vector(ep1, v2, o2, precision))
//                    {
//                        idata.p1 = ep1;
//                        idata.num++;
//                        if (!within_round(o1, ep1, precision) && point_on_vector(o1, v2, o2, precision))
//                        {
//                            idata.p2 = o1;
//                            idata.num++;
//                        }
//                        else if (!within_round(ep2, ep1, precision) && point_on_vector(ep2, v1, o1, precision))
//                        {
//                            idata.p2 = ep2;
//                            idata.num++;
//                        }
//                    }
//                }
//                else // ep1 not in v2 direction
//                {
//                    if (dot_product(v2, vo) >= 0) // o1 past o2 in v2 direction
//                    {
//                        if (point_on_vector(o1, v2, o2, precision))
//                        {
//                            idata.p1 = o1;
//                            idata.num++;
//                            if (!within_round(o2, o1, precision) && point_on_vector(o2, v1, o1, precision))
//                            {
//                                idata.p2 = o2;
//                                idata.num++;
//                            }
//                        }
//                        else // v1 completely covers v2, so test if o2 is on v1
//                        {
//                            if (point_on_vector(o2, v1, o1, precision))
//                            {
//                                idata.p1 = o2;
//                                idata.num++;
//                            }
//                            if (!within_round(ep2, o2, precision) && point_on_vector(ep2, v1, o1, precision))
//                            {
//                                if (idata.num == 0)
//                                    idata.p1 = ep2;
//                                else
//                                    idata.p2 = ep2;
//                                idata.num++;
//                            }
//                        }
//                    }
//                }
//            }
//        }
//        else // not in same line
//        {
//            // vectors cross, but may not intersect
//            // cross product:
//            // (v2X * v1Y) - (v2Y * v1X)
//            // compute cross product where vector(o2, iPoint) X v2 == 0
//            // iPoint = vo + a*v1
//            // voXv2Y + av1Xv2Y - voYv2X - av1Yv2X = 0
//            // a(v1Xv2Y - v1Yv2X) = voYv2X - voXv2Y
//            // a = (voYv2X - voXv2Y) / (v1Xv2Y - v1Yv2X)
//            
//            // maximum of one intersection points
//            Vector_2D::Measurement bottom = (v1.get_x() * v2.get_y()) - (v1.get_y() * v2.get_x());
//            if (!within_round(bottom, 0, precision))
//            {
//                Vector_2D::Measurement a = ((vo.get_y() * v2.get_x()) - (vo.get_x() * v2.get_y())) / bottom;
//                if (a >= 0 && a <= 1) // point is on v1, check if it is on v2
//                {
//                    Point_2D p(v1.get_x() * a + o1.get_x(), v1.get_y() * a + o1.get_y());
//                    vo = Vector_2D(o2, p); // reuse vo vector to check if p1 is on v2
//                    if (dot_product(v2, vo) >= 0 && vo.length() <= v2.length())
//                    {
//                        idata.p1 = p;
//                        idata.num++; // intersection point is on both vectors
//                        // check if intersection point is within round error of
//                        // origins and end points of the vectors
//                        if (within_round(idata.p1, o2, precision))
//                            idata.p1 = o2;
//                        else if (within_round(idata.p1, ep2, precision))
//                            idata.p1 = ep2;
//                        else if (within_round(idata.p1, o1, precision))
//                            idata.p1 = o1;
//                        else if (within_round(idata.p1, ep1, precision))
//                            idata.p1 = ep1;
//                    }
//                }
//            }
//        }
//        
//        return idata.num > 0;
    }

//    const bool point_on_vector(const Point_2D& point, const Vector_2D& v, const Point_2D& o, 
//            const Vector_2D::Measurement precision) 
//    {
//        Vector_2D vp(o, point);
//        if (cross_product_zero(v, vp, precision))
//        {
//            return (dot_product(v, vp) >= 0 && vp.length() <= v.length());
//        }
//        else
//            return false;
//    }
//
//    
//    const bool point_btwn_vectors(const Point_2D& p, const Vector_2D& v1, const Vector_2D& v2, 
//            const Point_2D& o, const Vector_2D::Measurement precision)
//    {
//        if (p == o)
//            return true;
//        
//        Vector_2D vop(o,p);
//        if (cross_product(v1, v2) > 0)
//        {
//            return ((cross_product(v1, vop) >= 0 || cross_product_zero(v1,vop, precision)) &&
//                    (cross_product(v2, vop) <= 0 || cross_product_zero(v2,vop, precision)));
//        }
//        else
//        {
//            return ((cross_product(v1, vop) <= 0 || cross_product_zero(v1,vop, precision)) &&
//                    (cross_product(v2, vop) >= 0 || cross_product_zero(v2,vop, precision)));
//        }
//    }
}
