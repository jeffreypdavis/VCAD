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
 * File:   Vector_3D.cpp
 * Author: Jeffrey Davis
 */

#include "Vector_3D.h"
#include <stdexcept>
#include <cfloat>

namespace VCAD_lib
{
    Vector_3D::Vector_3D(const Measurement x_val, const Measurement y_val, const Measurement z_val) 
            : x(x_val), y(y_val), z(z_val) {}
    
    Vector_3D::Vector_3D(const Point_3D& end_point) 
                : x(end_point.get_x()), y(end_point.get_y()), z(end_point.get_z()) {}
    
    Vector_3D::Vector_3D(const Point_3D& start_point, const Point_3D& end_point) :
                x(end_point.get_x() - start_point.get_x()),
                y(end_point.get_y() - start_point.get_y()),
                z(end_point.get_z() - start_point.get_z()) {}
    
    Vector_3D::Measurement Vector_3D::length() const
    {
        return sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
    }

    Vector_3D& Vector_3D::normalize()
    {
        Measurement length = this->length();
        if (length == 0)
            throw length_error("Cannot normalize a zero length vector");
        x /= length;
        y /= length;
        z /= length;
        return *this;
    }
    
    Vector_3D& Vector_3D::operator+=(const Vector_3D& v)
    {
        x += v.get_x();
        y += v.get_y();
        z += v.get_z();
        return *this;
    }
    
    Vector_3D& Vector_3D::operator-=(const Vector_3D& v)
    {
        x -= v.get_x();
        y -= v.get_y();
        z -= v.get_z();
        return *this;
    }
    
    Vector_3D& Vector_3D::operator*=(const Vector_3D::Measurement val)
    {
        x *= val;
        y *= val;
        z *= val;
        return *this;
    }

    const Vector_3D cross_product(const Vector_3D& v1, const Vector_3D& v2)
    {
        return Vector_3D(v1.get_y() * v2.get_z() - v1.get_z() * v2.get_y(),
                v1.get_z() * v2.get_x() - v1.get_x() * v2.get_z(),
                v1.get_x() * v2.get_y() - v1.get_y() * v2.get_x());
    }

    const Vector_3D::Measurement dot_product(const Vector_3D& v1, const Vector_3D& v2)
    {
        return v1.get_x() * v2.get_x() + v1.get_y() * v2.get_y() + v1.get_z() * v2.get_z();
    }

    const Vector_3D::Angle_Meas angle_between(const Vector_3D& v1, const Vector_3D& v2)
    {
        return acos(dot_product(v1, v2) / (v1.length() * v2.length()));
    }
    
    const Vector_3D orthogonal_projection(const Vector_3D& a, const Vector_3D& b)
    {
        Vector_3D::Measurement length = b.length();
        if (length == 0)
            throw length_error("Cannot orthogonal project a zero length vector");
        return b * (dot_product(a, b) / pow(length, 2));
    }
    
    const Vector_3D operator+(const Vector_3D& v1, const Vector_3D& v2)
    {
        return Vector_3D(v1.get_x() + v2.get_x(), v1.get_y() + v2.get_y(), v1.get_z() + v2.get_z());
    }
    
    const Vector_3D operator-(const Vector_3D& v1, const Vector_3D& v2)
    {
        return Vector_3D(v1.get_x() - v2.get_x(), v1.get_y() - v2.get_y(), v1.get_z() - v2.get_z());
    }
    
    const Vector_3D operator*(const Vector_3D& v, const Vector_3D::Measurement scalar)
    {
        return Vector_3D(v.get_x() * scalar, v.get_y() * scalar, v.get_z() * scalar);
    }
    
    const Vector_3D operator*(const Vector_3D::Measurement scalar, const Vector_3D& v)
    {
        return Vector_3D(v.get_x() * scalar, v.get_y() * scalar, v.get_z() * scalar);
    }
    
    const bool is_equal(const Vector_3D& vector1, const Vector_3D& vector2, 
            const Vector_3D::Measurement precision)
    {
#ifdef DEBUG_VECTOR_3D_IS_EQUAL
        cout << "Vector_3D is_equal begin\n";
        cout << "Vector_3D is_equal vector1 x: " << vector1.get_x() << " y: " << vector1.get_y() << " z: " << vector1.get_z() << "\n";
        cout << "vector_3D is_equal vector2 x: " << vector2.get_x() << " y: " << vector2.get_y() << " z: " << vector2.get_z() << "\n";
#endif
        Vector_3D::Measurement largest(fabs(vector1.get_x()));
        Vector_3D::Measurement error_bound(fabs(vector2.get_x()));
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
            
            if (fabs(vector1.get_y() - vector2.get_y()) <= error_bound)
            {
                largest = fabs(vector1.get_z());
                error_bound = fabs(vector2.get_z());
                if (error_bound > largest)
                    largest = error_bound;

                error_bound = largest > 1.0 ? (largest * precision) : precision;

#ifdef DEBUG_VECTOR_3D_IS_EQUAL
                if (fabs(vector1.get_z() - vector2.get_z()) > error_bound)
                {
                    // determine what multiple of DBL_EPSILON is needed to succeed
                    Vector_3D::Measurement val(fabs(vector1.get_z() - vector2.get_z()));
                    Vector_3D::Measurement divisor(DBL_EPSILON);
                    
                    if (largest > 1.0)
                        divisor *= largest;
                    cout << "Vector_3D is_equal vector z (" << val << ") FALSE needs a precision of at least: " << (val / divisor) << " * DBL_EPSILON\n";
                }
#endif                
                
                return fabs(vector1.get_z() - vector2.get_z()) <= error_bound;
            }
#ifdef DEBUG_VECTOR_3D_IS_EQUAL
            else
            {
                // determine what multiple of DBL_EPSILON is needed to succeed
                Vector_3D::Measurement val(fabs(vector1.get_y() - vector2.get_y()));
                Vector_3D::Measurement divisor(DBL_EPSILON);
                if (largest > 1.0)
                    divisor *= largest;
                cout << "Vector_3D is_equal vector y (" << val << ") FALSE needs a precision of at least: " << (val / divisor) << " * DBL_EPSILON\n";
            }
#endif
        }
#ifdef DEBUG_VECTOR_3D_IS_EQUAL
        else
        {
            // determine what multiple of DBL_EPSILON is needed to succeed
            Vector_3D::Measurement val(fabs(vector1.get_x() - vector2.get_x()));
            Vector_3D::Measurement divisor(DBL_EPSILON);
            if (largest > 1.0)
                divisor *= largest;
            cout << "Vector_3D is_equal vector x (" << val << ") FALSE needs a precision of at least: " << (val / divisor) << " * DBL_EPSILON\n";
        }
#endif
        return false;
    }
    
    const bool is_cp_zero(const Vector_3D& v1, const Vector_3D& v2, 
            const Vector_3D::Measurement precision)
    {
        Vector_3D::Measurement one(v1.get_y() * v2.get_z());
        Vector_3D::Measurement two(v1.get_z() * v2.get_y());
        
        Vector_3D::Measurement largest(fabs(one));
        Vector_3D::Measurement error_bound(fabs(two));
        if (error_bound > largest)
            largest = error_bound;
        
        error_bound = largest > 1.0 ? (largest * precision) : precision;
        
        if (fabs(one - two) <= error_bound) // zero
        {
            one = v1.get_x() * v2.get_z();
            two = v1.get_z() * v2.get_x();

            largest = fabs(one);
            error_bound = fabs(two);
            if (error_bound > largest)
                largest = error_bound;

            error_bound = largest > 1.0 ? (largest * precision) : precision;

            if (fabs(two - one) <= error_bound) // zero
            {
                one = v1.get_x() * v2.get_y();
                two = v1.get_y() * v2.get_x();

                largest = fabs(one);
                error_bound = fabs(two);
                if (error_bound > largest)
                    largest = error_bound;

                error_bound = largest > 1.0 ? (largest * precision) : precision;

#ifdef DEBUG_IS_CP_3D_ZERO
                if (fabs(one - two) > error_bound)
                {
                    Vector_3D::Measurement val(fabs(one - two));
                    Vector_3D::Measurement divisor(DBL_EPSILON);
                    if (largest > 1.0)
                        divisor *= largest;
                    cout << "is_cp_zero result z(" << val << ") compared to (" << error_bound << ") one (" << one << ") two (" << two << ") needs to have a precision of " << (val / divisor) << " * DBL_EPSILON\n";
                }
#endif
                return fabs(one - two) <= error_bound;
            }
#ifdef DEBUG_IS_CP_3D_ZERO
            else
            {
                Vector_3D::Measurement val(fabs(two - one));
                Vector_3D::Measurement divisor(DBL_EPSILON);
                if (largest > 1.0)
                    divisor *= largest;
                cout << "is_cp_zero result y(" << val << ") compared to (" << error_bound << ") one (" << one << ") two (" << two << ") needs to have a precision of " << (val / divisor) << " * DBL_EPSILON\n";
            }
#endif
        }
#ifdef DEBUG_IS_CP_3D_ZERO
        else
        {
            Vector_3D::Measurement val(fabs(one - two));
            Vector_3D::Measurement divisor(DBL_EPSILON);
            if (largest > 1.0)
                divisor *= largest;
            cout << "is_cp_zero result x(" << val << ") compared to (" << error_bound << ") one (" << one << ") two (" << two << ") needs to have a precision of " << (val / divisor) << " * DBL_EPSILON\n";
        }
        cout << "is_cp_zero returning false\n";
#endif        
        return false;
    }
    
    const bool is_same_line(const Point_3D& v1_start, const Point_3D& v1_end, 
            const Point_3D& v2_start, const Point_3D& v2_end, bool& same_direction,
            const Vector_3D::Measurement precision)
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
            
            if (is_cp_zero(Vector_3D(v2_start, v1_start), Vector_3D(v2_start, v2_end), precision))
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
            
            if (is_cp_zero(Vector_3D(v1_start, v2_start), Vector_3D(v1_start, v1_end), precision))
            {
                same_direction = true;
                return true;
            } // else return false
            same_direction = false;
            return false;
        }
        
        Vector_3D v1(v1_start, v1_end);
        Vector_3D v2(v2_start, v2_end);
        if (is_cp_zero(v1, v2, precision))
        {
            // same general line check if they start at the same origin
            if (is_equal(v1_start, v2_start, precision))
            {
                // start at same origin
                same_direction = dot_product(v1, v2) > 0;
                return true;
            } // else
            if (is_cp_zero(Vector_3D(v1_start, v2_start), v1, precision))
            {
                same_direction = dot_product(v1, v2) > 0;
                return true;
            } // else return false
        }
        
        same_direction = false;
        return false;
    }

    const bool is_pt_on_vector(const Point_3D& pt, const Point_3D& v_start, 
            const Point_3D& v_end, const Point_3D::Measurement precision)
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
        
        Vector_3D v_pt_f_start(v_start, pt);
        Vector_3D v_pt_f_end(v_end, pt);
        
#ifdef DEBUG_IS_PT_ON_VECTOR
        cout << "    is_pt_on_vector Testing from one end\n";
#endif
        Vector_3D v(v_start, v_end);

        Vector_3D::Measurement largest(fabs(v_start.get_x()));
        Vector_3D::Measurement error_bound(fabs(v_end.get_x()));
        if (error_bound > largest)
            largest = error_bound;

        error_bound = largest > 1.0 ? (largest * precision) : precision;
        if (fabs(v.get_x()) > error_bound) // if x is not zero
        {
            // ptX = oX + tX*vX
            // tX = (ptX - oX) / vX
            Vector_3D::Measurement tX((pt.get_x() - v_start.get_x()) / v.get_x());
//            Point_3D ip(v_start + v*tX);
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
            Vector_3D::Measurement tY((pt.get_y() - v_start.get_y()) / v.get_y());
//            Point_3D ip(v_start + v*tY);
#ifdef DEBUG_IS_PT_ON_VECTOR
            cout << "    is_pt_on_vector tY: " << tY << " " << (tY > 0 && tY < 1) << " ip x: " << ip.get_x() << " y: " << ip.get_y() << " z: " << ip.get_z() << " pt x: " << pt.get_x() << " y: " << pt.get_y() << " z: " << pt.get_z() << "\n";
#endif
            if (tY > 0 && tY < 1 && is_equal(v_start + v*tY, pt, precision))
                return true;
        }

        largest = fabs(v_start.get_z());
        error_bound = fabs(v_end.get_z());
        if (error_bound > largest)
            largest = error_bound;

        error_bound = largest > 1.0 ? (largest * precision) : precision;
        if (fabs(v.get_z()) > error_bound) // if z is not zero
        {
            // ptZ = oZ + tZ*vZ
            // tZ = (ptZ - oZ) / vZ
            Vector_3D::Measurement tZ((pt.get_z() - v_start.get_z()) / v.get_z());
//            Point_3D ip(v_start + v*tZ);
#ifdef DEBUG_IS_PT_ON_VECTOR
            cout << "    is_pt_on_vector tZ: " << tZ << " " << (tZ > 0 && tZ < 1) << " ip x: " << ip.get_x() << " y: " << ip.get_y() << " z: " << ip.get_z() << " pt x: " << pt.get_x() << " y: " << pt.get_y() << " z: " << pt.get_z() << "\n";
#endif
            if (tZ > 0 && tZ < 1 && is_equal(v_start + v*tZ, pt, precision))
                return true;
        }

#ifdef DEBUG_IS_PT_ON_VECTOR
        cout << "    is_pt_on_vector Testing from other end\n";
#endif
        v = Vector_3D(v_end, v_start);

        largest = fabs(v_start.get_x());
        error_bound = fabs(v_end.get_x());
        if (error_bound > largest)
            largest = error_bound;

        error_bound = largest > 1.0 ? (largest * precision) : precision;
        if (fabs(v.get_x()) > error_bound) // if x is not zero
        {
            // ptX = oX + tX*vX
            // tX = (ptX - oX) / vX
            Vector_3D::Measurement tX((pt.get_x() - v_end.get_x()) / v.get_x());
//            Point_3D ip(v_end + v*tX);
#ifdef DEBUG_IS_PT_ON_VECTOR
            cout << "    is_pt_on_vector tX: " << tX << " " << (tX > 0 && tX < 1) << " ip x: " << ip.get_x() << " y: " << ip.get_y() << " z: " << ip.get_z() << " pt x: " << pt.get_x() << " y: " << pt.get_y() << " z: " << pt.get_z() << "\n";
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
            Vector_3D::Measurement tY((pt.get_y() - v_end.get_y()) / v.get_y());
//            Point_3D ip(v_end + v*tY);
#ifdef DEBUG_IS_PT_ON_VECTOR
            cout << "    is_pt_on_vector tY: " << tY << " " << (tY > 0 && tY < 1) << " ip x: " << ip.get_x() << " y: " << ip.get_y() << " z: " << ip.get_z() << " pt x: " << pt.get_x() << " y: " << pt.get_y() << " z: " << pt.get_z() << "\n";
#endif
            if (tY > 0 && tY < 1 && is_equal(v_end + v*tY, pt, precision))
                return true;
        }

        largest = fabs(v_start.get_z());
        error_bound = fabs(v_end.get_z());
        if (error_bound > largest)
            largest = error_bound;

        error_bound = largest > 1.0 ? (largest * precision) : precision;
        if (fabs(v.get_z()) > error_bound) // if z is not zero
        {
            // ptZ = oZ + tZ*vZ
            // tZ = (ptZ - oZ) / vZ
            Vector_3D::Measurement tZ((pt.get_z() - v_end.get_z()) / v.get_z());
//            Point_3D ip(v_end + v*tZ);
#ifdef DEBUG_IS_PT_ON_VECTOR
            cout << "    is_pt_on_vector tZ: " << tZ << " " << (tZ > 0 && tZ < 1) << " ip x: " << ip.get_x() << " y: " << ip.get_y() << " z: " << ip.get_z() << " pt x: " << pt.get_x() << " y: " << pt.get_y() << " z: " << pt.get_z() << "\n";
#endif
            if (tZ > 0 && tZ < 1 && is_equal(v_end + v*tZ, pt, precision))
                return true;
        }
        
#ifdef DEBUG_IS_PT_ON_VECTOR
        cout << "    is_pt_on_vector returning false\n";
#endif
        return false;
    }
    
    const bool intersect_vectors_sl(const Point_3D& v1_start, const Point_3D& v1_end, 
            const Point_3D& v2_start, const Point_3D& v2_end, const bool same_direction,
            Vector_3D_idata& result, const Vector_3D::Measurement precision)
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
            if (dot_product(Vector_3D(v1_start, v1_end), Vector_3D(v1_start, v2_start)) > 0)
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
            if (dot_product(Vector_3D(v1_start, v1_end), Vector_3D(v1_start, v2_end)) > 0)
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
    
    const bool intersect_vectors_dl_t1(const Point_3D& o1, 
            const Vector_3D& v1, const Point_3D& o2, const Vector_3D& v2,
            const axis axis_a, const axis axis_b, Point_3D& i_point, 
            const Vector_3D::Measurement precision)
    {
#ifdef DEBUG_INTERSECT_VECTORS
        cout << "intersect_vectors_dl_t1 begin\n";
#endif
        // solve for t1 first
        // t2 = (v1B*(o2A - o1A) + v1A*(o1B - o2B)) / (v2B*v1A - V2A*v1B)
        // t1 = (o2B - o1B + t2*v2B) / v1B
        Vector_3D::Measurement o1A(0);
        Vector_3D::Measurement o1B(0);
        Vector_3D::Measurement v1A(0);
        Vector_3D::Measurement v1B(0);
        Vector_3D::Measurement o2A(0);
        Vector_3D::Measurement o2B(0);
        Vector_3D::Measurement v2A(0);
        Vector_3D::Measurement v2B(0);
        
        switch (axis_a) {
            case (X):
                o1A = o1.get_x();
                o2A = o2.get_x();
                v1A = v1.get_x();
                v2A = v2.get_x();
                break;
            case (Y):
                o1A = o1.get_y();
                o2A = o2.get_y();
                v1A = v1.get_y();
                v2A = v2.get_y();
                break;
            default:
                o1A = o1.get_z();
                o2A = o2.get_z();
                v1A = v1.get_z();
                v2A = v2.get_z();
                break;
        }
        
        switch (axis_b) {
            case (X):
                o1B = o1.get_x();
                o2B = o2.get_x();
                v1B = v1.get_x();
                v2B = v2.get_x();
                break;
            case (Y):
                o1B = o1.get_y();
                o2B = o2.get_y();
                v1B = v1.get_y();
                v2B = v2.get_y();
                break;
            default:
                o1B = o1.get_z();
                o2B = o2.get_z();
                v1B = v1.get_z();
                v2B = v2.get_z();
                break;
        }

        // try as provided and from opposite ends because bottom will be the same
        Vector_3D::Measurement bottom(v2B*v1A);
#ifdef DEBUG_INTERSECT_VECTORS
        cout << "intersect_vectors_dl_t1 bottom: " << bottom << "\n";
#endif
        Vector_3D::Measurement largest(fabs(bottom));
        Vector_3D::Measurement next(v2A*v1B);
#ifdef DEBUG_INTERSECT_VECTORS
        cout << "intersect_vectors_dl_t1 next: " << next << "\n";
#endif
        Vector_3D::Measurement error_bound(fabs(next));
        if (error_bound > largest)
            largest = error_bound;
        
        error_bound = largest > 1.0 ? (largest * precision) : precision;
        bottom -= next;
#ifdef DEBUG_INTERSECT_VECTORS
        cout << "intersect_vectors_dl_t1 bottom: " << bottom << "\n";
#endif

        if (fabs(bottom) > error_bound)
        {
            // try vectors as provided
            Vector_3D::Measurement t2(v1B * (o2A - o1A) + v1A * (o1B - o2B));
            t2 /= bottom;
#ifdef DEBUG_INTERSECT_VECTORS
            cout << "intersect_vectors_dl_t1 t2: " << t2 << "\n";
#endif
            if (t2 >= 0 && t2 <= 1)
            {
                Vector_3D::Measurement t1 = (o2B - o1B + t2 * v2B) / v1B;
#ifdef DEBUG_INTERSECT_VECTORS
                cout << "intersect_vectors_dl_t1 t1: " << t1 << "\n";
#endif
                if (t1 >= 0 && t1 <= 1)
                {
                    i_point = o1 + v1 * t1;
#ifdef DEBUG_INTERSECT_VECTORS
                    cout << "intersect_vectors_dl_t1 i_point_t1 x: " << i_point.get_x() << " y: " << i_point.get_y() << " z: " << i_point.get_z() << "\n";
                    cout << "intersect_vectors_dl_t1 is_pt_on_vector1: " << is_pt_on_vector(i_point, o1, o1 + v1, precision) << "\n";
                    cout << "intersect_vectors_dl_t1 is_pt_on_vector2: " << is_pt_on_vector(i_point, o2, o2 + v2, precision) << "\n";
#endif
                    if (is_pt_on_vector(i_point, o1, o1 + v1, precision) && is_pt_on_vector(i_point, o2, o2 + v2, precision))
                        return true;
                    i_point = o2 + v2 * t2;
#ifdef DEBUG_INTERSECT_VECTORS
                    cout << "intersect_vectors_dl_t1 i_point_t2 x: " << i_point.get_x() << " y: " << i_point.get_y() << " z: " << i_point.get_z() << "\n";
                    cout << "intersect_vectors_dl_t1 is_pt_on_vector1: " << is_pt_on_vector(i_point, o1, o1 + v1, precision) << "\n";
                    cout << "intersect_vectors_dl_t1 is_pt_on_vector2: " << is_pt_on_vector(i_point, o2, o2 + v2, precision) << "\n";
#endif
                    if (is_pt_on_vector(i_point, o1, o1 + v1, precision) && is_pt_on_vector(i_point, o2, o2 + v2, precision))
                        return true;
                }
            }
            
            // try vectors starting from end point
            // t2 = (v1B*(o2A - o1A) + v1A*(o1B - o2B)) / (v2B*v1A - V2A*v1B)
            // t1 = (o2B - o1B + t2*v2B) / v1B
            t2 = v1B * ((o1A + v1A) - (o2A + v2A)) + v1A * ((o2B + v2B) - (o1B + v1B));
            t2 /= bottom;
#ifdef DEBUG_INTERSECT_VECTORS
            cout << "intersect_vectors_dl_t1 t2: " << t2 << "\n";
#endif
            if (t2 >= 0 && t2 <= 1)
            {
                Vector_3D::Measurement t1 = (o1B + v1B - (o2B + v2B) + t2 * v2B) / v1B;
#ifdef DEBUG_INTERSECT_VECTORS
                cout << "intersect_vectors_dl_t1 t1: " << t1 << "\n";
#endif
                if (t1 >= 0 && t1 <= 1)
                {
                    i_point = o1 + v1 - v1 * t1;
#ifdef DEBUG_INTERSECT_VECTORS
                    cout << "intersect_vectors_dl_t1 i_point_v12_t1 x: " << i_point.get_x() << " y: " << i_point.get_y() << " z: " << i_point.get_z() << "\n";
                    cout << "intersect_vectors_dl_t1 is_pt_on_vector1: " << is_pt_on_vector(i_point, o1, o1 + v1, precision) << "\n";
                    cout << "intersect_vectors_dl_t1 is_pt_on_vector2: " << is_pt_on_vector(i_point, o2, o2 + v2, precision) << "\n";
#endif
                    if (is_pt_on_vector(i_point, o1, o1 + v1, precision) && is_pt_on_vector(i_point, o2, o2 + v2, precision))
                        return true;
                    i_point = o2 + v2 - v2 * t2;
#ifdef DEBUG_INTERSECT_VECTORS
                    cout << "intersect_vectors_dl_t1 i_point_v12_t2 x: " << i_point.get_x() << " y: " << i_point.get_y() << " z: " << i_point.get_z() << "\n";
                    cout << "intersect_vectors_dl_t1 is_pt_on_vector1: " << is_pt_on_vector(i_point, o1, o1 + v1, precision) << "\n";
                    cout << "intersect_vectors_dl_t1 is_pt_on_vector2: " << is_pt_on_vector(i_point, o2, o2 + v2, precision) << "\n";
#endif
                    if (is_pt_on_vector(i_point, o1, o1 + v1, precision) && is_pt_on_vector(i_point, o2, o2 + v2, precision))
                        return true;
                }
            }
        }
        else
            return false; // no need to calculate the next bottom because it will be zero too
        
        bottom = v2A*v1B;
#ifdef DEBUG_INTERSECT_VECTORS
        cout << "intersect_vectors_dl_t1 bottom: " << bottom << "\n";
#endif
        largest = fabs(bottom);
        next = v2B*v1A;
#ifdef DEBUG_INTERSECT_VECTORS
        cout << "intersect_vectors_dl_t1 next: " << next << "\n";
#endif
        error_bound = fabs(next);
        if (error_bound > largest)
            largest = error_bound;

        error_bound = largest > 1.0 ? (largest * precision) : precision;
        bottom -= next;
#ifdef DEBUG_INTERSECT_VECTORS
        cout << "intersect_vectors_dl_t1 bottom: " << bottom << "\n";
#endif
        
        if (fabs(bottom) > error_bound)
        {
            // try vectors with v1 starting from v1_end
            Vector_3D::Measurement t2(v1B * ((o1A + v1A) - o2A) + v1A * (o2B - (o1B + v1B)));
            t2 /= bottom;
#ifdef DEBUG_INTERSECT_VECTORS
            cout << "intersect_vectors_dl_t1 t2: " << t2 << "\n";
#endif
            if (t2 >= 0 && t2 <= 1)
            {
                Vector_3D::Measurement t1 = (o1B + v1B - o2B - t2 * v2B) / v1B;
#ifdef DEBUG_INTERSECT_VECTORS
                cout << "intersect_vectors_dl_t1 t1: " << t1 << "\n";
#endif
                if (t1 >= 0 && t1 <= 1)
                {
                    i_point = o1 + v1 - v1 * t1;
#ifdef DEBUG_INTERSECT_VECTORS
                    cout << "intersect_vectors_dl_t1 i_point_v1_t1 x: " << i_point.get_x() << " y: " << i_point.get_y() << " z: " << i_point.get_z() << "\n";
                    cout << "intersect_vectors_dl_t1 is_pt_on_vector1: " << is_pt_on_vector(i_point, o1, o1 + v1, precision) << "\n";
                    cout << "intersect_vectors_dl_t1 is_pt_on_vector2: " << is_pt_on_vector(i_point, o2, o2 + v2, precision) << "\n";
#endif
                    if (is_pt_on_vector(i_point, o1, o1 + v1, precision) && is_pt_on_vector(i_point, o2, o2 + v2, precision))
                        return true;
                    i_point = o2 + v2 * t2;
#ifdef DEBUG_INTERSECT_VECTORS
                    cout << "intersect_vectors_dl_t1 i_point_v1_t2 x: " << i_point.get_x() << " y: " << i_point.get_y() << " z: " << i_point.get_z() << "\n";
                    cout << "intersect_vectors_dl_t1 is_pt_on_vector1: " << is_pt_on_vector(i_point, o1, o1 + v1, precision) << "\n";
                    cout << "intersect_vectors_dl_t1 is_pt_on_vector2: " << is_pt_on_vector(i_point, o2, o2 + v2, precision) << "\n";
#endif
                    if (is_pt_on_vector(i_point, o1, o1 + v1, precision) && is_pt_on_vector(i_point, o2, o2 + v2, precision))
                        return true;
                }
            }
            
            // try vectors with v2 starting vrom v2_end
            // t2 = (v1B*(o2A - o1A) + v1A*(o1B - o2B)) / (v2B*v1A - V2A*v1B)
            // t1 = (o2B - o1B + t2*v2B) / v1B
            t2 = v1B * (o2A + v2A - o1A) + v1A * (o1B - (o2B + v2B));
            t2 /= bottom;
#ifdef DEBUG_INTERSECT_VECTORS
            cout << "intersect_vectors_dl_t1 t2: " << t2 << "\n";
#endif
            if (t2 >= 0 && t2 <= 1)
            {
                Vector_3D::Measurement t1 = (o2B + v2B - o1B - t2 * v2B) / v1B;
#ifdef DEBUG_INTERSECT_VECTORS
                cout << "intersect_vectors_dl_t1 t1: " << t1 << "\n";
#endif
                if (t1 >= 0 && t1 <= 1)
                {
                    i_point = o1 + v1 - v1 * t1;
#ifdef DEBUG_INTERSECT_VECTORS
                    cout << "intersect_vectors_dl_t1 i_point_v2_t1 x: " << i_point.get_x() << " y: " << i_point.get_y() << " z: " << i_point.get_z() << "\n";
                    cout << "intersect_vectors_dl_t1 is_pt_on_vector1: " << is_pt_on_vector(i_point, o1, o1 + v1, precision) << "\n";
                    cout << "intersect_vectors_dl_t1 is_pt_on_vector2: " << is_pt_on_vector(i_point, o2, o2 + v2, precision) << "\n";
#endif
                    if (is_pt_on_vector(i_point, o1, o1 + v1, precision) && is_pt_on_vector(i_point, o2, o2 + v2, precision))
                        return true;
                    i_point = o2 + v2 - v2 * t2;
#ifdef DEBUG_INTERSECT_VECTORS
                    cout << "intersect_vectors_dl_t1 i_point_v2_t2 x: " << i_point.get_x() << " y: " << i_point.get_y() << " z: " << i_point.get_z() << "\n";
                    cout << "intersect_vectors_dl_t1 is_pt_on_vector1: " << is_pt_on_vector(i_point, o1, o1 + v1, precision) << "\n";
                    cout << "intersect_vectors_dl_t1 is_pt_on_vector2: " << is_pt_on_vector(i_point, o2, o2 + v2, precision) << "\n";
#endif
                    if (is_pt_on_vector(i_point, o1, o1 + v1, precision) && is_pt_on_vector(i_point, o2, o2 + v2, precision))
                        return true;
                }
            }
        }
#ifdef DEBUG_INTERSECT_VECTORS
        else
            cout << "intersect_vectors_dl_t1 bottom is zero\n";
        cout << "intersect_vectors_dl_t1 returning false\n";
#endif
        return false;
    }
    
    const bool intersect_vectors_dl_t2(const Point_3D& o1, 
            const Vector_3D& v1, const Point_3D& o2, const Vector_3D& v2,
            const axis axis_a, const axis axis_b, Point_3D& i_point, 
            const Vector_3D::Measurement precision)
    {
#ifdef DEBUG_INTERSECT_VECTORS
        cout << "intersect_vectors_dl_t2 begin\n";
#endif
        // solve for t2 first
        // t1 = (v2B*(o2A - o1A) + v2A*(o1B - o2B) / (v1A*v2B - v1B*v2A)
        // t2 = (o1B - o2B + t1*v1B) / v2B
        Vector_3D::Measurement o1A(0);
        Vector_3D::Measurement o1B(0);
        Vector_3D::Measurement v1A(0);
        Vector_3D::Measurement v1B(0);
        Vector_3D::Measurement o2A(0);
        Vector_3D::Measurement o2B(0);
        Vector_3D::Measurement v2A(0);
        Vector_3D::Measurement v2B(0);
        
        switch (axis_a) {
            case (X):
                o1A = o1.get_x();
                o2A = o2.get_x();
                v1A = v1.get_x();
                v2A = v2.get_x();
                break;
            case (Y):
                o1A = o1.get_y();
                o2A = o2.get_y();
                v1A = v1.get_y();
                v2A = v2.get_y();
                break;
            default:
                o1A = o1.get_z();
                o2A = o2.get_z();
                v1A = v1.get_z();
                v2A = v2.get_z();
                break;
        }
        
        switch (axis_b) {
            case (X):
                o1B = o1.get_x();
                o2B = o2.get_x();
                v1B = v1.get_x();
                v2B = v2.get_x();
                break;
            case (Y):
                o1B = o1.get_y();
                o2B = o2.get_y();
                v1B = v1.get_y();
                v2B = v2.get_y();
                break;
            default:
                o1B = o1.get_z();
                o2B = o2.get_z();
                v1B = v1.get_z();
                v2B = v2.get_z();
                break;
        }

        // try as provided and from opposite ends because bottom will be the same
        Vector_3D::Measurement bottom(v2B*v1A);
#ifdef DEBUG_INTERSECT_VECTORS
        cout << "intersect_vectors_dl_t2 bottom: " << bottom << "\n";
#endif
        Vector_3D::Measurement largest(fabs(bottom));
        Vector_3D::Measurement next(v2A*v1B);
#ifdef DEBUG_INTERSECT_VECTORS
        cout << "intersect_vectors_dl_t2 next: " << next << "\n";
#endif
        Vector_3D::Measurement error_bound(fabs(next));
        if (error_bound > largest)
            largest = error_bound;
        
        error_bound = largest > 1.0 ? (largest * precision) : precision;
        bottom -= next;
#ifdef DEBUG_INTERSECT_VECTORS
        cout << "intersect_vectors_dl_t2 bottom: " << bottom << "\n";
#endif
        
        if (fabs(bottom) > error_bound)
        {
            // try vectors as provided
            // t1 = (v2B*(o2A - o1A) + v2A*(o1B - o2B) / (v1A*v2B - v1B*v2A)
            // t2 = (o1B - o2B + t1*v1B) / v2B
            Vector_3D::Measurement t1(v2B * (o2A - o1A) + v2A * (o1B - o2B));
            t1 /= bottom;
#ifdef DEBUG_INTERSECT_VECTORS
            cout << "intersect_vectors_dl_t2 t1: " << t1 << "\n";
#endif
            if (t1 >= 0 && t1 <= 1)
            {
                Vector_3D::Measurement t2((o1B - o2B + t1 * v1B) / v2B);
#ifdef DEBUG_INTERSECT_VECTORS
                cout << "intersect_vectors_dl_t2 t2: " << t2 << "\n";
#endif
                if (t2 >= 0 && t2 <= 1)
                {
                    i_point = o1 + v1 * t1;
#ifdef DEBUG_INTERSECT_VECTORS
                    cout << "intersect_vectors_dl_t2 i_point_t1 x: " << i_point.get_x() << " y: " << i_point.get_y() << " z: " << i_point.get_z() << "\n";
                    cout << "intersect_vectors_dl_t2 is_pt_on_vector1: " << is_pt_on_vector(i_point, o1, o1 + v1, precision) << "\n";
                    cout << "intersect_vectors_dl_t2 is_pt_on_vector2: " << is_pt_on_vector(i_point, o2, o2 + v2, precision) << "\n";
#endif
                    if (is_pt_on_vector(i_point, o1, o1 + v1, precision) && is_pt_on_vector(i_point, o2, o2 + v2, precision))
                        return true;
                    i_point = o2 + v2 * t2;
#ifdef DEBUG_INTERSECT_VECTORS
                    cout << "intersect_vectors_dl_t2 i_point_t2 x: " << i_point.get_x() << " y: " << i_point.get_y() << " z: " << i_point.get_z() << "\n";
                    cout << "intersect_vectors_dl_t2 is_pt_on_vector1: " << is_pt_on_vector(i_point, o1, o1 + v1, precision) << "\n";
                    cout << "intersect_vectors_dl_t2 is_pt_on_vector2: " << is_pt_on_vector(i_point, o2, o2 + v2, precision) << "\n";
#endif
                    if (is_pt_on_vector(i_point, o1, o1 + v1, precision) && is_pt_on_vector(i_point, o2, o2 + v2, precision))
                        return true;
                }
            }
            
            // try vectors starting from end point
            // t1 = (v2B*(o2A - o1A) + v2A*(o1B - o2B) / (v1A*v2B - v1B*v2A)
            // t2 = (o1B - o2B + t1*v1B) / v2B
            t1 = v2B * ((o1A + v1A) - (o2A + v2A)) + v2A * ((o2B + v2B) - (o1B + v1B));
            t1 /= bottom;
#ifdef DEBUG_INTERSECT_VECTORS
            cout << "intersect_vectors_dl_t2 t1: " << t1 << "\n";
#endif
            if (t1 >= 0 && t1 <= 1)
            {
                Vector_3D::Measurement t2((o1B + v1B - (o2B + v2B) + t1 * v1B) / v2B);
#ifdef DEBUG_INTERSECT_VECTORS
                cout << "intersect_vectors_dl_t2 t2: " << t2 << "\n";
#endif
                if (t2 >= 0 && t2 <= 1)
                {
                    i_point = o1 + v1 - v1 * t1;
#ifdef DEBUG_INTERSECT_VECTORS
                    cout << "intersect_vectors_dl_t2 i_point_v12_t1 x: " << i_point.get_x() << " y: " << i_point.get_y() << " z: " << i_point.get_z() << "\n";
                    cout << "intersect_vectors_dl_t2 is_pt_on_vector1: " << is_pt_on_vector(i_point, o1, o1 + v1, precision) << "\n";
                    cout << "intersect_vectors_dl_t2 is_pt_on_vector2: " << is_pt_on_vector(i_point, o2, o2 + v2, precision) << "\n";
#endif
                    if (is_pt_on_vector(i_point, o1, o1 + v1, precision) && is_pt_on_vector(i_point, o2, o2 + v2, precision))
                        return true;
                    i_point = o2 + v2 - v2 * t2;
#ifdef DEBUG_INTERSECT_VECTORS
                    cout << "intersect_vectors_dl_t2 i_point_v12_t2 x: " << i_point.get_x() << " y: " << i_point.get_y() << " z: " << i_point.get_z() << "\n";
                    cout << "intersect_vectors_dl_t2 is_pt_on_vector1: " << is_pt_on_vector(i_point, o1, o1 + v1, precision) << "\n";
                    cout << "intersect_vectors_dl_t2 is_pt_on_vector2: " << is_pt_on_vector(i_point, o2, o2 + v2, precision) << "\n";
#endif
                    if (is_pt_on_vector(i_point, o1, o1 + v1, precision) && is_pt_on_vector(i_point, o2, o2 + v2, precision))
                        return true;
                }
            }
        }
        else
            return false; // no need to calculate the next bottom because it will be zero too
        
        bottom = v2A*v1B;
#ifdef DEBUG_INTERSECT_VECTORS
        cout << "intersect_vectors_dl_t2 bottom: " << bottom << "\n";
#endif
        largest = fabs(bottom);
        next = v2B*v1A;
#ifdef DEBUG_INTERSECT_VECTORS
        cout << "intersect_vectors_dl_t2 next: " << next << "\n";
#endif
        error_bound = fabs(next);
        if (error_bound > largest)
            largest = error_bound;
        error_bound = largest > 1.0 ? (largest * precision) : precision;
        bottom -= next;
#ifdef DEBUG_INTERSECT_VECTORS
        cout << "intersect_vectors_dl_t2 bottom: " << bottom << "\n";
#endif
        
        if (fabs(bottom) > error_bound)
        {
            // try vectors with v1 starting from v1_end
            // t1 = (v2B*(o2A - (o1A + v1A)) + v2A*((o1B + v1B) - o2B) / (v1B*v2A - v1A*v2B)
            // t2 = (o1B + v1B - o2B - t1*v1B) / v2B
            Vector_3D::Measurement t1(v2B * (o2A - (o1A + v1A)) + v2A * (o1B + v1B - o2B));
            t1 /= bottom;
#ifdef DEBUG_INTERSECT_VECTORS
            cout << "intersect_vectors_dl_t2 t1: " << t1 << "\n";
#endif
            if (t1 >= 0 && t1 <= 1)
            {
                Vector_3D::Measurement t2((o1B + v1B - o2B - t1 * v1B) / v2B);
#ifdef DEBUG_INTERSECT_VECTORS
                cout << "intersect_vectors_dl_t2 t2: " << t2 << "\n";
#endif
                if (t2 >= 0 && t2 <= 1)
                {
                    i_point = o1 + v1 - v1 * t1;
#ifdef DEBUG_INTERSECT_VECTORS
                    cout << "intersect_vectors_dl_t2 i_point_v1_t1 x: " << i_point.get_x() << " y: " << i_point.get_y() << " z: " << i_point.get_z() << "\n";
                    cout << "intersect_vectors_dl_t2 is_pt_on_vector1: " << is_pt_on_vector(i_point, o1, o1 + v1, precision) << "\n";
                    cout << "intersect_vectors_dl_t2 is_pt_on_vector2: " << is_pt_on_vector(i_point, o2, o2 + v2, precision) << "\n";
#endif
                    if (is_pt_on_vector(i_point, o1, o1 + v1, precision) && is_pt_on_vector(i_point, o2, o2 + v2, precision))
                        return true;
                    i_point = o2 + v2 * t2;
#ifdef DEBUG_INTERSECT_VECTORS
                    cout << "intersect_vectors_dl_t2 i_point_v1_t2 x: " << i_point.get_x() << " y: " << i_point.get_y() << " z: " << i_point.get_z() << "\n";
                    cout << "intersect_vectors_dl_t2 is_pt_on_vector1: " << is_pt_on_vector(i_point, o1, o1 + v1, precision) << "\n";
                    cout << "intersect_vectors_dl_t2 is_pt_on_vector2: " << is_pt_on_vector(i_point, o2, o2 + v2, precision) << "\n";
#endif
                    if (is_pt_on_vector(i_point, o1, o1 + v1, precision) && is_pt_on_vector(i_point, o2, o2 + v2, precision))
                        return true;
                }
            }
            
            // try vectors with v2 starting from v2_end
            // t1 = (v2B*(o1A - (o2A + v2A)) + v2A*(o2B + v2B - o1B) / (v1A*v2B - v1B*v2A)
            // t2 = (o2B + v2B - o1B - t1*v1B) / v2B
            t1 = v2B * (o1A - (o2A + v2A)) + v2A * (o2B + v2B - o1B);
            t1 /= bottom;
#ifdef DEBUG_INTERSECT_VECTORS
            cout << "intersect_vectors_dl_t2 t1: " << t1 << "\n";
#endif
            if (t1 >= 0 && t1 <= 1)
            {
                Vector_3D::Measurement t2((o2B + v2B - o1B - t1 * v1B) / v2B);
#ifdef DEBUG_INTERSECT_VECTORS
                cout << "intersect_vectors_dl_t2 t2: " << t2 << "\n";
#endif
                if (t2 >= 0 && t2 <= 1)
                {
                    i_point = o1 + v1 - v1 * t1;
#ifdef DEBUG_INTERSECT_VECTORS
                    cout << "intersect_vectors_dl_t2 i_point_v2_t1 x: " << i_point.get_x() << " y: " << i_point.get_y() << " z: " << i_point.get_z() << "\n";
                    cout << "intersect_vectors_dl_t2 is_pt_on_vector1: " << is_pt_on_vector(i_point, o1, o1 + v1, precision) << "\n";
                    cout << "intersect_vectors_dl_t2 is_pt_on_vector2: " << is_pt_on_vector(i_point, o2, o2 + v2, precision) << "\n";
#endif
                    if (is_pt_on_vector(i_point, o1, o1 + v1, precision) && is_pt_on_vector(i_point, o2, o2 + v2, precision))
                        return true;
                    i_point = o2 + v2 - v2 * t2;
#ifdef DEBUG_INTERSECT_VECTORS
                    cout << "intersect_vectors_dl_t2 i_point_v2_t2 x: " << i_point.get_x() << " y: " << i_point.get_y() << " z: " << i_point.get_z() << "\n";
                    cout << "intersect_vectors_dl_t2 is_pt_on_vector1: " << is_pt_on_vector(i_point, o1, o1 + v1, precision) << "\n";
                    cout << "intersect_vectors_dl_t2 is_pt_on_vector2: " << is_pt_on_vector(i_point, o2, o2 + v2, precision) << "\n";
#endif
                    if (is_pt_on_vector(i_point, o1, o1 + v1, precision) && is_pt_on_vector(i_point, o2, o2 + v2, precision))
                        return true;
                }
            }
        }
#ifdef DEBUG_INTERSECT_VECTORS
        else
            cout << "intersect_vectors_dl_t2 bottom is zero\n";
        cout << "intersect_vectors_dl_t2 returning false\n";
#endif        
        return false;
    }

    const bool intersect_vectors_dl_cpz(const Point_3D& v1_start, const Point_3D& v1_end, 
            const Point_3D& v2_start, const Point_3D& v2_end, Point_3D& i_point, 
            const Vector_3D::Measurement precision)
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
        // pZ = o1Z + t*v1Z
        // where t is from 0 to 1
        //
        // cross product will be of V3 and v2 which should equal zero
        // v3 is defined as
        // v3 = p - o2
        // three equations: x, y, and z
        // each must be zero
        // cpX = (v2Y * v3Z) - (v2Z * v3Y) = 0
        // cpY = (v2X * v3Z) - (v2Z * v3X) = 0
        // cpZ = (v2X * v3Y) - (v2Y * v3X) = 0
        //
        // cpx
        // (v2Y * v3Z) - (v2Z * v3Y) = 0
        // (V2Y * (o1Z + t*v1Z - o2Z)) - (v2Z * (o1Y + t*v1Y - o2Y)) = 0
        // v2Y*o1Z + t*v1Z*v2Y - o2Z*v2Y - v2Z*o1Y - t*v1Y*v2Z + o2Y*v2Z = 0
        // t(v1Z*v2Y - v1Y*v2Z) = o2Z*v2Y + v2Z*o1Y - v2Y*o1Z - o2Y*v2Z
        // t = (o2Z*v2Y + v2Z*o1Y - v2Y*o1Z - o2Y*v2Z) / (v1Z*v2Y - v1Y*v2Z)
        // t = (v2Y*(o2Z - o1Z) + v2Z*(o1Y - o2Y)) / (v1Z*v2Y - v1Y*v2Z)
        //
        // cpy
        // (v2X * v3Z) - (v2Z * v3X) = 0
        // (v2X * (o1Z + t*v1Z - o2Z)) - (v2Z * (o1X + t*v1X - o2X)) = 0
        // v2X*o1Z + t*v1Z*v2X - o2Z*v2X - v2Z*o1X - t*v1X*v2Z + o2X*v2Z = 0
        // t(v1Z*v2X - v1X*v2Z) = o2Z*v2X + v2Z*o1X - v2X*o1Z - o2X*v2Z
        // t = (o2Z*v2X + v2Z*o1X - v2X*o1Z - o2X*v2Z) / (v1Z*v2X - v1X*v2Z)
        // t = (v2X*(o2Z - o1Z) + v2Z*(o1X - o2X)) / (v1Z*v2X - v1X*v2Z)
        //
        // cpz
        // (v2X * v3Y) - (v2Y * v3X) = 0
        // (v2X * (o1Y + t*v1Y - o2Y)) - (v2Y * (o1X + t*v1X - o2X)) = 0
        // v2X*o1Y + t*v1Y*v2X - o2Y*v2X - v2Y*o1X - t*v1X*v2Y + o2X*v2Y = 0
        // t(v1Y*v2X - v1X*v2Y) = o2Y*v2X + v2Y*o1X - v2X*o1Y - o2X*v2Y
        // t = (o2Y*v2X + v2Y*o1X - v2X*o1Y - o2X*v2Y) / (v1Y*v2X - v1X*v2Y)
        // t = (v2X*(o2Y - o1Y) + v2Y*(o1X - o2X)) / (v1Y*v2X - v1X*v2Y)
        //
        // try first equation
        // check for zero division
        // t = (v2Y*(o2Z - o1Z) + v2Z*(o1Y - o2Y)) / (v1Z*v2Y - v1Y*v2Z)
        Vector_3D v1(v1_start, v1_end);
        Vector_3D v2(v2_start, v2_end);
#ifdef DEBUG_INTERSECT_VECTORS
        cout << "intersect_vectors_dl_cpz v1 x: " << v1.get_x() << " y: " << v1.get_y() << " z: " << v1.get_z() << "\n";
        cout << "intersect_vectors_dl_cpz v2 x: " << v2.get_x() << " y: " << v2.get_y() << " z: " << v2.get_z() << "\n";
#endif
        Vector_3D::Measurement bottom(v1.get_z() * v2.get_y());
        Vector_3D::Measurement largest(fabs(bottom));
        Vector_3D::Measurement next(v1.get_y() * v2.get_z());
        bottom -= next;
        Vector_3D::Measurement error_bound(fabs(next));
        if (error_bound > largest)
            largest = error_bound;
        error_bound = largest > 1.0 ? (largest * precision) : precision;
//        cout << "intersect_vectors_dl_cpz fabs(bottom): " << fabs(bottom) << " error_bound: " << error_bound << "\n";
        if (fabs(bottom) > error_bound)
        {
            Vector_3D::Measurement t(v2.get_y() * (v2_start.get_z() - v1_start.get_z()) + v2.get_z() * (v1_start.get_y() - v2_start.get_y()));
            t /= bottom;
#ifdef DEBUG_INTERSECT_VECTORS
            cout << "intersect_vectors_dl_cpz first equation bottom: " << bottom << " t: " << t << "\n";
#endif
            if (t >= 0 && t <= 1)
            {
                // point is on v1, now check if it is in v2
                i_point = v1_start + v1 * t;
#ifdef DEBUG_INTERSECT_VECTORS
                cout << "intersect_vectors_dl_cpz i_point_cpz_x x: " << i_point.get_x() << " y: " << i_point.get_y() << " z: " << i_point.get_z() << "\n";
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
        
        // try second equation
        // t = (v2X*(o2Z - o1Z) + v2Z*(o1X - o2X)) / (v1Z*v2X - v1X*v2Z)
        bottom = v1.get_z() * v2.get_x();
        largest = fabs(bottom);
        next = v1.get_x() * v2.get_z();
        bottom -= next;
        error_bound = fabs(next);
        if (error_bound > largest)
            largest = error_bound;
        error_bound = largest > 1.0 ? (largest * precision) : precision;
//        cout << "intersect_vectors_dl_cpz fabs(bottom): " << fabs(bottom) << " error_bound: " << error_bound << "\n";
        // check for zero division
        if (fabs(bottom) > error_bound)
        {
            Vector_3D::Measurement t(v2.get_x() * (v2_start.get_z() - v1_start.get_z()) + v2.get_z() * (v1_start.get_x() - v2_start.get_x()));
            t /= bottom;
#ifdef DEBUG_INTERSECT_VECTORS
            cout << "intersect_vectors_dl_cpz second equation bottom: " << bottom << " t: " << t << "\n";
#endif
            if (t >= 0 && t <= 1)
            {
                // point is on v1, now check if it is in v2
                i_point = v1_start + v1 * t;
#ifdef DEBUG_INTERSECT_VECTORS
                cout << "intersect_vectors_dl_cpz i_point_cpz_y x: " << i_point.get_x() << " y: " << i_point.get_y() << " z: " << i_point.get_z() << "\n";
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
        
        // try third equation
        // t = (v2X*(o2Y - o1Y) + v2Y*(o1X - o2X)) / (v1Y*v2X - v1X*v2Y)
        // check for zero division
        bottom = v1.get_y() * v2.get_x();
        largest = fabs(bottom);
        next = v1.get_x() * v2.get_y();
        bottom -= next;
        error_bound = fabs(next);
        if (error_bound > largest)
            largest = error_bound;
        error_bound = largest > 1.0 ? (largest * precision) : precision;
//        cout << "intersect_vectors_dl_cpz fabs(bottom): " << fabs(bottom) << " error_bound: " << error_bound << "\n";
        if (fabs(bottom) > error_bound)
        {
            Vector_3D::Measurement t(v2.get_x() * (v2_start.get_y() - v1_start.get_y()) + v2.get_y() * (v1_start.get_x() - v2_start.get_x()));
            t /= bottom;
#ifdef DEBUG_INTERSECT_VECTORS
            cout << "intersect_vectors_dl_cpz third equation bottom: " << bottom << " t: " << t << "\n";
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
    
    const bool intersect_vectors_dl(const Point_3D& v1_start, 
            const Point_3D& v1_end, const Point_3D& v2_start, 
            const Point_3D& v2_end, Vector_3D_idata& result, 
            const Vector_3D::Measurement precision)
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
        
#ifndef INTERSECT_VECTORS_USE_LINE
#ifdef DEBUG_INTERSECT_VECTORS
        cout << "intersect_vectors_dl Testing cross product zero method\n";
#endif
        // try cross product equal to zero method
        Point_3D i_point(0,0,0);
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
#else
        Vector_3D::Measurement largest(fabs(v1_end.get_x()));
        Vector_3D::Measurement error_bound(fabs(v1_start.get_x()));
        if (error_bound > largest)
            largest = error_bound;
        
        error_bound = largest > 1.0 ? (largest * precision) : precision;
        if (fabs(v1.get_x()) > error_bound) // v1 x is not zero
        {
#ifdef DEBUG_INTERSECT_VECTORS
            cout << "intersect_vectors_dl testing v1 X Y\n";
#endif
            Point_3D i_point(0,0,0);
            if (intersect_vectors_dl_t1(v1_start, v1, v2_start, v2, X, Y, i_point, precision))
            {
                result.p1 = i_point;
                result.num = 1;
                return true;
            }
#ifdef DEBUG_INTERSECT_VECTORS
            cout << "intersect_vectors_dl testing v1 X Z\n";
#endif
            if (intersect_vectors_dl_t1(v1_start, v1, v2_start, v2, X, Z, i_point, precision))
            {
                result.p1 = i_point;
                result.num = 1;
                return true;
            }
        }
        
        largest = fabs(v1_end.get_y());
        error_bound = fabs(v1_start.get_y());
        if (error_bound > largest)
            largest = error_bound;
        
        error_bound = largest > 1.0 ? (largest * precision) : precision;
        
        if (fabs(v1.get_y()) > error_bound) // v1 y is not zero
        {
#ifdef DEBUG_INTERSECT_VECTORS
            cout << "intersect_vectors_dl testing v1 Y Z\n";
#endif
            Point_3D i_point(0,0,0);
            if (intersect_vectors_dl_t1(v1_start, v1, v2_start, v2, Y, Z, i_point, precision))
            {
                result.p1 = i_point;
                result.num = 1;
                return true;
            }
#ifdef DEBUG_INTERSECT_VECTORS
            cout << "intersect_vectors_dl testing v1 Y X\n";
#endif
            if (intersect_vectors_dl_t1(v1_start, v1, v2_start, v2, Y, X, i_point, precision))
            {
                result.p1 = i_point;
                result.num = 1;
                return true;
            }
        }
        
        largest = fabs(v1_end.get_z());
        error_bound = fabs(v1_start.get_z());
        if (error_bound > largest)
            largest = error_bound;
        
        error_bound = largest > 1.0 ? (largest * precision) : precision;
        if (fabs(v1.get_z()) > error_bound) // v1 z is not zero
        {
#ifdef DEBUG_INTERSECT_VECTORS
            cout << "intersect_vectors_dl testing v1 Z X\n";
#endif
            Point_3D i_point(0,0,0);
            if (intersect_vectors_dl_t1(v1_start, v1, v2_start, v2, Z, X, i_point, precision))
            {
                result.p1 = i_point;
                result.num = 1;
                return true;
            }
#ifdef DEBUG_INTERSECT_VECTORS
            cout << "intersect_vectors_dl testing v1 Z Y\n";
#endif
            if (intersect_vectors_dl_t1(v1_start, v1, v2_start, v2, Z, Y, i_point, precision))
            {
                result.p1 = i_point;
                result.num = 1;
                return true;
            }
        }
        
        largest = fabs(v2_end.get_x());
        error_bound = fabs(v2_start.get_x());
        if (error_bound > largest)
            largest = error_bound;
        
        error_bound = largest > 1.0 ? (largest * precision) : precision;
        if (fabs(v2.get_x()) > error_bound) // v2 x is not zero
        {
#ifdef DEBUG_INTERSECT_VECTORS
            cout << "intersect_vectors_dl testing v2 X Y\n";
#endif
            Point_3D i_point(0,0,0);
            if (intersect_vectors_dl_t2(v1_start, v1, v2_start, v2, X, Y, i_point, precision))
            {
                result.p1 = i_point;
                result.num = 1;
                return true;
            }
#ifdef DEBUG_INTERSECT_VECTORS
            cout << "intersect_vectors_dl testing v2 X Z\n";
#endif
            if (intersect_vectors_dl_t2(v1_start, v1, v2_start, v2, X, Z, i_point, precision))
            {
                result.p1 = i_point;
                result.num = 1;
                return true;
            }
        }
        
        largest = fabs(v2_end.get_y());
        error_bound = fabs(v2_start.get_y());
        if (error_bound > largest)
            largest = error_bound;
        
        error_bound = largest > 1.0 ? (largest * precision) : precision;
        if (fabs(v2.get_y()) > error_bound) // v2 y is not zero
        {
#ifdef DEBUG_INTERSECT_VECTORS
            cout << "intersect_vectors_dl testing v2 Y Z\n";
#endif
            Point_3D i_point(0,0,0);
            if (intersect_vectors_dl_t2(v1_start, v1, v2_start, v2, Y, Z, i_point, precision))
            {
                result.p1 = i_point;
                result.num = 1;
                return true;
            }
#ifdef DEBUG_INTERSECT_VECTORS
            cout << "intersect_vectors_dl testing v2 Y X\n";
#endif
            if (intersect_vectors_dl_t2(v1_start, v1, v2_start, v2, Y, X, i_point, precision))
            {
                result.p1 = i_point;
                result.num = 1;
                return true;
            }
        }
        
        largest = fabs(v2_end.get_z());
        error_bound = fabs(v2_start.get_z());
        if (error_bound > largest)
            largest = error_bound;
        
        error_bound = largest > 1.0 ? (largest * precision) : precision;
        if (fabs(v2.get_z()) > error_bound) // v2 z is not zero
        {
#ifdef DEBUG_INTERSECT_VECTORS
            cout << "intersect_vectors_dl testing v2 Z X\n";
#endif
            Point_3D i_point(0,0,0);
            if (intersect_vectors_dl_t2(v1_start, v1, v2_start, v2, Z, X, i_point, precision))
            {
                result.p1 = i_point;
                result.num = 1;
                return true;
            }
#ifdef DEBUG_INTERSECT_VECTORS
            cout << "intersect_vectors_dl testing v2 Z Y\n";
#endif
            if (intersect_vectors_dl_t2(v1_start, v1, v2_start, v2, Z, Y, i_point, precision))
            {
                result.p1 = i_point;
                result.num = 1;
                return true;
            }
        }
#endif
        
#ifdef DEBUG_INTERSECT_VECTORS
        cout << "intersect_vectors_dl returning false\n";
#endif
        return false;
    }
    
    const bool intersect_vectors(
            const Point_3D& v1_start, const Point_3D& v1_end,
            const Point_3D& v2_start, const Point_3D& v2_end, 
            Vector_3D_idata& result, const Vector_3D::Measurement precision)
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
    }
}
