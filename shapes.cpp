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
 * File:   shapes.cpp
 * Author: Jeffrey Davis
 */

#include "shapes.h"
#include <queue>
#include <cmath>
#include <memory>
#include "Point_2D.h"
#include "Point_3D.h"
#include "Facet_2D.h"
#include "Facet_3D.h"
#include "Mesh_2D.h"
#include "Mesh_2D.h"
#include "Mesh_3D.h"
#include "Mesh_3D.h"
#include <iostream>

namespace VCAD_lib
{
    void m_rectangle(Mesh_2D& mesh, const Point_2D::Measurement x_length, 
            const Point_2D::Measurement y_length, const bool center, 
            const Point_2D& origin)
    {
        Point_2D::Measurement orig_x = center ? origin.get_x() - (x_length / 2) : origin.get_x();
        Point_2D::Measurement orig_y = center ? origin.get_y() - (y_length / 2) : origin.get_y();

        shared_ptr<Point_2D> p1(new Point_2D(orig_x, orig_y));
        shared_ptr<Point_2D> p2(new Point_2D(orig_x, orig_y + y_length));
        shared_ptr<Point_2D> p3(new Point_2D(orig_x + x_length, orig_y));
        shared_ptr<Point_2D> p4(new Point_2D(orig_x + x_length, orig_y + y_length));
        
        mesh.push_back(Facet_2D(p1, p4, p2));
        mesh.push_back(Facet_2D(p1, p3, p4));
    }
    
    void m_circle(Mesh_2D& mesh, const Point_2D::Measurement radius, 
            const int steps_per_quarter, const Point_2D& origin)
    {
        const Point_2D::Angle_Meas two_pi = 6.283185307179586;
        const Point_2D::Angle_Meas step = two_pi / (4 * steps_per_quarter);
        // generate 2 dimensional circle
        Point_2D::Angle_Meas angle = step;
        shared_ptr<Point_2D> begin(new Point_2D(origin.get_x() + radius, origin.get_y()));
        shared_ptr<Point_2D> middle(new Point_2D(polar_point(radius, angle, origin)));
        angle += step;
        while (angle < two_pi)
        {
            shared_ptr<Point_2D> last(new Point_2D(polar_point(radius, angle, origin)));
            
            mesh.push_back(Facet_2D(begin, middle, last));
            middle = last;
            angle += step;
        }
    }
    
    void m_ellipse(Mesh_2D& mesh, const Point_2D::Measurement x_radius, 
            const Point_2D::Measurement y_radius, const int steps_per_quarter, 
            const Point_2D& origin)
    {
        // Ax^2 + By^2 = C
        // first equation = x_radius = sqrt(C / A)
        // second equation = y_radius = sqrt(C / B)
        // choose C to be A*B
        // x_radius = sqrt(B)
        // y_radius = sqrt(A)
        // A = y_radius^2
        // B = x_radius^2
        // C = A*B
        Point_2D::Measurement y_mult = pow(x_radius,2);
        Point_2D::Measurement x_mult = pow(y_radius,2);
        Point_2D::Measurement result = x_mult * y_mult;
        // Ax^2 + By^2 = C
        // y^2 = ( C - Ax^2 ) / B
        // y = +/- sqrt((C - Ax^2) / B)
        
        // find x range, set y = 0
        // Ax^2 = C
        // x = +/- sqrt(C / A)
        const Point_2D::Measurement x_max = sqrt(result / x_mult);
        const Point_2D::Measurement step = x_max / steps_per_quarter;
        const Point_2D::Measurement orig_x = origin.get_x();
        const Point_2D::Measurement orig_y = origin.get_y();
        Point_2D::Measurement x_pos = -x_max;
        shared_ptr<Point_2D> begin(new Point_2D(orig_x + x_pos,orig_y));
        x_pos += step;
        Point_2D::Measurement y_val = sqrt((result - x_mult * pow(x_pos, 2)) / y_mult);
        shared_ptr<Point_2D> t_middle(new Point_2D(orig_x + x_pos,orig_y + y_val));
        shared_ptr<Point_2D> b_middle(new Point_2D(orig_x + x_pos,orig_y - y_val));
        x_pos += step;
        while (x_pos < x_max)
        {
            y_val = sqrt((result - x_mult * pow(x_pos, 2)) / y_mult);
            shared_ptr<Point_2D> t_last(new Point_2D(orig_x + x_pos,orig_y + y_val));
            shared_ptr<Point_2D> b_last(new Point_2D(orig_x + x_pos,orig_y - y_val));
            // create facets
            mesh.push_back(Facet_2D(begin,t_last,t_middle));
            mesh.push_back(Facet_2D(begin,b_middle,b_last));
            // move last to middle
            t_middle = t_last;
            b_middle = b_last;
            // update x_pos
            x_pos += step;
        }
        
        // create last facets
        shared_ptr<Point_2D> end(new Point_2D(orig_x + x_max,orig_y));
        mesh.push_back(Facet_2D(begin,end,t_middle));
        mesh.push_back(Facet_2D(begin,b_middle,end));
    }
    
    void m_cuboid(Mesh_3D& mesh, const Point_3D::Measurement x_length, const Point_3D::Measurement y_length, 
            const Point_3D::Measurement z_length, const bool center, const Point_3D& origin)
    {
        const Point_3D::Measurement orig_x = center ? origin.get_x() - (x_length / 2) : origin.get_x();
        const Point_3D::Measurement orig_y = center ? origin.get_y() - (y_length / 2) : origin.get_y();
        const Point_3D::Measurement orig_z = center ? origin.get_z() - (z_length / 2) : origin.get_z();

        shared_ptr<Point_3D> p1(new Point_3D(orig_x,orig_y,orig_z));
        shared_ptr<Point_3D> p2(new Point_3D(orig_x,orig_y + y_length,orig_z));
        shared_ptr<Point_3D> p3(new Point_3D(orig_x + x_length,orig_y + y_length,orig_z));
        shared_ptr<Point_3D> p4(new Point_3D(orig_x + x_length,orig_y,orig_z));
        shared_ptr<Point_3D> p5(new Point_3D(orig_x,orig_y,orig_z + z_length));
        shared_ptr<Point_3D> p6(new Point_3D(orig_x,orig_y + y_length,orig_z + z_length));
        shared_ptr<Point_3D> p7(new Point_3D(orig_x + x_length,orig_y + y_length,orig_z + z_length));
        shared_ptr<Point_3D> p8(new Point_3D(orig_x + x_length,orig_y,orig_z + z_length));
        
        // bottom
        mesh.push_back(Facet_3D(p1, p2, p3));
        mesh.push_back(Facet_3D(p1, p3, p4));
        // top
        mesh.push_back(Facet_3D(p5, p7, p6));
        mesh.push_back(Facet_3D(p5, p8, p7));
        // sides
        mesh.push_back(Facet_3D(p1, p5, p6));
        mesh.push_back(Facet_3D(p1, p6, p2));
        mesh.push_back(Facet_3D(p4, p7, p8));
        mesh.push_back(Facet_3D(p4, p3, p7));
        mesh.push_back(Facet_3D(p1, p8, p5));
        mesh.push_back(Facet_3D(p1, p4, p8));
        mesh.push_back(Facet_3D(p2, p6, p7));
        mesh.push_back(Facet_3D(p2, p7, p3));
    }
    
    void m_cylinder(Mesh_3D& mesh, const Point_3D::Measurement b_radius, 
            const Point_3D::Measurement t_radius, const Point_3D::Measurement height, 
            const int steps_per_quarter, const bool center, const Point_3D& origin)
    {
        const Point_3D::Angle_Meas two_pi = 6.283185307179586;
        const Point_3D::Angle_Meas step = two_pi / (4 * steps_per_quarter);
        const Point_3D::Measurement orig_x = origin.get_x();
        const Point_3D::Measurement orig_y = origin.get_y();
        const Point_3D::Measurement orig_z = origin.get_z();
        const Point_3D::Measurement b_height = center ? -height / 2 : 0;
        const Point_3D::Measurement t_height = center ? height / 2 : height;
        Point_3D::Angle_Meas angle = step;
        // generate bottom circle
        if (fabs(b_radius) < mesh.get_precision()) // bottom single point
        {
            // single point
            shared_ptr<Point_3D> bottom(new Point_3D(orig_x,orig_y,orig_z + b_height));

            // Make top and side facets together
            shared_ptr<Point_3D> begin(new Point_3D(orig_x + t_radius, orig_y, orig_z + t_height));
            shared_ptr<Point_3D> middle(new Point_3D(cylindrical_point(t_radius, angle, t_height, origin)));
            mesh.push_back(Facet_3D(bottom, middle, begin)); // add side facet
            angle += step;
            while (angle < two_pi)
            {
                shared_ptr<Point_3D> last(new Point_3D(cylindrical_point(t_radius, angle, t_height, origin)));

                // create top facet
                mesh.push_back(Facet_3D(begin, middle, last)); 
                // create side facet
                mesh.push_back(Facet_3D(bottom, last, middle));

                middle = last;
                angle += step;
            }
            mesh.push_back(Facet_3D(bottom, begin, middle)); // add last side facet
        }
        else if (fabs(t_radius) < mesh.get_precision()) // top single point
        {
            // single point
            shared_ptr<Point_3D> top(new Point_3D(orig_x,orig_y,orig_z + t_height));

            // Make top and side facets together
            shared_ptr<Point_3D> begin(new Point_3D(orig_x + b_radius, orig_y, orig_z + b_height));
            shared_ptr<Point_3D> middle(new Point_3D(cylindrical_point(b_radius, angle, b_height, origin)));
            mesh.push_back(Facet_3D(top, begin, middle)); // add side facet
            angle += step;
            while (angle < two_pi)
            {
                shared_ptr<Point_3D> last(new Point_3D(cylindrical_point(b_radius, angle, b_height, origin)));

                // create bottom facet
                mesh.push_back(Facet_3D(begin, last, middle)); 
                // create side facet
                mesh.push_back(Facet_3D(top, middle, last));

                middle = last;
                angle += step;
            }
            mesh.push_back(Facet_3D(top, middle, begin)); // add last side facet
        }
        else
        {
            shared_ptr<Point_3D> b_begin(new Point_3D(orig_x + b_radius, orig_y, orig_z + b_height));
            shared_ptr<Point_3D> t_begin(new Point_3D(orig_x + t_radius, orig_y, orig_z + t_height));
            shared_ptr<Point_3D> b_middle(new Point_3D(cylindrical_point(b_radius, angle, b_height, origin)));
            shared_ptr<Point_3D> t_middle(new Point_3D(cylindrical_point(t_radius, angle, t_height, origin)));
            mesh.push_back(Facet_3D(b_begin, t_middle, t_begin));  // add side facet
            mesh.push_back(Facet_3D(b_begin, b_middle, t_middle)); // add side facet
            angle += step;
            while (angle < two_pi)
            {
                shared_ptr<Point_3D> b_last(new Point_3D(cylindrical_point(b_radius, angle, b_height, origin)));
                shared_ptr<Point_3D> t_last(new Point_3D(cylindrical_point(t_radius, angle, t_height, origin)));

                mesh.push_back(Facet_3D(b_begin, b_last, b_middle));  // add bottom facet
                mesh.push_back(Facet_3D(t_begin, t_middle, t_last));  // add top facet
                mesh.push_back(Facet_3D(b_middle, t_last, t_middle)); // add side facet
                mesh.push_back(Facet_3D(b_middle, b_last, t_last));   // add side facet
                b_middle = b_last;
                t_middle = t_last;
                angle += step;
            }
            mesh.push_back(Facet_3D(b_begin, t_begin, t_middle));  // add last side facets
            mesh.push_back(Facet_3D(b_begin, t_middle, b_middle)); // add last side facets
        }
    }
    
    void m_e_cylinder(Mesh_3D& mesh, const Point_3D::Measurement b_x_radius, 
            const Point_3D::Measurement b_y_radius, const Point_3D::Measurement t_x_radius, 
            const Point_3D::Measurement t_y_radius, const Point_3D::Measurement height, 
            const int steps_per_quarter, const bool center, const Point_3D& origin)
    {
        // bottom ellipse parameters
        Point_3D::Measurement b_y_mult = pow(b_x_radius,2);
        Point_3D::Measurement b_x_mult = pow(b_y_radius,2);
        Point_3D::Measurement b_result = b_x_mult * b_y_mult;

        // top ellipse parameters
        Point_3D::Measurement t_y_mult = pow(t_x_radius,2);
        Point_3D::Measurement t_x_mult = pow(t_y_radius,2);
        Point_3D::Measurement t_result = t_x_mult * t_y_mult;
        
        // cylinder parameters
        const Point_3D::Measurement orig_x = origin.get_x();
        const Point_3D::Measurement orig_y = origin.get_y();
        const Point_3D::Measurement orig_z = origin.get_z();
        const Point_3D::Measurement b_height = center ? -height / 2 : 0;
        const Point_3D::Measurement t_height = center ? height / 2 : height;
        // generate bottom circle
        if (fabs(b_x_radius) < mesh.get_precision() && 
                fabs(b_y_radius) < mesh.get_precision()) // bottom single point
        {
            // single point
            shared_ptr<Point_3D> bottom(new Point_3D(orig_x,orig_y,orig_z + b_height));

            // Make top and side facets together
            const Point_3D::Measurement x_max = sqrt(t_result / t_x_mult);
            const Point_3D::Measurement step = x_max / steps_per_quarter;
            Point_3D::Measurement x_pos = -x_max;
            shared_ptr<Point_3D> begin(new Point_3D(orig_x + x_pos,orig_y,orig_z + t_height));
            x_pos += step;
            Point_3D::Measurement y_val = sqrt((t_result - t_x_mult * pow(x_pos, 2)) / t_y_mult);
            shared_ptr<Point_3D> t_middle(new Point_3D(orig_x + x_pos,orig_y + y_val, orig_z + t_height));
            shared_ptr<Point_3D> b_middle(new Point_3D(orig_x + x_pos,orig_y - y_val, orig_z + t_height));
            // create side facets
            mesh.push_back(Facet_3D(bottom,begin,t_middle));
            mesh.push_back(Facet_3D(bottom,b_middle,begin));
            x_pos += step;
            while (x_pos < x_max)
            {
                y_val = sqrt((t_result - t_x_mult * pow(x_pos, 2)) / t_y_mult);
                shared_ptr<Point_3D> t_last(new Point_3D(orig_x + x_pos,orig_y + y_val, orig_z + t_height));
                shared_ptr<Point_3D> b_last(new Point_3D(orig_x + x_pos,orig_y - y_val, orig_z + t_height));
                // create top facets
                mesh.push_back(Facet_3D(begin,t_last,t_middle));
                mesh.push_back(Facet_3D(begin,b_middle,b_last));
                // create side facets
                mesh.push_back(Facet_3D(bottom,t_middle,t_last));
                mesh.push_back(Facet_3D(bottom,b_last,b_middle));
                // move last to middle
                t_middle = t_last;
                b_middle = b_last;
                // update x_pos
                x_pos += step;
            }

            // create last facets
            shared_ptr<Point_3D> end(new Point_3D(orig_x + x_max,orig_y,orig_z + t_height));
            mesh.push_back(Facet_3D(begin,end,t_middle));
            mesh.push_back(Facet_3D(begin,b_middle,end));
            mesh.push_back(Facet_3D(bottom,t_middle,end));
            mesh.push_back(Facet_3D(bottom,end,b_middle));
        }
        else if (fabs(t_x_radius) < mesh.get_precision() && 
                fabs(t_y_radius) < mesh.get_precision()) // top single point
        {
            // single point
            shared_ptr<Point_3D> top(new Point_3D(orig_x,orig_y,orig_z + t_height));

            // Make top and side facets together
            const Point_3D::Measurement x_max = sqrt(b_result / b_x_mult);
            const Point_3D::Measurement step = x_max / steps_per_quarter;
            Point_3D::Measurement x_pos = -x_max;
            shared_ptr<Point_3D> begin(new Point_3D(orig_x + x_pos,orig_y,orig_z + b_height));
            x_pos += step;
            Point_3D::Measurement y_val = sqrt((b_result - b_x_mult * pow(x_pos, 2)) / b_y_mult);
            shared_ptr<Point_3D> t_middle(new Point_3D(orig_x + x_pos,orig_y + y_val, orig_z + b_height));
            shared_ptr<Point_3D> b_middle(new Point_3D(orig_x + x_pos,orig_y - y_val, orig_z + b_height));
            // create side facets
            mesh.push_back(Facet_3D(top,t_middle,begin));
            mesh.push_back(Facet_3D(top,begin,b_middle));
            x_pos += step;
            while (x_pos < x_max)
            {
                y_val = sqrt((b_result - b_x_mult * pow(x_pos, 2)) / b_y_mult);
                shared_ptr<Point_3D> t_last(new Point_3D(orig_x + x_pos,orig_y + y_val, orig_z + b_height));
                shared_ptr<Point_3D> b_last(new Point_3D(orig_x + x_pos,orig_y - y_val, orig_z + b_height));
                // create bottom facets
                mesh.push_back(Facet_3D(begin,t_middle,t_last));
                mesh.push_back(Facet_3D(begin,b_last,b_middle));
                // create side facets
                mesh.push_back(Facet_3D(top,t_last,t_middle));
                mesh.push_back(Facet_3D(top,b_middle,b_last));
                // move last to middle
                t_middle = t_last;
                b_middle = b_last;
                // update x_pos
                x_pos += step;
            }

            // create last facets
            shared_ptr<Point_3D> end(new Point_3D(orig_x + x_max,orig_y,orig_z + b_height));
            mesh.push_back(Facet_3D(begin,t_middle,end));
            mesh.push_back(Facet_3D(begin,end,b_middle));
            mesh.push_back(Facet_3D(top,end,t_middle));
            mesh.push_back(Facet_3D(top,b_middle,end));
        }
        else
        {
            const Point_3D::Measurement b_x_max = sqrt(b_result / b_x_mult);
            const Point_3D::Measurement t_x_max = sqrt(t_result / t_x_mult);
            const Point_3D::Measurement b_step = b_x_max / steps_per_quarter;
            const Point_3D::Measurement t_step = t_x_max / steps_per_quarter;
            Point_3D::Measurement b_x_pos = -b_x_max;
            Point_3D::Measurement t_x_pos = -t_x_max;
            shared_ptr<Point_3D> b_begin(new Point_3D(orig_x + b_x_pos,orig_y,orig_z + b_height));
            shared_ptr<Point_3D> t_begin(new Point_3D(orig_x + t_x_pos,orig_y,orig_z + t_height));
            b_x_pos += b_step;
            t_x_pos += t_step;
            Point_3D::Measurement b_y_val = sqrt((b_result - b_x_mult * pow(b_x_pos, 2)) / b_y_mult);
            Point_3D::Measurement t_y_val = sqrt((t_result - t_x_mult * pow(t_x_pos, 2)) / t_y_mult);
            // TODO
            shared_ptr<Point_3D> b_p_middle(new Point_3D(orig_x + b_x_pos,orig_y + b_y_val, orig_z + b_height));
            shared_ptr<Point_3D> b_n_middle(new Point_3D(orig_x + b_x_pos,orig_y - b_y_val, orig_z + b_height));
            shared_ptr<Point_3D> t_p_middle(new Point_3D(orig_x + t_x_pos,orig_y + t_y_val, orig_z + t_height));
            shared_ptr<Point_3D> t_n_middle(new Point_3D(orig_x + t_x_pos,orig_y - t_y_val, orig_z + t_height));
            // create side facets
            mesh.push_back(Facet_3D(b_begin,t_begin,t_p_middle));
            mesh.push_back(Facet_3D(b_begin,t_p_middle,b_p_middle));
            mesh.push_back(Facet_3D(b_begin,t_n_middle,t_begin));
            mesh.push_back(Facet_3D(b_begin,b_n_middle,t_n_middle));
            b_x_pos += b_step;
            t_x_pos += t_step;
            while (b_x_pos < b_x_max)
            {
                b_y_val = sqrt((b_result - b_x_mult * pow(b_x_pos, 2)) / b_y_mult);
                t_y_val = sqrt((t_result - t_x_mult * pow(t_x_pos, 2)) / t_y_mult);
                shared_ptr<Point_3D> b_p_last(new Point_3D(orig_x + b_x_pos,orig_y + b_y_val, orig_z + b_height));
                shared_ptr<Point_3D> b_n_last(new Point_3D(orig_x + b_x_pos,orig_y - b_y_val, orig_z + b_height));
                shared_ptr<Point_3D> t_p_last(new Point_3D(orig_x + t_x_pos,orig_y + t_y_val, orig_z + t_height));
                shared_ptr<Point_3D> t_n_last(new Point_3D(orig_x + t_x_pos,orig_y - t_y_val, orig_z + t_height));
                // create bottom facets
                mesh.push_back(Facet_3D(b_begin,b_p_middle,b_p_last));
                mesh.push_back(Facet_3D(b_begin,b_n_last,b_n_middle));
                // create top facets
                mesh.push_back(Facet_3D(t_begin,t_p_last,t_p_middle));
                mesh.push_back(Facet_3D(t_begin,t_n_middle,t_n_last));
                // create side facets
                mesh.push_back(Facet_3D(b_p_middle,t_p_middle,t_p_last));
                mesh.push_back(Facet_3D(b_p_middle,t_p_last,b_p_last));
                mesh.push_back(Facet_3D(b_n_middle,t_n_last,t_n_middle));
                mesh.push_back(Facet_3D(b_n_middle,b_n_last,t_n_last));
                // move last to middle
                b_p_middle = b_p_last;
                b_n_middle = b_n_last;
                t_p_middle = t_p_last;
                t_n_middle = t_n_last;
                // update x_pos
                b_x_pos += b_step;
                t_x_pos += t_step;
            }

            // create last facets
            shared_ptr<Point_3D> b_end(new Point_3D(orig_x + b_x_max,orig_y,orig_z + b_height));
            shared_ptr<Point_3D> t_end(new Point_3D(orig_x + t_x_max,orig_y,orig_z + t_height));
            mesh.push_back(Facet_3D(b_begin,b_p_middle,b_end));
            mesh.push_back(Facet_3D(b_begin,b_end,b_n_middle));
            mesh.push_back(Facet_3D(t_begin,t_end,t_p_middle));
            mesh.push_back(Facet_3D(t_begin,t_n_middle,t_end));
            mesh.push_back(Facet_3D(b_end,t_p_middle,t_end));
            mesh.push_back(Facet_3D(b_end,b_p_middle,t_p_middle));
            mesh.push_back(Facet_3D(b_end,t_end,t_n_middle));
            mesh.push_back(Facet_3D(b_end,t_n_middle,b_n_middle));
        }
    }
    
    void m_sphere(Mesh_3D& mesh, const Point_3D::Measurement radius, 
            const int steps_per_quarter, const Point_3D& origin)
    {
        const Point_3D::Angle_Meas pi = 3.1415926535897932384;
        const Point_3D::Angle_Meas two_pi = pi * 2;
        const Point_3D::Angle_Meas step = pi / (2 * steps_per_quarter);
        const Point_3D::Measurement orig_x = origin.get_x();
        const Point_3D::Measurement orig_y = origin.get_y();
        const Point_3D::Measurement orig_z = origin.get_z();

        shared_ptr<Point_3D> top(new Point_3D(orig_x, orig_y, orig_z + radius));
        shared_ptr<Point_3D> bottom(new Point_3D(orig_x, orig_y, orig_z - radius));
        queue<shared_ptr<Point_3D>> init_vert_pts;
        queue<shared_ptr<Point_3D>> vert_pts; // vertical points
        // angular measuements
        Point_3D::Angle_Meas a_theta = 0;
        Point_3D::Angle_Meas a_phi = step;
        // fill queue with initial vertical set of points
        int num_phi_pts = (2 * steps_per_quarter) - 1;
        for (int count = 0; count < num_phi_pts; ++count, a_phi += step)
        {
            vert_pts.push(shared_ptr<Point_3D>(new Point_3D(spherical_point(radius, a_theta, a_phi, origin))));
            init_vert_pts.push(vert_pts.back());
        }
        // now go through and generate facets
        int num_theta_pts = (4 * steps_per_quarter) - 1;
        a_theta += step;
        for (int count = 0; count < num_theta_pts; ++count, a_theta += step)
        {
            a_phi = step;
            shared_ptr<Point_3D> p(new Point_3D(spherical_point(radius, a_theta, a_phi, origin)));
            vert_pts.push(p);
            mesh.push_back(Facet_3D(top, vert_pts.front(), p)); // add top level facet
            a_phi += step;
            // start count at 1 because created the initial point just above (p)
            for (int phi_count = 1; phi_count < num_phi_pts; ++phi_count, a_phi += step)
            {
                vert_pts.push(shared_ptr<Point_3D>(new Point_3D(spherical_point(radius, a_theta, a_phi, origin))));
                shared_ptr<Point_3D> p_p(vert_pts.front()); // previous p
                vert_pts.pop();
                mesh.push_back(Facet_3D(p, p_p, vert_pts.front()));
                mesh.push_back(Facet_3D(p, vert_pts.front(), vert_pts.back()));
                p = vert_pts.back();
            }
            mesh.push_back(Facet_3D(bottom, vert_pts.back(), vert_pts.front())); // add bottom level facet
            vert_pts.pop();
        }
        // add last remaining facets
        mesh.push_back(Facet_3D(top, vert_pts.front(), init_vert_pts.front())); // top level facet
        while (vert_pts.size() > 1)
        {
            shared_ptr<Point_3D> p(vert_pts.front());
            vert_pts.pop();
            shared_ptr<Point_3D> i_p(init_vert_pts.front());
            init_vert_pts.pop();
            mesh.push_back(Facet_3D(i_p, p, vert_pts.front()));
            mesh.push_back(Facet_3D(i_p, vert_pts.front(), init_vert_pts.front()));
        }
        mesh.push_back(Facet_3D(bottom, init_vert_pts.front(), vert_pts.front())); // bottom level facet
    }
    
    void m_ellipsoid(Mesh_3D& mesh, const Point_3D::Measurement x_radius, 
            const Point_3D::Measurement y_radius, const Point_3D::Measurement z_radius,
            const int steps_per_quarter, const Point_3D& origin)
    {
        // Ax^2 + By^2 + Cz^2 = D
        // x_radius = sqrt(D / A)
        // c_radius = sqrt(D / B)
        // z_radius = sqrt(D / C)
        // choose D to be A*B*C
        // x_radius = sqrt(B*C)
        // B*C = x_radius^2
        // y_radius = sqrt(A*C)
        // A*C = y_radius^2
        // z_radius = sqrt(A*B)
        // A*B = z_radius^2
        // solve for A,B, and C
        // A = z_radius^2 / B
        // substitute
        // (z_radius^2 / B) * C = y_radius^2
        // C = x_radius^2 / B
        // substitute
        // (z_radius^2 / B) * (x_radius^2 / B) = y_radius^2
        // B^2 = (z_radius^2 * x_radius^2) / y_radius^2
        // B = sqrt((z_radius^2 * x_radius^2) / y_radius^2)
        // C = x_radius^2 / B
        // A = z_radius^2 / B
        Point_3D::Measurement y_mult = sqrt((pow(z_radius,2) * pow(x_radius,2)) / pow(y_radius,2)); // B
        Point_3D::Measurement z_mult = pow(x_radius,2) / y_mult; // C
        Point_3D::Measurement x_mult = pow(z_radius,2) / y_mult; // A
        Point_3D::Measurement result = x_mult * y_mult * z_mult; // D
        // Ax^2 + By^2 + Cz^2 = D
        // find max z by setting x and y to zero
        // z = sqrt(D / C)
        const Point_3D::Measurement z_max = sqrt(result / z_mult);
        
        const Point_3D::Measurement z_step = z_max / steps_per_quarter;
        const Point_3D::Measurement orig_x = origin.get_x();
        const Point_3D::Measurement orig_y = origin.get_y();
        const Point_3D::Measurement orig_z = origin.get_z();
        
        // then starting at z=0, work to max z
        // for each level, start at -max_x and work to max x
        // Ax^2 + By^2 = D - z_step
        // y = +/- sqrt((D - z_step - Ax^2) / B)
        
        // keep track of previous layer points to create next layer of facets
        queue<shared_ptr<Point_3D>> pos_prev_layer; // positive side points
        queue<shared_ptr<Point_3D>> neg_prev_layer; // negative side points
        // start at bottom
        shared_ptr<Point_3D> z_max_pt(new Point_3D(orig_x,orig_y,orig_z - z_max)); // bottom most point
        Point_3D::Measurement z_pos = z_step - z_max; // z position of layer
        // find max x by setting y to zero and z to z_pos
        // Ax^2 + By^2 + Cz^2 = D
        // Ax^2 + Cz^2 = D
        // x_mult * x^2 + z_mult * z_pos^2 = result
        // x_mult * x^2 = result - z_mult * z_pos^2
        // x = sqrt((result - z_mult * z_pos^2) / x_mult)
//        Point_3D::Measurement x_max = sqrt((result - fabs(z_pos)) / x_mult);
        Point_3D::Measurement x_max = sqrt((result - z_mult * pow(z_pos, 2)) / x_mult);
        Point_3D::Measurement x_step = x_max / steps_per_quarter;
        Point_3D::Measurement x_pos = x_step - x_max; // x position
        shared_ptr<Point_3D> p_begin(new Point_3D(orig_x - x_max,orig_y,orig_z + z_pos));
        shared_ptr<Point_3D> p_pos_pt(p_begin);
        shared_ptr<Point_3D> p_neg_pt(p_begin);
        // create bottom most facets
        int num_layer_points = (2 * steps_per_quarter) - 1;
        for (int count=0; count < num_layer_points; x_pos += x_step, ++count)
        {
            // determine y value
            // Ax^2 + By^2 + Cz^2 = D
            // y_mult * y^2 = result - x_mult * x_pos^2 - z_mult * z_pos^2
            // y = sqrt((result - x_mult * x_pos^2 - z_mult * z_pos^2) / y_mult)
//            Point_3D::Measurement y = sqrt((result - fabs(z_pos) - x_mult * pow(x_pos, 2)) / y_mult);
            Point_3D::Measurement y = sqrt((result - x_mult * pow(x_pos, 2) - z_mult * pow(z_pos, 2)) / y_mult);
            pos_prev_layer.push(shared_ptr<Point_3D>(new Point_3D(orig_x + x_pos, orig_y + y, orig_z + z_pos)));
            neg_prev_layer.push(shared_ptr<Point_3D>(new Point_3D(orig_x + x_pos, orig_y - y, orig_z + z_pos)));
            mesh.push_back(Facet_3D(z_max_pt, p_pos_pt, pos_prev_layer.back()));
            mesh.push_back(Facet_3D(z_max_pt, neg_prev_layer.back(), p_neg_pt));
            p_pos_pt = pos_prev_layer.back();
            p_neg_pt = neg_prev_layer.back();
        }
        shared_ptr<Point_3D> p_end(new Point_3D(orig_x + x_max,orig_y,orig_z + z_pos));
        mesh.push_back(Facet_3D(z_max_pt, p_pos_pt, p_end));
        mesh.push_back(Facet_3D(z_max_pt, p_end, p_neg_pt));
        z_pos += z_step; // go to next level
        // now generate facets in the middle between the z top and z bottom
        for (int count=1;count < num_layer_points;z_pos += z_step,++count)
        {
            x_max = sqrt((result - z_mult * pow(z_pos, 2)) / x_mult); //sqrt((result - fabs(z_pos)) / x_mult);
            x_step = x_max / steps_per_quarter;
            x_pos = x_step - x_max;
            // generate first points, pop prev layer first points
            shared_ptr<Point_3D> begin(new Point_3D(orig_x - x_max, orig_y, orig_z + z_pos));
            shared_ptr<Point_3D> pos_pt(begin);
            shared_ptr<Point_3D> neg_pt(begin);
            p_pos_pt = p_begin;
            p_neg_pt = p_begin;
            for (int layer_count=0;layer_count < num_layer_points; x_pos += x_step,++layer_count)
            {
                // generate new points and facets
                Point_3D::Measurement y = sqrt((result - x_mult * pow(x_pos, 2) - z_mult * pow(z_pos, 2)) / y_mult); //sqrt((result - fabs(z_pos) - x_mult * pow(x_pos, 2)) / y_mult);
                pos_prev_layer.push(shared_ptr<Point_3D>(new Point_3D(orig_x + x_pos, orig_y + y, orig_z + z_pos)));
                neg_prev_layer.push(shared_ptr<Point_3D>(new Point_3D(orig_x + x_pos, orig_y - y, orig_z + z_pos)));
                // generate positive side facets
                mesh.push_back(Facet_3D(pos_prev_layer.back(),p_pos_pt,pos_pt));
                mesh.push_back(Facet_3D(pos_prev_layer.back(),pos_prev_layer.front(),p_pos_pt));
                // generate negative side facets
                mesh.push_back(Facet_3D(neg_prev_layer.back(),neg_pt,p_neg_pt));
                mesh.push_back(Facet_3D(neg_prev_layer.back(),p_neg_pt,neg_prev_layer.front()));
                // update point values for next cycle
                p_pos_pt = pos_prev_layer.front();
                pos_prev_layer.pop();
                pos_pt = pos_prev_layer.back();
                p_neg_pt = neg_prev_layer.front();
                neg_prev_layer.pop();
                neg_pt = neg_prev_layer.back();
            }
            // end point
            shared_ptr<Point_3D> end(new Point_3D(orig_x + x_max, orig_y, orig_z + z_pos));
            // generate end facets
            mesh.push_back(Facet_3D(end, p_pos_pt, pos_pt));
            mesh.push_back(Facet_3D(end, p_end, p_pos_pt));
            mesh.push_back(Facet_3D(end, neg_pt, p_neg_pt));
            mesh.push_back(Facet_3D(end, p_neg_pt, p_end));
            // prepare for next cycle through
            p_begin = begin;
            p_end = end; 
        }
        // final top most z point
        z_max_pt = shared_ptr<Point_3D>(new Point_3D(orig_x,orig_y,orig_z + z_max));
        // first points
        p_pos_pt = p_begin;
        p_neg_pt = p_begin;
        while (!pos_prev_layer.empty())
        {
            mesh.push_back(Facet_3D(z_max_pt,pos_prev_layer.front(),p_pos_pt));
            mesh.push_back(Facet_3D(z_max_pt,p_neg_pt,neg_prev_layer.front()));
            p_pos_pt = pos_prev_layer.front();
            pos_prev_layer.pop();
            p_neg_pt = neg_prev_layer.front();
            neg_prev_layer.pop();
        }
        // last facets
        mesh.push_back(Facet_3D(z_max_pt,p_end,p_pos_pt));
        mesh.push_back(Facet_3D(z_max_pt,p_neg_pt,p_end));
    }
}
