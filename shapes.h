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
 * File:   shapes.h
 * Author: Jeffrey Davis
 */

#ifndef SHAPES_H
#define SHAPES_H

#include "Point_2D.h"
#include "Point_3D.h"
using namespace std;

namespace VCAD_lib
{
    class Mesh_2D;
    class Mesh_3D;
    
    // 2D shapes
    void m_rectangle(Mesh_2D& mesh, const Point_2D::Measurement x_length, 
            const Point_2D::Measurement y_length, const bool center=false, 
            const Point_2D& origin=Point_2D(0,0));
    
    void m_circle(Mesh_2D& mesh, const Point_2D::Measurement radius, 
            const int steps_per_quarter, const Point_2D& origin=Point_2D(0,0));
    
    void m_ellipse(Mesh_2D& mesh, const Point_2D::Measurement x_radius, 
            const Point_2D::Measurement y_radius, const int steps_per_quarter, 
            const Point_2D& origin=Point_2D(0,0));
    
    // 3D shapes
    void m_cuboid(Mesh_3D& mesh, const Point_3D::Measurement x_length, 
            const Point_3D::Measurement y_length, const Point_3D::Measurement z_length, 
            const bool center=false, const Point_3D& origin=Point_3D(0,0,0));
    
    /*
     * cylinder
     * mesh: the mesh to add the cylinder to
     * b_radius: the bottom radius
     * t_radius: the top radius
     * height: the height of the cylinder
     * steps_per_quarter: the number of steps per quarter to take when generating
     *                    points
     * center: if true, the height will be centered on the origin
     * origin: the point that should be considered the origin when forming the
     *         cylinder
     */
    void m_cylinder(Mesh_3D& mesh, const Point_3D::Measurement b_radius, 
            const Point_3D::Measurement t_radius, const Point_3D::Measurement height, 
            const int steps_per_quarter, const bool center=false, 
            const Point_3D& origin=Point_3D(0,0,0));
    
    /*
     * elliptical cylinder
     * mesh: the mesh to add the elliptical cylinder to
     * b_x_radius: the bottom x radius
     * b_y_radius: the bottom y radius
     * t_x_radius: the top x radius
     * t_y_radius: the top y radius
     * height: the height of the cylinder
     * steps_per_quarter: the number of steps per quarter to take when generating
     *                    points
     * center: if true, the height will be centered on the origin
     * origin: the point that should be considered the origin when forming the
     *         elliptical cylinder
     */
    void m_e_cylinder(Mesh_3D& mesh, const Point_3D::Measurement b_x_radius, 
            const Point_3D::Measurement b_y_radius, const Point_3D::Measurement t_x_radius, 
            const Point_3D::Measurement t_y_radius, const Point_3D::Measurement height, 
            const int steps_per_quarter, const bool center=false, 
            const Point_3D& origin=Point_3D(0,0,0));
    
    void m_sphere(Mesh_3D& mesh, const Point_3D::Measurement radius, 
            const int steps_per_quarter, const Point_3D& origin=Point_3D(0,0,0));
    
    void m_ellipsoid(Mesh_3D& mesh, const Point_3D::Measurement x_radius, 
            const Point_3D::Measurement y_radius, const Point_3D::Measurement z_radius,
            const int steps_per_quarter, const Point_3D& origin=Point_3D(0,0,0));
}

#endif /* SHAPES_H */

