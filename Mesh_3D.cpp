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
 * File:   Mesh_3D.cpp
 * Author: Jeffrey Davis
 */

#include "Mesh_3D.h"
#include <cfloat>
#include <algorithm>
#include "Point_2D.h"
#include "Mesh_2D.h"

namespace VCAD_lib
{
    Mesh_3D::Facet_Find::Facet_Find(const Facet& facet) : facet_to_find(facet) {}
    
    bool Mesh_3D::Facet_Find::operator()(const Facet& facet)
    {
        return (facet.get_p1_index() == facet_to_find.get_p1_index() && facet.get_p2_index() == facet_to_find.get_p2_index() && facet.get_p3_index() == facet_to_find.get_p3_index()) || 
                (facet.get_p1_index() == facet_to_find.get_p2_index() && facet.get_p2_index() == facet_to_find.get_p3_index() && facet.get_p3_index() == facet_to_find.get_p1_index()) || 
                (facet.get_p1_index() == facet_to_find.get_p3_index() && facet.get_p2_index() == facet_to_find.get_p1_index() && facet.get_p3_index() == facet_to_find.get_p2_index());
    }
    
    Mesh_3D::const_iterator::const_iterator(const vector<shared_ptr<Point_3D>>::const_iterator point_it_begin, 
            const vector<Facet>::const_iterator facet_it_begin, const vector<Facet>::const_iterator facet_it_end, 
            const vector<Facet>::const_iterator position)
            : point_list_begin(point_it_begin), facets_begin(facet_it_begin), 
            facets_end(facet_it_end), current_facet(position), facet() 
    {
        if (current_facet != facets_end)
            update_facet();
    }
    
    void Mesh_3D::const_iterator::update_facet()
    {
        if (current_facet == facets_end)
            throw runtime_error("iterator is at end of mesh");
        else // update facet
            facet = Facet_3D(
                    *(point_list_begin + current_facet->get_p1_index()), 
                    *(point_list_begin + current_facet->get_p2_index()), 
                    *(point_list_begin + current_facet->get_p3_index()));
    }
    
    bool Mesh_3D::const_iterator::operator==(const const_iterator& other_it) const
    {
        return current_facet == other_it.current_facet;
    }
    
    bool Mesh_3D::const_iterator::operator!=(const const_iterator& other_it) const
    {
        return current_facet != other_it.current_facet;
    }
    
    Mesh_3D::const_iterator& Mesh_3D::const_iterator::operator++()
    {
        // ++iterator
        ++current_facet;
        if (current_facet != facets_end)
            update_facet();
        return *this;
    }
    
    Mesh_3D::const_iterator Mesh_3D::const_iterator::operator++(int)
    {
        // iterator++
        const_iterator orig = *this;
        ++(*this);  // do the actual increment
        return orig;        
    }
    
    Mesh_3D::const_iterator& Mesh_3D::const_iterator::operator--()
    {
        // --iterator
        if (current_facet == facets_begin)
            throw runtime_error("cannot decrement past beginning");
        --current_facet;
        update_facet();
        return *this;
    }
    
    Mesh_3D::const_iterator Mesh_3D::const_iterator::operator--(int)
    {
        // iterator--
        const_iterator orig = *this;
        --(*this);  // do the actual increment
        return orig;        
    }
    
    const Facet_3D& Mesh_3D::const_iterator::operator*() const
    {
        if (current_facet >= facets_end)
            throw runtime_error("iterator is past last facet of mesh");
        return facet;
    }
    
    const Facet_3D* const Mesh_3D::const_iterator::operator->() const
    {
        if (current_facet >= facets_end)
            throw runtime_error("iterator is past last facet of mesh");
        return &facet;
    }
    
    Mesh_3D::Mesh_3D() : precision(DBL_EPSILON * 21), point_list(), facet_list() {}
    
    Mesh_3D::Mesh_3D(const Measurement prec) : precision(prec), point_list(), facet_list() {}

    Mesh_3D::Mesh_3D(const Mesh_3D& orig) : precision(orig.precision), point_list(), facet_list(orig.facet_list)
    {
        
        for (vector<shared_ptr<Point_3D>>::const_iterator iter = orig.point_list.begin(); iter < orig.point_list.end(); ++iter)
        {
            shared_ptr<Point_3D> ptr(new Point_3D(**iter));
            point_list.push_back(ptr);
        }
    }
    
    Mesh_3D& Mesh_3D::operator=(const Mesh_3D& other)
    {
        precision = other.precision;
        facet_list = other.facet_list;
        point_list.clear();
        for (vector<shared_ptr<Point_3D>>::const_iterator iter = other.point_list.begin(); iter < other.point_list.end(); ++iter)
        {
            shared_ptr<Point_3D> ptr(new Point_3D(**iter));
            point_list.push_back(ptr);
        }
        return *this;
    }
    
    Mesh_3D::const_iterator Mesh_3D::begin() const 
    {
        return const_iterator(point_list.begin(), facet_list.begin(), facet_list.end(), facet_list.begin()); 
    }
    
    Mesh_3D::const_iterator Mesh_3D::cbegin() const
    {
        return const_iterator(point_list.begin(), facet_list.begin(), facet_list.end(), facet_list.begin()); 
    }
    
    Mesh_3D::const_iterator Mesh_3D::end() const
    {
        return const_iterator(point_list.begin(), facet_list.begin(), facet_list.end(), facet_list.end()); 
    }
    
    Mesh_3D::const_iterator Mesh_3D::cend() const
    {
        return const_iterator(point_list.begin(), facet_list.begin(), facet_list.end(), facet_list.end()); 
    }
    
    void Mesh_3D::push_back(const Facet_3D& facet)
    {
        // find points if mesh already contains them
        int p1_index(-1);
        int p2_index(-1);
        int p3_index(-1);
        int index(0);
        for (vector<shared_ptr<Point_3D>>::const_iterator point_it = point_list.begin(); point_it != point_list.end(); ++point_it)
        {
            // test p1
            if ((*point_it)->get_x() == facet.get_point1()->get_x() && (*point_it)->get_y() == facet.get_point1()->get_y() && (*point_it)->get_z() == facet.get_point1()->get_z())
                p1_index = index;
            else if ((*point_it)->get_x() == facet.get_point2()->get_x() && (*point_it)->get_y() == facet.get_point2()->get_y() && (*point_it)->get_z() == facet.get_point2()->get_z())
                p2_index = index;
            else if ((*point_it)->get_x() == facet.get_point3()->get_x() && (*point_it)->get_y() == facet.get_point3()->get_y() && (*point_it)->get_z() == facet.get_point3()->get_z())
                p3_index = index;
            ++index;
        }
        
        // add points if not found
        if (p1_index == -1)
        {
            point_list.push_back(shared_ptr<Point_3D>(new Point_3D(*facet.get_point1())));
            p1_index = index;
            ++index;
        }
        if (p2_index == -1)
        {
            point_list.push_back(shared_ptr<Point_3D>(new Point_3D(*facet.get_point2())));
            p2_index = index;
            ++index;
        }
        if (p3_index == -1)
        {
            point_list.push_back(shared_ptr<Point_3D>(new Point_3D(*facet.get_point3())));
            p3_index = index;
            ++index;
        }
        
        // add facet
        facet_list.push_back(Facet(p1_index, p2_index, p3_index));
    }
    
    void Mesh_3D::clear()
    {
        facet_list.clear();
        point_list.clear();
    }
    
    Mesh_3D::const_iterator Mesh_3D::erase(const_iterator it)
    {
        // first find point indices
        int p1_index(-1);
        vector<shared_ptr<Point_3D>>::const_iterator p1_it(point_list.end());
        int p2_index(-1);
        vector<shared_ptr<Point_3D>>::const_iterator p2_it(point_list.end());
        int p3_index(-1);
        vector<shared_ptr<Point_3D>>::const_iterator p3_it(point_list.end());
        int index(0);
        for (vector<shared_ptr<Point_3D>>::const_iterator point_it = point_list.begin(); point_it != point_list.end(); ++point_it)
        {
            // test p1
            if ((*point_it)->get_x() == it->get_point1()->get_x() && (*point_it)->get_y() == it->get_point1()->get_y() && (*point_it)->get_z() == it->get_point1()->get_z())
            {
                p1_index = index;
                p1_it = point_it;
            }
            else if ((*point_it)->get_x() == it->get_point2()->get_x() && (*point_it)->get_y() == it->get_point2()->get_y() && (*point_it)->get_z() == it->get_point2()->get_z())
            {
                p2_index = index;
                p2_it = point_it;
            }
            else if ((*point_it)->get_x() == it->get_point3()->get_x() && (*point_it)->get_y() == it->get_point3()->get_y() && (*point_it)->get_z() == it->get_point3()->get_z())
            {
                p3_index = index;
                p3_it = point_it;
            }
            ++index;
        }
        if (p1_index == -1 || p2_index == -1 || p3_index == -1)
            throw runtime_error("unable to locate facet points");
        // find facet in facet list
        vector<Facet>::const_iterator facet = find_if(facet_list.begin(), facet_list.end(), Facet_Find(Facet(p1_index, p2_index, p3_index)));
        if (facet == facet_list.end())
            throw runtime_error("unable to locate facet");
        
//        // find if any other facets are using the same points
//        bool p1_found(false);
//        bool p2_found(false);
//        bool p3_found(false);
//        for (vector<Facet>::const_iterator facet_it = facet_list.begin(); facet_it != facet_list.end(); ++facet_it)
//        {
//            if (facet_it == facet)
//                continue; // don't process the facet to remove
//            if (p1_index == facet_it->get_p1_index() || p1_index == facet_it->get_p2_index() || p1_index == facet_it->get_p3_index())
//                p1_found = true;
//            if (p2_index == facet_it->get_p1_index() || p2_index == facet_it->get_p2_index() || p2_index == facet_it->get_p3_index())
//                p2_found = true;
//            if (p3_index == facet_it->get_p1_index() || p3_index == facet_it->get_p2_index() || p3_index == facet_it->get_p3_index())
//                p3_found = true;
//        }
//
//        // if the points are not used by other facets remove the points
//        // remove the largest point indices first to not invalidate earlier point iterators
//        if (p1_index > p2_index && p1_index > p3_index)
//        {
//            if (!p1_found)
//                point_list.erase(p1_it);
//            if (p2_index > p3_index)
//            {
//                if (!p2_found)
//                    point_list.erase(p2_it);
//                if (!p3_found)
//                    point_list.erase(p3_it);
//            }
//            else
//            {
//                if (!p3_found)
//                    point_list.erase(p3_it);
//                if (!p2_found)
//                    point_list.erase(p2_it);
//            }
//        }
//        else if (p2_index > p1_index && p2_index > p3_index)
//        {
//            if (!p2_found)
//                point_list.erase(p2_it);
//            if (p1_index > p3_index)
//            {
//                if (!p1_found)
//                    point_list.erase(p1_it);
//                if (!p3_found)
//                    point_list.erase(p3_it);
//            }
//            else
//            {
//                if (!p3_found)
//                    point_list.erase(p3_it);
//                if (!p1_found)
//                    point_list.erase(p1_it);
//            }
//        }
//        else
//        {
//            if (!p3_found)
//                point_list.erase(p3_it);
//            if (p1_index > p2_index)
//            {
//                if (!p1_found)
//                    point_list.erase(p1_it);
//                if (!p2_found)
//                    point_list.erase(p2_it);
//            }
//            else
//            {
//                if (!p2_found)
//                    point_list.erase(p2_it);
//                if (!p1_found)
//                    point_list.erase(p1_it);
//            }
//        }
        
        // remove the facet
        facet_list.erase(facet);
    }
    
    Mesh_3D::const_iterator Mesh_3D::erase(const_iterator begin, const_iterator end)
    {
        const_iterator it(end);
        --it;
        while (it != begin)
        {
            this->erase(it);
            --it;
        }
        return this->erase(begin);
    }
    
//    Mesh_3D::const_iterator Mesh_3D::insert(const_iterator loc, const Facet_3D& facet)
//    {
//        // first find point indices
//        int loc_p1_index(-1);
//        int loc_p2_index(-1);
//        int loc_p3_index(-1);
//        int p1_index(-1);
//        int p2_index(-1);
//        int p3_index(-1);
//        int index(0);
//        for (vector<shared_ptr<Point_3D>>::const_iterator point_it = point_list.begin(); point_it != point_list.end(); ++point_it)
//        {
//            if ((*point_it)->get_x() == loc->get_point1()->get_x() && (*point_it)->get_y() == loc->get_point1()->get_y() && (*point_it)->get_z() == loc->get_point1()->get_z())
//                loc_p1_index = index;
//            else if ((*point_it)->get_x() == loc->get_point2()->get_x() && (*point_it)->get_y() == loc->get_point2()->get_y() && (*point_it)->get_z() == loc->get_point2()->get_z())
//                loc_p2_index = index;
//            else if ((*point_it)->get_x() == loc->get_point3()->get_x() && (*point_it)->get_y() == loc->get_point3()->get_y() && (*point_it)->get_z() == loc->get_point3()->get_z())
//                loc_p3_index = index;
//
//            if ((*point_it)->get_x() == facet.get_point1()->get_x() && (*point_it)->get_y() == facet.get_point1()->get_y() && (*point_it)->get_z() == facet.get_point1()->get_z())
//                p1_index = index;
//            else if ((*point_it)->get_x() == facet.get_point2()->get_x() && (*point_it)->get_y() == facet.get_point2()->get_y() && (*point_it)->get_z() == facet.get_point2()->get_z())
//                p2_index = index;
//            else if ((*point_it)->get_x() == facet.get_point3()->get_x() && (*point_it)->get_y() == facet.get_point3()->get_y() && (*point_it)->get_z() == facet.get_point3()->get_z())
//                p3_index = index;
//            
//            ++index;
//        }
//        if (loc_p1_index == -1 || loc_p2_index == -1 || loc_p3_index == -1)
//            throw runtime_error("unable to locate facet points");
//        // find facet in facet list
//        vector<Facet>::const_iterator loc_facet = find_if(facet_list.begin(), facet_list.end(), Facet_Find(Facet(loc_p1_index, loc_p2_index, loc_p3_index)));
//        if (loc_facet == facet_list.end())
//            throw runtime_error("unable to locate facet");
//
//        // add facet points that need to be added
//        if (p1_index == -1)
//        {
//            point_list.push_back(shared_ptr<Point_3D>(new Point_3D(*facet.get_point1())));
//            p1_index = index;
//            ++index;
//        }
//        if (p2_index == -1)
//        {
//            point_list.push_back(shared_ptr<Point_3D>(new Point_3D(*facet.get_point2())));
//            p2_index = index;
//            ++index;
//        }
//        if (p3_index == -1)
//        {
//            point_list.push_back(shared_ptr<Point_3D>(new Point_3D(*facet.get_point3())));
//            p3_index = index;
//        }
//        
//        // insert facet
//        vector<Facet>::const_iterator result = facet_list.insert(loc_facet, Facet(p1_index, p2_index, p3_index));
//
//        const_iterator it = begin();
//        for (vector<Facet>::const_iterator facet_it=facet_list.begin(); facet_it != facet_list.end(); ++facet_it)
//        {
//            if (result == facet_it)
//                break;
//            ++it;
//        }
//        return it;
//    }
//    
//    Mesh_3D::const_iterator Mesh_3D::insert(const_iterator loc, const_iterator from_facet, const_iterator to_facet)
//    {
//        const_iterator position(loc);
//        const_iterator it(to_facet);
//        --it;
//        while (it != from_facet)
//            position = this->insert(position, *it);
//        return this->insert(position, *from_facet);
//    }

    Mesh_3D& Mesh_3D::rotate(const Angle& angle)
    {
        for (vector<shared_ptr<Point_3D>>::const_iterator it = point_list.begin(); it != point_list.end(); ++it)
            (*it)->rotate(angle);
        return *this;
    }
    
    Mesh_3D& Mesh_3D::rotate(const Angle_Meas angle, const Vector_3D& axis)
    {
        for (vector<shared_ptr<Point_3D>>::const_iterator it = point_list.begin(); it != point_list.end(); ++it)
            (*it)->rotate(angle, axis);
        return *this;
    }
    
    Mesh_3D& Mesh_3D::rotate(const Angle& angle, const Point_3D& origin)
    {
        for (vector<shared_ptr<Point_3D>>::const_iterator it = point_list.begin(); it != point_list.end(); ++it)
            (*it)->rotate(angle, origin);
        return *this;
    }

    Mesh_3D& Mesh_3D::rotate(const Angle_Meas angle, const Vector_3D& axis, const Point_3D& origin)
    {
        for (vector<shared_ptr<Point_3D>>::const_iterator it = point_list.begin(); it != point_list.end(); ++it)
            (*it)->rotate(angle, axis, origin);
        return *this;
    }

    Mesh_3D& Mesh_3D::scale(const Measurement x_scalar, const Measurement y_scalar, 
            const Measurement z_scalar)
    {
        for (vector<shared_ptr<Point_3D>>::const_iterator it = point_list.begin(); it != point_list.end(); ++it)
            (*it)->scale(x_scalar, y_scalar, z_scalar);
        
        int num_neg(0);
        if (x_scalar < 0) { ++num_neg; }
        if (y_scalar < 0) { ++num_neg; }
        if (z_scalar < 0) { ++num_neg; }
        if (num_neg % 2 == 1) // odd number of negative multipliers
        {
            // invert the unit normal vector of each facet by swapping p2 and p3
            for (vector<Facet>::iterator it = facet_list.begin(); it != facet_list.end(); ++it)
                it->invert_unv();
        }
        
        return *this;
    }
    
    Mesh_3D& Mesh_3D::scale(const Measurement x_scalar, const Measurement y_scalar, 
            const Measurement z_scalar, const Point_3D& origin)
    {
        for (vector<shared_ptr<Point_3D>>::const_iterator it = point_list.begin(); it != point_list.end(); ++it)
            (*it)->scale(x_scalar, y_scalar, z_scalar, origin);
        
        int num_neg(0);
        if (x_scalar < 0) { ++num_neg; }
        if (y_scalar < 0) { ++num_neg; }
        if (z_scalar < 0) { ++num_neg; }
        if (num_neg % 2 == 1) // odd number of negative multipliers
        {
            // invert the unit normal vector of each facet by swapping p2 and p3
            for (vector<Facet>::iterator it = facet_list.begin(); it != facet_list.end(); ++it)
                it->invert_unv();
        }
        
        return *this;
    }

    Mesh_3D& Mesh_3D::translate(const Measurement x_val, const Measurement y_val, 
            const Measurement z_val)
    {
        for (vector<shared_ptr<Point_3D>>::const_iterator it = point_list.begin(); it != point_list.end(); ++it)
            (*it)->translate(x_val, y_val, z_val);
        return *this;
    }
    
    Mesh_3D& Mesh_3D::translate(const Vector_3D& v)
    {
        for (vector<shared_ptr<Point_3D>>::const_iterator it = point_list.begin(); it != point_list.end(); ++it)
            (*it)->translate(v);
        return *this;
    }
    
    Mesh_3D& Mesh_3D::move_x_pxy(const Point_3D& new_origin, const Vector_3D& x_axis, 
            const Point_3D& pt_xy_plane, const Point_3D& ref_origin)
    {
        for (vector<shared_ptr<Point_3D>>::const_iterator it = point_list.begin(); it != point_list.end(); ++it)
            (*it)->move_x_pxy(new_origin, x_axis, pt_xy_plane, ref_origin);
        return *this;
    }
    
    Mesh_3D& Mesh_3D::move_x_pxz(const Point_3D& new_origin, const Vector_3D& x_axis, 
            const Point_3D& pt_xz_plane, const Point_3D& ref_origin)
    {
        for (vector<shared_ptr<Point_3D>>::const_iterator it = point_list.begin(); it != point_list.end(); ++it)
            (*it)->move_x_pxz(new_origin, x_axis, pt_xz_plane, ref_origin);
        return *this;
    }
    
    Mesh_3D& Mesh_3D::move_y_pxy(const Point_3D& new_origin, const Vector_3D& y_axis, 
            const Point_3D& pt_xy_plane, const Point_3D& ref_origin)
    {
        for (vector<shared_ptr<Point_3D>>::const_iterator it = point_list.begin(); it != point_list.end(); ++it)
            (*it)->move_y_pxy(new_origin, y_axis, pt_xy_plane, ref_origin);
        return *this;
    }
    
    Mesh_3D& Mesh_3D::move_y_pyz(const Point_3D& new_origin, const Vector_3D& y_axis, 
            const Point_3D& pt_yz_plane, const Point_3D& ref_origin)
    {
        for (vector<shared_ptr<Point_3D>>::const_iterator it = point_list.begin(); it != point_list.end(); ++it)
            (*it)->move_y_pyz(new_origin, y_axis, pt_yz_plane, ref_origin);
        return *this;
    }
    
    Mesh_3D& Mesh_3D::move_z_pxz(const Point_3D& new_origin, const Vector_3D& z_axis, 
            const Point_3D& pt_xz_plane, const Point_3D& ref_origin)
    {
        for (vector<shared_ptr<Point_3D>>::const_iterator it = point_list.begin(); it != point_list.end(); ++it)
            (*it)->move_z_pxz(new_origin, z_axis, pt_xz_plane, ref_origin);
        return *this;
    }
    
    Mesh_3D& Mesh_3D::move_z_pyz(const Point_3D& new_origin, const Vector_3D& z_axis, 
            const Point_3D& pt_yz_plane, const Point_3D& ref_origin)
    {
        for (vector<shared_ptr<Point_3D>>::const_iterator it = point_list.begin(); it != point_list.end(); ++it)
            (*it)->move_z_pyz(new_origin, z_axis, pt_yz_plane, ref_origin);
        return *this;
    }
    
    Mesh_3D& Mesh_3D::operator+=(const Vector_3D& v)
    {
        return this->translate(v);
    }
    
    Mesh_3D& Mesh_3D::operator-=(const Vector_3D& v)
    {
        Vector_3D temp(v);
        return this->translate(-temp);
    }
    
    Mesh_3D& Mesh_3D::operator*=(const Measurement scalar)
    {
        return this->scale(scalar, scalar, scalar);
    }
    
    const bool mesh_contains_point(const Mesh_3D& mesh, const Point_3D& p, bool& pt_on_surface)
    {
        // determine if point is on or inside the mesh
        bool outside_mesh = false;
        for (Mesh_3D::const_iterator iter = mesh.begin(); iter != mesh.end(); ++iter)
        {
            // if point is on a facet, return true
            bool pt_on_side(false);
            if (iter->contains_point(p, pt_on_side, mesh.get_precision()))
            {
                pt_on_surface = true;
                return true;
            }
            else if (outside_mesh) // if already found to be outside, 
                continue;          // see if point is on a facet - no need to look if outside anymore
//            (*iter)->set_debug(false);
//            if (debug)
//                cout << "    after facet contains point\n";
            // determine if point is inside mesh by checking if 
            // point is behind each closest facet
            
            // vector from point to facet
            Point_3D inside_pt(iter->get_inside_point());
            Vector_3D v(p, inside_pt);
            bool closer_facet = false; // if there is a facet that is closer
            
            // go through mesh again looking for possible facets that are closer
            Point_3D i_point(0,0,0); // initialize here so it will be only assigned in loop
//            int count = 0;
            for (Mesh_3D::const_iterator it = mesh.begin(); it != mesh.end(); ++it)
            {
//                if (debug)
//                    cout << "    closer facet loop " << count++ << " of " << facets.size() << "\n";
                if (it == iter) // don't process the same facet
                    continue;
                // check if facet is in the same general direction and closer 
                // than the *iter facet
//                Vector_3D v2(p, (*it)->get_inside_point());
//                if (dot_product(v2, v) >= 0 && v2.length() < v.length())
//                {
                    // see if line v from p intersects *it facet
//                    if (intersect_line_facet_plane(v, p, **it, i_point, precision) && 
//                            (*it)->contains_point(i_point, precision))
                    if (intersect_line_facet_plane(v, p, *it, i_point, mesh.get_precision()) && 
                            is_pt_on_vector(i_point, p, inside_pt, mesh.get_precision()) && 
                            it->contains_point(i_point, pt_on_side, mesh.get_precision()))
                    {
                        // closer facet
                        closer_facet = true;
                        break;
                    }
//                }
            }
            
            if (!closer_facet) // this is the closest facet, so test
            {
//                if (debug)
//                    cout << "closest facet v.length zero? " << within_round(v.length(), 0, precision) << "\n";
                // if point is in front of facet
                if (dot_product(iter->get_unv(), v.normalize()) < 0)
                    outside_mesh = true;
            }
        }
        
        // if point is behind every closest facet, then the mesh contains the point
        if (!outside_mesh)
        {
            pt_on_surface = false;
            return true;
        }
        else
            return false;
    }
}
