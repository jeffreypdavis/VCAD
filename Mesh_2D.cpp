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
 * File:   Mesh_2D.cpp
 * Author: Jeffrey Davis
 */

#include "Mesh_2D.h"
#include <cfloat>
#include <algorithm>
#include <unordered_map>
#include <utility>
#include "Point_3D.h"
#include "Mesh_3D.h"

namespace VCAD_lib
{

    Mesh_2D::Facet_Find::Facet_Find(const Facet& facet) : facet_to_find(facet) {}
    
    bool Mesh_2D::Facet_Find::operator()(const Facet& facet)
    {
        return (facet.get_p1_index() == facet_to_find.get_p1_index() && facet.get_p2_index() == facet_to_find.get_p2_index() && facet.get_p3_index() == facet_to_find.get_p3_index()) || 
                (facet.get_p1_index() == facet_to_find.get_p2_index() && facet.get_p2_index() == facet_to_find.get_p3_index() && facet.get_p3_index() == facet_to_find.get_p1_index()) || 
                (facet.get_p1_index() == facet_to_find.get_p3_index() && facet.get_p2_index() == facet_to_find.get_p1_index() && facet.get_p3_index() == facet_to_find.get_p2_index());
    }
    
    Mesh_2D::const_iterator::const_iterator(const vector<shared_ptr<Point_2D>>::const_iterator point_it_begin, 
            const vector<Facet>::const_iterator facet_it_begin, const vector<Facet>::const_iterator facet_it_end, 
            const vector<Facet>::const_iterator position)
            : point_list_begin(point_it_begin), facets_begin(facet_it_begin), 
            facets_end(facet_it_end), current_facet(position), facet() 
    {
        if (current_facet != facets_end)
            update_facet();
    }
    
    void Mesh_2D::const_iterator::update_facet()
    {
        if (current_facet == facets_end)
            throw runtime_error("iterator is at end of mesh");
        else // update facet
            facet = Facet_2D(
                    *(point_list_begin + current_facet->get_p1_index()), 
                    *(point_list_begin + current_facet->get_p2_index()), 
                    *(point_list_begin + current_facet->get_p3_index()));
    }
    
    bool Mesh_2D::const_iterator::operator==(const const_iterator& other_it) const
    {
        return current_facet == other_it.current_facet;
    }
    
    bool Mesh_2D::const_iterator::operator!=(const const_iterator& other_it) const
    {
        return current_facet != other_it.current_facet;
    }
    
    Mesh_2D::const_iterator& Mesh_2D::const_iterator::operator++()
    {
        // ++iterator
        ++current_facet;
        if (current_facet != facets_end)
            update_facet();
        return *this;
    }
    
    Mesh_2D::const_iterator Mesh_2D::const_iterator::operator++(int)
    {
        // iterator++
        const_iterator orig = *this;
        ++(*this);  // do the actual increment
        return orig;        
    }
    
    Mesh_2D::const_iterator& Mesh_2D::const_iterator::operator--()
    {
        // --iterator
        if (current_facet == facets_begin)
            throw runtime_error("cannot decrement past beginning");
        --current_facet;
        update_facet();
        return *this;
    }
    
    Mesh_2D::const_iterator Mesh_2D::const_iterator::operator--(int)
    {
        // iterator--
        const_iterator orig = *this;
        --(*this);  // do the actual increment
        return orig;        
    }
    
    const Facet_2D& Mesh_2D::const_iterator::operator*() const
    {
        if (current_facet >= facets_end)
            throw runtime_error("iterator is past last facet of mesh");
        return facet;
    }
    
    const Facet_2D* const Mesh_2D::const_iterator::operator->() const
    {
        if (current_facet >= facets_end)
            throw runtime_error("iterator is past last facet of mesh");
        return &facet;
    }
    
    Mesh_2D::Mesh_2D() : precision(DBL_EPSILON * 21), point_list(), facet_list() {}
    
    Mesh_2D::Mesh_2D(const Measurement prec) : precision(prec), point_list(), facet_list() {}

    Mesh_2D::Mesh_2D(const Mesh_2D& orig) : precision(orig.precision), point_list(), facet_list(orig.facet_list)
    {
        
        for (vector<shared_ptr<Point_2D>>::const_iterator iter = orig.point_list.begin(); iter < orig.point_list.end(); ++iter)
        {
            shared_ptr<Point_2D> ptr(new Point_2D(**iter));
            point_list.push_back(ptr);
        }
    }
    
    Mesh_2D& Mesh_2D::operator=(const Mesh_2D& other)
    {
        precision = other.precision;
        facet_list = other.facet_list;
        point_list.clear();
        for (vector<shared_ptr<Point_2D>>::const_iterator iter = other.point_list.begin(); iter < other.point_list.end(); ++iter)
        {
            shared_ptr<Point_2D> ptr(new Point_2D(**iter));
            point_list.push_back(ptr);
        }
        return *this;
    }
    
    Mesh_2D::const_iterator Mesh_2D::begin() const 
    {
        return const_iterator(point_list.begin(), facet_list.begin(), facet_list.end(), facet_list.begin()); 
    }
    
    Mesh_2D::const_iterator Mesh_2D::cbegin() const
    {
        return const_iterator(point_list.begin(), facet_list.begin(), facet_list.end(), facet_list.begin()); 
    }
    
    Mesh_2D::const_iterator Mesh_2D::end() const
    {
        return const_iterator(point_list.begin(), facet_list.begin(), facet_list.end(), facet_list.end()); 
    }
    
    Mesh_2D::const_iterator Mesh_2D::cend() const
    {
        return const_iterator(point_list.begin(), facet_list.begin(), facet_list.end(), facet_list.end()); 
    }
    
    void Mesh_2D::push_back(const Facet_2D& facet)
    {
        // find points if mesh already contains them
        int p1_index(-1);
        int p2_index(-1);
        int p3_index(-1);
        int index(0);
        for (vector<shared_ptr<Point_2D>>::const_iterator point_it = point_list.begin(); point_it != point_list.end(); ++point_it)
        {
            // test p1
            if ((*point_it)->get_x() == facet.get_point1()->get_x() && (*point_it)->get_y() == facet.get_point1()->get_y())
                p1_index = index;
            else if ((*point_it)->get_x() == facet.get_point2()->get_x() && (*point_it)->get_y() == facet.get_point2()->get_y())
                p2_index = index;
            else if ((*point_it)->get_x() == facet.get_point3()->get_x() && (*point_it)->get_y() == facet.get_point3()->get_y())
                p3_index = index;
            ++index;
        }
        
        // add points if not found
        if (p1_index == -1)
        {
            point_list.push_back(shared_ptr<Point_2D>(new Point_2D(*facet.get_point1())));
            p1_index = index;
            ++index;
        }
        if (p2_index == -1)
        {
            point_list.push_back(shared_ptr<Point_2D>(new Point_2D(*facet.get_point2())));
            p2_index = index;
            ++index;
        }
        if (p3_index == -1)
        {
            point_list.push_back(shared_ptr<Point_2D>(new Point_2D(*facet.get_point3())));
            p3_index = index;
            ++index;
        }
        
        // add facet
        facet_list.push_back(Facet(p1_index, p2_index, p3_index));
    }
    
    void Mesh_2D::clear()
    {
        facet_list.clear();
        point_list.clear();
    }
    
    Mesh_2D::const_iterator Mesh_2D::erase(const_iterator it)
    {
        // first find point indices
        int p1_index(-1);
        vector<shared_ptr<Point_2D>>::const_iterator p1_it(point_list.end());
        int p2_index(-1);
        vector<shared_ptr<Point_2D>>::const_iterator p2_it(point_list.end());
        int p3_index(-1);
        vector<shared_ptr<Point_2D>>::const_iterator p3_it(point_list.end());
        int index(0);
        for (vector<shared_ptr<Point_2D>>::const_iterator point_it = point_list.begin(); point_it != point_list.end(); ++point_it)
        {
            // test p1
            if ((*point_it)->get_x() == it->get_point1()->get_x() && (*point_it)->get_y() == it->get_point1()->get_y())
            {
                p1_index = index;
                p1_it = point_it;
            }
            else if ((*point_it)->get_x() == it->get_point2()->get_x() && (*point_it)->get_y() == it->get_point2()->get_y())
            {
                p2_index = index;
                p2_it = point_it;
            }
            else if ((*point_it)->get_x() == it->get_point3()->get_x() && (*point_it)->get_y() == it->get_point3()->get_y())
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
        
        // find if any other facets are using the same points
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
    
    Mesh_2D::const_iterator Mesh_2D::erase(const_iterator begin, const_iterator end)
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
    
//    Mesh_2D::const_iterator Mesh_2D::insert(const_iterator loc, const Facet_2D& facet)
//    {
//        // first find point indices
//        int loc_p1_index(-1);
//        int loc_p2_index(-1);
//        int loc_p3_index(-1);
//        int p1_index(-1);
//        int p2_index(-1);
//        int p3_index(-1);
//        int index(0);
//        for (vector<shared_ptr<Point_2D>>::const_iterator point_it = point_list.begin(); point_it != point_list.end(); ++point_it)
//        {
//            if ((*point_it)->get_x() == loc->get_point1()->get_x() && (*point_it)->get_y() == loc->get_point1()->get_y())
//                loc_p1_index = index;
//            else if ((*point_it)->get_x() == loc->get_point2()->get_x() && (*point_it)->get_y() == loc->get_point2()->get_y())
//                loc_p2_index = index;
//            else if ((*point_it)->get_x() == loc->get_point3()->get_x() && (*point_it)->get_y() == loc->get_point3()->get_y())
//                loc_p3_index = index;
//
//            if ((*point_it)->get_x() == facet.get_point1()->get_x() && (*point_it)->get_y() == facet.get_point1()->get_y())
//                p1_index = index;
//            else if ((*point_it)->get_x() == facet.get_point2()->get_x() && (*point_it)->get_y() == facet.get_point2()->get_y())
//                p2_index = index;
//            else if ((*point_it)->get_x() == facet.get_point3()->get_x() && (*point_it)->get_y() == facet.get_point3()->get_y())
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
//            point_list.push_back(shared_ptr<Point_2D>(new Point_2D(*facet.get_point1())));
//            p1_index = index;
//            ++index;
//        }
//        if (p2_index == -1)
//        {
//            point_list.push_back(shared_ptr<Point_2D>(new Point_2D(*facet.get_point2())));
//            p2_index = index;
//            ++index;
//        }
//        if (p3_index == -1)
//        {
//            point_list.push_back(shared_ptr<Point_2D>(new Point_2D(*facet.get_point3())));
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
//    Mesh_2D::const_iterator Mesh_2D::insert(const_iterator loc, const_iterator from_facet, const_iterator to_facet)
//    {
//        const_iterator position(loc);
//        const_iterator it(to_facet);
//        --it;
//        while (it != from_facet)
//            position = this->insert(position, *it);
//        return this->insert(position, *from_facet);
//    }

    Mesh_2D& Mesh_2D::rotate(const Angle_Meas angle)
    {
        for (vector<shared_ptr<Point_2D>>::const_iterator it = point_list.begin(); it != point_list.end(); ++it)
            (*it)->rotate(angle);
        return *this;
    }
    
    Mesh_2D& Mesh_2D::rotate(const Angle_Meas angle, const Point_2D& origin)
    {
        for (vector<shared_ptr<Point_2D>>::const_iterator it = point_list.begin(); it != point_list.end(); ++it)
            (*it)->rotate(angle, origin);
        return *this;
    }

    Mesh_2D& Mesh_2D::scale(const Measurement x_scalar, const Measurement y_scalar)
    {
        for (vector<shared_ptr<Point_2D>>::const_iterator it = point_list.begin(); it != point_list.end(); ++it)
            (*it)->scale(x_scalar, y_scalar);
        
        int num_neg(0);
        if (x_scalar < 0) { ++num_neg; }
        if (y_scalar < 0) { ++num_neg; }
        if (num_neg % 2 == 1) // odd number of negative multipliers
        {
            // invert the unit normal vector of each facet by swapping p2 and p3
            for (vector<Facet>::iterator it = facet_list.begin(); it != facet_list.end(); ++it)
                it->invert_unv();
        }
        
        return *this;
    }
    
    Mesh_2D& Mesh_2D::scale(const Measurement x_scalar, const Measurement y_scalar, 
            const Point_2D& origin)
    {
        for (vector<shared_ptr<Point_2D>>::const_iterator it = point_list.begin(); it != point_list.end(); ++it)
            (*it)->scale(x_scalar, y_scalar, origin);
        
        int num_neg(0);
        if (x_scalar < 0) { ++num_neg; }
        if (y_scalar < 0) { ++num_neg; }
        if (num_neg % 2 == 1) // odd number of negative multipliers
        {
            // invert the unit normal vector of each facet by swapping p2 and p3
            for (vector<Facet>::iterator it = facet_list.begin(); it != facet_list.end(); ++it)
                it->invert_unv();
        }
        
        return *this;
    }

    Mesh_2D& Mesh_2D::translate(const Measurement x_val, const Measurement y_val)
    {
        for (vector<shared_ptr<Point_2D>>::const_iterator it = point_list.begin(); it != point_list.end(); ++it)
            (*it)->translate(x_val, y_val);
        return *this;
    }
    
    Mesh_2D& Mesh_2D::translate(const Vector_2D& v)
    {
        for (vector<shared_ptr<Point_2D>>::const_iterator it = point_list.begin(); it != point_list.end(); ++it)
            (*it)->translate(v);
        return *this;
    }
    
    Mesh_2D& Mesh_2D::move(const Point_2D& new_origin, const Vector_2D& axis,
            const bool is_x_axis, const Point_2D& ref_origin)
    {
        for (vector<shared_ptr<Point_2D>>::const_iterator it = point_list.begin(); it != point_list.end(); ++it)
            (*it)->move(new_origin, axis, is_x_axis, ref_origin);
        return *this;
    }
    
    Mesh_2D& Mesh_2D::operator+=(const Vector_2D& v)
    {
        return this->translate(v);
    }
    
    Mesh_2D& Mesh_2D::operator-=(const Vector_2D& v)
    {
        Vector_2D temp(v);
        return this->translate(-temp);
    }
    
    Mesh_2D& Mesh_2D::operator*=(const Measurement scalar)
    {
        return this->scale(scalar, scalar);
    }
    
    Mesh_3D& linear_extrude(const Mesh_2D& mesh_2d, Mesh_3D& mesh, 
            const Mesh_3D::Measurement height, const bool center)
    {
        if (height == 0)
        {
            for (Mesh_2D::const_iterator it = mesh_2d.begin(); it != mesh_2d.end(); ++it)
                mesh.push_back(Facet_3D(shared_ptr<Point_3D>(new Point_3D(it->get_point1()->get_x(), it->get_point1()->get_y(), 0)), 
                        shared_ptr<Point_3D>(new Point_3D(it->get_point2()->get_x(), it->get_point2()->get_y(), 0)), 
                        shared_ptr<Point_3D>(new Point_3D(it->get_point3()->get_x(), it->get_point3()->get_y(), 0))));
            
            return mesh;
        }
        
        // determine outside line segments
        struct Facet_Side {
            int point1;
            int point2;
            Facet_Side(const int pt1, const int pt2) : point1(pt1), point2(pt2) {}
            const bool operator==(const Facet_Side& side) const
            {
                return (point1 == side.point1 || point1 == side.point2) && 
                        (point2 == side.point2 || point2 == side.point1);
            }
        };

        struct Facet_Side_Hasher {
            const int operator()(const Facet_Side& side) const
            {
                return side.point1 ^ side.point2;
            }
        };
        
        unordered_map<Facet_Side,int,Facet_Side_Hasher> facet_sides;
        
        for (Mesh_2D::const_facet_iterator it = mesh_2d.facet_begin(); it != mesh_2d.facet_end(); ++it)
        {
            // p1p2
            Facet_Side side(it->get_p1_index(), it->get_p2_index());
            unordered_map<Facet_Side,int,Facet_Side_Hasher>::iterator side_it = facet_sides.find(side);
            if (side_it == facet_sides.end()) 
            {
                facet_sides.insert(pair<Facet_Side,int>(side, 1));
            }
            else
            {
                side_it->second = side_it->second + 1;
            }
            
            // p1p3
            side = Facet_Side(it->get_p1_index(), it->get_p3_index());
            side_it = facet_sides.find(side);
            if (side_it == facet_sides.end()) 
            {
                facet_sides.insert(pair<Facet_Side,int>(side, 1));
            }
            else
            {
                side_it->second = side_it->second + 1;
            }
            
            // p2p3
            side = Facet_Side(it->get_p2_index(), it->get_p3_index());
            side_it = facet_sides.find(side);
            if (side_it == facet_sides.end()) 
            {
                facet_sides.insert(pair<Facet_Side,int>(side, 1));
            }
            else
            {
                side_it->second = side_it->second + 1;
            }
        }
        
        // find outside segments
        vector<Facet_Side> outside_segments;
        for (unordered_map<Facet_Side,int,Facet_Side_Hasher>::const_iterator it = facet_sides.begin(); it != facet_sides.end(); ++it)
        {
            if (it->second == 1)
                outside_segments.push_back(it->first);
        }
        facet_sides.clear();
        
        Point_2D::Measurement lower_z = center ? -height / 2 : 0;
        Point_2D::Measurement upper_z = center ? height / 2 : height;
        if (height < 0)
        {
            lower_z = center ? -lower_z : height;
            upper_z = center ? -upper_z : 0;
        }
        
        // generate mesh facets
        for (Mesh_2D::const_iterator it = mesh_2d.begin(); it != mesh_2d.end(); ++it)
        {
            // bottom facet
            mesh.push_back(Facet_3D(shared_ptr<Point_3D>(new Point_3D(it->get_point1()->get_x(), it->get_point1()->get_y(), lower_z)), 
                    shared_ptr<Point_3D>(new Point_3D(it->get_point3()->get_x(), it->get_point3()->get_y(), lower_z)), 
                    shared_ptr<Point_3D>(new Point_3D(it->get_point2()->get_x(), it->get_point2()->get_y(), lower_z))));
            
            // top facet
            mesh.push_back(Facet_3D(shared_ptr<Point_3D>(new Point_3D(it->get_point1()->get_x(), it->get_point1()->get_y(), upper_z)), 
                    shared_ptr<Point_3D>(new Point_3D(it->get_point2()->get_x(), it->get_point2()->get_y(), upper_z)), 
                    shared_ptr<Point_3D>(new Point_3D(it->get_point3()->get_x(), it->get_point3()->get_y(), upper_z))));
        }
        
        // generate side facets
        for (vector<Facet_Side>::const_iterator it = outside_segments.begin(); it != outside_segments.end(); ++it)
        {
            Mesh_2D::const_point_iterator p1_it = mesh_2d.point_begin();
            Mesh_2D::const_point_iterator p2_it = mesh_2d.point_begin();
            advance(p1_it, it->point1);
            advance(p2_it, it->point2);
            
            mesh.push_back(Facet_3D(shared_ptr<Point_3D>(new Point_3D((*p1_it)->get_x(), (*p1_it)->get_y(), lower_z)), 
                    shared_ptr<Point_3D>(new Point_3D((*p2_it)->get_x(), (*p2_it)->get_y(), upper_z)), 
                    shared_ptr<Point_3D>(new Point_3D((*p1_it)->get_x(), (*p1_it)->get_y(), upper_z))));
            mesh.push_back(Facet_3D(shared_ptr<Point_3D>(new Point_3D((*p1_it)->get_x(), (*p1_it)->get_y(), lower_z)), 
                    shared_ptr<Point_3D>(new Point_3D((*p2_it)->get_x(), (*p2_it)->get_y(), lower_z)), 
                    shared_ptr<Point_3D>(new Point_3D((*p2_it)->get_x(), (*p2_it)->get_y(), upper_z))));
        }
        
        return mesh;
    }
}
