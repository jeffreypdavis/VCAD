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
 * File:   Valid_Mesh_2D.cpp
 * Author: Jeffrey Davis
 */

#include "Valid_Mesh_2D.h"

#include <algorithm>
#include <functional>
#include <unordered_set>
#include <utility>
#include "Mesh_2D.h"
#include "Point_3D.h"

namespace VCAD_lib
{

    Valid_Mesh_2D::Facet_Side::Facet_Side(const int p1, const int p2) : point1(p1), point2(p2) {}
    
    const bool Valid_Mesh_2D::Facet_Side::operator ==(const Facet_Side& side) const
    {
        return (point1 == side.point1 || point1 == side.point2) && (point2 == side.point2 || point2 == side.point1);
    }
    
    const int Valid_Mesh_2D::FS_Hasher::operator ()(const Facet_Side& side) const
    {
        return side.point1 ^ side.point2;
    }
    
    const int Valid_Mesh_2D::Facet_2D_Hasher::operator ()(const Facet_2D& facet) const
    {
        hash<Point_3D::Measurement> hasher;
        
        shared_ptr<Point_2D> pt = facet.get_point1();
        int hash_value = (31 * hasher(pt->get_x())) ^ (43 * hasher(pt->get_y()));
        pt = facet.get_point2();
        hash_value = hash_value ^ (31 * hasher(pt->get_x())) ^ (43 * hasher(pt->get_y()));
        pt = facet.get_point3();
        hash_value = hash_value ^ (31 * hasher(pt->get_x())) ^ (43 * hasher(pt->get_y()));
        return hash_value;
    }
    
    const bool Valid_Mesh_2D::Facet_2D_Predicate::operator ()(const Facet_2D& facet1, const Facet_2D& facet2) const
    {
        return (facet1.get_point1() == facet2.get_point1()) &&
                (facet1.get_point2() == facet2.get_point2()) &&
                (facet1.get_point3() == facet2.get_point3());
    }
    
    const Vector_2D::Measurement Valid_Mesh_2D::Facet_Sorter::area(const Facet_2D& f) const
    {
        Vector_2D p1p2(*f.get_point1(), *f.get_point2());
        Vector_2D p1p3(*f.get_point1(), *f.get_point3());
        Vector_2D p2p3(*f.get_point2(), *f.get_point3());
        
        if (p1p2.length() >= p1p3.length() && p1p2.length() > p2p3.length())
        {
            // p1p2 is the largest side
            if (p1p3.length() > p2p3.length())
            {
                // p1p3 is the next largest
                return fabs(0.5 * cross_product(p1p2, p1p3));
            }
            else
            {
                // p2p3 is the next largest
                return fabs(0.5 * cross_product(-p1p2, p2p3));
            }
        }
        else if (p1p3.length() >= p1p2.length() && p1p3.length() >= p2p3.length())
        {
            // p1p3 is the largest side
            if (p1p2.length() > p2p3.length())
            {
                // p1p2 is the next largest
                return fabs(0.5 * cross_product(p1p3, p1p2));
            }
            else
            {
                // p2p3 is the next largest
                return fabs(0.5 * cross_product(-p1p3, -p2p3));
            }
        }
        else
        {
            // p2p3 is the largest side
            if (p1p2.length() > p1p3.length())
            {
                // p1p2 is the next largest
                return fabs(0.5 * cross_product(p2p3, -p1p2));
            }
            else
            {
                // p1p3 is the next largest
                return fabs(0.5 * cross_product(-p2p3, -p1p3));
            }
        }
    }
    
    const bool Valid_Mesh_2D::Facet_Sorter::operator ()(const Facet_2D& f1, const Facet_2D& f2) const
    {
        // put larger facets before smaller facets
        return area(f1) >= area(f2);
    }
    
    Valid_Mesh_2D::Valid_Mesh_2D(const Mesh_2D& mesh) : precision(mesh.get_precision()), 
            all_points(), all_facets(), facets_inside_facets(),
            pts_on_facet_sides(), too_many_share_side() 
    {
        for (Mesh_2D::const_point_iterator it = mesh.point_begin(); it != mesh.point_end(); ++it)
        {
            all_points.push_back(*it);
        }
        for (Mesh_2D::const_facet_iterator it = mesh.facet_begin(); it != mesh.facet_end(); ++it)
        {
            all_facets.push_back(*it);
        }
    }

    const bool Valid_Mesh_2D::validate()
    {
        unordered_map<Facet_Side,int,FS_Hasher> sides;
        unordered_set<int> points;
        
        // go through facets and get unique points and count the number of side occurrences
        for (vector<Facet>::const_iterator facet_it = all_facets.begin(); facet_it != all_facets.end(); ++facet_it)
        {
            points.insert(facet_it->get_p1_index());
            points.insert(facet_it->get_p2_index());
            points.insert(facet_it->get_p3_index());
            
            // take each side of each facet and insert into the map of sides
            // if the side already exists, update the side count
            // p1p2
            Facet_Side side(facet_it->get_p1_index(), facet_it->get_p2_index());
            unordered_map<Facet_Side,int,FS_Hasher>::iterator it = sides.find(side);
            if (it == sides.end())
                sides.insert(pair<Facet_Side,int>(side,1));
            else
                it->second = it->second + 1;

            // p1p3
            side = Facet_Side(facet_it->get_p1_index(), facet_it->get_p3_index());
            it = sides.find(side);
            if (it == sides.end())
                sides.insert(pair<Facet_Side,int>(side,1));
            else
                it->second = it->second + 1;
            
            // p2p3
            side = Facet_Side(facet_it->get_p2_index(), facet_it->get_p3_index());
            it = sides.find(side);
            if (it == sides.end())
                sides.insert(pair<Facet_Side,int>(side,1));
            else
                it->second = it->second + 1;
        }
        
        // now take each unique side and look for points that are on the side
        // but are not end points
        for (unordered_map<Facet_Side,int,FS_Hasher>::const_iterator side_it = sides.begin(); side_it != sides.end(); ++side_it)
        {
            if (side_it->second > 2)
            {
                vector<Facet_2D> facets;
                for (vector<Facet>::const_iterator it = all_facets.begin(); it != all_facets.end(); ++it)
                {
                    int count = 0;
                    if (it->get_p1_index() == side_it->first.point1 || it->get_p1_index() == side_it->first.point2)
                        ++count;
                    if (it->get_p2_index() == side_it->first.point1 || it->get_p2_index() == side_it->first.point2)
                        ++count;
                    if (it->get_p3_index() == side_it->first.point1 || it->get_p3_index() == side_it->first.point2)
                        ++count;
                    if (count == 2) {
                        vector<shared_ptr<Point_2D>>::const_iterator p1_it = all_points.begin();
                        vector<shared_ptr<Point_2D>>::const_iterator p2_it = all_points.begin();
                        vector<shared_ptr<Point_2D>>::const_iterator p3_it = all_points.begin();
                        advance(p1_it, it->get_p1_index());
                        advance(p2_it, it->get_p2_index());
                        advance(p3_it, it->get_p3_index());
                        facets.push_back(Facet_2D(*p1_it, *p2_it, *p3_it));
                    }
                }
                
                too_many_share_side.push_back(facets);
            }
            
            vector<shared_ptr<Point_2D>>::const_iterator p1_it = all_points.begin();
            advance(p1_it, side_it->first.point1);

            vector<shared_ptr<Point_2D>>::const_iterator p2_it = all_points.begin();
            advance(p2_it, side_it->first.point2);
            
            for (unordered_set<int>::const_iterator it = points.begin(); it != points.end(); ++it)
            {
                // go to next side if point is an end point
                if (side_it->first.point1 == *it || side_it->first.point2 == *it)
                    continue;

                vector<shared_ptr<Point_2D>>::const_iterator point_it = all_points.begin();
                advance(point_it, *it);
                
                if (is_pt_on_vector(**point_it, **p1_it, **p2_it, precision) && pts_on_facet_sides.end() == find(pts_on_facet_sides.begin(), pts_on_facet_sides.end(), *point_it))
                    pts_on_facet_sides.push_back(*point_it);
            }
        }
        
        // create a vector of facets then sort them largest to smallest
        vector<Facet_2D> sorted_facets;
        for (vector<Facet>::const_iterator it = all_facets.begin(); it != all_facets.end(); ++it)
        {
            vector<shared_ptr<Point_2D>>::const_iterator p1_it = all_points.begin();
            vector<shared_ptr<Point_2D>>::const_iterator p2_it = all_points.begin();
            vector<shared_ptr<Point_2D>>::const_iterator p3_it = all_points.begin();
            advance(p1_it, it->get_p1_index());
            advance(p2_it, it->get_p2_index());
            advance(p3_it, it->get_p3_index());
            
            sorted_facets.push_back(Facet_2D(*p1_it, *p2_it, *p3_it));
        }
        sort(sorted_facets.begin(), sorted_facets.end(), Facet_Sorter());

        // now go through sorted facets and find all smaller facets that are 
        // inside the larger facet
        for (vector<Facet_2D>::const_iterator iter = sorted_facets.begin(); iter != sorted_facets.end(); ++iter)
        {
            Facet_2D larger_facet(*iter);
            
            vector<Facet_2D> f_list;
            for (vector<Facet_2D>::const_iterator smaller_iter = iter + 1; smaller_iter != sorted_facets.end(); ++smaller_iter)
            {
                bool pt_on_side(false);
                if (larger_facet.contains_point(smaller_iter->get_inside_point(), pt_on_side, precision))
                {
                    // found facet inside facet
                    f_list.push_back(*smaller_iter);
                }
            }
            
            if (!f_list.empty())
                facets_inside_facets.insert(pair<Facet_2D,vector<Facet_2D>>(larger_facet,f_list));
        }
        
        return pts_on_facet_sides.empty() && too_many_share_side.empty() && facets_inside_facets.empty();
    }
}
