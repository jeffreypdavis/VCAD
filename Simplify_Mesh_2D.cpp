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
 * File:   Simplify_Mesh_2D.cpp
 * Author: Jeffrey Davis
 * 
 * Available PreProcessor definitions
 * DEBUG_SIMPLIFY_MESH_2D
 * DEBUG_SIMPLIFY_MESH_2D_NEW_FACETS
 */

#include "Simplify_Mesh_2D.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include "Vector_2D.h"

namespace VCAD_lib
{
    
    Simplify_Mesh_2D::Pt_Remover::Segment::Segment(const int p1, const int p2, const bool internal) : point1(p1), point2(p2), is_internal(internal), used(false) {}
    
    const bool Simplify_Mesh_2D::Pt_Remover::Segment::operator==(const Segment& seg) const 
    {
        return (point1 == seg.point1 || point1 == seg.point2) && (point2 == seg.point1 || point2 == seg.point2);
    }
    
    const int Simplify_Mesh_2D::Pt_Remover::Segment::shares_pt(const Segment& seg) const 
    {
        
        if (point1 == seg.point1 || point1 == seg.point2)
            return point1;
        if (point2 == seg.point1 || point2 == seg.point2)
            return point2;
        return -1;
    }
    
    Simplify_Mesh_2D::Pt_Remover::Segments::Segment_Sort::Segment_Sort() {}

    const bool Simplify_Mesh_2D::Pt_Remover::Segments::Segment_Sort::operator ()(
            const Segment& seg1, const Segment& seg2) const
    {
        // move internal segments to the front
        return seg1.is_internal;
    }
    
    Simplify_Mesh_2D::Pt_Remover::Segments::Segments(const Mesh_2D& mesh) : 
            pt_begin(mesh.point_begin()), pt_end(mesh.point_end()), segments(), removed_segs() {}
    
    void Simplify_Mesh_2D::Pt_Remover::Segments::push_back(const int p1, const int p2, const bool internal)
    {
        Segment seg(p1,p2,internal);
        // there needs to be duplicates because that is how external segments
        // and internal segments are identified
        segments.push_back(seg);
        // sort segments so the internal segments move to the front
        sort(segments.begin(), segments.end(), Segment_Sort());
    }
    
    void Simplify_Mesh_2D::Pt_Remover::Segments::find_pts_to_remove(vector<int>& internal_pts, 
            vector<int>& perimeter_pts, vector<Facet>& same_plane_facets, const Mesh_2D& mesh) const
    {
        vector<Segment> internal_segs;
        vector<Segment> perimeter_segs;
        vector<int> all_perimeter_pts;
        // locate perimeter segment points and internal segments
        for (vector<Segment>::const_iterator it = segments.begin(); it != segments.end(); ++it)
        {
            if (internal_segs.end() != find(internal_segs.begin(), internal_segs.end(), *it))
                continue; // already processed internal segment, so move to next segment
            if (segments.end() == find(it + 1, segments.end(), *it))
            {
                // external segment because it is only found once
                if (all_perimeter_pts.end() == find(all_perimeter_pts.begin(), all_perimeter_pts.end(), it->point1))
                    all_perimeter_pts.push_back(it->point1);
                if (all_perimeter_pts.end() == find(all_perimeter_pts.begin(), all_perimeter_pts.end(), it->point2))
                    all_perimeter_pts.push_back(it->point2);
                perimeter_segs.push_back(*it);
            }
            else // internal segment because it appears twice
                internal_segs.push_back(*it);
        }
        
        // locate internal points that can be removed
        for (vector<Segment>::const_iterator it = internal_segs.begin(); it != internal_segs.end(); ++it)
        {
            if ((internal_pts.end() == find(internal_pts.begin(), internal_pts.end(), it->point1)) && 
                    (all_perimeter_pts.end() == find(all_perimeter_pts.begin(), all_perimeter_pts.end(), it->point1)))
                internal_pts.push_back(it->point1); // point is internal and can be removed
            if ((internal_pts.end() == find(internal_pts.begin(), internal_pts.end(), it->point2)) && 
                    (all_perimeter_pts.end() == find(all_perimeter_pts.begin(), all_perimeter_pts.end(), it->point2)))
                internal_pts.push_back(it->point2); // point is internal and can be removed
        }
        
        class Seg_Pt_Find {
        public:
            Seg_Pt_Find(int pt_index) : index(pt_index) {}
            const bool operator()(const Simplify_Mesh_2D::Pt_Remover::Segment& seg) const
            {
                return (seg.point1 == index) || (seg.point2 == index);
            }
        private:
            const int index;
        };
        
        // locate external points that might be able to be removed
        for (vector<int>::const_iterator it = all_perimeter_pts.begin(); it != all_perimeter_pts.end(); ++it)
        {
            // find two perimeter segments that share the same point
            Seg_Pt_Find seg_pt_find(*it);
            vector<Segment>::const_iterator seg_it1 = find_if(perimeter_segs.cbegin(), perimeter_segs.cend(), seg_pt_find);
            if (seg_it1 == perimeter_segs.cend())
                continue;
            vector<Segment>::const_iterator seg_it2 = find_if((seg_it1 + 1), perimeter_segs.cend(), seg_pt_find);
            if (seg_it2 == perimeter_segs.cend())
                continue;
            
            // test if segments are in a straight line
            
            int common_index = seg_it1->shares_pt(*seg_it2);
            if (common_index == -1)
                continue;
            Mesh_2D::const_point_iterator pt_iter(mesh.point_begin());
            advance(pt_iter, common_index);
            shared_ptr<Point_2D> p2(*pt_iter);
            pt_iter = mesh.point_begin();
            advance(pt_iter, (common_index == seg_it1->point1) ? seg_it1->point2 : seg_it1->point1);
            shared_ptr<Point_2D> p1(*pt_iter);
            pt_iter = mesh.point_begin();
            advance(pt_iter, (common_index == seg_it2->point1) ? seg_it2->point2 : seg_it2->point1);
            shared_ptr<Point_2D> p3(*pt_iter);
            
            bool same_direction(false);
            if (is_same_line(*p1, *p2, *p2, *p3, same_direction, mesh.get_precision()))
            {
                perimeter_pts.push_back(*it);
            }
        }
    }
    
    const Simplify_Mesh_2D::Pt_Remover::Segment* Simplify_Mesh_2D::Pt_Remover::Segments::get_next_segment(
            const Simplify_Mesh_2D::Pt_Remover::Segment* prev_segment) const
    {
        vector<Segment>::const_iterator it = segments.begin();
        if (prev_segment != 0)
        {
            while (it != segments.end())
            {
                if (*it == *prev_segment)
                {
                    ++it; // advance iterator to next segment
                    break;
                }

                ++it;
            }
        }
        if (it != segments.end())
            return &*it;
        
        return 0;
    }
    
    const Simplify_Mesh_2D::Pt_Remover::Segment* Simplify_Mesh_2D::Pt_Remover::Segments::find_connecting_segment(
            const Simplify_Mesh_2D::Pt_Remover::Segment& segment, const Simplify_Mesh_2D::Pt_Remover::Segment* prev_connecting_seg, 
            int& shared_pt, const Point_2D::Measurement precision) const
    {
#ifdef DEBUG_SIMPLIFY_MESH_2D_NEW_FACETS
        cout << "Simplify_Mesh_2D::Pt_Remover::Segments::find_connecting_segment begin\n";
#endif
        bool found(false);
        int shared_point;
#ifdef DEBUG_SIMPLIFY_MESH_2D_NEW_FACETS
        for (vector<Line_Segment>::const_iterator t_it = segments.begin(); t_it != segments.end(); ++t_it)
            cout << "Simplify_Mesh_2D::Pt_Remover::Segments::find_connecting_segment segments p1 x: " << (*t_it).point1->get_x() << " y: " <<
                    (*t_it).point1->get_y() << " p2 x: " << (*t_it).point2->get_x() << " y: " << (*t_it).point2->get_y() << "\n";
#endif
        vector<Segment>::const_iterator it = segments.begin();
        if (prev_connecting_seg != 0)
        {
#ifdef DEBUG_SIMPLIFY_MESH_2D_NEW_FACETS
            cout << "Simplify_Mesh_2D::Pt_Remover::Segments::find_conneting_segment locating last connecting segment in segments p1 x: " << 
                    prev_connecting_seg->point1->get_x() << " y: " << prev_connecting_seg->point1->get_y() << 
                    " p2 x: " << prev_connecting_seg->point2->get_x() << " y: " << prev_connecting_seg->point2->get_y() << "\n";
#endif
            while (it != segments.end())
            {
                if (*it == *prev_connecting_seg)
                {
#ifdef DEBUG_SIMPLIFY_MESH_2D_NEW_FACETS
                    cout << "Simplify_Mesh_2D::Pt_Remover::Segments::find_connecting_segment found segment\n";
#endif
                    found = true;
                    ++it; // move to next segment or end
                    break;
                }

                ++it;
            }
        }
        while (it != segments.end())
        {
#ifdef DEBUG_SIMPLIFY_MESH_2D_NEW_FACETS
            cout << "Simplify_Mesh_2D::Pt_Remover::Segments::find_connecting_segment checking segment p1 x: " << 
                    it->point1->get_x() << " y: " << it->point1->get_y() << " p2 x: " << it->point2->get_x() << 
                    " y: " << it->point2->get_y() << "\n";
#endif
            if (segment == *it) // do not return the same segment
            {
                ++it;
                continue;
            }

            int shared_point = segment.shares_pt(*it);
            if (shared_point != -1)
            {
#ifdef DEBUG_SIMPLIFY_MESH_2D_NEW_FACETS
                cout << "Simplify_Mesh_2D::Pt_Remover::Segments::find_connecting_segment segment shared pt x: " << 
                        shared_point->get_x() << " y: " << shared_point->get_y() << "\n";
#endif
                int index((segment.point1 == shared_point) ? segment.point2 : segment.point1);
                Mesh_2D::const_point_iterator pt_it(pt_begin);
                advance(pt_it, index);
                if (pt_it >= pt_end)
                    throw runtime_error("Invalid point index");
                shared_ptr<Point_2D> p1(*pt_it);
                pt_it = pt_begin;
                advance(pt_it, shared_point);
                if (pt_it >= pt_end)
                    throw runtime_error("Invalid point index");
                shared_ptr<Point_2D> p2(*pt_it);
                index = (it->point1 == shared_point) ? it->point2 : it->point1;
                pt_it = pt_begin;
                advance(pt_it, index);
                if (pt_it >= pt_end)
                    throw runtime_error("Invalid point index");
                shared_ptr<Point_2D> p3(*pt_it);

                bool same_direction(false);
                if (is_same_line(*p1, *p2, *p2, *p3, same_direction, precision))
                {
                    ++it; // same line so try next segment
                    continue;
                }

                // return segment
                shared_pt = shared_point;
                return &*it;
            }

            ++it;
        }
        
#ifdef DEBUG_SIMPLIFY_MESH_2D_NEW_FACETS
        cout << "Simplify_Mesh_2D::Pt_Remover::Segments::find_connecting_segment end returning 0\n";
#endif
        return 0;
    }

    const bool Simplify_Mesh_2D::Pt_Remover::Segments::is_seg_valid(const int p1, const int p2, 
            const vector<Facet_2D>& orig_facets, const Mesh_2D::const_point_iterator pt_begin, 
            const Mesh_2D::const_point_iterator pt_end, const Point_2D::Measurement precision) const
    {
#ifdef DEBUG_SIMPLIFY_MESH_2D_NEW_FACETS
        cout << "Simplify_Mesh_2D::Pt_Remover::Segments::is_seg_valid begin\n";
#endif
        // 1. is it an existing segment? is it removed already?
        if (segments.end() != find(segments.begin(), segments.end(), Segment(p1,p2,true)))
        {
#ifdef DEBUG_SIMPLIFY_MESH_2D_NEW_FACETS
            cout << "Simplify_Mesh_2D::Pt_Remover::Segments::is_seg_valid segment is an existing segment: returning true\n";
#endif
            return true;
        }
        if (removed_segs.end() != find(removed_segs.begin(), removed_segs.end(), Segment(p1,p2,true)))
        {
#ifdef DEBUG_SIMPLIFY_MESH_2D_NEW_FACETS
            cout << "Simplify_Mesh_2D::Pt_Remover::Segments::is_seg_valid segment is an already removed segment: returning false\n";
#endif
            return false; // segment has already been removed so it shouldn't be used again
        }
        // 2. does it intersect any existing or removed segments?
        // get segment p1 point value
        Segment seg(p1, p2, true);
        Mesh_2D::const_point_iterator pt_it(pt_begin);
        advance(pt_it, p1);
        if (pt_it >= pt_end)
            throw runtime_error("invalid point index");
        shared_ptr<Point_2D> p1_pt(*pt_it);
        
        // get segment p2 point value
        pt_it = pt_begin;
        advance(pt_it, p2);
        if (pt_it >= pt_end)
            throw runtime_error("invalid point index");
        shared_ptr<Point_2D> p2_pt(*pt_it);
        
#ifdef DEBUG_SIMPLIFY_MESH_2D_NEW_FACETS
        cout << "Simplify_Mesh_2D::Pt_Remover::Segments::is_seg_valid checking segment (" << 
                p1 << ", " << p2 << ") p1 x: " << p1_pt->get_x() << " y: " << p1_pt->get_y() << 
                " p2 x: " << p2_pt->get_x() << " y: " << p2_pt->get_y() << "\n";
#endif
        
        for (vector<Segment>::const_iterator it = segments.begin(); it != segments.end(); ++it)
        {
            // get segment points
            pt_it = pt_begin;
            advance(pt_it, it->point1);
            if (pt_it >= pt_end)
                throw runtime_error("invalid point index");
            shared_ptr<Point_2D> seg_p1(*pt_it);
            pt_it = pt_begin;
            advance(pt_it, it->point2);
            if (pt_it >= pt_end)
                throw runtime_error("invalid point index");
            shared_ptr<Point_2D> seg_p2(*pt_it);
            
#ifdef DEBUG_SIMPLIFY_MESH_2D_NEW_FACETS
            cout << "Simplify_Mesh_2D::Pt_Remover::Segments::is_seg_valid testing against existing segment (" << 
                    it->point1 << ", " << it->point2 << ") p1 x: " << seg_p1->get_x() << " y: " << seg_p1->get_y() << 
                    " p2 x: " << seg_p2->get_x() << " y: " << seg_p2->get_y() << "\n";
#endif
            
            // does the segment share a point
            int shared_pt(it->shares_pt(seg));
            if (shared_pt != -1) // shares a common point
            {
#ifdef DEBUG_SIMPLIFY_MESH_2D_NEW_FACETS
                cout << "Simplify_Mesh_2D::Pt_Remover::Segments::is_seg_valid shares common point: " << shared_pt << "\n";
#endif
                if (it->point1 == shared_pt)
                {
#ifdef DEBUG_SIMPLIFY_MESH_2D_NEW_FACETS
                    cout << "Simplify_Mesh_2D::Pt_Remover::Segments::is_seg_valid testing seg point1 is the shared point\n";
#endif
                    // test if non-common point is on vector
                    if (is_pt_on_vector(*seg_p2, *p1_pt, *p2_pt, precision))
                    {
#ifdef DEBUG_SIMPLIFY_MESH_2D_NEW_FACETS
                        cout << "Simplify_Mesh_2D::Pt_Remover::Segments::is_seg_valid testing segment point2 is on the main segment - returning false\n";
#endif
                        return false; // segment contains a smaller segment
                    }
                }
                else
                {
#ifdef DEBUG_SIMPLIFY_MESH_2D_NEW_FACETS
                    cout << "Simplify_Mesh_2D::Pt_Remover::Segments::is_seg_valid testing seg point2 is the shared point\n";
#endif
                    // test if non-common point is on vector
                    if (is_pt_on_vector(*seg_p1, *p1_pt, *p2_pt, precision))
                    {
#ifdef DEBUG_SIMPLIFY_MESH_2D_NEW_FACETS
                        cout << "Simplify_Mesh_2D::Pt_Remover::Segments::is_seg_valid testing segment point1 is on the main segment - returning false\n";
#endif
                        return false; // segment contains a smaller segment
                    }
                }
            }
            else // no common point, test if segments intersect
            {
#ifdef DEBUG_SIMPLIFY_MESH_2D_NEW_FACETS
                cout << "Simplify_Mesh_2D::Pt_Remover::Segments::is_seg_valid no shared points\n";
#endif
                Vector_2D_idata idata;
                if (intersect_vectors(*p1_pt, *p2_pt, *seg_p1, *seg_p2, idata, precision))
                {
#ifdef DEBUG_SIMPLIFY_MESH_2D_NEW_FACETS
                    cout << "Simplify_Mesh_2D::Pt_Remover::Segments::is_seg_valid segments intersect with no shared points - returning false\n";
#endif
                    return false; // segments intersect and do not share a point.  
                }
            }
        }
        for (vector<Segment>::const_iterator it = removed_segs.begin(); it != removed_segs.end(); ++it)
        {
            // get segment points
            pt_it = pt_begin;
            advance(pt_it, it->point1);
            if (pt_it >= pt_end)
                throw runtime_error("invalid point index");
            shared_ptr<Point_2D> seg_p1(*pt_it);
            pt_it = pt_begin;
            advance(pt_it, it->point2);
            if (pt_it >= pt_end)
                throw runtime_error("invalid point index");
            shared_ptr<Point_2D> seg_p2(*pt_it);
            
#ifdef DEBUG_SIMPLIFY_MESH_2D_NEW_FACETS
            cout << "Simplify_Mesh_2D::Pt_Remover::Segments::is_seg_valid testing against removed segment (" << 
                    it->point1 << ", " << it->point2 << ") p1 x: " << seg_p1->get_x() << " y: " << seg_p1->get_y() << 
                    " p2 x: " << seg_p2->get_x() << " y: " << seg_p2->get_y() << "\n";
#endif
            
            // does the segment share a point
            int shared_pt(it->shares_pt(seg));
            if (shared_pt != -1) // shares a common point
            {
#ifdef DEBUG_SIMPLIFY_MESH_2D_NEW_FACETS
                cout << "Simplify_Mesh_2D::Pt_Remover::Segments::is_seg_valid shares point: " << shared_pt << "\n";
#endif
                if (it->point1 == shared_pt)
                {
#ifdef DEBUG_SIMPLIFY_MESH_2D_NEW_FACETS
                    cout << "Simplify_Mesh_2D::Pt_Remover::Segments::is_seg_valid testing segment point1 is the shared point\n";
#endif
                    // test if non-common point is on vector
                    if (is_pt_on_vector(*seg_p2, *p1_pt, *p2_pt, precision))
                    {
#ifdef DEBUG_SIMPLIFY_MESH_2D_NEW_FACETS
                        cout << "Simplify_Mesh_2D::Pt_Remover::Segments::is_seg_valid testing segment point2 is on the original segment - returning false\n";
#endif
                        return false; // segment contains a smaller segment
                    }
                }
                else
                {
#ifdef DEBUG_SIMPLIFY_MESH_2D_NEW_FACETS
                    cout << "Simplify_Mesh_2D::Pt_Remover::Segments::is_seg_valid testing segment point2 is the shared point\n";
#endif
                    // test if non-common point is on vector
                    if (is_pt_on_vector(*seg_p1, *p1_pt, *p2_pt, precision))
                    {
#ifdef DEBUG_SIMPLIFY_MESH_2D_NEW_FACETS
                        cout << "Simplify_Mesh_2D::Pt_Remover::Segments::is_seg_valid testing segment point1 is on the original segment - returning false\n";
#endif
                        return false; // segment contains a smaller segment
                    }
                }
            }
            else // no common point, test if segments intersect
            {
#ifdef DEBUG_SIMPLIFY_MESH_2D_NEW_FACETS
                cout << "Simplify_Mesh_2D::Pt_Remover::Segments::is_seg_valid no shared points\n";
#endif
                Vector_2D_idata idata;
                if (intersect_vectors(*p1_pt, *p2_pt, *seg_p1, *seg_p2, idata, precision))
                {
#ifdef DEBUG_SIMPLIFY_MESH_2D_NEW_FACETS
                    cout << "Simplify_Mesh_2D::Pt_Remover::Segments::is_seg_valid no shared points, but segments intersect - returning false\n";
#endif
                    return false; // segments intersect and do not share a point.  
                }
            }
        }

        // 3. is a point halfway up segment inside the original facets?
        Vector_2D v(*p1_pt, *p2_pt);
        v *= 0.5;
        Point_2D half_way_pt(*p1_pt + v);
        
#ifdef DEBUG_SIMPLIFY_MESH_2D_NEW_FACETS
        cout << "Simplify_Mesh_2D::Pt_Remover::Segments::is_seg_valid half way pt x: " << half_way_pt.get_x() << " y: " << half_way_pt.get_y() << "\n";
#endif
        
        bool hwp_found(false); // half way point found
        for (vector<Facet_2D>::const_iterator it = orig_facets.begin(); it != orig_facets.end(); ++it)
        {
            bool pt_is_on_side(false);
            if (it->contains_point(half_way_pt, pt_is_on_side, precision))
                hwp_found = true;
        }
#ifdef DEBUG_SIMPLIFY_MESH_2D_NEW_FACETS
        if (hwp_found)
            cout << "Simplify_Mesh_2D::Pt_Remover::Segments::is_seg_valid original facet does contain half way point\n";
        else
            cout << "Simplify_Mesh_2D::Pt_Remover::Segments::is_seg_valid original facets does not contain half way point\n";
#endif
        return hwp_found; // segment is valid or not
    }
    
    void Simplify_Mesh_2D::Pt_Remover::Segments::process_segs(const Segment& seg1, const Segment& seg2, const int seg3_p1, const int seg3_p2)
    {
#ifdef DEBUG_SIMPLIFY_MESH_2D_NEW_FACETS
        cout << "Simplify_Mesh_2D::Pt_Remover::Segments::process_segs begin\n";
#endif
        Segment segment2(seg2);
        Segment segment3(seg3_p1, seg3_p2, true);
        
        vector<Segment>::iterator seg_it = find(segments.begin(), segments.end(), seg1);
        if (seg_it == segments.end())
            throw runtime_error("Unable to locate seg1");
        if (seg1.is_internal)
        {
#ifdef DEBUG_SIMPLIFY_MESH_2D_NEW_FACETS
        cout << "Simplify_Mesh_2D::Pt_Remover::Segments::process_segs seg1 is an internal segment\n";
#endif
            if (seg1.used)
            {
#ifdef DEBUG_SIMPLIFY_MESH_2D_NEW_FACETS
                cout << "Simplify_Mesh_2D::Pt_Remover::Segments::process_segs seg1 has been used before. removing\n";
#endif
                removed_segs.push_back(*seg_it);
                segments.erase(seg_it);
            }
            else
                seg_it->used = true;
        }
        else // seg1 is a perimeter segment
        {
#ifdef DEBUG_SIMPLIFY_MESH_2D_NEW_FACETS
            cout << "Simplify_Mesh_2D::Pt_Remover::Segments::process_segs seg1 is a perimeter segment. removing\n";
#endif
            removed_segs.push_back(*seg_it);
            segments.erase(seg_it);
        }
        
        seg_it = find(segments.begin(), segments.end(), segment2);
        if (seg_it == segments.end())
            throw runtime_error("Unable to locate segment2");
        if (seg_it->is_internal)
        {
#ifdef DEBUG_SIMPLIFY_MESH_2D_NEW_FACETS
            cout << "Simplify_Mesh_2D::Pt_Remover::Segments::process_segs segment2 is an internal segment\n";
#endif
            if (seg_it->used)
            {
#ifdef DEBUG_SIMPLIFY_MESH_2D_NEW_FACETS
                cout << "Simplify_Mesh_2D::Pt_Remover::Segments::process_segs segment2 has been used before. removing\n";
#endif
                removed_segs.push_back(*seg_it);
                segments.erase(seg_it);
            }
            else
                seg_it->used = true;
        }
        else // segment2 is a perimeter segment
        {
#ifdef DEBUG_SIMPLIFY_MESH_2D_NEW_FACETS
            cout << "Simplify_Mesh_2D::Pt_Remover::Segments::process_segs segment2 is a perimeter segment. removing\n";
#endif
            removed_segs.push_back(*seg_it);
            segments.erase(seg_it);
        }
        
        seg_it = find(segments.begin(), segments.end(), segment3);
        if (seg_it == segments.end()) // segment3 does not exist, so it must be internal
        {
#ifdef DEBUG_SIMPLIFY_MESH_2D_NEW_FACETS
            cout << "Simplify_Mesh_2D::Pt_Remover::Segments::process_segs segment3 is internal. setting used to true\n";
#endif
            // set used to true and add internal segment
            segment3.used = true;
            segments.push_back(segment3);
            sort(segments.begin(), segments.end(), Segment_Sort()); // sort segments
        }
        else // segment3 already exists
        {
            if (seg_it->is_internal)
            {
#ifdef DEBUG_SIMPLIFY_MESH_2D_NEW_FACETS
                cout << "Simplify_Mesh_2D::Pt_Remover::Segments::process_segs segment3 is an internal segment\n";
#endif
                if (seg_it->used)
                {
#ifdef DEBUG_SIMPLIFY_MESH_2D_NEW_FACETS
                    cout << "Simplify_Mesh_2D::Pt_Remover::Segments::process_segs segment3 has been used before. removing\n";
#endif
                    removed_segs.push_back(*seg_it);
                    segments.erase(seg_it);
                }
                else
                    seg_it->used = true;
            }
            else // segment3 is a perimeter segment
            {
#ifdef DEBUG_SIMPLIFY_MESH_2D_NEW_FACETS
                cout << "Simplify_Mesh_2D::Pt_Remover::Segments::process_segs segment3 is a perimeter segment. removing\n";
#endif
                removed_segs.push_back(*seg_it);
                segments.erase(seg_it);
            }
        }
    }

    Simplify_Mesh_2D::Pt_Remover::Pt_Remover(const Mesh_2D& mesh) : 
                orig_mesh(&mesh), internal_pts(), perimeter_pts(), same_plane_facets() 
    {
        // make a copy of the facets in the mesh
        for (Mesh_2D::const_facet_iterator it = mesh.facet_begin(); it != mesh.facet_end(); ++it)
            same_plane_facets.push_back(*it);
        // determine internal and perimeter points that can be removed
        Segments segments(mesh);
        for (vector<Facet>::const_iterator it = same_plane_facets.begin(); it != same_plane_facets.end(); ++it)
        {
            segments.push_back(it->get_p1_index(), it->get_p2_index(), false);
            segments.push_back(it->get_p1_index(), it->get_p3_index(), false);
            segments.push_back(it->get_p2_index(), it->get_p3_index(), false);
        }
        segments.find_pts_to_remove(internal_pts, perimeter_pts, same_plane_facets, mesh);
    }
    
    void Simplify_Mesh_2D::Pt_Remover::form_new_facets(const vector<Facet_2D>& orig_facets, 
            Segments& segments, vector<Facet>& new_facets) const
    {
#ifdef DEBUG_SIMPLIFY_MESH_2D_NEW_FACETS
        cout << "Simplify_Mesh_2D::Pt_Remover::form_new_facets begin\n";
        cout << "Simplify_Mesh_2D::Pt_Remover::form_new_facets Segments:\n";
        for (Segments::const_iterator it = segments.begin(); it != segments.end(); ++it)
        {
            Mesh_2D::const_point_iterator pt_it = orig_mesh->point_begin();
            advance(pt_it, it->point1);
            shared_ptr<Point_2D> p1(*pt_it);
            pt_it = orig_mesh->point_begin();
            advance(pt_it, it->point2);
            shared_ptr<Point_2D> p2(*pt_it);
            cout << "Simplify_Mesh_2D::Pt_Remover::form_new_facets (" << it->point1 << 
                    ", " << it->point2 << ") p1 x: " << p1->get_x() << " y: " << p1->get_y() << 
                    " p2 x: " << p2->get_x() << " y: " << p2->get_y() << "\n";
        }
#endif
        Vector_2D::Measurement plane_unv(cross_product(Vector_2D(*(orig_facets.begin()->get_point1()), *(orig_facets.begin()->get_point2())), 
                Vector_2D(*(orig_facets.begin()->get_point1()), *(orig_facets.begin()->get_point3()))));
        
        // form new facets
        const Segment* seg1 = segments.get_next_segment(0);
        while (seg1 != 0)
        {
            Segment segment1(*seg1);
#ifdef DEBUG_SIMPLIFY_MESH_2D_NEW_FACETS
            cout << "Simplify_Mesh_2D::Pt_Remover::form_new_facets segment1 p1: " << segment1.point1 << " p2: " << segment1.point2 << "\n";
#endif
            // get connecting segment
            int shared_pt(-1);
            const Segment* segment2 = segments.find_connecting_segment(*seg1, 0, shared_pt, orig_mesh->get_precision());
#ifdef DEBUG_SIMPLIFY_MESH_2D_NEW_FACETS
            cout << "Simplify_Mesh_2D::Pt_Remover::form_new_facets segment2 p1: " << segment2->point1 << " p2: " << segment2->point2 << "\n";
#endif
            
            while (segment2 != 0)
            {
                int seg3_p1 = seg1->point1 == shared_pt ? seg1->point2 : seg1->point1;
                int seg3_p2 = segment2->point1 == shared_pt ? segment2->point2 : segment2->point1;
#ifdef DEBUG_SIMPLIFY_MESH_2D_NEW_FACETS
                cout << "Simplify_Mesh_2D::Pt_Remover::form_new_facets formed segment3: (" << seg3_p1 << ", " << seg3_p2 << ")\n";
#endif
                if (segments.is_seg_valid(seg3_p1, seg3_p2, orig_facets, orig_mesh->point_begin(), orig_mesh->point_end(), orig_mesh->get_precision()))
                {
#ifdef DEBUG_SIMPLIFY_MESH_2D_NEW_FACETS
                    cout << "Simplify_Mesh_2D::Pt_Remover::form_new_facets seg3 is valid\n";
#endif
                    Facet facet(seg3_p1, shared_pt, seg3_p2);
                    // verify if facet unit normal is pointing in the right direction
                    Mesh_2D::const_point_iterator pt_it = orig_mesh->point_begin();
                    advance(pt_it, seg3_p1);
                    shared_ptr<Point_2D> p1(*pt_it);
                    pt_it = orig_mesh->point_begin();
                    advance(pt_it, shared_pt);
                    shared_ptr<Point_2D> p2(*pt_it);
                    pt_it = orig_mesh->point_begin();
                    advance(pt_it, seg3_p2);
                    shared_ptr<Point_2D> p3(*pt_it);
                    Vector_2D::Measurement unv(cross_product(Vector_2D(*p1, *p2), Vector_2D(*p1, *p3)));
                    if ((unv < 0 && plane_unv > 0) || (unv > 0 && plane_unv < 0)) // change facet point order
                        facet.invert_unv();
#ifdef DEBUG_SIMPLIFY_MESH_2D_NEW_FACETS
                    cout << "Simplify_Mesh_2D::Pt_Remover::form_new_facets Created new facet p1: " << facet.get_p1_index() << " p2: " << facet.get_p2_index() << " p3: " << facet.get_p3_index() << "\n";
#endif
                    new_facets.push_back(facet);
                    // process segments used
                    segments.process_segs(*seg1, *segment2, seg3_p1, seg3_p2);
                    break;
                }
                else // try to find another connecting segment
                {
#ifdef DEBUG_SIMPLIFY_MESH_2D_NEW_FACETS
                    cout << "Simplify_Mesh_2D::Pt_Remover::form_new_facets seg3 is not valid\n";
#endif
                    segment2 = segments.find_connecting_segment(*seg1, segment2, shared_pt, orig_mesh->get_precision());
                }
            }
            
            // get the next segment1
            Segments::const_iterator it = find(segments.begin(), segments.end(), segment1);
            if (it == segments.end())
                seg1 = segments.get_next_segment(0);
            else
                seg1 = segments.get_next_segment(&*it);
        }
    }
    
    void Simplify_Mesh_2D::Pt_Remover::rem_internal_pts()
    {
        for (vector<int>::const_iterator it = internal_pts.begin(); it != internal_pts.end(); ++it)
        {
#ifdef DEBUG_SIMPLIFY_MESH_2D
            Mesh_2D::const_point_iterator temp_pt_it = orig_mesh->point_begin();
            advance(temp_pt_it, *it);
            cout << "Simplify_Mesh_2D::Pt_Remover::rem_internal_pts removing internal point #" << *it << " x: " << (*temp_pt_it)->get_x() << " y: " << (*temp_pt_it)->get_y() << "\n";
#endif
            // form original facets surrounding internal point and the
            // perimeter segments surrounding the internal point
            vector<Facet> orig_facets;
            vector<Facet_2D> orig_facet_3ds;
            Segments perimeter_segs(*orig_mesh);
            vector<Facet>::const_iterator f_it = same_plane_facets.begin();
            while (f_it != same_plane_facets.end())
            {
                bool add(false);
                if (f_it->get_p1_index() == *it)
                {
                    perimeter_segs.push_back(f_it->get_p2_index(), f_it->get_p3_index(), false);
                    add = true;
                }
                else if (f_it->get_p2_index() == *it)
                {
                    perimeter_segs.push_back(f_it->get_p1_index(), f_it->get_p3_index(), false);
                    add = true;
                }
                else if (f_it->get_p3_index() == *it)
                {
                    perimeter_segs.push_back(f_it->get_p1_index(), f_it->get_p2_index(), false);
                    add = true;
                }
                if (add)
                {
#ifdef DEBUG_SIMPLIFY_MESH_2D
                    cout << "Simplify_Mesh_2D::Pt_Remover::rem_internal_pts adding facet p1: " << f_it->get_p1_index() << " p2: " << f_it->get_p2_index() << " p3: " << f_it->get_p3_index() << "\n";
#endif
                    orig_facets.push_back(*f_it);
                    shared_ptr<Point_2D> p1(*(orig_mesh->point_begin() + f_it->get_p1_index()));
                    shared_ptr<Point_2D> p2(*(orig_mesh->point_begin() + f_it->get_p2_index()));
                    shared_ptr<Point_2D> p3(*(orig_mesh->point_begin() + f_it->get_p3_index()));
                    orig_facet_3ds.push_back(Facet_2D(p1, p2, p3));
                    f_it = same_plane_facets.erase(f_it);
                }
                else
                    ++f_it;
            }
#ifdef DEBUG_SIMPLIFY_MESH_2D
            cout << "Simplify_Mesh_2D::Pt_Remover::rem_internal_pts forming new facets\n";
#endif            
            // remove internal point by forming new facets and replacing the original ones in same_plane_facets
            vector<Facet> new_facets;
            form_new_facets(orig_facet_3ds, perimeter_segs, new_facets);
            
            // add newly created facets back to the list of same plane facets
            for (vector<Facet>::const_iterator nf_it = new_facets.begin(); nf_it != new_facets.end(); ++nf_it)
            {
#ifdef DEBUG_SIMPLIFY_MESH_2D
                cout << "Simplify_Mesh_2D::Pt_Remover::rem_internal_pts adding newly created facet p1: " << nf_it->get_p1_index() << " p2: " << nf_it->get_p2_index() << " p3: " << nf_it->get_p3_index() << "\n";
#endif
                same_plane_facets.push_back(*nf_it);
            }
        }
    }
    
    void Simplify_Mesh_2D::Pt_Remover::rem_perimeter_pts()
    {
        for (vector<int>::const_iterator ppt_it = perimeter_pts.begin(); ppt_it != perimeter_pts.end(); ++ppt_it) 
        {
#ifdef DEBUG_SIMPLIFY_MESH_2D
            Mesh_2D::const_point_iterator temp_pt_it = orig_mesh->point_begin();
            advance(temp_pt_it, *ppt_it);
            cout << "Simplify_Mesh_2D::Pt_Remover::rem_perimeter_pts removing perimeter point #" << *ppt_it << " x: " << (*temp_pt_it)->get_x() << " y: " << (*temp_pt_it)->get_y() << "\n";
#endif
            // form original facets surrounding internal point and the
            // perimeter segments surrounding the perimeter point
            vector<Facet> orig_facets;
            vector<Facet_2D> orig_facet_2ds;
            Segments perimeter_segs(*orig_mesh);

            vector<Facet>::const_iterator f_it = same_plane_facets.begin();
            while (f_it != same_plane_facets.end())
            {
                bool add(false);
                if (f_it->get_p1_index() == *ppt_it)
                {
                    perimeter_segs.push_back(f_it->get_p2_index(), f_it->get_p3_index(), false);
                    add = true;
                }
                else if (f_it->get_p2_index() == *ppt_it)
                {
                    perimeter_segs.push_back(f_it->get_p1_index(), f_it->get_p3_index(), false);
                    add = true;
                }
                else if (f_it->get_p3_index() == *ppt_it)
                {
                    perimeter_segs.push_back(f_it->get_p1_index(), f_it->get_p2_index(), false);
                    add = true;
                }
                if (add)
                {
#ifdef DEBUG_SIMPLIFY_MESH_2D
                    cout << "Simplify_Mesh_2D::Pt_Remover::rem_perimeter_pts adding facet p1: " << f_it->get_p1_index() << " p2: " << f_it->get_p2_index() << " p3: " << f_it->get_p3_index() << "\n";
#endif
                    orig_facets.push_back(*f_it);
                    shared_ptr<Point_2D> p1(*(orig_mesh->point_begin() + f_it->get_p1_index()));
                    shared_ptr<Point_2D> p2(*(orig_mesh->point_begin() + f_it->get_p2_index()));
                    shared_ptr<Point_2D> p3(*(orig_mesh->point_begin() + f_it->get_p3_index()));
                    orig_facet_2ds.push_back(Facet_2D(p1, p2, p3));
                    f_it = same_plane_facets.erase(f_it);
                }
                else
                    ++f_it;
            }

            // locate segment end points and then form a new segment that connects them
            // this would remove that perimeter point
            vector<int> seg_pts;
            for (Segments::const_iterator seg_it = perimeter_segs.begin(); seg_it != perimeter_segs.end(); ++seg_it)
            {
                seg_pts.push_back(seg_it->point1);
                seg_pts.push_back(seg_it->point2);
            }
            vector<int>::const_iterator pt_it = seg_pts.cbegin();
            while (pt_it != seg_pts.cend())
            {
                int current_pt(*pt_it);
                vector<int>::const_iterator it = find(pt_it + 1, seg_pts.cend(), current_pt);
                if (it != seg_pts.cend())
                {
                    seg_pts.erase(it);
                    pt_it = find(seg_pts.cbegin(), seg_pts.cend(), current_pt);
                    pt_it = seg_pts.erase(pt_it);
                }
                else
                    ++pt_it;
            }
            // there should be two points left
            if (seg_pts.size() != 2)
                throw runtime_error("Invalid number of segment end points found: " + seg_pts.size());
            pt_it = seg_pts.begin();
#ifdef DEBUG_SIMPLIFY_MESH_2D
            cout << "Simplify_Mesh_2D::Pt_Remover::rem_perimeter_pts discovered end points #" << *pt_it << " and #" << *(pt_it + 1) << "\n";
#endif
            perimeter_segs.push_back(*pt_it, *(pt_it + 1), false);

#ifdef DEBUG_SIMPLIFY_MESH_2D
            cout << "Simplify_Mesh_2D::Pt_Remover::rem_perimeter_pts forming new facets\n";
#endif
            // remove perimeter point by forming new facets and replacing the original ones in same_plane_facets
            vector<Facet> new_facets;
            form_new_facets(orig_facet_2ds, perimeter_segs, new_facets);

            // add newly created facets back to the list of same plane facets
            for (vector<Facet>::const_iterator nf_it = new_facets.begin(); nf_it != new_facets.end(); ++nf_it)
            {
#ifdef DEBUG_SIMPLIFY_MESH_2D
                cout << "Simplify_Mesh_2D::Pt_Remover::rem_perimeter_pts adding new facet p1: " << nf_it->get_p1_index() << " p2: " << nf_it->get_p2_index() << " p3: " << nf_it->get_p3_index() << "\n";
#endif
                same_plane_facets.push_back(*nf_it);
            }
        }
    }

    Simplify_Mesh_2D::Simplify_Mesh_2D() {}
    
    const bool Simplify_Mesh_2D::operator()(Mesh_2D& mesh)
    {
        if (mesh.empty())
            return false;
        
#ifdef DEBUG_SIMPLIFY_MESH_2D
        cout << "Simplify_Mesh_2D::Operator() begin\n";
#endif
        // make a copy of mesh and do all operations on it
        Pt_Remover pt_remover(mesh);
#ifdef DEBUG_SIMPLIFY_MESH_2D
        int count = 0;
        for (Mesh_2D::const_point_iterator pt_iter = mesh.point_begin(); pt_iter != mesh.point_end(); ++pt_iter)
        {
            cout << "Simplify_Mesh_2D::operator() point #" << count++ << " x: " << (*pt_iter)->get_x() << " y: " << (*pt_iter)->get_y() << "\n";
        }
#endif        
        
#ifdef DEBUG_SIMPLIFY_MESH_2D
        cout << "Simplify_Mesh_2D::operator() removing internal points\n";
#endif
        // remove all internal points found
        pt_remover.rem_internal_pts();
        
#ifdef DEBUG_SIMPLIFY_MESH_2D
        cout << "Simplify_Mesh_2D::operator() removing perimeter points\n";
#endif
        pt_remover.rem_perimeter_pts();
        
        Mesh_2D temp_mesh(mesh.get_precision());
        for (Pt_Remover::const_iterator pr_it = pt_remover.begin(); pr_it != pt_remover.end(); ++pr_it)
        {
#ifdef DEBUG_SIMPLIFY_MESH_2D
            cout << "Simplify_Mesh_2D::operator() adding facet p1: " << pr_it->get_p1_index() << " p2: " << pr_it->get_p2_index() << " p3: " << pr_it->get_p3_index() << "\n";
#endif
            Mesh_2D::const_point_iterator pt_it(mesh.point_begin());
            advance(pt_it, pr_it->get_p1_index());
            shared_ptr<Point_2D> p1(*pt_it);
            pt_it = mesh.point_begin();
            advance(pt_it, pr_it->get_p2_index());
            shared_ptr<Point_2D> p2(*pt_it);
            pt_it = mesh.point_begin();
            advance(pt_it, pr_it->get_p3_index());
            shared_ptr<Point_2D> p3(*pt_it);
            temp_mesh.push_back(Facet_2D(p1, p2, p3));
        }
        
        if (temp_mesh.size() < mesh.size())
        {
            mesh.clear();
            for (Mesh_2D::const_iterator it = temp_mesh.begin(); it != temp_mesh.end(); ++it)
                mesh.push_back(*it);
            return true;
        }
        else
            return false;
    }

}
