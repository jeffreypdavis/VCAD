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
 * File:   Intersect_Meshes_3D.h
 * Author: Jeffrey Davis
 */

#ifndef INTERSECT_MESHES_3D_H
#define INTERSECT_MESHES_3D_H

#include <forward_list>
#include <memory>
#include "Point_3D.h"
#include "Facet.h"
#include "Facet_3D.h"
#include "Mesh_3D.h"

using namespace std;

namespace VCAD_lib
{

    class Intersect_Meshes_3D {
    private:
        /*
         * Intersect Point location
         * internal: inside facet
         * p1p2: on side from corner point 1 to corner point 2
         * p1p3: on side from corner point 1 to corner point 3
         * p2p3: on side from corner point 2 to corner point 3
         * p1: corner point 1
         * p2: corner point 2
         * p3: corner point 3
         */
        enum Location { internal, p1p2, p1p3, p2p3, p1, p2, p3 };
        
        /*
         * An intersection point common to two facets
         */
        struct Intersect_Point {
            shared_ptr<Point_3D> pt; // the intersection point
            Location f1_loc;         // the location of the intersection point on facet1
            Location f2_loc;         // the location of the intersection point on facet2
            Intersect_Point();       // null initialization
            /*
             * Constructor
             * 
             * Arguments:
             * i_point: the intersect point
             * f1_location: the intersect point location on facet1
             * f2_location: the intersect point location on facet2
             */
            Intersect_Point(const shared_ptr<Point_3D> i_point, const Location f1_location, 
                    const Location f2_location);
            /*
             * equality operator.  returns true if the points match as well as the locations.
             */
            const bool operator==(const Intersect_Point& ipt) const;
        };

        /*
         * A list of intersect points
         */
        class I_Pt_List {
        public:
            typedef vector<Intersect_Point>::size_type size_type;
            typedef vector<Intersect_Point>::iterator iterator;
            typedef vector<Intersect_Point>::const_iterator const_iterator;
            
            I_Pt_List();
            const bool empty() const { return i_points.empty(); }
            const size_type size() const { return i_points.size(); }
            void clear() { i_points.clear(); }
            /*
             * Intersect_Point iterators
             */
            iterator begin() { return i_points.begin(); }
            iterator end() { return i_points.end(); }
            const_iterator begin() const { return i_points.begin(); }
            const_iterator end() const { return i_points.end(); }
            const_iterator cbegin() const { return i_points.begin(); }
            const_iterator cend() const { return i_points.end(); }
            // add an intersect point
            void push_back(const Intersect_Point& ip);
            // erase an intersect point
            void erase(const_iterator it) { i_points.erase(it); }
            /*
             * Verify intersect points and try to correct any problems found
             * 
             * Arguments:
             * facet1: the facet1 used in the intersection
             * facet2: the facet2 used in the intersection
             */
            void validate(const Facet_3D& facet1, const Facet_3D& facet2);
        private:
            vector<Intersect_Point> i_points; // vector to hold intersect points
            /*
             * returns true if the points match, but ignores f1_loc and f2_loc
             */
            const bool matches(const Point_3D& p1, const Point_3D& p2) const;
            /*
             * If two intersect points have the same points, but different locations
             * then this method is called (from validate).  Determines which point 
             * to keep and which to remove.  In one case, both are points are
             * removed and a corner point is added instead.
             * 
             * Arguments:
             * matching_ip1: has the same intersect point as matching_ip2
             * matching_ip2: has the same intersect point as matching_ip1
             * facet1: the facet1 used in the intersection
             * facet2: the facet2 used in the intersection
             */
            const_iterator process_matching_i_pts(const Intersect_Point& matching_ip1, 
                    const const_iterator matching_ip2, const Facet_3D& facet1, 
                    const Facet_3D& facet2);
        };
        
        /*
         * Determines the intersect points common to two facets.
         * 
         * 1. Intersect each facet side of facet1 with each facet side of facet2
         * 2. If no intersection of facet1 side, check if the side intersects facet2
         * 3. if there is one intersection of facet1 side, check if either corner
         *    is in facet2.
         * 4. If a facet2 side was not intersected, check if it intersects facet1
         * 5. If a facet2 side has one intersection, check if either corner is inside facet1
         */
        class I_Pt_Locator {
        private:
            /*
             * An intersect point holder to facilitate processing
             */
            struct I_Pt_Data {
                int num;
                Intersect_Point ip1;
                Intersect_Point ip2;
                I_Pt_Data(); // sets num to zero and initializes ip1 and ip2 with null constructor
            };
            
            class Point_Find {
            public:
                Point_Find(const shared_ptr<Point_3D>& point, const Point_3D::Measurement prec);
                const bool operator()(const shared_ptr<Point_3D>& other_pt);
            private:
                const Point_3D::Measurement precision;
                const shared_ptr<Point_3D> pt;
            };
        public:
            /*
             * Constructor
             * 
             * Arguments:
             * f1: facet1 to intersect into facet2
             * f2: facet2 to intersect into facet1
             * prec: the precision to intersect the two facets
             */
            I_Pt_Locator(const Point_3D::Measurement prec);
            /*
             * find all intersect points for the two facets
             * 
             * Arguments:
             * intersect_points: filled with the intersect points found
             */
            const bool operator()(const Facet_3D& f1, const Facet_3D& f2, I_Pt_List& intersect_points);
        private:
            const Point_3D::Measurement precision;
            Facet_3D facet1;
            Facet_3D facet2;
            vector<shared_ptr<Point_3D>> generated_pts;
            
            const bool matches(const shared_ptr<Point_3D>& p1, const shared_ptr<Point_3D>& p2) const;
            
            /*
             * determines the location of the intersect point. On the side,
             * or is either corner.  If the location is a corner, i_pt is updated to
             * use the corner point location.
             * 
             * Arguments:
             * side_start: the facet side start point (p1 or p2)
             * side_end: the facet side end point (p2 or p3)
             * side_loc: the location of the side (p1p2, p1p3, or p2p3)
             * i_pt: the intersect point to check
             * loc: the determined i_pt location
             */
            void side_i_pt_loc(const shared_ptr<Point_3D> side_start, 
                    const shared_ptr<Point_3D> side_end, const Location side_loc, 
                    shared_ptr<Point_3D>& i_pt, Location& loc);

            /*
             * Intersects f1 side with the three sides of facet2.
             * Updates i_pt_data to have a maximum of two intersect points if
             * any are found.  Returns true if one or two intersect points were
             * found.
             * 
             * Arguments:
             * f1_side_start: the starting point of the f1 side (p1 or p2)
             * f1_Side_end: the ending point of the f1 side (p2 or p3)
             * f1_side: the location of the f1_side (p1p2, p1p3, or p2p3)
             * f2_side_start: the starting point of the f2 side (p1 or p2)
             * f2_side_end: the ending point of the f1_side (p2 or p3)
             * f2_side: the location of the f2 side (p1p2, p1p3, or p2p3)
             * i_pt_data: the intersect points found
             */
            const bool intersect_sides(const shared_ptr<Point_3D> f1_side_start, 
                    const shared_ptr<Point_3D> f1_side_end, const Location f1_side, 
                    const shared_ptr<Point_3D> f2_side_start, 
                    const shared_ptr<Point_3D> f2_side_end, const Location f2_side, 
                    I_Pt_Data& i_pt_data);

            /*
             * Intersects a vector with the three sides of a facet.  Fills p1 and p2
             * with intersect points and sets p1p2_intersected, p1p3_intersected, and/or 
             * p2p3_intersected to true if those facet sides were intersected by the 
             * vector.  Returns the number of intersect points found (0, 1, or 2)
             * 
             * Arguments:
             * f1_side: the side of facet1 to intersect into facet2
             * intersect_points: the list of intersect points to fill
             */
            void intersect_f1_side_to_f2(const Location f1_side, I_Pt_List& intersect_points);

            /*
             * Intersects a facet side into the intersecting_facet. Adds any internal
             * segments, updates side_points and single_i_points as needed.
             * 
             * Arguments:
             * f2_side: the side of facet2 to intersect into facet1
             * intersect_points: the list of intersect points to fill
             */
            void intersect_f2_side_to_f1(const Location f2_side, I_Pt_List& intersect_points);

            /*
             * Determine if a facet was intersected.  Returns the number of 
             * intersections (0, 1 or 2)
             * 
             * Arguments:
             * f2_side: the side of facet2 to test if it was intersected
             * intersect_points: the list of intersect points
             * corner1_found: set to true or false depending if the side start
             *                corner is found as an intersect point
             * corner2_found: set to true or false depending if the side end
             *                corner is found as an intersect point
             */
            const int is_f2_side_intersected(const Location f2_side, 
                const I_Pt_List& intersect_points, bool& corner1_found,
                bool& corner2_found);
        };
        
        /*
         * List of facets with functions to facilitate intersecting facets and 
         * maintain a consistent mesh
         */
        class Facets {
        public:
            typedef vector<Facet>::size_type size_type;
            typedef vector<Facet>::const_iterator const_iterator;
            typedef vector<shared_ptr<Point_3D>>::const_iterator pt_const_iterator;
            
            Facets();
            Facets(const Mesh_3D& mesh);
            const bool empty() const { return facet_list.empty(); }
            size_type size() const { return facet_list.size(); }
            const_iterator begin() const { return facet_list.begin(); }
            const_iterator cbegin() const { return facet_list.begin(); } 
            const_iterator end() const { return facet_list.end(); }
            const_iterator cend() const { return facet_list.end(); }
            pt_const_iterator pts_cbegin() const { return point_list.begin(); } 
            pt_const_iterator pts_cend() const { return point_list.end(); }
            void clear();
            /*
             * get the point associated with the index
             */
            const shared_ptr<Point_3D> get_point(int index) const;
            /*
             * Add a Facet
             */
            void push_back(const Facet_3D& facet);
            /*
             * Get the last facet in the list and remove it
             */
            const Facet_3D pop();
            /*
             * Clear the list and add new facets
             */
            void replace_all(const Facets& facets);
            /*
             * Checks if the facet list contains the facet
             */
            const bool contains(const Facet_3D& facet) const;
            /*
             * Finds the facet in the list.  Throws runtime_error if it cannot 
             * find the facet.
             */
            const Facet find_facet(const Facet_3D& facet) const;
            /*
             * Replaces facet with new_facets
             * returns a pointer to the first of the newly inserted facets
             * 
             * Arguments:
             * facet: the facet to replace
             * new_facets: the facets to replace facet with
             */
            const_iterator replace_facet(const Facet facet, const Facets& new_facets);
        private:
            vector<shared_ptr<Point_3D>> point_list; // list of points
            vector<Facet> facet_list; // list of Facet
            /*
             * Checks if the point values are the same
             */
            const bool matches(const Point_3D& p1, const Point_3D& p2) const;
        };

	/*
         * A Class to fracture a facet into new facets based on the intersect points
         */
        class Facet_Builder {
        private:
            /*
             * Facet is broken up into line segments.  The segments are used to form
             * new Facets.
             */
            struct Line_Segment {
                bool used; // internal line segments are used twice, external line segments only once
                Location location; // the location of the segment (internal or one of p1p2, p1p3, or p2p3)
                shared_ptr<Point_3D> point1; // the beginning of the line segment
                shared_ptr<Point_3D> point2; // the end of the line segment
                /*
                 * Constructor
                 *
                 * Arguments:
                 * p1: the start of the line segment
                 * p2: the end of the line segment
                 * location: the location of the line segment (internal or one of the sides p1p2, p1p3, or p2p3)
                 */
                Line_Segment(const shared_ptr<Point_3D>& p1, const shared_ptr<Point_3D>& p2, const Location location);
                /*
	         * determines if the line segments are equal (share same points.  Location is not considered)
                 */
                const bool operator==(const Line_Segment& seg) const;
                /*
                 * checks if the other line segment shares a point with this segment
                 * 
                 * Arguments:
                 * seg: the segment to check
                 * shared_pt: function sets this to the shared point if one is found
                 */
                const bool shares_pt(const Line_Segment& seg, shared_ptr<Point_3D>& shared_pt) const;
            };
            
            /*
             * A grouping of segment objects.  Ensures that any segment added does 
             * not intersect an already existing segment
             */
            class Segments {
                class Segment_Sort {
                public:
                    Segment_Sort();
                    const bool operator()(const Line_Segment& seg1, const Line_Segment& seg2) const;
                };
            public:
                typedef vector<Line_Segment>::size_type size_type;
                typedef vector<Line_Segment>::const_iterator const_iterator;
            
                Segments();
                const_iterator begin() const { return segments.begin(); }
                const_iterator end() const { return segments.end(); }
                const size_type size() const; // overall size of internal and external segments
                void clear(); // clear segments
                /*
                 * determines if the segment already exists in either external or internal segments
                 * 
                 * Arguments:
                 * p1: the start of the segment
                 * p2: the end of the segment
                 * 
                 * returns true if the segment was found, false otherwise
                 */
                const bool contains_segment(const shared_ptr<Point_3D>& p1, const shared_ptr<Point_3D>& p2);
                // add an internal segment
                void add_internal_segment(const shared_ptr<Point_3D>& p1, const shared_ptr<Point_3D>& p2);
                // add an external segment. side_loc specifies the side of the segment on the original facet
                void add_external_segment(const shared_ptr<Point_3D>& p1, const shared_ptr<Point_3D>& p2, const Location side_loc);
                /*
                 * get the next segment
                 * 
                 * Arguments:
                 * prev_Segment: the previous segment found.  Can be zero to indicate the first time used
                 * 
                 * returns a pointer to the next segment or zero if no segments are found
                 */
                const Line_Segment* get_next_segment(const Line_Segment* prev_segment) const;
                /*
                 * find the line segment and return a pointer to it
                 * 
                 * Arguments:
                 * p1: the start of the segment
                 * p2: the end of the segment
                 * 
                 * returns a pointer to the segment or throws runtime_error if 
                 * it could not find the segment
                 */
                const Line_Segment* find_segment(const shared_ptr<Point_3D>& p1, const shared_ptr<Point_3D>& p2) const;
                /*
                 * find the next connecting segment.  A connecting segment is one
                 * which shares a point with the segment
                 * 
                 * Arguments:
                 * segment: the segment to find a connecting segment for
                 * prev_connecting_seg: the previous connecting segment found.  Can be zero
                 * shared_pt: if a connecting segment is found, function sets this to the 
                 *            common point found between the two segments
                 * precision: the precision to find the connecting segment
                 * 
                 * returns a pointer to a connecting segment or zero if none were
                 * found.
                 */
                const Line_Segment* find_connecting_segment(const Line_Segment& segment, 
                        const Line_Segment* prev_connecting_seg, shared_ptr<Point_3D>& shared_pt, 
                        const Point_3D::Measurement precision) const;
                /* 
                 * determine if a segment is an external or internal segment
                 * 
                 * Arguments:
                 * seg3: the segment to test if it is an internal or external segment
                 * side_loc: if external, function sets this to the facet side it is on
                 * original_facet: the original facet being fractured
                 * p1p2_pts: p1p2 side intersect points
                 * p1p3_pts: p1p3 side intersect points
                 * p2p3_pts: p2p3 side intersect points
                 * internal_pts: the internal points found inside the facet
                 */
                const bool is_segment_external(const Line_Segment& seg3, Location& side_loc, 
                        const Facet_3D& original_facet, const vector<shared_ptr<Point_3D>>& p1p2_pts, 
                        const vector<shared_ptr<Point_3D>>& p1p3_pts, const vector<shared_ptr<Point_3D>>& p2p3_pts, 
                        const vector<shared_ptr<Point_3D>>& internal_pts) const;
                /* 
                 * check if the formed segment3 crosses other segments
                 * 
                 * Arguments:
                 * seg3: the segment to check
                 * orig_facet: the original facet that is being fractured
                 * p1p2_pts: p1p2 side intersect points
                 * p1p3_pts: p1p3 side intersect points
                 * p2p3_pts: p2p3 side intersect points
                 * internal_points: the original facet internal intersect points
                 * precision: the precision to test for intersects with
                 */
                const bool does_seg_intersect(const Line_Segment& seg3, const Facet_3D& orig_facet, 
                        const vector<shared_ptr<Point_3D>>& p1p2_pts, const vector<shared_ptr<Point_3D>>& p1p3_pts, 
                        const vector<shared_ptr<Point_3D>>& p2p3_pts, const vector<shared_ptr<Point_3D>>& internal_points, 
                        const Point_3D::Measurement precision) const;
                /* 
                 * updates used property of segment or removes them from the list of
                 * segments
                 * 
                 * Arguments:
                 * seg1: segment1 used to form a new facet
                 * seg2: segment2 used to form a new facet
                 */
                void process_used_segments(const Line_Segment& seg1, const Line_Segment& seg2);
            private:
                vector<Line_Segment> segments;
                vector<Line_Segment> removed_segs;
                /*
                 * processes segment
                 */
                void process_used_segment(vector<Line_Segment>::iterator seg);
            };
            
            /*
             * Intersect points located on the sides of a facet
             */
            struct Intersecting_Facet_Side_Pts {
                I_Pt_List p1p2_points;
                I_Pt_List p1p3_points;
                I_Pt_List p2p3_points;
                Intersecting_Facet_Side_Pts();
                void clear() { p1p2_points.clear(); p1p3_points.clear(); p2p3_points.clear(); }
            };
            
            /*
             * Find a segment that contains point
             */
            class Segment_Find {
            public:
                Segment_Find(const shared_ptr<Point_3D>& point);
                const bool operator()(const Intersect_Meshes_3D::Facet_Builder::Line_Segment& segment);
            private:
                const shared_ptr<Point_3D> pt;
            };
        public:
            /*
             * Constructor
             * 
             * Arguments:
             * for_f1: for facet1
             * orig_Facet: the original facet that will be fractured
             * prec: the precision to form new facets to
             */
            Facet_Builder(const bool for_f1, const Facet_3D& orig_facet, const Point_3D::Measurement prec);

            const Facet_3D& get_facet() const { return orig_facet; }
            
            /*
             * Add a new intersection.  Generates internal segments based on
             * intersection points.
             * 
             * Arguments
             * intersect_pts: all of the intersect points
             * side_points: the side intersect points for this facet
             * i_side_points: the side intersect points for the other facet
             * internal_points: the internal points found in the intersection
             */
            void add_intersection(I_Pt_List& intersect_pts);
            
            /*
             * creates new facets.  If no new facets are generated, new_facets
             * will be empty and function will return false.  If new facets were
             * generated, function will return true and the new facets will be in 
             * new_facets.
             * 
             * Arguments
             * new_facets: Facets list to store new facets in.
             */
            const bool form_new_facets(Facets& new_facets);
        private:
            const bool for_facet1;
            const Point_3D::Measurement precision;
            const Facet_3D orig_facet;
            vector<shared_ptr<Point_3D>> internal_pts;
            vector<shared_ptr<Point_3D>> p1p2_pts;
            vector<shared_ptr<Point_3D>> p1p3_pts;
            vector<shared_ptr<Point_3D>> p2p3_pts;
            Segments segments;
            
            /*
             * Generate an internal segment for facet based on a list of two 
             * intersect points. For example, if an intersecting facet side 
             * intersects orig_facet twice, then this would generate a segment
             * between those two points.
             * 
             * Arguments:
             * i_pts: a list of intersect points that has two points.
             */
            void process_two_i_pts(const I_Pt_List& i_pts);
            
            /*
             * Processes a list of intersecting facet side points.  If there is
             * only one point, it does nothing, if there are two, it calls
             * process_two_i_pts, and if there are more than two, it throws a
             * runtime_error because there should only be a maximum of two
             * side intersect points.  If process_two_pts was called, then the 
             * intersect points used are removed from t_intersect_pts
             * 
             * Arguments:
             * i_side_pts: the side intersect points of the intersecting facet side
             * t_intersect_pts: the temporary intersect point list
             */
            void process_i_side_pts(const I_Pt_List& i_side_pts, I_Pt_List& t_intersect_pts);
            
            /*
             * Generates internal segments based on the intersect points.
             * 
             * Arguments:
             * intersect_points: the total intersect points
             * intersecting_side_pts: the side intersect points
             * internal_pts: the internal intersect points
             */
            void gen_internal_segs(const I_Pt_List& intersect_points, 
                    const Intersecting_Facet_Side_Pts& intersecting_side_pts);
            
            /*
             * Create a link segment
             */
            const Line_Segment* create_link_segment(const shared_ptr<Point_3D> pt);
            
            /*
             * Look for any internal points that are not part of an internal
             * segment.  Creates a link segment for these or throws a runtime_error
             * if no link segment could be created.
             */
            void check_internal_pts();
            
            /*
             * Check if segment is connected to a facet corner or side point.
             * return true if it is, false if both end points are internal points
             */
            const bool check_internal_seg(const Line_Segment& seg);
            
            const bool find_path(const Segments::const_iterator& current_it, 
                    const shared_ptr<Point_3D>& point, vector<Line_Segment>& path, 
                    vector<Line_Segment>& checked_segs);
            
            const bool complete_path(vector<Line_Segment>& path, 
                    vector<Line_Segment>& checked_segs, const Line_Segment& segment);
            
            /*
             * Check for any internal points without segments or any internal
             * segments not connected to a corner or side point.
             */
            void verify_internal_paths();
            
            /*
             * Check for internal segments that are in a straight line without
             * another segment connecting to the shared point
             */
            void verify_internal_segs();
            
            /*
             * form the outside segments of the facet using side_point data. 
             * p1 and p2 are the corner points of the facet side
             * 
             * Arguments:
             * side_points: the side intersect points found
             * side: the side to form facets for (p1p2, p1p3, or p2p3)
             * p1: the side corner point (p1 or p2)
             * p2: the side end point (p2 or p3)
             */
            void gen_side_segs(vector<shared_ptr<Point_3D>>& side_points, 
                    const Location side, const shared_ptr<Point_3D>& p1, 
                    const shared_ptr<Point_3D>& p2);
        
            /*
             * Determine if a facet contains an internal point
             * 
             * Arguments:
             * p1: point1 of the facet
             * p2: point2 of the facet
             * p3: point3 of the facet
             * facet: the facet made up of p1, p2, and p3
             */
            const bool contains_internal_pt(const shared_ptr<Point_3D>& p1, 
                    const shared_ptr<Point_3D>& p2, const shared_ptr<Point_3D>& p3, 
                    const Facet_3D& facet);
            
            /*
             * Creates new facets, but keeps overlapping facets the same as other_new_facets
             * 
             * Arguments:
             * side_pts: the side intersect points found
             * internal_points: the internal intersect points found
             * new_facets: the function will put newly generated facets in this list
             */
            void build_facets(Facets& new_facets);
            
#ifdef DEBUG_INTERSECT_MESHES_3D_FACET_BUILDER
            const Point_3D::Measurement facet_area(const Facet_3D& facet);
#endif
        };
        
        /*
         * Sort intersected facets to know which facets are on the surface or 
         * inside the other mesh.
         */
        class Facet_Sorter {
        public:
            typedef vector<Facet>::const_iterator const_iterator;
            
            const_iterator f1_surface_begin() { return f1_on_surface_f2.begin(); }
            const_iterator f1_surface_end() { return f1_on_surface_f2.end(); }
            const_iterator f1_inside_begin() { return f1_inside_f2.begin(); }
            const_iterator f1_inside_end() { return f1_inside_f2.end(); }
            const_iterator f2_surface_begin() { return f2_on_surface_f1.begin(); }
            const_iterator f2_surface_end() { return f2_on_surface_f1.end(); }
            const_iterator f2_inside_begin() { return f2_inside_f1.begin(); }
            const_iterator f2_inside_end() { return f2_inside_f1.end(); }
            
            Facet_Sorter(const Point_3D::Measurement& prec);
            /*
             * determine the location of the facets in facets1 and facets2 in 
             * relation to each other.  each facet is determined if it is on the 
             * surface of or inside of the other mesh.
             * 
             * Arguments:
             * facets1: a mesh that was intersected by facets2
             * facets2: a mesh that was intersected by facets1
             */
            void sort(const Facets& facets1, const Facets& facets2);
            void clear();
        private:
            Point_3D::Measurement precision;
            vector<Facet> f1_on_surface_f2;
            vector<Facet> f1_inside_f2;
            vector<Facet> f2_on_surface_f1;
            vector<Facet> f2_inside_f1;
        };
    public:
        /*
         * Constructor.  Does not perform any actions because class does not 
         * contain any objects to initialize.
         */
        Intersect_Meshes_3D();
        /*
         * Intersect mesh1 into mesh2.  returns true if new facets were generated
         * because of the intersection.  mesh1_result and mesh2_result are only
         * updated if function returns true;
         * 
         * Arguments:
         * mesh1: mesh to intersect into mesh2
         * mesh2: mesh to intersect into mesh1
         * mesh1_result: the intersected mesh1 result fractured with facets aligned to mesh2_result
         * mesh2_result: the intersected mesh2 result fractured with facets aligned to mesh1_result
         */
        const bool operator()(const Mesh_3D& mesh1, const Mesh_3D& mesh2, Mesh_3D& mesh1_result, Mesh_3D& mesh2_result);
        /*
         * subtract mesh2 from mesh1 and store in result
         * 
         * Arguments:
         * mesh1: mesh to subtract mesh2 from
         * mesh2: mesh to subtract from mesh1
         * result: the result of mesh1 - mesh2
         */
        void difference(const Mesh_3D& mesh1, const Mesh_3D& mesh2, Mesh_3D& result);
        /* 
         * intersect mesh1 and mesh2 and keep only the facets that are in both.
         * 
         * Arguments:
         * mesh1: mesh to intersect into mesh2
         * mesh2: mesh to intersect into mesh1
         * result: the intersected mesh1 and mesh2 result
         */
        void intersection(const Mesh_3D& mesh1, const Mesh_3D& mesh2, Mesh_3D& result);
        /*
         * merge mesh1 and mesh2 together and store in result.
         * 
         * Arguments:
         * mesh1: mesh to combine into mesh2
         * mesh2: mesh to combine into mesh1
         * result: the combined mesh1 and mesh2
         */
        void merge(const Mesh_3D& mesh1, const Mesh_3D& mesh2, Mesh_3D& result);
    private:
        
        /*
         * Intersect two facets.  
         * 
         * Arguments:
         * facet1: a facet to intersect into facet2
         * facet2: a facet to intersect into facet1
         * precision: the precision to perform the intersection
         */
        void intersect_facets(Facets& facets1, Facets& facets2, const Point_3D::Measurement precision);
    };

}

#endif /* INTERSECT_MESHES5_3D_H */

