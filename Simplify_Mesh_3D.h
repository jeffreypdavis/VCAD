/* 
 * File:   Simplify_Mesh_3D.h
 * Author: Jeffrey Davis
 *
 */

#ifndef SIMPLIFY_MESH_3D_H
#define SIMPLIFY_MESH_3D_H

#include <memory>
#include <vector>
#include "Point_3D.h"
#include "Facet.h"
#include "Facet_3D.h"
#include "Mesh_3D.h"

namespace VCAD_lib
{
    class Simplify_Mesh_3D {
    private: // classes to help simplify mesh
        /*
         * Can create new facets as long as new segments generated have a point half way
         * along segment that is inside the facets and does not intersect other internal
         * segments
         * 
         * A perimeter point can only be removed if all of the connecting facets
         * form two planes where the perimeter point is common to both planes.
         */
        class Pt_Remover {
        private:
            struct Segment {
                bool is_internal;
                bool used;
                int point1;
                int point2;
                Segment(const int p1, const int p2, const bool internal);
                const bool operator==(const Segment& seg) const;
                const int shares_pt(const Segment& seg) const;
            };

            class Segments {
            private:
                class Segment_Sort {
                public:
                    Segment_Sort();
                    const bool operator()(const Segment& seg1, const Segment& seg2) const;
                };
            public:
                typedef vector<Segment>::const_iterator const_iterator;
                Segments(const Mesh_3D& mesh);
                const_iterator begin() const { return segments.begin(); }
                const_iterator end() const { return segments.end(); }
                /*
                 * add a segment
                 */
                void push_back(const int p1, const int p2, const bool internal);
                /*
                 * sort segments to locate internal and perimeter points that might
                 * be able to be removed
                 */
                void find_pts_to_remove(vector<int>& internal_pts, vector<int>& perimeter_pts, 
                        vector<Facet>& same_plane_facets, const Mesh_3D& mesh) const;
                /*
                 * Get next segment
                 */
                const Segment* get_next_segment(const Segment* prev_segment) const;
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
                const Segment* find_connecting_segment(const Segment& segment, 
                        const Segment* prev_connecting_seg, int& shared_pt, const Point_3D::Measurement precision) const;
                /*
                 * is_segment_valid
                 */
                const bool is_seg_valid(const int p1, const int p2, const vector<Facet_3D>& orig_facets, 
                        const Mesh_3D::const_point_iterator pt_begin, const Mesh_3D::const_point_iterator pt_end, 
                        const Point_3D::Measurement precision) const;
                /*
                 * process segments that have been used
                 */
                void process_segs(const Segment& seg1, const Segment& seg2, const int seg3_p1, const int seg3_p2);
            private:
                const Mesh_3D::const_point_iterator pt_begin;
                const Mesh_3D::const_point_iterator pt_end;
                vector<Segment> segments;
                vector<Segment> removed_segs;
            };
        public:
            typedef vector<int>::const_iterator perimeter_pt_iter;
            typedef vector<Facet>::const_iterator const_iterator;
            Pt_Remover(const vector<Facet>& f_data, const Mesh_3D& mesh);
            const_iterator begin() const { return same_plane_facets.begin(); }
            const_iterator end() const { return same_plane_facets.end(); }
            perimeter_pt_iter pt_begin() const { return perimeter_pts.begin(); }
            perimeter_pt_iter pt_end() const { return perimeter_pts.end(); }
            void rem_internal_pts();
            void rem_perimeter_pt(const int pt);
        private:
            const Mesh_3D* const orig_mesh;
            vector<int> internal_pts; // internal points that can be removed
            vector<int> perimeter_pts; // perimeter points that may or may not be removed
            vector<Facet> same_plane_facets;
            /*
             * form new facets
             */
            void form_new_facets(const vector<Facet_3D>& orig_facets, Segments& perimeter_segs, 
                    vector<Facet>& new_facets) const;
        };

        /*
         * List of same plane facets within Mesh
         */
        class Facet_Datas {
        private:
            class Facet_Find {
            public:
                Facet_Find(const Facet& facet);
                const bool operator()(const Facet& facet);
            private:
                Facet facet_to_find;
            };
        public:
            typedef vector<Pt_Remover>::iterator iterator;
            typedef vector<Pt_Remover>::const_iterator const_iterator;
            typedef vector<Pt_Remover>::size_type size_type;
            /*
             * get the point associated with the index
             */
            Facet_Datas();
            const size_type size() { return facet_datas.size(); }
            iterator begin() { return facet_datas.begin(); }
            iterator end() { return facet_datas.end(); }
            const_iterator cbegin() { return facet_datas.cbegin(); }
            const_iterator cend() { return facet_datas.cend(); }
            /*
             * sort facets into planes and locate internal and perimeter points to remove
             */
            void process_mesh(const Mesh_3D& mesh);
        private:
            vector<Pt_Remover> facet_datas;
            /*
             * find next same plane facets. 
             */
            void get_next_sp_facets(vector<Facet>& same_plane_facets, 
                    vector<Facet>& facets, const Mesh_3D& mesh) const;
        };
        
    public:
        // exception safety: strong guarantee
        Simplify_Mesh_3D();
        /*
         * Try to simplify mesh by removing unnecessary points
         */
        // exception safety: strong guarantee
        const bool operator()(Mesh_3D& mesh);
    };

}

#endif /* SIMPLIFY_MESH3_3D_H */

