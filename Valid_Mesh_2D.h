/* 
 * File:   Valid_Mesh_2D.h
 * Author: Jeffrey Davis
 *
 * Created on May 20, 2016, 12:42 PM
 */

#ifndef VALID_MESH_2D_H
#define VALID_MESH_2D_H

#include <vector>
#include <unordered_map>
#include <memory>
#include "Facet.h"
#include "Point_2D.h"
#include "Facet_2D.h"

namespace VCAD_lib
{
    class Mesh_2D;
    
    /*
     * Looks for the following problems:
     * 1. Facets that contain other facets
     * 2. Points in the middle of a facet side
     * 3. more than two facets that share the same side
     */
    class Valid_Mesh_2D {
    private:
        struct Facet_Side {
            int point1;
            int point2;
            Facet_Side(const int p1, const int p2);
            const bool operator==(const Facet_Side& side) const;
        };
        
        // Facet_Side hasher
        struct FS_Hasher {
            const int operator()(const Facet_Side& side) const;
        };
        
        // Facet hasher
        struct Facet_2D_Hasher {
            const int operator()(const Facet_2D& facet) const;
        };
        
        struct Facet_2D_Predicate {
            const bool operator()(const Facet_2D& facet1, const Facet_2D& facet2) const;
        };
        
        struct Facet_Sorter {
            const Vector_2D::Measurement area(const Facet_2D& f) const;
            const bool operator()(const Facet_2D& f1, const Facet_2D& f2) const;
        };
        
    public:
        typedef vector<shared_ptr<Point_2D>>::const_iterator pt_on_side_iterator;
        typedef vector<vector<Facet_2D>>::const_iterator too_many_share_side_iterator;
        // facets inside facets iterator (gives pair<Facet_2D,vector<Facet_2D>>
        typedef unordered_map<Facet_2D,vector<Facet_2D>,Facet_2D_Hasher,Facet_2D_Predicate>::const_iterator facets_inside_facet_iterator;
        Valid_Mesh_2D(const Mesh_2D& mesh);
        const bool validate();
        pt_on_side_iterator pts_on_side_begin() const { return pts_on_facet_sides.begin(); }
        pt_on_side_iterator pts_on_side_end() const { return pts_on_facet_sides.end(); }
        const bool pts_on_side_empty() const { return pts_on_facet_sides.empty(); }
        vector<shared_ptr<Point_2D>>::size_type pts_on_side_size() { return pts_on_facet_sides.size(); }
        too_many_share_side_iterator too_many_share_side_begin() const { return too_many_share_side.begin(); }
        too_many_share_side_iterator too_many_share_side_end() const { return too_many_share_side.end(); }
        const bool too_many_share_side_empty() const { return too_many_share_side.empty(); }
        vector<vector<Facet_2D>>::size_type too_many_share_side_size() { return too_many_share_side.size(); }
        facets_inside_facet_iterator facets_inside_facet_begin() const { return facets_inside_facets.begin(); }
        facets_inside_facet_iterator facets_inside_facet_end() const { return facets_inside_facets.end(); }
        const bool facets_inside_facet_empty() const { return facets_inside_facets.empty(); }
        unordered_map<Facet_2D,vector<Facet_2D>,Facet_2D_Hasher,Facet_2D_Predicate>::size_type facets_inside_facet_size() { return facets_inside_facets.size(); }
    private:
        const Point_2D::Measurement precision;
        vector<shared_ptr<Point_2D>> all_points; // stores all points from mesh
        vector<shared_ptr<Point_2D>> pts_on_facet_sides; // points that are in the middle of a facet side
        vector<Facet> all_facets; // stores all facets from mesh
        vector<vector<Facet_2D>> too_many_share_side; // more than two facets that share the same side
        // key is the facet that contains the facets in the vector value
        unordered_map<Facet_2D,vector<Facet_2D>,Facet_2D_Hasher,Facet_2D_Predicate> facets_inside_facets; 
    };


}

#endif /* VALIDATE_MESH_2D_H */

