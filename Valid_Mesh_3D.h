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
 * File:   Valid_Mesh_3D.h
 * Author: Jeffrey Davis
 */

#ifndef VALID_MESH_3D_H
#define VALID_MESH_3D_H

#include <vector>
#include <unordered_map>
#include <memory>
#include "Facet.h"
#include "Point_3D.h"
#include "Facet_3D.h"

using namespace std;

namespace VCAD_lib
{
    class Mesh_3D;
    
    /*
     * Looks for the following problems:
     * 1. Facets that contain other facets (can cause #3 below)
     * 2. Points in the middle of a facet side (causes mesh edges)
     * 3. facet sides belonging to more than two facets
     * 4. mesh edges (related to #3). For a valid 3D mesh, 
     *    there shouldn't be any unconnected sides
     */
    class Valid_Mesh_3D {
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
        struct Facet_3D_Hasher {
            const int operator()(const Facet_3D& facet) const;
        };
        
        struct Facet_3D_Predicate {
            const bool operator()(const Facet_3D& facet1, const Facet_3D& facet2) const;
        };
        
        class Facet_Find {
        public:
            Facet_Find(const Facet_3D& f);
            const bool operator()(const Facet_3D& f) const;
        private:
            const Facet_3D facet;
        };
        
        struct Facet_Sorter {
            const Vector_3D::Measurement area(const Facet_3D& f) const;
            const bool operator()(const Facet_3D& f1, const Facet_3D& f2) const;
        };
        
    public:
        typedef vector<shared_ptr<Point_3D>>::const_iterator pt_on_side_iterator;
        typedef vector<Facet_3D>::const_iterator edge_facet_iterator;
        typedef vector<vector<Facet_3D>>::const_iterator too_many_share_side_iterator;
        // facets inside facets iterator (gives pair<Facet_3D,vector<Facet_3D>>
        typedef unordered_map<Facet_3D,vector<Facet_3D>,Facet_3D_Hasher>::const_iterator facets_inside_facet_iterator;
        Valid_Mesh_3D(const Mesh_3D& mesh);
        // Maybe return an object that can be interpreted as a boolean and has all the data??
        const bool validate();
        pt_on_side_iterator pts_on_side_begin() const { return pts_on_facet_sides.begin(); }
        pt_on_side_iterator pts_on_side_end() const { return pts_on_facet_sides.end(); }
        const bool pts_on_side_empty() const { return pts_on_facet_sides.empty(); }
        vector<shared_ptr<Point_3D>>::size_type pts_on_side_size() { return pts_on_facet_sides.size(); }
        edge_facet_iterator edge_facets_begin() const { return edge_facets.begin(); }
        edge_facet_iterator edge_facets_end() const { return edge_facets.end(); }
        const bool edge_facets_empty() const { return edge_facets.empty(); }
        vector<Facet_3D>::size_type edge_facets_size() { return edge_facets.size(); }
        too_many_share_side_iterator too_many_share_side_begin() const { return too_many_share_side.begin(); }
        too_many_share_side_iterator too_many_share_side_end() const { return too_many_share_side.end(); }
        const bool too_many_share_side_empty() const { return too_many_share_side.empty(); }
        vector<vector<Facet_3D>>::size_type too_many_share_side_size() { return too_many_share_side.size(); }
        facets_inside_facet_iterator facets_inside_facet_begin() const { return facets_inside_facets.begin(); }
        facets_inside_facet_iterator facets_inside_facet_end() const { return facets_inside_facets.end(); }
        const bool facets_inside_facet_empty() const { return facets_inside_facets.empty(); }
        unordered_map<Facet_3D,vector<Facet_3D>,Facet_3D_Hasher,Facet_3D_Predicate>::size_type facets_inside_facet_size() { return facets_inside_facets.size(); }
    private:
        const Point_3D::Measurement precision;
        vector<shared_ptr<Point_3D>> all_points; // stores all points from mesh
        vector<shared_ptr<Point_3D>> pts_on_facet_sides; // points that are in the middle of a facet side
        vector<Facet> all_facets; // stores all facets from mesh
        vector<Facet_3D> edge_facets; // facets that form an edge
        vector<vector<Facet_3D>> too_many_share_side; // more than two facets that share the same side
        // key is the facet that contains the facets in the vector value
        unordered_map<Facet_3D,vector<Facet_3D>,Facet_3D_Hasher,Facet_3D_Predicate> facets_inside_facets;
        
        void find_same_plane_facets(vector<vector<Facet_3D>>& sp_facets) const;
    };


}

#endif /* VALIDATE_MESH_3D_H */

