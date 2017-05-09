/* 
 * File:   Mesh_3D.h
 * Author: jeffrey Davis
 *
 */

#ifndef MESH_3D_H
#define MESH_3D_H

#include <vector>
#include <memory>
#include <iterator>
#include "Point_3D.h"
#include "Facet.h"
#include "Facet_3D.h"

namespace VCAD_lib
{
    class Mesh_3D {
    private:
        class Facet_Find {
        public:
            Facet_Find(const Facet& facet);
            bool operator()(const Facet& facet);
        private:
            const Facet facet_to_find;
        };
    public:
        class const_iterator {
        public:
            typedef bidirectional_iterator_tag iterator_category;

            const_iterator(const vector<shared_ptr<Point_3D>>::const_iterator point_it_begin, 
                    const vector<Facet>::const_iterator facet_it_begin, 
                    const vector<Facet>::const_iterator facet_it_end, 
                    const vector<Facet>::const_iterator position);

            bool operator==(const const_iterator&) const;
            bool operator!=(const const_iterator&) const;

            const_iterator& operator++();
            const_iterator operator++(int);
            const_iterator& operator--();
            const_iterator operator--(int);

            const Facet_3D& operator*() const;
            const Facet_3D* const operator->() const;
        private:
            vector<shared_ptr<Point_3D>>::const_iterator point_list_begin;
            vector<Facet>::const_iterator facets_begin;
            vector<Facet>::const_iterator facets_end;
            vector<Facet>::const_iterator current_facet;
            Facet_3D facet;
            void update_facet();
        };        
        
        typedef vector<Facet>::size_type size_type;
        typedef Point_3D::Measurement Measurement;
        typedef Point_3D::Angle_Meas Angle_Meas;
        // to allow for users to get point list and facet point indices
        typedef vector<shared_ptr<Point_3D>>::const_iterator const_point_iterator;
        typedef vector<Facet>::const_iterator const_facet_iterator;
        
        // exception safety: strong guarantee
        Mesh_3D();
        // exception safety: strong guarantee - invalid_argument if precision is less than or equal to zero
        explicit Mesh_3D(const Measurement precision); // precision
        // exception safety: strong guarantee
        Mesh_3D(const Mesh_3D& orig);
        // exception safety: strong guarantee
        Mesh_3D& operator=(const Mesh_3D&);
        // exception safety: no throw
        Measurement get_precision() const { return precision; }
        // iterators
        const_iterator begin() const;
        const_iterator cbegin() const; 
        const_iterator end() const;
        const_iterator cend() const;
        
        const_point_iterator point_begin() const { return point_list.cbegin(); }
        const_point_iterator point_end() const { return point_list.cend(); }
        const_facet_iterator facet_begin() const { return facet_list.cbegin(); }
        const_facet_iterator facet_end() const { return facet_list.cend(); }
        
        void push_back(const Facet_3D& facet);
        size_type size() const { return facet_list.size(); }
        bool empty() const { return facet_list.empty(); }
        void clear();
        const_iterator erase(const_iterator it);
        const_iterator erase(const_iterator begin, const_iterator end);
//        const_iterator insert(const_iterator loc, const Facet_3D& facet);
//        const_iterator insert(const_iterator loc, const_iterator from_facet, const_iterator to_facet);
        // rotate
        // exception safety: basic guarantee
        Mesh_3D& rotate(const Angle& angle);
        // exception safety: basic guarantee
        Mesh_3D& rotate(const Angle_Meas angle, const Vector_3D& axis);
        // exception safety: basic guarantee
        Mesh_3D& rotate(const Angle& angle, const Point_3D& origin);
        // exception safety: basic guarantee
        Mesh_3D& rotate(const Angle_Meas angle, const Vector_3D& axis, const Point_3D& origin);
        // scale
        // exception safety: basic guarantee
        Mesh_3D& scale(const Measurement x_scalar, const Measurement y_scalar, 
                const Measurement z_scalar);
        // exception safety: basic guarantee
        Mesh_3D& scale(const Measurement x_scalar, const Measurement y_scalar, 
                const Measurement z_scalar, const Point_3D& origin);
        // translate
        // exception safety: no throw
        Mesh_3D& translate(const Measurement x_val, const Measurement y_val, 
                const Measurement z_val);
        // exception safety: no throw
        Mesh_3D& translate(const Vector_3D&);
        /*
         * move the facet to a different coordinate system and origin.
         * A new coordinate system must be supplied and this is done by 
         * specifying a new origin, a new axis, and one more point that is
         * not on the axis to define a plane.  ref_origin is an option to 
         * move based on a different origin than 0,0,0.
         * 
         * There are six methods to allow for specifying the new coordinate system
         * the names of the methods start with move.  Then next letter is the 
         * axis specified in the vector argument - x, y, or z axis.  The last part of the name
         * is identifying where the third point is located - xy plane, xz plane,
         * or yz plane.
         * 
         * exception safety: basic guarantee
         */
        Mesh_3D& move_x_pxy(const Point_3D& new_origin, const Vector_3D& x_axis, 
                const Point_3D& pt_xy_plane, const Point_3D& ref_origin=Point_3D(0,0,0));
        Mesh_3D& move_x_pxz(const Point_3D& new_origin, const Vector_3D& x_axis, 
                const Point_3D& pt_xz_plane, const Point_3D& ref_origin=Point_3D(0,0,0));
        Mesh_3D& move_y_pxy(const Point_3D& new_origin, const Vector_3D& y_axis, 
                const Point_3D& pt_xy_plane, const Point_3D& ref_origin=Point_3D(0,0,0));
        Mesh_3D& move_y_pyz(const Point_3D& new_origin, const Vector_3D& y_axis, 
                const Point_3D& pt_yz_plane, const Point_3D& ref_origin=Point_3D(0,0,0));
        Mesh_3D& move_z_pxz(const Point_3D& new_origin, const Vector_3D& z_axis, 
                const Point_3D& pt_xz_plane, const Point_3D& ref_origin=Point_3D(0,0,0));
        Mesh_3D& move_z_pyz(const Point_3D& new_origin, const Vector_3D& z_axis, 
                const Point_3D& pt_yz_plane, const Point_3D& ref_origin=Point_3D(0,0,0));
        Mesh_3D& operator+=(const Vector_3D&);
        // exception safety: no throw
        Mesh_3D& operator-=(const Vector_3D&);
        // exception safety: basic guarantee
        Mesh_3D& operator*=(const Measurement);
    private:
        Measurement precision;
        vector<shared_ptr<Point_3D>> point_list;
        vector<Facet> facet_list;
    };

    /*
     * Determine if a point is on or inside the mesh.  If the point is on or inside the
     * mesh, then pt_on_surface is set to true if the point is on the surface, or false
     * it the point is inside the mesh
     */
    const bool mesh_contains_point(const Mesh_3D& mesh, const Point_3D& p, bool& pt_on_surface);
}

#endif /* MESH2_3D_H */

