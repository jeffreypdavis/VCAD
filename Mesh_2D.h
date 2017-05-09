/* 
 * File:   Mesh_2D.h
 * Author: Jeffrey Davis
 *
 */

#ifndef MESH_2D_H
#define MESH_2D_H

#include <vector>
#include <memory>
#include <iterator>
#include "Point_2D.h"
#include "Facet.h"
#include "Facet_2D.h"

namespace VCAD_lib
{

    class Mesh_2D {
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

            const_iterator(const vector<shared_ptr<Point_2D>>::const_iterator point_it_begin, 
                    const vector<Facet>::const_iterator facet_it_begin, 
                    const vector<Facet>::const_iterator facet_it_end, 
                    const vector<Facet>::const_iterator position);

            bool operator==(const const_iterator&) const;
            bool operator!=(const const_iterator&) const;

            const_iterator& operator++();
            const_iterator operator++(int);
            const_iterator& operator--();
            const_iterator operator--(int);

            const Facet_2D& operator*() const;
            const Facet_2D* const operator->() const;
        private:
            vector<shared_ptr<Point_2D>>::const_iterator point_list_begin;
            vector<Facet>::const_iterator facets_begin;
            vector<Facet>::const_iterator facets_end;
            vector<Facet>::const_iterator current_facet;
            Facet_2D facet;
            void update_facet();
        };        
        
        typedef vector<Facet>::size_type size_type;
        typedef Point_2D::Measurement Measurement;
        typedef Point_2D::Angle_Meas Angle_Meas;
        // to allow for users to get point list and facet point indices
        typedef vector<shared_ptr<Point_2D>>::const_iterator const_point_iterator;
        typedef vector<Facet>::const_iterator const_facet_iterator;
        
        // exception safety: strong guarantee
        Mesh_2D();
        // exception safety: strong guarantee - invalid_argument if precision is less than or equal to zero
        explicit Mesh_2D(const Measurement precision); // precision
        // exception safety: strong guarantee
        Mesh_2D(const Mesh_2D& orig);
        // exception safety: strong guarantee
        Mesh_2D& operator=(const Mesh_2D&);
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
        
        void push_back(const Facet_2D& facet);
        size_type size() const { return facet_list.size(); }
        bool empty() const { return facet_list.empty(); }
        void clear();
        const_iterator erase(const_iterator it);
        const_iterator erase(const_iterator begin, const_iterator end);
//        const_iterator insert(const_iterator loc, const Facet_2D& facet);
//        const_iterator insert(const_iterator loc, const_iterator from_facet, const_iterator to_facet);
        // rotate
        // exception safety: basic guarantee
        Mesh_2D& rotate(const Angle_Meas angle);
        // exception safety: basic guarantee
        Mesh_2D& rotate(const Angle_Meas angle, const Point_2D& origin);
        // scale
        // exception safety: basic guarantee
        Mesh_2D& scale(const Measurement x_scalar, const Measurement y_scalar);
        // exception safety: basic guarantee
        Mesh_2D& scale(const Measurement x_scalar, const Measurement y_scalar, 
                const Point_2D& origin);
        // translate
        // exception safety: no throw
        Mesh_2D& translate(const Measurement x_val, const Measurement y_val);
        // exception safety: no throw
        Mesh_2D& translate(const Vector_2D&);
        Mesh_2D& move(const Point_2D& new_origin, const Vector_2D& axis,
                const bool is_x_axis, const Point_2D& ref_origin=Point_2D(0,0));
        Mesh_2D& operator+=(const Vector_2D&);
        // exception safety: no throw
        Mesh_2D& operator-=(const Vector_2D&);
        // exception safety: basic guarantee
        Mesh_2D& operator*=(const Measurement);
    private:
        Measurement precision;
        vector<shared_ptr<Point_2D>> point_list;
        vector<Facet> facet_list;
        const bool validate_facet(const shared_ptr<Point_2D> p1, const shared_ptr<Point_2D> p2, const shared_ptr<Point_2D> p3);
        const Vector_2D gen_facet_unv(const Facet& facet);
    };

    class Mesh_3D;
    
    Mesh_3D& linear_extrude(const Mesh_2D&, Mesh_3D&, const Point_2D::Measurement height, const bool center=false);
    
}

#endif /* MESH2_2D_H */

