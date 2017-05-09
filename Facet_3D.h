/* 
 * File:   Facet_3D.h
 * Author: Jeffrey Davis
 *
 */

#ifndef FACET_3D_H
#define FACET_3D_H

#include <memory>
#include "Point_3D.h"
#include "Vector_3D.h"

namespace VCAD_lib
{

    /*
     * Facets are maintained to have counter-clockwise point order.  This way
     * the unit normal vector will follow the right hand rule:
     * 
     * Counter-clockwise point order:
     * P3 -- P2
     *  \    /
     *    P1
     * Right hand rule:
     * Align your hand with fingers out along P1-P2.
     * Close your fingers moving from P1-P2 to P1-P3.  Thumb will point in
     * direction of unit normal vector
     */
    class Facet_3D {
    public:
        typedef Point_3D::Measurement Measurement;
        Facet_3D();
        Facet_3D(const shared_ptr<Point_3D> pt1, const shared_ptr<Point_3D> pt2, 
                const shared_ptr<Point_3D> pt3, const bool clockwise_order=false);
        Facet_3D(const shared_ptr<Point_3D> pt1, const shared_ptr<Point_3D> pt2, 
                const shared_ptr<Point_3D> pt3, const Vector_3D& un_vector, 
                const bool clockwise_order=false);
        // exception safety: no throw
        const shared_ptr<Point_3D> get_point1() const { return p1; }
        // exception safety: no throw
        const shared_ptr<Point_3D> get_point2() const { return p2; }
        // exception safety: no throw
        const shared_ptr<Point_3D> get_point3() const { return p3; }
        // exception safety: strong guarantee throws length_error if calculated vector has zero length
        const Vector_3D get_unv() const;
        // exception safety: no throw
        void invert_unv();
        // exception safety: no throw
        const Point_3D get_inside_point() const;
        const bool contains_point(const Point_3D& pt, bool& pt_is_on_side, const Measurement precision) const;
    private:
        shared_ptr<Point_3D> p1;
        shared_ptr<Point_3D> p2;
        shared_ptr<Point_3D> p3;
        void validate_points();
    };

    /* 
     * determines if two facets are within a rounding error of each other
     */
    const bool is_equal(const Facet_3D& f1, const Facet_3D& f2, const Facet_3D::Measurement precision);
//    bool within_round(const Facet_3D& v1, const Facet_3D& v2, const Point_3D::Measurement precision);

    // determine if the specified point is on the facet plane
    // exception safety: no throw
    const bool is_pt_on_facet_plane(const Point_3D& point, const Facet_3D& facet, 
            const Facet_3D::Measurement precision);
    
    // line is defined by the vector and point.  Facet plane is defined by f.
    // returns true if an intersect point was found.  Intersect point is stored
    // in the i_point argument
    // exception safety: no throw
    const bool intersect_line_facet_plane(const Vector_3D& v, const Point_3D& o, 
            const Facet_3D& f, Point_3D& i_point, const Facet_3D::Measurement precision);
}

#endif /* FACET2_3D_H */

