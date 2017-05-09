/* 
 * File:   Facet_2D.h
 * Author: Jeffrey Davis
 *
 */

#ifndef FACET_2D_H
#define FACET_2D_H

#include <memory>
#include "Point_2D.h"
#include "Vector_2D.h"

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
    class Facet_2D {
    public:
        typedef Point_2D::Measurement Measurement;
        Facet_2D();
        Facet_2D(const shared_ptr<Point_2D> pt1, const shared_ptr<Point_2D> pt2, 
                const shared_ptr<Point_2D> pt3);
        // exception safety: no throw
        const shared_ptr<Point_2D> get_point1() const { return p1; }
        // exception safety: no throw
        const shared_ptr<Point_2D> get_point2() const { return p2; }
        // exception safety: no throw
        const shared_ptr<Point_2D> get_point3() const { return p3; }
        // exception safety: no throw
        const Point_2D get_inside_point() const;
        const bool contains_point(const Point_2D& pt, bool& pt_is_on_side, const Measurement precision) const;
    private:
        shared_ptr<Point_2D> p1;
        shared_ptr<Point_2D> p2;
        shared_ptr<Point_2D> p3;
        void validate_points();
    };

    /* 
     * determines if two facets are within a rounding error of each other
     */
    const bool is_equal(const Facet_2D& f1, const Facet_2D& f2, const Facet_2D::Measurement precision);
}
#endif /* FACET2_2D_H */

