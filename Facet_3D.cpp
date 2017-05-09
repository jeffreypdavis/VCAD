/* 
 * File:   Facet_3D.cpp
 * Author: Jeffrey Davis
 * 
 */

#include "Facet_3D.h"

namespace VCAD_lib
{
    Facet_3D::Facet_3D() : p1(shared_ptr<Point_3D>()), p2(shared_ptr<Point_3D>()), p3(shared_ptr<Point_3D>()) {}
    
    Facet_3D::Facet_3D(const shared_ptr<Point_3D> pt1, const shared_ptr<Point_3D> pt2, 
            const shared_ptr<Point_3D> pt3, const bool clockwise_order) : p1(pt1), p2(pt2), p3(pt3)
    {
        validate_points();
        
        // change point order to counter-clockwise orientation
        if (clockwise_order)
            swap(p2, p3);
    }

    Facet_3D::Facet_3D(const shared_ptr<Point_3D> pt1, const shared_ptr<Point_3D> pt2, 
            const shared_ptr<Point_3D> pt3, const Vector_3D& un_vector, 
            const bool clockwise_order) : p1(pt1), p2(pt2), p3(pt3)
    { 
        validate_points();
        
        // change point order to counter-clockwise if clockwise
        if (clockwise_order)
            swap(p2, p3);

        if (un_vector.length() > 0)
        {
            if (dot_product(this->get_unv(), un_vector) < 0) // point order/unv error
                throw invalid_argument("Error: point order and unit normal vector do not align.");
        }
    }
    
    void Facet_3D::validate_points()
    {
        // check if p1 is the same as p2
        if (p1->get_x() == p2->get_x() && p1->get_y() == p2->get_y() && p1->get_z() == p2->get_z())
            throw runtime_error("invalid facet points: point 1 is the same as point 2");
        
        // check if p1 is the same as p3
        if (p1->get_x() == p3->get_x() && p1->get_y() == p3->get_y() && p1->get_z() == p3->get_z())
            throw runtime_error("invalid facet points: point 1 is the same as point 3");
        
        // check if p2 is the same as p3
        if (p2->get_x() == p3->get_x() && p2->get_y() == p3->get_y() && p2->get_z() == p3->get_z())
            throw runtime_error("invalid facet points: point 2 is the same as point 3");
        
        // check if the three points are in a straight line
        // test if the cross product is exactly zero
        Vector_3D cp = cross_product(Vector_3D(*p1, *p2), Vector_3D(*p1, *p3));
        if (cp.get_x() == 0 && cp.get_y() == 0 && cp.get_z() == 0)
            throw runtime_error("invalid facet points: points do not form a triangle, but a straight line");
    }

    const Vector_3D Facet_3D::get_unv() const
    {
        Vector_3D unv(cross_product(Vector_3D(*p1, *p2), Vector_3D(*p1, *p3)));
        unv.normalize();
        return unv;
    }
    
    void Facet_3D::invert_unv()
    {
        // reorder points
        swap(p2, p3);
    }
    
    const Point_3D Facet_3D::get_inside_point() const
    {
        Vector_3D v1(*p1,*p2);
        Vector_3D v2(*p2,*p3);
        return *p1 + v1 * 0.5 + v2 * 0.25;
    }

    const bool Facet_3D::contains_point(const Point_3D& pt, bool& pt_is_on_side, Measurement precision) const
    {
        // make sure point is on facet plane
        if (!is_pt_on_facet_plane(pt, *this, precision))
        {
//            cout << "contains_point: point is not on facet plane returning false\n";
            return false;
        }
        
        // if the point is on p1-p2 or p1-p3, then it is contained in the facet
        if (is_pt_on_vector(pt, *p1, *p2, precision) || 
                is_pt_on_vector(pt, *p1, *p3, precision) || 
                is_pt_on_vector(pt, *p2, *p3, precision))
        {
//            cout << "contains_point: point is on a vector returning true\n";
            pt_is_on_side = true;
            return true;
        }

        Vector_3D cp_v1v2(cross_product(Vector_3D(*p1, *p2), Vector_3D(*p1, *p3)));
        Vector_3D cp_v1vp(cross_product(Vector_3D(*p1, *p2), Vector_3D(*p1, pt)));
        Vector_3D cp_vpv2(cross_product(Vector_3D(*p1, pt), Vector_3D(*p1, *p3)));

//        cout << "contains_point: cp_v1v2 x: " << cp_v1v2.get_x() << " y: " << cp_v1v2.get_y() << " z: " << cp_v1v2.get_z() << "\n";
//        cout << "contains_point: cp_v1vp x: " << cp_v1vp.get_x() << " y: " << cp_v1vp.get_y() << " z: " << cp_v1vp.get_z() << "\n";
//        cout << "contains_point: cp_vpv2 x: " << cp_vpv2.get_x() << " y: " << cp_vpv2.get_y() << " z: " << cp_vpv2.get_z() << "\n";
//        cout << "contains_point: p1 dot_product(cp_v1v2, cp_v1vp): " << dot_product(cp_v1v2, cp_v1vp) << "\n";
//        cout << "contains_point: p1 dot_product(cp_v1v2, cp_vpv2): " << dot_product(cp_v1v2, cp_vpv2) << "\n";
//        bool from_p1 = dot_product(cp_v1v2, cp_v1vp) > 0 && dot_product(cp_v1v2, cp_vpv2) > 0;
        if (dot_product(cp_v1v2, cp_v1vp) > 0 && dot_product(cp_v1v2, cp_vpv2) > 0)
        {
        
            cp_v1v2 = cross_product(Vector_3D(*p2, *p3), Vector_3D(*p2, *p1));
            cp_v1vp = cross_product(Vector_3D(*p2, *p3), Vector_3D(*p2, pt));
            cp_vpv2 = cross_product(Vector_3D(*p2, pt), Vector_3D(*p2, *p1));

    //        cout << "contains_point: cp_v1v2 x: " << cp_v1v2.get_x() << " y: " << cp_v1v2.get_y() << " z: " << cp_v1v2.get_z() << "\n";
    //        cout << "contains_point: cp_v1vp x: " << cp_v1vp.get_x() << " y: " << cp_v1vp.get_y() << " z: " << cp_v1vp.get_z() << "\n";
    //        cout << "contains_point: cp_vpv2 x: " << cp_vpv2.get_x() << " y: " << cp_vpv2.get_y() << " z: " << cp_vpv2.get_z() << "\n";
    //        cout << "contains_point: p1 dot_product(cp_v1v2, cp_v1vp): " << dot_product(cp_v1v2, cp_v1vp) << "\n";
    //        cout << "contains_point: p1 dot_product(cp_v1v2, cp_vpv2): " << dot_product(cp_v1v2, cp_vpv2) << "\n";
//            bool from_p2 = dot_product(cp_v1v2, cp_v1vp) > 0 && dot_product(cp_v1v2, cp_vpv2) > 0;
            if (dot_product(cp_v1v2, cp_v1vp) > 0 && dot_product(cp_v1v2, cp_vpv2) > 0)
            {

                cp_v1v2 = cross_product(Vector_3D(*p3, *p1), Vector_3D(*p3, *p2));
                cp_v1vp = cross_product(Vector_3D(*p3, *p1), Vector_3D(*p3, pt));
                cp_vpv2 = cross_product(Vector_3D(*p3, pt), Vector_3D(*p3, *p2));

        //        cout << "contains_point: cp_v1v2 x: " << cp_v1v2.get_x() << " y: " << cp_v1v2.get_y() << " z: " << cp_v1v2.get_z() << "\n";
        //        cout << "contains_point: cp_v1vp x: " << cp_v1vp.get_x() << " y: " << cp_v1vp.get_y() << " z: " << cp_v1vp.get_z() << "\n";
        //        cout << "contains_point: cp_vpv2 x: " << cp_vpv2.get_x() << " y: " << cp_vpv2.get_y() << " z: " << cp_vpv2.get_z() << "\n";
        //        cout << "contains_point: p1 dot_product(cp_v1v2, cp_v1vp): " << dot_product(cp_v1v2, cp_v1vp) << "\n";
        //        cout << "contains_point: p1 dot_product(cp_v1v2, cp_vpv2): " << dot_product(cp_v1v2, cp_vpv2) << "\n";
//                bool from_p3 = dot_product(cp_v1v2, cp_v1vp) > 0 && dot_product(cp_v1v2, cp_vpv2) > 0;
                if (dot_product(cp_v1v2, cp_v1vp) > 0 && dot_product(cp_v1v2, cp_vpv2) > 0)
                {

            //        cout << "contains_point: from_p1: " << from_p1 << " from_p2: " << from_p2 << " from_p3: " << from_p3 << "\n";
            //        return (from_p1 && from_p2) || (from_p1 && from_p3) || (from_p2 && from_p3);
//                    return from_p1 && from_p2 && from_p3;
                    pt_is_on_side = false;
                    return true;
                }
            }
        }
        
        return false;
    }
    
    const bool is_equal(const Facet_3D& f1, const Facet_3D& f2, const Facet_3D::Measurement precision)
    {
        return (is_equal(*f1.get_point1(), *f2.get_point1(), precision) && is_equal(*f1.get_point2(), *f2.get_point2(), precision) && is_equal(*f1.get_point3(), *f2.get_point3(), precision)) || 
                (is_equal(*f1.get_point1(), *f2.get_point2(), precision) && is_equal(*f1.get_point2(), *f2.get_point3(), precision) && is_equal(*f1.get_point3(), *f2.get_point1(), precision)) || 
                (is_equal(*f1.get_point1(), *f2.get_point3(), precision) && is_equal(*f1.get_point2(), *f2.get_point1(), precision) && is_equal(*f1.get_point3(), *f2.get_point2(), precision));
    }
    
    const bool is_pt_on_facet_plane(const Point_3D& point, const Facet_3D& facet, 
            const Facet_3D::Measurement precision)
    {
        if (is_equal(point, *facet.get_point1(), precision))
            return true;
        
        // The equation of the plane P through (x0, y0, z0) that has a normal vector 
        // n=Ai + Bj + ck is:
        //    A(x - x0) + B(y - y0) + C(z - z0) = 0
        Vector_3D unv(facet.get_unv());
        Facet_3D::Measurement result(unv.get_x() * (point.get_x() - facet.get_point1()->get_x()));
        Facet_3D::Measurement largest(fabs(result));
        Facet_3D::Measurement next(unv.get_y() * (point.get_y() - facet.get_point1()->get_y()));
        result += next;
        Facet_3D::Measurement error_bound(fabs(next));
        if (error_bound > largest)
            largest = error_bound;
        next = unv.get_z() * (point.get_z() - facet.get_point1()->get_z());
        result += next;
        error_bound = fabs(next);
        if (error_bound > largest)
            largest = error_bound;
        
        error_bound = largest > 1.0 ? (largest * precision) : precision;
//        cout << "is_pt_on_facet_plane: result: " << result << " largest: " << largest << " (precision * largest): " << (precision * largest) << " needs a precision of DBL_EPSILON * " << (fabs(result) / (largest * DBL_EPSILON)) << "\n";
        return fabs(result) <= error_bound;
    }
    
    const bool intersect_line_facet_plane(const Vector_3D& v, const Point_3D& o, 
            const Facet_3D& facet, Point_3D& i_point, const Facet_3D::Measurement precision)
    {
        // facet - the facet to check if the line intersects
        // v - the vector forming the line
        // o - the vector's origin
        //
        // The parametric equations of a line l through points P=(x1,y1,z1) and Q=(x2,y2,z2)
        // are:
        //    x=x1 + (x2 - x1)t
        //    y=y1 + (y2 - y1)t
        //    z=z1 + (z2 - z1)t
        // where (x,y,z) is the general point of l, and the parameter t takes on all real values.
        //
        // line equations
        // x=p1X + (p2X - p1X)t
        // y=p1Y + (p2Y - p1Y)t
        // z=p1Z + (p2Z - p2Z)t
        //
        // plane equation
        // The equation of the plane P through (x0, y0, z0) that has a normal vector 
        // n=Ai + Bj + ck is:
        //    A(x - x0) + B(y - y0) + C(z - z0) = 0
        //
        // substitute line equations into equation
        // A(p1X + (p2X - p1X)t - x0) + B(p1Y + (p2Y - p1Y)t - y0) + C(p1Z + (p2Z - p1Z)t - z0) = 0
        // Ap1X + A(p2X - p1X)t - Ax0 + Bp1Y + B(p2Y - p1Y)t - By0 + Cp1Z + C(p2Z - p1Z)t - Cz0 = 0
        // A(p2X - p1X)t + B(p2Y - p1Y)t + C(p2Z - p1Z)t = -Ap1X + Ax0 - Bp1Y + By0 - Cp1Z + Cz0
        // t(A(p2X - p1X) + B(p2Y - p1Y) + C(p2Z - p1Z)) = Ax0  + By0 + Cz0 - Ap1X - Bp1Y - Cp1Z
        // t = (Ax0  + By0 + Cz0 - Ap1X - Bp1Y - Cp1Z) / (A(p2X - p1X) + B(p2Y - p1Y) + C(p2Z - p1Z))
        // t = (A(x0 - p1X) + B(y0 - p1Y) + C(z0 - p1Z)) / (A(p2X - p1X) + B(p2Y - p1Y) + C(p2Z - p1Z))
        //
        // solve x,y,z using line equations
        
        // A(x - x0) + B(y - y0) + C(z - z0) = 0
        // A(p1X + t*vX - x0) + B(p1Y + t*vY - y0) + C(p1Z + t*vZ - z0) = 0
        // A*p1X + t*A*vX - A*x0 + B*p1Y + t*B*vY - B*y0 + C*p1Z + t*C*vZ - C*z0 = 0
        // t*(A*vX + B*vY + C*vZ) = A*x0 - A*p1X + B*y0 - B*p1Y + C*z0 - C*p1Z
        // t = (A*x0 - A*p1X + B*y0 - B*p1Y + C*z0 - C*p1Z) / (A*vX + B*vY + C*vZ)
        //
        // possible points to use in plane equation
        // unit normal vector to use in plane equation
        Vector_3D unv = facet.get_unv();
        // the following two points are used for the line equations
        Facet_3D::Measurement bottom(unv.get_x() * v.get_x());
        Facet_3D::Measurement largest(fabs(bottom));
        Facet_3D::Measurement next(unv.get_y() * v.get_y());
        Facet_3D::Measurement error_bound(fabs(next));
        bottom += next;
        if (error_bound > largest)
            largest = error_bound;
        next = unv.get_z() * v.get_z();
        error_bound = fabs(next);
        bottom += next;
        if (error_bound > largest)
            largest = error_bound;

        error_bound = largest > 1.0 ? (largest * precision) : precision;
        if (fabs(bottom) > error_bound)
        {
//            cout << "ilfp: calculating from facet p1\n";
            // t = (A*x0 - A*p1X + B*y0 - B*p1Y + C*z0 - C*p1Z) / (A*vX + B*vY + C*vZ)
            Facet_3D::Measurement t(unv.get_x() * (facet.get_point1()->get_x() - o.get_x()) + unv.get_y() * (facet.get_point1()->get_y() - o.get_y()) + unv.get_z() * (facet.get_point1()->get_z() - o.get_z()));
            t /= bottom;
            // check if point does lie on the plane
            Point_3D p = o + v * t;

//            cout << "ilfp i_point x: " << p.get_x() << " y: " << p.get_y() << " z: " << p.get_z() << "\n";
            if (is_pt_on_facet_plane(p, facet, precision))
            {
                i_point = p;
                return true;
            }
            
//            cout << "ilfp: calculating from facet p2\n";
            // t = (A*x0 - A*p1X + B*y0 - B*p1Y + C*z0 - C*p1Z) / (A*vX + B*vY + C*vZ)
            t = unv.get_x() * (facet.get_point2()->get_x() - o.get_x()) + unv.get_y() * (facet.get_point2()->get_y() - o.get_y()) + unv.get_z() * (facet.get_point2()->get_z() - o.get_z());
            t /= bottom;
            // check if point does lie on the plane
            p = o + v * t;

//            cout << "ilfp i_point x: " << p.get_x() << " y: " << p.get_y() << " z: " << p.get_z() << "\n";
            if (is_pt_on_facet_plane(p, facet, precision))
            {
                i_point = p;
                return true;
            }
            
//            cout << "ilfp: calculating from facet p3\n";
            // t = (A*x0 - A*p1X + B*y0 - B*p1Y + C*z0 - C*p1Z) / (A*vX + B*vY + C*vZ)
            t = unv.get_x() * (facet.get_point3()->get_x() - o.get_x()) + unv.get_y() * (facet.get_point3()->get_y() - o.get_y()) + unv.get_z() * (facet.get_point3()->get_z() - o.get_z());
            t /= bottom;
            // check if point does lie on the plane
            p = o + v * t;

//            cout << "ilfp i_point x: " << p.get_x() << " y: " << p.get_y() << " z: " << p.get_z() << "\n";
            if (is_pt_on_facet_plane(p, facet, precision))
            {
                i_point = p;
                return true;
            }
        }
        else if (is_pt_on_facet_plane(o, facet, precision)) // vector is parallel, so see if origin is on facet plane
        {
            i_point = o;
            return true;
        }
        
        return false;
    }
}
