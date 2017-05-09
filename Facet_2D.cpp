/* 
 * File:   Facet_2D.cpp
 * Author: Jeffrey Davis
 * 
 */

#include "Facet_2D.h"

namespace VCAD_lib
{
    Facet_2D::Facet_2D() : p1(shared_ptr<Point_2D>()), p2(shared_ptr<Point_2D>()), p3(shared_ptr<Point_2D>()) {}
    
    Facet_2D::Facet_2D(const shared_ptr<Point_2D> pt1, const shared_ptr<Point_2D> pt2, 
            const shared_ptr<Point_2D> pt3) : p1(pt1), p2(pt2), p3(pt3)
    {
        validate_points();
        
        // swap points to make sure of counter clockwise point order
        if (cross_product(Vector_2D(*pt1, *pt2), Vector_2D(*pt1, *pt3)) < 0)
            swap(p2, p3);
    }

    void Facet_2D::validate_points()
    {
        // check if p1 is the same as p2
        if (p1->get_x() == p2->get_x() && p1->get_y() == p2->get_y())
            throw runtime_error("invalid facet points: point 1 is the same as point 2");
        
        // check if p1 is the same as p3
        if (p1->get_x() == p3->get_x() && p1->get_y() == p3->get_y())
            throw runtime_error("invalid facet points: point 1 is the same as point 3");
        
        // check if p2 is the same as p3
        if (p2->get_x() == p3->get_x() && p2->get_y() == p3->get_y())
            throw runtime_error("invalid facet points: point 2 is the same as point 3");
        
        // check if the three points are in a straight line
        // test if the cross product is exactly zero
        if (cross_product(Vector_2D(*p1, *p2), Vector_2D(*p1, *p3)) == 0)
            throw runtime_error("invalid facet points: points do not form a triangle, but a straight line");
    }
    
    const Point_2D Facet_2D::get_inside_point() const
    {
        Vector_2D v1(*p1,*p2);
        Vector_2D v2(*p2,*p3);
        return *p1 + v1 * 0.5 + v2 * 0.25;
    }

    const bool Facet_2D::contains_point(const Point_2D& pt, bool& pt_is_on_side, Measurement precision) const
    {
        // if the point is on p1-p2 or p1-p3, then it is contained in the facet
        if (is_pt_on_vector(pt, *p1, *p2, precision) || 
                is_pt_on_vector(pt, *p1, *p3, precision) || 
                is_pt_on_vector(pt, *p2, *p3, precision))
        {
//            cout << "contains_point: point is on a vector returning true\n";
            pt_is_on_side = true;
            return true;
        }

        Vector_2D::Measurement cp_v1v2(cross_product(Vector_2D(*p1, *p2), Vector_2D(*p1, *p3)));
        Vector_2D::Measurement cp_v1vp(cross_product(Vector_2D(*p1, *p2), Vector_2D(*p1, pt)));
        Vector_2D::Measurement cp_vpv2(cross_product(Vector_2D(*p1, pt), Vector_2D(*p1, *p3)));

//        cout << "contains_point: cp_v1v2 x: " << cp_v1v2.get_x() << " y: " << cp_v1v2.get_y() << "\n";
//        cout << "contains_point: cp_v1vp x: " << cp_v1vp.get_x() << " y: " << cp_v1vp.get_y() << "\n";
//        cout << "contains_point: cp_vpv2 x: " << cp_vpv2.get_x() << " y: " << cp_vpv2.get_y() << "\n";
//        cout << "contains_point: p1 dot_product(cp_v1v2, cp_v1vp): " << dot_product(cp_v1v2, cp_v1vp) << "\n";
//        cout << "contains_point: p1 dot_product(cp_v1v2, cp_vpv2): " << dot_product(cp_v1v2, cp_vpv2) << "\n";
        if (((cp_v1v2 > 0 && cp_v1vp > 0) || (cp_v1v2 < 0 && cp_v1vp < 0)) && ((cp_v1v2 > 0 && cp_vpv2 > 0) || (cp_v1v2 < 0 && cp_vpv2 < 0)))
        {
        
            cp_v1v2 = cross_product(Vector_2D(*p2, *p3), Vector_2D(*p2, *p1));
            cp_v1vp = cross_product(Vector_2D(*p2, *p3), Vector_2D(*p2, pt));
            cp_vpv2 = cross_product(Vector_2D(*p2, pt), Vector_2D(*p2, *p1));

//            cout << "contains_point: cp_v1v2 x: " << cp_v1v2.get_x() << " y: " << cp_v1v2.get_y() << "\n";
//            cout << "contains_point: cp_v1vp x: " << cp_v1vp.get_x() << " y: " << cp_v1vp.get_y() << "\n";
//            cout << "contains_point: cp_vpv2 x: " << cp_vpv2.get_x() << " y: " << cp_vpv2.get_y() << "\n";
//            cout << "contains_point: p1 dot_product(cp_v1v2, cp_v1vp): " << dot_product(cp_v1v2, cp_v1vp) << "\n";
//            cout << "contains_point: p1 dot_product(cp_v1v2, cp_vpv2): " << dot_product(cp_v1v2, cp_vpv2) << "\n";
            if (((cp_v1v2 > 0 && cp_v1vp > 0) || (cp_v1v2 < 0 && cp_v1vp < 0)) && ((cp_v1v2 > 0 && cp_vpv2 > 0) || (cp_v1v2 < 0 && cp_vpv2 < 0)))
            {

                cp_v1v2 = cross_product(Vector_2D(*p3, *p1), Vector_2D(*p3, *p2));
                cp_v1vp = cross_product(Vector_2D(*p3, *p1), Vector_2D(*p3, pt));
                cp_vpv2 = cross_product(Vector_2D(*p3, pt), Vector_2D(*p3, *p2));

//                cout << "contains_point: cp_v1v2 x: " << cp_v1v2.get_x() << " y: " << cp_v1v2.get_y() << "\n";
//                cout << "contains_point: cp_v1vp x: " << cp_v1vp.get_x() << " y: " << cp_v1vp.get_y() << "\n";
//                cout << "contains_point: cp_vpv2 x: " << cp_vpv2.get_x() << " y: " << cp_vpv2.get_y() << "\n";
//                cout << "contains_point: p1 dot_product(cp_v1v2, cp_v1vp): " << dot_product(cp_v1v2, cp_v1vp) << "\n";
//                cout << "contains_point: p1 dot_product(cp_v1v2, cp_vpv2): " << dot_product(cp_v1v2, cp_vpv2) << "\n";
                if (((cp_v1v2 > 0 && cp_v1vp > 0) || (cp_v1v2 < 0 && cp_v1vp < 0)) && ((cp_v1v2 > 0 && cp_vpv2 > 0) || (cp_v1v2 < 0 && cp_vpv2 < 0)))
                {
                    pt_is_on_side = false;
                    return true;
                }
            }
        }
        
        return false;
    }
    
    const bool is_equal(const Facet_2D& f1, const Facet_2D& f2, const Facet_2D::Measurement precision)
    {
        return (is_equal(*f1.get_point1(), *f2.get_point1(), precision) && is_equal(*f1.get_point2(), *f2.get_point2(), precision) && is_equal(*f1.get_point3(), *f2.get_point3(), precision)) || 
                (is_equal(*f1.get_point1(), *f2.get_point2(), precision) && is_equal(*f1.get_point2(), *f2.get_point3(), precision) && is_equal(*f1.get_point3(), *f2.get_point1(), precision)) || 
                (is_equal(*f1.get_point1(), *f2.get_point3(), precision) && is_equal(*f1.get_point2(), *f2.get_point1(), precision) && is_equal(*f1.get_point3(), *f2.get_point2(), precision));
    }
}
