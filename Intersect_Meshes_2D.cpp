/* 
 * File:   Intersect_Meshes_2D.cpp
 * Author: Jeffrey Davis
 * 
 * Available debug preprocessor definitions
 * DEBUG_INTERSECT_MESHES_2D
 * DEBUG_INTERSECT_MESHES_2D_I_PT_LIST
 * DEBUG_INTERSECT_MESHES_2D_I_PT_LOCATOR
 * DEBUG_INTERSECT_MESHES_2D_FACETS
 * DEBUG_INTERSECT_MESHES_2D_FACET_BUILDER
 * DEBUG_INTERSECT_MESHES_2D_FACET_BUILDER_SEGMENTS
 * DEBUG_INTERSECT_MESHES_2D_TRACE_FACETS
 * DEBUG_INTERSECT_MESHES_2D_FACET_SORTER
 * DEBUG_INTERSECT_MESHES_2D_DIFFERENCE
 * DEBUG_INTERSECT_MESHES_2D_INTERSECTION
 * DEBUG_INTERSECT_MESHES_2D_MERGE
 */

#include "Intersect_Meshes_2D.h"
#include <algorithm>
#include <vector>
#include <stack>
#include <cmath>
#include <cfloat>

namespace VCAD_lib
{

    Intersect_Meshes_2D::Intersect_Point::Intersect_Point() : 
            pt(), f1_loc(Location::internal), f2_loc(Location::internal) {}
    
    Intersect_Meshes_2D::Intersect_Point::Intersect_Point(const shared_ptr<Point_2D> i_point, 
            const Location f1_location, const Location f2_location) : pt(i_point), 
            f1_loc(f1_location), f2_loc(f2_location) {}
    
    const bool Intersect_Meshes_2D::Intersect_Point::operator==(const Intersect_Point& ipt) const
    {
        return (ipt.pt == pt || (ipt.pt->get_x() == pt->get_x() && ipt.pt->get_y() == pt->get_y())) && 
                ipt.f1_loc == f1_loc && ipt.f2_loc == f2_loc;
    }
    
    Intersect_Meshes_2D::I_Pt_List::I_Pt_List() : i_points() {}
    
    const bool Intersect_Meshes_2D::I_Pt_List::matches(const Point_2D& p1, const Point_2D& p2) const
    {
        return p1.get_x() == p2.get_x() && p1.get_y() == p2.get_y();
    }
    
    Intersect_Meshes_2D::I_Pt_List::const_iterator Intersect_Meshes_2D::I_Pt_List::process_matching_i_pts(
            const Intersect_Point& matching_ip1, const const_iterator matching_ip2, 
            const Facet_2D& facet1, const Facet_2D& facet2)
    {
#ifdef DEBUG_INTERSECT_MESHES_2D_I_PT_LIST
        cout << "Intersect_Meshes_2D::I_Pt_List::process_matching_i_pts begin\n";
        cout << "Intersect_Meshes_2D::I_Pt_List::process_matching_i_pts ip1 pt x: " << matching_ip1.pt->get_x() << " y: " << matching_ip1.pt->get_y() << " f1_loc: " << matching_ip1.f1_loc << " f2_loc: " << matching_ip1.f2_loc << "\n";
        cout << "Intersect_Meshes_2D::I_Pt_List::process_matching_i_pts ip2 pt x: " << (*matching_ip2).pt->get_x() << " y: " << (*matching_ip2).pt->get_y() << " f1_loc: " << (*matching_ip2).f1_loc << " f2_loc: " << (*matching_ip2).f2_loc << "\n";
#endif
        // only process if the side locations of either f1, or f2 are the same
        bool rem_ip1(false);
        bool rem_ip2(false);
        bool add_ip(false);
        shared_ptr<Point_2D> pt;
        Location f1_loc(Location::internal);
        Location f2_loc(Location::internal);
        
        bool f1_loc_same(false);
        Location ip1_loc(Location::internal);
        Location ip2_loc(Location::internal);
        Location common_loc(Location::internal);
        
        if (matching_ip1.f1_loc == (*matching_ip2).f1_loc)
        {
            f1_loc_same = true;
            ip1_loc = matching_ip1.f2_loc;
            ip2_loc = (*matching_ip2).f2_loc;
            common_loc = matching_ip1.f1_loc;
        }
        else if (matching_ip1.f2_loc == (*matching_ip2).f2_loc)
        {
            f1_loc_same = false;
            ip1_loc = matching_ip1.f1_loc;
            ip2_loc = (*matching_ip2).f1_loc;
            common_loc = matching_ip1.f2_loc;
        }
        else // no common locations - return
        {
            const_iterator it = find(i_points.begin(), i_points.end(), matching_ip1);
            if (it == end())
                throw runtime_error("Unable to locate existing intersect point");
            return ++it;
        }

        // process locations and determine which intersect point to use
        if (ip1_loc == Location::p1 || ip1_loc == Location::p2 || ip1_loc == Location::p3)
        {
            // prefer corner points over side or internal points
            if (ip2_loc == Location::p1p2 || ip2_loc == Location::p1p3 || ip2_loc == Location::p2p3 || ip2_loc == Location::internal)
                rem_ip2 = true;
        }
        else if (ip1_loc == Location::p1p2 || ip1_loc == Location::p1p3 || ip1_loc == Location::p2p3)
        {
            // prefer side points over internal points
            if (ip2_loc == Location::internal)
                rem_ip2 = true;
            else if (ip2_loc == Location::p1 || ip2_loc == Location::p2 || ip2_loc == Location::p3)
                rem_ip1 = true; // prefer corner points over side points
            else // same pt is on two sides... determine which one makes sense
            {
                // determine which location is the correct location...
#ifdef DEBUG_INTERSECT_MESHES_2D_I_PT_LIST
                cout << "Intersect_Meshes_2D::I_Pt_List::process_matching_i_pts matching points are both side points\n";
#endif
                // iterate through other intersect points and remove the side that generates too many side points
                vector<Intersect_Point> p1p2_pts;
                vector<Intersect_Point> p1p3_pts;
                vector<Intersect_Point> p2p3_pts;
                for (const_iterator it = i_points.begin(); it != i_points.end(); ++it)
                {
                    Location loc(f1_loc_same ? (*it).f2_loc : (*it).f1_loc);
                    if (loc == Location::p1)
                    {
                        p1p2_pts.push_back(*it);
                        p1p3_pts.push_back(*it);
                    }
                    else if (loc == Location::p2)
                    {
                        p1p2_pts.push_back(*it);
                        p2p3_pts.push_back(*it);
                    }
                    else if (loc == Location::p3)
                    {
                        p1p3_pts.push_back(*it);
                        p2p3_pts.push_back(*it);
                    }
                    else if (loc == Location::p1p2)
                        p1p2_pts.push_back(*it);
                    else if (loc == Location::p1p3)
                        p1p3_pts.push_back(*it);
                    else if (loc == Location::p2p3)
                        p2p3_pts.push_back(*it);
                }

                if (p1p2_pts.size() > 2)
                {
                    if (p1p2_pts.end() != find(p1p2_pts.begin(), p1p2_pts.end(), matching_ip1))
                        rem_ip1 = true;
                    else if (p1p2_pts.end() != find(p1p2_pts.begin(), p1p2_pts.end(), *matching_ip2))
                        rem_ip2 = true;
                }
                else if (p1p3_pts.size() > 2)
                {
                    if (p1p3_pts.end() != find(p1p3_pts.begin(), p1p3_pts.end(), matching_ip1))
                        rem_ip1 = true;
                    else if (p1p3_pts.end() != find(p1p3_pts.begin(), p1p3_pts.end(), *matching_ip2))
                        rem_ip2 = true;
                }
                else if (p2p3_pts.size() > 2)
                {
                    if (p2p3_pts.end() != find(p2p3_pts.begin(), p2p3_pts.end(), matching_ip1))
                        rem_ip1 = true;
                    else if (p2p3_pts.end() != find(p2p3_pts.begin(), p2p3_pts.end(), *matching_ip2))
                        rem_ip2 = true;
                }
                else
                {
                    shared_ptr<Point_2D> common_corner;
                    shared_ptr<Point_2D> f_diff_ep1;
                    shared_ptr<Point_2D> f_diff_ep2;
                    shared_ptr<Point_2D> f_same_ep1;
                    shared_ptr<Point_2D> f_same_ep2;

                    if (ip1_loc == Location::p1p2)
                    {
                        if (ip2_loc == Location::p1p3)
                        {
                            common_corner = f1_loc_same ? facet2.get_point1() : facet1.get_point1();
                            f_diff_ep1 = f1_loc_same ? facet2.get_point2() : facet1.get_point2();
                            f_diff_ep2 = f1_loc_same ? facet2.get_point3() : facet1.get_point3();
                        }
                        else // ip2_loc == p2p3
                        {
                            common_corner = f1_loc_same ? facet2.get_point2() : facet1.get_point2();
                            f_diff_ep1 = f1_loc_same ? facet2.get_point1() : facet1.get_point1();
                            f_diff_ep2 = f1_loc_same ? facet2.get_point3() : facet1.get_point3();
                        }
                    }
                    else if (ip1_loc == Location::p1p3)
                    {
                        if (ip2_loc == Location::p1p2)
                        {
                            common_corner = f1_loc_same ? facet2.get_point1() : facet1.get_point1();
                            f_diff_ep1 = f1_loc_same ? facet2.get_point3() : facet1.get_point3();
                            f_diff_ep2 = f1_loc_same ? facet2.get_point2() : facet1.get_point2();
                        }
                        else // ip2_loc == p2p3
                        {
                            common_corner = f1_loc_same ? facet2.get_point3() : facet1.get_point3();
                            f_diff_ep1 = f1_loc_same ? facet2.get_point1() : facet1.get_point1();
                            f_diff_ep2 = f1_loc_same ? facet2.get_point2() : facet1.get_point2();
                        }
                    }
                    else // ip1_loc == p2p3
                    {
                        if (ip2_loc == Location::p1p2)
                        {
                            common_corner = f1_loc_same ? facet2.get_point2() : facet1.get_point2();
                            f_diff_ep1 = f1_loc_same ? facet2.get_point3() : facet1.get_point3();
                            f_diff_ep2 = f1_loc_same ? facet2.get_point1() : facet1.get_point1();
                        }
                        else // ip2_loc == p1p3
                        {
                            common_corner = f1_loc_same ? facet2.get_point3() : facet1.get_point3();
                            f_diff_ep1 = f1_loc_same ? facet2.get_point2() : facet1.get_point2();
                            f_diff_ep2 = f1_loc_same ? facet2.get_point1() : facet1.get_point1();
                        }
                    }

                    // determine f1_ep1 and f1_ep2
                    switch (common_loc)
                    {
                        case (Location::p1):
                            f_same_ep1 = f1_loc_same ? facet1.get_point2() : facet2.get_point2();
                            f_same_ep2 = f1_loc_same ? facet1.get_point3() : facet2.get_point3();
                            break;
                        case (Location::p2):
                            f_same_ep1 = f1_loc_same ? facet1.get_point1() : facet2.get_point1();
                            f_same_ep2 = f1_loc_same ? facet1.get_point3() : facet2.get_point3();
                            break;
                        case (Location::p3):
                            f_same_ep1 = f1_loc_same ? facet1.get_point1() : facet2.get_point1();
                            f_same_ep2 = f1_loc_same ? facet1.get_point2() : facet2.get_point2();
                            break;
                        case (Location::p1p2):
                            f_same_ep1 = f1_loc_same ? facet1.get_point1() : facet2.get_point1();
                            f_same_ep2 = f1_loc_same ? facet1.get_point2() : facet2.get_point2();
                            break;
                        case (Location::p1p3):
                            f_same_ep1 = f1_loc_same ? facet1.get_point1() : facet2.get_point1();
                            f_same_ep2 = f1_loc_same ? facet1.get_point3() : facet2.get_point3();
                            break;
                        case (Location::p2p3):
                            f_same_ep1 = f1_loc_same ? facet1.get_point2() : facet2.get_point2();
                            f_same_ep2 = f1_loc_same ? facet1.get_point3() : facet2.get_point3();
                            break;
                        default: // Location::internal
                            throw runtime_error("unable to process matching point because location is internal");
                    }

                    // find the side that is closer to the f2 point and use that
                    // use cross product and dot product together to determine closest side
                    Vector_2D f_diff_v1(*common_corner, *f_diff_ep1);
                    Vector_2D f_diff_v2(*common_corner, *f_diff_ep2);
                    Vector_2D f_same_v1(*common_corner, *f_same_ep1);
                    Vector_2D f_same_v2(*common_corner, *f_same_ep2);

                    Vector_2D::Measurement cp_f_diff_v1v2(cross_product(f_diff_v1, f_diff_v2));
                    Vector_2D::Measurement cp_f_diff_v1_f_same_v1(cross_product(f_diff_v1, f_same_v1)); // could be a zero vector if end point is the same
                    Vector_2D::Measurement cp_f_diff_v1_f_same_v2(cross_product(f_diff_v1, f_same_v2)); // could be a zero vector if end point is the same
                    Vector_2D::Measurement cp_f_diff_v2_f_same_v1(cross_product(f_diff_v2, f_same_v1)); // could be a zero vector if end point is the same
                    Vector_2D::Measurement cp_f_diff_v2_f_same_v2(cross_product(f_diff_v2, f_same_v2)); // could be a zero vector if end point is the same

                    if (cp_f_diff_v1_f_same_v1 == 0) 
                    {
                        // cp_f_diff_v1_f_same_v1 is a zero vector
                        if ((cp_f_diff_v1_f_same_v2 > 0 && cp_f_diff_v2_f_same_v2 > 0) || (cp_f_diff_v1_f_same_v2 < 0 && cp_f_diff_v2_f_same_v2 < 0)) // f_same_v2 is on the same side of both f_diff_v1 and f_diff_v2
                        {
                            if ((cp_f_diff_v1v2 > 0 && cp_f_diff_v1_f_same_v2 > 0) || (cp_f_diff_v1v2 < 0 && cp_f_diff_v1_f_same_v2 < 0)) // f_diff_v2 is the closest side to f_same_v2, so choose matching_ip2
                                rem_ip1 = true;
                            else // f_diff_v1 is the closest side to f_same_vectors, so choose matching_ip1
                                rem_ip2 = true;
                        }
                        else // f_same_v2 is on different sides of f_diff_v1 and f_diff_v2, so choose matching_ip1
                            rem_ip2 = true;
                    }
                    else if (cp_f_diff_v1_f_same_v2 == 0) 
                    {
                        // cp_f_diff_v1_f_same_v2 is a zero vector
                        if ((cp_f_diff_v1_f_same_v1 > 0 && cp_f_diff_v2_f_same_v1 > 0) || (cp_f_diff_v1_f_same_v1 < 0 && cp_f_diff_v2_f_same_v1 < 0)) // f_same_v1 is on the same side of both f_diff_v1 and f_diff_v2
                        {
                            if ((cp_f_diff_v1v2 > 0 && cp_f_diff_v1_f_same_v1 > 0) || (cp_f_diff_v1v2 < 0 && cp_f_diff_v1_f_same_v1 < 0)) // f_diff_v2 is the closest side to f_same_v1, so choose matching_ip2
                                rem_ip1 = true;
                            else // f_diff_v1 is the closest side to f_same vectors, so choose matching_ip1
                                rem_ip2 = true;
                        }
                        else // f_same_v1 is on different sides of f_diff_v1 and f_diff_v2, so choose f_diff_v1
                            rem_ip2 = true;
                    }
                    else if (cp_f_diff_v2_f_same_v1 == 0) 
                    {
                        // cp_f_diff_v2_f_same_v1 is a zero vector
                        if ((cp_f_diff_v1_f_same_v2 > 0 && cp_f_diff_v2_f_same_v2 > 0) || (cp_f_diff_v1_f_same_v2 < 0 && cp_f_diff_v2_f_same_v2 < 0)) // f_same_v2 is on the same side of both f_diff_v1 and f_diff_v2
                        {
                            if ((cp_f_diff_v1v2 > 0 && cp_f_diff_v1_f_same_v2 > 0) || (cp_f_diff_v1v2 < 0 && cp_f_diff_v1_f_same_v2 < 0)) // f_diff_v2 is the closest side to f_same_v2, so choose matching_ip2
                                rem_ip1 = true;
                            else // f_diff_v1 is the closest side to f_same vectors, so choose matching_ip1
                                rem_ip2 = true;
                        }
                        else // f_same_v2 is on different sides of f_diff_v1 and f_diff_v2, so choose f_diff_v2
                            rem_ip1 = true;
                    }
                    else if (cp_f_diff_v2_f_same_v2 == 0) 
                    {
                        // cp_f_diff_v2_f_same_v2 is a zero vector
                        if ((cp_f_diff_v1_f_same_v1 > 0 && cp_f_diff_v2_f_same_v1 > 0) || (cp_f_diff_v1_f_same_v1 < 0 && cp_f_diff_v2_f_same_v1 < 0)) // f_same_v1 is on the same side of both f_diff_v1 and f_diff_v2
                        {
                            if ((cp_f_diff_v1v2 > 0 && cp_f_diff_v1_f_same_v1 > 0) || (cp_f_diff_v1v2 < 0 && cp_f_diff_v1_f_same_v1 < 0)) // f_diff_v2 is the closest side to f_same_v1, so choose matching_ip2
                                rem_ip1 = true;
                            else // f2_v1 is the closest side to f1_vectors, so choose matching_ip1
                                rem_ip2 = true;
                        }
                        else // f_same_v1 is on different sides of f_diff_v1 and f_diff_v2, so choose matching_ip2
                            rem_ip1 = true;
                    }
                    else if ((cp_f_diff_v1_f_same_v1 > 0 && cp_f_diff_v1_f_same_v2 > 0) || (cp_f_diff_v1_f_same_v1 < 0 && cp_f_diff_v1_f_same_v2 < 0)) // both f_same_v1 and f_same_v2 are on the same side of f_diff_v1
                    {
                        if ((cp_f_diff_v2_f_same_v1 > 0 && cp_f_diff_v2_f_same_v2 > 0) || (cp_f_diff_v2_f_same_v1 < 0 && cp_f_diff_v2_f_same_v2 < 0)) // both f_same_v1 and f_same_v2 are on the same side of f_diff_v2
                        {
                            if ((cp_f_diff_v1v2 > 0 && cp_f_diff_v1_f_same_v1 > 0) || (cp_f_diff_v1v2 < 0 && cp_f_diff_v1_f_same_v1 < 0)) // f_diff_v2 is the closest side to f_same vectors, so choose matching_ip2
                                rem_ip1 = true;
                            else // f_diff_v1 is the closest side to f_same vectors, so choose matching_ip1
                                rem_ip2 = true;
                        }
                        else // f_same_v1 and f_same_v2 are on different sides of f_diff_v2, so choose matching_ip2
                            rem_ip1 = true;
                    }
                    else // f_same_v1 and f_same_v2 are on different sides of f_diff_v1
                    {
                        if ((cp_f_diff_v2_f_same_v1 > 0 && cp_f_diff_v2_f_same_v2 > 0) || (cp_f_diff_v2_f_same_v1 < 0 && cp_f_diff_v2_f_same_v2 < 0)) // both f_same_v1 and f_same_v2 are on the same side of f_diff_v2, so choose matching_ip1
                            rem_ip2 = true;
                        else // f_same_v1 and f_same_v2 are on different sides of f_diff_v2
                        {
                            // unable to determine which side to choose, so remove both and use corner point
#ifdef DEBUG_INTERSECT_MESHES_2D_I_PT_LIST
                            cout << "Intersect_Meshes_2D::I_Pt_List::process_matching_i_pts using corner point\n";
#endif
                            rem_ip1 = true;
                            rem_ip2 = true;
                            add_ip = true;
                            pt = common_corner;
                            Location corner_loc = Location::internal;
                            if (ip1_loc == Location::p1p2)
                                corner_loc = ip2_loc == Location::p1p3 ? Location::p1 : Location::p2;
                            else if (ip1_loc == Location::p1p3)
                                corner_loc = ip2_loc == Location::p1p2 ? Location::p1 : Location::p3;
                            else // ip1_loc == Location::p2p3
                                corner_loc = ip2_loc == Location::p1p2 ? Location::p2 : Location::p3;
                            f1_loc = f1_loc_same ? matching_ip1.f1_loc : corner_loc;
                            f2_loc = f1_loc_same ? corner_loc : matching_ip1.f2_loc;
#ifdef DEBUG_INTERSECT_MESHES_2D_I_PT_LIST
                            cout << "Intersect_Meshes_2D::I_Pt_List::process_matching_i_pts add_ip: " << add_ip << " pt x: " << pt->get_x() << " y: " << pt->get_y() << " f1_loc: " << f1_loc << " f2_loc: " << f2_loc << "\n";
#endif
                        }
                    }
                }
            }
        } 
        else // ip1_loc == Location::internal
            rem_ip1 = true; // prefer side or corner points over internal points

        // remove the 'extra' intersect point if it was found
        if (rem_ip2)
        {
#ifdef DEBUG_INTERSECT_MESHES_2D_I_PT_LIST
            cout << "Intersect_Meshes_2D::I_Pt_List::process_matching_i_pts removing ip2\n";
#endif
            i_points.erase(matching_ip2);
        }

        if (add_ip)
        {
#ifdef DEBUG_INTERSECT_MESHES_2D_I_PT_LIST
            cout << "Intersect_Meshes_2D::I_Pt_List::process_matching_i_pts adding ip x: " << pt->get_x() << " y: " << pt->get_y() << " f1_loc: " << f1_loc << " f2_loc: " << f2_loc << "\n";
#endif
            Intersect_Point i_pt(pt, f1_loc, f2_loc);
            if (i_points.end() == find(i_points.begin(), i_points.end(), i_pt))
                i_points.push_back(i_pt);
#ifdef DEBUG_INTERSECT_MESHES_2D_I_PT_LIST
            else
                cout << "Intersect_Meshes_2D::I_Pt_List::process_matching_i_pts i_pt already exists\n";
#endif
        }
        
        if (rem_ip1)
        {
#ifdef DEBUG_INTERSECT_MESHES_2D_I_PT_LIST
            cout << "Intersect_Meshes_2D::I_Pt_List::process_matching_i_pts removing ip1\n";
#endif
            const_iterator it = find(i_points.begin(), i_points.end(), matching_ip1);
            if (it == end())
                throw runtime_error("Unable to locate existing intersect point");
            return i_points.erase(it);
        }
        
        const_iterator it = find(i_points.begin(), i_points.end(), matching_ip1);
        if (it == end())
            throw runtime_error("Unable to locate existing intersect point");
        return ++it;
    }
    
    const bool Intersect_Meshes_2D::I_Pt_List::validate(const Facet_2D& facet1, 
            const Facet_2D& facet2)
    {
#ifdef DEBUG_INTERSECT_MESHES_2D_I_PT_LIST
        cout << "Intersect_Meshes_2D::I_Pt_List::validate begin\n";
#endif
        // check if the point matches or is considered is_equal with another intersect point
        // and remove duplicate if possible
        const_iterator it = i_points.begin();
        while (it != i_points.end())
        {
            Intersect_Point ip(*it);
#ifdef DEBUG_INTERSECT_MESHES_2D_I_PT_LIST
            cout << "Intersect_Meshes_2D::I_Pt_List::validate checking intersect pt x: " << ip.pt->get_x() << " y: " << ip.pt->get_y() << " (" << ip.pt << ") f1_loc: " << ip.f1_loc << " f2_loc: " << ip.f2_loc << "\n";
#endif
            // increment it in the for loop definition, but update it in methods if necessary
            for (const_iterator it2 = ++it; it2 != i_points.end(); ++it2)
            {
#ifdef DEBUG_INTERSECT_MESHES_2D_I_PT_LIST
                cout << "Intersect_Meshes_2D::I_Pt_List::validate checking against i pt x: " << (*it2).pt->get_x() << " y: " << (*it2).pt->get_y() << " (" << (*it2).pt << ") f1_loc: " << (*it2).f1_loc << " f2_loc: " << (*it2).f2_loc << "\n";
#endif
                if (matches(*ip.pt, *(*it2).pt))
                {
#ifdef DEBUG_INTERSECT_MESHES_2D_I_PT_LIST
                    cout << "Intersect_Meshes_2D::I_Pt_List::validate pts match\n";
#endif
                    // points match
                    it = process_matching_i_pts(ip, it2, facet1, facet2);
                    break;
                }
            }
        }
        
#ifdef DEBUG_INTERSECT_MESHES_2D_I_PT_LIST
        cout << "Intersect_Meshes_2D::I_Pt_List::validate end\n";
#endif
    }
    
    void Intersect_Meshes_2D::I_Pt_List::push_back(const Intersect_Point& ip)
    {
        const_iterator it = find(i_points.begin(), i_points.end(), ip);
        if (it != i_points.end())
            return;
        
        i_points.push_back(ip);
    }
    
    Intersect_Meshes_2D::I_Pt_Locator::I_Pt_Locator(const Point_2D::Measurement prec) : precision(prec), 
            facet1(), facet2(), generated_pts() {}

    Intersect_Meshes_2D::I_Pt_Locator::I_Pt_Data::I_Pt_Data() : num(0), 
            ip1(), ip2() {}
    
    Intersect_Meshes_2D::I_Pt_Locator::Point_Find::Point_Find(const shared_ptr<Point_2D>& point, 
            const Point_2D::Measurement prec) : precision(prec), pt(point) {}
    
    const bool Intersect_Meshes_2D::I_Pt_Locator::Point_Find::operator()(const shared_ptr<Point_2D>& other_pt)
    {
        return is_equal(*pt, *other_pt, precision);
    }
    
    const bool Intersect_Meshes_2D::I_Pt_Locator::matches(const shared_ptr<Point_2D>& p1, 
            const shared_ptr<Point_2D>& p2) const
    {
        return p1->get_x() == p2->get_x() && p1->get_y() == p2->get_y();
    }
    
    void Intersect_Meshes_2D::I_Pt_Locator::side_i_pt_loc(const shared_ptr<Point_2D> side_start, 
            const shared_ptr<Point_2D> side_end, const Location side_loc, 
            shared_ptr<Point_2D>& i_pt, Location& loc)
    {
        if (is_equal(*i_pt, *side_start, precision))
        {
            i_pt = side_start;
            loc = (side_loc == Location::p2p3) ? Location::p2 : Location::p1;
        }
        else if (is_equal(*i_pt, *side_end, precision))
        {
            i_pt = side_end;
            loc = (side_loc == Location::p1p2) ? Location::p2 : Location::p3;
        }
        else
        {
            loc = side_loc;
            // check existing points for pt value
            vector<shared_ptr<Point_2D>>::const_iterator it = find_if(generated_pts.begin(), generated_pts.end(), Point_Find(i_pt, precision));
            if (it != generated_pts.end())
                i_pt = *it;
            else
                generated_pts.push_back(i_pt);
        }
    }

    const bool Intersect_Meshes_2D::I_Pt_Locator::intersect_sides(
            const shared_ptr<Point_2D> f1_side_start, 
            const shared_ptr<Point_2D> f1_side_end, const Location f1_side_loc, 
            const shared_ptr<Point_2D> f2_side_start, 
            const shared_ptr<Point_2D> f2_side_end, const Location f2_side_loc, 
            I_Pt_Data& i_pt_data)
    {
#ifdef DEBUG_INTERSECT_MESHES_2D_I_PT_LOCATOR
        cout << "Intersect_Meshes_2D::I_Pt_Locator::intersect_sides begin\n";
#endif
        Vector_2D_idata idata;
        if (intersect_vectors(*f1_side_start, *f1_side_end, *f2_side_start, *f2_side_end, idata, precision))
        {
#ifdef DEBUG_INTERSECT_MESHES_2D_I_PT_LOCATOR
            cout << "Intersect_Meshes_2D::I_Pt_Locator::intersect_sides found " << idata.num << " intersect points\n";
#endif
            shared_ptr<Point_2D> i_pt(new Point_2D(idata.p1));
            side_i_pt_loc(f1_side_start, f1_side_end, f1_side_loc, i_pt, i_pt_data.ip1.f1_loc);
            side_i_pt_loc(f2_side_start, f2_side_end, f2_side_loc, i_pt, i_pt_data.ip1.f2_loc);
            i_pt_data.ip1.pt = i_pt;
            ++i_pt_data.num;
#ifdef DEBUG_INTERSECT_MESHES_2D_I_PT_LOCATOR
            cout << "Intersect_Meshes_2D::I_Pt_Locator::intersect_sides ip1.pt x: " << i_pt_data.ip1.pt->get_x() << " y: " << i_pt_data.ip1.pt->get_y() << " f1_loc: " << i_pt_data.ip1.f1_loc << " f2_loc: " << i_pt_data.ip1.f2_loc << "\n";
#endif
            if (idata.num == 2)
            {
                i_pt = shared_ptr<Point_2D>(new Point_2D(idata.p2));
                side_i_pt_loc(f1_side_start, f1_side_end, f1_side_loc, i_pt, i_pt_data.ip2.f1_loc);
                side_i_pt_loc(f2_side_start, f2_side_end, f2_side_loc, i_pt, i_pt_data.ip2.f2_loc);
                i_pt_data.ip2.pt = i_pt;
                ++i_pt_data.num;
#ifdef DEBUG_INTERSECT_MESHES_2D_I_PT_LOCATOR
                cout << "Intersect_Meshes_2D::I_Pt_Locator::intersect_sides ip2.pt x: " << i_pt_data.ip2.pt->get_x() << " y: " << i_pt_data.ip2.pt->get_y() << " f1_loc: " << i_pt_data.ip2.f1_loc << " f2_loc: " << i_pt_data.ip2.f2_loc << "\n";
#endif
            }
        }
        return idata.num > 0;
    }
    
    void Intersect_Meshes_2D::I_Pt_Locator::intersect_f1_side_to_f2(
            const Location f1_side, I_Pt_List& intersect_points)
    {
        const shared_ptr<Point_2D> v1_start(f1_side == Location::p2p3 ? facet1.get_point2() : facet1.get_point1());
        const shared_ptr<Point_2D> v1_end(f1_side == Location::p1p2 ? facet1.get_point2() : facet1.get_point3());
        
        I_Pt_Data i_pt_data;
#ifdef DEBUG_INTERSECT_MESHES_2D_I_PT_LOCATOR
        cout << "Intersect_Meshes_2D::I_Pt_Locator::intersect_f1_side_to_f2 p1p2\n";
#endif
        if (intersect_sides(v1_start, v1_end, f1_side, facet2.get_point1(), facet2.get_point2(), Location::p1p2, i_pt_data))
        {
            intersect_points.push_back(i_pt_data.ip1);
            if (i_pt_data.num == 2)
            {
                intersect_points.push_back(i_pt_data.ip2);
                return; // no need to process further, found maximum of two intersection points
            }
        }
        
        I_Pt_Data t_i_pt_data;
#ifdef DEBUG_INTERSECT_MESHES_2D_I_PT_LOCATOR
        cout << "Intersect_Meshes_2D::I_Pt_Locator::intersect_f1_side_to_f2 p1p3\n";
#endif
        if (intersect_sides(v1_start, v1_end, f1_side, facet2.get_point1(), facet2.get_point3(), Location::p1p3, t_i_pt_data))
        {
//            if (i_pt_data.num == 0 || !is_equal(*i_pt_data.ip1.pt, *t_i_pt_data.ip1.pt, precision))
            if (i_pt_data.num == 0 || !matches(i_pt_data.ip1.pt, t_i_pt_data.ip1.pt))
            {
                intersect_points.push_back(t_i_pt_data.ip1);
                switch (i_pt_data.num) {
                    case (0):
                        i_pt_data.num = 1;
                        i_pt_data.ip1 = t_i_pt_data.ip1;
                        break;
                    case (1):
                        return; // found maximum of two intersect points, do not do any more processing
                }
            }
            if (t_i_pt_data.num == 2) // should not need to check point against i_pt_data.ip1
            {
                intersect_points.push_back(t_i_pt_data.ip2);
                return; // found maximum of two intersect points, do not do any more processing
            }
        }
        
#ifdef DEBUG_INTERSECT_MESHES_2D_I_PT_LOCATOR
        cout << "Intersect_Meshes_2D::I_Pt_Locator::intersect_f1_side_to_f2 p2p3\n";
#endif
        t_i_pt_data.num = 0;
        if (intersect_sides(v1_start, v1_end, f1_side, facet2.get_point2(), facet2.get_point3(), Location::p2p3, t_i_pt_data))
        {
//            if (i_pt_data.num == 0 || !is_equal(*i_pt_data.ip1.pt, *t_i_pt_data.ip1.pt, precision))
            if (i_pt_data.num == 0 || !matches(i_pt_data.ip1.pt, t_i_pt_data.ip1.pt))
            {
                intersect_points.push_back(t_i_pt_data.ip1);
                switch (i_pt_data.num) {
                    case (0):
                        i_pt_data.num = 1;
                        i_pt_data.ip1 = t_i_pt_data.ip1;
                        break;
                    case (1):
                        return; // found maximum of two intersect points, do not do any more processing
                }
            }
            if (t_i_pt_data.num == 2) // should not need to check point against i_pt_data.ip1
            {
                intersect_points.push_back(t_i_pt_data.ip2);
                return; // found maximum of two intersect points, do not do any more processing
            }
        }
        
        if (i_pt_data.num == 0)
        {
            bool pt_on_side(false);
            if (facet2.contains_point(*v1_start, pt_on_side, precision) && facet2.contains_point(*v1_end, pt_on_side, precision))
            {
#ifdef DEBUG_INTERSECT_MESHES_2D_I_PT_LOCATOR
                cout << "Intersect_Meshes_2D::I_Pt_Locator::intersect_f1_side_to_f2 contains v1_start and v1_end\n";
#endif
                intersect_points.push_back(Intersect_Point(v1_start, f1_side == Location::p2p3 ? Location::p2 : Location::p1, Location::internal));
                intersect_points.push_back(Intersect_Point(v1_end, f1_side == Location::p1p2 ? Location::p2 : Location::p3, Location::internal));
            }
//            else
//            {
//                // try to intersect with facet plane
//                Vector_2D i_vector(*v1_start, *v1_end);
//                Point_2D i_point(0,0);
//#ifdef DEBUG_INTERSECT_MESHES_2D_I_PT_LOCATOR
//                cout << "Intersect_Meshes_2D::I_Pt_Locator::intersect_f1_side_to_f2 intersect_line_facet_plane\n";
//#endif
//                if (intersect_line_facet_plane(i_vector, *v1_start, facet2, i_point, precision) &&
//                        is_pt_on_vector(i_point, *v1_start, *v1_end, precision) && 
//                        facet2.contains_point(i_point, pt_on_side, precision))
//                {
//                    Location loc(Location::internal);
//                    shared_ptr<Point_2D> i_pt(new Point_2D(i_point));
//                    side_i_pt_loc(v1_start, v1_end, f1_side, i_pt, loc);
//#ifdef DEBUG_INTERSECT_MESHES_2D_I_PT_LOCATOR
//                    cout << "Intersect_Meshes_2D::I_Pt_Locator::intersect_f1_side_to_f2 adding intersection point from intersect_line_facet_plane ip x: " << i_pt->get_x() << " y: " << i_pt->get_y() << " z: " << i_pt->get_z() << " f1_loc: " << loc << " f2_loc: internal\n";
//#endif
//                    intersect_points.push_back(Intersect_Point(i_pt, loc, Location::internal));
//                }
//            }
        }
        else // one intersect point
        {
            bool pt_on_side(false);
            // check if either end is inside the facet
//            if (!is_equal(*v1_start, *i_pt_data.ip1.pt, precision) && facet2.contains_point(*v1_start, pt_on_side, precision))
            if (i_pt_data.ip1.f1_loc != (f1_side == Location::p2p3 ? Location::p2 : Location::p1) && facet2.contains_point(*v1_start, pt_on_side, precision))
            {
#ifdef DEBUG_INTERSECT_MESHES_2D_I_PT_LOCATOR
                cout << "Intersect_Meshes_2D::I_Pt_Locator::intersect_f1_side_to_f2 adding v1_start intersection point\n";
#endif
                intersect_points.push_back(Intersect_Point(v1_start, f1_side == Location::p2p3 ? Location::p2 : Location::p1, Location::internal));
            }
//            else if (!is_equal(*v1_end, *i_pt_data.ip1.pt, precision) && facet2.contains_point(*v1_end, pt_on_side, precision))
            else if (i_pt_data.ip1.f1_loc != (f1_side == Location::p1p2 ? Location::p2 : Location::p3) && facet2.contains_point(*v1_end, pt_on_side, precision))
            {
#ifdef DEBUG_INTERSECT_MESHES_2D_I_PT_LOCATOR
                cout << "Intersect_Meshes_2D::I_Pt_Locator::intersect_f1_side_to_f2 adding v1_end intersection point\n";
#endif
                intersect_points.push_back(Intersect_Point(v1_end, f1_side == Location::p1p2 ? Location::p2 : Location::p3, Location::internal));
            }
        }
    }
    
    void Intersect_Meshes_2D::I_Pt_Locator::intersect_f2_side_to_f1(
            const Location f2_side, I_Pt_List& intersect_points)
    {
        const shared_ptr<Point_2D> v2_start(f2_side == Location::p2p3 ? facet2.get_point2() : facet2.get_point1());
        const shared_ptr<Point_2D> v2_end(f2_side == Location::p1p2 ? facet2.get_point2() : facet2.get_point3());
        
//        bool found_i_pt(false);
        bool pt_on_side(false);
        if (facet1.contains_point(*v2_start, pt_on_side, precision))
        {
#ifdef DEBUG_INTERSECT_MESHES_2D_I_PT_LOCATOR
            cout << "Intersect_Meshes_2D::I_Pt_Locator::intersect_f2_side_to_f1 adding v2_start pt x: " << v2_start->get_x() << " y: " << v2_start->get_y() << " f1_loc: internal f2_loc: " << (f2_side == Location::p2p3 ? Location::p2 : Location::p1) << "\n";
#endif
            intersect_points.push_back(Intersect_Point(v2_start, Location::internal, f2_side == Location::p2p3 ? Location::p2 : Location::p1));
//            found_i_pt = true;
        }

        if (facet1.contains_point(*v2_end, pt_on_side, precision))
        {
#ifdef DEBUG_INTERSECT_MESHES_2D_I_PT_LOCATOR
            cout << "Intersect_Meshes_2D::I_Pt_Locator::intersect_f2_side_to_f1 adding v2_end pt x: " << v2_end->get_x() << " y: " << v2_end->get_y() << " f1_loc: internal f2_loc: " << (f2_side == Location::p1p2 ? Location::p2 : Location::p3) << "\n";
#endif
            intersect_points.push_back(Intersect_Point(v2_end, Location::internal, f2_side == Location::p1p2 ? Location::p2 : Location::p3));
//            found_i_pt = true;
        }
        
//        if (found_i_pt)
//            return;
//
//        Point_2D i_point(0,0);
//        if (intersect_line_facet_plane(Vector_2D(*v2_start, *v2_end), *v2_start, facet1, i_point, precision) &&
//                is_pt_on_vector(i_point, *v2_start, *v2_end, precision) && 
//                facet1.contains_point(i_point, pt_on_side, precision))
//        {
//            // f1_loc should be internal
//            shared_ptr<Point_2D> i_pt(new Point_2D(i_point));
//            Location loc(f2_side);
//            side_i_pt_loc(v2_start, v2_end, f2_side, i_pt, loc);
//#ifdef DEBUG_INTERSECT_MESHES_2D_I_PT_LOCATOR
//            cout << "Intersect_Meshes_2D::I_Pt_Locator::intersect_f2_side_to_f1 adding i_point pt x: " << i_pt->get_x() << " y: " << i_pt->get_y() << " z: " << i_pt->get_z() << " f1_loc: internal f2_loc: " << loc << "\n";
//#endif
//            intersect_points.push_back(Intersect_Point(i_pt, Location::internal, loc));
//        }
    }
    
    const int Intersect_Meshes_2D::I_Pt_Locator::is_f2_side_intersected(const Location f2_side, 
            const I_Pt_List& intersect_points, bool& corner1_found, bool& corner2_found)
    {
        Location start_corner = f2_side == Location::p2p3 ? Location::p2 : Location::p1;
        Location end_corner = f2_side == Location::p1p2 ? Location::p2 : Location::p3;
        
        int count(0);
        for (I_Pt_List::const_iterator it = intersect_points.begin(); it != intersect_points.end(); ++it)
        {
            Location loc(it->f2_loc);
            if (loc == f2_side)
                ++count;
            else if (loc == start_corner)
            {
                ++count;
                corner1_found = true;
            }
            else if (loc == end_corner)
            {
                ++count;
                corner2_found = true;
            }
        }

        return count;
    }
    
    const bool Intersect_Meshes_2D::I_Pt_Locator::operator()(const Facet_2D& f1, 
            const Facet_2D& f2, I_Pt_List& intersect_points)
    {
        // assign values
        facet1 = f1;
        facet2 = f2;
#ifdef DEBUG_INTERSECT_MESHES_2D_I_PT_LOCATOR
        cout << "Intersect_Meshes_2D::I_Pt_Locator::locate facet1 p1 x: " << facet1.get_point1()->get_x() << 
                " y: " << facet1.get_point1()->get_y() << " p2 x: " << facet1.get_point2()->get_x() << 
                " y: " << facet1.get_point2()->get_y() << " p3 x: " << facet1.get_point3()->get_x() << 
                " y: " << facet1.get_point3()->get_y() << "\n";
        cout << "Intersect_Meshes_2D::I_Pt_Locator::locate facet2 p1 x: " << facet2.get_point1()->get_x() << 
                " y: " << facet2.get_point1()->get_y() << " p2 x: " << facet2.get_point2()->get_x() << 
                " y: " << facet2.get_point2()->get_y() << " p3 x: " << facet2.get_point3()->get_x() << 
                " y: " << facet2.get_point3()->get_y() << "\n";
#endif
        I_Pt_List intersect_pts;
        // intersect vector sides, form internal segments
        // intersect i_p1p2 to facet
#ifdef DEBUG_INTERSECT_MESHES_2D_I_PT_LOCATOR
        cout << "Intersect_Meshes_2D::I_Pt_Locator::locate intersecting i_p1p2 with facet\n";
#endif
        intersect_f1_side_to_f2(Location::p1p2, intersect_pts);
#ifdef DEBUG_INTERSECT_MESHES_2D_I_PT_LOCATOR
        int count = 0;
        for (I_Pt_List::const_iterator it = intersect_pts.begin(); it != intersect_pts.end(); ++it)
        {
            cout << "Intersect_Meshes_2D::I_Pt_Locator::locate intersect_pts[" << count++ << "] x: " << (*it).pt->get_x() << " y: " << (*it).pt->get_y() << " f1_loc: " << (*it).f1_loc << " f2_loc: " << (*it).f2_loc << "\n";
        }
#endif
        
        // intersect i_p1p3 to facet
#ifdef DEBUG_INTERSECT_MESHES_2D_I_PT_LOCATOR
        cout << "Intersect_Meshes_2D::I_Pt_Locator::locate intersecting i_p1p3 with facet. intersect_points size: " << intersect_pts.size() << "\n";
#endif
        intersect_f1_side_to_f2(Location::p1p3, intersect_pts);
#ifdef DEBUG_INTERSECT_MESHES_2D_I_PT_LOCATOR
        count = 0;
        for (I_Pt_List::const_iterator it = intersect_pts.begin(); it != intersect_pts.end(); ++it)
        {
            cout << "Intersect_Meshes_2D::I_Pt_Locator::locate intersect_pts[" << count++ << "] x: " << (*it).pt->get_x() << " y: " << (*it).pt->get_y() << " f1_loc: " << (*it).f1_loc << " f2_loc: " << (*it).f2_loc << "\n";
        }
#endif
        
        // intersect i_p2p3 to facet
#ifdef DEBUG_INTERSECT_MESHES_2D_I_PT_LOCATOR
        cout << "Intersect_Meshes_2D::I_Pt_Locator::locate intersecting i_p2p3 with facet. intersect_points size: " << intersect_pts.size() << "\n";
#endif
        intersect_f1_side_to_f2(Location::p2p3, intersect_pts);
#ifdef DEBUG_INTERSECT_MESHES_2D_I_PT_LOCATOR
        count = 0;
        for (I_Pt_List::const_iterator it = intersect_pts.begin(); it != intersect_pts.end(); ++it)
        {
            cout << "Intersect_Meshes_2D::I_Pt_Locator::locate intersect_pts[" << count++ << "] x: " << (*it).pt->get_x() << " y: " << (*it).pt->get_y() << " f1_loc: " << (*it).f1_loc << " f2_loc: " << (*it).f2_loc << "\n";
        }
#endif
        
        // if a facet side was not intersected, see if it intersects the intersecting_facet
        bool corner1_found(false);
        bool corner2_found(false);
        int intersection_ct(is_f2_side_intersected(Location::p1p2, intersect_pts, corner1_found, corner2_found));
        if (intersection_ct == 0)
        {
#ifdef DEBUG_INTERSECT_MESHES_2D_I_PT_LOCATOR
            cout << "Intersect_Meshes_2D::I_Pt_Locator::locate Intersecting f2 p1p2 side\n";
#endif
            intersect_f2_side_to_f1(Location::p1p2, intersect_pts);
        }
        else if (intersection_ct == 1)
        {
#ifdef DEBUG_INTERSECT_MESHES_2D_I_PT_LOCATOR
            cout << "Intersect_Meshes_2D::I_Pt_Locator::locate checking f2.p1 and f2.p2\n";
#endif
            bool pt_on_side(false);
            if (!corner1_found && facet1.contains_point(*facet2.get_point1(), pt_on_side, precision))
            {
#ifdef DEBUG_INTERSECT_MESHES_2D_I_PT_LOCATOR
                cout << "Intersect_Meshes_2D::I_Pt_Locator::locate adding f2.p1 pt x: " << facet2.get_point1()->get_x() << " y: " << facet2.get_point1()->get_y() << " f1_loc: internal f2_loc: p1\n";
#endif
                intersect_pts.push_back(Intersect_Point(facet2.get_point1(), Location::internal, Location::p1));
            }

            if (!corner2_found && facet1.contains_point(*facet2.get_point2(), pt_on_side, precision))
            {
#ifdef DEBUG_INTERSECT_MESHES_2D_I_PT_LOCATOR
                cout << "Intersect_Meshes_2D::I_Pt_Locator::locate adding f2.p2 pt x: " << facet2.get_point2()->get_x() << " y: " << facet2.get_point2()->get_y() << " f1_loc: internal f2_loc: p2\n";
#endif
                intersect_pts.push_back(Intersect_Point(facet2.get_point2(), Location::internal, Location::p2));
            }
        }
        
        corner1_found = false;
        corner2_found = false;
        intersection_ct = is_f2_side_intersected(Location::p1p3, intersect_pts, corner1_found, corner2_found);
        if (intersection_ct == 0)
        {
#ifdef DEBUG_INTERSECT_MESHES_2D_I_PT_LOCATOR
            cout << "Intersect_Meshes_2D::I_Pt_Locator::locate Intersecting f2 p1p3 side\n";
#endif
            intersect_f2_side_to_f1(Location::p1p3, intersect_pts);
        }
        else if (intersection_ct == 1)
        {
#ifdef DEBUG_INTERSECT_MESHES_2D_I_PT_LOCATOR
            cout << "Intersect_Meshes_2D::I_Pt_Locator::locate checking f2.p3\n";
#endif
            bool pt_on_side(false);
            if (!corner2_found && facet1.contains_point(*facet2.get_point3(), pt_on_side, precision))
            {
#ifdef DEBUG_INTERSECT_MESHES_2D_I_PT_LOCATOR
                cout << "Intersect_Meshes_2D::I_Pt_Locator::locate adding f2.p3 pt x: " << facet2.get_point3()->get_x() << " y: " << facet2.get_point3()->get_y() << " f1_loc: internal f2_loc: p3\n";
#endif
                intersect_pts.push_back(Intersect_Point(facet2.get_point3(), Location::internal, Location::p3));
            }
        }
        
        corner1_found = false;
        corner2_found = false;
        intersection_ct = is_f2_side_intersected(Location::p2p3, intersect_pts, corner1_found, corner2_found);
        if (intersection_ct == 0)
        {
#ifdef DEBUG_INTERSECT_MESHES_2D_I_PT_LOCATOR
            cout << "Intersect_Meshes_2D::I_Pt_Locator::locate Intersecting f2 p2p3 side\n";
#endif
            intersect_f2_side_to_f1(Location::p2p3, intersect_pts);
        }
#ifdef DEBUG_INTERSECT_MESHES_2D_I_PT_LOCATOR
        cout << "Intersect_Meshes_2D::I_Pt_Locator::locate after f2 to f1 intersect_pts size: " << intersect_pts.size() << "\n";
#endif
        if (intersect_pts.empty())
            return false;

        // validate intersect points
        intersect_pts.validate(facet1, facet2);
        
        for (I_Pt_List::const_iterator it = intersect_pts.begin(); it != intersect_pts.end(); ++it)
            intersect_points.push_back(*it);
        return true;
    }
    
    Intersect_Meshes_2D::Facets::Facets() : point_list(), facet_list() {}
    
    Intersect_Meshes_2D::Facets::Facets(const Mesh_2D& mesh) : point_list(), facet_list() 
    {
#ifdef DEBUG_INTERSECT_MESHES_2D_FACETS
        cout << "Intersect_Meshes_2D::Facets begin\n";
#endif
        for (Mesh_2D::const_iterator it = mesh.begin(); it != mesh.end(); ++it)
        {
#ifdef DEBUG_INTERSECT_MESHES_2D_FACETS
            cout << "Intersect_Meshes_2D::Facets adding facet p1 x: " << it->get_point1()->get_x() << 
                    " y: " << it->get_point1()->get_y() << " p2 x: " << it->get_point2()->get_x() << 
                    " y: " << it->get_point2()->get_y() << " p3 x: " << it->get_point3()->get_x() << 
                    " y: " << it->get_point3()->get_y() << "\n";
#endif
            int p1_index(-1);
            int p2_index(-1);
            int p3_index(-1);
            int index(0);
            for (vector<shared_ptr<Point_2D>>::const_iterator p_it = point_list.begin(); p_it != point_list.end(); ++p_it)
            {
                if (matches(*it->get_point1(), **p_it))
                    p1_index = index;
                else if (matches(*it->get_point2(), **p_it))
                    p2_index = index;
                else if (matches(*it->get_point3(), **p_it))
                    p3_index = index;
                ++index;
            }
            if (p1_index == -1)
            {
                point_list.push_back(it->get_point1());
                p1_index = index;
                ++index;
            }
            if (p2_index == -1)
            {
                point_list.push_back(it->get_point2());
                p2_index = index;
                ++index;
            }
            if (p3_index == -1)
            {
                point_list.push_back(it->get_point3());
                p3_index = index;
            }
            facet_list.push_back(Facet(p1_index, p2_index, p3_index));
        }
#ifdef DEBUG_INTERSECT_MESHES_2D_FACETS
        cout << "Intersect_Meshes_2D::Facets end\n";
#endif
    }

    void Intersect_Meshes_2D::Facets::clear()
    {
        facet_list.clear();
        point_list.clear();
    }
    
    const bool Intersect_Meshes_2D::Facets::matches(const Point_2D& p1, const Point_2D& p2) const
    {
        return p1.get_x() == p2.get_x() && p1.get_y() == p2.get_y();
    }
    
    const shared_ptr<Point_2D> Intersect_Meshes_2D::Facets::get_point(int index) const
    {
        return *(point_list.begin() + index);
    }
    
    void Intersect_Meshes_2D::Facets::push_back(const shared_ptr<Point_2D>& p1, 
            const shared_ptr<Point_2D>& p2, const shared_ptr<Point_2D>& p3)
    {
#ifdef DEBUG_INTERSECT_MESHES_2D_FACETS
        cout << "Intersect_Meshes_2D::Facets::push_back(p1, p2, p3) begin\n";
        cout << "Intersect_Meshes_2D::Facets::push_back p1 x: " << p1->get_x() << " y: " << 
                p1->get_y() << " p2 x: " << p2->get_x() << " y: " << p2->get_y() << 
                " p3 x: " << p3->get_x() << " y: " << p3->get_y() << "\n";
#endif
        int p1_index(-1);
        int p2_index(-1);
        int p3_index(-1);
        int index(0);
        for (vector<shared_ptr<Point_2D>>::const_iterator p_it = point_list.begin(); p_it != point_list.end(); ++p_it)
        {
            if (matches(*p1, **p_it))
                p1_index = index;
            else if (matches(*p2, **p_it))
                p2_index = index;
            else if (matches(*p3, **p_it))
                p3_index = index;
            ++index;
        }
        if (p1_index == -1)
        {
            point_list.push_back(p1);
            p1_index = index;
            ++index;
        }
        if (p2_index == -1)
        {
            point_list.push_back(p2);
            p2_index = index;
            ++index;
        }
        if (p3_index == -1)
        {
            point_list.push_back(p3);
            p3_index = index;
        }
#ifdef DEBUG_INTERSECT_MESHES_2D_FACETS
        cout << "Intersect_Meshes_2D::Facets::push_back(p1, p2, p3) adding facet p1: " << p1_index << " p2: " << p2_index << " p3: " << p3_index << "\n";
#endif
        facet_list.push_back(Facet(p1_index, p2_index, p3_index));
    }
    
    void Intersect_Meshes_2D::Facets::push_back(const Facet_2D& facet)
    {
#ifdef DEBUG_INTERSECT_MESHES_2D_FACETS
        cout << "Intersect_Meshes_2D::Facets::push_back(facet) begin\n";
        cout << "Intersect_Meshes_2D::Facets::push_back p1 x: " << facet.get_point1()->get_x() << 
                " y: " << facet.get_point1()->get_y() << " p2 x: " << facet.get_point2()->get_x() << 
                " y: " << facet.get_point2()->get_y() << " p3 x: " << facet.get_point3()->get_x() << 
                " y: " << facet.get_point3()->get_y() << "\n";
#endif
        int p1_index(-1);
        int p2_index(-1);
        int p3_index(-1);
        int index(0);
        for (vector<shared_ptr<Point_2D>>::const_iterator p_it = point_list.begin(); p_it != point_list.end(); ++p_it)
        {
            if (matches(*facet.get_point1(), **p_it))
                p1_index = index;
            else if (matches(*facet.get_point2(), **p_it))
                p2_index = index;
            else if (matches(*facet.get_point3(), **p_it))
                p3_index = index;
            ++index;
        }
        if (p1_index == -1)
        {
            point_list.push_back(facet.get_point1());
            p1_index = index;
            ++index;
        }
        if (p2_index == -1)
        {
            point_list.push_back(facet.get_point2());
            p2_index = index;
            ++index;
        }
        if (p3_index == -1)
        {
            point_list.push_back(facet.get_point3());
            p3_index = index;
        }
#ifdef DEBUG_INTERSECT_MESHES_2D_FACETS
        cout << "Intersect_Meshes_2D::Facets::push_back(facet) adding facet p1: " << p1_index << " p2: " << p2_index << " p3: " << p3_index << "\n";
#endif
        facet_list.push_back(Facet(p1_index, p2_index, p3_index));
    }
    
    const Facet_2D Intersect_Meshes_2D::Facets::pop()
    {
        Facet last_facet = *(--facet_list.end());
        // do not check for points to remove and update other facet point indices.  
        // Leave unused points in list to reduce processing time
        facet_list.pop_back();
        
        return Facet_2D(*(point_list.begin() + last_facet.get_p1_index()), 
                *(point_list.begin() + last_facet.get_p2_index()),
                *(point_list.begin() + last_facet.get_p3_index()));
    }
    
    void Intersect_Meshes_2D::Facets::replace_all(const Facets& facets)
    {
        facet_list.clear();
        point_list.clear();
        point_list.insert(point_list.begin(), facets.point_list.begin(),facets.point_list.end());
        facet_list.insert(facet_list.begin(), facets.facet_list.begin(),facets.facet_list.end());
    }
    
    const bool Intersect_Meshes_2D::Facets::contains(const Facet_2D& facet) const
    {
        int p1_index(-1);
        int p2_index(-1);
        int p3_index(-1);
        
        int index(0);
        for (vector<shared_ptr<Point_2D>>::const_iterator it = point_list.begin(); it != point_list.end(); ++it)
        {
            if (*it == facet.get_point1())
                p1_index = index;
            else if (*it == facet.get_point2())
                p2_index = index;
            else if (*it == facet.get_point3())
                p3_index = index;
            
            ++index;
        }
        
        if (p1_index == -1 || p2_index == -1 || p3_index == -1)
            return false;
        Facet f(p1_index, p2_index, p3_index);
        return facet_list.end() != find(facet_list.begin(), facet_list.end(), f);
    }
    
    Intersect_Meshes_2D::Facets::const_iterator Intersect_Meshes_2D::Facets::replace_facet(
            const Facet facet, const Facets& new_facets)
    {
#ifdef DEBUG_INTERSECT_MESHES_2D_FACETS
        cout << "Intersect_Meshes_2D::Facets::replace_facet begin\n";
        cout << "Intersect_Meshes_2D::Facets::replace_facet facet p1(" << facet.get_p1_index() << ") x: " << 
                this->get_point(facet.get_p1_index())->get_x() << " y: " << this->get_point(facet.get_p1_index())->get_y() << 
                " p2(" << facet.get_p2_index() << ") x: " << this->get_point(facet.get_p2_index())->get_x() << " y: " << 
                this->get_point(facet.get_p2_index())->get_y() << " p3(" << facet.get_p3_index() << ") x: " << 
                this->get_point(facet.get_p3_index())->get_x() << " y: " << this->get_point(facet.get_p3_index())->get_y() << "\n"; 
        int count = 0;
        for (vector<Facet>::const_iterator it = facet_list.begin(); it != facet_list.end(); ++it)
        {
            cout << "Intersect_Meshes_2D::Facets::replace_facet before facets[" << count << "] p1: " << it->get_p1_index() << 
                    " p2: " << it->get_p2_index() << " p3: " << it->get_p3_index() << "\n";
            ++count;
        }
#endif        
        vector<Facet>::const_iterator f_loc = find(facet_list.begin(), facet_list.end(), facet);
        
        if (f_loc == facet_list.end())
            throw runtime_error("Unable to locate facet");
        
        // remove Facet at f_loc
#ifdef DEBUG_INTERSECT_MESHES_2D_FACETS
        cout << "Intersect_Meshes_2D::Facets::replace_facet f_loc p1: " << f_loc->get_p1_index() << " p2: " << 
                f_loc->get_p2_index() << " p3: " << f_loc->get_p3_index() << "\n";
        cout << "Intersect_Meshes_2D::Facets::replace_facet erasing f_loc facets size: " << facet_list.size() << "\n";
#endif
            
        f_loc = facet_list.erase(f_loc); // f_loc now points to the facet just after the one removed

#ifdef DEBUG_INTERSECT_MESHES_2D_FACETS
        cout << "Intersect_Meshes_2D::Facets::replace_facet erased f_loc facets size: " << facet_list.size() << "\n";
        count = 0;
        for (vector<Facet>::const_iterator it = facet_list.begin(); it != facet_list.end(); ++it)
        {
            cout << "Intersect_Meshes_2D::Facets::replace_facet after erase facets[" << count << "] p1: " << it->get_p1_index() << 
                    " p2: " << it->get_p2_index() << " p3: " << it->get_p3_index() << "\n";
            ++count;
        }
        if (f_loc != facet_list.end())
            cout << "Intersect_Meshes_2D::Facets::replace_facet f_loc p1: " << f_loc->get_p1_index() << 
                    " p2: " << f_loc->get_p2_index() << " p3: " << f_loc->get_p3_index() << "\n";
        else
            cout << "Intersect_Meshes_2D::Facets::replace_facet f_loc == facets.end\n";
#endif
        
        for (Facets::const_iterator it = --new_facets.end(); it != new_facets.begin(); --it)
        {
            int p1_index(-1);
            int p2_index(-1);
            int p3_index(-1);
            
            int index(0);
            for (vector<shared_ptr<Point_2D>>::const_iterator p_it = point_list.begin(); p_it != point_list.end(); ++p_it)
            {
                if (matches(*new_facets.get_point(it->get_p1_index()), **p_it))
                    p1_index = index;
                else if (matches(*new_facets.get_point(it->get_p2_index()), **p_it))
                    p2_index = index;
                else if (matches(*new_facets.get_point(it->get_p3_index()), **p_it))
                    p3_index = index;
                ++index;
            }

            if (p1_index == -1)
            {
                point_list.push_back(new_facets.get_point(it->get_p1_index()));
                p1_index = index;
                ++index;
            }
            if (p2_index == -1)
            {
                point_list.push_back(new_facets.get_point(it->get_p2_index()));
                p2_index = index;
                ++index;
            }
            if (p3_index == -1)
            {
                point_list.push_back(new_facets.get_point(it->get_p3_index()));
                p3_index = index;
            }
#ifdef DEBUG_INTERSECT_MESHES_2D_FACETS
            cout << "Intersect_Meshes_2D::Facets::replace_facet inserting facet p1: " << p1_index << " p2: " << p2_index << " p3: " << p3_index << "\n";
#endif
            f_loc = facet_list.insert(f_loc, Facet(p1_index, p2_index, p3_index));
        }

        int p1_index(-1);
        int p2_index(-1);
        int p3_index(-1);

        int index(0);
        for (vector<shared_ptr<Point_2D>>::const_iterator p_it = point_list.begin(); p_it != point_list.end(); ++p_it)
        {
            if (matches(*new_facets.get_point(new_facets.begin()->get_p1_index()), **p_it))
                p1_index = index;
            else if (matches(*new_facets.get_point(new_facets.begin()->get_p2_index()), **p_it))
                p2_index = index;
            else if (matches(*new_facets.get_point(new_facets.begin()->get_p3_index()), **p_it))
                p3_index = index;
            ++index;
        }

        if (p1_index == -1)
        {
            point_list.push_back(new_facets.get_point(new_facets.begin()->get_p1_index()));
            p1_index = index;
            ++index;
        }
        if (p2_index == -1)
        {
            point_list.push_back(new_facets.get_point(new_facets.begin()->get_p2_index()));
            p2_index = index;
            ++index;
        }
        if (p3_index == -1)
        {
            point_list.push_back(new_facets.get_point(new_facets.begin()->get_p3_index()));
            p3_index = index;
        }
#ifdef DEBUG_INTERSECT_MESHES_2D_FACETS
        cout << "Intersect_Meshes_2D::Facets::replace_facet inserting facet p1: " << p1_index << " p2: " << p2_index << " p3: " << p3_index << "\n";
#endif
        f_loc = facet_list.insert(f_loc, Facet(p1_index, p2_index, p3_index));
        
#ifdef DEBUG_INTERSECT_MESHES_2D_FACETS
        count = 0;
        for (vector<Facet>::const_iterator it = facet_list.begin(); it != facet_list.end(); ++it)
        {
            cout << "Intersect_Meshes_2D::Facets::replace_facet after facets[" << count << "] p1: " << it->get_p1_index() << 
                    " p2: " << it->get_p2_index() << " p3: " << it->get_p3_index() << "\n";
            ++count;
        }
        
        if (f_loc != facet_list.end())
            cout << "Intersect_Meshes_2D::Facets::replace_facet f_loc p1: " << f_loc->get_p1_index() << 
                    " p2: " << f_loc->get_p2_index() << " p3: " << f_loc->get_p3_index() << "\n";
        else
            cout << "Intersect_Meshes_2D::Facets::replace_facet f_loc == facets.end\n";
        cout << "Intersect_Meshes_2D::Facets::replace_facet end. returning f_loc\n";
#endif
        return f_loc;
    }
    
    const Facet Intersect_Meshes_2D::Facets::find_facet(const Facet_2D& facet) const
    {
#ifdef DEBUG_INTERSECT_MESHES_2D_FACETS
        cout << "Intersect_Meshes_2D::Facets::find_facet begin\n";
        cout << "Intersect_Meshes_2D::Facets::find_facet facet p1 x: " << facet.get_point1()->get_x() << 
                " y: " << facet.get_point1()->get_y() << " p2 x: " << facet.get_point2()->get_x() << 
                " y: " << facet.get_point2()->get_y() << " p3 x: " << facet.get_point3()->get_x() << 
                " y: " << facet.get_point3()->get_y() << "\n";
#endif
        int p1_index(-1);
        int p2_index(-1);
        int p3_index(-1);
        
        int index(0);
        for (vector<shared_ptr<Point_2D>>::const_iterator it = point_list.begin(); it != point_list.end(); ++it)
        {
            if (matches(*facet.get_point1(), **it))
                p1_index = index;
            else if (matches(*facet.get_point2(), **it))
                p2_index = index;
            else if (matches(*facet.get_point3(), **it))
                p3_index = index;
            ++index;
        }
        
        if (p1_index == -1 || p2_index == -1 || p3_index == -1)
        {
#ifdef DEBUG_INTERSECT_MESHES_2D_FACETS
            cout << "Intersect_Meshes_2D::Facets::find_facet unable to locate point index p1_index: " << 
                    p1_index << " p2_index: " << p2_index << " p3_index: " << p3_index << "\n";
#endif
            throw runtime_error("Unable to locate facet");
        }
        
        Facet f(p1_index, p2_index, p3_index);
        if (facet_list.end() == find(facet_list.begin(), facet_list.end(), f))
        {
#ifdef DEBUG_INTERSECT_MESHES_2D_FACETS
            cout << "Intersect_Meshes_2D::Facets::find_facet unable to locate facet in facet list\n";
#endif
            throw runtime_error("Unable to locate facet");
        }
        
#ifdef DEBUG_INTERSECT_MESHES_2D_FACETS
        cout << "Intersect_Meshes_2D::Facets::find_facet found facet p1: " << f.get_p1_index() << 
                " p2: " << f.get_p2_index() << " p3: " << f.get_p3_index() << "\n";
#endif
        return f;
    }
    
    Intersect_Meshes_2D::Facet_Builder::Line_Segment::Line_Segment(const shared_ptr<Point_2D>& p1, 
            const shared_ptr<Point_2D>& p2, const Location loc) : used(false), location(loc), point1(p1), point2(p2) {}
    
    const bool Intersect_Meshes_2D::Facet_Builder::Line_Segment::operator==(const Line_Segment& seg) const
    {
        return (point1 == seg.point1 || point1 == seg.point2) && 
                (point2 == seg.point2 || point2 == seg.point1);
    }
    
    const bool Intersect_Meshes_2D::Facet_Builder::Line_Segment::shares_pt(
            const Line_Segment& seg, shared_ptr<Point_2D>& shared_pt) const
    {
        if (point1 == seg.point1 || point1 == seg.point2)
        {
            shared_pt = point1;
            return true;
        }
        else if (point2 == seg.point1 || point2 == seg.point2)
        {
            shared_pt = point2;
            return true;
        }
        
        return false;
    }
    
    Intersect_Meshes_2D::Facet_Builder::Segments::Segment_Sort::Segment_Sort() {}

    const bool Intersect_Meshes_2D::Facet_Builder::Segments::Segment_Sort::operator ()(
            const Line_Segment& seg1, const Line_Segment& seg2) const
    {
        // move internal segments to the front
        return seg1.location != Location::internal && seg2.location == Location::internal;
    }
    
    Intersect_Meshes_2D::Facet_Builder::Segments::Segments() : segments(), removed_segs() {}
    
    const Intersect_Meshes_2D::Facet_Builder::Segments::size_type Intersect_Meshes_2D::Facet_Builder::Segments::size() const
    {
        return segments.size();
    }
    
    void Intersect_Meshes_2D::Facet_Builder::Segments::clear() 
    {
        segments.clear();
        removed_segs.clear();
    }

    const bool Intersect_Meshes_2D::Facet_Builder::Segments::contains_segment(
            const shared_ptr<Point_2D>& p1, const shared_ptr<Point_2D>& p2)
    {
        // location does not matter, so use p1 since it isn't a valid segment location
        Line_Segment seg(p1, p2, Location::p1);
        return segments.end() != find(segments.begin(), segments.end(), seg);
    }
    
    void Intersect_Meshes_2D::Facet_Builder::Segments::add_external_segment(
            const shared_ptr<Point_2D>& p1, const shared_ptr<Point_2D>& p2, const Location side_loc)
    {
        // don't check if it already exists because it shouldn't
        segments.push_back(Line_Segment(p1, p2, side_loc));
        // don't sort because is is already in the back
    }
    
    void Intersect_Meshes_2D::Facet_Builder::Segments::add_internal_segment(
            const shared_ptr<Point_2D>& p1, const shared_ptr<Point_2D>& p2)
    {
        Line_Segment seg(p1, p2, Location::internal);
        if (segments.end() == find(segments.begin(), segments.end(), seg))
        {
            segments.push_back(Line_Segment(p1, p2, Location::internal));
            sort(segments.begin(), segments.end(), Segment_Sort());
        }
    }
    
    const Intersect_Meshes_2D::Facet_Builder::Line_Segment* Intersect_Meshes_2D::Facet_Builder::Segments::get_next_segment(
            const Line_Segment* prev_segment) const
    {
        if (!segments.empty())
        {
            vector<Line_Segment>::const_iterator it = segments.begin();
            if (prev_segment != 0)
            {
                while (it != segments.end())
                {
                    if (*it == *prev_segment)
                    {
                        ++it; // advance iterator to next segment
                        break;
                    }
                    
                    ++it;
                }
            }
            if (it != segments.end())
                return &*it;
        }
        
        return 0;
    }
    
    const Intersect_Meshes_2D::Facet_Builder::Line_Segment* Intersect_Meshes_2D::Facet_Builder::Segments::find_segment(
            const shared_ptr<Point_2D>& p1, const shared_ptr<Point_2D>& p2) const
    {
        // location does not matter, so use p1 since it is not a valid segment location
        Line_Segment seg(p1, p2, Location::p1);
        vector<Line_Segment>::const_iterator it = find(segments.begin(), segments.end(), seg);
        if (it != segments.end())
            return &*it;
        return 0;
    }
    
    const Intersect_Meshes_2D::Facet_Builder::Line_Segment* Intersect_Meshes_2D::Facet_Builder::Segments::find_connecting_segment(
            const Line_Segment& segment, const Line_Segment* prev_connecting_seg, shared_ptr<Point_2D>& shared_pt, 
            const Point_2D::Measurement precision) const
    {
#ifdef DEBUG_INTERSECT_MESHES_2D_FACET_BUILDER_SEGMENTS
        cout << "Intersect_Meshes_2D::Facet_Builder::Segments::find_connecting_segment begin\n";
#endif
        bool found(false);
        shared_ptr<Point_2D> shared_point;
        if (!segments.empty())
        {
#ifdef DEBUG_INTERSECT_MESHES_2D_FACET_BUILDER_SEGMENTS
            for (vector<Line_Segment>::const_iterator t_it = segments.begin(); t_it != segments.end(); ++t_it)
                cout << "Intersect_Meshes_2D::Facet_Builder::Segments::find_connecting_segment segments p1 x: " << 
                        (*t_it).point1->get_x() << " y: " << (*t_it).point1->get_y() << " p2 x: " << 
                        (*t_it).point2->get_x() << " y: " << (*t_it).point2->get_y() << "\n";
#endif
            vector<Line_Segment>::const_iterator it = segments.begin();
            if (prev_connecting_seg != 0)
            {
#ifdef DEBUG_INTERSECT_MESHES_2D_FACET_BUILDER_SEGMENTS
                cout << "Intersect_Meshes_2D::Facet_Builder::Segments::find_conneting_segment locating last connecting segment in segments p1 x: " << 
                        prev_connecting_seg->point1->get_x() << " y: " << prev_connecting_seg->point1->get_y() << 
                        " p2 x: " << prev_connecting_seg->point2->get_x() << " y: " << prev_connecting_seg->point2->get_y() << "\n";
#endif
                while (it != segments.end())
                {
                    if (*it == *prev_connecting_seg)
                    {
#ifdef DEBUG_INTERSECT_MESHES_2D_FACET_BUILDER_SEGMENTS
                        cout << "Intersect_Meshes_2D::Facet_Builder::Segments::find_connecting_segment found segment\n";
#endif
                        found = true;
                        ++it; // move to next segment or end
                        break;
                    }
                    
                    ++it;
                }
            }
            while (it != segments.end())
            {
#ifdef DEBUG_INTERSECT_MESHES_2D_FACET_BUILDER_SEGMENTS
                cout << "Intersect_Meshes_2D::Facet_Builder::Segments::find_connecting_segment checking segment p1 x: " << 
                        it->point1->get_x() << " y: " << it->point1->get_y() << " p2 x: " << it->point2->get_x() << " y: " << 
                        it->point2->get_y() << "\n";
#endif
                if (segment == *it) // do not return the same segment
                {
                    ++it;
                    continue;
                }
                
                if (segment.shares_pt(*it, shared_point))
                {
#ifdef DEBUG_INTERSECT_MESHES_2D_FACET_BUILDER_SEGMENTS
                    cout << "Intersect_Meshes_2D::Facet_Builder::Segments::find_connecting_segment segment shared pt x: " << 
                            shared_point->get_x() << " y: " << shared_point->get_y() << "\n";
#endif
                    if (segment.location == Location::internal)
                    {
                        if (it->location == Location::internal)
                        {
                            // both segments are internal, so check if they form a straight line
                            bool same_direction(false);
                            shared_ptr<Point_2D> seg1_pt(segment.point1 == shared_point ? segment.point2 : segment.point1);
                            shared_ptr<Point_2D> seg2_pt((*it).point1 == shared_point ? (*it).point2 : (*it).point1);
                            if (is_same_line(*shared_point, *seg1_pt, *shared_point, *seg2_pt, same_direction, precision) && !same_direction)
                            {
                                ++it; // same line so try next segment
                                continue;
                            }
#ifdef DEBUG_INTERSECT_MESHES_2D_FACET_BUILDER_SEGMENTS
                            cout << "Intersect_Meshes_2D::Facet_Builder::Segments::find_connecting_segment found connecting segment\n";
#endif
                            shared_pt = shared_point;
                            return &*it;
                        }
                        
                        // connecting segment is a perimeter segment
#ifdef DEBUG_INTERSECT_MESHES_2D_FACET_BUILDER_SEGMENTS
                        cout << "Intersect_Meshes_2D::Facet_Builder::Segments::find_connecting_segment found connecting segment\n";
#endif
                        // no need to check if it is on the same side
                        shared_pt = shared_point;
                        return &*it;
                    }
                    
                    // segment is a perimeter segment
                    if (it->location == segment.location) // connecting segment is on the same side
                    {
                        ++it; // do not return a segment from the same side
                        continue;
                    }
#ifdef DEBUG_INTERSECT_MESHES_2D_FACET_BUILDER_SEGMENTS
                    cout << "Intersect_Meshes_2D::Facet_Builder::Segments::find_connecting_segment found connecting segment\n";
#endif
                    // segment is external, so do not need to check if it is in a straight line
                    shared_pt = shared_point;
                    return &*it;
                }
                
                ++it;
            }
        }
        
#ifdef DEBUG_INTERSECT_MESHES_2D_FACET_BUILDER_SEGMENTS
        cout << "Intersect_Meshes_2D::Facet_Builder::Segments::find_connecting_segment end returning 0\n";
#endif
        return 0;
    }
    
    const bool Intersect_Meshes_2D::Facet_Builder::Segments::is_segment_external(
            const Line_Segment& seg3, Location& side_loc, const Facet_2D& original_facet,
            const vector<shared_ptr<Point_2D>>& p1p2_pts, const vector<shared_ptr<Point_2D>>& p1p3_pts, 
            const vector<shared_ptr<Point_2D>>& p2p3_pts, const vector<shared_ptr<Point_2D>>& internal_pts) const
    {
        for (vector<shared_ptr<Point_2D>>::const_iterator it = internal_pts.begin(); it != internal_pts.end(); ++it)
        {
            if (*it == seg3.point1 || *it == seg3.point2)
                return false; // if either end of seg3 is an internal point, then it is an internal segment
        }
        
        bool found_p1(false);
        bool found_p2(false);
        
        // check p1p2 side
        if (seg3.point1 == original_facet.get_point1() || seg3.point1 == original_facet.get_point2())
            found_p1 = true;
        if (seg3.point2 == original_facet.get_point1() || seg3.point2 == original_facet.get_point2())
            found_p2 = true;
        for (vector<shared_ptr<Point_2D>>::const_iterator it = p1p2_pts.begin(); it != p1p2_pts.end(); ++it)
        {
            if (*it == seg3.point1)
                found_p1 = true;
            else if (*it == seg3.point2)
                found_p2 = true;
        }
        
        if (found_p1 && found_p2)
        {
            side_loc = Location::p1p2;
            return true;
        }
        
        found_p1 = false;
        found_p2 = false;
        
        // check p1p3 side
        if (seg3.point1 == original_facet.get_point1() || seg3.point1 == original_facet.get_point3())
            found_p1 = true;
        if (seg3.point2 == original_facet.get_point1() || seg3.point2 == original_facet.get_point3())
            found_p2 = true;
        for (vector<shared_ptr<Point_2D>>::const_iterator it = p1p3_pts.begin(); it != p1p3_pts.end(); ++it)
        {
            if (*it == seg3.point1)
                found_p1 = true;
            else if (*it == seg3.point2)
                found_p2 = true;
        }
        
        if (found_p1 && found_p2)
        {
            side_loc = Location::p1p3;
            return true;
        }
        
        found_p1 = false;
        found_p2 = false;
        
        // check p2p3 side
        if (seg3.point1 == original_facet.get_point2() || seg3.point1 == original_facet.get_point3())
            found_p1 = true;
        if (seg3.point2 == original_facet.get_point2() || seg3.point2 == original_facet.get_point3())
            found_p2 = true;
        for (vector<shared_ptr<Point_2D>>::const_iterator it = p2p3_pts.begin(); it != p2p3_pts.end(); ++it)
        {
            if (*it == seg3.point1)
                found_p1 = true;
            else if (*it == seg3.point2)
                found_p2 = true;
        }
        
        if (found_p1 && found_p2)
        {
            side_loc = Location::p2p3;
            return true;
        }
            
        return false;
    }

    const bool Intersect_Meshes_2D::Facet_Builder::Segments::does_seg_intersect(
            const Line_Segment& segment, const Facet_2D& orig_facet, 
            const vector<shared_ptr<Point_2D>>& p1p2_pts, const vector<shared_ptr<Point_2D>>& p1p3_pts, 
            const vector<shared_ptr<Point_2D>>& p2p3_pts, const vector<shared_ptr<Point_2D>>& internal_points, 
            const Point_2D::Measurement precision) const
    {
#ifdef DEBUG_INTERSECT_MESHES_2D_FACET_BUILDER_SEGMENTS
        cout << "Intersect_Meshes_2D::Facet_Builder::Segments::does_seg_intersect checking segment p1 x: " << 
                segment.point1->get_x() << " y: " << segment.point1->get_y() << " ptr: " << segment.point1.get() <<  
                " p2 x: " << segment.point2->get_x() << " y: " << segment.point2->get_y() << 
                " ptr: " << segment.point2.get() << " location: " << segment.location << "\n";
#endif
        if (removed_segs.end() != find(removed_segs.begin(), removed_segs.end(), segment))
        {
#ifdef DEBUG_INTERSECT_MESHES_2D_FACET_BUILDER_SEGMENTS
            cout << "Intersect_Meshes_2D::Facet_Builder::Segments::does_seg_intersect segment matches removed segment. returning true\n";
#endif
            return true; // segment has already been processed, return true
        }
        if (segments.end() != find(segments.begin(), segments.end(), segment))
        {
#ifdef DEBUG_INTERSECT_MESHES_2D_FACET_BUILDER_SEGMENTS
            cout << "Intersect_Meshes_2D::Facet_Builder::Segments::does_seg_intersect segment matches segment. returning true\n";
#endif
            return false; // segment is still being processed, return false
        }
        
        Location side_loc(Location::internal);
        const bool segment_is_ext(is_segment_external(segment, side_loc, orig_facet, p1p2_pts, p1p3_pts, p2p3_pts, internal_points));
#ifdef DEBUG_INTERSECT_MESHES_2D_FACET_BUILDER_SEGMENTS
        cout << "Intersect_Meshes_2D::Facet_Builder::Segments::does_seg_intersect segment_is_ext: " << (segment_is_ext ? "TRUE" : "FALSE") << " side_loc: " << side_loc << "\n";
#endif
        
        // check if segment crosses any removed segments
        for (vector<Line_Segment>::const_iterator it = removed_segs.begin(); it != removed_segs.end(); ++it)
        {
#ifdef DEBUG_INTERSECT_MESHES_2D_FACET_BUILDER_SEGMENTS
            cout << "Intersect_Meshes_2D::Facet_Builder::Segments::does_seg_intersect testing removed segment p1 x: " << 
                    it->point1->get_x() << " y: " << it->point1->get_y() << " ptr: " << it->point1.get() <<  " p2 x: " << 
                    it->point2->get_x() << " y: " << it->point2->get_y() << " ptr: " << it->point2.get() << " location: " << 
                    it->location << "\n";
#endif
            if (segment_is_ext && it->location == side_loc) // if both are external and on the same side
            {
                shared_ptr<Point_2D> shared_pt;
                if (segment.shares_pt(*it, shared_pt))
                {
#ifdef DEBUG_INTERSECT_MESHES_2D_FACET_BUILDER_SEGMENTS
                    cout << "Intersect_Meshes_2D::Facet_Builder::Segments::does_seg_intersect removed external segment shares point with segment\n";
#endif

                    if ((*it).location == side_loc)
                    {
#ifdef DEBUG_INTERSECT_MESHES_2D_FACET_BUILDER_SEGMENTS
                        cout << "Intersect_Meshes_2D::Facet_Builder::Segments::does_seg_intersect testing segment end points\n";
#endif
                        shared_ptr<Point_2D> non_common_pt((*it).point1 == shared_pt ? (*it).point2 : (*it).point1);
                        if (is_pt_on_vector(*non_common_pt, *segment.point1, *segment.point2, precision))
                        {
#ifdef DEBUG_INTERSECT_MESHES_2D_FACET_BUILDER_SEGMENTS
                            cout << "Intersect_Meshes_2D::Facet_Builder::Segments::does_seg_intersect segment intersects removed external segment\n";
#endif
                            return true;
                        }
                        non_common_pt = segment.point1 == shared_pt ? segment.point2 : segment.point1;
                        if (is_pt_on_vector(*non_common_pt, *(*it).point1, *(*it).point2, precision))
                        {
#ifdef DEBUG_INTERSECT_MESHES_2D_FACET_BUILDER_SEGMENTS
                            cout << "Intersect_Meshes_2D::Facet_Builder::Segments::does_seg_intersect segment intersects removed external segment\n";
#endif
                            return true;
                        }
                    }
                }
            }
            else if (!segment_is_ext && it->location == Location::internal) // both segments are internal
            {
                shared_ptr<Point_2D> shared_pt;
                if (!segment.shares_pt(*it, shared_pt))
                {
#ifdef DEBUG_INTERSECT_MESHES_2D_FACET_BUILDER_SEGMENTS
                    cout << "Intersect_Meshes_2D::Facet_Builder::Segments::does_seg_intersect intersecting vectors\n";
#endif
                    // intersect
                    Vector_2D_idata idata;
                    if (intersect_vectors(*segment.point1, *segment.point2, *(*it).point1, *(*it).point2, idata, precision))
                    {
#ifdef DEBUG_INTERSECT_MESHES_2D_FACET_BUILDER_SEGMENTS
                        cout << "Intersect_Meshes_2D::Facet_Builder::Segments::does_seg_intersect segment intersects a removed segment\n";
#endif
                        return true; // vectors do not share a common point, but did intersect
                    }
                }
                else // segments share point
                {
                    Vector_2D_idata idata;
                    if (intersect_vectors(*segment.point1, *segment.point2, *(it->point1), *(it->point2), idata, precision) && idata.num == 2)
                    {
#ifdef DEBUG_INTERSECT_MESHES_2D_FACET_BUILDER_SEGMENTS
                        cout << "Intersect_Meshes_2D::Facet_Builder::Segments::does_seg_intersect segment intersects a removed segment\n";
#endif
                        return true; // vectors shared a common point and an additional intersection
                    }
                }
            }
        }
        
        // check if segment crosses any segments
        for (vector<Line_Segment>::const_iterator it = segments.begin(); it != segments.end(); ++it)
        {
#ifdef DEBUG_INTERSECT_MESHES_2D_FACET_BUILDER_SEGMENTS
            cout << "Intersect_Meshes_2D::Facet_Builder::Segments::does_seg_intersect testing segment p1 x: " << 
                    it->point1->get_x() << " y: " << it->point1->get_y() << " ptr: " << it->point1.get() <<  
                    " p2 x: " << it->point2->get_x() << " y: " << it->point2->get_y() << " ptr: " << it->point2.get() << "\n";
#endif
            if (segment_is_ext && side_loc == it->location) // both segments are on the same side
            {
                shared_ptr<Point_2D> shared_pt;
                if (segment.shares_pt(*it, shared_pt))
                {
#ifdef DEBUG_INTERSECT_MESHES_2D_FACET_BUILDER_SEGMENTS
                    cout << "Intersect_Meshes_2D::Facet_Builder::Segments::does_seg_intersect segment shares point with existing segment\n";
                    cout << "Intersect_Meshes_2D::Facet_Builder::Segments::does_seg_intersect testing segment end points\n";
#endif
                    shared_ptr<Point_2D> non_common_pt(it->point1 == shared_pt ? it->point2 : it->point1);
                    if (is_pt_on_vector(*non_common_pt, *segment.point1, *segment.point2, precision))
                    {
#ifdef DEBUG_INTERSECT_MESHES_2D_FACET_BUILDER_SEGMENTS
                        cout << "Intersect_Meshes_2D::Facet_Builder::Segments::does_seg_intersect segment intersects external segment\n";
#endif
                        return true;
                    }
                    if (is_pt_on_vector(*non_common_pt, *(*it).point1, *(*it).point2, precision))
                    {
#ifdef DEBUG_INTERSECT_MESHES_2D_FACET_BUILDER_SEGMENTS
                        cout << "Intersect_Meshes_2D::Facet_Builder::Segments::does_seg_intersect segment intersects external segment\n";
#endif
                        return true;
                    }
                } // don't test if they don't share a common point - shouldn't form an intersect pt because *it is external
            }
            else if (!segment_is_ext && it->location == Location::internal) // both are internal
            {
                shared_ptr<Point_2D> shared_pt;
                if (!segment.shares_pt(*it, shared_pt))
                {
#ifdef DEBUG_INTERSECT_MESHES_2D_FACET_BUILDER_SEGMENTS
                    cout << "Intersect_Meshes_2D::Facet_Builder::Segments::does_seg_intersect intersecting vectors\n";
#endif
                    // intersect
                    Vector_2D_idata idata;
                    if (intersect_vectors(*segment.point1, *segment.point2, *(*it).point1, *(*it).point2, idata, precision))
                    {
#ifdef DEBUG_INTERSECT_MESHES_2D_FACET_BUILDER_SEGMENTS
                        cout << "Intersect_Meshes_2D::Facet_Builder::Segments::does_seg_intersect segment intersects an internal segment\n";
#endif
                        return true; // vectors do not share a common point, but did intersect
                    }
                }
                else // segments share point
                {
                    Vector_2D_idata idata;
                    if (intersect_vectors(*segment.point1, *segment.point2, *(it->point1), *(it->point2), idata, precision) && idata.num == 2)
                    {
#ifdef DEBUG_INTERSECT_MESHES_2D_FACET_BUILDER_SEGMENTS
                        cout << "Intersect_Meshes_2D::Facet_Builder::Segments::does_seg_intersect segment intersects an internal segment\n";
#endif
                        return true; // vectors shared a common point and an additional intersection
                    }
                }
            }
        }

        return false;
    }
    
    void Intersect_Meshes_2D::Facet_Builder::Segments::process_used_segment(vector<Line_Segment>::iterator seg)
    {
        if (seg->location == Location::internal)
        {
#ifdef DEBUG_INTERSECT_MESHES_2D_FACET_BUILDER_SEGMENTS
            cout << "Intersect_Meshes_2D::Facet_Builder::Segments::process_used_segments segment is an internal segment\n";
#endif
            if (seg->used)
            {
#ifdef DEBUG_INTERSECT_MESHES_2D_FACET_BUILDER_SEGMENTS
                cout << "Intersect_Meshes_2D::Facet_Builder::Segments::process_used_segments erasing segment\n";
#endif
                removed_segs.push_back(*seg);
                segments.erase(seg);
            }
            else // update used to true
            {
#ifdef DEBUG_INTERSECT_MESHES_2D_FACET_BUILDER_SEGMENTS
                cout << "Intersect_Meshes_2D::Facet_Builder::Segments::process_used_segments setting segment used to true\n";
#endif
                seg->used = true;
            }
        }
        else // external segment
        {
#ifdef DEBUG_INTERSECT_MESHES_2D_FACET_BUILDER_SEGMENTS
            cout << "Intersect_Meshes_2D::Facet_Builder::Segments::process_used_segments segment is an external segment. erasing\n";
#endif
            removed_segs.push_back(*seg);
            segments.erase(seg);
        }
    }
    
    // update used or remove segments from finished facets
    void Intersect_Meshes_2D::Facet_Builder::Segments::process_used_segments(
            const Line_Segment& seg1, const Line_Segment& seg2)
    {
#ifdef DEBUG_INTERSECT_MESHES_2D_FACET_BUILDER_SEGMENTS
        cout << "Intersect_Meshes_2D::Facet_Builder::Segments::process_used_segments begin\n";
        cout << "Intersect_Meshes_2D::Facet_Builder::Segments::process_used_segments seg1 p1 x: " << 
                seg1.point1->get_x() << " y: " << seg1.point1->get_y() << " p2 x: " << seg1.point2->get_x() << 
                " y: " << seg1.point2->get_y() << " location: " << seg1.location << " used: " << seg1.used << "\n";
        cout << "Intersect_Meshes_2D::Facet_Builder::Segments::process_used_segments seg2 p1 x: " << 
                seg2.point1->get_x() << " y: " << seg2.point1->get_y() << " p2 x: " << seg2.point2->get_x() << 
                " y: " << seg2.point2->get_y() << " location: " << seg2.location << " used: " << seg2.used << "\n";
#endif
        // determine segment 3
        shared_ptr<Point_2D> shared_pt;
        if (!seg1.shares_pt(seg2, shared_pt))
            throw runtime_error("segment1 and segment2 do not share a common point");
        // location should not matter: use location p1 since it is not a valid segment location
        Line_Segment seg3(seg1.point1 == shared_pt ? seg1.point2 : seg1.point1, seg2.point1 == shared_pt ? seg2.point2 : seg2.point1, Location::p1);
#ifdef DEBUG_INTERSECT_MESHES_2D_FACET_BUILDER_SEGMENTS
        cout << "Intersect_Meshes_2D::Facet_Builder::Segments::process_used_segments seg3 p1 x: " << 
                seg3.point1->get_x() << " y: " << seg3.point1->get_y() << " p2 x: " << 
                seg3.point2->get_x() << " y: " << seg3.point2->get_y() << "\n";
#endif

        vector<Line_Segment>::iterator it = find(segments.begin(), segments.end(), seg1);
        if (it == segments.end())
            throw runtime_error("Unable to locate segment1");
        process_used_segment(it);
        
        it = find(segments.begin(), segments.end(), seg2);
        if (it == segments.end())
            throw runtime_error("Unable to locate segment2");
        process_used_segment(it);
        
        it = find(segments.begin(), segments.end(), seg3);
        if (it == segments.end())
        {
            // segment3 does not exist
#ifdef DEBUG_INTERSECT_MESHES_2D_FACET_BUILDER_SEGMENTS
            cout << "Intersect_Meshes_2D::Facet_Builder::Segments::process_used_segments seg3 does not exist. adding\n";
#endif
            // add line segment 3 because it was formed in the build algorithm
            this->add_internal_segment(seg3.point1, seg3.point2);
            it = find(segments.begin(), segments.end(), seg3);
        }
        process_used_segment(it);
        
#ifdef DEBUG_INTERSECT_MESHES_2D_FACET_BUILDER_SEGMENTS
        cout << "Intersect_Meshes_2D::Facet_Builder::Segments::process_used_segments end\n";
#endif
    }
    
    Intersect_Meshes_2D::Facet_Builder::Intersecting_Facet_Side_Pts::Intersecting_Facet_Side_Pts() : p1p2_points(), p1p3_points(), p2p3_points() {}
    
    Intersect_Meshes_2D::Facet_Builder::Segment_Find::Segment_Find(const shared_ptr<Point_2D>& point) : pt(point) {}
    
    const bool Intersect_Meshes_2D::Facet_Builder::Segment_Find::operator ()(const Intersect_Meshes_2D::Facet_Builder::Line_Segment& segment)
    {
        return pt == segment.point1 || pt == segment.point2;
    }
    
    Intersect_Meshes_2D::Facet_Builder::Facet_Builder(const bool for_f1, 
            const Facet_2D& original_facet, const Point_2D::Measurement prec) : 
            for_facet1(for_f1), precision(prec), orig_facet(original_facet), 
            internal_pts(), p1p2_pts(), p1p3_pts(), p2p3_pts(), segments() {}

    void Intersect_Meshes_2D::Facet_Builder::process_two_i_pts(const I_Pt_List& i_pts)
    {
        I_Pt_List::const_iterator it = i_pts.begin();
        Intersect_Point ip1(*it);
        ++it;
        Intersect_Point ip2(*it);
        
        Location ip1_loc = for_facet1 ? ip1.f1_loc : ip1.f2_loc;
        Location ip2_loc = for_facet1 ? ip2.f1_loc : ip2.f2_loc;
        
#ifdef DEBUG_INTERSECT_MESHES_2D_FACET_BUILDER
        cout << "Intersect_Meshes_2D::Facet_Builder::process_two_i_pts ip1 x: " << 
                ip1.pt->get_x() << " y: " << ip1.pt->get_y() << " ptr: " << (ip1.pt) << 
                " loc " << ip1_loc << "\n";
        cout << "Intersect_Meshes_2D::Facet_Builder::process_two_i_pts ip2 x: " << 
                ip2.pt->get_x() << " y: " << ip2.pt->get_y() << " ptr: " << (ip2.pt) << 
                " loc " << ip2_loc << "\n";
#endif
        
        switch (ip1_loc) 
        {
            case (Location::internal):
#ifdef DEBUG_INTERSECT_MESHES_2D_FACET_BUILDER
                cout << "Intersect_Meshes_2D::Facet_Builder::process_two_i_pts adding segment p1 x: " << 
                        ip1.pt->get_x() << " y: " << ip1.pt->get_y() << " p2 x: " << ip2.pt->get_x() << 
                        " y: " << ip2.pt->get_y() << "\n";
#endif
                segments.add_internal_segment(ip1.pt, ip2.pt);
                break;
            case (Location::p1p2):
                if (ip2_loc != Location::p1p2 && ip2_loc != Location::p1 && 
                        ip2_loc != Location::p2)
                {
#ifdef DEBUG_INTERSECT_MESHES_2D_FACET_BUILDER
                    cout << "Intersect_Meshes_2D::Facet_Builder::process_two_i_pts adding segment p1 x: " << 
                            ip1.pt->get_x() << " y: " << ip1.pt->get_y() << " p2 x: " << ip2.pt->get_x() << 
                            " y: " << ip2.pt->get_y() << "\n";
#endif
                    segments.add_internal_segment(ip1.pt, ip2.pt);
                }
                break;
            case (Location::p1p3):
                if (ip2_loc != Location::p1p3 && ip2_loc != Location::p1 && 
                        ip2_loc != Location::p3)
                {
#ifdef DEBUG_INTERSECT_MESHES_2D_FACET_BUILDER
                    cout << "Intersect_Meshes_2D::Facet_Builder::process_two_i_pts adding segment p1 x: " << 
                            ip1.pt->get_x() << " y: " << ip1.pt->get_y() << " p2 x: " << ip2.pt->get_x() << 
                            " y: " << ip2.pt->get_y() << "\n";
#endif
                    segments.add_internal_segment(ip1.pt, ip2.pt);
                }
                break;
            case (Location::p2p3):
                if (ip2_loc != Location::p2p3 && ip2_loc != Location::p2 && 
                        ip2_loc != Location::p3)
                {
#ifdef DEBUG_INTERSECT_MESHES_2D_FACET_BUILDER
                    cout << "Intersect_Meshes_2D::Facet_Builder::process_two_i_pts adding segment p1 x: " << 
                            ip1.pt->get_x() << " y: " << ip1.pt->get_y() << " p2 x: " << ip2.pt->get_x() << 
                            " y: " << ip2.pt->get_y() << "\n";
#endif
                    segments.add_internal_segment(ip1.pt, ip2.pt);
                }
                break;
            case (Location::p1):
                if (ip2_loc != Location::p1p2 && ip2_loc != Location::p1p3 && 
                        ip2_loc != Location::p2 && ip2_loc != Location::p3)
                {
#ifdef DEBUG_INTERSECT_MESHES_2D_FACET_BUILDER
                    cout << "Intersect_Meshes_2D::Facet_Builder::process_two_i_pts adding segment p1 x: " << 
                            ip1.pt->get_x() << " y: " << ip1.pt->get_y() << " p2 x: " << ip2.pt->get_x() << 
                            " y: " << ip2.pt->get_y() << "\n";
#endif
                    segments.add_internal_segment(ip1.pt, ip2.pt);
                }
                break;
            case (Location::p2):
                if (ip2_loc != Location::p1p2 && ip2_loc != Location::p2p3 && 
                        ip2_loc != Location::p1 && ip2_loc != Location::p3)
                {
#ifdef DEBUG_INTERSECT_MESHES_2D_FACET_BUILDER
                    cout << "Intersect_Meshes_2D::Facet_Builder::process_two_i_pts adding segment p1 x: " << 
                            ip1.pt->get_x() << " y: " << ip1.pt->get_y() << " p2 x: " << ip2.pt->get_x() << 
                            " y: " << ip2.pt->get_y() << "\n";
#endif
                    segments.add_internal_segment(ip1.pt, ip2.pt);
                }
                break;
            case (Location::p3):
                if (ip2_loc != Location::p1p3 && ip2_loc != Location::p2p3 && 
                        ip2_loc != Location::p1 && ip2_loc != Location::p2)
                {
#ifdef DEBUG_INTERSECT_MESHES_2D_FACET_BUILDER
                    cout << "Intersect_Meshes_2D::Facet_Builder::process_two_i_pts adding segment p1 x: " << 
                            ip1.pt->get_x() << " y: " << ip1.pt->get_y() << " p2 x: " << ip2.pt->get_x() << 
                            " y: " << ip2.pt->get_y() << "\n";
#endif
                    segments.add_internal_segment(ip1.pt, ip2.pt);
                }
                break;
            default:
                throw runtime_error("undefined point location");
        }
    }
    
    void Intersect_Meshes_2D::Facet_Builder::process_i_side_pts(const I_Pt_List& i_side_pts, 
            I_Pt_List& t_intersect_pts)
    {
        if (i_side_pts.size() == 2)
        {
            process_two_i_pts(i_side_pts);
            // remove from temp intersect pts
            I_Pt_List::const_iterator t_it = find(t_intersect_pts.begin(), t_intersect_pts.end(), *i_side_pts.begin());
            if (t_it != t_intersect_pts.end())
                t_intersect_pts.erase(t_it);
            t_it = find(t_intersect_pts.begin(), t_intersect_pts.end(), *(++i_side_pts.begin()));
            if (t_it != t_intersect_pts.end())
                t_intersect_pts.erase(t_it);
        }
        else if (i_side_pts.size() > 2)
        {
#ifdef DEBUG_INTERSECT_MESHES_2D_FACET_BUILDER
            int count = 0;
            for (I_Pt_List::const_iterator it = i_side_pts.begin(); it != i_side_pts.end(); ++it)
            {
                cout << "Intersect_Meshes_2D::Facet_Builder::process_i_side_pts: i_side_pts[" << count++ << "] pt x: " << 
                        (*it).pt->get_x() << " y: " << (*it).pt->get_y() << " f1_loc: " << (*it).f1_loc << " f2_loc: " << 
                        (*it).f2_loc << "\n";
            }
#endif
            throw runtime_error("invalid number of facet intersecting side points");
        }
    }
    
    void Intersect_Meshes_2D::Facet_Builder::gen_internal_segs(const I_Pt_List& intersect_points, 
            const Intersecting_Facet_Side_Pts intersecting_side_pts)
    {
#ifdef DEBUG_INTERSECT_MESHES_2D_FACET_BUILDER
        cout << "Intersect_Meshes_2D::Facet_Builder::gen_internal_segs begin for_facet1: " << for_facet1 << "\n";
#endif
        I_Pt_List t_intersect_pts;
        for (I_Pt_List::const_iterator it = intersect_points.begin(); it != intersect_points.end(); ++it)
            t_intersect_pts.push_back(*it);
        
#ifdef DEBUG_INTERSECT_MESHES_2D_FACET_BUILDER
        cout << "Intersect_Meshes_2D::Facet_Builder::gen_internal_segs intersecting_side_pts.p1p2_points.size() is " << 
                intersecting_side_pts.p1p2_points.size() << "\n";
#endif
        process_i_side_pts(intersecting_side_pts.p1p2_points, t_intersect_pts);
#ifdef DEBUG_INTERSECT_MESHES_2D_FACET_BUILDER
        cout << "Intersect_Meshes_2D::Facet_Builder::gen_internal_segs intersecting_side_pts.p1p3_points.size() is " << 
                intersecting_side_pts.p1p3_points.size() << "\n";
#endif
        process_i_side_pts(intersecting_side_pts.p1p3_points, t_intersect_pts);
#ifdef DEBUG_INTERSECT_MESHES_2D_FACET_BUILDER
        cout << "Intersect_Meshes_2D::Facet_Builder::gen_internal_segs intersecting_side_pts.p2p3_points.size() is " << 
                intersecting_side_pts.p2p3_points.size() << "\n";
#endif
        process_i_side_pts(intersecting_side_pts.p2p3_points, t_intersect_pts);

        // process remaining points
#ifdef DEBUG_INTERSECT_MESHES_2D_FACET_BUILDER
        cout << "Intersect_Meshes_2D::Facet_Builder::gen_internal_segs facet leftover intersect points: " << 
                t_intersect_pts.size() << "\n";
        int count = 0;
        for (I_Pt_List::const_iterator it = t_intersect_pts.begin(); it != t_intersect_pts.end(); ++it)
        {
            cout << "Intersect_Meshes_2D::Facet_Builder::gen_internal_segs t_intersect_pts[" << 
                    count++ << "] x: " << it->pt->get_x() << " y: " << it->pt->get_y() << 
                    " f1_loc: " << it->f1_loc << " f2_loc: " << it->f2_loc << "\n";
        }
#endif
        
        // process remaining single_i_points to form any further internal segments
        if (t_intersect_pts.size() == 1)
        {
            const Intersect_Point* ip1(&*t_intersect_pts.begin());
#ifdef DEBUG_INTERSECT_MESHES_2D_FACET_BUILDER
            cout << "Intersect_Meshes_2D::Facet_Builder::gen_internal_segs t_intersect_pts ip1 x: " << 
                    ip1->pt->get_x() << " y: " << ip1->pt->get_y() << " f1_loc: " << ip1->f1_loc << 
                    " f2_loc: " << ip1->f2_loc << "\n";
#endif
        }
        else if (t_intersect_pts.size() == 2)
        {
            // process points
            I_Pt_List leftover_pts;
            I_Pt_List::const_iterator t_it = t_intersect_pts.begin();
            leftover_pts.push_back(*t_it);
            ++t_it;
            leftover_pts.push_back(*t_it);
            process_two_i_pts(leftover_pts);
        }
        else if (t_intersect_pts.size() == 3)
        {
            I_Pt_List::const_iterator it = t_intersect_pts.begin();
            const Intersect_Point* ip1(&*it);
            ++it;
            const Intersect_Point* ip2(&*it);
            ++it;
            const Intersect_Point* ip3(&*it);
            
#ifdef DEBUG_INTERSECT_MESHES_2D_FACET_BUILDER
            cout << "Intersect_Meshes_2D::Facet_Builder::gen_internal_segs 3 intersect points remaining\n";
            cout << "Intersect_Meshes_2D::Facet_Builder::gen_internal_segs ip1 x: " << ip1->pt->get_x() << 
                    " y: " << ip1->pt->get_y() << " f1_loc: " << ip1->f1_loc << " f2_loc: " << ip1->f2_loc << "\n";
            cout << "Intersect_Meshes_2D::Facet_Builder::gen_internal_segs ip2 x: " << ip2->pt->get_x() << 
                    " y: " << ip2->pt->get_y() << " f1_loc: " << ip2->f1_loc << " f2_loc: " << ip2->f2_loc << "\n";
            cout << "Intersect_Meshes_2D::Facet_Builder::gen_internal_segs ip3 x: " << ip3->pt->get_x() << 
                    " y: " << ip3->pt->get_y() << " f1_loc: " << ip3->f1_loc << " f2_loc: " << ip3->f2_loc << "\n";
#endif
            
            // if all points are internal
            if (for_facet1)
            {
                if (ip1->f1_loc == Location::internal && ip2->f1_loc == Location::internal && ip3->f1_loc == Location::internal)
                {
#ifdef DEBUG_INTERSECT_MESHES_2D_FACET_BUILDER
                    cout << "Intersect_Meshes_2D::Facet_Builder::gen_internal_segs all three intersect points are internal points - adding segments\n";
#endif
                    segments.add_internal_segment(ip1->pt, ip2->pt);
                    segments.add_internal_segment(ip2->pt, ip3->pt);
                    segments.add_internal_segment(ip3->pt, ip1->pt);
                }
                else if (!((ip1->f1_loc == Location::p1 || ip1->f1_loc == Location::p2 || ip1->f1_loc == Location::p3) &&
                        (ip2->f1_loc == Location::p1 || ip2->f1_loc == Location::p2 || ip2->f1_loc == Location::p3) &&
                        (ip3->f1_loc == Location::p1 || ip3->f1_loc == Location::p2 || ip3->f1_loc == Location::p3)))
                    throw runtime_error("invalid number of single intersect points: 3");
            }
            else
            {
                if (ip1->f2_loc == Location::internal && ip2->f2_loc == Location::internal && ip3->f2_loc == Location::internal)
                {
#ifdef DEBUG_INTERSECT_MESHES_2D_FACET_BUILDER
                    cout << "Intersect_Meshes_2D::Facet_Builder::gen_internal_segs all three points are internal points - adding segments\n";
#endif
                    segments.add_internal_segment(ip1->pt, ip2->pt);
                    segments.add_internal_segment(ip2->pt, ip3->pt);
                    segments.add_internal_segment(ip3->pt, ip1->pt);
                }
                else if (!((ip1->f2_loc == Location::p1 || ip1->f2_loc == Location::p2 || ip1->f2_loc == Location::p3) &&
                        (ip2->f2_loc == Location::p1 || ip2->f2_loc == Location::p2 || ip2->f2_loc == Location::p3) &&
                        (ip3->f2_loc == Location::p1 || ip3->f2_loc == Location::p2 || ip3->f2_loc == Location::p3)))
                    throw runtime_error("invalid number of single intersect points: 3");
            }
        }
        else if (t_intersect_pts.size() != 0)
            throw runtime_error("invalid number of facet intersect_points");
        
#ifdef DEBUG_INTERSECT_MESHES_2D_FACET_BUILDER
        cout << "Intersect_Meshes_2D::Facet_Builder::gen_internal_segs end. segments size: " << segments.size() << "\n";
#endif
    } 

    void Intersect_Meshes_2D::Facet_Builder::add_intersection(I_Pt_List& intersect_pts)
    {
        class Point_Find {
        public:
            Point_Find(const shared_ptr<Point_2D>& point) : pt(point) {}
            const bool operator()(const shared_ptr<Point_2D>& list_pt)
            {
                return pt->get_x() == list_pt->get_x() && pt->get_y() == list_pt->get_y();
            }
        private:
            const shared_ptr<Point_2D> pt;
        };
     
#ifdef DEBUG_INTERSECT_MESHES_2D_FACET_BUILDER
        cout << "Intersect_Meshes_2D::Facet_Builder::add_intersection Adding a new intersection\n";
        cout.flush();
        cout << "Intersect_Meshes_2D::Facet_Builder::add_intersection Current internal_points:\n";
        for (vector<shared_ptr<Point_2D>>::const_iterator it = internal_pts.begin(); it != internal_pts.end(); ++it)
        {
            cout << "Intersect_Meshes_2D::Facet_Builder::add_intersection existing internal point x: " << 
                    (*it)->get_x() << " y: " << (*it)->get_y() << "\n";
            cout.flush();
        }
        cout << "Intersect_Meshes_2D::Facet_Builder::add_intersection Current p1p2_pts:\n";
        cout.flush();
        for (vector<shared_ptr<Point_2D>>::const_iterator it = p1p2_pts.begin(); it != p1p2_pts.end(); ++it)
        {
            cout << "Intersect_Meshes_2D::Facet_Builder::add_intersection existing p1p2 point x: " << 
                    (*it)->get_x() << " y: " << (*it)->get_y() << "\n";
            cout.flush();
        }
        cout << "Intersect_Meshes_2D::Facet_Builder::add_intersection Current p1p3_pts:\n";
        cout.flush();
        for (vector<shared_ptr<Point_2D>>::const_iterator it = p1p3_pts.begin(); it != p1p3_pts.end(); ++it)
        {
            cout << "Intersect_Meshes_2D::Facet_Builder::add_intersection existing p1p3 point x: " << 
                    (*it)->get_x() << " y: " << (*it)->get_y() << "\n";
            cout.flush();
        }
        cout << "Intersect_Meshes_2D::Facet_Builder::add_intersection Current p2p3_pts:\n";
        cout.flush();
        for (vector<shared_ptr<Point_2D>>::const_iterator it = p2p3_pts.begin(); it != p2p3_pts.end(); ++it)
        {
            cout << "Intersect_Meshes_2D::Facet_Builder::add_intersection existing p2p3 point x: " << 
                    (*it)->get_x() << " y: " << (*it)->get_y() << "\n";
            cout.flush();
        }
#endif
        
        Intersecting_Facet_Side_Pts i_side_pts;
        for (I_Pt_List::iterator ip_it = intersect_pts.begin(); ip_it != intersect_pts.end(); ++ip_it)
        {
            Location loc(for_facet1 ? ip_it->f1_loc : ip_it->f2_loc); // the location on this facet
            Location i_loc(for_facet1 ? ip_it->f2_loc : ip_it->f1_loc); // the location on the intersecting facet
            
            if (loc == Location::internal)
            {
                vector<shared_ptr<Point_2D>>::const_iterator internal_it = find_if(internal_pts.begin(), internal_pts.end(), Point_Find(ip_it->pt));
                if (internal_it == internal_pts.end())
                {
                    internal_pts.push_back(ip_it->pt);
#ifdef DEBUG_INTERSECT_MESHES_2D_FACET_BUILDER
                    cout << "Intersect_Meshes_2D::Facet_Builder::add_intersection adding new internal point x: " << 
                            ip_it->pt->get_x() << " y: " << ip_it->pt->get_y() << "\n";
                    cout.flush();
#endif
                }
                else // update intersect point to use this point instead
                    ip_it->pt = *internal_it;
            }
            else if (loc == Location::p1p2)
            {
                vector<shared_ptr<Point_2D>>::const_iterator side_it = find_if(p1p2_pts.begin(), p1p2_pts.end(), Point_Find(ip_it->pt));
                if (side_it == p1p2_pts.end())
                {
                    p1p2_pts.push_back(ip_it->pt);
#ifdef DEBUG_INTERSECT_MESHES_2D_FACET_BUILDER
                    cout << "Intersect_Meshes_2D::Facet_Builder::add_intersection adding new p1p2 point x: " << 
                            ip_it->pt->get_x() << " y: " << ip_it->pt->get_y() << "\n";
                    cout.flush();
#endif
                }
                else // update intersect point to use this point instead
                    ip_it->pt = *side_it;
            }
            else if (loc == Location::p1p3)
            {
                vector<shared_ptr<Point_2D>>::const_iterator side_it = find_if(p1p3_pts.begin(), p1p3_pts.end(), Point_Find(ip_it->pt));
                if (side_it == p1p3_pts.end())
                {
                    p1p3_pts.push_back(ip_it->pt);
#ifdef DEBUG_INTERSECT_MESHES_2D_FACET_BUILDER
                    cout << "Intersect_Meshes_2D::Facet_Builder::add_intersection adding new p1p3 point x: " << 
                            ip_it->pt->get_x() << " y: " << ip_it->pt->get_y() << "\n";
                    cout.flush();
#endif
                }
                else // update intersect point to use this point instead
                    ip_it->pt = *side_it;
            }
            else if (loc == Location::p2p3)
            {
                vector<shared_ptr<Point_2D>>::const_iterator side_it = find_if(p2p3_pts.begin(), p2p3_pts.end(), Point_Find(ip_it->pt));
                if (side_it == p2p3_pts.end())
                {
                    p2p3_pts.push_back(ip_it->pt);
#ifdef DEBUG_INTERSECT_MESHES_2D_FACET_BUILDER
                    cout << "Intersect_Meshes_2D::Facet_Builder::add_intersection adding new p2p3 point x: " << 
                            ip_it->pt->get_x() << " y: " << ip_it->pt->get_y() << "\n";
                    cout.flush();
#endif
                }
                else // update intersect point to use this point instead
                    ip_it->pt = *side_it;
            }
            else if (loc == Location::p1)
                ip_it->pt = orig_facet.get_point1();
            else if (loc == Location::p2)
                ip_it->pt = orig_facet.get_point2();
            else // (loc == Location::p3)
                ip_it->pt = orig_facet.get_point3();
            
            switch (i_loc)
            {
                case (Location::p1):
                    i_side_pts.p1p2_points.push_back(*ip_it);
                    i_side_pts.p1p3_points.push_back(*ip_it);
                    break;
                case (Location::p2):
                    i_side_pts.p1p2_points.push_back(*ip_it);
                    i_side_pts.p2p3_points.push_back(*ip_it);
                    break;
                case (Location::p3):
                    i_side_pts.p1p3_points.push_back(*ip_it);
                    i_side_pts.p2p3_points.push_back(*ip_it);
                    break;
                case (Location::p1p2):
                    i_side_pts.p1p2_points.push_back(*ip_it);
                    break;
                case (Location::p1p3):
                    i_side_pts.p1p3_points.push_back(*ip_it);
                    break;
                case (Location::p2p3):
                    i_side_pts.p2p3_points.push_back(*ip_it);
                    break;
                default:
                    break;
            }
            
        }

#ifdef DEBUG_INTERSECT_MESHES_2D_FACET_BUILDER
        cout << "Intersect_Meshes_2D::Facet_Builder::add_intersection internal_points after processing:\n";
        cout.flush();
        for (vector<shared_ptr<Point_2D>>::const_iterator it = internal_pts.begin(); it != internal_pts.end(); ++it)
        {
            cout << "Intersect_Meshes_2D::Facet_Builder::add_intersection current internal point x: " << 
                    (*it)->get_x() << " y: " << (*it)->get_y() << "\n";
            cout.flush();
        }
        cout << "Intersect_Meshes_2D::Facet_Builder::add_intersection p1p2_pts after processing:\n";
        cout.flush();
        for (vector<shared_ptr<Point_2D>>::const_iterator it = p1p2_pts.begin(); it != p1p2_pts.end(); ++it)
        {
            cout << "Intersect_Meshes_2D::Facet_Builder::add_intersection current p1p2 point x: " << 
                    (*it)->get_x() << " y: " << (*it)->get_y() << "\n";
            cout.flush();
        }
        cout << "Intersect_Meshes_2D::Facet_Builder::add_intersection p1p3_pts after processing:\n";
        cout.flush();
        for (vector<shared_ptr<Point_2D>>::const_iterator it = p1p3_pts.begin(); it != p1p3_pts.end(); ++it)
        {
            cout << "Intersect_Meshes_2D::Facet_Builder::add_intersection current p1p3 point x: " << 
                    (*it)->get_x() << " y: " << (*it)->get_y() << "\n";
            cout.flush();
        }
        cout << "Intersect_Meshes_2D::Facet_Builder::add_intersection p2p3_pts after processing:\n";
        cout.flush();
        for (vector<shared_ptr<Point_2D>>::const_iterator it = p2p3_pts.begin(); it != p2p3_pts.end(); ++it)
        {
            cout << "Intersect_Meshes_2D::Facet_Builder::add_intersection current p2p3 point x: " << 
                    (*it)->get_x() << " y: " << (*it)->get_y() << "\n";
            cout.flush();
        }
        cout << "Intersect_Meshes_2D::Facet_Builder::add_intersection generating internal segments\n";
        cout.flush();
#endif
        
        // now generate internal segments
        this->gen_internal_segs(intersect_pts, i_side_pts);
    }
    
    const Intersect_Meshes_2D::Facet_Builder::Line_Segment* Intersect_Meshes_2D::Facet_Builder::create_link_segment(const shared_ptr<Point_2D> pt)
    {
        Line_Segment seg(pt, orig_facet.get_point1(), Location::internal);
        if (segments.does_seg_intersect(seg, orig_facet, p1p2_pts, p1p3_pts, p2p3_pts, internal_pts, precision))
        {
            seg.point2 = orig_facet.get_point2();
            if (segments.does_seg_intersect(seg, orig_facet, p1p2_pts, p1p3_pts, p2p3_pts, internal_pts, precision))
            {
                seg.point2 = orig_facet.get_point3();
                if (segments.does_seg_intersect(seg, orig_facet, p1p2_pts, p1p3_pts, p2p3_pts, internal_pts, precision))
                {
                    // if none of the corner points can link to the point, then try for side points
                    for (vector<shared_ptr<Point_2D>>::const_iterator side_pt_it = p1p2_pts.begin(); side_pt_it != p1p2_pts.end(); ++side_pt_it)
                    {
                        seg.point2 = *side_pt_it;
                        if (!segments.does_seg_intersect(seg, orig_facet, p1p2_pts, p1p3_pts, p2p3_pts, internal_pts, precision))
                        {
                            segments.add_internal_segment(seg.point1, seg.point2);
                            return &*find(segments.begin(), segments.end(), Line_Segment(seg.point1, seg.point2, Location::internal));
//                            return &*(--segments.internal_segs_end());
                        }
                    }
                    for (vector<shared_ptr<Point_2D>>::const_iterator side_pt_it = p1p3_pts.begin(); side_pt_it != p1p3_pts.end(); ++side_pt_it)
                    {
                        seg.point2 = *side_pt_it;
                        if (!segments.does_seg_intersect(seg, orig_facet, p1p2_pts, p1p3_pts, p2p3_pts, internal_pts, precision))
                        {
                            segments.add_internal_segment(seg.point1, seg.point2);
                            return &*find(segments.begin(), segments.end(), Line_Segment(seg.point1, seg.point2, Location::internal));
//                            return &*(--segments.internal_segs_end());
                        }
                    }
                    for (vector<shared_ptr<Point_2D>>::const_iterator side_pt_it = p2p3_pts.begin(); side_pt_it != p2p3_pts.end(); ++side_pt_it)
                    {
                        seg.point2 = *side_pt_it;
                        if (!segments.does_seg_intersect(seg, orig_facet, p1p2_pts, p1p3_pts, p2p3_pts, internal_pts, precision))
                        {
                            segments.add_internal_segment(seg.point1, seg.point2);
                            return &*find(segments.begin(), segments.end(), Line_Segment(seg.point1, seg.point2, Location::internal));
//                            return &*(--segments.internal_segs_end());
                        }
                    }
                    // try other internal points last chance
                    for (vector<shared_ptr<Point_2D>>::const_iterator internal_pt_it = p1p2_pts.begin(); internal_pt_it != p1p2_pts.end(); ++internal_pt_it)
                    {
                        if ((*internal_pt_it)->get_x() == pt->get_x() && (*internal_pt_it)->get_y() == pt->get_y()) // don't process the same internal point
                            continue;
                        seg.point2 = *internal_pt_it;
                        if (!segments.does_seg_intersect(seg, orig_facet, p1p2_pts, p1p3_pts, p2p3_pts, internal_pts, precision))
                        {
                            segments.add_internal_segment(seg.point1, seg.point2);
                            return &*find(segments.begin(), segments.end(), Line_Segment(seg.point1, seg.point2, Location::internal));
//                            return &*(--segments.internal_segs_end());
                        }
                    }
                    // else if no link segment could be created,
                    return 0;
                }
                else
                    segments.add_internal_segment(seg.point1, seg.point2);
            }
            else
                segments.add_internal_segment(seg.point1, seg.point2);
        }
        else
            segments.add_internal_segment(seg.point1, seg.point2);
        return &*find(segments.begin(), segments.end(), Line_Segment(seg.point1, seg.point2, Location::internal));
//        return &*(--segments.internal_segs_end());
    }
    
    void Intersect_Meshes_2D::Facet_Builder::check_internal_pts()
    {
        // go through internal points and look for any not in a segment
        for (vector<shared_ptr<Point_2D>>::const_iterator it = internal_pts.begin(); it != internal_pts.end(); ++it)
        {
            bool found(false);
            for (Segments::const_iterator seg_it = segments.begin(); seg_it != segments.end(); ++seg_it)
            {
                if (seg_it->location != Location::internal)
                    break;
                if (seg_it->point1 == *it || seg_it->point2 == *it)
                {
                    found = true;
                    break;
                }
            }
            if (!found) // create a link segment to the point
            {
#ifdef DEBUG_INTERSECT_MESHES_2D_FACET_BUILDER
                cout << "Intersect_Meshes_2D::Facet_Builder::check_internal_pts found isolated intersect point. Creating link segment to internal point x: " << 
                        (*it)->get_x() << " y: " << (*it)->get_y() << "\n";
#endif
                const Line_Segment* link_seg = create_link_segment(*it);
#ifdef DEBUG_INTERSECT_MESHES_2D_FACET_BUILDER
                if (link_seg == 0)
                    cout << "Intersect_Meshes_2D::Facet_Builder::check_internal_pts Could not create link segment to isolated internal point x: " << 
                            (*it)->get_x() << " y: " << (*it)->get_y() << "\n";
                else
                    cout << "Intersect_Meshes_2D::Facet_Builder::check_internal_pts created link segment p1 x: " << 
                            link_seg->point1->get_x() << " y: " << link_seg->point1->get_y() << " p2 x: " << 
                            link_seg->point2->get_x() << " y: " << link_seg->point2->get_y() << "\n";
#endif
                if (link_seg == 0)
                    throw runtime_error("Unable to generate link segment for internal point");
            }
        }
    }
    
    const bool Intersect_Meshes_2D::Facet_Builder::check_internal_seg(const Line_Segment& seg)
    {
        // look if points are internal points
        bool p1_found(false); // found matching internal point for p1
        bool p2_found(false); // found matching internal point for p2
        for (vector<shared_ptr<Point_2D>>::const_iterator internal_pt_it = internal_pts.begin(); internal_pt_it != internal_pts.end(); ++internal_pt_it)
        {
            if (seg.point1 == *internal_pt_it)
                p1_found = true;
            else if (seg.point2 == *internal_pt_it)
                p2_found = true;
        }
        return !(p1_found && p2_found);
    }
    
    const bool Intersect_Meshes_2D::Facet_Builder::find_path(const Segments::const_iterator& current_it, 
            const shared_ptr<Point_2D>& point, vector<Line_Segment>& path, 
            vector<Line_Segment>& checked_segs)
    {
        // create a stack to process all the connecting segments found on the path
        stack<shared_ptr<Point_2D>> pt_stack; // non shared point
        // create path from point
        Segments::const_iterator prev_seg = current_it;
        shared_ptr<Point_2D> pt(point);
        Segments::const_iterator it = find_if(prev_seg + 1, segments.end(), Segment_Find(pt));
        while (it != segments.end() && it->location == Location::internal)
        {
            // check if segment is already in the path - forms a circle
            if (it == prev_seg || path.end() != find(path.begin(), path.end(), *it))
            {
                // move to next segment
                it = find_if(it + 1, segments.end(), Segment_Find(pt));
                if (it == segments.end() || it->location != Location::internal)
                {
                    // look for a segment in the stack
                    if (!pt_stack.empty())
                    {
                        pt = pt_stack.top();
                        pt_stack.pop();
                        it = find_if(segments.begin(), segments.end(), Segment_Find(pt));
                    }
                }
            }
            else if (checked_segs.end() != find(checked_segs.begin(), checked_segs.end(), *it)) // check if segment is a checked segment
                return true; // found segment path to corner or side of facet
            else
            {
                // add line segment to path
                path.push_back(*it);
                if (check_internal_seg(*it)) // check if the segment does connect with the perimeter
                    return true;
                else
                {
                    // add non shared point to stack
                    pt_stack.push((pt == it->point1) ? it->point2 : it->point1);
                    // get next connecting segment
                    it = find_if(it + 1, segments.end(), Segment_Find(pt));
                    if (it == segments.end() || it->location != Location::internal) // no more connecting segments
                    {
                        // look for a segment in the stack
                        if (!pt_stack.empty())
                        {
                            pt = pt_stack.top();
                            pt_stack.pop();
                            it = find_if(segments.begin(), segments.end(), Segment_Find(pt));
                        }
                    }
                }
            }
        }
        return false;
    }

    const bool Intersect_Meshes_2D::Facet_Builder::complete_path(vector<Line_Segment>& path, 
            vector<Line_Segment>& checked_segs, const Line_Segment& segment)
    {
        bool found(false);
        // create a link segment from point1
        const Line_Segment* link_seg(create_link_segment(segment.point1));
        // check link segment
        if (link_seg != 0)
        {
            if (check_internal_seg(*link_seg))
            {
                checked_segs.push_back(*link_seg);
                checked_segs.push_back(segment);
                // add all path elements to checked_segs
                for (vector<Line_Segment>::const_iterator path_it = path.begin(); path_it != path.end(); ++path_it)
                    checked_segs.push_back(*path_it);
                return true;
            } 
            else
                path.push_back(*link_seg);
        }
        // try from point2
        link_seg = create_link_segment(segment.point2);
        // check link segment
        if (link_seg != 0)
        {
            if (check_internal_seg(*link_seg))
            {
                checked_segs.push_back(*link_seg);
                checked_segs.push_back(segment);
                // add all path elements to checked_segs
                for (vector<Line_Segment>::const_iterator path_it = path.begin(); path_it != path.end(); ++path_it)
                    checked_segs.push_back(*path_it);
                return true;
            }
            else
                path.push_back(*link_seg);
        }
        return false;
    }
    
    void Intersect_Meshes_2D::Facet_Builder::verify_internal_paths()
    {
        vector<Line_Segment> checked_segs;
        // Go through internal segments and verify that they are connected to
        // a side or corner point
        Segments::const_iterator seg_it = segments.begin();
        while (seg_it != segments.end() && seg_it->location == Location::internal)
        {
#ifdef DEBUG_INTERSECT_MESHES_2D_FACET_BUILDER
            cout << "Intersect_Meshes_2D::Facet_Builder::verify_internal_segs checking internal segment p1 x: " << 
                    seg_it->point1->get_x() << " y: " << seg_it->point1->get_y() << " p2 x: " << 
                    seg_it->point2->get_x() << " y: " << seg_it->point2->get_y() << "\n";
#endif
            if (checked_segs.end() != find(checked_segs.begin(), checked_segs.end(), *seg_it))
            {
#ifdef DEBUG_INTERSECT_MESHES_2D_FACET_BUILDER
                cout << "Intersect_Meshes_2D::Facet_Builder::verify_internal_segs segment has already been verified\n";
#endif
                ++seg_it;
                continue; // already checked segment, so go to next one
            }
            
            if (check_internal_seg(*seg_it)) // if one of the points is not an internal point
            {
#ifdef DEBUG_INTERSECT_MESHES_2D_FACET_BUILDER
                cout << "Intersect_Meshes_2D::Facet_Builder::verify_internal_segs segment contains a perimeter point\n";
#endif
                checked_segs.push_back(*seg_it); // segment connects to facet perimeter
            }
            else // segment has both points inside the facet
            {
#ifdef DEBUG_INTERSECT_MESHES_2D_FACET_BUILDER
                cout << "Intersect_Meshes_2D::Facet_Builder::verify_internal_segs segment point1 and point2 are internal points\n";
#endif
                // segment is fully inside facet, try locating linking segments on either side of the segment and follow
                // the trail to a side or corner point
                
                // check if there is an already processed segment that shares point1
                vector<Line_Segment>::const_iterator it = find_if(checked_segs.begin(), checked_segs.end(), Segment_Find(seg_it->point1));
                if (it != checked_segs.end())
                {
#ifdef DEBUG_INTERSECT_MESHES_2D_FACET_BUILDER
                    cout << "Intersect_Meshes_2D::Facet_Builder::verify_internal_segs segment shares point1 with a segment connected to the perimeter\n";
#endif
                    checked_segs.push_back(*seg_it);
                    ++seg_it;
                    continue;
                }
                // check if there is an already processed segment that shares point2
                it = find_if(checked_segs.begin(), checked_segs.end(), Segment_Find(seg_it->point2));
                if (it != checked_segs.end())
                {
#ifdef DEBUG_INTERSECT_MESHES_2D_FACET_BUILDER
                    cout << "Intersect_Meshes_2D::Facet_Builder::verify_internal_segs segment shares point2 with a segment connected to the perimeter\n";
#endif
                    checked_segs.push_back(*seg_it);
                    ++seg_it;
                    continue;
                }
#ifdef DEBUG_INTERSECT_MESHES_2D_FACET_BUILDER
                cout << "Intersect_Meshes_2D::Facet_Builder::verify_internal_segs Finding path to perimeter\n";
#endif
                // no connecting segment linked to a side or corner point was found
                // create a stack of connecting segments looking for a segment that is connected to 
                vector<Line_Segment> path;
                path.push_back(*seg_it);
                if (find_path(seg_it, seg_it->point1, path, checked_segs))
                {
#ifdef DEBUG_INTERSECT_MESHES_2D_FACET_BUILDER
                    cout << "Intersect_Meshes_2D::Facet_Builder::verify_internal_segs Found path to perimeter from point1.  Adding path to checked_segs\n";
#endif
                    // add all path elements to checked_segs
                    for (vector<Line_Segment>::const_iterator path_it = path.begin(); path_it != path.end(); ++path_it)
                        checked_segs.push_back(*path_it);
                }
                else if (find_path(seg_it, seg_it->point2, path, checked_segs))
                {
#ifdef DEBUG_INTERSECT_MESHES_2D_FACET_BUILDER
                    cout << "Intersect_Meshes_2D::Facet_Builder::verify_internal_segs Found path to perimeter from point2.  Adding path to checked_segs\n";
#endif
                    // add all path elements to checked_segs
                    for (vector<Line_Segment>::const_iterator path_it = path.begin(); path_it != path.end(); ++path_it)
                        checked_segs.push_back(*path_it);
                }
                else // create a link segment
                {
#ifdef DEBUG_INTERSECT_MESHES_2D_FACET_BUILDER
                    cout << "Intersect_Meshes_2D::Facet_Builder::verify_internal_segs Could not find path to perimeter.  Generating a link segment\n";
#endif
                    // try to link segment with perimeter of facet
                    if (this->complete_path(path, checked_segs, *path.begin()))
                    {
#ifdef DEBUG_INTERSECT_MESHES_2D_FACET_BUILDER
                        cout << "Intersect_Meshes_2D::Facet_Builder::verify_internal_segs Generated a link segment to the perimeter\n";
#endif
                        seg_it = find(segments.begin(), segments.end(), *path.begin());
                        ++seg_it;
                        continue;
                    }
                    else if (path.size() > 1)
                    {
#ifdef DEBUG_INTERSECT_MESHES_2D_FACET_BUILDER
                        cout << "Intersect_Meshes_2D::Facet_Builder::verify_internal_segs Unable to generate a link segment to the perimeter\n";
#endif
                        // try from other segments in the path
                        bool found(false);
                        vector<Line_Segment>::const_iterator path_it(++path.begin());
                        while (path_it != path.end())
                        {
#ifdef DEBUG_INTERSECT_MESHES_2D_FACET_BUILDER
                            cout << "Intersect_Meshes_2D::Facet_Builder::verify_internal_segs Attempting to generate a link segment from segment p1 x: " << 
                                    path_it->point1->get_x() << " y: " << path_it->point1->get_y() << " p2 x: " << 
                                    path_it->point2->get_x() << " y: " << path_it->point2->get_y() << "\n";
#endif
                            Line_Segment path_seg(path_it->point1, path_it->point2, path_it->location);
                            if (this->complete_path(path, checked_segs, path_seg))
                            {
                                found = true;
                                break;
                            }
                            else
                            {
                                path_it = find(path.begin(), path.end(), path_seg);
                                ++path_it;
                            }
                        }
                        if (found)
                        {
#ifdef DEBUG_INTERSECT_MESHES_2D_FACET_BUILDER
                            cout << "Intersect_Meshes_2D::Facet_Builder::verify_internal_segs Generated a link segment to the perimeter\n";
#endif
                            seg_it = find(segments.begin(), segments.end(), *path.begin());
                            ++seg_it;
                            continue;
                        }
                        else
                            throw runtime_error("Unable to link segment with facet perimeter");
                    }
                }
            }
            ++seg_it;
        }
    }

    void Intersect_Meshes_2D::Facet_Builder::verify_internal_segs()
    {
        // check each internal segment if it forms a straight line with 
        // another internal segment
        Segments::const_iterator seg_it = segments.begin();
        while (seg_it != segments.end() && seg_it->location == Location::internal)
        {
            // check if segment point1 internal point
            if (internal_pts.end() != find(internal_pts.begin(), internal_pts.end(), seg_it->point1))
            {
                // try to find a connecting segment to the internal point
                vector<Line_Segment> connecting_segs;
                Segments::const_iterator it = find_if(seg_it + 1, segments.end(), Segment_Find(seg_it->point1));
                while (it != segments.end()) // found connecting segment
                {
                    connecting_segs.push_back(*it);
                    it = find_if(it + 1, segments.end(), Segment_Find(seg_it->point1));
                }
                if (connecting_segs.size() == 1)
                {
                    bool same_direction(false);
                    if (is_same_line(*seg_it->point1, *seg_it->point2, 
                            *connecting_segs.begin()->point1, *connecting_segs.begin()->point2, 
                            same_direction, precision))
                    {
                        Line_Segment seg(*seg_it);
                        // segments are in a straight line
                        // add a link segment to the common point
                        create_link_segment(seg_it->point1);
                        seg_it = find(segments.begin(), segments.end(), seg);
                    }
                }
            }
            // check if segment point2 internal point
            if (internal_pts.end() != find(internal_pts.begin(), internal_pts.end(), seg_it->point2))
            {
                // try to find a connecting segment to the internal point
                vector<Line_Segment> connecting_segs;
                Segments::const_iterator it = find_if(seg_it + 1, segments.end(), Segment_Find(seg_it->point2));
                while (it != segments.end()) // found connecting segment
                {
                    connecting_segs.push_back(*it);
                    it = find_if(it + 1, segments.end(), Segment_Find(seg_it->point1));
                }
                if (connecting_segs.size() == 1)
                {
                    bool same_direction(false);
                    if (is_same_line(*seg_it->point1, *seg_it->point2, 
                            *connecting_segs.begin()->point1, *connecting_segs.begin()->point2, 
                            same_direction, precision))
                    {
                        Line_Segment seg(*seg_it);
                        // segments are in a straight line
                        // add a link segment to the common point
                        create_link_segment(seg_it->point2);
                        seg_it = find(segments.begin(), segments.end(), seg);
                    }
                }
            }
            ++seg_it;
        }
    }
    
    void Intersect_Meshes_2D::Facet_Builder::gen_side_segs(vector<shared_ptr<Point_2D>>& side_points, 
            const Location side, const shared_ptr<Point_2D>& p1, const shared_ptr<Point_2D>& p2)
    {
#ifdef DEBUG_INTERSECT_MESHES_2D_FACET_BUILDER
        cout << "Intersect_Meshes_2D::Facet_Builder::gen_side_segs begin\n";
        cout << "Intersect_Meshes_2D::Facet_Builder::gen_side_segs p1 x: " << p1->get_x() << " y: " << p1->get_y() << "\n";
        cout << "Intersect_Meshes_2D::Facet_Builder::gen_side_segs p2 x: " << p2->get_x() << " y: " << p2->get_y() << "\n";
        for (vector<shared_ptr<Point_2D>>::const_iterator it = side_points.begin(); it != side_points.end(); ++it)
            cout << "Intersect_Meshes_2D::Facet_Builder::gen_side_segs side_point x: " << (*it)->get_x() << 
                    " y: " << (*it)->get_y() << "\n";
        cout << "Intersect_Meshes_2D::Facet_Builder::gen_side_segs removing corner points\n";
#endif
        
        // remove any corner points
//        vector<shared_ptr<Point_2D>> side_pts;
//        for (I_Pt_List::const_iterator it = side_points.begin(); it != side_points.end(); ++it)
//        {
//            if ((for_facet1 && (*it).f1_loc == side) || 
//                    (!for_facet1 && (*it).f2_loc == side))
//                side_pts.push_back((*it).pt);
//        }
#ifdef DEBUG_INTERSECT_MESHES_2D_FACET_BUILDER
        cout << "Intersect_Meshes_2D::Facet_Builder::gen_side_segs side_points size: " << side_points.size() << "\n";
        for (vector<shared_ptr<Point_2D>>::const_iterator it = side_points.begin(); it != side_points.end(); ++it)
            cout << "Intersect_Meshes_2D::Facet_Builder::gen_side_segs side_points x: " << (*it)->get_x() << " y: " << 
                    (*it)->get_y() << "\n";
#endif        
        
        if (side_points.empty())
        {
#ifdef DEBUG_INTERSECT_MESHES_2D_FACET_BUILDER
            cout << "Intersect_Meshes_2D::Facet_Builder::gen_side_segs adding segment: p1 x: " << p1->get_x() << 
                    " y: " << p1->get_y() << " ptr: " << p1.get() << " p2 x: " << 
                    p2->get_x() << " y: " << p2->get_y() << " ptr: " << p2.get()  << "\n";
            int size = segments.size();
#endif
            segments.add_external_segment(p1, p2, side);
#ifdef DEBUG_INTERSECT_MESHES_2D_FACET_BUILDER
            if (segments.size() > size)
                cout << "Intersect_Meshes_2D::Facet_Builder::gen_side_segs added segment\n";
            else
                cout << "Intersect_Meshes_2D::Facet_Builder::gen_side_segs did not add segment\n";
#endif
        }
        else
        {
            class Point_Sort {
            public:
                Point_Sort(const shared_ptr<Point_2D> cp) : corner_pt(cp) {}
                const bool operator()(const shared_ptr<Point_2D> sp1, const shared_ptr<Point_2D> sp2) const
                {
                    return Vector_2D(*corner_pt, *sp1).length() < Vector_2D(*corner_pt, *sp2).length();
                }
            private:
                const shared_ptr<Point_2D> corner_pt;
            };
            
            sort(side_points.begin(), side_points.end(), Point_Sort(p1));

            shared_ptr<Point_2D> prev_pt(p1); // start with the corner point
            for (vector<shared_ptr<Point_2D>>::const_iterator it = side_points.begin(); it != side_points.end(); ++it)
            {
#ifdef DEBUG_INTERSECT_MESHES_2D_FACET_BUILDER
                cout << "Intersect_Meshes_2D::Facet_Builder::gen_side_segs adding segment: p1 x: " << prev_pt->get_x() << 
                        " y: " << prev_pt->get_y() << " ptr: " << prev_pt << " p2 x: " << 
                        (*it)->get_x() << " y: " << (*it)->get_y() << " ptr: " << (*it) << "\n";
                int size = segments.size();
#endif
                segments.add_external_segment(prev_pt, *it, side);
#ifdef DEBUG_INTERSECT_MESHES_2D_FACET_BUILDER
                if (segments.size() > size)
                    cout << "Intersect_Meshes_2D::Facet_Builder::gen_side_segs added segment\n";
                else
                    cout << "Intersect_Meshes_2D::Facet_Builder::gen_side_segs did not add segment\n";
#endif
                prev_pt = *it;
            }
            // add last segment
#ifdef DEBUG_INTERSECT_MESHES_2D_FACET_BUILDER
            cout << "Intersect_Meshes_2D::Facet_Builder::gen_side_segs add last segment: p1 x: " << prev_pt->get_x() << 
                    " y: " << prev_pt->get_y() << " ptr: " << prev_pt.get() << " p2 x: " << 
                    p2->get_x() << " y: " << p2->get_y() << " ptr: " << p2.get() << "\n";
            int size = segments.size();
#endif
            segments.add_external_segment(prev_pt, p2, side);
#ifdef DEBUG_INTERSECT_MESHES_2D_FACET_BUILDER
            if (segments.size() > size)
                cout << "Intersect_Meshes_2D::Facet_Builder::gen_side_segs added segment\n";
            else
                cout << "Intersect_Meshes_2D::Facet_Builder::gen_side_segs did not add segment\n";
#endif
        }
#ifdef DEBUG_INTERSECT_MESHES_2D_FACET_BUILDER
        cout << "Intersect_Meshes_2D::Facet_Builder::gen_side_segments end\n";
#endif
    }

    const bool Intersect_Meshes_2D::Facet_Builder::contains_internal_pt(
            const shared_ptr<Point_2D>& p1, const shared_ptr<Point_2D>& p2, 
            const shared_ptr<Point_2D>& p3, const Facet_2D& facet)
    {
        for (vector<shared_ptr<Point_2D>>::const_iterator internal_pt_iter = internal_pts.begin(); 
                internal_pt_iter != internal_pts.end(); ++internal_pt_iter)
        {
//            const Intersect_Point* ip = &*internal_point_iter;
            
            bool pt_on_side(false);
            if (!is_equal(**internal_pt_iter, *p1, precision) && 
                    !is_equal(**internal_pt_iter, *p2, precision) && 
                    !is_equal(**internal_pt_iter, *p3, precision) && 
                    facet.contains_point(**internal_pt_iter, pt_on_side, precision))
                return true;
        }
        
        return false;
    }
    
    void Intersect_Meshes_2D::Facet_Builder::build_facets(Facets& new_facets)
    {
#ifdef DEBUG_INTERSECT_MESHES_2D_FACET_BUILDER
        cout << "Intersect_Meshes_2D::Facet_Builder::build_facets orig_facet p1 x: " << 
                orig_facet.get_point1()->get_x() << " y: " << orig_facet.get_point1()->get_y() << 
                " ptr: " << (orig_facet.get_point1().get()) << " p2 x: " << 
                orig_facet.get_point2()->get_x() << " y: " << orig_facet.get_point2()->get_y() << 
                " ptr: " << (orig_facet.get_point2().get()) << " p3 x: " << 
                orig_facet.get_point3()->get_x() << " y: " << orig_facet.get_point3()->get_y() << 
                " ptr: " << (orig_facet.get_point3().get()) << "\n";
        cout << "Intersect_Meshes_2D::Facet_Builder::build_facets: internal_segments:\n";
        cout.flush();
        Segments::const_iterator temp_it(segments.begin());
        while (temp_it != segments.end() && temp_it->location == Location::internal)
        {
            cout << "Intersect_Meshes_2D::Facet_Builder::build_facets segment p1 x: " << temp_it->point1->get_x() << 
                    " y: " << temp_it->point1->get_y() << " ptr: " << 
                    temp_it->point1 << " p2 x: " << temp_it->point2->get_x() << " y: " << temp_it->point2->get_y() << 
                    " ptr: " << temp_it->point2 << " used: " << temp_it->used << "\n";
            ++temp_it;
        }
        cout << "Intersect_Meshes_2D::Facet_Builder::build_facets: external_segments:\n";
        cout.flush();
        while (temp_it != segments.end())
        {
            cout << "Intersect_Meshes_2D::Facet_Builder::build_facets segment p1 x: " << temp_it->point1->get_x() << 
                    " y: " << temp_it->point1->get_y() << " p2 x: " << temp_it->point2->get_x() << 
                    " y: " << temp_it->point2->get_y() << " loc: " << temp_it->location << "\n";
            ++temp_it;
        }
        cout.flush();
        for (vector<shared_ptr<Point_2D>>::const_iterator it = internal_pts.begin(); it != internal_pts.end(); ++it)
        {
            cout << "Intersect_Meshes_2D::Facet_Builder::build_facets internal point x: " << (*it)->get_x() << 
                    " y: " << (*it)->get_y() << "\n";
        }
        cout.flush();
        cout << "Intersect_Meshes_2D::Facet_Builder::build_facets begin\n";
        cout.flush();
#endif
        const Line_Segment* segp = segments.get_next_segment(0); // get initial segment
        while (segp != 0)
        {
            Line_Segment segment1(*segp);
#ifdef DEBUG_INTERSECT_MESHES_2D_FACET_BUILDER
            cout << "Intersect_Meshes_2D::Facet_Builder::build_facets segment1 p1 x=" << 
                    segment1.point1->get_x() << ", y=" << segment1.point1->get_y() << 
                    " ptr: " << segment1.point1 << " p2(x=" << segment1.point2->get_x() << 
                    ", y=" << segment1.point2->get_y() << " ptr: " << segment1.point2 << "\n";
            cout.flush();
#endif
            shared_ptr<Point_2D> shared_pt;
            shared_ptr<Point_2D> p2;
            shared_ptr<Point_2D> p3;
            
            segp = segments.find_connecting_segment(segment1, 0, shared_pt, precision);
//            cout << "segp=" << segp << "\n";
            cout.flush();
            while (segp != 0)
            {
                Line_Segment segment2(*segp);
#ifdef DEBUG_INTERSECT_MESHES_2D_FACET_BUILDER
                cout << "Intersect_Meshes_2D::Facet_Builder::build_facets segment2 p1 x=" << 
                        segment2.point1->get_x() << ", y=" << segment2.point1->get_y() << 
                        " ptr: " << segment2.point2 << " p2 x=" << segment2.point2->get_x() << 
                        ", y=" << segment2.point2->get_y() << " ptr: " << segment2.point2 << "\n";
                cout.flush();
#endif

                // shared_pt is p1
                p2 = segment1.point1 == shared_pt ? segment1.point2 : segment1.point1;
                p3 = segment2.point1 == shared_pt ? segment2.point2 : segment2.point1;

                // do not know location now, so use invalid segment location p1
                Line_Segment seg3(p2, p3, Location::p1);
#ifdef DEBUG_INTERSECT_MESHES_2D_FACET_BUILDER
                cout << "Intersect_Meshes_2D::Facet_Builder::build_facets segment3 p1 x=" << 
                        seg3.point1->get_x() << ", y=" << seg3.point1->get_y() << " ptr: " << 
                        seg3.point1 << " p2 x=" << seg3.point2->get_x() << ", y=" << 
                        seg3.point2->get_y() << " ptr: " << seg3.point2 << "\n";
                cout << "Intersect_Meshes_2D::Facet_Builder::build_facets checking if segment3 intersects any other segments\n";
                cout.flush();
#endif
                if (segments.does_seg_intersect(seg3, orig_facet, p1p2_pts, p1p3_pts, p2p3_pts, internal_pts, precision))
                {
#ifdef DEBUG_INTERSECT_MESHES_2D_FACET_BUILDER
                    cout << "Intersect_Meshes_2D::Facet_Builder::build_facets segment3 intersects another segment\n";
                    cout.flush();
#endif
                    segp = segments.find_connecting_segment(segment1, segp, shared_pt, precision);
                    continue;
                }
                
                // create facet
#ifdef DEBUG_INTERSECT_MESHES_2D_FACET_BUILDER
                cout << "Intersect_Meshes_2D::Facet_Builder::build_facets forming facet p1(x=" << 
                        shared_pt->get_x() << ", y=" << shared_pt->get_y() << "), p2(x=" << p2->get_x() << 
                        ", y=" << p2->get_y() << "), p3(x=" << p3->get_x() << ", y=" << p3->get_y() << ")\n";
                cout.flush();
#endif
                Facet_2D facet(shared_pt, p2, p3);

#ifdef DEBUG_INTERSECT_MESHES_2D_FACET_BUILDER
                cout << "Intersect_Meshes_2D::Facet_Builder::build_facets checking if facet has already been created\n";
                cout.flush();
#endif
                // make sure facet unit normal is pointing in the same direction
                // as orig_facet
//                if (dot_product(facet.get_unv(), orig_facet.get_unv()) < 0)
//                {
//                    facet.invert_unv();
//                    swap(p2, p3);
//                }
                
                // check if facet already exists
                if (new_facets.contains(facet)) // new_facets.end() != find_if(new_facets.begin(), new_facets.end(), Facet_find(facet, precision)))
                {
#ifdef DEBUG_INTERSECT_MESHES_2D_FACET_BUILDER
                    cout << "Intersect_Meshes_2D::Facet_Builder::build_facets facet has already been created\n";
                    cout.flush();
#endif
                    segp = segments.find_connecting_segment(segment1, segp, shared_pt, precision);
                    continue;
                }

#ifdef DEBUG_INTERSECT_MESHES_2D_FACET_BUILDER
                cout << "Intersect_Meshes_2D::Facet_Builder::build_facets checking if facet contains an internal point\n";
                cout.flush();
#endif
                // check if facet contains an internal point
                if (contains_internal_pt(shared_pt, p2, p3, facet))
                {
#ifdef DEBUG_INTERSECT_MESHES_2D_FACET_BUILDER
                    cout << "Intersect_Meshes_2D::Facet_Builder::build_facets facet contains an internal point\n";
                    cout.flush();
#endif
                    segp = segments.find_connecting_segment(segment1, segp, shared_pt, precision);
                }
                else
                {
#ifdef DEBUG_INTERSECT_MESHES_2D_FACET_BUILDER
                    cout << "Intersect_Meshes_2D::Facet_Builder::build_facets facet does not contain an internal point.  found facet to add\n";
                    cout.flush();
#endif
                    break;
                }
            }
            
            if (segp == 0)
            {
                segp = segments.get_next_segment(&segment1);
                if (segp != 0 && *segp == segment1)
                    throw runtime_error("next segment is the same segment");
                continue;
            }
            
            // update / remove segments as necessary
            Line_Segment segment2(*segp);
            segments.process_used_segments(segment1, segment2);

#ifdef DEBUG_INTERSECT_MESHES_2D_FACET_BUILDER
            cout << "Intersect_Meshes_2D::Facet_Builder::build_facets adding facet: " << new_facets.size() << "\n";
            cout.flush();
#endif
            new_facets.push_back(Facet_2D(shared_pt, p2, p3));
#ifdef DEBUG_INTERSECT_MESHES_2D_FACET_BUILDER
            cout << "Intersect_Meshes_2D::Facet_Builder::build_facets added facet: " << new_facets.size() << "\n";
            cout.flush();
#endif
            segp = segments.get_next_segment(0);
        }
#ifdef DEBUG_INTERSECT_MESHES_2D_FACET_BUILDER
        cout << "Intersect_Meshes_2D::Facet_Builder::build_facets end\n";
        cout.flush();
#endif
    }
    
    const Point_2D::Measurement Intersect_Meshes_2D::Facet_Builder::facet_area(const Facet_2D& facet)
    {
        Vector_2D p1p2(*facet.get_point1(), *facet.get_point2());
        Vector_2D p1p3(*facet.get_point1(), *facet.get_point3());
        Vector_2D p2p3(*facet.get_point2(), *facet.get_point3());
        
        if (p1p2.length() >= p1p3.length() && p1p2.length() > p2p3.length())
        {
            // p1p2 is the hypotenuse
            if (p1p3.length() > p2p3.length())
            {
                // use p1p3 is the base
                double angle = angle_between(p1p2,p1p3);
                double height = p1p2.length() * sin(angle);
                return 0.5 * height * p1p3.length();
            }
            else
            {
                // use p2p3 as the base
                double angle = angle_between(-p1p2,p2p3);
                double height = p1p2.length() * sin(angle);
                return 0.5 * height * p2p3.length();
            }
        }
        else if (p1p3.length() >= p1p2.length() && p1p3.length() >= p2p3.length())
        {
            // p1p3 is the hypotenuse
            if (p1p2.length() > p2p3.length())
            {
                // use p1p2 as the base
                double angle = angle_between(p1p2,p1p3);
                double height = p1p3.length() * sin(angle);
                return 0.5 * height * p1p2.length();
            }
            else
            {
                // use p2p3 as the base
                double angle = angle_between(-p1p3,-p2p3);
                double height = p1p3.length() * sin(angle);
                return 0.5 * height * p2p3.length();
            }
        }
        else
        {
            // p2p3 is the hypotenuse
            if (p1p2.length() > p1p3.length())
            {
                // use p1p2 as the base
                double angle = angle_between(-p1p2,p2p3);
                double height = p2p3.length() * sin(angle);
                return 0.5 * height * p1p2.length();
            }
            else
            {
                // use p1p3 as the base
                double angle = angle_between(-p1p3,-p2p3);
                double height = p2p3.length() * sin(angle);
                return 0.5 * height * p1p3.length();
            }
        }
    }
    
    const bool Intersect_Meshes_2D::Facet_Builder::form_new_facets(Facets& new_facets)
    {
        // if there are internal points
        if (!internal_pts.empty())
        {
            // check for any link segments that may need to be added for lone internal points
            check_internal_pts();

            // check if any link facets need to be added for internal segments
            // not connected to the facet perimeter
            verify_internal_paths();
            
            // check for any internal segments that form a straight line
            // and add a link segment to allow facets to be formed
            verify_internal_segs();
        }
        
        // generate perimeter line segments (external)
        shared_ptr<Point_2D> p1(orig_facet.get_point1());
        shared_ptr<Point_2D> p2(orig_facet.get_point2());
        shared_ptr<Point_2D> p3(orig_facet.get_point3());
        
        // generate side segments
#ifdef DEBUG_INTERSECT_MESHES_2D_FACET_BUILDER
        cout << "Intersect_Meshes_2D::Facet_Builder::form_new_facets creating side segments for p1p2\n";
#endif
        gen_side_segs(p1p2_pts, Location::p1p2, p1, p2);
#ifdef DEBUG_INTERSECT_MESHES_2D_FACET_BUILDER
        cout << "Intersect_Meshes_2D::Facet_Builder::form_new_facets creating side segments for p1p3\n";
#endif
        gen_side_segs(p1p3_pts, Location::p1p3, p1, p3);
#ifdef DEBUG_INTERSECT_MESHES_2D_FACET_BUILDER
        cout << "Intersect_Meshes_2D::Facet_Builder::form_new_facets creating side segments for p2p3\n";
#endif
        gen_side_segs(p2p3_pts, Location::p2p3, p2, p3);
        
        if (segments.size() > 3)
        {
#ifdef DEBUG_INTERSECT_MESHES_2D_FACET_BUILDER
            cout << "Intersect_Meshes_2D::Facet_Builder::form_new_facets Forming Facets\n";
#endif
            Facets temp;
            this->build_facets(temp);
            
#ifdef DEBUG_INTERSECT_MESHES_2D_FACET_BUILDER
            cout << "Intersect_Meshes_2D::Facet_Builder::form_new_facets replacing Facet with new facets\n";
            cout << "Intersect_Meshes_2D::Facet_Builder::form_new_facets original Facet p1 x: " << orig_facet.get_point1()->get_x() << 
                    " y: " << orig_facet.get_point1()->get_y() << " p2 x: " << orig_facet.get_point2()->get_x() << 
                    " y: " << orig_facet.get_point2()->get_y() << " p3 x: " << orig_facet.get_point3()->get_x() << 
                    " y: " << orig_facet.get_point3()->get_y() << "\n";
            double orig_area(facet_area(orig_facet));
            double sum(0);
            for (Facets::const_iterator it = temp.begin(); it != temp.end(); ++it)
            {
                shared_ptr<Point_2D> p1(temp.get_point(it->get_p1_index()));
                shared_ptr<Point_2D> p2(temp.get_point(it->get_p2_index()));
                shared_ptr<Point_2D> p3(temp.get_point(it->get_p3_index()));
                cout << "        Intersect_Meshes_2D::Facet_Builder::form_new_facets new Facet p1 x: " << p1->get_x() << 
                        " y: " << p1->get_y() << " p2 x: " << p2->get_x() << " y: " << 
                        p2->get_y() << " p3 x: " << p3->get_x() << " y: " << p3->get_y() << "\n";
                sum += facet_area(Facet_2D(p1,p2,p3));
            }
            if (orig_area - sum > 0.1)
                cout << "MISSING FACETS!!! area difference: " << (orig_area - sum) << "\n";
#endif
            new_facets.replace_all(temp);
            return true;
        }
        else // no intersection
        {
#ifdef DEBUG_INTERSECT_MESHES_2D_FACET_BUILDER
            cout << "Intersect_Meshes_2D::Facet_Builder::form_new_facets no intersection\n";
#endif
            return false;
        }
    }
    
    Intersect_Meshes_2D::Facet_Sorter::Facet_Sorter(const Point_2D::Measurement& prec) : precision(prec),  
            f1_inside_f2(), f2_inside_f1() {}
    
    void Intersect_Meshes_2D::Facet_Sorter::sort(const Facets& facets1, 
            const Facets& facets2)
    {
#ifdef DEBUG_INTERSECT_MESHES_2D_FACET_SORTER
        cout << "Intersect_Meshes_2D::Facet_Sorter::sort begin\n";
#endif
//        // consider all f1 facets inside f2 and remove them if they are found to be outside
//        for (Facets::const_iterator it = facets1.begin(); it != facets1.end(); ++it)
//        {
//            f1_inside_f2.push_back(*it);
//        }
//        
        vector<Facet> unchecked_f2_facets;
        // consider all f2 facets inside f1 and remove them if they are found to be outside
        for (Facets::const_iterator it = facets2.begin(); it != facets2.end(); ++it)
        {
//            f2_inside_f1.push_back(*it);
            unchecked_f2_facets.push_back(*it);
        }
        
#ifdef DEBUG_INTERSECT_MESHES_2D_FACET_SORTER
        cout << "Intersect_Meshes_2D::Facet_Sorter::sort beginning first round\n";
#endif
        for (Facets::const_iterator f1_it = facets1.begin(); f1_it != facets1.end(); ++f1_it)
        {
            const Facet_2D f1(facets1.get_point(f1_it->get_p1_index()), facets1.get_point(f1_it->get_p2_index()), facets1.get_point(f1_it->get_p3_index()));
            const Point_2D f1_ip(f1.get_inside_point());
#ifdef DEBUG_INTERSECT_MESHES_2D_FACET_SORTER
            cout << "Intersect_Meshes_2D::Facet_Sorter::sort sorting facets1 facet (p1: " << 
                    f1_it->get_p1_index() << " p2: " << f1_it->get_p2_index() << " p3: " << 
                    f1_it->get_p3_index() << ") p1 x: " << f1.get_point1()->get_x() << " y: " << 
                    f1.get_point1()->get_y() << " p2 x: " << f1.get_point2()->get_x() << " y: " << 
                    f1.get_point2()->get_y() << " p3 x: " << f1.get_point3()->get_x() << " y: " << 
                    f1.get_point3()->get_y() << " internal point is x: " << f1_ip.get_x() << " y: " << 
                    f1_ip.get_y() << "\n";
#endif
            // determine if point is on or inside the mesh
            for (Facets::const_iterator f2_it = facets2.begin(); f2_it != facets2.end(); ++f2_it)
            {
                const Facet_2D f2(facets2.get_point(f2_it->get_p1_index()), facets2.get_point(f2_it->get_p2_index()), facets2.get_point(f2_it->get_p3_index()));
                // if point is on a facet, return true
                bool pt_on_side(false);
                if (f2.contains_point(f1_ip, pt_on_side, precision))
                {
#ifdef DEBUG_INTERSECT_MESHES_2D_FACET_SORTER
                    cout << "Intersect_Meshes_2D::Facet_Sorter::sort sorting facets1 facet p1: " << 
                            f1_it->get_p1_index() << " p2: " << f1_it->get_p2_index() << " p3: " << 
                            f1_it->get_p3_index() << " internal point is on surface of f2 p1: " << 
                            f2_it->get_p1_index() << " p2: " << f2_it->get_p2_index() << " p3: " << 
                            f2_it->get_p3_index() << "\n";
#endif
                    // both f1 and f2 overlap
                    f1_inside_f2.push_back(*f1_it);
                    f2_inside_f1.push_back(*f2_it);
                    // remove f2 facet from unchecked
                    vector<Facet>::const_iterator it = find(unchecked_f2_facets.begin(), unchecked_f2_facets.end(), *f2_it);
                    if (it != unchecked_f2_facets.end())
                        unchecked_f2_facets.erase(it); // remove from list since it is on the f1 mesh
                    break; // go to next f1
                }
            }
        }
        
#ifdef DEBUG_INTERSECT_MESHES_2D_FACET_SORTER
        cout << "Intersect_Meshes_2D::Facet_Sorter::sort beginning second round with " << unchecked_f2_facets.size() << " facets2 facets to sort\n";
#endif
        for (vector<Facet>::const_iterator f2_it = unchecked_f2_facets.begin(); f2_it != unchecked_f2_facets.end(); ++f2_it)
        {
            const Facet_2D f2(facets2.get_point(f2_it->get_p1_index()), facets2.get_point(f2_it->get_p2_index()), facets2.get_point(f2_it->get_p3_index()));
            const Point_2D f2_ip(f2.get_inside_point());
#ifdef DEBUG_INTERSECT_MESHES_2D_FACET_SORTER
            cout << "Intersect_Meshes_2D::Facet_Sorter::sort sorting facets2 facet (p1: " << 
                    f2_it->get_p1_index() << " p2: " << f2_it->get_p2_index() << " p3: " << 
                    f2_it->get_p3_index() << ") p1 x: " << f2.get_point1()->get_x() << " y: " << 
                    f2.get_point1()->get_y() << " p2 x: " << f2.get_point2()->get_x() << " y: " << 
                    f2.get_point2()->get_y() << " p3 x: " << f2.get_point3()->get_x() << " y: " << 
                    f2.get_point3()->get_y() << " internal point is x: " << f2_ip.get_x() << " y: " << 
                    f2_ip.get_y() << "\n";
#endif
            // determine if point is inside the mesh
            for (Facets::const_iterator f1_it = facets1.begin(); f1_it != facets1.end(); ++f1_it)
            {
                const Facet_2D f1(facets1.get_point(f1_it->get_p1_index()), facets1.get_point(f1_it->get_p2_index()), facets1.get_point(f1_it->get_p3_index()));
                // if point is on a facet, return true
                bool pt_on_side(false);
                if (f1.contains_point(f2_ip, pt_on_side, precision))
                {
#ifdef DEBUG_INTERSECT_MESHES_2D_FACET_SORTER
                    cout << "Intersect_Meshes_2D::Facet_Sorter::sort sorting facets2 facet p1: " << 
                            f2_it->get_p1_index() << " p2: " << f2_it->get_p2_index() << " p3: " << 
                            f2_it->get_p3_index() << " internal point is on surface of f1 p1: " << 
                            f1_it->get_p1_index() << " p2: " << f1_it->get_p2_index() << " p3: " << 
                            f1_it->get_p3_index() << "\n";
#endif
                    // both f1 and f2 are on the surface of each other
                    f2_inside_f1.push_back(*f2_it);
                    break; // go to next f1
                }
            }
        }
#ifdef DEBUG_INTERSECT_MESHES_2D_FACET_SORTER
        cout << "Intersect_Meshes_2D::Facet_Sorter::sort end\n";
#endif
    }
    
    void Intersect_Meshes_2D::Facet_Sorter::clear() 
    {
        f1_inside_f2.clear();
        f2_inside_f1.clear();
    }
    
    Intersect_Meshes_2D::Intersect_Meshes_2D() {}
    
    void Intersect_Meshes_2D::intersect_facets(Facets& facets1, Facets& facets2, 
            const Point_2D::Measurement precision)
    {
#ifdef DEBUG_INTERSECT_MESHES_2D
        cout << "intersect_meshes_2D::intersect_facets begin\n";
#endif
        // create a list of facet builders for facets1
        vector<Facet_Builder> f1_builders;
        for (Facets::const_iterator it = facets1.begin(); it != facets1.end(); ++it)
        {
            f1_builders.push_back(Facet_Builder(true, 
                    Facet_2D(facets1.get_point(it->get_p1_index()), 
                            facets1.get_point(it->get_p2_index()), 
                            facets1.get_point(it->get_p3_index())), 
                    precision));
        }
        
        I_Pt_Locator i_pt_locator(precision);
        
        // intersect all facets together
        Facets::const_iterator f2_it = facets2.begin();
        while (f2_it != facets2.end())
        {
            Facet_Builder f2_builder(false, 
                    Facet_2D(facets2.get_point(f2_it->get_p1_index()), 
                            facets2.get_point(f2_it->get_p2_index()),
                            facets2.get_point(f2_it->get_p3_index())), 
                    precision);
            
            // intersect all of facets1 with facets2 facet
            for (vector<Facet_Builder>::iterator f1_it = f1_builders.begin(); f1_it != f1_builders.end(); ++f1_it)
            {
#ifdef DEBUG_INTERSECT_MESHES_2D
                cout << "intersect_meshes_2D::intersect_facets facet1 p1 x: " << f1_it->get_facet().get_point1()->get_x() << 
                        " y: " << f1_it->get_facet().get_point1()->get_y() << " p2 x: " << f1_it->get_facet().get_point2()->get_x() << 
                        " y: " << f1_it->get_facet().get_point2()->get_y() << " p3 x: " << f1_it->get_facet().get_point3()->get_x() << 
                        " y: " << f1_it->get_facet().get_point3()->get_y() << "\n";
                cout << "intersect_meshes_2D::intersect_facets facet2 p1 x: " << f2_builder.get_facet().get_point1()->get_x() << " y: " << 
                        f2_builder.get_facet().get_point1()->get_y() << " p2 x: " << f2_builder.get_facet().get_point2()->get_x() << " y: " << 
                        f2_builder.get_facet().get_point2()->get_y() << " p3 x: " << f2_builder.get_facet().get_point3()->get_x() << " y: " << 
                        f2_builder.get_facet().get_point3()->get_y() << "\n";
                cout.flush();
#endif
                I_Pt_List intersect_pts;
                if (i_pt_locator(f1_it->get_facet(), f2_builder.get_facet(), intersect_pts))
                {
#ifdef DEBUG_INTERSECT_MESHES_2D
                    cout << "intersect_meshes_2D::intersect_facets intersect_pts size: " << intersect_pts.size() << "\n";
                    cout.flush();
                    for (I_Pt_List::const_iterator ip_it = intersect_pts.begin(); ip_it != intersect_pts.end(); ++ip_it)
                    {
                        cout << "intersect_meshes_2D::intersect_facets intersect point x: " << ip_it->pt->get_x() << " y: " << 
                                ip_it->pt->get_y() << " f1_loc: " << ip_it->f1_loc << " f2_loc: " << ip_it->f2_loc << "\n";
                        cout.flush();
                    }
                    cout << "intersect_meshes_2D::intersect_facets adding intersection to facets1 facet builder\n";
                    cout.flush();
#endif                    
                    // add intersection to f1 facet
                    f1_it->add_intersection(intersect_pts);
#ifdef DEBUG_INTERSECT_MESHES_2D
                    cout << "intersect_meshes_2D::intersect_facets adding intersection to facets2 facet builder\n";
                    cout.flush();
#endif                    
                    // add intersection to f2_facet
                    f2_builder.add_intersection(intersect_pts);
                }
            }
            
#ifdef DEBUG_INTERSECT_MESHES_2D
            cout << "intersect_meshes_2D::intersect_facets forming new facets for facets2 facet\n";
            cout.flush();
#endif                    
            // build facets2 facet
            Facets new_facets;
            if (f2_builder.form_new_facets(new_facets))
            {
#ifdef DEBUG_INTERSECT_MESHES_2D
                cout << "intersect_meshes_2D::intersect_facets formed new facets for facets2 facet\n";
                cout.flush();
                for (Facets::const_iterator it = new_facets.begin(); it != new_facets.end(); ++it)
                {
                    shared_ptr<Point_2D> p1(new_facets.get_point(it->get_p1_index()));
                    shared_ptr<Point_2D> p2(new_facets.get_point(it->get_p2_index()));
                    shared_ptr<Point_2D> p3(new_facets.get_point(it->get_p3_index()));
                    cout << "        intersect_meshes_2D::intersect_facets new facet p1 x: " << p1->get_x() <<
                            " y: " << p1->get_y() << " p2 x: " << p2->get_x() << " y: " << p2->get_y() << 
                            " p3 x: " << p3->get_x() << " y: " << p3->get_y() << "\n";
                    cout.flush();
                }
#endif                    
                f2_it = facets2.replace_facet(*f2_it, new_facets);
                advance(f2_it, new_facets.size());
            }
            else
            {
#ifdef DEBUG_INTERSECT_MESHES_2D
                cout << "intersect_meshes_2D::intersect_facets no intersection\n";
                cout.flush();
#endif
                ++f2_it;
            }
        }
        
        // form facets1 facets
        Facets::const_iterator f1_it = facets1.begin();
        for (vector<Facet_Builder>::iterator fb1_it = f1_builders.begin(); fb1_it != f1_builders.end(); ++fb1_it)
        {
#ifdef DEBUG_INTERSECT_MESHES_2D
            cout << "intersect_meshes_2D::intersect_facets forming new facets for facets1 facet p1 x: " << 
                    fb1_it->get_facet().get_point1()->get_x() << " y: " << fb1_it->get_facet().get_point1()->get_y() << 
                    " p2 x: " << fb1_it->get_facet().get_point2()->get_x() << " y: " << 
                    fb1_it->get_facet().get_point2()->get_y() << " p3 x: " << 
                    fb1_it->get_facet().get_point3()->get_x() << " y: " << 
                    fb1_it->get_facet().get_point3()->get_y() << "\n";
            cout.flush();
#endif                    
            // build facets1 facet
            Facets new_facets;
            if (fb1_it->form_new_facets(new_facets))
            {
#ifdef DEBUG_INTERSECT_MESHES_2D
                cout << "intersect_meshes_2D::intersect_facets formed new facets for facets1 facet\n";
                cout.flush();
                for (Facets::const_iterator it = new_facets.begin(); it != new_facets.end(); ++it)
                {
                    shared_ptr<Point_2D> p1(new_facets.get_point(it->get_p1_index()));
                    shared_ptr<Point_2D> p2(new_facets.get_point(it->get_p2_index()));
                    shared_ptr<Point_2D> p3(new_facets.get_point(it->get_p3_index()));
                    cout << "        intersect_meshes_2D::intersect_facets new facet p1 x: " << p1->get_x() <<
                            " y: " << p1->get_y() << " p2 x: " << p2->get_x() <<
                            " y: " << p2->get_y() << " p3 x: " << p3->get_x() << 
                            " y: " << p3->get_y() << "\n";
                    cout.flush();
                }
#endif                    
                f1_it = facets1.replace_facet(*f1_it, new_facets);
                advance(f1_it, new_facets.size());
            }
            else
            {
#ifdef DEBUG_INTERSECT_MESHES_2D
                cout << "intersect_meshes_2D::intersect_facets no intersection\n";
                cout.flush();
#endif
                ++f1_it;
            }
        }
    }
    
    const bool Intersect_Meshes_2D::operator()(const Mesh_2D& mesh1, const Mesh_2D& mesh2, Mesh_2D& mesh1_result, Mesh_2D& mesh2_result)
    {
#ifdef INTERSECT_MESHES_2D
        cout << "Intersect_Meshes_2D::operator() begin\n";
#endif
        Facets facets1(mesh1);
        Facets facets2(mesh2);
        
#ifdef INTERSECT_MESHES_2D
        cout << "Intersect_Meshes_2D::operator() calling intersect_meshes\n";
#endif
        this->intersect_facets(facets1, facets2, mesh2_result.get_precision());
#ifdef INTERSECT_MESHES_2D
        cout << "Intersect_Meshes_2D::operator() after intersect_meshes\n";
#endif
        if (facets1.size() > mesh1.size() || facets2.size() > mesh2.size())
        {
#ifdef INTERSECT_MESHES_2D
            cout << "Intersect_Meshes_2D::operator() new facets were generated.  Updating mesh1_result\n";
#endif
            mesh1_result.clear();
            for (Facets::const_iterator it = facets1.begin(); it != facets1.end(); ++it)
            {
                mesh1_result.push_back(Facet_2D(facets1.get_point(it->get_p1_index()), facets1.get_point(it->get_p2_index()), facets1.get_point(it->get_p3_index())));
            }
#ifdef INTERSECT_MESHES_2D
            cout << "Intersect_Meshes_2D::operator() Updating mesh2_result\n";
#endif
            mesh2_result.clear();
            for (Facets::const_iterator it = facets2.begin(); it != facets2.end(); ++it)
            {
                mesh2_result.push_back(Facet_2D(facets2.get_point(it->get_p1_index()), facets2.get_point(it->get_p2_index()), facets2.get_point(it->get_p3_index())));
            }
#ifdef INTERSECT_MESHES_2D
            cout << "Intersect_Meshes_2D::operator() returning true\n";
#endif
            return true;
        }
        else
        {
#ifdef INTERSECT_MESHES_2D
            cout << "Intersect_Meshes_2D::operator() no new facets were generated.  Returning false\n";
#endif
            return false;
        }
    }
    
    void Intersect_Meshes_2D::difference(const Mesh_2D& mesh1, const Mesh_2D& mesh2, Mesh_2D& result)
    {
#ifdef DEBUG_INTERSECT_MESHES_2D_DIFFERENCE
        cout << "Intersect_Meshes_2D::difference begin\n";
#endif
        Facets facets1(mesh1);
        Facets facets2(mesh2);
        
        this->intersect_facets(facets2, facets1, result.get_precision());
#ifdef DEBUG_INTERSECT_MESHES_2D_DIFFERENCE
        cout << "Intersect_Meshes_2D::difference intersected meshes size facets1: " << facets1.size() << " facets2: " << facets2.size() << "\n";
        cout << "Intersect_Meshes_2D::difference facets1: polyhedron(points=[";
        for (Facets::pt_const_iterator it = facets1.pts_cbegin(); it != facets1.pts_cend(); ++it)
        {
            if (it != facets1.pts_cbegin())
                cout << ',';
            cout << '[' << (*it)->get_x() << ',' << (*it)->get_y() << ']';
        }
        cout << "], faces=[";
        for (Facets::const_iterator it = facets1.cbegin(); it != facets1.cend(); ++it)
        {
            if (it != facets1.cbegin())
                cout << ',';
            // write points in clockwise order because openscad prefers it
            cout << '[' << it->get_p1_index() << ',' << it->get_p3_index() << ',' << it->get_p2_index() << ']';
        }
        cout << "], convexity=4);\n";
        cout << "Intersect_Meshes_2D::difference facets2: polyhedron(points=[";
        for (Facets::pt_const_iterator it = facets2.pts_cbegin(); it != facets2.pts_cend(); ++it)
        {
            if (it != facets2.pts_cbegin())
                cout << ',';
            cout << '[' << (*it)->get_x() << ',' << (*it)->get_y() << ']';
        }
        cout << "], faces=[";
        for (Facets::const_iterator it = facets2.cbegin(); it != facets2.cend(); ++it)
        {
            if (it != facets2.cbegin())
                cout << ',';
            // write points in clockwise order because openscad prefers it
            cout << '[' << it->get_p1_index() << ',' << it->get_p3_index() << ',' << it->get_p2_index() << ']';
        }
        cout << "], convexity=4);\n";
        cout << "Intersect_Meshes_2D::difference sorting facets\n";
        cout.flush();
#endif
        Facet_Sorter facet_sorter(result.get_precision());
        facet_sorter.sort(facets1, facets2);
        result.clear();
#ifdef DEBUG_INTERSECT_MESHES_2D_DIFFERENCE
        cout << "Intersect_Meshes_2D::difference adding facets from facets1 that are outside of facets2\n";
#endif
        // add facets from facets1 that are not in facets2
        for (Facets::const_iterator it = facets1.begin(); it != facets1.end(); ++it)
        {
            if (facet_sorter.f1_inside_end() == find(facet_sorter.f1_inside_begin(), facet_sorter.f1_inside_end(), *it))
                result.push_back(Facet_2D(facets1.get_point(it->get_p1_index()), facets1.get_point(it->get_p2_index()), facets1.get_point(it->get_p3_index())));
        }
#ifdef DEBUG_INTERSECT_MESHES_2D_DIFFERENCE
            cout << "Intersect_Meshes_2D::difference end\n";
#endif
    }
    
    void Intersect_Meshes_2D::intersection(const Mesh_2D& mesh1, const Mesh_2D& mesh2, Mesh_2D& result)
    {
#ifdef DEBUG_INTERSECT_MESHES_2D_INTERSECTION
        cout << "Intersect_Meshes_2D::intersection begin\n";
#endif
        Facets facets1(mesh1);
        Facets facets2(mesh2);
        
        this->intersect_facets(facets2, facets1, result.get_precision());
#ifdef DEBUG_INTERSECT_MESHES_2D_INTERSECTION
        cout << "Intersect_Meshes_2D::intersection intersected meshes size facets1: " << facets1.size() << " facets2: " << facets2.size() << "\n";
        cout << "Intersect_Meshes_2D::intersection facets1: polyhedron(points=[";
        for (Facets::pt_const_iterator it = facets1.pts_cbegin(); it != facets1.pts_cend(); ++it)
        {
            if (it != facets1.pts_cbegin())
                cout << ',';
            cout << '[' << (*it)->get_x() << ',' << (*it)->get_y() << ']';
        }
        cout << "], faces=[";
        for (Facets::const_iterator it = facets1.cbegin(); it != facets1.cend(); ++it)
        {
            if (it != facets1.cbegin())
                cout << ',';
            // write points in clockwise order because openscad prefers it
            cout << '[' << it->get_p1_index() << ',' << it->get_p3_index() << ',' << it->get_p2_index() << ']';
        }
        cout << "], convexity=4);\n";
        cout << "Intersect_Meshes_2D::intersection facets2: polyhedron(points=[";
        for (Facets::pt_const_iterator it = facets2.pts_cbegin(); it != facets2.pts_cend(); ++it)
        {
            if (it != facets2.pts_cbegin())
                cout << ',';
            cout << '[' << (*it)->get_x() << ',' << (*it)->get_y() << ']';
        }
        cout << "], faces=[";
        for (Facets::const_iterator it = facets2.cbegin(); it != facets2.cend(); ++it)
        {
            if (it != facets2.cbegin())
                cout << ',';
            // write points in clockwise order because openscad prefers it
            cout << '[' << it->get_p1_index() << ',' << it->get_p3_index() << ',' << it->get_p2_index() << ']';
        }
        cout << "], convexity=4);\n";
        cout << "Intersect_Meshes_2D::intersection sorting facets\n";
        cout.flush();
#endif
        Facet_Sorter facet_sorter(result.get_precision());
        facet_sorter.sort(facets1, facets2);
        result.clear();
#ifdef DEBUG_INTERSECT_MESHES_2D_INTERSECTION
        cout << "Intersect_Meshes_2D::intersection adding facets from facets1 that are inside of facets2\n";
#endif
        for (Facet_Sorter::const_iterator it = facet_sorter.f1_inside_begin(); it != facet_sorter.f1_inside_end(); ++it)
        {
            result.push_back(Facet_2D(facets1.get_point(it->get_p1_index()), facets1.get_point(it->get_p2_index()), facets1.get_point(it->get_p3_index())));
        }
#ifdef DEBUG_INTERSECT_MESHES_2D_INTERSECTION
        cout << "Intersect_Meshes_2D::intersection end\n";
#endif
    }
    
    void Intersect_Meshes_2D::merge(const Mesh_2D& mesh1, const Mesh_2D& mesh2, Mesh_2D& result)
    {
#ifdef DEBUG_INTERSECT_MESHES_2D_MERGE
        cout << "Intersect_Meshes_2D::merge begin\n";
#endif
        Facets facets1(mesh1);
        Facets facets2(mesh2);
        
        this->intersect_facets(facets2, facets1, mesh1.get_precision());
#ifdef DEBUG_INTERSECT_MESHES_2D_MERGE
        cout << "Intersect_Meshes_2D::merge intersected meshes size facets1: " << facets1.size() << " facets2: " << facets2.size() << "\n";
        cout << "Intersect_Meshes_2D::merge facets1: polyhedron(points=[";
        for (Facets::pt_const_iterator it = facets1.pts_cbegin(); it != facets1.pts_cend(); ++it)
        {
            if (it != facets1.pts_cbegin())
                cout << ',';
            cout << '[' << (*it)->get_x() << ',' << (*it)->get_y() << ']';
        }
        cout << "], faces=[";
        for (Facets::const_iterator it = facets1.cbegin(); it != facets1.cend(); ++it)
        {
            if (it != facets1.cbegin())
                cout << ',';
            // write points in clockwise order because openscad prefers it
            cout << '[' << it->get_p1_index() << ',' << it->get_p3_index() << ',' << it->get_p2_index() << ']';
        }
        cout << "], convexity=4);\n";
        cout << "Intersect_Meshes_2D::merge facets2: polyhedron(points=[";
        for (Facets::pt_const_iterator it = facets2.pts_cbegin(); it != facets2.pts_cend(); ++it)
        {
            if (it != facets2.pts_cbegin())
                cout << ',';
            cout << '[' << (*it)->get_x() << ',' << (*it)->get_y() << ']';
        }
        cout << "], faces=[";
        for (Facets::const_iterator it = facets2.cbegin(); it != facets2.cend(); ++it)
        {
            if (it != facets2.cbegin())
                cout << ',';
            // write points in clockwise order because openscad prefers it
            cout << '[' << it->get_p1_index() << ',' << it->get_p3_index() << ',' << it->get_p2_index() << ']';
        }
        cout << "], convexity=4);\n";
        cout << "Intersect_Meshes_2D::merge sorting facets\n";
        cout.flush();
#endif
        Facet_Sorter facet_sorter(result.get_precision());
        facet_sorter.sort(facets1, facets2);
        result.clear();
#ifdef DEBUG_INTERSECT_MESHES_2D_MERGE
        cout << "Intersect_Meshes_2D::merge adding all facets from facets1\n";
#endif
        // add this mesh facets that are not inside mesh
        for (Facets::const_iterator it = facets1.begin(); it != facets1.end(); ++it)
        {
            result.push_back(Facet_2D(facets1.get_point(it->get_p1_index()), facets1.get_point(it->get_p2_index()), facets1.get_point(it->get_p3_index())));
        }
            
#ifdef DEBUG_INTERSECT_MESHES_2D_MERGE
        cout << "Intersect_Meshes_2D::merge adding facets from facets2 that are outside facets1\n";
#endif
        // add mesh facet that is not inside or on this mesh
        for (Facets::const_iterator it = facets2.begin(); it != facets2.end(); ++it)
        {
            if (facet_sorter.f2_inside_end() == find(facet_sorter.f2_inside_begin(), facet_sorter.f2_inside_end(), *it))
                result.push_back(Facet_2D(facets2.get_point(it->get_p1_index()), facets2.get_point(it->get_p2_index()), facets2.get_point(it->get_p3_index())));
        }
#ifdef DEBUG_INTERSECT_MESHES_2D_MERGE
        cout << "Intersect_Meshes_2D::merge end\n";
#endif
    }

}
