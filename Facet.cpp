/* 
 * File:   Facet.cpp
 * Author: Jeffrey Davis
 * 
 * Created on March 23, 2017, 10:47 AM
 */

#include "Facet.h"

namespace VCAD_lib
{

    Facet::Facet(const int point1_index, const int point2_index, const int point3_index) : 
            p1_index(point1_index), p2_index(point2_index), p3_index(point3_index) {}

    void Facet::invert_unv()
    {
        // swap
        int temp = p2_index;
        p2_index = p3_index;
        p3_index = temp;
    }
    
    const bool Facet::operator==(const Facet& facet) const
    {
        return (facet.get_p1_index() == p1_index && facet.get_p2_index() == p2_index && facet.get_p3_index() == p3_index) || 
                (facet.get_p2_index() == p1_index && facet.get_p3_index() == p2_index && facet.get_p1_index() == p3_index) || 
                (facet.get_p3_index() == p1_index && facet.get_p1_index() == p2_index && facet.get_p2_index() == p3_index);
    }
    
}
