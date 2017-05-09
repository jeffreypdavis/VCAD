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
 * File:   Facet.cpp
 * Author: Jeffrey Davis
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
