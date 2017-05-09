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
 * File:   Facet.h
 * Author: Jeffrey Davis
 */

#ifndef FACET_H
#define FACET_H

namespace VCAD_lib
{

    class Facet {
    public:
        Facet(const int point1_index, const int point2_index, const int point3_index);
        const int get_p1_index() const { return p1_index; }
        const int get_p2_index() const { return p2_index; }
        const int get_p3_index() const { return p3_index; }
        void invert_unv();
        const bool operator==(const Facet& facet) const;
    private:
        int p1_index;
        int p2_index;
        int p3_index;
    };

}

#endif /* FACET_H */

