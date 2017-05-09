/* 
 * File:   Facet.h
 * Author: Jeffrey Davis
 *
 * Created on March 23, 2017, 10:47 AM
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

