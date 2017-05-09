/* 
 * File:   VSCAD_Error.h
 * Author: Jeffrey Davis
 *
 * Created on May 3, 2016, 9:44 PM
 */

#ifndef VSCAD_ERROR_H
#define VSCAD_ERROR_H

#include <stdexcept>
using namespace std;

namespace VCAD_lib
{

    class VSCAD_Error : public runtime_error {
    public:
        explicit VSCAD_Error (const string& what_arg);
        explicit VSCAD_Error (const char* what_arg);
    };

}

#endif /* VSCAD_ERROR_H */

