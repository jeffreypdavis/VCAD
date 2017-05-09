/* 
 * File:   VSCAD_Error.cpp
 * Author: Jeffrey Davis
 * Version: 1.0
 * 
 * Created on May 3, 2016, 9:44 PM
 */

#include "VSCAD_Error.h"

namespace VCAD_lib
{
    VSCAD_Error::VSCAD_Error (const string& what_arg) : runtime_error(what_arg) {}
    VSCAD_Error::VSCAD_Error (const char* what_arg) : runtime_error(what_arg) {}
}
