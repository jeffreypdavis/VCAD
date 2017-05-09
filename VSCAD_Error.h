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
 * File:   VSCAD_Error.h
 * Author: Jeffrey Davis
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

