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
 * File:   stl.h
 * Author: Jeffrey Davis
 */

#ifndef STL_H
#define STL_H

#include <string>
#include <stdexcept>
using namespace std;

namespace VCAD_lib
{
    class STL_Error : public runtime_error {
    public:
        explicit STL_Error (const string& what_arg);
        explicit STL_Error (const char* what_arg);
    };
    
    class Mesh_3D;
    
    const int read_stl(Mesh_3D& mesh, const string& filename, string& comment,
            const bool clockwise_order=false, const bool ignore_unv=false);
    
    const int read_stl_ascii(Mesh_3D& mesh, const string& filename, string& comment,
            const bool clockwise_order=false, const bool ignore_unv=false);
    
    /*
     * Test if architecture is little endian.  returns true if little endian, 
     * false if big endian
     */
    const bool is_little_endian();

    const short read_short(ifstream& ifs);
    
    const short read_short_cbo(ifstream& ifs);

    void write_short(ofstream& ofs, const short val);

    void write_short_cbo(ofstream& ofs, const short val);
    
    const unsigned int read_uint(ifstream& ifs);
    
    const unsigned int read_uint_cbo(ifstream& ifs);
    
    void write_uint(ofstream& ofs, const unsigned int val);
    
    void write_uint_cbo(ofstream& ofs, const unsigned int val);
    
    const float read_float(ifstream& ifs);
    
    const float read_float_cbo(ifstream& ifs);
    
    void write_float(ofstream& ofs, const float val);
    
    void write_float_cbo(ofstream& ofs, const float val);
    
    const int read_stl_bin(Mesh_3D& mesh, const string& filename, string& comment,
            const bool clockwise_order=false, const bool ignore_unv=false);
    
    /*
     * read stl file using normal byte order (little or big endian depending
     * on machine architecture)
     */
    const int read_stl_bin_nbo(Mesh_3D& mesh, const string& filename, string& comment,
            const bool clockwise_order=false, const bool ignore_unv=false);
    
    /*
     * read stl file and change byte order (little to big endian or big to little
     * endian)
     */
    const int read_stl_bin_cbo(Mesh_3D& mesh, const string& filename, string& comment,
            const bool clockwise_order=false, const bool ignore_unv=false);
    
    // write stl as an ascii stl file
    void write_stl(const Mesh_3D& mesh, const string& filename, 
            const string& comment, const bool clockwise_order=false, 
            const bool zero_unv=false);
    
    // write stl as a binary stl file
    void write_stl_bin(const Mesh_3D& mesh, const string& filename, 
            const string& comment, const short attribute=0, 
            const bool clockwise_order=false, const bool zero_unv=false);
    
    // write stl as a binary stl file using normal machine byte order
    void write_stl_bin_nbo(const Mesh_3D& mesh, const string& filename, 
            const string& comment, const short attribute=0, 
            const bool clockwise_order=false, const bool zero_unv=false);
    
    // write stl as a binary stl file by changing byte order
    void write_stl_bin_cbo(const Mesh_3D& mesh, const string& filename, 
            const string& comment, const short attribute=0, 
            const bool clockwise_order=false, const bool zero_unv=false);
}

#endif /* STL_H */

