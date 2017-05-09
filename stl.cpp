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
 * File:   stl.cpp
 * Author: Jeffrey Davis
 */

#include "stl.h"
#include <fstream>
#include <cstring>
#include <memory>
#include <sstream>
#include "Point_3D.h"
#include "Vector_3D.h"
#include "Facet_3D.h"
#include "Mesh_3D.h"

namespace VCAD_lib
{
    STL_Error::STL_Error (const string& what_arg) : runtime_error(what_arg) {}
    STL_Error::STL_Error (const char* what_arg) : runtime_error(what_arg) {}

    // taken from http://www.cplusplus.com/faq/sequences/strings/trim/
    std::string& trim_right_inplace(std::string& s, const std::string& delimiters = " \f\n\r\t\v" )
    {
      return s.erase( s.find_last_not_of( delimiters ) + 1 );
    }

    // taken from http://www.cplusplus.com/faq/sequences/strings/trim/
    std::string& trim_left_inplace(std::string& s, const std::string& delimiters = " \f\n\r\t\v" )
    {
      return s.erase( 0, s.find_first_not_of( delimiters ) );
    }

    // taken from http://www.cplusplus.com/faq/sequences/strings/trim/
    std::string& trim_inplace(std::string& s, const std::string& delimiters = " \f\n\r\t\v" )
    {
      return trim_left_inplace( trim_right_inplace( s, delimiters ), delimiters );
    }

    std::string& to_lower(std::string& s)
    {
        for (int i = 0; i < s.size(); ++i)
            s[i] = tolower(s[i]);
        return s;
    }
    
    const int read_stl(Mesh_3D& mesh, const string& filename, string& comment,
            const bool clockwise_order, const bool ignore_unv)
    {
        ifstream ifs;
        ifs.open(filename);
        if (!ifs.is_open())
            throw STL_Error("Unable to open file '" + filename + "'");
        
        char test[7] = {0,0,0,0,0,0,0};
        ifs.get(test,6);
        if (ifs.eof() || strlen(test) < 5)
            throw STL_Error("Invalid STL file: " + filename);
        else // successful
        {
            ifs.close();
            if (tolower(test[0]) == 's' && tolower(test[1]) == 'o' && 
                    tolower(test[2]) == 'l' & tolower(test[3]) == 'i' &&
                    tolower(test[4]) == 'd')
                // ascii stl file
                return read_stl_ascii(mesh, filename, comment, clockwise_order, ignore_unv);
            else
                // binary stl file
                return read_stl_bin(mesh, filename, comment, clockwise_order, ignore_unv);
        }
    }
    
    const int read_stl_ascii(Mesh_3D& mesh, const string& filename, string& comment,
            const bool clockwise_order, const bool ignore_unv)
    {
        ifstream ifs;
        ifs.open(filename);
        if (!ifs.is_open())
            throw STL_Error("Unable to open file '" + filename + "'");
        string line;
        getline(ifs, line);
        // test for errors
        if (ifs.eof() || ifs.fail())
            throw STL_Error("Invalid STL file: " + filename);

        if (line.size() < 5)
            throw STL_Error("Invalid STL file: " + filename);

        // test if line starts with solid
        string solid = line.substr(0,5);
        if (to_lower(solid).compare("solid") != 0)
            throw STL_Error("Invalid ascii STL file: " + filename);

        // set comment
        comment = line.substr(5,line.size());
        trim_inplace(comment);

        Mesh_3D temp(mesh.get_precision()); // temporary mesh to store facets in case of error

        // read facets
        while (ifs >> line)
        {
            if (to_lower(line).compare("facet") != 0)
            {
                if (line.compare("endsolid") == 0)
                    break; // end of stl file
                // handle error
                stringstream ss;
                ss << "expected 'facet', found '" << line << "' on facet " << temp.size();
                STL_Error e(string(ss.str()));
                throw e;
            }
            if (!(ifs >> line) || to_lower(line).compare("normal") != 0)
            {
                // handle error
                stringstream ss;
                ss << "expected 'normal', found '" << line << "' on facet " << temp.size();
                STL_Error e(string(ss.str()));
                throw e;
            }
            double unv_x(0);
            double unv_y(0);
            double unv_z(0);
            if (!(ifs >> unv_x))
            {
                // handle error
                stringstream ss;
                ss << "unable to read the first facet normal decimal value on facet " << temp.size();
                STL_Error e(string(ss.str()));
                throw e;
            }
            if (!(ifs >> unv_y))
            {
                // handle error
                stringstream ss;
                ss << "unable to read the second facet normal decimal value on facet " << temp.size();
                STL_Error e(string(ss.str()));
                throw e;
            }
            if (!(ifs >> unv_z))
            {
                // handle error
                stringstream ss;
                ss << "unable to read the third facet normal decimal value on facet " << temp.size();
                STL_Error e(string(ss.str()));
                throw e;
            }
            if (!(ifs >> line) || to_lower(line).compare("outer") != 0)
            {
                // handle error
                stringstream ss;
                ss << "expected 'outer', found '" << line << "' on facet " << temp.size();
                STL_Error e(string(ss.str()));
                throw e;
            }
            if (!(ifs >> line) || to_lower(line).compare("loop") != 0)
            {
                // handle error
                stringstream ss;
                ss << "expected 'loop', found '" << line << "' on facet " << temp.size();
                STL_Error e(string(ss.str()));
                throw e;
            }
            if (!(ifs >> line) || to_lower(line).compare("vertex") != 0)
            {
                // handle error
                stringstream ss;
                ss << "expected 'vertex', found '" << line << "' on facet " << temp.size();
                STL_Error e(string(ss.str()));
                throw e;
            }
            double x;
            double y;
            double z;
            if (!(ifs >> x))
            {
                // handle error
                stringstream ss;
                ss << "unable to read the first decimal value on facet " << temp.size() << " vertex #1";
                STL_Error e(string(ss.str()));
                throw e;
            }
            if (!(ifs >> y))
            {
                // handle error
                stringstream ss;
                ss << "unable to read the second decimal value on facet " << temp.size() << " vertex #1";
                STL_Error e(string(ss.str()));
                throw e;
            }
            if (!(ifs >> z))
            {
                // handle error
                stringstream ss;
                ss << "unable to read the third decimal value on facet " << temp.size() << " vertex #1";
                STL_Error e(string(ss.str()));
                throw e;
            }
            shared_ptr<Point_3D> p1(new Point_3D(x,y,z));
            
            if (!(ifs >> line) || to_lower(line).compare("vertex") != 0)
            {
                // handle error
                stringstream ss;
                ss << "expected 'vertex', found '" << line << "' on facet " << temp.size();
                STL_Error e(string(ss.str()));
                throw e;
            }
            if (!(ifs >> x))
            {
                // handle error
                stringstream ss;
                ss << "unable to read the first decimal value on facet " << temp.size() << " vertex #2";
                STL_Error e(string(ss.str()));
                throw e;
            }
            if (!(ifs >> y))
            {
                // handle error
                stringstream ss;
                ss << "unable to read the second decimal value on facet " << temp.size() << " vertex #2";
                STL_Error e(string(ss.str()));
                throw e;
            }
            if (!(ifs >> z))
            {
                // handle error
                stringstream ss;
                ss << "unable to read the third decimal value on facet " << temp.size() << " vertex #2";
                STL_Error e(string(ss.str()));
                throw e;
            }
            shared_ptr<Point_3D> p2(new Point_3D(x,y,z));

            if (!(ifs >> line) || to_lower(line).compare("vertex") != 0)
            {
                // handle error
                stringstream ss;
                ss << "expected 'vertex', found '" << line << "' on facet " << temp.size();
                STL_Error e(string(ss.str()));
                throw e;
            }
            if (!(ifs >> x))
            {
                // handle error
                stringstream ss;
                ss << "unable to read the first decimal value on facet " << temp.size() << " vertex #3";
                STL_Error e(string(ss.str()));
                throw e;
            }
            if (!(ifs >> y))
            {
                // handle error
                stringstream ss;
                ss << "unable to read the second decimal value on facet " << temp.size() << " vertex #3";
                STL_Error e(string(ss.str()));
                throw e;
            }
            if (!(ifs >> z))
            {
                // handle error
                stringstream ss;
                ss << "unable to read the third decimal value on facet " << temp.size() << " vertex #3";
                STL_Error e(string(ss.str()));
                throw e;
            }
            shared_ptr<Point_3D> p3(new Point_3D(x,y,z));
            
            if (!(ifs >> line) || to_lower(line).compare("endloop") != 0)
            {
                // handle error
                stringstream ss;
                ss << "expected 'endloop', found '" << line << "' on facet " << temp.size();
                STL_Error e(string(ss.str()));
                throw e;
            }
            if (!(ifs >> line) || to_lower(line).compare("endfacet") != 0)
            {
                // handle error
                stringstream ss;
                ss << "expected 'endfacet', found '" << line << "' on facet " << temp.size();
                STL_Error e(string(ss.str()));
                throw e;
            }
            
            // form facet
            if (ignore_unv)
                temp.push_back(Facet_3D(p1, p2, p3, clockwise_order));
            else
                temp.push_back(Facet_3D(p1, p2, p3, Vector_3D(unv_x,unv_y,unv_z),clockwise_order));
        }
        
        ifs.close();
        
        mesh.clear();
        for (Mesh_3D::const_iterator it = temp.cbegin(); it != temp.cend(); ++it)
            mesh.push_back(*it);
        
        return temp.size();
    }
    
    const bool is_little_endian()
    {
        // this would be little endian 1
        stringstream ss;
        ss.put(0x01);
        ss.put(0x00);
        ss.put(0x00);
        ss.put(0x00);
        int temp = 0;
        ss.read((char *) &temp, 4);
        return temp == 1;
    }

    const short read_short(ifstream& ifs)
    {
        short result = 0;
        ifs.read((char *)&result, 2);
        if (ifs.eof() || ifs.fail())
            throw STL_Error("Invalid STL file.  Unable to read short");
        return result;
    }

    const short read_short_cbo(ifstream& ifs)
    {
        char c1 = ifs.get();
        if (ifs.eof())
            throw STL_Error("Invalid STL file.  Unable to read short");
        char c2 = ifs.get();
        if (ifs.eof())
            throw STL_Error("Invalid STL file.  Unable to read short");
        stringstream ss;
        ss.put(c2);
        ss.put(c1);
        short result = 0;
        ss.read((char *)&result, 2);
        if (ss.fail() || ss.eof())
            throw STL_Error("Error changing byte order for a short");
        return result;
    }

    void write_short(ofstream& ofs, const short val)
    {
        ofs.write((char *)&val, 2);
        if (ofs.fail())
            throw STL_Error("Error writing short to stl file");
    }

    void write_short_cbo(ofstream& ofs, const short val)
    {
        stringstream ss;
        ss.write((char *)&val, 2);
        char c1 = ss.get();
        if (ss.eof())
            throw STL_Error("Error changing byte order for a short");
        char c2 = ss.get();
        if (ss.eof())
            throw STL_Error("Error changing byte order for a short");
        ofs.put(c2);
        if (ofs.fail())
            throw STL_Error("Error writing short to binary stl file");
        ofs.put(c1);
        if (ofs.fail())
            throw STL_Error("Error writing short to binary stl file");
    }

    const unsigned int read_uint(ifstream& ifs)
    {
        unsigned int result = 0;
        ifs.read((char *)&result, 4);
        if (ifs.eof() || ifs.fail())
            throw STL_Error("Invalid STL file.  Unable to read unsigned int");
        return result;
    }
    
    const unsigned int read_uint_cbo(ifstream& ifs)
    {
        char c1 = ifs.get();
        if (ifs.eof())
            throw STL_Error("Invalid STL file.  Unable to read unsigned int");
        char c2 = ifs.get();
        if (ifs.eof())
            throw STL_Error("Invalid STL file.  Unable to read unsigned int");
        char c3 = ifs.get();
        if (ifs.eof())
            throw STL_Error("Invalid STL file.  Unable to read unsigned int");
        char c4 = ifs.get();
        if (ifs.eof())
            throw STL_Error("Invalid STL file.  Unable to read unsigned int");
        stringstream ss;
        ss.put(c4);
        ss.put(c3);
        ss.put(c2);
        ss.put(c1);
        unsigned int result = 0;
        ss.read((char *)&result, 4);
        if (ss.fail() || ss.eof())
            throw STL_Error("Error changing byte order for an unsigned int");
        return result;
    }

    void write_uint(ofstream& ofs, const unsigned int val)
    {
        ofs.write((char *)&val, 4);
        if (ofs.fail())
            throw STL_Error("Error writing unsigned int to stl file");
    }

    void write_uint_cbo(ofstream& ofs, const unsigned int val)
    {
        stringstream ss;
        ss.write((char *)&val, 4);
        char c1 = ss.get();
        if (ss.eof())
            throw STL_Error("Error changing byte order for an unsigned int");
        char c2 = ss.get();
        if (ss.eof())
            throw STL_Error("Error changing byte order for an unsigned int");
        char c3 = ss.get();
        if (ss.eof())
            throw STL_Error("Error changing byte order for an unsigned int");
        char c4 = ss.get();
        if (ss.eof())
            throw STL_Error("Error changing byte order for an unsigned int");
        ofs.put(c4);
        if (ofs.fail())
            throw STL_Error("Error writing unsigned int to binary stl file");
        ofs.put(c3);
        if (ofs.fail())
            throw STL_Error("Error writing unsigned int to binary stl file");
        ofs.put(c2);
        if (ofs.fail())
            throw STL_Error("Error writing unsigned int to binary stl file");
        ofs.put(c1);
        if (ofs.fail())
            throw STL_Error("Error writing unsigned int to binary stl file");
    }

    const float read_float(ifstream& ifs)
    {
        float result = 0;
        ifs.read((char *)&result, 4);
        if (ifs.eof() || ifs.fail())
            throw STL_Error("Invalid STL file.  Unable to read float");
        return result;
    }
    
    const float read_float_cbo(ifstream& ifs)
    {
        char c1 = ifs.get();
        if (ifs.eof())
            throw STL_Error("Invalid STL file.  Unable to read float");
        char c2 = ifs.get();
        if (ifs.eof())
            throw STL_Error("Invalid STL file.  Unable to read float");
        char c3 = ifs.get();
        if (ifs.eof())
            throw STL_Error("Invalid STL file.  Unable to read float");
        char c4 = ifs.get();
        if (ifs.eof())
            throw STL_Error("Invalid STL file.  Unable to read float");
        stringstream ss;
        ss.put(c4);
        ss.put(c3);
        ss.put(c2);
        ss.put(c1);
        float result = 0;
        ss.read((char *)&result, 4);
        if (ss.fail() || ss.eof())
            throw STL_Error("Error changing byte order for a float");
        return result;
    }

    void write_float(ofstream& ofs, const float val)
    {
        ofs.write((char *)&val, 4);
        if (ofs.fail())
            throw STL_Error("Error writing float to stl file");
    }

    void write_float_cbo(ofstream& ofs, const float val)
    {
        stringstream ss;
        ss.write((char *)&val, 4);
        char c1 = ss.get();
        if (ss.eof())
            throw STL_Error("Error changing byte order for a float");
        char c2 = ss.get();
        if (ss.eof())
            throw STL_Error("Error changing byte order for a float");
        char c3 = ss.get();
        if (ss.eof())
            throw STL_Error("Error changing byte order for a float");
        char c4 = ss.get();
        if (ss.eof())
            throw STL_Error("Error changing byte order for a float");
        ofs.put(c4);
        if (ofs.fail())
            throw STL_Error("Error writing float to binary stl file");
        ofs.put(c3);
        if (ofs.fail())
            throw STL_Error("Error writing float to binary stl file");
        ofs.put(c2);
        if (ofs.fail())
            throw STL_Error("Error writing float to binary stl file");
        ofs.put(c1);
        if (ofs.fail())
            throw STL_Error("Error writing float to binary stl file");
    }
    
    const int read_stl_bin(Mesh_3D& mesh, const string& filename, string& comment,
            const bool clockwise_order, const bool ignore_unv)
    {
        // binary stl files are usually little endian
        if (is_little_endian())
            return read_stl_bin_nbo(mesh, filename, comment, clockwise_order, ignore_unv);
        else
            return read_stl_bin_cbo(mesh, filename, comment, clockwise_order, ignore_unv);
    }
    
    const int read_stl_bin_nbo(Mesh_3D& mesh, const string& filename, string& comment,
            const bool clockwise_order, const bool ignore_unv)
    {
        ifstream ifs;
        ifs.open(filename, ios::binary);
        if (!ifs.is_open())
            throw STL_Error("Unable to open file '" + filename + "'");

        // read 80 bytes for the comment
        char buf [81];
        buf[80] = 0;
        ifs.read(buf,80);
        if (ifs.eof() || ifs.fail())
            throw STL_Error("Invalid STL file: " + filename);
        comment = string(buf);

        // read 4 byte int to find out how many facets there are
        unsigned int num_facets = read_uint(ifs);
        
        Mesh_3D temp(mesh.get_precision()); // temporary mesh to hold facets
        ifs.peek();
        while (!ifs.eof())
        {
            // read unit normal vector
            float x = read_float(ifs);
            float y = read_float(ifs);
            float z = read_float(ifs);
            Vector_3D unv(x, y, z);
            
            // read point 1
            x = read_float(ifs);
            y = read_float(ifs);
            z = read_float(ifs);
            shared_ptr<Point_3D> p1(new Point_3D(x,y,z));
            
            // read point 2
            x = read_float(ifs);
            y = read_float(ifs);
            z = read_float(ifs);
            shared_ptr<Point_3D> p2(new Point_3D(x,y,z));
            
            // read point 3
            x = read_float(ifs);
            y = read_float(ifs);
            z = read_float(ifs);
            shared_ptr<Point_3D> p3(new Point_3D(x,y,z));

            // read attribute
            short attribute = read_short(ifs);
            ifs.peek(); // for next iteration

            // form facet
            if (ignore_unv)
                temp.push_back(Facet_3D(p1, p2, p3, clockwise_order));
            else
                temp.push_back(Facet_3D(p1, p2, p3, unv, clockwise_order));
        }
        
        ifs.close();
        
        if (temp.size() != num_facets)
        {
            stringstream ss;
            ss << "Invalid STL file. Expected " << num_facets << " facets and found " << temp.size();
            STL_Error e(ss.str());
            throw e;
        }
        
        mesh.clear();
        for (Mesh_3D::const_iterator it = temp.cbegin(); it != temp.cend(); ++it)
            mesh.push_back(*it);
        return temp.size();
    }
    
    const int read_stl_bin_cbo(Mesh_3D& mesh, const string& filename, string& comment,
            const bool clockwise_order, const bool ignore_unv)
    {
        ifstream ifs;
        ifs.open(filename, ios::binary);
        if (!ifs.is_open())
            throw STL_Error("Unable to open file '" + filename + "'");

        // read 80 bytes for the comment
        char buf [81];
        buf[80] = 0;
        ifs.read(buf,80);
        if (ifs.eof() || ifs.fail())
            throw STL_Error("Invalid STL file: " + filename);
        comment = string(buf);

        // read 4 byte int to find out how many facets there are
        unsigned int num_facets = read_uint_cbo(ifs);

        Mesh_3D temp(mesh.get_precision()); // temporary mesh to hold facets
        ifs.peek();
        while (!ifs.eof())
        {
            // read unit normal vector
            float x = read_float_cbo(ifs);
            float y = read_float_cbo(ifs);
            float z = read_float_cbo(ifs);
            Vector_3D unv(x, y, z);
            
            // read point 1
            x = read_float_cbo(ifs);
            y = read_float_cbo(ifs);
            z = read_float_cbo(ifs);
            shared_ptr<Point_3D> p1(new Point_3D(x, y, z));

            // read point 2
            x = read_float_cbo(ifs);
            y = read_float_cbo(ifs);
            z = read_float_cbo(ifs);
            shared_ptr<Point_3D> p2(new Point_3D(x, y, z));

            // read point 3
            x = read_float_cbo(ifs);
            y = read_float_cbo(ifs);
            z = read_float_cbo(ifs);
            shared_ptr<Point_3D> p3(new Point_3D(x, y, z));

            // read attribute
            short attribute = read_short_cbo(ifs);
            ifs.peek(); // for next iteration

            // form facet
            if (ignore_unv)
                temp.push_back(Facet_3D(p1, p2, p3, clockwise_order));
            else
                temp.push_back(Facet_3D(p1, p2, p3, unv, clockwise_order));
        }
        
        ifs.close();
        
        if (temp.size() != num_facets)
        {
            stringstream ss;
            ss << "Invalid STL file. Expected " << num_facets << " facets and found " << temp.size();
            STL_Error e(ss.str());
            throw e;
        }

        mesh.clear();
        for (Mesh_3D::const_iterator it = temp.cbegin(); it != temp.cend(); ++it)
            mesh.push_back(*it);
        return temp.size();
    }
    
    void write_stl(const Mesh_3D& mesh, const string& filename, 
            const string& comment, const bool clockwise_order, 
            const bool zero_unv)
    {
        ofstream ofs;
        ofs.open(filename);
        if (!ofs.is_open())
        {
            STL_Error e("Unable to open file '" + filename + "'");
            throw e;
        }
        
        // erase any newlines
        string cmnt(comment);
        size_t pos = string::npos;
        while ((pos = cmnt.find('\r')) != string::npos)
            cmnt.erase(pos,1);
        while ((pos = cmnt.find('\n')) != string::npos)
            cmnt.erase(pos,1);
        
        // write comment
        if (cmnt.size() > 0)
        {
            ofs << "solid ";
            if (ofs.fail())
                throw STL_Error("Error writing to file: " + filename);
            ofs << cmnt;
            if (ofs.fail())
                throw STL_Error("Error writing to file: " + filename);
            ofs << endl;
            if (ofs.fail())
                throw STL_Error("Error writing to file: " + filename);
        }
        else
        {
            ofs << "solid";
            if (ofs.fail())
                throw STL_Error("Error writing to file: " + filename);
            ofs << endl;
            if (ofs.fail())
                throw STL_Error("Error writing to file: " + filename);
        }
        
        for (Mesh_3D::const_iterator it = mesh.begin(); it != mesh.end(); ++it)
        {
            // write unit normal
            ofs << "  facet normal ";
            if (ofs.fail())
                throw STL_Error("Error writing to file: " + filename);
            if (zero_unv)
            {
                ofs << "0 0 0";
                if (ofs.fail())
                    throw STL_Error("Error writing to file: " + filename);
                ofs << endl;
                if (ofs.fail())
                    throw STL_Error("Error writing to file: " + filename);
            }
            else
            {
                ofs << it->get_unv().get_x();
                if (ofs.fail())
                    throw STL_Error("Error writing to file: " + filename);
                ofs << ' ';
                if (ofs.fail())
                    throw STL_Error("Error writing to file: " + filename);
                ofs << it->get_unv().get_y();
                if (ofs.fail())
                    throw STL_Error("Error writing to file: " + filename);
                ofs << ' ';
                if (ofs.fail())
                    throw STL_Error("Error writing to file: " + filename);
                ofs << it->get_unv().get_z();
                if (ofs.fail())
                    throw STL_Error("Error writing to file: " + filename);
                ofs << endl;
                if (ofs.fail())
                    throw STL_Error("Error writing to file: " + filename);
            }
            ofs << "    outer loop";
            if (ofs.fail())
                throw STL_Error("Error writing to file: " + filename);
            ofs << endl;
            if (ofs.fail())
                throw STL_Error("Error writing to file: " + filename);
            ofs << "      vertex ";
            if (ofs.fail())
                throw STL_Error("Error writing to file: " + filename);
            ofs << it->get_point1()->get_x();
            if (ofs.fail())
                throw STL_Error("Error writing to file: " + filename);
            ofs << ' ';
            if (ofs.fail())
                throw STL_Error("Error writing to file: " + filename);
            ofs << it->get_point1()->get_y();
            if (ofs.fail())
                throw STL_Error("Error writing to file: " + filename);
            ofs << ' ';
            if (ofs.fail())
                throw STL_Error("Error writing to file: " + filename);
            ofs << it->get_point1()->get_z();
            if (ofs.fail())
                throw STL_Error("Error writing to file: " + filename);
            ofs << endl;
            if (ofs.fail())
                throw STL_Error("Error writing to file: " + filename);
            ofs << "      vertex ";
            if (ofs.fail())
                throw STL_Error("Error writing to file: " + filename);
            ofs << it->get_point2()->get_x();
            if (ofs.fail())
                throw STL_Error("Error writing to file: " + filename);
            ofs << ' ';
            if (ofs.fail())
                throw STL_Error("Error writing to file: " + filename);
            ofs << it->get_point2()->get_y();
            if (ofs.fail())
                throw STL_Error("Error writing to file: " + filename);
            ofs << ' '; 
            if (ofs.fail())
                throw STL_Error("Error writing to file: " + filename);
            ofs << it->get_point2()->get_z();
            if (ofs.fail())
                throw STL_Error("Error writing to file: " + filename);
            ofs << endl;
            if (ofs.fail())
                throw STL_Error("Error writing to file: " + filename);
            ofs << "      vertex ";
            if (ofs.fail())
                throw STL_Error("Error writing to file: " + filename);
            ofs << it->get_point3()->get_x();
            if (ofs.fail())
                throw STL_Error("Error writing to file: " + filename);
            ofs << ' ';
            if (ofs.fail())
                throw STL_Error("Error writing to file: " + filename);
            ofs << it->get_point3()->get_y();
            if (ofs.fail())
                throw STL_Error("Error writing to file: " + filename);
            ofs << ' ';
            if (ofs.fail())
                throw STL_Error("Error writing to file: " + filename);
            ofs << it->get_point3()->get_z();
            if (ofs.fail())
                throw STL_Error("Error writing to file: " + filename);
            ofs << endl;
            if (ofs.fail())
                throw STL_Error("Error writing to file: " + filename);
            ofs << "    endloop";
            if (ofs.fail())
                throw STL_Error("Error writing to file: " + filename);
            ofs << endl;
            if (ofs.fail())
                throw STL_Error("Error writing to file: " + filename);
            ofs << "  endfacet";
            if (ofs.fail())
                throw STL_Error("Error writing to file: " + filename);
            ofs << endl;
            if (ofs.fail())
                throw STL_Error("Error writing to file: " + filename);
        }
        if (cmnt.size() > 0)
        {
            ofs << "endsolid ";
            if (ofs.fail())
                throw STL_Error("Error writing to file: " + filename);
            ofs << cmnt;
            if (ofs.fail())
                throw STL_Error("Error writing to file: " + filename);
            ofs << endl;
            if (ofs.fail())
                throw STL_Error("Error writing to file: " + filename);
        }
        else
        {
            ofs << "endsolid";
            if (ofs.fail())
                throw STL_Error("Error writing to file: " + filename);
            ofs << endl;
            if (ofs.fail())
                throw STL_Error("Error writing to file: " + filename);
        }
        
        ofs.close();
    }
    
    void write_stl_bin(const Mesh_3D& mesh, const string& filename, 
            const string& comment, const short attribute, 
            const bool clockwise_order, const bool zero_unv)
    {
        // binary stl files are usually little endian
        // with counter clockwise point order and unv pointing outside
        if (is_little_endian())
            write_stl_bin_nbo(mesh, filename, comment, attribute, clockwise_order, zero_unv);
        else
            write_stl_bin_cbo(mesh, filename, comment, attribute, clockwise_order, zero_unv);
    }
    
    // write binary stl file with 'normal' byte order - normal depends on machine
    // architecture and is either little endian or big endian.  Little endian
    // is the typical expected format of stl files
    void write_stl_bin_nbo(const Mesh_3D& mesh, const string& filename, 
            const string& comment, const short attribute, 
            const bool clockwise_order, const bool zero_unv)
    {
        ofstream ofs;
        ofs.open(filename, ios::binary);
        if (!ofs.is_open())
            throw STL_Error("Unable to open file '" + filename + "'");
        
        // write comment
        const char * c_comment = comment.c_str();
        
        if (strlen(c_comment) >= 80)
        {
            ofs.write(c_comment, 80);
            if (ofs.fail())
                throw STL_Error("Unable to write file '" + filename + "'");
        }
        else
        {
            int length = strlen(c_comment);
            ofs.write(c_comment, length);
            if (ofs.fail())
                throw STL_Error("Unable to write file '" + filename + "'");
            length = 80 - length;
            for (int i = 0; i < length; ++i)
            {
                ofs.put(0);
                if (ofs.fail())
                    throw STL_Error("Unable to write file '" + filename + "'");
            }
        }
        
        // write number of facets
        write_uint(ofs, mesh.size());
        
        for (Mesh_3D::const_iterator it = mesh.begin(); it != mesh.end(); ++it)
        {
            // write unit normal vector
            if (zero_unv)
            {
                write_float(ofs, 0);
                write_float(ofs, 0);
                write_float(ofs, 0);
            }
            else
            {
                write_float(ofs, float(it->get_unv().get_x()));
                write_float(ofs, float(it->get_unv().get_y()));
                write_float(ofs, float(it->get_unv().get_z()));
            }
            // write point 1
            write_float(ofs, float(it->get_point1()->get_x()));
            write_float(ofs, float(it->get_point1()->get_y()));
            write_float(ofs, float(it->get_point1()->get_z()));
            
            if (clockwise_order) // swap points two and three
            {
                // write point 3
                write_float(ofs, float(it->get_point3()->get_x()));
                write_float(ofs, float(it->get_point3()->get_y()));
                write_float(ofs, float(it->get_point3()->get_z()));
                // write point 2
                write_float(ofs, float(it->get_point2()->get_x()));
                write_float(ofs, float(it->get_point2()->get_y()));
                write_float(ofs, float(it->get_point2()->get_z()));
            }
            else
            {
                // write point 2
                write_float(ofs, float(it->get_point2()->get_x()));
                write_float(ofs, float(it->get_point2()->get_y()));
                write_float(ofs, float(it->get_point2()->get_z()));
                // write point 3
                write_float(ofs, float(it->get_point3()->get_x()));
                write_float(ofs, float(it->get_point3()->get_y()));
                write_float(ofs, float(it->get_point3()->get_z()));
            }
            // write attribute
            write_short(ofs, attribute);
        }
        
        ofs.close();
    }
    
    void write_stl_bin_cbo(const Mesh_3D& mesh, const string& filename, 
            const string& comment, const short attribute, 
            const bool clockwise_order, const bool zero_unv)
    {
        ofstream ofs;
        ofs.open(filename, ios::binary);
        if (!ofs.is_open())
            throw STL_Error("Unable to open file '" + filename + "'");
        
        // write comment
        const char * c_comment = comment.c_str();
        
        if (strlen(c_comment) >= 80)
        {
            ofs.write(c_comment, 80);
            if (ofs.fail())
                throw STL_Error("Unable to write file '" + filename + "'");
        }
        else
        {
            int length = strlen(c_comment);
            ofs.write(c_comment, length);
            if (ofs.fail())
                throw STL_Error("Unable to write file '" + filename + "'");
            length = 80 - length;
            for (int i = 0; i < length; ++i)
            {
                ofs.put(0);
                if (ofs.fail())
                    throw STL_Error("Unable to write file '" + filename + "'");
            }
        }
        
        // write number of facets
        write_uint_cbo(ofs, mesh.size());
        
        for (Mesh_3D::const_iterator it = mesh.begin(); it != mesh.end(); ++it)
        {
            // write unit normal vector
            if (zero_unv)
            {
                write_float_cbo(ofs, 0);
                write_float_cbo(ofs, 0);
                write_float_cbo(ofs, 0);
            }
            else
            {
                write_float_cbo(ofs, float(it->get_unv().get_x()));
                write_float_cbo(ofs, float(it->get_unv().get_y()));
                write_float_cbo(ofs, float(it->get_unv().get_z()));
            }
            // write point 1
            write_float_cbo(ofs, float(it->get_point1()->get_x()));
            write_float_cbo(ofs, float(it->get_point1()->get_y()));
            write_float_cbo(ofs, float(it->get_point1()->get_z()));
            
            if (clockwise_order) // swap points two and three
            {
                // write point 3
                write_float_cbo(ofs, float(it->get_point3()->get_x()));
                write_float_cbo(ofs, float(it->get_point3()->get_y()));
                write_float_cbo(ofs, float(it->get_point3()->get_z()));
                // write point 2
                write_float_cbo(ofs, float(it->get_point2()->get_x()));
                write_float_cbo(ofs, float(it->get_point2()->get_y()));
                write_float_cbo(ofs, float(it->get_point2()->get_z()));
            }
            else
            {
                // write point 2
                write_float_cbo(ofs, float(it->get_point2()->get_x()));
                write_float_cbo(ofs, float(it->get_point2()->get_y()));
                write_float_cbo(ofs, float(it->get_point2()->get_z()));
                // write point 3
                write_float_cbo(ofs, float(it->get_point3()->get_x()));
                write_float_cbo(ofs, float(it->get_point3()->get_y()));
                write_float_cbo(ofs, float(it->get_point3()->get_z()));
            }
            // write attribute
            write_short_cbo(ofs, attribute);
        }
        
        ofs.close();
    }
}
