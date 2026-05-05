/*
  Copyright 2008-2011 The University of Texas at Austin

        Authors: Joe Rivera <transfix@ices.utexas.edu>
        Advisor: Chandrajit Bajaj <bajaj@cs.utexas.edu>

  This file is part of VolumeRover.

  VolumeRover is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License version 2.1 as published by the Free Software Foundation.

  VolumeRover is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/

#ifndef __CVCRAW_GEOMETRY_IO_H__
#define __CVCRAW_GEOMETRY_IO_H__

/*
 * CVC Raw Geometry file I/O functions:
 *  read, write
 *
 * Change Log:
 * 04/02/2010 - Moved this code into it's own header from cvcraw_geometry.h
 */

#include <stdexcept>
#include <fstream>
#include <iostream>
#include <vector>
#include <string>

#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lambda/lambda.hpp>

namespace 
{
  static inline bool ends_with(const std::string& haystack, const std::string& needle)
  {
    return haystack.rfind(needle) == haystack.size() - needle.size();
  }
}

namespace cvcraw_geometry
{
  template <class geometry_map_container>
  static inline void write(const geometry_map_container& geometry,
			   const std::string& filename)
  {
    using namespace std;
    using namespace boost; using boost::format;

    typedef typename geometry_map_container::points_t    points_t;
    typedef typename geometry_map_container::boundary_t  boundary_t;
    typedef typename geometry_map_container::normals_t   normals_t;
    typedef typename geometry_map_container::colors_t    colors_t;
    typedef typename geometry_map_container::lines_t     lines_t;
    typedef typename geometry_map_container::triangles_t triangles_t;
    typedef typename geometry_map_container::quads_t     quads_t;

    const points_t& points = geometry.points;
    const boundary_t& boundary = geometry.boundary;
    const normals_t& normals = geometry.normals;
    const colors_t& colors = geometry.colors;
    const lines_t& lines = geometry.lines;
    const triangles_t& tris = geometry.tris;
    const quads_t& quads = geometry.quads;
        
    ofstream outf(filename.c_str());
    if(!outf)
      throw runtime_error(string("Could not open ") + filename);
    
    unsigned int num_elems;
    if(boundary.size() != points.size()) //if we don't have boundary information
      {
	num_elems = 
	  lines.size() > 0 ? lines.size() :
	  tris.size() > 0 ? tris.size()   :
	  quads.size() > 0 ? quads.size() :
	  0;
      }
    else
      {
	num_elems =
	  tris.size() > 0 ? tris.size()/4   :
	  quads.size() > 0 ? quads.size()/6 :
	  0;
      }
    outf << points.size() << " " << num_elems << endl;
    if(!outf)
      throw runtime_error("Could not write number of points or number of tris to file!");

    // arand: 4-21-2011
    // changed to check file extension
    // normals and colors are only printed if .rawn/.rawnc/.rawc
    // have been specified...
    bool haveNormals = (normals.size() == points.size());
    bool printNormals = (ends_with(filename,".rawn") || ends_with(filename,".rawnc"));
    bool haveColors = (colors.size() == points.size());
    bool printColors = (ends_with(filename,".rawc") || ends_with(filename,".rawnc"));

    if (printNormals && !haveNormals) {
      std::cout << "WARNING: file with normals requested but not available." << std::endl;
    }

    if (printColors && !haveColors) {
      std::cout << "WARNING: file with normals requested but not available." << std::endl;
    }
    
    for(typename points_t::const_iterator i = points.begin();
	i != points.end();
	i++)
      {
	outf << (*i)[0] << " " << (*i)[1] << " " << (*i)[2];

	if(haveNormals && printNormals)
	  outf << " " 
	       << normals[std::distance(points.begin(),i)][0] << " " 
	       << normals[std::distance(points.begin(),i)][1] << " " 
	       << normals[std::distance(points.begin(),i)][2];
	if(haveColors && printColors)
	  outf << " " 
	       << colors[std::distance(points.begin(),i)][0] << " " 
	       << colors[std::distance(points.begin(),i)][1] << " " 
	       << colors[std::distance(points.begin(),i)][2];
	if(boundary.size() == points.size())
	  outf << " " << boundary[std::distance(points.begin(),i)];
	outf << endl;
	if(!outf)
	  throw runtime_error(str(format("Error writing vertex %1%") % std::distance(points.begin(),i)));
      }
    
    if(lines.size() != 0)
      {
	for(typename lines_t::const_iterator i = lines.begin(); i != lines.end(); i++)
	  {
	    typedef typename lines_t::value_type cell_type;
	    for(typename cell_type::const_iterator j = i->begin(); j != i->end(); j++)
	      {
		outf << *j;
		if(std::next(j) == i->end()) outf << endl;
		else outf << " ";
	      }

	    if(!outf)
	      throw runtime_error(str(format("Error writing line %1%") % std::distance(lines.begin(),i)));
	  }
      }
    else if(tris.size() != 0)
      {
	if(boundary.size() != points.size()) //if we don't have boundary info, don't treat these tris as part of a tet
	  {
	    for(typename triangles_t::const_iterator i = tris.begin(); i != tris.end(); i++)
	      {
		typedef typename triangles_t::value_type cell_type;
		for(typename cell_type::const_iterator j = i->begin(); j != i->end(); j++)
		  {
		    outf << *j;
		    if(std::next(j) == i->end()) outf << endl;
		    else outf << " ";
		  }
		
		if(!outf)
		  throw runtime_error(str(format("Error writing triangle %1%") % std::distance(tris.begin(),i)));
	      }
	  }
	else
	  {
	    for(unsigned int i = 0; i < tris.size()/4; i++)
	      {
		outf << tris[4*i][0] << " " << tris[4*i][1] << " " 
		     << tris[4*i][2] << " " << tris[4*i+1][2] << endl;
		if(!outf)
		  throw runtime_error(str(format("Error writing tetrahedron %1%") % i));
	      }
	  }
      }
    else if(quads.size() != 0)
      {
	if(boundary.size() != points.size()) //if we don't have boundary info, don't tread these quads as part of a hexa
	  {
	    for(typename quads_t::const_iterator i = quads.begin(); i != quads.end(); i++)
	      {
		typedef typename quads_t::value_type cell_type;
		for(typename cell_type::const_iterator j = i->begin(); j != i->end(); j++)
		  {
		    outf << *j;
		    if(std::next(j) == i->end()) outf << endl;
		    else outf << " ";
		  }
		
		if(!outf)
		  throw runtime_error(str(format("Error writing quad %1%") % std::distance(quads.begin(),i)));
	      }
	  }
	else
	  {
	    for(unsigned int i = 0; i < quads.size()/6; i++)
	      {
#if 0
		outf << quads[6*i][0] << " " << quads[6*i][1] << " "
		     << quads[6*i][2] << " " << quads[6*i][3] << " "
		     << quads[6*i+1][1] << " " << quads[6*i+1][0] << " "
		     << quads[6*i+1][3] << " " << quads[6*i+1][2] << endl;
#endif
		outf << quads[6*i][0] << " " << quads[6*i][1] << " "
		     << quads[6*i][2] << " " << quads[6*i][3] << " "
		     << quads[6*i+1][0] << " " << quads[6*i+1][1] << " "
		     << quads[6*i+1][2] << " " << quads[6*i+1][3] << endl;
		

		if(!outf)
		  throw runtime_error(str(format("Error writing hexahedron %1%") % i));		  
	      }
	  }
      }
  }

  template <class geometry_map_container>
  static inline geometry_map_container& read(geometry_map_container& geometry,
					     const std::string& filename)
  {
    using namespace std;
    using namespace boost; using boost::format;
    using namespace boost::lambda;

    typedef typename geometry_map_container::points_t    points_t;
    typedef typename geometry_map_container::boundary_t  boundary_t;
    typedef typename geometry_map_container::normals_t   normals_t;
    typedef typename geometry_map_container::colors_t    colors_t;
    typedef typename geometry_map_container::lines_t     lines_t;
    typedef typename geometry_map_container::triangles_t triangles_t;
    typedef typename geometry_map_container::quads_t     quads_t;

    typedef typename points_t::value_type    point_t;
    typedef typename normals_t::value_type   normal_t;
    typedef typename colors_t::value_type    color_t;
    typedef typename lines_t::value_type     line_t;
    typedef typename triangles_t::value_type triangle_t;
    typedef typename quads_t::value_type     quad_t;

    points_t points;
    boundary_t boundary;
    normals_t normals;
    colors_t colors;
    lines_t lines;
    triangles_t tris;
    quads_t quads;
        
    ifstream inf(filename.c_str());
    if(!inf)
      throw runtime_error(string("Could not open ") + filename);
    
    unsigned int line_num = 0;
    string line;
    vector<string> split_line;
    unsigned int num_verts, num_elems;
    
    getline(inf, line); line_num++;
    if(!inf)
      throw runtime_error(str(format("Error reading file %1%, line %2%")
				   % filename
				   % line_num));
    trim(line);
    split(split_line,
	  line,
	  is_any_of(" "),
	  token_compress_on);
    if(split_line.size() != 1 &&
       split_line.size() != 2)
      throw runtime_error(str(format("Not a cvc-raw file (wrong number of tokens: [%1%])")
			      % split_line.size()));

    try
      {
	num_verts = lexical_cast<unsigned int>(split_line[0]);
	num_elems = split_line.size() > 1 ?
	  lexical_cast<unsigned int>(split_line[1]) : 0;
    
	for(unsigned int vt = 0; vt < num_verts; vt++)
	  {
	    getline(inf, line); line_num++;
	    if(!inf)
	      throw runtime_error(str(format("Error reading file %1%, line %2%")
				      % filename
				      % line_num));
	    trim(line);
	    split(split_line,
		  line,
		  is_any_of(" "),
		  token_compress_on);
	    switch(split_line.size())
	      {
	      case 3: // raw
	      case 4: // raw with boundary
		{
		  point_t point;
		  for(int i = 0; i < 3; i++)
		    point[i] = lexical_cast<double>(split_line[i]);
		  points.push_back(point);
		  if(split_line.size() == 4)
		    boundary.push_back(lexical_cast<int>(split_line[3]));
		}
		break;
	      case 6: // rawn or rawc
	      case 7: // rawn or rawc with boundary
		{
		  bool is_rawn = ends_with(filename,".rawn");
		  point_t point;
		  for(int i = 0; i < 3; i++)
		    point[i] = lexical_cast<double>(split_line[i]);
		  points.push_back(point);

		  if(is_rawn)
		    {
		      normal_t normal;
		      for(int i = 3; i < 6; i++)
			normal[i-3] = lexical_cast<double>(split_line[i]);
		      normals.push_back(normal);
		    }
		  else
		    {
		      color_t color;
		      for(int i = 3; i < 6; i++)
			color[i-3] = lexical_cast<double>(split_line[i]);
		      colors.push_back(color);
		    }

		  if(split_line.size() == 7)
		    boundary.push_back(lexical_cast<int>(split_line[6]));
		}
		break;
	      case 9:  // rawnc
	      case 10: // rawnc with boundary
		{
		  point_t point;
		  normal_t normal;
		  color_t color;

		  for(int i = 0; i < 3; i++)
		    point[i-0] = lexical_cast<double>(split_line[i]);
		  points.push_back(point);

		  for(int i = 3; i < 6; i++)
		    normal[i-3] = lexical_cast<double>(split_line[i]);
		  normals.push_back(normal);

		  for(int i = 6; i < 9; i++)
		    color[i-6] = lexical_cast<double>(split_line[i]);
		  colors.push_back(color);

		  if(split_line.size() == 10)
		    boundary.push_back(lexical_cast<int>(split_line[9]));
		}
		break;
	      default:
		{
		  throw runtime_error(str(format("Not a cvc-raw file (wrong number of tokens: [%1%])")
					  % split_line.size()));
		}
		break;
	      }
	  }
    
	for(unsigned int tri = 0; tri < num_elems; tri++)
	  {
	    getline(inf, line); line_num++;
	    if(!inf)
	      throw runtime_error(str(format("Error reading file %1%, line %2%")
				      % filename
				      % line_num));
	    trim(line);
	    split(split_line,
		  line,
		  is_any_of(" "),
		  token_compress_on);
	    switch(split_line.size())
	      {
	      case 2: // lines
		{
		  line_t line;
		  for(unsigned int i = 0; i < split_line.size(); i++)
		    line[i] = lexical_cast<unsigned int>(split_line[i]);
		  lines.push_back(line);
		}
		break;
	      case 3: // tris
		{
		  triangle_t triangle;
		  for(unsigned int i = 0; i < split_line.size(); i++)
		    triangle[i] = lexical_cast<unsigned int>(split_line[i]);
		  tris.push_back(triangle);
		}
		break;
	      case 4: // quads or tetrahedrons
		{
		  //if we didn't collect boundary information earlier, then we have quads
		  if(boundary.size() != points.size())
		    {
		      quad_t quad;
		      for(unsigned int i = 0; i < split_line.size(); i++)
			quad[i] = lexical_cast<unsigned int>(split_line[i]);
		      quads.push_back(quad);
		    }
		  else
		    {
		      unsigned int t[4];
		      for(unsigned int i = 0; i < 4; i++)
			t[i] = lexical_cast<unsigned int>(split_line[i]);
		      
		      triangle_t tet_tris[4];
		      tet_tris[0][0] = t[0]; tet_tris[0][1] = t[2]; tet_tris[0][2] = t[1];
		      tet_tris[1][0] = t[1]; tet_tris[1][1] = t[2]; tet_tris[1][2] = t[3];
		      tet_tris[2][0] = t[0]; tet_tris[2][1] = t[3]; tet_tris[2][2] = t[2];
		      tet_tris[3][0] = t[0]; tet_tris[3][1] = t[1]; tet_tris[3][2] = t[3];

		      for(unsigned int i = 0; i < 4; i++)
			tris.push_back(tet_tris[i]);
		    }
		}
		break;
	      case 8: // hexahedrons
		{
		  //if we don't have boundary information, then something is wrong!
		  if(boundary.size() != points.size())
		    throw runtime_error("Incorrect cvc-raw file: missing boundary info for hex verts");

		  unsigned int t[8];
		  for(unsigned int i = 0; i < 8; i++)
		    t[i] = lexical_cast<unsigned int>(split_line[i]);
		  // TODO: investigate hex and tet index ordering... somehow I had to do this different than geoframe
		  quad_t hex_quads[6];
#if 0
		  hex_quads[0][0] = t[0]; hex_quads[0][1] = t[3]; hex_quads[0][2] = t[2]; hex_quads[0][3] = t[1];
		  hex_quads[1][0] = t[4]; hex_quads[1][1] = t[5]; hex_quads[1][2] = t[6]; hex_quads[1][3] = t[7];
		  hex_quads[2][0] = t[0]; hex_quads[2][1] = t[4]; hex_quads[2][2] = t[7]; hex_quads[2][3] = t[3];
		  hex_quads[3][0] = t[1]; hex_quads[3][1] = t[2]; hex_quads[3][2] = t[6]; hex_quads[3][3] = t[5];
		  hex_quads[4][0] = t[0]; hex_quads[4][1] = t[1]; hex_quads[4][2] = t[5]; hex_quads[4][3] = t[4];
		  hex_quads[5][0] = t[2]; hex_quads[5][1] = t[3]; hex_quads[5][2] = t[7]; hex_quads[5][3] = t[6];
#endif

		  hex_quads[0][0] = t[0]; hex_quads[0][1] = t[1]; hex_quads[0][2] = t[2]; hex_quads[0][3] = t[3];
		  hex_quads[1][0] = t[4]; hex_quads[1][1] = t[5]; hex_quads[1][2] = t[6]; hex_quads[1][3] = t[7];
		  hex_quads[2][0] = t[0]; hex_quads[2][1] = t[4]; hex_quads[2][2] = t[7]; hex_quads[2][3] = t[3];
		  hex_quads[3][0] = t[1]; hex_quads[3][1] = t[2]; hex_quads[3][2] = t[6]; hex_quads[3][3] = t[5];
		  hex_quads[4][0] = t[0]; hex_quads[4][1] = t[1]; hex_quads[4][2] = t[5]; hex_quads[4][3] = t[4];
		  hex_quads[5][0] = t[2]; hex_quads[5][1] = t[3]; hex_quads[5][2] = t[7]; hex_quads[5][3] = t[6];

		  for(unsigned int i = 0; i < 6; i++)
		    quads.push_back(hex_quads[i]);
		}
		break;
	      default:
		{
		  throw runtime_error(str(format("Not a cvc-raw file {num components found: %1%}") % split_line.size()));   
		}
		break;
	      }
	  }
      }
    catch(std::exception& e)
      {
	throw runtime_error(str(format("Error reading file %1%, line %2%, contents: '%3%', reason: %4%")
				% filename
				% line_num
				% line
				% string(e.what())));
      }

    //adjust indices if they start from 1 rather than 0.  (i.e. search for a 0 in each index
    //list. If there are no zeros, decrement each index by one)
#ifdef CVCRAW_GEOMETRY__CORRECT_INDEX_START
    {
      bool lines_has_zero = false;
      bool tris_has_zero = false;
      bool quads_has_zero = false;

      for(typename lines_t::const_iterator i = lines.begin();
	  i != lines.end();
	  i++)
	if(find(i->begin(), i->end(), 0) != i->end())
	  {
	    lines_has_zero = true;
	    break;
	  }

      for(typename triangles_t::const_iterator i = tris.begin();
	  i != tris.end();
	  i++)
	if(find(i->begin(), i->end(), 0) != i->end())
	  {
	    tris_has_zero = true;
	    break;
	  }
      
      for(typename quads_t::const_iterator i = quads.begin();
	  i != quads.end();
	  i++)
	if(find(i->begin(), i->end(), 0) != i->end())
	  {
	    quads_has_zero = true;
	    break;
	  }

      if(!lines_has_zero)
	for(typename lines_t::iterator i = lines.begin();
	    i != lines.end();
	    i++)
	  for_each(i->begin(), i->end(), --_1);

      if(!tris_has_zero)
	for(typename triangles_t::iterator i = tris.begin();
	    i != tris.end();
	    i++)
	  for_each(i->begin(), i->end(), --_1);

      if(!quads_has_zero)
	for(typename quads_t::iterator i = quads.begin();
	    i != quads.end();
	    i++)
	  for_each(i->begin(), i->end(), --_1);
    }
#endif

    geometry.points = points;
    geometry.boundary = boundary;
    geometry.normals = normals;
    geometry.colors = colors;
    geometry.lines = lines;
    geometry.tris = tris;
    geometry.quads = quads;
    
    return geometry;
  }

}

#endif
