/*
  Copyright 2008 The University of Texas at Austin
  
	Authors: Jose Rivera <transfix@ices.utexas.edu>
	Advisor: Chandrajit Bajaj <bajaj@cs.utexas.edu>

  This file is part of Volume Rover.

  Volume Rover is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  Volume Rover is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with iotree; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

/* $Id: cvcraw_geometry.h,v 1.2 2008/10/05 17:45:00 transfix Exp $ */

/*
  Joe R. - transfix@ices.utexas.edu - 10/08/2008 - v0.1

  This single header library provides some useful facilities for CVC geometry
  I/O.  It provides a geometry type called cvcraw_geometry::geomery_t to collect
  geometry information.  All information is contained in STL conforming containers,
  except the boundary information which uses boost::dynamic_bitset<> (I don't think
  dynamic_bitset<> can be used like an STL container exactly).

  Currently, the geometry_t class supplies 3 operations:

  calculate_extents() - simply calculates the smallest bounding box such that
                        all vertices fit inside that box.
  merge()             - merges *this with another geometry_t object
  tri_surface()       - extracts the boundary of *this and returns it as a
                        triangulated surface.  This is useful for code that
			only works with tri surfaces.

  Note: I am considering moving these operations outside of geometry_t and instead
  implementing them as template functions so they aren't dependent on the
  cvcraw_geometry::geometry_t type exactly.  However I need to spend some more time
  thinking about what an appropriate set of traits would be for generic geometry types.
                  
  The library also provides the following I/O template functions:
  read()              - reads the geometry from the filename provided and stores it
                        into the container object provided.  It figures out what kind
			of geometry is present in the file based on it's structure.
			It supports all CVC mesh types including those used by LBIE.
  write()             - writes the geometry to the filename provided.  It will write
                        the correct kind of geometry file depending on what information
			is present in the geometry container.

  If there is an error during reading or writing, an exception will be thrown.
*/

#ifndef __CVCRAW_GEOMETRY_H__
#define __CVCRAW_GEOMETRY_H__

#include <fstream>
#include <vector>
#include <string>
#include <stdexcept>
#include <algorithm>
#include <boost/format.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/array.hpp>
#include <boost/utility.hpp>
#include <boost/dynamic_bitset.hpp>

namespace 
{
  static inline bool ends_with(const std::string& haystack, const std::string& needle)
  {
    return haystack.rfind(needle) == haystack.size() - needle.size();
  }
}

namespace cvcraw_geometry
{
  class geometry_t
  {
  public:
    typedef boost::array<double,3> point_t;
    typedef boost::array<double,3> normal_t;
    typedef boost::array<double,3> color_t;
    typedef boost::array<int,2>    line_t;
    typedef boost::array<int,3>    triangle_t;
    typedef boost::array<int,4>    quad_t;

    typedef std::vector<point_t>    points_t;
    typedef boost::dynamic_bitset<> boundary_t;
    typedef std::vector<normal_t>   normals_t;
    typedef std::vector<color_t>    colors_t;
    typedef std::vector<line_t>     lines_t;
    typedef std::vector<triangle_t> triangles_t;
    typedef std::vector<quad_t>     quads_t;

    geometry_t() {}
    geometry_t(const geometry_t& copy)
      : points(copy.points), normals(copy.normals),
      colors(copy.colors), lines(copy.lines),
      tris(copy.tris), quads(copy.quads),
      min_ext(copy.min_ext), max_ext(copy.max_ext)
      {}
    ~geometry_t() {}

    geometry_t& operator=(const geometry_t& copy)
      {
	points = copy.points;
	normals = copy.normals;
	colors = copy.colors;
	lines = copy.lines;
	tris = copy.tris;
	quads = copy.quads;
	min_ext = copy.min_ext;
	max_ext = copy.max_ext;
	return *this;
      }

    points_t    points;
    boundary_t  boundary;
    normals_t   normals;
    colors_t    colors;
    lines_t     lines;
    triangles_t tris;
    quads_t     quads;

    point_t min_ext;
    point_t max_ext;

    point_t min() const { return min_ext; }
    point_t max() const { return max_ext; }

    void calculate_extents()
    {
      using namespace std;

      if(points.empty())
	{
	  fill(min_ext.begin(),min_ext.end(),0.0);
	  fill(max_ext.begin(),max_ext.end(),0.0);
	  return;
	}

      min_ext = max_ext = points[0];
      for(points_t::const_iterator i = points.begin();
	  i != points.end();
	  i++)
	{
	  if(min_ext[0] > (*i)[0]) min_ext[0] = (*i)[0];
	  if(min_ext[1] > (*i)[1]) min_ext[1] = (*i)[1];
	  if(min_ext[2] > (*i)[2]) min_ext[2] = (*i)[2];
	  if(max_ext[0] < (*i)[0]) max_ext[0] = (*i)[0];
	  if(max_ext[1] < (*i)[1]) max_ext[1] = (*i)[1];
	  if(max_ext[2] < (*i)[2]) max_ext[2] = (*i)[2];
	}
    }

    geometry_t& merge(const geometry_t& geom)
      {
	using namespace std;

	//append vertex info
	points.insert(points.end(),
		      geom.points.begin(),
		      geom.points.end());
	normals.insert(normals.end(),
		       geom.normals.begin(),
		       geom.normals.end());
	colors.insert(colors.end(),
		      geom.colors.begin(),
		      geom.colors.end());
	
	//we apparently don't have normal iterators for boundary_t :/
	unsigned int old_size = boundary.size();
	boundary.resize(boundary.size() + geom.boundary.size());
	for(unsigned int i = 0; i < geom.boundary.size(); i++)
	  boundary[old_size + i] = geom.boundary[i];

	//append and modify index info
	lines.insert(lines.end(),
		     geom.lines.begin(),
		     geom.lines.end());
	{
	  lines_t::iterator i = lines.begin();
	  advance(i, lines.size() - geom.lines.size());
	  while(i != lines.end())
	    {
	      for(line_t::iterator j = i->begin();
		  j != i->end();
		  j++)
		*j += points.size() - geom.points.size();
	      i++;
	    }
	}

	tris.insert(tris.end(),
		    geom.tris.begin(),
		    geom.tris.end());
	{
	  triangles_t::iterator i = tris.begin();
	  advance(i, tris.size() - geom.tris.size());
	  while(i != tris.end())
	    {
	      for(triangle_t::iterator j = i->begin();
		  j != i->end();
		  j++)
		*j += points.size() - geom.points.size();
	      i++;
	    }
	}

	quads.insert(quads.end(),
		     geom.quads.begin(),
		     geom.quads.end());
	{
	  quads_t::iterator i = quads.begin();
	  advance(i, quads.size() - geom.quads.size());
	  while(i != quads.end())
	    {
	      for(quad_t::iterator j = i->begin();
		  j != i->end();
		  j++)
		*j += points.size() - geom.points.size();
	      i++;
	    }
	}

	calculate_extents();
	return *this;
      }

    //returns a simple tri surface for the boundary
    //doesn't remove extra non boundary points
    geometry_t tri_surface() const
    {
      geometry_t new_geom;

      for(points_t::const_iterator i = points.begin();
	  i != points.end();
	  i++)
	new_geom.points.push_back(*i);

      //if we don't have boundary info, just insert the tris directly
      if(boundary.size() != points.size())
	{
	  for(triangles_t::const_iterator i = tris.begin();
	      i != tris.end();
	      i++)
	    new_geom.tris.push_back(*i);

	  for(quads_t::const_iterator i = quads.begin();
	      i != quads.end();
	      i++)
	    {
	      triangle_t quad_tris[2];
	      quad_tris[0][0] = (*i)[0];
	      quad_tris[0][1] = (*i)[1];
	      quad_tris[0][2] = (*i)[3];

	      quad_tris[1][0] = (*i)[1];
	      quad_tris[1][1] = (*i)[2];
	      quad_tris[1][2] = (*i)[3];

	      new_geom.tris.push_back(quad_tris[0]);
	      new_geom.tris.push_back(quad_tris[1]);
	    }
	}
      else
	{
	  for(triangles_t::const_iterator i = tris.begin();
	      i != tris.end();
	      i++)
	    if(boundary[(*i)[0]] &&
	       boundary[(*i)[1]] &&
	       boundary[(*i)[2]])
	      new_geom.tris.push_back(*i);

	  for(quads_t::const_iterator i = quads.begin();
	      i != quads.end();
	      i++)
	    if(boundary[(*i)[0]] &&
	       boundary[(*i)[1]] &&
	       boundary[(*i)[2]] &&
	       boundary[(*i)[3]])
	      {
		triangle_t quad_tris[2];
		quad_tris[0][0] = (*i)[0];
		quad_tris[0][1] = (*i)[1];
		quad_tris[0][2] = (*i)[3];
		
		quad_tris[1][0] = (*i)[1];
		quad_tris[1][1] = (*i)[2];
		quad_tris[1][2] = (*i)[3];
		
		new_geom.tris.push_back(quad_tris[0]);
		new_geom.tris.push_back(quad_tris[1]);
	      }
	}
      
      return new_geom;
    }
  };
  

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
    
    for(typename points_t::const_iterator i = points.begin();
	i != points.end();
	i++)
      {
	outf << (*i)[0] << " " << (*i)[1] << " " << (*i)[2];
	if(normals.size() == points.size())
	  outf << " " 
	       << normals[distance(points.begin(),i)][0] << " " 
	       << normals[distance(points.begin(),i)][1] << " " 
	       << normals[distance(points.begin(),i)][2];
	if(colors.size() == points.size())
	  outf << " " 
	       << colors[distance(points.begin(),i)][0] << " " 
	       << colors[distance(points.begin(),i)][1] << " " 
	       << colors[distance(points.begin(),i)][2];
	if(boundary.size() == points.size())
	  outf << " " << boundary[distance(points.begin(),i)];
	outf << endl;
	if(!outf)
	  throw runtime_error(str(format("Error writing vertex %1%") % distance(points.begin(),i)));
      }
    
    if(lines.size() != 0)
      {
	for(typename lines_t::const_iterator i = lines.begin(); i != lines.end(); i++)
	  {
	    typedef typename lines_t::value_type cell_type;
	    for(typename cell_type::const_iterator j = i->begin(); j != i->end(); j++)
	      {
		outf << *j;
		if(next(j) == i->end()) outf << endl;
		else outf << " ";
	      }

	    if(!outf)
	      throw runtime_error(str(format("Error writing line %1%") % distance(lines.begin(),i)));
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
		    if(next(j) == i->end()) outf << endl;
		    else outf << " ";
		  }
		
		if(!outf)
		  throw runtime_error(str(format("Error writing triangle %1%") % distance(tris.begin(),i)));
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
		    if(next(j) == i->end()) outf << endl;
		    else outf << " ";
		  }
		
		if(!outf)
		  throw runtime_error(str(format("Error writing quad %1%") % distance(quads.begin(),i)));
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
    split(split_line,
		 line,
		 is_any_of(" "));
    if(split_line.size() != 1 &&
       split_line.size() != 2)
      throw runtime_error("Not a cvc-raw file");

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
	    split(split_line,
		  line,
		  is_any_of(" "));
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
		  throw runtime_error("Not a cvc-raw file");
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
	    split(split_line,
		  line,
		  is_any_of(" "));
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
#pragma message("TODO: investigate hex and tet index ordering... somehow I had to do this different than geoframe")
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
	throw runtime_error(str(format("Error reading file %1%, line %2%, reason: %3%")
				% filename
				% line_num
				% string(e.what())));
      }

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
