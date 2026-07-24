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

#include <cvcraw_geometry/cvcgeom.h>
#include <cvcraw_geometry/utility.h>
#ifndef DISABLE_CONVERSION
#include <cvcraw_geometry/cvcraw_geometry.h>
#endif

#ifdef CVCRAW_GEOMETRY_ENABLE_PROJECT
#include <cvcraw_geometry/project_verts.h>
#endif

#ifdef CVCRAW_GEOMETRY_ENABLE_BUNNY
#include <cvcraw_geometry/bunny.h>
#endif

#include <map>
#include <list>
#include <limits>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/regex.hpp>

namespace CVCGEOM_NAMESPACE
{
  cvcgeom_t::cvcgeom_t()
    : _extents_set(false)
  {
    for(uint64_t i = 0; i < 3; i++)
      _min[i] = _max[i] = 0.0;

    init_ptrs();
  }

  cvcgeom_t::cvcgeom_t(const cvcgeom_t& geom)
    : _points(geom._points), _boundary(geom._boundary),
      _normals(geom._normals), _colors(geom._colors),
      _lines(geom._lines), _triangles(geom._triangles),
      _quads(geom._quads), _extents_set(geom._extents_set),
      _min(geom._min), _max(geom._max) 
  {
    //make sure all our pointers are valid
    init_ptrs();
  }

  cvcgeom_t::~cvcgeom_t() {}

  void cvcgeom_t::copy(const cvcgeom_t &geom)
  {
    _points = geom._points;
    _boundary = geom._boundary;
    _normals = geom._normals;
    _colors = geom._colors;
    _lines = geom._lines;
    _triangles = geom._triangles;
    _quads = geom._quads;
    _extents_set = geom._extents_set;
    _min = geom._min;
    _max = geom._max;

    //make sure all our pointers are valid
    init_ptrs();
  }

  cvcgeom_t& cvcgeom_t::operator=(const cvcgeom_t& geom)
  {
    copy(geom); return *this;
  }

#ifndef DISABLE_CONVERSION
  cvcgeom_t::cvcgeom_t(const geometry_t& geom)
    : _extents_set(false)
  {
    for(uint64_t i = 0; i < 3; i++)
      _min[i] = _max[i] = 0.0;
    copy(geom);
  }

  void cvcgeom_t::copy(const geometry_t& geom)
  {
    init_ptrs();
    *_points = geom.points;
    *_boundary = geom.boundary;
    *_normals = geom.normals;
    *_colors = geom.colors;
    *_lines = geom.lines;
    *_triangles = geom.tris;
    *_quads = geom.quads;
  }

  cvcgeom_t& cvcgeom_t::operator=(const geometry_t& geom)
  {
    copy(geom); return *this;
  }

  cvcgeom_t::operator geometry_t() const
  {
    geometry_t geom;
    //init_ptrs();
    if(_points)
      geom.points = const_points();
    if(_boundary)
      geom.boundary = const_boundary();
    if(_normals)
      geom.normals = const_normals();
    if(_colors)
      geom.colors = const_colors();
    if(_lines)
      geom.lines = const_lines();
    if(_triangles)
      geom.tris = const_triangles();
    if(_quads)
      geom.quads = const_quads();
    geom.min_ext = min_point();
    geom.max_ext = max_point();
    return geom;
  }
#endif

  point_t cvcgeom_t::min_point() const
  {
    if(!_extents_set) calc_extents();
    return _min;
  }

  point_t cvcgeom_t::max_point() const
  {
    if(!_extents_set) calc_extents();
    return _max;
  }

  uint64_t cvcgeom_t::num_points() const
  {
    //TODO: throw an exception if _points.size() != _boundary.size() 
    //!= _normals.size() != _colors.size()
    return _points ? _points->size() : 0;
  }

  uint64_t cvcgeom_t::num_lines() const
  {
    return _lines ? _lines->size() : 0;
  }

  uint64_t cvcgeom_t::num_triangles() const
  {
    return _triangles ? _triangles->size() : 0;
  }

  uint64_t cvcgeom_t::num_quads() const
  {
    return _quads ? _quads->size() : 0;
  }

  bool cvcgeom_t::empty() const
  {
    if(num_points() == 0) // ||   (num_lines() == 0 && num_triangles() == 0 && num_quads() == 0))
      return true;
    return false;
  }

  cvcgeom_t& cvcgeom_t::merge(const cvcgeom_t& geom)
  {
    using namespace std;

    //append vertex info
    points().insert(points().end(),
                    geom.points().begin(),
                    geom.points().end());
    normals().insert(normals().end(),
                     geom.normals().begin(),
                     geom.normals().end());
    colors().insert(colors().end(),
                    geom.colors().begin(),
                    geom.colors().end());
	
    //we apparently don't have normal iterators for boundary_t :/
    unsigned int old_size = boundary().size();
    boundary().resize(boundary().size() + geom.boundary().size());
    for(unsigned int i = 0; i < geom.boundary().size(); i++)
      boundary()[old_size + i] = geom.boundary()[i];

    //append and modify index info
    lines().insert(lines().end(),
                   geom.lines().begin(),
                   geom.lines().end());
    {
      lines_t::iterator i = lines().begin();
      advance(i, lines().size() - geom.lines().size());
      while(i != lines().end())
	{
	  for(line_t::iterator j = i->begin();
	      j != i->end();
	      j++)
	    *j += points().size() - geom.points().size();
	  i++;
	}
    }
    
    triangles().insert(triangles().end(),
                       geom.triangles().begin(),
                       geom.triangles().end());
    {
      triangles_t::iterator i = triangles().begin();
      advance(i, triangles().size() - geom.triangles().size());
      while(i != triangles().end())
	{
	  for(triangle_t::iterator j = i->begin();
	      j != i->end();
	      j++)
	    *j += points().size() - geom.points().size();
	  i++;
	}
    }

    quads().insert(quads().end(),
                   geom.quads().begin(),
                   geom.quads().end());
    {
      quads_t::iterator i = quads().begin();
      advance(i, quads().size() - geom.quads().size());
      while(i != quads().end())
	{
	  for(quad_t::iterator j = i->begin();
	      j != i->end();
	      j++)
	    *j += points().size() - geom.points().size();
	  i++;
	}
    }

    return *this;
  }

  cvcgeom_t cvcgeom_t::tri_surface() const
  {
    cvcgeom_t new_geom;

    new_geom._points = _points;
    new_geom._colors = _colors;

    //if we don't have boundary info, just insert the tris directly
    if(boundary().size() != points().size())
      {
        new_geom._triangles = _triangles;

	for(quads_t::const_iterator i = quads().begin();
	    i != quads().end();
	    i++)
	  {
	    triangle_t quad_tris[2];
	    quad_tris[0][0] = (*i)[0];
	    quad_tris[0][1] = (*i)[1];
	    quad_tris[0][2] = (*i)[3];

	    quad_tris[1][0] = (*i)[1];
	    quad_tris[1][1] = (*i)[2];
	    quad_tris[1][2] = (*i)[3];

	    new_geom.triangles().push_back(quad_tris[0]);
	    new_geom.triangles().push_back(quad_tris[1]);
	  }
      }
    else
      {
	for(triangles_t::const_iterator i = triangles().begin();
	    i != triangles().end();
	    i++)
	  if(boundary()[(*i)[0]] &&
	     boundary()[(*i)[1]] &&
	     boundary()[(*i)[2]])
	    new_geom.triangles().push_back(*i);

	for(quads_t::const_iterator i = quads().begin();
	    i != quads().end();
	    i++)
	  if(boundary()[(*i)[0]] &&
	     boundary()[(*i)[1]] &&
	     boundary()[(*i)[2]] &&
	     boundary()[(*i)[3]])
	    {
	      triangle_t quad_tris[2];
	      quad_tris[0][0] = (*i)[0];
	      quad_tris[0][1] = (*i)[1];
	      quad_tris[0][2] = (*i)[3];
		
	      quad_tris[1][0] = (*i)[1];
	      quad_tris[1][1] = (*i)[2];
	      quad_tris[1][2] = (*i)[3];
		
	      new_geom.triangles().push_back(quad_tris[0]);
	      new_geom.triangles().push_back(quad_tris[1]);
	    }
      }
      
    return new_geom;
  }

  cvcgeom_t& cvcgeom_t::calculate_surf_normals()
  {
    using namespace std;
    cvcgeom_t tri_geom = this->tri_surface();
      
    //find all the triangles that each point is a part of
    map<point_t, list<unsigned int> > neighbor_tris;
    for(triangles_t::const_iterator i = tri_geom.const_triangles().begin();
	i != tri_geom.const_triangles().end();
	i++)
      for(triangle_t::const_iterator j = i->begin();
	  j != i->end();
	  j++)
	neighbor_tris[tri_geom.const_points()[*j]]
	  .push_back(distance(tri_geom.const_triangles().begin(),i));

    //now iterate over our points, and for each point push back a new normal vector
    normals().clear();
    for(points_t::const_iterator i = const_points().begin();
	i != const_points().end();
	i++)
      {
	vector_t norm = {{ 0.0, 0.0, 0.0 }};
	list<unsigned int> &neighbors = neighbor_tris[*i];
	for(list<unsigned int>::iterator j = neighbors.begin();
	    j != neighbors.end();
	    j++)
	  {
	    vector_t cur_tri_norm;
	    const triangle_t &tri = tri_geom.const_triangles()[*j];
	    point_t p[3] = { tri_geom.const_points()[tri[0]],
			     tri_geom.const_points()[tri[1]],
			     tri_geom.const_points()[tri[2]] };
	    vector_t v1 = {{ p[1][0] - p[0][0],
			     p[1][1] - p[0][1],
			     p[1][2] - p[0][2] }};
	    vector_t v2 = {{ p[2][0] - p[0][0], 
			     p[2][1] - p[0][1],
			     p[2][2] - p[0][2] }};
	    utility::cross(cur_tri_norm,v1,v2);
	    utility::normalize(cur_tri_norm);
	    for(int k = 0; k < 3; k++)
	      norm[k] += cur_tri_norm[k];
	  }
	if(!neighbors.empty())
	  for(int j = 0; j < 3; j++)
	    norm[j] /= neighbors.size();
	normals().push_back(norm);
      }

    return *this;
  }
  
  cvcgeom_t cvcgeom_t::generate_wire_interior() const
  {
    cvcgeom_t new_geom(*this);

    if(!boundary().empty() && 
       boundary().size() == points().size())
      {
	//tet mesh
	if(!triangles().empty())
	  {
	    new_geom.triangles().clear();
	    for(triangles_t::const_iterator j = triangles().begin();
		j != triangles().end();
		j++)
	      {
		//add the tri face boundary lines if both verts
		//of the potential line to be added do not lie 
		//on the boundary
		if(!boundary()[(*j)[0]] || !boundary()[(*j)[1]])
		  {
		    line_t line = {{ (*j)[0], (*j)[1] }};
		    new_geom.lines().push_back(line);
		  }
		if(!boundary()[(*j)[1]] || !boundary()[(*j)[2]])
		  {
		    line_t line = {{ (*j)[1], (*j)[2] }};
		    new_geom.lines().push_back(line);
		  }
		if(!boundary()[(*j)[2]] || !boundary()[(*j)[0]])
		  {
		    line_t line = {{ (*j)[2], (*j)[0] }};
		    new_geom.lines().push_back(line);
		  }
		  
		//add the triangle to new_geom if all verts are on the boundary
		if(boundary()[(*j)[0]] && boundary()[(*j)[1]] && boundary()[(*j)[2]])
		  new_geom.triangles().push_back(*j);
	      }
	  }
	//hex mesh
	if(!quads().empty())
	  {
	    new_geom.quads().clear();
	    for(quads_t::const_iterator j = quads().begin();
		j != quads().end();
		j++)
	      {
		//same as comment above except for quads...
		if(!boundary()[(*j)[0]] || !boundary()[(*j)[1]])
		  {
		    line_t line = {{ (*j)[0], (*j)[1] }};
		    new_geom.lines().push_back(line);
		  }
		if(!boundary()[(*j)[1]] || !boundary()[(*j)[2]])
		  {
		    line_t line = {{ (*j)[1], (*j)[2] }};
		    new_geom.lines().push_back(line);
		  }
		if(!boundary()[(*j)[2]] || !boundary()[(*j)[3]])
		  {
		    line_t line = {{ (*j)[2], (*j)[3] }};
		    new_geom.lines().push_back(line);
		  }
		if(!boundary()[(*j)[3]] || !boundary()[(*j)[0]])
		  {
		    line_t line = {{ (*j)[3], (*j)[0] }};
		    new_geom.lines().push_back(line);
		  }
		  
		if(boundary()[(*j)[0]] && boundary()[(*j)[1]] && 
		   boundary()[(*j)[2]] && boundary()[(*j)[3]])
		  new_geom.quads().push_back(*j);
	      }
	  }
      }

    return new_geom;
  }
  
  cvcgeom_t& cvcgeom_t::invert_normals()
  {
    for(normals_t::iterator i = normals().begin();
	i != normals().end();
	i++)
      for(int j = 0; j < 3; j++)
	(*i)[j] *= -1.0;
    return *this;
  }

  cvcgeom_t& cvcgeom_t::reorient()
  {
    using namespace std;

    //TODO: make this work for quads

    if(!triangles().empty())
      {
	//find all the triangles that each point is a part of
	map<point_t, list<unsigned int> > neighbor_tris;
	for(triangles_t::iterator i = triangles().begin();
	    i != triangles().end();
	    i++)
	  for(triangle_t::iterator j = i->begin();
	      j != i->end();
	      j++)
	    neighbor_tris[points()[*j]]
	      .push_back(distance(triangles().begin(),i));
	  
	/*
	  IGNORE THIS
	    
	  for each triangle (i),
	  for each point (j) in that triangle (i)
	  for each triangle (k) that point (j) is a part of
	  if the angle between the normal of k and i is obtuse,
	  reverse the vertex order of k
	*/
	  
	//if we don't have normals, lets calculate them
	if(points().size() != normals().size())
	  calculate_surf_normals();
	  
	for(points_t::iterator i = points().begin();
	    i != points().end();
	    i++)
	  {
	    //skip non boundary verts because their normals are null
	    if(!boundary()[distance(points().begin(),i)]) continue;
	      
	    list<unsigned int> &neighbors = neighbor_tris[*i];
	    for(list<unsigned int>::iterator j = neighbors.begin();
		j != neighbors.end();
		j++)
	      {
		triangle_t &tri = triangles()[*j];
		for(int k=0; k<3; k++)
		  if(utility::dot(normals()[tri[k]],
				  normals()[distance(points().begin(),i)]) < 0)
		    for(int l=0; l<3; l++)
		      normals()[tri[k]][l] *= -1.0;
	      }
	  }
      }
      
    return *this;
  }

  cvcgeom_t& cvcgeom_t::clear()
  {
    return (*this = cvcgeom_t());
  }

  cvcgeom_t& cvcgeom_t::project(const cvcgeom_t& input)
  {
#ifdef CVCRAW_GEOMETRY_ENABLE_PROJECT
    cvcgeom_t ref = input.tri_surface();
    project_verts::project(points().begin(),
                           points().end(),
                           ref.const_points().begin(),
                           ref.const_points().end(),
                           ref.const_triangles().begin(),
                           ref.const_triangles().end());
#else
# pragma message("cvcgeom_t::project() disabled")
#endif
    return *this;
  }

  void cvcgeom_t::init_ptrs()
  {
    if(!_points)
      _points.reset(new points_t);
    if(!_boundary)
      _boundary.reset(new boundary_t);
    if(!_normals)
      _normals.reset(new normals_t);
    if(!_colors)
      _colors.reset(new colors_t);
    if(!_lines)
      _lines.reset(new lines_t);
    if(!_triangles)
      _triangles.reset(new triangles_t);
    if(!_quads)
      _quads.reset(new quads_t);
  }

  void cvcgeom_t::calc_extents() const
  {
    using namespace std;
    
    if(empty())
      {
	fill(_min.begin(),_min.end(),numeric_limits<scalar_t>::max());
	fill(_max.begin(),_max.end(),-numeric_limits<scalar_t>::max());
	return;
      }
    
    _min = _max = points()[0];

    for(points_t::const_iterator i = points().begin();
	i != points().end();
	i++)
      {
	if(_min[0] > (*i)[0]) _min[0] = (*i)[0];
	if(_min[1] > (*i)[1]) _min[1] = (*i)[1];
	if(_min[2] > (*i)[2]) _min[2] = (*i)[2];
	if(_max[0] < (*i)[0]) _max[0] = (*i)[0];
	if(_max[1] < (*i)[1]) _max[1] = (*i)[1];
	if(_max[2] < (*i)[2]) _max[2] = (*i)[2];
      }

    _extents_set = true;

    return;
  }

  template<class container_ptr_t>
  void make_unique(container_ptr_t& cp)
  {
    if(!cp.unique())
      {
        container_ptr_t tmp(cp);
        cp.reset(new typename container_ptr_t::element_type(*tmp));
      }
  }

  void cvcgeom_t::pre_write(ARRAY_TYPE at)
  {
    //invalidate calculated extents
    _extents_set = false;

    switch(at)
      {
      case POINTS: make_unique(_points); break;
      case BOUNDARY: make_unique(_boundary); break;
      case NORMALS: make_unique(_normals); break;
      case COLORS: make_unique(_colors); break;
      case LINES: make_unique(_lines); break;
      case TRIANGLES: make_unique(_triangles); break;
      case QUADS: make_unique(_quads); break;
      }
  }

  // arand: written 4-11-2011
  //        directly read a cvc-raw type file into the data structure
  cvcgeom_t::cvcgeom_t(const std::string & filename) {


    std::string errors;
    boost::regex file_extension("^(.*)(\\.\\S*)$");
    boost::smatch what;

    if(boost::regex_match(filename, what, file_extension)) {
      
      //using namespace std;
      //cout << what[0] << " ::: " << what[1] << " ::: " << what[2] << endl;
      

      if (what[2].compare(".raw") == 0 ||
	  what[2].compare(".rawn") == 0 ||
	  what[2].compare(".rawnc") == 0 ||
	  what[2].compare(".rawc") == 0 ) { // raw data types
	read_raw(filename);
      } else if (what[2].compare(".off") == 0) { // off files
	read_off(filename);	
      } else {	 

	// code should never get here...
	// this isn't working... but eventually set it up
	//throw CVC_NAMESPACE::UnsupportedGeometryFileType(std::string(BOOST_CURRENT_FUNCTION) + 
	//						 std::string(": Cannot read ") + filename);
      }
    }    
  }

  // arand, 11-16-2011: added off reader
  //        this is not fully functional but it does handle the most
  //        common variants of off files...

  void cvcgeom_t::read_off(const std::string & filename) {

    using namespace std;
    using namespace boost; using boost::format;


    init_ptrs();
    _extents_set = false;


    std::ifstream inf(filename.c_str());
    if(!inf)
      throw runtime_error(string("Could not open ") + filename);


    unsigned int num_verts, num_elems;
    string line;
    vector<std::string> split_line;
    int line_num = 0;


    getline(inf, line); line_num++; // junk the first line
    if(!inf)
      throw runtime_error(str(format("Error reading file %1%, line %2%")
			      % filename
			      % line_num));


    //string headword;
    //inf >> headword;    
    //if (headword.compare("OFF") != 0 && 
    //	headword.compare("COFF") != 0) {
    //  throw runtime_error(str(format("Error reading header for file %1%")
    //		      % filename));
    //}

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
    if(split_line.size() != 3)
      throw runtime_error(str(format("Not an OFF file (wrong number of tokens in line 2: [%1%])")
			      % split_line.size()));


    try {
      num_verts = lexical_cast<unsigned int>(split_line[0]);
      num_elems = split_line.size() > 1 ?
	lexical_cast<unsigned int>(split_line[1]) : 0;
      
      for(unsigned int vt = 0; vt < num_verts; vt++) {
	
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

	switch(split_line.size()) {

	case 6: // colors
	case 7: // colors with alpha
	  // arand: ignoring transparency...
	  color_t color;
	  for(int i = 3; i < 6; i++)
	    color[i-3] = lexical_cast<double>(split_line[i]);
	  _colors->push_back(color);
	case 3: // no colors
	  point_t point;
	  for(int i = 0; i < 3; i++) {
	    point[i] = lexical_cast<double>(split_line[i]);
	  }
	  _points->push_back(point);

	  break;
	default:
	  throw runtime_error(str(format("Error reading file %1%, line %2%")
				  % filename
				  % line_num));
	}
      }

    }catch(std::exception& e) {
	  throw runtime_error(str(format("Error reading file %1%, line %2%, contents: '%3%', reason: %4%")
				  % filename
				  % line_num
				  % line
				  % string(e.what())));
	}
	

    for(unsigned int tri = 0; tri < num_elems; tri++) {
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


	int size = lexical_cast<int>(split_line[0]);

	// don't handle segments yet... just planar facets
	if (size < 3) {
	  throw runtime_error(str(format("Error reading file %1%, line %2%")
				  % filename
				  % line_num));
	}

	// arand: just triangulate every surface so that
	//        we don't have special cases...
	for (int i=2; i<split_line.size()-1; i++) {
	  triangle_t triangle;
	  triangle[0] = lexical_cast<unsigned int>(split_line[1]);
	  triangle[1] = lexical_cast<unsigned int>(split_line[i]);
	  triangle[2] = lexical_cast<unsigned int>(split_line[i+1]);
	  _triangles->push_back(triangle);
	}
    }


	

  }

  void cvcgeom_t::read_raw(const std::string & filename) {
    using namespace std;
    using namespace boost; using boost::format;
    
    init_ptrs();

    // arand: bug fix, 8-23-2011
    _extents_set = false;

    /*
    points_ptr_t    _points;
    boundary_ptr_t  _boundary;
    normals_ptr_t   _normals;
    colors_ptr_t    _colors;
    lines_ptr_t     _lines;
    triangles_ptr_t _triangles;
    quads_ptr_t     _quads;
    */

    std::ifstream inf(filename.c_str());
    if(!inf)
      throw runtime_error(string("Could not open ") + filename);
    
    unsigned int line_num = 0;
    string line;
    vector<std::string> split_line;
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
		  _points->push_back(point);
		  if(split_line.size() == 4)
		    _boundary->push_back(lexical_cast<int>(split_line[3]));
		}
		break;
	      case 6: // rawn or rawc
	      case 7: // rawn or rawc with boundary
		{
		  bool is_rawn = ends_with(filename,".rawn");
		  point_t point;
		  for(int i = 0; i < 3; i++)
		    point[i] = lexical_cast<double>(split_line[i]);
		  _points->push_back(point);

		  if(is_rawn)
		    {
		      vector_t normal;
		      for(int i = 3; i < 6; i++)
			normal[i-3] = lexical_cast<double>(split_line[i]);
		      _normals->push_back(normal);
		    }
		  else
		    {
		      color_t color;
		      for(int i = 3; i < 6; i++)
			color[i-3] = lexical_cast<double>(split_line[i]);
		      _colors->push_back(color);
		    }

		  if(split_line.size() == 7)
		    _boundary->push_back(lexical_cast<int>(split_line[6]));
		}
		break;
	      case 9:  // rawnc
	      case 10: // rawnc with boundary
		{
		  point_t point;
		  vector_t normal;
		  color_t color;

		  for(int i = 0; i < 3; i++)
		    point[i-0] = lexical_cast<double>(split_line[i]);
		  _points->push_back(point);

		  for(int i = 3; i < 6; i++)
		    normal[i-3] = lexical_cast<double>(split_line[i]);
		  _normals->push_back(normal);

		  for(int i = 6; i < 9; i++)
		    color[i-6] = lexical_cast<double>(split_line[i]);
		  _colors->push_back(color);

		  if(split_line.size() == 10)
		    _boundary->push_back(lexical_cast<int>(split_line[9]));
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
		  _lines->push_back(line);
		}
		break;
	      case 3: // tris
		{
		  triangle_t triangle;
		  for(unsigned int i = 0; i < split_line.size(); i++)
		    triangle[i] = lexical_cast<unsigned int>(split_line[i]);
		  _triangles->push_back(triangle);
		}
		break;
	      case 4: // quads or tetrahedrons
		{
		  //if we didn't collect boundary information earlier, then we have quads
		  if(_boundary->size() != _points->size())
		    {
		      quad_t quad;
		      for(unsigned int i = 0; i < split_line.size(); i++)
			quad[i] = lexical_cast<unsigned int>(split_line[i]);
		      _quads->push_back(quad);
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
			_triangles->push_back(tet_tris[i]);
		    }
		}
		break;
	      case 8: // hexahedrons
		{
		  //if we don't have boundary information, then something is wrong!
		  if(_boundary->size() != _points->size())
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
		    _quads->push_back(hex_quads[i]);
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


  }

  // 06/01/2012 -- initial implementation
  cvcgeom_t bunny()
  {
    cvcgeom_t geom;
#ifdef CVCRAW_GEOMETRY_ENABLE_BUNNY
    for(size_t i = 0; i < 34835; i++)
      {
        point_t point;
        vector_t normal;
        for(int j = 0; j < 3; j++)
          point[j] = BUNNY_VERTS[i*3+j];
        for(int j = 0; j < 3; j++)
          normal[j] = BUNNY_NORMS[i*3+j];
        geom.points().push_back(point);
        geom.normals().push_back(normal);
      }

    for(size_t i = 0; i < 69473; i++)
      {
        triangle_t tri;
        for(int j = 0; j < 3; j++)
          tri[j] = BUNNY_TRIS[i*3+j];
        geom.triangles().push_back(tri);
      }
#endif
    return geom;
  }
}
