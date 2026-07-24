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

#ifdef USING_TILING

#include <cvcraw_geometry/contours.h>
#include <cvcraw_geometry/utility.h>

#include <ContourTiler/Slice.h>
#include <ContourTiler/print_utils.h>

#ifndef DISABLE_CONVERSION
#include <cvcraw_geometry/cvcraw_geometry.h>
#endif

#ifdef CVCRAW_GEOMETRY_ENABLE_PROJECT
#include <cvcraw_geometry/project_verts.h>
#endif

#include <log4cplus/logger.h>
#include <log4cplus/loggingmacros.h>

#include <map>
#include <list>
#include <limits>
#include <fstream>
#include <iostream>
#include <stdexcept>

#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/foreach.hpp>


namespace CVCGEOM_NAMESPACE
{
  typedef cvcgeom_t::point_t     point_t;
  typedef cvcgeom_t::points_t    points_t;
  typedef cvcgeom_t::boundary_t  boundary_t;
  typedef cvcgeom_t::normals_t   normals_t;
  typedef cvcgeom_t::colors_t    colors_t;
  typedef cvcgeom_t::lines_t     lines_t;
  typedef cvcgeom_t::triangles_t triangles_t;
  typedef cvcgeom_t::quads_t     quads_t;

  contours_t::contours_t() {}

  contours_t::~contours_t() {}

  contours_t::contours_t(const std::list<std::string>& components, 
			 const std::vector<CONTOURTILER_NAMESPACE::Slice>& slices, 
			 int z_first, int z_last, const std::string& name) {

    using namespace std;
    using namespace boost; using boost::format;
    using namespace CONTOURTILER_NAMESPACE;

    _components = vector<string>(components.begin(), components.end());
    _slices = slices;
    _z_first = z_first;
    _z_last = z_last;
    _z_scale = 0.05;
    _name = name;
    
    refresh_geom();
  }

  void contours_t::refresh_geom()
  {
    log4cplus::Logger logger = log4cplus::Logger::getInstance("cvcraw_geometry.contours_t.refresh_geom");

    using namespace std;
    using namespace boost; using boost::format;
    using namespace CONTOURTILER_NAMESPACE;

    LOG4CPLUS_TRACE(logger, "1");
    _geom.clear();
    LOG4CPLUS_TRACE(logger, "2");

    std::unordered_map<Point_3, int> p2i;
    size_t nextIdx = 0;
    for (const auto& slice : _slices) {
      list<string> components;
      LOG4CPLUS_TRACE(logger, "3");
      slice.components(back_inserter(components));
      LOG4CPLUS_TRACE(logger, "4");
      for (const auto& component : components) {
	LOG4CPLUS_TRACE(logger, "5");
	for(Slice::Contour_const_iterator it = slice.begin(component); it != slice.end(component); ++it) {
	  const Polygon_2& p = (*it)->polygon();
	  for (Polygon_2::Edge_const_iterator eit = p.edges_begin(); eit != p.edges_end(); ++eit) {
	    Point_3 pts[2];
	    pts[0] = eit->source();
	    pts[1] = eit->target();
	    cvcgeom_t::line_t line;
	    // Create the line
	    for (int i = 0; i < 2; ++i) {
	      if (p2i.find(pts[i]) == p2i.end()) {
		p2i[pts[i]] = nextIdx++;
		point_t p;
		p[0] = pts[i].x();
		p[1] = pts[i].y();
		p[2] = (pts[i].z()+0.5) * _z_scale;
		_geom.points().push_back(p);
	      }
	      size_t idx = p2i[pts[i]];
	      line[i] = idx;
	    }
	    // add the line
	    _geom.lines().push_back(line);
	  }
	}
      }

    }
  }

}

#endif
