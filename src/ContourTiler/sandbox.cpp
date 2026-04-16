#include <algorithm>
#include <vector>
#include <utility>
#include <iostream>
#include <fstream>

#include <ContourTiler/common.h>
#include <ContourTiler/Contour.h>
#include <ContourTiler/reader_ser.h>
#include <ContourTiler/reader_gnuplot.h>
#include <ContourTiler/print_utils.h>
#include <ContourTiler/tiler_operations.h>
#include <ContourTiler/Segment_3_undirected.h>
#include <ContourTiler/augment.h>
#include <ContourTiler/polygon_utils.h>
#include <ContourTiler/sweep_line_visitors.h>
#include <ContourTiler/remove_contour_intersections.h>
#include <ContourTiler/offset_polygon.h>
#include <ContourTiler/polygon_difference.h>
#include <ContourTiler/interp.h>
#include <ContourTiler/Contour2.h>
#include <ContourTiler/Slice2.h>

#include <log4cplus/logger.h>
#include <log4cplus/loggingmacros.h>
#include <log4cplus/fileappender.h>
#include <log4cplus/loglevel.h>
#include <log4cplus/configurator.h>

#include <boost/unordered_set.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/foreach.hpp>

#include <CGAL/Boolean_set_operations_2.h>

using namespace std;
using namespace CONTOURTILER_NAMESPACE;

int main(int argc, char** argv) {
  log4cplus::PropertyConfigurator::doConfigure("log4cplus.properties");
  log4cplus::Logger logger = log4cplus::Logger::getInstance("sandbox");

  Polygon_2 P;
  ifstream in(argv[1]);
  double x, y;
  in >> x >> y;
  while (!in.eof()) {
    P.push_back(Point_2(x, y));
    in >> x >> y;
  }
  in.close();

  map<Point_2, Point_2> old2new;
  P = adjust_nonsimple_polygon(P, boost::lexical_cast<double>(argv[2]), old2new);
  LOG4CPLUS_INFO(logger, pp(P));
  // vector<Polygon_2> polygons;
  // split_nonsimple(P, back_inserter(polygons));
  // for (const auto& p : polygons) {
  //   LOG4CPLUS_INFO(logger, pp(p));
  // }

  return 0;
}
