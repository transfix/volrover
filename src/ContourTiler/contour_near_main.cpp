#include <ContourTiler/main.h>

#include <iostream>
#include <fstream>
#include <list>
#include <stdexcept>

#include <boost/lexical_cast.hpp>
#include <boost/foreach.hpp>

#include <ContourTiler/reader_gnuplot.h>
#include <ContourTiler/reader_ser.h>
#include <ContourTiler/contour_graph.h>
#include <ContourTiler/tiler.h>
#include <ContourTiler/tiler_operations.h>
#include <ContourTiler/print_utils.h>
#include <ContourTiler/Slice2.h>
#include <ContourTiler/cl_options.h>

#include <log4cplus/logger.h>
#include <log4cplus/loggingmacros.h>
#include <log4cplus/fileappender.h>
#include <log4cplus/loglevel.h>
#include <log4cplus/configurator.h>

using namespace std;
using namespace boost; using boost::format;
using namespace CONTOURTILER_NAMESPACE;

struct cn_options
{
  cn_options() : num(10) {}
  string fn;
  string data_dir;
  string component;
  int slice;
  int num;
};

cn_options parse_cn_options(int argc, char** argv)
{
  const string usage("Usage: ContourNear [options] file\n"
		     "Options:\n"
                     "  -n arg               data directory (required)\n"
                     "  -c arg               name of component to check (required)\n"
                     "  -s arg               slice number (required)\n"
                     "  -m arg               number of components to put in command (default=10)\n"
		     );

  cn_options o;
  int c;
  while ((c = getopt(argc, argv, "n:c:s:m:")) != -1) {
    switch(c) {
    case 'n':
      o.data_dir = string(optarg);
      break;
    case 'c':
      o.component = string(optarg);
      break;
    case 's':
      o.slice = lexical_cast<int>(optarg);
      break;
    case 'm':
      o.num = lexical_cast<int>(optarg);
      break;
    case '?':
      cout << usage << endl;
      exit(1);
      break;
    }
  }
  if (optind > argc-1) {
    cout << usage << endl;
    exit(1);
  }

  o.fn = argv[optind];

  return o;
}

bool file_exists(const string& fn)
{
  ifstream testfile(fn.c_str());
  if (testfile) {
    testfile.close();
    return true;
  }
  return false;
}

struct Component
{
  Component(string name, Number_type dist) : _name(name), _dist(dist) {}
  bool operator<(const Component& rhs) const {
    return _dist < rhs._dist;
  }
  string name() const { return _name; }
  Number_type dist() const { return _dist; }
  string _name;
  Number_type _dist;
};

Point_2 centroid(const Polygon_2& P)
{
  Point_2 sum(0,0);
  for (auto p_it = P.vertices_begin(); p_it != P.vertices_end(); ++p_it) {
    const Point_2& p = *p_it;
    sum = Point_2(sum.x() + p.x(), sum.y() + p.y());
  }
  return Point_2(sum.x()/P.size(), sum.y()/P.size());
}

int main(int argc, char** argv)
{
  if (file_exists("log4cplus.properties")) {
    log4cplus::PropertyConfigurator::doConfigure("log4cplus.properties");
  }
  else {
    log4cplus::BasicConfigurator::doConfigure();
  }

  log4cplus::Logger logger = log4cplus::Logger::getInstance("tiler.main");
  
  cn_options o = parse_cn_options(argc, argv);

  list<Contour_handle> contours;
  list<Contour_exception> exceptions;
  list<string> components, components_skip;
  // components.push_back(o.component);
  vector<string> empty_components;
  read_contours_ser(o.data_dir + "/" + o.fn, back_inserter(contours), back_inserter(exceptions),
                    o.slice, o.slice, -1,
                    components.begin(), components.end(),
                    components_skip.begin(), components_skip.end());

  for (const auto& e : exceptions) {
    LOG4CPLUS_WARN(logger, "Error in reading contours: " << e.what());
  }

  list<Point_2> ccents;
  for (const auto& contour : contours) {
    if (contour->info().object_name() == o.component) {
      ccents.push_back(centroid(contour->polygon()));
    }
  }

  list<Component> comps;
  for (const auto& contour : contours) {
    const string component = contour->info().object_name();
    const Point_2 cent = centroid(contour->polygon());
    for (const auto& ccent : ccents) {
      Component c(component, sqrt(CGAL::squared_distance(cent, ccent)));
      comps.push_back(c);
    }
  }
  comps.sort();
  comps.reverse();
  for (const auto& c : comps) {
    cout << c.name() << " dist=" << c.dist() << endl;
  }
  
  comps.reverse();
  set<string> include;
  for (const auto& c : comps) {
    if (include.size() < o.num && c.name() != o.component) {
      include.insert(c.name());
    }
  }

  cout << "./ContourTilerBin -e 1e-15 -C 0.01 -d output -o raw -z 0.05 -n data -s " << o.slice << " " << o.slice+1 << " ";
  cout << "-c " << o.component << " ";
  for (const auto& c : include) {
    cout << "-c " << c << " ";
  }
  cout << "AxonsDendrites" << endl;
}
