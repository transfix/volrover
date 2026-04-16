#include <ContourTiler/main.h>

#include <iostream>
#include <fstream>
#include <list>
#include <stdexcept>

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
using namespace CONTOURTILER_NAMESPACE;

void print_usage()
{
  cout << endl;
  cout << "Usage: contour_tiler [options] file1 [file2...]" << endl;
  cout << endl;
  cout << " I/O options" << endl;
  cout << "  -f, --format (dat | ser)        input file format" << endl;
  cout << "  -n, --indir directory           directory for input files" << endl;
  cout << "  -b, --basename basename         base name for output files" << endl;
  cout << "  -d, --outdir directory          directory for output files" << endl;
  cout << "  -c, --componentname name        name of component to tile" << endl;
  cout << "  -x, --componentskip name        name of component to skip" << endl;
  cout << "  -s, --slices start end          start and end slices (inclusive)" << endl;
  cout << "  -z, --zscale scale              scale factor for z" << endl;
  cout << "  -m, --multi                     multiple surfaces to tile" << endl;
  cout << "      --printcomponents           don't tile; output names of available components" << endl;
  cout << "      --color r g b               component color 0-255" << endl;
  cout << endl;
  cout << " Tiling options" << endl;
  cout << "  -r, --removeintersections       remove 3D intersections between constructed" << endl;
  cout << "                                  components.  Must also specify a value > 0" << endl;
  cout << "                                  for the -C option." << endl;
  cout << "  -e, --collinearepsilon          min area for a triangle formed by 3" << endl;
  cout << "                                  consecutive points in a contour to have" << endl;
  cout << "                                  to avoid discarding of the middle point" << endl;
  cout << "                                  as a collinearity.  Recommended value for" << endl;
  cout << "                                  neuropil data is 1e-18." << endl;
  cout << "  -C, --curationdelta             min distance between contours in the same slice" << endl;
  cout << "  -I, --2dintremove               only remove 2D contour intersections.  Output will go to" << endl;
  cout << "                                  identically-named ser files in output directory" << endl;
  cout << "  -v, --correspondoverlap overlap required overlap for contours to correspond" << endl;
  cout << endl;
  cout << " Debugging options" << endl;
  cout << "  -o, --outformat (raw | gnuplot) tile output format" << endl;
  cout << "  -i, --intformat (raw | gnuplot) tile intermediate output format" << endl << endl;
  cout << "  -p, --phaseformat (raw | gnuplot) tile phases output format" << endl << endl;
  cout << "  -B, --outputbottom              output bottom slice" << endl;
  cout << "  -T, --outputtop                 output top slice" << endl;
  cout << "  -O, --outputotv                 output OTV" << endl;
  cout << "  -P, --outputotvpairs            output OTV pairs" << endl;
  cout << "  -V, --outputvertexlabels        output vertex labels" << endl;
  cout << "  -U, --outputuntiledregions      output untiled regions" << endl;
  cout << "  -R, --defaultprecision prec     default debug output precision" << endl;
  cout << "  -S, --smoothingfactor prec      contour smoothing factor (-1 - no smoothing)" << endl;
  // cout << "  -F, --nointerpur                don't interpolate untiled regions" << endl;
  cout << "      --close componentname e     find all components that have a point " << endl;
  cout << "                                  within e of componentname" << endl;
}

int main2(int argc, char** argv)
{
  ifstream testfile("log4cplus.properties");
  if (testfile) {
    testfile.close();
    log4cplus::PropertyConfigurator::doConfigure("log4cplus.properties");
  }
  else {
    log4cplus::BasicConfigurator::doConfigure();
  }

  vector<Contour_handle> slice;
  read_contours_gnuplot2("../test_data/temp1.dat", back_inserter(slice), 61);
  Segment_2 seg(slice[1]->polygon()[0], slice[1]->polygon()[1]);
  std::cout << intersects_proper(seg, slice[0]->polygon()) << std::endl;
}

// These are here just temporarily to get 2d intersection removal working
// The canonical version of this function is found in tiler.cpp
Number_type my_round(Number_type d, int dec) {
  d *= pow((Number_type)10, (Number_type)dec);
  d = (d < 0) ? ceil(d-0.5) : floor(d+0.5);
  d /= pow((Number_type)10, (Number_type)dec);
  return d;
}

// The canonical version of this function is found in tiler.cpp
void my_round(Polygon_2& p, int dec) {
  Polygon_2::Vertex_iterator vit;
  for (vit = p.vertices_begin(); vit != p.vertices_end(); ++vit) {
    Point_2 pnt = *vit;
    pnt = Point_25_<Kernel>(my_round(pnt.x(), dec), my_round(pnt.y(), dec), pnt.z(), pnt.id());
    p.set(vit, pnt);
  }

  p = Polygon_2(p.vertices_begin(), unique(p.vertices_begin(), p.vertices_end()));
}

// The canonical version of this function is found in tiler.cpp
void my_round(Slice& slice, int dec) {
  list<string> components;
  slice.components(back_inserter(components));
  for (list<string>::const_iterator it = components.begin();
       it != components.end();
       ++it)
  {
    for (Slice::Contour_iterator cit = slice.begin(*it);
	 cit != slice.end(*it);
	 ++cit)
    {
      my_round((*cit)->polygon(), dec);
    }
  }
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

int main(int argc, char** argv)
{
  if (file_exists("log4cplus.properties")) {
    log4cplus::PropertyConfigurator::doConfigure("log4cplus.properties");
  }
  else {
    log4cplus::BasicConfigurator::doConfigure();
  }

  log4cplus::Logger logger = log4cplus::Logger::getInstance("tiler.main");

  srand(5);

  cl_options options;
  try {
    options = cl_parse(argc, argv);
  }
  catch (std::runtime_error& e) {
    print_usage();
    return 1;
  }

  if (options.options.contour_curation_delta() > options.options.z_scale()) {
    cout << endl << "Error: contour separation distance (" << options.options.contour_curation_delta() << 
      ") cannot exceed slice spacing (" << options.options.z_scale() << ")" << endl << endl;
    return 1;
  }

  boost::unordered_map<string, Color> comp2color;
  for (int i = 1; i < 5000; ++i)
  {
    stringstream ss;
    ss << setfill('0');
    ss << setw(3) << i;
    comp2color["a" + ss.str()] = Color(2*rand() / (double)RAND_MAX, 2*rand() / (double)RAND_MAX, 2*rand() / (double)RAND_MAX);
  }

  if (options.print_components)
  {
    list<Contour_handle> contours;
    list<Contour_exception> exceptions;
    // read_contours_ser(options.files[0], back_inserter(contours), back_inserter(exceptions),
    //     	      options.start, options.end, options.smoothing_factor, 
    //     	      options.components.begin(), options.components.end());
    read_contours_ser_new(options.files[0], back_inserter(contours), back_inserter(exceptions),
                          options.start, options.end, options.smoothing_factor, 
                          options.components.begin(), options.components.end(),
                          options.components_skip.begin(), options.components_skip.end());

    set<string> components;
    for (list<Contour_handle>::iterator it = contours.begin(); it != contours.end(); ++it)
    {
      const Contour_handle contour = *it;
      string component = contour->info().object_name();
      components.insert(component);
    }

    cout << "Components:" << endl;
    set<string>::const_iterator it = components.begin();
    for (; it != components.end(); ++it) {
      cout << "   " << *it << endl;
    }
  }
  else if (options.ser2)
  {
    throw logic_error("not supported");
    // list<Contour2_handle> contours;
    // list<Contour_exception> exceptions;
    // vector<string> empty_components;
    // read_contours_ser2(options.files[0], back_inserter(contours), back_inserter(exceptions),
    // 		       options.start, options.end, options.smoothing_factor, 
    // 		       options.components.begin(), options.components.end(),
    // 		       options.components_skip.begin(), options.components_skip.end());

    // for (list<Contour_exception>::const_iterator it = exceptions.begin(); it != exceptions.end(); ++it) {
    //   LOG4CPLUS_WARN(logger, "Error in reading contours: " << it->what());
    // }

    // for (list<Contour2_handle>::iterator it = contours.begin(); it != contours.end(); ++it) {
    //   LOG4CPLUS_TRACE(logger, "Read in contour: ");
    //   for (Contour2::Polygon_iterator it2 = (*it)->begin(); it2 != (*it)->end(); ++it2) {
    // 	LOG4CPLUS_TRACE(logger, "  " << pp(*it2));
    //   }
    // }
    // int slice_begin = options.start;
    // int slice_end = options.end + 1;
    // vector<Slice2> slices(slice_end - slice_begin);
    // for (list<Contour2_handle>::iterator it = contours.begin(); it != contours.end(); ++it) {
    //   Contour2_handle contour = *it;
    //   const int z = (int) contour->info().slice();
    //   const string component = contour->info().object_name();
    //   LOG4CPLUS_TRACE(logger, "Contour: " << z << " " << component);
    //   Slice2& slice = slices[z - slice_begin];
    //   slice[component] = contour;
    // }

    // if (options.td_int_remove_only) {
    //   typedef vector<Slice2>::iterator Slice_iter;
    //   for (Slice_iter it = slices.begin(); it != slices.end(); ++it) {
    // 	Slice2& s = *it;
    // 	validate(s);
    // 	remove_collinear(s, options.options.collinear_epsilon());
    // 	validate(s);

    // 	if (options.options.contour_curation_delta() > 0) {
    // 	  remove_intersections(s, options.options.contour_curation_delta());
    // 	  validate(s);
    // 	}

    // 	// cout << options.options.output_dir() << "/" << options.options.base_name() << "." << s.z() << endl;
    // 	// cout << options.options.output_dir() << "/Volumejosef." << setw(3) << setfill('0') << s.z() << endl;
    // 	stringstream ss;
    // 	// ss << options.options.output_dir() << "/Volumejosef." << setw(3) << setfill('0') << s.z();
    // 	ss << options.options.output_dir() << "/Volumejosef." << s.begin()->second->info().slice();
    // 	string fn = ss.str();
    // 	if (!file_exists(fn)) {
    // 	  ofstream out(fn.c_str());
    // 	  print_ser(out, s, 0.05);
    // 	  cout << "Wrote to " << fn << endl;
    // 	}
    // 	else {
    // 	  cout << "WARNING: " << fn << " exists!" << endl;
    // 	}
    //   }
    // }
    // else { // tile
    //   // clock_t start = clock();
    //   // tile(slices.begin(), slices.end(), comp2color, options.options);
    //   // cout << "Time: " << (clock() - start)/(double)CLOCKS_PER_SEC << endl;
    // }
  }
  else if (options.ser)
  {
    list<Contour_handle> contours;
    list<Contour_exception> exceptions;
    vector<string> empty_components;
    // read_contours_ser(options.files[0], back_inserter(contours), back_inserter(exceptions),
    //     	       options.start, options.end, options.smoothing_factor, 
    //     	       options.components.begin(), options.components.end(),
    //     	       options.components_skip.begin(), options.components_skip.end());
    read_contours_ser_new(options.files[0], back_inserter(contours), back_inserter(exceptions),
                          options.start, options.end, options.smoothing_factor, 
                          options.components.begin(), options.components.end(),
                          options.components_skip.begin(), options.components_skip.end());

    for (list<Contour_exception>::const_iterator it = exceptions.begin(); it != exceptions.end(); ++it) {
      LOG4CPLUS_WARN(logger, "Error in reading contours: " << it->what());
    }
    int slice_begin = options.start;
    int slice_end = options.end + 1;
    vector<Slice> slices(slice_end - slice_begin);
    for (list<Contour_handle>::iterator it = contours.begin(); it != contours.end(); ++it)
    {
      const Contour_handle contour = *it;
      int z = (int) contour->polygon()[0].z();
      string component = contour->info().object_name();
      LOG4CPLUS_TRACE(logger, "Read " << component << ": " << pp(contour->polygon()));
      Slice& slice = slices[z - slice_begin];
      slice.push_back(component, contour);
    }

    if (options.td_int_remove_only) {
      typedef vector<Slice>::iterator Slice_iter;
      for (Slice_iter it = slices.begin(); it != slices.end(); ++it) {
	Slice& s = *it;
	// if (!s.empty()) {
	//   LOG4CPLUS_DEBUG(logger, "Removing 2D contour intersections in slice " << s.z());
	// }
	s.validate();
	s.remove_collinear(options.options.collinear_epsilon());
	s.validate();
	// LOG4CPLUS_TRACE(logger, "Removed collinear points in slice " << s.z());

	if (options.options.contour_curation_delta() > 0) {
	  s.remove_intersections(options.options.contour_curation_delta());
	  s.validate();
	  // LOG4CPLUS_TRACE(logger, "Removed intersections in slice " << s.z());

	  my_round(s, 5);
	  s.validate();
	  // LOG4CPLUS_TRACE(logger, "Rounded in slice " << s.z());
	}

	// cout << options.options.output_dir() << "/" << options.options.base_name() << "." << s.z() << endl;
	// cout << options.options.output_dir() << "/Volumejosef." << setw(3) << setfill('0') << s.z() << endl;
	stringstream ss;
	// ss << options.options.output_dir() << "/Volumejosef." << setw(3) << setfill('0') << s.z();
	ss << options.options.output_dir() << "/Volumejosef." << s.z();
	string fn = ss.str();
	if (!file_exists(fn)) {
	  ofstream out(fn.c_str());
	  print_ser(out, s, 0.05);
	  cout << "Wrote to " << fn << endl;
	}
	else {
	  cout << "WARNING: " << fn << " exists!" << endl;
	}
      }
    }
    else { // tile
      clock_t start = clock();
      tile(slices.begin(), slices.end(), comp2color, options.options);
      cout << "Time: " << (clock() - start)/(double)CLOCKS_PER_SEC << endl;
    }
  }
  else if (options.options.remove_intersections())
  {
    string names[] = {"a001", "a020"};
    vector<vector<Contour_handle> > slices;
    for (int i = 0; i < options.files.size(); ++i)
    {
      static log4cplus::Logger logger = log4cplus::Logger::getInstance("tiler");
      LOG4CPLUS_INFO(logger, "reading " << options.files[i]);
      
      vector<Contour_handle> slice;
      Number_type z = i/2 + options.start;
      read_contours_gnuplot2(options.files[i], back_inserter(slice), z);
      vector<Contour_handle>::iterator it = slice.begin();
      for (; it != slice.end(); ++it) {
	(*it)->info().object_name() = names[i%2];
	(*it)->info().slice() = z;
      }
      slices.push_back(slice);
    }

    vector<Slice> slices_(2);
    slices_[0].insert(names[0], slices[0].begin(), slices[0].end());
    slices_[0].insert(names[1], slices[1].begin(), slices[1].end());
    slices_[1].insert(names[0], slices[2].begin(), slices[2].end());
    slices_[1].insert(names[1], slices[3].begin(), slices[3].end());
    tile(slices_.begin(), slices_.end(), comp2color, options.options);
  }
  else
  {
    vector<vector<Contour_handle> > slices;
    vector<Slice> slices_;
    for (int i = 0; i < options.files.size(); ++i)
    {
      static log4cplus::Logger logger = log4cplus::Logger::getInstance("tiler");
      LOG4CPLUS_INFO(logger, "reading " << options.files[i]);
      
      vector<Contour_handle> slice;
      Number_type z = i + options.start;
      if (options.options.multi())
	z = (i % 2) + 1;
      read_contours_gnuplot2(options.files[i], back_inserter(slice), z);
      slices.push_back(slice);
      Slice s;
//       s.insert("a001", slice.begin(), slice.end());
      s.insert("", slice.begin(), slice.end());
      slices_.push_back(s);
    }

//     vector<Slice> slices_(2);
//     slices_[0].insert("a002", slices[0].begin(), slices[0].end());
//     slices_[1].insert("a002", slices[1].begin(), slices[1].end());
    tile(slices_.begin(), slices_.end(), comp2color, options.options);

//     vector<Contour_handle> bac, tac;
//     augment(slices[0].begin(), slices[0].end(), slices[1].begin(), slices[1].end(), back_inserter(bac), back_inserter(tac));

// //     bottom.replace(component, bac.begin(), bac.end());
// //     top.replace(component, tac.begin(), tac.end());

// //     tile(slices[0].begin(), slices[0].end(), slices[1].begin(), slices[1].end(), options.options);
//     tile(bac.begin(), bac.end(), tac.begin(), tac.end(), options.options);
  }
  

//   test_graph();

}
