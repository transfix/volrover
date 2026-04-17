#include <ContourTiler/cl_options.h>
#include <ContourTiler/print_utils.h>
#include <log4cplus/logger.h>
#include <log4cplus/loggingmacros.h>

#include <boost/lexical_cast.hpp>

using namespace std;

CONTOURTILER_BEGIN_NAMESPACE

cl_options cl_parse(int argc, char** argv)
{
  log4cplus::Logger logger = log4cplus::Logger::getInstance("tiler.cl_parse");

  LOG4CPLUS_TRACE(logger, "1");
  vector<string> args;
  for (int i = 0; i < argc; ++i) {
    LOG4CPLUS_TRACE(logger, "2");
    args.push_back(argv[i]);
  }
  LOG4CPLUS_TRACE(logger, "3");
  return cl_parse(args);
  LOG4CPLUS_TRACE(logger, "4");
}

cl_options cl_parse(const vector<string>& argv)
{
  int argc = argv.size();
 
  cl_options options;

  if (argc < 2)
    throw std::runtime_error("Illegal usage");

  string input_dir;

  for (int i = 1; i < argc; ++i)
  {
    string arg(argv[i]);
    if (arg == "-b" || arg == "--basename")
    {
      options.options.base_name() = argv[++i];
    }
    else if (arg == "-f" || arg == "--format")
    {
      arg = argv[++i];
      if (arg == "dat")
	options.ser = false;
      else if (arg == "ser")
	options.ser = true;
      else if (arg == "ser2")
	options.ser2 = true;
    }
    else if (arg == "--color")
    {
      options.options.color() = true;
      arg = argv[++i];
      options.options.color_r() = boost::lexical_cast<int>(arg);
      arg = argv[++i];
      options.options.color_g() = boost::lexical_cast<int>(arg);
      arg = argv[++i];
      options.options.color_b() = boost::lexical_cast<int>(arg);
    }
    else if (arg == "-d" || arg == "--outdir")
    {
      options.options.output_dir() = argv[++i];
    }
    else if (arg == "-n" || arg == "--indir")
    {
      input_dir = argv[++i];
    }
    else if (arg == "-o" || arg == "--outformat")
    {
      arg = argv[++i];
      if (arg == "raw")
	options.options.output_raw() = true;
      else if (arg == "gnuplot")
	options.options.output_gnuplot() = true;
    }
    else if (arg == "-i" || arg == "--intformat")
    {
      arg = argv[++i];
      if (arg == "raw")
	options.options.output_intermediate_raw() = true;
      else if (arg == "gnuplot")
	options.options.output_intermediate_gnuplot() = true;
    }
    else if (arg == "-p" || arg == "--phaseformat")
    {
      arg = argv[++i];
      if (arg == "raw")
	options.options.output_phases_raw() = true;
      else if (arg == "gnuplot")
	options.options.output_phases_gnuplot() = true;
    }
    else if (arg == "-c" || arg == "--componentname")
    {
      arg = argv[++i];
      options.components.insert(arg);
    }
    else if (arg == "-x" || arg == "--componentskip")
    {
      arg = argv[++i];
      options.components_skip.insert(arg);
    }
    else if (arg == "-s" || arg == "--slices")
    {
      arg = argv[++i];
      options.start = boost::lexical_cast<int>(arg);
      arg = argv[++i];
      options.end = boost::lexical_cast<int>(arg);
    }
    else if (arg == "-z" || arg == "--zscale")
    {
      arg = argv[++i];
      options.options.z_scale() = boost::lexical_cast<double>(arg);
    }
    else if (arg == "-v" || arg == "--correspondoverlap")
    {
      arg = argv[++i];
      options.options.correspondence_overlap() = boost::lexical_cast<double>(arg);
    }
    else if (arg == "-e" || arg == "--collinearepsilon")
    {
      arg = argv[++i];
      options.options.collinear_epsilon() = boost::lexical_cast<Number_type>(arg);
    }
    else if (arg == "-C" || arg == "--curationdelta")
    {
      arg = argv[++i];
      options.options.contour_curation_delta() = boost::lexical_cast<Number_type>(arg);
    }
    else if (arg == "-I" || arg == "--2dintremove")
    {
      options.td_int_remove_only = true;
    }
    else if (arg == "-m" || arg == "--multi")
    {
      options.options.multi() = true;
    }
    else if (arg == "--printcomponents")
    {
      options.print_components = true;
    }
    else if (arg == "-r" || arg == "--removeintersections")
    {
      options.options.remove_intersections() = true;
    }
    else if (arg == "-B" || arg == "--outputbottom")
      options.options.output_bottom() = true;
    else if (arg == "-T" || arg == "--outputtop")
      options.options.output_top() = true;
    else if (arg == "-O" || arg == "--outputotv")
      options.options.output_otv() = true;
    else if (arg == "-P" || arg == "--outputotvpairs")
      options.options.output_otv_pairs() = true;
    else if (arg == "-V" || arg == "--outputvertexlabels")
      options.options.output_vertex_labels() = true;
    else if (arg == "-U" || arg == "--outputuntiledregions")
      options.options.output_untiled_regions_gnuplot() = true;
    else if (arg == "-R" || arg == "--defaultprecision")
    {
      set_default_pp_precision(boost::lexical_cast<int>(argv[++i]));
    }
    else if (arg == "-S" || arg == "--smoothingfactor")
      options.smoothing_factor = boost::lexical_cast<int>(argv[++i]);
    else if (arg == "-F" || arg == "--nointerpur")
      options.options.interp_untiled_regions() = false;
    else if (arg[0] == '-')
    {
      throw std::runtime_error("Unrecognized option");
    }
    else
    {
      if (!input_dir.empty())
	arg = input_dir + "/" + arg;
      options.files.push_back(arg);
    }
  }

  return options;
}

CONTOURTILER_END_NAMESPACE
