/*
  Copyright 2007-2025 The University of Texas at Austin

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

// volutils - Unified volume utility CLI
// Consolidates the ~60 separate VolUtils programs from VolumeRover
// into a single tool with subcommands.

#include <cvc/volume.h>
#include <cvc/volume_file_info.h>
#include <cvc/volume_file_io.h>
#include <cvc/volume_ops.h>
#include <cvc/geometry.h>
#include <cvc/geometry_file_io.h>
#include <cvc/algorithm.h>
#include <cvc/types.h>

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <stdexcept>
#include <iomanip>
#include <cmath>
#include <sstream>

namespace po = boost::program_options;
namespace fs = boost::filesystem;

// ── Helpers ──

static std::string type_to_string(cvc::data_type t)
{
  return std::string(cvc::data_type_strings[t]);
}

static cvc::data_type string_to_type(const std::string& s)
{
  // Map user-friendly names to internal type strings
  static const struct { const char* alias; cvc::data_type dt; } type_map[] = {
    { "UChar",    cvc::UChar },   { "unsigned char",  cvc::UChar },
    { "UShort",   cvc::UShort },  { "unsigned short", cvc::UShort },
    { "UInt",     cvc::UInt },    { "unsigned int",   cvc::UInt },
    { "Float",    cvc::Float },   { "float",          cvc::Float },
    { "Double",   cvc::Double },  { "double",         cvc::Double },
  };
  for (const auto& t : type_map)
    if (s == t.alias)
      return t.dt;
  throw std::runtime_error("Unknown voxel type: " + s +
    ". Valid types: UChar, UShort, UInt, Float, Double");
}

// ── info ──

static int cmd_info(int argc, char** argv)
{
  po::options_description desc("volutils info - display volume metadata");
  desc.add_options()
    ("help,h", "show help")
    ("input,i", po::value<std::string>()->required(), "input volume file");

  po::positional_options_description pos;
  pos.add("input", 1);

  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(desc).positional(pos).run(), vm);
  if (vm.count("help")) { std::cout << desc << "\n"; return 0; }
  po::notify(vm);

  cvc::volume_file_info vfi(vm["input"].as<std::string>());

  std::cout << "File:       " << vfi.filename() << "\n"
            << "Dimensions: " << vfi.XDim() << " x " << vfi.YDim() << " x " << vfi.ZDim() << "\n"
            << "BBox:       [" << vfi.XMin() << ", " << vfi.YMin() << ", " << vfi.ZMin() << "] - ["
                               << vfi.XMax() << ", " << vfi.YMax() << ", " << vfi.ZMax() << "]\n"
            << "Span:       " << vfi.XSpan() << " x " << vfi.YSpan() << " x " << vfi.ZSpan() << "\n"
            << "Variables:  " << vfi.numVariables() << "\n"
            << "Timesteps:  " << vfi.numTimesteps() << "\n";

  for (unsigned v = 0; v < vfi.numVariables(); ++v)
  {
    std::cout << "  Var " << v << ": type=" << vfi.voxelTypeStr(v)
              << " min=" << vfi.min(v,0) << " max=" << vfi.max(v,0) << "\n";
  }
  return 0;
}

// ── stats ──

static int cmd_stats(int argc, char** argv)
{
  po::options_description desc("volutils stats - compute volume statistics");
  desc.add_options()
    ("help,h", "show help")
    ("input,i", po::value<std::string>()->required(), "input volume file");

  po::positional_options_description pos;
  pos.add("input", 1);

  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(desc).positional(pos).run(), vm);
  if (vm.count("help")) { std::cout << desc << "\n"; return 0; }
  po::notify(vm);

  cvc::volume vol;
  vol.read(vm["input"].as<std::string>());
  cvc::volume_stats s = cvc::compute_stats(vol);

  std::cout << std::setprecision(12)
            << "Min:     " << s.min << "\n"
            << "Max:     " << s.max << "\n"
            << "Mean:    " << s.mean << "\n"
            << "StdDev:  " << s.std_dev << "\n"
            << "Voxels:  " << s.num_voxels << "\n";
  return 0;
}

// ── convert ──

static int cmd_convert(int argc, char** argv)
{
  po::options_description desc("volutils convert - convert between volume formats or voxel types");
  desc.add_options()
    ("help,h", "show help")
    ("input,i", po::value<std::string>()->required(), "input volume file")
    ("output,o", po::value<std::string>()->required(), "output volume file")
    ("type,t", po::value<std::string>(), "output voxel type (UChar, UShort, UInt, Float, Double)");

  po::positional_options_description pos;
  pos.add("input", 1);
  pos.add("output", 1);

  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(desc).positional(pos).run(), vm);
  if (vm.count("help")) { std::cout << desc << "\n"; return 0; }
  po::notify(vm);

  cvc::volume vol;
  vol.read(vm["input"].as<std::string>());

  if (vm.count("type"))
    vol.voxelType(string_to_type(vm["type"].as<std::string>()));

  vol.write(vm["output"].as<std::string>());
  return 0;
}

// ── add ──

static int cmd_add(int argc, char** argv)
{
  po::options_description desc("volutils add - add two volumes element-wise");
  desc.add_options()
    ("help,h", "show help")
    ("input,i", po::value<std::vector<std::string>>()->multitoken()->required(), "two input volume files")
    ("output,o", po::value<std::string>()->required(), "output volume file");

  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(desc).run(), vm);
  if (vm.count("help")) { std::cout << desc << "\n"; return 0; }
  po::notify(vm);

  auto inputs = vm["input"].as<std::vector<std::string>>();
  if (inputs.size() != 2)
    throw std::runtime_error("add requires exactly 2 input files");

  cvc::volume a, b;
  a.read(inputs[0]);
  b.read(inputs[1]);
  cvc::volume result = cvc::vol_add(a, b);
  result.write(vm["output"].as<std::string>());
  return 0;
}

// ── subtract ──

static int cmd_subtract(int argc, char** argv)
{
  po::options_description desc("volutils subtract - subtract second volume from first");
  desc.add_options()
    ("help,h", "show help")
    ("input,i", po::value<std::vector<std::string>>()->multitoken()->required(), "two input volume files")
    ("output,o", po::value<std::string>()->required(), "output volume file");

  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(desc).run(), vm);
  if (vm.count("help")) { std::cout << desc << "\n"; return 0; }
  po::notify(vm);

  auto inputs = vm["input"].as<std::vector<std::string>>();
  if (inputs.size() != 2)
    throw std::runtime_error("subtract requires exactly 2 input files");

  cvc::volume a, b;
  a.read(inputs[0]);
  b.read(inputs[1]);
  cvc::volume result = cvc::vol_subtract(a, b);
  result.write(vm["output"].as<std::string>());
  return 0;
}

// ── scale ──

static int cmd_scale(int argc, char** argv)
{
  po::options_description desc("volutils scale - multiply volume by scalar");
  desc.add_options()
    ("help,h", "show help")
    ("input,i", po::value<std::string>()->required(), "input volume file")
    ("output,o", po::value<std::string>()->required(), "output volume file")
    ("factor,f", po::value<double>()->required(), "scale factor");

  po::positional_options_description pos;
  pos.add("input", 1);
  pos.add("output", 1);

  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(desc).positional(pos).run(), vm);
  if (vm.count("help")) { std::cout << desc << "\n"; return 0; }
  po::notify(vm);

  cvc::volume vol;
  vol.read(vm["input"].as<std::string>());
  cvc::volume result = cvc::vol_scale(vol, vm["factor"].as<double>());
  result.write(vm["output"].as<std::string>());
  return 0;
}

// ── normalize ──

static int cmd_normalize(int argc, char** argv)
{
  po::options_description desc("volutils normalize - remap voxel values to [min, max]");
  desc.add_options()
    ("help,h", "show help")
    ("input,i", po::value<std::string>()->required(), "input volume file")
    ("output,o", po::value<std::string>()->required(), "output volume file")
    ("min", po::value<double>()->default_value(0.0), "new minimum")
    ("max", po::value<double>()->default_value(1.0), "new maximum");

  po::positional_options_description pos;
  pos.add("input", 1);
  pos.add("output", 1);

  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(desc).positional(pos).run(), vm);
  if (vm.count("help")) { std::cout << desc << "\n"; return 0; }
  po::notify(vm);

  cvc::volume vol;
  vol.read(vm["input"].as<std::string>());
  cvc::volume result = cvc::vol_normalize(vol, vm["min"].as<double>(), vm["max"].as<double>());
  result.write(vm["output"].as<std::string>());
  return 0;
}

// ── clip ──

static int cmd_clip(int argc, char** argv)
{
  po::options_description desc("volutils clip - zero voxels above threshold");
  desc.add_options()
    ("help,h", "show help")
    ("input,i", po::value<std::string>()->required(), "input volume file")
    ("output,o", po::value<std::string>()->required(), "output volume file")
    ("threshold,t", po::value<double>()->required(), "clipping threshold");

  po::positional_options_description pos;
  pos.add("input", 1);
  pos.add("output", 1);

  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(desc).positional(pos).run(), vm);
  if (vm.count("help")) { std::cout << desc << "\n"; return 0; }
  po::notify(vm);

  cvc::volume vol;
  vol.read(vm["input"].as<std::string>());
  cvc::volume result = cvc::vol_clip(vol, vm["threshold"].as<double>());
  result.write(vm["output"].as<std::string>());
  return 0;
}

// ── negate ──

static int cmd_negate(int argc, char** argv)
{
  po::options_description desc("volutils negate - negate all voxel values");
  desc.add_options()
    ("help,h", "show help")
    ("input,i", po::value<std::string>()->required(), "input volume file")
    ("output,o", po::value<std::string>()->required(), "output volume file");

  po::positional_options_description pos;
  pos.add("input", 1);
  pos.add("output", 1);

  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(desc).positional(pos).run(), vm);
  if (vm.count("help")) { std::cout << desc << "\n"; return 0; }
  po::notify(vm);

  cvc::volume vol;
  vol.read(vm["input"].as<std::string>());
  cvc::volume result = cvc::vol_negate(vol);
  result.write(vm["output"].as<std::string>());
  return 0;
}

// ── mask ──

static int cmd_mask(int argc, char** argv)
{
  po::options_description desc("volutils mask - zero voxels where mask is nonzero");
  desc.add_options()
    ("help,h", "show help")
    ("input,i", po::value<std::string>()->required(), "input volume file")
    ("mask,m", po::value<std::string>()->required(), "mask volume file")
    ("output,o", po::value<std::string>()->required(), "output volume file")
    ("inverse", "use inverse mask (zero where mask IS zero)");

  po::positional_options_description pos;
  pos.add("input", 1);
  pos.add("output", 1);

  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(desc).positional(pos).run(), vm);
  if (vm.count("help")) { std::cout << desc << "\n"; return 0; }
  po::notify(vm);

  cvc::volume vol, mask_vol;
  vol.read(vm["input"].as<std::string>());
  mask_vol.read(vm["mask"].as<std::string>());

  cvc::volume result = vm.count("inverse")
    ? cvc::vol_inverse_mask(vol, mask_vol)
    : cvc::vol_mask(vol, mask_vol);
  result.write(vm["output"].as<std::string>());
  return 0;
}

// ── downsample ──

static int cmd_downsample(int argc, char** argv)
{
  po::options_description desc("volutils downsample - reduce volume resolution");
  desc.add_options()
    ("help,h", "show help")
    ("input,i", po::value<std::string>()->required(), "input volume file")
    ("output,o", po::value<std::string>()->required(), "output volume file")
    ("factor,f", po::value<unsigned int>()->default_value(2), "downsample factor");

  po::positional_options_description pos;
  pos.add("input", 1);
  pos.add("output", 1);

  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(desc).positional(pos).run(), vm);
  if (vm.count("help")) { std::cout << desc << "\n"; return 0; }
  po::notify(vm);

  cvc::volume vol;
  vol.read(vm["input"].as<std::string>());
  unsigned int f = vm["factor"].as<unsigned int>();
  cvc::volume result = cvc::vol_downsample(vol, f, f, f);
  result.write(vm["output"].as<std::string>());
  return 0;
}

// ── bunny ──

static int cmd_bunny(int argc, char** argv)
{
  po::options_description desc("volutils bunny - output Stanford bunny geometry or SDF volume");
  desc.add_options()
    ("help,h", "show help")
    ("output,o", po::value<std::string>()->required(), "output file (.off for geometry, .rawiv/.mrc for volume)")
    ("volume", "output SDF volume instead of geometry")
    ("dims,d", po::value<unsigned int>()->default_value(64), "volume dimensions (cube)")
    ("padding,p", po::value<double>()->default_value(0.1), "bounding box padding factor");

  po::positional_options_description pos;
  pos.add("output", 1);

  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(desc).positional(pos).run(), vm);
  if (vm.count("help")) { std::cout << desc << "\n"; return 0; }
  po::notify(vm);

  std::string output = vm["output"].as<std::string>();

  // Load built-in bunny geometry (the .bunny extension triggers bunny_io)
  cvc::geometry bunny;
  bunny.read("builtin.bunny");

  if (vm.count("volume"))
  {
    // Compute SDF from bunny geometry
    unsigned int d = vm["dims"].as<unsigned int>();
    double pad = vm["padding"].as<double>();

    // Compute padded bounding box
    auto pmin = bunny.min_point();
    auto pmax = bunny.max_point();
    double ext[3];
    for (int i = 0; i < 3; ++i) ext[i] = pmax[i] - pmin[i];
    double max_ext = std::max({ext[0], ext[1], ext[2]});
    double half = max_ext * (1.0 + pad) * 0.5;
    double cx = (pmin[0] + pmax[0]) * 0.5;
    double cy = (pmin[1] + pmax[1]) * 0.5;
    double cz = (pmin[2] + pmax[2]) * 0.5;
    cvc::bounding_box bbox(cx - half, cy - half, cz - half,
                           cx + half, cy + half, cz + half);

    // Use sdf_library directly via the public API
#ifdef CVC_ENABLE_SDF
    cvc::volume sdf_vol = cvc::sdf(bunny, cvc::dimension(d, d, d), bbox);
    sdf_vol.write(output);
#else
    throw std::runtime_error("SDF support not enabled (CVC_ENABLE_SDF=OFF)");
#endif
    std::cout << "Wrote bunny SDF volume " << d << "^3 to " << output << "\n";
  }
  else
  {
    bunny.write(output);
    std::cout << "Wrote bunny geometry (" << bunny.num_points() << " verts, "
              << bunny.num_tris() << " tris) to " << output << "\n";
  }
  return 0;
}

// ── rotate ──

static int cmd_rotate(int argc, char** argv)
{
  po::options_description desc("volutils rotate - rotate volume around Z-axis");
  desc.add_options()
    ("help,h", "show help")
    ("input,i", po::value<std::string>()->required(), "input volume file")
    ("output,o", po::value<std::string>()->required(), "output volume file")
    ("angle,a", po::value<double>()->required(), "rotation angle in degrees")
    ("count,n", po::value<int>()->default_value(1), "number of evenly-spaced rotations");

  po::positional_options_description pos;
  pos.add("input", 1);
  pos.add("output", 1);

  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(desc).positional(pos).run(), vm);
  if (vm.count("help")) { std::cout << desc << "\n"; return 0; }
  po::notify(vm);

  cvc::volume vol;
  vol.read(vm["input"].as<std::string>());
  std::string output = vm["output"].as<std::string>();
  int count = vm["count"].as<int>();

  if (count <= 1)
  {
    double angle_rad = vm["angle"].as<double>() * M_PI / 180.0;
    cvc::volume result = cvc::vol_rotate_z(vol, angle_rad);
    result.write(output);
    std::cout << "Rotated " << vm["angle"].as<double>() << " degrees -> " << output << "\n";
  }
  else
  {
    // Generate count rotations evenly spaced around 360 degrees
    for (int n = 0; n < count; n++)
    {
      double angle_rad = 2.0 * M_PI * n / count;
      cvc::volume result = cvc::vol_rotate_z(vol, angle_rad);
      std::ostringstream ss;
      ss << output << "." << std::setw(4) << std::setfill('0') << n << ".rawiv";
      result.write(ss.str());
    }
    std::cout << "Wrote " << count << " rotations\n";
  }
  return 0;
}

// ── ssim ──

static int cmd_ssim(int argc, char** argv)
{
  po::options_description desc("volutils ssim - compute Structural Similarity Index (SSIM)");
  desc.add_options()
    ("help,h", "show help")
    ("input,i", po::value<std::vector<std::string>>()->multitoken()->required(), "two input volume files")
    ("output,o", po::value<std::string>(), "output SSIM map volume file")
    ("window,w", po::value<int>()->default_value(11), "Gaussian window size (odd)")
    ("sigma,s", po::value<double>()->default_value(1.5), "Gaussian sigma");

  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(desc).run(), vm);
  if (vm.count("help")) { std::cout << desc << "\n"; return 0; }
  po::notify(vm);

  auto inputs = vm["input"].as<std::vector<std::string>>();
  if (inputs.size() != 2)
    throw std::runtime_error("ssim requires exactly 2 input files");

  cvc::volume a, b;
  a.read(inputs[0]);
  b.read(inputs[1]);

  cvc::ssim_result result = cvc::vol_ssim(a, b,
    vm["window"].as<int>(), vm["sigma"].as<double>());

  std::cout << std::setprecision(12) << "Mean SSIM: " << result.mean_ssim << "\n";

  if (vm.count("output"))
  {
    result.ssim_map.write(vm["output"].as<std::string>());
    std::cout << "Wrote SSIM map to " << vm["output"].as<std::string>() << "\n";
  }
  return 0;
}

// ── project ──

static int cmd_project(int argc, char** argv)
{
  po::options_description desc("volutils project - forward ray projection of volume");
  desc.add_options()
    ("help,h", "show help")
    ("input,i", po::value<std::string>()->required(), "input volume file")
    ("output,o", po::value<std::string>()->required(), "output projection volume")
    ("angles,a", po::value<std::string>()->required(), "file with angles in degrees (one per line)")
    ("step", po::value<double>()->default_value(0.5), "ray step size");

  po::positional_options_description pos;
  pos.add("input", 1);
  pos.add("output", 1);

  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(desc).positional(pos).run(), vm);
  if (vm.count("help")) { std::cout << desc << "\n"; return 0; }
  po::notify(vm);

  // Read angles
  std::vector<double> angles;
  {
    std::ifstream f(vm["angles"].as<std::string>());
    if (!f) throw std::runtime_error("Cannot open angles file");
    double deg;
    while (f >> deg) angles.push_back(deg * M_PI / 180.0);
  }

  cvc::volume vol;
  vol.read(vm["input"].as<std::string>());
  cvc::volume result = cvc::vol_project(vol, angles, vm["step"].as<double>());
  result.write(vm["output"].as<std::string>());
  std::cout << "Projected " << angles.size() << " angles -> " << vm["output"].as<std::string>() << "\n";
  return 0;
}

// ── backproject ──

static int cmd_backproject(int argc, char** argv)
{
  po::options_description desc("volutils backproject - filtered back-projection (tomographic reconstruction)");
  desc.add_options()
    ("help,h", "show help")
    ("input,i", po::value<std::string>()->required(), "input projection volume")
    ("output,o", po::value<std::string>()->required(), "output reconstructed volume")
    ("angles,a", po::value<std::string>()->required(), "file with angles in degrees (one per line)")
    ("dim,d", po::value<unsigned int>()->required(), "output cube dimension")
    ("no-filter", "disable FFT ramp filter");

  po::positional_options_description pos;
  pos.add("input", 1);
  pos.add("output", 1);

  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(desc).positional(pos).run(), vm);
  if (vm.count("help")) { std::cout << desc << "\n"; return 0; }
  po::notify(vm);

  std::vector<double> angles;
  {
    std::ifstream f(vm["angles"].as<std::string>());
    if (!f) throw std::runtime_error("Cannot open angles file");
    double deg;
    while (f >> deg) angles.push_back(deg * M_PI / 180.0);
  }

  cvc::volume proj;
  proj.read(vm["input"].as<std::string>());
  bool filter = !vm.count("no-filter");
  cvc::volume result = cvc::vol_back_project(proj, angles,
    vm["dim"].as<unsigned int>(), filter);
  result.write(vm["output"].as<std::string>());
  std::cout << "Reconstructed " << vm["dim"].as<unsigned int>() << "^3 volume -> "
            << vm["output"].as<std::string>() << "\n";
  return 0;
}

// ── vol2img ──

static int cmd_vol2img(int argc, char** argv)
{
  po::options_description desc("volutils vol2img - export volume slices as images");
  desc.add_options()
    ("help,h", "show help")
    ("input,i", po::value<std::string>()->required(), "input volume file")
    ("dir,d", po::value<std::string>()->required(), "output directory")
    ("format,f", po::value<std::string>()->default_value("slice_%05d.png"), "filename pattern (printf-style)");

  po::positional_options_description pos;
  pos.add("input", 1);
  pos.add("dir", 1);

  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(desc).positional(pos).run(), vm);
  if (vm.count("help")) { std::cout << desc << "\n"; return 0; }
  po::notify(vm);

  cvc::volume vol;
  vol.read(vm["input"].as<std::string>());

  std::string dir = vm["dir"].as<std::string>();
  fs::create_directories(dir);

  cvc::vol_to_slices(vol, dir, vm["format"].as<std::string>());
  std::cout << "Exported " << vol.ZDim() << " slices to " << dir << "/\n";
  return 0;
}

// ── img2vol ──

static int cmd_img2vol(int argc, char** argv)
{
  po::options_description desc("volutils img2vol - import image stack into volume");
  desc.add_options()
    ("help,h", "show help")
    ("input,i", po::value<std::vector<std::string>>()->multitoken()->required(), "input image files (in Z order)")
    ("output,o", po::value<std::string>()->required(), "output volume file");

  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(desc).run(), vm);
  if (vm.count("help")) { std::cout << desc << "\n"; return 0; }
  po::notify(vm);

  auto inputs = vm["input"].as<std::vector<std::string>>();
  cvc::volume result = cvc::slices_to_volume(inputs);
  result.write(vm["output"].as<std::string>());
  std::cout << "Imported " << inputs.size() << " images -> " << vm["output"].as<std::string>() << "\n";
  return 0;
}

// ── rgba-merge ──

static int cmd_rgba_merge(int argc, char** argv)
{
  po::options_description desc("volutils rgba-merge - merge 4 volumes into RGBA");
  desc.add_options()
    ("help,h", "show help")
    ("input,i", po::value<std::vector<std::string>>()->multitoken()->required(), "4 input volume files (R G B A)")
    ("output,o", po::value<std::string>()->required(), "output volume file");

  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(desc).run(), vm);
  if (vm.count("help")) { std::cout << desc << "\n"; return 0; }
  po::notify(vm);

  auto inputs = vm["input"].as<std::vector<std::string>>();
  if (inputs.size() != 4)
    throw std::runtime_error("rgba-merge requires exactly 4 input files (R G B A)");

  cvc::volume r, g, b, a;
  r.read(inputs[0]); g.read(inputs[1]); b.read(inputs[2]); a.read(inputs[3]);
  cvc::volume result = cvc::vol_rgba_merge(r, g, b, a);
  result.write(vm["output"].as<std::string>());
  return 0;
}

// ── extract ──

static int cmd_extract(int argc, char** argv)
{
  po::options_description desc("volutils extract - extract variable/timestep from multi-variable volume");
  desc.add_options()
    ("help,h", "show help")
    ("input,i", po::value<std::string>()->required(), "input volume file")
    ("output,o", po::value<std::string>()->required(), "output volume file")
    ("var,v", po::value<unsigned int>()->default_value(0), "variable index")
    ("time,t", po::value<unsigned int>()->default_value(0), "timestep index");

  po::positional_options_description pos;
  pos.add("input", 1);
  pos.add("output", 1);

  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(desc).positional(pos).run(), vm);
  if (vm.count("help")) { std::cout << desc << "\n"; return 0; }
  po::notify(vm);

  cvc::volume vol;
  vol.read(vm["input"].as<std::string>(),
           vm["var"].as<unsigned int>(),
           vm["time"].as<unsigned int>());
  vol.write(vm["output"].as<std::string>());
  std::cout << "Extracted var=" << vm["var"].as<unsigned int>()
            << " time=" << vm["time"].as<unsigned int>()
            << " -> " << vm["output"].as<std::string>() << "\n";
  return 0;
}

// ── clamp-min ──

static int cmd_clamp_min(int argc, char** argv)
{
  po::options_description desc("volutils clamp-min - clamp voxels to a minimum value");
  desc.add_options()
    ("help,h", "show help")
    ("input,i", po::value<std::string>()->required(), "input volume file")
    ("output,o", po::value<std::string>()->required(), "output volume file")
    ("value", po::value<double>()->required(), "minimum value to clamp to");

  po::positional_options_description pos;
  pos.add("input", 1);
  pos.add("output", 1);

  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(desc).positional(pos).run(), vm);
  if (vm.count("help")) { std::cout << desc << "\n"; return 0; }
  po::notify(vm);

  cvc::volume vol;
  vol.read(vm["input"].as<std::string>());
  cvc::volume result = cvc::vol_clamp_min(vol, vm["value"].as<double>());
  result.write(vm["output"].as<std::string>());
  return 0;
}

// ── difference ──

static int cmd_difference(int argc, char** argv)
{
  po::options_description desc("volutils difference - absolute difference between two volumes");
  desc.add_options()
    ("help,h", "show help")
    ("input,i", po::value<std::vector<std::string>>()->multitoken()->required(), "two input volume files")
    ("output,o", po::value<std::string>()->required(), "output volume file");

  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(desc).run(), vm);
  if (vm.count("help")) { std::cout << desc << "\n"; return 0; }
  po::notify(vm);

  auto inputs = vm["input"].as<std::vector<std::string>>();
  if (inputs.size() != 2)
    throw std::runtime_error("difference requires exactly 2 input files");

  cvc::volume a, b;
  a.read(inputs[0]);
  b.read(inputs[1]);
  cvc::volume result = cvc::vol_difference(a, b);
  result.write(vm["output"].as<std::string>());
  return 0;
}

// ── average ──

static int cmd_average(int argc, char** argv)
{
  po::options_description desc("volutils average - compute element-wise average of N volumes");
  desc.add_options()
    ("help,h", "show help")
    ("input,i", po::value<std::vector<std::string>>()->multitoken()->required(), "input volume files (2 or more)")
    ("output,o", po::value<std::string>()->required(), "output volume file");

  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(desc).run(), vm);
  if (vm.count("help")) { std::cout << desc << "\n"; return 0; }
  po::notify(vm);

  auto inputs = vm["input"].as<std::vector<std::string>>();
  if (inputs.size() < 2)
    throw std::runtime_error("average requires at least 2 input files");

  cvc::volume sum;
  sum.read(inputs[0]);
  for (size_t i = 1; i < inputs.size(); ++i)
  {
    cvc::volume tmp;
    tmp.read(inputs[i]);
    sum = cvc::vol_add(sum, tmp);
  }
  cvc::volume result = cvc::vol_scale(sum, 1.0 / static_cast<double>(inputs.size()));
  result.write(vm["output"].as<std::string>());
  std::cout << "Averaged " << inputs.size() << " volumes -> " << vm["output"].as<std::string>() << "\n";
  return 0;
}

// ── sum ──

static int cmd_sum(int argc, char** argv)
{
  po::options_description desc("volutils sum - sum N volumes element-wise");
  desc.add_options()
    ("help,h", "show help")
    ("input,i", po::value<std::vector<std::string>>()->multitoken()->required(), "input volume files (2 or more)")
    ("output,o", po::value<std::string>()->required(), "output volume file");

  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(desc).run(), vm);
  if (vm.count("help")) { std::cout << desc << "\n"; return 0; }
  po::notify(vm);

  auto inputs = vm["input"].as<std::vector<std::string>>();
  if (inputs.size() < 2)
    throw std::runtime_error("sum requires at least 2 input files");

  cvc::volume result;
  result.read(inputs[0]);
  for (size_t i = 1; i < inputs.size(); ++i)
  {
    cvc::volume tmp;
    tmp.read(inputs[i]);
    result = cvc::vol_add(result, tmp);
  }
  result.write(vm["output"].as<std::string>());
  std::cout << "Summed " << inputs.size() << " volumes -> " << vm["output"].as<std::string>() << "\n";
  return 0;
}

// ── interpolate ──

static int cmd_interpolate(int argc, char** argv)
{
  po::options_description desc("volutils interpolate - linear interpolation between two volumes");
  desc.add_options()
    ("help,h", "show help")
    ("input,i", po::value<std::vector<std::string>>()->multitoken()->required(), "two input volume files")
    ("output,o", po::value<std::string>()->required(), "output file prefix")
    ("steps,n", po::value<int>()->default_value(5), "number of interpolation steps");

  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(desc).run(), vm);
  if (vm.count("help")) { std::cout << desc << "\n"; return 0; }
  po::notify(vm);

  auto inputs = vm["input"].as<std::vector<std::string>>();
  if (inputs.size() != 2)
    throw std::runtime_error("interpolate requires exactly 2 input files");

  cvc::volume a, b;
  a.read(inputs[0]);
  b.read(inputs[1]);

  int steps = vm["steps"].as<int>();
  std::string prefix = vm["output"].as<std::string>();

  for (int n = 1; n < steps; ++n)
  {
    double t = static_cast<double>(n) / static_cast<double>(steps - 1);
    cvc::volume sa = cvc::vol_scale(a, 1.0 - t);
    cvc::volume sb = cvc::vol_scale(b, t);
    cvc::volume result = cvc::vol_add(sa, sb);
    std::ostringstream ss;
    ss << prefix << n << ".rawiv";
    result.write(ss.str());
  }
  std::cout << "Wrote " << (steps - 1) << " interpolated volumes\n";
  return 0;
}

// ── compare ──

static int cmd_compare(int argc, char** argv)
{
  po::options_description desc("volutils compare - compare two volumes voxel-by-voxel");
  desc.add_options()
    ("help,h", "show help")
    ("input,i", po::value<std::vector<std::string>>()->multitoken()->required(), "two input volume files")
    ("threshold,t", po::value<double>()->default_value(0.01), "error threshold for comparison");

  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(desc).run(), vm);
  if (vm.count("help")) { std::cout << desc << "\n"; return 0; }
  po::notify(vm);

  auto inputs = vm["input"].as<std::vector<std::string>>();
  if (inputs.size() != 2)
    throw std::runtime_error("compare requires exactly 2 input files");

  cvc::volume a, b;
  a.read(inputs[0]);
  b.read(inputs[1]);

  double thresh = vm["threshold"].as<double>();

  if (a.XDim() != b.XDim() || a.YDim() != b.YDim() || a.ZDim() != b.ZDim())
  {
    std::cerr << "DIFFER: dimensions mismatch ("
              << a.XDim() << "x" << a.YDim() << "x" << a.ZDim() << " vs "
              << b.XDim() << "x" << b.YDim() << "x" << b.ZDim() << ")\n";
    return 1;
  }

  uint64_t total = a.XDim() * a.YDim() * a.ZDim();
  uint64_t diff_count = 0;
  double max_diff = 0.0;
  for (uint64_t i = 0; i < total; ++i)
  {
    double d = std::fabs(a(i) - b(i));
    if (d > thresh) ++diff_count;
    if (d > max_diff) max_diff = d;
  }

  if (diff_count == 0)
  {
    std::cout << "MATCH: volumes are identical within threshold " << thresh << "\n";
    return 0;
  }
  else
  {
    std::cout << "DIFFER: " << diff_count << " of " << total
              << " voxels differ (max diff=" << max_diff << ")\n";
    return 1;
  }
}

// ── bbox-shift ──

static int cmd_bbox_shift(int argc, char** argv)
{
  po::options_description desc("volutils bbox-shift - shift volume bounding box by offset");
  desc.add_options()
    ("help,h", "show help")
    ("input,i", po::value<std::string>()->required(), "input volume file")
    ("output,o", po::value<std::string>()->required(), "output volume file")
    ("dx", po::value<double>()->required(), "X offset")
    ("dy", po::value<double>()->required(), "Y offset")
    ("dz", po::value<double>()->required(), "Z offset");

  po::positional_options_description pos;
  pos.add("input", 1);
  pos.add("output", 1);

  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(desc).positional(pos).run(), vm);
  if (vm.count("help")) { std::cout << desc << "\n"; return 0; }
  po::notify(vm);

  cvc::volume vol;
  vol.read(vm["input"].as<std::string>());

  double dx = vm["dx"].as<double>();
  double dy = vm["dy"].as<double>();
  double dz = vm["dz"].as<double>();

  cvc::bounding_box bb = vol.boundingBox();
  bb.minx += dx; bb.maxx += dx;
  bb.miny += dy; bb.maxy += dy;
  bb.minz += dz; bb.maxz += dz;
  vol.boundingBox(bb);
  vol.write(vm["output"].as<std::string>());
  return 0;
}

// ── bbox-scale ──

static int cmd_bbox_scale(int argc, char** argv)
{
  po::options_description desc("volutils bbox-scale - scale volume bounding box");
  desc.add_options()
    ("help,h", "show help")
    ("input,i", po::value<std::string>()->required(), "input volume file")
    ("output,o", po::value<std::string>()->required(), "output volume file")
    ("sx", po::value<double>()->required(), "X scale factor")
    ("sy", po::value<double>()->required(), "Y scale factor")
    ("sz", po::value<double>()->required(), "Z scale factor");

  po::positional_options_description pos;
  pos.add("input", 1);
  pos.add("output", 1);

  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(desc).positional(pos).run(), vm);
  if (vm.count("help")) { std::cout << desc << "\n"; return 0; }
  po::notify(vm);

  cvc::volume vol;
  vol.read(vm["input"].as<std::string>());

  double sx = vm["sx"].as<double>();
  double sy = vm["sy"].as<double>();
  double sz = vm["sz"].as<double>();

  cvc::bounding_box bb = vol.boundingBox();
  bb.minx *= sx; bb.maxx *= sx;
  bb.miny *= sy; bb.maxy *= sy;
  bb.minz *= sz; bb.maxz *= sz;
  vol.boundingBox(bb);
  vol.write(vm["output"].as<std::string>());
  return 0;
}

// ── bbox-set ──

static int cmd_bbox_set(int argc, char** argv)
{
  po::options_description desc("volutils bbox-set - set volume bounding box explicitly");
  desc.add_options()
    ("help,h", "show help")
    ("input,i", po::value<std::string>()->required(), "input volume file")
    ("output,o", po::value<std::string>()->required(), "output volume file")
    ("xmin", po::value<double>()->required(), "minimum X")
    ("ymin", po::value<double>()->required(), "minimum Y")
    ("zmin", po::value<double>()->required(), "minimum Z")
    ("xmax", po::value<double>()->required(), "maximum X")
    ("ymax", po::value<double>()->required(), "maximum Y")
    ("zmax", po::value<double>()->required(), "maximum Z");

  po::positional_options_description pos;
  pos.add("input", 1);
  pos.add("output", 1);

  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(desc).positional(pos).run(), vm);
  if (vm.count("help")) { std::cout << desc << "\n"; return 0; }
  po::notify(vm);

  cvc::volume vol;
  vol.read(vm["input"].as<std::string>());

  cvc::bounding_box bb(
    vm["xmin"].as<double>(), vm["ymin"].as<double>(), vm["zmin"].as<double>(),
    vm["xmax"].as<double>(), vm["ymax"].as<double>(), vm["zmax"].as<double>());
  vol.boundingBox(bb);
  vol.write(vm["output"].as<std::string>());
  return 0;
}

// ── fill ──

static int cmd_fill(int argc, char** argv)
{
  po::options_description desc("volutils fill - create volume filled with constant value");
  desc.add_options()
    ("help,h", "show help")
    ("output,o", po::value<std::string>()->required(), "output volume file")
    ("xdim", po::value<unsigned int>()->default_value(64), "X dimension")
    ("ydim", po::value<unsigned int>()->default_value(64), "Y dimension")
    ("zdim", po::value<unsigned int>()->default_value(64), "Z dimension")
    ("value", po::value<double>()->default_value(0.0), "fill value")
    ("type,t", po::value<std::string>()->default_value("Float"), "voxel type");

  po::positional_options_description pos;
  pos.add("output", 1);

  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(desc).positional(pos).run(), vm);
  if (vm.count("help")) { std::cout << desc << "\n"; return 0; }
  po::notify(vm);

  unsigned int xd = vm["xdim"].as<unsigned int>();
  unsigned int yd = vm["ydim"].as<unsigned int>();
  unsigned int zd = vm["zdim"].as<unsigned int>();
  double val = vm["value"].as<double>();

  cvc::volume vol(cvc::dimension(xd, yd, zd),
                  string_to_type(vm["type"].as<std::string>()));
  uint64_t total = static_cast<uint64_t>(xd) * yd * zd;
  for (uint64_t i = 0; i < total; ++i)
    vol(i, val);

  vol.write(vm["output"].as<std::string>());
  std::cout << "Created " << xd << "x" << yd << "x" << zd
            << " volume filled with " << val << "\n";
  return 0;
}

// ── edge ──

static int cmd_edge(int argc, char** argv)
{
  po::options_description desc("volutils edge - Laplacian of Gaussian edge detection");
  desc.add_options()
    ("help,h", "show help")
    ("input,i", po::value<std::string>()->required(), "input volume file")
    ("output,o", po::value<std::string>()->required(), "output volume file")
    ("sigma,s", po::value<double>()->default_value(1.0), "Gaussian sigma");

  po::positional_options_description pos;
  pos.add("input", 1);
  pos.add("output", 1);

  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(desc).positional(pos).run(), vm);
  if (vm.count("help")) { std::cout << desc << "\n"; return 0; }
  po::notify(vm);

  cvc::volume vol;
  vol.read(vm["input"].as<std::string>());

  double sigma = vm["sigma"].as<double>();

  // Build 5x5x5 LoG kernel
  double kernel[125];
  int idx = 0;
  double kernel_sum = 0.0;
  double sigma2 = sigma * sigma;
  double norm = 1.0 / (std::pow(2.0 * M_PI, 1.5) * std::pow(sigma, 5));
  for (int di = -2; di <= 2; ++di)
    for (int dj = -2; dj <= 2; ++dj)
      for (int dk = -2; dk <= 2; ++dk)
      {
        double r2 = di*di + dj*dj + dk*dk;
        double val = norm * (r2 / sigma2 - 3.0) * std::exp(-r2 / (2.0 * sigma2));
        kernel[idx++] = val;
        kernel_sum += val;
      }
  // Zero-center the kernel
  for (int n = 0; n < 125; ++n)
    kernel[n] -= kernel_sum / 125.0;

  uint64_t xd = vol.XDim(), yd = vol.YDim(), zd = vol.ZDim();
  cvc::volume result(cvc::dimension(xd, yd, zd), vol.voxelType(), vol.boundingBox());

  for (uint64_t x = 0; x < xd; ++x)
    for (uint64_t y = 0; y < yd; ++y)
      for (uint64_t z = 0; z < zd; ++z)
      {
        double val = 0.0;
        int ki = 0;
        for (int di = -2; di <= 2; ++di)
          for (int dj = -2; dj <= 2; ++dj)
            for (int dk = -2; dk <= 2; ++dk)
            {
              int64_t nx = static_cast<int64_t>(x) + di;
              int64_t ny = static_cast<int64_t>(y) + dj;
              int64_t nz = static_cast<int64_t>(z) + dk;
              if (nx >= 0 && nx < static_cast<int64_t>(xd) &&
                  ny >= 0 && ny < static_cast<int64_t>(yd) &&
                  nz >= 0 && nz < static_cast<int64_t>(zd))
                val += vol(nx, ny, nz) * kernel[ki];
              ++ki;
            }
        result(x, y, z, val);
      }

  result.write(vm["output"].as<std::string>());
  std::cout << "Edge detection (sigma=" << sigma << ") -> " << vm["output"].as<std::string>() << "\n";
  return 0;
}

// ── Main dispatcher ──

struct command_entry {
  const char* name;
  const char* help;
  int (*func)(int, char**);
};

static const command_entry commands[] = {
  { "info",       "display volume file metadata",         cmd_info },
  { "stats",      "compute volume statistics",            cmd_stats },
  { "convert",    "convert format or voxel type",         cmd_convert },
  { "extract",    "extract variable/timestep",            cmd_extract },
  { "add",        "add two volumes element-wise",         cmd_add },
  { "subtract",   "subtract two volumes element-wise",    cmd_subtract },
  { "difference", "absolute difference between volumes",  cmd_difference },
  { "sum",        "sum N volumes element-wise",           cmd_sum },
  { "average",    "average N volumes element-wise",       cmd_average },
  { "scale",      "multiply volume by scalar",            cmd_scale },
  { "normalize",  "remap voxel values to [min, max]",     cmd_normalize },
  { "interpolate","linear interpolation between volumes", cmd_interpolate },
  { "clip",       "zero voxels above threshold",          cmd_clip },
  { "clamp-min",  "clamp voxels to minimum value",       cmd_clamp_min },
  { "negate",     "negate all voxel values",              cmd_negate },
  { "mask",       "apply mask volume",                    cmd_mask },
  { "fill",       "create constant-filled volume",        cmd_fill },
  { "downsample", "reduce volume resolution",             cmd_downsample },
  { "rotate",     "rotate volume around Z-axis",          cmd_rotate },
  { "edge",       "LoG edge detection",                   cmd_edge },
  { "compare",    "compare two volumes voxel-by-voxel",   cmd_compare },
  { "ssim",       "compute SSIM between two volumes",     cmd_ssim },
  { "bbox-shift", "shift bounding box by offset",         cmd_bbox_shift },
  { "bbox-scale", "scale bounding box",                   cmd_bbox_scale },
  { "bbox-set",   "set bounding box explicitly",          cmd_bbox_set },
  { "project",    "forward ray projection",               cmd_project },
  { "backproject","filtered back-projection (FBP)",       cmd_backproject },
  { "vol2img",    "export slices as images (ImageMagick)", cmd_vol2img },
  { "img2vol",    "import image stack as volume",         cmd_img2vol },
  { "rgba-merge", "merge 4 volumes into RGBA",            cmd_rgba_merge },
  { "bunny",      "output Stanford bunny geometry or SDF", cmd_bunny },
};

static void print_usage()
{
  std::cout << "Usage: volutils <command> [options]\n\n"
            << "Commands:\n";
  for (const auto& c : commands)
    std::cout << "  " << std::left << std::setw(14) << c.name << c.help << "\n";
  std::cout << "\nRun 'volutils <command> --help' for command-specific options.\n";
}

int main(int argc, char** argv)
{
  if (argc < 2)
  {
    print_usage();
    return 1;
  }

  std::string cmd = argv[1];
  if (cmd == "--help" || cmd == "-h")
  {
    print_usage();
    return 0;
  }

  for (const auto& c : commands)
  {
    if (cmd == c.name)
    {
      try
      {
        return c.func(argc - 1, argv + 1);
      }
      catch (const std::exception& e)
      {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
      }
    }
  }

  std::cerr << "Unknown command: " << cmd << "\n";
  print_usage();
  return 1;
}
