#include <iostream>

#include <log4cplus/logger.h>
#include <log4cplus/loggingmacros.h>
#include <log4cplus/fileappender.h>
#include <log4cplus/loglevel.h>
#include <log4cplus/configurator.h>

#include <boost/filesystem.hpp>
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/foreach.hpp>

#include <ContourTiler/common.h>
#include <CGAL/Plane_3.h>
#include <ContourTiler/print_utils.h>

// For triangulation
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangle_2.h>

#include <ContourTiler/ecs.h>

using namespace std;
// using namespace CONTOURTILER_NAMESPACE;

CONTOURTILER_BEGIN_NAMESPACE

// typedef CGAL::Bbox_3 Bbox_3;
typedef Kernel::Plane_3 Plane_3;
typedef Kernel::Vector_3 Vector_3;
typedef Kernel::Direction_3 Direction_3;
typedef Kernel::Triangle_2 Triangle_2;
typedef Kernel::Triangle_3 Triangle_3;

typedef CGAL::Triangulation_vertex_base_2<Kernel>                     Vb;
typedef CGAL::Triangulation_vertex_base_with_info_2<bool, Kernel, Vb>     Info;
typedef CGAL::Constrained_triangulation_face_base_2<Kernel>           Fb;
typedef CGAL::Triangulation_data_structure_2<Info,Fb>              TDS;
typedef CGAL::Exact_predicates_tag                               Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<Kernel, TDS, Itag> CDT;

typedef vector<int> iline;
typedef iline tri;
typedef vector<double> dline;

// typedef boost::unordered_map<Point_3, list<Point_3> > CMap;
// typedef boost::unordered_map<int, list<int> > CMap;
// typedef boost::unordered_map<int, boost::unordered_set<int> > CMap;

// This maps segments to their original order in the clipped triangle
// typedef boost::unordered_map<pair<int,int>, pair<int,int> > SegOrderMap;
typedef boost::unordered_set<Segment_2> Segments;
// typedef boost::unordered_set<pair<int,int> > Segments;

typedef boost::unordered_map<Point_3, int> Vertex2Index;

typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
const boost::char_separator<char> sep(" ");

enum Keep {
  KEEP,
  TOSS,
  UNKNOWN
};

void ecs_Bbox_3::add(const Point_3& p) {
  for (int i = 0; i < 3; ++i) {
    mins[i] = ::min(mins[i], p[i]);
    maxs[i] = ::max(maxs[i], p[i]);
  }
}

template <typename T>
vector<T> myreadline(istream& in)
{
  string line;
  getline(in, line);
  tokenizer tok(line, sep);
  vector<T> v;
  for (const auto& s : tok) {
    v.push_back(boost::lexical_cast<T>(s));
  }
  return v;
}

iline readiline(istream& in)
{
  return myreadline<int>(in);
}

iline readiline(istream& in, int offset)
{
  iline ret = myreadline<int>(in);
  for (int i = 0; i < ret.size(); ++i) {
    ret[i] += offset;
  }
  return ret;
}

dline readdline(istream& in)
{
  return myreadline<double>(in);
}

tri new_tri(int i, int j, int k)
{
  tri t(3);
  t[0] = i;
  t[1] = j;
  t[2] = k;
  return t;
}

bool is_inside(const Point_3& p, const ecs_Bbox_3& bb)
{
  const double e = -0.0000001;
  for (int i = 0; i < 3; ++i) {
    if ((bb.apply[i] && p[i] < bb.min(i)+e) || 
	(bb.apply[i+3] && p[i] > bb.max(i)-e)) 
      return false;
  }
  return true;
}

double intersection(const Plane_3& p, const Segment_3& seg, Point_3& ip)
{
  log4cplus::Logger logger = log4cplus::Logger::getInstance("ecs.intersection");

  Point_3 p0 = seg.source();
  Point_3 p1 = seg.target();
  Vector_3 v = p1 - p0;
  double t = -(p.a()*p0.x() + p.b()*p0.y() + p.c()*p0.z() + p.d())/(p.a()*v.x()+p.b()*v.y()+p.c()*v.z());
  ip = Point_3(p0.x()+t*v.x(), p0.y()+t*v.y(), p0.z()+t*v.z());
  LOG4CPLUS_TRACE(logger, "t = " << t);
  return t;
}

Plane_3 get_plane(int idx, const ecs_Bbox_3& bb)
{
  return Plane_3((idx%3)==0?1:0,
  		 (idx%3)==1?1:0,
  		 (idx%3)==2?1:0,
  		 (idx/3)==0?-bb.min(idx%3):-bb.max(idx%3));
}

Point_3 operator-(const Point_3& p0, const Point_3& p1) {
  return Point_3(p1.x()-p0.x(), p1.y()-p0.y(), p1.z()-p0.z());
}

Point_3 operator+(const Point_3& v0, const Point_3& v1) {
  return Point_3(v1.x()+v0.x(), v1.y()+v0.y(), v1.z()+v0.z());
}

Point_3 operator*(double d, const Point_3& v) {
  return Point_3(d*v.x(), d*v.y(), d*v.z());
}

Point_3 operator/(const Point_3& v, double d) {
  return Point_3(v.x()/d, v.y()/d, v.z()/d);
}

double dot(const Point_3& v0, const Point_3& v1) {
  return v0[0]*v1[0] + v0[1]*v1[1] + v0[2]*v1[2];
}

Point_3 cross(const Point_3& v0, const Point_3& v1) {
  double x = (v0[1]*v1[2]-v1[1]*v0[2]);
  // double y = v1[0]*v0[2]-v0[0]*v1[2];
  // double y = v0[0]*v1[2]-v1[0]*v0[2];
  double y = (v0[0]*v1[2]-v1[0]*v0[2])*-1;
  double z = (v0[0]*v1[1]-v1[0]*v0[1]);
  return Point_3(x, y, z);
}

double length(const Point_3& v) {
  return sqrt(dot(v, v));
}

Point_3 normalized(const Point_3& v) {
  return v/length(v);
}

Point_3 normal(const Triangle_3& t) {
  return cross(t[1]-t[0], t[2]-t[0]);
}

// Intersect two bounding box planes with a triangle
Point_3 intersect(const ecs_Bbox_3& bb, int p1, int p2, const tri& t, const vector<Point_3>& verts)
{
  log4cplus::Logger logger = log4cplus::Logger::getInstance("ecs.intersect2");

  // Find the plane equation for the triangle
  Triangle_3 t3(verts[t[0]], verts[t[1]], verts[t[2]]);
  Point_3 n = normalized(normal(t3));
  const double a = n.x();
  const double b = n.y();
  const double c = n.z();
  const double d = -(a*t3[0].x() + b*t3[0].y() + c*t3[0].z());

  LOG4CPLUS_TRACE(logger, "d = " << d << "  plane = " << a*t3[1].x() + b*t3[1].y() + c*t3[1].z());
  assert(fabs(d+(a*t3[1].x() + b*t3[1].y() + c*t3[1].z())) < .0000001);
  assert(fabs(d+(a*t3[2].x() + b*t3[2].y() + c*t3[2].z())) < .0000001);
  assert(p1%3 != p2%3);

  if (p2%3 < p1%3) {
    swap(p1, p2);
  }

  // Comment (1)
  // It is possible that the intersection point be on the outside of
  // the remaining plane.  If this is the case, bring it to the intersection
  // of all three relevant planes.

  double x, y, z;
  if (p1%3 == 0) {
    x = bb.val(p1);
    if (p2%3 == 1) {
      y = bb.val(p2);
      z = -(a*x + b*y + d)/c;
      // See comment (1)
      // z = max(z, bb.min(2));
      // z = min(z, bb.max(2));
    }
    else {
      z = bb.val(p2);
      y = -(a*x + c*z + d)/b;
      // See comment (1)
      // y = max(y, bb.min(1));
      // y = min(y, bb.max(1));
    }
  }
  else {
    y = bb.val(p1);
    z = bb.val(p2);
    x = -(b*y + c*z + d)/a;
    // See comment (1)
    // x = max(x, bb.min(0));
    // x = min(x, bb.max(0));
  }

  if (!is_inside(Point_3(x,y,z), bb)) {
    LOG4CPLUS_ERROR(logger, "Not inside bb: " << pp(Point_3(x,y,z)));
    throw logic_error("Not inside");
  }
  return Point_3(x,y,z);
}

// s.source() must be inside the bounding box
Point_3 intersect(const ecs_Bbox_3& bb, const Segment_3& s, int& iplane)
{
  log4cplus::Logger logger = log4cplus::Logger::getInstance("ecs.intersect");

  if (!is_inside(s.source(), bb)) {
    throw logic_error("Source must be inside bounding box");
  }
  if (is_inside(s.target(), bb)) {
    throw logic_error("Target must be outside bounding box");
  }

  double tmin = 2;
  Point_3 pmin;
  for (int i = 0; i < 6; ++i) {
    // const Plane_3& p = planes[i];
    Plane_3 p = get_plane(i, bb);
    Point_3 ip;
    double t = intersection(p, s, ip);
    if (t >= 0 && t <= 1) {
      // iplanes.push_back(i);
      if (t < tmin) {
	iplane = i;
	tmin = t;
	pmin = ip;
      }
    }
  }
  
  LOG4CPLUS_TRACE(logger, "intersected point = " << pp(pmin));
  assert(tmin != 2);
  if (!is_inside(pmin, bb)) {
    throw logic_error("Not inside");
  }

  return pmin;
}

double flip(int plane_idx) {
  if (plane_idx < 3)
    return plane_idx%3 == 1?-1:1;
  return plane_idx%3 != 1?-1:1;
}

Point_3 to3(const Point_2& p, ecs_Bbox_3 bb, int plane_idx)
{
  double plane_loc = bb.min(plane_idx%3);
  if (plane_idx > 2) {
    plane_loc = bb.max(plane_idx%3);
  }
  Point_2 p2(flip(plane_idx)*p.x(), p.y());
  double ret[3];
  int p2i = 0;
  for (int i = 0; i < 3; ++i) {
    if (i == plane_idx%3) {
      ret[i] = plane_loc;
    }
    else {
      ret[i] = p2[p2i++];
    }
  }
  if (!is_inside(Point_3(ret[0], ret[1], ret[2]), bb)) {
    throw logic_error("Not inside - to3");
  }
  return Point_3(ret[0], ret[1], ret[2]);
}

Point_2 to2(const Point_3& p3, int plane_idx)
{
  double ret[2];
  int p2i = 0;
  for (int i = 0; i < 3; ++i) {
    if (i != plane_idx%3) {
      ret[p2i++] = p3[i];
    }
  }
  ret[0] = flip(plane_idx) * ret[0];
  return Point_2(ret[0], ret[1]);
}

void create_boundaries(Polygon_2* boundaries, ecs_Bbox_3 bb)
{
  static log4cplus::Logger logger = log4cplus::Logger::getInstance("create_boundaries");

  double xmin = bb.min(0);
  double ymin = bb.min(1);
  double zmin = bb.min(2);
  double xmax = bb.max(0);
  double ymax = bb.max(1);
  double zmax = bb.max(2);
  double f;
  
  f = flip(0);
  boundaries[0].push_back(Point_2(f*ymin, zmin));
  boundaries[0].push_back(Point_2(f*ymax, zmin));
  boundaries[0].push_back(Point_2(f*ymax, zmax));
  boundaries[0].push_back(Point_2(f*ymin, zmax));

  f = flip(1);
  boundaries[1].push_back(Point_2(f*xmin, zmin));
  boundaries[1].push_back(Point_2(f*xmax, zmin));
  boundaries[1].push_back(Point_2(f*xmax, zmax));
  boundaries[1].push_back(Point_2(f*xmin, zmax));

  f = flip(2);
  boundaries[2].push_back(Point_2(f*xmin, ymin));
  boundaries[2].push_back(Point_2(f*xmax, ymin));
  boundaries[2].push_back(Point_2(f*xmax, ymax));
  boundaries[2].push_back(Point_2(f*xmin, ymax));

  f = flip(3);
  boundaries[3].push_back(Point_2(f*ymin, zmin));
  boundaries[3].push_back(Point_2(f*ymax, zmin));
  boundaries[3].push_back(Point_2(f*ymax, zmax));
  boundaries[3].push_back(Point_2(f*ymin, zmax));

  f = flip(4);
  boundaries[4].push_back(Point_2(f*xmin, zmin));
  boundaries[4].push_back(Point_2(f*xmax, zmin));
  boundaries[4].push_back(Point_2(f*xmax, zmax));
  boundaries[4].push_back(Point_2(f*xmin, zmax));

  f = flip(5);
  boundaries[5].push_back(Point_2(f*xmin, ymin));
  boundaries[5].push_back(Point_2(f*xmax, ymin));
  boundaries[5].push_back(Point_2(f*xmax, ymax));
  boundaries[5].push_back(Point_2(f*xmin, ymax));

  for (int i = 0; i < 6; ++i) {
    if (!boundaries[i].is_counterclockwise_oriented()) {
      LOG4CPLUS_TRACE(logger, "Reversing " << i);
      boundaries[i].reverse_orientation();
    }
  }
}

Keep keep(const Triangle_2& t, Segments& segments)
{
  static log4cplus::Logger logger = log4cplus::Logger::getInstance("keep");

  vector<Segment_2> shared;
  vector<Point_2> others;

  for (int i = 0; i < 3; ++i) {
    for (int j = i+1; j < 3; ++j) {
      int indices[] = {i, j};
      assert((j+(j-i))%3 != i && (j+(j-i))%3 != j);
      Point_2 other = t[(j+(j-i))%3];
      // Check segment to see if it's in the segment set.  If so, check
      // orientation and toss or keep accordingly.
      for (int k = 0; k < 2; ++k) {
	Segment_2 s(t[indices[k]], t[indices[1-k]]);
	if (segments.find(s) != segments.end()) {
	  if (CGAL::left_turn(s.source(), s.target(), other)) {
	    shared.push_back(s);
	    others.push_back(other);
	  }
	  else {
	    return KEEP;
	  }
	}
      }
    }
  }

  if (shared.empty())
    return UNKNOWN;

  int size = shared.size();
  if (size == 1) {
    Segment_2 s(shared[0]);
    segments.insert(Segment_2(others[0], s.target()));
    segments.insert(Segment_2(s.source(), others[0]));
  }
  else if (size == 2) {
    if (shared[0].source() == shared[1].target()) {
      segments.insert(Segment_2(shared[1].source(), shared[0].target()));
    }
    else {
      segments.insert(Segment_2(shared[0].source(), shared[1].target()));
    }
  }
  return TOSS;
}

// Returns the index of the verts
// int add_vertex(const Point_3& v, vector<Point_3>& verts, bool force_new = false)
int add_vertex(const Point_3& v, vector<Point_3>& verts, vector<int>& i2i, Vertex2Index& v2i, const ecs_Bbox_3& bb, bool force_new = false)
{
  static log4cplus::Logger logger = log4cplus::Logger::getInstance("ecs.add_vertex");

  if (!force_new) {
    if (!is_inside(v, bb)) {
      LOG4CPLUS_ERROR(logger, "Vertex outside: " << pp(v));
      throw logic_error("Vertex outside");
    }
    if (v2i.find(v) != v2i.end()) {
      return v2i[v];
    }
  }
  int index = verts.size();
  verts.push_back(v);
  i2i.push_back(-1);
  v2i[v] = index;
  return index;
}

void add_tri(const tri& t, vector<tri>& tris, const vector<Point_3>& verts, vector<int>& i2i)
{
  tris.push_back(t);
  for (int i = 0; i < 3; ++i) {
    i2i[t[i]] = t[i];
  }
}

void add_connection(int a, int b, Segments& segments, int plane, const vector<Point_3>& verts) {
  if (a != b) {
    segments.insert(Segment_2(to2(verts[a],plane), to2(verts[b],plane)));
  }
}

void triangulate(Polygon_2* boundaries, int plane_idx, vector<Point_3>& verts, vector<tri>& tris, 
		 vector<int>& i2i, Vertex2Index& v2i, 
		 ecs_Bbox_3 bb, 
		 Segments& segments, bool keep_all)
{
  static log4cplus::Logger logger = log4cplus::Logger::getInstance("ecs.triangulate");

  LOG4CPLUS_TRACE(logger, "Creating triangulation");
  CDT cdt;
  LOG4CPLUS_TRACE(logger, "  Boundary");
  cdt.insert(boundaries[plane_idx].vertices_begin(), boundaries[plane_idx].vertices_end());

  LOG4CPLUS_TRACE(logger, "  Adding constraints");
  for (const auto& s : segments) {
    cdt.insert_constraint(s.source(), s.target());
  }

  LOG4CPLUS_TRACE(logger, "Storing triangle vertices");
  // Loop through the triangulation and store the vertices of each triangle
  list<Triangle_2> triangles, temp, keepers;
  for (CDT::Finite_faces_iterator ffi = cdt.finite_faces_begin();
       ffi != cdt.finite_faces_end();
       ++ffi)
  {
    Triangle_2 t(ffi->vertex(0)->point(), 
		 ffi->vertex(1)->point(), 
		 ffi->vertex(2)->point());
    triangles.push_back(t);
  }
  
  if (!keep_all) {
    bool keep_going = true;
    while (keep_going) {
      // while (!triangles.empty()) {
      LOG4CPLUS_TRACE(logger, "Triangles to check: " << triangles.size());
      keep_going = false;
      for (const auto& t : triangles) {
	Keep k = keep(t, segments);
	if (k == KEEP) {
	  keepers.push_back(t);
	  keep_going = true;
	}
	else if (k == UNKNOWN) {
	  temp.push_back(t);
	}
	else {
	  keep_going = true;
	}
      }
      triangles = temp;
      temp.clear();
    }
  }
  keepers.insert(keepers.end(), triangles.begin(), triangles.end());

  for (const auto& t : keepers) {
    tri tr(3);
    for (int i = 0; i < 3; ++i) {
      Point_2 p2 = t[i];
      Point_3 p3 = to3(p2, bb, plane_idx);
      tr[i] = add_vertex(p3, verts, i2i, v2i, bb);
      // tr[i] = verts.size();
      // verts.push_back(p3);
    }
    add_tri(tr, tris, verts, i2i);
    // tris.push_back(tr);
  }
}

void process(tri& t, vector<Point_3>& verts, vector<tri>& tris, const ecs_Bbox_3& bb, Segments* segments, vector<int>& i2i, Vertex2Index& v2i)
{
  log4cplus::Logger logger = log4cplus::Logger::getInstance("ecs.process");

  // inside[i] = (0|1|2)
  vector<int> inside, outside;

  // Initialize inside (outside) vector.  These vectors contain indices
  // to the vertices inside (outside) the bounding box.
  for (int j = 0; j < 3; ++j) {
    if (is_inside(verts[t[j]], bb)) inside.push_back(j);
    else outside.push_back(j);
  }
  if (inside.size() > 0) {
    if (outside.empty()) {
      add_tri(t, tris, verts, i2i);
    }
    else if (outside.size() == 1) {
      //----------------------------
      // 1 outside
      // Will now have two triangles
      //----------------------------
      int outidx = t[outside[0]];
      Point_3 out = verts[outidx];
      int iplanes[2];
      for (int j = 0; j < 2; ++j) {
	Point_3 in1 = verts[t[inside[j]]];
	Point_3 in2 = verts[t[inside[(j+1)%2]]];
	Point_3 in3 = intersect(bb, Segment_3(in1, out), iplanes[j]);
	// Point_3 in3 = intersect(bb, Segment_3(in1, out), iplanes);
	// t[outside[0]] = add_vertex(in3, verts, i2i);
	t[outside[0]] = add_vertex(in3, verts, i2i, v2i, bb);
	add_tri(t, tris, verts, i2i);
	t[inside[j]] = t[outside[0]];
      }
      bool swap = !((inside[1] - inside[0]) > 1);
      if (iplanes[0] == iplanes[1]) {
	int plane = iplanes[0];
	if (swap)
	  add_connection(t[inside[1]], t[inside[0]], segments[plane], plane, verts);
	else
	  add_connection(t[inside[0]], t[inside[1]], segments[plane], plane, verts);
      }
      else {
	// triangle/plane/plane intersection
	try {
	Point_3 tppi = intersect(bb, iplanes[0], iplanes[1], t, verts);
	int tppiv = add_vertex(tppi, verts, i2i, v2i, bb);
	if (swap) {
	  add_connection(t[inside[1]], tppiv, segments[iplanes[1]], iplanes[1], verts);
	  add_connection(tppiv, t[inside[0]], segments[iplanes[0]], iplanes[0], verts);
	  add_tri(new_tri(t[inside[1]], tppiv, t[inside[0]]), tris, verts, i2i);
	}
	else {
	  add_connection(t[inside[0]], tppiv, segments[iplanes[0]], iplanes[0], verts);
	  add_connection(tppiv, t[inside[1]], segments[iplanes[1]], iplanes[1], verts);
	  add_tri(new_tri(t[inside[0]], tppiv, t[inside[1]]), tris, verts, i2i);
	}
	} catch(logic_error& e) {
	}
      }
    }
    else if (outside.size() == 2) {
      //----------------------------
      // 2 outside
      // Still have 1 triangle
      //----------------------------
      int inidx = t[inside[0]];
      Point_3 in = verts[inidx];
      int iplanes[2];
      for (int j = 0; j < 2; ++j) {
	Point_3 out = verts[t[outside[j]]];
	out = intersect(bb, Segment_3(in, out), iplanes[j]);
	// t[outside[j]] = add_vertex(out, verts, i2i);
	t[outside[j]] = add_vertex(out, verts, i2i, v2i, bb);
      }
      add_tri(t, tris, verts, i2i);
      bool swap = (outside[1] - outside[0]) > 1;
      LOG4CPLUS_TRACE(logger, "iplanes[0] = " << iplanes[0]);
      LOG4CPLUS_TRACE(logger, "iplanes[1] = " << iplanes[1]);
      if (iplanes[0] == iplanes[1]) {
	int plane = iplanes[0];
	if (swap)
	  add_connection(t[outside[1]], t[outside[0]], segments[plane], plane, verts);
	else
	  add_connection(t[outside[0]], t[outside[1]], segments[plane], plane, verts);
      }
      else {
	// triangle/plane/plane intersection
	try {
	Point_3 tppi = intersect(bb, iplanes[0], iplanes[1], t, verts);
	int tppiv = add_vertex(tppi, verts, i2i, v2i, bb);
	if (swap) {
	  add_connection(t[outside[1]], tppiv, segments[iplanes[1]], iplanes[1], verts);
	  add_connection(tppiv, t[outside[0]], segments[iplanes[0]], iplanes[0], verts);
	  add_tri(new_tri(t[outside[1]], tppiv, t[outside[0]]), tris, verts, i2i);
	}
	else {
	  add_connection(t[outside[0]], tppiv, segments[iplanes[0]], iplanes[0], verts);
	  add_connection(tppiv, t[outside[1]], segments[iplanes[1]], iplanes[1], verts);
	  add_tri(new_tri(t[outside[0]], tppiv, t[outside[1]]), tris, verts, i2i);
	}
	} catch(logic_error& e) {
	}
      }
    }
  }
}

void process_ecs(const vector<string>& filenames, string outfn, ecs_Bbox_3 bb, bool bb_init, bool crop_only)
{
  log4cplus::Logger logger = log4cplus::Logger::getInstance("ecs.main");

  const bool box_only = false;
  const bool keep_all = false;

  Segments segments[6];
  Polygon_2 boundaries[6];

  vector<Point_3> verts;
  // index to index
  // i2i[idx] = -1; // if it is not used by a triangle
  vector<int> i2i;
  Vertex2Index v2i;
  vector<tri> tris;
  int offset = 0;

  // Read each file and process triangles
  list<tri> to_process;
  for (const auto& filename : filenames) {
    bool success = true;
    try {
      LOG4CPLUS_INFO(logger, "Reading " << filename);
      LOG4CPLUS_TRACE(logger, "  Offset = " << offset);
      LOG4CPLUS_TRACE(logger, "  verts.size() = " << verts.size());
      ifstream in(filename.c_str());

      // Read header
      iline first = readiline(in);
      int num_verts = first[0];
      int num_tris = first[1];

      // Read vertices
      for (int i = 0; i < num_verts; ++i) {
	dline line = readdline(in);
	if (in.eof()) throw logic_error("End of file reached prematurely");
	Point_3 p(line[0], line[1], line[2]);
	// add_vertex(p, verts, i2i, true);
	add_vertex(p, verts, i2i, v2i, bb, true);
	if (!bb_init) {
	  bb.add(p);
	}
      }

      // Read triangles
      for (int j = 0; j < num_tris; ++j) {
	tri t = readiline(in, offset);
	if (in.eof()) throw logic_error("End of file reached prematurely");
	LOG4CPLUS_TRACE(logger, "Adding triangle: " << t[0] << " " << t[1] << " " << t[2]);
	to_process.push_back(t);
      }
    }
    catch(exception& e) {
      LOG4CPLUS_WARN(logger, "Failed to read " << filename << ": " << e.what());
      success = false;
    }

    offset = verts.size();
  }

  if (!bb_init) {
    bb.contract();
    LOG4CPLUS_INFO(logger, "Using bounding box: min = (" << 
		   bb.min(0) << ", " << bb.min(1) << ", " << bb.min(2) << ") max = (" <<
		   bb.max(0) << ", " << bb.max(1) << ", " << bb.max(2) << ")");
  }
  create_boundaries(boundaries, bb);

  LOG4CPLUS_INFO(logger, "Processing meshes");
  for (auto t : to_process) {
    LOG4CPLUS_TRACE(logger, "Processing triangle: " << t[0] << " " << t[1] << " " << t[2]);
    process(t, verts, tris, bb, segments, i2i, v2i);
  }

  // Triangulate the boundaries
  if (!crop_only) {
    if (box_only) {
      tris.clear();
    }
    for (int i = 0; i < 6; ++i) {
      if (bb.apply[i]) {
	LOG4CPLUS_TRACE(logger, "Triangulating boundary " << i);
	triangulate(boundaries, i, verts, tris, i2i, v2i, bb, segments[i], keep_all);
      }
    }
  }

  //-----------------------------
  // Write results to disk
  //-----------------------------
  LOG4CPLUS_INFO(logger, "Writing " << outfn);
  string ext = outfn.substr(outfn.length() - 3, 3);
  transform(ext.begin(), ext.end(), ext.begin(), ::tolower);
  bool off = (ext == "off");

  // Count the number of vertices
  int num_write_verts = 0;
  int num_write_tris = 0;
  for (int i = 0; i < verts.size(); ++i) {
    // If any of the coordinates are nan, then update i2i
    if (verts[i].x() != verts[i].x() ||
  	verts[i].y() != verts[i].y() ||
  	verts[i].z() != verts[i].z()) {
      i2i[i] = -2;
    }
    if (i2i[i] >= 0) {
      num_write_verts++;
    }
  }
  for (const auto& t : tris) {
    // Skip if i2i is -2 (vertex is nan)
    if (i2i[t[0]] >= 0 && i2i[t[1]] >= 0 && i2i[t[2]] >= 0) {
      num_write_tris++;
    }
  }
  stringstream ss;
  if (off) {
    ss << "OFF" << endl;
    ss << num_write_verts << " " << num_write_tris << " 0" << endl;
  }
  else {
    ss << num_write_verts << " " << num_write_tris << endl;
  }

  int cur_idx = 0;
  for (int i = 0; i < verts.size(); ++i) {
    if (i2i[i] >= 0) {
      i2i[i] = cur_idx++;
      const Point_3& v = verts[i];
      ss << v.x() << " " << v.y() << " " << v.z() << endl;
    }
  }
  for (const auto& t : tris) {
    // Skip if i2i is -2 (vertex is nan)
    if (i2i[t[0]] >= 0 && i2i[t[1]] >= 0 && i2i[t[2]] >= 0) {
      if (crop_only) {
  	if (off) {
  	  ss << "3 " << i2i[t[0]] << " " << i2i[t[1]] << " " << i2i[t[2]] << endl;
  	}
  	else {
  	  ss << i2i[t[0]] << " " << i2i[t[1]] << " " << i2i[t[2]] << endl;
  	}
      }
      // Flip normals
      else {
  	if (off) {
  	  ss << "3 " << i2i[t[1]] << " " << i2i[t[0]] << " " << i2i[t[2]] << endl;
  	}
  	else {
  	  ss << i2i[t[1]] << " " << i2i[t[0]] << " " << i2i[t[2]] << endl;
  	}
      }
    }
  }

  ofstream out(outfn.c_str());
  out << ss.str();
  out.close();
}

CONTOURTILER_END_NAMESPACE
