// References
// 1. R. A. Jarvis. On the identification of the convex hull of a finite set of points in the plane. Inform. Process. Lett., 2:18-21, 1973
// 2. J. O'Rourke. Computational Geometry in C. Cambridge University Press, 1994.
// 3. F. P. Preparata and M. I. Shamos. Computational Geometry: An Introduction. Springer-Verlag, New York, 1985.

#include <iostream>
#include <fstream>
#include <deque>

#include <log4cplus/logger.h>
#include <log4cplus/loggingmacros.h>
#include <log4cplus/fileappender.h>
#include <log4cplus/loglevel.h>
#include <log4cplus/configurator.h>

#include <boost/lexical_cast.hpp>

#include <ContourTiler/reader_gnuplot.h>
#include <ContourTiler/print_utils.h>
#include <ContourTiler/offset_polygon.h>
#include <ContourTiler/kernel_utils.h>

CONTOURTILER_BEGIN_NAMESPACE

using namespace std;

// typedef Polygon_2::Vertex_circulator Circ;
typedef CGAL::Direction_2<Kernel> Direction_2;
typedef CGAL::Ray_2<Kernel> Ray_2;

typedef Bso_polygon_2::Vertex_circulator Circ;
typedef Bso_kernel::Ray_2 Bso_ray_2;
//typedef Bso_kernel::Segment_2 Bso_segment_2;
typedef Bso_kernel::Direction_2 Bso_direction_2;

// CGAL::Object intersection(const Ray_2& r_, const Segment_2& s_)
// {
//   typedef Bso_kernel K;
//   typedef K::Ray_2 Bso_ray_2;
//   typedef K::Segment_2 Bso_segment_2;

//   Bso_ray_2 r = change_kernel<K>(r_);
//   Bso_segment_2 s = change_kernel<K>(s_);
//   CGAL::Object o = CGAL::intersection(r, s);
//   const Bso_point_2* ip = CGAL::object_cast<Bso_point_2>(&o);
//   const Bso_segment_2* is = CGAL::object_cast<Bso_segment_2>(&o);
//   CGAL::Object ret;
//   if (ip)
//     ret = CGAL::make_object(to_common<K>(*ip, 0));
//   else if (is)
//     ret = CGAL::make_object(to_common<K>(*is, 0));
//   return ret;
// }

// CGAL::Object intersection(const Segment_2& s1_, const Segment_2& s2_)
// {
//   typedef Bso_kernel K;
//   typedef K::Segment_2 Bso_segment_2;

//   Bso_segment_2 s1 = change_kernel<K>(s1_);
//   Bso_segment_2 s2 = change_kernel<K>(s2_);
//   CGAL::Object o = CGAL::intersection(s1, s2);
//   const Bso_point_2* ip = CGAL::object_cast<Bso_point_2>(&o);
//   const Bso_segment_2* is = CGAL::object_cast<Bso_segment_2>(&o);
//   CGAL::Object ret;
//   if (ip)
//     ret = CGAL::make_object(to_common<K>(*ip, 0));
//   else if (is)
//     ret = CGAL::make_object(to_common<K>(*is, 0));
//   return ret;
// }

Bso_direction_2 direction(const Bso_point_2& from, const Bso_point_2& to)
{
  Bso_segment_2 seg(from, to);
  return Bso_direction_2(seg);
}

bool outer_is_max(const Bso_ray_2& ray, const Bso_point_2& p1, const Bso_point_2& p2, const Bso_point_2& p3, const Bso_direction_2& dmin)
{
  const Bso_point_2& p0 = ray.source();
  const Bso_direction_2& d0 = ray.direction();
  const Bso_direction_2 d1 = direction(p0, p1);  // d1 is prev
  const Bso_direction_2 d2 = direction(p0, p2);  // d2 is the cur point direction
  const Bso_direction_2 d3 = direction(p0, p3);  // d3 is next

  // return ((d1 == d0 || d2.counterclockwise_in_between(d1, d0)) && 
  return (d2.counterclockwise_in_between(d1, d0) && 
	  d2.counterclockwise_in_between(d3, d0) &&
	  (dmin == d0 || d2.counterclockwise_in_between(dmin, d0)) &&
	  CGAL::orientation(p1, p2, p3) == CGAL::RIGHT_TURN);
}

bool inner_is_max(const Bso_ray_2& ray, const Bso_point_2& p1, const Bso_point_2& p2, const Bso_point_2& p3, const Bso_direction_2& dmax)
{
  const Bso_point_2& p0 = ray.source();
  const Bso_direction_2& d0 = ray.direction();
  const Bso_direction_2 d1 = direction(p0, p1);  // d1 is prev
  const Bso_direction_2 d2 = direction(p0, p2);  // d2 is the cur point direction
  const Bso_direction_2 d3 = direction(p0, p3);  // d3 is next

   return (d1.counterclockwise_in_between(d2, d0) && 
	   d3.counterclockwise_in_between(d2, d0) &&
	   (dmax == d2 || dmax.counterclockwise_in_between(d2, d0)) &&
	   CGAL::orientation(p1, p2, p3) == CGAL::LEFT_TURN);
}

template <typename F, typename Out_iter>
void find_local_max(Circ start, Circ end, const Bso_point_2& p, const Bso_point_2& q, F is_max, Out_iter points)
{
  Circ cur = start;
  Circ next = cur+1;
  const Bso_direction_2 d0 = direction(p, q);
  const Bso_ray_2 ray(p, d0);
  Bso_point_2 p1 = *(cur-1);
  Bso_point_2 p2 = *cur;
  Bso_point_2 p3 = *next;
  Bso_direction_2 dmax = direction(p, p1);
  do {
    if (is_max(ray, p1, p2, p3, dmax)) {
      *points++ = cur;
      dmax = direction(p, *cur);
    }
    ++cur;
    ++next;
    p1 = p2;
    p2 = p3;
    p3 = *next;
  } while (cur != end);
}

class Hull_checker
{
public:
  // Hull_checker(const Bso_polygon_2& P) : _P(change_kernel<Bso_kernel>(P)) {}
  Hull_checker(const Bso_polygon_2& P) : _P(P) {}
  bool operator()(Circ p, Circ q) {
    // typedef Bso_kernel::Segment_2 Bso_segment_2;
    // Bso_segment_2 s = change_kernel<Bso_kernel>(Segment_2(*p, *q));
    Bso_segment_2 s(*p, *q);
    for (Bso_polygon_2::Edge_const_iterator it = _P.edges_begin(); it != _P.edges_end(); ++it) {
      if (CGAL::do_intersect(s, *it)) return false;
    }
    return true;
  }
private:
  Bso_polygon_2 _P;
};

// Returns the pair (constrained hull point, unconstrained hull point)
template <typename F>
pair<Circ, Circ> hull_point_outer(const Circ cur, const Circ prev, const Circ* stop, F constraints)
{
  static log4cplus::Logger logger = log4cplus::Logger::getInstance("decimate.hull_point_outer");

  // TODO: Take out stop stuff

  // This is troublesome because we may be able to get to a vertex that is beyond (prev,cur).
  // This is possible if the previous check was constrained.
  Bso_ray_2 r0(*prev, *cur);
  // if (CGAL::orientation(*prev, *cur, *(cur+1)) == CGAL::RIGHT_TURN) {
    r0 = Bso_ray_2(*cur, *prev);
  //   LOG4CPLUS_TRACE(logger, "  r = cur->prev (turn)");
  // }
  // else if (stop) {
  //   r0 = Bso_ray_2(*cur, *prev);
  //   LOG4CPLUS_TRACE(logger, "  r = cur->prev (stop)");
  // }
  // else
  //   LOG4CPLUS_TRACE(logger, "  r = prev->cur");

  if (stop)
    LOG4CPLUS_TRACE(logger, "  prev: " << pp(*prev) << " cur: " << pp(*cur) << " stop: " << pp(**stop));
  else
    LOG4CPLUS_TRACE(logger, "  prev: " << pp(*prev) << " cur: " << pp(*cur) << " stop: null");

  Bso_ray_2 urmax(*cur, *(cur+1)); // the unconstrained max
  Circ c = cur+2;
  Circ hullpt = cur+1;
  Circ uhullpt = cur+1;
  // true if we've failed a constraint check
  bool failure_mode = false;
  do {
    Bso_ray_2 r(*cur, *c);
    LOG4CPLUS_TRACE(logger, "  checking " << pp(*c));
    const bool between = (r.direction() == urmax.direction()) ||
      r.direction().counterclockwise_in_between(r0.direction(), urmax.direction());
    const bool passed = constraints(cur, c);
    if (between) {
      LOG4CPLUS_TRACE(logger, "    uhullpt updated");
      urmax = r;
      uhullpt = c;
      if (!failure_mode) {
	if (passed) {
	  Bso_polygon_2 tempp(cur, c+1);
	  if (tempp.is_simple() && tempp.is_clockwise_oriented()) {
	    LOG4CPLUS_TRACE(logger, "    hullpt updated");
	    hullpt = c;
	  }
	}
	else {
	  failure_mode = true;
	}
      }
    }
    ++c;
  // } while (c != cur && (!stop || c-1 != *stop));
  } while (c != cur && !failure_mode);//(!stop || c-1 != *stop));

  // if (stop && uhullpt != *stop) {
  //   LOG4CPLUS_TRACE(logger, "    updating uhullpt to stop");
  //   uhullpt = *stop;
  // }

  // return make_pair(hullpt, uhullpt);
  return make_pair(hullpt, hullpt);
}

template <typename F>
Circ hull_point_inner(Circ cur, Circ prev, Circ stop, F constraints)
{
  Bso_ray_2 r0(*cur, *prev);
  // if (CGAL::orientation(*prev, *cur, *(cur+1)) == CGAL::RIGHT_TURN)
    // r0 = Bso_ray_2(*prev, *cur);

  Bso_ray_2 rmin(*cur, *(cur+1));
  Circ c = cur+2;
  Circ hullpt = cur+1;
  do {
    Bso_ray_2 r(*cur, *c);
    const bool between = (r.direction() == rmin.direction()) ||
      r.direction().counterclockwise_in_between(rmin.direction(), r0.direction());
    if (constraints(cur, c)) {
      if (between) {
	rmin = r;
	hullpt = c;
      }
    }
    else {
      break;
    }
    ++c;
  } while (c != cur && c != stop);

  return hullpt;
}

// Find a polygon's constrained convex hull using a variant of Jarvis's march [1,2,3].
template <typename F>
Bso_polygon_2 convex_hull_outer(const Bso_polygon_2& polygon, F constraints)
{
  static log4cplus::Logger logger = log4cplus::Logger::getInstance("decimate.convex_hull_outer");

  Bso_polygon_2 hull;

  // Find an initial point on the convex hull
  Circ begin = polygon.vertices_circulator();
  Circ hullpt, uhullpt;
  boost::tie(hullpt, uhullpt) = hull_point_outer(begin, begin-1, 0, constraints);
  Circ first = hullpt;
  Circ prev = begin;
  // Now loop through finding all hull points
  LOG4CPLUS_TRACE(logger, "convex hull outer");
  do {
    hull.push_back(*hullpt);
    LOG4CPLUS_TRACE(logger, "hullpt: " << pp(*hullpt) << " first: " << pp(*first) << " begin: " << pp(*begin));
    Circ temp;
    if (hullpt != uhullpt) {
      boost::tie(temp, uhullpt) = hull_point_outer(hullpt, prev, &uhullpt, constraints);
    }
    else {
      boost::tie(temp, uhullpt) = hull_point_outer(hullpt, prev, 0, constraints);
    }
    prev = hullpt;
    hullpt = temp;
  } while (hullpt != first);

  return hull;
}

// Find a polygon's constrained inner convex hull using a variant of Jarvis's march [1,2,3].
template <typename F>
Bso_polygon_2 convex_hull_inner(const Bso_polygon_2& polygon, F constraints)
{
  static log4cplus::Logger logger = log4cplus::Logger::getInstance("decimate.convex_hull_inner");

  Bso_polygon_2 hull;

  // Find an initial point on the convex hull
  Circ begin = polygon.vertices_circulator();
  Circ hullpt = hull_point_inner(begin, begin-1, begin-1, constraints);
  Circ first = hullpt;
  Circ prev = begin;
  // Now loop through finding all hull points
  LOG4CPLUS_TRACE(logger, "convex hull inner");
  do {
    hull.push_back(*hullpt);
    LOG4CPLUS_TRACE(logger, "  hullpt: " << pp(*hullpt) << " first: " << pp(*first) << " begin: " << pp(*begin));
    Circ temp = hull_point_inner(hullpt, prev, first, constraints);
    prev = hullpt;
    hullpt = temp;
  } while (hullpt != first);

  return hull;
}

template <typename F, typename Out_iter>
void find_local_max(Circ start, const Bso_point_2& p, const Bso_point_2& q, F is_max, Out_iter points)
{
  find_local_max(start, start, p, q, is_max, points);
}

CGAL::Orientation orientation(Circ p)
{
  return CGAL::orientation(*(p-1), *p, *(p+1));
}

bool right_turn(Circ p) {
  return orientation(p) == CGAL::RIGHT_TURN;
}

bool left_turn(Circ p) {
  return orientation(p) == CGAL::LEFT_TURN;
}

Bso_point_2 intersection(const Bso_ray_2& r, Circ& p)
{
  static log4cplus::Logger logger = log4cplus::Logger::getInstance("decimate.intersection");

  // typedef Bso_kernel K;
  // typedef K::Ray_2 Bso_ray_2;
  // typedef K::Segment_2 Bso_segment_2;

  // Bso_ray_2 r = change_kernel<K>(r_);
  Bso_segment_2 s(*p, *(p+1));
  // Bso_segment_2 s = change_kernel<K>(Segment_2(*p, *(p+1)));
  // LOG4CPLUS_TRACE(logger, "Ray: " << pp(r.source()) << " " << pp(r.point(1)) << " Segment: " << pp(s));
  // CGAL::Object o = CGAL::intersection(r, s);
  CGAL::Object o = CGAL::intersection(r, s);
  // const Bso_point_2* ip = CGAL::object_cast<Bso_point_2>(&o);
  // const Bso_segment_2* is = CGAL::object_cast<Bso_segment_2>(&o);
  const Bso_point_2* ip = CGAL::object_cast<Bso_point_2>(&o);
  const Bso_segment_2* is = CGAL::object_cast<Bso_segment_2>(&o);
  Circ orig = p;
  while (!ip && !is) {
    p++;
    s = Bso_segment_2(*p, *(p+1));
    // s = change_kernel<K>(Segment_2(*p, *(p+1)));
    o = CGAL::intersection(r, s);
    // ip = CGAL::object_cast<Bso_point_2>(&o);
    // is = CGAL::object_cast<Bso_segment_2>(&o);
    ip = CGAL::object_cast<Bso_point_2>(&o);
    is = CGAL::object_cast<Bso_segment_2>(&o);
    if (p == orig) {
      throw logic_error("Found no intersections");
    }
  }
  // if (ip) return to_common<K>(*ip, 0);
  if (ip) return *ip;
  return *(p+1);
}

Bso_polygon_2 completed(const deque<Bso_point_2>& points)
{
  if (points.size() < 4) return Bso_polygon_2();
  typedef deque<Bso_point_2>::const_iterator Iter;
  typedef deque<Bso_point_2>::const_reverse_iterator rIter;
  Iter first = points.begin();
  Iter second = first; second++;
  rIter rfirst = points.rbegin();
  rIter rsecond = rfirst; rsecond++;
  Bso_segment_2 seg1(*first, *second);
  Bso_segment_2 seg2(*rfirst, *rsecond);
  if (!do_intersect(seg1, seg2)) {
    // Check first point
    if (CGAL::squared_distance(*first, seg2) < 0.000001) {
      deque<Bso_point_2> new_points(points.begin(), points.end());
      new_points.pop_back();
      return Bso_polygon_2(new_points.begin(), new_points.end());
    }
    return Bso_polygon_2();
  }
  
  CGAL::Object o = CGAL::intersection(seg1, seg2);
  const Bso_point_2* ip = CGAL::object_cast<Bso_point_2>(&o);
  const Bso_segment_2* is = CGAL::object_cast<Bso_segment_2>(&o);
  Bso_point_2 firstpt;
  if (ip) {
    firstpt = *ip;
  }
  else if (is) {
    firstpt = is->source();
  }
  else {
    throw logic_error("Should have had intersection");
  }
  deque<Bso_point_2> new_points(points.begin(), points.end());
  new_points.pop_front();
  new_points.pop_back();
  new_points.push_front(firstpt);
  return Bso_polygon_2(new_points.begin(), new_points.end());
}

bool local_max_inner(const Bso_point_2& p, Circ q)
{
  return CGAL::right_turn(*q, p, *(q+1));
}

bool local_max_outer(const Bso_point_2& p, Circ q)
{
  return CGAL::left_turn(*q, p, *(q+1));
}

typedef boost::shared_ptr<Circ> HCirc;

pair<HCirc, Bso_point_2> first_intersection(const Bso_ray_2& r, Circ begin, Circ end)
{
  for (Circ it = begin; it != end; ++it) {
    CGAL::Object o = CGAL::intersection(r, Bso_segment_2(*it, *(it+1)));
    const Bso_point_2* ip = CGAL::object_cast<Bso_point_2>(&o);
    const Bso_segment_2* is = CGAL::object_cast<Bso_segment_2>(&o);
    if (ip) {
      return make_pair(HCirc(new Circ(it)), *ip);
    }
    if (is) {
      return make_pair(HCirc(new Circ(it)), is->target());
    }
  }
  return make_pair(HCirc(), Bso_point_2());
}

Polygon_2 decimate_2(Polygon_2 inner_, Polygon_2 outer_)
{
  static log4cplus::Logger logger = log4cplus::Logger::getInstance("decimate.decimate_2");

  Bso_polygon_2 inner = change_kernel<Bso_kernel>(inner_);
  Bso_polygon_2 outer = change_kernel<Bso_kernel>(outer_);

  // Get the constrained convex hull of the inner polygon and update
  // the inner variable
  inner = convex_hull_outer(inner, Hull_checker(outer));
  gnuplot_print_polygon("output/hull1.dat", to_common<Bso_kernel>(inner, 0));

  outer = convex_hull_inner(outer, Hull_checker(inner));
  gnuplot_print_polygon("output/hull2.dat", to_common<Bso_kernel>(outer, 0));

  // Initialization
  Circ pcirc = inner.vertices_circulator();
  Circ pi = pcirc+1;

  // Find po
  Bso_ray_2 r(*(pcirc-1), *pcirc);
  Bso_direction_2 dir = r.direction();
  dir = dir.transform(CGAL::Aff_transformation_2<Bso_kernel>(CGAL::Rotation(), -1, 0));
  r = Bso_ray_2(*pcirc, dir);
  Circ cur = outer.vertices_circulator();
  Circ prev = cur - 1;
  Circ po;
  Bso_point_2 curp;
  bool initialized = false;
  LOG4CPLUS_TRACE(logger, "init po -- pcirc-1: " << pp(*(pcirc-1)) << " pcirc: " << pp(*pcirc));
  do {
    Bso_segment_2 s(*prev, *cur);
    LOG4CPLUS_TRACE(logger, "  segment: " << pp(s));
    CGAL::Object o = CGAL::intersection(r, s);
    const Bso_point_2* ip = CGAL::object_cast<Bso_point_2>(&o);
    const Bso_segment_2* is = CGAL::object_cast<Bso_segment_2>(&o);
    if (ip || is) {
      Bso_point_2 testp = *cur;
      if (ip) testp = *ip;
      LOG4CPLUS_TRACE(logger, "  found: " << pp(testp));
      if (!initialized || CGAL::has_smaller_distance_to_point(*pcirc, testp, curp)) {
	LOG4CPLUS_TRACE(logger, "  smaller");
  	po = cur;
  	curp = testp;
  	initialized = true;
      }
    } 
    prev = cur;
    cur++;
  } while (cur != outer.vertices_circulator());

  if (!initialized) {
    throw logic_error("Unable to find po");
  }

  // Now loop through
  deque<Bso_point_2> points;
  Bso_point_2 p = *pcirc;
  points.push_back(p);
  Bso_polygon_2 decimated;
  enum STEP { TO_UNKNOWN, TO_INSIDE, TO_OUTSIDE };
  STEP last_step = TO_INSIDE;
  int count = 0;
  while (decimated.size() == 0) {
    count++;
    LOG4CPLUS_TRACE(logger, "p: " << pp(p) << " (iteration " << count << ")");

    Circ maxi = pi;
    Circ maxo = po;
    STEP next_step = TO_UNKNOWN;
    while (!local_max_inner(p, maxi)) {
      maxi++;
    }
    // while (!local_max_outer(p, maxo)) {
    //   maxo++;
    //   if (left_turn(maxo)) {
    // 	// next_step = TO_OUTSIDE;
    // 	break;
    //   }
    // }
    maxo = right_turn(po) ? po : po+1;

    LOG4CPLUS_TRACE(logger, "  pi: " << pp(*pi) << " maxi: " << pp(*maxi));
    LOG4CPLUS_TRACE(logger, "  po: " << pp(*po) << " maxo: " << pp(*maxo));

    Bso_point_2 outer_int, inner_int;
    HCirc outer_int_vertex, inner_int_vertex;
    Bso_direction_2 outer_dir(Bso_ray_2(p, *maxi));
    // boost::tie(outer_int_vertex, outer_int) = first_intersection(Bso_ray_2(p, *maxi), po-1, po-2);
    boost::tie(outer_int_vertex, outer_int) = first_intersection(Bso_ray_2(*maxi, outer_dir), po-1, po-2);
    if (right_turn(maxo)) {
      if (last_step == TO_OUTSIDE) {
	LOG4CPLUS_TRACE(logger, "  pi-1 -> pi-2");
	boost::tie(inner_int_vertex, inner_int)  = first_intersection(Bso_ray_2(p, *maxo), pi-1, pi-2);
      }
      else {
	LOG4CPLUS_TRACE(logger, "  pi -> pi-1");
	boost::tie(inner_int_vertex, inner_int)  = first_intersection(Bso_ray_2(p, *maxo), pi, pi-1);
      }
    }

    // Logging
    if (logger.getLogLevel() == log4cplus::TRACE_LOG_LEVEL) {
      LOG4CPLUS_TRACE(logger, "  outer intersection: " << pp(outer_int) << " vertex: " << pp(**outer_int_vertex));
      if (inner_int_vertex)
	LOG4CPLUS_TRACE(logger, "  inner intersection: " << pp(inner_int) << " vertex: " << pp(**inner_int_vertex));

      string fn("output/temp" + boost::lexical_cast<string>(count) + ".dat");
      ofstream tout(fn.c_str());
      for (deque<Bso_point_2>::iterator it = points.begin(); it != points.end(); ++it)
	tout << to_common<Bso_kernel>(*it, 0) << endl;
      tout << endl << endl;
      tout.close();

      string infn("output/to-in" + boost::lexical_cast<string>(count) + ".dat");
      ofstream ins(infn.c_str());
      ins << to_common<Bso_kernel>(*po, 0) << endl << endl << endl;
      ins << to_common<Bso_kernel>(p, 0) << endl;
      ins << to_common<Bso_kernel>(inner_int, 0) << endl << endl << endl;
      // ins << to_common<Bso_kernel>(*maxo, 0) << endl << endl << endl;
      ins.close();

      string outfn("output/to-out" + boost::lexical_cast<string>(count) + ".dat");
      ofstream outs(outfn.c_str());
      outs << to_common<Bso_kernel>(*pi, 0) << endl << endl << endl;
      outs << to_common<Bso_kernel>(p, 0) << endl;
      outs << to_common<Bso_kernel>(outer_int, 0) << endl << endl << endl;
      // outs << to_common<Bso_kernel>(*maxi, 0) << endl << endl << endl;
      outs.close();
    }

    if (inner_int_vertex && has_larger_distance_to_point(p, inner_int, *maxo) && has_smaller_distance_to_point(p, inner_int, outer_int)) {
      LOG4CPLUS_TRACE(logger, "  to inner");
      last_step = TO_INSIDE;
      p = inner_int;
      po = maxo+1;
      pi = *inner_int_vertex + 1;
    }
    else {
      LOG4CPLUS_TRACE(logger, "  to outer");
      last_step = TO_OUTSIDE;
      p = outer_int;
      po = *outer_int_vertex + 1;
      pi = maxi+1;
    }

    // // p is on outside
    // if (last_step == TO_OUTSIDE) {
    //   if (right_turn(po) && right_turn(po+1)) {
    // 	next_step = TO_INSIDE;
    // 	maxo = po;
    //   }
    //   else if (right_turn(po) && left_turn(po+1)) {
    // 	next_step = right_turn(*po, p, *maxi) ? TO_OUTSIDE : TO_INSIDE;
    // 	maxo = po;
    //   }
    //   else if (left_turn(po+1)) {
    // 	next_step = TO_OUTSIDE;
    //   }
    //   else if (CGAL::do_intersect(Bso_ray_2(p, *maxi), Bso_segment_2(*po, *(po+1)))) {
    // 	next_step = TO_OUTSIDE;
    //   }
    //   else {
    // 	next_step = right_turn(*(po+1), p, *maxi) ? TO_OUTSIDE : TO_INSIDE;
    //   }
    // }
    // // p is on inside
    // else {
    //   if (left_turn(pi)) {
    // 	next_step = TO_OUTSIDE;
    //   }
    //   else if (left_turn(po)) {
    // 	next_step = TO_OUTSIDE;
    //   }
    //   else if (CGAL::do_intersect(Bso_ray_2(p, *maxi), Bso_segment_2(*(po-1), *po)) &&
    // 	       !CGAL::do_intersect(Bso_segment_2(p, *maxi), Bso_segment_2(*(po-1), *po))) {
    // 	next_step = TO_OUTSIDE;
    //   }
    //   else {
    // 	next_step = right_turn(*po, p, *maxi) ? TO_OUTSIDE : TO_INSIDE;
    //   }
    // }

    // if (next_step == TO_INSIDE && left_turn(maxo)) {
    //   maxi++;
    //   next_step = TO_OUTSIDE;
    // }

    // if (next_step == TO_UNKNOWN) {
    //   LOG4CPLUS_TRACE(logger, "  second check");
    //   // If we're inside, first try to go to the adjacent point.  If the adjacent
    //   // point is a right-hand turn, then we have to go to the outside.
    //   if (last_step == TO_INSIDE) {
    // 	LOG4CPLUS_TRACE(logger, "  last inside");
    // 	next_step = left_turn(pi) ? TO_OUTSIDE : TO_INSIDE;
    //   }
    //   // If we're outside, do checks
    //   else {
    // 	LOG4CPLUS_TRACE(logger, "  last outside");
    // 	next_step = (left_turn(maxo) || left_turn(*maxo, p, *maxi)) ? TO_OUTSIDE : TO_INSIDE;
    //   }
    // }


    // if (next_step == TO_OUTSIDE) {
    //   LOG4CPLUS_TRACE(logger, "  to outside");
    //   Circ temp = po;
    //   if (last_step == TO_INSIDE)
    //   	temp--;
    //   p = intersection(Bso_ray_2(p, *maxi), temp);
    //   // p = CGAL::midpoint(*maxi, p);
    //   po = temp;
    //   pi = maxi;
    // }
    // else {
    //   LOG4CPLUS_TRACE(logger, "  to inside");
    //   Circ temp = pi;
    //   if (last_step == TO_OUTSIDE)
    // 	temp--;
    //   // try {
    // 	p = intersection(Bso_ray_2(p, *maxo), temp);
    // 	// p = CGAL::midpoint(*maxo, p);
    // 	pi = temp;
    //   // }
    //   // catch (logic_error& e) {
    //   // 	next_step = TO_OUTSIDE;
    //   // 	po = maxo+1;
    //   // 	p = intersection(Bso_ray_2(p, *maxo), po);
    //   // 	// p = CGAL::midpoint(*maxo, p);
    //   // 	pi = maxi;
    //   // }
    //   po = maxo;
    // }
    // pi++;
    // po++;
    // last_step = next_step;
    points.push_back(p);

    decimated = completed(points);
  }

  return to_common<Bso_kernel>(decimated, 0);
}

CONTOURTILER_END_NAMESPACE

using namespace CONTOURTILER_NAMESPACE;

void arbitrary(int argc, char** argv)
{
  // Read inner
  list<Polygon_2> polygons;
  read_polygons_gnuplot2("test-data/decimate/p" + string(argv[1]) + ".dat", back_inserter(polygons), 0);
  Polygon_2 inner = *polygons.begin();
  gnuplot_print_polygon("output/orig1.dat", inner);

  // Read outer
  polygons.clear();
  read_polygons_gnuplot2("test-data/decimate/q" + string(argv[1]) + ".dat", back_inserter(polygons), 0);
  Polygon_2 outer = *polygons.begin();
  gnuplot_print_polygon("output/orig2.dat", outer);

  Polygon_2 p = decimate_2(inner, outer);
  gnuplot_print_polygon("output/out1.dat", p);
}

void single(int argc, char** argv)
{
  static log4cplus::Logger logger = log4cplus::Logger::getInstance("decimate.single");

  // Read polygon
  list<Polygon_2> polygons;
  read_polygons_gnuplot2("test-data/decimate/" + string(argv[1]), back_inserter(polygons), 0);
  const Polygon_2 p = *polygons.begin();
  gnuplot_print_polygon("output/orig.dat", p);

  double epsilon = boost::lexical_cast<double>(argv[2]);
  list<Polygon_with_holes_2> all_inner, all_outer;
  ContourTiler::offset_polygon_negative(p, -epsilon/2, back_inserter(all_inner));
  all_outer.push_back(ContourTiler::offset_polygon_positive(p, epsilon/2));
  Polygon_2 outer = all_outer.begin()->outer_boundary();
  gnuplot_print_polygon("output/outer.dat", outer);

  ofstream out_inner("output/inner.dat");
  ofstream out_out("output/out.dat");
  // list<Polygon_with_holes_2>::iterator end = all_inner.begin();
  // end++;
  // for (list<Polygon_with_holes_2>::iterator it = all_inner.begin(); it != end; ++it) {
    // const Polygon_2& inner = it->outer_boundary();
  const Polygon_2& inner = all_inner.begin()->outer_boundary();
    gnuplot_print_polygon(out_inner, inner);

    Polygon_2 q = decimate_2(inner, outer);
    gnuplot_print_polygon(out_out, q);

    out_inner << endl << endl;
    out_out << endl << endl;
  // }
  out_inner.close();
  out_out.close();

  LOG4CPLUS_INFO(logger, "Original polygon: " << p.size() << "  Decimated polygon: " << q.size());
}

int main(int argc, char** argv)
{
  log4cplus::PropertyConfigurator::doConfigure("log4cplus.properties");

  // arbitrary(argc, argv);
  single(argc, argv);

  return 0;
}

