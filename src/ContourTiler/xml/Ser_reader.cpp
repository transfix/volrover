#include <iostream>
#include <stdexcept>
#include <algorithm>
#include <vector>

#include <boost/tokenizer.hpp>
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/format.hpp>

#include <ContourTiler/xml/tinyxml.h>
#include <ContourTiler/xml/Ser_reader.h>

using namespace std;
using namespace boost; using boost::format;

//------------------------------------------------------------
// class Ser_transform
//------------------------------------------------------------
class Ser_transform
{
// coefficients are
//   1 x y xy xx yy
public:
  Ser_transform() {
    double xI[] = {0, 1, 0, 0, 0, 0};
    double yI[] = {0, 0, 1, 0, 0, 0};
    _xcoef = std::vector<double>(xI, xI+6);
    _ycoef = std::vector<double>(yI, yI+6);
  }
  Ser_transform(std::vector<double> xcoef, std::vector<double> ycoef) : _xcoef(xcoef), _ycoef(ycoef) {}
  Ser_point operator()(const Ser_point& point) const {
    double x = point.x();
    double y = point.y();
    x = -_xcoef[0] + x*_xcoef[1] + y*_xcoef[2] + x*y*_xcoef[3] + x*x*_xcoef[4] + y*y*_xcoef[5];
    y = -_ycoef[0] + x*_ycoef[1] + y*_ycoef[2] + x*y*_ycoef[3] + x*x*_ycoef[4] + y*y*_ycoef[5];
    return Ser_point(x,y);
  }
  bool is_I() const {
    bool I = true;
    for (int i = 0; i < 6; ++i) {
      if (i != 1 && _xcoef[i] != 0) {
        I = false;
      }
      if (i != 2 && _ycoef[i] != 0) {
        I = false;
      }
    }
    return I;
  }

  friend std::ostream& operator<<(std::ostream& out, const Ser_transform& transform) {
    out << "xcoef = [";
    for (int i = 0; i < 6; ++i) {
      out << transform._xcoef[i] << " ";
    }
    out << "] ycoef = [";
    for (int i = 0; i < 6; ++i) {
      out << transform._ycoef[i] << " ";
    }
    out << "]";
    return out;
  }

private:
  std::vector<double> _xcoef;
  std::vector<double> _ycoef;
};

//------------------------------------------------------------
// utilities
//------------------------------------------------------------

// trim from start
static inline std::string &ltrim(std::string &s) {
  s.erase(s.begin(), std::find_if(s.begin(), s.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
  return s;
}

// trim from end
static inline std::string &rtrim(std::string &s) {
  s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
  return s;
}

// trim from both ends
static inline std::string &trim(std::string &s) {
  return ltrim(rtrim(s));
}

//------------------------------------------------------------
// processing
//------------------------------------------------------------

Ser_transform process_transform(TiXmlElement* transform)
{
  if (!transform) throw Ser_exception("transform cannot be null");
  string xcoef_str, ycoef_str;

  if (transform->QueryValueAttribute("xcoef", &xcoef_str) != TIXML_SUCCESS) {
    throw Ser_exception("Failed to get xcoef");
  }

  if (transform->QueryValueAttribute("ycoef", &ycoef_str) != TIXML_SUCCESS) {
    throw Ser_exception("Failed to get ycoef");
  }

  // coefficients are
  //   1 x y xy xx yy
  char_separator<char> sep(" \t\n");
  tokenizer<char_separator<char> > xtok(xcoef_str, sep), ytok(ycoef_str, sep);
  vector<double> xcoef, ycoef;
  for (const auto& token : xtok) {
    xcoef.push_back(lexical_cast<double>(token));
  }
  for (const auto& token : ytok) {
    ycoef.push_back(lexical_cast<double>(token));
  }

  return Ser_transform(xcoef, ycoef);
}

Ser_contour process_contour(TiXmlElement* contour, const Ser_transform& transform)
{
  if (!contour) throw Ser_exception("contour cannot be null");
  string name, points_str;
  if (contour->QueryValueAttribute("name", &name) != TIXML_SUCCESS) {
    throw Ser_exception("Failed to get contour name");
  }

  if (contour->QueryValueAttribute("points", &points_str) != TIXML_SUCCESS) {
    throw Ser_exception("Failed to get contour points");
  }

  // Parse points
  typedef boost::char_separator<char> Separator;
  typedef tokenizer<Separator> Tokenizer;

  vector<Ser_point> points;
  Separator sep(" \t\n,");
  Tokenizer tok(points_str, sep);
  for (Tokenizer::iterator it = tok.begin(); 
       it != tok.end();
       ++it)
  {
    const double x = lexical_cast<double>(*it++);
    const double y = lexical_cast<double>(*it);
    points.push_back(transform(Ser_point(x,y)));
  }

  return Ser_contour(name, points.begin(), points.end());
}

Ser_section process_section(const string& filebase, int index, vector<string>& warnings, bool apply_transforms)
{
  string filename = filebase + "." + lexical_cast<string>(index);
  
  TiXmlDocument doc(filename);
  if (!doc.LoadFile()) {
    throw Ser_exception("Section " + lexical_cast<string>(index) 
                        + " doesn't exist for base " + filebase);
  }

  vector<Ser_contour> contours;

  TiXmlHandle docHandle(&doc);
  TiXmlElement* section = docHandle.FirstChild("Section").ToElement();
  if (!section) throw Ser_exception("No section node");
  string thickness_str;
  if (section->QueryValueAttribute("thickness", &thickness_str) != TIXML_SUCCESS) {
    warnings.push_back("Failed to get section thickness in section " + lexical_cast<string>(index));
    thickness_str = "0.05";
  }
  const double thickness = lexical_cast<double>(thickness_str);
  
  TiXmlElement* transform_xml = section->FirstChildElement("Transform");
  while (transform_xml) {
    Ser_transform transform;
    if (apply_transforms) {
      transform = process_transform(transform_xml);
    }
    TiXmlElement* contour_xml = transform_xml->FirstChildElement("Contour");
    while (contour_xml) {
      Ser_contour contour = process_contour(contour_xml, transform);
      if (!transform.is_I()) {
        stringstream ss;
        ss << "Contour " << contour.name() << " has non-identity transform in section " << index << " (" << transform << ")";
        // warnings.push_back(ss.str());
      }
      contours.push_back(contour);
      contour_xml = contour_xml->NextSiblingElement();
    }

    transform_xml = transform_xml->NextSiblingElement();
  }
  return Ser_section(index, thickness, contours.begin(), contours.end());
}

