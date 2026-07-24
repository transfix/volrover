#include <algorithm> //for min/max
#include <iterator>
#include <string>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <functional> //compose2, etc
#include <boost/lambda/lambda.hpp>
#include <boost/utility.hpp> //for prior()
#include <boost/format.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <ColorTable2/ColorTable.h>
#include <ColorTable2/Table.h>
#include <qlayout.h>

#if QT_VERSION < 0x040000
#include <ColorTable2/XoomedOut.h>
#endif

#ifndef COLORTABLE2_DISABLE_CONTOUR_TREE
#include <VolMagick/VolMagick.h>
#endif

namespace CVCColorTable
{
  template <typename T> static inline T clamp(T val)
  {
    return std::max(MIN_RANGE,std::min(MAX_RANGE,val));
  }

  ColorTable::ColorTable(QWidget *parent, 
#if QT_VERSION < 0x040000
                      const char *name
#else
                      Qt::WindowFlags flags
#endif
                      )
    : QFrame( parent, 
#if QT_VERSION < 0x040000
              name
#else
              flags
#endif
              )
  {
    //    setMaximumHeight(150);
    //    setMinimumHeight(80);
    setFrameStyle( QFrame::Panel | QFrame::Raised );
    QBoxLayout *layout = 
#if QT_VERSION < 0x040000
      new QBoxLayout( this, QBoxLayout::Down );
#else
    new QBoxLayout( QBoxLayout::Down, this );
#endif
    layout->setContentsMargins(3, 3, 3, 3);
    layout->setSpacing(3);

#if QT_VERSION < 0x040000
    _xoomedOut = new XoomedOut(info(),this,"XoomedOut");
#else
    _xoomedOut = new Table(Table::BACKGROUND |
                           Table::RANGE_BARS,
                           info(),this);
#endif

    layout->addWidget(_xoomedOut);
    _xoomedOut->setFixedHeight(20);

    _xoomedIn = new Table(Table::BACKGROUND |
                          Table::COLOR_BARS |
                          Table::ISOCONTOUR_BARS |
                          Table::OPACITY_NODES,
                          info(),this
#if QT_VERSION < 0x040000
                          ,"Table"
#endif
                          );
    layout->addWidget(_xoomedIn);

#if QT_VERSION < 0x040000
    connect(_xoomedOut,SIGNAL(minChanged(double)),_xoomedIn,SLOT(setMin(double)));
    connect(_xoomedOut,SIGNAL(minExploring(double)),_xoomedIn,SLOT(setMin(double)));
    connect(_xoomedOut,SIGNAL(maxChanged(double)),_xoomedIn,SLOT(setMax(double)));
    connect(_xoomedOut,SIGNAL(maxExploring(double)),_xoomedIn,SLOT(setMax(double)));
#else
    connect(_xoomedOut,SIGNAL(rangeMinChanged(double)),_xoomedIn,SLOT(setMin(double)));
    connect(_xoomedOut,SIGNAL(rangeMaxChanged(double)),_xoomedIn,SLOT(setMax(double)));
#endif

#if QT_VERSION < 0x040000
    connect(_xoomedIn,SIGNAL(changed()),_xoomedOut,SLOT(update()));
#else
    connect(_xoomedIn,SIGNAL(changed()),_xoomedOut,SLOT(update()));
#endif
    connect(_xoomedIn,SIGNAL(changed()),SIGNAL(changed()));

    info() = default_transfer_function();

    _opacityCubed = false;
  }

  ColorTable::~ColorTable()
  {
  }

  QSize ColorTable::sizeHint() const
  {
    return QSize(150, 150);
  }

  QSize ColorTable::minimumSizeHint() const
  {
    return QSize(150, 150);
  }

  boost::shared_array<unsigned char> ColorTable::getTable(unsigned int size) const
  {
    if(size==0) return boost::shared_array<unsigned char>();
    boost::shared_array<unsigned char> table(new unsigned char[size*4]);

    //if for some reason we don't have enough nodes, return a simple grayscale ramp
    if(info().colorNodes().size() < 2 || info().opacityNodes().size() < 2)
      {
	for(unsigned int i = 0; i < size; i++)
	  for(unsigned int j = 0; j < 4; j++)
	    table[i*4+j] = static_cast<unsigned char>((double(i)/double(size-1))*255.0);
	return table;
      }

    //TODO: i'd like to do this on the graphics hardware for speed,
    //but lets just do a software implementation for now

    //    double inc = (MAX_RANGE-MIN_RANGE)/(size-1 > 0 ? size-1 : 1);
    //    for(double pos = MIN_RANGE; pos <= MAX_RANGE; pos += inc)
#if 0
    for(unsigned int i = 0; i < size; i++)
      {
	double pos = MIN_RANGE + (MAX_RANGE-MIN_RANGE)*(double(i)/(double(size)-1 > 0 ? double(size)-1 : 1));
	color_nodes::const_iterator high_itr = info().colorNodes().lower_bound(color_node(pos));
	if(high_itr == info().colorNodes().end()) //we're past the defined nodes! just return what we have 
	  return table;
	color_nodes::const_iterator low_itr = 
	  high_itr == info().colorNodes().begin() ? high_itr : boost::prior(high_itr);
	opacity_nodes::const_iterator high_opac_itr = info().opacityNodes().lower_bound(opacity_node(pos));
	if(high_opac_itr == info().opacityNodes().end()) //past defined nodes!
	  return table;
	opacity_nodes::const_iterator low_opac_itr =
	  high_opac_itr == info().opacityNodes().begin() ? high_opac_itr : boost::prior(high_opac_itr);

	color_node high = *high_itr;
	color_node low = *low_itr;
	opacity_node high_opac = *high_opac_itr;
	opacity_node low_opac = *low_opac_itr;
	
	double interval_pos = (pos - low.position)/(high.position - low.position);
	qDebug("i == %d",i);
	qDebug("low.position = %f, high.position = %f",low.position,high.position);
	qDebug("pos = %f,interval_pos = %f",pos,interval_pos);
	
	table[i*4 + 0] = static_cast<unsigned char>(clamp(low.r + (high.r - low.r)*interval_pos)*255.0);
	table[i*4 + 1] = static_cast<unsigned char>(clamp(low.g + (high.g - low.g)*interval_pos)*255.0);
	table[i*4 + 2] = static_cast<unsigned char>(clamp(low.b + (high.b - low.b)*interval_pos)*255.0);
	table[i*4 + 3] = 
	  static_cast<unsigned char>(clamp(low_opac.value + (high_opac.value - low_opac.value)*interval_pos)*255.0);
      }
#endif

    for(color_nodes::const_iterator cur = info().colorNodes().begin();
	cur != boost::prior(info().colorNodes().end());
	cur++)
      {
	color_nodes::const_iterator next = boost::next(cur);
	unsigned int cur_idx = static_cast<unsigned int>(((cur->position - MIN_RANGE)/(MAX_RANGE - MIN_RANGE))*255.0);
	//unsigned int next_idx = cur_idx + ((next->position - cur->position)/(MAX_RANGE - MIN_RANGE))*255.0;
	unsigned int next_idx = static_cast<unsigned int>(((next->position - MIN_RANGE)/(MAX_RANGE - MIN_RANGE))*255.0);
	for(unsigned int i = cur_idx; i <= next_idx; i++)
	  {
	    double u = (double(i-cur_idx)/double(next_idx-cur_idx));
	    table[i*4 + 0] = static_cast<unsigned char>(clamp(cur->r + (next->r - cur->r)*u)*255.0);
	    table[i*4 + 1] = static_cast<unsigned char>(clamp(cur->g + (next->g - cur->g)*u)*255.0);
	    table[i*4 + 2] = static_cast<unsigned char>(clamp(cur->b + (next->b - cur->b)*u)*255.0);
	  }
      }

    unsigned int last_idx = -1;
    for(opacity_nodes::const_iterator cur = info().opacityNodes().begin();
	cur != boost::prior(info().opacityNodes().end());
	cur++)
      {
	opacity_nodes::const_iterator next = boost::next(cur);
	unsigned int cur_idx = static_cast<unsigned int>(((cur->position - MIN_RANGE)/(MAX_RANGE - MIN_RANGE))*255.0);
	unsigned int next_idx = cur_idx + 
	  static_cast<unsigned int>(((next->position - cur->position)/(MAX_RANGE - MIN_RANGE))*255.0);
	for(unsigned int i = cur_idx; i <= next_idx; i++)
	  {
	    double u = (double(i-cur_idx)/double(next_idx-cur_idx));

	    double v = clamp(cur->value + (next->value - cur->value)*u);
	    // arand: trick to make opacity work better, 4-11-2011	    
            if(_opacityCubed) v = v*v*v;
	    table[i*4 + 3] = static_cast<unsigned char>(v*255.0);
	  }
	last_idx = next_idx;
      }

  
  // arand, 6-8-2011
  // fix up the boundary cases which were causing a long standing bug
  color_nodes::const_iterator cur = info().colorNodes().begin();
  table[0] = static_cast<unsigned char>(clamp(cur->r)*255.0);
  table[1] = static_cast<unsigned char>(clamp(cur->g)*255.0);
  table[2] = static_cast<unsigned char>(clamp(cur->b)*255.0);

  cur = info().colorNodes().end();
  cur--;
  table[4*(size-1)+0] = static_cast<unsigned char>(clamp(cur->r)*255.0);
  table[4*(size-1)+1] = static_cast<unsigned char>(clamp(cur->g)*255.0);
  table[4*(size-1)+2] = static_cast<unsigned char>(clamp(cur->b)*255.0);

  opacity_nodes::const_iterator cur1 = info().opacityNodes().begin();
  double v = clamp(cur1->value);
  if(_opacityCubed) v = v*v*v;
  table[3] = static_cast<unsigned char>(v*255.0);

  cur1 = info().opacityNodes().end();
  cur1--;
  v = clamp(cur1->value);
  if(_opacityCubed) v = v*v*v;	  
  table[4*(size-1)+3] = static_cast<unsigned char>(v*255.0);  
 
    return table;
  }


  unsigned char* ColorTable::getCharTable(unsigned int size) const
  {
    if(size==0) return NULL;
    unsigned char* table = new unsigned char[size * 4];

    //if for some reason we don't have enough nodes, return a simple grayscale ramp
    if(info().colorNodes().size() < 2 || info().opacityNodes().size() < 2)
      {
	for(unsigned int i = 0; i < size; i++)
	  for(unsigned int j = 0; j < 4; j++)
	    table[i*4+j] = static_cast<unsigned char>((double(i)/double(size-1))*255.0);
	return table;
      }

    for(color_nodes::const_iterator cur = info().colorNodes().begin();
	cur != boost::prior(info().colorNodes().end());
	cur++)
      {
	color_nodes::const_iterator next = boost::next(cur);
	unsigned int cur_idx = static_cast<unsigned int>(((cur->position - MIN_RANGE)/(MAX_RANGE - MIN_RANGE))*255.0);
	//unsigned int next_idx = cur_idx + ((next->position - cur->position)/(MAX_RANGE - MIN_RANGE))*255.0;
	unsigned int next_idx = static_cast<unsigned int>(((next->position - MIN_RANGE)/(MAX_RANGE - MIN_RANGE))*255.0);
	for(unsigned int i = cur_idx; i <= next_idx; i++)
	  {
	    double u = (double(i-cur_idx)/double(next_idx-cur_idx));
	    table[i*4 + 0] = static_cast<unsigned char>(clamp(cur->r + (next->r - cur->r)*u)*255.0);
	    table[i*4 + 1] = static_cast<unsigned char>(clamp(cur->g + (next->g - cur->g)*u)*255.0);
	    table[i*4 + 2] = static_cast<unsigned char>(clamp(cur->b + (next->b - cur->b)*u)*255.0);
	  }
      }

    for(opacity_nodes::const_iterator cur = info().opacityNodes().begin();
	cur != boost::prior(info().opacityNodes().end());
	cur++)
      {
	opacity_nodes::const_iterator next = boost::next(cur);
	unsigned int cur_idx = static_cast<unsigned int>(((cur->position - MIN_RANGE)/(MAX_RANGE - MIN_RANGE))*255.0);
	unsigned int next_idx = cur_idx + 
	  static_cast<unsigned int>(((next->position - cur->position)/(MAX_RANGE - MIN_RANGE))*255.0);
	for(unsigned int i = cur_idx; i <= next_idx; i++)
	  {
	    double u = (double(i-cur_idx)/double(next_idx-cur_idx));

	    // arand: trick to make opacity work better, 4-11-2011
            if(_opacityCubed)
              {
                double v = clamp(cur->value + (next->value - cur->value)*u);
                v = v*v*v;
                table[i*4 + 3] = static_cast<unsigned char>(v*255.0);
              }
            else
              table[i*4 + 3] = static_cast<unsigned char>(clamp(cur->value + (next->value - cur->value)*u)*255.0);
	  }
      }

  // arand, 6-8-2011
  // fix up the boundary cases which were causing a long standing bug
  color_nodes::const_iterator cur = info().colorNodes().begin();
  table[0] = static_cast<unsigned char>(clamp(cur->r)*255.0);
  table[1] = static_cast<unsigned char>(clamp(cur->g)*255.0);
  table[2] = static_cast<unsigned char>(clamp(cur->b)*255.0);

  cur = info().colorNodes().end();
  cur--;
  table[4*(size-1)+0] = static_cast<unsigned char>(clamp(cur->r)*255.0);
  table[4*(size-1)+1] = static_cast<unsigned char>(clamp(cur->g)*255.0);
  table[4*(size-1)+2] = static_cast<unsigned char>(clamp(cur->b)*255.0);

  opacity_nodes::const_iterator cur1 = info().opacityNodes().begin();
  double v = clamp(cur1->value);
  if(_opacityCubed) v = v*v*v;	  
  table[3] = static_cast<unsigned char>(v*255.0);

  cur1 = info().opacityNodes().end();
  cur1--;
  v = clamp(cur1->value);
  if(_opacityCubed) v = v*v*v;	  
  table[4*(size-1)+3] = static_cast<unsigned char>(v*255.0);  


    return table;
  }

  // arand, 10-7-2011, new Float table for eventually better Volume Rendering
  boost::shared_array<float> ColorTable::getFloatTable(unsigned int size) const
  {
    if(size==0) return boost::shared_array<float>();
    boost::shared_array<float> table(new float[size*4]);

    //if for some reason we don't have enough nodes, return a simple grayscale ramp
    if(info().colorNodes().size() < 2 || info().opacityNodes().size() < 2)
      {
	for(unsigned int i = 0; i < size; i++)
	  for(unsigned int j = 0; j < 4; j++)
	    table[i*4+j] =float(i)/float(size-1);
	return table;
      }

    //TODO: i'd like to do this on the graphics hardware for speed,
    //but lets just do a software implementation for now

    //    double inc = (MAX_RANGE-MIN_RANGE)/(size-1 > 0 ? size-1 : 1);
    //    for(double pos = MIN_RANGE; pos <= MAX_RANGE; pos += inc)
#if 0
    for(unsigned int i = 0; i < size; i++)
      {
	double pos = MIN_RANGE + (MAX_RANGE-MIN_RANGE)*(double(i)/(double(size)-1 > 0 ? double(size)-1 : 1));
	color_nodes::const_iterator high_itr = info().colorNodes().lower_bound(color_node(pos));
	if(high_itr == info().colorNodes().end()) //we're past the defined nodes! just return what we have 
	  return table;
	color_nodes::const_iterator low_itr = 
	  high_itr == info().colorNodes().begin() ? high_itr : boost::prior(high_itr);
	opacity_nodes::const_iterator high_opac_itr = info().opacityNodes().lower_bound(opacity_node(pos));
	if(high_opac_itr == info().opacityNodes().end()) //past defined nodes!
	  return table;
	opacity_nodes::const_iterator low_opac_itr =
	  high_opac_itr == info().opacityNodes().begin() ? high_opac_itr : boost::prior(high_opac_itr);

	color_node high = *high_itr;
	color_node low = *low_itr;
	opacity_node high_opac = *high_opac_itr;
	opacity_node low_opac = *low_opac_itr;
	
	double interval_pos = (pos - low.position)/(high.position - low.position);
	qDebug("i == %d",i);
	qDebug("low.position = %f, high.position = %f",low.position,high.position);
	qDebug("pos = %f,interval_pos = %f",pos,interval_pos);
	
	table[i*4 + 0] = clamp(low.r + (high.r - low.r)*interval_pos);
	table[i*4 + 1] = clamp(low.g + (high.g - low.g)*interval_pos);
	table[i*4 + 2] = clamp(low.b + (high.b - low.b)*interval_pos);
	table[i*4 + 3] = clamp(low_opac.value + (high_opac.value - low_opac.value)*interval_pos);
      }
#endif


    // arand: BUG warning.. I think the table size of 256 may be hardcoded below...
    // the general case is above?
    for(color_nodes::const_iterator cur = info().colorNodes().begin();
	cur != boost::prior(info().colorNodes().end());
	cur++)
      {
	color_nodes::const_iterator next = boost::next(cur);
	unsigned int cur_idx = static_cast<unsigned int>(((cur->position - MIN_RANGE)/(MAX_RANGE - MIN_RANGE))*255.0);

	unsigned int next_idx = static_cast<unsigned int>(((next->position - MIN_RANGE)/(MAX_RANGE - MIN_RANGE))*255.0);
	for(unsigned int i = cur_idx; i <= next_idx; i++)
	  {
	    double u = (double(i-cur_idx)/double(next_idx-cur_idx));
	    table[i*4 + 0] = clamp(cur->r + (next->r - cur->r)*u);
	    table[i*4 + 1] = clamp(cur->g + (next->g - cur->g)*u);
	    table[i*4 + 2] = clamp(cur->b + (next->b - cur->b)*u);
	  }
      }

    unsigned int last_idx = -1;
    for(opacity_nodes::const_iterator cur = info().opacityNodes().begin();
	cur != boost::prior(info().opacityNodes().end());
	cur++)
      {
	opacity_nodes::const_iterator next = boost::next(cur);
	unsigned int cur_idx = static_cast<unsigned int>(((cur->position - MIN_RANGE)/(MAX_RANGE - MIN_RANGE))*255.0);
	unsigned int next_idx = cur_idx + 
	  static_cast<unsigned int>(((next->position - cur->position)/(MAX_RANGE - MIN_RANGE))*255.0);
	for(unsigned int i = cur_idx; i <= next_idx; i++)
	  {
	    double u = (double(i-cur_idx)/double(next_idx-cur_idx));

	    float v = clamp(cur->value + (next->value - cur->value)*u);
	    // arand: trick to make opacity work better, 4-11-2011	    
            if(_opacityCubed) v = v*v*v;
	    table[i*4 + 3] = v;
	  }
	last_idx = next_idx;
      }

  color_nodes::const_iterator cur = info().colorNodes().begin();
  table[0] = clamp(cur->r);
  table[1] = clamp(cur->g);
  table[2] = clamp(cur->b);

  cur = info().colorNodes().end();
  cur--;
  table[4*(size-1)+0] = clamp(cur->r);
  table[4*(size-1)+1] = clamp(cur->g);
  table[4*(size-1)+2] = clamp(cur->b);

  opacity_nodes::const_iterator cur1 = info().opacityNodes().begin();
  float v = clamp(cur1->value);
  if(_opacityCubed) v = v*v*v;
  table[3] = v;

  cur1 = info().opacityNodes().end();
  cur1--;
  v = clamp(cur1->value);
  if(_opacityCubed) v = v*v*v;	  
  table[4*(size-1)+3] = v;
  
  return table;
  }








  bool ColorTable::interactiveUpdates() const
  {
    return  _xoomedIn->interactiveUpdates();
  }

  bool ColorTable::opacityCubed() const
  {
    return _opacityCubed;
  }

  ColorTable::color_table_info& ColorTable::color_table_info::normalize()
  {
    using namespace std;
    //using namespace boost::lambda;
    //using boost::lambda::_1;

    {
      color_nodes::iterator bound;
      bound = colorNodes().lower_bound(color_node(MIN_RANGE));
      if(bound == colorNodes().end())
	colorNodes().insert(color_node(MIN_RANGE,0.0,0.0,0.0));
      else if(bound->position > MIN_RANGE) //just duplicate the node closest to MIN_RANGE
	colorNodes().insert(color_node(MIN_RANGE,bound->r,bound->g,bound->b));
      
      bound = colorNodes().lower_bound(color_node(MAX_RANGE));
      if(bound == colorNodes().end())
	colorNodes().insert(color_node(MAX_RANGE,1.0,1.0,1.0));
      else if(bound->position > MAX_RANGE)
	colorNodes().insert(color_node(MAX_RANGE,bound->r,bound->g,bound->b));
    }

    {
      opacity_nodes::iterator bound;
      bound = opacityNodes().lower_bound(opacity_node(MIN_RANGE));
      if(bound == opacityNodes().end())
	opacityNodes().insert(opacity_node(MIN_RANGE,0.0));
      else if(bound->position > MIN_RANGE) //just duplicate the node closest to MIN_RANGE
	opacityNodes().insert(opacity_node(MIN_RANGE,bound->value));
      
      bound = opacityNodes().lower_bound(opacity_node(MAX_RANGE));
      if(bound == opacityNodes().end())
	opacityNodes().insert(opacity_node(MAX_RANGE,1.0));
      else if(bound->position > MAX_RANGE)
	opacityNodes().insert(opacity_node(MAX_RANGE,bound->value));
    }


    // arand, 2/20 added boost::lambda namespace below to fix a compile issue
    //               maybe a boost version changed to cause the
    //               using boost::lambda above to break?

    //Now that we're sure we have a range of nodes within [MIN_RANGE,MAX_RANGE],
    //lets clear out everything not in that range.
    {
      color_nodes pruned;
      remove_copy_if(colorNodes().begin(),colorNodes().end(),
		     insert_iterator<color_nodes>(pruned,pruned.begin()),
		     (boost::lambda::_1 < color_node(MIN_RANGE)) || (color_node(MAX_RANGE) < boost::lambda::_1));
      colorNodes() = pruned;
    }
    
    {
      opacity_nodes pruned;
      remove_copy_if(opacityNodes().begin(),opacityNodes().end(),
		     insert_iterator<opacity_nodes>(pruned,pruned.begin()),
		     (boost::lambda::_1 < opacity_node(MIN_RANGE)) || (opacity_node(MAX_RANGE) < boost::lambda::_1));
      opacityNodes() = pruned;
    }

    return *this;
  }

  void ColorTable::update()
  {
    _xoomedIn->update();
    _xoomedOut->update();
  }

  void ColorTable::reset()
  {
    info() = default_transfer_function();
    update();
  }

  void ColorTable::read(const std::string& filename)
  {
    info() = read_transfer_function(filename);
    update();
  }

  void ColorTable::write(const std::string& filename)
  {
       write_transfer_function(filename,info());
  }

#if !defined(COLORTABLE2_DISABLE_CONTOUR_TREE) || !defined(COLORTABLE2_DISABLE_CONTOUR_SPECTRUM)
  void ColorTable::setContourVolume(const VolMagick::Volume& vol)
  {
    _xoomedIn->setContourVolume(vol);
  }
#endif

  ColorTable::color_table_info ColorTable::default_transfer_function()
  {
    color_table_info cti;


   //*
    cti.colorNodes().insert(color_node(MIN_RANGE,0.0,0.0,0.0));
    cti.colorNodes().insert(color_node(MIN_RANGE + 0.5*(MAX_RANGE - MIN_RANGE),
				       1.0,0.0,0.0));
    cti.colorNodes().insert(color_node(MIN_RANGE + 0.75*(MAX_RANGE - MIN_RANGE),
				       1.0,1.0,0.0));
    cti.colorNodes().insert(color_node(MAX_RANGE,0.0,1.0,0.0));

    cti.opacityNodes().insert(opacity_node(MIN_RANGE,0.0));
    cti.opacityNodes().insert(opacity_node(MIN_RANGE + 0.25*(MAX_RANGE - MIN_RANGE),
					      0.75));
    cti.opacityNodes().insert(opacity_node(MIN_RANGE + 0.75*(MAX_RANGE - MIN_RANGE),
					      0.25));
    cti.opacityNodes().insert(opacity_node(MAX_RANGE,1.0));
    /*/

    // arand: new default which looks good with white and black backgrounds...
    cti.colorNodes().insert(color_node(MIN_RANGE,1.0,0.0,0.0));
    cti.colorNodes().insert(color_node(MIN_RANGE + 0.5*(MAX_RANGE - MIN_RANGE),
				       1.0,1.0,0.0));
    cti.colorNodes().insert(color_node(MAX_RANGE,0.0,0.0,1.0));    
    
    cti.opacityNodes().insert(opacity_node(MIN_RANGE,0.0));
    cti.opacityNodes().insert(opacity_node(MIN_RANGE + 0.25*(MAX_RANGE - MIN_RANGE),
					      0.5));
    cti.opacityNodes().insert(opacity_node(MIN_RANGE + 0.5*(MAX_RANGE - MIN_RANGE),
					      0.12));
    cti.opacityNodes().insert(opacity_node(MAX_RANGE,0.75)); 
	*/

    return cti;
  }

  std::string ColorTable::transfer_function_filename = "";

  ColorTable::color_table_info ColorTable::read_transfer_function(const std::string& filename)
  {
    using namespace std;
    using namespace boost; using boost::format;

    color_table_info cti;

    ifstream inf(filename.c_str());
    if(!inf)
      throw runtime_error(string("Could not open ") + filename);

    unsigned int line_num = 0;
    string line;
    vector<string> split_line;

#define CHECK_LINE(check_str)						\
    {									\
      getline(inf, line); line_num++;					\
      if(!inf)								\
	std::cout << str(format("Error reading file %1%, line %2%")	\
				% filename				\
		    % line_num) << std::endl;				\
      if(line != check_str)						\
	std::cout << str(format("Error reading file %1%, line %2%: "	\
				       "string '%3%' not found")	\
				% filename				\
				% line_num				\
			 % string(check_str)) << endl;			\
    }

    /* // arand: old version crashed when a bad vinay is loaded
      if(!inf)								\
	throw runtime_error(str(format("Error reading file %1%, line %2%") \
				% filename				\
				% line_num));				\
      if(line != check_str)						\
	throw runtime_error(str(format("Error reading file %1%, line %2%: " \
				       "string '%3%' not found")	\
				% filename				\
				% line_num				\
				% string(check_str)));			\
    */

    CHECK_LINE("Anthony and Vinay are Great.");

    while(!getline(inf,line).eof())
      {
	line_num++;
	if(!inf)
	  throw runtime_error(str(format("Error reading file %1%, line %2%")
				  % filename
				  % line_num));
	
	if(line == "Alphamap") //process alphamap
	  {
	    CHECK_LINE("Number of nodes");
	    
	    getline(inf, line); line_num++;
	    if(!inf)
	      throw runtime_error(str(format("Error reading file %1%, line %2%")
				      % filename
				      % line_num));
	    unsigned int num_nodes = lexical_cast<unsigned int>(line);
	    
	    CHECK_LINE("Position and opacity");

	    for(unsigned int i = 0; i < num_nodes; i++)
	      {
		getline(inf, line); line_num++;
		if(!inf)
		  throw runtime_error(str(format("Error reading file %1%, line %2%")
					  % filename
					  % line_num));
		split(split_line,line,is_any_of(" "));
		if(split_line.size() != 2)
		  throw runtime_error(str(format("Error reading file %1%, line %2%: "
						 "Invalid position and opacity")
					  % filename
					  % line_num));
		cti.opacityNodes().insert(ColorTable::opacity_node(lexical_cast<double>(split_line[0]),
								   lexical_cast<double>(split_line[1])));
	      }
	  }
	else if(line == "ColorMap")
	  {
	    CHECK_LINE("Number of nodes");
	    
	    getline(inf, line); line_num++;
	    if(!inf)
	      throw runtime_error(str(format("Error reading file %1%, line %2%")
				      % filename
				      % line_num));
	    unsigned int num_nodes = lexical_cast<unsigned int>(line);
	    
	    CHECK_LINE("Position and RGB");
	    
	    for(unsigned int i = 0; i < num_nodes; i++)
	      {
		getline(inf, line); line_num++;
		if(!inf)
		  throw runtime_error(str(format("Error reading file %1%, line %2%")
					  % filename
					  % line_num));
		split(split_line,line,is_any_of(" "));
		if(split_line.size() != 4)
		  throw runtime_error(str(format("Error reading file %1%, line %2%: "
						 "Invalid position and RGB values")
					  % filename
					  % line_num));
		cti.colorNodes().insert(ColorTable::color_node(lexical_cast<double>(split_line[0]),
							       lexical_cast<double>(split_line[1]),
							       lexical_cast<double>(split_line[2]),
							       lexical_cast<double>(split_line[3])));
	      }
	  }
	else if(line == "IsocontourMap")
	  {
	    CHECK_LINE("Number of nodes");
	    
	    getline(inf, line); line_num++;
	    if(!inf)
	      throw runtime_error(str(format("Error reading file %1%, line %2%")
				      % filename
				      % line_num));
	    unsigned int num_nodes = lexical_cast<unsigned int>(line);
	    
	    CHECK_LINE("Position");
	    
	    for(unsigned int i = 0; i < num_nodes; i++)
	      {
		getline(inf, line); line_num++;
		if(!inf)
		  throw runtime_error(str(format("Error reading file %1%, line %2%")
					  % filename
					  % line_num));
		cti.isocontourNodes().insert(ColorTable::isocontour_node(lexical_cast<double>(line)));
	      }
	  }
      }

#undef CHECK_LINE

    transfer_function_filename = filename;

    return cti.normalize();
  }

  void ColorTable::write_transfer_function(const std::string& filename,
					   const ColorTable::color_table_info& cti)
  {
    using namespace std;
    using namespace boost; using boost::format;

    color_table_info local_cti = cti;
    local_cti.normalize();

    ofstream outf(filename.c_str());
    if(!outf)
      throw runtime_error(string("Could not open ") + filename);

    outf << "Anthony and Vinay are Great." << endl;

    outf << "Alphamap" << endl;
    outf << "Number of nodes" << endl;
    outf << local_cti.opacityNodes().size() << endl;
    outf << "Position and opacity" << endl;
    //write out the beginning and ending first
    outf << local_cti.opacityNodes().begin()->position << " "
	 << local_cti.opacityNodes().begin()->value << endl;
    outf << prior(local_cti.opacityNodes().end())->position << " "
	 << prior(local_cti.opacityNodes().end())->value << endl;
    if(local_cti.opacityNodes().size() > 2)
      for(opacity_nodes::iterator i = std::next(local_cti.opacityNodes().begin());
	  i != prior(local_cti.opacityNodes().end());
	  i++)
	outf << i->position << " " << i->value << endl;
    
    outf << "ColorMap" << endl;
    outf << "Number of nodes" << endl;
    outf << local_cti.colorNodes().size() << endl;
    outf << "Position and RGB" << endl;
    //write out the beginning and ending first
    outf << local_cti.colorNodes().begin()->position << " "
	 << local_cti.colorNodes().begin()->r << " "
	 << local_cti.colorNodes().begin()->g << " "
	 << local_cti.colorNodes().begin()->b << endl;
    outf << prior(local_cti.colorNodes().end())->position << " "
	 << prior(local_cti.colorNodes().end())->r << " "
	 << prior(local_cti.colorNodes().end())->g << " "
	 << prior(local_cti.colorNodes().end())->b << endl;
    if(local_cti.colorNodes().size() > 2)
      for(color_nodes::iterator i = std::next(local_cti.colorNodes().begin());
	  i != prior(local_cti.colorNodes().end());
	  i++)
	outf << i->position << " "
	     << i->r << " "
	     << i->g << " "
	     << i->b << endl;
    
    outf << "IsocontourMap" << endl;
    outf << "Number of nodes" << endl;
    outf << local_cti.isocontourNodes().size() << endl;
    outf << "Position" << endl;
    for(isocontour_nodes::iterator i = local_cti.isocontourNodes().begin();
	i != local_cti.isocontourNodes().end();
	i++)
      outf << i->position << endl;
  }

  void ColorTable::write_full_color_table( const std::string& filename, const ColorTable::color_table_info& cti )
  {
    color_table_info local_cti = cti;

    FILE *ofp = fopen( filename.c_str(), "w");
    if( !ofp ) {
      std::cout << "Could not open " << filename << std::endl;
      return;
    }
    
    /*  
    ofstream outf(filename.c_str());
    if(!outf)
      throw runtime_error(string("Could not open ") + filename);
    */

     unsigned int size = 256;
     //=======================================================
     unsigned char* table = new unsigned char[size * 4];

     //if for some reason we don't have enough nodes, return a simple grayscale ramp
     if(local_cti.colorNodes().size() < 2 || local_cti.opacityNodes().size() < 2)
      {
	for(unsigned int i = 0; i < size; i++)
	  for(unsigned int j = 0; j < 4; j++)
	    table[i*4+j] = static_cast<unsigned char>((double(i)/double(size-1))*255.0);
      }
     else
      {
         for(color_nodes::const_iterator cur = local_cti.colorNodes().begin();
	    cur != boost::prior(local_cti.colorNodes().end());
	    cur++)
          {
	    color_nodes::const_iterator next = boost::next(cur);
	    unsigned int cur_idx = static_cast<unsigned int>(((cur->position - MIN_RANGE)/(MAX_RANGE - MIN_RANGE))*255.0);
	    //unsigned int next_idx = cur_idx + ((next->position - cur->position)/(MAX_RANGE - MIN_RANGE))*255.0;
	    unsigned int next_idx = static_cast<unsigned int>(((next->position - MIN_RANGE)/(MAX_RANGE - MIN_RANGE))*255.0);
	    for(unsigned int i = cur_idx; i <= next_idx; i++)
	      {
	        double u = (double(i-cur_idx)/double(next_idx-cur_idx));
	        table[i*4 + 0] = static_cast<unsigned char>(clamp(cur->r + (next->r - cur->r)*u)*255.0);
	        table[i*4 + 1] = static_cast<unsigned char>(clamp(cur->g + (next->g - cur->g)*u)*255.0);
	        table[i*4 + 2] = static_cast<unsigned char>(clamp(cur->b + (next->b - cur->b)*u)*255.0);
	      }
          }

         bool opacityCubed = true; // default setting for test, this must be fixed

         for(opacity_nodes::const_iterator cur = local_cti.opacityNodes().begin();
   	    cur != boost::prior(local_cti.opacityNodes().end());
	    cur++)
          {
	    opacity_nodes::const_iterator next = boost::next(cur);
	    unsigned int cur_idx = static_cast<unsigned int>(((cur->position - MIN_RANGE)/(MAX_RANGE - MIN_RANGE))*255.0);
	    unsigned int next_idx = cur_idx + 
	      static_cast<unsigned int>(((next->position - cur->position)/(MAX_RANGE - MIN_RANGE))*255.0);
	    for(unsigned int i = cur_idx; i <= next_idx; i++)
	      {
	        double u = (double(i-cur_idx)/double(next_idx-cur_idx));

	        // arand: trick to make opacity work better, 4-11-2011
                if(opacityCubed)
                  {
                    double v = clamp(cur->value + (next->value - cur->value)*u);
                    v = v*v*v;
                    table[i*4 + 3] = static_cast<unsigned char>(v*255.0);
                  }
                else
                  table[i*4 + 3] = static_cast<unsigned char>(clamp(cur->value + (next->value - cur->value)*u)*255.0);
   	     }
          }

         // arand, 6-8-2011
         // fix up the boundary cases which were causing a long standing bug
         color_nodes::const_iterator cur = local_cti.colorNodes().begin();
         table[0] = static_cast<unsigned char>(clamp(cur->r)*255.0);
         table[1] = static_cast<unsigned char>(clamp(cur->g)*255.0);
         table[2] = static_cast<unsigned char>(clamp(cur->b)*255.0);

         cur = local_cti.colorNodes().end();
         cur--;
         table[4*(size-1)+0] = static_cast<unsigned char>(clamp(cur->r)*255.0);
         table[4*(size-1)+1] = static_cast<unsigned char>(clamp(cur->g)*255.0);
         table[4*(size-1)+2] = static_cast<unsigned char>(clamp(cur->b)*255.0);

         opacity_nodes::const_iterator cur1 = local_cti.opacityNodes().begin();
         double v = clamp(cur1->value);
         if(opacityCubed) v = v*v*v;	  
         table[3] = static_cast<unsigned char>(v*255.0);

         cur1 = local_cti.opacityNodes().end();
         cur1--;
         v = clamp(cur1->value);
         if(opacityCubed) v = v*v*v;	  
         table[4*(size-1)+3] = static_cast<unsigned char>(v*255.0);  
      }

      // write out into file
      fprintf( ofp, "%d\n", size );
      //outf << size << endl;
      for(int i = 0; i < size; i++ )
         fprintf( ofp, "%d %d %d %d\n", table[ i*4 ], table[ i*4 +1 ], table[ i*4 +2 ], table[ i*4 +3 ] );
         //outf << table[ i*4 ] << " " << table[ i*4 + 1 ] << " " << table[ i*4 + 2 ] << " " << table[ i*4 + 3 ] << endl;
      fclose( ofp );
  }

  void ColorTable::interactiveUpdates(bool b)
  {
    _xoomedIn->interactiveUpdates(b);
  }

  void ColorTable::opacityCubed(bool b)
  {
    _opacityCubed = b;
  }
}
