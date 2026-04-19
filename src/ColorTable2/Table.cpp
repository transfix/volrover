#include <iostream>
#include <map>
#include <algorithm>
#include <utility>
#include <boost/utility.hpp>
#include <boost/current_function.hpp>
#include <boost/algorithm/minmax_element.hpp>
#include <boost/foreach.hpp>
#include <cmath>
#include <cstdlib>
#include <ColorTable2/Table.h>
#include <QGLViewer/manipulatedCameraFrame.h>

#if QT_VERSION < 0x040000
#include <qpopupmenu.h>
#include <qfiledialog.h>
#include <qcolordialog.h>
#else
#include <QMouseEvent>
#include <QContextMenuEvent>
#include <QMenu>
#include <QFileDialog>
#include <QColorDialog>
#endif

#ifndef COLORTABLE2_DISABLE_CONTOUR_TREE
#include <contourtree/computeCT.h>
#endif

#ifndef COLORTABLE2_DISABLE_CONTOUR_SPECTRUM
#include <Contour/contour.h>
#include <Contour/datasetreg3.h>
#endif


#ifdef __GNUC__
#pragma message("TODO: implement vertical rendering, and make sure everything depends on MAX_RANGE/MIN_RANGE")
#endif

namespace CVCColorTable
{
  const float HANDLE_SIZE = 5.0; //size of draggable handles

  //in the screen coordinate system, 0.0 is the near plane, and 1.0 is far plane
  const float COLOR_BACKGROUND_LAYER = 0.9;
  const float HISTOGRAM_LAYER = 0.85;
  const float COLOR_BAR_LAYER = 0.8;
  const float RANGE_BAR_LAYER = 0.7;
  const float CONTOUR_TREE_LAYER = 0.6;
  const float CONTOUR_SPECTRUM_LAYER = 0.5;
  const float ISOCONTOUR_BAR_LAYER = 0.4;
  const float OPACITY_NODE_LAYER = 0.3;
  // arand, 8-24-2011: reordered layers above so that obacity nodes are 
  //                   always easy to move...


  Table::Table( ColorTable::color_table_info& cti, QWidget *parent, 
#if QT_VERSION < 0x040000 || defined QT3_SUPPORT
                      const char *name
#else
                      Qt::WindowFlags flags
#endif
                      )
    : QGLViewer(parent, flags), _cti(cti),
      _min(MIN_RANGE), _max(MAX_RANGE),
      _constraint(new qglviewer::WorldConstraint()),
      _selectedObj(-1), _interactiveUpdates(true),
      _visibleComponents(BACKGROUND |
                         COLOR_BARS |
                         ISOCONTOUR_BARS |
                         OPACITY_NODES),
      _rangeMin(MIN_RANGE),
      _rangeMax(MAX_RANGE),
      _dirtyContourTree(true),
      _dirtyContourSpectrum(true),
      _dirtyHistogram(true)
  {
      allocateInformDialg();
  }

  Table::Table( boost::uint64_t components,
                      ColorTable::color_table_info& cti, QWidget *parent, 
#if QT_VERSION < 0x040000 || defined QT3_SUPPORT
                      const char *name
#else
                      Qt::WindowFlags flags
#endif
                      )
    : QGLViewer(parent, flags), _cti(cti),
      _min(MIN_RANGE), _max(MAX_RANGE),
      _constraint(new qglviewer::WorldConstraint()),
      _selectedObj(-1), _interactiveUpdates(true),
      _visibleComponents(components),
      _rangeMin(MIN_RANGE),
      _rangeMax(MAX_RANGE),
      _dirtyContourTree(true),
      _dirtyContourSpectrum(true),
      _dirtyHistogram(true)
  {
      allocateInformDialg();
  }

  Table::~Table()
  {
     if( m_DrawInformDialog ) delete m_DrawInformDialog;
  }

  void Table::rangeMin(double val)
  {
    _rangeMin = val;
    emit rangeMinChanged(val);
    update();
  }

  void Table::rangeMax(double val)
  {
    _rangeMax = val;
    emit rangeMaxChanged(val);
    update();
  }

  void Table::interactiveUpdates(bool b)
  { 
    _interactiveUpdates = b; 
  }

  void Table::visibleComponents(boost::uint64_t components)
  {
    _visibleComponents = components;
    update();
  }

  void Table::setContourVolume(const VolMagick::Volume& vol)
  {
    _contourVolume = vol;
    _dirtyContourTree = true;
    _dirtyContourSpectrum = true;
    _dirtyHistogram = true;
    update();
  }

  void Table::setMin(double min)
  {
    _min = std::max(MIN_RANGE,std::min(MAX_RANGE,min)); //silently clamp
    //ensure min <= max
    double tmp_min = std::min(_min,_max);
    double tmp_max = std::max(_min,_max);
    _min = tmp_min;
    _max = tmp_max;
    update();
  }
 
  void Table::setMax(double max)
  {
    _max = std::max(MIN_RANGE,std::min(MAX_RANGE,max)); //silently clamp
    //ensure min <= max
    double tmp_min = std::min(_min,_max);
    double tmp_max = std::max(_min,_max);
    _min = tmp_min;
    _max = tmp_max;
    update();
  }

  void Table::showOpacityFunction(bool b)
  {
    visibleComponents(b ?
                      visibleComponents() | OPACITY_NODES :
                      visibleComponents() & ~OPACITY_NODES);
    update();
  }

  void Table::showTransferFunction(bool b)
  {
    visibleComponents(b ?
                      visibleComponents() | BACKGROUND :
                      visibleComponents() & ~BACKGROUND);

    visibleComponents(b ?
                      visibleComponents() | COLOR_BARS :
                      visibleComponents() & ~COLOR_BARS);
    update();
  }

  void Table::showContourTree(bool b)
  {
    visibleComponents(b ?
                      visibleComponents() | CONTOUR_TREE :
                      visibleComponents() & ~CONTOUR_TREE);
    update();
  }

  void Table::allocateInformDialg()
  {
     m_DrawInformDialog = new InfoDialog();
  }

  void Table::showInformDialog(bool b)
  {
     updateInformDialog( -1, 0.0, SHOW );
     b?m_DrawInformDialog->show():m_DrawInformDialog->hide();
  }

  void Table::showContourSpectrum(bool b)
  {
    visibleComponents(b ?
                      visibleComponents() | CONTOUR_SPECTRUM :
                      visibleComponents() & ~CONTOUR_SPECTRUM);
    update();
  }

  void Table::showHistogram(bool b)
  {
    visibleComponents(b ?
                      visibleComponents() | HISTOGRAM :
                      visibleComponents() & ~HISTOGRAM);
    update();
  }

  void Table::init()
  {
    //this is probably not necessary because of the start/stop screen coordinates calls
    _constraint->setRotationConstraintType(qglviewer::AxisPlaneConstraint::FORBIDDEN);
    _constraint->setTranslationConstraintType(qglviewer::AxisPlaneConstraint::FORBIDDEN);
    // Qt6: camera() returns Camera*, call its manipulatedFrame() method
    if (camera() && camera()->frame())
      camera()->frame()->setConstraint(_constraint.get());
  }

  void Table::draw()
  {
    startScreenCoordinatesSystem();
    drawTable();
    stopScreenCoordinatesSystem();
  }

  void Table::drawWithNames()
  {
    //begin selection already set up the screen coord matrix
    drawTable(true);
    //qDebug("%s: glGetError: %d",BOOST_CURRENT_FUNCTION,glGetError());
  }

  void Table::drawTable(bool withNames)
  {
    using namespace boost; //for std::next() and prior()
    const GLfloat color_bar_node_color[3] = { 1.0, 0.0, 0.0 };
    const GLfloat isocontour_bar_node_color[3] = { 0.0, 1.0, 0.0 };
    const GLfloat opacity_node_color[3] = { 0.0, 0.0, 1.0 };
    const GLfloat range_bar_node_color[3] = { 0.5, 0.5, 0.5 };
    const GLfloat contour_tree_color[3] = { 0.0, 0.0, 0.0 };
    int name = 1;

    if(withNames)
      {
	_nameMap.clear();
	_nameMap.resize(_cti.colorNodes().size() +
			_cti.isocontourNodes().size() +
			_cti.opacityNodes().size() +
                        2 /* range bars */);
      }

    glDisable(GL_LIGHTING);

    //qDebug("glGetError: %d",glGetError());
    
    /*********** Color Background ************/
    if(visibleComponents() & BACKGROUND)
      {
        //if we don't have enough  color nodes, just draw a grayscale thing
        if(_cti.colorNodes().size() < 2)
          {
            glBegin(GL_QUADS);
            glColor3f(_min,_min,_min);
            glVertex3f(0,0, COLOR_BACKGROUND_LAYER);
            glColor3f(_max,_max,_max);
            glVertex3f(width()-1,0, COLOR_BACKGROUND_LAYER);
            glColor3f(_max,_max,_max);
            glVertex3f(width()-1,height()-1, COLOR_BACKGROUND_LAYER);
            glColor3f(_min,_min,_min);
            glVertex3f(0,height()-1, COLOR_BACKGROUND_LAYER);
            glEnd();
          }
        else
          {
            glBegin(GL_QUADS);

            for(ColorTable::color_nodes::const_iterator i = _cti.colorNodes().begin();
                i != prior(_cti.colorNodes().end());
                i++)
              {
                double x0 = ((i->position - _min)/(_max - _min))*(width()-1);
                double x1 = ((std::next(i)->position - _min)/(_max - _min))*(width()-1);

		//std::cout << "x0 == " << x0 << std::endl;
		//std::cout << "x1 == " << x1 << std::endl;

                glColor3f(i->r,i->g,i->b);
                glVertex3f(x0,0.0, COLOR_BACKGROUND_LAYER);
                glColor3f(std::next(i)->r,std::next(i)->g,std::next(i)->b);
                glVertex3f(x1,0.0, COLOR_BACKGROUND_LAYER);
                glColor3f(std::next(i)->r,std::next(i)->g,std::next(i)->b);
                glVertex3f(x1,height()-1.0, COLOR_BACKGROUND_LAYER);
                glColor3f(i->r,i->g,i->b);
                glVertex3f(x0,height()-1.0, COLOR_BACKGROUND_LAYER);
              }

            glEnd();
          }
      } else if (visibleComponents()) 
      {
	// arand: draw a solid background color
	//        for now, a light grey is ok...
	// FUTURE: allow the user to select this color	
	glBegin(GL_QUADS);
	glColor3f(.5, .5, .5);
	glVertex3f(0,0, COLOR_BACKGROUND_LAYER);
	glColor3f(.5, .5, .5);
	glVertex3f(width()-1,0, COLOR_BACKGROUND_LAYER);
	glColor3f(.5, .5, .5);
	glVertex3f(width()-1,height()-1, COLOR_BACKGROUND_LAYER);
	glColor3f(.5, .5, .5);
	glVertex3f(0,height()-1, COLOR_BACKGROUND_LAYER);
	glEnd();	  
      }

    //qDebug("glGetError: %d",glGetError());
#if 1
    /*********** Color Bars ************/
    if(visibleComponents() & COLOR_BARS)
      {
        for(ColorTable::color_nodes::const_iterator i = _cti.colorNodes().begin();
            i != _cti.colorNodes().end();
            i++)
          {
            drawBar(((i->position - _min)/(_max - _min))*(width()-1),
                    COLOR_BAR_LAYER,
                    color_bar_node_color,
                    !withNames ? -1 : name);
            if(withNames) _nameMap[name-1] = i;
            name++;
          }
      }
#endif
    //qDebug("glGetError: %d",glGetError());
    
#if 1
    /*********** Isocontour Bars ************/
    if(visibleComponents() & ISOCONTOUR_BARS)
      {
        for(ColorTable::isocontour_nodes::const_iterator i = _cti.isocontourNodes().begin();
            i != _cti.isocontourNodes().end();
            i++)
          {
            drawBar(((i->position - _min)/(_max - _min))*(width()-1),
                    ISOCONTOUR_BAR_LAYER,
                    isocontour_bar_node_color,
                    !withNames ? -1 : name);
            if(withNames) _nameMap[name-1] = i;
            name++;
          }
      }
#endif
    //qDebug("glGetError: %d",glGetError());

#if 1
    /*********** Range Bars ************/
    if(visibleComponents() & RANGE_BARS)
      {
        drawBar(((_rangeMin - _min)/(_max - _min))*(width()-1),
                RANGE_BAR_LAYER,
                range_bar_node_color,
                !withNames ? -1 : name);
        if(withNames) _nameMap[name-1] = &_rangeMin;
        name++;

        drawBar(((_rangeMax - _min)/(_max - _min))*(width()-1),
                RANGE_BAR_LAYER,
                range_bar_node_color,
                !withNames ? -1 : name);
        if(withNames) _nameMap[name-1] = &_rangeMax;
        name++;
      }
#endif
    //qDebug("glGetError: %d",glGetError());

#if 1
    /*********** Opacity Nodes ************/
    if(visibleComponents() & OPACITY_NODES)
      {
        //if we don't have enough, lets pretend we have 2 nodes, making a linear ramp from MIN_RANGE to MAX_RANGE
        if(_cti.opacityNodes().size() < 2)
          {
            //draw ramp
            glBegin(GL_LINES);
            /*glColor3f(1.0-opacity_node_color[0],
              1.0-opacity_node_color[1],
              1.0-opacity_node_color[2]);*/
            glColor3f(0.0,0.0,0.0); //black line
            glVertex3f(0.0, height()-1.0, OPACITY_NODE_LAYER);
            glVertex3f(width()-1.0, 0.0, OPACITY_NODE_LAYER);
            glEnd();

            glBegin(GL_QUADS);
            glColor3fv(opacity_node_color);
	
            //opacity is 0.0 at 0.0
            glVertex3f(0.0 - HANDLE_SIZE,height()-1.0 + HANDLE_SIZE, OPACITY_NODE_LAYER-0.01);
            glVertex3f(0.0 + HANDLE_SIZE,height()-1.0 + HANDLE_SIZE, OPACITY_NODE_LAYER-0.01);
            glVertex3f(0.0 + HANDLE_SIZE,height()-1.0 - HANDLE_SIZE, OPACITY_NODE_LAYER-0.01);
            glVertex3f(0.0 - HANDLE_SIZE,height()-1.0 - HANDLE_SIZE, OPACITY_NODE_LAYER-0.01);

            //opacity is 1.0 at 1.0
            glVertex3f(width()-1.0 - HANDLE_SIZE,0.0 + HANDLE_SIZE, OPACITY_NODE_LAYER-0.01);
            glVertex3f(width()-1.0 + HANDLE_SIZE,0.0 + HANDLE_SIZE, OPACITY_NODE_LAYER-0.01);
            glVertex3f(width()-1.0 + HANDLE_SIZE,0.0 - HANDLE_SIZE, OPACITY_NODE_LAYER-0.01);
            glVertex3f(width()-1.0 - HANDLE_SIZE,0.0 - HANDLE_SIZE, OPACITY_NODE_LAYER-0.01);
            glEnd();
          }
        else
          {
            for(ColorTable::opacity_nodes::const_iterator i = _cti.opacityNodes().begin();
                i != prior(_cti.opacityNodes().end());
                i++)
              {
                double x0 = ((i->position - _min)/(_max - _min))*(width()-1);
                double x1 = ((std::next(i)->position - _min)/(_max - _min))*(width()-1);
                double y0 = (1.0 - i->value)*(height()-1);
                double y1 = (1.0 - std::next(i)->value)*(height()-1);

                //draw line between these nodes
                glBegin(GL_LINES);
                /*glColor3f(1.0-opacity_node_color[0],
                  1.0-opacity_node_color[1],
                  1.0-opacity_node_color[2]);*/
                glColor3f(0.0,0.0,0.0); //black line
                glVertex3f(x0,y0, OPACITY_NODE_LAYER);
                glVertex3f(x1,y1, OPACITY_NODE_LAYER);
                glEnd();

                if(withNames)
                  {
                    //glPushName(name++);
                    //glPushName(8);
                    glLoadName(name);
                    _nameMap[name-1] = i;
                    name++;
                  }
                glBegin(GL_QUADS);
                glColor3fv(opacity_node_color);
                glVertex3f(x0 - HANDLE_SIZE,y0 + HANDLE_SIZE, OPACITY_NODE_LAYER-0.01);
                glVertex3f(x0 + HANDLE_SIZE,y0 + HANDLE_SIZE, OPACITY_NODE_LAYER-0.01);
                glVertex3f(x0 + HANDLE_SIZE,y0 - HANDLE_SIZE, OPACITY_NODE_LAYER-0.01);
                glVertex3f(x0 - HANDLE_SIZE,y0 - HANDLE_SIZE, OPACITY_NODE_LAYER-0.01);
                glEnd();
                if(withNames)
                  {
                    //glPopName();
                    glLoadName(0);
                  }
              }

            //draw the last node
            {
              ColorTable::opacity_nodes::const_iterator i = prior(_cti.opacityNodes().end());
              double x0 = ((i->position - _min)/(_max - _min))*(width()-1);
              double y0 = (1.0 - i->value)*(height()-1);
              if(withNames)
                {
                  //glPushName(name++);
                  glLoadName(name);
                  _nameMap[name-1] = i;
                  name++;
                }
              glBegin(GL_QUADS);
              glColor3fv(opacity_node_color);
              glVertex3f(x0 - HANDLE_SIZE,y0 + HANDLE_SIZE, OPACITY_NODE_LAYER-0.01);
              glVertex3f(x0 + HANDLE_SIZE,y0 + HANDLE_SIZE, OPACITY_NODE_LAYER-0.01);
              glVertex3f(x0 + HANDLE_SIZE,y0 - HANDLE_SIZE, OPACITY_NODE_LAYER-0.01);
              glVertex3f(x0 - HANDLE_SIZE,y0 - HANDLE_SIZE, OPACITY_NODE_LAYER-0.01);
              glEnd();
              if(withNames)
                {
                  //glPopName();
                  glLoadName(0);
                }
            }
          }
      }

    if(visibleComponents() & CONTOUR_TREE)
      {
        computeContourTree();

        for(std::vector<int>::iterator i = _contourTreeEdges.begin();
            i != _contourTreeEdges.end();
            i++)
          {
            if(boost::next(i) == _contourTreeEdges.end()) break;

            int v1 = *i;
            int v2 = *boost::next(i);

            double p1x,p1y, p2x,p2y;
            double clampedP1x,clampedP1y,clampedP2x,clampedP2y,mag;
            
            // get the vertex coords
            p1x = _contourTreeVertices[v1 * 2 + 0]; //func_val
            p1y = _contourTreeVertices[v1 * 2 + 1]; //norm_x
            p2x = _contourTreeVertices[v2 * 2 + 0]; //func_val
            p2y = _contourTreeVertices[v2 * 2 + 1]; //norm_x;

            //qDebug("%s :: p1x = %f, p1y = %f, p2x = %f, p2y = %f", BOOST_CURRENT_FUNCTION, p1x, p1y, p2x, p2y);            

            // don't add edges that are outside of the zoomed in region
            if (p1x < _min && p2x < _min) continue;
            if (p1x > _max && p2x > _max) continue;

            // find the cases where edges span the edge of our bounding box
            if (p1x < p2x) {
              if (p1x < _min) {
                mag = (_min-p1x)/(p2x-p1x);
                clampedP1y = p1y * (1.0-mag) + p2y * mag;
                clampedP1x = _min;
              }
              else {
                clampedP1x = p1x;
                clampedP1y = p1y;
              }
              if (p2x > _max) {
                mag = (_max-p1x)/(p2x-p1x);
                clampedP2y = p1y * (1.0-mag) + p2y * mag;
                clampedP2x = _max;
              }
              else {
                clampedP2x = p2x;
                clampedP2y = p2y;
              }
            }
            else { // p2x <= p1x (but probably p2x < p1x)
              if (p2x < _min) {
                mag = (_min-p2x)/(p1x-p2x);
                clampedP1y = p2y * (1.0-mag) + p1y * mag;
                clampedP1x = _min;
              }
              else {
                clampedP1x = p2x;
                clampedP1y = p2y;
              }
              if (p1x > _max) {
                mag = (_max-p2x)/(p1x-p2x);
                clampedP2y = p2y * (1.0-mag) + p1y * mag;
                clampedP2x = _max;
              }
              else {
                clampedP2x = p1x;
                clampedP2y = p1y;
              }
            }

            //qDebug("%s :: x0 = %f, y0 = %f, x1 = %f, y1 = %f", BOOST_CURRENT_FUNCTION, 
            //       clampedP1x, clampedP1y, clampedP2x, clampedP2y);

            double x0 = ((clampedP1x - _min)/(_max - _min))*(width()-1);
            double x1 = ((clampedP2x - _min)/(_max - _min))*(width()-1);
            double y0 = (1.0 - clampedP1y)*(height()-1)*.8+height()*.1;
            double y1 = (1.0 - clampedP2y)*(height()-1)*.8+height()*.1;
	    // arand, 8-25-2011: added .8 above so that contour tree lines don't get lost
	    //                   in the window boundary

            glBegin(GL_LINES);
            glColor3fv(contour_tree_color);
            glVertex3f(x0,y0, CONTOUR_TREE_LAYER);
            glVertex3f(x1,y1, CONTOUR_TREE_LAYER);
            glEnd();
          }
      }



    if(visibleComponents() & HISTOGRAM) {
      // arand, 9-2-2011: initial histogram implementation
      // FIXME: the histogram does not scale as with the "zoom" feature
      //        of the transfer function...

      if (_dirtyHistogram) {
	// grab the histogram
	_histogram = _contourVolume.histogram();
	_dirtyHistogram = false;
      }
      // draw the histogram
     
      // loop over datapoints...

      const boost::uint64_t* hist = get<0>(_histogram);
      VolMagick::uint64 len = _histogram.get<1>();
      
      long histmax = 0;
      long histmin = 100000000;
      long histsecond = 0;
      for (int i=0; i<len; i++) {
	if (hist[i] >= histmax) {
	  histsecond = histmax;
	  histmax = hist[i];
	} else if (hist[i] > histsecond) {
	  histsecond = hist[i];
	}
	if (hist[i] < histmin) histmin = hist[i];       
      }

      if (histmax > 1.5*histsecond) {
	histmax = 1.5*histsecond;
      }
      
      if (histmax != histmin) { 
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
      
	glBegin(GL_QUADS);
	if (visibleComponents() & BACKGROUND)
	  glColor4d(0.25,0.25,0.25,.4);
	else
	  glColor4d(0,0,0,1);

	for (int i=0; i<len; i++) {
	  
	  double x0 = ((double)i)/((double)len)*(width()-1);
	  double x1 = ((double)i+1.0)/((double)len)*(width()-1);
	  double y = (1.0 - (hist[i]-histmin)/((double)(histmax-histmin)))*(height()-1);
	  
	  glVertex3f(x0,height(), HISTOGRAM_LAYER);
	  glVertex3f(x1,height(), HISTOGRAM_LAYER);
	  glVertex3f(x1,y, HISTOGRAM_LAYER);
	  glVertex3f(x0,y, HISTOGRAM_LAYER);
	  
	}	
	glDisable(GL_BLEND);

	glEnd();
      }

    }

    if(visibleComponents() & CONTOUR_SPECTRUM)
      {
        computeContourSpectrum();

        for(function_map::iterator curFunc = _contourSpectrum.begin();
            curFunc != _contourSpectrum.end();
            curFunc++)
          {
            color_t &specFuncColor = curFunc->second.first;
            function_t &specFunc = curFunc->second.second;

            double curPos,nextPos, clampedCurPos,clampedNextPos,mag;
            double func1, func2;

            glBegin(GL_LINES);
            glColor3dv(&(specFuncColor[0]));
	
            // loop through the function table
            for (size_t i=0; i < specFunc.size(); i++) 
              {
                // current and next index normalized X coords
                curPos = i / double(specFunc.size());
                nextPos = (i+1) / double(specFunc.size());
	
                if(i+1 >= specFunc.size()) break;

                // make sure the points are within range
                if (nextPos < _min) continue;
                if (curPos > _max) break;
                
                // if needed find the intersection of the line with the min and max
#if 0
                if( curPos < _min ) 
                  {
                    mag = (_min-curPos)/(nextPos-curPos);
                    func1 = specFunc[i] * (1.0-mag) + specFunc[i+1] * (mag);
                    clampedCurPos = _min;
                  }
                else 
                  {
                    func1 = specFunc[i];
                    clampedCurPos = curPos;
                  }
                if( nextPos > _max ) 
                  {
                    mag = (_max-curPos)/(nextPos-curPos);
                    func2 = specFunc[i] * (1.0-mag) + specFunc[i+1] * (mag);
                    clampedNextPos = _max;
                  }
                else 
                  {
                    func2 = specFunc[i+1];
                    clampedNextPos = nextPos;
                  }
#endif

                clampedCurPos = curPos;
                clampedNextPos = nextPos;
                func1 = specFunc[i];
                func2 = specFunc[i+1];

                double x0 = ((clampedCurPos - _min)/(_max - _min))*(width()-1);
                double x1 = ((clampedNextPos - _min)/(_max - _min))*(width()-1);
                double y0 = (1.0 - func1)*(height()-1);
                double y1 = (1.0 - func2)*(height()-1);
                
                glVertex3f(x0,y0, CONTOUR_SPECTRUM_LAYER);
                glVertex3f(x1,y1, CONTOUR_SPECTRUM_LAYER);
              }
            
            glEnd();
          }
      }
#endif
    //qDebug("glGetError: %d",glGetError());
  }

  void Table::computeContourTree()
  {
     if( !_dirtyContourTree ) return;
  
     int dim[3] = { static_cast<int>(_contourVolume.XDim()),
                    static_cast<int>(_contourVolume.YDim()),
                    static_cast<int>(_contourVolume.ZDim()) };

     _contourTreeVertices.clear();
     _contourTreeEdges.clear();

#ifndef COLORTABLE2_DISABLE_CONTOUR_TREE
     if(_contourVolume.voxelType() != cvc::UChar)
     {
          _contourVolume.map(0.0,255.0);
          _contourVolume.voxelType(cvc::UChar);
     }

     CTVTX* verts = NULL;
     CTEDGE* edges = NULL;
     int no_vtx = 0, no_edge = 0;
            
     computeCT(*_contourVolume,dim,no_vtx,no_edge,&verts,&edges);
 
     if(verts)
       {
        //normalize the vertices
        float xmin = verts[0].norm_x, xmax=xmin,
        ymin = verts[0].func_val, ymax=ymin;
	int i;
		
	for (i=0; i<no_vtx; i++)
          {
            if (verts[i].func_val < ymin) ymin = verts[i].func_val;
            else if (verts[i].func_val > ymax) ymax = verts[i].func_val;
                    
            if (verts[i].norm_x < xmin) xmin = verts[i].norm_x;
            else if (verts[i].norm_x > xmax) xmax = verts[i].norm_x;
          }
	for (i=0; i<no_vtx; i++) 
          {
             verts[i].func_val = (verts[i].func_val-ymin) / (ymax-ymin);
             verts[i].norm_x = (verts[i].norm_x-xmin) / (xmax-xmin);
          }
        for(int i = 0; i < no_vtx; i++)
          {
             _contourTreeVertices.push_back(verts[i].func_val);
             _contourTreeVertices.push_back(verts[i].norm_x);
          }
       }
     if(edges)
       {
         for(int i = 0; i < no_edge; i++)
          {
             _contourTreeEdges.push_back(edges[i].v1);
             _contourTreeEdges.push_back(edges[i].v2);
          }
       }

     std::free(verts);
     std::free(edges);
#endif
     _dirtyContourTree = false;
  }

  void Table::computeContourSpectrum()
  {
     if( !_dirtyContourSpectrum ) return;

#ifndef COLORTABLE2_DISABLE_CONTOUR_SPECTRUM
     float span[3], orig[3];
     ConDataset* the_data;
     Signature	*sig;
     int dim[3] = 
     {
        _contourVolume.XDim(),
        _contourVolume.YDim(),
        _contourVolume.ZDim()
     };

     orig[0] = _contourVolume.XMin();
     orig[1] = _contourVolume.YMin();
     orig[2] = _contourVolume.ZMin();
            
     span[0] = _contourVolume.XSpan();
     span[1] = _contourVolume.YSpan();
     span[2] = _contourVolume.ZSpan();

     if(_contourVolume.voxelType() != cvc::UChar)
     {
        _contourVolume.map(0.0,255.0);
        _contourVolume.voxelType(cvc::UChar);
     }

     // make a libcontour variable out of dataBuffer
     the_data = newDatasetReg(CONTOUR_UCHAR, CONTOUR_REG_3D, 1, 1, dim, *_contourVolume);
     ((Datareg3 *)the_data->data->getData(0))->setOrig(orig);
     ((Datareg3 *)the_data->data->getData(0))->setSpan(span);
     // compute the contour spectrum
     sig = getSignatureFunctions(the_data, 0,0);

     color_t isovalue_color = { 0.0, 0.0, 0.0 };
     color_t area_color     = { 1.0, 0.0, 0.0 };
     color_t min_vol_color  = { 0.0, 1.0, 0.0 };
     color_t max_vol_color  = { 0.0, 0.0, 1.0 };
     color_t gradient_color = { 1.0, 1.0, 0.0 };
     _contourSpectrum["isovalue"].first = isovalue_color;
     _contourSpectrum["area"].first = area_color;
     _contourSpectrum["min_vol"].first = min_vol_color;
     _contourSpectrum["max_vol"].first = max_vol_color;
     _contourSpectrum["gradient"].first = gradient_color;

     _contourSpectrum["isovalue"].second.clear();
     _contourSpectrum["area"].second.clear();
     _contourSpectrum["min_vol"].second.clear();
     _contourSpectrum["max_vol"].second.clear();
     _contourSpectrum["gradient"].second.clear();

     for(int i = 0; i < 256; i++)
     {
         _contourSpectrum["isovalue"].second.push_back(sig[0].fx[i]);
         _contourSpectrum["area"].second.push_back(sig[0].fy[i]);
         _contourSpectrum["min_vol"].second.push_back(sig[1].fy[i]);
         _contourSpectrum["max_vol"].second.push_back(sig[2].fy[i]);
         _contourSpectrum["gradient"].second.push_back(sig[3].fy[i]);
     }

     bool tmp = true;

     //normalize the functions
     for(function_map::iterator i = _contourSpectrum.begin();
        i != _contourSpectrum.end();
        i++)
     {
        function_t &func = i->second.second;
        if(func.empty()) continue;
        typedef function_t::iterator iterator;
        std::pair<iterator,iterator> result =
        boost::minmax_element(func.begin(),
                                func.end());

	float mymin = *(result.first);
	float mymax = *(result.second);

        for (auto& val : func) {
          val = (val-mymin)/(mymax-mymin);
        }
     }

     delete the_data;
     delete sig;
#endif
     _dirtyContourSpectrum = false;
  }

  void Table::drawBar(double x_pos, double depth, const GLfloat *color_3f, GLint name)
  {
    if(name != -1)
      //glPushName(name);
      glLoadName(name);

    GLint params[2];
    //back up current setting
    //glGetIntegerv(GL_POLYGON_MODE,params);

    //draw bar extending from node with negative color
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    glBegin(GL_QUADS);
    glColor3f(1.0-color_3f[0],
	      1.0-color_3f[1],
	      1.0-color_3f[2]);
    glVertex3f(x_pos - 1.0, height()-1, depth);
    glVertex3f(x_pos + 1.0, height()-1, depth);
    glVertex3f(x_pos + 1.0, 0.0, depth);
    glVertex3f(x_pos - 1.0, 0.0, depth);
    glEnd();

    //draw filled in node
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    glBegin(GL_QUADS);
    glColor3f(color_3f[0],
	      color_3f[1],
	      color_3f[2]);
    glVertex3f(x_pos - HANDLE_SIZE, height()-1, depth-0.01);
    glVertex3f(x_pos + HANDLE_SIZE, height()-1, depth-0.01);
    glVertex3f(x_pos + HANDLE_SIZE, height()-1 - HANDLE_SIZE*2, depth-0.01);
    glVertex3f(x_pos - HANDLE_SIZE, height()-1 - HANDLE_SIZE*2, depth-0.01);
    glEnd();

    //draw outline of node with negative color
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    glBegin(GL_QUADS);
    glColor3f(1.0-color_3f[0],
	      1.0-color_3f[1],
	      1.0-color_3f[2]);
    glVertex3f(x_pos - HANDLE_SIZE, height()-1, depth-0.02);
    glVertex3f(x_pos + HANDLE_SIZE, height()-1, depth-0.02);
    glVertex3f(x_pos + HANDLE_SIZE, height()-1 - HANDLE_SIZE*2, depth-0.02);
    glVertex3f(x_pos - HANDLE_SIZE, height()-1 - HANDLE_SIZE*2, depth-0.02);
    glEnd();

    //restore previous setting for polygon mode
    //glPolygonMode(GL_FRONT,params[0]);
    //glPolygonMode(GL_BACK,params[1]);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

    if(name != -1)
    //glPopName();
      glLoadName(0);
  }

  void Table::beginSelection(const QPoint& point)
  {
#if 1
    // Make OpenGL context current (may be needed with several viewers ?)
    makeCurrent();

    // Prepare the selection mode
    glSelectBuffer(selectBufferSize(), selectBuffer());
    glRenderMode(GL_SELECT);
    glInitNames();
    glPushName(0);

    // Loads the matrices
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    static GLint viewport[4];
    camera()->getViewport(viewport);
    gluPickMatrix(point.x(), point.y(), selectRegionWidth(), selectRegionHeight(), viewport);

    glOrtho(0, width(), height(), 0, 0.0, -1.0);
    glMatrixMode(GL_MODELVIEW);
    //glPushMatrix();
    glLoadIdentity();
#endif
  }

  void Table::postSelection(const QPoint& point)
  {
#if 1
    //only apply selection if there isn't something already selected
    if(_selectedObj != -1) return;

    qDebug("Selected: %d", selectedName());
    _selectedPoint = point;
    _selectedObj = selectedName();
    _selectedObj = _selectedObj == 0 ? -1 : _selectedObj; //0 also means selected nothing
#endif
  }
  
  void Table::mousePressEvent(QMouseEvent *e)
  {
#if QT_VERSION < 0x040000
    if(e->button() == QMouseEvent::LeftButton)
#else
    if(e->button() == Qt::LeftButton)
#endif
      {
	select(e->pos());
	if(_selectedObj != -1)
	  _selectedPoint = e->pos();
      }

    QGLViewer::mousePressEvent(e);
  }

  void Table::mouseMoveEvent(QMouseEvent *e)
  {
    using namespace boost; //for any_cast<>() and prior()

    if(_selectedObj != -1)
      {
	double newpos = _min + e->pos().x()*((_max-_min)/(width()-1));
	double newopac = 1.0 - e->pos().y()*(1.0/(height()-1));

	//clamp
	newpos = std::max(MIN_RANGE,std::min(MAX_RANGE,newpos));
	newopac = std::max(MIN_RANGE,std::min(MAX_RANGE,newopac));

	//qDebug("newpos = %f, newopac = %f",newpos,newopac);

	if(_nameMap[_selectedObj-1].type() == 
	   typeid(ColorTable::isocontour_nodes::const_iterator))
	  {
	    ColorTable::isocontour_nodes::const_iterator node_itr = 
	      any_cast<ColorTable::isocontour_nodes::const_iterator>(_nameMap[_selectedObj-1]);

	    ColorTable::isocontour_node node = *node_itr;
	    node.position = newpos;

            //Update contour spectrum/tree info. dialog
            updateInformDialog( node.id(), newpos, MOVE );

	    //Only change this node if it wont consume a nearby node via unsuccessful re-insertion.
	    if(_cti.isocontourNodes().find(node) != _cti.isocontourNodes().end())
	      return;
	    _cti.isocontourNodes().erase(node_itr);
	    _cti.isocontourNodes().insert(node);
	  }
	else if(_nameMap[_selectedObj-1].type() ==
		typeid(ColorTable::color_nodes::const_iterator))
	  {
	    ColorTable::color_nodes::const_iterator node_itr =
	      any_cast<ColorTable::color_nodes::const_iterator>(_nameMap[_selectedObj-1]);

	    //We cannot allow the first and last color nodes to change position!
	    if(node_itr == _cti.colorNodes().begin())
	      newpos = MIN_RANGE;
	    else if(node_itr == prior(_cti.colorNodes().end()))
	      newpos = MAX_RANGE;

	    ColorTable::color_node node = *node_itr;
	    node.position = newpos;
	    //Only change this node if it wont consume a nearby node via unsuccessful re-insertion.
	    //This doesn't happen for the first and last nodes, so check for them.
	    if(node_itr != _cti.colorNodes().begin() &&
	       node_itr != prior(_cti.colorNodes().end()) &&
	       _cti.colorNodes().find(node) != _cti.colorNodes().end())
	      return;
	    _cti.colorNodes().erase(node_itr);
	    _cti.colorNodes().insert(node);
	  }
	else if(_nameMap[_selectedObj-1].type() ==
		typeid(ColorTable::opacity_nodes::const_iterator))
	  {
	    ColorTable::opacity_nodes::const_iterator node_itr =
	      any_cast<ColorTable::opacity_nodes::const_iterator>(_nameMap[_selectedObj-1]);

	    //We cannot allow the first and last opacity nodes to change position!
	    //However, opacity change for those nodes is ok...
	    if(node_itr == _cti.opacityNodes().begin())
	      newpos = MIN_RANGE;
	    else if(node_itr == prior(_cti.opacityNodes().end()))
	      newpos = MAX_RANGE;
	    
	    ColorTable::opacity_node node = *node_itr;
	    node.position = newpos;
	    node.value = newopac;
	    //Only change this node if it wont consume a nearby node via unsuccessful re-insertion.
	    //This doesn't happen for the first and last nodes, so check for them.
	    if(node_itr != _cti.opacityNodes().begin() &&
	       node_itr != prior(_cti.opacityNodes().end()) &&
	       _cti.opacityNodes().find(node) != _cti.opacityNodes().end())
	      return;
	    _cti.opacityNodes().erase(node_itr);
	    _cti.opacityNodes().insert(node);
	  }
        else if(_nameMap[_selectedObj-1].type() == typeid(double*))
          {
            double *rangeVar = any_cast<double*>(_nameMap[_selectedObj-1]);
            
            *rangeVar = newpos;

            if(rangeVar == &_rangeMin && newpos > _rangeMax) //make sure min is less than max
              _rangeMax = newpos;
            if(rangeVar == &_rangeMax && newpos < _rangeMin) //make sure max is greater than min
              _rangeMin = newpos;

            if(rangeVar == &_rangeMin)
              emit rangeMinChanged(newpos);
            if(rangeVar == &_rangeMax)
              emit rangeMaxChanged(newpos);
          }

	//recalculate _nameMap because select id's probably changed if node order is different
	select(e->pos());

	update();
	if(_interactiveUpdates) emit changed();
      }

    QGLViewer::mouseMoveEvent(e);
  }

  void Table::mouseReleaseEvent(QMouseEvent *e)
  {
    std::cout << BOOST_CURRENT_FUNCTION << ": called!" << std::endl;

    if(_selectedObj != -1)
      {
	if(!_interactiveUpdates) emit changed();
	_selectedObj = -1; //de-select
      }

    QGLViewer::mouseReleaseEvent(e);
  }

  void Table::contextMenuEvent(QContextMenuEvent *e)
  {
    using namespace boost; //for std::next(), prior() & any_cast<>()

    std::cout << BOOST_CURRENT_FUNCTION << ": called!" << std::endl;

    bool modified = false;
    POPUPSELECTION selection = showPopup(e->globalPos());
    double newpos = _min + e->pos().x()*((_max-_min)/(width()-1));
    double newopac = 1.0 - e->pos().y()*(1.0/(height()-1));

    //do a selection to see if we're over an object
    select(e->pos());

    switch(selection)
      {
      case LOAD_MAP:
	{
	  QString filename;
#if QT_VERSION < 0x040000
	  filename = QFileDialog::getOpenFileName("", "Vinay files (*.vinay)", 
						  this, "open file dialog", "Choose a file");
#else
          filename = QFileDialog::getOpenFileName(this, tr("open file dialog"), 
                                                  "", tr("Vinay files (*.vinay)"));
#endif
	  if(filename.isNull()) return;
          const char *c_filename =
#if QT_VERSION < 0x040000
            filename.ascii()
#else
            filename.toUtf8().constData()
#endif
            ;
	  info() = ColorTable::read_transfer_function(c_filename);
	  modified = true;
	}
	break;
      case SAVE_MAP:
	{
	  QString filename;
          //SimpleFileSaveDialog dialog;
          QFileDialog dialog( this, tr("Save Color Map Dialog"), "./", tr("Vinay files (*.vinay);;Ctable files( *.ctbl)") );
          dialog.setFileMode( QFileDialog::AnyFile );
          dialog.setViewMode(QFileDialog::Detail);
	  dialog.setAcceptMode(QFileDialog::AcceptSave);

          QString filter;
          //filename = dialog.getSaveFileName( this, tr("Save Color Map"), "", tr("Vinay files (*.vinay);;Ctable files( *.ctbl)") );
          //filter = dialog.selectedNameFilter();
          
          QStringList filenames; 
          if( dialog.exec() ) {
             filter = dialog.selectedNameFilter();
             filename = dialog.selectedFiles().at(0);
          }
         
	  if(filename.isNull()) return;	 

          if( filter.compare( QString("Vinay files (*.vinay)") ) == 0 )
	    ColorTable::write_transfer_function(filename.toStdString(),info());
          else if( filter.compare( QString("Ctable files( *.ctbl)") ) == 0 )
	    ColorTable::write_full_color_table(filename.toStdString(), info());
          else
             std::cout << "Error: wrong file type: " << filter.toStdString() << std::endl;
	}
	break;
      case ADD_COLOR:
	{
	  ColorTable::color_nodes::iterator prev_iter, next_iter;
	  next_iter = info().colorNodes().upper_bound(ColorTable::color_node(newpos));
	  if(next_iter == info().colorNodes().end())
	    {
	      qDebug("ColorTable::Table::contextMenuEvent(): "
		     "WARNING: attempted to add color node past the end node, ignoring");
	      return;
	    }

	  if(next_iter == info().colorNodes().begin())
	    {
	      qDebug("ColorTable::Table::contextMenuEvent(): "
		     "WARNING: attempted to add color node before the start node, ignoring");
	      return;
	    }

	  prev_iter = prior(next_iter);

	  //linearly interpolate to get color between nodes at new position
	  ColorTable::color_node prev_node = *prev_iter;
	  ColorTable::color_node next_node = *next_iter;
	  double factor = (newpos - prev_node.position)/
	    (next_node.position - prev_node.position);
	  info().colorNodes().insert(
	    ColorTable::color_node(newpos,
				   prev_node.r + factor*(next_node.r - prev_node.r),
				   prev_node.g + factor*(next_node.g - prev_node.g),
				   prev_node.b + factor*(next_node.b - prev_node.b))
	    );

	  modified = true;
	}
	break;
      case ADD_ISOCONTOUR:
	{
          int _id = info().isocontourNodes().size();
	  info().isocontourNodes().insert(ColorTable::isocontour_node(newpos, _id));

          //Update contour spectrum/tree info. dialog
          updateInformDialog( _id, newpos, NEW );
 
	  modified = true;
	}
	break;
      case ADD_ALPHA:
	{
	  info().opacityNodes().insert(ColorTable::opacity_node(newpos,newopac));
	  modified = true;
	}
	break;
      case DELETE_SELECTION:
	{
	  if(_selectedObj < 1)
	    {
	      qDebug("ColorTable::Table::contextMenuEvent(): "
		     "WARNING: no object selected for deletion");
	      return;
	    }

	  if(_nameMap[_selectedObj-1].type() == 
	     typeid(ColorTable::isocontour_nodes::const_iterator)) {
            ColorTable::isocontour_nodes::const_iterator node_itr = 
	      any_cast<ColorTable::isocontour_nodes::const_iterator>(_nameMap[_selectedObj-1]);

	    ColorTable::isocontour_node node = *node_itr;

            //Update contour spectrum/tree info. dialog
            updateInformDialog( node.id(), 0.0, REMOVE );

	    info().isocontourNodes().erase(node_itr);
          }
	  else if(_nameMap[_selectedObj-1].type() == 
		  typeid(ColorTable::color_nodes::const_iterator))
	    {
	      ColorTable::color_nodes::const_iterator node_itr =
		any_cast<ColorTable::color_nodes::const_iterator>(_nameMap[_selectedObj-1]);
	      if(node_itr != info().colorNodes().begin() && //cannot delete first and last
		 node_itr != prior(info().colorNodes().end()))
		info().colorNodes().erase(node_itr);
	      else
		qDebug("ColorTable::Table::contextMenuEvent(): "
		       "WARNING: cannot delete first and last color nodes!");
	    }
	  else if(_nameMap[_selectedObj-1].type() ==
		  typeid(ColorTable::opacity_nodes::const_iterator))
	    {
	      ColorTable::opacity_nodes::const_iterator node_itr =
		any_cast<ColorTable::opacity_nodes::const_iterator>(_nameMap[_selectedObj-1]);
	      if(node_itr != info().opacityNodes().begin() && //cannot delete first and last
		 node_itr != prior(info().opacityNodes().end()))
		info().opacityNodes().erase(node_itr);
	      else
		qDebug("ColorTable::Table::contextMenuEvent(): "
		       "WARNING: cannot delete first and last opacity nodes!");
	    }

	  modified = true;
	}
	break;
      case EDIT_SELECTION:
	{
	  if(_selectedObj < 1)
	    {
	      qDebug("ColorTable::Table::contextMenuEvent(): "
		     "WARNING: no object selected for editing.");
	      return;
	    }

	  if(_nameMap[_selectedObj-1].type() == 
		  typeid(ColorTable::color_nodes::const_iterator))
	    {
	      ColorTable::color_nodes::const_iterator node_itr =
		any_cast<ColorTable::color_nodes::const_iterator>(_nameMap[_selectedObj-1]);

	      ColorTable::color_node node = *node_itr;
	      QColor color(int(node.r * 255.0),
			   int(node.g * 255.0),
			   int(node.b * 255.0));
	      QColor c = QColorDialog::getColor(color);
	      node.r = double(c.red())/255.0;
	      node.g = double(c.green())/255.0;
	      node.b = double(c.blue())/255.0;
	      info().colorNodes().erase(node_itr);
	      info().colorNodes().insert(node);

	      modified = true;
	    }
	  else
	    qDebug("ColorTable::Table::contextMenuEvent(): "
		   "WARNING: can only edit color nodes in this manner.");
	}
	break;
      case RESET:
	{
	  info() = ColorTable::default_transfer_function();
	  modified = true;
	}
	break;
      }

    if(modified)
      {
	update();
	emit changed();
      }
    
    _selectedObj = -1; //de-select

    e->accept();
  }

  Table::POPUPSELECTION Table::showPopup(QPoint point)
  {
#if QT_VERSION < 0x040000
    QPopupMenu popup;
    QPopupMenu add;
    QPopupMenu display;
    
    add.insertItem( "Alpha Node", ADD_ALPHA);
    add.insertItem( "Color Node", ADD_COLOR);
    add.insertItem( "Isocontour Node", ADD_ISOCONTOUR);
    
    display.insertItem( "Contour Spectrum", DISP_CONTOUR_SPECTRUM);
    display.setItemChecked(DISP_CONTOUR_SPECTRUM, /*m_DrawContourSpec*/ false);
    display.insertItem( "Contour Tree", DISP_CONTOUR_TREE);
    display.setItemChecked(DISP_CONTOUR_TREE, /*m_DrawContourTree*/ false);
    display.insertItem( "Contour Information", DISP_INFORM_DIALOG);
    display.setItemChecked(DISP_INFORM_DIALOG, /*m_DrawInfomDialog*/ false);
    display.insertItem( "Opacity Function", DISP_ALPHA_MAP);
    display.setItemChecked(DISP_ALPHA_MAP, /*m_DrawAlphaMap*/ false);
    display.insertItem( "Color Function", DISP_TRANS_MAP);
    display.setItemChecked(DISP_TRANS_MAP, /*m_DrawTransMap*/ false);
    display.insertItem( "Histogram", DISP_HISTOGRAM);
    display.setItemChecked(DISP_HISTOGRAM, /*m_DrawHistogram*/ false);
    
    popup.insertItem( "Open", LOAD_MAP );
    popup.insertItem( "Save", SAVE_MAP );
    popup.insertItem( "Add", &add, ADD_MENU );
    popup.insertItem( "Display", &display, DISP_MENU);
    popup.insertItem( "Edit", EDIT_SELECTION);
    popup.insertItem( "Delete", DELETE_SELECTION );
    popup.insertItem( "Reset", RESET);
    
    return (POPUPSELECTION)popup.exec(point);
#else
    QMenu popup;
    QMenu add("Add");
    QMenu display("Display");
    std::map<QAction*,POPUPSELECTION> actionmap;
    QAction *action = NULL;

    actionmap[add.addAction("Alpha Node")] = ADD_ALPHA;
    actionmap[add.addAction("Color Node")] = ADD_COLOR;
    actionmap[add.addAction("Isocontour Node")] = ADD_ISOCONTOUR;

#ifndef COLORTABLE2_DISABLE_CONTOUR_SPECTRUM
    actionmap[ action = display.addAction("Contour Spectrum") ] = DISP_CONTOUR_SPECTRUM;
    action->setCheckable(true);
    action->setChecked(visibleComponents() & CONTOUR_SPECTRUM);
    connect(action,SIGNAL(triggered(bool)),SLOT(showContourSpectrum(bool)));
#endif
#ifndef COLORTABLE2_DISABLE_CONTOUR_TREE
    actionmap[ action = display.addAction("Contour Tree") ] = DISP_CONTOUR_TREE;
    action->setCheckable(true);
    action->setChecked(visibleComponents() & CONTOUR_TREE);
    connect(action,SIGNAL(triggered(bool)),SLOT(showContourTree(bool)));
#endif
    actionmap[ action = display.addAction("Contour Information") ] = DISP_INFORM_DIALOG;
    action->setCheckable(true);
    action->setChecked(false);
    connect(action,SIGNAL(triggered(bool)),SLOT(showInformDialog(bool)));

    actionmap[ action = display.addAction("Opacity Function") ] = DISP_ALPHA_MAP;
    action->setCheckable(true);
    action->setChecked(visibleComponents() & OPACITY_NODES);
    connect(action,SIGNAL(triggered(bool)),SLOT(showOpacityFunction(bool)));

    actionmap[ action = display.addAction("Transfer Function") ] = DISP_TRANS_MAP;
    action->setCheckable(true);
    action->setChecked(visibleComponents() & BACKGROUND);
    connect(action,SIGNAL(triggered(bool)),SLOT(showTransferFunction(bool)));

    actionmap[ action = display.addAction("Histogram") ] = DISP_HISTOGRAM;
    action->setCheckable(true);
    action->setChecked(visibleComponents() & HISTOGRAM);
    connect(action,SIGNAL(triggered(bool)),SLOT(showHistogram(bool)));
    
    actionmap[popup.addAction("Open")] = LOAD_MAP;
    actionmap[popup.addAction("Save")] = SAVE_MAP;
    actionmap[popup.addMenu(&add)] = ADD_MENU;
    actionmap[popup.addMenu(&display)] = DISP_MENU;
    actionmap[popup.addAction("Edit")] = EDIT_SELECTION;
    actionmap[popup.addAction("Delete")] = DELETE_SELECTION;
    actionmap[popup.addAction("Reset")] = RESET;

    return actionmap[popup.exec(point)];
#endif
  }

  void Table::updateInformDialog( const int _id, const double _newpos, CONTOURSTATUS _status )
  {
     if( (_newpos < 0.0 ) || (_newpos > 1.0 ) ) return;

     switch( _status ) {
     case SHOW:
      {
        // if it's not computed yet, compute.
        computeContourSpectrum();
        computeContourTree();

        // update information if there are isocontours already
        for(ColorTable::isocontour_nodes::const_iterator i = _cti.isocontourNodes().begin();
            i != _cti.isocontourNodes().end();
            i++)
          {
             float isoval, area, minvol, maxvol, grad;
             int ncomp;
             computeInformation( i->position, &isoval, &area, &minvol, &maxvol, &grad, &ncomp);
             m_DrawInformDialog->updateInfo( i->_id, i->position, isoval, area, minvol, maxvol, grad, ncomp );
          }
      }
     break;
     case NEW:
     case MOVE:
      {
         float isoval, area, minvol, maxvol, grad;
         int ncomp;
         computeInformation( _newpos, &isoval, &area, &minvol, &maxvol, &grad, &ncomp);
         m_DrawInformDialog->updateInfo( _id, _newpos, isoval, area, minvol, maxvol, grad, ncomp );
      }
     break;
     case REMOVE:
        m_DrawInformDialog->remove( _id );
     break;
     }
  }

  void Table::computeInformation( const double pos, float *isoval, float *area, float *minvol, float *maxvol, float *grad, int *ncomp)
  {
     *isoval = 0.0f;
     *area = 0.0f;
     *minvol = 0.0f;
     *maxvol = 0.0f;
     *grad = 0.0f;
     *ncomp = 0;

     if( !_dirtyContourSpectrum )
     {
        // compute information
        double tpos = pos * 256.0;
        int i  = (int)(tpos);
        int pi = ((i+1)>255)?255:(i+1);
        double dx = tpos - i;
        double mdx = 1.0 - dx;
   
        float val[10] = {0.0};
   
        val[0] = _contourSpectrum["isovalue"].second[i];
        val[1] = _contourSpectrum["isovalue"].second[pi];
        val[2] = _contourSpectrum["area"].second[i];
        val[3] = _contourSpectrum["area"].second[pi];
        val[4] = _contourSpectrum["min_vol"].second[i];
        val[5] = _contourSpectrum["min_vol"].second[pi];
        val[6] = _contourSpectrum["max_vol"].second[i];
        val[7] = _contourSpectrum["max_vol"].second[pi];
        val[8] = _contourSpectrum["gradient"].second[i];
        val[9] = _contourSpectrum["gradient"].second[pi];
   
        for( int j = 0;  j < 5; j++ )
           val[j*2] = val[j*2]*mdx + val[j*2+1]*dx;
   
        // map iso value
        *isoval = val[0] * (_contourVolume.max() - _contourVolume.min()) + _contourVolume.min();
        *area = val[2];
        *minvol = val[4];
        *maxvol = val[6];
        *grad = val[8];
     }
 
     if( !_dirtyContourTree )
     {
        int ncomponents = 0; 
        for(std::vector<int>::iterator i = _contourTreeEdges.begin();
            i != _contourTreeEdges.end();
            i++)
          {
            if(boost::next(i) == _contourTreeEdges.end()) break;
   
            int v1 = *i;
            int v2 = *boost::next(i);
   
            double p1x,p1y, p2x,p2y;
            
            // get the vertex coords
            p1x = _contourTreeVertices[v1 * 2 + 0]; //func_val
            //p1y = _contourTreeVertices[v1 * 2 + 1]; //norm_x
            p2x = _contourTreeVertices[v2 * 2 + 0]; //func_val
            //p2y = _contourTreeVertices[v2 * 2 + 1]; //norm_x;
   
            // don't add edges that are outside of the zoomed in region
            if (p1x < _min && p2x < _min) continue;
            if (p1x > _max && p2x > _max) continue;
            if ( ( p1x < pos ) && ( p2x > pos ) )
               ncomponents++;
         }
       *ncomp = ncomponents;
     }
  }
}
