/*
  Copyright 2002-2003 The University of Texas at Austin
  
  Authors: Anthony Thane <thanea@ices.utexas.edu>
  Vinay Siddavanahalli <skvinay@cs.utexas.edu>
  Advisor: Chandrajit Bajaj <bajaj@cs.utexas.edu>

  This file is part of Volume Rover.

  Volume Rover is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License version 2.1 as published by the Free Software Foundation.

  Volume Rover is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/

#include <ColorTable2/XoomedOut.h>
#include <cstdlib>
#include <boost/utility.hpp>
#include <boost/next_prior.hpp>

#if QT_VERSION < 0x040000
#include <qpixmap.h>
#else
#include <QPixmap>
#endif

#if QT_VERSION > 0x040000
#include <QMouseEvent>
#endif

namespace CVCColorTable
{
  XoomedOut::XoomedOut( ColorTable::color_table_info& cti, QWidget *parent, 
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
              ), m_ColorTableInformation(cti)
  {
    setFrameStyle( QFrame::Panel | QFrame::Sunken );
    m_RangeMin = MIN_RANGE;
    m_RangeMax = MAX_RANGE;
#if QT_VERSION < 0x040000
    setBackgroundMode( Qt::NoBackground );
#else
    setAttribute(Qt::WA_OpaquePaintEvent);
#endif
  }

  XoomedOut::~XoomedOut()
  {
  }

  static int MapToPixel(double input, int start, int end)
  {
    if (input <= 1.0 && input >= 0.0) {
      int width = end-start;
      int offset = (int)(input * width);
      return start + offset;
    }
    else if (input > 1.0) {
      return end;
    }
    else {
      return start;
    }
  }

  static double MapToDouble(int input, int start, int end)
  {
    int width = end-start;
    int offset = input-start;

    return (double)offset/(double)width;
  }


  void XoomedOut::drawContents( QPainter* painter )
  {
    drawContentsHorizontal(painter);
  }

  void XoomedOut::drawContentsVertical( QPainter* painter )
  {
    Q_UNUSED(painter);

    QRect rect = contentsRect();
    QRect reducedRect(rect.x(), rect.y()+5, rect.width(), rect.height()-10);
	
    QPixmap  pix( rect.size() );       // Pixmap for double-buffering

#if QT_VERSION < 0x040000
    QColor test = paletteBackgroundColor();
#else
    QColor test = palette().color(backgroundRole());
#endif

    //pix.setBrush(backgroundBursh());
    pix.fill( test ); //this, rect.topLeft() );  // fill with widget background

    QPainter p( &pix );

    drawColorMapVertical(&p, reducedRect);
    drawControlBarVertical(&p, reducedRect, m_RangeMin);
    drawControlBarVertical(&p, reducedRect, m_RangeMax);

#if QT_VERSION < 0x040000
    bitBlt( this, rect.topLeft(), &pix );
#else
    QPainter widgetPainter(this);
    widgetPainter.drawPixmap(rect.topLeft(), pix);
#endif
    p.end();
  }

  void XoomedOut::drawContentsHorizontal( QPainter* painter )
  {
    Q_UNUSED(painter);

    QRect rect = contentsRect();
    QRect reducedRect(rect.x()+5, rect.y(), rect.width()-10, rect.height());
	
    QPixmap  pix( rect.size() );       // Pixmap for double-buffering

#if QT_VERSION < 0x040000
    QColor test = paletteBackgroundColor();
#else
    QColor test = palette().color(backgroundRole());
#endif

    //pix.setBrush(backgroundBursh());
    pix.fill( test ); //this, rect.topLeft() );  // fill with widget background

    QPainter p( &pix );

    drawColorMapHorizontal(&p, reducedRect);
    drawControlBarHorizontal(&p, reducedRect, m_RangeMin);
    drawControlBarHorizontal(&p, reducedRect, m_RangeMax);

#if QT_VERSION < 0x040000
    bitBlt( this, rect.topLeft(), &pix );
#else
    QPainter widgetPainter(this);
    widgetPainter.drawPixmap(rect.topLeft(), pix);
#endif
    p.end();
  }

  void XoomedOut::drawColorMapVertical( QPainter* painter, QRect rect )
  {
    using namespace boost;

    int miny = rect.y();
    int maxy = rect.y() + rect.height() - 1;
    int minx = rect.x();
    int maxx = rect.x() + rect.width() - 1;
    double mag;
    int pos1, pos2;
    double cr1, cr2;
    double cg1, cg2;
    double cb1, cb2;

    //int size=rect.height();

    //if we dont have enough color nodes, just pretend we have 2 nodes for grayscale
    if(m_ColorTableInformation.colorNodes().size() < 2)
      {
	pos1 = MapToPixel(MIN_RANGE, miny, maxy);
	pos2 = MapToPixel(MAX_RANGE, miny, maxy);
	cr1 = 0.0; cr2 = 1.0;
	cg1 = 0.0; cg2 = 1.0;
	cb1 = 0.0; cb2 = 1.0;
	for (int y=pos1; y<=pos2; y++) {
	  mag = MapToDouble(y, pos1, pos2);
	  /*pMap[y*4] = cr1*(1.0f-mag) + cr2*(mag);
	    pMap[y*4+1] = cg1*(1.0f-mag) + cg2*(mag);
	    pMap[y*4+2] = cb1*(1.0f-mag) + cb2*(mag);*/
	  painter->setPen(QColor(MapToPixel(cr1*(1.0f-mag) + cr2*(mag), 0, 255), 
				 MapToPixel(cg1*(1.0f-mag) + cg2*(mag), 0, 255), 
				 MapToPixel(cb1*(1.0f-mag) + cb2*(mag), 0, 255)));
	  painter->drawLine(minx, y,maxx,y);
	}
      }
    else
      {
	for(ColorTable::color_nodes::const_iterator i = m_ColorTableInformation.colorNodes().begin();
	    i != boost::prior(m_ColorTableInformation.colorNodes().end());
	    i++)
	  {
	    pos1 = MapToPixel(i->position, miny, maxy);
	    pos2 = MapToPixel(std::next(i)->position, miny, maxy);
	    cr1 = i->r;
	    cr2 = std::next(i)->r;
	    cg1 = i->g;
	    cg2 = std::next(i)->g;
	    cb1 = i->b;
	    cb2 = std::next(i)->b;
	    for (int y=pos1; y<=pos2; y++)
	      {
		mag = MapToDouble(y, pos1, pos2);
		/*pMap[y*4] = cr1*(1.0f-mag) + cr2*(mag);
		  pMap[y*4+1] = cg1*(1.0f-mag) + cg2*(mag);
		  pMap[y*4+2] = cb1*(1.0f-mag) + cb2*(mag);*/
		painter->setPen(QColor(MapToPixel(cr1*(1.0f-mag) + cr2*(mag), 0, 255), 
				       MapToPixel(cg1*(1.0f-mag) + cg2*(mag), 0, 255), 
				       MapToPixel(cb1*(1.0f-mag) + cb2*(mag), 0, 255)));
		painter->drawLine(minx, y,maxx,y);
	      }
	  }
      }
  }

  void XoomedOut::drawControlBarVertical( QPainter* painter, QRect rect, double doublePos )
  {
    int pos;
    int miny = rect.y();
    int maxy = rect.y() + rect.height() - 1;
    int minx = rect.x();
    int maxx = rect.x() + rect.width() - 1;
    pos = MapToPixel(doublePos, miny, maxy);

    painter->setPen(QColor(255, 255, 255));
    painter->drawLine(minx, pos-1,maxx,pos-1);
    painter->setPen(QColor(192, 192, 192));
    painter->drawLine(minx, pos,maxx,pos);
    painter->setPen(QColor(128, 128, 128));
    painter->drawLine(minx, pos+1,maxx,pos+1);

  }

  void XoomedOut::drawColorMapHorizontal( QPainter* painter, QRect rect )
  {
    using namespace boost;

    int miny = rect.y();
    int maxy = rect.y() + rect.height() - 1;
    int minx = rect.x();
    int maxx = rect.x() + rect.width() - 1;
    double mag;
    int pos1, pos2;
    double cr1, cr2;
    double cg1, cg2;
    double cb1, cb2;

    //int size=rect.height();

    //if we dont have enough color nodes, just pretend we have 2 nodes for grayscale
    if(m_ColorTableInformation.colorNodes().size() < 2)
      {
	pos1 = MapToPixel(MIN_RANGE, minx, maxx);
	pos2 = MapToPixel(MAX_RANGE, minx, maxx);
	cr1 = 0.0; cr2 = 1.0;
	cg1 = 0.0; cg2 = 1.0;
	cb1 = 0.0; cb2 = 1.0;
	for (int x=pos1; x<=pos2; x++) {
	  mag = MapToDouble(x, pos1, pos2);
	  /*pMap[y*4] = cr1*(1.0f-mag) + cr2*(mag);
	    pMap[y*4+1] = cg1*(1.0f-mag) + cg2*(mag);
	    pMap[y*4+2] = cb1*(1.0f-mag) + cb2*(mag);*/
	  painter->setPen(QColor(MapToPixel(cr1*(1.0f-mag) + cr2*(mag), 0, 255), 
				 MapToPixel(cg1*(1.0f-mag) + cg2*(mag), 0, 255), 
				 MapToPixel(cb1*(1.0f-mag) + cb2*(mag), 0, 255)));
	  painter->drawLine(x, miny,x,maxy);
	}
      }
    else
      {
	for(ColorTable::color_nodes::const_iterator i = m_ColorTableInformation.colorNodes().begin();
	    i != boost::prior(m_ColorTableInformation.colorNodes().end());
	    i++)
	  {
	    pos1 = MapToPixel(i->position, minx, maxx);
	    pos2 = MapToPixel(std::next(i)->position, minx, maxx);
	    cr1 = i->r;
	    cr2 = std::next(i)->r;
	    cg1 = i->g;
	    cg2 = std::next(i)->g;
	    cb1 = i->b;
	    cb2 = std::next(i)->b;
	    for (int x=pos1; x<=pos2; x++)
	      {
		mag = MapToDouble(x, pos1, pos2);
		/*pMap[y*4] = cr1*(1.0f-mag) + cr2*(mag);
		  pMap[y*4+1] = cg1*(1.0f-mag) + cg2*(mag);
		  pMap[y*4+2] = cb1*(1.0f-mag) + cb2*(mag);*/
		painter->setPen(QColor(MapToPixel(cr1*(1.0f-mag) + cr2*(mag), 0, 255), 
				       MapToPixel(cg1*(1.0f-mag) + cg2*(mag), 0, 255), 
				       MapToPixel(cb1*(1.0f-mag) + cb2*(mag), 0, 255)));
		painter->drawLine(x, miny,x,maxy);
	      }
	  }
      }
  }

  void XoomedOut::drawControlBarHorizontal( QPainter* painter, QRect rect, double doublePos )
  {
    int pos;
    int miny = rect.y();
    int maxy = rect.y() + rect.height() - 1;
    int minx = rect.x();
    int maxx = rect.x() + rect.width() - 1;
    pos = MapToPixel(doublePos, minx, maxx);

    painter->setPen(QColor(255, 255, 255));
    painter->drawLine(pos-1, miny,pos-1,maxy);
    painter->setPen(QColor(192, 192, 192));
    painter->drawLine(pos, miny,pos,maxy);
    painter->setPen(QColor(128, 128, 128));
    painter->drawLine(pos+1, miny,pos+1,maxy);
  }

  void XoomedOut::mouseMoveEvent( QMouseEvent* q )
  {
    mouseMoveEventHorizontal(q);
  }

  void XoomedOut::mousePressEvent( QMouseEvent* q )
  {
    mousePressEventHorizontal(q);
  }

  void XoomedOut::mouseReleaseEvent( QMouseEvent* q )
  {
    mouseReleaseEventHorizontal(q);
  }

  void XoomedOut::mouseMoveEventVertical( QMouseEvent* q )
  {
    //int x = q->x();
    int y = qRound(q->position().y());
    QRect rect = contentsRect();
    QRect reducedRect(rect.x(), rect.y()+5, rect.width(), rect.height()-10);

    if( m_SelectedHandle == MIN_HANDLE ) {
      m_RangeMin = getDoublePosVertical(y, reducedRect);
      if (m_RangeMin > MAX_RANGE) m_RangeMin = MAX_RANGE;
      if (m_RangeMin < MIN_RANGE) m_RangeMin = MIN_RANGE;
      if (getIntPosVertical(m_RangeMax, reducedRect) - getIntPosVertical(m_RangeMin, reducedRect) <=3) {
	m_RangeMin = getDoublePosVertical(getIntPosVertical(m_RangeMax, reducedRect)-3, reducedRect);
      }

      emit minExploring( m_RangeMin );
    }
    if( m_SelectedHandle == MAX_HANDLE ) {
      m_RangeMax = getDoublePosVertical(y, reducedRect);
      if (m_RangeMax > MAX_RANGE) m_RangeMax = MAX_RANGE;
      if (m_RangeMax < MIN_RANGE) m_RangeMax = MIN_RANGE;
      if (getIntPosVertical(m_RangeMax, reducedRect) - getIntPosVertical(m_RangeMin, reducedRect) <=3) {
	m_RangeMax = getDoublePosVertical(getIntPosVertical(m_RangeMin, reducedRect)+3, reducedRect);
      }

      emit maxExploring( m_RangeMax );
    }

    repaint();
  }

  void XoomedOut::mousePressEventVertical( QMouseEvent* q )
  {
    //int x = q->x();
    int y = qRound(q->position().y());
    QRect rect = contentsRect();
    QRect reducedRect(rect.x(), rect.y()+5, rect.width(), rect.height()-10);
    int intMin, intMax;

    intMin = getIntPosVertical( m_RangeMin, reducedRect );
    intMax = getIntPosVertical( m_RangeMax, reducedRect );
    if( ((y <= intMin + 2) && (y >= intMin-2)) && ((y <= intMax + 2) && (y >= intMax - 2))) {
      m_SelectedHandle = ( abs(y-intMin) <= abs(y-intMax)? MIN_HANDLE: MAX_HANDLE  );
    }
    else if( (y <= intMin + 2) && (y >= intMin-2) ) {
      m_SelectedHandle = MIN_HANDLE;
    }
    else if( (y <= intMax + 2) && (y >= intMax - 2) ) {
      m_SelectedHandle = MAX_HANDLE;
    }
    else {
      m_SelectedHandle = NO_HANDLE;
    }

    update();
  }

  void XoomedOut::mouseReleaseEventVertical( QMouseEvent* q )
  {
    Q_UNUSED(q);
    //int x = q->x();
    //int y = q->y();
    update();
  }


  void XoomedOut::mouseMoveEventHorizontal( QMouseEvent* q )
  {
    int x = qRound(q->position().x());
    //int y = q->y();
    QRect rect = contentsRect();
    QRect reducedRect(rect.x()+5, rect.y(), rect.width()-10, rect.height());

    if( m_SelectedHandle == MIN_HANDLE ) {
      m_RangeMin = getDoublePosHorizontal(x, reducedRect);
      if (m_RangeMin > MAX_RANGE) m_RangeMin = MAX_RANGE;
      if (m_RangeMin < MIN_RANGE) m_RangeMin = MIN_RANGE;
      if (getIntPosHorizontal(m_RangeMax, reducedRect) - getIntPosHorizontal(m_RangeMin, reducedRect) <=3) {
	m_RangeMin = getDoublePosHorizontal(getIntPosHorizontal(m_RangeMax, reducedRect)-3, reducedRect);
      }

      emit minExploring( m_RangeMin );
    }
    if( m_SelectedHandle == MAX_HANDLE ) {
      m_RangeMax = getDoublePosHorizontal(x, reducedRect);
      if (m_RangeMax > MAX_RANGE) m_RangeMax = MAX_RANGE;
      if (m_RangeMax < MIN_RANGE) m_RangeMax = MIN_RANGE;
      if (getIntPosHorizontal(m_RangeMax, reducedRect) - getIntPosHorizontal(m_RangeMin, reducedRect) <=3) {
	m_RangeMax = getDoublePosHorizontal(getIntPosHorizontal(m_RangeMin, reducedRect)+3, reducedRect);
      }

      emit maxExploring( m_RangeMax );
    }

    repaint();
  }

  void XoomedOut::mousePressEventHorizontal( QMouseEvent* q )
  {
    int x = qRound(q->position().x());
    //int y = q->y();
    QRect rect = contentsRect();
    QRect reducedRect(rect.x()+5, rect.y(), rect.width()-10, rect.height());
    int intMin, intMax;

    intMin = getIntPosHorizontal( m_RangeMin, reducedRect );
    intMax = getIntPosHorizontal( m_RangeMax, reducedRect );
    if( ((x <= intMin + 2) && (x >= intMin-2)) && ((x <= intMax + 2) && (x >= intMax - 2))) {
      m_SelectedHandle = ( abs(x-intMin) <= abs(x-intMax)? MIN_HANDLE: MAX_HANDLE  );
    }
    else if( (x <= intMin + 2) && (x >= intMin-2) ) {
      m_SelectedHandle = MIN_HANDLE;
    }
    else if( (x <= intMax + 2) && (x >= intMax - 2) ) {
      m_SelectedHandle = MAX_HANDLE;
    }
    else {
      m_SelectedHandle = NO_HANDLE;
    }

    update();
  }

  void XoomedOut::mouseReleaseEventHorizontal( QMouseEvent* q )
  {
    Q_UNUSED(q);
    //int x = q->x();
    //int y = q->y();
    update();
  }

  int XoomedOut::getIntPosVertical(double doublePos, QRect rect)
  {
    return (int)(doublePos*(rect.height()-1))+rect.y();
  }

  double XoomedOut::getDoublePosVertical(int intPos, QRect rect)
  {
    return (double)(intPos-rect.y())/(double)(rect.height()-1);
  }

  int XoomedOut::getIntPosHorizontal(double doublePos, QRect rect)
  {
    return (int)(doublePos*(rect.width()-1))+rect.x();
  }

  double XoomedOut::getDoublePosHorizontal(int intPos, QRect rect)
  {
    return (double)(intPos-rect.x())/(double)(rect.width()-1);
  }

}
