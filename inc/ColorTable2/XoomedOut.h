/*
  Copyright 2002-2003 The University of Texas at Austin
  
    Authors: Anthony Thane <thanea@ices.utexas.edu>
             Vinay Siddavanahalli <skvinay@cs.utexas.edu>
	     Jose Rivera <transfix@ices.utexas.edu>
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

#ifndef __CVCCOLORTABLE__XOOMEDOUT_H__
#define __CVCCOLORTABLE__XOOMEDOUT_H__

#if QT_VERSION >= 0x040000
# include <QFrame>
# include <QPainter>
#else
# include <qframe.h>
# include <qpainter.h>
#endif

#include <ColorTable2/ColorTable.h>

namespace CVCColorTable
{
class XoomedOut : public QFrame 
{
  Q_OBJECT
public:
  XoomedOut(ColorTable::color_table_info& cti, QWidget *parent,
#if QT_VERSION < 0x040000
            const char *name=0
#else
            Qt::WindowFlags flags={}
#endif
            );

  virtual ~XoomedOut();

protected:
	virtual void drawContents( QPainter* painter );
///\fn virtual void drawContentsVertical( QPainter* painter )
///\brief Draws the widget in a vertical orientation
///\param painter A QPainter object to use for drawing
	virtual void drawContentsVertical( QPainter* painter );
///\fn virtual void drawContentsHorizontal( QPainter* painter )
///\brief Draws the widget in a horizontal orientation
///\param painter A QPainter object to use for drawing
	virtual void drawContentsHorizontal( QPainter* painter );

///\fn void drawColorMapVertical( QPainter* painter, QRect rect )
///\brief Draws the ColorMap in a vertical orientation
///\param painter A QPainter object to use for drawing
///\param rect The bounding rectangle
	void drawColorMapVertical( QPainter* painter, QRect rect );
///\fn void drawControlBarVertical( QPainter* painter, QRect rect, double doublePos )
///\brief Draws a control bar in a vertical orientation
///\param painter A QPainter object to use for drawing
///\param rect The bounding rectangle
///\param doublePos A position along the widget
	void drawControlBarVertical( QPainter* painter, QRect rect, double doublePos );

///\fn void drawColorMapHorizontal( QPainter* painter, QRect rect )
///\brief Draws the ColorMap in a horizontal orientation
///\param painter A QPainter object to use for drawing
///\param rect The bounding rectangle
	void drawColorMapHorizontal( QPainter* painter, QRect rect );
///\fn void drawControlBarHorizontal( QPainter* painter, QRect rect, double doublePos )
///\brief Draws a control bar in a horizontal orientation
///\param painter A QPainter object to use for drawing
///\param rect The bounding rectangle
///\param doublePos A position along the widget
	void drawControlBarHorizontal( QPainter* painter, QRect rect, double doublePos );

	ColorTable::color_table_info& m_ColorTableInformation;
	double m_RangeMin;
	double m_RangeMax;


	enum HANDLE {MIN_HANDLE, MAX_HANDLE, NO_HANDLE};

	HANDLE m_SelectedHandle ;

	void mouseMoveEvent( QMouseEvent* q );
	void mousePressEvent( QMouseEvent* q );
	void mouseReleaseEvent( QMouseEvent* q );

///\fn void mouseMoveEventVertical( QMouseEvent* q )
///\brief Handles mouse move events for a vertically oriented widget
	void mouseMoveEventVertical( QMouseEvent* q );
///\fn void mousePressEventVertical( QMouseEvent* q )
///\brief Handles mouse press events for a vertically oriented widget
	void mousePressEventVertical( QMouseEvent* q );
///\fn void mouseReleaseEventVertical( QMouseEvent* q )
///\brief Handles mouse release events for a vertically oriented widget
	void mouseReleaseEventVertical( QMouseEvent* q );

///\fn void mouseMoveEventHorizontal( QMouseEvent* q )
///\brief Handles mouse move events for a horizontally oriented widget
	void mouseMoveEventHorizontal( QMouseEvent* q );
///\fn void mousePressEventHorizontal( QMouseEvent* q )
///\brief Handles mouse press events for a horizontally oriented widget
	void mousePressEventHorizontal( QMouseEvent* q );
///\fn void mouseReleaseEventHorizontal( QMouseEvent* q )
///\brief Handles mouse release events for a horizontally oriented widget
	void mouseReleaseEventHorizontal( QMouseEvent* q );

///\fn int getIntPosVertical(double doublePos, QRect rect)
///\brief Converts from a double in the range [0,1] to an int
///\param doublePos A position along the widget
///\param rect The bounding rectangle
///\return An int
	int getIntPosVertical(double doublePos, QRect rect);
///\fn double getDoublePosVertical(int intPos, QRect rect)
///\brief Converts from an int to a double in the range [0,1]
///\param intPos A position along the widget
///\param rect The bounding rectangle
///\return A double
	double getDoublePosVertical(int intPos, QRect rect);

///\fn int getIntPosHorizontal(double doublePos, QRect rect)
///\brief Converts from a double in the range [0,1] to an int
///\param doublePos A position along the widget
///\param rect The bounding rectangle
///\return An int
	int getIntPosHorizontal(double doublePos, QRect rect);
///\fn double getDoublePosHorizontal(int intPos, QRect rect)
///\brief Converts from an int to a double in the range [0,1]
///\param intPos A position along the widget
///\param rect The bounding rectangle
///\return A double
	double getDoublePosHorizontal(int intPos, QRect rect);

signals:
///\fn void minChanged( double value )
///\brief Signals when the minimum range density has changed
	void minChanged( double value );
///\fn void minExploring( double value )
///\brief Signals while the minimum range density is changing
	void minExploring( double value );
///\fn void maxChanged( double value )
///\brief Signals when the maximum range density has changed
	void maxChanged( double value );
///\fn void maxExploring( double value )
///\brief Signals while the maximum range density is changing
	void maxExploring( double value );
};
}
#endif
