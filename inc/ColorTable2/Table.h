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

#ifndef __CVCCOLORTABLE__TABLE_H__
#define __CVCCOLORTABLE__TABLE_H__

#define NOMINMAX
#include <QGLViewer/qglviewer.h>
#include <ColorTable2/ColorTable.h>
#include <ColorTable2/InfoDialog.h>

#include <vector>
#include <map>
#include <boost/scoped_ptr.hpp>
#include <boost/any.hpp>
#include <boost/cstdint.hpp>
#include <boost/array.hpp>
#include <boost/tuple/tuple.hpp>

#if !defined(COLORTABLE2_DISABLE_CONTOUR_TREE) || !defined(COLORTABLE_DISABLE_CONTOUR_SPECTRUM)    
#include <VolMagick/VolMagick.h>
#endif

namespace CVCColorTable
{
  class Table : public QGLViewer
  {
    Q_OBJECT
  public:
    Table( ColorTable::color_table_info& cti, QWidget *parent, 
#if QT_VERSION < 0x040000 || defined QT3_SUPPORT
              const char *name=0
#else
              Qt::WindowFlags flags={}
#endif
              );

    Table( boost::uint64_t components,
              ColorTable::color_table_info& cti, QWidget *parent, 
#if QT_VERSION < 0x040000 || defined QT3_SUPPORT
              const char *name=0
#else
              Qt::WindowFlags flags={}
#endif
              );

    virtual ~Table();

    double min() const { return _min; }
    double max() const { return _max; }

    double rangeMin() const { return _rangeMin; }
    void rangeMin(double val);
    double rangeMax() const { return _rangeMax; }
    void rangeMax(double val);

    const ColorTable::color_table_info& info() const { return _cti; }
    ColorTable::color_table_info& info() { return _cti; }

    bool interactiveUpdates() const { return _interactiveUpdates; }

    //The different components of the table
    static const boost::uint64_t BACKGROUND       = 1 << 0;
    static const boost::uint64_t COLOR_BARS       = 1 << 1;
    static const boost::uint64_t ISOCONTOUR_BARS  = 1 << 2;
    static const boost::uint64_t OPACITY_NODES    = 1 << 3;
    static const boost::uint64_t RANGE_BARS       = 1 << 4;
    static const boost::uint64_t CONTOUR_TREE     = 1 << 5;
    static const boost::uint64_t CONTOUR_SPECTRUM = 1 << 6;
    static const boost::uint64_t HISTOGRAM        = 1 << 7;
    
    boost::uint64_t visibleComponents() const { return _visibleComponents; }
    void visibleComponents(boost::uint64_t components);

#if !defined(COLORTABLE2_DISABLE_CONTOUR_TREE) || !defined(COLORTABLE_DISABLE_CONTOUR_SPECTRUM)    
    void setContourVolume(const VolMagick::Volume& vol);
#endif

  public slots:
    void interactiveUpdates(bool b);
    void setMin(double min);
    void setMax(double max);
    void showOpacityFunction(bool b);
    void showTransferFunction(bool b);
    void showContourTree(bool b);
    void showContourSpectrum(bool b);
    void showHistogram(bool b);
    void showInformDialog(bool b);

  signals:
    void changed();

    void rangeMinChanged(double);
    void rangeMaxChanged(double);

  protected:
    virtual void init();
    virtual void draw();
    virtual void drawWithNames();

    virtual void drawTable(bool withNames = false);
    void drawBar(double x_pos, double depth, const GLfloat *color_3f, GLint name = -1);

    virtual void beginSelection(const QPoint& point); //overloading this because we want to select with ortho proj
    virtual void postSelection(const QPoint& point);

    virtual void mousePressEvent(QMouseEvent* e);
    virtual void mouseMoveEvent(QMouseEvent* e);
    virtual void mouseReleaseEvent(QMouseEvent* e);
    virtual void contextMenuEvent(QContextMenuEvent* e);

    void computeContourTree();
    void computeContourSpectrum();
    void computeInformation( const double pos, float *isoval, float *area, float *minvol, float *maxvol, float *grad, int *nComp);

    void allocateInformDialg(void);
    void updateInformDialog( const int _id, const double _newpos, CONTOURSTATUS _status );

    // arand, 8-24, 2011: added DO_NOTHING as the default element in the list
    //                    this fixes a bug that caused alpha nodes to be inserted
    //                    when the users just tried to close the popup menu
    enum POPUPSELECTION { DO_NOTHING, ADD_ALPHA, ADD_COLOR, ADD_ISOCONTOUR,
			  DISP_CONTOUR_SPECTRUM, DISP_CONTOUR_TREE, DISP_INFORM_DIALOG,
                          DISP_ALPHA_MAP,
			  DISP_TRANS_MAP, DISP_HISTOGRAM,
			  DELETE_SELECTION,EDIT_SELECTION, SAVE_MAP, LOAD_MAP,
			  ADD_MENU, DISP_MENU, RESET };
    POPUPSELECTION showPopup(QPoint point);

    ColorTable::color_table_info& _cti;
    
    //density range to show
    double _min;
    double _max;

    boost::scoped_ptr<qglviewer::WorldConstraint> _constraint;

    int _selectedObj;
    QPoint _selectedPoint;

    std::vector<boost::any> _nameMap;

    bool _interactiveUpdates; //if this is false, we only emit on 'mouse up'

    //bitfield for enabling/disabling the rendered components of the table
    boost::uint64_t _visibleComponents;

    //Range bars min/max.  Used for setting the min/max of an accompanying
    //table.
    double _rangeMin;
    double _rangeMax;

#if !defined(COLORTABLE2_DISABLE_CONTOUR_TREE) || !defined(COLORTABLE2_DISABLE_CONTOUR_SPECTRUM)
    VolMagick::Volume  _contourVolume;
#endif

    // TODO: remove the need for ColorTable2 needing to know about contour trees and spectrums.
    // Use generic 2D geometry layer and 1D function layers instead, and have both be computed
    // from outside.

    //Contour tree pointers
    std::vector<float> _contourTreeVertices;
    std::vector<int>   _contourTreeEdges;
    bool _dirtyContourTree;

    //Contour spectrum functions
    typedef boost::array<double,3> color_t;
    typedef std::vector<float> function_t;
    typedef std::map<std::string, std::pair<color_t,function_t> > function_map;
    function_map _contourSpectrum;
    bool _dirtyContourSpectrum;

    //Histogram stuff
    boost::tuple<const VolMagick::uint64 *,VolMagick::uint64> _histogram;
    bool _dirtyHistogram;

    //Information dialog
    InfoDialog *m_DrawInformDialog;
  };
}

#endif
