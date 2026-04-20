/*
  Copyright 2012 The University of Texas at Austin

        Authors: Joe Rivera <transfix@ices.utexas.edu>
        Advisor: Chandrajit Bajaj <bajaj@cs.utexas.edu>

  This file is part of libCVC.

  libCVC is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License version 2.1 as published by the Free Software Foundation.

  libCVC is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/

/* $Id: GeometryRenderer.cpp 5883 2012-07-20 19:52:38Z transfix $ */

// GLEW must be included before any OpenGL headers
#include <GL/glew.h>

#include <GeometryRenderer/GeometryRenderer.h>

#include <cvcraw_geometry/cvcraw_geometry.h>
#include <cvcraw_geometry/cvcgeom.h>

#include <QGLViewer/frame.h>
#include <QGLViewer/camera.h>

namespace CVC_NAMESPACE
{
  void GeometryRenderer::defaultConstructor()
  {
    getState("rendering_mode")
      .value("filled")
      .comment("The overall geometry rendering mode: 'filled', 'wireframe', 'boundingbox'");
  }

  void GeometryRenderer::render(const std::string& sceneRoot)
  {
    cvcapp.log(2,str(boost::format("%s :: rendering: %s\n")
                     % BOOST_CURRENT_FUNCTION
                     % sceneRoot));

    //main scene rendering via renderState
    getState(sceneRoot).traverse(
      boost::bind(&GeometryRenderer::renderState, 
                  boost::ref(*this),
                  sceneRoot,
                  boost::placeholders::_1));
  }

  void GeometryRenderer::renderState(const std::string& root,
                                     const std::string& s)
  {
    //TODO: handle sub volume selector and names
    //TODO: VBO
    //TODO: texture mapping to make volume rendering easier

    cvcapp.log(4,
               str(
                   boost::format("%s :: state: %s\n")
                   % BOOST_CURRENT_FUNCTION
                   % s));

    cvcraw_geometry::cvcgeom_t geom;
    if(cvcstate(s).isData<cvcraw_geometry::cvcgeom_t>())
      {
        geom = cvcstate(s).data<cvcraw_geometry::cvcgeom_t>();
      }
    else if(cvcstate(s).isData<cvcraw_geometry::geometry_t>())
      {
        geom = cvcstate(s).data<cvcraw_geometry::geometry_t>();
      }
    else if(cvcstate(s).value() == "glPushMatrix")
      {
        cvcapp.log(4,
                   str(
                       boost::format("%s :: glPushMatrix %s\n")
                       % BOOST_CURRENT_FUNCTION
                       % s));

        glPushMatrix();
      }
    else if(cvcstate(s).value() == "glPopMatrix")
      {
        cvcapp.log(4,
                   str(
                       boost::format("%s :: glPopMatrix %s\n")
                       % BOOST_CURRENT_FUNCTION
                       % s));

        glPopMatrix();
      }
    else if(cvcstate(s).isData<qglviewer::Frame>())
      {
        cvcapp.log(4,
                   str(
                       boost::format("%s :: Frame %s\n")
                       % BOOST_CURRENT_FUNCTION
                       % s));

        qglviewer::Frame fr = cvcstate(s).data<qglviewer::Frame>();
        glMultMatrixd(fr.matrix());
      }
    else if(cvcstate(s).value() == "glPushAttrib")
      {
        cvcapp.log(4,
                   str(
                       boost::format("%s :: glPushAttrib %s\n")
                       % BOOST_CURRENT_FUNCTION
                       % s));

        //if no mask set, default to all bits
        //GLbitfield mask = GL_LIGHTING_BIT | GL_ENABLE_BIT | GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT;
        GLbitfield mask = GL_ALL_ATTRIB_BITS;
        if(cvcstate(s).isData<GLbitfield>())
          mask = cvcstate(s).data<GLbitfield>();
        glPushAttrib(mask);
      }
    else if(cvcstate(s).value() == "glPopAttrib")
      {
        cvcapp.log(4,
                   str(
                       boost::format("%s :: glPopAttrib %s\n")
                       % BOOST_CURRENT_FUNCTION
                       % s));

        glPopAttrib();
      }
    else if(cvcstate(s).isData<Callback>())
      {
        Callback c = cvcstate(s).data<Callback>();
        c();
      }
    else if(cvcstate(s).isData<BoundingBox>())
      {
        cvcapp.log(4,
                   str(
                       boost::format("%s :: bounding box %s\n")
                       % BOOST_CURRENT_FUNCTION
                       % s));

        doDrawBoundingBox(s);
      }

    //nothing to do if empty
    if(geom.empty()) return;

    //handle naming if a name id is provided
    bool hasName = cvcstate(s)("name").initialized();
    if(hasName) glPushName(cvcstate(s)("name").value<int>());

    //do draw geometry
#if 0
    {
      using namespace qglviewer;
      const Vec camPos = camera()->position();
      //const GLfloat pos[4] = {camPos[0]+70.0f,camPos[1]+50.0f,camPos[2]+100.0f,1.0};
      const GLfloat pos[4] = {camPos[0],camPos[1],camPos[2],1.0};

      // arand: old lights
      //GLfloat diffuseColor[] = {1.0f, 1.0f, 1.0f, 1.0f};
      //GLfloat specularColor[] = {1.0f, 1.0f, 1.0f, 1.0f};
      //GLfloat ambientColor[] = {0.0f, 0.0f, 0.0f, 1.0f};

      GLfloat diffuseColor[] = {0.90f, 0.90f, 0.90f, 1.0f};
      GLfloat specularColor[] = {0.60f, 0.60f, 0.60f, 1.0f};
      GLfloat ambientColor[] = {0.0f, 0.0f, 0.0f, 1.0f};

      // arand, 7-19-2011
      // I am not sure if this code setting up the lighting options 
      // actually does anything... I just tried to copy TexMol
      glLightfv(GL_LIGHT0, GL_POSITION, pos);      
      glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuseColor);
      glLightfv(GL_LIGHT0, GL_SPECULAR, specularColor);
      glLightfv(GL_LIGHT0, GL_AMBIENT, ambientColor);

      glEnable(GL_LIGHTING);
      glEnable(GL_DEPTH_TEST);
    
      glEnable(GL_LIGHT0);
      //glEnable(GL_LIGHT1);
      glEnable(GL_NORMALIZE);
    
      //// arand: added to render both sides of the surface... 4-12-2011
      glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
    }
#endif

    GLint params[2];
    //back up current setting
    glGetIntegerv(GL_POLYGON_MODE,params);


    glEnable(GL_LIGHTING);
    
    //make sure we have normals!
    if(geom.const_points().size() != geom.const_normals().size()) {
      geom.calculate_surf_normals();
    }

    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_NORMAL_ARRAY);

#if 0
    if(usingVBO())
      {
        glBindBufferARB(GL_ARRAY_BUFFER_ARB,scene_geom.vboArrayBufferID);
        glVertexPointer(3, GL_DOUBLE, 0, (GLvoid*)scene_geom.vboArrayOffsets[0]);
        glNormalPointer(GL_DOUBLE, 0, (GLvoid*)scene_geom.vboArrayOffsets[1]);
      }
    else
#endif
      {
        glVertexPointer(3, GL_DOUBLE, 0, &(geom.const_points()[0]));
        glNormalPointer(GL_DOUBLE, 0, &(geom.const_normals()[0]));
      }

    if(geom.const_colors().size() == geom.const_points().size())
      {
        glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
        glEnable(GL_COLOR_MATERIAL);
        glEnableClientState(GL_COLOR_ARRAY);
#if 0
        if(usingVBO())
          glColorPointer(3, GL_DOUBLE, 0, (GLvoid*)scene_geom.vboArrayOffsets[2]);
        else
#endif
          glColorPointer(3, GL_DOUBLE, 0, &(geom.const_colors()[0])); 
      } 
    else {
      glColor3f(1.0,1.0,1.0);
    }

    if(cvcstate(s)("render_mode").value()=="points")
      {
        glPushAttrib(GL_LIGHTING_BIT | GL_POINT_BIT);
        glDisable(GL_LIGHTING);
        glPointSize(2.0);
        glEnable(GL_POINT_SMOOTH);
        
        //no need for normals when point rendering
        glDisableClientState(GL_NORMAL_ARRAY);

#if 0
        if(usingVBO())
          glDrawArrays(GL_POINTS, 0, scene_geom.vboVertSize/3);
        else
#endif
          glDrawArrays(GL_POINTS, 0, geom.const_points().size());
        
        glPopAttrib();
      }
    else if(cvcstate(s)("render_mode").value()=="lines")
      {
        glLineWidth(cvcstate(s)("line_width").value<float>());
        glDisable(GL_LIGHTING);
        glEnable(GL_LINE_SMOOTH);
        //glHint(GL_LINE_SMOOTH,GL_NICEST); // arand: this doesn't seem to improve anything
        glDisableClientState(GL_NORMAL_ARRAY);
        
        glPolygonMode(GL_FRONT, GL_LINE);

#if 0        
        if(usingVBO())
          {
            
            glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER_ARB,
                            scene_geom.vboLineElementArrayBufferID);
            glDrawElements(GL_LINES, 
                           geom.const_lines().size()*2,
                           GL_UNSIGNED_INT, 0);
          }
        else
#endif
          glDrawElements(GL_LINES, 
                         geom.const_lines().size()*2,
                         GL_UNSIGNED_INT, &(geom.const_lines()[0]));
      }
    else if(cvcstate(s)("render_mode").value()=="triangles")
      {
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
#if 0
        if(usingVBO())
          {
            glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER_ARB,
                            scene_geom.vboTriElementArrayBufferID);
            glDrawElements(GL_TRIANGLES, 
                           geom.const_triangles().size()*3,
                           GL_UNSIGNED_INT, 0);
          }
        else
#endif
          glDrawElements(GL_TRIANGLES, 
                         geom.const_triangles().size()*3,
                         GL_UNSIGNED_INT, &(geom.const_triangles()[0]));
        
#if 0
        if (_drawGeometryNormals) {
          // draw normals...
          glLineWidth(2.0*cvcstate(s)("line_width").value<float>());
          glBegin(GL_LINES);
          for (int i=0; i< geom.points().size(); i++) {
            double len = sqrt(geom.const_normals()[i][0]*geom.const_normals()[i][0]+geom.const_normals()[i][1]*geom.const_normals()[i][1]+geom.const_normals()[i][2]*geom.const_normals()[i][2]);
            len *= 2.0;
            glVertex3d(geom.const_points()[i][0],geom.const_points()[i][1],geom.const_points()[i][2]);
            glVertex3d(geom.const_points()[i][0]+geom.const_normals()[i][0]/len,
                       geom.const_points()[i][1]+geom.const_normals()[i][1]/len,
                       geom.const_points()[i][2]+geom.const_normals()[i][2]/len);	     	      
          }
          glEnd();
          
        }
#endif
      }
    else if(cvcstate(s)("render_mode").value()=="triangles_flat")
      {
        // arand, 9-12-2011: added "flat" mode. This is probably much slower
        //                   but displays un-oriented triangulations much nicer

        glBegin(GL_TRIANGLES);
        for (int i=0; i<geom.triangles().size(); i++) {
          int t1,t2,t0;
          t0 = geom.const_triangles()[i][0];
          t1 = geom.const_triangles()[i][1];
          t2 = geom.const_triangles()[i][2];

          double nx,ny,nz;		
          double v1x,v2x,v1y,v2y,v1z,v2z;
		
          v1x = geom.const_points()[t1][0] - geom.const_points()[t0][0];
          v1y = geom.const_points()[t1][1] - geom.const_points()[t0][1];
          v1z = geom.const_points()[t1][2] - geom.const_points()[t0][2];
          v2x = geom.const_points()[t2][0] - geom.const_points()[t0][0];
          v2y = geom.const_points()[t2][1] - geom.const_points()[t0][1];
          v2z = geom.const_points()[t2][2] - geom.const_points()[t0][2];

          double lv1 = sqrt(v1x*v1x+v1y*v1y+v1z*v1z);
          double lv2 = sqrt(v2x*v2x+v2y*v2y+v2z*v2z);

          nx = (v1y*v2z-v1z*v2y)/(lv1*lv2);
          ny = (v1z*v2x-v1x*v2z)/(lv1*lv2);
          nz = (v1x*v2y-v1y*v2x)/(lv1*lv2);

          glNormal3d(nx,ny,nz);
          if (!geom.const_colors().empty())
            glColor3d(geom.const_colors()[t0][0],geom.const_colors()[t0][1],geom.const_colors()[t0][2]);
          glVertex3d(geom.const_points()[t0][0],geom.const_points()[t0][1],geom.const_points()[t0][2]);
          glVertex3d(geom.const_points()[t1][0],geom.const_points()[t1][1],geom.const_points()[t1][2]);
          glVertex3d(geom.const_points()[t2][0],geom.const_points()[t2][1],geom.const_points()[t2][2]);
        }
        glEnd();	      

#if 0
        if (_drawGeometryNormals) {
          // draw normals...
          glLineWidth(2.0*cvcstate(s)("line_width").value<float>());
          glBegin(GL_LINES);
          for (int i=0; i< geom.points().size(); i++) {
            double len = sqrt(geom.const_normals()[i][0]*geom.const_normals()[i][0]+geom.const_normals()[i][1]*geom.const_normals()[i][1]+geom.const_normals()[i][2]*geom.const_normals()[i][2]);
            len *= 2.0;
            glVertex3d(geom.const_points()[i][0],geom.const_points()[i][1],geom.const_points()[i][2]);
            glVertex3d(geom.const_points()[i][0]+geom.const_normals()[i][0]/len,
                       geom.const_points()[i][1]+geom.const_normals()[i][1]/len,
                       geom.const_points()[i][2]+geom.const_normals()[i][2]/len);	     	      
          }
          glEnd();
        }
#endif
      }
    else if(cvcstate(s)("render_mode").value()=="triangle_wireframe")
      {
        glLineWidth(cvcstate(s)("line_width").value<float>());
        
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
#if 0
        if(usingVBO())
          {
            glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER_ARB,
                            scene_geom.vboTriElementArrayBufferID);
            glDrawElements(GL_TRIANGLES, 
                           geom.const_triangles().size()*3,
                           GL_UNSIGNED_INT, 0);
          }
        else
#endif
          glDrawElements(GL_TRIANGLES, 
                         geom.const_triangles().size()*3,
                         GL_UNSIGNED_INT, &(geom.const_triangles()[0]));
      }
    else if(cvcstate(s)("render_mode").value()=="triangle_filled_wire")
      {
        glLineWidth(cvcstate(s)("line_width").value<float>());

        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
        glEnable(GL_POLYGON_OFFSET_FILL);
        glPolygonOffset(1.0,1.0);
#if 0
        if(usingVBO())
          {
            glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER_ARB,
                            scene_geom.vboTriElementArrayBufferID);
            glDrawElements(GL_TRIANGLES, 
                           geom.const_triangles().size()*3,
                           GL_UNSIGNED_INT, 0);
          }
        else
#endif
          glDrawElements(GL_TRIANGLES, 
                         geom.const_triangles().size()*3,
                         GL_UNSIGNED_INT, &(geom.const_triangles()[0]));

#if 0
        if (_drawGeometryNormals) {
          // draw normals...
          glLineWidth(2.0*cvcstate(s)("line_width").value<float>());
          glBegin(GL_LINES);
          for (int i=0; i< geom.points().size(); i++) {
            double len = sqrt(geom.const_normals()[i][0]*geom.const_normals()[i][0]+geom.const_normals()[i][1]*geom.const_normals()[i][1]+geom.const_normals()[i][2]*geom.const_normals()[i][2]);
            len *= 2.0;
            glVertex3d(geom.const_points()[i][0],geom.const_points()[i][1],geom.const_points()[i][2]);
            glVertex3d(geom.const_points()[i][0]+geom.const_normals()[i][0]/len,
                       geom.const_points()[i][1]+geom.const_normals()[i][1]/len,
                       geom.const_points()[i][2]+geom.const_normals()[i][2]/len);	     	      
          }
          glEnd();

        }
#endif

        glPolygonOffset(0.0,0.0);
        glDisable(GL_POLYGON_OFFSET_FILL);
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
        glDisable(GL_LIGHTING);
        glDisableClientState(GL_COLOR_ARRAY);
        glColor3f(0.0,0.0,0.0); //black wireframe
#if 0
        if(usingVBO())
          {
            glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER_ARB,
                            scene_geom.vboTriElementArrayBufferID);
            glDrawElements(GL_TRIANGLES, 
                           geom.const_triangles().size()*3,
                           GL_UNSIGNED_INT, 0);
          }
        else
#endif
          glDrawElements(GL_TRIANGLES, 
                         geom.const_triangles().size()*3,
                         GL_UNSIGNED_INT, &(geom.const_triangles()[0]));
      }
    else if(cvcstate(s)("render_mode").value()=="triangle_flat_filled_wire")
      {
        glLineWidth(cvcstate(s)("line_width").value<float>());

        // new hack mode...
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
        glEnable(GL_POLYGON_OFFSET_FILL);
        glPolygonOffset(1.0,1.0);
	    
        glBegin(GL_TRIANGLES);
        for (int i=0; i<geom.triangles().size(); i++) {
          int t1,t2,t0;
          t0 = geom.const_triangles()[i][0];
          t1 = geom.const_triangles()[i][1];
          t2 = geom.const_triangles()[i][2];
	      
          double nx,ny,nz;
          double v1x,v2x,v1y,v2y,v1z,v2z;
	      
          v1x = geom.const_points()[t1][0] - geom.const_points()[t0][0];
          v1y = geom.const_points()[t1][1] - geom.const_points()[t0][1];
          v1z = geom.const_points()[t1][2] - geom.const_points()[t0][2];
          v2x = geom.const_points()[t2][0] - geom.const_points()[t0][0];
          v2y = geom.const_points()[t2][1] - geom.const_points()[t0][1];
          v2z = geom.const_points()[t2][2] - geom.const_points()[t0][2];
		
          double lv1 = sqrt(v1x*v1x+v1y*v1y+v1z*v1z);
          double lv2 = sqrt(v2x*v2x+v2y*v2y+v2z*v2z);
	      
          nx = (v1y*v2z-v1z*v2y)/(lv1*lv2);
          ny = (v1z*v2x-v1x*v2z)/(lv1*lv2);
          nz = (v1x*v2y-v1y*v2x)/(lv1*lv2);
	      
          glNormal3d(nx,ny,nz);
          if (!geom.const_colors().empty())
            glColor3d(geom.const_colors()[t0][0],geom.const_colors()[t0][1],geom.const_colors()[t0][2]);
          glVertex3d(geom.const_points()[t0][0],geom.const_points()[t0][1],geom.const_points()[t0][2]);
          glVertex3d(geom.const_points()[t1][0],geom.const_points()[t1][1],geom.const_points()[t1][2]);
          glVertex3d(geom.const_points()[t2][0],geom.const_points()[t2][1],geom.const_points()[t2][2]);
        }
        glEnd();

#if 0
        if (_drawGeometryNormals) {
          // draw normals...
          glLineWidth(2.0*cvcstate(s)("line_width").value<float>());
          glBegin(GL_LINES);
          for (int i=0; i< geom.points().size(); i++) {
            double len = sqrt(geom.const_normals()[i][0]*geom.const_normals()[i][0]+geom.const_normals()[i][1]*geom.const_normals()[i][1]+geom.const_normals()[i][2]*geom.const_normals()[i][2]);
            len *= 2.0;
            glVertex3d(geom.const_points()[i][0],geom.const_points()[i][1],geom.const_points()[i][2]);
            glVertex3d(geom.const_points()[i][0]+geom.const_normals()[i][0]/len,
                       geom.const_points()[i][1]+geom.const_normals()[i][1]/len,
                       geom.const_points()[i][2]+geom.const_normals()[i][2]/len);	     	      
          }
          glEnd();
        }
#endif
	      
        glPolygonOffset(0.0,0.0);
        glDisable(GL_POLYGON_OFFSET_FILL);
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
        glDisable(GL_LIGHTING);
        glDisableClientState(GL_COLOR_ARRAY);
        glColor3f(0.0,0.0,0.0); //black wireframe

        // draw wireframe...
        glBegin(GL_TRIANGLES);
        for (int i=0; i<geom.triangles().size(); i++) {
          int t1,t2,t0;
          t0 = geom.const_triangles()[i][0];
          t1 = geom.const_triangles()[i][1];
          t2 = geom.const_triangles()[i][2];
	      
          double nx,ny,nz;
          double v1x,v2x,v1y,v2y,v1z,v2z;
	      
          v1x = geom.const_points()[t1][0] - geom.const_points()[t0][0];
          v1y = geom.const_points()[t1][1] - geom.const_points()[t0][1];
          v1z = geom.const_points()[t1][2] - geom.const_points()[t0][2];
          v2x = geom.const_points()[t2][0] - geom.const_points()[t0][0];
          v2y = geom.const_points()[t2][1] - geom.const_points()[t0][1];
          v2z = geom.const_points()[t2][2] - geom.const_points()[t0][2];
		
          double lv1 = sqrt(v1x*v1x+v1y*v1y+v1z*v1z);
          double lv2 = sqrt(v2x*v2x+v2y*v2y+v2z*v2z);
	      
          nx = (v1y*v2z-v1z*v2y)/(lv1*lv2);
          ny = (v1z*v2x-v1x*v2z)/(lv1*lv2);
          nz = (v1x*v2y-v1y*v2x)/(lv1*lv2);
	      
          glNormal3d(nx,ny,nz);
          glVertex3d(geom.const_points()[t0][0],geom.const_points()[t0][1],geom.const_points()[t0][2]);
          glVertex3d(geom.const_points()[t1][0],geom.const_points()[t1][1],geom.const_points()[t1][2]);
          glVertex3d(geom.const_points()[t2][0],geom.const_points()[t2][1],geom.const_points()[t2][2]);
        }
        glEnd();
      }
    else if(cvcstate(s)("render_mode").value()=="quads")
      {
#if 0
        if(usingVBO())
          {
            glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER_ARB,
                            scene_geom.vboQuadElementArrayBufferID);
            glDrawElements(GL_QUADS, 
                           geom.const_quads().size()*4,
                           GL_UNSIGNED_INT, 0);
          }
        else
#endif
          glDrawElements(GL_QUADS, 
                         geom.const_quads().size()*4,
                         GL_UNSIGNED_INT, &(geom.const_quads()[0]));
      }
    else if(cvcstate(s)("render_mode").value()=="quad_wireframe")
      {
        glLineWidth(cvcstate(s)("line_width").value<float>());

        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
#if 0
        if(usingVBO())
          {
            glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER_ARB,
                            scene_geom.vboQuadElementArrayBufferID);
            glDrawElements(GL_QUADS, 
                           geom.const_quads().size()*4,
                           GL_UNSIGNED_INT, 0);
          }
        else
#endif
          glDrawElements(GL_QUADS, 
                         geom.const_quads().size()*4,
                         GL_UNSIGNED_INT, &(geom.const_quads()[0]));
      }
    else if(cvcstate(s)("render_mode").value()=="quad_filled_wire")
      {
        glLineWidth(cvcstate(s)("line_width").value<float>());

        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
        glEnable(GL_POLYGON_OFFSET_FILL);
        glPolygonOffset(1.0,1.0);
#if 0
        if(usingVBO())
          {
            glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER_ARB,
                            scene_geom.vboQuadElementArrayBufferID);
            glDrawElements(GL_QUADS, 
                           geom.const_quads().size()*4,
                           GL_UNSIGNED_INT, 0);
          }
        else
#endif
          glDrawElements(GL_QUADS, 
                         geom.const_quads().size()*4,
                         GL_UNSIGNED_INT, &(geom.const_quads()[0]));
        glPolygonOffset(0.0,0.0);
        glDisable(GL_POLYGON_OFFSET_FILL);
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
        glDisable(GL_LIGHTING);
        glDisableClientState(GL_COLOR_ARRAY);
        glColor3f(0.0,0.0,0.0); //black wireframe
#if 0
        if(usingVBO())
          {
            glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER_ARB,
                            scene_geom.vboQuadElementArrayBufferID);
            glDrawElements(GL_QUADS, 
                           geom.const_quads().size()*4,
                           GL_UNSIGNED_INT, 0);
          }
        else
#endif
          glDrawElements(GL_QUADS, 
                         geom.const_quads().size()*4,
                         GL_UNSIGNED_INT, &(geom.const_quads()[0]));
      }
    else if(cvcstate(s)("render_mode").value()=="tetra")
      {
        glLineWidth(cvcstate(s)("line_width").value<float>());

#if 0
        if(usingVBO())
          {
            glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER_ARB,
                            scene_geom.vboTriElementArrayBufferID);
            glDrawElements(GL_TRIANGLES, 
                           geom.const_triangles().size()*3,
                           GL_UNSIGNED_INT, 0);
          }
        else
#endif
          glDrawElements(GL_TRIANGLES,
                         geom.const_triangles().size()*3,
                         GL_UNSIGNED_INT, &(geom.const_triangles()[0]));

#if 0
        if(usingVBO())
          {
            glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER_ARB,
                            scene_geom.vboLineElementArrayBufferID);
            glDrawElements(GL_LINES, 
                           geom.const_lines().size()*2,
                           GL_UNSIGNED_INT, 0);
          }
        else
#endif
          glDrawElements(GL_LINES, 
                         geom.const_lines().size()*2,
                         GL_UNSIGNED_INT, &(geom.const_lines()[0]));
      }
    else if(cvcstate(s)("render_mode").value()=="hexa")
      {
        glLineWidth(cvcstate(s)("line_width").value<float>());

#if 0
        if(usingVBO())
          {
            glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER_ARB,
                            scene_geom.vboQuadElementArrayBufferID);
            glDrawElements(GL_QUADS, 
                           geom.const_quads().size()*4,
                           GL_UNSIGNED_INT, 0);
          }
        else
#endif
          glDrawElements(GL_QUADS, 
                         geom.const_quads().size()*4,
                         GL_UNSIGNED_INT, &(geom.const_quads()[0]));
        
#if 0
        if(usingVBO())
          {
            glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER_ARB,
                            scene_geom.vboLineElementArrayBufferID);
            glDrawElements(GL_LINES, 
                           geom.const_lines().size()*2,
                           GL_UNSIGNED_INT, 0);
          }
        else
#endif
          glDrawElements(GL_LINES, 
                         geom.const_lines().size()*2,
                         GL_UNSIGNED_INT, &(geom.const_lines()[0]));
      }

    glDisableClientState(GL_VERTEX_ARRAY);
    glDisableClientState(GL_NORMAL_ARRAY);
    glDisableClientState(GL_TEXTURE_COORD_ARRAY);
    glDisableClientState(GL_COLOR_ARRAY);

#if 0    
    if(usingVBO())
      {
        glBindBufferARB(GL_ARRAY_BUFFER_ARB,0);
        glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER_ARB,0);
      }
#endif
    
    //restore previous setting for polygon mode
    glPolygonMode(GL_FRONT,params[0]);
    glPolygonMode(GL_BACK,params[1]);

    if(hasName) glPopName();
  }

  void GeometryRenderer::handleStateChanged(const std::string& childState)
  {
    cvcapp.log(2,str(boost::format("%s :: state changed: %s\n")
                     % BOOST_CURRENT_FUNCTION
                     % childState));
  }

  void GeometryRenderer::doDrawBoundingBox(const std::string& s)
  {
    BoundingBox bbox;
    bbox = cvcstate(s).data<BoundingBox>();
    doDrawBoundingBox(bbox);
  }

  void GeometryRenderer::doDrawBoundingBox(const BoundingBox& bbox)
  {
    double minx, miny, minz;
    double maxx, maxy, maxz;
    float bgcolor[4];

    minx = bbox.XMin();
    miny = bbox.YMin();
    minz = bbox.ZMin();
    maxx = bbox.XMax();
    maxy = bbox.YMax();
    maxz = bbox.ZMax();

    glPushAttrib(GL_CURRENT_BIT | GL_ENABLE_BIT);
    glLineWidth(1.0);
    glEnable(GL_LINE_SMOOTH);
    glDisable(GL_LIGHTING);
    glGetFloatv(GL_COLOR_CLEAR_VALUE, bgcolor);
    glColor4f(1.0-bgcolor[0],1.0-bgcolor[1],1.0-bgcolor[2],1.0-bgcolor[3]);
    /* front face */
    glBegin(GL_LINE_LOOP);
    glVertex3d(minx,miny,minz);
    glVertex3d(maxx,miny,minz);
    glVertex3d(maxx,maxy,minz);
    glVertex3d(minx,maxy,minz);
    glEnd();
    /* back face */
    glBegin(GL_LINE_LOOP);
    glVertex3d(minx,miny,maxz);
    glVertex3d(maxx,miny,maxz);
    glVertex3d(maxx,maxy,maxz);
    glVertex3d(minx,maxy,maxz);
    glEnd();
    /* connecting lines */
    glBegin(GL_LINES);
    glVertex3d(minx,maxy,minz);
    glVertex3d(minx,maxy,maxz);
    glVertex3d(minx,miny,minz);
    glVertex3d(minx,miny,maxz);
    glVertex3d(maxx,maxy,minz);
    glVertex3d(maxx,maxy,maxz);
    glVertex3d(maxx,miny,minz);
    glVertex3d(maxx,miny,maxz);
    glEnd();
    glDisable(GL_LINE_SMOOTH);
    glPopAttrib();
  }
}
