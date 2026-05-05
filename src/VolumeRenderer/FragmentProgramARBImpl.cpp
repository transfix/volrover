/*
  Copyright 2002-2003 The University of Texas at Austin

	Authors: Anthony Thane <thanea@ices.utexas.edu>
        Advisor: Chandrajit Bajaj <bajaj@cs.utexas.edu>

  This file is part of VolumeLibrary.

  VolumeLibrary is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License version 2.1 as published by the Free Software Foundation.

  VolumeLibrary is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/

// FragmentProgramARBImpl.cpp: implementation of the FragmentProgramARBImpl class.
//
//////////////////////////////////////////////////////////////////////

#include <VolumeRenderer/FragmentProgramARBImpl.h>
#include <stdio.h>
#include <string.h>

#include <boost/current_function.hpp>
#include <boost/format.hpp>
#include <iostream>

#ifdef UPLOAD_DATA_RESIZE_HACK
#include <VolMagick/VolMagick.h>
#endif

#ifdef GLERROR_OUTPUT
#define PRINT_GLERROR                              \
        {                                          \
          using namespace std;                     \
          using namespace boost; using boost::format;                   \
          cerr<<str(format("%1%,%2%,%3%,glGetError()==%4%\n")        \
                    % BOOST_CURRENT_FUNCTION       \
                    % __FILE__                     \
                    % __LINE__                     \
                    % glGetError()                 \
                    );                             \
        }
#else
#define PRINT_GLERROR
#endif


using namespace OpenGLVolumeRendering;

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

FragmentProgramARBImpl::FragmentProgramARBImpl()
{
	m_Initialized = false;
	m_Width = -1;
	m_Height = -1;
	m_Depth = -1;
	m_GL_VERSION_1_2 = true;
}

FragmentProgramARBImpl::~FragmentProgramARBImpl()
{

}

// Initializes the renderer.  Should be called again if the renderer is
// moved to a different openGL context.  If this returns false, do not try
// to use it to do volumeRendering
bool FragmentProgramARBImpl::initRenderer()
{
	if (!UnshadedBase::initRenderer() || !initExtensions() || !initTextureNames() || !initFragmentProgram()) {
		m_Initialized = false;
		m_Width = -1;
		m_Height = -1;
		m_Depth = -1;
		return false;
	}
	else {
		m_Initialized = true;
		return true;
	}
}

// Makes the check necessary to determine if this renderer is 
// compatible with the hardware its running on
bool FragmentProgramARBImpl::checkCompatibility() const
{
  //so here's the deal, I've tried this on both mac osx and win32 with the same graphics hardware (Radeon9800).
  //This code requires 3D textures so we either need GL_VERSION_1_2 or GL_EXT_texture3D...
  //On win32 GL_EXT_texture3D is supported, but GL_VERSION_1_2 is not supported.  And the opposite is true on mac osx.
  //It's the same hardware, what gives?? :(
  return (glewIsSupported("GL_VERSION_1_2") ||
	  glewIsSupported("GL_EXT_texture3D")) &&
    glewIsSupported("GL_ARB_vertex_program") &&
    glewIsSupported("GL_ARB_fragment_program") &&
    glewIsSupported("GL_ARB_multitexture");
}

// Uploads colormapped data
bool FragmentProgramARBImpl::uploadColormappedData(const GLubyte* data, int width, int height, int depth)
{
  const GLubyte * upload_data = data;

	// bail if we haven't been initialized properly
	if (!m_Initialized) return false;

	// clear previous errors
	GLenum error = glGetError();
        PRINT_GLERROR;

#ifdef UPLOAD_DATA_RESIZE_HACK
	VolMagick::Volume vol(data,
			      VolMagick::Dimension(width,height,depth),
			      VolMagick::UChar);
	vol.resize(VolMagick::Dimension(128,128,128)); //force 128^3 to hopefully skirt an issue on intel macs!
	upload_data = *vol;
	width = vol.XDim();
	height = vol.YDim();
	depth = vol.ZDim();
#endif

	printf("FragmentProgramARBImpl::uploadColormappedData: m_DataTextureName: %d, width: %d, height: %d, depth: %d\n",m_DataTextureName,width,height,depth);

	if(m_GL_VERSION_1_2)
	  glBindTexture(GL_TEXTURE_3D, m_DataTextureName);
	else
	  glBindTexture(GL_TEXTURE_3D_EXT, m_DataTextureName);

	if (width!=m_Width || height!=m_Height || depth!=m_Depth) {
	  if(m_GL_VERSION_1_2)
           {
             //first check for errors with our parameters
             glTexImage3D(GL_PROXY_TEXTURE_3D, 0, GL_LUMINANCE, width, height,
                             depth, 0, GL_LUMINANCE, GL_UNSIGNED_BYTE, upload_data);

	     int i;
             GLint params[4];
             for (i=0; i<4; i++)
               params[i] = -1;
             glGetTexLevelParameteriv(GL_PROXY_TEXTURE_3D, 0, GL_TEXTURE_WIDTH, params);
             glGetTexLevelParameteriv(GL_PROXY_TEXTURE_3D, 0, GL_TEXTURE_HEIGHT, params + 1);
             glGetTexLevelParameteriv(GL_PROXY_TEXTURE_3D, 0, GL_TEXTURE_DEPTH, params + 2);
             glGetTexLevelParameteriv(GL_PROXY_TEXTURE_3D, 0, GL_TEXTURE_INTERNAL_FORMAT, params + 3);
             for(i=0; i<4; i++)
               printf("%d\n",params[i]);
 
	    glTexImage3D(GL_TEXTURE_3D, 0, GL_LUMINANCE, width, height,
			    depth, 0, GL_LUMINANCE, GL_UNSIGNED_BYTE, upload_data);
           }
	  else
           {
	    glTexImage3DEXT(GL_TEXTURE_3D_EXT, 0, GL_LUMINANCE, width, height,
			       depth, 0, GL_LUMINANCE, GL_UNSIGNED_BYTE, upload_data);
           }
	}
	else {
	  if(m_GL_VERSION_1_2)
	    glTexSubImage3D(GL_TEXTURE_3D, 0, 0,0,0, width, height,
			       depth, GL_LUMINANCE, GL_UNSIGNED_BYTE, upload_data);
	  else
	    glTexSubImage3DEXT(GL_TEXTURE_3D_EXT, 0, 0,0,0, width, height,
				  depth, GL_LUMINANCE, GL_UNSIGNED_BYTE, upload_data);
		
	}

	if(m_GL_VERSION_1_2)
	  {
	    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);
	    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	  }
	else
	  {
	    glTexParameteri(GL_TEXTURE_3D_EXT, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	    glTexParameteri(GL_TEXTURE_3D_EXT, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	    glTexParameteri(GL_TEXTURE_3D_EXT, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);
	    glTexParameteri(GL_TEXTURE_3D_EXT, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	    glTexParameteri(GL_TEXTURE_3D_EXT, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	  }

	// save the width height and depth
	m_Width = width;   m_HintDimX = width;
	m_Height = height; m_HintDimY = height;
	m_Depth = depth;   m_HintDimZ = depth;

	// test for error
	error = glGetError();
        PRINT_GLERROR;
	if (error == GL_NO_ERROR) {
		return true;
	}
	else {
		return false;
	}

}

// Tests to see if the given parameters would return an error
bool FragmentProgramARBImpl::testColormappedData(int width, int height, int depth)
{
	// bail if we haven't been initialized properly
	if (!m_Initialized) return false;

	// nothing above 512
	if (width>512 || height>512 || depth>512) {
		return false;
	}

	// clear previous errors
	GLenum error;
	int c =0;
	while (glGetError()!=GL_NO_ERROR && c<10) c++;
        PRINT_GLERROR;
	if(m_GL_VERSION_1_2)
	  glTexImage3D(GL_PROXY_TEXTURE_3D, 0, GL_LUMINANCE, width, height,
			  depth, 0, GL_LUMINANCE, GL_UNSIGNED_BYTE, 0);
	else
	  glTexImage3DEXT(GL_PROXY_TEXTURE_3D_EXT, 0, GL_LUMINANCE, width, height,
			     depth, 0, GL_LUMINANCE, GL_UNSIGNED_BYTE, 0);
	
	// test for error
	error = glGetError();
        PRINT_GLERROR;
	if (error == GL_NO_ERROR) {
		return true;
	}
	else {
		return false;
	}
}

// Uploads the transfer function for the colormapped data
bool FragmentProgramARBImpl::uploadColorMap(const GLubyte* colorMap)
{
	// bail if we haven't been initialized properly
	if (!m_Initialized) return false;

	// clear previous errors
	GLenum error = glGetError();
        PRINT_GLERROR;	

	//printf("FragmentProgramARBImpl::uploadColorMap: m_TransferTextureName: %d\n",m_TransferTextureName);

	glBindTexture(GL_TEXTURE_1D, m_TransferTextureName);

	glTexImage1D(GL_TEXTURE_1D, 0, 4, 256, 0, GL_RGBA, GL_UNSIGNED_BYTE, colorMap);

	glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);

	// test for error
	error = glGetError();
        PRINT_GLERROR;
	if (error == GL_NO_ERROR) {
		return true;
	}
	else {
		return false;
	}
}

// Uploads the transfer function for the colormapped data
bool FragmentProgramARBImpl::uploadColorMap(const GLfloat* colorMap)
{
	// bail if we haven't been initialized properly
	if (!m_Initialized) return false;

	// clear previous errors
	GLenum error = glGetError();
        PRINT_GLERROR;	

	//printf("FragmentProgramARBImpl::uploadColorMap: m_TransferTextureName: %d\n",m_TransferTextureName);

	glBindTexture(GL_TEXTURE_1D, m_TransferTextureName);

	glTexImage1D(GL_TEXTURE_1D, 0, 4, 256, 0, GL_RGBA, GL_FLOAT, colorMap);

	glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);

	// test for error
	error = glGetError();
        PRINT_GLERROR;
	if (error == GL_NO_ERROR) {
		return true;
	}
	else {
		return false;
	}
}

// Performs the actual rendering.
bool FragmentProgramARBImpl::renderVolume()
{
	// bail if we haven't been initialized properly
	if (!m_Initialized) return false;

        PRINT_GLERROR;

	// set up the state
	glPushAttrib( GL_ENABLE_BIT | GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
	glColor4f(1.0, 1.0, 1.0, 1.0);
	glDisable(GL_CULL_FACE);
	glDisable(GL_LIGHTING);
	glEnable(GL_BLEND);

        initFragmentProgram();

	// get the fragment program ready
	glEnable(GL_FRAGMENT_PROGRAM_ARB);
	glBindProgramARB(GL_FRAGMENT_PROGRAM_ARB, m_FragmentProgramName);

	//glBlendFunc( GL_ONE, GL_ONE_MINUS_SRC_ALPHA );
	//glBlendFunc( GL_SRC_ALPHA, GL_ONE );
	glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );
	glDepthMask( GL_FALSE );

	// bind the transfer function
	glActiveTextureARB(GL_TEXTURE1_ARB);
	glEnable(GL_TEXTURE_1D);
	glBindTexture(GL_TEXTURE_1D, m_TransferTextureName);

	// bind the data texture
	glActiveTextureARB(GL_TEXTURE0_ARB);
	if(m_GL_VERSION_1_2)
	  {
	    glEnable(GL_TEXTURE_3D);
	    glBindTexture(GL_TEXTURE_3D, m_DataTextureName);
	  }
	else
	  {
	    glEnable(GL_TEXTURE_3D_EXT);
	    glBindTexture(GL_TEXTURE_3D_EXT, m_DataTextureName);
	  }

	computePolygons();

	convertToTriangles();

	renderTriangles();

	// unbind the fragment program
	glBindProgramARB(GL_FRAGMENT_PROGRAM_ARB, 0);

	// restore the state
	glPopAttrib();

        PRINT_GLERROR;

	return true;

}

// Initializes the necessary extensions.
bool FragmentProgramARBImpl::initExtensions()
{
  return ((m_GL_VERSION_1_2=glewIsSupported("GL_VERSION_1_2")) ||
	  glewIsSupported("GL_EXT_texture3D")) &&
    glewIsSupported("GL_ARB_vertex_program") &&
    glewIsSupported("GL_ARB_fragment_program") &&
    glewIsSupported("GL_ARB_multitexture");
}

// Gets the opengl texture IDs
bool FragmentProgramARBImpl::initTextureNames()
{
	// clear previous errors
	GLenum error = glGetError();
        PRINT_GLERROR;

	// get the names
	glGenTextures(1, &m_DataTextureName);
	glGenTextures(1, &m_TransferTextureName);

	// test for error
	error = glGetError();
        PRINT_GLERROR;
	if (error==GL_NO_ERROR) {
		return true;
	}
	else {
		return false;
	}
}

// Gets the fragment program ready
bool FragmentProgramARBImpl::initFragmentProgram()
{
	// clear previous errors
	GLenum error = glGetError();
        PRINT_GLERROR;

	glGenProgramsARB(1, &m_FragmentProgramName);

	const GLubyte program[] = 
          "!!ARBfp1.0\n"
          "PARAM c0 = {0.5, 1, 2.7182817, 0};\n"
          "TEMP R0;\n"
          "TEX R0.x, fragment.texcoord[0].xyzx, texture[0], 3D;\n"
          "TEX result.color, R0.x, texture[1], 1D;\n"
          "END\n";
	int len = strlen((const char*)program);
	glBindProgramARB(GL_FRAGMENT_PROGRAM_ARB, m_FragmentProgramName);
	glProgramStringARB(GL_FRAGMENT_PROGRAM_ARB, GL_PROGRAM_FORMAT_ASCII_ARB, len, program);

	// test for error
	error = glGetError();
        PRINT_GLERROR;
	if (error==GL_NO_ERROR) {
		return true;
	}
	else {
		return false;
	}
}

// Render the actual triangles
void FragmentProgramARBImpl::renderTriangles()
{
  PRINT_GLERROR;

	// set up the client render state
	glEnableClientState(GL_VERTEX_ARRAY);
	glEnableClientState(GL_TEXTURE_COORD_ARRAY);

	glTexCoordPointer(3, GL_FLOAT, 0, m_TextureArray.get());
	glVertexPointer(3, GL_FLOAT, 0, m_VertexArray.get());

	// render the triangles
	glDrawElements(GL_TRIANGLES, m_NumTriangles*3, GL_UNSIGNED_INT, m_TriangleArray.get());
	
	glDisableClientState(GL_VERTEX_ARRAY);
	glDisableClientState(GL_TEXTURE_COORD_ARRAY);


  PRINT_GLERROR;
}
