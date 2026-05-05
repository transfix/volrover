/*
  Copyright 2002-2003,2011 The University of Texas at Austin

	Authors: Anthony Thane <thanea@ices.utexas.edu>
                 Joe Rivera <transfix@ices.utexas.edu>
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

// FragmentProgramGLSLImpl.cpp: implementation of the FragmentProgramGLSLImpl class.
//
//////////////////////////////////////////////////////////////////////

#include <VolumeRenderer/FragmentProgramGLSLImpl.h>
#include <cstdio>
#include <cstring>

#include <string>

#include <boost/current_function.hpp>
#include <boost/format.hpp>
#include <iostream>

#ifdef UPLOAD_DATA_RESIZE_HACK
#include <VolMagick/VolMagick.h>
#endif

#include <cvc_compat.h>

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


using namespace OpenGLVolumeRendering;

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

FragmentProgramGLSLImpl::FragmentProgramGLSLImpl() :
  _initialized(false),
  _isGL2_0(false),
  _glslAvailable(false),
  _program(0),
  _vertexShader(0),
  _fragmentShader(0),
  _lightsTex(0),
  _lightColorsTex(0)
{
  _GLSL_UNIFORM_TYPE_MAP[GL_FLOAT] = "GL_FLOAT";
  _GLSL_UNIFORM_TYPE_MAP[GL_FLOAT_VEC2] = "GL_FLOAT_VEC2";
  _GLSL_UNIFORM_TYPE_MAP[GL_FLOAT_VEC3] = "GL_FLOAT_VEC3";
  _GLSL_UNIFORM_TYPE_MAP[GL_FLOAT_VEC4] = "GL_FLOAT_VEC4";
  _GLSL_UNIFORM_TYPE_MAP[GL_INT] = "GL_INT";
  _GLSL_UNIFORM_TYPE_MAP[GL_INT_VEC2] = "GL_INT_VEC2";
  _GLSL_UNIFORM_TYPE_MAP[GL_INT_VEC3] = "GL_INT_VEC3";
  _GLSL_UNIFORM_TYPE_MAP[GL_INT_VEC4] = "GL_INT_VEC4";
  _GLSL_UNIFORM_TYPE_MAP[GL_BOOL] = "GL_BOOL";
  _GLSL_UNIFORM_TYPE_MAP[GL_BOOL_VEC2] = "GL_BOOL_VEC2";
  _GLSL_UNIFORM_TYPE_MAP[GL_BOOL_VEC3] = "GL_BOOL_VEC3";
  _GLSL_UNIFORM_TYPE_MAP[GL_BOOL_VEC4] = "GL_BOOL_VEC4";
  _GLSL_UNIFORM_TYPE_MAP[GL_FLOAT_MAT2] = "GL_FLOAT_MAT2";
  _GLSL_UNIFORM_TYPE_MAP[GL_FLOAT_MAT3] = "GL_FLOAT_MAT3";
  _GLSL_UNIFORM_TYPE_MAP[GL_FLOAT_MAT4] = "GL_FLOAT_MAT4";
#if 0
  _GLSL_UNIFORM_TYPE_MAP[GL_FLOAT_MAT2x3] = "GL_FLOAT_MAT2x3";
  _GLSL_UNIFORM_TYPE_MAP[GL_FLOAT_MAT2x4] = "GL_FLOAT_MAT2x4";
  _GLSL_UNIFORM_TYPE_MAP[GL_FLOAT_MAT3x2] = "GL_FLOAT_MAT3x2";
  _GLSL_UNIFORM_TYPE_MAP[GL_FLOAT_MAT3x4] = "GL_FLOAT_MAT3x4";
  _GLSL_UNIFORM_TYPE_MAP[GL_FLOAT_MAT4x2] = "GL_FLOAT_MAT4x2";
  _GLSL_UNIFORM_TYPE_MAP[GL_FLOAT_MAT4x3] = "GL_FLOAT_MAT4x3";
#endif
  _GLSL_UNIFORM_TYPE_MAP[GL_SAMPLER_1D] = "GL_SAMPLER_1D";
  _GLSL_UNIFORM_TYPE_MAP[GL_SAMPLER_2D] = "GL_SAMPLER_2D";
  _GLSL_UNIFORM_TYPE_MAP[GL_SAMPLER_3D] = "GL_SAMPLER_3D";
  _GLSL_UNIFORM_TYPE_MAP[GL_SAMPLER_CUBE] = "GL_SAMPLER_CUBE";
  _GLSL_UNIFORM_TYPE_MAP[GL_SAMPLER_1D_SHADOW] = "GL_SAMPLER_1D_SHADOW";
  _GLSL_UNIFORM_TYPE_MAP[GL_SAMPLER_2D_SHADOW] = "GL_SAMPLER_2D_SHADOW";

  m_Width = -1;
  m_Height = -1;
  m_Depth = -1;
}

FragmentProgramGLSLImpl::~FragmentProgramGLSLImpl()
{

}

// Initializes the renderer.  Should be called again if the renderer is
// moved to a different openGL context.  If this returns false, do not try
// to use it to do volumeRendering
bool FragmentProgramGLSLImpl::initRenderer()
{
  if (!UnshadedBase::initRenderer() || !initExtensions() || !initTextureNames() || !initFragmentProgram()) {
    _initialized = false;
    m_Width = -1;
    m_Height = -1;
    m_Depth = -1;
    return false;
  }
  else {
    _initialized = true;
    return true;
  }
}

// Makes the check necessary to determine if this renderer is 
// compatible with the hardware its running on
bool FragmentProgramGLSLImpl::checkCompatibility() const
{
  return glewIsSupported("GL_VERSION_2_0") &&
    glewIsSupported("GL_ARB_vertex_shader") &&
    glewIsSupported("GL_ARB_fragment_shader") &&
    glewIsSupported("GL_ARB_multitexture");
}

// Uploads colormapped data
bool FragmentProgramGLSLImpl::uploadColormappedData(const GLubyte* data, int width, int height, int depth)
{
  const GLubyte * upload_data = data;

  // bail if we haven't been initialized properly
  if (!initialized()) return false;

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

  printf("FragmentProgramGLSLImpl::uploadColormappedData: m_DataTextureName: %d, width: %d, height: %d, depth: %d\n",m_DataTextureName,width,height,depth);

	
  glBindTexture(GL_TEXTURE_3D, m_DataTextureName);
	
  if (width!=m_Width || height!=m_Height || depth!=m_Depth) {
	
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
  else {
    glTexSubImage3D(GL_TEXTURE_3D, 0, 0,0,0, width, height,
                    depth, GL_LUMINANCE, GL_UNSIGNED_BYTE, upload_data);
  }

  glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
  glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);
  glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	
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
bool FragmentProgramGLSLImpl::testColormappedData(int width, int height, int depth)
{
  // bail if we haven't been initialized properly
  if (!initialized()) return false;

  // nothing above 512
  if (width>512 || height>512 || depth>512) {
    return false;
  }

  // clear previous errors
  GLenum error;
  int c =0;
  while (glGetError()!=GL_NO_ERROR && c<10) c++;
  PRINT_GLERROR;
  glTexImage3D(GL_PROXY_TEXTURE_3D, 0, GL_LUMINANCE, width, height,
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
bool FragmentProgramGLSLImpl::uploadColorMap(const GLubyte* colorMap)
{
  // bail if we haven't been initialized properly
  if (!initialized()) return false;

  // clear previous errors
  GLenum error = glGetError();
  PRINT_GLERROR;	

  //printf("FragmentProgramGLSLImpl::uploadColorMap: m_TransferTextureName: %d\n",m_TransferTextureName);

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

// Performs the actual rendering.
bool FragmentProgramGLSLImpl::renderVolume()
{
  // bail if we haven't been initialized properly
  if (!initialized()) return false;

  // set up the state
  glPushAttrib( GL_ENABLE_BIT | GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
  glColor4f(1.0, 1.0, 1.0, 1.0);
  glDisable(GL_CULL_FACE);
  glDisable(GL_LIGHTING);
  glEnable(GL_BLEND);

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
  glEnable(GL_TEXTURE_3D);
  glBindTexture(GL_TEXTURE_3D, m_DataTextureName);

  //use GLSL programs if they've been loaded, else default to fixed functionality
  bool using_prog = false;
  GLint num_shaders = 0;
  if(_isGL2_0)
  {
      glGetProgramiv(_program,GL_ATTACHED_SHADERS,&num_shaders);
      if(num_shaders > 0)
      {
          glUseProgram(_program);
          using_prog = true;
      }
  }
	
  computePolygons();

  convertToTriangles();

  renderTriangles();

  // unbind the GLSL program
  if(using_prog)
    glUseProgram(0);

  // restore the state
  glPopAttrib();

  return true;
}

// Initializes the necessary extensions.
bool FragmentProgramGLSLImpl::initExtensions()
{
  return (_isGL2_0 = glewIsSupported("GL_VERSION_2_0")) &&
    glewIsSupported("GL_ARB_vertex_shader") &&
    glewIsSupported("GL_ARB_fragment_shader") &&
    glewIsSupported("GL_ARB_multitexture");
}

// Gets the opengl texture IDs
bool FragmentProgramGLSLImpl::initTextureNames()
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
bool FragmentProgramGLSLImpl::initFragmentProgram()
{
  using namespace boost; using boost::format;
  using namespace std;

  _isGL2_0 = glewIsSupported("GL_VERSION_2_0");
  _glslAvailable =
    (GLEW_ARB_vertex_shader && GLEW_ARB_fragment_shader) ||
    _isGL2_0;
        
  cerr<<str(format("GLSL available: %s\n") % (_glslAvailable ? "yes" : "no"));
  cerr<<str(format("OpenGL 2.0: %s\n") % (_isGL2_0 ? "yes" : "no"));
        
  if(_isGL2_0)
    _program = glCreateProgram();
  else
    _program = glCreateProgramObjectARB();

  _initialized = true;

  //TODO: write up and load actual fragment shader! :)

  return _isGL2_0;
}

// Render the actual triangles
void FragmentProgramGLSLImpl::renderTriangles()
{
  //PRINT_GLERROR;

  // set up the client render state
  glEnableClientState(GL_VERTEX_ARRAY);
  glEnableClientState(GL_TEXTURE_COORD_ARRAY);

  glTexCoordPointer(3, GL_FLOAT, 0, m_TextureArray.get());
  glVertexPointer(3, GL_FLOAT, 0, m_VertexArray.get());

  // render the triangles
  glDrawElements(GL_TRIANGLES, m_NumTriangles*3, GL_UNSIGNED_INT, m_TriangleArray.get());
	
  glDisableClientState(GL_VERTEX_ARRAY);
  glDisableClientState(GL_TEXTURE_COORD_ARRAY);

  //  PRINT_GLERROR;
}


void FragmentProgramGLSLImpl::setVertexShader(const char *prog)
{
  using namespace std;

  if(!initialized())
    throw UninitializedException("Could not set vertex shader!");
    
  unsetVertexShader();
    
  if(_isGL2_0)
    {
      _vertexShader = glCreateShader(GL_VERTEX_SHADER);

      const char *progList[1] = { prog };
      glShaderSource(_vertexShader,1,progList,NULL);
      glCompileShader(_vertexShader);

      GLint compile_status = 0;
      glGetShaderiv(_vertexShader,GL_COMPILE_STATUS,&compile_status);
      if(compile_status == GL_FALSE) //if there was a compile error...
        {
          GLint log_length = 0;
          glGetShaderiv(_vertexShader,GL_INFO_LOG_LENGTH,&log_length);
          if(log_length > 0)
            {
              string log_string(log_length,0);
              glGetShaderInfoLog(_vertexShader,log_length,NULL,&(log_string[0]));
              throw ShaderCompileException("Error compiling vertex shader: " + log_string);
            }
          else
            throw ShaderCompileException("Error compiling vertex shader");
        }

      glAttachShader(_program, _vertexShader);
      glLinkProgram(_program);

      GLint link_status = 0;
      glGetProgramiv(_program,GL_LINK_STATUS,&link_status);
      if(link_status == GL_FALSE) //if there was a link error
        {
          GLint log_length = 0;
          glGetProgramiv(_program,GL_INFO_LOG_LENGTH,&log_length);
          if(log_length > 0)
            {
              string log_string(log_length,0);
              glGetProgramInfoLog(_program,log_length,NULL,&(log_string[0]));
              throw ShaderLinkException("Error linking GLSL program: " + log_string);
            }
          else
            throw ShaderLinkException("Error linking GLSL program");
        }
    }
  else
    {
      _vertexShader = glCreateShaderObjectARB(GL_VERTEX_SHADER_ARB);

      const char *progList[1] = { prog };
      glShaderSourceARB(_vertexShader,1,progList,NULL);
      glCompileShaderARB(_vertexShader);

      GLint compile_status = 0;
      glGetObjectParameterivARB(_vertexShader,GL_OBJECT_COMPILE_STATUS_ARB,&compile_status);
      if(compile_status == GL_FALSE)
        {
          GLint log_length = 0;
          glGetObjectParameterivARB(_vertexShader,GL_OBJECT_INFO_LOG_LENGTH_ARB,&log_length);
          if(log_length > 0)
            {
              string log_string(log_length,0);
              glGetInfoLogARB(_vertexShader,log_length,NULL,&(log_string[0]));
              throw ShaderCompileException("Error compiling vertex shader: " + log_string);
            }
          else
            throw ShaderCompileException("Error compiling vertex shader");
        }

      glAttachObjectARB(_program, _vertexShader);
      glLinkProgramARB(_program);

      GLint link_status = 0;
      glGetObjectParameterivARB(_program,GL_OBJECT_LINK_STATUS_ARB,&link_status);
      if(link_status == GL_FALSE)
        {
          GLint log_length = 0;
          glGetObjectParameterivARB(_program,GL_OBJECT_INFO_LOG_LENGTH_ARB,&log_length);
          if(log_length > 0)
            {
              string log_string(log_length,0);
              glGetInfoLogARB(_program,log_length,NULL,&(log_string[0]));
              throw ShaderCompileException("Error linking GLSL program: " + log_string);
            }
          else
            throw ShaderCompileException("Error linking GLSL program");
        }
    }

  updateUniforms();
}

void FragmentProgramGLSLImpl::setFragmentShader(const char *prog)
{
  using namespace std;

  if(!initialized())
    throw UninitializedException("Could not set fragment shader!");
    
  unsetFragmentShader();
    
  if(_isGL2_0)
    {
      _fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);

      const char *progList[1] = { prog };
      glShaderSource(_fragmentShader,1,progList,NULL);
      glCompileShader(_fragmentShader);

      GLint compile_status = 0;
      glGetShaderiv(_fragmentShader,GL_COMPILE_STATUS,&compile_status);
      if(compile_status == GL_FALSE) //if there was a compile error...
        {
          GLint log_length = 0;
          glGetShaderiv(_fragmentShader,GL_INFO_LOG_LENGTH,&log_length);
          if(log_length > 0)
            {
              string log_string(log_length,0);
              glGetShaderInfoLog(_fragmentShader,log_length,NULL,&(log_string[0]));
              throw ShaderCompileException("Error compiling fragment shader: " + log_string);
            }
          else
            throw ShaderCompileException("Error compiling fragment shader");
        }

      glAttachShader(_program, _fragmentShader);
      glLinkProgram(_program);

      GLint link_status = 0;
      glGetProgramiv(_program,GL_LINK_STATUS,&link_status);
      if(link_status == GL_FALSE) //if there was a link error
        {
          GLint log_length = 0;
          glGetProgramiv(_program,GL_INFO_LOG_LENGTH,&log_length);
          if(log_length > 0)
            {
              string log_string(log_length,0);
              glGetProgramInfoLog(_program,log_length,NULL,&(log_string[0]));
              throw ShaderLinkException("Error linking GLSL program: " + log_string);
            }
          else
            throw ShaderLinkException("Error linking GLSL program");
        }
    }
  else
    {
      _fragmentShader = glCreateShaderObjectARB(GL_FRAGMENT_SHADER_ARB);

      const char *progList[1] = { prog };
      glShaderSourceARB(_fragmentShader,1,progList,NULL);
      glCompileShaderARB(_fragmentShader);

      GLint compile_status = 0;
      glGetObjectParameterivARB(_fragmentShader,GL_OBJECT_COMPILE_STATUS_ARB,&compile_status);
      if(compile_status == GL_FALSE)
        {
          GLint log_length = 0;
          glGetObjectParameterivARB(_fragmentShader,GL_OBJECT_INFO_LOG_LENGTH_ARB,&log_length);
          if(log_length > 0)
            {
              string log_string(log_length,0);
              glGetInfoLogARB(_fragmentShader,log_length,NULL,&(log_string[0]));
              throw ShaderCompileException("Error compiling fragment shader: " + log_string);
            }
          else
            throw ShaderCompileException("Error compiling fragment shader");
        }

      glAttachObjectARB(_program, _fragmentShader);
      glLinkProgramARB(_program);

      GLint link_status = 0;
      glGetObjectParameterivARB(_program,GL_OBJECT_LINK_STATUS_ARB,&link_status);
      if(link_status == GL_FALSE)
        {
          GLint log_length = 0;
          glGetObjectParameterivARB(_program,GL_OBJECT_INFO_LOG_LENGTH_ARB,&log_length);
          if(log_length > 0)
            {
              string log_string(log_length,0);
              glGetInfoLogARB(_program,log_length,NULL,&(log_string[0]));
              throw ShaderCompileException("Error linking GLSL program: " + log_string);
            }
          else
            throw ShaderCompileException("Error linking GLSL program");
        }
    }

  updateUniforms();
}

void FragmentProgramGLSLImpl::unsetVertexShader()
{
  if(!initialized())
    throw UninitializedException("Could not un-set vertex shader!");
    
    
  if(_isGL2_0)
    {
      if(glIsShader(_vertexShader))
        {
          glDetachShader(_program, _vertexShader);
          glDeleteShader(_vertexShader);
        }
    }
  else
    {
      glDetachObjectARB(_program, _vertexShader);
      glDeleteObjectARB(_vertexShader);
    }
  _vertexShader = 0;
}

void FragmentProgramGLSLImpl::unsetFragmentShader()
{
  if(!initialized())
    throw UninitializedException("Could not un-set fragment shader!");
        
  if(_isGL2_0)
    {
      if(glIsShader(_fragmentShader))
        {
          glDetachShader(_program, _fragmentShader);
          glDeleteShader(_fragmentShader);
        }
    }
  else
    {
      glDetachObjectARB(_program, _fragmentShader);
      glDeleteObjectARB(_fragmentShader);
    }
  _fragmentShader = 0;
}

void FragmentProgramGLSLImpl::updateUniforms()
{
  using namespace boost; using boost::format;
  using namespace std;
 
  if(!initialized())
    throw UninitializedException("renderer not yet initialized");
  
  if(_isGL2_0)
    {
      //make sure the program is ok
      glValidateProgram(_program);
      GLint status;
      glGetProgramiv(_program,GL_VALIDATE_STATUS,&status);
      if(status == GL_FALSE)
        throw InvalidProgramException("Error updating uniforms: " + getGLSLLogString());          
      
      //Make _program the current program.
      //glUniform*() calls apply to whatever the current program is
      glUseProgram(_program);
        
      //update light
      GLint light_loc = glGetUniformLocation(_program,"lightDir");
      if(light_loc != -1) //if the program has a lightDir...
        {
          if(_lights.empty())
            {
              std::cerr<<BOOST_CURRENT_FUNCTION<<": no lights defined"<<std::endl;
            }
          else
            {
              //output some debug info about the uniform lightDir
              printUniformInfo("lightDir");

              glUniform3fv(light_loc,1,_lights[0].c_array());
              GLenum err = glGetError();
              if(err == GL_INVALID_OPERATION)
                throw InvalidProgramException("Error updating lightDir: " + getGLSLLogString());
              
              GLfloat u_light[3];
              glGetUniformfv(_program,light_loc,u_light);
              cerr<<str(format("Uniform var lightDir values: (%f,%f,%f)\n")
                        % u_light[0]
                        % u_light[1]
                        % u_light[2]);
            }
        }
      else
        cerr<<BOOST_CURRENT_FUNCTION<<": program doesn't have a lightDir"<<endl;

      //update sampler texture ids

      glUseProgram(0); //turn off _program to return state to what's expected
    }
  else
    {
      cerr<<BOOST_CURRENT_FUNCTION<<": Cannot update uniforms: only opengl 2.0 support for now"<<endl;
    }
}

void FragmentProgramGLSLImpl::printUniformInfo(const std::string& uniformName)
{
  using namespace std;
  using namespace boost; using boost::format;

  if(!initialized())
    throw UninitializedException("renderer not yet initialized");

  if(!_isGL2_0)
    throw InsufficientHardwareException("opengl 2.0 support for now");
  
  GLint uniform_id = glGetUniformLocation(_program,uniformName.c_str());
  if(uniform_id == -1)
    throw InvalidProgramException("program doesn't have a uniform called " + uniformName);

  GLint max_length, actual_length, uniform_size;
  GLenum uniform_type;
  glGetProgramiv(_program,GL_ACTIVE_UNIFORM_MAX_LENGTH,&max_length);
  if(max_length == 0)
    throw InvalidProgramException("no active uniform variables exist!");
  string namebuf(max_length,0);
  glGetActiveUniform(_program,
                     uniform_id,
                     max_length,
                     &actual_length,
                     &uniform_size,
                     &uniform_type,
                     &(namebuf[0]));
  namebuf.resize(actual_length);
  cerr<<str(format("Uniform var: %s\n\tsize: %d\n\ttype: %s\n")
            % namebuf.c_str()
            % uniform_size
            % _GLSL_UNIFORM_TYPE_MAP[uniform_type].c_str());
}

std::string FragmentProgramGLSLImpl::getGLSLLogString() const
{
  using namespace std;

  if(!initialized())
    throw UninitializedException("renderer not yet initialized");

  if(!_isGL2_0)
    throw InsufficientHardwareException("opengl 2.0 support for now");

  GLint log_length = 0;
  glGetProgramiv(_program,GL_INFO_LOG_LENGTH,&log_length);
  if(log_length > 0)
    {
      string log_string(log_length,0);
      glGetProgramInfoLog(_program,log_length,NULL,&(log_string[0]));
      return log_string;
    }

  return string();
}
