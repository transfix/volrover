/*
  Copyright 2002-2003 The University of Texas at Austin
  
	Authors: Anthony Thane <thanea@ices.utexas.edu>
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

#include <SimpleExample/SimpleExample.h>
#include <stdio.h>
#include <stdlib.h>

// header file for the volume renderer
#include <VolumeRenderer/VolumeRenderer.h>
#include <GL/glut.h>

/* Global Variables */
int W=640;
int H=480;

int mousex;
int mousey;
int cdx = 0;
int cdy = 0;
int whichbutton;
int state = ROTATE;

VolumeRenderer volumeRenderer;

GLdouble rotMatrix[16];
GLdouble transMatrix[16];

/// the following 3 functions show the interaction with the volume renderer

void InitVolumeRenderer() {
	if (!volumeRenderer.initRenderer()) {
		printf("Warning, there was an error initializing the volume renderer\n");
	}
}

void InitData() {
	// load the colormap from a file
	unsigned char byte_map[256*4];
	FILE* fp;
	// load the data from a file
	fp = fopen( "vh4.rawiv", "rb" );
	if (!fp) {
		printf("There was an error loading the data file\n");
		exit(1);
	}
	unsigned char* data = new unsigned char[128*128*128];
	unsigned char* rgbaData = new unsigned char[128*128*128*4];
	// skip the header
	fread(data, sizeof(unsigned char), 68, fp);

	printf("Reading the data file\n");
	// load the data
	fread(data, sizeof(unsigned char), 128*128*128, fp);
	fclose(fp);

	fp = fopen( "colormap.map", "rb" );
	if (!fp) {
		printf("There was an error loading the color map\n");
		exit(1);
	}
	fread(byte_map, sizeof(unsigned char), 256*4, fp);
	fclose(fp);

	unsigned int c, i;
	for (c = 0; c<128*128*128; c++) {
		i = data[c];
		rgbaData[c*4+0] = byte_map[i*4+0];
		rgbaData[c*4+1] = byte_map[i*4+1];
		rgbaData[c*4+2] = byte_map[i*4+2];
		rgbaData[c*4+3] = byte_map[i*4+3];
	}

	// send the data to the volume renderer
	if (!volumeRenderer.uploadColorMappedData(data, 128, 128, 128)) {
		// woops, maybe color mapped data wasnt supported

	}
	volumeRenderer.uploadColorMap(byte_map);
	//volumeRenderer.uploadRGBAData(rgbaData, 128, 128, 128);
	

	// commenting this out since shaded rendering isnt finished yet
	//printf("Calculating gradients\n");
	//volumeRenderer.calculateGradientsFromDensities(data, 128, 128, 128);

	delete [] data;


	printf("done uploading volume data\n");



}

void Display() {
	
	static int accumCnt=0;
	accumCnt++;
	
	// Clear the buffer
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	
	
	// Set up the modelview matrix
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glTranslated(0.0, 0.0, -4.0);
	glMultMatrixd(transMatrix);
	glMultMatrixd(rotMatrix);
	

	// the volume is rendered
	volumeRenderer.renderVolume();

	//glAccum(GL_MULT, 0.25);
	//glAccum(GL_ADD, 0.1);
	/*glAccum(GL_ACCUM, 0.1);
	glAccum(GL_RETURN, 1.0);

	if (accumCnt == 20) {
		glClear(GL_ACCUM_BUFFER_BIT);
		accumCnt = 0;
	}
	else if (accumCnt == 1) {
		glAccum(GL_LOAD, 1.0);
	}*/
		
		
	glutSwapBuffers();
}


int BeginGraphics (int* argc, char** argv, const char* name) {
	int win;
	
	glutInit(argc,argv);
	glutInitWindowSize(W,H);
	glutInitWindowPosition(100,100);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
	win = glutCreateWindow(name);
	glutSetWindow(win);
	init();
	
	RegisterCallbacks();
	
	glutMainLoop();
	return 0;
}

void RegisterCallbacks() {
	glutDisplayFunc(Display);
	glutIdleFunc(Idle);
	glutMouseFunc(MouseButton);
	glutMotionFunc(MouseMove);
	glutKeyboardFunc(Keyboard);
}

void init() {
	InitProjection();
	InitModelMatrix();
	InitLighting();
	InitState();
	InitVolumeRenderer();
	InitData();
	
	
}

void InitProjection() {
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective( 30.0,
		(GLdouble) W / (GLdouble) H,
		0.01,
		20.0
		);
}

void InitModelMatrix() {
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	/* initialize rotMatrix and transMatrix to the identity */
	glGetDoublev(GL_MODELVIEW_MATRIX, rotMatrix);
	glGetDoublev(GL_MODELVIEW_MATRIX, transMatrix);
}

void InitLighting() {
	GLfloat lightdir[] = {-1.0, 1.0, 1.0, 0.0};
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glLightfv(GL_LIGHT0, GL_POSITION, lightdir);
	glEnable(GL_NORMALIZE);
}

void InitState() {
	glClearColor(0.2, 0.2, 0.2, 1.0);
	glClearAccum(0.0,0.0,0.0,0.0);
	glColor4d(1.0, 1.0, 1.0, 1.0);
	glPolygonMode( GL_FRONT_AND_BACK, GL_FILL );
	glEnable(GL_DEPTH_TEST);
	glDisable(GL_CULL_FACE);
}

void Idle() {
	Rotate(cdx, cdy);
	glutPostRedisplay();
}

/* Callback function for the mouse button */
void MouseButton(int button, int mstate, int x, int y) {
	/* sets the initial positions for mousex and mousey */
	if (mstate == GLUT_DOWN) {
		mousex = x;
		mousey = y;
		whichbutton = button;
		cdx = 0;
		cdy = 0;
	}
	
}

/* Callback function for the mouse motion */
void MouseMove(int x, int y) {
	int dx = x-mousex;
	int dy = y-mousey;
	
	switch (state) {
	case ROTATE:
		Rotate(dx, dy);
		break;
	case TRANSLATE:
		if (whichbutton == GLUT_LEFT_BUTTON)
			Translate(dx, -dy,0);
		else
			Translate(0,0,-dy);
		break;
	}
	mousex = x;
	mousey = y;
	glutPostRedisplay();
}

void Keyboard(unsigned char key, int x, int y) {
	
	switch (key) {
	case 'r':
	case 'R':
		state = ROTATE;
		break;
	case 't':
	case 'T':
		state = TRANSLATE;
		break;
	}
}


void Rotate(int dx, int dy) {
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();
	glRotated(((double)dy), 1.0, 0.0, 0.0);
	glRotated(((double)dx), 0.0, 1.0, 0.0);
	glMultMatrixd(rotMatrix);
	glGetDoublev(GL_MODELVIEW_MATRIX, rotMatrix);
	glPopMatrix();
}

void Translate(int dx, int dy, int dz) {
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();
	glTranslated(((double)dx)/100.0, ((double)dy)/100.0, ((double)dz)/10.0);
	glMultMatrixd(transMatrix);
	glGetDoublev(GL_MODELVIEW_MATRIX, transMatrix);
	glPopMatrix();
}


int main(int argc, char** argv) {
	BeginGraphics(&argc, argv, "Simple Example");
	
	return 0;
}
