/*
  Copyright 2000-2002 The University of Texas at Austin

	Authors: Sanghun Park <hun@ices.utexas.edu>
        Advisor: Chandrajit Bajaj <bajaj@cs.utexas.edu>

  This file is part of volren.

  volren is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License version 2.1 as published by the Free Software Foundation.

  volren is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/

#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>

#include <volren/vr.h>

ImageBuffer* vrRayTracing(MultiVolRenEnv* menv,int quiet)
{
  ImageBuffer* img;
  Viewing* view;
  int ntx, nty, tid;
  unsigned char tile_image[TILE_RES*RGBA_SIZE];

  view = menv->env[0]->view;

  img = (ImageBuffer *)calloc(1, sizeof(ImageBuffer));
  vrCreateImage(img, view->pix_width, view->pix_height);
  ntx = img->width / TILE_SIZE; nty = img->height / TILE_SIZE;

  tid = 0;
  while (tid < ntx*nty) {
    vrTracingTile(menv, tid, tile_image);
    vrCopyTile(img, tid, tile_image);
    if(!quiet)
      printf("tile %d \n", tid);
    fflush(stdout);
    tid ++;
  }

  return img;
}

void vrTracingTile(MultiVolRenEnv *menv, int tid, unsigned char* tile_image)
{
  int i, j, n_tile_x, point[2], offset;
  Viewing* view;
	
  view = menv->env[0]->view; /* use the first volume's viewing parameters */

  n_tile_x = view->pix_width / TILE_SIZE;
  point[0] = (tid % n_tile_x) * TILE_SIZE;
  point[1] = (tid / n_tile_x) * TILE_SIZE;
  offset = 0;
  for (j = 0; j < TILE_SIZE; j++) {
    for (i = 0; i < TILE_SIZE; i++, offset++) {
      vrFireOneRay(menv, point[0]+i, point[1]+j, tile_image, offset);
    }
  }
}

void vrFireOneRay(MultiVolRenEnv* menv, int nx, int ny, unsigned char* tile_image, int offset)
{
  VolRenEnv *env;
  Volume* vol;
  Rendering* rend;
  Viewing* view;
  RayMisc* misc;
  Point3d ray_org, start_pnt, end_pnt;
  float min_x, min_y, min_z, max_x, max_y, max_z;
  int rt, i, j, curidx[MAX_ENV][3] /* the current cell */;
  float opac, r, g, b;
  unsigned char rendflag, *ucptr;
  float den;
  //float s_x, s_y, s_z;
  float orig_x, orig_y, orig_z;
  float span_x, span_y, span_z;
  float pnt[3], spn[3];
  float color[3] = {0, 0, 0};
  float new_opac, ratio;
  Vector3d norm;
  iRay ray;
  Cell cell[MAX_ENV];
  float min_val, max_val, vals[8];

  env = menv->env[0]; /* use the viewing parameters for the first volume */
  rend = env->rend;
  view = env->view;
  misc = env->misc;
  vol = env->vol;
  
  vrGetPixCoord(view, nx, ny, &ray_org);
  if (view->persp) {
    /* ray direction is toward eye */
    view->raydir.x = view->eye.x - ray_org.x;
    view->raydir.y = view->eye.y - ray_org.y;
    view->raydir.z = view->eye.z - ray_org.z;
    vrNormalize(&(view->raydir));
    rt = vrComputeIntersection(menv->metavol, &(view->eye), &(view->raydir), &start_pnt, &end_pnt);
  } else {
    rt = vrComputeIntersection(menv->metavol, &ray_org, &(view->vpn), &start_pnt, &end_pnt);
  }

  if (!rt) {		 /* not intersected */
    for (i = 0; i < 3; i++) {
      *(tile_image+RGBA_SIZE*offset+i) =  misc->back_color[i];
    }
    *(tile_image+RGBA_SIZE*offset+3) = 0;
    return;
  }
	
  vrSetRayMisc(env);
	
  /* Accumulate along the ray */
  min_x = menv->metavol->minb.x; min_y = menv->metavol->minb.y; min_z = menv->metavol->minb.z;
  max_x = menv->metavol->maxb.x; max_y = menv->metavol->maxb.y; max_z = menv->metavol->maxb.z;
  opac = 0.0f;
  //curidx[0] = curidx[1] = curidx[2] = -1;
  for(i=0; i<MAX_ENV; i++)
    for(j=0; j<3; j++)
      curidx[i][j] = -1;
  r = g = b = 0.0f;
  spn[0] = misc->step_x; spn[1] = misc->step_y; spn[2] = misc->step_z; /* Every iteration, step this much in this direction */
  pnt[0] = start_pnt.x; pnt[1] = start_pnt.y; pnt[2] = start_pnt.z; /* point on the subvolume in which we're starting */

  for(i=0; i<menv->n_env; i++)
    {
      if (ISOSURFACE(menv->env[i]->rend->render_mode))
	{
	  vrSetContourRay(&ray, view, pnt);
	  cell[i].orig[0] = vol->minb.x;
	  cell[i].orig[1] = vol->minb.y;
	  cell[i].orig[2] = vol->minb.z;
	  cell[i].span[0] = vol->span[0];
	  cell[i].span[1] = vol->span[1];
	  cell[i].span[2] = vol->span[2];
	}
    }
	
  while (opac < THRESHOLD_OPC && /* while we're below the opacity threshold and within the subvolume... */
	 min_x <= pnt[0] && pnt[0] <= max_x && min_y <= pnt[1] && pnt[1] <= max_y &&
	 min_z <= pnt[2] && pnt[2] <= max_z) {
    int idx[MAX_ENV][3];
	
    for(i=0; i<menv->n_env; i++)
      {
	env = menv->env[i];
	vol = env->vol;
	rend = env->rend;

	//orig_x = vol->orig.x; orig_y = vol->orig.y; orig_z = vol->orig.z;
	orig_x = vol->minb.x; orig_y = vol->minb.y; orig_z = vol->minb.z;
	span_x = vol->span[0]; span_y = vol->span[1]; span_z = vol->span[2];

	idx[i][0] = (int)((pnt[0] - orig_x) / span_x); /* the volume and cell we're in */
	idx[i][1] = (int)((pnt[1] - orig_y) / span_y);
	idx[i][2] = (int)((pnt[2] - orig_z) / span_z);
	  
	//printf("idx[i] == { %d, %d, %d } ",idx[i][0],idx[i][1],idx[i][2]);
	//printf("spn[0] == %f, spn[1] == %f, spn[2] == %f \n",spn[0],spn[1],spn[2]);

	new_opac = 0;

	if (idx[i][0] < 0 || idx[i][0] > vol->sub_ext[0]-2 || 
	    idx[i][1] < 0 || idx[i][1] > vol->sub_ext[1]-2 ||
	    idx[i][2] < 0 || idx[i][2] > vol->sub_ext[2]-2)
	  continue;

	if (idx[i][0] != curidx[i][0] || idx[i][1] != curidx[i][1] || idx[i][2] != curidx[i][2]) {
	  /* we're in a new cell */
	  
	  if (ISOSURFACE(rend->render_mode) || COLDEN(rend->render_mode)) {
	    vrGetVertDensities(vol, idx[i], vals);
	    min_val = max_val = vals[0];
	    for (j = 1; j < 8; j++) {
	      min_val = MIN2(min_val, vals[j]);
	      max_val = MAX2(max_val, vals[j]);
	    }
	    if (ISOSURFACE(rend->render_mode)) {
	      cell[i].id[0] = idx[i][0];
	      cell[i].id[1] = idx[i][1];
	      cell[i].id[2] = idx[i][2];
					
	      cell[i].func[0] = vals[0]; cell[i].func[1] = vals[1];
	      cell[i].func[2] = vals[5]; cell[i].func[3] = vals[4];
	      cell[i].func[4] = vals[2]; cell[i].func[5] = vals[3];
	      cell[i].func[6] = vals[7]; cell[i].func[7] = vals[6]; 
	      for (j = 0; j < rend->n_surface; j++) {
		if (rend->surf[j].value <= max_val && rend->surf[j].value >= min_val) {
		  vrComputeColorOpacitySurf(env, pnt, idx[i], j, &ray, &(cell[i]), vals, color, &new_opac);

		  ratio = new_opac * (1-opac); /* opacity increases as we go further into the subvolume, */
		  r += color[0] * ratio;       /* so color contribution decreases */
		  g += color[1] * ratio;
		  b += color[2] * ratio;
		  opac += ratio;  
		}
	      }
	    }
	    if (COLDEN(rend->render_mode)) {
	      if(min_val <= rend->max_den && max_val >= rend->min_den) {
		vrComputeColorOpacityInterp(env, pnt, idx[i], vals, color, &new_opac);
		
		ratio = new_opac * (1-opac);
		r += color[0] * ratio;
		g += color[1] * ratio;
		b += color[2] * ratio;
		opac += ratio;
	      }
	    }
	  }

	  if (RAYSHADING(rend->render_mode)) {
	    rendflag = 0x01;
	    vrGetDenNorm(env, pnt, idx[i], &den, &norm, &rendflag, vals);
	    if (rendflag) {
	      vrComputeColorOpacity(env, pnt, idx[i], den, &norm, vals, color, &new_opac);
	      assert(new_opac < 1.1);

	      ratio = new_opac * (1-opac);
	      r += color[0] * ratio;
	      g += color[1] * ratio;
	      b += color[2] * ratio;
	      opac += ratio;              
	      //printf("den: %f, new_opac: %f, opac: %f\n",den, new_opac, opac);
	    }
	  }

	  curidx[i][0] = idx[i][0]; 
	  curidx[i][1] = idx[i][1]; 
	  curidx[i][2] = idx[i][2];


	}
	/*else {
	// use the same value in the cell
	ratio = new_opac * (1-opac);
	r += color[0] * ratio;
	g += color[1] * ratio;
	b += color[2] * ratio;
	opac += ratio;
	}*/
      }

    pnt[0] += spn[0];   pnt[1] += spn[1];   pnt[2] += spn[2];
  }

  if (opac > THRESHOLD_OPC) {
    float invert;
    
    invert = 1.0f / opac;
    r *= invert;
    g *= invert;
    b *= invert;
  } else {
    float rest;

    rest = 1.0f - opac;
    r += misc->back_color[0] * rest;
    g += misc->back_color[1] * rest;
    b += misc->back_color[2] * rest;
  }

  ucptr = tile_image + RGBA_SIZE*offset;
  *ucptr = (unsigned char)MIN2(MAX_COLOR, r);
  *(ucptr+1) = (unsigned char)MIN2(MAX_COLOR, g);
  *(ucptr+2) = (unsigned char)MIN2(MAX_COLOR, b);
  *(ucptr+3) = (unsigned char)(opac*255);
}

