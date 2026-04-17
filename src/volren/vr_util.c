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

#include <math.h>
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#ifdef WIN32
#else
#include <sys/time.h>
#endif

#include <volren/vr.h>

float vrNormalize(Vector3d *vector)
{
  float len;

  len = (float) sqrt(vector->x*vector->x + vector->y*vector->y + vector->z*vector->z);
  if (VERYSMALL(len)) {
    vector->x = vector->y = vector->z = 0;
  } else {
    vector->x /= len;
    vector->y /= len;
    vector->z /= len;
  }

  return len;
}

void vrCopyShading(Shading* dest, Shading* src)
{
  dest->ambient.r = src->ambient.r;
  dest->ambient.g = src->ambient.g;
  dest->ambient.b = src->ambient.b;
  dest->diffuse.r = src->diffuse.r;
  dest->diffuse.g = src->diffuse.g;
  dest->diffuse.b = src->diffuse.b;
  dest->specular.r = src->specular.r;
  dest->specular.g = src->specular.g;
  dest->specular.b = src->specular.b;
  dest->shining = src->shining;
}

void vrGetPixCoord(Viewing* view, int nx, int ny, Point3d* pnt)
{
  float win_pix_x, win_pix_y;

  win_pix_x = view->pix_sz_x;    win_pix_y = view->pix_sz_y;

  pnt->x = view->win_sp.x + (float)nx * win_pix_x * view->vpu.x
    + (float)ny * win_pix_y * view->vup.x;
  pnt->y = view->win_sp.y + (float)nx * win_pix_x * view->vpu.y
    + (float)ny * win_pix_y * view->vup.y;
  pnt->z = view->win_sp.z + (float)nx * win_pix_x * view->vpu.z
    + (float)ny * win_pix_y * view->vup.z;
}

int vrComputeIntersection(Volume* vol, Point3d* eye_p, Vector3d* dir_p, 
			  Point3d* start_pnt_p, Point3d* end_pnt_p)
{
  //Volume* vol;
  float   x1, y1, z1, xdir, ydir, zdir;
  float   xb1, yb1, zb1, xb2, yb2, zb2;
  float   t1x, t1y, t1z, t2x, t2y, t2z, tmp;
  float   first_px, first_py, first_pz, second_px, second_py, second_pz;
  char    tnear, tfar;
  float   max_near, min_far;

  //vol = env->vol;

  x1   = eye_p->x;        y1   = eye_p->y;        z1   = eye_p->z; /* eye point */
  xb1  = vol->minb.x;     yb1  = vol->minb.y;     zb1  = vol->minb.z; /* subvolume lower bound */
  xb2  = vol->maxb.x;     yb2  = vol->maxb.y;     zb2  = vol->maxb.z; /* subvolume upper bound */
  xdir = -(dir_p->x);     ydir = -(dir_p->y);     zdir = -(dir_p->z); /* ray direction (dir_p) is toward eye */
  if (VERYSMALL(xdir)) {
    if (x1 < xb1 || x1 > xb2) return(0); /* the ray misses the subvolume */
    xdir = SMALLEST;
    t1x = -LARGEST; t2x = LARGEST;
  } else {
    first_px = xb1;     second_px = xb2;
    t1x = (xb1 - x1) / xdir;    t2x = (xb2 - x1) / xdir; /* max t values for the intersection of the subvolume */

    if (t1x > t2x) {
      tmp = t1x;  t1x = t2x;  t2x = tmp;
      first_px = xb2;     second_px = xb1;
    }
  }

  if (VERYSMALL(ydir)) {
    if (y1 < yb1 || y1 > yb2) return(0);
    ydir = SMALLEST;
    t1y = -LARGEST; t2y = LARGEST;
  } else {
    first_py = yb1;     second_py = yb2;
    t1y = (yb1 - y1) / ydir;    t2y = (yb2 - y1) / ydir;

    if (t1y > t2y) {
      tmp = t1y;  t1y = t2y;  t2y = tmp;
      first_py = yb2;     second_py = yb1;
    }
  }

  if (VERYSMALL(zdir)) {
    if (z1 < zb1 || z1 > zb2) return(0);
    zdir = SMALLEST;
    t1z = -LARGEST; t2z = LARGEST;
  } else {
    first_pz = zb1;     second_pz = zb2;
    t1z = (zb1 - z1) / zdir;     t2z = (zb2 - z1) / zdir; 

    if (t1z > t2z) {
      tmp = t1z;      t1z = t2z;      t2z = tmp;
      first_pz = zb2;         second_pz = zb1;
    }
  }

  /* Calaulate Max t_near	*/   
  if (t1x < t1y) {
    if (t1y < t1z) {
      tnear = 'z'; max_near = t1z;
    } else {
      tnear = 'y'; max_near = t1y;
    }
  } else {
    if (t1x < t1z) {
      tnear = 'z'; max_near = t1z;
    } else {
      tnear = 'x'; max_near = t1x;     
    }
  }

  /* Calculate Min t_far */
  if (t2x > t2y) {
    if (t2y > t2z) {
      tfar = 'z';  min_far = t2z;
    } else {
      tfar = 'y';  min_far = t2y;
    }
  } else {
    if (t2x > t2z) {
      tfar = 'z';  min_far = t2z;
    } else {
      tfar = 'x';  min_far = t2x;
    }
  }

  /* If Don't Intersect, Return 0 */
  if (max_near < 0.0)	return(0);
  if (max_near > min_far)	return(0);

  /* Calculate Intersection Points for Each Plane */
  switch (tnear) {
  case 'x' :  vrCalcIntersectTx(eye_p, dir_p, first_px, start_pnt_p);
    break;
  case 'y' :  vrCalcIntersectTy(eye_p, dir_p, first_py, start_pnt_p);
    break;
  case 'z' :  vrCalcIntersectTz(eye_p, dir_p, first_pz, start_pnt_p);
    break;
  }

  switch (tfar) {
  case 'x' :  vrCalcIntersectTx(eye_p, dir_p, second_px, end_pnt_p);
    break;
  case 'y' :  vrCalcIntersectTy(eye_p, dir_p, second_py, end_pnt_p);
    break;
  case 'z' :  vrCalcIntersectTz(eye_p, dir_p, second_pz, end_pnt_p);
    break;
  }

  return(1);
}

void vrCalcIntersectTx(Point3d *ray_org, Vector3d *dir, float px, Point3d *tx)
{
  float           tmp;

  tmp = (px - ray_org->x) / dir->x;
  tx->x = px;
  tx->y = dir->y * tmp + ray_org->y;
  tx->z = dir->z * tmp + ray_org->z;
}

void vrCalcIntersectTy(Point3d *ray_org, Vector3d *dir, float py, Point3d *ty)
{
  float           tmp;

  tmp = (py - ray_org->y) / dir->y;
  ty->x = dir->x * tmp + ray_org->x;
  ty->y = py;
  ty->z = dir->z * tmp + ray_org->z;
}


void vrCalcIntersectTz(Point3d *ray_org, Vector3d *dir, float pz, Point3d *tz)
{
  float           tmp;

  tmp = (pz - ray_org->z) / dir->z;
  tz->x = dir->x * tmp + ray_org->x;
  tz->y = dir->y * tmp + ray_org->y;
  tz->z = pz;
}

void vrGetVertDensities(Volume* vol, int idx[3], float val[8])
{
  int i;
  int off[8];

  /* the offsets specify the verticies clockwise, front face first (front meaning lower z coordinate) */

  off[0] = 0;                         off[1] = 1;
  off[2] = vol->sub_ext[0];           off[3] = off[2] + 1;
  off[4] = vol->slc_size;             off[5] = off[4] + 1;
  off[6] = off[4] + vol->sub_ext[0];  off[7] = off[6] + 1;

  switch (vol->type) {
  case UCHAR:
    {
      unsigned char* ucptr = vol->den.uc_ptr + idx[2]*vol->slc_size + idx[1]*vol->sub_ext[0]+idx[0];
      for (i = 0; i < 8; i++) {
	val[i] = *(ucptr + off[i]);
      }
      break;
    }
  case USHORT:
    {
      unsigned short* usptr = vol->den.us_ptr + idx[2]*vol->slc_size + idx[1]*vol->sub_ext[0]+idx[0];
      for(i=0;i<8;i++)
	val[i] = *(usptr + off[i]);
      break;
    }
  case FLOAT:
    {
      float* fptr = vol->den.f_ptr + idx[2]*vol->slc_size + idx[1]*vol->sub_ext[0]+idx[0];
      assert(0);
      for (i = 0; i < 8; i++) {
	val[i] = *(fptr + off[i]);
      }
      break;
    }
  case RAWV:
    {
      int sizeOfDenVal = vrRawVSizeOf(vol->den.rawv_var_types[vol->den.den_var]);
      unsigned char *rawvptr;

      /* advance pointer to the cell we want in density data */
      rawvptr = (unsigned char *)vrRawVGetValue(vol,0,vol->den.den_var,idx);

#define RAWV_DEN(t) for(i=0;i<8;i++) val[i] = (float)*((t *)(rawvptr + off[i]*sizeOfDenVal))
      
      RAWV_CAST(vol->den.rawv_var_types[vol->den.den_var],RAWV_DEN);

      break;
    }
  default:
    fprintf(stderr, "ERROR: unsupported data type: %d\n", vol->type);
  }
}

float vrTriInterp(float w[3], float val[8])
{
  float tmp0, tmp1, tmp2;

  tmp0 = (1-w[2])*val[0] + w[2]*val[4];
  tmp1 = (1-w[2])*val[2] + w[2]*val[6];
  tmp2 = (1-w[1])*tmp0 + w[1]*tmp1;

  tmp0 = (1-w[2])*val[1] + w[2]*val[5];
  tmp1 = (1-w[2])*val[3] + w[2]*val[7];
  tmp0 = (1-w[1])*tmp0 + w[1]*tmp1;

  //if(tmp0 != 0 || tmp1 != 0 || tmp2 != 0)
  //  printf("tmp: %f %f %f\n",tmp0,tmp1,tmp2);

  return((1-w[0])*tmp2 + w[0]*tmp0);
}

void vrGetDenNorm(VolRenEnv* env, float pnt[3], int idx[3], float *den, Vector3d *norm, 
		  unsigned char *flag, float val[8])
{
  float w[3];
  Volume *vol;
  Rendering *rend;

  vol = env->vol;
  rend = env->rend;

  vrGetVertDensities(vol, idx, val);

  w[0] = (pnt[0] - (vol->minb.x + idx[0]*vol->span[0]))/vol->span[0];
  w[1] = (pnt[1] - (vol->minb.y + idx[1]*vol->span[1]))/vol->span[1];
  w[2] = (pnt[2] - (vol->minb.z + idx[2]*vol->span[2]))/vol->span[2];

  *den = vrTriInterp(w, val); 
  if (*den < rend->min_den || *den > rend->max_den) {
    *flag &= 0;		// 11111010 
    return;
  }

#ifdef FAST_NORMAL
  norm->x = (1-w[2])*((1-w[1])*(val[1]-val[0])+w[1]*(val[3]-val[2]))
    +w[2]*((1-w[1])*(val[5]-val[4])+w[1]*(val[7]-val[6]));

  norm->y = (1-w[2])*((1-w[0])*(val[2]-val[0])+w[0]*(val[3]-val[1]))
    +w[2]*((1-w[0])*(val[6]-val[4])+w[0]*(val[7]-val[5]));

  norm->z = (1-w[0])*((1-w[1])*(val[4]-val[0])+w[0]*(val[6]-val[2]))
    +w[0]*((1-w[1])*(val[5]-val[1])+w[0]*(val[7]-val[3]));
#else
  vrSplineNorm(env, w, idx, norm);
#endif
}

void vrSplineNorm(VolRenEnv* env, float w[3], int idx[3], Vector3d* norm)
{
  //static float val[4][4][4];
  //static int curidx[3] = {-1, -1, -1};
  // deBoor Helpers
  float delta;				  		// the weights of deBoor algo
  //static float Derivative[4][4][3][3];				// leaves of the deBoor tree
  float D[4][4][2];							// knots of the deBoor tree
  int   i, j, k, l, m;
  float normal[3];
  
  /* this is ugly but better than the alternative :( */
#define val env->val
#define curidx env->curidx
#define Derivative env->Derivative
  //float *val[][];
  //float 

  Volume *vol = env->vol;

  //val = env->val;
  //curidx = env->curidx;
  //Derivative = env->Derivative;

  if(!env->notfirst)
    {
      curidx[0] = -1; curidx[1] = -1; curidx[2] = -1;
      env->notfirst = 1;
    }

  if (idx[0] != curidx[0] || idx[1] != curidx[1] || idx[2] != curidx[2]) {
    curidx[0] = idx[0]; curidx[1] = idx[1]; curidx[2] = idx[2];
    for (k = -1; k <= 2; k++) {
      for (j = -1; j <= 2; j++) {
	for (i = -1; i <= 2; i++) {
	  int vi, vj, vk, n;
	  vi = MIN2(vol->sub_ext[0]-1, MAX2(idx[0]+i, 0));
	  vj = MIN2(vol->sub_ext[1]-1, MAX2(idx[1]+j, 0));
	  vk = MIN2(vol->sub_ext[2]-1, MAX2(idx[2]+k, 0));
	  //n = vk*vol->slc_size + vj*vol->dim[0] + vi;
	  n = vol->yinc[vk] + vol->xinc[vj] + vi;

	  switch (vol->type) {
	  case UCHAR: 
	    val[k+1][j+1][i+1] = vol->den.uc_ptr[n];
	    break;
	  case USHORT:
	    val[k+1][j+1][i+1] = vol->den.us_ptr[n];
	    break;
	  case FLOAT:
	    val[k+1][j+1][i+1] = vol->den.f_ptr[n];
	    break;
	  case RAWV:
#define RAWV_DEN_SPLINE(t) \
            { \
            unsigned char *rawvptr; \
            int _i; \
            rawvptr = vol->den.rawv_ptr; \
	    for(_i=0;_i<vol->den.rawv_num_vars-1;_i++) rawvptr+=vol->vol_size*vol->den.rawv_num_timesteps*vrRawVSizeOf(vol->den.rawv_var_types[_i]); \
	    val[k+1][j+1][i+1] = (float)*((t *)(rawvptr + n*vrRawVSizeOf(vol->den.rawv_var_types[vol->den.rawv_num_vars-1]))); \
            }

	    RAWV_CAST(vol->den.rawv_var_types[vol->den.rawv_num_vars-1],RAWV_DEN_SPLINE);
	    
	    break;
	  }
	}
      }
    }

    // d(main)
    for (k = 0; k < 3; k++)
      for (j = 1; j < 3; j++)
	for (i = 1; i < 3; i++)
	  Derivative[i][j][k][0] = val[i][j][k+1] - val[i][j][k];

    // d(main+1)
    for (k = 1; k < 3; k++)
      for (j = 0; j < 3; j++)
	for (i = 1; i < 3; i++)
	  Derivative[k][i][j][1] = val[i][j+1][k] - val[i][j][k];

    // d(main+2)
    for (k = 1; k < 3; k++)
      for (j = 1; j < 3; j++)
	for (i = 0; i < 3; i++)
	  Derivative[j][k][i][2] = val[i+1][j][k] - val[i][j][k];
  }

  // deBoor for the 3 Derivatives (degree=2x1x1)
  for (l = 0; l < 3; l++) {
    // collapse the derivative-direction first
    for (k = 1; k < 3; k++) {
      for (j = 1; j < 3; j++) {
	delta = 0.5f * (1.0f - w[l]);
	D[k][j][0] = delta * Derivative[k][j][0][l] + (1.0f-delta) * Derivative[k][j][1][l];
	delta = 0.5f * (2.0f - w[l]);
	D[k][j][1] = delta * Derivative[k][j][1][l] + (1.0f-delta) * Derivative[k][j][2][l];

	delta = 1.0f - w[l];
	D[k][j][0] = delta * D[k][j][0] + (1.0f-delta) * D[k][j][1];
      }
    }
    m = (l+1)%3;

    // now the next direction
    for (k = 1; k < 3; k++) {
      delta = 1.0f - w[m];
      D[k][0][0] = delta * D[k][1][0] + (1.0f-delta) * D[k][2][0];
    }

    m = (l+2)%3;

    // and finally the last one
    delta = 1.0f - w[m];
    normal[l] = delta * D[1][0][0] + (1.0f-delta) * D[2][0][0];
  }
  norm->x = normal[0];
  norm->y = normal[1];
  norm->z = normal[2];

#undef val
#undef curidx
#undef Derivative
}

void vrComputeColorOpacityInterp(VolRenEnv* env, float pnt[3], int idx[3], float val[8], 
				 float color[3], float *opac)
{
  float w[3], den;
  Volume *vol;
  Rendering *rend;
  unsigned char * cbptr;
  float colors[3][8];
  //Vector3d *norm;

  vol = env->vol;
  rend = env->rend;

  w[0] = (pnt[0] - (vol->minb.x + idx[0]*vol->span[0]))/vol->span[0];
  w[1] = (pnt[1] - (vol->minb.y + idx[1]*vol->span[1]))/vol->span[1];
  w[2] = (pnt[2] - (vol->minb.z + idx[2]*vol->span[2]))/vol->span[2];

  if(vol->type == RAWV)
    {
      vrRawVGetVertColor(vol,0,idx,colors);
      color[0] = vrTriInterp(w,colors[0]);
      color[1] = vrTriInterp(w,colors[1]);
      color[2] = vrTriInterp(w,colors[2]);
      den = vrTriInterp(w, val);
      //if(!VERYSMALL(den))
      //printf("den == %f\n",den);

      *opac = den / rend->max_den;
      //if(!VERYSMALL(*opac)) printf("*opac == %f\n",*opac);
      return;
    }
  else
    {
      den = vrTriInterp(w, val);

      cbptr = rend->coldentbl + (int)(den) * RGBA_SIZE;

      color[0] = *cbptr;
      color[1] = *(cbptr + 1);
      color[2] = *(cbptr + 2);
      
      *opac = *(cbptr + 3) / 255.0f; /* alpha */
    }
}

void vrComputeColorOpacity(VolRenEnv* env, float pnt[3], int idx[3], float den, Vector3d* norm,
			   float vals[8], float color[3], float* opac)
{
  int len;
  unsigned char * cbptr;
  Rendering *rend;
  RayMisc* misc;
  Shading shade;
  Volume *vol;
  //float dens;

  rend = env->rend;
  misc = env->misc;
  vol = env->vol;

  if(env->vol->type == RAWV)
    {
      float colors[3][8];
      float w[3];
      float c[3], avg, val;
      int i,j,k;
      float tmp0, tmp1, tmp2;
      
      vrRawVGetVertColor(vol,0,idx,colors);
      w[0] = (pnt[0] - (vol->minb.x + idx[0]*vol->span[0]))/vol->span[0];
      w[1] = (pnt[1] - (vol->minb.y + idx[1]*vol->span[1]))/vol->span[1];
      w[2] = (pnt[2] - (vol->minb.z + idx[2]*vol->span[2]))/vol->span[2];
      //c[0] = vrTriInterp(w,colors[0]);
      //c[1] = vrTriInterp(w,colors[1]);
      //c[2] = vrTriInterp(w,colors[2]);
      //den = vrTriInterp(w,val);

      /* do color interpolation only with cells that have a density greater than zero */
      val = 0.0;
      avg = c[0] = c[1] = c[2] = 0;
      for(k=0; k<3; k++)
	{
	  for(i=0,j=0,avg=0; i<8; i++)
	    if(vals[i] > val)
	      {
		avg += colors[k][i];
		j++;
	      }
	  if(j>0) avg /= j;
	  
	  tmp0 = (1-w[2])*(vals[0]>val?colors[k][0]:avg) + w[2]*(vals[4]>val?colors[k][4]:avg);
	  tmp1 = (1-w[2])*(vals[2]>val?colors[k][2]:avg) + w[2]*(vals[6]>val?colors[k][6]:avg);
	  tmp2 = (1-w[1])*tmp0 + w[1]*tmp1;
	  
	  tmp0 = (1-w[2])*(vals[1]>val?colors[k][1]:avg) + w[2]*(vals[5]>val?colors[k][5]:avg);
	  tmp1 = (1-w[2])*(vals[3]>val?colors[k][3]:avg) + w[2]*(vals[7]>val?colors[k][7]:avg);;
	  tmp0 = (1-w[1])*tmp0 + w[1]*tmp1;
	  
	  c[k] = (1-w[0])*tmp2 + w[0]*tmp0;
	}

      len = (int) vrNormalize(norm);
      if (len >= MAX_GRADIENT) len = MAX_GRADIENT-1;
      *opac = den * misc->gradtbl[len] / rend->max_den; 

      if (VERYSMALL(*opac)) {
	color[0] = color[1] = color[2] = 0;
      } else {
	shade.ambient.r = (unsigned char)c[0];
	shade.ambient.g = (unsigned char)c[1];
	shade.ambient.b = (unsigned char)c[2];
	shade.diffuse.r = (unsigned char)c[0];
	shade.diffuse.g = (unsigned char)c[1];
	shade.diffuse.b = (unsigned char)c[2];
	shade.specular.r = (unsigned char)c[0];
	shade.specular.g = (unsigned char)c[1];
	shade.specular.b = (unsigned char)c[2];
	shade.shining = 15;

	vrPhongShading(env, color, &shade, norm);
      }
    }
  else
    {
      /* find out material */
      cbptr = rend->coldentbl + (int)(den) * RGBA_SIZE;
      len = (int) vrNormalize(norm);
      if (len >= MAX_GRADIENT) len = MAX_GRADIENT-1;
      *opac = *(cbptr + 3) * misc->gradtbl[len] / 255.0f;

      /*
      if(den < rend->min_den || den > rend->max_den) {
	printf("den: %f, opac: %f, color : %d %d %d\n", den, *opac, *cbptr, *(cbptr+1), *(cbptr+2));
	}
      */
      if (VERYSMALL(*opac)) {
	color[0] = color[1] = color[2] = 0;
      } else {
	shade.ambient.r = *cbptr;
	shade.ambient.g = *(cbptr+1);
	shade.ambient.b = *(cbptr+2);
	shade.diffuse.r = *cbptr;
	shade.diffuse.g = *(cbptr+1);
	shade.diffuse.b = *(cbptr+2);
	shade.specular.r = *cbptr;
	shade.specular.g = *(cbptr+1);
	shade.specular.b = *(cbptr+2);
	shade.shining = 15;

	vrPhongShading(env, color, &shade, norm);
      }
    }
}

void vrPhongShading(VolRenEnv *env, float color[3], Shading* shade, Vector3d* norm)
{
  Vector3d  light_vec, view_vec, half_vec;
  Viewing *view;
  Rendering *rend;
  int i;
  float xl, xh;

  view = env->view;
  rend = env->rend;

  if (view->persp) {
    view_vec.x = view->raydir.x;
    view_vec.y = view->raydir.y;
    view_vec.z = view->raydir.z;
  } else {
    view_vec.x = view->vpn.x;
    view_vec.y = view->vpn.y;
    view_vec.z = view->vpn.z;
  }

  //color[0] = (float) shade->ambient.r; 
  //color[1] = (float) shade->ambient.g; 
  //color[2] = (float) shade->ambient.b;

  color[0] = 0; 
  color[1] = 0; 
  color[2] = 0;
  //return;		// no shading

  for (i = 0; i < rend->n_light; i++) {
    light_vec.x = rend->light[i].dir.x;
    light_vec.y = rend->light[i].dir.y;
    light_vec.z = rend->light[i].dir.z;

    half_vec.x = light_vec.x + view_vec.x;
    half_vec.y = light_vec.y + view_vec.y;
    half_vec.z = light_vec.z + view_vec.z;
    vrNormalize(&half_vec);

    xl = INNERPROD(*norm, light_vec);
    xh = INNERPROD(*norm, half_vec);

    if (rend->light_both) {
      xl = (float)fabs(xl);
      xh = (float)fabs(xh);
    } else {
      //xv = INNERPROD(*norm, view_vec);
      if (xl >= 0) {
	// xl = (xl > 0.0) ? xl : 0.0f;
	xh = (xh > 0.0) ? xh : 0.0f;
      } else {
	xl = 0.0f;
	xh = 0.0f;
      }
    }
    //if (xh != 0.0) xh = (float)pow(xh, shade->shining);
    xh = spec_val[(int)(xh*QUANTIZE_SIZE)];
    color[0] = shade->diffuse.r * xl * rend->light[i].color.r / 256 
      + rend->light[i].color.r * xh;
    color[1] = shade->diffuse.g * xl * rend->light[i].color.g / 256 
      + rend->light[i].color.g * xh;
    color[2] = shade->diffuse.b * xl * rend->light[i].color.g / 256 
      + rend->light[i].color.b * xh;
  }

  color[0] = MIN2(0.9f * color[0], 255);
  color[1] = MIN2(0.9f * color[1], 255);
  color[2] = MIN2(0.9f * color[2], 255);
}

void vrComputeColorOpacitySurf(VolRenEnv* env, float pnt[3], int idx[3], int n, iRay *ray, Cell* cell,  
			       float vals[8], float color[3], float* opac)
{
  int i;
  float val;
  float w[3];
  Vector3d normal;
  Surface* surf;
  Volume *vol;

  vol = env->vol;

  *opac = 0;
  color[0] = color[1] = color[2] = 0;

  if (n >= env->rend->n_surface) return;

  surf = &(env->rend->surf[n]);
  val = surf->value;
  if (!iso_intersectW(*ray, val, cell, w)) return;

  *opac = env->rend->surf[n].opacity;
  vrSplineNorm(env, w, cell->id, &normal);
  vrNormalize(&normal);
  if(env->vol->type == RAWV)
    {
      /*float c[3];
      vrRawVGetColor(env->vol,0,idx,c);
      if(c[0] == 0.0 && c[1] == 0.0 && c[2] == 0.0)
	{
	  *opac = 0;
	  color[0] = color[1] = color[2] = 0;
	  return;
	  }*/
      float colors[3][8];
      float c[3], avg;
      int j,k;

      float tmp0, tmp1, tmp2;

      //vrRawVGetColor(env->vol,0,idx,c);
      
      vrRawVGetVertColor(vol,0,idx,colors);
      w[0] = (pnt[0] - (vol->minb.x + idx[0]*vol->span[0]))/vol->span[0];
      w[1] = (pnt[1] - (vol->minb.y + idx[1]*vol->span[1]))/vol->span[1];
      w[2] = (pnt[2] - (vol->minb.z + idx[2]*vol->span[2]))/vol->span[2];
      //c[0] = vrTriInterp(w,colors[0]);
      //c[1] = vrTriInterp(w,colors[1]);
      //c[2] = vrTriInterp(w,colors[2]);

      /* do color interpolation only with cells that contribute to the surface */
      avg = c[0] = c[1] = c[2] = 0;
      for(k=0; k<3; k++)
	{
	  for(i=0,j=0,avg=0; i<8; i++)
	    if(vals[i] >= val)
	      {
		avg += colors[k][i];
		j++;
	      }
	  if(j>0) avg /= j;

	  tmp0 = (1-w[2])*(vals[0]>=val?colors[k][0]:avg) + w[2]*(vals[4]>=val?colors[k][4]:avg);
	  tmp1 = (1-w[2])*(vals[2]>=val?colors[k][2]:avg) + w[2]*(vals[6]>=val?colors[k][6]:avg);
	  tmp2 = (1-w[1])*tmp0 + w[1]*tmp1;
	  
	  tmp0 = (1-w[2])*(vals[1]>=val?colors[k][1]:avg) + w[2]*(vals[5]>=val?colors[k][5]:avg);
	  tmp1 = (1-w[2])*(vals[3]>=val?colors[k][3]:avg) + w[2]*(vals[7]>=val?colors[k][7]:avg);;
	  tmp0 = (1-w[1])*tmp0 + w[1]*tmp1;

	  c[k] = (1-w[0])*tmp2 + w[0]*tmp0;
	}
      

      surf->shading.ambient.r = (unsigned char)c[0];
      surf->shading.ambient.g = (unsigned char)c[1];
      surf->shading.ambient.b = (unsigned char)c[2];
      surf->shading.diffuse.r = (unsigned char)c[0];
      surf->shading.diffuse.g = (unsigned char)c[1];
      surf->shading.diffuse.b = (unsigned char)c[2];
      surf->shading.specular.r = (unsigned char)c[0];
      surf->shading.specular.g = (unsigned char)c[1];
      surf->shading.specular.b = (unsigned char)c[2];
      surf->shading.shining = 15;
    }
  vrPhongShading(env, color, &(surf->shading), &normal);
}

void vrCopyTile(ImageBuffer* img, int tid, unsigned char* tile)
{
  int ntx, i, j, k;
  int pos[2];
  unsigned char* buf;

  

  ntx = img->width / TILE_SIZE;
  buf = img->buffer;
  pos[0] = (tid % ntx);
  pos[1] = (tid / ntx);
  for (i = 0; i < TILE_SIZE; i++) {
    for (j = 0; j < TILE_SIZE; j++) {
      for (k = 0; k < RGBA_SIZE; k++)
	buf[RGBA_SIZE*((i+pos[1]*TILE_SIZE)*img->width + (j+pos[0]*TILE_SIZE)) + k] = tile[RGBA_SIZE*(i*TILE_SIZE + j) + k];
    } 
  } 
}

double vrGetTime()
{
#ifdef WIN32
  clock_t curtime = clock();
  return(double)(curtime) / CLOCKS_PER_SEC;
#else
  struct timeval t;
  gettimeofday(&t, NULL);
  return t.tv_sec + t.tv_usec / 1000000.0;
#endif
}

void *vrRawVGetValue(Volume *vol, int timestep, unsigned int var, int idx[3])
{
  unsigned int i,size=0;
  for(i=0; i<var; i++) size += vrRawVSizeOf(vol->den.rawv_var_types[i]);
  return vol->den.rawv_ptr + vol->den.rawv_num_timesteps*vol->vol_size*size + timestep*vol->vol_size*vrRawVSizeOf(vol->den.rawv_var_types[var]) + 
    (idx[2]*vol->slc_size + idx[1]*vol->dim[0] + idx[0])*vrRawVSizeOf(vol->den.rawv_var_types[var]);
}

void vrRawVGetColor(Volume *vol, int timestep, int idx[3], float color[3])
{
#define RAWV_R(t) color[0] = (float)*((t *)vrRawVGetValue(vol,timestep,vol->den.r_var,idx))
#define RAWV_G(t) color[1] = (float)*((t *)vrRawVGetValue(vol,timestep,vol->den.g_var,idx))
#define RAWV_B(t) color[2] = (float)*((t *)vrRawVGetValue(vol,timestep,vol->den.b_var,idx))

  RAWV_CAST(vol->den.rawv_var_types[vol->den.r_var],RAWV_R);
  RAWV_CAST(vol->den.rawv_var_types[vol->den.g_var],RAWV_G);
  RAWV_CAST(vol->den.rawv_var_types[vol->den.b_var],RAWV_B);
}

void vrRawVGetDensity(Volume *vol, int timestep, int idx[3], float *dens)
{
#define RAWV_DENS(t) *dens = (float)*((t *)vrRawVGetValue(vol,timestep,vol->den.den_var,idx))
  RAWV_CAST(vol->den.rawv_var_types[vol->den.den_var],RAWV_DENS);
}

void vrRawVGetVertColor(Volume *vol, int timestep, int idx[3], float color[3][8])
{
  int i;
  int off[8];
  unsigned char *rawvptr;
  int size_r, size_g, size_b;

  off[0] = 0;                     off[1] = 1;
  off[2] = vol->dim[0];           off[3] = off[2] + 1;
  off[4] = vol->slc_size;         off[5] = off[4] + 1;
  off[6] = off[4] + vol->dim[0];  off[7] = off[6] + 1;

  /* red */
  size_r = vrRawVSizeOf(vol->den.rawv_var_types[vol->den.r_var]);
  rawvptr = (unsigned char *)vrRawVGetValue(vol,timestep,vol->den.r_var,idx);
#define RAWV_VERT_R(t) for(i=0;i<8;i++) color[0][i] = (float)*((t *)(rawvptr + off[i]*size_r))
  RAWV_CAST(vol->den.rawv_var_types[vol->den.r_var],RAWV_VERT_R);
  /* green */
  size_g = vrRawVSizeOf(vol->den.rawv_var_types[vol->den.g_var]);
  rawvptr = (unsigned char *)vrRawVGetValue(vol,timestep,vol->den.g_var,idx);
#define RAWV_VERT_G(t) for(i=0;i<8;i++) color[1][i] = (float)*((t *)(rawvptr + off[i]*size_g))
  RAWV_CAST(vol->den.rawv_var_types[vol->den.g_var],RAWV_VERT_G);
  /* blue */
  size_b = vrRawVSizeOf(vol->den.rawv_var_types[vol->den.b_var]);
  rawvptr = (unsigned char *)vrRawVGetValue(vol,timestep,vol->den.b_var,idx);
#define RAWV_VERT_B(t) for(i=0;i<8;i++) color[2][i] = (float)*((t *)(rawvptr + off[i]*size_b))
  RAWV_CAST(vol->den.rawv_var_types[vol->den.b_var],RAWV_VERT_B);
}

MultiVolRenEnv *vrCopyMenv(MultiVolRenEnv *menv)
{
  int i;
  MultiVolRenEnv *ret=NULL;

  ret = (MultiVolRenEnv *)calloc(1,sizeof(MultiVolRenEnv));
  ret->n_env = menv->n_env;
  for(i=0; i<ret->n_env; i++)
    {
      ret->metavol = (Volume *)calloc(1,sizeof(Volume));
      memcpy(ret->metavol,menv->metavol,sizeof(Volume));
      ret->env[i] = vrCopyEnv(menv->env[i]);
    }

  return ret;
}

VolRenEnv *vrCopyEnv(VolRenEnv *env)
{
  VolRenEnv *ret = (VolRenEnv *)calloc(1,sizeof(VolRenEnv));

  /* This is not a deep copy, so be sure to only free this stuff in one thread */
  ret->vol = (Volume *)calloc(1,sizeof(Volume));
  memcpy(ret->vol,env->vol,sizeof(Volume));
  ret->view = (Viewing *)calloc(1,sizeof(Viewing));
  memcpy(ret->view,env->view,sizeof(Viewing));
  ret->rend = (Rendering *)calloc(1,sizeof(Rendering));
  memcpy(ret->rend,env->rend,sizeof(Rendering));
  ret->misc = (RayMisc *)calloc(1,sizeof(RayMisc));
  memcpy(ret->misc,env->misc,sizeof(RayMisc));

  return ret;
}
