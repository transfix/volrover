#include <cvc_compat.h>
#include <cvc_compat.h>
/***************************************************************************
 *   Copyright (C) 2009 by Bharadwaj Subramanian   *
 *   bharadwajs@axon.ices.utexas.edu   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
#include <VolMagick/VolMagick.h>
#include "Mesh.h"
#include <iostream>
#include <set>
#include <math.h>
#include "Grid.h"

#define EPSILON 0.02

void ReSamplingCubicSpline1OrderPartials(float *coeff, int nx, int ny, int nz, float *dxyz, float *minExtent, float *zeroToTwoPartialValue,  float *p);
void geometricFlow(Mesh &mesh);

int projectMeshToIsosurface(VolMagick::Volume *vol,Mesh *mesh,float isovalue)
{
  vector<VolMagick::Volume> gradient; // Remove when zqyork's gradient computation is ready.
  VolMagick::calcGradient(cvcapp, gradient,*vol);
//   float *coefficients;
//  int count=0;
//  int dim[3];
// 
//  dim[0]=coeff->XDim(); dim[1]=coeff->YDim(); dim[2]=coeff->ZDim();
// 
//   coefficients = new float[ coeff->XDim() * coeff->YDim() * coeff->ZDim() ];
/*
    int  index, index1;
    for(int p=0;p<coeff->XDim();p++)
        for(int q=0;q<coeff->YDim();q++)
            for(int r=0;r<coeff->ZDim();r++)
            {

               index1 = (r*dim[1]+q)*dim[0]+p;
               coefficients[index1]= (*coeff)(p,q,r);
            }*/


  int numverts=mesh->vertices.size();

  int numDirtyVerts=0;
  int fixed=0;

  for(int i=0;i<numverts;i++)
  {
    float x=mesh->vertices[i].x;
    float y=mesh->vertices[i].y;
    float z=mesh->vertices[i].z;
    
    float value=vol->interpolate(x,y,z);

    if(value<=isovalue*(1-EPSILON)||value>=isovalue*(1+EPSILON)) // We have a point away fromt he isosurface!
    {
      float newvalue;
      int itercount=0;
     // cout<<"Before:  "<<value<<" After: ";
      do
      {
        float d, dist=isovalue-value;
        float abs_grad_hls;
        
        float d_hls[10]; // The partial derivatives in x,y and z directions.
        float p[3];
        float dxyz[3];
        float minExt[3];
        
        d_hls[0]=value;
//         p[0]=x; dxyz[2]=coeff->XSpan(); minExt[0]=coeff->XMin();
//         p[1]=y; dxyz[1]=coeff->YSpan(); minExt[1]=coeff->YMin();
//         p[2]=z; dxyz[0]=coeff->ZSpan(); minExt[2]=coeff->ZMin();
        
        
        float absnormal=0;
            
        // Compute the absolute value of the gradient at that point.
         d_hls[0]=gradient[0].interpolate(x,y,z);
         d_hls[1]=gradient[1].interpolate(x,y,z);
         d_hls[2]=gradient[2].interpolate(x,y,z);
//         ReSamplingCubicSpline1OrderPartials(coefficients,coeff->XDim(),coeff->YDim(),coeff->ZDim(),dxyz, minExt, d_hls, p);

        abs_grad_hls=d_hls[0]*d_hls[0]+d_hls[1]*d_hls[1]+d_hls[2]*d_hls[2];
//        abs_grad_hls=d_hls[3]*d_hls[3]+d_hls[1]*d_hls[1]+d_hls[2]*d_hls[2];

        if(abs_grad_hls<=0.001) abs_grad_hls=1.0;

        d=dist/abs_grad_hls; // Get distance

        // Update the distance with the EPSILON.
        //d*=(1+EPSILON/2);

        //Move the point.
//         mesh->vertices[i].x+=  d * d_hls[3] ;
//         mesh->vertices[i].y+=  d * d_hls[2] ;
//         mesh->vertices[i].z+=  d * d_hls[1] ;
        mesh->vertices[i].x+=  d * d_hls[0] ;
        mesh->vertices[i].y+=  d * d_hls[1] ;
        mesh->vertices[i].z+=  d * d_hls[2] ;

       // cout<<"x: "<<mesh->vertices[i].x<<" y: "<<mesh->vertices[i].y<<" z: "<<mesh->vertices[i].z<<endl;

        // Bad idea. Delete.
//         if(mesh->vertices[i].x<vol->XMin()) mesh->vertices[i].x=vol->XMin();
//         if(mesh->vertices[i].y<vol->YMin()) mesh->vertices[i].y=vol->YMin();
//         if(mesh->vertices[i].z<vol->ZMin()) mesh->vertices[i].z=vol->ZMin();
//         if(mesh->vertices[i].x>vol->XMax()) mesh->vertices[i].x=vol->XMax();
//         if(mesh->vertices[i].y>vol->YMax()) mesh->vertices[i].y=vol->YMax();
//         if(mesh->vertices[i].z>vol->ZMax()) mesh->vertices[i].z=vol->ZMax();
        
        newvalue=vol->interpolate(mesh->vertices[i].x,mesh->vertices[i].y,mesh->vertices[i].z);

       // cout<<"Here2"<<endl;
        itercount++;

        if(newvalue>=isovalue*(1-EPSILON)&&newvalue<=isovalue*(1+EPSILON))
          fixed++;
        else
          value=newvalue;

        
      }while((newvalue<=isovalue*(1-EPSILON)||newvalue>=isovalue*(1+EPSILON))&&itercount<50);
   //   cout<<itercount<<" "<<newvalue<<endl;

      numDirtyVerts++;
    }
  }

  cout<<"NumDirtyVerts = "<<numDirtyVerts<<", fixed = "<<fixed<<endl;

//   delete[] coefficients;
  
  return numDirtyVerts-fixed;
}


int correctMeshByProjection(VolMagick::Volume *vol,Mesh *mesh,float isovalue,float tolerance,bool &dirty)
{
  vector<VolMagick::Volume> gradient; // Remove when zqyork's gradient computation is ready.
  VolMagick::calcGradient(cvcapp, gradient,*vol);
/*
  float *coefficients;
  int count=0;
  int dim[3];

  dim[0]=coeff->XDim(); dim[1]=coeff->YDim(); dim[2]=coeff->ZDim();

  coefficients = new float[ coeff->XDim() * coeff->YDim() * coeff->ZDim() ];

  int  index, index1;
  for(int p=0;p<coeff->XDim();p++)
    for(int q=0;q<coeff->YDim();q++)
      for(int r=0;r<coeff->ZDim();r++)
  {

    index1 = (r*dim[1]+q)*dim[0]+p;
    coefficients[index1]= (*coeff)(p,q,r);
  }*/
  
  int numverts=mesh->vertices.size();

  int numDirtyVerts=0;
  int fixed=0;
  tolerance/=2; // To be fair to both parties.
  isovalue=isovalue*(1+tolerance);



  for(int i=0;i<numverts;i++)
  {
    double x=mesh->vertices[i].x;
    double y=mesh->vertices[i].y;
    double z=mesh->vertices[i].z;
    
    double value=vol->interpolate(x,y,z);

    if(value<=isovalue*(1-EPSILON)) // We have a point inside!
    {
      double newvalue;
      int itercount=0;
     // cout<<"Before:  "<<value<<" After: ";
     do
     {
        double d, dist=abs(value-isovalue);
        double abs_grad_hls;
        
        float d_hls[10]; // The partial derivatives in x,y and z directions.
        float p[3];
        float dxyz[3];
        float minExt[3];
        
        d_hls[0]=value;
//         p[0]=x; dxyz[2]=coeff->XSpan(); minExt[0]=coeff->XMin();
//         p[1]=y; dxyz[1]=coeff->YSpan(); minExt[1]=coeff->YMin();
//         p[2]=z; dxyz[0]=coeff->ZSpan(); minExt[2]=coeff->ZMin();
            
        // Compute the absolute value of the gradient at that point.
        d_hls[0]=gradient[0].interpolate(x,y,z);
        d_hls[1]=gradient[1].interpolate(x,y,z);
        d_hls[2]=gradient[2].interpolate(x,y,z);

//         ReSamplingCubicSpline1OrderPartials(coefficients,coeff->XDim(),coeff->YDim(),coeff->ZDim(),dxyz, minExt, d_hls, p);
        
        abs_grad_hls=d_hls[0]*d_hls[0]+d_hls[1]*d_hls[1]+d_hls[2]*d_hls[2];
//         abs_grad_hls=d_hls[3]*d_hls[3]+d_hls[1]*d_hls[1]+d_hls[2]*d_hls[2];

        if(abs_grad_hls<=0.00000001) abs_grad_hls=1.0;

        d=dist/abs_grad_hls; // Get distance

        // Update the distance with the tolerance.
        //d*=(1+tolerance/2);

        //Move the point.
        mesh->vertices[i].x+=  d * d_hls[0] ;
        mesh->vertices[i].y+=  d * d_hls[1] ;
        mesh->vertices[i].z+=  d * d_hls[2] ;

//         mesh->vertices[i].x+=  d * d_hls[3] ;
//         mesh->vertices[i].y+=  d * d_hls[2] ;
//         mesh->vertices[i].z+=  d * d_hls[1] ;

        
        newvalue=vol->interpolate(mesh->vertices[i].x,mesh->vertices[i].y,mesh->vertices[i].z);
        itercount++;

        if(newvalue>=isovalue*(1-EPSILON)&&newvalue<=isovalue*(1+EPSILON))
            fixed++;
        else
          value=newvalue;

        
        
     }while(newvalue<isovalue*(1-EPSILON)&&itercount<50);
   // cout<<itercount<<" "<<newvalue<<endl;

      numDirtyVerts++;
    }
  }

  if(numDirtyVerts>0) dirty=true;
  else dirty=false;
  cout<<"NumDirtyVerts = "<<numDirtyVerts<<", fixed = "<<fixed<<endl;
//   delete[] coefficients;
  return numDirtyVerts-fixed;
}

bool testBBoxIntersect(const BoundingBox &bb1,const BoundingBox &bb2)
// bool testBBoxIntersect(const BoundingBox &b1,const VolMagick::Volume &v,const float isovalue,const float tolerance)
{
//   bool x_b1_test = (b1.minX <= b2.minX && b2.minX <= b1.maxX) || (b1.minX <= b2.maxX && b2.maxX <= b1.maxX),
//         y_b1_test = (b1.minY <= b2.minY && b2.minY <= b1.maxY) || (b1.minY <= b2.maxY && b2.maxY <= b1.maxY),
//         z_b1_test = (b1.minZ <= b2.minZ && b2.minZ <= b1.maxZ) || (b1.minZ <= b2.maxZ && b2.maxZ <= b1.maxZ);
//   bool x_b2_test = (b2.minX <= b1.minX && b1.minX <= b2.maxX) || (b2.minX <= b1.maxX && b1.maxX <= b2.maxX),
//         y_b2_test = (b2.minY <= b1.minY && b1.minY <= b2.maxY) || (b2.minY <= b1.maxY && b1.maxY <= b2.maxY),
//         z_b2_test = (b2.minZ <= b1.minZ && b1.minZ <= b2.maxZ) || (b2.minZ <= b1.maxZ && b1.maxZ <= b2.maxZ);
//   
//    if((x_b1_test && y_b1_test && z_b1_test) || (x_b2_test && y_b2_test && z_b2_test))
//     return true;

//   if((b1.maxX >= b2.minX && b1.maxX <= b2.maxX) && (b1.maxY >= b2.minY && b1.maxY <= b2.maxY) &&(b1.maxZ >= b2.minZ && b1.maxZ <= b2.maxZ)) return true;
//   if((b1.maxX >= b2.minX && b1.maxX <= b2.maxX) && (b1.maxY >= b2.minY && b1.maxY <= b2.maxY) &&(b1.minZ >= b2.minZ && b1.minZ <= b2.maxZ)) return true;
//   if((b1.maxX >= b2.minX && b1.maxX <= b2.maxX) && (b1.minY >= b2.minY && b1.minY <= b2.maxY) &&(b1.maxZ >= b2.minZ && b1.maxZ <= b2.maxZ)) return true;
//   if((b1.minX >= b2.minX && b1.minX <= b2.maxX) && (b1.maxY >= b2.minY && b1.maxY <= b2.maxY) &&(b1.maxZ >= b2.minZ && b1.maxZ <= b2.maxZ)) return true;
//   if((b1.maxX >= b2.minX && b1.maxX <= b2.maxX) && (b1.minY >= b2.minY && b1.minY <= b2.maxY) &&(b1.minZ >= b2.minZ && b1.minZ <= b2.maxZ)) return true;
//   if((b1.minX >= b2.minX && b1.minX <= b2.maxX) && (b1.maxY >= b2.minY && b1.maxY <= b2.maxY) &&(b1.minZ >= b2.minZ && b1.minZ <= b2.maxZ)) return true;
//   if((b1.minX >= b2.minX && b1.minX <= b2.maxX) && (b1.minY >= b2.minY && b1.minY <= b2.maxY) &&(b1.maxZ >= b2.minZ && b1.maxZ <= b2.maxZ)) return true;
//   if((b1.minX >= b2.minX && b1.minX <= b2.maxX) && (b1.minY >= b2.minY && b1.minY <= b2.maxY) &&(b1.minZ >= b2.minZ && b1.minZ <= b2.maxZ)) return true;
//   
//   
//   if((b2.maxX >= b1.minX && b2.maxX <= b1.maxX) && (b2.maxY >= b1.minY && b2.maxY <= b1.maxY) &&(b2.maxZ >= b1.minZ && b2.maxZ <= b1.maxZ)) return true;
//   if((b2.maxX >= b1.minX && b2.maxX <= b1.maxX) && (b2.maxY >= b1.minY && b2.maxY <= b1.maxY) &&(b2.minZ >= b1.minZ && b2.minZ <= b1.maxZ)) return true;
//   if((b2.maxX >= b1.minX && b2.maxX <= b1.maxX) && (b2.minY >= b1.minY && b2.minY <= b1.maxY) &&(b2.maxZ >= b1.minZ && b2.maxZ <= b1.maxZ)) return true;
//   if((b2.minX >= b1.minX && b2.minX <= b1.maxX) && (b2.maxY >= b1.minY && b2.maxY <= b1.maxY) &&(b2.maxZ >= b1.minZ && b2.maxZ <= b1.maxZ)) return true;
//   if((b2.maxX >= b1.minX && b2.maxX <= b1.maxX) && (b2.minY >= b1.minY && b2.minY <= b1.maxY) &&(b2.minZ >= b1.minZ && b2.minZ <= b1.maxZ)) return true;
//   if((b2.minX >= b1.minX && b2.minX <= b1.maxX) && (b2.maxY >= b1.minY && b2.maxY <= b1.maxY) &&(b2.minZ >= b1.minZ && b2.minZ <= b1.maxZ)) return true;
//   if((b2.minX >= b1.minX && b2.minX <= b1.maxX) && (b2.minY >= b1.minY && b2.minY <= b1.maxY) &&(b2.maxZ >= b1.minZ && b2.maxZ <= b1.maxZ)) return true;
//   if((b2.minX >= b1.minX && b2.minX <= b1.maxX) && (b2.minY >= b1.minY && b2.minY <= b1.maxY) &&(b2.minZ >= b1.minZ && b2.minZ <= b1.maxZ)) return true;
//     float ptsA[8][3];
//     
//     ptsA[0][0]=b1.minX; ptsA[0][1]=b1.minY; ptsA[0][2]=b1.minZ;
//     ptsA[1][0]=b1.maxX; ptsA[1][1]=b1.minY; ptsA[1][2]=b1.minZ;
//     ptsA[2][0]=b1.maxX; ptsA[2][1]=b1.maxY; ptsA[2][2]=b1.minZ;
//     ptsA[3][0]=b1.minX; ptsA[3][1]=b1.maxY; ptsA[3][2]=b1.minZ;
//     ptsA[4][0]=b1.minX; ptsA[4][1]=b1.minY; ptsA[4][2]=b1.maxZ;
//     ptsA[5][0]=b1.maxX; ptsA[5][1]=b1.minY; ptsA[5][2]=b1.maxZ;
//     ptsA[6][0]=b1.maxX; ptsA[6][1]=b1.maxY; ptsA[6][2]=b1.maxZ;
//     ptsA[7][0]=b1.minX; ptsA[7][1]=b1.maxY; ptsA[7][2]=b1.maxZ;
    
//     float ptsB[8][3];
    
//     ptsB[0][0]=b2.minX; ptsB[0][1]=b2.minY; ptsB[0][2]=b2.minZ;
//     ptsB[1][0]=b1.maxX; ptsB[1][1]=b1.minY; ptsB[1][2]=b1.minZ;
//     ptsB[2][0]=b1.maxX; ptsB[2][1]=b1.maxY; ptsB[2][2]=b1.minZ;
//     ptsB[3][0]=b1.minX; ptsB[3][1]=b1.maxY; ptsB[3][2]=b1.minZ;
//     ptsB[4][0]=b1.minX; ptsB[4][1]=b1.minY; ptsB[4][2]=b1.maxZ;
//     ptsB[5][0]=b1.maxX; ptsB[5][1]=b1.minY; ptsB[5][2]=b1.maxZ;
//     ptsB[6][0]=b2.maxX; ptsB[6][1]=b2.maxY; ptsB[6][2]=b2.maxZ;
//     ptsB[7][0]=b1.minX; ptsB[7][1]=b1.maxY; ptsB[7][2]=b1.maxZ;
//     
//     for(int i=0;i<8;i++)
//     {
//       float value=v(pts[i][0],pts[i][1],pts[i][2]);
//       
//       if(value<isovalue*(1+tolerance/2)*(1.0-EPSILON))
//         return true;
//     }
//   float k[6],pts[6][3];
//   bool test[3];
//   for(int i=0;i<3;i++)
//   {
//     k[i]=(ptsB[6][i]-ptsA[0][i])/(ptsB[6][i]-ptsB[0][i]);
//     if(k[i]>=0 && k[i]<=1) return true;
//     {
//       for(int j=0; j<3; j++)
//       {
//         pts[i][j]=k[i]*ptsB[0][j]+(1.0-k[i])*ptsB[6][j];
//         test[j]=pts[i][j]<=ptsA[6][j] &&pts[i][j]>=ptsA[0][j];
//       }
//      if (test[0] || test[1] || test[2]) return true;
//     }
//   }
//   for(int i=3;i<6;i++)
//   {
//       k[i]=(ptsB[6][i-3]-ptsA[6][i-3])/(ptsB[6][i-3]-ptsB[0][i-3]);
//       if (k[i] >= 0 && k[i]<=1) return true;
//       {
//         for(int j=0; j<3; j++)
//         {
//             pts[i][j]=k[i]*ptsB[0][j]+(1.0-k[i])*ptsB[6][j];
//             test[j]=  pts[i][j]<=ptsA[6][j] &&pts[i][j]>=ptsA[0][j];
//         }
//       }
//       if (test[0] || test[1] || test[2]) return true;
//   }
  
 /* 
  for(int i=0;i<6;i++) cout<<k[i]<<" ";
  cout<<endl;
  
  for(int i=0;i<6;i++)
    if(0<=k[i]&&k[i]<=1) return true;*/
  if(bb1.maxX<bb2.minX || bb1.minX > bb2. maxX) return false;
  if(bb1.maxY<bb2.minY || bb1.minY > bb2. maxY) return false;
  if(bb1.maxZ<bb2.minZ || bb1.minZ > bb2. maxZ) return false;
  
  return true;
}

void correctMeshes(vector<string> &vols,vector<string> &procvols,vector<Mesh> &meshes,vector<float> &isovalues,float tolerance)
{
  int dirtyVertCount=0,iterations=0;
 
 // int **seenArray;
 // seenArray=new int *[vols.size()];

 // for(int i=0;i<meshes.size();i++)
 //   seenArray[i]=new int [meshes.size()];

  // Perform simple normal correction
  cout<<"Performing simple normal correction."<<endl;
  
  for(unsigned int i=0;i<procvols.size();i++)
  {
    VolMagick::Volume v;
    VolMagick::readVolumeFile(cvcapp, v,vols[i]);
    projectMeshToIsosurface(&v,&meshes[i],isovalues[i]); // Really tight.
    cout<<"Peforming geometric flow for #"<<i<<endl;
    geometricFlow(meshes[i]);
    meshes[i].computeBoundingBox();
  }

  do
  {
    dirtyVertCount=0;
    for(unsigned int i=0;i<vols.size();i++)
     for(unsigned int j=0;j<meshes.size();j++)
        if(j!=i)
        { 
//           if(i==5&&j==25)
//           {
//             cout<<"Mesh #5: BB- ("<<meshes[i].bb.minX<<","<<meshes[i].bb.minY<<","<<meshes[i].bb.minZ<<")-("<<meshes[i].bb.maxX<<","<<meshes[i].bb.maxY<<","<<meshes[i].bb.maxZ<<")"<<endl; 
//             cout<<"Mesh #5: BB- ("<<meshes[j].bb.minX<<","<<meshes[j].bb.minY<<","<<meshes[j].bb.minZ<<")-("<<meshes[j].bb.maxX<<","<<meshes[j].bb.maxY<<","<<meshes[j].bb.maxZ<<")"<<endl;
//           }
          if(testBBoxIntersect(meshes[j].bb,meshes[i].bb))
          {
            VolMagick::Volume v;
            VolMagick::readVolumeFile(cvcapp, v,vols[i]);
            int numDirtyVerts;
            bool dirt;
            cout<<"Processing mesh #"<<j<<" with volume #"<<i<<"."<<endl;
            
            numDirtyVerts=correctMeshByProjection(&v,&meshes[j],isovalues[i],tolerance,dirt);
            dirtyVertCount+=numDirtyVerts;
    
            if(dirt)
            {
                for(int g=0;g<1;g++)
                {
                    cout<<"Peforming geometric flow iter# "<<g<<" for #"<<j<<endl;
                    geometricFlow(meshes[j]);
                }
                 meshes[j].computeBoundingBox();
            }
        }
        else
        {
            cout<<"Bounding box non-intersecting for meshes #"<<j<<" and #"<<i<<endl;
        }
        }

    iterations++;

    cout<<"Iteration #"<<iterations<<"; "<<dirtyVertCount<<" dirty vertices left."<<endl;
    
  }while((dirtyVertCount>0)&&iterations<3);

 // for(int i=0;i<vols.size();i++)
 //   delete[] seenArray[i];
 // delete[] seenArray;
      
}
