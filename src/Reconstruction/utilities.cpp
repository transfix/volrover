/* 
All right are reserved by The National Key Lab of Scientific and Engineering Computing,
Chinese Academy of Sciences.

        Author: Ming Li <ming@lsec.cc.ac.cn>
        Advisor: Guoliang Xu <xuguo@lsec.cc.ac.cn>

This file is part of CIMOR. CIMOR stands for:
        Computational Inverse Methods of Reconstruction

*/


#include <stdio.h>
#include <cmath>
#include <Reconstruction/utilities.h>
#include <Reconstruction/Reconstruction.h>
#include <cstdio>



float *VertArray;
int   *FaceArray, numbpts, numbtris;





void ObtainRotationMatFromViews(float* Rmat, struct Views* view)
{
 int i, j;
 Views* v;
 float rotmat[9];

 for ( v=view, i=0; v; v = v->next, i++ ) 
    {
     ObtainRotationMatrixAroundOrigin(rotmat,v);

     for ( j = 0; j < 9; j++ ) Rmat[i*9 + j] = rotmat[j];

    }
}


void InverseRotationMatrix(float *Rmat,int nv)
{
  float rotmat[9], rotmatinv[9];
  int i,j;
  for( i = 0; i < nv; i++ )
	{
    for ( j = 0; j < 9; j++ ) 
	  rotmat[j] = Rmat[i*9+j];
  
    Matrix3Transpose(rotmat, rotmatinv);

    for ( j = 0; j < 9; j++ ) 
      Rmat[i*9+j] = rotmatinv[j];             

    }



}


float  angle_set_negPI_to_PI(float angle)
{

        if ( !std::isfinite(angle) ) return(0);
        angle -= TWOPI*((int) (angle/TWOPI));
        while ( angle <= -M_PI ) angle += TWOPI;
        while ( angle >   M_PI ) angle -= TWOPI;

        return(angle);
}


float    bcos(float x)
{
        x = angle_set_negPI_to_PI (x);

        if( fabs(x - M_PI_2) < TRIGPRECISION || fabs(x + M_PI_2) < TRIGPRECISION )      //x=+90 or -90 degrees
                return (0.0);

        return (cos(x));
}

float    bsin(float x)
{
        x = angle_set_negPI_to_PI (x);

        if( fabs(x - M_PI) < TRIGPRECISION || fabs(x + M_PI) < TRIGPRECISION )  //x=180 or -180 degrees
                return (0.0);

        return(sin(x));
}



/*****************************************************************************
 Descriptions:
      Obtain a rotation matrix which rotate a volume data around coordinate 
      origin(also is the volume origin) to make the view vector be the z-axis.  

Arguments:
      rotmat     : The output rotation matrix.
      [v1,v2,v3] : The view unit vector(line of sight). 
      
 *****************************************************************************/
void ObtainRotationMatrixAroundOrigin(float* rotmat,struct Views* v)
{
 float p,norm;
 norm = v->x*v->x + v->y*v->y + v->z*v->z;
 if(fabs(norm-1) > 0.000001) 
    {
     printf("\n view vector is not a unit vector. ");

    }
   
int      j;
float    a = v->a/2.0;
float    ca = bcos(a);
float    sa = bsin(a);
float    z1 = v->z + 1;
float    f1, f2;
float    s, x,y,z;


if( v->z == 1.0 ) 
{
 s= ca;
 x = 0;
 y = 0;
 z = sa;
//printf("\n s = %f x = %f ,y = %f ,z = %f ", s ,x,y,z);
//getchar();

}



if ( v->z < 1.0 - SMALLFLOAT ) {
  if ( v->z < SMALLFLOAT - 1.0 ) {   // z= -1.
     s = 0;
     x = sa;
     y = ca;
     z =  0;
    
     }

 else {                 
                                            // -1 < z < 1
     f1 = sqrt(z1/2.0);
     f2 = sqrt(1.0/(2*z1));
      s = f1*ca;
      x = f2*(v->x*sa - v->y*ca);
      y = f2*(v->x*ca + v->y*sa);
      z = f1*sa;
     }

}
 rotmat[0] = s*s + x*x - y*y - z*z;
 rotmat[1] = 2*x*y - 2*s*z;
 rotmat[2] = 2*x*z + 2*s*y;
 rotmat[3] = 2*x*y + 2*s*z;
 rotmat[4] = s*s - x*x + y*y - z*z;
 rotmat[5] = 2*y*z - 2*s*x;
 rotmat[6] = 2*x*z - 2*s*y;
 rotmat[7] = 2*y*z + 2*s*x;
 rotmat[8] = s*s - x*x - y*y + z*z;
 
 
 //RotateMatrix_z(v->x, v->y, v->z, rotmat);
 

//my rotation matrix.
/* 
if( v->x == 0.0 && v->y ==0.0 && v->z > 0.0) 
    {
     printf("\n view vector is z-axis.don't need rotation.");
     rotmat[1] = rotmat[2] = rotmat[3] = rotmat[5] = rotmat[6] = rotmat[7] = 0.0; 
     rotmat[0] = rotmat[4] = rotmat[8] = 1.0;
    }

else
   {
    p = sqrt(1-v->z*v->z);

    rotmat[0] = -v->x*v->z/p;             // v->y/p; 
    rotmat[1] = -v->y*v->z/p;             // -v->x/p; 
    rotmat[2] = p;                       // 0.0;
    rotmat[3] = v->y/p;                  // -v->x*v->z/p; 
    rotmat[4] = -v->x/p;                // -v->y*v->z/p;
    rotmat[5] =  0.0;                    // p; 
    rotmat[6] =  v->x;
    rotmat[7] =  v->y;
    rotmat[8] =  v->z;
   }

*/



 /*
 float oldx[3]={1.0,0.0,0.0}, rx[3],newx[3], angle, rotmatinv[9];
 int j;
 for ( j = 0; j < 9; j++ )
   printf("\n rotmat = %f ", rotmat[j]);

 p = sqrt(1-v->z*v->z);
 
 // if(v->x <0 || (v->x == 0 && v->y >0) ||(v->x ==0 && v->y ==0 && v->z > 0) ) {
 rx[0] = -v->x*v->z/p;
 rx[1] = -v->y*v->z/p;
 rx[2] = p;
 //}

 

 Matrix3Transpose(rotmat, rotmatinv); 
 float nv[3]={v->x, v->y,v->z}, ry[3] = {v->y/p, -v->x/p, 0.0};

 printf("\n rx = %f %f %f ", rx[0], rx[1], rx[2]);
 printf("\n ry = %f %f %f ",v->y/p, -v->x/p, 0.0);
 
  MatrixMultiply(rotmatinv,3,3,rx,3,1,newx);

 printf("\n newx = %f %f %f ", newx[0], newx[1], newx[2]);

 // MatrixMultiply(rotmatinv,3,3,nv,3,1,newx);
 // printf("\n vz = %f %f %f ", newx[0], newx[1], newx[2]);

 //MatrixMultiply(rotmatinv,3,3,ry,3,1,newx);
 //printf("\n vy = %f %f %f ", newx[0], newx[1], newx[2]);


 angle =  Angle_Of_Two_Vectors(oldx, newx);
 printf("\n angle = %f ", angle);
 

 */


}


/*
int main(int argc,char **argv)
{
float v1,v2,v3, norm, v[3];
float rotmat[9];
v1 = 0;
v2 = 0.00001;
v3 = 1;

norm = sqrt(v1*v1 + v2*v2 + v3*v3);
v1 = v1/norm;
v2 = v2/norm;
v3 = v3/norm;
v[0] = v1; v[1] = v2; v[2] = v3;


ObtainRotationAroundOrigin(rotmat, v1,v2,v3);
float t[3];
t[0] = t[1] = t[2] = 0;

for(int i = 0; i < 3; i++)
   for(int j = 0; j < 3; j++)
      {
       t[i]=t[i] + rotmat[i*3+j]*v[j];
      }

printf("result vector = %f %f %f ", t[0],t[1],t[2]);


}

*/
 
void Matrix3Transpose(float* matrix, float* trmatrix) 
{
int   i, j;
for ( i=0; i<3; i++ )
    for ( j=0; j<3; j++ )
        trmatrix[3*i+j] = matrix[3*j+i];


}

void MatrixTranspose(float *A,int m,int n)
{
 int i,j;
 float *At=(float *)malloc(m*n*sizeof(float));
 for(i=0;i<m;i++)
         for(j=0;j<n;j++)
                 At[j*m+i]=A[i*n+j];
 for(i=0;i<m*n;i++)
         A[i]=At[i];
 free(At);

}




void  MatrixMultiply(float* A1,int m1,int n1,float* A2,int m2,int n2,float* Affi)
{
  int i,j,k, in1, in2;
 float sum;
 if(n1!=m2) {printf("A cannot multiply B.");return;}
 for(i=0;i<m1;i++)
   {
	 in1 = i*n1;
	 in2 = i*n2;

   for(j=0;j<n2;j++)
	 {
	 sum=0.0;

	 for(k=0;k<n1;k++)
	   sum= sum+A1[in1+k]*A2[k*n2+j];

	 Affi[in2+j]=sum;

	 }
   }
}






void  ReadRawFile(char *filename)
{
 int i,j,k,u,v,w;
 int npts,nbtris;
 float  x=0.0,y=0.0,z=0.0;
 FILE *fp;

 fp=fopen(filename,"r");
 if ( (fp==NULL) || (fp==0) )
         {printf("\nCould not open the file surface.raw for writing!\n");return;}


    fscanf(fp, "%d %d", &numbpts, &numbtris);
    VertArray=(float *)malloc(3*numbpts*sizeof(float));
    FaceArray=(int *)malloc(3*numbtris*sizeof(int));
    for(i=0;i<3*numbpts;i++) VertArray[i]=0.0;
    for(i=0;i<3*numbtris;i++) FaceArray[i]=0;

    for (i=0; i<numbpts; i++)
    {
        fscanf(fp, "%f %f %f", &x, &y, &z);
        VertArray[3*i+0] = x;
        VertArray[3*i+1] = y;
        VertArray[3*i+2] = z;
    }
    for (i=0; i<numbtris; i++)
    {
        fscanf(fp, "%d %d %d", &u, &v, &w);
        FaceArray[3*i+0] = u;
        FaceArray[3*i+1] = v;
        FaceArray[3*i+2] = w;
    }
 fclose(fp);
  

}

struct Views* ObtainViewVectorFromTriangularSphere(char *filename)
{
 
float       norm,vec1[3],vec2[3];
//int      *faceArray,numbpts,numbtris;
int i;
Views* view, *head;

view = (Views *)malloc(sizeof(Views));
head = view;

ReadRawFile(filename);

for( i = 0; i < numbpts; i++)
   {
    vec1[0] =  VertArray[3*i+0];
    vec1[1] =  VertArray[3*i+1];
    vec1[2] =  VertArray[3*i+2];

    norm = sqrt(vec1[0]*vec1[0]+vec1[1]*vec1[1]+vec1[2]*vec1[2]);
    
    vec1[0] = vec1[0]/norm;
    vec1[1] = vec1[1]/norm;
    vec1[2] = vec1[2]/norm;

    view->x = vec1[0];
    view->y = vec1[1];
    view->z = vec1[2];

    view->a = 0.0;
    if( i < numbpts-1) 
      {
       view->next =  (Views *)malloc(sizeof(Views));
       view = view->next;
      }
    else view->next = NULL;

 


  }


/*  //Find the parallel vector.
for( i = 0; i < numbpts-1; i++ )
   {
    vec1[0] =  VertArray[3*i+0];
    vec1[1] =  VertArray[3*i+1];
    vec1[2] =  VertArray[3*i+2];
    for( j = i+1; j < numbpts; j++ )
       {
        vec2[0] =  VertArray[3*j+0];
        vec2[1] =  VertArray[3*j+1];
        vec2[2] =  VertArray[3*j+2];  
        if(fabs(InnerProduct(vec1,vec2))<0.001)




       } 


   } 

*/

view->next = NULL;
return head;


}


/*************************************************************************/
/* InnerProduct                                                          */
/*************************************************************************/
float InnerProduct(float u[],float v[])
{
 float w;
 w = u[0]*v[0] + u[1]*v[1] + u[2]*v[2];
 return(w);
}

long double InnerProductl(long double *u,long double *v,int n)
{
long double w;
int i;

 w = 0.0L;
 for ( i = 0; i < n; i++ )
   w = w+u[i]*v[i];

return(w);
}




/*************************************************************************/
/* GaussInverse                                                          */
/* - Inverse the matrix A by Gauss elimination methods                   */
/*************************************************************************/
void gaussinverse (long double* a, int n, long double eps,int* message)
//int             n, *message;
///long double         *a, eps;
///*   n        -integer, the dimension of the two dimensional matrix a[n][n].
//a        -pointer, point to a matrix to be inversed and also the
//computing result
//eps      -controllor, that control the computation by testing the
//pivot element
//message, -pointer, point to an integer, if pivot element less than eps
//message = 0,
//otherwise message = 1.
//*/

{
        long double        max;
        int           k,ik,jk,i,j,*z;

        z =  (int *)malloc(2*n*sizeof(int));
        *message = 1;
        for (k=0; k<n; k++) {
                max = 0.0L;
                for (i = k;i<n; i++)
                        for (j = k;j<n; j++)
                                if (fabsl(*(a+i*n+j)) > max) {
                                        ik = i;
                                        jk = j;
                                        max = fabsl(*(a+i*n+j));
                                };
                                if (max < eps     || max == 0.0L) {
                                        *message = 0;
                                        printf("The matrix in gaussinverse is singular, %Lf\n",max);
                                        return;
                                };
                                max = 1.0L/ *(a+ik*n+jk);
                                *(a+ik*n+jk) = 1.0L;
                                z[2*k] = ik;
                                z[2*k+1] = jk;
                                Exchangerowcolumn(a,n,k,ik,jk);
                                for (j = 0;j<n; j++)
                                        *(a+k*n+j)  =  *(a+k*n+j) * max;
                                for (i = 0;i<n; i++)
                                        if (i != k) {
                                                max = *(a+i*n+k);
                                                *(a+i*n+k) = 0.0L;
                                                for (j = 0;j<n; j++)
                                                        *(a+i*n+j) = *(a+i*n+j)- max * *(a+k*n+j);
                                        }
        };
        for (k = n-2; k > -1; k--){
                ik = z[2*k+1];
                jk = z[2*k];
                Exchangerowcolumn(a,n,k,ik,jk);
        }
        free(z);
}


void Exchangerowcolumn(long double* a, int n,int k, int ik, int jk)
//int             n,k,ik,jk;
////long double          *a;
///*  exchange the (k,ik) rows
//(k,ik) columns.
//*/
{
 long double         b;
 int            j;
 /*  exchange the (k,ik) rows */
 if (ik != k){
	 for  (j=0; j<n; j++) {
		 b = *(a+ik*n+j);
		 *(a+ik*n+j) =  *(a+k*n+j);
		 *(a+k*n+j) = b;
		 
	 }
	 
 }       /*  exchange the (k,ik) column  */
 if (jk != k){
	 for  (j=0; j<n; j++) {
		 b = *(a+j*n+jk);
		 *(a+j*n+jk) =  *(a+j*n+k);
		 *(a+j*n+k) = b;
		  
	 }
 }

}




Views*  quaternion_from_view_vector(Views* view)
{
  //  view_normalize(&view);

  //  Quaternion      q(1,0,0,0);                             // z = 1

		float  s,x,y,z;

        float  factor;

		s = 1.0; x = y = z = 0.0;   // z = 1;
     

        if ( view->z < 1.0 - SMALLFLOAT ) {
                if ( view->z < SMALLFLOAT - 1.0 ) {     // z = -1
                        s = 0;
                        y = 1;
                } else {                                                        // -1 < z < 1
                        s      = sqrt((1.0 + view->z)/2.0);
                        factor = sqrt(0.5/(1 + view->z));
                        x      = -view->y*factor;
                        y      = view->x*factor;
                }
        }

		//        q.check();

		// if ( verbose & VERB_DEBUG )
                printf("DEBUG quaternion_from_view_vector: %f %f %f %f\n",
                                s, x, y, z);
				float rotmat[9];
				int i;
 rotmat[0] = s*s + x*x - y*y - z*z;
 rotmat[1] = 2*x*y - 2*s*z;
 rotmat[2] = 2*x*z + 2*s*y;
 rotmat[3] = 2*x*y + 2*s*z;
 rotmat[4] = s*s - x*x + y*y - z*z;
 rotmat[5] = 2*y*z - 2*s*x;
 rotmat[6] = 2*x*z - 2*s*y;
 rotmat[7] = 2*y*z + 2*s*x;
 rotmat[8] = s*s - x*x - y*y + z*z;	
 for (i = 0; i < 9; i++ ) printf("\n rotmat = %f ", rotmat[i]);
	
				//   return(q);

				
        float           a = atan2(z, s);
        float          ca = bcos(a);
        float          sa = bsin(a);
        float          q2 = s*s + z*z;
        float          f2 = 2*sqrt(q2);


 if ( q2 < SMALLFLOAT ) a = atan2(x, y);

 struct   Views*  v;
 v = (Views*)malloc(sizeof(Views));
 v->next  = NULL;
 v->x = f2*(y*ca + x*sa);
 v->y = f2*(y*sa - x*ca);
 v->z = 2*q2 - 1;
 v->a = (float)2*a;
        //view_normalize(&view);

		// if ( verbose & VERB_DEBUG )
                printf("DEBUG view_from_quaternion: %f %f %f %f\n",
                                v->x, v->y, v->z, v->a);

				


		

        return(v);



}


/*---------------------------------------------------------------------------
Author : Prof. Xu.
Angle_Of_Three_Points -- compute the angle of three points p1,p2,p3. That
                         is the angle between p1-p2 and p3-p2
----------------------------------------------------------------------------*/
float  Angle_Of_Two_Vectors(float* p12, float* p32)
{
float  result, length1, length3;

//p12[0] = p1[0] - p2[0];
//p12[1] = p1[1] - p2[1];
//p12[2] = p1[2] - p2[2];

//p32[0] = p3[0] - p2[0];
//p32[1] = p3[1] - p2[1];
//p32[2] = p3[2] - p2[2];

length1 = InnerProduct(p12,p12);
if (length1 < 0.0000001) return(0.0);
length3 = InnerProduct(p32,p32);
if (length3 < 0.0000001) return(0.0);

length1 = sqrt(length1);
length3 = sqrt(length3);


result = InnerProduct(p12,p32)/(length1*length3);

if (result > 1.0) result = 1.0;
if (result < -1.0) result = -1.0;

result = acos(result);

return(result);
}



// Author : Prof. Xu.
/************************************************************************/
void RotateMatrix_z(float nx, float ny, float nz, float *matrix)
//float  nx,ny,nz;                  /* the  normal  at the given point    */
//float  matrix[3][3];              /* the rotating matrix                */
{
float  c1,c2,s1,s2, normal,normalz;

normal = sqrt(nx*nx + ny*ny + nz*nz);
normalz = sqrt(nx*nx + ny*ny);
c1 = nz/normal;       c2 = -ny/normalz;
s1 = normalz/normal;  s2 = nx/normalz;

if (normalz < 0.001) {
   matrix[0] = 1.0;    matrix[1] = 0.0;    matrix[2] = 0.0;
   matrix[3] = 0.0;    matrix[4] = 1.0;    matrix[5] = 0.0;
   matrix[6] = 0.0;    matrix[7] = 0.0;    matrix[8] = 1.0;
}
if (normalz >= 0.001) {
   matrix[0] = c2;     matrix[1] = -c1*s2; matrix[2] = s1*s2;
   matrix[3] = s2;     matrix[4] = c1*c2;  matrix[5] = -s1*c2;
   matrix[6] = 0.0;    matrix[7] = s1;     matrix[8] = c1;
}
}


float TrilinearInterpolation8(float xd, float yd, float zd, float x1, float x2, float x3, float x4, float x5, float x6, float x7, float x8)
{
  float C00, C10, C01, C11;
  float C0, C1, C;

  C00 = (1.0-xd) * x1 + xd * x2;
  C10 = (1.0-xd) * x3 + xd * x4;
  C01 = (1.0-xd) * x5 + xd * x6;
  C11 = (1.0-xd) * x7 + xd * x8;

  C0 = (1.0-yd) * C00 + yd * C10;
  C1 = (1.0-yd) * C01 + yd * C11;

  C = (1.0-zd) * C0 + zd * C1;


  return C;


}

bool  fft1D_shift(fftw_complex *Vec,int n)
{


  // if(fabs(n/2.0-n/2) >0.0 ) return false;
  
  float temp, temp1, tmp, tmp1;
  int i, n_2;

  n_2 = n/2;
  printf("\nn=%d ", n);getchar();
  if(n%2 == 0)
	{
  for ( i = 0; i <n_2; i++ )
	{
	  temp   = Vec[i][0];
	  temp1  = Vec[i][1]; 

	  Vec[i][0] = Vec[i+n_2][0];
	  Vec[i][1] = Vec[i+n_2][1];

	  Vec[i+n_2][0] = temp;
	  Vec[i+n_2][1] = temp1;
	
	}
	}
  else
	{ 
	  tmp  = Vec[n_2][0];
	  tmp1 = Vec[n_2][1];
 

  for ( i = 0; i <n_2; i++ )
	{
	  temp   = Vec[i][0];
	  temp1  = Vec[i][1]; 

	  Vec[i][0] = Vec[i+1+n_2][0];
	  Vec[i][1] = Vec[i+1+n_2][1];


	  Vec[i+n_2][0] = temp;
	  Vec[i+n_2][1] = temp1;
	
	}
  Vec[n-1][0] = tmp;
  Vec[n-1][1] = tmp1;

	}

  return true;

}

bool fft2D_shift(fftw_complex *Vec,int m, int n)
{

  int i,j, n_2, m_2, ii, jj, size;
  float temp1, temp2;
  fftw_complex *v;

  size = m * n;

  v = (fftw_complex *)malloc(m*n*sizeof(fftw_complex));
  for ( i = 0; i < m; i++)
	for ( j = 0; j < n; j++ )
	  {
		v[i*n+j][0] = Vec[i*n+j][0];
		v[i*n+j][1] = Vec[i*n+j][1];
	  }
  if(m != n ) return false;

  n_2 = n/2;
  m_2 = m/2;


 //if(fabs(n/2.0-n/2) >0.0 || fabs(m/2.0-m/2) ) return false;
  if(n%2 ==0 && m%2 == 0)
	{
	  for ( i = 0; i < m_2; i++ )
		for ( j = 0; j < n_2; j++ )
		  {
			temp1 =  Vec[i*n+j][0];
			temp2 =  Vec[i*n+j][1];
			
			Vec[i*n+j][0] = Vec[(i+m_2)*n+j+n_2][0];
			Vec[i*n+j][1] = Vec[(i+m_2)*n+j+n_2][1];
			
			Vec[(i+m_2)*n+j+n_2][0] = temp1;
			Vec[(i+m_2)*n+j+n_2][1] = temp2;
			
			
			temp1 =  Vec[i*n+j+n_2][0];
			temp2 =  Vec[i*n+j+n_2][1];
			
			Vec[i*n+j+n_2][0] = Vec[(i+m_2)*n+j][0];  
			Vec[i*n+j+n_2][1] = Vec[(i+m_2)*n+j][1];
			
			Vec[(i+m_2)*n+j][0] = temp1;
			Vec[(i+m_2)*n+j][1] = temp2;
			
		  }
	}
  else
	{
	  temp1 = Vec[n_2*n+m_2][0];
	  temp2 = Vec[n_2*n+m_2][1];

	  //Vec[n_2*n+m_2][0] = Vec[0][0];
	  //Vec[n_2*n+m_2][1] = Vec[0][1];

	  //Vec[0][0] = Vec[size-1][0];
	  //Vec[0][1] = Vec[size-1][1];

	  Vec[size-1][0] = temp1;
	  Vec[size-1][1] = temp2;


	  for ( i = 0; i < m; i++ )
		for ( j = 0; j < n; j++ )
		  {


			if(i < m_2) ii = m_2 + 1 + i;
			else ii = i -  m_2;

			if(j < n_2) jj =  n_2 + 1 + j;
			else jj = j - n_2;

			//			temp1 = v[i*n+j][0];
			//temp2 = v[i*n+j][1];


			Vec[i*n+j][0] = v[ii*n+jj][0];
			Vec[i*n+j][1] = v[ii*n+jj][1];


			//	Vec[(ii-1)*n+jj-1][0] = temp1;
			//Vec[(ii-1)*n+jj-1][1] = temp2;

		  }

	}
  fftw_free(v);
}


void EulerMatrice(float *Rmat,struct EulerAngles *Eulers)
{

 int i, j;
 EulerAngles* v;
 float rot, tilt, psi, rotmat[9];

 for ( v=Eulers, i=0; v; v = v->next, i++ ) 
    {
	  rot = v->rot;
	  tilt = v->tilt;
	  psi = v->psi;

	  euler2matrix(rot, tilt, psi,rotmat);

     for ( j = 0; j < 9; j++ ) Rmat[i*9 + j] = rotmat[j];

    }



}





void euler2matrix(float alpha, float beta, float gamma,
                         float A[9])
{
    float ca, sa, cb, sb, cg, sg;
    float cc, cs, sc, ss;

	//printf("\nalpha beta gamma=%f %f %f ", alpha, beta, gamma);
    alpha = alpha * M_PI/180.0; //DEG2RAD(alpha);
	beta  = beta * M_PI/180.0;  //DEG2RAD(beta);
	gamma = gamma * M_PI/180.0;  //DEG2RAD(gamma);

	//printf("\nrad alpha beta gamma=%f %f %f ", alpha, beta, gamma);

    ca = cos(alpha);
    cb = cos(beta);
    cg = cos(gamma);
    sa = sin(alpha);
    sb = sin(beta);
    sg = sin(gamma);
    cc = cb * ca;
    cs = cb * sa;
    sc = sb * ca;
    ss = sb * sa;

    A[0] =  cg * cc - sg * sa;
    A[1] =  cg * cs + sg * ca;
    A[2] = -cg * sb;
    A[3] = -sg * cc - cg * sa;
    A[4] = -sg * cs + cg * ca;
    A[5] =  sg * sb;
    A[6] =  sc;
    A[7] =  ss;
    A[8] =  cb;

	/*
	for ( int i = 0; i < 9; i++ )
	  printf("\nA=%f ", A[i]); 
	*/
}






void euler2view(float rot, float tilt, float psi, float view[3])
{
 

  view[0] = cos(rot) * sin(tilt);
  view[1] = sin(rot) * sin(tilt);
  view[2] = cos(tilt);
        

}


/*
Euler Angles ranges.
rot : 0-2PI;
tilt: 0-PI/2;
psi : 0-2PI;
 */

/******************************************************
p1: number of samples between 0 and 2PI for rot.
p2: number of samples between 0 and PI/2 for tilt.
p3: the same in-plane rotation angle for each rot and tilt.
***********************************************************/    
EulerAngles* phantomEulerAngles(float p1, float p2, float p3)
{
  int i, j, n1, n2, num;
  EulerAngles *Eulers, *head;

  n1  = (int)p1;
  n2  = (int)p2;
  num = n1 * n2;

  Eulers = (struct EulerAngles *)malloc(sizeof(EulerAngles));
  head = Eulers;

  printf("\nn1 n2 =%d %d", n1, n2);
  for ( i = 0; i < n1; i++ )
	{

	  for ( j = 0; j < n2; j++ )
		{
	
		  Eulers->rot = i * 1.0/n1 * 360 + 0.0;//PI2;
		  //Eulers->tilt = j * 1.0/n2 * 90 + 0;//PI/2.0;
		  Eulers->tilt = j * 1.0/n2 * 180 + 0.0;
		  Eulers->psi  = p3;
                  //printf("\neuler=%f %f %f ", Eulers->rot, Eulers->tilt, Eulers->psi);
		  // printf("\nrot tilt psi=%f %f %f ", Eulers->rot, Eulers->tilt, Eulers->psi);
		  if(i*n2 + j < num-1) 
			{
			  Eulers->next = (EulerAngles*)malloc(sizeof(EulerAngles));
			  Eulers = Eulers->next;
 
			}
		  else Eulers->next = NULL; 
	
		}

	  
	}


  return head;




}

EulerAngles* phantomEulerLimitedAngles(float p1, float p2, float p3)
{

  int i, j, n1, n2, num;
  EulerAngles *Eulers, *head;

  n1  = (int)p1;
  n2  = (int)p2;

  if(n1 == 0) n1 = 1;

  num = n1 * n2;


 
  Eulers = (struct EulerAngles *)malloc(sizeof(EulerAngles));
  head = Eulers;


  for ( i = 0; i < n1; i++ )
	{

	  for ( j = 0; j < n2; j++ )
		{
	
		  Eulers->rot  = i * 1.0/n1 * 360 + 0;//  rot = 0;
		  Eulers->tilt = -60 + j * 1.0/n2 * 120 + 0;// tilt \in[-60, 60];
		  Eulers->psi  = p3;
		  //printf("\nrot tilt psi=%f %f %f ", Eulers->rot, Eulers->tilt, Eulers->psi);
		  if(i*n2 + j < num-1) 
			{
			  Eulers->next = (EulerAngles*)malloc(sizeof(EulerAngles));
			  Eulers = Eulers->next;
 
			}
		  else Eulers->next = NULL; 
	
		}

	  
	}


  return head;



}





/*
Project a point to plane e1 X e2.
output:
point[0] = <point, e1>
point[1] = <point, e2>
point[2] = 0.0;

*/
void ProjectPoint2Plane(float point[3], float e1[3], float e2[3])
{

  float u, v;
  u = InnerProduct(point, e1);
  v = InnerProduct(point, e2);

  point[0] = u;
  point[1] = v;
  point[2] = 0.0;


/*
point[0] = u * e1[0] + v * e2[0];
point[1] = u * e1[1] + v * e2[1];
point[2] = u * e1[2] + v * e2[2];
*/
  

}

float MaxError(fftw_complex * ar1, fftw_complex * ar2, int size)
{

  int i,j;
  float r1, r2, s1, s2, max;
  max = 0.0;

  for ( i = 0; i < size; i++ )
	{
	  r1 = ar1[i][0];
	  r2 = ar2[i][0];

	  r1 = r1 - r2;

	  s1 = ar1[i][1];
	  s2 = ar2[i][1];

	  s1 = s1 - s2;

	  r1 = (fabs(r1)>fabs(s1))?fabs(r1):fabs(s1);
	  max = (max > r1)?max:r1;
	}
  return max;

}

float MaxError_L2(fftw_complex * ar1, fftw_complex * ar2, int size)
{
  int i,j;
  float r1, r2, s1, s2, max, temp;
  max = 0.0;

 for ( i = 0; i < size; i++ )
	{
	  r1 = ar1[i][0];
	  r2 = ar2[i][0];

	  r1 = r1 - r2;

	  s1 = ar1[i][1];
	  s2 = ar2[i][1];

	  s1 = s1 - s2;

	  temp = sqrt(r1*r1+ s1*s1);

      max = (temp > max)?temp:max;
 
	}

 return max;

}

/*----------------------------------------------------------------------------*/
void Volume_Projection(Oimage *Volume, float rotmat[9], int sample_num, float *prjimg)
{
  int     nx, ny, nz, nynz, partition_t, half, size, gM, ix, iy, iz;  
  int     i, j, k, numb_in, ii;
  float   e1[3], e2[3], d[3], a[3], R, R2, center[3];
  float   A, B, C, t1, t2, sum, delt, temp, values_t;
  float   xx, x, y, z, kot, weit, xd, yd, zd;
  long double Kot[4], Weit[4]; 
  float   r1, r2, r3, r4, r5, r6, r7, r8;
  

  Kot[0] = (-0.8611363115940526L + 1.0L)/2.0L;
  Kot[1] = (-0.3399810435848563L + 1.0L)/2.0L;
  Kot[2] =  (0.3399810435848563L + 1.0L)/2.0L;
  Kot[3] =  (0.8611363115940526L + 1.0L)/2.0L;
  
  Weit[0] = 0.3478548451374539L/2.0L;
  Weit[1] = 0.6521451548625461L/2.0L;
  Weit[2] = 0.6521451548625461L/2.0L;
  Weit[3] = 0.3478548451374539L/2.0L;

  nx = Volume->nx;
  ny = Volume->ny;
  nz = Volume->nz;

  nynz = ny * nz;

  //  for ( i = 0; i < nx * ny * nz; i++ )
  //printf("\nVolume=%f ", Volume->data[i]);
  //getchar();


  //printf("\nbegin compute volume projection.\n");
  
  partition_t = 4;
  partition_t = 1000;   //Xu
  gM          = 4 * partition_t; 

  half      = sample_num/2;  //sample_num=nx+1; half=nx/2;
  size      = sample_num*sample_num;
  center[0] = 0.0;
  center[1] = 0.0;
  center[2] = 0.0;

  R  = (Volume->nx-1)/2.0;
  R2 = R * R;
 
  for ( i = 0; i < size; i++ )
	prjimg[i] = 0.0;

  for ( i = 0; i < 3; i++ ) {
	  e1[i] = rotmat[i];
	  e2[i] = rotmat[3+i];
	  d[i]  = rotmat[6+i];

          if(fabs(e1[i])<SMALLFLOAT) e1[i] = 0.0;
          if(fabs(e2[i])<SMALLFLOAT) e2[i] = 0.0;


	  //printf("\ne1=%f e2=%f d=%f ", e1[i], e2[i], d[i]);
	}


  A = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];

// for ( i = -half; i < half; i++ )
//        for ( j = -half; j < half; j++ )
  for ( i = -half; i <= half; i++ )
	for ( j = -half; j <= half; j++ )
	  {
/*		a[0] = (i+0.5)*e1[0] + (j+0.5)*e2[0];      // centers are (0,0,0).
		a[1] = (i+0.5)*e1[1] + (j+0.5)*e2[1];
		a[2] = (i+0.5)*e1[2] + (j+0.5)*e2[2];
*/
 
                a[0] = i*e1[0] + j*e2[0];
                a[1] = i*e1[1] + j*e2[1];
                a[2] = i*e1[2] + j*e2[2];


		B = 2 * d[0] * (a[0] - center[0]) + 
		    2 * d[1] * (a[1] - center[1]) +
		    2 * d[2] * (a[2] - center[2]);
		C = (a[0] - center[0]) * (a[0] - center[0]) + 
		    (a[1] - center[1]) * (a[1] - center[1]) +
		    (a[2] - center[2]) * (a[2] - center[2]) - R2;
		

		if(B*B-4*A*C <= 0.0 ) continue;

		C  = sqrt(B*B-4*A*C);
						
		t1 = 0.5*(-B-C)/A;
		t2 = 0.5*(-B+C)/A;
		//printf("\nt1=%f t2=%f ", t1, t2);
		if(t1 > t2 ) {temp = t1; t1 = t2; t2 = temp;}
		//delt = (t2-t1)/4.0;
		  delt = (t2-t1)/partition_t;
		sum = 0.0;

/*
		for ( k = 0; k < gM; k++ )
		  {
			numb_in = k/4;
			ii      = k - 4 * numb_in; 
			kot     = (float)Kot[ii]; 
			xx      = t1 + delt * (numb_in + kot);
			x = a[0] + xx * d[0] + half;  // - 0.5;
			y = a[1] + xx * d[1] + half;  // - 0.5;
			z = a[2] + xx * d[2] + half;  // - 0.5;
			
			ix = x; xd = x - ix;
			iy = y; yd = y - iy;
			iz = z; zd = z - iz;

			//printf("\nix iy iz = %d %d %d ", ix, iy, iz);
			//values_t =  Volume->data[ix*nynz+iy*ny+iz];
			values_t = TrilinearInterpolation8(xd, yd, zd, 
			Volume->data[ix*nynz+iy*ny+iz],
		   	Volume->data[(ix+1)*nynz+iy*ny+iz],
            Volume->data[ix*nynz+(iy+1)*ny+iz],
		    Volume->data[(ix+1)*nynz+(iy+1)*ny+iz],
            Volume->data[ix*nynz+iy*ny+iz+1],
			Volume->data[(ix+1)*nynz+iy*ny+iz+1],
			Volume->data[ix*nynz+(iy+1)*ny+iz+1],
			Volume->data[(ix+1)*nynz+(iy+1)*ny+iz+1]);





			//printf("\nvalues_t = %f ", values_t);
			weit = (float)Weit[ii];
			sum = sum + weit*values_t;  //Gaussian quadrature.   

		  }

       sum = sum * delt ;

*/

                for ( k = 0; k < partition_t; k++ )
                   {
                     xx = t1+ delt * k;
                     x = a[0] + xx * d[0] + half; 
                     y = a[1] + xx * d[1] + half; 
                     z = a[2] + xx * d[2] + half;

                     ix = (int)floor(x); xd = x - ix;
                     iy = (int)floor(y); yd = y - iy;
                     iz = (int)floor(z); zd = z - iz;
/*
                     if( ix <0 || ix >= nx-1 || iy < 0 || iy >= ny-1 || iz < 0 || iz >= nz-1) 
                        {
			  printf("\nvolume projection error.ix iy iz=%d %d %d out ranged.", ix, iy, iz);
                          
                        }
*/
                        r1 = (ix<0 || iy<0 || iz<0 || ix > nx-1 || iy > nx-1 || iz > nx-1)?0.0:Volume->data[ix    * nynz + iy    * ny+iz];
                        r2 = (ix<-1|| iy<0 || iz<0 || ix > nx-2 || iy > nx-1 || iz > nx-1)?0.0:Volume->data[(ix+1)* nynz + iy    * ny+iz];
                        r3 = (ix<0 || iy<-1|| iz<0 || ix > nx-1 || iy > nx-2 || iz > nx-1)?0.0:Volume->data[ix    * nynz + (iy+1)* ny+iz];
                        r4 = (ix<-1|| iy<-1|| iz<0 || ix > nx-2 || iy > nx-2 || iz > nx-1)?0.0:Volume->data[(ix+1)* nynz + (iy+1)* ny+iz];
                        r5 = (ix<0 || iy<0 || iz<-1|| ix > nx-1 || iy > nx-1 || iz > nx-2)?0.0:Volume->data[ix    * nynz + iy    * ny+iz+1];
                        r6 = (ix<-1|| iy<0 || iz<-1|| ix > nx-2 || iy > nx-1 || iz > nx-2)?0.0:Volume->data[(ix+1)* nynz + iy    * ny+iz+1];
                        r7 = (ix<0 || iy<-1|| iz<-1|| ix > nx-1 || iy > nx-2 || iz > nx-2)?0.0:Volume->data[ix    * nynz + (iy+1)* ny+iz+1];
                        r8 = (ix<-1|| iy<-1|| iz<-1|| ix > nx-2 || iy > nx-2 || iz > nx-2)?0.0:Volume->data[(ix+1)* nynz + (iy+1)* ny+iz+1];


/*
if (ix<0 || iy<0 || iz<0 || ix > nx-1 || iy > nx-1 || iz > nx-1) printf("Signal of outlier1\n");
if (ix<-1|| iy<0 || iz<0 || ix > nx-2 || iy > nx-1 || iz > nx-1) printf("Signal of outlier2\n");
if (ix<0 || iy<-1|| iz<0 || ix > nx-1 || iy > nx-2 || iz > nx-1) printf("Signal of outlier3\n");
if (ix<-1|| iy<-1|| iz<0 || ix > nx-2 || iy > nx-2 || iz > nx-1) printf("Signal of outlier4\n");
if (ix<0 || iy<0 || iz<-1|| ix > nx-1 || iy > nx-1 || iz > nx-2) printf("Signal of outlier5\n");
if (ix<-1|| iy<0 || iz<-1|| ix > nx-2 || iy > nx-1 || iz > nx-2) printf("Signal of outlier6\n");
if (ix<0 || iy<-1|| iz<-1|| ix > nx-1 || iy > nx-2 || iz > nx-2) printf("Signal of outlier7\n");
if (ix<-1|| iy<-1|| iz<-1|| ix > nx-2 || iy > nx-2 || iz > nx-2) printf("Signal of outlier8\n");
*/

                     values_t = TrilinearInterpolation8(xd, yd, zd, r1, r2, r3, r4, r5, r6, r7, r8);

  
		     sum = sum + values_t * delt ;
                    }
		//getchar();

		prjimg[(i+half)*sample_num + j+half] = sum;
	  }

  /*
  for ( i = 0; i < sample_num; i++ )
	{
	  printf("\n");

	for ( j = 0; j < sample_num; j++ )
	  printf("%f ", prjimg[i * sample_num + j]);
	}
  */
  //getchar();



}


void Volume_GridProjection(Oimage* Volume, float rotmat[9], int sample_num, float *prjimg)
{
  int nx, ny, nz, nynz, i, j, k, x, y, z, ix, iy, iz, ii;
  float ox, oy, oz, old[3], Trotmat[9], X[3], xf, yf, weight;
  float value;

  nx = Volume->nx;
  ny = Volume->ny;
  nz = Volume->nz;

  nynz = ny * nz;

  ox = (nx-1)/2.0;
  oy = (ny-1)/2.0;
  oz = (nz-1)/2.0;
  //printf("\nnx ny nz=%d %d %d ox oy oz=%f %f %f ", nx, ny, nz, ox, oy, oz);
  Matrix3Transpose(rotmat, Trotmat); 
  for (i = 0; i < nynz; i++ ) prjimg[i] = 0.0;

  for ( z = 0; z < nz; z++ )
	{
	  old[2] = z - oz;
	  for ( x = 0; x < nx; x++ )
		{
		  old[0] = x - ox;
		  for ( y = 0; y < ny; y++ )
			{
			  old[1] = y - oy;

			  MatrixMultiply(Trotmat,3,3,old,3,1,X);
			  X[0] = X[0] + ox;
			  X[1] = X[1] + oy;
			  X[2] = X[2] + oz;
			  /*
			  //Nearest neighbor interpo.

			  ix = floor(X[0]+0.5);
			  iy = floor(X[1]+0.5);
			  iz = floor(X[2]+0.5);

			  if( ix >= 0 && ix <= nx-1 && iy >= 0 && iy <= ny-1 )
				{
				  xf = ix - X[0];
				  yf = iy - X[1];
				  weight = 1-sqrt(xf*xf+yf*yf);
				  if( weight > 0 )
					{
					  i  = x*nynz+ y*nz+z;
					  value = Volume->data[i];
					  ii = ix * ny + iy;
					  prjimg[ii] += weight*value;
					}
				}
				/// end Nearest neighbor interpo.
			  */
			  
			  ///Bilinear interpolation.
			  ix = X[0];
			  iy = X[1];
			  iz = X[2];

			  if( ix >= 0 && ix < nx-1 && iy >= 0 && iy < ny-1 )
				{
				  xf = X[0] - ix;
				  yf = X[1] - iy;
				  i  = x*nynz+ y*nz+z;

				  value = Volume->data[i];
				  //printf("\nxf yf=%f %f value=%f ", xf, yf, value);
				  for (j = ix; j < ix+2; j++ )
					{
					  xf = 1 - xf;
					  for ( k = iy; k < iy+2; k++ )
						{
						  yf = 1 - yf;
						  ii = j * ny + k;
						  prjimg[ii] += value*xf*yf;
						}
					}
				}
			  //end Bilinear interpolation.

			}
		}
	} 
}


void Imageinterpo_BiLinear(float *data, int nx, int ny)
{

  int i, j;
  int rnx, rny, rsize;
  float a, b, c, d;

  rnx   = nx + 1;
  rny   = ny + 1;
  rsize = rnx * rny;
 
  float *result = (float *)malloc(rnx * rny * sizeof(float));
  memset(result, 0, rsize*sizeof(float));

  
  for ( i = 0; i < rnx; i++ )
	for ( j = 0; j < rny; j++ )
	  {
		
		a =(i<1 || j < 1 || i >nx || j> nx  )?0:data[(i-1)*nx+j-1];
		b =(i<1 || j < 0 || i >nx || j> nx-1)?0:data[(i-1)*nx+j];
 		c =(i<0 || j < 1 || i >nx-1 || j>nx )?0:data[i*nx+j-1];
		d =(i<0 || j < 0 || i >nx-1 || j>nx-1)?0:data[i*nx+j];
        result[i*rny + j] =  0.25*(a+b+c+d);
 
      
	  }


  free(data);
  data = result;

}

void Imageinterpo_BiLinear(float *data, int nx, int ny, float *result)
{

  int i, j;
  int rnx, rny, rsize;
  float a, b, c, d;

  rnx   = nx + 1;
  rny   = ny + 1;
  rsize = rnx * rny;
for ( i = 0; i < rsize; i++ ) result[i] = 0.0;
for ( i = 0; i < nx; i++ )
        for ( j = 0; j < ny; j++ )
          {
           //result[i*rny + j] = data[i*nx+j];

		a =(i<1 || j < 1 || i >nx || j> nx  )?0:data[(i-1)*nx+j-1];
		b =(i<1 || j < 0 || i >nx || j> nx-1)?0:data[(i-1)*nx+j];
 		c =(i<0 || j < 1 || i >nx-1 || j>nx )?0:data[i*nx+j-1];
		d =(i<0 || j < 0 || i >nx-1 || j>nx-1)?0:data[i*nx+j];
        result[i*rny + j] =  0.25*(a+b+c+d);
          }
}



void WriteSpiFile(char* filename, float* data, int nx, int ny, int nz, float rot, float tilt, float psi)
{

   float min = 0.0, max = 0.0;

   SpiderHeader* sh = new SpiderHeader();


    sh->fNslice = (float)nz;
    sh->fNrow   = (float)ny;
    sh->fNcol   = (float)nx;

    sh->fIform  = (sh->fNslice > 1)?3:1;

    printf("\nsh->fNslice=%f , sh->fIform=%f sh->fNrow sh->fNcol=%f %f  x=%d ny=%d", sh->fNslice, sh->fIform, sh->fNrow, sh->fNcol, nx, ny);
    sh->fFmax   = max;
    sh->fFmin   = min;
    sh->fAv     = 0.0;

    sh->fSig    = -1;
    sh->fIhist  = 0 ;
    sh->fLabrec = (float) ceil((float)(256 / (float)sh->fNcol));

    sh->fPhi    = rot;
    sh->fTheta  = tilt;
    sh->fPsi    = psi;

    sh->fXoff   = 0;
    sh->fYoff   = 0;
    sh->fZoff   = 0;


    int imageSize =nx * ny* nz * sizeof(float);


    FILE* fp = fopen(filename, "wb");
    fwrite(sh, sizeof(SpiderHeader), 1, fp);
    fseek(fp, sizeof(SpiderHeader), SEEK_SET);
    fwrite(data, imageSize, 1, fp);


    fclose(fp);
    delete sh;


}

/*
void  readSpiFile(float* data, long memsize, char* spifile, EulerAngles* eulers)
{

 FILE* fp;
 if( (fp = fopen(spifile,"rb")) == NULL) return;
 SpiderHeader* sh = new SpiderHeader();
 fread(sh, sizeof(SpiderHeader), 1, fp);
 //printf("\nspiderheadersize=%d ", sizeof(SpiderHeader));

 eulers->rot   = sh->fPhi;
 eulers->tilt  = sh->fTheta;
 eulers->psi   = sh->fPsi;
 fseek(fp, sizeof(SpiderHeader), SEEK_SET);

 fread(data, memsize, 1, fp);


  delete sh;

  fclose(fp);

}
*/


size_t FREAD(void *dest, size_t size, size_t nitems, FILE *&fp, int reverse)
{
    size_t retval;
    if (!reverse)
        retval = fread(dest, size, nitems, fp);
    else
    {
        char *ptr = (char *)dest;
        int end = 0;
        retval = 0;
        for (int n = 0; n < nitems; n++)
        {
            for (int i = size - 1; i >= 0; i--)
            {
                if (fread(ptr + i, 1, 1, fp) != 1)
                {
                    end = 1;
                    break;
                }
            }
            if (end)
                break;
            else
                retval++;
            ptr += size;
        }
    }
    return retval;
}

/* This function reads a Spider volume. Call it as
 *
 *     int dim[3];
 *         float *data=NULL;
 *             readSpiderFile("myfile.vol",'V',dim,&data);
 *                 readSpiderFile("myfile.img",'I',dim,&data);
 *                     
 *                         Code errors:
 *                             0 - OK
 *                                 1 - Cannot open file
 *                                     2 - File is not a Spider volume
 *                                         3 - Problem when computing the size of the file
 *                                             4 - filetype is not 'V' or 'I'
 *                                             */
int readSpiderFile(const char *filename, char filetype,
    int dim[], float **data, EulerAngles* eulers)
{
    FILE* fp=NULL;
    union {
        unsigned char c[4];
        float         f;
    } file_type;
    int fileReversed=0, machineReversed=0, reversed=0;
    SpiderHeader header;
    unsigned long usfNcol, usfNrow, usfNslice, usfHeader, size, tmpSize,
        volSize, n, currentPosition;
    struct stat info;
    float *ptr;

    /* Set dimensions to 0, just in case it cannot be read */
    dim[0]=0;
    dim[1]=0;
    dim[2]=0;    

    /* Check that the filetype is correct */
    if (filetype!='V' and filetype!='I')
        return 4;

    /* Open file */
    if ((fp = fopen(filename, "rb")) == NULL)
        return 1;

    /* Check that the input file is really a Spider volume */
    currentPosition=ftell(fp);

    /* Check file type */
    fseek(fp, currentPosition+16, SEEK_SET);
    for (int i = 0; i < 4; i++)
        fread(&(file_type.c[i]), sizeof(unsigned char), 1, fp);
    fseek(fp,  currentPosition+0, SEEK_SET);
    switch (filetype)
    {
        case 'V':
            if (file_type.c[0]==  0 && file_type.c[1]== 0 &&
                file_type.c[2]== 64 && file_type.c[3]==64)
            {
                fileReversed=0;
            }
            else if (file_type.c[0]==64 && file_type.c[1]==64 &&
                     file_type.c[2]== 0 && file_type.c[3]== 0)
            {
                fileReversed=1;
            }
            else
            {
                fclose(fp);
                return 2;
            }
            break;
        case 'I':
            if (file_type.c[0]==  0 && file_type.c[1]== 0 &&
                file_type.c[2]==128 && file_type.c[3]==63)
            {
                fileReversed=0;
            }
            else if (file_type.c[0]==63 && file_type.c[1]==128 &&
                     file_type.c[2]== 0 && file_type.c[3]==  0)
            {
                fileReversed=1;
            }
            else
            {
                fclose(fp);
                return 2;
            }
            break;
    }

    /* Check machine type */
    file_type.f = 1;
    if (file_type.c[0]==63 && file_type.c[1]==128 && file_type.c[2]==0 &&
        file_type.c[3]==0)
        machineReversed=1;

    /* Read header */
    reversed=fileReversed ^ machineReversed;
    if (!reversed)
        fread(&header, sizeof(SpiderHeader), 1, fp);
    else
    {
        FREAD(&header,             sizeof(float),  36, fp, true);
        FREAD(&header.fGeo_matrix, sizeof(double),  9, fp, true);
        FREAD(&header.fAngle1,     sizeof(float),  13, fp, true);
        FREAD(&header.fNada2,      sizeof(char),  756, fp, true);
    }

    /* Compute file size, header size and volume dimensions */
    usfNcol = (unsigned long) header.fNcol;
    usfNrow = (unsigned long) header.fNrow;
    usfNslice = (unsigned long) header.fNslice;
    usfHeader = (unsigned long) header.fNcol *
                (unsigned long) header.fLabrec * sizeof(float);
    if (fstat(fileno(fp), &info))
    {
        fclose(fp);
        return 3;
    }

    /* Readjust the number of rows in "aberrant" images*/
    if (filetype=='I' || header.fIform == 1)
        if (usfNcol*usfNrow*sizeof(float) == info.st_size)
        {
            usfNrow = (unsigned long)(--header.fNrow);
            --header.fNrec;
        }

    /* Check that the file size is correct */
    switch (filetype)
    {
        case 'I':
            size = usfHeader + usfNcol * usfNrow * sizeof(float);
            if ((size != info.st_size) || (header.fIform != 1))
            {
                fclose(fp);
                return 2;
            }
            break;
        case 'V':
            size = usfHeader + usfNslice * usfNcol * usfNrow * sizeof(float);
            if ((size != info.st_size) || (header.fIform != 3))
            {
                fclose(fp);
                return 2;
            }
            break;
    }
    
    /* Read the extra filling header space */
    tmpSize = (unsigned long)(header.fNcol * header.fLabrec * 4);
    tmpSize -= sizeof(SpiderHeader);
    for (unsigned int i = 0; i < tmpSize / 4; i++)
    {
        float tmp;
        fread(&tmp, sizeof(float), 1, fp);
    }
    currentPosition=ftell(fp);

    /* Fill the dimensions */
    dim[0]=(int)header.fNcol;
    dim[1]=(int)header.fNrow;
    if (filetype=='V') dim[2]=(int)header.fNslice;
    else               dim[2]=1;
    volSize = (unsigned long) dim[0]*dim[1]*dim[2];


    if(filetype == 'I' ) {
      eulers->rot   = header.fPhi;
      eulers->tilt  = header.fTheta;
      eulers->psi   = header.fPsi;
    }

    /* Read the whole file */
    //if (*data!=NULL) free(data);

    *data=(float *) malloc(volSize*sizeof(float));
    ptr=*data;
    for (n=0; n<volSize; n++)
        FREAD(ptr++, sizeof(float), 1, fp, reversed);

    /* Close file */
    fclose(fp);
    return 0;
}




void ReadRawIVFile(const char* filename, const char* path, float *data)
{
/*
 RawIVHeader *rawivHeader=new RawIVHeader();
 FILE* fp;
 unsigned int memsize,nx, ny, nz;
  //QDir dataDir;
  //dataDir.setCurrent(path);
  //printf("\nnewcurrentpath=%s ", dataDir.current().path().ascii());
  printf("\nfilename =%s ", filename);

  if( (fp = fopen(filename,"rb")) == NULL)
        {printf("\nerror read files. ");return ;}

 fread(rawivHeader, sizeof(RawIVHeader), 1, fp);

 fseek(fp, sizeof(RawIVHeader), SEEK_SET);

 memsize = uint64(rawivHeader->dim[0]) * uint64(rawivHeader->dim[1]) * uint64(rawivHeader->dim[2])*sizeof(float);
 nx = uint64(rawivHeader->dim[0]);
 ny = uint64(rawivHeader->dim[1]);
 nz = uint64(rawivHeader->dim[2]);   
 printf("\nnx ny nz = %d %d %d RawIV  memsize ================%d ", nx, ny, nz, memsize);
 fread(data, memsize, 1, fp);

 delete rawivHeader;
 fclose(fp);
*/
}


