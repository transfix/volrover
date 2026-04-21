/*
All right are reserved by The National Key Lab of Scientific and Engineering Computing,
Chinese Academy of Sciences.

        Author: Ming Li <liming@lsec.cc.ac.cn>
        Advisor: Guoliang Xu <xuguo@lsec.cc.ac.cn>

This file is part of CIMOR. CIMOR stands for:
        Computational Inverse Methods of Reconstruction

*/

//Declaration.==================================================================================
/*Bscale:  only right for 1 or 2.
*/ 

#include <stdio.h>
#include <stdlib.h>

#include <Reconstruction/Reconstruction.h>
#include <ctype.h>
#include <math.h>
#include <fftw3.h>
#include <iostream>
#include <sys/time.h>
#include <qdir.h>

using namespace std;


void CubeCoefficients();
void Sortingfloat(float *point,int *index, int number);

// global variables, Xu added
float Cube_Coeff[216];
int PRO_LENGTH, PRO_LENGTH_SUB, SUB, THICK, BStartZ, BFinishZ;

Reconstruction::Reconstruction()
{
  tempf        = NULL;
  slicemat     = NULL;
  ReBFT        = NULL;
  bgrids       = NULL;
  bspline      = NULL;
  BsplineFT    = NULL;
  Matrix       = NULL;
  coefficients = NULL;
  gd           = NULL;
  VolImg       = NULL;
  //VolF         = NULL;
  OrthoOrdered = NULL;
  
  Rmat         = NULL;
  InvRmat      = NULL;
  gd_FFT       = NULL;
  gdxdphi      = NULL;
  Mat          = NULL;
  O_Kcoef      = NULL;
}

Reconstruction::~Reconstruction()
{
  //delete bspline;
}

// Xu cahnged this code 1004
void Reconstruction::setThick(int thickness)
{
  THICK = thickness;
  THICK = THICK/2;
  THICK = THICK*2 + 1;
  
  if (THICK > BFinishX - BStartX + 1) THICK = BFinishX - BStartX + 1;
  
  BStartZ  = -THICK/2; 
  BFinishZ =  THICK/2;
  
  printf("\nnew THICK==================%d, BStartZ = %d, BFinishZ = %d, BStartX = %d, BFinishX = %d ", THICK, BStartZ, BFinishZ, BStartX, BFinishX);
}

/***************************************************************************
 Descriptions:
     The main program for reconstruction.

****************************************************************************/
Oimage* Reconstruction::Reconstruction3D_Ordered(PrjImg* gd, int iter, float tau, 
                         EulerAngles* eulers, Views* nview, int m_Itercounts, int object)
{

  float *coef, *O_coef_Result, *Non_o_coef, *f, rotmat[9];
  int i, ui;
  clock_t t,t_end ;
  float time,time_end ;
  
  // InitialFunction(ImgNx, ImgNy, ImgNz);
  
  coef   = (float *)malloc(usedtolbsp * sizeof(float));
  O_coef_Result  = (float *)malloc(usedtolbsp * sizeof(float));
  Non_o_coef  = (float *)malloc(usedtolbsp * sizeof(float));
  f      = (float *)malloc(tolbsp * sizeof(float));
  
  if(m_Itercounts == 1 )
    { 
      // ObtainRotationMatFromViews(Rmat, nview);
      //EulerMatrice(Rmat,eulers);
      //for ( i = 0; i < nv * 9; i++ ) InvRmat[i] = Rmat[i];	  
      //InverseRotationMatrix(InvRmat, nv);	
      
      printf("\nCreate the Gram Schmidt processing matrix.\n");
      bspline->GramSchmidtofBsplineBaseFunctions2();
      
      printf("\nEvaluate Ortho Bspline Basis at volume grids.\n ");
      bspline->Evaluate_OrthoBSpline_Basis_AtVolGrid2();
      
      //FFT_gdMatrix(gd);
      
      compare_Xdphi_Phi_pd(nview);
      printf("\ncompare finish.");getchar();
      
      // EvaluateBsplineBaseFT();
      // TestXdF_Gd();
      // testxdf_gd(nview);
      
      // printf("test OK."); getchar();
      
      //printf("\nFFT_gdMatrix over.");
      //EvaluateFFT_GdXdPhiMatrix();
      //printf("\nFFT_GdXdPhiMatrix over.");
      
      //EvaluateFFT_GdXdOrthoPhiMatrix();
      //printf("\nFFT_GdXdOrthoPhiMatrix over.");
      
      //printf("\nEvaluate BsplineFT integral ijk_pqr");
      //t=clock();
      //time=(float)t/(float)CLOCKS_PER_SEC;
      //Evaluate_BsplineFT_ijk_pqr();
      
      //t_end=clock();
      //time_end=(float)t_end/(float)CLOCKS_PER_SEC;
      //time=time_end-time;
      //printf("\n Time = %f \n",time);
      
      //InitialFunction(ImgNx, ImgNy, ImgNz);
      
      for ( ui = 0; ui <usedtolbsp; ui++ )
	O_Kcoef[ui]   = coefficients[ui];
      
      free(coefficients);
      // free(gd_FFT);
      // free(slicemat);
      
    }
  
  
  
  Oimage* Object = ReconstructionIters(iter,tau,coef, O_coef_Result, Non_o_coef, f);
  
  /*
    free(f);
    free(coef);
    free(O_coef_Result);
    free(Non_o_coef);
    
    return Object;
  */
  
}

// Main program
Oimage* Reconstruction::Reconstruction3D(int ReconManner, int iter, float tau, EulerAngles* eulers, int m_Itercounts, int object)
{
  
  float         *coef=NULL, *O_coef_Result=NULL, *Non_o_coef=NULL, *f=NULL, rotmat[9];
  int           i, ui;
  clock_t       t,t_end ;
  float        time,time_end ;
  Oimage        *image = NULL, *vol=NULL;
  
  
  printf("\nVolImgSize = %d ", VolImgSize);
  coef   = (float *)malloc(usedtolbsp * sizeof(float));    // Bspline coeffieicent previous step
  O_coef_Result  = (float *)malloc(usedtolbsp * sizeof(float));    // Result 
  Non_o_coef  = (float *)malloc(usedtolbsp * sizeof(float));    // Non-orthogonal coeff
  f      = (float *)malloc(VolImgSize * sizeof(float));    // outpot volume data
  
  //vol    = InitImageParameters(1, ImgNx,ImgNy,ImgNz, 1);
  
  if (m_Itercounts == 1 ) { 
    
    //bspline setting--------------------------------------
    
    printf("\nStart the Gram Schmidt processing matrix.");
    bspline->GramSchmidtofBsplineBaseFunctions2();        // othogonalization
    printf("Finished GramSchmidtofBsplineBaseFunctions. \n");
    /*
      
      printf("4 Test by Xu----------------------------------------------------------------\n");
      printf("\nEvaluate Ortho BsplineBasis at volume grid. ");
      bspline->Evaluate_OrthoBSpline_Basis_AtVolGrid2();    // OBasis function values on grid points
      printf("Finished. \n");
      
      printf("5 Test by Xu----------------------------------------------------------------\n");
      //spline->Evaluate_BSpline_Basis_AtGaussNodes();
      printf("\nEvaluate Ortho BsplineBasis at Gauss Nodes. ");
      bspline->Evaluate_OrthoBSpline_Basis_AtGaussNodes();  //  OBasis function values on Gaussian points
      printf("Finished. \n");
      
    */
    
    printf("\nStart Evaluate Ortho BsplineBasis at Image Grid. ");
    bspline->Evaluate_BSpline_Basis_AtImgGrid();          // ???
    printf("Finished. \n");
    //getchar();
    
    printf("\nStart Evaluate Ortho BsplineBasis at sub Image Grid. ");
    bspline->Evaluate_BSpline_Basis_AtImgGrid_sub();          // ???
    printf("Finished. \n");
    
    //set transform matrix-----------------------------------------
    //   if(ReconManner == 1 ) EulerMatrice(Rmat,eulers);                            // Compute ratation matroces, 36 3*3 matrices
    /*
      for ( i = 0; i < nv * 9; i++ ) InvRmat[i] = Rmat[i];
      printf("\n Start InverseRotationMatrix----------------------------------------------------------------\n");
      InverseRotationMatrix(InvRmat, nv);                   // Inverse all the ratation matrices
      printf("\n End InverseRotationMatrix----------------------------------------------------------------\n");
    */
    
    // set initial function-----------------------------------
    //printf("\n Start InitialFunction----------------------------------------------------\n");
    //InitialFunction(0); 
    //printf("\n Finish InitialFunction----------------------------------------------------\n");
    //VolF = InitialFunction(3); 
    // Get ortho coefficients of initial function bspline.
    
    printf("\n Start ConvertToOrthoCoeffs----------------------------------------------------\n");
    ConvertToOrthoCoeffs();
    printf("\n Finish ConvertToOrthoCoeffs----------------------------------------------------\n");
    
    //set object phantom--------------------------------------
    
    //image = Phantoms_Volume(object); // object to be constructed, 0,cylinder, 
    //Phantoms_gd(image);     
    
    if(ReconManner == 1 ) {
      EulerMatrice(Rmat,eulers); 
      printf("\n Start Phantoms_Volume----------------------------------------------------\n");
      VolImg = Phantoms_Volume(object);                                                    // 1, big shpere + small sphere
      printf("\n End  Phantoms_Volume----------------------------------------------------\n");
      //VolImg = VolF;
      //Test_GridProjectioin(eulers);
      
      printf("\n Start Phantoms_gd----------------------------------------------------\n");
      Phantoms_gd(VolImg, eulers); 
      printf("\n Finish Phantoms_gd----------------------------------------------------\n");
    }        
    
    if (OrderManner == 1 ) {
      OrthoOrdered = (int *)malloc(nv*sizeof(int));
      setOrders(nv);
      CurrentId = 0;
    }
    
    /*  
	printf("\nStart FFT_gdMatrix----------------------------------------------------------------\n");
	FFT_gdMatrix();                                       // FFT of 2d images
	printf("\nFinish  FFT_gdMatrix");
    */	  
    t=clock();
    time=(double)t/(double)CLOCKS_PER_SEC;
    
    //if (Recon_Method == 2 )     //error. for Recon_Method =1. for ComputeTau need VolBsplineProjections.
    //{
    printf("\n Start VolBsplineProjections----------------------------------------------------\n");
    VolBsplineProjections();
    //VolBsplineProjections_FA();
    //getchar();
    printf("\n Finish VolBsplineProjections. This code is slow--------------------------------\n");
    //}
    
    if(Recon_Method == 1 )
      {
	printf("\nStart FFT_gdMatrix----------------------------------------------------------------\n");
	FFT_gdMatrix();                                       // FFT of 2d images
	printf("\nFinish  FFT_gdMatrix");
	
	printf("\nFFT_GdXdPhiMatrix.");
	EvaluateFFT_GdXdPhiMatrix();                           // compute FFT of g_dX_di --  a slice page 11 
	printf("Finished. \n");
	printf("8 Test by Xu----------------------------------------------------------------\n");
	
	
	//printf("\nFFT_GdXdOrthoPhiMatrix.");
	//EvaluateFFT_GdXdOrthoPhiMatrix();                      // zhengjiao Hua
	//printf("Finished. \n");	  
	printf("\nEvaluate BsplineFT integral ijk_pqr");       
	Evaluate_BsplineFT_ijk_pqr();                          // (4.14)
	printf("Finished. \n");
	
	printf("9 Test by Xu----------------------------------------------------------------\n");
      }
    
    t_end=clock();
    time_end=(double)t_end/(double)CLOCKS_PER_SEC;
    time=time_end-time;
    printf("\n Time = %f \n",time);
    
    // InitialFunction(1);                                    // initial function, 0--sphere, 1 -- ??
    
    //for ( ui = 0; ui <usedtolbsp; ui++ )
    //	O_Kcoef[ui]   = coefficients[ui];                  // O_Kcoeff is current coeff 
    
    //  free(coefficients);
    // free(gd_FFT);
    // free(slicemat);
    //kill_all_but_main_img(image);
    //free(image);
    
    
    if(newnv > nv ) newnv = nv;
    for ( i = 0; i < nv; i++ ) { EJ1_k[i] = 0.0; EJ1_k1[i] = 0.0; } 
    //free(VolImg->data);
    //free(VolImg);
  }
  
  // iteration in reconstruction
  printf("\n Start ReconstructionIter----------------------------------------------------\n");
  Oimage* Object = ReconstructionIters(iter,tau,coef, O_coef_Result, Non_o_coef, f);
  printf("\n Finish ReconstructionIter----------------------------------------------------\n");
  
  free(f);
  free(coef);
  free(O_coef_Result);
  free(Non_o_coef);
  
  f     = NULL;
  coef  = NULL;
  O_coef_Result = NULL;
  Non_o_coef = NULL;
  
  return Object;
}




Oimage* Reconstruction::Reconstruction3D_FBP(int iter, float tau, EulerAngles* eulers, int m_Itercounts, int object)
{
  int i;
  Oimage  *image = NULL, *vol=NULL;


  image = Phantoms_Volume(object);
  
  Phantoms_gd(image, NULL);

  //printf("\nmemmory");
  //getchar();

  kill_all_but_main_img(image);

  free(image);
  //printf("\nmemmory");
  //getchar();

  EulerMatrice(Rmat,eulers);
  for ( i = 0; i < nv * 9; i++ ) InvRmat[i] = Rmat[i];
  
  InverseRotationMatrix(InvRmat, nv);
  
  vol    = InitImageParameters(1, ImgNx,ImgNy,ImgNz, 1);

  backprojection(vol);
  
  return vol;
   
}


/***************************************************************************
Descriptions:
    The main iter process for reconstruction. 

****************************************************************************/
Oimage* Reconstruction::ReconstructionIters(int iter, float tau, float *coef, 
                        float *O_coef_Result, float *Non_o_coef, float *f )
// coef -- k step orthogonal B-spline coefficients
// O_coef_Result-- k+1 step orthogonal B-spline coefficients.
// Non_o_coef-- Non-Othogonal B-spline coefficients
// f    -- volume data
{
float         taual,taube, tauga, taula, IsoC;
int           it, i, j, k, i1, j1, k1,n, usedN2, half;
float         ai, bj, ck, aibj, aibjck, tau2, Xdfphi, *dcoefs, *SchmidtMat;
int           ii,jj, iii, jjj,  usedN1, *index, *subset;
vector<int>   suppi , outliers;
vector<CUBE>  cubes, cubes1;

//fftw_complex  *in=NULL, *out=NULL;   // complex data for 3d FFT

clock_t t,t_end, t1, t1_end;
float time,time_end, time1, time1_end;
float maxf, minf;
int    sub, subImgNx, subVolImgSize, subimg2dsize, subhalf, a, b;
float *subf = NULL, weigth, fill_value, comp_tau, total;
int    padfactor, padVolImgSize, padImgNx;
float *J234phi,  *J5phi;
//int   flows;
 //float *J234phi, *J3phi, *J4phi, *J5phi;

// set intial values
half          = ImgNx/2; 
sub           = 1;
// sub           = 2;    // Our methods.  subdivision volume f for compute central slice.

padfactor     = 1;    // xmipp,  bsoft method. for interpolate in fourier space. 1: do not interpolate. 
 //padfactor     = 2;  
if(padfactor == 2 ) sub = 1;

printf("\n Xu-- ImgNx = %d, sub = %d, N = %d fac=%f\n", ImgNx, sub, N, fac);

subImgNx      = sub * ImgNx;
subimg2dsize  = (subImgNx+1) * (subImgNx+1);
subVolImgSize = (subImgNx+1) * (subImgNx+1) * (subImgNx+1);
subhalf       = subImgNx/2;
a             = -half;
b             = half;
weigth        = (b-a)*1.0/subImgNx; 
weigth        = weigth * weigth * weigth;

padVolImgSize = padfactor*ImgNx * padfactor * ImgNx * padfactor * ImgNx;
padImgNx      = padfactor*ImgNx - 1;

//flows = 1;
//flows = 2;
//flows = 3;


// malloc memeries
// in   = (fftw_complex *)fftw_malloc(subVolImgSize * sizeof(fftw_complex));
//out  = (fftw_complex *)fftw_malloc(subVolImgSize * sizeof(fftw_complex));
//subf = (float *)malloc(subVolImgSize * sizeof(float));
//float* padf = (float *)malloc(padVolImgSize * sizeof(float));

/* 
//Compute J1 by using FFT and Central slice theorem.
in   = (fftw_complex *)fftw_malloc(subImgNx*subImgNx*subImgNx*sizeof(fftw_complex));
out  = (fftw_complex *)fftw_malloc(subImgNx*subImgNx*subImgNx*sizeof(fftw_complex));
*/


//fftw_complex* in_sub   = (fftw_complex *)fftw_malloc(subImgNx*subImgNx*subImgNx*sizeof(fftw_complex));
//fftw_complex* out_sub  = (fftw_complex *)fftw_malloc(subImgNx*subImgNx*subImgNx*sizeof(fftw_complex));

xdgdphi   = (float *)malloc(usedtolbsp*sizeof(float));
J234phi   = (float *)malloc(usedtolbsp*sizeof(float));
//J3phi   = (float *)malloc(usedtolbsp*sizeof(float));
//J4phi   = (float *)malloc(usedtolbsp*sizeof(float));
J5phi     = (float *)malloc(usedtolbsp*sizeof(float));
dcoefs    = (float *)malloc(usedtolbsp*sizeof(float));
index     = (int *)malloc(nv*sizeof(int));
subset    = (int *)malloc(newnv*sizeof(int));

 
SchmidtMat= (float *)malloc((usedN+1)*(usedN+1)*sizeof(float));
maxf = 0.0;
fill_value = 0.0;


for ( i = 0; i < (usedN+1)*(usedN+1); i++ ) SchmidtMat[i] = (float)bspline->SchmidtMat[i];

// sortf not used. 
/*
taual = tau * al/tolbsp;
taube = tau * be/tolbsp;
tauga = tau * ga/tolbsp;
taula = tau * la/tolbsp;
tau2  = 2.0 * tau * reconj1/tolbsp;  
*/

//tau   = 1.0;      //for new temporal step-size computing. 8-31-09
taual = tau * al;
taube = tau * be;
tauga = tau * ga;
taula = tau * la;
//tau2  = 2.0 * tau * reconj1;
tau2  = 2.0 * reconj1;
total = al + be + ga + la; 

 
printf("\nreconj1 = %f taual=%f taube=%f tauga=%f taula=%f alpha=%f", reconj1, al, be, ga, la, alpha);  
printf("\nTotal view number==============%d ", nv);


usedN1 = usedN+1; 
usedN2 = (usedN + 1) * (usedN + 1); // The used B-splines for 2d == (N-3)^2


//Initialize coef.
printf("\nInitialize the coefficients of Bspline Basis.\n");
for ( i = 0; i < usedtolbsp; i++ ) {
   coef[i]= O_Kcoef[i];   // Orthogonal coefficients
   Non_o_coef[i]  = 0.0;
   //Matrix[i] = 0.0;                                                               //Compute variation J1. Method 0.
   O_coef_Result[i]  = 0.0;
}
   
printf("\nBegin the iteration processing for reconstruction.\n ");
t=clock();
time=(double)t/(double)CLOCKS_PER_SEC;
 printf("++++ t = %f, time = %f, (float)CLOCKS_PER_SEC = %e\n", (double)t, time, (double)CLOCKS_PER_SEC);
printf("\nnewnv =========%d bandwidth=%d", newnv, bandwidth);

// Start iterations    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
for ( it = 0; it < iter; it++ ) {
   printf("\n");
   printf("\niteration step = %d, flows = %d, reconj1 = %f, taube = %f\n", it + 1, flows, reconj1, taube);
   
	 
   for ( i = 0; i < VolImgSize; i++ ) f[i] = 0.0;
//   for ( i = 0; i < subVolImgSize; i++ ) subf[i] = 0.0;

    
   ConvertCoeffsToNewCoeffs(coef, Non_o_coef, bspline->SchmidtMat);//convert to non-ortho bspline basis coefficients.

//   for ( i = 0; i < usedtolbsp; i++ ) {
//      if (Non_o_coef[i] < 0.0) Non_o_coef[i] = 0.0;   // Xu Cut small coefficients 0912, the result is not good
//   } 	 

   // obtain volume data from non-ortho bspline function
   if( Bscale == 1.0 && sub == 1 ) bspline->ObtainObjectFromNonOrthoCoeffs(Non_o_coef,f);
   if( Bscale == 2.0 && sub == 1 ) bspline->ObtainObjectFromNonOrthoCoeffs_FA(Non_o_coef,f);
   //if( sub == 2 ) bspline->ObtainObjectFromNonOrthoCoeffs_sub(Non_o_coef,subf, sub);
	 
	  ////////////////////////pad methods.
          /////////////////////////ImgPad(f, padf, ImgNx+1, padImgNx+1,fill_value);
	  /////////////////////////FFT_padf(padf, padImgNx+1);
	  ////////////////////////// end pad methods.

   for ( i = 0; i < usedtolbsp; i++ )  {
      xdgdphi[i] = 0.0;
      dcoefs[i]  = 0.0;    // Xu added 0912
      J234phi[i] = 0.0;    // Xu moved 0912
      J5phi[i]   = 0.0;    // Xu movwd 0912
   }

   //Compute variation J1. Method 2.  by using FFT and Central slice theorem.------------------------------------------------
   /*
	 if(sub == 1 ) ComputeIFSFf_gd(f, padf, padfactor, ImgNx+1, padImgNx+1, sub, fill_value, Non_o_coef, in, out);
	 if(sub == 2 ) ComputeIFSFf_gd(subf, padf, padfactor, ImgNx+1, padImgNx+1, sub, fill_value, Non_o_coef, in, out);
   */


   //Compute variation J1. Method 0.-----------------------------------------------------------------------------------------
   //EvaluateCubicBsplineBaseProjectionIntegralMatrix(Non_o_coef);

   //Compute variaition J1. Method 1.----------------------------------------------------------------------------------------
   if(reconj1 != 0.0 ) {
      t1=clock();
      time1=(double)t1/(double)CLOCKS_PER_SEC;
      //printf("Start ComputeXdf_gd\n");

      // Xu added this line for test 0918
      //ComputeEJ1_Acurate(f, Non_o_coef);
      if (newnv == nv ) 
       {
         if(Recon_Method == 1) EvaluateCubicBsplineBaseProjectionIntegralMatrix(Non_o_coef);
         if(Recon_Method == 2) ComputeXdf_gd(f, Non_o_coef);

       }
  
      else {
         /*select order sub set from all projection sets.-------------------------------------*/ 
         // method 1. Compute energy J1.
        if (OrderManner == 0 ) {
            printf("Start ComputeEJ1.\n");
            ComputeEJ1(f, Non_o_coef, index);
            
            // Test 
            t1_end=clock();
            time1_end=(double)t1_end/(double)CLOCKS_PER_SEC;
            printf("End ComputeEJ1.                                   time=%f\n", time1_end-time1);
         }

         //method 2. employ ortho ordered subset.
         if (OrderManner == 1 ) {
            for ( i = 0; i < newnv; i++ ) subset[i] = i;
            getSubSetFromOrthoSets(subset);
            for ( i = 0; i < newnv; i++) index[nv-newnv+i] = subset[i]; 
            printf("End Ortho Order.                                   time=%f\n", time1_end-time1);
         }

         t1=clock();
         time1=(double)t1/(double)CLOCKS_PER_SEC;
         if(Recon_Method == 2) ComputeXdf_gd(newnv, index, f, Non_o_coef);
      }  // end if (newnv == nv ) 

      t1_end=clock();
      time1_end=(double)t1_end/(double)CLOCKS_PER_SEC;
      printf("End ComputeXdf_gd.                                time=%f, newnv, nv = %d, %d, Start and End %d, %d\n", time1_end-time1, newnv, nv,
                      BStartX, BFinishX);
   }

//   ConvertCoeffsToNewCoeffs(xdgdphi, J5phi, bspline->InvSchmidtMat);
//   for ( i = 0; i < usedtolbsp; i++ ) printf("J5phi=%f ", J5phi[i]);
//   for ( i = 0; i < usedtolbsp; i++ ) xdgdphi[i] = J5phi[i];
   //Finish Compute variation J1.-----------------------------------------------------------------------------------------------

   // Xu comment out 0912, move forward 
   /*
   for ( i = 0; i < usedtolbsp; i++ ) {
      J234phi[i] = 0.0; 
      J5phi[i] = 0.0;
   }
   */

   if(taual != 0.0 || tauga != 0.0 || taula != 0.0 ) {
      EestimateIsoValueC(f, &IsoC,VolImgSize, cubes,cubes1, fac);  
      n = cubes.size();
      printf("\ncubes.size=%d IsoC =%f ", n, IsoC);
   

      //compute J4 part 2. outliers. 
      if(tauga != 0.0 ) {
         ComputeHc_f(f, &IsoC, outliers);  
         ComputeJ4outliers(J234phi, f, outliers, tauga); 
      }

      CubeCoefficients();
      //subdivide narrow band.
      GetSubNarrowBand(cubes, cubes1, &IsoC, Non_o_coef, subf, fac);
      //printf("\nEnd GetSubNarrowBand.\n");
      n = cubes1.size();    
      printf("\ncubes1.size=%d ", n);
      ComputeJ2J4J5(Non_o_coef, &IsoC, n, cubes1, V0, taual, tauga, taula, J234phi, J5phi);

//      ConvertCoeffsToNewCoeffs(J234phi, J5phi, bspline->InvSchmidtMat);
//      for ( i = 0; i < usedtolbsp; i++ ) {J234phi[i] = J5phi[i];J5phi[i] = 0.0;}
   }


   // Compute J_3
   if(taube != 0.0 ) {

      t1=clock();
      time1=(double)t1/(double)CLOCKS_PER_SEC;

      //ComputeJ3_Xu(Non_o_coef, taube, J234phi);    
      if (flows == 1) ComputeJ3_HTF(Non_o_coef, 1.0, J234phi);   // Xu changed taube to one 0920 
      if (flows == 2) ComputeJ3_MCF(Non_o_coef, 1.0, J234phi);   // Xu changed taube to one 0920 
      if (flows == 3) ComputeJ3_WMF(Non_o_coef, 1.0, J234phi);   // Xu changed taube to one 0920 

      t1_end=clock();
      time1_end=(double)t1_end/(double)CLOCKS_PER_SEC;
      printf("End Compute J3.                                   time=%f\n", time1_end-time1);
   }


   // Xu added the following 20 lines for automatically choosing taube, 0914
   
   float inner1, inner2, inner3, alpha;
   if (reconj1 != 0.0 && taube != 0.0) {

      inner1 = 0.0;
      inner2 = 0.0;
      inner3 = 0.0;
      for (i = 0; i < usedtolbsp; i++) { 
         inner1 = inner1 + xdgdphi[i]*xdgdphi[i];
         inner2 = inner2 + J234phi[i]*J234phi[i];
         inner3 = inner3 + xdgdphi[i]*J234phi[i];
      }

      
      inner1 = sqrt(inner1);
      inner2 = sqrt(inner2);
      alpha = inner1/inner2;
      //if (alpha > 2.0) alpha = 2.0;
      
      printf("Inner1, inner2, %f, %f, alpha = %f, Angle = %f\n", inner1, inner2/taube, alpha,
              acos(inner3/(inner1*inner2))*180/3.1415926);
      for (i = 0; i < usedtolbsp; i++) { 
         J234phi[i] = alpha*J234phi[i];
      }
     


      /*
      if (inner3 >= 0) {
         for (i = 0; i < usedtolbsp; i++) { 
            xdgdphi[i] = 0.0; 
         }
      }  else {
         alpha = - inner1/inner3;
         for (i = 0; i < usedtolbsp; i++) { 
            J234phi[i] = alpha*J234phi[i];
         }
         printf("Inner1, inner2, inner3 = %f, %f, %f alpha = %f\n", inner1, inner2, inner3, alpha);
      }
      */
   }
  

   int i1u, j1u, k1u, ui, i2, j2, k2;//,bandwidth; 
   float sum,sum1;

   //bandwidth = 3;  //The result is already good
   //bandwidth = 4;
   //bandwidth = 5;
   //bandwidth = 6;  
   //bandwidth = 7;  
   ////////////////bandwidth = 8;  // good enough 
   //bandwidth = 10; 
   //bandwidth = 12; 
   //bandwidth = 14; 
   //bandwidth = 16; 
   //bandwidth = 18;  // almost exact 
   //bandwidth = 20;  // almost exact 
   //printf("\n\n");
   // compute J1
   if (reconj1 != 0.0) { 
      // Loop for basis 
      t1=clock();
      time1=(double)t1/(double)CLOCKS_PER_SEC;

      for ( i1 = 0; i1 <= usedN; i1++ ) {
         i1u = i1*usedN1;
         //printf("\ni1=%d ", i1);
         for ( j1 = 0; j1 <= usedN; j1++ ) {
            j1u = j1*usedN1;
            for ( k1 = 0; k1 <= usedN; k1++ ) {
               
               k1u = k1*usedN1;
               jj  = i1*usedN2+j1*(usedN+1) + k1;
	       //O_coef_Result[jj] = coef[jj];   // Xu  0912
               //Xdfphi = 0.0;
               
               // Orthogonalize
               i2 = i1 - bandwidth; 
               if (i2 < 0) i2 = 0;
               //for ( i = 0; i <= i1; i++ ) {
               for ( i = i2; i <= i1; i++ ) {
	          //ai  = (float)bspline->SchmidtMat[i1u+i] ;
	          ai  = SchmidtMat[i1u+i] ;
//printf("\n i1 = %d, j1 = %d, k1 = %d, i = %d, ai = %e\n", i1, j1, k1, i, ai);
	          iii = i*usedN2;

                  j2 = j1 - bandwidth; 
                  if (j2 < 0) j2 = 0;
                  //for ( j = 0; j <= j1; j++ ) {
                  for ( j = j2; j <= j1; j++ ) {
                     //bj = (float)bspline->SchmidtMat[j1u+j];
                     bj = SchmidtMat[j1u+j];
		     aibj = ai * bj;
		     jjj  = j * usedN1;
				
                     k2 = k1 - bandwidth; 
                     if (k2 < 0) k2 = 0;
	             //for ( k = 0; k <= k1 ; k++ ) {
	             for ( k = k2; k <= k1 ; k++ ) {
                        //ck = (float)bspline->SchmidtMat[k1u+k]; 
                        ck = SchmidtMat[k1u+k]; 
		        aibjck = aibj*ck;
	                ii =  iii + jjj + k;
                        //O_coef_Result[jj]   = O_coef_Result[jj]   - aibjck*(tau2 *xdgdphi[ii]);        //for J1
                        // O_coef_Result[jj]   = O_coef_Result[jj]   - aibjck*xdgdphi[ii];   // Xu changed 0912     //for J1
                        dcoefs[jj]   = dcoefs[jj]   - aibjck*xdgdphi[ii];   // Xu changed 0912     //for J1
                        //Xdfphi += aibjck*xdgdphi[ii];
		     }
                  }
               }
               //printf("Xdfphi=%f ", Xdfphi);
               //O_coef_Result[jj]  = O_coef_Result[jj] - xdgdphi[jj];   
               // dcoefs[jj] = O_coef_Result[jj] - coef[jj]; // Xu comment out 0912  // difference	
            }
         }
      }

      t1_end=clock();
      time1_end=(double)t1_end/(double)CLOCKS_PER_SEC;
      printf("End update coefficients using J1.                 time=%f\n", time1_end-time1);
   }  // end of J1


   // compute J2--J5
   if (total != 0) { 
      // Loop for basis 
      for ( i1 = 0; i1 <= usedN; i1++ ) {
         i1u = i1*usedN1;
         for ( j1 = 0; j1 <= usedN; j1++ ) {
            j1u = j1*usedN1;
            for ( k1 = 0; k1 <= usedN; k1++ ) {
               
               k1u = k1*usedN1;
               jj  = i1*usedN2+j1*(usedN+1) + k1;
	       //O_coef_Result[jj] = coef[jj];  // Xu 0912
 

               // Orthogonalize
               i2 = i1 - bandwidth; 
               if (i2 < 0) i2 = 0;
               //for ( i = 0; i <= i1; i++ ) {
               for ( i = i2; i <= i1; i++ ) {
	          //ai  = (float)bspline->SchmidtMat[i1u+i] ;
	          ai  = SchmidtMat[i1u+i] ;
	          iii = i*usedN2;

                  j2 = j1 - bandwidth; 
                  if (j2 < 0) j2 = 0;
                  //for ( j = 0; j <= j1; j++ ) {
                  for ( j = j2; j <= j1; j++ ) {
                     //bj = (float)bspline->SchmidtMat[j1u+j];
                     bj = SchmidtMat[j1u+j];
		     aibj = ai * bj;
		     jjj  = j * usedN1;
				
                     k2 = k1 - bandwidth; 
                     if (k2 < 0) k2 = 0;
	             //for ( k = 0; k <= k1 ; k++ ) {
	             for ( k = k2; k <= k1 ; k++ ) {
                        //ck = (float)bspline->SchmidtMat[k1u+k]; 
                        ck = SchmidtMat[k1u+k]; 
		        aibjck = aibj*ck;
	                ii =  iii + jjj + k;
                        // O_coef_Result[jj]   = O_coef_Result[jj]   - aibjck*(J234phi[ii] +  J5phi[ii]); // Xu 0912 //for J2--J5
                        dcoefs[jj]   = dcoefs[jj]   - aibjck*(J234phi[ii] +  J5phi[ii]);   // Xu changed 0912     //for J1
		     }
                  }
               }

//              O_coef_Result[jj]   = O_coef_Result[jj]   - (J234phi[jj] +  J5phi[jj]);
            }
         }
      }
   }      // end of J2--J5


   
   //Compute temporal step-size.
   //comp_tau = 1.0;  // Xu 0912
   comp_tau = tau;  // Xu 0924
   
   if (reconj1 != 0.0 || taube != 0.0) { 
      t1=clock();
      time1=(double)t1/(double)CLOCKS_PER_SEC;

      // dcoefs is orthogonal difference, O_coef_Result is non-orto, tempararaly
      ConvertCoeffsToNewCoeffs(dcoefs, O_coef_Result, bspline->SchmidtMat);//convert to non-ortho bspline basis coefficients.
      t1_end=clock();
      time1_end=(double)t1_end/(double)CLOCKS_PER_SEC;
      printf("End ConvertCoeffs.                                time=%f\n", time1_end-time1);

      //  compute tau
      t1=clock();
      time1=(double)t1/(double)CLOCKS_PER_SEC;
      // Non_o_coef, O_coef_Result are non-orthogonal coeffients, currently
      if (newnv == nv ) {
         ComputeTau(Non_o_coef, O_coef_Result, &comp_tau,reconj1, taube, flows);    // Xu added 
      }  else {
         ComputeTau(newnv, index, Non_o_coef, O_coef_Result, &comp_tau, reconj1, taube, flows);   // Xu added
      }
       
      t1_end=clock();
      time1_end=(double)t1_end/(double)CLOCKS_PER_SEC;

      // Xu added 
      sum = 0.0;
      sum1 = 0.0;
      for ( ui = 0; ui < usedtolbsp; ui++ ) {
          sum = sum + O_coef_Result[ui]*O_coef_Result[ui];
          sum1 = sum1 + coef[ui]*coef[ui];
      }
      sum = sum/usedtolbsp;
      sum = sqrt(sum);
      printf("End computed tau====================%f      time=%f, residual = %e, function = %f, %e\n", 
                                                    comp_tau, time1_end-time1, sum, sum1, sum/sum1);
   }

   //Get updated ortho Bspline coeffs.      
   //comp_tau = comp_tau*0.618;
   // In the following, O_coef_Result, coef and dcoefs are orthogonal coeffients
   for ( i = 0; i < usedtolbsp; i++ ) {
      O_coef_Result[i] = coef[i] + comp_tau * dcoefs[i];
      coef[i] = O_coef_Result[i];
   } 	 

   cubes.clear();
   cubes1.clear();
   outliers.clear();
 
}  // End iteration    
 
// save variables for next iters.
// O_Kcoef and  coef are ortho Bspline coeffs
for ( i = 0; i <usedtolbsp; i++ )  {
   O_Kcoef[i] = coef[i];   // is O_Kcoef orth for the first time??  // Xu 0912
   //printf("\ncoef=%f ", coef[i]);
}
 
t_end=clock();
time_end=(double)t_end/(double)CLOCKS_PER_SEC;
printf("\n++++ t = %f, time = %f, CLOCKS_PER_SEC = %e\n", t_end, time, (double)CLOCKS_PER_SEC);
time=time_end-time;
Iters = Iters + iter;
 printf("\nReconstruction Time = %f \n",(double)time);
printf("\nToltal iter number = %d \n", Iters);



//get reconstructed volume data.
Oimage* Object = InitImageParameters(1, ImgNx+1, ImgNy+1,ImgNz+1, 1); // 6.6 ImgNx, ImgNy, ImgNz changed to ImgNx+1, ....
 
Object->nx = Object->ny = Object->nz = ImgNx+1;   //6.6 ImgNx to ImgNx +1

//convert to non-ortho bspline basis coefficients.
ConvertCoeffsToNewCoeffs(coef, Non_o_coef, bspline->SchmidtMat);

// Xu added this for test J1, 0918 
//ComputeEJ1_Acurate(f, Non_o_coef);

//for ( i = 0; i < usedtolbsp; i++ ) {
//   if (Non_o_coef[i] < 0.0) Non_o_coef[i] = 0.0;   // Xu Cut small coefficients 0912
//} 	 

if(Bscale == 1.0 ) bspline->ObtainObjectFromNonOrthoCoeffs(Non_o_coef,Object->data);
if(Bscale == 2.0)  bspline->ObtainObjectFromNonOrthoCoeffs_FA(Non_o_coef,Object->data);

maxf = Object->data[0];
minf = Object->data[0];
for ( i = 0; i <= ImgNx; i++ ) {
   for ( j = 0; j <= ImgNx; j++ ) {
      for ( k = 0; k <= ImgNx; k++ ) {
         i1 = i*img2dsize + j * (ImgNx+1) + k;
	 if (maxf < Object->data[i1])  maxf = Object->data[i1];
         if (minf > Object->data[i1])  minf = Object->data[i1];
      }  
   }
}	 
printf("\nmaxf========= %f minf=========%f  ", maxf, minf);
 

/*
 ttt = 113;
 cut = minf + (maxf - minf)*ttt/255;
 printf(" Cut ==== %f\n", cut);
 //for ( i = 0; i <usedtolbsp; i++ )  {
 //   if (O_Kcoef[i] < cut) O_Kcoef[i] = minf;
 //}
 */

//free(subf);   subf  = NULL;
//free(padf);   padf  = NULL;
 free(J234phi); J234phi = NULL;
// free(J3phi); J3phi = NULL;
// free(J4phi); J4phi = NULL;
 free(J5phi); J5phi = NULL;
//fftw_free(in); fftw_free(out); 
//in = NULL;     out = NULL;	 
free(dcoefs);
free(index);
free(subset);

free(SchmidtMat);
return Object;

}


float Reconstruction::ConvertCoeffsToNewCoeffs(float *coef, float *Non_o_coef, long double *TransformMat)
{
  int i, j, k, s, id,  iN2, jN1, usedN2;
  float *F, *F0, *A, *AT, Fjk;
  float *ai;
  int    N2, usedn;
  N2 = (N+1)*(N+1);
 
  usedN2 = (usedN + 1)*(usedN + 1);

  usedn = usedN;

  ai = (float *)malloc((N-3) * sizeof(float));

  F  = (float *)malloc(usedN2*sizeof(float));
  F0 = (float *)malloc(usedN2*sizeof(float));
  A  = (float *)malloc(usedN2*sizeof(float));
  AT = (float *)malloc(usedN2*sizeof(float));

  for ( i = 0; i < usedN + 1; i++ ) ai[i] = 0.0;

  for ( i = 0; i < (usedN+1)*(usedN+1)*(usedN+1); i++ )
	Non_o_coef[i] = 0.0;

  for ( i = 0; i < usedN2; i++ )
	{
	  F[i]  = 0.0;
	  F0[i] = 0.0;
	  A[i]  = 0.0;
	  AT[i] = 0.0;
	 
	}


  for ( i = 0; i < usedN + 1; i++ )
      {
        jN1 = i * (usedN+1);
   	for ( j = 0; j < usedN + 1; j++ )
	  {
	   A[jN1 + j] = (float)TransformMat[jN1 + j];
           AT[jN1 + j] = A[jN1 + j];
	  }
       }

  MatrixTranspose(AT,usedn+1,usedn+1);

  for ( i = 0; i < usedN+1; i++ )
     {
         iN2 = i * usedN2;

	  for( j = 0; j < usedN + 1; j++ )
              {
                jN1 = j * (usedN+1);
		for ( k = 0; k < usedN + 1; k++ )
		  F0[jN1+k] = coef[iN2 + jN1 + k];
              } 
      MatrixMultiply(F0, usedn+1, usedn+1, A, usedn+1, usedn+1, F);
      MatrixMultiply(AT, usedn+1, usedn+1, F, usedn+1, usedn+1, F0);

      for ( j = 0; j < usedN+1; j++)
		ai[j] = A[i*(usedN+1)+j];

      for ( s = 0; s <= usedN; s++ )
          {
           iN2 = s * usedN2;

       	   for ( j = 0; j <= usedN; j++ )
	       {
                jN1 = j * (usedN+1);

  	        for ( k = 0; k <= usedN; k++ )
	 	   Non_o_coef[iN2 + jN1 + k] =Non_o_coef[iN2 + jN1 + k]+ ai[s]*F0[jN1 + k];
               }
           }

     }   

  free(ai); ai = NULL;
  free(F);F = NULL;
  free(F0); F0 = NULL;
  free(A); A = NULL;
  free(AT); AT = NULL;

  /*
	  for ( s = 0; s <= N; s++ )
		for ( j = 0; j <= N; j++ )
		  for ( k = 0; k <= N; k++ )
			printf("\ncoef=%f        ", Non_o_coef[s*N2 + j*(N+1) + k]);
getchar();

  */
}




/************************************************************************
 Descriptions:
     Convert the one dimension index to ijk index.

 ************************************************************************/
void Reconstruction:: Converttoijk(vector<int> &suppi,int *ijk)
{
  int i,n,NN2, NN1,remainder;

  //NN2 = (ImgNx+1) * (ImgNx+1); //(N-1)*(N-1);
  NN1 = (ImgNx+1);         //(N-1);

  n = suppi.size();

  for ( i = 0; i < n; i++ )
	{
	  ijk[3*i+0] = suppi[i]/img2dsize;
	  remainder  = suppi[i]%img2dsize;
    ijk[3*i+1] = remainder/NN1;
    ijk[3*i+2] = remainder%NN1;

//	printf("\ni j k = %d %d %d ", ijk[3*i+0], ijk[3*i+1], ijk[3*i+2]);
	}
  
}


void Reconstruction::GetSubNarrowBand(vector<CUBE> &cubes, vector<CUBE> &cubes1, 
                                        float *IsoC, float *coefs, float *f, float param) //, vector<float> &DeltafIsoC,  vector<float> &PdeltafIsoC)
{
  int s, i,j, k, ii, x, y, z, p, q, r, scale, lx, rx, ly, ry, lz, rz, id, jd, kd, id4, jd2, usedN2, ii8;
  float X, Y, Z, del, phi_p, phi_q, phi_r, phi_pq, value, rscale, sum, sum8;
  float cubef[8], subcubef[27], df[3],s8[8], palpha, ox, oy, oz, length, length25;
  int x9, x19, y3, y13, z1, index[8];
  CUBE  onecube;
  scale  = (int)Bscale;
  rscale = 1.0/Bscale;
  //del    = 1.0/factor;
  usedN2 = (usedN+1)*(usedN+1);  
  palpha = param * alpha;

/*
  //for test.
  int subimg2dsize, subImgNx;
  subImgNx = factor * ImgNx;
  subimg2dsize = (subImgNx+1) * (subImgNx+1);
// end test.
*/

//using stack.


while(!cubes.empty()) {
/*
//old methods.
while(!ijk.at(0).empty())
{

 i = ijk.at(0).back();
 j = ijk.at(1).back();
 k = ijk.at(2).back();

 ijk.at(0).pop_back();
 ijk.at(1).pop_back();
 ijk.at(2).pop_back();
*/

   onecube = cubes.back();  //Acess last element.
   cubes.pop_back();        // Delete last element.

   ox = onecube.ox;
   oy = onecube.oy;
   oz = onecube.oz;
   length = onecube.l;
   length25 = 0.25*length;

   cubef[0] = onecube.cubef8[0]; //f[i  *  img2dsize+j    *(ImgNx+1)+k  ];
   cubef[1] = onecube.cubef8[1]; //f[i  *  img2dsize+j    *(ImgNx+1)+k+1];
   cubef[2] = onecube.cubef8[2]; //f[i  *  img2dsize+(j+1)*(ImgNx+1)+k  ];
   cubef[3] = onecube.cubef8[3]; //f[i  *  img2dsize+(j+1)*(ImgNx+1)+k+1];
   cubef[4] = onecube.cubef8[4]; //f[(i+1)*img2dsize+j    *(ImgNx+1)+k  ];
   cubef[5] = onecube.cubef8[5]; //f[(i+1)*img2dsize+j    *(ImgNx+1)+k+1];
   cubef[6] = onecube.cubef8[6]; //f[(i+1)*img2dsize+(j+1)*(ImgNx+1)+k  ];
   cubef[7] = onecube.cubef8[7]; //f[(i+1)*img2dsize+(j+1)*(ImgNx+1)+k+1];

   sum8 = 0.125*(cubef[0]+cubef[1]+cubef[2]+cubef[3]+cubef[4]+cubef[5]+cubef[6]+cubef[7]);
   // printf("\n%f %f %f %f %f %f %f %f", cubef[0], 	cubef[1],cubef[2], cubef[3], cubef[4],cubef[5],cubef[6],cubef[7]);  

   //split the current cube into eight sub-cubes.
   //compute 27 sub cube vertex values.
   ii = 0;
   for ( x = 0; x <= 2; x++ ) {
      //id    = x * 0.5;
      //id4   = 4 * id;
      //df[0] = x * 0.5;

      for ( y = 0; y <= 2; y++ ) {
         //jd    = y * 0.5;
         //jd2   = 2 * jd;
         //df[1] = y * 0.5;

         for ( z = 0; z <= 2; z++, ii++) {
            //kd    = z * 0.5;
            //df[2] = z * 0.5;
             //trilinear interpolation.         
/*
             subcubef[ii] = cubef[0] * (1-df[2])*(1-df[1])*(1-df[0])+
                            cubef[1] * df[2]    *(1-df[1])*(1-df[0])+
                            cubef[2] * (1-df[2])*df[1]    *(1-df[0])+
                            cubef[3] * df[2]    *df[1]    *(1-df[0])+  
                            cubef[4] * (1-df[2])*(1-df[1])*df[0]    +
                            cubef[5] * df[2]    *(1-df[1])*df[0]    +
                            cubef[6] * (1-df[2])*df[1]    *df[0]    +
                            cubef[7] * df[2]    *df[1]    *df[0];
*/

             // This is the fastest way among the four. 
             ii8 = ii*8;
             subcubef[ii] = cubef[0] * Cube_Coeff[ii8    ] +
                            cubef[1] * Cube_Coeff[ii8 + 1] +
                            cubef[2] * Cube_Coeff[ii8 + 2] +
                            cubef[3] * Cube_Coeff[ii8 + 3] +
                            cubef[4] * Cube_Coeff[ii8 + 4] +
                            cubef[5] * Cube_Coeff[ii8 + 5] +
                            cubef[6] * Cube_Coeff[ii8 + 6] +
                            cubef[7] * Cube_Coeff[ii8 + 7];


/*
             if((x==0||x==2) && (y==0||y==2) && (z==0||z==2))
                subcubef[ii] = cubef[id4+jd2+kd];
             else if((x==1) && (y==0||y==2) && (z==0||z==2)) 
                subcubef[ii] = 0.5*(cubef[jd2+kd]+ cubef[4+jd2+kd]);  
             else if((y==1) && (x==0||x==2) && (z==0||z==2))
                subcubef[ii] = 0.5*(cubef[id4+kd]+ cubef[id4+2+kd]);
             else if((z==1) && (x==0||x==2) && (y==0||y==2))
                subcubef[ii] = 0.5*(cubef[id4+jd2]+ cubef[id4+jd2+1]);
             else if(x==1 && y==1 && z==1)
                subcubef[ii] = sum8;
             else if(x==0 || x==2)
                subcubef[ii] = 0.25*(cubef[id4]+cubef[id4+1]+cubef[id4+2]+cubef[id4+3]); 
             else if(y==0 || y==2)
                subcubef[ii] = 0.25*(cubef[jd2]+cubef[jd2+1]+cubef[4+jd2]+cubef[4+jd2+1]);    
             else
                subcubef[ii] = 0.25*(cubef[kd]+cubef[2+kd]+cubef[4+kd]+cubef[4+2+kd]);    
*/             
             //printf("\ni j k = %d %d %d subcubef=%f ", x, y, z, subcubef[ii]);
         }
      }
   }

   for ( x = 0; x < 2; x++ ) {
      x9 = x*9;
      x19 = x9 + 9;
      for ( y = 0; y < 2; y++ ) {
         y3 = y*3;
         y13 = y3 + 3;
         for ( z = 0; z < 2; z++ ) {
            z1 = z + 1;
            s8[0] = subcubef[x9+y3+z  ];
            s8[1] = subcubef[x9+y3+z1];
            s8[2] = subcubef[x9+y13+z  ];
            s8[3] = subcubef[x9+y13+z1];
            s8[4] = subcubef[x19+y3+z  ];
            s8[5] = subcubef[x19+y3+z1];
            s8[6] = subcubef[x19+y13+z  ];
            s8[7] = subcubef[x19+y13+z1];
            
            // Xu changed
            // SortDecreasing(s8, 8);
            Sortingfloat(s8, index, 8);
            //printf("\ns8-s1=%f ", s8[0]-s8[7]);

            // Xu changed
            //if(s8[0]-*IsoC < -alpha || s8[7]-*IsoC > alpha ) continue;
            if(s8[index[7]]-*IsoC < -alpha || s8[index[0]]-*IsoC > alpha ) continue;
            //printf("\n%f %f %f %f %f %f %f %f ", s8[0], s8[1], s8[2], s8[3],s8[4], s8[5], s8[6], s8[7]);
           
            onecube.ox = ox + (-1+2*x)*length25;
            onecube.oy = oy + (-1+2*y)*length25;
            onecube.oz = oz + (-1+2*z)*length25;
            onecube.l  = length*0.5;

            /*
            // Xu comment these out 
            onecube.cubef8[0] = s8[0];
            onecube.cubef8[1] = s8[1];
            onecube.cubef8[2] = s8[2];
            onecube.cubef8[3] = s8[3];
            onecube.cubef8[4] = s8[4];
            onecube.cubef8[5] = s8[5];
            onecube.cubef8[6] = s8[6];
            onecube.cubef8[7] = s8[7];
            */

            // Xu changed
            //if (s8[0] - s8[7] < palpha) {
            if (s8[index[7]] - s8[index[0]] < palpha) {
                //printf("\ns8-s1=%f ", s8[0]-s8[7]);
                // Xu asked. Are these following 5 lines useful? 
                //sum  = 0.125*(s8[0]+s8[1]+s8[2]+s8[3]+s8[4]+s8[5]+s8[6]+s8[7]);
                //onecube.cubef8[0] = sum;                 // the first element saves the middle point value.
                //sum  = sum - *IsoC;
                //DeltafIsoC.push_back(Deltafunc(sum));
                //PdeltafIsoC.push_back(DeltafuncPartials(sum));
     
                onecube.l = onecube.l * onecube.l * onecube.l;
                cubes1.push_back(onecube); 
                //printf("\nonecube.oxyz=%f %f %f l=%f ", onecube.ox, onecube.oy, onecube.oz, onecube.l);
            } else { 

               // Xu put these here
               onecube.cubef8[0] = s8[0];
               onecube.cubef8[1] = s8[1];
               onecube.cubef8[2] = s8[2];
               onecube.cubef8[3] = s8[3];
               onecube.cubef8[4] = s8[4];
               onecube.cubef8[5] = s8[5];
               onecube.cubef8[6] = s8[6];
               onecube.cubef8[7] = s8[7];

               cubes.push_back(onecube);
            }
         }
      } 
   }
}
 
//old methods.
/*
  for ( s = 0; s < ijk.at(0).size(); s++ )
     {
      i = ijk.at(0).at(s);
      j = ijk.at(1).at(s);
      k = ijk.at(2).at(s);

      lx = i/scale - 3;
      rx = lx + 3;
      ly = j/scale - 3;
      ry = ly + 3;
      lz = k/scale - 3;
      rz = lz + 3;
 
      for ( x = 0; x <= factor; x++ )
        {
         X = (i + x * del)*rscale;
         for ( y = 0; y <= factor; y++ )
            {
             Y = (j + y * del)*rscale;
             for ( z = 0; z <= factor; z++ )
               {
                Z = (k + z * del)*rscale;
                sum = 0.0;
                //printf("\nX Y Z = %f %f %f ", X, Y, Z);
              
                for ( p = lx; p <= rx; p++ ) 
                  {
                   if( p < 0 || p > usedN ) continue;
                   bspline->Spline_N_Base(X-p-2.0, &value);
                   phi_p = value;
                   id    = p * usedN2;
                   for ( q = ly; q <= ry; q++ )
                      {
                       if( q < 0 || q > usedN ) continue;
                       bspline->Spline_N_Base(Y-q-2.0, &value);
                       phi_q  = value;
                       phi_pq = phi_p * phi_q;
                       jd     = q * (usedN+1);

                       for ( r = lz; r <= rz; r++ )
                          {
                           if( r < 0 || r > usedN ) continue;
                           bspline->Spline_N_Base(Z-r-2.0, &value);
                           phi_r = value;
                           ii    = id + jd + r;
                           sum  =sum  +  coefs[ii] * phi_pq * phi_r;
                           //printf("\nii = %d B=%1.15e coef=%1.15e ji=%1.18e", ii, phi_pq*phi_r, coefs[ii],
                                     //coefs[ii] * phi_pq * phi_r);

                          }
                       }
                    }
                //if ((x == 0 || x == factor) && (y == 0 || y == factor) && (z == 0 || z == factor) ) 
                   //printf("\nsum-f=%e ", sum-f[(i+x/factor)*img2dsize+(j+y/factor)*(ImgNx+1)+k+z/factor]);
                //printf("\nsum-subf=%e ", sum-f[(factor*i+x)*subimg2dsize+(factor*j+y)*(subImgNx+1)+factor*k+z]);
                //printf("\nlx rx ly ry lz rz = %d %d %d %d %d %d ", lx, rx, ly, ry, lz, rz);
                 //printf("sum = %1.15e ", sum);
                //printf("\nX Y Z = %f %f %f ", X, Y, Z);
                sum = sum - *IsoC;
                if(fabs(sum) > alpha ) continue;
                subijk.at(0).push_back(factor*i+x);
                subijk.at(1).push_back(factor*j+y);
                subijk.at(2).push_back(factor*k+z);
                //printf("\nx y z = %d %d %d ", factor*i+x, factor*j+y, factor*k+z);
                //printf("\n ix iy iz = %d %d %d new ix iy iz=%d %d %d ", i, j, k, (factor*i+x)/factor, (factor*j+y)/factor, (factor*k+z)/factor);

                DeltafIsoC.push_back(Deltafunc(sum));
                PdeltafIsoC.push_back(DeltafuncPartials(sum));


                }
              }
           }
 
     }
*/
}

/************************************************************************
 Descriptions:

 ************************************************************************/
void Reconstruction::ComputeDeltafIsoC(float *f, float *IsoC, int n, 
                                       vector<vector<int> > &ijk, float *DeltafIsoC, float *PdeltafIsoC)
{
  int i,ii, ix,iy,iz;
  //int N2;

  //N2 = (ImgNx+1) * (ImgNx+1);

  for ( i = 0; i < n; i++ )
	{
/*
	  ix = ijk[3*i + 0];
          iy = ijk[3*i + 1];
	  iz = ijk[3*i + 2];
*/
          ix = ijk.at(0).at(i);
          iy = ijk.at(1).at(i);
          iz = ijk.at(2).at(i);

	  ii = ix * img2dsize + iy * (ImgNx+1) + iz;

      DeltafIsoC[i] = Deltafunc(f[ii]-*IsoC);
      PdeltafIsoC[i]= DeltafuncPartials(f[ii]-*IsoC);
	  //printf("\nDeltafIsoC=%f ", DeltafIsoC[i]);
	}

}





/************************************************************************
 Descriptions:
     Compute the Hessian Matrix of Volume data f in the narrow band.

*************************************************************************/
void Reconstruction::ComputeHfGradientf(float *coefs, float *IsoC, int n, vector<vector<int> > &subijk, 
            vector<CUBE> &cubes1, int subfactor, float *DeltafIsoC, float *PdeltafIsoC, float *Hf, 
            float *Gradientf, float *Normgradf, float *HGf, float *GHGf,  float *XtN, float *NNt, float *XtP, float *f)
{
  int   id, id9, i, j, k, ii, id3, lix, rix, liy, riy, liz, riz, scale, half;
  int   ix, iy, iz, N2;
  float partials[20], temp, epsilon, X[3], del, rscale, values[3], value, normal[3], m[9], ngradf, rngradf;
  float Bx, Bx1, Bx2, By[10], By1[10], By2[10], Bz[10], Bz1[10], Bz2[10];
  float Bxy,Bxy1, Bx1y, Bx1y1, Bx2y, Bxy2, Imatrix[9], sum;  
  int   id91, id92, id93, id94, id95, id96, id97, id98, id31, id32, jly, klz;
  float rscale2;
  CUBE  onecube;

  epsilon = 0.000001;  //If Norm gradient is zero.
  scale   = (int)Bscale;
  rscale  = 1.0/Bscale;
  rscale2 = rscale*rscale;   // Xu Added
  N2      = (usedN+1)*(usedN+1);
  del     = 1.0/subfactor * rscale;
  half    = ImgNx/2;

  for ( i = 0; i < 9; i++ ) Imatrix[i] = 0.0;
  Imatrix[0] = Imatrix[4] = Imatrix[8] = 1.0;


/*
//for test.
int subimg2dsize , subImgNx, tx, ty, tz;
subImgNx     = subfactor*ImgNx;
subimg2dsize = (subImgNx+1)*(subImgNx+1);
*/


// Loop for all the cubes in the narrow band
for ( id = 0; id < n; id++ ) {
   id3  = 3 * id;
   id31 = id3 + 1;
   id32 = id3 + 2;

   id9  = 9 * id;
   id91 = id9 + 1;
   id92 = id9 + 2;
   id93 = id9 + 3;
   id94 = id9 + 4;
   id95 = id9 + 5;
   id96 = id9 + 6;
   id97 = id9 + 7;
   id98 = id9 + 8;
 
   onecube = cubes1.at(id);
   //printf("\nonecube.oxyz=%f %f %f l=%f ", onecube.ox, onecube.oy, onecube.oz, onecube.l);

  
 /*
   ix = (int)onecube.ox;
   iy = (int)onecube.oy;
   iz = (int)onecube.oz;
   */
   // Xu changed
   ix = floor(onecube.ox);
   iy = floor(onecube.oy);
   iz = floor(onecube.oz);

   lix = ix/scale-3;  
   rix = lix+3;           
   if (lix < 0) lix = 0;
   if (rix > usedN) rix = usedN;

   liy = iy/scale-3;
   riy = liy+3;          
   if (liy < 0) liy = 0;
   if (riy > usedN) riy = usedN;

   liz = iz/scale-3;
   riz = liz+3;         
   if (liz < 0) liz = 0;
   if (riz > usedN) riz = usedN;

   if(lix > usedN || liy > usedN || liz > usedN ) {printf("\nerror. Lix liy liz is bigger than usedN.line 1235.\n ");getchar();}
   X[0] = onecube.ox*rscale - 2.0;
   X[1] = onecube.oy*rscale - 2.0;
   X[2] = onecube.oz*rscale - 2.0;
 
   // Xu added the following 14 lines
   for ( j = liy; j <= riy; j++ ) {
      jly = j - liy;
      bspline->Spline_N_Base_2(X[1]-j, values);
      By[jly]   = values[0];
      By1[jly]  = values[1]*rscale;
      By2[jly]  = values[2]*rscale2;
   }

   for ( k = liz; k <= riz; k++ ) {
      klz = k - liz; 
      bspline->Spline_N_Base_2(X[2]-k, values);
      Bz[klz]  = values[0];
      Bz1[klz] = values[1]*rscale;
      Bz2[klz] = values[2]*rscale2;
   }

//old methods.
/*
         ix = subijk.at(0).at(id);
         iy = subijk.at(1).at(id);
         iz = subijk.at(2).at(id);

         ix = ix/subfactor;
         iy = iy/subfactor;
         iz = iz/subfactor;
         //printf("\n ix iy iz = %d %d %d ", ix, iy, iz);

        ////////tx = subijk.at(0).at(id);
        ////////ty = subijk.at(1).at(id);
        ///////tz = subijk.at(2).at(id);
        //printf("\nx y z = %d %d %d ", tx, ty, tz); 

        lix = ix/scale-3;  
        rix = lix+3;           

        liy = iy/scale-3;
        riy = liy+3;          

        liz = iz/scale-3;
        riz = liz+3;         

        //printf("\nlx rx ly ry lz rz = %d %d %d %d %d %d ", lix, rix, liy, riy, liz, riz);
        X = subijk.at(0).at(id)*del;     
        Y = subijk.at(1).at(id)*del;    
        Z = subijk.at(2).at(id)*del;
        //printf("\nX Y Z = %f %f %f ", X, Y, Z);
*/  
   sum = 0.0;
   // Loop for bsisis functions 
   for ( i = lix; i<= rix; i++ ) {
      //if( i < 0 || i > usedN ) continue;    // Xu comment out
      bspline->Spline_N_Base_2(X[0]-i, values);
      /*
      Bx  = values[0];
      Bx1 = values[1];
      Bx2 = values[2];  
      */

      // Xu Changed
      Bx  = values[0];
      Bx1 = values[1]*rscale;
      Bx2 = values[2]*rscale2;  
      for ( j = liy; j <= riy; j++ ) {
         //if ( j < 0 || j > usedN ) continue;  // Xu comment out 
         //bspline->Spline_N_Base_3(X[1]-j, values);
         /*
         By  = values[0];
         By1 = values[1];
         By2 = values[2];
         */

         // Xu changed
         //By  = values[0];
         //By1 = values[1]*rscale;
         //By2 = values[2]*rscale2;

         jly = j - liy;
         Bxy  = Bx  * By[jly];
         Bx1y = Bx1 * By[jly];
         Bxy1 = Bx  * By1[jly];
       
         Bx2y = Bx2 * By[jly];
         Bx1y1= Bx1 * By1[jly];
         //Bx1y = Bx1 * By[jly];
         Bxy2 = Bx  * By2[jly];

         for ( k = liz; k <= riz; k++ ) {
            //if ( k < 0 || k > usedN ) continue;   // Xu comment out 
            //bspline->Spline_N_Base_3(X[2]-k, values);
            /*
            Bz  = values[0];
            Bz1 = values[1];
            Bz2 = values[2];
            */

            // Xu changed 
            //Bz  = values[0];
            //Bz1 = values[1]*rscale;
            //Bz2 = values[2]*rscale2;


            klz = k - liz;  
            ii = i*N2 + j*(usedN+1) + k; // 
	    sum = sum + coefs[ii] * Bxy * Bz[klz];
                      //printf("\nii = %d B=%1.15e coef=%1.15e ji=%1.18e", ii, Bxy*Bz, coefs[ii]
                        //      ,coefs[ii] * Bxy * Bz);	

            // compute gradient
	    Gradientf[id3 ] = Gradientf[id3 ] + coefs[ii]*Bx1y*Bz[klz];
	    Gradientf[id31] = Gradientf[id31] + coefs[ii]*Bxy1*Bz[klz];
	    Gradientf[id32] = Gradientf[id32] + coefs[ii]*Bxy *Bz1[klz];
			
			
            // compute Hessian 
	    Hf[id9 ] = Hf[id9 ] + coefs[ii]*Bx2y*Bz[klz]; //partials[4];
	    Hf[id91] = Hf[id91] + coefs[ii]*Bx1y1*Bz[klz]; //partials[5];
	    Hf[id92] = Hf[id92] + coefs[ii]*Bx1y*Bz1[klz]; //partials[6];
	    //Hf[id93] = Hf[id93] + coefs[ii]*Bx1y1*Bz[klz]; //partials[5];  // Xu comment out because symmetricty
	    Hf[id94] = Hf[id94] + coefs[ii]*Bxy2*Bz[klz];  //partials[7];
	    Hf[id95] = Hf[id95] + coefs[ii]*Bxy1*Bz1[klz]; //partials[8];
	    //Hf[id96] = Hf[id96] + coefs[ii]*Bx1y*Bz1[klz]; //partials[6];
	    //Hf[id97] = Hf[id97] + coefs[ii]*Bxy1*Bz1[klz]; //partials[8];
            Hf[id98] = Hf[id98] + coefs[ii]*Bxy*Bz2[klz]; //partials[9];
	 }		
      }
   }

   // Xu added 
   Hf[id93] = Hf[id91];
   Hf[id96] = Hf[id92];
   Hf[id97] = Hf[id95];

   // compute the length of gradient
   temp = Gradientf[id3]*Gradientf[id3] + Gradientf[id31]*Gradientf[id31] + Gradientf[id32]* Gradientf[id32];	
   Normgradf[id]    = sqrt(temp);
  
   DeltafIsoC[id] = Deltafunc(sum - *IsoC);
   PdeltafIsoC[id]= DeltafuncPartials(sum - *IsoC);

 
        //if(fabs(Gradientf[id3+0]-f[tx*subimg2dsize+ty*(subImgNx+1)+tz]) != 0.0 ) 
        //if(fabs(Hf[id9+0]-f[tx*subimg2dsize+ty*(subImgNx+1)+tz]) != 0.0 ) 
        //printf("\ngradient-subf=%e ", Gradientf[id3+0]-f[tx*subimg2dsize+ty*(subImgNx+1)+tz]);
        //printf("\nHf-subf=%e ", Hf[id9+0]-f[tx*subimg2dsize+ty*(subImgNx+1)+tz]);
        //printf("sum = %1.15e ", sum);
   // compute H*G
   for ( j = 0; j < 3; j++ ) {
      // HGf[id3 + j] = 0.0; // Xu have set to zero before

      for ( k = 0; k < 3; k++ ) {
         HGf[id3+j] = HGf[id3+j] + Hf[id9 + 3*j + k] * Gradientf[id3 + k];
      }
   }

   // compute G^T*H*G
   for ( j = 0; j < 3; j++ ) {
       GHGf[id] = GHGf[id] + Gradientf[id3 + j] * HGf[id3 + j];
   }

   normal[0] = Gradientf[id3];
   normal[1] = Gradientf[id31];
   normal[2] = Gradientf[id32];


   ngradf = Normgradf[id];
   // Xu changed
   if( ngradf <= epsilon ) ngradf= epsilon;
   rngradf   = 1.0/ngradf;

   // surface normal N
   normal[0] = normal[0]*rngradf;
   normal[1] = normal[1]*rngradf;
   normal[2] = normal[2]*rngradf;
  
   X[0] = onecube.ox - half;
   X[1] = onecube.oy - half;
   X[2] = onecube.oz - half;

   // Compute X^T* N
   XtN[id] = X[0]*normal[0]+X[1]*normal[1]+X[2]*normal[2];

   // Compute N*N^T
   MatrixMultiply(normal,3,1,normal,1,3,m);
   for ( j = 0; j < 9; j++ ) {
       NNt[id9+j] = m[j]; 
   }
   
   // Compute P = I - N*N^T
   for ( j = 0; j < 9; j++ ) {
       m[j] = Imatrix[j] - m[j];
   }

   // Compute X^T*P
   MatrixMultiply(X, 1, 3, m, 3, 3, values);
   for ( j = 0; j < 3; j++ ) {
      XtP[id3+j] = values[j];
   }



   // Testing the data 
   float tem, Div, tem2;
   /*
   // Test function values, The test result shows that it is OK
   tem = 36.0 - (X[0]*X[0]+X[1]*X[1]+X[2]*X[2]);
   printf("\n sum = %1.15e, X*X = %f, error = %e ", sum, tem, sum-tem);
   */

   /*
   // Test gradient, The test result shows that it is OK. 
   tem = sqrt(InnerProduct(Gradientf+id3,Gradientf+id3));
   tem2 = 2*sqrt(X[0]*X[0]+X[1]*X[1]+X[2]*X[2]);
   printf("\n ecact gradient = %f, computed = %f", tem2, tem);
   */

   // Test Hessian, the result shows it is not OK , the exact result is -2*I_3
   //printf("\n H = %e,%e,%e,%e,%e,%e,%e,%e,%e", Hf[id9],Hf[id91],Hf[id92],Hf[id93],Hf[id94],Hf[id95],Hf[id96],Hf[id97],Hf[id98]); 


   // test curvature, it is not ok. Because Hessian is not OK.
   tem = 1.0/sqrt(X[0]*X[0]+X[1]*X[1]+X[2]*X[2]);
   
   Div =  ((Hf[id9] + Hf[id9+4] + Hf[id9+8])*rngradf -
             - (GHGf[i])*rngradf*rngradf*rngradf);  // Xu changed
   //printf("\n  exact cur = %f, computed = %f", tem, -0.5*Div);
   //printf("\n H = %e,%e,%e,%e,%e,%e,%e,%e,%e", Hf[id9],Hf[id91],Hf[id92],Hf[id93],Hf[id94],Hf[id95],Hf[id96],Hf[id97],Hf[id98]); 


}  // end id loop
}


//                                                                              
void Reconstruction::ComputeG_HfGf(int n, vector<vector<int> > &ijk, float *coef, float *Hf, float *Gradientf, float *Hfx, float *Hfy, float *Hfz, float *G_HfGf)
{

  int   id, id9, i, j, k, ii, id3, lix, rix, liy, riy, liz, riz, scale;
  int   ix, iy, iz, N2;
  float partials[20], gf[3], gfx[3], gfy[3], gfz[3], hf[9], g_hfgf[9];
  float tmp1[9], tmp2[9];

  scale = (int)Bscale;

  N2 = (usedN+1)*(usedN+1);

  for ( id = 0; id < n; id++ )
	{
	  id3 = 3 * id;
	  id9 = 9 * id;

/*
	ix = ijk[id3 + 0];
        iy = ijk[id3 + 1];
	iz = ijk[id3 + 2];
*/

         ix = ijk.at(0).at(id);
         iy = ijk.at(1).at(id);
         iz = ijk.at(2).at(id);



/*
	for ( i = 0; i <= usedN; i++ )
	  for ( j = 0; j <= usedN; j++ )
		for ( k = 0; k <= usedN; k++ )
          {
*/


        lix = ix/scale-3;
        rix = ix/scale;

        liy = iy/scale-3;
        riy = iy/scale;

        liz = iz/scale-3;
        riz = iz/scale;

        for ( i = lix; i<= rix; i++ )
            {
              if( i < 0 || i > usedN ) continue;

              for ( j = liy; j <= riy; j++ )
                 {
                   if ( j < 0 || j > usedN ) continue;

                 for ( k = liz; k <= riz; k++ )
                    {
                       if ( k < 0 || k > usedN ) continue;


			ii = i*N2 + j*(usedN+1) + k;
			bspline->Phi_ijk_Partials_ImgGrid(ix, iy, iz, i, j, k, partials);
			
			Hfx[0] = Hfx[0] + coef[ii]*partials[10];
			Hfx[1] = Hfx[1] + coef[ii]*partials[11];
			Hfx[2] = Hfx[2] + coef[ii]*partials[12];
			Hfx[3] = Hfx[3] + coef[ii]*partials[11];
			Hfx[4] = Hfx[4] + coef[ii]*partials[13];
			Hfx[5] = Hfx[5] + coef[ii]*partials[14];
			Hfx[6] = Hfx[6] + coef[ii]*partials[12];
			Hfx[7] = Hfx[7] + coef[ii]*partials[14];
			Hfx[8] = Hfx[8] + coef[ii]*partials[15];


			Hfy[0] = Hfy[0] + coef[ii]*partials[11];
			Hfy[1] = Hfy[1] + coef[ii]*partials[13];
			Hfy[2] = Hfy[2] + coef[ii]*partials[14];
			Hfy[3] = Hfy[3] + coef[ii]*partials[13];
			Hfy[4] = Hfy[4] + coef[ii]*partials[16];
			Hfy[5] = Hfy[5] + coef[ii]*partials[17];
			Hfy[6] = Hfy[6] + coef[ii]*partials[14];
			Hfy[7] = Hfy[7] + coef[ii]*partials[17];
			Hfy[8] = Hfy[8] + coef[ii]*partials[18];
			



			Hfz[0] = Hfz[0] + coef[ii]*partials[12];
			Hfz[1] = Hfz[1] + coef[ii]*partials[14];
			Hfz[2] = Hfz[2] + coef[ii]*partials[15];
			Hfz[3] = Hfz[3] + coef[ii]*partials[14];
			Hfz[4] = Hfz[4] + coef[ii]*partials[17];
			Hfz[5] = Hfz[5] + coef[ii]*partials[18];
			Hfz[6] = Hfz[6] + coef[ii]*partials[15];
			Hfz[7] = Hfz[7] + coef[ii]*partials[18];
			Hfz[8] = Hfz[8] + coef[ii]*partials[19];
						
		     }
               }
         }

	gfx[0] = hf[0] = Hf[id9+0];     
	gfx[1] = hf[1] = Hf[id9+1];
	gfx[2] = hf[2] = Hf[id9+2];
	gfy[0] = hf[3] = Hf[id9+3];
	gfy[1] = hf[4] = Hf[id9+4];
	gfy[2] = hf[5] = Hf[id9+5];
	gfz[0] = hf[6] = Hf[id9+6];
	gfz[1] = hf[7] = Hf[id9+7];
	gfz[2] = hf[8] = Hf[id9+8];

	gf[0]  = Gradientf[id3+0];
	gf[1]  = Gradientf[id3+1];
	gf[2]  = Gradientf[id3+2];

	MatrixMultiply(Hfx,3,3,gf,3,1,tmp1);
	MatrixMultiply(hf,3,3,gfx,3,1,tmp2);
	G_HfGf[id9+0] = tmp1[0] + tmp2[0];
	G_HfGf[id9+3] = tmp1[1] + tmp2[1]; 
	G_HfGf[id9+6] = tmp1[2] + tmp2[2];


	MatrixMultiply(Hfy,3,3,gf,3,1,tmp1);
	MatrixMultiply(hf,3,3,gfy,3,1,tmp2);
	G_HfGf[id9+1] = tmp1[0] + tmp2[0];
	G_HfGf[id9+4] = tmp1[1] + tmp2[1]; 
	G_HfGf[id9+7] = tmp1[2] + tmp2[2];

	MatrixMultiply(Hfz,3,3,gf,3,1,tmp1);
	MatrixMultiply(hf,3,3,gfz,3,1,tmp2);
	G_HfGf[id9+2] = tmp1[0] + tmp2[0];
	G_HfGf[id9+5] = tmp1[1] + tmp2[1]; 
	G_HfGf[id9+8] = tmp1[2] + tmp2[2];

//printf("\nG_HfGf=%f %f %f %f %f %f %f %f %f ", G_HfGf[id9+0],G_HfGf[id9+1],G_HfGf[id9+2], G_HfGf[id9+3],G_HfGf[id9+4],G_HfGf[id9+5],G_HfGf[id9+6],G_HfGf[id9+7], G_HfGf[id9+8]); 
    }
}


/************************************************************************
 Descriptions:
     
*************************************************************************/
void Reconstruction:: ComputeHGfGHGf(int n, vector<vector<int> > &ijk, float *Hf, float *Gradientf,
                                     float *HGf, float *GHGf)
{
  int i, j, k, ix, iy, iz;

  for ( i = 0; i < n; i++ )
	{
	  GHGf[i] = 0.0;

/*
	  ix = ijk[3*i + 0];
	  iy = ijk[3*i + 1];
	  iz = ijk[3*i + 2];
*/

         ix = ijk.at(0).at(i);
         iy = ijk.at(1).at(i);
         iz = ijk.at(2).at(i);


	for ( j = 0; j < 3; j++ )
	  {
		HGf[3*i + j] = 0.0;

	  for ( k = 0; k < 3; k++ )
		HGf[3*i+j] = HGf[3*i+j] + Hf[9*i + 3*j + k] * Gradientf[3*i + k];
	  }

	for ( j = 0; j < 3; j++ )
	  GHGf[i] = GHGf[i] + Gradientf[3*i + j] * HGf[3*i + j];

	//if(fabs(GHGf[i]) >0.0)
	// printf("\nix iy iz=%d %d %d", ix, iy, iz);

	}



}



/************************************************************************
 Descriptions:
     Compute the energy functional variational J2.

 ************************************************************************/
float Reconstruction:: ComputeJ2(int i1,int j1,int k1, int n, vector<vector<int> >  &ijk, 
								 float *Gradientf, float *Hf, float *HGf, float *GHGf,
                                 float *Ngradf, float *G_HfGf, float *DeltafIsoC, 
                                 float *PdeltafIsoC)
{
  float J2f, cJ2f,J2=0.0;
  float Gradphi[3], Hphi[9], GHphiGf, HphiGf[3], Phi_ijk;
  int   i, j, k, ii, ix, iy, iz, s1, s2, id9;
  float partials[10], hf[9], tmp1[9], tmp2[9], gf[3], vec1[3], vec2[3], g_hfgf[9], term1, term2, term3, div, ngradf, epsilon;
 
  epsilon = 0.0001;
 

  for ( i = 0; i < n; i++ )
	{
	  s1 = 3*i;  
	  s2 = 9*i;
	  
/*
	  ix = ijk[s1 + 0];
	  iy = ijk[s1 + 1];
	  iz = ijk[s1 + 2];
*/

          ix = ijk.at(0).at(i);
          iy = ijk.at(1).at(i);
          iz = ijk.at(2).at(i);
	  
	  GHphiGf = 0.0;
	  bspline->Ortho_Phi_ijk_Partials_ImgGrid(ix, iy, iz, i1, j1, k1, partials);

	  Phi_ijk     = partials[0];	  
	  Gradphi[0] = partials[1];
	  Gradphi[1] = partials[2];
	  Gradphi[2] = partials[3];
	  
	  Hphi[0] =  partials[4];
	  Hphi[1] =  partials[5];
	  Hphi[2] =  partials[6];
	  Hphi[3] =  partials[5];
	  Hphi[4] =  partials[7];
	  Hphi[5] =  partials[8];
	  Hphi[6] =  partials[6];
	  Hphi[7] =  partials[8];
	  Hphi[8] =  partials[9];
	  
	  
	  for ( j = 0; j < 3; j++ )
		{
		  HphiGf[j] = 0.0;
		  
		  for ( k = 0; k < 3; k++ )
			HphiGf[j] = HphiGf[j] + Hphi[3*j + k] * Gradientf[s1 + k];
		}
	  
	  for ( j = 0; j < 3; j++ )
		GHphiGf = GHphiGf + Gradientf[s1 + j] * HphiGf[j];
	  
	//compute term 2.
	  hf[0] = Hf[s2+0];     
	  hf[1] = Hf[s2+1];
	  hf[2] = Hf[s2+2];
	  hf[3] = Hf[s2+3];
	  hf[4] = Hf[s2+4];
	  hf[5] = Hf[s2+5];
	  hf[6] = Hf[s2+6];
	  hf[7] = Hf[s2+7];
	  hf[8] = Hf[s2+8];
	  
	  gf[0]  = Gradientf[s1+0];
	  gf[1]  = Gradientf[s1+1];
	  gf[2]  = Gradientf[s1+2];
	  //printf("gf=%f %f %f ", gf[0], gf[1], gf[2]);

	  MatrixMultiply(hf,3,3,hf,3,3,tmp1);
	  MatrixMultiply(tmp1,3,3,gf,3,1,vec1);  

	  g_hfgf[0] = G_HfGf[s2+0];
 	  g_hfgf[1] = G_HfGf[s2+1];
 	  g_hfgf[2] = G_HfGf[s2+2];
 	  g_hfgf[3] = G_HfGf[s2+3];
 	  g_hfgf[4] = G_HfGf[s2+4];
 	  g_hfgf[5] = G_HfGf[s2+5];
 	  g_hfgf[6] = G_HfGf[s2+6];
 	  g_hfgf[7] = G_HfGf[s2+7];
 	  g_hfgf[8] = G_HfGf[s2+8];
 
	  MatrixMultiply(g_hfgf,3,3,gf,3,1,vec2);

	  vec1[0] = vec1[0] + vec2[0];
	  vec1[1] = vec1[1] + vec2[1];
	  vec1[2] = vec1[2] + vec2[2];

	  if(Ngradf[i] == 0.0 ) ngradf = epsilon;
          else ngradf = Ngradf[i];
 
   	  term2 = -InnerProduct(vec1,gf)/ngradf*Phi_ijk;

	  J2f = term2;
     // printf("\nterm2=%f ", term2);
      //if(fabs(Ngradf[i])<SMALLFLOAT) {printf("\n%e ", Ngradf[i]); getchar();}
	  //compute term 1 delta(g) ||delta(f)||.
	  //printf("\nJ2=%f ",J2);
	  
	  term1 = (2*(Gradphi[0] * HGf[s1+0] + Gradphi[1]*HGf[s1+1] + Gradphi[2]*HGf[s1+2])+ GHphiGf)*Ngradf[i] ;
	  J2f =  J2f + term1;
	  //printf("term1=%f ", term1);


	  // J2f[ii] = (J2f[ii] + GHphiGf[i]) * Ngradf[i] * DeltafIsoC[i];
	  //compute term 3.
	  //printf("\nJ2=%f ",J2);
	  //if(fabs(Ngradf[i]) <= SMALLFLOAT ) div = 0.0;
	  //else 
		div =  ((hf[0] + hf[4] + hf[8])/ngradf -
				(HGf[s1+0]+HGf[s1+1]+HGf[s1+2])/(ngradf*ngradf*ngradf));
		term3 = div * GHGf[i] * Phi_ijk; 
		//printf("div=%f g=%f term3=%f",div, GHGf[i], term3);

	  J2f =  J2f - term3;
	  //if(ix== 2 && iy== 2 && iz==1) 
	  //printf("Ngradf=%f term1 2 3=%f %f %f div=%f gf=%f %f %f",Ngradf[i], term1, term2, term3, div, gf[0], gf[1], gf[2]);
	  // compute integral.
	 

	  J2 = J2 + J2f*DeltafIsoC[i]*Vcube;
	  //printf("J2=%f DeltafisoC=%f Vcube=%f ", J2, DeltafIsoC[i], Vcube);
	}
 
  /*
  for ( i = 0; i < N-2; i++ )
	for ( j = 0; j < N-2; j++ )
	  for ( k = 0; k < N-2; k++ )
		{
		  cJ2f = EvaluateCubeCenterValue(J2f,i,j,k);
		  J2 = J2 + cJ2f*Vcube;
		  
		}          
  */
  //getchar();

  //free(J2f);         J2f = NULL;
  //free(Gradphi); Gradphi = NULL;
  //free(Hphi);    Hphi    = NULL;
  //free(GHphiGf); GHphiGf = NULL;
  //free(HphiGf);  HphiGf  = NULL;
  //free(Phi_ijk); Phi_ijk = NULL;
  return J2;
  
}

float Reconstruction::ComputeJ3(int i1, int j1, int k1, int n, vector<vector<int> > &ijk, float *Gradientf, float *Ngradf, float *DeltafIsoC, float *PdeltafIsoC)
{

  float J3f , cJ3f,J3=0.0, phi_ijk;
  float Gradphi[3], GfGphi;
  int   i, j, k, ix, iy, iz, s1;
  float partials[10], epsilon, ngradf;

  epsilon = 0.000001;

  printf("/n, n ============== %d/n", n);
  return 0.0;

  for ( i = 0; i < n; i++ )
	{
         s1 = 3*i;  
/*
 	ix = ijk[s1 + 0];
        iy = ijk[s1 + 1];
	iz = ijk[s1 + 2];
*/

          ix = ijk.at(0).at(i);
          iy = ijk.at(1).at(i);
          iz = ijk.at(2).at(i);



         GfGphi = 0.0;

	bspline->Ortho_Phi_ijk_Partials_ImgGrid(ix, iy, iz, i1, j1, k1, partials);

        Gradphi[0] = partials[1];
	Gradphi[1] = partials[2];
        Gradphi[2] = partials[3];

	phi_ijk    = partials[0];
   
	//printf("\n Gradphi = %f %f %f ",partials[1],  partials[2], partials[3]);
	for ( j = 0; j < 3; j++ )
	  GfGphi = GfGphi + Gradientf[s1 + j] * Gradphi[j];

	if(fabs(Ngradf[i]) == 0.0 ) ngradf = Ngradf[i] + epsilon;
        else ngradf = Ngradf[i];

	J3f = PdeltafIsoC[i] * Ngradf[i] * phi_ijk + GfGphi/ngradf * DeltafIsoC[i]; 

	//if(fabs(PdeltafIsoC[i] * Ngradf[i] * phi_ijk) > 0.001 )
	  //printf("\nterm1=%f Pdelta=%f phi_ijk=%f ", PdeltafIsoC[i] * Ngradf[i] * phi_ijk, PdeltafIsoC[i], phi_ijk);

	J3 = J3 + J3f*Vcube;

	}
  
  /*
  for ( i = 0; i < ImgNx; i++ )
	for ( j = 0; j < ImgNy; j++ )
	  for ( k = 0; k < ImgNz; k++ )
		{
		  ii = i * img2dsize + j * ImgNz + k;

		  cJ3f = J3f[ii]; //EvaluateCubeCenterValue(J3f,i,j,k);
		  J3 = J3 + cJ3f*Vcube;
		  //printf("\ncJ3f=%f Vcube=%f", cJ3f, Vcube);
		}
  */          
  //getchar();

  return J3;
  
}


void Reconstruction::ComputeHc_f(float* f, float* IsoC, vector<int> &outliers)
{
  for (int i = 0; i < VolImgSize; i++ )
	if ((*IsoC - f[i])>= 0.0 )  outliers.push_back(i);



}



float Reconstruction::ComputeJ4(int i1, int j1, int k1, int n, vector<vector<int> > &ijk, float* f, float *DeltafIsoC, vector<int> *outliers)
{
  int        i, s, ix, iy, iz, id, n2, remainder;
  float      partials[10], term1, term2;

  //compute term1.

  term1 = 0.0;
  for ( i = 0; i < n; i++ )
	{
      s = i + i + i;

/*	  ix = ijk[s + 0];
	  iy = ijk[s + 1];
	  iz = ijk[s + 2];
*/

          ix = ijk.at(0).at(i);
          iy = ijk.at(1).at(i);
          iz = ijk.at(2).at(i);


	  id = ix * img2dsize + iy * ImgNy + iz;

	  bspline->Ortho_Phi_ijk_Partials_ImgGrid(ix, iy, iz, i1, j1, k1, partials);

	  term1 += DeltafIsoC[i] * f[id] * f[id] * partials[0];
	}

  term1 = term1 * Vcube;

  //compute term2.
  n2 = outliers->size();
  term2 = 0.0;
  for ( i = 0; i < n2; i++ )
	{
	  s = outliers->at(i);
	  ix        = s /img2dsize;
	  remainder = s%img2dsize;
 	  iy        = remainder/ImgNz;
      iz        = remainder%ImgNz;
	  id = ix * img2dsize + iy * ImgNz + iz;
	  bspline->Ortho_Phi_ijk_Partials_ImgGrid(ix, iy, iz, i1, j1, k1, partials);

	  term2 += f[id] * partials[0];
	}
  term2 = term2 * Vcube;
 
  return (term1 + term2);
}


void Reconstruction::ComputeJ4outliers(float *J234phi, float *f, vector<int> &outliers, float tauga)
{
 int i, j, k, s, i1, ii, jj, n, ix, iy, iz, remainder, id, jd;
 int lx, rx, ly, ry, lz, rz, scale, usedN2;
 float X[3], Bx, By, Bz, Bxy,  rscale, cubev;
  
  scale  = (int)Bscale;
  rscale = 1.0 / Bscale;
  cubev  = 1.0;
  cubev = cubev*tauga;

  usedN2  = (usedN+1)*(usedN+1);
  //compute term2.
  n = outliers.size();
  for ( i1 = 0; i1 < n; i1++ )
	{
	  s = outliers.at(i1);
	  ix        = s /img2dsize;
	  remainder = s%img2dsize;
 	  iy        = remainder/(ImgNx+1);
          iz        = remainder%(ImgNx+1);
	  ii = ix * img2dsize + iy * (ImgNx+1) + iz;

/*
          lx = ix/scale-3;  
          rx = lx+3;           
          if (lx < 0) lx = 0;
          if (rx > usedN) rx = usedN;

          ly = iy/scale-3;
          ry = ly+3;          
          if (ly < 0) ly = 0;
          if (ry > usedN) ry = usedN;

          lz = iz/scale-3;
          rz = lz+3;         
          if (lz < 0) lz = 0;
          if (rz > usedN) rz = usedN;

          if(lx > usedN || ly > usedN || lz > usedN ) {printf("\nerror. Lx ly lz is bigger than usedN.line 1993.\n ");getchar();}
*/
          X[0] = ix*rscale - 2.0;
          X[1] = iy*rscale - 2.0;
          X[2] = iz*rscale - 2.0;
 
         // Loop for bsisis functions.  
         for ( i = lx; i<= rx; i++ ) {
           if( i < 0 || i > usedN ) continue;           
            bspline->Spline_N_Base(X[0]-i, &Bx);
           id  = i * usedN2;

           for ( j = ly; j <= ry; j++ ) {
               if( j < 0 || j > usedN ) continue;
               bspline->Spline_N_Base(X[1]-j, &By);
               jd = j *(usedN+1);
               Bxy = Bx*By;
               for ( k = lz; k <= rz; k++ ) {
                 if( k < 0 || k > usedN ) continue;
                 bspline->Spline_N_Base(X[2]-k, &Bz);
                 jj = id + jd + k;
	         J234phi[jj] += 2*f[ii] * Bxy*Bz*cubev;
	   }
        }
     }
  }

}



float Reconstruction::ComputeJ5(int i1, int j1, int k1, int n, vector<vector<int> > &ijk, vector<float> &DeltafIsoC, float *HGf,float* Ngradf , float *Gradientf, float *Hf, float V0)
{
  int i, j, ix, iy, iz, s1, id9;
  float result1,result2, X[3], normal[3], half, temp, temp1[3], m[9], Imatrix[9];
  float partials[10], gphi[3], div, r2, epsilon, ngradf;
  result1 = 0.0;
  result2 = 0.0;
  half = ImgNx/2.0;

  epsilon = 0.000001;

  for ( i = 0; i < 9; i++ ) Imatrix[i] = 0.0;
  Imatrix[0] = Imatrix[4] = Imatrix[8] = 1.0;

  for ( i = 0; i < n; i++ )
	{
	  s1 = 3*i;  
	  id9 = 9 * i;
/*
	  ix = ijk[s1 + 0];
	  iy = ijk[s1 + 1];
	  iz = ijk[s1 + 2];
*/
          ix = ijk.at(0).at(i);
          iy = ijk.at(1).at(i);
          iz = ijk.at(2).at(i);



	  X[0] = ix - half;
	  X[1] = iy - half; 
	  X[2] = iz - half;

	  normal[0] = Gradientf[s1+0];
	  normal[1] = Gradientf[s1+1];
	  normal[2] = Gradientf[s1+2];

          ngradf = Ngradf[i];
         if( Ngradf[i] == 0.0 ) ngradf= epsilon;

	  //compute part1.
	  result1 += InnerProduct(X, normal) * DeltafIsoC[i];

	  //compute part2.
	  normal[0] = normal[0]/ngradf;
	  normal[1] = normal[1]/ngradf;
	  normal[2] = normal[2]/ngradf;

	  MatrixMultiply(normal,3,1,normal,1,3,m);
	  for ( j = 0; j < 9; j++ )
		m[j] = Imatrix[j] - m[j];
	  
	  MatrixMultiply(X, 1, 3, m, 3, 3, temp1);
	  bspline->Ortho_Phi_ijk_Partials_ImgGrid(ix, iy, iz, i1, j1, k1, partials);
	  gphi[0] = partials[1];
	  gphi[1] = partials[2];
	  gphi[2] = partials[3];

 	
	  r2 = InnerProduct(temp1, gphi) - partials[0];
	  div =  ((Hf[id9+0] + Hf[id9+4] + Hf[id9+8])/ngradf -
				(HGf[s1+0]+HGf[s1+1]+HGf[s1+2])/(ngradf*ngradf*ngradf));
	  
      div = div * InnerProduct(X, normal) * partials[0];
	  r2  = r2 - div;
	  result2 += r2 * DeltafIsoC[i];


	}

  result1 = result1 * 1.0/3 * Vcube;
  result1 = result1 - V0;

  result2 = result2 * Vcube;

  return(result1 * result2 * 2.0/3);



}

void Reconstruction::ComputeJ5_Xu(float *J5phi, float taula, int n, vector<vector<int> > &subijk, vector<CUBE> &cubes1, 
                                   //int subfactor, float *DeltafIsoC, float *HGf, float *Ngradf , float *Gradientf, 
                                   int subfactor, float *DeltafIsoC, float *GHGf, float *Ngradf , float *Gradientf, 
                                   float *Hf, float *XtN,  float *NNt, float *XtP, float V0)
{

  int    i, j, i1, j1, k1, ix, iy, iz, jj, id, jd, half, s1, id9, usedN1, usedN2;
  float  Imatrix[9], result1, result2, epsilon, r2;
  float  partials[4], X[3], normal[3], gphi[3], div, Div, ngradf, rngradf, temp, temp1[3], del;
  int    lx, rx, ly, ry, lz, rz, scale, jly, klz;
  float  Bx, Bx1, By[10], By1[10], Bz[10], Bz1[10], Bx1y, Bxy1, Bxy, Bxyz, values[3], rscale, subcubeV;
  //float  Bx, Bx1, By, By1, Bz, Bz1, Bx1y, Bxy1, Bxy, Bxyz, values[3], rscale, subcubeV;

  CUBE   onecube;
  float  OX, OY, OZ, R;

  OX = ImgNx/2.0;
  OY = ImgNx/2.0;
  OZ = ImgNx/2.0;

  scale   = (int)Bscale;
  rscale  = 1.0/Bscale;
  usedN1 = usedN + 1;
  usedN2  = usedN1*usedN1;
  half    = ImgNx/2;
  del     = 1.0/subfactor;
  epsilon = 0.000001;

  for ( i = 0; i < 9; i++ ) Imatrix[i] = 0.0;
  Imatrix[0] = Imatrix[4] = Imatrix[8] = 1.0;

 //compute part1. 
  result1 = 0.0;
  for ( i = 0; i < n; i++ ) {
     s1 = 3*i;  
     id9 = 9 * i;

     onecube = cubes1.at(i);

     // Xu comment out
     /*
     ix = (int)onecube.ox;
     iy = (int)onecube.oy;
     iz = (int)onecube.oz;
     */

     X[0] = onecube.ox - half;
     X[1] = onecube.oy - half;
     X[2] = onecube.oz - half;

     subcubeV = onecube.l;    // * onecube.l * onecube.l;  //volume of sub cube.

     normal[0] = Gradientf[s1  ];
     normal[1] = Gradientf[s1+1];
     normal[2] = Gradientf[s1+2];
			
     result1 += InnerProduct(X, normal) * DeltafIsoC[i] * subcubeV;
  }
  
  result1 = result1 * 0.33333333333 ;
  result1 = result1 - V0;
  
  // New method 
  // Loop for cubes in the narrow band
  for ( i = 0; i < n; i++ ) {
     s1 = 3*i;  
     id9 = 9 * i;

     onecube = cubes1.at(i);

     ix = (int)onecube.ox;
     iy = (int)onecube.oy;
     iz = (int)onecube.oz;

     X[0] = onecube.ox*rscale-2.0;
     X[1] = onecube.oy*rscale-2.0;
     X[2] = onecube.oz*rscale-2.0;
     subcubeV = onecube.l;  

     lx = ix/scale-3;  
     rx = lx+3;           
     if (lx < 0) lx = 0;
     if (rx > usedN) rx = usedN;

     ly = iy/scale-3;
     ry = ly+3;          
     if (ly < 0) ly = 0;
     if (ry > usedN) ry = usedN;

     lz = iz/scale-3;
     rz = lz+3;
     if (lz < 0) lz = 0;
     if (rz > usedN) rz = usedN;
      
     // Xu added the following 12 lines 
     for ( j1 = ly; j1 <= ry; j1++ ) {
        jly = j1 - ly;
        bspline->Spline_N_Base_2(X[1]-j1, values);
        By[jly]  = values[0];
	By1[jly]  = values[1]*rscale; 
     }

     for ( k1 = lz; k1 <= rz; k1++ ) {
        klz = k1 - lz;
        bspline->Spline_N_Base_2(X[2]-k1, values);
        Bz[klz]  = values[0];
        Bz1[klz] = values[1]*rscale; 
     }


     ngradf = Ngradf[i];
     if( Ngradf[i] <= epsilon ) ngradf= epsilon;
     rngradf   = 1.0/ngradf;

     //for ( j = 0; j < 9; j++ ) 
     //   m[j] = Imatrix[j] - NNt[id9+j];

     temp1[0] = XtP[s1]; 
     temp1[1] = XtP[s1+1]; 
     temp1[2] = XtP[s1+2];
     	
     Div =  ((Hf[id9] + Hf[id9+4] + Hf[id9+8])*rngradf -
              //(HGf[s1]+HGf[s1+1]+HGf[s1+2])*rngradf*rngradf*rngradf);
             - (GHGf[i])*rngradf*rngradf*rngradf);  // Xu changed

     /*
     //Test for compute R and Mean curvature.
     R    = sqrt((onecube.ox-OX)*(onecube.ox-OX)+(onecube.oy-OY)*(onecube.oy-OY)+(onecube.oz-OZ)*(onecube.oz-OZ));
     if( R < (ImgNx-4.0)/2.0 ) {
          printf("\nH=====%f ", -Div*0.5);
          printf(" rR======%f R=%f", 1.0/R, R);
     }
     */

	  // Loop for bsisis functions 
     for ( i1 = lx; i1 <= rx; i1++ ) {
        //if( i1 < 0 || i1 > usedN ) continue;  // Xu comment this line
        bspline->Spline_N_Base_2(X[0]-i1, values);
        Bx  = values[0];
        Bx1 = values[1]*rscale; 
        id  = i1*usedN2;

	for ( j1 = ly; j1 <= ry; j1++ ) {
	   //if( j1 < 0 || j1 > usedN ) continue;  // Xu comment this line
	   //bspline->Spline_N_Base_2(X[1]-j1, values);
	   //By  = values[0];
	   // By1 = values[1]*rscale; 
	   jly = j1-ly;
	   Bxy = Bx * By[jly];
           Bx1y= Bx1* By[jly];
           Bxy1= Bx * By1[jly];
	   jd  = j1 * usedN1;
 
	   for ( k1 = lz; k1 <= rz; k1++ ) {
	      //if( k1 < 0 || k1 > usedN ) continue; // Xu comment this line

	      jj  = id + jd + k1;
 
	      //bspline->Spline_N_Base_2(X[2]-k1, values);
              //Bz  = values[0];
	      //Bz1 = values[1]*rscale; 

              klz = k1 - lz;
	      Bxyz = Bxy*Bz[klz];

              gphi[0] = Bx1y*Bz[klz]; 
              gphi[1] = Bxy1*Bz[klz]; 
              gphi[2] = Bxy*Bz1[klz];

              r2 =  temp1[0]*gphi[0]+temp1[1]*gphi[1]+temp1[2]*gphi[2] - Bxyz;   //partials[0];
			
              div = Div * XtN[i] * Bxyz; 
	      r2  = r2 - div;

	      J5phi[jj] += r2 * DeltafIsoC[i] * subcubeV;
           }
	}
     }
  }

  
  for ( i = 0; i < usedtolbsp; i++ ) {
     //if (J5phi[i] != 0.0 && taula != 0.0) {
	J5phi[i] = taula*result1*J5phi[i]*0.6666666667;
     //}
  }
  // end of new method


  
  /*
  // Old method
  // Loop for the basis functions
  for ( i1 = 0; i1 <= usedN; i1++ ) {
     for ( j1 = 0; j1 <= usedN; j1++ ) {
        for ( k1 = 0; k1 <= usedN; k1++ ) {
		
           jj  = i1*usedN2+j1*(usedN+1) + k1;
		
	   //result1 = 0.0;
           result2 = 0.0;

           lx = i1 *scale;   
	   rx = (lx + 4)*scale;

	   ly = j1 *scale;
	   ry = (ly + 4)*scale;

           lz = k1 *scale;
	   rz = (lz + 4)*scale;

           // Loop for cubes in the narrow band 
	   for ( i = 0; i < n; i++ ) {
              s1 = 3*i;  
	      id9 = 9 * i;

              onecube = cubes1.at(i);
              // Xu asked, are these three lines need to be changed to floor
              ix = (int)onecube.ox;
              iy = (int)onecube.oy;
              iz = (int)onecube.oz;

			  //                        X[0] = onecube.ox - half;
			  //         X[1] = onecube.oy - half;
			  //        X[2] = onecube.oz - half;

			  //        subcubeV = onecube.l ;    // * onecube.l * onecube.l;  //volume of sub cube.
//old methods.
//                        ix = subijk.at(0).at(i);
//                        iy = subijk.at(1).at(i);
//                        iz = subijk.at(2).at(i);
 			
//			X[0] = ix*del - half;
//			X[1] = iy*del - half; 
//			X[2] = iz*del - half;
			
//			normal[0] = Gradientf[s1+0];
			  //		normal[1] = Gradientf[s1+1];
			  //			normal[2] = Gradientf[s1+2];
			
			//compute part1.
			  //			result1 += InnerProduct(X, normal) * DeltafIsoC[i] * subcubeV;



	      if (ix >= lx && ix < rx && iy >= ly && iy < ry && iz >= lz && iz < rz){

		 //compute part2.
                 X[0] = onecube.ox*rscale-2.0;
                 X[1] = onecube.oy*rscale-2.0;
                 X[2] = onecube.oz*rscale-2.0;
                 subcubeV = onecube.l;    // * onecube.l * onecube.l;  //volume of sub cube.

                               //normal[0] = Gradientf[s1+0];
                               //normal[1] = Gradientf[s1+1];
                               //normal[2] = Gradientf[s1+2];


		 ngradf = Ngradf[i];
		 if( Ngradf[i] <= epsilon ) ngradf= epsilon;
                 rngradf   = 1.0/ngradf;

				//normal[0] = normal[0]*rngradf;
				//normal[1] = normal[1]*rngradf;
				//normal[2] = normal[2]*rngradf;
				
		 //MatrixMultiply(normal,3,1,normal,1,3,m);
		 for ( j = 0; j < 9; j++ ) {
                    //m[j] = Imatrix[j] - m[j];
		    m[j] = Imatrix[j] - NNt[id9+j];
                 }
		                		
		 //MatrixMultiply(X, 1, 3, m, 3, 3, temp1);
                 temp1[0] = XtP[s1]; temp1[1] = XtP[s1+1]; temp1[2] = XtP[s1+2];
		 bspline->Spline_N_Base_2(X[0]-i1, values);
                 Bx  = values[0];
                 Bx1 = values[1]*rscale;  // Xu asked, rscale need to be mutiplied???? Yes.
	
                 bspline->Spline_N_Base_2(X[1]-j1, values);
                 By  = values[0];
                 By1 = values[1]*rscale;  // Xu asked, rscale need to be mutiplied???? Yes.

                 bspline->Spline_N_Base_2(X[2]-k1, values);
                 Bz  = values[0];
                 Bz1 = values[1]*rscale; //  Xu asked, rscale need to be mutiplied???? Yes.

                 Bxyz = Bx*By*Bz;
 
                 gphi[0] = Bx1*By*Bz; //partials[1];
		 gphi[1] = Bx*By1*Bz; //partials[2];
		 gphi[2] = Bx*By*Bz1; //partials[3];
				
				
		 //r2 = InnerProduct(temp1, gphi) - Bxyz;   //partials[0];
		 r2 = temp1[0]*gphi[0]+temp1[1]*gphi[1]+temp1[2]*gphi[2] - Bxyz;   //partials[0];
		 div =  ((Hf[id9+0] + Hf[id9+4] + Hf[id9+8])*rngradf -
						(HGf[s1+0]+HGf[s1+1]+HGf[s1+2])*rngradf*rngradf*rngradf);
				
		 //div = div * InnerProduct(X, normal) * Bxyz;  //partials[0];
	         div = div * XtN[i] * Bxyz; 
		 r2  = r2 - div;
                 result2 += r2 * DeltafIsoC[i] * subcubeV;

	      }
	   }

	   if (result2 != 0.0 && taula != 0.0) {
	      J5phi[jj] = taula * result1 * result2 * 0.6666666667;
              //printf("\nJ5phi=%e ", J5phi[jj]);
	   }
        }
     }
  }
  // end of old method
  */
 

 for ( i = 0; i < usedtolbsp; i++ ) printf("\nJ5phi=%e ", J5phi[i]);
  
}

void Reconstruction::ComputeJ2J4J5(float *coefs, float *IsoC, int n, vector<CUBE> &cubes1, float V0,  float taual, float tauga, float taula, 
                                   float *J234phi, float *J5phi)
{
  int   id, i, j, k, ii, jj, kk, jd, id3, lx, rx, ly, ry, lz, rz, scale, half;
  int   ix, iy, iz, usedN1, usedN2;
  float partials[20], temp, epsilon, X[3], del, rscale, values[3], value, normal[3], m[9], ngradf, rngradf;
  float Bx[10], Bx1[10], Bx2[10], Bx3[10], By[10], By1[10], By2[10], By3[10], Bz[10], Bz1[10], Bz2[10], Bz3[10];
  float Bxy,Bxy1, Bx1y, Bx1y1, Bx2y, Bxy2, Bx3y, Bx2y1, Bx1y2, Bxy3, Imatrix[9], sum;
  int   ilx, jly, klz;
  float rscale2, rscale3, Div, div, result1, subcubeV, DeltasubV, DeltasubV_J2, DeltasubV_J4, DeltasubV_J5, DivXtN;
  CUBE  onecube;
  float DeltafIsoC, PdeltafIsoC, Hf[9], Gradientf[3], Ngradf, HGf[3], GHGf, XtN, NNt[9], XtP[3], gphi[3];
  float Hphi[9], Hfx[9], Hfy[9], Hfz[9], G_HfGf[9], Gfx[3], Gfy[3], Gfz[3], tmp[3], grad_g[3]; //, Bx3, By3, Bz3;
  float Func_phi[64], Grad_phi[192], Hsn_phi[384];
  int   kk30, kk31, kk32, kk60, kk61, kk62, kk63, kk64, kk65;

  epsilon = 0.000001;  //If Norm gradient is zero.
  scale   = (int)Bscale;
  rscale  = 1.0/Bscale;
  rscale2 = rscale*rscale;   // Xu Added
  rscale3 = rscale*rscale2;
  usedN1  = usedN+1;
  usedN2  = usedN1*usedN1;
  half    = ImgNx/2;
  result1 = 0.0;

  for ( i = 0; i < 9; i++ ) Imatrix[i] = 0.0;
  Imatrix[0] = Imatrix[4] = Imatrix[8] = 1.0;


// Loop for all the cubes in the narrow band
for ( id = 0; id < n; id++ ) {
   onecube  = cubes1.at(id);
   subcubeV = onecube.l;
 
   // Xu changed
   ix = floor(onecube.ox);
   iy = floor(onecube.oy);
   iz = floor(onecube.oz);

   lx = ix/scale-3;  
   rx = lx+3;           
   if (lx < 0) lx = 0;
   if (rx > usedN) rx = usedN;

   ly = iy/scale-3;
   ry = ly+3;          
   if (ly < 0) ly = 0;
   if (ry > usedN) ry = usedN;

   lz = iz/scale-3;
   rz = lz+3;         
   if (lz < 0) lz = 0;
   if (rz > usedN) rz = usedN;

   if(lx > usedN || ly > usedN || lz > usedN ) {printf("\nerror. Lx ly lz is bigger than usedN.line 1235.\n ");getchar();}
   X[0] = onecube.ox*rscale - 2.0;
   X[1] = onecube.oy*rscale - 2.0;
   X[2] = onecube.oz*rscale - 2.0;
 
   // Xu added the following 14 lines
   for ( i = lx; i<= rx; i++ ) {
      ilx = i - lx; 
      bspline->Spline_N_Base_2(X[0]-i, values);
      Bx[ilx]  = values[0];
      Bx1[ilx] = values[1]*rscale;
      Bx2[ilx] = values[2]*rscale2;  
   }

   for ( j = ly; j <= ry; j++ ) {
      jly = j - ly;
      bspline->Spline_N_Base_2(X[1]-j, values);
      By[jly]   = values[0];
      By1[jly]  = values[1]*rscale;
      By2[jly]  = values[2]*rscale2;
   }

   for ( k = lz; k <= rz; k++ ) {
      klz = k - lz; 
      bspline->Spline_N_Base_2(X[2]-k, values);
      Bz[klz]  = values[0];
      Bz1[klz] = values[1]*rscale;
      Bz2[klz] = values[2]*rscale2;
   }

   if (taual != 0.0 ) {
      for ( i = lx; i <= rx; i++ ) {
         ilx = i - lx;
         bspline->Spline_N_Base_3(X[0]-i, values);
         Bx3[ilx] = values[0]*rscale3;
      }

      for ( j = ly; j <= ry; j++ ) {
         jly = j - ly;
         bspline->Spline_N_Base_3(X[1]-j, values);
         By3[jly] = values[0]*rscale3;
      }

      for ( k = lz; k <= rz; k++ ) {
         klz = k - lz;
         bspline->Spline_N_Base_3(X[2]-k, values);
         Bz3[klz] = values[0]*rscale3;
      }
   }

   // set all the initial values to zero
   sum = 0.0;
   for ( i = 0; i < 3; i++ ) 
      Gradientf[i] = 0.0;
   for ( i = 0; i < 9; i++ ) 
      Hf[i]  = 0.0;

   if (taual != 0.0 ) { 
      for ( i = 0; i < 9; i++ ) {
         Hfx[i] = Hfy[i] = Hfz[i] = 0.0;
      }
   }

   kk = 0;
   // Loop for basis functions. Compute Gradientf,Hf...., up to the third order partials of f at the cube's center 
   for ( i = lx; i<= rx; i++ ) {
      ilx = i - lx; 
      jj  = i * usedN2;

      for ( j = ly; j <= ry; j++ ) {
         jly = j - ly;
         Bxy  = Bx[ilx]  * By[jly];
         Bx1y = Bx1[ilx] * By[jly];
         Bxy1 = Bx[ilx]  * By1[jly];
       
         Bx2y = Bx2[ilx] * By[jly];
         Bx1y1= Bx1[ilx] * By1[jly];
         Bxy2 = Bx[ilx]  * By2[jly];
         jd   = j * usedN1;

         if (taual != 0.0) {
            Bx3y  = Bx3[ilx]*By[jly] ;  //partials[10];
            Bx2y1 = Bx2[ilx]*By1[jly];  //partials[11];
            Bx1y2 = Bx1[ilx]*By2[jly];  //partials[13];
            Bxy3  = Bx[ilx]*By3[jly];   //partials[16];
         }

         for ( k = lz; k <= rz; k++ ) {
            klz = k - lz;  
            ii = jj + jd + k; // 

            // compute function value of phi
            Func_phi[kk] = Bxy * Bz[klz];

            // compute function value of f
	     sum = sum + coefs[ii] * Func_phi[kk];

            // compute gradient of phi
            kk30 = 3*kk;
            kk31 = kk30 + 1;
            kk32 = kk30 + 2;

            Grad_phi[kk30] = Bx1y*Bz[klz];
            Grad_phi[kk31] = Bxy1*Bz[klz];
            Grad_phi[kk32] = Bxy *Bz1[klz];

            // compute gradient of f
	    Gradientf[0] = Gradientf[0] + coefs[ii]*Grad_phi[kk30];
	    Gradientf[1] = Gradientf[1] + coefs[ii]*Grad_phi[kk31];
	    Gradientf[2] = Gradientf[2] + coefs[ii]*Grad_phi[kk32];
			
            // compute Hessian of phi 
            kk60 = 6*kk;
            kk61 = kk60 + 1;
            kk62 = kk60 + 2;
            kk63 = kk60 + 3;
            kk64 = kk60 + 4;
            kk65 = kk60 + 5;

            Hsn_phi[kk60] = Bx2y*Bz[klz];
            Hsn_phi[kk61] = Bx1y1*Bz[klz];
            Hsn_phi[kk62] = Bx1y*Bz1[klz];
            Hsn_phi[kk63] = Bxy2*Bz[klz];
            Hsn_phi[kk64] = Bxy1*Bz1[klz];
            Hsn_phi[kk65] = Bxy*Bz2[klz];

            // compute Hessian 
	    Hf[0] = Hf[0] + coefs[ii]*Hsn_phi[kk60];  //partials[4];
	    Hf[1] = Hf[1] + coefs[ii]*Hsn_phi[kk61];  //partials[5];
	    Hf[2] = Hf[2] + coefs[ii]*Hsn_phi[kk62];  //partials[6];
	    Hf[4] = Hf[4] + coefs[ii]*Hsn_phi[kk63];  //partials[7];
	    Hf[5] = Hf[5] + coefs[ii]*Hsn_phi[kk64];  //partials[8];
            Hf[8] = Hf[8] + coefs[ii]*Hsn_phi[kk65];  //partials[9];

            // Compute J2
            if(taual != 0.0 ) {

               //third order partials that is required by J2
               Hfx[0] = Hfx[0] + coefs[ii]*Bx3y *Bz[klz];  //partials[10];
               Hfx[1] = Hfx[1] + coefs[ii]*Bx2y1*Bz[klz];  //partials[11];
               Hfx[2] = Hfx[2] + coefs[ii]*Bx2y *Bz1[klz]; //partials[12];
               Hfx[4] = Hfx[4] + coefs[ii]*Bx1y2*Bz[klz];  //partials[13];
               Hfx[5] = Hfx[5] + coefs[ii]*Bx1y1*Bz1[klz]; //partials[14];
               Hfx[8] = Hfx[8] + coefs[ii]*Bx1y *Bz2[klz]; //partials[15];

               Hfy[4] = Hfy[4] + coefs[ii]*Bxy3*Bz[klz];   //partials[16];
               Hfy[5] = Hfy[5] + coefs[ii]*Bxy2*Bz1[klz];  //partials[17];
               Hfy[8] = Hfy[8] + coefs[ii]*Bxy1*Bz2[klz];  //partials[18];

               Hfz[8] = Hfz[8] + coefs[ii]*Bxy*Bz3[klz];   //partials[19];
            }  // end J2 

            kk = kk + 1;
	 }		
      }
   }    // end for loop basis functions

   // Xu added 
   Hf[3] = Hf[1];
   Hf[6] = Hf[2];
   Hf[7] = Hf[5];

   // for J2
   if(taual != 0.0 ) {
      Hfx[3] = Hfx[1];  
      Hfx[6] = Hfx[2];  
      Hfx[7] = Hfx[5];

      Hfy[0] = Hfx[1];  
      Hfy[3] = Hfy[1] = Hfx[4];  
      Hfy[6] = Hfy[2] = Hfx[5]; 
      Hfy[7] = Hfy[5];

      Hfz[0] = Hfx[2];  
      Hfz[3] = Hfz[1] = Hfx[5];  
      Hfz[6] = Hfz[2] = Hfx[8]; 
      Hfz[4] = Hfy[5];  
      Hfz[7] = Hfz[5] = Hfy[8]; 
 
      // grad of fx, fy, fz
      Gfx[0] = Hf[0]; Gfx[1] = Hf[1]; Gfx[2] = Hf[2];
      Gfy[0] = Hf[3]; Gfy[1] = Hf[4]; Gfy[2] = Hf[5];
      Gfz[0] = Hf[6]; Gfz[1] = Hf[7]; Gfz[2] = Hf[8];

      /* // Xu comment out
      MatrixMultiply(Hfx,3,3,Gradientf,3,1,values);
      MatrixMultiply(Hf,3,3,Gfx,3,1,tmp);
      G_HfGf[0] = values[0] + tmp[0];
      G_HfGf[3] = values[1] + tmp[1];
      G_HfGf[6] = values[2] + tmp[2];

      MatrixMultiply(Hfy,3,3,Gradientf,3,1,values);
      MatrixMultiply(Hf,3,3,Gfy,3,1,tmp);
      G_HfGf[1] = values[0] + tmp[0];
      G_HfGf[4] = values[1] + tmp[1];
      G_HfGf[7] = values[2] + tmp[2];

      MatrixMultiply(Hfz,3,3,Gradientf,3,1,values);
      MatrixMultiply(Hf,3,3,Gfz,3,1,tmp);
      G_HfGf[2] = values[0] + tmp[0];
      G_HfGf[5] = values[1] + tmp[1];
      G_HfGf[8] = values[2] + tmp[2];
      */

      // Xu changed , compute grad(g)
      // compute g_x
      MatrixMultiply(Hfx,3,3,Gradientf,3,1,values);
      MatrixMultiply(Hf,3,3,Gfx,3,1,tmp);
      G_HfGf[0] = values[0] + 2.0*tmp[0];
      G_HfGf[1] = values[1] + 2.0*tmp[1];
      G_HfGf[2] = values[2] + 2.0*tmp[2];
      grad_g[0] = InnerProduct(Gradientf, G_HfGf);

      // compute g_y
      MatrixMultiply(Hfy,3,3,Gradientf,3,1,values);
      MatrixMultiply(Hf,3,3,Gfy,3,1,tmp);
      G_HfGf[0] = values[0] + 2.0*tmp[0];
      G_HfGf[1] = values[1] + 2.0*tmp[1];
      G_HfGf[2] = values[2] + 2.0*tmp[2];
      grad_g[1] = InnerProduct(Gradientf, G_HfGf);

      // compute g_z
      MatrixMultiply(Hfz,3,3,Gradientf,3,1,values);
      MatrixMultiply(Hf,3,3,Gfz,3,1,tmp);
      G_HfGf[0] = values[0] + 2.0*tmp[0];
      G_HfGf[1] = values[1] + 2.0*tmp[1];
      G_HfGf[2] = values[2] + 2.0*tmp[2];
      grad_g[2] = InnerProduct(Gradientf, G_HfGf);

   }  // end for J2
 
   // compute the length of gradient
   temp = Gradientf[0]*Gradientf[0] + Gradientf[1]*Gradientf[1] + Gradientf[2]* Gradientf[2];	
   Ngradf    = sqrt(temp);

   //printf("Ngrad = %e\n", Ngradf);
  
   DeltafIsoC = Deltafunc(sum - *IsoC);
   PdeltafIsoC= DeltafuncPartials(sum - *IsoC);
   DeltasubV  = DeltafIsoC * subcubeV; 
   DeltasubV_J2  = taual*DeltasubV; 
   DeltasubV_J4  = tauga*DeltasubV; 

   // compute H*G
   HGf[0] = 0.0; 
   HGf[1] = 0.0; 
   HGf[2] = 0.0; 
   for ( j = 0; j < 3; j++ ) {
      for ( k = 0; k < 3; k++ ) {
         HGf[j] = HGf[j] + Hf[3*j + k] * Gradientf[k];
      }
   }

   // compute G^T*H*G
   GHGf = 0.0;
   for ( j = 0; j < 3; j++ ) {
       GHGf = GHGf + Gradientf[j] * HGf[j];
   }
   
   // Why this X is different from that X
   X[0] = onecube.ox - half;
   X[1] = onecube.oy - half;
   X[2] = onecube.oz - half;

  // compute J5 part1.
  if(taula != 0.0 ) result1 += InnerProduct(X, Gradientf) * DeltafIsoC * subcubeV;


   ngradf = Ngradf;
   // Xu changed
   if( ngradf <= epsilon ) ngradf= epsilon;
   rngradf   = 1.0/ngradf;

   // surface normal N
   normal[0] = Gradientf[0]*rngradf;
   normal[1] = Gradientf[1]*rngradf;
   normal[2] = Gradientf[2]*rngradf;
  
   // Compute X^T* N
   XtN = X[0]*normal[0]+X[1]*normal[1]+X[2]*normal[2];

   // Compute N*N^T
   MatrixMultiply(normal,3,1,normal,1,3,NNt);
   
   // Compute P = I - N*N^T
   for ( j = 0; j < 9; j++ ) {
       m[j] = Imatrix[j] - NNt[j];
   }

   // Compute X^T*P
   MatrixMultiply(X, 1, 3, m, 3, 3, XtP);
     	
   Div =  ((Hf[0] + Hf[4] + Hf[8])*rngradf -
           - GHGf*rngradf*rngradf*rngradf);  // Xu changed

   DivXtN = Div*XtN;

   kk = 0;
   // Loop for bsisis functions. compute J2, and J5 part 2.
   for ( i = lx; i <= rx; i++ ) {
      ii  = i*usedN2;
      for ( j = ly; j <= ry; j++ ) {
	 jd  = j * usedN1;
	 for ( k = lz; k <= rz; k++ ) {
	    jj  = ii + jd + k;

            kk30 = 3*kk;

            //compute J2.----------------------------------------------------------------------
            if (taual != 0.0 ) {
               kk60 = kk*6;
              
               // Hessian of phi
	       Hphi[0] = Hsn_phi[kk60]; 
	       Hphi[1] = Hsn_phi[kk60+1];
	       Hphi[2] = Hsn_phi[kk60+2];
	       Hphi[4] = Hsn_phi[kk60+3]; 
	       Hphi[5] = Hsn_phi[kk60+4];
               Hphi[8] = Hsn_phi[kk60+5];

               Hphi[3] = Hphi[1];
               Hphi[6] = Hphi[2];
               Hphi[7] = Hphi[5];

               // 2 (\nabla \phi)^T * H(f) * \nabla f
               temp = 2*InnerProduct(Grad_phi+kk30, HGf);

               MatrixMultiply(Hphi, 3, 3, Gradientf, 3, 1, values);
               temp = (temp + InnerProduct(Gradientf, values))*Ngradf;   // the first term of (2.7)

               /* Xu changed
               // compute gradient of g.
               MatrixMultiply(Hf, 3, 3, Hf, 3, 3, m);
               MatrixMultiply(m, 3, 3, Gradientf, 3, 1, tmp);
               MatrixMultiply(G_HfGf, 3, 3, Gradientf, 3, 1, values);
               tmp[0] = tmp[0] + values[0];
               tmp[1] = tmp[1] + values[1];
               tmp[2] = tmp[2] + values[2];    // tmp is gradient of g
               */

               //////temp = temp - InnerProduct(tmp, Gradientf)*rngradf*Func_phi[kk]; // minus the second term            

               temp = temp - InnerProduct(grad_g, Gradientf)*rngradf*Func_phi[kk]; // minus the second term            
               temp = temp - Div * GHGf * Func_phi[kk];                            // minus the third term
               J234phi[jj] += temp * DeltasubV_J2; 
               //printf("id = %d, temp = %f, DeltasubV = %e\n", id, temp, DeltasubV);
            }

            //compute J4 part1.------------------------------------------------------------------------

            if(tauga != 0.0 ) J234phi[jj] += sum * sum * Func_phi[kk] * DeltasubV_J4;



            //compute J5 part2.------------------------------------------------------------------
            if (taula != 0.0 ) {
               temp =  XtP[0]*Grad_phi[kk30] + XtP[1]*Grad_phi[kk30+1]+XtP[2]*Grad_phi[kk30+2] - Func_phi[kk];   
               div = DivXtN * Func_phi[kk]; 
	       temp = temp - div;
	       J5phi[jj] += temp * DeltasubV;
               //printf("\ngphi=%f %f %f phi=%f temp=%f DeltasubV=%f", Grad_phi[kk30], Grad_phi[kk30+1], Grad_phi[kk30+2], Func_phi[kk], temp, DeltasubV);
            }
      
            kk = kk + 1;
         }  // end k  loop
      }     // end j  loop
   }        // end i  loop for basis
}           // end id loop for cubes, which is the out loop

if (taula != 0.0 ) {
   result1 = result1 * 0.33333333333 ;
   result1 = result1 - V0;
   result1 = result1*0.6666666667*taula;
   for ( i = 0; i < usedtolbsp; i++ ) {
      J5phi[i] = result1*J5phi[i];
      //printf("\nJ5phi=%e ", J5phi[i]);
   }
}


}


/************************************************************************
 ************************************************************************/
float Reconstruction::EvaluateCubeCenterValue(float *J2f,int i, int j,int k)
{
  float result = 0.0;

  result = J2f[i*img2dsize+j*(N-1)+k] + J2f[i*img2dsize+j*(N-1)+k+1] +
	J2f[i*img2dsize+(j+1)*(N-1)+k] + J2f[i*img2dsize+(j+1)*(N-1)+k+1] +      
    J2f[(i+1)*img2dsize+j*(N-1)+k] + J2f[(i+1)*img2dsize+j*(N-1)+k+1] +
    J2f[(i+1)*img2dsize+(j+1)*(N-1)+k] + J2f[(i+1)*img2dsize+(j+1)*(N-1)+k+1];

  result = 0.125 * result;

  return result;
}





/***************************************************************************
Descriptions:

Arguments:
     N   : 0,1...,N are the total Cubic Bspline Base Functions along each 
           Volume axis.
     M   : 0,1,..,M are the total Gauss nodes for Gauss quadrature numerical 
           integral. M = (N -2)*4;

****************************************************************************/
void  Reconstruction::EvaluateCubicBsplineBaseProjectionIntegralMatrix(float *Non_o_coef)
{
  int     i,j,k,i1,j1,k1, a, b, c, ai, bj, ss;
 float   res; 
 char    filename[50];
 int     ii,jj, iii,jjj,kkk,jjjj,  id,jd;

 int      usedN2, usedn2, Begin, End;

 Begin = BStartX - BFinishX;
 End   = BFinishX - BStartX;

/*
for ( ss = 0; ss <usedtolbsp; ss++ )
   { Matrix[ss] = 0.0;
   }
*/
 usedN2 = (usedN+1) *(usedN+1);
 usedn2 = (2*usedN+1) * (2*usedN+1);

// printf("\nMatrix"); 

 //for ( i1 = 0; i1 <= usedN; i1++ )
 for ( i1 = BStartX; i1 <= BFinishX; i1++ )
   {
	 iii = (i1-BStartX) * usedN2;
    //for ( j1 = 0; j1 <= usedN; j1++ )
	for ( j1 = BStartX; j1 <= BFinishX; j1++ )
      {
		jjj = (j1-BStartX)*(usedN+1);
       //for ( k1 = 0; k1 <= usedN; k1++ )
	   for ( k1 = BStartX; k1 <= BFinishX; k1++ )
          {
           jj = iii + jjj + k1 - BStartX;
		   
		   //for ( i = 0; i <= usedN; i++ )
		   for ( i = BStartX; i <= BFinishX; i++ )
              {
				id = (i-BStartX) * usedN2;
			   //a  = i1 - i + usedN;
			   a = i1 - i - Begin;
			   ai = a * usedn2;
 
			   //for ( j = 0; j <= usedN; j++ )
			   for ( j = BStartX; j <= BFinishX; j++ )
                {
				  jd = (j-BStartX) * (usedN+1) ;
				 //b  = j1 - j + usedN;
				 b = j1 - j - Begin;
				 bj = b * (2*usedN+1);

                 //for ( k = 0; k <= usedN; k++ )
				 for ( k = BStartX; k <= BFinishX; k++ )
                    {
                     ii   = id + jd + k - BStartX;
					 //c    = k1 - k + usedN;
					 c =   k1 - k - Begin;
                     jjjj = ai + bj + c;
					 /*jjjj = (i1-BStartX)*(usedN+1)*(usedN+1)*(usedN+1)*(usedN+1)*(usedN+1)
					   +(j1-BStartX)*(usedN+1)*(usedN+1)*(usedN+1)*(usedN+1)+
					   (k1-BStartX)*(usedN+1)*(usedN+1)*(usedN+1)
					   + (i-BStartX)*(usedN+1)*(usedN+1) + 
					   (j-BStartX)*(usedN+1) + k-BStartX;*/
                     res = Mat[jjjj];

                     //Matrix[jj] = Matrix[jj] + Non_o_coef[ii]*res; 
                     xdgdphi[jj] = xdgdphi[jj] + Non_o_coef[ii]*res; 
		//			  if(res != 0 ) printf("res=%f Non_o_coef=%f    ", res, Non_o_coef[ii]);
					} 
                }
             }
		   //printf("\ni1 j1 k1 = %d %d %d Matrix=%f gdxdphi=%f cha=========%e ", i1, j1, k1, Matrix[jj], gdxdphi[jj], (Matrix[jj] - gdxdphi[jj])/Matrix[jj]);

		   //if(jj == 0 ) getchar();
		  
          }
       }
    }

 for ( i = 0; i < usedtolbsp; i++ ) xdgdphi[i] = xdgdphi[i] - gdxdphi[i];
//getchar();
}





void Reconstruction::EvaluateFFT_GdXdPhiMatrix()
{
  int i, j, i1, j1, k1, N1, half, img2dsize3, imgsize;
  int  ii, jj, iii,jjj, N2, isize, jj3;
  fftw_complex *slice, *B_k;
  Oimage *image;
  float  rotmat[9], oldx[3], X1[3],  X_2, Y_2, Z_2, a, b;
  float X, Y, Z, s, t, sX, sY, sZ, wk_x, wk_y; 
  float e1[3], e2[3], d[3], scale;

	N1         = usedN+1;
	N2         = N1*N1;
	img2dsize3 = img2dsize * 3;
	half       = ImgNx/2;
	scale      = PI2/(ImgNx+1);
	 
	a          = -ImgNx/2.0;
	b          =  ImgNx/2.0;
	

  slice=(fftw_complex *)fftw_malloc(img2dsize*sizeof(fftw_complex));
  B_k = (fftw_complex *)fftw_malloc(sizeof(fftw_complex));
  
  slicemat=(float *)malloc(nv * 3 * img2dsize * sizeof(float));
  ReBFT = (float *)malloc(nv * img2dsize * 2 * sizeof(float));
  
  for ( i = 0; i < nv*3*img2dsize; i++ ) slicemat[i] = 0.0;
  for ( i = 0; i < nv * img2dsize; i++ ) ReBFT[i] = 0.0;
  
  
  for (i = 0; i < nv; i++ )
	{
	  isize   = i * img2dsize3;
	  imgsize = 2* i * img2dsize;
	  
	  for ( j = 0; j < 9; j++ ) {
		rotmat[j] = Rmat[i*9 + j];
		//	   printf("\nrotmat=%f ", rotmat[j]);
	  }
	  //	 getchar();
	  
	  for ( i1 = 0; i1 < 3; i1++ )
		{
		  e1[i1] = rotmat[i1];
		  e2[i1] = rotmat[3+i1];  
		  d[i1]  = rotmat[6+i1];
		}
	  
	  for ( i1 = -half; i1 <= half; i1++ )
		{
		  
		 //printf("\nPI2=%f ", PI2); getchar();
		  ii = (i1 + half)  * (ImgNx+1);
		 
		  for ( j1 = -half; j1 <= half; j1++ )
			{
			  // oldx[1] = j1*PI2/(ny*1.0);
			  
			  jj = ii + j1 + half;
			  jj3 = 3*jj;
			  


			  wk_x = i1*PI2/((b-a)+1.0);
			  wk_y = j1*PI2/((b-a)+1.0);

/*			  wk_x  = i1 * scale; 
			  wk_y  = j1 * scale;
*/
			  X1[0] = wk_x * e1[0] + wk_y * e2[0];
			  X1[1] = wk_x * e1[1] + wk_y * e2[1];
			  X1[2] = wk_x * e1[2] + wk_y * e2[2];

			  			  
			  slicemat[isize+jj3+0] = X1[0];
			  slicemat[isize+jj3+1] = X1[1];
			  slicemat[isize+jj3+2] = X1[2];

			 //Bspline FT central slice Real parts.
			 X = Bscale * X1[0]; 
			 Y = Bscale * X1[1]; 
			 Z = Bscale * X1[2];

			 X_2 = X/2.0;
			 Y_2 = Y/2.0;
			 Z_2 = Z/2.0;
			 
			 if( fabs(X_2) == 0.0) sX = 1.0;
			 else sX = sin(X_2)/X_2;
			 
			 if( fabs(Y_2) == 0.0 ) sY = 1.0;
			 else sY = sin(Y_2)/Y_2;
			 
			 
			 if( fabs(Z_2) == 0.0 ) sZ = 1.0;
			 else sZ = sin(Z_2)/Z_2;
			 //printf("\nsZ=%f", sZ ); 
			 s = sX * sY * sZ;
			 
			 t = s * s;
			 s = t * t;
			 ReBFT[imgsize + 2*jj+0] = fabs(Bscale*Bscale*Bscale)* s;
			 ReBFT[imgsize + 2*jj+1] = fabs(Bscale*Bscale*Bscale)* fabs(Bscale*Bscale*Bscale)*s*s;

			 gd_FFT[i * img2dsize +jj][0] = gd_FFT[i * img2dsize +jj][0] * ReBFT[imgsize + 2*jj+0];
			 gd_FFT[i * img2dsize +jj][1] = gd_FFT[i * img2dsize +jj][1] * ReBFT[imgsize + 2*jj+0];

		   }
	   }
   }

 
 image  = InitImageParameters(1, ImgNx+1, ImgNx+1, 1, 1);

	 image->ux = scale;
	 image->uy = scale;
	 image->uz = 1.0;
	 image->nx = ImgNx+1; 
	 image->ny = ImgNx+1;
	 image->nz = 1;

for ( i1 = 0; i1 < img2dsize; i1++ )
   {slice[i1][0] = 0.0;slice[i1][1] = 0.0;}

 gdxdphi = (float *)malloc(usedtolbsp*sizeof(float));
 for ( jj = 0; jj < usedtolbsp; jj++ ) gdxdphi[jj] = 0.0;

 printf("\nBegin compute GdXdPhiIntegral\n");
 for ( i1 = BStartX; i1 <= BFinishX; i1++ )
   {
	 iii = (i1 - BStartX) * N2;
	 
	 for ( j1 = BStartX; j1 <= BFinishX; j1++ )
	   {
		 jjj = (j1 - BStartX) * N1;
		 for ( k1 = BStartX; k1 <= BFinishX; k1++ )
		   {
			 
			 jj = iii + jjj + k1 - BStartX;
			 //printf("\ni1 j1 k1=%d %d %d ", i1, j1, k1);

			 gdxdphi[jj] = EvaluateFFT_GdXdPhiIntegral(i1, j1, k1, slice, B_k, image);
			 
			 //printf("%f  ", gdxdphi[jj]);
		   } 
	   }
   }
 printf("...Finished. \n");
 //getchar();
 
 fftw_free(slice);   slice = NULL;
 fftw_free(B_k);     B_k   = NULL;
 
 kill_all_but_main_img(image);    
 free(image); image = NULL;
 // free(slicemat);

//getchar();
}


void Reconstruction::EvaluateFFT_GdXdOrthoPhiMatrix()
{
  int i, j, k, i1, j1, k1, ii, jj, kk, ii1, jj1, kk1;
  int  N2, s;
  float *gdxdOphi, ai, bj, ck, aibj, aibjck;


  gdxdOphi = (float *)malloc(usedtolbsp * sizeof(float));
  for ( s = 0; s < usedtolbsp; s++ )
	gdxdOphi[s] = 0.0;

 
 
  N2 = (usedN + 1) * (usedN + 1);

  for (i = 0; i <= usedN; i++ )
	{
	  ii = i * N2;
	  for ( j = 0; j <= usedN; j++ )
		{
		  jj = j * ( usedN+1);
		  for ( k = 0; k <= usedN; k++ )
			{
			  kk = ii + jj + k;

			  for ( i1 = 0; i1 <= i; i1++ )
				{
                    ai  = (float)bspline->SchmidtMat[i*(usedN+1)+i1] ;
					ii1 = i1 * N2;

					for ( j1 = 0; j1 <= j; j1++ )
					  {
						bj = (float)bspline->SchmidtMat[j*(usedN+1)+j1];
						aibj = ai * bj;

						jj1 = j1 * (usedN + 1);

						for ( k1 = 0; k1 <= k; k1++ )
						  {
							ck = (float)bspline->SchmidtMat[k*(usedN+1)+k1]; 

							aibjck = aibj * ck;
							kk1 = ii1 + jj1 + k1;
							gdxdOphi[kk] = gdxdOphi[kk] + aibjck * gdxdphi[kk1];
						  }  
					  }
				}
			}
		}
	}

  for ( s = 0; s < usedtolbsp; s++ )
	{	gdxdphi[s] = gdxdOphi[s];
	  //printf("xdgdphi%f ", gdxdphi[s]);
	}

  free(gdxdOphi); gdxdOphi = NULL;
}



/*******************************************************************************

********************************************************************************/
/*float Reconstruction::PickfromMat(unsigned long ii,unsigned long jj)
{
 FILE*  fp;
 float djj,dii,size,dmn;
 int    j,mn, in;
 unsigned long i;

 char   filename[50];

 djj   = (float)jj;
 dii   = (float)ii;

 size  = (float)djj*tolbsp + dii ;

 dmn   = size/msize;
 mn    = dmn;
 
 in    = size - mn*msize;
 printf("\nsize=%f ", size);

  printf("ii=%u, jj=%u ,djj = %f dii = %f size = %f dmn = %f mn = %d in = %d", ii,jj,djj, dii, size,dmn,mn, in);

  if( dmn - mn == 0.0 )
    {
	  sprintf(filename, "matrix%c%d.txt", ch,mn);
     fp = fopen(filename,"r");
	 printf("\n open the %s \n ", filename );
	 // printf("ii=%u, jj=%u ,djj = %f dii = %f size = %f dmn = %f mn = %d in = %d", ii,jj,djj, dii, size,dmn,mn, in);
	 //getchar();
     for ( i = 0; i < msize; i++)
	   {
       fscanf(fp, "%f  ", &Mat[i]);  
	   // printf("\nMat=%f ",Mat[i]);getchar();
	   }

     fclose(fp);
    }
 
return Mat[in];
}
*/

/********************************************************************************
Descriptions:
     

*********************************************************************************/
 /*void Reconstruction::InserttoMatrix(unsigned long i, unsigned long j,float res)
{
 
 fprintf(fp, "%f  ", res);

}
 */


/********************************************************************


*********************************************************************/
/*void Reconstruction::writetofile(float *Mat, int  size,char* filename)
{
unsigned long i;
FILE* fp;

fp = fopen(filename, "w");

//fprintf(fp, "%d \n", size);

for ( i = 0; i < size; i++)
   fprintf(fp, "%f  ", Mat[i]);


fclose(fp);


}

*/



float Reconstruction::EvaluateFFT_GdXdPhiIntegral(int i1, int j1, int k1, fftw_complex *slice, fftw_complex *B_k, Oimage *image)
{

  int      i,j, ii1, jj1, ii, jj, jj3, half;
  Views* v;
  float    sum = 0.0, a, b, t, area, integral, iB, jB, kB;
  float   rotmat[9],translate[3]= {0.0,0.0,0.0}, X[3];
  fftw_complex *gd_i, *oneslice;//, *Xdf;
  int     isize, isize3;
  int p,q,r, s;


  // image->ux = PI2/(ImgNx+0.0); image->uy = PI2/(ImgNy+0.0);

  area = image->ux*image->uy;
  iB   = i1*Bscale;
  jB   = j1*Bscale;
  kB   = k1*Bscale;

   half = ImgNx/2;

  a    = -half;
  b    = half;
 
  for ( i=0; i < nv; i++ ) 
   {
	 isize  = i * img2dsize;
	 isize3 = i * img2dsize * 3;

	 //for ( j = 0; j < 9; j++ ) 
	 // rotmat[j] = Rmat[i*9 + j];

	 /*
	 BsplineCentralSlice(i1, j1, k1,i, rotmat, slice, B_k,1);  //Evaluate X_d(\phi_pqr).

	 	 
	 for ( j = 0; j < img2dsize; j++ )
	   {
		 image->data[j] = gd_FFT[isize+j][0] * slice[j][0] + gd_FFT[isize+j][1] * slice[j][1];
		 //printf("%f ", image->data[j]);
		 
	   }
	 */

	 integral = 0.0;

	 for ( ii1 = -half; ii1 <= half; ii1++ )
	   {
		 ii = (ii1 + half) * (ImgNx+1);
		 
		 for ( jj1 = -half; jj1 <= half; jj1++ )
		   {
			 jj = ii + jj1 + half;
			 
			 jj3 = jj + jj + jj;
			 
			 X[0] = slicemat[isize3 + jj3 + 0];
			 X[1] = slicemat[isize3 + jj3 + 1];
			 X[2] = slicemat[isize3 + jj3 + 2];	
			 
			 t  = iB*X[0] + jB*X[1] + kB*X[2];
			 
			 image->data[jj] = gd_FFT[isize+jj][0] * cos(t) - gd_FFT[isize+jj][1] * sin(t); //6-14-09
			 //printf("%f ", image->data[j]);
			 //integral = integral + image->data[jj];
		   }
	   }
	 //sum = sum + integral; //SimpsonIntegrationOn2DImage(image); 
	   sum = sum + SimpsonIntegrationOn2DImage(image);
	 
   }
  //printf("sum=%f ", sum);
  //sum = 1.0/nv * sum*area;
  sum = 1.0/(nv+0.0) * sum;
  //getchar(); 
  
  return sum;
  

}




/********************************************************************************/
Oimage* Reconstruction::GdiFromGd(Oimage* gd, int i)
{

int i1,j1,ii ;

Oimage* gd_i    = InitImageParameters(1, N-1, N-1, 1, 1);
//gd_i->image->ox = gd->image->ox;
//gd_i->image->oy = gd->image->oy;
//gd_i->image->oz = gd->image->oz;


float*                  newdata = (float *) gd_i->data;
float                  background = 0.0;

for ( i1 = 0; i1 < N-1; i1++ )
   for ( j1 = 0; j1 < N-1; j1++ )
       {
        ii   =  i1*(N-1) + j1;

        newdata[ii] = gd->data[i*img2dsize+ ii]; 
       }

return gd_i;
}


void Reconstruction::FFT_gdMatrix()
{
  int i,j, i1, j1;
  int    isize, size;
  fftw_complex *in, *out;


  in = (fftw_complex *)fftw_malloc(img2dsize * sizeof(fftw_complex));
									
  out = (fftw_complex *)fftw_malloc(img2dsize * sizeof(fftw_complex));
  gd_FFT = (fftw_complex *)fftw_malloc(nv * img2dsize * sizeof(fftw_complex));

 									

  for ( i = 0; i < nv * img2dsize; i++ )
	{
	  gd_FFT[i][0] = 0.0;
	  gd_FFT[i][1] = 0.0;
	}



  for ( i = 0; i < nv; i++ )
	{
	  isize = i * img2dsize;
	  for ( j = 0; j < img2dsize; j++ )
		{
		  in[j][0] = 0.0; in[j][1] = 0.0;
		  out[j][0] = 0.0; out[j][1] = 0.0; 
		 
		}

	   FFT_gdi(i, in, out);

	  for ( j = 0; j < img2dsize; j++ )
		{
		  gd_FFT[isize+j][0] = out[j][0];
		  gd_FFT[isize+j][1] = out[j][1];
		  //printf("%f     ",sqrt(out[j][0]*out[j][0]+out[j][1]*out[j][1])); 
          //printf("\ngdFFT=%f %f  ", gd_FFT[isize+j][0], gd_FFT[isize+j][1]);
		}
	  //fftw_free(out);
	  //	  out = NULL;

	}

  //getchar();
  fftw_free(in);     in = NULL;
  fftw_free(out);   out = NULL;

}

void Reconstruction::FFT_gdi(int i, fftw_complex * in, fftw_complex *out)
{

  int i1,j1,ii, jj, half;
  fftw_plan fft_Gd;
  int    isize, imgnx;
  float wk, r1, r2;

  imgnx = ImgNx;
 

  isize = i * img2dsize;
  half  = ImgNx/2;
 
  for ( i1 = 0;  i1 <= ImgNx; i1++ )
	{
	  ii = i1 * (ImgNx+1);
	  
	  for ( j1 = 0; j1 <= ImgNx; j1++ )
		{
		  
		  jj =  ii + j1;

		  //wk =  PI2*(i1*0.5+j1*0.5);// matlab fftshift();

		  wk =  M_PI*(i1+j1)*ImgNx/(ImgNx+1.0);  // 6.10-09
          //wk = 0.0;
		  in[jj][0] = gd->data[isize + jj] * cos(wk);
		  in[jj][1] = gd->data[isize + jj] * sin(wk);
		} 
	}
  
  //fft2D_shift(out, ImgNx+1, ImgNx+1);
  fft_Gd =  fftw_plan_dft_2d(ImgNx+1, ImgNx+1, in, out,FFTW_FORWARD,FFTW_ESTIMATE);
  
  fftw_execute(fft_Gd);
  
  //fft2D_shift(out, ImgNx+1, ImgNx+1);


  //modified 6.15.09.
  //for ( i1 = 0;  i1 <= ImgNx; i1++ )
  for ( i1 = -half;  i1 <= half; i1++ )
	{
	  ii = (i1+half) * (ImgNx+1);
	  //printf("\n");	  
	  //	  for ( j1 = 0; j1 <= ImgNy; j1++ )
	  for ( j1 = -half; j1 <= half; j1++ )
		{
		  
		  jj =  ii + j1+half;
		  //wk =  PI2*(i1*0.5+j1*0.5);
		  wk =  M_PI*(i1+j1)*ImgNx/(ImgNx+1.0);// matlab fftshift();
		  //wk = 0.0;
		  r1 = cos(wk)*out[jj][0]-sin(wk)*out[jj][1];
		  r2 = cos(wk)*out[jj][1]+sin(wk)*out[jj][0];

		  out[jj][0] = r1;
		  out[jj][1] = r2;
		  //printf("%f %fi ", out[ii+j1][0], out[ii+j1][1]);
		}
	}
  
  fftw_destroy_plan(fft_Gd);


  //getchar();

}



void Reconstruction::Convertgd2BsplineCoefs(float *gddata)
{
  int i, j, v, id, iv, N2, iv1;
  N2 = (usedN+1) * (usedN+1);
  gd_coefs = (float *)malloc(nv * N2 *sizeof(float));

  for ( v = 0; v < nv; v++ )
	{
	  iv  = v * N2;
	  iv1 = v * img2dsize;

	  for ( i = 2; i <= N-2; i++)
		for ( j = 2; j <= N-2; j++ )
		  {
			id = i * (ImgNx+1) + j;
			gd_coefs[iv + (i-2)*(usedN+1)+j-2] = gddata[iv1+id];
		  }
	  bspline->ConvertToInterpolationCoefficients_2D(gd_coefs+iv, usedN+1, usedN+1,CVC_DBL_EPSILON);
	}

  /*
  float *data = (float*)malloc(img2dsize * sizeof(float));
  int i1, j1, ii;
  float s1, s2;

  for ( v = 0; v < nv; v++ )
	{
	  iv  = v * N2;
	  iv1 = v * img2dsize;
	  for ( i = 0; i < img2dsize; i++ ) data[i] = 0.0; 

	  for ( i = 0; i <= ImgNx; i++ )
		for ( j = 0; j <= ImgNx; j++ )
		  {
			for (i1 = BStartX; i1 <= BFinishX; i1++ )
			  {
				for ( j1 = BStartX; j1 <= BFinishX; j1++ )
				  {
					ii = (i1 - BStartX)*(usedN+1) + j1 - BStartX;
					s1 = bspline->Cubic_Bspline_Interpo_kernel(i-ImgNx/2.0, i1);
					s2 = bspline->Cubic_Bspline_Interpo_kernel(j-ImgNx/2.0, j1);
					
					data[i*(ImgNx+1)+j] += gd_coefs[iv+ii] * s1 * s2;

				  }
			  }
		  
			printf("\ni j =%d %d cha=======%e ", i, j, data[i*(ImgNx+1)+j]- gd->data[iv1+i*(ImgNx+1)+j]); 
		  }
	}
  */

}


/************************************************************************/
float* Reconstruction::EvaluateBsplineBaseFT()
{
  float         R, delR, s, X, B, *B_k;

  fftw_complex  *in, *out;

  fftw_plan     ifft_BsplineFT, fft_Bspline;

  int           i, j,  N3;


  //bspline->GramSchmidtofBsplineBaseFunctions(N,M);

  // bspline->Evaluate_OrthoBSpline_Basis_AtGaussNodes();


  N3 =  3*(N+1);


  
  in = (fftw_complex *)fftw_malloc(ImgNx*sizeof(fftw_complex));
  out= (fftw_complex *)fftw_malloc(ImgNx*sizeof(fftw_complex));

  BsplineFT = (fftw_complex *)fftw_malloc((N+1)*ImgNx*sizeof(fftw_complex));
  
  for ( i = 0; i <= N; i++ )
	{
   for ( j = 0; j < ImgNx; j++ )
       {
		 in[j][0] = bspline->BernBaseGrid[j*N3+3*i];//*cos(PI*(ImgNx-1.0)*j/ImgNx);
		//printf("\nBernbaseGrid=%f ", bspline->BernBaseGrid[j*N3+3*i]);
		 in[j][1] =0.0;// bspline->BernBaseGrid[j*N3+3*i]*sin(PI*(ImgNx-1.0)*j/ImgNx);
	   }
   //   getchar();
   //	for ( j = 0; j <= M; j++ )
   //   printf("\n in=%f %f ", in[j][0], in[j][1]);


    fft_Bspline = fftw_plan_dft_1d(ImgNx,in,out,FFTW_FORWARD,FFTW_ESTIMATE);

  fftw_execute(fft_Bspline);

  fft1D_shift(out ,ImgNx);

	for ( j = 0; j < ImgNx; j++ )
	  {
		BsplineFT[i*ImgNx+j][0] = out[j][0];
		BsplineFT[i*ImgNx+j][1] = out[j][1];


		// printf("\nnew out=%f %f ", out[j][0], out[j][1]);

	  }
    //getchar();

    fftw_destroy_plan(fft_Bspline);
	}

  fftw_free(in);   in = NULL;
  fftw_free(out);  out = NULL;



  return NULL;
}
// p , q, r are alpha , beta, gamma in paper.
/************************************************************************/
void Reconstruction::Evaluate_BsplineFT_ijk_pqr()
{
  int a, b, c, i, j, k, p, q, r, v, s, size;
 
  float centers[3], centers1[3], rotmat[9], x[3], x1[3], Bradius, dist, sum, ijk[3], pqr[3];
 
  int    ss, ii, jj, iii, jjj, kkk, id, jd, jjjj, N2, tol;
 char filename[50];
 Oimage *oneproj = NULL, *image=NULL;
 fftw_complex *oneproj_j, *slice, *B_k;
 int Begin, End;
 float e1[3], e2[3], d[3], scale;

 scale = PI2/(ImgNx+1);

 N2 = (2*usedN+1) *(2*usedN+1);

 slice=(fftw_complex *)fftw_malloc(img2dsize*sizeof(fftw_complex));
B_k = (fftw_complex *)fftw_malloc(sizeof(fftw_complex));

 image  = InitImageParameters(1, ImgNx+1, ImgNx+1, 1, 1);
 image->ux = scale;
 image->uy = scale;
 image->uz = 1.0;
 image->nx = ImgNx+1;
 image->ny = ImgNx+1;
 image->nz = 1;



for ( i = 0; i < img2dsize; i++ )
   {slice[i][0] = 0.0;slice[i][1] = 0.0;}
  
 delx = Bscale;

  //delx = 1.0;

  Bradius = delx * 2.0 * sqrt(3.0);

  Bradius = 4.0*Bradius*Bradius;
  
  Mat = (float *)malloc((2*usedN+1) * (2*usedN+1) *(2*usedN+1)*sizeof(float));


  Begin = BStartX - BFinishX;
  End   = BFinishX - BStartX;
  
  printf("\nBegin===========%d  End=%d", Begin, End); 
    for ( a = Begin ; a <= End; a++ )
	{

	  iii = (a-Begin) *N2;
	  //if(a <= 0) {p = BStartX; i = p - a;}
	  //else { i = BStartX; p = a + i;} 

	  //ijk[0] = i * Bscale;
	  //pqr[0] = p * Bscale;

	  for ( b = Begin; b <= End; b++ )
	  {

		jjj = (b-Begin) * (2*usedN+1);
		//if(b <= 0) {q = BStartX; j = q - b;}
		//else { j = BStartX; q = b + j;} 

		//ijk[1] = j * Bscale;
		//pqr[1] = q * Bscale;

		for ( c = Begin; c <= End; c++ )
		  {

			jj = iii + jjj + c - Begin;
			//if(c <= 0) {r = BStartX; k = r - c;}
			//else { k = BStartX; r = c + k;} 
			//ijk[2] = k * Bscale;
			//pqr[2] = r * Bscale;


			sum = 0.0;
			for ( v=0; v < nv;  v++ )
			  {
				for ( s = 0; s < 9; s++ ) rotmat[s] = Rmat[v*9 + s];
				for ( s = 0; s < 3; s++ )
				  {
					e1[s] = rotmat[s];
					e2[s] = rotmat[3+s];
					d[s]  = rotmat[6+s];
					
				  }

				//printf("\na b c = %d %d %d ijk=%f %f %f pqr=%f %f %f ", a, b, c, ijk[0], ijk[1], ijk[2], pqr[0], pqr[1], pqr[2]);		
				
				//ProjectPoint2Plane(ijk, e1, e2);
				//ProjectPoint2Plane(pqr, e1, e2);

				//dist = (ijk[0] - pqr[0])* (ijk[0] - pqr[0])+(ijk[1]-pqr[1])*(ijk[1]-pqr[1]);
//if (dist- Bradius != 0) printf("\n dist =====================  %e,  %e \n", dist,  Bradius); 

				//if(dist <=  Bradius )
				//{

				BsplineCentralSlice(a, b, c, v, rotmat, slice, B_k, 2);
				
				for ( s = 0; s < img2dsize; s++) 
				  image->data[s] = slice[s][0];
				sum = sum + SimpsonIntegrationOn2DImage(image);
			  
                                //sum += BsplineCentralSlice1(a, b, c, rotmat); //new algorithm.

 	
			  }
			Mat[jj] = 1.0/(nv+0.0)*sum;
			//printf("\nMat=%f ",Mat[jj]);if(jj==0) getchar();	
		  } 
	  }
	}    
  
	/*
  
  int i1, j1, half, isize;
  half = ImgNx/2;
  Mat = (float *)malloc((usedN+1) * (usedN+1) *(usedN+1)*(usedN+1) * (usedN+1)*(usedN+1)*sizeof(float));



  fftw_complex *B_a = (fftw_complex *)fftw_malloc(sizeof(fftw_complex));
  for ( i = BStartX; i<=BFinishX; i++ )
	for ( j = BStartX; j <= BFinishX; j++ )
	  for ( k = BStartX; k <= BFinishX; k++ )
		{
		  for ( a = BStartX; a <= BFinishX; a++ ) 
			for ( b = BStartX; b <= BFinishX; b++ )
			  for ( c = BStartX; c <= BFinishX; c++ )
				{
				  sum = 0.0;
				  jjj = (i-BStartX)*(usedN+1)*(usedN+1)*(usedN+1)*(usedN+1)*(usedN+1)
					+(j-BStartX)*(usedN+1)*(usedN+1)*(usedN+1)*(usedN+1)+
					(k-BStartX)*(usedN+1)*(usedN+1)*(usedN+1)
					+ (a-BStartX)*(usedN+1)*(usedN+1) + 
					(b-BStartX)*(usedN+1) + c-BStartX;
				  for ( v = 0; v < nv; v++ )
					{
					  isize = v * img2dsize * 3;

					  for ( s = 0; s < 9; s++ ) rotmat[s] = Rmat[v*9+s];
					  for ( s = 0; s < 3; s++ )
						{
						  e1[s] = rotmat[s]; 
						  e2[s] = rotmat[3+s];
						  d[s]  = rotmat[6+s];
						}

					  
					  for ( i1 = -half; i1 <= half; i1++ )
						for ( j1 = -half; j1 <= half; j1++ )
						  {
							jj = (i1+half)*(ImgNx+1) + j1+half;
							x[0] = slicemat[isize+jj*3+0];
							x[1] = slicemat[isize+jj*3+1];
							x[2] = slicemat[isize+jj*3+2];
   							VolBsplineBaseFT(i,j, k, x[0], x[1], x[2], Bscale, B_k);
							VolBsplineBaseFT(a,b, c, x[0], x[1], x[2], Bscale, B_a);
							slice[jj][0] = B_k[0][0] * B_a[0][0] + B_k[0][1] * B_a[0][1];
							slice[jj][1] = B_k[0][0] * B_a[0][1] - B_k[0][1] * B_a[0][0];

						  }
					  for ( s = 0; s < img2dsize; s++ )
						image->data[s] = slice[s][0];
					  image->ux = PI2/ImgNx; 
					  image->uy = PI2/ImgNx;
					  sum += SimpsonIntegrationOn2DImage(image);
					}
				  Mat[jjj] = 1.0/nv * sum;
				}
		}


		fftw_free(B_a);       B_a  = NULL;
  


	*/




  fftw_free(slice);      slice = NULL;
  fftw_free(B_k);        B_k   = NULL;

  kill_all_but_main_img(image);    
  free(image);           image = NULL;
  
}



void  Reconstruction::BsplineCentralSlice(int i, int j, int k, int iv, float rotmat[9], fftw_complex *slice, fftw_complex *B_k, int index) 
 {

   int i1, j1, ii, jj, half;
   float oldx[3],  X[3], e1[3], e2[3], d[3], a, b, wk_x, wk_y;
   int    isize, jj3;

   
    half = ImgNx/2;
   //half = FN/2;
/*
   a    = -half+0.0;
   b    = half+0.0;
*/


  isize = iv * img2dsize * 3;

 for ( i1 = -half; i1 <= half; i1++ )
	{
	  ii = (i1 + half) * (ImgNx+1);

	for ( j1 = -half; j1 <= half; j1++ )
	  {
        jj = ii + j1 + half;

		jj3 = jj + jj + jj;
	
		X[0] = slicemat[isize + jj3 + 0];
		X[1] = slicemat[isize + jj3 + 1];
		X[2] = slicemat[isize + jj3 + 2];
	
		//VolBsplineBaseFT(i, j, k, X[0], X[1], X[2], Bscale, B_k);//, index);
		VolBsplineBaseFT2(i, j, k, X[0], X[1], X[2], iv, jj, B_k, index);
	    slice[jj][0] = B_k[0][0];  
        slice[jj][1] = B_k[0][1];  
		//if(i == 1 && j == 1 && k == 1) 
		// printf("\ni j k = %d %d %d slice = %f %f ", i, j, k,B_k[0][0],  B_k[0][1]);
		// fftw_free(B_k); B_k = NULL;
	  }
	}

 }


void Reconstruction::GetCentralSlice(const fftw_complex *Vol, float *COEFS, int sub)
{
  int i, j, k, s, t, v, ii, jj, ix, iy, iz, nx, ny, nz, nynz, half, lx, ly, rx, ry, id, id1;
  float rotmat[9], e1[3], e2[3], wk_x, wk_y, X[3], area, ftscale, wk, rnv;
  float xd, yd, zd, real, im, r1, r2, r3, r4, r5, r6, r7, r8, sum, value, *result;
  int    proj_length, proj_size, N2, scale, proj_imgnx;
  fftw_complex *slice, *out;
  fftw_plan     invfft_slice;
  int subImgNx, subhalf;

  subImgNx = sub * ImgNx;
  
  nx = subImgNx;
  ny = subImgNx;
  nz = subImgNx;
  nynz = ny * nz;
  half = ImgNx/2;
  subhalf = subImgNx/2;
  //  printf("\nnx ny nz subhalf=%d %d %d %d ", nx, ny, nz, subhalf);
  scale = (int)Bscale;
  ////ftscale = PI2/(ImgNx+1.0);
  //ftscale = PI2/ImgNx;
  ftscale    = 1.0/ImgNx;
  //ftscale      = 1.0;

  proj_length = 8*scale;
  proj_imgnx  = 8 * scale -1;

  proj_size   = proj_length * proj_length;
  N2          = (usedN+1)*(usedN+1);
  area        = 1.0*1.0;
  rnv         = 1.0/nv;

  slice = (fftw_complex*)fftw_malloc(subImgNx*subImgNx*sizeof(fftw_complex));
  out   = (fftw_complex*)fftw_malloc(subImgNx*subImgNx*sizeof(fftw_complex));
  result =(float*)malloc(ImgNx*ImgNx*sizeof(float));

  //xdgdphi = (float *)malloc(usedtolbsp*sizeof(float));
 

  /*  
//compute exact fourier volume.
int i1, j1, k1, iii, subVolImgSize = nx * ny * nz;
float X1[3];
fftw_complex *B_k;
B_k = (fftw_complex*)fftw_malloc(sizeof(fftw_complex));
//ftscale = PI2/(ImgNx+1);
//subhalf = subImgNx/2;

fftw_complex *eVol = (fftw_complex*)fftw_malloc(subVolImgSize*sizeof(fftw_complex));
for ( i = 0; i < subVolImgSize; i++ ) {eVol[i][0] = 0.0; eVol[i][1] = 0.0;}

for(i = -subhalf; i <= subhalf; i++ )
   for ( j = -subhalf; j <= subhalf; j++ )
      for ( k = -subhalf; k <= subhalf; k++ )
          {
           id    = (i+subhalf)*nynz + (j+subhalf)*(subImgNx+1) + k + subhalf;
           X1[0] =  i*ftscale;
           X1[1] =  j*ftscale;
           X1[2] =  k*ftscale;

           for (i1 = BStartX; i1 <= BFinishX; i1++ )
             {
              for ( j1 = BStartX; j1 <= BFinishX; j1++ )
                 {
                  for ( k1 = BStartX; k1 <= BFinishX; k1++ )
                     {
                      iii = (i1-BStartX)*(usedN+1) * (usedN+1) + (j1-BStartX)*(usedN+1) + k1-BStartX;
                      VolBsplineBaseFT(i1, j1, k1, X1[0], X1[1], X1[2], Bscale, B_k);

                      eVol[id][0] += COEFS[iii]*B_k[0][0];
                      eVol[id][1] += COEFS[iii]*B_k[0][1];
                      }
                 }
              }
                    
            if(fabs(eVol[id][0]-Vol[id][0]) >SMALLFLOAT)  printf("\ni j k =%d %d %d exact =%f %f    FTVol=%f %f", i, j, k, eVol[id][0], eVol[id][1], Vol[id][0], Vol[id][1]);
                     
           }

fftw_free(eVol); eVol = NULL;
fftw_free(B_k);  B_k = NULL;

getchar();

//end.
*/

  for (i = 0; i < usedtolbsp; i++ ) xdgdphi[i] = 0.0;
 
 for ( v = 0; v < nv; v++ )
	{
	  for ( i = 0; i < 9; i++ )
		rotmat[i] = Rmat[v*9+i];

	  for (i = 0; i <3; i++ )
		{
		  e1[i] = rotmat[i];
		  e2[i] = rotmat[3+i];

		  if(fabs(e1[i])<SMALLFLOAT) e1[i] = 0.0;
		  if(fabs(e2[i])<SMALLFLOAT) e2[i] = 0.0;


		  //printf("\ne1 e2 d=%1.15e %1.15e %1.15e ", e1[i], e2[i], rotmat[6+i]);
		}
	  
	  for ( i = 0; i < subImgNx*subImgNx; i++ ) 
		{
		  slice[i][0] = 0.0;
		  slice[i][1] = 0.0;
		} 

	  for ( i = -subhalf; i < subhalf; i++ )
		for ( j = -subhalf; j < subhalf; j++ )
		  {
			ii = (i+subhalf)*subImgNx+j+subhalf;
			
			X[0] = i * e1[0] + j * e2[0];
			X[1] = i * e1[1] + j * e2[1];
			X[2] = i * e1[2] + j * e2[2];
			 //printf("\ni j=%d %d  X=%2.16e %1.16e %2.16e ", i, j, X[0], X[1], X[2]);
			/*
			//Nearest neigbor interpolation.
	        ix = floor(X[0]+0.5);
			iy = floor(X[1]+0.5);
			iz = floor(X[2]+0.5);
			//printf("\nix iy iz=%d %d %d ", ix, iy, iz);

			
			ix = ix + subhalf; iy = iy + subhalf; iz = iz + subhalf;

			real = im = 0.0;

			if( ix >= 0 && ix < ImgNx && iy >= 0 && iy < ImgNx && iz >= 0 && iz < ImgNx)
			  {
				real = Vol[ix*ImgNx*ImgNx+iy*ImgNx+iz][0];
				im   = Vol[ix*ImgNx*ImgNx+iy*ImgNx+iz][1];

			  }
			*/
			//ix   = (fabs(X[0])<SMALLFLOAT)?0:(int)floor(X[0]); 
			//iy   = (fabs(X[1])<SMALLFLOAT)?0:(int)floor(X[1]);
			//iz   = (fabs(X[2])<SMALLFLOAT)?0:(int)floor(X[2]);
			//printf("\nX=%f %f %f ", X[0], X[1], X[2]);

			ix     = (int)floor(X[0]);
			iy     = (int)floor(X[1]);
			iz     = (int)floor(X[2]);

                        //printf("ix iy iz =%d %d %d half=%d ", ix, iy, iz, half);

			xd   = (X[0] - ix);
			yd   = (X[1] - iy);
			zd   = (X[2] - iz);
			//printf("\nxd yd zd=%1.15e %1.15e %1.15e ", xd, yd, zd);            

			ix = ix + subhalf; iy = iy + subhalf; iz = iz + subhalf;

                        //if(v==3) {ix -=1; iy-=1;}
                         //printf("ix iy iz=%d %d %d ", ix, iy, iz);
			//if(ix < 0 || iy < 0 || iz < 0 ) 
			//{printf("\nix iy iz=%d %d %d", ix, iy, iz);getchar();}
                        //Trilinear interpolation.
			r1 = (ix<0 || iy<0 || iz<0 || ix > nx-1 || iy > nx-1 || iz > nx-1)?0.0:Vol[ix    * nynz + iy    * ny+iz][0];
			r2 = (ix<-1|| iy<0 || iz<0 || ix > nx-2 || iy > nx-1 || iz > nx-1)?0.0:Vol[(ix+1)* nynz + iy    * ny+iz][0];
			r3 = (ix<0 || iy<-1|| iz<0 || ix > nx-1 || iy > nx-2 || iz > nx-1)?0.0:Vol[ix    * nynz + (iy+1)* ny+iz][0];
			r4 = (ix<-1|| iy<-1|| iz<0 || ix > nx-2 || iy > nx-2 || iz > nx-1)?0.0:Vol[(ix+1)* nynz + (iy+1)* ny+iz][0];
			r5 = (ix<0 || iy<0 || iz<-1|| ix > nx-1 || iy > nx-1 || iz > nx-2)?0.0:Vol[ix    * nynz + iy    * ny+iz+1][0];
			r6 = (ix<-1|| iy<0 || iz<-1|| ix > nx-2 || iy > nx-1 || iz > nx-2)?0.0:Vol[(ix+1)* nynz + iy    * ny+iz+1][0];
			r7 = (ix<0 || iy<-1|| iz<-1|| ix > nx-1 || iy > nx-2 || iz > nx-2)?0.0:Vol[ix    * nynz + (iy+1)* ny+iz+1][0];
			r8 = (ix<-1|| iy<-1|| iz<-1|| ix > nx-2 || iy > nx-2 || iz > nx-2)?0.0:Vol[(ix+1)* nynz + (iy+1)* ny+iz+1][0];
                        //printf("\nr1 r2 r3 r4 r5 r6 r7 r8=%e %e %e %e %e %e %e %e ", r1, r2, r3, r4, r5, r6, r7, r8);
			real = TrilinearInterpolation8(xd, yd, zd,r1, r2, r3, r4, r5, r6, r7, r8);

			//if(fabs(real-r1) > SMALLFLOAT ) printf("\nreal-r1=%e ", real-r1);												

			r1 = (ix<0 || iy<0 || iz<0 || ix > nx-1 || iy > nx-1 || iz > nx-1)?0.0:Vol[ix    * nynz + iy    * ny+iz][1];
			r2 = (ix<-1|| iy<0 || iz<0 || ix > nx-2 || iy > nx-1 || iz > nx-1)?0.0:Vol[(ix+1)* nynz + iy    * ny+iz][1];
			r3 = (ix<0 || iy<-1|| iz<0 || ix > nx-1 || iy > nx-2 || iz > nx-1)?0.0:Vol[ix    * nynz + (iy+1)* ny+iz][1];
			r4 = (ix<-1|| iy<-1|| iz<0 || ix > nx-2 || iy > nx-2 || iz > nx-1)?0.0:Vol[(ix+1)* nynz + (iy+1)* ny+iz][1];
			r5 = (ix<0 || iy<0 || iz<-1|| ix > nx-1 || iy > nx-1 || iz > nx-2)?0.0:Vol[ix    * nynz + iy    * ny+iz+1][1];
			r6 = (ix<-1|| iy<0 || iz<-1|| ix > nx-2 || iy > nx-1 || iz > nx-2)?0.0:Vol[(ix+1)* nynz + iy    * ny+iz+1][1];
			r7 = (ix<0 || iy<-1|| iz<-1|| ix > nx-1 || iy > nx-2 || iz > nx-2)?0.0:Vol[ix    * nynz + (iy+1)* ny+iz+1][1];
			r8 = (ix<-1|| iy<-1|| iz<-1|| ix > nx-2 || iy > nx-2 || iz > nx-2)?0.0:Vol[(ix+1)* nynz + (iy+1)* ny+iz+1][1];				

			//printf("\ni1 i2 i3 i4 i5 i6 i7 i8=%e %e %e %e %e %e %e %e ", r1, r2, r3, r4, r5, r6, r7, r8);
			im = TrilinearInterpolation8(xd, yd, zd,r1, r2, r3, r4, r5, r6, r7, r8);
			//if(fabs(im-r1) > SMALLFLOAT ) printf("im-r1=%e ", im - r1);

			slice[ii][0] = real;
			slice[ii][1] = im ;
//printf("\nslice=%f %f ", slice[ii][0], slice[ii][1]);

/*
			if(v==3) 
                            {
                             slice[ii][0] = Vol[ix    * nynz + iy    * ny+3-iz][0]; 
                             slice[ii][1] = Vol[ix    * nynz + iy    * ny+3-iz][1];
                             printf("\nslice=%f %f real im=%f %f ", slice[ii][0], slice[ii][1], real, im);
                            }
*/
                        //printf("\niix iy iz=%d %d %d real im=%e %e xf yf zf=%e %e %e", ix, iy, iz, real, im, xd, yd, zd);

			//if( fabs(real) > 1000000 || fabs(im) > 10000) 
			 // {printf("\ni j = %d %d slice=%f %f gd=%f %f cha======%f %f ", i, j, real, im, gd_FFT[v*img2dsize+ii][0],gd_FFT[v*img2dsize+(i+half)*(ImgNx+1)+j+half][1], slice[ii][0], slice[ii][1]);
			//	getchar();
			 // }
		  }
	  //getchar();

	  for (i = -subhalf; i < subhalf; i++ )
		for ( j = -subhalf;  j < subhalf; j++ )
		  {
			ii = (i+subhalf)*subImgNx+j+subhalf;
			wk = -M_PI*(i+j);
			//wk = 0.0;

			real = cos(wk)*slice[ii][0] - sin(wk)*slice[ii][1];
			im   = cos(wk)*slice[ii][1] + sin(wk)*slice[ii][0]; 

			slice[ii][0] = real;
			slice[ii][1] = im;
			//printf("\ni j = %d %d slice=%f %f gd=%f %f ", i, j, real, im, gd_FFT[v*img2dsize+(i+half)*(ImgNx+1)+j+half][0],gd_FFT[v*img2dsize+(i+half)*(ImgNx+1)+j+half][1]);
          	        //printf("\ni j = %d %d slice=%f %f", i, j, real, im);

		  }
	  //getchar();

	  invfft_slice = fftw_plan_dft_2d(subImgNx,subImgNx,slice,out,FFTW_BACKWARD,FFTW_ESTIMATE);
	  fftw_execute(invfft_slice);
	  //out = 1/[(ImgNx+1)^2] * out.

	  for (ii = i = 0; i < subImgNx; i++ )
		for ( j = 0;  j < subImgNx; j++, ii++ )
		  {
			//ii = i*ImgNx+j;
			wk = -M_PI*(i+j);
			//wk = 0.0;
			
			real = cos(wk)*out[ii][0] - sin(wk)*out[ii][1];
			im   = cos(wk)*out[ii][1] + sin(wk)*out[ii][0];

			out[ii][0] = real;
			out[ii][1] = im;
 
		  }

	  //fftw_destroy_plan(invfft_slice);
	  //printf("\n\n"); 			  
/*
	    for ( ii = i = 0; i < ImgNx; i++ )
	   for ( j = 0; j < ImgNx; j++, ii++ )
	  	{
		  	   printf("\ni j = %d %d inverse=%f  gd=%f ", i, j, 1.0/(ImgNx*ImgNx)*out[ii][0], gd->data[v*img2dsize+i*(ImgNx+1)+j]);
	  	   printf("cha=========%e ",  1.0/(ImgNx*ImgNx)*out[ii][0] - gd->data[v*img2dsize+i*(ImgNx+1)+j]);
	  	   //printf("\ni j = %d %d inverse=%f  gd=%f ", i, j, 1.0/(ImgNx*ImgNx)*out[ii][0], gd->data[v*ImgNx*ImgNx+ii]);

	  	}
	  getchar();
*/	      
	  for ( jj = i = 0; i < ImgNx; i++ )
		for ( j = 0; j < ImgNx; j++, jj++ )
		{
			ii = sub*i* subImgNx + sub*j;
			
			result[jj] =  1.0/(subImgNx*subImgNx) * out[ii][0] - gd->data[v*img2dsize+i*(ImgNx+1)+j];
			//printf("\nsub=%d slice=%f  gd=%f", sub, 1.0/(subImgNx*subImgNx) * out[ii][0] ,gd->data[v*img2dsize+i*(ImgNx+1)+j]); 
		  }
	  //getchar();	  
	  
	  id  = v * usedtolbsp * 2;
      id1 = v * usedtolbsp*proj_size;
	  
	  for ( i = BStartX; i <= BFinishX; i++ )
	  	for ( j = BStartX; j <= BFinishX; j++ )
		  for( k = BStartX; k <= BFinishX; k++ )
			{ 
			  ii = (i-BStartX)*N2 + (j-BStartX)*(usedN+1)+k-BStartX;
			  lx = startp[id+2*ii+0];
			  ly = startp[id+2*ii+1];
			  rx = ((lx+proj_imgnx)>=half)?(half-1):(lx+proj_imgnx); 
			  ry = ((ly+proj_imgnx)>=half)?(half-1):(ly+proj_imgnx);
			  
			  sum = 0.0;
			  if(lx <-half || ly < -half ) 
				{printf("\nlx rx ly ry=%d %d %d %d are too small.", lx, rx, ly, ry);
				  getchar();}
			  for(s = lx; s <= rx; s++ )
				for(t = ly; t <= ry; t++ )
				  {
					//value = out[(s+half)*(ImgNx+1)+t+half][0] * proj_VolB[id1+ii*proj_size+(s-lx)*proj_length+t-ly];
					//value = out[(s+half)*ImgNx+t+half][0] * proj_VolB[id1+ii*proj_size+(s-lx)*proj_length+t-ly];
					value = result[(s+half)*ImgNx+t+half] * proj_VolB[id1+ii*proj_size+(s-lx)*proj_length+t-ly];
					
					//printf("\nout=%f  proj_img=%f ", out[(s+half)*(ImgNx+1)+t+half][0], proj_VolB[id1+ii*proj_size+(s-lx)*proj_length+t-ly]);
					if( (t == -half  &&  (s == -half || s == half-1)) || (t == half-1 &&(s == -half || s == half-1)))
					  sum += 0.25 * value;
					else if(s > -half && s < half-1 && t > -half && t < half-1)
					  sum += value;
					else
					  sum += 0.5 *value;
					/*                    if(fabs(sum)>100000000) 
					  {
					  printf("\nsum=%f out=%f proj=%f s t =%d %d lx ly rx ry=%d %d %d %d", sum,out[(s+half)*ImgNx+t+half][0] , proj_VolB[id1+ii*proj_size+(s-lx)*proj_length+t-ly], s, t, lx, ly, rx, ry);
					  getchar();
					  }
					*/
				  }
			  //printf("sum=%f ", sum);
			  xdgdphi[ii] += sum * area * rnv;
			  //printf("\nxgdphi=%f nv=%d", xdgdphi[ii], nv);
			}
			  
	}
 
  fftw_destroy_plan(invfft_slice);
  fftw_free(slice); fftw_free(out);
  slice = NULL;  out = NULL;
  free(result);
  //free(xdgdphi);
}





float Reconstruction::BsplineCentralSlice1(int i, int j, int k, float rotmat[9])
{
int i1, j1, s, thres, max_i, max_j, max_k, ri, rj, rk, temp, IMGNX, IMGNY, id;
float e1[3], e2[3], d[3], e1x, e2x, e1y, e2y, del_X, del_Y, ci, dj;
float u, v, A, X[3], sum, value, area, u1, v1, u2, v2, u3, v3, a[3], b[3], c[3], inn1, inn2, inn3, inn4;
fftw_complex *value1, *value2, *value3;

value1 = (fftw_complex*)malloc(sizeof(fftw_complex));
value2 = (fftw_complex*)malloc(sizeof(fftw_complex));
value3  = (fftw_complex*)malloc(sizeof(fftw_complex));

thres = 10;
// reorder i,j, k to max_i > max_j > max_k.
max_i = (i>j)?i:j;
max_j = (j>k)?j:k;
max_k = (temp=i<j?i:j)<k?temp:k;
if(max_i < max_j) {temp = max_i; max_i = max_j; max_j = temp;}
  
max_i = (max_i < thres)?thres:max_i;
max_j = (max_j < thres)?thres:max_j;
//printf("\nmax_i=%d max_j=%d max_k=%d i j k = %d %d %d ", max_i, max_j, max_k, i, j, k);
del_X = M_PI/(2.0 * max_i * Bscale);
del_Y = M_PI/(2.0 * max_j * Bscale);

//del_X = FW / (ImgNx1);
//del_Y = FW / (ImgNx+1);

IMGNX = FW / del_X;
IMGNY = FW / del_Y;
//IMGNX = ImgNx;
//IMGNY = ImgNx;

//printf("\nIMGNX =%d IMGNY=%d ", IMGNX, IMGNY);

for ( s = 0; s < 3; s++ )
   {
    e1[s] = rotmat[s];
    e2[s] = rotmat[3+s];
    d[s]  = rotmat[6+s];
   }
//printf("\ne1=%f %f %f e2=%f %f %f ", e1[0], e1[1], e1[2], e2[0], e2[1], e2[2]);

inn1 = sqrt(e1[0]*e1[0]+e1[1]*e1[1]);
inn2 = sqrt(e2[0]*e2[0]+e2[1]*e2[1]);

inn3 = sqrt(e1[1]*e1[1]+e1[2]*e1[2]);
inn4 = sqrt(e2[1]*e2[1]+e2[2]*e2[2]);


if(inn1 >= SMALLFLOAT && inn2 >= SMALLFLOAT && (e1[0]*e2[0]+e1[1]*e2[1])/(inn1 * inn2) != 1.0 ) 
   {
    e1x = e1[0]; e1y = e1[1]; e2x = e2[0]; e2y = e2[1]; 
    ri  = i; rj = j; rk = k;  id = 2;
   }
  
else if(inn3 >= SMALLFLOAT && inn4 >= SMALLFLOAT && (e1[1]*e2[1]+e1[2]*e2[2])/(inn3 * inn4) != 1.0 ) 
   {
    e1x = e1[1]; e1y = e1[2]; e2x = e2[1]; e2y = e2[2];
    ri  = j;  rj = k;  rk = i; id = 0;
   }

else {
      e1x = e1[0]; e1y = e1[2]; e2x = e2[0]; e2y = e2[2];
      ri  = i;  rj = k;  rk = j; id = 1;
     }
//if(ri != i || rj != j || rk != k ) {printf("\ne1x e1y=%e %e e2x e2y=%e %e  ri rj rk =%d %d %d ", e1x, e1y, e2x, e2y, ri, rj, rk);getchar();}
A  = e1x * e2y - e2x * e1y;

ci  = -IMGNX * del_X;
dj  = -IMGNY * del_Y;
//if(fabs(A) < 0.001)  printf("\nA=%e ", A);
u1  = (ci * e2y - dj * e2x) / A;
v1  = (dj * e1x - ci * e1y) / A;

dj  = (-IMGNY+1) * del_Y;

u2  = (ci * e2y - dj * e2x) / A;
v2  = (dj * e1x - ci * e1y) / A;

ci  = (-IMGNX+1) * del_X;
dj  = -IMGNY * del_Y;

u3  = (ci * e2y - dj * e2x) / A;
v3  = (dj * e1x - ci * e1y) / A;

a[0] = (u2-u1)*e1[0] + (v2-v1)*e2[0];
a[1] = (u2-u1)*e1[1] + (v2-v1)*e2[1];
a[2] = (u2-u1)*e1[2] + (v2-v1)*e2[2];
  
b[0] = (u3-u1)*e1[0] + (v3-v1)*e2[0];
b[1] = (u3-u1)*e1[1] + (v3-v1)*e2[1];
b[2] = (u3-u1)*e1[2] + (v3-v1)*e2[2];

c[0] = a[1]*b[2]-a[2]*b[1];
c[1] = a[2]*b[0]-a[0]*b[2];
c[2] = a[0]*b[1]-a[1]*b[0];

area = sqrt(c[0]*c[0]+c[1]*c[1]+c[2]*c[2]);
//printf("\ndel_X = %f del_Y=%f ", del_X, del_Y);
//printf("\nA = %f u1 v1=%f %f u2 v2=%f %f u3 v3=%f %f ", A, u1, v1, u2, v2, u3, v3);
//printf("\na=%f %f %f b=%f %f %f c=%f %f %f ", a[0], a[1], a[2], b[0], b[1], b[2], c[0], c[1], c[2]);
//printf("\narea=%f ", area);

sum = 0.0;
for ( i1 = -IMGNX; i1 <= IMGNX; i1++ )
  {
   ci     = i1 * del_X;
  
    BsplineFT1D(ri, ci, Bscale,value1);
    
   for ( j1 = -IMGNY; j1 <= IMGNY; j1++ )
      {
       dj = j1 * del_Y;
       BsplineFT1D(rj, dj, Bscale,value2);

       u  = (ci * e2y - dj * e2x) / A;
       v  = (dj * e1x - ci * e1y) / A;
/*
       X[0] = u * e1[0] + v * e2[0];
       X[1] = u * e1[1] + v * e2[1];
       X[2] = u * e1[2] + v * e2[2];
*/  
       X[id]  = u * e1[id] + v * e2[id];
//printf("\nci dj=%f %f X Y Z = %f %f %f ", ci, dj, X[0], X[1], X[2]);
       BsplineFT1D(rk, X[id], Bscale,value3);
       value = (value1[0][0] * value2[0][0] - value1[0][1] * value2[0][1] ) * value3[0][0] 
               - (value1[0][0] * value2[0][1] + value2[0][0]*value1[0][1]) * value3[0][1];

       if( (j1 == -IMGNY  &&  (i1 == -IMGNX || i1 == IMGNX)) || (j1 == IMGNY &&(i1 == -IMGNX || i1 == IMGNX)))
            sum += 0.25 * value;
         else if(i1 > -IMGNX && i1 < IMGNX && j1 > -IMGNY && j1 < IMGNY)
            sum += value;
         else
            sum += 0.5 *value;

      } 
  }

//getchar();
fftw_free(value1);
fftw_free(value2);
fftw_free(value3);

sum = area * sum;
return sum;

}



void Reconstruction::compare_Xdphi_Phi_pd(Views* nview)
{
  /*  
  int i,j,k,i1,j1,k1;
  float *xdB, *prjimg, rotmat[9];
  Oimage*   p, *Xdphi;
  fftw_complex *B_k=NULL, *in0, *out0, *in, *out, *IMAGE, *Phi_pd, *OUT0, *SLICE;
  fftw_plan     fft_Xdphi, fft_phi;
  int N3 = 3*(N+1), NX, NY;
  float e1[3],e2[3], d[3], X, Y, Z, r1, r2, s1, s2, max, MAX, wk_x, wk_y, wk;
  fftw_complex *slice;
  int a, b, ii, jj, v, dn; 
  float xa,xb, ya, yb, za,zb;

  // p = InitImageParameters(1, nx,ny,nz, nv);
  

  // NX = nx +2; NY = ny +2;
  //NX = (nx+2)*sqrt(3)/2;

  NX = ImgNx * sqrt(3)/2; //change nx+2 to nx.

  NX = 2*NX+2;
  NY = NX;
  MAX= 0.0;
  dn = NX/2-ImgNx/2;  

  in0 =(fftw_complex*)fftw_malloc(NX*NY*sizeof(fftw_complex));
  out0=(fftw_complex*)fftw_malloc(NX*NY*sizeof(fftw_complex));
  OUT0=(fftw_complex*)fftw_malloc(NX*NY*sizeof(fftw_complex));


  //IMAGE=(fftw_complex*)fftw_malloc((nx+1)*(ny+1)*sizeof(fftw_complex));
  //SLICE=(fftw_complex*)fftw_malloc((nx+1)*(ny+1)*sizeof(fftw_complex));
  
IMAGE=(fftw_complex*)fftw_malloc(ImgNx*ImgNy*sizeof(fftw_complex));
SLICE=(fftw_complex*)fftw_malloc(ImgNx*ImgNy*sizeof(fftw_complex));


  prjimg = (float*)malloc(NX*NY*sizeof(float));

  /// slice = (fftw_complex *)malloc((nx+1) * (ny+1) * sizeof(fftw_complex));
  slice = (fftw_complex *)malloc(NX * NY * sizeof(fftw_complex));
  B_k = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)); 


  //choose one bspline base.compute projection.
  for ( v = 0; v < nv; v++ )
	{
	  
	  for ( i1 = 0; i1 < 9; i1++ )
		{
			  rotmat[i1] = Rmat[9*v+i1];//printf("\nrotmat=%f ", rotmat[i1]);
		}
	  
	  
	 
	  printf("\nbegin new view v=%d ", v);
	  for ( i = BStartX; i <= BFinishX; i++ )
	  	for ( j = BStartX; j <= BFinishX; j++ )
		  for( k = BStartX; k <= BFinishX; k++ )
			{
  		  //  i = 0; j = 0; k = 0;
		  printf("\nbegin new spline base i j k =%d %d %d ", i,j,k);
		  bspline->Bspline_Projection(i, j, k,  rotmat,NX, 0.0, prjimg);
		  
		  for ( i1 = 0; i1 < NX; i1++ )
			for(j1 = 0; j1 < NY; j1++ )
			  {
				wk = PI2*(i1*0.5+j1*0.5);
				in0[i1*NY+j1][0] = prjimg[i1*NY+j1]*cos(wk);
				in0[i1*NY+j1][1] = prjimg[i1*NY+j1]*sin(wk);
			  }
		  //  getchar();
		  
		  //compute fft of projection.

		  // fft2D_shift(in0,NX,NY);

		  fft_Xdphi = fftw_plan_dft_2d(NX,NY,in0,out0,FFTW_FORWARD,FFTW_ESTIMATE);
		  
		  fftw_execute(fft_Xdphi);
		  ///fft2D_shift(out0,NX,NY);
		  
		  //printf("\n\nfft2dshift.");

		 
		  for( i1 = 0; i1 < NX; i1++ )
		  //for ( i1 = 1; i1 < nx+2; i1++ )
		  ///for ( i1 = 0; i1 < nx+2; i1++ )
			{//printf("\n");
			  for ( j1 = 0; j1 < NY; j1++ )
			 // for ( j1 = 1; j1 < ny+2; j1++ )
			  ///for ( j1 = 0; j1 < ny+2; j1++ )
				{
				  wk = PI2*(i1*0.5+j1*0.5);

				  r1 = cos(wk)*out0[i1*NX+j1][0]-sin(wk)*out0[i1*NY+j1][1];
				  r2 = cos(wk)*out0[i1*NX+j1][1]+sin(wk)*out0[i1*NY+j1][0];
				  out0[i1*NX+j1][0] = r1;
				  out0[i1*NX+j1][1] = r2;

				  printf("\n%f %f ",out0[i1*NX+j1][0],  out0[i1*NX+j1][1] );
				}
			}
		  getchar();

		  for ( i1 = 0; i1 < ImgNx; i1++)
			for ( j1 = 0; j1 < ImgNy; j1++ )
			  {

				IMAGE[i1*ImgNx+j1][0] = out0[(i1+dn)*NX+j1+dn][0];
				IMAGE[i1*ImgNx+j1][1] = out0[(i1+dn)*NX+j1+dn][1];
			  }
		  
  //  Evaluate central slice.
		  
		  ///a = -nx/2-1;
		  ///b =  nx/2+1;
		  a = -NX/2;
		  b = NX/2;


		  //printf("\na=%d b=%d ",a, b);
	
		  for ( i1 = 0; i1 < 3; i1++ )
			{
			  e1[i1] = rotmat[i1];
			  e2[i1] = rotmat[3+i1];  
			  d[i1]  = rotmat[6+i1];
			}
		  // printf("\ne1=%f %f %f e2=%f %f %f d=%f %f %f ", e1[0], e1[1], e1[2], e2[0], e2[1], e2[2], d[0], d[1], d[2]);
		  // printf("\n central slice.");
		  
		  max = 0.0;
		  
		  for ( i1 = -NX/2; i1 < NX/2; i1++ )
		  ///for ( i1 = -nx/2-1; i1 <= nx/2; i1++ )
		  //for ( i1 = 0; i1 <= nx+1; i1++ )
			{
			  // printf("\n");
			  for ( j1 = -NY/2; j1 < NY/2; j1++ )
			  /// for ( j1 = -ny/2-1; j1 <= ny/2; j1++ )
				//for ( j1 = 0; j1 <= ny+1; j1++ )
				{
				  
				  
				  wk_x = i1*PI2/(b-a);
				  wk_y = j1*PI2/(b-a);
				  
				  X = wk_x * e1[0] + wk_y * e2[0];
				  Y = wk_x * e1[1] + wk_y * e2[1];
				  Z = wk_x * e1[2] + wk_y * e2[2];
				  
				  //printf("\nX Y Z = %f %f %f ", X, Y,Z);	
				  VolBsplineBaseFT(i, j, k, X, Y, Z,Bscale, B_k);
				  


				  ii = (i1+NX/2)*NX+j1+NY/2;
				  ///ii = (i1+nx/2+1)*(ny+2)+j1+ny/2+1;
				  //ii = i1*NY+j1;

				
				
				  
				  	
				  slice[ii][0] = B_k[0][0];
				  slice[ii][1] = B_k[0][1];
				  
				
			
				  printf("\n%f %f ",slice[ii][0] ,slice[ii][1]); 
				  
				  
				}

			}
		  //fft2D_shift(slice,nx+1,ny+1);
		  getchar();


		  for ( i1 = 0; i1 < ImgNx; i1++)
			for ( j1 = 0; j1 < ImgNy; j1++ )
			  {

				SLICE[i1*ImgNx+j1][0] = slice[(i1+dn)*NX+j1+dn][0];
				SLICE[i1*ImgNx+j1][1] = slice[(i1+dn)*NX+j1+dn][1];
			  }
		 
		 
		  max = MaxError(IMAGE, SLICE,ImgNx*ImgNy);		  
		  //max = MaxError(out0, slice,NX*NY);
	
		  //max = MaxError_L2(out0, slice,NX*NY);

		  MAX = (max > MAX)?max:MAX;

		  // if(i==0 && j== 1  && k == 0) getchar();

		  printf("\nmax=%f i j k =%d %d %d v=%d", max, i,j,k,v);
		  //printf("\ne1=%f %f %f e2=%f %f %f d=%f %f %f", e1[0], e1[1], e1[2], e2[0], e2[1], e2[2], d[0], d[1], d[2]);

		  //if(max > 0.2)
		  //{
		  //   getchar();
		  //}
		}

	}

  printf("\nMAX = %f ", MAX); getchar();
  fftw_free(in0);
  fftw_free(out0);
  fftw_free(IMAGE);
  fftw_free(SLICE);
  fftw_free(OUT0);
  fftw_free(slice);
  fftw_free(B_k);
  free(prjimg);
*/
  	
	// 1D compare.
	/*
	fftw_complex *data1D, *outdata1D, *FTdata1D; 
  
  float x;long double values[6];
  
  B_k = (fftw_complex *)fftw_malloc(sizeof(fftw_complex));
 
		  // test 1D ok.
  		  	
printf("\nBspline base FT ");

 int I=4;
 int  n;
 float a, b, T, scal;
 a = 0.0; b = nx; 

 // a = -2.0; b = nx-2;
 n = nx ; T=(b-a)/n;

  scal = 0.8;
 //scal = 1.0;

 data1D = (fftw_complex*)malloc(n*sizeof(fftw_complex));
 outdata1D = (fftw_complex*)malloc(n*sizeof(fftw_complex));
 FTdata1D = (fftw_complex*)malloc(n*sizeof(fftw_complex));


 for ( i1 = 0; i1 < n; i1++ )
   // for ( i1 = 0; i1 < n; i1++ )
   {
	 //x =-2.0 + i1*4.0/nx;
	 //x = a +  i1-i*1.0;
	 
	 x = 1.0/scal*(a + i1*T)-I*1.0 ;
	 
	 bspline->Spline_N_Base(x,n-1,values);
	 // bspline->Spline_N_Base(x,n-1,values);
	 
	 data1D[i1][0] = (float)values[0]; //BernBaseGrid[i1*N3+3*i];
	 data1D[i1][1] = 0.0;
	 printf("\ndata1D=%f ", data1D[i1][0]);
   }
 
 fft_Xdphi = fftw_plan_dft_1d(n,data1D,outdata1D,FFTW_FORWARD,FFTW_ESTIMATE);
 
 //fft_Xdphi = fftw_plan_dft_1d(n,data1D,outdata1D,FFTW_FORWARD,FFTW_ESTIMATE);
 fftw_execute(fft_Xdphi);
 
 fft1D_shift(outdata1D ,n);
 
 //fft1D_shift(outdata1D ,n);
 
 for ( i1 = 0; i1 < n; i1++ )
   //for ( i1 = 0; i1 < n; i1++ )
   {
	 printf("\noutdata1D=%1.8f %1.8f ", outdata1D[i1][0],outdata1D[i1][1]);
   }
 
 printf("\n \n");
 float wk;float r1, r2,s1,s2, tem, max;
 max = 0.0;
 
 
 for ( i1 = 0; i1 < n; i1++ )
   {
	 wk = (i1-n/2) *1.0*PI2/(b-a);//ok. 
	 
	 BsplineFT1D(I, wk, scal, B_k);
	 
	 FTdata1D[i1][0] = 1.0/T*(cos(a*wk)*B_k[0][0]-sin(a*wk)*B_k[0][1]);
	 FTdata1D[i1][1] = 1.0/T*(cos(a*wk)*B_k[0][1]+sin(a*wk)*B_k[0][0]);
	 
	 printf("\nFTdata1D=%1.8f %1.8f ", FTdata1D[i1][0],FTdata1D[i1][1]);
	 r1 = outdata1D[i1][0];
	 r2 = FTdata1D[i1][0];
	 
	 r1 = r1 - r2;
	 
	 s1 = outdata1D[i1][1];
	 s2 = FTdata1D[i1][1];
	 
	 s1 = s1 - s2;
	 //printf("\nr1=%f s1=%f ",r1, s1);

	 r1 = (fabs(r1)>fabs(s1))?fabs(r1):fabs(s1);
	 max = (max > r1)?max:r1;
	 
	 
   }
 printf("\nmax=%f ", max);
 getchar();
  */		   
  }




//scale = 1.0;
void  Reconstruction::VolBsplineBaseFT2(int i, int j, int k, float X, float Y, float Z, int iv, int jj, fftw_complex *B_k, int index)
{

 float s, sX, sY, sZ,X_2, Y_2, Z_2, t;
  
 
 if(index == 1)
   s =  ReBFT[2*iv*img2dsize + 2*jj+0];

 else if (index == 2) 
   s = ReBFT[2*iv*img2dsize + 2*jj+1];

  t = i*Bscale*X + j*Bscale*Y + k*Bscale*Z;
  B_k[0][0] = s * cos(t);

  B_k[0][1] = -s * sin(t);
 



}


/************************************************************************/
void  Reconstruction::VolBsplineBaseFT(int i, int j, int k, float X, float Y, float Z, float scale, fftw_complex *B_k)
{

  float s, sX, sY, sZ,X_2, Y_2, Z_2, x, y, z, t;
  //  fftw_complex *B_k;
  // B_k = (fftw_complex *)fftw_malloc(sizeof(fftw_complex));

  x = X;
  y = Y;
  z = Z;


  X   = scale * X;
  Y   = scale * Y;
  Z   = scale * Z;



  X_2 = X/2.0;

  Y_2 = Y/2.0;

  Z_2 = Z/2.0;

  if( fabs(X_2) == 0.0 ) sX = 1.0;
  else sX = sin(X_2)/X_2;

  if( fabs(Y_2) == 0.0 ) sY = 1.0;
  else sY = sin(Y_2)/Y_2;


  if( fabs(Z_2) == 0.0 ) sZ = 1.0;
  else sZ = sin(Z_2)/Z_2;

  s = sX * sY * sZ;

  t = s * s;
  s = t * t;

  t = i * X + j * Y + k * Z;
  B_k[0][0] =  fabs(scale*scale*scale)* s * cos(t);

  B_k[0][1] = -fabs(scale*scale*scale) * s * sin(t);
 
  // printf("\n ok= %f %f X_2=%f Y_2= %f Z_2=%f s=%f", B_k[0][0], B_k[0][1], X_2, Y_2, Z_2, s); 
  //return B_k;


}

void Reconstruction::BsplineBaseFT2D(int i, int j, float X, float Y, float scale, fftw_complex* B_k)
{
  float x, y, X_2, Y_2, sX, sY, s, t;
  x = X;
  y = Y;
 
  X = scale * X;
  Y = scale * Y;

  X_2 = X/2.0;
  Y_2 = Y/2.0;

  if( fabs(X_2) == 0.0 ) sX = 1.0;
  else sX = sin(X_2)/X_2;

  if( fabs(Y_2) == 0.0 ) sY = 1.0;
  else sY = sin(Y_2)/Y_2;

  s = sX * sY;

  t = s * s;
  s = t * t;

  t = i * X + j * Y;

  B_k[0][0] =  fabs(scale*scale) * s * cos(t);
 
  B_k[0][1] = -fabs(scale*scale) * s * sin(t);
 
}

/************************************************************************/
/* float Reconstruction::testMatrix()
 {
   int i,j;
   fftw_complex *s1, *s2, *s;

   for ( i = 0; i <= N; i++ )
	 {
	   printf("\n");

	 for ( j = 0; j <= N; j++ )
	   {

		 s1 = BsplineFT(i,0.5, Bscale);
		 s2 =  BsplineFT(j, 0.5,Bscale);
		 s[0][0] = s1[0][0] * s2[0][0] - s1[0][1] * s2[0][1];
		 s[0][1] = s1[0][0] * s2[0][1] + s1[0][1] * s2[0][0];


		 printf("%f+i%f  ", s[0][0], s[0][1]);
	   }
	 }



   
 }
*/

/************************************************************************/

void Reconstruction::BsplineFT1D(int k, float omega, float scale,fftw_complex *B_k)
{
  float s,omega_2, Omega;
  
 
  Omega = omega;

  omega = omega * scale;
  omega_2 = omega/2.0;
  
  if(fabs(omega_2) == 0.0) s = 1.0;

  else
	s = (sin(omega_2)/(omega_2));

  s = s * s * s * s;
  s = s * s;

  B_k[0][0] = fabs(scale) * fabs(scale) * s * cos(k * omega);

  B_k[0][1] = -fabs(scale) * fabs(scale) * s * sin(k * omega);

}







/************************************************************************
Description:
        General image structure initialization with data memory allocation.
@Algorithm:
        Whenever a new image is read or created, this function should be
        called first. It allocates memory and sets up a number of defaults.
        The Projimg structures and data block are also allocated.

Arguments:
//        DataType datatype       the new data type.
        int  c         number of channels.
        int  x         x-dimension.
        int  y         y-dimension.
        int  z         z-dimension.
        int  n         number of images.               
@Returns:
        Oimage*                         the new image structure, NULL if initialization failed.
**************************************************************************/
Oimage* Reconstruction::InitImageParameters(int c, int  x,
                                int y, int  z, int n)
{
  int  i;

  if ( n < 1 ) return(NULL);

  // Call the more primitive initialization
  
  Oimage*         p = (Oimage*)malloc(sizeof(Oimage));
  //InitImageHeader(c, x, y, z, n);
  
  if ( !p ) return(p);
  
  int  datasize = c*x*y*z*n;
  
  printf("\ndatasize=%u \n", datasize);
  
  // Allocate data memory
  
  p->data = (float *)malloc(datasize*sizeof(float));
  
  for ( i = 0; i < datasize; i++ )
		  p->data[i] = 0.0;
	

        return(p);
}



/************************************************************************
Description:
        General image header structure initialization.
Algorithm:
        Whenever a new image is read or created, this function should be
        called first. It allocates memory and sets up a number of defaults.
        The Projimg structures are also allocated, but not the data block.
Arguments:
//        DataType datatype       the new data type.
        int  c         number of channels.
        int  x         x-dimension.
        int  y         y-dimension.
        int  z         z-dimension.
        int  n         number of images.               
Returns:
        Oimage*                 the new image structure, NULL if initialization failed.
**************************************************************************/
Oimage* Reconstruction::InitImageHeader(int c, int x,
                                int y, int z, int n)
{
        if ( n < 1 ) return(NULL);

        // Call the more primitive initialization
        Oimage* p =  InitImage();

        if ( !p ) return(p);

        p->nx = x;
        p->ny = y;
        p->nz = z;
        p->c = c;
        
		// p->image = (Projimg *)malloc(n * sizeof(Projimg));

//        if ( verbose & VERB_DEBUG )
            //    printf("DEBUG init_img_header: Image created\n");

        return(p);
}


/************************************************************************
Description:
        General image structure initialization.
Algorithm:
        Whenever a new image is read or created, this function should be
        called first. It allocates memory and sets up a number of defaults.
Arguments:
        .
Returns:
        Oimage*                         new image structure, NULL if initialization failed.
**************************************************************************/
Oimage*  Reconstruction::InitImage()
{
        // Allocate memory for the image parameter structure
        Oimage*  p = (Oimage*)malloc(sizeof(Oimage));
        if ( !p ) return(p);

        // Set parameter defaults
        p->nx = p->ny = p->nz = ImgNx;
        p->c = 1;
        p->ux = p->uy = p->uz = 0.0;
   
        return(p);
}

/************************************************************************/
int Reconstruction::CountList(char* list)
{
        if ( !list ) return(0);

        int              n = 0;
        char**          item;

        for ( item = (char **) list; item; item = (char **) *item ) n++;

        return(n);
}



VolMagick::Volume * Reconstruction::GetVolume(Oimage *p) {

  VolMagick::Volume * result = new VolMagick::Volume;

 result->voxelType(VolMagick::Float);
 printf("\np->nx=%d %d %d ", p->nx,p->ny, p->nz);

 // for volume data.
 result->dimension(VolMagick::Dimension(p->nx,p->ny,p->nz));
 printf("\np->nx ny nz = %d %d %d", p->nx, p->ny, p->nz); 

 float minExt[3];//={0.0,0.0,0.0};
 float maxExt[3];
 
 minExt[0] = StartXYZ[0];//bgrids->StartXYZ[0];
 minExt[1] = StartXYZ[1];//bgrids->StartXYZ[1];
 minExt[2] = StartXYZ[2];//bgrids->StartXYZ[2];


 maxExt[0] = FinishXYZ[0];//bgrids->FinishXYZ[0]-1;
 maxExt[1] = FinishXYZ[1];//bgrids->FinishXYZ[1]-1;
 maxExt[2] = FinishXYZ[2];//bgrids->FinishXYZ[2]-1;
 
 printf("\nbounder %f %f %f %f %f %f", StartXYZ[0], StartXYZ[1], StartXYZ[2], FinishXYZ[0], FinishXYZ[1], FinishXYZ[2]);
 
 result->boundingBox(VolMagick::BoundingBox(minExt[0],minExt[1],minExt[2], maxExt[0],maxExt[1],maxExt[2]));
 

 for(VolMagick::uint64 i = 0; i < p->nx; i++)
   for(VolMagick::uint64 j = 0; j < p->ny; j++)
	 for(VolMagick::uint64 k = 0; k < p->nz; k++)
	   {
		 //result(i,j,k, p->data[k*p->ny*p->nx+j*p->nx+i]);
	     (*result)(i,j,k, p->data[i*p->ny*p->nx+j*p->nx+k]);
			  // printf("result = %f    ",  p->data[k*p->x*p->y+j*p->x+i]);
	   }
 
 return result;

}

boost::tuple<bool,VolMagick::Volume> Reconstruction::ConvertToVolume(Oimage *p)
{

 VolMagick::Volume result;

 result.voxelType(VolMagick::Float);
 printf("\np->nx=%d %d %d ", p->nx,p->ny, p->nz);

 // for volume data.
 result.dimension(VolMagick::Dimension(p->nx,p->ny,p->nz));
 printf("\np->nx ny nz = %d %d %d", p->nx, p->ny, p->nz); 

 float minExt[3];//={0.0,0.0,0.0};
 float maxExt[3];
 
 minExt[0] = StartXYZ[0];//bgrids->StartXYZ[0];
 minExt[1] = StartXYZ[1];//bgrids->StartXYZ[1];
 minExt[2] = StartXYZ[2];//bgrids->StartXYZ[2];


 maxExt[0] = FinishXYZ[0];//bgrids->FinishXYZ[0]-1;
 maxExt[1] = FinishXYZ[1];//bgrids->FinishXYZ[1]-1;
 maxExt[2] = FinishXYZ[2];//bgrids->FinishXYZ[2]-1;
 
 printf("\nbounder %f %f %f %f %f %f", StartXYZ[0], StartXYZ[1], StartXYZ[2], FinishXYZ[0], FinishXYZ[1], FinishXYZ[2]);
 
 result.boundingBox(VolMagick::BoundingBox(minExt[0],minExt[1],minExt[2], maxExt[0],maxExt[1],maxExt[2]));
 

 for(VolMagick::uint64 i = 0; i < p->nx; i++)
   for(VolMagick::uint64 j = 0; j < p->ny; j++)
	 for(VolMagick::uint64 k = 0; k < p->nz; k++)
	   {
		 //result(i,j,k, p->data[k*p->ny*p->nx+j*p->nx+i]);
		 result(i,j,k, p->data[i*p->ny*p->nx+j*p->nx+k]);
			  // printf("result = %f    ",  p->data[k*p->x*p->y+j*p->x+i]);
	   }
 
 return boost::make_tuple(true,result);
 
}


/*****************************************************************************/
void  Reconstruction::SaveVolume(Oimage *p)
{

 VolMagick::Volume result;



 result.voxelType(VolMagick::Float);
 //printf("\np->nx=%d %d %d ", p->nx,p->ny, p->nz);

 // for volume data.
 result.dimension(VolMagick::Dimension(p->nx,p->ny,p->nz));
 //printf("\np->nx ny nz = %d %d %d", p->nx, p->ny, p->nz); 

 float minExt[3];//={0.0,0.0,0.0};
 float maxExt[3];
 
 minExt[0] = StartXYZ[0];//bgrids->StartXYZ[0];
 minExt[1] = StartXYZ[1];//bgrids->StartXYZ[1];
 minExt[2] = StartXYZ[2];//bgrids->StartXYZ[2];


 maxExt[0] = FinishXYZ[0];//bgrids->FinishXYZ[0]-1;
 maxExt[1] = FinishXYZ[1];//bgrids->FinishXYZ[1]-1;
 maxExt[2] = FinishXYZ[2];//bgrids->FinishXYZ[2]-1;
 
 printf("\nbounder %f %f %f %f %f %f\n", StartXYZ[0], StartXYZ[1], StartXYZ[2], FinishXYZ[0], FinishXYZ[1], FinishXYZ[2]);
 printf("\nBegin new reconstruction process.");
 
 result.boundingBox(VolMagick::BoundingBox(minExt[0],minExt[1],minExt[2], maxExt[0],maxExt[1],maxExt[2]));

std::cout<<"Filling in output volume with data\n"; 

 for(VolMagick::uint64 i = 0; i < p->nx; i++)
   for(VolMagick::uint64 j = 0; j < p->ny; j++)
	 for(VolMagick::uint64 k = 0; k < p->nz; k++)
	   {
		 result(i,j,k, p->data[i*p->ny*p->nx+j*p->nx+k]);
	   }

  //volfilename.clear();
  //volfilename("testlm.rawiv"); 
  //volfilename = "testlm.rawiv";
 
std::cout<<"Creating volume file\n";

     VolMagick::createVolumeFile(cvcapp, volfilename,
                                  result.boundingBox(),
                                  result.dimension(),
                                  std::vector<VolMagick::VoxelType>(1, result.voxelType()));


 std::cout<<"Writing out output file "<<volfilename<<"\n";


  VolMagick::writeVolumeFile(cvcapp, result, volfilename);
  

}



boost::tuple<bool,VolMagick::Volume> Reconstruction::gd2Volume(float *data)
{
//for gd projection.
  
  VolMagick::Volume result;
  result.voxelType(VolMagick::Float);

 result.dimension(VolMagick::Dimension(ImgNx,ImgNy,4));

 float minExt[3];
 float maxExt[3];


 minExt[0] = bgrids->StartXYZ[0]+0.5;
 minExt[1] = bgrids->StartXYZ[1]+0.5;
 minExt[2] = -0.01;//bgrids->StartXYZ[2]+0.5;
 

 maxExt[0] = bgrids->FinishXYZ[0]-0.5;
 maxExt[1] = bgrids->FinishXYZ[1]-0.5;
 maxExt[2] = 0.01;//bgrids->FinishXYZ[2];

 result.boundingBox(VolMagick::BoundingBox(minExt[0],minExt[1],minExt[2], maxExt[0],maxExt[1],maxExt[2]));
 

 for(VolMagick::uint64 i = 0; i < ImgNx; i++)
      for(VolMagick::uint64 j = 0; j < ImgNy; j++)
         for(VolMagick::uint64 k = 0; k < 4; k++)
             {
			   if(k==0 || k==3 )continue;
			    result(i,j,k, data[j*ImgNy+i]);
			  //printf("result = %f    ",  data[j*ImgNy+i]);
             }
 //getchar();

 return boost::make_tuple(true,result);
 
 
 //free(boundingbox);


}




/*************************************************************************
Descriptions:
     Test function:  f(x,y,z) = exp{-(|x-0.5| + |y-0.5|+|z-0.5|)};
     Cylinder     :  (x-0.5)^2+(y-0.5)^2 = 0.3^2;   0.2 =< z <= 0.8;                

**************************************************************************/
// m * m * m spline volume.
Oimage *Reconstruction::Phantoms( int dimx, int dimy, int dimz, int object)
{
 
  int i,j,k,index, size;
 float x,y,z, width, length, R, oz; 
 
 Oimage*   p ;

 p    = InitImageParameters(1, dimx,dimy,dimz, nv);
 size = dimx * dimy * dimz;

//for Cylinder.
 switch(object)
   {
   case 0:
	 oz = ox;
	 //R = (nx-1)/2.0;

	 R = ImgNx/2.0;

	 //length = (nx-1)/2.0;

	 length = ImgNx/2.0;
	 length = length/2.0;
	 length = length * length;
	 //printf("\nox oy oz=%f %f %f length=%f Bscale=%f", ox, oy, oz, length, Bscale);
	 // getchar();
	 width  = sqrt(R*R - length);

	 for ( i = 0; i < dimx; i++ )
	   for ( j = 0; j < dimy; j++ )
		 for ( k = 0; k < dimz; k++ )
		   {
			 x = i*Bscale;
			 y = j*Bscale;
			 z = k*Bscale;
			 //printf("\nx y z = %f %f %f ", x, y, z);			 
			 index = i*dimx*dimx + j*dimx + k;
			 if( ((x-ox)*(x-ox) + (y-oy)*(y-oy)) <= length &&                 fabs(z - oz)<= width) 
			   p->data[index] = 1.0;
			 else 
			   p->data[index] = 0.0;
			 //printf("p->data=%f   ", p->data[index]);
		   }
	 printf("\ncylinder");
	 break;
	 

   case 1: printf("\ntest");
	 break;
   }

 //getchar();



 // convert to coefficients.
 /*   
 float *coefs, *coeffs, *Ocoeffs;
 coefs        = (float *)malloc(size * sizeof(float));
 coeffs       = (float *)malloc(size * sizeof(float));
 coefficients = (float *)malloc(size * sizeof(float));
 Ocoeffs      = (float *)malloc(size * sizeof(float));

 for ( i = 0; i < size; i++ )
   {
	 // printf("\np->data=%f ", p->data[i]);
	 coefs[i] = p->data[i];
   }


 for ( i = 0; i < size; i++ )
   {
	 coeffs[i] = 0.0;coefficients[i] = 0.0;
   }

  bspline->ConvertToInterpolationCoefficients_3D(coefs, dimx, dimy,dimz,CVC_DBL_EPSILON);




 bspline->GramSchmidtofBsplineBaseFunctions2();
bspline->Evaluate_OrthoBSpline_Basis_AtVolGrid2();


// for test.

 for ( i = 0; i < size; i++ )
   {
	 coefficients[i] = coefs[i];
   }

 

  ConvertCoeffsToNewCoeffs(coefs, coefficients, bspline->InvSchmidtMat);//convert to ortho bspline basis coefficients.
 */
// convert to coefficients over.




 // printf("\northo coeffs.\n");
 //for ( i = 0; i < (N+1)*(N+1)*(N+1); i++ )
 //printf("coeffs=%f orthNon_o_coeffs=%f                   ", coefs[i], coefficients[i]);

// getchar();

 /*
 ConvertCoeffsToNewCoeffs(coefficients, Ocoeffs, bspline->InvSchmidtMat);
  printf("\nnew origin coeffs.\n");
 for ( i = 0; i < (N+1)*(N+1)*(N+1); i++ )
  if(fabs(coefs[i]-Ocoeffs[i]) > SMALLFLOAT) printf("%f ", Ocoeffs[i]);
 printf("\nOK!");
 getchar();
 */
 
 //over.

 //p->image->ox = ox;  
 //p->image->oy = oy;  
 //p->image->oz = oy; 

p->nx         = dimx;  
p->ny         = dimy;  
p->nz         = dimz;  

p->ux        = Bscale;   
p->uy        = Bscale;   
p->uz        = Bscale;   


// V0 = PI * 0.09 * 0.6;
//p->image->vx = 1;  
//p->image->vy = 1;  
//p->image->vz = 1;  


//gd = ProjectionImageFromView(p,  nview);


//free(coefs);
// free(coeffs);
//free(Ocoeffs);


//p->image->angle = 0;
//p->image->background = 0; 

return(p);



}

//n * n * n image Volume.
Oimage *Reconstruction::Phantoms_Volume(int object)
{
int i,j,k,l;
int  size, N2, index;

float x,y,z, width, length, a, b, c;
float X, Y, Z, R, OX, OY, OZ,Ox, Oy, Oz, r; 
 
Oimage*   p ;

// p    = InitImageParameters(1, ImgNx,ImgNy,ImgNz, nv);
/* size = (ImgNx * ImgNy * ImgNz;
  printf("\nsize=%u ", size);

*/
p    = InitImageParameters(1, ImgNx+1,ImgNy+1,ImgNz+1, 1);  //6.6  change to ImgNx +1.
  
N2   = (ImgNy+1) * (ImgNz+1);
 
//printf("\nobject=%d ", object);getchar();
//for Cylinder.
switch(object) {
   case 0:
      OX = (ImgNx)/2.0;
      OY = OZ = OX;

      R = (ImgNx-2)/2.0;
      //R = ImgNx/2.0;
      //R = nx/2.0;
      //length = (nx-1)/2.0;
      length = (ImgNx-1)/2.0;
      //length = (ImgNx)/2.0;

      length = length/2.0;
      length = length * length;
	 //printf("\nox oy oz=%f %f %f length=%f Bscale=%f", ox, oy, oz, length, Bscale);
	 //getchar();
       width  = sqrt(R*R - length); // image lies in a sphere x^2 + y^2 + z^2 <= R^2.

       for ( i = 0; i <= ImgNx; i++ ) {
	  for ( j = 0; j <= ImgNy; j++ ) {
	     for ( k = 0; k <= ImgNz; k++ ) {
	        x = i-OX;
	        y = j-OY;  
	        z = k-OZ;  

	        //printf("\ImgNx y z = %f %f %f ", x, y, z);			 
	        index = i*N2 + j*(ImgNz+1) + k;
	        if( (x*x + y*y) <= length && fabs(z)<= width) 
	           p->data[index] = 1.0;
	        else 
	           p->data[index] = 0.0;
	            //if ( fabs(z) <= 0.00001)  printf("x=%f y=%f z=%f  =p->data=%f   ", x,y,z,p->data[index]);
            }
         }
      }

      V0 = M_PI * length * 2*width;
      volfilename = "Cylinder.rawiv";
      //printf("\nCylinder Volume = %f ", V0);
      break;
	 

   case 1:   // Large Ellisoid Sphere and Sphere." 
      a  = (ImgNx-1)/3.0;
      b  = (ImgNx-1)/4.0;
      c  = b;

      Ox = (ImgNx)/2.0;//0.0;//(ImgNx-1)/2.0 - 2.0;
      Oy = (ImgNx)/2.0;//0.0;//(ImgNx-1)/4.0;
      Oz = (ImgNx)/4.0;

      //r  = (ImgNx-1)/3.0;
       //r  = r * r


      OX = (ImgNx)/2.0;//0.0;//(ImgNx-1)/2.0;
      OY = (ImgNx)/2.0;//0.0;//3*(ImgNx-1)/4.0;//(ImgNy-1)/2.0 + 1.0;
      OZ = 3*(ImgNx)/4.0;      //OY;
      R  = (ImgNx-1)/4.0-1.0;//(ImgNx-1)/2.0 - 2.0;
      R  = R * R;
      printf("\nOx Oy Oz r=%f %f %f %f OX OY OZ R=%f %f %f %f", Ox, Oy, Oz, r, OX, OY, OZ, R);

      for ( i = 0; i <= ImgNx; i++ ) {
         for ( j = 0; j <= ImgNy; j++ ) {
	    for ( k = 0; k <= ImgNz; k++ ) {
	       x = (i - Ox)/a;
	       y = (j - Oy)/b;
	       z = (k - Oz)/c;

	       X = i - OX;
	       Y = j - OY;
	       Z = k - OZ;
			 
	       index = i*N2 + j*(ImgNz+1) + k;
	       if ( (x*x + y*y+ z*z)<= 1 || (X*X+Y*Y+Z*Z)<= R) 
		  p->data[index]= 1.0;
	       else 
	          p->data[index] = 0.0;
	          // printf("Object->data[ii] = %f ", p->data[index]);		 
	    }
         }
      }
	 
      printf("\nEllisoid Sphere and Sphere.");
      volfilename = "Ellisoid_sphere.rawiv";
      break;

   case 2:   //Sphere.

      OX = ImgNx/2.0;//0.0;//(ImgNx-1)/2.0;
      OY = ImgNx/2.0;//0.0;//3*(ImgNx-1)/4.0;//(ImgNy-1)/2.0 + 1.0;
      OZ = ImgNx/2.0;//0.0;//3*(ImgNx-1)/4.0;//(ImgNy-1)/2.0 + 1.0;
      //OZ = 3*(ImgNx-1)/4.0;      //OY;

      printf("OXYZ ============ %f, %f, %f\n", OX, OY, OZ);
      printf("OXYZ ============ %f, %f, %f\n", OX, OY, OZ);

      R  = ImgNx/4.0-1.0;//(ImgNx-1)/2.0 - 2.0;

      R  = R * R;

      for ( i = 0; i <= ImgNx; i++ ) {   //6.6
         for ( j = 0; j <= ImgNy; j++ ) {
	    for ( k = 0; k <= ImgNz; k++ ) {
	       X = i - OX;
	       Y = j - OY;
	       Z = k - OZ;
			 
	       index = i*N2 + j*(ImgNz+1) + k; //6.6
	       if ((X*X+Y*Y+Z*Z)<= R) 
		  p->data[index]=  1.0;
	       else 
	          p->data[index] = 0.0;
		  //if ( Z == 0.0)  printf("\nX=%f Y=%f Z=%f  =p->data=%f   ", X,Y,Z,p->data[index]);

                if (fabs(X) < 4 && fabs(Y) < 4 && fabs(Z) < 4) 
                   p->data[index]=  1.0;
	        else 
	           p->data[index] = 0.0;           

		if (fabs(X) +  fabs(Y) + fabs(Z) < 4) 
	           p->data[index]=  1.0;
	        else 
	           p->data[index] = 0.0;    

	        // printf("Object->data[ii] = %f ", p->data[index]);		 
		// p->data[index]= exp(-(X*X + Y*Y + Z*Z)/R);
                //p->data[index] = 0.0;
                //p->data[index] = 1.0;
            }
         }
      }
      volfilename = "sphere.rawiv";
      break;
	 
   case 3:   // Xu added, a set of  Spheres." 
      float cx[10], cy[10], cz[10], w,alf, radius;
      int l, lx;

      alf = 1.0;
      radius = 4.0;

      for ( i = 0; i <= ImgNx; i++ ) {
         for ( j = 0; j <= ImgNy; j++ ) {
           for ( k = 0; k <= ImgNz; k++ ) {
              index = i*N2 + j*(ImgNz+1) + k;
              p->data[index] = 0.0;
           }
         } 
      }

      OX = ImgNx/2.0;//0.0;//(ImgNx-1)/2.0;
      OY = ImgNx/2.0;//0.0;//3*(ImgNx-1)/4.0;//(ImgNy-1)/2.0 + 1.0;
      OZ = ImgNx/2.0;//0.0;//3*(ImgNx-1)/4.0;//(ImgNy-1)/2.0 + 1.0;

      for (l = 0; l < 10; l++) {
         lx = 2*l/6;
         lx = 2*l - 6*lx; 
         cx[l] = lx + OX-3.0; 
         
         lx = 4*l/10;
         lx = 4*l - 10*lx; 
         cy[l] = lx + OY -5.0; 

         lx = 6*l/13;
         lx = 6*l - 13*lx; 
         cz[l] = lx  + OZ - 6.0; 

         for ( i = 0; i <= ImgNx; i++ ) {
            for ( j = 0; j <= ImgNy; j++ ) {
   	       for ( k = 0; k <= ImgNz; k++ ) {

	          index = i*N2 + j*(ImgNz+1) + k;
                  w = alf*((i-cx[l])*(i-cx[l])+(j-cy[l])*(j-cy[l])+(k-cz[l])*(k-cz[l]) - radius);
	          p->data[index] = p->data[index] + exp(-w);
	          // printf("Object->data[ii] = %f ", p->data[index]);		 
	       }
            }
         }
         for ( i = 0; i <= ImgNx; i++ ) {
            for ( j = 0; j <= ImgNy; j++ ) {
   	       for ( k = 0; k <= ImgNz; k++ ) {
	          index = i*N2 + j*(ImgNz+1) + k;
	          if (p->data[index] >= 2.0) p->data[index] = 2.0;
	       }
            }
         }
      }
      volfilename = "Xu_sphere.rawiv";	 
      printf("\n Ganssian Map.");
      break;
}

p->nx         = ImgNx+1;  
p->ny         = ImgNy+1;  
p->nz         = ImgNz+1;  

p->ux        = 1.0;   
p->uy        = 1.0;   
p->uz        = 1.0;   

return (p);
}

//-----------------------------------------------------------------------------
void Reconstruction::Phantoms_gd(Oimage *Volume, EulerAngles* eulers)
{
  int v, i,j, size, sample_num;
  float rotmat[9], *data = NULL, *prjimg=NULL; 
  char filename[100];
  sample_num = ImgNx+1;
  size = sample_num * sample_num;
  EulerAngles* home = NULL;
 
  //home = eulers;

  // data = (float *)malloc(nv * size *sizeof(float));
  prjimg = (float *)malloc(size *sizeof(float));

  
  for ( v = 0; v < nv; v++ )
	{
	  for ( i = 0; i < 9; i++ )
		{
                 rotmat[i] = Rmat[v*9+i];
                 //printf("\nrotmat=%f ", rotmat[i]);
		}	  
	  Volume_Projection(Volume, rotmat,  sample_num, prjimg);
	  //Volume_GridProjection(Volume, rotmat,  sample_num, prjimg);
	  //printf("\nBspline projection.");
	  //for ( i = 0; i < size; i++ ) 
	  for ( i = 0; i <= ImgNx; i++ )
             {
			   //printf("\n");
              for ( j = 0; j <= ImgNx; j++ ) 
		{
		  gd->data[v*size + i*(ImgNx+1) + j] = prjimg[i*(ImgNx+1) + j];
		  //if (v == 0)  printf("\ni=%d  j=%d  gd=%f ", i, j, prjimg[i*(ImgNx+1) + j]);
		  //printf("%f ", prjimg[i*(ImgNx+1) + j]);
		}
             }
          
          //if(v+1 < 10 )  sprintf(filename, "column1600000%d.spi",v+1);
          //else sprintf(filename, "column160000%d.spi",v+1);

          //WriteSpiFile(filename, gd->data+v*size, ImgNx+1, ImgNx+1, 1, home->rot, home->tilt, home->psi);
          //home = home->next;
	  
	}


 free(prjimg); prjimg = NULL;

}


Oimage* Reconstruction::InitialFunction(int function, const char* filename, const char *path)
{
  
 int i,j,k, index, index2, N2, usedN2, scale, ddim[3];
 Oimage *p ;
 
 scale  = (int)Bscale;
 N2     = (ImgNy+1)*(ImgNz+1);
 usedN2 = (usedN+1)*(usedN+1);
 p      = InitImageParameters(1, ImgNx+1,ImgNy+1,ImgNz+1, 1);
 O_Kcoef  = (float *)malloc(usedtolbsp * sizeof(float));
 

 if(filename == NULL || path == NULL ) DesignInitFunction(function, p);
 else if ( filename != NULL && path != NULL ) 
 {
//  VolMagick::VolumeFileInfo volinfo;
//  volinfo.read(filename);
//  cout << "Dimension: " << volinfo.XDim() << "x" << volinfo.YDim() << "x" << volinfo.ZDim() << endl;

  VolMagick::Volume Vol;
  VolMagick::readVolumeFile(cvcapp, Vol,filename,0,0);
  printf("\nVolMagic dim========%d %d %d ",(int)Vol.XDim(), (int)Vol.YDim(),(int)Vol.ZDim());

  ddim[0] = Vol.XDim();
  ddim[1] = Vol.YDim();
  ddim[2] = Vol.ZDim();

   for(VolMagick::uint64 i = 0; i < ddim[0]; i++)
    for(VolMagick::uint64 j = 0; j < ddim[1]; j++)
      for(VolMagick::uint64 k = 0; k < ddim[2]; k++) {
         p->data[((i*ddim[1]+j)*ddim[2]+k)] = Vol(i,j,k);
         }

 } 
   
 printf("\ninit name=====%s reslut.path===%s \n", filename, path) ;
 
   p->nx         = ImgNx+1;  
   p->ny         = ImgNy+1;  
   p->nz         = ImgNz+1;  

   p->ux        = Bscale;   
   p->uy        = Bscale;   
   p->uz        = Bscale;   
  
  

   /*
   for ( i = 2; i <= N - 2; i++ )
	for ( j = 2; j <= N - 2; j++ )
	  for ( k = 2; k <= N - 2; k++ )
		{
		  index = i*N2 + j*(ImgNz+1) + k;

	  // printf("\np->data=%f ", p->data[i]);
		  O_Kcoef[(i-2)*usedN2+(j-2)*(usedN+1)+k-2] = p->data[index];
		}
 
  */

   for ( i = 0; i <=usedN; i++ ) {
      for ( j = 0; j <= usedN; j++ ) {
         for ( k = 0; k <= usedN; k++ ) {
            index  = i*usedN2 + j*(usedN+1) + k;
            index2 = (2+i)*scale*N2 + (j+2)*scale*(ImgNz+1) + (k+2)*scale;
            O_Kcoef[index] = p->data[index2];
         } 
      }
   }


  
   bspline->ConvertToInterpolationCoefficients_3D(O_Kcoef, usedN+1, usedN+1,usedN+1,CVC_DBL_EPSILON);

  
  /*
for ( i = 0; i <=usedN; i++ )
   for ( j = 0; j <= usedN; j++ )
      for ( k = 0; k <= usedN; k++ )
		{
		  if(i = 0 && j == 0 && j == 1) O_Kcoef[index] = 1.0;
		  else O_Kcoef[index] = 0.0;
		}
  */


  if(Bscale == 1.0 ) bspline->ObtainObjectFromNonOrthoCoeffs(O_Kcoef,p->data);
  if(Bscale == 2.0 ) bspline->ObtainObjectFromNonOrthoCoeffs_FA(O_Kcoef,p->data);

  /*
	 for ( i = 0; i <= ImgNx; i++ )
	   for ( j = 0; j <= ImgNx; j++ )
		 for ( k = 0; k <= ImgNx; k++ )
		    printf("\ni =%d j =%d k = %d p=%f   ", i, j, k, p->data[i*img2dsize+ j*(ImgNx+1) + k]);getchar();
  */


 
  //for ( i = 0; i < usedtolbsp; i++ ) printf("\nO_Kcoef=%f ", O_Kcoef[i]);getchar();
  return p;

}
 
void Reconstruction::DesignInitFunction(int function, Oimage* p)
{
 int i,j,k, index, N2, scale;

 float x,y,z, oz, R;

 float OX, OY, OZ, length, width;

 scale = (int)Bscale;
 N2 = (ImgNy+1)*(ImgNz+1);


 oz = Ox;
 printf("\nOx Oy Oz=%f  %f %f in InitialFunction and function = %d", Ox, Oy,oz, function);


//memset(O_Kcoef, 0, usedtolbsp*sizeof(float));
switch(function) {
   case 0:                       //for sphere.
   OX = ImgNx/2.0;
   OY = ImgNx/2.0;
   OZ = ImgNx/2.0;
	  
   R = (ImgNx-4.0)/2.0;
   //printf("\nR=============%f ", R);
   //R = (ImgNx-3.0)/5.0;        // Xu 
 
   //R  = ImgNx/4.0-1.0;           // Xu
   R = R * R;
   printf("\nold R^2 = %f \n", R);  
   for ( i = 0; i <= ImgNx; i++ )     //6.6 
      for ( j = 0; j <= ImgNy; j++ ) //6.6
         for ( k = 0; k <= ImgNz; k++ ) {
	    x = i - OX;
	    y = j - OY;
	    z = k - OZ;
	
            index = i*N2 + j*(ImgNz+1) + k;
	    if( (x*x + y*y+ z*z)<= R ) { 
	       //p->data[index] = R - (x*x + y*y+ z*z) ;//+ rand()%10;  // continuous f
	       p->data[index] = 1.0;      // Xu         dis-continuous f
	    } else {
	       p->data[index] = 0.0;
	    }
			  
	    //                          p->data[index] = 1.0;   // Changed
	    //		p->data[index] = R - (x*x + y*y+ z*z) ;// Xu Chnage;
	    //	p->data[index] = exp(p->data[index]) - 1.0;;
	 
	    p->data[index] = 0.0;  // Xu 0915 //0924.  
   }
           
	  
   //V0 = 4.0/3 * PI * R * sqrt(R);
   printf("\nSphere V0 = %f ", V0);
  
   break;
   
   case 1:
      printf("\ninitial function is zero.");
      break;

   case 2:
      OX = (ImgNx)/2.0;
      OY = OZ = OX;

      R = (ImgNx-2)/2.0;
        //R = ImgNx/2.0;
      //R = nx/2.0;

      //length = (nx-1)/2.0;
      length = (ImgNx-1)/2.0;
      //length = (ImgNx)/2.0;

      length = length/2.0;
      length = length * length;
      //printf("\nox oy oz=%f %f %f length=%f Bscale=%f", ox, oy, oz, length, Bscale);
      //getchar();
      width  = sqrt(R*R - length); // image lies in a sphere x^2 + y^2 + z^2 <= R^2.

      for ( i = 0; i <= ImgNx; i++ ) {
         for ( j = 0; j <= ImgNy; j++ ) {
	    for ( k = 0; k <= ImgNz; k++ ) {
               x = i-OX;
               y = j-OY;  
               z = k-OZ;  

	       //printf("\ImgNx y z = %f %f %f ", x, y, z);			 
	       index = i*N2 + j*(ImgNz+1) + k;
	       if( (x*x + y*y) <= length && fabs(z)<= width) {
                  p->data[index] = 1.0;
	       } else  {
	          p->data[index] = 0.0;
	         //if ( fabs(z) <= 0.00001)  printf("x=%f y=%f z=%f  =p->data=%f   ", x,y,z,p->data[index]);
	       }
            }
         }
      }
      //V0 = PI * length * 2*width;
      printf("\nCylinder Volume = %f ", V0);
      break;
 
   case 3:   // Xu added, a set of  Spheres." 
      float cx[10], cy[10], cz[10], w,alf, radius;
      int l, lx;

      alf = 1.0;
      radius = 4.0;

      for ( i = 0; i <= ImgNx; i++ ) {
         for ( j = 0; j <= ImgNy; j++ ) {
           for ( k = 0; k <= ImgNz; k++ ) {
              index = i*N2 + j*(ImgNz+1) + k;
              p->data[index] = 0.0;
           }
         } 
      }

      OX = ImgNx/2.0;//0.0;//(ImgNx-1)/2.0;
      OY = ImgNx/2.0;//0.0;//3*(ImgNx-1)/4.0;//(ImgNy-1)/2.0 + 1.0;
      OZ = ImgNx/2.0;//0.0;//3*(ImgNx-1)/4.0;//(ImgNy-1)/2.0 + 1.0;

      for (l = 0; l < 10; l++) {
         lx = 2*l/6;
         lx = 2*l - 6*lx; 
         cx[l] = lx + OX-3.0; 
         
         lx = 4*l/10;
         lx = 4*l - 10*lx; 
         cy[l] = lx + OY -5.0; 

         lx = 6*l/13;
         lx = 6*l - 13*lx; 
         cz[l] = lx  + OZ - 6.0; 

         for ( i = 0; i <= ImgNx; i++ ) {
            for ( j = 0; j <= ImgNy; j++ ) {
   	       for ( k = 0; k <= ImgNz; k++ ) {

	          index = i*N2 + j*(ImgNz+1) + k;
                  w = alf*((i-cx[l])*(i-cx[l])+(j-cy[l])*(j-cy[l])+(k-cz[l])*(k-cz[l]) - radius);
	          p->data[index] = p->data[index] + exp(-w);
	          // printf("Object->data[ii] = %f ", p->data[index]);		 
	       }
            }
         }
         for ( i = 0; i <= ImgNx; i++ ) {
            for ( j = 0; j <= ImgNy; j++ ) {
   	       for ( k = 0; k <= ImgNz; k++ ) {
	          index = i*N2 + j*(ImgNz+1) + k;
	          if (p->data[index] >= 2.0) p->data[index] = 2.0;
	       }
            }
         }
      }
	 
      printf("\n Ganssian Map.");
      break;
   }


}

void Reconstruction::ConvertToOrthoCoeffs()
{ 

  float *Ocoeffs;
  int   memsize, i;

  Ocoeffs      = (float *)malloc(usedtolbsp * sizeof(float));
  memsize      = usedtolbsp * sizeof(float);
 
  for ( i = 0; i < usedtolbsp; i++ )
	{
	  Ocoeffs[i] = 0.0; 
	  //printf("\nusedtolbsp=%d Ocoeffs=%f ", usedtolbsp, Ocoeffs[i]);
	}
  //memset(Ocoeffs, 0.0, memsize);


  ConvertCoeffsToNewCoeffs(O_Kcoef, Ocoeffs, bspline->InvSchmidtMat);//convert to ortho bspline basis coefficients.

  /*

 float *temp =  (float *)malloc(usedtolbsp * sizeof(float));
 //memset(temp, 0.0, memsize);
 ConvertCoeffsToNewCoeffs(Ocoeffs, temp, bspline->SchmidtMat);


 for ( i = 0; i < usedtolbsp; i++ )
   if(fabs(O_Kcoef[i] - temp[i] ) != 0.0 ) printf("\n%e ", (O_Kcoef[i]-temp[i])/O_Kcoef[i]);
 getchar();
  */

  for ( i = 0; i < usedtolbsp; i++ )
	O_Kcoef[i] = Ocoeffs[i];




  //free(temp);     temp    = NULL;
  free(Ocoeffs);  Ocoeffs = NULL;
  
    
  //int i1;

//float *data = (float *)malloc(VolImgSize*sizeof(float));

//ConvertCoeffsToNewCoeffs(coefficients, coefs, bspline->SchmidtMat);

//bspline->ObtainObjectFromNonOrthoCoeffs(coefs,p->data);

// VolImg = p;

/*
for ( i = 0; i <= ImgNx; i++ )
   for ( j = 0; j <= ImgNx; j++ )
     for ( k = 0; k <= ImgNx; k++ )
           {
                 i1 = i*img2dsize + j * (ImgNx+1) + k;
                 //if(fabs(VolImg->data[i1] - VolF->data[i1]) != 0.0) printf("\np=%e VolImg=%e  i=%d j=%d k=%d", VolF->data[i1], VolImg->data[i1], i,j,k);            
                            }


 free(data);
  
*/
  /*

  printf("\nnew coef");
  for ( i = 0; i < usedtolbsp; i++ )
	{
	  printf("\ncoeff[ii] = %f ",Ocoeffs[i]);
	}  
   getchar();
  */
 


/*
 float a, b,c,d;
 float *s, f;
 long double values[6];

 int numb = 13;

 s = (float *)malloc(numb*sizeof(float));
 a = 1.0; b = 2.0; c= -1.5; d = 4.0;

 //s[0] = 0.5*(a+b);
 //s[1] = 0.5*(b+c);
 //s[2] = 0.5*(c+d);


 s[0] = s[1] = s[2] = s[10] = s[11] = s[12] = 0.0;
 s[3] = s[4] = s[5] = s[6] = s[7] = s[8] = s[9] = 1.0; 

 bspline->ConvertToInterpolationCoefficients_1D(s, numb,CVC_DBL_EPSILON ) ;

 for ( i = 0; i < numb; i++ ) printf("\ns = %f ", s[i] ); getchar();


*/
 /*
 for ( j = 0; j < 4; j++ )
   {
	 y = -0.5 + j + 0.0;
	 f = 0.0;
  for ( i = 0; i < 3; i++ )
	{
	  x = y - i;
	  bspline->Spline_N_Base(x,2,values);
	  f = f+ s[i] * values[0];
	}
  printf("\nf=%f ",f);
  getchar();
   }
 */

  //VolImg = p;

}




/************************************************************************/
float Reconstruction::SimpsonIntegrationOn2DImage(Oimage* p)
{
 if(p->nz > 1) {printf("Error.Not a plane image."); return -1;}

 int  i,j,k;
 float area, midvalue;
 float result = 0.0; 

 area = p->ux*p->uy;
 // printf("\n area=%f ", area);getchar();

 for ( i = 0; i < p->nx-1; i++ )
   {
     k = i * p->ny;

     for( j = 0; j < p->ny-1; j++ )
        {
          
         midvalue = 0.25*(p->data[k+j]+p->data[k+p->ny+j]
                             +p->data[k+j+1]+p->data[k+p->ny+j+1]);
                  
		 // printf("\ndata=%f %f %f %f %d %d ", p->data[k+j],p->data[k+p->y+j],p->data[k+j+1],p->data[k+p->y+j+1], p->x, p->y);
        result = result + area*midvalue;


       } 
   }
 //printf("\narea=%f result=%f ", area, result);
 //getchar();
return result;
 
}

//2D gaussian .
/*
float Reconstruction::GaussianIntegrationOn2DImage(Oimage* p)
{
  int i, j, ii, jj, gM, numb_in_i, numb_in_j, id;
  float sum1, sum, delt;

  gM  = 4 * FN - 1;
  sum = 0.0;
  delt= p->ux;

  for ( i = 0; i <= gM; i++ )
	{
	  id        = i * (gM + 1);
	  numb_in_i = i/4;
	  ii      = i - 4 * numb_in_i;
	  sum1    = 0.0; 
	  for ( j = 0; j <= gM; j++ )
		{
		  numb_in_j = j/4;
		  jj        = j - 4 * numb_in_j;
		  sum1     += Weit[j] * delt * p->data[id + j];
		}

	  sum += Weit[i] * delt * sum1;  
}
*/
/************************************************************************/
Oimage* Reconstruction::ImageForTestSimpsonIntegration()
{

Oimage*   p = InitImageParameters(1, 4,4,2, 1);
for(int i = 0; i < 4; i++ )
   for( int j = 0; j < 4; j++ )
      p->data[i*4+j] = i*4+j;
p->ux = p->uy = p->uz =0.5;
}


/*************************************************************************
 Descriptions:
     Initialize the variables for reconstruction.
 ************************************************************************/
bool Reconstruction::Initialize(int dimx,int dimy,int dimz, int m, int NewNv, int BandWidth,
      float Alpha, float Fac, float Volume, int Flow , int Recon_method)
{

  int i,j, ImgDim[3];
  //int StartXYZ[3],  FinishXYZ[3];
  long double h;
 
 bgrids      = (BGrids *)malloc(sizeof(BGrids)); 
 bspline     = new Bspline();
 

 ImgNx = dimx;  // inputed image dimension.  needed to be interpolated to ImgNx+1  in order to make the center at the center of Volume grids.
 ImgNy = dimy;  // inputed image dimension.
 ImgNz = dimz;  // inputed image dimension.
 //VolImgSize = ImgNx * ImgNy * ImgNz;
VolImgSize = (ImgNx+1) * (ImgNy+1) * (ImgNz+1); ///6.6

 ImgDim[0] = ImgNx;
 ImgDim[1] = ImgNy;
 ImgDim[2] = ImgNz;

 h = ImgNx*1.0L/m;
 Bscale =(float)h;
 printf("\nImgNx =%d, Bscale=%f ", ImgNx, Bscale);


//Define scale. 
//delx = 1.0;    //delx   = 1.0/(nx-1);



 printf("\nInitialize the new 3D volume.\n");
 ox         = ImgNx/2.0;
 oy         = ImgNy/2.0; 
 
 Ox = ImgNx*1.0/2.0;  //6.6
 Oy = ImgNy*1.0/2.0;  //6.6
 
/*   //6.6
 gox = ImgNx/2.0;
 goy = ImgNy/2.0;
 goz = ImgNz/2.0;
*/
 if(ImgNx != ImgNy || ImgNy!= ImgNz) 
    { 
     printf("\n Error in initial data. the three dimension along x,y,z are not equal.");
     return false;
    }

 //numb_2dimg = number;
 img2dsize  = (ImgNx+1) * (ImgNy+1);
// M          = 4*(ImgNx-1)-1;
 N          = m;     // The number B-spline intervals used. 
 usedN      = N - 4; // Why it is not N - 3? 

 tolbsp     = (N+1)*(N+1)*(N+1); // The total number of Bspline
 //msize  = 64*64*64*64*16;
 

 
 ////StartXYZ[0] = -(ImgNx/2) - 1;  
 ////StartXYZ[1] = -(ny/2) - 1;
 ////StartXYZ[2] = -(nz/2) - 1;


 ////FinishXYZ[0] = ImgNx/2 + 1;
 ////FinishXYZ[1] = ny/2 + 1;
 ////FinishXYZ[2] = nz/2 + 1;

 StartXYZ[0] = 0.0;//-(ImgNx/2);
 StartXYZ[1] = 0.0;//-(ImgNy/2);
 StartXYZ[2] = 0.0;//-(ImgNz/2);

 FinishXYZ[0] = ImgNx;
 FinishXYZ[1] = ImgNy;
 FinishXYZ[2] = ImgNz;

 ////h = (ImgNx+2.0L)/m;



 BStartX  = -(m/2-2);
 BFinishX =  m/2-2;
 printf("\nBstart=====================%d  BFinshX=%d ", BStartX, BFinishX); 
 usedtolbsp = (N - 3) * ( N - 3) * ( N-3); // (N+1 - 4)^3, the used B-spline number 
 
 bgrids->dim[0] = m+1;
 bgrids->dim[1] = m+1;
 bgrids->dim[2] = m+1;

 bgrids->size = bgrids->dim[0] * bgrids->dim[1] * bgrids->dim[2];

 bgrids->StartXYZ[0] = -(ImgNx/2);//StartXYZ[0];
 bgrids->StartXYZ[1] = -(ImgNx/2);//StartXYZ[1];
 bgrids->StartXYZ[2] = -(ImgNx/2);//StartXYZ[2];

 bgrids->FinishXYZ[0] = (ImgNx/2);//FinishXYZ[0];
 bgrids->FinishXYZ[1] = (ImgNx/2);//FinishXYZ[1];
 bgrids->FinishXYZ[2] = (ImgNx/2);//FinishXYZ[2];

 bgrids->scale = h;

 bspline->BsplineSetting(ImgDim, bgrids);


 //nv     = CountList((char *)eulers);

 pixelscale = 1.0 * 1.0;
 Vcube  = 1.0 * pixelscale;
 //V0     = Volume;

 newnv     = NewNv;
 bandwidth = BandWidth;
 alpha     = Alpha;    
 fac       = Fac; 
 flows     = Flow;

 Recon_Method = Recon_method;

 printf("\nfixed parameters----------------------------------------------------\n");
 printf("Recon_method=%d alpha=%f fac=%f\n ", Recon_Method, alpha, fac);
 printf("--------------------------------------------------------------------\n");
 Iters = 0;


 printf("\nvivew number  = %d  usedtolbsp=%d ", nv,  usedtolbsp);
 
 //O_Kcoef   = (float *)malloc(usedtolbsp * sizeof(float));
	  
	  //Rmat   = (float *)malloc(nv *9*sizeof(float));
	  //InvRmat = (float *)malloc(nv *9*sizeof(float));

	  //gd       = (PrjImg *)malloc(sizeof(PrjImg));
	  //gd->data = (float *)malloc(nv * (ImgNx+1) * (ImgNy+1) * sizeof(float));
	  ///Matrix = (float *)malloc(usedtolbsp * sizeof(float));                     //Compute variation J1. Method 0.

	  //for test. 
	  

	  //EulerMatrice(Rmat,eulers);


 // return true;

}

void Reconstruction::SetJ12345Coeffs(float ReconJ1, float Al, float Be, float Ga, float La)
{
 reconj1 = ReconJ1;
 al     = Al;
 be     = Be;
 ga     = Ga;
 la     = La;

}

void Reconstruction::setTolNvRmatGdSize(const int TolNv)
{
 nv    = TolNv;
 Rmat     = (float *)malloc(nv *9*sizeof(float));
 InvRmat  = (float *)malloc(nv *9*sizeof(float));

 gd       = (PrjImg *)malloc(sizeof(PrjImg));
 gd->data = (float *)malloc(nv * img2dsize * sizeof(float));

 EJ1_k    = (float*)malloc(nv*sizeof(float));
 EJ1_k1   = (float*)malloc(nv*sizeof(float));
}



/************************************************************************
 * @Function: kill_all_but_main_img
 * @Description:
 *         Deallocate all pointers in an image structure.
 *         @Algorithm:
 *                 All allocated blocks within a structure is freed, but not the structure.
 *                 @Arguments:
 Bimage*                         the image structure to be deallocated.
 @Returns:    
 int                                     error code (<0 means failure).
 *************************************************************************/
int Reconstruction::kill_all_but_main_img(Oimage* p)
{   
        if ( p == NULL ) return(0);
    
        //if ( p->image ) free(p->image);

        if ( p->data ) free(p->data);
        return(0); 
}   
  

/*************************************************************************
 Descriptions:
     Define a Delta Function, where alpha > 0 is a given parameter, which
     control the support of Delta Function.
 *************************************************************************/
float Reconstruction::Deltafunc(float x)
{
 
  if (fabs(x) > alpha)
	return 0;
  else 
	return (1.0/(2*alpha)*(1+cos(M_PI*x/alpha)));



} 

float Reconstruction::Deltafunc3d(float x, float y, float z)
{
  float deltaxyz;

  deltaxyz = Deltafunc(x) * Deltafunc(y) * Deltafunc(z);

  return deltaxyz; 

}



float Reconstruction::DeltafuncPartials(float x)
{
  if (fabs(x) > alpha) return 0;
  else 
	return (-0.5*M_PI/(alpha*alpha)*sin(M_PI*x/alpha));

}



/*************************************************************************
 Descriptions:
     Sort the array in the decreasing order.
**************************************************************************/
void Reconstruction::SortDecreasing(float *arr,int *oldindex,int n)
{
int maxIndex;int pass,j;float temp;

for ( pass = 0; pass < n; pass++ )
    {
     oldindex[pass] = pass;
   }


for(pass=0;pass<n-1;pass++)
   {

    maxIndex=pass;
 
    for(j=pass+1;j<n;j++)
       if(arr[j]>arr[maxIndex])
          maxIndex=j;

    if(maxIndex!=pass)
      {
       temp          = arr[pass];
       arr[pass]     = arr[maxIndex];
       arr[maxIndex] = temp;

       temp               = oldindex[pass];
       oldindex[pass]     = oldindex[maxIndex];
       oldindex[maxIndex] = temp;
      }
   }

}

void Reconstruction::SortDecreasing(float *arr,int n)
{
int maxIndex;int pass,j;float temp;

for(pass=0;pass<n-1;pass++)
   {

    maxIndex=pass;

    for(j=pass+1;j<n;j++)
       if(arr[j]>arr[maxIndex])
          maxIndex=j;

    if(maxIndex!=pass)
      {
       temp          = arr[pass];
       arr[pass]     = arr[maxIndex];
       arr[maxIndex] = temp;

      }
   }

}


/*---------------------------------------------------------------------------
 * Sortingfloat -- Sorting a float array from small to big
 * ----------------------------------------------------------------------------*/
void Sortingfloat(float *point, int *index, int number)
{
int    n,temp,gap,i,j;
n = number;
/*
 * if (n > 256) n = 256;
 * */
for (i = 0; i< n; i++) {
    index[i] = i;
}
for (gap = n/2; gap > 0; gap /= 2)
    for (i = gap; i <n; i++) {
        for (j = i - gap; j >= 0 &&
            point[index[j]] > point[index[j+gap]]; j-=gap){
            temp = index[j];
            index[j] =  index[j+gap];
            index[j+gap] = temp;
        }
    }
}




/*************************************************************************
 Descriptions:
     Estimate Iso-value c (IsoC) from volume data f.
     Return   the non-zero volume index set of Delta Function (f-IsoC).
 *************************************************************************/
void Reconstruction::EestimateIsoValueC(float *f, float *IsoC,int n,vector<int> &suppi)
{
  int i,*oldindex=NULL, cindex;
  float *sortf = (float*)malloc(n*sizeof(float));

  for ( i = 0; i < n ;i++ ) sortf[i] = f[i];

  
//  SortDecreasing(f, oldindex,n);
    SortDecreasing(sortf,n);


  cindex = V0/Vcube;
  *IsoC   = sortf[cindex-1];


//old methods.
  for ( i = 0; i <n; i++ )
     if( fabs(f[i] - *IsoC) <= alpha )
	  {
	   suppi.push_back(i);
          }

  free(sortf);   sortf = NULL;  
}

//new methods.
void Reconstruction::EestimateIsoValueC(float *f, float *IsoC, int n, vector<CUBE> &cubes, 
                             vector<CUBE> &cubes1, float param)  //, vector<float> &DeltafIsoC,  vector<float> &PdeltafIsoC)
{
  int i,cindex;
  float *sortf = (float*)malloc(n*sizeof(float));
  CUBE onecube;
  float palpha, sum;
  int   index[8];
  
  palpha = param * alpha;

  for ( i = 0; i < n ;i++ ) sortf[i] = f[i];
  SortDecreasing(sortf,n);
  cindex = V0/Vcube;
  *IsoC   = sortf[cindex-1];
  printf("\nmax = %f , mini= %f, isoC = %f, iso_in_display = %f", sortf[0], sortf[n-1], *IsoC, 255*(*IsoC-sortf[n-1])/(sortf[0]-sortf[n-1]));

  int j, k, p,q,r, ii;
  float s[8];
  for ( i = 0; i < ImgNx; i++ )
     for ( j = 0; j < ImgNx; j++ ) 
        for ( k = 0; k < ImgNx; k++ )
           {

             for ( p = 0; p <2; p++ )
             for ( q = 0; q <2; q++ )
             for ( r = 0; r <2; r++ )
                s[p*4+q*2+r] = f[(i+p)*img2dsize+(j+q)*(ImgNx+1)+k+r];

             // Xu change sorting, because sorting make function values disorderd 
             //SortDecreasing(s, 8);
             Sortingfloat(s, index, 8);


             //printf("\nds8-s1=%f ", s[0]-s[7]);
             // Xu changed 
             //if(s[0]-*IsoC < -alpha || s[7]-*IsoC > alpha ) continue;
             if(s[index[7]]-*IsoC < -alpha || s[index[0]]-*IsoC > alpha ) continue;
             onecube.ox = i+0.5;
             onecube.oy = j+0.5;
             onecube.oz = k+0.5;
             onecube.l  = 1.0;
             /* 
             // Xu comment out 
             onecube.cubef8[0] = s[0];
             onecube.cubef8[1] = s[1];
             onecube.cubef8[2] = s[2];
             onecube.cubef8[3] = s[3];
             onecube.cubef8[4] = s[4];
             onecube.cubef8[5] = s[5];
             onecube.cubef8[6] = s[6];
             onecube.cubef8[7] = s[7];
             */
              
             //if (s[0] - s[7] < palpha) {
             if (s[index[7]] - s[index[0]] < palpha) { // Xu changed

                // Xu asked: Are the following 5 lines useful?
                //sum  = 0.125*(s[0]+s[1]+s[2]+s[3]+s[4]+s[5]+s[6]+s[7]);
                //onecube.cubef8[0] = sum;                 // the first element saves the middle point value.
                //sum  = sum - *IsoC;
                //DeltafIsoC.push_back(Deltafunc(sum));
                //PdeltafIsoC.push_back(DeltafuncPartials(sum));

                onecube.l = onecube.l * onecube.l * onecube.l;
                cubes1.push_back(onecube);
// printf("\nonecube.oxyz=%f %f %f l=%f ", onecube.ox, onecube.oy, onecube.oz, onecube.l);

            } else {

                // Xu move to here
                onecube.cubef8[0] = s[0];
                onecube.cubef8[1] = s[1];
                onecube.cubef8[2] = s[2];
                onecube.cubef8[3] = s[3];
                onecube.cubef8[4] = s[4];
                onecube.cubef8[5] = s[5];
                onecube.cubef8[6] = s[6];
                onecube.cubef8[7] = s[7];

                cubes.push_back(onecube);
            }


             //cubes.push_back(onecube);

/*
             printf("\n");
             for ( p = 0; p <2; p++ )
             for ( q = 0; q <2; q++ )
             for ( r = 0; r <2; r++ )
                printf("%e ", f[(i+p)*img2dsize+(j+q)*(ImgNx+1)+k+r]-*IsoC);
*/

            
/*
//old methods.
                pointxyz.at(0).push_back(i);   //pointxyz.at(0) is x.
                pointxyz.at(1).push_back(j); 
                pointxyz.at(2).push_back(k); 
*/

           }

  free(sortf);   sortf = NULL;  

}



//test 1D schmidt .
float Reconstruction::Test1DSchmidtMatrix()
{
  /*
  float *f, *coefs, *Non_o_coefs;
  int i,j;

  f     = (float *)malloc(ImgNx*sizeof(float));
  coefs = (float *)malloc(ImgNx*sizeof(float));
  Non_o_coefs = (float *)malloc(ImgNx*sizeof(float));

  for (i = 0; i < ImgNx; i++ )
	{f[i] = 1.0;coefs[i] = f[i];}


  bspline->ConvertToInterpolationCoefficients_1D(coefs, ImgNx,CVC_DBL_EPSILON);


*/

}



void Reconstruction::ObtainObjectFromCoeffs(float *O_coef_Result,float *f)
{

  int   i, j, k, usedN, usedN2, usedN3;
  int   ii,jj, iii,jjj,kkk,xN3,id, jd, yN3,zN3,NN, N2;
 float    Phi_i, Phi_j, Phi_k, Phi_ij;
 /*
 usedN = N - 4;

 NN = ny*nz;
 usedN2 = (usedN+1)*(usedN+1);
 usedN3 = 3*(usedN+1);


 //3 = 3*(N+1);

 int x, y,z;

 for ( x = 0; x < nx; x++ ) 
   {
	 iii = x * NN;
     xN3 = x * usedN3;

    for ( y = 0; y < ny; y++ )
	  {
		jjj = y * nz;
		yN3 = y * usedN3;
 
       for ( z = 0; z < nz; z++ )
          {
           ii  = iii + jjj + z;
		   zN3 = z * usedN3;
		   //printf("\nii=%u ", ii);
          for ( i = 0; i <= usedN; i++ )
			{
			   Phi_i =  OrthoBBaseGrid2[xN3+i+i+i];
			  //Phi_i = bspline->BernBaseGrid[xN3+i+i+i]; 
			  id = i*usedN2;
			  //printf("\nx=%d,Phi_i = %f ", x, Phi_i);
             for ( j = 0; j <= usedN; j++ )
			   {
				 Phi_j  = OrthoBBaseGrid2[yN3+j+j+j];
				 // Phi_j = bspline->BernBaseGrid[yN3+j+j+j];  
                 Phi_ij = Phi_i * Phi_j;
				 // printf("y=%d, Phi_j %f ", y, Phi_j);
				 jd = j*(usedN+1);

                for ( k = 0; k <= usedN; k++ )
                   {
                    jj = id+ jd + k;

                    Phi_k = OrthoBBaseGrid2[zN3+k+k+k];
					//Phi_k = bspline->BernBaseGrid[zN3+k+k+k]; 

					//					printf("z=%d, Phi_k %f ", z, Phi_k);  
					//	printf("\nx = %d y = %d z = %d i = %d j=%d k=%d Phi_i=%f, Phi_j=%f, Phi_k=%f,Phi_ijk=%f ", x, y, z,i, j, k, Phi_i , Phi_j, Phi_k,Phi_ij*Phi_k);

                    f[ii] = f[ii] + O_coef_Result[jj] * Phi_ij * Phi_k;

                   }
				//	printf("\nii = %u f=%f ", ii, f[ii]);
				//	getchar();
 
               }
			}
		  }
	  }
   }
 */
}





void Reconstruction::simple_backprojection(float *oneprj, float *invrotmat, Oimage *vol, int diameter)
{
  int i, j, k, l, m, half;
  int  N2, index;
  
  
  float dim2, x, y, z, xp, yp;
  float value1, value2, scalex, scaley, scale1, value;
  float radius2, x2, y2, z2, z2_plus_y2;
  
  // Use minus-tilt, because code copied from OldXmipp
  //Euler_angles2matrix(img.rot(), -img.tilt(), img.psi(), A);
  //A = A.inv();

  N2 = ImgNx * ImgNy;

  radius2 = diameter / 2.;
  radius2 = radius2 * radius2;
  dim2 = ImgNx / 2;
  
  for (i = 0; i < ImgNz; i++)
    {
	  z = -i + dim2;   /*** Z points upwards ***/
	  z2 = z * z;
	  for (j = 0; j < ImgNy; j++)
		{
		  y = j - dim2;
		  y2 = y * y;
		  z2_plus_y2 = z2 + y2;
		  x = 0 - dim2;   /***** X for k == 0 *****/
		  //xp = x * A(0, 0) + y * A(1, 0) + z * A(2, 0) + dim2;
		  //yp = x * A(0, 1) + y * A(1, 1) + z * A(2, 1) + dim2;
		  
		  xp = x * invrotmat[0] + y * invrotmat[3] + z * invrotmat[6] + dim2;
		  yp = x * invrotmat[1] + y * invrotmat[4] + z * invrotmat[7] + dim2;
		  
		  for (k = 0; k < ImgNx; k++, xp += invrotmat[0] , yp += invrotmat[1], x++)
            {
			  x2 = x * x;
			  //if (x2 + z2_plus_y2 > radius2)
			  //  continue;
			  if ((xp >= (ImgNx - 1) || xp < 0) || (yp >= (ImgNy - 1) || yp < 0))
				continue;
			  
			  /**** interpolation ****/
			  l = (int)yp;
			  m = (int)xp;
			  scalex = xp - m;
			  scaley = yp - l;
			  scale1 = 1. - scalex;
			  value1 = scalex * oneprj[(m+1)*ImgNx + l] + scale1 * oneprj[m*ImgNx + l];
			  value2 = scalex * oneprj[(m+1)*ImgNx + l+1] + scale1 * oneprj[m*ImgNx + l+1];
			  value  = scaley * value2 + (1. - scaley) * value1;
	
			  index = k * N2 + j * ImgNz + i;
			  
			  vol->data[index] += value;
			  
			  
			}
		}
	}
}


void Reconstruction::backprojection(Oimage *vol)
{

  //  Oimage*   vol ;
  float *oneprj, invrotmat[9], diameter=0.0;
  
  int    i,j,k, size;

  size = ImgNx * ImgNy;

  vol->nx = vol->ny = vol->nz = ImgNx;

  oneprj = (float *)malloc(size * sizeof(float));

  //vol    = InitImageParameters(1, ImgNx,ImgNy,ImgNz, 1);


 for ( i = 0; i < nv; i++ )
   {
	 
	 for ( j = 0; j < 9; j++ )
	   invrotmat[j] = InvRmat[i*9+j];


	 for ( j = 0; j < size; j++ ) 
		{
		  oneprj[j] = gd->data[i*size + j] ;

		}

	 simple_backprojection(oneprj, invrotmat, vol, diameter);

   }

 
 free(oneprj); oneprj = NULL;
//return vol;
  
}


float Reconstruction::MeanError(float*data1  , float* data0)
{

  float sum = 0.0;
  int   i;
  for (  i = 0 ; i < img2dsize; i++ )
    {
      sum = sum + (data1[i]-data0[i])*(data1[i]-data0[i]);

    }

  return sum;

}



float Reconstruction::GlobalMeanError(Oimage* object)
{
 PrjImg *Gd1 = (PrjImg *)malloc(sizeof(PrjImg));
 Gd1->data = (float *)malloc(nv * (ImgNx+1) * (ImgNy+1) * sizeof(float));

 //PrjImg *Gd0 = (PrjImg *)malloc(sizeof(PrjImg));
 //Gd0->data = (float *)malloc(nv * (ImgNx+1) * (ImgNy+1) * sizeof(float));

 printf("\ntolnv=%d ", nv);

 //EulerAngles* euler;// = new EulerAngles();
//euler = readFiles(filename, path, Gd0);


 //float *rmat =  (float *)malloc(tolnv *9*sizeof(float));

 //EulerMatrice(rmat,euler);

 int   i, j, v, size, sample_num;
 float rotmat[9];


/*
 for( i = 0; i < nv; i++ ){
   printf("\neuler=%f %f %f ", euler->rot, euler->tilt, euler->psi);

   euler=euler->next;
 }
*/
 sample_num = ImgNx+1;
 size = sample_num * sample_num;

 float *data1 = (float *)malloc(size *sizeof(float));
 float *data0 = (float *)malloc(size *sizeof(float));

 for ( i = 0; i < size; i++ )
  {
    data1[i] = 0.0;
    data0[i] = 0.0;   
  }

 for ( v = 0; v < nv; v++ )
   {
     for ( i = 0; i < 9; i++ )
       {rotmat[i] = Rmat[v*9+i];
        //printf("\nrotmat=%f ", rotmat[i]);
       }
      
	 //memset(data1, 0, size*sizeof(float));
     for ( i = 0; i < size; i++ ) data1[i] = 0.0;
     Volume_Projection(object, rotmat,  sample_num, data1);
     //Volume_GridProjection(object, rotmat,  sample_num, data1);
 
     for ( i = 0; i <= ImgNx; i++ )
	   {// printf("\n");
		 for ( j = 0; j <= ImgNx; j++ )
		   {
			 Gd1->data[v*size + i*(ImgNx+1) + j] = data1[i*(ImgNx+1) + j];			 //Gd0->data[v*size + i*(ImgNx+1) + j] = data0[i*(ImgNx+1) + j];
			 
			// printf("xdf=%f  ",Gd1->data[v*size + i*(ImgNx+1) + j] );
		   }
	   }
	 
   }
 
 

 float glbalError = 0.0, temp;

 for (  v = 0; v < nv; v++ )
   {
     for (  i = 0; i <= ImgNx; i++ )
	   {
		 
	   for ( j = 0; j <= ImgNx; j++ )
		 {
		   data1[i*(ImgNx+1) + j] = Gd1->data[v*size+i*(ImgNx+1) + j];
		   data0[i*(ImgNx+1) + j] = gd->data[v*img2dsize+i*(ImgNx+1) + j];
		   
		 }
	   }

     temp = MeanError(data1  , data0);
     //printf("\nv=%d meanerror=%e ", v, temp);
     glbalError = glbalError + temp;
   }

 printf("\nglbalError=%e ", glbalError);
 glbalError = sqrt(glbalError)/nv;

 printf("\nGlobal Mean Error ================================= %20.15e \n\n", glbalError);

 free(Gd1->data); 
 //free(Gd0->data);
 
 free(Gd1); Gd1 = NULL;
 //free(Gd0);
 free(data1); data1 = NULL;
 free(data0);  data0 = NULL;
 //free(rmat);
// delete euler;

 return glbalError;
}

/*
// old .
void Reconstruction::VolBsplineProjections()
{
  int i, j, k, t, s, v, proj_length, start_point0[2],start_point[2], id, id1, ii, jj, N2, proj_size, scale;
  float rotmat[9], *prjimg, *prjimg_sub, e1[3], e2[3];
  int sub, proj_length_sub, proj_size_sub;
  sub         = 2;

  scale       = (int)Bscale;   

  if(sub == 2)
	proj_length_sub = sub*8*scale-1;

  if (sub == 1)
	proj_length_sub = 8 * scale;

  proj_size_sub = proj_length_sub * proj_length_sub; 

  proj_length   = 8 * scale;
  proj_size     = proj_length * proj_length; 

  proj_VolB = (float *)malloc(nv*usedtolbsp*proj_size*sizeof(float));
  prjimg    = (float *)malloc(proj_size*sizeof(float));
  prjimg_sub= (float *)malloc(proj_size_sub*sizeof(float));

  startp    = (int    *)malloc(nv*usedtolbsp*2*sizeof(int)); 

  N2        = (usedN+1)*(usedN+1);

  for ( v = 0; v < nv; v++ )
	{
	  id  = v * usedtolbsp*2;
	  id1 = v * usedtolbsp*proj_size;
	  for ( i = 0; i < 9; i++ )
		{
		  rotmat[i] = Rmat[9*v+i];//printf("\nrotmat=%f ", rotmat[i1]);
		}
	  
	  
	  for (i = 0; i <3; i++ )
		{
		  e1[i] = rotmat[i];
		  e2[i] = rotmat[3+i];
		}
	  //printf("\nv=%d ", v);
	  //first compute Phi_000.
	  bspline->Bspline_Projection(0, 0, 0, sub, rotmat,ImgNx+1, 0.0, prjimg_sub, start_point0);
	  ii = (0-BStartX)*N2 + (0-BStartX)*(usedN+1)+0-BStartX;
	  	
	  for ( jj = s = 0; s < proj_length; s++ )
		{
		  //printf("\n");
		  for ( t = 0; t < proj_length; t++, jj++ )
			{
			  proj_VolB[id1+ii*proj_size+jj]=prjimg_sub[sub*s*proj_length_sub+sub*t]; 
			  prjimg[jj] = prjimg_sub[sub*s*proj_length_sub+sub*t]; 
			  
			  //printf("%f ", prjimg_sub[sub*s*proj_length_sub+sub*t]);
			}
		}

	  startp[id+2*ii+0] = start_point0[0];
	  startp[id+2*ii+1] = start_point0[1];
			    
	
	  for ( i = BStartX; i <= BFinishX; i++ )
	  	for ( j = BStartX; j <= BFinishX; j++ )
		  for( k = BStartX; k <= BFinishX; k++ )
			{ 
			  ii = (i-BStartX)*N2 + (j-BStartX)*(usedN+1)+k-BStartX;
	
			  		  
			  // old projection methods.			
			  //bspline->Bspline_Projection(i, j, k, sub, rotmat,ImgNx+1, 0.0, prjimg_sub, start_point);
			  //for ( s = 0; s < proj_size; s++ )
			//	proj_VolB[id1+ii*proj_size+s]=prjimg_sub[s]; 

			  //startp[id+2*ii+0] = start_point[0];
			  //startp[id+2*ii+1] = start_point[1];			  
			   //end old projection.
			   
			  			  	  
			  //quick projection methods by shift.
			  if(i != 0 || j != 0 || k != 0 )
				{ 
				  ComputePhi_ijk(i, j, k, e1, e2, start_point0,start_point, proj_length, sub, prjimg, prjimg_sub);

				for ( s = 0; s < proj_size; s++ )
				  proj_VolB[id1+ii*proj_size+s]=prjimg[s]; 

				startp[id+2*ii+0] = start_point[0];
				startp[id+2*ii+1] = start_point[1];
			  
				}
			  //end quick projection.
			  
			  
			  // startp[id+2*ii+0] = start_point[0];
			  // startp[id+2*ii+1] = start_point[1];
			  
			  
			} 
	}


  //orthognalize Vol Bspline projections.
  //getchar();
  free(prjimg);
  free(prjimg_sub);
}
*/


//Modified by Prof. Xu. 
void Reconstruction::VolBsplineProjections()
{
  int i, j, k, t, s, v, start_point0[2],start_point[2], id, id1, ii, jj, N2, proj_size, proj_size_sub, scale, tem;
  float rotmat[9], *prjimg, *prjimg_sub, e1[3], e2[3];

  //float *prjimg_sub_test; // Xu Added
  int   start_point0_test[2];

  //sub         = 1;
  //sub         = 2;
  //sub         = 4;
  //sub         = 8;
  SUB         = 16;
  //SUB         = 32;

  scale       = (int)Bscale;   

  /* Xu commoned out 
  if(sub == 2)
	proj_length_sub = sub*8*scale-1;

  if (sub == 1)
	proj_length_sub = 8 * scale;
  */

  // Xu added this line
  //proj_length_sub = sub*8*scale;
  // PRO_LENGTH and PRO_LENGTH_SUB need to be reconsidered, to make them as small as possible
  PRO_LENGTH_SUB = SUB*6*scale+1;   // Xu changed the box from [-4,4]^3 to [-3,3]^3
  proj_size_sub = PRO_LENGTH_SUB * PRO_LENGTH_SUB; 

  PRO_LENGTH   = 6 * scale + 1;     // Xu changed the box from [-4,4]^3 to [-3,3]^3
  proj_size     = PRO_LENGTH*PRO_LENGTH; 

  //proj_VolB = (float *)malloc(nv*usedtolbsp*proj_size*sizeof(float));
  proj_VolB = (float *)malloc(nv*proj_size_sub*sizeof(float));   // Xu changed to smaller data set
  //prjimg    = (float *)malloc(proj_size*sizeof(float));
  //prjimg_sub     = (float *)malloc(proj_size_sub*sizeof(float));
  //prjimg_sub_test= (float *)malloc(proj_size_sub*sizeof(float));
  
  //startp    = (int    *)malloc(nv*usedtolbsp*2*sizeof(int)); 
  startp    = (int    *)malloc(nv*2*sizeof(int));    // Xu changed 

  printf("Which is large %d, %d, %d, %d, nv, bs, proj_size = %d, %d, %d\n", nv*usedtolbsp*proj_size, proj_size, proj_size_sub,nv*usedtolbsp, nv, usedtolbsp, proj_size);
  N2        = (usedN+1)*(usedN+1);

   // Loop for projection directions
   for ( v = 0; v < nv; v++ ) {
      id  = v * usedtolbsp*2;
      id1 = v * usedtolbsp*proj_size;
      for ( i = 0; i < 9; i++ ) {
        rotmat[i] = Rmat[9*v+i];//printf("\nrotmat=%f ", rotmat[i1]);
      }
	  
      // coordinate directions	  
      for (i = 0; i <3; i++ ) {
	  e1[i] = rotmat[i];
	  e2[i] = rotmat[3+i];
      }
      ///printf("\nv==============================================%d, proj_length_sub=%d ,proj_length=%d", 
         ///        v,proj_length_sub,proj_length);

      //printf("\ne1 = %f, %f, %f, e2 = %f, %f, %f", e1[0], e1[1],e1[2],e2[0],e2[1],e2[2]);

      //first compute Phi_000. which is in prjimg_sub on fine grid, start_point0 = [-3,3] curently
      //bspline->Bspline_Projection(0, 0, 0, sub, rotmat,ImgNx+1, 0.0, prjimg_sub, start_point0);
      bspline->Bspline_Projection(0, 0, 0, SUB, rotmat,ImgNx+1, 0.0, proj_VolB+v*proj_size_sub, startp+v*2); // Xu changed

/*   // Xu comment    
    
      ii = (0-BStartX)*N2 + (0-BStartX)*(usedN+1)+0-BStartX;
	  	
      startp[id+2*ii+0] = start_point0[0];
      startp[id+2*ii+1] = start_point0[1];
			    
      //printf("BStartX = %d, BFinishX = %d, scale = %d \n", BStartX, BFinishX, scale);
	
      for ( i = BStartX; i <= BFinishX; i++ )
         for ( j = BStartX; j <= BFinishX; j++ )
	    for( k = BStartX; k <= BFinishX; k++ ) { 
               ii = (i-BStartX)*N2 + (j-BStartX)*(usedN+1)+k-BStartX;
	
               //quick projection methods.
	       //if (i != 0 || j != 0 || k != 0 ) {  // Xu commoned out this if
                  ComputePhi_ijk(i, j, k, e1, e2, start_point0,start_point, sub, prjimg, prjimg_sub);

		  for ( s = 0; s < proj_size; s++ ) {
		     proj_VolB[id1+ii*proj_size+s]=prjimg[s]; // First compute 
                     //printf("\n New 000image %d error= %f ", s, prjimg_sub[sub*s]-prjimg[s]);
                  }
		  startp[id+2*ii+0] = start_point[0];
		  startp[id+2*ii+1] = start_point[1];
			  
	       //}
	       //end quick projection.
            } 
  */
   }  // end v loop


  //orthognalize Vol Bspline projections.
  
  /*
  free(prjimg);
  free(prjimg_sub);
  */
}




void Reconstruction::VolBsplineProjections_FA()
{
  int i, j, k, t, s, v, proj_length, start_point[2], id, id1, ii, N2, proj_size, scale;
  float rotmat[9], *prjimg, *prjimg000;
  int    proj_size_000, sub;

  sub           = 2;

  scale         = (int)Bscale;   
  proj_length   = 8 * scale;
  proj_size     = proj_length * proj_length;
  proj_size_000 = (2*proj_length-1)*(2*proj_length-1);

	  proj_VolB   = (float *)malloc(nv*usedtolbsp*proj_size*sizeof(float));
	  prjimg      = (float *)malloc(proj_size*sizeof(float));
  prjimg000   = (float *)malloc(proj_size_000*sizeof(float));

  startp    = (int    *)malloc(nv*usedtolbsp*2*sizeof(int)); 

  N2        = (usedN+1)*(usedN+1);

  for ( v = 0; v < nv; v++ )
	{
	  id  = v * usedtolbsp*2;
	  id1 = v * usedtolbsp*proj_size;
	  for ( i = 0; i < 9; i++ )
		{
		  rotmat[i] = Rmat[9*v+i];//printf("\nrotmat=%f ", rotmat[i1]);
		}
	  
	  for ( i = BStartX; i <= BFinishX; i++ )
	  	for ( j = BStartX; j <= BFinishX; j++ )
		  for( k = BStartX; k <= BFinishX; k++ )
			{ 
			  ii = (i-BStartX)*N2 + (j-BStartX)*(usedN+1)+k-BStartX;

			  if( i == 0 && j== 0 && k == 0 )			 
				bspline->Bspline_Projection(i, j, k, sub, rotmat,ImgNx+1, 0.0, prjimg000, start_point);



			  //bspline->Bspline_GridProjection(i, j, k,  rotmat,ImgNx+1, 0.0, prjimg, start_point);

			  for ( s = 0; s < proj_size; s++ )
				proj_VolB[id1+ii*proj_size+s]=prjimg[s]; 
			  /*
			  for ( s = 0; s <= ImgNx; s++ )
				{//printf("\n");
				for (t = 0; t <= ImgNx; t++ )
				  {
					printf("\ncha====%e ", prjimg[s*proj_length+t]-gd->data[v*img2dsize+s*(ImgNx+1)+t]);
				  }
				}
			  getchar();
			  */
 			  startp[id+2*ii+0] = start_point[0];
			  startp[id+2*ii+1] = start_point[1];
			  
			  
			} 
	}


  //orthognalize Vol Bspline projections.
  //getchar();
  free(prjimg);
}



void Reconstruction::ComputeXdf (float* coef)
{
}

void Reconstruction::ComputeFFT_gd()
{
}

void Reconstruction:: subdivision_gd()
{
}


void Reconstruction::testFFT_sub()
{
  float *x, step, wk, weight, step_sub, real, im, a, b;
  int i, half, nx, subnx, sub, k,  subhalf;
  fftw_complex *in, *out, *conj_out, *sample_out;
 fftw_plan    test_fft;
  
 sub = 2;
 nx =  100;
 subnx = nx * sub;
 subhalf = subnx/2;
 step = 0.01;

 half = nx/2;
 a = -half*step; 
 b =  half*step;


 step_sub = step*1.0/sub;

 weight = (b-a)*1.0/subnx;



 x = (float *)malloc((subnx+1) * sizeof(float));

 for (i = -subhalf; i <= subhalf; i++ )
    x[i+subhalf] = exp(step_sub*i);


in = (fftw_complex*)fftw_malloc(subnx * sizeof(fftw_complex));
out = (fftw_complex*)fftw_malloc(subnx * sizeof(fftw_complex));
 conj_out =  (fftw_complex*)fftw_malloc((subnx-1) * sizeof(fftw_complex));
sample_out = (fftw_complex*)fftw_malloc(subnx * sizeof(fftw_complex));




for ( i = 0; i < subnx; i++ ) {
   in[i][0] = x[i]; in[i][1] = 0.0;
   out[i][0] = 0.0; out[i][1] = 0.0;
}

 
test_fft = fftw_plan_dft_1d(subnx,in,out,FFTW_FORWARD,FFTW_ESTIMATE);
fftw_execute(test_fft);

for ( k = 0; k < subnx; k++ ) 
   {
    wk = -PI2 * k/(b-a) * a ;
    out[k][0] = weight*(out[k][0]*cos(wk) - out[k][1] * sin(wk));
    out[k][0] = weight*(out[k][1]*cos(wk) + out[k][0] * sin(wk));

	//if(k != 0 )
	//{
	//conj_out[k][0] = out[k][0];
	//	conj_out[k][1] = -out[k][1];
	//}

   }

 for ( k = 0; k < subhalf; k++ ) 
   {
	 sample_out[subhalf+k][0] = out[k][0];
	 sample_out[subhalf+k][1] = out[k][1];

   }
 for ( k = 1; k < subhalf; k++ )
   {
	 sample_out[k][0] = sample_out[subnx-k][0];
	 sample_out[k][1] = -sample_out[subnx-k][1];
        
   } 
 sample_out[0][0] = out[subnx][0];
 sample_out[0][1] = -out[subnx][1];


 test_fft = fftw_plan_dft_1d(subnx,sample_out,out,FFTW_BACKWARD,FFTW_ESTIMATE);
 fftw_execute(test_fft);

 for ( k = 0; k < subnx; k++ )
   {
	 wk = M_PI*k;
	 wk = 0.0;
	 real = cos(wk)*out[k][0]-sin(wk)*out[k][1];
	 im   = cos(wk)*out[k][1]-sin(wk)*out[k][0];
 
     out[k][0] = 1.0/subnx * PI2/ImgNx * real;
     out[k][1] = 1.0/subnx * PI2/ImgNx * im;

   }	


for ( i = 0; i <= nx; i++ ) {
  printf("\ni = %d x=%f out=%f  %f ", i, x[2*i], out[half+i][0], out[half+i][1]);
}
/*

for ( i = 0; i <= half; i=i++) {
printf("\ni = %d out=%f  %f  ", i, weight*out[i][0], weight*out[i][1]);

}
*/
}


// This piece of code is modified by Xu for speed up
void Reconstruction::ComputeXdf_gd(float *f, float *coefs)
{
int i, j, v, v2, vp, ii, lx, rx, ly, ry, s, t, proj_imgnx, scale;
float *prjimg=NULL, *prjimg_ijk =NULL, sum, value, area, e1[3], e2[3];
int k, N2, proj_size, proj_size_sub, half, halfI, usedN1;
int i1, j1, k1, tem1, tem2, start_point[2],v9, ImgNx1, BNB, lxpl;
//Oimage *vol = InitImageParameters(1, ImgNx+1,ImgNy+1,ImgNz+1, 1);  
  
scale = (int)Bscale;
//proj_imgnx = 8 * scale - 1;
proj_imgnx = PRO_LENGTH - 1;
proj_size   = PRO_LENGTH * PRO_LENGTH;

proj_size_sub = PRO_LENGTH_SUB * PRO_LENGTH_SUB;
area   = 1.0;
area   = area * 1.0/nv;

usedN1 = usedN+1;
N2     = usedN1*usedN1;

prjimg  = (float *)malloc(img2dsize *sizeof(float));
prjimg_ijk  = (float *)malloc(proj_size*sizeof(float));
// xdgdphi = (float *)malloc(usedtolbsp*sizeof(float));

ImgNx1 = ImgNx+1; 
half = ImgNx/2; 
halfI = half*ImgNx1+half;

//for ( i = 0; i < VolImgSize; i++ ) vol->data[i] = f[i];
//vol->nx   = ImgNx1;
//vol->ny   = ImgNx1;
//vol->nz   = ImgNx1;

BNB = BStartX*N2 + BStartX; 

//vol->ux = vol->uy = vol->uz = 1.0;

for ( v = 0; v < nv; v++ ) {
   v2 = v + v; 
   vp = v*proj_size_sub;

   for (i = 0; i <3; i++ ) {
      v9 = v*9;
      e1[i] = Rmat[v9+i];
      e2[i] = Rmat[v9+3+i];
   }

   //Project the volume in the given directions, which is X_df.
   //Volume_Projection(vol, rotmat,  ImgNx+1,  prjimg);    /// Direct projection, using volume data.
   Volume_Projection_FromBasis(v, ImgNx+1,  prjimg, coefs); // Direct projection, using the projection of basis
   //Volume_GridProjection(vol, rotmat,  ImgNx+1,  prjimg);// Direct projection using bilinear interpolation.

   //	  for ( i = 0; i < img2dsize; i++ )
   //		printf("\nxdf=%f  gd=%f ", prjimg[i], gd->data[v*img2dsize+i]);


   for ( i = 0; i < img2dsize; i++ ) {
      prjimg[i] = prjimg[i] - gd->data[v*img2dsize+i];
   }

   //id  = v * usedtolbsp * 2;
   //id1 = v * usedtolbsp*proj_size;

   for ( i = BStartX; i <= BFinishX; i++ ) {
      i1 = i*N2 - BNB;
      for ( j = BStartX; j <= BFinishX; j++ ) {
         j1 = (j-BStartX)*usedN1; 
         //for ( k = BStartX; k <= BFinishX; k++ ) { 
         for ( k = BStartZ; k <= BFinishZ; k++ ) {    // Xu changed 1004
            //ii = (i-BStartX)*N2 + (j-BStartX)*(usedN+1)+k-BStartX;
            ii = i1 + j1 + k;

            //ComputePhi_ijk_Bilinear(i, j, k, e1, e2, startp + v2, start_point, SUB, prjimg_ijk, proj_VolB+vp);
            ComputePhi_ijk_Constant(i, j, k, e1, e2, startp + v2, start_point, SUB, prjimg_ijk, proj_VolB+vp);

            //tem4 = id + ii + ii;
            lx = start_point[0];
            ly = start_point[1];

            rx = lx+proj_imgnx;
            if (rx > half) rx = half;

            ry = ly+proj_imgnx;
            if (ry > half) ry = half;

            lxpl = lx*PRO_LENGTH + ly; 
	    sum = 0.0;
	    for (s = lx; s <= rx; s++ ) {

               tem1 = s*ImgNx1+halfI;
               tem2 = s*PRO_LENGTH - lxpl;

	       for (t = ly; t <= ry; t++ ) {
                  sum = sum + prjimg[tem1+t] * prjimg_ijk[tem2+t];;

                  // Xu comment out 
                  /*
                  if( t == -half) {
                     printf("\n boundary value = %e",    value);
                     if (s == -half || s == half) {
                        sum = sum -  0.75 * value;
                     } else {
                        sum = sum -  0.5 * value;
                     }
                     continue;
                  }

                  if( t == half) {
                     printf("\n boundary value = %e",    value);
                     if (s == -half || s == half) {
                        sum = sum -  0.75 * value;
                     } else {
                        sum = sum -  0.5 * value;
                     }
                     continue;
                  }

                  if ((s == -half) || (s == half) ) {
                     sum = sum - 0.5 * value;
                  }
                  */

                  /*
	          if( (t == -half  &&  (s == -half || s == half)) || (t == half &&(s == -half || s == half)))
	             sum += 0.25 * value;
	          else if(s > -half && s < half && t > -half && t < half)
	             sum += value;
	          else
	             sum += 0.5 *value;	
                  */
               }
            }
            xdgdphi[ii] += sum * area;
         } // end k loop
      }    // end j loop
   }       // end i loop
}          // end v loop

  //getchar();
//  free(vol->data);
//  free(vol);
  free(prjimg);
  free(prjimg_ijk);
  
}

//-----------------------------------------------------------------------------
void Reconstruction::ComputeXdf_gd(int newnv, int *index, float *f, float *coefs)
{
int i, j, v, v2, vp, ii, lx, rx, ly, ry, s, t, proj_imgnx, scale, newv;
float *prjimg=NULL, *prjimg_ijk =NULL, sum, value, area, e1[3], e2[3];
int k, N2, proj_size, proj_size_sub, half, halfI, usedN1;
int i1, j1, k1, tem1, tem2, start_point[2],v9, ImgNx1, BNB, lxpl;
  
scale = (int)Bscale;
proj_imgnx = PRO_LENGTH - 1;
proj_size   = PRO_LENGTH * PRO_LENGTH;

proj_size_sub = PRO_LENGTH_SUB * PRO_LENGTH_SUB;
area   = 1.0;
area   = area * 1.0/newnv;

usedN1 = usedN+1;
N2     = usedN1*usedN1;

prjimg  = (float *)malloc(img2dsize *sizeof(float));
prjimg_ijk  = (float *)malloc(proj_size*sizeof(float));

ImgNx1 = ImgNx+1; 
half = ImgNx/2; 
halfI = half*ImgNx1+half;

BNB = BStartX*N2 + BStartX; 

for ( v = 0; v < newnv; v++ ) {
//   v2 = v + v; 
//   vp = v*proj_size_sub;
    
     newv = index[nv-newnv+v];
     v2 = 2*newv; 
     vp = newv * proj_size_sub;

   for (i = 0; i <3; i++ ) {
      //v9 = v*9;
      v9 = newv * 9;
      e1[i] = Rmat[v9+i];
      e2[i] = Rmat[v9+3+i];
   }

   //Project the volume in the given directions, which is X_df.
   //Volume_Projection_FromBasis(v, ImgNx+1,  prjimg, coefs); // Direct projection, using the projection of basis
   Volume_Projection_FromBasis(newv, ImgNx+1,  prjimg, coefs); // Direct projection, using the projection of basis


   for ( i = 0; i < img2dsize; i++ ) {
      prjimg[i] = prjimg[i] - gd->data[newv*img2dsize+i];
   }

   for ( i = BStartX; i <= BFinishX; i++ ) {
      i1 = i*N2 - BNB;
      for ( j = BStartX; j <= BFinishX; j++ ) {
         j1 = (j-BStartX)*usedN1; 
         //for ( k = BStartX; k <= BFinishX; k++ ) { 
         for ( k = BStartZ; k <= BFinishZ; k++ ) {    // Xu changed 1004
            //ii = (i-BStartX)*N2 + (j-BStartX)*(usedN+1)+k-BStartX;
            ii = i1 + j1 + k;

            //ComputePhi_ijk_Bilinear(i, j, k, e1, e2, startp + v2, start_point, SUB, prjimg_ijk, proj_VolB+vp);
            ComputePhi_ijk_Constant(i, j, k, e1, e2, startp + v2, start_point, SUB, prjimg_ijk, proj_VolB+vp);

            //tem4 = id + ii + ii;
            lx = start_point[0];
            ly = start_point[1];

            rx = lx+proj_imgnx;
            if (rx > half) rx = half;

            ry = ly+proj_imgnx;
            if (ry > half) ry = half;

            lxpl = lx*PRO_LENGTH + ly; 
	    sum = 0.0;
	    for (s = lx; s <= rx; s++ ) {

               tem1 = s*ImgNx1+halfI;
               tem2 = s*PRO_LENGTH - lxpl;

	       for (t = ly; t <= ry; t++ ) {
                  sum = sum + prjimg[tem1+t] * prjimg_ijk[tem2+t];;

               }
            }
            xdgdphi[ii] += sum * area;
         } // end k loop
      }    // end j loop
   }       // end i loop
}          // end v loop

  free(prjimg);
  free(prjimg_ijk);
  
}

//-----------------------------------------------------------------------------
void Reconstruction::ComputeEJ1(float *f, float *coefs, int *index)
{
int i, v;
float *prjimg=NULL, ww;

prjimg  = (float *)malloc(img2dsize *sizeof(float));

for ( v = 0; v < nv; v++ ) {
   //Project the volume in the given directions, which is X_df.
   //Volume_Projection_FromBasis(v, ImgNx+1,  prjimg, coefs); // Direct projection, using the projection of basis
   Volume_Projection_FromBasis_Simplify(v, ImgNx+1,  prjimg, coefs); // Direct projection, using the projection of basis

   EJ1_k1[v] = 0.0;
   for ( i = 0; i < img2dsize; i++ ) {
      ww = prjimg[i] - gd->data[v*img2dsize+i];
      EJ1_k1[v] += ww * ww;
   }
}    // end v loop

for ( i = 0; i < nv; i++ ) EJ1_k[i] = fabs(EJ1_k1[i] - EJ1_k[i]);  

Sortingfloat(EJ1_k, index, nv);

//for ( i = 0; i < nv; i++ )  printf("EJ1=%f , index = %d\n", EJ1_k[index[i]], index[i]);
printf("EJ1=%f , index = %d\n", EJ1_k[index[nv-1]], index[nv-1]);

for ( i = 0; i < nv; i++ ) EJ1_k[i] = EJ1_k1[i];

free(prjimg);
}


// Xu added this code for test J1,  0918
//-----------------------------------------------------------------------------
void Reconstruction::ComputeEJ1_Acurate(float *f, float *coefs)
{
int i, v;
float *prjimg=NULL, ww, totalE;

prjimg  = (float *)malloc(img2dsize *sizeof(float));

totalE = 0.0;
for ( v = 0; v < nv; v++ ) {
   //Project the volume in the given directions, which is X_df.
   Volume_Projection_FromBasis(v, ImgNx+1,  prjimg, coefs); // Direct projection, using the projection of basis
   //Volume_Projection_FromBasis_Simplify(v, ImgNx+1,  prjimg, coefs); // Direct projection, using the projection of basis

   for ( i = 0; i < img2dsize; i++ ) {
      ww = prjimg[i] - gd->data[v*img2dsize+i];
      totalE  = totalE + ww * ww;
   }
}    // end v loop
printf("Total EJ1 ============ %f\n", totalE);

free(prjimg);
}


//-----------------------------------------------------------------------------
//Compute Temporal step-size tau.   This code is modified by Xu 0918
void Reconstruction::ComputeTau(float *coefs, float *diff_coefs, float *tau,
                                float reconj1, float taube, int flows)
{
int i, v, a, b, c;
float *prjimg=NULL, *dprjimg =NULL, denominator, numerator;
float numerator1, denominator1, epsilon, tau0;
float lengthf, lengthdf, slength, inner;
  
// if (flows == 3 && reconj1 == 0.0) return; // do nothing

denominator = 0.0;
numerator = 0.0;
denominator1 = 0.0;
numerator1 = 0.0;

tau0 = 0.1;
epsilon = 0.00001;

if (reconj1 != 0.0) {
   prjimg  = (float *)malloc(img2dsize *sizeof(float));
   dprjimg  = (float *)malloc(img2dsize *sizeof(float));

   for ( v = 0; v < nv; v++ ) {

      Volume_Projection_FromBasis_2(v, ImgNx+1,  prjimg, coefs, dprjimg, diff_coefs); // Direct projection, using the projection of basis
      for ( i = 0; i < img2dsize; i++ ) {
         prjimg[i] = prjimg[i] - gd->data[v*img2dsize+i];
      }

      for ( i = 0; i < img2dsize; i++ ) {
         numerator += prjimg[i] * dprjimg[i];
         denominator += dprjimg[i] * dprjimg[i];

      }
   }       // end v loop

   free(prjimg);
   free(dprjimg);
   tau0 = -numerator/denominator;
   *tau = tau0;
 }
 
J3_Tau(tau0, coefs, diff_coefs, taube, flows, &numerator1, &denominator1);


//printf("fenzi = %f, %f, fenmu = %f, %f\n", numerator,  numerator1, denominator, denominator1);
denominator = denominator + denominator1;
if(fabs(denominator) <  epsilon) denominator = epsilon;
*tau = -(numerator + numerator1)/denominator;
  
//printf("tau0 = %f,  tau =  %f\n", tau0, *tau);
printf("tau0 = %f,  total tau =  %f, n d = %e, %e\n", tau0, *tau, numerator1,denominator1);

//if (*tau < 0.0) *tau = 0.0;

//if (*tau > tau0) *tau = tau0;       // Xu Test 0927

}

//-----------------------------------------------------------------------------
//Compute Temporal step-size tau.   This code is modified by Xu 0918
void Reconstruction::J3_Tau(float tau0, float *coefs, float *diff_coefs, float taube, 
                            int flows, float *numerator1, float *denominator1)
{
int  i, j, k, l, a, b, c, *index, size, size2;
float gradient[10], gradient1[10], epsilon;
float length, lengthf, lengthdf, slength, inner;
float *grad_leng_f, dg, ww, depsilon, depsilon2, depsilon_sq;

*denominator1 = 0.0;
*numerator1 = 0.0;
epsilon = 0.00001*0.00001;

// Xu added the following for J3,  0918
if (taube != 0.0 && flows == 1) {
 
   for (a = 1; a < ImgNx; a++) {
      for (b = 1; b < ImgNy; b++) {
         for (c = 1; c < ImgNz; c++) {

            B_spline_Function_Gradient_Grid(coefs, gradient, a, b, c);
            B_spline_Function_Gradient_Grid(diff_coefs, gradient1, a, b, c);

            // for minimizing \!\nabla f\|^2 
            *numerator1 = *numerator1 + InnerProduct(gradient, gradient1);
            *denominator1 = *denominator1 + InnerProduct(gradient1, gradient1);
            
         }     // end c loop
      }        // end b loop
   }           // end a loop
   return;
}              // end if 


// Xu added the following for J3, 0920
if (taube != 0.0 && flows == 2) {
   depsilon = 0.001;
   depsilon2 = depsilon + depsilon;
   depsilon_sq = depsilon*depsilon;
   dg = 0.0;
   for (a = 1; a < ImgNx; a++) {
      for (b = 1; b < ImgNy; b++) {
         for (c = 1; c < ImgNz; c++) {

            B_spline_Function_Gradient_Grid(coefs, gradient, a, b, c);
            B_spline_Function_Gradient_Grid(diff_coefs, gradient1, a, b, c);

            lengthf = InnerProduct(gradient, gradient);
            if(lengthf  <  epsilon) lengthf = epsilon;
            slength = sqrt(lengthf);


            inner    = InnerProduct(gradient, gradient1);
            lengthdf = InnerProduct(gradient1, gradient1);

            
            *numerator1 = *numerator1 + inner/slength;  
            /*
            ww = lengthdf/slength  - inner*inner/(lengthf*slength);
            *denominator1 = *denominator1 + ww;
            */

            dg = dg + (inner + depsilon*lengthdf)/sqrt(lengthf + depsilon2*inner + depsilon_sq*lengthdf);
//if (lengthdf > 10.0) printf("slength = %f, innr = %f, de = %f, %e, a, b, c = %d, %d, %d\n", slength, inner/slength, lengthdf, lengthdf, a , b, c);

            //denominator1 = denominator1 + lengthdf/slength;
         }     // end c loop
      }        // end b loop
   }           // end a loop

   *denominator1 = (dg - *numerator1)/depsilon;

   //printf("old fen mu = %f, new fenmu = %f\n",  *denominator1, ww);

   *numerator1 = 0.5* *numerator1;
   return;
}              // end if 


// Xu modified the following for J3, 0925
if (taube != 0.0 && flows == 3) {
   float  gradient[10],gradient1[10], length, tem, epsilon, grad, grad2, 
          normal[3], gradient2[3],rgrad, rgrad2, rgrad4,H, H2;
   float  Hf[9], HGf[3], sum, sum1;
   float Hnormal[3], H2normal[3], Hfn[3], nHfn, nHfn_normal[3], H_phi[9], nH_phi[3];
   float ww, c1, c2, c3, c4, c5, c6, c7, c8, c9;

   depsilon = 0.0001;
   sum = 0.0; 
   sum1 = 0.0; 
   for (a = 1; a < ImgNx; a++) {
      for (b = 1; b < ImgNy; b++) {
         for (c = 1; c < ImgNz; c++) {

            // at point epsilon = 0;
            B_spline_Function_Hessian_Grid(coefs, gradient, a, b, c); // Xu 0924 
            B_spline_Function_Hessian_Grid(diff_coefs, gradient1, a, b, c); // Xu 0924 

            gradient2[0] = gradient[1] * gradient[1];
            gradient2[1] = gradient[2] * gradient[2];
            gradient2[2] = gradient[3] * gradient[3];

            grad2 = gradient2[0]+ gradient2[1] + gradient2[2];
            if (grad2 < epsilon) grad2 = epsilon;

            grad = sqrt(grad2);
            rgrad = 1.0/grad;
            rgrad2 = 1.0/grad2;

            //Compute mean curvature H = - 0.5 * div(normal);
            H =   gradient[1]*(gradient[4]*gradient[1] + gradient[5]*gradient[2] + gradient[6]*gradient[3])
                + gradient[2]*(gradient[5]*gradient[1] + gradient[7]*gradient[2] + gradient[8]*gradient[3])
                + gradient[3]*(gradient[6]*gradient[1] + gradient[8]*gradient[2] + gradient[9]*gradient[3]);

            H = H*rgrad2 - gradient[4] - gradient[7] - gradient[9];
            H = 0.5 * H * rgrad;
            H2 = H * H;
            normal[0] = gradient[1] * rgrad ;
            normal[1] = gradient[2] * rgrad ;
            normal[2] = gradient[3] * rgrad ;
        
            // Hessian of f 
            Hf[0] = gradient[4];
            Hf[1] = gradient[5];
            Hf[2] = gradient[6];
            Hf[3] = gradient[5];
            Hf[4] = gradient[7];
            Hf[5] = gradient[8];
            Hf[6] = gradient[6];
            Hf[7] = gradient[8];
            Hf[8] = gradient[9];
 
            // H*n
            MatrixMultiply(Hf, 3,3, normal, 3, 1, Hfn);
  
            // Hn * 2H / |\nabla f\| 
            ww = rgrad*H*2;
            Hfn[0] = Hfn[0]*ww;
            Hfn[1] = Hfn[1]*ww;
            Hfn[2] = Hfn[2]*ww;

            // n^T * H * n
            nHfn  = InnerProduct(normal, Hfn);
 
            Hnormal[0] = normal[0] * H;
            Hnormal[1] = normal[1] * H;
            Hnormal[2] = normal[2] * H;

            H2normal[0] = normal[0] * H2;
            H2normal[1] = normal[1] * H2;
            H2normal[2] = normal[2] * H2;


            nHfn_normal[0] = normal[0] * nHfn;
            nHfn_normal[1] = normal[1] * nHfn;
            nHfn_normal[2] = normal[2] * nHfn;
          
            // Xu add these c's for speedup 0925
            c1 =  - H2normal[0] + Hfn[0] - nHfn_normal[0];
            c2 =  - H2normal[1] + Hfn[1] - nHfn_normal[1];
            c3 =  - H2normal[2] + Hfn[2] - nHfn_normal[2];

            c4 =   Hnormal[0]*normal[0];
            c5 = 2*Hnormal[0]*normal[1];
            c6 = 2*Hnormal[0]*normal[2];
            c7 =   Hnormal[1]*normal[1];
            c8 = 2*Hnormal[1]*normal[2];
            c9 =   Hnormal[2]*normal[2];

            sum = sum - H * (gradient1[4] +  gradient1[7] + gradient1[9]);
  
            sum = sum + c1*gradient1[1] + c2*gradient1[2] + c3*gradient1[3] 
                      + c4*gradient1[0] + c5*gradient1[5] + c6*gradient1[6]
                      + c7*gradient1[7] + c8*gradient1[8] + c9*gradient1[9]; 


            // at the other point 
            for (i = 0; i < 10; i++) {
                gradient[i] = gradient[i] + depsilon * gradient1[i];
            }

            gradient2[0] = gradient[1] * gradient[1];
            gradient2[1] = gradient[2] * gradient[2];
            gradient2[2] = gradient[3] * gradient[3];

            grad2 = gradient2[0]+ gradient2[1] + gradient2[2];
            if (grad2 < epsilon) grad2 = epsilon;

            grad = sqrt(grad2);
            rgrad = 1.0/grad;
            rgrad2 = 1.0/grad2;

            //Compute mean curvature H = - 0.5 * div(normal);
            H =   gradient[1]*(gradient[4]*gradient[1] + gradient[5]*gradient[2] + gradient[6]*gradient[3])
                + gradient[2]*(gradient[5]*gradient[1] + gradient[7]*gradient[2] + gradient[8]*gradient[3])
                + gradient[3]*(gradient[6]*gradient[1] + gradient[8]*gradient[2] + gradient[9]*gradient[3]);

            H = H*rgrad2 - gradient[4] - gradient[7] - gradient[9];
            H = 0.5 * H * rgrad;
            H2 = H * H;
            normal[0] = gradient[1] * rgrad ;
            normal[1] = gradient[2] * rgrad ;
            normal[2] = gradient[3] * rgrad ;
        
            // Hessian of f 
            Hf[0] = gradient[4];
            Hf[1] = gradient[5];
            Hf[2] = gradient[6];
            Hf[3] = gradient[5];
            Hf[4] = gradient[7];
            Hf[5] = gradient[8];
            Hf[6] = gradient[6];
            Hf[7] = gradient[8];
            Hf[8] = gradient[9];
 
            // H*n
            MatrixMultiply(Hf, 3,3, normal, 3, 1, Hfn);
  
            // Hn * 2H / |\nabla f\| 
            ww = rgrad*H*2;
            Hfn[0] = Hfn[0]*ww;
            Hfn[1] = Hfn[1]*ww;
            Hfn[2] = Hfn[2]*ww;

            // n^T * H * n
            nHfn  = InnerProduct(normal, Hfn);

            Hnormal[0] = normal[0] * H;
            Hnormal[1] = normal[1] * H;
            Hnormal[2] = normal[2] * H;

            H2normal[0] = normal[0] * H2;
            H2normal[1] = normal[1] * H2;
            H2normal[2] = normal[2] * H2;


            nHfn_normal[0] = normal[0] * nHfn;
            nHfn_normal[1] = normal[1] * nHfn;
            nHfn_normal[2] = normal[2] * nHfn;
          
            // Xu add these c's for speedup 0925
            c1 =  - H2normal[0] + Hfn[0] - nHfn_normal[0];
            c2 =  - H2normal[1] + Hfn[1] - nHfn_normal[1];
            c3 =  - H2normal[2] + Hfn[2] - nHfn_normal[2];

            c4 =   Hnormal[0]*normal[0];
            c5 = 2*Hnormal[0]*normal[1];
            c6 = 2*Hnormal[0]*normal[2];
            c7 =   Hnormal[1]*normal[1];
            c8 = 2*Hnormal[1]*normal[2];
            c9 =   Hnormal[2]*normal[2];

            sum1 = sum1 - H * (gradient1[4] +  gradient1[7] + gradient1[9]);

            sum1 = sum1 + c1*gradient1[1] + c2*gradient1[2] + c3*gradient1[3] 
                        + c4*gradient1[0] + c5*gradient1[5] + c6*gradient1[6]
                        + c7*gradient1[7] + c8*gradient1[8] + c9*gradient1[9]; 

         }     // end c loop
      }        // end b loop
   }           // end a loop

   *numerator1 = sum;
   *denominator1 = (sum1 - sum)/depsilon;
   
   if (sum > 0.0)  {
      printf("Warning, The Energy is not decreasing at this moment\n"); // 0927 
   }
   if (*denominator1 < 0)  {
      printf("Warning, nagative denominator1\n");
      *denominator1 = - sum/depsilon;
   }

    printf("Flow = 3, N sum1 = %e, D = %e, %e, local tau = %f \n", *numerator1, sum1, *denominator1, - *numerator1/ *denominator1 );
    *numerator1 = 0.5* *numerator1;
}

}


//-----------------------------------------------------------------------------
//Compute Temporal step-size tau.
// coefs, diff_coefs are non-orthogonal 
void Reconstruction::ComputeTau(int newnv, int *index, float *coefs, float *diff_coefs, 
                                float *tau, float reconj1, float taube, int flows)
{
int i, v, a, b, c, newv;
float denominator, numerator, numerator1, denominator1;
float *prjimg=NULL, *dprjimg =NULL, epsilon, length, lengthf, lengthdf, slength, inner, tau0;
  
//if (flows == 3 && reconj1 == 0.0) return; // do nothing
 
denominator = 0.0;
numerator = 0.0;

tau0 = 0.1;
epsilon = 0.00001;

if (reconj1 != 0.0) {
   prjimg  = (float *)malloc(img2dsize *sizeof(float));
   dprjimg  = (float *)malloc(img2dsize *sizeof(float));

   for ( v = 0; v < newnv; v++ ) {

      newv = index[nv-newnv+v];
      Volume_Projection_FromBasis_2(newv, ImgNx+1,  prjimg, coefs, dprjimg, diff_coefs); 

      for ( i = 0; i < img2dsize; i++ ) {
         prjimg[i] = prjimg[i] - gd->data[newv*img2dsize+i];
      }

      for ( i = 0; i < img2dsize; i++ ) {
         numerator += prjimg[i] * dprjimg[i];
         denominator += dprjimg[i] * dprjimg[i];

      }
   }       // end v loop
   free(prjimg);
   free(dprjimg);
   tau0 = -numerator/denominator;
   *tau = tau0;
}

J3_Tau(tau0, coefs, diff_coefs, taube, flows, &numerator1, &denominator1);

printf("fenzi = %f, %f, fenmu = %f, %f\n", numerator,  numerator1, denominator, denominator1);
denominator = denominator + denominator1;
if(fabs(denominator) <  epsilon) denominator = epsilon;
*tau = -(numerator + numerator1)/denominator;
  
printf("tau0 = %f,  tau =  %f, n d = %e, %e\n", tau0, *tau, numerator1,denominator1);

//if (*tau < 0.0) *tau = 0.0;
//if (*tau > tau0) *tau = tau0;  // Xu Test 0927
}


//------------------------------------------------------------------------------
void Reconstruction::ComputeIFSFf_gd(const float* f, float* padf, int padfactor, int nx, int padnx, int sub, float fill_value, 
                                     float* Non_o_coef, fftw_complex *in, fftw_complex *out)
{
  int i, j, k, ii, id, jd, wki, wkj;
  int subImgNx, subimg2dsize, subhalf;
  fftw_plan fft_Volf;
  float ai, bj, wk;


//pad volume f.
if (padfactor == 2 )
{
  ImgPad(f, padf, nx, padnx,fill_value);
  FFT_padf(padf, padnx);
}

else if(padfactor == 1) 
{
  //sub = 1;
  subImgNx     = sub * ImgNx;
  subimg2dsize = (subImgNx+1) * (subImgNx+1); 
  subhalf      = subImgNx/2;

  
  // transfer data to in for FFT.
  for ( ii = i = 0; i < subImgNx; i++ ) 
	{
	  id  = i * subimg2dsize;
	  wki = M_PI * i;
	  
	  for ( j = 0; j < subImgNx; j++ )
		{
		  jd  = j *(subImgNx+1);
		  wkj = M_PI*j;
		  
		  for ( k = 0;  k < subImgNx; k++, ii++ )
			{
			  wk = wki + wkj + M_PI*k;
                          wk = M_PI*(i+j+k);
			 			  
			  in[ii][0] = f[id+jd+k] * cos(wk); 
			  in[ii][1] = f[id+jd+k] * sin(wk);
			  //printf("\nin=%f %f", in[ii][0], in[ii][1]);		  
			}
		}
	}
  
  //getchar();
  //printf("\nsubImgNx=%d subhalf=%d", subImgNx, subhalf);
  //getchar();
  fft_Volf =  fftw_plan_dft_3d(subImgNx, subImgNx, subImgNx, in, out,FFTW_FORWARD,FFTW_ESTIMATE);
  
  fftw_execute(fft_Volf);
  
  for ( i = -subhalf; i < subhalf; i++ ) 
     {
      id  = (i+subhalf)*subImgNx*subImgNx;
      //wki = PI*i;

      for ( j = -subhalf; j < subhalf; j++ )
		{
          jd  = (j+subhalf)*subImgNx;
          //wkj = PI*j;
		  
		  for ( k = -subhalf;  k < subhalf; k++ )
			{
			  ii = id + jd + k + subhalf;
			  //wk = wki + wkj + PI*k;   //increase the errors.

			  wk = M_PI*(i+j+k);
			  
			  ai = cos(wk)*out[ii][0]-sin(wk)*out[ii][1];
			  bj = cos(wk)*out[ii][1]+sin(wk)*out[ii][0];
			  
			  //out[ii][0] = 0.125*ai; // for subimage.
			  //out[ii][1] = 0.125*bj;
			  out[ii][0] = ai; 
			  out[ii][1] = bj;
			  //if(k == 0 ) printf("\ni  j k =%d %d %d ii=%d sub=%d out=%f %f ", i, j, k, ii, sub, out[ii][0], out[ii][1]);
			}
		}
	 }
  //getchar();
  fftw_destroy_plan(fft_Volf);

  GetCentralSlice(out, Non_o_coef, sub); 
}

}


void Reconstruction::ImgPad(const float *f, float* padf, int size, int padsize, float fill_value)
{
  int padVolImgSize,padimg2dsize,  i, j, k, id, id0, idi, idj, id0i, id0j;

  padVolImgSize = padsize * padsize * padsize;
  padimg2dsize  = padsize * padsize;

  
  for (i = 0; i < padsize; i++ )
     {
        idi = i*padimg2dsize;
        id0i = i*img2dsize;

	for ( j = 0; j  < padsize; j++ )
           {
            idj = j * padsize;
            id0j = j * size;

	    for ( k = 0; k < padsize; k++ )
		{
		  id  = idi + idj + k;

		  if(i >=  size || j >= size || k >= size )// size = ImgNx+1. for odd volume.
			padf[id] = fill_value;
		  else
			{
			  id0 = id0i + id0j + k;
			  padf[id] = f[id0]; 

			}

		  // printf("\ni j k = %d %d %d padf=%f ", i, j, k, padf[id]);
		}
           }
      }
    
  //for ( i = 0; i < padVolImgSize; i++ ) 
  //printf("\n sf=%f ",  padf[i]);

}

void Reconstruction::FFT_padf(const float* padf, int padsize)
{
  int padVolImgSize, v, i, j, k, x, y, hx, hy, ix, iy, iz, padhalf, half;
  float rotmat[9], e1[3], e2[3], X[3], d2, inv, wk, real, im;
  int id, id1, ii, lx, rx, ly, ry, s, t;
  float sum, area, value;
  int    proj_length, proj_size, N2, scale, proj_imgnx;

  scale = (int)Bscale;

  proj_length = 8*scale;
  proj_imgnx  = 8 * scale -1;

  proj_size   = proj_length * proj_length;
  N2          = (usedN+1)*(usedN+1);
  area        = 1.0*1.0;


  padVolImgSize = padsize * padsize *padsize;
  padhalf = padsize/2;
  inv     = 1.0/(2*ImgNx);
  half    = ImgNx/2;

  fftw_plan fft_padf, inv_slice;
  fftw_complex *in, *out, *slice, *out1;

  //in   =  (fftw_complex *) padf;


  in   = (fftw_complex *)fftw_malloc(padVolImgSize*sizeof(fftw_complex));
  out  = (fftw_complex *)fftw_malloc(padVolImgSize*sizeof(fftw_complex));
  slice = (fftw_complex *)fftw_malloc(ImgNx*ImgNx*sizeof(fftw_complex));
  out1 = (fftw_complex *)fftw_malloc(ImgNx*ImgNx*sizeof(fftw_complex));

   
 for (v = i = 0; i < padsize; i++ ) 
	for (j = 0; j < padsize; j++ )
	  for ( k = 0; k < padsize; k++ , v++ )
		{
		  wk = M_PI * ( i + j + k);

		  in[v][0] = padf[v]*cos(wk);
		  in[v][1] = padf[v]*sin(wk); 
		  //printf("\ni j k = %d %d %d in=%f %f  padf=%f ", i, j, k, in[v][0], in[v][1], padf[v]);
		}

 //getchar();

  fft_padf = fftw_plan_dft_3d(padsize,padsize, padsize, in, out,FFTW_FORWARD,FFTW_ESTIMATE);
  fftw_execute(fft_padf);

  fftw_destroy_plan(fft_padf);
  
  ImgPhaseShift2Origin(out, padsize);
  
  //for (i = 0; i < padVolImgSize; i++ ) printf("\nout=% f %f ", out[i][0], out[i][1]);
  //getchar();
  for (i = 0; i < usedtolbsp; i++ ) xdgdphi[i] = 0.0;
 
  for ( v = 0; v < nv; v++ )
	{
	  for ( i = 0; i < 9; i++ )
		rotmat[i] = Rmat[v*9+i];
	  
	  for (i = 0; i <3; i++ )
		{
		  e1[i] = rotmat[i];
		  e2[i] = rotmat[3+i];
		}
	  for ( i = 0; i < ImgNx*ImgNx; i++ ) 
		{
		  slice[i][0] = 0.0;
		  slice[i][1] = 0.0;
		}

	  for ( x = -half; x < half; x++ )
		{
		  //hx = x;
		  //if ( hx > padhalf ) hx -= padsize;
		  for ( y = -half;  y < half; y++ )
			{
			  i = (x+half)*ImgNx+y+half;
			  //hy = y;
			  //if ( hy > padhalf ) hy -= padsize;

			  X[0] = x * e1[0] + y * e2[0];
			  X[1] = x * e1[1] + y * e2[1];
			  X[2] = x * e1[2] + y * e2[2];

			  //d2 = inv * inv * (X[0]*X[0] + X[1]*X[1] + X[2]*X[2]);

			  //if(d2 < 0.25 )
			  //{
				  iz = floor(2*X[2] + 0.5);   // Nearest neighbour
				  iy = floor(2*X[1] + 0.5);
				  ix = floor(2*X[0] + 0.5);

				  iz += padhalf;
				  iy += padhalf;
				  ix += padhalf;

				  j = ix*padsize*padsize + iy*padsize + iz;
				  //printf("\nix iy iz = %d %d  %d ", ix, iy, iz);
				  if(ix >= 0 && ix < 2*ImgNx && iy >= 0 && iy < 2*ImgNx && iz >= 0 && iz < 2*ImgNx)
					{
					  slice[i][0] = out[j][0];
					  slice[i][1] = out[j][1];
					}
				  
				  
				  //}

			}
		}

	  //for ( i = 0; i < img2dsize; i++ ) printf("\nslice=%f %f ", slice[i][0], slice[i][1]);
	  //getchar();

	  ImgPhaseShift2Origin2D(slice, ImgNx);
	  //for ( i = 0; i < img2dsize; i++ ) printf("\nslice=%f %f ", slice[i][0], slice[i][1]);
	  //getchar();
	  fft_padf = fftw_plan_dft_2d(ImgNx,ImgNx,slice,out1,FFTW_BACKWARD,FFTW_ESTIMATE);
	  fftw_execute(fft_padf);
	  
	  for (x = i = 0; i < ImgNx; i++ )
		for ( j = 0;  j < ImgNx; j++, x++ )
		  {
			//ii = i*ImgNx+j;
			wk = -M_PI*(i+j);
			
			real = cos(wk)*out1[x][0] - sin(wk)*out1[x][1];
			im   = cos(wk)*out1[x][1] + sin(wk)*out1[x][0];

			out1[x][0] = 1.0/(ImgNx*ImgNx)*real;
			out1[x][1] = 1.0/(ImgNx*ImgNx)*im;
 
		  }
	 /* 
	  for ( x = i = 0; i < ImgNx; i++ )
		for (j = 0; j < ImgNx; j++, x++ )
		  printf("\ninverse=%f  gd=%f  ", out1[x][0], gd->data[v*img2dsize+i*(ImgNx+1)+j]);
	  getchar();
	  */

	  for ( x = i = 0; i < ImgNx; i++ )
		for ( j = 0; j < ImgNx; j++, x++ )
		  out1[x][0] = out1[x][0] - gd->data[v*img2dsize+i*(ImgNx+1)+j];
	  

	  id  = v * usedtolbsp * 2;
      id1 = v * usedtolbsp*proj_size;
	  
	  for ( i = BStartX; i <= BFinishX; i++ )
	  	for ( j = BStartX; j <= BFinishX; j++ )
		  for( k = BStartX; k <= BFinishX; k++ )
			{ 
			  ii = (i-BStartX)*N2 + (j-BStartX)*(usedN+1)+k-BStartX;
			  lx = startp[id+2*ii+0];
			  ly = startp[id+2*ii+1];
			  rx = ((lx+proj_imgnx)>=half)?(half-1):(lx+proj_imgnx); 
			  ry = ((ly+proj_imgnx)>=half)?(half-1):(ly+proj_imgnx);
			  
			  sum = 0.0;
			  //if(lx <-4 || ly < -4 ) {
			  //printf("\nlx rx ly ry=%d %d %d %d ", lx, rx, ly, ry);getchar();}
			  for(s = lx; s <= rx; s++ )
				for(t = ly; t <= ry; t++ )
				  {
					//value = out[(s+half)*(ImgNx+1)+t+half][0] * proj_VolB[id1+ii*proj_size+(s-lx)*proj_length+t-ly];
					value = out1[(s+half)*ImgNx+t+half][0] * proj_VolB[id1+ii*proj_size+(s-lx)*proj_length+t-ly];

					//printf("\nout=%f  proj_img=%f ", out[(s+half)*(ImgNx+1)+t+half][0], proj_VolB[id1+ii*proj_size+(s-lx)*proj_length+t-ly]);

					if( (t == -half  &&  (s == -half || s == half-1)) || (t == half-1 &&(s == -half || s == half-1)))
					  sum += 0.25 * value;
					else if(s > -half && s < half-1 && t > -half && t < half-1)
					  sum += value;
					else
					  sum += 0.5 *value;
					
				  }
			  //printf("sum=%f ", sum);
			  xdgdphi[ii] += sum * area * 1.0/nv;
			  //printf("xgdphi=%f ", xdgdphi[ii]);
			}
			  
			  

	}
   
  fftw_destroy_plan(fft_padf);
	  
  fftw_free(in);
  fftw_free(slice);
  fftw_free(out);
  fftw_free(out1);

}


void Reconstruction::ImgPhaseShift2Origin(fftw_complex* ftdata, int padsize)
{
  int x, y, z, l, k, h, i, padhalf;
  float shift, sl, skl, w, real, im;

  shift = 0.25;//(ImgNx+1)/2.0;
  padhalf = padsize/2;


  for ( x= -padhalf; x < padhalf; x++) 
	{
	  //l = x;
	  //if ( x > padhalf ) l -= padsize;
	  sl = x*shift;

	  for ( y=-padhalf; y < padhalf; y++) 
		{
		  //k = y;
		  //if ( y > padhalf ) k -= padsize;
		  skl = sl + y*shift;
		  for ( z=-padhalf; z < padhalf; z++ )
			{
			  i = (x+padhalf)*padsize*padsize+(y+padhalf)*padsize+z+padhalf;
			  //h = z;
			  //if( z > padhalf ) h -= padsize;
			  w = PI2 * (skl + z*shift);
 
			  real = cos(w)*ftdata[i][0]-sin(w)*ftdata[i][1];
			  im   = cos(w)*ftdata[i][1]+sin(w)*ftdata[i][0];

			  ftdata[i][0] = real;
			  ftdata[i][1] = im;

			}
		}
	}




}

void Reconstruction::ImgPhaseShift2Origin2D(fftw_complex* ftdata, int padsize)
{
  int x, y, z, l, k, h, i, padhalf;
  float shift, sl, skl, w, real, im;

  shift = -0.5;//-(ImgNx+1)/2.0;
  padhalf = padsize/2;


  for ( x=-padhalf; x < padhalf; x++) 
	{
	  //l = x;
	  //if ( x > padhalf ) l -= padsize;
	  sl = x*shift;

	  for ( y=-padhalf; y < padhalf; y++) 
		{
		  i = (x+padhalf)*padsize+(y+padhalf);
		  //k = y;
		  //if ( y > padhalf ) k -= padsize;
		  skl = sl + y*shift;
		  
		  w = PI2 * skl;
 
		  real = cos(w)*ftdata[i][0]-sin(w)*ftdata[i][1];
		  im   = cos(w)*ftdata[i][1]+sin(w)*ftdata[i][0];
		  
		  ftdata[i][0] = real;
		  ftdata[i][1] = im;
		  
		  
		}
	}
  



}


void Reconstruction::Test_GridProjectioin(EulerAngles* eulers)
{
Oimage* p = InitImageParameters(1, ImgNx,ImgNy,ImgNz, 1);  //6.6  change to ImgNx +1.
  
 float OX, OY, OZ, R, length, width, x, y, z;
 int i, j, k, N2, index;

 OX = (ImgNx)/2.0;
 OY = OZ = OX;

 R = (ImgNx-2)/2.0;
 N2   = ImgNy * ImgNz;
 
 length = (ImgNx-1)/2.0;
 
 length = length/2.0;
 length = length * length;
 width  = sqrt(R*R - length); 
 
 for ( i = 0; i < ImgNx; i++ )
   for ( j = 0; j < ImgNy; j++ )
	 for ( k = 0; k < ImgNz; k++ )
	   {
		 x = i-OX;
		 y = j-OY;  
		 z = k-OZ;  
		 

		 index = i*N2 + j*ImgNz + k;
		 if( (x*x + y*y) <= length && fabs(z)<= width) 
		   p->data[index] = 1.0;
		 else 
		   p->data[index] = 0.0;

		 //printf("\ni j k = %d %d %d f=%f ", i, j, k, p->data[index]);
	   }

 p->nx = p->ny = p->nz = ImgNx;

 // do projections.
 int size = ImgNx * ImgNx;

 float* prjimg = (float *)malloc(size *sizeof(float));
 int v;
 float rotmat[9];
 EulerAngles* home = NULL;
 char filename[100];
 home = eulers;

 for ( v = 0; v < nv; v++ )
   {
	 for ( i = 0; i < 9; i++ )
	   {
		 rotmat[i] = Rmat[v*9+i];
				 
	   }
	 //Volume_Projection(p, rotmat,  ImgNx, prjimg);
	 Volume_GridProjection(p, rotmat,  ImgNx, prjimg);
	 for ( k = i = 0; i < ImgNx; i++ )
	   {
		 printf("\n");
		 for ( j = 0; j < ImgNx; j++ , k++)
		   printf("%f ", prjimg[k]);
	   }
	 getchar();

	 //if(v+1 < 10 )  sprintf(filename, "column800000%d.spi",v+1);
	 //else sprintf(filename, "column80000%d.spi",v+1);
	 
	 //WriteSpiFile(filename, prjimg, ImgNx, ImgNx, 1, home->rot, home->tilt, home->psi);
	 //home = home->next;

   }


free(prjimg); prjimg = NULL;



}





/*
//old.
void Reconstruction::ComputePhi_ijk(int i, int j, int k, float e1[3], float e2[3], int start_point0[2], int start_point[2], int proj_length, int sub, float* prjimg, float* prjimg_sub)
{
 
  int i1, j1, ii, lx, rx, ly, ry, x0, y0, half;
  float s1, s2, x, y, xf, yf, r1, r2, r3, r4;
  int    proj_length_sub, size;

  s1 = i * e1[0] + j * e1[1] + k * e1[2];
  s2 = i * e2[0] + j * e2[1] + k * e2[2];

  //printf("\nstartpoint 0 = %d %d e1=%f %f %f e2=%f %f %f", start_point0[0], start_point0[1], e1[0], e1[1], e1[2], e2[0], e2[1], e2[2]);
  lx = (int)floor(start_point0[0] + s1);
  ly = (int)floor(start_point0[1] + s2);
  
  half = ImgNx*0.5;
  proj_length_sub = sub * proj_length - 1;
 

  lx = (-half <= lx)?lx:(-half); //image range: [-half, half]X[-half, half];
  ly = (-half <= ly)?ly:(-half); 

  rx = lx + proj_length - 1;
  ry = ly + proj_length - 1;


  rx = (rx>half)?half:rx; 
  ry = (ry>half)?half:ry;
			  
  size = proj_length * proj_length;
  for ( i1 = 0; i1 < size; i1++ )	prjimg[i1] = 0.0;

  // printf("\ni j k = %d %d %d lx rx ly ry = %d %d %d %d half=%d", i, j, k, lx, rx, ly, ry, half);
  for ( i1 = 0; i1 <= rx-lx; i1++ ) 
	{
	  //printf("\n");
	  for ( j1 = 0; j1 <= ry-ly; j1++ )
		{
		  x = lx + i1 - s1;
		  y = ly + j1 - s2;
		  //printf("\nold x y = %f %f ", x, y);
		  x = x + x;
		  y = y + y;
		  
		  x0 = (int)floor(x);
		  y0 = (int)floor(y);
		
		  xf = x - x0;
		  yf = y - y0;
		  
		  //printf("\ni1 j1 = %d %d xf yf = %f %f 2x=%f 2y=%f s1 s2=%f %f ", i1, j1, xf, yf, x, y, s1, s2);
		  x0 = ((x0-2*start_point0[0])>=0)?(x0-2*start_point0[0]):0;
		  y0 = ((y0-2*start_point0[1])>=0)?(y0-2*start_point0[1]):0;
		  //printf("x0 y0 = %d %d ", x0, y0);
		  
		  r1 = (x0>=proj_length_sub   || y0>=proj_length_sub)?0.0:prjimg_sub[x0*proj_length_sub+y0];
		  r2 = (x0>=proj_length_sub-1 || y0>=proj_length_sub)?0.0:prjimg_sub[(x0+1)*proj_length_sub+y0];
		  r3 = (x0>=proj_length_sub   || y0>=proj_length_sub-1)?0.0:prjimg_sub[x0*proj_length_sub+y0+1];
		  r4 = (x0>=proj_length_sub-1 || y0>=proj_length_sub-1)?0.0:prjimg_sub[(x0+1)*proj_length_sub+y0+1];
		  
		  //printf("r1 r2 r3 r4=%f %f %f %f ", r1, r2, r3, r4);
		  prjimg[i1*proj_length+j1] = 
			(1-xf)*(1-yf)*r1+xf*(1-yf)*r2+(1-xf)*yf*r3+xf*yf*r4;    

		  //printf("%f ", prjimg[i1*proj_length+j1]);
		}
	}
  //getchar();

  start_point[0] = lx;
  start_point[1] = ly;

}
*/

// Modified by Prof. Xu.
// prjimg_sub = X_d\Phi_{000}, it is on fine grid, which is input
// prjimg is the result after translation, its size is proj_length^2
// start_point0 is low-left corner of its domain
void Reconstruction::ComputePhi_ijk_Bilinear(int i, int j, int k, float e1[3], float e2[3], int start_point0[2], 
                                int start_point[2], int sub, float* prjimg, float* prjimg_sub)
{
 
int i1, j1, ii, kk, lx, rx, ly, ry, x0, y0, half;
float s1, s2, x, y, xf, yf, xf1, yf1, r1, r2, r3, r4, lx1, ly1, xi;
int    x0p, x0py0, x0py01;

s1 = Bscale*(i * e1[0] + j * e1[1] + k * e1[2]);
s2 = Bscale*(i * e2[0] + j * e2[1] + k * e2[2]);

//printf("\nstartpoint 0 = %d %d e1=%f %f %f e2=%f %f %f", start_point0[0], start_point0[1], e1[0], e1[1], e1[2], e2[0], e2[1], e2[2]);
lx = (int)floor(start_point0[0] + s1);
ly = (int)floor(start_point0[1] + s2);
  
half = ImgNx*0.5;

lx = (-half <= lx)?lx:(-half); //image range: [-half, half]X[-half, half];
ly = (-half <= ly)?ly:(-half); 

rx = lx + PRO_LENGTH - 1;
ry = ly + PRO_LENGTH - 1;

rx = (rx>half)?half:rx; 
ry = (ry>half)?half:ry;
			  
lx1 = lx - s1;
ly1 = ly - s2;

for ( i1 = 0; i1 < PRO_LENGTH; i1++ ) {

   //printf("\n");
   x = lx1 + i1;
   x = sub*x;  // Xu added
   x0 = floor(x);   // for bilinear 
   xf = x - x0;
   x0 = x0-sub*start_point0[0];
   kk = i1*PRO_LENGTH; 

   if (x0 >= 0 && x0 < PRO_LENGTH_SUB-1) {
      x0p = x0*PRO_LENGTH_SUB;
      for ( j1 = 0; j1 < PRO_LENGTH; j1++ ) {

         y = ly1 + j1;
         y = sub*y;
         y0 = floor(y);  // for bilinear
         yf = y - y0;
         y0 = y0-sub*start_point0[1];

         if (y0 >= 0 && y0 < PRO_LENGTH_SUB-1) {

            // bilinear interpolation
            // Four neighbour values
            /*
            r1 = (x0>=PRO_LENGTH_SUB   || y0>=PRO_LENGTH_SUB)?0.0:prjimg_sub[x0p+y0];
            r2 = (x0>=PRO_LENGTH_SUB-1 || y0>=PRO_LENGTH_SUB)?0.0:prjimg_sub[x0p+PRO_LENGTH_SUB+y0];
            r3 = (x0>=PRO_LENGTH_SUB   || y0>=PRO_LENGTH_SUB-1)?0.0:prjimg_sub[x0p+y0+1];
            r4 = (x0>=PRO_LENGTH_SUB-1 || y0>=PRO_LENGTH_SUB-1)?0.0:prjimg_sub[x0p+PRO_LENGTH_SUB+y0+1];
            */
            
            x0py0 = x0p+y0;
            x0py01 = x0p+y0+1;
            r1 = prjimg_sub[x0py0];
            r2 = prjimg_sub[x0py0 + PRO_LENGTH_SUB];
            r3 = prjimg_sub[x0py01];
            r4 = prjimg_sub[x0py01 + PRO_LENGTH_SUB];
            
		  
            //printf("r1 r2 r3 r4=%f %f %f %f ", r1, r2, r3, r4);
            // Bilinear Interpolation

            //prjimg[kk + j1] = (1-xf)*(r1 + yf*(r3 - r1)) + xf*(r2 + yf*(r4 - r2));    // method 1, middle 
            //prjimg[kk + j1] = r1 + yf*(r3 - r1) + xf*(r2 - r1 + yf*(r4 - r2 - r3 + r1) );// method 2 slowest 
            prjimg[kk + j1] = (1-xf)*(1-yf)*r1+xf*(1-yf)*r2+(1-xf)*yf*r3+xf*yf*r4;    // method 3, more operations, but fastest, ???!!!
            //xf1 = 1 - xf;
            //yf1 = 1 - yf;
            //prjimg[kk + j1] = xf1*yf1*r1+xf*yf1*r2+xf1*yf*r3+xf*yf*r4;                  // method 4, less operations, but not fastest, ???!!!
      
         }  else {
            prjimg[kk + j1] = 0.0;
            //printf("\n, out of the reange-----------------------");
         }

      }
   }  else {
      for ( j1 = 0; j1 < PRO_LENGTH; j1++ ) {
         prjimg[kk + j1] = 0.0;
         //printf("\n, out of the reange-----------------------");
      }
   }

}

start_point[0] = lx;
start_point[1] = ly;

  //printf("\n Start0 = %d, %d, Start = %d, %d\n", start_point0[0],start_point0[1], lx, ly);
}


// Modified by Prof. Xu.
// prjimg_sub = X_d\Phi_{000}, it is on fine grid, which is input
// prjimg is the result after translation, its size is proj_length^2
// start_point0 is low-left corner of its domain
void Reconstruction::ComputePhi_ijk_Constant(int i, int j, int k, float e1[3], float e2[3], int start_point0[2], 
                                int start_point[2], int sub, float* prjimg, float* prjimg_sub)
{
 
int   i1, j1, ii, kk, lx, rx, ly, ry, x0, y0, half;
float s1, s2, x, y, xf, yf, xf1, yf1, r1, r2, r3, r4, lx1, ly1, xi;
int   x0p, x0py0, x0py01,Y0[10], subp0, subp1;

s1 = Bscale*(i * e1[0] + j * e1[1] + k * e1[2]);
s2 = Bscale*(i * e2[0] + j * e2[1] + k * e2[2]);

lx = floor(start_point0[0] + s1);
ly = floor(start_point0[1] + s2);
  
half = ImgNx*0.5;
/*
lx = (-half <= lx)?lx:(-half); //image range: [-half, half]X[-half, half];
ly = (-half <= ly)?ly:(-half); 
*/
if(lx < -half) lx = -half;
if(ly < -half) ly = -half;

rx = lx + PRO_LENGTH - 1;
ry = ly + PRO_LENGTH - 1;

/*
rx = (rx>half)?half:rx; 
ry = (ry>half)?half:ry;
*/
if(rx>half) rx = half;
if(ry>half) ry = half;

lx1 = lx - s1;
ly1 = ly - s2;

subp0 = sub*start_point0[0];
subp1 = sub*start_point0[1];

for ( j1 = 0; j1 < PRO_LENGTH; j1++ ) {
   y = ly1 + j1;
   y = sub*y + 0.5;
   y0 = floor(y);
   //Y0[j1] = y0-sub*start_point0[1];
   Y0[j1] = y0-subp1;
}
			  
//kk = 0;
for ( i1 = 0; i1 < PRO_LENGTH; i1++ ) {

   x = lx1 + i1;
   x = sub*x + 0.5; 
   x0 = floor(x);    
   x0 = x0-subp0;
   kk = i1*PRO_LENGTH; 


   if (x0 >= 0 && x0 < PRO_LENGTH_SUB-1) {
      x0p = x0*PRO_LENGTH_SUB;
      for ( j1 = 0; j1 < PRO_LENGTH; j1++ ) {
         y0 = Y0[j1];
         if (y0 >= 0 && y0 < PRO_LENGTH_SUB-1) {
            prjimg[kk + j1] = prjimg_sub[x0p+y0];
         }  else {
            prjimg[kk + j1] = 0.0;
          }
      }
   }  else {
      for ( j1 = 0; j1 < PRO_LENGTH; j1++ ) {
         prjimg[kk + j1] = 0.0;
      }
   }   // end if 

   //kk = kk + PRO_LENGTH;
}      // end i1 loop

start_point[0] = lx;
start_point[1] = ly;

}

// start_point0 is low-left corner of its domain
void Reconstruction::ComputePhi_ijk_Constant_Simplify(int i, int j, int k, float e1[3], float e2[3], int start_point0[2], 
                                int start_point[2], int sub, float* prjimg, float* prjimg_sub)
{
 
int   i1, j1, ii, kk, lx, ly, x0, y0, half;
float s1, s2, x, y, xf, yf, xf1, yf1, r1, r2, r3, r4, lx1, ly1, xi;
int   x0p, x0py0, x0py01,Y0[10], subp0, subp1;

s1 = Bscale*(i * e1[0] + j * e1[1] + k * e1[2]);
s2 = Bscale*(i * e2[0] + j * e2[1] + k * e2[2]);

lx = floor(start_point0[0] + s1);
ly = floor(start_point0[1] + s2);
  
half = ImgNx*0.5;
if(lx < -half) lx = -half;
if(ly < -half) ly = -half;

lx1 = lx - s1;
ly1 = ly - s2;

subp0 = sub*start_point0[0];
subp1 = sub*start_point0[1];

//for ( j1 = 0; j1 < PRO_LENGTH; j1++ ) {
for ( j1 = 2; j1 < PRO_LENGTH - 2; j1++ ) {
   y = ly1 + j1;
   y = sub*y + 0.5;
   y0 = floor(y);
   Y0[j1] = y0-subp1;
}
			  
//for ( i1 = 0; i1 < PRO_LENGTH; i1++ ) {
for ( i1 = 2; i1 < PRO_LENGTH -2; i1++ ) {

   x = lx1 + i1;
   x = sub*x + 0.5; 
   x0 = floor(x);    
   x0 = x0-subp0;
   kk = i1*PRO_LENGTH; 


   if (x0 >= 0 && x0 < PRO_LENGTH_SUB-1) {
      x0p = x0*PRO_LENGTH_SUB;
      //for ( j1 = 0; j1 < PRO_LENGTH; j1++ ) {
      for ( j1 = 2; j1 < PRO_LENGTH - 2; j1++ ) {
         y0 = Y0[j1];
         if (y0 >= 0 && y0 < PRO_LENGTH_SUB-1) {
            prjimg[kk + j1] = prjimg_sub[x0p+y0];
         }  
      }
   }  
}      // end i1 loop

start_point[0] = lx + 2;
start_point[1] = ly + 2;

}

/*----------------------------------------------------------------------------*/
void Reconstruction::Volume_Projection_FromBasis(int v, int sample_num, float *prjimg, float *coefs)
{
int i, j, k, ii, jj, kk, lx, rx, ly, ry, s, t, jd, proj_imgnx, scale;
int  proj_size, proj_size_sub, half, halfI, start_point[2];
int  usedN1,  usedN2, tt, tem1, tem2, tem3, ImgNx1, v2, v9, vp, lxpl;
float   e1[3], e2[3], d[3], a[3], rotmat[9],*prjimg_ijk =NULL;

scale = (int)Bscale;
//proj_imgnx = 8 * scale - 1;
proj_imgnx = PRO_LENGTH - 1;
proj_size     = PRO_LENGTH * PRO_LENGTH;
proj_size_sub = PRO_LENGTH_SUB * PRO_LENGTH_SUB;
half = ImgNx/2;
ImgNx1 = ImgNx+1;
halfI = half*ImgNx1 + half;


prjimg_ijk  = (float *)malloc(proj_size*sizeof(float));

//usedN = N - 4;
usedN1 = usedN+1;
usedN2 = usedN1*usedN1;

//id1 = v * usedtolbsp*proj_size;
//id  = v * usedtolbsp*2;

v2 = v + v;
v9 = v*9;
vp = v*proj_size_sub;

for ( i = 0; i < img2dsize; i++ ) {
   prjimg[i] = 0.0;
}

//for ( i = 0; i < 9; i++ ) {
//   rotmat[i] = Rmat[v9+i];
//}
// coordinate directions
for (i = 0; i <3; i++ ) {
   e1[i] = Rmat[v9 + i];
   e2[i] = Rmat[v9 + 3 + i];
}

// Speed up
for ( i = BStartX; i <= BFinishX; i++ ) {
   ii = (i-BStartX)*usedN2 - BStartX;
   for ( j = BStartX; j <= BFinishX; j++ ) {
     jj = (j-BStartX)*usedN1;
     //for ( k = BStartX; k <= BFinishX; k++ ) {
     for ( k = BStartZ; k <= BFinishZ; k++ ) {    // Xu changed 1004
         kk = ii+ jj + k;

         //ComputePhi_ijk_Bilinear(i, j, k, e1, e2, startp + v2, start_point, SUB, prjimg_ijk, proj_VolB+vp);
         ComputePhi_ijk_Constant(i, j, k, e1, e2, startp + v2, start_point, SUB, prjimg_ijk, proj_VolB+vp);

         //ii = (i-BStartX)*usedN2 + (j-BStartX)*(usedN+1)+k-BStartX;
         lx = start_point[0];
         ly = start_point[1];


         rx = lx+proj_imgnx;
         if (rx > half) rx = half;

         ry = ly+proj_imgnx;
         if (ry > half) ry = half; 

         lxpl = lx*PRO_LENGTH + ly;

         for (s = lx; s <= rx; s++ ) {

            tem2 = s*PRO_LENGTH - lxpl;
            tem3 = s*ImgNx1 + halfI;

            for (t = ly; t <= ry; t++ ) {
               tt = tem3 + t;
               prjimg[tt] = prjimg[tt] + coefs[kk] *  prjimg_ijk[tem2+t]; // Second use of proj_VolB
            }
         }
      }
   }
}



/*
// Old
for ( i = BStartX; i <= BFinishX; i++ ) {
  //id = (i-BStartX)*usedN2;
  for ( j = BStartX; j <= BFinishX; j++ ) {
      //jd = (j-BStartX)*(usedN+1);
      for ( k = BStartX; k <= BFinishX; k++ ) {
         //jj = id+ jd + k-BStartX;

         ii = (i-BStartX)*usedN2 + (j-BStartX)*(usedN+1)+k-BStartX;
         lx = startp[id+2*ii+0];
         ly = startp[id+2*ii+1];
         rx = ((lx+proj_imgnx)>half)?half:(lx+proj_imgnx);
         ry = ((ly+proj_imgnx)>half)?half:(ly+proj_imgnx);

         for (s = lx; s <= rx; s++ ) {
            for (t = ly; t <= ry; t++ ) {
               tt = (s+half)*(ImgNx+1)+t+half;
               prjimg[tt] = prjimg[tt] + coefs[ii]*proj_VolB[id1+ii*proj_size+(s-lx)*proj_length+t-ly];
            }
         }
      }
   }
}
*/

free(prjimg_ijk);

}


/*----------------------------------------------------------------------------*/
void Reconstruction::Volume_Projection_FromBasis_Simplify(int v, int sample_num, float *prjimg, float *coefs)
{
int i, j, k, ii, jj, kk, lx, rx, ly, ry, s, t, jd, proj_imgnx, scale;
int  proj_size_sub, half, halfI, start_point[2], half_i;
int  usedN1,  usedN2, tt, tem1, tem2, tem3, ImgNx1, v2, v9, vp, lxpl;
float   e1[3], e2[3], d[3], a[3], rotmat[9],*prjimg_ijk =NULL, s1, s2;

scale = (int)Bscale;
proj_imgnx = PRO_LENGTH - 1;
//proj_imgnx = PRO_LENGTH - 5;
proj_size_sub = PRO_LENGTH_SUB * PRO_LENGTH_SUB;
half = ImgNx/2;
ImgNx1 = ImgNx+1;
halfI = half*ImgNx1 + half;
half_i = (BStartX + BFinishX)/2;


prjimg_ijk  = (float *)malloc(PRO_LENGTH* PRO_LENGTH*sizeof(float));

for ( i = 0; i < PRO_LENGTH* PRO_LENGTH; i++ ) {
   prjimg_ijk[i] = 0.0;
}

usedN1 = usedN+1;
usedN2 = usedN1*usedN1;

v2 = v + v;
v9 = v*9;
vp = v*proj_size_sub;

for ( i = 0; i < img2dsize; i++ ) {
   prjimg[i] = 0.0;
}

// coordinate directions
for (i = 0; i <3; i++ ) {
   e1[i] = Rmat[v9 + i];
   e2[i] = Rmat[v9 + 3 + i];
}

 ComputePhi_ijk_Constant(half_i, half_i, half_i, e1, e2, startp + v2, start_point, SUB, prjimg_ijk, proj_VolB+vp);

start_point[0] = startp[v2];
start_point[1] = startp[v2+1];

// Speed up
for ( i = BStartX; i <= BFinishX; i++ ) {
   ii = (i-BStartX)*usedN2 - BStartX;
   for ( j = BStartX; j <= BFinishX; j++ ) {
     jj = (j-BStartX)*usedN1;
     //for ( k = BStartX; k <= BFinishX; k++ ) {
     for ( k = BStartZ; k <= BFinishZ; k++ ) {   // Xu Changed 1004
         kk = ii+ jj + k;

         //ComputePhi_ijk_Constant_Simplify(i, j, k, e1, e2, startp + v2, start_point, SUB, prjimg_ijk, proj_VolB+vp);

         //lx = start_point[0];
         //ly = start_point[1];

         
         s1 = Bscale*(i * e1[0] + j * e1[1] + k * e1[2]);
         s2 = Bscale*(i * e2[0] + j * e2[1] + k * e2[2]);

         lx = floor(start_point[0] + s1);
         ly = floor(start_point[1] + s2);
         

         rx = lx+proj_imgnx;
         if (rx > half) rx = half;

         ry = ly+proj_imgnx;
         if (ry > half) ry = half; 

         if(lx < -half) lx = -half;
         if(ly < -half) ly = -half;

         lxpl = lx*PRO_LENGTH + ly;

         //for (s = lx; s <= rx; s++ ) {
         //for (s = lx+1; s < rx; s++ ) {
         for (s = lx+2; s < rx-1; s++ ) {

            tem2 = s*PRO_LENGTH - lxpl;
            tem3 = s*ImgNx1 + halfI;

            //for (t = ly; t <= ry; t++ ) {
            //for (t = ly+1; t < ry; t++ ) {
            for (t = ly+2; t < ry-1; t++ ) {
               tt = tem3 + t;
               prjimg[tt] = prjimg[tt] + coefs[kk] *  prjimg_ijk[tem2+t]; // Second use of proj_VolB
            }
         }
      }
   }
}

free(prjimg_ijk);

}

/*----------------------------------------------------------------------------*/
void Reconstruction::Volume_Projection_FromBasis_2(int v, int sample_num, float *prjimg, float *coefs, float *prjimg1, float *coefs1)
{
int i, j, k, ii, jj, kk, lx, rx, ly, ry, s, t, jd, proj_imgnx, scale;
int  proj_size, proj_size_sub, half, halfI, start_point[2];
int  usedN1,  usedN2, tt, tem1, tem2, tem3, ImgNx1, v2, v9, vp, lxpl;
float   e1[3], e2[3], d[3], a[3], rotmat[9],*prjimg_ijk =NULL;

scale = (int)Bscale;
//proj_imgnx = 8 * scale - 1;
proj_imgnx = PRO_LENGTH - 1;
proj_size     = PRO_LENGTH * PRO_LENGTH;
proj_size_sub = PRO_LENGTH_SUB * PRO_LENGTH_SUB;
half = ImgNx/2;
ImgNx1 = ImgNx+1;
halfI = half*ImgNx1 + half;

prjimg_ijk  = (float *)malloc(proj_size*sizeof(float));

usedN1 = usedN+1;
usedN2 = usedN1*usedN1;

v2 = v + v;
v9 = v*9;
vp = v*proj_size_sub;

for ( i = 0; i < img2dsize; i++ ) {
   prjimg[i] = 0.0;
   prjimg1[i] = 0.0;
}

// coordinate directions
for (i = 0; i <3; i++ ) {
   e1[i] = Rmat[v9 + i];
   e2[i] = Rmat[v9 + 3 + i];
}

// Speed up
for ( i = BStartX; i <= BFinishX; i++ ) {
   ii = (i-BStartX)*usedN2 - BStartX;
   for ( j = BStartX; j <= BFinishX; j++ ) {
     jj = (j-BStartX)*usedN1;
     //for ( k = BStartX; k <= BFinishX; k++ ) {
     for ( k = BStartZ; k <= BFinishZ; k++ ) { // Xu changed 1004
         kk = ii+ jj + k;

         //ComputePhi_ijk_Bilinear(i, j, k, e1, e2, startp + v2, start_point, SUB, prjimg_ijk, proj_VolB+vp);
         ComputePhi_ijk_Constant(i, j, k, e1, e2, startp + v2, start_point, SUB, prjimg_ijk, proj_VolB+vp);

         //ii = (i-BStartX)*usedN2 + (j-BStartX)*(usedN+1)+k-BStartX;
         lx = start_point[0];
         ly = start_point[1];

         rx = lx+proj_imgnx;
         if (rx > half) rx = half;

         ry = ly+proj_imgnx;
         if (ry > half) ry = half; 

         lxpl = lx*PRO_LENGTH + ly;

         for (s = lx; s <= rx; s++ ) {

            tem2 = s*PRO_LENGTH - lxpl;
            tem3 = s*ImgNx1 + halfI;

            for (t = ly; t <= ry; t++ ) {
               tt = tem3 + t;
               prjimg[tt] = prjimg[tt] + coefs[kk] *  prjimg_ijk[tem2+t]; // Second use of proj_VolB
               prjimg1[tt] = prjimg1[tt] + coefs1[kk] *  prjimg_ijk[tem2+t]; // Second use of proj_VolB
            }
         }
      }
   }
}

free(prjimg_ijk);

}


// Cimpoute J_3, see book page 72, example 3.5.1
/*----------------------------------------------------------------------------*/
void  Reconstruction::ComputeJ3_HTF(float *coefs, float taube, float *J234phi)
{
int i, j, k, i1, j1, k1, ii, jj, kk, aa, bb, cc, abc, ilx, jly, klz;
int a, b, c, al, bl, cl, N2, lx, rx, ly, ry, lz, rz, scale, size, usedN1;
float sum, gradient[10],length,partials[20], tem, epsilon, *partials_000;
float Bx[10], Bx1[10],  By[10], By1[10], Bz[10], Bz1[10], values[2], rscale;
float Bxy,Bxy1, Bx1y, Imatrix[9];

//printf("Begin Compute J3.\n");
epsilon = 0.000001;
//epsilon = 0.1;
scale   = (int)Bscale;
rscale  = 1.0/Bscale;
usedN1 = usedN+1;
N2 = (usedN+1)*(usedN+1);
//size = (4*scale-1)*(4*scale-1);

//partials_000 = (float*)malloc(size*(4*scale-1)*4*sizeof(float)); 

//new methods.
 for (a = 0; a <= ImgNx; a++) {
    for (b = 0; b <= ImgNy; b++) {
       for (c = 0; c <= ImgNz; c++) {

          B_spline_Function_Gradient_Grid(coefs, gradient, a, b, c); 

          lx = a/scale-3;
          rx = lx+3;
          if (lx < 0) lx = 0;
          if (rx > usedN) rx = usedN;
          if (lx > usedN ) lx = usedN;

          ly = b/scale-3;
          ry = ly+3;
          if (ly < 0) ly = 0;
          if (ry > usedN) ry = usedN;
          if (ly > usedN) ly = usedN;
          
          lz = c/scale-3;
          rz = lz+3;
          if (lz < 0) lz = 0;
          if (rz > usedN) rz = usedN;
          if (lz > usedN) lz = usedN;
 

          for ( i = lx; i<= rx; i++ ) {
             ilx = i - lx;
             bspline->Spline_N_Base_1(a-2-i, values);
             Bx[ilx]  = values[0];
             Bx1[ilx] = values[1]*rscale;
          }

          for ( j = ly; j <= ry; j++ ) {
             jly = j - ly;
             bspline->Spline_N_Base_1(b-2-j, values);
             By[jly]   = values[0];
             By1[jly]  = values[1]*rscale;
          }

          for ( k = lz; k <= rz; k++ ) {
             klz = k - lz;
             bspline->Spline_N_Base_1(c-2-k, values);
             Bz[klz]  = values[0];
             Bz1[klz] = values[1]*rscale;
          }

          for  (i = lx; i <= rx; i++) {
             ii = i*N2;
             ilx = i - lx; 
             for (j = ly; j <= ry; j++) {
                jj = j * usedN1;
                jly = j - ly;

                Bxy  = Bx[ilx]  * By[jly];
                Bx1y = Bx1[ilx] * By[jly];
                Bxy1 = Bx[ilx]  * By1[jly];


                for (k = lz; k <= rz; k++) {
                   klz = k - lz; 
                   kk = ii + jj + k;

                   partials[0] = Bx1y*Bz[klz];
                   partials[1] = Bxy1*Bz[klz];
                   partials[2] = Bxy *Bz1[klz];


                   /*
                   // using these for minimizing \|\nabla f \|, comment out these for heat flow
                   length = InnerProduct(gradient, gradient);
                   if( length < epsilon ) length = epsilon;
                   length = 1.0/sqrt(length);

                   partials[0] = partials[0]*length;
                   partials[1] = partials[1]*length;
                   partials[2] = partials[2]*length;
                   */

                   sum = InnerProduct(gradient, partials);
                   J234phi[kk] = J234phi[kk] +  taube*sum;

                }
             }
          }
       }
    }
 }


/*
//old methods.
 // Loop for B-spline basis
//for ( i = BStartX; i <= BFinishX; i++ ) {
for ( i = 0; i <= usedN; i++ ) {

   // the range of grid points 
   //aa = i - BStartX + 2;

   i1 = i*N2;
   lx = i *scale;  // ??? why 
   //rx = (lx + 4)*scale;
   rx = lx + 4*scale;

//   for ( j = BStartX; j <= BFinishX; j++ ) {
     for ( j = 0; j <= usedN; j++ ) {

      //bb = j - BStartX + 2;

      j1 = j*(usedN+1);
      ly = j *scale;
      //ry = (ly + 4)*scale;
      ry = ly + 4*scale;

      //for ( k = BStartX; k <= BFinishX; k++ ) {
       for ( k = 0; k <= usedN; k++ ) {

//         cc = k - BStartX + 2;
         k1 = k;
         ii = i1 + j1 + k1;
         lz = k *scale;
         //rz = (lz + 4)*scale;
         rz = lz + 4*scale;

         sum = 0.0;
         // Loop for grid point within the support of B-splines
         for  (a = lx+1; a < rx; a++) {
            al = (a - lx - 1)*4*size;

            for (b = ly+1; b < ry; b++) {
               bl = (b - ly -1)*4*(rz-lz-1);

               for (c = lz+1; c < rz; c++) {
                  cl = (c - lz - 1)*4;

                  abc = al+bl+cl;
                  if( i == 0 && j == 0 && k == 0 ) {
                       bspline->Phi_ijk_Partials_ImgGrid_3(a, b, c, i, j, k, partials);
                       partials_000[abc]   = partials[0];
                       partials_000[abc+1] = partials[1];
                       partials_000[abc+2] = partials[2];
                       partials_000[abc+3] = partials[3];
                  } else {
                     partials[0] = partials_000[abc];
                     partials[1] = partials_000[abc+1];
                     partials[2] = partials_000[abc+2];
                     partials[3] = partials_000[abc+3];
                  }

                  B_spline_Function_Gradient_Grid(coefs, gradient, a, b, c); // this guy is called 27 times for each point, 0912
                  length = InnerProduct(gradient, gradient);
                  if( length < epsilon ) length = epsilon;

//tem = 0.0;
                  
                  ///
                  ////tem = -partials[0] /length;
                  ////partials[0] = partials[1] + (gradient[0]*gradient[3]+gradient[1]*gradient[4]+gradient[2]*gradient[5]) *  tem;
                  ////partials[1] = partials[2] + (gradient[0]*gradient[4]+gradient[1]*gradient[6]+gradient[2]*gradient[7]) *  tem;
                  ////partials[2] = partials[3] + (gradient[0]*gradient[5]+gradient[1]*gradient[7]+gradient[2]*gradient[8]) *  tem;
                  

                  // without using gradient length
                  length = 1.0/sqrt(length);
                  partials[0] = partials[1]*length;
                  partials[1] = partials[2]*length;
                  partials[2] = partials[3]*length; 
                  
                  sum = sum + InnerProduct(gradient, partials);

               } 
            }
         }
         J234phi[ii] = J234phi[ii] +  taube*sum;
      }
   }
}
free(partials_000);
//printf("\nJ3=%f ", J234phi[ii]);
*/
}

// Cimpoute J_3, see book page 72, example 3.5.1
/*----------------------------------------------------------------------------*/
void  Reconstruction::ComputeJ3_MCF(float *coefs, float taube, float *J234phi)
{
int i, j, k, i1, j1, k1, ii, jj, kk, aa, bb, cc, abc, ilx, jly, klz;
int a, b, c, al, bl, cl, N2, lx, rx, ly, ry, lz, rz, scale, size, usedN1;
float sum, gradient[10],length,partials[20], tem, epsilon, *partials_000;
float Bx[10], Bx1[10], By[10], By1[10], Bz[10], Bz1[10], values[3], rscale;
float Bxy,Bxy1, Bx1y, Imatrix[9];

//printf("Begin Compute J3.\n");
epsilon = 0.000001;
epsilon = epsilon*epsilon; 
scale   = (int)Bscale;
rscale  = 1.0/Bscale;
usedN1 = usedN+1;
N2 = (usedN+1)*(usedN+1);
//size = (4*scale-1)*(4*scale-1);

//partials_000 = (float*)malloc(size*(4*scale-1)*4*sizeof(float)); 

//new methods.
 for (a = 0; a <= ImgNx; a++) {
    for (b = 0; b <= ImgNy; b++) {
       for (c = 0; c <= ImgNz; c++) {

          //B_spline_Function_Hessian_Grid(coefs, gradient, a, b, c); // this guy is called 27 times for each point, 0912
          B_spline_Function_Gradient_Grid(coefs, gradient, a, b, c); 

          length = InnerProduct(gradient, gradient);
          if( length < epsilon ) length = epsilon;
          length = 1.0/sqrt(length);

          gradient[0] = gradient[0]*length;
          gradient[1] = gradient[1]*length;
          gradient[2] = gradient[2]*length;

          lx = a/scale-3;
          rx = lx+3;
          if (lx < 0) lx = 0;
          if (rx > usedN) rx = usedN;
          if (lx > usedN ) lx = usedN;

          ly = b/scale-3;
          ry = ly+3;
          if (ly < 0) ly = 0;
          if (ry > usedN) ry = usedN;
          if (ly > usedN) ly = usedN;
          
          lz = c/scale-3;
          rz = lz+3;
          if (lz < 0) lz = 0;
          if (rz > usedN) rz = usedN;
          if (lz > usedN) lz = usedN;
 

          for ( i = lx; i<= rx; i++ ) {
             ilx = i - lx;
             bspline->Spline_N_Base_2(a-2-i, values);
             Bx[ilx]  = values[0];
             Bx1[ilx] = values[1]*rscale;
          }

          for ( j = ly; j <= ry; j++ ) {
             jly = j - ly;
             bspline->Spline_N_Base_2(b-2-j, values);
             By[jly]   = values[0];
             By1[jly]  = values[1]*rscale;
          }

          for ( k = lz; k <= rz; k++ ) {
             klz = k - lz;
             bspline->Spline_N_Base_2(c-2-k, values);
             Bz[klz]  = values[0];
             Bz1[klz] = values[1]*rscale;
          }

          for  (i = lx; i <= rx; i++) {
             ii = i*N2;
             ilx = i - lx; 
             for (j = ly; j <= ry; j++) {
                jj = j * usedN1;
                jly = j - ly;

                Bxy  = Bx[ilx]  * By[jly];
                Bx1y = Bx1[ilx] * By[jly];
                Bxy1 = Bx[ilx]  * By1[jly];


                for (k = lz; k <= rz; k++) {
                   klz = k - lz; 
                   kk = ii + jj + k;

                   partials[0] = Bx1y*Bz[klz];
                   partials[1] = Bxy1*Bz[klz];
                   partials[2] = Bxy *Bz1[klz];

                   /*
                   tem = Bxy*Bz[klz] /length;
                   partials[0] = partials[0] + (gradient[0]*gradient[3]+gradient[1]*gradient[4]+gradient[2]*gradient[5]) *  tem;
                   partials[1] = partials[1] + (gradient[0]*gradient[4]+gradient[1]*gradient[6]+gradient[2]*gradient[7]) *  tem;
                   partials[2] = partials[2] + (gradient[0]*gradient[5]+gradient[1]*gradient[7]+gradient[2]*gradient[8]) *  tem;
                   */

                   
                   //partials[0] = partials[0]*length;
                   //partials[1] = partials[1]*length;
                   //partials[2] = partials[2]*length;
                   

                   sum = InnerProduct(gradient, partials);
                   J234phi[kk] = J234phi[kk] +  taube*sum;

                }
             }
          }
       }
    }
 }


}


// Cimpoute J_3, see book page 68, (3.3.7)
/*----------------------------------------------------------------------------*/
void  Reconstruction::ComputeJ3_WMF(float *coefs, float taube, float *J234phi)
{
int i, j, k, i1, j1, k1, ii, jj, kk, aa, bb, cc, abc, ilx, jly, klz;
int a, b, c, al, bl, cl, N2, lx, rx, ly, ry, lz, rz, scale, size, usedN1;
float sum, gradient[20],length,partials[20], tem, epsilon, *partials_000, 
      grad, grad2, normal[3], gradient2[3],rgrad, rgrad2, rgrad4,H, H2;
float Bx[10], Bx1[10], Bx2[10], By[10], By1[10], By2[10], Bz[10], Bz1[10], Bz2[10], values[3], rscale, rscale2;
float Bxy,Bxy1, Bx1y, Bx1y1, Bx2y, Bxy2;  //omega = H * norm(f).
float div_Pgphi1, div_Pgphi2, Hf[9], HGf[3];
float Hnormal[3], H2normal[3], Hfn[3], nHfn, nHfn_normal[3], H_phi[9], nH_phi[3];
float ww, c1, c2, c3, c4, c5, c6, c7, c8, c9;

//printf("Begin Compute J3.\n");
epsilon = 0.000001;
//epsilon = epsilon*epsilon;
scale   = (int)Bscale;
rscale  = 1.0/Bscale;
rscale2 = rscale*rscale;  

usedN1 = usedN+1;
N2 = (usedN+1)*(usedN+1);
//size = (4*scale-1)*(4*scale-1);

//new methods.
 for (a = 0; a <= ImgNx; a++) {
    for (b = 0; b <= ImgNy; b++) {
       for (c = 0; c <= ImgNz; c++) {

          B_spline_Function_Hessian_Grid(coefs, gradient, a, b, c); // Xu 0924 
          //B_spline_Function_Gradient_Grid_All(coefs, gradient, a, b, c); 

          gradient2[0] = gradient[1] * gradient[1];
          gradient2[1] = gradient[2] * gradient[2];
          gradient2[2] = gradient[3] * gradient[3];

          grad2 = gradient2[0]+ gradient2[1] + gradient2[2];
          if (grad2 < epsilon) grad2 = epsilon;

          grad = sqrt(grad2);
          rgrad = 1.0/grad;
          rgrad2 = 1.0/grad2;

          //Compute mean curvature H = - 0.5 * div(normal);
          H =   gradient[1]*(gradient[4]*gradient[1] + gradient[5]*gradient[2] + gradient[6]*gradient[3])
              + gradient[2]*(gradient[5]*gradient[1] + gradient[7]*gradient[2] + gradient[8]*gradient[3])
              + gradient[3]*(gradient[6]*gradient[1] + gradient[8]*gradient[2] + gradient[9]*gradient[3]);

          H = H*rgrad2 - gradient[4] - gradient[7] - gradient[9];
          H = 0.5 * H * rgrad;
          H2 = H * H;
          normal[0] = gradient[1] * rgrad ;
          normal[1] = gradient[2] * rgrad ;
          normal[2] = gradient[3] * rgrad ;
        
          // Hessian of f 
          Hf[0] = gradient[4];
          Hf[1] = gradient[5];
          Hf[2] = gradient[6];
          Hf[3] = gradient[5];
          Hf[4] = gradient[7];
          Hf[5] = gradient[8];
          Hf[6] = gradient[6];
          Hf[7] = gradient[8];
          Hf[8] = gradient[9];
 
          // H*n
          MatrixMultiply(Hf, 3,3, normal, 3, 1, Hfn);

          // Hn * 2H / |\nabla f\| 
          ww = rgrad*H*2;
          Hfn[0] = Hfn[0]*ww;
          Hfn[1] = Hfn[1]*ww;
          Hfn[2] = Hfn[2]*ww;

          // n^T * H * n
          nHfn  = InnerProduct(normal, Hfn);

          Hnormal[0] = normal[0] * H;
          Hnormal[1] = normal[1] * H;
          Hnormal[2] = normal[2] * H;

          H2normal[0] = normal[0] * H2;
          H2normal[1] = normal[1] * H2;
          H2normal[2] = normal[2] * H2;


          nHfn_normal[0] = normal[0] * nHfn;
          nHfn_normal[1] = normal[1] * nHfn;
          nHfn_normal[2] = normal[2] * nHfn;
          
          // Xu add these c's for speedup 0925
          c1 =  - H2normal[0] + Hfn[0] - nHfn_normal[0];
          c2 =  - H2normal[1] + Hfn[1] - nHfn_normal[1];
          c3 =  - H2normal[2] + Hfn[2] - nHfn_normal[2];

          c4 =   Hnormal[0]*normal[0];
          c5 = 2*Hnormal[0]*normal[1];
          c6 = 2*Hnormal[0]*normal[2];
          c7 =   Hnormal[1]*normal[1];
          c8 = 2*Hnormal[1]*normal[2];
          c9 =   Hnormal[2]*normal[2];

          lx = a/scale-3;
          rx = lx+3;
          if (lx < 0) lx = 0;
          if (rx > usedN) rx = usedN;
          if (lx > usedN ) lx = usedN;

          ly = b/scale-3;
          ry = ly+3;
          if (ly < 0) ly = 0;
          if (ry > usedN) ry = usedN;
          if (ly > usedN) ly = usedN;
          
          lz = c/scale-3;
          rz = lz+3;
          if (lz < 0) lz = 0;
          if (rz > usedN) rz = usedN;
          if (lz > usedN) lz = usedN;
 

          for ( i = lx; i<= rx; i++ ) {
             ilx = i - lx;
             bspline->Spline_N_Base_2(a-2-i, values);
             Bx[ilx]  = values[0];
             Bx1[ilx] = values[1]*rscale;
             Bx2[ilx] = values[2]*rscale2;
          }

          for ( j = ly; j <= ry; j++ ) {
             jly = j - ly;
             bspline->Spline_N_Base_2(b-2-j, values);
             By[jly]   = values[0];
             By1[jly]  = values[1]*rscale;
             By2[jly]  = values[2]*rscale2;
          }

          for ( k = lz; k <= rz; k++ ) {
             klz = k - lz;
             bspline->Spline_N_Base_2(c-2-k, values);
             Bz[klz]  = values[0];
             Bz1[klz] = values[1]*rscale;
             Bz2[klz] = values[2]*rscale2;
          }

          for  (i = lx; i <= rx; i++) {
             ii = i*N2;
             ilx = i - lx; 
             for (j = ly; j <= ry; j++) {
                jj = j * usedN1;
                jly = j - ly;

                Bxy  = Bx[ilx]  * By[jly];
                Bx1y = Bx1[ilx] * By[jly];
                Bxy1 = Bx[ilx]  * By1[jly];
                Bx1y1= Bx1[ilx] * By1[jly];

                Bx2y = Bx2[ilx] * By[jly];
                Bxy2 = Bx[ilx]  * By2[jly];

                for (k = lz; k <= rz; k++) {
                   klz = k - lz; 
                   kk = ii + jj + k;

                   // grad(\phi)
                   partials[0] = Bx1y*Bz[klz];
                   partials[1] = Bxy1*Bz[klz];
                   partials[2] = Bxy *Bz1[klz];

                   // Hessian(\phi)
                   H_phi[0]= Bx2y *Bz[klz];
                   H_phi[1]= Bx1y1*Bz[klz];
                   H_phi[2]= Bx1y *Bz1[klz];

                   //H_phi[3]= H_phi[1]; 
                   H_phi[4]= Bxy2 *Bz[klz];
                   H_phi[5]= Bxy1  *Bz1[klz];

                   //H_phi[6]= H_phi[2];
                   //H_phi[7]= H_phi[5];
                   H_phi[8]= Bxy  *Bz2[klz];


                   /* 
                   //       H^2 * n^T * grad(\phi)               - H \Delta phi  
                   sum = -InnerProduct(H2normal, partials) - H * (H_phi[0] +  H_phi[4] + H_phi[8]);

                   //       2*H* n^T * G(f) * grad(\phi)/|grad(f)| 
                   sum += InnerProduct(Hfn, partials);

                   //        H * n^T *G(\phi) 
                   MatrixMultiply(Hnormal,1,3,H_phi,3,3, nH_phi);

                   //         n^T * G(\phi) * n               -2n^TG(f)*n*n^T * grad(\phi)
                   sum += InnerProduct(nH_phi, normal) - InnerProduct(nHfn_normal, partials); 
                  */ 


                   // Xu replaced above 4 sentences wit the following for speedup, 0925
                   sum = - H * (H_phi[0] +  H_phi[4] + H_phi[8]);

                   sum = sum + c1*partials[0] + c2*partials[1] + c3*partials[2] 
                             + c4*H_phi[0]    + c5*H_phi[1]    + c6*H_phi[2]
                             + c7*H_phi[4]    + c8*H_phi[5]    + c9*H_phi[8]; 
                       

                   J234phi[kk] += taube * sum;

                    //printf("\nsum=%f inner=%f H=%f  coef=%f %f %f %f %f %f %f %f %f ", sum, InnerProduct(partials, HGf), 
                    //      H, coef_phix, coef_phiy, coef_phiz, coef_phixx, coef_phixy, coef_phixz, coef_phiyy, coef_phiyz, coef_phizz);          
                }
             }
          }
       }
    }
 }

}

//------------------------------------------------------------------------------
void Reconstruction::B_spline_Function_Gradient_Grid(float *coefs, float *gradient, int ix, int iy, int iz)
{
int i,j,k, l, ii, scale, lx, rx, ly, ry, lz, rz, N2;
float partials[20], t[10];
int i1, j1, usedN1;

scale   = (int)Bscale;
usedN1 = usedN+1;
N2 = usedN1*usedN1;

rx = ix/scale;
lx = rx-3;

ry = iy/scale;
ly = ry-3;

rz = iz/scale;
lz = rz-3;


//for (l = 0; l < 10; l++) {
for (l = 0; l < 3; l++) {
   gradient[l] = 0.0;
   //t[l] = 0.0;
}

for ( i = lx; i<= rx; i++ ) {
   if( i < 0 || i > usedN ) continue;

   i1 = i*N2;
   for ( j = ly; j <= ry; j++ ) {
      if ( j < 0 || j > usedN ) continue;

      j1 = j*usedN1; 
      for ( k = lz; k <= rz; k++ ) {
         if ( k < 0 || k > usedN ) continue;
         //printf("\ni j k =%d %d %d ", i , j, k );
         //ii = i*N2 + j*(usedN+1) + k;
         ii = i1 + j1 + k;

         //bspline->Phi_ijk_Partials_ImgGrid(ix, iy, iz, i, j, k, partials);
         bspline->Phi_ijk_Partials_ImgGrid_3(ix, iy, iz, i, j, k, partials);
//coefs[ii] = 1.0;
//coefs[ii] = i+j+k;


         //for (l = 0; l < 9; l++) {
         for (l = 0; l < 3; l++) {
            gradient[l] = gradient[l] + coefs[ii]*partials[l+1];
            //t[l] = t[l]  + coefs[ii]*partials[l];
         }
         
      }
   }
}

 //  printf("ix ,iy iz = %d, %d, %d, grad = %f,  %f,%f,%f,  %f, %f, %f, %f,%f\n", ix, iy, iz, t[0],t[1],t[2],t[3],t[4],t[5],t[6],t[7],t[8]);
}

//0921. added.
void Reconstruction::B_spline_Function_Gradient_Grid_All(float *coefs, float *gradient, int ix, int iy, int iz)
{
int i,j,k, l, ii, scale, lx, rx, ly, ry, lz, rz, N2;
float partials[20], t[10];
int i1, j1, usedN1;

scale   = (int)Bscale;
usedN1 = usedN+1;
N2 = usedN1*usedN1;

rx = ix/scale;
lx = rx-3;

ry = iy/scale;
ly = ry-3;

rz = iz/scale;
lz = rz-3;


//for (l = 0; l < 10; l++) {
for (l = 0; l < 20; l++) {
   gradient[l] = 0.0;
   //t[l] = 0.0;
}
for ( i = lx; i<= rx; i++ ) {
   if( i < 0 || i > usedN ) continue;

   i1 = i*N2;
   for ( j = ly; j <= ry; j++ ) {
      if ( j < 0 || j > usedN ) continue;

      j1 = j*usedN1;
      for ( k = lz; k <= rz; k++ ) {
         if ( k < 0 || k > usedN ) continue;
         ii = i1 + j1 + k;

         bspline->Phi_ijk_Partials_ImgGrid(ix, iy, iz, i, j, k, partials);
         for (l = 0; l < 20; l++) {
            gradient[l] = gradient[l] + coefs[ii]*partials[l];
         }

      }
   }
}

}




void Reconstruction::B_spline_Function_Hessian_Grid(float *coefs, float *gradient, int ix, int iy, int iz)
{
int i,j,k, l, ii, scale, lx, rx, ly, ry, lz, rz, N2;
float partials[10];
int i1, j1, usedN1;

scale   = (int)Bscale;
usedN1 = usedN+1;
N2 = usedN1*usedN1;

rx = ix/scale;
lx = rx-3;

ry = iy/scale;
ly = ry-3;

rz = iz/scale;
lz = rz-3;


for (l = 0; l < 10; l++) {
   gradient[l] = 0.0;
}

for ( i = lx; i<= rx; i++ ) {
   if( i < 0 || i > usedN ) continue;

   i1 = i*N2;
   for ( j = ly; j <= ry; j++ ) {
      if ( j < 0 || j > usedN ) continue;

      j1 = j*usedN1; 
      for ( k = lz; k <= rz; k++ ) {
         if ( k < 0 || k > usedN ) continue;
         ii = i1 + j1 + k;

         bspline->Phi_ijk_Partials_ImgGrid_9(ix, iy, iz, i, j, k, partials);
         for (l = 0; l < 10; l++) {
            //gradient[l] = gradient[l] + (i+j+k)*partials[l]; // Xu Test plane
            //gradient[l] = gradient[l] + (i*i+j*j+k*k)*partials[l]; // Xu Test sphere

            gradient[l] = gradient[l] + coefs[ii]*partials[l];
         }
      }
   }
}

}


// ----------------------------------------------------------------------------
void CubeCoefficients()
{
int ii, i, j, k, x, y, z;
float df[3];

ii = 0;
for ( x = 0; x <= 2; x++ ) {
   df[0] = x * 0.5;

   for ( y = 0; y <= 2; y++ ) {
      df[1] = y * 0.5;

      for ( z = 0; z <= 2; z++, ii++) {
          df[2] = z * 0.5;
          
          Cube_Coeff[8*ii + 0] =  (1-df[2])*(1-df[1])*(1-df[0]);
          Cube_Coeff[8*ii + 1] =  df[2]    *(1-df[1])*(1-df[0]);
          Cube_Coeff[8*ii + 2] =  (1-df[2])*df[1]    *(1-df[0]);
          Cube_Coeff[8*ii + 3] =  df[2]    *df[1]    *(1-df[0]); 
          Cube_Coeff[8*ii + 4] =  (1-df[2])*(1-df[1])*df[0];
          Cube_Coeff[8*ii + 5] =  df[2]    *(1-df[1])*df[0];
          Cube_Coeff[8*ii + 6] =  (1-df[2])*df[1]    *df[0];
          Cube_Coeff[8*ii + 7] =  df[2]    *df[1]    *df[0];
      }
   }
}

}

void Reconstruction::readFiles(const char* filename, const char* path, int dimN)   //0928.
{
  char prjfile[50];
  //char* p;
    int i, j, numb, id, imgsize, imgnx;
  FILE* fp;
  EulerAngles  *home, *eulers;
  int          dim[3];

  if(dimN%2 == 0 ) {                    //0928.
    imgnx = ImgNx;
    imgsize = ImgNx * ImgNx;
  }
  else {
    imgnx = ImgNx + 1;
    imgsize = imgnx * imgnx;
  }


  prjimg = (PrjImg *)malloc(sizeof(PrjImg));

  prjimg->data = (float *)malloc(nv * imgsize * sizeof(float));

  float *data=NULL;

  QDir dataDir;
  dataDir.setCurrent(path);
  printf("\nnewcurrentpath=%s ", dataDir.current().path().ascii());
  printf("\nfilename =%s ", filename);

  if( (fp = fopen(filename,"r")) == NULL)
        {printf("\nerror read files. ");return ;}

  id = 0;

  eulers =  new EulerAngles();
  volfilename = "newvolume.rawiv";

  home = eulers;
  while(id < nv)
    {

     fscanf(fp, "%s %d", &prjfile, &numb);
     readSpiderFile(prjfile,'I',dim, &data, eulers);

     for ( i = 0; i < imgnx; i++ )
        for ( j = 0; j < imgnx; j++ )
             prjimg->data[id*imgsize+i*imgnx+j ] = data[j*imgnx+i];


      printf("\nrot =%f tilt=%f psi=%f ", eulers->rot, eulers->tilt, eulers->psi);
      if(id < nv - 1)
        {
           eulers->next =  new EulerAngles();
           eulers = eulers->next;
         }
       else
         eulers->next = NULL;
      id++;
 
    free(data); data = NULL; 
    }

  EulerMatrice(Rmat,home);

  fclose(fp);

  if ( dimN%2 != 0 )  {                           //0928
     free(gd);
     gd = prjimg;gd->data = prjimg->data;

  }

  delete home;  home = NULL;

}

void Reconstruction::imageInterpolation()
{


    //      PrjImg *imgp = (PrjImg*)malloc(sizeof(PrjImg));
    //      imgp->data = (float*)malloc(nv*img2dsize*sizeof(float));

          int v, i, j, imgsize, size;
          imgsize = ImgNx * ImgNx;

          float *data  = (float*)malloc(imgsize * sizeof(float));
          float *data1 = (float*)malloc(img2dsize * sizeof(float));


          for ( v = 0; v <nv; v++ )
                {
                  for (  i = 0; i < imgsize; i++ )
                        data[i] = prjimg->data[v*imgsize+i];

                  size = v*img2dsize;
                  Imageinterpo_BiLinear(data, ImgNx, ImgNx, data1);
                  for (  i = 0; i < img2dsize; i++ )
                        gd->data[size+i] = data1[i];
                        //imgp->data[size+i] = data1[i];
               }

          free(data1);  data1 = NULL;
          free(prjimg->data); free(prjimg); prjimg = NULL;
          //prjimg = imgp;
          //gd = imgp;
          free(data); data = NULL;

}

void Reconstruction::SetNewNv(int NewNv)
{
 newnv = NewNv;
}

void Reconstruction::SetBandWidth(int BandWidth)
{
 bandwidth = BandWidth;
}

void Reconstruction::SetFlow(int Flow)
{
 flows = Flow;
}


void Reconstruction::setOrders(const int tolnv)
{


     //orderSet = new OrderSet();
     int     i,j;

     //orderSet->nv = SubSetN;

     //orderSet->id = (int *)malloc(tolnv*sizeof(int));


     int     *chosen, minprod_id, listi,listN, NN;
     float   *product, minprod, z1[3], z2[3], z3[3];


     NN = 2;
     chosen  = (int *)malloc(tolnv * sizeof(int));
     product = (float *)malloc(tolnv * sizeof(float));

     for (j = 0; j < tolnv; j++ )
        {
          chosen[j]       = 0;
          OrthoOrdered[j] = 0;
          product[j]      = 0.0;
        }

     i = 0;  // Pick first projection as the first one.
     chosen[i] = 1;
      OrthoOrdered[0] = i;

      std::cerr<<"sorting projections as orthogonal orders.\n";;

      for( i = 1; i < tolnv; i++ )
         {
           minprod = MAXFLOAT;
           listi =OrthoOrdered[i-1];
           z1[0] = Rmat[listi*9 + 6 ];
           z1[1] = Rmat[listi*9 + 7 ];
           z1[2] = Rmat[listi*9 + 8 ];
           //printf("\nz1=%f z2=%f z3=%f ", z1[0], z1[1], z1[2]);

           if( i > NN )
              {
                listN = OrthoOrdered[i-NN-1];
                z3[0] = Rmat[listN*9 + 6 ];
                z3[1] = Rmat[listN*9 + 7 ];
                z3[2] = Rmat[listN*9 + 8 ];
              }
           for ( j = 0; j < tolnv; j++ )
               if(chosen[j] == 0 )
                   {
                     z2[0] = Rmat[j*9 + 6 ];
                     z2[1] = Rmat[j*9 + 7 ];
                     z2[2] = Rmat[j*9 + 8 ];

                     product[j] += fabs(InnerProduct(z1, z2));
                     if(i > NN )
                        product[j] -= fabs(InnerProduct(z3, z2));
                    if( product[j] < minprod )
                         {
                           minprod = product[j];
                           minprod_id = j;
                         }
                   }

            OrthoOrdered[i] = minprod_id;
            chosen[minprod_id] = 1;

         }

     printf("\nOrtho Order \n");
     for ( i = 0; i < tolnv; i++ ) printf("%d  ", OrthoOrdered[i]);
         printf("\nfinished.");
      free(chosen);        chosen = NULL;
      free(product);      product = NULL;
}


void Reconstruction::getSubSetFromOrthoSets(int *subset)
{

 int i, newnv_1;
 newnv_1 = newnv;

 if(CurrentId == nv ) CurrentId = 0;
 else if(CurrentId + newnv > nv ) newnv_1 = nv - CurrentId;
for ( i = 0; i < newnv_1; i++ )
   {subset[i] = OrthoOrdered[CurrentId+i];
    printf("\nsubSet=%d ", subset[i]);}
           
 
 CurrentId += newnv_1;


}

void Reconstruction::setOrderManner(int ordermanner)
{
OrderManner = ordermanner;
}
