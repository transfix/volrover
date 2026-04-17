/*
  Copyright 2006 The University of Texas at Austin

        Authors: Sangmin Park <smpark@ices.utexas.edu>
        Advisor: Chandrajit Bajaj <bajaj@cs.utexas.edu>

  This file is part of PEDetection.

  PEDetection is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License version 2.1 as published by the Free Software Foundation.

  PEDetection is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/

#include <stdlib.h>
#include <stdio.h>
#include <iostream.h>
#include <math.h>
#include <float.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#ifndef _WIN32
#include <unistd.h>
#endif

#include <PEDetection/Skeleton.h>
#include <PEDetection/CompileOptions.h>

#define	MIN(X, Y) { (((X)<(Y))? (return (X)): (return Y)); }
#define	SWAP(A, B, Temp) { (Temp)=(A); (A)=(B); (B)=(Temp); }

//----------------------------------------------------------------------------
// cSkeleton Class Member Functions
//----------------------------------------------------------------------------

template <class _DataType>
cSkeleton<_DataType>::cSkeleton()
{
	InitVolume_muc = NULL;
	Distance_mi = NULL;
	VoxelFlags_muc = NULL;
	CCVolume_muc = NULL;
	GVFDistance_mf = NULL;
	DistanceFromSkeletons_mi = NULL;
	Skeletons_muc = NULL;
}


// destructor
template <class _DataType>
cSkeleton<_DataType>::~cSkeleton()
{
	delete [] InitVolume_muc;
	delete [] Distance_mi;
	delete [] VoxelFlags_muc;
	delete [] CCVolume_muc;
	delete [] GVFDistance_mf;
	delete [] DistanceFromSkeletons_mi;
	delete [] Skeletons_muc;
}


template <class _DataType>
void cSkeleton<_DataType>::setData(_DataType *Data, float Minf, float Maxf)
{
	MinData_mf = Minf;
	MaxData_mf = Maxf;
	Data_mT = Data;
}


template<class _DataType>
void cSkeleton<_DataType>::setWHD(int W, int H, int D)
{
	Width_mi = W;
	Height_mi = H;
	Depth_mi = D;
	
	WtimesH_mi = W*H;
	WHD_mi = W*H*D;
	
}

template<class _DataType>
void cSkeleton<_DataType>::setProbabilityHistogram(float *Prob, int NumMaterial, int *Histo, float HistoF)
{
	MaterialProb_mf = Prob;
	Histogram_mi = Histo;
	HistogramFactorI_mf = HistoF;
	HistogramFactorG_mf = HistoF;
	NumMaterial_mi = NumMaterial;
}

template <class _DataType>
void cSkeleton<_DataType>::setGradient(float *Grad, float Minf, float Maxf)
{
	GradientMag_mf = Grad;
	MinGrad_mf = Minf;
	MaxGrad_mf = Maxf;
}

template<class _DataType>
void cSkeleton<_DataType>::setSecondDerivative(float *SecondD, float Min, float Max)
{ 
	SecondDerivative_mf = SecondD; 
	MinSecond_mf = Min;
	MaxSecond_mf = Max;
}

template<class _DataType>
void cSkeleton<_DataType>::setGradientVectors(float *GVec)
{ 
	GradientVec_mf = GVec;
}

template<class _DataType>
unsigned char *cSkeleton<_DataType>::getFlaggedVoxelVolume()
{ 
	return VoxelFlags_muc;
}

template<class _DataType>
int *cSkeleton<_DataType>::getDistanceVolume()
{
	return Distance_mi;
}


template<class _DataType>
void cSkeleton<_DataType>::Skeletonize(char *OutFileName, _DataType MatMin, _DataType MatMax)
{

	BinarySegment(MatMin, MatMax); // Using only Min & Max intensity values
//	BinarySegment2(MatMin, MatMax); // Removing the zero 2nd D values
//	SaveInitVolume(MatMin, MatMax);
	

	ComputeDistance(); // Compute the distance from bondary field
//	SaveDistanceVolume();
	

	ComputeGVF(); // Allocate memory to "GVFDistance_mf" and Compute GVF
	
	FlagNonUniformGradient(); // Marking nonuniform voxels
	SaveVolume(VoxelFlags_muc, (float)0.0, (float)255.0, "Flagged");
	printf ("The flagged voxel volume is saved\n"); fflush(stdout);

	ConnectingFlaggedVoxels();	// Original
//	ConnectingFlaggedVoxels2(); // Compute Ave. Radius
	SaveVolume(VoxelFlags_muc, (float)0.0, (float)255.0, "Connected");
	printf ("The Connected voxel volume is saved\n"); fflush(stdout);
	
/*	
	int		MaxCCLoc_Ret;
	ConnectedComponents("CC", MaxCCLoc_Ret);
	SaveVolume(CCVolume_muc, (float)0.0, (float)255.0, "CCVolume");
	printf ("The Connected Component volume is saved\n"); fflush(stdout);


	int		MaxCCLoc = MaxCCLoc_Ret, RootLoc_Ret, EndLoc_Ret;
	FindingRootAndEndLoc(MaxCCLoc, RootLoc_Ret, EndLoc_Ret);

	int		Threshold_Distance = 170;
	ComputeSkeletons(RootLoc_Ret, EndLoc_Ret, Threshold_Distance);
	SaveVolume(Skeletons_muc, (float)0.0, (float)255.0, "Skeletons");
	printf ("The skeleton volume dataset is saved\n"); fflush(stdout);
*/


}


// Using only min & max intensity values
template<class _DataType>
void cSkeleton<_DataType>::BinarySegment(_DataType MatMin, _DataType MatMax)
{
	int		i, NumSegmentedVoxels;

	
	printf ("Binary Segmentation ... \n");

	delete [] Distance_mi;
	Distance_mi = new int[WHD_mi];
	delete [] VoxelFlags_muc;
	VoxelFlags_muc = new unsigned char[WHD_mi];
	
	for (i=0; i<WHD_mi; i++) {
		Distance_mi[i] = 0;
		VoxelFlags_muc[i] = (unsigned char)FLAG_EMPTY;
	}

	NumSegmentedVoxels = 0;
	for (i=0; i<WHD_mi; i++) {
		if (MatMin <= Data_mT[i] && Data_mT[i] <= MatMax) {
			Distance_mi[i] = (int)255;
			VoxelFlags_muc[i] = FLAG_SEGMENTED;
			NumSegmentedVoxels++;
		}
		else Distance_mi[i] = (int)0;
	}
	
	// It is for the Index() function
	Distance_mi[0] = (int)0;
	VoxelFlags_muc[0] = (unsigned char)FLAG_EMPTY;
	NumSegmentedVoxels_mi = NumSegmentedVoxels;



	delete [] InitVolume_muc;
	InitVolume_muc = new unsigned char [WHD_mi];
	
	NumSegmentedVoxels = 0;
	for (i=0; i<WHD_mi; i++) {
		if (Distance_mi[i]==255) InitVolume_muc[i] = (unsigned char)255;
		else InitVolume_muc[i] = (unsigned char)0;
	}



	
	printf ("Using only intensity values\n");
	printf ("Num. Segmented Voxels = %d\n", NumSegmentedVoxels_mi);
	printf ("Total Num. Voxels = %d\n", WHD_mi);
	printf ("Num. Segmented / Total Num. Voxels = %f %%\n", (double)NumSegmentedVoxels_mi/WHD_mi*100.0);
	printf ("\n");
	fflush(stdout);
	
}


// Using second derivative
// Removing segmented voxels that have zero-crossing second derivatives
// Using FindZeroCrossingLocation() to get the zero 2nd derivative locations
template<class _DataType>
void cSkeleton<_DataType>::BinarySegment2(_DataType MatMin, _DataType MatMax)
{
	int		i, j, k, l, m, n, loc[8], CubeXYZ_i[8][3], NumSegmentedVoxels;
	int		NumPosSigns, Idx;
	double	GradVec_d[3], DataLoc_d[3], Length_d;
	double	AveGradM_d, Step_d, GradM_d[15];
	
	
	printf ("Binary Segmentation 2 ... \n"); fflush(stdout);

	delete [] Distance_mi;
	Distance_mi = new int[WHD_mi];
	delete [] VoxelFlags_muc;
	VoxelFlags_muc = new unsigned char[WHD_mi];
	

	for (i=0; i<WHD_mi; i++) {
		if (MatMin <= Data_mT[i] && Data_mT[i] <= MatMax) {
			Distance_mi[i] = (int)255;
			VoxelFlags_muc[i] = FLAG_SEGMENTED;
		}
		else {
			Distance_mi[i] = (int)0;
			VoxelFlags_muc[i] = (unsigned char)FLAG_EMPTY;
		}
	}
	// It is for the Index() function
	Distance_mi[0] = (int)0;
	VoxelFlags_muc[0] = (unsigned char)FLAG_EMPTY;
	

	for (k=1; k<Depth_mi-1; k++) {
		for (j=1; j<Height_mi-1; j++) {
			for (i=1; i<Width_mi-1; i++) {
				
				loc[0] = Index(i, j, k);
				if (Distance_mi[loc[0]]==0) continue;
				
				Idx = 0;
				NumPosSigns = 0;
				AveGradM_d = 0.0;
				for (n=k; n<=k+1; n++) {
					for (m=j; m<=j+1; m++) {
						for (l=i; l<=i+1; l++) {
							loc[Idx] = Index(l, m, n);
							if (SecondDerivative_mf[loc[Idx]]>=0) NumPosSigns++;
							CubeXYZ_i[Idx][0] = l;
							CubeXYZ_i[Idx][1] = m;
							CubeXYZ_i[Idx][2] = n;
							AveGradM_d += (double)GradientMag_mf[loc[Idx]];
							Idx++;
						}
					}
				}
				AveGradM_d /= 8.0;
				
				if (NumPosSigns<8 && NumPosSigns>0 && AveGradM_d>5.0) {

					for (n=0; n<8; n++) {
						
						for (l=0; l<3; l++) GradVec_d[l] = (double)GradientVec_mf[loc[n]*3 + l];
						Length_d = 0.0;
						for (l=0; l<3; l++) Length_d += GradVec_d[l]*GradVec_d[l];
						Length_d = sqrt (Length_d);
						for (l=0; l<3; l++) GradVec_d[l] /= Length_d;
						
						Idx = 0;
						Step_d = -1.5;
						for (l=0; l<3; l++) DataLoc_d[l] = (double)CubeXYZ_i[n][l] + GradVec_d[l]*Step_d;
						GradM_d[Idx] = GradientInterpolation(DataLoc_d);
						Idx++;
						
						Step_d = -1.5 + 0.25;
						for (l=0; l<3; l++) DataLoc_d[l] = (double)CubeXYZ_i[n][l] + GradVec_d[l]*Step_d;
						GradM_d[Idx] = GradientInterpolation(DataLoc_d);
						Idx++;
						
						for (Step_d = -1.5 + 0.25*2.0; Step_d<=1.5; Step_d+=0.25) {
							for (l=0; l<3; l++) DataLoc_d[l] = (double)CubeXYZ_i[n][l] + GradVec_d[l]*Step_d;
							GradM_d[Idx] = GradientInterpolation(DataLoc_d);

							
							if (GradM_d[Idx-2] < GradM_d[Idx-1] && GradM_d[Idx-1] > GradM_d[Idx]) {
								Distance_mi[loc[n]] = (int)0;
								VoxelFlags_muc[loc[n]] = (unsigned char)FLAG_EMPTY;
								break;
							}
							
							
							Idx++;
						}
						/*
						printf ("\n");
						printf ("(%3d %3d %3d): ", i, j, k);
						printf ("Vec = (%6.3f %6.3f %6.3f), ", GradVec_d[0], GradVec_d[1], GradVec_d[2]);
						printf ("Second D = ");
						for (l=0; l<Idx; l++) {
							printf ("%8.3f ", SecondD_d[l]);
						}
						printf ("\n");
						*/
						
					}
					
				}

			}
		}
	}


	delete [] InitVolume_muc;
	InitVolume_muc = new unsigned char [WHD_mi];
	

	NumSegmentedVoxels = 0;
	for (i=0; i<WHD_mi; i++) {

		if (Distance_mi[i]==255) {
			if (VoxelFlags_muc[i] != FLAG_SEGMENTED) {
				printf ("Error!!!, VoxelFlags_muc[i] should be FLAG_SEGMENTED\n");
				exit(1);
			}
			NumSegmentedVoxels++;
			InitVolume_muc[i] = (unsigned char)255;
		}
		else {
			if (VoxelFlags_muc[i] != FLAG_EMPTY) {
				printf ("Error!!!, VoxelFlags_muc[i] should be FLAG_EMPTY\n");
				exit(1);
			}
			InitVolume_muc[i] = (unsigned char)0;
		}

	}
	NumSegmentedVoxels_mi = NumSegmentedVoxels;
	
	printf ("The voxels of zero second derivative are removed\n");
	printf ("Num. Segmented Voxels = %d\n", NumSegmentedVoxels_mi);
	printf ("Total Num. Voxels = %d\n", WHD_mi);
	printf ("Num. Segmented / Total Num. Voxels = %f %%\n", (double)NumSegmentedVoxels_mi/WHD_mi*100.0);
	printf ("\n");
	fflush(stdout);
	
}


template<class _DataType>
void cSkeleton<_DataType>::ComputeDistance()
{
	int		i, j, k, n, df_i, db_i,  d_i, Tempi;
	int		MaxRes, *Buffer_i;
	

	printf ("Computing Distance ... \n"); fflush(stdout);
	if (Distance_mi == NULL) { printf ("Distance_mi is NULL\n"); exit(1); }

	MaxRes = (Width_mi>Height_mi)? Width_mi : Height_mi;
	MaxRes = (MaxRes>Depth_mi)? MaxRes : Depth_mi;
	Buffer_i = new int [MaxRes];
	
	
	// Step 1: X-axis
	// Forward scan
	for (k=0; k<Depth_mi; k++) {
		for (j=0; j<Height_mi; j++) {
			df_i = Width_mi-1;
			for (i=0; i<Width_mi; i++) {
				if (Distance_mi[Index(i, j, k)]>0) df_i++;
				else df_i = 0;
				Distance_mi[Index(i, j, k)] = df_i*df_i;
			}
		}
	}

	// Backward scan
	for (k=0; k<Depth_mi; k++) {
		for (j=0; j<Height_mi; j++) {
			db_i = Width_mi-1;
			for (i=Width_mi-1; i>=0; i--) {
				if (Distance_mi[Index(i, j, k)]>0) db_i++;
				else db_i = 0;
				Tempi = Distance_mi[Index(i, j, k)];
				Distance_mi[Index(i, j, k)] = (Tempi < db_i*db_i)? Tempi : db_i*db_i;
			}
		}
	}


	// Step 2: Y-axis
	int		w_i, rStart, rMax, rEnd;

	for (k=0; k<Depth_mi; k++) {
		for (i=0; i<Width_mi; i++) {
			for (j=0; j<Height_mi; j++) {
				Buffer_i[j] = Distance_mi[Index(i, j, k)];
			}
			
			for (j=0; j<Height_mi; j++) {
				d_i = Buffer_i[j];
				if (d_i>0) {
					rMax = (int)(sqrt((double)d_i)) + 1;
					rStart = (rMax<(j-1))? rMax : (j-1);
					rEnd = (rMax<(Height_mi-1 - j))? rMax : (Height_mi-1 - j);
					for (n=-rStart; n<=rEnd; n++) {
						w_i = Buffer_i[j+n] + n*n;
						if (w_i<d_i) d_i = w_i;
					}
				}
				
				Distance_mi[Index(i, j, k)] = d_i;
			}
		}
	}

	// Step 3: Z-axis

	for (j=0; j<Height_mi; j++) {
		for (i=0; i<Width_mi; i++) {
			for (k=0; k<Depth_mi; k++) {
				Buffer_i[k] = Distance_mi[Index(i, j, k)];
			}
			
			for (k=0; k<Depth_mi; k++) {
				d_i = Buffer_i[k];
				if (d_i>0) {
					rMax = (int)(sqrt((double)d_i)) + 1;
					rStart = (rMax<(k-1))? rMax : (k-1);
					rEnd = (rMax<(Depth_mi-1 - k))? rMax : (Depth_mi-1 - k);
					for (n=-rStart; n<=rEnd; n++) {
						w_i = Buffer_i[k+n] + n*n;
						if (w_i<d_i) d_i = w_i;
					}
				}
				
				Distance_mi[Index(i, j, k)] = d_i;
			}
		}
	}

#ifdef	DEBUG_DIST_TF
	printf ("\nStep 3:\n");
	Display_Distance(76);

	printf ("\nStep 3:\n");
	Display_Distance(77);

	printf ("\nStep 3:\n");
	Display_Distance(78);

	printf ("\nStep 3:\n");
	Display_Distance(79);

	printf ("\nStep 3:\n");
	Display_Distance(80);

	printf ("\nStep 3:\n");
	Display_Distance(81);

	printf ("\nStep 3:\n");
	Display_Distance(82);

	printf ("\nStep 3:\n");
	Display_Distance(83);

	printf ("\nStep 3:\n");
	Display_Distance(84);

	printf ("\nStep 3:\n");
	Display_Distance(85);

	printf ("\nStep 3:\n");
	Display_Distance(86);
#endif




	delete [] Buffer_i;
	
}



template<class _DataType>
void cSkeleton<_DataType>::ComputeDistanceVolume(double GradThreshold)
{
	int		i, k, DataCoor[3], FoundZeroCrossingLoc_i;
	double	GradVec_d[3], CurrDataLoc_d[3], ZeroCrossingLoc_d[3], LocalMaxGradient_d;
	double	DataPosFromZeroCrossingLoc_d, Tempd;
	int		CurrZ = 0, Distance_i;
	unsigned char	*DistanceVolume_uc = new unsigned char [WHD_mi];
	

	for (i=0; i<WHD_mi; i++) {
	
		if (GradientMag_mf[i] < GradThreshold) {
			DistanceVolume_uc[i] = (unsigned char)255;
			continue;
		}

		DataCoor[2] = i/WtimesH_mi;
		DataCoor[1] = (i - DataCoor[2]*WtimesH_mi)/Width_mi;
		DataCoor[0] = i % Width_mi;

		if (CurrZ==DataCoor[2]) {
			printf ("Z = %3d, ", DataCoor[2]);
			printf ("\n");
			fflush(stdout);
			CurrZ++;
		}
		
		//-------------------------------------------------------------------------------------------
		// Finding the local maxima of gradient magnitudes along the gradient direction. 
		// It climbs the mountain of gradient magnitudes to find the zero-crossing second derivative 
		// Return ZeroCrossingLoc_d, LocalMaxGradient, DataPosFromZeroCrossingLoc_d
		//-------------------------------------------------------------------------------------------
		// Getting the gradient vector at the position
		for (k=0; k<3; k++) GradVec_d[k] = (double)GradientVec_mf[i*3 + k];
		Tempd = sqrt (GradVec_d[0]*GradVec_d[0] + GradVec_d[1]*GradVec_d[1] + GradVec_d[2]*GradVec_d[2]);
		if (fabs(Tempd)<1e-6) {
			DistanceVolume_uc[i] = (unsigned char)255;
			continue; // To skip zero-length vectors
		}
		for (k=0; k<3; k++) GradVec_d[k] /= Tempd; // Normalize the gradient vector
		for (k=0; k<3; k++) CurrDataLoc_d[k] = (double)DataCoor[k];
		FoundZeroCrossingLoc_i = this->PED_m->FindZeroCrossingLocation(CurrDataLoc_d, GradVec_d, 
									ZeroCrossingLoc_d, LocalMaxGradient_d, DataPosFromZeroCrossingLoc_d);
		if (FoundZeroCrossingLoc_i) {
			// Riversed Distance 0(opaque) - 254(transparent)
			Distance_i = (int)(fabs(DataPosFromZeroCrossingLoc_d)/15.0*254.0);
			DistanceVolume_uc[i] = (unsigned char)Distance_i;
		}
		else DistanceVolume_uc[i] = (unsigned char)255;

	}

	printf ("Computing Distance is done\n");
	printf ("\n");
	fflush(stdout);

	SaveVolume(DistanceVolume_uc, (float)0.0, (float)255.0, "DistByVec");
	printf ("Distance volume (DistByVec) is saved\n"); fflush (stdout);
	
	delete [] DistanceVolume_uc;

}




template<class _DataType>
void cSkeleton<_DataType>::ComputeGVF()
{
	int		i, j, k, loc[3];
	float	Vec_f[3];
	double	Length_d;
	

	printf ("Computing GVF ... \n"); fflush(stdout);
	
	delete [] GVFDistance_mf;
	GVFDistance_mf = new float[WHD_mi*3];

	for (k=0; k<Depth_mi; k++) {
		for (j=0; j<Height_mi; j++) {
			for (i=0; i<Width_mi; i++) {

				loc[0] = Index (i, j, k);

				Vec_f[0] = Distance_mi[Index(i+1, j, k)] - Distance_mi[Index(i-1, j, k)];
				Vec_f[1] = Distance_mi[Index(i, j+1, k)] - Distance_mi[Index(i, j-1, k)];
				Vec_f[2] = Distance_mi[Index(i, j, k+1)] - Distance_mi[Index(i, j, k-1)];

				// Normalize
				Length_d = sqrt((double)Vec_f[0]*Vec_f[0] + Vec_f[1]*Vec_f[1] + Vec_f[2]*Vec_f[2]);
				if (fabs(Length_d)<1e-6) {
					GVFDistance_mf[loc[0]*3 + 0] = (float)0.0; 
					GVFDistance_mf[loc[0]*3 + 1] = (float)0.0;
					GVFDistance_mf[loc[0]*3 + 2] = (float)0.0;
					VoxelFlags_muc[loc[0]] = FLAG_EMPTY; // A zero-length vector
				}
				else {
					GVFDistance_mf[loc[0]*3 + 0] = (float)((double)Vec_f[0]/Length_d);
					GVFDistance_mf[loc[0]*3 + 1] = (float)((double)Vec_f[1]/Length_d);
					GVFDistance_mf[loc[0]*3 + 2] = (float)((double)Vec_f[2]/Length_d);
				}
			}
		}
	}

	// For the Index() function
	VoxelFlags_muc[0] = FLAG_EMPTY;

	GVFDistance_mf[0] = (float)(0.0);
	GVFDistance_mf[1] = (float)(0.0);
	GVFDistance_mf[2] = (float)(0.0);
	
}

// Checking nonuniform voxels by computing dot product with ave. vectors of 2^3 voxels
template<class _DataType>
void cSkeleton<_DataType>::FlagNonUniformGradient()
{
	int				i, j, k, l, m, n, Idx, loc[8];
	int				NumNonunifromVoxels, IsNonuniform;
	float			AveVec_f[3], DotProduct_f, Tempf;


	printf ("\n");
	printf ("Flag NonUniform Gradient() ... \n");
	fflush(stdout);


	// Flagging non-uniform voxels using GVF.
	// When the vector length is zero, then mark the voxels as nonuniform voxels
	for (k=1; k<Depth_mi-1; k++) {
		for (j=1; j<Height_mi-1; j++) {
			for (i=1; i<Width_mi-1; i++) {

				loc[0] = Index (i, j, k);
				if (Distance_mi[loc[0]]<=1) continue;

				if (VoxelFlags_muc[loc[0]]==FLAG_EMPTY || 
					VoxelFlags_muc[loc[0]]==FLAG_NONUNIFORM) continue;


				AveVec_f[0] = AveVec_f[1] = AveVec_f[2] = 0.0;
				Idx = 0;
				for (n=k; n<=k+1; n++) {
					for (m=j; m<=j+1; m++) {
						for (l=i; l<=i+1; l++) {
							loc[Idx] = Index (l, m, n);
							AveVec_f[0] += GVFDistance_mf[loc[Idx]*3 + 0];
							AveVec_f[1] += GVFDistance_mf[loc[Idx]*3 + 1];
							AveVec_f[2] += GVFDistance_mf[loc[Idx]*3 + 2];
							Idx++;
						}
					}
				}
				for (l=0; l<3; l++) AveVec_f[l] /= 8.0;
				// When the vector length is zero, then mark the voxels as nonuniform voxels
				Tempf = sqrtf (AveVec_f[0]*AveVec_f[0] + AveVec_f[1]*AveVec_f[1] + AveVec_f[2]*AveVec_f[2]);
				if (Tempf<1e-6) {
					for (m=0; m<8; m++) VoxelFlags_muc[loc[m]] = 0;
					continue;
				}
				
				IsNonuniform = false;
				for (m=0; m<8; m++) {
					DotProduct_f = 0;
					for (l=0; l<3; l++) {
						DotProduct_f += AveVec_f[l]*GVFDistance_mf[loc[m]*3 + l];
					}
					if (DotProduct_f<=0) {
						IsNonuniform = true;
						break;
					}
				}
				
				
				if (IsNonuniform) {
					for (m=0; m<8; m++) {
						if (Distance_mi[loc[m]]>=1) VoxelFlags_muc[loc[m]] = FLAG_NONUNIFORM;
					}
				}

			}
		}
	}

	NumNonunifromVoxels = 0;
	for (i=0; i<WHD_mi; i++) {
		if (VoxelFlags_muc[i]==FLAG_NONUNIFORM) NumNonunifromVoxels++;

	}
	
	printf ("Num. Nonuniform Voxels = %d\n", NumNonunifromVoxels);
	printf ("Num. Segmented Voxels = %d\n", NumSegmentedVoxels_mi);
	printf ("Num. Nonuniform / NumSegmented = %f %%\n", (float)NumNonunifromVoxels/NumSegmentedVoxels_mi*100.0);
	printf ("Num. Nonuniform / Total Num. Voxels = %f %%\n", (double)NumNonunifromVoxels/WHD_mi*100.0);
	fflush (stdout);
}


template<class _DataType>
void cSkeleton<_DataType>::ConnectingFlaggedVoxels()
{
	int				i, j, k, loc[3], NextPos_i[3], NumContinuousVoxels;
	unsigned char	VoxelFlag_uc;


	printf ("Connecting Flagged Voxels ... \n");
	fflush(stdout);


#ifdef	DEBUG_CONNECTING
//	int		MaxNumContinuousVoxels = 0;
#endif

	for (k=1; k<Depth_mi-1; k++) {
		for (j=1; j<Height_mi-1; j++) {
			for (i=1; i<Width_mi-1; i++) {
	
				loc[0] = Index(i, j, k);
				if (VoxelFlags_muc[loc[0]]==FLAG_NONUNIFORM) {


#ifdef	DEBUG_CONNECTING
					printf ("Start (%d %d %d)\n", i, j, k);
#endif
					NextPos_i[0] = i;
					NextPos_i[1] = j;
					NextPos_i[2] = k;

					NumContinuousVoxels = 0;
					do {

						NextPos_i[0] = (int)((float)NextPos_i[0] + GVFDistance_mf[loc[0]*3 + 0] + 0.5);
						NextPos_i[1] = (int)((float)NextPos_i[1] + GVFDistance_mf[loc[0]*3 + 1] + 0.5);
						NextPos_i[2] = (int)((float)NextPos_i[2] + GVFDistance_mf[loc[0]*3 + 2] + 0.5);

						loc[0] = Index (NextPos_i[0], NextPos_i[1], NextPos_i[2]);
						NumContinuousVoxels++;
						
#ifdef	DEBUG_CONNECTING
						printf ("(%d %d %d)", NextPos_i[0], NextPos_i[1], NextPos_i[2]);
						DisplayFlag(VoxelFlags_muc[loc[0]]);
						printf (",  ");
#endif
						
						VoxelFlag_uc = VoxelFlags_muc[loc[0]];
						if (VoxelFlag_uc == FLAG_CONNECTED || 
							VoxelFlag_uc == FLAG_NONUNIFORM ||
							VoxelFlag_uc == FLAG_EMPTY
							) break;
						
						VoxelFlags_muc[loc[0]] = FLAG_CONNECTED;
						
					} while (VoxelFlag_uc != FLAG_NONUNIFORM);

				}
			}
		}
	}

	for (i=0; i<WHD_mi; i++) {
		if (VoxelFlags_muc[i]==FLAG_NONUNIFORM) VoxelFlags_muc[i] = FLAG_CONNECTED;
	}

	int NumConnectedVoxels = 0;
	for (i=0; i<WHD_mi; i++) {
		if (VoxelFlags_muc[i]==FLAG_CONNECTED) NumConnectedVoxels++;
	}
	printf ("Num. Connected Voxels = %d\n", NumConnectedVoxels);
	printf ("Num. Connected / NumSegmented = %f %%\n", (float)NumConnectedVoxels/NumSegmentedVoxels_mi*100.0);
	printf ("Num. Connected / Total Num. Voxels = %f %%\n", (double)NumConnectedVoxels/WHD_mi*100.0);
	fflush (stdout);
	
}


// Verifying whether the flagged voxels is in blood vessels
// by computing the radii at the flagged voxel locations
template<class _DataType>
void cSkeleton<_DataType>::ConnectingFlaggedVoxels2()
{
	int				i, j, k, loc[3], NextPos_i[3], NumContinuousVoxels;
	unsigned char	VoxelFlag_uc;
	double			StartPt_d[3], Direction_d[3], AveRadius_d;
	double			NewStartPt_Ret[3];


	printf ("Connecting Flagged Voxels 2... \n");
	fflush(stdout);

	for (k=1; k<Depth_mi-1; k++) {
		for (j=1; j<Height_mi-1; j++) {
			for (i=1; i<Width_mi-1; i++) {
	
				loc[0] = Index(i, j, k);
				if (VoxelFlags_muc[loc[0]]==FLAG_NONUNIFORM && 
					Distance_mi[loc[0]]>=2 &&
					Distance_mi[loc[0]]<=4		) {

					StartPt_d[0] = (double)i;
					StartPt_d[1] = (double)j;
					StartPt_d[2] = (double)k;
					
					Direction_d[0] = (double)GVFDistance_mf[loc[0]*3 + 0];
					Direction_d[1] = (double)GVFDistance_mf[loc[0]*3 + 1];
					Direction_d[2] = (double)GVFDistance_mf[loc[0]*3 + 2];
					if (Normalize(Direction_d)<1e-6) continue;


					AveRadius_d = ComputeAveRadius(StartPt_d, Direction_d, NewStartPt_Ret);
		

#ifdef	DEBUG_CONNECTING2
					if (AveRadius_d>0.0) {
						printf ("Start (%3d %3d %3d) --> ", i, j, k);
						printf ("Ave. Radius = %8.3f", AveRadius_d);
						printf ("\n"); fflush (stdout);
					}
					else printf ("\n"); fflush (stdout);
#endif
				

					// If the radius is bigger than the threshold, then remove the flagged voxel
					if (AveRadius_d >= 10.0 || AveRadius_d <= 0.0) {
						VoxelFlags_muc[loc[0]] = FLAG_SEGMENTED;
						continue;
					}


					NextPos_i[0] = i;
					NextPos_i[1] = j;
					NextPos_i[2] = k;

					NumContinuousVoxels = 0;
					do {

						NextPos_i[0] = (int)((float)NextPos_i[0] + GVFDistance_mf[loc[0]*3 + 0] + 0.5);
						NextPos_i[1] = (int)((float)NextPos_i[1] + GVFDistance_mf[loc[0]*3 + 1] + 0.5);
						NextPos_i[2] = (int)((float)NextPos_i[2] + GVFDistance_mf[loc[0]*3 + 2] + 0.5);

						loc[0] = Index (NextPos_i[0], NextPos_i[1], NextPos_i[2]);
						NumContinuousVoxels++;
						
#ifdef	DEBUG_CONNECTING2
//						printf ("(%d %d %d)", NextPos_i[0], NextPos_i[1], NextPos_i[2]);
//						DisplayFlag(VoxelFlags_muc[loc[0]]);
//						printf (",  ");
#endif
						
						VoxelFlag_uc = VoxelFlags_muc[loc[0]];
						if (VoxelFlag_uc == FLAG_CONNECTED || 
							VoxelFlag_uc == FLAG_NONUNIFORM || 
							VoxelFlag_uc == FLAG_EMPTY ||
							Distance_mi[loc[0]] >=10 
							) break;
						
						
						StartPt_d[0] = (double)NextPos_i[0];
						StartPt_d[1] = (double)NextPos_i[1];
						StartPt_d[2] = (double)NextPos_i[2];

						Direction_d[0] = (double)GVFDistance_mf[loc[0]*3 + 0];
						Direction_d[1] = (double)GVFDistance_mf[loc[0]*3 + 1];
						Direction_d[2] = (double)GVFDistance_mf[loc[0]*3 + 2];
						if (Normalize(Direction_d)<1e-6) continue;
						
						AveRadius_d = ComputeAveRadius(StartPt_d, Direction_d, NewStartPt_Ret);

#ifdef	DEBUG_CONNECTING2
						if (AveRadius_d>0.0) {
							printf ("Start (%3d %3d %3d) == ", i, j, k);
							printf ("Ave. Radius = %8.3f", AveRadius_d);
							printf ("\n"); fflush (stdout);
						}
						else printf ("\n"); fflush (stdout);
#endif

						// If the radius is bigger than the threshold, then remove the flagged voxel
						if (AveRadius_d>=10.0 || AveRadius_d <= 0.0) {
							VoxelFlags_muc[loc[0]] = FLAG_SEGMENTED;
							break;
						}
						
						
						VoxelFlags_muc[loc[0]] = FLAG_CONNECTED;
						
					} while (VoxelFlag_uc != FLAG_NONUNIFORM);

				}
			}
		}
	}

	for (i=0; i<WHD_mi; i++) {
		if (VoxelFlags_muc[i]==FLAG_NONUNIFORM && Distance_mi[loc[0]]>=5) VoxelFlags_muc[i] = FLAG_SEGMENTED;
	}

	int NumConnectedVoxels = 0;
	for (i=0; i<WHD_mi; i++) {
		if (VoxelFlags_muc[i]==FLAG_CONNECTED) NumConnectedVoxels++;
	}
	printf ("Num. Connected Voxels = %d\n", NumConnectedVoxels);
	printf ("Num. Connected / NumSegmented = %f %%\n", (float)NumConnectedVoxels/NumSegmentedVoxels_mi*100.0);
	printf ("Num. Connected / Total Num. Voxels = %f %%\n", (double)NumConnectedVoxels/WHD_mi*100.0);
	fflush (stdout);
	
}

template<class _DataType>
double cSkeleton<_DataType>::ComputeAveRadius(double *StartPt_d, double *Direction_d, double *NewStartPt_Ret)
{
	int		l, m;
	double	Rays_d[8*3], HitLocations_d[8*3], Radius_d[8], AveRadius_d;
	double	EndPt_d[3];
	

	//-----------------------------------------------------------------------
	// Computing the new starting point with the flagged voxel location
	// Computing the average radius at the starting point
	ComputePerpendicular8Rays(StartPt_d, Direction_d, Rays_d);
	
	
#ifdef	DEBUG_COMP_AVE_R
	printf ("StartPt = (%7.2f, %7.2f, %7.2f), ", StartPt_d[0], StartPt_d[1], StartPt_d[2]);
	printf ("Direction = (%7.3f, %7.3f, %7.3f), ", Direction_d[0], Direction_d[1], Direction_d[2]);
	fflush (stdout);
#endif
	
	
	// Computing the 8 hit locations with the blood vessel boundary and
	// Returning 8 radii at the starting point to Radius_d and
	ComputeRadius(StartPt_d, Rays_d, HitLocations_d, Radius_d);
	
	

#ifdef	DEBUG_COMP_AVE_R
	printf ("Radii = ");
	for (m=0; m<8; m++) {
		printf ("%8.3f ", Radius_d[m]);
	}
	fflush (stdout);
#endif
	
	
	
	for (l=0; l<3; l++) StartPt_d[l] = 0.0;
	for (m=0; m<8; m++) {
		for (l=0; l<3; l++) StartPt_d[l] += HitLocations_d[m*3 + l];
	}
	for (l=0; l<3; l++) StartPt_d[l] /= 8.0;


	//-----------------------------------------------------------------------
	// Computing the new end point with the flagged voxel location
	for (l=0; l<3; l++) EndPt_d[l] = StartPt_d[l] + Direction_d[l];
	// Computing the average radius at the starting point
	ComputePerpendicular8Rays(EndPt_d, Direction_d, Rays_d);
	// Computing the 8 hit locations with the blood vessel boundary and
	// Returning 8 radii at the starting point to Radius_d and
	ComputeRadius(EndPt_d, Rays_d, HitLocations_d, Radius_d);
	for (l=0; l<3; l++) EndPt_d[l] = 0.0;
	for (m=0; m<8; m++) {
		for (l=0; l<3; l++) EndPt_d[l] += HitLocations_d[m*3 + l];
	}
	for (l=0; l<3; l++) EndPt_d[l] /= 8.0;


	//-----------------------------------------------------------------------
	// Recomputing the hit locations and the average radius with the new 
	// calculated starting point
	for (l=0; l<3; l++) Direction_d[l] = EndPt_d[l] - StartPt_d[l];
	if (Normalize(Direction_d)<1e-6) {
		AveRadius_d = -1.0;
		for (l=0; l<3; l++) NewStartPt_Ret[l] = -1.0;
		return AveRadius_d;
	}
	
	// Computing the average radius at the starting point
	ComputePerpendicular8Rays(StartPt_d, Direction_d, Rays_d);
	// Computing the 8 hit locations with the blood vessel boundary and
	// Returning 8 radii at the starting point to Radius_d and
	ComputeRadius(StartPt_d, Rays_d, HitLocations_d, Radius_d);
	
	
	int		NumPositiveRadii = 0;
	
	AveRadius_d = 0.0;
	for (m=0; m<8; m++) {
		if (Radius_d[m]>0) {
			AveRadius_d += Radius_d[m];
			NumPositiveRadii++;
		}
	}
	// The blood vessel should not be open
	// When the number of the positive radii is small (< 5),
	// I say that the blood vessel is open.
	// Note: Negative radius is computed, when there is no nearest zero f''
	if (NumPositiveRadii>=5) AveRadius_d /= (double)NumPositiveRadii;
	else AveRadius_d = -1.0;

	for (l=0; l<3; l++) NewStartPt_Ret[l] = StartPt_d[l];


#ifdef	DEBUG_COMP_AVE_R
	printf ("Average = %8.3f", AveRadius_d);
	printf ("\n"); fflush (stdout);
#endif

	
	return AveRadius_d;
	
}



template<class _DataType>
void cSkeleton<_DataType>::ComputeRadius(double *StartPt, double *Rays8, double *HitLocs8, double *Radius8)
{
	int		i, k, Found;
	double	ZeroCrossingLoc_Ret[3], FirstDAtTheLoc_Ret, DataPosFromZeroCrossingLoc_Ret;

	
	for (i=0; i<8; i++) {
		Found = getANearestZeroCrossing(StartPt, &Rays8[i*3], ZeroCrossingLoc_Ret, 
									FirstDAtTheLoc_Ret, DataPosFromZeroCrossingLoc_Ret);
		if (Found) {
			for (k=0; k<3; k++) HitLocs8[i*3 + k] = ZeroCrossingLoc_Ret[k];
			Radius8[i] = DataPosFromZeroCrossingLoc_Ret;
		}
		else {
			for (k=0; k<3; k++) HitLocs8[i*3 + k] = -1.0;
			Radius8[i] = -1.0;
		}
	}

}





template<class _DataType>
void cSkeleton<_DataType>::ConnectedComponents(char *OutFileName, int &MaxCCLoc_Ret)
{
	int 			i, j, k, X_i, Y_i, Z_i, loc[3], NumCCVoxels;
	map<int, int>				CCList_m;
	map<int, int>::iterator		CCList_it;
	map<int, unsigned char>				CC_m;
	map<int, unsigned char>::iterator	CC_it;
	

	CCVolume_muc = new unsigned char[WHD_mi];
	
	// Copying the voxel flag volume
	for (i=0; i<WHD_mi; i++) {
		if (VoxelFlags_muc[i]==FLAG_CONNECTED) CCVolume_muc[i] = VoxelFlags_muc[i];
		else CCVolume_muc[i] = (unsigned char)FLAG_EMPTY;
	}
	
	//------------------------------------------------------------------------------
	// Extracting CC, or 
	// grouping voxels with the connected information
	//------------------------------------------------------------------------------
	CCList_m.clear();
	CC_m.clear();
	for (k=1; k<Depth_mi-1; k++) {
		for (j=1; j<Height_mi-1; j++) {
			for (i=1; i<Width_mi-1; i++) {
	
				loc[0] = k*WtimesH_mi + j*Width_mi + i;

				if (CCVolume_muc[loc[0]]==FLAG_CONNECTED) {
					NumCCVoxels = MarkingCC(loc[0], (unsigned char)FLAG_EMPTY, CCVolume_muc);
					CCList_m[loc[0]] = NumCCVoxels;
				}

			}
		}
	}


	//------------------------------------------------------------------------------
	// Finding a CC that has the max number of voxels and
	// saving the coordinate to MaxLocX, MaxLocY, and MaxLocZ.
	//------------------------------------------------------------------------------
	int		TotalNumCCVoxels = 0;
	int		MinNumVoxels=99999, MaxNumVoxels=-99999;
	int		MaxLocX, MaxLocY, MaxLocZ;
	
	CCList_it = CCList_m.begin();
	for (i=0; i<CCList_m.size(); i++, CCList_it++) {

		loc[0] =(*CCList_it).first;
		NumCCVoxels = (*CCList_it).second;
		TotalNumCCVoxels += NumCCVoxels;
		
		Z_i = loc[0]/WtimesH_mi;
		Y_i = (loc[0] - Z_i*WtimesH_mi)/Width_mi;
		X_i = loc[0] % Width_mi;

		if (MinNumVoxels > NumCCVoxels) MinNumVoxels = NumCCVoxels;
		if (MaxNumVoxels < NumCCVoxels) { 
			MaxNumVoxels = NumCCVoxels; 
			MaxCCLoc_Ret = loc[0];
			MaxLocX = X_i;
			MaxLocY = Y_i;
			MaxLocZ = Z_i;
		}

		if (NumCCVoxels>10) {
			printf ("(%3d %3d %3d): ", X_i, Y_i, Z_i);
			printf ("Num. Connected Voxels = %d\n", NumCCVoxels);
		}
	
	}

	printf ("Min & Max Num. Voxels = %d %d\n", MinNumVoxels, MaxNumVoxels);
	printf ("Max Loc xyz = (%d, %d %d)\n", MaxLocX, MaxLocY, MaxLocZ);
	printf ("Num. of CC = %d\n", (int)CCList_m.size());
	printf ("Total Num. CC Voxles = %d\n", TotalNumCCVoxels);
	fflush (stdout);


	// Copying the voxel flag volume
	for (i=0; i<WHD_mi; i++) {
		if (VoxelFlags_muc[i]==FLAG_CONNECTED) CCVolume_muc[i] = VoxelFlags_muc[i];
		else CCVolume_muc[i] = (unsigned char)FLAG_EMPTY;
	}

	
	//------------------------------------------------------------------------------
	// Generating the CC volume, CCVolume_uc
	// CCList_m is cleared at this stage.
	//------------------------------------------------------------------------------
	int				CurrMaxLoc;
	unsigned char	GreyColor_uc;
	
	GreyColor_uc = 255;
	
	do {
	
		MaxNumVoxels=-99999;
		CCList_it = CCList_m.begin();
		for (i=0; i<CCList_m.size(); i++, CCList_it++) {

			loc[0] =(*CCList_it).first;
			NumCCVoxels = (*CCList_it).second;
			if (MaxNumVoxels < NumCCVoxels) {
				MaxNumVoxels = NumCCVoxels;
				CurrMaxLoc = loc[0];
			}
		}
		
		CCList_m.erase(CurrMaxLoc);
		NumCCVoxels = MarkingCC(CurrMaxLoc, GreyColor_uc, CCVolume_muc);
		Z_i = CurrMaxLoc/WtimesH_mi;
		Y_i = (CurrMaxLoc - Z_i*WtimesH_mi)/Width_mi;
		X_i = CurrMaxLoc % Width_mi;

		printf ("Loc = (%3d %3d %3d), ", X_i, Y_i, Z_i);
		printf ("Grey Color = %3d, ", GreyColor_uc);
		printf ("Num. Voxels = %6d, ", MaxNumVoxels);
		printf ("Num. Voxels (rechecked) = %6d\n", NumCCVoxels);
		fflush (stdout);
		
		GreyColor_uc -= 2;
		if (GreyColor_uc==FLAG_CONNECTED) GreyColor_uc -=2;
		
	} while (CCList_m.size()>0 && GreyColor_uc>5 );

	
}


template<class _DataType>
int cSkeleton<_DataType>::MarkingCC(int CCLoc, unsigned char MarkingNum, unsigned char *CCVolume_uc)
{
	int			i, j, k, X_i, Y_i, Z_i, loc[3], NumCCVoxels;
	map<int, unsigned char>				CC_m;
	map<int, unsigned char>::iterator	CC_it;

	
	CC_m.clear();
	loc[0] = CCLoc;
	CC_m[loc[0]] = (unsigned char)1;
	CCVolume_uc[loc[0]] = MarkingNum;
	NumCCVoxels = 1;
	
	do {

		CC_it = CC_m.begin();
		loc[0] = (*CC_it).first;
		CC_m.erase(loc[0]);
		Z_i = loc[0]/WtimesH_mi;
		Y_i = (loc[0] - Z_i*WtimesH_mi)/Width_mi;
		X_i = loc[0] % Width_mi;

		for (k=Z_i-1; k<=Z_i+1; k++) {
			for (j=Y_i-1; j<=Y_i+1; j++) {
				for (i=X_i-1; i<=X_i+1; i++) {

					loc[0] = Index(i, j, k);
					if (CCVolume_uc[loc[0]]==FLAG_CONNECTED) {
						CC_m[loc[0]] = (unsigned char)1;
						CCVolume_uc[loc[0]] = MarkingNum;
						NumCCVoxels++;
					}

				}
			}
		}

	} while (CC_m.size()>0);
	CC_m.clear();
	
	return NumCCVoxels;
}


template<class _DataType>
void cSkeleton<_DataType>::FindingRootAndEndLoc(int MaxCCLoc, int &RootLoc_Ret, int &EndLoc_Ret)
{
	int			i, j, k, X_i, Y_i, Z_i, loc[3];
	int			Root_i[3];


	printf ("Finding the root and end voxels...\n"); fflush (stdout);
	
	//-----------------------------------------------------------------------
	// Finding the root voxel of skeleton
	// The root voxel has the smallest Y coordinate, which means that
	// it is the front of the heart
	Root_i[0] = 0; 			// X Coordinate
	Root_i[1] = Height_mi;	// Y Coordinate
	Root_i[2] = 0; 			// Z coordinate
	
	for (j=0; j<Height_mi; j++) {
		for (k=0; k<Depth_mi; k++) {
			for (i=0; i<Width_mi; i++) {

				if (CCVolume_muc[Index(i, j, k)]==255) {
					Root_i[0] = i;
					Root_i[1] = j;
					Root_i[2] = k;
					i = j = k = WHD_mi; // Exiting this three overlapped loops
				}
				
			}
		}
	}
	
	RootLoc_Ret = Index (Root_i[0], Root_i[1], Root_i[2]);
	printf ("The root voxel = (%d, %d, %d)\n", Root_i[0], Root_i[1], Root_i[2]); fflush (stdout);
	//-----------------------------------------------------------------------
	

	//-----------------------------------------------------------------------
	// Finding the end voxels of skeleton
	// The end voxel is the farthest voxel from the root voxel
	map<int, unsigned char>				CurrBoundary_m;
	
	
	DistanceFromSkeletons_mi = new int [WHD_mi];

	loc[0] = RootLoc_Ret;
	CurrBoundary_m.clear();
	CurrBoundary_m[loc[0]] = (unsigned char)1; // Adding the current location
	EndLoc_Ret = ComputeDistanceFromCurrSkeletons (CurrBoundary_m);
	CurrBoundary_m.clear();

	Z_i = EndLoc_Ret/WtimesH_mi;
	Y_i = (EndLoc_Ret - Z_i*WtimesH_mi)/Width_mi;
	X_i = EndLoc_Ret % Width_mi;
	printf ("The end voxel = (%d, %d, %d)\n", X_i, Y_i, Z_i);
	//-----------------------------------------------------------------------
	

	int		Maxi=-9999999, Mini=9999999;
	for (i=0; i<WHD_mi; i++) {
		if (Mini > DistanceFromSkeletons_mi[i]) Mini = DistanceFromSkeletons_mi[i];
		if (Maxi < DistanceFromSkeletons_mi[i]) Maxi = DistanceFromSkeletons_mi[i];
	}
	printf ("Min & Max distance from the root = %d, %d\n", Mini, Maxi); fflush (stdout);
	
}


template<class _DataType>
void cSkeleton<_DataType>::ComputeSkeletons(int RootLoc, int EndLoc, int Threshold_Dist)
{
	int		i, j, k, X_i, Y_i, Z_i, loc[3];
	int		FurthestVoxelLoc, FurthestVoxelDistance_i;
	map<int, unsigned char>		CurrBoundary_m;
	cStack<int>					EndLocs_stack;
	

	// Initializing the skeleton volume data
	Skeletons_muc = new unsigned char [WHD_mi];
	for (i=0; i<WHD_mi; i++) Skeletons_muc[i] = (unsigned char)0;
	Skeletons_muc[RootLoc] = (unsigned char)255;
	Skeletons_muc[EndLoc] = (unsigned char)255;


	FurthestVoxelDistance_i = SkeletonTracking_Marking(EndLoc);
	printf ("Furthest Voxel Distance = %d\n", FurthestVoxelDistance_i); fflush (stdout);


	//--------------------------------------------------------------------------------------
	// Computing other branches of the skeletons
	EndLocs_stack.Clear();
	do {
	
		CurrBoundary_m.clear();
		for (i=0; i<WHD_mi; i++) {
			if (Skeletons_muc[i]==255) CurrBoundary_m[i] = (unsigned char)0;
		}
		FurthestVoxelLoc = ComputeDistanceFromCurrSkeletons(CurrBoundary_m);
		FurthestVoxelDistance_i = SkeletonTracking_Marking(FurthestVoxelLoc);
		EndLocs_stack.Push(FurthestVoxelLoc);
		
		printf ("Furthest Voxel Distance = %d\n", FurthestVoxelDistance_i); fflush (stdout);
		
	} while (FurthestVoxelDistance_i >= Threshold_Dist);
	
	printf ("Num. end points = %d\n", EndLocs_stack.Size()); fflush (stdout);
	//--------------------------------------------------------------------------------------
	

	// Marking the root voxel
	Z_i = RootLoc/WtimesH_mi;
	Y_i = (RootLoc - Z_i*WtimesH_mi)/Width_mi;
	X_i = RootLoc % Width_mi;
	for (k=Z_i-5; k<=Z_i+5; k++) {
		for (j=Y_i-5; j<=Y_i+5; j++) {
			for (i=X_i-5; i<=X_i+5; i++) {

				loc[0] = Index (i, j, k);
				Skeletons_muc[loc[0]] = (unsigned char)245;
			}
		}
	}

	// Marking the end voxel
	Z_i = EndLoc/WtimesH_mi;
	Y_i = (EndLoc - Z_i*WtimesH_mi)/Width_mi;
	X_i = EndLoc % Width_mi;
	for (k=Z_i-3; k<=Z_i+3; k++) {
		for (j=Y_i-3; j<=Y_i+3; j++) {
			for (i=X_i-3; i<=X_i+3; i++) {

				loc[0] = Index (i, j, k);
				Skeletons_muc[loc[0]] = (unsigned char)240;
			}
		}
	}


	// Marking the end voxles of the skeletons
	int		GreyColor_i = 235;
		
	do {

		EndLocs_stack.Pop(loc[0]);
		Z_i = loc[0]/WtimesH_mi;
		Y_i = (loc[0] - Z_i*WtimesH_mi)/Width_mi;
		X_i = loc[0] % Width_mi;
		for (k=Z_i-2; k<=Z_i+2; k++) {
			for (j=Y_i-2; j<=Y_i+2; j++) {
				for (i=X_i-2; i<=X_i+2; i++) {

					loc[0] = Index (i, j, k);
					Skeletons_muc[loc[0]] = (unsigned char)GreyColor_i;
				}
			}
		}
		GreyColor_i -= 1;
		printf ("The grey colors for the end voxels = %d\n", GreyColor_i); fflush (stdout);
		if (GreyColor_i<60) GreyColor_i = 60;
		
	} while (!EndLocs_stack.IsEmpty());


	// Combining the initial segmented volume with the skeletons
	for (i=0; i<WHD_mi; i++) {
		if (InitVolume_muc[i]==255 && Skeletons_muc[i]==0) {
			Skeletons_muc[i] = 50;
		}
	}
	
	
}


template<class _DataType>
int cSkeleton<_DataType>::SkeletonTracking_Marking(int EndLoc)
{
	int		i, j, k, X_i, Y_i, Z_i, loc[3];
	int		MinDistLoc, CurrDistanceFromSkeletons_i;
	int		FurthestDistance_i;


	Skeletons_muc[EndLoc] = (unsigned char)255;
	FurthestDistance_i = DistanceFromSkeletons_mi[EndLoc];
	CurrDistanceFromSkeletons_i = DistanceFromSkeletons_mi[EndLoc];
	MinDistLoc = EndLoc;
	
	do {
		Z_i = MinDistLoc/WtimesH_mi;
		Y_i = (MinDistLoc - Z_i*WtimesH_mi)/Width_mi;
		X_i = MinDistLoc % Width_mi;

		for (k=Z_i-1; k<=Z_i+1; k++) {
			for (j=Y_i-1; j<=Y_i+1; j++) {
				for (i=X_i-1; i<=X_i+1; i++) {
				
					loc[0] = Index (i, j, k);
					if (CurrDistanceFromSkeletons_i > DistanceFromSkeletons_mi[loc[0]] && 
						DistanceFromSkeletons_mi[loc[0]]>=0) {
						
						CurrDistanceFromSkeletons_i = DistanceFromSkeletons_mi[loc[0]];
						MinDistLoc = loc[0];
						Skeletons_muc[loc[0]] = (unsigned char)255;

						i = j = k = Width_mi; // To exit the three overlapped loops
						
					}
				
				}
			}
		}
	
	} while(CurrDistanceFromSkeletons_i > 1);


	return FurthestDistance_i;

}


template<class _DataType>
int cSkeleton<_DataType>::ComputeDistanceFromCurrSkeletons(map<int, unsigned char> &InitialBoundary_m)
{
	int		i, j, k, l, X_i, Y_i, Z_i, loc[3];
	map<int, unsigned char>				CurrBoundary_m, NextBoundary_m;
	map<int, unsigned char>::iterator	Boundary_it;
	int		CurrDistance_i;


	CurrBoundary_m.clear();
	NextBoundary_m.clear();
	for (i=0; i<WHD_mi; i++) DistanceFromSkeletons_mi[i] = -1;
	CurrDistance_i = 0;
	
	Boundary_it = InitialBoundary_m.begin();
	for (i=0; i<InitialBoundary_m.size(); i++, Boundary_it++) {
		loc[0] = (*Boundary_it).first;
		CurrBoundary_m[loc[0]] = (unsigned char)1;
		DistanceFromSkeletons_mi[loc[0]] = CurrDistance_i;
	}

	do {
		
		CurrDistance_i++;
		Boundary_it = CurrBoundary_m.begin();
		for (l=0; l<CurrBoundary_m.size(); l++, Boundary_it++) {

			loc[0] = (*Boundary_it).first;
			Z_i = loc[0]/WtimesH_mi;
			Y_i = (loc[0] - Z_i*WtimesH_mi)/Width_mi;
			X_i = loc[0] % Width_mi;

			for (k=Z_i-1; k<=Z_i+1; k++) {
				for (j=Y_i-1; j<=Y_i+1; j++) {
					for (i=X_i-1; i<=X_i+1; i++) {

						loc[0] = Index(i, j, k);
						if (CCVolume_muc[loc[0]]==255 && DistanceFromSkeletons_mi[loc[0]]<0) {
							NextBoundary_m[loc[0]] = (unsigned char)1; // Adding the loc to the curr map
							DistanceFromSkeletons_mi[loc[0]] = CurrDistance_i;
						}

					}
				}
			}

		}
		
		if (NextBoundary_m.size()==0) break;
		CurrBoundary_m.clear();
		
		Boundary_it = NextBoundary_m.begin();
		for (l=0; l<NextBoundary_m.size(); l++, Boundary_it++) {
			CurrBoundary_m[(*Boundary_it).first] = (unsigned char)1;
		}
		NextBoundary_m.clear();
		
	} while (1);


	// returning one of the furthest voxel locations
	Boundary_it = CurrBoundary_m.begin();
	loc[0] = (*Boundary_it).first;
	return loc[0];

}




template <class _DataType>
void cSkeleton<_DataType>::ComputePerpendicular8Rays(double *StartPt, double *Direction, double *Rays)
{
	int		i, k;
	double	PointOnPlane_d[3], Weight_d, Tempd;


	// Compute a plane equation and a perpendicular ray direction to the vessel direction
	Weight_d = -(StartPt[0]*Direction[0] + StartPt[1]*Direction[1] + StartPt[2]*Direction[2]);
	
	if (fabs(Direction[2])>1e-5) {
		PointOnPlane_d[0] = StartPt[0] + 10.0;
		PointOnPlane_d[1] = StartPt[1];
	}
	else if (fabs(Direction[0])>1e-5) {
		PointOnPlane_d[0] = StartPt[0];
		PointOnPlane_d[1] = StartPt[1] + 10.0;
	} 
	else if (fabs(Direction[1])>1e-5) {
		PointOnPlane_d[0] = StartPt[0] + 10.0;
		PointOnPlane_d[1] = StartPt[1];
	} 
	
	if (fabs(Direction[2])<1e-5) {
		PointOnPlane_d[2] = 0.0;
	}
	else {
		PointOnPlane_d[2] = (-PointOnPlane_d[0]*Direction[0] - 
								PointOnPlane_d[1]*Direction[1] - Weight_d)/Direction[2];
	}

	// Compute Ray0
	for (k=0; k<3; k++) Rays[k] = PointOnPlane_d[k] - StartPt[k];

	// Compute Ray1 with the cross product
	Rays[1*3 + 0] = Direction[1]*Rays[0*3 + 2] - Direction[2]*Rays[0*3 + 1];
	Rays[1*3 + 1] = Direction[2]*Rays[0*3 + 0] - Direction[0]*Rays[0*3 + 2];
	Rays[1*3 + 2] = Direction[0]*Rays[0*3 + 1] - Direction[1]*Rays[0*3 + 0];

	// Compute Ray2
	for (k=0; k<3; k++) Rays[2*3 + k] = (Rays[k] + Rays[1*3 + k])/2.0;

	// Compute Ray3 with the cross product
	Rays[3*3 + 0] = Direction[1]*Rays[2*3 + 2]-Direction[2]*Rays[2*3 + 1];
	Rays[3*3 + 1] = Direction[2]*Rays[2*3 + 0]-Direction[0]*Rays[2*3 + 2];
	Rays[3*3 + 2] = Direction[0]*Rays[2*3 + 1]-Direction[1]*Rays[2*3 + 0];

	// Compute Ray4 with the cross product
	for (k=0; k<3; k++) Rays[4*3 + k] = (Rays[1*3 + k] + Rays[2*3 + k])/2.0; // R4 = (R1+R2)/2
	for (k=0; k<3; k++) Rays[5*3 + k] = (Rays[2*3 + k] + Rays[0*3 + k])/2.0; // R5 = (R0+R2)/2
	for (k=0; k<3; k++) Rays[6*3 + k] = (Rays[0*3 + k] - Rays[3*3 + k])/2.0; // R6 = (R0+R3)/2
	for (k=0; k<3; k++) Rays[7*3 + k] = (Rays[3*3 + k] + Rays[1*3 + k])/2.0; // R7 = (R3-R1)/2

	for (i=0; i<=7; i++) {
		Tempd = sqrt(Rays[i*3 + 0]*Rays[i*3 + 0] + Rays[i*3 + 1]*Rays[i*3 + 1] + Rays[i*3 + 2]*Rays[i*3 + 2]);
		for (k=0; k<3; k++) Rays[i*3 + k] /= Tempd; // Normalize the ray from 4-7
	}


}


template <class _DataType>
int cSkeleton<_DataType>::getANearestZeroCrossing(double *CurrLoc, double *DirVec, double *ZeroCrossingLoc_Ret, 
									double& FirstDAtTheLoc_Ret, double& DataPosFromZeroCrossingLoc_Ret)
{
	int			k, OutOfVolume;
	double		FirstD_d[3], SecondD_d[3], Step;
	double		StartPos, EndPos, MiddlePos, LocAlongGradDir[3], Increase_d;
	double 		ThresholdGM_d = 5.0;

	Increase_d = 0.2;
	//--------------------------------------------------------------------------------------
	// The Positive and Negative Direction of a Gradient Vector: Repeat twice
	// for the Positive and Negative Directions
	//--------------------------------------------------------------------------------------

	for (k=0; k<3; k++) LocAlongGradDir[k] = CurrLoc[k] + DirVec[k]*Increase_d;
	FirstD_d[0] = GradientInterpolation(LocAlongGradDir);
	SecondD_d[0] = SecondDInterpolation(LocAlongGradDir);
	if (fabs(SecondD_d[0])<1e-5 && FirstD_d[0]>=ThresholdGM_d) {
		for (k=0; k<3; k++) ZeroCrossingLoc_Ret[k] = LocAlongGradDir[k];
		FirstDAtTheLoc_Ret = FirstD_d[0];
		DataPosFromZeroCrossingLoc_Ret = Increase_d;
		return true;
	}

	for (k=0; k<3; k++) LocAlongGradDir[k] = CurrLoc[k] + DirVec[k]*Increase_d*2.0;
	FirstD_d[1] = GradientInterpolation(LocAlongGradDir);
	SecondD_d[1] = SecondDInterpolation(LocAlongGradDir);
	if (fabs(SecondD_d[1])<1e-5 && FirstD_d[1]>=ThresholdGM_d) {
		for (k=0; k<3; k++) ZeroCrossingLoc_Ret[k] = LocAlongGradDir[k];
		FirstDAtTheLoc_Ret = FirstD_d[1];
		DataPosFromZeroCrossingLoc_Ret = Increase_d*2.0;
		return true;
	}

	if (SecondD_d[0]*SecondD_d[1]<0.0 && FirstD_d[0]>=ThresholdGM_d/2.0 && FirstD_d[1]>=ThresholdGM_d/2.0) {
		StartPos = Increase_d;
		EndPos = Increase_d*2.0;
		MiddlePos = (StartPos + EndPos)/2.0;
		SecondD_d[2] = SecondD_d[1];
		// Binary Search of the zero-crossing location
		do {
			for (k=0; k<3; k++) LocAlongGradDir[k] = (CurrLoc[k]) + DirVec[k]*MiddlePos;
			FirstD_d[1] = GradientInterpolation(LocAlongGradDir);
			SecondD_d[1] = SecondDInterpolation(LocAlongGradDir);
			if (fabs(SecondD_d[1])<1e-5) {
				if (FirstD_d[1]>=ThresholdGM_d) {
					for (k=0; k<3; k++) ZeroCrossingLoc_Ret[k] = LocAlongGradDir[k];
					FirstDAtTheLoc_Ret = FirstD_d[1];
					DataPosFromZeroCrossingLoc_Ret = MiddlePos;
					return true;
				}
				else break;
			}
			if (SecondD_d[1]*SecondD_d[2]<0.0) {
				StartPos = MiddlePos;
				MiddlePos = (StartPos + EndPos)/2.0;
				SecondD_d[0] = SecondD_d[1];
			}
			if (SecondD_d[0]*SecondD_d[1]<0.0) {
				EndPos = MiddlePos;
				MiddlePos = (StartPos + EndPos)/2.0;
				SecondD_d[2] = SecondD_d[1];
			}
		} while (fabs(StartPos-EndPos)>1e-5);

		// When the gradient magnitude is less than 10.0,
		// recompute f' and f''.
		for (k=0; k<3; k++) LocAlongGradDir[k] = CurrLoc[k] + DirVec[k]*Increase_d;
		FirstD_d[0] = GradientInterpolation(LocAlongGradDir);
		SecondD_d[0] = SecondDInterpolation(LocAlongGradDir);

		for (k=0; k<3; k++) LocAlongGradDir[k] = CurrLoc[k] + DirVec[k]*Increase_d*2.0;
		FirstD_d[1] = GradientInterpolation(LocAlongGradDir);
		SecondD_d[1] = SecondDInterpolation(LocAlongGradDir);

	}

	for (Step=Increase_d*3.0; fabs(Step)<=15.0+1e-5; Step+=Increase_d) {
		OutOfVolume = false;
		for (k=0; k<3; k++) {
			LocAlongGradDir[k] = CurrLoc[k] + DirVec[k]*Step;
			if (LocAlongGradDir[k]<0.0) OutOfVolume = true;
		}
		if (OutOfVolume) break;
		if (LocAlongGradDir[0]>=(double)Width_mi) break;
		if (LocAlongGradDir[1]>=(double)Height_mi) break;
		if (LocAlongGradDir[2]>=(double)Depth_mi) break;
		FirstD_d[2] = GradientInterpolation(LocAlongGradDir);
		SecondD_d[2] = SecondDInterpolation(LocAlongGradDir);

		if (fabs(SecondD_d[2])<1e-5 && FirstD_d[2]>=ThresholdGM_d) {
			for (k=0; k<3; k++) ZeroCrossingLoc_Ret[k] = LocAlongGradDir[k];
			FirstDAtTheLoc_Ret = FirstD_d[2];
			DataPosFromZeroCrossingLoc_Ret = Step;
			return true;
		}
		if (SecondD_d[1]*SecondD_d[2]<0.0 && FirstD_d[1]>=ThresholdGM_d/2.0 && FirstD_d[2]>=ThresholdGM_d/2.0) {
			StartPos = Step - Increase_d;
			EndPos = Step;
			MiddlePos = (StartPos + EndPos)/2.0;
			SecondD_d[0] = SecondD_d[1];
			// Binary Search of the zero-crossing location
			do {
				for (k=0; k<3; k++) LocAlongGradDir[k] = CurrLoc[k] + DirVec[k]*MiddlePos;
				FirstD_d[1] = GradientInterpolation(LocAlongGradDir);
				SecondD_d[1] = SecondDInterpolation(LocAlongGradDir);
				if (fabs(SecondD_d[1])<1e-5) {
					if (FirstD_d[1]>=ThresholdGM_d) {
						for (k=0; k<3; k++) ZeroCrossingLoc_Ret[k] = LocAlongGradDir[k];
						FirstDAtTheLoc_Ret = FirstD_d[1];
						DataPosFromZeroCrossingLoc_Ret = MiddlePos;
						return true;
					}
					else break;
				}
				
				if (SecondD_d[1]*SecondD_d[2]<0.0) {
					StartPos = MiddlePos;
					MiddlePos = (StartPos + EndPos)/2.0;
					SecondD_d[0] = SecondD_d[1];
				}
				if (SecondD_d[0]*SecondD_d[1]<0.0) {
					EndPos = MiddlePos;
					MiddlePos = (StartPos + EndPos)/2.0;
					SecondD_d[2] = SecondD_d[1];
				}
			} while (fabs(StartPos-EndPos)>1e-5);

		}

		for (k=0; k<=1; k++) FirstD_d[k] = FirstD_d[k+1];
		for (k=0; k<=1; k++) SecondD_d[k] = SecondD_d[k+1];


	}

	
	return false;
}



template <class _DataType>
double cSkeleton<_DataType>::Normalize(double	*Vec)
{
	double	Length = sqrt(Vec[0]*Vec[0] + Vec[1]*Vec[1] + Vec[2]*Vec[2]);
	for (int k=0; k<3; k++) Vec[k] /= Length;
	return Length;
}


// Trilinear Interpolation
/*
Vxyz =  V000 (1 - Vx) (1 - Vy) (1 - Vz) +
		V100 Vx (1 - Vy) (1 - Vz) + 
		V010 (1 - Vx) Vy (1 - Vz) + 
		V110 Vx Vy (1 - Vz) + 
		V001 (1 - Vx) (1 - Vy) Vz +
		V101 Vx (1 - Vy) Vz + 
		V011 (1 - Vx) Vy Vz + 
		V111 Vx Vy Vz  
*/
// The vector (Vx, Vy, Vz) should have unit length or 1.
// Tri-linear Interpolation
// (Vx, Vy, Vz) is a location in the volume dataset
template<class _DataType>
double cSkeleton<_DataType>::GradientInterpolation(double* LocXYZ)
{
	return GradientInterpolation(LocXYZ[0], LocXYZ[1], LocXYZ[2]);
}

template<class _DataType>
double cSkeleton<_DataType>::GradientInterpolation(double LocX, double LocY, double LocZ)
{
	int		i, j, k, loc[2], X, Y, Z;
	double	RetGradM, GradM[8], Vx, Vy, Vz;


	X = (int)floor(LocX+1e-8);
	Y = (int)floor(LocY+1e-8);
	Z = (int)floor(LocZ+1e-8);
	Vx = LocX - (double)X;
	Vy = LocY - (double)Y;
	Vz = LocZ - (double)Z;

	for (i=0; i<8; i++) GradM[i] = 0.0;
	loc[1] = 0;
	for (k=Z; k<=Z+1; k++) {
		for (j=Y; j<=Y+1; j++) {
			for (i=X; i<=X+1; i++) {
				if (i<0 || j<0 || k<0 || i>=Width_mi || j>=Height_mi ||  k>=Depth_mi) loc[1]++;
				else {
					loc[0] = k*WtimesH_mi + j*Width_mi + i;
					GradM[loc[1]] = (double)GradientMag_mf[loc[0]];
					loc[1]++;
				}
			}
		}
	}

	loc[1] = 0;
	RetGradM = 0.0;
	for (k=0; k<=1; k++) {
		for (j=0; j<=1; j++) {
			for (i=0; i<=1; i++) {
				RetGradM += GradM[loc[1]]*	((double)(1-i) - Vx*pow((double)-1.0, (double)i))*
											((double)(1-j) - Vy*pow((double)-1.0, (double)j))*
											((double)(1-k) - Vz*pow((double)-1.0, (double)k));
				loc[1] ++;
			}
		}
	}

	return RetGradM;
}



// Tri-linear Interpolation
// (Vx, Vy, Vz) is a location in the volume dataset
template<class _DataType>
double cSkeleton<_DataType>::SecondDInterpolation(double* LocXYZ)
{
	return SecondDInterpolation(LocXYZ[0], LocXYZ[1], LocXYZ[2]);
}

template<class _DataType>
double cSkeleton<_DataType>::SecondDInterpolation(double LocX, double LocY, double LocZ)
{
	int		i, j, k, loc[3], X, Y, Z;
	double	RetSecondD, SecondD[8], Vx, Vy, Vz;


	X = (int)floor(LocX+1e-8);
	Y = (int)floor(LocY+1e-8);
	Z = (int)floor(LocZ+1e-8);
	Vx = LocX - (double)X;
	Vy = LocY - (double)Y;
	Vz = LocZ - (double)Z;

	for (i=0; i<8; i++) SecondD[i] = 0.0;
	loc[1] = 0;
	for (k=Z; k<=Z+1; k++) {
		for (j=Y; j<=Y+1; j++) {
			for (i=X; i<=X+1; i++) {
				if (i<0 || j<0 || k<0 || i>=Width_mi || j>=Height_mi ||  k>=Depth_mi) loc[1]++;
				else {
					loc[0] = k*WtimesH_mi + j*Width_mi + i;
					SecondD[loc[1]] = (double)SecondDerivative_mf[loc[0]];
					loc[1]++;
				}
			}
		}
	}

	loc[1] = 0;
	RetSecondD = 0.0;
	for (k=0; k<=1; k++) {
		for (j=0; j<=1; j++) {
			for (i=0; i<=1; i++) {
				RetSecondD += SecondD[loc[1]]*	((double)(1-i) - Vx*pow((double)-1.0, (double)i))*
												((double)(1-j) - Vy*pow((double)-1.0, (double)j))*
												((double)(1-k) - Vz*pow((double)-1.0, (double)k));
				loc[1] ++;
			}
		}
	}

	return RetSecondD;
}



template<class _DataType>
int	cSkeleton<_DataType>::Index(int X, int Y, int Z, int ith, int NumElements)
{
	if (X<0 || Y<0 || Z<0 || X>=Width_mi || Y>=Height_mi || Z>=Depth_mi) return 0;
	else return (Z*WtimesH_mi*NumElements + Y*Width_mi*NumElements + X*NumElements + ith);
}

template<class _DataType>
int	cSkeleton<_DataType>::Index(int X, int Y, int Z)
{
	if (X<0 || Y<0 || Z<0 || X>=Width_mi || Y>=Height_mi || Z>=Depth_mi) return 0;
	else return (Z*WtimesH_mi + Y*Width_mi + X);
}

	

template<class _DataType>
void cSkeleton<_DataType>::SaveInitVolume(_DataType Min, _DataType Max)
{
	char	Postfix[256];
	
	sprintf (Postfix, "InitDF_%03d_%03d", (int)Min, (int)Max);
	SaveVolume(Distance_mi, (float)0.0, (float)255.0, Postfix);
	printf ("The init binary volume is saved, Postfix = %s\n", Postfix);
	fflush(stdout);
}

template<class _DataType>
void cSkeleton<_DataType>::SaveDistanceVolume()
{
	int		i, Mini, Maxi;
	
	
	Mini = 9999999;
	Maxi = -9999999;
	
	for (i=0; i<WHD_mi; i++) {
		if (Mini > Distance_mi[i]) Mini = Distance_mi[i];
		if (Maxi < Distance_mi[i]) Maxi = Distance_mi[i];
	}
	printf ("Min & Max Distance = (%d, %d)\n", Mini, Maxi);
	fflush(stdout);
	
	SaveVolume(Distance_mi, (float)Mini, (float)Maxi, "Dist");
	

/*
	//-----------------------------------------------------
	// For Debugging
	unsigned char	*CovertedDist = new unsigned char[WHD_mi];
	
	for (i=0; i<WHD_mi; i++) {
		if (Distance_mi[i]>=1 && Distance_mi[i]<=4) CovertedDist[i] = (unsigned char)(255-Distance_mi[i]);
		else CovertedDist[i] = (unsigned char)0;
	}
	SaveVolume(CovertedDist, (float)Mini, (float)Maxi, "Dist_Converted");
	delete [] CovertedDist;
*/	
	
}



template<class _DataType>
void cSkeleton<_DataType>::DisplayFlag(int Flag)
{
	switch (Flag) {
		case FLAG_CONNECTED : printf ("CON"); break;
		case FLAG_LOCAL_MAX	: printf ("MAX"); break;
		case FLAG_LOCAL_MIN	: printf ("MIN"); break;
		case FLAG_SEGMENTED	: printf ("SEG"); break;
		case FLAG_EMPTY     : printf ("EPT"); break;
		default: printf ("ERR"); break;
	}

}


template<class _DataType>
void cSkeleton<_DataType>::Display_Distance(int ZPlane)
{
	int		i, j;
	

	for (i=0; i<108; i++) printf ("Z Plane = %d       ", ZPlane);
	printf ("\n");

	printf ("    ");
	for (i=0; i<Width_mi; i++) printf ("%4d", i);
	printf ("\n");
	
	for (j=0; j<Height_mi; j++) {
		printf ("%4d", j);
		for (i=0; i<Width_mi; i++) {
			printf ("%4d", Distance_mi[Index(i, j, ZPlane)]);
		}
		printf ("\n");
	}
	fflush(stdout);
}

template<class _DataType>
void cSkeleton<_DataType>::Display_Distance2(int ZPlane)
{
	int		i, j;
	

	for (i=0; i<108; i++) printf ("Z Plane = %d       ", ZPlane);
	printf ("\n");

	printf ("    ");
	for (i=0; i<Width_mi; i++) printf ("%4d", i);
	printf ("\n");
	
	for (j=0; j<Height_mi; j++) {
		printf ("%4d", j);
		for (i=0; i<Width_mi; i++) {
			if (i>=110 && i<=175 && j>=211 && j<=280) {
				printf ("%4d", Distance_mi[Index(i, j, ZPlane)]);
			}
			else printf ("    ");
		}
		printf ("\n");
	}
	fflush(stdout);
}



template<class _DataType>
void cSkeleton<_DataType>::Destroy()
{
	delete [] InitVolume_muc;
	InitVolume_muc = NULL;

	delete [] Distance_mi;
	Distance_mi = NULL;

	delete [] CCVolume_muc;
	CCVolume_muc = NULL;

	delete [] GVFDistance_mf;
	GVFDistance_mf = NULL;

	delete [] DistanceFromSkeletons_mi;
	DistanceFromSkeletons_mi = NULL;	

	delete [] Skeletons_muc;
	Skeletons_muc = NULL;
}


cSkeleton<unsigned char>	__SkeletonValue0;
//cSkeleton<unsigned short>	__SkeletonValue1;
//cSkeleton<int>				__SkeletonValue2;
//cSkeleton<float>			__SkeletonValue3;


