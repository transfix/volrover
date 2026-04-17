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

#include <PEDetection/ConnV_AmbiF_Table.h>
#include <PEDetection/TFGeneration.h>
#include <PEDetection/CompileOptions.h>
#include <PEDetection/Stack.h>

#define TRUE	1
#define FALSE	0

// Neighbor Index Table for the Connectivity
int	NeighborIdxTable[27][2] = {
	{-1, -1,},		// 0
	{4, 10,},		// 1 
	{-1, -1,},		// 2
	{4, 12,},		// 3 
	{-1, -1,},		// 4
	{4, 14,},		// 5
	{-1, -1,},		// 6
	{4, 16,}, 		// 7
	{-1, -1,},		// 8
	{10, 12,},		// 9
	{-1, -1,},		// 10
	{10, 14,},	 	// 11
	{-1, -1,},		// 12
	{-1, -1,},		// 13
	{-1, -1,},		// 14
	{12, 16,},		// 15
	{-1, -1,},		// 16
	{14, 16,},		// 17
	{-1, -1,},		// 18
	{10, 22,},		// 19
	{-1, -1,},		// 20
	{12, 22,},		// 21
	{-1, -1,},		// 22
	{14, 22,},		// 23
	{-1, -1,},		// 24
	{16, 22,},		// 25
	{-1, -1,},		// 26
};


template<class T>
T MIN(T X, T Y) { if (X <= Y) { return X; } else { return Y;} }


template<class T>
void SaveVolume(T *data, float Minf, float Maxf, char *Name);

//----------------------------------------------------------------------------
// cTFGeneration Class Member Functions
//----------------------------------------------------------------------------

template <class _DataType>
cTFGeneration<_DataType>::cTFGeneration()
{
	H_vg_mf = NULL;
	VoxelStatus_muc = NULL;
	SecondDerivative_mf = NULL;
	VoxelVolume_muc = NULL;
	
	// Kindlmann's sigma computation
	H_vg_md = NULL;
	G_v_md = NULL;
	
	// Connected Surface Computation
	CCIndexTable_mi = NULL;

	// Saving the distance from a voxel to a hit location
	DistanceToZeroCrossingLoc_mf = NULL;
	
	HistogramVolume_mi = NULL;
}


// destructor
template <class _DataType>
cTFGeneration<_DataType>::~cTFGeneration()
{
	delete [] VoxelVolume_muc;
	delete [] H_vg_md;
	delete [] G_v_md;
	delete [] CCIndexTable_mi;
	delete [] DistanceToZeroCrossingLoc_mf;
	delete [] HistogramVolume_mi;
}


template <class _DataType>
void cTFGeneration<_DataType>::setData(_DataType *Data, float Minf, float Maxf)
{
	MinData_mf = Minf;
	MaxData_mf = Maxf;
	Data_mT = Data;
}

template <class _DataType>
void cTFGeneration<_DataType>::setGradient(float *Grad, float Minf, float Maxf)
{
	MinGrad_mf = Minf;
	MaxGrad_mf = Maxf;
	GradientMag_mf = Grad;
}

template<class _DataType>
void cTFGeneration<_DataType>::setGradientVectors(float *GradVec)
{
	GradientVec_mf = GradVec;
}

template<class _DataType>
void cTFGeneration<_DataType>::setWHD(int W, int H, int D)
{
	Width_mi = W;
	Height_mi = H;
	Depth_mi = D;
	
	WtimesH_mi = W*H;
	WHD_mi = W*H*D;
	
	VoxelStatus_muc = new unsigned char[WHD_mi];
	AlphaVolume_mf = new float[WHD_mi];
}

template<class _DataType>
void cTFGeneration<_DataType>::setProbabilityHistogram(float *Prob, int NumMaterial, int *Histo, float HistoF)
{
	MaterialProb_mf = Prob;
	Histogram_mi = Histo;
	HistogramFactorI_mf = HistoF;
	HistogramFactorG_mf = HistoF;
	NumMaterial_mi = NumMaterial;
}


template<class _DataType>
void cTFGeneration<_DataType>::ComputeDistanceHistoAt(double Intensity, double GradMag,
									float *DistanceVolume_f, float MinDist, float MaxDist)
{
	int		i, NumI_i, NumIG_i;
	int		*Histogram = new int[(int)(MaxDist)+1];
	
	
	printf ("MinDist = %f, MaxDist = %f\n", MinDist, MaxDist);
	for (i=0; i<(int)(MaxDist)+1; i++) {
		printf ("%3d, ", i);
	}
	printf ("\n");

	for (float GM=GradMag-90.0; GM<GradMag+90.0; GM+=1.0) {
	
		for (i=0; i<(int)(MaxDist)+1; i++) Histogram[i]=0;

		NumI_i = 0;
		NumIG_i = 0;
		for (i=0; i<WHD_mi; i++) {
			if (fabs(Intensity - Data_mT[i])<1e-5) NumI_i++;
			if (fabs(Intensity - Data_mT[i])<1e-5 && fabs(GM - GradientMag_mf[i])<0.5) NumIG_i++;

			if (fabs(Intensity - Data_mT[i])<1e-5 && fabs(GM - GradientMag_mf[i])<0.5) {
				Histogram[(int)(DistanceVolume_f[i])]++;
			}
		}

		printf ("Histogram of the Distance Volume at Intensity = %f, Gradient = %f\n", Intensity, GM);
		for (i=0; i<(int)(MaxDist)+1; i++) {
			printf ("%3d, ", Histogram[i]);
		}
		printf ("\n\n");
		printf ("The number of voxels at the Intensity & GradM = %d / %d (Only intensity)\n", NumIG_i, NumI_i);
	}
	printf ("\n");
	fflush(stdout);
		


}


template<class _DataType>
void cTFGeneration<_DataType>::ComputeDistanceVolume(double GradThreshold, char *TargetFileName)
{
	int		i, k, DataCoor[3], FoundZeroCrossingLoc_i;
	double	GradVec_d[3], CurrDataLoc_d[3], ZeroCrossingLoc_d[3], LocalMaxGradient_d;
	double	DataPosFromZeroCrossingLoc_d, Tempd;
	int		CurrZ = 0, Distance_i;
	unsigned char	*DistanceVolume_uc = new unsigned char [WHD_mi];
	
	
	printf ("TargetFileName = %s\n", TargetFileName);
	for (i=0; i<WHD_mi; i++) {
	
		if (GradientMag_mf[i] < GradThreshold) {
			DistanceVolume_uc[i] = (unsigned char)255;
			continue;
		}

		DataCoor[2] = i/WtimesH_mi;
		DataCoor[1] = (i - DataCoor[2]*WtimesH_mi)/Width_mi;
		DataCoor[0] = i % Width_mi;

		if (CurrZ<=DataCoor[2]) {
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
		if (fabs(Tempd)<1e-5) {
			DistanceVolume_uc[i] = (unsigned char)255;
			continue; // To skip zero-length vectors
		}
		for (k=0; k<3; k++) GradVec_d[k] /= Tempd; // Normalize the gradient vector
		for (k=0; k<3; k++) CurrDataLoc_d[k] = (double)DataCoor[k];
		FoundZeroCrossingLoc_i = FindZeroCrossingLocation(CurrDataLoc_d, GradVec_d, ZeroCrossingLoc_d, 
														LocalMaxGradient_d, DataPosFromZeroCrossingLoc_d);
		if (FoundZeroCrossingLoc_i) {
			// Riversed Distance 0(opaque) - 254(transparent)
			Distance_i = (int)(fabs(DataPosFromZeroCrossingLoc_d)/15.0*255.0);
			DistanceVolume_uc[i] = (unsigned char)Distance_i;
		}
		else DistanceVolume_uc[i] = (unsigned char)255;

	}

	printf ("Computing Distance is done\n");
	printf ("\n");
	fflush(stdout);

	SaveVolume(DistanceVolume_uc, (float)0.0, (float)255.0, "Distance");
	printf ("Distance volume is saved\n"); fflush (stdout);
	
	delete [] DistanceVolume_uc;

}

// Computing Distance_uc volume 
//
// Drawing Histograms
//		IG_General, IGAtHit, IDAtHit Graphs Computation
//
// Volume Saving
//		DistanceVolume_f, GMVolume_uc Saving
//
template<class _DataType>
float* cTFGeneration<_DataType>::ComputeDistanceVolume2(double GradThreshold, float& Minf_ret, float& Maxf_ret,
												 char *TargetFileName, int NumBins)
{
	int		i, j, k, loc[3], DataCoor[3], MaxAxisResolution_i, FoundZeroCrossingLoc_i;
	double	GradVec_d[3], CurrDataLoc_d[3], ZeroCrossingLoc_d[3], LocalMaxGradient_d;
	double	DataPosFromZeroCrossingLoc_d, Tempd;
	float	ComputedMaxDistance_f;
	int		l, m, n;
	
	float	*DistanceVolume_f = new float [WHD_mi];
	unsigned char	*DistanceVolume_uc = new unsigned char [WHD_mi];
	unsigned char	*GMVolume_uc = new unsigned char [WHD_mi];
	


	//----------------------------------------------------------------------------
	// Initializing the I-G histogram
	int		*HistogramGI; // Histogram of (Intensity, Gradient Magnitude)
	int		Intensity, GradMag;
	HistogramGI = new int [(NumBins+1)*(NumBins+1)];
	for(i=0; i<(NumBins+1)*(NumBins+1); i++) HistogramGI[i] = 0;
	//----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	// Initializing the I-D histogram
	int		*HistogramDI; // Histogram of (Intensity, Gradient Magnitude)
	int		Distance_i;
	HistogramDI = new int [(NumBins+1)*(NumBins+1)*3];
	for(i=0; i<(NumBins+1)*(NumBins+1)*3; i++) HistogramDI[i] = 0;
	//----------------------------------------------------------------------------

	int		CurrZ = 0;
	float	MinData_f, MaxData_f, MinGM_f, MaxGM_f, MinSD_f, MaxSD_f;

	printf ("Min & Max Value Computation\n");
	fflush (stdout);
	// Computing distance and save it into DistanceToZeroCrossingLoc_mf[]
	ComputeAndSaveDistanceFloat(GradThreshold, MinData_f, MaxData_f, MinGM_f, MaxGM_f, 
								MinSD_f, MaxSD_f, TargetFileName);

	MinGM_f = MinGrad_mf;
//	MaxGM_f = MaxGrad_mf + 10.0;
	MaxGM_f = MaxGrad_mf;
	ComputeIGGraphGeneral("IG", 256, MinGM_f, MaxGM_f);
	
	Minf_ret = FLT_MAX;
	Maxf_ret = -FLT_MAX;
	ComputedMaxDistance_f = -FLT_MAX;
	
	MaxAxisResolution_i = Width_mi;
	if (MaxAxisResolution_i < Height_mi) MaxAxisResolution_i = Height_mi;
	if (MaxAxisResolution_i < Depth_mi) MaxAxisResolution_i = Depth_mi;

	ZeroCrossingCells_mf = new float[WHD_mi];
	for (i=0; i<WHD_mi; i++) ZeroCrossingCells_mf[i] = MinSecMag_mf - 10.0;
	

#ifdef	ZERO_CELLS_GRADMAG							
	//----------------------------------------------------------------------------
	// To save zero crossing cells with Gradient Magnitude values
	unsigned char *ZeroCrossingCells_uc = new unsigned char[WHD_mi];
	for (i=0; i<WHD_mi; i++) ZeroCrossingCells_uc[i] = (unsigned char)0;

	// Initializing the I-G histogram at Zero Cells
	int		*Histogram_ZeroCell_GI; // Histogram of (Intensity, Gradient Magnitude)
	Histogram_ZeroCell_GI = new int [(NumBins+1)*(NumBins+1)];
	for(i=0; i<(NumBins+1)*(NumBins+1); i++) Histogram_ZeroCell_GI[i] = 0;
	//----------------------------------------------------------------------------
#endif


	for (i=0; i<WHD_mi; i++) {
	
		if (GradientMag_mf[i] < GradThreshold) {
			DistanceVolume_f[i] = (float)MaxAxisResolution_i;
			DistanceVolume_uc[i] = (unsigned char)255;
			continue;
		}

		DataCoor[2] = i/WtimesH_mi;
		DataCoor[1] = (i - DataCoor[2]*WtimesH_mi)/Width_mi;
		DataCoor[0] = i % Width_mi;

		if (CurrZ<=DataCoor[2]) {
			printf ("Z = %3d, ", DataCoor[2]);
			printf ("Min & Max Distance = %7.2f %7.2f, ", Minf_ret, Maxf_ret);
			printf ("Computed Max Distance = %.5f, ", ComputedMaxDistance_f);
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
			DistanceVolume_f[i] = (float)MaxAxisResolution_i;
			DistanceVolume_uc[i] = (unsigned char)255;
			continue; // To skip zero-length vectors
		}
		for (k=0; k<3; k++) GradVec_d[k] /= Tempd; // Normalize the gradient vector
		for (k=0; k<3; k++) CurrDataLoc_d[k] = (double)DataCoor[k];
//		FoundZeroCrossingLoc_i = FindZeroCrossingLocation(CurrDataLoc_d, GradVec_d, ZeroCrossingLoc_d, 
//														LocalMaxGradient_d, DataPosFromZeroCrossingLoc_d);

		for (k=0; k<3; k++) ZeroCrossingLoc_d[k] = CurrDataLoc_d[k] + GradVec_d[k]*DistanceToZeroCrossingLoc_mf[i];
		LocalMaxGradient_d = GradientInterpolation(ZeroCrossingLoc_d);
		DataPosFromZeroCrossingLoc_d = (double)DistanceToZeroCrossingLoc_mf[i];
		if (DataPosFromZeroCrossingLoc_d>=MaxAxisResolution_i) FoundZeroCrossingLoc_i = false;
		else FoundZeroCrossingLoc_i = true;
		
		if (FoundZeroCrossingLoc_i) {
			if (LocalMaxGradient_d >= GradThreshold) {
				DistanceVolume_f[i] = (float)DataPosFromZeroCrossingLoc_d;

				for (n=(int)ZeroCrossingLoc_d[2]; n<=(int)ZeroCrossingLoc_d[2]+1; n++) {
					for (m=(int)ZeroCrossingLoc_d[1]; m<=(int)ZeroCrossingLoc_d[1]+1; m++) {
						for (l=(int)ZeroCrossingLoc_d[0]; l<=(int)ZeroCrossingLoc_d[0]+1; l++) {
							if (l>=Width_mi || m>=Height_mi || n>=Depth_mi) continue;
							loc[0] = n*WtimesH_mi + m*Width_mi + l;

							ZeroCrossingCells_mf[loc[0]] = SecondDerivative_mf[loc[0]];
							
#ifdef	ZERO_CELLS_GRADMAG							
							// Zero Cells with Gradient Magnitudes
							//----------------------------------------------------
							// Making a volume dataset
							double	Tempd = GradientMag_mf[loc[0]]/100.0*255.0;
							int		Tempi;
							Tempi = (int)Tempd;
							if (Tempd<=0) Tempi = 0;
							if (Tempd>=255) Tempi = 255;
							ZeroCrossingCells_uc[loc[0]] = (unsigned char)Tempi;
							
						#ifdef	TWO_DIM_R3_IMAGE
							if (n==1) 
						#else
							if (1) 
						#endif
							{
								// Computing a Histogram At the Hit Locations
								Intensity = (int)(DataInterpolation(ZeroCrossingLoc_d)/255.0*NumBins);
//								Tempd = (LocalMaxGradient_d-MinGrad_mf)/(MaxGrad_mf-MinGrad_mf)*NumBins;
								Tempd = (LocalMaxGradient_d-MinGM_f)/(MaxGM_f-MinGM_f)*NumBins;
								GradMag = (int)(Tempd);

								if (Intensity<0) Intensity = 0;
								if (Intensity>NumBins) Intensity = NumBins;
								if (GradMag<0) GradMag = 0;
								if (GradMag>NumBins) GradMag = NumBins;
								Histogram_ZeroCell_GI[GradMag*(NumBins+1) + Intensity]++;
							}
							//-----------------------------------------------------
#endif

						}
					}
				}

			}
			else DistanceVolume_f[i] = (float)MaxAxisResolution_i;

			if (DistanceVolume_f[i] > (float)MaxAxisResolution_i) {
				DistanceVolume_f[i] = (float)MaxAxisResolution_i;
			}
		}
		else DistanceVolume_f[i] = (float)MaxAxisResolution_i;
		
		if (Minf_ret > DistanceVolume_f[i]) Minf_ret = DistanceVolume_f[i];
		if (Maxf_ret < DistanceVolume_f[i]) Maxf_ret = DistanceVolume_f[i];
		
		if (ComputedMaxDistance_f < DistanceVolume_f[i] && 
			DistanceVolume_f[i] <= (float)MaxAxisResolution_i - 1e-4) {
			ComputedMaxDistance_f = DistanceVolume_f[i];
		}

		if (FoundZeroCrossingLoc_i) {
//			Intensity = (int)(DataInterpolation(ZeroCrossingLoc_d)/255.0*NumBins);
			Intensity = (int)(Data_mT[i]);
//			Tempd = (LocalMaxGradient_d-MinGrad_mf)/(MaxGrad_mf-MinGrad_mf)*NumBins;	// Orginal Min Max
			Tempd = (LocalMaxGradient_d-MinGM_f)/(MaxGM_f-MinGM_f)*NumBins;	// New Min Max
			GradMag = (int)(Tempd);
			
			if (Intensity<0) Intensity = 0;
			if (Intensity>NumBins) Intensity = NumBins;
			if (GradMag<0) GradMag = 0;
			if (GradMag>NumBins) GradMag = NumBins;
			
			#ifdef	TWO_DIM_R3_IMAGE
				if (DataCoor[2]==1)	HistogramGI[GradMag*(NumBins+1) + Intensity]++;
			#else
				HistogramGI[GradMag*(NumBins+1) + Intensity]++;
			#endif
			
			// Gradient Magnitudes at Hit Locations
			GMVolume_uc[i] = (unsigned char)GradMag;

			// Distance 0 is opaque & 255 is transparent
			Distance_i = (int)(fabs(DataPosFromZeroCrossingLoc_d)/15.0*255.0);
			if (Distance_i > 255) Distance_i = (unsigned char)255;
			DistanceVolume_uc[i] = (unsigned char)Distance_i;


			// Histogram
//			Tempd = (LocalMaxGradient_d-MinGrad_mf)/(MaxGrad_mf-MinGrad_mf)*(NumBins-15) + 15.0; // 15 - 240
			Tempd = (LocalMaxGradient_d-MinGM_f)/(MaxGM_f-MinGM_f)*(NumBins-15) + 15.0; // 15 - 240
			Distance_i = (int)(Tempd + DataPosFromZeroCrossingLoc_d); // 0 - 255
			if (Distance_i<0) Distance_i = 0;
			if (Distance_i>NumBins) Distance_i = NumBins;
			#ifdef	TWO_DIM_R3_IMAGE
				if (DataPosFromZeroCrossingLoc_d<0) {
					if (DataCoor[2]==1)	HistogramDI[Distance_i*(NumBins+1)*3 + Intensity*3]++; // Red
				}
				else if (DataCoor[2]==1) HistogramDI[Distance_i*(NumBins+1)*3 + Intensity*3 + 2]++; // Blue
			#else
				if (DataPosFromZeroCrossingLoc_d<0) HistogramDI[Distance_i*(NumBins+1)*3 + Intensity*3]++; // Red
				else HistogramDI[Distance_i*(NumBins+1)*3 + Intensity*3 + 2]++; // Blue
			#endif
		}
		else {
			DistanceVolume_uc[i] = (unsigned char)255;
			GMVolume_uc[i] = (unsigned char)0;
		}



	}
	printf ("Min Max Distance = %f %f\n", Minf_ret, Maxf_ret);
	printf ("Computing Distance is done\n");
	printf ("\n");
	fflush(stdout);


	SaveVolume(DistanceVolume_uc, (float)0.0, (float)255.0, "Distance");
	SaveVolume(GMVolume_uc, (float)0.0, (float)255.0, "GMHit");
	delete [] DistanceVolume_uc;
	delete [] GMVolume_uc;

	Maxf_ret = ComputedMaxDistance_f;

	// Display G-I Histogram
	printf ("IG Graph at Hit Locations\n");
	printf ("%d %d\n", NumBins, NumBins);
	for (j=0; j<NumBins; j++) {
		for (i=0; i<NumBins; i++) {
			printf ("%5d ", HistogramGI[j*(NumBins+1) + i]);
		}
		printf ("\n");
	}
	printf ("\n\n");
	fflush(stdout);
	delete [] HistogramGI;
	
	
	// Display D-I Histogram
	printf ("ID Graph at Hit Locations\n");
	printf ("%d %d\n", NumBins, NumBins);
	for (j=0; j<NumBins; j++) {
		for (i=0; i<NumBins; i++) {
			printf ("%5d ", HistogramDI[j*(NumBins+1)*3 + i*3 + 0]);
			printf ("%1d ", HistogramDI[j*(NumBins+1)*3 + i*3 + 1]);
			printf ("%5d ", HistogramDI[j*(NumBins+1)*3 + i*3 + 2]);
		}
		printf ("\n");
	}
	printf ("\n\n");
	fflush(stdout);
	delete [] HistogramDI;

	printf ("Zero Cells \n");
	SaveVolume(ZeroCrossingCells_mf, MinSecMag_mf, MaxSecMag_mf, "ZeroCells");
	fflush (stdout);

#ifdef	ZERO_CELLS_GRADMAG							

	// Display G-I Histogram
	printf ("IG Graph at four surrounded voxels at the Hit Locations\n");
	printf ("%d %d\n", NumBins, NumBins);
	for (j=0; j<NumBins; j++) {
		for (i=0; i<NumBins; i++) {
			printf ("%5d ", Histogram_ZeroCell_GI[j*(NumBins+1) + i]);
		}
		printf ("\n");
	}
	printf ("\n\n");
	fflush(stdout);
	delete [] Histogram_ZeroCell_GI;


	// Zero Cells with Gradient Magnitudes
	SaveVolume(ZeroCrossingCells_uc, MinSecMag_mf, MaxSecMag_mf, "ZeroCells_GradM");
#endif

//	Connected_Positive_SecondD_Voxels_RunSave(false);

	
	
	return DistanceVolume_f;

}


template<class _DataType>
void cTFGeneration<_DataType>::Compute_Dist_GM(double GradThreshold, float& Minf_ret, float& Maxf_ret,
				unsigned char *DistanceVolume_uc, unsigned char *GMVolume_uc, char *TargetFileName)
{
	int		i, k, DataCoor[3], MaxAxisResolution_i, FoundZeroCrossingLoc_i;
	double	GradVec_d[3], CurrDataLoc_d[3], ZeroCrossingLoc_d[3], LocalMaxGradient_d;
	double	DataPosFromZeroCrossingLoc_d, Tempd;
	float	ComputedMaxDistance_f;	
	float	*DistanceVolume_f = new float [WHD_mi];
//	unsigned char	*DistanceVolume_uc = new unsigned char [WHD_mi];
//	unsigned char	*GMVolume_uc = new unsigned char [WHD_mi];
	
	int		CurrZ = 0;
	float	MinData_f, MaxData_f, MinGM_f, MaxGM_f, MinSD_f, MaxSD_f;
	int		Intensity_i, GradMag_i, Distance_i;

	// Computing distance and save it into DistanceToZeroCrossingLoc_mf[]
	ComputeAndSaveDistanceFloat(GradThreshold, MinData_f, MaxData_f, MinGM_f, MaxGM_f, 
								MinSD_f, MaxSD_f, TargetFileName);

	MinGM_f = MinGrad_mf;
	MaxGM_f = MaxGrad_mf;
	
	Minf_ret = FLT_MAX;
	Maxf_ret = -FLT_MAX;
	ComputedMaxDistance_f = -FLT_MAX;
	
	MaxAxisResolution_i = Width_mi;
	if (MaxAxisResolution_i < Height_mi) MaxAxisResolution_i = Height_mi;
	if (MaxAxisResolution_i < Depth_mi) MaxAxisResolution_i = Depth_mi;

#ifdef TF_GENERATION_COMPUTE_DOT_PRODUCT
	double	Threshold_GM_CurrLoc = 5.0;
	double	Threshold_GM_BoundaryLoc = 20.0;
	double	GradVecAtHitLoc_d[3];
	double	Total_Dot_Product_d = 0.0, DotP_d;
	double	MinDotP_d = 10.0, MaxDotP_d = -10.0;
	double	Total_GradMag_d = 0.0;
	int		TotalNumVoxels = 0, NumNegVoxels = 0, NumPosVoxels = 0;
	int		NumVoxels[100];
	double	GradMag_d[100];
	double	SumDistance_d[100];
	double	DotProduct_d[100];
	for (i=0; i<100; i++) {
		NumVoxels[i] = 0;
		GradMag_d[i] = 0;
		SumDistance_d[i] = 0;
		DotProduct_d[i] = 0;
	}
#endif

	for (i=0; i<WHD_mi; i++) {
	
		if (GradientMag_mf[i] < GradThreshold) {
			DistanceVolume_f[i] = (float)MaxAxisResolution_i;
			DistanceVolume_uc[i] = (unsigned char)255;
			continue;
		}

		DataCoor[2] = i/WtimesH_mi;
		DataCoor[1] = (i - DataCoor[2]*WtimesH_mi)/Width_mi;
		DataCoor[0] = i % Width_mi;

		if (CurrZ<=DataCoor[2]) {
			printf ("Z = %3d, ", DataCoor[2]);
			printf ("Min & Max Distance = %7.2f %7.2f, ", Minf_ret, Maxf_ret);
			printf ("Computed Max Distance = %.5f, ", ComputedMaxDistance_f);
			printf ("\n");
			fflush(stdout);
			CurrZ++;
		}

		//-------------------------------------------------------------------------------------------
		// Finding the local maxima of gradient magnitudes along the gradient direction. 
		// It climbs the mountain of gradient magnitudes to find the zero-crossing second derivative 
		// Return ZeroCrossingLoc_d, LocalMaxGradient, DataPosFromZeroCrossingLoc_d
		// GradVec_d[] = the gradient vector at the current position
		//-------------------------------------------------------------------------------------------
		// Getting the gradient vector at the position
		for (k=0; k<3; k++) GradVec_d[k] = (double)GradientVec_mf[i*3 + k];
		Tempd = sqrt (GradVec_d[0]*GradVec_d[0] + GradVec_d[1]*GradVec_d[1] + GradVec_d[2]*GradVec_d[2]);
		if (fabs(Tempd)<1e-6) {
			DistanceVolume_f[i] = (float)MaxAxisResolution_i;
			DistanceVolume_uc[i] = (unsigned char)255;
			continue; // To skip zero-length vectors
		}
		for (k=0; k<3; k++) GradVec_d[k] /= Tempd; // Normalize the gradient vector
		for (k=0; k<3; k++) CurrDataLoc_d[k] = (double)DataCoor[k];
		for (k=0; k<3; k++) ZeroCrossingLoc_d[k] = CurrDataLoc_d[k] + GradVec_d[k]*DistanceToZeroCrossingLoc_mf[i];
		LocalMaxGradient_d = GradientInterpolation(ZeroCrossingLoc_d);
		DataPosFromZeroCrossingLoc_d = (double)DistanceToZeroCrossingLoc_mf[i];
		if (DataPosFromZeroCrossingLoc_d>=MaxAxisResolution_i) FoundZeroCrossingLoc_i = false;
		else FoundZeroCrossingLoc_i = true;
		
		if (FoundZeroCrossingLoc_i) {
			if (LocalMaxGradient_d >= GradThreshold) {
				DistanceVolume_f[i] = (float)DataPosFromZeroCrossingLoc_d;
			}
			else DistanceVolume_f[i] = (float)MaxAxisResolution_i;

			if (DistanceVolume_f[i] > (float)MaxAxisResolution_i) {
				DistanceVolume_f[i] = (float)MaxAxisResolution_i;
			}
		}
		else DistanceVolume_f[i] = (float)MaxAxisResolution_i;
		
		if (Minf_ret > DistanceVolume_f[i]) Minf_ret = DistanceVolume_f[i];
		if (Maxf_ret < DistanceVolume_f[i]) Maxf_ret = DistanceVolume_f[i];
		
		if (ComputedMaxDistance_f < DistanceVolume_f[i] && 
			DistanceVolume_f[i] <= (float)MaxAxisResolution_i - 1e-4) {
			ComputedMaxDistance_f = DistanceVolume_f[i];
		}

		if (FoundZeroCrossingLoc_i) {
			Intensity_i = (int)(Data_mT[i]);
			Tempd = (LocalMaxGradient_d-MinGM_f)/(MaxGM_f-MinGM_f)*255;	// New Min Max
			GradMag_i = (int)(Tempd);
			
			if (Intensity_i<0) Intensity_i = 0;
			if (Intensity_i>255) Intensity_i = 255;
			if (GradMag_i<0) GradMag_i = 0;
			if (GradMag_i>255) GradMag_i = 255;
						
			// Gradient Magnitudes at Hit Locations
			GMVolume_uc[i] = (unsigned char)GradMag_i;

			// Distance 0 is opaque & 255 is transparent
			Distance_i = (int)(fabs(DataPosFromZeroCrossingLoc_d)/15.0*255.0);
			if (Distance_i > 255) Distance_i = (unsigned char)255;
			DistanceVolume_uc[i] = (unsigned char)Distance_i;

#ifdef TF_GENERATION_COMPUTE_DOT_PRODUCT
			GradVecInterpolation(ZeroCrossingLoc_d, GradVecAtHitLoc_d);
			NormalizeVector(GradVecAtHitLoc_d);
			DotP_d = GradVec_d[0]*GradVecAtHitLoc_d[0] + GradVec_d[1]*GradVecAtHitLoc_d[1] + 
						GradVec_d[2]*GradVecAtHitLoc_d[2];
			if (GradientMag_mf[i]>Threshold_GM_CurrLoc && 
				LocalMaxGradient_d>Threshold_GM_BoundaryLoc && DotP_d >= 0) {
				Total_GradMag_d += LocalMaxGradient_d;
				Total_Dot_Product_d += DotP_d;
				if (MinDotP_d > DotP_d) MinDotP_d = DotP_d;
				if (MaxDotP_d < DotP_d) MaxDotP_d = DotP_d;
				if (DotP_d > 0) NumPosVoxels++;
				if (DotP_d < 0) NumNegVoxels++;

				int	DotP_i = (int)(DotP_d*100.0);
				if (DotP_i>=100) DotP_i = 99;
				NumVoxels[DotP_i]++;
				GradMag_d[DotP_i]+=LocalMaxGradient_d;

				SumDistance_d[DotP_i] += fabs(DataPosFromZeroCrossingLoc_d);

				DotProduct_d[DotP_i] += DotP_d;

				TotalNumVoxels++;
			}
#endif
		}
		else {
			DistanceVolume_uc[i] = (unsigned char)255;
			GMVolume_uc[i] = (unsigned char)0;
		}
	}

#ifdef TF_GENERATION_COMPUTE_DOT_PRODUCT
	printf ("\n");
	printf ("Threshold GM for the voxel current location = %f\n", Threshold_GM_CurrLoc);
	printf ("Threshold GM for the boundary location = %f\n", Threshold_GM_BoundaryLoc);
	printf ("Num Voxels = %d/%d\n", TotalNumVoxels, WHD_mi);
	printf ("Dot Product = %f\n", Total_Dot_Product_d);
	printf ("Ave. Dot Product = %f\n", Total_Dot_Product_d/TotalNumVoxels);
	printf ("Ave. Grad. Mag = %f\n", Total_GradMag_d/TotalNumVoxels);
	printf ("Min & Max Dot Product = %f, %f\n", MinDotP_d, MaxDotP_d);
	printf ("Num Negative Voxels = %d, Num Positive Voxels = %d\n", NumNegVoxels, NumPosVoxels);
	
	double	AccumulatedPercent_d[100];
	for (i=0; i<100; i++) AccumulatedPercent_d[i] = 0.0;
	AccumulatedPercent_d[99] = 100.0*NumVoxels[99]/TotalNumVoxels;
	for (i=98; i>=0; i--) {
		AccumulatedPercent_d[i] = 100.0*NumVoxels[i]/TotalNumVoxels;
		AccumulatedPercent_d[i] += AccumulatedPercent_d[i+1];
	}
	for (i=0; i<100; i++) {
		printf ("NumVoxels[%3d] = %8d, ", i, NumVoxels[i]);
		if (NumVoxels[i]>0) printf ("Ave. DotP = %6.4f, ", 1.0*DotProduct_d[i]/NumVoxels[i]);
		else printf ("Ave. DotP =      0 ");
		if (NumVoxels[i]>0) printf ("GradMag = %8.4f ", GradMag_d[i]/NumVoxels[i]);
		else printf ("GradMag =        0 ");
		printf ("Percent = %6.2f, ", 100.0*NumVoxels[i]/TotalNumVoxels);
		printf ("Acc. Percent = %6.2f, ", AccumulatedPercent_d[i]);
		if (NumVoxels[i]>0)printf ("Ave. Distance = %f, ", SumDistance_d[i]/NumVoxels[i]);
		else printf ("Ave. Distance = 0, ");
		printf ("\n"); fflush (stdout);
	}
#endif
	
	printf ("Computing Distance is done\n");
	printf ("\n");
	fflush(stdout);

	Maxf_ret = ComputedMaxDistance_f;

}


template<class _DataType>
void cTFGeneration<_DataType>::ComputeConnectedSurfaceVolume2(double GradThreshold, char *TargetFileName)
{
	int		i;
	char	LocalMaxGM_FileName[512];
	int		LocalMaxGM_fd;
	float	MinGM_f, MaxGM_f;
	
	
	sprintf (LocalMaxGM_FileName, "%s_LocalMaxGMVolume.rawiv", TargetFileName);
	printf ("LocalMaxGM FileName = %s\n", LocalMaxGM_FileName);
	fflush (stdout);

	if ((LocalMaxGM_fd = open (LocalMaxGM_FileName, O_RDONLY)) < 0) {

		printf ("The file, %s, does not exist\n", LocalMaxGM_FileName);
		fflush (stdout);

		// Making the volume, CSurfaceVolume_mi[], which contains
		// local maximum voxels
		MarkingLocalMaxVoxels(GradThreshold);

		unsigned char	*TempVolume_uc = new unsigned char [WHD_mi];
		for (i=0; i<WHD_mi; i++) TempVolume_uc[i] = (unsigned char)CCSurfaceVolume_mi[i];
		SaveVolume(TempVolume_uc, 0.0, 255.0, "LocalMaxGMVolume");
		delete [] TempVolume_uc;
		
	}
	else {

		printf ("Reading the file, %s\n", LocalMaxGM_FileName);
		fflush (stdout);


		unsigned char	TempHeader[68];
		if (read(LocalMaxGM_fd, TempHeader, 68) != (unsigned int)68) {
			cout << "The file could not be read " << LocalMaxGM_FileName << endl;
			close (LocalMaxGM_fd);
			exit(1);
		}
		unsigned char	*LocalMaxVolume_uc = new unsigned char [WHD_mi];
		if (read(LocalMaxGM_fd, LocalMaxVolume_uc, WHD_mi) != (unsigned int)WHD_mi) {
			cout << "The file could not be read " << LocalMaxGM_FileName << endl;
			close (LocalMaxGM_fd);
			exit(1);
		}
		MinGM_f = 99999.0;
		MaxGM_f = -99999.0;
		CCSurfaceVolume_mi = new int [WHD_mi];
		for (i=0; i<WHD_mi; i++) {
			if (LocalMaxVolume_uc[i]==255) CCSurfaceVolume_mi[i] = 0; // Positive Voxels
			else CCSurfaceVolume_mi[i] = -1; // Negative & Empty Voxels

			if (LocalMaxVolume_uc[i]>0) {
				if (MinGM_f > GradientMag_mf[i]) MinGM_f = GradientMag_mf[i];
				if (MaxGM_f < GradientMag_mf[i]) MaxGM_f = GradientMag_mf[i];
			}
		}
		delete [] LocalMaxVolume_uc;

		printf ("Min & Max Gradient Magnitudes at the hit boundaries = %f, %f\n", MinGM_f, MaxGM_f);
		fflush (stdout);

	}
	
	unsigned char *SurfaceIndexVolume_uc;
	SurfaceIndexVolume_uc = ConnectedPositiveVoxels();
	SaveVolume(SurfaceIndexVolume_uc, 0.0, 255.0, "SurfaceIndexVolume");
	delete [] SurfaceIndexVolume_uc;
	

}


template<class _DataType>
void cTFGeneration<_DataType>::MarkingLocalMaxVoxels(double GradThreshold)
{
	int		i, k, loc[3], DataCoor_i[3], FoundZeroCrossingLoc_i;
	double	GradVec_d[3], CurrDataLoc_d[3], ZeroCrossingLoc_d[3], LocalMaxGradient_d;
	double	DataPosFromZeroCrossingLoc_d, Tempd;
	int		l, m, n, CurrZ = 0;
	int		ZeroLoc_i[3];
	float	MinGM_f, MaxGM_f;
	

	int	Empty_i = 0;
	int	Negative_i = 50;
	int	Positive_i = 255;
	
	CCSurfaceVolume_mi = new int [WHD_mi];
	for (i=0; i<WHD_mi; i++) CCSurfaceVolume_mi[i] = Empty_i;


	MinGM_f = 99999.0;
	MaxGM_f = -99999.0;
	
	for (i=0; i<WHD_mi; i++) {
	
		if (GradientMag_mf[i] < GradThreshold) continue;
		if (CCSurfaceVolume_mi[i]!=Empty_i) continue;
		
		DataCoor_i[2] = i/WtimesH_mi;
		DataCoor_i[1] = (i - DataCoor_i[2]*WtimesH_mi)/Width_mi;
		DataCoor_i[0] = i % Width_mi;

		if (CurrZ==DataCoor_i[2]) {
			printf ("Z = %3d, ", DataCoor_i[2]);
			printf ("\n"); fflush(stdout);
			CurrZ++;
		}


		//-------------------------------------------------------------------------------------------
		// Finding the local maxima of gradient magnitudes along the gradient direction. 
		// It climbs the mountain of gradient magnitudes to find the zero-crossing second derivative 
		// It returns ZeroCrossingLoc_d, LocalMaxGradient, DataPosFromZeroCrossingLoc_d
		//-------------------------------------------------------------------------------------------
		// Getting the gradient vector at the position
		for (k=0; k<3; k++) GradVec_d[k] = (double)GradientVec_mf[i*3 + k];
		Tempd = sqrt (GradVec_d[0]*GradVec_d[0] + GradVec_d[1]*GradVec_d[1] + GradVec_d[2]*GradVec_d[2]);
		if (fabs(Tempd)<1e-6) continue; // To skip zero-length vectors
		for (k=0; k<3; k++) GradVec_d[k] /= Tempd; // Normalizing the gradient vector
		for (k=0; k<3; k++) CurrDataLoc_d[k] = (double)DataCoor_i[k];
		FoundZeroCrossingLoc_i = FindZeroCrossingLocation(CurrDataLoc_d, GradVec_d, ZeroCrossingLoc_d, 
														LocalMaxGradient_d, DataPosFromZeroCrossingLoc_d);
		if (FoundZeroCrossingLoc_i) {
			if (LocalMaxGradient_d >= GradThreshold) {
			
				ZeroLoc_i[0] = (int)floor(ZeroCrossingLoc_d[0]);
				ZeroLoc_i[1] = (int)floor(ZeroCrossingLoc_d[1]);
				ZeroLoc_i[2] = (int)floor(ZeroCrossingLoc_d[2]);
				
				for (n=ZeroLoc_i[2]; n<=ZeroLoc_i[2]+1; n++) {
					for (m=ZeroLoc_i[1]; m<=ZeroLoc_i[1]+1; m++) {
						for (l=ZeroLoc_i[0]; l<=ZeroLoc_i[0]+1; l++) {
							
							if (l>=Width_mi || m>=Height_mi || n>=Depth_mi) continue;
							
							loc[0] = n*WtimesH_mi + m*Width_mi + l;
							if (CCSurfaceVolume_mi[loc[0]]!=Empty_i) continue;
							if (SecondDerivative_mf[loc[0]]>=0) CCSurfaceVolume_mi[loc[0]] = Positive_i;
							if (SecondDerivative_mf[loc[0]]<0) CCSurfaceVolume_mi[loc[0]] = Negative_i;

							if (MinGM_f > GradientMag_mf[loc[0]]) MinGM_f = GradientMag_mf[loc[0]];
							if (MaxGM_f < GradientMag_mf[loc[0]]) MaxGM_f = GradientMag_mf[loc[0]];

						}
					}
				}
				// Treating the cube
			}
		}
		
	}
	
	printf ("Min & Max Gradient Magnitudes at the hit boundaries = %f, %f\n", MinGM_f, MaxGM_f);
	fflush (stdout);
	
}


template<class _DataType>
void cTFGeneration<_DataType>::ComputeAndSaveDistanceFloat(double GradThreshold, 
								float &MinData_Ret, float &MaxData_Ret, 
								float &MinGM_Ret, 	float &MaxGM_Ret, 
								float &MinSD_Ret, float &MaxSD_Ret, char *TargetFileName)
{
	int		i, k;
	float	MinDist_f, MaxDist_f;
	int		DistanceToZero_fd;
	

	int MaxAxisResolution_i = Width_mi;
	if (MaxAxisResolution_i < Height_mi) MaxAxisResolution_i = Height_mi;
	if (MaxAxisResolution_i < Depth_mi) MaxAxisResolution_i = Depth_mi;
	DistanceToZeroCrossingLoc_mf = new float [WHD_mi];
	for (i=0; i<WHD_mi; i++) DistanceToZeroCrossingLoc_mf[i] = MaxAxisResolution_i;
	
	char	DistanceFileName[512];
	sprintf (DistanceFileName, "%s_DistanceToZero.rawiv", TargetFileName);
	printf ("Distance to Zero Crossing Locs FileName = %s\n", DistanceFileName);
	fflush (stdout);

	if ((DistanceToZero_fd = open (DistanceFileName, O_RDONLY)) < 0) {
		printf ("ComputeAndSaveDistanceFloat: The file, %s, does not exist\n", DistanceFileName);
		fflush (stdout);

		// Saving distance to zero-crossing locations into DistanceToZeroCrossingLoc_mf[]
		MinMaxComputationAtHitLocations(GradThreshold, MinData_Ret, MaxData_Ret, 
										MinGM_Ret, MaxGM_Ret, MinSD_Ret, MaxSD_Ret, MinDist_f, MaxDist_f);

		printf ("ComputeAndSaveDistanceFloat: Saving Volume ... \n"); fflush (stdout);

		SaveVolume(DistanceToZeroCrossingLoc_mf, MinDist_f, MaxDist_f, "DistanceToZero");

		printf ("Saving is done\n"); fflush (stdout);

	}
	else {
	
		printf ("Reading the file, %s\n", DistanceFileName);
		fflush (stdout);

		unsigned char	TempHeader_uc[68];
		if (read(DistanceToZero_fd, TempHeader_uc, 68) != (unsigned int)68) {
			cout << "The file could not be read " << DistanceFileName << endl;
			close (DistanceToZero_fd);
			exit(1);
		}
		if (read(DistanceToZero_fd, DistanceToZeroCrossingLoc_mf, sizeof(float)*WHD_mi)
			!= (unsigned int)sizeof(float)*WHD_mi) {
			cout << "The file could not be read " << DistanceFileName << endl;
			close (DistanceToZero_fd);
			exit(1);
		}

		int		l, m, n, loc[3], DataCoor_i[3], ZeroLoc_i[3];
		double	GradVec_d[3], CurrDataLoc_d[3], ZeroCrossingLoc_d[3];
		double	Tempd;

		MinData_Ret = MinGM_Ret = MinSD_Ret = 99999.0;
		MaxData_Ret = MaxGM_Ret = MaxSD_Ret = -99999.0;

		CCSurfaceVolume_mi = new int [WHD_mi];
		for (i=0; i<WHD_mi; i++) {

			DataCoor_i[2] = i/WtimesH_mi;
			DataCoor_i[1] = (i - DataCoor_i[2]*WtimesH_mi)/Width_mi;
			DataCoor_i[0] = i % Width_mi;

			// Getting the gradient vector at the position
			for (k=0; k<3; k++) GradVec_d[k] = (double)GradientVec_mf[i*3 + k];
			Tempd = sqrt (GradVec_d[0]*GradVec_d[0] + GradVec_d[1]*GradVec_d[1] + GradVec_d[2]*GradVec_d[2]);
			if (fabs(Tempd)<1e-6) continue; // To skip zero-length vectors
			for (k=0; k<3; k++) GradVec_d[k] /= Tempd; // Normalizing the gradient vector
			for (k=0; k<3; k++) CurrDataLoc_d[k] = (double)DataCoor_i[k];
			for (k=0; k<3; k++) {
				ZeroCrossingLoc_d[k] = CurrDataLoc_d[k] + GradVec_d[k]*DistanceToZeroCrossingLoc_mf[i];
			}
			for (k=0; k<3; k++) ZeroLoc_i[k] = (int)ZeroCrossingLoc_d[k];
			for (n=ZeroLoc_i[2]; n<=ZeroLoc_i[2]+1; n++) {
				for (m=ZeroLoc_i[1]; m<=ZeroLoc_i[1]+1; m++) {
					for (l=ZeroLoc_i[0]; l<=ZeroLoc_i[0]+1; l++) {
						if (l<0 || m<0 || n<0 || l>=Width_mi || m>=Height_mi || n>=Depth_mi) continue;
						loc[0] = n*WtimesH_mi + m*Width_mi + l;
						if (MinData_Ret > (float)Data_mT[loc[0]]) MinData_Ret = (float)Data_mT[loc[0]];
						if (MaxData_Ret < (float)Data_mT[loc[0]]) MaxData_Ret = (float)Data_mT[loc[0]];
						if (MinGM_Ret > GradientMag_mf[loc[0]]) MinGM_Ret = GradientMag_mf[loc[0]];
						if (MaxGM_Ret < GradientMag_mf[loc[0]]) MaxGM_Ret = GradientMag_mf[loc[0]];
						if (MinSD_Ret > SecondDerivative_mf[loc[0]]) MinSD_Ret = SecondDerivative_mf[loc[0]];
						if (MaxSD_Ret < SecondDerivative_mf[loc[0]]) MaxSD_Ret = SecondDerivative_mf[loc[0]];
					}
				}
			}
		}

		printf ("Min & Max Data Values at the hit boundaries = %f, %f\n", MinData_Ret, MaxData_Ret);
		printf ("Min & Max Gradient Magnitudes at the hit boundaries = %f, %f\n", MinGM_Ret, MaxGM_Ret);
		printf ("Min & Max Second Derivatives at the hit boundaries = %f, %f\n", MinSD_Ret, MaxSD_Ret);
		fflush (stdout);

	}
}


template<class _DataType>
void cTFGeneration<_DataType>::MinMaxComputationAtHitLocations(double GradThreshold, 
								float &MinData_Ret, float &MaxData_Ret, 
								float &MinGM_Ret, 	float &MaxGM_Ret, 
								float &MinSD_Ret, float &MaxSD_Ret,
								float &MinDist_Ret, float &MaxDist_Ret)
{
	int		i, k, loc[3], DataCoor_i[3], FoundZeroCrossingLoc_i;
	double	GradVec_d[3], CurrDataLoc_d[3], ZeroCrossingLoc_d[3], LocalMaxGradient_d;
	double	DataPosFromZeroCrossingLoc_d, Tempd;
	int		l, m, n, CurrZ = 0, MaxAxisResolution_i;
	int		ZeroLoc_i[3];
	float	MinData_f, MaxData_f;
	float	MinGM_f, MaxGM_f;
	float	MinSD_f, MaxSD_f;
	float	MinDist_f, MaxDist_f;
	
	
	MinData_f = MinGM_f = MinSD_f = MinDist_f = 99999.0;
	MaxData_f = MaxGM_f = MaxSD_f = MaxDist_f = -99999.0;

	MaxAxisResolution_i = Width_mi;
	if (MaxAxisResolution_i < Height_mi) MaxAxisResolution_i = Height_mi;
	if (MaxAxisResolution_i < Depth_mi) MaxAxisResolution_i = Depth_mi;

	for (i=0; i<WHD_mi; i++) {
	
		if (GradientMag_mf[i] < GradThreshold) {
			DistanceToZeroCrossingLoc_mf[i] = (float)MaxAxisResolution_i;
			continue;
		}
		DataCoor_i[2] = i/WtimesH_mi;
		DataCoor_i[1] = (i - DataCoor_i[2]*WtimesH_mi)/Width_mi;
		DataCoor_i[0] = i % Width_mi;
		if (CurrZ<=DataCoor_i[2]) {
			printf ("Z = %3d, ", DataCoor_i[2]);
			printf ("\n"); fflush(stdout);
			CurrZ++;
		}
		//-------------------------------------------------------------------------------------------
		// Finding the local maxima of gradient magnitudes along the gradient direction. 
		// It climbs the mountain of gradient magnitudes to find the zero-crossing second derivative 
		// It returns ZeroCrossingLoc_d, LocalMaxGradient, DataPosFromZeroCrossingLoc_d
		//-------------------------------------------------------------------------------------------
		// Getting the gradient vector at the position
		for (k=0; k<3; k++) GradVec_d[k] = (double)GradientVec_mf[i*3 + k];
		Tempd = sqrt (GradVec_d[0]*GradVec_d[0] + GradVec_d[1]*GradVec_d[1] + GradVec_d[2]*GradVec_d[2]);
		if (fabs(Tempd)<1e-6) {
			DistanceToZeroCrossingLoc_mf[i] = (float)MaxAxisResolution_i;
			continue; // To skip zero-length vectors
		}
		
		for (k=0; k<3; k++) GradVec_d[k] /= Tempd; // Normalizing the gradient vector
		for (k=0; k<3; k++) CurrDataLoc_d[k] = (double)DataCoor_i[k];
		FoundZeroCrossingLoc_i = FindZeroCrossingLocation(CurrDataLoc_d, GradVec_d, ZeroCrossingLoc_d, 
														LocalMaxGradient_d, DataPosFromZeroCrossingLoc_d);
		
		if (FoundZeroCrossingLoc_i) {
			if (LocalMaxGradient_d >= GradThreshold) {
			
		
//				printf ("MinMaxComputationAtHitLocations -- 2 -- 1\n"); fflush (stdout);
		
				DistanceToZeroCrossingLoc_mf[i] = (float)DataPosFromZeroCrossingLoc_d;
			
				if (MinDist_f > (float)DataPosFromZeroCrossingLoc_d) MinDist_f = (float)DataPosFromZeroCrossingLoc_d;
				if (MaxDist_f < (float)DataPosFromZeroCrossingLoc_d) MaxDist_f = (float)DataPosFromZeroCrossingLoc_d;
				
		
//				printf ("MinMaxComputationAtHitLocations -- 2 -- 2\n"); fflush (stdout);
		
				ZeroLoc_i[0] = (int)floor(ZeroCrossingLoc_d[0]);
				ZeroLoc_i[1] = (int)floor(ZeroCrossingLoc_d[1]);
				ZeroLoc_i[2] = (int)floor(ZeroCrossingLoc_d[2]);
				
		
//				printf ("MinMaxComputationAtHitLocations -- 2 -- 3\n"); fflush (stdout);
		
//				printf ("ZeroLoc_i = %d %d %d\n", ZeroLoc_i[0], ZeroLoc_i[1], ZeroLoc_i[2]); fflush (stdout);
		
				for (n=ZeroLoc_i[2]; n<=ZeroLoc_i[2]+1; n++) {
					for (m=ZeroLoc_i[1]; m<=ZeroLoc_i[1]+1; m++) {
						for (l=ZeroLoc_i[0]; l<=ZeroLoc_i[0]+1; l++) {
							if (l<0 || m<0 || n<0 || l>=Width_mi || m>=Height_mi || n>=Depth_mi) continue;
							loc[0] = n*WtimesH_mi + m*Width_mi + l;
							if (MinData_f > (float)Data_mT[loc[0]]) MinData_f = (float)Data_mT[loc[0]];
							if (MaxData_f < (float)Data_mT[loc[0]]) MaxData_f = (float)Data_mT[loc[0]];
							if (MinGM_f > GradientMag_mf[loc[0]]) MinGM_f = GradientMag_mf[loc[0]];
							if (MaxGM_f < GradientMag_mf[loc[0]]) MaxGM_f = GradientMag_mf[loc[0]];
							if (MinSD_f > SecondDerivative_mf[loc[0]]) MinSD_f = SecondDerivative_mf[loc[0]];
							if (MaxSD_f < SecondDerivative_mf[loc[0]]) MaxSD_f = SecondDerivative_mf[loc[0]];
						}
					}
				}
				// Treating the cube
		
//				printf ("MinMaxComputationAtHitLocations -- 2 -- 4\n"); fflush (stdout);
		

			}
			else DistanceToZeroCrossingLoc_mf[i] = (float)MaxAxisResolution_i;
			
			if (DistanceToZeroCrossingLoc_mf[i] > (float)MaxAxisResolution_i) {
				DistanceToZeroCrossingLoc_mf[i] = (float)MaxAxisResolution_i;
			}
		}
		else DistanceToZeroCrossingLoc_mf[i] = (float)MaxAxisResolution_i;
		
//		printf ("MinMaxComputationAtHitLocations -- 3\n"); fflush (stdout);
		

	}
	MinData_Ret = MinData_f;
	MaxData_Ret = MaxData_f;
	MinGM_Ret =   MinGM_f;
	MaxGM_Ret =   MaxGM_f;
	MinSD_Ret =   MinSD_f;
	MaxSD_Ret =   MaxSD_f;
	MinDist_Ret = MinDist_f;
	MaxDist_Ret = MaxDist_f;

	printf ("Min & Max Data Values at the hit boundaries = %f, %f\n", MinData_f, MaxData_f);
	printf ("Min & Max Gradient Magnitudes at the hit boundaries = %f, %f\n", MinGM_f, MaxGM_f);
	printf ("Min & Max Second Derivatives at the hit boundaries = %f, %f\n", MinSD_f, MaxSD_f);
	printf ("Min & Max Distance from a voxel to the hit boundary = %f, %f\n", MinDist_f, MaxDist_f);
	fflush (stdout);
	
}



// Removing phantom edges in the second derivative
template<class _DataType>
unsigned char *cTFGeneration<_DataType>::ConnectedPositiveVoxels()
{
	int				i, j, k, l, m, n, Xi, Yi, Zi, loc[28];
	int				CCIndex_i, Idx, DataLoc_i;
	float			T_f[4];
	cStack<int>		PositiveVoxels_stack;
	int				NumVoxels;
	
	
	printf ("Connected Positive Voxels ... \n");
	fflush (stdout);
	
	Initialize_CCSurfaceIndexTable();
	
	CCIndex_i = 0;
	PositiveVoxels_stack.Clear();
	for (i=0; i<WHD_mi; i++) {

		if (CCSurfaceVolume_mi[i]!=0) continue;

		PositiveVoxels_stack.Push(i);

		NumVoxels = 0;
		CCIndex_i++;
		
		do {

			PositiveVoxels_stack.Pop(DataLoc_i);
			CCSurfaceVolume_mi[DataLoc_i] = CCIndex_i;
			
			NumVoxels++;
			
			if (NumVoxels%100000==0) {
				printf ("Num CC Voxels = %d\n", NumVoxels);
				fflush (stdout);
			}
			
			
			Zi = DataLoc_i/WtimesH_mi;
			Yi = (DataLoc_i - Zi*WtimesH_mi)/Width_mi;
			Xi = DataLoc_i % Width_mi;

			Idx = -1;
			for (n=Zi-1; n<=Zi+1; n++) {
				for (m=Yi-1; m<=Yi+1; m++) {
					for (l=Xi-1; l<=Xi+1; l++) {
						Idx++;
						if (l<0 || m<0 || n<0 || l>=Width_mi || m>=Height_mi || n>=Depth_mi) {
							loc[Idx] = -1;
						}
						else loc[Idx] = n*WtimesH_mi + m*Width_mi + l;
					}
				}
			}

			for (j=0; j<27; j++) {

				// if the current location has a positive voxel, then ...
				if (loc[j]<0) continue;
				if (CCSurfaceVolume_mi[loc[j]]==0) {
			
					switch (j) {

						// Diagonal Elements
						case 1: case 3: case 5: case 7: 
						case 9: case 11: case 15: case 17:
						case 19: case 21: case 23: case 25:
							T_f[0] = SecondDerivative_mf[loc[13]];
							T_f[1] = SecondDerivative_mf[loc[j]];
							if (loc[NeighborIdxTable[j][0]]<0 || loc[NeighborIdxTable[j][1]]<0) break;
							T_f[2] = SecondDerivative_mf[loc[NeighborIdxTable[j][0]]];
							T_f[3] = SecondDerivative_mf[loc[NeighborIdxTable[j][1]]];
							if (T_f[0]*T_f[1] > T_f[2]*T_f[3]) {
								PositiveVoxels_stack.Push(loc[j]);
							}
							break;

						// Directly Connected Elements
						case 4: case 10: case 12: case 14: case 16: case 22:
							PositiveVoxels_stack.Push(loc[j]);
							break;
						default: break;
					}
				}
			}

		} while (PositiveVoxels_stack.Size()>0);
		
		
		put_CCSurfaceIndex(CCIndex_i, NumVoxels);
		
//		printf ("Num Connected Voxels = %d, Index = %d\n", NumVoxels, CCIndex_i);
//		fflush (stdout);

	}

	QuickSort3Elements(CCSurfaceIndex_mi, CCIndex_i+1, 'Y');

	int	NewIndex[256*2], MaxNumNewIndex;
	
	for (k=0, j=255, i=CCIndex_i; i>=0; i--, j-=5, k++) {
		printf ("j = %5d, Index = %5d, ", j, CCSurfaceIndex_mi[i*3]);
		printf ("Num Voxels = %d, ", CCSurfaceIndex_mi[i*3+1]);
		if (j<=1) break;
		CCSurfaceIndex_mi[i*3+2] = j;
		printf ("New Index = %d, ", CCSurfaceIndex_mi[i*3+2]);
		printf ("\n"); fflush (stdout);
		
		NewIndex[k*2] = CCSurfaceIndex_mi[i*3]; // Old Index
		NewIndex[k*2+1] = CCSurfaceIndex_mi[i*3+2];	// New Index
		if (k>=255) break;
	}
	MaxNumNewIndex = k;
	printf ("\n");
	printf ("Max Num New Index = %d\n", MaxNumNewIndex);
	fflush (stdout);
	
	unsigned char	*SurfaceIndexVolume_uc = new unsigned char[WHD_mi];

	
	for (i=0; i<WHD_mi; i++) {
		
		for (j=0; j<MaxNumNewIndex; j++) {
			if (CCSurfaceVolume_mi[i]==NewIndex[j*2]) {
				SurfaceIndexVolume_uc[i] = NewIndex[j*2+1];
				break;
			}
		}
		if (j==MaxNumNewIndex) SurfaceIndexVolume_uc[i] = 0;
	}
	
	
	return SurfaceIndexVolume_uc;
}


template<class _DataType>
void cTFGeneration<_DataType>::Initialize_CCSurfaceIndexTable()
{
	delete [] CCSurfaceIndex_mi;
	MaxSize_CCSurfaceIndexTable_mi = 100;

	CCSurfaceIndex_mi = new int [MaxSize_CCSurfaceIndexTable_mi*3];
}

template<class _DataType>
void cTFGeneration<_DataType>::put_CCSurfaceIndex(int Index, int NumVoxels)
{
	if (MaxSize_CCSurfaceIndexTable_mi<=Index) {
		IncreaseSize_CCSurfaceIndexTable(Index);
		CCSurfaceIndex_mi[Index*3+0] = Index;
		CCSurfaceIndex_mi[Index*3+1] = NumVoxels;
		CCSurfaceIndex_mi[Index*3+2] = Index;
	}
	else {
		CCSurfaceIndex_mi[Index*3+0] = Index;
		CCSurfaceIndex_mi[Index*3+1] = NumVoxels;
		CCSurfaceIndex_mi[Index*3+2] = Index;
	}
}

template<class _DataType>
void cTFGeneration<_DataType>::IncreaseSize_CCSurfaceIndexTable(int CurrMaxCCIndex_i)
{
	int	i, *CCSurfaceIndexTable_tempi;


	if (MaxSize_CCSurfaceIndexTable_mi<=CurrMaxCCIndex_i) {
		
		MaxSize_CCSurfaceIndexTable_mi*=2;
		CCSurfaceIndexTable_tempi = new int [MaxSize_CCSurfaceIndexTable_mi*3];
		for (i=0; i<CurrMaxCCIndex_i; i++) {
			CCSurfaceIndexTable_tempi[i*3+0] = CCSurfaceIndex_mi[i*3+0];
			CCSurfaceIndexTable_tempi[i*3+1] = CCSurfaceIndex_mi[i*3+1];
			CCSurfaceIndexTable_tempi[i*3+2] = CCSurfaceIndex_mi[i*3+2];
		}
		for (i=CurrMaxCCIndex_i; i<MaxSize_CCSurfaceIndexTable_mi; i++) {
			CCSurfaceIndexTable_tempi[i*3+0] = -1;
			CCSurfaceIndexTable_tempi[i*3+1] = -1;
			CCSurfaceIndexTable_tempi[i*3+2] = -1;
		}
		delete [] CCSurfaceIndex_mi;
		CCSurfaceIndex_mi = CCSurfaceIndexTable_tempi;
	}
}



template<class _DataType>
void cTFGeneration<_DataType>::MarkingLocalMaxVoxels2(double GradThreshold)
{
	int			i, j, k, loc[27];
	int			Xi, Yi, Zi;
	double		ZeroLoc_d[3], SD_d[4], GradM_d[3], GradVec_d[3];
	double		Length_d, FrontLoc_d[3], BackLoc_d[3];
	int			Positive_i, Negative_i, Empty_i;
	
	
	printf ("GradThreshold = %f\n", GradThreshold);
	Empty_i = 0;
	Negative_i = 1;
	Positive_i = 255;
	
	CCSurfaceVolume_mi = new int [WHD_mi];
	for (i=0; i<WHD_mi; i++) CCSurfaceVolume_mi[i] = Empty_i;
	
	for (i=0; i<WHD_mi; i++) {
		
		Zi = i/WtimesH_mi;
		Yi = (i - Zi*WtimesH_mi)/Width_mi;
		Xi = i % Width_mi;

		if (Xi>=Width_mi-1 || Yi>=Height_mi-1 || Zi>=Depth_mi-1) continue;

		loc[0] = i;
		loc[1] = Zi*WtimesH_mi + Yi*Width_mi + Xi+1;
		loc[2] = Zi*WtimesH_mi + (Yi+1)*Width_mi + Xi;
		loc[3] = (Zi+1)*WtimesH_mi + Yi*Width_mi + Xi;
		for (k=0; k<4; k++) SD_d[k] = SecondDerivative_mf[loc[k]];

		for (j=1; j<=3; j++) {

			if (SD_d[0]*SD_d[j]<0) {
				switch (j) {
					case 1 : 
						ZeroLoc_d[0] = SD_d[0]/(SD_d[0] - SD_d[j]) + Xi;
						ZeroLoc_d[1] = (double)Yi;
						ZeroLoc_d[2] = (double)Zi;
						break;
					case 2 : 
						ZeroLoc_d[0] = (double)Xi;
						ZeroLoc_d[1] = SD_d[0]/(SD_d[0] - SD_d[j]) + Yi;
						ZeroLoc_d[2] = (double)Zi;
						break;
					case 3 : 
						ZeroLoc_d[0] = (double)Xi;
						ZeroLoc_d[1] = (double)Yi;
						ZeroLoc_d[2] = SD_d[0]/(SD_d[0] - SD_d[j]) + Zi;
						break;
					default: break;
				}

				GradVecInterpolation(ZeroLoc_d, GradVec_d);
				Length_d = sqrt (GradVec_d[0]*GradVec_d[0] + GradVec_d[1]*GradVec_d[1]
								 + GradVec_d[2]*GradVec_d[2]);
				if (fabs(Length_d)<1e-5) continue;

				for (k=0; k<3; k++) GradVec_d[k] /= Length_d;
				for (k=0; k<3; k++) {
					FrontLoc_d[k] = ZeroLoc_d[k] + GradVec_d[k]*0.1;
					BackLoc_d[k] = ZeroLoc_d[k] + GradVec_d[k]*(-0.1);
				}
				GradM_d[0] = GradientInterpolation(FrontLoc_d);
				GradM_d[1] = GradientInterpolation(ZeroLoc_d);
				GradM_d[2] = GradientInterpolation(BackLoc_d);
				// Local Maximum
				if (GradM_d[1]>GradM_d[0] && GradM_d[1]>GradM_d[2]) {

	//				printf ("<-- Local Max ");

					if (SD_d[0]>=0) {
						if (CCSurfaceVolume_mi[loc[0]]==Empty_i) CCSurfaceVolume_mi[loc[0]] = Positive_i;// Positive
						if (CCSurfaceVolume_mi[loc[j]]==Empty_i) CCSurfaceVolume_mi[loc[j]] = Negative_i;// Negative
					} else 
					if (SD_d[0]<0) {
						if (CCSurfaceVolume_mi[loc[0]]==Empty_i) CCSurfaceVolume_mi[loc[0]] = Negative_i;// Negative
						if (CCSurfaceVolume_mi[loc[j]]==Empty_i) CCSurfaceVolume_mi[loc[j]] = Positive_i;// Positive
					}
				}

	//			printf ("\n"); fflush (stdout);
			}
		}

	}


	int		NumPositive;
	for (i=0; i<WHD_mi; i++) {
		
		Zi = i/WtimesH_mi;
		Yi = (i - Zi*WtimesH_mi)/Width_mi;
		Xi = i % Width_mi;

		if (Xi<=0 || Yi<=0 || Zi<=0 || Xi>=Width_mi-1 || Yi>=Height_mi-1 || Zi>=Depth_mi-1) continue;
		if (CCSurfaceVolume_mi[i]!=Empty_i) continue;
		if (SecondDerivative_mf[i]<0) continue;

		// When SecondDerivative_mf[i]>0 and CCSurfaceVolume_mi[i]==Empty
		loc[0] = i;
		loc[1] = Zi*WtimesH_mi + Yi*Width_mi + Xi-1;
		loc[2] = Zi*WtimesH_mi + Yi*Width_mi + Xi+1;
		loc[3] = Zi*WtimesH_mi + (Yi-1)*Width_mi + Xi;
		loc[4] = Zi*WtimesH_mi + (Yi+1)*Width_mi + Xi;
		loc[5] = (Zi-1)*WtimesH_mi + Yi*Width_mi + Xi;
		loc[6] = (Zi+1)*WtimesH_mi + Yi*Width_mi + Xi;

		NumPositive = 0;
		for (j=1; j<=6; j++) {
			if (CCSurfaceVolume_mi[loc[j]]==Positive_i) NumPositive++;
		}
		if (NumPositive>=2) CCSurfaceVolume_mi[i] = Positive_i-100;
	}

	
}



template<class _DataType>
void cTFGeneration<_DataType>::ComputeConnectedSurfaceVolume(double GradThreshold, char *TargetFileName)
{
	int		i, k, loc[3], DataCoor[3], FoundZeroCrossingLoc_i;
	double	GradVec_d[3], CurrDataLoc_d[3], ZeroCrossingLoc_d[3], LocalMaxGradient_d;
	double	DataPosFromZeroCrossingLoc_d, Tempd;
	int		l, m, n, CubeIndex_i, MinCCIndex_i, CurrZ = 0;
	int		CCIndex_i, ZeroLoc_i[3], TempLoc_i[3];
	int		*CCVolume_i;
	

	printf ("Computing Connected Surface Volume ... \n");
	printf ("TargetFileName = %s\n", TargetFileName);
	fflush (stdout);
	
	Initialize_CCIndexTable();
	CCVolume_i = new int [WHD_mi];
	for (i=0; i<WHD_mi; i++) CCVolume_i[i] = -1;


	CCIndex_i = 1;
	for (i=0; i<WHD_mi; i++) {
	
		if (GradientMag_mf[i] < GradThreshold) continue;

		DataCoor[2] = i/WtimesH_mi;
		DataCoor[1] = (i - DataCoor[2]*WtimesH_mi)/Width_mi;
		DataCoor[0] = i % Width_mi;

		if (CurrZ==DataCoor[2]) {
			printf ("Z = %3d, ", DataCoor[2]);
			printf ("\n"); fflush(stdout);
			CurrZ++;
		}

		if (DataCoor[2]>=60 && DataCoor[2]<=65) {  }
		else continue;


		//-------------------------------------------------------------------------------------------
		// Finding the local maxima of gradient magnitudes along the gradient direction. 
		// It climbs the mountain of gradient magnitudes to find the zero-crossing second derivative 
		// It returns ZeroCrossingLoc_d, LocalMaxGradient, DataPosFromZeroCrossingLoc_d
		//-------------------------------------------------------------------------------------------
		// Getting the gradient vector at the position
		for (k=0; k<3; k++) GradVec_d[k] = (double)GradientVec_mf[i*3 + k];
		Tempd = sqrt (GradVec_d[0]*GradVec_d[0] + GradVec_d[1]*GradVec_d[1] + GradVec_d[2]*GradVec_d[2]);
		if (fabs(Tempd)<1e-6) continue; // To skip zero-length vectors
		for (k=0; k<3; k++) GradVec_d[k] /= Tempd; // Normalizing the gradient vector
		for (k=0; k<3; k++) CurrDataLoc_d[k] = (double)DataCoor[k];
		FoundZeroCrossingLoc_i = FindZeroCrossingLocation(CurrDataLoc_d, GradVec_d, ZeroCrossingLoc_d, 
														LocalMaxGradient_d, DataPosFromZeroCrossingLoc_d);
		if (FoundZeroCrossingLoc_i) {
			if (LocalMaxGradient_d >= GradThreshold) {
			

				printf ("Curr Data Loc = %.2f, %.2f, %.2f, ", CurrDataLoc_d[0], CurrDataLoc_d[1], CurrDataLoc_d[2]);
				printf ("Grad Vec = %.4f, %.4f, %.4f, ", GradVec_d[0], GradVec_d[1], GradVec_d[2]);
				printf ("Zero Loc = %.4f, %.4f, %.4f, ", ZeroCrossingLoc_d[0], ZeroCrossingLoc_d[1], ZeroCrossingLoc_d[2]);
				printf ("Max Grad Magnitude = %.4f, ", LocalMaxGradient_d);
				printf ("Data Position = %.4f, ", DataPosFromZeroCrossingLoc_d);
				printf ("\n"); fflush (stdout);

				
				ZeroLoc_i[0] = (int)floor(ZeroCrossingLoc_d[0]);
				ZeroLoc_i[1] = (int)floor(ZeroCrossingLoc_d[1]);
				ZeroLoc_i[2] = (int)floor(ZeroCrossingLoc_d[2]);
				
				CubeIndex_i = -1;
				for (n=ZeroLoc_i[2]; n<=ZeroLoc_i[2]+1; n++) {
					for (m=ZeroLoc_i[1]; m<=ZeroLoc_i[1]+1; m++) {
						for (l=ZeroLoc_i[0]; l<=ZeroLoc_i[0]+1; l++) {
							
							CubeIndex_i++;
							if (l>=Width_mi || m>=Height_mi || n>=Depth_mi) continue;
							TempLoc_i[0] = l;
							TempLoc_i[1] = m;
							TempLoc_i[2] = n;
							loc[0] = n*WtimesH_mi + m*Width_mi + l;
							if (SecondDerivative_mf[loc[0]]<0) continue;
							
							if (IsConnectedVoxel_To_ZeroCrossing(ZeroCrossingLoc_d, TempLoc_i)) {
							
								MinCCIndex_i = FindMinIndex_18Adjacency(CCVolume_i, TempLoc_i);
								
								
								if (MinCCIndex_i>0) printf ("MinCCIndex_i = %d\n", MinCCIndex_i);
								
								
								if (MinCCIndex_i>0) {
									if (CCIndexTable_mi[CCVolume_i[loc[0]]]>0)
										CCIndexTable_mi[CCVolume_i[loc[0]]] = MinCCIndex_i;
									if (CCIndexTable_mi[CCVolume_i[i]]>0)
										CCIndexTable_mi[CCVolume_i[i]] = MinCCIndex_i;
									CCVolume_i[loc[0]] = MinCCIndex_i;
									CCVolume_i[i] = MinCCIndex_i;
								}
								else {
									CCVolume_i[loc[0]] = CCIndex_i;
									CCVolume_i[i] = CCIndex_i;
									CCIndex_i++;
									if (MaxSize_CCIndexTable_mi<=CCIndex_i) {
									
										printf ("Increasing MaxSize_CCIndexTable_mi = %d\n", MaxSize_CCIndexTable_mi*2);
										fflush (stdout);
										
										IncreaseSize_CCIndexTable(CCIndex_i);
									}
								}
							}

						}
					}
				}
				// Treating the cube
			}
		}
		
	}


	unsigned char	*CCVolume_uc;
	CCVolume_uc = Rearrange_CCVolume(CCVolume_i);
//	SaveVolume(CCVolume_uc, 0.0, 255.0, "CCVolume");
	delete [] CCVolume_uc;

	
	printf ("Computing Connected Surface Volume is done \n");
	fflush (stdout);

}



template<class _DataType>
void cTFGeneration<_DataType>::Initialize_CCIndexTable()
{
	delete [] CCIndexTable_mi;
	MaxSize_CCIndexTable_mi = 100;

	CCIndexTable_mi = new int [MaxSize_CCIndexTable_mi];
	
}


template<class _DataType>
void cTFGeneration<_DataType>::IncreaseSize_CCIndexTable(int CurrMaxCCIndex_i)
{
	int	i, *CCIndexTable_tempi;

	if (MaxSize_CCIndexTable_mi<=CurrMaxCCIndex_i) {
		
		MaxSize_CCIndexTable_mi*=2;
		CCIndexTable_tempi = new int [MaxSize_CCIndexTable_mi];
		
		for (i=0; i<CurrMaxCCIndex_i; i++) {
			CCIndexTable_tempi[i] = CCIndexTable_mi[i];
		}
		for (i=CurrMaxCCIndex_i; i<MaxSize_CCIndexTable_mi; i++) {
			CCIndexTable_tempi[i] = -1;
		}
		delete [] CCIndexTable_mi;
		CCIndexTable_mi = CCIndexTable_tempi;
	}
}


template<class _DataType>
unsigned char *cTFGeneration<_DataType>::Rearrange_CCVolume(int *CCVolume_i)
{
	int				i, j, SurfaceIdx, NewIdx;
	int				*Freq_Index_i, *ReducedFreq_Index_i;
	unsigned char	*CCVolume_uc;
	

	for (i=0; i<WHD_mi; i++) {
		SurfaceIdx = CCVolume_i[i];
		if (SurfaceIdx>0) {
			if (CCIndexTable_mi[SurfaceIdx]>0) {
				CCVolume_i[i] = CCIndexTable_mi[SurfaceIdx];
			}
		}
	}

	Freq_Index_i = new int [MaxSize_CCIndexTable_mi*3];
	for (i=0; i<MaxSize_CCIndexTable_mi; i++) {
		Freq_Index_i[i*3 + 0] = 0;		// Frequencies
		Freq_Index_i[i*3 + 1] = i+1;	// Index
		Freq_Index_i[i*3 + 2] = i+1;	// Index
	}

	// Computing Frequencies
	for (i=0; i<WHD_mi; i++) {
		SurfaceIdx = CCVolume_i[i];
		if (SurfaceIdx>0) Freq_Index_i[SurfaceIdx*3]++;
	}

#ifdef	DEBUG_TF_CCVOLUME
		printf ("Before the quick sorting \n");
		fflush (stdout);
#endif							

	int		NumBins = 0;
	for (i=0; i<MaxSize_CCIndexTable_mi; i++) {
		if (Freq_Index_i[i*3 + 0]>1) {
			NumBins++;
			printf ("i = %d, ", i);
			printf ("Freq = %d, ", Freq_Index_i[i*3 + 0]);
			printf ("Index = %d, ", Freq_Index_i[i*3 + 1]);
			printf ("\n"); fflush (stdout);
		}
	}
	
	ReducedFreq_Index_i = new int [NumBins*3];
	for (j=0, i=0; i<MaxSize_CCIndexTable_mi; i++) {
		if (Freq_Index_i[i*3 + 0]>1) {
			ReducedFreq_Index_i[j*3] = Freq_Index_i[i*3];
			ReducedFreq_Index_i[j*3+1] = Freq_Index_i[i*3+1];
			ReducedFreq_Index_i[j*3+2] = Freq_Index_i[i*3+2];
			j++;
		}
	}	

	QuickSort3Elements (ReducedFreq_Index_i, NumBins, 'X', 'Y', 'Z');
	
#ifdef	DEBUG_TF_CCVOLUME
		printf ("After the sorting \n");
		fflush (stdout);
#endif							

	for (i=0; i<NumBins; i++) {
		if (ReducedFreq_Index_i[i*3 + 0]>1) {
			printf ("i = %d, ", i);
			printf ("Freq = %d, ", ReducedFreq_Index_i[i*3 + 0]);
			printf ("Index = %d, ", ReducedFreq_Index_i[i*3 + 1]);
			printf ("Index2 = %d, ", ReducedFreq_Index_i[i*3 + 2]);
			printf ("\n"); fflush (stdout);
		}
	}

	for (i=0; i<WHD_mi; i++) {
		SurfaceIdx = CCVolume_i[i];
		NewIdx = ReducedFreq_Index_i[SurfaceIdx*3 + 1];
		CCVolume_i[i] = NewIdx;
	}
	
	CCVolume_uc = new unsigned char [WHD_mi];
	for (i=0; i<WHD_mi; i++) {
		if (CCVolume_i[i]<255) CCVolume_uc[i] = CCVolume_i[i];
		else CCVolume_uc[i] = 255;
	}


	delete [] Freq_Index_i;
	delete [] ReducedFreq_Index_i;
	
	return CCVolume_uc;
	return NULL;
}


template<class _DataType>
int cTFGeneration<_DataType>::FindMinIndex_18Adjacency(int *CCVolume_i, int *TempLoc3_i)
{
	int		Idx, l, m, n, loc[27], CubeIdx;
	int		MinIdx;
	float	T_f[5];
	
	
	CubeIdx = -1;
	for (n=TempLoc3_i[2]-1; n<=TempLoc3_i[2]+1; n++) {
		for (m=TempLoc3_i[1]-1; m<=TempLoc3_i[1]+1; m++) {
			for (l=TempLoc3_i[0]-1; l<=TempLoc3_i[0]+1; l++) {
				CubeIdx++;
				if (l<0 || m<0 || n<0 || l>=Width_mi || m>=Height_mi || n>=Depth_mi) loc[CubeIdx] = -1;
				else loc[CubeIdx] = n*WtimesH_mi + m*Width_mi + l;
			}
		}
	}
	
	MinIdx = -1;
	for (Idx=0; Idx<27; Idx++) {
		
		if (loc[CubeIdx]<0) continue;
		if (SecondDerivative_mf[loc[Idx]]<0) continue;
		
		switch (Idx) {

			// Diagonal Elements
			case 1: case 3: case 5: case 7: 
			case 9: case 11: case 15: case 17:
			case 19: case 21: case 23: case 25:
				T_f[0] = SecondDerivative_mf[loc[13]];
				T_f[1] = SecondDerivative_mf[loc[Idx]];
				T_f[2] = SecondDerivative_mf[loc[NeighborIdxTable[Idx][0]]];
				T_f[3] = SecondDerivative_mf[loc[NeighborIdxTable[Idx][1]]];
				if (T_f[2]<0 && T_f[3]<0) {
					if (T_f[0]*T_f[1] > T_f[2]*T_f[3]) {
						if (MinIdx > CCVolume_i[loc[Idx]] && CCVolume_i[loc[Idx]]>0)
							MinIdx = CCVolume_i[loc[Idx]];
					}
				}
				break;

			// Directly Connected Elements
			case 4: case 10: case 12: case 14: case 16: case 22:
				T_f[0] = SecondDerivative_mf[loc[13]];
				T_f[1] = SecondDerivative_mf[loc[Idx]];
				if (MinIdx > CCVolume_i[loc[Idx]] && CCVolume_i[loc[Idx]]>0)
					MinIdx = CCVolume_i[loc[Idx]];
				break;

			default: break;
		}
	}

	return 	MinIdx;
	
}


template<class _DataType>
int cTFGeneration<_DataType>::IsConnectedVoxel_To_ZeroCrossing(double *ZeroCrossingLoc3_d, int *TempLoc3_i)
{
	int		k;
	double	VectorTowardVoxel_d[3], Length_d;
	double	FrontLoc_d[3], FrontSD_d;
	double	Front2Loc_d[3], Front2SD_d;


	for (k=0; k<3; k++) {
		VectorTowardVoxel_d[k] = (double)TempLoc3_i[k] - ZeroCrossingLoc3_d[k];
	}
	Length_d = sqrt(VectorTowardVoxel_d[0]*VectorTowardVoxel_d[0] +
					VectorTowardVoxel_d[1]*VectorTowardVoxel_d[1] +
					VectorTowardVoxel_d[2]*VectorTowardVoxel_d[2]);
	// Normalizing the vector
	for (k=0; k<3; k++) VectorTowardVoxel_d[k] /= Length_d;
	for (k=0; k<3; k++) {
		FrontLoc_d[k] = ZeroCrossingLoc3_d[k] + VectorTowardVoxel_d[k]*0.1;
		Front2Loc_d[k] = ZeroCrossingLoc3_d[k] + VectorTowardVoxel_d[k]*0.2;
	}

	FrontSD_d = SecondDInterpolation(FrontLoc_d);
	Front2SD_d = SecondDInterpolation(Front2Loc_d);

	if (FrontSD_d>0 && Front2SD_d>0) return true;
	else return false;
	
	return false;	// To suppress the compile warning
}


template<class _DataType>
void cTFGeneration<_DataType>::HistoVolumeSigma_Evaluation(double GradThreshold, int NumElements, 
												int DoSigmaEvaluation, char *TargetFileName)
{
	int		i, loc[3], DataCoor[3], ithSp;
	int		CurrZ;
	double	Sigma_d[5], SigmaFromHistogram_d;
	double	CenterLoc_d[5][3], ZeroDist_d[5], Distance_d;
	double	HistogramFactorI_d, HistogramFactorG_d, HistogramFactorH_d;
	double	Den_d, FD_d, GradMag_d;
	double	P_vg_d, H_vg_d;
	int		NumBinsI_i, NumBinsG_i, NumBinsH_i, MaxNumBins_i;
	int		v_i, g_i, NumAveRealDist_i[256], NumAveP_vg_i[256], NumAveDistToZero_i[256];
	int		NumRealDist_Pvg_Diff_i[256], NumRealDist_DistToZero_Diff_i[256];
	double	AveRealDist_d[256], AveP_vg_d[256], AveDistToZero_d[256];
	double	RealDist_Pvg_Diff_d[256], RealDist_DistToZero_Diff_d[256];



	delete [] H_vg_md;
	// Computing G_v[] & H_vg_md[]
	ComputeHistogramVolume(NumElements, SigmaFromHistogram_d);
	printf ("Computed Sigma = %.6f\n", SigmaFromHistogram_d);

	float	MinData_f, MaxData_f, MinGM_f, MaxGM_f, MinSD_f, MaxSD_f;
	// Computing distance and save it into DistanceToZeroCrossingLoc_mf[]
	ComputeAndSaveDistanceFloat(GradThreshold, MinData_f, MaxData_f, MinGM_f, MaxGM_f, 
								MinSD_f, MaxSD_f, TargetFileName);

	if (DoSigmaEvaluation) {	}	// Continue this function
	else return;

	HistogramFactorI_d = (double)NumElements/(MaxData_mf - MinData_mf);
	HistogramFactorG_d = (double)NumElements/(MaxGrad_mf - MinGrad_mf);
	HistogramFactorH_d = (double)NumElements/(MaxSecMag_mf - MinSecMag_mf);

	NumBinsI_i = (int)((MaxData_mf - MinData_mf)*HistogramFactorI_d);
	NumBinsG_i = (int)((MaxGrad_mf - MinGrad_mf)*HistogramFactorG_d);
	NumBinsH_i = (int)((MaxSecMag_mf - MinSecMag_mf)*HistogramFactorH_d);

	MaxNumBins_i = NumBinsI_i;
	if (MaxNumBins_i < NumBinsG_i) MaxNumBins_i = NumBinsG_i;
	if (MaxNumBins_i < NumBinsH_i) MaxNumBins_i = NumBinsH_i;

	printf ("NumBinsI_i = %d, ", NumBinsI_i);
	printf ("NumBinsG_i = %d, ", NumBinsG_i);
	printf ("NumBinsH_i = %d, ", NumBinsH_i);
	printf ("MaxNumBins_i = %d, ", MaxNumBins_i);
	printf ("\n"); fflush (stdout);

	for (i=0; i<256; i++) {
		AveRealDist_d[i] = 0.0;			NumAveRealDist_i[i] = 0;
		AveP_vg_d[i] = 0.0;				NumAveP_vg_i[i] = 0;
		AveDistToZero_d[i] = 0.0;		NumAveDistToZero_i[i] = 0;	
		RealDist_Pvg_Diff_d[i] = 0.0;			NumRealDist_Pvg_Diff_i[i] = 0;	
		RealDist_DistToZero_Diff_d[i] = 0.0;	NumRealDist_DistToZero_Diff_i[i] = 0;	
	}

	ithSp = 0;
	CenterLoc_d[ithSp][0] = 40.0;
//	CenterLoc_d[ithSp][0] = 64.0;
	CenterLoc_d[ithSp][1] = 64.0;
	CenterLoc_d[ithSp][2] = 64.0;
	Sigma_d[ithSp] = 3.0;
	ZeroDist_d[ithSp] = 20.0;
	
	ithSp = 1;
	CenterLoc_d[ithSp][0] = 85.0;
//	CenterLoc_d[ithSp][0] = 64.0;
	CenterLoc_d[ithSp][1] = 64.0;
	CenterLoc_d[ithSp][2] = 64.0;
	Sigma_d[ithSp] = 1.0;
	ZeroDist_d[ithSp] = 20.0;

	int		PrintCount_i = 0;
	CurrZ = 0;
	for (i=0; i<WHD_mi; i++) {
	
		if (GradientMag_mf[i] < GradThreshold) continue;

		DataCoor[2] = i/WtimesH_mi;
		DataCoor[1] = (i - DataCoor[2]*WtimesH_mi)/Width_mi;
		DataCoor[0] = i % Width_mi;
		if (CurrZ<=DataCoor[2]) {
			printf ("Z = %3d, ", DataCoor[2]);
			printf ("\n"); fflush (stdout);
			CurrZ++;
		}
		if (DataCoor[0]<64) ithSp = 0;
		else ithSp = 1;
		Distance_d = sqrt(	(CenterLoc_d[ithSp][0]-DataCoor[0])*(CenterLoc_d[ithSp][0]-DataCoor[0]) + 
							(CenterLoc_d[ithSp][1]-DataCoor[1])*(CenterLoc_d[ithSp][1]-DataCoor[1]) +
							(CenterLoc_d[ithSp][2]-DataCoor[2])*(CenterLoc_d[ithSp][2]-DataCoor[2])	);
		Den_d = ((double)Data_mT[i]-MinData_mf)*HistogramFactorI_d;			// Intensity Values
		FD_d = ((double)GradientMag_mf[i]-MinGrad_mf)*HistogramFactorG_d;		// Gradient Magnitudes

		g_i = (int)(FD_d+0.5);
		v_i = (int)(Den_d+0.5);

		GradMag_d = G_v_md[v_i];
		H_vg_d = (double)H_vg_md[g_i*(MaxNumBins_i+1) + v_i];
		if (GradMag_d - 0.001 <= 1e-5) P_vg_d = 255.0;
		else P_vg_d = -Sigma_d[ithSp]*Sigma_d[ithSp]*H_vg_d/(GradMag_d-0.001);
//		else P_vg_d = -SigmaFromHistogram_d*SigmaFromHistogram_d*H_vg_d/(GradMag_d-0.001);

		loc[0] = (int)(fabs(ZeroDist_d[ithSp] - Distance_d) + 0.5);
		AveRealDist_d[loc[0]] += fabs(ZeroDist_d[ithSp] - Distance_d);
		NumAveRealDist_i[loc[0]]++;

		AveP_vg_d[loc[0]] += P_vg_d;
		NumAveP_vg_i[loc[0]]++;

		if (fabs(DistanceToZeroCrossingLoc_mf[i])<20.0) {
			AveDistToZero_d[loc[0]] += (double)fabs(DistanceToZeroCrossingLoc_mf[i]);
			NumAveDistToZero_i[loc[0]]++;
		}
		
		RealDist_Pvg_Diff_d[loc[0]] += fabs(fabs(P_vg_d) - fabs(ZeroDist_d[ithSp]-Distance_d));
		NumRealDist_Pvg_Diff_i[loc[0]]++;
		
		RealDist_DistToZero_Diff_d[loc[0]] += fabs(fabs(DistanceToZeroCrossingLoc_mf[i])-fabs((ZeroDist_d[ithSp]-Distance_d)));
		NumRealDist_DistToZero_Diff_i[loc[0]]++;
		if (DataCoor[1]==64 && DataCoor[2]==64) {
			PrintCount_i++;
			printf ("%3d,%3d,%3d, ", DataCoor[0], DataCoor[1], DataCoor[2]);
			printf ("DistFromC = %7.2f ", Distance_d);
			printf ("RealDistTo0 = %9.4f, ", fabs(ZeroDist_d[ithSp] - Distance_d));
			printf ("CompDistTo0 = %9.5f, ", DistanceToZeroCrossingLoc_mf[i]);
			printf ("SD = %9.5f, ", SecondDerivative_mf[i]);
			printf ("FD = %9.5f, ", GradientMag_mf[i]);
			printf ("P_vg = %9.5f, ", P_vg_d);
			printf ("Diff P_vg = %9.5f, ", fabs(fabs(P_vg_d) - fabs((ZeroDist_d[ithSp] - Distance_d))));
			printf ("Diff DistTo0 = %9.5f, ", fabs(fabs(DistanceToZeroCrossingLoc_mf[i])-fabs((ZeroDist_d[ithSp]-Distance_d))));
			printf ("\n"); fflush (stdout);
		}
	}
	
	for (i=0; i<256; i++) {
		if (NumAveRealDist_i[i]>0) AveRealDist_d[i] /= (double)NumAveRealDist_i[i];
		if (NumAveDistToZero_i[i]>0) AveDistToZero_d[i] /= (double)NumAveDistToZero_i[i];
		if (NumAveP_vg_i[i]>0) AveP_vg_d[i] /= (double)NumAveP_vg_i[i];
		if (NumRealDist_Pvg_Diff_i[i]>0) RealDist_Pvg_Diff_d[i] /= (double)NumRealDist_Pvg_Diff_i[i];
		if (NumRealDist_DistToZero_Diff_i[i]>0) RealDist_DistToZero_Diff_d[i] /= (double)NumRealDist_DistToZero_Diff_i[i];
	}
	
	for (i=0; i<256; i++) {
		if (NumAveRealDist_i[i]>0 || NumAveP_vg_i[i]>0 || NumAveDistToZero_i[i]>0 ||
			NumRealDist_Pvg_Diff_i[i]>0 || NumRealDist_DistToZero_Diff_i[i]>0) {
			printf ("Dist = %3d ", i);
			printf ("RealDistTo0 = %9.5f, #%5d ", AveRealDist_d[i], NumAveRealDist_i[i]);
			printf ("CompDistTo0 = %9.5f, #%5d ", AveDistToZero_d[i], NumAveDistToZero_i[i]);
			printf ("P_vg = %9.5f, #%5d ", AveP_vg_d[i], NumAveP_vg_i[i]);
			printf ("Diff P_vg = %9.5f, ", RealDist_Pvg_Diff_d[i]);
			printf ("Diff DistTo0 = %9.5f, #%5d ", RealDist_DistToZero_Diff_d[i], NumRealDist_DistToZero_Diff_i[i]);
			printf ("HitRatio = %5.2f %%", (double)NumAveDistToZero_i[i]/NumAveRealDist_i[i]*100.0);
			printf ("\n"); fflush (stdout);
		}
	}
	
	printf ("\n");
	for (i=0; i<256; i++) {
		if (NumAveRealDist_i[i]>0 || NumAveP_vg_i[i]>0 || NumAveDistToZero_i[i]>0 ||
			NumRealDist_Pvg_Diff_i[i]>0 || NumRealDist_DistToZero_Diff_i[i]>0) {
			printf ("Dist = %3d ", i);
			printf ("RealDistTo0 = %9.5f, #%5d ", AveRealDist_d[i], NumAveRealDist_i[i]);
			printf ("CompDistTo0 = %9.5f, #%5d ", AveDistToZero_d[i], NumAveDistToZero_i[i]);
			printf ("P_vg = %9.5f, #%5d ", AveP_vg_d[i], NumAveP_vg_i[i]);
			printf ("Diff P_vg = %9.5f, ", fabs(fabs(AveRealDist_d[i])-fabs(AveP_vg_d[i])));
			printf ("Diff DistTo0 = %9.5f, ", fabs(fabs(AveRealDist_d[i])-fabs(AveDistToZero_d[i])));
			printf ("HitRatio = %5.2f %%", (double)NumAveDistToZero_i[i]/NumAveRealDist_i[i]*100.0);
			printf ("\n"); fflush (stdout);
		}
	}
	

}

template<class _DataType>
void cTFGeneration<_DataType>::HistoVolumeSigma_Evaluation2(double GradThreshold, int NumElements, char *TargetName,
								int DoSigmaEvaluation, char *TargetFileName)
{
	int		i, loc[3], DataCoor[3];
	int		CurrZ;
	double	Sigma_d;
	double	CenterLoc_d[3], Distance_d;
	double	HistogramFactorI_d, HistogramFactorG_d, HistogramFactorH_d;
	double	Den_d, FD_d;
	double	P_vg_d, H_vg_d;
	int		NumBinsI_i, NumBinsG_i, NumBinsH_i, MaxNumBins_i;
	int		v_i, g_i, NumAveRealDist_i[256], NumAveP_vg_i[256], NumAveDistToZero_i[256];
	int		NumRealDist_Pvg_Diff_i[256], NumRealDist_DistToZero_Diff_i[256];
	double	AveRealDist_d[256], AveP_vg_d[256], AveDistToZero_d[256];
	double	RealDist_Pvg_Diff_d[256], RealDist_DistToZero_Diff_d[256];
	double	ZeroDist_d;


	printf ("TargetName = %s\n", TargetName);
	delete [] H_vg_md;
	// Computing G_v[] & H_vg_md[]
	ComputeHistogramVolume(NumElements, Sigma_d);

	printf ("Computed Sigma = %.6f\n", Sigma_d);
	Sigma_d = 3.0;
	printf ("Assigned New Sigma = %.6f\n", Sigma_d);

	float	MinData_f, MaxData_f, MinGM_f, MaxGM_f, MinSD_f, MaxSD_f;
	// Computing distance and save it into DistanceToZeroCrossingLoc_mf[]
	ComputeAndSaveDistanceFloat(GradThreshold, MinData_f, MaxData_f, MinGM_f, MaxGM_f, 
								MinSD_f, MaxSD_f, TargetFileName);

	if (DoSigmaEvaluation) {	}	// Continue this function
	else return;

	HistogramFactorI_d = (double)NumElements/(MaxData_mf - MinData_mf);
	HistogramFactorG_d = (double)NumElements/(MaxGrad_mf - MinGrad_mf);
	HistogramFactorH_d = (double)NumElements/(MaxSecMag_mf - MinSecMag_mf);

	NumBinsI_i = (int)((MaxData_mf - MinData_mf)*HistogramFactorI_d);
	NumBinsG_i = (int)((MaxGrad_mf - MinGrad_mf)*HistogramFactorG_d);
	NumBinsH_i = (int)((MaxSecMag_mf - MinSecMag_mf)*HistogramFactorH_d);

	MaxNumBins_i = NumBinsI_i;
	if (MaxNumBins_i < NumBinsG_i) MaxNumBins_i = NumBinsG_i;
	if (MaxNumBins_i < NumBinsH_i) MaxNumBins_i = NumBinsH_i;

	printf ("NumBinsI_i = %d, ", NumBinsI_i);
	printf ("NumBinsG_i = %d, ", NumBinsG_i);
	printf ("NumBinsH_i = %d, ", NumBinsH_i);
	printf ("MaxNumBins_i = %d, ", MaxNumBins_i);
	printf ("\n"); fflush (stdout);

	for (i=0; i<256; i++) {
		AveRealDist_d[i] = 0.0;			NumAveRealDist_i[i] = 0;
		AveP_vg_d[i] = 0.0;				NumAveP_vg_i[i] = 0;
		AveDistToZero_d[i] = 0.0;		NumAveDistToZero_i[i] = 0;	
		RealDist_Pvg_Diff_d[i] = 0.0;			NumRealDist_Pvg_Diff_i[i] = 0;	
		RealDist_DistToZero_Diff_d[i] = 0.0;	NumRealDist_DistToZero_Diff_i[i] = 0;	
	}

	CenterLoc_d[0] = 64.0;
	CenterLoc_d[1] = 64.0;
	CenterLoc_d[2] = 64.0;
	ZeroDist_d = 30.0;
	
	CurrZ = 0;
	for (i=0; i<WHD_mi; i++) {
	
		if (GradientMag_mf[i] < GradThreshold) continue;

		DataCoor[2] = i/WtimesH_mi;
		DataCoor[1] = (i - DataCoor[2]*WtimesH_mi)/Width_mi;
		DataCoor[0] = i % Width_mi;
		if (CurrZ<=DataCoor[2]) {
			printf ("Z = %3d, ", DataCoor[2]);
			printf ("\n"); fflush (stdout);
			CurrZ++;
		}
		
		Distance_d = sqrt(	(CenterLoc_d[0]-DataCoor[0])*(CenterLoc_d[0]-DataCoor[0]) + 
							(CenterLoc_d[1]-DataCoor[1])*(CenterLoc_d[1]-DataCoor[1]) +
							(CenterLoc_d[2]-DataCoor[2])*(CenterLoc_d[2]-DataCoor[2])	);
							
		Den_d = ((double)Data_mT[i]-MinData_mf)*HistogramFactorI_d;			// Intensity Values
		FD_d = ((double)GradientMag_mf[i]-MinGrad_mf)*HistogramFactorG_d;		// Gradient Magnitudes
//			SD_d = ((double)SecondDerivative_mf[i]-MinSecMag_mf)*HistogramFactorH_d;// 2nd derivative
		v_i = (int)(Den_d+0.5);
		g_i = (int)(FD_d+0.5);

		H_vg_d = (double)H_vg_md[g_i*(MaxNumBins_i+1) + v_i];
		P_vg_d = -Sigma_d*Sigma_d*H_vg_d/GradientMag_mf[i];

		loc[0] = (int)(Distance_d + 0.5);
		AveRealDist_d[loc[0]] += (ZeroDist_d - Distance_d);
		NumAveRealDist_i[loc[0]]++;

		AveP_vg_d[loc[0]] += P_vg_d;
		NumAveP_vg_i[loc[0]]++;

		if (fabs(DistanceToZeroCrossingLoc_mf[i])<20.0) {
			AveDistToZero_d[loc[0]] += (double)fabs(DistanceToZeroCrossingLoc_mf[i]);
			NumAveDistToZero_i[loc[0]]++;
		}
		RealDist_Pvg_Diff_d[loc[0]] += fabs(fabs(P_vg_d) - fabs((ZeroDist_d - Distance_d)));
		NumRealDist_Pvg_Diff_i[loc[0]]++;
		RealDist_DistToZero_Diff_d[loc[0]] += fabs(fabs(DistanceToZeroCrossingLoc_mf[i])-fabs((ZeroDist_d-Distance_d)));
		NumRealDist_DistToZero_Diff_i[loc[0]]++;

		if (fabs(ZeroDist_d - Distance_d)<10.0 && DistanceToZeroCrossingLoc_mf[i]>100.0) {
			printf ("%3d,%3d,%3d, ", DataCoor[0], DataCoor[1], DataCoor[2]);
//			printf ("(%5.1f,", 	(CenterLoc_d[0]-DataCoor[0])*(CenterLoc_d[0]-DataCoor[0]));
//			printf ("%5.1f,", 	(CenterLoc_d[1]-DataCoor[1])*(CenterLoc_d[1]-DataCoor[1]));
//			printf ("%5.1f),", 	(CenterLoc_d[2]-DataCoor[2])*(CenterLoc_d[2]-DataCoor[2]));
			printf ("DistFromC = %7.2f ", Distance_d);
			printf ("RealDistTo0 = %9.4f, ", fabs(ZeroDist_d - Distance_d));
			printf ("CompDistTo0 = %9.5f, ", DistanceToZeroCrossingLoc_mf[i]);
			printf ("P_vg = %9.5f, ", P_vg_d);
			printf ("Diff P_vg = %9.5f, ", fabs(fabs(P_vg_d) - fabs((ZeroDist_d - Distance_d))));
			printf ("Diff DistTo0 = %9.5f, ", fabs(fabs(DistanceToZeroCrossingLoc_mf[i])-fabs((ZeroDist_d-Distance_d))));
			printf ("\n"); fflush (stdout);
		}
		
	}
	
	for (i=0; i<256; i++) {
		if (NumAveRealDist_i[i]>0) AveRealDist_d[i] /= (double)NumAveRealDist_i[i];
		if (NumAveDistToZero_i[i]>0) AveDistToZero_d[i] /= (double)NumAveDistToZero_i[i];
		if (NumAveP_vg_i[i]>0) AveP_vg_d[i] /= (double)NumAveP_vg_i[i];
		if (NumRealDist_Pvg_Diff_i[i]>0) RealDist_Pvg_Diff_d[i] /= (double)NumRealDist_Pvg_Diff_i[i];
		if (NumRealDist_DistToZero_Diff_i[i]>0) RealDist_DistToZero_Diff_d[i] /= (double)NumRealDist_DistToZero_Diff_i[i];
	}
	
	for (i=0; i<256; i++) {
		if (NumAveRealDist_i[i]>0 || NumAveP_vg_i[i]>0 || NumAveDistToZero_i[i]>0 ||
			NumRealDist_Pvg_Diff_i[i]>0 || NumRealDist_DistToZero_Diff_i[i]>0) {
			printf ("Dist = %3d ", i);
			printf ("RealDistTo0 = %9.5f, #%5d ", AveRealDist_d[i], NumAveRealDist_i[i]);
			printf ("CompDistTo0 = %9.5f, #%5d ", AveDistToZero_d[i], NumAveDistToZero_i[i]);
			printf ("P_vg = %9.5f, #%5d ", AveP_vg_d[i], NumAveP_vg_i[i]);
			printf ("Diff P_vg = %9.5f, ", RealDist_Pvg_Diff_d[i]);
			printf ("Diff DistTo0 = %9.5f, #%5d ", RealDist_DistToZero_Diff_d[i], NumRealDist_DistToZero_Diff_i[i]);
			printf ("HitRatio = %5.2f %%", (double)NumAveDistToZero_i[i]/NumAveRealDist_i[i]*100.0);
			printf ("\n"); fflush (stdout);
		}
	}
	
	printf ("\n");
	for (i=0; i<256; i++) {
		if (NumAveRealDist_i[i]>0 || NumAveP_vg_i[i]>0 || NumAveDistToZero_i[i]>0 ||
			NumRealDist_Pvg_Diff_i[i]>0 || NumRealDist_DistToZero_Diff_i[i]>0) {
			printf ("Dist = %3d ", i);
			printf ("RealDistTo0 = %9.5f, #%5d ", AveRealDist_d[i], NumAveRealDist_i[i]);
			printf ("CompDistTo0 = %9.5f, #%5d ", AveDistToZero_d[i], NumAveDistToZero_i[i]);
			printf ("P_vg = %9.5f, #%5d ", AveP_vg_d[i], NumAveP_vg_i[i]);
			printf ("Diff P_vg = %9.5f, ", fabs(fabs(AveRealDist_d[i])-fabs(AveP_vg_d[i])) );
			printf ("Diff DistTo0 = %9.5f, ", fabs(fabs(AveRealDist_d[i])-fabs(AveDistToZero_d[i])));
			printf ("HitRatio = %5.2f %%", (double)NumAveDistToZero_i[i]/NumAveRealDist_i[i]*100.0);
			printf ("\n"); fflush (stdout);
		}
	}
	

}


template<class _DataType>
void cTFGeneration<_DataType>::ComputeHistogramVolume(int NumElements, double& Sigma_Ret)
{
	int		i, NumBinsI, NumBinsG, NumBinsH, MaxNumBins;
	double	HistogramFactorI, HistogramFactorG, HistogramFactorH;
	double	Den_d, FD_d, SD_d;
	int		v_i, g_i, h_i, Freq_i;
	double	MaxFD_d, MaxSD_d, MinSD_d;
	

	HistogramFactorI = (double)NumElements/(MaxData_mf-MinData_mf);
	HistogramFactorG = (double)NumElements/(MaxGrad_mf-MinGrad_mf);
	HistogramFactorH = (double)NumElements/(MaxSecMag_mf-MinSecMag_mf);
	
	NumBinsI = (int)((MaxData_mf-MinData_mf)*HistogramFactorI);
	NumBinsG = (int)((MaxGrad_mf-MinGrad_mf)*HistogramFactorG);
	NumBinsH = (int)((MaxSecMag_mf-MinSecMag_mf)*HistogramFactorH);
	
	MaxNumBins = NumBinsI;
	if (MaxNumBins < NumBinsG) MaxNumBins = NumBinsG;
	if (MaxNumBins < NumBinsH) MaxNumBins = NumBinsH;

	// Kindlmann's Histogram Volume
	HistogramVolume_mi = new int [(MaxNumBins+1)*(MaxNumBins+1)*(MaxNumBins+1)];
	for (i=0; i<(MaxNumBins+1)*(MaxNumBins+1)*(MaxNumBins+1); i++) {
		HistogramVolume_mi[i] = 0;
	}

	printf ("Num Bins I, G, H = %d, %d, %d ", NumBinsI, NumBinsG, NumBinsH);
	printf ("Max Num Bins = %d ", MaxNumBins);
	printf ("\n"); fflush (stdout);

	int	*NumVoxelsHvg_i = new int[(MaxNumBins+1)*(MaxNumBins+1)];
	int	*NumVoxelsGv_i = new int[(MaxNumBins+1)];

	H_vg_md = new double[(MaxNumBins+1)*(MaxNumBins+1)];
	G_v_md = new double[(MaxNumBins+1)];
	
	for (i=0; i<(MaxNumBins+1)*(MaxNumBins+1); i++) {
		H_vg_md[i] = 0.0;
		NumVoxelsHvg_i[i] = (int)0;
	}
	for (i=0; i<(MaxNumBins+1); i++) {
		G_v_md[i] = (double)0.0;
		NumVoxelsGv_i[i] = (int)0;
	}

	MaxFD_d = 0.0;
	MaxSD_d = 0.0;
	MinSD_d = 0.0;
	for (i=0; i<WHD_mi; i++) {
		Den_d = ((double)Data_mT[i]-MinData_mf)*HistogramFactorI;			// Intensity Values
		FD_d = ((double)GradientMag_mf[i]-MinGrad_mf)*HistogramFactorG;		// Gradient Magnitudes
		SD_d = ((double)SecondDerivative_mf[i]-MinSecMag_mf)*HistogramFactorH;// 2nd derivative
		v_i = (int)(Den_d+0.5);
		g_i = (int)(FD_d+0.5);
		h_i = (int)(SD_d+0.5);
		
		HistogramVolume_mi[h_i*(MaxNumBins+1)*(MaxNumBins+1) + g_i*(MaxNumBins+1) + v_i]++;

		if (MaxFD_d < GradientMag_mf[i]) MaxFD_d = GradientMag_mf[i];
		if (MinSD_d > SecondDerivative_mf[i]) MinSD_d = SecondDerivative_mf[i];
		if (MaxSD_d < SecondDerivative_mf[i]) MaxSD_d = SecondDerivative_mf[i];
	}

	Sigma_Ret = MaxFD_d/(sqrt(exp(1.0))*MaxSD_d);
	printf ("\n");
	printf ("Sigma Computation from the Real Min & Max Values\n");
	printf ("exp(1) = %.15f\n", exp(1.0));
	printf ("sqrt(exp(1)) = %.15f\n", sqrt(exp(1.0)));
	printf ("MaxFD_d = %.5f\n", MaxFD_d);
	printf ("MaxSD_d = %.5f\n", MaxSD_d);
	printf ("MinSD_d = %.5f\n", MinSD_d);
	printf ("Sigma_Ret = %.5f\n", Sigma_Ret);


	// Computing H_vg[]
	double	SumHvg_d;
	int		FreqHvg_i;
	for (g_i=0; g_i<(MaxNumBins+1); g_i++) {
		for (v_i=0; v_i<(MaxNumBins+1); v_i++) {
			SumHvg_d = 0;
			FreqHvg_i = 0;
			for (h_i=0; h_i<(MaxNumBins+1); h_i++) {
				Freq_i = HistogramVolume_mi[h_i*(MaxNumBins+1)*(MaxNumBins+1) + g_i*(MaxNumBins+1) + v_i];
				if (Freq_i>0) {
					SumHvg_d += ((double)MinSecMag_mf + (double)h_i/HistogramFactorH)*Freq_i;
					FreqHvg_i += Freq_i;
				}
			}
			if (FreqHvg_i>0) H_vg_md[g_i*(MaxNumBins+1) + v_i] = SumHvg_d/FreqHvg_i;
			else H_vg_md[g_i*(MaxNumBins+1) + v_i] = 0.0;
		}
	}


	// Computing G_v[]
	double	SumGv_d;
	int		FreqGv_i;
	for (v_i=0; v_i<(MaxNumBins+1); v_i++) {
		SumGv_d = 0;
		FreqGv_i = 0;
		for (h_i=0; h_i<(MaxNumBins+1); h_i++) {
			for (g_i=0; g_i<(MaxNumBins+1); g_i++) {
				Freq_i = HistogramVolume_mi[h_i*(MaxNumBins+1)*(MaxNumBins+1) + g_i*(MaxNumBins+1) + v_i];
				if (Freq_i>0) {
					SumGv_d += ((double)MinGrad_mf + (double)g_i/HistogramFactorG)*Freq_i;
					FreqGv_i += Freq_i;
				}
			}
			if (FreqGv_i>0) G_v_md[v_i] = SumGv_d/(double)FreqGv_i;
			else G_v_md[v_i] = 0.0;
		}
	}


	printf ("G_v = \n");
	for (v_i=0; v_i<(MaxNumBins+1); v_i++) {
		printf ("%8.3f ", G_v_md[v_i]);
		printf ("\n"); fflush (stdout);
	}
	printf ("\n");
	fflush (stdout);

	printf ("H_vg = \n");
	for (g_i=0; g_i<(MaxNumBins+1); g_i++) {
		for (v_i=0; v_i<(MaxNumBins+1); v_i++) {
			printf ("%8.3f ", H_vg_md[g_i*(MaxNumBins+1) + v_i]);
		}
		printf ("\n");
		fflush (stdout);
	}
	printf ("\n");
	fflush (stdout);

	int		Sum_i;
	printf ("v g Graph = \n");
	printf ("%d %d\n", MaxNumBins, MaxNumBins);
	for (g_i=0; g_i<MaxNumBins; g_i++) {
		for (v_i=0; v_i<MaxNumBins; v_i++) {
			Sum_i = 0;
			for (h_i=0; h_i<(MaxNumBins+1); h_i++) {
				Sum_i += HistogramVolume_mi[h_i*(MaxNumBins+1)*(MaxNumBins+1) + g_i*(MaxNumBins+1) + v_i];
			}
			printf ("%6d ", Sum_i);
		}
		printf ("\n");
		fflush (stdout);
	}
	printf ("\n");
	fflush (stdout);


	int		j, k;
	int		ZeroHIdx_i = (int)(-MinSecMag_mf*HistogramFactorH+0.5);
	int		*vh_Old_i = new int [(MaxNumBins+1)*(MaxNumBins+1)];
	int		*vh_New_i = new int [(MaxNumBins+1)*(MaxNumBins+1)];
	
	printf ("ZeroHIdx_i = %d\n", ZeroHIdx_i);

	for (i=0; i<(MaxNumBins+1)*(MaxNumBins+1); i++) {
		vh_Old_i[i] = 0;
		vh_New_i[i] = 0;
	}
	for (h_i=0; h_i<(MaxNumBins+1); h_i++) {
		for (v_i=0; v_i<(MaxNumBins+1); v_i++) {
			Sum_i = 0;
			for (g_i=0; g_i<(MaxNumBins+1); g_i++) {
				Sum_i += HistogramVolume_mi[h_i*(MaxNumBins+1)*(MaxNumBins+1) + g_i*(MaxNumBins+1) + v_i];
			}
			vh_Old_i[h_i*(MaxNumBins+1) + v_i] = Sum_i;
		}
	}

	for (j=MaxNumBins/2-1, h_i=ZeroHIdx_i-1; h_i>=0; h_i--, j--) {
		for (k=0, v_i=0; v_i<(MaxNumBins+1); v_i++, k++) {
			if (j<0 || k<0) break;
			vh_New_i[j*(MaxNumBins+1) + k] = vh_Old_i[h_i*(MaxNumBins+1) + v_i];
		}
	}
	for (j=MaxNumBins/2, h_i=ZeroHIdx_i; h_i<(MaxNumBins+1); h_i++, j++) {
		for (k=0, v_i=0; v_i<(MaxNumBins+1); v_i++, k++) {
			if (j>=(MaxNumBins+1) || k>=(MaxNumBins+1)) break;
			vh_New_i[j*(MaxNumBins+1) + k] = vh_Old_i[h_i*(MaxNumBins+1) + v_i];
		}
	}

	
	printf ("v h Graph = \n");
	printf ("%d %d\n", MaxNumBins, MaxNumBins);
	for (h_i=0; h_i<MaxNumBins; h_i++) {
		for (v_i=0; v_i<MaxNumBins; v_i++) {
			printf ("%6d ", vh_New_i[h_i*(MaxNumBins+1) + v_i]);
		}
		printf ("\n");
		fflush (stdout);
	}
	printf ("\n");
	fflush (stdout);


	int		MaxFDIdx_i = 0, MaxSDIdx_i = 0;
	for (h_i=0; h_i<(MaxNumBins+1); h_i++) {
		for (g_i=0; g_i<(MaxNumBins+1); g_i++) {
			for (v_i=0; v_i<(MaxNumBins+1); v_i++) {
				Freq_i = HistogramVolume_mi[h_i*(MaxNumBins+1)*(MaxNumBins+1) + g_i*(MaxNumBins+1) + v_i];
				if (Freq_i>0) {
					if (MaxFDIdx_i < g_i) MaxFDIdx_i = g_i;
					if (MaxSDIdx_i < h_i) MaxSDIdx_i = h_i;
				}
			}
		}
	}

	MaxFD_d = (double)MinGrad_mf + (double)MaxFDIdx_i/HistogramFactorG;
	MaxSD_d = (double)MinSecMag_mf + (double)MaxSDIdx_i/HistogramFactorH;

	Sigma_Ret = MaxFD_d/(sqrt(exp(1.0))*MaxSD_d);
	printf ("\n");
	printf ("Sigma Computation from the Volume Histogram\n");
	printf ("exp(1) = %.15f\n", exp(1.0));
	printf ("sqrt(exp(1)) = %.15f\n", sqrt(exp(1.0)));
	printf ("MaxFD_d = %.5f, MaxFDIdx = %d\n", MaxFD_d, MaxFDIdx_i);
	printf ("MaxSD_d = %.5f, MaxSDIdx = %d\n", MaxSD_d, MaxSDIdx_i);
	printf ("Sigma_Ret = %.5f\n", Sigma_Ret);
	fflush (stdout);


	delete [] vh_Old_i;
	delete [] vh_New_i;
	delete [] NumVoxelsHvg_i;
	delete [] NumVoxelsGv_i;
	
}



template<class _DataType>
void cTFGeneration<_DataType>::Connected_Positive_SecondD_Voxels_RunSave(int ReadFile)
{
	int		i, j, NumVoxels, VoxelIndex, TotalNumCC, CCNum;
	int		ZeroCells_fd, ZeroCellFileExist;
	float	MinCC_f, MaxCC_f, SecondD_f;
	char	CCVolumeFileName[512], ZeroCellFileName[512];
	double	AveGradM;
	

	if (ReadFile) {
		ZeroCrossingCells_mf = new float[WHD_mi];
		for (i=0; i<WHD_mi; i++) ZeroCrossingCells_mf[i] = MinSecMag_mf - 1.0;
		sprintf (ZeroCellFileName, "/work/smpark/data-Engine/C_Engine/Engine_ZeroCells.rawiv");
		if ((ZeroCells_fd = open (ZeroCellFileName, O_RDONLY)) < 0) ZeroCellFileExist = 0;
		else ZeroCellFileExist = 1;

		if (ZeroCellFileExist && ReadFile) {
			cout << "Read the Zero-Crossing Cell File. " <<  ZeroCellFileName << endl;

			if (read(ZeroCells_fd, ZeroCrossingCells_mf, sizeof(float)*WHD_mi) != (unsigned int)sizeof(float)*WHD_mi) {
				cout << "The file could not be read " << ZeroCellFileName << endl;
				close (ZeroCells_fd);
				exit(1);
			}
		}
		else {
			cout << "The Zero-Crossing Cell File does not exist" << endl;
			cout.flush();
			exit(1);
		}
	}
	
	
	ConnectedPositiveSecondDVoxels();


	if (fabsf(MinSecMag_mf)<fabsf(MaxSecMag_mf)) {
		MinCC_f = MinSecMag_mf;
		MaxCC_f = MinSecMag_mf*(-1.0);
	}
	else {
		MinCC_f = MaxSecMag_mf*(-1.0);
		MaxCC_f = MaxSecMag_mf;
	}


	TotalNumCC = 0;
	for (i=0; i<MAX_NUM_CC_STACKS; i++) {
		if (ConnectedVoxel_Stack[i].Size()<100)	ConnectedVoxel_Stack[i].Clear();
		else {
			printf ("CC# = %d, Stack# = %d, Num Voxels = %d\n", TotalNumCC, i, ConnectedVoxel_Stack[i].Size());
			TotalNumCC++;
			}
	}
	
	printf ("Total Num of CC = %d\n", TotalNumCC);
	float	*CCVolume_f = new float [WHD_mi];
	for (i=0; i<WHD_mi; i++) CCVolume_f[i] = MinCC_f;
	
	CCNum=0;
	for (i=0; i<MAX_NUM_CC_STACKS; i++) {
	
		if (!ConnectedVoxel_Stack[i].IsEmpty()) {
	
			NumVoxels = ConnectedVoxel_Stack[i].Size();
			AveGradM = 0.0;
			do {
				ConnectedVoxel_Stack[i].Pop(VoxelIndex);
				SecondD_f = SecondDerivative_mf[VoxelIndex];
				if (SecondD_f<MinCC_f) SecondD_f = MinCC_f;
				if (SecondD_f>MaxCC_f) SecondD_f = MaxCC_f;
				
				CCVolume_f[VoxelIndex] = SecondD_f;
				AveGradM += (double)GradientMag_mf[VoxelIndex];
			} while (!ConnectedVoxel_Stack[i].IsEmpty());
			ConnectedVoxel_Stack[i].Clear();

			AveGradM /= (double)NumVoxels;
			
			if (AveGradM>=5.0) {

				printf ("CC# = %d, Stack# = %d, Num Voxels = %d, ", CCNum, i, NumVoxels);
				printf ("Ave. GradM = %f\n", AveGradM);
				fflush(stdout);
				
				CCVolume_f[WHD_mi-1] = MaxCC_f;
				sprintf (CCVolumeFileName, "CCVolume_%03d", CCNum);
				CCNum++;
				SaveVolume(CCVolume_f, MinCC_f, MaxCC_f, CCVolumeFileName);
			}
			printf ("\n");
			
			for (j=0; j<WHD_mi; j++) CCVolume_f[j] = MinCC_f;
	
		}
		
	}
	delete [] CCVolume_f;


}


template<class _DataType>
void cTFGeneration<_DataType>::ConnectedPositiveSecondDVoxels()
{
	int		i, j, k, l, n, loc[8], ConfigIndex;
	int		VertexI[8];
	float	Cube_f[8];
	

	ConnectedComponentV_mi = new int[WHD_mi];
	for (i=0; i<WHD_mi; i++) ConnectedComponentV_mi[i] = -1;

	for (k=0; k<Depth_mi-1; k++) {
	
		printf ("Connected Positive SecondD Voxels(): Progress Z = %d, ", k);
	
		int		NumStacks = 0;
		for (n=0; n<MAX_NUM_CC_STACKS; n++) {
			if (!ConnectedVoxel_Stack[n].IsEmpty()) NumStacks++;
		}
		printf ("Num Used Stacks = %d\n", NumStacks);
		fflush(stdout);
	
		for (j=0; j<Height_mi-1; j++) {
	
			for (i=0; i<Width_mi-1; i++) {
			
				loc[0] = k*WtimesH_mi + (j+1)*Width_mi + i;		// Vertex Index = 0
				loc[1] = k*WtimesH_mi + (j+1)*Width_mi + i+1;	// Vertex Index = 1
				loc[2] = k*WtimesH_mi + j*Width_mi + i;			// Vertex Index = 2
				loc[3] = k*WtimesH_mi + j*Width_mi + i+1;		// Vertex Index = 3

				loc[4] = (k+1)*WtimesH_mi + (j+1)*Width_mi + i;	// Vertex Index = 4
				loc[5] = (k+1)*WtimesH_mi + (j+1)*Width_mi + i+1;//Vertex Index = 5
				loc[6] = (k+1)*WtimesH_mi + j*Width_mi + i;		// Vertex Index = 6
				loc[7] = (k+1)*WtimesH_mi + j*Width_mi + i+1;	// Vertex Index = 7
				
				ConfigIndex = 0;

				for (n=7; n>=0; n--) {
					ConfigIndex <<= 1;
					Cube_f[n] = ZeroCrossingCells_mf[loc[n]];
					if (ZeroCrossingCells_mf[loc[n]]>=0.0) ConfigIndex |= 1;

//					Cube_f[n] = SecondDerivative_mf[loc[n]];
//					if (SecondDerivative_mf[loc[n]]>=0.0) ConfigIndex |= 1;
				}
				

				if (ConfigIndex==255) continue;

				// C4 Config Index = 24, 36, 66, 129
				if (ConfigIndex==24 || ConfigIndex==36 || ConfigIndex==66 || ConfigIndex==129) {

#ifdef	DEBUG_CC_ZERO_CELLS
//					printf ("C4 \n");
#endif

					if (ConnectedC4(Cube_f)) {
						VertexI[0] = ConnectedVerticesTable_gi[ConfigIndex][0][0];
						VertexI[1] = ConnectedVerticesTable_gi[ConfigIndex][1][0];
						ConnectedVoxels(loc[VertexI[0]], loc[VertexI[1]]);
					}
				}
				else {
				
					// Connected Vertices
					for (l=0; l<7; l++) {
						VertexI[0] = ConnectedVerticesTable_gi[ConfigIndex][l][0];
						VertexI[1] = ConnectedVerticesTable_gi[ConfigIndex][l][1];
						if (VertexI[0]<0) break;
						else {
							ConnectedVoxels(loc[VertexI[0]], loc[VertexI[1]]);
						}
					}
					
					// Ambiguous Faces
					for (l=0; l<6; l++) {
						VertexI[0] = AmbiguousFaceTable_gi[ConfigIndex][l][0];
						VertexI[1] = AmbiguousFaceTable_gi[ConfigIndex][l][1];
						VertexI[2] = AmbiguousFaceTable_gi[ConfigIndex][l][2];
						VertexI[3] = AmbiguousFaceTable_gi[ConfigIndex][l][3];
						if (VertexI[0]<0) break;
						else {
							if (Cube_f[VertexI[0]]*Cube_f[VertexI[1]]>Cube_f[VertexI[2]]*Cube_f[VertexI[3]]) {
								ConnectedVoxels(loc[VertexI[0]], loc[VertexI[1]]);
							}
							else {
								ConnectedVoxels(loc[VertexI[0]], loc[VertexI[0]]);
								ConnectedVoxels(loc[VertexI[1]], loc[VertexI[1]]);
							}
						}
					}
				}
				
		
			}
		}
	}


}


template<class _DataType>
void cTFGeneration<_DataType>::ConnectedVoxels(int VertexI1, int VertexI2)
{
	int		i, CCNum1, CCNum2;


	CCNum1 = ConnectedComponentV_mi[VertexI1];
	CCNum2 = ConnectedComponentV_mi[VertexI2];


	if (CCNum1<0 && CCNum2<0) {
		for (i=0; i<MAX_NUM_CC_STACKS; i++) {
			if (ConnectedVoxel_Stack[i].IsEmpty()) {
				if (VertexI1==VertexI2) {
					ConnectedVoxel_Stack[i].Push(VertexI1);
					ConnectedComponentV_mi[VertexI1] = i;
				}
				else {
					ConnectedVoxel_Stack[i].Push(VertexI1);
					ConnectedVoxel_Stack[i].Push(VertexI2);
					ConnectedComponentV_mi[VertexI1] = i;
					ConnectedComponentV_mi[VertexI2] = i;
				}
				break;
			}
		}
		if (i>=MAX_NUM_CC_STACKS) {
			printf ("Error!: the stack size is too small: ");
			printf ("Stack Size = %d\n", i);
			exit(1);
		}
		return;
	}
	
	if (CCNum1>=0 && CCNum2<0) {
		ConnectedVoxel_Stack[CCNum1].Push(VertexI2);
		ConnectedComponentV_mi[VertexI2] = CCNum1;
		return;
	}
	
	if (CCNum1<0 && CCNum2>=0) {
		ConnectedVoxel_Stack[CCNum2].Push(VertexI1);
		ConnectedComponentV_mi[VertexI1] = CCNum2;
		return;
	}
	
	if (CCNum1>=0 && CCNum2>=0) {
	
		if (VertexI1==VertexI2 || CCNum1==CCNum2) return;

//		printf ("CCNum %d & %d are combined, ", CCNum1, CCNum2);
//		printf ("Size = %d %d\n", ConnectedVoxel_Stack[CCNum1].Size(), ConnectedVoxel_Stack[CCNum2].Size());

		int		VoxelIndex;
		if (ConnectedVoxel_Stack[CCNum1].Size() > ConnectedVoxel_Stack[CCNum2].Size()) {
			do {
				ConnectedVoxel_Stack[CCNum2].Pop(VoxelIndex);
				ConnectedVoxel_Stack[CCNum1].Push(VoxelIndex);
				ConnectedComponentV_mi[VoxelIndex] = CCNum1;
			} while (!ConnectedVoxel_Stack[CCNum2].IsEmpty());
			ConnectedVoxel_Stack[CCNum2].Clear();
		}
		else {
			do {
				ConnectedVoxel_Stack[CCNum1].Pop(VoxelIndex);
				ConnectedVoxel_Stack[CCNum2].Push(VoxelIndex);
				ConnectedComponentV_mi[VoxelIndex] = CCNum2;
			} while (!ConnectedVoxel_Stack[CCNum1].IsEmpty());
			ConnectedVoxel_Stack[CCNum1].Clear();
		}
		return;
	}
	
}

template<class _DataType>
int cTFGeneration<_DataType>::ConnectedC4(float *Cube8)
{
	double	a, b, c, d, e, f, g, h, a_x, a_y, a_z;
	double	DisC_T, G_T;
	

	a = (-Cube8[0] + Cube8[1] + Cube8[2] - Cube8[3] + Cube8[4] - Cube8[5] - Cube8[6] + Cube8[7]); // x*y*z
	b = (+Cube8[0] - Cube8[1] - Cube8[2] + Cube8[3]); // x*y
	c = (+Cube8[0] - Cube8[2] - Cube8[4] + Cube8[6]); // y*z
	d = (+Cube8[0] - Cube8[1] - Cube8[4] + Cube8[5]); // x*z
	e = (-Cube8[0] + Cube8[1]); // x
	f = (-Cube8[0] + Cube8[2]); // y
	g = (-Cube8[0] + Cube8[4]); // z
	h = Cube8[0];

	a_x = (a*e - b*d);
	a_y = (a*f - b*c);
	a_z = (a*g - c*d);
	
	DisC_T = a*h*a*h + b*g*b*g + c*e*c*e + d*f*d*f - 2*a*b*g*h - 2*a*c*e*h
			- 2*a*d*f*h - 2*b*c*e*g - 2*b*d*f*g - 2*c*d*e*f + 4*a*e*f*g + 4*b*c*d*h;
	G_T = a_x * a_y * a_z;

	if (DisC_T<=0 && G_T<0) return false;
	else return true;
}

template<class _DataType>
void cTFGeneration<_DataType>::InitTF(int NumClusters)
{
	int		i, j;

	if (MAX_NUM_MATERIALS <= NumClusters) {
		cout << "cTFGeneration: Error!" << endl;
		cout << "cTFGeneration: The number of clusters is bigger than the give number" << endl;
		exit(1);
	}
	for (i=0; i<NumClusters; i++) TF_mi[i] = new int [256*256*4]; // RGBA
	for (i=0; i<NumClusters; i++) {
		for (j=0; j<256*256*4; j++) TF_mi[i][j]=0;
	}
}


template<class _DataType>
float* cTFGeneration<_DataType>::ComputeSecondD()
{
	int     i, j, k, loc[3];
	double	Dx, Dy, Dz, SDx, SDy, SDz, SecondD;
	double	Tempd;


	float *SecondDerivative = new float [WHD_mi];

	if (Depth_mi==1) { // for 2D
		for (j=0; j<Height_mi; j++) {
			for (k=0; k<Width_mi; k++) {
				if (k==0) SDx = (double)(GradientMag_mf[j*Width_mi + k+1] - 0)/2.0;
				else if (k==Width_mi-1) SDx = (double)(0 - GradientMag_mf[j*Width_mi + k-1])/2.0;
				else SDx = (double)(GradientMag_mf[j*Width_mi + k+1] - GradientMag_mf[j*Width_mi + k-1])/2.0;

				if (j==0) SDy = (double)(GradientMag_mf[(j+1)*Width_mi + k] - 0)/2.0;
				else if (j==Height_mi-1) SDy = (double)(0 - GradientMag_mf[(j-1)*Width_mi + k])/2.0;
				else SDy = (double)(GradientMag_mf[(j+1)*Width_mi + k] - GradientMag_mf[(j-1)*Width_mi + k])/2.0;

				// The Gradient Vector
				Dx = GradientVec_mf[j*Width_mi*3 + k*3 + 0];
				Dy = GradientVec_mf[j*Width_mi*3 + k*3 + 1];
				Dz = GradientVec_mf[j*Width_mi*3 + k*3 + 2];
				Tempd = sqrt(Dx*Dx + Dy*Dy + Dz*Dz);
				if (fabs(Tempd)<1e-6) {
					Dx = 0.0; Dy = 0.0;
				}
				else {
					// Normalize the Gradient Vector
					Dx /= Tempd;  Dy /= Tempd;
				}
				SecondD = SDx*Dx + SDy*Dy; // Dot Product
				SecondDerivative[j*Width_mi + k] = (float)SecondD;
			}
		}
	}
	else { // if (Depth_mi > 1) for 3D
		for (i=0; i<Depth_mi; i++) {
			for (j=0; j<Height_mi; j++) {
				for (k=0; k<Width_mi; k++) {
					loc[0] = i*WtimesH_mi + j*Width_mi + k+1;
					loc[1] = i*WtimesH_mi + j*Width_mi + k-1;
					if (k==0) SDx = (double)(GradientMag_mf[loc[0]] - 0)/2.0;
					else if (k==Width_mi-1) SDx = (double)(0 - GradientMag_mf[loc[1]])/2.0;
					else SDx = (double)(GradientMag_mf[loc[0]] - GradientMag_mf[loc[1]])/2.0;

					loc[0] = i*WtimesH_mi + (j+1)*Width_mi + k;
					loc[1] = i*WtimesH_mi + (j-1)*Width_mi + k;
					if (j==0) SDy = (double)(GradientMag_mf[loc[0]] - 0)/2.0;
					else if (j==Height_mi-1) SDy = (double)(0 - GradientMag_mf[loc[1]])/2.0;
					else SDy = (double)(GradientMag_mf[loc[0]] - GradientMag_mf[loc[1]])/2.0;

					loc[0] = (i+1)*WtimesH_mi + j*Width_mi + k;
					loc[1] = (i-1)*WtimesH_mi + j*Width_mi + k;
					if (i==0) SDz = (double)(GradientMag_mf[loc[0]] - 0)/2.0;
					else if (i==Depth_mi-1) SDz = (double)(0 - GradientMag_mf[loc[1]])/2.0;
					else SDz = (double)(GradientMag_mf[loc[0]] - GradientMag_mf[loc[1]])/2.0;

					// The Gradient Vector
					loc[2] = i*WtimesH_mi*3 + j*Width_mi*3 + k*3;
					Dx = GradientVec_mf[loc[2] + 0];
					Dy = GradientVec_mf[loc[2] + 1];
					Dz = GradientVec_mf[loc[2] + 2];
					Tempd = sqrt(Dx*Dx + Dy*Dy + Dz*Dz);
					if (fabs(Tempd)<1e-6) {
						Dx = 0.0; Dy = 0.0; Dz = 0.0;
					}
					else {
						// Normalize the Gradient Vector
						Dx /= Tempd;  Dy /= Tempd;	Dz /= Tempd;
					}

					SecondD = SDx*Dx + SDy*Dy + SDz*Dz; // Dot Product
					SecondDerivative[i*WtimesH_mi + j*Width_mi + k] = (float)SecondD;
				}
			}
		}
	}
	
	return SecondDerivative;
}


// Directional Second Derivative
template<class _DataType>
void cTFGeneration<_DataType>::ComputeSecondDerivative(char *TargetName)
{
	int		i;
	
		
#ifdef	BILATERAL_FILTER_FOR_SECOND_DERIVATIVE
	double	MinSecondD_d = FLT_MAX, MaxSecondD_d = FLT_MIN;
	double	Sigma_Domain = 50.0, Sigman_Range = 400.0;
	float	*Filtered_SecondD_f;
	char	BilateralFileName[500];
	int		Bilateral_fd1, Bilateral_fd2, SecondDFileExist;

	cout << "Apply Bilateral Filter to the Second Derivatives" << endl;
	sprintf (BilateralFileName, "%s_FilteredSecondD.bilateral", TargetName);
	if ((Bilateral_fd1 = open (BilateralFileName, O_RDONLY)) < 0) SecondDFileExist = 0;
	else SecondDFileExist = 1;
	
	if (SecondDFileExist) {
		cout << "Read the filtered second derivative file. " <<  BilateralFileName << endl;
		Filtered_SecondD_f = new float[WHD_mi];
		if (read(Bilateral_fd1, Filtered_SecondD_f, sizeof(float)*WHD_mi) != (unsigned int)sizeof(float)*WHD_mi) {
			cout << "The file could not be read " << BilateralFileName << endl;
			close (Bilateral_fd1);
			exit(1);
		}
		SecondDerivative_mf = Filtered_SecondD_f;
		Filtered_SecondD_f = NULL;
	}
	else {

		SecondDerivative_mf = ComputeSecondD();

		cout << "Create a New One: " <<  BilateralFileName << endl;
		Filtered_SecondD_f = Bilateral_Filter(SecondDerivative_mf, 5, MinSecondD_d, MaxSecondD_d, Sigma_Domain, Sigman_Range);
		delete [] SecondDerivative_mf;
		SecondDerivative_mf = Filtered_SecondD_f;
		Filtered_SecondD_f = NULL;

		if ((Bilateral_fd2 = open (BilateralFileName, O_CREAT | O_WRONLY)) < 0) {
			cout << "could not open " << BilateralFileName << endl;
			exit(1);
		}
		if (write(Bilateral_fd2, SecondDerivative_mf, sizeof(float)*WHD_mi)!=
			(unsigned int)sizeof(float)*WHD_mi) {
			cout << "The file could not be written " << BilateralFileName << endl;
			close (Bilateral_fd2);
			exit(1);
		}
		if (chmod(BilateralFileName, S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH) < 0) {
			cout << "chmod was not worked to file " << BilateralFileName << endl;
			exit(1);
		}
	}

#endif



#ifdef	ANISOTROPIC_DIFFUSION_FOR_SECONDD
	char	AnisotropicFileName[500];
	int		Anisotropic_fd1, SecondDFileExist;
	int		NumDiffusionRepeat = 3;

	cout << "Apply Anisotropic Diffusion to the Second Derivatives" << endl;
	sprintf (AnisotropicFileName, "%s_DiffusionR%02d.second", TargetName, NumDiffusionRepeat);
	if ((Anisotropic_fd1 = open (AnisotropicFileName, O_RDONLY)) < 0) SecondDFileExist = 0;
	else SecondDFileExist = 1;
	
	if (SecondDFileExist) {
		cout << "Read the filtered second derivative file. " <<  AnisotropicFileName << endl;
		delete [] SecondDerivative_mf;
		SecondDerivative_mf = new float[WHD_mi];
		if (read(Anisotropic_fd1, SecondDerivative_mf, sizeof(float)*WHD_mi) != (unsigned int)sizeof(float)*WHD_mi) {
			cout << "The file could not be read " << AnisotropicFileName << endl;
			close (Anisotropic_fd1);
			exit(1);
		}
	}
	else {
		SecondDerivative_mf = ComputeSecondD();

		cout << "Create a New One: " <<  AnisotropicFileName << endl;
		Anisotropic_Diffusion_3DScalar(SecondDerivative_mf, NumDiffusionRepeat);
		
#ifdef	SAVE_VOLUME_SecondD
		int		Anisotropic_fd2;
		if ((Anisotropic_fd2 = open (AnisotropicFileName, O_CREAT | O_WRONLY)) < 0) {
			cout << "could not open " << AnisotropicFileName << endl;
			exit(1);
		}
		if (write(Anisotropic_fd2, SecondDerivative_mf, sizeof(float)*WHD_mi)!=
			(unsigned int)sizeof(float)*WHD_mi) {
			cout << "The file could not be written " << AnisotropicFileName << endl;
			close (Anisotropic_fd2); exit(1);
		}
		if (chmod(AnisotropicFileName, S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH) < 0) {
			cout << "chmod was not worked to file " << AnisotropicFileName << endl;
			exit(1);
		}
#endif
	}

#endif


#ifndef	BILATERAL_FILTER_FOR_SECOND_DERIVATIVE
#ifndef ANISOTROPIC_DIFFUSION_FOR_SECONDD
	SecondDerivative_mf = ComputeSecondD();
#endif	
#endif	

	
	MinSecMag_mf = 99999.99999;
	MaxSecMag_mf = -99999.99999;
	for (i=0; i<WHD_mi; i++) {
		if (SecondDerivative_mf[i]<MinSecMag_mf) MinSecMag_mf = SecondDerivative_mf[i];
		if (SecondDerivative_mf[i]>MaxSecMag_mf) MaxSecMag_mf = SecondDerivative_mf[i];
	}
	
	cout << "Computing the second derivative is done" << endl;
	cout << "Min & Max Second Derivative = " << MinSecMag_mf << ", " << MaxSecMag_mf << endl;
	
}




// Trilinear Interpolation
//Vxyz =  V000 (1 - Vx) (1 - Vy) (1 - Vz) +
//		V100 Vx (1 - Vy) (1 - Vz) + 
//		V010 (1 - Vx) Vy (1 - Vz) + 
//		V110 Vx Vy (1 - Vz) + 
//		V001 (1 - Vx) (1 - Vy) Vz +
//		V101 Vx (1 - Vy) Vz + 
//		V011 (1 - Vx) Vy Vz + 
//		V111 Vx Vy Vz  

// The vector (Vx, Vy, Vz) should have unit length or 1.
// Tri-linear Interpolation
// (Vx, Vy, Vz) is a location in the volume dataset
template<class _DataType>
double cTFGeneration<_DataType>::GradientInterpolation(double* LocXYZ)
{
	return GradientInterpolation(LocXYZ[0], LocXYZ[1], LocXYZ[2]);
}

template<class _DataType>
double cTFGeneration<_DataType>::GradientInterpolation2(double LocX, double LocY, double LocZ)
{
	int		i, j, k, loc[8];
	double	Vx, Vy, Vz, RetGradM;
	float	Cubef[8];


	i = (int)floor(LocX+1e-8);
	j = (int)floor(LocY+1e-8);
	k = (int)floor(LocZ+1e-8);


	loc[0] = k*WtimesH_mi + (j+1)*Width_mi + i;		// Vertex Index = 0
	loc[1] = k*WtimesH_mi + (j+1)*Width_mi + i+1;	// Vertex Index = 1
	loc[2] = k*WtimesH_mi + j*Width_mi + i;			// Vertex Index = 2
	loc[3] = k*WtimesH_mi + j*Width_mi + i+1;		// Vertex Index = 3

	loc[4] = (k+1)*WtimesH_mi + (j+1)*Width_mi + i;	// Vertex Index = 4
	loc[5] = (k+1)*WtimesH_mi + (j+1)*Width_mi + i+1;//Vertex Index = 5
	loc[6] = (k+1)*WtimesH_mi + j*Width_mi + i;		// Vertex Index = 6
	loc[7] = (k+1)*WtimesH_mi + j*Width_mi + i+1;	// Vertex Index = 7

	Vx = LocX - (double)i;
	Vy = 1.0 - (LocY - (double)j);
	Vz = LocZ - (double)k;


	if (i>=Width_mi || j>=Height_mi || k>=Depth_mi) {
		for (i=0; i<8; i++) loc[i] = 0;
	}
	else {
		if (i+1>=Width_mi) {
			loc[1] = 0; loc[3] = 0; loc[5] = 0; loc[7] = 0; 
		}
		if (j+1>=Height_mi) {
			loc[0] = 0; loc[1] = 0; loc[4] = 0; loc[5] = 0; 
		}
		if (k+1>=Depth_mi) {
			loc[4] = 0; loc[5] = 0; loc[6] = 0; loc[7] = 0; 
		}
	}


	for (i=0; i<8; i++) {
		if (loc[i]==0) Cubef[i] = 0.0;
		else Cubef[i] = GradientMag_mf[loc[i]];
	}


	double	a, b, c, d, e, f, g, h;

	a = (-Cubef[0] + Cubef[1] + Cubef[2] - Cubef[3] + Cubef[4] - Cubef[5] - Cubef[6] + Cubef[7]); // x*y*z
	b = (+Cubef[0] - Cubef[1] - Cubef[2] + Cubef[3]); // x*y
	c = (+Cubef[0] - Cubef[2] - Cubef[4] + Cubef[6]); // y*z
	d = (+Cubef[0] - Cubef[1] - Cubef[4] + Cubef[5]); // x*z
	e = (-Cubef[0] + Cubef[1]); // x
	f = (-Cubef[0] + Cubef[2]); // y
	g = (-Cubef[0] + Cubef[4]); // z
	h = Cubef[0];

	RetGradM = a*Vx*Vy*Vz + b*Vx*Vy + c*Vy*Vz + d*Vz*Vx + e*Vx + f*Vy + g*Vz + h;

	return RetGradM;

}


template<class _DataType>
double cTFGeneration<_DataType>::GradientInterpolation(double LocX, double LocY, double LocZ)
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
int cTFGeneration<_DataType>::GradVecInterpolation(double* LocXYZ, double* GradVec_Ret)
{
	return GradVecInterpolation(LocXYZ[0], LocXYZ[1], LocXYZ[2], GradVec_Ret);
}

template<class _DataType>
int cTFGeneration<_DataType>::GradVecInterpolation2(double LocX, double LocY, double LocZ, double* GradVec_Ret)
{
	int		i, j, k, loc[8];
	double	Vx, Vy, Vz;
	float	Cubef[8];


	i = (int)floor(LocX+1e-8);
	j = (int)floor(LocY+1e-8);
	k = (int)floor(LocZ+1e-8);

	loc[0] = k*WtimesH_mi + (j+1)*Width_mi + i;		// Vertex Index = 0
	loc[1] = k*WtimesH_mi + (j+1)*Width_mi + i+1;	// Vertex Index = 1
	loc[2] = k*WtimesH_mi + j*Width_mi + i;			// Vertex Index = 2
	loc[3] = k*WtimesH_mi + j*Width_mi + i+1;		// Vertex Index = 3

	loc[4] = (k+1)*WtimesH_mi + (j+1)*Width_mi + i;	// Vertex Index = 4
	loc[5] = (k+1)*WtimesH_mi + (j+1)*Width_mi + i+1;//Vertex Index = 5
	loc[6] = (k+1)*WtimesH_mi + j*Width_mi + i;		// Vertex Index = 6
	loc[7] = (k+1)*WtimesH_mi + j*Width_mi + i+1;	// Vertex Index = 7

	Vx = LocX - (double)i;
	Vy = 1.0 - (LocY - (double)j);
	Vz = LocZ - (double)k;

	if (i>=Width_mi || j>=Height_mi || k>=Depth_mi) {
		for (i=0; i<8; i++) loc[i] = 0;
	}
	else {
		if (i+1>=Width_mi) {
			loc[1] = 0; loc[3] = 0; loc[5] = 0; loc[7] = 0; 
		}
		if (j+1>=Height_mi) {
			loc[0] = 0; loc[1] = 0; loc[4] = 0; loc[5] = 0; 
		}
		if (k+1>=Depth_mi) {
			loc[4] = 0; loc[5] = 0; loc[6] = 0; loc[7] = 0; 
		}
	}

	for (i=0; i<8; i++) {
		if (loc[0]==0) Cubef[i] = 0.0;
		else Cubef[i] = GradientMag_mf[loc[i]];
	}


	double	a, b, c, d, e, f, g;

	a = (-Cubef[0] + Cubef[1] + Cubef[2] - Cubef[3] + Cubef[4] - Cubef[5] - Cubef[6] + Cubef[7]); // x*y*z
	b = (+Cubef[0] - Cubef[1] - Cubef[2] + Cubef[3]); // x*y
	c = (+Cubef[0] - Cubef[2] - Cubef[4] + Cubef[6]); // y*z
	d = (+Cubef[0] - Cubef[1] - Cubef[4] + Cubef[5]); // x*z
	e = (-Cubef[0] + Cubef[1]); // x
	f = (-Cubef[0] + Cubef[2]); // y
	g = (-Cubef[0] + Cubef[4]); // z

	// Partial Derivative of Gradient Magnitude
	GradVec_Ret[0] = a*Vy*Vz + b*Vy + d*Vz + e;
	GradVec_Ret[1] = a*Vx*Vz + b*Vx + c*Vz + f;
	GradVec_Ret[2] = a*Vx*Vy + c*Vy + d*Vx + g;

	return true;
}


template<class _DataType>
int cTFGeneration<_DataType>::GradVecInterpolation(double LocX, double LocY, double LocZ, double* GradVec_Ret)
{
	int		i, j, k, loc[3], X, Y, Z;
	double	RetVec[3], GradVec[8*3], Vx, Vy, Vz, Weight_d;



	X = (int)floor(LocX+1e-8);
	Y = (int)floor(LocY+1e-8);
	Z = (int)floor(LocZ+1e-8);
	Vx = LocX - (double)X;
	Vy = LocY - (double)Y;
	Vz = LocZ - (double)Z;

	if (LocX<0.0 || LocX>=(double)Width_mi) return -1;
	if (LocY<0.0 || LocY>=(double)Height_mi) return -1;
	if (LocZ<0.0 || LocZ>=(double)Depth_mi) return -1;

	for (i=0; i<8*3; i++) GradVec[i] = 0.0;
	loc[1] = 0;
	for (k=Z; k<=Z+1; k++) {
		for (j=Y; j<=Y+1; j++) {
			for (i=X; i<=X+1; i++) {
				if (i<0 || j<0 || k<0 || i>=Width_mi || j>=Height_mi ||  k>=Depth_mi) loc[1]++;
				else {
					loc[0] = (k*WtimesH_mi + j*Width_mi + i)*3;
					GradVec[loc[1]*3 + 0] = (double)GradientVec_mf[loc[0] + 0];
					GradVec[loc[1]*3 + 1] = (double)GradientVec_mf[loc[0] + 1];
					GradVec[loc[1]*3 + 2] = (double)GradientVec_mf[loc[0] + 2];
					loc[1]++;
				}
			}
		}
	}

	loc[1] = 0;
	for (k=0; k<3; k++) RetVec[k] = 0.0;
	for (k=0; k<=1; k++) {
		for (j=0; j<=1; j++) {
			for (i=0; i<=1; i++) {
				Weight_d =	((double)(1-i) - Vx*pow((double)-1.0, (double)i))*
							((double)(1-j) - Vy*pow((double)-1.0, (double)j))*
							((double)(1-k) - Vz*pow((double)-1.0, (double)k));
				RetVec[0] += GradVec[loc[1]*3+0]*Weight_d;
				RetVec[1] += GradVec[loc[1]*3+1]*Weight_d;
				RetVec[2] += GradVec[loc[1]*3+2]*Weight_d;
				loc[1] ++;
			}
		}
	}

	for (k=0; k<3; k++) GradVec_Ret[k] = RetVec[k];

	return true;
}



// Tri-linear Interpolation
// (Vx, Vy, Vz) is a location in the volume dataset
template<class _DataType>
double cTFGeneration<_DataType>::DataInterpolation(double* LocXYZ)
{
	return DataInterpolation(LocXYZ[0], LocXYZ[1], LocXYZ[2]);
}

template<class _DataType>
double cTFGeneration<_DataType>::DataInterpolation(double LocX, double LocY, double LocZ)
{
	int		i, j, k, loc[3], X, Y, Z;
	double	RetData, Data[8], Vx, Vy, Vz;




	X = (int)floor(LocX+1e-8);
	Y = (int)floor(LocY+1e-8);
	Z = (int)floor(LocZ+1e-8);
	Vx = LocX - (double)X;
	Vy = LocY - (double)Y;
	Vz = LocZ - (double)Z;

	for (i=0; i<8; i++) Data[i] = 0.0;
	loc[1] = 0;
	for (k=Z; k<=Z+1; k++) {
		for (j=Y; j<=Y+1; j++) {
			for (i=X; i<=X+1; i++) {
				if (i<0 || j<0 || k<0 || i>=Width_mi || j>=Height_mi ||  k>=Depth_mi) loc[1]++;
				else {
					loc[0] = k*WtimesH_mi + j*Width_mi + i;
					Data[loc[1]] = (double)Data_mT[loc[0]];
					loc[1]++;
				}
			}
		}
	}

	loc[1] = 0;
	RetData = 0.0;
	for (k=0; k<=1; k++) {
		for (j=0; j<=1; j++) {
			for (i=0; i<=1; i++) {
				RetData += Data[loc[1]]*	((double)(1-i) - Vx*pow((double)-1.0, (double)i))*
											((double)(1-j) - Vy*pow((double)-1.0, (double)j))*
											((double)(1-k) - Vz*pow((double)-1.0, (double)k));
				loc[1] ++;
			}
		}
	}

	return RetData;

}




// Tri-linear Interpolation
// (Vx, Vy, Vz) is a location in the volume dataset
template<class _DataType>
double cTFGeneration<_DataType>::SecondDInterpolation(double* LocXYZ)
{
	return SecondDInterpolation(LocXYZ[0], LocXYZ[1], LocXYZ[2]);
}

template<class _DataType>
double cTFGeneration<_DataType>::SecondDInterpolation(double LocX, double LocY, double LocZ)
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
void cTFGeneration<_DataType>::ComputeIGGraph(char *OutFileName, int MaterialNum, int NumElements)
{
	int		i, NumBinsI, NumBinsG, NumBins, MinFreq, MaxFreq;
	int		loc[3], *Histogram, DataValue;
	double	HistogramFactorI, HistogramFactorG;
	

	HistogramFactorI = (double)NumElements/(MaxData_mf-MinData_mf);
	HistogramFactorG = (double)NumElements/(MaxGrad_mf-MinGrad_mf);
	
	NumBinsI = (int)((MaxData_mf-MinData_mf)*HistogramFactorI)+1;
	NumBinsG = (int)((MaxGrad_mf-MinGrad_mf)*HistogramFactorG)+1;
	if (NumBinsI > NumBinsG) NumBins = NumBinsI;
	else NumBins = NumBinsG;

	cout << "Num Bins I & G = " << NumBinsI << " " << NumBinsG << endl;
	cout << "Num Bins = " << NumBins << endl;

	Histogram = new int[NumBins*NumBins];
	for (i=0; i<NumBins*NumBins; i++) Histogram[i]=0;
	for (i=0; i<WHD_mi; i++) {
		loc[0] = (int)(((double)Data_mT[i]-MinData_mf)*HistogramFactorI);		// Intensity Values
		loc[1] = (int)(((double)GradientMag_mf[i]-MinGrad_mf)*HistogramFactorG);// Gradient Magnitudes

		DataValue = (int)(((double)Data_mT[i]-MinData_mf)*HistogramFactorI_mf);	// A Data Value in the Histogram
		loc[2] = DataValue*NumMaterial_mi + MaterialNum;
		if (MaterialProb_mf[loc[2]]>0.1) {
			Histogram[loc[1]*NumBins + loc[0]]++;
		}
	}

	int		*SortFreq = new int [NumBins*NumBins];
	for (i=0; i<NumBins*NumBins; i++) SortFreq[i] = Histogram[i];
	QuickSort(SortFreq, 0, NumBins*NumBins-1);

	loc[0] = (int)((double)NumBins*NumBins*0.0001);
	loc[1] = (int)((double)NumBins*NumBins*0.9999);
	MinFreq = SortFreq[loc[0]];
	MaxFreq = SortFreq[loc[1]];
	cout << "100% Min & 100% Max Freq of 2D Histogram = " << SortFreq[0] << " " << SortFreq[NumBins*NumBins-1] << endl;
	cout << "0.01% Min & 0.01% Max Freq of 2D Histogram = " << MinFreq << " " << MaxFreq << endl;
	
	while (MaxFreq==0 && loc[1]<NumBins*NumBins) {
		MaxFreq = SortFreq[loc[1]++];
	}
	delete [] SortFreq;

	cout << "0.01% Min & New Max (Greater than 0) Freq of 2D Histogram = " << MinFreq << " " << MaxFreq << endl;
	
	char	ImageFile[200];
	sprintf (ImageFile, "%s_Mat%d_IG-Graph-Individual.ppm", OutFileName, MaterialNum);
	SaveIGGraph(ImageFile, Histogram, NumBins, MinFreq, MaxFreq);
}

template<class _DataType>
void cTFGeneration<_DataType>::ComputeIGGraphAllMaterials(char *OutFileName, int NumElements)
{
	int		i, j, k, NumBinsI, NumBinsG, NumBins, MinFreq, MaxFreq;
	int		loc[3], *Histogram[50], DataValue, *HistogramAll;
	double	HistogramFactorI, HistogramFactorG;
	

	HistogramFactorI = (double)NumElements/(MaxData_mf-MinData_mf);
	HistogramFactorG = (double)NumElements/(MaxGrad_mf-MinGrad_mf);
	
	NumBinsI = (int)((MaxData_mf-MinData_mf)*HistogramFactorI)+1;
	NumBinsG = (int)((MaxGrad_mf-MinGrad_mf)*HistogramFactorG)+1;
	if (NumBinsI > NumBinsG) NumBins = NumBinsI;
	else NumBins = NumBinsG;

	cout << "Num Bins I & G = " << NumBinsI << " " << NumBinsG << endl;
	cout << "Num Bins = " << NumBins << endl;

	for (k=0; k<50; k++) Histogram[k] = NULL;
	for (k=0; k<NumMaterial_mi; k++) Histogram[k] = new int[NumBins*NumBins];
	for (k=0; k<NumMaterial_mi; k++) {
		for (i=0; i<NumBins*NumBins; i++) {
			Histogram[k][i]=0;
		}
	}

	HistogramAll = new int[NumBins*NumBins];
	for (i=0; i<NumBins*NumBins; i++) HistogramAll[i]=0;

	for (i=0; i<WHD_mi; i++) {
		loc[0] = (int)(((double)Data_mT[i]-MinData_mf)*HistogramFactorI);		// Intensity Values
		loc[1] = (int)(((double)GradientMag_mf[i]-MinGrad_mf)*HistogramFactorG);// Gradient Magnitudes

		DataValue = (int)(((double)Data_mT[i]-MinData_mf)*HistogramFactorI_mf);	// A Data Value in the Histogram
		for (j=0; j<NumMaterial_mi; j++) {
			loc[2] = DataValue*NumMaterial_mi + j;
			if (MaterialProb_mf[loc[2]]>0.1) {
				Histogram[j][loc[1]*NumBins + loc[0]]++;
			}
		}
		HistogramAll[loc[1]*NumBins + loc[0]]++;
	}

	int		*SortFreq = new int [NumBins*NumBins];
	for (i=0; i<NumBins*NumBins; i++) SortFreq[i] = HistogramAll[i];
	QuickSort(SortFreq, 0, NumBins*NumBins-1);

	double		Lower;
	for (Lower=0.1; Lower>1e-6; Lower*=0.1) {

		loc[0] = (int)((double)NumBins*NumBins*Lower);
		loc[1] = (int)((double)NumBins*NumBins*(1.0-Lower));
		MinFreq = SortFreq[loc[0]];
		MaxFreq = SortFreq[loc[1]];
		cout << "100% Min & 100% Max Freq of 2D Histogram = ";
		cout << SortFreq[0] << " " << SortFreq[NumBins*NumBins-1] << endl;
		cout << Lower << "% Min & Max Freq of 2D Histogram = " << MinFreq << " " << MaxFreq << endl;
		if (MaxFreq<=0) continue;

		char	ImageFile[200];
		for (k=0; k<NumMaterial_mi; k++) {
			sprintf (ImageFile, "%s_Mat%d_IG-Graph%f.ppm", OutFileName, k, Lower);
			SaveIGGraph(ImageFile, Histogram[k], NumBins, MinFreq, MaxFreq);
		}
		sprintf (ImageFile, "%s_AllMat_IG-Graph%f.ppm", OutFileName, Lower);
		SaveIGGraph(ImageFile, HistogramAll, NumBins, MinFreq, MaxFreq);
	}
	
	delete [] SortFreq;
	delete [] HistogramAll;
	for (k=0; k<50; k++) delete [] Histogram[k];
}


template<class _DataType>
void cTFGeneration<_DataType>::ComputeIGGraphGeneral(char *OutFileName, int NumElements)
{

	ComputeIGGraphGeneral(OutFileName, NumElements, MinGrad_mf, MaxGrad_mf);

}

template<class _DataType>
void cTFGeneration<_DataType>::ComputeIGGraphGeneral(char *OutFileName, int NumElements, float MinGM, float MaxGM)
{
	int		i, j, NumBinsI, NumBinsG, NumBins;
	int		loc[3], *HistogramAll;
	double	HistogramFactorI, HistogramFactorG;
	double	DataLoc_d[3], Data_d, GM_d;
	int		DataCoor_i[3];


	HistogramFactorI = (double)NumElements/(MaxData_mf-MinData_mf);
	HistogramFactorG = (double)NumElements/(MaxGM-MinGM);
	
	NumBinsI = (int)((MaxData_mf-MinData_mf)*HistogramFactorI);
	NumBinsG = (int)((MaxGM-MinGM)*HistogramFactorG);
	if (NumBinsI > NumBinsG) NumBins = NumBinsI;
	else NumBins = NumBinsG;

	cout << "Num Bins I & G = " << NumBinsI << " " << NumBinsG << endl;
	cout << "Num Bins = " << NumBins << endl;

	HistogramAll = new int[(NumBins+1)*(NumBins+1)];
	for (i=0; i<NumBins*NumBins; i++) HistogramAll[i]=0;

	for (i=0; i<WHD_mi; i++) {
	
		DataCoor_i[2] = i/WtimesH_mi;
		DataCoor_i[1] = (i - DataCoor_i[2]*WtimesH_mi)/Width_mi;
		DataCoor_i[0] = i % Width_mi;

#ifdef	TWO_DIM_R3_IMAGE
		loc[0] = (int)(((double)Data_mT[i]-MinData_mf)*HistogramFactorI);		// Intensity Values
		loc[1] = (int)(((double)GradientMag_mf[i]-MinGM)*HistogramFactorG);// Gradient Magnitudes
		if (loc[1]<0) loc[1] = 0;
		if (loc[1]>NumBins) loc[1] = NumBins;
		if (DataCoor_i[2]==1)	HistogramAll[loc[1]*NumBins + loc[0]]++;
#else
		DataLoc_d[0] = MIN((double)DataCoor_i[0]+0.5, (double)Width_mi-1);
		DataLoc_d[1] = MIN((double)DataCoor_i[1]+0.5, (double)Height_mi-1);
		DataLoc_d[2] = MIN((double)DataCoor_i[2]+0.5, (double)Depth_mi-1);
		
		Data_d = DataInterpolation(DataLoc_d);
		GM_d = GradientInterpolation(DataLoc_d);
		loc[0] = (int)((Data_d - MinData_mf)*HistogramFactorI);		// Intensity Values
		loc[1] = (int)((GM_d - MinGM)*HistogramFactorG);// Gradient Magnitudes
		if (loc[0]<0) loc[0] = 0;
		if (loc[0]>NumBins) loc[0] = NumBins;
		if (loc[1]<0) loc[1] = 0;
		if (loc[1]>NumBins) loc[1] = NumBins;
		HistogramAll[loc[1]*NumBins + loc[0]]++;
#endif
	}

	printf ("OutFileName = %s\n", OutFileName);
	printf ("IG Graph General with Min & Max Grad Mag. %f, %f\n", MinGM, MaxGM);
	printf ("%d %d\n", NumBins, NumBins);
	for (j=0; j<NumBins; j++) {
		for (i=0; i<NumBins; i++) {
			loc[0] = j*NumBins + i;
			printf ("%5d ", HistogramAll[loc[0]]);
		}
		printf ("\n");
	}
	printf ("\n");
	fflush(stdout);
	
	delete [] HistogramAll;

}


template<class _DataType>
void cTFGeneration<_DataType>::ComputeIGGraphBoundaries(char *OutFileName, int NumElements)
{
	int		i, j, k, NumBinsI, NumBinsG, NumBins, MinFreq, MaxFreq;
	int		loc[3], *Histogram[50], DataValue, *HistogramAll, *HistogramAllBoundaries;
	int		TrueBoundary;
	double	HistogramFactorI, HistogramFactorG;
	

	HistogramFactorI = (double)NumElements/(MaxData_mf-MinData_mf);
	HistogramFactorG = (double)NumElements/(MaxGrad_mf-MinGrad_mf);
	
	NumBinsI = (int)((MaxData_mf-MinData_mf)*HistogramFactorI)+1;
	NumBinsG = (int)((MaxGrad_mf-MinGrad_mf)*HistogramFactorG)+1;
	if (NumBinsI > NumBinsG) NumBins = NumBinsI;
	else NumBins = NumBinsG;

	cout << "Num Bins I & G = " << NumBinsI << " " << NumBinsG << endl;
	cout << "Num Bins = " << NumBins << endl;

	for (k=0; k<50; k++) Histogram[k] = NULL;
	for (k=0; k<NumMaterial_mi; k++) Histogram[k] = new int[NumBins*NumBins];
	for (k=0; k<NumMaterial_mi; k++) {
		for (i=0; i<NumBins*NumBins; i++) {
			Histogram[k][i]=0;
		}
	}

	HistogramAll = new int[NumBins*NumBins];
	HistogramAllBoundaries = new int[NumBins*NumBins];
	for (i=0; i<NumBins*NumBins; i++) HistogramAll[i]=0;
	for (i=0; i<NumBins*NumBins; i++) HistogramAllBoundaries[i]=0;
	
	TrueBoundary = false;
	for (i=0; i<WHD_mi; i++) {
		loc[0] = (int)(((double)Data_mT[i]-MinData_mf)*HistogramFactorI);		// Intensity Values
		loc[1] = (int)(((double)GradientMag_mf[i]-MinGrad_mf)*HistogramFactorG);// Gradient Magnitudes

		DataValue = (int)(((double)Data_mT[i]-MinData_mf)*HistogramFactorI_mf);	// A Data Value in the Histogram
		for (j=0; j<NumMaterial_mi; j++) {
			loc[2] = DataValue*NumMaterial_mi + j;
			if (MaterialProb_mf[loc[2]]>0.1) {
				TrueBoundary = IsMaterialBoundary(i, j); // i = Data Location, j = Material Number
				if (TrueBoundary) Histogram[j][loc[1]*NumBins + loc[0]]++;
			}
		}
		if (TrueBoundary) HistogramAllBoundaries[loc[1]*NumBins + loc[0]]++;
		TrueBoundary = false;
		HistogramAll[loc[1]*NumBins + loc[0]]++;
	}

	int		*SortFreq = new int [NumBins*NumBins];
	for (i=0; i<NumBins*NumBins; i++) SortFreq[i] = HistogramAll[i];
	QuickSort(SortFreq, 0, NumBins*NumBins-1);

	double		Lower;
	for (Lower=0.1; Lower>1e-6; Lower*=0.1) {

		loc[0] = (int)((double)NumBins*NumBins*Lower);
		loc[1] = (int)((double)NumBins*NumBins*(1.0-Lower));
		MinFreq = SortFreq[loc[0]];
		MaxFreq = SortFreq[loc[1]];
		cout << "100% Min & 100% Max Freq of 2D Histogram = ";
		cout << SortFreq[0] << " " << SortFreq[NumBins*NumBins-1] << endl;
		cout << Lower << "% Min & Max Freq of 2D Histogram = " << MinFreq << " " << MaxFreq << endl;

		if (MaxFreq<=0) continue;

		char	ImageFile[200];

		for (k=0; k<NumMaterial_mi; k++) {
			sprintf (ImageFile, "%s_Mat%d_IG-Graph-Boundary%f.ppm", OutFileName, k, Lower);
			SaveIGGraph(ImageFile, Histogram[k], NumBins, MinFreq, MaxFreq);
		}
		sprintf (ImageFile, "%s_AllMat_IG-Graph-Boundary%f.ppm", OutFileName, Lower);
		SaveIGGraph(ImageFile, HistogramAllBoundaries, NumBins, MinFreq, MaxFreq);
	}
	
	delete [] SortFreq;
	delete [] HistogramAll;
	delete [] HistogramAllBoundaries;
	for (k=0; k<50; k++) delete [] Histogram[k];

}


// Graph Between Intensities and Second Derivatives
template<class _DataType>
void cTFGeneration<_DataType>::ComputeIHGraphAllMaterials(char *OutFileName, int NumElements)
{
	int		i, j, k, NumBinsI, NumBinsH, NumBins, MinFreq, MaxFreq;
	int		loc[3], *Histogram[50], DataValue, *HistogramAll;
	double	HistogramFactorI, HistogramFactorH;
	

	HistogramFactorI = (double)NumElements/(MaxData_mf-MinData_mf);
	HistogramFactorH = (double)NumElements/(MaxSecMag_mf-MinSecMag_mf);
	
	NumBinsI = (int)((MaxData_mf-MinData_mf)*HistogramFactorI)+1;
	NumBinsH = (int)((MaxSecMag_mf-MinSecMag_mf)*HistogramFactorH)+1;
	if (NumBinsI > NumBinsH) NumBins = NumBinsI;
	else NumBins = NumBinsH;

	cout << "Num Bins I & H = " << NumBinsI << " " << NumBinsH << endl;
	cout << "Num Bins = " << NumBins << endl;

	for (k=0; k<50; k++) Histogram[k] = NULL;
	for (k=0; k<NumMaterial_mi; k++) Histogram[k] = new int[NumBins*NumBins];
	for (k=0; k<NumMaterial_mi; k++) {
		for (i=0; i<NumBins*NumBins; i++) Histogram[k][i]=0;
	}

	HistogramAll = new int[NumBins*NumBins];
	for (i=0; i<NumBins*NumBins; i++) HistogramAll[i]=0;

	for (i=0; i<WHD_mi; i++) {
		loc[0] = (int)(((double)Data_mT[i]-MinData_mf)*HistogramFactorI);		// Intensity Values
		loc[1] = (int)(((double)SecondDerivative_mf[i]-MinSecMag_mf)*HistogramFactorH);

		DataValue = (int)(((double)Data_mT[i]-MinData_mf)*HistogramFactorI_mf);	// A Data Value in the Histogram

		for (j=0; j<NumMaterial_mi; j++) {
			loc[2] = DataValue*NumMaterial_mi + j;
			if (MaterialProb_mf[loc[2]]>0.1) {
				Histogram[j][loc[1]*NumBins + loc[0]]++;
			}
		}
		HistogramAll[loc[1]*NumBins + loc[0]]++;

	}

	int		*SortFreq = new int [NumBins*NumBins];
	for (i=0; i<NumBins*NumBins; i++) SortFreq[i] = HistogramAll[i];
	QuickSort(SortFreq, 0, NumBins*NumBins-1);

	double		Lower=0.01;
	for (Lower=0.1; Lower>1e-6; Lower*=0.1) {

		loc[0] = (int)((double)NumBins*NumBins*Lower);
		loc[1] = (int)((double)NumBins*NumBins*(1.0-Lower));
		MinFreq = SortFreq[loc[0]];
		MaxFreq = SortFreq[loc[1]];
		cout << "100% Min & 100% Max Freq of 2D Histogram = ";
		cout << SortFreq[0] << " " << SortFreq[NumBins*NumBins-1] << endl;
		cout << Lower << "% Min & Max Freq of 2D Histogram = " << MinFreq << " " << MaxFreq << endl;
		if (MaxFreq<=0) continue;

		char	ImageFile[200];
		for (k=0; k<NumMaterial_mi; k++) {
			sprintf (ImageFile, "%s_Mat%d_IH-Graph%f.ppm", OutFileName, k, Lower);
			SaveIGGraph(ImageFile, Histogram[k], NumBins, MinFreq, MaxFreq);
		}
		sprintf (ImageFile, "%s_AllMat_IH-Graph%f.ppm", OutFileName, Lower);
		SaveIGGraph(ImageFile, HistogramAll, NumBins, MinFreq, MaxFreq);
	}
	
	delete [] SortFreq;
	delete [] HistogramAll;
	for (k=0; k<50; k++) delete [] Histogram[k];
}



template<class _DataType>
void cTFGeneration<_DataType>::SaveIGGraph(char *OutFileName, int *IGHisto, int NumBins, int MinFreq, int MaxFreq)
{
	int		i, j, loc[3], Tempi;
	FILE	*PPM_fs;
	

	cout << "Graph File Name: " << OutFileName << endl;
	PPM_fs = fopen (OutFileName, "wt");
	fprintf (PPM_fs, "P3\n%d %d\n 255\n", NumBins, NumBins);
	for (j=NumBins-1; j>=0; j--) {
		for (i=0; i<NumBins; i++) {
			loc[0] = j*NumBins + i;
			Tempi = (int)(((double)IGHisto[loc[0]]-MinFreq)/((double)MaxFreq-MinFreq)*255.0);
			if (Tempi < 0) Tempi = 0;
			if (Tempi > 255) Tempi = 255;
			fprintf (PPM_fs, "%d %d %d\n", Tempi, Tempi, Tempi);
		}
	}
	fclose (PPM_fs);

}


template<class _DataType>
void cTFGeneration<_DataType>::Save_H_vg(char *FileName, float *Hvg, int NumBins, float MinSD, float MaxSD)
{
	int		i, j, loc[3], Tempi;
	FILE	*PPM_fs;
	

	cout << "Graph File Name: " << FileName << endl;
	PPM_fs = fopen (FileName, "wt");
	fprintf (PPM_fs, "P3\n%d %d\n 255\n", NumBins, NumBins);
	for (j=NumBins-1; j>=0; j--) {
		for (i=0; i<NumBins; i++) {
			loc[0] = j*NumBins + i;
			Tempi = (int)(((double)Hvg[loc[0]]-MinSD)/((double)MaxSD-MinSD)*255.0);
			if (Tempi < 0) Tempi = 0;
			if (Tempi > 255) Tempi = 255;
			fprintf (PPM_fs, "%d %d %d\n", Tempi, Tempi, Tempi);
		}
	}
	fclose (PPM_fs);

}

template<class _DataType>
void cTFGeneration<_DataType>::Save_P_vg_Color(char *FileName, float *Pvg, int NumBins, float Minf, float Maxf)
{
	int		i, j, loc[3], R, B;
	FILE	*PPM_fs;
	

	cout << "Graph File Name: " << FileName << endl;
	PPM_fs = fopen (FileName, "wt");
	fprintf (PPM_fs, "P3\n%d %d\n 255\n", NumBins, NumBins);
	for (j=NumBins-1; j>=0; j--) {
		for (i=0; i<NumBins; i++) {
			loc[0] = j*NumBins + i;
			if (Pvg[loc[0]]>0) {
				R = (int)(((double)Pvg[loc[0]]-0)/((double)Maxf-0)*255.0);
				B = 0;
			}
			else {
				B = (int)((fabs((double)Pvg[loc[0]])-0)/(fabs((double)Minf)-0)*255.0);
				R=0;
			}
				
			if (R < 0) R = 0;
			if (R > 255) R = 255;
			if (B < 0) B = 0;
			if (B > 255) B = 255;
			fprintf (PPM_fs, "%d %d %d\n", R, 255, B);
		}
	}
	fclose (PPM_fs);

}



template <class _DataType>
void cTFGeneration<_DataType>::ComputeH_vg(char *OutFileName, int NumElements)
{
	int			i, loc[3], *Freq_vg;
	double		HistogramFactorI, HistogramFactorG;
	int			DataValue;
	

	delete [] H_vg_mf;
	H_vg_mf = new float[NumElements*NumElements];
	Freq_vg = new int[NumElements*NumElements];

	HistogramFactorI = (double)NumElements/(MaxData_mf-MinData_mf);
	HistogramFactorG = (double)NumElements/(MaxGrad_mf-MinGrad_mf);
	
	cout << "Compute H(v,g) Graph = " <<endl;
	cout << "Num NumElements = " << NumElements << endl;

	for (i=0; i<NumElements*NumElements; i++) {
		H_vg_mf[i] = 0;
		Freq_vg[i] = 0;
	}

	int		TrueBoundary=false;
	int		Mini=99999, Maxi=-99999;
	int		*ZeroSecondD = new int [NumElements*NumElements];
	for (i=0; i<NumElements*NumElements; i++) ZeroSecondD[i] = 0;

	cout << "Compute the sum of the second derivatives" << endl;
	for (i=0; i<WHD_mi; i++) {
		loc[0] = (int)(((double)Data_mT[i]-MinData_mf)*HistogramFactorI);		// Intensity Values
		loc[1] = (int)(((double)GradientMag_mf[i]-MinGrad_mf)*HistogramFactorG);// Gradient Magnitudes
		DataValue = (int)(((double)Data_mT[i]-MinData_mf)*HistogramFactorI_mf);	// A Data Value in the Histogram
		loc[2] = DataValue*NumMaterial_mi + 2;

		if (MaterialProb_mf[loc[2]]>0.1) {
//			TrueBoundary = IsMaterialBoundary(i, 2); // i = Data Location, j = Material Number
			TrueBoundary = IsMaterialBoundary(i, 2, 4); // i = Data Location, j = Material Number
		}

		if (200<loc[0] && loc[0]<250 && 200<loc[1] && loc[1]<250 && TrueBoundary) 
			cout << "2nd D = " << loc[0] << " " << loc[1] << " " << SecondDerivative_mf[i] << endl;

		if (-1.0<=SecondDerivative_mf[i] && SecondDerivative_mf[i]<1.0 && GradientMag_mf[i]>0 && TrueBoundary)
			ZeroSecondD[loc[1]*NumElements + loc[0]]++;
		TrueBoundary = false;


		H_vg_mf[loc[1]*NumElements + loc[0]] += SecondDerivative_mf[i];
		Freq_vg[loc[1]*NumElements + loc[0]] ++;
//		if (i%WtimesH_mi==0) cout << "Depth = " << i/WtimesH_mi << endl;
	}


	cout << "Compute Average of H_vg" << endl;
	for (i=0; i<NumElements*NumElements; i++) {
		if (Freq_vg[i]<=0) H_vg_mf[i] = 0.0;
		else H_vg_mf[i] /= (float)Freq_vg[i];
	}

	if (OutFileName==NULL) { // If OutFileName == NULL, then this function is called by Compute TF()
		delete [] Freq_vg;
		return;
	}


	for (i=0; i<NumElements*NumElements; i++) {
		if (Mini>ZeroSecondD[i]) Mini = ZeroSecondD[i];
		if (Maxi<ZeroSecondD[i]) Maxi = ZeroSecondD[i];
	}
	cout << "Min & Max ZeroSecondD = " << Mini << " " << Maxi << endl;
	char	ImageFile[200];
	sprintf (ImageFile, "%s_Zero2ndD_Mat2-Boundary.ppm", OutFileName);
	SaveIGGraph(ImageFile, ZeroSecondD, NumElements, Mini, 20);

}


// Finding the best boundaries while it minimizes the sum of the second derivatives
template <class _DataType>
void cTFGeneration<_DataType>::MakeMaterialNumVolume(int MatNum)
{
	int			i, loc[2];
	int			DataValue;


	delete [] VoxelStatus_muc;
	VoxelStatus_muc = new unsigned char [WHD_mi];
	for (i=0; i<WHD_mi; i++) VoxelStatus_muc[i] = (unsigned char)0;

	for (i=0; i<WHD_mi; i++) {
		DataValue = (int)(((double)Data_mT[i]-MinData_mf)*HistogramFactorI_mf);	// A Data Value in the Histogram
		loc[0] = DataValue*NumMaterial_mi + MatNum;
		if (MaterialProb_mf[loc[0]]>0.1) {
			VoxelStatus_muc[i] = (unsigned char)1; // Inside of Boundaries
			if (IsMaterialBoundary(i, MatNum)) VoxelStatus_muc[i] = (unsigned char)2; // Material Boundaries
		}
	}
}


 // i = Data Location, j = Material Number
template <class _DataType>
int cTFGeneration<_DataType>::IsMaterialBoundaryInMaterialVolume(int DataLoc, int MaterialNum)
{
	int			i, j, k, loc[3];
	int			XCoor, YCoor, ZCoor;
	
	printf ("MaterialNum = %d\n", MaterialNum);
	ZCoor = DataLoc / WtimesH_mi;
	YCoor = (DataLoc - ZCoor*WtimesH_mi) / Height_mi;
	XCoor = DataLoc % Width_mi;

	// Checking all 26 neighbors, whether at least one of them is a different material
	for (k=ZCoor-1; k<=ZCoor+1; k++) {
		if (k<0 || k>=Depth_mi) continue;
		for (j=YCoor-1; j<=YCoor+1; j++) {
			if (j<0 || j>=Height_mi) continue;
			for (i=XCoor-1; i<=XCoor+1; i++) {
				if (i<0 || i>=Width_mi) continue;
				
				// DataValue = A Data Value in the Histogram
				loc[0] = k*WtimesH_mi + j*Width_mi + i;
				if (VoxelStatus_muc[loc[0]]==(unsigned char)0) return true;
			}
		}
	}
	return false;
}

// Finding the zero second directional derivatives, which have local maximum gradients
template <class _DataType>
void cTFGeneration<_DataType>::FindZeroSecondD(char *OutFileName)
{
	int			i, j, k, loc[3];
	int			DataCoor[3];
	double		GVec[3], GradientDir[3], Step, Tempd;
	double		GradM[3];


	printf ("OutFileName = %s\n", OutFileName);
//	for (j=0; j<=NumMaterial_mi; j++)
	for (j=3; j<=3; j++) {

		for (i=0; i<WHD_mi; i++) {

			for (k=0; k<3; k++) GVec[k] = (double)GradientVec_mf[i*3 + k];
			Tempd = sqrt (GVec[0]*GVec[0] + GVec[1]*GVec[1] + GVec[2]*GVec[2]);
			if (fabs(Tempd)<1e-5) continue;
			for (k=0; k<3; k++) GVec[k] /= Tempd; // Normalize the gradient vector

			DataCoor[2] = i / WtimesH_mi;
			DataCoor[1] = (i - DataCoor[2]*WtimesH_mi) / Height_mi;
			DataCoor[0] = i % Width_mi;

			loc[0] = 0;
			for (Step=-0.2; Step<=0.2;Step+=0.2) {

				for (k=0; k<3; k++) GradientDir[k] = (double)DataCoor[k] + GVec[k]*Step;

				GradM[loc[0]] = GradientInterpolation(GradientDir[0], GradientDir[1], GradientDir[2]);

				Tempd = SecondDInterpolation(GradientDir[0], GradientDir[1], GradientDir[2]);

			}

		}
	}

}


// Marking Boundary Voxels
template <class _DataType>
void cTFGeneration<_DataType>::InitializeVoxelStatus()
{
	if (VoxelStatus_muc==NULL) {
		cout << "VoxelStatus_muc is NULL" << endl;
		exit(1);
	}
	for (int i=0; i<WHD_mi; i++) {
		VoxelStatus_muc[i] = (unsigned char)0;
	}
}


// Marking Boundary Voxels
template <class _DataType>
void cTFGeneration<_DataType>::InitializeAlphaVolume()
{
	if (AlphaVolume_mf==NULL) {
		cout << "AlphaVolume_muc is NULL" << endl;
		exit(1);
	}
	for (int i=0; i<WHD_mi; i++) {
		AlphaVolume_mf[i] = (float)0.0;
	}
}


// Marking Boundary Voxels
template <class _DataType>
void cTFGeneration<_DataType>::MarkingGradDirectedVoxels(_DataType MatMin, _DataType MatMax)
{
	int			i, k, loc[3], DataCoor[3], NextDataCoor[3];
	double		GradVec[3], LocAlongGradDir[3], Tempd;
	double		Step, FirstD_d[3], SecondD_d[3];
	int			NumNeighbors, Neighbor26_i[26];


	map<int, unsigned char> MaterialBoundary_m;
	MaterialBoundary_m.clear();

	int		NumMatBoundary=0;
	for (i=0; i<WHD_mi; i++) {

		// Marking the material inside
		if (MatMin<=Data_mT[i] && Data_mT[i]<=MatMax && VoxelStatus_muc[i]==0) {
			VoxelStatus_muc[i] = MAT_INSIDE_BOUNDARY; // Mark the Material Volume
		}
		if (IsMaterialBoundaryUsingMinMax(i, MatMin, MatMax)) {
			MaterialBoundary_m[i]=(unsigned char)0; // Add it to the map
			NumMatBoundary++;
		}
	}
	cout << endl << "Num Mat Boundary = " << NumMatBoundary << endl;
	

#ifdef	DEBUG_TF
	int		PrintSliceNum = 25;
#endif

	map<int, unsigned char>::iterator	Boundary_it;
	
	while ((int)MaterialBoundary_m.size()>0) {
	

#ifdef	DEBUG_TF
		if ((int)MaterialBoundary_m.size()%10000==0)
			cout << "MaterialBoundary_m.size() = " << (int)MaterialBoundary_m.size() << endl;
#endif

		Boundary_it = MaterialBoundary_m.begin();
		loc[0] = (*Boundary_it).first;
		DataCoor[2] = loc[0]/WtimesH_mi;
		DataCoor[1] = (loc[0] - DataCoor[2]*WtimesH_mi)/Width_mi;
		DataCoor[0] = loc[0] % Width_mi;
		MaterialBoundary_m.erase(loc[0]);
		// If the location is already checked, then skip it.

		if (VoxelStatus_muc[loc[0]]>0) continue;
		

		for (k=0; k<3; k++) GradVec[k] = (double)GradientVec_mf[loc[0]*3 + k];
		Tempd = sqrt (GradVec[0]*GradVec[0] + GradVec[1]*GradVec[1] + GradVec[2]*GradVec[2]);
		if (fabs(Tempd)<1e-5) continue;
		for (k=0; k<3; k++) GradVec[k] /= Tempd; // Normalize the gradient vector

		//--------------------------------------------------------------------------------------
		// The Positive Direction of a Gradient Vector
		//--------------------------------------------------------------------------------------
		for (k=0; k<3; k++) LocAlongGradDir[k] = ((double)DataCoor[k])+ GradVec[k]*(-0.2);
		FirstD_d[0] = GradientInterpolation(LocAlongGradDir);
		SecondD_d[0] = SecondDInterpolation(LocAlongGradDir);
		for (k=0; k<3; k++) LocAlongGradDir[k] = ((double)DataCoor[k]);
		FirstD_d[1] = GradientInterpolation(LocAlongGradDir);
		SecondD_d[1] = SecondDInterpolation(LocAlongGradDir);
		for (Step=0.2; Step<=1.0; Step+=0.2) {
			for (k=0; k<3; k++) LocAlongGradDir[k] = ((double)DataCoor[k]) + GradVec[k]*Step;
			FirstD_d[2] = GradientInterpolation(LocAlongGradDir);
			SecondD_d[2] = SecondDInterpolation(LocAlongGradDir);

			if (FirstD_d[1] < FirstD_d[2]) { // Gradient Magnitudes should increase to be MAT_BOUNDARY
				for (k=0; k<3; k++) NextDataCoor[k] = (int)ceil(LocAlongGradDir[k]);
				if (DataCoor[0]==NextDataCoor[0] && DataCoor[1]==NextDataCoor[1] && DataCoor[2]==NextDataCoor[2]) { }
				else {
					loc[1] = NextDataCoor[2]*WtimesH_mi + NextDataCoor[1]*Width_mi + NextDataCoor[0];
					if (VoxelStatus_muc[loc[1]]==0 || VoxelStatus_muc[loc[1]]==MAT_INSIDE_BOUNDARY) {
						VoxelStatus_muc[loc[1]] = (unsigned char)MAT_BOUNDARY; // Marking the boundary
						MaterialBoundary_m[loc[1]]=(unsigned char)0; // Adding it to the map
					}
				}
			}
			else {
				for (k=0; k<3; k++) LocAlongGradDir[k] = ((double)DataCoor[k]) + GradVec[k]*(Step-0.2);
				for (k=0; k<3; k++) NextDataCoor[k] = (int)floor(LocAlongGradDir[k]);
				loc[1] = NextDataCoor[2]*WtimesH_mi + NextDataCoor[1]*Width_mi + NextDataCoor[0];

				// If it is zero-crossing point, then mark it.
				if (SecondD_d[0]*SecondD_d[2]<0.0) {
					VoxelStatus_muc[loc[1]] = (unsigned char)MAT_ZERO_CROSSING; // Marking the zero-crossing voxel
					if (VoxelStatus_muc[loc[1]]==0 || VoxelStatus_muc[loc[1]]==MAT_INSIDE_BOUNDARY) {
						MaterialBoundary_m[loc[1]]=(unsigned char)0; // Adding it to the map
					}

					NumNeighbors = FindNeighbors8Or26(loc[1], Neighbor26_i);
					for (k=0; k<NumNeighbors; k++) {
						if (VoxelStatus_muc[Neighbor26_i[k]]==0 || 
							VoxelStatus_muc[Neighbor26_i[k]]==MAT_INSIDE_BOUNDARY) {
							MaterialBoundary_m[Neighbor26_i[k]]=(unsigned char)0; // Adding it to the map
						}
					}
				}
				else {
					VoxelStatus_muc[loc[1]] = MAT_OTHER_MATERIALS;
					break;
				}

			}
			for (k=0; k<=1; k++) FirstD_d[k] = FirstD_d[k+1];
			for (k=0; k<=1; k++) SecondD_d[k] = SecondD_d[k+1];
		}

		//--------------------------------------------------------------------------------------
		// The Negative Direction of a Gradient Vector
		//--------------------------------------------------------------------------------------
		for (k=0; k<3; k++) LocAlongGradDir[k] = ((double)DataCoor[k]) + GradVec[k]*0.2;
		FirstD_d[0] = GradientInterpolation(LocAlongGradDir);
		SecondD_d[0] = SecondDInterpolation(LocAlongGradDir);
		for (k=0; k<3; k++) LocAlongGradDir[k] = ((double)DataCoor[k]);
		FirstD_d[1] = GradientInterpolation(LocAlongGradDir);
		SecondD_d[1] = SecondDInterpolation(LocAlongGradDir);
		for (Step=-0.2; Step>=-1.0; Step-=0.2) {
			for (k=0; k<3; k++) LocAlongGradDir[k] = ((double)DataCoor[k]) + GradVec[k]*Step;
			FirstD_d[2] = GradientInterpolation(LocAlongGradDir);
			SecondD_d[2] = SecondDInterpolation(LocAlongGradDir);

			if (FirstD_d[1] < FirstD_d[2]) { // Gradient Magnitudes should increase
				for (k=0; k<3; k++) NextDataCoor[k] = (int)ceil(LocAlongGradDir[k]);
				if (DataCoor[0]==NextDataCoor[0] && DataCoor[1]==NextDataCoor[1] && DataCoor[2]==NextDataCoor[2]) { }
				else {
					loc[1] = NextDataCoor[2]*WtimesH_mi + NextDataCoor[1]*Width_mi + NextDataCoor[0];
					if (VoxelStatus_muc[loc[1]]==0 || VoxelStatus_muc[loc[1]]==MAT_INSIDE_BOUNDARY) {
						VoxelStatus_muc[loc[1]] = (unsigned char)MAT_BOUNDARY; // Marking the Boundary
						MaterialBoundary_m[loc[1]]=(unsigned char)0;  // Adding it to the map
					}
				}
			}
			else {
				for (k=0; k<3; k++) LocAlongGradDir[k] = ((double)DataCoor[k]) + GradVec[k]*(Step+0.2);
				for (k=0; k<3; k++) NextDataCoor[k] = (int)floor(LocAlongGradDir[k]);
				loc[1] = NextDataCoor[2]*WtimesH_mi + NextDataCoor[1]*Width_mi + NextDataCoor[0];

				if (SecondD_d[0]*SecondD_d[2]<0.0) {
					VoxelStatus_muc[loc[1]] = (unsigned char)MAT_ZERO_CROSSING; // Marking the zero-crossing voxel
					if (VoxelStatus_muc[loc[1]]==0) {
						MaterialBoundary_m[loc[1]]=(unsigned char)0; // Adding it to the map
					}
					NumNeighbors = FindNeighbors8Or26(loc[1], Neighbor26_i);
					for (k=0; k<NumNeighbors; k++) {
						if (VoxelStatus_muc[Neighbor26_i[k]]==0 || 
							VoxelStatus_muc[Neighbor26_i[k]]==MAT_INSIDE_BOUNDARY) {
							MaterialBoundary_m[Neighbor26_i[k]]=(unsigned char)0; // Adding it to the map
						}
					}
				}
				else {
					VoxelStatus_muc[loc[1]] = MAT_OTHER_MATERIALS;
					break;
				}
			}
			for (k=0; k<=1; k++) FirstD_d[k] = FirstD_d[k+1];
			for (k=0; k<=1; k++) SecondD_d[k] = SecondD_d[k+1];
		}
	}

}

//-----------------------------------------------------------------------------------------------
// Finding a Zero-Crossing Location along Gradient Direction
// Input:
//	CurrLoc: Start Location double[3] (x, y, z)
//	GradVec: Gradient Vector at CurrLoc. It should be normalized
//
// Output:
// ZeroCrossingLoc_Ret: (double[3])
//		Zero Crossing Location from the given location with "CurrLoc" along gradient direction
// FirstDAtTheLoc_Ret:  (double[3])
//		The First Derivative Value (or a gradient magnitude) at the zero-crossing location.
// DataPosFromZeroCrossingLoc_Ret: (double)
//		Zero-crossing location form the given data location. The given data location is 0
//
// Return:
// if it found a zero-crossing location, then return true. Otherwise return false
//
//-----------------------------------------------------------------------------------------------
template <class _DataType>
int cTFGeneration<_DataType>::FindZeroCrossingLocation(double *CurrLoc, double *GradVec, double *ZeroCrossingLoc_Ret, 
									double& FirstDAtTheLoc_Ret, double& DataPosFromZeroCrossingLoc_Ret)
{
	return FindZeroCrossingLocation(CurrLoc, GradVec, ZeroCrossingLoc_Ret, 
									FirstDAtTheLoc_Ret, DataPosFromZeroCrossingLoc_Ret, (double)0.2);
}


template <class _DataType>
int cTFGeneration<_DataType>::FindZeroCrossingLocation(double *CurrLoc, double *GradVec, double *ZeroCrossingLoc_Ret, 
									double& FirstDAtTheLoc_Ret, double& DataPosFromZeroCrossingLoc_Ret, double Interval)
{
	int			i, k, OutOfVolume, FoundZero_i;
	double		FirstD_d[7], SecondD_d[7], Step_d;
	double		StartPos, EndPos, MiddlePos, LocAlongGradDir[3], Increase_d;


	Increase_d = Interval;

	for (i=0, Step_d=-Increase_d*2.0; fabs(Step_d)<=Increase_d*2.0 + 1e-5; Step_d+=Increase_d, i++) {
		for (k=0; k<3; k++) LocAlongGradDir[k] = CurrLoc[k]+ GradVec[k]*Step_d;
		FirstD_d[i] = GradientInterpolation(LocAlongGradDir);
		SecondD_d[i] = SecondDInterpolation(LocAlongGradDir);
		if (fabs(SecondD_d[i])<1e-5) {
			for (k=0; k<3; k++) ZeroCrossingLoc_Ret[k] = LocAlongGradDir[k];
			FirstDAtTheLoc_Ret = FirstD_d[i];
			DataPosFromZeroCrossingLoc_Ret = Step_d;
			return true;
		}
	}

	FoundZero_i = false;
	if (SecondD_d[2]*SecondD_d[3]<0.0 && SecondD_d[2]*SecondD_d[1]<0.0) {
		FoundZero_i = true;
		if (FirstD_d[1] < FirstD_d[3]) {
			StartPos = 0.0;
			EndPos = Increase_d;
		}
		else {
			StartPos = -Increase_d;
			EndPos = 0.0;
		}
	} else if (SecondD_d[2]*SecondD_d[3]<0.0 && SecondD_d[2]*SecondD_d[1]>0.0) {
		FoundZero_i = true;
		StartPos = 0.0;
		EndPos = Increase_d;
	} else if (SecondD_d[2]*SecondD_d[3]>0.0 && SecondD_d[2]*SecondD_d[1]<0.0) {
		FoundZero_i = true;
		StartPos = -Increase_d;
		EndPos = 0.0;
	}
	else {
		FoundZero_i = false;
	}


	if (FoundZero_i) {
		for (k=0; k<3; k++) LocAlongGradDir[k] = (CurrLoc[k]) + GradVec[k]*StartPos;
		SecondD_d[0] = SecondDInterpolation(LocAlongGradDir);
		for (k=0; k<3; k++) LocAlongGradDir[k] = (CurrLoc[k]) + GradVec[k]*EndPos;
		SecondD_d[2] = SecondDInterpolation(LocAlongGradDir);
		MiddlePos = (StartPos + EndPos)/2.0;
		// Binary Search of the zero-crossing location
		do {
			for (k=0; k<3; k++) LocAlongGradDir[k] = (CurrLoc[k]) + GradVec[k]*MiddlePos;
			SecondD_d[1] = SecondDInterpolation(LocAlongGradDir);
			if (fabs(SecondD_d[1])<1e-5) {
				for (k=0; k<3; k++) ZeroCrossingLoc_Ret[k] = LocAlongGradDir[k];
				FirstDAtTheLoc_Ret = FirstD_d[1];
				DataPosFromZeroCrossingLoc_Ret = MiddlePos;
				
#ifdef	DEBUG_TF_ZC
		printf ("Found Zero Loc 1: ");
		printf ("loc = (%7.2f %7.2f %7.2f) ", ZeroCrossingLoc_Ret[0], ZeroCrossingLoc_Ret[1], ZeroCrossingLoc_Ret[2]);
		printf ("DataPosFromZeroCrossingLoc_Ret = %f ", DataPosFromZeroCrossingLoc_Ret);
		Dist_d = sqrt ( (ZeroCrossingLoc_Ret[0]-64.0)*(ZeroCrossingLoc_Ret[0]-64.0) +
						(ZeroCrossingLoc_Ret[1]-64.0)*(ZeroCrossingLoc_Ret[1]-64.0) +
						(ZeroCrossingLoc_Ret[2]-64.0)*(ZeroCrossingLoc_Ret[2]-64.0) );
		printf ("Real Dist = %f ", Dist_d);
		printf ("\n"); fflush (stdout);
#endif
				return true;
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
		for (k=0; k<3; k++) ZeroCrossingLoc_Ret[k] = LocAlongGradDir[k];
		FirstDAtTheLoc_Ret = GradientInterpolation(LocAlongGradDir);
		DataPosFromZeroCrossingLoc_Ret = MiddlePos;
		
#ifdef	DEBUG_TF_ZC
		printf ("Found Zero Loc 2: ");
		printf ("loc = (%7.2f %7.2f %7.2f) ", ZeroCrossingLoc_Ret[0], ZeroCrossingLoc_Ret[1], ZeroCrossingLoc_Ret[2]);
		printf ("DataPosFromZeroCrossingLoc_Ret = %f ", DataPosFromZeroCrossingLoc_Ret);
		Dist_d = sqrt ( (ZeroCrossingLoc_Ret[0]-64.0)*(ZeroCrossingLoc_Ret[0]-64.0) +
						(ZeroCrossingLoc_Ret[1]-64.0)*(ZeroCrossingLoc_Ret[1]-64.0) +
						(ZeroCrossingLoc_Ret[2]-64.0)*(ZeroCrossingLoc_Ret[2]-64.0) );
		printf ("Real Dist = %f ", Dist_d);
		printf ("\n"); fflush (stdout);
#endif

		return true;
	}

	if ((FirstD_d[0]+FirstD_d[1]) < (FirstD_d[3]+FirstD_d[4])) {
		Increase_d = Interval;
		FirstD_d[0] = FirstD_d[2];
		FirstD_d[1] = FirstD_d[3];
	} 
	else {
		if ((FirstD_d[0]+FirstD_d[1]) > (FirstD_d[3]+FirstD_d[4])) {
			Increase_d = Interval*(-1.0);
			FirstD_d[0] = FirstD_d[2];
			FirstD_d[1] = FirstD_d[1];
		}
		else {
			printf ("CurrLoc = %f %f %f ", CurrLoc[0], CurrLoc[1], CurrLoc[2]);
			printf ("\n"); fflush (stdout);
			printf ("(FirstD_d[0]+FirstD_d[1]) == (FirstD_d[3]+FirstD_d[4]) = true ");
			printf ("%f + %f = %f + %f ", FirstD_d[0], FirstD_d[1], FirstD_d[3], FirstD_d[4]);
			printf ("\n"); fflush (stdout);
			for (i=0; i<7; i++) printf ("%f ", FirstD_d[i]);
			printf ("\n\n"); fflush (stdout);

			Increase_d = Interval;
			FirstD_d[0] = FirstD_d[2];
			FirstD_d[1] = FirstD_d[3];
		}
	}
	

//	int		CurrVoxelLoc_i[3];
//	for (k=0; k<3; k++) CurrVoxelLoc_i[k] = (int)LocAlongGradDir[k];

	int		Runtime_i = 0;
	do {
		for (Step_d=Increase_d*2.0; fabs(Step_d)<=15.0+1e-5; Step_d+=Increase_d) {
			OutOfVolume = false;
			for (k=0; k<3; k++) {
				LocAlongGradDir[k] = CurrLoc[k] + GradVec[k]*Step_d;
				if (LocAlongGradDir[k]<0.0) OutOfVolume = true;
			}
			if (OutOfVolume) break;
			if (LocAlongGradDir[0]>=(double)Width_mi) break;
			if (LocAlongGradDir[1]>=(double)Height_mi) break;
			if (LocAlongGradDir[2]>=(double)Depth_mi) break;
			FirstD_d[2] = GradientInterpolation(LocAlongGradDir);
			SecondD_d[2] = SecondDInterpolation(LocAlongGradDir);

	#ifdef	DEBUG_TF_ZC
		printf ("(1st, 2nd) = (%f %f) ", FirstD_d[2], SecondD_d[2]);
		printf ("loc = (%7.2f %7.2f %7.2f) ", LocAlongGradDir[0], LocAlongGradDir[1], LocAlongGradDir[2]);
		printf ("Step = %f ", Step_d);
		Dist_d = sqrt ( (LocAlongGradDir[0]-64.0)*(LocAlongGradDir[0]-64.0) +
						(LocAlongGradDir[1]-64.0)*(LocAlongGradDir[1]-64.0) +
						(LocAlongGradDir[2]-64.0)*(LocAlongGradDir[2]-64.0) );
		printf ("Dist = %f ", Dist_d);
		printf ("\n"); fflush (stdout);
	#endif

			 // Climbing the gradient magnitude mountain
	//			if (FirstD_d[1] < FirstD_d[2] && SecondD_d[1]*SecondD_d[2]>0.0) {
			if (SecondD_d[1]*SecondD_d[2]>0.0) {
				for (k=0; k<=1; k++) FirstD_d[k] = FirstD_d[k+1];
				for (k=0; k<=1; k++) SecondD_d[k] = SecondD_d[k+1];
				continue; // Keep going to the gradient direction
			}
			else {
				if (fabs(SecondD_d[2])<1e-5) {
					for (k=0; k<3; k++) ZeroCrossingLoc_Ret[k] = LocAlongGradDir[k];
					FirstDAtTheLoc_Ret = FirstD_d[2];
					DataPosFromZeroCrossingLoc_Ret = Step_d;
	#ifdef	DEBUG_TF_ZC
			printf ("Found Zero Loc 3: ");
			printf ("loc = (%7.2f %7.2f %7.2f) ", ZeroCrossingLoc_Ret[0], ZeroCrossingLoc_Ret[1], ZeroCrossingLoc_Ret[2]);
			printf ("DataPosFromZeroCrossingLoc_Ret = %f ", DataPosFromZeroCrossingLoc_Ret);
			Dist_d = sqrt ( (ZeroCrossingLoc_Ret[0]-64.0)*(ZeroCrossingLoc_Ret[0]-64.0) +
							(ZeroCrossingLoc_Ret[1]-64.0)*(ZeroCrossingLoc_Ret[1]-64.0) +
							(ZeroCrossingLoc_Ret[2]-64.0)*(ZeroCrossingLoc_Ret[2]-64.0) );
			printf ("Real Dist = %f ", Dist_d);
			printf ("\n"); fflush (stdout);
	#endif
					return true;
				}
				if (SecondD_d[1]*SecondD_d[2]<0.0) {
					StartPos = Step_d - Increase_d;
					EndPos = Step_d;
					MiddlePos = (StartPos + EndPos)/2.0;
					SecondD_d[0] = SecondD_d[1];
					// Binary Search of the zero-crossing location
					do {
						for (k=0; k<3; k++) LocAlongGradDir[k] = (CurrLoc[k]) + GradVec[k]*MiddlePos;
						FirstD_d[1] = GradientInterpolation(LocAlongGradDir);
						SecondD_d[1] = SecondDInterpolation(LocAlongGradDir);
						if (fabs(SecondD_d[1])<1e-5) {
							for (k=0; k<3; k++) ZeroCrossingLoc_Ret[k] = LocAlongGradDir[k];
							FirstDAtTheLoc_Ret = FirstD_d[1];
							DataPosFromZeroCrossingLoc_Ret = MiddlePos;
	#ifdef	DEBUG_TF_ZC
			printf ("Found Zero Loc 4: ");
			printf ("loc = (%7.2f %7.2f %7.2f) ", ZeroCrossingLoc_Ret[0], ZeroCrossingLoc_Ret[1], ZeroCrossingLoc_Ret[2]);
			printf ("DataPosFromZeroCrossingLoc_Ret = %f ", DataPosFromZeroCrossingLoc_Ret);
			Dist_d = sqrt ( (ZeroCrossingLoc_Ret[0]-64.0)*(ZeroCrossingLoc_Ret[0]-64.0) +
							(ZeroCrossingLoc_Ret[1]-64.0)*(ZeroCrossingLoc_Ret[1]-64.0) +
							(ZeroCrossingLoc_Ret[2]-64.0)*(ZeroCrossingLoc_Ret[2]-64.0) );
			printf ("Real Dist = %f ", Dist_d);
			printf ("\n"); fflush (stdout);
	#endif
							return true;
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
					for (k=0; k<3; k++) ZeroCrossingLoc_Ret[k] = LocAlongGradDir[k];
					FirstDAtTheLoc_Ret = FirstD_d[1];
					DataPosFromZeroCrossingLoc_Ret = MiddlePos;
	#ifdef	DEBUG_TF_ZC
			printf ("Found Zero Loc 5: ");
			printf ("loc = (%7.2f %7.2f %7.2f) ", ZeroCrossingLoc_Ret[0], ZeroCrossingLoc_Ret[1], ZeroCrossingLoc_Ret[2]);
			printf ("DataPosFromZeroCrossingLoc_Ret = %f ", DataPosFromZeroCrossingLoc_Ret);
			Dist_d = sqrt ( (ZeroCrossingLoc_Ret[0]-64.0)*(ZeroCrossingLoc_Ret[0]-64.0) +
							(ZeroCrossingLoc_Ret[1]-64.0)*(ZeroCrossingLoc_Ret[1]-64.0) +
							(ZeroCrossingLoc_Ret[2]-64.0)*(ZeroCrossingLoc_Ret[2]-64.0) );
			printf ("Real Dist = %f ", Dist_d);
			printf ("\n"); fflush (stdout);
	#endif
					return true;
				}

			}
		}
		Increase_d *= (-1.0);
		Runtime_i++;
	} while (Runtime_i<=1);

	return false;
}

template <class _DataType>
int cTFGeneration<_DataType>::FindZeroCrossingLocation2(double *CurrLoc, double *GradVec, double *ZeroCrossingLoc_Ret, 
									double& FirstDAtTheLoc_Ret, double& DataPosFromZeroCrossingLoc_Ret, double Interval)
{
	int			k, OutOfVolume;
	double		FirstD_d[3], SecondD_d[3], Step;
	double		StartPos, EndPos, MiddlePos, LocAlongGradDir[3], Increase_d;


	Increase_d = Interval;
	//--------------------------------------------------------------------------------------
	// The Positive and Negative Direction of a Gradient Vector: Repeat twice
	// for the Positive and Negative Directions
	//--------------------------------------------------------------------------------------
#ifdef	DEBUG_TF_ZC
		printf ("FindZeroCrossingLocation = ");
		printf ("Curr loc = (%7.2f %7.2f %7.2f) ", CurrLoc[0], CurrLoc[1], CurrLoc[2]);
		printf ("GradVec = (%7.4f %7.4f %7.4f) ", GradVec[0], GradVec[1], GradVec[2]);
		double		OrgGradVec[3], Tempd, Dist_d;
		
		for (k=0; k<3; k++) OrgGradVec[k] = 64.0 - CurrLoc[k];
		Tempd = sqrt (OrgGradVec[0]*OrgGradVec[0] + OrgGradVec[1]*OrgGradVec[1] + OrgGradVec[2]*OrgGradVec[2]);
		for (k=0; k<3; k++) OrgGradVec[k] /= Tempd;
		for (k=0; k<3; k++) GradVec[k] = OrgGradVec[k];
		printf ("Org. GradVec = (%7.4f %7.4f %7.4f) ", OrgGradVec[0], OrgGradVec[1], OrgGradVec[2]);
		printf ("\n"); fflush (stdout);
#endif

	do {
		for (k=0; k<3; k++) LocAlongGradDir[k] = CurrLoc[k]+ GradVec[k]*(-Increase_d);
		FirstD_d[0] = GradientInterpolation(LocAlongGradDir);
		SecondD_d[0] = SecondDInterpolation(LocAlongGradDir);
		if (fabs(SecondD_d[0])<1e-5) {
			for (k=0; k<3; k++) ZeroCrossingLoc_Ret[k] = LocAlongGradDir[k];
			FirstDAtTheLoc_Ret = FirstD_d[0];
			DataPosFromZeroCrossingLoc_Ret = (-Increase_d);
			return true;
		}
		
#ifdef	DEBUG_TF_ZC
		printf ("(1st, 2nd) = (%f %f) ", FirstD_d[0], SecondD_d[0]);
		printf ("loc = (%7.2f %7.2f %7.2f) ", LocAlongGradDir[0], LocAlongGradDir[1], LocAlongGradDir[2]);
		printf ("GradVec = (%7.4f %7.4f %7.4f) ", GradVec[0], GradVec[1], GradVec[2]);
		printf ("Step = %f ", -Increase_d);
		Dist_d = sqrt ( (LocAlongGradDir[0]-64.0)*(LocAlongGradDir[0]-64.0) +
						(LocAlongGradDir[1]-64.0)*(LocAlongGradDir[1]-64.0) +
						(LocAlongGradDir[2]-64.0)*(LocAlongGradDir[2]-64.0) );
		printf ("Dist = %f ", Dist_d);
		printf ("\n"); fflush (stdout);
#endif

		for (k=0; k<3; k++) LocAlongGradDir[k] = CurrLoc[k];
		FirstD_d[1] = GradientInterpolation(LocAlongGradDir);
		SecondD_d[1] = SecondDInterpolation(LocAlongGradDir);
		if (fabs(SecondD_d[1])<1e-5) {
			for (k=0; k<3; k++) ZeroCrossingLoc_Ret[k] = LocAlongGradDir[k];
			FirstDAtTheLoc_Ret = FirstD_d[1];
			DataPosFromZeroCrossingLoc_Ret = 0.0;
			return true;
		}

#ifdef	DEBUG_TF_ZC
		printf ("(1st, 2nd) = (%f %f) ", FirstD_d[1], SecondD_d[1]);
		printf ("loc = (%7.2f %7.2f %7.2f) ", LocAlongGradDir[0], LocAlongGradDir[1], LocAlongGradDir[2]);
		printf ("Step = %f ", 0.0);
		Dist_d = sqrt ( (LocAlongGradDir[0]-64.0)*(LocAlongGradDir[0]-64.0) +
						(LocAlongGradDir[1]-64.0)*(LocAlongGradDir[1]-64.0) +
						(LocAlongGradDir[2]-64.0)*(LocAlongGradDir[2]-64.0) );
		printf ("Dist = %f ", Dist_d);
		printf ("\n"); fflush (stdout);
#endif


		if (SecondD_d[0]*SecondD_d[1]<0.0) {
			StartPos = -Increase_d;
			EndPos = 0.0;
			MiddlePos = (StartPos + EndPos)/2.0;
			SecondD_d[2] = SecondD_d[1];
			// Binary Search of the zero-crossing location
			do {
				for (k=0; k<3; k++) LocAlongGradDir[k] = (CurrLoc[k]) + GradVec[k]*MiddlePos;
				FirstD_d[1] = GradientInterpolation(LocAlongGradDir);
				SecondD_d[1] = SecondDInterpolation(LocAlongGradDir);
				if (fabs(SecondD_d[1])<1e-5) {
					for (k=0; k<3; k++) ZeroCrossingLoc_Ret[k] = LocAlongGradDir[k];
					FirstDAtTheLoc_Ret = FirstD_d[1];
					DataPosFromZeroCrossingLoc_Ret = MiddlePos;
					return true;
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
			for (k=0; k<3; k++) ZeroCrossingLoc_Ret[k] = LocAlongGradDir[k];
			FirstDAtTheLoc_Ret = FirstD_d[1];
			DataPosFromZeroCrossingLoc_Ret = MiddlePos;
			return true;
		}


		for (Step=Increase_d; fabs(Step)<=15.0+1e-5; Step+=Increase_d) {
			OutOfVolume = false;
			for (k=0; k<3; k++) {
				LocAlongGradDir[k] = CurrLoc[k] + GradVec[k]*Step;
				if (LocAlongGradDir[k]<0.0) OutOfVolume = true;
			}
			if (OutOfVolume) break;
			if (LocAlongGradDir[0]>=(double)Width_mi) break;
			if (LocAlongGradDir[1]>=(double)Height_mi) break;
			if (LocAlongGradDir[2]>=(double)Depth_mi) break;
			FirstD_d[2] = GradientInterpolation(LocAlongGradDir);
			SecondD_d[2] = SecondDInterpolation(LocAlongGradDir);

#ifdef	DEBUG_TF_ZC
		printf ("(1st, 2nd) = (%f %f) ", FirstD_d[2], SecondD_d[2]);
		printf ("loc = (%7.2f %7.2f %7.2f) ", LocAlongGradDir[0], LocAlongGradDir[1], LocAlongGradDir[2]);
		printf ("Step = %f ", Step);
		Dist_d = sqrt ( (LocAlongGradDir[0]-64.0)*(LocAlongGradDir[0]-64.0) +
						(LocAlongGradDir[1]-64.0)*(LocAlongGradDir[1]-64.0) +
						(LocAlongGradDir[2]-64.0)*(LocAlongGradDir[2]-64.0) );
		printf ("Dist = %f ", Dist_d);
		printf ("\n"); fflush (stdout);
#endif

			 // Climbing the gradient magnitude mountain
//			if (FirstD_d[1] < FirstD_d[2] && SecondD_d[1]*SecondD_d[2]>0.0) {
			if (SecondD_d[1]*SecondD_d[2]>0.0) {
				for (k=0; k<=1; k++) FirstD_d[k] = FirstD_d[k+1];
				for (k=0; k<=1; k++) SecondD_d[k] = SecondD_d[k+1];
				continue; // Keep going to the gradient direction
			}
			else {
				if (fabs(SecondD_d[2])<1e-5) {
					for (k=0; k<3; k++) ZeroCrossingLoc_Ret[k] = LocAlongGradDir[k];
					FirstDAtTheLoc_Ret = FirstD_d[2];
					DataPosFromZeroCrossingLoc_Ret = Step;
					return true;
				}
				if (SecondD_d[1]*SecondD_d[2]<0.0) {
					StartPos = Step - Increase_d;
					EndPos = Step;
					MiddlePos = (StartPos + EndPos)/2.0;
					SecondD_d[0] = SecondD_d[1];
					// Binary Search of the zero-crossing location
					do {
						for (k=0; k<3; k++) LocAlongGradDir[k] = (CurrLoc[k]) + GradVec[k]*MiddlePos;
						FirstD_d[1] = GradientInterpolation(LocAlongGradDir);
						SecondD_d[1] = SecondDInterpolation(LocAlongGradDir);
						if (fabs(SecondD_d[1])<1e-5) {
							for (k=0; k<3; k++) ZeroCrossingLoc_Ret[k] = LocAlongGradDir[k];
							FirstDAtTheLoc_Ret = FirstD_d[1];
							DataPosFromZeroCrossingLoc_Ret = MiddlePos;
							return true;
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
					for (k=0; k<3; k++) ZeroCrossingLoc_Ret[k] = LocAlongGradDir[k];
					FirstDAtTheLoc_Ret = FirstD_d[1];
					DataPosFromZeroCrossingLoc_Ret = MiddlePos;
					return true;
				}

				for (k=0; k<=1; k++) FirstD_d[k] = FirstD_d[k+1];
				for (k=0; k<=1; k++) SecondD_d[k] = SecondD_d[k+1];
//				break;
			}
		}

#ifdef	DEBUG_TF_ZC
		printf ("\n");
#endif
		
		Increase_d -= Interval*2.0;	// Step 1: Increase = 0.2, Step 2: Increase = -0.2, Step 3: Increase = -0.6
	} while (Increase_d >= -(Interval*2.0 + Interval/2.0)); // Repeat twice for positive and negative directions



	int		Found_i;
	double	ZeroLoc_d[3];
	Found_i = SearchAgain_WithVectorAdjustment(CurrLoc, GradVec, ZeroLoc_d);
	if (Found_i) {
		DataPosFromZeroCrossingLoc_Ret = sqrt (	(CurrLoc[0]-ZeroLoc_d[0])*(CurrLoc[0]-ZeroLoc_d[0]) +
												(CurrLoc[1]-ZeroLoc_d[1])*(CurrLoc[1]-ZeroLoc_d[1]) +
												(CurrLoc[2]-ZeroLoc_d[2])*(CurrLoc[2]-ZeroLoc_d[2]));
		for (k=0; k<3; k++) ZeroCrossingLoc_Ret[k] = ZeroLoc_d[k];
		FirstDAtTheLoc_Ret = GradientInterpolation(ZeroLoc_d);
		return true;
	}


	return false;
}


template <class _DataType>
int cTFGeneration<_DataType>::SearchAgain_WithVectorAdjustment(double *CurrLoc, double *GradVec, 
								double *ZeroLoc_Ret)
{
	int			k, OutOfVolume;
	double		FirstD_d[5], SecondD_d[5], Step;
	double		StartPos, EndPos, MiddlePos, LocAlongGradDir[3], Increase_d;
	double		AjdGradVec_d[3], Tempd;
	
	

#ifdef	DEBUG_ZEROLOC_SEARCHAGAIN
		printf ("SearchAgain with VectorAdjustment = ");
		printf ("loc = (%7.2f %7.2f %7.2f) ", CurrLoc[0], CurrLoc[1], CurrLoc[2]);
		printf ("GradVec = %7.4f, %7.4f, %7.4f ", GradVec[0], GradVec[1], GradVec[2]);
		printf ("\n"); fflush (stdout);
#endif
	
	for (k=0; k<3; k++) LocAlongGradDir[k] = CurrLoc[k]+ GradVec[k]*(-0.2);
	FirstD_d[3] = GradientInterpolation(LocAlongGradDir);
	SecondD_d[3] = SecondDInterpolation(LocAlongGradDir);


	for (k=0; k<3; k++) LocAlongGradDir[k] = CurrLoc[k]+ GradVec[k]*(-0.1);
	FirstD_d[0] = GradientInterpolation(LocAlongGradDir);
	SecondD_d[0] = SecondDInterpolation(LocAlongGradDir);

	for (k=0; k<3; k++) LocAlongGradDir[k] = CurrLoc[k];
	FirstD_d[1] = GradientInterpolation(LocAlongGradDir);
	SecondD_d[1] = SecondDInterpolation(LocAlongGradDir);
	
	for (k=0; k<3; k++) LocAlongGradDir[k] = CurrLoc[k]+ GradVec[k]*(0.1);
	FirstD_d[2] = GradientInterpolation(LocAlongGradDir);
	SecondD_d[2] = SecondDInterpolation(LocAlongGradDir);


	for (k=0; k<3; k++) LocAlongGradDir[k] = CurrLoc[k]+ GradVec[k]*(0.2);
	FirstD_d[4] = GradientInterpolation(LocAlongGradDir);
	SecondD_d[4] = SecondDInterpolation(LocAlongGradDir);


	if (FirstD_d[1]<FirstD_d[0]) {
		Increase_d = -0.5;
	}
	else if (FirstD_d[1]<FirstD_d[2]) {
		Increase_d = 0.5;
	}
	else if (FirstD_d[3]+FirstD_d[0] < FirstD_d[2]+FirstD_d[4]) {
		Increase_d = 0.5;
	}
	else if (FirstD_d[3]+FirstD_d[0] > FirstD_d[2]+FirstD_d[4]) {
		Increase_d = -0.5;
	}
	else {
#ifdef	DEBUG_ZEROLOC_SEARCHAGAIN
		printf ("loc = (%7.2f %7.2f %7.2f) ", CurrLoc[0], CurrLoc[1], CurrLoc[2]);
		printf ("FD = %7.4f, %7.4f, %7.4f, %7.4f, %7.4f ", FirstD_d[3], FirstD_d[0], FirstD_d[1], FirstD_d[2], FirstD_d[4]);
		printf ("SD = %7.4f, %7.4f, %7.4f, %7.4f, %7.4f ", SecondD_d[3], SecondD_d[0], SecondD_d[1], SecondD_d[2], SecondD_d[4]);
		printf ("\n"); fflush (stdout);
#endif
		return false;
	}
	
	FirstD_d[0] = FirstD_d[1];
	SecondD_d[0] = SecondD_d[1];

#ifdef	DEBUG_ZEROLOC_SEARCHAGAIN
		printf ("(1st, 2nd) = (%f %f) ", FirstD_d[0], SecondD_d[0]);
		printf ("loc = (%7.2f %7.2f %7.2f) ", CurrLoc[0], CurrLoc[1], CurrLoc[2]);
		printf ("\n"); fflush (stdout);
#endif

	for (k=0; k<3; k++) LocAlongGradDir[k] = CurrLoc[k]+ GradVec[k]*Increase_d;
	FirstD_d[1] = GradientInterpolation(LocAlongGradDir);
	SecondD_d[1] = SecondDInterpolation(LocAlongGradDir);

#ifdef	DEBUG_ZEROLOC_SEARCHAGAIN
	double		Dist_d;
	printf ("(1st, 2nd) = (%f %f) ", FirstD_d[1], SecondD_d[1]);
	printf ("loc = (%7.2f %7.2f %7.2f) ", LocAlongGradDir[0], LocAlongGradDir[1], LocAlongGradDir[2]);
	printf ("Org. GradVec = (%7.4f %7.4f %7.4f) ", GradVec[0], GradVec[1], GradVec[2]);
	printf ("\n"); fflush (stdout);
#endif


	for (Step=Increase_d*2.0; fabs(Step)<=15.0+1e-5; Step+=Increase_d) {

		if (GradVecInterpolation(LocAlongGradDir, AjdGradVec_d)<0) return false;
		Tempd = sqrt (AjdGradVec_d[0]*AjdGradVec_d[0] + AjdGradVec_d[1]*AjdGradVec_d[1] + AjdGradVec_d[2]*AjdGradVec_d[2]);
		for (k=0; k<3; k++) AjdGradVec_d[k] /= Tempd;
		OutOfVolume = false;
		for (k=0; k<3; k++) {
			LocAlongGradDir[k] = CurrLoc[k] + AjdGradVec_d[k]*Step;
			if (LocAlongGradDir[k]<0.0) OutOfVolume = true;
		}
		if (OutOfVolume) break;
		if (LocAlongGradDir[0]>=(double)Width_mi) break;
		if (LocAlongGradDir[1]>=(double)Height_mi) break;
		if (LocAlongGradDir[2]>=(double)Depth_mi) break;
		FirstD_d[2] = GradientInterpolation(LocAlongGradDir);
		SecondD_d[2] = SecondDInterpolation(LocAlongGradDir);

#ifdef	DEBUG_ZEROLOC_SEARCHAGAIN
		printf ("(1st, 2nd) = (%f %f) ", FirstD_d[2], SecondD_d[2]);
		printf ("loc = (%7.2f %7.2f %7.2f) ", LocAlongGradDir[0], LocAlongGradDir[1], LocAlongGradDir[2]);
		printf ("GradVec = (%7.4f %7.4f %7.4f) ", AjdGradVec_d[0], AjdGradVec_d[1], AjdGradVec_d[2]);
		printf ("Step = %f ", Step);
		printf ("\n"); fflush (stdout);
#endif


		 // Climbing the gradient magnitude mountain
		if (FirstD_d[1] < FirstD_d[2] && SecondD_d[1]*SecondD_d[2]>0.0) {
			for (k=0; k<=1; k++) FirstD_d[k] = FirstD_d[k+1];
			for (k=0; k<=1; k++) SecondD_d[k] = SecondD_d[k+1];
			continue; // Keep going to the gradient direction
		}
		else {
			if (fabs(SecondD_d[2])<1e-5) {
				for (k=0; k<3; k++) ZeroLoc_Ret[k] = LocAlongGradDir[k];

#ifdef	DEBUG_ZEROLOC_SEARCHAGAIN
		printf ("Found Zero Loc "); fflush (stdout);
		printf ("loc = (%7.2f %7.2f %7.2f) ", ZeroLoc_Ret[0], ZeroLoc_Ret[1], ZeroLoc_Ret[2]); fflush (stdout);
		Dist_d = sqrt ( (ZeroLoc_Ret[0]-CurrLoc[0])*(ZeroLoc_Ret[0]-CurrLoc[0]) + 
						(ZeroLoc_Ret[1]-CurrLoc[1])*(ZeroLoc_Ret[1]-CurrLoc[1]) + 
						(ZeroLoc_Ret[2]-CurrLoc[2])*(ZeroLoc_Ret[2]-CurrLoc[2]) );
		printf ("Dist = %f ", Dist_d); fflush (stdout);
		printf ("\n"); fflush (stdout);
#endif

				return true;
			}
			if (SecondD_d[1]*SecondD_d[2]<0.0) {
				StartPos = Step - Increase_d;
				EndPos = Step;
				MiddlePos = (StartPos + EndPos)/2.0;
				SecondD_d[0] = SecondD_d[1];
				// Binary Search of the zero-crossing location
				do {
					for (k=0; k<3; k++) LocAlongGradDir[k] = CurrLoc[k] + AjdGradVec_d[k]*MiddlePos;
					FirstD_d[1] = GradientInterpolation(LocAlongGradDir);
					SecondD_d[1] = SecondDInterpolation(LocAlongGradDir);
					if (fabs(SecondD_d[1])<1e-5) {
						for (k=0; k<3; k++) ZeroLoc_Ret[k] = LocAlongGradDir[k];

#ifdef	DEBUG_ZEROLOC_SEARCHAGAIN
		printf ("Found Zero Loc "); fflush (stdout);
		printf ("loc = (%7.2f %7.2f %7.2f) ", ZeroLoc_Ret[0], ZeroLoc_Ret[1], ZeroLoc_Ret[2]); fflush (stdout);
		Dist_d = sqrt ( (ZeroLoc_Ret[0]-CurrLoc[0])*(ZeroLoc_Ret[0]-CurrLoc[0]) + 
						(ZeroLoc_Ret[1]-CurrLoc[1])*(ZeroLoc_Ret[1]-CurrLoc[1]) + 
						(ZeroLoc_Ret[2]-CurrLoc[2])*(ZeroLoc_Ret[2]-CurrLoc[2]) );
		printf ("Dist = %f ", Dist_d); fflush (stdout);
		printf ("\n"); fflush (stdout);
#endif

						return true;
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
				for (k=0; k<3; k++) ZeroLoc_Ret[k] = LocAlongGradDir[k];

#ifdef	DEBUG_ZEROLOC_SEARCHAGAIN
		printf ("Found Zero Loc "); fflush (stdout);
		printf ("loc = (%7.2f %7.2f %7.2f) ", ZeroLoc_Ret[0], ZeroLoc_Ret[1], ZeroLoc_Ret[2]); fflush (stdout);
		Dist_d = sqrt ( (ZeroLoc_Ret[0]-CurrLoc[0])*(ZeroLoc_Ret[0]-CurrLoc[0]) + 
						(ZeroLoc_Ret[1]-CurrLoc[1])*(ZeroLoc_Ret[1]-CurrLoc[1]) + 
						(ZeroLoc_Ret[2]-CurrLoc[2])*(ZeroLoc_Ret[2]-CurrLoc[2]) );
		printf ("Dist = %f ", Dist_d); fflush (stdout);
		printf ("\n"); fflush (stdout);
#endif
				return true;
			}
			break;
		}
	}

	return false;
}



// Marking Boundary Voxels
template <class _DataType>
void cTFGeneration<_DataType>::MarkingZeroCrossingLoc(_DataType MatMin, _DataType MatMax)
{
	int			i, k, loc[3], DataCoor[3], ZeroCrossingLoc_i[3];
	double		GradVec[3], ZeroCrossingLoc_d[3], Tempd;
	double		LocalMaxGradient, DataPosFromZeroCrossingLoc_d;
	int			NumNeighbors, Neighbor26_i[26], FoundZeroCrossingLoc;


	map<int, unsigned char> InitialBoundaryLocs_map;
	map<int, unsigned char> ZeroCrossingLocs_map;
	InitialBoundaryLocs_map.clear();
	ZeroCrossingLocs_map.clear();

	int		NumMatBoundary=0;

	for (i=0; i<WHD_mi; i++) {

//		if (MatMin<=Data_mT[i] && Data_mT[i]<=MatMax && VoxelStatus_muc[i]==0) {
//			VoxelStatus_muc[i] = MAT_INSIDE_BOUNDARY; // Mark the Material Volume
//		}

		if (IsMaterialBoundaryUsingMinMax(i, MatMin, MatMax)) {
			InitialBoundaryLocs_map[i]=(unsigned char)0; // Add it to the map
			NumMatBoundary++;
		}
	}
	cout << endl << "Num Mat Boundary = " << NumMatBoundary << " ";
	cout << "Material Range = " << (int)MatMin << " " << (int)MatMax << endl;


//(105 110 25)
//	i = 25*WtimesH_mi + 110*Width_mi + 105;
//	InitialBoundaryLocs_map[i]=(unsigned char)0; // Add it to the map
	

	map<int, unsigned char>::iterator	Boundary_it;

#ifdef	DEBUG_TF
	int		NumRepeat=0;
#endif

	while ((int)InitialBoundaryLocs_map.size()>0 || (int)ZeroCrossingLocs_map.size()>0) {


#ifdef	DEBUG_TF
		NumRepeat++;
		if (NumRepeat%1000==0) {
			printf ("InitialBoundaryLocs_map_size = %d ", (int)InitialBoundaryLocs_map.size());
			printf ("# Repeat = %d\n", NumRepeat);
			fflush(stdout);
		}
#endif


		if ((int)ZeroCrossingLocs_map.size()>0) {
			Boundary_it = ZeroCrossingLocs_map.begin();
			loc[0] = (*Boundary_it).first;


			ZeroCrossingLocs_map.erase(loc[0]);
		}
		else {
			Boundary_it = InitialBoundaryLocs_map.begin();
			loc[0] = (*Boundary_it).first;
			InitialBoundaryLocs_map.erase(loc[0]);
		}
		DataCoor[2] = loc[0]/WtimesH_mi;
		DataCoor[1] = (loc[0] - DataCoor[2]*WtimesH_mi)/Width_mi;
		DataCoor[0] = loc[0] % Width_mi;

		// If the location is already checked, then skip it.
		// This is the same as if (VoxelStatus_muc[loc[0]]>0) continue;
		if ((VoxelStatus_muc[loc[0]] & MAT_CHECKED_VOXEL)) continue;

		// Getting the gradient vector at the position
		for (k=0; k<3; k++) GradVec[k] = (double)GradientVec_mf[loc[0]*3 + k];
		Tempd = sqrt (GradVec[0]*GradVec[0] + GradVec[1]*GradVec[1] + GradVec[2]*GradVec[2]);
		if (fabs(Tempd)<1e-6) {

#ifdef	DEBUG_TF
			printf ("Length is 0 at (%d %d %d) %X\n", DataCoor[0], DataCoor[1], DataCoor[2], VoxelStatus_muc[loc[0]]);
#endif
			VoxelStatus_muc[loc[0]] |= MAT_CHECKED_VOXEL;
			continue;
		}
		for (k=0; k<3; k++) GradVec[k] /= Tempd; // Normalize the gradient vector

		//-------------------------------------------------------------------------------------------
		// Finding the local maxima of gradient magnitudes along the gradient direction. 
		// It climbs the mountain of gradient magnitudes to find the zero-crossing second derivative 
		// Return ZeroCrossingLoc_d, LocalMaxGradient, DataPosFromZeroCrossingLoc_d
		//-------------------------------------------------------------------------------------------
		int			FoundSigma, IsCorrectBoundary;
		double		CurrDataLoc_d[3], DataValuePos_d, DataValueNeg_d;
		double		SigmaPosDir_Ret, SigmaNegDir_Ret, SecondDPosDir_Ret, SecondDNegDir_Ret;
		for (k=0; k<3; k++) CurrDataLoc_d[k] = (double)DataCoor[k];


		FoundZeroCrossingLoc = FindZeroCrossingLocation(CurrDataLoc_d, GradVec, ZeroCrossingLoc_d, 
														LocalMaxGradient, DataPosFromZeroCrossingLoc_d);
		


		for (k=0; k<3; k++) ZeroCrossingLoc_i[k] = (int)floor(ZeroCrossingLoc_d[k]);
		loc[1] = ZeroCrossingLoc_i[2]*WtimesH_mi + ZeroCrossingLoc_i[1]*Width_mi + ZeroCrossingLoc_i[0];


#ifdef	DEBUG_TF

		if (FoundZeroCrossingLoc) {
			printf ("DataCoor = (%d %d %d), %X ", DataCoor[0], DataCoor[1], DataCoor[2], VoxelStatus_muc[loc[0]]);
			printf ("Zero2nd = %8.4f %8.4f %8.4f, ", ZeroCrossingLoc_d[0], ZeroCrossingLoc_d[1], ZeroCrossingLoc_d[2]);
			printf ("MaxGradient = %7.2f, ", LocalMaxGradient);
			printf ("DataPosFromZero2nd = %8.4f ", DataPosFromZeroCrossingLoc_d);
			printf ("At Zero = %X\n", VoxelStatus_muc[loc[1]]);
			fflush(stdout);
		}
		else {
			printf ("DataCoor = (%d %d %d), %X ", DataCoor[0], DataCoor[1], DataCoor[2], VoxelStatus_muc[loc[0]]);
			printf ("No Zero Crossing Loc\n");
//			FindZeroCrossingLocation(CurrDataLoc_d, GradVec, ZeroCrossingLoc_d, 
//											LocalMaxGradient, DataPosFromZeroCrossingLoc_d);
		}

#endif

		if ((VoxelStatus_muc[loc[1]] & MAT_ZERO_CROSSING)) continue;

		if (FoundZeroCrossingLoc) {

			IsCorrectBoundary = false;
			double	GradVec_d[3], SamplingInterval_d;
			GradVecInterpolation(ZeroCrossingLoc_d, GradVec_d);
			Tempd = sqrt (GradVec_d[0]*GradVec_d[0] + GradVec_d[1]*GradVec_d[1] + GradVec_d[2]*GradVec_d[2]);
			if (fabs(Tempd)<1e-6) continue;
			for (k=0; k<3; k++) GradVec_d[k] /= Tempd; // Normalize the gradient vector
			

#ifdef	DEBUG_TF_SIGMA
			printf ("Compute Positive Direction Sigma\n");
#endif

			FoundSigma = ComputeSigma(ZeroCrossingLoc_d, GradVec_d, SigmaPosDir_Ret, SecondDPosDir_Ret, 0.1);
			if (FoundSigma) {
				// 2.5759 = 99.00%, 3.89 = 99.99%
				for (SamplingInterval_d=1.0; SamplingInterval_d<3.89; SamplingInterval_d+=0.5) {
					DataValuePos_d = getDataValue(ZeroCrossingLoc_d, GradVec_d, SigmaPosDir_Ret*SamplingInterval_d);
					if (((double)MatMin-1e-5) <= DataValuePos_d && 
						((double)MatMax+1e-5) >= DataValuePos_d ) IsCorrectBoundary = true;
				}
			}

#ifdef	DEBUG_TF_SIGMA
			printf ("Compute Negative Direction Sigma\n");
#endif


			FoundSigma = ComputeSigma(ZeroCrossingLoc_d, GradVec_d, SigmaNegDir_Ret, SecondDNegDir_Ret, -0.1);
			if (FoundSigma) {
				// 2.5759 = 99.00%, 3.89 = 99.99%
				for (SamplingInterval_d=1.0; SamplingInterval_d<3.89; SamplingInterval_d+=0.5) {
					DataValueNeg_d = getDataValue(ZeroCrossingLoc_d, GradVec_d, SigmaNegDir_Ret*SamplingInterval_d);
					if (((double)MatMin-1e-5) <= DataValueNeg_d && 
						((double)MatMax+1e-5) >= DataValueNeg_d ) IsCorrectBoundary = true;
				}
			}
			

#ifdef	DEBUG_TF_SIGMA
			printf ("\n");
			fflush(stdout);
#endif


			if (!IsCorrectBoundary) {

#ifdef	DEBUG_TF
				printf ("Incorrect Boundary: Found Sigma = %d  ", FoundSigma);
				printf ("DataPos = %f SigmaPos = %f, ", DataValuePos_d, SigmaPosDir_Ret);
				printf ("(%f ", ZeroCrossingLoc_d[0] + GradVec_d[0]*SigmaPosDir_Ret);
				printf ("%f ",  ZeroCrossingLoc_d[1] + GradVec_d[1]*SigmaPosDir_Ret);
				printf ("%f) ", ZeroCrossingLoc_d[2] + GradVec_d[2]*SigmaPosDir_Ret);
				printf ("DataNeg = %f SigmaNeg = %f, ", DataValueNeg_d, SigmaNegDir_Ret);
				printf ("(%f ", ZeroCrossingLoc_d[0] + GradVec_d[0]*SigmaNegDir_Ret);
				printf ("%f ",  ZeroCrossingLoc_d[1] + GradVec_d[1]*SigmaNegDir_Ret);
				printf ("%f) ", ZeroCrossingLoc_d[2] + GradVec_d[2]*SigmaNegDir_Ret);
				printf ("GradVec = ");
				printf ("(%f ",    GradVec_d[0]);
				printf ("%f ",     GradVec_d[1]);
				printf ("%f)\n\n", GradVec_d[2]);
				fflush (stdout);
#endif
				VoxelStatus_muc[loc[0]] |= MAT_CHECKED_VOXEL; // Marking the Checked Voxel
				continue;
			}
			else {
				VoxelStatus_muc[loc[1]] |= MAT_ZERO_CROSSING;
				NumNeighbors = FindZeroConnectedNeighbors(loc[1], Neighbor26_i);
				for (k=0; k<NumNeighbors; k++) {
					if (VoxelStatus_muc[Neighbor26_i[k]] == 0) {
						ZeroCrossingLocs_map[Neighbor26_i[k]]=(unsigned char)0; // Add it to the map
					}
				}

#ifdef	DEBUG_TF

				if (NumNeighbors > 0) {
					printf ("NumNeighbors = %d   ", NumNeighbors);
					printf ("ZeroCrossing = (%d %d %d) ", ZeroCrossingLoc_i[0], ZeroCrossingLoc_i[1], ZeroCrossingLoc_i[2]);
					printf ("%x ", VoxelStatus_muc[loc[1]]);
					printf ("Neigbors = ");
					for (k=0; k<NumNeighbors; k++) {
						if (VoxelStatus_muc[Neighbor26_i[k]] == 0) {
							int Z = Neighbor26_i[k]/WtimesH_mi;
							int Y = (Neighbor26_i[k] - Z*WtimesH_mi)/Width_mi;
							int X = Neighbor26_i[k] % Width_mi;
							printf ("(%d %d %d) ", X, Y, Z);
						}
					}
					printf ("\n\n");
					fflush(stdout);
				}
				else {
					printf ("No Zero Crossing Neighbors\n\n");
				}

#endif


			}

		}
		VoxelStatus_muc[loc[0]] |= MAT_CHECKED_VOXEL; // Marking the Checked Voxel

		
	} // while ((int)InitialBoundaryLocs_map.size()>0)
}


template <class _DataType>
double cTFGeneration<_DataType>::getDataValue(double *Loc_d, double *NormalizedGradVec, double Position_d)
{
	int		k;
	double	NewLoc[3];
	
	for (k=0; k<3; k++) NewLoc[k] = Loc_d[k] + NormalizedGradVec[k]*Position_d;
	return DataInterpolation(NewLoc);
}


// Marking Boundary Voxels
template <class _DataType>
int cTFGeneration<_DataType>::ComputeSigma(double *CurrLoc, double *NormalizedGradVec, double& SigmaPosDir_Ret,
								double& SigmaNegDir_Ret, double& SecondDPosDir_Ret, double& SecondDNegDir_Ret)
{
	int		FoundPosDir, FoundNegDir;



	// Positive Direction
	FoundPosDir = ComputeSigma(CurrLoc, NormalizedGradVec, SigmaPosDir_Ret, SecondDPosDir_Ret, 0.1);

	// Negative Direction
	FoundNegDir = ComputeSigma(CurrLoc, NormalizedGradVec, SigmaNegDir_Ret, SecondDNegDir_Ret, -0.1);

	if (FoundPosDir && FoundNegDir) return 2;
	else if (FoundPosDir) return 1; // Found positive directin only
	else if (FoundNegDir) return -1;// Found negative directin only
	else return 0;
}


template <class _DataType>
int cTFGeneration<_DataType>::ComputeSigma(double *CurrLoc, double *NormalizedGradVec, 
											double& Sigma_Ret, double& SecondD_Ret, double Increase_d)
{
	int			k;
	double		LocAlongGradDir[3], Step_d, FirstD_d[3], SecondD_d[3];
	double		HalfStep_d, CurrStep_d;


	for (k=0; k<3; k++) LocAlongGradDir[k] = CurrLoc[k];
	FirstD_d[0] = GradientInterpolation(LocAlongGradDir);
	SecondD_d[0] = SecondDInterpolation(LocAlongGradDir);


#ifdef	DEBUG_TF_SIGMA
	int		Print_i;
	if (((int)CurrLoc[0])==80) Print_i = true;
	else Print_i = false;
	if (Print_i) {
		printf ("(1st, 2nd) = (%7.4f, %7.4f) ", FirstD_d[0], SecondD_d[0]);
		printf ("Loc=(%5.2f, %5.2f, %5.2f) ", LocAlongGradDir[0], LocAlongGradDir[1], LocAlongGradDir[2]);
		printf ("Position = %7.4f", 0.0);
		printf ("\n"); fflush (stdout);
	}
#endif

	for (k=0; k<3; k++) LocAlongGradDir[k] = CurrLoc[k] + NormalizedGradVec[k]*Increase_d;
	FirstD_d[1] = GradientInterpolation(LocAlongGradDir);
	SecondD_d[1] = SecondDInterpolation(LocAlongGradDir);

#ifdef	DEBUG_TF_SIGMA
	if (Print_i) {
		printf ("(1st, 2nd) = (%7.4f, %7.4f) ", FirstD_d[1], SecondD_d[1]);
		printf ("Loc=(%5.2f, %5.2f, %5.2f) ", LocAlongGradDir[0], LocAlongGradDir[1], LocAlongGradDir[2]);
		printf ("Position = %7.4f", Increase_d);
		printf ("\n"); fflush (stdout);
	}
#endif


	for (Step_d=Increase_d*2.0; fabs(Step_d)<=5.0+1e-5; Step_d+=Increase_d) {

		// When the location goes out of volume, then break the loop
		for (k=0; k<3; k++) {
			LocAlongGradDir[k] = CurrLoc[k] + NormalizedGradVec[k]*Step_d;
			if (LocAlongGradDir[k]<0.0) return false;
		}
		FirstD_d[2] = GradientInterpolation(LocAlongGradDir);
		SecondD_d[2] = SecondDInterpolation(LocAlongGradDir);

		if (fabs(SecondD_d[0]) < fabs(SecondD_d[1]) && 
			fabs(SecondD_d[2]) < fabs(SecondD_d[1])) {
			
			CurrStep_d = Step_d - Increase_d;
			do {
				
				Increase_d /= 2.0;
				HalfStep_d = CurrStep_d - Increase_d;
				for (k=0; k<3; k++) LocAlongGradDir[k] = CurrLoc[k] + NormalizedGradVec[k]*HalfStep_d;
				SecondD_d[0] = SecondDInterpolation(LocAlongGradDir);

				HalfStep_d = CurrStep_d + Increase_d;
				for (k=0; k<3; k++) LocAlongGradDir[k] = CurrLoc[k] + NormalizedGradVec[k]*HalfStep_d;
				SecondD_d[2] = SecondDInterpolation(LocAlongGradDir);

#ifdef	DEBUG_TF_SIGMA
	if (Print_i) {
		printf ("F''=(%8.5f, %8.5f, %8.5f) ", SecondD_d[0], SecondD_d[1], SecondD_d[2]);
		printf ("Loc=(%5.2f, %5.2f, %5.2f) ", LocAlongGradDir[0], LocAlongGradDir[1], LocAlongGradDir[2]);
		printf ("Position = %9.6f ", CurrStep_d);
		printf ("Biggest = ");
		if (fabs(SecondD_d[0])>fabs(SecondD_d[1])) printf (" 0 ");
		else if (fabs(SecondD_d[2])>fabs(SecondD_d[1])) printf (" 2 ");
		else printf (" 1 ");
		printf ("\n"); fflush (stdout);
	}
#endif

				if (fabs(SecondD_d[0])>fabs(SecondD_d[1])) {
					SecondD_d[1] = SecondD_d[0];
					CurrStep_d -= Increase_d;
				}
				else if (fabs(SecondD_d[2])>fabs(SecondD_d[1])) {
					SecondD_d[1] = SecondD_d[2];
					CurrStep_d += Increase_d;
				}
				
			} while (fabs(Increase_d)>1e-5);
			

			for (k=0; k<3; k++) LocAlongGradDir[k] = CurrLoc[k] + NormalizedGradVec[k]*CurrStep_d;
			FirstD_d[1] = GradientInterpolation(LocAlongGradDir);
			Sigma_Ret = FirstD_d[1]/(fabs(SecondD_d[1])*sqrt(exp((double)1.0)));
			SecondD_Ret = SecondD_d[1];

#ifdef	DEBUG_TF_SIGMA
		if (Print_i) {
			printf ("Sigma = %f ", Sigma_Ret);
			printf ("(1stD, 2ndD) = (%9.6f, %9.6f) ", FirstD_d[1], SecondD_d[1]);
			printf ("Position = %9.6f ", CurrStep_d);
			printf ("\n"); fflush (stdout);
		}
#endif

			return true;
		}
		for (k=0; k<=1; k++) FirstD_d[k] = FirstD_d[k+1];
		for (k=0; k<=1; k++) SecondD_d[k] = SecondD_d[k+1];
	}
	return false;
}





// Marking Boundary Voxels
template <class _DataType>
void cTFGeneration<_DataType>::ComputeAlpha(_DataType MatMin, _DataType MatMax)
{
	int			i, k, loc[3], DataCoor[3], NextDataCoor[3];
	double		GradVec[3], LocAlongGradDir[3], Tempd;
	double		Step, FirstD_d[3], SecondD_d[3], ZeroCrossingLoc_d[3];
	double		LocalMaxGradient, LocalExtremeSecondD, LocalExtremeSecondDPos;
	double		Sigma_d, PositionAlongGrad_d, Alpha_d, Increase;
	double		DataPosFromZeroCrossingLoc_d, DataValue_d;
	int			NumNeighbors, Neighbor26_i[26], FoundZeroCrossingLoc;
	int			OutOfVolume, FoundLocalExtremeSecondD, LocalMinDirection_i;


	map<int, unsigned char> MaterialBoundary_m;
	MaterialBoundary_m.clear();

	int		NumMatBoundary=0;
	for (i=0; i<WHD_mi; i++) {

//		if (MatMin<=Data_mT[i] && Data_mT[i]<=MatMax && VoxelStatus_muc[i]==0) {
//			VoxelStatus_muc[i] = MAT_INSIDE_BOUNDARY; // Mark the Material Volume
//		}
		if (IsMaterialBoundaryUsingMinMax(i, MatMin, MatMax)) {
			MaterialBoundary_m[i]=(unsigned char)0; // Add it to the map
			NumMatBoundary++;
		}
	}
	cout << endl << "Num Mat Boundary = " << NumMatBoundary << endl;
	

	map<int, unsigned char>::iterator	Boundary_it;
	
#ifdef	DEBUG_TF
	int		NumRepeat = 0;
	map<int, unsigned char> Checked_map;
	map<int, unsigned char>::iterator	Checked_it;
	Checked_map.clear();
#endif
	
	while ((int)MaterialBoundary_m.size()>0) {


#ifdef	DEBUG_TF
		NumRepeat++;
		if (NumRepeat%1000==0) {
			printf ("Num Repeat = %d, ", NumRepeat);
			printf ("Size of Map = %d\n", (int)MaterialBoundary_m.size());
		}
//		if (NumRepeat++ > 1000000) exit(1);
#endif		
	
		Boundary_it = MaterialBoundary_m.begin();
		loc[0] = (*Boundary_it).first;
		DataCoor[2] = loc[0]/WtimesH_mi;
		DataCoor[1] = (loc[0] - DataCoor[2]*WtimesH_mi)/Width_mi;
		DataCoor[0] = loc[0] % Width_mi;
		MaterialBoundary_m.erase(loc[0]);


#ifdef	DEBUG_TF
		unsigned char	Count;
		Checked_it = Checked_map.find(loc[0]);
		if (Checked_it == Checked_map.end()) Checked_map[loc[0]] = (unsigned char)0;
		else {
			Count = (*Checked_it).second + 1;
			Checked_map[loc[0]] = Count;
			if (Count>1) {
				printf ("DataCoor = %3d %3d %3d is re-checked. ", DataCoor[0], DataCoor[1], DataCoor[2]);
				printf ("# = %d, Marking = %d\n", Count, (int)VoxelStatus_muc[loc[0]]);
			}
		}
#endif

		// If the location is already checked, then skip it.
		// This is the same as if (VoxelStatus_muc[loc[0]]>0) continue;
		if (VoxelStatus_muc[loc[0]]==MAT_BOUNDARY_HAS_ALPHA ||
			VoxelStatus_muc[loc[0]]==MAT_OTHER_MATERIALS) continue;

		// Getting the gradient vector at the position
		for (k=0; k<3; k++) GradVec[k] = (double)GradientVec_mf[loc[0]*3 + k];
		Tempd = sqrt (GradVec[0]*GradVec[0] + GradVec[1]*GradVec[1] + GradVec[2]*GradVec[2]);
		if (fabs(Tempd)<1e-5) {
			if (MatMin<=Data_mT[i] && Data_mT[i]<=MatMax) {
				VoxelStatus_muc[loc[0]] = MAT_BOUNDARY_HAS_ALPHA;
				AlphaVolume_mf[loc[0]] = (float)1.0;
			}
			else {
				VoxelStatus_muc[loc[0]] = MAT_OTHER_MATERIALS;
				AlphaVolume_mf[loc[0]] = (float)0.0;
			}
			continue;
		}
		for (k=0; k<3; k++) GradVec[k] /= Tempd; // Normalize the gradient vector


		//-------------------------------------------------------------------------------------------
		// Finding the local maxima of gradient magnitudes along the gradient direction. 
		// It climbs the mountain of gradient magnitudes to find the zero-crossing second derivative 
		// Return ZeroCrossingLoc_d, LocalMaxGradient, DataPosFromZeroCrossingLoc_d
		//-------------------------------------------------------------------------------------------
		double		CurrDataLoc_d[3];
		for (k=0; k<3; k++) CurrDataLoc_d[k] = (double)DataCoor[k];
		FoundZeroCrossingLoc = FindZeroCrossingLocation(CurrDataLoc_d, GradVec, ZeroCrossingLoc_d, 
														LocalMaxGradient, DataPosFromZeroCrossingLoc_d);

		if (!FoundZeroCrossingLoc) {
			VoxelStatus_muc[loc[0]] = MAT_OTHER_MATERIALS;
			continue;
		}

		// Inverse the direction
		// From the Zero-Crossing Location to the given data position
		if (DataPosFromZeroCrossingLoc_d<0.0) Increase = 0.1;
		else Increase = -0.1;
		LocalExtremeSecondD = 0.0;

		//--------------------------------------------------------------------------------------
		// The Positive or Negative Direction of a Gradient Vector
		//--------------------------------------------------------------------------------------
		for (k=0; k<3; k++) LocAlongGradDir[k] = ZeroCrossingLoc_d[k];
		FirstD_d[0] = GradientInterpolation(LocAlongGradDir);
		SecondD_d[0] = SecondDInterpolation(LocAlongGradDir);
		for (k=0; k<3; k++) LocAlongGradDir[k] = ZeroCrossingLoc_d[k] + GradVec[k]*Increase;
		FirstD_d[1] = GradientInterpolation(LocAlongGradDir);
		SecondD_d[1] = SecondDInterpolation(LocAlongGradDir);

		if (SecondD_d[1]<0.0) LocalMinDirection_i = 1;
		else LocalMinDirection_i = 0;

		FoundLocalExtremeSecondD = 0;
		for (Step=Increase*2.0; fabs(Step)<=5.0+1e-5; Step+=Increase) {
			// When the location goes out of volume, then break the loop
			OutOfVolume = 0;
			for (k=0; k<3; k++) {
				LocAlongGradDir[k] = ZeroCrossingLoc_d[k] + GradVec[k]*Step;
				if (LocAlongGradDir[k]<0.0) OutOfVolume = 1;
			}
			if (OutOfVolume) break;

			FirstD_d[2] = GradientInterpolation(LocAlongGradDir);
			SecondD_d[2] = SecondDInterpolation(LocAlongGradDir);

			FoundLocalExtremeSecondD = 0;
			if (LocalMinDirection_i) {
				if (SecondD_d[0] > SecondD_d[1] && 
//					SecondD_d[2] > SecondD_d[1] && SecondD_d[1]<0.0) {
					SecondD_d[2] > SecondD_d[1]) {
					LocalExtremeSecondD = SecondD_d[1];
					LocalExtremeSecondDPos = fabs(Step);
					FoundLocalExtremeSecondD = 1;
				}
			}
			else { // Local Max Direction
				if (SecondD_d[0] < SecondD_d[1] &&
//					SecondD_d[2] < SecondD_d[1] && SecondD_d[1]>0.0) {
					SecondD_d[2] < SecondD_d[1]) {
					LocalExtremeSecondD = SecondD_d[1];
					LocalExtremeSecondDPos = fabs(Step);
					FoundLocalExtremeSecondD = 1;
				}
			}

			if (FoundLocalExtremeSecondD) break;
			for (k=0; k<=1; k++) FirstD_d[k] = FirstD_d[k+1];
			for (k=0; k<=1; k++) SecondD_d[k] = SecondD_d[k+1];
		}

		// Assertions
		if (fabs(Step)>=5.0) {
			printf ("Too much steps \n");
//			exit(1);
		}

		if (FoundLocalExtremeSecondD) {

			Sigma_d = LocalMaxGradient/LocalExtremeSecondD/exp((double)1.0);

			// Adding Neighbor Voxels
			if (DataPosFromZeroCrossingLoc_d<0.0) Increase = 1.0001;
			else Increase = -1.0001;
			for (Step=0.0; fabs(Step)<=fabs(Sigma_d); Step+=Increase) {
				// Adding neighbors
				for (k=0; k<3; k++) LocAlongGradDir[k] = ZeroCrossingLoc_d[k] + GradVec[k]*Step;
				if (LocAlongGradDir[0]<0.0 || LocAlongGradDir[1]<0.0 || LocAlongGradDir[2]<0.0) break;

				for (k=0; k<3; k++) NextDataCoor[k] = (int)floor(LocAlongGradDir[k]);
				loc[1] = NextDataCoor[2]*WtimesH_mi + NextDataCoor[1]*Width_mi + NextDataCoor[0];
				
//-------------------------------
// Make sure that it will be removed 
				MaterialBoundary_m[loc[1]]=(unsigned char)0; // Adding it to the map
				NumNeighbors = FindNeighbors8Or26(loc[1], Neighbor26_i);
				for (k=0; k<NumNeighbors; k++) {
					if (VoxelStatus_muc[Neighbor26_i[k]] != MAT_BOUNDARY_HAS_ALPHA &&
						VoxelStatus_muc[Neighbor26_i[k]] != MAT_OTHER_MATERIALS &&
						GradientMag_mf[Neighbor26_i[k]] > FirstD_d[1] ) {
						MaterialBoundary_m[Neighbor26_i[k]]=(unsigned char)0; // Adding it to the map
					}
				}
//-------------------------------

			}
			// Increase = -1.0001, or 1.0001
			for (k=0; k<3; k++) LocAlongGradDir[k] = ZeroCrossingLoc_d[k] + GradVec[k]*(fabs(Sigma_d)*2.0*Increase);
			DataValue_d = DataInterpolation(LocAlongGradDir);
			if (MatMin <= DataValue_d && DataValue_d <= MatMax) {
				AlphaVolume_mf[loc[0]] = (float)1.0;
				VoxelStatus_muc[loc[0]] = MAT_BOUNDARY_HAS_ALPHA;
			}
			else {
				PositionAlongGrad_d = -(Sigma_d*Sigma_d) * SecondDerivative_mf[loc[0]] / GradientMag_mf[loc[0]];
				Alpha_d = (1.0-fabs(PositionAlongGrad_d));
				if (Alpha_d<0.0) Alpha_d = 0.0;
				AlphaVolume_mf[loc[0]] = (float)Alpha_d;
				VoxelStatus_muc[loc[0]] = MAT_BOUNDARY_HAS_ALPHA;
			}
			

			MaterialBoundary_m.erase(loc[0]);

		}
		else {
			VoxelStatus_muc[loc[0]] = MAT_OTHER_MATERIALS;
			MaterialBoundary_m.erase(loc[0]);
		}
	
	}
	
		
		
}


// Marking Boundary Voxels
template <class _DataType>
void cTFGeneration<_DataType>::ComputeSigmaAtBoundary(_DataType MatMin, _DataType MatMax)
{
	int			i, k, loc[3];
	int			DataCoor[3];
	double		GradVec[3], LocAlongGradDir[3], Step, Tempd;
	double		ConsecutiveGradM_d[3], ConsecutiveSecondD_d[3];
	double		DataValue_d, FirstD_d, SecondD_d;
	map<int, float> 	GradPositiveDirBoundary_map;
	map<int, float> 	GradNegativeDirBoundary_map;
	map<int, float>::iterator currXYZ_it;
	int			CurrPointi, FoundLocalMinSecondD, FoundLocalMaxSecondD, PositiveDirHasSecondMax;
	int			PositiveDirIsInsideBoundary=0, NegativeDirIsInsideBoundary=0, NextDataCoor[3];
	double		LocalMaxGradM=-9999999, LocalMinSecondM=9999999, LocalMaxSecondM=-9999999, Sigma;
	double		SamplingDistance = 0.001, ZeroCrossingLoc_d[3];


	MaterialBoundary_mm.clear(); // Remove All Saved Boundaries

	for (i=0; i<WHD_mi; i++) {

		// if the voxel has the zero-crossing second derivative value, then ...
		if (VoxelStatus_muc[i]==MAT_ZERO_CROSSING) {

			for (k=0; k<3; k++) GradVec[k] = (double)GradientVec_mf[i*3 + k];
			Tempd = sqrt (GradVec[0]*GradVec[0] + GradVec[1]*GradVec[1] + GradVec[2]*GradVec[2]);
			if (fabs(Tempd)<1e-5) continue;
			for (k=0; k<3; k++) GradVec[k] /= Tempd; // Normalize the gradient vector

			DataCoor[2] = i / WtimesH_mi;
			DataCoor[1] = (i - DataCoor[2]*WtimesH_mi) / Height_mi;
			DataCoor[0] = i % Width_mi;

			printf ("Z Coor = %d\n", DataCoor[2]);
			
#ifdef 	DEBUG_TF
			int		PrintSliceNum = 25;
			if (DataCoor[2]==PrintSliceNum) {
				printf ("\nOriginal Position (LocationXYZ), data, 1stD, 2ndD\n");
				printf ("Loc = (%3d %3d %3d) ", DataCoor[0], DataCoor[1], DataCoor[2]);
				printf ("%3d %7.2f %7.2f ", (int)Data_mT[i], GradientMag_mf[i], SecondDerivative_mf[i]);
				printf ("GVec = (%5.2f %5.2f %5.2f)\n", GradVec[0], GradVec[1], GradVec[2]);
			}
#endif


			PositiveDirIsInsideBoundary=0;
			NegativeDirIsInsideBoundary=0;
			LocalMaxGradM=-9999999;
			LocalMinSecondM=9999999;
			LocalMaxSecondM=-9999999;


			CurrPointi=0; // Position --> -1.0
			for (k=0; k<3; k++) LocAlongGradDir[k] = (double)DataCoor[k] + GradVec[k]*(-1.0);
			FirstD_d = GradientInterpolation(LocAlongGradDir[0], LocAlongGradDir[1], LocAlongGradDir[2]);
			SecondD_d= SecondDInterpolation(LocAlongGradDir[0], LocAlongGradDir[1], LocAlongGradDir[2]);
			ConsecutiveGradM_d[CurrPointi] = FirstD_d;
			ConsecutiveSecondD_d[CurrPointi] = SecondD_d;

			CurrPointi=1; // Position --> -1.0 + SamplingDistance
			for (k=0; k<3; k++) LocAlongGradDir[k] = (double)DataCoor[k] + GradVec[k]*(-1.0 + SamplingDistance);
			FirstD_d = GradientInterpolation(LocAlongGradDir[0], LocAlongGradDir[1], LocAlongGradDir[2]);
			SecondD_d= SecondDInterpolation(LocAlongGradDir[0], LocAlongGradDir[1], LocAlongGradDir[2]);
			ConsecutiveGradM_d[CurrPointi] = FirstD_d;
			ConsecutiveSecondD_d[CurrPointi] = SecondD_d;

			for (Step=-0.8; Step<=1.0; Step+=SamplingDistance) {
			
				for (k=0; k<3; k++) LocAlongGradDir[k] = (double)DataCoor[k] + GradVec[k]*Step;
				FirstD_d = GradientInterpolation(LocAlongGradDir[0], LocAlongGradDir[1], LocAlongGradDir[2]);
				SecondD_d= SecondDInterpolation(LocAlongGradDir[0], LocAlongGradDir[1], LocAlongGradDir[2]);
				ConsecutiveGradM_d[2] = FirstD_d;
				ConsecutiveSecondD_d[2] = SecondD_d;

				// Check local maxima of gradient magnitudes and Zero-crossing Second derivatives
				if (ConsecutiveGradM_d[0] < ConsecutiveGradM_d[1] && 
					ConsecutiveGradM_d[2] < ConsecutiveGradM_d[1] && 
					ConsecutiveSecondD_d[0]*ConsecutiveSecondD_d[2] < 0.0) {
					LocalMaxGradM = ConsecutiveGradM_d[1];
					break;
				}
				for (k=0; k<=1; k++) ConsecutiveGradM_d[k] = ConsecutiveGradM_d[k+1];
				for (k=0; k<=1; k++) ConsecutiveSecondD_d[k] = ConsecutiveSecondD_d[k+1];
			}
			for (k=0; k<3; k++) ZeroCrossingLoc_d[k] = LocAlongGradDir[k];


			//---------------------------------------------------------------------------------------
			// Positive Direction
			//---------------------------------------------------------------------------------------
			double		Increase = 0.01;

			FoundLocalMinSecondD = 0;
			FoundLocalMaxSecondD = 0;
			PositiveDirHasSecondMax = -1;
			GradPositiveDirBoundary_map.clear();
			
			CurrPointi=0; // Current Position - 0.2
			for (k=0; k<3; k++) LocAlongGradDir[k] = (double)DataCoor[k] - GradVec[k]*Increase;
			FirstD_d = GradientInterpolation(LocAlongGradDir[0], LocAlongGradDir[1], LocAlongGradDir[2]);
			SecondD_d= SecondDInterpolation(LocAlongGradDir[0], LocAlongGradDir[1], LocAlongGradDir[2]);
			ConsecutiveGradM_d[CurrPointi] = FirstD_d;
			ConsecutiveSecondD_d[CurrPointi] = SecondD_d;

			CurrPointi=1; // Current Position
			for (k=0; k<3; k++) LocAlongGradDir[k] = (double)DataCoor[k];
			FirstD_d = GradientInterpolation(LocAlongGradDir[0], LocAlongGradDir[1], LocAlongGradDir[2]);
			SecondD_d= SecondDInterpolation(LocAlongGradDir[0], LocAlongGradDir[1], LocAlongGradDir[2]);
			ConsecutiveGradM_d[CurrPointi] = FirstD_d;
			ConsecutiveSecondD_d[CurrPointi] = SecondD_d;

			for (Step=Increase; Step<=5.0; Step+=Increase) {
			
				for (k=0; k<3; k++) LocAlongGradDir[k] = (double)DataCoor[k] + GradVec[k]*Step;
				FirstD_d = GradientInterpolation(LocAlongGradDir[0], LocAlongGradDir[1], LocAlongGradDir[2]);
				SecondD_d= SecondDInterpolation(LocAlongGradDir[0], LocAlongGradDir[1], LocAlongGradDir[2]);
				ConsecutiveGradM_d[2] = FirstD_d;
				ConsecutiveSecondD_d[2] = SecondD_d;

				// Insert the currect location to the map, GradPositiveDirBoundary_map
				// If the same location is inserted, then it will be overwritten
				for (k=0; k<3; k++) NextDataCoor[k] = (int)floor(LocAlongGradDir[k]+0.5);
				loc[1] = NextDataCoor[2]*WtimesH_mi + NextDataCoor[1]*Width_mi + NextDataCoor[0];
				GradPositiveDirBoundary_map[loc[1]] = (float)0.0; // Add the boundary to the map

				if (ConsecutiveSecondD_d[0] > ConsecutiveSecondD_d[1] && 
					ConsecutiveSecondD_d[2] > ConsecutiveSecondD_d[1] && !FoundLocalMinSecondD) {
//					if (ConsecutiveGradM_d[0] > ConsecutiveGradM_d[1] && 
//						ConsecutiveGradM_d[1] > ConsecutiveGradM_d[2]) {
						LocalMinSecondM = ConsecutiveSecondD_d[1];
						FoundLocalMinSecondD = 1;
						PositiveDirHasSecondMax = 0;
						break;
//					}
				}
				if (ConsecutiveSecondD_d[0] < ConsecutiveSecondD_d[1] && 
					ConsecutiveSecondD_d[2] < ConsecutiveSecondD_d[1] && !FoundLocalMaxSecondD) {
//					if (ConsecutiveGradM_d[0] > ConsecutiveGradM_d[1] && 
//						ConsecutiveGradM_d[1] > ConsecutiveGradM_d[2]) {
						LocalMaxSecondM = ConsecutiveSecondD_d[1];
						FoundLocalMaxSecondD = 1;
						PositiveDirHasSecondMax = 1;
						break;
//					}
				}
				for (k=0; k<=1; k++) ConsecutiveGradM_d[k] = ConsecutiveGradM_d[k+1];
				for (k=0; k<=1; k++) ConsecutiveSecondD_d[k] = ConsecutiveSecondD_d[k+1];

			}

			for (k=0; k<3; k++) NextDataCoor[k] = (int)ceil(LocAlongGradDir[k]+0.5);
			DataValue_d = DataInterpolation(NextDataCoor[0], NextDataCoor[1], NextDataCoor[2]);
			if (MatMin <= DataValue_d && DataValue_d <= MatMax) PositiveDirIsInsideBoundary = 1;


			// Assertion
			if (FoundLocalMinSecondD && FoundLocalMaxSecondD) { 
				printf ("Error in searching extremum of second derivatives\n");
				exit(1);
			}
			if (Step>=5.0) {
				printf ("Error in finding local extrem in positive direction\n");
				exit(1);
			}
			if (PositiveDirHasSecondMax==-1) {
				printf ("Error in finding second derivative maximum\n");
				exit(1);
			}

			printf ("FoundLocal Min & Max SecondD = %d, %d\n", FoundLocalMinSecondD, FoundLocalMaxSecondD);
			printf ("Local Min & Max SecondD = %7.2f, %7.2f\n", LocalMinSecondM, LocalMaxSecondM);

			//---------------------------------------------------------------------------------------
			// Negative Direction
			//---------------------------------------------------------------------------------------
			double		Decrease = -0.01;

			GradNegativeDirBoundary_map.clear();
			CurrPointi=0;
			for (k=0; k<3; k++) LocAlongGradDir[k] = (double)DataCoor[k] - GradVec[k]*Decrease;
			FirstD_d = GradientInterpolation(LocAlongGradDir[0], LocAlongGradDir[1], LocAlongGradDir[2]);
			SecondD_d= SecondDInterpolation(LocAlongGradDir[0], LocAlongGradDir[1], LocAlongGradDir[2]);
			ConsecutiveGradM_d[CurrPointi] = FirstD_d;
			ConsecutiveSecondD_d[CurrPointi] = SecondD_d;

			CurrPointi=1;
			for (k=0; k<3; k++) LocAlongGradDir[k] = (double)DataCoor[k];
			FirstD_d = GradientInterpolation(LocAlongGradDir[0], LocAlongGradDir[1], LocAlongGradDir[2]);
			SecondD_d= SecondDInterpolation(LocAlongGradDir[0], LocAlongGradDir[1], LocAlongGradDir[2]);
			ConsecutiveGradM_d[CurrPointi] = FirstD_d;
			ConsecutiveSecondD_d[CurrPointi] = SecondD_d;

			for (Step=Decrease; Step>=-5.0;Step+=Decrease) {

				for (k=0; k<3; k++) LocAlongGradDir[k] = (double)DataCoor[k] + GradVec[k]*Step;
				FirstD_d = GradientInterpolation(LocAlongGradDir[0], LocAlongGradDir[1], LocAlongGradDir[2]);
				SecondD_d= SecondDInterpolation(LocAlongGradDir[0], LocAlongGradDir[1], LocAlongGradDir[2]);
				ConsecutiveGradM_d[2] = FirstD_d;
				ConsecutiveSecondD_d[2] = SecondD_d;

				// Insert the currect location to the map, GradNegativeDirBoundary_map
				// If the same location is inserted, then it will be overwritten
				for (k=0; k<3; k++) NextDataCoor[k] = (int)ceil(LocAlongGradDir[k]+0.5);
				loc[1] = NextDataCoor[2]*WtimesH_mi + NextDataCoor[1]*Width_mi + NextDataCoor[0];
				GradNegativeDirBoundary_map[loc[1]] = (float)0.0; // Add the boundary to the map

				if (ConsecutiveSecondD_d[0] > ConsecutiveSecondD_d[1] && 
					ConsecutiveSecondD_d[2] > ConsecutiveSecondD_d[1] && !FoundLocalMinSecondD) {
//					if (ConsecutiveGradM_d[0] > ConsecutiveGradM_d[1] && 
//						ConsecutiveGradM_d[1] > ConsecutiveGradM_d[2]) {
						LocalMinSecondM = ConsecutiveSecondD_d[1];
						FoundLocalMinSecondD = 1;
						break;
//					}
				}
				if (ConsecutiveSecondD_d[0] < ConsecutiveSecondD_d[1] && 
					ConsecutiveSecondD_d[2] < ConsecutiveSecondD_d[1] && !FoundLocalMaxSecondD) {
//					if (ConsecutiveGradM_d[0] > ConsecutiveGradM_d[1] && 
//						ConsecutiveGradM_d[1] > ConsecutiveGradM_d[2]) {
						LocalMaxSecondM = ConsecutiveSecondD_d[1];
						FoundLocalMaxSecondD = 1;
						break;
//					}
				}
				for (k=0; k<=1; k++) ConsecutiveGradM_d[k] = ConsecutiveGradM_d[k+1];
				for (k=0; k<=1; k++) ConsecutiveSecondD_d[k] = ConsecutiveSecondD_d[k+1];
			}

			for (k=0; k<3; k++) NextDataCoor[k] = (int)ceil(LocAlongGradDir[k]+0.5);
			DataValue_d = DataInterpolation(NextDataCoor[0], NextDataCoor[1], NextDataCoor[2]);
			if (MatMin <= DataValue_d && DataValue_d <= MatMax) NegativeDirIsInsideBoundary = 1;


			// Assertion
			if (Step<=-5.0) {
				printf ("Error in finding local extrem in negative direction\n");

				printf ("FoundLocal Min & Max SecondD = %d, %d\n", FoundLocalMinSecondD, FoundLocalMaxSecondD);
				printf ("Local Min & Max SecondD = %7.2f, %7.2f\n", LocalMinSecondM, LocalMaxSecondM);
				
				for (Step=-1.5; Step<=1.5; Step+=0.01) {

					for (k=0; k<3; k++) LocAlongGradDir[k] = (double)DataCoor[k] + GradVec[k]*Step;
					DataValue_d = DataInterpolation(LocAlongGradDir);
					FirstD_d = GradientInterpolation(LocAlongGradDir);
					SecondD_d= SecondDInterpolation(LocAlongGradDir);
					printf ("Loc = (%5.2f %5.2f %5.2f), Data = %6.2f, 1st & 2nd = (%7.2f %7.2f) ", 
							LocAlongGradDir[0], LocAlongGradDir[1], LocAlongGradDir[2], DataValue_d, FirstD_d, SecondD_d);
					if (fabs(Step)<1e-6) printf ("--> Current Point\n");
					else printf ("\n");
				}


				exit(1);
			}
			if (FoundLocalMinSecondD && FoundLocalMaxSecondD) { }
			else {
				printf ("Error in finding local min and max second derivatives\n");
				exit(1);
			}
//			if (PositiveDirIsInsideBoundary && NegativeDirIsInsideBoundary) {
//				printf ("Error in finding inside boundary\n");
//				exit(1);
//			}
			

			// Isolated Voxels are removed
			if (!PositiveDirIsInsideBoundary && !NegativeDirIsInsideBoundary) {
				currXYZ_it  = GradPositiveDirBoundary_map.begin();
				for (k=0; k<(int)GradPositiveDirBoundary_map.size(); k++, currXYZ_it++) {
					loc[0] = (*currXYZ_it).first;
					VoxelStatus_muc[loc[0]]=(unsigned char)0; // Remove the voxel from the material volume
				}
				currXYZ_it  = GradNegativeDirBoundary_map.begin();
				for (k=0; k<(int)GradNegativeDirBoundary_map.size(); k++, currXYZ_it++) {
					loc[0] = (*currXYZ_it).first;
					VoxelStatus_muc[loc[0]]=(unsigned char)0; // Remove the voxel from the material volume
				}
				continue;
			}
			
			currXYZ_it  = GradPositiveDirBoundary_map.begin();
			if (PositiveDirIsInsideBoundary) {
				for (k=0; k<(int)GradPositiveDirBoundary_map.size(); k++, currXYZ_it++) {
					loc[0] = (*currXYZ_it).first;
					AddBoundary(loc[0], (float)0.0); // Location, Sigma
				}
			}
			else {
				if (PositiveDirHasSecondMax) {
					for (k=0; k<(int)GradPositiveDirBoundary_map.size(); k++, currXYZ_it++) {
						loc[0] = (*currXYZ_it).first;
						Sigma = LocalMaxGradM/LocalMaxSecondM/exp((double)1.0);
						Sigma *= Sigma;
						AddBoundary(loc[0], (float)(Sigma)); // Location, Sigma
					}
				}
				else {
					for (k=0; k<(int)GradPositiveDirBoundary_map.size(); k++, currXYZ_it++) {
						loc[0] = (*currXYZ_it).first;
						Sigma = LocalMaxGradM/LocalMinSecondM/exp((double)1.0);
						Sigma *= Sigma;
						AddBoundary(loc[0], (float)(Sigma)); // Location, Sigma
					}
				}
			}
			
			currXYZ_it  = GradNegativeDirBoundary_map.begin();
			if (NegativeDirIsInsideBoundary) {
				for (k=0; k<(int)GradNegativeDirBoundary_map.size(); k++, currXYZ_it++) {
					loc[0] = (*currXYZ_it).first;
					AddBoundary(loc[0], (float)0.0); // Location, Sigma
				}
			}
			else {
				if (PositiveDirHasSecondMax) {
					for (k=0; k<(int)GradNegativeDirBoundary_map.size(); k++, currXYZ_it++) {
						loc[0] = (*currXYZ_it).first;
						Sigma = LocalMaxGradM/LocalMinSecondM/exp((double)1.0);
						Sigma *= Sigma;
						AddBoundary(loc[0], (float)(Sigma)); // Location, Sigma
					}
				}
				else {
					for (k=0; k<(int)GradNegativeDirBoundary_map.size(); k++, currXYZ_it++) {
						loc[0] = (*currXYZ_it).first;
						Sigma = LocalMaxGradM/LocalMaxSecondM/exp((double)1.0);
						Sigma *= Sigma;
						AddBoundary(loc[0], (float)(Sigma)); // Location, Sigma
					}
				}
			}
		} // if (VoxelStatus_muc[i]==MAT_ZERO_CROSSING)
	} // for (i=0; i<WHD_mi; i++)
}


template<class _DataType>
void cTFGeneration<_DataType>::ComputeMinMax(_DataType MatMin, _DataType MatMax)
{
	int			i, loc[2];
	float		GradMin, GradMax, SecondMin, SecondMax;

	GradMin = 99999;
	GradMax = -99999;
	SecondMin = 99999;
	SecondMax = -99999;
	
	printf ("MatMin = %f, MatMax = %f\n", (float)MatMin, (float)MatMax);
	printf ("Size of MaterialBoundary_mm = %d\n", (int)MaterialBoundary_mm.size());
	
	map<int, float>::iterator Boundary_it  = MaterialBoundary_mm.begin();
	for (i=0; i<(int)MaterialBoundary_mm.size(); i++, Boundary_it++) {
		loc[0] = (*Boundary_it).first;
		if (GradMin > GradientMag_mf[loc[0]]) GradMin = GradientMag_mf[loc[0]];
		if (GradMax < GradientMag_mf[loc[0]]) GradMax = GradientMag_mf[loc[0]];
		if (SecondMin > SecondDerivative_mf[loc[0]]) SecondMin = SecondDerivative_mf[loc[0]];
		if (SecondMax < SecondDerivative_mf[loc[0]]) SecondMax = SecondDerivative_mf[loc[0]];
	}

	printf ("Grad Min & Max = (%f, %f)\n", GradMin, GradMax);
	printf ("Second Min & Max = (%f, %f)\n", SecondMin, SecondMax);
}


template<class _DataType>
void cTFGeneration<_DataType>::ComputeTF(char *OutFileName, _DataType MatMin, _DataType MatMax, int *TransFunc)
{

	// Step 1
	// Marking boundary voxels along gradient direction
	InitializeVoxelStatus();

	MarkingZeroCrossingLoc(MatMin, MatMax);
	
//	MarkingGradDirectedVoxels(MatMin, MatMax);
	SaveBoundaryVolume(OutFileName, MatMin, MatMax);


//	ComputeSigmaAtBoundary(MatMin, MatMax);
	
	
	
	// Step 2
	// Connectivity Check (Betti Number 1)
	
	// Step 3
	// Compute alpha at each boundary voxel
//	InitializeAlphaVolume();
//	ComputeAlpha(MatMin, MatMax);
//	SaveAlphaVolume(OutFileName, MatMin, MatMax);
	
	// Step 4
	// Converting the alpha values to the 2D transfer function

	printf ("TransFunc[0] = %d\n", TransFunc[0]);

}


// Three Value Location: X, Y, and Z
template<class _DataType>
void cTFGeneration<_DataType>::AddBoundary(int BoundaryLoc, float Sigma)
{
	map<int, float>::iterator Boundary_it  = MaterialBoundary_mm.find(BoundaryLoc);
	
	if (Boundary_it==MaterialBoundary_mm.end()) {
		// Add the current location to the map
		MaterialBoundary_mm[BoundaryLoc] = Sigma;
	}
}


template<class _DataType>
int cTFGeneration<_DataType>::getNumBoundary()
{
	return (int)MaterialBoundary_mm.size();
}


template <class _DataType>
int cTFGeneration<_DataType>::IsMaterialBoundaryUsingMinMax(int DataLoc, _DataType MatMin, _DataType MatMax)
{
	int			i, j, k, loc[3];
	int			XCoor, YCoor, ZCoor;
	
	
	// The given location should be between the min and max values
	if (Data_mT[DataLoc] < MatMin || MatMax < Data_mT[DataLoc]) return false;
	
	ZCoor = DataLoc / WtimesH_mi;
	YCoor = (DataLoc - ZCoor*WtimesH_mi) / Height_mi;
	XCoor = DataLoc % Width_mi;

	// Checking all 26 neighbors, whether at least one of them is a different material
	for (k=ZCoor-1; k<=ZCoor+1; k++) {
		if (k<0 || k>=Depth_mi) continue;
		for (j=YCoor-1; j<=YCoor+1; j++) {
			if (j<0 || j>=Height_mi) continue;
			for (i=XCoor-1; i<=XCoor+1; i++) {
				if (i<0 || i>=Width_mi) continue;

				loc[0] = k*WtimesH_mi + j*Width_mi + i;
				if (Data_mT[loc[0]] < MatMin || MatMax < Data_mT[loc[0]]) return true;

			}
		}
	}
	
	return false;
}

 // i = Data Location, j = Material Number
template <class _DataType>
int cTFGeneration<_DataType>::IsMaterialBoundary(int DataLoc, int MaterialNum)
{
	_DataType	DataValue;
	int			i, j, k, loc[3];
	int			XCoor, YCoor, ZCoor;
	
	
	ZCoor = DataLoc / WtimesH_mi;
	YCoor = (DataLoc - ZCoor*WtimesH_mi) / Height_mi;
	XCoor = DataLoc % Width_mi;

//	cout << "Data Location = " << DataLoc << " ";
//	cout << "X Y Z Coordinates = " << XCoor << ", " << YCoor << ", " << ZCoor << endl;

	// Checking all 26 neighbors, whether at least one of them is a different material
	for (k=ZCoor-1; k<=ZCoor+1; k++) {
		if (k<0 || k>=Depth_mi) return true;
		for (j=YCoor-1; j<=YCoor+1; j++) {
			if (j<0 || j>=Height_mi) return true;
			for (i=XCoor-1; i<=XCoor+1; i++) {
				if (i<0 || i>=Width_mi) return true;
				
				// DataValue = A Data Value in the Histogram
				loc[0] = k*WtimesH_mi + j*Width_mi + i;
				DataValue = (int)(((double)Data_mT[loc[0]]-MinData_mf)*HistogramFactorI_mf);

				loc[1] = (int)(DataValue*NumMaterial_mi + MaterialNum);
				if (MaterialProb_mf[loc[1]]<0.1) return true;
				
			}
		}
	}
	
	return false;
}


template <class _DataType>
int cTFGeneration<_DataType>::IsMaterialBoundary(int DataLoc, int MaterialNum, int Neighbor)
{
	_DataType	DataValue;
	int			i, j, k, loc[3];
	int			XCoor, YCoor, ZCoor;
	
	
	ZCoor = DataLoc / WtimesH_mi;
	YCoor = (DataLoc - ZCoor*WtimesH_mi) / Height_mi;
	XCoor = DataLoc % Width_mi;

//	cout << "Data Location = " << DataLoc << " ";
//	cout << "X Y Z Coordinates = " << XCoor << ", " << YCoor << ", " << ZCoor << endl;

	// Checking all 26 neighbors, whether at least one of them is a different material
	for (k=ZCoor-1; k<=ZCoor+1; k++) {
		if (k<0 || k>=Depth_mi) continue;
		for (j=YCoor-1; j<=YCoor+1; j++) {
			if (j<0 || j>=Height_mi) continue;
			for (i=XCoor-1; i<=XCoor+1; i++) {
				if (i<0 || i>=Width_mi) continue;
				
				// DataValue = A Data Value in the Histogram
				loc[0] = k*WtimesH_mi + j*Width_mi + i;
				DataValue = (int)(((double)Data_mT[loc[0]]-MinData_mf)*HistogramFactorI_mf);

				loc[1] = (int)(DataValue*NumMaterial_mi + MaterialNum);
				if (MaterialProb_mf[loc[1]]<0.1) {
					loc[2] = (int)(DataValue*NumMaterial_mi + Neighbor);
					if (MaterialProb_mf[loc[2]]>=0.1) return true;
				}
				
			}
		}
	}
	
	return false;
}


template<class _DataType>
void cTFGeneration<_DataType>::InitZeroCrossingVolume(int X, int Y, int Z)
{
	printf ("XYZ = %d, %d, %d\n", X, Y, Z);
#ifdef		SAVE_ZERO_CROSSING_VOLUME

	ZCVolumeGridSize_md = (double)1.0/10.0;
	int	GridSize_i = (int)(floor(1.0/ZCVolumeGridSize_md+0.5));
	VoxelVolume_muc = new unsigned char [GridSize_i*GridSize_i*GridSize_i];

	
	StartLoc_gi[0] = X;
	StartLoc_gi[1] = Y;
	StartLoc_gi[2] = Z;

	GridSize_i = (int)(floor(1.0/ZCVolumeGridSize_md+0.5));
	ZeroCrossingVoxels_muc = new unsigned char [64*64*64*GridSize_i*GridSize_i*GridSize_i];
	for (int i=0; i<64*64*64*GridSize_i*GridSize_i*GridSize_i; i++) {
		ZeroCrossingVoxels_muc[i] = (unsigned char)0;
	}

#endif
}

template<class _DataType>
void cTFGeneration<_DataType>::SaveZeroCrossingVolume(char *Prefix)
{
	int		binfile_fd2;
	char	FileName[200];


	sprintf (FileName, "%s_Voxel_%03d_%03d_%03d.raw", Prefix, StartLoc_gi[0], StartLoc_gi[1], StartLoc_gi[2]);
	cout << "Zero Crossing Volume File Name = " << FileName << endl;
	cout.flush();

	if ((binfile_fd2 = open (FileName, O_CREAT | O_WRONLY)) < 0) {
		printf ("could not open %s\n", FileName);
		exit(1);
	}
	if (write(binfile_fd2, ZeroCrossingVoxels_muc, sizeof(unsigned char)*640*640*640)
				!=(unsigned int)sizeof(unsigned char)*640*640*640) {
		printf ("The file could not be written : %s\n", FileName);
		close (binfile_fd2);
		exit(1);
	}
	close (binfile_fd2);

	if (chmod(FileName, S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH) < 0) {
		printf ("chmod was not worked to file %s\n", FileName);
		exit(1);
	}

}

template<class _DataType>
void cTFGeneration<_DataType>::InitBoundaryVolume()
{
	BV_GridSize_md = 1.0/10.0; // Boundary Volume Grid Size
	ZCVolumeGridSize_md = 1.0/10.0;
}


// Find 26 Neighbors
//
// Inputs:
//		CenterLoc: The Start Point to find zero-crossing (Z*WH + Y*W + X)
//		
// Returns:
//
// 		*Neighbors_Ret: Zero-crossing neighbors
//
// Notes:
// 		Grid Resolution = BV_GridSize_md
template<class _DataType>
int cTFGeneration<_DataType>::FindBoundaryConnectedNeighbors(int CenterLoc, int *Neighbors_Ret)
{
	int			i, j, k, l, loc[3], XYZCoor_i[3];
	int			FoundZeroCrossingLoc, XNeighbor[3], YNeighbor[3], ZNeighbor[3];
	double		XCoor_d, YCoor_d, ZCoor_d, GradVec_d[3], CurrLoc_d[3], Tempd;
	double		ZeroCrossingLoc_d[3], FirstDAtTheLoc_d, DataPosFromZeroCrossingLoc_d;
	double		Increase_d;


	XYZCoor_i[2] = CenterLoc / WtimesH_mi;
	XYZCoor_i[1] = (CenterLoc - XYZCoor_i[2]*WtimesH_mi) / Height_mi;
	XYZCoor_i[0] = CenterLoc % Width_mi;

	for (l=0; l<=2; l++) {
		XNeighbor[l] = 0;
		YNeighbor[l] = 0;
		ZNeighbor[l] = 0;
	}

	Increase_d = ZCVolumeGridSize_md;
	for (ZCoor_d=(double)XYZCoor_i[2]; ZCoor_d<(double)XYZCoor_i[2]+1.0; ZCoor_d+=Increase_d) {
		for (YCoor_d=(double)XYZCoor_i[1]; YCoor_d<(double)XYZCoor_i[1]+1.0; YCoor_d+=Increase_d) {
			for (XCoor_d=(double)XYZCoor_i[0]; XCoor_d<(double)XYZCoor_i[0]+1.0; XCoor_d+=Increase_d) {
			
				if (ZCVolumeGridSize_md-1e-5 <= fabs(XCoor_d-XYZCoor_i[0]) && 
					fabs(XCoor_d-XYZCoor_i[0]) <= 1.0-ZCVolumeGridSize_md-1e-5 &&
					ZCVolumeGridSize_md-1e-5 <= fabs(YCoor_d-XYZCoor_i[1]) && 
					fabs(YCoor_d-XYZCoor_i[1]) <= 1.0-ZCVolumeGridSize_md-1e-5 &&
					ZCVolumeGridSize_md-1e-5 <= fabs(ZCoor_d-XYZCoor_i[2]) && 
					fabs(ZCoor_d-XYZCoor_i[2]) <= 1.0-ZCVolumeGridSize_md-1e-5) continue;
			
				CurrLoc_d[0] = XCoor_d + ZCVolumeGridSize_md/2.0;
				CurrLoc_d[1] = YCoor_d + ZCVolumeGridSize_md/2.0;
				CurrLoc_d[2] = ZCoor_d + ZCVolumeGridSize_md/2.0;
				GradVecInterpolation(CurrLoc_d, GradVec_d);
				Tempd = sqrt(GradVec_d[2]*GradVec_d[2] + GradVec_d[1]*GradVec_d[1] + GradVec_d[0]*GradVec_d[0]);
				if (fabs(Tempd)<1e-5) continue;
				for (k=0; k<3; k++) GradVec_d[k] /= Tempd; // Normalize the gradient vector
				//Find Zero Connected Neighbors
				FoundZeroCrossingLoc = FindZeroCrossingLocation(CurrLoc_d, GradVec_d, ZeroCrossingLoc_d, 
											FirstDAtTheLoc_d, DataPosFromZeroCrossingLoc_d, 0.2);

				if (FoundZeroCrossingLoc && fabs(DataPosFromZeroCrossingLoc_d)<=ZCVolumeGridSize_md) {
					if (fabs(XCoor_d-XYZCoor_i[0])<1e-5) XNeighbor[0] = -1; // Negative X Direction
					if (fabs(XCoor_d-XYZCoor_i[0])>=(1.0-ZCVolumeGridSize_md-1e-5)) XNeighbor[1] = 1; // Positive X Direction
					
					if (fabs(YCoor_d-XYZCoor_i[1])<1e-5) YNeighbor[0] = -1; // Negative Y Direction
					if (fabs(YCoor_d-XYZCoor_i[1])>=(1.0-ZCVolumeGridSize_md-1e-5)) YNeighbor[1] = 1; // Positive Y Direction
					
					if (fabs(ZCoor_d-XYZCoor_i[2])<1e-5) ZNeighbor[0] = -1; // Negative Z Direction
					if (fabs(ZCoor_d-XYZCoor_i[2])>=(1.0-ZCVolumeGridSize_md-1e-5)) ZNeighbor[1] = 1; // Positive Z Direction

				}
				
			}
		}
	}

	int		XCoor_i, YCoor_i, ZCoor_i;
	map<int, unsigned char> Neighbor_map;
	Neighbor_map.clear();
	Neighbor_map[CenterLoc] = (unsigned char)0; // The Given Location
	
// Considering only 26 neighbors except the center location
	for (i=0; i<=2; i++) {
		for (j=0; j<=2; j++) {
			for (k=0; k<=2; k++) {
				XCoor_i = XYZCoor_i[0] + XNeighbor[k];
				YCoor_i = XYZCoor_i[1] + YNeighbor[j];
				ZCoor_i = XYZCoor_i[2] + ZNeighbor[i];
				if (XCoor_i<0 || XCoor_i>=Width_mi) continue;
				if (YCoor_i<0 || YCoor_i>=Height_mi) continue;
				if (ZCoor_i<0 || ZCoor_i>=Depth_mi) continue;

				loc[0] = ZCoor_i*WtimesH_mi + YCoor_i*Width_mi + XCoor_i;
				Neighbor_map[loc[0]] = (unsigned char)0;

			}
		}
	}

	Neighbor_map.erase(CenterLoc);
		
	map<int, unsigned char>::iterator Neighbors_it  = Neighbor_map.begin();
	
	for (i=0; i<(int)Neighbor_map.size(); Neighbors_it++, i++) {
		Neighbors_Ret[i] = 	(*Neighbors_it).first;
	}
	
	int NumNeighbors = (int)Neighbor_map.size();
	Neighbor_map.clear();
	return NumNeighbors;

}



// Find 8 Neighbors for 2D Case, 26 Neighbors for 3D Case
template<class _DataType>
int cTFGeneration<_DataType>::FindZeroConnectedNeighbors(int CenterLoc, int *Neighbors)
{

	int			i, j, k, l, loc[3], XYZCoor_i[3];
	int			FoundZeroCrossingLoc, XNeighbor[3], YNeighbor[3], ZNeighbor[3];
	double		XCoor_d, YCoor_d, ZCoor_d, GradVec_d[3], CurrLoc_d[3], Tempd;
	double		ZeroCrossingLoc_d[3], FirstDAtTheLoc_d, DataPosFromZeroCrossingLoc_d;
	double		Increase_d;


	XYZCoor_i[2] = CenterLoc / WtimesH_mi;
	XYZCoor_i[1] = (CenterLoc - XYZCoor_i[2]*WtimesH_mi) / Height_mi;
	XYZCoor_i[0] = CenterLoc % Width_mi;

	// Threshold Values
	// 1/8  --> 1/64
	// 1/32 --> 1/1024

	for (l=0; l<=2; l++) {
		XNeighbor[l] = 0;
		YNeighbor[l] = 0;
		ZNeighbor[l] = 0;
	}
	
	Increase_d = ZCVolumeGridSize_md;
	for (ZCoor_d=(double)XYZCoor_i[2]; ZCoor_d<(double)XYZCoor_i[2]+1.0; ZCoor_d+=Increase_d) {
		for (YCoor_d=(double)XYZCoor_i[1]; YCoor_d<(double)XYZCoor_i[1]+1.0; YCoor_d+=Increase_d) {
			for (XCoor_d=(double)XYZCoor_i[0]; XCoor_d<(double)XYZCoor_i[0]+1.0; XCoor_d+=Increase_d) {
			
				if (ZCVolumeGridSize_md-1e-5 <= fabs(XCoor_d-XYZCoor_i[0]) && 
					fabs(XCoor_d-XYZCoor_i[0]) <= 1.0-ZCVolumeGridSize_md-1e-5 &&
					ZCVolumeGridSize_md-1e-5 <= fabs(YCoor_d-XYZCoor_i[1]) && 
					fabs(YCoor_d-XYZCoor_i[1]) <= 1.0-ZCVolumeGridSize_md-1e-5 &&
					ZCVolumeGridSize_md-1e-5 <= fabs(ZCoor_d-XYZCoor_i[2]) && 
					fabs(ZCoor_d-XYZCoor_i[2]) <= 1.0-ZCVolumeGridSize_md-1e-5) continue;
			
				CurrLoc_d[0] = XCoor_d+ZCVolumeGridSize_md/2.0;
				CurrLoc_d[1] = YCoor_d+ZCVolumeGridSize_md/2.0;
				CurrLoc_d[2] = ZCoor_d+ZCVolumeGridSize_md/2.0;
				GradVecInterpolation(CurrLoc_d, GradVec_d);
				Tempd = sqrt(GradVec_d[2]*GradVec_d[2] + GradVec_d[1]*GradVec_d[1] + GradVec_d[0]*GradVec_d[0]);
				if (fabs(Tempd)<1e-5) continue;
				for (k=0; k<3; k++) GradVec_d[k] /= Tempd; // Normalize the gradient vector
				//Find Zero Connected Neighbors
				FoundZeroCrossingLoc = FindZeroCrossingLocation(CurrLoc_d, GradVec_d, ZeroCrossingLoc_d, 
											FirstDAtTheLoc_d, DataPosFromZeroCrossingLoc_d, 0.2);
				if (FoundZeroCrossingLoc && fabs(DataPosFromZeroCrossingLoc_d)<=ZCVolumeGridSize_md) {
					if (fabs(XCoor_d-XYZCoor_i[0])<1e-5) XNeighbor[0] = -1; // Negative Direction
					if (fabs(XCoor_d-XYZCoor_i[0])>=(1.0-ZCVolumeGridSize_md-1e-5)) XNeighbor[1] = 1; // Positive Direction
					
					if (fabs(YCoor_d-XYZCoor_i[1])<1e-5) YNeighbor[0] = -1;
					if (fabs(YCoor_d-XYZCoor_i[1])>=(1.0-ZCVolumeGridSize_md-1e-5)) YNeighbor[1] = 1;
					
					if (fabs(ZCoor_d-XYZCoor_i[2])<1e-5) ZNeighbor[0] = -1;
					if (fabs(ZCoor_d-XYZCoor_i[2])>=(1.0-ZCVolumeGridSize_md-1e-5)) ZNeighbor[1] = 1;

				}
			}
		}
	}

	int		XCoor_i, YCoor_i, ZCoor_i;
	map<int, unsigned char> Neighbor_map;
	Neighbor_map.clear();
	Neighbor_map[CenterLoc] = (unsigned char)0; // The Given Location
	
// Considering only 26 neighbors
	for (i=0; i<=2; i++) {
		for (j=0; j<=2; j++) {
			for (k=0; k<=2; k++) {
				XCoor_i = XYZCoor_i[0] + XNeighbor[k];
				YCoor_i = XYZCoor_i[1] + YNeighbor[j];
				ZCoor_i = XYZCoor_i[2] + ZNeighbor[i];
				if (XCoor_i<0 || XCoor_i>=Width_mi) continue;
				if (YCoor_i<0 || YCoor_i>=Height_mi) continue;
				if (ZCoor_i<0 || ZCoor_i>=Depth_mi) continue;
				

#ifdef		SAVE_ZERO_CROSSING_VOLUME

				if (StartLoc_gi[0]-20<=XCoor_i && StartLoc_gi[0]+64+20 > XCoor_i && 
					StartLoc_gi[1]-20<=YCoor_i && StartLoc_gi[1]+64+20 > YCoor_i && 
					StartLoc_gi[2]-20<=ZCoor_i && StartLoc_gi[2]+64+20 > ZCoor_i	) {
					loc[0] = ZCoor_i*WtimesH_mi + YCoor_i*Width_mi + XCoor_i;
					Neighbor_map[loc[0]] = (unsigned char)0;
				}
				else continue;

#else				
				loc[0] = ZCoor_i*WtimesH_mi + YCoor_i*Width_mi + XCoor_i;
				Neighbor_map[loc[0]] = (unsigned char)0;
#endif



			
			}
		}
	}


	Neighbor_map.erase(CenterLoc);
		
	map<int, unsigned char>::iterator Neighbors_it  = Neighbor_map.begin();
	
	for (i=0; i<(int)Neighbor_map.size(); Neighbors_it++, i++) {
		Neighbors[i] = 	(*Neighbors_it).first;
	}
	

#ifdef		SAVE_ZERO_CROSSING_VOLUME
	int 	GridSize_i;

	if (StartLoc_gi[0]<=XYZCoor_i[0] && StartLoc_gi[1]<=XYZCoor_i[1] && StartLoc_gi[2]<=XYZCoor_i[2] &&
		StartLoc_gi[0]+64 > XYZCoor_i[0] && StartLoc_gi[1]+64 > XYZCoor_i[1] && StartLoc_gi[2]+64 > XYZCoor_i[2]) {

		GridSize_i = (int)(floor(1.0/ZCVolumeGridSize_md+0.5));
		for (i=0; i<GridSize_i*GridSize_i*GridSize_i; i++) {
			VoxelVolume_muc[i] = (unsigned char)0;
		}


		cout << " TFG: Save Zero Crossing Loc at: XYZCoor_i = ";
		cout << "(" << XYZCoor_i[0] << " " << XYZCoor_i[1] << " " << XYZCoor_i[2] << ")  (";
		cout << (XYZCoor_i[0] - StartLoc_gi[0])*GridSize_i << " ";
		cout << (XYZCoor_i[1] - StartLoc_gi[1])*GridSize_i << " ";
		cout << (XYZCoor_i[2] - StartLoc_gi[2])*GridSize_i << ")" << endl;
		cout.flush();

	
		for (k=0,ZCoor_d=(double)XYZCoor_i[2];ZCoor_d<(double)XYZCoor_i[2]+1.0-1e-5;ZCoor_d+=ZCVolumeGridSize_md,k++){
			for (j=0,YCoor_d=(double)XYZCoor_i[1];YCoor_d<(double)XYZCoor_i[1]+1.0-1e-5;YCoor_d+=ZCVolumeGridSize_md,j++){
				for (i=0,XCoor_d=(double)XYZCoor_i[0];XCoor_d<(double)XYZCoor_i[0]+1.0-1e-5;XCoor_d+=ZCVolumeGridSize_md,i++){

					CurrLoc_d[0] = XCoor_d + ZCVolumeGridSize_md/2.0;
					CurrLoc_d[1] = YCoor_d + ZCVolumeGridSize_md/2.0;
					CurrLoc_d[2] = ZCoor_d + ZCVolumeGridSize_md/2.0;
					GradVecInterpolation(CurrLoc_d, GradVec_d);
					Tempd = sqrt(GradVec_d[2]*GradVec_d[2] + GradVec_d[1]*GradVec_d[1] + GradVec_d[0]*GradVec_d[0]);
					if (fabs(Tempd)<1e-5) continue;
					for (int l=0; l<3; l++) GradVec_d[l] /= Tempd; // Normalize the gradient vector
					
					//Find Zero Crossing Neighbors
					FoundZeroCrossingLoc = FindZeroCrossingLocation(CurrLoc_d, GradVec_d, ZeroCrossingLoc_d, 
												FirstDAtTheLoc_d, DataPosFromZeroCrossingLoc_d, 0.2);

					if (FoundZeroCrossingLoc) {
						loc[0] = k*GridSize_i*GridSize_i + j*GridSize_i + i;
						if  (fabs(DataPosFromZeroCrossingLoc_d)<=ZCVolumeGridSize_md*6) 
							VoxelVolume_muc[loc[0]] = (unsigned char)100;
						if  (fabs(DataPosFromZeroCrossingLoc_d)<=ZCVolumeGridSize_md*4) 
							VoxelVolume_muc[loc[0]] = (unsigned char)150;
						if  (fabs(DataPosFromZeroCrossingLoc_d)<=ZCVolumeGridSize_md*2) 
							VoxelVolume_muc[loc[0]] = (unsigned char)200;
						if  (fabs(DataPosFromZeroCrossingLoc_d)<=ZCVolumeGridSize_md) 
							VoxelVolume_muc[loc[0]] = (unsigned char)255;
					}
				}
			}
		}

		int		InsideVoxelLoc;
		
		for (k=0; k<GridSize_i; k++) {
			for (j=0; j<GridSize_i; j++) {
				for (i=0; i<GridSize_i; i++) {

					InsideVoxelLoc = k*GridSize_i*GridSize_i + j*GridSize_i + i;
					loc[0] = 	(XYZCoor_i[2] - StartLoc_gi[2])*64*64*GridSize_i*GridSize_i*GridSize_i +
									 k*64*64*GridSize_i*GridSize_i + 
								(XYZCoor_i[1] - StartLoc_gi[1])*64*GridSize_i*GridSize_i + j*64*GridSize_i +
								(XYZCoor_i[0] - StartLoc_gi[0])*GridSize_i + i;
					ZeroCrossingVoxels_muc[loc[0]] = VoxelVolume_muc[InsideVoxelLoc];
				}
			}
		}
	}
	else {
		GridSize_i = (int)(floor(1.0/ZCVolumeGridSize_md+0.5));
		cout << " The voxel is out of range in the zero-crossing volume: XYZCoor_i = (";
		cout << XYZCoor_i[0] << " " << XYZCoor_i[1] << " " << XYZCoor_i[2] << ")  (";
		cout << (XYZCoor_i[0] - StartLoc_gi[0])*GridSize_i << " ";
		cout << (XYZCoor_i[1] - StartLoc_gi[1])*GridSize_i << " ";
		cout << (XYZCoor_i[2] - StartLoc_gi[2])*GridSize_i << ")" << endl;
		cout.flush();
	}	

#endif



#ifdef		SAVE_ZERO_CROSSING_SECOND_DERIVATIVE
		cout << "Save Zero Crossing Second Derivative Voxel at ";
		cout << ": XYZCoor_i = " << XYZCoor_i[0] << " " << XYZCoor_i[1] << " " << XYZCoor_i[2] << endl;

		int 	j, binfile_fd2;
		char	VoxelFileName[500];
		int GridSize_i = (int)(floor(1.0/ZCVolumeGridSize_md+0.5));
		unsigned char	*VoxelVolume2_uc = new unsigned char [GridSize_i*GridSize_i*GridSize_i];
		for (i=0; i<GridSize_i*GridSize_i*GridSize_i; i++) {
			VoxelVolume2_uc[i] = (unsigned char)0;
		}

		for (k=0,ZCoor_d=(double)XYZCoor_i[2];ZCoor_d<(double)XYZCoor_i[2]+1.0-1e-5;ZCoor_d+=ZCVolumeGridSize_md,k++) {
			for (j=0,YCoor_d=(double)XYZCoor_i[1];YCoor_d<(double)XYZCoor_i[1]+1.0-1e-5;YCoor_d+=ZCVolumeGridSize_md,j++) {
				for (i=0,XCoor_d=(double)XYZCoor_i[0];XCoor_d<(double)XYZCoor_i[0]+1.0-1e-5;XCoor_d+=ZCVolumeGridSize_md,i++) {

					CurrLoc_d[0] = XCoor_d + ZCVolumeGridSize_md/2.0;
					CurrLoc_d[1] = YCoor_d + ZCVolumeGridSize_md/2.0;
					CurrLoc_d[2] = ZCoor_d + ZCVolumeGridSize_md/2.0;
					GradVecInterpolation(CurrLoc_d, GradVec_d);
					Tempd = sqrt(GradVec_d[2]*GradVec_d[2] + GradVec_d[1]*GradVec_d[1] + GradVec_d[0]*GradVec_d[0]);
					if (fabs(Tempd)<1e-5) continue;
					for (int l=0; l<3; l++) GradVec_d[l] /= Tempd; // Normalize the gradient vector
					
					//Find Zero Crossing Neighbors
					FoundZeroCrossingLoc = FindZeroCrossingLocation(CurrLoc_d, GradVec_d, ZeroCrossingLoc_d, 
												FirstDAtTheLoc_d, DataPosFromZeroCrossingLoc_d, 0.2);

					if (FoundZeroCrossingLoc) {
						loc[0] = k*GridSize_i*GridSize_i + j*GridSize_i + i;
						if  (fabs(DataPosFromZeroCrossingLoc_d)<=ZCVolumeGridSize_md*16) 
							VoxelVolume2_uc[loc[0]] = (unsigned char)50;
						if  (fabs(DataPosFromZeroCrossingLoc_d)<=ZCVolumeGridSize_md*8) 
							VoxelVolume2_uc[loc[0]] = (unsigned char)100;
						if  (fabs(DataPosFromZeroCrossingLoc_d)<=ZCVolumeGridSize_md*4) 
							VoxelVolume2_uc[loc[0]] = (unsigned char)150;
						if  (fabs(DataPosFromZeroCrossingLoc_d)<=ZCVolumeGridSize_md*2) 
							VoxelVolume2_uc[loc[0]] = (unsigned char)200;
						if  (fabs(DataPosFromZeroCrossingLoc_d)<=ZCVolumeGridSize_md) 
							VoxelVolume2_uc[loc[0]] = (unsigned char)255;
					}
				}
			}
		}

		sprintf (VoxelFileName, "Voxel_%03d_%03d_%03d.raw", XYZCoor_i[0], XYZCoor_i[1], XYZCoor_i[2]);
		cout << "Voxel File Name = " << VoxelFileName << endl;
		cout.flush();


		if ((binfile_fd2 = open (VoxelFileName, O_CREAT | O_WRONLY)) < 0) {
			printf ("could not open %s\n", VoxelFileName);
			exit(1);
		}
		if (write(binfile_fd2, VoxelVolume2_uc, sizeof(unsigned char)*GridSize_i*GridSize_i*GridSize_i)
					!=(unsigned int)sizeof(unsigned char)*GridSize_i*GridSize_i*GridSize_i) {
			printf ("The file could not be written : %s\n", VoxelFileName);
			close (binfile_fd2);
			exit(1);
		}
		close (binfile_fd2);
		if (chmod(VoxelFileName, S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH) < 0) {
			printf ("chmod was not worked to file %s\n", VoxelFileName);
			exit(1);
		}
		delete [] VoxelVolume2_uc;
#endif


	int NumNeighbors = (int)Neighbor_map.size();
	Neighbor_map.clear();
	return NumNeighbors;

}

// Find 8 Neighbors for 2D Case, 26 Neighbors for 3D Case
template<class _DataType>
int cTFGeneration<_DataType>::FindZeroConnectedNeighbors2(int CenterLoc, int *Neighbors)
{

	int			i, k, l, loc[3], XYZCoor_i[3];
	int			FoundZeroCrossingLoc, XNeighbor[2], YNeighbor[2], ZNeighbor[2];
	double		XCoor_d, YCoor_d, ZCoor_d, GradVec_d[3], CurrLoc_d[3], Tempd;
	double		ZeroCrossingLoc_d[3], FirstDAtTheLoc_d, DataPosFromZeroCrossingLoc_d;
	double		Increase_d;


	XYZCoor_i[2] = CenterLoc / WtimesH_mi;
	XYZCoor_i[1] = (CenterLoc - XYZCoor_i[2]*WtimesH_mi) / Height_mi;
	XYZCoor_i[0] = CenterLoc % Width_mi;

	// Threshold Values
	// 1/8  --> 1/64
	// 1/32 --> 1/1024
	for (l=0; l<=1; l++) {
		XNeighbor[l] = 0;
		YNeighbor[l] = 0;
		ZNeighbor[l] = 0;
	}
	Increase_d = ZCVolumeGridSize_md;
	for (ZCoor_d=(double)XYZCoor_i[2]; ZCoor_d<(double)XYZCoor_i[2]+1.0; ZCoor_d+=Increase_d) {
		for (YCoor_d=(double)XYZCoor_i[1]; YCoor_d<(double)XYZCoor_i[1]+1.0; YCoor_d+=Increase_d) {
			for (XCoor_d=(double)XYZCoor_i[0]; XCoor_d<(double)XYZCoor_i[0]+1.0; XCoor_d+=Increase_d) {
			
				if (ZCVolumeGridSize_md-1e-5 <= fabs(XCoor_d-XYZCoor_i[0]) && 
					fabs(XCoor_d-XYZCoor_i[0]) <= 1.0-ZCVolumeGridSize_md-1e-5 &&
					ZCVolumeGridSize_md-1e-5 <= fabs(YCoor_d-XYZCoor_i[1]) && 
					fabs(YCoor_d-XYZCoor_i[1]) <= 1.0-ZCVolumeGridSize_md-1e-5 &&
					ZCVolumeGridSize_md-1e-5 <= fabs(ZCoor_d-XYZCoor_i[2]) && 
					fabs(ZCoor_d-XYZCoor_i[2]) <= 1.0-ZCVolumeGridSize_md-1e-5) continue;
			
				CurrLoc_d[0] = XCoor_d+ZCVolumeGridSize_md/2.0;
				CurrLoc_d[1] = YCoor_d+ZCVolumeGridSize_md/2.0;
				CurrLoc_d[2] = ZCoor_d+ZCVolumeGridSize_md/2.0;
				GradVecInterpolation(CurrLoc_d, GradVec_d);
				Tempd = sqrt(GradVec_d[2]*GradVec_d[2] + GradVec_d[1]*GradVec_d[1] + GradVec_d[0]*GradVec_d[0]);
				if (fabs(Tempd)<1e-5) continue;
				for (k=0; k<3; k++) GradVec_d[k] /= Tempd; // Normalize the gradient vector
				//Find Zero Connected Neighbors
				FoundZeroCrossingLoc = FindZeroCrossingLocation(CurrLoc_d, GradVec_d, ZeroCrossingLoc_d, 
											FirstDAtTheLoc_d, DataPosFromZeroCrossingLoc_d, 0.2);
				if (FoundZeroCrossingLoc && fabs(DataPosFromZeroCrossingLoc_d)<=1/ZCVolumeGridSize_md/32.0) {
					if (fabs(XCoor_d-XYZCoor_i[0])<1e-5) XNeighbor[0] = -1; // Negative Direction
					if (fabs(XCoor_d-XYZCoor_i[0])>=(1.0-ZCVolumeGridSize_md-1e-5)) XNeighbor[1] = 1; // Positive Direction
					
					if (fabs(YCoor_d-XYZCoor_i[1])<1e-5) YNeighbor[0] = -1;
					if (fabs(YCoor_d-XYZCoor_i[1])>=(1.0-ZCVolumeGridSize_md-1e-5)) YNeighbor[1] = 1;
					
					if (fabs(ZCoor_d-XYZCoor_i[2])<1e-5) ZNeighbor[0] = -1;
					if (fabs(ZCoor_d-XYZCoor_i[2])>=(1.0-ZCVolumeGridSize_md-1e-5)) ZNeighbor[1] = 1;

				}
			}
		}
	}

	int		XCoor_i, YCoor_i, ZCoor_i;
	map<int, unsigned char> Neighbor_map;
	Neighbor_map.clear();
	Neighbor_map[CenterLoc] = (unsigned char)0; // The Given Location
	
// Considering only 6 neighbors
	for (k=0; k<=1; k++) {
		XCoor_i = XYZCoor_i[0] + XNeighbor[k];
		YCoor_i = XYZCoor_i[1];
		ZCoor_i = XYZCoor_i[2];
		loc[0] = ZCoor_i*WtimesH_mi + YCoor_i*Width_mi + XCoor_i;
		Neighbor_map[loc[0]] = (unsigned char)0;
	}
	for (k=0; k<=1; k++) {
		XCoor_i = XYZCoor_i[0];
		YCoor_i = XYZCoor_i[1] + YNeighbor[k];
		ZCoor_i = XYZCoor_i[2];
		loc[0] = ZCoor_i*WtimesH_mi + YCoor_i*Width_mi + XCoor_i;
		Neighbor_map[loc[0]] = (unsigned char)0;
	}
	for (k=0; k<=1; k++) {
		XCoor_i = XYZCoor_i[0];
		YCoor_i = XYZCoor_i[1];
		ZCoor_i = XYZCoor_i[2] + ZNeighbor[k];
		loc[0] = ZCoor_i*WtimesH_mi + YCoor_i*Width_mi + XCoor_i;
		Neighbor_map[loc[0]] = (unsigned char)0;
	}

	Neighbor_map.erase(CenterLoc);
		
	map<int, unsigned char>::iterator Neighbors_it  = Neighbor_map.begin();
	
	for (i=0; i<(int)Neighbor_map.size(); Neighbors_it++, i++) {
		Neighbors[i] = 	(*Neighbors_it).first;
	}
	

#ifdef		SAVE_ZERO_CROSSING_SECOND_DERIVATIVE

//	int		n, m;
//	if (XYZCoor_i[2]==25) {
	
		cout << "Save Zero Crossing Second Derivative: XYZCoor_i = ";
		cout << XYZCoor_i[0] << " " << XYZCoor_i[1] << " " << XYZCoor_i[2] << endl;

		int 	j, GridSize_i, binfile_fd2;
		char	FileName[500];
		GridSize_i = (int)(floor(1.0/ZCVolumeGridSize_md+0.5));
		for (i=0; i<GridSize_i*GridSize_i*GridSize_i; i++) {
			VoxelVolume_muc[i] = (unsigned char)0;
		}

		for (k=0,ZCoor_d=(double)XYZCoor_i[2];ZCoor_d<(double)XYZCoor_i[2]+1.0-1e-5;ZCoor_d+=ZCVolumeGridSize_md,k++){
			for (j=0,YCoor_d=(double)XYZCoor_i[1];YCoor_d<(double)XYZCoor_i[1]+1.0-1e-5;YCoor_d+=ZCVolumeGridSize_md,j++){
				for (i=0,XCoor_d=(double)XYZCoor_i[0];XCoor_d<(double)XYZCoor_i[0]+1.0-1e-5;XCoor_d+=ZCVolumeGridSize_md,i++){

					CurrLoc_d[0] = XCoor_d + ZCVolumeGridSize_md/2.0;
					CurrLoc_d[1] = YCoor_d + ZCVolumeGridSize_md/2.0;
					CurrLoc_d[2] = ZCoor_d + ZCVolumeGridSize_md/2.0;
					GradVecInterpolation(CurrLoc_d, GradVec_d);
					Tempd = sqrt(GradVec_d[2]*GradVec_d[2] + GradVec_d[1]*GradVec_d[1] + GradVec_d[0]*GradVec_d[0]);
					if (fabs(Tempd)<1e-5) continue;
					for (int l=0; l<3; l++) GradVec_d[l] /= Tempd; // Normalize the gradient vector
					
					//Find Zero Crossing Neighbors
					FoundZeroCrossingLoc = FindZeroCrossingLocation(CurrLoc_d, GradVec_d, ZeroCrossingLoc_d, 
												FirstDAtTheLoc_d, DataPosFromZeroCrossingLoc_d, 0.2);
//												FirstDAtTheLoc_d, DataPosFromZeroCrossingLoc_d, ZCVolumeGridSize_md/20.0);


#ifdef		DEBUG_TF

					printf ("ijk = %3d %3d %3d ", i, j, k);
					printf ("CurrLoc = %8.4f %8.4f %8.4f ", CurrLoc_d[0], CurrLoc_d[1], CurrLoc_d[2]);
					printf ("GradVec = %6.2f %6.2f %6.2f ", GradVec_d[0], GradVec_d[1], GradVec_d[2]);
					printf ("ZeroLoc = %8.4f %8.4f %8.4f ", ZeroCrossingLoc_d[0], ZeroCrossingLoc_d[1], ZeroCrossingLoc_d[2]);
					if  (fabs(DataPosFromZeroCrossingLoc_d)<=ZCVolumeGridSize_md) printf ("C=255 ");
					else if  (fabs(DataPosFromZeroCrossingLoc_d)<=ZCVolumeGridSize_md*2) printf ("C=200 ");
					else if  (fabs(DataPosFromZeroCrossingLoc_d)<=ZCVolumeGridSize_md*4) printf ("C=150 ");
					else if  (fabs(DataPosFromZeroCrossingLoc_d)<=ZCVolumeGridSize_md*8) printf ("C=100 ");
					else if  (fabs(DataPosFromZeroCrossingLoc_d)<=ZCVolumeGridSize_md*16) printf ("C=50 ");
					printf ("ZeroPos = %20.16f\n", DataPosFromZeroCrossingLoc_d);

#endif

					if (FoundZeroCrossingLoc) {
						loc[0] = k*GridSize_i*GridSize_i + j*GridSize_i + i;
						if  (fabs(DataPosFromZeroCrossingLoc_d)<=ZCVolumeGridSize_md*16) 
							VoxelVolume_muc[loc[0]] = (unsigned char)50;
						if  (fabs(DataPosFromZeroCrossingLoc_d)<=ZCVolumeGridSize_md*8) 
							VoxelVolume_muc[loc[0]] = (unsigned char)100;
						if  (fabs(DataPosFromZeroCrossingLoc_d)<=ZCVolumeGridSize_md*4) 
							VoxelVolume_muc[loc[0]] = (unsigned char)150;
						if  (fabs(DataPosFromZeroCrossingLoc_d)<=ZCVolumeGridSize_md*2) 
							VoxelVolume_muc[loc[0]] = (unsigned char)200;
						if  (fabs(DataPosFromZeroCrossingLoc_d)<=ZCVolumeGridSize_md) 
							VoxelVolume_muc[loc[0]] = (unsigned char)255;
					}
				}
			}
		}

		sprintf (FileName, "Voxel_%03d_%03d_%03d.raw", XYZCoor_i[0], XYZCoor_i[1], XYZCoor_i[2]);
		cout << "Voxel File Name = " << FileName << endl;
		cout.flush();


		if ((binfile_fd2 = open (FileName, O_CREAT | O_WRONLY)) < 0) {
			printf ("could not open %s\n", FileName);
			exit(1);
		}
		if (write(binfile_fd2, VoxelVolume_muc, sizeof(unsigned char)*GridSize_i*GridSize_i*GridSize_i)
					!=(unsigned int)sizeof(unsigned char)*GridSize_i*GridSize_i*GridSize_i) {
			printf ("The file could not be written : %s\n", FileName);
			close (binfile_fd2);
			exit(1);
		}
		close (binfile_fd2);
		if (chmod(FileName, S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH) < 0) {
			printf ("chmod was not worked to file %s\n", FileName);
			exit(1);
		}
//	}


#endif
	
	int NumNeighbors = (int)Neighbor_map.size();
	Neighbor_map.clear();
	return NumNeighbors;

}



// Find 8 Neighbors for 2D Case, 26 Neighbors for 3D Case
template<class _DataType>
int cTFGeneration<_DataType>::FindNeighbors8Or26(int CenterLoc, int *Neighbors)
{
	int		i, XYZCoor_i[3];
	int		Neighbors26[26*3];
	
	XYZCoor_i[2] = CenterLoc / WtimesH_mi;
	XYZCoor_i[1] = (CenterLoc - XYZCoor_i[2]*WtimesH_mi) / Height_mi;
	XYZCoor_i[0] = CenterLoc % Width_mi;
	
	int NumNeighbor = FindNeighbors8Or26(XYZCoor_i, Neighbors26);
	
	for (i=0; i<NumNeighbor; i++) {
		Neighbors[i] = 	Neighbors26[i*3 + 2]*WtimesH_mi + 
						Neighbors26[i*3 + 1]*Width_mi + 
						Neighbors26[i*3 + 0];
	}
	
	return NumNeighbor;
}

// Find 8 Neighbors for 2D Case, 26 Neighbors for 3D Case
template<class _DataType>
int cTFGeneration<_DataType>::FindNeighbors8Or26(int *CenterLoc, int *Neighbors)
{
	int		i, j, k, NumNeighbor;


	if (Depth_mi<=1) {
		NumNeighbor=0;
		for (j=CenterLoc[1]-1; j<=CenterLoc[1]+1; j++) {
			if (j<0 || j>=Height_mi) continue;
			for (k=CenterLoc[0]-1; k<=CenterLoc[0]+1; k++) {
				if (k<0 || k>=Width_mi) continue;
				if (k==CenterLoc[0] && j==CenterLoc[1]) continue;
				Neighbors[NumNeighbor*3 + 0] = CenterLoc[0]+k;	// X
				Neighbors[NumNeighbor*3 + 1] = CenterLoc[1]+j;	// Y
				Neighbors[NumNeighbor*3 + 2] = 0;					// Z
				NumNeighbor++;
			}
		}
	}
	else {
		NumNeighbor=0;
		for (i=CenterLoc[2]-1; i<=CenterLoc[2]+1; i++) {
			if (i<0 || i>=Depth_mi) continue;
			for (j=CenterLoc[1]-1; j<=CenterLoc[1]+1; j++) {
				if (j<0 || j>=Height_mi) continue;
				for (k=CenterLoc[0]-1; k<=CenterLoc[0]+1; k++) {
					if (k<0 || k>=Width_mi) continue;
					if (k==CenterLoc[0] && j==CenterLoc[1] && i==CenterLoc[2]) continue;
					Neighbors[NumNeighbor*3 + 0] = k;    // X
					Neighbors[NumNeighbor*3 + 1] = j;    // Y
					Neighbors[NumNeighbor*3 + 2] = i;    // Z
					NumNeighbor++;
				}
			}
		}
	}
	
	return NumNeighbor;
}


// Find 4 Neighbors for 2D Case, 6 Neighbors for 3D Case
template<class _DataType>
void cTFGeneration<_DataType>::FindNeighbors4Or6(int *CenterLoc, int *Neighbors)
{
	int		l;


	if (Depth_mi<=1) {
		l=0;
		Neighbors[l*3 + 0] = CenterLoc[0]-1;	// X
		Neighbors[l*3 + 1] = CenterLoc[1];		// Y
		Neighbors[l*3 + 2] = 0;					// Z
		l++;
		Neighbors[l*3 + 0] = CenterLoc[0]+1;	// X
		Neighbors[l*3 + 1] = CenterLoc[1];		// Y
		Neighbors[l*3 + 2] = 0;					// Z
		l++;
		Neighbors[l*3 + 0] = CenterLoc[0];		// X
		Neighbors[l*3 + 1] = CenterLoc[1]-1;	// Y
		Neighbors[l*3 + 2] = 0;					// Z
		l++;
		Neighbors[l*3 + 0] = CenterLoc[0];		// X
		Neighbors[l*3 + 1] = CenterLoc[1]+1;	// Y
		Neighbors[l*3 + 2] = 0;					// Z
	}
	else {
		l=0;
		Neighbors[l*3 + 0] = CenterLoc[0]-1;	// X
		Neighbors[l*3 + 1] = CenterLoc[1];		// Y
		Neighbors[l*3 + 2] = CenterLoc[2];		// Z
		l++;
		Neighbors[l*3 + 0] = CenterLoc[0]+1;	// X
		Neighbors[l*3 + 1] = CenterLoc[1];		// Y
		Neighbors[l*3 + 2] = CenterLoc[2];		// Z
		l++;
		Neighbors[l*3 + 0] = CenterLoc[0];		// X
		Neighbors[l*3 + 1] = CenterLoc[1]-1;	// Y
		Neighbors[l*3 + 2] = CenterLoc[2];		// Z
		l++;
		Neighbors[l*3 + 0] = CenterLoc[0];		// X
		Neighbors[l*3 + 1] = CenterLoc[1]+1;	// Y
		Neighbors[l*3 + 2] = CenterLoc[2];		// Z
		l++;
		Neighbors[l*3 + 0] = CenterLoc[0];		// X
		Neighbors[l*3 + 1] = CenterLoc[1];		// Y
		Neighbors[l*3 + 2] = CenterLoc[2]-1;	// Z
		l++;
		Neighbors[l*3 + 0] = CenterLoc[0];		// X
		Neighbors[l*3 + 1] = CenterLoc[1];		// Y
		Neighbors[l*3 + 2] = CenterLoc[2]+1;	// Z
	}
}



template<class _DataType>
void cTFGeneration<_DataType>::SaveSecondDerivative(char *OutFileName, int NumLeaping)
{
	FILE			*out_st;
	char			SecondDFileName[200];
	int				i, k, l, loc[5], RColor, GColor, BColor, DataCoor[3];
	double			GVec[3], GradientDir[3], GradM[3], Step, Length_d, SecondD[3];


	if (NumLeaping > 1) k=NumLeaping;
	else k=0;
	for (; k<Depth_mi; k+=NumLeaping) {	// Jump to next slice

		//-----------------------------------------------------------------------
		// Second Derivative
		if (Depth_mi==1 && NumLeaping==1) {
			sprintf (SecondDFileName, "%sSecondD.ppm", OutFileName);
		}
		else sprintf (SecondDFileName, "%s%04dSecondD.ppm", OutFileName, k);

		if ((k==0 || k==Depth_mi-1) && Depth_mi>1) {  }
		else {
			printf ("Second Derivative File = %s\n", SecondDFileName); fflush (stdout);
			if ((out_st=fopen(SecondDFileName, "w"))==NULL) {
				printf ("Could not open %s\n", SecondDFileName);
				exit(1);
			}
			fprintf (out_st, "P3\n%d %d\n", Width_mi, Height_mi);
			fprintf (out_st, "%d\n", 255);
			
			for (i=0; i<WtimesH_mi; i++) {
				loc[0] = k*WtimesH_mi + i;
				if (SecondDerivative_mf[loc[0]]>0) {
					RColor = (int)(((double)SecondDerivative_mf[loc[0]]-0)/((double)MaxSecMag_mf-0)*255.0);
					BColor = 0;
				}
				else {
					BColor = (int)((fabs((double)SecondDerivative_mf[loc[0]])-0)/(fabs((double)MinSecMag_mf)-0)*255.0);
					RColor=0;
				}

				if (RColor < 0) RColor = 0;
				if (RColor > 255) RColor = 255;
				if (BColor < 0) BColor = 0;
				if (BColor > 255) BColor = 255;
				
				// Finding zero second derivatives and local maximum gradient locations				
				loc[1] = k*WtimesH_mi + i;
				for (l=0; l<3; l++) GVec[l] = (double)GradientVec_mf[loc[1]*3 + l];
				Length_d = sqrt (GVec[0]*GVec[0] + GVec[1]*GVec[1] + GVec[2]*GVec[2]);
				if (fabs(Length_d)<1e-5) {
					GColor = 0;
				}
				else {
					for (l=0; l<3; l++) GVec[l] /= Length_d; // Normalize the gradient vector

					DataCoor[2] = k;
					DataCoor[1] = i / Width_mi;
					DataCoor[0] = i % Width_mi;

					GColor=0;

					Step=-1.0;
					for (l=0; l<3; l++) GradientDir[l] = (double)DataCoor[l] + GVec[l]*Step;
					GradM[0] = GradientInterpolation(GradientDir[0], GradientDir[1], GradientDir[2]);
					SecondD[0] = SecondDInterpolation(GradientDir[0], GradientDir[1], GradientDir[2]);
					
					Step=-0.75;
					for (l=0; l<3; l++) GradientDir[l] = (double)DataCoor[l] + GVec[l]*Step;
					GradM[1] = GradientInterpolation(GradientDir[0], GradientDir[1], GradientDir[2]);
					SecondD[1] = SecondDInterpolation(GradientDir[0], GradientDir[1], GradientDir[2]);
					for (Step=-0.50; Step<=1.0;Step+=0.25) {

						for (l=0; l<3; l++) GradientDir[l] = (double)DataCoor[l] + GVec[l]*Step;
						GradM[2] = GradientInterpolation(GradientDir[0], GradientDir[1], GradientDir[2]);
						SecondD[2] = SecondDInterpolation(GradientDir[0], GradientDir[1], GradientDir[2]);
						
						// Threshold Value = 10 for gradient magnitudes
						if (SecondD[0]*SecondD[2]<0 && GradM[1]>2) {
							GColor=155; // Only Zero Crossing Second Derivatives with Threshold Values
							if (GradM[1]>=GradM[0] && GradM[1]>=GradM[2]) {
								GColor=255; // Local Maximum Gradient Magnitude with Zero Crossing Second Derivatives
								break;
							}
						}
						for (l=0; l<=1; l++) {
							GradM[l] = GradM[l+1];
							SecondD[l] = SecondD[l+1];
						}
					}
				}


//				loc[1] = k*WtimesH_mi + i;
//				if (fabs((double)SecondDerivative_mf[loc[1]])<1e-1 && GradientMag_mf[loc[1]]>0) G=255;
//				else G=0;

				fprintf (out_st, "%d %d %d\n", RColor, GColor, BColor);
			}
			fclose(out_st);
		}
	}
	
}


template <class _DataType>
void cTFGeneration<_DataType>::SaveBoundaryVolume(char *OutFileName, _DataType MatMin, _DataType MatMax)
{
	int     		i, binfile_fd2;
	char			MaterialVolumeName[200];
	unsigned char	MinData=255, MaxData=0;


	sprintf (MaterialVolumeName, "%s_%03d-%03d.raw", OutFileName, (int)MatMin, (int)MatMax);
	cout << "Boundary Volume File Name = " << MaterialVolumeName << endl;

	for (i=0; i<WHD_mi; i++) {
		if ((VoxelStatus_muc[i] & MAT_ZERO_CROSSING)) {
			VoxelStatus_muc[i] = (unsigned char)Data_mT[i];
			if (MinData > VoxelStatus_muc[i]) MinData = VoxelStatus_muc[i];
			if (MaxData < VoxelStatus_muc[i]) MaxData = VoxelStatus_muc[i];
		}
		else VoxelStatus_muc[i] = (unsigned char)0;
	}

	if ((binfile_fd2 = open (MaterialVolumeName, O_CREAT | O_WRONLY)) < 0) {
		printf ("could not open %s\n", MaterialVolumeName);
		exit(1);
	}
	if (write(binfile_fd2, VoxelStatus_muc, sizeof(unsigned char)*WHD_mi)!=
		(unsigned int)sizeof(unsigned char)*WHD_mi) {
		printf ("The file could not be written : %s\n", MaterialVolumeName);
		close (binfile_fd2);
		exit(1);
	}
	close (binfile_fd2);
	if (chmod(MaterialVolumeName, S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH) < 0) {
		printf ("chmod was not worked to file %s\n", MaterialVolumeName);
		exit(1);
	}
	
	cout << endl;
	cout << "Min & Max Value of Extracted Material Volume = " << (int)MinData << " " << (int)MaxData << endl;
	cout << endl;

}

template <class _DataType>
void cTFGeneration<_DataType>::SaveAlphaVolume(char *OutFileName, _DataType MatMin, _DataType MatMax)
{
	int     		i, binfile_fd2;
	char			MaterialVolumeName[200];
	unsigned char	MinData=255, MaxData=0;
	unsigned char	*AlphaVolume_muc = new unsigned char[WHD_mi];
	

	sprintf (MaterialVolumeName, "%s_%03d-%03dAlpha.raw", OutFileName, (int)MatMin, (int)MatMax);
	cout << "Alpha Volume File Name = " << MaterialVolumeName << endl;

	for (i=0; i<WHD_mi; i++) {

		AlphaVolume_muc[i] = (unsigned char)(AlphaVolume_mf[i]*255.0);
		if (VoxelStatus_muc[i]==0 && MatMin<=Data_mT[i] && Data_mT[i]<=MatMax) {
			AlphaVolume_muc[i] = (unsigned char)255;
		}
	}
	
	if ((binfile_fd2 = open (MaterialVolumeName, O_CREAT | O_WRONLY)) < 0) {
		printf ("could not open %s\n", MaterialVolumeName);
		exit(1);
	}
	if (write(binfile_fd2, AlphaVolume_muc, sizeof(unsigned char)*WHD_mi)!=
		(unsigned int)sizeof(unsigned char)*WHD_mi) {
		printf ("The file could not be written : %s\n", MaterialVolumeName);
		close (binfile_fd2);
		exit(1);
	}
	close (binfile_fd2);
	if (chmod(MaterialVolumeName, S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH) < 0) {
		printf ("chmod was not worked to file %s\n", MaterialVolumeName);
		exit(1);
	}
	
	cout << endl;
	cout << "The Min & Max Values of the Alpha Volume = " << (int)MinData << " " << (int)MaxData << endl;
	cout << endl;

	delete [] AlphaVolume_muc;
}


template <class _DataType>
void cTFGeneration<_DataType>::RemoveSecondD()
{
	delete [] SecondDerivative_mf;
	SecondDerivative_mf = NULL;
}

template <class _DataType>
void cTFGeneration<_DataType>::Destroy()
{
	for (int i=0; i<NumMaterial_mi; i++){
		 delete [] TF_mi[i];
		 TF_mi[i] = NULL;
	}
	delete [] H_vg_md;
	H_vg_md = NULL;
	delete [] G_v_md;
	G_v_md = NULL;
	delete [] CCIndexTable_mi;
	CCIndexTable_mi = NULL;
}


template<class _DataType>
void cTFGeneration<_DataType>::DisplayIGS(int X1, int X2, int Y1, int Y2, int Z1, int Z2, float *DistanceV)
{
	int		i, j, k, l, m, n, loc[8];
	int		Width, Height, Depth, WtimesH, WHD;
	float	X_f, Y_f, XG_f, YG_f, Tempf;
	int		Min_i, Max_i;
	float	MinGradM_f, MaxGradM_f, MinSecondD_f, MaxSecondD_f;
	float	MinDist_f, MaxDist_f;
	

	Width = X2 - X1 + 1;
	Height = Y2 - Y1 + 1;
	Depth = Z2 - Z1 + 1;
	WtimesH = Width*Height;
	WHD = Width*Height*Depth;
	
	int		*Intensity_i = new int [WHD];
	float	*Gradient_Mag_f = new float [WHD];
	float	*SecondD_f = new float [WHD];
	

	printf ("Ranges X Y Z = (%d, %d), (%d, %d), (%d, %d)\n", X1, X2, Y1, Y2, Z1, Z2);

	Min_i = 9999999;
	Max_i = -9999999;
	printf ("Intensities\n");
	for (n=0, k=Z1; k<=Z2; k++, n++) {
		printf ("Z = %d\n", k);
		for (m=0, j=Y1; j<=Y2; j++, m++) {
			for (l=0, i=X1; i<=X2; i++, l++) {
			
				loc[0] = k*WtimesH_mi + j*Width_mi + i;
				printf ("%4d ", (int)Data_mT[loc[0]]);

				loc[1] = n*WtimesH + m*Width + l;
				Intensity_i[loc[1]] = (int)Data_mT[loc[0]];
				
				if (Min_i > Intensity_i[loc[1]]) Min_i = Intensity_i[loc[1]];
				if (Max_i < Intensity_i[loc[1]]) Max_i = Intensity_i[loc[1]];
			}
			printf ("\n");
		}
		printf ("\n");
	}
	printf ("Min Max Data Value = %d %d\n", Min_i, Max_i);
	printf ("\n");
	fflush(stdout);
	Intensity_i[0] = 0;
	
	printf ("Gradient Magnitudes\n");
	MinGradM_f = 9999999.9999;
	MaxGradM_f = -9999999.9999;
	for (n=0; n<Depth; n++) {
		printf ("Z = %d\n", Z1);
		for (m=0; m<Height; m++) {
			for (l=0; l<Width; l++) {
			
				loc[0] = n*WtimesH + m*Width + l;
				
				if (l<=0) loc[1] = 0;
				else loc[1] = n*WtimesH + m*Width + l-1;
				if (l>=Width-1) loc[2] = 0;
				else loc[2] = n*WtimesH + m*Width + l+1;
				
				if (m<=0) loc[3] = 0;
				else loc[3] = n*WtimesH + (m-1)*Width + l;
				if (m>=Height-1) loc[4] = 0;
				else loc[4] = n*WtimesH + (m+1)*Width + l;
				
				X_f = Intensity_i[loc[2]] - Intensity_i[loc[1]];
				Y_f = Intensity_i[loc[4]] - Intensity_i[loc[3]];
				Gradient_Mag_f[loc[0]] = sqrt(X_f*X_f + Y_f*Y_f);
				if (MinGradM_f > Gradient_Mag_f[loc[0]]) MinGradM_f = Gradient_Mag_f[loc[0]];
				if (MaxGradM_f < Gradient_Mag_f[loc[0]]) MaxGradM_f = Gradient_Mag_f[loc[0]];

				printf ("%7.2f ", Gradient_Mag_f[loc[0]]);
				
			}
			printf ("\n");
		}
		printf ("\n");
	}
	printf ("Min Max Gradient Magnitude = %f %f\n", MinGradM_f, MaxGradM_f);
	printf ("\n");
	fflush(stdout);
	Gradient_Mag_f[0] = (float)0.0;
	
	
	printf ("Second Derivatives\n");
	MinSecondD_f = 999999.9999;
	MaxSecondD_f = -999999.9999;
	for (n=0; n<Depth; n++) {
		printf ("Z = %d\n", Z1);
		for (m=0; m<Height; m++) {
			for (l=0; l<Width; l++) {
			
				loc[0] = n*WtimesH + m*Width + l;
				
				if (l<=0) loc[1] = 0;
				else loc[1] = n*WtimesH + m*Width + l-1;
				if (l>=Width-1) loc[2] = 0;
				else loc[2] = n*WtimesH + m*Width + l+1;
				
				if (m<=0) loc[3] = 0;
				else loc[3] = n*WtimesH + (m-1)*Width + l;
				if (m>=Height-1) loc[4] = 0;
				else loc[4] = n*WtimesH + (m+1)*Width + l;
				
				XG_f = Intensity_i[loc[2]] - Intensity_i[loc[1]];
				YG_f = Intensity_i[loc[4]] - Intensity_i[loc[3]];

				Tempf = sqrt(XG_f*XG_f + YG_f*YG_f);
				if (fabs(Tempf)<1e-7) {
					XG_f = 0.0;
					YG_f = 0.0;
				}
				else {
					XG_f /= Tempf;
					YG_f /= Tempf;
				}

				X_f = (Gradient_Mag_f[loc[2]] - Gradient_Mag_f[loc[1]])*XG_f;
				Y_f = (Gradient_Mag_f[loc[4]] - Gradient_Mag_f[loc[3]])*YG_f;
				SecondD_f[loc[0]] = X_f + Y_f;
				if (MinSecondD_f > SecondD_f[loc[0]]) MinSecondD_f = SecondD_f[loc[0]];
				if (MaxSecondD_f < SecondD_f[loc[0]]) MaxSecondD_f = SecondD_f[loc[0]];
				
				if (Gradient_Mag_f[loc[0]]>25) printf ("%9.4f ", SecondD_f[loc[0]]);
				else printf ("%9.4f ", 0.0);
				
			}
			printf ("\n");
		}
		printf ("\n");
	}
	printf ("Min Max Second D = %f %f\n", MinSecondD_f, MaxSecondD_f);
	printf ("\n");
	fflush(stdout);


	printf ("Distance From Zero-Crossing Locations\n");
	MinDist_f = 9999999.9999;
	MaxDist_f = -9999999.9999;
	for (n=0; n<Depth; n++) {
		for (m=0; m<Height; m++) {
			for (l=0; l<Width; l++) {
				loc[0] = n*WtimesH + m*Width + l;
				if (MinDist_f>DistanceV[loc[0]]) MinDist_f = DistanceV[loc[0]];
				if (MaxDist_f<DistanceV[loc[0]]) MaxDist_f = DistanceV[loc[0]];
			}
		}
	}
	for (n=0; n<Depth; n++) {
		printf ("Z = %d\n", Z1);
		for (m=0; m<Height; m++) {
			for (l=0; l<Width; l++) {
			
				loc[0] = n*WtimesH + m*Width + l;
				if (fabsf(DistanceV[loc[0]])>15.0) printf ("%7.2f ", 0.0);
				else printf ("%7.2f ", 15.0 - fabs(DistanceV[loc[0]]));

//				else if (DistanceV[loc[0]]>=0.0) printf ("%7.2f ", 256.0 - DistanceV[loc[0]]*256.0/16.0);
//				else printf ("%7.2f ", (256.0 - fabsf(DistanceV[loc[0]])*256.0/16.0)*(-1.0));
				
			}
			printf ("\n");
		}
		printf ("\n");
	}
	printf ("Min Max Gradient Magnitude = %f %f\n", MinGradM_f, MaxGradM_f);
	printf ("\n");
	fflush(stdout);
	Gradient_Mag_f[0] = (float)0.0;


	int		NumBinsI, NumBinsG, NumBinsS, NumBins;
	int		*HistogramAll;
	double	HistogramFactorI, HistogramFactorG, HistogramFactorS;
	

	HistogramFactorI = (double)256.0/(Max_i - Min_i);
	HistogramFactorG = (double)256.0/(MaxGradM_f - MinGradM_f);
	HistogramFactorS = (double)256.0/(MaxSecondD_f - MinSecondD_f);
	
	NumBinsI = (int)((Max_i - Min_i)*HistogramFactorI);
	NumBinsG = (int)((MaxGradM_f - MinGradM_f)*HistogramFactorG);
	NumBinsS = (int)((MaxSecondD_f - MinSecondD_f)*HistogramFactorS);
	
	if (NumBinsI > NumBinsG) NumBins = NumBinsI;
	else NumBins = NumBinsG;
	if (NumBins < NumBinsS) NumBins = NumBinsS;

	printf ("Num Bins I G S = %d %d %d\n", NumBinsI, NumBinsG, NumBinsS);
	printf ("Num Bins = %d\n", NumBins);
	fflush(stdout);

	HistogramAll = new int[NumBins*NumBins];

	// Making the I-G Graph
	for (i=0; i<NumBins*NumBins; i++) HistogramAll[i]=0;
	for (i=0; i<WHD; i++) {
		loc[0] = (int)(((double)Intensity_i[i] - Min_i)*HistogramFactorI);		// Intensity Values
		loc[1] = (int)(((double)Gradient_Mag_f[i] - MinGradM_f)*HistogramFactorG);
		
		HistogramAll[loc[1]*NumBins + loc[0]]++;
	}

	printf ("I-G Graph General\n");
	for (j=0; j<NumBins; j++) {
		for (i=0; i<NumBins; i++) {
			loc[0] = j*NumBins + i;
			printf ("%5d ", HistogramAll[loc[0]]);
		}
		printf ("\n");
	}
	printf ("\n");
	fflush(stdout);
	
	// Making the I-S Graph
	for (i=0; i<NumBins*NumBins; i++) HistogramAll[i]=0;
	for (i=0; i<WHD; i++) {
		loc[0] = (int)(((double)Intensity_i[i] - Min_i)*HistogramFactorI);		// Intensity Values
		loc[1] = (int)(((double)SecondD_f[i] - MinSecondD_f)*HistogramFactorS);// Gradient Magnitudes
		HistogramAll[loc[1]*NumBins + loc[0]]++;
	}

	printf ("I-S Graph General\n");
	for (j=0; j<NumBins; j++) {
		for (i=0; i<NumBins; i++) {
			loc[0] = j*NumBins + i;
			printf ("%5d ", HistogramAll[loc[0]]);
		}
		printf ("\n");
	}
	printf ("\n");
	fflush(stdout);
	
	delete [] HistogramAll;

}


template<class _DataType>
int	cTFGeneration<_DataType>::Index(int X, int Y, int Z, int ith, int NumElements)
{
	if (X<0 || Y<0 || Z<0 || X>=Width_mi || Y>=Height_mi || Z>=Depth_mi) return 0;
	else return (Z*WtimesH_mi*NumElements + Y*Width_mi*NumElements + X*NumElements + ith);
}

template<class _DataType>
int	cTFGeneration<_DataType>::Index(int X, int Y, int Z)
{
	if (X<0 || Y<0 || Z<0 || X>=Width_mi || Y>=Height_mi || Z>=Depth_mi) return 0;
	else return (Z*WtimesH_mi + Y*Width_mi + X);
}


template<class _DataType>
void cTFGeneration<_DataType>::NormalizeVector(double *GradVector3)
{
	double	Length_d;
	Length_d = sqrt (GradVector3[0]*GradVector3[0] + GradVector3[1]*GradVector3[1] + GradVector3[2]*GradVector3[2]);
	GradVector3[0] /= Length_d;
	GradVector3[1] /= Length_d;
	GradVector3[2] /= Length_d;
}


//---------------------------------------------------------------------------------------
// Quick Sort Locations
//---------------------------------------------------------------------------------------

// Sort first Axis1, second Axis2, and third Axis3
template<class _DataType>
void cTFGeneration<_DataType>::QuickSort3Elements(int *DataArray3, int NumLocs, char Axis1, char Axis2, char Axis3)
{
	int		i, Curr_Loc=0, ith_Element=0, NumNextElements=0;


	if (NumLocs <= 1) return;

	// Sorting with Axis 1
	QuickSort3Elements(DataArray3, NumLocs, Axis1);

	// Sorting with Axis2
	switch (Axis1) { // Check whether the first Axis1 element is the same or not
		case 'X' : ith_Element = 0; break;
		case 'Y' : ith_Element = 1; break;
		case 'Z' : ith_Element = 2; break;
		default : break;
	}
	Curr_Loc = 0;
	do {
		NumNextElements = 0;
		for (i=Curr_Loc; i<NumLocs; i++) {
			if (DataArray3[Curr_Loc*3 + ith_Element]==DataArray3[i*3 + ith_Element]) NumNextElements++;
			else break;
		}
		QuickSort3Elements(&DataArray3[Curr_Loc*3], NumNextElements, Axis2);
		Curr_Loc += NumNextElements;
	} while (Curr_Loc < NumLocs);

	// Sorting with Axis3
	switch (Axis2) { // Check whether the first Axis1 element is the same or not
		case 'X' : ith_Element = 0; break;
		case 'Y' : ith_Element = 1; break;
		case 'Z' : ith_Element = 2; break;
		default : break;
	}
	Curr_Loc = 0;
	do {
		NumNextElements = 0;
		for (i=Curr_Loc; i<NumLocs; i++) {
			if (DataArray3[Curr_Loc*3 + ith_Element]==DataArray3[i*3 + ith_Element]) NumNextElements++;
			else break;
		}
		QuickSort3Elements(&DataArray3[Curr_Loc*3], NumNextElements, Axis3);
		Curr_Loc += NumNextElements;
	} while (Curr_Loc < NumLocs);

}


char	AxisTF_gc;
int		ith_ElementTF_gi;

template<class _DataType> 
void cTFGeneration<_DataType>::QuickSort3Elements(int *DataArray3, int NumLocs, char Axis)
{
	if (NumLocs <= 1) return;
	AxisTF_gc = Axis;
	switch (Axis) {
		case 'X' : ith_ElementTF_gi = 0; break;
		case 'Y' : ith_ElementTF_gi = 1; break;
		case 'Z' : ith_ElementTF_gi = 2; break;
		default : break;
	}
//	cout << "Axis = " << AxisTF_gc << endl;
//	cout << "Start Quick Sort " << endl;
	QuickSortLocs(DataArray3, 0, NumLocs-1);

}

template<class _DataType> 
void cTFGeneration<_DataType>::QuickSortLocs(int *DataArray3, int p, int r)
{
	int q;

	if (p<r) {
		q = PartitionLocs(DataArray3, p, r);
		QuickSortLocs(DataArray3, p, q-1);
		QuickSortLocs(DataArray3, q+1, r);
	}
}


template<class _DataType> 
void cTFGeneration<_DataType>::SwapLocs(int& x, int& y)
{
	int		Temp;

	Temp = x;
	x = y;
	y = Temp;
}


template<class _DataType> 
int cTFGeneration<_DataType>::PartitionLocs(int *DataArray3, int low, int high)
{
	int		left, right;
	int 	pivot_item, pivot_itemX, pivot_itemY, pivot_itemZ;
	

	pivot_item = DataArray3[low*3 + ith_ElementTF_gi];
	pivot_itemX = DataArray3[low*3 + 0];
	pivot_itemY = DataArray3[low*3 + 1];
	pivot_itemZ = DataArray3[low*3 + 2];
	
	left = low;
	right = high;

	while ( left < right ) {
		while( DataArray3[left*3 + ith_ElementTF_gi] <= pivot_item  && left<=high) left++;
		while( DataArray3[right*3 + ith_ElementTF_gi] > pivot_item && right>=low) right--;
		if ( left < right ) {
			SwapLocs(DataArray3[left*3 + 0], DataArray3[right*3 + 0]);
			SwapLocs(DataArray3[left*3 + 1], DataArray3[right*3 + 1]);
			SwapLocs(DataArray3[left*3 + 2], DataArray3[right*3 + 2]);
		}
	}

	DataArray3[low*3+0] = DataArray3[right*3+0];
	DataArray3[low*3+1] = DataArray3[right*3+1];
	DataArray3[low*3+2] = DataArray3[right*3+2];
	DataArray3[right*3+0] = pivot_itemX;
	DataArray3[right*3+1] = pivot_itemY;
	DataArray3[right*3+2] = pivot_itemZ;

	return right;
}









void _QuickSort_FunctionGenerator()
{
	int				data_i[] = {10, 20, 30};
	char			data_c[] = {10, 20, 30};
	unsigned char	data_uc[] = {10, 20, 30};
	unsigned short	data_us[] = {10, 20, 30};
	float			data_f[] = {10., 20., 30.};
	double			data_d[] = {10., 20., 30.};
	
	QuickSort(data_i, 0, 2);
	QuickSort(data_c, 0, 2);
	QuickSort(data_uc, 0, 2);
	QuickSort(data_us, 0, 2);
	QuickSort(data_f, 0, 2);
	QuickSort(data_d, 0, 2);
}


template<class _DataType>
void QuickSort(_DataType* data, int p, int r)
{
	if (p==0 && r==0) return;
	int q;

	if (p<r) {
		q = Partition(data, p, r);
		QuickSort(data, p, q-1);
		QuickSort(data, q+1, r);
	}
}


template<class _DataType> 
void Swap(_DataType& x, _DataType& y)
{
	_DataType	Temp;

	Temp = x;
	x = y;
	y = Temp;
}


template<class _DataType> 
int Partition(_DataType* data, int low, int high)
{
	int 		left, right;
	_DataType 	pivot_item;
	
	pivot_item = data[low];
	
	left = low;
	right = high;
	while ( left < right ) {
		while( data[left] <= pivot_item && left<=high) left++;
		while( data[right] > pivot_item && right>=low) right--;
		if ( left < right ) {
			Swap(data[left], data[right]);
		}
	}

	data[low] = data[right];
	data[right] = pivot_item;

	return right;
}


template class cTFGeneration<unsigned char>;
template class cTFGeneration<unsigned short>;
template class cTFGeneration<int>;
template class cTFGeneration<float>;










/*
// Marking Boundary Voxels
template <class _DataType>
void cPEDetection<_DataType>::VesselTracking(_DataType MatMin, _DataType MatMax, float *StartPt, float *EndPt)
{
	int			k, loc[3], DataCoor[3], ZeroCrossingLoc_i[3];
	double		GradVec[3], ZeroCrossingLoc_d[3], Tempd;
	double		LocalMaxGradient, DataPosFromZeroCrossingLoc_d;
	int			NumNeighbors, Neighbor26_i[26], FoundZeroCrossingLoc;
	
	double		VesselDirection_d[3], Weight_d, PointOnPlane_d[3];
	double		Ray1_d[3], Ray2_d[3], Ray3_d[3], Ray4_d[3];
	double		StartPt_d[3], Increase_d;


	// Compute a plane equation and a perpendicular ray direction to the vessel direction
	for (k=0; k<3; k++) VesselDirection_d[k] = (double)EndPt[k] - StartPt[k];
	Tempd = sqrt(	VesselDirection_d[0]*VesselDirection_d[0] + 
					VesselDirection_d[1]*VesselDirection_d[1] + 
					VesselDirection_d[2]*VesselDirection_d[2]	);
	for (k=0; k<3; k++) VesselDirection_d[k] /= Tempd; // Normalize the vessel direction
	Weight_d = -(StartPt[0]*VesselDirection_d[0] + StartPt[1]*VesselDirection_d[1] + StartPt[2]*VesselDirection_d[2]);
	
	PointOnPlane_d[0] = (double)StartPt[0] + 10.0;
	PointOnPlane_d[1] = (double)StartPt[1];
	if (fabs(VesselDirection_d[2])<1e-5) {
		PointOnPlane_d[2] = 0.0;
	}
	else {
		PointOnPlane_d[2] = (-PointOnPlane_d[0]*VesselDirection_d[0] - 
								PointOnPlane_d[1]*VesselDirection_d[1] - Weight_d)/VesselDirection_d[2];
	}

#ifdef	DEBUG_PED
	printf ("Start & End Points = (%5.2f %5.2f %5.2f) --> (%5.2f %5.2f %5.2f) ", StartPt[0], StartPt[1], StartPt[2], 
				EndPt[0], EndPt[1], EndPt[2]);
	printf ("Direction = (%f %f %f) ", VesselDirection_d[0], VesselDirection_d[1], VesselDirection_d[2]);
	printf ("W = %f\n", Weight_d);
	printf ("A Point on the Plane = (%f %f %f)\n", PointOnPlane_d[0], PointOnPlane_d[1], PointOnPlane_d[2]);
	fflush(stdout);
#endif

	// Compute Ray1
	for (k=0; k<3; k++) Ray1_d[k] = PointOnPlane_d[k] - StartPt[k];
	Tempd = sqrt(	Ray1_d[0]*Ray1_d[0] + Ray1_d[1]*Ray1_d[1] + Ray1_d[2]*Ray1_d[2]	);
	for (k=0; k<3; k++) Ray1_d[k] /= Tempd; // Normalize the ray

	// Compute Ray2 with the cross product
	Ray2_d[0] = VesselDirection_d[1]*Ray1_d[2]-VesselDirection_d[2]*Ray1_d[1];
	Ray2_d[1] = VesselDirection_d[2]*Ray1_d[0]-VesselDirection_d[0]*Ray1_d[2];
	Ray2_d[2] = VesselDirection_d[0]*Ray1_d[1]-VesselDirection_d[1]*Ray1_d[0];

	// Compute Ray3
	for (k=0; k<3; k++) Ray3_d[k] = (Ray1_d[k] + Ray2_d[k])/2.0;
	Tempd = sqrt(	Ray3_d[0]*Ray3_d[0] + Ray3_d[1]*Ray3_d[1] + Ray3_d[2]*Ray3_d[2]	);
	for (k=0; k<3; k++) Ray3_d[k] /= Tempd; // Normalize the ray

	// Compute Ray4 with the cross product
	Ray4_d[0] = VesselDirection_d[1]*Ray3_d[2]-VesselDirection_d[2]*Ray3_d[1];
	Ray4_d[1] = VesselDirection_d[2]*Ray3_d[0]-VesselDirection_d[0]*Ray3_d[2];
	Ray4_d[2] = VesselDirection_d[0]*Ray3_d[1]-VesselDirection_d[1]*Ray3_d[0];

#ifdef	DEBUG_PED
	printf ("Ray1-4: (%f %f %f), (%f %f %f), (%f %f %f), (%f %f %f)\n", 
			Ray1_d[0], Ray1_d[1], Ray1_d[2], Ray2_d[0], Ray2_d[1], Ray2_d[2], 
			Ray3_d[0], Ray3_d[1], Ray3_d[2], Ray4_d[0], Ray4_d[1], Ray4_d[2]);
	fflush(stdout);
#endif	


	map<int, unsigned char> InitialBoundaryLocs_map;
	map<int, unsigned char> ZeroCrossingLocs_map;
	InitialBoundaryLocs_map.clear();
	ZeroCrossingLocs_map.clear();

	int		DataLoc;

	Increase_d = 0.2;
	for (k=0; k<3; k++) StartPt_d[k] = (double)StartPt[k];
	// Ray 1
	DataLoc = getANearestBoundary(StartPt_d, Ray1_d, Increase_d, MatMin, MatMax);
	if (DataLoc>0) InitialBoundaryLocs_map[DataLoc]=(unsigned char)0; // Add it to the map
	DataLoc = getANearestBoundary(StartPt_d, Ray1_d, -Increase_d, MatMin, MatMax);
	if (DataLoc>0) InitialBoundaryLocs_map[DataLoc]=(unsigned char)0; // Add it to the map
	
	// Ray 2
	DataLoc = getANearestBoundary(StartPt_d, Ray2_d, Increase_d, MatMin, MatMax);
	if (DataLoc>0) InitialBoundaryLocs_map[DataLoc]=(unsigned char)0; // Add it to the map
	DataLoc = getANearestBoundary(StartPt_d, Ray2_d, -Increase_d, MatMin, MatMax);
	if (DataLoc>0) InitialBoundaryLocs_map[DataLoc]=(unsigned char)0; // Add it to the map

	// Ray 3
	DataLoc = getANearestBoundary(StartPt_d, Ray3_d, Increase_d, MatMin, MatMax);
	if (DataLoc>0) InitialBoundaryLocs_map[DataLoc]=(unsigned char)0; // Add it to the map
	DataLoc = getANearestBoundary(StartPt_d, Ray3_d, -Increase_d, MatMin, MatMax);
	if (DataLoc>0) InitialBoundaryLocs_map[DataLoc]=(unsigned char)0; // Add it to the map
	
	// Ray 4
	DataLoc = getANearestBoundary(StartPt_d, Ray4_d, Increase_d, MatMin, MatMax);
	if (DataLoc>0) InitialBoundaryLocs_map[DataLoc]=(unsigned char)0; // Add it to the map
	DataLoc = getANearestBoundary(StartPt_d, Ray4_d, -Increase_d, MatMin, MatMax);
	if (DataLoc>0) InitialBoundaryLocs_map[DataLoc]=(unsigned char)0; // Add it to the map
	
	
#ifdef	DEBUG_PED
	printf ("\n");
	fflush(stdout);
#endif	


	map<int, unsigned char>::iterator	Boundary_it;

#ifdef	DEBUG_PED
	int		NumRepeat=0;
#endif

	while ((int)InitialBoundaryLocs_map.size()>0 || (int)ZeroCrossingLocs_map.size()>0) {


#ifdef	DEBUG_PED
		NumRepeat++;
		if (NumRepeat%100==0) {
			printf ("InitialBoundaryLocs_map_size = %d ", (int)InitialBoundaryLocs_map.size());
			printf ("# Repeat = %d\n", NumRepeat);
			fflush(stdout);
		}
		if (NumRepeat>1000) break;
#endif


		if ((int)ZeroCrossingLocs_map.size()>0) {
			Boundary_it = ZeroCrossingLocs_map.begin();
			loc[0] = (*Boundary_it).first;


#ifdef	DEBUG_PED

			printf ("Zero Crossing Loc = ");
			for (k=0; k<(int)ZeroCrossingLocs_map.size(); k++, Boundary_it++) {
				loc[2] = (*Boundary_it).first;
				int Z = loc[2]/WtimesH_mi;
				int Y = (loc[2] - Z*WtimesH_mi)/Width_mi;
				int X = loc[2] % Width_mi;
				printf ("(%d %d %d) %X ", X, Y, Z, VoxelStatus_muc[loc[2]]);
			}
			printf ("\n");

#endif




			ZeroCrossingLocs_map.erase(loc[0]);
		}
		else {
			Boundary_it = InitialBoundaryLocs_map.begin();
			loc[0] = (*Boundary_it).first;
			InitialBoundaryLocs_map.erase(loc[0]);
		}
		DataCoor[2] = loc[0]/WtimesH_mi;
		DataCoor[1] = (loc[0] - DataCoor[2]*WtimesH_mi)/Width_mi;
		DataCoor[0] = loc[0] % Width_mi;

		// If the location is already checked, then skip it.
		// This is the same as if (VoxelStatus_muc[loc[0]]>0) continue;
		if ((VoxelStatus_muc[loc[0]] & MAT_CHECKED_VOXEL)) continue;

		// Getting the gradient vector at the position
		for (k=0; k<3; k++) GradVec[k] = (double)GradientVec_mf[loc[0]*3 + k];
		Tempd = sqrt (GradVec[0]*GradVec[0] + GradVec[1]*GradVec[1] + GradVec[2]*GradVec[2]);
		if (fabs(Tempd)<1e-6) {

#ifdef	DEBUG_PED
			printf ("Length is 0 at (%d %d %d) %X\n", DataCoor[0], DataCoor[1], DataCoor[2], VoxelStatus_muc[loc[0]]);
#endif
			VoxelStatus_muc[loc[0]] |= MAT_CHECKED_VOXEL;
			continue;
		}
		for (k=0; k<3; k++) GradVec[k] /= Tempd; // Normalize the gradient vector

		//-------------------------------------------------------------------------------------------
		// Finding the local maxima of gradient magnitudes along the gradient direction. 
		// It climbs the mountain of gradient magnitudes to find the zero-crossing second derivative 
		// Return ZeroCrossingLoc_d, LocalMaxGradient, DataPosFromZeroCrossingLoc_d
		//-------------------------------------------------------------------------------------------
		int			FoundSigma, IsCorrectBoundary;
		double		CurrDataLoc_d[3], DataValuePos_d, DataValueNeg_d;
		double		SigmaPosDir_Ret, SigmaNegDir_Ret, SecondDPosDir_Ret, SecondDNegDir_Ret;
		for (k=0; k<3; k++) CurrDataLoc_d[k] = (double)DataCoor[k];



		FoundZeroCrossingLoc = FindZeroCrossingLocation(CurrDataLoc_d, GradVec, ZeroCrossingLoc_d, 
														LocalMaxGradient, DataPosFromZeroCrossingLoc_d);
		


		for (k=0; k<3; k++) ZeroCrossingLoc_i[k] = (int)floor(ZeroCrossingLoc_d[k]);
		loc[1] = ZeroCrossingLoc_i[2]*WtimesH_mi + ZeroCrossingLoc_i[1]*Width_mi + ZeroCrossingLoc_i[0];


#ifdef	DEBUG_PED

		if (FoundZeroCrossingLoc) {
			printf ("Found Zero Crossing Loc\n");
			printf ("DataCoor = (%d %d %d), %X ", DataCoor[0], DataCoor[1], DataCoor[2], VoxelStatus_muc[loc[0]]);
			printf ("Zero2nd = %8.4f %8.4f %8.4f, ", ZeroCrossingLoc_d[0], ZeroCrossingLoc_d[1], ZeroCrossingLoc_d[2]);
			printf ("MaxGradient = %7.2f, ", LocalMaxGradient);
			printf ("DataPosFromZero2nd = %8.4f ", DataPosFromZeroCrossingLoc_d);
			printf ("At Zero = %X\n", VoxelStatus_muc[loc[1]]);
			fflush(stdout);
		}
		else {
			printf ("Did Not Find  Zero Crossing Loc\n");
			printf ("DataCoor = (%d %d %d), %X ", DataCoor[0], DataCoor[1], DataCoor[2], VoxelStatus_muc[loc[0]]);
			printf ("No Zero Crossing Loc\n");
//			FindZeroCrossingLocation(CurrDataLoc_d, GradVec, ZeroCrossingLoc_d, 
//											LocalMaxGradient, DataPosFromZeroCrossingLoc_d);
		}

#endif

		if ((VoxelStatus_muc[loc[1]] & MAT_ZERO_CROSSING)) continue;

		if (FoundZeroCrossingLoc) {

			IsCorrectBoundary = false;
			double	GradVec_d[3], SamplingInterval_d;
			GradVecInterpolation(ZeroCrossingLoc_d, GradVec_d);
			Tempd = sqrt (GradVec_d[0]*GradVec_d[0] + GradVec_d[1]*GradVec_d[1] + GradVec_d[2]*GradVec_d[2]);
			if (fabs(Tempd)<1e-6) continue;
			for (k=0; k<3; k++) GradVec_d[k] /= Tempd; // Normalize the gradient vector
			

#ifdef	DEBUG_PED
			printf ("Compute Positive Direction Sigma\n");
#endif

			FoundSigma = ComputeSigma(ZeroCrossingLoc_d, GradVec_d, SigmaPosDir_Ret, SecondDPosDir_Ret, 0.1);
			if (FoundSigma) {
				// 2.5759 = 99.00%, 3.89 = 99.99%
				for (SamplingInterval_d=1.0; SamplingInterval_d<3.89; SamplingInterval_d+=0.5) {
					DataValuePos_d = getDataValue(ZeroCrossingLoc_d, GradVec_d, SigmaPosDir_Ret*SamplingInterval_d);
					if (((double)MatMin-1e-5) <= DataValuePos_d && 
						((double)MatMax+1e-5) >= DataValuePos_d ) IsCorrectBoundary = true;
				}
			}

#ifdef	DEBUG_PED
			printf ("Compute Negative Direction Sigma\n");
#endif


			FoundSigma = ComputeSigma(ZeroCrossingLoc_d, GradVec_d, SigmaNegDir_Ret, SecondDNegDir_Ret, -0.1);
			if (FoundSigma) {
				// 2.5759 = 99.00%, 3.89 = 99.99%
				for (SamplingInterval_d=1.0; SamplingInterval_d<3.89; SamplingInterval_d+=0.5) {
					DataValueNeg_d = getDataValue(ZeroCrossingLoc_d, GradVec_d, SigmaNegDir_Ret*SamplingInterval_d);
					if (((double)MatMin-1e-5) <= DataValueNeg_d && 
						((double)MatMax+1e-5) >= DataValueNeg_d ) IsCorrectBoundary = true;
				}
			}
			

#ifdef	DEBUG_PED
			printf ("\n");
			fflush(stdout);
#endif


			if (!IsCorrectBoundary) {

#ifdef	DEBUG_PED
				printf ("Incorrect Boundary: Found Sigma = %d  ", FoundSigma);
				printf ("DataPos = %f SigmaPos = %f, ", DataValuePos_d, SigmaPosDir_Ret);
				printf ("(%f ", ZeroCrossingLoc_d[0] + GradVec_d[0]*SigmaPosDir_Ret);
				printf ("%f ",  ZeroCrossingLoc_d[1] + GradVec_d[1]*SigmaPosDir_Ret);
				printf ("%f) ", ZeroCrossingLoc_d[2] + GradVec_d[2]*SigmaPosDir_Ret);
				printf ("DataNeg = %f SigmaNeg = %f, ", DataValueNeg_d, SigmaNegDir_Ret);
				printf ("(%f ", ZeroCrossingLoc_d[0] + GradVec_d[0]*SigmaNegDir_Ret);
				printf ("%f ",  ZeroCrossingLoc_d[1] + GradVec_d[1]*SigmaNegDir_Ret);
				printf ("%f) ", ZeroCrossingLoc_d[2] + GradVec_d[2]*SigmaNegDir_Ret);
				printf ("GradVec = ");
				printf ("(%f ",    GradVec_d[0]);
				printf ("%f ",     GradVec_d[1]);
				printf ("%f)\n\n", GradVec_d[2]);
				fflush (stdout);
#endif
				VoxelStatus_muc[loc[0]] |= MAT_CHECKED_VOXEL; // Marking the Checked Voxel
				continue;
			}
			else {
				VoxelStatus_muc[loc[1]] |= MAT_ZERO_CROSSING;
				NumNeighbors = FindZeroConnectedNeighbors(loc[1], Neighbor26_i);
				for (k=0; k<NumNeighbors; k++) {
					if (VoxelStatus_muc[Neighbor26_i[k]] == 0) {
						ZeroCrossingLocs_map[Neighbor26_i[k]]=(unsigned char)0; // Add it to the map
					}
				}

#ifdef	DEBUG_PED

				if (NumNeighbors > 0) {
					printf ("NumNeighbors = %d   ", NumNeighbors);
					printf ("ZeroCrossing = (%d %d %d) ", ZeroCrossingLoc_i[0], ZeroCrossingLoc_i[1], ZeroCrossingLoc_i[2]);
					printf ("%x ", VoxelStatus_muc[loc[1]]);
					printf ("Neigbors = ");
					for (k=0; k<NumNeighbors; k++) {
						if (VoxelStatus_muc[Neighbor26_i[k]] == 0) {
							int Z = Neighbor26_i[k]/WtimesH_mi;
							int Y = (Neighbor26_i[k] - Z*WtimesH_mi)/Width_mi;
							int X = Neighbor26_i[k] % Width_mi;
							printf ("(%d %d %d) ", X, Y, Z);
						}
					}
					printf ("\n\n");
					fflush(stdout);
				}
				else {
					printf ("No Zero Crossing Neighbors\n\n");
				}

#endif


			}

		}
		VoxelStatus_muc[loc[0]] |= MAT_CHECKED_VOXEL; // Marking the Checked Voxel

		
	} // while ((int)InitialBoundaryLocs_map.size()>0)
}
*/

