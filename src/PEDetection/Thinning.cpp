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

#include <PEDetection/Thinning.h>
#include <PEDetection/ThinningTables.h>
#include <PEDetection/CompileOptions.h>

#define	INDEX(X, Y, Z) ((Z)*WtimesH_mi + (Y)*Width_mi + (X))
#define	SWAP(A, B, Temp) { (Temp)=(A); (A)=(B); (B)=(Temp); }

//----------------------------------------------------------------------------
// cThinning Class Member Functions
//----------------------------------------------------------------------------

template <class _DataType>
cThinning<_DataType>::cThinning()
{

}


// destructor
template <class _DataType>
cThinning<_DataType>::~cThinning()
{
	delete []  ThreeColorVolume_muc;
	ThreeColorVolume_muc = NULL;
	
	delete [] ConnectedComponentVolume_muc;
	ConnectedComponentVolume_muc = NULL;
}


template <class _DataType>
void cThinning<_DataType>::setData(_DataType *Data, float Minf, float Maxf)
{
	MinData_mf = Minf;
	MaxData_mf = Maxf;
	Data_mT = Data;
}


template<class _DataType>
void cThinning<_DataType>::setWHD(int W, int H, int D)
{
	Width_mi = W;
	Height_mi = H;
	Depth_mi = D;
	
	WtimesH_mi = W*H;
	WHD_mi = W*H*D;
	
	ThreeColorVolume_muc = new unsigned char[WHD_mi];
	ConnectedComponentVolume_muc = new unsigned char[WHD_mi];
}

template<class _DataType>
void cThinning<_DataType>::setProbabilityHistogram(float *Prob, int NumMaterial, int *Histo, float HistoF)
{
	MaterialProb_mf = Prob;
	Histogram_mi = Histo;
	HistogramFactorI_mf = HistoF;
	HistogramFactorG_mf = HistoF;
	NumMaterial_mi = NumMaterial;
}

template<class _DataType>
unsigned char* cThinning<_DataType>::getCCVolume()
{
	return ConnectedComponentVolume_muc;
}

template<class _DataType>
void cThinning<_DataType>::InitThreeColorVolume(_DataType MatMin, _DataType MatMax)
{
	int		i;
	

	for (i=0; i<WHD_mi; i++) {
		if (MatMin <= Data_mT[i] && Data_mT[i] <= MatMax) ThreeColorVolume_muc[i] = (unsigned char)1;
		else ThreeColorVolume_muc[i] = (unsigned char)0;
	}
	// It is for the Index() function
	ThreeColorVolume_muc[0] = (unsigned char)0;
}


template<class _DataType>
void cThinning<_DataType>::Thinning4Subfields(char *OutFileName, _DataType MatMin, _DataType MatMax)
{
	int		SK_fd, i;
	char	Sk_FileName[500];


	sprintf (Sk_FileName, "%s_Skeletons_%03d_%03d.raw", OutFileName, (int)MatMin, (int)MatMax);
	if ((SK_fd = open (Sk_FileName, O_RDONLY)) < 0) {
		cout << "Could Not Open the File, " << Sk_FileName << endl;
		cout << "Create a New Skeleton File " << Sk_FileName << endl;

		Skeletonize(OutFileName, MatMin, MatMax);

	}
	else {
		cout << "Read the Skeleton File: " << Sk_FileName << endl;
		if (read(SK_fd, ThreeColorVolume_muc, sizeof(unsigned char)*WHD_mi) != sizeof(unsigned char)*WHD_mi) {
			cout << "File Reading Error: " << Sk_FileName << endl;
			close (SK_fd);
			exit(1);
		}
		for (i=0; i<WHD_mi; i++) {
			if (ThreeColorVolume_muc[i]) ThreeColorVolume_muc[i] = (unsigned char)1;
			else ThreeColorVolume_muc[i] = (unsigned char)0;
		}
		
	#ifdef	SAVE_THREE_COLOR_VOLUME
		sprintf (Sk_FileName, "%s_SK_Read", OutFileName);
		SaveThreeColorVolume(Sk_FileName, MatMin, MatMax);
	#endif
	}
}


template<class _DataType>
void cThinning<_DataType>::ConnectedComponents(char *OutFileName, _DataType MatMin, _DataType MatMax)
{
	int				i, j, k, l, m, n, X_i, Y_i, Z_i, Num1Voxels;
	int				loc[3], NumTotal1Voxels, CubeIndex_i[27];
	map<int, unsigned char> EndVoxelStack_m, StackTemp_m, CCStack_m, CCStackBackup_m;
	map<int, unsigned char>::iterator EndVoxelStack_it, StackTemp_it, CCStack_it, CCStackBackup_it;
	

#ifdef	DEBUG_THINNING
		cout << "Removing Isolated Voxels" << endl;
		cout.flush();
#endif
	int		NumOrg1Voxels, NumEachBranch[27];

	for (i=0; i<WHD_mi; i++) ConnectedComponentVolume_muc[i] = ThreeColorVolume_muc[i];

	// Removing all isolated voxels and 6-branch and more-branch voxels, or
	// keeping only blood vessels

	EndVoxelStack_mm.clear();
	EndVoxelStack_m.clear();
	for (i=0; i<27; i++) NumEachBranch[i] = 0;
	NumOrg1Voxels = 0;
	NumTotal1Voxels = 0;
	for (k=0; k<Depth_mi; k++) {
		for (j=0; j<Height_mi; j++) {
			for (i=0; i<Width_mi; i++) {

				// Isolation Check
				loc[0] = k*WtimesH_mi + j*Width_mi + i;
				if (ConnectedComponentVolume_muc[loc[0]]==1) {
					NumOrg1Voxels++;
					Num1Voxels = 0;
					for (n=-1; n<=1; n++) {
						for (m=-1; m<=1; m++) {
							for (l=-1; l<=1; l++) {
								loc[1] = Index(i+l, j+m, k+n);
								if (ConnectedComponentVolume_muc[loc[1]]==1) {
									CubeIndex_i[Num1Voxels] = loc[1];
									Num1Voxels++;
								}
							}
						}
					}
					// Removing Isolated Voxels
					if (Num1Voxels==1) ConnectedComponentVolume_muc[loc[0]] = (unsigned char)0;
					if (Num1Voxels==2) {
						EndVoxelStack_m[loc[0]] = (unsigned char)1;
						if (loc[0]==CubeIndex_i[0])	EndVoxelStack_mm[loc[0]] = CubeIndex_i[1];
						else EndVoxelStack_mm[loc[0]] = CubeIndex_i[0];
					}
					NumEachBranch[Num1Voxels]++;
				}
				
			}
		}
	}

#ifdef	DEBUG_THINNING
	for (i=1; i<27; i++) {
		cout << "The Total Number of Each Branch [" << i << "]" << " = " << NumEachBranch[i] << endl;
	}
	cout << "The Total Number of Original One-Voxels = " << NumOrg1Voxels << endl;
	cout << "The Total Number of One-Voxels = " << NumTotal1Voxels << endl;
	cout << endl;
	cout.flush();
#endif


	// Finding all connected components from the end voxels
	do {

		EndVoxelStack_it = EndVoxelStack_m.begin();
		loc[2] = (*EndVoxelStack_it).first;
		EndVoxelStack_m.erase(loc[2]);
		if (ConnectedComponentVolume_muc[loc[2]]>=2) continue;

		CCStack_m.clear();
		CCStack_m[loc[2]] = (unsigned char)1;
		CCStackBackup_m.clear();
		CCStackBackup_m[loc[2]] = (unsigned char)1;
		do {

			CCStack_it = CCStack_m.begin();
			loc[0] = (*CCStack_it).first;
			CCStack_m.erase(loc[0]);
			Z_i = loc[0]/WtimesH_mi;
			Y_i = (loc[0] - Z_i*WtimesH_mi)/Width_mi;
			X_i = loc[0] % Width_mi;


			// ----------------------------------------------------------------------
			// Note:
			// This part should be replaced with a blood vessel recognition program
			// The next if statement is for CV (D1-SRS01) only
			//
			if (X_i<=80 || X_i>=450 || Y_i<=115 || Y_i>=360 || Z_i<=15) {
				CCStack_m.clear();
				CCStackBackup_m.clear();
				EndVoxelStack_mm.erase(loc[2]);
				break;
			}
			// ----------------------------------------------------------------------



			StackTemp_m.clear();
			for (n=-1; n<=1; n++) {
				for (m=-1; m<=1; m++) {
					for (l=-1; l<=1; l++) {
						loc[1] = Index(X_i+l, Y_i+m, Z_i+n);
						if (ConnectedComponentVolume_muc[loc[1]]==1)
							StackTemp_m[loc[1]] = (unsigned char)1;
					}
				}
			}
			
			StackTemp_m.erase(loc[0]); // Removing the center point
			if (StackTemp_m.size()==1 || StackTemp_m.size()==2 || StackTemp_m.size()==3) {
				StackTemp_it = StackTemp_m.begin();
				for (i=0; i<StackTemp_m.size(); i++, StackTemp_it++) {
					CCStackBackup_it = CCStackBackup_m.find((*StackTemp_it).first);
					if (CCStackBackup_it==CCStackBackup_m.end()) {
						CCStack_m[(*StackTemp_it).first] = (unsigned char)1;
						CCStackBackup_m[(*StackTemp_it).first] = (unsigned char)1;
					}
				}
			}

		} while (CCStack_m.size()>0);
		

#ifdef	DEBUG_THINNING
		printf ("An End Point, ");
		int	TempX, TempY, TempZ;
		TempZ = loc[2]/WtimesH_mi;
		TempY = (loc[2] - TempZ*WtimesH_mi)/Width_mi;
		TempX = loc[2] % Width_mi;
		printf ("XYZ = %3d %3d %3d: ", TempX, TempY, TempZ);
		if ((int)CCStackBackup_m.size()==0) printf ("<-- Removed (Out of Range)     ");
		else printf ("# Connected Components = %4d, ", (int)CCStackBackup_m.size());
		printf ("# End Voxels = %4d\n", (int)EndVoxelStack_m.size());
		fflush(stdout);
#endif


		// Coloring all connected voxels with 255 grey colors
		unsigned char GreyColor;
		if (CCStackBackup_m.size()>=255) GreyColor = 255;
		else if (CCStackBackup_m.size()>=2) GreyColor = CCStackBackup_m.size();
		else GreyColor = 0;
		
		CCStackBackup_it = CCStackBackup_m.begin();
		for (i=0; i<CCStackBackup_m.size(); i++, CCStackBackup_it++) {
			loc[1] = (*CCStackBackup_it).first;
			ConnectedComponentVolume_muc[loc[1]] = GreyColor;
		}
		
	} while (EndVoxelStack_m.size()>0);


#ifdef	DEBUG_THINNING
		printf ("# New End Voxels = %4d\n", (int)EndVoxelStack_mm.size());
		fflush(stdout);
#endif

#ifdef	SAVE_THREE_COLOR_VOLUME
	char	CCFileName[500];
	sprintf (CCFileName, "%s_SK_CC", OutFileName);
	SaveVolume(CCFileName, MatMin, MatMax, ConnectedComponentVolume_muc);
	SaveThickVolume(CCFileName, MatMin, MatMax, ConnectedComponentVolume_muc, (unsigned char)3);
#endif


}


template<class _DataType>
void cThinning<_DataType>::Skeletonize(char *OutFileName, _DataType MatMin, _DataType MatMax)
{
	int		i, j, NumTotalThinningRepeat_i, NumRepeat_i[9];
	int		Reduce_One_More_Time;


	InitThreeColorVolume(MatMin, MatMax);

#ifdef	SAVE_THREE_COLOR_VOLUME
	char	InitRawFileName[500];
	sprintf (InitRawFileName, "%s_Init", OutFileName);
	SaveThreeColorVolume(InitRawFileName, MatMin, MatMax);
#endif

	for (i=0; i<9; i++) NumRepeat_i[i] = 0;

	int NumSubfields;

	do {
	
	#ifdef	DEBUG_THINNING
			cout << endl;
			cout << "----- Thinning Step -----" << endl;
			cout.flush();
	#endif
		NumSubfields = 0;
		NumTotalThinningRepeat_i = 0;
		do {


	#ifdef	DEBUG_THINNING
			cout << "Reduce ORTH or DIAG Nonend Voxels of Subfield 0, 1, 2, and 3" << endl;
			cout.flush();
	#endif






			Reduce_ORTH_DIAG_Nonend_Voxels(0, 0, 0, NumRepeat_i[0]); // Subfield 0
			if (NumRepeat_i[0]>0) NumSubfields++;


			Reduce_ORTH_DIAG_Nonend_Voxels(0, 0, 1, NumRepeat_i[1]); // Subfield 1
			if (NumRepeat_i[1]>0) NumSubfields++;

			Reduce_ORTH_DIAG_Nonend_Voxels(0, 1, 0, NumRepeat_i[2]); // Subfield 2
			if (NumRepeat_i[2]>0) NumSubfields++;


			Reduce_ORTH_DIAG_Nonend_Voxels(0, 1, 1, NumRepeat_i[3]); // Subfield 3
			if (NumRepeat_i[3]>0) NumSubfields++;

	#ifdef	DEBUG_THINNING
			for (i=0; i<4; i++) {
				cout << "Num Repeat Subfield " << i << "= " << NumRepeat_i[i] << endl;
			}
			cout.flush();
	#endif


			Reduce_One_More_Time = false;
			for (i=0; i<4; i++) {
				if (NumRepeat_i[i]>0) {
					Reduce_One_More_Time = true;
					break;
				}
			}

			for (i=0; i<4; i++) {
				NumTotalThinningRepeat_i += NumRepeat_i[i];
			}

	#ifdef	DEBUG_THINNING
			cout << "---> The Total Number of Reduced Voxels of the Thinning Step = ";
			cout << NumTotalThinningRepeat_i << endl;
			cout.flush();
	#endif


		} while (Reduce_One_More_Time);


//------------------------------------------------------------------------------------
//
// Cleaning Step
//
//------------------------------------------------------------------------------------


	#ifdef	DEBUG_THINNING
			cout << "NumSubfields = " << NumSubfields << endl;
			cout << endl;
			cout << "----- Cleaning Step -----" << endl;
			cout.flush();
	#endif
		

		for (i=0; i<=(NumSubfields+1)/2; ) {



	#ifdef	DEBUG_THINNING
			cout << "i = " << i << " / " << (NumSubfields+1)/2 << endl;
			cout.flush();
	#endif

			Remove_All_Weakly_End_Voxels(0, 0, 0, NumRepeat_i[0]); // Subfield 0
			All_Upper_End_Voxels_Even	(0, 0, 0, NumRepeat_i[4]); // Subfield 0
			i++;

			if (i>(NumSubfields+1)/2) break;
			Remove_All_Weakly_End_Voxels(0, 0, 1, NumRepeat_i[1]); // Subfield 1
			All_Lower_End_Voxels_Odd	(0, 0, 1, NumRepeat_i[5]); // Subfield 1
			i++;

			if (i>(NumSubfields+1)/2) break;
			Remove_All_Weakly_End_Voxels(0, 1, 0, NumRepeat_i[2]); // Subfield 2
			All_Upper_End_Voxels_Even	(0, 1, 0, NumRepeat_i[6]); // Subfield 2
			i++;

			if (i>(NumSubfields+1)/2) break;
			Remove_All_Weakly_End_Voxels(0, 1, 1, NumRepeat_i[3]); // Subfield 3
			All_Lower_End_Voxels_Odd	(0, 1, 1, NumRepeat_i[7]); // Subfield 3
			i++;

			if (i>(NumSubfields+1)/2) break;
//			i += NumRepeat_i[3] + NumRepeat_i[7];
//			if (i>(NumTotalThinningRepeat_i+1)/2) break;
			
	#ifdef	DEBUG_THINNING
			int		SubTotal = 0;
			for (j=0; j<4; j++) {
				cout << "Num Repeat Remove All_Weakly_End_Voxels Subfield " << j << " = " << NumRepeat_i[j] << endl;
				SubTotal += NumRepeat_i[j];
			}
			for (j=0; j<4; j++) {
				if (j%2==0) cout << "Num Repeat Upper End_Voxels Subfield " << j << " = " << NumRepeat_i[j+4] << endl;
				if (j%2==1) cout << "Num Repeat Lower End_Voxels Subfield " << j << " = " << NumRepeat_i[j+4] << endl;
				SubTotal += NumRepeat_i[j+4];
			}
			cout << "Sub-Total Num Cleaned Voxels = " << SubTotal << endl;
			cout << "Total Num Cleaned Voxels = " << i << " / " << (NumSubfields+1)/2 << endl;
			cout.flush();
	#endif

			Reduce_One_More_Time = false;
			for (j=0; j<8; j++) {
				if (NumRepeat_i[j]>0) {
					Reduce_One_More_Time = true;
					break;
				}
			}
			if (!Reduce_One_More_Time) break;

		} // The End of Cleaning
		
	} while (NumTotalThinningRepeat_i > 0);
	
#ifdef	SAVE_THREE_COLOR_VOLUME
	char	SK_FileName[500];
	sprintf (SK_FileName, "%s_Skeletons", OutFileName);
	SaveThreeColorVolume(SK_FileName, MatMin, MatMax);
	cout << SK_FileName << "--> Final" << endl << endl;
	SaveThickSkeletons(SK_FileName, MatMin, MatMax);
#endif

}

template<class _DataType>
void cThinning<_DataType>::All_Lower_End_Voxels_Odd(int Px, int Py, int Pz, int& NumRepeat_Ret)
{
	int		i, j, k, l, m, n, NumRepeat, loc[3], CubeIndex_i, Num1Voxels, NeighborIndex;
//	int		StartX, StartY;
	int		Modxy, Modyz;

#ifdef	DEBUG_THINNING
	cout << "All Lower_End_Voxels_Odd: Px & Py = " << Px << " " << Py << endl;
	cout.flush();
#endif

	Modxy = (Px + Py) % 2;
	Modyz = (Py + Pz) % 2;

	NumRepeat = 0;
	for (k=0; k<Depth_mi; k++) {

//		StartY = (Py + (k % 2)) % 2;
//		StartX = (Px + (k % 2)) % 2;
//		for (j=StartY; j<Height_mi; j+=2) {
//			for (i=StartX; i<Width_mi; i+=2) {

		for (j=0; j<Height_mi; j++) {
			for (i=0; i<Width_mi; i++) {

				if ((i+j)%2==Modxy && (j+k)%2==Modyz) { }
				else continue;

				loc[0] = k*WtimesH_mi + j*Width_mi + i;
				if (ThreeColorVolume_muc[loc[0]]==0) continue;

				Num1Voxels = 0;
				CubeIndex_i = -1;
				for (n=-1; n<=1; n++) {
					for (m=-1; m<=1; m++) {
						for (l=-1; l<=1; l++) {
							CubeIndex_i++;
							if (ThreeColorVolume_muc[Index(i+l, j+m, k+n)]==1) {
								if (l==0 && m==0 && n==0) continue;
								NeighborIndex = CubeIndex_i;
								Num1Voxels++;
							}
						}
					}
				}

				if (Num1Voxels==1 && NeighborIndex!=6 && NeighborIndex!=8 && 
					NeighborIndex!=24 && NeighborIndex!=26) {
					ThreeColorVolume_muc[loc[0]] = (unsigned char)0;
					NumRepeat++;
				}


			}
		}
	}

	NumRepeat_Ret = NumRepeat;


}

template<class _DataType>
void cThinning<_DataType>::All_Upper_End_Voxels_Even(int Px, int Py, int Pz, int& NumRepeat_Ret)
{
	int		i, j, k, l, m, n, NumRepeat, loc[3], CubeIndex_i, Num1Voxels, NeighborIndex;
//	int		StartX, StartY;
	int		Modxy, Modyz;

#ifdef	DEBUG_THINNING
	cout << "All Upper_End_Voxels_Even: (Px, Py, Pz) = " << Px << " " << Py << " " << Pz << endl;
	cout.flush();
#endif

	Modxy = (Px + Py) % 2;
	Modyz = (Py + Pz) % 2;

	NumRepeat = 0;
	for (k=0; k<Depth_mi; k++) {
		for (j=0; j<Height_mi; j++) {
			for (i=0; i<Width_mi; i++) {

				if ((i+j)%2==Modxy && (j+k)%2==Modyz) { }
				else continue;


				loc[0] = k*WtimesH_mi + j*Width_mi + i;
				if (ThreeColorVolume_muc[loc[0]]==0) continue;

				Num1Voxels = 0;
				CubeIndex_i = -1;
				for (n=-1; n<=1; n++) {
					for (m=-1; m<=1; m++) {
						for (l=-1; l<=1; l++) {
							CubeIndex_i++;
							if (ThreeColorVolume_muc[Index(i+l, j+m, k+n)]==1) {
								if (l==0 && m==0 && n==0) continue;
								NeighborIndex = CubeIndex_i;
								Num1Voxels++;
							}
						}
					}
				}

				if (Num1Voxels==1 && NeighborIndex!=0 && NeighborIndex!=2 && 
					NeighborIndex!=18 && NeighborIndex!=20) {
					ThreeColorVolume_muc[loc[0]] = (unsigned char)0;
					NumRepeat++;
				}

			}
		}
	}

	NumRepeat_Ret = NumRepeat;


}

template<class _DataType>
void cThinning<_DataType>::Remove_All_Weakly_End_Voxels(int Px, int Py, int Pz, int& NumRepeat_Ret)
{
	int		i, j, k, NumRepeat, loc[3];
	int		Modxy, Modyz;
	

#ifdef	DEBUG_THINNING
	cout << "Remove All_Weakly_End_Voxels: (Px, Py, Pz) = " << Px << " " << Py << " " << Pz << endl;
	cout.flush();
#endif


	Modxy = (Px + Py) % 2;
	Modyz = (Py + Pz) % 2;

	NumRepeat = 0;
	for (k=0; k<Depth_mi; k++) {
		for (j=0; j<Height_mi; j++) {
			for (i=0; i<Width_mi; i++) {

				if ((i+j)%2==Modxy && (j+k)%2==Modyz) { }
				else continue;
				

				loc[0] = k*WtimesH_mi + j*Width_mi + i;
				if (ThreeColorVolume_muc[loc[0]]==0) continue;

/*
				if (Is_Adjacent_Two_1_Voxels(i, j, k)) {
					cout << "Is Adjacent_Two_1_Voxels = true, (i, j, k) = ";
					cout << i << " " << j << " " << k << endl;
				}
				if (Is_Adjacent_Three_1_Voxels(i, j, k)) {
					cout << "Is Adjacent_Three_1_Voxels = true, (i, j, k) = ";
					cout << i << " " << j << " " << k << endl;
				}
*/
				
				
				
				if (Is_Adjacent_Two_1_Voxels(i, j, k) || Is_Adjacent_Three_1_Voxels(i, j, k)) {
					ThreeColorVolume_muc[loc[0]] = (unsigned char)0;
					NumRepeat++;
				}

			}
		}
	}

	NumRepeat_Ret = NumRepeat;


}

// p is adjacent to exactly two 1-voxels that are 26-adjacent to each other
template<class _DataType>
int cThinning<_DataType>::Is_Adjacent_Two_1_Voxels(int Px, int Py, int Pz)
{
	int				i, j, k, Num1Voxels, CubeIndex_i, OneVoxelIndex[27], FirstSecond;
	unsigned char	VoxelValue;


	Num1Voxels = 0;
	CubeIndex_i = -1;
	FirstSecond = 0;
	for (k=-1; k<=1; k++) {
		for (j=-1; j<=1; j++) {
			for (i=-1; i<=1; i++) {

				CubeIndex_i++;
				VoxelValue = ThreeColorVolume_muc[Index(Px+i, Py+j, Pz+k)];
				if (VoxelValue==1) {
					if (i==0 && j==0 && k==0) continue; // Skip the center voxel
					else {
						OneVoxelIndex[FirstSecond++] = CubeIndex_i;
						Num1Voxels++;
					}
				}
			}
		}
	}

	if (Num1Voxels==2) return Two_OneVoxels_26Adjacent[OneVoxelIndex[0]][OneVoxelIndex[1]];
	else return false;
}

// p is adjacent to exactly three 1-voxels that are in the same boundary layer of N(p)
template<class _DataType>
int cThinning<_DataType>::Is_Adjacent_Three_1_Voxels(int Px, int Py, int Pz)
{
	int				i, j, k, CubeIndex_i, OneVoxelIndex[27], FirstSecondThird, Num1Voxels;
	unsigned char	VoxelValue;


	CubeIndex_i = -1;
	FirstSecondThird = 0;
	Num1Voxels = 0;
	for (k=-1; k<=1; k++) {
		for (j=-1; j<=1; j++) {
			for (i=-1; i<=1; i++) {

				CubeIndex_i++;
				VoxelValue = ThreeColorVolume_muc[Index(Px+i, Py+j, Pz+k)];
				if (VoxelValue==1) {
					if (i==0 && j==0 && k==0) continue; // Skip the center voxel
					else {
						OneVoxelIndex[FirstSecondThird++] = CubeIndex_i;
						Num1Voxels++;
					}
				}
			}
		}
	}


	if (Num1Voxels==3 && (Three_OneVoxels_SameBoundary[OneVoxelIndex[0]] &
						  Three_OneVoxels_SameBoundary[OneVoxelIndex[1]] &
						  Three_OneVoxels_SameBoundary[OneVoxelIndex[2]]) > 0 ) return true;
	else return false;

	return false;
}


template<class _DataType>
void cThinning<_DataType>::Reduce_ORTH_DIAG_Nonend_Voxels(int Px, int Py, int Pz, int& NumRepeat_Ret)
{
	int				i, j, k, l, m, n, NumRepeat, loc[3], CubeIndex_i, Num1Voxels;
	unsigned char	Cube_uc1[27], Cube_uc2[27];
	int				Modxy, Modyz;

	Modxy = (Px + Py) % 2;
	Modyz = (Py + Pz) % 2;


	NumRepeat = 0;
	for (k=0; k<Depth_mi; k++) {
		for (j=0; j<Height_mi; j++) {
			for (i=0; i<Width_mi; i++) {

				if ((i+j)%2==Modxy && (j+k)%2==Modyz) { }
				else continue;

/*
				if (
					(i==10 && j==20 && k==32) || 
					(i==12 && j==20 && k==32) || 
					(i==10 && j==22 && k==32) ||
					(k==1000)
					

				) {
					int 	SubfieldNum = Modxy+Modyz+1;
					if (Modxy==1 && Modyz==0) SubfieldNum = 4;
					printf ("Subfield # = %d, (%d %d %d), ", SubfieldNum, i, j, k);
					printf ("Modxy = %d, Modyz = %d, ", Modxy, Modyz);
					printf ("Px, Py, Pz = %d %d %d ", Px, Py, Pz);
					
					loc[0] = k*WtimesH_mi + j*Width_mi + i;
					printf ("Color = %d", ThreeColorVolume_muc[loc[0]]);
					printf ("\n");
					
					fflush(stdout);
				}
*/







				loc[0] = k*WtimesH_mi + j*Width_mi + i;
				if (ThreeColorVolume_muc[loc[0]]==0) continue;

				// Extracing All 26 Neighbors
				CubeIndex_i = -1;
				Num1Voxels = 0;
				for (n=-1; n<=1; n++) {
					for (m=-1; m<=1; m++) {
						for (l=-1; l<=1; l++) {
							CubeIndex_i++;
							if (ThreeColorVolume_muc[Index(i+l, j+m, k+n)]==1) {
								Cube_uc1[CubeIndex_i] = (unsigned char)1; // for ORTH_Reducible
								Cube_uc2[CubeIndex_i] = (unsigned char)1; // A copy for DIAG_Reducible
								Num1Voxels++;
							}
							else {
								Cube_uc1[CubeIndex_i] = (unsigned char)0; // for ORTH_Reducible
								Cube_uc2[CubeIndex_i] = (unsigned char)0; // A copy for DIAG_Reducible
							}
						}
					}
				}




				if (
					(i==10 && j==20 && k==32) || 
					(i==12 && j==20 && k==32) || 
					(i==10 && j==22 && k==32) ||
					(k==10000) 
					

				) {
					printf ("ORTH_DIAG: Num1Voxels = %d, XYZ=(%d %d %d)\n", Num1Voxels, i, j, k);
					fflush(stdout);
				}






				if (Num1Voxels==2) continue; // Do not reduce the end voxel, since it is an end voxel
/*
				if (Num1Voxels==3) {
					if (Is_DiametricallyAdjacent(Cube_uc1)) continue;
//						if (Is_DiametricallyAdjacent(i, j, k, Cube_uc1)) continue;
				}
*/

/*
				if (i==14 && j==25 && k==33) {
					printf ("14 25 33!!!\n");
					exit(1);
				
				}
				printf ("Before ORTH DIAG: i j k = %d %d %d, ", i, j, k);
				printf ("Px & Py = (%d %d), ", Px, Py);
				printf ("\n");
				fflush(stdout);
*/

				if (Is_ORTH_Reducible(i, j, k, Cube_uc1) || Is_DIAG_Reducible(i, j, k, Cube_uc2)) {

/*
					printf ("Px & Py = (%d %d), ", Px, Py);
					printf ("Reducible = (%d %d %d) ------------------ ", i, j, k);
					printf ("\n");
*/

					ThreeColorVolume_muc[loc[0]] = (unsigned char)0;
					NumRepeat++;
				}
				else {
				
/*
					printf ("Px & Py = (%d %d), ", Px, Py);
					printf ("Not Reducible = (%d %d %d)\n", i, j, k);
*/
				
				}


			}
		}
	}

	NumRepeat_Ret = NumRepeat;
}

//
//  X  Y  Z            X  Y  Z
// +1 -1 -1  lne-usw  -1 +1 +1
// -1 -1 -1  lnw-use  +1 +1 +1
// +1 -1 +1  lse-unw  -1 +1 -1
// -1 -1 +1  lsw-une  +1 +1 -1
template<class _DataType>
int	cThinning<_DataType>::Is_DiametricallyAdjacent(int Xi, int Yi, int Zi, unsigned char *Cube_uc)
{
	int		l, m, n, Num1Voxels[2];
	int		lne_i, usw_i, lnw_i, use_i, lse_i, unw_i, lsw_i, une_i;
	
	
	lne_i = 20;	lnw_i = 18;	lsw_i =  0;	lse_i =  2;
	une_i = 26;	unw_i = 24;	usw_i =  6;	use_i =  8;
	
//  X  Y  Z            X  Y  Z
// +1 -1 -1  lne-usw  -1 +1 +1
	if (Cube_uc[lne_i]==1 && Cube_uc[usw_i]==1) {
		Num1Voxels[0] = 0;
		Num1Voxels[1] = 0;
		for (n=-1; n<=1; n++) {
			for (m=-1; m<=1; m++) {
				for (l=-1; l<=1; l++) {
					if (ThreeColorVolume_muc[Index(Xi+l+2, Yi+m-2, Zi+n-2)]==1) Num1Voxels[0]++;
				}
			}
		}
		for (n=-1; n<=1; n++) {
			for (m=-1; m<=1; m++) {
				for (l=-1; l<=1; l++) {
					if (ThreeColorVolume_muc[Index(Xi+l-2, Yi+m+2, Zi+n+2)]==1) Num1Voxels[1]++;
				}
			}
		}
		if (Num1Voxels[0]==2 && Num1Voxels[1]==2) return true;
	}
//  X  Y  Z            X  Y  Z
// -1 -1 -1  lnw-use  +1 +1 +1
	if (Cube_uc[lnw_i]==1 && Cube_uc[use_i]==1) {
		Num1Voxels[0] = 0;
		Num1Voxels[1] = 0;
		for (n=-1; n<=1; n++) {
			for (m=-1; m<=1; m++) {
				for (l=-1; l<=1; l++) {
					if (ThreeColorVolume_muc[Index(Xi+l-2, Yi+m-2, Zi+n-2)]==1) Num1Voxels[0]++;
				}
			}
		}
		for (n=-1; n<=1; n++) {
			for (m=-1; m<=1; m++) {
				for (l=-1; l<=1; l++) {
					if (ThreeColorVolume_muc[Index(Xi+l+2, Yi+m+2, Zi+n+2)]==1) Num1Voxels[1]++;
				}
			}
		}
		if (Num1Voxels[0]==2 && Num1Voxels[1]==2) return true;
	}
//  X  Y  Z            X  Y  Z
// +1 -1 +1  lse-unw  -1 +1 -1
	if (Cube_uc[lse_i]==1 && Cube_uc[unw_i]==1) {
		Num1Voxels[0] = 0;
		Num1Voxels[1] = 0;
		for (n=-1; n<=1; n++) {
			for (m=-1; m<=1; m++) {
				for (l=-1; l<=1; l++) {
					if (ThreeColorVolume_muc[Index(Xi+l+2, Yi+m-2, Zi+n+2)]==1) Num1Voxels[0]++;
				}
			}
		}
		for (n=-1; n<=1; n++) {
			for (m=-1; m<=1; m++) {
				for (l=-1; l<=1; l++) {
					if (ThreeColorVolume_muc[Index(Xi+l-2, Yi+m+2, Zi+n-2)]==1) Num1Voxels[1]++;
				}
			}
		}
		if (Num1Voxels[0]==2 && Num1Voxels[1]==2) return true;
	}
//  X  Y  Z            X  Y  Z
// -1 -1 +1  lsw-une  +1 +1 -1
	if (Cube_uc[lsw_i]==1 && Cube_uc[une_i]==1) {
		Num1Voxels[0] = 0;
		Num1Voxels[1] = 0;
		for (n=-1; n<=1; n++) {
			for (m=-1; m<=1; m++) {
				for (l=-1; l<=1; l++) {
					if (ThreeColorVolume_muc[Index(Xi+l-2, Yi+m-2, Zi+n+2)]==1) Num1Voxels[0]++;
				}
			}
		}
		for (n=-1; n<=1; n++) {
			for (m=-1; m<=1; m++) {
				for (l=-1; l<=1; l++) {
					if (ThreeColorVolume_muc[Index(Xi+l+2, Yi+m+2, Zi+n-2)]==1) Num1Voxels[1]++;
				}
			}
		}
		if (Num1Voxels[0]==2 && Num1Voxels[1]==2) return true;
	}

	return false;	
}


template<class _DataType>
int	cThinning<_DataType>::Is_DiametricallyAdjacent(unsigned char *Cube_uc)
{
	int		lne_i, usw_i, lnw_i, use_i, lse_i, unw_i, lsw_i, une_i;
	
	
	lne_i = 20;	lnw_i = 18;	lsw_i =  0;	lse_i =  2;
	une_i = 26;	unw_i = 24;	usw_i =  6;	use_i =  8;
	
	if (Cube_uc[lne_i]==1 && Cube_uc[usw_i]==1) return true;
	if (Cube_uc[lnw_i]==1 && Cube_uc[use_i]==1) return true;
	if (Cube_uc[lse_i]==1 && Cube_uc[unw_i]==1) return true;
	if (Cube_uc[lsw_i]==1 && Cube_uc[une_i]==1) return true;

	return false;	
}


/*
	ai[0] = 17; bi[0] = 14;		ci[0] =  8; di[0] =  5; ei[0] =  4; 
	ai[1] =  7; bi[1] =  4; 	ci[1] =  6; di[1] = 12; ei[1] =  3; 
	ai[2] = 15; bi[2] = 12; 	ci[2] = 24; di[2] = 22; ei[2] = 21; 
	ai[3] = 25; bi[3] = 22; 	ci[3] = 26; di[3] = 14; ei[3] = 23; 
*/

int ai_g[4] = {17, 7, 15, 25};
int bi_g[4] = {14, 4, 12, 22};
int ci_g[4] = {8, 6, 24, 26};
int di_g[4] = {5, 12, 22, 14};
int ei_g[4] = {4, 3, 21, 23};


template<class _DataType>
int	cThinning<_DataType>::Is_ORTH_Reducible(int Xi, int Yi, int Zi, unsigned char *Cube_uc)
{
	int				RotationNum_i;
	int				upi, lpi;


	upi = 16; lpi = 10;
	RotationNum_i = -1;
	while (RotationNum_i<=5) {
	
		RotateCube_ORTH(RotationNum_i++, Cube_uc);

/*
		if (Xi==15 && Yi==7 && Zi==33) {
			printf ("Rotation Num = %d, (%d %d %d)\n", RotationNum_i, Xi, Yi, Zi);
			DisplayCube(Cube_uc);
		}
*/

		if ((Cube_uc[upi]==0 && Cube_uc[lpi]==1) &&
			(!(Cube_uc[ai_g[0]]==1) || (Cube_uc[bi_g[0]]==1)) &&
			(!(Cube_uc[ci_g[0]]==1) || (Cube_uc[bi_g[0]]==1 || Cube_uc[di_g[0]]==1 || Cube_uc[ei_g[0]]==1)) &&
			(!(Cube_uc[ai_g[1]]==1) || (Cube_uc[bi_g[1]]==1)) &&
			(!(Cube_uc[ci_g[1]]==1) || (Cube_uc[bi_g[1]]==1 || Cube_uc[di_g[1]]==1 || Cube_uc[ei_g[1]]==1)) &&
			(!(Cube_uc[ai_g[2]]==1) || (Cube_uc[bi_g[2]]==1)) &&
			(!(Cube_uc[ci_g[2]]==1) || (Cube_uc[bi_g[2]]==1 || Cube_uc[di_g[2]]==1 || Cube_uc[ei_g[2]]==1)) &&
			(!(Cube_uc[ai_g[3]]==1) || (Cube_uc[bi_g[3]]==1)) &&
			(!(Cube_uc[ci_g[3]]==1) || (Cube_uc[bi_g[3]]==1 || Cube_uc[di_g[3]]==1 || Cube_uc[ei_g[3]]==1))

			) {
			
				return true;
			}
	
	}
	
	return false;
}

template<class _DataType>
void cThinning<_DataType>::DisplayCube(unsigned char *Cube27)
{
	int		i;
	
	for (i=6; i<=8; i++) printf ("%1d ", Cube27[i]); printf ("\n");
	for (i=3; i<=5; i++) printf ("%1d ", Cube27[i]); printf ("\n");
	for (i=0; i<=2; i++) printf ("%1d ", Cube27[i]); printf ("\n");
	printf ("\n");

	for (i=15; i<=17; i++) printf ("%1d ", Cube27[i]); printf ("\n");
	for (i=12; i<=14; i++) printf ("%1d ", Cube27[i]); printf ("\n");
	for (i=9; i<=11; i++) printf ("%1d ", Cube27[i]); printf ("\n");
	printf ("\n");
	
	for (i=24; i<=26; i++) printf ("%1d ", Cube27[i]); printf ("\n");
	for (i=21; i<=23; i++) printf ("%1d ", Cube27[i]); printf ("\n");
	for (i=18; i<=20; i++) printf ("%1d ", Cube27[i]); printf ("\n");
	printf ("\n");

/*
	for (int i=0; i<27; i++) {
		printf ("%1d ", Cube27[i]);
		if ((i%3)==2) printf ("\n");
		if ((i%9)==8) printf ("\n");
	}
*/
}

template<class _DataType>
int	cThinning<_DataType>::Is_ORTH_Edge(int Xi, int Yi, int Zi, unsigned char *Cube_uc)
{
	int				RotationNum_i, New[6];
	int				ai[8], bi[8], ci[8], upi, lpi;
	
	upi = 16; lpi = 10;
	New[0] =  1;	New[1] =  4;	New[2] =  7;	New[3] = 25;	New[4] = 22;	New[5] = 29;
	ai[1] =  8;	ai[2] =  5;	ai[3] =  2;
	bi[1] = 17;	bi[2] = 14;	bi[3] =  11;
	ci[1] = 26;	ci[2] = 23;	ci[3] =  20;

	ai[4] =  8;	ai[5] =  5;	ai[6] =  2;
	bi[4] = 17;	bi[5] = 14;	bi[6] =  11;
	ci[4] = 26;	ci[5] = 23;	ci[6] =  20;


	RotationNum_i = 0;
	do {
	
		if ((Cube_uc[upi]==0 && Cube_uc[lpi]==1) &&
			(Cube_uc[New[0]]==0 && Cube_uc[New[1]]==0 && Cube_uc[New[2]]==0 && 
			 Cube_uc[New[3]]==0 && Cube_uc[New[4]]==0 && Cube_uc[New[5]]==0) &&

			(  (Cube_uc[ai[1]]==0 && Cube_uc[ai[2]]==0 && Cube_uc[ai[3]]==0 && 
				Cube_uc[ci[1]]==0 && Cube_uc[ci[2]]==0 && Cube_uc[ci[3]]==0 && 
				(!(Cube_uc[bi[1]]==1) || (Cube_uc[bi[2]]==1))) ||
			   ((Cube_uc[bi[1]]==0 && Cube_uc[bi[2]]==0 && Cube_uc[bi[3]]==0) && 
			     ((!(Cube_uc[ai[1]]==1) || (Cube_uc[ai[2]]==1)) || (!(Cube_uc[ci[1]]==1) || (Cube_uc[ci[2]]==1))) ) ||
			   (Cube_uc[ai[1]]==0 && Cube_uc[ai[2]]==0 && Cube_uc[ai[3]]==0 && 
				Cube_uc[bi[1]]==0 && Cube_uc[bi[2]]==0 && Cube_uc[bi[3]]==0 && 
				Cube_uc[ci[1]]==0 && Cube_uc[ci[2]]==0 && Cube_uc[ci[3]]==0) ||

			   (Cube_uc[ai[4]]==0 && Cube_uc[ai[5]]==0 && Cube_uc[ai[6]]==0 && 
				Cube_uc[ci[4]]==0 && Cube_uc[ci[5]]==0 && Cube_uc[ci[6]]==0 && 
				(!(Cube_uc[bi[4]]==1) || (Cube_uc[bi[5]]==1))) ||
			   ((Cube_uc[bi[4]]==0 && Cube_uc[bi[5]]==0 && Cube_uc[bi[6]]==0) && 
			     ((!(Cube_uc[ai[4]]==1) || (Cube_uc[ai[5]]==1)) || (!(Cube_uc[ci[4]]==1) || (Cube_uc[ci[5]]==1))) ) ||
			   (Cube_uc[ai[4]]==0 && Cube_uc[ai[5]]==0 && Cube_uc[ai[6]]==0 && 
				Cube_uc[bi[4]]==0 && Cube_uc[bi[5]]==0 && Cube_uc[bi[6]]==0 && 
				Cube_uc[ci[4]]==0 && Cube_uc[ci[5]]==0 && Cube_uc[ci[6]]==0)

			)
			
			) return true;

		RotateCube_ORTH(RotationNum_i++, Cube_uc);
		
	} while (RotationNum_i<=5);
	
	return false;
}


/*
	ai[0] = 4;	ai[1] = 22;
	bi[0] = 7;	bi[1] = 25;
	ci[0] = 5;	ci[1] = 23;
	di[0] = 3;	di[1] = 19;
	ei[0] = 1;	ei[1] = 21;
	Nsn3x3[0] = 10;	Nsn3x3[1] = 11;	Nsn3x3[2] = 12;	Nsn3x3[3] = 14;
	Nsn3x3[4] = 15;	Nsn3x3[5] = 16;	Nsn3x3[6] = 17;
*/

int Dai_g[2] = {4, 22};
int Dbi_g[2] = {7, 25};
int Dci_g[2] = {5, 23};
int Ddi_g[2] = {3, 19};
int Dei_g[2] = {1, 21};
int Nsn3x3_gi[7] = {10, 11, 12, 14, 15, 16, 17};


template<class _DataType>
int	cThinning<_DataType>::Is_DIAG_Reducible(int Xi, int Yi, int Zi, unsigned char *Cube_uc)
{
	int		RotationNum_i;
	int		lwpi;


	//----------------------------------------------------------
	// Checking whether 3x3 neighborhood containing lw(p) are 0s
	int	EdgeI1[8] = {1, 10, 19, 22, 25, 16, 7, 4};
	int	EdgeI2[8] = {3, 4, 5, 14, 23, 22, 21, 12};
	int	EdgeI3[8] = {9, 10, 11, 14, 17, 16, 15, 12};
	int	i, Print_b, NumZeros;
	
	Print_b = false;
	NumZeros = 0;
	for (i=0; i<8; i++) {
		if (Cube_uc[EdgeI1[i]]==0) NumZeros++;
	}
	if (NumZeros==7 && (Cube_uc[1]==0 || Cube_uc[19]==0 || Cube_uc[25]==0 || Cube_uc[7]==0)) 
		Print_b = true;
	
	if (Print_b==false) {
		NumZeros = 0;
		for (i=0; i<8; i++) {
			if (Cube_uc[EdgeI2[i]]==0) NumZeros++;
		}
		if (NumZeros==7 && (Cube_uc[3]==0 || Cube_uc[5]==0 || Cube_uc[23]==0 || Cube_uc[21]==0)) 
			Print_b = true;
	}
	
	if (Print_b==false) {
		NumZeros = 0;
		for (i=0; i<8; i++) {
			if (Cube_uc[EdgeI3[i]]==0) NumZeros++;
		}
		if (NumZeros==7 && (Cube_uc[9]==0 || Cube_uc[11]==0 || Cube_uc[17]==0 || Cube_uc[15]==0)) 
			Print_b = true;
	}
	
	if (Xi==14 && Yi==25) Print_b = true;
	
	//----------------------------------------------------------


	lwpi = 9;
	RotationNum_i = -1;
	while (RotationNum_i<=11) {
	
		RotateCube_DIAG(RotationNum_i++, Cube_uc);
		
		
		
/*		
		//--------------------------------------------------------------------------------
		if (Print_b) {
			printf ("X Y Z = (%d %d %d), Rotation Num = %d\n", Xi, Yi, Zi, RotationNum_i);
			DisplayCube(Cube_uc);
			fflush(stdout);
		}
		//--------------------------------------------------------------------------------
*/
		
		
		
	
		if ((Cube_uc[lwpi]==1 && Cube_uc[Nsn3x3_gi[0]]==0 && Cube_uc[Nsn3x3_gi[1]]==0 && 
				Cube_uc[Nsn3x3_gi[2]]==0 && Cube_uc[Nsn3x3_gi[3]]==0 && Cube_uc[Nsn3x3_gi[4]]==0 && 
				Cube_uc[Nsn3x3_gi[5]]==0 && Cube_uc[Nsn3x3_gi[6]]==0) &&

			 (!(Cube_uc[Dai_g[0]]==0) || (Cube_uc[Dbi_g[0]]==0 || Cube_uc[Dci_g[0]]==0 || 
			 	Cube_uc[Ddi_g[0]]==0 || Cube_uc[Dei_g[0]]==0)) &&
			 (!(Cube_uc[Dai_g[1]]==0) || (Cube_uc[Dbi_g[1]]==0 || Cube_uc[Dci_g[1]]==0 || 
			 	Cube_uc[Ddi_g[1]]==0 || Cube_uc[Dei_g[1]]==0)) &&

			 ConnectedTo_lwp(Cube_uc) ) {

			 
/*
				//--------------------------------------------------------------------------------
			 	printf ("DIAG: Reducible Voxle! --------------------------------------------------------- ");
				printf ("\n\n\n");
			 	fflush(stdout);
	 			//--------------------------------------------------------------------------------
*/
			 
			 	return true;
			}
	}

	return false;
}

// 26-neighbor connect
template<class _DataType>
int cThinning<_DataType>::ConnectedTo_lwp2(unsigned char *Cube27)
{
	int				i, NumOneVoxels_south, NumOneVoxels_north;
	unsigned char	ConnectingIndex, Connected;

	
	NumOneVoxels_south = 0;
	NumOneVoxels_north = 0;
	for (i=0; i<=8; i++) if (Cube27[i]==1) NumOneVoxels_south++;
	for (i=18; i<=26; i++) if (Cube27[i]==1) NumOneVoxels_north++;
	if ((NumOneVoxels_south + NumOneVoxels_north)==0) return true;

	Connected = true;
	if (NumOneVoxels_south>0) {
		if (Cube27[4]==1) Connected = true;
		if (Cube27[4]==0) {
			ConnectingIndex = Cube27[0] | Cube27[1]<<1 | Cube27[2]<<2 | Cube27[3]<<3 | 
							Cube27[5]<<4 | Cube27[6]<<5 | Cube27[7]<<6 | Cube27[8]<<7;
			Connected = ConnectedTo_lwpTable[ConnectingIndex];				
			if (Connected==0) return false;
		}
	}

	if (NumOneVoxels_north>0) {
		if (Cube27[22]==1) Connected = true;
		if (Cube27[22]==0) {
			ConnectingIndex = Cube27[18] | Cube27[19]<<1 | Cube27[20]<<2 | Cube27[21]<<3 | 
							Cube27[23]<<4 | Cube27[24]<<5 | Cube27[25]<<6 | Cube27[26]<<7;
			Connected = ConnectedTo_lwpTable[ConnectingIndex];				
			if (Connected==0) return false;
		}
	}
	
	return true;
}

// 18-neighbor connect
template<class _DataType>
int cThinning<_DataType>::ConnectedTo_lwp(unsigned char *Cube27)
{
	int				i, NumOneVoxels_south, NumOneVoxels_north;
	unsigned char	ConnectingIndex, Connected;

	
	NumOneVoxels_south = 0;
	NumOneVoxels_north = 0;
	for (i=0; i<9; i++) if (Cube27[i]==1) NumOneVoxels_south++;
	for (i=18; i<27; i++) if (Cube27[i]==1) NumOneVoxels_north++;
	if ((NumOneVoxels_south + NumOneVoxels_north)==0) return true;

	if (NumOneVoxels_south>0) {
		if (Cube27[4]==1) {
			if (Cube27[0]==1 || Cube27[1]==1 || Cube27[3]==1) Connected = true;
			else Connected = false;
		}
		else { // if (Cube27[4]==0)
			ConnectingIndex = Cube27[0] | Cube27[1]<<1 | Cube27[2]<<2 | Cube27[3]<<3 | 
							Cube27[5]<<4 | Cube27[6]<<5 | Cube27[7]<<6 | Cube27[8]<<7;
			Connected = ConnectedTo_lwpTable[ConnectingIndex];				
		}
		if (!Connected) return false;
	}

	if (NumOneVoxels_north>0) {
		if (Cube27[22]==1) {
			if (Cube27[18]==1 || Cube27[19]==1 || Cube27[21]==1) Connected = true;
			else Connected = false;
		}
		else { // if (Cube27[22]==0)
			ConnectingIndex = Cube27[18] | Cube27[19]<<1 | Cube27[20]<<2 | Cube27[21]<<3 | 
							Cube27[23]<<4 | Cube27[24]<<5 | Cube27[25]<<6 | Cube27[26]<<7;
			Connected = ConnectedTo_lwpTable[ConnectingIndex];				
		}
		if (!Connected) return false;
	}
	
	return true;
}


template<class _DataType>
void cThinning<_DataType>::RotateCube_ORTH(int RotationNumber, unsigned char *Cube27)
{
	switch (RotationNumber) {
		case 0: break;
		case 1: case 2: case 3: RotateAroundZ(Cube27); break;

		case 4: RotateAroundY(Cube27); break;
		case 5: RotateAroundX(Cube27); RotateAroundX(Cube27); break;
		default : return;
	}
}

template<class _DataType>
void cThinning<_DataType>::RotateCube_DIAG(int RotationNumber, unsigned char *Cube27)
{
	switch (RotationNumber) {
		case 0: break;
		case 1: case 2: case 3: RotateAroundY(Cube27); break;

		case 4: RotateAroundX(Cube27); break;
		case 5: case 6: case 7: RotateAroundY(Cube27); break;
		
		case 8: RotateAroundX(Cube27); break;
		case 9: case 10: case 11: RotateAroundY(Cube27); break;
		
		default : return;
	}
}


int	RotationAroundXTable[3][8] = {
	{2, 5, 8, 17, 26, 23, 20, 11},
	{1, 4, 7, 16, 25, 22, 19, 10},
	{0, 3, 6, 15, 24, 21, 18,  9}
};

int	RotationAroundYTable[3][8] = {
	{6, 7, 8, 17, 26, 25, 24, 15},
	{3, 4, 5, 14, 23, 22, 21, 12},
	{0, 1, 2, 11, 20, 19, 18,  9}
};

int	RotationAroundZTable[3][8] = {
	{ 0,  1,  2,  5,  8,  7,  6,  3},
	{ 9, 10, 11, 14, 17, 16, 15, 12},
	{18, 19, 20, 23, 26, 25, 24, 21},
};


template<class _DataType>
void cThinning<_DataType>::RotateAroundX(unsigned char *Cube27)
{
	int				i, j;
	unsigned char	Edge_uc[3][2];
	
	for (j=0; j<3; j++) {
		for (i=0; i<=1; i++) {
			Edge_uc[j][i] = Cube27[RotationAroundXTable[j][i]];
		}
	}
	for (j=0; j<3; j++) {
		for (i=0; i<=5; i++) {
			Cube27[RotationAroundXTable[j][i]] = Cube27[RotationAroundXTable[j][i+2]];
		}
	}
	for (j=0; j<3; j++) {
		for (i=6; i<=7; i++) {
			Cube27[RotationAroundXTable[j][i]] = Edge_uc[j][i-6];
		}
	}
}

template<class _DataType>
void cThinning<_DataType>::RotateAroundY(unsigned char *Cube27)
{
	int				i, j;
	unsigned char	Edge_uc[3][2];
	
	for (j=0; j<3; j++) {
		for (i=0; i<2; i++) {
			Edge_uc[j][i] = Cube27[RotationAroundYTable[j][i]];
		}
	}
	for (j=0; j<3; j++) {
		for (i=0; i<=5; i++) {
			Cube27[RotationAroundYTable[j][i]] = Cube27[RotationAroundYTable[j][i+2]];
		}
	}
	for (j=0; j<3; j++) {
		for (i=6; i<=7; i++) {
			Cube27[RotationAroundYTable[j][i]] = Edge_uc[j][i-6];
		}
	}
}

template<class _DataType>
void cThinning<_DataType>::RotateAroundZ(unsigned char *Cube27)
{
	int				i, j;
	unsigned char	Edge_uc[3][2];
	

	for (j=0; j<3; j++) {
		for (i=0; i<2; i++) {
			Edge_uc[j][i] = Cube27[RotationAroundZTable[j][i]];
		}
	}
	for (j=0; j<3; j++) {
		for (i=0; i<=5; i++) {
			Cube27[RotationAroundZTable[j][i]] = Cube27[RotationAroundZTable[j][i+2]];
		}
	}
	for (j=0; j<3; j++) {
		for (i=6; i<=7; i++) {
			Cube27[RotationAroundZTable[j][i]] = Edge_uc[j][i-6];
		}
	}
}


// 	ai = 17; bi = 14; ci = 8; di = 5; ei = 4; upi = 16; lpi = 10;
template<class _DataType>
void cThinning<_DataType>::SymmetryAlongX(unsigned char *Cube27)
{
	unsigned char Tempuc;
	SWAP(Cube27[0], Cube27[2], Tempuc);
	SWAP(Cube27[3], Cube27[5], Tempuc);
	SWAP(Cube27[6], Cube27[8], Tempuc);
	
	SWAP(Cube27[9], Cube27[11], Tempuc);
	SWAP(Cube27[12], Cube27[14], Tempuc);
	SWAP(Cube27[15], Cube27[17], Tempuc);
	
	SWAP(Cube27[18], Cube27[20], Tempuc);
	SWAP(Cube27[21], Cube27[23], Tempuc);
	SWAP(Cube27[24], Cube27[26], Tempuc);
}

template<class _DataType>
void cThinning<_DataType>::SymmetryAlongY(unsigned char *Cube27)
{
	unsigned char Tempuc;
	SWAP(Cube27[0], Cube27[6], Tempuc);
	SWAP(Cube27[1], Cube27[7], Tempuc);
	SWAP(Cube27[2], Cube27[8], Tempuc);
	
	SWAP(Cube27[9], Cube27[15], Tempuc);
	SWAP(Cube27[10], Cube27[16], Tempuc);
	SWAP(Cube27[11], Cube27[17], Tempuc);
	
	SWAP(Cube27[18], Cube27[24], Tempuc);
	SWAP(Cube27[19], Cube27[25], Tempuc);
	SWAP(Cube27[20], Cube27[26], Tempuc);
}

template<class _DataType>
void cThinning<_DataType>::SymmetryAlongZ(unsigned char *Cube27)
{
	int		i;
	unsigned char Tempuc;
	
	for (i=0; i<9; i++) {
		SWAP(Cube27[i], Cube27[i+18], Tempuc);
	}
}

template<class _DataType>
int	cThinning<_DataType>::Index(int X, int Y, int Z)
{
	if (X<0 || Y<0 || Z<0 || X>=Width_mi || Y>=Height_mi || Z>=Depth_mi) return 0;
	else return (Z*WtimesH_mi + Y*Width_mi + X);
}

template<class _DataType>
void cThinning<_DataType>::SaveThreeColorVolume(char *OutFileName, _DataType MatMin, _DataType MatMax)
{
	SaveThreeColorVolume(OutFileName, MatMin, MatMax, ThreeColorVolume_muc);
}


template<class _DataType>
void cThinning<_DataType>::SaveThreeColorVolume(char *OutFileName, _DataType MatMin, _DataType MatMax,
							unsigned char	*SkeletonVolume)
{
	int		i, binfile_fd2, Num1Voxels;
	char	FileName[200];
	unsigned char	*ThickSkeletons = new unsigned char [WHD_mi];


	sprintf (FileName, "%s_%03d_%03d.raw", OutFileName, (int)MatMin, (int)MatMax);
	cout << "Skeleton Volume File Name = " << FileName << endl;
	cout.flush();

	Num1Voxels = 0;
	for (i=0; i<WHD_mi; i++) {
		if (SkeletonVolume[i]) {
			ThickSkeletons[i] = (unsigned char)255;
			Num1Voxels++;
		}
		else ThickSkeletons[i] = (unsigned char)0;
	}
	cout << "The number of 1-voxels = " << Num1Voxels << endl;
	cout.flush();


	if ((binfile_fd2 = open (FileName, O_CREAT | O_WRONLY)) < 0) {
		printf ("could not open %s\n", FileName);
		exit(1);
	}
	if (write(binfile_fd2, ThickSkeletons, sizeof(unsigned char)*WHD_mi)
				!=sizeof(unsigned char)*WHD_mi) {
		printf ("The file could not be written : %s\n", FileName);
		close (binfile_fd2);
		exit(1);
	}
	close (binfile_fd2);

	if (chmod(FileName, S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH) < 0) {
		printf ("chmod was not worked to file %s\n", FileName);
		exit(1);
	}
	delete [] ThickSkeletons;
}
	
template<class _DataType>
void cThinning<_DataType>::SaveVolume(char *OutFileName, _DataType MatMin, _DataType MatMax,
							unsigned char	*VolumeInput)
{
	int		i, binfile_fd2, Num1Voxels;
	char	FileName[200];


	sprintf (FileName, "%s_%03d_%03d.raw", OutFileName, (int)MatMin, (int)MatMax);
	cout << "Volume File Name = " << FileName << endl;
	cout.flush();

	Num1Voxels = 0;
	for (i=0; i<WHD_mi; i++) {
		if (VolumeInput[i]>1) Num1Voxels++;
	}
	cout << "The number of 1-voxels = " << Num1Voxels << endl;
	cout.flush();


	if ((binfile_fd2 = open (FileName, O_CREAT | O_WRONLY)) < 0) {
		printf ("could not open %s\n", FileName);
		exit(1);
	}
	if (write(binfile_fd2, VolumeInput, sizeof(unsigned char)*WHD_mi)
				!=sizeof(unsigned char)*WHD_mi) {
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
void cThinning<_DataType>::SaveThickVolume(char *OutFileName, _DataType MatMin, _DataType MatMax, 
											unsigned char *VolumeInput)
{
	SaveThickVolume (OutFileName, MatMin, MatMax, (unsigned char)1);
}


template<class _DataType>
void cThinning<_DataType>::SaveThickVolume(char *OutFileName, _DataType MatMin, _DataType MatMax, 
											unsigned char *VolumeInput, unsigned char Threshold)
{
	int		i, j, k, l, m, n, loc[3], binfile_fd2, Num1Voxels;
	char	FileName[200];
	unsigned char	*ThickSkeletons = new unsigned char [WHD_mi];

	sprintf (FileName, "%s_Thick_%03d_%03d.raw", OutFileName, (int)MatMin, (int)MatMax);
	cout << "Thick Skeleton Volume File Name = " << FileName << endl;
	cout.flush();


	for (i=0; i<WHD_mi; i++) ThickSkeletons[i] = (unsigned char)0;


	Num1Voxels = 0;
	for (k=1; k<Depth_mi-1; k++) {
		for (j=1; j<Height_mi-1; j++) {
			for (i=1; i<Width_mi-1; i++) {

				// Drawing Thick Skeletons
				loc[0] = k*WtimesH_mi + j*Width_mi + i;
				if (VolumeInput[loc[0]]>=Threshold) {
					Num1Voxels++;
					for (n=-1; n<=1; n++) {
						for (m=-1; m<=1; m++) {
							for (l=-1; l<=1; l++) {
								loc[1] = (k+n)*WtimesH_mi + (j+m)*Width_mi + (i+l);
								ThickSkeletons[loc[1]] = (unsigned char)255;
							}
						}
					}
				}
				// The End of Drawing

			}
		}
	}
	cout << "The number of 1-voxels = " << Num1Voxels << endl;
	cout.flush();


	if ((binfile_fd2 = open (FileName, O_CREAT | O_WRONLY)) < 0) {
		printf ("could not open %s\n", FileName);
		exit(1);
	}
	if (write(binfile_fd2, ThickSkeletons, sizeof(unsigned char)*WHD_mi)
				!=sizeof(unsigned char)*WHD_mi) {
		printf ("The file could not be written : %s\n", FileName);
		close (binfile_fd2);
		exit(1);
	}
	close (binfile_fd2);

	if (chmod(FileName, S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH) < 0) {
		printf ("chmod was not worked to file %s\n", FileName);
		exit(1);
	}

	delete [] ThickSkeletons;

}
	
template<class _DataType>
void cThinning<_DataType>::SaveThickSkeletons(char *OutFileName, _DataType MatMin, _DataType MatMax)
{
	int		i, j, k, l, m, n, loc[3], binfile_fd2, Isolation, Num1Voxels;
	char	FileName[200];
	unsigned char	*ThickSkeletons = new unsigned char [WHD_mi];

	sprintf (FileName, "%s_Thick_%03d_%03d.raw", OutFileName, (int)MatMin, (int)MatMax);
	cout << "Thick Skeleton Volume File Name = " << FileName << endl;
	cout.flush();


	for (i=0; i<WHD_mi; i++) ThickSkeletons[i] = (unsigned char)0;

	for (k=1; k<Depth_mi-1; k++) {
		for (j=1; j<Height_mi-1; j++) {
			for (i=1; i<Width_mi-1; i++) {

				// Isolatione Check
				loc[0] = k*WtimesH_mi + j*Width_mi + i;
				Isolation = true;
				if (ThreeColorVolume_muc[loc[0]]==1) {
					for (n=-1; n<=1; n++) {
						for (m=-1; m<=1; m++) {
							for (l=-1; l<=1; l++) {
								loc[1] = (k+n)*WtimesH_mi + (j+m)*Width_mi + (i+l);
								if (ThreeColorVolume_muc[loc[1]]==1) {
									if (l==0 && m==0 && n==0) continue;
									else {
										Isolation = false;
										n=2; m=2; l=2; // break the three loops
									}
								}
							}
						}
					}
					if (Isolation) ThreeColorVolume_muc[loc[0]]=0;
				}
				
			}
		}
	}


	Num1Voxels = 0;
	for (k=1; k<Depth_mi-1; k++) {
		for (j=1; j<Height_mi-1; j++) {
			for (i=1; i<Width_mi-1; i++) {

				// Drawing Thick Skeletons
				loc[0] = k*WtimesH_mi + j*Width_mi + i;
				if (ThreeColorVolume_muc[loc[0]]) {
					Num1Voxels++;
					for (n=-1; n<=1; n++) {
						for (m=-1; m<=1; m++) {
							for (l=-1; l<=1; l++) {
								loc[1] = (k+n)*WtimesH_mi + (j+m)*Width_mi + (i+l);
								ThickSkeletons[loc[1]] = (unsigned char)255;
							}
						}
					}
				}
				// The End of Drawing

			}
		}
	}
	cout << "The number of 1-voxels = " << Num1Voxels << endl;
	cout.flush();


	if ((binfile_fd2 = open (FileName, O_CREAT | O_WRONLY)) < 0) {
		printf ("could not open %s\n", FileName);
		exit(1);
	}
	if (write(binfile_fd2, ThickSkeletons, sizeof(unsigned char)*WHD_mi)
				!=sizeof(unsigned char)*WHD_mi) {
		printf ("The file could not be written : %s\n", FileName);
		close (binfile_fd2);
		exit(1);
	}
	close (binfile_fd2);

	if (chmod(FileName, S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH) < 0) {
		printf ("chmod was not worked to file %s\n", FileName);
		exit(1);
	}

	delete [] ThickSkeletons;

}
	

template<class _DataType>
void cThinning<_DataType>::Destroy()
{
	delete []  ThreeColorVolume_muc;
	ThreeColorVolume_muc = NULL;
}

cThinning<unsigned char>	__ThinningValue0;
//cThinning<unsigned short>	__ThinningValue1;
//cThinning<int>				__ThinningValue2;
//cThinning<float>			__ThinningValue3;


