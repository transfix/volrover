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
#include <map.h>


#include <PEDetection/STree.h>


template<class _DataType>
cSTree_Eval<_DataType>::cSTree_Eval()
{
	NumWrongArteries_mi = 0;
	NumWrongVeins_mi = 0;
	for (int i=0; i<MAX_HISTO_ELEMENTS; i++) {
		HistoWrongA_mi[i] = 0;
		HistoWrongV_mi[i] = 0;
		HistoWrongM_mi[i] = 0;
		HistoWrongT_mi[i] = 0;
	}
}

template<class _DataType>
cSTree_Eval<_DataType>::~cSTree_Eval()
{

}

template<class _DataType>
void cSTree_Eval<_DataType>::setWHD(int W, int H, int D)
{
	Width_mi = W;
	Height_mi = H;
	Depth_mi = D;
	
	WtimesH_mi = W*H;
	WHD_mi = W*H*D;
}

template <class _DataType>
void cSTree_Eval<_DataType>::setData(_DataType *Data, float Minf, float Maxf)
{
	MinData_mf = Minf;
	MaxData_mf = Maxf;
	Data_mT = Data;
}

template <class _DataType>
void cSTree_Eval<_DataType>::setMuscleRange(int M1, int M2)
{
	Range_Muscles_mi[0] = M1;
	Range_Muscles_mi[1] = M2;
}

template <class _DataType>
void cSTree_Eval<_DataType>::setBVRange(int BV1, int BV2)
{
	Range_BloodVessels_mi[0] = BV1;
	Range_BloodVessels_mi[1] = BV2;
}


template<class _DataType>
void cSTree_Eval<_DataType>::ComputingPixelSpacing()
{

	int		NameLength_i = strlen(TargetName_gc);
	int		PatientNum_i;

	PatientNum_i = (TargetName_gc[NameLength_i-2] - '0')*10;
	PatientNum_i += (TargetName_gc[NameLength_i-1] - '0');
	printf ("cSTree Eval: Patient Num = %d\n", PatientNum_i);
	printf ("\n"); fflush (stdout);

	switch (PatientNum_i) {
		case 20: VoxelSize_gf = 0.742188; break;
		case 33: VoxelSize_gf = 0.742188; break;
		case 35: VoxelSize_gf = 0.781250; break;
		case 41: VoxelSize_gf = 0.742188; break;
		case 56: VoxelSize_gf = 0.703125; break;
		case 68: VoxelSize_gf = 0.644531; break;
		case 70: VoxelSize_gf = 0.644531; break;
		case 77: VoxelSize_gf = 0.703125; break;
		case 94: VoxelSize_gf = 0.742188; break;
		default: VoxelSize_gf = 0.742188; break;
	}


}


template<class _DataType>
void cSTree_Eval<_DataType>::AVSeparation_Evaluate(unsigned char *BV)
{
	int				i, StartLoc_i;
	cSkeletonVoxel	*Start_SVoxel;
	map<int, cSkeletonVoxel *>::iterator	StartLoc_it;


	ComputingPixelSpacing();


	BVVolume_muc = BV;
	
	TotalNumBranches_mi = 0;
	NumWrongClasss_mi = 0;
	NumMissing_mi = 0;
	for (i=0; i<20; i++) NumWrongBranch_Generations_mi[i] = 0;	
	
	if ((int)Skeleton_mmap.size()==0) {
		printf ("Evalution: Skeleton gmap size = 0\n"); fflush (stdout);
	}
	else {
		StartLoc_it = StartLoc_mmap.begin();
		for (i=0; i<(int)StartLoc_mmap.size(); i++, StartLoc_it++) {
			StartLoc_i = (*StartLoc_it).first;
			Start_SVoxel = (*StartLoc_it).second;
			if (Start_SVoxel->AV_c=='A') CurrAVType_mc = 'A';
			else if (Start_SVoxel->AV_c=='V') CurrAVType_mc = 'V';

			Evaluation(StartLoc_i);
			
		}
	}
	printf ("Total # Branches = %d\n", 	TotalNumBranches_mi);
	
	printf ("# Wrong Branches = %d, ", 	NumWrongClasss_mi);
	printf ("%.2f%%\n", (float)NumWrongClasss_mi/TotalNumBranches_mi*100);
	
	printf ("# Wrong Arteries = %d, ", 	NumWrongArteries_mi);
	printf ("%.2f%%\n", (float)NumWrongArteries_mi/TotalNumBranches_mi*100);
	
	printf ("# Wrong Veins = %d, ", 	NumWrongVeins_mi);
	printf ("%.2f%%\n", (float)NumWrongVeins_mi/TotalNumBranches_mi*100);
	
	printf ("# Missing Branches = %d, ", 	NumMissing_mi);
	printf ("%.2f%%\n", (float)NumMissing_mi/TotalNumBranches_mi*100);
	
	printf ("# Wrong + MIssing Branches = %d, ", 	NumWrongClasss_mi+NumMissing_mi);
	printf ("%.2f%%\n", (float)(NumMissing_mi+NumWrongClasss_mi)/TotalNumBranches_mi*100);
	
	printf ("Total # Dead Ends = %d\n", (int)SkeletonDeadEnds_mmap.size());
	printf ("Total # Voxels = %d\n", (int)Skeleton_mmap.size());

	printf ("Wrong Branch Generations:\n"); 
	double	AveGenerations_d;
	int		TotalGenerations_i = 0, NumWrongBranches_i = 0;
	
	
	for (i=0; i<20; i++) {
		printf ("Ith Generation = %2d, # Branches = %d\n", i, NumWrongBranch_Generations_mi[i]);

		TotalGenerations_i += i*NumWrongBranch_Generations_mi[i];
		NumWrongBranches_i += NumWrongBranch_Generations_mi[i];
	}
	printf ("TotalGenerations = %d ", TotalGenerations_i);
	printf ("# Wrong Branches = %d\n", NumWrongBranches_i);
	AveGenerations_d = (double)TotalGenerations_i/NumWrongBranches_i;
	printf ("Ave. Generations = %.4f\n", AveGenerations_d);
	fflush (stdout);

	int		MaxR_i = 12;
	printf ("\nWrong Branch Histogram ... \n");
	for (i=0; i<MaxR_i; i++) {
		printf ("Wrong Branch [%2d] = %3d\n", i, HistoWrongT_mi[i]);
	}

	printf ("\nMissing Branch Histogram ... \n");
	for (i=0; i<MaxR_i; i++) {
		printf ("Missing Branch [%2d] = %3d\n", i, HistoWrongM_mi[i]);
	}

	printf ("\nWrong Artery Branch Histogram ... \n");
	for (i=0; i<MaxR_i; i++) {
		printf ("Wrong Artery Branch [%2d] = %3d\n", i, HistoWrongA_mi[i]);
	}

	printf ("\nWrong Vein Branch Histogram ... \n");
	for (i=0; i<MaxR_i; i++) {
		printf ("Wrong Vein Branch [%2d] = %3d\n", i, HistoWrongV_mi[i]);
	}
}

template<class _DataType>
void cSTree_Eval<_DataType>::Evaluation(int CurrLoc)
{
	int				i, loc[3], CurrLoc_i, NextLoc_i, NumNextVoxels_i;
	int				CurrPt_i[3], NumWrongVoxels_i, NumTotalVoxels_i;
	int				NumWrongA_i, NumWrongV_i;
	int				MaxBVR_i, MaxBVMR_i;
	int				HistoIndex_i;
	cSkeletonVoxel	*CurrSVoxel=NULL;
	map<int, cSkeletonVoxel *>::iterator	Skeleton_it;
	cStack<int>		WrongVoxels_stack;
	
	
	CurrLoc_i = CurrLoc;
	NumWrongVoxels_i = 0;
	NumWrongA_i = 0;
	NumWrongV_i = 0;
	NumTotalVoxels_i = 0;

	TotalNumBranches_mi++;
	
	do {

		if (NumRepeat_Eval_mi++>50000) break;	// To prevent infinite loops
		

	 	Skeleton_it = Skeleton_mmap.find(CurrLoc_i);
		if (Skeleton_it==Skeleton_mmap.end()) break;
		
		CurrSVoxel = (*Skeleton_it).second;
		CurrSVoxel->getXYZ(CurrPt_i[0], CurrPt_i[1], CurrPt_i[2]);
		// Flip X
		loc[0] = CurrPt_i[2]*Width_mi*Height_mi + CurrPt_i[1]*Width_mi + (Width_mi - CurrPt_i[0] - 1);
		if (CurrAVType_mc=='A') {
			if (BVVolume_muc[loc[0]]!=255) {
				NumWrongVoxels_i++;
				NumWrongV_i++;
				WrongVoxels_stack.Push(CurrLoc_i);
			}
		}
		else if (CurrAVType_mc=='V') {
			if (BVVolume_muc[loc[0]]!=230 && BVVolume_muc[loc[0]]!=240) {
				NumWrongA_i++;
				NumWrongVoxels_i++;
				WrongVoxels_stack.Push(CurrLoc_i);
			}
		}
		NumTotalVoxels_i++;
		
		
		
		NumNextVoxels_i = CurrSVoxel->getNumNext();
		if (NumNextVoxels_i==0) break;
		else if (NumNextVoxels_i==1) {
			CurrLoc_i = CurrSVoxel->NextVoxels_i[0];
		}
		else {
			for (i=0; i<NumNextVoxels_i; i++) {
				NextLoc_i = CurrSVoxel->NextVoxels_i[i];
				if (NextLoc_i<0) break;
				Evaluation(NextLoc_i);
			}
			break;
		}
		
	} while (1);
	


	double	AveDiameterBV_d = 0.0;
	double	AveDiameterBVM_d = 0.0;
	
	if (NumWrongVoxels_i==NumTotalVoxels_i) {
		printf ("\nFound Wrong Branch: # Wrong Voxels = %d, ", NumWrongVoxels_i);
		printf ("# Wrong Arteries = %d, ", NumWrongA_i);
		printf ("# Wrong Veins = %d\n", NumWrongV_i);
		
		for (i=0; i<WrongVoxels_stack.Size(); i++) {
			WrongVoxels_stack.IthValue(i, CurrLoc_i);

			printf ("\tID = %8d ", CurrLoc_i);

			IndexInverse(CurrLoc_i, CurrPt_i[0], CurrPt_i[1], CurrPt_i[2]);
			loc[0] = CurrPt_i[2]*Width_mi*Height_mi + CurrPt_i[1]*Width_mi + (Width_mi - CurrPt_i[0] - 1);
			printf ("AV Value = %3d: ", BVVolume_muc[loc[0]]);
			
			Skeleton_it = Skeleton_mmap.find(CurrLoc_i);
			if (Skeleton_it==Skeleton_mmap.end()) {
				printf ("%3d %3d %3d\n", CurrPt_i[0], CurrPt_i[1], CurrPt_i[2]); fflush (stdout);
				continue;
			}

			MaxBVR_i = RecomputingMaxR_InBloodVessels((Width_mi-CurrPt_i[0]-1), CurrPt_i[1], CurrPt_i[2]);
			printf ("BVR = %2d ", MaxBVR_i);
			AveDiameterBV_d += (2*MaxBVR_i + 1)*VoxelSize_gf;
			
			MaxBVMR_i = RecomputingMaxR_InBV_Muscles((Width_mi-CurrPt_i[0]-1), CurrPt_i[1], CurrPt_i[2]);
			printf ("BVMR = %2d ", MaxBVMR_i);
			AveDiameterBVM_d += (2*MaxBVMR_i + 1)*VoxelSize_gf;


			CurrSVoxel = (*Skeleton_it).second;
			CurrSVoxel->Display();
			
		}
		printf ("# Wrong Arteries = %d, ", NumWrongA_i);
		printf ("# Wrong Veins = %d, ", NumWrongV_i);
		AveDiameterBV_d /= (double)NumWrongVoxels_i;
		AveDiameterBVM_d /= (double)NumWrongVoxels_i;
		printf ("Ave. Diameter BV  = %8.4f\n", AveDiameterBV_d);
		printf ("Ave. Diameter BVM = %8.4f\n", AveDiameterBVM_d);
		printf ("\n"); fflush (stdout);

		NumWrongBranch_Generations_mi[CurrSVoxel->Generation_i] += 1;
		HistoIndex_i = (int)((AveDiameterBV_d + AveDiameterBVM_d)/2.0 + 0.5);
		
		
		
		if (BVVolume_muc[loc[0]]>0) {
			NumWrongClasss_mi++;
			HistoWrongT_mi[HistoIndex_i]++;
			if (NumWrongA_i>NumWrongV_i) {
				NumWrongArteries_mi++;
				HistoWrongA_mi[HistoIndex_i]++;
			}
			else {
				NumWrongVeins_mi++;
				HistoWrongV_mi[HistoIndex_i]++;
			}
		}
		if (BVVolume_muc[loc[0]]==0) {
			NumMissing_mi++;
			HistoWrongM_mi[HistoIndex_i]++;
		}
	}
	
}


template<class _DataType>
int cSTree_Eval<_DataType>::RecomputingMaxR_InBV_Muscles(int Xi, int Yi, int Zi)
{
//	int				i, CurrSpR_i, CurrSpID_i;
	int				j, l, loc[3], DX, DY, DZ, *SphereIndex_i;
	int				SpR_i, NewCenter_i[3], NewSpR_i;
	int				NumVoxels, MaxBVMR_i;
	cStack<int>		VoxelLocs_stack;
	

	SpR_i = 3;
		
	VoxelLocs_stack.setDataPointer(0);
	for (j=0; j<=SpR_i; j++) {
		SphereIndex_i = getSphereIndex(j, NumVoxels);
		if (NumVoxels==0) break;

		for (l=0; l<NumVoxels; l++) {
			DX = SphereIndex_i[l*3 + 0];
			DY = SphereIndex_i[l*3 + 1];
			DZ = SphereIndex_i[l*3 + 2];
			loc[0] = Index (Xi+DX, Yi+DY, Zi+DZ);
			if (Data_mT[loc[0]]>=Range_Muscles_mi[0] &&
				Data_mT[loc[0]]<=Range_BloodVessels_mi[1]) {
				VoxelLocs_stack.Push(loc[0]);
			}
		}
	}

	ComputingTheBiggestSphereAt_DataRange(VoxelLocs_stack, 
						Range_Muscles_mi[0], Range_BloodVessels_mi[1], &NewCenter_i[0], NewSpR_i);

	if (NewSpR_i<0) MaxBVMR_i = 0;
	else MaxBVMR_i = NewSpR_i;

	return MaxBVMR_i;

}


template<class _DataType>
int cSTree_Eval<_DataType>::RecomputingMaxR_InBloodVessels(int Xi, int Yi, int Zi)
{
//	int				i, CurrSpR_i, CurrSpID_i;
	int				j, l, loc[3], DX, DY, DZ, *SphereIndex_i;
	int				SpR_i, NewCenter_i[3], NewSpR_i;
	int				NumVoxels, MaxBVR_i;
	cStack<int>		VoxelLocs_stack;
	

	SpR_i = 3;
		
	VoxelLocs_stack.setDataPointer(0);
	for (j=0; j<=SpR_i; j++) {
		SphereIndex_i = getSphereIndex(j, NumVoxels);
		if (NumVoxels==0) break;

		for (l=0; l<NumVoxels; l++) {
			DX = SphereIndex_i[l*3 + 0];
			DY = SphereIndex_i[l*3 + 1];
			DZ = SphereIndex_i[l*3 + 2];
			loc[0] = Index (Xi+DX, Yi+DY, Zi+DZ);
			if (Data_mT[loc[0]]>=Range_BloodVessels_mi[0] &&
				Data_mT[loc[0]]<=Range_BloodVessels_mi[1]) {
				VoxelLocs_stack.Push(loc[0]);
			}
		}
	}

	ComputingTheBiggestSphereAt_DataRange(VoxelLocs_stack, 
						Range_BloodVessels_mi[0], Range_BloodVessels_mi[1], &NewCenter_i[0], NewSpR_i);

	if (NewSpR_i<0) MaxBVR_i = 0;
	else MaxBVR_i = NewSpR_i;

	return MaxBVR_i;
}


template<class _DataType>
int cSTree_Eval<_DataType>::ComputingTheBiggestSphereAt_DataRange(cStack<int> &VoxelLocs, 
								_DataType Lower_Th, _DataType Upper_Th, 
								int *Center3_ret, int &SphereR_ret)
{
	int			i, l, loc[3], Xi, Yi, Zi, DX, DY, DZ, *SphereIndex_i;
	int			FoundSphere_i, SphereR_i, NumVoxels, ReturnNew_i;


	ReturnNew_i = false;
	SphereR_ret = -1;
	for (i=0; i<VoxelLocs.Size(); i++) {
		VoxelLocs.IthValue(i, loc[0]);
		if ((Data_mT[loc[0]] < Lower_Th || Data_mT[loc[0]] > Upper_Th)) continue;
		IndexInverse(loc[0], Xi, Yi, Zi);

		FoundSphere_i = false;
		SphereR_i = -1;
		
		do {
		
			SphereR_i++;
			SphereIndex_i = getSphereIndex(SphereR_i, NumVoxels);
			if (NumVoxels==0) {
				SphereR_i--;
				break;
			}
			for (l=0; l<NumVoxels; l++) {
				DX = SphereIndex_i[l*3 + 0];
				DY = SphereIndex_i[l*3 + 1];
				DZ = SphereIndex_i[l*3 + 2];
				loc[1] = Index (Xi+DX, Yi+DY, Zi+DZ);
				if (Data_mT[loc[1]] < Lower_Th || Data_mT[loc[1]] > Upper_Th) {
					SphereR_i--;
					FoundSphere_i = true;
					break;
				}
			}
			
		} while (FoundSphere_i==false);

		if (SphereR_ret < SphereR_i) {
			SphereR_ret = SphereR_i;
			Center3_ret[0] = Xi;
			Center3_ret[1] = Yi;
			Center3_ret[2] = Zi;
			ReturnNew_i = true;
		}
		
	}

	return ReturnNew_i;
}



template<class _DataType>
void cSTree_Eval<_DataType>::LoadTreeStructure(char *SkeletonFileName_c)
{
	int				i, j, CurrLoc_i, NumSVoxels_i, NumNext_i;
	int				NumStartLocs_i=0, StartLoc_i;
	float			Dist_f;
	cSkeletonVoxel	*SVoxel, *DeadEnd_SVoxel;
	map<int, cSkeletonVoxel *>::iterator	Skeleton_it;
	

	FILE	*Skeleton_fs;
	Skeleton_fs = fopen (SkeletonFileName_c, "r");
	if (Skeleton_fs==NULL) {
		printf("Cannot open the skeleton file %s\n", SkeletonFileName_c); 
		fflush (stdout);
	}
	else {
		printf ("Loading ... %s\n", SkeletonFileName_c);

		fscanf (Skeleton_fs, "%d", &NumStartLocs_i);
		printf ("Num Start Locs = %d\n", NumStartLocs_i); fflush (stdout);
		StartLoc_mmap.clear();
		for (i=0; i<NumStartLocs_i; i++) {
			fscanf (Skeleton_fs, "%d", &StartLoc_i);
			SVoxel = new cSkeletonVoxel;
			fscanf (Skeleton_fs, "%d %d %d", &SVoxel->Xi, &SVoxel->Yi, &SVoxel->Zi);
			fscanf (Skeleton_fs, "%d ", &SVoxel->Generation_i);
			fscanf (Skeleton_fs, "%c", &SVoxel->LR_c);
			fscanf (Skeleton_fs, "%c", &SVoxel->AV_c);
			fscanf (Skeleton_fs, "%c", &SVoxel->End_c);
			
			printf ("StartLoc = %d ", StartLoc_i);
			printf ("XYZ = %3d %3d %3d ", SVoxel->Xi, SVoxel->Yi, SVoxel->Zi);
			printf ("Generation = %d ", SVoxel->Generation_i);
			printf ("LR = %c ", SVoxel->LR_c);
			printf ("AV = %c ", SVoxel->AV_c);
			printf ("\n"); fflush (stdout);
			
			StartLoc_mmap[StartLoc_i] = SVoxel;
		}


		fscanf (Skeleton_fs, "%d", &NumSVoxels_i);
		printf ("Total num of SVoxels = %d\n", NumSVoxels_i); fflush (stdout);
		
		Dist_f = 0;
		Skeleton_mmap.clear();
		for (i=0; i<NumSVoxels_i; i++) {
			SVoxel = new cSkeletonVoxel;
			fscanf (Skeleton_fs, "%d", &CurrLoc_i);
			fscanf (Skeleton_fs, "%d %d %d", &SVoxel->Xi, &SVoxel->Yi, &SVoxel->Zi);
			fscanf (Skeleton_fs, "%d ", &SVoxel->Generation_i);
			fscanf (Skeleton_fs, "%c", &SVoxel->LR_c);
			fscanf (Skeleton_fs, "%c", &SVoxel->AV_c);
			fscanf (Skeleton_fs, "%d", &SVoxel->PrevVoxel_i);
			fscanf (Skeleton_fs, "%d", &NumNext_i);
			if (NumNext_i==0) SVoxel->End_c = 'E';
			
			for (j=0; j<NumNext_i; j++) {
				fscanf (Skeleton_fs, "%d", &SVoxel->NextVoxels_i[j]);
			}
			Skeleton_mmap[CurrLoc_i] = SVoxel;
			/*
			printf ("Curr Loc = %d ", CurrLoc_i);
			printf ("XYZ = %3d %3d %3d ", SVoxel->Xi, SVoxel->Yi, SVoxel->Zi);
			printf ("Generation = %2d ", SVoxel->Generation_i);
			printf ("%c%c%c ", SVoxel->LR_c,SVoxel->AV_c, SVoxel->End_c);
			printf ("Prev = %d ", SVoxel->PrevVoxel_i);
			printf ("#Next = %d ", NumNext_i);
			printf ("\n"); fflush (stdout);
			*/

			if (NumNext_i==0) {
				Dist_f += 0.1;
				DeadEnd_SVoxel = new cSkeletonVoxel;
				DeadEnd_SVoxel->Copy(SVoxel);
				SkeletonDeadEnds_mmap[Dist_f] = DeadEnd_SVoxel;
			}
		}
		fclose (Skeleton_fs);
		printf ("SkeletonDeadEnds gmap size = %d\n", (int)SkeletonDeadEnds_mmap.size());
		printf ("Skeleton gmap size = %d\n", (int)Skeleton_mmap.size());
		printf ("Loading is done\n"); fflush (stdout);
	}
}


template<class _DataType>
int *cSTree_Eval<_DataType>::getSphereIndex(int SphereRadius, int &NumVoxels_ret)
{
	int 	*SphereIndex_ret;
	switch (SphereRadius) {
		case 0: SphereIndex_ret = &SphereR00_gi[0]; NumVoxels_ret = NumSphereR00_gi; break;
		case 1: SphereIndex_ret = &SphereR01_gi[0]; NumVoxels_ret = NumSphereR01_gi; break;
		case 2: SphereIndex_ret = &SphereR02_gi[0]; NumVoxels_ret = NumSphereR02_gi; break;
		case 3: SphereIndex_ret = &SphereR03_gi[0]; NumVoxels_ret = NumSphereR03_gi; break;
		case 4: SphereIndex_ret = &SphereR04_gi[0]; NumVoxels_ret = NumSphereR04_gi; break;
		case 5: SphereIndex_ret = &SphereR05_gi[0]; NumVoxels_ret = NumSphereR05_gi; break;
		case 6: SphereIndex_ret = &SphereR06_gi[0]; NumVoxels_ret = NumSphereR06_gi; break;
		case 7: SphereIndex_ret = &SphereR07_gi[0]; NumVoxels_ret = NumSphereR07_gi; break;
		case 8: SphereIndex_ret = &SphereR08_gi[0]; NumVoxels_ret = NumSphereR08_gi; break;
		case 9: SphereIndex_ret = &SphereR09_gi[0]; NumVoxels_ret = NumSphereR09_gi; break;
		case 10:SphereIndex_ret = &SphereR10_gi[0]; NumVoxels_ret = NumSphereR10_gi; break;
		case 11:SphereIndex_ret = &SphereR11_gi[0]; NumVoxels_ret = NumSphereR11_gi; break;
		case 12:SphereIndex_ret = &SphereR12_gi[0]; NumVoxels_ret = NumSphereR12_gi; break;
		case 13:SphereIndex_ret = &SphereR13_gi[0]; NumVoxels_ret = NumSphereR13_gi; break;
		case 14:SphereIndex_ret = &SphereR14_gi[0]; NumVoxels_ret = NumSphereR14_gi; break;
		case 15:SphereIndex_ret = &SphereR15_gi[0]; NumVoxels_ret = NumSphereR15_gi; break;
		case 16:SphereIndex_ret = &SphereR16_gi[0]; NumVoxels_ret = NumSphereR16_gi; break;
		case 17:SphereIndex_ret = &SphereR17_gi[0]; NumVoxels_ret = NumSphereR17_gi; break;
		case 18:SphereIndex_ret = &SphereR18_gi[0]; NumVoxels_ret = NumSphereR18_gi; break;
		case 19:SphereIndex_ret = &SphereR19_gi[0]; NumVoxels_ret = NumSphereR19_gi; break;
		case 20:SphereIndex_ret = &SphereR20_gi[0]; NumVoxels_ret = NumSphereR20_gi; break;
		case 21:SphereIndex_ret = &SphereR21_gi[0]; NumVoxels_ret = NumSphereR21_gi; break;
		case 22:SphereIndex_ret = &SphereR22_gi[0]; NumVoxels_ret = NumSphereR22_gi; break;
		case 23:SphereIndex_ret = &SphereR23_gi[0]; NumVoxels_ret = NumSphereR23_gi; break;
		case 24:SphereIndex_ret = &SphereR24_gi[0]; NumVoxels_ret = NumSphereR24_gi; break;
		case 25:SphereIndex_ret = &SphereR25_gi[0]; NumVoxels_ret = NumSphereR25_gi; break;
		case 26:SphereIndex_ret = &SphereR26_gi[0]; NumVoxels_ret = NumSphereR26_gi; break;
		default: SphereIndex_ret = NULL; NumVoxels_ret = 0; break;
	}
	return SphereIndex_ret;
}



template<class _DataType>
int cSTree_Eval<_DataType>::Index(int X, int Y, int Z)
{
	if (X<0 || Y<0 || Z<0 || X>=Width_mi || Y>=Height_mi || Z>=Depth_mi) return 0;
	else return (Z*WtimesH_mi + Y*Width_mi + X);
}

template<class _DataType>
int cSTree_Eval<_DataType>::Index(int *Loc3)
{
	if (Loc3[0]<0 || Loc3[1]<0 || Loc3[2]<0 || Loc3[0]>=Width_mi || 
		Loc3[1]>=Height_mi || Loc3[2]>=Depth_mi) return 0;
	else return (Loc3[2]*WtimesH_mi + Loc3[1]*Width_mi + Loc3[0]);
}

template<class _DataType>
void cSTree_Eval<_DataType>::IndexInverse(int Loc, int &X, int &Y, int &Z)
{
	if (Loc<0 || Loc>=WHD_mi) {
		X = Y = Z = 0;
	}
	else {
		Z = Loc/WtimesH_mi;
		Y = (Loc - Z*WtimesH_mi)/Width_mi;
		X = Loc % Width_mi;
	}
}

template<class _DataType>
void cSTree_Eval<_DataType>::IndexInverse(int Loc, int *Center3_ret)
{
	if (Loc<0 || Loc>=WHD_mi) {
		Center3_ret[0] = Center3_ret[1] = Center3_ret[2] = 0;
	}
	else {
		Center3_ret[2] = Loc/WtimesH_mi;
		Center3_ret[1] = (Loc - Center3_ret[2]*WtimesH_mi)/Width_mi;
		Center3_ret[0] = Loc % Width_mi;
	}
}


template class cSTree_Eval<unsigned char>;
//template class cSTree_Eval<unsigned short>;
//template class cSTree_Eval<int>;
//template class cSTree_Eval<float>;
