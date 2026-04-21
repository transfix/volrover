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
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <stdlib.h>
#include <fstream>
#include <cassert>
#include <VolMagick/VolMagick.h>
#include <Mesher.h>
#include <UniformGrid.h>
#include <Mesh.h>
#include <MMHLS/generateMesh.h>
#include <FastContouring/FastContouring.h>
#include <QFileInfo>
#include <QString>

using namespace std;
void correctMeshes(vector<string> &vols,vector<string> &procvols,vector<Mesh> &meshes,vector<float> &isovalues,float tolerance);
void mopUpDirt(Mesh &m, float threshold);
void ReSamplingCubicSpline1OrderPartials(float *coeff, int nx, int ny, int nz, float *dxyz, float *minExtent, float *zeroToTwoPartialValue,  float *p) ;

struct ManifestEntry
{
  unsigned int materialId;
  float isovalue;
  float color[3];
};

void triSurfToMesh(const FastContouring::TriSurf &triSurf, Mesh &mesh)
{
    for(int i=0;i<triSurf.verts.size();i+=3)
    {
      Point p;
      p.x=triSurf.verts[i]; p.y=triSurf.verts[i+1]; p.z=triSurf.verts[i+2];
      mesh.vertices.push_back(p);
    }
    
    for(int i=0;i<triSurf.tris.size();i+=3)
    {
      vector<unsigned int> t;
      t.push_back(triSurf.tris[i]);t.push_back(triSurf.tris[i+1]);t.push_back(triSurf.tris[i+2]);
      mesh.triangles.push_back(t);
    }
}
void readManifestFile(char *fileName,map<string,ManifestEntry> &manifest)
{
  fstream manifestFile;
  manifestFile.open(fileName,ios::in);
  int n;
  manifestFile>>n;

  manifestFile.ignore(1); // To chuck that pesky \n from the stream
  for(int i=0;i<n;i++)
  {
    char volFileName[300],isovalue[20],matId[10],red[20],green[20],blue[20];
    ManifestEntry mEntry;
    manifestFile.getline(volFileName,299,',');
    manifestFile.getline(matId,9,',');
    manifestFile.getline(isovalue,20,',');
    manifestFile.getline(red,19,',');
    manifestFile.getline(green,19,',');
    manifestFile.getline(blue,19);

    cout<<"File "<<volFileName<<endl
        <<"Material: "<<matId<<endl
        <<"Isovalue: "<<isovalue<<endl
        <<"Color: ("<<red<<","<<green<<","<<blue<<")"<<endl;

    mEntry.materialId=atoi(matId);
    mEntry.isovalue=atof(isovalue);
    mEntry.color[0]=atof(red);
    mEntry.color[1]=atof(green);
    mEntry.color[2]=atof(blue);

    manifest[string(volFileName)]=mEntry;
  }

  manifestFile.close();
}


/*

int main(int argc,char *argv[])
{
  if(argc<8)
  {
    cout<<"generateMesh <manifest-file> <iso-value> <tolerance> <volume-thresh> <begin-mesh> <end-mesh> <mesh-name-prefix>"<<endl;
    return 0;
  }

  map<string,ManifestEntry> manifest;
  cout<<"Reading manifest file."<<endl;
  readManifestFile(argv[1],manifest);

  float scalingFactor=atof(argv[2]);

  float tolerance=atof(argv[3]);

  int begin=atoi(argv[5]);
  int end=atoi(argv[6]);
  
  vector<float> isovalues;

  cout<<"Scaling Factor is: "<<scalingFactor<<" tolerance: "<<tolerance<<endl;

  vector<Mesh> meshes;
  vector<string> volumeNames;
  vector<string> processedVolumes;
  vector<unsigned int> matIds;
  int count=-1;
  // For every volume indicated in manifest, generate mesh.
  for(map<string,ManifestEntry>::iterator mIter=manifest.begin();mIter!=manifest.end();mIter++)
  {
    //cout<<mIter->second.color[0]<<","<<mIter->second.color[1]<<","<<mIter->second.color[2]<<endl;
    volumeNames.push_back(mIter->first);
    count++;
    if(count<begin&&count>end) continue;
    Mesh mesh(mIter->second.color);

    processedVolumes.push_back(mIter->first);

    cout<<"Processing volume "<<mIter->first<<endl;
    
    matIds.push_back(mIter->second.materialId);

    VolMagick::Volume vol;
    VolMagick::readVolumeFile(cvcapp, vol,mIter->first);
    
    isovalues.push_back(mIter->second.isovalue*scalingFactor);

//     Mesher mesher;
//     UniformGrid ugrid(1,1,1);
//     ugrid.importVolume(vol);
// 
//     mesher.setGrid((Grid &)(ugrid));
//     mesher.setMesh(&mesh);
    cout<<"Generating mesh for isovalue "<<mIter->second.isovalue*scalingFactor<<endl;
//     mesher.generateMesh(mIter->second.isovalue*scalingFactor);
    //mesh.reorientMesh();
     FastContouring::ContourExtractor *contexta;
     contexta=new FastContouring::ContourExtractor;

     contexta->setVolume(vol);
     FastContouring::TriSurf result=contexta->extractContour(mIter->second.isovalue*scalingFactor);
     triSurfToMesh(result,mesh);
    meshes.push_back(mesh);
    //scount++;
     delete contexta;
  }

  cout<<"Finished generating meshes"<<endl;
  cout<<"Number of meshes: "<<meshes.size()<<endl;
  // compute bounding boxes.
  for(int i=0;i<meshes.size();i++) 
  {
    cout<<"Computing bb for #"<<i<<endl;
    cout<<"Num verts: "<<meshes[i].vertices.size()<<endl;
    cout<<"Num tris: "<<meshes[i].triangles.size()<<endl;
    cout<<"Corresponding volume: "<<volumeNames[i]<<endl;
    meshes[i].computeBoundingBox();
  }
  
  cout<<"Computed bounding boxes."<<endl;
  //Cleaning meshes acc threshold.
  cout<<"Begin mesh mop-up."<<endl;
  float threshold=atof(argv[4]);
  for(int i=0;i<meshes.size();i++) 
  {
    cout<<"Cleaning mesh with matId: "<<matIds[i]<<endl;
    mopUpDirt(meshes[i],threshold);
  }
  
//   cout<<"Saving meshes before correction."<<endl;
//   for(int i=0;i<meshes.size();i++)
//   {
//     stringstream s;
//     s<<argv[7]<<"-before-"<<matIds[i]<<".rawc";
//     meshes[i].saveMesh(s.str());
//   }
  // Now we have generated the meshes. Call the appropriate correct mesh function - modify correct meshes.

  cout<<"Begin mesh correction."<<endl;
  correctMeshes(volumeNames,processedVolumes,meshes,isovalues,tolerance);

  cout<<"Finished mesh correction."<<endl;
  // Now we have to save the meshes.
  cout<<"Saving meshes."<<endl;
  for(int i=0;i<meshes.size();i++)
  {
    stringstream s;
    s<<argv[7]<<"-"<<matIds[i]<<".rawc";
    meshes[i].saveMesh(s.str());
  }

  cout<<"All done."<<endl;
  return 0;
}

*/

void generateMesh(string manifestFile, float isoratio, float tolerance, float volthresh, int meshStart, int meshEnd, string outPref)
{


  map<string,ManifestEntry> manifest;
  cout<<"Reading manifest file."<<endl;
  readManifestFile((char*)manifestFile.c_str(), manifest);

  float scalingFactor=isoratio;


  int begin=meshStart;
  int end=meshEnd;
  
  vector<float> isovalues;

  cout<<"Scaling Factor is: "<<scalingFactor<<" tolerance: "<<tolerance<<endl;

  vector<Mesh> meshes;
  vector<string> volumeNames;
  vector<string> processedVolumes;
  vector<unsigned int> matIds;
  int count=0;
  // For every volume indicated in manifest, generate mesh.
  for(map<string,ManifestEntry>::iterator mIter=manifest.begin();mIter!=manifest.end();mIter++)
  {
    //cout<<mIter->second.color[0]<<","<<mIter->second.color[1]<<","<<mIter->second.color[2]<<endl;
    volumeNames.push_back(mIter->first);
    count++;
    if(count<begin || count>end) continue;
    Mesh mesh(mIter->second.color);

    processedVolumes.push_back(mIter->first);

    cout<<"Processing volume "<<mIter->first<<endl;
    
    matIds.push_back(mIter->second.materialId);

    VolMagick::Volume vol;
    VolMagick::readVolumeFile(cvcapp, vol,mIter->first);
    
    isovalues.push_back(mIter->second.isovalue*scalingFactor);

//     Mesher mesher;
//     UniformGrid ugrid(1,1,1);
//     ugrid.importVolume(vol);
// 
//     mesher.setGrid((Grid &)(ugrid));
//     mesher.setMesh(&mesh);
    cout<<"Generating mesh for isovalue "<<mIter->second.isovalue*scalingFactor<<endl;
//     mesher.generateMesh(mIter->second.isovalue*scalingFactor);
    //mesh.reorientMesh();
     FastContouring::ContourExtractor *contexta;
     contexta=new FastContouring::ContourExtractor;

     contexta->setVolume(vol);
     FastContouring::TriSurf result=contexta->extractContour(mIter->second.isovalue*scalingFactor);
     triSurfToMesh(result,mesh);
    meshes.push_back(mesh);
    //scount++;
     delete contexta;
  }

  cout<<"Finished generating meshes"<<endl;
  cout<<"Number of meshes: "<<meshes.size()<<endl;
  // compute bounding boxes.
  for(int i=0;i<meshes.size();i++) 
  {
    cout<<"Computing bb for #"<<i<<endl;
    cout<<"Num verts: "<<meshes[i].vertices.size()<<endl;
    cout<<"Num tris: "<<meshes[i].triangles.size()<<endl;
    cout<<"Corresponding volume: "<<volumeNames[i]<<endl;
    meshes[i].computeBoundingBox();
  }
  
  cout<<"Computed bounding boxes."<<endl;
  //Cleaning meshes acc threshold.
  cout<<"Begin mesh mop-up."<<endl;
  float threshold=volthresh;
  for(int i=0;i<meshes.size();i++) 
  {
    cout<<"Cleaning mesh with matId: "<<matIds[i]<<endl;
    mopUpDirt(meshes[i],threshold);
  }
  
//   cout<<"Saving meshes before correction."<<endl;
//   for(int i=0;i<meshes.size();i++)
//   {
//     stringstream s;
//     s<<argv[7]<<"-before-"<<matIds[i]<<".rawc";
//     meshes[i].saveMesh(s.str());
//   }
  // Now we have generated the meshes. Call the appropriate correct mesh function - modify correct meshes.

  cout<<"Begin mesh correction."<<endl;
  correctMeshes(volumeNames,processedVolumes,meshes,isovalues,tolerance);

  cout<<"Finished mesh correction."<<endl;
  // Now we have to save the meshes.
  cout<<"Saving meshes."<<endl;

  QFileInfo fi = QString::fromStdString(manifestFile);
  QString name = fi.fileName();

  string outDir = manifestFile;
  outDir.erase(outDir.length()-name.toStdString().length(), name.toStdString().length());


  for(int i=0;i<meshes.size();i++)
  {
    stringstream s;
    s<<outDir<<outPref<<"-"<<(meshStart+i) /*matIds[i]*/<<".rawc";
    meshes[i].saveMesh(s.str());
  }

  cout<<"All done."<<endl;
}
