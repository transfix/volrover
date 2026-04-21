/******************************************************************************
				Copyright   

This code is developed within the Computational Visualization Center at The 
University of Texas at Austin.

This code has been made available to you under the auspices of a Lesser General 
Public License (LGPL) (http://www.ices.utexas.edu/cvc/software/license.html) 
and terms that you have agreed to.

Upon accepting the LGPL, we request you agree to acknowledge the use of use of 
the code that results in any published work, including scientific papers, 
films, and videotapes by citing the following references:

C. Bajaj, Z. Yu, M. Auer
Volumetric Feature Extraction and Visualization of Tomographic Molecular Imaging
Journal of Structural Biology, Volume 144, Issues 1-2, October 2003, Pages 
132-143.

If you desire to use this code for a profit venture, or if you do not wish to 
accept LGPL, but desire usage of this code, please contact Chandrajit Bajaj 
(bajaj@ices.utexas.edu) at the Computational Visualization Center at The 
University of Texas at Austin for a different license.
******************************************************************************/

#include <iostream>
#include <string>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include <xmlrpc/XmlRpc.h>

#include <Segmentation/GenSeg/genseg.h>
#include <Segmentation/SegCapsid/segcapsid.h>
#include <Segmentation/SegMed/segmed.h>
#include <Segmentation/SegMonomer/segmonomer.h>
#include <Segmentation/SegSubunit/segsubunit.h>
#include <Segmentation/SecStruct/secstruct.h>

using namespace XmlRpc;
using namespace std;

XmlRpcServer s;

class GenSegmentation : public XmlRpcServerMethod
{
public:
  GenSegmentation(XmlRpcServer* s) : XmlRpcServerMethod("GenSegmentation", s) {}
  
  void execute(XmlRpcValue& params, XmlRpcValue& result)
  {
    if(!generalSegmentation(params,result)) cout << "Error running general segmentation." << endl;
    else cout << "General segmentation complete." << endl;
  }
 
  string help() { return string("Segment volume data, producing a set of subunit volumes."); } 
} GenSegmentation(&s);


class SegmentCapsid : public XmlRpcServerMethod
{
public:
  SegmentCapsid(XmlRpcServer* s) : XmlRpcServerMethod("SegmentCapsid", s) {}
	
  void execute(XmlRpcValue& params, XmlRpcValue& result)
  {
    if(!virusSegCapsid(params,result)) cout << "Error segmenting capsid." << endl;
    else cout << "Capsid Segmentation complete." << endl;
  }
	
  string help() { return string("Segment virus capsid layer."); }
} SegmentCapsid(&s);

class SegmentMedical : public XmlRpcServerMethod
{
public:
  SegmentMedical(XmlRpcServer* s) : XmlRpcServerMethod("SegmentMedical", s) {}

  void execute(XmlRpcValue& params, XmlRpcValue& result)
  {
    if(!medicalSegmentation(params,result)) cout << "Error running medical segmentation." << endl;
    else cout << "Medical segmentation complete." << endl;
  }

  string help() { return string("Segmentation for medical data."); }
} SegmentMedical(&s);

class SegmentMonomer : public XmlRpcServerMethod
{
public:
  SegmentMonomer(XmlRpcServer* s) : XmlRpcServerMethod("SegmentMonomer", s) {}
  
  void execute(XmlRpcValue& params, XmlRpcValue& result)
  {
    if(!virusSegMonomer(params, result)) cout << "Error segmenting monomer." << endl;
    else cout << "Monomer segmentation complete." << endl;
  }
  
  string help() { return string("Segment monomer."); }
} SegmentMonomer(&s);

class SegmentSubunit : public XmlRpcServerMethod
{
public:
  SegmentSubunit(XmlRpcServer* s) : XmlRpcServerMethod("SegmentSubunit", s) {}

  void execute(XmlRpcValue& params, XmlRpcValue& result)
  {
    if(!virusSegSubunit(params, result)) cout << "Error segmenting subunit." << endl;
    else cout << "Subunit segmentation complete." << endl;
  }
  
  string help() { return string("Segment subunit."); }
} SegmentSubunit(&s);

class SecondaryStructureDetection : public XmlRpcServerMethod
{
public:
  SecondaryStructureDetection(XmlRpcServer* s) : XmlRpcServerMethod("SecondaryStructureDetection", s) {}
  
  void execute(XmlRpcValue& params, XmlRpcValue& result)
  {
    if(!secondaryStructureDetection(params, result)) cout << "Error running secondary structure detection." << endl;
    else cout << "Secondary structure detection complete." << endl;
  }
  
  string help() { return string("Secondary Structure Detection!"); }
};

int main(int argc, char **argv)
{
  if(argc != 2)
    {
      cerr<< "Usage: " << argv[0] << " serverPort\n";
      return -1;
    }
  
  int port = atoi(argv[1]);
  XmlRpc::setVerbosity(5);
  s.bindAndListen(port);
  s.enableIntrospection(true);
  s.work(-1.0);
  
  return 0;
}
