#include <iostream>

using namespace std;
#include <dglib/DgBoundedIDGG.h>
#include <dglib/DgCell.h>
#include <dglib/DgIDGG.h>
#include <dglib/DgIDGGBase.h>
#include <dglib/DgIDGGSBase.h>
#include <dglib/DgInputStream.h>
#include <dglib/DgOutAIGenFile.h>
#include <dglib/DgOutputStream.h>
#include <memory>
using namespace std;
////////////////////////////////////////////////////////////////////////////////
void doTransform (string inFileName,string outFileName,char inputDelimiter, 
char outputDelimiter,const DgIDGGBase& dgg ,string inAddType = "Q2DI",string outAddType = "SEQNUM",bool isApSeq =false,int precision=7
)
{

   // set-up to convert to degrees
   DgGeoSphDegRF::makeRF(dgg.geoRF(), dgg.geoRF().name() + "Deg");
   // set-up the input reference frame
   const char* calValstr = "test remainder,test ,";
   bool inSeqNum = false;
   const DgRFBase* pInRF = NULL;
   if (inAddType == "GEO") pInRF = &dgg.geoRF();
   else if (inAddType == "PROJTRI") pInRF = &dgg.projTriRF();
   else if (inAddType == "VERTEX2DD") pInRF = &dgg.vertexRF();
   else if (inAddType == "Q2DD") pInRF = &dgg.q2ddRF();
   else if (inAddType == "INTERLEAVE") pInRF = &dgg.intRF();
   else if (inAddType == "PLANE") pInRF = &dgg.planeRF();
   else if (inAddType == "Q2DI") pInRF = &dgg;
   else if (inAddType == "SEQNUM")
   {
      if (isApSeq)
         ::report("input_address_type of SEQNUM not supported for dggs_aperture_type of SEQUENCE",
                  DgBase::Fatal);

      inSeqNum = true;
      pInRF = &dgg;
   }
   const DgRFBase& inRF = *pInRF;

   // set-up the output reference frame

   bool outSeqNum = false;
   const DgRFBase* pOutRF = NULL;
   if (outAddType == "GEO") pOutRF = &dgg.geoRF();
   else if (outAddType == "PROJTRI") pOutRF = &dgg.projTriRF();
   else if (outAddType == "VERTEX2DD") pOutRF = &dgg.vertexRF();
   else if (outAddType == "Q2DD") pOutRF = &dgg.q2ddRF();
   else if (outAddType == "INTERLEAVE") {
      if (isApSeq)
         ::report("output_address_type of INTERLEAVE not supported for dggs_aperture_type of SEQUENCE",
                  DgBase::Fatal);

      pOutRF = &dgg.intRF();
   } else if (outAddType == "PLANE") pOutRF = &dgg.planeRF();
   else if (outAddType == "Q2DI") pOutRF = &dgg;
   else if (outAddType == "SEQNUM")
   {
      outSeqNum = true;
      pOutRF = &dgg;
   }
   const DgRFBase& outRF = *pOutRF;

   // set the precision

   const_cast<DgRFBase&>(outRF).setPrecision(precision);

   // now process the addresses in the input file

   const int maxLine = 1000;
   char buff[maxLine];
   

   dgcout << "transforming values..." << endl;

   DgInputStream inFile(inFileName, "", DgBase::Fatal);
   ofstream* pOutFile;
   pOutFile = new DgOutputStream(outFileName, "", DgBase::Fatal);

   ofstream& outFile = *pOutFile;

   char delimStr[2];
   delimStr[0] = inputDelimiter;
   delimStr[1] = '\0';

   while (1)
   {
      // get the next line

      inFile.getline(buff, maxLine);//一次读取一行
      if (inFile.eof()) break;

      // parse the address

      DgLocation* loc = NULL;
      if (inSeqNum)
      {
         char* snStr;
         // 分割字符串函数
         snStr = strtok(buff, delimStr);
         unsigned long int sNum;
         if (sscanf(snStr, "%lu", &sNum) != 1)
         {
            ::report("doTransform(): invalid SEQNUM " + string(snStr),
                     DgBase::Fatal);
         }
         loc = static_cast<const DgIDGGBase&>(inRF).bndRF().locFromSeqNum(sNum);
      }
      else
      {
         loc = new DgLocation(inRF);
         string snStr1 =string(buff);
         loc->fromString(buff, inputDelimiter);
      }

      // convert the address
      outRF.convert(loc);

      // output the converted line
  
      if (outSeqNum)
      {
         outFile << static_cast<const DgIDGGBase&>(outRF).bndRF().seqNum(*loc);
      }
      else
      {
         outFile << loc->asString(outputDelimiter);
      }
 
      outFile << outputDelimiter << calValstr << endl;
 
      outFile << endl;

      delete loc;

   }

   inFile.close();
   outFile.close();
   delete pOutFile;

} // void doTransform

vector<DgQ2DICoord>  trans2QDI (string inFileName,  const DgIDGGBase& dgg ,string inAddType = "SEQNUM",string outAddType = "Q2DI",bool isApSeq =false,int precision=7
)
{
      // set-up to convert to degrees
   DgGeoSphDegRF::makeRF(dgg.geoRF(), dgg.geoRF().name() + "Deg");
   // set-up the input reference frame
   // const char* calValstr = "test remainder,test ,";
   bool inSeqNum = false;
   const DgRFBase* pInRF = NULL;
   if (inAddType == "GEO") pInRF = &dgg.geoRF();
   else if (inAddType == "PROJTRI") pInRF = &dgg.projTriRF();
   else if (inAddType == "VERTEX2DD") pInRF = &dgg.vertexRF();
   else if (inAddType == "Q2DD") pInRF = &dgg.q2ddRF();
   else if (inAddType == "INTERLEAVE") pInRF = &dgg.intRF();
   else if (inAddType == "PLANE") pInRF = &dgg.planeRF();
   else if (inAddType == "Q2DI") pInRF = &dgg;
   else if (inAddType == "SEQNUM")
   {
      if (isApSeq)
         ::report("input_address_type of SEQNUM not supported for dggs_aperture_type of SEQUENCE",
                  DgBase::Fatal);

      inSeqNum = true;
      pInRF = &dgg;
   }
   const DgRFBase& inRF = *pInRF;

   const DgRFBase& outRF =dgg;

   // set the precision

   const_cast<DgRFBase&>(outRF).setPrecision(precision);

   // now process the addresses in the input file

   const int maxLine = 1000;
   char buff[maxLine];
   

   dgcout << "transforming values..." << endl;

   DgInputStream inFile(inFileName, "", DgBase::Fatal);

   char delimStr[2];
   char inputDelimiter = ' ';
   delimStr[0] =inputDelimiter;
   delimStr[1] = '\0';
   vector<DgQ2DICoord> adds ;
   // 跳过第一行 属于列名
   inFile.getline(buff, maxLine);//一次读取一行
   while (1)
   {
      // get the next line

      inFile.getline(buff, maxLine);//一次读取一行
      if (inFile.eof()) break;

      // parse the address

      DgLocation* loc = NULL;
      if (inSeqNum)
      {
         char* snStr;
         // 分割字符串函数
         snStr = strtok(buff, delimStr);
         unsigned long int sNum;
         if (sscanf(snStr, "%lu", &sNum) != 1)
         {
            ::report("doTransform(): invalid SEQNUM " + string(snStr),
                     DgBase::Fatal);
         }
         loc = static_cast<const DgIDGGBase&>(inRF).bndRF().locFromSeqNum(sNum);
      }
      else
      {
         loc = new DgLocation(inRF);
         string snStr1 =string(buff);
         loc->fromString(buff, inputDelimiter);
      }

      // convert the address
      outRF.convert(loc);
      const DgQ2DICoord add = static_cast<const DgAddress<DgQ2DICoord>*>(loc->address())->address();
      adds.push_back(add);
      delete loc;

   }
   inFile.close();

   return adds;
}


