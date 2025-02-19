/*******************************************************************************
    Copyright (C) 2021 Kevin Sahr

    This file is part of DGGRID.

    DGGRID is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    DGGRID is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
*******************************************************************************/
////////////////////////////////////////////////////////////////////////////////
//
// binpres.cpp: perform presence/absence binning of point values into grid
//
// Created by Kevin Sahr, October 1, 2001
//    version 3.1b: November 11, 2001
//    version 4.3b: June 21, 2003
//
////////////////////////////////////////////////////////////////////////////////

#include <iostream>

using namespace std;

#include "dggrid.h"
#include <dglib/DgIDGGBase.h>
#include <dglib/DgIDGGSBase.h>
#include <dglib/DgBoundedIDGG.h>
#include <dglib/DgInputStream.h>

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
BinPresenceParam::BinPresenceParam (DgParamList& plist)
      : MainParam(plist), wholeEarth (true), outFile (0), outSeqNum (false),
        inputDelimiter (' '), outputDelimiter (' '), outputAllCells (true)
{
      /////// fill state variables from the parameter list //////////

      string dummy;
      getParamValue(plist, "bin_coverage", dummy, false);
      if (dummy == "GLOBAL")
         wholeEarth = true;
      else
         wholeEarth = false;

      //// input file names

      string inFileStr;
      getParamValue(plist, "input_files", inFileStr, false);

      char* names = new char[inFileStr.length() + 1];
      inFileStr.copy(names, string::npos);
      names[inFileStr.length()] = 0;

      char* name = strtok(names, " ");
      while (name != NULL)
      {
         inputFiles.push_back(string(name));
         name = strtok(NULL, " ");
      }
      delete [] names;

      // input delimiter

      getParamValue(plist, "input_delimiter", dummy, false);
      if (dummy.length() != 3 || dummy.c_str()[0] != '"' ||
          dummy.c_str()[2] != '"')
      {
         ::report(
          "invalid input_delimiter; must be a single char in double quotes",
          DgBase::Fatal);
      }
      inputDelimiter = dummy.c_str()[1];

      // output file name

      getParamValue(plist, "output_file_name", outFileNameBase, false);

      // output address type

      getParamValue(plist, "output_address_type", outAddType, false);

      // output delimiter

      getParamValue(plist, "output_delimiter", dummy, false);
      if (dummy.length() != 3 || dummy.c_str()[0] != '"' ||
          dummy.c_str()[2] != '"')
      {
         ::report(
          "invalid output_delimiter; must be a single char in double quotes",
          DgBase::Fatal);
      }
      outputDelimiter = dummy.c_str()[1];

      // cell_output_control

      getParamValue(plist, "cell_output_control", dummy, false);
      if (dummy == "OUTPUT_ALL") outputAllCells = true;
      else outputAllCells = false;

      // output counts?

      getParamValue(plist, "output_count", outputCount, false);

} // BinPresenceParam::BinPresenceParam

////////////////////////////////////////////////////////////////////////////////
BinPresenceParam::~BinPresenceParam () {}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
void BinPresenceParam::dump (void)
{
   MainParam::dump();

   dgcout << "BEGIN BINPRESENCE PARAMETER DUMP" << endl;

   dgcout << " wholeEarth: " << wholeEarth << endl;
   dgcout << " outFileNameBase: " << outFileNameBase << endl;
   dgcout << " outFileName: " << outFileName << endl;

   dgcout << " *outFile: ";
   if (outFile)
      dgcout << "(allocated)" << endl;
   else
      dgcout << "null" << endl;

   dgcout << " inputFiles: " << endl;
   for (unsigned int i = 0; i < inputFiles.size(); i++)
      dgcout << "  " << i << " " << inputFiles[i] << endl;

   dgcout << " outAddType: " << outAddType << endl;
   dgcout << " outSeqNum: " << outSeqNum << endl;
   dgcout << " inputDelimiter: " << inputDelimiter << endl;
   dgcout << " outputDelimiter: " << outputDelimiter << endl;
   dgcout << " inFormatStr: " << inFormatStr << endl;
   dgcout << " outputAllCells: " << outputAllCells << endl;
   dgcout << " outputCount: " << ((outputCount) ? "true" : "false") << endl;

   dgcout << "END BINPRESENCE PARAMETER DUMP" << endl;

} // void BinPresenceParam::dump

////////////////////////////////////////////////////////////////////////////////
void binPresGlobal (BinPresenceParam& dp)
{
   ////// create the reference frames ////////

   DgRFNetwork net0;
   const DgGeoSphRF& geoRF = *(DgGeoSphRF::makeRF(net0, dp.datum, dp.earthRadius));
   const DgIDGGSBase *idggs = DgIDGGSBase::makeRF(net0, geoRF, dp.vert0,
             dp.azimuthDegs, dp.aperture, dp.actualRes+1, dp.gridTopo,
             dp.gridMetric, "IDGGS", dp.projType, dp.isMixed43, dp.numAp4,
             dp.isSuperfund, dp.isApSeq, dp.apSeq);
   const DgIDGGBase& dgg = idggs->idggBase(dp.actualRes);

   dgcout << "Res " << dgg.outputRes() << " " << dgg.gridStats() << endl;

   // set-up to convert to degrees
   DgGeoSphDegRF::makeRF(geoRF, geoRF.name() + "Deg");

   // set-up the output reference frame

   dp.outSeqNum = false;
   const DgRFBase* pOutRF = NULL;
   if (dp.outAddType == "PROJTRI") pOutRF = &dgg.projTriRF();
   else if (dp.outAddType == "VERTEX2DD") pOutRF = &dgg.vertexRF();
   else if (dp.outAddType == "Q2DD") pOutRF = &dgg.q2ddRF();
   else if (dp.outAddType == "INTERLEAVE") pOutRF = &dgg.intRF();
   else if (dp.outAddType == "PLANE") pOutRF = &dgg.planeRF();
   else if (dp.outAddType == "Q2DI") pOutRF = &dgg;
   else if (dp.outAddType == "SEQNUM")
   {
      dp.outSeqNum = true;
      pOutRF = &dgg;
   }
   else
   {
      ::report("binPresGlobal(): invalid output_address_type " +
               dp.outAddType, DgBase::Fatal);
   }

   const DgRFBase& outRF = *pOutRF;

   // create an array to store the values

   int nClasses = (int) dp.inputFiles.size();
   //bool** vals = new (bool (*[dgg.bndRF().size()]));
   bool** vals = new bool*[dgg.bndRF().size()];
   for (unsigned long int i = 0; i < dgg.bndRF().size(); i++)
   {
      vals[i] = new bool[nClasses];
      for (int j = 0; j < nClasses; j++) vals[i][j] = false;
   }

   // now process the points in each input file

   dgcout << "binning points..." << endl;

   const int maxLine = 100;
   char buff[maxLine];
   double lon, lat;
   for (unsigned int fc = 0; fc < dp.inputFiles.size(); fc++)
   {
      DgInputStream inFile(dp.inputFiles[fc], "", DgBase::Fatal);

      while (1)
      {
         inFile.getline(buff, maxLine);
         if (inFile.eof()) break;

         int result;
         result = sscanf(buff, dp.inFormatStr.c_str(), &lon, &lat);
         if (result != 2)
         {
            ::report("binPresGlobal(): invalid format in file " +
                     dp.inputFiles[fc], DgBase::Fatal);
         }

         DgLocation* tloc = geoRF.makeLocation(DgGeoCoord(lon, lat, false));
         dgg.convert(tloc);
         long int sNum = dgg.bndRF().seqNum(*tloc);
         long int ndx = sNum - 1;
         delete tloc;

         vals[ndx][fc] = true;
      }

      inFile.close();
   }

   ///// output the values /////

   for (unsigned long int i = 0; i < dgg.bndRF().size(); i++)
   {
      int count = 0;
      for (int j = 0; j < nClasses; j++) if (vals[i][j]) count++;

      if (!dp.outputAllCells && count == 0) continue;

      unsigned long int sNum = i + 1;
      if (dp.outSeqNum)
         *dp.outFile << sNum << dp.outputDelimiter;
      else
      {
         DgLocation* tloc = dgg.bndRF().locFromSeqNum(sNum);
         outRF.convert(tloc);
         *dp.outFile << tloc->asString(dp.outputDelimiter)
                     << dp.outputDelimiter;
         delete tloc;
      }

      if (dp.outputCount) *dp.outFile << count << dp.outputDelimiter;
      for (int j = 0; j < nClasses; j++)
         *dp.outFile << ((vals[i][j]) ? 1 : 0);
      *dp.outFile << endl;
   }

   ///// clean-up /////

   for (unsigned long int i = 0; i < dgg.bndRF().size(); i++) delete [] vals[i];
   delete [] vals;

} // void binPresGlobal

////////////////////////////////////////////////////////////////////////////////
class QuadVals {

   public:

      bool isUsed;
      DgIVec2D offset;     // offset of min (i, j)
      DgIVec2D upperRight; // (maxi, maxj) relative to offset
      int numI;
      int numJ;

      bool*** vals;

};

////////////////////////////////////////////////////////////////////////////////
void binPresPartial (BinPresenceParam& dp)
{
   ////// create the reference frames ////////
   DgRFNetwork net0;
   const DgGeoSphRF& geoRF = *(DgGeoSphRF::makeRF(net0, dp.datum, dp.earthRadius));
      const DgIDGGSBase *idggs = DgIDGGSBase::makeRF(net0, geoRF, dp.vert0,
             dp.azimuthDegs, dp.aperture, dp.actualRes+1, dp.gridTopo,
             dp.gridMetric, "IDGGS", dp.projType, dp.isMixed43, dp.numAp4,
             dp.isSuperfund, dp.isApSeq, dp.apSeq);
      const DgIDGGBase& dgg = idggs->idggBase(dp.actualRes);

   dgcout << "Res " << dgg.outputRes() << " " << dgg.gridStats() << endl;

   // set-up to convert to degrees
   DgGeoSphDegRF::makeRF(geoRF, geoRF.name() + "Deg");

   // set-up the output reference frame

   dp.outSeqNum = false;
   const DgRFBase* pOutRF = NULL;
   if (dp.outAddType == "PROJTRI") pOutRF = &dgg.projTriRF();
   else if (dp.outAddType == "VERTEX2DD") pOutRF = &dgg.vertexRF();
   else if (dp.outAddType == "Q2DD") pOutRF = &dgg.q2ddRF();
   else if (dp.outAddType == "INTERLEAVE") pOutRF = &dgg.intRF();
   else if (dp.outAddType == "PLANE") pOutRF = &dgg.planeRF();
   else if (dp.outAddType == "Q2DI") pOutRF = &dgg;
   else if (dp.outAddType == "SEQNUM")
   {
      dp.outSeqNum = true;
      pOutRF = &dgg;
   }
   else
   {
      ::report("binPresPartial(): invalid output_address_type " +
               dp.outAddType, DgBase::Fatal);
   }

   const DgRFBase& outRF = *pOutRF;

   // create a place to store the values by quad

   QuadVals qvals[12];
   for (int q = 0; q < 12; q++)
   {
      qvals[q].isUsed = false;
      qvals[q].offset = DgIVec2D(dgg.maxI() + 1, dgg.maxJ() + 1);
      qvals[q].upperRight = DgIVec2D(-1, -1);
      qvals[q].numI = 0;
      qvals[q].numJ = 0;
      qvals[q].vals = 0;
   }

   // now make a first pass through the input files and determine what
   // cells are represented

   dgcout << "determing quad bounds..." << endl;

   const int maxLine = 100;
   char buff[maxLine];
   double lon, lat;
   for (unsigned int fc = 0; fc < dp.inputFiles.size(); fc++)
   {
      DgInputStream inFile(dp.inputFiles[fc], "", DgBase::Fatal);

      while (1)
      {
         inFile.getline(buff, maxLine);
         if (inFile.eof()) break;

         int result;
         result = sscanf(buff, dp.inFormatStr.c_str(), &lon, &lat);
         if (result != 2)
         {
            ::report("binPresPartial(): invalid format in file " +
                      dp.inputFiles[fc], DgBase::Fatal);
         }

         DgLocation* tloc = geoRF.makeLocation(DgGeoCoord(lon, lat, false));
         dgg.convert(tloc);

         int q = dgg.getAddress(*tloc)->quadNum();
         const DgIVec2D& coord = dgg.getAddress(*tloc)->coord();
         QuadVals& qv = qvals[q];

         qv.isUsed = true;
         if (coord.i() < qv.offset.i()) qv.offset.setI(coord.i());
         if (coord.i() > qv.upperRight.i()) qv.upperRight.setI(coord.i());
         if (coord.j() < qv.offset.j()) qv.offset.setJ(coord.j());
         if (coord.j() > qv.upperRight.j()) qv.upperRight.setJ(coord.j());

         delete tloc;
      }

      inFile.close();
   }

   // now initialize the vals storage in the quads which are used

   int nClasses = (int) dp.inputFiles.size();
   for (int q = 0; q < 12; q++)
   {
      QuadVals& qv = qvals[q];
      if (!qv.isUsed) continue;

      qv.upperRight = qv.upperRight - qv.offset; // make relative

      qv.numI = (int) qv.upperRight.i() + 1;
      qv.numJ = (int) qv.upperRight.j() + 1;
      qv.vals = new bool**[qv.numI];
      for (int i = 0; i < qv.numI; i++)
      {
         qv.vals[i] = new bool*[qv.numJ];

         for (int j = 0; j < qv.numJ; j++)
         {
            qv.vals[i][j] = new bool[nClasses];
            for (int k = 0; k < nClasses; k++) qv.vals[i][j][k] = false;
         }
      }
   }

   // now process the points in each input file

   dgcout << "binning points..." << endl;

   for (unsigned int fc = 0; fc < dp.inputFiles.size(); fc++)
   {
      DgInputStream inFile(dp.inputFiles[fc], "", DgBase::Fatal);

      while (1)
      {
         inFile.getline(buff, maxLine);
         if (inFile.eof()) break;

         int result;
         result = sscanf(buff, dp.inFormatStr.c_str(), &lon, &lat);
         if (result != 2)
         {
            ::report("binPresPartial(): invalid format in file " +
                      dp.inputFiles[fc], DgBase::Fatal);
         }

         DgLocation* tloc = geoRF.makeLocation(DgGeoCoord(lon, lat, false));
         dgg.convert(tloc);

         int q = dgg.getAddress(*tloc)->quadNum();
         QuadVals& qv = qvals[q];
         DgIVec2D coord = dgg.getAddress(*tloc)->coord() - qv.offset;
         delete tloc;

         qv.vals[coord.i()][coord.j()][fc] = true;
      }

      inFile.close();
   }

   ///// output the values /////

   if (dp.outputAllCells)
   {
      int count;
      bool* dummy = new bool[nClasses];
      for (int i = 0; i < nClasses; i++) dummy[i] = false;

      bool* vvec;
      for (unsigned long int i = 0; i < dgg.bndRF().size(); i++)
      {
         unsigned long int sNum = i + 1;
         DgLocation* tloc = dgg.bndRF().locFromSeqNum(sNum);

         // check to see if there is a value for this cell
         count = 0;
         vvec = dummy;
         int q = dgg.getAddress(*tloc)->quadNum();
         QuadVals& qv = qvals[q];
         if (qv.isUsed)
         {
            DgIVec2D coord = dgg.getAddress(*tloc)->coord() - qv.offset;
            if (coord.i() >= 0 && coord.j() >= 0 &&
                coord.i() <= qv.upperRight.i() &&
                coord.j() <= qv.upperRight.j())
            {
              vvec = qv.vals[coord.i()][coord.j()];
              for (int k = 0; k < nClasses; k++) if (vvec[k]) count++;
            }
         }

         // output the value

         if (dp.outSeqNum)
            *dp.outFile << sNum << dp.outputDelimiter;
         else
         {
            DgLocation* tloc = dgg.bndRF().locFromSeqNum(sNum);
            outRF.convert(tloc);
            *dp.outFile << tloc->asString(dp.outputDelimiter)
                        << dp.outputDelimiter;
            delete tloc;
         }

         if (dp.outputCount) *dp.outFile << count << dp.outputDelimiter;
         for (int k = 0; k < nClasses; k++)
            *dp.outFile << ((vvec[k]) ? 1 : 0);
         *dp.outFile << endl;

         delete tloc;
      }
      delete [] dummy;
   }
   else
   {
      for (int q = 0; q < 12; q++)
      {
         QuadVals& qv = qvals[q];
         if (!qv.isUsed) continue;

         for (int i = 0; i < qv.numI; i++)
         {
            for (int j = 0; j < qv.numJ; j++)
            {
               int count = 0;
               for (int k = 0; k < nClasses; k++) if (qv.vals[i][j][k]) count++;
               if (count == 0) continue;

               DgIVec2D coord(qv.offset.i() + i, qv.offset.j() + j);
               DgLocation* tloc = dgg.makeLocation(DgQ2DICoord(q, coord));

               if (dp.outSeqNum)
               {
                  unsigned long int sNum = dgg.bndRF().seqNum(*tloc);
                  *dp.outFile << sNum << dp.outputDelimiter;
               }
               else
               {
                  outRF.convert(tloc);
                  *dp.outFile << tloc->asString(dp.outputDelimiter)
                              << dp.outputDelimiter;
               }
               delete tloc;

               if (dp.outputCount) *dp.outFile << count << dp.outputDelimiter;
               for (int k = 0; k < nClasses; k++)
                  *dp.outFile << ((qv.vals[i][j][k]) ? 1 : 0);
               *dp.outFile << endl;
            }
         }
      }
   }


   ///// clean-up /////

   for (int q = 0; q < 12; q++)
   {
      QuadVals& qv = qvals[q];
      if (!qv.isUsed) continue;

      for (int i = 0; i < qv.numI; i++)
      {
         for (int j = 0; j < qv.numJ; j++) delete [] qv.vals[i][j];
         delete [] qv.vals[i];
      }

      delete [] qv.vals;
   }

} // void binPresPartial

////////////////////////////////////////////////////////////////////////////////
void doBinPresence (BinPresenceParam& dp, DgGridPList& plist)
{
   char tmpStr[DgRFBase::maxFmtStr];
   snprintf(tmpStr, DgRFBase::maxFmtStr, "%%lf%c%%lf%c", dp.inputDelimiter, dp.inputDelimiter);
   dp.inFormatStr = tmpStr;

   for (dp.curGrid = 1; dp.curGrid <= dp.numGrids; dp.curGrid++)
   {
      // first get the grid placement

      dp.outFileName = dp.outFileNameBase;
      dp.metaOutFileName = dp.metaOutFileNameBase;

      orientGrid(dp, plist);

      if (dp.numGrids > 1)
      {
         string suffix = string(".") + dgg::util::to_string(dp.curGrid, 4);
         dp.metaOutFileName = dp.metaOutFileName + suffix;
         dp.outFileName = dp.outFileName + suffix;
      }

      /////// open the output file as applicable //////

      dp.outFile = new ofstream();
      dp.outFile->open(dp.outFileName.c_str());
      dp.outFile->setf(ios::fixed, ios::floatfield);
      dp.outFile->precision(dp.precision);

      if (dp.numGrids > 1 || dp.placeRandom)
      {
         ofstream metaOutFile;
         metaOutFile.open(dp.metaOutFileName.c_str());
         metaOutFile.setf(ios::fixed, ios::floatfield);
         metaOutFile.precision(12);
         metaOutFile << plist;
         metaOutFile.close();
      }

      // now do it

      if (dp.curGrid == dp.numGrids) dp.lastGrid = true;

      if (dp.wholeEarth) binPresGlobal(dp);
      else binPresPartial(dp);

      dgcout << endl;

      // close the output files

      dp.outFile->close();
      delete dp.outFile;
   }

   dgcout << "** presence binning complete **" << endl;

} // void doBinPres

////////////////////////////////////////////////////////////////////////////////
