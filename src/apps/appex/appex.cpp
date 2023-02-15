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
// appex.cpp: simple application that demonstrates using the dglib library to
//            manipulate DGG cells.
//
////////////////////////////////////////////////////////////////////////////////

#include <iostream>

using namespace std;

#include <dglib/DgIDGGS4H.h>
#include <dglib/DgIDGGS4D.h>
#include <dglib/DgIDGG.h>
#include "dglib/DgBoundedIDGG.h"

////////////////////////////////////////////////////////////////////////////////
int main (int, char**)
{
   ///// create the DGG /////

   // create the reference frame (RF) conversion network
   DgRFNetwork net0;

   // create the geodetic reference frame
   // reference frames must be created dynamically using makeRF
   // they will be deleted by the Network
   const DgGeoSphRF& geoRF = *(DgGeoSphRF::makeRF(net0, "GS0"));

   // create the ISEA4H grid system with resolutions 0-9; requires a
   // fixed icosahedron vertex and edge azimuth
   // DgGeoCoord vert0(11.25L, 58.28252559L, false); // args: lon, lat, isRadians
   DgGeoCoord vert0(0.0L, 90.0L, false); // args: lon, lat, isRadians

   long double azimuth = 0.0L;

   // all DGGS's must be created using a factory makeRF method
   // the DGGS is memory managed by the DgRFNetwork
   const DgIDGGS4D * idggsDPtr = DgIDGGS4D::makeRF(net0, geoRF, vert0, azimuth, 10,"ISEA4D_D8","ISEA",dgg::topo::D8);
//    const DgIDGGS4H* idggsPtr = DgIDGGS4H::makeRF(net0, geoRF, vert0, azimuth, 10);
//   const DgIDGGS4H& idggs = *idggsPtr;
    const DgIDGGS4D& idggs = *idggsDPtr;

   // get the resolution 7 dgg from the dggs
   const DgIDGG& dgg = idggs.idgg(3);
   cout << dgg.gridStats() << endl;

   //////// now use the DGG /////////

   ///// given a point in lon/lat, get the cell it lies within /////

   // first create a DgLocation in geoRF coordinates
   // DgGeoCoord geoAddress(-122.7083, 42.1947, false);
   DgGeoCoord geoAddress(177.5, 26.5, false);
    //保证DgGeoCoord是geoRF模板类的一个才能makelocation
   DgLocation* thePt = geoRF.makeLocation(geoAddress);
   cout << "the point " << *thePt << endl;

    const DgBoundedIDGG  loc1 = static_cast<const DgIDGG&>(dgg).bndRF();
    DgLocation* loc2 =  loc1.locFromSeqNum(2);
    cout << "the loc2 " << *loc2 << endl;

   // converting the point location to the dgg RF determines which cell it's in
   dgg.convert(thePt);
//   cout << "* lies in cell " << *thePt << endl;

   // we can get the cell's vertices, which are defined in geoRF
   DgPolygon verts;
   int ptsPerEdgeDensify = 3;
   dgg.setVertices(*thePt, verts, ptsPerEdgeDensify);

//   cout << "* with densified cell boundary:\n" << verts << endl;
   
   DgLocVector neighbors;

   dgg.setNeighbors(*thePt, neighbors);
//   cout << "*neighbors:" << neighbors << endl;
   // we can get the cell's center point by converting the cell back to geoRF
//   geoRF.convert(thePt);
//   cout << "* and cell center point:" << *thePt << endl;
//
//   // we can extract the coordinates into primitive data types
//   const DgGeoCoord& centCoord = *geoRF.getAddress(*thePt);
//   double latRads = centCoord.lat();
//   double lonRads = centCoord.lon();
//   cout << "* center point lon,lat in radians: "
//        << lonRads << ", " << latRads << endl;
//
//   const DgGeoCoord& firstVert = *geoRF.getAddress(verts[0]);
//   double latDegs = firstVert.latDegs();
//   double lonDegs = firstVert.lonDegs();
//   cout << "* first boundary vertex lon,lat in degrees: "
//        << lonDegs << ", " << latDegs << endl;

   delete thePt;

   return 0;

} // main 

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
