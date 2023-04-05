#include <iostream>

using namespace std;

#include <dglib/DgIDGGS4H.h>
#include <dglib/DgIDGGS4D.h>
#include <dglib/DgIDGG.h>
#include "dglib/DgBoundedIDGG.h"
#include "dglib/DgBoundedIDGGS.h"
#include "calDgg.h"
#include "testfun.h"
int main (int, char**)
{
//    testIDGGbounds();
    // testcalArea();
    // testcalange();
    // testreadadd();
    // testfileio();
    // testtrans2QDI();
    // testiteratefile();
    testresult2file();
    // printf("end");
    // DgRFNetwork net0;
    // const DgGeoSphRF& geoRF = *(DgGeoSphRF::makeRF(net0, "GS0"));
    // DgGeoCoord vert0(0.0L, 90.0L, false); // args: lon, lat, isRadians

    // long double azimuth = 0.0L;
    // const DgIDGGS4D * idggsDPtr = DgIDGGS4D::makeRF(net0, geoRF, vert0, azimuth, 10,"ISEA4D_D4","ISEA",dgg::topo::D4);
    // const DgIDGGS4D& idggs = *idggsDPtr;
    // const DgIDGGBase& dgg = idggs.idgg(6);


    // dgcerr<<dgg.gridTopo()<< " setAddNeighbors not yet implemented"<<endl;

    // dgcerr<<dgg.gridTopo()<< " setAddNeighbors not yet implemented"<<endl;
}


