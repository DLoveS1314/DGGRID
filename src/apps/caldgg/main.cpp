#include <iostream>

using namespace std;

#include <dglib/DgIDGGS4H.h>
#include <dglib/DgIDGGS4D.h>
#include <dglib/DgIDGG.h>
#include "dglib/DgBoundedIDGG.h"
#include "dglib/DgBoundedIDGGS.h"

void testIDGGbounds()
{
    DgRFNetwork net0;
    const DgGeoSphRF& geoRF = *(DgGeoSphRF::makeRF(net0, "GS0"));
     DgGeoCoord vert0(0.0L, 90.0L, false); // args: lon, lat, isRadians

    long double azimuth = 0.0L;
    const DgIDGGS4D * idggsDPtr = DgIDGGS4D::makeRF(net0, geoRF, vert0, azimuth, 10,"ISEA4D_D4","ISEA",dgg::topo::D4);
    const DgIDGGS4D& idggs = *idggsDPtr;
    const DgBoundedIDGGS dgBoundedIdggs = DgBoundedIDGGS(*idggsDPtr);
    DgResAdd<DgQ2DICoord>  add = DgResAdd ( DgQ2DICoord (1, DgIVec2D(3,5)) ,  3);
    unsigned long long int seqnum = dgBoundedIdggs.seqNumAddress(add);
    cout<<"seqnum " << seqnum<<endl;
}
int main (int, char**)
{
    testIDGGbounds();
}

