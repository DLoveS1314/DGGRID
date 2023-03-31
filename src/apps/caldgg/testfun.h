
#ifndef TEST_FUN_H
#define TEST_FUN_H

#include "dglib/DgIDGG.h"
#include <dglib/DgIDGGutil.h>
#include <memory>
#include <iostream>

using namespace std;

#include <dglib/DgIDGGS4H.h>
#include <dglib/DgIDGGS4D.h>
#include <dglib/DgIDGG.h>
#include "dglib/DgBoundedIDGG.h"
#include "dglib/DgBoundedIDGGS.h"
#include "calDgg.h"
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

void testcalArea()
{
    DgRFNetwork net0;
    const DgGeoSphRF& geoRF = *(DgGeoSphRF::makeRF(net0, "GS0"));
    DgGeoCoord vert0(0.0L, 90.0L, false); // args: lon, lat, isRadians

    long double azimuth = 0.0L;
    const DgIDGGS4D * idggsDPtr = DgIDGGS4D::makeRF(net0, geoRF, vert0, azimuth, 10,"ISEA4D_D4","ISEA",dgg::topo::D4);
    const DgIDGGS4D& idggs = *idggsDPtr;
    const DgIDGGBase& dgg = idggs.idgg(6);
//    static_cast<const DgIDGG&>(idgg);
    calDgg cdgg = calDgg(dgg);
    DgQ2DICoord dggadress = DgQ2DICoord(2,DgIVec2D(10,10));
    cdgg.calarea(dggadress);
//    const DgBoundedIDGGS dgBoundedIdggs = DgBoundedIDGGS(*idggsDPtr);
//    DgResAdd<DgQ2DICoord>  add = DgResAdd ( DgQ2DICoord (1, DgIVec2D(3,5)) ,  3);
//    unsigned long long int seqnum = dgBoundedIdggs.seqNumAddress(add);
//    cout<<"seqnum " << seqnum<<endl;
}
void testsphTriSolve()
{
    cout<<"测试sphTriSolve 的面积和geoTriArea 结果是否相同"<<endl;
    // 结果是相同的 但是使用了不同的球面三角形面积算法
    SphTri tri;
    // long double M_PI_180 =  0.0174532925199432957692369076848861271111L;
    // long double M_180_PI =  57.29577951308232087679815481410517033240547L;
//    long double  lat1 = 20.5*M_PI_180;
    struct GeoCoord a={lat:20.5*M_PI_180,lon:40.5*M_PI_180};
    struct GeoCoord b={lat:21.5*M_PI_180,lon:41.5*M_PI_180};
    struct GeoCoord c={lat:19.5*M_PI_180,lon:42.5*M_PI_180};

    tri.verts[0]= a;
    tri.verts[1]= b;
    tri.verts[2]= c;
    printSphTri(tri);
    sphTriSolve(&tri);
    printSphTri(tri);
    DgGeoCoord v1 = DgGeoCoord(a.lon,a.lat);
    DgGeoCoord v2 = DgGeoCoord(b.lon,b.lat);
    DgGeoCoord v3 = DgGeoCoord(c.lon,c.lat);
 long double area =   DgGeoCoord::geoTriArea(v1, v2, v3);
 std::cout << "* area  " << area*DgGeoSphRF::earthRadiusKM()*DgGeoSphRF::earthRadiusKM() << std::endl;
}

void testcalange ()
{
    /*比较种生成角度的方式 结果是否相同*/
    DgRFNetwork net0;
    const DgGeoSphRF& geoRF = *(DgGeoSphRF::makeRF(net0, "GS0"));
    DgGeoCoord vert0(0.0L, 90.0L, false); // args: lon, lat, isRadians

    long double azimuth = 0.0L;
    const DgIDGGS4D * idggsDPtr = DgIDGGS4D::makeRF(net0, geoRF, vert0, azimuth, 10,"ISEA4D_D4","ISEA",dgg::topo::D4);
    const DgIDGGS4D& idggs = *idggsDPtr;
    const DgIDGGBase& dgg = idggs.idgg(6);
//    static_cast<const DgIDGG&>(idgg);
    calDgg cdgg = calDgg(dgg);
    DgQ2DICoord dggadress = DgQ2DICoord(2,DgIVec2D(10,10));
    vector<DgQ2DICoord> dggadress_nei = cdgg.dmdneicell(dggadress);
    // DgPolygon poly;
    // cdgg.setAddVertices(dggadress,poly);
    cout<< "complete"<<endl;
    // cout<< poly<<endl;

    cdgg.calangle(dggadress,dggadress_nei[0],dggadress_nei[1]);
//    const DgBoundedIDGGS dgBoundedIdggs = DgBoundedIDGGS(*idggsDPtr);
//    DgResAdd<DgQ2DICoord>  add = DgResAdd ( DgQ2DICoord (1, DgIVec2D(3,5)) ,  3);
//    unsigned long long int seqnum = dgBoundedIdggs.seqNumAddress(add);
//    cout<<"seqnum " << seqnum<<endl;
}
#endif //DGGRID_CALDGG_H
