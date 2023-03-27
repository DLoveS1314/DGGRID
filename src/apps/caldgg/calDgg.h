//
// Created by dls on 23-3-15.
//

#ifndef DGGRID_CALDGG_H
#define DGGRID_CALDGG_H

#include "dglib/DgIDGG.h"
#include <dglib/DgIDGGutil.h>
#include <memory>
class calDgg {
public:

    calDgg( const DgIDGGBase &idgg);


    const DgIDGGBase& idgg (void) const { return IDGG_; }
    // 计算面积
    long double calarea(DgQ2DICoord& add,int ptsPerEdgeDensify=0)const;
    // 得到多边形顶点(经纬度坐标)
    void setAddVertices (const DgQ2DICoord& add,
                                   DgPolygon& vec, int densify=0) const;
 
    DgGeoCoord    getGeocoord(const DgQ2DICoord& add)const;
    DgQ2DICoord    getQ2DI(const DgGeoCoord& add)const;

    long double caldistance( const DgQ2DICoord& add,const DgQ2DICoord& add1,bool useearthRadius=true) const;
    DgGeoCoord calmidpoint(const DgQ2DICoord& add,const DgQ2DICoord& add1 ) const;
    // 计算交点
    DgGeoCoord calintersect(const GeoCoord& sv11, const GeoCoord& sv12,
                     const GeoCoord& sv21, const GeoCoord& sv22, int sign) const;
    
    // 计算临近 后续添加所有格网的临近 这里只实现了菱形的八邻近
    void calnei (const DgQ2DICoord& add,
                                   DgLocVector& vec) const;
    // 计算以当前编码为中心的两个临近编码练成的二面角
//    long double  calangle (const DgQ2DICoord& center,
//                                   const DgQ2DICoord& add1,const DgQ2DICoord& add2) const;
    
    long double  calangle (const DgQ2DICoord& center,
                           const DgQ2DICoord& add1,const DgQ2DICoord& add2) const;
    //    生成菱形四邻域和八邻域
    vector<DgQ2DICoord> dmdfourneicell(DgQ2DICoord& add1 ) const;
    vector<DgQ2DICoord> dmdneicell(const DgQ2DICoord& add1 ) const;

    long double maxval(const long double val1, const long double val2) const
    {
    /*
        return the maxmum of two variables
    */
        long double maxx;
        if (val1>val2) maxx=val1;
        else maxx=val2;
        return maxx;
    } /* long double maxval */
    long double minval(const long double val1, const long double val2) const
    {
    /*
    return the minmum of two variables
    */
        long double minn;
        if (val1<val2) minn=val1;
        else minn=val2;
        return minn;
    } /* long double minval */
    //   DgGeoCoord (const DgDVec2D& coord, bool rads = true)
    //     { if (rads) *this = coord; else *this = coord * M_PI_180; }
 

private:
   // 得到多边形面积
    long double geoPolyArea(const DgPolygon &poly, const DgGeoCoord& center) const;
    const DgIDGGBase& IDGG_;

};


#endif //DGGRID_CALDGG_H
