//
// Created by dls on 23-3-15.
//

#ifndef DGGRID_CALDGG_H
#define DGGRID_CALDGG_H

#include "dglib/DgIDGG.h"
#include <dglib/DgIDGGutil.h>
#include <memory>
#include <math.h>


class calDgg {
public:

    calDgg( const DgIDGGBase &idgg);

    static bool useearthRadius; //统一规定是否使用地球半径 涉及面积、周长、边长等函数
    const DgIDGGBase& idgg (void) const { return IDGG_; }
    // 计算面积
    long double calarea(const DgQ2DICoord &add, int ptsPerEdgeDensify = 0) const;
    // 计算周长
    long double calper(const DgQ2DICoord& add,int ptsPerEdgeDensify=0) const;
    //计算紧致度
    long double calzsc(const DgQ2DICoord &add, int ptsPerEdgeDensify = 0) const;
    //计算邻近格点最大最小边比 和距离均值 先返回最大最小比 再返回均值
    vector<long double>  calmaxmintotaldis(const DgQ2DICoord& add) const;

    //计算邻近格点最大最小角度比 和角度均值 先返回最大最小比 再返回均值 单位弧度制
    vector<long double>  calmaxmintotalangle(const DgQ2DICoord& add) const;

    /* 计算交点和中点的距离最大最小比（格网边和对对偶边中点距离 ）和 距离均值 先返回最大最小比 再返回均值 单位弧度制    
    计算的是两中点重合度 A comparison of intercell metrics on discrete global grid systems fig2
    */
    vector<long double>  calmaxmintotalcsd(const DgQ2DICoord& add) const;
   
    // 得到多边形顶点(经纬度坐标)
    void setAddVertices (const DgQ2DICoord& add,
                                   DgPolygon& vec, int densify=0) const;
 
    DgGeoCoord    getGeocoord(const DgQ2DICoord& add)const;
    DgQ2DICoord   getQ2DI(const DgGeoCoord& add)const;

    long double caldistance( const DgQ2DICoord& add,const DgQ2DICoord& add1 ) const;
    long double caldistance( const DgGeoCoord& add,const DgGeoCoord& add1 ) const;


    // 计算 编码代表的经纬度中心连线的中点坐标
    DgGeoCoord calmidpoint(const DgQ2DICoord& add,const DgQ2DICoord& add1 ) const;
    DgGeoCoord calmidpoint(const DgGeoCoord& add,const DgGeoCoord& add1 ) const;

    // 计算交点
    DgGeoCoord calintersect(const GeoCoord& sv11, const GeoCoord& sv12,
                     const GeoCoord& sv21, const GeoCoord& sv22, int sign=1) const;

    DgGeoCoord calintersect(const DgGeoCoord& sv11, const DgGeoCoord& sv12,
                     const DgGeoCoord& sv21, const DgGeoCoord& sv22, int sign=1) const;
    // 计算临近 后续添加所有格网的临近 这里只实现了菱形的八邻近
    void calnei (const DgQ2DICoord& add,
                                   DgLocVector& vec) const;
    // 计算以当前编码为中心的两个临近编码练成的二面角
//    long double  calangle (const DgQ2DICoord& center,
//                                   const DgQ2DICoord& add1,const DgQ2DICoord& add2) const;
    // 计算周长  long double calarea()const;
    long double  calangle (const DgQ2DICoord& center,
                           const DgQ2DICoord& add1,const DgQ2DICoord& add2) const;
    // 四邻域在求中心点时需要
    vector<DgQ2DICoord> edgecell(const DgQ2DICoord& add1 ) const;

    // 生成菱形四邻域和八邻域
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
   // 得到多边形面积 不带球面半径的面积
    long double geoPolyArea(const DgPolygon &poly, const DgGeoCoord& center) const;

    // 得到多边形周长 不带球面半径的面积 点要顺序存储（感觉库里应该是顺序存储的 没有查验）
    long double geoPolyper(const DgPolygon &poly ) const;
    bool isEqual(const long double a, const long double b) const  {
        return fabs(a - b) < std::numeric_limits<double>::epsilon();
    }
    const DgIDGGBase& IDGG_;

};
// 初始化静态变量
bool calDgg::useearthRadius = true;

#endif //DGGRID_CALDGG_H
