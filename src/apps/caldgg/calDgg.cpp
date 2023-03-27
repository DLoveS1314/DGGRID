//
// Created by dls on 23-3-15.
//

#include "calDgg.h"
#include "dglib/DgEllipsoidRF.h"
long double calDgg::calarea(DgQ2DICoord &add, int ptsPerEdgeDensify) const
{
    // cout<<"res "<<idgg().res()<<endl;
    //        获取边界点
    DgPolygon verts;
    // 此处指针无人删除
    // cout << "* lies in cell " << add << endl;

    // DgLocation * loc = idgg().makeLocation(add);
    shared_ptr<DgLocation> loc(idgg().makeLocation(add) ) ;
    // cout << "* lies in cell " << *loc << endl;
    // idgg().setVertices(*loc, verts, ptsPerEdgeDensify);
    this->setAddVertices(add, verts, ptsPerEdgeDensify);
//    // 此处指针不需要删除 他只是指向了一个adress 里面的内容 相应的类会去删除，
//    // 但是后续要至于直接把这个指针指向null，这个是局部变量可以不指向
    const DgGeoCoord* center  = idgg().geoRF().getAddress(*idgg().geoRF().convert(loc.get()));
    // shared_ptr<DgGeoCoord> center(  );
    //loc最终的结果会被改变
    // DgGeoCoord center = static_cast<const DgAddress<DgGeoCoord>*>(idgg().geoRF().convert(loc)->address())->address();
    // DgGeoCoord center = *(idgg().geoRF().getAddress(*idgg().geoRF().convert(loc)));
    // cout << "* lies in cell " << *loc << endl;

    double area = (this->geoPolyArea(verts, *center))*DgGeoSphRF::earthRadiusKM()*DgGeoSphRF::earthRadiusKM();
    // cout << "* area  " << area*DgGeoSphRF::earthRadiusKM()*DgGeoSphRF::earthRadiusKM() << endl;
    // delete loc;
    return  area;
}
long double calDgg::geoPolyArea (const DgPolygon& poly, const DgGeoCoord& center)  const

{
    //
// returns area of spherical polygon in radians;
//
// assumes DgAddressBase for poly is DgGeoCoord
//
// Assumes poly is a "reasonably" convex polygon (i.e. a cell bouncary) and
// that center is an interior point.
//
    long double totArea = 0.0L;

    const DgGeoSphRF* geoRF = dynamic_cast<const DgGeoSphRF*>(&poly.rf());
    if (geoRF == 0)
        report("DgGeoCoord::geoPolyArea() non-geo polygon", DgBase::Fatal);

    // do each sub-triangle
    for (int i = 0; i < poly.size(); i++)
    {
        const DgGeoCoord& v1 = *geoRF->getAddress(poly[i]);
        const DgGeoCoord& v2 = *geoRF->getAddress(poly[(i + 1) % poly.size()]);

        totArea += DgGeoCoord::geoTriArea(center, v1, v2);
    }

    return totArea;

}

DgGeoCoord  calDgg::getGeocoord(const DgQ2DICoord& add) const
{
    // invQuantify idgg的这个函数暂且没有启用 
    // idgg().invQuantify(add);
    shared_ptr<DgLocation> loc(idgg().makeLocation(add) ) ;
    const DgGeoCoord* center  = idgg().geoRF().getAddress(*idgg().geoRF().convert(loc.get()));
    return * center;
}

DgQ2DICoord calDgg::getQ2DI(const DgGeoCoord &add) const
{
    // Quantify idgg的这个函数暂且没有启用 
    // idgg().Quantify(add);
    shared_ptr<DgLocation> loc(idgg().backFrame() .makeLocation(add) ) ;
    const DgQ2DICoord* qij = idgg().getAddress( *(idgg().convert(loc.get())) ) ;
    return  *qij;
}
// useearthRadius 默认是true
long double calDgg::caldistance(const DgQ2DICoord &add,const  DgQ2DICoord &add1,bool useearthRadius ) const
{
    DgGeoCoord addg = this->getGeocoord(add);
    DgGeoCoord addg1 = this->getGeocoord(add1);
    long double distance =   DgGeoCoord::gcDist(addg, addg1); 
    if (useearthRadius)//如果使用地球半径 那就乘上
      {
        distance = DgGeoSphRF::earthRadiusKM() * distance;
       } 
    return distance;
}

DgGeoCoord calDgg::calmidpoint(const DgQ2DICoord& add,const DgQ2DICoord& add1 ) const
{
    DgGeoCoord addg = this->getGeocoord(add);
    DgGeoCoord addg1 = this->getGeocoord(add1);
    DgGeoCoord midpoint =   DgGeoSphRF::midPoint(addg, addg1); 
 
    return midpoint;
}

DgGeoCoord calDgg::calintersect(const GeoCoord &sv11, const GeoCoord &sv12, const GeoCoord &sv21, const GeoCoord &sv22, int sign) const
{   
  /*
  参考的GCintersect
        sign=1：两个大圆弧段，两端分别为 sv11 和 sv12，
            分别为 sv21 和 sv22
        sign=0: 两个完整的大圆，一次通过 sv11 和 sv12，一次通过
             sv21 和 sv22，返回北半球的交点
  */     
    GeoCoord pt;
    Vec3D nn1,nn2,
    p11,p12,
    p21,p22,
    pp,pp2;
    long double a,b,maxlon,minlon;

    /* calculate the intersect point of two great circle */
    p11=llxyz(sv11);
    p12=llxyz(sv12);
    p21=llxyz(sv21);
    p22=llxyz(sv22);
    nn1.x=p11.y*p12.z-p12.y*p11.z;
    nn1.y=-p11.x*p12.z+p12.x*p11.z;
    nn1.z=p11.x*p12.y-p12.x*p11.y;
    nn2.x=p21.y*p22.z-p22.y*p21.z;
    nn2.y=-p21.x*p22.z+p22.x*p21.z;
    nn2.z=p21.x*p22.y-p22.x*p21.y;
    if ((nn2.z*nn1.y-nn1.z*nn2.y)!= 0.0L) 
    {
    b=(nn1.x*nn2.y-nn2.x*nn1.y)/(nn2.z*nn1.y-nn1.z*nn2.y);
    a=(nn2.x*nn1.z-nn1.x*nn2.z)/(nn1.y*nn2.z-nn2.y*nn1.z);
    pp.x=1/sqrtl(a*a+b*b+1);
    pp.y=a*pp.x;
    pp.z=b*pp.x;
    }
     else if (((nn2.z*nn1.y-nn1.z*nn2.y)==0.0L) &&
            ((nn1.x*nn2.y-nn2.x*nn1.y)==0.0L) && ((nn1.x*nn2.z-nn2.x*nn1.z)==0.0L)) 
    {
        report("Error in GCintersect: the two great circle planes are parallel.\n",
                DgBase::Fatal);
    } 
    else if (((nn2.z*nn1.y-nn1.z*nn2.y)==0.0L) && (nn1.z!=0.0L)) 
    {
        pp.x=0.0L;
        pp.y=1.0L/sqrtl(1+nn1.y*nn1.y/nn1.z/nn1.z);
        pp.z=-nn1.y/nn1.z*pp.y;
    }
     else if (((nn2.z*nn1.y-nn1.z*nn2.y)==0.0L) && (nn2.z!=0.0L)) 
     {
        pp.x=0.0L;
        pp.y=1.0L/sqrtl(1.0L+nn2.y*nn2.y/nn2.z/nn2.z);
        pp.z=-nn2.y/nn2.z*pp.y;
    } 
    else if (((nn2.z*nn1.y-nn1.z*nn2.y)==0.0L) && (nn1.y!=0.0L)) 
    {
        pp.x=0.0L;
        pp.z=1/sqrtl(1+nn1.z*nn1.z/nn1.y/nn1.y);
        pp.y=-nn1.z/nn1.y*pp.z;
    } 
    else if (((nn2.z*nn1.y-nn1.z*nn2.y)==0.0L) && (nn2.y!=0.0L)) 
    {
        pp.x=0.0L;
        pp.z=1.0L/sqrtl(1.0L+nn2.z*nn2.z/nn2.y/nn2.y);
        pp.y=-nn2.z/nn2.y*pp.z;
    }

    if (sign==0) 
    {
        if (pp.z<0.0L) {
            pp.x=0.0L-pp.x;
            pp.y=-pp.y;
            pp.z=-pp.z;
    }
    pt=xyzll(pp);
    return pt;
    } 
    else
     {
        /* judge if the point is on the two great circle segment */

        pt=xyzll(pp);
        maxlon=this->maxval(sv11.lon,sv12.lon);
        minlon=this->minval(sv11.lon,sv12.lon);
        if ((pt.lon<=maxlon) && (pt.lon>=minlon))
        return pt;
        else 
        {
            pp2.x=-pp.x;
            pp2.y=-pp.y;
            pp2.z=-pp.z;
            pt=xyzll(pp2);
            if ((pt.lon<=maxlon) && (pt.lat>=minlon))
                return pt;
            else 
            {
                dgcerr << "Error of GCintersect: the point is not on great circle segment.\n";
                pt.lat=UNDEFVAL; pt.lon=UNDEFVAL;
                return pt;
            }
        }
    }
}

/// @brief Calculate the adjacent grid of idgg, only realize the rhombus
/// @param add Piece number row number column number
/// @param vec saved structure
void calDgg::calnei(const DgQ2DICoord &add, DgLocVector &vec) const
{

    // //重写了 DgDiscRF setAddNeighbors
    //DgLocVector在 DgDiscRF setNeighbors 已经转化为了DGIDGG空间 也就是Q i j
   if (this->idgg().gridTopo()== Diamond)
   {
       if (this->idgg().gridMetric() == D4)//现在不存在D8的菱形 所以只能D4菱形 也就是只能边临近和角临近选择一样
       {
           //八临近
           vector <DgQ2DICoord> coords =  this->dmdneicell(add);
           vec.clearAddress();
           for (int i = 0; i < coords.size(); i++)
               vec.push_back(*(this->idgg().makeLocation(coords[i])));
       }
   }
   else
   {
    dgcerr<<"setAddNeighbors "<<this->idgg().gridTopo()<< " not yet implemented"<<endl;
   }
}

long double calDgg::calangle(const DgQ2DICoord &center, const DgQ2DICoord &add1, const DgQ2DICoord &add2) const
{
    // 改编 自DgGeoCoord::geoTriArea
    DgGeoCoord centerg = this->getGeocoord(center);
    DgGeoCoord addg1 = this->getGeocoord(add1);
    DgGeoCoord addg2 = this->getGeocoord(add2);
      // determine the edges

   long double a = DgGeoCoord::gcDist(addg1, addg2);
   long double b = DgGeoCoord::gcDist(centerg, addg2);
   long double c = DgGeoCoord::gcDist(centerg, addg1);

   // determine the angles using half-angle formulas

   long double s = (a + b + c) / 2.0L;

   long double sinsa = sinl(s - a);
   long double sinsb = sinl(s - b);
   long double sinsc = sinl(s - c);

   long double k = sqrtl(sinsa * sinsb * sinsc / sinl(s));

   long double bigA = 2.0L * atanl(k / sinsa);

    long double l1[2],l2[2],l3[2],edges[2];
    cout<<"夹角"<<bigA<<endl;

	
    l1[0]=centerg.lat(); l1[1]=centerg.lon();
    l2[0]=addg1.lat(); l2[1]=addg1.lon();
    l3[0]=addg2.lat(); l3[1]=addg2.lon();
    edges[0]=acosl(cosl(M_PI_2-l2[0])*cosl(M_PI_2-l3[0])+
            sinl(M_PI_2-l2[0])*sinl(M_PI_2-l3[0])*cosl(l2[1]-l3[1]));
    edges[1]=acosl(cosl(M_PI_2-l1[0])*cosl(M_PI_2-l3[0])+
            sinl(M_PI_2-l1[0])*sinl(M_PI_2-l3[0])*cosl(l1[1]-l3[1]));
    edges[2]=acosl(cosl(M_PI_2-l2[0])*cosl(M_PI_2-l1[0])+
            sinl(M_PI_2-l2[0])*sinl(M_PI_2-l1[0])*cosl(l2[1]-l1[1]));
    long double angles =acosl((cosl(edges[0])-cosl(edges[1])*
               cosl(edges[2]))/(sinl(edges[1])*sinl(edges[2])));
    cout<<"夹角"<<angles<<endl;
   
  
 
}

 

 

vector<DgQ2DICoord> calDgg::dmdfourneicell(DgQ2DICoord &add1) const
{
   return vector<DgQ2DICoord>();
}

vector<DgQ2DICoord> calDgg::dmdneicell(const DgQ2DICoord &add1) const
{
   vector<DgQ2DICoord> neiadd;
    //maxI和MaxJ 对于菱形 是相等的所以下面出现了混用的情况
    long long int maxI =this->idgg().maxI();
    long long int maxJ =this->idgg().maxJ();
	// 点在内部的情况
    if(add1.coord().i() != 0 && add1.coord().i() != maxI && add1.coord().j() != 0 && add1.coord().j() != maxJ)
	{
		neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() + 1, add1.coord().j())));
		neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() + 1, add1.coord().j() + 1)));
		neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i(), add1.coord().j() + 1)));
		neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() - 1, add1.coord().j() + 1)));
		neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() - 1, add1.coord().j())));
		neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() - 1, add1.coord().j() - 1)));
		neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i(), add1.coord().j() - 1)));
		neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() + 1, add1.coord().j() - 1)));
	}
	else if(add1.coord().i() == 0 && add1.coord().j() != 0 && add1.coord().j() != maxJ)//第0行的情况 但是列不是0或maxJ
	{
		if(add1.quadNum() == 1)
		{
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() + 1, add1.coord().j())));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() + 1, add1.coord().j() + 1)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i(), add1.coord().j() + 1)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() + 4, DgIVec2D(maxJ - add1.coord().j() - 1, maxJ)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() + 4, DgIVec2D(maxJ - add1.coord().j(), maxJ)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() + 4, DgIVec2D(maxJ - add1.coord().j() + 1,maxJ)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i(), add1.coord().j() - 1)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() + 1, add1.coord().j() - 1)));
		}
		else if(add1.quadNum() > 1 && add1.quadNum() <= 5)
		{
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() + 1, add1.coord().j())));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() + 1, add1.coord().j() + 1)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i(), add1.coord().j() + 1)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() - 1, DgIVec2D(maxJ - add1.coord().j() - 1, maxJ)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() - 1, DgIVec2D(maxJ - add1.coord().j(), maxJ)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() - 1, DgIVec2D(maxJ - add1.coord().j() + 1,maxJ)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i(), add1.coord().j() - 1)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() + 1, add1.coord().j() - 1)));
		}
		else if(add1.quadNum() >= 6 && add1.quadNum() <= 10)
		{
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() + 1, add1.coord().j())));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() + 1, add1.coord().j() + 1)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i(), add1.coord().j() + 1)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() - 5, DgIVec2D(maxI,add1.coord().j() + 1)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() - 5, DgIVec2D(maxI,add1.coord().j())));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() - 5, DgIVec2D(maxI,add1.coord().j() - 1)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i(), add1.coord().j() - 1)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() + 1, add1.coord().j() - 1)));
		}

	}
	else if(add1.coord().i() == maxI && add1.coord().j() != 0 && add1.coord().j() != maxJ)//最大行 但是列不是0或maxJ
	{
		if(add1.quadNum() >= 1 && add1.quadNum() <= 5)
		{
			neiadd.push_back(DgQ2DICoord(add1.quadNum() + 5, DgIVec2D(0, add1.coord().j())));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() + 5, DgIVec2D(0, add1.coord().j() + 1)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i(), add1.coord().j() + 1)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() - 1, add1.coord().j() + 1)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() - 1, add1.coord().j())));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() - 1, add1.coord().j() - 1)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i(), add1.coord().j() - 1)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() + 5, DgIVec2D(0, add1.coord().j() - 1)));
		}
		else if(add1.quadNum() >= 6 && add1.quadNum() < 10)
		{
			neiadd.push_back(DgQ2DICoord(add1.quadNum() + 1, DgIVec2D(maxJ - add1.coord().j(), 0)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() + 1, DgIVec2D(maxJ - add1.coord().j() - 1, 0)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i(), add1.coord().j() + 1)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() - 1, add1.coord().j() + 1)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() - 1, add1.coord().j())));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() - 1, add1.coord().j() - 1)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i(), add1.coord().j() - 1)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() + 1, DgIVec2D(maxJ - add1.coord().j() + 1, 0)));
		}
		else if(add1.quadNum() == 10)
		{
			neiadd.push_back(DgQ2DICoord(add1.quadNum() - 4, DgIVec2D(maxJ - add1.coord().j(), 0)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() - 4, DgIVec2D(maxJ - add1.coord().j() - 1, 0)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i(), add1.coord().j() + 1)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() - 1, add1.coord().j() + 1)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() - 1, add1.coord().j())));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() - 1, add1.coord().j() - 1)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i(), add1.coord().j() - 1)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() - 4, DgIVec2D(maxJ - add1.coord().j() + 1, 0)));
		}
	}
	else if(add1.coord().j() == 0 && add1.coord().i() != 0 && add1.coord().i() != maxJ) //第0列 但是行不是0 或最大行 
	{
		if(add1.quadNum() == 1)
		{
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() + 1, add1.coord().j())));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() + 1, add1.coord().j() + 1)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i(), add1.coord().j() + 1)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() - 1, add1.coord().j() + 1)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() - 1, add1.coord().j())));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() + 9, DgIVec2D(add1.coord().i() - 1, maxJ)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() + 9, DgIVec2D(add1.coord().i(), maxJ)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() + 9, DgIVec2D(add1.coord().i() + 1, maxJ)));
		}
		else if(add1.quadNum() > 1 && add1.quadNum() <= 5)
		{
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() + 1, add1.coord().j())));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() + 1, add1.coord().j() + 1)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i(), add1.coord().j() + 1)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() - 1, add1.coord().j() + 1)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() - 1, add1.coord().j())));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() + 4, DgIVec2D(add1.coord().i() - 1, maxJ)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() + 4, DgIVec2D(add1.coord().i(), maxJ)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() + 4, DgIVec2D(add1.coord().i() + 1, maxJ)));
		}
		else if(add1.quadNum() == 6)
		{
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() + 1, add1.coord().j())));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() + 1, add1.coord().j() + 1)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i(), add1.coord().j() + 1)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() - 1, add1.coord().j() + 1)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() - 1, add1.coord().j())));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() + 4, DgIVec2D(maxI, maxI - add1.coord().i() + 1)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() + 4, DgIVec2D(maxI, maxI - add1.coord().i())));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() + 4, DgIVec2D(maxI, maxI - add1.coord().i() - 1)));
		}
		else if(add1.quadNum() > 6 && add1.quadNum() <= 10)
		{
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() + 1, add1.coord().j())));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() + 1, add1.coord().j() + 1)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i(), add1.coord().j() + 1)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() - 1, add1.coord().j() + 1)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() - 1, add1.coord().j())));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() - 1, DgIVec2D(maxI, maxI - add1.coord().i() + 1)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() - 1, DgIVec2D(maxI, maxI - add1.coord().i())));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() - 1, DgIVec2D(maxI, maxI - add1.coord().i() - 1)));
		}
	}
	else if(add1.coord().j() == maxJ && add1.coord().i() != 0 && add1.coord().i() != maxJ)//最大列 但是行不是0 或最大行 
	{
		if(add1.quadNum() >= 1 && add1.quadNum() < 5)
		{
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() + 1, add1.coord().j())));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() + 1, DgIVec2D(0, maxI - add1.coord().i() - 1)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() + 1, DgIVec2D(0, maxI - add1.coord().i())));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() + 1, DgIVec2D(0, maxI - add1.coord().i() + 1)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() - 1, add1.coord().j())));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() - 1, add1.coord().j() - 1)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i(), add1.coord().j() - 1)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() + 1, add1.coord().j() - 1)));
		}
		else if(add1.quadNum() == 5)
		{
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() + 1, add1.coord().j())));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() - 4, DgIVec2D(0, maxI - add1.coord().i() - 1)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() - 4, DgIVec2D(0, maxI - add1.coord().i())));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() - 4, DgIVec2D(0, maxI - add1.coord().i() + 1)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() - 1, add1.coord().j())));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() - 1, add1.coord().j() - 1)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i(), add1.coord().j() - 1)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() + 1, add1.coord().j() - 1)));
		}
		else if(add1.quadNum() >= 6 && add1.quadNum() < 10)
		{
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() + 1, add1.coord().j())));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() - 4, DgIVec2D(add1.coord().i() + 1, 0)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() - 4, DgIVec2D(add1.coord().i(), 0)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() - 4, DgIVec2D(add1.coord().i() - 1, 0)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() - 1, add1.coord().j())));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() - 1, add1.coord().j() - 1)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i(), add1.coord().j() - 1)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() + 1, add1.coord().j() - 1)));
		}
		else if(add1.quadNum() == 10)
		{
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() + 1, add1.coord().j())));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() - 9, DgIVec2D(add1.coord().i() + 1, 0)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() - 9, DgIVec2D(add1.coord().i(), 0)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() - 9, DgIVec2D(add1.coord().i() - 1, 0)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() - 1, add1.coord().j())));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() - 1, add1.coord().j() - 1)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i(), add1.coord().j() - 1)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() + 1, add1.coord().j() - 1)));
		}
	}
	
	else if(add1.coord().j() == maxJ && add1.coord().i() == 0)//剩余的特殊情况 最大列 最小行 
	{
		//对于1~5属于极点共点，是锐角 所以多一个邻近 7~11 其余是钝角 所以少一个邻近
		if(add1.quadNum() == 1)//九个邻域
		{
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() + 1, add1.coord().j())));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() + 1, DgIVec2D(0,maxJ - 1)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() + 1, DgIVec2D(0,maxJ)));
            
			neiadd.push_back(DgQ2DICoord(add1.quadNum() + 2, DgIVec2D(0,maxJ)));

			neiadd.push_back(DgQ2DICoord(add1.quadNum() + 3, DgIVec2D(0,maxJ)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() + 4, DgIVec2D(0, maxJ)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() + 4, DgIVec2D(1, maxJ)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i(), add1.coord().j() - 1)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() + 1, add1.coord().j() - 1)));
		}
		else if(add1.quadNum() == 2)
		{
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() + 1, add1.coord().j())));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() + 1, DgIVec2D(0,maxJ - 1)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() + 1, DgIVec2D(0,maxJ)));

			neiadd.push_back(DgQ2DICoord(add1.quadNum() + 2, DgIVec2D(0,maxJ)));

			neiadd.push_back(DgQ2DICoord(add1.quadNum() + 3, DgIVec2D(0,maxJ)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() - 1, DgIVec2D(0, maxJ)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() - 1, DgIVec2D(1, maxJ)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i(), add1.coord().j() - 1)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() + 1, add1.coord().j() - 1)));
		}
		else if(add1.quadNum() == 3)
		{
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() + 1, add1.coord().j())));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() + 1, DgIVec2D(0,maxJ - 1)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() + 1, DgIVec2D(0,maxJ)));

			neiadd.push_back(DgQ2DICoord(add1.quadNum() + 2, DgIVec2D(0,maxJ)));

			neiadd.push_back(DgQ2DICoord(add1.quadNum() - 2, DgIVec2D(0,maxJ)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() - 1, DgIVec2D(0, maxJ)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() - 1, DgIVec2D(1, maxJ)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i(), add1.coord().j() - 1)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() + 1, add1.coord().j() - 1)));
		}
		else if(add1.quadNum() == 4)
		{
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() + 1, add1.coord().j())));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() + 1, DgIVec2D(0,maxJ - 1)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() + 1, DgIVec2D(0,maxJ)));

			neiadd.push_back(DgQ2DICoord(add1.quadNum() - 3, DgIVec2D(0,maxJ)));

			neiadd.push_back(DgQ2DICoord(add1.quadNum() - 2, DgIVec2D(0,maxJ)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() - 1, DgIVec2D(0, maxJ)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() - 1, DgIVec2D(1, maxJ)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i(), add1.coord().j() - 1)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() + 1, add1.coord().j() - 1)));
		}
		else if(add1.quadNum() == 5)
		{
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() + 1, add1.coord().j())));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() - 4, DgIVec2D(0,maxJ - 1)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() - 4, DgIVec2D(0,maxJ)));

			neiadd.push_back(DgQ2DICoord(add1.quadNum() - 3, DgIVec2D(0,maxJ)));

			neiadd.push_back(DgQ2DICoord(add1.quadNum() - 2, DgIVec2D(0,maxJ)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() - 1, DgIVec2D(0, maxJ)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() - 1, DgIVec2D(1, maxJ)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i(), add1.coord().j() - 1)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() + 1, add1.coord().j() - 1)));
		}
		else if(add1.quadNum() >= 6 && add1.quadNum() < 10)
		{
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() + 1, add1.coord().j())));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() - 4, DgIVec2D(1,0)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() - 4, DgIVec2D(0,0)));
// 重复第三个 保持八个邻近
			// neiadd.push_back(DgQ2DICoord(add1.quadNum() - 4, DgIVec2D(0,0)));

			neiadd.push_back(DgQ2DICoord(add1.quadNum() - 5, DgIVec2D(maxI, maxJ)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() - 5, DgIVec2D(maxI,maxJ - 1)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i(), add1.coord().j() - 1)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() + 1, add1.coord().j() - 1)));
		}
		else if(add1.quadNum() == 10)
		{
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() + 1, add1.coord().j())));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() - 9, DgIVec2D(1,0)));//源程序错误
			neiadd.push_back(DgQ2DICoord(add1.quadNum() - 9, DgIVec2D(0,0)));

			// 重复第三个 保持八个邻近	
			// neiadd.push_back(DgQ2DICoord(add1.quadNum() - 9, DgIVec2D(0,0)));


			neiadd.push_back(DgQ2DICoord(add1.quadNum() - 5, DgIVec2D(maxI, maxJ)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() - 5, DgIVec2D(maxI,maxJ - 1)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i(), add1.coord().j() - 1)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() + 1, add1.coord().j() - 1)));
		}
	}
	else if(add1.coord().j() == maxJ && add1.coord().i() == maxI)//最大列 最大行 七个邻域
	{
		//都是是钝角 所以少一个邻近
		if(add1.quadNum() >= 1 && add1.quadNum() < 5)
		{
			neiadd.push_back(DgQ2DICoord(add1.quadNum() + 5, DgIVec2D(0, maxJ)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() + 1, DgIVec2D(0,0)));
// 重复第二个 保持八个编码
			// neiadd.push_back(DgQ2DICoord(add1.quadNum() + 1, DgIVec2D(0,0)));

			neiadd.push_back(DgQ2DICoord(add1.quadNum() + 1, DgIVec2D(0,1)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() - 1, add1.coord().j())));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() - 1, add1.coord().j() - 1)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i(), add1.coord().j() - 1)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() + 5, DgIVec2D(0,maxJ - 1)));
		}
		else if(add1.quadNum() == 5)
		{
			neiadd.push_back(DgQ2DICoord(add1.quadNum() + 5, DgIVec2D(0, maxJ)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() - 4, DgIVec2D(0,0)));

// 重复第二个 保持八个编码
			// neiadd.push_back(DgQ2DICoord(add1.quadNum() - 4, DgIVec2D(0,0)));

			neiadd.push_back(DgQ2DICoord(add1.quadNum() - 4, DgIVec2D(0,1)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() - 1, add1.coord().j())));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() - 1, add1.coord().j() - 1)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i(), add1.coord().j() - 1)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() + 5, DgIVec2D(0,maxJ - 1)));
		}
		else if(add1.quadNum() >= 6 && add1.quadNum() < 10)
		{
			neiadd.push_back(DgQ2DICoord(add1.quadNum() + 1, DgIVec2D(0,0)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() - 4, DgIVec2D(maxI,0)));

// 重复第二个 保持八个编码
			// neiadd.push_back(DgQ2DICoord(add1.quadNum() - 4, DgIVec2D(maxI,0)));

			neiadd.push_back(DgQ2DICoord(add1.quadNum() - 4, DgIVec2D(maxI - 1,0)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() - 1, add1.coord().j())));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() - 1, add1.coord().j() - 1)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i(), add1.coord().j() - 1)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() + 1, DgIVec2D(1,0)));
		}
		else if(add1.quadNum() == 10)
		{
			neiadd.push_back(DgQ2DICoord(add1.quadNum() - 4, DgIVec2D(0,0)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() - 9, DgIVec2D(maxI,0)));

// 重复第二个 保持八个编码
			// neiadd.push_back(DgQ2DICoord(add1.quadNum() - 9, DgIVec2D(maxI,0)));

			neiadd.push_back(DgQ2DICoord(add1.quadNum() - 9, DgIVec2D(maxI - 1,0)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() - 1, add1.coord().j())));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() - 1, add1.coord().j() - 1)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i(), add1.coord().j() - 1)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() - 4, DgIVec2D(1,0)));
		}
	}
	else if(add1.coord().j() == 0 && add1.coord().i() == 0)//最小列 最小行 七个邻域
	{
		if(add1.quadNum() == 1)
		{
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() + 1, add1.coord().j())));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() + 1, add1.coord().j() + 1)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i(), add1.coord().j() + 1)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() + 4, DgIVec2D(maxI - 1, maxJ)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() + 4, DgIVec2D(maxI,maxJ)));
			
			// 重复第五个 保持八个编码	
			// neiadd.push_back(DgQ2DICoord(add1.quadNum() + 4, DgIVec2D(maxI,maxJ)));

			neiadd.push_back(DgQ2DICoord(add1.quadNum() + 9, DgIVec2D(0,maxJ)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() + 9, DgIVec2D(1,maxJ)));
		}
		else if(add1.quadNum() > 1 && add1.quadNum() <= 5)
		{
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() + 1, add1.coord().j())));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() + 1, add1.coord().j() + 1)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i(), add1.coord().j() + 1)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() - 1, DgIVec2D(maxI - 1, maxJ)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() - 1, DgIVec2D(maxI,maxJ)));
			// 重复第五个 保持八个编码	
			// neiadd.push_back(DgQ2DICoord(add1.quadNum() - 1, DgIVec2D(maxI,maxJ)));

			neiadd.push_back(DgQ2DICoord(add1.quadNum() + 4, DgIVec2D(0,maxJ)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() + 4, DgIVec2D(1,maxJ)));
		}
		else if(add1.quadNum() == 6)
		{
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() + 1, add1.coord().j())));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() + 1, add1.coord().j() + 1)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i(), add1.coord().j() + 1)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() - 5, DgIVec2D(maxI,1)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() - 5, DgIVec2D(maxI,0)));
			// 重复第五个 保持八个编码	
			// neiadd.push_back(DgQ2DICoord(add1.quadNum() - 5, DgIVec2D(maxI,0)));

			neiadd.push_back(DgQ2DICoord(add1.quadNum() + 4, DgIVec2D(maxI,maxJ)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() + 4, DgIVec2D(maxI,maxJ - 1)));
		}
		else if(add1.quadNum() > 6 && add1.quadNum() <= 10)
		{
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() + 1, add1.coord().j())));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() + 1, add1.coord().j() + 1)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i(), add1.coord().j() + 1)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() - 5, DgIVec2D(maxI,1)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() - 5, DgIVec2D(maxI,0)));
			// 重复第五个 保持八个编码	
			// neiadd.push_back(DgQ2DICoord(add1.quadNum() - 5, DgIVec2D(maxI,0)));

			neiadd.push_back(DgQ2DICoord(add1.quadNum() - 1, DgIVec2D(maxI,maxJ)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() - 1, DgIVec2D(maxI,maxJ - 1)));
		}
		
	}
	else if(add1.coord().j() == 0 && add1.coord().i() == maxI)		//最小列 最大行 与最大列和最小行是相反的
		{
			//最小列 最大行 与最大列和最小行是相反的
			if(add1.quadNum() == 1)
			{
				neiadd.push_back(DgQ2DICoord(add1.quadNum() + 5, DgIVec2D(0,0)));
				neiadd.push_back(DgQ2DICoord(add1.quadNum() + 5, DgIVec2D(0,1)));
				neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i(), add1.coord().j() + 1)));
				neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() - 1, add1.coord().j() + 1)));
				neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() - 1, add1.coord().j())));
				neiadd.push_back(DgQ2DICoord(add1.quadNum() + 9, DgIVec2D(maxI - 1,maxJ)));
				neiadd.push_back(DgQ2DICoord(add1.quadNum() + 9, DgIVec2D(maxI,maxJ)));
			
			// 重复第八个 保持八个编码	
				// neiadd.push_back(DgQ2DICoord(add1.quadNum() + 9, DgIVec2D(maxI,maxJ)));

			}
			else if(add1.quadNum() > 1 && add1.quadNum() <= 5)
			{
				neiadd.push_back(DgQ2DICoord(add1.quadNum() + 5, DgIVec2D(0,0)));
				neiadd.push_back(DgQ2DICoord(add1.quadNum() + 5, DgIVec2D(0,1)));
				neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i(), add1.coord().j() + 1)));
				neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() - 1, add1.coord().j() + 1)));
				neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() - 1, add1.coord().j())));
				neiadd.push_back(DgQ2DICoord(add1.quadNum() + 4, DgIVec2D(maxI - 1,maxJ)));
				neiadd.push_back(DgQ2DICoord(add1.quadNum() + 4, DgIVec2D(maxI,maxJ)));
			
			// 重复第八个 保持八个编码	
				// neiadd.push_back(DgQ2DICoord(add1.quadNum() + 4, DgIVec2D(maxI,maxJ)));

			}
			else if(add1.quadNum() == 6)
			{
				neiadd.push_back(DgQ2DICoord(add1.quadNum() + 1, DgIVec2D(maxI,0)));
				neiadd.push_back(DgQ2DICoord(add1.quadNum() + 1, DgIVec2D(maxI - 1,0)));
				neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i(), add1.coord().j() + 1)));
				neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() - 1, add1.coord().j() + 1)));
				neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() - 1, add1.coord().j())));
				neiadd.push_back(DgQ2DICoord(add1.quadNum() + 4, DgIVec2D(maxI,1)));
				neiadd.push_back(DgQ2DICoord(add1.quadNum() + 4, DgIVec2D(maxI,0)));

				neiadd.push_back(DgQ2DICoord(add1.quadNum() + 3, DgIVec2D(maxI,0)));// 去除一个保持八邻近

				neiadd.push_back(DgQ2DICoord(add1.quadNum() + 2, DgIVec2D(maxI,0)));
			}
			else if(add1.quadNum() == 7)
			{
				neiadd.push_back(DgQ2DICoord(add1.quadNum() + 1, DgIVec2D(maxI,0)));
				neiadd.push_back(DgQ2DICoord(add1.quadNum() + 1, DgIVec2D(maxI - 1,0)));
				neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i(), add1.coord().j() + 1)));
				neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() - 1, add1.coord().j() + 1)));
				neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() - 1, add1.coord().j())));
				neiadd.push_back(DgQ2DICoord(add1.quadNum() - 1, DgIVec2D(maxI,1)));
				neiadd.push_back(DgQ2DICoord(add1.quadNum() - 1, DgIVec2D(maxI,0)));

				neiadd.push_back(DgQ2DICoord(add1.quadNum() + 3, DgIVec2D(maxI,0))); //去除一个保持八邻近

				neiadd.push_back(DgQ2DICoord(add1.quadNum() + 2, DgIVec2D(maxI,0)));
			}
			else if(add1.quadNum() == 8)
			{
				neiadd.push_back(DgQ2DICoord(add1.quadNum() + 1, DgIVec2D(maxI,0)));
				neiadd.push_back(DgQ2DICoord(add1.quadNum() + 1, DgIVec2D(maxI - 1,0)));
				neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i(), add1.coord().j() + 1)));
				neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() - 1, add1.coord().j() + 1)));
				neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() - 1, add1.coord().j())));
				neiadd.push_back(DgQ2DICoord(add1.quadNum() - 1, DgIVec2D(maxI,1)));
				neiadd.push_back(DgQ2DICoord(add1.quadNum() - 1, DgIVec2D(maxI,0)));

				neiadd.push_back(DgQ2DICoord(add1.quadNum() - 2, DgIVec2D(maxI,0)));//去除一个保持八邻近

				neiadd.push_back(DgQ2DICoord(add1.quadNum() + 2, DgIVec2D(maxI,0)));
			}
			else if(add1.quadNum() == 9)
			{
				neiadd.push_back(DgQ2DICoord(add1.quadNum() + 1, DgIVec2D(maxI,0)));
				neiadd.push_back(DgQ2DICoord(add1.quadNum() + 1, DgIVec2D(maxI - 1,0)));
				neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i(), add1.coord().j() + 1)));
				neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() - 1, add1.coord().j() + 1)));
				neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() - 1, add1.coord().j())));
				neiadd.push_back(DgQ2DICoord(add1.quadNum() - 1, DgIVec2D(maxI,1)));
				neiadd.push_back(DgQ2DICoord(add1.quadNum() - 1, DgIVec2D(maxI,0)));

				neiadd.push_back(DgQ2DICoord(add1.quadNum() - 2, DgIVec2D(maxI,0)));//去除一个保持八邻近

				neiadd.push_back(DgQ2DICoord(add1.quadNum() - 3, DgIVec2D(maxI,0)));
			}
			else if(add1.quadNum() == 10)
			{
				neiadd.push_back(DgQ2DICoord(add1.quadNum() - 4, DgIVec2D(maxI,0)));
				neiadd.push_back(DgQ2DICoord(add1.quadNum() - 4, DgIVec2D(maxI - 1,0)));
				neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i(), add1.coord().j() + 1)));
				neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() - 1, add1.coord().j() + 1)));
				neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() - 1, add1.coord().j())));
				neiadd.push_back(DgQ2DICoord(add1.quadNum() - 1, DgIVec2D(maxI,1)));
				neiadd.push_back(DgQ2DICoord(add1.quadNum() - 1, DgIVec2D(maxI,0)));

				neiadd.push_back(DgQ2DICoord(add1.quadNum() - 2, DgIVec2D(maxI,0))); //去除一个保持八邻近

				neiadd.push_back(DgQ2DICoord(add1.quadNum() - 3, DgIVec2D(maxI,0)));
			}
		}
    return neiadd;
}

calDgg::calDgg(  const DgIDGGBase &idgg) : IDGG_(idgg) 
{

}

void calDgg::setAddVertices(const DgQ2DICoord &add, DgPolygon &vec, int densify) const 
{
    /*得到的直接是 Dggecoord*/
    idgg().setAddVertices(add,vec,densify);

}

// long double DgGeoCoord::geoPolyArea
