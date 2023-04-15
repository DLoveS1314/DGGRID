//
// Created by dls on 23-3-15.
//

#include "calDgg.h"
#include "dglib/DgEllipsoidRF.h"
#include "dglib/DgBoundedIDGG.h"
#include <dglib/DgOutputStream.h>
#include <dglib/DgInputStream.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <dirent.h>
#include <chrono>
#include <omp.h>


bool calDgg::useearthRadius = true;
string getbasename(string filename);
long double calDgg::calarea(const DgQ2DICoord &add, int ptsPerEdgeDensify) const
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
	double r2 =1.0;
	if(calDgg::useearthRadius)
	{
		r2=DgGeoSphRF::earthRadiusKM()*DgGeoSphRF::earthRadiusKM() ;
	}
    double area =  this->geoPolyArea(verts, *center) * r2 ;
    // cout << "* area  " << area*DgGeoSphRF::earthRadiusKM()*DgGeoSphRF::earthRadiusKM() << endl;
    // delete loc;
    return  area;
}

long double calDgg::calper( const DgQ2DICoord &add, int ptsPerEdgeDensify) const
{
	DgPolygon verts;
    // shared_ptr<DgLocation> loc(idgg().makeLocation(add) ) ;/
    this->setAddVertices(add, verts, ptsPerEdgeDensify);
	double r = 1.0;
    if (calDgg::useearthRadius)//如果使用地球半径 那就乘上
	{
        r = DgGeoSphRF::earthRadiusKM() ;
	} 
	double per =  this->geoPolyper(verts ) * r ;
   return per;
}

long double calDgg::calzsc(const DgQ2DICoord &add, int ptsPerEdgeDensify) const
{
	// 分母
	long double per = this->calper(add,ptsPerEdgeDensify);


	long double area = this->calarea(add,ptsPerEdgeDensify);
	double r2 =1.0;
	if(calDgg::useearthRadius)
	{
		r2=DgGeoSphRF::earthRadiusKM()*DgGeoSphRF::earthRadiusKM() ;
	}
	double const pi = 3.141592653589793238;
	long double sqrtdata = 4*pi*area -area*area/r2 ; //参考罗师姐论文
	// 分子
	long double Numerator =  sqrt(sqrtdata);
	
	long double zsc =  Numerator/per;
   	return zsc;
}

long double calDgg::calzsc(const double area, const double per) const
{
	double r2 =1.0;
	if(calDgg::useearthRadius)
	{
		r2=DgGeoSphRF::earthRadiusKM()*DgGeoSphRF::earthRadiusKM() ;
	}
	long double sqrtdata = 4*M_PI*area -area*area/r2 ; //参考罗师姐论文
	// 分子
	long double Numerator =  sqrt(sqrtdata);
	
	long double zsc =  Numerator/per;
   	return zsc;
}

long double calDgg::geoPolyArea (const DgPolygon& poly, const DgGeoCoord& center)  const

{

    // 保证参考系已经是DgGeoSphRF空间才能进行计算
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


long double calDgg::geoPolyper (const DgPolygon& poly )  const

{
 
    long double totper = 0.0L;

    const DgGeoSphRF* geoRF = dynamic_cast<const DgGeoSphRF*>(&poly.rf());
    if (geoRF == 0)
        report("DgGeoCoord::geoPolyArea() non-geo polygon", DgBase::Fatal);

    // do each sub-triangle
    for (int i = 0; i < poly.size(); i++)
    {
        const DgGeoCoord& v1 = *geoRF->getAddress(poly[i]);
        const DgGeoCoord& v2 = *geoRF->getAddress(poly[(i + 1) % poly.size()]); //通过取余实现 roll

        totper += DgGeoCoord::gcDist( v1, v2);
    }

    return totper;

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

vector<long double> calDgg::calmaxmintotaldis(const DgQ2DICoord &add) const
{
	shared_ptr<DgLocation> loc(idgg().makeLocation(add) ) ;

    // const DgGeoCoord* center  = idgg().geoRF().getAddress(*idgg().geoRF().convert(loc.get()));
	vector<long double> returnval ;
	DgLocVector locv(this->idgg());
	this->calnei(add,locv);//不需要重新排序 
	long double minVal = LDBL_MAX ;
	long double maxVal = 0; //加个大写避免与类中函数名minval重复 防止出现意外情况
	long double totalval = 0;

	// do each point
    for (int i = 0; i < locv.size(); i++)
    {
		long double dis = this->caldistance(add,*idgg().getAddress(locv[i]));
       	totalval+=dis;
		minVal =this->minval(minVal,dis);
		maxVal =  this->maxval(maxVal,dis);

    }
	long double max_min_ratio = maxVal/minVal;
	long double mean = totalval/locv.size();
	returnval.push_back(max_min_ratio);
	returnval.push_back(mean);
	// returnval.push_back(maxval);

    return returnval;
}

// calmaxmintotalcsd 的辅助函数 不再放在类里 不再计算公共边
// vector<DgGeoCoord> Getmindis(vector<DgAddressBase*>add1,vector<DgAddressBase*>add2)  
// {
// 	//先转化为
// 		// DgAddress<DgVertex2DDCoord>* fullAdd =
//     //                dynamic_cast< DgAddress<DgVertex2DDCoord>* >(v[i]);
// 	vector<DgGeoCoord> coords;
// 	return coords

// }
vector<long double>  calDgg::calmaxmintotalcsd(const DgQ2DICoord& add) const 
{

	//求邻域
	vector<long double> mid_dis_vec;//按照默认的neis顺序记录每个邻近的两中点间距离

	shared_ptr<DgLocation> loc( idgg().makeLocation(add) ) ;

    // const DgGeoCoord* center  = idgg().geoRF().getAddress(*idgg().geoRF().convert(loc.get()));

	// 先求得格网本身的边界中点
	DgPolygon verts_c; //rf是 georf
	this->setAddVertices(add,verts_c,0);
	// vector<DgAddressBase*>& v_c = verts_c.addressVec();
	vector<DgGeoCoord> midpoints_c;

	//求面邻近：暂且无法求三角形的 后续把师姐给的三角形面代码加上
	vector<DgQ2DICoord> neis = this->edgecell(add);
	//邻域循环 得到每一个面csd 两中心点距离
	for (int nei = 0; nei < neis.size(); nei++)
    {

		//计算两个格网间中心点位置

		DgGeoCoord nei_midpoint = this->calmidpoint(add,neis[nei]);

		//循环计算 与格网本身的边界中点的距离 求的的最小距离就是真正的两中点间距离
		long double min_mid_dis = LDBL_MAX;
		for (int mp=0;mp <midpoints_c.size();mp++)
		{
			long double mid_dis = this->caldistance(nei_midpoint,midpoints_c[mp]);
			min_mid_dis =this->minval(min_mid_dis,mid_dis);
		}
		mid_dis_vec.push_back(min_mid_dis);

    }
	vector<long double> returnval ;
 
	long double minVal = LDBL_MAX ;
	long double maxVal = 0; //加个大写避免与类中函数名minval重复 防止出现意外情况
	long double totalval = 0;
	for(int m_v =0;m_v <mid_dis_vec.size();m_v++ ) 
	{
		long double dis =mid_dis_vec[m_v];
       	totalval+=dis;
		minVal =this->minval(minVal,dis);
		maxVal =  this->maxval(maxVal,dis);
	}
  
	long double max_min_ratio = maxVal/minVal;
	long double mean = totalval/mid_dis_vec.size();
	returnval.push_back(max_min_ratio);
	returnval.push_back(mean);
	// returnval.push_back(maxval);

    return returnval;

}

 
long double calDgg::calAzimuth(const GeoCoord &pt1, const GeoCoord &pt2) const
{
	long double az;
	az=atan2l(cosl(pt2.lat)*sinl(pt2.lon-pt1.lon), cosl(pt1.lat)*sinl(pt2.lat)-
		sinl(pt1.lat)*cosl(pt2.lat)*cosl(pt2.lon-pt1.lon));
	az=az+M_PI;//转为0~360度
	return az ;
}

long double calDgg::calAzimuth(const DgGeoCoord &pt1, const DgGeoCoord &pt2) const
{   
	// 计算p2为终点 的方位角
   	long double az;
	az=atan2l(cosl(pt2.lat())*sinl(pt2.lon()-pt1.lon()), cosl(pt1.lat())*sinl(pt2.lat())-
		sinl(pt1.lat())*cosl(pt2.lat())*cosl(pt2.lon()-pt1.lon()));
	az=az+M_PI;//转为0~360度
	return az ;
}
// // 自定义结构体，用于保存double类型的值和对应的索引
// struct IndexedValue {
//     long double value;
//     int index;
// };

// // 比较函数，用于std::sort排序时按value从小到大排序
// bool cmp(const IndexedValue &a, const IndexedValue &b) {
//     return a.value < b.value;
// }

vector<long double > calDgg::sortByAzimuth(const DgQ2DICoord& add,const DgLocVector locs) const
{
	// 先计算方位角
	vector<long double > azi;
	// vector< IndexedValue > azi( locs.size());
	DgGeoCoord add_c = this->getGeocoord(add);
	// for (int i = 0; i < locs.size(); i++)
	// {
	// 	DgGeoCoord add_e = this->getGeocoord( *idgg().getAddress(locs[i]));
	// 	azi[i].value = this->calAzimuth(add_c,add_e);
	// 	azi[i].index=i;
	// }
		for (int i = 0; i < locs.size(); i++)
	{
		DgGeoCoord add_e = this->getGeocoord( *idgg().getAddress(locs[i]));
		azi.push_back(this->calAzimuth(add_c,add_e)) ;
	} 
	sort(azi.begin(), azi.end());//默认升序排列 必须从小到大排序 通过方位角计算的夹角才是正确的
    return azi;
}

vector<DgQ2DICoord> calDgg::readadd(string path) const
{
	vector<DgQ2DICoord> add;
	//先读取只包含编码的文档
	// 打开文件
	std::ifstream file(path);
	// 检查文件是否成功打开
	if (!file.is_open()) {
		std::cout << "Failed to open file." << std::endl;
		return vector<DgQ2DICoord>();
	}
	// 逐行读取文件内容
	std::vector<int> nums;
	std::string line;
	while (std::getline(file, line)) {
		// 将字符串转换为整数并存储在容器中
		nums.push_back(std::stoull(line));
	}

	file.close();
	//把编码转化为DgQ2DICoord
	for(auto num : nums)
	{
		add.push_back(idgg().bndRF().addFromSeqNum(num)) ;
	}
    return add;
}

vector<long double> calDgg::calmaxmintotalangle(const DgQ2DICoord &add) const
{
	shared_ptr<DgLocation> loc(idgg().makeLocation(add) ) ;

    // const DgGeoCoord* center  = idgg().geoRF().getAddress(*idgg().geoRF().convert(loc.get()));
	vector<long double> returnval ;
	DgLocVector locv(this->idgg());
	this->calnei(add,locv);
	//三角形的邻近顺序不固定怎么办 只能先算方位角 然后按照方位角排序
	vector<long double> azi = this->sortByAzimuth(add,locv);
	long double minVal = LDBL_MAX ;
	long double maxVal = 0; //加个大写避免与类中函数名minval重复 防止出现意外情况
	long double totalval = 0;
	long double angle =0;
	// do each point  通过方位角计算夹角最后一个夹角单独计算
    for (int i = 0; i < azi.size()-1; i++)
    {
		 
		angle = azi[i+1]-azi[i];
		if (angle>M_PI) //防止钝角出现 虽然应该不可能
		{
			angle =M_2_PI- angle;
		}
       	totalval+=angle;
		minVal =this->minval(minVal,angle);
		maxVal =  this->maxval(maxVal,angle);
    }
	// 最后一个角一般是钝角
	angle = azi[azi.size()-1]-azi[0];//一般是钝角
	if (angle>M_PI) //防止钝角出现 虽然应该不可能
	{
		angle =M_2_PI- angle;
	}
	totalval+=angle;
	minVal =this->minval(minVal,angle);
	maxVal =  this->maxval(maxVal,angle);

	long double max_min_ratio = maxVal/minVal;
	long double mean = totalval/locv.size();

	returnval.push_back(max_min_ratio);
	returnval.push_back(mean);
	// returnval.push_back(maxval);

    return returnval;
}
// useearthRadius 默认是true
long double calDgg::caldistance(const DgQ2DICoord &add,const  DgQ2DICoord &add1  ) const
{
    DgGeoCoord addg = this->getGeocoord(add);
    DgGeoCoord addg1 = this->getGeocoord(add1);
    long double distance =   DgGeoCoord::gcDist(addg, addg1); 
    if (calDgg::useearthRadius)//如果使用地球半径 那就乘上
      {
        distance = DgGeoSphRF::earthRadiusKM() * distance;
       } 
    return distance;
}

long double calDgg::caldistance(const DgGeoCoord &add, const DgGeoCoord &add1 ) const
{
	long double distance =   DgGeoCoord::gcDist(add, add1); 
    if (calDgg::useearthRadius)//如果使用地球半径 那就乘上
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

DgGeoCoord calDgg::calmidpoint(const DgGeoCoord &add, const DgGeoCoord &add1) const
{
	DgGeoCoord midpoint =   DgGeoSphRF::midPoint(add, add1); 
    return midpoint;
}

DgGeoCoord calDgg::calintersect(const GeoCoord &sv11, const GeoCoord &sv12, const GeoCoord &sv21, const GeoCoord &sv22, int sign) const
{   
  /*
  参考的 GCintersect
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

DgGeoCoord calDgg::calintersect(const DgGeoCoord &g11, const DgGeoCoord &g12, const DgGeoCoord &g21, const DgGeoCoord &g22, int sign) const
{
//    GeoCoord sv11, sv12, sv21,sv22;
//    sv11.lon = g11.lon(); sv11.lat = g11.lat();
//    sv12.lon = g12.lon(); sv12.lat = g12.lat();
//    sv21.lon = g21.lon(); sv21.lat = g21.lat();

//    sv22.lon = g22.lon(); sv22.lat = g22.lat();
    return DgGeoCoord();
}

/// @brief Calculate the adjacent grid of idgg, only realize the rhombus
/// @param add Piece number row number column number
/// @param vec saved structure 保存生成的邻域编码 编码是QIj 的location 
void calDgg::calnei(const DgQ2DICoord &add, DgLocVector &vec) const
{
	//三角形的邻近顺序不是按照逆时针顺序来的一定谨记！
	// 直接使用dggride的源程序 这个里面 菱形经过了特殊处理 去 除了九邻近，和七邻近 都补充为了八邻近
	// 三角形 输出的是九邻近 六边形输出的是六邻近或者五邻近， 按照卷积核的大小输出来处理
    //DgLocVector在 DgDiscRF setNeighbors 已经转化为了DGIDGG空间 也就是Q i j空间
   if (this->idgg().gridTopo()== Diamond) 
   //菱形的在源程序里菱形经过了特殊处理 去 除了九邻近，和七邻近 只包含了八邻近 所以不能使用 这里重新写一下
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
		idgg().setNeighbors(add,vec);	
   }
//    else
//    {
//     dgcerr<<"setAddNeighbors "<<this->idgg().gridTopo()<< " not yet implemented"<<endl;
//    }
}


long double calDgg::calangle(const DgQ2DICoord &center, const DgQ2DICoord &add1, const DgQ2DICoord &add2) const
{
    // 改编 自DgGeoCoord::geoTriArea 返回的是弧度 
	//  也可以计算方位角 然后返回锐角
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

    // long double l1[2],l2[2],l3[2],edges[2];
    // cout<<"夹角"<<bigA<<endl;

	// 经过测试两个结果相同
    // l1[0]=centerg.lat(); l1[1]=centerg.lon();
    // l2[0]=addg1.lat(); l2[1]=addg1.lon();
    // l3[0]=addg2.lat(); l3[1]=addg2.lon();
    // edges[0]=acosl(cosl(M_PI_2-l2[0])*cosl(M_PI_2-l3[0])+
    //         sinl(M_PI_2-l2[0])*sinl(M_PI_2-l3[0])*cosl(l2[1]-l3[1]));
    // edges[1]=acosl(cosl(M_PI_2-l1[0])*cosl(M_PI_2-l3[0])+
    //         sinl(M_PI_2-l1[0])*sinl(M_PI_2-l3[0])*cosl(l1[1]-l3[1]));
    // edges[2]=acosl(cosl(M_PI_2-l2[0])*cosl(M_PI_2-l1[0])+
    //         sinl(M_PI_2-l2[0])*sinl(M_PI_2-l1[0])*cosl(l2[1]-l1[1]));
    // long double angles =acosl((cosl(edges[0])-cosl(edges[1])*
    //            cosl(edges[2]))/(sinl(edges[1])*sinl(edges[2])));
    // cout<<"夹角"<<angles<<endl;
	return  bigA ;
}
vector<DgQ2DICoord> calDgg::edgecell(const DgQ2DICoord &add) const
	{
		// 改写自 DgIDGGBase::setAddNeighbors 只能处理菱形和六边形 三角形不可以 所以增加一个if 
		if (this->idgg().gridTopo()== Triangle)
		{
				vector <DgQ2DICoord> coords = this->idgg().triedgenei(add);
				return  coords;
		}
		else
		{	
			vector<DgQ2DICoord> vec;
			// //重写了 DgDiscRF setAddNeighbors
			//DgLocVector在 DgDiscRF setNeighbors 已经转化为了DGIDGG空间 也就是Q i j
			//空的vector d但是 参考系变成了grid2D() 先求在一个面片里的行列号 边界处有错误的话 在后面改写
			DgLocVector ngh2d(idgg().grid2D());
			idgg().grid2D().setAddNeighbors(add.coord(), ngh2d);
			int q = add.quadNum();
			DgLocVector ngh2dNoDup(idgg());//开始处理跨面
			DgIVec2D c;
			for (int i = 0; i < ngh2d.size(); i++)
			{
				DgQ2DICoord c2di(q, *idgg().grid2D().getAddress(ngh2d[i])); //这个时候是有错误的 因为他包含了边界点
				const DgBoundedIDGG& bndRF =idgg().bndRF();
				c2di = bndRF.q2dixToQ2di(c2di);//纠正错误的边界编码
				// check for duplicates
				bool keeper = true;
				if (!idgg().isCongruent() && add.coord() == DgIVec2D(0, 0))
				{
					for (int i = 0; i < ngh2dNoDup.size(); i++)
					{
						const DgQ2DICoord& veci = *idgg().getAddress(ngh2dNoDup[i]);
						if (c2di == veci) 
						{
							keeper = false;
							break;
						}
					}
				}
				if (keeper)
				{

					vec.push_back(c2di);
				}
			}
			return vec;

		}
	
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



// void outchunk(string outFileName ,vector<results> resus,vector<string> names )
// {
// 		// vector<DgQ2DICoord> adds =getQDIfromfile(inFileName,inAddType,outAddType);
// 		// calonefile(adds);
//         // shared_ptr<ofstream> pOutFile(new DgOutputStream(outFileName, "", DgBase::Fatal));
//         // ofstream& outFile = *pOutFile;
// 		std::ofstream outFile(outFileName,std::ios_base::app); // 创建 ofstream 并打开文件
// 		if (!outFile.is_open()) 
// 		{ // 检查文件是否成功打开
//         	std::cerr << " result2file Failed to open file " <<outFileName<< std::endl;
// 			exit(0);
// 		}
// 		//首先输入标题
// 		// const char* title = " name , area , per , zsc , disminmax ,disavg , angleminmax , angleavg , csdminmax , csdavg ";
// 		// outFile<<title<<endl;
// 		for(int i =0;i<resus.size() ;i++)
// 		{
// 			outFile<<names[i];
// 			outFile<<outputDelimiter<<resus[i].area;
// 			outFile<<outputDelimiter<<resus[i].per;
// 			outFile<<outputDelimiter<<resus[i].zsc;
// 			outFile<<outputDelimiter<<resus[i].disminmax;
// 			outFile<<outputDelimiter<<resus[i].disavg;
// 			outFile<<outputDelimiter<<resus[i].angleminmax;
// 			outFile<<outputDelimiter<<resus[i].angleavg;
// 			outFile<<outputDelimiter<<resus[i].csdminmax;
// 			outFile<<outputDelimiter<<resus[i].csdavg;
// 			outFile<<endl;
// 		}
// 		// for(int index =0;index<resus.size() ;index++)
// 		// {
// 		// 	outFile<<resu.name;
// 		// 	outFile<<outputDelimiter<<resu.area;
// 		// 	outFile<<outputDelimiter<<resu.per;
// 		// 	outFile<<outputDelimiter<<resu.zsc;
// 		// 	outFile<<outputDelimiter<<resu.disminmax;
// 		// 	outFile<<outputDelimiter<<resu.disavg;
// 		// 	outFile<<outputDelimiter<<resu.angleminmax;
// 		// 	outFile<<outputDelimiter<<resu.angleavg;
// 		// 	outFile<<outputDelimiter<<resu.csdminmax;
// 		// 	outFile<<outputDelimiter<<resu.csdavg;
// 		// 	outFile<<endl;
// 		// }
// 		//循环输入文件
//    		outFile.close();
// }


/// @brief 
/// @param outFileName  输出的文件名 
/// @param outputDelimiter 输出所用的间隔符
/// @param inFilepath 输入文件所在的路径 
/// @param inAddType 输入的编码类型 见DGGRID 能够作为输入的编码
/// @param outAddType 输出的编码类型 默认且必须是"Q2DI"
/// @param suffix 输入文件的后缀名
void calDgg::result2file(string outFileName, char outputDelimiter,string inFilePath, string inAddType,string outAddType,string suffix,int ptsPerEdgeDensify  ) const
{

		vector <string> paths =this->getfilName(inFilePath,suffix);
		vector<results> resus;
		vector<string> names;
		vector<string> filenames;
		// vector<DgQ2DICoord> adds;
		// for(int i =0;i< paths.size(); i++)
		// {
		// 	fs::path fspath(paths[i]);
		// 	std::string filename = fspath.filename().string();
		// 	string name = getbasename(filename);
		// 	names.push_back(name);
		// 	filenames.push_back(paths[i]);
		// }
		// int num_threads = 4;

		// #pragma omp parallel num_threads(num_threads)
		// #pragma omp parallel 

		{
			int thread_id = omp_get_thread_num();
			int num_threads = omp_get_num_threads();
			// #pragma omp for
			for(int i =0;i< filenames.size(); i++)
			{
				// #pragma omp critical
				// {
				std::cout << "result2file Thread " << thread_id << " of " << num_threads   << " "<<filenames[i]<<std::endl;
					// adds =this->getQDIfromfile(filenames[i]);
					// adds.empty();
				string path =filenames[i];
				// this->getQDIfromfile(path);
				// }
				results resu =this->calonefile(this->getQDIfromfile(path), ptsPerEdgeDensify);
				// #pragma omp critical
				{
					resus.push_back(resu);
				}
			}
		}
		std::ofstream outFile(outFileName); // 创建 ofstream 并打开文件
		if (!outFile.is_open()) 
		{ // 检查文件是否成功打开
        	std::cerr << " result2file Failed to open file " <<outFileName<< std::endl;
			exit(0);
		}

		//首先输入标题
 
		const char* title = " name , area , per , zsc , disminmax ,disavg , angleminmax , angleavg , csdminmax , csdavg ";
		outFile<<title<<endl;
		for(int i =0;i<resus.size() ;i++)
		{
			outFile<<names[i];
			outFile<<outputDelimiter<<resus[i].area;
			outFile<<outputDelimiter<<resus[i].per;
			outFile<<outputDelimiter<<resus[i].zsc;

			outFile<<outputDelimiter<<resus[i].disminmax;
			outFile<<outputDelimiter<<resus[i].disavg;

			outFile<<outputDelimiter<<resus[i].angleminmax;
			outFile<<outputDelimiter<<resus[i].angleavg;

			outFile<<outputDelimiter<<resus[i].csdminmax;
			outFile<<outputDelimiter<<resus[i].csdavg;
			outFile<<endl;
		}
   		outFile.close();

}

void calDgg::result2file1(string outFileName, char outputDelimiter, string inFilePath, string inAddType, string outAddType, string suffix, int ptsPerEdgeDensify) const
{
		vector <string> paths =this->getfilName(inFilePath,suffix);
		vector<cloumnsdata> cdatas(paths.size());
		// vector<DgQ2DICoord> adds;
		for(auto ii =0;ii< paths.size(); ii++)
		{
			cdatas[ii].setname(paths[ii]);
			cdatas[ii].adds =this->getQDIfromfile(paths[ii]);
		}
  		// omp_set_num_threads(10);
		// #pragma omp parallel  
		{
			int thread_id = omp_get_thread_num();
			int num_threads = omp_get_num_threads();
			// #pragma omp for  
			for(auto i =0; i<cdatas.size(); i++)
			// for(auto i =0; i<10; i++)
			{
				// #pragma omp critical
				{
					std::cout << i<<" result2file Thread " << thread_id << " of " << num_threads   << " "<<cdatas[i].name<<std::endl;
				}
				cdatas[i].resu =this->calonefile(cdatas[i].adds, ptsPerEdgeDensify);
			}

		} 
		std::ofstream outFile (outFileName); // 创建 ofstream 并打开文件
		if (!outFile.is_open()) 
		{ // 检查文件是否成功打开
        	std::cerr << " result2file Failed to open file " <<outFileName<< std::endl;
			exit(0);
		}
		const char* title = " name , area , per , zsc , disminmax ,disavg , angleminmax , angleavg , csdminmax , csdavg ";
		outFile<<title<<endl;
		for(int i =0;i<cdatas.size() ;i++)
		{
			outFile<<cdatas[i].name;
			outFile<<outputDelimiter<<cdatas[i].resu.area;
			outFile<<outputDelimiter<<cdatas[i].resu.per;
			outFile<<outputDelimiter<<cdatas[i].resu.zsc;

			outFile<<outputDelimiter<<cdatas[i].resu.disminmax;
			outFile<<outputDelimiter<<cdatas[i].resu.disavg;

			outFile<<outputDelimiter<<cdatas[i].resu.angleminmax;
			outFile<<outputDelimiter<<cdatas[i].resu.angleavg;

			outFile<<outputDelimiter<<cdatas[i].resu.csdminmax;
			outFile<<outputDelimiter<<cdatas[i].resu.csdavg;
			outFile<<endl;
		}
   		outFile.close();
		// cdatas.empty();
}

vector<DgQ2DICoord> calDgg::getQDIfromfile (string inFileName,    string inAddType ,string outAddType  )  const
{
	const DgIDGGBase&  dgg= this->idgg();
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
		inSeqNum = true;
		pInRF = &dgg;
	}
	const DgRFBase& inRF = *pInRF;

	const DgRFBase& outRF =dgg;

	// set the precision

	// const_cast<DgRFBase&>(outRF).setPrecision(precision);

	// now process the addresses in the input file

	const int maxLine = 1000;
	char buff[maxLine];

	// dgcout << "transforming values..." << endl;

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
		loc = NULL;

	}
	inFile.close();

	return adds;
}
vector<string> calDgg::getfilName(string inFilePath, string suffix) const
{
	vector<string> absfilenames;
    DIR *dir = opendir(inFilePath.c_str());
    if (dir)
    {
        struct dirent *entry;
        while ((entry = readdir(dir)) != nullptr)
        {
            if (entry->d_type == DT_REG && std::string(entry->d_name).substr(std::string(entry->d_name).find_last_of(".") + 1) == suffix)
            {
				absfilenames.push_back(inFilePath + "/" + entry->d_name);
				// string filename = string(path);
                // std::cout << "Found CSV file: " << inFilePath << "/" << entry->d_name << std::endl;
            }
        }
        closedir(dir);
		return absfilenames;
    }
    else
    {
        std::cerr << "Could not open directory: " << inFilePath << std::endl;
        return absfilenames;
    }
}
long double calmean(vector<long double> vecs)
{
	long double sum = 0;
	for (int i = 0; i < vecs.size(); i++) {
        sum += vecs[i];
    }
 	double mean = sum / vecs.size();
	return mean;
}
results calDgg::calonefile(vector<DgQ2DICoord> adds,int ptsPerEdgeDensify) const
{
	results resu;

	vector<  long double > areas;
	vector<  long double > pers;
	vector<  long double > zscs;

	vector<  long double > disminmaxs;
	vector<  long double > disavgs;

	vector<  long double > angleminmaxs;
	vector<  long double > angleavgs;

	vector<  long double > csdminmaxs;
	vector<  long double > csdavgs;

	 for(auto add : adds)
	 {
		long double area =this->calarea(add,ptsPerEdgeDensify);
		long double per = this->calper(add,ptsPerEdgeDensify);
		areas.push_back( area);
		pers.push_back(per);
		zscs.push_back(this->calzsc(area,per));
		// cout<<"1"<<endl;

		vector<long double> disvec = this->calmaxmintotaldis(add);
		disminmaxs.push_back(disvec[0]);
		disavgs.push_back(disvec[1]);
		// cout<<"2"<<endl;
		
		vector<long double> anglevec = this->calmaxmintotalangle(add);
		angleminmaxs.push_back(anglevec[0]);
		angleavgs.push_back(anglevec[1]);
		// cout<<"3"<<endl;

		vector<long double> csdvec = this->calmaxmintotalangle(add);
		csdminmaxs.push_back(csdvec[0]);
		csdavgs.push_back(csdvec[1]);
		// cout<<"4"<<endl;

	 }

	 resu.area = calmean(areas);
	 resu.per = calmean(pers);
	 resu.zsc = calmean(zscs);

	 resu.angleminmax  = calmean(angleminmaxs);
	 resu.angleavg  = calmean(angleavgs);

	 resu.disminmax = calmean(disminmaxs);
	 resu.disavg  = calmean(disavgs);

	 resu.csdminmax  = calmean(csdminmaxs);
	 resu.csdavg = calmean(csdavgs);


    return resu;
}

void calDgg::setAddVertices(const DgQ2DICoord &add, DgPolygon &vec, int densify) const 
{
    /*得到的直接是 Dggecoord*/
    idgg().setAddVertices(add,vec,densify);

}

// long double DgGeoCoord::geoPolyArea

unsigned long long int calDgg::getseqnum(const DgLocation loc) const
{
    return static_cast<const DgIDGGBase&>(this->idgg()).bndRF().seqNum(loc);
}

DgQ2DICoord calDgg::locFromSeqNum(unsigned long long int sNum) const
{
	shared_ptr<DgLocation> loc (static_cast<const DgIDGGBase&>(this->idgg()).bndRF().locFromSeqNum(sNum));
	const DgQ2DICoord add = static_cast<const DgAddress<DgQ2DICoord>*>(loc->address())->address();
    return add;
}
 