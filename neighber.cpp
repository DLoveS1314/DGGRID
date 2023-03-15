#include <iostream>
#include <fstream>
//#include "DgIDGGS4H.h"
#include "DgIDGGS4D.h"
#include "DgIDGG.h"
#include "DgBoundedIDGG.h"
#include "DgEllipsoidRF.h"
#include "ogrsf_frmts.h"
#include <gdal.h>
#include <gdal_priv.h>
#include<ogr_spatialref.h>
#include <cpl_string.h>
#include <cpl_conv.h>

using namespace std;

struct polygon
{
	vector<GeoCoord> Vertices;
};

class MyLatLng //经纬度为数值
{ 
public:
	double m_Longitude, m_Latitude;  
	MyLatLng()
	{
		m_Longitude=0;
		m_Latitude=0;
	}
	MyLatLng(double longitude,double latitude)
	{  
		m_Longitude=longitude;
		m_Latitude=latitude;
	}

public:
	void Set(double longitude,double latitude)
	{
		m_Longitude=longitude;
		m_Latitude=latitude;
	}

	MyLatLng& operator =( MyLatLng right)
	{
		this->m_Longitude = right.m_Longitude;
		this->m_Latitude = right.m_Latitude;
		return *this;
	}
};


//计算大圆弧对应的方位角
double getAngle(MyLatLng A,MyLatLng B)
{
	double result = 0.0;
	if((A.m_Latitude == B.m_Latitude)&&(A.m_Longitude == B.m_Longitude))
	{
		return result;
	}
	else if(A.m_Longitude == B.m_Longitude)
	{
		if(A.m_Latitude > B.m_Latitude)
			result = M_PI;
	}
	else
	{
		double pp = sin(B.m_Latitude) * sin(A.m_Latitude) + cos(B.m_Latitude) * cos(A.m_Latitude) * cos(B.m_Longitude - A.m_Longitude);
		if(pp > 1)
			pp = 1.0;
		if(pp < -1)
			pp = -1.0;
		double c = acos(pp);
		double qq = cos(B.m_Latitude) * sin(B.m_Longitude - A.m_Longitude) / sin(c);
		if(qq > 1)
			qq = 1.0;
		if(qq < -1)
			qq = -1.0;
		double Angle = asin(qq);
		result = Angle;
		if((B.m_Latitude > A.m_Latitude) && (B.m_Longitude > A.m_Longitude))
		{}
		else if((B.m_Latitude < A.m_Latitude) && (B.m_Longitude < A.m_Longitude))
		{
			result = M_PI - result;
		}
		else if((B.m_Latitude < A.m_Latitude) && (B.m_Longitude > A.m_Longitude))
		{
			result = M_PI - result;
		}
		else if((B.m_Latitude > A.m_Latitude) && (B.m_Longitude < A.m_Longitude))
		{
			result = result + 2 * M_PI;
		}
	}
	return result;
}


double distanceofpoint2plane(Vec3D point, Vec3D P , Vec3D Q)//P，Q为平面经过的两个点的坐标
{
	double A = P.y * Q.z - P.z * Q.y;
	double B = Q.x * P.z - P.x * Q.z;
	double C = P.x * Q.y - Q.x * P.y;
	double distance = abs(A * point.x + B * point.y + C * point.z) / sqrt(A * A + B * B + C * C);
	return distance;
}

//计算跨面两个单元与基础三角形边的交点
MyLatLng bigcirintersectpt(DgQ2DICoord& add1,DgQ2DICoord& endadd1,DgQ2DICoord& Dgtrigcintersectpointrc,const DgGeoCoord& centCoord1,const DgGeoCoord& centCoord2, const DgIDGG& dgg,DgGeoSphRF& geoRF)
{
	MyLatLng intersectcPt;
	GeoCoord add1trivertex0;
	GeoCoord add1trivertex1;
	if(endadd1.quadNum() - add1.quadNum() == 1)
	{
		if(add1.quadNum() >= 1 && add1.quadNum() < 5)
		{
			//三角形两个顶点
			add1trivertex0 = dgg.projTriRF().sphIcosa().sphIcosa().icotri[add1.quadNum() - 1][0];
			add1trivertex1 = dgg.projTriRF().sphIcosa().sphIcosa().icotri[add1.quadNum() - 1][2];
		}
		if(add1.quadNum() >= 6 && add1.quadNum() < 10)
		{
			//三角形两个顶点
			add1trivertex0 = dgg.projTriRF().sphIcosa().sphIcosa().icotri[add1.quadNum() + 9][0];
			add1trivertex1 = dgg.projTriRF().sphIcosa().sphIcosa().icotri[add1.quadNum() + 9][1];
		}

	}
	else if(endadd1.quadNum() - add1.quadNum() == -1)
	{
		if(add1.quadNum() > 1 && add1.quadNum() <= 5)
		{
			//三角形两个顶点
			add1trivertex0 = dgg.projTriRF().sphIcosa().sphIcosa().icotri[endadd1.quadNum() - 1][0];
			add1trivertex1 = dgg.projTriRF().sphIcosa().sphIcosa().icotri[endadd1.quadNum() - 1][2];
		}
		if(add1.quadNum() > 6 && add1.quadNum() <= 10)
		{
			//三角形两个顶点
			add1trivertex0 = dgg.projTriRF().sphIcosa().sphIcosa().icotri[endadd1.quadNum() + 9][0];
			add1trivertex1 = dgg.projTriRF().sphIcosa().sphIcosa().icotri[endadd1.quadNum() + 9][1];
		}
	}
	else if(endadd1.quadNum() - add1.quadNum() == 5)
	{
		if(add1.quadNum() >= 1 && add1.quadNum() <= 5)
		{
			//三角形两个顶点
			add1trivertex0 = dgg.projTriRF().sphIcosa().sphIcosa().icotri[add1.quadNum() + 4][0];
			add1trivertex1 = dgg.projTriRF().sphIcosa().sphIcosa().icotri[add1.quadNum() + 4][1];
		}

	}
	else if(endadd1.quadNum() - add1.quadNum() == -5)
	{
		if(add1.quadNum() >= 6 && add1.quadNum() <= 10)
		{
			//三角形两个顶点
			add1trivertex0 = dgg.projTriRF().sphIcosa().sphIcosa().icotri[endadd1.quadNum() + 4][0];
			add1trivertex1 = dgg.projTriRF().sphIcosa().sphIcosa().icotri[endadd1.quadNum() + 4][1];
		}
	}
	else if(endadd1.quadNum() - add1.quadNum() == 4)
	{
		if(add1.quadNum() > 1 && add1.quadNum() <= 5)
		{
			//三角形两个顶点
			add1trivertex0 = dgg.projTriRF().sphIcosa().sphIcosa().icotri[add1.quadNum() + 4][0];
			add1trivertex1 = dgg.projTriRF().sphIcosa().sphIcosa().icotri[add1.quadNum() + 4][2];
		}

		if(add1.quadNum() == 1)
		{
			//三角形两个顶点
			add1trivertex0 = dgg.projTriRF().sphIcosa().sphIcosa().icotri[add1.quadNum() - 1][0];
			add1trivertex1 = dgg.projTriRF().sphIcosa().sphIcosa().icotri[add1.quadNum() - 1][1];
		}
		if(add1.quadNum() == 6)
		{
			//三角形两个顶点
			add1trivertex0 = dgg.projTriRF().sphIcosa().sphIcosa().icotri[add1.quadNum() + 9][0];
			add1trivertex1 = dgg.projTriRF().sphIcosa().sphIcosa().icotri[add1.quadNum() + 9][2];
		}
	}
	else if(endadd1.quadNum() - add1.quadNum() == -4)
	{
		if(add1.quadNum() >= 6 && add1.quadNum() < 10)
		{
			//三角形两个顶点
			add1trivertex0 = dgg.projTriRF().sphIcosa().sphIcosa().icotri[endadd1.quadNum() + 4][0];
			add1trivertex1 = dgg.projTriRF().sphIcosa().sphIcosa().icotri[endadd1.quadNum() + 4][2];
		}

		if(add1.quadNum() == 5)
		{
			//三角形两个顶点
			add1trivertex0 = dgg.projTriRF().sphIcosa().sphIcosa().icotri[endadd1.quadNum() - 1][0];
			add1trivertex1 = dgg.projTriRF().sphIcosa().sphIcosa().icotri[endadd1.quadNum() - 1][1];
		}
		if(add1.quadNum() == 10)
		{
			//三角形两个顶点
			add1trivertex0 = dgg.projTriRF().sphIcosa().sphIcosa().icotri[endadd1.quadNum() + 9][0];
			add1trivertex1 = dgg.projTriRF().sphIcosa().sphIcosa().icotri[endadd1.quadNum() + 9][2];
		}
	}
	else if(endadd1.quadNum() - add1.quadNum() == 9)
	{
		if(add1.quadNum() == 1)
		{
			//三角形两个顶点
			add1trivertex0 = dgg.projTriRF().sphIcosa().sphIcosa().icotri[add1.quadNum() + 4][0];
			add1trivertex1 = dgg.projTriRF().sphIcosa().sphIcosa().icotri[add1.quadNum() + 4][2];
		}
	}
	else if(endadd1.quadNum() - add1.quadNum() == -9)
	{
		if(add1.quadNum() == 10)
		{
			//三角形两个顶点
			add1trivertex0 = dgg.projTriRF().sphIcosa().sphIcosa().icotri[endadd1.quadNum() + 4][0];
			add1trivertex1 = dgg.projTriRF().sphIcosa().sphIcosa().icotri[endadd1.quadNum() + 4][2];
		}
	}

	//跨边界的两个单元
	GeoCoord add1degree ,end1degree;
	add1degree.lon = centCoord1.lon();
	add1degree.lat = centCoord1.lat();
	end1degree.lat = centCoord2.lat();
	end1degree.lon = centCoord2.lon();
	//计算两个大圆弧交点
	GeoCoord trigcintersectpoint = GCintersect(add1trivertex0,add1trivertex1,add1degree,end1degree,1);
	DgGeoCoord Dgtrigcintersectpointadd(trigcintersectpoint.lon,trigcintersectpoint.lat,true);
	DgLocation* Dgtrigcintersectpoint = geoRF.makeLocation(Dgtrigcintersectpointadd);
	//两个大圆弧交点所在六边形单元的行列
	dgg.convert(Dgtrigcintersectpoint);
	const DgQ2DICoord* intersectpoint = dgg.getAddress(*Dgtrigcintersectpoint);
	Dgtrigcintersectpointrc = const_cast<DgQ2DICoord&>(*intersectpoint);
	//得出经纬度
	geoRF.convert(Dgtrigcintersectpoint);
	/* cout << "* and cell center point:" << *thePt2 << endl;*/
	const DgGeoCoord& centintersect = *geoRF.getAddress(*Dgtrigcintersectpoint);
	intersectcPt.Set(centintersect.lon(),centintersect.lat());
	return intersectcPt;

}
//八邻近
vector<DgQ2DICoord> neicell(DgQ2DICoord& add1,int maxI,int maxJ)
{
	vector<DgQ2DICoord> neiadd;
	//maxI和MaxJ 对于菱形 是相等的所以下面出现了混用的情况
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
			// neiadd.push_back(DgQ2DICoord(add1.quadNum() + 2, DgIVec2D(0,maxJ)));
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
			// neiadd.push_back(DgQ2DICoord(add1.quadNum() + 2, DgIVec2D(0,maxJ)));
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
			// neiadd.push_back(DgQ2DICoord(add1.quadNum() + 2, DgIVec2D(0,maxJ)));
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
			// neiadd.push_back(DgQ2DICoord(add1.quadNum() - 3, DgIVec2D(0,maxJ)));
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
			// neiadd.push_back(DgQ2DICoord(add1.quadNum() - 3, DgIVec2D(0,maxJ)));
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
			neiadd.push_back(DgQ2DICoord(add1.quadNum() - 4, DgIVec2D(0,0)));

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
			neiadd.push_back(DgQ2DICoord(add1.quadNum() - 9, DgIVec2D(0,0)));


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
			neiadd.push_back(DgQ2DICoord(add1.quadNum() + 1, DgIVec2D(0,0)));

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
			neiadd.push_back(DgQ2DICoord(add1.quadNum() - 4, DgIVec2D(0,0)));

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
			neiadd.push_back(DgQ2DICoord(add1.quadNum() - 4, DgIVec2D(maxI,0)));

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
			neiadd.push_back(DgQ2DICoord(add1.quadNum() - 9, DgIVec2D(maxI,0)));

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
			neiadd.push_back(DgQ2DICoord(add1.quadNum() + 4, DgIVec2D(maxI,maxJ)));

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
			neiadd.push_back(DgQ2DICoord(add1.quadNum() - 1, DgIVec2D(maxI,maxJ)));

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
			neiadd.push_back(DgQ2DICoord(add1.quadNum() - 5, DgIVec2D(maxI,0)));

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
			neiadd.push_back(DgQ2DICoord(add1.quadNum() - 5, DgIVec2D(maxI,0)));

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
				neiadd.push_back(DgQ2DICoord(add1.quadNum() + 9, DgIVec2D(maxI,maxJ)));

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
				neiadd.push_back(DgQ2DICoord(add1.quadNum() + 4, DgIVec2D(maxI,maxJ)));

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
				// neiadd.push_back(DgQ2DICoord(add1.quadNum() + 3, DgIVec2D(maxI,0))); 去除一个保持八邻近
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
				// neiadd.push_back(DgQ2DICoord(add1.quadNum() + 3, DgIVec2D(maxI,0))); 去除一个保持八邻近
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
				// neiadd.push_back(DgQ2DICoord(add1.quadNum() - 2, DgIVec2D(maxI,0)));去除一个保持八邻近
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
				// neiadd.push_back(DgQ2DICoord(add1.quadNum() - 2, DgIVec2D(maxI,0)));去除一个保持八邻近
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
				// neiadd.push_back(DgQ2DICoord(add1.quadNum() - 2, DgIVec2D(maxI,0))); 去除一个保持八邻近
				neiadd.push_back(DgQ2DICoord(add1.quadNum() - 3, DgIVec2D(maxI,0)));
			}
		}
	return neiadd;
}

//四邻近
vector<DgQ2DICoord> fourneicell(DgQ2DICoord& add1,int maxI,int maxJ)
{
	vector<DgQ2DICoord> neiadd;
	if(add1.coord().i() != 0 && add1.coord().i() != maxI &&
		add1.coord().j() != 0 && add1.coord().j() != maxJ)
	{
		neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() + 1, add1.coord().j())));
		neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i(), add1.coord().j() + 1)));
		neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() - 1, add1.coord().j())));
		neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i(), add1.coord().j() - 1)));

	}
	else if(add1.coord().i() == 0 && add1.coord().j() != 0 && add1.coord().j() != maxJ)
	{
		if(add1.quadNum() == 1)
		{
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() + 1, add1.coord().j())));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i(), add1.coord().j() + 1)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() + 4, DgIVec2D(maxJ - add1.coord().j(), maxJ)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i(), add1.coord().j() - 1)));

		}
		else if(add1.quadNum() > 1 && add1.quadNum() <= 5)
		{
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() + 1, add1.coord().j())));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i(), add1.coord().j() + 1)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() - 1, DgIVec2D(maxJ - add1.coord().j(), maxJ)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i(), add1.coord().j() - 1)));

		}
		else if(add1.quadNum() >= 6 && add1.quadNum() <= 10)
		{
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() + 1, add1.coord().j())));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i(), add1.coord().j() + 1)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() - 5, DgIVec2D(maxI,add1.coord().j())));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i(), add1.coord().j() - 1)));
		}

	}
	else if(add1.coord().i() == maxI && add1.coord().j() != 0 && add1.coord().j() != maxJ)
	{
		if(add1.quadNum() >= 1 && add1.quadNum() <= 5)
		{
			neiadd.push_back(DgQ2DICoord(add1.quadNum() + 5, DgIVec2D(0, add1.coord().j())));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i(), add1.coord().j() + 1)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() - 1, add1.coord().j())));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i(), add1.coord().j() - 1)));

		}
		else if(add1.quadNum() >= 6 && add1.quadNum() < 10)
		{
			neiadd.push_back(DgQ2DICoord(add1.quadNum() + 1, DgIVec2D(maxJ - add1.coord().j(), 0)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i(), add1.coord().j() + 1)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() - 1, add1.coord().j())));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i(), add1.coord().j() - 1)));

		}
		else if(add1.quadNum() == 10)
		{
			neiadd.push_back(DgQ2DICoord(add1.quadNum() - 4, DgIVec2D(maxJ - add1.coord().j(), 0)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i(), add1.coord().j() + 1)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() - 1, add1.coord().j())));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i(), add1.coord().j() - 1)));
		}
	}
	else if(add1.coord().j() == 0 && add1.coord().i() != 0 && add1.coord().i() != maxJ)
	{
		if(add1.quadNum() == 1)
		{
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() + 1, add1.coord().j())));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i(), add1.coord().j() + 1)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() - 1, add1.coord().j())));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() + 9, DgIVec2D(add1.coord().i(), maxJ)));

		}
		else if(add1.quadNum() > 1 && add1.quadNum() <= 5)
		{
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() + 1, add1.coord().j())));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i(), add1.coord().j() + 1)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() - 1, add1.coord().j())));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() + 4, DgIVec2D(add1.coord().i(), maxJ)));
		}
		else if(add1.quadNum() == 6)
		{
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() + 1, add1.coord().j())));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i(), add1.coord().j() + 1)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() - 1, add1.coord().j())));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() + 4, DgIVec2D(maxI, maxI - add1.coord().i())));

		}
		else if(add1.quadNum() > 6 && add1.quadNum() <= 10)
		{
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() + 1, add1.coord().j())));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i(), add1.coord().j() + 1)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() - 1, add1.coord().j())));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() - 1, DgIVec2D(maxI, maxI - add1.coord().i())));

		}
	}
	else if(add1.coord().j() == maxJ && add1.coord().i() != 0 && add1.coord().i() != maxJ)
	{
		if(add1.quadNum() >= 1 && add1.quadNum() < 5)
		{
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() + 1, add1.coord().j())));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() + 1, DgIVec2D(0, maxI - add1.coord().i())));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() - 1, add1.coord().j())));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i(), add1.coord().j() - 1)));

		}
		else if(add1.quadNum() == 5)
		{
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() + 1, add1.coord().j())));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() - 4, DgIVec2D(0, maxI - add1.coord().i())));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() - 1, add1.coord().j())));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i(), add1.coord().j() - 1)));

		}
		else if(add1.quadNum() >= 6 && add1.quadNum() < 10)
		{
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() + 1, add1.coord().j())));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() - 4, DgIVec2D(add1.coord().i(), 0)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() - 1, add1.coord().j())));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i(), add1.coord().j() - 1)));

		}
		else if(add1.quadNum() == 10)
		{
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() + 1, add1.coord().j())));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() - 9, DgIVec2D(add1.coord().i(), 0)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() - 1, add1.coord().j())));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i(), add1.coord().j() - 1)));

		}
	}
	else if(add1.coord().j() == maxJ && add1.coord().i() == 0)
	{
		if(add1.quadNum() == 1)
		{
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() + 1, add1.coord().j())));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() + 1, DgIVec2D(0,maxJ)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() + 4, DgIVec2D(0, maxJ)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i(), add1.coord().j() - 1)));
		}
		else if(add1.quadNum() == 2)
		{
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() + 1, add1.coord().j())));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() + 1, DgIVec2D(0,maxJ)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() - 1, DgIVec2D(0, maxJ)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i(), add1.coord().j() - 1)));

		}
		else if(add1.quadNum() == 3)
		{
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() + 1, add1.coord().j())));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() + 1, DgIVec2D(0,maxJ)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() - 1, DgIVec2D(0, maxJ)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i(), add1.coord().j() - 1)));

		}
		else if(add1.quadNum() == 4)
		{
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() + 1, add1.coord().j())));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() + 1, DgIVec2D(0,maxJ)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() - 1, DgIVec2D(0, maxJ)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i(), add1.coord().j() - 1)));
		}
		else if(add1.quadNum() == 5)
		{
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() + 1, add1.coord().j())));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() - 4, DgIVec2D(0,maxJ)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() - 1, DgIVec2D(0, maxJ)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i(), add1.coord().j() - 1)));

		}
		else if(add1.quadNum() >= 6 && add1.quadNum() < 10)
		{
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() + 1, add1.coord().j())));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() - 4, DgIVec2D(0,0)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() - 5, DgIVec2D(maxI, maxJ)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i(), add1.coord().j() - 1)));

		}
		else if(add1.quadNum() == 10)
		{
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() + 1, add1.coord().j())));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() + 9, DgIVec2D(0,0)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() - 5, DgIVec2D(maxI, maxJ)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i(), add1.coord().j() - 1)));

		}
	}
	else if(add1.coord().j() == maxJ && add1.coord().i() == maxI)
	{
		if(add1.quadNum() >= 1 && add1.quadNum() < 5)
		{
			neiadd.push_back(DgQ2DICoord(add1.quadNum() + 5, DgIVec2D(0, maxJ)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() + 1, DgIVec2D(0,0)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() - 1, add1.coord().j())));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i(), add1.coord().j() - 1)));

		}
		else if(add1.quadNum() == 5)
		{
			neiadd.push_back(DgQ2DICoord(add1.quadNum() + 5, DgIVec2D(0, maxJ)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() - 4, DgIVec2D(0,0)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() - 1, add1.coord().j())));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i(), add1.coord().j() - 1)));

		}
		else if(add1.quadNum() >= 6 && add1.quadNum() < 10)
		{
			neiadd.push_back(DgQ2DICoord(add1.quadNum() + 1, DgIVec2D(0,0)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() - 4, DgIVec2D(maxI,0)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() - 1, add1.coord().j())));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i(), add1.coord().j() - 1)));

		}
		else if(add1.quadNum() == 10)
		{
			neiadd.push_back(DgQ2DICoord(add1.quadNum() - 4, DgIVec2D(0,0)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() - 9, DgIVec2D(maxI,0)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() - 1, add1.coord().j())));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i(), add1.coord().j() - 1)));
		}
	}
	else if(add1.coord().j() == 0 && add1.coord().i() == 0)
	{
		if(add1.quadNum() == 1)
		{
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() + 1, add1.coord().j())));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i(), add1.coord().j() + 1)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() + 4, DgIVec2D(maxI,maxJ)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() + 9, DgIVec2D(0,maxJ)));
		}
		else if(add1.quadNum() > 1 && add1.quadNum() <= 5)
		{
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() + 1, add1.coord().j())));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i(), add1.coord().j() + 1)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() - 1, DgIVec2D(maxI,maxJ)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() + 4, DgIVec2D(0,maxJ)));

		}
		else if(add1.quadNum() == 6)
		{
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() + 1, add1.coord().j())));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i(), add1.coord().j() + 1)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() - 5, DgIVec2D(maxI,0)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() + 4, DgIVec2D(maxI,maxJ)));

		}
		else if(add1.quadNum() > 6 && add1.quadNum() <= 10)
		{
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() + 1, add1.coord().j())));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i(), add1.coord().j() + 1)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() - 5, DgIVec2D(maxI,0)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() - 1, DgIVec2D(maxI,maxJ)));
		}
		else if(add1.coord().j() == 0 && add1.coord().i() == maxI)
		{
			if(add1.quadNum() == 1)
			{
				neiadd.push_back(DgQ2DICoord(add1.quadNum() + 5, DgIVec2D(0,0)));
				neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i(), add1.coord().j() + 1)));
				neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() - 1, add1.coord().j())));
				neiadd.push_back(DgQ2DICoord(add1.quadNum() + 9, DgIVec2D(maxI,maxJ)));
			}
			else if(add1.quadNum() > 1 && add1.quadNum() <= 5)
			{
				neiadd.push_back(DgQ2DICoord(add1.quadNum() + 5, DgIVec2D(0,0)));
				neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i(), add1.coord().j() + 1)));
				neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() - 1, add1.coord().j())));
				neiadd.push_back(DgQ2DICoord(add1.quadNum() + 4, DgIVec2D(maxI,maxJ)));
			}
			else if(add1.quadNum() == 6)
			{
				neiadd.push_back(DgQ2DICoord(add1.quadNum() + 1, DgIVec2D(maxI,0)));
				neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i(), add1.coord().j() + 1)));
				neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() - 1, add1.coord().j())));
				neiadd.push_back(DgQ2DICoord(add1.quadNum() + 4, DgIVec2D(maxI,0)));
			}
			else if(add1.quadNum() == 7)
			{
				neiadd.push_back(DgQ2DICoord(add1.quadNum() + 1, DgIVec2D(maxI,0)));
				neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i(), add1.coord().j() + 1)));
				neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() - 1, add1.coord().j())));
				neiadd.push_back(DgQ2DICoord(add1.quadNum() - 1, DgIVec2D(maxI,0)));
			}
			else if(add1.quadNum() == 8)
			{
				neiadd.push_back(DgQ2DICoord(add1.quadNum() + 1, DgIVec2D(maxI,0)));
				neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i(), add1.coord().j() + 1)));
				neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() - 1, add1.coord().j())));
				neiadd.push_back(DgQ2DICoord(add1.quadNum() - 1, DgIVec2D(maxI,0)));

			}
			else if(add1.quadNum() == 9)
			{
				neiadd.push_back(DgQ2DICoord(add1.quadNum() + 1, DgIVec2D(maxI,0)));
				neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i(), add1.coord().j() + 1)));
				neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() - 1, add1.coord().j())));
				neiadd.push_back(DgQ2DICoord(add1.quadNum() - 1, DgIVec2D(maxI,0)));
			}
			else if(add1.quadNum() == 10)
			{
				neiadd.push_back(DgQ2DICoord(add1.quadNum() - 4, DgIVec2D(maxI,0)));
				neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i(), add1.coord().j() + 1)));
				neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() - 1, add1.coord().j())));
				neiadd.push_back(DgQ2DICoord(add1.quadNum() - 1, DgIVec2D(maxI,0)));
			}
		}
	}

	return neiadd;
}
vector<polygon> discreteline(vector<vector<MyLatLng>> LinePoints)
{
	// create the reference frame (RF) conversion network
	DgRFNetwork net0;

	// create the geodetic reference frame
	DgGeoSphRF geoRF(net0, "GS0");

	// create the ISEA4H grid system with resolutions 0-9; requires a
	// fixed icosahedron vertex and edge azimuth
	/*DgGeoCoord vert0(11.25L, 58.28252559L, false); */// args: lon, lat, isRadians
	DgGeoCoord vert0(112L, 33L, false);
	/*DgGeoCoord vert0(123.2628L, 35.7635L, false);*/
	long double azimuth = 0.0L;
	const DgIDGGS4D idggs(net0, geoRF, vert0, azimuth, 25); // last argument is
	// number of resolutions

	// get the resolution 7 dgg from the dggs
	const DgIDGG& dgg = idggs.idgg(9);

	int a = dgg.maxI();
	int b = dgg.maxJ();
	cout<<a<<endl;
	cout<<b<<endl;
	//总的离散线单元
	vector<DgQ2DICoord> lisanxianadds;

	for(int j = 0;j < LinePoints.size();j++)
	{

		for(int i = 0;i < LinePoints[j].size()-1;i++)
		{
			//cout<<LinePoints[j][i].m_Longitude<<" "<<LinePoints[j][i].m_Latitude<<endl;
			DgGeoCoord geoAddress(LinePoints[j][i].m_Longitude, LinePoints[j][i].m_Latitude, false);
			//DgGeoCoord geoAddress(112.06122998404138, 35.28261899, false);//112.06122998404138, 35.282618990033619

			DgLocation* thePt1 = geoRF.makeLocation(geoAddress);
			//得出行列
			dgg.convert(thePt1);
			const DgQ2DICoord* add = dgg.getAddress(*thePt1 );
			DgQ2DICoord add1 = const_cast<DgQ2DICoord&>(*add);

			//得出经纬度
			geoRF.convert(thePt1);
			const DgGeoCoord& centCoord1 = *geoRF.getAddress(*thePt1);
			/*MyLatLng firstcPt(centCoord1.lonDegs(),centCoord1.latDegs());*/
			MyLatLng firstcPtRad(centCoord1.lon(),centCoord1.lat());
			//cout<<LinePoints[j][i + 1].m_Longitude<<" "<<LinePoints[j][i + 1].m_Latitude<<endl;
			DgGeoCoord geoAddress1(LinePoints[j][i + 1].m_Longitude, LinePoints[j][i + 1].m_Latitude, false);
			//cout<<geoAddress1.lon()<<" "<<geoAddress1.lat()<<endl;
			DgLocation* thePt2 = geoRF.makeLocation(geoAddress1);
			//cout << "the point " << *thePt2 << endl;
			//得出行列
			dgg.convert(thePt2);
			//cout << "* lies in cell " << *thePt2 << endl;

			const DgQ2DICoord* endadd = dgg.getAddress(*thePt2 );
			DgQ2DICoord endadd1 = const_cast<DgQ2DICoord&>(*endadd);
			//得出经纬度
			geoRF.convert(thePt2);
			//cout << "* and cell center point:" << *thePt2 << endl;
			const DgGeoCoord& centCoord2 = *geoRF.getAddress(*thePt2);
			/*MyLatLng endcPt(centCoord2.lonDegs(),centCoord2.latDegs());*/
			MyLatLng endcPtRad(centCoord2.lon(),centCoord2.lat());

			double firstendPt = getAngle(firstcPtRad, endcPtRad);  //算出起点和终点的方位角
			//double firstendPt = getAngle(firstcPt, endcPt);  //算出起点和终点的方位角

			GeoCoord fcpt;
			GeoCoord ecpt;
			fcpt.lon = firstcPtRad.m_Longitude;
			fcpt.lat = firstcPtRad.m_Latitude;
			ecpt.lon = endcPtRad.m_Longitude;
			ecpt.lat = endcPtRad.m_Latitude;
			Vec3D firstcPtxyz = llxyz(fcpt);
			Vec3D endcPtxyz = llxyz(ecpt);

			vector<DgQ2DICoord> neiadd;
			/*vector<MyLatLng> neiaddlonlat;*/
			vector<MyLatLng> neiaddlonlatRad;
			vector<double> neiaddAnglepaixu;
			vector<double> neiaddangle;
			vector<int> bestcode;
			//离散线单元的行列
			vector<DgQ2DICoord> lisanxianadd;
			DgQ2DICoord Dgtrigcintersectpointrc;
			DgQ2DICoord Dgtrigcintersectpointrc1;
			DgQ2DICoord Dgtrigcintersectpointrc2;

			if(add1.quadNum() != endadd1.quadNum())
			{
				MyLatLng intersectcPt;
				MyLatLng intersectcPt1;
				MyLatLng intersectcPt2;
				intersectcPt = bigcirintersectpt(add1,endadd1,Dgtrigcintersectpointrc,geoAddress,geoAddress1,dgg,geoRF);

				if(add1.quadNum() != Dgtrigcintersectpointrc.quadNum())
				{
					vector<DgQ2DICoord> neiinterptrc = fourneicell(Dgtrigcintersectpointrc,a,b);
					for(int m = 0; m < neiinterptrc.size();m++)
					{
						if(add1.quadNum() == neiinterptrc[m].quadNum())
						{
							Dgtrigcintersectpointrc1 = neiinterptrc[m];
							Dgtrigcintersectpointrc2 = Dgtrigcintersectpointrc;
							DgLocation* neiintercenadd =dgg.makeLocation(Dgtrigcintersectpointrc1) ;
							geoRF.convert(neiintercenadd);
							const DgGeoCoord& DGeoneiintercenter = *geoRF.getAddress(*neiintercenadd);

							intersectcPt1 = MyLatLng(DGeoneiintercenter.lon(),DGeoneiintercenter.lat());
							intersectcPt2 = intersectcPt;
						}
					}
				}
				else
				{
					vector<DgQ2DICoord> neiinterptrc = fourneicell(Dgtrigcintersectpointrc,a,b);
					for(int m = 0; m < neiinterptrc.size();m++)
					{
						if(endadd1.quadNum() == neiinterptrc[m].quadNum())
						{
							Dgtrigcintersectpointrc1 = Dgtrigcintersectpointrc;
							Dgtrigcintersectpointrc2 = neiinterptrc[m];
							DgLocation* neiintercenadd =dgg.makeLocation(Dgtrigcintersectpointrc2) ;
							geoRF.convert(neiintercenadd);
							const DgGeoCoord& DGeoneiintercenter = *geoRF.getAddress(*neiintercenadd);

							intersectcPt2 = MyLatLng(DGeoneiintercenter.lon(),DGeoneiintercenter.lat());
							intersectcPt1 = intersectcPt;
						}
					}
				}


				GeoCoord incpt1;

				incpt1.lon = intersectcPt1.m_Longitude;
				incpt1.lat = intersectcPt1.m_Latitude;

				Vec3D intersectcPtxyz1 = llxyz(incpt1);

				GeoCoord incpt2;

				incpt2.lon = intersectcPt2.m_Longitude;
				incpt2.lat = intersectcPt2.m_Latitude;

				Vec3D intersectcPtxyz2 = llxyz(incpt2);

				//菱形的八个边邻近单元
				neiadd = neicell(add1,a,b);
				double add1intersectPt1 =  getAngle(firstcPtRad,intersectcPt1);
				//由行列转经纬度

				for(int a = 0;a < neiadd.size();a++)
				{
					DgLocation* neicenadd =dgg.makeLocation(neiadd[a]) ;
					geoRF.convert(neicenadd);
					const DgGeoCoord& DGeoneicenter = *geoRF.getAddress(*neicenadd);
					/*MyLatLng Gneicenter(DGeoneicenter.lonDegs(),DGeoneicenter.latDegs());*/
					MyLatLng GneicenterRad(DGeoneicenter.lon(),DGeoneicenter.lat());
					neiaddlonlatRad.push_back(GneicenterRad);
					/*neiaddlonlat.push_back(Gneicenter);*/
				}



				vector<double> copyneiaddangle;
				//得到八个相邻单元与AB单元中心连接的大圆弧之间的夹角
				for(int b = 0;b < neiadd.size();b++)
				{
					double angle = abs(getAngle(firstcPtRad , neiaddlonlatRad[b]) - add1intersectPt1);
					//cout<<angle<<endl;
					if(angle > M_PI)
						angle = abs(2 * M_PI - angle);
					//cout<<angle<<endl;

					neiaddangle.push_back(angle);
					copyneiaddangle.push_back(angle);

				}

				sort(copyneiaddangle.begin(),copyneiaddangle.end());

				for(int c = 0;c < neiaddangle.size();c++)
				{


					if(neiaddangle[c] == copyneiaddangle[0])
					{

						bestcode.push_back(c);    
						///*cout<<i<<endl;*/
					}

					if(neiaddangle[c] == copyneiaddangle[1])
					{

						bestcode.push_back(c);   
						/*cout<<i<<endl;*/
					}

					if(neiaddangle[c] == copyneiaddangle[2])
					{

						bestcode.push_back(c);    
						/*cout<<i<<endl;*/
					}
					if(neiaddangle[c] == copyneiaddangle[3])
					{

						bestcode.push_back(c);    
						///*cout<<i<<endl;*/
					}



				}

				neiadd.clear();
				neiaddlonlatRad.clear();
				while (add1 != Dgtrigcintersectpointrc1)
				{
					/* lisanxianadd.push_back(add1);*/
					lisanxianadds.push_back(add1);

					neiadd = neicell(add1,a,b);

					//由行列转经纬度

					for(int d = 0;d < neiadd.size();d++)
					{
						DgLocation* neicenadd =dgg.makeLocation(neiadd[d]) ;
						geoRF.convert(neicenadd);
						const DgGeoCoord& DGeoneicenter = *geoRF.getAddress(*neicenadd);
						/*MyLatLng Gneicenter(DGeoneicenter.lonDegs(),DGeoneicenter.latDegs());*/
						MyLatLng GneicenterRad(DGeoneicenter.lon(),DGeoneicenter.lat());
						neiaddlonlatRad.push_back(GneicenterRad);
						/*neiaddlonlat.push_back(Gneicenter);*/
					}
					//判断最小球面距离单元




					vector<Vec3D> Geoneiaddxyz;
					for(int k = 0;k < bestcode.size();k++)
					{
						GeoCoord lonlatPt;
						Vec3D xyzPt;

						lonlatPt.lon = neiaddlonlatRad[bestcode[k]].m_Longitude;
						lonlatPt.lat = neiaddlonlatRad[bestcode[k]].m_Latitude;
						xyzPt = llxyz(lonlatPt);
						Geoneiaddxyz.push_back(xyzPt);

					}


					vector<double> distance;
					vector<double> copydistance;
					for(int r = 0;r < bestcode.size();r++)
					{
						double dist = distanceofpoint2plane(Geoneiaddxyz[r],firstcPtxyz,intersectcPtxyz1);
						distance.push_back(dist);
						/*copydistance.push_back(dist);*/
					}



					if(distance[0] < distance[1] && distance[0] < distance[2] && distance[0] < distance[3])
					{
						add1 =  neiadd[bestcode[0]];
					}
					if(distance[1] < distance[0] && distance[1] < distance[2] && distance[1] < distance[3])
					{
						add1 =  neiadd[bestcode[1]];
					}
					if(distance[2] < distance[0] && distance[2] < distance[1] && distance[2] < distance[3])
					{
						add1 =  neiadd[bestcode[2]];
					}

					if(distance[3] < distance[0] && distance[3] < distance[1] && distance[3] < distance[2])
					{
						add1 =  neiadd[bestcode[3]];
					}

					neiadd.clear();
					neiaddlonlatRad.clear();
					neiaddangle.clear();
					Geoneiaddxyz.clear();

				}//离散线两个点之间的格网单元
				//         /* lisanxianadd.push_back(endadd1);*/

				bestcode.clear();
				copyneiaddangle.clear();
				add1 = Dgtrigcintersectpointrc2;

				neiadd = neicell(add1,a,b);
				double intersectendPt =  getAngle(intersectcPt2,endcPtRad);	



				//由行列转经纬度

				for(int a = 0;a < neiadd.size();a++)
				{
					DgLocation* neicenadd =dgg.makeLocation(neiadd[a]) ;
					geoRF.convert(neicenadd);
					const DgGeoCoord& DGeoneicenter = *geoRF.getAddress(*neicenadd);
					/*MyLatLng Gneicenter(DGeoneicenter.lonDegs(),DGeoneicenter.latDegs());*/
					MyLatLng GneicenterRad(DGeoneicenter.lon(),DGeoneicenter.lat());
					neiaddlonlatRad.push_back(GneicenterRad);
				}




				//得到六个相邻单元与AB单元中心连接的大圆弧之间的夹角
				for(int b = 0;b < neiadd.size();b++)
				{
					double angle = abs(getAngle(intersectcPt2, neiaddlonlatRad[b]) - intersectendPt);
					//cout<<angle<<endl;
					if(angle > M_PI)
						angle = abs(2 * M_PI - angle);
					//cout<<angle<<endl;

					neiaddangle.push_back(angle);
					copyneiaddangle.push_back(angle);

				}

				sort(copyneiaddangle.begin(),copyneiaddangle.end());

				for(int c = 0;c < neiaddangle.size();c++)
				{


					if(neiaddangle[c] == copyneiaddangle[0])
					{

						bestcode.push_back(c);    
						///*cout<<i<<endl;*/
					}

					if(neiaddangle[c] == copyneiaddangle[1])
					{

						bestcode.push_back(c);   
						/*cout<<i<<endl;*/
					}

					if(neiaddangle[c] == copyneiaddangle[2])
					{

						bestcode.push_back(c);    
						/*cout<<i<<endl;*/
					}
					if(neiaddangle[c] == copyneiaddangle[3])
					{

						bestcode.push_back(c);    
						///*cout<<i<<endl;*/
					}



				}

				neiadd.clear();
				neiaddlonlatRad.clear();
				while (add1 != endadd1)
				{
					/* lisanxianadd.push_back(add1);*/
					lisanxianadds.push_back(add1);

					neiadd = neicell(add1,a,b);

					//由行列转经纬度

					for(int d = 0;d < neiadd.size();d++)
					{
						DgLocation* neicenadd =dgg.makeLocation(neiadd[d]) ;
						geoRF.convert(neicenadd);
						const DgGeoCoord& DGeoneicenter = *geoRF.getAddress(*neicenadd);

						MyLatLng GneicenterRad(DGeoneicenter.lon(),DGeoneicenter.lat());
						neiaddlonlatRad.push_back(GneicenterRad);
						/*MyLatLng Gneicenter(DGeoneicenter.lonDegs(),DGeoneicenter.latDegs());

						neiaddlonlat.push_back(Gneicenter);*/
					}



					vector<Vec3D> Geoneiaddxyz;
					for(int k = 0;k < bestcode.size();k++)
					{
						GeoCoord lonlatPt;
						Vec3D xyzPt;

						lonlatPt.lon = neiaddlonlatRad[bestcode[k]].m_Longitude;
						lonlatPt.lat = neiaddlonlatRad[bestcode[k]].m_Latitude;
						xyzPt = llxyz(lonlatPt);
						Geoneiaddxyz.push_back(xyzPt);

					}


					vector<double> distance;
					vector<double> copydistance;
					for(int r = 0;r < bestcode.size();r++)
					{
						double dist = distanceofpoint2plane(Geoneiaddxyz[r],intersectcPtxyz2,endcPtxyz);
						distance.push_back(dist);
						/*copydistance.push_back(dist);*/
					}



					if(distance[0] < distance[1] && distance[0] < distance[2] && distance[0] < distance[3])
					{
						add1 =  neiadd[bestcode[0]];
					}
					if(distance[1] < distance[0] && distance[1] < distance[2] && distance[1] < distance[3])
					{
						add1 =  neiadd[bestcode[1]];
					}
					if(distance[2] < distance[0] && distance[2] < distance[1] && distance[2] < distance[3])
					{
						add1 =  neiadd[bestcode[2]];
					}

					if(distance[3] < distance[0] && distance[3] < distance[1] && distance[3] < distance[2])
					{
						add1 =  neiadd[bestcode[3]];
					}

					neiadd.clear();
					neiaddlonlatRad.clear();
					neiaddangle.clear();
					Geoneiaddxyz.clear();

				}//离散线两个点之间的格网单元


				lisanxianadds.push_back(endadd1);


			}
			else
			{
				//菱形的八个边邻近单元
				neiadd = neicell(add1,a,b);

				//由行列转经纬度

				for(int a = 0;a < neiadd.size();a++)
				{
					DgLocation* neicenadd =dgg.makeLocation(neiadd[a]) ;
					geoRF.convert(neicenadd);
					const DgGeoCoord& DGeoneicenter = *geoRF.getAddress(*neicenadd);
					/*MyLatLng Gneicenter(DGeoneicenter.lonDegs(),DGeoneicenter.latDegs());*/
					MyLatLng GneicenterRad(DGeoneicenter.lon(),DGeoneicenter.lat());
					neiaddlonlatRad.push_back(GneicenterRad);
				}



				vector<double> copyneiaddangle;
				//得到六个相邻单元与AB单元中心连接的大圆弧之间的夹角
				for(int b = 0;b < neiadd.size();b++)
				{
					double angle = abs(getAngle(firstcPtRad, neiaddlonlatRad[b]) - firstendPt);
					//cout<<angle<<endl;
					if(angle > M_PI)
						angle = abs(2 * M_PI - angle);
					//cout<<angle<<endl;

					neiaddangle.push_back(angle);
					copyneiaddangle.push_back(angle);

				}

				sort(copyneiaddangle.begin(),copyneiaddangle.end());

				for(int c = 0;c < neiaddangle.size();c++)
				{


					if(neiaddangle[c] == copyneiaddangle[0])
					{

						bestcode.push_back(c);    
						///*cout<<i<<endl;*/
					}

					if(neiaddangle[c] == copyneiaddangle[1])
					{

						bestcode.push_back(c);   
						/*cout<<i<<endl;*/
					}

					if(neiaddangle[c] == copyneiaddangle[2])
					{

						bestcode.push_back(c);    
						/*cout<<i<<endl;*/
					}
					if(neiaddangle[c] == copyneiaddangle[3])
					{

						bestcode.push_back(c);    
						///*cout<<i<<endl;*/
					}



				}

				neiadd.clear();
				neiaddlonlatRad.clear();
				while (add1 != endadd1)
				{
					/* lisanxianadd.push_back(add1);*/
					lisanxianadds.push_back(add1);

					neiadd = neicell(add1,a,b);

					//由行列转经纬度

					for(int d = 0;d < neiadd.size();d++)
					{
						DgLocation* neicenadd =dgg.makeLocation(neiadd[d]) ;
						geoRF.convert(neicenadd);
						const DgGeoCoord& DGeoneicenter = *geoRF.getAddress(*neicenadd);

						MyLatLng GneicenterRad(DGeoneicenter.lon(),DGeoneicenter.lat());
						neiaddlonlatRad.push_back(GneicenterRad);
						/*MyLatLng Gneicenter(DGeoneicenter.lonDegs(),DGeoneicenter.latDegs());

						neiaddlonlat.push_back(Gneicenter);*/
					}



					vector<Vec3D> Geoneiaddxyz;
					for(int k = 0;k < bestcode.size();k++)
					{
						GeoCoord lonlatPt;
						Vec3D xyzPt;

						lonlatPt.lon = neiaddlonlatRad[bestcode[k]].m_Longitude;
						lonlatPt.lat = neiaddlonlatRad[bestcode[k]].m_Latitude;
						xyzPt = llxyz(lonlatPt);
						Geoneiaddxyz.push_back(xyzPt);

					}


					vector<double> distance;
					vector<double> copydistance;
					for(int r = 0;r < bestcode.size();r++)
					{
						double dist = distanceofpoint2plane(Geoneiaddxyz[r],firstcPtxyz,endcPtxyz);
						distance.push_back(dist);
						/*copydistance.push_back(dist);*/
					}



					if(distance[0] < distance[1] && distance[0] < distance[2] && distance[0] < distance[3])
					{
						add1 =  neiadd[bestcode[0]];
					}
					if(distance[1] < distance[0] && distance[1] < distance[2] && distance[1] < distance[3])
					{
						add1 =  neiadd[bestcode[1]];
					}
					if(distance[2] < distance[0] && distance[2] < distance[1] && distance[2] < distance[3])
					{
						add1 =  neiadd[bestcode[2]];
					}

					if(distance[3] < distance[0] && distance[3] < distance[1] && distance[3] < distance[2])
					{
						add1 =  neiadd[bestcode[3]];
					}

					neiadd.clear();
					neiaddlonlatRad.clear();
					neiaddangle.clear();
					Geoneiaddxyz.clear();

				}//离散线两个点之间的格网单元
				//         /* lisanxianadd.push_back(endadd1);*/
				lisanxianadds.push_back(endadd1);
			}

		}//一条线上的格网单元
	}//所有线上的格网单元

	vector<polygon>lisanxiancells;
	DgPolygon verts;
	int ptsPerEdgeDensify = 0;
	for(int n = 0;n < lisanxianadds.size();n++)
	{

		DgLocation* lisanxiancenadd =dgg.bndRF().discRF().makeLocation(lisanxianadds[n]) ;
		dgg.setVertices(*lisanxiancenadd, verts, ptsPerEdgeDensify);
		polygon cell;
		for(int r = 0;r < verts.size();r++)
		{
			const DgGeoCoord& Vert = *geoRF.getAddress(verts[r]);
			GeoCoord sv;

			sv.lat = Vert.latDegs();
			sv.lon = Vert.lonDegs();
			cell.Vertices.push_back(sv);
		}
		lisanxiancells.push_back(cell);


	}

	return lisanxiancells;
}

OGRGeometry* CreateGeometry(polygon cell)
{
	OGRPolygon *pGeometry = (OGRPolygon*)OGRGeometryFactory::createGeometry(wkbPolygon);
	OGRLinearRing *pRing = (OGRLinearRing *)OGRGeometryFactory::createGeometry(wkbLinearRing);
	OGRPoint pt;

	double R = 6371007.180918475;
	for(int i = 0;i < cell.Vertices.size();i++)
	{
		double x = cell.Vertices[i].lon;
		double y = cell.Vertices[i].lat;


		pt.setX(x);
		pt.setY(y);

		pRing->addPoint(&pt);
	}

	pt.setX(cell.Vertices[0].lon);
	pt.setY(cell.Vertices[0].lat);

	pRing->addPoint(&pt);
	pRing->closeRings();
	pGeometry->addRing(pRing);
	OGRGeometry* pGeo = (OGRGeometry*)pGeometry;
	return pGeo;
}

int main(int argc, char* argv[])
{
	GDALAllRegister();
	GDALDataset   *poDS1;

	//为了支持中文路径，请添加下面这句代码
	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8","NO");
	CPLSetConfigOption("SHAPE_ENCODING","");  //解决中文乱码问题
	//读取shp文件
	poDS1 = (GDALDataset*) GDALOpenEx("henanline.shp", GDAL_OF_VECTOR, NULL, NULL, NULL );
	if( poDS1 == NULL )
	{
		printf( "Open failed.\n%s" );
		return 0;
	}

	OGRLayer  *poLayer1;
	poLayer1 = poDS1->GetLayer(0); //读取层
	OGRFeature *poFeature1;
	int iFeatureCount1=poLayer1 -> GetFeatureCount();   //获得要素个数
	cout<<"要素个数"<<iFeatureCount1<<endl;
	vector<MyLatLng> shpPoint;                                     //用于存点数据
	vector<vector<MyLatLng>> shpPointS;
	poLayer1 -> ResetReading();
	int i=0;


	while( (poFeature1 = poLayer1->GetNextFeature()) != NULL )
	{

		///获取要素中的几何体
		OGRGeometry* poGeometry = poFeature1->GetGeometryRef();
		OGRwkbGeometryType geotype;
		geotype = poGeometry -> getGeometryType();        //获取该文件类型，点状、线状or面状 

		OGRLineString *poLine = NULL;
		poLine = (OGRLineString*)poGeometry;

		int pointnums = poLine->getNumPoints();
		for(int i = 0;i < pointnums;i++)
		{
			MyLatLng Point(poLine -> getX(i),poLine -> getY(i));

			shpPoint.push_back(Point);

		}
		shpPointS.push_back(shpPoint);
		shpPoint.clear();
	}

	GDALClose( poDS1 );

	vector<polygon> lsxcells;
	lsxcells = discreteline(shpPointS);

	OGRSpatialReference spatialReference;
	spatialReference.SetWellKnownGeogCS("WGS84");
	const char *pszDriverName = "ESRI Shapefile";
	GDALDriver *poDriver;
	poDriver = GetGDALDriverManager()->GetDriverByName(pszDriverName );
	if( poDriver == NULL )
	{
		printf( "%s driver not available.\n", pszDriverName ); 
		return 0; 
	}
	GDALDataset *poDS;

	poDS = poDriver->Create("lisanxianoutShp9.shp", 0, 0, 0, GDT_Unknown, NULL ); //创建shp文件
	if( poDS == NULL )
	{
		printf( "Creation of output file failed.\n" ); 
		return 0;
	}
	OGRLayer *poLayer;
	poLayer = poDS->CreateLayer( "testpolygon", &spatialReference, wkbPolygon, NULL );
	if( poLayer == NULL )
	{ 
		printf( "Layer creation failed.\n" );
		return 0;
	}


	//下面创建属性表
	//先创建一个叫FieldID整形属性
	OGRFieldDefn oFieldID("FieldID",OFTInteger);
	poLayer->CreateField(&oFieldID);

	//再创建一个叫FeatureName的字符型属性，字符长度为
	OGRFieldDefn oFieldName("FieldName",OFTString);
	oFieldName.SetWidth(100);
	poLayer->CreateField(&oFieldName);


	OGRFeatureDefn *poDefn = poLayer->GetLayerDefn();


	for(int i = 0;i < lsxcells.size();i++)
	{

		//创建菱形
		OGRFeature *poFeaturePentagon = OGRFeature::CreateFeature(poDefn);
		poFeaturePentagon->SetField(0,2);
		poFeaturePentagon->SetField(1,"菱形");
		poFeaturePentagon->SetGeometry(CreateGeometry(lsxcells[i]));
		poLayer->CreateFeature(poFeaturePentagon);
		OGRFeature::DestroyFeature(poFeaturePentagon);

	}

	/*OGRDataSource::DestroyDataSource(poDS);*/
	printf("\n数据集创建完成！\n");
	return 0;
}	 

