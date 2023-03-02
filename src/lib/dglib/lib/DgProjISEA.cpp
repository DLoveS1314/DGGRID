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
// DgProjISEA.cpp: DgProjISEA class implementation
//
////////////////////////////////////////////////////////////////////////////////

#include <cmath>
#include <climits>

#include <dglib/DgProjISEA.h>

////////////////////////////////////////////////////////////////////////////////
DgProjISEAInv::DgProjISEAInv (const DgRF<DgProjTriCoord, long double>& from,
                       const DgRF<DgGeoCoord, long double>& to)
         : DgConverter<DgProjTriCoord, long double, DgGeoCoord, long double>(from, to),
           pProjTriRF_ (0)
{
   pProjTriRF_ = dynamic_cast<const DgProjTriRF*>(&fromFrame());

   if (!pProjTriRF_)
   {
      report("DgProjISEAInv::DgProjISEAInv(): "
        " fromFrame not of type DgProjTriRF", DgBase::Fatal);
   }

} // DgProjISEAInv::DgProjISEAInv

////////////////////////////////////////////////////////////////////////////////
DgGeoCoord
DgProjISEAInv::convertTypedAddress (const DgProjTriCoord& addIn) const
{
//cout << "***DgProjISEAInv: DgProjTriCoord: " << addIn << endl;
   IcosaGridPt gridpt;
   gridpt.pt.x = addIn.coord().x();
   gridpt.pt.y = addIn.coord().y();
   gridpt.triangle = addIn.triNum();

//cout << "    gridpt.triangle .x .y: " << gridpt.triangle << ", " <<
//      gridpt.pt.x << ", " << gridpt.pt.y << endl;

   GeoCoord ll = snyderInv(gridpt, projTriRF().sphIcosa().sphIcosa());

//cout << " ll.lon, ll.lat: " << ll.lon << ", " <<
//ll.lat << endl;
   DgGeoCoord geoPt(ll.lon, ll.lat);
   geoPt.normalize();

//cout << "    geoPt: " << geoPt << endl;
   return geoPt;

} // DgGeoCoord DgProjISEAInv::convertTypedAddress

//////////////////////////////////isea前向 属于是 球面投到正多面体//////////////////////////////////////////////
DgProjISEAFwd::DgProjISEAFwd (const DgRF<DgGeoCoord, long double>& from,
                    const DgRF<DgProjTriCoord, long double>& to)
         : DgConverter<DgGeoCoord, long double, DgProjTriCoord, long double>(from, to)
{
//    DgRF<DgProjTriCoord, long double> 的子类有一个是DgProjTriRF
   pProjTriRF_= dynamic_cast<const DgProjTriRF*>(&toFrame());

   if (!pProjTriRF_)
   {
      report("DgProjISEAFwd::DgProjISEAFwd(): "
        " toFrame not of type DgProjTriRF", DgBase::Fatal);
   }

} // DgProjISEAFwd::DgProjISEAFwd

///////////////////////////////////实现基类的虚函数完成地址转换/////////////////////////////////////////////
DgProjTriCoord
DgProjISEAFwd::convertTypedAddress (const DgGeoCoord& addIn) const
{

//cout << "***DgProjISEAFwd: geoPt: " << addIn << endl;
   GeoCoord ll;

   ll.lon = addIn.lon();
   ll.lat = addIn.lat();

//cout << "   ll.lon, ll.lat: " << ll.lon << ", " << ll.lat << endl;
//前向计算 是一个c语言形式的函数
   IcosaGridPt gridpt = snyderFwd(ll, projTriRF().sphIcosa());
//cout << "    gridpt.triangle .x .y: " << gridpt.triangle << ", " <<
//gridpt.pt.x << ", " << gridpt.pt.y << endl;

//cout << "DgProjTriCoord: " << DgProjTriCoord(gridpt.triangle,
//                               DgDVec2D(gridpt.pt.x, gridpt.pt.y)) << endl;

   return DgProjTriCoord(gridpt.triangle, DgDVec2D(gridpt.pt.x, gridpt.pt.y));

} // DgProjTriCoord DgProjISEAFwd::convertTypedAddress

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/*
   The C++ methods above are wrappers for the C functions below, which were
   written by Lian Song and Kevin Sahr.
*/
////////////////////////////////////////////////////////////////////////////////
/**
 * 投影相关的常数
 *
 * 根据投影条件A, 有(1)式=(3)式，因此：
 * \f$ \frac{(G-\theta) \pi R^2}{180} = \frac{1}{2} (R' \tan g)^2 \sin\theta \cos\theta \f$，
 * 所以\f$ R'= \frac{2R}{\tan g} \sqrt{\frac{\pi}{15 \sqrt{3}}} \f$，将二十面体的G（36.0度）、
 * \f$\theta\f$（37.37736814度）带入即得\f$R'\f$=0.9103832814.
 */
static const long double R1 = 0.9103832815L;
static const long double R1S = R1 * R1;

static const long double DH    = 37.37736814L * M_PI_180;
static const long double GH    = 36.0L * M_PI_180;
static const long double cot30 = 1.0L / tanl(30.0L * M_PI_180);
static const long double tanDH = tanl(DH);
static const long double cosDH = cosl(DH);
static const long double sinGH = sinl(GH);
static const long double cosGH = cosl(GH);
/**
 * 设置坐标原点为二十面体三角面中心，X坐标的偏移量。
 *
 * 假设球体半径为1（单位球），则与之等积的正二十面体的边长\f$ a = \frac{2}{15}\sqrt{15 \sqrt{3} \pi}\f$，（求a的公式在贲进教员文章中第7页表1中）
 * 即1.204591006。若以三角面中心为原点，X轴平移量为\f$ \frac{1}{2}a \f$，即0.6022955029。
 */
static const long double originXOff = 0.6022955029L;
/*开始是以正二十面体的某个顶点的位置来确定正而是面体的，现在将坐标的原点移动到了三角形的中心，此处的x，y偏移量即为二者相差*/
/**
 * 设置坐标原点为二十面体三角面中心，Y坐标的偏移量。
 *
 * 假设球体半径为1（单位球），则与之等积的正二十面体的边长\f$ a = \frac{2}{15}\sqrt{15 \sqrt{3} \pi}\f$，
 * 即1.204591006。若以三角面中心为原点，Y轴平移量为\f$ \frac{\sqrt{3}}{6}a \f$，即0.3477354707。
 */
static const long double originYOff = 0.3477354707L;
/** 二十面体展开后的实际边长
 *
 * 由等积投影条件a得到，严格表达式为：
 * \f$ a = \frac{2}{15}\sqrt{15 \sqrt{3} \pi} R \f$
 */
static const long double icosaEdge = 2.0L * originXOff;

////////////////////////////////////////////////////////////////////////////////
Vec2D sllxy (const GeoCoord& geoVect, SphIcosa& sphico, int nTri)
/*
   ISEA projection: from sphere to plane
   //参数为某一点的地理坐标，记录二十面体的全部信息的结构体，及所在的三角形的编号 ，返回值为对应的平面二维坐标

*/
{
   long double ph, fh, azh, azh1, dazh, h, dz, z, azh0, ag, cosAzh, sinAzh;
   Vec2D  Planevect;
    //找到了编号对应的三角形的信息
   const PreCompGeo& cent = sphico.triCen[nTri];

   long double cosLat = cosl(geoVect.lat);
   long double sinLat = sinl(geoVect.lat);

    // 球面三角面中心与上（下）顶点连接的大弧的方位角Az
   dazh = sphico.dazh[nTri];

   long double tmp = cent.sinLat * sinLat + cent.cosLat * cosLat *
       cosl(geoVect.lon - cent.pt.lon);
//   判断一下是否大雨1？为什么
   if (tmp > M_ONE) tmp = M_ONE;
   if (tmp < -M_ONE) tmp = -M_ONE;
    /*
  * 待计算点geoVect到三角面中心cent的球面距离z（单位：弧度），
  * 计算公式见文献[1]的(13)式。//参考施耐德原文
  */
   z = acosl(tmp);
    // z超过限差g(37.37736814),退出.
   if (z > DH + 0.00000005L)
   {
      dgcout << "nTri: " << nTri << "  z: " << z
             << "  DH+: " << DH + 0.00000005L << "  diff: "
             << ((DH + 0.00000005L) - z) << endl;
      dgcout << "1: The point: ";
      printGeoCoord(geoVect);
      dgcout << " is located on another polygon." << endl;
      report("Unable to continue.", DgBase::Fatal);
   }
    // 计算geoVect的方位角Az"与三角面中心cent方位角daz之差，这实际是在计算geoVect与三角面中心cent所形成的弧与中线之间的夹角
    //azh = azhr - dazh;//dazh为中线的方位角，azhr为AT的方位角
   azh = atan2l(cosLat * sinl(geoVect.lon - cent.pt.lon),
         cent.cosLat * sinLat - cent.sinLat * cosLat * cosl(geoVect.lon -
         cent.pt.lon)) - dazh;
    // 调整geoVect的等效点到基本单元内
   if (azh < 0.0) azh = azh + 2.0 * M_PI;
   azh0 = azh;//azh0记录的是原始的azh角
   if ((azh >= 120.0 * M_PI_180) && (azh <= 240.0 * M_PI_180)) azh -= 120.0 * M_PI_180;
   if (azh > 240.0 * M_PI_180) azh -= 240.0 * M_PI_180;

   cosAzh = cosl(azh);
   sinAzh = sinl(azh);
    // 计算球面上方位角为Az的大弧A'D'的长度q，计算公式见文献[1]的(9)式。
   dz = atan2l(tanDH, cosAzh + cot30 * sinAzh);
    // 检查是否超限
   if (z > dz + 0.00000005) {
      dgcout << "1: The point: ";
      printGeoCoord(geoVect);
      dgcout << " is located on another polygon." << endl;
      report("Unable to continue.", DgBase::Fatal);
   }
    // 公式(6)中H
   h = acosl(sinAzh * sinGH * cosDH - cosAzh * cosGH);
    // 公式(7)中AG
   ag = azh + GH + h - 180.0L * M_PI_180;
    // 计算球面方位角为Az的大弧投影到平面后对应的方位角Az'，计算公式见文献[1]的(8)式。
   azh1 = atan2l(2.0L * ag, R1S * tanDH * tanDH - 2.0L * ag * cot30);
    // 公式(11)中f即比例系数，将公式(10)代入(11)式即得
   fh = tanDH / (2.0L * (cosl(azh1) + cot30 * sinl(azh1)) * sinl(dz / 2.0L));
    // 公式(12)中的ρ
   ph = 2.0L * R1 * fh * sinl(z / 2.0L);
    // 将坐标还原
   if ((azh0 >= 120.0L * M_PI_180) && (azh0 < 240.0L * M_PI_180)) azh1 += 120.0L * M_PI_180;
   if (azh0 >= 240.0L * M_PI_180) azh1 += 240L * M_PI_180;
    // 坐标归"1"化, 坐标原点为三角形角点   为何要归一化？
   Planevect.x = (ph * sinl(azh1) + originXOff) / icosaEdge;
   Planevect.y = (ph * cosl(azh1) + originYOff) / icosaEdge;

   return (Planevect);

} /* Vec2D sllxy */

////////////////////////////////////////////////////////////////////////////////
IcosaGridPt snyderFwd (const GeoCoord& ll, DgSphIcosa& sphicosa)
/*
   project the point ll (lat, lon in radius) to the
   plane and return relevant info in IcosaGridPt
*/
{
   IcosaGridPt gridpt;

   gridpt.triangle = sphicosa.whichIcosaTri(ll);

   if (gridpt.triangle < 0)
   {
      dgcout << "ERROR: point in no triangle:";
      printGeoCoord(ll);
      dgcout << endl;

      gridpt.pt.x = M_ZERO;
      gridpt.pt.y = M_ZERO;

      return gridpt;
   }

   gridpt.pt = sllxy(ll, sphicosa.sphIcosa(), gridpt.triangle);

   return gridpt;

} /* IcosaGridPt snyderFwd */

////////////////////////////////////////////////////////////////////////////////
GeoCoord snyderInv (const IcosaGridPt& icosaPt, SphIcosa& sphicosa)
/*
    project the point icosaPt.pt (x, y in ISEA in the coordinate
    system of the specified triangle) to lon, lat in radians

    this replaces the old snyderInv/sxyll combo.
*/
{
  long double ddazh,ph,fh,azh,azh1,dazh,h,fazh,flazh,dz,z;
  long double sinlat,sinlon;
  long double azh0;
  GeoCoord Geovect;
  const PreCompGeo& cent = sphicosa.triCen[icosaPt.triangle];

  Vec2D pt;
  pt.x = icosaPt.pt.x * icosaEdge - originXOff;
  pt.y = icosaPt.pt.y * icosaEdge - originYOff;

  ddazh = sphicosa.dazh[icosaPt.triangle];

  if ((fabsl(pt.x) < PRECISION) && (fabsl(pt.y) < PRECISION))
  {
    Geovect.lat=cent.pt.lat; Geovect.lon=cent.pt.lon;
  }
  else
  {
    ph=sqrtl(pt.x*pt.x+pt.y*pt.y);
    azh1=atan2l(pt.x,pt.y);

    if (azh1<0.0L) azh1=azh1+2*M_PI;
    azh0=azh1;
    if ((azh1>120.0L*M_PI_180) && (azh1<=240.0L*M_PI_180)) azh1=azh1-120.0L*M_PI_180;
    if (azh1>240.0L*M_PI_180) azh1=azh1-240.0L*M_PI_180;

    azh=azh1;

    if (fabsl(azh1) > PRECISION)
    {
       long double agh=R1S*tanDH*tanDH/(2.0L*(1.0L/tanl(azh1)+cot30));

       //cout << "agh: " << agh << endl;

       dazh=1.0;
       while (fabsl(dazh) > PRECISION)
        {
         h=acosl(sinl(azh)*sinGH*cosDH-cosl(azh)*cosGH);
         fazh=agh-azh-GH-h+M_PI;
         flazh=((cosl(azh)*sinGH*cosDH+sinl(azh)*cosGH)/sinl(h))-1.0;
         dazh=-fazh/flazh;
         azh=azh+dazh;
         //cout << "loop: h: " << h << "  fazh: " << fazh*(180.0/M_PI) <<
         //  "  flazh: " << flazh*(180.0/M_PI) <<
         //  "  dazh: " << dazh*(180.0/M_PI) <<
         //  "  azh: " << azh*(180.0/M_PI) << endl;
        }
    }
    else azh = azh1 = 0.0;

    dz=atan2l(tanDH,cosl(azh)+cot30*sinl(azh));
    fh=tanDH/(2.0L*(cosl(azh1)+cot30*sinl(azh1))*sinl(dz/2.0L));
    z=2.0*asinl(ph/(2.0L*R1*fh));
    if ((azh0>=120*M_PI_180) && (azh0<240.0L*M_PI_180)) azh=azh+120*M_PI_180;
    if (azh0>=240.0L*M_PI_180) azh=azh+240.0L*M_PI_180;

    // now reposition to the actual triangle

    azh += ddazh;

    while (azh <= -M_PI) azh += M_2PI;
    while (azh > M_PI) azh -= M_2PI;

    sinlat=cent.sinLat * cosl(z) + cent.cosLat * sinl(z) * cosl(azh);
    if (sinlat > M_ONE) sinlat = M_ONE;
    if (sinlat < -M_ONE) sinlat = -M_ONE;
    Geovect.lat = asinl(sinlat);

    if (fabsl(fabsl(Geovect.lat) - M_PI_2) < M_EPSILON)
    {
       Geovect.lat = (Geovect.lat > M_ZERO) ? M_PI_2 : -M_PI_2;
       Geovect.lon = M_ZERO;
    }
    else
    {
      sinlon = sinl(azh)*sinl(z)/cosl(Geovect.lat);
      long double coslon = (cosl(z) - cent.sinLat * sinl(Geovect.lat)) /
              cent.cosLat/cosl(Geovect.lat);
      if (sinlon > M_ONE) sinlon = M_ONE;
      if (sinlon < -M_ONE) sinlon = -M_ONE;
      if (coslon > M_ONE) coslon = M_ONE;
      if (coslon < -M_ONE) coslon =-M_ONE;
      Geovect.lon = cent.pt.lon+asinl(sinlon);
      Geovect.lon = cent.pt.lon+atan2l(sinlon, coslon);
      if (Geovect.lon <= -M_PI) Geovect.lon += M_2PI;
      if (Geovect.lon >= M_PI) Geovect.lon -= M_2PI;
    }
  }
  return Geovect;

} /* GeoCoord snyderInv */

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
