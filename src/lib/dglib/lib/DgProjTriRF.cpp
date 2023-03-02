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
// DgProjTriRF.cpp: DgProjTriRF class implementation
//
////////////////////////////////////////////////////////////////////////////////

#include <dglib/DgProjTriRF.h>

////////////////////////////////////////////////////////////////////////////////
/////////////////////////////每次创建一层格网的时候都需要创建一个新的球面参数///////////////////////////////////////////////////
DgSphIcosa::DgSphIcosa (const DgGeoCoord& vert0, long double azimuthDegs)
{
   sphIcosa_.pt.lon = vert0.lon();
   sphIcosa_.pt.lat = vert0.lat();
   sphIcosa_.azimuth = azimuthDegs * M_PI/180;

   ico12verts();

} // DgSphIcosa::DgSphIcosa

////////////////////////////////////////////////////////////////////////////////
ostream& operator<< (ostream& str, const DgSphIcosa& dgsi)
{
   const SphIcosa& si = dgsi.sphIcosa_;

   DgGeoCoord tmp(si.pt);
   str << "vert0: " << tmp << endl;
   str << "az0: " << si.azimuth * 180.0 / M_PI << endl;

   str << "vertices:\n";
   str << "{\n";
   for (int i = 0; i < 12; i++)
   {
      DgGeoCoord tmp(si.icoverts[i]);
      str << "  " << i << ": " << tmp << endl;
   }
   str << "}\n";

   str << "triangles:\n";
   str << "{\n";
   for (int i = 0; i < 20; i++)
   {
      str << "  " << i;
      for (int j = 0; j < 3; j++)
      {
         DgGeoCoord tmp(si.icotri[i][j]);
         str << " " << tmp;
      }
      str << endl;
   }
   str << "}\n";

   return str;

} // ostream& operator<<

////////////////////////////////////////////////////////////////////////////////
GeoCoord coordtrans(const GeoCoord& newNPold, const GeoCoord& ptold,
                    long double lon0)
/*   
返回以newNPold为极点的经纬度坐标return使用方位投影公式 第3章- 第3讲方位投影及其应用 最后怎么返回真正的经纬度你
the new coordinates of any point in orginal coordinate system.
   Define a point (newNPold) in orginal coordinate system as the North
   Pole in new coordinate system, and the great circle connect the original
   and new North Pole as the lon0 longitude in new coordinate system, given
   any point in orginal coordinate system, this function return the new
   coordinates.*/

 {
  long double cosptnewlat, cosptnewlon;
  GeoCoord ptnew;

  cosptnewlat = sinl(newNPold.lat)*sinl(ptold.lat) +
                cosl(newNPold.lat)*cosl(ptold.lat)*cosl(newNPold.lon-ptold.lon);
  if (cosptnewlat>1.0L) cosptnewlat=1.0L;
  if (cosptnewlat<-1.0L) cosptnewlat=-1.0L;
  ptnew.lat = acosl(cosptnewlat);
  if (fabsl(ptnew.lat-0.) < PRECISION*100000)
      ptnew.lon=0.;
  else if (fabsl(ptnew.lat-M_PI) < PRECISION*100000)
      ptnew.lon=0.;
  else
   {
    cosptnewlon = (sinl(ptold.lat)*cosl(newNPold.lat) - cosl(ptold.lat)*
                  sinl(newNPold.lat)*cosl(newNPold.lon-ptold.lon))/sinl(ptnew.lat);
    if (cosptnewlon>1.0L) cosptnewlon=1.0L;
    if (cosptnewlon<-1.0L) cosptnewlon=-1.0L;
    ptnew.lon = acosl(cosptnewlon);
    if ((ptold.lon-newNPold.lon)>=0 && (ptold.lon-newNPold.lon) < M_PI)
      ptnew.lon=-ptnew.lon+lon0;
    else ptnew.lon=ptnew.lon+lon0;
    if (ptnew.lon>M_PI) ptnew.lon -= 2*M_PI;
    if (ptnew.lon<-M_PI) ptnew.lon += 2*M_PI;
   }
  ptnew.lat = M_PI/2-ptnew.lat;
  return ptnew;
 }

////////////////////////////////////////////////////////////////////////////////
int
DgSphIcosa::whichIcosaTri (const GeoCoord& pt)
/*
   Return the index of the icosahedron triangle in which pt occurs.
*/
{
   // start by assuming the first face is the minimum
   int minFace = 0;
   // 计算与三角形中心点大圆弧长 来比较离得最近的面是哪个
   long double minDist = DgGeoCoord::gcDist(sphIcosa().triCen[0].pt, pt);

   // compare against the rest of the faces
   for (int i = 1; i < 20; i++)
   {
      const DgGeoCoord& triCen = sphIcosa().triCen[i].pt;
      long double dist = DgGeoCoord::gcDist(triCen, pt);
      if (dist < minDist) {
         minDist = dist;
         minFace = i;
      }
   }

   return minFace;

} /* int DgSphIcosa::whichIcosaTri */

////////////////////////////////////////////////////////////////////////////////
void
DgSphIcosa::ico12verts (void)
/*
   //初始化正二十面体的参数 Fill in the icosahedron particulars given one point and one edge's azimuth.
*/
{
   GeoCoord newnpold, vertsnew[12];
   int i, j;
   Vec3D ptri[3];
  /* 构成三角形的顶点索引（符合右手法则，首个一定是上或下顶点） */

   static int verts[20][3] = {

	   {  0,  1,  2 },   /*  0 */
	   {  0,  2,  3 },   /*  1 */
	   {  0,  3,  4 },   /*  2 */
           {  0,  4,  5 },   /*  3 */
	   {  0,  5,  1 },   /*  4 */
	   {  6,  2,  1 },   /*  5 */
           {  7,  3,  2 },   /*  6 */
           {  8,  4,  3 },   /*  7 */
           {  9,  5,  4 },   /*  8 */
           { 10,  1,  5 },   /*  9 */
           {  2,  6,  7 },   /* 10 */
           {  3,  7,  8 },   /* 11 */
           {  4,  8,  9 },   /* 12 */
           {  5,  9, 10 },   /* 13 */
           {  1, 10,  6 },   /* 14 */
           { 11,  7,  6 },   /* 15 */
           { 11,  8,  7 },   /* 16 */
           { 11,  9,  8 },   /* 17 */
           { 11, 10,  9 },   /* 18 */
           { 11,  6, 10 }    /* 19 */
   };
  /* 标准二十面体顶点的坐标初始化 */ 
  //使用方位投影公式 第3章- 第3讲方位投影及其应用 意思是计算新极点出的经纬度
  // 新坐标系北极点的纬度
   newnpold.lat = sphIcosa().pt.lat;
     // 新坐标系北极点的经度是任意角度, 指定为0.
   newnpold.lon = 0.0; //极点处的角度是任意的
   for (i = 1; i <= 5; i++)
   {//先计算以真正的南北两级为顶点的坐标经纬度
          //纬度为26.565051177的5个点
     vertsnew[i].lat = 26.565051177 * M_PI / 180.0;
         //考虑到方位角的影响
     vertsnew[i].lon = -sphIcosa().azimuth + 72 * (i - 1) * M_PI / 180.0;
//     转换为【-pi pi 】区间
     if (vertsnew[i].lon > M_PI-PRECISION) vertsnew[i].lon -= 2 * M_PI;
     if (vertsnew[i].lon < -(M_PI+PRECISION)) vertsnew[i].lon += 2 * M_PI;
     vertsnew[i+5].lat = -26.565051177 * M_PI / 180;
     vertsnew[i+5].lon =
              -sphIcosa().azimuth + (36.0 + 72.0 * (i - 1)) * M_PI / 180.0;
     if (vertsnew[i + 5].lon > M_PI - PRECISION)
     {
        vertsnew[i + 5].lon -= 2 * M_PI;
     }
     if (vertsnew[i+5].lon < -(M_PI+PRECISION)) vertsnew[i+5].lon += 2 * M_PI;
   }
   vertsnew[11].lat = -90.0 * M_PI / 180.0;
   vertsnew[11].lon = 0.0 * M_PI / 180.0;

     //0#点是用户直接给定的点,不参与坐标变换.（此处涉及到关于正二十面体的定义的位置在问题中应该很重要 杜灵瑀 2016.08.04） 
       sphIcosa().icoverts[0].lat = sphIcosa().pt.lat;
       sphIcosa().icoverts[0].lon = sphIcosa().pt.lon;
/***** hardwire for bug test ******/


vertsnew[0].lat = 90.0 * M_PI / 180.0;
vertsnew[0].lon = 0.0 * M_PI / 180.0;
//for (i = 0; i < 12; i++)
//{
////   sphIcosa().icoverts[i].lat =  vertsnew[i];
////   sphIcosa().icoverts[i].lon = vertsnew[i].lon;
//    dgcout << "(" << i<<" "<< vertsnew[i].lon * M_180_PI << ", " <<  vertsnew[i].lat * M_180_PI << ")"<<"\n";
//}

  /* 根据用户的定义对二十面体的顶点进行坐标变换 */

/***************这个结果输出的是基于新极点的经纬度，并不是真正的经纬度，暂且不知道真正经纬度在哪转换的(输出的事真正的经纬度 但是按照第3章- 第3讲方位投影及其应用公式应该是
 关于新极点的新经纬度呀 而且coordtrans的解释也是这样说的 就是新极点下的经纬度 但是输出的值是在老坐标系下度经纬度值* )*******************/
   for (i = 0; i < 12; i++)
   {
      sphIcosa().icoverts[i] =
                   coordtrans(newnpold, vertsnew[i], sphIcosa().pt.lon);
//       dgcout << "(nnn" << i<<" "<<sphIcosa().icoverts[i].lon * M_180_PI << ", " << sphIcosa().icoverts[i].lat * M_180_PI << ")"<<"\n";

   }
   /* calculate the triangle-specific values */
  /* 计算以描述三角形的各项值 */

   for (i = 0; i < 20; i++)
   {
      /* set the vertices */

      for (j = 0; j < 3; j++)
      {
         sphIcosa().icotri[i][j] = sphIcosa().icoverts[verts[i][j]];
      }
    /* 预先计算与判断点在中心相关的数值 */

      /* pre-calculate the center point values */

      sphIcosa().triCen[i].pt = sphTricenpoint(sphIcosa().icotri[i]);
      sphIcosa().triCen[i].sinLat = sinl(sphIcosa().triCen[i].pt.lat);
      sphIcosa().triCen[i].sinLon = sinl(sphIcosa().triCen[i].pt.lon);
      sphIcosa().triCen[i].cosLat = cosl(sphIcosa().triCen[i].pt.lat);
      sphIcosa().triCen[i].cosLon = cosl(sphIcosa().triCen[i].pt.lon);

      /* pre-calculate the dazh's */
    /*
     * 预先计算球面三角形中心到上（下）顶点连线与经线（地理北）的夹角，详见文献[1]公式（14）。
     * 多面体与地球的位置关系可以任意定义，而Snyder投影的计算必须化归到基本单元中，因此必须
     * 记录这个方位角。
     */
      const GeoCoord* t = sphIcosa().icotri[i];
      const PreCompGeo& c = sphIcosa().triCen[i];
      sphIcosa().dazh[i] = atan2l(cosl(t[0].lat) * sinl(t[0].lon - c.pt.lon),
                  c.cosLat * sinl(t[0].lat) - c.sinLat * cosl(t[0].lat) *
                  cosl(t[0].lon - c.pt.lon));
      // 避免计算误差（贲进添加）
      if (fabs(sphIcosa().dazh[i]) < PRECISION)
         sphIcosa().dazh[i] = 0.0;
/*
      sphIcosa().dazh[i] = atan2l(cosl(sphIcosa().icotri[i][0].lat) *
                  sinl(sphIcosa().icotri[i][0].lon - sphIcosa().triCen[i].pt.lon),
                  sphIcosa().triCen[i].cosLat * sinl(sphIcosa().icotri[i][0].lat) -
                  sphIcosa().triCen[i].sinLat * cosl(sphIcosa().icotri[i][0].lat) *
                  cosl(sphIcosa().icotri[i][0].lon - sphIcosa().triCen[i].pt.lon));
*/

      /* pre-calculate the pt-in-triangle constants */
    /* 预先计算判断点在三角形中所需参数 */

      for (j = 0; j < 3; j++) ptri[j] = llxyz(sphIcosa().icotri[i][j]);

      sphIcosa().ptin[i].p0x0 = ptri[1].y * ptri[2].z - ptri[2].y * ptri[1].z;
      sphIcosa().ptin[i].p0y0 = ptri[1].x * ptri[2].z - ptri[2].x * ptri[1].z;
      sphIcosa().ptin[i].p0z0 = ptri[1].x * ptri[2].y - ptri[2].x * ptri[1].y;

      sphIcosa().ptin[i].t00 = ptri[0].x * sphIcosa().ptin[i].p0x0 -
                              ptri[0].y * sphIcosa().ptin[i].p0y0 +
                              ptri[0].z * sphIcosa().ptin[i].p0z0;

      sphIcosa().ptin[i].p0x1 = ptri[0].y * ptri[2].z - ptri[2].y * ptri[0].z;
      sphIcosa().ptin[i].p0y1 = ptri[0].x * ptri[2].z - ptri[2].x * ptri[0].z;
      sphIcosa().ptin[i].p0z1 = ptri[0].x * ptri[2].y - ptri[2].x * ptri[0].y;

      sphIcosa().ptin[i].t01 = ptri[1].x * sphIcosa().ptin[i].p0x1 -
                              ptri[1].y * sphIcosa().ptin[i].p0y1 +
                              ptri[1].z * sphIcosa().ptin[i].p0z1;

      sphIcosa().ptin[i].p0x2 = ptri[0].y * ptri[1].z - ptri[1].y * ptri[0].z;
      sphIcosa().ptin[i].p0y2 = ptri[0].x * ptri[1].z - ptri[1].x * ptri[0].z;
      sphIcosa().ptin[i].p0z2 = ptri[0].x * ptri[1].y - ptri[1].x * ptri[0].y;

      sphIcosa().ptin[i].t02 = ptri[2].x * sphIcosa().ptin[i].p0x2 -
                              ptri[2].y * sphIcosa().ptin[i].p0y2 +
                              ptri[2].z * sphIcosa().ptin[i].p0z2;

   }

/*
   for (int i = 0; i < 20; i++)
   {
      printf("%d\n", i);
      for (int j = 0; j < 3; j++) printGeoCoord(sphIcosa().icotri[i][j]);
      printGeoCoord(sphIcosa().triCen[i].pt);
      printf("\n");
   }
*/

} /* void ico12verts */

////////////////////////////////////////////////////////////////////////////////
const char*
DgProjTriRF::str2add (DgProjTriCoord* add, const char* str, char delimiter)
                                                                          const
{
   if (!add) add = new DgProjTriCoord();

   char delimStr[2];
   delimStr[0] = delimiter;
   delimStr[1] = '\0';

   char* tmpStr = new char[strlen(str) + 1];
   strcpy(tmpStr, str);

   // get the triNum

   char* tok = strtok(tmpStr, delimStr);
   int t;
   if (sscanf(tok, "%d", &t) != 1)
   {
      ::report("DgProjTriRF::fromString() invalid value in string " +
               string(tok), DgBase::Fatal);
   }

   const char* tmp = &(str[strlen(tok) + 1]);
   DgDVec2D vec;
   tmp = vec.fromString(tmp, delimiter);

   *add = DgProjTriCoord(t, vec);

   return tmp;

} // const char* DgProjTriRF::str2add
