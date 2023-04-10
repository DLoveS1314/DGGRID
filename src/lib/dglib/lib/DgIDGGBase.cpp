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
// DgIDGGBase.cpp: DgIDGGBase class implementation
//
// Version 7.0 - Kevin Sahr, 11/15/14
// Version 6.1 - Kevin Sahr, 5/23/13
//
////////////////////////////////////////////////////////////////////////////////

#include <dglib/DgBoundedIDGG.h>
#include <dglib/DgDmdD4Grid2DS.h>
#include <dglib/DgHexGrid2DS.h>
#include <dglib/DgIDGGBase.h>
#include <dglib/DgIDGGSBase.h>
#include <dglib/DgProjFuller.h>
#include <dglib/DgProjISEA.h>
#include <dglib/DgRadixString.h>
#include <dglib/DgSeriesConverter.h>
#include <dglib/DgTriGrid2DS.h>

#include <cfloat>
#include <climits>
#include <cmath>
#include <set>

////////////////////////////////////////////////////////////////////////////////
const DgGeoSphRF& DgIDGGBase::geoRF      (void) const { return dggs()->geoRF(); }
const DgGeoCoord& DgIDGGBase::vert0      (void) const { return dggs()->vert0(); }
long double       DgIDGGBase::azDegs     (void) const { return dggs()->azDegs(); }
const string&     DgIDGGBase::projType   (void) const { return dggs()->projType(); }
DgGridTopology    DgIDGGBase::gridTopo   (void) const { return dggs()->gridTopo(); }
DgGridMetric      DgIDGGBase::gridMetric (void) const { return dggs()->gridMetric(); }

////////////////////////////////////////////////////////////////////////////////
const DgQuadEdgeCells DgIDGGBase::edgeTable_[12] = {

   DgQuadEdgeCells(0,  true,  0,  0, 0), // quad 0 should never occur
   DgQuadEdgeCells(1,  true,  0,  2, 6),
   DgQuadEdgeCells(2,  true,  0,  3, 7),
   DgQuadEdgeCells(3,  true,  0,  4, 8),
   DgQuadEdgeCells(4,  true,  0,  5, 9),
   DgQuadEdgeCells(5,  true,  0,  1, 10),
   DgQuadEdgeCells(6,  false, 11, 2, 7),
   DgQuadEdgeCells(7,  false, 11, 3, 8),
   DgQuadEdgeCells(8,  false, 11, 4, 9),
   DgQuadEdgeCells(9,  false, 11, 5, 10),
   DgQuadEdgeCells(10, false, 11, 1, 6),
   DgQuadEdgeCells(11, false, 11, 0, 0)  // quad 11 should never occur

};

////////////////////////////////////////////////////////////////////////////////
const char*
DgIDGGBase::str2add (DgQ2DICoord* add, const char* str, char delimiter) const
{
   if (!add) add = new DgQ2DICoord();

   char delimStr[2];
   delimStr[0] = delimiter;
   delimStr[1] = '\0';

   char* tmpStr = new char[strlen(str) + 1];
   strcpy(tmpStr, str);

   // get the quadNum

   char* tok = strtok(tmpStr, delimStr);
   int q;
   if (sscanf(tok, "%d", &q) != 1)
   {
      ::report("DgQ2DIRF::fromString() invalid value in string " +
               string(tok), DgBase::Fatal);
   }

   const char* tmp = &(str[strlen(tok) + 1]);
   DgIVec2D vec;
   tmp = vec.fromString(tmp, delimiter);

   *add = DgQ2DICoord(q, vec);

   return tmp;

} // const char* DgIDGGBase::str2add

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
DgIDGGBase::DgIDGGBase (const DgIDGGSBase* dggs, const DgGeoSphRF& geoRF,
             unsigned int aperture, int res, const string& name,
             DgGridTopology gridTopo, DgGridMetric gridMetric,
             unsigned int precision)
   : DgDiscRF<DgQ2DICoord, DgGeoCoord, long double>
          (geoRF.network(), geoRF, name, gridTopo, gridMetric),
     dggs_ (dggs), sphIcosa_(0), aperture_(aperture), res_(res),
     precision_(precision), grid2D_(0), grid2DS_(0), ccFrame_(0),
     projTriRF_(0), vertexRF_(0), q2ddRF_(0), bndRF_(0), intRF_(0), planeRF_(0)
{//参数在相应的子类中进行初始化 如DgDmdIDGG::initialize
   //initialize();

} // DgIDGGBase::DgIDGGBase

////////////////////////////////////////////////////////////////////////////////
/*
DgIDGGBase::DgIDGGBase (const DgIDGGBase& rfIn)
   : DgDiscRF<DgQ2DICoord, DgGeoCoord, long double> (rfIn),
        dggs_ (NULL), sphIcosa_(0), aperture_(rfIn.aperture()),
        res_(rfIn.res()), precision_(rfIn.precision()),
        grid2D_(0), grid2DS_(0), ccFrame_(0), projTriRF_(0),
        vertexRF_(0), q2ddRF_(0), bndRF_(0), intRF_(0), planeRF_(0)
{
   //initialize();

} // DgIDGGBase::DgIDGGBase
*/

////////////////////////////////////////////////////////////////////////////////
DgIDGGBase::~DgIDGGBase()
{
 delete sphIcosa_;
 delete bndRF_;
}

////////////////////////////////////////////////////////////////////////////////
void
DgIDGGBase::createConverters (void)
{

//   grid2DS()在子类中进行赋值 grid2DS().grids()[res()] 这个得到是一个 vector<const DgDiscRF<A, B, DB>*>结构 强制转化为子类 这个是一个空数组
//grid2DS 有很多层级 只取一级即可表示当前的grid2D_
   grid2D_ = dynamic_cast<const DgDiscRF2D*>(grid2DS().grids()[res()]);
   //cout << "== GRID2D: " << string(*grid2D_);
//    bndRF_ 赋值位置
   bndRF_ = new DgBoundedIDGG(*this);
   //cout << "== BNDRF:: " << string(*bndRF_) << endl;

   // create the intermediate RFs 是全局的参考系统

   projTriRF_ = DgProjTriRF::makeRF(network(), name() + string("projTri"),
                sphIcosa_);//投影空间
   vertexRF_ = DgVertex2DDRF::makeRF(network(), name() + string("vertex"));//顶点空间
   q2ddRF_ = DgQ2DDRF::makeRF(network(), name() + string("q2dd"));
   intRF_ = DgInterleaveRF::makeRF(network(), name() + string("int"));
   planeRF_ = DgPlaneTriRF::makeRF(network(), name() + string("plane"));

   // create the converters; for convenience use where they are in overall
   // sequence for name

   DgIcosaProj* icosaProj = NULL;
   if (projType() == "ISEA")
      icosaProj = new DgProjISEA(geoRF(), projTriRF());
   else if (projType() == "FULLER")
      icosaProj = new DgProjFuller(geoRF(), projTriRF());
   else
      report("DgIDGGBase::initialize(): invalid projection type " + projType(),
             DgBase::Fatal);
    //正算，从经纬度 （ DgGeoSphRF）到离散格网空间 q ,i j 标识
   const DgConverterBase* c1to2 = &(icosaProj->forward());
   const DgConverterBase* c2to3 = new DgProjTriToVertex2DD(projTriRF(), vertexRF());
   const DgConverterBase* c3to4 = new DgVertex2DDToQ2DDConverter(vertexRF(), q2ddRF());
   const DgConverterBase* c4to5 = new DgQ2DDtoIConverter(q2ddRF(), *this);
    //反算，从离散格网空间 q ,i j 到标识经纬度 （ DgGeoSphRF）
   const DgConverterBase* c5to4 = new DgQ2DItoDConverter(*this, q2ddRF()); //离散空间到归算空间
   const DgConverterBase* c4to3 = new DgQ2DDtoVertex2DDConverter(q2ddRF(), vertexRF());//归算空间到顶点空间
   const DgConverterBase* c3to2 = new DgVertex2DDtoProjTri(vertexRF(), projTriRF());//顶点空间到投影空间
   const DgConverterBase* c2to1 = &(icosaProj->inverse());//投影空间逆转换为经纬度

   // done with icosaProj; the fwd/inv converters are in the RFNetwork
   delete icosaProj;

   DgConverterBase* toInt = new DgQ2DItoInterleaveConverter(*this, intRF());
   DgConverterBase* toPlane = new DgPlaneTriProj(projTriRF(), planeRF());

   // create the series converters that will replace the default DgDiscRF
   // converters

   vector<const DgConverterBase*> sc;
   sc.push_back(c1to2);
   sc.push_back(c2to3);
   sc.push_back(c3to4);
   sc.push_back(c4to5);
   new DgSeriesConverter(sc, true);
   sc.resize(0);

   sc.push_back(c5to4);
   sc.push_back(c4to3);
   sc.push_back(c3to2);
   sc.push_back(c2to1);
   new DgSeriesConverter(sc, true);
   sc.resize(0);

   //// now fill in the other connections; begin with from vertexRF and then
   //// connect to it as needed

   // vertexRF -> geoRF
   sc.push_back(c3to2);
   sc.push_back(c2to1);
   new DgSeriesConverter(sc, true);
   sc.resize(0);

   // vertexRF -> Q2DI
   sc.push_back(c3to4);
   sc.push_back(c4to5);
   new DgSeriesConverter(sc, true);
   sc.resize(0);

   // vertexRF -> projTriRF is c3to2 above

   // vertexRF -> planeRF
   sc.push_back(c3to2);
   sc.push_back(toPlane);
   new DgSeriesConverter(sc, true);
   sc.resize(0);

   // vertexRF -> Q2DD is c3to4 above

   // vertexRF -> intRF
   sc.push_back(c3to4);
   sc.push_back(c4to5);
   sc.push_back(toInt);
   new DgSeriesConverter(sc, true);
   sc.resize(0);

   /// now do from projTriRF

   // projTriRF -> geoRF is c2to1 above

   // projTriRF -> Q2DI
   sc.push_back(c2to3);
   sc.push_back(network().getConverter(vertexRF(), *this));
   new DgSeriesConverter(sc, true);
   sc.resize(0);

   // projTriRF -> vertexRF is c2to3 above
   // projTriRF -> planeRF is toPlane above

   // projTriRF -> Q2DD
   sc.push_back(c2to3);
   sc.push_back(c3to4);
   new DgSeriesConverter(sc, true);
   sc.resize(0);

   // projTriRF -> intRF
   sc.push_back(c2to3);
   sc.push_back(network().getConverter(vertexRF(), intRF()));
   new DgSeriesConverter(sc, true);
   sc.resize(0);

   /// do from Q2DD

   // Q2DD -> geoRF
   sc.push_back(c4to3);
   sc.push_back(network().getConverter(vertexRF(), geoRF()));
   new DgSeriesConverter(sc, true);
   sc.resize(0);

   // Q2DD -> Q2DI is c4to5 above
   // Q2DD -> vertexRF is c4to3 above

   // Q2DD -> projTriRF
   sc.push_back(c4to3);
   sc.push_back(network().getConverter(vertexRF(), projTriRF()));
   new DgSeriesConverter(sc, true);
   sc.resize(0);

   // Q2DD -> planeRF
   sc.push_back(c4to3);
   sc.push_back(network().getConverter(vertexRF(), planeRF()));
   new DgSeriesConverter(sc, true);
   sc.resize(0);

   // Q2DD -> intRF
   sc.push_back(c4to3);
   sc.push_back(network().getConverter(vertexRF(), intRF()));
   new DgSeriesConverter(sc, true);
   sc.resize(0);

   /// do from Q2DI

   // Q2DI -> geoRF is series converter given above

   // Q2DI -> vertexRF
   sc.push_back(c5to4);
   sc.push_back(c4to3);
   new DgSeriesConverter(sc, true);
   sc.resize(0);

   // Q2DI -> projTriRF
   sc.push_back(c5to4);
   sc.push_back(network().getConverter(q2ddRF(), projTriRF()));
   new DgSeriesConverter(sc, true);
   sc.resize(0);

   // Q2DI -> planeRF
   sc.push_back(c5to4);
   sc.push_back(network().getConverter(q2ddRF(), planeRF()));
   new DgSeriesConverter(sc, true);
   sc.resize(0);

   // Q2DI -> Q2DD is c5to4 above
   // Q2DI -> intRF is toInt above

   /// finally from geoRF

   // geoRF -> Q2DI is series converter given above

   // geoRF -> vertexRF
   sc.push_back(c1to2);
   sc.push_back(c2to3);
   new DgSeriesConverter(sc, true);
   sc.resize(0);

   // geoRF -> projTriRF is c1to2 above

   // geoRF -> planeRF
   sc.push_back(network().getConverter(geoRF(), vertexRF()));
   sc.push_back(network().getConverter(vertexRF(), planeRF()));
   new DgSeriesConverter(sc, true);
   sc.resize(0);

   // geoRF -> Q2DD
   sc.push_back(c1to2);
   sc.push_back(c2to3);
   sc.push_back(c3to4);
   new DgSeriesConverter(sc, true);
   sc.resize(0);

   // geoRF -> intRF
   sc.push_back(network().getConverter(geoRF(), *this));
   sc.push_back(toInt);
   new DgSeriesConverter(sc, true);
   sc.resize(0);

} // DgIDGGBase::createConverters

//////////////////////////获取边界点//////////////////////////////////////////////////////
void
DgIDGGBase::setVertices (const DgLocation& loc, DgPolygon& vec,
                     int densify) const
{
   vec.clearAddress();
   backFrame().convert(vec);

   DgLocation tLoc(loc);
//cout << "*** " << loc << endl;
   convert(&tLoc);
//cout << "**** " << tLoc << endl;

   setAddVertices(*getAddress(tLoc), vec, densify);

} // void DgDiscRF::setVertices

////////////////////////////////////////////////////////////////////////////////
void
DgIDGGBase::setAddVertices (const DgQ2DICoord& add, DgPolygon& vec,
                        int densify) const
{
    //转到不带面编码的行列号模式，并且参考系转成grid2D对应的参考系
   DgLocation* tmpLoc = grid2D().makeLocation(add.coord());
//cout << "a: " << *tmpLoc << endl;
    DgPolygon dummy(ccFrame());
    vec = dummy;  // force empty RF to allow for network change
    grid2D().setVertices(*tmpLoc, vec);//得到过度坐标系下坐标
   delete tmpLoc;

//cout << "A: " << vec << endl;
   ccFrame().convert(vec);
//cout << "B: " << vec << endl;

   // densify
   vec.densify(densify);
//cout << "C: " << vec << endl;

   // kludge to jump nets and add the quad number

   DgPolygon tmpVec(q2ddRF());
   vector<DgAddressBase*>& v = tmpVec.addressVec();
   for (int i = 0; i < vec.size(); i++)
   {
      v.push_back(new DgAddress<DgQ2DDCoord>(DgQ2DDCoord(add.quadNum(),
                       *(ccFrame().getAddress(vec[i])))));
   }
   vec = tmpVec;

//cout << "D: " << vec << endl;

   vertexRF().convert(vec);

//cout << "E: " << vec << endl;

   if (!isCongruent() && add.coord() == DgIVec2D(0, 0))
   {
      // we need to explicitly go to vertexRF to look for non-keepers
      // to clip
      vector<DgAddressBase*>& v = vec.addressVec();
      vector<DgAddressBase*> newV;
      for (unsigned long i = 0; i < v.size(); i++) {
         DgAddress<DgVertex2DDCoord>* fullAdd =
                   dynamic_cast< DgAddress<DgVertex2DDCoord>* >(v[i]);
         DgVertex2DDCoord& add = fullAdd->address();
         if (add.keep()) {
            newV.push_back(v[i]);
         } else {
            delete v[i];
         }
      }

      v.resize(0);
      for (unsigned long i = 0; i < newV.size(); i++) {
         v.push_back(newV[i]);
      }
   }

   // now convert to the geoRF

//cout << "F: " << vec << endl;
   geoRF().convert(&vec);
//cout << "G: " << vec << endl;

   // Release the Kraken... I mean, the vector's pointers:
   dgg::util::release(v);

} // DgIDGGBase::setAddVertices

////////////////////////////可以吧这个删了加上自己的内容 不需要返回值 因为是引用传递////////////////////////////////////////////////////
void
DgIDGGBase::setAddNeighbors (const DgQ2DICoord& add,
                                   DgLocVector& vec) const
{
// //重写了 DgDiscRF setAddNeighbors
//DgLocVector在 DgDiscRF setNeighbors 已经转化为了DGIDGG空间 也就是Q i j
   if (this->gridTopo()== Diamond)
   {
       if (this->gridMetric() == D4)//现在不存在D8的菱形 所以只能D4菱形 也就是只能边临近和角临近选择一样
       {
           //八临近
           vector <DgQ2DICoord> coords =  this->neicell(add);
           vec.clearAddress();
           for (int i = 0; i < coords.size(); i++)
               vec.push_back(*makeLocation(coords[i]));
       }
   }
   else if (this->gridTopo()== Triangle)
   {
		//九个编码 ,二阶邻近。
		vector <DgQ2DICoord> coords =  this->trisecnei(add);
		vec.clearAddress();
		for (int i = 0; i < coords.size(); i++)
			vec.push_back(*makeLocation(coords[i]));
   }
   //源程序
   else
   {
       //空的vector d但是 参考系变成了grid2D()
       DgLocVector ngh2d(grid2D());
       grid2D().setAddNeighbors(add.coord(), ngh2d);
//cout << " >> DgIDGGBase::setAddNeighbors center: " << add << endl;
//cout << "  ngh2d: " << ngh2d << endl;
//cout << " isCongruent: " << (isCongruent() ? "yes" : "no");

       int q = add.quadNum();
       DgLocVector ngh2dNoDup(*this);
       vec.clearAddress();
       DgIVec2D c;
       for (int i = 0; i < ngh2d.size(); i++)
       {
           DgQ2DICoord c2di(q, *grid2D().getAddress(ngh2d[i]));
//cout << "*** i: " << i << " " << c2di;
           c2di = bndRF().q2dixToQ2di(c2di);
//cout << " -> " << c2di << endl;

           // check for duplicates
           bool keeper = true;
           if (!isCongruent() && add.coord() == DgIVec2D(0, 0))
           {
               for (int i = 0; i < ngh2dNoDup.size(); i++)
               {
                   const DgQ2DICoord& veci = *this->getAddress(ngh2dNoDup[i]);
//cout << "   " << i << " " << veci << " -> " << (c2di == veci) <<  endl;
                   if (c2di == veci)
                   {
                       keeper = false;
                       break;
                   }
               }
           }

           if (keeper)
           {
               DgLocation* tmpLoc = this->makeLocation(c2di);
               ngh2dNoDup.push_back(*tmpLoc);
               delete tmpLoc;
           }
       }

//cout << "ngh2dNoDup: " << ngh2dNoDup << endl;
       // now build the vector; the push_back will take care of converting
       for (int i = 0; i < ngh2dNoDup.size(); i++)
           vec.push_back(ngh2dNoDup[i]);

//cout << "final neigh vec for add: " << add << endl;
//cout << vec << endl;
//cout << "-------" << endl;

   }

}


////////////////////////////////////////////////////////////////////////////////
void
DgIDGGBase::setAddNeighborsBdry2 (const DgQ2DICoord& add,
                                   DgLocVector& vec) const
{
   DgLocVector ngh2d(grid2D());
   grid2D().setAddNeighborsBdry2(add.coord(), ngh2d);
//cout << " >> DgIDGGBase::setAddNeighborsBdry2:  ngh2d: " << endl;
//cout << ngh2d << endl;

   int q = add.quadNum();
   DgLocVector ngh2dNoDup(*this);
   vec.clearAddress();
   DgIVec2D c;
   for (int i = 0; i < ngh2d.size(); i++)
   {
      DgQ2DICoord c2di(q, *grid2D().getAddress(ngh2d[i]));
//cout << "*** i: " << i << " " << c2di;
      c2di = bndRF().q2dixToQ2di(c2di);
//cout << " -> " << c2di << endl;

      // check for duplicates
      bool keeper = true;
      if (!isCongruent() && add.coord() == DgIVec2D(0, 0))
      {
         for (int i = 0; i < ngh2dNoDup.size(); i++)
         {
            const DgQ2DICoord& veci = *this->getAddress(ngh2dNoDup[i]);
//cout << "   " << i << " " << veci << " -> " << (c2di == veci) <<  endl;
            if (c2di == veci)
            {
               keeper = false;
               break;
            }
         }
      }

      if (keeper)
      {
         DgLocation* tmpLoc = this->makeLocation(c2di);
         ngh2dNoDup.push_back(*tmpLoc);
         delete tmpLoc;
      }
   }

//cout << "ngh2dNoDup: " << ngh2dNoDup << endl;
   // now build the vector; the push_back will take care of converting
   for (int i = 0; i < ngh2dNoDup.size(); i++)
      vec.push_back(ngh2dNoDup[i]);

//cout << "vec: " << vec << endl;

}
//特殊处理的八邻近 每个格网都值生成八个格网 七个的补充一个 九个的删去一个 因为七个和九个的极少 所以不影响结果，按照卷积顺序排列好
vector<DgQ2DICoord> DgIDGGBase::neicell(const DgQ2DICoord &add1) const 
{
    vector<DgQ2DICoord> neiadd;
    //maxI和MaxJ 对于菱形 是相等的所以下面出现了混用的情况
    long long int maxI =this->maxI();
    long long int maxJ =this->maxJ();
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

vector<DgQ2DICoord> DgIDGGBase::triedgenei(const DgQ2DICoord &add1) const
{
	
	//对于三角形的行列 行对应的就是菱形的行 列分两种 首先基础列j是菱形所在列 包含的两个三角形列为 2*j 2*j +1 分别是下三角形
	//和上三角形
	// 总结跨面三角形 一共有第0行的上三角形 最大行的下三角形  第0列和最大列  跨两个面的0行最大列 0列最大行
	vector<DgQ2DICoord> neiadd;
	long long int maxI =this->maxI();
    long long int maxJ =this->maxJ();
	//行列均不为0和最大值 表示在三角形内部 属于是菱形面片的内部
	if(add1.coord().i() != 0 && add1.coord().i() != maxI &&
		add1.coord().j() != 0 && add1.coord().j() != maxJ) 
	{
		//面内三角形邻近查找
		if(add1.coord().j() % 2 != 0) //列是奇数列的全部是上三角形
		{
			//面内上三角形的三个边邻近单元
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() + 0, add1.coord().j() - 1)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() + 0, add1.coord().j() + 1)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() - 1, add1.coord().j() - 1)));

		}
		else //下三角形
		{
			//面内下三角形的三个边邻近单元 
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() + 0, add1.coord().j() + 1)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() - 0, add1.coord().j() - 1)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() + 1, add1.coord().j() + 1)));
		
		}
	}
	//第0行 的下三角形 且不是第0列的下三角形 下三角形的第0列是跨面的特殊处理 
	// 下三角形列不会是最大列 但不代表最大列不是跨面的
	else if(add1.coord().i() == 0 && add1.coord().j() % 2 == 0 	&& add1.coord().j() != 0)
	
	{
		//面内下三角形的三个边邻近单元 
		neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() + 0, add1.coord().j() + 1)));
		neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() - 0, add1.coord().j() - 1)));
		neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() + 1, add1.coord().j() + 1)));

	}
	//第大行的上三角形 且不是最大列的上三角形 上三角形的最大列是跨面的特殊处理 上三角形列不会是第0列
	else if(add1.coord().i() == maxI && add1.coord().j() % 2 == 1 && add1.coord().j() != maxJ)
	{
		//面内上三角形的三个边邻近单元
		neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() + 0, add1.coord().j() - 1)));
		neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() + 0, add1.coord().j() + 1)));
		neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() - 1, add1.coord().j() - 1)));
	}

	//跨面的情况

	// 行为0  最小行 列不是最大列的奇数列（上三角形）  跨一个面 上三角形
	else if(add1.coord().i() == 0 && add1.coord().j() % 2 == 1 && add1.coord().j() != maxJ)
	{
		if(add1.quadNum() == 1)
		{
			
			//跨面三角形邻近查找
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() + 0, add1.coord().j() - 1)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() + 0, add1.coord().j() + 1)));
			// if((maxJ / 2 - add1.coord().j() / 2) = (maxI -add1.coord().i()) )
			// { cout<<"results same"<<endl;}
			// maxJ / 2 - add1.coord().j() / 2 不就是最大行减去当前行？和菱形的计算是一样的
			neiadd.push_back(DgQ2DICoord(add1.quadNum() + 4, DgIVec2D(maxJ / 2 - add1.coord().j() / 2, maxJ)));
		}
		else if(add1.quadNum() > 1 && add1.quadNum() <= 5)
		{
			//跨面三角形邻近查找
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() + 0, add1.coord().j() - 1)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() + 0, add1.coord().j() + 1)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() - 1, DgIVec2D(maxJ / 2 - add1.coord().j() / 2, maxJ)));
		}
		else if(add1.quadNum() >= 6 && add1.quadNum() <= 10)
		{
			//跨面三角形邻近查找
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() + 0, add1.coord().j() - 1)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() + 0, add1.coord().j() + 1)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() - 5, DgIVec2D(maxI, add1.coord().j() - 1)));
		}

	}
	
	// 行为0 最小行 列是最大列    跨两个面  上三角形
	else if(add1.coord().i() == 0 && add1.coord().j() % 2 == 1 && add1.coord().j() == maxJ)
	{
		if(add1.quadNum() == 1)
		{
			//跨面上三角形邻近查找
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() + 0, add1.coord().j() - 1)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() + 1, DgIVec2D( 0, maxJ)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() + 4, DgIVec2D(0, maxJ)));
		}
		else if(add1.quadNum() > 1 && add1.quadNum() < 5)
		{
			//跨面三角形邻近查找
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() + 0, add1.coord().j() - 1)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() + 1, DgIVec2D( 0, maxJ)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() - 1, DgIVec2D(0, maxJ)));
		}
		else if(add1.quadNum() == 5)
		{
			//跨面三角形邻近查找
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() + 0, add1.coord().j() - 1)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() - 4, DgIVec2D( 0, maxJ)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() - 1, DgIVec2D(0, maxJ)));
		}
		else if(add1.quadNum() >= 6 && add1.quadNum() < 10)//源程序代码给的有错误 10的情况不一样
		{
			//跨面三角形邻近查找
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() + 0, add1.coord().j() - 1)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() - 4, DgIVec2D( 0, 0)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() - 5, DgIVec2D( maxI, maxJ - 1)));

		}
		else if(add1.quadNum()== 10)
		{
			//跨面三角形邻近查找 
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() + 0, add1.coord().j() - 1)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() - 9, DgIVec2D( 0, 0)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() - 5, DgIVec2D( maxI, maxJ - 1)));

		}

	}
	// 最大行 不等于0的偶数列 跨一个面 下三角形 进行修改 让下三角形从顶开始逆时针输出
	else if(add1.coord().i() == maxI && add1.coord().j() % 2 == 0 && add1.coord().j() != 0)
	{
		neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() + 0, add1.coord().j() + 1)));
		neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() - 0, add1.coord().j() - 1)));
		if(add1.quadNum() >= 1 && add1.quadNum() <= 5)
		{
			neiadd.push_back(DgQ2DICoord(add1.quadNum() + 5, DgIVec2D(0, add1.coord().j() + 1)));
		}
		else if(add1.quadNum() >= 6 && add1.quadNum() < 10)
		{
			neiadd.push_back(DgQ2DICoord(add1.quadNum() + 1, DgIVec2D(maxJ / 2 - add1.coord().j() / 2, 0)));
		}
		else if(add1.quadNum() == 10)
		{
			neiadd.push_back(DgQ2DICoord(add1.quadNum() - 4, DgIVec2D(maxJ / 2 - add1.coord().j() / 2, 0)));
		}
	}
	// 最大行 等于0的偶数列 （最小列）  跨两个面 下三角形
	else if(add1.coord().i() == maxI && add1.coord().j() == 0)
	{
		neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() + 0, add1.coord().j() + 1)));

		if(add1.quadNum() == 1)
		{
			//跨面下三角形的三个边邻近单元 
			neiadd.push_back(DgQ2DICoord(add1.quadNum() + 9 , DgIVec2D(maxI, maxJ)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() + 5, DgIVec2D(0, add1.coord().j() + 1)));
		}
		else if(add1.quadNum() > 1 && add1.quadNum() <= 5)
		{
			//跨面下三角形的三个边邻近单元 

			neiadd.push_back(DgQ2DICoord(add1.quadNum() + 4, DgIVec2D(maxI, maxJ)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() + 5, DgIVec2D(0, add1.coord().j() + 1)));
		}
		else if(add1.quadNum() == 6)
		{
			//跨面下三角形的三个边邻近单元 
			neiadd.push_back(DgQ2DICoord(add1.quadNum() + 4, DgIVec2D(maxI, 0)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() + 1, DgIVec2D(maxI, 0)));
		}
		else if(add1.quadNum() > 6 && add1.quadNum() < 10)
		{
			//跨面下三角形的三个边邻近单元 
			neiadd.push_back(DgQ2DICoord(add1.quadNum() - 1, DgIVec2D(maxI, 0)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() + 1, DgIVec2D(maxI, 0)));
		}
		else if(add1.quadNum() == 10)
		{
			//跨面下三角形的三个边邻近单元 
			neiadd.push_back(DgQ2DICoord(add1.quadNum() - 1, DgIVec2D(maxI, 0)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum() - 4, DgIVec2D(maxI, 0)));
		}
	}
	
	// 最小列 行不是最大行 下三角形 跨一面
	else if( add1.coord().j() == 0 &&  add1.coord().i() != maxI)
	{
		neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() + 0, add1.coord().j() + 1)));

		if(add1.quadNum() == 1)
		{
			 
			neiadd.push_back(DgQ2DICoord(add1.quadNum() + 9, DgIVec2D(add1.coord().i(), maxJ)));
		}
		else if(add1.quadNum() > 1 && add1.quadNum() <= 5)
		{
			neiadd.push_back(DgQ2DICoord(add1.quadNum() + 4, DgIVec2D(add1.coord().i(), maxJ)));
		}
		else if(add1.quadNum() == 6)
		{
			neiadd.push_back(DgQ2DICoord(add1.quadNum() + 4, DgIVec2D(maxI, (maxI - add1.coord().i()) * 2)));
		}
		else if(add1.quadNum() > 6 && add1.quadNum() <= 10)
		{
			neiadd.push_back(DgQ2DICoord(add1.quadNum() - 1, DgIVec2D(maxI, (maxI - add1.coord().i()) * 2)));
		}
		neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() + 1, add1.coord().j() + 1)));
	}
	// 最大列 行不是最小行 跨一个面  上三角形
	else if( add1.coord().j() == maxJ &&  add1.coord().i() != 0)
	{
		neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() + 0, add1.coord().j() - 1)));
		if(add1.quadNum() >= 1 && add1.quadNum() < 5)
		{
			//跨面上三角形的三个边邻近单元
			
			neiadd.push_back(DgQ2DICoord(add1.quadNum() + 1, DgIVec2D(0, (maxI - add1.coord().i()) * 2 + 1)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() - 1, add1.coord().j() - 1)));
		}
		else if(add1.quadNum() == 5)
		{
			//跨面上三角形的三个边邻近单元
			
			neiadd.push_back(DgQ2DICoord(add1.quadNum() - 4, DgIVec2D(0, (maxI - add1.coord().i()) * 2 + 1)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() - 1, add1.coord().j() - 1)));
		}
		else if(add1.quadNum() >= 6 && add1.quadNum() < 10)
		{
			//跨面上三角形的三个边邻近单元
			neiadd.push_back(DgQ2DICoord(add1.quadNum() - 4, DgIVec2D(add1.coord().i(), 0)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() - 1, add1.coord().j() - 1)));
		}
		else if(add1.quadNum() == 10)
		{
			//跨面上三角形的三个边邻近单元
			neiadd.push_back(DgQ2DICoord(add1.quadNum() - 9, DgIVec2D(add1.coord().i(), 0)));
			neiadd.push_back(DgQ2DICoord(add1.quadNum(), DgIVec2D(add1.coord().i() - 1, add1.coord().j() - 1)));
		}
	}
	return neiadd;

}
 
vector<DgQ2DICoord> DgIDGGBase::trisecnei(const DgQ2DICoord &add1) const
{
	vector<DgQ2DICoord> firstnei = this->triedgenei(add1);
	std::vector<DgQ2DICoord> secnei_vec;
	secnei_vec.insert(secnei_vec.end(), firstnei.begin(), firstnei.end());//首先加入一阶邻近的编码
	// secnei_vec.assign(firstnei.begin(), firstnei.end());//讲内容复制到vec，但是会删除以前的内容
	for (const auto& ele : firstnei) 
	{
		 vector<DgQ2DICoord> secneiadd  = this->triedgenei(ele);
		 for(const auto& secele : secneiadd)
		 {
			 if (secele == add1)
			 {
				continue;
			 }
			else
			{
				secnei_vec.push_back(secele);
			}
		 }
	}
	// secnei_vec.assign(secnei_set.begin(), secnei_set.end());
	return secnei_vec;
}
// 四邻近源程序有实现不再实现
vector<DgQ2DICoord> DgIDGGBase::fourneicell(DgQ2DICoord &add1) {
    return vector<DgQ2DICoord>();
}
// DgIDGGBase::setAddNeighborsBdry2

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
