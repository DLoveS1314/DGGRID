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
// DgDmdIDGG.cpp: DgDmdIDGG class implementation
//
// Kevin Sahr, 8/12/20
//
////////////////////////////////////////////////////////////////////////////////

#include <dglib/DgBoundedIDGG.h>
#include <dglib/DgDmdD4Grid2D.h>
#include <dglib/DgDmdD4Grid2DS.h>
#include <dglib/DgDmdD8Grid2D.h>
#include <dglib/DgDmdD8Grid2DS.h>
#include <dglib/DgDmdIDGG.h>
#include <dglib/DgIDGGS4D.h>
#include <dglib/DgRadixString.h>
#include <dglib/DgSeriesConverter.h>

#include <cfloat>
#include <climits>
#include <cmath>

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
DgDmdIDGG::DgDmdIDGG (const DgIDGGS4D& dggs, unsigned int aperture,
              int res, const string& name, DgGridMetric gridMetric,
              unsigned int precision)
   : DgIDGGBase (&dggs, dggs.geoRF(), aperture, res, name, Diamond, gridMetric,
                 precision),
	   scaleFac_ (1.0L)
{
   initialize();

} // DgDmdIDGG::DgDmdIDGG

////////////////////////////////////////////////////////////////////////////////
DgDmdIDGG::DgDmdIDGG (const DgDmdIDGG& rfIn)
   : DgIDGGBase (rfIn.dggs(), rfIn.geoRF(), rfIn.aperture(), rfIn.res(),
                 rfIn.name(), Diamond, rfIn.gridMetric(), rfIn.precision()),
	scaleFac_ (rfIn.scaleFac())
{
   initialize();

} // DgDmdIDGG::DgDmdIDGG

////////////////////////////////////////////////////////////////////////////////
DgDmdIDGG::~DgDmdIDGG (void) { }

////////////////////////////////////////////////////////////////////////////////
const DgIDGGS4D&
DgDmdIDGG::dmdDggs (void) const
{ return *(static_cast<const DgIDGGS4D*>(dggs())); }

////////////////////////////连着DGIDGGBase里的参数一起初始化////////////////////////////////////////////////////
void
DgDmdIDGG::initialize (void)
{
   // verify parameter validity

   string apErrStr = string("DgDmdIDGG::initialize(): invalid aperture " +
             dgg::util::to_string(aperture()) + string(" for grid topo ") +
             to_string(gridTopo()));

   if (gridTopo() != Diamond)
      report("DgDmdIDGG::initialize(): invalid grid topo " +
             to_string(gridTopo()), DgBase::Fatal);

   if (aperture() != 4) report(apErrStr, DgBase::Fatal);

   // create some internal data structures
   setUndefLoc(makeLocation(undefAddress()));
   //每一层都有一个
   sphIcosa_ = new DgSphIcosa(vert0(), azDegs());//根据 方位角和指定的经纬度确定球面的方向 并提前计算一些投影转换参数

   isAligned_ = false;
   isCongruent_ = false;

   // initialize parent values as if this is grid res 0 0层无法获取父格网 所以要放在循环外面
   long double parentScaleFac = 1.0;
   unsigned long long int parentNCells = 1;

   // get actual parent values if there is a parent grid 循环创建dgg 需要知道上一层格网信息
   if (res() > 0) {
      const DgDmdIDGG& parentIDGG = dmdDggs().dmdIdgg(res() - 1);

      parentScaleFac = parentIDGG.scaleFac();
      parentNCells = parentIDGG.gridStats().nCells();
   }

   // set-up local network to scale so that quad (and consequently dmd) edge
   // length is 1.0 信息保存在了局部路由里locNet_ ccFrame_单位是 double x double y 属于过渡坐标系坐标
//   DgContCartRF这个类里也没有什么内容 就是求一下距离
   ccFrame_ = DgContCartRF::makeRF(locNet_, name() + "CC1");
   if (gridMetric() == D4)
      grid2DS_ = DgDmdD4Grid2DS::makeRF(locNet_, ccFrame(), res() + 1, 4, true, false, name() + string("D4H2DS"));
   else // must be D8
      grid2DS_ = DgDmdD8Grid2DS::makeRF(locNet_, ccFrame(), res() + 1, 4, true, false, name() + string("D8H2DS"));
   //cout << "== NEW GRID2DS:" << endl;
   //cout << *grid2DS_;

   if (res() == 0)
      maxD_ = 0;
   else {
      double factor = parentScaleFac * 2.0L; // aperture 4 每次二分

      scaleFac_ = factor;
      maxD_ = factor - 1.0L;//行从0开始计算

      //cout << res() << " " << aperture();
      //cout << " f: " << factor << " maxD: " << maxD_ << endl;
   }

   maxI_ = maxD();//返回值是int类型 maxD_是float类型
   maxJ_ = maxD();//菱形的行列大小事一样的 不同网格 算法不同
   mag_ = maxD() + 1;//未发现有什么用 g代表什么意思？在三角形格网有用 计算DgBoundedIDGG的时候offsetPerQuad_ 有用
   firstAdd_ = DgQ2DICoord(1, DgIVec2D(0, 0));
   lastAdd_ = DgQ2DICoord(10, DgIVec2D(maxI(), maxJ()));

   if (res() == 0)
      gridStats_.setNCells(10);
   else
      gridStats_.setNCells(parentNCells * 4);
    //创建内部的转换器，父类的函数 所有的IDGG都要这样实现
   createConverters();

   ///// calculate the statistics /////

   gridStats_.setPrecision(precision());

   long double tmpLen = DgGeoSphRF::icosaEdgeKM();
////// NEEDS UPDATING 每一层参数计算
   gridStats_.setCellDistKM(tmpLen / pow(sqrt((long double) aperture()), res()));

      // a = globeArea / (#cells - 2); 总的地球面积 除以格网数量
      gridStats_.setCellAreaKM(DgGeoSphRF::totalAreaKM() / gridStats_.nCells());

   gridStats_.setCLS(2.0L * 2.0L * DgGeoSphRF::earthRadiusKM() *
                     asinl(sqrt(gridStats_.cellAreaKM() / M_PI) /
                     (2.0L * DgGeoSphRF::earthRadiusKM())));

} // DgDmdIDGG::initialize

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
