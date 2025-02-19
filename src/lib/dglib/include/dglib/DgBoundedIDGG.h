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
// DgBoundedIDGG.h: DgBoundedIDGG class definitions
//
// Version 7.0 - Kevin Sahr, 12/14/14
// Version 6.1 - Kevin Sahr, 5/23/13
//
////////////////////////////////////////////////////////////////////////////////

#ifndef DGBOUNDEDIDGG_H
#define DGBOUNDEDIDGG_H

#include <dglib/DgBoundedRF.h>
#include <dglib/DgBoundedRF2D.h>
#include <dglib/DgGeoSphRF.h>
#include <dglib/DgIDGGBase.h>
#include <dglib/DgIDGGutil.h>

using namespace std;

////////////////////////////////////////////////////////////////////////////////
class DgBoundedIDGG : public DgBoundedRF<DgQ2DICoord, DgGeoCoord, long double> {

   public:

      DgBoundedIDGG (const DgIDGGBase& IDGGin);
     ~DgBoundedIDGG (void) { delete bnd2D_; }
//    先把列加1 因为是行优先 然后行加1 载然后 qnum加1
      virtual DgQ2DICoord& incrementAddress (DgQ2DICoord& add) const;
      virtual DgQ2DICoord& decrementAddress (DgQ2DICoord& add) const;

      virtual bool validAddress (const DgQ2DICoord& add) const;

      const DgQ2DICoord& invalidAdd (void) const
                         { return idgg().undefAddress(); }

      const DgIDGGBase& idgg (void) const { return IDGG_; }
      // 通过记录这个数量 免去记录res，其实记录res更好 但是不符合 Dgbounded 的逻辑，应该Dgboundeds才有res
      unsigned long long int offsetPerQuad (void) const { return offsetPerQuad_; }
//        行列号转seqnum
      virtual unsigned long long int seqNumAddress (const DgQ2DICoord& add) const;
//        seqnum转行列号
      virtual DgQ2DICoord addFromSeqNum (unsigned long long int sNum) const;

      virtual DgQ2DICoord q2dixToQ2di (const DgQ2DICoord& add) const;

      virtual operator string (void) const
        {
           string s = "=== DgBoundedIDGG: " + DgBoundedRF::operator string();
           s += "\n offsetPerQuad: " + dgg::util::to_string(offsetPerQuad());
           s += "\n BND2D: " + string(*bnd2D_);

           return s;
        }


   //protected:

      const DgBoundedRF2D& bnd2D (void) const { return *bnd2D_; }

   private:

      const DgIDGGBase& IDGG_;
      DgBoundedRF2D* bnd2D_;
      unsigned long long int offsetPerQuad_;
};

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
#endif
