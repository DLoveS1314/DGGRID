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
// DgAddress.h: DgAddress class definitions
//
////////////////////////////////////////////////////////////////////////////////

#ifndef DGADDRESS_H
#define DGADDRESS_H

#include <dglib/DgAddressBase.h>

class DgDistanceBase;

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////吧坐标类A如DgGeoCoord统一包装成 Dgaddress
/// 实现一些统一的接口如 buildLocation功能中接受这个address参数 可以实现所有对应坐标类的转换//
/// DgRF<A, D>::makeLocation
/// ////////////////////////////////////////////////
template <class A> class DgAddress : public DgAddressBase {

   public:

      DgAddress (void) {}
      DgAddress (const DgAddress<A>& add) : address_ (add.address()) {}
      DgAddress (const A& address) : address_ (address) {}

      void setAddress (const A& address) { address_ = address; }

      A& address (void) { return address_; }

      const A& address (void) const { return address_; }

      DgAddress<A>& operator= (const DgAddress<A>& add) 
                                    { address_ = add.address(); return *this; }

   protected:

      virtual ostream& writeTo (ostream& stream) const
                                    { return stream << address_; }

   private:

      A address_;

};

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

#endif
