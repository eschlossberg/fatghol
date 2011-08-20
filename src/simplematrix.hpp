/**
 * @file   simple.hpp
 *
 * Interface of the matrix class.
 * Make Python wrappper with::
 *   gcc -pthread -fno-strict-aliasing -DNDEBUG -g -fwrapv -O2 -Wall -fPIC -I/usr/include/python2.5 -c simple_wrap.cxx -o simple_wrap.o
 *   g++ -pthread -shared  -llinbox -llinboxsage -Wl,-O1 -Wl,-Bsymbolic-functions simple_wrap.o -o _simple.so
 *
 * @author  riccardo.murri@gmail.com
 * @version $Revision$
 */
/*
 * Copyright (c) 2005, 2006, 2009 riccardo.murri@gmail.com
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA
 *
 */

#ifndef SIMPLE_HPP
#define SIMPLE_HPP


//#include <linbox/blackbox/triplesbb.h>
#include <linbox/field/PID-integer.h>
#include <linbox/blackbox/sparse.h>
#include <linbox/solutions/rank.h>
#include <linbox/solutions/methods.h>


class SimpleMatrix {
public:
  /** Default ctor: make null sparse matrix. */
  SimpleMatrix(size_t const m, size_t const n);  

  /** Add an entry. */
  void addToEntry(size_t const i, size_t const j, 
                  int const value);

  /** Return rank of this matrix. */
  unsigned long rank(void);

  size_t const num_rows;
  size_t const num_columns;

private:
  typedef LinBox::PID_integer _CoefficientRingType;
  _CoefficientRingType ZZ;
  LinBox::SparseMatrix<_CoefficientRingType> m;
};


// --- inline methods ---

inline
SimpleMatrix::SimpleMatrix(size_t const m, size_t const n):
  num_rows(m), num_columns(n), ZZ(), m(ZZ, m, n)
{
  // nothing to do
}


#include <iostream>
inline 
void
SimpleMatrix::addToEntry(size_t const i, size_t const j,
                      int const value)
{
  _CoefficientRingType::Element a;
  ZZ.init(a, value);
  m.refEntry(i, j) += a;
}


inline
unsigned long
SimpleMatrix::rank()
{
  unsigned long r;
  return LinBox::rank(r, m);
}


#endif // SIMPLE_HPP
