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
 * Copyright (c) 2009-2011 riccardo.murri@gmail.com
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

#include <cassert>
#include <fstream>
#include <string>
#include <iostream>

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

  /** Return the value of entry at row @a i and column @j */
  int getEntry(size_t const i, size_t const j) const;

  /** Return rank of this matrix. */
  unsigned long rank(void);

  size_t const num_rows;
  size_t const num_columns;

  /** Return `true` if product of given matrices is null */
  friend bool is_null_product(const SimpleMatrix &M1, 
                              const SimpleMatrix &M2);

  /** Dump entries to named file */
  void save(const char *const filename);

  /** Load entries from named file.  Return true on successful load, and false on error. */
  bool load(const char *const filename);

private:
  typedef LinBox::PID_integer _CoefficientRingType;
  typedef LinBox::SparseMatrix<_CoefficientRingType> _MatrixType;
  _CoefficientRingType ZZ;
  _MatrixType m;
};


// --- inline methods ---

inline
SimpleMatrix::SimpleMatrix(size_t const m, size_t const n):
  num_rows(m), num_columns(n), ZZ(), m(ZZ, m, n)
{
  // nothing to do
}


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
int
SimpleMatrix::getEntry(size_t const i, size_t const j) const
{
  if (num_rows == 0 or num_columns == 0)
    return 0;
  else {
    double result = 0;
    ZZ.convert(result, m.getEntry(i, j));
    return result;
  };
}


inline
unsigned long
SimpleMatrix::rank()
{
  if (num_rows == 0 or num_columns == 0)
    return 0;
  else {
#if defined(FATGHOL_USE_LINBOX_ELIMINATION_PIVOT_LINEAR)
    LinBox::Method::SparseElimination se;
    se.strategy(LinBox::Specifier::PIVOT_LINEAR);
#elif defined(FATGHOL_USE_LINBOX_ELIMINATION_PIVOT_NONE)
    LinBox::Method::SparseElimination se;
    se.strategy(LinBox::Specifier::PIVOT_NONE);
#endif

    unsigned long r;
#if defined(FATGHOL_USE_LINBOX_ELIMINATION_PIVOT_LINEAR) || defined(FATGHOL_USE_LINBOX_ELIMINATION_PIVOT_NONE)
    // use the prescribed elimination strategy
    return LinBox::rank(r, m, se);
#elif defined(FATGHOL_USE_LINBOX_DEFAULT)
    // use LinBox' default strategy (blackbox, as of 1.1.7)
    return LinBox::rank(r, m);
#else
# ifndef SWIG
#  error Rheinfall support not yet implemented.  Please define FATGHOL_USE_LINBOX_ELIMINATION_PIVOT_LINEAR.
# endif
#endif  
  }
}


inline
bool
is_null_product(const SimpleMatrix &A, const SimpleMatrix &B)
{
  assert (A.m.coldim() == B.m.rowdim());
  SimpleMatrix::_CoefficientRingType ZZ;
  for(size_t i=0; i<A.m.rowdim(); i++) {
    for (size_t j=0; j<B.m.coldim(); j++) {
      SimpleMatrix::_CoefficientRingType::Element x;
      ZZ.init(x, 0);
      for (size_t k=0; k<A.m.coldim(); k++) 
        x += A.m.getEntry(i,k) * B.m.getEntry(k,j);
      if (not (x == 0))
        return false;
    }
  }
  return true;
}


inline
bool
SimpleMatrix::load(const char *const filename)
{
  std::ifstream input(filename);
  if ((not input.is_open()) or input.bad()) 
    return false;

  m.read(input, LinBox::FORMAT_DETECT);

  // all done
  input.close();
  return true;
}


inline
void
SimpleMatrix::save(const char *const filename)
{
  std::ofstream output(filename);
  m.write(output, LinBox::FORMAT_GUILLAUME);

  // all done
  output.close();
}


#endif // SIMPLE_HPP
