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
#include <cstdint>
#include <fstream>
#include <string>
#include <iostream>

#ifdef FATGHOL_USE_RHEINFALL
# include <rheinfall/rank.hpp>
# include <map>  // std::map
# include <util> // std::pair
#else // use LinBox
# include <linbox/field/PID-integer.h>
# include <linbox/blackbox/sparse.h>
# include <linbox/solutions/rank.h>
# include <linbox/solutions/methods.h>
#endif


class SimpleMatrix {
public:
  /** Default ctor: make null sparse matrix. */
  SimpleMatrix(int32_t const m, int32_t const n);  

  /** Add an entry. */
  void addToEntry(int32_t const i, int32_t const j, 
                  int const value);

  /** Return the value of entry at row @a i and column @j */
  int getEntry(int32_t const i, int32_t const j) const;

  /** Return rank of this matrix. */
  unsigned int32_t rank(void);

  int32_t const num_rows;
  int32_t const num_columns;

  /** Return `true` if product of given matrices is null */
  friend bool is_null_product(const SimpleMatrix &M1, 
                              const SimpleMatrix &M2);

  /** Dump entries to named file */
  void save(const char *const filename);

  /** Load entries from named file.  Return true on successful load, and false on error. */
  bool load(const char *const filename);

private:
#ifdef FATGHOL_USE_RHEINFALL
  typedef std::pair< const int32_t,int64_t >                  _coord_and_val;
  typedef std::map< int32_t, int64_t, std::less<int32_t> >    _simplerow;
  typedef std::pair< const int32_t,_simplerow >               _coord_and_simplerow;
  typedef std::map< int32_t, _simplerow, std::less<int32_t> > _simplerows;
  _simplerows m;
#else
  typedef LinBox::PID_integer _CoefficientRingType;
  typedef LinBox::SparseMatrix<_CoefficientRingType> _MatrixType;
  _CoefficientRingType ZZ;
  _MatrixType m;
#endif
};


// --- inline methods ---

inline
SimpleMatrix::SimpleMatrix(int32_t const m, int32_t const n):
  num_rows(m), num_columns(n), 
#ifdef FATGHOL_USE_RHEINFALL
#  error Rheinfall support not yet implemented.  Please define FATGHOL_USE_LINBOX_ELIMINATION_PIVOT_LINEAR.
#else
  ZZ(), m(ZZ, m, n)
#endif
{
  // nothing to do
}


inline 
void
SimpleMatrix::addToEntry(int32_t const i, int32_t const j,
                      int const value)
{
#ifdef FATGHOL_USE_RHEINFALL
    if (0 == m.count(i))
      m[i][j] = value;
    else {
      if (0 == m[i].count(j))
        m[i][j] = value;
      else
        m[i][j] += value;
#else
  _CoefficientRingType::Element a;
  ZZ.init(a, value);
  m.refEntry(i, j) += a;
#endif 
}


inline 
int
SimpleMatrix::getEntry(int32_t const i, int32_t const j) const
{
  if (num_rows == 0 or num_columns == 0)
    return 0;
  else {
#ifdef FATGHOL_USE_RHEINFALL
    if (0 == m.count(i))
      return 0;
    else {
      if (0 == m[i].count(j))
        return 0;
      else
        return m[i][j];
    };
#else // use LinBox
    double result = 0;
    ZZ.convert(result, m.getEntry(i, j));
    return result;
#endif
  };
}


inline
unsigned int32_t
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

    unsigned int32_t r;
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
  assert (A.num_columns == B.num_rows);
#ifdef FATGHOL_USE_RHEINFALL
  for(int32_t i=0; i<A.num_rows; i++) {
    for (int32_t j=0; j<B.num_columns; j++) {
      int64_t x = 0;
      for (int32_t k=0; k<A.num_columns; k++) 
        x += A.m.getEntry(i,k) * B.m.getEntry(k,j);
      if (not (x == 0))
        return false;
    };
  };
#else // use LinBox
  assert (A.m.coldim() == B.m.rowdim());
  SimpleMatrix::_CoefficientRingType ZZ;
  for(int32_t i=0; i<A.m.rowdim(); i++) {
    for (int32_t j=0; j<B.m.coldim(); j++) {
      SimpleMatrix::_CoefficientRingType::Element x;
      ZZ.init(x, 0);
      for (int32_t k=0; k<A.m.coldim(); k++) 
        x += A.m.getEntry(i,k) * B.m.getEntry(k,j);
      if (not (x == 0))
        return false;
    }
  }
#endif
  return true;
}


inline
bool
SimpleMatrix::load(const char *const filename)
{
  std::ifstream input(filename);
  if ((not input.is_open()) or input.bad()) 
    return false;

#ifdef FATGHOL_USE_RHEINFALL
  // XXX: this is basically ripped off Rheinfall's "rank.hpp"
  int32_t nrows, ncols; 
  char M;
  input >> nrows >> ncols >> M;
  if (input.fail() or 'M' != M)
    throw std::runtime_error("Cannot read SMS header");
  if (nrows != num_rows)
    throw std::runtime_error("Number of rows in SMS header does not match"
                             " number of rows passed to constructor");
  if (ncols != num_columns)
    throw std::runtime_error("Number of columns in SMS header does not match"
                             " number of columns passed to constructor");
  int32_t i, j;
  int64_t value;
  while (not input.eof()) {
    input >> i >> j >> value;
    if (0 == i and 0 == j and 0 == value)
      break; // end of matrix stream
    // SMS indices are 1-based
    --i;
    --j;
    // ignore zero entries in matrix -- they shouldn't be here in the first place
    if (0 == value) 
      continue; 
    m[i][j] = value;
  }; // while(not eof)
#else // use LinBox
  m.read(input, LinBox::FORMAT_DETECT);
#endif

  // all done
  input.close();
  return true;
}


inline
void
SimpleMatrix::save(const char *const filename)
{
  std::ofstream output(filename);

#ifdef FATGHOL_USE_RHEINFALL
  output << num_rows <<" "<< num_cols <<" "<< "M" << std::endl;
  for(_simplerows::const_iterator r = m.begin(); r != m.end(); ++r)
    for(_simplerow::const_iterator c = r->second.begin(); c != r->second.end(); ++c)
      output << r->first <<" "<< c->first <<" "<< c->second << std::endl;
  output << "0 0 0" << std::endl;
#else // use LinBox
  m.write(output, LinBox::FORMAT_GUILLAUME);
#endif

  // all done
  output.close();
}


#endif // SIMPLE_HPP
