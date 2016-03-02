// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2008 Benoit Jacob <jacob.benoit.1@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "main.h"

template<typename MatrixType> void matrixVisitor(const MatrixType& p)
{
  typedef typename MatrixType::Scalar Scalar;

  int rows = p.rows();
  int cols = p.cols();

  // construct a random matrix where all coefficients are different
  MatrixType m;
  m = MatrixType::Random(rows, cols);
  for(int i = 0; i < m.size(); i++)
    for(int i2 = 0; i2 < i; i2++)
      while(m(i) == m(i2)) // yes, ==
        m(i) = ei_random<Scalar>();
  
  Scalar minc = Scalar(1000), maxc = Scalar(-1000);
  int minrow=0,mincol=0,maxrow=0,maxcol=0;
  for(int j = 0; j < cols; j++)
  for(int i = 0; i < rows; i++)
  {
    if(m(i,j) < minc)
    {
      minc = m(i,j);
      minrow = i;
      mincol = j;
    }
    if(m(i,j) > maxc)
    {
      maxc = m(i,j);
      maxrow = i;
      maxcol = j;
    }
  }
  int eigen_minrow, eigen_mincol, eigen_maxrow, eigen_maxcol;
  Scalar eigen_minc, eigen_maxc;
  eigen_minc = m.minCoeff(&eigen_minrow,&eigen_mincol);
  eigen_maxc = m.maxCoeff(&eigen_maxrow,&eigen_maxcol);
  VERIFY(minrow == eigen_minrow);
  VERIFY(maxrow == eigen_maxrow);
  VERIFY(mincol == eigen_mincol);
  VERIFY(maxcol == eigen_maxcol);
  VERIFY_IS_APPROX(minc, eigen_minc);
  VERIFY_IS_APPROX(maxc, eigen_maxc);
  VERIFY_IS_APPROX(minc, m.minCoeff());
  VERIFY_IS_APPROX(maxc, m.maxCoeff());
}

template<typename VectorType> void vectorVisitor(const VectorType& w)
{
  typedef typename VectorType::Scalar Scalar;

  int size = w.size();

  // construct a random vector where all coefficients are different
  VectorType v;
  v = VectorType::Random(size);
  for(int i = 0; i < size; i++)
    for(int i2 = 0; i2 < i; i2++)
      while(v(i) == v(i2)) // yes, ==
        v(i) = ei_random<Scalar>();
  
  Scalar minc = Scalar(1000), maxc = Scalar(-1000);
  int minidx=0,maxidx=0;
  for(int i = 0; i < size; i++)
  {
    if(v(i) < minc)
    {
      minc = v(i);
      minidx = i;
    }
    if(v(i) > maxc)
    {
      maxc = v(i);
      maxidx = i;
    }
  }
  int eigen_minidx, eigen_maxidx;
  Scalar eigen_minc, eigen_maxc;
  eigen_minc = v.minCoeff(&eigen_minidx);
  eigen_maxc = v.maxCoeff(&eigen_maxidx);
  VERIFY(minidx == eigen_minidx);
  VERIFY(maxidx == eigen_maxidx);
  VERIFY_IS_APPROX(minc, eigen_minc);
  VERIFY_IS_APPROX(maxc, eigen_maxc);
  VERIFY_IS_APPROX(minc, v.minCoeff());
  VERIFY_IS_APPROX(maxc, v.maxCoeff());
}

void test_eigen2_visitor()
{
  for(int i = 0; i < g_repeat; i++) {
    CALL_SUBTEST_1( matrixVisitor(Matrix<float, 1, 1>()) );
    CALL_SUBTEST_2( matrixVisitor(Matrix2f()) );
    CALL_SUBTEST_3( matrixVisitor(Matrix4d()) );
    CALL_SUBTEST_4( matrixVisitor(MatrixXd(8, 12)) );
    CALL_SUBTEST_5( matrixVisitor(Matrix<double,Dynamic,Dynamic,RowMajor>(20, 20)) );
    CALL_SUBTEST_6( matrixVisitor(MatrixXi(8, 12)) );
  }
  for(int i = 0; i < g_repeat; i++) {
    CALL_SUBTEST_7( vectorVisitor(Vector4f()) );
    CALL_SUBTEST_4( vectorVisitor(VectorXd(10)) );
    CALL_SUBTEST_4( vectorVisitor(RowVectorXd(10)) );
    CALL_SUBTEST_8( vectorVisitor(VectorXf(33)) );
  }
}
