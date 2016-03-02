// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2008 Gael Guennebaud <gael.guennebaud@inria.fr>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef EIGEN_NO_ASSERTION_CHECKING
#define EIGEN_NO_ASSERTION_CHECKING
#endif

static int nb_temporaries;

#define EIGEN_DENSE_STORAGE_CTOR_PLUGIN { if(size!=0) nb_temporaries++; }

#include "main.h"
#include <Eigen/Cholesky>
#include <Eigen/QR>

#define VERIFY_EVALUATION_COUNT(XPR,N) {\
    nb_temporaries = 0; \
    XPR; \
    if(nb_temporaries!=N) std::cerr << "nb_temporaries == " << nb_temporaries << "\n"; \
    VERIFY( (#XPR) && nb_temporaries==N ); \
  }

template<typename MatrixType,template <typename,int> class CholType> void test_chol_update(const MatrixType& symm)
{
  typedef typename MatrixType::Scalar Scalar;
  typedef typename MatrixType::RealScalar RealScalar;
  typedef Matrix<Scalar, MatrixType::RowsAtCompileTime, 1> VectorType;

  MatrixType symmLo = symm.template triangularView<Lower>();
  MatrixType symmUp = symm.template triangularView<Upper>();
  MatrixType symmCpy = symm;

  CholType<MatrixType,Lower> chollo(symmLo);
  CholType<MatrixType,Upper> cholup(symmUp);

  for (int k=0; k<10; ++k)
  {
    VectorType vec = VectorType::Random(symm.rows());
    RealScalar sigma = internal::random<RealScalar>();
    symmCpy += sigma * vec * vec.adjoint();

    // we are doing some downdates, so it might be the case that the matrix is not SPD anymore
    CholType<MatrixType,Lower> chol(symmCpy);
    if(chol.info()!=Success)
      break;

    chollo.rankUpdate(vec, sigma);
    VERIFY_IS_APPROX(symmCpy, chollo.reconstructedMatrix());

    cholup.rankUpdate(vec, sigma);
    VERIFY_IS_APPROX(symmCpy, cholup.reconstructedMatrix());
  }
}

template<typename MatrixType> void cholesky(const MatrixType& m)
{
  typedef typename MatrixType::Index Index;
  /* this test covers the following files:
     LLT.h LDLT.h
  */
  Index rows = m.rows();
  Index cols = m.cols();

  typedef typename MatrixType::Scalar Scalar;
  typedef typename NumTraits<Scalar>::Real RealScalar;
  typedef Matrix<Scalar, MatrixType::RowsAtCompileTime, MatrixType::RowsAtCompileTime> SquareMatrixType;
  typedef Matrix<Scalar, MatrixType::RowsAtCompileTime, 1> VectorType;

  MatrixType a0 = MatrixType::Random(rows,cols);
  VectorType vecB = VectorType::Random(rows), vecX(rows);
  MatrixType matB = MatrixType::Random(rows,cols), matX(rows,cols);
  SquareMatrixType symm =  a0 * a0.adjoint();
  // let's make sure the matrix is not singular or near singular
  for (int k=0; k<3; ++k)
  {
    MatrixType a1 = MatrixType::Random(rows,cols);
    symm += a1 * a1.adjoint();
  }

  // to test if really Cholesky only uses the upper triangular part, uncomment the following
  // FIXME: currently that fails !!
  //symm.template part<StrictlyLower>().setZero();

  {
    SquareMatrixType symmUp = symm.template triangularView<Upper>();
    SquareMatrixType symmLo = symm.template triangularView<Lower>();
    
    LLT<SquareMatrixType,Lower> chollo(symmLo);
    VERIFY_IS_APPROX(symm, chollo.reconstructedMatrix());
    vecX = chollo.solve(vecB);
    VERIFY_IS_APPROX(symm * vecX, vecB);
    matX = chollo.solve(matB);
    VERIFY_IS_APPROX(symm * matX, matB);

    // test the upper mode
    LLT<SquareMatrixType,Upper> cholup(symmUp);
    VERIFY_IS_APPROX(symm, cholup.reconstructedMatrix());
    vecX = cholup.solve(vecB);
    VERIFY_IS_APPROX(symm * vecX, vecB);
    matX = cholup.solve(matB);
    VERIFY_IS_APPROX(symm * matX, matB);

    MatrixType neg = -symmLo;
    chollo.compute(neg);
    VERIFY(chollo.info()==NumericalIssue);

    VERIFY_IS_APPROX(MatrixType(chollo.matrixL().transpose().conjugate()), MatrixType(chollo.matrixU()));
    VERIFY_IS_APPROX(MatrixType(chollo.matrixU().transpose().conjugate()), MatrixType(chollo.matrixL()));
    VERIFY_IS_APPROX(MatrixType(cholup.matrixL().transpose().conjugate()), MatrixType(cholup.matrixU()));
    VERIFY_IS_APPROX(MatrixType(cholup.matrixU().transpose().conjugate()), MatrixType(cholup.matrixL()));
    
    // test some special use cases of SelfCwiseBinaryOp:
    MatrixType m1 = MatrixType::Random(rows,cols), m2(rows,cols);
    m2 = m1;
    m2 += symmLo.template selfadjointView<Lower>().llt().solve(matB);
    VERIFY_IS_APPROX(m2, m1 + symmLo.template selfadjointView<Lower>().llt().solve(matB));
    m2 = m1;
    m2 -= symmLo.template selfadjointView<Lower>().llt().solve(matB);
    VERIFY_IS_APPROX(m2, m1 - symmLo.template selfadjointView<Lower>().llt().solve(matB));
    m2 = m1;
    m2.noalias() += symmLo.template selfadjointView<Lower>().llt().solve(matB);
    VERIFY_IS_APPROX(m2, m1 + symmLo.template selfadjointView<Lower>().llt().solve(matB));
    m2 = m1;
    m2.noalias() -= symmLo.template selfadjointView<Lower>().llt().solve(matB);
    VERIFY_IS_APPROX(m2, m1 - symmLo.template selfadjointView<Lower>().llt().solve(matB));
  }

  // LDLT
  {
    int sign = internal::random<int>()%2 ? 1 : -1;

    if(sign == -1)
    {
      symm = -symm; // test a negative matrix
    }

    SquareMatrixType symmUp = symm.template triangularView<Upper>();
    SquareMatrixType symmLo = symm.template triangularView<Lower>();

    LDLT<SquareMatrixType,Lower> ldltlo(symmLo);
    VERIFY_IS_APPROX(symm, ldltlo.reconstructedMatrix());
    vecX = ldltlo.solve(vecB);
    VERIFY_IS_APPROX(symm * vecX, vecB);
    matX = ldltlo.solve(matB);
    VERIFY_IS_APPROX(symm * matX, matB);

    LDLT<SquareMatrixType,Upper> ldltup(symmUp);
    VERIFY_IS_APPROX(symm, ldltup.reconstructedMatrix());
    vecX = ldltup.solve(vecB);
    VERIFY_IS_APPROX(symm * vecX, vecB);
    matX = ldltup.solve(matB);
    VERIFY_IS_APPROX(symm * matX, matB);

    VERIFY_IS_APPROX(MatrixType(ldltlo.matrixL().transpose().conjugate()), MatrixType(ldltlo.matrixU()));
    VERIFY_IS_APPROX(MatrixType(ldltlo.matrixU().transpose().conjugate()), MatrixType(ldltlo.matrixL()));
    VERIFY_IS_APPROX(MatrixType(ldltup.matrixL().transpose().conjugate()), MatrixType(ldltup.matrixU()));
    VERIFY_IS_APPROX(MatrixType(ldltup.matrixU().transpose().conjugate()), MatrixType(ldltup.matrixL()));

    if(MatrixType::RowsAtCompileTime==Dynamic)
    {
      // note : each inplace permutation requires a small temporary vector (mask)

      // check inplace solve
      matX = matB;
      VERIFY_EVALUATION_COUNT(matX = ldltlo.solve(matX), 0);
      VERIFY_IS_APPROX(matX, ldltlo.solve(matB).eval());


      matX = matB;
      VERIFY_EVALUATION_COUNT(matX = ldltup.solve(matX), 0);
      VERIFY_IS_APPROX(matX, ldltup.solve(matB).eval());
    }

    // restore
    if(sign == -1)
      symm = -symm;

    // check matrices coming from linear constraints with Lagrange multipliers
    if(rows>=3)
    {
      SquareMatrixType A = symm;
      int c = internal::random<int>(0,rows-2);
      A.bottomRightCorner(c,c).setZero();
      // Make sure a solution exists:
      vecX.setRandom();
      vecB = A * vecX;
      vecX.setZero();
      ldltlo.compute(A);
      VERIFY_IS_APPROX(A, ldltlo.reconstructedMatrix());
      vecX = ldltlo.solve(vecB);
      VERIFY_IS_APPROX(A * vecX, vecB);
    }
    
    // check non-full rank matrices
    if(rows>=3)
    {
      int r = internal::random<int>(1,rows-1);
      Matrix<Scalar,Dynamic,Dynamic> a = Matrix<Scalar,Dynamic,Dynamic>::Random(rows,r);
      SquareMatrixType A = a * a.adjoint();
      // Make sure a solution exists:
      vecX.setRandom();
      vecB = A * vecX;
      vecX.setZero();
      ldltlo.compute(A);
      VERIFY_IS_APPROX(A, ldltlo.reconstructedMatrix());
      vecX = ldltlo.solve(vecB);
      VERIFY_IS_APPROX(A * vecX, vecB);
    }
    
    // check matrices with a wide spectrum
    if(rows>=3)
    {
      RealScalar s = (std::min)(16,std::numeric_limits<RealScalar>::max_exponent10/8);
      Matrix<Scalar,Dynamic,Dynamic> a = Matrix<Scalar,Dynamic,Dynamic>::Random(rows,rows);
      Matrix<RealScalar,Dynamic,1> d =  Matrix<RealScalar,Dynamic,1>::Random(rows);
      for(int k=0; k<rows; ++k)
        d(k) = d(k)*std::pow(RealScalar(10),internal::random<RealScalar>(-s,s));
      SquareMatrixType A = a * d.asDiagonal() * a.adjoint();
      // Make sure a solution exists:
      vecX.setRandom();
      vecB = A * vecX;
      vecX.setZero();
      ldltlo.compute(A);
      VERIFY_IS_APPROX(A, ldltlo.reconstructedMatrix());
      vecX = ldltlo.solve(vecB);
      VERIFY_IS_APPROX(A * vecX, vecB);
    }
  }

  // update/downdate
  CALL_SUBTEST(( test_chol_update<SquareMatrixType,LLT>(symm)  ));
  CALL_SUBTEST(( test_chol_update<SquareMatrixType,LDLT>(symm) ));
}

template<typename MatrixType> void cholesky_cplx(const MatrixType& m)
{
  // classic test
  cholesky(m);

  // test mixing real/scalar types

  typedef typename MatrixType::Index Index;

  Index rows = m.rows();
  Index cols = m.cols();

  typedef typename MatrixType::Scalar Scalar;
  typedef typename NumTraits<Scalar>::Real RealScalar;
  typedef Matrix<RealScalar, MatrixType::RowsAtCompileTime, MatrixType::RowsAtCompileTime> RealMatrixType;
  typedef Matrix<Scalar, MatrixType::RowsAtCompileTime, 1> VectorType;

  RealMatrixType a0 = RealMatrixType::Random(rows,cols);
  VectorType vecB = VectorType::Random(rows), vecX(rows);
  MatrixType matB = MatrixType::Random(rows,cols), matX(rows,cols);
  RealMatrixType symm =  a0 * a0.adjoint();
  // let's make sure the matrix is not singular or near singular
  for (int k=0; k<3; ++k)
  {
    RealMatrixType a1 = RealMatrixType::Random(rows,cols);
    symm += a1 * a1.adjoint();
  }

  {
    RealMatrixType symmLo = symm.template triangularView<Lower>();

    LLT<RealMatrixType,Lower> chollo(symmLo);
    VERIFY_IS_APPROX(symm, chollo.reconstructedMatrix());
    vecX = chollo.solve(vecB);
    VERIFY_IS_APPROX(symm * vecX, vecB);
//     matX = chollo.solve(matB);
//     VERIFY_IS_APPROX(symm * matX, matB);
  }

  // LDLT
  {
    int sign = internal::random<int>()%2 ? 1 : -1;

    if(sign == -1)
    {
      symm = -symm; // test a negative matrix
    }

    RealMatrixType symmLo = symm.template triangularView<Lower>();

    LDLT<RealMatrixType,Lower> ldltlo(symmLo);
    VERIFY_IS_APPROX(symm, ldltlo.reconstructedMatrix());
    vecX = ldltlo.solve(vecB);
    VERIFY_IS_APPROX(symm * vecX, vecB);
//     matX = ldltlo.solve(matB);
//     VERIFY_IS_APPROX(symm * matX, matB);
  }
}

// regression test for bug 241
template<typename MatrixType> void cholesky_bug241(const MatrixType& m)
{
  eigen_assert(m.rows() == 2 && m.cols() == 2);

  typedef typename MatrixType::Scalar Scalar;
  typedef Matrix<Scalar, MatrixType::RowsAtCompileTime, 1> VectorType;

  MatrixType matA;
  matA << 1, 1, 1, 1;
  VectorType vecB;
  vecB << 1, 1;
  VectorType vecX = matA.ldlt().solve(vecB);
  VERIFY_IS_APPROX(matA * vecX, vecB);
}

// LDLT is not guaranteed to work for indefinite matrices, but happens to work fine if matrix is diagonal.
// This test checks that LDLT reports correctly that matrix is indefinite. 
// See http://forum.kde.org/viewtopic.php?f=74&t=106942 and bug 736
template<typename MatrixType> void cholesky_definiteness(const MatrixType& m)
{
  eigen_assert(m.rows() == 2 && m.cols() == 2);
  MatrixType mat;
  LDLT<MatrixType> ldlt(2);
  
  {
    mat << 1, 0, 0, -1;
    ldlt.compute(mat);
    VERIFY(!ldlt.isNegative());
    VERIFY(!ldlt.isPositive());
  }
  {
    mat << 1, 2, 2, 1;
    ldlt.compute(mat);
    VERIFY(!ldlt.isNegative());
    VERIFY(!ldlt.isPositive());
  }
  {
    mat << 0, 0, 0, 0;
    ldlt.compute(mat);
    VERIFY(ldlt.isNegative());
    VERIFY(ldlt.isPositive());
  }
  {
    mat << 0, 0, 0, 1;
    ldlt.compute(mat);
    VERIFY(!ldlt.isNegative());
    VERIFY(ldlt.isPositive());
  }
  {
    mat << -1, 0, 0, 0;
    ldlt.compute(mat);
    VERIFY(ldlt.isNegative());
    VERIFY(!ldlt.isPositive());
  }
}

template<typename MatrixType> void cholesky_verify_assert()
{
  MatrixType tmp;

  LLT<MatrixType> llt;
  VERIFY_RAISES_ASSERT(llt.matrixL())
  VERIFY_RAISES_ASSERT(llt.matrixU())
  VERIFY_RAISES_ASSERT(llt.solve(tmp))
  VERIFY_RAISES_ASSERT(llt.solveInPlace(&tmp))

  LDLT<MatrixType> ldlt;
  VERIFY_RAISES_ASSERT(ldlt.matrixL())
  VERIFY_RAISES_ASSERT(ldlt.permutationP())
  VERIFY_RAISES_ASSERT(ldlt.vectorD())
  VERIFY_RAISES_ASSERT(ldlt.isPositive())
  VERIFY_RAISES_ASSERT(ldlt.isNegative())
  VERIFY_RAISES_ASSERT(ldlt.solve(tmp))
  VERIFY_RAISES_ASSERT(ldlt.solveInPlace(&tmp))
}

void test_cholesky()
{
  int s = 0;
  for(int i = 0; i < g_repeat; i++) {
    CALL_SUBTEST_1( cholesky(Matrix<double,1,1>()) );
    CALL_SUBTEST_3( cholesky(Matrix2d()) );
    CALL_SUBTEST_3( cholesky_bug241(Matrix2d()) );
    CALL_SUBTEST_3( cholesky_definiteness(Matrix2d()) );
    CALL_SUBTEST_4( cholesky(Matrix3f()) );
    CALL_SUBTEST_5( cholesky(Matrix4d()) );
    s = internal::random<int>(1,EIGEN_TEST_MAX_SIZE);
    CALL_SUBTEST_2( cholesky(MatrixXd(s,s)) );
    s = internal::random<int>(1,EIGEN_TEST_MAX_SIZE/2);
    CALL_SUBTEST_6( cholesky_cplx(MatrixXcd(s,s)) );
  }

  CALL_SUBTEST_4( cholesky_verify_assert<Matrix3f>() );
  CALL_SUBTEST_7( cholesky_verify_assert<Matrix3d>() );
  CALL_SUBTEST_8( cholesky_verify_assert<MatrixXf>() );
  CALL_SUBTEST_2( cholesky_verify_assert<MatrixXd>() );

  // Test problem size constructors
  CALL_SUBTEST_9( LLT<MatrixXf>(10) );
  CALL_SUBTEST_9( LDLT<MatrixXf>(10) );
  
  TEST_SET_BUT_UNUSED_VARIABLE(s)
  TEST_SET_BUT_UNUSED_VARIABLE(nb_temporaries)
}
