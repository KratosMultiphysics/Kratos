//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:       Massimo Petracca $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:           September 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_SHELLQ4_LOCAL_COORDINATE_SYSTEM_H_INCLUDED)
#define  KRATOS_SHELLQ4_LOCAL_COORDINATE_SYSTEM_H_INCLUDED


namespace Kratos
{

/** \brief ShellQ4_LocalCoordinateSystem
 *
 * This class represent the local coordinate system of any element whose geometry
 * is a QUAD 4 in 3D space
 */
class ShellQ4_LocalCoordinateSystem
{

 public:

  KRATOS_CLASS_POINTER_DEFINITION( ShellQ4_LocalCoordinateSystem );

  typedef double RealType;

  typedef array_1d<RealType, 3> Vector3Type;

  typedef std::vector<Vector3Type> Vector3ContainerType;

  typedef Matrix MatrixType;

  typedef matrix_row< const MatrixType > ConstMatrixRowType;

 private:

  Vector3ContainerType mP;
  Vector3Type mCenter;
  MatrixType mOrientation;
  RealType mArea;

 private:

  template< class TVector3 >
  inline static double NormalizeVector3(TVector3 & V)
  {
    double vn = V(0) * V(0) + V(1) * V(1) + V(2) * V(2);
    if(vn != 0.0 && vn != 1.0) {
      vn = std::sqrt(vn);
      V /= vn;
    }
    return vn;
  }

 public:

  ShellQ4_LocalCoordinateSystem(const Vector3Type & P1global,
                                const Vector3Type & P2global,
                                const Vector3Type & P3global,
                                const Vector3Type & P4global)
      : mP(4)
      , mOrientation(3, 3)
  {
    // Form the basis vectors alligning the local X direction
    // with the 1-2 side

    // compute the central point
    noalias( mCenter )  = P1global;
    noalias( mCenter ) += P2global;
    noalias( mCenter ) += P3global;
    noalias( mCenter ) += P4global;
    mCenter *= 0.25;

    // compute the diagonal vectors
    Vector3Type d13( P3global - P1global );
    Vector3Type d24( P4global - P2global );

    // compute the Normal vector at the element center
    // as the cross product of the 2 diagonals.
    // While normalizing the normal vector save its norm to compute the area.
    // Note: the norm should be divided by 2 because its computed from the
    // cross product of the diagonals, which gives twice the area!
    Vector3Type e3;
    MathUtils<RealType>::CrossProduct(e3,        d13, d24);
    mArea = NormalizeVector3( e3 );
    mArea /= 2.0;

    // compute the local X direction parallel to the projection of the side 1-2 onto
    // the local xy plane.
    Vector3Type e1( P2global - P1global );
    RealType e1_dot_e3 = MathUtils<RealType>::Dot(e1, e3);
    noalias( e1 ) -= e1_dot_e3 * e3;
    NormalizeVector3( e1 );

    // finally compute the local Y direction to be orthogonal to both X and Z local directions
    Vector3Type e2;
    MathUtils<RealType>::CrossProduct(e2,    e3, e1);
    NormalizeVector3( e2 );

    // form the 3x3 transformation matrix
    for(int i = 0; i < 3; i++)
    {
      mOrientation(0, i) = e1(i);
      mOrientation(1, i) = e2(i);
      mOrientation(2, i) = e3(i);
    }

    // transform global coordinates to the local coordinate system
    for(int i = 0; i < 3; i++)
    {
      mP[0](i) = mOrientation(i, 0) * (P1global(0) - mCenter(0)) + mOrientation(i, 1) * (P1global(1) - mCenter(1)) + mOrientation(i, 2) * (P1global(2) - mCenter(2));
      mP[1](i) = mOrientation(i, 0) * (P2global(0) - mCenter(0)) + mOrientation(i, 1) * (P2global(1) - mCenter(1)) + mOrientation(i, 2) * (P2global(2) - mCenter(2));
      mP[2](i) = mOrientation(i, 0) * (P3global(0) - mCenter(0)) + mOrientation(i, 1) * (P3global(1) - mCenter(1)) + mOrientation(i, 2) * (P3global(2) - mCenter(2));
      mP[3](i) = mOrientation(i, 0) * (P4global(0) - mCenter(0)) + mOrientation(i, 1) * (P4global(1) - mCenter(1)) + mOrientation(i, 2) * (P4global(2) - mCenter(2));
    }
  }

  ShellQ4_LocalCoordinateSystem(const Vector3Type & P1global,
                                const Vector3Type & P2global,
                                const Vector3Type & P3global,
                                const Vector3Type & P4global,
                                double alpha)
      : mP(4)
      , mOrientation(3, 3)
  {
    // Form the basis vectors alligning the local X direction
    // with the 1-2 side

    // compute the central point
    noalias( mCenter )  = P1global;
    noalias( mCenter ) += P2global;
    noalias( mCenter ) += P3global;
    noalias( mCenter ) += P4global;
    mCenter *= 0.25;

    // compute the diagonal vectors
    Vector3Type d13( P3global - P1global );
    Vector3Type d24( P4global - P2global );

    // compute the Normal vector at the element center
    // as the cross product of the 2 diagonals.
    // While normalizing the normal vector save its norm to compute the area.
    // Note: the norm should be divided by 2 because it's computed from the
    // cross product of the diagonals, which gives twice the area!
    Vector3Type e3;
    MathUtils<RealType>::CrossProduct(e3,        d13, d24);
    mArea = NormalizeVector3( e3 );
    mArea /= 2.0;

    // compute the local X direction parallel to the projection of the side 1-2 onto
    // the local xy plane.
    Vector3Type e1( P2global - P1global );
    RealType e1_dot_e3 = MathUtils<RealType>::Dot(e1, e3);
    noalias( e1 ) -= e1_dot_e3 * e3;
    Quaternion<RealType>::FromAxisAngle(e3(0), e3(1), e3(2), alpha).RotateVector3(e1);
    NormalizeVector3( e1 );

    // finally compute the local Y direction to be orthogonal to both X and Z local directions
    Vector3Type e2;
    MathUtils<RealType>::CrossProduct(e2,    e3, e1);
    NormalizeVector3( e2 );

    // form the 3x3 transformation matrix
    for(int i = 0; i < 3; i++)
    {
      mOrientation(0, i) = e1(i);
      mOrientation(1, i) = e2(i);
      mOrientation(2, i) = e3(i);
    }

    // transform global coordinates to the local coordinate system
    for(int i = 0; i < 3; i++)
    {
      mP[0](i) = mOrientation(i, 0) * (P1global(0) - mCenter(0)) + mOrientation(i, 1) * (P1global(1) - mCenter(1)) + mOrientation(i, 2) * (P1global(2) - mCenter(2));
      mP[1](i) = mOrientation(i, 0) * (P2global(0) - mCenter(0)) + mOrientation(i, 1) * (P2global(1) - mCenter(1)) + mOrientation(i, 2) * (P2global(2) - mCenter(2));
      mP[2](i) = mOrientation(i, 0) * (P3global(0) - mCenter(0)) + mOrientation(i, 1) * (P3global(1) - mCenter(1)) + mOrientation(i, 2) * (P3global(2) - mCenter(2));
      mP[3](i) = mOrientation(i, 0) * (P4global(0) - mCenter(0)) + mOrientation(i, 1) * (P4global(1) - mCenter(1)) + mOrientation(i, 2) * (P4global(2) - mCenter(2));
    }
  }

 public:

  inline const Vector3ContainerType & Nodes()const { return mP; }

  inline const Vector3Type & P1()const { return mP[0]; }
  inline const Vector3Type & P2()const { return mP[1]; }
  inline const Vector3Type & P3()const { return mP[2]; }
  inline const Vector3Type & P4()const { return mP[3]; }
  inline const Vector3Type & Center()const { return mCenter; }

  inline const RealType X1()const { return mP[0][0]; }
  inline const RealType X2()const { return mP[1][0]; }
  inline const RealType X3()const { return mP[2][0]; }
  inline const RealType X4()const { return mP[3][0]; }

  inline const RealType Y1()const { return mP[0][1]; }
  inline const RealType Y2()const { return mP[1][1]; }
  inline const RealType Y3()const { return mP[2][1]; }
  inline const RealType Y4()const { return mP[3][1]; }

  inline const RealType Z1()const { return mP[0][2]; }
  inline const RealType Z2()const { return mP[1][2]; }
  inline const RealType Z3()const { return mP[2][2]; }
  inline const RealType Z4()const { return mP[3][2]; }

  inline const RealType Area()const { return mArea; }

  inline const MatrixType & Orientation()const { return mOrientation; }

  inline const ConstMatrixRowType Vx()const { return row( mOrientation, 0 ); }
  inline const ConstMatrixRowType Vy()const { return row( mOrientation, 1 ); }
  inline const ConstMatrixRowType Vz()const { return row( mOrientation, 2 ); }

  inline const double WarpageFactor()const { return this->Z1(); }
  inline const bool IsWarped()const { return std::abs(this->WarpageFactor()) > 0.0; }

  inline void ComputeTotalRotationMatrix( MatrixType & R )const
  {
    size_t mat_size = 24;
    if(R.size1() != mat_size || R.size2() != mat_size)
      R.resize(mat_size, mat_size, false);

    noalias( R ) = ZeroMatrix(mat_size, mat_size);

    for(size_t k = 0; k < 8; k++)
    {
      size_t i = k * 3;
      R(i  , i  ) = mOrientation(0, 0);   R(i  , i+1) = mOrientation(0, 1);   R(i  , i+2) = mOrientation(0, 2);
      R(i+1, i  ) = mOrientation(1, 0);   R(i+1, i+1) = mOrientation(1, 1);   R(i+1, i+2) = mOrientation(1, 2);
      R(i+2, i  ) = mOrientation(2, 0);   R(i+2, i+1) = mOrientation(2, 1);   R(i+2, i+2) = mOrientation(2, 2);
    }
  }

  inline void ComputeTotalWarpageMatrix( MatrixType & W, RealType wf )const
  {
    size_t mat_size = 24;
    if(W.size1() != mat_size || W.size2() != mat_size)
      W.resize(mat_size, mat_size, false);

    noalias( W ) = IdentityMatrix(mat_size, mat_size);

    W( 0,  4) = -wf;
    W( 1,  3) =  wf;

    W( 6, 10) =  wf;
    W( 7,  9) = -wf;

    W(12, 16) = -wf;
    W(13, 15) =  wf;

    W(18, 22) =  wf;
    W(19, 21) = -wf;
  }

  inline void ComputeTotalWarpageMatrix( MatrixType & W )const
  {
    ComputeTotalWarpageMatrix( W, this->WarpageFactor() );
  }

  inline void ComputeLocalToGlobalTransformationMatrix( MatrixType & R )
  {
    size_t mat_size = 24;
    if(R.size1() != mat_size || R.size2() != mat_size)
      R.resize(mat_size, mat_size, false);

    noalias( R ) = ZeroMatrix(mat_size, mat_size);

    // form the global transformation matrix
    for(size_t k = 0; k < 8; k++)
    {
      size_t i = k * 3;
      R(i  , i  ) = mOrientation(0, 0);   R(i  , i+1) = mOrientation(0, 1);   R(i  , i+2) = mOrientation(0, 2);
      R(i+1, i  ) = mOrientation(1, 0);   R(i+1, i+1) = mOrientation(1, 1);   R(i+1, i+2) = mOrientation(1, 2);
      R(i+2, i  ) = mOrientation(2, 0);   R(i+2, i+1) = mOrientation(2, 1);   R(i+2, i+2) = mOrientation(2, 2);
    }

    // if needed, include the warpage correction.
    if(this->IsWarped())
    {
      MatrixType W( IdentityMatrix(24, 24) );
      RealType z1 = this->Z1();
      RealType z2 = this->Z2();
      RealType z3 = this->Z3();
      RealType z4 = this->Z4();
      W( 0,  4) = -z1;
      W( 1,  3) =  z1;
      W( 6, 10) = -z2;
      W( 7,  9) =  z2;
      W(12, 16) = -z3;
      W(13, 15) =  z3;
      W(18, 22) = -z4;
      W(19, 21) =  z4;
      R = prod( W, R );
    }
  }

  template<int Tid>
  inline const Vector3Type & P()const {
    KRATOS_ERROR_IF_NOT(Tid > 0 && Tid < 5) << "The index should be 1-based, from 1 to 4" << std::endl;
    return mP[Tid - 1];
  }

  template<int TComponent>
  inline const double Coord_ij(int Tid1, int Tid2)const {
    KRATOS_ERROR_IF_NOT(TComponent >= 0 && TComponent < 3) <<
        "The component index should be 0-based, from 0 to 2" <<std::endl;
    double ci, cj;
    if(Tid1 < 5) {
      int i = Tid1 - 1;
      ci = mP[i][TComponent];
    }
    else  {
      int i = Tid1 - 1 - 4;
      ci = mP[i][TComponent] + mP[i == 3 ? 0 : i + 1][TComponent];
    }
    if(Tid2 < 5) {
      int i = Tid2 - 1;
      cj = mP[i][TComponent];
    }
    else  {
      int i = Tid2 - 1 - 4;
      cj = mP[i][TComponent] + mP[i == 3 ? 0 : i + 1][TComponent];
    }
    return cj - ci;
  }

  inline const double Xij(int i, int j)const { return Coord_ij<0>(i,j); }

  inline const double Yij(int i, int j)const { return Coord_ij<1>(i,j); }

  inline const double Zij(int i, int j)const { return Coord_ij<2>(i,j); }

};

}


#endif // KRATOS_SHELLQ4_LOCAL_COORDINATE_SYSTEM_H_INCLUDED
