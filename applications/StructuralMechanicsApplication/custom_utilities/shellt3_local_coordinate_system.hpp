// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Massimo Petracca
//
#if !defined(KRATOS_SHELLT3_LOCAL_COORDINATE_SYSTEM_H_INCLUDED )
#define  KRATOS_SHELLT3_LOCAL_COORDINATE_SYSTEM_H_INCLUDED

#include "utilities/quaternion.h"

namespace Kratos
{

/** \brief ShellT3_LocalCoordinateSystem
*
* This class represent the local coordinate system of any element whose geometry
* is a TRIANGLE 3 in 3D space
*/
class ShellT3_LocalCoordinateSystem
{

public:

    KRATOS_CLASS_POINTER_DEFINITION( ShellT3_LocalCoordinateSystem );

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
        if(vn != 0.0 && vn != 1.0)
        {
            vn = std::sqrt(vn);
            V /= vn;
        }
        return vn;
    }

public:

    ShellT3_LocalCoordinateSystem(const Vector3Type & P1global,
                                  const Vector3Type & P2global,
                                  const Vector3Type & P3global)
        : mP(3)
        , mOrientation(3, 3)
    {
        // Form the basis vectors alligning the local X direction
        // with the 1-2 side

        // compute the central point
        noalias( mCenter )  = P1global;
        noalias( mCenter ) += P2global;
        noalias( mCenter ) += P3global;
        mCenter /= 3.0;

        // compute the local X direction parallel to the side 1-2
        Vector3Type e1( P2global - P1global );

        // compute the temporary local Y direction
        Vector3Type e2( P3global - P1global );

        // compute the Normal vector at the element center
        // While normalizing the normal vector save its norm to compute the area.
        Vector3Type e3;
        MathUtils<RealType>::CrossProduct(e3,    e1, e2);
        mArea = NormalizeVector3( e3 );
        mArea /= 2.0;

        // finally compute the local Y direction to be orthogonal to both X and Z local directions
        MathUtils<RealType>::CrossProduct(e2,    e3, e1);

        // normalize the X and Y directions
        NormalizeVector3( e1 );
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
        }
    }

    ShellT3_LocalCoordinateSystem(const Vector3Type & P1global,
                                  const Vector3Type & P2global,
                                  const Vector3Type & P3global,
                                  double alpha)
        : mP(3)
        , mOrientation(3, 3)
    {
        // Form the basis vectors alligning the local X direction
        // with the 1-2 side

        // compute the central point
        noalias( mCenter )  = P1global;
        noalias( mCenter ) += P2global;
        noalias( mCenter ) += P3global;
        mCenter /= 3.0;

        // compute the local X direction parallel to the side 1-2
        Vector3Type e1( P2global - P1global );

        // compute the temporary local Y direction
        Vector3Type e2( P3global - P1global );

        // compute the Normal vector at the element center
        // While normalizing the normal vector save its norm to compute the area.
        Vector3Type e3;
        MathUtils<RealType>::CrossProduct(e3,    e1, e2);
        mArea = NormalizeVector3( e3 );
        mArea /= 2.0;

        // apply the extra in-plane rotation
        Quaternion<RealType>::FromAxisAngle(e3(0), e3(1), e3(2), alpha).RotateVector3(e1);

        // finally compute the local Y direction to be orthogonal to both X and Z local directions
        MathUtils<RealType>::CrossProduct(e2,    e3, e1);

        // normalize the X and Y directions
        NormalizeVector3( e1 );
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
        }
    }

public:

    inline const Vector3ContainerType & Nodes()const
    {
        return mP;
    }

    inline const Vector3Type & P1()const
    {
        return mP[0];
    }
    inline const Vector3Type & P2()const
    {
        return mP[1];
    }
    inline const Vector3Type & P3()const
    {
        return mP[2];
    }
    inline const Vector3Type & Center()const
    {
        return mCenter;
    }

    inline RealType X1()const
    {
        return mP[0][0];
    }
    inline RealType X2()const
    {
        return mP[1][0];
    }
    inline RealType X3()const
    {
        return mP[2][0];
    }

    inline RealType Y1()const
    {
        return mP[0][1];
    }
    inline RealType Y2()const
    {
        return mP[1][1];
    }
    inline RealType Y3()const
    {
        return mP[2][1];
    }

    inline RealType Z1()const
    {
        return mP[0][2];
    }
    inline RealType Z2()const
    {
        return mP[1][2];
    }
    inline RealType Z3()const
    {
        return mP[2][2];
    }

    inline RealType Area()const
    {
        return mArea;
    }

    inline const MatrixType & Orientation()const
    {
        return mOrientation;
    }

    inline const ConstMatrixRowType Vx()const
    {
        return row( mOrientation, 0 );
    }
    inline const ConstMatrixRowType Vy()const
    {
        return row( mOrientation, 1 );
    }
    inline const ConstMatrixRowType Vz()const
    {
        return row( mOrientation, 2 );
    }

    inline void ComputeTotalRotationMatrix( MatrixType & R )const
    {
        size_t mat_size = 18;
        if(R.size1() != mat_size || R.size2() != mat_size)
            R.resize(mat_size, mat_size, false);

        noalias( R ) = ZeroMatrix(mat_size, mat_size);

        for(size_t k = 0; k < 6; k++)
        {
            size_t i = k * 3;
            R(i  , i  ) = mOrientation(0, 0);
            R(i  , i+1) = mOrientation(0, 1);
            R(i  , i+2) = mOrientation(0, 2);
            R(i+1, i  ) = mOrientation(1, 0);
            R(i+1, i+1) = mOrientation(1, 1);
            R(i+1, i+2) = mOrientation(1, 2);
            R(i+2, i  ) = mOrientation(2, 0);
            R(i+2, i+1) = mOrientation(2, 1);
            R(i+2, i+2) = mOrientation(2, 2);
        }
    }

};

}


#endif // KRATOS_SHELLT3_LOCAL_COORDINATE_SYSTEM_H_INCLUDED
