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
//

#if !defined(KRATOS_SHAPE_SHELLQ4_COORDINATE_TRANSFORMATION_H_INCLUDED )
#define  KRATOS_SHAPE_SHELLQ4_COORDINATE_TRANSFORMATION_H_INCLUDED

#include "shape_shellq4_local_coordinate_system.hpp"

namespace Kratos
{
/** \brief ShapeShellQ4_CoordinateTransformation
*
* This class represents a basic (linear) coordinate transformation that can be used
* by any element whose geometry is a QUAD 4 in 3D space, with 6 D.O.F.s per node.
* It's main aim is to:
* 1) Create the local coordinate system
* 2) Transform the incoming global displacements in local coordinate system
* 3) Transform the outgoing matrices and vectors in global coordinate system
*/
class ShapeShellQ4_CoordinateTransformation
{

public:

    KRATOS_CLASS_POINTER_DEFINITION(ShapeShellQ4_CoordinateTransformation);

    typedef Element::GeometryType GeometryType;

    typedef Vector VectorType;

    typedef Matrix MatrixType;

public:

    ShapeShellQ4_CoordinateTransformation(const GeometryType::Pointer& pGeometry)
        : mpGeometry(pGeometry)
    {
    }

    virtual ~ShapeShellQ4_CoordinateTransformation()
    {
    }

private:

    ShapeShellQ4_CoordinateTransformation(const ShapeShellQ4_CoordinateTransformation& other);

    ShapeShellQ4_CoordinateTransformation& operator = (const ShapeShellQ4_CoordinateTransformation& other);

public:

    virtual ShapeShellQ4_CoordinateTransformation::Pointer Create(GeometryType::Pointer pGeometry)const
    {
        return ShapeShellQ4_CoordinateTransformation::Pointer(new ShapeShellQ4_CoordinateTransformation(pGeometry));
    }

    virtual void Initialize()
    {
    }

    virtual void InitializeSolutionStep()
    {
    }

    virtual void FinalizeSolutionStep()
    {
    }

    virtual void InitializeNonLinearIteration()
    {
    }

    virtual void FinalizeNonLinearIteration()
    {
    }

    virtual ShapeShellQ4_LocalCoordinateSystem CreateReferenceCoordinateSystem()const
    {
        const GeometryType& geom = GetGeometry();
        return ShapeShellQ4_LocalCoordinateSystem(geom[0].GetInitialPosition(),
                                             geom[1].GetInitialPosition(),
                                             geom[2].GetInitialPosition(),
                                             geom[3].GetInitialPosition());
    }

    virtual ShapeShellQ4_LocalCoordinateSystem CreateLocalCoordinateSystem()const
    {
        return CreateReferenceCoordinateSystem();
    }

    virtual Vector CalculateLocalDisplacements(const ShapeShellQ4_LocalCoordinateSystem& LCS,
            const VectorType& globalDisplacements)
    {
        MatrixType R(24, 24);
        LCS.ComputeTotalRotationMatrix(R);
        if (LCS.IsWarped()) {
            MatrixType W(24, 24);
            LCS.ComputeTotalWarpageMatrix(W);
            R = prod(W, R);
        }
        return prod(R, globalDisplacements);
    }

    virtual void FinalizeCalculations(const ShapeShellQ4_LocalCoordinateSystem& LCS,
                                      const VectorType& globalDisplacements,
                                      const VectorType& localDisplacements,
                                      MatrixType& rLeftHandSideMatrix,
                                      VectorType& rRightHandSideVector,
                                      const bool RHSrequired,
                                      const bool LHSrequired)
    {
        MatrixType R(24, 24);
        LCS.ComputeTotalRotationMatrix(R);
        if (LCS.IsWarped()) {
            MatrixType W(24, 24);
            LCS.ComputeTotalWarpageMatrix(W);
            R = prod(W, R);
        }

        if (LHSrequired) {
            MatrixType temp(24, 24);
            noalias(temp) = prod(trans(R), rLeftHandSideMatrix);
            noalias(rLeftHandSideMatrix) = prod(temp, R);
        }

        if (RHSrequired) {
            rRightHandSideVector = prod(trans(R), rRightHandSideVector);
        }
    }

    virtual MatrixType GetNodalDeformationalRotationTensor(const ShapeShellQ4_LocalCoordinateSystem& LCS,
            const Vector& globalDisplacements,
            size_t nodeid)
    {
        return IdentityMatrix(3);
    }

    virtual MatrixType GetNodalDeformationalRotationTensor(const ShapeShellQ4_LocalCoordinateSystem& LCS,
            const Vector& globalDisplacements,
            const Vector& N)
    {
        return IdentityMatrix(3);
    }

public:

    inline const GeometryType& GetGeometry()const
    {
        return *mpGeometry;
    }

protected:

    ShapeShellQ4_CoordinateTransformation()
        : mpGeometry(GeometryType::Pointer())
    {
    }

private:

    GeometryType::Pointer mpGeometry;

private:

    friend class Serializer;

    virtual void save(Serializer& rSerializer) const
    {
        rSerializer.save("pGeom", mpGeometry);
    }

    virtual void load(Serializer& rSerializer)
    {
        rSerializer.load("pGeom", mpGeometry);
    }

};

}


#endif // KRATOS_SHAPE_SHELLQ4_COORDINATE_TRANSFORMATION_H_INCLUDED
