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

#if !defined(KRATOS_SHAPE_SHELLT3_COORDINATE_TRANSFORMATION_H_INCLUDED )
#define  KRATOS_SHAPE_SHELLT3_COORDINATE_TRANSFORMATION_H_INCLUDED

#include "shape_shellt3_local_coordinate_system.hpp"

namespace Kratos
{
/** \brief ShapeShellT3_CoordinateTransformation
*
* This class represents a basic (linear) coordinate transformation that can be used
* by any element whose geometry is a TRIANGLE 3 in 3D space, with 6 D.O.F.s per node.
* It's main aim is to:
* 1) Create the local coordinate system
* 2) Transform the incoming global displacements in local coordinate system
* 3) Transform the outgoing matrices and vectors in global coordinate system
*/
class ShapeShellT3_CoordinateTransformation
{

public:

    KRATOS_CLASS_POINTER_DEFINITION(ShapeShellT3_CoordinateTransformation);

    typedef Element::GeometryType GeometryType;

    typedef Vector VectorType;

    typedef Matrix MatrixType;

public:

    ShapeShellT3_CoordinateTransformation(const GeometryType::Pointer& pGeometry)
        : mpGeometry(pGeometry)
    {
    }

    virtual ~ShapeShellT3_CoordinateTransformation()
    {
    }

private:

    ShapeShellT3_CoordinateTransformation(const ShapeShellT3_CoordinateTransformation& other);

    ShapeShellT3_CoordinateTransformation& operator = (const ShapeShellT3_CoordinateTransformation& other);

public:

    virtual ShapeShellT3_CoordinateTransformation::Pointer Create(GeometryType::Pointer pGeometry)const
    {
        return ShapeShellT3_CoordinateTransformation::Pointer(new ShapeShellT3_CoordinateTransformation(pGeometry));
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

    virtual ShapeShellT3_LocalCoordinateSystem CreateReferenceCoordinateSystem()const
    {
        const GeometryType& geom = GetGeometry();
        return ShapeShellT3_LocalCoordinateSystem(geom[0].GetInitialPosition(),
                                             geom[1].GetInitialPosition(),
                                             geom[2].GetInitialPosition());
    }

    virtual ShapeShellT3_LocalCoordinateSystem CreateLocalCoordinateSystem()const
    {
        return CreateReferenceCoordinateSystem();
    }

    virtual Vector CalculateLocalDisplacements(const ShapeShellT3_LocalCoordinateSystem& LCS,
            const VectorType& globalDisplacements)
    {
        MatrixType R(18, 18);
        LCS.ComputeTotalRotationMatrix(R);
        return prod(R, globalDisplacements);
    }

    virtual void FinalizeCalculations(const ShapeShellT3_LocalCoordinateSystem& LCS,
                                      const VectorType& globalDisplacements,
                                      const VectorType& localDisplacements,
                                      MatrixType& rLeftHandSideMatrix,
                                      VectorType& rRightHandSideVector,
                                      const bool RHSrequired,
                                      const bool LHSrequired)
    {
        MatrixType R(18, 18);
        LCS.ComputeTotalRotationMatrix(R);

        if (LHSrequired) {
            MatrixType temp(18, 18);
            noalias(temp) = prod(trans(R), rLeftHandSideMatrix);
            noalias(rLeftHandSideMatrix) = prod(temp, R);
        }

        if (RHSrequired) {
            rRightHandSideVector = prod(trans(R), rRightHandSideVector);
        }
    }

    virtual MatrixType GetNodalDeformationalRotationTensor(const ShapeShellT3_LocalCoordinateSystem& LCS,
            const Vector& globalDisplacements,
            size_t nodeid)
    {
        return IdentityMatrix(3);
    }

    virtual MatrixType GetNodalDeformationalRotationTensor(const ShapeShellT3_LocalCoordinateSystem& LCS,
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

    ShapeShellT3_CoordinateTransformation()
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


#endif // KRATOS_SHAPE_SHELLT3_COORDINATE_TRANSFORMATION_H_INCLUDED
