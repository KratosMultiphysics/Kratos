//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Andrea Gorgi
//

#if !defined(KRATOS_QUADRATURE_POINT_COUPLING_GEOMETRY_2D_H_INCLUDED )
#define KRATOS_QUADRATURE_POINT_COUPLING_GEOMETRY_2D_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/variables.h"
#include "geometries/quadrature_point_curve_on_surface_geometry.h"

namespace Kratos
{
/**
 * @class QuadraturePointCouplingGeometry2D
 * @ingroup KratosCore
 * @brief A single quadrature point, that can be used for geometries without
 *        a predefined integration scheme, i.e. they can handle material point elements,
 *        isogeometric analysis elements or standard finite elements which are defined
 *        at a single quadrature point.
 *        This point defines a line segment described on an underlying surface.
 *        Shape functions and integration types have to be precomputed and are set from
 *        outside.
 *        The parent pointer can provide the address to the owner of this quadrature point.
 */

template<class TPointType>
class QuadraturePointCouplingGeometry2D : public Geometry<TPointType> {
public:
    KRATOS_CLASS_POINTER_DEFINITION(QuadraturePointCouplingGeometry2D);

    typedef Geometry<TPointType> BaseType;
    typedef Geometry<TPointType> GeometryType;

    typedef typename GeometryType::Pointer GeometryPointer;
    typedef std::vector<GeometryPointer> GeometryPointerVector;

    /// Pointer definition of CouplingGeometry

    typedef TPointType PointType;

    typedef typename BaseType::IndexType IndexType;
    typedef typename BaseType::SizeType SizeType;

    typedef typename BaseType::IntegrationPointsArrayType IntegrationPointsArrayType;
    typedef typename BaseType::PointsArrayType PointsArrayType;
    typedef typename BaseType::CoordinatesArrayType CoordinatesArrayType;
    typedef std::vector<CoordinatesArrayType> CoordinatesArrayVectorType;
    typedef PointerVector<GeometryType> GeometriesArrayType;

    ///@}
    ///@name Public Static Members
    ///@{

    static constexpr IndexType Master = 0;
    static constexpr IndexType Slave = 1;

    static constexpr IndexType CouplingGeometry = 2;

    /// Constructor with points and geometry shape function container
    // QuadraturePointCouplingGeometry2D(
    //     QuadraturePointCurveOnSurfaceGeometryPointer pMasterQuadraturePoint,
    //     QuadraturePointCurveOnSurfaceGeometryPointer pSlaveQuadraturePoint,
    //     const PointsArrayType& ThisPoints,
    //     GeometryShapeFunctionContainerType& ThisGeometryShapeFunctionContainer,
    //     GeometryType* pGeometryParent
    //     ): BaseType(ThisPoints, ThisGeometryShapeFunctionContainer, pGeometryParent)
    //     , mpQuadraturePointMaster(pMasterQuadraturePoint)
    //     , mpQuadraturePointSlave(pSlaveQuadraturePoint)
    // {

    // }

    // QuadraturePointCouplingGeometry2D(
    //     QuadraturePointCurveOnSurfaceGeometryPointer pSlaveQuadraturePoint,
    //     const PointsArrayType& ThisPoints,
    //     GeometryShapeFunctionContainerType& ThisGeometryShapeFunctionContainer,
    //     GeometryType* pGeometryParent
    //     ): BaseType(ThisPoints, ThisGeometryShapeFunctionContainer, pGeometryParent)
    //     , mpQuadraturePointSlave(pSlaveQuadraturePoint)
    // {

    // }

    QuadraturePointCouplingGeometry2D(
        GeometryPointer pMasterQuadraturePoint,
        GeometryPointer pSlaveQuadraturePoint,
        GeometryPointer pCouplingGeometry): 
        BaseType()
    {
        mpGeometries.resize(3);

        mpGeometries[0] = pMasterQuadraturePoint;
        mpGeometries[1] = pSlaveQuadraturePoint;

        mpGeometries[2] = pCouplingGeometry;
    }


    void SetGeometryPart(
        const IndexType Index,
        GeometryPointer pGeometry
        ) override
    {
        KRATOS_DEBUG_ERROR_IF(mpGeometries.size() <= Index) << "Index out of range: "
            << Index << " composite contains only of: "
            << mpGeometries.size() << " geometries." << std::endl;

        if (Index == 0){
            this->SetGeometryData(&(pGeometry->GetGeometryData()));
        }

        mpGeometries[Index] = pGeometry;
    }

    GeometryPointer pGetGeometryPart(const IndexType Index) override
    {
        return mpGeometries[Index];
    }

    const GeometryPointer pGetGeometryPart(const IndexType Index) const override
    {
        return mpGeometries[Index];
    }

    std::string Info() const
    {
        return "Quadrature point for two contact-paired points in curve on surface.";
    }

    /// Print information about this object.
    void PrintInfo( std::ostream& rOStream ) const
    {
        rOStream << "Quadrature point for two contact-paired points in curve on surface.";
    }

    /// Print object's data.
    void PrintData( std::ostream& rOStream ) const
    {
    }

    ///@}

private:

    GeometryPointerVector mpGeometries;


    // typename QuadraturePointCurveOnSurfaceGeometry<TPointType>::Pointer mpQuadraturePointMaster;
    // typename QuadraturePointCurveOnSurfaceGeometry<TPointType>::Pointer mpQuadraturePointSlave;

    ///@}
}; // Class QuadraturePointCouplingGeometry2D

}  // namespace Kratos.

#endif // KRATOS_QUADRATURE_POINT_COUPLING_GEOMETRY_2D_H_INCLUDED
