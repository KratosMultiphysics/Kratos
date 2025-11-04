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

#pragma once

// System includes

// External includes

// Project includes
#include "includes/variables.h"
#include "geometries/quadrature_point_geometry.h"

namespace Kratos
{
/**
 * @class QuadraturePointGeometry
 * @ingroup KratosCore
 * @brief A single quadrature point, that can be used for geometries without
 *        a predefined integration scheme, i.e. they can handle material point elements,
 *        isogeometric analysis elements or standard finite elements which are defined
 *        at a single quadrature point.
 *        This point defines a line segment described on a underlying surface.
 *        Shape functions and integration types have to be precomputed and are set from
 *        from outside.
 *        The parent pointer can provide the address to the owner of this quadrature point.
 */
template<class TPointType, int TWorkingSpaceDimension = 3>
class QuadraturePointCurveGeometry
    : public QuadraturePointGeometry<TPointType, TWorkingSpaceDimension, 2, 1>
{
public:

    /// Pointer definition of QuadraturePointGeometry
    KRATOS_CLASS_POINTER_DEFINITION(QuadraturePointCurveGeometry);

    using BaseType = QuadraturePointGeometry<TPointType, TWorkingSpaceDimension, 2, 1>;
    using GeometryType = Geometry<TPointType>;

    using IndexType = typename GeometryType::IndexType;
    using SizeType = typename GeometryType::SizeType;

    using PointsArrayType = typename GeometryType::PointsArrayType;
    using CoordinatesArrayType = typename GeometryType::CoordinatesArrayType;

    using IntegrationPointType = typename GeometryType::IntegrationPointType;
    using IntegrationPointsArrayType = typename GeometryType::IntegrationPointsArrayType;

    using ShapeFunctionsGradientsType = typename GeometryData::ShapeFunctionsGradientsType;

    using GeometryShapeFunctionContainerType = GeometryShapeFunctionContainer<GeometryData::IntegrationMethod>;

    using IntegrationPointsContainerType = typename GeometryType::IntegrationPointsContainerType;
    using ShapeFunctionsValuesContainerType = typename GeometryType::ShapeFunctionsValuesContainerType;
    using ShapeFunctionsLocalGradientsContainerType = typename GeometryType::ShapeFunctionsLocalGradientsContainerType;

    /// using base class functions
    using BaseType::Jacobian;
    using BaseType::DeterminantOfJacobian;
    using BaseType::ShapeFunctionsValues;
    using BaseType::ShapeFunctionsLocalGradients;
    using BaseType::InverseOfJacobian;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor with points and geometry shape function container
    QuadraturePointCurveGeometry(
        const PointsArrayType& ThisPoints,
        GeometryShapeFunctionContainerType& ThisGeometryShapeFunctionContainer,
        double LocalTangentsU,
        double LocalTangentsV)
        : BaseType(ThisPoints, ThisGeometryShapeFunctionContainer)
        , mLocalTangentsU(LocalTangentsU)
        , mLocalTangentsV(LocalTangentsV)
    {
    }

    /// Constructor with points, geometry shape function container, parent
    QuadraturePointCurveGeometry(
        const PointsArrayType& ThisPoints,
        GeometryShapeFunctionContainerType& ThisGeometryShapeFunctionContainer,
        double LocalTangentsU,
        double LocalTangentsV,
        GeometryType* pGeometryParent)
        : BaseType(ThisPoints, ThisGeometryShapeFunctionContainer, pGeometryParent)
        , mLocalTangentsU(LocalTangentsU)
        , mLocalTangentsV(LocalTangentsV)
    {
    }

    /// Destructor.
    ~QuadraturePointCurveGeometry() override = default;

    /// Copy constructor.
    QuadraturePointCurveGeometry(QuadraturePointCurveGeometry const& rOther )
        : BaseType( rOther )
        , mLocalTangentsU(rOther.mLocalTangentsU)
        , mLocalTangentsV(rOther.mLocalTangentsV)
    {
    }

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    QuadraturePointCurveGeometry& operator=( const QuadraturePointCurveGeometry& rOther )
    {
        BaseType::operator=( rOther );

        mLocalTangentsU = rOther.mLocalTangentsU;
        mLocalTangentsV = rOther.mLocalTangentsV;

        return *this;
    }

    ///@}
    ///@name Dynamic access to internals
    ///@{

    /// Assign with array_1d<double, 3>
    void Assign(
        const Variable<array_1d<double, 3>>& rVariable,
        const array_1d<double, 3>& rInput) override
    {
        if (rVariable == LOCAL_TANGENT)
        {
            mLocalTangentsU = rInput[0];
            mLocalTangentsV = rInput[1];
        }
    }

    /// Calculate with array_1d<double, 3>
    void Calculate(
        const Variable<array_1d<double, 3>>& rVariable,
        array_1d<double, 3>& rOutput) const override
    {
        if (rVariable == LOCAL_TANGENT)
        {
            rOutput[0] = mLocalTangentsU;
            rOutput[1] = mLocalTangentsV;
            rOutput[2] = 0.0;
        }
    }

    ///@}
    ///@name Normal
    ///@{

    /*
    * @brief computes the normal of the curve
    *        laying on the underlying surface.
    * @param rResult Normal to the curve lying on the surface.
    * @param IntegrationPointIndex index should be always 0.
    */
    CoordinatesArrayType Normal(
        IndexType IntegrationPointIndex,
        GeometryData::IntegrationMethod ThisMethod) const override
    {
        KRATOS_DEBUG_ERROR_IF(IntegrationPointIndex != 0)
            << "Trying to access Normal of QuadraturePointCurveGeometryOnSurface "
            << "with an integration point index != 0." << std::endl;

        CoordinatesArrayType normal = ZeroVector(3);
        normal[0] = mLocalTangentsV;
        normal[1] = -mLocalTangentsU;

        return normal;
    }

    ///@}
    ///@name Kratos Geometry Families
    ///@{

    /**
     * @brief Gets the geometry family.
     * @details This function returns the family type of the geometry. The geometry family categorizes the geometry into a broader classification, aiding in its identification and processing.
     * @return GeometryData::KratosGeometryFamily The geometry family.
     */
    GeometryData::KratosGeometryFamily GetGeometryFamily() const override
    {
        return GeometryData::KratosGeometryFamily::Kratos_Quadrature_Geometry;
    }

    /**
     * @brief Gets the geometry type.
     * @details This function returns the specific type of the geometry. The geometry type provides a more detailed classification of the geometry.
     * @return GeometryData::KratosGeometryType The specific geometry type.
     */
    GeometryData::KratosGeometryType GetGeometryType() const override
    {
        return GeometryData::KratosGeometryType::Kratos_Quadrature_Point_Curve_On_Surface_Geometry;
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "Quadrature point for a curve on surface.";
    }

    /// Print information about this object.
    void PrintInfo( std::ostream& rOStream ) const override
    {
        rOStream << "Quadrature point for a curve on surface.";
    }

    /// Print object's data.
    void PrintData( std::ostream& rOStream ) const override
    {
    }
    ///@}

private:
    ///@name Static Member Variables
    ///@{

    double mLocalTangentsU;
    double mLocalTangentsV;

    ///@}
    ///@name Serialization
    ///@{

    /// Default constructor for serializer
    QuadraturePointCurveGeometry()
        : BaseType()
    {
    }

    friend class Serializer;

    void save( Serializer& rSerializer ) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseType);
        rSerializer.save("LocalTangentsU", mLocalTangentsU);
        rSerializer.save("LocalTangentsV", mLocalTangentsV);
    }

    void load( Serializer& rSerializer ) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BaseType );
        rSerializer.load("LocalTangentsU", mLocalTangentsU);
        rSerializer.load("LocalTangentsV", mLocalTangentsV);
    }

    ///@}
}; // Class Geometry

}  // namespace Kratos.
