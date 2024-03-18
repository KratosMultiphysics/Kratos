//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Tobias Teschemacher
//                   Veronika Singer
//                   Philipp Bucher
//

#if !defined(KRATOS_QUADRATURE_POINT_CURVE_ON_SURFACE_GEOMETRY_H_INCLUDED )
#define  KRATOS_QUADRATURE_POINT_CURVE_ON_SURFACE_GEOMETRY_H_INCLUDED

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
 * @brief A sinlge quadrature point, that can be used for geometries without
 *        a predefined integration scheme, i.e. they can handle material point elements,
 *        isogeometric analysis elements or standard finite elements which are defined
 *        at a single quadrature point.
 *        This point defines a line segment described on a underlying surface.
 *        Shape functions and integration types have to be precomputed and are set from
 *        from outside.
 *        The parent pointer can provide the adress to the owner of this quadrature point.
 */
template<class TPointType>
class QuadraturePointCurveOnSurfaceGeometry
    : public QuadraturePointGeometry<TPointType, 3, 2, 1>
{
public:

    /// Pointer definition of QuadraturePointGeometry
    KRATOS_CLASS_POINTER_DEFINITION(QuadraturePointCurveOnSurfaceGeometry);

    typedef QuadraturePointGeometry<TPointType, 3, 2, 1> BaseType;
    typedef Geometry<TPointType> GeometryType;

    typedef typename GeometryType::IndexType IndexType;
    typedef typename GeometryType::SizeType SizeType;

    typedef typename GeometryType::PointsArrayType PointsArrayType;
    typedef typename GeometryType::CoordinatesArrayType CoordinatesArrayType;

    typedef typename GeometryType::IntegrationPointType IntegrationPointType;
    typedef typename GeometryType::IntegrationPointsArrayType IntegrationPointsArrayType;

    typedef typename GeometryData::ShapeFunctionsGradientsType ShapeFunctionsGradientsType;

    typedef GeometryShapeFunctionContainer<GeometryData::IntegrationMethod> GeometryShapeFunctionContainerType;

    typedef typename GeometryType::IntegrationPointsContainerType IntegrationPointsContainerType;
    typedef typename GeometryType::ShapeFunctionsValuesContainerType ShapeFunctionsValuesContainerType;
    typedef typename GeometryType::ShapeFunctionsLocalGradientsContainerType ShapeFunctionsLocalGradientsContainerType;

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
    QuadraturePointCurveOnSurfaceGeometry(
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
    QuadraturePointCurveOnSurfaceGeometry(
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
    ~QuadraturePointCurveOnSurfaceGeometry() override = default;

    /// Copy constructor.
    QuadraturePointCurveOnSurfaceGeometry(QuadraturePointCurveOnSurfaceGeometry const& rOther )
        : BaseType( rOther )
        , mLocalTangentsU(rOther.mLocalTangentsU)
        , mLocalTangentsV(rOther.mLocalTangentsV)
    {
    }

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    QuadraturePointCurveOnSurfaceGeometry& operator=( const QuadraturePointCurveOnSurfaceGeometry& rOther )
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

    /// Calculate with Vector
    void Calculate(
        const Variable<Vector>& rVariable,
        Vector& rOutput) const override
    {
        if (rVariable == DETERMINANTS_OF_JACOBIAN_PARENT)
        {
            DeterminantOfJacobianParent(rOutput);
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
            << "Trying to access Normal of QuadraturePointCurveOnSurface "
            << "with an integration point index != 0." << std::endl;

        Matrix J;
        this->Jacobian(J, IntegrationPointIndex, ThisMethod);

        array_1d<double, 3> a_1 = column(J, 0);
        array_1d<double, 3> a_2 = column(J, 1);

        CoordinatesArrayType normal = a_2 * mLocalTangentsU - a_1 * mLocalTangentsV;

        return normal;
    }

    ///@}
    ///@name Jacobian
    ///@{

    /*
    * @brief returns the respective curve length of this
    *        quadrature point.
    * @param IntegrationPointIndex index should be always 0.
    */
    double DeterminantOfJacobian(
        IndexType IntegrationPointIndex,
        GeometryData::IntegrationMethod ThisMethod) const override
    {
        KRATOS_DEBUG_ERROR_IF(IntegrationPointIndex != 0)
            << "Trying to access DeterminantOfJacobian of QuadraturePointCurveOnSurface "
            << "with an integration point index != 0." << std::endl;

        Matrix J;
        this->Jacobian(J, IntegrationPointIndex, ThisMethod);

        array_1d<double, 3> a_1 = column(J, 0);
        array_1d<double, 3> a_2 = column(J, 1);

        return norm_2(a_1 * mLocalTangentsU + a_2 * mLocalTangentsV);
    }

    /* @brief returns the respective segment length of this
     *        quadrature point. Length of vector always 1.
     * @param rResult vector of results of this quadrature point.
     */
    Vector& DeterminantOfJacobian(
        Vector& rResult,
        GeometryData::IntegrationMethod ThisMethod) const override
    {
        if (rResult.size() != 1)
            rResult.resize(1, false);

        rResult[0] = this->DeterminantOfJacobian(0, ThisMethod);

        return rResult;
    }

    /* @brief returns the respective segment length of this
     *        quadrature point, computed on the parent of this geometry.
     *        Required for reduced quadrature point geometries (Not all
     *        nodes are part of this geometry - used for mapping).
     * @param rResult vector of results of this quadrature point.
     */
    Vector& DeterminantOfJacobianParent(
        Vector& rResult) const
    {
        if (rResult.size() != 1)
            rResult.resize(1, false);

        Matrix J;
        this->GetGeometryParent(0).Jacobian(J, this->IntegrationPoints()[0]);

        rResult[0] = norm_2(column(J, 0) * mLocalTangentsU + column(J, 1) * mLocalTangentsV);

        return rResult;
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
    QuadraturePointCurveOnSurfaceGeometry()
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

#endif // KRATOS_QUADRATURE_POINT_CURVE_ON_SURFACE_GEOMETRY_H_INCLUDED  defined
