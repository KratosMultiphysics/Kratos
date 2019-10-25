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

    /// Constructor with points and all shape function containers separately
    QuadraturePointCurveOnSurfaceGeometry(
        const PointsArrayType& ThisPoints,
        const IntegrationPointsContainerType& rIntegrationPoints,
        ShapeFunctionsValuesContainerType& rShapeFunctionValues,
        ShapeFunctionsLocalGradientsContainerType& rShapeFunctionsDerivativesVector,
        array_1d<double, 2> LocalTangents2d)
        : BaseType(ThisPoints, rIntegrationPoints, rShapeFunctionValues, rShapeFunctionsDerivativesVector)
        , mLocalTangents2d(LocalTangents2d)
    {
    }

    /// Constructor with points and all shape function containers separately including the parent
    QuadraturePointCurveOnSurfaceGeometry(
        const PointsArrayType& ThisPoints,
        const IntegrationPointsContainerType& rIntegrationPoints,
        const ShapeFunctionsValuesContainerType& rShapeFunctionValues,
        const ShapeFunctionsLocalGradientsContainerType& rShapeFunctionsDerivativesVector,
        array_1d<double, 2> LocalTangents2d,
        GeometryType* pGeometryParent)
        : BaseType(ThisPoints, rIntegrationPoints, rShapeFunctionValues, rShapeFunctionsDerivativesVector, pGeometryParent)
        , mLocalTangents2d(LocalTangents2d)
    {
    }

    /// Constructor with points and geometry shape function container
    QuadraturePointCurveOnSurfaceGeometry(
        const PointsArrayType& ThisPoints,
        GeometryShapeFunctionContainerType& ThisGeometryShapeFunctionContainer,
        array_1d<double, 2> LocalTangents2d)
        : BaseType(ThisPoints, ThisGeometryShapeFunctionContainer)
        , mLocalTangents2d(LocalTangents2d)
    {
    }
    /// Constructor with points, geometry shape function container, parent
    QuadraturePointCurveOnSurfaceGeometry(
        const PointsArrayType& ThisPoints,
        GeometryShapeFunctionContainerType& ThisGeometryShapeFunctionContainer,
        array_1d<double, 2> LocalTangents2d,
        GeometryType* pGeometryParent)
        : BaseType(ThisPoints, ThisGeometryShapeFunctionContainer, pGeometryParent)
        , mLocalTangents2d(LocalTangents2d)
    {
    }

    explicit QuadraturePointCurveOnSurfaceGeometry(const PointsArrayType& ThisPoints)
        : BaseType(ThisPoints)
        , mLocalTangents2d(ZeroVector(2))
    {
    }

    /**
     * Copy constructor.
     * Constructs this geometry as a copy of given geometry.
     *
     * @note This copy constructor does not copy the points, thus,
     * the new geometry shares points with the source geometry.
     * Any changes to the new geometry points affect the source
     * geometry points too.
     */
    QuadraturePointCurveOnSurfaceGeometry(QuadraturePointCurveOnSurfaceGeometry const& rOther )
        : BaseType( rOther )
        , mLocalTangents2d(rOther.mLocalTangents2d)
    {
    }

    /**
     * Copy constructor from a geometry with other point type.
     * Construct this geometry as a copy of given geometry which
     * has different type of points. The given goemetry's
     * TOtherPointType* must be implicity convertible to this
     * geometry PointType.
     *
     * @note This copy constructor does not copy the points, thus,
     * the new geometry shares points with the source geometry.
     * Any changes to the new geometry points affect the source
     * geometry points too.
     */
    template<class TOtherPointType>
    QuadraturePointCurveOnSurfaceGeometry(
        QuadraturePointCurveOnSurfaceGeometry<TOtherPointType> const& rOther )
        : BaseType( rOther )
        , mLocalTangents2d(rOther.mLocalTangents2d)
    {
    }

    /// Destructor. Does nothing!!!
    ~QuadraturePointCurveOnSurfaceGeometry() override = default;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    QuadraturePointCurveOnSurfaceGeometry& operator=( const QuadraturePointCurveOnSurfaceGeometry& rOther )
    {
        BaseType::operator=( rOther );

        mLocalTangents2d = rOther.mLocalTangents2d;

        return *this;
    }

    /// Assignment operator for geometries with different point type.
    template<class TOtherPointType>
    QuadraturePointCurveOnSurfaceGeometry& operator=(
        QuadraturePointCurveOnSurfaceGeometry<TOtherPointType> const & rOther )
    {
        BaseType::operator=( rOther );

        mLocalTangents2d = rOther.mLocalTangents2d;

        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    typename GeometryType::Pointer Create( PointsArrayType const& ThisPoints ) const override
    {
        return typename GeometryType::Pointer( new QuadraturePointCurveOnSurfaceGeometry<TPointType>( ThisPoints ) );
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
        IntegrationMethod ThisMethod) const override
    {
        Matrix J;
        this->Jacobian(J, IntegrationPointIndex, ThisMethod);

        array_1d<double, 3> a_1 = column(J, 0);
        array_1d<double, 3> a_2 = column(J, 1);

        CoordinatesArrayType normal = a_2 * mLocalTangents2d[0] - a_1 * mLocalTangents2d[1];

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
        IntegrationMethod ThisMethod) const override
    {
        Matrix J;
        this->Jacobian(J, IntegrationPointIndex, ThisMethod);

        array_1d<double, 3> a_1 = column(J, 0);
        array_1d<double, 3> a_2 = column(J, 1);

        return norm_2(a_1 * mLocalTangents2d[0] + a_2 * mLocalTangents2d[1]);
    }

    ///@}
    ///@name Access to Member Variables
    ///@{

    array_1d<double, 2>& GetLocalTangents2d()
    {
        return mLocalTangents2d;
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

    array_1d<double, 2> mLocalTangents2d;

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save( Serializer& rSerializer ) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseType);
        rSerializer.save("LocalTangents2d", mLocalTangents2d);
    }

    void load( Serializer& rSerializer ) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BaseType );
        rSerializer.load("LocalTangents2d", mLocalTangents2d);
    }

    QuadraturePointCurveOnSurfaceGeometry()
        : BaseType()
        , mLocalTangents2d(
            ZeroVector(2))
    {
    }

    ///@}
}; // Class Geometry

}  // namespace Kratos.

#endif // KRATOS_QUADRATURE_POINT_CURVE_ON_SURFACE_GEOMETRY_H_INCLUDED  defined
