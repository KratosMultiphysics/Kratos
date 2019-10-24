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
//

#if !defined(KRATOS_QUADRATURE_POINT_GEOMETRY_H_INCLUDED )
#define  KRATOS_QUADRATURE_POINT_GEOMETRY_H_INCLUDED

// System includes

// External includes

// Project includes
#include "geometries/geometry.h"
#include "geometries/geometry_dimension.h"


namespace Kratos
{
/**
 * @class QuadraturePointGeometry
 * @ingroup KratosCore
 * @brief A sinlge quadrature point, that can be used for geometries without
 *        a predefined integration scheme, i.e. they can handle material point elements,
 *        isogeometric analysis elements or standard finite elements which are defined
 *        at a single quadrature point.
 *        Shape functions and integration types are precomputed and are set from
 *        from outside.
 *        The parent pointer can provide the adress to the owner of this quadrature point.
 */
template<class TPointType,
    int TWorkingSpaceDimension,
    int TLocalSpaceDimension = TWorkingSpaceDimension,
    int TDimension = TLocalSpaceDimension>
class QuadraturePointGeometry
    : public Geometry<TPointType>
{
public:

    /// Pointer definition of QuadraturePointGeometry
    KRATOS_CLASS_POINTER_DEFINITION( QuadraturePointGeometry );

    typedef Geometry<TPointType> BaseType;
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
    QuadraturePointGeometry(
        const PointsArrayType& ThisPoints,
        const IntegrationPointsContainerType& rIntegrationPoints,
        ShapeFunctionsValuesContainerType& rShapeFunctionValues,
        ShapeFunctionsLocalGradientsContainerType& rShapeFunctionsDerivativesVector)
        : BaseType(ThisPoints, &mGeometryData)
        , mGeometryData(
            &msGeometryDimension,
            GeometryData::GI_GAUSS_1,
            rIntegrationPoints,
            rShapeFunctionValues,
            rShapeFunctionsDerivativesVector)
    {
    }

    /// Constructor with points and all shape function containers separately including the parent
    QuadraturePointGeometry(
        const PointsArrayType& ThisPoints,
        const IntegrationPointsContainerType& rIntegrationPoints,
        const ShapeFunctionsValuesContainerType& rShapeFunctionValues,
        const ShapeFunctionsLocalGradientsContainerType& rShapeFunctionsDerivativesVector,
        GeometryType* pGeometryParent)
        : BaseType(ThisPoints, &mGeometryData)
        , mGeometryData(
            &msGeometryDimension,
            GeometryData::GI_GAUSS_1,
            rIntegrationPoints,
            rShapeFunctionValues,
            rShapeFunctionsDerivativesVector)
        , mpGeometryParent(pGeometryParent)
    {
    }

    /// Constructor with points and geometry shape function container
    QuadraturePointGeometry(
        const PointsArrayType& ThisPoints,
        GeometryShapeFunctionContainerType& ThisGeometryShapeFunctionContainer)
        : BaseType(ThisPoints, &mGeometryData)
        , mGeometryData(
            &msGeometryDimension,
            ThisGeometryShapeFunctionContainer)
    {
    }
    /// Constructor with points, geometry shape function container, parent
    QuadraturePointGeometry(
        const PointsArrayType& ThisPoints,
        GeometryShapeFunctionContainerType& ThisGeometryShapeFunctionContainer,
        GeometryType* pGeometryParent)
        : BaseType(ThisPoints, &mGeometryData)
        , mGeometryData(
            &msGeometryDimension,
            ThisGeometryShapeFunctionContainer)
        , mpGeometryParent(pGeometryParent)
    {
    }

    explicit QuadraturePointGeometry(const PointsArrayType& ThisPoints)
        : BaseType(ThisPoints, &mGeometryData)
        , mGeometryData(
            &msGeometryDimension,
            GeometryData::GI_GAUSS_1,
            {}, {}, {})
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
    QuadraturePointGeometry( QuadraturePointGeometry const& rOther )
        : BaseType( rOther )
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
    QuadraturePointGeometry( QuadraturePointGeometry<TOtherPointType, TWorkingSpaceDimension, TLocalSpaceDimension, TDimension> const& rOther )
        : BaseType( rOther )
    {
    }

    /// Destructor. Does nothing!!!
    ~QuadraturePointGeometry() override {}

    ///@}
    ///@name Operators
    ///@{

    /**
     * Assignment operator.
     *
     * @note This copy constructor does not copy the points, thus,
     * the new geometry shares points with the source geometry.
     * Any changes to the new geometry points affect the source
     * geometry points too.
     *
     * @see Clone
     * @see ClonePoints
     */
    QuadraturePointGeometry& operator=( const QuadraturePointGeometry& rOther )
    {
        BaseType::operator=( rOther );

        return *this;
    }

    /**
     * Assignment operator for geometries with different point type.
     *
     * @note This copy constructor does not copy the points, thus,
     * the new geometry shares points with the source geometry.
     * Any changes to the new geometry points affect the source
     * geometry points too.
     *
     * @see Clone
     * @see ClonePoints
     */
    template<class TOtherPointType>
    QuadraturePointGeometry& operator=( QuadraturePointGeometry<TOtherPointType, TWorkingSpaceDimension, TLocalSpaceDimension, TDimension> const & rOther )
    {
        BaseType::operator=( rOther );

        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    typename BaseType::Pointer Create( PointsArrayType const& ThisPoints ) const override
    {
        return typename BaseType::Pointer( new QuadraturePointGeometry( ThisPoints ) );
    }

    ///@}
    ///@name Parent
    ///@{

    GeometryType& GetGeometryParent(IndexType Index) const override
    {
        return *mpGeometryParent;
    }

    void SetGeometryParent(GeometryType* pGeometryParent) override
    {
        mpGeometryParent = pGeometryParent;
    }

    ///@}
    ///@name Geometrical Operations
    ///@{

    /// Returns the domain size of this quadrature point.
    double DomainSize() const override
    {
        Vector temp;
        temp = this->DeterminantOfJacobian(temp, GeometryData::GI_GAUSS_1);
        const IntegrationPointsArrayType& r_integration_points = this->IntegrationPoints();
        double domain_size = 0.0;

        for (std::size_t i = 0; i < r_integration_points.size(); ++i) {
            domain_size += temp[i] * r_integration_points[i].Weight();
        }
        return domain_size;
    }

    /**
    * @brief Calculates global location of this integration point.
    * \f[
    * c_i = \sum_j^N(x_j)*x_i
    * \f]
    * j is the index of the node and i the global direction (x,y,z).
    *
    * @return Point which is the location of this quadrature point.
    */
    Point Center() const override
    {
        const std::size_t node_number = this->PointsNumber();

        Point point(0.0, 0.0, 0.0);
        const Matrix& r_N = this->ShapeFunctionsValues();

        for (IndexType point_number = 0; point_number < this->IntegrationPointsNumber(); ++point_number) {
            for (IndexType i = 0; i < node_number; ++i) {
                point += (*this)[i] * r_N(point_number, i);
            }
        }
        return point;
    }

    ///@}
    ///@name Coordinates
    ///@{

    /**
    * This method provides the global coordinates corresponding to the
    * local coordinates. The mapping is done on the parent.
    * Error if Parent is not assigned is only thrown in debug mode.
    */
    CoordinatesArrayType& GlobalCoordinates(
        CoordinatesArrayType& rResult,
        CoordinatesArrayType const& LocalCoordinates
    ) const override
    {
        KRATOS_DEBUG_ERROR_IF(mpGeometryParent == nullptr)
            << "Trying to call GlobalCoordinates(LocalCoordinates) from quadrature point. "
            << "Pointer to parent is not assigned." << std::endl;

        return mpGeometryParent->GlobalCoordinates(rResult, LocalCoordinates);
    }

    CoordinatesArrayType& GlobalCoordinates(
        CoordinatesArrayType& rResult,
        CoordinatesArrayType const& LocalCoordinates,
        Matrix& DeltaPosition
    ) const override
    {
        KRATOS_DEBUG_ERROR_IF(mpGeometryParent == nullptr)
            << "Trying to call GlobalCoordinates(LocalCoordinates, DeltaPosition) from quadrature point. "
            << "Pointer to parent is not assigned." << std::endl;

        return mpGeometryParent->GlobalCoordinates(rResult, LocalCoordinates, DeltaPosition);
    }

    ///@}
    ///@name Jacobian
    ///@{

    /**
    * @brief Jacobian in given point. Computed on parent
    *        geometry.
    * Error if Parent is not assigned is only thrown in debug mode.
    * @param rCoordinates point which jacobians has to
    *        be calculated in it.
    * @return jacobian matrix \f$ J \f$ in given point.
    */
    Matrix& Jacobian(
        Matrix& rResult,
        const CoordinatesArrayType& rCoordinates
    ) const override
    {
        KRATOS_DEBUG_ERROR_IF(mpGeometryParent == nullptr)
            << "Trying to call Jacobian(LocalCoordinates) from quadrature point. "
            << "Pointer to parent is not assigned." << std::endl;

        return mpGeometryParent->Jacobian(rResult, rCoordinates);
    }

    /** Determinant of jacobian in given point. Computed on parent
    *        geometry.
    * Error if Parent is not assigned is only thrown in debug mode.
    * @param rPoint point which determinant of jacobians has to
    *        be calculated in it.
    * @return Determinamt of jacobian matrix \f$ |J| \f$ in given
    *         point.
    */
    double DeterminantOfJacobian(
        const CoordinatesArrayType& rPoint
    ) const override
    {
        KRATOS_DEBUG_ERROR_IF(mpGeometryParent == nullptr)
            << "Trying to call DeterminantOfJacobian(rPoint) from quadrature point. "
            << "Pointer to parent is not assigned." << std::endl;

        return mpGeometryParent->DeterminantOfJacobian(rPoint);
    }

    Matrix& InverseOfJacobian(
        Matrix& rResult,
        const CoordinatesArrayType& rCoordinates
    ) const override
    {
        KRATOS_DEBUG_ERROR_IF(mpGeometryParent == nullptr)
            << "Trying to call InverseOfJacobian(rPoint) from quadrature point. "
            << "Pointer to parent is not assigned." << std::endl;

        return mpGeometryParent->InverseOfJacobian(rResult, rCoordinates);
    }

    ///@}
    ///@name Shape Function
    ///@{

    Vector& ShapeFunctionsValues(
        Vector &rResult,
        const CoordinatesArrayType& rCoordinates
    ) const override
    {
        KRATOS_DEBUG_ERROR_IF(mpGeometryParent == nullptr)
            << "Trying to call ShapeFunctionsValues(rCoordinates) from quadrature point. "
            << "Pointer to parent is not assigned." << std::endl;

        return mpGeometryParent->ShapeFunctionsValues(rResult, rCoordinates);
    }

    virtual Matrix& ShapeFunctionsLocalGradients(
        Matrix& rResult,
        const CoordinatesArrayType& rPoint
    ) const override
    {
        KRATOS_DEBUG_ERROR_IF(mpGeometryParent == nullptr)
            << "Trying to call ShapeFunctionsLocalGradients(rPoint) from quadrature point. "
            << "Pointer to parent is not assigned." << std::endl;

        return mpGeometryParent->ShapeFunctionsLocalGradients(rResult, rPoint);
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "Quadrature point templated by local space dimension and working space dimension.";
    }

    /// Print information about this object.
    void PrintInfo( std::ostream& rOStream ) const override
    {
        rOStream << "Quadrature point templated by local space dimension and working space dimension.";
    }

    /// Print object's data.
    void PrintData( std::ostream& rOStream ) const override
    {
    }
    ///@}

private:
    ///@name Static Member Variables
    ///@{

    static const GeometryDimension msGeometryDimension;

    ///@}
    ///@name Member Variables
    ///@{

    GeometryData mGeometryData;

    // quatrature point can be related to a parent geometry. To keep the connection,
    // this geometry is related to the integration point.
    GeometryType* mpGeometryParent;

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save( Serializer& rSerializer ) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BaseType );
        rSerializer.save("pGeometryParent", mpGeometryParent);
        rSerializer.save("GeometryData", mGeometryData);
    }

    void load( Serializer& rSerializer ) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BaseType );
        rSerializer.load("pGeometryParent", mpGeometryParent);

        GeometryDimension local(3, 3, 3);
        GeometryDimension *temp = &local;
        rSerializer.load("pGeometryDimension", temp);

        rSerializer.load("GeometryData", mGeometryData);

        mGeometryData.SetGeometryDimension(&msGeometryDimension);

    }

    QuadraturePointGeometry()
        : BaseType(
            PointsArrayType(),
            &mGeometryData)
        , mGeometryData(
            &msGeometryDimension,
            GeometryData::GI_GAUSS_1,
            {}, {}, {})
    {
    }

    ///@}
}; // Class Geometry

///@name Input and output
///@{

/// input stream function
template<class TPointType,
    int TWorkingSpaceDimension,
    int TLocalSpaceDimension,
    int TDimension>
inline std::istream& operator >> (
    std::istream& rIStream,
    QuadraturePointGeometry<TPointType, TWorkingSpaceDimension, TLocalSpaceDimension, TDimension>& rThis );

/// output stream function
template<class TPointType,
    int TWorkingSpaceDimension,
    int TLocalSpaceDimension,
    int TDimension>
inline std::ostream& operator << (
    std::ostream& rOStream,
    const QuadraturePointGeometry<TPointType, TWorkingSpaceDimension, TLocalSpaceDimension, TDimension>& rThis )
{
    rThis.PrintInfo( rOStream );
    rOStream << std::endl;
    rThis.PrintData( rOStream );

    return rOStream;
}

///@}
///@name Type Dimension Definition
///@{

template<class TPointType,
    int TWorkingSpaceDimension,
    int TLocalSpaceDimension,
    int TDimension>
const GeometryDimension QuadraturePointGeometry<
    TPointType,
    TWorkingSpaceDimension,
    TLocalSpaceDimension,
    TDimension>::msGeometryDimension(
        TLocalSpaceDimension,
        TWorkingSpaceDimension,
        TLocalSpaceDimension);

///@}

}  // namespace Kratos.

#endif // KRATOS_QUADRATURE_POINT_GEOMETRY_H_INCLUDED  defined
