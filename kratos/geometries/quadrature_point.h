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

#if !defined(KRATOS_INTEGRATION_POINT_2D_H_INCLUDED )
#define  KRATOS_INTEGRATION_POINT_2D_H_INCLUDED

// System includes

// External includes

// Project includes
#include "geometries/geometry.h"
#include "geometries/geometry_dimension.h"


namespace Kratos
{
/**
 * @class QuadraturePoint
 * @ingroup KratosCore
 * @brief A sinlge quadrature point, that can be used for geometries without
 *        a predefined integration scheme, i.e. they can handle material point elements,
 *        isogeometric analysis elements or standard finite elements which are defined
 *        at a single quadrature point.
 *        Shape functions and integration types are precomputed and are set from
 *        from outside.
 *        The parent pointer can provide the adress to the owner of this quadrature point.
 */
template<int TWorkingSpaceDimension, int TLocalSpaceDimension, class TPointType>
class QuadraturePoint
    : public Geometry<TPointType>
{
public:

    /// Pointer definition of QuadraturePoint
    KRATOS_CLASS_POINTER_DEFINITION( QuadraturePoint );

    typedef Geometry<TPointType> BaseType;
    typedef Geometry<TPointType> GeometryType;

    typedef typename GeometryType::PointsArrayType PointsArrayType;

    typedef typename GeometryType::IntegrationPointsArrayType IntegrationPointsArrayType;

    typedef typename GeometryData::ShapeFunctionsGradientsType ShapeFunctionsGradientsType;

    typedef GeometryShapeFunctionContainer<IntegrationMethod> GeometryShapeFunctionContainerType;

    typedef typename GeometryType::IntegrationMethod IntegrationMethod;
    typedef typename GeometryType::IntegrationPointsContainerType IntegrationPointsContainerType;
    typedef typename GeometryType::ShapeFunctionsValuesContainerType ShapeFunctionsValuesContainerType;
    typedef typename GeometryType::ShapeFunctionsLocalGradientsContainerType ShapeFunctionsLocalGradientsContainerType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor with points and all shape function containers separately
    QuadraturePoint(
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
    QuadraturePoint(
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
    QuadraturePoint(
        const PointsArrayType& ThisPoints,
        GeometryShapeFunctionContainerType& ThisGeometryShapeFunctionContainer)
        : BaseType(ThisPoints, &mGeometryData)
        , mGeometryData(
            &msGeometryDimension,
            ThisGeometryShapeFunctionContainer)
    {
    }
    /// Constructor with points, geometry shape function container, parent
    QuadraturePoint(
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

    explicit QuadraturePoint(const PointsArrayType& ThisPoints)
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
    QuadraturePoint( QuadraturePoint const& rOther )
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
    QuadraturePoint( QuadraturePoint<TWorkingSpaceDimension, TLocalSpaceDimension, TOtherPointType> const& rOther )
        : BaseType( rOther )
    {
    }

    /// Destructor. Does nothing!!!
    ~QuadraturePoint() override {}

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
    QuadraturePoint& operator=( const QuadraturePoint& rOther )
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
    QuadraturePoint& operator=( QuadraturePoint<TWorkingSpaceDimension, TLocalSpaceDimension, TOtherPointType> const & rOther )
    {
        BaseType::operator=( rOther );

        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    typename BaseType::Pointer Create( PointsArrayType const& ThisPoints ) const override
    {
        return typename BaseType::Pointer( new QuadraturePoint( ThisPoints ) );
    }

    ///@}
    ///@name Parent
    ///@{

    //GeometryType& GetGeometryParent(IndexType Index) const override
    //{
    //    return *mpGeometryParent;
    //}

    //void SetGeometryParent(GeometryType* pGeometryParent) override
    //{
    //    mpGeometryParent = pGeometryParent;
    //}

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
        const std::size_t points_number = this->PointsNumber();

        Point point(0.0, 0.0, 0.0);
        const Matrix& N = this->ShapeFunctionsValues();

        for (IndexType point_number = 0; point_number < this->IntegrationPointsNumber(); ++point_number) {
            for (IndexType i = 0; i < this->PointsNumber(); ++i) {
                point += (*this)[i] * N(point_number, i);
            }
        }
        return point;
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
    }

    void load( Serializer& rSerializer ) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BaseType );
    }

    QuadraturePoint()
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
    ///@name Private Friend
    ///@{

    //template<class TOtherPointType> friend class QuadraturePoint<TWorkingSpaceDimension, TLocalSpaceDimension, TPointType>;

    ///@}
}; // Class Geometry

///@name Input and output
///@{

/// input stream function
template<int TWorkingSpaceDimension, int TLocalSpaceDimension, class TPointType>
inline std::istream& operator >> (
    std::istream& rIStream,
    QuadraturePoint<TWorkingSpaceDimension, TLocalSpaceDimension, TPointType>& rThis );

/// output stream function
template<int TWorkingSpaceDimension, int TLocalSpaceDimension, class TPointType>
inline std::ostream& operator << (
    std::ostream& rOStream,
    const QuadraturePoint<TWorkingSpaceDimension, TLocalSpaceDimension, TPointType>& rThis )
{
    rThis.PrintInfo( rOStream );
    rOStream << std::endl;
    rThis.PrintData( rOStream );

    return rOStream;
}

///@}
///@name Type Dimension Definition
///@{

template<int TWorkingSpaceDimension, int TLocalSpaceDimension, class TPointType>
const GeometryDimension QuadraturePoint<TWorkingSpaceDimension, TLocalSpaceDimension, TPointType>::msGeometryDimension(
    TLocalSpaceDimension,
    TWorkingSpaceDimension,
    TLocalSpaceDimension);

///@}

}  // namespace Kratos.

#endif // KRATOS_INTEGRATION_POINT_2D_H_INCLUDED  defined
