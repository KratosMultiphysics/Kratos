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

#if !defined(KRATOS_INTEGRATION_POINT_SURFACE_3D_H_INCLUDED )
#define  KRATOS_INTEGRATION_POINT_SURFACE_3D_H_INCLUDED

// System includes

// External includes

// Project includes
#include "geometries/geometry.h"
#include "geometries/geometry_dimension.h"


namespace Kratos
{
/**
 * @class IntegrationPointSurface3d
 * @ingroup KratosCore
 * @brief A sinlge integration point, that can be used for geometries without
 *        a predefined integration scheme, i.e. they can handle material point elements,
 *        isogeometric analysis elements or standard finite elements which are defined
 *        at a single integration point.
 *        Shape functions and integration types have to be precomputed and are set from
 *        from outside.
 */
template<class TPointType> class IntegrationPointSurface3d
    : public Geometry<TPointType>
{
public:

    /**
     * Pointer definition of IntegrationPointSurface3d
     */
    KRATOS_CLASS_POINTER_DEFINITION( IntegrationPointSurface3d );

    typedef Geometry<TPointType> BaseType;
    typedef Geometry<TPointType> GeometryType;

    typedef PointerVector<TPointType> PointsArrayType;

    typedef IntegrationPoint<2> IntegrationPointType;
    typedef typename GeometryType::IntegrationPointsArrayType IntegrationPointsArrayType;

    typedef typename GeometryData::ShapeFunctionsGradientsType ShapeFunctionsGradientsType;

    typedef GeometryShapeFunctionContainer<IntegrationMethod> GeometryShapeFunctionContainerType;

    typedef typename GeometryType::IntegrationMethod IntegrationMethod;
    typedef typename GeometryType::IntegrationPointsContainerType IntegrationPointsContainerType;
    typedef typename GeometryType::ShapeFunctionsValuesContainerType ShapeFunctionsValuesContainerType;
    typedef typename GeometryType::ShapeFunctionsLocalGradientsContainerType ShapeFunctionsLocalGradientsContainerType;

    // Typedef for generic amount of derivatives.
    typedef typename GeometryType::ShapeFunctionsType ShapeFunctionsType;
    typedef typename GeometryType::ShapeFunctionsIntegrationPointsType ShapeFunctionsIntegrationPointsType;
    typedef typename GeometryType::ShapeFunctionsContainerType ShapeFunctionsContainerType;

    ///@}
    ///@name Life Cycle
    ///@{

    IntegrationPointSurface3d(
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

    IntegrationPointSurface3d(
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

    IntegrationPointSurface3d(
        const PointsArrayType& ThisPoints,
        const GeometryShapeFunctionContainerType& ThisGeometryShapeFunctionContainer)
        : BaseType(ThisPoints, &mGeometryData)
        , mGeometryData(
            &msGeometryDimension,
            ThisGeometryShapeFunctionContainer)
    {
    }

    IntegrationPointSurface3d(
        const PointsArrayType& ThisPoints,
        const GeometryShapeFunctionContainerType& ThisGeometryShapeFunctionContainer,
        GeometryType* pGeometryParent)
        : BaseType(ThisPoints, &mGeometryData)
        , mGeometryData(
            &msGeometryDimension,
            ThisGeometryShapeFunctionContainer)
        , mpGeometryParent(pGeometryParent)
    {
    }

    explicit IntegrationPointSurface3d(const PointsArrayType& ThisPoints)
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
    IntegrationPointSurface3d( IntegrationPointSurface3d const& rOther )
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
    template<class TOtherPointType> IntegrationPointSurface3d( IntegrationPointSurface3d<TOtherPointType> const& rOther )
        : BaseType( rOther )
    {
    }

    /**
     * Destructor. Does nothing!!!
     */
    ~IntegrationPointSurface3d() override {}

    GeometryData::KratosGeometryFamily GetGeometryFamily() const override
    {
        return GeometryData::Kratos_generic_family;
    }

    GeometryData::KratosGeometryType GetGeometryType() const override
    {
        return GeometryData::Kratos_generic_type;
    }

    /**
     * Operators
     */

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
    IntegrationPointSurface3d& operator=( const IntegrationPointSurface3d& rOther )
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
    IntegrationPointSurface3d& operator=( IntegrationPointSurface3d<TOtherPointType> const & rOther )
    {
        BaseType::operator=( rOther );

        return *this;
    }

    ///@}
    ///@name Operations
    ///@{
    typename BaseType::Pointer Create( PointsArrayType const& ThisPoints ) const override
    {
        return typename BaseType::Pointer( new IntegrationPointSurface3d( ThisPoints ) );
    }

    ///@}

    GeometryType& GetGeometryParent(IndexType Index) const override
    {
        return *mpGeometryParent;
    }

    /** Function returns the respective curve length on
    the underlying surface
    */
    double Area() const override
    {
        Vector temp;
        temp = this->DeterminantOfJacobian(temp, GeometryData::GI_GAUSS_1);
        const IntegrationPointsArrayType& r_integration_points = this->IntegrationPoints();
        double area = 0.0;

        for (std::size_t i = 0; i < r_integration_points.size(); ++i) {
            area += temp[i] * r_integration_points[i].Weight();
        }

        return area;
    }


    /** Returns the length of this curve segment.
    */
    double DomainSize() const override
    {
        return this->Area();
    }

    /** Calculates global location of this integration point.

    \f[
    c_i = \sum_j^n(x_j)*x_i
    \f]

    j is the index of the node and i the global direction (x,y,z).

    @return Point which is the location of this geometry.
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
    ///@name Information
    ///@{
    std::string Info() const override
    {
        return "2 dimensional single integration point defined in 3D space.";
    }

    /**
     * Print information about this object.
     * @param rOStream Stream to print into it.
     * @see PrintData()
     * @see Info()
     */
    void PrintInfo( std::ostream& rOStream ) const override
    {
        rOStream << "2 dimensional single integration point defined in 3D space.";
    }

    /**
     * Print geometry's data into given stream.
     * Prints it's points by the order they stored in
     * the geometry and then center point of geometry.
     *
     * @param rOStream Stream to print into it.
     * @see PrintInfo()
     * @see Info()
     */
    void PrintData( std::ostream& rOStream ) const override
    {
    }
    ///@}

protected:

private:
    ///@name Static Member Variables
    ///@{

    static const GeometryDimension msGeometryDimension;

    ///@}
    ///@name Member Variables
    ///@{

    GeometryData mGeometryData;

    // Integration point can be related to a parent geometry. To keep the connection,
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

    IntegrationPointSurface3d()
        : BaseType(
            PointsArrayType(),
            &mGeometryData)
        , mGeometryData(
            &msGeometryDimension,
            GeometryData::GI_GAUSS_1,
            {}, {}, {})
    {
    }
    /**
     * Private Friends
     */

    template<class TOtherPointType> friend class IntegrationPointSurface3d;

    /**
     * Un accessible methods
     */

}; // Class Geometry

/**
 * Input and output
 */
/**
 * input stream function
 */
template< class TPointType > inline std::istream& operator >> (
    std::istream& rIStream,
    IntegrationPointSurface3d<TPointType>& rThis );

/**
         * output stream function
 */
template<class TPointType> inline std::ostream& operator << (
    std::ostream& rOStream,
    const IntegrationPointSurface3d<TPointType>& rThis )
{
    rThis.PrintInfo( rOStream );
    rOStream << std::endl;
    rThis.PrintData( rOStream );
    return rOStream;
}

template<class TPointType>
const GeometryDimension IntegrationPointSurface3d<TPointType>::msGeometryDimension(
    2,
    3,
    2);

}  // namespace Kratos.

#endif // KRATOS_INTEGRATION_POINT_SURFACE_3D_H_INCLUDED  defined
