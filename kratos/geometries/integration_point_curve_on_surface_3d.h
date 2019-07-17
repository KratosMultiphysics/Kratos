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

#if !defined(KRATOS_INTEGRATION_POINT_CURVE_ON_SURFACE_3D_H_INCLUDED )
#define  KRATOS_INTEGRATION_POINT_CURVE_ON_SURFACE_3D_H_INCLUDED

// System includes

// External includes

// Project includes
#include "geometries/geometry.h"
#include "geometries/geometry_dimension.h"


namespace Kratos
{
/**
 * @class IntegrationPointCurveOnSurface3d
 * @ingroup KratosCore
 * @brief A sinlge integration point, that can be used for geometries without
 *        a predefined integration scheme, i.e. they can handle material point elements,
 *        isogeometric analysis elements or standard finite elements which are defined
 *        at a single integration point.
 *        This point defines a line segment described on a underlying surface.
 *        Shape functions and integration types have to be precomputed and are set from
 *        from outside.
 */
template<class TPointType> class IntegrationPointCurveOnSurface3d
    : public Geometry<TPointType>
{
public:

    /**
     * Pointer definition of IntegrationPointCurveOnSurface3d
     */
    KRATOS_CLASS_POINTER_DEFINITION( IntegrationPointCurveOnSurface3d );

    typedef Geometry<TPointType> BaseType;
    typedef Geometry<TPointType> GeometryType;

    typedef PointerVector<TPointType> PointsArrayType;
    typedef array_1d<double, 2> LocalCoordinatesArray2dType;

    typedef IntegrationPoint<2> IntegrationPointType;
    typedef typename GeometryType::IntegrationPointsArrayType IntegrationPointsArrayType;

    typedef typename GeometryData::ShapeFunctionsGradientsType ShapeFunctionsGradientsType;

    typedef array_1d<double, 2> LocalTangentsArray2dType;

    typedef typename GeometryType::IntegrationMethod IntegrationMethod;
    typedef typename GeometryType::IntegrationPointsContainerType IntegrationPointsContainerType;
    typedef typename GeometryType::ShapeFunctionsValuesContainerType ShapeFunctionsValuesContainerType;
    typedef typename GeometryType::ShapeFunctionsLocalGradientsContainerType ShapeFunctionsLocalGradientsContainerType;

    ///@}
    ///@name Life Cycle
    ///@{

    IntegrationPointCurveOnSurface3d(
        const PointsArrayType& ThisPoints,
        const IntegrationPointsContainerType& rIntegrationPoints,
        const LocalTangentsArray2dType& rLocalTangents2d,
        ShapeFunctionsValuesContainerType& rShapeFunctionValues,
        ShapeFunctionsLocalGradientsContainerType& rShapeFunctionsDerivativesVector)
        : BaseType(ThisPoints, &mGeometryData)
        , mLocalTangents2d(rLocalTangents2d)
        , mGeometryData(
            &msGeometryDimension,
            GeometryData::GI_GAUSS_1,
            rIntegrationPoints,
            rShapeFunctionValues,
            rShapeFunctionsDerivativesVector)
    {
    }

    explicit IntegrationPointCurveOnSurface3d(const PointsArrayType& ThisPoints)
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
    IntegrationPointCurveOnSurface3d( IntegrationPointCurveOnSurface3d const& rOther )
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
    template<class TOtherPointType> IntegrationPointCurveOnSurface3d(
        IntegrationPointCurveOnSurface3d<TOtherPointType> const& rOther )
        : BaseType( rOther )
    {
    }

    /**
     * Destructor. Does nothing!!!
     */
    ~IntegrationPointCurveOnSurface3d() override {}

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
    IntegrationPointCurveOnSurface3d& operator=( const IntegrationPointCurveOnSurface3d& rOther )
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
    IntegrationPointCurveOnSurface3d& operator=( IntegrationPointCurveOnSurface3d<TOtherPointType> const & rOther )
    {
        BaseType::operator=( rOther );

        return *this;
    }

    ///@}
    ///@name Operations
    ///@{
    typename BaseType::Pointer Create( PointsArrayType const& ThisPoints ) const override
    {
        return typename BaseType::Pointer( new IntegrationPointCurveOnSurface3d( ThisPoints ) );
    }

    ///@}
    ///@name Informations
    ///@{

    /** Function returns the respective curve length on
    the underlying surface
    */
    double Length() const override
    {
        Vector temp;
        temp = DeterminantOfJacobian(temp, GeometryData::GI_GAUSS_1);
        const IntegrationPointsArrayType& r_integration_points = this->IntegrationPoints();
        double length = 0.0;

        for (std::size_t i = 0; i < r_integration_points.size(); ++i) {
            length += temp[i] * r_integration_points[i].Weight();
        }

        return length;
    }


    /** Returns the length of this curve segment.
    */
    double DomainSize() const override
    {
        return this->Length();
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


    /** Determinant of jacobians for given integration method. This
    method calculate determinant of jacobian in all
    integrations points of given integration method.

    @return Vector of double which is vector of determinants of
    jacobians \f$ |J|_i \f$ where \f$ i=1,2,...,n \f$ is the
    integration point index of given integration method.

    @see Jacobian
    @see InverseOfJacobian
    */
    Vector& DeterminantOfJacobian(Vector& rResult, IntegrationMethod ThisMethod) const override
    {
        if (rResult.size() != this->IntegrationPointsNumber(ThisMethod))
            rResult.resize(this->IntegrationPointsNumber(ThisMethod), false);

        for (unsigned int pnt = 0; pnt < this->IntegrationPointsNumber(ThisMethod); pnt++)
        {
            rResult[pnt] = DeterminantOfJacobian(
                0, ThisMethod);
        }
        return rResult;
    }

    /** Determinant of jacobian in specific integration point of
    given integration method. This method calculate determinant
    of jacobian in given integration point of given integration
    method.

    The tangential integration weight is already applied to the
    length of the line segment.

    @param IntegrationPointIndex index of integration point which jacobians has to
    be calculated in it.

    @param IntegrationPointIndex index of integration point
    which determinant of jacobians has to be calculated in it.

    @return Determinamt of jacobian matrix \f$ |J|_i \f$ where \f$
    i \f$ is the given integration point index of given
    integration method.

    @see Jacobian
    @see InverseOfJacobian
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

    /**
    * @brief TO BE DONE
    */
    array_1d<double, 3> Normal(
        array_1d<double, 3>& rResult,
        IndexType IntegrationPointIndex,
        IntegrationMethod ThisMethod) const override
    {
        Matrix J;
        this->Jacobian(J, IntegrationPointIndex, ThisMethod);

        array_1d<double, 3> a_1 = column(J, 0);
        array_1d<double, 3> a_2 = column(J, 1);

        rResult = a_2 * mLocalTangents2d[0] - a_1 * mLocalTangents2d[1];

        return rResult;
    }

    /** Tangents in global space of the curve defined on the surface.

    The tangential integration weight is already applied to the
    length of the line segment.

    @see Jacobian
    @see InverseOfJacobian
    */
    virtual array_1d<double, 3> Tangent(
        array_1d<double, 3>& rResult,
        IndexType IntegrationPointIndex,
        IntegrationMethod ThisMethod) const
    {
        Matrix J;
        this->Jacobian(J, IntegrationPointIndex, ThisMethod);

        array_1d<double, 3> a_1 = column(J, 0);
        array_1d<double, 3> a_2 = column(J, 1);

        rResult = a_1 * mLocalTangents2d[0] + a_2 * mLocalTangents2d[1];

        return rResult;
    }

    ///@}
    ///@name Information
    ///@{
    std::string Info() const override
    {
        return "2 dimensional single line integration point defined in 3D space.";
    }

    /**
     * Print information about this object.
     * @param rOStream Stream to print into it.
     * @see PrintData()
     * @see Info()
     */
    void PrintInfo( std::ostream& rOStream ) const override
    {
        rOStream << "2 dimensional single line integration point defined in 3D space.";
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
        BaseType::PrintData( rOStream );
        std::cout << std::endl;
        Matrix jacobian;
        Jacobian( jacobian, TPointType() );
        rOStream << "    Jacobian in the integration point\t : " << jacobian;
    }
    ///@}

protected:

    /**
    * there are no protected class members
     */

private:
    ///@name Static Member Variables
    ///@{

    static const GeometryDimension msGeometryDimension;

    ///@}
    ///@name Private Member Variables
    ///@{

    GeometryData mGeometryData;

    LocalTangentsArray2dType mLocalTangents2d;

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

    IntegrationPointCurveOnSurface3d()
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

    template<class TOtherPointType> friend class IntegrationPointCurveOnSurface3d;

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
    IntegrationPointCurveOnSurface3d<TPointType>& rThis );

/**
         * output stream function
 */
template<class TPointType> inline std::ostream& operator << (
    std::ostream& rOStream,
    const IntegrationPointCurveOnSurface3d<TPointType>& rThis )
{
    rThis.PrintInfo( rOStream );
    rOStream << std::endl;
    rThis.PrintData( rOStream );
    return rOStream;
}

template<class TPointType>
const GeometryDimension IntegrationPointCurveOnSurface3d<TPointType>::msGeometryDimension(1,
    3,
    2);

}  // namespace Kratos.

#endif // KRATOS_INTEGRATION_POINT_CURVE_ON_SURFACE_3D_H_INCLUDED  defined
