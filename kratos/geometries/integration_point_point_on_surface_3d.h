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

#if !defined(KRATOS_INTEGRATION_POINT_POINT_ON_SURFACE_3D_H_INCLUDED )
#define  KRATOS_INTEGRATION_POINT_POINT_ON_SURFACE_3D_H_INCLUDED

// System includes

// External includes

// Project includes
#include "geometries/geometry.h"
#include "geometries/geometry_dimension.h"


namespace Kratos
{
/**
 * @class IntegrationPointPointOnSurface3d
 * @ingroup KratosCore
 * @brief A sinlge integration point, that can be used for geometries without 
 *        a predefined integration scheme, i.e. they can handle material point elements,
 *        isogeometric analysis elements or standard finite elements which are defined
 *        at a single integration point.
 *        This point defines a line segment described on a underlying surface.
 *        Shape functions and integration types have to be precomputed and are set from 
 *        from outside.
 */
template<class TPointType> class IntegrationPointPointOnSurface3d
    : public Geometry<TPointType>
{
public:

    /**
     * Pointer definition of IntegrationPointPointOnSurface3d
     */
    KRATOS_CLASS_POINTER_DEFINITION( IntegrationPointPointOnSurface3d );

    typedef Geometry<TPointType> BaseType;
    typedef Geometry<TPointType> GeometryType;

    typedef PointerVector<TPointType> PointsArrayType;
    typedef typename PointType::CoordinatesArrayType CoordinatesArrayType;

    typedef IntegrationPoint<2> IntegrationPointType;
    typedef std::vector<IntegrationPointType> IntegrationPointsArrayType;

    typedef GeometryData::ShapeFunctionsGradientsType ShapeFunctionsGradientsType;

    ///@}
    ///@name Life Cycle
    ///@{

    IntegrationPointPointOnSurface3d(
        const PointsArrayType& ThisPoints,
        const CoordinatesArrayType& rLocalCoordinates,
        Matrix& rShapeFunctionValues,
        ShapeFunctionsGradientsType& rShapeFunctionsDerivativesVector)
        : BaseType(ThisPoints, &mGeometryData)
    {
        IntegrationPoint(rLocalCoordinates[0], rLocalCoordinates[1], 1.0);

        IntegrationPointsArrayType IntegrationPoints = IntegrationPointsContainerType(1);
        IntegrationPoints[0] = IntegrationPoint;

        mGeometryData = GeometryData(
            msGeometryDimension,
            GeometryData::GI_GAUSS_1,
            IntegrationPoints,
            rShapeFunctionValues,
            rShapeFunctionsDerivativesVector);
    }

    IntegrationPointPointOnSurface3d(
        const PointsArrayType& ThisPoints,
        const IntegrationPointsArrayType& rIntegrationPoints,
        Matrix& rShapeFunctionValues,
        ShapeFunctionsGradientsType& rShapeFunctionsDerivativesVector)
        : BaseType(ThisPoints, &mGeometryData)
    {
        mGeometryData = GeometryData(
            msGeometryDimension,
            GeometryData::GI_GAUSS_1,
            rIntegrationPoints,
            rShapeFunctionValues,
            rShapeFunctionsDerivativesVector);
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
    IntegrationPointPointOnSurface3d( IntegrationPointPointOnSurface3d const& rOther )
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
    template<class TOtherPointType> IntegrationPointPointOnSurface3d(
        IntegrationPointPointOnSurface3d<TOtherPointType> const& rOther )
        : BaseType( rOther )
    {
    }

    /**
     * Destructor. Does nothing!!!
     */
    ~IntegrationPointPointOnSurface3d() override {}

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
    IntegrationPointPointOnSurface3d& operator=( const IntegrationPointPointOnSurface3d& rOther )
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
    IntegrationPointPointOnSurface3d& operator=( IntegrationPointPointOnSurface3d<TOtherPointType> const & rOther )
    {
        BaseType::operator=( rOther );

        return *this;
    }

    ///@}
    ///@name Operations
    ///@{
    typename BaseType::Pointer Create( PointsArrayType const& ThisPoints ) const override
    {
        return typename BaseType::Pointer( new IntegrationPointPointOnSurface3d( ThisPoints ) );
    }

    ///@}
    /** Calculates global location of this integration point.

    \f[
    c_i = \sum_j^n(x_j)*x_i
    \f]

    j is the index of the node and i the global direction (x,y,z).

    @return Point which is the location of this geometry.
    */
    Point Center() const override
    {
        const SizeType points_number = PointsNumber();

        Point point(0.0, 0.0, 0.0);
        const Matrix& N = this->ShapeFunctionsValues();

        for (IndexType point_number = 0; point_number < IntegrationPointsNumber(); ++point_number) {
            for (IndexType i = 0; i < PointsNumber(); ++i) {
                point += (*this)[i] * N(point_number, i);
            }
        }
        return point;
    }


    /** Determinant of point is always 1.0, independent on the deformation.
    */
    Vector& DeterminantOfJacobian(Vector& rResult, IntegrationMethod ThisMethod) const override
    {
        if (rResult.size() != this->IntegrationPointsNumber(ThisMethod))
            rResult.resize(this->IntegrationPointsNumber(ThisMethod), false);

        for (unsigned int pnt = 0; pnt < this->IntegrationPointsNumber(ThisMethod); pnt++)
        {
            rResult[pnt] = 1.0;
        }
        return rResult;
    }

    /** Determinant of point is always 1.0, independent on the deformation.
    */
    double DeterminantOfJacobian(
        IndexType IntegrationPointIndex,
        IntegrationMethod ThisMethod) const override
    {
        return 1.0;
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
        Jacobian( jacobian, PointType() );
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

    array_1d<double, 2> mTangents;

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

    IntegrationPointPointOnSurface3d(): BaseType( PointsArrayType(), &mGeometryData ) {}

    /**
     * Private Friends
     */

    template<class TOtherPointType> friend class IntegrationPointPointOnSurface3d;

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
    IntegrationPointPointOnSurface3d<TPointType>& rThis );

/**
         * output stream function
 */
template<class TPointType> inline std::ostream& operator << (
    std::ostream& rOStream,
    const IntegrationPointPointOnSurface3d<TPointType>& rThis )
{
    rThis.PrintInfo( rOStream );
    rOStream << std::endl;
    rThis.PrintData( rOStream );
    return rOStream;
}

template<class TPointType>
const GeometryDimension IntegrationPointPointOnSurface3d<TPointType>::msGeometryDimension(
    0,
    3,
    2);

}  // namespace Kratos.

#endif // KRATOS_INTEGRATION_POINT_POINT_ON_SURFACE_3D_H_INCLUDED  defined 
