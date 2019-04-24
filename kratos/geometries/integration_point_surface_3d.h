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
     * Type Definitions
     */
    /**
     * Base Type: Geometry
     */
    typedef Geometry<TPointType> BaseType;

    /**
     * Type of edge geometry
     */
    typedef Line3D3<TPointType> EdgeType;

    /**
     * Pointer definition of IntegrationPointSurface3d
     */
    KRATOS_CLASS_POINTER_DEFINITION( IntegrationPointSurface3d );

    /**
     * Integration methods implemented in geometry.
     */
    typedef GeometryData::IntegrationMethod IntegrationMethod;

    /**
     * A Vector of counted pointers to Geometries.
     * Used for returning edges of the geometry.
     */
    typedef typename BaseType::GeometriesArrayType GeometriesArrayType;

    /**
     * Redefinition of template parameter TPointType.
     */
    typedef TPointType PointType;

    /**
     * Array of coordinates. Can be Nodes, Points or IntegrationPointSurface3ds
     */
    typedef typename BaseType::CoordinatesArrayType CoordinatesArrayType;


    /**
     * Type used for indexing in geometry class.
     * std::size_t used for indexing
     * point or integration point access methods and also all other
     * methods which need point or integration point index.
     */
    typedef typename BaseType::IndexType IndexType;

    /**
     * This typed used to return size or dimension in
     * geometry. Dimension, WorkingDimension, PointsNumber and
     * ... return this type as their results.
     */
    typedef typename BaseType::SizeType SizeType;

    /**
     * Array of counted pointers to point.
     * This type used to hold geometry's points.
     */
    typedef  typename BaseType::PointsArrayType PointsArrayType;

    /**
     * This type used for representing an integration point in
     * geometry. This integration point is a point with an
     * additional weight component.
     */
    typedef typename BaseType::IntegrationPointType IntegrationPointType;

    /**
     * A Vector of IntegrationPointType which used to hold
     * integration points related to an integration
     * method. IntegrationPointSurface3ds functions used this type to return
     * their results.
     */
    typedef typename BaseType::IntegrationPointsArrayType IntegrationPointsArrayType;

    /**
     * A Vector of IntegrationPointsArrayType which used to hold
     * integration points related to different integration method
     * implemented in geometry.
     */
    typedef typename BaseType::IntegrationPointsContainerType IntegrationPointsContainerType;

    /**
     * A third order tensor used as shape functions' values
     * container.
     */
    typedef typename BaseType::ShapeFunctionsValuesContainerType
    ShapeFunctionsValuesContainerType;

    /**
             * A fourth order tensor used as shape functions' local
             * gradients container in geometry.
     */
    typedef typename BaseType::ShapeFunctionsLocalGradientsContainerType
    ShapeFunctionsLocalGradientsContainerType;

    /**
             * A third order tensor to hold jacobian matrices evaluated at
             * integration points. Jacobian and InverseOfJacobian functions
             * return this type as their result.
     */
    typedef typename BaseType::JacobiansType JacobiansType;

    /**
     * A third order tensor to hold shape functions' local
     * gradients. ShapefunctionsLocalGradients function return this
     * type as its result.
     */
    typedef typename BaseType::ShapeFunctionsGradientsType ShapeFunctionsGradientsType;

    /**
     * A third order tensor to hold shape functions' local second derivatives.
     * ShapefunctionsLocalGradients function return this
     * type as its result.
     */
    typedef typename BaseType::ShapeFunctionsSecondDerivativesType
        ShapeFunctionsSecondDerivativesType;


    /**
    * A fourth order tensor to hold shape functions' local third derivatives.
    * ShapefunctionsLocalGradients function return this
    * type as its result.
     */
    typedef typename BaseType::ShapeFunctionsThirdDerivativesType
        ShapeFunctionsThirdDerivativesType;

    /**
    * Type of the normal vector used for normal to edges in geomety.
     */
    typedef typename BaseType::NormalType NormalType;

    ///@}
    ///@name Life Cycle
    ///@{

    IntegrationPointSurface3d(
        const PointsArrayType& ThisPoints,
        const CoordinatesArrayType& rLocalCoordinates,
        const double& rIntegrationWeight,
        Matrix& rShapeFunctions)
        : BaseType( ThisPoints, nullptr )
    {
        IntegrationPoint(rLocalCoordinates[0], rLocalCoordinates[1], rIntegrationWeight);

        mIntegrationPoints = IntegrationPointsContainerType(1);
        mIntegrationPoints[0] = IntegrationPoint;

        KRATOS_ERROR_IF(rShapeFunctions.size1() != ThisPoints.size()) <<
            "Length of shape function array in matrix does not fit to number of points. Length of point vector: " << ThisPoints.size()
            << " Length of corresponding shape function vector: " << rShapeFunctions.size2() << std::endl;
        mShapeFunctions = rShapeFunctions;
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
    ///@name Integration Points
    ///@{

    /** IntegrationPointSurface3d is a geometry which is based on a single
    gauss point. Thus, the size is always 1. Independent on the integration
    method.

    @return SizeType which is always 1.
    */
    SizeType IntegrationPointsNumber() const override
    {
        return 1;
    }

    /** IntegrationPointSurface3d is a geometry which is based on a single
    gauss point. Thus, the size is always 1. Independent on the integration
    method.

    @return SizeType which is always 1.
    */
    SizeType IntegrationPointsNumber(IntegrationMethod ThisMethod) const override
    {
        return 1;
    }


    /** This method returns the member variable of the precomputed integration point.

    @return const IntegrationPointsArrayType which is Vector of size 1
    of this integration points
    for default integrating method.
    */
    const IntegrationPointsArrayType& IntegrationPoints() const override
    {
        return mIntegrationPoints;
    }

    /** Integtation points for given integration
    method. This method use integration points data base to
    obtain integration points Vector respected to
    given method.

    @return const IntegrationPointsArrayType which is Vector of integration points
    for default integrating method.
    */
    const IntegrationPointsArrayType& IntegrationPoints(IntegrationMethod ThisMethod) const override
    {
        return mIntegrationPoints;
    }

    ///@}
    ///@name Informations
    ///@{

    /** Dimension of the Integration Point, is always set to 3.

    @return SizeType, dimension of this geometry.
    @see WorkingSpaceDimension()
    @see LocalSpaceDimension()
    */
    inline SizeType Dimension() const override
    {
        return 3;
    }

    /** Working Space Dimension of the Integration Point, is always set to 3.

    @return SizeType, working space dimension of this geometry.
    @see Dimension()
    @see LocalSpaceDimension()
    */
    inline SizeType WorkingSpaceDimension() const override
    {
        return 3();
    }

    /** Local Space Dimension of surface integration point are always the 
    two coordinates of the surface.

    @return SizeType, local space dimension of this geometry.
    @see Dimension()
    @see WorkingSpaceDimension()
    */
    inline SizeType LocalSpaceDimension() const override
    {
        return 2;
    }

    ///@}

    /** Calculates global location of this integration point.

    \f[
    c_i = \sum_j^n(x_j)*x_i
    \f]

    j is the index of the node and i the global direction (x,y,z).

    @return Point which is the location of this geometry.
    */
    virtual Point Center() const
    {
        const SizeType points_number = PointsNumber();

        Point location = ZeroVector(3);
        const Matrix& ShapeFunctionValues();

        for (IndexType i = 0; i < PointsNumber(); ++i) {
            location.Coordinates() += (*this)[i].Coordinates()* Matrix(0, i);
        }
        return location;
    }


    ///@name Shape Function
    ///@{
    /** This method gives all non-zero shape functions values
    evaluated at the rCoordinates provided

    @return Vector of values of shape functions \f$ F_{i} \f$
    where i is the shape function index (for NURBS it is the index
    of the local enumeration in the element).

    @see ShapeFunctionValue
    @see ShapeFunctionsLocalGradients
    @see ShapeFunctionLocalGradient
    */
    double ShapeFunctionValue(
        IndexType IntegrationPointIndex,
        IndexType ShapeFunctionIndex
    ) const override
    {
        mShapeFunctionValues(0, ShapeFunctionIndex)

        return 0;
    }

    /** This method gives the shape functions values corresponding to
    this integration point evaluated at each nodes. There is no 
    calculation of the shape functions it just returns it from 
    member container.

    @return Matrix of values of shape functions \f$ F_{ij} \f$
    where i is the integration point index and j is the shape
    function index. In other word component \f$ f_{ij} \f$ is value
    of the shape function corresponding to node j evaluated in
    integration point i of default integration method.

    @see ShapeFunctionValue
    @see ShapeFunctionsLocalGradients
    @see ShapeFunctionLocalGradient
    */
    const Matrix& ShapeFunctionsValues() const
    {
        Matrix ShapeFunctionsValues = ZeroMatrix(0, PointsNumber());

        for (IndexType i = 0; i < PointsNumber(); ++i)
        {
            ShapeFunctionsValues(0, i) = mShapeFunctions(0, i);
        }

        return ShapeFunctionsValues();
    }

    /** This method gives the shape functions values corresponding to
    this integration point evaluated at each nodes.

    \note The integration at one point is independent on the method.

    @see ShapeFunctionsValues
    */
    const Matrix& ShapeFunctionsValues(IntegrationMethod ThisMethod)  const
    {

        return ShapeFunctionsValues();
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
    ///@name Members
    ///@{

    IntegrationPointsContainerType mIntegrationPoints;

    /** This member holds the shape functions and its evaluated
    *   derivative  */
    Matrix mShapeFunctions;

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

    IntegrationPointSurface3d(): BaseType( PointsArrayType(), nullptr ) {}

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

}  // namespace Kratos.

#endif // KRATOS_INTEGRATION_POINT_SURFACE_3D_H_INCLUDED  defined 
