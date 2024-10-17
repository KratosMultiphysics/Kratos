//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//                   Janosch Stascheit
//                   Felix Nagel
//  contributors:    Hoang Giang Bui
//                   Josep Maria Carbonell
//

#pragma once

// System includes

// External includes

// Project includes
#include "geometries/geometry.h"
#include "integration/line_gauss_legendre_integration_points.h"
#include "integration/line_collocation_integration_points.h"
#include "utilities/intersection_utilities.h"
#include "utilities/geometry_utilities.h"
namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/**
 * @class Line3D2
 * @ingroup KratosCore
 * @brief An two node 3D line geometry with linear shape functions
 * @details The node ordering corresponds with:
 *      0----------1 --> u
 * @author Riccardo Rossi
 * @author Janosch Stascheit
 * @author Felix Nagel
 */
template<class TPointType>

class Line3D2 : public Geometry<TPointType>
{

public:
    ///@}
    ///@name Type Definitions
    ///@{

    /// Geometry as base class.
    typedef Geometry<TPointType> BaseType;
    using Geometry<TPointType>::ShapeFunctionsValues;

    /// Pointer definition of Line3D2
    KRATOS_CLASS_POINTER_DEFINITION( Line3D2 );

    /// Type of edge geometry
    typedef Line3D2<TPointType> EdgeType;

    /** Integration methods implemented in geometry.
    */
    typedef GeometryData::IntegrationMethod IntegrationMethod;

    /** A Vector of counted pointers to Geometries. Used for
    returning edges of the geometry.
     */
    typedef typename BaseType::GeometriesArrayType GeometriesArrayType;

    /** Redefinition of template parameter TPointType.
     */
    typedef TPointType PointType;

    /** Type used for indexing in geometry class.std::size_t used for indexing
    point or integration point access methods and also all other
    methods which need point or integration point index.
    */
    typedef typename BaseType::IndexType IndexType;


    /** This typed used to return size or dimension in
    geometry. Dimension, WorkingDimension, PointsNumber and
    ... return this type as their results.
    */
    typedef typename BaseType::SizeType SizeType;

    /** Array of counted pointers to point. This type used to hold
    geometry's points.
    */
    typedef  typename BaseType::PointsArrayType PointsArrayType;

    /** This type used for representing an integration point in
    geometry. This integration point is a point with an
    additional weight component.
    */
    typedef typename BaseType::IntegrationPointType IntegrationPointType;

    /** A Vector of IntegrationPointType which used to hold
    integration points related to an integration
    method. IntegrationPoints functions used this type to return
    their results.
    */
    typedef typename BaseType::IntegrationPointsArrayType IntegrationPointsArrayType;

    /** A Vector of IntegrationPointsArrayType which used to hold
    integration points related to different integration method
    implemented in geometry.
    */
    typedef typename BaseType::IntegrationPointsContainerType IntegrationPointsContainerType;

    /** A third order tensor used as shape functions' values
    continer.
    */
    typedef typename BaseType::ShapeFunctionsValuesContainerType ShapeFunctionsValuesContainerType;

    /** A fourth order tensor used as shape functions' local
    gradients container in geometry.
    */
    typedef typename BaseType::ShapeFunctionsLocalGradientsContainerType ShapeFunctionsLocalGradientsContainerType;

    /** A third order tensor to hold jacobian matrices evaluated at
    integration points. Jacobian and InverseOfJacobian functions
    return this type as their result.
    */
    typedef typename BaseType::JacobiansType JacobiansType;

    /** A third order tensor to hold shape functions' local
    gradients. ShapefunctionsLocalGradients function return this
    type as its result.
    */
    typedef typename BaseType::ShapeFunctionsGradientsType ShapeFunctionsGradientsType;

    /** Type of the normal vector used for normal to edges in geomety.
     */
    typedef typename BaseType::NormalType NormalType;

    /**
     * Type of coordinates array
     */
    typedef typename BaseType::CoordinatesArrayType CoordinatesArrayType;

    ///@}
    ///@name Life Cycle
    ///@{

//     Line3D2( const PointType& FirstPoint, const PointType& SecondPoint )
//         : BaseType( PointsArrayType(), &msGeometryData )
//     {
//         this->Points().push_back( typename PointType::Pointer( new PointType( FirstPoint ) ) );
//         this->Points().push_back( typename PointType::Pointer( new PointType( SecondPoint ) ) );
// //         BaseType::Points().push_back( typename PointType::Pointer( new PointType( FirstPoint ) ) );
// //         BaseType::Points().push_back( typename PointType::Pointer( new PointType( SecondPoint ) ) );
//     }

    Line3D2( typename PointType::Pointer pFirstPoint, typename PointType::Pointer pSecondPoint )
        : BaseType( PointsArrayType(), &msGeometryData )
    {
        BaseType::Points().push_back( pFirstPoint );
        BaseType::Points().push_back( pSecondPoint );
    }

    explicit Line3D2( const PointsArrayType& ThisPoints )
        : BaseType( ThisPoints, &msGeometryData )
    {
        if ( BaseType::PointsNumber() != 2 )
            KRATOS_ERROR << "Invalid points number. Expected 2, given " << BaseType::PointsNumber() << std::endl;
    }

    /// Constructor with Geometry Id
    explicit Line3D2(
        const IndexType GeometryId,
        const PointsArrayType& ThisPoints
    ) : BaseType(GeometryId, ThisPoints, &msGeometryData)
    {
        KRATOS_ERROR_IF( this->PointsNumber() != 2 ) << "Invalid points number. Expected 2, given " << this->PointsNumber() << std::endl;
    }

    /// Constructor with Geometry Name
    explicit Line3D2(
        const std::string& rGeometryName,
        const PointsArrayType& rThisPoints
    ) : BaseType(rGeometryName, rThisPoints, &msGeometryData)
    {
        KRATOS_ERROR_IF(this->PointsNumber() != 2) << "Invalid points number. Expected 2, given " << this->PointsNumber() << std::endl;
    }

    /** Copy constructor.
    Construct this geometry as a copy of given geometry.

    @note This copy constructor don't copy the points and new
    geometry shares points with given source geometry. It's
    obvious that any change to this new geometry's point affect
    source geometry's points too.
    */
    Line3D2( Line3D2 const& rOther )
        : BaseType( rOther )
    {
    }


    /** Copy constructor from a geometry with other point type.
    Construct this geometry as a copy of given geometry which
    has different type of points. The given goemetry's
    TOtherPointType* must be implicity convertible to this
    geometry PointType.

    @note This copy constructor don't copy the points and new
    geometry shares points with given source geometry. It's
    obvious that any change to this new geometry's point affect
    source geometry's points too.
    */
    template<class TOtherPointType> explicit Line3D2( Line3D2<TOtherPointType> const& rOther )
        : BaseType( rOther )
    {
    }

    /// Destructor. Do nothing!!!
    ~Line3D2() override {}

    GeometryData::KratosGeometryFamily GetGeometryFamily() const override
    {
        return GeometryData::KratosGeometryFamily::Kratos_Linear;
    }

    GeometryData::KratosGeometryType GetGeometryType() const override
    {
        return GeometryData::KratosGeometryType::Kratos_Line3D2;
    }

    ///@}
    ///@name Operators
    ///@{

    /** Assignment operator.

    @note This operator don't copy the points and this
    geometry shares points with given source geometry. It's
    obvious that any change to this geometry's point affect
    source geometry's points too.

    @see Clone
    @see ClonePoints
    */
    Line3D2& operator=( const Line3D2& rOther )
    {
        BaseType::operator=( rOther );

        return *this;
    }

    /** Assignment operator for geometries with different point type.

    @note This operator don't copy the points and this
    geometry shares points with given source geometry. It's
    obvious that any change to this geometry's point affect
    source geometry's points too.

    @see Clone
    @see ClonePoints
    */
    template<class TOtherPointType>
    Line3D2& operator=( Line3D2<TOtherPointType> const & rOther )
    {
        BaseType::operator=( rOther );

        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Creates a new geometry pointer
     * @param NewGeometryId the ID of the new geometry
     * @param rThisPoints the nodes of the new geometry
     * @return Pointer to the new geometry
     */
    typename BaseType::Pointer Create(
        const IndexType NewGeometryId,
        PointsArrayType const& rThisPoints
        ) const override
    {
        return typename BaseType::Pointer( new Line3D2( NewGeometryId, rThisPoints ) );
    }

    /**
     * @brief Creates a new geometry pointer
     * @param NewGeometryId the ID of the new geometry
     * @param rGeometry reference to an existing geometry
     * @return Pointer to the new geometry
     */
    typename BaseType::Pointer Create(
        const IndexType NewGeometryId,
        const BaseType& rGeometry
    ) const override
    {
        auto p_geometry = typename BaseType::Pointer( new Line3D2( NewGeometryId, rGeometry.Points() ) );
        p_geometry->SetData(rGeometry.GetData());
        return p_geometry;
    }

    /**
     * @brief Lumping factors for the calculation of the lumped mass matrix
     * @param rResult Vector containing the lumping factors
     * @param LumpingMethod The lumping method considered. The three methods available are:
     *      - The row sum method
     *      - Diagonal scaling
     *      - Evaluation of M using a quadrature involving only the nodal points and thus automatically yielding a diagonal matrix for standard element shape function
     */
    Vector& LumpingFactors(
        Vector& rResult,
        const typename BaseType::LumpingMethods LumpingMethod = BaseType::LumpingMethods::ROW_SUM
        )  const override
    {
	if(rResult.size() != 2)
           rResult.resize( 2, false );
        rResult[0] = 0.5;
        rResult[1] = 0.5;
        return rResult;
    }

    ///@}
    ///@name Informations
    ///@{

    /** This method calculate and return Length or charactereistic
    length of this geometry depending to it's dimension. For one
    dimensional geometry for example Line it returns length of it
    and for the other geometries it gives Characteristic length
    otherwise.

    @return double value contains length or Characteristic
    length
    @see Area()
    @see Volume()
    @see DomainSize()
    */
    double Length() const override
    {
        const TPointType& point0 = BaseType::GetPoint(0);
        const TPointType& point1 = BaseType::GetPoint(1);
        const double lx = point0.X() - point1.X();
        const double ly = point0.Y() - point1.Y();
        const double lz = point0.Z() - point1.Z();

        const double length = lx * lx + ly * ly + lz * lz;

        return std::sqrt( length );
    }

    /** This method calculate and return area or surface area of
    this geometry depending to it's dimension. For one dimensional
    geometry it returns length, for two dimensional it gives area
    and for three dimensional geometries it gives surface area.

    @return double value contains area or surface
    area.
    @see Length()
    @see Volume()
    @see DomainSize()
    */
    double Area() const override
    {
        return Length();
    }


    /** This method calculate and return length, area or volume of
    this geometry depending to it's dimension. For one dimensional
    geometry it returns its length, for two dimensional it gives area
    and for three dimensional geometries it gives its volume.

    @return double value contains length, area or volume.
    @see Length()
    @see Area()
    @see Volume()
    */
    double DomainSize() const override
    {
        return Length();
    }


    ///@}
    ///@name Intersections
    ///@{

    /**
     * @brief Test if this geometry intersects with other geometry
     * @param  ThisGeometry Geometry to intersect with
     * @return True if the geometries intersect, False in any other case.
     * @details We always check the intersection from the higher LocalSpaceDimension to the lower one
     */
    bool HasIntersection(const BaseType& rThisGeometry) const override
    {
        const BaseType& r_geom = *this;
        if (rThisGeometry.LocalSpaceDimension() > r_geom.LocalSpaceDimension()) {
            return rThisGeometry.HasIntersection(r_geom);
        }
        // Both objects are lines
        Point intersection_point;
        return IntersectionUtilities::ComputeLineLineIntersection(r_geom, rThisGeometry[0], rThisGeometry[1], intersection_point);
    }

    ///@}
    ///@name Spatial Operations
    ///@{

    /**
    * @brief Computes the distance between an point in
    *        global coordinates and the closest point
    *        of this geometry.
    *        If projection fails, double::max will be returned.
    * @param rPointGlobalCoordinates the point to which the
    *        closest point has to be found.
    * @param Tolerance accepted orthogonal error.
    * @return Distance to geometry.
    *         positive -> outside of to the geometry (for 2D and solids)
    *         0        -> on/ in the geometry.
    */
    double CalculateDistance(
        const CoordinatesArrayType& rPointGlobalCoordinates,
        const double Tolerance = std::numeric_limits<double>::epsilon()
        ) const override
    {
        const Point point(rPointGlobalCoordinates);
        return GeometryUtils::PointDistanceToLineSegment3D(this->GetPoint(0), this->GetPoint(1), point);
    }

    ///@}
    ///@name Jacobian
    ///@{

    /** Jacobians for given  method. This method
    calculate jacobians matrices in all integrations points of
    given integration method.

    @param ThisMethod integration method which jacobians has to
    be calculated in its integration points.

    @return JacobiansType a Vector of jacobian
    matrices \f$ J_i \f$ where \f$ i=1,2,...,n \f$ is the integration
    point index of given integration method.

    @see DeterminantOfJacobian
    @see InverseOfJacobian
    */
    JacobiansType& Jacobian( JacobiansType& rResult, IntegrationMethod ThisMethod ) const override
    {
        Matrix jacobian( 3, 1 );
        jacobian( 0, 0 ) = ( this->GetPoint( 1 ).X() - this->GetPoint( 0 ).X() ) * 0.5; //on the Gauss points (J is constant at each element)
        jacobian( 1, 0 ) = ( this->GetPoint( 1 ).Y() - this->GetPoint( 0 ).Y() ) * 0.5;
        jacobian( 2, 0 ) = ( this->GetPoint( 1 ).Z() - this->GetPoint( 0 ).Z() ) * 0.5;

        if ( rResult.size() != BaseType::IntegrationPointsNumber( ThisMethod ) )
        {
            // KLUDGE: While there is a bug in ublas vector resize, I have to put this beside resizing!!
            JacobiansType temp( BaseType::IntegrationPointsNumber( ThisMethod ) );
            rResult.swap( temp );
        }

        std::fill( rResult.begin(), rResult.end(), jacobian );

        return rResult;
    }

    /** Jacobians for given  method. This method
    calculate jacobians matrices in all integrations points of
    given integration method.

    @param ThisMethod integration method which jacobians has to
    be calculated in its integration points.

    @return JacobiansType a Vector of jacobian
    matrices \f$ J_i \f$ where \f$ i=1,2,...,n \f$ is the integration
    point index of given integration method.

    @param DeltaPosition Matrix with the nodes position increment which describes
    the configuration where the jacobian has to be calculated.

    @see DeterminantOfJacobian
    @see InverseOfJacobian
    */
    JacobiansType& Jacobian( JacobiansType& rResult, IntegrationMethod ThisMethod, Matrix & DeltaPosition ) const override
    {
        Matrix jacobian( 3, 1 );
        jacobian( 0, 0 ) = ( (this->GetPoint( 1 ).X() - DeltaPosition(1,0)) - (this->GetPoint( 0 ).X() - DeltaPosition(0,0)) ) * 0.5; //on the Gauss points (J is constant at each element)
	jacobian( 1, 0 ) = ( (this->GetPoint( 1 ).Y() - DeltaPosition(1,1)) - (this->GetPoint( 0 ).Y() - DeltaPosition(0,1)) ) * 0.5;
        jacobian( 2, 0 ) = ( (this->GetPoint( 1 ).Z() - DeltaPosition(1,2)) - (this->GetPoint( 0 ).Z() - DeltaPosition(0,2)) ) * 0.5;

        if ( rResult.size() != BaseType::IntegrationPointsNumber( ThisMethod ) )
        {
            // KLUDGE: While there is a bug in ublas vector resize, I have to put this beside resizing!!
            JacobiansType temp( BaseType::IntegrationPointsNumber( ThisMethod ) );
            rResult.swap( temp );
        }

        std::fill( rResult.begin(), rResult.end(), jacobian );

        return rResult;
    }

    /** Jacobian in specific integration point of given integration
    method. This method calculate jacobian matrix in given
    integration point of given integration method.

    @param IntegrationPointIndex index of integration point which jacobians has to
    be calculated in it.

    @param ThisMethod integration method which jacobians has to
    be calculated in its integration points.

    @return Matrix<double> Jacobian matrix \f$ J_i \f$ where \f$
    i \f$ is the given integration point index of given
    integration method.

    @see DeterminantOfJacobian
    @see InverseOfJacobian
    */
    Matrix& Jacobian( Matrix& rResult, IndexType IntegrationPointIndex, IntegrationMethod ThisMethod ) const override
    {
        rResult.resize( 3, 1, false );
        //on the Gauss points (J is constant at each element)
        rResult( 0, 0 ) = ( this->GetPoint( 1 ).X() - this->GetPoint( 0 ).X() ) * 0.5;
        rResult( 1, 0 ) = ( this->GetPoint( 1 ).Y() - this->GetPoint( 0 ).Y() ) * 0.5;
        rResult( 2, 0 ) = ( this->GetPoint( 1 ).Z() - this->GetPoint( 0 ).Z() ) * 0.5;
        return rResult;
    }

    /** Jacobian in given point. This method calculate jacobian
    matrix in given point.

    @param rPoint point which jacobians has to
    be calculated in it.

    @return Matrix of double which is jacobian matrix \f$ J \f$ in given point.

    @see DeterminantOfJacobian
    @see InverseOfJacobian
    */
    Matrix& Jacobian( Matrix& rResult, const CoordinatesArrayType& rPoint ) const override
    {
        rResult.resize( 3, 1, false );
        //on the Gauss points (J is constant at each element)
        rResult( 0, 0 ) = ( this->GetPoint( 1 ).X() - this->GetPoint( 0 ).X() ) * 0.5;
        rResult( 1, 0 ) = ( this->GetPoint( 1 ).Y() - this->GetPoint( 0 ).Y() ) * 0.5;
        rResult( 2, 0 ) = ( this->GetPoint( 1 ).Z() - this->GetPoint( 0 ).Z() ) * 0.5;
        return rResult;
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
    Vector& DeterminantOfJacobian( Vector& rResult, IntegrationMethod ThisMethod ) const override
    {
        const unsigned int integration_points_number = msGeometryData.IntegrationPointsNumber( ThisMethod );
        if(rResult.size() != integration_points_number)
        {
            rResult.resize(integration_points_number,false);
        }

        const double detJ = 0.5*(this->Length());

        for ( unsigned int pnt = 0; pnt < integration_points_number; pnt++ )
        {
            rResult[pnt] = detJ;
        }
        return rResult;
    }

    /** Determinant of jacobian in specific integration point of
    given integration method. This method calculate determinant
    of jacobian in given integration point of given integration
    method.

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
    double DeterminantOfJacobian( IndexType IntegrationPointIndex, IntegrationMethod ThisMethod ) const override
    {
        return 0.5*(this->Length());
    }

    /** Determinant of jacobian in given point. This method calculate determinant of jacobian
    matrix in given point.

    @param rPoint point which determinant of jacobians has to
    be calculated in it.

    @return Determinamt of jacobian matrix \f$ |J| \f$ in given
    point.

    @see DeterminantOfJacobian
    @see InverseOfJacobian
    */
    double DeterminantOfJacobian( const CoordinatesArrayType& rPoint ) const override
    {
        return 0.5*(this->Length());
    }


    /** Inverse of jacobians for given integration method. This method
    calculate inverse of jacobians matrices in all integrations points of
    given integration method.

    @param ThisMethod integration method which inverse of jacobians has to
    be calculated in its integration points.

    @return Inverse of jacobian
    matrices \f$ J^{-1}_i \f$ where \f$ i=1,2,...,n \f$ is the integration
    point index of given integration method.

    @see Jacobian
    @see DeterminantOfJacobian
    */
    JacobiansType& InverseOfJacobian( JacobiansType& rResult, IntegrationMethod ThisMethod ) const override
    {
        rResult[0] = ZeroMatrix( 1, 1 );
        rResult[0]( 0, 0 ) = 2.0 * MathUtils<double>::Norm3(( this->GetPoint( 1 ) ) - ( this->GetPoint( 0 ) ) );
        return rResult;
    }

    /** Inverse of jacobian in specific integration point of given integration
    method. This method calculate Inverse of jacobian matrix in given
    integration point of given integration method.

    @param IntegrationPointIndex index of integration point which inverse of jacobians has to
    be calculated in it.

    @param ThisMethod integration method which inverse of jacobians has to
    be calculated in its integration points.

    @return Inverse of jacobian matrix \f$ J^{-1}_i \f$ where \f$
    i \f$ is the given integration point index of given
    integration method.

    @see Jacobian
    @see DeterminantOfJacobian
    */
    Matrix& InverseOfJacobian( Matrix& rResult, IndexType IntegrationPointIndex, IntegrationMethod ThisMethod ) const override
    {
        rResult = ZeroMatrix( 1, 1 );
        rResult( 0, 0 ) = 2.0 * MathUtils<double>::Norm3(( this->GetPoint( 1 ) ) - ( this->GetPoint( 0 ) ) );
        return( rResult );
    }

    /** Inverse of jacobian in given point. This method calculate inverse of jacobian
    matrix in given point.

    @param rPoint point which inverse of jacobians has to
    be calculated in it.

    @return Inverse of jacobian matrix \f$ J^{-1} \f$ in given point.

    @see DeterminantOfJacobian
    @see InverseOfJacobian
    */
    Matrix& InverseOfJacobian( Matrix& rResult, const CoordinatesArrayType& rPoint ) const override
    {
        rResult = ZeroMatrix( 1, 1 );
        rResult( 0, 0 ) = 2.0 * MathUtils<double>::Norm3(( this->GetPoint( 1 ) ) - ( this->GetPoint( 0 ) ) );
        return( rResult );
    }

    ///@}
    ///@name Edge
    ///@{

    /**
     * @brief This method gives you number of all edges of this geometry.
     * @details For example, for a hexahedron, this would be 12
     * @return SizeType containes number of this geometry edges.
     * @see EdgesNumber()
     * @see Edges()
     * @see GenerateEdges()
     * @see FacesNumber()
     * @see Faces()
     * @see GenerateFaces()
     */
    SizeType EdgesNumber() const override
    {
        return 1;
    }

    /**
     * @brief This method gives you all edges of this geometry.
     * @details This method will gives you all the edges with one dimension less than this geometry.
     * For example a triangle would return three lines as its edges or a tetrahedral would return four triangle as its edges but won't return its six edge lines by this method.
     * @return GeometriesArrayType containes this geometry edges.
     * @see EdgesNumber()
     * @see Edge()
     */
    GeometriesArrayType GenerateEdges() const override
    {
        GeometriesArrayType edges = GeometriesArrayType();
        edges.push_back( Kratos::make_shared<EdgeType>( this->pGetPoint( 0 ), this->pGetPoint( 1 ) ) );
        return edges;
    }

    ///@}
    ///@name Face
    ///@{

    /**
     * @brief Returns the number of faces of the current geometry.
     * @details This is only implemented for 3D geometries, since 2D geometries only have edges but no faces
     * @see EdgesNumber
     * @see Edges
     * @see Faces
     */
    SizeType FacesNumber() const override
    {
        return 0;
    }

    ///@}
    ///@name Shape Function
    ///@{

    /**
     * @brief This method gives all non-zero shape functions values evaluated at the rCoordinates provided
     * @note There is no control if the return vector is empty or not!
     * @return Vector of values of shape functions \f$ F_{i} \f$ where i is the shape function index (for NURBS it is the inde of the local enumeration in the element).
     * @see ShapeFunctionValue
     * @see ShapeFunctionsLocalGradients
     * @see ShapeFunctionLocalGradient
     */
    Vector& ShapeFunctionsValues (Vector &rResult, const CoordinatesArrayType& rCoordinates) const override
    {
        if(rResult.size() != 2) {
            rResult.resize(2, false);
        }

        rResult[0] =  0.5 * ( 1.0 - rCoordinates[0]);
        rResult[1] =  0.5 * ( 1.0 + rCoordinates[0]);

        return rResult;
    }

    /**
     * @brief This method gives value of given shape function evaluated in given point.
     * @param rPoint Point of evaluation of the shape function. This point must be in local coordinate.
     * @param ShapeFunctionIndex index of node which correspounding shape function evaluated in given integration point.
     * @return Value of given shape function in given point.
     * @see ShapeFunctionsValues
     * @see ShapeFunctionsLocalGradients
     * @see ShapeFunctionLocalGradient
     */
    double ShapeFunctionValue( IndexType ShapeFunctionIndex,
                                       const CoordinatesArrayType& rPoint ) const override
    {
        switch ( ShapeFunctionIndex )
        {
        case 0:
            return( 0.5*( 1.0 - rPoint[0] ) );

        case 1:
            return( 0.5*( 1.0 + rPoint[0] ) );

        default:
            KRATOS_ERROR << "Wrong index of shape function!" << *this << std::endl;
        }

        return 0;
    }

    Matrix& ShapeFunctionsLocalGradients( Matrix& rResult,
            const CoordinatesArrayType& rPoint ) const override
    {
        rResult = ZeroMatrix( 2, 1 );
        rResult( 0, 0 ) = -0.5;
        rResult( 1, 0 ) =  0.5;
        return( rResult );
    }

    /**
     * It computes the unit normal of the geometry, if possible
     * @return The normal of the geometry
     */
    array_1d<double, 3> Normal(const CoordinatesArrayType& rPointLocalCoordinates) const override
    {
        KRATOS_ERROR << "ERROR: Line3D2 can not define a normal. Please, define the normal in your implementation" << std::endl;
        return ZeroVector(3);
    }

    /**
     * @brief Returns whether given arbitrary point is inside the Geometry and the respective
     * local point for the given global point
     * @param rPoint The point to be checked if is inside o note in global coordinates
     * @param rResult The local coordinates of the point
     * @param Tolerance The  tolerance that will be considered to check if the point is inside or not
     * @return True if the point is inside, false otherwise
     */
    bool IsInside(
        const CoordinatesArrayType& rPoint,
        CoordinatesArrayType& rResult,
        const double Tolerance = std::numeric_limits<double>::epsilon()
        ) const override
    {
        PointLocalCoordinates( rResult, rPoint );

        if ( std::abs( rResult[0] ) <= (1.0 + Tolerance) ) {
            return true;
        }

        return false;
    }

    /** Test intersection of the geometry with a box (AABB)
     * Tests the intersection of the geometry with
     * a 3D box defined by rLowPoint and rHighPoint
     *
     * @param  rLowPoint  Lower point of the box to test the intersection
     * @param  rHighPoint Higher point of the box to test the intersection
     * @return            True if the geometry intersects the box, False in any other case.
     */
    bool HasIntersection(const Point& rLowPoint, const Point& rHighPoint) const override
    {
        return IntersectionUtilities::ComputeLineBoxIntersection(
            rLowPoint, rHighPoint, this->GetPoint(0), this->GetPoint(1));
    }

    /**
     * @brief Returns the local coordinates of a given arbitrary point
     * @param rResult The vector containing the local coordinates of the point
     * @param rPoint The point in global coordinates
     * @return The vector containing the local coordinates of the point
     */
    CoordinatesArrayType& PointLocalCoordinates(
        CoordinatesArrayType& rResult,
        const CoordinatesArrayType& rPoint
        ) const override
    {
        rResult.clear();

        const TPointType& r_first_point  = BaseType::GetPoint(0);
        const TPointType& r_second_point = BaseType::GetPoint(1);

        // Project point
        const double tolerance = 1e-14; // Tolerance

        const double length = Length();

        const double length_1 = std::sqrt( std::pow(rPoint[0] - r_first_point[0], 2)
                    + std::pow(rPoint[1] - r_first_point[1], 2) + std::pow(rPoint[2] - r_first_point[2], 2));

        const double length_2 = std::sqrt( std::pow(rPoint[0] - r_second_point[0], 2)
                    + std::pow(rPoint[1] - r_second_point[1], 2) + std::pow(rPoint[2] - r_second_point[2], 2));

        if (length_1 <= (length + tolerance) && length_2 <= (length + tolerance)) {
            rResult[0] = 2.0 * length_1/(length + tolerance) - 1.0;
        } else if (length_1 > (length + tolerance)) {
            rResult[0] = 2.0 * length_1/(length + tolerance) - 1.0; // NOTE: The same value as before, but it will be > than 1
        } else if (length_2 > (length + tolerance)) {
            rResult[0] = 1.0 - 2.0 * length_2/(length + tolerance);
        } else {
            rResult[0] = 2.0; // Out of the line!!!
        }

        return rResult ;
    }

    ///@}
    ///@name Input and output
    ///@{

    /** Turn back information as a string.

    @return String contains information about this geometry.
    @see PrintData()
    @see PrintInfo()
    */
    std::string Info() const override
    {
        return "1 dimensional line with 2 nodes in 3D space";
    }

    /** Print information about this object.

    @param rOStream Stream to print into it.
    @see PrintData()
    @see Info()
    */
    void PrintInfo( std::ostream& rOStream ) const override
    {
        rOStream << "1 dimensional line with 2 nodes in 3D space";
    }

    /** Print geometry's data into given stream. Prints it's points
    by the order they stored in the geometry and then center
    point of geometry.

    @param rOStream Stream to print into it.
    @see PrintInfo()
    @see Info()
    */
    void PrintData( std::ostream& rOStream ) const override
    {
        // Base Geometry class PrintData call
        BaseType::PrintData( rOStream );
        std::cout << std::endl;

        // If the geometry has valid points, calculate and output its data
        if (this->AllPointsAreValid()) {
            Matrix jacobian;
            this->Jacobian( jacobian, PointType() );
            rOStream << "    Jacobian\t : " << jacobian;
        }
    }

    ///@}
    ///@name Friends
    ///@{


    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{


    ///@}
    ///@name Protected member Variables
    ///@{


    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{


    ///@}
    ///@name Protected  Access
    ///@{


    ///@}
    ///@name Protected Inquiry
    ///@{


    ///@}
    ///@name Protected LifeCycle
    ///@{


    ///@}

private:
    ///@name Static Member Variables
    ///@{

    static const GeometryData msGeometryData;

    static const GeometryDimension msGeometryDimension;

    ///@}
    ///@name Member Variables
    ///@{

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

    Line3D2(): BaseType( PointsArrayType(), &msGeometryData ) {}


    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    static Matrix CalculateShapeFunctionsIntegrationPointsValues( typename BaseType::IntegrationMethod ThisMethod )
    {
        const IntegrationPointsContainerType& all_integration_points = AllIntegrationPoints();
        const IntegrationPointsArrayType& IntegrationPoints = all_integration_points[static_cast<int>(ThisMethod)];
        int integration_points_number = IntegrationPoints.size();
        Matrix N( integration_points_number, 2 );

        for ( int it_gp = 0; it_gp < integration_points_number; it_gp++ )
        {
            double e = IntegrationPoints[it_gp].X();
            N( it_gp, 0 ) = 0.5 * ( 1 - e );
            N( it_gp, 1 ) = 0.5 * ( 1 + e );
        }

        return N;
    }

    static ShapeFunctionsGradientsType CalculateShapeFunctionsIntegrationPointsLocalGradients( typename BaseType::IntegrationMethod ThisMethod )
    {
        const IntegrationPointsContainerType& all_integration_points = AllIntegrationPoints();
        const IntegrationPointsArrayType& IntegrationPoints = all_integration_points[static_cast<int>(ThisMethod)];
        ShapeFunctionsGradientsType DN_De( IntegrationPoints.size() );

        for ( unsigned int it_gp = 0; it_gp < IntegrationPoints.size(); it_gp++ )
        {
            Matrix aux_mat = ZeroMatrix(2,1);
            aux_mat(0,0) = -0.5;
            aux_mat(1,0) =  0.5;
            DN_De[it_gp] = aux_mat;
        }

        return DN_De;

    }

    static const IntegrationPointsContainerType AllIntegrationPoints()
    {
        IntegrationPointsContainerType integration_points = {{
                Quadrature<LineGaussLegendreIntegrationPoints1, 1, IntegrationPoint<3> >::GenerateIntegrationPoints(),
                Quadrature<LineGaussLegendreIntegrationPoints2, 1, IntegrationPoint<3> >::GenerateIntegrationPoints(),
                Quadrature<LineGaussLegendreIntegrationPoints3, 1, IntegrationPoint<3> >::GenerateIntegrationPoints(),
                Quadrature<LineGaussLegendreIntegrationPoints4, 1, IntegrationPoint<3> >::GenerateIntegrationPoints(),
                Quadrature<LineGaussLegendreIntegrationPoints5, 1, IntegrationPoint<3> >::GenerateIntegrationPoints(),
                Quadrature<LineCollocationIntegrationPoints1, 1, IntegrationPoint<3> >::GenerateIntegrationPoints(),
                Quadrature<LineCollocationIntegrationPoints2, 1, IntegrationPoint<3> >::GenerateIntegrationPoints(),
                Quadrature<LineCollocationIntegrationPoints3, 1, IntegrationPoint<3> >::GenerateIntegrationPoints(),
                Quadrature<LineCollocationIntegrationPoints4, 1, IntegrationPoint<3> >::GenerateIntegrationPoints(),
                Quadrature<LineCollocationIntegrationPoints5, 1, IntegrationPoint<3> >::GenerateIntegrationPoints()
            }
        };
        return integration_points;
    }

    static const ShapeFunctionsValuesContainerType AllShapeFunctionsValues()
    {
        ShapeFunctionsValuesContainerType shape_functions_values = {{
                Line3D2<TPointType>::CalculateShapeFunctionsIntegrationPointsValues( GeometryData::IntegrationMethod::GI_GAUSS_1 ),
                Line3D2<TPointType>::CalculateShapeFunctionsIntegrationPointsValues( GeometryData::IntegrationMethod::GI_GAUSS_2 ),
                Line3D2<TPointType>::CalculateShapeFunctionsIntegrationPointsValues( GeometryData::IntegrationMethod::GI_GAUSS_3 ),
                Line3D2<TPointType>::CalculateShapeFunctionsIntegrationPointsValues( GeometryData::IntegrationMethod::GI_GAUSS_4 ),
                Line3D2<TPointType>::CalculateShapeFunctionsIntegrationPointsValues( GeometryData::IntegrationMethod::GI_GAUSS_5 ),
                Line3D2<TPointType>::CalculateShapeFunctionsIntegrationPointsValues( GeometryData::IntegrationMethod::GI_EXTENDED_GAUSS_1 ),
                Line3D2<TPointType>::CalculateShapeFunctionsIntegrationPointsValues( GeometryData::IntegrationMethod::GI_EXTENDED_GAUSS_2 ),
                Line3D2<TPointType>::CalculateShapeFunctionsIntegrationPointsValues( GeometryData::IntegrationMethod::GI_EXTENDED_GAUSS_3 ),
                Line3D2<TPointType>::CalculateShapeFunctionsIntegrationPointsValues( GeometryData::IntegrationMethod::GI_EXTENDED_GAUSS_4 ),
                Line3D2<TPointType>::CalculateShapeFunctionsIntegrationPointsValues( GeometryData::IntegrationMethod::GI_EXTENDED_GAUSS_5 )
            }
        };
        return shape_functions_values;
    }

    static const ShapeFunctionsLocalGradientsContainerType AllShapeFunctionsLocalGradients()
    {
        ShapeFunctionsLocalGradientsContainerType shape_functions_local_gradients = {{
                Line3D2<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients( GeometryData::IntegrationMethod::GI_GAUSS_1 ),
                Line3D2<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients( GeometryData::IntegrationMethod::GI_GAUSS_2 ),
                Line3D2<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients( GeometryData::IntegrationMethod::GI_GAUSS_3 ),
                Line3D2<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients( GeometryData::IntegrationMethod::GI_GAUSS_4 ),
                Line3D2<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients( GeometryData::IntegrationMethod::GI_GAUSS_5 ),
                Line3D2<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients( GeometryData::IntegrationMethod::GI_EXTENDED_GAUSS_1 ),
                Line3D2<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients( GeometryData::IntegrationMethod::GI_EXTENDED_GAUSS_2 ),
                Line3D2<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients( GeometryData::IntegrationMethod::GI_EXTENDED_GAUSS_3 ),
                Line3D2<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients( GeometryData::IntegrationMethod::GI_EXTENDED_GAUSS_4 ),
                Line3D2<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients( GeometryData::IntegrationMethod::GI_EXTENDED_GAUSS_5 )

            }
        };
        return shape_functions_local_gradients;
    }

    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Private Friends
    ///@{

    template<class TOtherPointType> friend class Line3D2;

    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}

}; // Class Geometry

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template<class TPointType>
inline std::istream& operator >> ( std::istream& rIStream,
                                   Line3D2<TPointType>& rThis );

/// output stream function
template<class TPointType>
inline std::ostream& operator << ( std::ostream& rOStream,
                                   const Line3D2<TPointType>& rThis )
{
    rThis.PrintInfo( rOStream );
    rOStream << std::endl;
    rThis.PrintData( rOStream );

    return rOStream;
}

///@}


template<class TPointType>
const GeometryData Line3D2<TPointType>::msGeometryData(
        &msGeometryDimension,
        GeometryData::IntegrationMethod::GI_GAUSS_1,
        Line3D2<TPointType>::AllIntegrationPoints(),
        Line3D2<TPointType>::AllShapeFunctionsValues(),
        AllShapeFunctionsLocalGradients() );

template<class TPointType>
const GeometryDimension Line3D2<TPointType>::msGeometryDimension(3, 1);

}  // namespace Kratos.