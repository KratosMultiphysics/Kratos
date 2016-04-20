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


#if !defined(KRATOS_QUADRILATERAL_2D_9_H_INCLUDED )
#define  KRATOS_QUADRILATERAL_2D_9_H_INCLUDED

// System includes

// External includes

// Project includes
#include "geometries/line_2d_3.h"
#include "integration/quadrilateral_gauss_legendre_integration_points.h"


namespace Kratos
{
/**
 * A nine node quadrilateral geometry. While the shape functions are only defined in
 * 2D it is possible to define an arbitrary orientation in space. Thus it can be used for
 * defining surfaces on 3D elements.
 */

template<class TPointType> class Quadrilateral2D9 : public Geometry<TPointType>
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
    typedef Line2D3<TPointType> EdgeType;

    /**
     * Pointer definition of Quadrilateral2D9
     */
    KRATOS_CLASS_POINTER_DEFINITION( Quadrilateral2D9 );

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
     * Array of coordinates. Can be Nodes, Points or IntegrationPoints
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
     * This type is used to return size or dimension in
     * geometry. Dimension, WorkingDimension, PointsNumber and
     * ... return this type as their results.
     */
    typedef typename BaseType::SizeType SizeType;

    /**
     * Array of counted pointers to point.
     * This type is used to hold geometry's points.
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
     * method. IntegrationPoints functions used this type to return
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

    /**
     * Life Cycle
     */

    Quadrilateral2D9( const PointType& Point1, const PointType& Point2,
                      const PointType& Point3, const PointType& Point4,
                      const PointType& Point5, const PointType& Point6,
                      const PointType& Point7, const PointType& Point8,
                      const PointType& Point9 )
        : BaseType( PointsArrayType(), &msGeometryData )
    {
        this->Points().push_back( typename PointType::Pointer( new PointType( Point1 ) ) );
        this->Points().push_back( typename PointType::Pointer( new PointType( Point2 ) ) );
        this->Points().push_back( typename PointType::Pointer( new PointType( Point3 ) ) );
        this->Points().push_back( typename PointType::Pointer( new PointType( Point4 ) ) );
        this->Points().push_back( typename PointType::Pointer( new PointType( Point5 ) ) );
        this->Points().push_back( typename PointType::Pointer( new PointType( Point6 ) ) );
        this->Points().push_back( typename PointType::Pointer( new PointType( Point7 ) ) );
        this->Points().push_back( typename PointType::Pointer( new PointType( Point8 ) ) );
        this->Points().push_back( typename PointType::Pointer( new PointType( Point9 ) ) );
    }

    Quadrilateral2D9( typename PointType::Pointer pPoint1,
                      typename PointType::Pointer pPoint2,
                      typename PointType::Pointer pPoint3,
                      typename PointType::Pointer pPoint4,
                      typename PointType::Pointer pPoint5,
                      typename PointType::Pointer pPoint6,
                      typename PointType::Pointer pPoint7,
                      typename PointType::Pointer pPoint8,
                      typename PointType::Pointer pPoint9 )
        : BaseType( PointsArrayType(), &msGeometryData )
    {
        this->Points().push_back( pPoint1 );
        this->Points().push_back( pPoint2 );
        this->Points().push_back( pPoint3 );
        this->Points().push_back( pPoint4 );
        this->Points().push_back( pPoint5 );
        this->Points().push_back( pPoint6 );
        this->Points().push_back( pPoint7 );
        this->Points().push_back( pPoint8 );
        this->Points().push_back( pPoint9 );
    }

    Quadrilateral2D9( const PointsArrayType& ThisPoints )
        : BaseType( ThisPoints, &msGeometryData )
    {
        if ( this->PointsNumber() != 9 )
            KRATOS_THROW_ERROR( std::invalid_argument,
                                "Invalid points number. Expected 9, given " , this->PointsNumber() );
    }

    /**
     * Copy constructor.
     * Constructs this geometry as a copy of given geometry.
     *
     * @note This copy constructor don't copy the points and new
     * geometry shares points with given source geometry. It's
     * obvious that any change to this new geometry's point affect
     * source geometry's points too.
     */
    Quadrilateral2D9( Quadrilateral2D9 const& rOther )
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
     * @note This copy constructor don't copy the points and new
     * geometry shares points with given source geometry. It's
     * obvious that any change to this new geometry's point affect
     * source geometry's points too.
     */
    template<class TOtherPointType> Quadrilateral2D9(
        Quadrilateral2D9<TOtherPointType> const& rOther )
        : BaseType( rOther )
    {
    }

    /**
     * Destructor. Does nothing
     */
    virtual ~Quadrilateral2D9() {}

    GeometryData::KratosGeometryFamily GetGeometryFamily()
    {
        return GeometryData::Kratos_Quadrilateral;
    }

    GeometryData::KratosGeometryType GetGeometryType()
    {
        return GeometryData::Kratos_Quadrilateral2D9;
    }

    /**
     * Operators
     */

    /**
     * Assignment operator.
     *
     * @note This operator don't copy the points and this
     * geometry shares points with given source geometry. It's
     * obvious that any change to this geometry's point affect
     * source geometry's points too.
     *
     * @see Clone
     * @see ClonePoints
     */
    Quadrilateral2D9& operator=( const Quadrilateral2D9& rOther )
    {
        BaseType::operator=( rOther );

        return *this;
    }

    /**
     * Assignment operator for geometries with different point type.
     *
     * @note This operator don't copy the points and this
     * geometry shares points with given source geometry. It's
     * obvious that any change to this geometry's point affect
     * source geometry's points too.
     *
     * @see Clone
     * @see ClonePoints
     */
    template<class TOtherPointType>
    Quadrilateral2D9& operator=( Quadrilateral2D9<TOtherPointType> const & rOther )
    {
        BaseType::operator=( rOther );

        return *this;
    }

    /**
     * Operations
     */
    typename BaseType::Pointer Create( PointsArrayType const& ThisPoints ) const
    {
        return typename BaseType::Pointer( new Quadrilateral2D9( ThisPoints ) );
    }

    virtual Geometry< Point<3> >::Pointer Clone() const
    {
        Geometry< Point<3> >::PointsArrayType NewPoints;

        //making a copy of the nodes TO POINTS (not Nodes!!!)
        for ( IndexType i = 0 ; i < this->size() ; i++ )
        {
                NewPoints.push_back(boost::make_shared< Point<3> >(( *this )[i]));
        }

        //creating a geometry with the new points
        Geometry< Point<3> >::Pointer p_clone( new Quadrilateral2D9< Point<3> >( NewPoints ) );

        return p_clone;
    }

    /**
     * lumping factors for the calculation of the lumped mass matrix
     */
    virtual Vector& LumpingFactors( Vector& rResult ) const
    {
        if(rResult.size() != 9)
            rResult.resize( 9, false );

        for ( int i = 0; i < 4; i++ ) rResult[i] = 1.00 / 36.00;

        for ( int i = 4; i < 8; i++ ) rResult[i] = 1.00 / 9.00;

        rResult[8] = 4.00 / 9.00;

        return rResult;
    }

    /**
     * Informations
     */

    /**
     * :TODO: the charactereistic sizes have to be reviewed
     * by the one who is willing to use them!
     */
    /** This method calculates and returns Length or charactereistic
     * length of this geometry depending on it's dimension. For one
     * dimensional geometry for example Line it returns length of it
     * and for the other geometries it gives Characteristic length
     * otherwise.
     *
     * @return double value contains length or Characteristic
     * length
     * @see Area()
     * @see Volume()
     * @see DomainSize()
     */
    virtual double Length() const
    {
        return sqrt( fabs( this->DeterminantOfJacobian( PointType() ) ) );
    }

    /**
     * :TODO: the charactereistic sizes have to be reviewed
     * by the one who is willing to use them!
     */
    /**
     * This method calculates and returns area or surface area of
     * this geometry depending to it's dimension. For one dimensional
     * geometry it returns zero, for two dimensional it gives area
     * and for three dimensional geometries it gives surface area.
     * @return double value contains area or surface
     * area.
     * @see Length()
     * @see Volume()
     * @see DomainSize()
     */
    virtual double Area() const
    {

        Vector temp;
        this->DeterminantOfJacobian( temp, msGeometryData.DefaultIntegrationMethod() );
        const IntegrationPointsArrayType& integration_points = this->IntegrationPoints( msGeometryData.DefaultIntegrationMethod() );
        double Area = 0.00;

        for ( unsigned int i = 0; i < integration_points.size(); i++ )
        {
            Area += temp[i] * integration_points[i].Weight();
        }

        //KRATOS_WATCH(temp)
        return Area;

        //return fabs(DeterminantOfJacobian(PointType())) * 0.5;
    }

    /**
     * :TODO: the charactereistic sizes have to be reviewed
     * by the one who is willing to use them!
     */
    /**
     * This method calculate and return length, area or volume of
     * this geometry depending to it's dimension. For one dimensional
     * geometry it returns its length, for two dimensional it gives area
     * and for three dimensional geometries it gives its volume.
     *
     * @return double value contains length, area or volume.
     * @see Length()
     * @see Area()
     * @see Volume()
     */
    virtual double DomainSize() const
    {
        return fabs( this->DeterminantOfJacobian( PointType() ) ) * 0.5;
    }

    /**
     * Returns whether given arbitrary point is inside the Geometry
     */
    virtual bool IsInside( const CoordinatesArrayType& rPoint, CoordinatesArrayType& rResult )
    {
        this->PointLocalCoordinates( rResult, rPoint );

        if ( rResult[0] >= -1.0 && rResult[0] <= 1.0 )
            if ( rResult[1] >= -1.0 && rResult[1] <= 1.0 )
                return true;

        return false;
    }

    /** This method gives you number of all edges of this
    geometry. This method will gives you number of all the edges
    with one dimension less than this geometry. for example a
    triangle would return three or a tetrahedral would return
    four but won't return nine related to its six edge lines.

    @return SizeType containes number of this geometry edges.
    @see Edges()
    @see Edge()
     */
    virtual SizeType EdgesNumber() const
    {
        return 4;
    }

    /** This method gives you all edges of this geometry. This
     * method will gives you all the edges with one dimension less
     * than this geometry. for example a triangle would return
     * three lines as its edges or a tetrahedral would return four
     * triangle as its edges but won't return its six edge
     * lines by this method.
     *
     * @return GeometriesArrayType containes this geometry edges.
     * @see EdgesNumber()
     * @see Edge()
     */
    virtual GeometriesArrayType Edges( void )
    {
        GeometriesArrayType edges = GeometriesArrayType();
        edges.push_back( boost::make_shared<EdgeType>( this->pGetPoint( 0 ), this->pGetPoint( 4 ), this->pGetPoint( 1 ) ) );
        edges.push_back( boost::make_shared<EdgeType>( this->pGetPoint( 1 ), this->pGetPoint( 5 ), this->pGetPoint( 2 ) ) );
        edges.push_back( boost::make_shared<EdgeType>( this->pGetPoint( 2 ), this->pGetPoint( 6 ), this->pGetPoint( 3 ) ) );
        edges.push_back( boost::make_shared<EdgeType>( this->pGetPoint( 3 ), this->pGetPoint( 7 ), this->pGetPoint( 0 ) ) );
        return edges;
    }

    /**
     * Shape Function
     */

    /**
     * :TODO: implemented but not yet tested
     */
    /**
     * Calculates the value of a given shape function at a given point.
     *
     * @param ShapeFunctionIndex The number of the desired shape function
     * @param rPoint the given point in local coordinates at which the
     * value of the shape function is calculated
     *
     * @return the value of the shape function at the given point
     */
    virtual double ShapeFunctionValue( IndexType ShapeFunctionIndex,
                                       const CoordinatesArrayType& rPoint ) const
    {
        double fx1 = 0.5 * ( rPoint[0] - 1 ) * rPoint[0];
        double fx2 = 0.5 * ( rPoint[0] + 1 ) * rPoint[0];
        double fx3 = 1 - rPoint[0] * rPoint[0];
        double fy1 = 0.5 * ( rPoint[1] - 1 ) * rPoint[1];
        double fy2 = 0.5 * ( rPoint[1] + 1 ) * rPoint[1];
        double fy3 = 1 - rPoint[1] * rPoint[1];

        switch ( ShapeFunctionIndex )
        {
        case 0:
            return( fx1*fy1 );
        case 1:
            return( fx2*fy1 );
        case 2:
            return( fx2*fy2 );
        case 3:
            return( fx1*fy2 );
        case 4:
            return( fx3*fy1 );
        case 5:
            return( fx2*fy3 );
        case 6:
            return( fx3*fy2 );
        case 7:
            return( fx1*fy3 );
        case 8:
            return( fx3*fy3 );
        default:
            KRATOS_THROW_ERROR( std::logic_error,
                                "Wrong index of shape function!" , *this );
        }

        return 0;
    }
    
        /** This method gives all non-zero shape functions values
    evaluated at the rCoordinates provided

    @return Vector of values of shape functions \f$ F_{i} \f$
    where i is the shape function index (for NURBS it is the index
    of the local enumeration in the element).

    @see ShapeFunctionValue
    @see ShapeFunctionsLocalGradients
    @see ShapeFunctionLocalGradient
    */
    virtual Vector& ShapeFunctionsValues (Vector &rResult, const CoordinatesArrayType& rCoordinates) const
    {
        if(rResult.size() != 9) rResult.resize(9,false);

        double fx1 = 0.5 * ( rCoordinates[0] - 1 ) * rCoordinates[0];
        double fx2 = 0.5 * ( rCoordinates[0] + 1 ) * rCoordinates[0];
        double fx3 = 1 - rCoordinates[0] * rCoordinates[0];
        double fy1 = 0.5 * ( rCoordinates[1] - 1 ) * rCoordinates[1];
        double fy2 = 0.5 * ( rCoordinates[1] + 1 ) * rCoordinates[1];
        double fy3 = 1 - rCoordinates[1] * rCoordinates[1];

        rResult[0] =   fx1*fy1 ;
        rResult[1] =   fx2*fy1 ;
        rResult[2] =   fx2*fy2 ;
        rResult[3] =   fx1*fy2 ;
        rResult[4] =   fx3*fy1 ;
        rResult[5] =   fx2*fy3 ;
        rResult[6] =   fx3*fy2 ;
        rResult[7] =   fx1*fy3 ;
        rResult[8] =   fx3*fy3 ;

        return rResult;
    }


    /**
     * Input and Output
     */


    /**
     * Turn back information as a string.
     * @return String contains information about this geometry.
     * @see PrintData()
     * @see PrintInfo()
     */
    virtual std::string Info() const
    {
        return "2 dimensional quadrilateral with nine nodes in 2D space";
    }

    /**
     * Print information about this object.
     * @param rOStream Stream to print into it.
     * @see PrintData()
     * @see Info()
     */
    virtual void PrintInfo( std::ostream& rOStream ) const
    {
        rOStream << "2 dimensional quadrilateral with nine nodes in 2D space";
    }

    /**
     * Print geometry's data into given stream. Prints it's points
     * by the order they stored in the geometry and then center
     * point of geometry.
     *
     * @param rOStream Stream to print into it.
     * @see PrintInfo()
     * @see Info()
     */
    virtual void PrintData( std::ostream& rOStream ) const
    {
        BaseType::PrintData( rOStream );
        std::cout << std::endl;
        Matrix jacobian;
        this->Jacobian( jacobian, PointType() );
        rOStream << "    Jacobian in the origin\t : " << jacobian;
    }

    /**
     * Calculates the local gradients for all integration points for
     * given integration method
     */
    virtual ShapeFunctionsGradientsType ShapeFunctionsLocalGradients(
        IntegrationMethod ThisMethod )
    {
        ShapeFunctionsGradientsType localGradients
            = CalculateShapeFunctionsIntegrationPointsLocalGradients( ThisMethod );
        const int integration_points_number
            = msGeometryData.IntegrationPointsNumber( ThisMethod );
        ShapeFunctionsGradientsType Result( integration_points_number );

        for ( int pnt = 0; pnt < integration_points_number; pnt++ )
        {
            Result[pnt] = localGradients[pnt];
        }

        return Result;
    }

    /**
     * Calculates the local gradients for all integration points for the
     * default integration method
     */
    virtual ShapeFunctionsGradientsType ShapeFunctionsLocalGradients()
    {
        IntegrationMethod ThisMethod = msGeometryData.DefaultIntegrationMethod();
        ShapeFunctionsGradientsType localGradients
            = CalculateShapeFunctionsIntegrationPointsLocalGradients( ThisMethod );
        const int integration_points_number
            = msGeometryData.IntegrationPointsNumber( ThisMethod );
        ShapeFunctionsGradientsType Result( integration_points_number );

        for ( int pnt = 0; pnt < integration_points_number; pnt++ )
        {
            Result[pnt] = localGradients[pnt];
        }

        return Result;
    }

    /**
     * Calculates the gradients in terms of local coordinateds of
     * all shape functions in a given point.
     *
     * @param rPoint the current point at which the gradients are calculated
     * @return the gradients of all shape functions
     * \f$ \frac{\partial N^i}{\partial \xi_j} \f$
     */
    virtual Matrix& ShapeFunctionsLocalGradients( Matrix& rResult,
            const CoordinatesArrayType& rPoint ) const
    {
        double fx1 = 0.5 * ( rPoint[0] - 1 ) * rPoint[0];
        double fx2 = 0.5 * ( rPoint[0] + 1 ) * rPoint[0];
        double fx3 = 1 - rPoint[0] * rPoint[0];
        double fy1 = 0.5 * ( rPoint[1] - 1 ) * rPoint[1];
        double fy2 = 0.5 * ( rPoint[1] + 1 ) * rPoint[1];
        double fy3 = 1 - rPoint[1] * rPoint[1];

        double gx1 = 0.5 * ( 2 * rPoint[0] - 1 );
        double gx2 = 0.5 * ( 2 * rPoint[0] + 1 );
        double gx3 = -2.0 * rPoint[0];
        double gy1 = 0.5 * ( 2 * rPoint[1] - 1 );
        double gy2 = 0.5 * ( 2 * rPoint[1] + 1 );
        double gy3 = -2.0 * rPoint[1];

        rResult.resize( 9, 2, false );
        noalias( rResult ) = ZeroMatrix( 9, 2 );
        rResult( 0, 0 ) = gx1 * fy1;
        rResult( 0, 1 ) = fx1 * gy1;
        rResult( 1, 0 ) = gx2 * fy1;
        rResult( 1, 1 ) = fx2 * gy1;
        rResult( 2, 0 ) = gx2 * fy2;
        rResult( 2, 1 ) = fx2 * gy2;
        rResult( 3, 0 ) = gx1 * fy2;
        rResult( 3, 1 ) = fx1 * gy2;
        rResult( 4, 0 ) = gx3 * fy1;
        rResult( 4, 1 ) = fx3 * gy1;
        rResult( 5, 0 ) = gx2 * fy3;
        rResult( 5, 1 ) = fx2 * gy3;
        rResult( 6, 0 ) = gx3 * fy2;
        rResult( 6, 1 ) = fx3 * gy2;
        rResult( 7, 0 ) = gx1 * fy3;
        rResult( 7, 1 ) = fx1 * gy3;
        rResult( 8, 0 ) = gx3 * fy3;
        rResult( 8, 1 ) = fx3 * gy3;

        return rResult;
    }

    /**
     * returns the local coordinates of all nodes of the current geometry
     * @param rResult a Matrix object that will be overwritten by the result
     * @return the local coordinates of all nodes
     */
    virtual Matrix& PointsLocalCoordinates( Matrix& rResult ) const
    {
        rResult.resize( 9, 2, false );
        noalias( rResult ) = ZeroMatrix( 9, 2 );
        rResult( 0, 0 ) = -1.0;
        rResult( 0, 1 ) = -1.0;
        rResult( 1, 0 ) =  1.0;
        rResult( 1, 1 ) = -1.0;
        rResult( 2, 0 ) =  1.0;
        rResult( 2, 1 ) =  1.0;
        rResult( 3, 0 ) = -1.0;
        rResult( 3, 1 ) =  1.0;
        rResult( 4, 0 ) =  0.0;
        rResult( 4, 1 ) = -1.0;
        rResult( 5, 0 ) =  1.0;
        rResult( 5, 1 ) =  0.0;
        rResult( 6, 0 ) =  0.0;
        rResult( 6, 1 ) =  1.0;
        rResult( 7, 0 ) = -1.0;
        rResult( 7, 1 ) =  0.0;
        rResult( 8, 0 ) =  0.0;
        rResult( 8, 1 ) =  0.0;
        return rResult;
    }

    /**
     * returns the shape function gradients in an arbitrary point,
     * given in local coordinates
     *
     * @param rResult the matrix of gradients, will be overwritten
     * with the gradients for all shape functions in given point
     * @param rPoint the given point the gradients are calculated in
     */
    virtual Matrix& ShapeFunctionsGradients( Matrix& rResult, PointType& rPoint )
    {
        double fx1 = 0.5 * ( rPoint.X() - 1 ) * rPoint.X();
        double fx2 = 0.5 * ( rPoint.X() + 1 ) * rPoint.X();
        double fx3 = 1 - rPoint.X() * rPoint.X();
        double fy1 = 0.5 * ( rPoint.Y() - 1 ) * rPoint.Y();
        double fy2 = 0.5 * ( rPoint.Y() + 1 ) * rPoint.Y();
        double fy3 = 1 - rPoint.Y() * rPoint.Y();

        double gx1 = 0.5 * ( 2 * rPoint.X() - 1 );
        double gx2 = 0.5 * ( 2 * rPoint.X() + 1 );
        double gx3 = -2.0 * rPoint.X();
        double gy1 = 0.5 * ( 2 * rPoint.Y() - 1 );
        double gy2 = 0.5 * ( 2 * rPoint.Y() + 1 );
        double gy3 = -2.0 * rPoint.Y();

        rResult.resize( 9, 2, false );
        noalias( rResult ) = ZeroMatrix( 9, 2 );
        rResult( 0, 0 ) = gx1 * fy1;
        rResult( 0, 1 ) = fx1 * gy1;
        rResult( 1, 0 ) = gx2 * fy1;
        rResult( 1, 1 ) = fx2 * gy1;
        rResult( 2, 0 ) = gx2 * fy2;
        rResult( 2, 1 ) = fx2 * gy2;
        rResult( 3, 0 ) = gx1 * fy2;
        rResult( 3, 1 ) = fx1 * gy2;
        rResult( 4, 0 ) = gx3 * fy1;
        rResult( 4, 1 ) = fx3 * gy1;
        rResult( 5, 0 ) = gx2 * fy3;
        rResult( 5, 1 ) = fx2 * gy3;
        rResult( 6, 0 ) = gx3 * fy2;
        rResult( 6, 1 ) = fx3 * gy2;
        rResult( 7, 0 ) = gx1 * fy3;
        rResult( 7, 1 ) = fx1 * gy3;
        rResult( 8, 0 ) = gx3 * fy3;
        rResult( 8, 1 ) = fx3 * gy3;

        return rResult;
    }

    /**
     * returns the second order derivatives of all shape functions
     * in given arbitrary points
     * @param rResult a third order tensor which contains the second derivatives
     * @param rPoint the given point the second order derivatives are calculated in
     */
    virtual ShapeFunctionsSecondDerivativesType& ShapeFunctionsSecondDerivatives( ShapeFunctionsSecondDerivativesType& rResult, const CoordinatesArrayType& rPoint ) const
    {
        if ( rResult.size() != this->PointsNumber() )
        {
            // KLUDGE: While there is a bug in ublas vector resize, I have to put this beside resizing!!
            ShapeFunctionsGradientsType temp( this->PointsNumber() );
            rResult.swap( temp );
        }

        for ( unsigned int i = 0; i < this->PointsNumber(); i++ )
        {
            rResult[i].resize( 2, 2, false );
            noalias( rResult[i] ) = ZeroMatrix( 2, 2 );
        }

        double fx1 = 0.5 * ( rPoint[0] - 1 ) * rPoint[0];
        double fx2 = 0.5 * ( rPoint[0] + 1 ) * rPoint[0];
        double fx3 = 1 - rPoint[0] * rPoint[0];
        double fy1 = 0.5 * ( rPoint[1] - 1 ) * rPoint[1];
        double fy2 = 0.5 * ( rPoint[1] + 1 ) * rPoint[1];
        double fy3 = 1 - rPoint[1] * rPoint[1];

        double gx1 = 0.5 * ( 2 * rPoint[0] - 1 );
        double gx2 = 0.5 * ( 2 * rPoint[0] + 1 );
        double gx3 = -2.0 * rPoint[0];
        double gy1 = 0.5 * ( 2 * rPoint[1] - 1 );
        double gy2 = 0.5 * ( 2 * rPoint[1] + 1 );
        double gy3 = -2.0 * rPoint[1];

        double hx1 = 1.0;
        double hx2 = 1.0;
        double hx3 = -2.0;
        double hy1 = 1.0;
        double hy2 = 1.0;
        double hy3 = -2.0;

        rResult[0]( 0, 0 ) = hx1 * fy1;
        rResult[0]( 0, 1 ) = gx1 * gy1;
        rResult[0]( 1, 0 ) = gx1 * gy1;
        rResult[0]( 1, 1 ) = fx1 * hy1;

        rResult[1]( 0, 0 ) = hx2 * fy1;
        rResult[1]( 0, 1 ) = gx2 * gy1;
        rResult[1]( 1, 0 ) = gx2 * gy1;
        rResult[1]( 1, 1 ) = fx2 * hy1;

        rResult[2]( 0, 0 ) = hx2 * fy2;
        rResult[2]( 0, 1 ) = gx2 * gy2;
        rResult[2]( 1, 0 ) = gx2 * gy2;
        rResult[2]( 1, 1 ) = fx2 * hy2;

        rResult[3]( 0, 0 ) = hx1 * fy2;
        rResult[3]( 0, 1 ) = gx1 * gy2;
        rResult[3]( 1, 0 ) = gx1 * gy2;
        rResult[3]( 1, 1 ) = fx1 * hy2;

        rResult[4]( 0, 0 ) = hx3 * fy1;
        rResult[4]( 0, 1 ) = gx3 * gy1;
        rResult[4]( 1, 0 ) = gx3 * gy1;
        rResult[4]( 1, 1 ) = fx3 * hy1;

        rResult[5]( 0, 0 ) = hx2 * fy3;
        rResult[5]( 0, 1 ) = gx2 * gy3;
        rResult[5]( 1, 0 ) = gx2 * gy3;
        rResult[5]( 1, 1 ) = fx2 * hy3;

        rResult[6]( 0, 0 ) = hx3 * fy2;
        rResult[6]( 0, 1 ) = gx3 * gy2;
        rResult[6]( 1, 0 ) = gx3 * gy2;
        rResult[6]( 1, 1 ) = fx3 * hy2;

        rResult[7]( 0, 0 ) = hx1 * fy3;
        rResult[7]( 0, 1 ) = gx1 * gy3;
        rResult[7]( 1, 0 ) = gx1 * gy3;
        rResult[7]( 1, 1 ) = fx1 * hy3;

        rResult[8]( 0, 0 ) = hx3 * fy3;
        rResult[8]( 0, 1 ) = gx3 * gy3;
        rResult[8]( 1, 0 ) = gx3 * gy3;
        rResult[8]( 1, 1 ) = fx3 * hy3;

        return rResult;
    }


    /**
    * returns the third order derivatives of all shape functions
    * in given arbitrary points
    * @param rResult a fourth order tensor which contains the third derivatives
    * @param rPoint the given point the third order derivatives are calculated in
     */
    virtual ShapeFunctionsThirdDerivativesType& ShapeFunctionsThirdDerivatives( ShapeFunctionsThirdDerivativesType& rResult, const CoordinatesArrayType& rPoint ) const
    {

        if ( rResult.size() != this->PointsNumber() )
        {
            // KLUDGE: While there is a bug in
            // ublas vector resize, I have to put this beside resizing!!
//                 ShapeFunctionsGradientsType
            ShapeFunctionsThirdDerivativesType temp( this->PointsNumber() );
            rResult.swap( temp );
        }

        for ( IndexType i = 0; i < rResult.size(); i++ )
        {
            boost::numeric::ublas::vector<Matrix> temp( this->PointsNumber() );
            rResult[i].swap( temp );
        }

        for ( unsigned int i = 0; i < this->PointsNumber(); i++ )
        {
            for ( unsigned int j = 0; j < 2; j++ )
            {
                rResult[i][j].resize( 2, 2, false );
                noalias( rResult[i][j] ) = ZeroMatrix( 2, 2 );
            }
        }

//        double fx1 = 0.5 * ( rPoint[0] - 1 ) * rPoint[0];

//        double fx2 = 0.5 * ( rPoint[0] + 1 ) * rPoint[0];
//        double fx3 = 1 - rPoint[0] * rPoint[0];
//        double fy1 = 0.5 * ( rPoint[1] - 1 ) * rPoint[1];
//        double fy2 = 0.5 * ( rPoint[1] + 1 ) * rPoint[1];
//        double fy3 = 1 - rPoint[1] * rPoint[1];

        double gx1 = 0.5 * ( 2 * rPoint[0] - 1 );
        double gx2 = 0.5 * ( 2 * rPoint[0] + 1 );
        double gx3 = -2.0 * rPoint[0];
        double gy1 = 0.5 * ( 2 * rPoint[1] - 1 );
        double gy2 = 0.5 * ( 2 * rPoint[1] + 1 );
        double gy3 = -2.0 * rPoint[1];

        double hx1 = 1.0;
        double hx2 = 1.0;
        double hx3 = -2.0;
        double hy1 = 1.0;
        double hy2 = 1.0;
        double hy3 = -2.0;

        rResult[0][0]( 0, 0 ) = 0.0;
        rResult[0][0]( 0, 1 ) = hx1 * gy1;
        rResult[0][0]( 1, 0 ) = hx1 * gy1;
        rResult[0][0]( 1, 1 ) = gx1 * hy1;
        rResult[0][1]( 0, 0 ) = hx1 * gy1;
        rResult[0][1]( 0, 1 ) = gx1 * hy1;
        rResult[0][1]( 1, 0 ) = gx1 * hy1;
        rResult[0][1]( 1, 1 ) = 0.0;

        rResult[1][0]( 0, 0 ) = 0.0;
        rResult[1][0]( 0, 1 ) = hx2 * gy1;
        rResult[1][0]( 1, 0 ) = hx2 * gy1;
        rResult[1][0]( 1, 1 ) = gx2 * hy1;
        rResult[1][1]( 0, 0 ) = hx2 * gy1;
        rResult[1][1]( 0, 1 ) = gx2 * hy1;
        rResult[1][1]( 1, 0 ) = gx2 * hy1;
        rResult[1][1]( 1, 1 ) = 0.0;

        rResult[2][0]( 0, 0 ) = 0.0;
        rResult[2][0]( 0, 1 ) = hx2 * gy2;
        rResult[2][0]( 1, 0 ) = hx2 * gy2;
        rResult[2][0]( 1, 1 ) = gx2 * hy2;
        rResult[2][1]( 0, 0 ) = hx2 * gy2;
        rResult[2][1]( 0, 1 ) = gx2 * hy2;
        rResult[2][1]( 1, 0 ) = gx2 * hy2;
        rResult[2][1]( 1, 1 ) = 0.0;

        rResult[3][0]( 0, 0 ) = 0.0;
        rResult[3][0]( 0, 1 ) = hx1 * gy2;
        rResult[3][0]( 1, 0 ) = hx1 * gy2;
        rResult[3][0]( 1, 1 ) = gx1 * hy2;
        rResult[3][1]( 0, 0 ) = hx1 * gy2;
        rResult[3][1]( 0, 1 ) = gx1 * hy2;
        rResult[3][1]( 1, 0 ) = gx1 * hy2;
        rResult[3][1]( 1, 1 ) = 0.0;

        rResult[4][0]( 0, 0 ) = 0.0;
        rResult[4][0]( 0, 1 ) = hx3 * gy1;
        rResult[4][0]( 1, 0 ) = hx3 * gy1;
        rResult[4][0]( 1, 1 ) = gx3 * hy1;
        rResult[4][1]( 0, 0 ) = hx3 * gy1;
        rResult[4][1]( 0, 1 ) = gx3 * hy1;
        rResult[4][1]( 1, 0 ) = gx3 * hy1;
        rResult[4][1]( 1, 1 ) = 0.0;

        rResult[5][0]( 0, 0 ) = 0.0;
        rResult[5][0]( 0, 1 ) = hx2 * gy3;
        rResult[5][0]( 1, 0 ) = hx2 * gy3;
        rResult[5][0]( 1, 1 ) = gx2 * hy3;
        rResult[5][1]( 0, 0 ) = hx2 * gy3;
        rResult[5][1]( 0, 1 ) = gx2 * hy3;
        rResult[5][1]( 1, 0 ) = gx2 * hy3;
        rResult[5][1]( 1, 1 ) = 0.0;

        rResult[6][0]( 0, 0 ) = 0.0;
        rResult[6][0]( 0, 1 ) = hx3 * gy2;
        rResult[6][0]( 1, 0 ) = hx3 * gy2;
        rResult[6][0]( 1, 1 ) = gx3 * hy2;
        rResult[6][1]( 0, 0 ) = hx3 * gy2;
        rResult[6][1]( 0, 1 ) = gx3 * hy2;
        rResult[6][1]( 1, 0 ) = gx3 * hy2;
        rResult[6][1]( 1, 1 ) = 0.0;

        rResult[7][0]( 0, 0 ) = 0.0;
        rResult[7][0]( 0, 1 ) = hx1 * gy3;
        rResult[7][0]( 1, 0 ) = hx1 * gy3;
        rResult[7][0]( 1, 1 ) = gx1 * hy3;
        rResult[7][1]( 0, 0 ) = hx1 * gy3;
        rResult[7][1]( 0, 1 ) = gx1 * hy3;
        rResult[7][1]( 1, 0 ) = gx1 * hy3;
        rResult[7][1]( 1, 1 ) = 0.0;

        rResult[8][0]( 0, 0 ) = 0.0;
        rResult[8][0]( 0, 1 ) = hx3 * gy3;
        rResult[8][0]( 1, 0 ) = hx3 * gy3;
        rResult[8][0]( 1, 1 ) = gx3 * hy3;
        rResult[8][1]( 0, 0 ) = hx3 * gy3;
        rResult[8][1]( 0, 1 ) = gx3 * hy3;
        rResult[8][1]( 1, 0 ) = gx3 * hy3;
        rResult[8][1]( 1, 1 ) = 0.0;


        return rResult;
    }

protected:


    /**
     * there are no protected class members
     */


private:

    /**
     * Static Member Variables
     */

    static const GeometryData msGeometryData;

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    virtual void save( Serializer& rSerializer ) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, PointsArrayType );
    }

    virtual void load( Serializer& rSerializer )
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, PointsArrayType );
    }

    Quadrilateral2D9(): BaseType( PointsArrayType(), &msGeometryData ) {}

    /**
     * Private Operations
     */




    /**
     * :TODO: implemented but not yet tested
     */
    /**
     * Calculates the values of all shape function in all integration points.
     * Integration points are expected to be given in local coordinates
     *
     * @param ThisMethod the current integration method
     * @return the matrix of values of every shape function in each integration point
     *
     */
    static Matrix CalculateShapeFunctionsIntegrationPointsValues(
        typename BaseType::IntegrationMethod ThisMethod )
    {
        IntegrationPointsContainerType all_integration_points = AllIntegrationPoints();
        IntegrationPointsArrayType integration_points = all_integration_points[ThisMethod];
        //number of integration points
        const int integration_points_number = integration_points.size();
        //number of nodes in current geometry
        const int points_number = 9;
        //setting up return matrix
        Matrix shape_function_values( integration_points_number, points_number );
        //loop over all integration points

        for ( int pnt = 0; pnt < integration_points_number; pnt++ )
        {
            double fx1 = 0.5 * ( integration_points[pnt].X() - 1 ) * integration_points[pnt].X();
            double fx2 = 0.5 * ( integration_points[pnt].X() + 1 ) * integration_points[pnt].X();
            double fx3 = 1 - integration_points[pnt].X() * integration_points[pnt].X();
            double fy1 = 0.5 * ( integration_points[pnt].Y() - 1 ) * integration_points[pnt].Y();
            double fy2 = 0.5 * ( integration_points[pnt].Y() + 1 ) * integration_points[pnt].Y();
            double fy3 = 1 - integration_points[pnt].Y() * integration_points[pnt].Y();

            shape_function_values( pnt, 0 ) = ( fx1 * fy1 );
            shape_function_values( pnt, 1 ) = ( fx2 * fy1 );
            shape_function_values( pnt, 2 ) = ( fx2 * fy2 );
            shape_function_values( pnt, 3 ) = ( fx1 * fy2 );
            shape_function_values( pnt, 4 ) = ( fx3 * fy1 );
            shape_function_values( pnt, 5 ) = ( fx2 * fy3 );
            shape_function_values( pnt, 6 ) = ( fx3 * fy2 );
            shape_function_values( pnt, 7 ) = ( fx1 * fy3 );
            shape_function_values( pnt, 8 ) = ( fx3 * fy3 );
        }

        return shape_function_values;
    }

    /**
     * :TODO: implemented but not yet tested
     */
    /**
     * Calculates the local gradients of all shape functions in
     * all integration points.
     * Integration points are expected to be given in local coordinates
     *
     * @param ThisMethod the current integration method
     * @return the vector of the gradients of all shape functions
     * in each integration point
     */
    static ShapeFunctionsGradientsType CalculateShapeFunctionsIntegrationPointsLocalGradients(
        typename BaseType::IntegrationMethod ThisMethod )
    {
        IntegrationPointsContainerType all_integration_points = AllIntegrationPoints();
        IntegrationPointsArrayType integration_points = all_integration_points[ThisMethod];
        //number of integration points
        const int integration_points_number = integration_points.size();
        ShapeFunctionsGradientsType d_shape_f_values( integration_points_number );
        //initialising container
        //std::fill(d_shape_f_values.begin(), d_shape_f_values.end(), Matrix(4,2));
        //loop over all integration points

        for ( int pnt = 0; pnt < integration_points_number; pnt++ )
        {
            double fx1 = 0.5 * ( integration_points[pnt].X() - 1 ) * integration_points[pnt].X();
            double fx2 = 0.5 * ( integration_points[pnt].X() + 1 ) * integration_points[pnt].X();
            double fx3 = 1 - integration_points[pnt].X() * integration_points[pnt].X();
            double fy1 = 0.5 * ( integration_points[pnt].Y() - 1 ) * integration_points[pnt].Y();
            double fy2 = 0.5 * ( integration_points[pnt].Y() + 1 ) * integration_points[pnt].Y();
            double fy3 = 1 - integration_points[pnt].Y() * integration_points[pnt].Y();

            double gx1 = 0.5 * ( 2 * integration_points[pnt].X() - 1 );
            double gx2 = 0.5 * ( 2 * integration_points[pnt].X() + 1 );
            double gx3 = -2.0 * integration_points[pnt].X();
            double gy1 = 0.5 * ( 2 * integration_points[pnt].Y() - 1 );
            double gy2 = 0.5 * ( 2 * integration_points[pnt].Y() + 1 );
            double gy3 = -2.0 * integration_points[pnt].Y();

            Matrix result( 9, 2 );
            result( 0, 0 ) = gx1 * fy1;
            result( 0, 1 ) = fx1 * gy1;
            result( 1, 0 ) = gx2 * fy1;
            result( 1, 1 ) = fx2 * gy1;
            result( 2, 0 ) = gx2 * fy2;
            result( 2, 1 ) = fx2 * gy2;
            result( 3, 0 ) = gx1 * fy2;
            result( 3, 1 ) = fx1 * gy2;
            result( 4, 0 ) = gx3 * fy1;
            result( 4, 1 ) = fx3 * gy1;
            result( 5, 0 ) = gx2 * fy3;
            result( 5, 1 ) = fx2 * gy3;
            result( 6, 0 ) = gx3 * fy2;
            result( 6, 1 ) = fx3 * gy2;
            result( 7, 0 ) = gx1 * fy3;
            result( 7, 1 ) = fx1 * gy3;
            result( 8, 0 ) = gx3 * fy3;
            result( 8, 1 ) = fx3 * gy3;

            d_shape_f_values[pnt] = result;
        }

        return d_shape_f_values;
    }

    /**
     * :TODO: testing
     */
    static const IntegrationPointsContainerType AllIntegrationPoints()
    {
        IntegrationPointsContainerType integration_points =
        {
            {
                Quadrature < QuadrilateralGaussLegendreIntegrationPoints1,
                2, IntegrationPoint<3> >::GenerateIntegrationPoints(),
                Quadrature < QuadrilateralGaussLegendreIntegrationPoints2,
                2, IntegrationPoint<3> >::GenerateIntegrationPoints(),
                Quadrature < QuadrilateralGaussLegendreIntegrationPoints3,
                2, IntegrationPoint<3> >::GenerateIntegrationPoints(),
                Quadrature < QuadrilateralGaussLegendreIntegrationPoints4,
                2, IntegrationPoint<3> >::GenerateIntegrationPoints()
            }
        };
        return integration_points;
    }

    /**
     * :TODO: testing
     */
    static const ShapeFunctionsValuesContainerType AllShapeFunctionsValues()
    {
        ShapeFunctionsValuesContainerType shape_functions_values =
        {
            {
                Quadrilateral2D9<TPointType>::CalculateShapeFunctionsIntegrationPointsValues(
                    GeometryData::GI_GAUSS_1 ),
                Quadrilateral2D9<TPointType>::CalculateShapeFunctionsIntegrationPointsValues(
                    GeometryData::GI_GAUSS_2 ),
                Quadrilateral2D9<TPointType>::CalculateShapeFunctionsIntegrationPointsValues(
                    GeometryData::GI_GAUSS_3 ),
                Quadrilateral2D9<TPointType>::CalculateShapeFunctionsIntegrationPointsValues(
                    GeometryData::GI_GAUSS_4 ),
            }
        };
        return shape_functions_values;
    }

    /**
     * :TODO: testing
     */
    static const ShapeFunctionsLocalGradientsContainerType AllShapeFunctionsLocalGradients()
    {
        ShapeFunctionsLocalGradientsContainerType shape_functions_local_gradients =
        {
            {
                Quadrilateral2D9<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients
                ( GeometryData::GI_GAUSS_1 ),
                Quadrilateral2D9<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients
                ( GeometryData::GI_GAUSS_2 ),
                Quadrilateral2D9<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients
                ( GeometryData::GI_GAUSS_3 ),
                Quadrilateral2D9<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients
                ( GeometryData::GI_GAUSS_4 ),
            }
        };
        return shape_functions_local_gradients;
    }

    /**
     * Private Friends
     */

    template<class TOtherPointType> friend class Quadrilateral2D9;

    /**
     * Un accessible methods
     */

}; // Class Quadrilateral2D9

/**
 * Input and output
 */

/**
 * input stream function
 */
template< class TPointType > inline std::istream& operator >> (
    std::istream& rIStream,
    Quadrilateral2D9<TPointType>& rThis );

/**
 * output stream function
 */
template< class TPointType > inline std::ostream& operator << (
    std::ostream& rOStream,
    const Quadrilateral2D9<TPointType>& rThis )
{
    rThis.PrintInfo( rOStream );
    rOStream << std::endl;
    rThis.PrintData( rOStream );
    return rOStream;
}

template<class TPointType>
const GeometryData Quadrilateral2D9<TPointType>::msGeometryData(
    2, 2, 2,
    GeometryData::GI_GAUSS_3,
    Quadrilateral2D9<TPointType>::AllIntegrationPoints(),
    Quadrilateral2D9<TPointType>::AllShapeFunctionsValues(),
    AllShapeFunctionsLocalGradients()
);

}  // namespace Kratos.

#endif // KRATOS_QUADRILATERAL_2D_9_H_INCLUDED  defined 
