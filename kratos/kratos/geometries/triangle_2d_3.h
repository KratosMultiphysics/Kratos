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



#if !defined(KRATOS_TRIANGLE_2D_3_H_INCLUDED )
#define  KRATOS_TRIANGLE_2D_3_H_INCLUDED


// System includes

// External includes

// Project includes
#include "geometries/line_2d_2.h"
#include "integration/triangle_gauss_legendre_integration_points.h"
#include "integration/triangle_collocation_integration_points.h"

//#include  "utilities/triangle_triangle_intersection.h"


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
 * A three node element geometry. While the shape functions are only defined in
 * 2D it is possible to define an arbitrary orientation in space. Thus it can be used for
 * defining surfaces on 3D elements.
 */

template<class TPointType> class Triangle2D3
    : public Geometry<TPointType>
{
public:
    ///@}
    ///@name Type Definitions
    ///@{

    /**
     * Geometry as base class.
     */
    typedef Geometry<TPointType> BaseType;

    /**
     * Type of edge geometry
     */
    typedef Line2D2<TPointType> EdgeType;

    /**
     * Pointer definition of Triangle2D3
     */
    KRATOS_CLASS_POINTER_DEFINITION( Triangle2D3 );

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
     * This type used to hold geometry's points.
     */
    typedef  typename BaseType::PointsArrayType PointsArrayType;

    /**
     * Array of coordinates. Can be Nodes, Points or IntegrationPoints
     */
    typedef typename BaseType::CoordinatesArrayType CoordinatesArrayType;

    /**
     * This type used for representing an integration point in geometry.
     * This integration point is a point with an additional weight component.
     */
    typedef typename BaseType::IntegrationPointType IntegrationPointType;

    /**
     * A Vector of IntegrationPointType which used to hold
     * integration points related to an integration
     * method.
     * IntegrationPoints functions used this type to return
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
    typedef typename BaseType::ShapeFunctionsValuesContainerType ShapeFunctionsValuesContainerType;

    /**
     * A fourth order tensor used as shape functions' local
     * gradients container in geometry.
     */
    typedef typename BaseType::ShapeFunctionsLocalGradientsContainerType ShapeFunctionsLocalGradientsContainerType;

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
    * A third order tensor to hold shape functions' local third derivatives.
    * ShapefunctionsLocalGradients function return this
    * type as its result.
    */
    typedef typename BaseType::ShapeFunctionsThirdDerivativesType
    ShapeFunctionsThirdDerivativesType;

    /**
     * Type of the normal vector used for normal to edges in geometry.
     */
    typedef typename BaseType::NormalType NormalType;

    ///@}
    ///@name Life Cycle
    ///@{

//     Triangle2D3( const PointType& FirstPoint,
//                  const PointType& SecondPoint,
//                  const PointType& ThirdPoint )
//         : BaseType( PointsArrayType(), &msGeometryData )
//     {
//         this->Points().push_back( typename PointType::Pointer( new PointType( FirstPoint ) ) );
//         this->Points().push_back( typename PointType::Pointer( new PointType( SecondPoint ) ) );
//         this->Points().push_back( typename PointType::Pointer( new PointType( ThirdPoint ) ) );
//     }

    Triangle2D3( typename PointType::Pointer pFirstPoint,
                 typename PointType::Pointer pSecondPoint,
                 typename PointType::Pointer pThirdPoint )
        : BaseType( PointsArrayType(), &msGeometryData )
    {
        this->Points().push_back( pFirstPoint );
        this->Points().push_back( pSecondPoint );
        this->Points().push_back( pThirdPoint );
    }

    Triangle2D3( const PointsArrayType& ThisPoints )
        : BaseType( ThisPoints, &msGeometryData )
    {
        if ( this->PointsNumber() != 3 )
            KRATOS_THROW_ERROR( std::invalid_argument,
                                "Invalid points number. Expected 3, given " , this->PointsNumber() );
    }

    /**
     * Copy constructor.
     * Construct this geometry as a copy of given geometry.
     *
     * @note This copy constructor does not copy the points and new
     * geometry shares points with given source geometry. It is
     * obvious that any change to this new geometry's point affect
     * source geometry's points too.
     */
    Triangle2D3( Triangle2D3 const& rOther )
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
     * @note This copy constructor does not copy the points and new
     * geometry shares points with given source geometry. It is
     * obvious that any change to this new geometry's point affect
     * source geometry's points too.
     */
    template<class TOtherPointType> Triangle2D3( Triangle2D3<TOtherPointType> const& rOther )
        : BaseType( rOther )
    {
    }

    /**
     * Destructor. Does nothing!!!
     */
    virtual ~Triangle2D3() {}

    GeometryData::KratosGeometryFamily GetGeometryFamily()
    {
        return GeometryData::Kratos_Triangle;
    }

    GeometryData::KratosGeometryType GetGeometryType()
    {
        return GeometryData::Kratos_Triangle2D3;
    }

    ///@}
    ///@name Operators
    ///@{

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
    Triangle2D3& operator=( const Triangle2D3& rOther )
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
    Triangle2D3& operator=( Triangle2D3<TOtherPointType> const & rOther )
    {
        BaseType::operator=( rOther );
        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    typename BaseType::Pointer Create( PointsArrayType const& ThisPoints ) const
    {
        return typename BaseType::Pointer( new Triangle2D3( ThisPoints ) );
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
        Geometry< Point<3> >::Pointer p_clone( new Triangle2D3< Point<3> >( NewPoints ) );

        return p_clone;
    }

    /**
     * returns the local coordinates of all nodes of the current geometry
     * @param rResult a Matrix object that will be overwritten by the result
     * @return the local coordinates of all nodes
     */
    virtual Matrix& PointsLocalCoordinates( Matrix& rResult ) const
    {
        rResult = ZeroMatrix( 3, 2 );
        rResult( 0, 0 ) =  0.0;
        rResult( 0, 1 ) =  0.0;
        rResult( 1, 0 ) =  1.0;
        rResult( 1, 1 ) =  0.0;
        rResult( 2, 0 ) =  0.0;
        rResult( 2, 1 ) =  1.0;
        return rResult;
    }

    //lumping factors for the calculation of the lumped mass matrix
    virtual Vector& LumpingFactors( Vector& rResult ) const
    {
        if(rResult.size() != 3)
            rResult.resize( 3, false );
        std::fill( rResult.begin(), rResult.end(), 1.00 / 3.00 );
        return rResult;
    }

    ///@}
    ///@name Information
    ///@{

    /**
     * This method calculates and returns Length or charactereistic
     * length of this geometry depending on it's dimension.
     * For one dimensional geometry for example Line it returns
     * length of it and for the other geometries it gives Characteristic
     * length otherwise.
     * In the current geometry this function returns the determinant of
     * jacobian
     *
     * @return double value contains length or Characteristic
     * length
     * @see Area()
     * @see Volume()
     * @see DomainSize()
     */
    /**
     * :TODO: could be replaced by something more suitable
     * (comment by janosch)
     */
    virtual double Length() const
    {
        //return sqrt(fabs( DeterminantOfJacobian(PointType()))*0.5);
        double length = 0.000;
        length = 1.1283791670955 * sqrt( fabs( Area() ) );  // (4xA/Pi)^(1/2)
        return length;
    }


    /** This method calculates and returns area or surface area of
     * this geometry depending to it's dimension. For one dimensional
     * geometry it returns zero, for two dimensional it gives area
     * and for three dimensional geometries it gives surface area.
     *
     * @return double value contains area or surface
     * area.
     * @see Length()
     * @see Volume()
     * @see DomainSize()
     */
    /**
     * :TODO: could be replaced by something more suitable
     * (comment by janosch)
     */
    virtual double Area() const
    {
        const PointType& p0 = this->operator [](0);
        const PointType& p1 = this->operator [](1);
        const PointType& p2 = this->operator [](2);

        double x10 = p1.X() - p0.X();
        double y10 = p1.Y() - p0.Y();

        double x20 = p2.X() - p0.X();
        double y20 = p2.Y() - p0.Y();

        double detJ = x10 * y20-y10 * x20;
        return 0.5*detJ;
    }


    /// detect if two triangle are intersected
    virtual bool HasIntersection( const BaseType& rThisGeometry )
    {
        const BaseType& geom_1 = *this;
        const BaseType& geom_2 = rThisGeometry;
        return  NoDivTriTriIsect(geom_1[0].Coordinates(), geom_1[1].Coordinates() , geom_1[2].Coordinates(),
                                 geom_2[0].Coordinates(), geom_2[1].Coordinates(),  geom_2[2].Coordinates());
    }


    /// detect if  triangle and box are intersected
    virtual bool HasIntersection( const Point<3, double>& rLowPoint, const Point<3, double>& rHighPoint )
    {
        //return true;
        const BaseType& geom_1 = *this;
        //std::size_t dim        =  geom_1.WorkingSpaceDimension();

        Point<3, double> boxcenter;
        Point<3, double> boxhalfsize;

        boxcenter[0]   = 0.50 * (rLowPoint[0] + rHighPoint[0]);
        boxcenter[1]   = 0.50 * (rLowPoint[1] + rHighPoint[1]);
        boxcenter[2]   = 0.00;


        boxhalfsize[0] = 0.50 * (rHighPoint[0] - rLowPoint[0]);
        boxhalfsize[1] = 0.50 * (rHighPoint[1] - rLowPoint[1]);
        boxhalfsize[2] = 0.00;

        std::size_t size = geom_1.size();
        std::vector<Point<3, double> > triverts;
        triverts.resize(size);
        for(unsigned int i = 0; i< size; i++ )
            triverts[i] =  geom_1.GetPoint(i);

        bool result = false;
        result      = TriBoxOverlap(boxcenter, boxhalfsize, triverts);
        return result;
    }



    /** This method calculates and returns length, area or volume of
     * this geometry depending to it's dimension. For one dimensional
     * geometry it returns its length, for two dimensional it gives area
     * and for three dimensional geometries it gives its volume.
     *
     * @return double value contains length, area or volume.
     * @see Length()
     * @see Area()
     * @see Volume()
     */
    /**
     * :TODO: could be replaced by something more suitable
     * (comment by janosch)
     */
    virtual double DomainSize() const
    {
        return this->Area();
    }

    /**
     * Returns whether given arbitrary point is inside the Geometry
     */
    virtual bool IsInside( const CoordinatesArrayType& rPoint, CoordinatesArrayType& rResult )
    {
        const double zero = 1E-14;
        this->PointLocalCoordinates( rResult, rPoint );
        if( ( rResult[0] >= (0.0-zero) ) && ( rResult[0] <= 1.0 + zero ) )
            if( ( rResult[1] >= 0.0-zero ) && (rResult[1] <= 1.0 + zero ) )
                if(((1.0-(rResult[0] + rResult[1])) >= 0.0-zero) &&  ((1.0-(rResult[0] + rResult[1])) <= 1.0 + zero))
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
        return 3;
    }


    virtual SizeType FacesNumber() const
    {
        return 3;
    }


    /** This method gives you all edges of this geometry. This
    method will gives you all the edges with one dimension less
    than this geometry. for example a triangle would return
    three lines as its edges or a tetrahedral would return four
    triangle as its edges but won't return its six edge
    lines by this method.

    @return GeometriesArrayType containes this geometry edges.
    @see EdgesNumber()
    @see Edge()
     */
    virtual GeometriesArrayType Edges( void )
    {
        GeometriesArrayType edges = GeometriesArrayType();

        edges.push_back( boost::make_shared<EdgeType>( this->pGetPoint( 0 ), this->pGetPoint( 1 ) ) );
        edges.push_back( boost::make_shared<EdgeType>( this->pGetPoint( 1 ), this->pGetPoint( 2 ) ) );
        edges.push_back( boost::make_shared<EdgeType>( this->pGetPoint( 2 ), this->pGetPoint( 0 ) ) );
        return edges;
    }



    //Connectivities of faces required
    virtual void NumberNodesInFaces (boost::numeric::ublas::vector<unsigned int>& NumberNodesInFaces) const
    {
        if(NumberNodesInFaces.size() != 3 )
            NumberNodesInFaces.resize(3,false);
        // Linear Triangles have elements of 2 nodes as faces
        NumberNodesInFaces[0]=2;
        NumberNodesInFaces[1]=2;
        NumberNodesInFaces[2]=2;

    }

    virtual void NodesInFaces (boost::numeric::ublas::matrix<unsigned int>& NodesInFaces) const
    {
        if(NodesInFaces.size1() != 3 || NodesInFaces.size2() != 3)
            NodesInFaces.resize(3,3,false);

        NodesInFaces(0,0)=0;//face or other node
        NodesInFaces(1,0)=1;
        NodesInFaces(2,0)=2;

        NodesInFaces(0,1)=1;//face or other node 
        NodesInFaces(1,1)=2;
        NodesInFaces(2,1)=0;

        NodesInFaces(0,2)=2;//face or other node
        NodesInFaces(1,2)=0;
        NodesInFaces(2,2)=1;

    }



    ///@}
    ///@name Shape Function
    ///@{

    /**
     * TODO: implemented but not yet tested
     */
    /**
     * Calculates the value of a given shape function at a given point.
     *
     * @param ShapeFunctionIndex The number of the desired shape function
     * @param rPoint the given point in local coordinates at which the value of the shape
     * function is calculated
     *
     * @return the value of the shape function at the given point
     */
    virtual double ShapeFunctionValue( IndexType ShapeFunctionIndex,
                                       const CoordinatesArrayType& rPoint ) const
    {
        switch ( ShapeFunctionIndex )
        {
        case 0:
            return( 1.0 -rPoint[0] - rPoint[1] );
        case 1:
            return( rPoint[0] );
        case 2:
            return( rPoint[1] );
        default:
            KRATOS_THROW_ERROR( std::logic_error,
                                "Wrong index of shape function!" ,
                                *this );
        }

        return 0;
    }
    
        /** This method gives all non-zero shape functions values
    evaluated at the rCoordinates provided

    \note There is no control if the return vector is empty or not!

    @return Vector of values of shape functions \f$ F_{i} \f$
    where i is the shape function index (for NURBS it is the index
    of the local enumeration in the element).

    @see ShapeFunctionValue
    @see ShapeFunctionsLocalGradients
    @see ShapeFunctionLocalGradient
    */

    virtual Vector& ShapeFunctionsValues (Vector &rResult, const CoordinatesArrayType& rCoordinates) const
    {
      if(rResult.size() != 3) rResult.resize(3,false);
      rResult[0] =  1.0 -rCoordinates[0] - rCoordinates[1];
      rResult[1] =  rCoordinates[0] ;
      rResult[2] =  rCoordinates[1] ;
        
        return rResult;
    }


    /**
     * Calculates the Gradients of the shape functions.
     * Calculates the gradients of the shape functions with regard to
     * the global coordinates in all
     * integration points (\f$ \frac{\partial N^i}{\partial X_j} \f$)
     *
     * @param rResult a container which takes the calculated gradients
     * @param ThisMethod the given IntegrationMethod
     *
     * @return the gradients of all shape functions with regard to the global coordinates

    */
    virtual ShapeFunctionsGradientsType& ShapeFunctionsIntegrationPointsGradients(
        ShapeFunctionsGradientsType& rResult,
        IntegrationMethod ThisMethod ) const
    {
        const unsigned int integration_points_number =
            msGeometryData.IntegrationPointsNumber( ThisMethod );

        boost::numeric::ublas::bounded_matrix<double,3,2> DN_DX;
        double x10 = this->Points()[1].X() - this->Points()[0].X();
        double y10 = this->Points()[1].Y() - this->Points()[0].Y();

        double x20 = this->Points()[2].X() - this->Points()[0].X();
        double y20 = this->Points()[2].Y() - this->Points()[0].Y();

        //Jacobian is calculated:
        //  |dx/dxi  dx/deta|	|x1-x0   x2-x0|
        //J=|				|=	|			  |
        //  |dy/dxi  dy/deta|	|y1-y0   y2-y0|


        double detJ = x10 * y20-y10 * x20;

        DN_DX(0,0) = -y20 + y10;
        DN_DX(0,1) =  x20 - x10;
        DN_DX(1,0) =  y20;
        DN_DX(1,1) = -x20;
        DN_DX(2,0) = -y10;
        DN_DX(2,1) =  x10;

        DN_DX /= detJ;

        //workaround by riccardo
        if(rResult.size() != integration_points_number)
        {
            rResult.resize(integration_points_number,false);
        }
        for(unsigned int i=0; i<integration_points_number; i++)
            rResult[i] = DN_DX;

        return rResult;
    }

    virtual ShapeFunctionsGradientsType& ShapeFunctionsIntegrationPointsGradients(
        ShapeFunctionsGradientsType& rResult,
        Vector& determinants_of_jacobian,
        IntegrationMethod ThisMethod ) const
    {
        const unsigned int integration_points_number =
            msGeometryData.IntegrationPointsNumber( ThisMethod );

        boost::numeric::ublas::bounded_matrix<double,3,2> DN_DX;
        double x10 = this->Points()[1].X() - this->Points()[0].X();
        double y10 = this->Points()[1].Y() - this->Points()[0].Y();

        double x20 = this->Points()[2].X() - this->Points()[0].X();
        double y20 = this->Points()[2].Y() - this->Points()[0].Y();

        //Jacobian is calculated:
        //  |dx/dxi  dx/deta|	|x1-x0   x2-x0|
        //J=|				|=	|			  |
        //  |dy/dxi  dy/deta|	|y1-y0   y2-y0|


        double detJ = x10 * y20-y10 * x20;

        DN_DX(0,0) = -y20 + y10;
        DN_DX(0,1) =  x20 - x10;
        DN_DX(1,0) =  y20;
        DN_DX(1,1) = -x20;
        DN_DX(2,0) = -y10;
        DN_DX(2,1) =  x10;

        DN_DX /= detJ;

        //workaround by riccardo
        if(rResult.size() != integration_points_number)
        {
            rResult.resize(integration_points_number,false);
        }
        for(unsigned int i=0; i<integration_points_number; i++)
            rResult[i] = DN_DX;

        if(determinants_of_jacobian.size() != integration_points_number )
            determinants_of_jacobian.resize(integration_points_number,false);

        for(unsigned int i=0; i<integration_points_number; i++)
            determinants_of_jacobian[i] = detJ;

        return rResult;
    }


    ///@}
    ///@name Input and output
    ///@{

    /**
     * Turn back information as a string.
     *
     * @return String contains information about this geometry.
     * @see PrintData()
     * @see PrintInfo()
     */
    virtual std::string Info() const
    {
        return "2 dimensional triangle with three nodes in 2D space";
    }

    /**
     * Print information about this object.
     * @param rOStream Stream to print into it.
     * @see PrintData()
     * @see Info()
     */
    virtual void PrintInfo( std::ostream& rOStream ) const
    {
        rOStream << "2 dimensional triangle with three nodes in 2D space";
    }

    /**
     * Print geometry's data into given stream.
     * Prints it's points
     * by the order they stored in the geometry and then center
     * point of geometry.
     *
     * @param rOStream Stream to print into it.
     * @see PrintInfo()
     * @see Info()
     */
    /**
     * :TODO: needs to be reviewed because it is not properly implemented yet
     * (comment by janosch)
     */
    virtual void PrintData( std::ostream& rOStream ) const
    {
        PrintInfo( rOStream );
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
     * Calculates the gradients in terms of local coordinates
     * of all shape functions in a given point.
     *
     * @param rPoint the current point at which the gradients are calculated in local
     * coordinates
     * @return the gradients of all shape functions
     * \f$ \frac{\partial N^i}{\partial \xi_j} \f$
     */
    virtual Matrix& ShapeFunctionsLocalGradients( Matrix& rResult,
            const CoordinatesArrayType& rPoint ) const
    {
        rResult = ZeroMatrix( 3, 2 );
        rResult( 0, 0 ) = -1.0;
        rResult( 0, 1 ) = -1.0;
        rResult( 1, 0 ) =  1.0;
        rResult( 1, 1 ) =  0.0;
        rResult( 2, 0 ) =  0.0;
        rResult( 2, 1 ) =  1.0;
        return rResult;
    }



    /**
     * returns the shape function gradients in an arbitrary point,
     * given in local coordinates
     *
     * @param rResult the matrix of gradients,
     * will be overwritten with the gradients for all
     * shape functions in given point
     * @param rPoint the given point the gradients are calculated in
     */
    virtual Matrix& ShapeFunctionsGradients( Matrix& rResult, PointType& rPoint )
    {
        rResult = ZeroMatrix( 3, 2 );
        rResult( 0, 0 ) = -1.0;
        rResult( 0, 1 ) = -1.0;
        rResult( 1, 0 ) =  1.0;
        rResult( 1, 1 ) =  0.0;
        rResult( 2, 0 ) =  0.0;
        rResult( 2, 1 ) =  1.0;
        return rResult;
    }

    /**
     * returns the second order derivatives of all shape functions
     * in given arbitrary pointers
     * @param rResult a third order tensor which contains the second derivatives
     * @param rPoint the given point the second order derivatives are calculated in
     */
    virtual ShapeFunctionsSecondDerivativesType& ShapeFunctionsSecondDerivatives( ShapeFunctionsSecondDerivativesType& rResult, const CoordinatesArrayType& rPoint ) const
    {
        if ( rResult.size() != this->PointsNumber() )
        {
            // KLUDGE: While there is a bug in
            // ublas vector resize, I have to put this beside resizing!!
            ShapeFunctionsGradientsType temp( this->PointsNumber() );
            rResult.swap( temp );
        }

        if(rResult[0].size1() != 2 || rResult[0].size2() != 2 )
            rResult[0].resize( 2, 2,false );

        if(rResult[1].size1() != 2 || rResult[1].size2() != 2 )
            rResult[1].resize( 2, 2,false );

        if(rResult[2].size1() != 2 || rResult[2].size2() != 2 )
            rResult[2].resize( 2, 2,false );

        rResult[0]( 0, 0 ) = 0.0;
        rResult[0]( 0, 1 ) = 0.0;
        rResult[0]( 1, 0 ) = 0.0;
        rResult[0]( 1, 1 ) = 0.0;
        rResult[1]( 0, 0 ) = 0.0;
        rResult[1]( 0, 1 ) = 0.0;
        rResult[1]( 1, 0 ) = 0.0;
        rResult[1]( 1, 1 ) = 0.0;
        rResult[2]( 0, 0 ) = 0.0;
        rResult[2]( 0, 1 ) = 0.0;
        rResult[2]( 1, 0 ) = 0.0;
        rResult[2]( 1, 1 ) = 0.0;
        return rResult;
    }

    /**
     * returns the third order derivatives of all shape functions
     * in given arbitrary pointers
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

        rResult[0][0].resize( 2, 2,false );

        rResult[0][1].resize( 2, 2,false );
        rResult[1][0].resize( 2, 2,false );
        rResult[1][1].resize( 2, 2,false );
        rResult[2][0].resize( 2, 2,false );
        rResult[2][1].resize( 2, 2,false );

        for ( int i = 0; i < 3; i++ )
        {

            rResult[i][0]( 0, 0 ) = 0.0;
            rResult[i][0]( 0, 1 ) = 0.0;
            rResult[i][0]( 1, 0 ) = 0.0;
            rResult[i][0]( 1, 1 ) = 0.0;
            rResult[i][1]( 0, 0 ) = 0.0;
            rResult[i][1]( 0, 1 ) = 0.0;
            rResult[i][1]( 1, 0 ) = 0.0;
            rResult[i][1]( 1, 1 ) = 0.0;
        }

        return rResult;
    }

    ///@}
    ///@name Friends
    ///@{

    ///@}

protected:
    /**
     * There are no protected members in class Triangle2D3
     */

private:
    ///@name Static Member Variables
    ///@{
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

    Triangle2D3():BaseType( PointsArrayType(), &msGeometryData ) {}

    ///@}
    ///@name Member Variables
    ///@{


    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    /**
     * TODO: implemented but not yet tested
     */
    /**
     * Calculates the values of all shape function in all integration points.
     * Integration points are expected to be given in local coordinates
     * @param ThisMethod the current integration method
     * @return the matrix of values of every shape function in each integration point
     * :KLUDGE: number of points is hard-coded -> be careful if you want to copy and paste!
     */
    static Matrix CalculateShapeFunctionsIntegrationPointsValues(
        typename BaseType::IntegrationMethod ThisMethod )
    {
        IntegrationPointsContainerType all_integration_points =
            AllIntegrationPoints();
        IntegrationPointsArrayType integration_points =
            all_integration_points[ThisMethod];
        //number of integration points
        const int integration_points_number = integration_points.size();
        //number of nodes in current geometry
        const int points_number = 3;
        //setting up return matrix
        Matrix shape_function_values( integration_points_number, points_number );
        //loop over all integration points

        for ( int pnt = 0; pnt < integration_points_number; pnt++ )
        {
            shape_function_values( pnt, 0 ) = 1.0
                                              - integration_points[pnt].X()
                                              - integration_points[pnt].Y();
            shape_function_values( pnt, 1 ) = integration_points[pnt].X();
            shape_function_values( pnt, 2 ) = integration_points[pnt].Y();
        }

        return shape_function_values;
    }

    /**
     * TODO: implemented but not yet tested
     */
    /**
     * Calculates the local gradients of all shape functions
     * in all integration points.
     * Integration points are expected to be given in local coordinates
     *
     * @param ThisMethod the current integration method
     * @return the vector of the gradients of all shape functions
     * in each integration point
     */
    static ShapeFunctionsGradientsType
    CalculateShapeFunctionsIntegrationPointsLocalGradients(
        typename BaseType::IntegrationMethod ThisMethod )
    {
        IntegrationPointsContainerType all_integration_points =
            AllIntegrationPoints();
        IntegrationPointsArrayType integration_points =
            all_integration_points[ThisMethod];
        //number of integration points
        const int integration_points_number = integration_points.size();
        ShapeFunctionsGradientsType d_shape_f_values( integration_points_number );
        //initialising container
        //std::fill(d_shape_f_values.begin(), d_shape_f_values.end(), Matrix(4,2));
        //loop over all integration points

        for ( int pnt = 0; pnt < integration_points_number; pnt++ )
        {
            Matrix result( 3, 2 );
            result( 0, 0 ) = -1.0;
            result( 0, 1 ) = -1.0;
            result( 1, 0 ) =  1.0;
            result( 1, 1 ) =  0.0;
            result( 2, 0 ) =  0.0;
            result( 2, 1 ) =  1.0;
            d_shape_f_values[pnt] = result;
        }

        return d_shape_f_values;
    }

    /**
     * TODO: testing
     */
    static const IntegrationPointsContainerType AllIntegrationPoints()
    {
        IntegrationPointsContainerType integration_points =
        {
            {
                Quadrature<TriangleGaussLegendreIntegrationPoints1, 2, IntegrationPoint<3> >::GenerateIntegrationPoints(),
                Quadrature<TriangleGaussLegendreIntegrationPoints2, 2, IntegrationPoint<3> >::GenerateIntegrationPoints(),
                Quadrature<TriangleGaussLegendreIntegrationPoints3, 2, IntegrationPoint<3> >::GenerateIntegrationPoints(),
                Quadrature<TriangleGaussLegendreIntegrationPoints4, 2, IntegrationPoint<3> >::GenerateIntegrationPoints(),
                Quadrature<TriangleGaussLegendreIntegrationPoints5, 2, IntegrationPoint<3> >::GenerateIntegrationPoints(),
                Quadrature<TriangleCollocationIntegrationPoints1, 2, IntegrationPoint<3> >::GenerateIntegrationPoints(),
                Quadrature<TriangleCollocationIntegrationPoints2, 2, IntegrationPoint<3> >::GenerateIntegrationPoints(),
                Quadrature<TriangleCollocationIntegrationPoints3, 2, IntegrationPoint<3> >::GenerateIntegrationPoints(),
                Quadrature<TriangleCollocationIntegrationPoints4, 2, IntegrationPoint<3> >::GenerateIntegrationPoints(),
                Quadrature<TriangleCollocationIntegrationPoints5, 2, IntegrationPoint<3> >::GenerateIntegrationPoints()
            }
        };
        return integration_points;
    }

    /**
     * TODO: testing
     */
    static const ShapeFunctionsValuesContainerType AllShapeFunctionsValues()
    {
        ShapeFunctionsValuesContainerType shape_functions_values =
        {
            {
                Triangle2D3<TPointType>::CalculateShapeFunctionsIntegrationPointsValues(
                    GeometryData::GI_GAUSS_1 ),
                Triangle2D3<TPointType>::CalculateShapeFunctionsIntegrationPointsValues(
                    GeometryData::GI_GAUSS_2 ),
                Triangle2D3<TPointType>::CalculateShapeFunctionsIntegrationPointsValues(
                    GeometryData::GI_GAUSS_3 ),
                Triangle2D3<TPointType>::CalculateShapeFunctionsIntegrationPointsValues(
                    GeometryData::GI_GAUSS_4 ),
                Triangle2D3<TPointType>::CalculateShapeFunctionsIntegrationPointsValues(
                    GeometryData::GI_GAUSS_5 ),
                Triangle2D3<TPointType>::CalculateShapeFunctionsIntegrationPointsValues(
                    GeometryData::GI_EXTENDED_GAUSS_1 ),
                Triangle2D3<TPointType>::CalculateShapeFunctionsIntegrationPointsValues(
                    GeometryData::GI_EXTENDED_GAUSS_2 ),
                Triangle2D3<TPointType>::CalculateShapeFunctionsIntegrationPointsValues(
                    GeometryData::GI_EXTENDED_GAUSS_3 ),
                Triangle2D3<TPointType>::CalculateShapeFunctionsIntegrationPointsValues(
                    GeometryData::GI_EXTENDED_GAUSS_4 ),
                Triangle2D3<TPointType>::CalculateShapeFunctionsIntegrationPointsValues(
                    GeometryData::GI_EXTENDED_GAUSS_5 ),
            }
        };
        return shape_functions_values;
    }

    /**
     * TODO: testing
     */
    static const ShapeFunctionsLocalGradientsContainerType
    AllShapeFunctionsLocalGradients()
    {
        ShapeFunctionsLocalGradientsContainerType shape_functions_local_gradients =
        {
            {
                Triangle2D3<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients( GeometryData::GI_GAUSS_1 ),
                Triangle2D3<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients( GeometryData::GI_GAUSS_2 ),
                Triangle2D3<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients( GeometryData::GI_GAUSS_3 ),
                Triangle2D3<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients( GeometryData::GI_GAUSS_4 ),
                Triangle2D3<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients( GeometryData::GI_GAUSS_5 ),
                Triangle2D3<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients( GeometryData::GI_EXTENDED_GAUSS_1 ),
                Triangle2D3<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients( GeometryData::GI_EXTENDED_GAUSS_2 ),
                Triangle2D3<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients( GeometryData::GI_EXTENDED_GAUSS_3 ),
                Triangle2D3<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients( GeometryData::GI_EXTENDED_GAUSS_4 ),
                Triangle2D3<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients( GeometryData::GI_EXTENDED_GAUSS_5 ),
            }
        };
        return shape_functions_local_gradients;
    }


//*************************************************************************************
//*************************************************************************************

    /* Triangle/triangle intersection test routine,
    * by Tomas Moller, 1997.
    * See article "A Fast Triangle-Triangle Intersection Test",
    * Journal of Graphics Tools, 2(2), 1997
    *
    * Updated June 1999: removed the divisions -- a little faster now!
    * Updated October 1999: added {} to CROSS and SUB macros
    *
    * int NoDivTriTriIsect(float V0[3],float V1[3],float V2[3],
    *                      float U0[3],float U1[3],float U2[3])
    *
    * parameters: vertices of triangle 1: V0,V1,V2
    *             vertices of triangle 2: U0,U1,U2
    * result    : returns 1 if the triangles intersect, otherwise 0
    *
    */

    bool NoDivTriTriIsect( const Point<3,double>& V0,
                           const Point<3,double>& V1,
                           const Point<3,double>& V2,
                           const Point<3,double>& U0,
                           const Point<3,double>& U1,
                           const Point<3,double>& U2)
    {
        short index;
        double d1,d2;
        double du0,du1,du2,dv0,dv1,dv2;
        double du0du1,du0du2,dv0dv1,dv0dv2;
        double vp0,vp1,vp2;
        double up0,up1,up2;
        double bb,cc,max;
        array_1d<double,2> isect1, isect2;
        array_1d<double,3> D;
        array_1d<double,3> E1,E2;
        array_1d<double,3> N1,N2;


        const double epsilon =  1E-6;

// compute plane equation of triangle(V0,V1,V2) //
        noalias(E1) = V1-V0;
        noalias(E2) = V2-V0;
        MathUtils<double>::CrossProduct(N1, E1, E2);
        d1=-inner_prod(N1,V0);
// plane equation 1: N1.X+d1=0 //

// put U0,U1,U2 into plane equation 1 to compute signed distances to the plane//
        du0=inner_prod(N1,U0)+d1;
        du1=inner_prod(N1,U1)+d1;
        du2=inner_prod(N1,U2)+d1;

// coplanarity robustness check //
        if(fabs(du0)<epsilon) du0=0.0;
        if(fabs(du1)<epsilon) du1=0.0;
        if(fabs(du2)<epsilon) du2=0.0;

        du0du1=du0*du1;
        du0du2=du0*du2;

        if(du0du1>0.00 && du0du2>0.00)// same sign on all of them + not equal 0 ? //
            return false;                   // no intersection occurs //

// compute plane of triangle (U0,U1,U2) //
        noalias(E1) = U1 - U0;
        noalias(E2) = U2 - U0;
        MathUtils<double>::CrossProduct(N2, E1, E2);
        d2=-inner_prod(N2,U0);
// plane equation 2: N2.X+d2=0 //

// put V0,V1,V2 into plane equation 2 //
        dv0=inner_prod(N2,V0)+d2;
        dv1=inner_prod(N2,V1)+d2;
        dv2=inner_prod(N2,V2)+d2;

        if(fabs(dv0)<epsilon) dv0=0.0;
        if(fabs(dv1)<epsilon) dv1=0.0;
        if(fabs(dv2)<epsilon) dv2=0.0;

        dv0dv1=dv0*dv1;
        dv0dv2=dv0*dv2;

        if(dv0dv1>0.00 && dv0dv2>0.00)// same sign on all of them + not equal 0 ? //
            return false;                   // no intersection occurs //

// compute direction of intersection line //
        MathUtils<double>::CrossProduct(D, N1, N2);

// compute and index to the largest component of D //
        max=(double)fabs(D[0]);
        index=0;
        bb=(double)fabs(D[1]);
        cc=(double)fabs(D[2]);
        if(bb>max)
        {
            max=bb,index=1;
        }
        if(cc>max)
        {
            max=cc,index=2;
        }

// this is the simplified projection onto L//
        vp0=V0[index];
        vp1=V1[index];
        vp2=V2[index];

        up0=U0[index];
        up1=U1[index];
        up2=U2[index];


// compute interval for triangle 1 //
        double a,b,c,x0,x1;
        if( New_Compute_Intervals(vp0,vp1,vp2,dv0, dv1,dv2,dv0dv1,dv0dv2, a,b,c,x0,x1)==true )
        {
            return coplanar_tri_tri(N1,V0,V1,V2,U0,U1,U2);
        }

// compute interval for triangle 2 //
        double d,e,f,y0,y1;
        if( New_Compute_Intervals(up0, up1, up2, du0, du1, du2, du0du1, du0du2, d, e, f, y0, y1)==true )
        {
            return coplanar_tri_tri(N1,V0,V1,V2,U0,U1,U2);
        }


        double xx,yy,xxyy,tmp;
        xx=x0*x1;
        yy=y0*y1;
        xxyy=xx*yy;

        tmp=a*xxyy;
        isect1[0]=tmp+b*x1*yy;
        isect1[1]=tmp+c*x0*yy;

        tmp=d*xxyy;
        isect2[0]=tmp+e*xx*y1;
        isect2[1]=tmp+f*xx*y0;


        Sort(isect1[0],isect1[1]);
        Sort(isect2[0],isect2[1]);

        if(isect1[1]<isect2[0] || isect2[1]<isect1[0]) return false;
        return true;
    }


//*************************************************************************************
//*************************************************************************************


// sort so that a<=b //
    void Sort(double& a, double& b)
    {
        if(a>b)
        {
            double c;
            c=a;
            a=b;
            b=c;
        }
    }

//*************************************************************************************
//*************************************************************************************

    bool New_Compute_Intervals( double& VV0,
                                double& VV1,
                                double& VV2,
                                double& D0,
                                double& D1,
                                double& D2,
                                double& D0D1,
                                double& D0D2,
                                double& A,
                                double& B,
                                double& C,
                                double& X0,
                                double& X1
                              )
    {
        if(D0D1>0.00)
        {
            // here we know that D0D2<=0.0 //
            // that is D0, D1 are on the same side, D2 on the other or on the plane //
            A=VV2;
            B=(VV0-VV2)*D2;
            C=(VV1-VV2)*D2;
            X0=D2-D0;
            X1=D2-D1;
        }
        else if(D0D2>0.00)
        {
            // here we know that d0d1<=0.0 //
            A=VV1;
            B=(VV0-VV1)*D1;
            C=(VV2-VV1)*D1;
            X0=D1-D0;
            X1=D1-D2;
        }
        else if(D1*D2>0.00 || D0!=0.00)
        {
            // here we know that d0d1<=0.0 or that D0!=0.0 //
            A=VV0;
            B=(VV1-VV0)*D0;
            C=(VV2-VV0)*D0;
            X0=D0-D1;
            X1=D0-D2;
        }
        else if(D1!=0.00)
        {
            A=VV1;
            B=(VV0-VV1)*D1;
            C=(VV2-VV1)*D1;
            X0=D1-D0;
            X1=D1-D2;
        }
        else if(D2!=0.00)
        {
            A=VV2;
            B=(VV0-VV2)*D2;
            C=(VV1-VV2)*D2;
            X0=D2-D0;
            X1=D2-D1;
        }
        else
        {
            ///Triangles are coplanar
            return true;
        }

        return false;

    }

//*************************************************************************************
//*************************************************************************************

    bool coplanar_tri_tri( const array_1d<double, 3>& N,
                           const Point<3,double>& V0,
                           const Point<3,double>& V1,
                           const Point<3,double>& V2,
                           const Point<3,double>& U0,
                           const Point<3,double>& U1,
                           const Point<3,double>& U2)
    {
        array_1d<double, 3 > A;
        short i0,i1;

        // first project onto an axis-aligned plane, that maximizes the area //
        // of the triangles, compute indices: i0,i1. //
        A[0]=fabs(N[0]);
        A[1]=fabs(N[1]);
        A[2]=fabs(N[2]);
        if(A[0]>A[1])
        {
            if(A[0]>A[2])
            {
                i0=1;      // A[0] is greatest //
                i1=2;
            }
            else
            {
                i0=0;      // A[2] is greatest //
                i1=1;
            }
        }
        else   // A[0]<=A[1] //
        {
            if(A[2]>A[1])
            {
                i0=0;      // A[2] is greatest //
                i1=1;
            }
            else
            {
                i0=0;      // A[1] is greatest //
                i1=2;
            }
        }

        // test all edges of triangle 1 against the edges of triangle 2 //
        //std::cout<< "Proof One " << std::endl;
        if ( Edge_Against_Tri_Edges(i0, i1, V0,V1,U0,U1,U2)==true) return true;

        //std::cout<< "Proof Two " << std::endl;
        if ( Edge_Against_Tri_Edges(i0, i1, V1,V2,U0,U1,U2)==true) return true;

        //std::cout<< "Proof Three " << std::endl;
        if ( Edge_Against_Tri_Edges(i0, i1, V2,V0,U0,U1,U2)==true) return true;

        // finally, test if tri1 is totally contained in tri2 or vice versa //
        if (Point_In_Tri( i0, i1, V0, U0, U1, U2)==true) return true;
        if (Point_In_Tri( i0, i1, U0, V0, V1, V2)==true) return true;

        return false;
    }

//*************************************************************************************
//*************************************************************************************

    bool Edge_Against_Tri_Edges(const short& i0,
                                const short& i1,
                                const Point<3,double>& V0,
                                const Point<3,double>& V1,
                                const Point<3,double>&U0,
                                const Point<3,double>&U1,
                                const Point<3,double>&U2)
    {

        double Ax,Ay,Bx,By,Cx,Cy,e,d,f;
        Ax=V1[i0]-V0[i0];
        Ay=V1[i1]-V0[i1];
        // test edge U0,U1 against V0,V1 //

        //std::cout<< "Proof One B " << std::endl;
        if(Edge_Edge_Test(Ax, Ay, Bx, By, Cx, Cy, e, d, f, i0, i1, V0, U0, U1)==true) return true;
        // test edge U1,U2 against V0,V1 //
        //std::cout<< "Proof Two B " << std::endl;
        if(Edge_Edge_Test(Ax, Ay, Bx, By, Cx, Cy, e, d, f, i0, i1, V0, U1, U2)==true) return true;
        // test edge U2,U1 against V0,V1 //
        if(Edge_Edge_Test(Ax, Ay, Bx, By, Cx, Cy, e, d, f, i0, i1, V0, U2, U0)==true) return true;

        return false;
    }


//*************************************************************************************
//*************************************************************************************

//   this edge to edge test is based on Franlin Antonio's gem:
//   "Faster Line Segment Intersection", in Graphics Gems III,
//   pp. 199-202
    bool Edge_Edge_Test(double& Ax,
                        double& Ay,
                        double& Bx,
                        double& By,
                        double& Cx,
                        double& Cy,
                        double& e,
                        double& d,
                        double& f,
                        const short& i0,
                        const short& i1,
                        const Point<3,double>&V0,
                        const Point<3,double>&U0,
                        const Point<3,double>&U1)
    {
        Bx=U0[i0]-U1[i0];
        By=U0[i1]-U1[i1];
        Cx=V0[i0]-U0[i0];
        Cy=V0[i1]-U0[i1];
        f=Ay*Bx-Ax*By;
        d=By*Cx-Bx*Cy;

        if(std::fabs(f)<1E-10) f = 0.00;
        if(std::fabs(d)<1E-10) d = 0.00;


        if((f>0.00 && d>=0.00 && d<=f) || (f<0.00 && d<=0.00 && d>=f))
        {
            e=Ax*Cy-Ay*Cx;
            //std::cout<< "e =  "<< e << std::endl;
            if(f>0.00)
            {
                if(e>=0.00 && e<=f) return true;
            }
            else
            {
                if(e<=0.00 && e>=f) return true;
            }
        }
        return false;
    }

//*************************************************************************************
//*************************************************************************************


    bool Point_In_Tri(const short& i0,
                      const short& i1,
                      const Point<3,double>& V0,
                      const Point<3,double>& U0,
                      const Point<3,double>& U1,
                      const Point<3,double>& U2)
    {
        double a,b,c,d0,d1,d2;
        // is T1 completly inside T2? //
        // check if V0 is inside tri(U0,U1,U2) //
        a=U1[i1]-U0[i1];
        b=-(U1[i0]-U0[i0]);
        c=-a*U0[i0]-b*U0[i1];
        d0=a*V0[i0]+b*V0[i1]+c;

        a=U2[i1]-U1[i1];
        b=-(U2[i0]-U1[i0]);
        c=-a*U1[i0]-b*U1[i1];
        d1=a*V0[i0]+b*V0[i1]+c;

        a=U0[i1]-U2[i1];
        b=-(U0[i0]-U2[i0]);
        c=-a*U2[i0]-b*U2[i1];
        d2=a*V0[i0]+b*V0[i1]+c;
        if(d0*d1>0.0)
        {
            if(d0*d2>0.0) return true;
        }

        return false;
    }

//*************************************************************************************
//*************************************************************************************


    inline bool TriBoxOverlap(Point<3, double>& boxcenter, Point<3, double>& boxhalfsize, std::vector< Point<3, double> >& triverts)
    {

        /*    use separating axis theorem to test overlap between triangle and box */
        /*    need to test for overlap in these directions: */
        /*    1) the {x,y,z}-directions (actually, since we use the AABB of the triangle */
        /*       we do not even need to test these) */
        /*    2) normal of the triangle */
        /*    3) crossproduct(edge from tri, {x,y,z}-directin) */
        /*       this gives 3x3=9 more tests */

        double min,max,d,p0,p1,p2,rad,fex,fey,fez;
        array_1d<double,3 > v0,v1,v2;
        array_1d<double,3 > axis;
        array_1d<double,3 > normal, e0, e1 ,e2;
//
// 		  /* This is the fastest branch on Sun */
// 		  /* move everything so that the boxcenter is in (0,0,0) */
        noalias(v0) = triverts[0]- boxcenter;
        noalias(v1) = triverts[1]- boxcenter;
        noalias(v2) = triverts[2]- boxcenter;
//
// 		  /* compute triangle edges */
        noalias(e0) = v1 - v0;      /* tri edge 0 */
        noalias(e1) = v2 - v1;      /* tri edge 1 */
        noalias(e2) = v0 - v2;      /* tri edge 2 */
//
// 		  /* Bullet 3:  */
// 		  /*  test the 9 tests first (this was faster) */
        fex = fabs(e0[0]);
        fey = fabs(e0[1]);
        fez = fabs(e0[2]);
        //AXISTEST_X01(e0[2], e0[1], fez, fey);
        if (AxisTest_X01(e0[2], e0[1], fez, fey, p0, p2, min,max, rad, v0, v2, boxhalfsize)==0) return false;
        //AXISTEST_Y02(e0[2], e0[0], fez, fex);
        if(AxisTest_Y02( e0[2], e0[0], fez, fex, p0, p2, min,max, rad, v0, v2, boxhalfsize)==0) return false;
        //AXISTEST_Z12(e0[1], e0[0], fey, fex);
        //if(AxisTest_Z12(e0[1], e0[0], fey, fex, p1, p2, min, max, rad, v1,v2, boxhalfsize )==0) return false;



        fex = fabs(e1[0]);
        fey = fabs(e1[1]);
        fez = fabs(e1[2]);
        //AXISTEST_X01(e1[2], e1[1], fez, fey);
        if( AxisTest_X01(e1[2], e1[1], fez, fey, p0, p2, min,max, rad, v0, v2, boxhalfsize)==0) return false;
        //AXISTEST_Y02(e1[2], e1[0], fez, fex);
        if(AxisTest_Y02( e1[2], e1[0], fez, fex, p0, p2, min,max, rad, v0, v2, boxhalfsize)==0) return false;
        //AXISTEST_Z0(e1[1], e1[0], fey, fex);
        //if( AxisTest_Z0(e1[1], e1[0], fey, fex, p0,  p1, min, max, rad, v0, v1,boxhalfsize)==0) return false;


        fex = fabs(e2[0]);
        fey = fabs(e2[1]);
        fez = fabs(e2[2]);
        //AXISTEST_X2(e2[2], e2[1], fez, fey);
        if (AxisTest_X2(e2[2], e2[1], fez, fey, p0, p1, min, max, rad, v0, v1, boxhalfsize )==0) return false;
        //AXISTEST_Y1(e2[2], e2[0], fez, fex);
        if (AxisTest_Y1(e2[2], e2[0], fez, fex, p0, p1, min, max, rad, v0, v1, boxhalfsize )==0) return false;
        //AXISTEST_Z12(e2[1], e2[0], fey, fex);
        //if(AxisTest_Z12(e2[1], e2[0], fey, fex, p1, p2, min, max, rad, v1,v2, boxhalfsize ) ==0) return false;


        /* Bullet 1: */
        /*  first test overlap in the {x,y,z}-directions */
        /*  find min, max of the triangle each direction, and test for overlap in */
        /*  that direction -- this is equivalent to testing a minimal AABB around */
        /*  the triangle against the AABB */

        /* test in X-direction */
        FindMinMax(v0[0],v1[0],v2[0],min,max);
        if(min>boxhalfsize[0] || max<-boxhalfsize[0]) return false;

        /* test in Y-direction */
        FindMinMax(v0[1],v1[1],v2[1],min,max);
        if(min>boxhalfsize[1] || max<-boxhalfsize[1]) return false;

        /* test in Z-direction */
        FindMinMax(v0[2],v1[2],v2[2],min,max);
        if(min>boxhalfsize[2] || max<-boxhalfsize[2]) return false;

        /* Bullet 2: */
        /*  test if the box intersects the plane of the triangle */
        /*  compute plane equation of triangle: normal*x+d=0 */
        MathUtils<double>::CrossProduct(normal, e0, e1);
        d=-inner_prod(normal,v0);  /* plane eq: normal.x+d=0 */
        if(!planeBoxOverlap(normal,d,boxhalfsize)) return false;

        return true;   /* box and triangle overlaps */
    }



//*************************************************************************************
//*************************************************************************************

    void FindMinMax(const double& x0,
                    const double& x1,
                    const double& x2,
                    double& min,
                    double& max)
    {
        min = max = x0;
        if(x1<min) min=x1;
        if(x1>max) max=x1;
        if(x2<min) min=x2;
        if(x2>max) max=x2;
    }

//*************************************************************************************
//*************************************************************************************

    bool planeBoxOverlap(const array_1d<double,3 >& normal,  const double& d, const array_1d<double,3 >& maxbox)
    {
        int q;
        array_1d<double,3 >  vmin,vmax;
        for(q=0; q<=2; q++)
        {
            if(normal[q]>0.00)
            {
                vmin[q]=-maxbox[q];
                vmax[q]= maxbox[q];
            }
            else
            {
                vmin[q]=maxbox[q];
                vmax[q]=-maxbox[q];
            }
        }
        if(inner_prod(normal,vmin)+d>0.00) return false;
        if(inner_prod(normal,vmax)+d>=0.00) return true;

        return false;
    }
//*************************************************************************************
//*************************************************************************************

    /*======================== X-tests ========================*/
    unsigned int  AxisTest_X01(double& a,   double& b,
                               double& fa,  double& fb,
                               double& p0,  double& p2,
                               double& min, double& max, double& rad,
                               array_1d<double,3 >& v0,
                               array_1d<double,3 >& v2,
                               Point<3, double>& boxhalfsize
                              )
    {
        p0 = a*v0[1] - b*v0[2];
        p2 = a*v2[1] - b*v2[2];
        if(p0<p2)
        {
            min=p0;
            max=p2;
        }
        else
        {
            min=p2;
            max=p0;
        }
        rad = fa * boxhalfsize[1] + fb * boxhalfsize[2];
        if(min>rad || max<-rad) return 0;
        else return 1;
    }

    unsigned int  AxisTest_X2(double& a,   double& b,
                              double& fa,  double& fb,
                              double& p0,  double& p1,
                              double& min, double& max, double& rad,
                              array_1d<double,3 >& v0,
                              array_1d<double,3 >& v1,
                              Point<3, double>& boxhalfsize
                             )
    {
        p0 = a*v0[1] - b*v0[2];
        p1 = a*v1[1] - b*v1[2];
        if(p0<p1)
        {
            min=p0;
            max=p1;
        }
        else
        {
            min=p1;
            max=p0;
        }
        rad = fa * boxhalfsize[1] + fb * boxhalfsize[2];
        if(min>rad || max<-rad) return 0;
        else return 1;
    }
//*************************************************************************************
//*************************************************************************************

    /*======================== Y-tests ========================*/
    unsigned int  AxisTest_Y02(double& a, double& b,
                               double& fa,  double& fb,
                               double& p0,  double& p2,
                               double& min, double& max, double& rad,
                               array_1d<double,3 >& v0,
                               array_1d<double,3 >& v2,
                               Point<3, double>& boxhalfsize
                              )
    {

        p0 = -a*v0[0] + b*v0[2];
        p2 = -a*v2[0] + b*v2[2];
        if(p0<p2)
        {
            min=p0;
            max=p2;
        }
        else
        {
            min=p2;
            max=p0;
        }
        rad = fa * boxhalfsize[0] + fb * boxhalfsize[2];
        if(min>rad || max<-rad) return 0;
        else return 1;
    }

    unsigned int  AxisTest_Y1(double& a,   double& b,
                              double& fa,  double& fb,
                              double& p0,  double& p1,
                              double& min, double& max, double& rad,
                              array_1d<double,3 >& v0,
                              array_1d<double,3 >& v1,
                              Point<3, double>& boxhalfsize
                             )

    {
        p0 = -a*v0[0] + b*v0[2];
        p1 = -a*v1[0] + b*v1[2];
        if(p0<p1)
        {
            min=p0;
            max=p1;
        }
        else
        {
            min=p1;
            max=p0;
        }
        rad = fa * boxhalfsize[0] + fb * boxhalfsize[2];
        if(min>rad || max<-rad) return 0;
        else return 1;
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

    template<class TOtherPointType> friend class Triangle2D3;

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
/**
 * input stream functions
 */
template<class TPointType> inline std::istream& operator >> (
    std::istream& rIStream,
    Triangle2D3<TPointType>& rThis );
/**
 * output stream functions
 */
template<class TPointType> inline std::ostream& operator << (
    std::ostream& rOStream,
    const Triangle2D3<TPointType>& rThis )
{
    rThis.PrintInfo( rOStream );
    rOStream << std::endl;
    rThis.PrintData( rOStream );
    return rOStream;
}

///@}

template<class TPointType> const
GeometryData Triangle2D3<TPointType>::msGeometryData(
    2, 2, 2,
    GeometryData::GI_GAUSS_1,
    Triangle2D3<TPointType>::AllIntegrationPoints(),
    Triangle2D3<TPointType>::AllShapeFunctionsValues(),
    AllShapeFunctionsLocalGradients()
);
}// namespace Kratos.

#endif // KRATOS_TRIANGLE_2D_3_H_INCLUDED  defined 

