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



#if !defined(KRATOS_TETRAHEDRA_3D_10_H_INCLUDED )
#define  KRATOS_TETRAHEDRA_3D_10_H_INCLUDED

// System includes

// External includes

// Project includes
#include "geometries/triangle_3d_6.h"
#include "integration/tetrahedron_gauss_legendre_integration_points.h"


namespace Kratos
{
/**
 * An eight node hexahedra geometry with linear shape functions
 */

template<class TPointType> class Tetrahedra3D10 : public Geometry<TPointType>
{
public:
    /**
     * Type Definitions
     */

    /**
     * Geometry as base class.
     */
    typedef Geometry<TPointType> BaseType;

    /**
     * Type of edge and face geometries
     */
    typedef Line3D3<TPointType> EdgeType;
    typedef Triangle3D6<TPointType> FaceType;

    /**
     * Pointer definition of Tetrahedra3D10
     */
    KRATOS_CLASS_POINTER_DEFINITION( Tetrahedra3D10 );

    /**
     * Integration methods implemented in geometry.
     */
    typedef GeometryData::IntegrationMethod IntegrationMethod;

    /**
     * A Vector of counted pointers to Geometries. Used for
     * returning edges of the geometry.
     */
    typedef typename BaseType::GeometriesArrayType GeometriesArrayType;

    /**
     * Redefinition of template parameter TPointType.
     */
    typedef TPointType PointType;

    /**
     * Type used for indexing in geometry class.std::size_t used for indexing
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
     * Array of counted pointers to point. This type used to hold
     * geometry's points.
     */
    typedef typename BaseType::PointsArrayType PointsArrayType;

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
     * Type of the normal vector used for normal to edges in geomety.
     */
    typedef typename BaseType::NormalType NormalType;

    /**
     * Type of coordinates array
     */
    typedef typename BaseType::CoordinatesArrayType CoordinatesArrayType;

    /**
     * Type of Matrix
     */
    typedef Matrix MatrixType;


    /**
     * Life Cycle
     */

//     Tetrahedra3D10( const PointType& Point1, const PointType& Point2,
//                     const PointType& Point3, const PointType& Point4 ,
//                     const PointType& Point5, const PointType& Point6,
//                     const PointType& Point7, const PointType& Point8 ,
//                     const PointType& Point9, const PointType& Point10
//                   )
//         : BaseType( PointsArrayType(), &msGeometryData )
//     {
//         this->Points().reserve( 10 );
//         this->Points().push_back( typename PointType::Pointer( new PointType( Point1 ) ) );
//         this->Points().push_back( typename PointType::Pointer( new PointType( Point2 ) ) );
//         this->Points().push_back( typename PointType::Pointer( new PointType( Point3 ) ) );
//         this->Points().push_back( typename PointType::Pointer( new PointType( Point4 ) ) );
//         this->Points().push_back( typename PointType::Pointer( new PointType( Point5 ) ) );
//         this->Points().push_back( typename PointType::Pointer( new PointType( Point6 ) ) );
//         this->Points().push_back( typename PointType::Pointer( new PointType( Point7 ) ) );
//         this->Points().push_back( typename PointType::Pointer( new PointType( Point8 ) ) );
//         this->Points().push_back( typename PointType::Pointer( new PointType( Point9 ) ) );
//         this->Points().push_back( typename PointType::Pointer( new PointType( Point10 ) ) );
//     }

    Tetrahedra3D10( typename PointType::Pointer pPoint1,
                    typename PointType::Pointer pPoint2,
                    typename PointType::Pointer pPoint3,
                    typename PointType::Pointer pPoint4,
                    typename PointType::Pointer pPoint5,
                    typename PointType::Pointer pPoint6,
                    typename PointType::Pointer pPoint7,
                    typename PointType::Pointer pPoint8,
                    typename PointType::Pointer pPoint9,
                    typename PointType::Pointer pPoint10
                  )
        : BaseType( PointsArrayType(), &msGeometryData )
    {
        this->Points().reserve( 10 );
        this->Points().push_back( pPoint1 );
        this->Points().push_back( pPoint2 );
        this->Points().push_back( pPoint3 );
        this->Points().push_back( pPoint4 );
        this->Points().push_back( pPoint5 );
        this->Points().push_back( pPoint6 );
        this->Points().push_back( pPoint7 );
        this->Points().push_back( pPoint8 );
        this->Points().push_back( pPoint9 );
        this->Points().push_back( pPoint10 );
    }

    Tetrahedra3D10( const PointsArrayType& ThisPoints )
        : BaseType( ThisPoints, &msGeometryData )
    {
        if ( this->PointsNumber() != 10 )
            KRATOS_THROW_ERROR( std::invalid_argument,
                                "Invalid points number. Expected 10, given " ,
                                this->PointsNumber() );
    }

    /**
     * Copy constructor.
     * Construct this geometry as a copy of given geometry.
     *
     * @note This copy constructor don't copy the points and new
     * geometry shares points with given source geometry. It's
     * obvious that any change to this new geometry's point affect
     * source geometry's points too.
     */
    Tetrahedra3D10( Tetrahedra3D10 const& rOther )
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
    template<class TOtherPointType> Tetrahedra3D10( Tetrahedra3D10<TOtherPointType> const& rOther )
        : BaseType( rOther )
    {
    }

    /// Destructor. Does nothing!!!
    virtual ~Tetrahedra3D10() {}

    GeometryData::KratosGeometryFamily GetGeometryFamily()
    {
        return GeometryData::Kratos_Tetrahedra;
    }

    GeometryData::KratosGeometryType GetGeometryType()
    {
        return GeometryData::Kratos_Tetrahedra3D10;
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
    Tetrahedra3D10& operator=( const Tetrahedra3D10& rOther )
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
    Tetrahedra3D10& operator=( Tetrahedra3D10<TOtherPointType> const & rOther )
    {
        BaseType::operator=( rOther );

        return *this;
    }


    /**
     * Operations
     */

    typename BaseType::Pointer Create( PointsArrayType const& ThisPoints ) const
    {
        return typename BaseType::Pointer( new Tetrahedra3D10( ThisPoints ) );
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
        Geometry< Point<3> >::Pointer p_clone( new Tetrahedra3D10< Point<3> >( NewPoints ) );

        return p_clone;
    }


    //lumping factors for the calculation of the lumped mass matrix
    virtual Vector& LumpingFactors( Vector& rResult ) const
    {
        if(rResult.size() != 10)
            rResult.resize( 10, false );
        std::fill( rResult.begin(), rResult.end(), 1.00 / 10.00 );
        return rResult;
    }


    /**
     * Informations
     */

    /**
     * This method calculates and returns Length or charactereistic
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
     *
     * :TODO: might be necessary to reimplement
     */
    virtual double Length() const
    {
        return sqrt( fabs( this->DeterminantOfJacobian( PointType() ) ) );
    }

    /**
     * This method calculates and returns area or surface area of
     * this geometry depending to it's dimension. For one dimensional
     * geometry it returns zero, for two dimensional it gives area
     * and for three dimensional geometries it gives surface area.
     *
     *
     * @return double value contains area or surface area.
     * @see Length()
     * @see Volume()
     * @see DomainSize()
     *
     * :TODO: might be necessary to reimplement
     */
    virtual double Area() const
    {
        return Volume();
    }


    virtual double Volume() const //Not a closed formula for a quadratic tetrahedra
    {

        Vector temp;
        this->DeterminantOfJacobian( temp, msGeometryData.DefaultIntegrationMethod() );
        const IntegrationPointsArrayType& integration_points = this->IntegrationPoints( msGeometryData.DefaultIntegrationMethod() );
        double Volume = 0.00;

        for ( unsigned int i = 0; i < integration_points.size(); i++ )
        {
            Volume += temp[i] * integration_points[i].Weight();
        }

        //KRATOS_WATCH(temp)
        return Volume;
    }


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
     *
     * :TODO: might be necessary to reimplement
     */
    virtual double DomainSize() const
    {
        return  Volume();
    }


    virtual bool IsInside( const CoordinatesArrayType& rPoint, CoordinatesArrayType& rResult )
    {
        this->PointLocalCoordinates( rResult, rPoint );

        if ( rResult[0] >= 0.0 - 1.0e-8 && rResult[0] <= 1.0 + 1.0e-8 )
            if ( rResult[1] >= 0.0 - 1.0e-8 && rResult[1] <= 1.0 + 1.0e-8 )
                if ( rResult[2] >= 0.0 - 1.0e-8 && rResult[2] <= 1.0 + 1.0e-8 )
                    if ((( 1.0 - ( rResult[0] + rResult[1] + rResult[2] ) ) >= 0.0 - 1.0e-8 ) && (( 1.0 - ( rResult[0] + rResult[1] + rResult[2] ) ) <= 1.0 + 1.0e-8 ) )
                        return true;

        return false;
    }


    /** This method gives you number of all edges of this
    geometry.
    @return SizeType containes number of this geometry edges.
    @see Edges()
    @see Edge()
     */
    // will be used by refinement algorithm, thus uncommented. janosch.
    virtual SizeType EdgesNumber() const
    {
        return 6;
    }

    virtual SizeType FacesNumber() const
    {
        return 4;
    }

    /** This method gives you all edges of this geometry.

    @return GeometriesArrayType containes this geometry edges.
    @see EdgesNumber()
    @see Edge()
    */
    virtual GeometriesArrayType Edges( void )
    {
        GeometriesArrayType edges = GeometriesArrayType();
        typedef typename Geometry<TPointType>::Pointer EdgePointerType;
        edges.push_back( EdgePointerType( new EdgeType(
                                              this->pGetPoint( 0 ),
                                              this->pGetPoint( 4 ),
                                              this->pGetPoint( 1 ) ) ) );
        edges.push_back( EdgePointerType( new EdgeType(
                                              this->pGetPoint( 1 ),
                                              this->pGetPoint( 5 ),
                                              this->pGetPoint( 2 ) ) ) );
        edges.push_back( EdgePointerType( new EdgeType(
                                              this->pGetPoint( 2 ),
                                              this->pGetPoint( 6 ),
                                              this->pGetPoint( 0 ) ) ) );
        edges.push_back( EdgePointerType( new EdgeType(
                                              this->pGetPoint( 0 ),
                                              this->pGetPoint( 7 ),
                                              this->pGetPoint( 3 ) ) ) );
        edges.push_back( EdgePointerType( new EdgeType(
                                              this->pGetPoint( 1 ),
                                              this->pGetPoint( 8 ),
                                              this->pGetPoint( 3 ) ) ) );
        edges.push_back( EdgePointerType( new EdgeType(
                                              this->pGetPoint( 2 ),
                                              this->pGetPoint( 9 ),
                                              this->pGetPoint( 3 ) ) ) );

        return edges;
    }

    virtual GeometriesArrayType Faces( void )
    {
        GeometriesArrayType faces = GeometriesArrayType();
        typedef typename Geometry<TPointType>::Pointer FacePointerType;
        faces.push_back( FacePointerType( new FaceType(
                                              this->pGetPoint( 0 ),
                                              this->pGetPoint( 2 ),
                                              this->pGetPoint( 1 ),
                                              this->pGetPoint( 6 ),
                                              this->pGetPoint( 5 ),
                                              this->pGetPoint( 4 ) ) ) );
        faces.push_back( FacePointerType( new FaceType(
                                              this->pGetPoint( 0 ),
                                              this->pGetPoint( 3 ),
                                              this->pGetPoint( 2 ),
                                              this->pGetPoint( 7 ),
                                              this->pGetPoint( 9 ),
                                              this->pGetPoint( 6 ) ) ) );
        faces.push_back( FacePointerType( new FaceType(
                                              this->pGetPoint( 0 ),
                                              this->pGetPoint( 1 ),
                                              this->pGetPoint( 3 ),
                                              this->pGetPoint( 4 ),
                                              this->pGetPoint( 8 ),
                                              this->pGetPoint( 7 ) ) ) );
        faces.push_back( FacePointerType( new FaceType(
                                              this->pGetPoint( 2 ),
                                              this->pGetPoint( 3 ),
                                              this->pGetPoint( 1 ),
                                              this->pGetPoint( 9 ),
                                              this->pGetPoint( 8 ),
                                              this->pGetPoint( 5 ) ) ) );
        return faces;
    }



    /**
     * Shape Function
     */

    /**
     * Calculates the value of a given shape function at a given point.
     *
     * @param ShapeFunctionIndex The number of the desired shape function
     * @param rPoint the given point in local coordinates at which the
     * value of the shape function is calculated
     *
     * @return the value of the shape function at the given point
     * TODO: TO BE VERIFIED
     */
    virtual double ShapeFunctionValue( IndexType ShapeFunctionIndex,
                                       const CoordinatesArrayType& rPoint ) const
    {
        double fourthCoord = 1.0 - ( rPoint[0] + rPoint[1] + rPoint[2] );

        switch ( ShapeFunctionIndex )
        {
        case 0:
            return( fourthCoord*( 2*fourthCoord - 1 ) );
        case 1:
            return( rPoint[0]*( 2*rPoint[0] - 1 ) );
        case 2:
            return( rPoint[1]*( 2*rPoint[1] - 1 ) );
        case 3:
            return( rPoint[2]*( 2*rPoint[2] - 1 ) );
        case 4:
            return( 4*fourthCoord*rPoint[0] );
        case 5:
            return( 4*rPoint[0]*rPoint[1] );
        case 6:
            return( 4*fourthCoord*rPoint[1] );
        case 7:
            return( 4*fourthCoord*rPoint[2] );
        case 8:
            return( 4*rPoint[0]*rPoint[2] );
        case 9:
            return( 4*rPoint[1]*rPoint[2] );
        default:
            KRATOS_THROW_ERROR( std::logic_error,
                                "Wrong index of shape function!" , *this );
        }

        return 0;
    }


    /**
     * Input and output
     */

    /**
     * Turn back information as a string.
     *
     * @return String contains information about this geometry.
     * @see PrintData()
     * @see PrintInfo()
     */
    virtual std::string Info() const
    {
        return "3 dimensional tetrahedra with ten nodes in 3D space";
    }

    /**
     * Print information about this object.
     *
     * @param rOStream Stream to print into it.
     * @see PrintData()
     * @see Info()
     */
    virtual void PrintInfo( std::ostream& rOStream ) const
    {
        rOStream << "3 dimensional tetrahedra with ten nodes in 3D space";
    }

    /**
     * Print geometry's data into given stream.
     * Prints it's points by the order they stored in the geometry
     * and then center point of geometry.
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

    Tetrahedra3D10(): BaseType( PointsArrayType(), &msGeometryData ) {}


    /**
     * Private Operations
     */

    /**
     * TODO: TO BE VERIFIED
     */
    /**
     * Calculates the gradients in terms of local coordinateds
     * of all shape functions in a given point.
     *
     * @param rPoint the current point at which the gradients are calculated
     * @return the gradients of all shape functions
     * \f$ \frac{\partial N^i}{\partial \xi_j} \f$
     */
    static Matrix& ShapeFunctionsLocalGradientsPrivate( Matrix& result,
            const CoordinatesArrayType& rPoint )
    {
        double fourthCoord = 1.0 - ( rPoint[0] + rPoint[1] + rPoint[2] );
        double fourthCoord_DX = -1.0;
        double fourthCoord_DY = -1.0;
        double fourthCoord_DZ = -1.0;

        if ( result.size1() != 10 || result.size2() != 3 )
            result.resize( 10, 3, false );

        result( 0, 0 ) = ( 4 * fourthCoord - 1.0 ) * fourthCoord_DX;

        result( 0, 1 ) = ( 4 * fourthCoord - 1.0 ) * fourthCoord_DY;

        result( 0, 2 ) = ( 4 * fourthCoord - 1.0 ) * fourthCoord_DZ;

        result( 1, 0 ) =  4 * rPoint[0] - 1.0;

        result( 1, 1 ) =  0.0;

        result( 1, 2 ) =  0.0;

        result( 2, 0 ) =  0.0;

        result( 2, 1 ) =  4 * rPoint[1] - 1.0;

        result( 2, 2 ) =  0.0;

        result( 3, 0 ) =  0.0;

        result( 3, 1 ) =  0.0 ;

        result( 3, 2 ) =  4 * rPoint[2] - 1.0 ;

        result( 4, 0 ) =  4 * fourthCoord_DX * rPoint[0] + 4 * fourthCoord;

        result( 4, 1 ) =  4 * fourthCoord_DY * rPoint[0];

        result( 4, 2 ) =  4 * fourthCoord_DZ * rPoint[0];

        result( 5, 0 ) =  4 * rPoint[1];

        result( 5, 1 ) =  4 * rPoint[0];

        result( 5, 2 ) =  0.0;

        result( 6, 0 ) =  4 * fourthCoord_DX * rPoint[1];

        result( 6, 1 ) =  4 * fourthCoord_DY * rPoint[1] + 4 * fourthCoord;

        result( 6, 2 ) =  4 * fourthCoord_DZ * rPoint[1];

        result( 7, 0 ) =  4 * fourthCoord_DX * rPoint[2];

        result( 7, 1 ) =  4 * fourthCoord_DY * rPoint[2];

        result( 7, 2 ) =  4 * fourthCoord_DZ * rPoint[2] + 4 * fourthCoord;

        result( 8, 0 ) =  4 * rPoint[2];

        result( 8, 1 ) =  0.0;

        result( 8, 2 ) =  4 * rPoint[0];

        result( 9, 0 ) =  0.0;

        result( 9, 1 ) =  4 * rPoint[2];

        result( 9, 2 ) =  4 * rPoint[1];

        return result;
    }

    /**
     * TODO: TO BE VERIFIED
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
        IntegrationPointsContainerType all_integration_points =
            AllIntegrationPoints();
        IntegrationPointsArrayType integration_points =
            all_integration_points[ThisMethod];
        //number of integration points
        const int integration_points_number = integration_points.size();
        //number of nodes in current geometry
        const int points_number = 10;
        //setting up return matrix
        Matrix shape_function_values( integration_points_number, points_number );
        //loop over all integration points

        for ( int pnt = 0; pnt < integration_points_number; pnt++ )
        {
            double fourthCoord = 1.0 - ( integration_points[pnt].X() + integration_points[pnt].Y() + integration_points[pnt].Z() );

            shape_function_values( pnt, 0 ) = fourthCoord * ( 2.0 * fourthCoord - 1.0 );

            shape_function_values( pnt, 1 ) = integration_points[pnt].X() * ( 2.0 * integration_points[pnt].X() - 1.0 );

            shape_function_values( pnt, 2 ) = integration_points[pnt].Y() * ( 2.0 * integration_points[pnt].Y() - 1.0 );

            shape_function_values( pnt, 3 ) = integration_points[pnt].Z() * ( 2.0 * integration_points[pnt].Z() - 1.0 );

            shape_function_values( pnt, 4 ) = 4.0 * fourthCoord * integration_points[pnt].X();

            shape_function_values( pnt, 5 ) = 4.0 * integration_points[pnt].X() * integration_points[pnt].Y();

            shape_function_values( pnt, 6 ) = 4.0 * fourthCoord * integration_points[pnt].Y();

            shape_function_values( pnt, 7 ) = 4.0 * fourthCoord * integration_points[pnt].Z();

            shape_function_values( pnt, 8 ) = 4.0 * integration_points[pnt].X() * integration_points[pnt].Z();

            shape_function_values( pnt, 9 ) = 4.0 * integration_points[pnt].Y() * integration_points[pnt].Z();

        }

        return shape_function_values;
    }

    /**
     * TODO: TO BE VERIFIED
     */
    /**
     * Calculates the local gradients of all shape functions in all integration points.
     * Integration points are expected to be given in local coordinates
     *
     * @param ThisMethod the current integration method
     * @return the vector of the gradients of all shape functions
     * in each integration point
     *
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
        //loop over all integration points

        for ( int pnt = 0; pnt < integration_points_number; pnt++ )
        {
            double fourthCoord = 1.0 - ( integration_points[pnt].X() + integration_points[pnt].Y() + integration_points[pnt].Z() );
            double fourthCoord_DX = -1.0;
            double fourthCoord_DY = -1.0;
            double fourthCoord_DZ = -1.0;

            Matrix result = ZeroMatrix( 10, 3 );
            result( 0, 0 ) = ( 4 * fourthCoord - 1.0 ) * fourthCoord_DX;
            result( 0, 1 ) = ( 4 * fourthCoord - 1.0 ) * fourthCoord_DY;
            result( 0, 2 ) = ( 4 * fourthCoord - 1.0 ) * fourthCoord_DZ;
            result( 1, 0 ) =  4 * integration_points[pnt].X() - 1.0;
            result( 1, 1 ) =  0.0;
            result( 1, 2 ) =  0.0;
            result( 2, 0 ) =  0.0;
            result( 2, 1 ) =  4 * integration_points[pnt].Y() - 1.0;
            result( 2, 2 ) =  0.0;
            result( 3, 0 ) =  0.0;
            result( 3, 1 ) =  0.0;
            result( 3, 2 ) =  4 * integration_points[pnt].Z() - 1.0 ;
            result( 4, 0 ) =  4 * fourthCoord_DX * integration_points[pnt].X() + 4 * fourthCoord;
            result( 4, 1 ) =  4 * fourthCoord_DY * integration_points[pnt].X();
            result( 4, 2 ) =  4 * fourthCoord_DZ * integration_points[pnt].X();
            result( 5, 0 ) =  4 * integration_points[pnt].Y();
            result( 5, 1 ) =  4 * integration_points[pnt].X();
            result( 5, 2 ) =  0.0;
            result( 6, 0 ) =  4 * fourthCoord_DX * integration_points[pnt].Y();
            result( 6, 1 ) =  4 * fourthCoord_DY * integration_points[pnt].Y() + 4 * fourthCoord;
            result( 6, 2 ) =  4 * fourthCoord_DZ * integration_points[pnt].Y();
            result( 7, 0 ) =  4 * fourthCoord_DX * integration_points[pnt].Z();
            result( 7, 1 ) =  4 * fourthCoord_DY * integration_points[pnt].Z();
            result( 7, 2 ) =  4 * fourthCoord_DZ * integration_points[pnt].Z() + 4 * fourthCoord;
            result( 8, 0 ) =  4 * integration_points[pnt].Z();
            result( 8, 1 ) =  0.0;
            result( 8, 2 ) =  4 * integration_points[pnt].X();
            result( 9, 0 ) =  0.0;
            result( 9, 1 ) =  4 * integration_points[pnt].Z();
            result( 9, 2 ) =  4 * integration_points[pnt].Y();

            d_shape_f_values[pnt] = result;
        }

        return d_shape_f_values;
    }

    static const IntegrationPointsContainerType AllIntegrationPoints()
    {
        IntegrationPointsContainerType integration_points =
        {
            {
                Quadrature < TetrahedronGaussLegendreIntegrationPoints1,
                3, IntegrationPoint<3> >::GenerateIntegrationPoints(),
                Quadrature < TetrahedronGaussLegendreIntegrationPoints2,
                3, IntegrationPoint<3> >::GenerateIntegrationPoints(),
                Quadrature < TetrahedronGaussLegendreIntegrationPoints3,
                3, IntegrationPoint<3> >::GenerateIntegrationPoints(),
                Quadrature < TetrahedronGaussLegendreIntegrationPoints4,
                3, IntegrationPoint<3> >::GenerateIntegrationPoints(),
                Quadrature < TetrahedronGaussLegendreIntegrationPoints5,
                3, IntegrationPoint<3> >::GenerateIntegrationPoints(),
            }
        };
        return integration_points;
    }

    static const ShapeFunctionsValuesContainerType AllShapeFunctionsValues()
    {
        ShapeFunctionsValuesContainerType shape_functions_values =
        {
            {
                Tetrahedra3D10<TPointType>::CalculateShapeFunctionsIntegrationPointsValues(
                    GeometryData::GI_GAUSS_1 ),
                Tetrahedra3D10<TPointType>::CalculateShapeFunctionsIntegrationPointsValues(
                    GeometryData::GI_GAUSS_2 ),
                Tetrahedra3D10<TPointType>::CalculateShapeFunctionsIntegrationPointsValues(
                    GeometryData::GI_GAUSS_3 ),
                Tetrahedra3D10<TPointType>::CalculateShapeFunctionsIntegrationPointsValues(
                    GeometryData::GI_GAUSS_4 ),
                Tetrahedra3D10<TPointType>::CalculateShapeFunctionsIntegrationPointsValues(
                    GeometryData::GI_GAUSS_5 )
            }
        };
        return shape_functions_values;
    }

    /**
     * TODO: TO BE VERIFIED
     */
    static const ShapeFunctionsLocalGradientsContainerType
    AllShapeFunctionsLocalGradients()
    {
        ShapeFunctionsLocalGradientsContainerType shape_functions_local_gradients =
        {
            {
                Tetrahedra3D10<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients(
                    GeometryData::GI_GAUSS_1 ),
                Tetrahedra3D10<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients(
                    GeometryData::GI_GAUSS_2 ),
                Tetrahedra3D10<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients(
                    GeometryData::GI_GAUSS_3 ),
                Tetrahedra3D10<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients(
                    GeometryData::GI_GAUSS_4 ),
                Tetrahedra3D10<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients(
                    GeometryData::GI_GAUSS_5 )
            }
        };
        return shape_functions_local_gradients;
    }


    /**
     * Private Friends
     */

    template<class TOtherPointType> friend class Tetrahedra3D10;


    /**
     * Un accessible methods
     */


};// Class Tetrahedra3D10


/**
 * Input and output
 */

/**
 * input stream function
 */
template<class TPointType> inline std::istream& operator >> (
    std::istream& rIStream, Tetrahedra3D10<TPointType>& rThis );

/**
 * output stream function
 */
template<class TPointType> inline std::ostream& operator << (
    std::ostream& rOStream, const Tetrahedra3D10<TPointType>& rThis )
{
    rThis.PrintInfo( rOStream );
    rOStream << std::endl;
    rThis.PrintData( rOStream );

    return rOStream;
}


template<class TPointType> const
GeometryData Tetrahedra3D10<TPointType>::msGeometryData(
    3, 3, 3, GeometryData::GI_GAUSS_2,
    Tetrahedra3D10<TPointType>::AllIntegrationPoints(),
    Tetrahedra3D10<TPointType>::AllShapeFunctionsValues(),
    AllShapeFunctionsLocalGradients()
);

}// namespace Kratos.

#endif // KRATOS_TETRAHEDRA_3D_4_H_INCLUDED  defined 
