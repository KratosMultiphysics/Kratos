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

#if !defined(KRATOS_LINE_2D_2_H_INCLUDED )
#define  KRATOS_LINE_2D_2_H_INCLUDED

// System includes

// External includes

// Project includes
#include "geometries/geometry.h"
#include "integration/line_gauss_legendre_integration_points.h"
#include "integration/line_collocation_integration_points.h"


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
 * @class Line2D2
 * @ingroup KratosCore
 * @brief An two node 2D line geometry with linear shape functions
 * @details The node ordering corresponds with: 
 *      0----------1 --> u  
 * @author Riccardo Rossi
 * @author Janosch Stascheit
 * @author Felix Nagel
 */
template<class TPointType>

class Line2D2 : public Geometry<TPointType>
{
public:
    ///@}
    ///@name Type Definitions
    ///@{

    /// Geometry as base class.
    typedef Geometry<TPointType> BaseType;

    /// Pointer definition of Line2D2
    KRATOS_CLASS_POINTER_DEFINITION( Line2D2 );

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



//     Line2D2( const PointType& FirstPoint, const PointType& SecondPoint )
//         : BaseType( PointsArrayType(), &msGeometryData )
//     {
//         BaseType::Points().push_back( typename PointType::Pointer( new PointType( FirstPoint ) ) );
//         BaseType::Points().push_back( typename PointType::Pointer( new PointType( SecondPoint ) ) );
//     }

    Line2D2( typename PointType::Pointer pFirstPoint, typename PointType::Pointer pSecondPoint )
        : BaseType( PointsArrayType(), &msGeometryData )
    {
        BaseType::Points().push_back( pFirstPoint );
        BaseType::Points().push_back( pSecondPoint );
    }


    Line2D2( const PointsArrayType& ThisPoints )
        : BaseType( ThisPoints, &msGeometryData )
    {
        if ( BaseType::PointsNumber() != 2 )
            KRATOS_ERROR << "Invalid points number. Expected 2, given " << BaseType::PointsNumber() << std::endl;
    }

    /** Copy constructor.
    Construct this geometry as a copy of given geometry.

    @note This copy constructor don't copy the points and new
    geometry shares points with given source geometry. It's
    obvious that any change to this new geometry's point affect
    source geometry's points too.
    */
    Line2D2( Line2D2 const& rOther )
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
    template<class TOtherPointType> Line2D2( Line2D2<TOtherPointType> const& rOther )
        : BaseType( rOther )
    {
    }

    /// Destructor. Do nothing!!!
    ~Line2D2() override {}

    GeometryData::KratosGeometryFamily GetGeometryFamily() const override
    {
        return GeometryData::Kratos_Linear;
    }

    GeometryData::KratosGeometryType GetGeometryType() const override
    {
        return GeometryData::Kratos_Line2D2;
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
    Line2D2& operator=( const Line2D2& rOther )
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
    Line2D2& operator=( Line2D2<TOtherPointType> const & rOther )
    {
        BaseType::operator=( rOther );

        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    typename BaseType::Pointer Create( PointsArrayType const& ThisPoints ) const override
    {
        return typename BaseType::Pointer( new Line2D2( ThisPoints ) );
    }

    // Geometry< Point<3> >::Pointer Clone() const override
    // {
    //     Geometry< Point<3> >::PointsArrayType NewPoints;

    //     //making a copy of the nodes TO POINTS (not Nodes!!!)
    //     for ( IndexType i = 0 ; i < this->size() ; i++ )
    //     {
    //         NewPoints.push_back(Kratos::make_shared< Point<3> >((*this)[i]));
    //     }

    //     //creating a geometry with the new points
    //     Geometry< Point<3> >::Pointer p_clone( new Line2D2< Point<3> >( NewPoints ) );

    //     return p_clone;
    // }

    //lumping factors for the calculation of the lumped mass matrix
    Vector& LumpingFactors( Vector& rResult ) const override
    {
        if(rResult.size() != 2)
        {
               rResult.resize( 2, false );
        }

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
        const TPointType& FirstPoint  = BaseType::GetPoint(0);
        const TPointType& SecondPoint = BaseType::GetPoint(1);
        const double lx = FirstPoint.X() - SecondPoint.X();
        const double ly = FirstPoint.Y() - SecondPoint.Y();

        const double length = lx * lx + ly * ly;

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
        const TPointType& FirstPoint = BaseType::GetPoint(0);
        const TPointType& SecondPoint = BaseType::GetPoint(1);
        const double lx = FirstPoint.X() - SecondPoint.X();
        const double ly = FirstPoint.Y() - SecondPoint.Y();

        const double length = lx * lx + ly * ly;

        return std::sqrt( length );
    }

//      virtual void Bounding_Box(BoundingBox<TPointType, BaseType>& rResult) const
//              {
//                 //rResult.Geometry() = *(this);
//                 BaseType::Bounding_Box(rResult.LowPoint(), rResult.HighPoint());
//              }





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
        Matrix jacobian( 2, 1 );
        jacobian( 0, 0 ) = ( BaseType::GetPoint( 1 ).X() - BaseType::GetPoint( 0 ).X() ) * 0.5; //on the Gauss points (J is constant at each element)
        jacobian( 1, 0 ) = ( BaseType::GetPoint( 1 ).Y() - BaseType::GetPoint( 0 ).Y() ) * 0.5;

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
        Matrix jacobian( 2, 1 );
        jacobian( 0, 0 ) = ( (BaseType::GetPoint( 1 ).X() - DeltaPosition(1,0)) - (BaseType::GetPoint( 0 ).X() - DeltaPosition(0,0)) ) * 0.5; //on the Gauss points (J is constant at each element)
        jacobian( 1, 0 ) = ( (BaseType::GetPoint( 1 ).Y() - DeltaPosition(1,1)) - (BaseType::GetPoint( 0 ).Y() - DeltaPosition(0,1)) ) * 0.5;

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
        rResult.resize( 2, 1, false );
        //on the Gauss points (J is constant at each element)
        rResult( 0, 0 ) = ( BaseType::GetPoint( 1 ).X() - BaseType::GetPoint( 0 ).X() ) * 0.5;
        rResult( 1, 0 ) = ( BaseType::GetPoint( 1 ).Y() - BaseType::GetPoint( 0 ).Y() ) * 0.5;
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
        rResult.resize( 2, 1, false );
        //on the Gauss points (J is constant at each element)
        rResult( 0, 0 ) = ( BaseType::GetPoint( 1 ).X() - BaseType::GetPoint( 0 ).X() ) * 0.5;
        rResult( 1, 0 ) = ( BaseType::GetPoint( 1 ).Y() - BaseType::GetPoint( 0 ).Y() ) * 0.5;
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
        KRATOS_ERROR << "Jacobian is not square" << std::endl;
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
        KRATOS_ERROR << "Jacobian is not square" << std::endl;
        return rResult;
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
        KRATOS_ERROR << "Jacobian is not square" << std::endl;
        return rResult;
    }

    /** EdgesNumber
    @return SizeType containes number of this geometry edges.
    */
    SizeType EdgesNumber() const override
    {
        return 2;
    }


    /** FacesNumber
    @return SizeType containes number of this geometry edges/faces.
    */
    SizeType FacesNumber() const override
    {
      return EdgesNumber();
    }



    //Connectivities of faces required
    void NumberNodesInFaces (DenseVector<unsigned int>& NumberNodesInFaces) const override
    {
        if(NumberNodesInFaces.size() != 2 )
            NumberNodesInFaces.resize(2,false);

        // Lines have 1 node in edges/faces
        NumberNodesInFaces[0]=1;
        NumberNodesInFaces[1]=1;

    }

    void NodesInFaces (DenseMatrix<unsigned int>& NodesInFaces) const override
    {
        // faces in columns
        if(NodesInFaces.size1() != 2 || NodesInFaces.size2() != 2)
            NodesInFaces.resize(2,2,false);

        //face 1
        NodesInFaces(0,0)=0;//contrary node to the face
        NodesInFaces(1,0)=1;

        //face 2
        NodesInFaces(0,1)=1;//contrary node to the face
        NodesInFaces(1,1)=0;
    }

    ///@}
    ///@name Shape Function
    ///@{

    Vector& ShapeFunctionsValues (Vector &rResult, const CoordinatesArrayType& rCoordinates) const override
    {
        if(rResult.size() != 2)
        {
            rResult.resize(2, false);
        }

        rResult[0] =  0.5 * ( 1.0 - rCoordinates[0]);
        rResult[1] =  0.5 * ( 1.0 + rCoordinates[0]);

        return rResult;
    }

    double ShapeFunctionValue( IndexType ShapeFunctionIndex,
                                       const CoordinatesArrayType& rPoint ) const override
    {
        switch ( ShapeFunctionIndex )
        {
        case 0:
            return( 0.5 * ( 1.0 - rPoint[0] ) );
        case 1:
            return( 0.5 * ( 1.0 + rPoint[0] ) );
        default:
            KRATOS_ERROR << "Wrong index of shape function!" << *this << std::endl;
        }

        return 0;
    }

    ///@}
    ///@name Shape Function Integration Points Gradient
    ///@{

    ShapeFunctionsGradientsType& ShapeFunctionsIntegrationPointsGradients( ShapeFunctionsGradientsType& rResult, IntegrationMethod ThisMethod ) const override
    {
        KRATOS_ERROR << "Jacobian is not square" << std::endl;
        return rResult;
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
        return "1 dimensional line in 2D space";
    }

    /** Print information about this object.

    @param rOStream Stream to print into it.
    @see PrintData()
    @see Info()
    */
    void PrintInfo( std::ostream& rOStream ) const override
    {
        rOStream << "1 dimensional line in 2D space";
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
        BaseType::PrintData( rOStream );
        std::cout << std::endl;
        Matrix jacobian;
        Jacobian( jacobian, PointType() );
        rOStream << "    Jacobian\t : " << jacobian;
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
     * @param rPoint the current point at which the gradients are calculated
     * @return the gradients of all shape functions
     * \f$ \frac{\partial N^i}{\partial \xi_j} \f$
     */
    Matrix& ShapeFunctionsLocalGradients( Matrix& rResult,
            const CoordinatesArrayType& rPoint ) const override
    {
        // Setting up result matrix
        if(rResult.size1() != 2 || rResult.size2() != 1)
        {
            rResult.resize( 2, 1, false );
        }
        noalias( rResult ) = ZeroMatrix( 2, 1 );
        rResult( 0, 0 ) = - 0.5;
        rResult( 1, 0 ) =   0.5;

        return( rResult );
    }

    /**
     * It computes the area normal of the geometry
     * @param rPointLocalCoordinates Refernce to the local coordinates of the
     * point in where the unit normal is to be computed
     * @return The area normal in the given point
     */
    array_1d<double, 3> AreaNormal(const CoordinatesArrayType& rPointLocalCoordinates) const override
    {
        // We define the normal
        array_1d<double,3> normal;

        // We get the local points
        const TPointType& first_point  = BaseType::GetPoint(0);
        const TPointType& second_point = BaseType::GetPoint(1);

        // We compute the normal
        normal[0] = second_point[1] -  first_point[1];
        normal[1] =  first_point[0] - second_point[0];
        normal[2] = 0.0;

        return normal;
    }

    /**
     * returns the local coordinates of all nodes of the current geometry
     * @param rResult a Matrix object that will be overwritten by the result
     * @return the local coordinates of all nodes
     */
    Matrix& PointsLocalCoordinates( Matrix& rResult ) const override
    {
        if(rResult.size1() != 2 || rResult.size2() != 1)
        {
            rResult.resize( 2, 1, false );
        }
        noalias( rResult ) = ZeroMatrix( 2, 1 );
        rResult( 0, 0 ) = -1.0;
        rResult( 1, 0 ) =  1.0;
        return rResult;
    }

    /**
     * returns the shape function gradients in an arbitrary point,
     * given in local coordinates
     * @param rResult the matrix of gradients, will be overwritten
     * with the gradients for all
     * shape functions in given point
     * @param rPoint the given point the gradients are calculated in
     */
    virtual Matrix& ShapeFunctionsGradients( Matrix& rResult, CoordinatesArrayType& rPoint )
    {
        if(rResult.size1() != 2 || rResult.size2() != 1)
        {
            rResult.resize( 2, 1, false );
        }
        noalias( rResult ) = ZeroMatrix( 2, 1 );

        rResult( 0, 0 ) = - 0.5;
        rResult( 1, 0 ) =   0.5;
        return rResult;
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
        ) override
    {
        PointLocalCoordinates( rResult, rPoint );

        if ( std::abs( rResult[0] ) <= (1.0 + Tolerance) ) {
            return true;
        }

        return false;
    }

    /** Test the intersection with another geometry
     *  Test if this geometry intersects with other line_2d_2
     *
     * @param  rOtherGeometry Geometry to intersect with
     * @return True if the geometries intersect, False in any other case.
     */
    bool HasIntersection(const BaseType& rOtherGeometry) override
    {
        const double tolerance = std::numeric_limits<double>::epsilon();
        // We get the local points
        const TPointType& first_point  = BaseType::GetPoint(0); //p1
        const TPointType& second_point = BaseType::GetPoint(1); //p2

        // We get the other line's points
        const TPointType& first_point_other  = *rOtherGeometry(0); //p3
        const TPointType& second_point_other = *rOtherGeometry(1); //p4

        // parametric coordinate of intersection on current line
        const double numerator   = ( (first_point[0]-first_point_other[0])*(first_point_other[1] - second_point_other[1]) - (first_point[1]-first_point_other[1])*(first_point_other[0]-second_point_other[0]) );
        const double denominator = ( (first_point[0]-second_point[0])*(first_point_other[1] - second_point_other[1]) - (first_point[1]-second_point[1])*(first_point_other[0]-second_point_other[0]) );
        if (std::abs(denominator) < tolerance) // this means parallel lines.
            return false;
        const double t = numerator  /  denominator;

        return (0.0-tolerance<=t) && (t<=1.0+tolerance);
    }

    /** Test intersection of the geometry with a box (AABB)
     * Tests the intersection of the geometry with
     * a 3D box defined by rLowPoint and rHighPoint
     *
     * @param  rLowPoint  Lower point of the box to test the intersection
     * @param  rHighPoint Higher point of the box to test the intersection
     * @return            True if the geometry intersects the box, False in any other case.
     */
    bool HasIntersection(const Point& rLowPoint, const Point& rHighPoint) override
    {
        const double tolerance = std::numeric_limits<double>::epsilon();
        // We get the local points
        const TPointType& first_point  = BaseType::GetPoint(0);
        const TPointType& second_point = BaseType::GetPoint(1);

        if (    // If one of the point is inside the box then there is an intersection. If none is inside then we check further.
                ( (first_point[0] >= rLowPoint[0] && first_point[0] <= rHighPoint[0])
                    && (first_point[1] >= rLowPoint[1] && first_point[1] <= rHighPoint[1]) ) // IF the first point is inside the box
                ||
                  ( (second_point[0] >= rLowPoint[0] && second_point[0] <= rHighPoint[0])
                    && (second_point[1] >= rLowPoint[1] && second_point[1] <= rHighPoint[1]) ) // IF the second point is inside the box
            )
            return true;

        const double high_x = rHighPoint[0];
        const double high_y = rHighPoint[1];
        const double low_x = rLowPoint[0];
        const double low_y = rLowPoint[1];

        const double denominator = ( second_point[0] - first_point[0] );
        const double numerator = (second_point[1] - first_point[1]);
        const double slope = std::abs(denominator) > tolerance ? std::abs(numerator) > tolerance ? numerator / denominator : 1.0e-12 : 1.0e12;

        // Intersection with left vertical line of the box that is x = low_x
        const double y_1 = slope*( low_x - first_point[0] ) + first_point[1];
        if(y_1 >= low_y - tolerance && y_1 <= high_y+tolerance) // If y intersection is between two y bounds there is an intersection
            return true;
        // Intersection with right vertical line of the box that is x = high_x
        const double y_2 = slope*( high_x - first_point[0] ) + first_point[1];
        if(y_2 >= low_y - tolerance && y_2 <= high_y+tolerance) // If y intersection is between two y bounds there is an intersection
            return true;
        // Intersection with bottom horizontal line of the box that is y = low_y
        const double x_1 = first_point[0] + ( (low_y - first_point[1]) / slope );
        if(x_1 >= low_x-tolerance && x_1 <= high_x+tolerance) // If x intersection is between two x bounds there is an intersection
            return true;
        // Intersection with top horizontal line of the box that is y = high_y
        const double x_2 = first_point[0] + ( (high_y - first_point[1]) / slope );
        if(x_2 >= low_x-tolerance && x_2 <= high_x+tolerance) // If x intersection is between two x bounds there is an intersection
            return true;


        return false;
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

        const TPointType& FirstPoint  = BaseType::GetPoint(0);
        const TPointType& SecondPoint = BaseType::GetPoint(1);

        // Project point
        const double Tolerance = 1e-14; // Tolerance

        const double Length = std::sqrt((SecondPoint[0] - FirstPoint[0]) * (SecondPoint[0] - FirstPoint[0])
                    + (SecondPoint[1] - FirstPoint[1]) * (SecondPoint[1] - FirstPoint[1]));

        const double Length1 = std::sqrt((rPoint[0] - FirstPoint[0]) * (rPoint[0] - FirstPoint[0])
                    + (rPoint[1] - FirstPoint[1]) * (rPoint[1] - FirstPoint[1]));

        const double Length2 = std::sqrt((rPoint[0] - SecondPoint[0]) * (rPoint[0] - SecondPoint[0])
                    + (rPoint[1] - SecondPoint[1]) * (rPoint[1] - SecondPoint[1]));

        if (Length1 <= (Length + Tolerance) && Length2 <= (Length + Tolerance))
        {
            rResult[0] = 2.0 * Length1/(Length + Tolerance) - 1.0;
        }
        else if (Length1 > (Length + Tolerance))
        {
            rResult[0] = 2.0 * Length1/(Length + Tolerance) - 1.0; // NOTE: The same value as before, but it will be > than 1
        }
        else if (Length2 > (Length + Tolerance))
        {
            rResult[0] = 1.0 - 2.0 * Length2/(Length + Tolerance);
        }
        else
        {
            rResult[0] = 2.0; // Out of the line!!!
        }

        return( rResult );
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

    Line2D2(): BaseType( PointsArrayType(), &msGeometryData ) {}

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    static Matrix CalculateShapeFunctionsIntegrationPointsValues( typename BaseType::IntegrationMethod ThisMethod )
    {
        const IntegrationPointsContainerType& all_integration_points = AllIntegrationPoints();
        const IntegrationPointsArrayType& IntegrationPoints = all_integration_points[ThisMethod];
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
        const IntegrationPointsArrayType& IntegrationPoints = all_integration_points[ThisMethod];
        ShapeFunctionsGradientsType DN_De( IntegrationPoints.size() );
        std::fill( DN_De.begin(), DN_De.end(), Matrix( 2, 1 ) );

        for ( unsigned int it_gp = 0; it_gp < IntegrationPoints.size(); it_gp++ )
        {
            // double e = IntegrationPoints[it_gp].X();
            DN_De[it_gp]( 0, 0 ) = -0.5;
            DN_De[it_gp]( 1, 0 ) =  0.5;
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
                Line2D2<TPointType>::CalculateShapeFunctionsIntegrationPointsValues( GeometryData::GI_GAUSS_1 ),
                Line2D2<TPointType>::CalculateShapeFunctionsIntegrationPointsValues( GeometryData::GI_GAUSS_2 ),
                Line2D2<TPointType>::CalculateShapeFunctionsIntegrationPointsValues( GeometryData::GI_GAUSS_3 ),
                Line2D2<TPointType>::CalculateShapeFunctionsIntegrationPointsValues( GeometryData::GI_GAUSS_4 ),
                Line2D2<TPointType>::CalculateShapeFunctionsIntegrationPointsValues( GeometryData::GI_GAUSS_5 ),
                Line2D2<TPointType>::CalculateShapeFunctionsIntegrationPointsValues( GeometryData::GI_EXTENDED_GAUSS_1 ),
                Line2D2<TPointType>::CalculateShapeFunctionsIntegrationPointsValues( GeometryData::GI_EXTENDED_GAUSS_2 ),
                Line2D2<TPointType>::CalculateShapeFunctionsIntegrationPointsValues( GeometryData::GI_EXTENDED_GAUSS_3 ),
                Line2D2<TPointType>::CalculateShapeFunctionsIntegrationPointsValues( GeometryData::GI_EXTENDED_GAUSS_4 ),
                Line2D2<TPointType>::CalculateShapeFunctionsIntegrationPointsValues( GeometryData::GI_EXTENDED_GAUSS_5 )
            }
        };
        return shape_functions_values;
    }

    static const ShapeFunctionsLocalGradientsContainerType AllShapeFunctionsLocalGradients()
    {
        ShapeFunctionsLocalGradientsContainerType shape_functions_local_gradients = {{
                Line2D2<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients( GeometryData::GI_GAUSS_1 ),
                Line2D2<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients( GeometryData::GI_GAUSS_2 ),
                Line2D2<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients( GeometryData::GI_GAUSS_3 ),
                Line2D2<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients( GeometryData::GI_GAUSS_4 ),
                Line2D2<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients( GeometryData::GI_GAUSS_5 ),
                Line2D2<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients( GeometryData::GI_EXTENDED_GAUSS_1 ),
                Line2D2<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients( GeometryData::GI_EXTENDED_GAUSS_2 ),
                Line2D2<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients( GeometryData::GI_EXTENDED_GAUSS_3 ),
                Line2D2<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients( GeometryData::GI_EXTENDED_GAUSS_4 ),
                Line2D2<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients( GeometryData::GI_EXTENDED_GAUSS_5 )
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

    template<class TOtherPointType> friend class Line2D2;

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
                                   Line2D2<TPointType>& rThis );

/// output stream function
template<class TPointType>
inline std::ostream& operator << ( std::ostream& rOStream,
                                   const Line2D2<TPointType>& rThis )
{
    rThis.PrintInfo( rOStream );
    rOStream << std::endl;
    rThis.PrintData( rOStream );

    return rOStream;
}

///@}


template<class TPointType>
const GeometryData Line2D2<TPointType>::msGeometryData( 2,
        2,
        1,
        GeometryData::GI_GAUSS_1,
        Line2D2<TPointType>::AllIntegrationPoints(),
        Line2D2<TPointType>::AllShapeFunctionsValues(),
        AllShapeFunctionsLocalGradients() );

}  // namespace Kratos.

#endif // KRATOS_LINE_2D_H_INCLUDED  defined
