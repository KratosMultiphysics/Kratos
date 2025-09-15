//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Klaus B. Sautter


#if !defined(KRATOS_LINE_3D_N_H_INCLUDED )
#define  KRATOS_LINE_3D_N_H_INCLUDED

// System includes

// External includes

// Project includes
#include "geometries/geometry.h"
#include "integration/line_gauss_legendre_integration_points.h"


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
 * @class Line3DN
 * @ingroup KratosCore
 * @brief A arbitrary node 3D line geometry with quadratic shape functions
 */
template<class TPointType>

class Line3DN : public Geometry<TPointType>
{
public:
    ///@}
    ///@name Type Definitions
    ///@{

    /// Geometry as base class.
    typedef Geometry<TPointType> BaseType;

    /// Pointer definition of Line3DN
    KRATOS_CLASS_POINTER_DEFINITION( Line3DN );

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

//     Line3DN( const PointType& FirstPoint, const PointType& SecondPoint, const PointType& ThirdPoint )
//         : BaseType( PointsArrayType(), &msGeometryData )
//     {
//         BaseType::Points().push_back( typename PointType::Pointer( new PointType( FirstPoint ) ) );
//         BaseType::Points().push_back( typename PointType::Pointer( new PointType( SecondPoint ) ) );
//         BaseType::Points().push_back( typename PointType::Pointer( new PointType( ThirdPoint ) ) );
//     }

    /* Line3DN( typename PointType::Pointer pFirstPoint, typename PointType::Pointer pSecondPoint,
             typename PointType::Pointer pThirdPoint )
        : BaseType( PointsArrayType(), &msGeometryData )
    {
        BaseType::Points().push_back( pFirstPoint );
        BaseType::Points().push_back( pSecondPoint );
        BaseType::Points().push_back( pThirdPoint );
    } */

    Line3DN( const PointsArrayType& ThisPoints )
        : BaseType( ThisPoints, &msGeometryData )
    {
        //if ( BaseType::PointsNumber() != 3 )
        //    KRATOS_ERROR << "Invalid points number. Expected 3, given " << BaseType::PointsNumber() << std::endl;
    }

    /** Copy constructor.
     * Construct this geometry as a copy of given geometry.
     * @note This copy constructor don't copy the points and new
     * geometry shares points with given source geometry. It's
     * obvious that any change to this new geometry's point affect
     * source geometry's points too.
     */
    Line3DN( Line3DN const& rOther )
        : BaseType( rOther )
    {
    }

    /** Copy constructor from a geometry with other point type.
     * Construct this geometry as a copy of given geometry which
     * has different type of points. The given goemetry's
     * TOtherPointType* must be implicity convertible to this
     * geometry PointType.
     * @note This copy constructor don't copy the points and new
     * geometry shares points with given source geometry. It's
     * obvious that any change to this new geometry's point affect
     * source geometry's points too.
     */
    template<class TOtherPointType> Line3DN( Line3DN<TOtherPointType> const& rOther )
        : BaseType( rOther )
    {
    }

    /// Destructor. Do nothing!!!
    ~Line3DN() override {}

/*     GeometryData::KratosGeometryFamily GetGeometryFamily() const override
    {
        return GeometryData::Kratos_Linear;
    }

    GeometryData::KratosGeometryType GetGeometryType() const override
    {
        return GeometryData::Kratos_Line3DN;
    } */

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
    Line3DN& operator=( const Line3DN& rOther )
    {
        BaseType::operator=( rOther );
        return *this;
    }

    /** Assignment operator for geometries with different point type.
     * @note This operator don't copy the points and this
     * geometry shares points with given source geometry. It's
     * obvious that any change to this geometry's point affect
     * source geometry's points too.
     * @see Clone
     * @see ClonePoints
     */
    template<class TOtherPointType>
    Line3DN& operator=( Line3DN<TOtherPointType> const & rOther )
    {
        BaseType::operator=( rOther );
        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    typename BaseType::Pointer Create( PointsArrayType const& ThisPoints ) const override
    {
      // line 3 conectivities order 1-3-2
      return typename BaseType::Pointer( new Line3DN( ThisPoints ) );
    }

    // Geometry< Point<3> >::Pointer Clone() const override
    // {
    //     Geometry< Point<3> >::PointsArrayType NewPoints;

    //     //making a copy of the nodes TO POINTS (not Nodes!!!)
    //     for ( IndexType i = 0 ; i < this->size() ; i++ )
    //     {
    //             NewPoints.push_back(Kratos::make_shared< Point<3> >(( *this )[i]));
    //     }

    //     //creating a geometry with the new points
    //     Geometry< Point<3> >::Pointer p_clone( new Line3DN< Point<3> >( NewPoints ) );

    //     return p_clone;
    // }

    //lumping factors for the calculation of the lumped mass matrix
    Vector& LumpingFactors( Vector& rResult,
        const typename BaseType::LumpingMethods LumpingMethod) const override
    {
        KRATOS_ERROR << "'LumpingFactors' not available for arbitrarty noded line" << std::endl;
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
        const int number_points  = BaseType::PointsNumber();
        const TPointType& point0 = BaseType::GetPoint(0);
        const TPointType& point1 = BaseType::GetPoint(number_points-1);
        const double lx = point0.X() - point1.X();
        const double ly = point0.Y() - point1.Y();
        const double lz = point0.Z() - point1.Z();

        const double length = lx * lx + ly * ly + lz * lz;

        return sqrt( length );
    }

    /** This method calculate and return area or surface area of
    this geometry depending to it's dimension. For one dimensional
    geometry it returns lentgh, for two dimensional it gives area
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
        KRATOS_ERROR << "'DomainSize' not available for arbitrarty noded line" << std::endl;
        return 0.0;
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
        KRATOS_ERROR << "'IsInside' not available for arbitrarty noded line" << std::endl;
        return false;
    }

    ///@}
    ///@name Jacobian
    ///@{

    /** Jacobians for given  method. This method
     * calculate jacobians matrices in all integrations points of
     * given integration method.
     *
     * @param ThisMethod integration method which jacobians has to
     * be calculated in its integration points.
     *
     * @return JacobiansType a Vector of jacobian
     * matrices \f$ J_i \f$ where \f$ i=1,2,...,n \f$ is the integration
     * point index of given integration method.
     *
     * @see DeterminantOfJacobian
     * @see InverseOfJacobian
     */
    JacobiansType& Jacobian( JacobiansType& rResult, IntegrationMethod ThisMethod ) const override
    {
        KRATOS_ERROR << "'Jacobian' not available for arbitrarty noded line" << std::endl;
        return rResult;
    }


   /** Jacobians for given  method. This method
     * calculate jacobians matrices in all integrations points of
     * given integration method.
     *
     * @param ThisMethod integration method which jacobians has to
     * be calculated in its integration points.
     *
     * @return JacobiansType a Vector of jacobian
     * matrices \f$ J_i \f$ where \f$ i=1,2,...,n \f$ is the integration
     * point index of given integration method.
     *
     * @param DeltaPosition Matrix with the nodes position increment which describes
     * the configuration where the jacobian has to be calculated.
     *
     * @see DeterminantOfJacobian
     * @see InverseOfJacobian
     */
    JacobiansType& Jacobian( JacobiansType& rResult, IntegrationMethod ThisMethod, Matrix & DeltaPosition ) const override
    {
        KRATOS_ERROR << "'Jacobian' not available for arbitrarty noded line" << std::endl;
        return rResult;
    }

    /** Jacobian in specific integration point of given integration
     * method. This method calculate jacobian matrix in given
     * integration point of given integration method.
     *
     * @param IntegrationPointIndex index of integration point which jacobians has to
     * be calculated in it.
     *
     * @param ThisMethod integration method which jacobians has to
     * be calculated in its integration points.
     *
     * @return Matrix<double> Jacobian matrix \f$ J_i \f$ where \f$
     * i \f$ is the given integration point index of given
     * integration method.
     * @see DeterminantOfJacobian
     * @see InverseOfJacobian
     */
    Matrix& Jacobian( Matrix& rResult, IndexType IntegrationPointIndex, IntegrationMethod ThisMethod ) const override
    {
        KRATOS_ERROR << "'Jacobian' not available for arbitrarty noded line" << std::endl;
        return rResult;
    }

    /** Jacobian in given point. This method calculate jacobian
     * matrix in given point.
     *
     * @param rPoint point which jacobians has to
     * be calculated in it.
     *
     * @return Matrix of double which is jacobian matrix \f$ J \f$ in given point.
     *
     * @see DeterminantOfJacobian
     * @see InverseOfJacobian
     */
    Matrix& Jacobian( Matrix& rResult, const CoordinatesArrayType& rPoint ) const override
    {
        KRATOS_ERROR << "'Jacobian' not available for arbitrarty noded line" << std::endl;
        return rResult;
    }

    /** Determinant of jacobians for given integration method. This
     * method calculate determinant of jacobian in all
     * integrations points of given integration method.
     *
     * @return Vector of double which is vector of determinants of
     * jacobians \f$ |J|_i \f$ where \f$ i=1,2,...,n \f$ is the
     * integration point index of given integration method.
     *
     * @see Jacobian
     * @see InverseOfJacobian
     */
    Vector& DeterminantOfJacobian( Vector& rResult, IntegrationMethod ThisMethod ) const override
    {
        KRATOS_ERROR << "Jacobian is not square" << std::endl;
        return rResult;
    }

    /** Determinant of jacobian in specific integration point of
     * given integration method. This method calculate determinant
     * of jacobian in given integration point of given integration
     * method.
     *
     * @param IntegrationPointIndex index of integration point which jacobians has to
     * be calculated in it.
     *
     * @param IntegrationPointIndex index of integration point
     * which determinant of jacobians has to be calculated in it.
     *
     * @return Determinamt of jacobian matrix \f$ |J|_i \f$ where \f$
     * i \f$ is the given integration point index of given
     * integration method.

    * Inverse of jacobians for given integration method. This method
     * calculate inverse of jacobians matrices in all integrations points of
     * given integration method.
     *
     * @param ThisMethod integration method which inverse of jacobians has to
     * be calculated in its integration points.
     *
     * @return Inverse of jacobian
     * matrices \f$ J^{-1}_i \f$ where \f$ i=1,2,...,n \f$ is the integration
     * point index of given integration method.
     *
     * @see Jacobian
     * @see DeterminantOfJacobian
     */
    JacobiansType& InverseOfJacobian( JacobiansType& rResult, IntegrationMethod ThisMethod ) const override
    {
        KRATOS_ERROR << "Jacobian is not square" << std::endl;
        return rResult;
    }

    /** Inverse of jacobian in specific integration point of given integration
     * method. This method calculate Inverse of jacobian matrix in given
     * integration point of given integration method.
     *
     * @param IntegrationPointIndex index of integration point which inverse of jacobians has to
     * be calculated in it.
     *
     * @param ThisMethod integration method which inverse of jacobians has to
     * be calculated in its integration points.
     *
     * @return Inverse of jacobian matrix \f$ J^{-1}_i \f$ where \f$
     * i \f$ is the given integration point index of given
     * integration method.
     *
     * @see Jacobian
     * @see DeterminantOfJacobian
     */
    Matrix& InverseOfJacobian( Matrix& rResult, IndexType IntegrationPointIndex,
                                       IntegrationMethod ThisMethod ) const override
    {
        KRATOS_ERROR << "Jacobian is not square" << std::endl;
        return rResult;
    }

    /** Inverse of jacobian in given point. This method calculate inverse of jacobian
     * matrix in given point.
     *
     * @param rPoint point which inverse of jacobians has to
     * be calculated in it.
     *
     * @return Inverse of jacobian matrix \f$ J^{-1} \f$ in given point.
     *
     * @see DeterminantOfJacobian
     * @see InverseOfJacobian
     */
    Matrix& InverseOfJacobian( Matrix& rResult, const CoordinatesArrayType& rPoint ) const override
    {
        KRATOS_ERROR << "Jacobian is not square" << std::endl;
        return rResult;
    }


    /** EdgesNumber
    @return SizeType contains number of this geometry edges.
    */
    SizeType EdgesNumber() const override
    {
        return 0;
    }

    SizeType FacesNumber() const override
    {
        return 0;
    }

    ///@}
    ///@name Shape Function
    ///@{

    double ShapeFunctionValue( IndexType ShapeFunctionIndex,
                                       const CoordinatesArrayType& rPoint ) const override
    {
        KRATOS_ERROR << "'ShapeFunctionValue' not available for arbitrarty noded line" << std::endl;
        return 0;
    }



    void ShapeFunctionsIntegrationPointsGradients( ShapeFunctionsGradientsType& rResult, IntegrationMethod ThisMethod ) const override
    {
        KRATOS_ERROR << "Jacobian is not square" << std::endl;
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
        return "1 dimensional line with n nodes in 3D space";
    }

    /** Print information about this object.

    @param rOStream Stream to print into it.
    @see PrintData()
    @see Info()
    */
    void PrintInfo( std::ostream& rOStream ) const override
    {
        rOStream << "1 dimensional line with n nodes in 3D space";
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
    }

    /**
     * Calculates the local gradients for all integration points for
     * given integration method
     */
    virtual ShapeFunctionsGradientsType ShapeFunctionsLocalGradients(
        IntegrationMethod ThisMethod )
    {
        KRATOS_ERROR << "'ShapeFunctionsLocalGradients' not available for arbitrarty noded line" << std::endl;
        ShapeFunctionsGradientsType Result( 0 );
        return Result;
    }

    /**
     * Calculates the local gradients for all integration points for the
     * default integration method
     */
    virtual ShapeFunctionsGradientsType ShapeFunctionsLocalGradients()
    {
        KRATOS_ERROR << "'ShapeFunctionsLocalGradients' not available for arbitrarty noded line" << std::endl;
        ShapeFunctionsGradientsType Result( 0 );
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
        KRATOS_ERROR << "'ShapeFunctionsLocalGradients' not available for arbitrarty noded line" << std::endl;
        return( rResult );
    }

    /**
     * returns the local coordinates of all nodes of the current geometry
     * @param rResult a Matrix object that will be overwritten by the result
     * @return the local coordinates of all nodes
     */
    Matrix& PointsLocalCoordinates( Matrix& rResult ) const override
    {
        KRATOS_ERROR << "'PointsLocalCoordinates' not available for arbitrarty noded line" << std::endl;
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
        KRATOS_ERROR << "'ShapeFunctionsGradients' not available for arbitrarty noded line" << std::endl;
        return rResult;
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

    Line3DN(): BaseType( PointsArrayType(), &msGeometryData ) {}


    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    static Matrix CalculateShapeFunctionsIntegrationPointsValues( typename BaseType::IntegrationMethod ThisMethod )
    {
        KRATOS_ERROR << "'CalculateShapeFunctionsIntegrationPointsValues' not available for arbitrarty noded line" << std::endl;
        Matrix N( 0, 3 );
        return N;
    }

    static ShapeFunctionsGradientsType CalculateShapeFunctionsIntegrationPointsLocalGradients(
        typename BaseType::IntegrationMethod ThisMethod )
    {
        //KRATOS_ERROR << "'CalculateShapeFunctionsIntegrationPointsLocalGradients' not available for arbitrarty noded line" << std::endl;
        ShapeFunctionsGradientsType DN_De( 0 );
        return DN_De;
    }

    static const IntegrationPointsContainerType AllIntegrationPoints()
    {
        //KRATOS_ERROR << "'AllIntegrationPoints' not available for arbitrarty noded line" << std::endl;
        IntegrationPointsContainerType integration_points;
        return integration_points;
    }

    static const ShapeFunctionsValuesContainerType AllShapeFunctionsValues()
    {
        //KRATOS_ERROR << "'AllIntegrationPoints' not available for arbitrarty noded line" << std::endl;
        ShapeFunctionsValuesContainerType shape_functions_values;
        return shape_functions_values;
    }

    static const ShapeFunctionsLocalGradientsContainerType AllShapeFunctionsLocalGradients()
    {
        //KRATOS_ERROR << "'AllShapeFunctionsLocalGradients' not available for arbitrarty noded line" << std::endl;
        ShapeFunctionsLocalGradientsContainerType shape_functions_local_gradients;
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

    template<class TOtherPointType> friend class Line3DN;

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
                                   Line3DN<TPointType>& rThis );

/// output stream function
template<class TPointType>
inline std::ostream& operator << ( std::ostream& rOStream,
                                   const Line3DN<TPointType>& rThis )
{
    rThis.PrintInfo( rOStream );
    rOStream << std::endl;
    rThis.PrintData( rOStream );

    return rOStream;
}

///@}


template<class TPointType>
const GeometryData Line3DN<TPointType>::msGeometryData(
        &msGeometryDimension,
        GeometryData::IntegrationMethod::GI_GAUSS_2,
        Line3DN<TPointType>::AllIntegrationPoints(),
        Line3DN<TPointType>::AllShapeFunctionsValues(),
        AllShapeFunctionsLocalGradients() );


template<class TPointType>
const GeometryDimension Line3DN<TPointType>::msGeometryDimension(3, 1);

}  // namespace Kratos.

#endif // KRATOS_LINE_3D_3_H_INCLUDED  defined
