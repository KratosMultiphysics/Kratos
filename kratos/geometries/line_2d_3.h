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

#if !defined(KRATOS_LINE_2D_3_H_INCLUDED )
#define  KRATOS_LINE_2D_3_H_INCLUDED

// System includes

// External includes

// Project includes
#include "geometries/geometry.h"
#include "integration/line_gauss_legendre_integration_points.h"
#include "utilities/integration_utilities.h"

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
 * @class Line2D3
 * @ingroup KratosCore
 * @brief An three node 2D line geometry with quadratic shape functions
 * @details The node ordering corresponds with:
 *      0-----2----1
 * @author Riccardo Rossi
 * @author Janosch Stascheit
 * @author Felix Nagel
 */
template<class TPointType>

class Line2D3 : public Geometry<TPointType>
{
public:
    ///@}
    ///@name Type Definitions
    ///@{

    /// Geometry as base class.
    typedef Geometry<TPointType> BaseType;

    /// Pointer definition of Line2D3
    KRATOS_CLASS_POINTER_DEFINITION( Line2D3 );

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

    Line2D3( const PointType& FirstPoint,
             const PointType& SecondPoint,
             const PointType& ThirdPoint )
        : BaseType( PointsArrayType(), &msGeometryData )
    {
        BaseType::Points().push_back( typename PointType::Pointer( new PointType( FirstPoint ) ) );
        BaseType::Points().push_back( typename PointType::Pointer( new PointType( SecondPoint ) ) );
        BaseType::Points().push_back( typename PointType::Pointer( new PointType( ThirdPoint ) ) );
    }

    Line2D3( typename PointType::Pointer pFirstPoint,
             typename PointType::Pointer pSecondPoint,
             typename PointType::Pointer pThirdPoint )
        : BaseType( PointsArrayType(), &msGeometryData )
    {
        BaseType::Points().push_back( pFirstPoint );
        BaseType::Points().push_back( pSecondPoint );
        BaseType::Points().push_back( pThirdPoint );
    }

    explicit Line2D3( const PointsArrayType& rThisPoints )
        : BaseType( rThisPoints, &msGeometryData )
    {
        KRATOS_ERROR_IF( BaseType::PointsNumber() != 3 ) << "Invalid points number. Expected 3, given " << BaseType::PointsNumber() << std::endl;
    }

    /// Constructor with Geometry Id
    explicit Line2D3(
        const IndexType GeometryId,
        const PointsArrayType& rThisPoints
        ) : BaseType( GeometryId, rThisPoints, &msGeometryData )
    {
        KRATOS_ERROR_IF( this->PointsNumber() != 3 ) << "Invalid points number. Expected 3, given " << this->PointsNumber() << std::endl;
    }

    /// Constructor with Geometry Name
    explicit Line2D3(
        const std::string& rGeometryName,
        const PointsArrayType& rThisPoints
    ) : BaseType(rGeometryName, rThisPoints, &msGeometryData)
    {
        KRATOS_ERROR_IF(this->PointsNumber() != 3) << "Invalid points number. Expected 3, given " << this->PointsNumber() << std::endl;
    }

    /** Copy constructor.
     * Construct this geometry as a copy of given geometry.
     * @note This copy constructor don't copy the points and new
     * geometry shares points with given source geometry. It's
     * obvious that any change to this new geometry's point affect
     * source geometry's points too.
     */
    Line2D3( Line2D3 const& rOther )
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
    template<class TOtherPointType> explicit Line2D3( Line2D3<TOtherPointType> const& rOther )
        : BaseType( rOther )
    {
    }

    /// Destructor. Do nothing!!!
    ~Line2D3() override {}

    GeometryData::KratosGeometryFamily GetGeometryFamily() const override
    {
        return GeometryData::Kratos_Linear;
    }

    GeometryData::KratosGeometryType GetGeometryType() const override
    {
        return GeometryData::Kratos_Line2D3;
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
    Line2D3& operator=( const Line2D3& rOther )
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
    Line2D3& operator=( Line2D3<TOtherPointType> const & rOther )
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
        return typename BaseType::Pointer( new Line2D3( NewGeometryId, rThisPoints ) );
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
        auto p_geometry = typename BaseType::Pointer( new Line2D3( NewGeometryId, rGeometry.Points() ) );
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
	if(rResult.size() != 3)
           rResult.resize( 3, false );
        rResult[0] = 0.25;
        rResult[2] = 0.5;
        rResult[1] = 0.25;
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
        Vector temp;
        const IntegrationMethod integration_method = IntegrationUtilities::GetIntegrationMethodForExactMassMatrixEvaluation(*this);
        this->DeterminantOfJacobian( temp, integration_method );
        const IntegrationPointsArrayType& r_integration_points = this->IntegrationPoints( integration_method );
        double length = 0.0;

        for (std::size_t i = 0; i < r_integration_points.size(); ++i) {
            length += temp[i] * r_integration_points[i].Weight();
        }

        return length;
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
        BoundedMatrix<double,3,3> X;
        BoundedMatrix<double,3,1> DN;
        for(IndexType i=0; i<this->size(); ++i) {
            const auto& r_node = this->GetPoint(i);
            X(0, i) = r_node.X();
            X(1, i) = r_node.Y();
            X(2, i) = r_node.Z();
        }

        static constexpr double MaxNormPointLocalCoordinates = 300.0;
        static constexpr std::size_t MaxIteratioNumberPointLocalCoordinates = 500;
        static constexpr double MaxTolerancePointLocalCoordinates = 1.0e-8;

        Matrix J = ZeroMatrix( 1, 1 );
        Matrix invJ = ZeroMatrix( 1, 1 );

        // Starting with xi = 0
        if (rResult.size() != 3)
            rResult.resize(3, false);
        rResult = ZeroVector( 3 );
        double delta_xi = 0.0;
        const array_1d<double, 3> zero_array = ZeroVector(3);
        array_1d<double, 3> current_global_coords;

        //Newton iteration:
        for ( IndexType k = 0; k < MaxIteratioNumberPointLocalCoordinates; k++ ) {
            noalias(current_global_coords) = zero_array;
            this->GlobalCoordinates( current_global_coords, rResult );

            noalias( current_global_coords ) = rPoint - current_global_coords;

            // Derivatives of shape functions
            Matrix shape_functions_gradients;
            shape_functions_gradients = ShapeFunctionsLocalGradients(shape_functions_gradients, rResult );
            noalias(DN) = prod(X, shape_functions_gradients);

            noalias(J) = prod(trans(DN), DN);
            const array_1d<double, 1> res = prod(trans(DN), current_global_coords);

            // The inverted jacobian matrix
            invJ(0, 0) = 1.0/J( 0, 0 );

            delta_xi = invJ(0, 0) * res[0];

            rResult[0] += delta_xi;

            if ( delta_xi > MaxNormPointLocalCoordinates ) {
                KRATOS_WARNING_IF("Line2D3", k > 0) << "detJ =\t" << J( 0, 0 ) << " DeltaX =\t" << delta_xi << " stopping calculation. Iteration:\t" << k << std::endl;
                break;
            }

            if ( delta_xi < MaxTolerancePointLocalCoordinates )
                break;
        }

        return rResult;
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
        // Getting derivatives of shape functions
        const ShapeFunctionsGradientsType shape_functions_gradients = CalculateShapeFunctionsIntegrationPointsLocalGradients( ThisMethod );

        const std::size_t number_of_integration_points = this->IntegrationPointsNumber( ThisMethod );
        if (rResult.size() !=  number_of_integration_points) {
            JacobiansType temp( number_of_integration_points );
            rResult.swap( temp );
        }

        // Loop over all integration points
        for (std::size_t pnt = 0; pnt < number_of_integration_points; ++pnt) {
            // Initializing jacobian matrix
            noalias(rResult[pnt]) = ZeroMatrix( 2, 1 );

            // Loop over all nodes
            for (std::size_t i = 0; i < this->PointsNumber(); ++i) {
                const auto& r_node = this->GetPoint(i);
                rResult[pnt](0, 0) += r_node.X() * shape_functions_gradients[pnt](i, 0);
                rResult[pnt](1, 0) += r_node.Y() * shape_functions_gradients[pnt](i, 0);
            }
        } // End of loop over all integration points

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
     * @param rDeltaPosition Matrix with the nodes position increment which describes
     * the configuration where the jacobian has to be calculated.
     *
     * @see DeterminantOfJacobian
     * @see InverseOfJacobian
     */
    JacobiansType& Jacobian( JacobiansType& rResult, IntegrationMethod ThisMethod, Matrix& rDeltaPosition ) const override
    {
        // Getting derivatives of shape functions
        ShapeFunctionsGradientsType shape_functions_gradients =  CalculateShapeFunctionsIntegrationPointsLocalGradients( ThisMethod );
        const std::size_t number_of_integration_points = this->IntegrationPointsNumber( ThisMethod );

        // Getting values of shape functions
        Matrix shape_functions_values = CalculateShapeFunctionsIntegrationPointsValues( ThisMethod );

        if ( rResult.size() != number_of_integration_points ) {
            JacobiansType temp( number_of_integration_points );
            rResult.swap( temp );
        }

        // Loop over all integration points
        for (std::size_t pnt = 0; pnt < number_of_integration_points; ++pnt ) {
            // Initializing jacobian matrix
            noalias(rResult[pnt]) = ZeroMatrix( 2, 1 );

            // Loop over all nodes
            for (std::size_t i = 0; i < this->PointsNumber(); ++i ) {
                const auto& r_node = this->GetPoint(i);
                rResult[pnt](0, 0) += (r_node.X() - rDeltaPosition(i,0)) * shape_functions_gradients[pnt](i, 0);
                rResult[pnt](1, 0) += (r_node.Y() - rDeltaPosition(i,1)) * shape_functions_gradients[pnt](i, 0);
            }
        }// End of loop over all integration points

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
        // Setting up size of jacobian matrix
        rResult.resize(2, 1, false);
        noalias(rResult) = ZeroMatrix( 2, 1 );

        // Derivatives of shape functions
        ShapeFunctionsGradientsType shape_functions_gradients = CalculateShapeFunctionsIntegrationPointsLocalGradients( ThisMethod );
        Matrix shape_function_gradient_in_integration_point = shape_functions_gradients( IntegrationPointIndex );

        // Values of shape functions in integration points
        DenseVector<double> ShapeFunctionsValuesInIntegrationPoint = ZeroVector( 3 );
        ShapeFunctionsValuesInIntegrationPoint = row( CalculateShapeFunctionsIntegrationPointsValues( ThisMethod ), IntegrationPointIndex );

        // Elements of jacobian matrix (e.g. J(1,1) = dX1/dXi1)
        // Loop over all nodes
        for (std::size_t i = 0; i < this->PointsNumber(); ++i ) {
            const auto& r_node = this->GetPoint(i);
            rResult(0, 0) += r_node.X() * shape_function_gradient_in_integration_point(i, 0);
            rResult(1, 0) += r_node.Y() * shape_function_gradient_in_integration_point(i, 0);
        }

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
        // Setting up size of jacobian matrix
        rResult.resize( 2, 1, false );
        noalias(rResult) = ZeroMatrix( 2, 1 );

        // Derivatives of shape functions
        Matrix shape_functions_gradients;
        shape_functions_gradients = ShapeFunctionsLocalGradients( shape_functions_gradients, rPoint );

        // Elements of jacobian matrix (e.g. J(1,1) = dX1/dXi1)
        // Loop over all nodes
        for (std::size_t i = 0; i < this->PointsNumber(); ++i) {
            const auto& r_node = this->GetPoint(i);
            rResult(0, 0) += r_node.X() * shape_functions_gradients(i, 0);
            rResult(1, 0) += r_node.Y() * shape_functions_gradients(i, 0);
        }

        return rResult;
    }

    /**
     * @brief Determinant of jacobians for given integration method.
     * @details This method calculate determinant of jacobian in all integrations points of given integration method.
     * @return Vector of double which is vector of determinants of jacobians \f$ |J|_i \f$ where \f$ i=1,2,...,n \f$ is the integration point index of given integration method.
     * @see Jacobian
     * @see InverseOfJacobian
     */
    Vector& DeterminantOfJacobian( Vector& rResult, IntegrationMethod ThisMethod ) const override
    {
        const std::size_t number_of_integration_points = this->IntegrationPointsNumber( ThisMethod );
        if( rResult.size() != number_of_integration_points)
            rResult.resize( number_of_integration_points, false );

        Matrix J(2, 1);
        for (std::size_t pnt = 0; pnt < number_of_integration_points; ++pnt) {
            this->Jacobian( J, pnt, ThisMethod);
            rResult[pnt] = std::sqrt(std::pow(J(0,0), 2) + std::pow(J(1,0), 2));
        }
        return rResult;
    }

    /**
     * @brief Determinant of jacobian in specific integration point of given integration method. This method calculate determinant of jacobian in given integration point of given integration method.
     * @param IntegrationPointIndex index of integration point which jacobians has to be calculated in it.
     * @param IntegrationPointIndex index of integration point which determinant of jacobians has to be calculated in it.
     * @return Determinamt of jacobian matrix \f$ |J|_i \f$ where \f$ i \f$ is the given integration point index of given integration method.
     * @see Jacobian
     * @see InverseOfJacobian
     */
    double DeterminantOfJacobian( IndexType IntegrationPointIndex, IntegrationMethod ThisMethod ) const override
    {
        Matrix J(2, 1);
        this->Jacobian( J, IntegrationPointIndex, ThisMethod);
        return std::sqrt(std::pow(J(0,0), 2) + std::pow(J(1,0), 2));
    }

    /**
     * @brief Determinant of jacobian in given point. This method calculate determinant of jacobian matrix in given point.
     * @param rPoint point which determinant of jacobians has to be calculated in it.
     * @return Determinamt of jacobian matrix \f$ |J| \f$ in given point.
     * @see DeterminantOfJacobian
     * @see InverseOfJacobian
     */
    double DeterminantOfJacobian( const CoordinatesArrayType& rPoint ) const override
    {
        Matrix J(2, 1);
        this->Jacobian( J, rPoint);
        return std::sqrt(std::pow(J(0,0), 2) + std::pow(J(1,0), 2));
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
     Vector& ShapeFunctionsValues(
        Vector& rResult,
        const CoordinatesArrayType& rCoordinates
        ) const override
    {
        if(rResult.size() != 3) {
            rResult.resize(3, false);
        }

        rResult[0] = 0.5 * (rCoordinates[0] - 1.0) * rCoordinates[0];
        rResult[1] = 0.5 * (rCoordinates[0] + 1.0) * rCoordinates[0];
        rResult[2] = 1.0 - rCoordinates[0] * rCoordinates[0];

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
    double ShapeFunctionValue(
        IndexType ShapeFunctionIndex,
        const CoordinatesArrayType& rPoint
        ) const override
    {
        switch ( ShapeFunctionIndex )
        {
        case 0:
            return( 0.5*( rPoint[0] - 1.0 )*rPoint[0] );
        case 1:
            return( 0.5*( rPoint[0] + 1.0 )*rPoint[0] );
        case 2:
            return( 1.0 -rPoint[0]*rPoint[0] );
        default:
            KRATOS_ERROR << "Wrong index of shape function!" << *this << std::endl;
        }

        return 0;
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
        return "1 dimensional line with 3 nodes in 2D space";
    }

    /** Print information about this object.

    @param rOStream Stream to print into it.
    @see PrintData()
    @see Info()
    */
    void PrintInfo( std::ostream& rOStream ) const override
    {
        rOStream << "1 dimensional line with 3 nodes in 2D space";
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
     * @param rPoint the current point at which the gradients are calculated
     * @return the gradients of all shape functions
     * \f$ \frac{\partial N^i}{\partial \xi_j} \f$
     */
    Matrix& ShapeFunctionsLocalGradients( Matrix& rResult,
            const CoordinatesArrayType& rPoint ) const override
    {
        // Setting up result matrix
        if(rResult.size1() != 3 || rResult.size2() != 1) {
            rResult.resize( 3, 1, false );
        }

        noalias( rResult ) = ZeroMatrix( 3, 1 );
        rResult( 0, 0 ) =  rPoint[0] - 0.5;
        rResult( 1, 0 ) =  rPoint[0] + 0.5;
        rResult( 2, 0 ) = -rPoint[0] * 2.0;
        return( rResult );
    }

    /**
     * returns the local coordinates of all nodes of the current geometry
     * @param rResult a Matrix object that will be overwritten by the result
     * @return the local coordinates of all nodes
     */
    Matrix& PointsLocalCoordinates( Matrix& rResult ) const override
    {
        if(rResult.size1() != 3 || rResult.size2() != 1)
        {
            rResult.resize( 3, 1, false );
        }
        noalias( rResult ) = ZeroMatrix( 3, 1 );
        rResult( 0, 0 ) = -1.0;
        rResult( 1, 0 ) =  1.0;
        rResult( 2, 0 ) =  0.0;
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
        if(rResult.size1() != 3 || rResult.size2() != 1)
        {
            rResult.resize( 3, 1, false );
        }
        noalias( rResult ) = ZeroMatrix( 3, 1 );

        rResult( 0, 0 ) =  rPoint[0] - 0.5;
        rResult( 1, 0 ) =  rPoint[0] + 0.5;
        rResult( 2, 0 ) = -rPoint[0] * 2.0;
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

    Line2D3(): BaseType( PointsArrayType(), &msGeometryData ) {}

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
        Matrix N( integration_points_number, 3 );

        for ( int it_gp = 0; it_gp < integration_points_number; it_gp++ )
        {
            double e = IntegrationPoints[it_gp].X();
            N( it_gp, 0 ) = 0.5 * ( e - 1 ) * e;
            N( it_gp, 2 ) = 1.0 - e * e;
            N( it_gp, 1 ) = 0.5 * ( 1 + e ) * e;
        }

        return N;
    }

    static ShapeFunctionsGradientsType CalculateShapeFunctionsIntegrationPointsLocalGradients(
        typename BaseType::IntegrationMethod ThisMethod )
    {
        const IntegrationPointsContainerType& all_integration_points = AllIntegrationPoints();
        const IntegrationPointsArrayType& IntegrationPoints = all_integration_points[ThisMethod];
        ShapeFunctionsGradientsType DN_De( IntegrationPoints.size() );
        std::fill( DN_De.begin(), DN_De.end(), Matrix( 3, 1 ) );

        for ( unsigned int it_gp = 0; it_gp < IntegrationPoints.size(); it_gp++ )
        {
            double e = IntegrationPoints[it_gp].X();
            DN_De[it_gp]( 0, 0 ) = e - 0.5;
            DN_De[it_gp]( 2, 0 ) = -2.0 * e;
            DN_De[it_gp]( 1, 0 ) = e + 0.5;
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
                Quadrature<LineGaussLegendreIntegrationPoints5, 1, IntegrationPoint<3> >::GenerateIntegrationPoints()
            }
        };
        return integration_points;
    }

    static const ShapeFunctionsValuesContainerType AllShapeFunctionsValues()
    {
        ShapeFunctionsValuesContainerType shape_functions_values = {{
                Line2D3<TPointType>::CalculateShapeFunctionsIntegrationPointsValues( GeometryData::GI_GAUSS_1 ),
                Line2D3<TPointType>::CalculateShapeFunctionsIntegrationPointsValues( GeometryData::GI_GAUSS_2 ),
                Line2D3<TPointType>::CalculateShapeFunctionsIntegrationPointsValues( GeometryData::GI_GAUSS_3 ),
                Line2D3<TPointType>::CalculateShapeFunctionsIntegrationPointsValues( GeometryData::GI_GAUSS_4 ),
                Line2D3<TPointType>::CalculateShapeFunctionsIntegrationPointsValues( GeometryData::GI_GAUSS_5 )
            }
        };
        return shape_functions_values;
    }

    static const ShapeFunctionsLocalGradientsContainerType AllShapeFunctionsLocalGradients()
    {
        ShapeFunctionsLocalGradientsContainerType shape_functions_local_gradients = {{
                Line2D3<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients( GeometryData::GI_GAUSS_1 ),
                Line2D3<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients( GeometryData::GI_GAUSS_2 ),
                Line2D3<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients( GeometryData::GI_GAUSS_3 ),
                Line2D3<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients( GeometryData::GI_GAUSS_4 ),
                Line2D3<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients( GeometryData::GI_GAUSS_5 )
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

    template<class TOtherPointType> friend class Line2D3;

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
                                   Line2D3<TPointType>& rThis );

/// output stream function
template<class TPointType>
inline std::ostream& operator << ( std::ostream& rOStream,
                                   const Line2D3<TPointType>& rThis )
{
    rThis.PrintInfo( rOStream );
    rOStream << std::endl;
    rThis.PrintData( rOStream );

    return rOStream;
}

///@}

template<class TPointType>
const GeometryData Line2D3<TPointType>::msGeometryData(
        &msGeometryDimension,
        GeometryData::GI_GAUSS_2,
        Line2D3<TPointType>::AllIntegrationPoints(),
        Line2D3<TPointType>::AllShapeFunctionsValues(),
        AllShapeFunctionsLocalGradients() );

template<class TPointType>
const GeometryDimension Line2D3<TPointType>::msGeometryDimension(
    2, 2, 1);

}  // namespace Kratos.

#endif // KRATOS_LINE_2D_3_H_INCLUDED  defined
