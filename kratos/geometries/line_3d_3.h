//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
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
#include "geometries/line_3d_2.h"
#include "integration/line_gauss_legendre_integration_points.h"
#include "integration/line_collocation_integration_points.h"
#include "utilities/integration_utilities.h"
#include "utilities/polynomial_utilities.h"

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
 * @class Line3D3
 * @ingroup KratosCore
 * @brief An three node 3D line geometry with quadratic shape functions
 * @details The node ordering corresponds with:
 *      0-----2----1
 * @author Riccardo Rossi
 * @author Janosch Stascheit
 * @author Felix Nagel
 */
template<class TPointType>

class Line3D3 : public Geometry<TPointType>
{
public:
    ///@name Type Definitions
    ///@{

    /// Geometry as base class.
    typedef Geometry<TPointType> BaseType;

    /// Pointer definition of Line3D3
    KRATOS_CLASS_POINTER_DEFINITION( Line3D3 );

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
    container.
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

    /** Type of the normal vector used for normal to edges in geometry.
     */
    typedef typename BaseType::NormalType NormalType;

    /**
     * Type of coordinates array
     */
    typedef typename BaseType::CoordinatesArrayType CoordinatesArrayType;

    ///@}
    ///@name Life Cycle
    ///@{

    Line3D3( typename PointType::Pointer pFirstPoint, typename PointType::Pointer pSecondPoint,
             typename PointType::Pointer pThirdPoint )
        : BaseType( PointsArrayType(), &msGeometryData )
    {
        BaseType::Points().push_back( pFirstPoint );
        BaseType::Points().push_back( pSecondPoint );
        BaseType::Points().push_back( pThirdPoint );
    }

    Line3D3( const PointsArrayType& ThisPoints )
        : BaseType( ThisPoints, &msGeometryData )
    {
        if ( BaseType::PointsNumber() != 3 )
            KRATOS_ERROR << "Invalid points number. Expected 3, given " << BaseType::PointsNumber() << std::endl;
    }

    /// Constructor with Geometry Id
    explicit Line3D3(
        const IndexType GeometryId,
        const PointsArrayType& rThisPoints
    ) : BaseType(GeometryId, rThisPoints, &msGeometryData)
    {
        KRATOS_ERROR_IF( this->PointsNumber() != 3 ) << "Invalid points number. Expected 3, given " << this->PointsNumber() << std::endl;
    }

    /// Constructor with Geometry Name
    explicit Line3D3(
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
    Line3D3( Line3D3 const& rOther )
        : BaseType( rOther )
    {
    }

    /** Copy constructor from a geometry with other point type.
     * Construct this geometry as a copy of given geometry which
     * has different type of points. The given goemetry's
     * TOtherPointType* must be implicitly convertible to this
     * geometry PointType.
     * @note This copy constructor don't copy the points and new
     * geometry shares points with given source geometry. It's
     * obvious that any change to this new geometry's point affect
     * source geometry's points too.
     */
    template<class TOtherPointType> Line3D3( Line3D3<TOtherPointType> const& rOther )
        : BaseType( rOther )
    {
    }

    /// Destructor. Do nothing!!!
    ~Line3D3() override {}

    /**
     * @brief Gets the geometry family.
     * @details This function returns the family type of the geometry. The geometry family categorizes the geometry into a broader classification, aiding in its identification and processing.
     * @return GeometryData::KratosGeometryFamily The geometry family.
     */
    GeometryData::KratosGeometryFamily GetGeometryFamily() const override
    {
        return GeometryData::KratosGeometryFamily::Kratos_Linear;
    }

    /**
     * @brief Gets the geometry type.
     * @details This function returns the specific type of the geometry. The geometry type provides a more detailed classification of the geometry.
     * @return GeometryData::KratosGeometryType The specific geometry type.
     */
    GeometryData::KratosGeometryType GetGeometryType() const override
    {
        return GeometryData::KratosGeometryType::Kratos_Line3D3;
    }

    /**
     * @brief Gets the geometry order type.
     * @details This function returns the order type of the geometry. The order type relates to the polynomial degree of the geometry.
     * @return GeometryData::KratosGeometryOrderType The geometry order type.
     */
    GeometryData::KratosGeometryOrderType GetGeometryOrderType() const override
    {
        return GeometryData::KratosGeometryOrderType::Kratos_Linear_Order;
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
    Line3D3& operator=( const Line3D3& rOther )
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
    Line3D3& operator=( Line3D3<TOtherPointType> const & rOther )
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
        return typename BaseType::Pointer( new Line3D3(NewGeometryId, rThisPoints ) );
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
        auto p_geometry = typename BaseType::Pointer( new Line3D3( NewGeometryId, rGeometry.Points() ) );
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
        rResult[0] = 1.0/6.0;
        rResult[2] = 2.0/3.0;
        rResult[1] = 1.0/6.0;
        return rResult;
    }

    ///@}
    ///@name Informations
    ///@{

    /** This method calculates and returns Length or charactereistic
    length of this geometry depending on its dimension. For one
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
        const IntegrationMethod integration_method = IntegrationUtilities::GetIntegrationMethodForExactMassMatrixEvaluation(*this);
        return IntegrationUtilities::ComputeDomainSize(*this, integration_method);
    }

    /** This method calculates and returns area or surface area of
    this geometry depending on its dimension. For one dimensional
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


    /** This method calculates and returns length, area or volume of
    this geometry depending on its dimension. For one dimensional
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
        this->PointLocalCoordinates( rResult, rPoint );

        if ( std::abs( rResult[0] ) <= (1.0 + Tolerance) )
            return true;

        return false;
    }

    /**
     * @brief Determinant of jacobians for given integration method.
     * @details This method calculates determinant of jacobian in all integrations points of given integration method.
     * @return Vector of double which is vector of determinants of jacobians \f$ |J|_i \f$ where \f$ i=1,2,...,n \f$ is the integration point index of given integration method.
     * @see Jacobian
     * @see InverseOfJacobian
     */
    Vector& DeterminantOfJacobian( Vector& rResult, IntegrationMethod ThisMethod ) const override
    {
        const std::size_t number_of_integration_points = this->IntegrationPointsNumber( ThisMethod );
        if( rResult.size() != number_of_integration_points)
            rResult.resize( number_of_integration_points, false );

        Matrix J(3, 1);
        for (std::size_t pnt = 0; pnt < number_of_integration_points; ++pnt) {
            this->Jacobian( J, pnt, ThisMethod);
            rResult[pnt] = std::sqrt(std::pow(J(0,0), 2) + std::pow(J(1,0), 2) + std::pow(J(2,0), 2));
        }
        return rResult;
    }

    /**
     * @brief Determinant of jacobian in specific integration point of given integration method. This method calculates determinant of jacobian in given integration point of given integration method.
     * @param IntegrationPointIndex index of integration point which jacobians has to be calculated in it.
     * @param IntegrationPointIndex index of integration point which determinant of jacobians has to be calculated in it.
     * @return Determinamt of jacobian matrix \f$ |J|_i \f$ where \f$ i \f$ is the given integration point index of given integration method.
     * @see Jacobian
     * @see InverseOfJacobian
     */
    double DeterminantOfJacobian(
        IndexType IntegrationPointIndex,
        IntegrationMethod ThisMethod ) const override
    {
        Matrix J(3, 1);
        this->Jacobian( J, IntegrationPointIndex, ThisMethod);
        return std::sqrt(std::pow(J(0,0), 2) + std::pow(J(1,0), 2) + std::pow(J(2,0), 2));
    }

    /**
     * @brief Determinant of jacobian in given point. This method calculates determinant of jacobian matrix in given point.
     * @param rPoint point which determinant of jacobians has to be calculated in it.
     * @return Determinamt of jacobian matrix \f$ |J| \f$ in given point.
     * @see DeterminantOfJacobian
     * @see InverseOfJacobian
     */
    double DeterminantOfJacobian( const CoordinatesArrayType& rPoint ) const override
    {
        Matrix J(3, 1);
        this->Jacobian( J, rPoint);
        return std::sqrt(std::pow(J(0,0), 2) + std::pow(J(1,0), 2) + std::pow(J(2,0), 2));
    }

    ///@}
    ///@name Edges and faces
    ///@{

    /** EdgesNumber
    @return SizeType contains number of this geometry edges.
    */
    SizeType EdgesNumber() const override
    {
        return 2;
    }

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
    double ShapeFunctionValue( IndexType ShapeFunctionIndex,
                                       const CoordinatesArrayType& rPoint ) const override
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
        if(rResult.size1() != 3 || rResult.size2() != 1) {
            rResult.resize( 3, 1, false );
        }
        noalias( rResult ) = ZeroMatrix( 3, 1 );
        rResult( 0, 0 ) = -1.0;
        rResult( 1, 0 ) =  1.0;
        rResult( 2, 0 ) =  0.0;
        return rResult;
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
        // Define distance_objective as the gradient of ||point - global_coordinate(xi)||**2,
        // that is (point - local_coordinate(xi))*coordinate_derivative(xi) (times -2, which we ignore)
        // The zeros of the distance_objective are local extremes of the distance between point and line
        // We will get the zeros contained in the interval [-1, 1] (if any). If, for any of these zeros,
        // the point is on the line, the local coordinate is our result.
        rResult.clear();

        constexpr double TOLERANCE = 1e-12;

        const auto& r_p0 = this->GetPoint(0);
        const auto& r_p1 = this->GetPoint(1);
        const auto& r_p2 = this->GetPoint(2);

        // The method is prone to numerical instability close to the ends of the search interval
        // exit early in that case
        const array_1d<double,3> d0 = r_p0 - rPoint;
        if (MathUtils<double>::Dot3(d0, d0) < TOLERANCE) {
            rResult[0] = -1.0;
            return rResult;
        }
        const array_1d<double,3> d1 = r_p1 - rPoint;
        if (MathUtils<double>::Dot3(d1, d1) < TOLERANCE) {
            rResult[0] = 1.0;
            return rResult;
        }

        array_1d<double, 3> c1 = r_p0 + r_p1 - 2*r_p2;
        array_1d<double, 3> c2 = r_p1 - r_p0;
        array_1d<double, 3> c3 = r_p2 - rPoint;

        const double aux1 = MathUtils<double>::Dot3(c1, c1);
        if (aux1 < TOLERANCE) {
            // The geometry is a straight line, fall back to Line3D2.h
            auto line = Line3D2<TPointType>(
                this->pGetPoint(0), this->pGetPoint(1));
            return line.PointLocalCoordinates(rResult, rPoint);
        }

        const double aux2 = MathUtils<double>::Dot3(c1, c3);
        if (std::abs(aux2) < TOLERANCE) {
            // r_p2 == rPoint (we got the center of the line),
            //the local coordinate is 0
            return rResult;
        }

        PolynomialUtilities::PolynomialType distance_objective{
            0.5 * aux1,
            0.75 * MathUtils<double>::Dot3(c1, c2),
            0.25 * MathUtils<double>::Dot3(c2, c2) + aux2,
            0.5 * MathUtils<double>::Dot3(c2, c3)
        };

        std::vector<PolynomialUtilities::IntervalType> root_ranges;
        PolynomialUtilities::IsolateRoots(
            root_ranges, distance_objective,
            PolynomialUtilities::IntervalType{-1,1});

        Vector shape;
        for (const auto& interval: root_ranges) {
            double root = PolynomialUtilities::FindRoot(distance_objective, interval);
            // if point == line(root) we found our coordinate;
            rResult[0] = root;
            this->ShapeFunctionsValues(shape, rResult);
            array_1d<double,3> d = (shape[0]*r_p0 + shape[1]*r_p1 + shape[2]*r_p2) - rPoint;
            if (MathUtils<double>::Dot3(d, d) < TOLERANCE) {
                return rResult;
            }
        }

        // No points in the interval [-1,1] correspond to our local coordinate
        rResult[0] = 2.0;
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
        if(rResult.size1() != 3 || rResult.size2() != 1) {
            rResult.resize( 3, 1, false );
        }
        noalias( rResult ) = ZeroMatrix( 3, 1 );

        rResult( 0, 0 ) =  rPoint[0] - 0.5;
        rResult( 1, 0 ) =  rPoint[0] + 0.5;
        rResult( 2, 0 ) = -rPoint[0] * 2.0;
        return rResult;
    }

    ///@}
    ///@name Input and output
    ///@{

    /// @copydoc Geometry::Name
    std::string Name() const override
    {
        return "Line3D3N";
    }

    /** Turn back information as a string.

    @return String contains information about this geometry.
    @see PrintData()
    @see PrintInfo()
    */
    std::string Info() const override
    {
        return "1 dimensional line with 3 nodes in 3D space";
    }

    /** Print information about this object.

    @param rOStream Stream to print into it.
    @see PrintData()
    @see Info()
    */
    void PrintInfo( std::ostream& rOStream ) const override
    {
        rOStream << "1 dimensional line with 3 nodes in 3D space";
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

    Line3D3(): BaseType( PointsArrayType(), &msGeometryData ) {}

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
        const IntegrationPointsArrayType& IntegrationPoints = all_integration_points[static_cast<int>(ThisMethod)];
        ShapeFunctionsGradientsType DN_De( IntegrationPoints.size() );
        std::fill( DN_De.begin(), DN_De.end(), Matrix( 3, 1 ) );

        for ( unsigned int it_gp = 0; it_gp < IntegrationPoints.size(); it_gp++ )
        {
            Matrix aux_mat = ZeroMatrix(3,1);
            const double e = IntegrationPoints[it_gp].X();
            aux_mat(0,0) = e - 0.5;
            aux_mat(2,0) = -2.0 * e;
            aux_mat(1,0) = e + 0.5;
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
                Line3D3<TPointType>::CalculateShapeFunctionsIntegrationPointsValues( GeometryData::IntegrationMethod::GI_GAUSS_1 ),
                Line3D3<TPointType>::CalculateShapeFunctionsIntegrationPointsValues( GeometryData::IntegrationMethod::GI_GAUSS_2 ),
                Line3D3<TPointType>::CalculateShapeFunctionsIntegrationPointsValues( GeometryData::IntegrationMethod::GI_GAUSS_3 ),
                Line3D3<TPointType>::CalculateShapeFunctionsIntegrationPointsValues( GeometryData::IntegrationMethod::GI_GAUSS_4 ),
                Line3D3<TPointType>::CalculateShapeFunctionsIntegrationPointsValues( GeometryData::IntegrationMethod::GI_GAUSS_5 ),
                Line3D3<TPointType>::CalculateShapeFunctionsIntegrationPointsValues( GeometryData::IntegrationMethod::GI_EXTENDED_GAUSS_1 ),
                Line3D3<TPointType>::CalculateShapeFunctionsIntegrationPointsValues( GeometryData::IntegrationMethod::GI_EXTENDED_GAUSS_2 ),
                Line3D3<TPointType>::CalculateShapeFunctionsIntegrationPointsValues( GeometryData::IntegrationMethod::GI_EXTENDED_GAUSS_3 ),
                Line3D3<TPointType>::CalculateShapeFunctionsIntegrationPointsValues( GeometryData::IntegrationMethod::GI_EXTENDED_GAUSS_4 ),
                Line3D3<TPointType>::CalculateShapeFunctionsIntegrationPointsValues( GeometryData::IntegrationMethod::GI_EXTENDED_GAUSS_5 )
            }
        };
        return shape_functions_values;
    }

    static const ShapeFunctionsLocalGradientsContainerType AllShapeFunctionsLocalGradients()
    {
        ShapeFunctionsLocalGradientsContainerType shape_functions_local_gradients = {{
                Line3D3<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients( GeometryData::IntegrationMethod::GI_GAUSS_1 ),
                Line3D3<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients( GeometryData::IntegrationMethod::GI_GAUSS_2 ),
                Line3D3<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients( GeometryData::IntegrationMethod::GI_GAUSS_3 ),
                Line3D3<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients( GeometryData::IntegrationMethod::GI_GAUSS_4 ),
                Line3D3<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients( GeometryData::IntegrationMethod::GI_GAUSS_5 ),
                Line3D3<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients( GeometryData::IntegrationMethod::GI_EXTENDED_GAUSS_1 ),
                Line3D3<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients( GeometryData::IntegrationMethod::GI_EXTENDED_GAUSS_2 ),
                Line3D3<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients( GeometryData::IntegrationMethod::GI_EXTENDED_GAUSS_3 ),
                Line3D3<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients( GeometryData::IntegrationMethod::GI_EXTENDED_GAUSS_4 ),
                Line3D3<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients( GeometryData::IntegrationMethod::GI_EXTENDED_GAUSS_5 )
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

    template<class TOtherPointType> friend class Line3D3;

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
                                   Line3D3<TPointType>& rThis );

/// output stream function
template<class TPointType>
inline std::ostream& operator << ( std::ostream& rOStream,
                                   const Line3D3<TPointType>& rThis )
{
    rThis.PrintInfo( rOStream );
    rOStream << std::endl;
    rThis.PrintData( rOStream );

    return rOStream;
}

///@}


template<class TPointType>
const GeometryData Line3D3<TPointType>::msGeometryData(
        &msGeometryDimension,
        GeometryData::IntegrationMethod::GI_GAUSS_2,
        Line3D3<TPointType>::AllIntegrationPoints(),
        Line3D3<TPointType>::AllShapeFunctionsValues(),
        AllShapeFunctionsLocalGradients() );

template<class TPointType>
const GeometryDimension Line3D3<TPointType>::msGeometryDimension(3, 1);

}  // namespace Kratos.
