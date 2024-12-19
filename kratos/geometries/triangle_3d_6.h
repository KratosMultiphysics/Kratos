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
#include "geometries/line_3d_3.h"
#include "integration/triangle_gauss_legendre_integration_points.h"
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
 * @class Triangle3D6
 * @ingroup KratosCore
 * @brief A six node 3D triangular geometry with quadratic shape functions
 * @details While the shape functions are only defined in 2D it is possible to define an arbitrary orientation in space. Thus it can be used for defining surfaces on 3D elements.
 * The node ordering corresponds with:
 *          2
 *          |`\
 *          |  `\
 *          5    `4
 *          |      `\
 *          |        `\
 *          0-----3----1
 * @author Riccardo Rossi
 * @author Janosch Stascheit
 * @author Felix Nagel
 */
template<class TPointType> class Triangle3D6
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
    typedef Line3D3<TPointType> EdgeType;

    /**
     * Pointer definition of Triangle3D6
     */
    KRATOS_CLASS_POINTER_DEFINITION( Triangle3D6 );

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
     * Type of the normal vector used for normal to edges in geomety.
     */
    typedef typename BaseType::NormalType NormalType;

    ///@}
    ///@name Life Cycle
    ///@{

//     Triangle3D6( const PointType& FirstPoint,
//                  const PointType& SecondPoint,
//                  const PointType& ThirdPoint,
//                  const PointType& FourthPoint,
//                  const PointType& FifthPoint,
//                  const PointType& SixthPoint
//                )
//         : BaseType( PointsArrayType(), &msGeometryData )
//     {
//         this->Points().push_back( typename PointType::Pointer( new PointType( FirstPoint ) ) );
//         this->Points().push_back( typename PointType::Pointer( new PointType( SecondPoint ) ) );
//         this->Points().push_back( typename PointType::Pointer( new PointType( ThirdPoint ) ) );
//         this->Points().push_back( typename PointType::Pointer( new PointType( FourthPoint ) ) );
//         this->Points().push_back( typename PointType::Pointer( new PointType( FifthPoint ) ) );
//         this->Points().push_back( typename PointType::Pointer( new PointType( SixthPoint ) ) );
//     }

    Triangle3D6( typename PointType::Pointer pFirstPoint,
                 typename PointType::Pointer pSecondPoint,
                 typename PointType::Pointer pThirdPoint,
                 typename PointType::Pointer pFourthPoint,
                 typename PointType::Pointer pFifthPoint,
                 typename PointType::Pointer pSixthPoint
               )
        : BaseType( PointsArrayType(), &msGeometryData )
    {
        this->Points().push_back( pFirstPoint );
        this->Points().push_back( pSecondPoint );
        this->Points().push_back( pThirdPoint );
        this->Points().push_back( pFourthPoint );
        this->Points().push_back( pFifthPoint );
        this->Points().push_back( pSixthPoint );
    }

    explicit Triangle3D6(
        const PointsArrayType& ThisPoints
    ) : BaseType( ThisPoints, &msGeometryData )
    {
        KRATOS_ERROR_IF( this->PointsNumber() != 6 ) << "Invalid points number. Expected 6, given " << this->PointsNumber() << std::endl;
    }

    /// Constructor with Geometry Id
    explicit Triangle3D6(
        const IndexType GeometryId,
        const PointsArrayType& rThisPoints
    ) : BaseType(GeometryId, rThisPoints, &msGeometryData)
    {
        KRATOS_ERROR_IF( this->PointsNumber() != 6 ) << "Invalid points number. Expected 6, given " << this->PointsNumber() << std::endl;
    }

    /// Constructor with Geometry Name
    explicit Triangle3D6(
        const std::string& rGeometryName,
        const PointsArrayType& rThisPoints
    ) : BaseType(rGeometryName, rThisPoints, &msGeometryData)
    {
        KRATOS_ERROR_IF(this->PointsNumber() != 20) << "Invalid points number. Expected 20, given " << this->PointsNumber() << std::endl;
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
    Triangle3D6( Triangle3D6 const& rOther )
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
    template<class TOtherPointType> Triangle3D6( Triangle3D6<TOtherPointType> const& rOther )
        : BaseType( rOther )
    {
    }

    /**
     * Destructor. Does nothing!!!
     */
    ~Triangle3D6() override {}

    /**
     * @brief Gets the geometry family.
     * @details This function returns the family type of the geometry. The geometry family categorizes the geometry into a broader classification, aiding in its identification and processing.
     * @return GeometryData::KratosGeometryFamily The geometry family.
     */
    GeometryData::KratosGeometryFamily GetGeometryFamily() const override
    {
        return GeometryData::KratosGeometryFamily::Kratos_Triangle;
    }

    /**
     * @brief Gets the geometry type.
     * @details This function returns the specific type of the geometry. The geometry type provides a more detailed classification of the geometry.
     * @return GeometryData::KratosGeometryType The specific geometry type.
     */
    GeometryData::KratosGeometryType GetGeometryType() const override
    {
        return GeometryData::KratosGeometryType::Kratos_Triangle3D6;
    }

    /**
     * @brief Gets the geometry order type.
     * @details This function returns the order type of the geometry. The order type relates to the polynomial degree of the geometry.
     * @return GeometryData::KratosGeometryOrderType The geometry order type.
     */
    GeometryData::KratosGeometryOrderType GetGeometryOrderType() const override
    {
        return GeometryData::KratosGeometryOrderType::Kratos_Quadratic_Order;
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
    Triangle3D6& operator=( const Triangle3D6& rOther )
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
    Triangle3D6& operator=( Triangle3D6<TOtherPointType> const & rOther )
    {
        BaseType::operator=( rOther );
        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Creates a new geometry pointer
     * @param rThisPoints the nodes of the new geometry
     * @return Pointer to the new geometry
     */
    typename BaseType::Pointer Create(
        PointsArrayType const& rThisPoints
        ) const override
    {
        return typename BaseType::Pointer( new Triangle3D6( rThisPoints ) );
    }

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
        return typename BaseType::Pointer( new Triangle3D6( NewGeometryId, rThisPoints ) );
    }

    /**
     * @brief Creates a new geometry pointer
     * @param rGeometry reference to an existing geometry
     * @return Pointer to the new geometry
     */
    typename BaseType::Pointer Create(
        const BaseType& rGeometry
        ) const override
    {
        auto p_geometry = typename BaseType::Pointer( new Triangle3D6( rGeometry.Points() ) );
        p_geometry->SetData(rGeometry.GetData());
        return p_geometry;
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
        auto p_geometry = typename BaseType::Pointer( new Triangle3D6( NewGeometryId, rGeometry.Points() ) );
        p_geometry->SetData(rGeometry.GetData());
        return p_geometry;
    }

    /**
     * returns the local coordinates of all nodes of the current geometry
     * @param rResult a Matrix object that will be overwritten by the result
     * @return the local coordinates of all nodes
     */
    Matrix& PointsLocalCoordinates( Matrix& rResult ) const override
    {
        rResult.resize( 6, 2,false );
        noalias( rResult ) = ZeroMatrix( 6, 2 );
        rResult( 0, 0 ) =  0.0;
        rResult( 0, 1 ) =  0.0;
        rResult( 1, 0 ) =  1.0;
        rResult( 1, 1 ) =  0.0;
        rResult( 2, 0 ) =  0.0;
        rResult( 2, 1 ) =  1.0;
        rResult( 3, 0 ) =  0.5;
        rResult( 3, 1 ) =  0.0;
        rResult( 4, 0 ) =  0.5;
        rResult( 4, 1 ) =  0.5;
        rResult( 5, 0 ) =  0.0;
        rResult( 5, 1 ) =  0.5;
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
     * @todo Could be replaced by something more suitable
     * @return double value contains length or Characteristic
     * length
     * @see Area()
     * @see Volume()
     * @see DomainSize()
     */
    double Length() const override
    {
        return std::sqrt(std::abs(this->DeterminantOfJacobian( PointType() ) ) );
    }

    /**
     * @brief This method calculates and returns area or surface area of this geometry depending to it's dimension.
     * @details For one dimensional geometry it returns zero, for two dimensional it gives area
     * and for three dimensional geometries it gives surface area.
     * @return double value contains area or surface area
     * @see Length()
     * @see Volume()
     * @see DomainSize()
     */
    double Area() const override
    {
        const IntegrationMethod integration_method = msGeometryData.DefaultIntegrationMethod();
        return IntegrationUtilities::ComputeDomainSize(*this, integration_method);
    }

    // TODO: Code activated in June 2023
    // /**
    //  * @brief This method calculates and returns the volume of this geometry.
    //  * @return Error, the volume of a 2D geometry is not defined
    //  * @see Length()
    //  * @see Area()
    //  * @see Volume()
    //  */
    // double Volume() const override
    // {
    //     KRATOS_ERROR << "Triangle3D6:: Method not well defined. Replace with DomainSize() instead" << std::endl;
    //     return 0.0;
    // }

    /**
     * @brief This method calculates and returns length, area or volume of this geometry depending to it's dimension.
     * @details For one dimensional geometry it returns its length, for two dimensional it gives area and for three dimensional geometries it gives its volume.
     * @return double value contains length, area or volume.
     * @see Length()
     * @see Area()
     * @see Volume()
     */
    double DomainSize() const override
    {
        return Area();
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
        )const  override
    {
        this->PointLocalCoordinates( rResult, rPoint );

        if ( (rResult[0] >= (0.0-Tolerance)) && (rResult[0] <= (1.0+Tolerance)) )
        {
            if ( (rResult[1] >= (0.0-Tolerance)) && (rResult[1] <= (1.0+Tolerance)) )
            {
                if ( (rResult[0] + rResult[1]) <= (1.0+Tolerance) )
                {
                    return true;
                }
            }
        }

        return false;
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
    double ShapeFunctionValue( IndexType ShapeFunctionIndex,
                                       const CoordinatesArrayType& rPoint ) const override
    {

        double thirdCoord = 1 - rPoint[0] - rPoint[1];

        switch ( ShapeFunctionIndex )
        {
        case 0:
            return( thirdCoord*( 2*thirdCoord - 1 ) );
        case 1:
            return( rPoint[0]*( 2*rPoint[0] - 1 ) );
        case 2:
            return( rPoint[1]*( 2*rPoint[1] - 1 ) );
        case 3:
            return( 4*thirdCoord*rPoint[0] );
        case 4:
            return( 4*rPoint[0]*rPoint[1] );
        case 5:
            return( 4*rPoint[1]*thirdCoord );

        default:
            KRATOS_ERROR << "Wrong index of shape function!" << *this << std::endl;
        }

        return 0;
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
        return GeometryUtils::PointDistanceToTriangle3D(this->GetPoint(0), this->GetPoint(1), this->GetPoint(2), this->GetPoint(3), this->GetPoint(4), this->GetPoint(5), point);
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
    std::string Info() const override
    {
        return "2 dimensional triangle with six nodes in 3D space";
    }

    /**
     * Print information about this object.
     * @param rOStream Stream to print into it.
     * @see PrintData()
     * @see Info()
     */
    void PrintInfo( std::ostream& rOStream ) const override
    {
        rOStream << "2 dimensional triangle with six nodes in 3D space";
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
    void PrintData( std::ostream& rOStream ) const override
    {
        // Base Geometry class PrintData call
        BaseType::PrintData( rOStream );
        std::cout << std::endl;

        // If the geometry has valid points, calculate and output its data
        if (this->AllPointsAreValid()) {
            Matrix jacobian;
            this->Jacobian( jacobian, PointType() );
            rOStream << "    Jacobian in the origin\t : " << jacobian;
        }
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
        return 3;
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

        edges.push_back( Kratos::make_shared<EdgeType>( this->pGetPoint( 0 ), this->pGetPoint( 1 ), this->pGetPoint( 3 ) ) );
        edges.push_back( Kratos::make_shared<EdgeType>( this->pGetPoint( 1 ), this->pGetPoint( 2 ), this->pGetPoint( 4 ) ) );
        edges.push_back( Kratos::make_shared<EdgeType>( this->pGetPoint( 2 ), this->pGetPoint( 0 ), this->pGetPoint( 5 ) ) );
        return edges;
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
    Matrix& ShapeFunctionsLocalGradients( Matrix& rResult,
            const CoordinatesArrayType& rPoint ) const override
    {
        rResult.resize( 6, 2,false );
        double thirdCoord = 1 - rPoint[0] - rPoint[1];
        double thirdCoord_DX = -1;
        double thirdCoord_DY = -1;

        noalias( rResult ) = ZeroMatrix( 6, 2 );
        rResult( 0, 0 ) = ( 4 * thirdCoord - 1 ) * thirdCoord_DX;
        rResult( 0, 1 ) = ( 4 * thirdCoord - 1 ) * thirdCoord_DY;
        rResult( 1, 0 ) =  4 * rPoint[0] - 1;
        rResult( 1, 1 ) =  0;
        rResult( 2, 0 ) =  0;
        rResult( 2, 1 ) =  4 * rPoint[1] - 1;
        rResult( 3, 0 ) =  4 * thirdCoord_DX * rPoint[0] + 4 * thirdCoord;
        rResult( 3, 1 ) =  4 * thirdCoord_DY * rPoint[0];
        rResult( 4, 0 ) =  4 * rPoint[1];
        rResult( 4, 1 ) =  4 * rPoint[0];
        rResult( 5, 0 ) =  4 * rPoint[1] * thirdCoord_DX;
        rResult( 5, 1 ) =  4 * rPoint[1] * thirdCoord_DY + 4 * thirdCoord;
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
    virtual Matrix& ShapeFunctionsGradients( Matrix& rResult, CoordinatesArrayType& rPoint )
    {
        rResult.resize( 6, 2 ,false);
        double thirdCoord = 1 - rPoint[0] - rPoint[1];
        double thirdCoord_DX = -1;
        double thirdCoord_DY = -1;

        noalias( rResult ) = ZeroMatrix( 6, 2 );
        rResult( 0, 0 ) = ( 4 * thirdCoord - 1 ) * thirdCoord_DX;
        rResult( 0, 1 ) = ( 4 * thirdCoord - 1 ) * thirdCoord_DY;
        rResult( 1, 0 ) =  4 * rPoint[0] - 1;
        rResult( 1, 1 ) =  0;
        rResult( 2, 0 ) =  0;
        rResult( 2, 1 ) =  4 * rPoint[1] - 1;
        rResult( 3, 0 ) =  4 * thirdCoord_DX * rPoint[0] + 4 * thirdCoord;
        rResult( 3, 1 ) =  4 * thirdCoord_DY * rPoint[0];
        rResult( 4, 0 ) =  4 * rPoint[1];
        rResult( 4, 1 ) =  4 * rPoint[0];
        rResult( 5, 0 ) =  4 * rPoint[1] * thirdCoord_DX;
        rResult( 5, 1 ) =  4 * rPoint[1] * thirdCoord_DY + 4 * thirdCoord;
        return rResult;
    }

    /**
     * returns the second order derivatives of all shape functions
     * in given arbitrary pointers
     * @param rResult a third order tensor which contains the second derivatives
     * @param rPoint the given point the second order derivatives are calculated in
     */
    ShapeFunctionsSecondDerivativesType& ShapeFunctionsSecondDerivatives( ShapeFunctionsSecondDerivativesType& rResult, const CoordinatesArrayType& rPoint ) const override
    {
        if ( rResult.size() != this->PointsNumber() )
        {
            // KLUDGE: While there is a bug in
            // ublas vector resize, I have to put this beside resizing!!
            ShapeFunctionsGradientsType temp( this->PointsNumber() );
            rResult.swap( temp );
        }

        rResult[0].resize( 2, 2 ,false);
        rResult[1].resize( 2, 2 ,false);
        rResult[2].resize( 2, 2 ,false);
        rResult[3].resize( 2, 2 ,false);
        rResult[4].resize( 2, 2 ,false);
        rResult[5].resize( 2, 2 ,false);

        rResult[0]( 0, 0 ) = 4.0;
        rResult[0]( 0, 1 ) = 4.0;
        rResult[0]( 1, 0 ) = 4.0;
        rResult[0]( 1, 1 ) = 4.0;
        rResult[1]( 0, 0 ) = 4.0;
        rResult[1]( 0, 1 ) = 0.0;
        rResult[1]( 1, 0 ) = 0.0;
        rResult[1]( 1, 1 ) = 0.0;
        rResult[2]( 0, 0 ) = 0.0;
        rResult[2]( 0, 1 ) = 0.0;
        rResult[2]( 1, 0 ) = 0.0;
        rResult[2]( 1, 1 ) = 4.0;
        rResult[3]( 0, 0 ) = -8.0;
        rResult[3]( 0, 1 ) = -4.0;
        rResult[3]( 1, 0 ) = -4.0;
        rResult[3]( 1, 1 ) = 0.0;
        rResult[4]( 0, 0 ) = 0.0;
        rResult[4]( 0, 1 ) = 4.0;
        rResult[4]( 1, 0 ) = 4.0;
        rResult[4]( 1, 1 ) = 0.0;
        rResult[5]( 0, 0 ) = 0.0;
        rResult[5]( 0, 1 ) = -4.0;
        rResult[5]( 1, 0 ) = -4.0;
        rResult[5]( 1, 1 ) = -8.0;

        return rResult;
    }

    /**
     * returns the third order derivatives of all shape functions
     * in given arbitrary pointers
     * @param rResult a fourth order tensor which contains the third derivatives
     * @param rPoint the given point the third order derivatives are calculated in
     */
    ShapeFunctionsThirdDerivativesType& ShapeFunctionsThirdDerivatives( ShapeFunctionsThirdDerivativesType& rResult, const CoordinatesArrayType& rPoint ) const override
    {
        if ( rResult.size() != this->PointsNumber() )
        {
            rResult.resize( this->PointsNumber() );
        }

        for ( IndexType i = 0; i < rResult.size(); i++ )
        {
            rResult[i].resize( this->PointsNumber() );
        }

        rResult[0][0].resize( 2, 2,false );

        rResult[0][1].resize( 2, 2,false );
        rResult[1][0].resize( 2, 2,false );
        rResult[1][1].resize( 2, 2 ,false);
        rResult[2][0].resize( 2, 2 ,false);
        rResult[2][1].resize( 2, 2 ,false);
        rResult[3][0].resize( 2, 2 ,false);
        rResult[3][1].resize( 2, 2 ,false);
        rResult[4][0].resize( 2, 2 ,false);
        rResult[4][1].resize( 2, 2 ,false);
        rResult[5][0].resize( 2, 2 ,false);
        rResult[5][1].resize( 2, 2 ,false);


        for ( int i = 0; i < 6; i++ )
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
        if (this->EdgesAreStraight()) {
            return GeometryUtils::PointLocalCoordinatesStraightEdgesTriangle(*this, rResult, rPoint);
        } else {
            return BaseType::PointLocalCoordinates( rResult, rPoint );
        }
    }

    ///@}
    ///@name Friends
    ///@{

    ///@}

protected:
    /**
     * There are no protected members in class Triangle3D6
     */

private:
    ///@name Static Member Variables
    ///@{

    static const GeometryData msGeometryData;

    static const GeometryDimension msGeometryDimension;

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

    Triangle3D6(): BaseType( PointsArrayType(), &msGeometryData ) {}

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
        IntegrationPointsArrayType integration_points = all_integration_points[static_cast<int>(ThisMethod)];
        //number of integration points
        const int integration_points_number = integration_points.size();
        //number of nodes in current geometry
        const int points_number = 6;
        //setting up return matrix
        Matrix shape_function_values( integration_points_number, points_number );
        //loop over all integration points

        for ( int pnt = 0; pnt < integration_points_number; pnt++ )
        {
            double thirdCoord = 1 - integration_points[pnt].X() - integration_points[pnt].Y();

            shape_function_values( pnt, 0 ) = thirdCoord * ( 2 * thirdCoord - 1 ) ;
            shape_function_values( pnt, 1 ) = integration_points[pnt].X() * ( 2 * integration_points[pnt].X() - 1 ) ;
            shape_function_values( pnt, 2 ) =  integration_points[pnt].Y() * ( 2 * integration_points[pnt].Y() - 1 ) ;
            shape_function_values( pnt, 3 ) =  4 * thirdCoord * integration_points[pnt].X();
            shape_function_values( pnt, 4 ) =  4 * integration_points[pnt].X() * integration_points[pnt].Y();
            shape_function_values( pnt, 5 ) =  4 * integration_points[pnt].Y() * thirdCoord;

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
        IntegrationPointsArrayType integration_points = all_integration_points[static_cast<int>(ThisMethod)];
        //number of integration points
        const int integration_points_number = integration_points.size();
        ShapeFunctionsGradientsType d_shape_f_values( integration_points_number );
        //initialising container
        //std::fill(d_shape_f_values.begin(), d_shape_f_values.end(), Matrix(4,2));
        //loop over all integration points

        for ( int pnt = 0; pnt < integration_points_number; pnt++ )
        {
            Matrix result( 6, 2 );
            double thirdCoord = 1 - integration_points[pnt].X() - integration_points[pnt].Y();
            double thirdCoord_DX = -1;
            double thirdCoord_DY = -1;

            noalias( result ) = ZeroMatrix( 6, 2 );
            result( 0, 0 ) = ( 4 * thirdCoord - 1 ) * thirdCoord_DX;
            result( 0, 1 ) = ( 4 * thirdCoord - 1 ) * thirdCoord_DY;
            result( 1, 0 ) =  4 * integration_points[pnt].X() - 1;
            result( 1, 1 ) =  0;
            result( 2, 0 ) =  0;
            result( 2, 1 ) =  4 * integration_points[pnt].Y() - 1;
            result( 3, 0 ) =  4 * thirdCoord_DX * integration_points[pnt].X() + 4 * thirdCoord;
            result( 3, 1 ) =  4 * thirdCoord_DY * integration_points[pnt].X();
            result( 4, 0 ) =  4 * integration_points[pnt].Y();
            result( 4, 1 ) =  4 * integration_points[pnt].X();
            result( 5, 0 ) =  4 * integration_points[pnt].Y() * thirdCoord_DX;
            result( 5, 1 ) =  4 * integration_points[pnt].Y() * thirdCoord_DY + 4 * thirdCoord;

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
                Quadrature<TriangleGaussLegendreIntegrationPoints3, 2, IntegrationPoint<3> >::GenerateIntegrationPoints()
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
                Triangle3D6<TPointType>::CalculateShapeFunctionsIntegrationPointsValues(
                    GeometryData::IntegrationMethod::GI_GAUSS_1 ),
                Triangle3D6<TPointType>::CalculateShapeFunctionsIntegrationPointsValues(
                    GeometryData::IntegrationMethod::GI_GAUSS_2 ),
                Triangle3D6<TPointType>::CalculateShapeFunctionsIntegrationPointsValues(
                    GeometryData::IntegrationMethod::GI_GAUSS_3 )
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
                Triangle3D6<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients( GeometryData::IntegrationMethod::GI_GAUSS_1 ),
                Triangle3D6<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients( GeometryData::IntegrationMethod::GI_GAUSS_2 ),
                Triangle3D6<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients( GeometryData::IntegrationMethod::GI_GAUSS_3 )
            }
        };
        return shape_functions_local_gradients;
    }

    /** This method gives all shape functions values
     * evaluated at the rCoordinates provided
     * @return Vector of values of shape functions \f$ F_{i} \f$
     * where i is the shape function index (for NURBS it is the index
     * of the local enumeration in the element).
     *
     * @see ShapeFunctionValue
     * @see ShapeFunctionsLocalGradients
     * @see ShapeFunctionLocalGradient
     */
    Vector& ShapeFunctionsValues (
        Vector &rResult,
        const CoordinatesArrayType& rCoordinates) const override
    {
        if(rResult.size() != 6) {
            rResult.resize(6, false);
        }

        const double xi = rCoordinates[0];
        const double eta = rCoordinates[1];
        const double zeta = 1.0 - xi - eta;

        rResult[0] = zeta * (2.0 * zeta - 1.0);
        rResult[1] = xi * (2.0 * xi - 1.0);
        rResult[2] = eta * (2.0 * eta - 1.0);
        rResult[3] = 4.0 * xi * zeta;
        rResult[4] = 4.0 * xi * eta;
        rResult[5] = 4.0 * eta * zeta;

        return rResult;
    }


    /**
     * @brief Checks if edges are straight. We iterate though all edges
     * and check that the sum of 0-2 and 2-1 segments is no bigger than 0-1.
     * @return bool edges are straight or not
     */
    bool EdgesAreStraight() const
    {
        constexpr double tol = 1e-6;
        constexpr std::array<std::array<size_t, 3>, 6> edges{
            {{0, 1, 3}, {1, 2, 4}, {2, 0, 5}}};
        const auto& r_points = this->Points();
        for (const auto& r_edge : edges) {
            const double a = MathUtils<double>::Norm3(r_points[r_edge[0]] - r_points[r_edge[1]]);
            const double b = MathUtils<double>::Norm3(r_points[r_edge[1]] - r_points[r_edge[2]]);
            const double c = MathUtils<double>::Norm3(r_points[r_edge[2]] - r_points[r_edge[0]]);
            if (b + c > a*(1.0+tol) ) {
                return false;
            }
        }
        return true;
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

    template<class TOtherPointType> friend class Triangle3D6;

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
    Triangle3D6<TPointType>& rThis );
/**
 * output stream functions
 */
template<class TPointType> inline std::ostream& operator << (
    std::ostream& rOStream,
    const Triangle3D6<TPointType>& rThis )
{
    rThis.PrintInfo( rOStream );
    rOStream << std::endl;
    rThis.PrintData( rOStream );
    return rOStream;
}

///@}

template<class TPointType> const
GeometryData Triangle3D6<TPointType>::msGeometryData(
    &msGeometryDimension,
    GeometryData::IntegrationMethod::GI_GAUSS_2,
    Triangle3D6<TPointType>::AllIntegrationPoints(),
    Triangle3D6<TPointType>::AllShapeFunctionsValues(),
    AllShapeFunctionsLocalGradients()
);

template<class TPointType> const
GeometryDimension Triangle3D6<TPointType>::msGeometryDimension(3, 2);

}// namespace Kratos.
