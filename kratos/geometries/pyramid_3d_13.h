//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//	                 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher (https://github.com/philbucher)
//                   Ashish Darekar

#pragma once

// System includes
#include <cmath> // std::abs for double

// External includes

// Project includes
#include "includes/define.h"
#include "geometries/geometry.h"
#include "utilities/integration_utilities.h"
#include "integration/pyramid_gauss_legendre_integration_points.h"
#include "integration/pyramid_gauss_lobatto_integration_points.h"

namespace Kratos {

///@name Kratos Classes
///@{

/**
 * @class Pyramid3D13
 * @ingroup KratosCore
 * @brief A 13 node pyramid geometry with quadratic shape functions
 * @details The node ordering corresponds with:
 *                     4
 *                   ,/|\
 *                 ,/ .'|\
 *               ,/   μ | \
 *             ,/    .^ | `.
 *           ,9      || 12  \
 *         ,/       .'|  |   \
 *       ,/        10 +--|----\-----> η
 *      0-------8--.'-`\ -3    11
 *       `\        |    `\ `\    \
 *         `5     .'      `\ `7   \
 *           `\   |         ξ  `\  \
 *             `\.'              `\`\
 *                1--------6---------2
 *
 *
 * @author Philipp Bucher, Ashish Darekar
 */
template<class TPointType>
class Pyramid3D13 : public Geometry<TPointType>
{
public:
    ///@}
    ///@name Type Definitions
    ///@{

    /// Geometry as base class.
    typedef Geometry<TPointType> BaseType;

    /// Pointer definition of Pyramid3D13
    KRATOS_CLASS_POINTER_DEFINITION(Pyramid3D13);

    /** Integration methods implemented in geometry.
     */
    typedef GeometryData::IntegrationMethod IntegrationMethod;

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
    typedef typename BaseType::PointsArrayType PointsArrayType;

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

    /**
     * A third order tensor to hold shape functions' local
     * gradients. ShapefunctionsLocalGradients function return this
     * type as its result.
     */
    typedef typename BaseType::ShapeFunctionsGradientsType ShapeFunctionsGradientsType;

    /**
     * A fourth order tensor used as shape functions' local
     * gradients container in geometry.
     */
    typedef typename BaseType::ShapeFunctionsLocalGradientsContainerType
    ShapeFunctionsLocalGradientsContainerType;

    /**
    * Type of coordinates array
    */
    typedef typename BaseType::CoordinatesArrayType CoordinatesArrayType;

    /**
     * A Vector of counted pointers to Geometries. Used for
     * returning edges of the geometry.
     */
    typedef typename BaseType::GeometriesArrayType GeometriesArrayType;


    ///@}
    ///@name Life Cycle
    ///@{

    explicit Pyramid3D13(
        typename PointType::Pointer pPoint1,
        typename PointType::Pointer pPoint2,
        typename PointType::Pointer pPoint3,
        typename PointType::Pointer pPoint4,
        typename PointType::Pointer pPoint5,
        typename PointType::Pointer pPoint6,
        typename PointType::Pointer pPoint7,
        typename PointType::Pointer pPoint8,
        typename PointType::Pointer pPoint9,
        typename PointType::Pointer pPoint10,
        typename PointType::Pointer pPoint11,
        typename PointType::Pointer pPoint12,
        typename PointType::Pointer pPoint13)
        : BaseType( PointsArrayType(), &msGeometryData )
    {
        this->Points().reserve(13);
        this->Points().push_back(pPoint1);
        this->Points().push_back(pPoint2);
        this->Points().push_back(pPoint3);
        this->Points().push_back(pPoint4);
        this->Points().push_back(pPoint5);
        this->Points().push_back(pPoint6);
        this->Points().push_back(pPoint7);
        this->Points().push_back(pPoint8);
        this->Points().push_back(pPoint9);
        this->Points().push_back(pPoint10);
        this->Points().push_back(pPoint11);
        this->Points().push_back(pPoint12);
        this->Points().push_back(pPoint13);
    }

    explicit Pyramid3D13( const PointsArrayType& ThisPoints )
        : BaseType( ThisPoints, &msGeometryData )
    {
        KRATOS_ERROR_IF( this->PointsNumber() != 13 ) << "Invalid points number. Expected 13, given " << this->PointsNumber() << std::endl;
    }

    /// Constructor with Geometry Id
    explicit Pyramid3D13(
        const IndexType GeometryId,
        const PointsArrayType& rThisPoints
    ) : BaseType( GeometryId, rThisPoints, &msGeometryData)
    {
        KRATOS_ERROR_IF( this->PointsNumber() != 13 ) << "Invalid points number. Expected 13, given " << this->PointsNumber() << std::endl;
    }

    /// Constructor with Geometry Name
    explicit Pyramid3D13(
        const std::string& rGeometryName,
        const PointsArrayType& rThisPoints
    ) : BaseType(rGeometryName, rThisPoints, &msGeometryData)
    {
        KRATOS_ERROR_IF(this->PointsNumber() != 13) << "Invalid points number. Expected 13, given " << this->PointsNumber() << std::endl;
    }

    /** Copy constructor.
     Construct this geometry as a copy of given geometry.

    @note This copy constructor don't copy the points and new
    geometry shares points with given source geometry. It's
    obvious that any change to this new geometry's point affect
    source geometry's points too.
    */
    Pyramid3D13(Pyramid3D13 const& rOther)
    : BaseType(rOther)
    {
    }

    /** Copy constructor from a geometry with other point type.
     Construct this geometry as a copy of given geometry which
    has different type of points. The given goemetry's
    TOtherPointType* must be implicitly convertible to this
    geometry PointType.

    @note This copy constructor don't copy the points and new
    geometry shares points with given source geometry. It's
    obvious that any change to this new geometry's point affect
    source geometry's points too.
    */
    template<class TOtherPointType> Pyramid3D13(Pyramid3D13<TOtherPointType> const& rOther)
    : BaseType(rOther)
    {
    }

    /**
     * @brief Gets the geometry family.
     * @details This function returns the family type of the geometry. The geometry family categorizes the geometry into a broader classification, aiding in its identification and processing.
     * @return GeometryData::KratosGeometryFamily The geometry family.
     */
    GeometryData::KratosGeometryFamily GetGeometryFamily() const override
    {
        return GeometryData::KratosGeometryFamily::Kratos_Pyramid;
    }

    /**
     * @brief Gets the geometry type.
     * @details This function returns the specific type of the geometry. The geometry type provides a more detailed classification of the geometry.
     * @return GeometryData::KratosGeometryType The specific geometry type.
     */
    GeometryData::KratosGeometryType GetGeometryType() const override
    {
        return GeometryData::KratosGeometryType::Kratos_Pyramid3D13;
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

    /** Assignment operator.

    @note This operator don't copy the points and this
    geometry shares points with given source geometry. It's
    obvious that any change to this geometry's point affect
    source geometry's points too.

    @see Clone
    @see ClonePoints
    */
    Pyramid3D13& operator=(const Pyramid3D13& rOther)
    {
        BaseType::operator=(rOther);

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
    Pyramid3D13& operator=(Pyramid3D13<TOtherPointType> const & rOther)
    {
        BaseType::operator=(rOther);

        return *this;
    }

    ///@}
    ///@name Create Methods
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
        return typename BaseType::Pointer( new Pyramid3D13( NewGeometryId, rThisPoints ) );
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
        auto p_geometry = typename BaseType::Pointer( new Pyramid3D13( NewGeometryId, rGeometry.Points() ) );
        p_geometry->SetData(rGeometry.GetData());
        return p_geometry;
    }

    ///@}
    ///@name Informations
    ///@{

    /**
     * @brief This method gives you number of all edges of this geometry.
     * @return SizeType contains number of this geometry edges.
     * @see EdgesNumber()
     * @see Edges()
     * @see GenerateEdges()
     * @see FacesNumber()
     * @see Faces()
     * @see GenerateFaces()
     */
    SizeType EdgesNumber() const override
    {
        return 8;
    }

    /**
     * @brief Returns the number of faces of the current geometry.
     * @see EdgesNumber
     * @see Edges
     * @see Faces
     */
    SizeType FacesNumber() const override
    {
        return 5;
    }

    /**
     * @brief This method calculates and returns volume of this geometry.
     * @details For one and two dimensional geometry it returns zero and for three dimensional it gives volume of geometry.
     * @return double value contains volume.
     * @see Length()
     * @see Area()
     * @see DomainSize()
     */
    double Volume() const override
    {
        return IntegrationUtilities::ComputeVolume3DGeometry(*this);
    }

    /**
     * This method calculates and returns length, area or volume of
     * this geometry depending on its dimension. For one dimensional
     * geometry it returns its length, for two dimensional it gives area
     * and for three dimensional geometries it gives its volume.
     *
     * @return double value contains length, area or volume.
     * @see Length()
     * @see Area()
     * @see Volume()
     */
    double DomainSize() const override
    {
        return Volume();
    }

    /**
    * @brief Checks if given point in local space coordinates of this geometry
    *        is inside the geometry boundaries.
    * @param rPointLocalCoordinates the point on the geometry,
    *        which shall be checked if it lays within
    *        the boundaries.
    * @param Tolerance the tolerance to the boundary.
    * @return  0 -> failed
    *          1 -> inside /boundary/ vertex
    */
    int IsInsideLocalSpace(
        const CoordinatesArrayType& rPointLocalCoordinates,
        const double Tolerance = std::numeric_limits<double>::epsilon()
    ) const override
    {
        if ( std::abs( rPointLocalCoordinates[0] ) <= (1.0 + Tolerance) ) {
            if ( std::abs( rPointLocalCoordinates[1] ) <= (1.0 + Tolerance) ) {
                if ( std::abs( rPointLocalCoordinates[2] ) <= (1.0 + Tolerance) ) {
                    if ( (std::abs(rPointLocalCoordinates[0]) +
                          std::abs(rPointLocalCoordinates[1]) +
                          rPointLocalCoordinates[2]) <= (1.0 + Tolerance) ) {
                        return 1;
                    }
                }
            }
        }

        return 0;
    }

    /**
    * Returns a matrix of the local coordinates of all points
    * @param rResult a Matrix that will be overwritten by the results
    * @return the coordinates of all points of the current geometry
    */
    Matrix& PointsLocalCoordinates( Matrix& rResult ) const override
    {
        if ( rResult.size1() != 13 || rResult.size2() != 3 )
            rResult.resize( 13, 3, false );

        rResult( 0, 0 ) = -1.0;
        rResult( 0, 1 ) = -1.0;
        rResult( 0, 2 ) = -1.0;

        rResult( 1, 0 ) = +1.0;
        rResult( 1, 1 ) = -1.0;
        rResult( 1, 2 ) = -1.0;

        rResult( 2, 0 ) = +1.0;
        rResult( 2, 1 ) = +1.0;
        rResult( 2, 2 ) = -1.0;

        rResult( 3, 0 ) = -1.0;
        rResult( 3, 1 ) = +1.0;
        rResult( 3, 2 ) = -1.0;

        rResult( 4, 0 ) =  0.0;
        rResult( 4, 1 ) =  0.0;
        rResult( 4, 2 ) = +1.0;

        rResult( 5, 0 ) =  0.0;
        rResult( 5, 1 ) = -0.5;
        rResult( 5, 2 ) = -1.0;

        rResult( 6, 0 ) = +0.5;
        rResult( 6, 1 ) =  0.0;
        rResult( 6, 2 ) = -1.0;

        rResult( 7, 0 ) =  0.0;
        rResult( 7, 1 ) = +0.5;
        rResult( 7, 2 ) = -1.0;

        rResult( 8, 0 ) = +0.5;
        rResult( 8, 1 ) =  0.0;
        rResult( 8, 2 ) = -1.0;

        rResult( 9, 0 ) = -0.5;
        rResult( 9, 1 ) = -0.5;
        rResult( 9, 2 ) =  0.0;

        rResult( 10, 0 ) = +0.5;
        rResult( 10, 1 ) = -0.5;
        rResult( 10, 2 ) = 0.0;

        rResult( 11, 0 ) = +0.5;
        rResult( 11, 1 ) = +0.5;
        rResult( 11, 2 ) = 0.0;

        rResult( 12, 0 ) = -0.5;
        rResult( 12, 1 ) = +0.5;
        rResult( 12, 2 ) = 0.0;

        return rResult;
    }

    /**
     * Shape Function
     */

    /** This method gives all non-zero shape functions values
    evaluated at the rCoordinates provided

    @return Vector of values of shape functions \f$ F_{i} \f$
    where i is the shape function index (for NURBS it is the index
    of the local enumeration in the element).

    @see ShapeFunctionValue
    @see ShapeFunctionsLocalGradients
    @see ShapeFunctionLocalGradient
    */
    Vector& ShapeFunctionsValues(Vector &rResult, const CoordinatesArrayType& rCoordinates) const override
    {
        if(rResult.size() != 13) rResult.resize(13,false);

        for (std::size_t i=0; i<13; ++i) {
            rResult[i] = ShapeFunctionValue(i, rCoordinates);
        }

        return rResult;
    }

    /**
     * Calculates the value of a given shape function at a given point.
     *
     * @param ShapeFunctionIndex The number of the desired shape function
     * @param rPoint the given point in local coordinates at which the
     * value of the shape function is calculated
     *
     * @return the value of the shape function at the given point
     */
    double ShapeFunctionValue( IndexType ShapeFunctionIndex,
                                       const CoordinatesArrayType& rPoint ) const override
    {
        return ShapeFunctionValueImpl(ShapeFunctionIndex, rPoint);
    }

    /**
     * Calculates the values of all shape function in all integration points.
     * Integration points are expected to be given in local coordinates
     *
     * @param ThisMethod the current integration method
     * @return the matrix of values of every shape function in each integration point
     *
     */

    static Matrix CalculateShapeFunctionsIntegrationPointsValues(typename BaseType::IntegrationMethod ThisMethod)
    {
        IntegrationPointsContainerType all_integration_points = AllIntegrationPoints();
        IntegrationPointsArrayType integration_points = all_integration_points[static_cast<int>(ThisMethod)];
        //number of integration points
        const std::size_t integration_points_number = integration_points.size();
        //number of nodes in current geometry
        const std::size_t points_number = 13;
        //setting up return matrix
        Matrix shape_function_values( integration_points_number, points_number );

        //loop over all integration points
        for (std::size_t pnt = 0; pnt<integration_points_number; ++pnt) {
            for (std::size_t i=0; i<points_number; ++i) {
                shape_function_values( pnt, i ) = ShapeFunctionValueImpl(i, integration_points[pnt]);
            }
        }

        return shape_function_values;
    }

        /**
     * Calculates the gradients in terms of local coordinateds
     * of all shape functions in a given point.
     * @param rPoint the current point at which the gradients are calculated
     * @return the gradients of all shape functions
     * \f$ \frac{\partial N^i}{\partial \xi_j} \f$
     */
    Matrix& ShapeFunctionsLocalGradients(
        Matrix& rResult,
        const CoordinatesArrayType& rPoint
        ) const override
    {
        if(rResult.size1() != this->PointsNumber() || rResult.size2() != this->LocalSpaceDimension())
            rResult.resize(this->PointsNumber(),this->LocalSpaceDimension(),false);

        CalculateShapeFunctionsLocalGradients(rResult, rPoint);

        return rResult;
    }

    /**
     * Calculates the gradients in terms of local coordinateds
     * of all shape functions in a given point.
     *
     * @param rPoint the current point at which the gradients are calculated
     * @return the gradients of all shape functions
     * \f$ \frac{\partial N^i}{\partial \xi_j} \f$
     */
    static Matrix& CalculateShapeFunctionsLocalGradients(
        Matrix& rResult,
        const CoordinatesArrayType& rPoint
        )
    {
        rResult.resize( 13, 3, false );
        noalias( rResult ) = ZeroMatrix( 13, 3 );

        rResult( 0, 0 ) = (+0.0625) * (1 - rPoint[1]) * (1 - rPoint[2]) * (1 + 6*rPoint[0] + rPoint[1] + 4*rPoint[0]*rPoint[1] + rPoint[2] + 2*rPoint[0]*rPoint[2] - rPoint[1]*rPoint[2] + 4*rPoint[0]*rPoint[1]*rPoint[2]) ;
        rResult( 0, 1 ) = (+0.0625) * (1 - rPoint[0]) * (1 - rPoint[2]) * (1 + rPoint[0] + 6*rPoint[1] + 4*rPoint[0]*rPoint[1] + rPoint[2] - rPoint[0]*rPoint[2] + 2*rPoint[1]*rPoint[2] + 4*rPoint[0]*rPoint[1]*rPoint[2]) ;
        rResult( 0, 2 ) = (+0.125) * (1 - rPoint[0]) * (1 - rPoint[1]) * (1 + rPoint[0] + rPoint[1] + 2*rPoint[2] + rPoint[0]*rPoint[2] + rPoint[1]*rPoint[2] + 2*rPoint[0]*rPoint[1]*rPoint[2]) ;

        rResult( 1, 0 ) = (-0.0625) * (1 - rPoint[1]) * (1 - rPoint[2]) * (1 - 6*rPoint[0] + rPoint[1] - 4*rPoint[0]*rPoint[1] + rPoint[2] - 2*rPoint[0]*rPoint[2] - rPoint[1]*rPoint[2] - 4*rPoint[0]*rPoint[1]*rPoint[2]) ;
        rResult( 1, 1 ) = (+0.0625) * (1 + rPoint[0]) * (1 - rPoint[2]) * (1 - rPoint[0] + 6*rPoint[1] - 4*rPoint[0]*rPoint[1] + rPoint[2] + rPoint[0]*rPoint[2] + 2*rPoint[1]*rPoint[2] - 4*rPoint[0]*rPoint[1]*rPoint[2]) ;
        rResult( 1, 2 ) = (+0.125) * (1 + rPoint[0]) * (1 - rPoint[1]) * (1 - rPoint[0] + rPoint[1] + 2*rPoint[2] - rPoint[0]*rPoint[2] + rPoint[1]*rPoint[2] - 2*rPoint[0]*rPoint[1]*rPoint[2]) ;

        rResult( 2, 0 ) = (-0.0625) * (1 + rPoint[1]) * (1 - rPoint[2]) * (1 - 6*rPoint[0] - rPoint[1] + 4*rPoint[0]*rPoint[1] + rPoint[2] - 2*rPoint[0]*rPoint[2] + rPoint[1]*rPoint[2] + 4*rPoint[0]*rPoint[1]*rPoint[2]) ;
        rResult( 2, 1 ) = (-0.0625) * (1 + rPoint[0]) * (1 - rPoint[2]) * (1 - rPoint[0] - 6*rPoint[1] + 4*rPoint[0]*rPoint[1] + rPoint[2] + rPoint[0]*rPoint[2] - 2*rPoint[1]*rPoint[2] + 4*rPoint[0]*rPoint[1]*rPoint[2]) ;
        rResult( 2, 2 ) = (+0.125) * (1 + rPoint[0]) * (1 + rPoint[1]) * (1 - rPoint[0] - rPoint[1] + 2*rPoint[2] - rPoint[0]*rPoint[2] - rPoint[1]*rPoint[2] + 2*rPoint[0]*rPoint[1]*rPoint[2]) ;

        rResult( 3, 0 ) = (+0.0625) * (1 + rPoint[1]) * (1 - rPoint[2]) * (1 + 6*rPoint[0] - rPoint[1] - 4*rPoint[0]*rPoint[1] + rPoint[2] + 2*rPoint[0]*rPoint[2] + rPoint[1]*rPoint[2] - 4*rPoint[0]*rPoint[1]*rPoint[2]) ;
        rResult( 3, 1 ) = (-0.0625) * (1 - rPoint[0]) * (1 - rPoint[2]) * (1 + rPoint[0] - 6*rPoint[1] - 4*rPoint[0]*rPoint[1] + rPoint[2] - rPoint[0]*rPoint[2] - 2*rPoint[1]*rPoint[2] - 4*rPoint[0]*rPoint[1]*rPoint[2]) ;
        rResult( 3, 2 ) = (+0.125) * (1 - rPoint[0]) * (1 + rPoint[1]) * (1 + rPoint[0] - rPoint[1] + 2*rPoint[2] + rPoint[0]*rPoint[2] - rPoint[1]*rPoint[2] - 2*rPoint[0]*rPoint[1]*rPoint[2]) ;

        rResult( 4, 0 ) = 0.00 ;
        rResult( 4, 1 ) = 0.00 ;
        rResult( 4, 2 ) = (0.5) + rPoint[2] ;

        rResult( 5, 0 ) = (-0.25) * (rPoint[0]) * (1 - rPoint[1]) * (1 - rPoint[2]) * (2 + rPoint[1] + rPoint[1]*rPoint[2]) ;
        rResult( 5, 1 ) = (-0.125) * (1 - std::pow(rPoint[0],2.0)) * (1 - rPoint[2]) * (1 + 2*rPoint[1] - rPoint[2] + 2*rPoint[1]*rPoint[2]) ;
        rResult( 5, 2 ) = (-0.25) * (1 - std::pow(rPoint[0],2.0)) * (1 - rPoint[1]) * (1 + rPoint[1]*rPoint[2]) ;

        rResult( 6, 0 ) = (0.125) * (1 - std::pow(rPoint[1],2.0)) * (1 - rPoint[2]) * (1 - 2*rPoint[0] - rPoint[2] - 2*rPoint[0]*rPoint[2]) ;
        rResult( 6, 1 ) = (-0.25) * (1 + rPoint[0]) * (rPoint[1]) *  (1 - rPoint[2]) * (2 - rPoint[0] - rPoint[0]*rPoint[2]) ;
        rResult( 6, 2 ) = (-0.25) * (1 + rPoint[0]) * (1 - std::pow(rPoint[1],2.0)) *  (1 - rPoint[0]*rPoint[2]) ;

        rResult( 7, 0 ) = (-0.25) * (rPoint[0]) * (1 + rPoint[1]) * (1 - rPoint[2]) * (2 - rPoint[1] - rPoint[1]*rPoint[2]) ;
        rResult( 7, 1 ) = (+0.125) * (1 - std::pow(rPoint[0],2.0)) * (1 - rPoint[2]) * (1 - 2*rPoint[1] - rPoint[2] - 2*rPoint[1]*rPoint[2]) ;
        rResult( 7, 2 ) = (-0.25) * (1 - std::pow(rPoint[0],2.0)) * (1 + rPoint[1]) * (1 - rPoint[1]*rPoint[2]) ;

        rResult( 8, 0 ) = (-0.125) * (1 - std::pow(rPoint[1],2.0)) * (1 - rPoint[2]) * (1 + 2*rPoint[0] - rPoint[2] + 2*rPoint[0]*rPoint[2]) ;
        rResult( 8, 1 ) = (-0.25) * (1 - rPoint[0]) * (rPoint[1]) *  (1 - rPoint[2]) * (2 + rPoint[0] + rPoint[0]*rPoint[2]) ;
        rResult( 8, 2 ) = (-0.25) * (1 - rPoint[0]) * (1 - std::pow(rPoint[1],2.0)) *  (1 + rPoint[0]*rPoint[2]) ;

        rResult( 9, 0 ) = (-0.25) * (1 - rPoint[1]) * (1 - std::pow(rPoint[2],2.0)) ;
        rResult( 9, 1 ) = (-0.25) * (1 - rPoint[0]) * (1 - std::pow(rPoint[2],2.0)) ;
        rResult( 9, 2 ) = (-0.5) * (1 - rPoint[0]) * (1 - rPoint[1]) * (rPoint[2]) ;

        rResult( 10, 0 ) = (+0.25) * (1 - rPoint[1]) * (1 - std::pow(rPoint[2],2.0)) ;
        rResult( 10, 1 ) = (-0.25) * (1 + rPoint[0]) * (1 - std::pow(rPoint[2],2.0)) ;
        rResult( 10, 2 ) = (-0.5) * (1 + rPoint[0]) * (1 - rPoint[1]) * (rPoint[2]) ;

        rResult( 11, 0 ) = (+0.25) * (1 + rPoint[1]) * (1 - std::pow(rPoint[2],2.0)) ;
        rResult( 11, 1 ) = (+0.25) * (1 + rPoint[0]) * (1 - std::pow(rPoint[2],2.0)) ;
        rResult( 11, 2 ) = (-0.5) * (1 + rPoint[0]) * (1 + rPoint[1]) * (rPoint[2]) ;

        rResult( 12, 0 ) = (-0.25) * (1 + rPoint[1]) * (1 - std::pow(rPoint[2],2.0)) ;
        rResult( 12, 1 ) = (+0.25) * (1 - rPoint[0]) * (1 - std::pow(rPoint[2],2.0)) ;
        rResult( 12, 2 ) = (-0.5) * (1 - rPoint[0]) * (1 + rPoint[1]) * (rPoint[2]) ;

        return rResult;
    }

    /**
     * Calculates the local gradients of all shape functions in all integration points.
     * Integration points are expected to be given in local coordinates
     *
     * @param ThisMethod the current integration method
     * @return the vector of the gradients of all shape functions
     * in each integration point
     *
     */
    static ShapeFunctionsGradientsType CalculateShapeFunctionsIntegrationPointsLocalGradients(typename BaseType::IntegrationMethod ThisMethod)
    {
        IntegrationPointsContainerType all_integration_points = AllIntegrationPoints();
        IntegrationPointsArrayType integration_points = all_integration_points[static_cast<int>(ThisMethod)]; //number of integration points

        const std::size_t integration_points_number = integration_points.size();
        ShapeFunctionsGradientsType d_shape_f_values(integration_points_number); //initialising container

        Matrix result;

        //loop over all integration points
        for (std::size_t pnt = 0; pnt<integration_points_number; ++pnt) {
            d_shape_f_values[pnt] = CalculateShapeFunctionsLocalGradients(result, integration_points[pnt]);
        }

        return d_shape_f_values;
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
        return "3 dimensional pyramid with 13 nodes in 3D space";
    }

    /** Print information about this object.

    @param rOStream Stream to print into it.
    @see PrintData()
    @see Info()
    */
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
    }

    /** Print geometry's data into given stream. Prints it's points
     by the order they stored in the geometry and then center
    point of geometry.

    @param rOStream Stream to print into it.
    @see PrintInfo()
    @see Info()
    */
    void PrintData(std::ostream& rOStream) const override
    {
        BaseType::PrintData(rOStream);
    }

    ///@}

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

    Pyramid3D13() : BaseType( PointsArrayType(), &msGeometryData ) {}

    ///@}
    ///@name Private Operations
    ///@{

    // see ShapeFunctionValue
    // special function such that it can also be used in the static methods
    static double ShapeFunctionValueImpl(
        IndexType ShapeFunctionIndex,
        const CoordinatesArrayType& rPoint)
    {
        switch ( ShapeFunctionIndex )
        {
        case 0:
            return( (-0.0625) * (1 - rPoint[0]) * (1 - rPoint[1]) * (1 - rPoint[2]) * (4 + 3*rPoint[0] + 3*rPoint[1] + 2*rPoint[0]*rPoint[1] + 2*rPoint[2] + rPoint[0]*rPoint[2] + rPoint[1]*rPoint[2] + 2*rPoint[0]*rPoint[1]*rPoint[2]) );
        case 1:
            return( (-0.0625) * (1 + rPoint[0]) * (1 - rPoint[1]) * (1 - rPoint[2]) * (4 - 3*rPoint[0] + 3*rPoint[1] - 2*rPoint[0]*rPoint[1] + 2*rPoint[2] - rPoint[0]*rPoint[2] + rPoint[1]*rPoint[2] - 2*rPoint[0]*rPoint[1]*rPoint[2]) );
        case 2:
            return( (-0.0625) * (1 + rPoint[0]) * (1 + rPoint[1]) * (1 - rPoint[2]) * (4 - 3*rPoint[0] - 3*rPoint[1] + 2*rPoint[0]*rPoint[1] + 2*rPoint[2] - rPoint[0]*rPoint[2] - rPoint[1]*rPoint[2] + 2*rPoint[0]*rPoint[1]*rPoint[2]) );
        case 3:
            return( (-0.0625) * (1 - rPoint[0]) * (1 + rPoint[1]) * (1 - rPoint[2]) * (4 + 3*rPoint[0] - 3*rPoint[1] - 2*rPoint[0]*rPoint[1] + 2*rPoint[2] + rPoint[0]*rPoint[2] - rPoint[1]*rPoint[2] - 2*rPoint[0]*rPoint[1]*rPoint[2]) );
        case 4:
            return( (0.5) * (rPoint[2]) * (1 + rPoint[2]) );
        case 5:
            return( (0.125) * (1 - std::pow(rPoint[0],2.0)) * (1 - rPoint[1]) * (1 - rPoint[2]) * (2 + rPoint[1] + rPoint[1]*rPoint[2]) );
        case 6:
            return( (0.125) * (1 + rPoint[0]) * (1 - std::pow(rPoint[1],2.0)) * (1 - rPoint[2]) * (2 - rPoint[0] - rPoint[0]*rPoint[2]) );
        case 7:
            return( (0.125) * (1 - std::pow(rPoint[0],2.0)) * (1 + rPoint[1]) * (1 - rPoint[2]) * (2 - rPoint[1] - rPoint[1]*rPoint[2]) );
        case 8:
            return( (0.125) * (1 - rPoint[0]) * (1 - std::pow(rPoint[1],2.0)) * (1 - rPoint[2]) * (2 + rPoint[0] + rPoint[0]*rPoint[2]) );
        case 9:
            return( (0.25) * (1 - rPoint[0]) * (1 - rPoint[1]) * (1 - std::pow(rPoint[2],2.0)) );
        case 10:
            return( (0.25) * (1 + rPoint[0]) * (1 - rPoint[1]) * (1 - std::pow(rPoint[2],2.0)) );
        case 11:
            return( (0.25) * (1 + rPoint[0]) * (1 + rPoint[1]) * (1 - std::pow(rPoint[2],2.0)) );
        case 12:
            return( (0.25) * (1 - rPoint[0]) * (1 + rPoint[1]) * (1 - std::pow(rPoint[2],2.0)) );
        default:
            KRATOS_ERROR << "Wrong index of shape function:" << ShapeFunctionIndex  << std::endl;
        }

        return 0;
    }

    static const IntegrationPointsContainerType AllIntegrationPoints()
    {
        IntegrationPointsContainerType integration_points =
        {
            {
                Quadrature < PyramidGaussLegendreIntegrationPoints1,
                3, IntegrationPoint<3> >::GenerateIntegrationPoints(),
                Quadrature < PyramidGaussLegendreIntegrationPoints2,
                3, IntegrationPoint<3> >::GenerateIntegrationPoints(),
                Quadrature < PyramidGaussLegendreIntegrationPoints3,
                3, IntegrationPoint<3> >::GenerateIntegrationPoints(),
                Quadrature < PyramidGaussLegendreIntegrationPoints4,
                3, IntegrationPoint<3> >::GenerateIntegrationPoints(),
                Quadrature < PyramidGaussLegendreIntegrationPoints5,
                3, IntegrationPoint<3> >::GenerateIntegrationPoints(),
                Quadrature < PyramidGaussLobattoIntegrationPoints1,
                3, IntegrationPoint<3> >::GenerateIntegrationPoints()
            }
        };
        return integration_points;
    }

    static const ShapeFunctionsValuesContainerType AllShapeFunctionsValues()
    {
        ShapeFunctionsValuesContainerType shape_functions_values =
        {
            {
                Pyramid3D13<TPointType>::CalculateShapeFunctionsIntegrationPointsValues(GeometryData::IntegrationMethod::GI_GAUSS_1),
                Pyramid3D13<TPointType>::CalculateShapeFunctionsIntegrationPointsValues(GeometryData::IntegrationMethod::GI_GAUSS_2),
                Pyramid3D13<TPointType>::CalculateShapeFunctionsIntegrationPointsValues(GeometryData::IntegrationMethod::GI_GAUSS_3),
                Pyramid3D13<TPointType>::CalculateShapeFunctionsIntegrationPointsValues(GeometryData::IntegrationMethod::GI_GAUSS_4),
                Pyramid3D13<TPointType>::CalculateShapeFunctionsIntegrationPointsValues(GeometryData::IntegrationMethod::GI_GAUSS_5),
                Pyramid3D13<TPointType>::CalculateShapeFunctionsIntegrationPointsValues(GeometryData::IntegrationMethod::GI_LOBATTO_1)
            }
        };
        return shape_functions_values;
    }

    static const ShapeFunctionsLocalGradientsContainerType AllShapeFunctionsLocalGradients()
    {
        ShapeFunctionsLocalGradientsContainerType shape_functions_local_gradients =
        {
            {
                Pyramid3D13<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients(GeometryData::IntegrationMethod::GI_GAUSS_1),
                Pyramid3D13<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients(GeometryData::IntegrationMethod::GI_GAUSS_2),
                Pyramid3D13<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients(GeometryData::IntegrationMethod::GI_GAUSS_3),
                Pyramid3D13<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients(GeometryData::IntegrationMethod::GI_GAUSS_4),
                Pyramid3D13<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients(GeometryData::IntegrationMethod::GI_GAUSS_5),
                Pyramid3D13<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients(GeometryData::IntegrationMethod::GI_LOBATTO_1)
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

    template<class TOtherPointType> friend class Pyramid3D13;

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
inline std::istream& operator >> (std::istream& rIStream,
                    Pyramid3D13<TPointType>& rThis);

/// output stream function
template<class TPointType>
inline std::ostream& operator << (std::ostream& rOStream,
                    const Pyramid3D13<TPointType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

template<class TPointType> const
GeometryData Pyramid3D13<TPointType>::msGeometryData(
    &msGeometryDimension,
    GeometryData::IntegrationMethod::GI_GAUSS_2,
    Pyramid3D13<TPointType>::AllIntegrationPoints(),
    Pyramid3D13<TPointType>::AllShapeFunctionsValues(),
    AllShapeFunctionsLocalGradients()
);

template<class TPointType> const
GeometryDimension Pyramid3D13<TPointType>::msGeometryDimension(3, 3);

}  // namespace Kratos.
