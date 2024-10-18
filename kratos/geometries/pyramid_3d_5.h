//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher (https://github.com/philbucher)
//                   Ashish Darekar
//

#pragma once

// System includes
#include <cmath> // std::abs for double

// External includes

// Project includes
#include "includes/define.h"
#include "geometries/geometry.h"
#include "geometries/triangle_3d_3.h"
#include "geometries/quadrilateral_3d_4.h"
#include "integration/pyramid_gauss_legendre_integration_points.h"
#include "utilities/geometry_utilities.h"

namespace Kratos {

///@name Kratos Classes
///@{

/**
 * @class Pyramid3D5
 * @ingroup KratosCore
 * @brief A five node pyramid geometry with linear shape functions
 * @details The node ordering corresponds with:
 *                     4
 *                   ,/|\
 *                 ,/ .'|\
 *               ,/   μ | \
 *             ,/    .^ | `.
 *           ,/      || '.  \
 *         ,/       .'|  |   \
 *       ,/         | +--|----\-----> η
 *      0----------.'-`\ -3    `.
 *       `\        |    `\ `\    \
 *         `\     .'      `\ `\   \
 *           `\   |         ξ  `\  \
 *             `\.'              `\`\
 *                1------------------2
 *
 *
 * @author Philipp Bucher, Ashish Darekar
 */
template<class TPointType>
class Pyramid3D5 : public Geometry<TPointType>
{
public:
    ///@}
    ///@name Type Definitions
    ///@{

    /// Geometry as base class.
    typedef Geometry<TPointType> BaseType;

    /**
     * Type of edge and face geometries
     */
    typedef Line3D2<TPointType> EdgeType;
    typedef Triangle3D3<TPointType> FaceType1;
    typedef Quadrilateral3D4<TPointType> FaceType2;

    /// Pointer definition of Pyramid3D5
    KRATOS_CLASS_POINTER_DEFINITION(Pyramid3D5);

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

    explicit Pyramid3D5(
        typename PointType::Pointer pPoint1,
        typename PointType::Pointer pPoint2,
        typename PointType::Pointer pPoint3,
        typename PointType::Pointer pPoint4,
        typename PointType::Pointer pPoint5)
        : BaseType( PointsArrayType(), &msGeometryData )
    {
        this->Points().reserve(5);
        this->Points().push_back(pPoint1);
        this->Points().push_back(pPoint2);
        this->Points().push_back(pPoint3);
        this->Points().push_back(pPoint4);
        this->Points().push_back(pPoint5);
    }

    explicit Pyramid3D5( const PointsArrayType& ThisPoints )
        : BaseType( ThisPoints, &msGeometryData )
    {
        KRATOS_ERROR_IF( this->PointsNumber() != 5 ) << "Invalid points number. Expected 5, given " << this->PointsNumber() << std::endl;
    }

    /// Constructor with Geometry Id
    explicit Pyramid3D5(
        const IndexType GeometryId,
        const PointsArrayType& rThisPoints
    ) : BaseType( GeometryId, rThisPoints, &msGeometryData)
    {
        KRATOS_ERROR_IF( this->PointsNumber() != 5 ) << "Invalid points number. Expected 5, given " << this->PointsNumber() << std::endl;
    }

    /// Constructor with Geometry Name
    explicit Pyramid3D5(
        const std::string& rGeometryName,
        const PointsArrayType& rThisPoints
    ) : BaseType(rGeometryName, rThisPoints, &msGeometryData)
    {
        KRATOS_ERROR_IF(this->PointsNumber() != 5) << "Invalid points number. Expected 5, given " << this->PointsNumber() << std::endl;
    }

    /** Copy constructor.
     Construct this geometry as a copy of given geometry.

    @note This copy constructor don't copy the points and new
    geometry shares points with given source geometry. It's
    obvious that any change to this new geometry's point affect
    source geometry's points too.
    */
    Pyramid3D5(Pyramid3D5 const& rOther)
    : BaseType(rOther)
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
    template<class TOtherPointType> Pyramid3D5(Pyramid3D5<TOtherPointType> const& rOther)
    : BaseType(rOther)
    {
    }

    GeometryData::KratosGeometryFamily GetGeometryFamily() const override
    {
        return GeometryData::KratosGeometryFamily::Kratos_Pyramid;
    }

    GeometryData::KratosGeometryType GetGeometryType() const override
    {
        return GeometryData::KratosGeometryType::Kratos_Pyramid3D5;
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
    Pyramid3D5& operator=(const Pyramid3D5& rOther)
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
    Pyramid3D5& operator=(Pyramid3D5<TOtherPointType> const & rOther)
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
        return typename BaseType::Pointer( new Pyramid3D5( NewGeometryId, rThisPoints ) );
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
        auto p_geometry = typename BaseType::Pointer( new Pyramid3D5( NewGeometryId, rGeometry.Points() ) );
        p_geometry->SetData(rGeometry.GetData());
        return p_geometry;
    }

    ///@}
    ///@name Informations
    ///@{

    /**
     * @brief This method gives you number of all edges of this geometry.
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
        return 8;
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
        typedef typename Geometry<TPointType>::Pointer EdgePointerType;
        edges.push_back( EdgePointerType( new EdgeType(
                                              this->pGetPoint( 0 ),
                                              this->pGetPoint( 1 ) ) ) );
        edges.push_back( EdgePointerType( new EdgeType(
                                              this->pGetPoint( 1 ),
                                              this->pGetPoint( 2 ) ) ) );
        edges.push_back( EdgePointerType( new EdgeType(
                                              this->pGetPoint( 2 ),
                                              this->pGetPoint( 3 ) ) ) );
        edges.push_back( EdgePointerType( new EdgeType(
                                              this->pGetPoint( 3 ),
                                              this->pGetPoint( 0 ) ) ) );
        edges.push_back( EdgePointerType( new EdgeType(
                                              this->pGetPoint( 0 ),
                                              this->pGetPoint( 4 ) ) ) );
        edges.push_back( EdgePointerType( new EdgeType(
                                              this->pGetPoint( 1 ),
                                              this->pGetPoint( 4 ) ) ) );
        edges.push_back( EdgePointerType( new EdgeType(
                                              this->pGetPoint( 2 ),
                                              this->pGetPoint( 4 ) ) ) );
        edges.push_back( EdgePointerType( new EdgeType(
                                              this->pGetPoint( 3 ),
                                              this->pGetPoint( 4 ) ) ) );
        return edges;
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
     * @brief Returns all faces of the current geometry.
     * @details This is only implemented for 3D geometries, since 2D geometries only have edges but no faces
     * @return GeometriesArrayType containes this geometry faces.
     * @see EdgesNumber
     * @see GenerateEdges
     * @see FacesNumber
     */
    GeometriesArrayType GenerateFaces() const override
    {
        GeometriesArrayType faces = GeometriesArrayType();
        typedef typename Geometry<TPointType>::Pointer FacePointerType;
        faces.push_back( FacePointerType( new FaceType1(
                                              this->pGetPoint( 0 ),
                                              this->pGetPoint( 1 ),
                                              this->pGetPoint( 4 ) ) ) );
        faces.push_back( FacePointerType( new FaceType1(
                                              this->pGetPoint( 1 ),
                                              this->pGetPoint( 2 ),
                                              this->pGetPoint( 4 ) ) ) );
        faces.push_back( FacePointerType( new FaceType2(
                                              this->pGetPoint( 0 ),
                                              this->pGetPoint( 1 ),
                                              this->pGetPoint( 2 ),
                                              this->pGetPoint( 3 ) ) ) );
        faces.push_back( FacePointerType( new FaceType1(
                                              this->pGetPoint( 2 ),
                                              this->pGetPoint( 3 ),
                                              this->pGetPoint( 4 ) ) ) );
        faces.push_back( FacePointerType( new FaceType1(
                                              this->pGetPoint( 3 ),
                                              this->pGetPoint( 0 ),
                                              this->pGetPoint( 4 ) ) ) );
        return faces;
    }

    /** This method calculate and return volume of this
     geometry. For one and two dimensional geometry it returns
    zero and for three dimensional it gives volume of geometry.

    @return double value contains volume.
    @see Length()
    @see Area()
    @see DomainSize()
    */
    double Volume() const override
    {
        return IntegrationUtilities::ComputeVolume3DGeometry(*this);
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
        if ( rResult.size1() != 5 || rResult.size2() != 3 )
            rResult.resize( 5, 3, false );

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
        if(rResult.size() != 5) rResult.resize(5,false);

        for (std::size_t i=0; i<5; ++i) {
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
        const std::size_t points_number = 5;
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
        rResult.resize( 5, 3, false );
        noalias( rResult ) = ZeroMatrix( 5, 3 );

        rResult( 0, 0 ) =  (-0.125) * ( 1 - rPoint[1] ) * ( 1 - rPoint[2] ) ;
        rResult( 0, 1 ) =  (-0.125) * ( 1 - rPoint[0] ) * ( 1 - rPoint[2] ) ;
        rResult( 0, 2 ) =  (-0.125) * ( 1 - rPoint[0] ) * ( 1 - rPoint[1] ) ;

        rResult( 1, 0 ) =  (+0.125) * ( 1 - rPoint[1] ) * ( 1 - rPoint[2] ) ;
        rResult( 1, 1 ) =  (-0.125) * ( 1 + rPoint[0] ) * ( 1 - rPoint[2] ) ;
        rResult( 1, 2 ) =  (-0.125) * ( 1 + rPoint[0] ) * ( 1 - rPoint[1] ) ;

        rResult( 2, 0 ) =  (+0.125) * ( 1 + rPoint[1] ) * ( 1 - rPoint[2] ) ;
        rResult( 2, 1 ) =  (+0.125) * ( 1 + rPoint[0] ) * ( 1 - rPoint[2] ) ;
        rResult( 2, 2 ) =  (-0.125) * ( 1 + rPoint[0] ) * ( 1 + rPoint[1] ) ;

        rResult( 3, 0 ) =  (-0.125) * ( 1 + rPoint[1] ) * ( 1 - rPoint[2] ) ;
        rResult( 3, 1 ) =  (+0.125) * ( 1 - rPoint[0] ) * ( 1 - rPoint[2] ) ;
        rResult( 3, 2 ) =  (-0.125) * ( 1 - rPoint[0] ) * ( 1 + rPoint[1] ) ;

        rResult( 4, 0 ) =   0.00 ;
        rResult( 4, 1 ) =   0.00 ;
        rResult( 4, 2 ) =  +0.50 ;

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
        // Point to compute distance to
        const Point point(rPointGlobalCoordinates);

        // Check if point is inside
        CoordinatesArrayType aux_coordinates;
        if (this->IsInside(rPointGlobalCoordinates, aux_coordinates, Tolerance)) {
            return 0.0;
        }

        // Compute distance to faces
        std::array<double, 5> distances;
        distances[0]  = GeometryUtils::PointDistanceToTriangle3D(this->GetPoint(0), this->GetPoint(1), this->GetPoint(4), point);
        distances[1]  = GeometryUtils::PointDistanceToTriangle3D(this->GetPoint(1), this->GetPoint(2), this->GetPoint(4), point);
        distances[2]  = GeometryUtils::PointDistanceToQuadrilateral3D(this->GetPoint(0), this->GetPoint(1), this->GetPoint(2), this->GetPoint(3), point);
        distances[3]  = GeometryUtils::PointDistanceToTriangle3D(this->GetPoint(2), this->GetPoint(3), this->GetPoint(4), point);
        distances[4]  = GeometryUtils::PointDistanceToTriangle3D(this->GetPoint(3), this->GetPoint(0), this->GetPoint(4), point);
        const auto min = std::min_element(distances.begin(), distances.end());
        return *min;
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
        return "3 dimensional pyramid with 5 nodes in 3D space";
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

    Pyramid3D5() : BaseType( PointsArrayType(), &msGeometryData ) {}

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
            return( (0.125) * (1 - rPoint[0]) * (1 - rPoint[1]) * (1 - rPoint[2]) );
        case 1:
            return( (0.125) * (1 + rPoint[0]) * (1 - rPoint[1]) * (1 - rPoint[2]) );
        case 2:
            return( (0.125) * (1 + rPoint[0]) * (1 + rPoint[1]) * (1 - rPoint[2]) );
        case 3:
            return( (0.125) * (1 - rPoint[0]) * (1 + rPoint[1]) * (1 - rPoint[2]) );
        case 4:
            return( (0.5) * (1 + rPoint[2]) );
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
                Pyramid3D5<TPointType>::CalculateShapeFunctionsIntegrationPointsValues(GeometryData::IntegrationMethod::GI_GAUSS_1),
                Pyramid3D5<TPointType>::CalculateShapeFunctionsIntegrationPointsValues(GeometryData::IntegrationMethod::GI_GAUSS_2),
                Pyramid3D5<TPointType>::CalculateShapeFunctionsIntegrationPointsValues(GeometryData::IntegrationMethod::GI_GAUSS_3),
                Pyramid3D5<TPointType>::CalculateShapeFunctionsIntegrationPointsValues(GeometryData::IntegrationMethod::GI_GAUSS_4),
                Pyramid3D5<TPointType>::CalculateShapeFunctionsIntegrationPointsValues(GeometryData::IntegrationMethod::GI_GAUSS_5)
            }
        };
        return shape_functions_values;
    }

    static const ShapeFunctionsLocalGradientsContainerType AllShapeFunctionsLocalGradients()
    {
        ShapeFunctionsLocalGradientsContainerType shape_functions_local_gradients =
        {
            {
                Pyramid3D5<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients(GeometryData::IntegrationMethod::GI_GAUSS_1),
                Pyramid3D5<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients(GeometryData::IntegrationMethod::GI_GAUSS_2),
                Pyramid3D5<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients(GeometryData::IntegrationMethod::GI_GAUSS_3),
                Pyramid3D5<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients(GeometryData::IntegrationMethod::GI_GAUSS_4),
                Pyramid3D5<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients(GeometryData::IntegrationMethod::GI_GAUSS_5)
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

    template<class TOtherPointType> friend class Pyramid3D5;

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
                    Pyramid3D5<TPointType>& rThis);

/// output stream function
template<class TPointType>
inline std::ostream& operator << (std::ostream& rOStream,
                    const Pyramid3D5<TPointType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

template<class TPointType> const
GeometryData Pyramid3D5<TPointType>::msGeometryData(
    &msGeometryDimension,
    GeometryData::IntegrationMethod::GI_GAUSS_2,
    Pyramid3D5<TPointType>::AllIntegrationPoints(),
    Pyramid3D5<TPointType>::AllShapeFunctionsValues(),
    AllShapeFunctionsLocalGradients()
);

template<class TPointType> const
GeometryDimension Pyramid3D5<TPointType>::msGeometryDimension(3, 3);

}  // namespace Kratos.
