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
//                   Bodhinanda Chandra
//

#pragma once

// System includes

// External includes

// Project includes
#include "geometries/quadrilateral_3d_4.h"
#include "utilities/integration_utilities.h"
#include "utilities/geometry_utilities.h"
#include "integration/hexahedron_gauss_legendre_integration_points.h"
#include "integration/hexahedron_gauss_lobatto_integration_points.h"

namespace Kratos
{
/**
 * @class Hexahedra3D8
 * @ingroup KratosCore
 * @brief An eight node hexahedra geometry with linear shape functions
 * @details The node ordering corresponds with:
 *             v
 *      3----------2
 *      |\     ^   |\
 *      | \    |   | \
 *      |  \   |   |  \
 *      |   7------+---6
 *      |   |  +-- |-- | -> u
 *      0---+---\--1   |
 *       \  |    \  \  |
 *        \ |     \  \ |
 *         \|      w  \|
 *          4----------5
 * @author Riccardo Rossi
 * @author Janosch Stascheit
 * @author Felix Nagel
 */
template<class TPointType> class Hexahedra3D8 : public Geometry<TPointType>
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
    typedef Line3D2<TPointType> EdgeType;
    typedef Quadrilateral3D4<TPointType> FaceType;

    /**
     * Pointer definition of Hexahedra3D8
     */
    KRATOS_CLASS_POINTER_DEFINITION( Hexahedra3D8 );

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
     * A third order tensor to hold shape functions' local second derivatives.
     * ShapefunctionsLocalGradients function return this
     * type as its result.
     */
    typedef typename BaseType::ShapeFunctionsSecondDerivativesType ShapeFunctionsSecondDerivativesType;

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
    Hexahedra3D8( const PointType& Point1, const PointType& Point2,
                  const PointType& Point3, const PointType& Point4,
                  const PointType& Point5, const PointType& Point6,
                  const PointType& Point7, const PointType& Point8 )
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
    }

    Hexahedra3D8( typename PointType::Pointer pPoint1,
                  typename PointType::Pointer pPoint2,
                  typename PointType::Pointer pPoint3,
                  typename PointType::Pointer pPoint4,
                  typename PointType::Pointer pPoint5,
                  typename PointType::Pointer pPoint6,
                  typename PointType::Pointer pPoint7,
                  typename PointType::Pointer pPoint8 )
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
    }

    explicit Hexahedra3D8( const PointsArrayType& ThisPoints )
        : BaseType( ThisPoints, &msGeometryData )
    {
        if ( this->PointsNumber() != 8 )
            KRATOS_ERROR << "Invalid points number. Expected 8, given " << this->PointsNumber() << std::endl;
    }

    /// Constructor with Geometry Id
    explicit Hexahedra3D8(
        const IndexType GeometryId,
        const PointsArrayType& rThisPoints
        ) : BaseType( GeometryId, rThisPoints, &msGeometryData )
    {
        KRATOS_ERROR_IF( this->PointsNumber() != 8 ) << "Invalid points number. Expected 8, given " << this->PointsNumber() << std::endl;
    }

    /// Constructor with Geometry Name
    explicit Hexahedra3D8(
        const std::string& rGeometryName,
        const PointsArrayType& rThisPoints
    ) : BaseType( rGeometryName, rThisPoints, &msGeometryData)
    {
        KRATOS_ERROR_IF(this->PointsNumber() != 8) << "Invalid points number. Expected 8, given " << this->PointsNumber() << std::endl;
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
    Hexahedra3D8( Hexahedra3D8 const& rOther )
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
    template<class TOtherPointType> explicit Hexahedra3D8( Hexahedra3D8<TOtherPointType> const& rOther )
        : BaseType( rOther )
    {
    }

    /// Destructor. Does nothing!!!
    ~Hexahedra3D8() override {}

    /**
     * @brief Gets the geometry family.
     * @details This function returns the family type of the geometry. The geometry family categorizes the geometry into a broader classification, aiding in its identification and processing.
     * @return GeometryData::KratosGeometryFamily The geometry family.
     */
    GeometryData::KratosGeometryFamily GetGeometryFamily() const override
    {
        return GeometryData::KratosGeometryFamily::Kratos_Hexahedra;
    }

    /**
     * @brief Gets the geometry type.
     * @details This function returns the specific type of the geometry. The geometry type provides a more detailed classification of the geometry.
     * @return GeometryData::KratosGeometryType The specific geometry type.
     */
    GeometryData::KratosGeometryType GetGeometryType() const override
    {
        return GeometryData::KratosGeometryType::Kratos_Hexahedra3D8;
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
    Hexahedra3D8& operator=( const Hexahedra3D8& rOther )
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
    Hexahedra3D8& operator=( Hexahedra3D8<TOtherPointType> const & rOther )
    {
        BaseType::operator=( rOther );

        return *this;
    }


    /**
     * Operations
     */

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
        return typename BaseType::Pointer( new Hexahedra3D8( NewGeometryId, rThisPoints ) );
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
        auto p_geometry = typename BaseType::Pointer( new Hexahedra3D8( NewGeometryId, rGeometry.Points() ) );
        p_geometry->SetData(rGeometry.GetData());
        return p_geometry;
    }

    /**
     * Information
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
    double Length() const override
    {
        return std::sqrt( std::abs( this->DeterminantOfJacobian( PointType() ) ) );
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
    double Area() const override
    {
        return Volume();
    }

    /**
     * @brief This method calculate and return volume of this geometry.
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
    double DomainSize() const override
    {
        return Volume();
    }

    /**
     * Returns a matrix of the local coordinates of all points
     * @param rResult a Matrix that will be overwritten by the results
     * @return the coordinates of all points of the current geometry
     */
    Matrix& PointsLocalCoordinates( Matrix& rResult ) const override
    {
        if ( rResult.size1() != 8 || rResult.size2() != 3 )
            rResult.resize( 8, 3, false );

        rResult( 0, 0 ) = -1.0;
        rResult( 0, 1 ) = -1.0;
        rResult( 0, 2 ) = -1.0;

        rResult( 1, 0 ) = 1.0;
        rResult( 1, 1 ) = -1.0;
        rResult( 1, 2 ) = -1.0;

        rResult( 2, 0 ) = 1.0;
        rResult( 2, 1 ) = 1.0;
        rResult( 2, 2 ) = -1.0;

        rResult( 3, 0 ) = -1.0;
        rResult( 3, 1 ) = 1.0;
        rResult( 3, 2 ) = -1.0;

        rResult( 4, 0 ) = -1.0;
        rResult( 4, 1 ) = -1.0;
        rResult( 4, 2 ) = 1.0;

        rResult( 5, 0 ) = 1.0;
        rResult( 5, 1 ) = -1.0;
        rResult( 5, 2 ) = 1.0;

        rResult( 6, 0 ) = 1.0;
        rResult( 6, 1 ) = 1.0;
        rResult( 6, 2 ) = 1.0;

        rResult( 7, 0 ) = -1.0;
        rResult( 7, 1 ) = 1.0;
        rResult( 7, 2 ) = 1.0;

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
        ) const override
    {
        this->PointLocalCoordinates( rResult, rPoint );

        if ( std::abs( rResult[0] ) <= (1.0 + Tolerance) ) {
            if ( std::abs( rResult[1] ) <= (1.0 + Tolerance) ) {
                if ( std::abs( rResult[2] ) <= (1.0 + Tolerance) ) {
                    return true;
                }
            }
        }

        return false;
    }


    /** This method gives you number of all edges of this
    geometry.
    @return SizeType containes number of this geometry edges.
    @see Edges()
    @see Edge()
     */
    // will be used by refinement algorithm, thus uncommented. janosch.
    SizeType EdgesNumber() const override
    {
        return 12;
    }

    /** This method gives you all edges of this geometry.

    @return GeometriesArrayType containes this geometry edges.
    @see EdgesNumber()
    @see Edge()
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
                                              this->pGetPoint( 4 ),
                                              this->pGetPoint( 5 ) ) ) );
        edges.push_back( EdgePointerType( new EdgeType(
                                              this->pGetPoint( 5 ),
                                              this->pGetPoint( 6 ) ) ) );
        edges.push_back( EdgePointerType( new EdgeType(
                                              this->pGetPoint( 6 ),
                                              this->pGetPoint( 7 ) ) ) );
        edges.push_back( EdgePointerType( new EdgeType(
                                              this->pGetPoint( 7 ),
                                              this->pGetPoint( 4 ) ) ) );
        edges.push_back( EdgePointerType( new EdgeType(
                                              this->pGetPoint( 0 ),
                                              this->pGetPoint( 4 ) ) ) );
        edges.push_back( EdgePointerType( new EdgeType(
                                              this->pGetPoint( 1 ),
                                              this->pGetPoint( 5 ) ) ) );
        edges.push_back( EdgePointerType( new EdgeType(
                                              this->pGetPoint( 2 ),
                                              this->pGetPoint( 6 ) ) ) );
        edges.push_back( EdgePointerType( new EdgeType(
                                              this->pGetPoint( 3 ),
                                              this->pGetPoint( 7 ) ) ) );
        return edges;
    }

    /** This method calculates and returns the average edge length of the geometry
    *
    * @return double value with the average edge length
    *
    */
    double AverageEdgeLength() const override
    {
        const TPointType& p0 = this->GetPoint(0);
        const TPointType& p1 = this->GetPoint(1);
        const TPointType& p2 = this->GetPoint(2);
        const TPointType& p3 = this->GetPoint(3);
        const TPointType& p4 = this->GetPoint(4);
        const TPointType& p5 = this->GetPoint(5);
        const TPointType& p6 = this->GetPoint(6);
        const TPointType& p7 = this->GetPoint(7);
        return (MathUtils<double>::Norm3(p0-p1) +
            MathUtils<double>::Norm3(p1-p2) +
            MathUtils<double>::Norm3(p2-p3) +
            MathUtils<double>::Norm3(p3-p0) +
            MathUtils<double>::Norm3(p4-p5) +
            MathUtils<double>::Norm3(p5-p6) +
            MathUtils<double>::Norm3(p6-p7) +
            MathUtils<double>::Norm3(p7-p4) +
            MathUtils<double>::Norm3(p0-p4) +
            MathUtils<double>::Norm3(p1-p5) +
            MathUtils<double>::Norm3(p2-p6) +
            MathUtils<double>::Norm3(p3-p7)) /12.0;
    }

    /** This method calculates and returns the maximum edge
     * length of the geometry
     *
     * @return double value with the maximum edge length
     *
     * @see MinEdgeLength()
     * @see AverageEdgeLength()
     */
    double MaxEdgeLength() const override {
        const auto edges = GenerateEdges();
        double max_edge_length = 0.0;
        for (const auto& r_edge: edges) {
            max_edge_length = std::max(max_edge_length, r_edge.Length());
        }
        return max_edge_length;
    }

    /** This method calculates and returns the minimum edge
     * length of the geometry
     *
     * @return double value with the maximum edge length
     *
     * @see MaxEdgeLength()
     * @see AverageEdgeLength()
     */
    double MinEdgeLength() const override {
        const auto edges = GenerateEdges();
        double min_edge_length = std::numeric_limits<double>::max();
        for (const auto& r_edge: edges) {
            min_edge_length = std::min(min_edge_length, r_edge.Length());
        }
        return min_edge_length;
    }

    ///@}
    ///@name Face
    ///@{

    /**
     * @brief Returns the number of faces of the current geometry.
     * @details This is only implemented for 3D geometries, since 2D geometries only have edges but no faces
     * @see EdgesNumber
     * @see Edges
     * @see Faces
     */
    SizeType FacesNumber() const override
    {
        return 6;
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
        faces.push_back( FacePointerType( new FaceType(
                                              this->pGetPoint( 3 ),
                                              this->pGetPoint( 2 ),
                                              this->pGetPoint( 1 ),
                                              this->pGetPoint( 0 ) ) ) );
        faces.push_back( FacePointerType( new FaceType(
                                              this->pGetPoint( 0 ),
                                              this->pGetPoint( 1 ),
                                              this->pGetPoint( 5 ),
                                              this->pGetPoint( 4 ) ) ) );
        faces.push_back( FacePointerType( new FaceType(
                                              this->pGetPoint( 2 ),
                                              this->pGetPoint( 6 ),
                                              this->pGetPoint( 5 ),
                                              this->pGetPoint( 1 ) ) ) );
        faces.push_back( FacePointerType( new FaceType(
                                              this->pGetPoint( 7 ),
                                              this->pGetPoint( 6 ),
                                              this->pGetPoint( 2 ),
                                              this->pGetPoint( 3 ) ) ) );
        faces.push_back( FacePointerType( new FaceType(
                                              this->pGetPoint( 7 ),
                                              this->pGetPoint( 3 ),
                                              this->pGetPoint( 0 ),
                                              this->pGetPoint( 4 ) ) ) );
        faces.push_back( FacePointerType( new FaceType(
                                              this->pGetPoint( 4 ),
                                              this->pGetPoint( 5 ),
                                              this->pGetPoint( 6 ),
                                              this->pGetPoint( 7 ) ) ) );
        return faces;
    }


    bool HasIntersection( const Point& rLowPoint, const Point& rHighPoint ) const override
    {
        using Quadrilateral3D4Type = Quadrilateral3D4<TPointType>;
        // Check if faces have intersection
        if(Quadrilateral3D4Type(this->pGetPoint(3),this->pGetPoint(2), this->pGetPoint(1), this->pGetPoint(0)).HasIntersection(rLowPoint, rHighPoint))
            return true;
        if(Quadrilateral3D4Type(this->pGetPoint(0),this->pGetPoint(1), this->pGetPoint(5), this->pGetPoint(4)).HasIntersection(rLowPoint, rHighPoint))
            return true;
        if(Quadrilateral3D4Type(this->pGetPoint(2),this->pGetPoint(6), this->pGetPoint(5), this->pGetPoint(1)).HasIntersection(rLowPoint, rHighPoint))
            return true;
        if(Quadrilateral3D4Type(this->pGetPoint(7),this->pGetPoint(6), this->pGetPoint(2), this->pGetPoint(3)).HasIntersection(rLowPoint, rHighPoint))
            return true;
        if(Quadrilateral3D4Type(this->pGetPoint(7),this->pGetPoint(3), this->pGetPoint(0), this->pGetPoint(4)).HasIntersection(rLowPoint, rHighPoint))
            return true;
        if(Quadrilateral3D4Type(this->pGetPoint(4),this->pGetPoint(5), this->pGetPoint(6), this->pGetPoint(7)).HasIntersection(rLowPoint, rHighPoint))
            return true;

        CoordinatesArrayType local_coordinates;
        // if there are no faces intersecting the box then or the box is inside the hexahedron or it does not have intersection
        if(IsInside(rLowPoint,local_coordinates))
            return true;

        return false;
    }

    /**
     * @brief Implements the calculus of the 8 solid angles of the hexa
     * @details Implements the calculus of the 8 solid angles of the hexa
     * @param rSolidAngles The solid angles of the geometry
     */
    void ComputeSolidAngles(Vector& rSolidAngles) const override
    {
        if(rSolidAngles.size() != 8) {
          rSolidAngles.resize(8, false);
        }

        Vector dihedral_angles(24);
        ComputeDihedralAngles(dihedral_angles);

        for (unsigned int i = 0; i < 8; ++i) {
            rSolidAngles[i] = dihedral_angles[3*i]
              + dihedral_angles[3*i + 1]
              + dihedral_angles[3*i + 2]
              - Globals::Pi;
        }
    }

    /**
     * @brief Implements the calculus of the 24 dihedral angles of the hexa
     * @details Implements the calculus of the 24 dihedral angles of the hexa.
     * Each edge has two different dihedral angles in each extreme.
     * @param rDihedralAngles The dihedral angles
     */
    void ComputeDihedralAngles(Vector& rDihedralAngles) const override
    {
        if(rDihedralAngles.size() != 24) {
          rDihedralAngles.resize(24, false);
        }
        const auto faces = this->GenerateFaces();
        // The three faces that contain the node i can be obtained by doing:
        // faces[faces_0[i]], faces[faces_1[i]], faces[faces_2[i]]
        const std::array<unsigned int, 8> faces_0 = {0,0,0,0,5,5,5,5};
        const std::array<unsigned int, 8> faces_1 = {1,1,3,3,1,1,3,3};
        const std::array<unsigned int, 8> faces_2 = {4,2,2,4,4,2,2,4};

        array_1d<double, 3> normal_0, normal_1, normal_2;
        double dihedral_angle_0, dihedral_angle_1, dihedral_angle_2;
        for (unsigned int i = 0; i < 8; ++i) {
            const TPointType& r_point_i = this->GetPoint(i);
            noalias(normal_0) = faces[faces_0[i]].UnitNormal(r_point_i);
            noalias(normal_1) = faces[faces_1[i]].UnitNormal(r_point_i);
            noalias(normal_2) = faces[faces_2[i]].UnitNormal(r_point_i);
            dihedral_angle_0 = std::acos(inner_prod(normal_0, -normal_1));
            dihedral_angle_1 = std::acos(inner_prod(normal_0, -normal_2));
            dihedral_angle_2 = std::acos(inner_prod(normal_2, -normal_1));
            rDihedralAngles[i*3] = dihedral_angle_0;
            rDihedralAngles[i*3 + 1] = dihedral_angle_1;
            rDihedralAngles[i*3 + 2] = dihedral_angle_2;
        }
    }

    /**
     * @brief Calculates the min dihedral angle quality metric.
     * @details Calculates the min dihedral angle quality metric.
     * The min dihedral angle is min angle between two faces of the element
     * in radians
     * @return [description]
     */
    double MinDihedralAngle() const override {
      Vector dihedral_angles(24);
      ComputeDihedralAngles(dihedral_angles);
      double min_dihedral_angle = dihedral_angles[0];
      for (unsigned int i = 1; i < 24; i++) {
        min_dihedral_angle = std::min(dihedral_angles[i], min_dihedral_angle);
      }
      return min_dihedral_angle;
    }

    /**
     * @brief Calculates the max dihedral angle quality metric.
     * @details Calculates the max dihedral angle quality metric.
     * The max dihedral angle is max angle between two faces of the element
     * in radians
     * @return [description]
     */
    double MaxDihedralAngle() const override {
        Vector dihedral_angles(24);
        ComputeDihedralAngles(dihedral_angles);
        double max_dihedral_angle = dihedral_angles[0];
        for (unsigned int i = 1; i < 24; i++) {
            max_dihedral_angle = std::max(dihedral_angles[i], max_dihedral_angle);
        }
        return max_dihedral_angle;
    }

    /**
     * @brief Calculates the volume to average edge length quality metric.
     * @details The average edge length is calculated using the RMS.
     * This metric is bounded by the interval (0,1) being:
     *  1 -> Optimal value
     *  0 -> Worst value
     *
     * \f$ \frac{V}{\sqrt{\frac{1}{6}\sum{A_{i}^{2}}}} \f$
     *
     * @return [description]
     */
    double VolumeToRMSEdgeLength() const override {
        const auto edges = GenerateEdges();
        double sum_squared_lengths = 0.0;
        for (const auto& r_edge : edges) {
            const double length = r_edge.Length();
            sum_squared_lengths += length*length;
        }

        const double rms_edge = std::sqrt(1.0/12.0 * sum_squared_lengths);

        return Volume() / std::pow(rms_edge, 3.0);
    }

    /**
     * @brief Calculates the shortest to longest edge quality metric.
     * Calculates the shortest to longest edge quality metric.
     * This metric is bounded by the interval (0,1) being:
     *  1 -> Optimal value
     *  0 -> Worst value
     *
     * \f$ \frac{l}{L} \f$
     *
     * @return [description]
     */
    double ShortestToLongestEdgeQuality() const override {
        const auto edges = GenerateEdges();
        double min_edge_length = std::numeric_limits<double>::max();
        double max_edge_length = -std::numeric_limits<double>::max();
        for (const auto& r_edge: edges) {
            min_edge_length = std::min(min_edge_length, r_edge.Length());
            max_edge_length = std::max(max_edge_length, r_edge.Length());
        }
        return min_edge_length / max_edge_length;
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
    double ShapeFunctionValue( IndexType ShapeFunctionIndex,
                                       const CoordinatesArrayType& rPoint ) const override
    {
        switch ( ShapeFunctionIndex )
        {
        case 0:
            return( 0.125*( 1.0 - rPoint[0] )*( 1.0 - rPoint[1] )*( 1.0 - rPoint[2] ) );
        case 1:
            return( 0.125*( 1.0 + rPoint[0] )*( 1.0 - rPoint[1] )*( 1.0 - rPoint[2] ) );
        case 2:
            return( 0.125*( 1.0 + rPoint[0] )*( 1.0 + rPoint[1] )*( 1.0 - rPoint[2] ) );
        case 3:
            return( 0.125*( 1.0 - rPoint[0] )*( 1.0 + rPoint[1] )*( 1.0 - rPoint[2] ) );
        case 4:
            return( 0.125*( 1.0 - rPoint[0] )*( 1.0 - rPoint[1] )*( 1.0 + rPoint[2] ) );
        case 5:
            return( 0.125*( 1.0 + rPoint[0] )*( 1.0 - rPoint[1] )*( 1.0 + rPoint[2] ) );
        case 6:
            return( 0.125*( 1.0 + rPoint[0] )*( 1.0 + rPoint[1] )*( 1.0 + rPoint[2] ) );
        case 7:
            return( 0.125*( 1.0 - rPoint[0] )*( 1.0 + rPoint[1] )*( 1.0 + rPoint[2] ) );
        default:
            KRATOS_ERROR << "Wrong index of shape function!" << *this << std::endl;
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
    Vector& ShapeFunctionsValues (Vector &rResult, const CoordinatesArrayType& rCoordinates) const override
    {
      if(rResult.size() != 8) rResult.resize(8,false);
      rResult[0] =  0.125*( 1.0 - rCoordinates[0] )*( 1.0 - rCoordinates[1] )*( 1.0 - rCoordinates[2] ) ;
      rResult[1] =  0.125*( 1.0 + rCoordinates[0] )*( 1.0 - rCoordinates[1] )*( 1.0 - rCoordinates[2] ) ;
      rResult[2] =  0.125*( 1.0 + rCoordinates[0] )*( 1.0 + rCoordinates[1] )*( 1.0 - rCoordinates[2] ) ;
      rResult[3] =  0.125*( 1.0 - rCoordinates[0] )*( 1.0 + rCoordinates[1] )*( 1.0 - rCoordinates[2] ) ;
      rResult[4] =  0.125*( 1.0 - rCoordinates[0] )*( 1.0 - rCoordinates[1] )*( 1.0 + rCoordinates[2] ) ;
      rResult[5] =  0.125*( 1.0 + rCoordinates[0] )*( 1.0 - rCoordinates[1] )*( 1.0 + rCoordinates[2] ) ;
      rResult[6] =  0.125*( 1.0 + rCoordinates[0] )*( 1.0 + rCoordinates[1] )*( 1.0 + rCoordinates[2] ) ;
      rResult[7] =  0.125*( 1.0 - rCoordinates[0] )*( 1.0 + rCoordinates[1] )*( 1.0 + rCoordinates[2] ) ;
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
    Matrix& ShapeFunctionsLocalGradients( Matrix& rResult, const CoordinatesArrayType& rPoint ) const override
    {
        if ( rResult.size1() != 8 || rResult.size2() != 3 )
            rResult.resize( 8, 3, false );

        rResult( 0, 0 ) = -0.125 * ( 1.0 - rPoint[1] ) * ( 1.0 - rPoint[2] );
        rResult( 0, 1 ) = -0.125 * ( 1.0 - rPoint[0] ) * ( 1.0 - rPoint[2] );
        rResult( 0, 2 ) = -0.125 * ( 1.0 - rPoint[0] ) * ( 1.0 - rPoint[1] );
        rResult( 1, 0 ) =  0.125 * ( 1.0 - rPoint[1] ) * ( 1.0 - rPoint[2] );
        rResult( 1, 1 ) = -0.125 * ( 1.0 + rPoint[0] ) * ( 1.0 - rPoint[2] );
        rResult( 1, 2 ) = -0.125 * ( 1.0 + rPoint[0] ) * ( 1.0 - rPoint[1] );
        rResult( 2, 0 ) =  0.125 * ( 1.0 + rPoint[1] ) * ( 1.0 - rPoint[2] );
        rResult( 2, 1 ) =  0.125 * ( 1.0 + rPoint[0] ) * ( 1.0 - rPoint[2] );
        rResult( 2, 2 ) = -0.125 * ( 1.0 + rPoint[0] ) * ( 1.0 + rPoint[1] );
        rResult( 3, 0 ) = -0.125 * ( 1.0 + rPoint[1] ) * ( 1.0 - rPoint[2] );
        rResult( 3, 1 ) =  0.125 * ( 1.0 - rPoint[0] ) * ( 1.0 - rPoint[2] );
        rResult( 3, 2 ) = -0.125 * ( 1.0 - rPoint[0] ) * ( 1.0 + rPoint[1] );
        rResult( 4, 0 ) = -0.125 * ( 1.0 - rPoint[1] ) * ( 1.0 + rPoint[2] );
        rResult( 4, 1 ) = -0.125 * ( 1.0 - rPoint[0] ) * ( 1.0 + rPoint[2] );
        rResult( 4, 2 ) =  0.125 * ( 1.0 - rPoint[0] ) * ( 1.0 - rPoint[1] );
        rResult( 5, 0 ) =  0.125 * ( 1.0 - rPoint[1] ) * ( 1.0 + rPoint[2] );
        rResult( 5, 1 ) = -0.125 * ( 1.0 + rPoint[0] ) * ( 1.0 + rPoint[2] );
        rResult( 5, 2 ) =  0.125 * ( 1.0 + rPoint[0] ) * ( 1.0 - rPoint[1] );
        rResult( 6, 0 ) =  0.125 * ( 1.0 + rPoint[1] ) * ( 1.0 + rPoint[2] );
        rResult( 6, 1 ) =  0.125 * ( 1.0 + rPoint[0] ) * ( 1.0 + rPoint[2] );
        rResult( 6, 2 ) =  0.125 * ( 1.0 + rPoint[0] ) * ( 1.0 + rPoint[1] );
        rResult( 7, 0 ) = -0.125 * ( 1.0 + rPoint[1] ) * ( 1.0 + rPoint[2] );
        rResult( 7, 1 ) =  0.125 * ( 1.0 - rPoint[0] ) * ( 1.0 + rPoint[2] );
        rResult( 7, 2 ) =  0.125 * ( 1.0 - rPoint[0] ) * ( 1.0 + rPoint[1] );

        return rResult;
    }


    /**
     * Calculates the second derivatives in terms of local coordinateds
     * of all shape functions in a given point.
     *
     * @param rResult the given container will be overwritten by the results
     * @param rPoint the given local coordinates the derivatives will be evaluated for.
     * @return a third order tensor containing the second order derivatives for each shape function
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

        for ( unsigned int i = 0; i < this->PointsNumber(); ++i )
        {
            rResult[i].resize(3, 3, false);
        }

        rResult[0]( 0, 0 ) = 0.0;
        rResult[0]( 0, 1 ) = 0.125 * ( 1.0 - rPoint[2] );
        rResult[0]( 0, 2 ) = 0.125 * ( 1.0 - rPoint[1] );
        rResult[0]( 1, 0 ) = 0.125 * ( 1.0 - rPoint[2] );
        rResult[0]( 1, 1 ) = 0.0;
        rResult[0]( 1, 2 ) = 0.125 * ( 1.0 - rPoint[0] );
        rResult[0]( 2, 0 ) = 0.125 * ( 1.0 - rPoint[1] );
        rResult[0]( 2, 1 ) = 0.125 * ( 1.0 - rPoint[0] );
        rResult[0]( 2, 2 ) = 0.0;

        rResult[1]( 0, 0 ) = 0.0;
        rResult[1]( 0, 1 ) = -0.125 * ( 1.0 - rPoint[2] );
        rResult[1]( 0, 2 ) = -0.125 * ( 1.0 - rPoint[1] );
        rResult[1]( 1, 0 ) = -0.125 * ( 1.0 - rPoint[2] );
        rResult[1]( 1, 1 ) = 0.0;
        rResult[1]( 1, 2 ) = 0.125 * ( 1.0 + rPoint[0] );
        rResult[1]( 2, 0 ) = -0.125 * ( 1.0 - rPoint[1] );
        rResult[1]( 2, 1 ) = 0.125 * ( 1.0 + rPoint[0] );
        rResult[1]( 2, 2 ) = 0.0;

        rResult[2]( 0, 0 ) = 0.0;
        rResult[2]( 0, 1 ) = 0.125 * ( 1.0 - rPoint[2] );
        rResult[2]( 0, 2 ) = -0.125 * ( 1.0 + rPoint[1] );
        rResult[2]( 1, 0 ) = 0.125 * ( 1.0 - rPoint[2] );
        rResult[2]( 1, 1 ) = 0.0;
        rResult[2]( 1, 2 ) = -0.125 * ( 1.0 + rPoint[0] );
        rResult[2]( 2, 0 ) = -0.125 * ( 1.0 + rPoint[1] );
        rResult[2]( 2, 1 ) = -0.125 * ( 1.0 + rPoint[0] );
        rResult[2]( 2, 2 ) = 0.0;

        rResult[3]( 0, 0 ) = 0.0;
        rResult[3]( 0, 1 ) = -0.125 * ( 1.0 - rPoint[2] );
        rResult[3]( 0, 2 ) = 0.125 * ( 1.0 + rPoint[1] );
        rResult[3]( 1, 0 ) = -0.125 * ( 1.0 - rPoint[2] );
        rResult[3]( 1, 1 ) = 0.0;
        rResult[3]( 1, 2 ) = -0.125 * ( 1.0 - rPoint[0] );
        rResult[3]( 2, 0 ) = 0.125 * ( 1.0 + rPoint[1] );
        rResult[3]( 2, 1 ) = -0.125 * ( 1.0 - rPoint[0] );
        rResult[3]( 2, 2 ) = 0.0;

        rResult[4]( 0, 0 ) = 0.0;
        rResult[4]( 0, 1 ) = 0.125 * ( 1.0 + rPoint[2] );
        rResult[4]( 0, 2 ) = -0.125 * ( 1.0 - rPoint[1] );
        rResult[4]( 1, 0 ) = 0.125 * ( 1.0 + rPoint[2] );
        rResult[4]( 1, 1 ) = 0.0;
        rResult[4]( 1, 2 ) = -0.125 * ( 1.0 - rPoint[0] );
        rResult[4]( 2, 0 ) = -0.125 * ( 1.0 - rPoint[1] );
        rResult[4]( 2, 1 ) = -0.125 * ( 1.0 - rPoint[0] );
        rResult[4]( 2, 2 ) = 0.0;

        rResult[5]( 0, 0 ) = 0.0;
        rResult[5]( 0, 1 ) = -0.125 * ( 1.0 + rPoint[2] );
        rResult[5]( 0, 2 ) = 0.125 * ( 1.0 - rPoint[1] );
        rResult[5]( 1, 0 ) = -0.125 * ( 1.0 + rPoint[2] );
        rResult[5]( 1, 1 ) = 0.0;
        rResult[5]( 1, 2 ) = -0.125 * ( 1.0 + rPoint[0] );
        rResult[5]( 2, 0 ) = 0.125 * ( 1.0 - rPoint[1] );
        rResult[5]( 2, 1 ) = -0.125 * ( 1.0 + rPoint[0] );
        rResult[5]( 2, 2 ) = 0.0;

        rResult[6]( 0, 0 ) = 0.0;
        rResult[6]( 0, 1 ) = 0.125 * ( 1.0 + rPoint[2] );
        rResult[6]( 0, 2 ) = 0.125 * ( 1.0 + rPoint[1] );
        rResult[6]( 1, 0 ) = 0.125 * ( 1.0 + rPoint[2] );
        rResult[6]( 1, 1 ) = 0.0;
        rResult[6]( 1, 2 ) = 0.125 * ( 1.0 + rPoint[0] );
        rResult[6]( 2, 0 ) = 0.125 * ( 1.0 + rPoint[1] );
        rResult[6]( 2, 1 ) = 0.125 * ( 1.0 + rPoint[0] );
        rResult[6]( 2, 2 ) = 0.0;

        rResult[7]( 0, 0 ) = 0.0;
        rResult[7]( 0, 1 ) = -0.125 * ( 1.0 + rPoint[2] );
        rResult[7]( 0, 2 ) = -0.125 * ( 1.0 + rPoint[1] );
        rResult[7]( 1, 0 ) = -0.125 * ( 1.0 + rPoint[2] );
        rResult[7]( 1, 1 ) = 0.0;
        rResult[7]( 1, 2 ) = 0.125 * ( 1.0 - rPoint[0] );
        rResult[7]( 2, 0 ) = -0.125 * ( 1.0 + rPoint[1] );
        rResult[7]( 2, 1 ) = 0.125 * ( 1.0 - rPoint[0] );
        rResult[7]( 2, 2 ) = 0.0;

        return rResult;
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
        std::array<double, 6> distances;
        distances[0] = GeometryUtils::PointDistanceToQuadrilateral3D(this->GetPoint(3), this->GetPoint(2), this->GetPoint(1), this->GetPoint(0), point);
        distances[1] = GeometryUtils::PointDistanceToQuadrilateral3D(this->GetPoint(0), this->GetPoint(1), this->GetPoint(5), this->GetPoint(4), point);
        distances[2] = GeometryUtils::PointDistanceToQuadrilateral3D(this->GetPoint(2), this->GetPoint(6), this->GetPoint(5), this->GetPoint(1), point);
        distances[3] = GeometryUtils::PointDistanceToQuadrilateral3D(this->GetPoint(7), this->GetPoint(6), this->GetPoint(2), this->GetPoint(3), point);
        distances[4] = GeometryUtils::PointDistanceToQuadrilateral3D(this->GetPoint(7), this->GetPoint(3), this->GetPoint(0), this->GetPoint(4), point);
        distances[5] = GeometryUtils::PointDistanceToQuadrilateral3D(this->GetPoint(4), this->GetPoint(5), this->GetPoint(6), this->GetPoint(7), point);
        const auto min = std::min_element(distances.begin(), distances.end());
        return *min;
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
        return "3 dimensional hexahedra with eight nodes in 3D space";
    }

    /**
     * Print information about this object.
     *
     * @param rOStream Stream to print into it.
     * @see PrintData()
     * @see Info()
     */
    void PrintInfo( std::ostream& rOStream ) const override
    {
        rOStream << "3 dimensional hexahedra with eight nodes in 3D space";
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

protected:

    /**
     * there are no protected class members
     */

private:

    /**
     * Static Member Variables
     */

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

    Hexahedra3D8(): BaseType( PointsArrayType(), &msGeometryData ) {}

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
    static Matrix ShapeFunctionsLocalGradients( const CoordinatesArrayType& rPoint )
    {
        Matrix result = ZeroMatrix( 8, 3 );
        result( 0, 0 ) = -0.125 * ( 1.0 - rPoint[1] ) * ( 1.0 - rPoint[2] );
        result( 0, 1 ) = -0.125 * ( 1.0 - rPoint[0] ) * ( 1.0 - rPoint[2] );
        result( 0, 2 ) = -0.125 * ( 1.0 - rPoint[0] ) * ( 1.0 - rPoint[1] );
        result( 1, 0 ) =  0.125 * ( 1.0 - rPoint[1] ) * ( 1.0 - rPoint[2] );
        result( 1, 1 ) = -0.125 * ( 1.0 + rPoint[0] ) * ( 1.0 - rPoint[2] );
        result( 1, 2 ) = -0.125 * ( 1.0 + rPoint[0] ) * ( 1.0 - rPoint[1] );
        result( 2, 0 ) =  0.125 * ( 1.0 + rPoint[1] ) * ( 1.0 - rPoint[2] );
        result( 2, 1 ) =  0.125 * ( 1.0 + rPoint[0] ) * ( 1.0 - rPoint[2] );
        result( 2, 2 ) = -0.125 * ( 1.0 + rPoint[0] ) * ( 1.0 + rPoint[1] );
        result( 3, 0 ) = -0.125 * ( 1.0 + rPoint[1] ) * ( 1.0 - rPoint[2] );
        result( 3, 1 ) =  0.125 * ( 1.0 - rPoint[0] ) * ( 1.0 - rPoint[2] );
        result( 3, 2 ) = -0.125 * ( 1.0 - rPoint[0] ) * ( 1.0 + rPoint[1] );
        result( 4, 0 ) = -0.125 * ( 1.0 - rPoint[1] ) * ( 1.0 + rPoint[2] );
        result( 4, 1 ) = -0.125 * ( 1.0 - rPoint[0] ) * ( 1.0 + rPoint[2] );
        result( 4, 2 ) =  0.125 * ( 1.0 - rPoint[0] ) * ( 1.0 - rPoint[1] );
        result( 5, 0 ) =  0.125 * ( 1.0 - rPoint[1] ) * ( 1.0 + rPoint[2] );
        result( 5, 1 ) = -0.125 * ( 1.0 + rPoint[0] ) * ( 1.0 + rPoint[2] );
        result( 5, 2 ) =  0.125 * ( 1.0 + rPoint[0] ) * ( 1.0 - rPoint[1] );
        result( 6, 0 ) =  0.125 * ( 1.0 + rPoint[1] ) * ( 1.0 + rPoint[2] );
        result( 6, 1 ) =  0.125 * ( 1.0 + rPoint[0] ) * ( 1.0 + rPoint[2] );
        result( 6, 2 ) =  0.125 * ( 1.0 + rPoint[0] ) * ( 1.0 + rPoint[1] );
        result( 7, 0 ) = -0.125 * ( 1.0 + rPoint[1] ) * ( 1.0 + rPoint[2] );
        result( 7, 1 ) =  0.125 * ( 1.0 - rPoint[0] ) * ( 1.0 + rPoint[2] );
        result( 7, 2 ) =  0.125 * ( 1.0 - rPoint[0] ) * ( 1.0 + rPoint[1] );
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
        IntegrationPointsContainerType all_integration_points = AllIntegrationPoints();
        IntegrationPointsArrayType& integration_points = all_integration_points[static_cast<int>(ThisMethod)];

        //number of integration points
        const int integration_points_number = integration_points.size();
        //number of nodes in current geometry
        const int points_number = 8;
        //setting up return matrix
        Matrix shape_function_values( integration_points_number, points_number );
        //loop over all integration points

        for ( int pnt = 0; pnt < integration_points_number; pnt++ )
        {
            shape_function_values( pnt, 0 ) =
                0.125 * ( 1.0 - integration_points[pnt].X() )
                * ( 1.0 - integration_points[pnt].Y() )
                * ( 1.0 - integration_points[pnt].Z() );
            shape_function_values( pnt, 1 ) =
                0.125 * ( 1.0 + integration_points[pnt].X() )
                * ( 1.0 - integration_points[pnt].Y() )
                * ( 1.0 - integration_points[pnt].Z() );
            shape_function_values( pnt, 2 ) =
                0.125 * ( 1.0 + integration_points[pnt].X() )
                * ( 1.0 + integration_points[pnt].Y() )
                * ( 1.0 - integration_points[pnt].Z() );
            shape_function_values( pnt, 3 ) =
                0.125 * ( 1.0 - integration_points[pnt].X() )
                * ( 1.0 + integration_points[pnt].Y() )
                * ( 1.0 - integration_points[pnt].Z() );
            shape_function_values( pnt, 4 ) =
                0.125 * ( 1.0 - integration_points[pnt].X() )
                * ( 1.0 - integration_points[pnt].Y() )
                * ( 1.0 + integration_points[pnt].Z() );
            shape_function_values( pnt, 5 ) =
                0.125 * ( 1.0 + integration_points[pnt].X() )
                * ( 1.0 - integration_points[pnt].Y() )
                * ( 1.0 + integration_points[pnt].Z() );
            shape_function_values( pnt, 6 ) =
                0.125 * ( 1.0 + integration_points[pnt].X() )
                * ( 1.0 + integration_points[pnt].Y() )
                * ( 1.0 + integration_points[pnt].Z() );
            shape_function_values( pnt, 7 ) =
                0.125 * ( 1.0 - integration_points[pnt].X() )
                * ( 1.0 + integration_points[pnt].Y() )
                * ( 1.0 + integration_points[pnt].Z() );
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
        IntegrationPointsArrayType integration_points = all_integration_points[static_cast<int>(ThisMethod)];
        //number of integration points
        const int integration_points_number = integration_points.size();
        ShapeFunctionsGradientsType d_shape_f_values( integration_points_number );
        //initialising container
        //loop over all integration points

        for ( int pnt = 0; pnt < integration_points_number; pnt++ )
        {
            Matrix& result = d_shape_f_values[pnt];
            result = ZeroMatrix( 8, 3 );
            result( 0, 0 ) =
                -0.125 * ( 1.0 - integration_points[pnt].Y() )
                * ( 1.0 - integration_points[pnt].Z() );
            result( 0, 1 ) =
                -0.125 * ( 1.0 - integration_points[pnt].X() )
                * ( 1.0 - integration_points[pnt].Z() );
            result( 0, 2 ) =
                -0.125 * ( 1.0 - integration_points[pnt].X() )
                * ( 1.0 - integration_points[pnt].Y() );
            result( 1, 0 ) =
                0.125 * ( 1.0 - integration_points[pnt].Y() )
                * ( 1.0 - integration_points[pnt].Z() );
            result( 1, 1 ) =
                -0.125 * ( 1.0 + integration_points[pnt].X() )
                * ( 1.0 - integration_points[pnt].Z() );
            result( 1, 2 ) =
                -0.125 * ( 1.0 + integration_points[pnt].X() )
                * ( 1.0 - integration_points[pnt].Y() );
            result( 2, 0 ) =
                0.125 * ( 1.0 + integration_points[pnt].Y() )
                * ( 1.0 - integration_points[pnt].Z() );
            result( 2, 1 ) =
                0.125 * ( 1.0 + integration_points[pnt].X() )
                * ( 1.0 - integration_points[pnt].Z() );
            result( 2, 2 ) =
                -0.125 * ( 1.0 + integration_points[pnt].X() )
                * ( 1.0 + integration_points[pnt].Y() );
            result( 3, 0 ) =
                -0.125 * ( 1.0 + integration_points[pnt].Y() )
                * ( 1.0 - integration_points[pnt].Z() );
            result( 3, 1 ) =
                0.125 * ( 1.0 - integration_points[pnt].X() )
                * ( 1.0 - integration_points[pnt].Z() );
            result( 3, 2 ) =
                -0.125 * ( 1.0 - integration_points[pnt].X() )
                * ( 1.0 + integration_points[pnt].Y() );
            result( 4, 0 ) =
                -0.125 * ( 1.0 - integration_points[pnt].Y() )
                * ( 1.0 + integration_points[pnt].Z() );
            result( 4, 1 ) =
                -0.125 * ( 1.0 - integration_points[pnt].X() )
                * ( 1.0 + integration_points[pnt].Z() );
            result( 4, 2 ) =
                0.125 * ( 1.0 - integration_points[pnt].X() )
                * ( 1.0 - integration_points[pnt].Y() );
            result( 5, 0 ) =
                0.125 * ( 1.0 - integration_points[pnt].Y() )
                * ( 1.0 + integration_points[pnt].Z() );
            result( 5, 1 ) =
                -0.125 * ( 1.0 + integration_points[pnt].X() )
                * ( 1.0 + integration_points[pnt].Z() );
            result( 5, 2 ) =
                0.125 * ( 1.0 + integration_points[pnt].X() )
                * ( 1.0 - integration_points[pnt].Y() );
            result( 6, 0 ) =
                0.125 * ( 1.0 + integration_points[pnt].Y() )
                * ( 1.0 + integration_points[pnt].Z() );
            result( 6, 1 ) =
                0.125 * ( 1.0 + integration_points[pnt].X() )
                * ( 1.0 + integration_points[pnt].Z() );
            result( 6, 2 ) =
                0.125 * ( 1.0 + integration_points[pnt].X() )
                * ( 1.0 + integration_points[pnt].Y() );
            result( 7, 0 ) =
                -0.125 * ( 1.0 + integration_points[pnt].Y() )
                * ( 1.0 + integration_points[pnt].Z() );
            result( 7, 1 ) =
                0.125 * ( 1.0 - integration_points[pnt].X() )
                * ( 1.0 + integration_points[pnt].Z() );
            result( 7, 2 ) =
                0.125 * ( 1.0 - integration_points[pnt].X() )
                * ( 1.0 + integration_points[pnt].Y() );
        }

        return d_shape_f_values;
    }

    static const IntegrationPointsContainerType AllIntegrationPoints()
    {
        IntegrationPointsContainerType integration_points =
        {
            {
                Quadrature < HexahedronGaussLegendreIntegrationPoints1, 3, IntegrationPoint<3> >::GenerateIntegrationPoints(),
                Quadrature < HexahedronGaussLegendreIntegrationPoints2, 3, IntegrationPoint<3> >::GenerateIntegrationPoints(),
                Quadrature < HexahedronGaussLegendreIntegrationPoints3, 3, IntegrationPoint<3> >::GenerateIntegrationPoints(),
                Quadrature < HexahedronGaussLegendreIntegrationPoints4, 3, IntegrationPoint<3> >::GenerateIntegrationPoints(),
                Quadrature < HexahedronGaussLegendreIntegrationPoints5, 3, IntegrationPoint<3> >::GenerateIntegrationPoints(),
                Quadrature < HexahedronGaussLobattoIntegrationPoints1, 3, IntegrationPoint<3> >::GenerateIntegrationPoints(),
                Quadrature < HexahedronGaussLobattoIntegrationPoints2, 3, IntegrationPoint<3> >::GenerateIntegrationPoints(),
                Quadrature < HexahedronGaussLobattoIntegrationPoints2, 3, IntegrationPoint<3> >::GenerateIntegrationPoints()
            }
        };
        return integration_points;
    }

    static const ShapeFunctionsValuesContainerType AllShapeFunctionsValues()
    {
        ShapeFunctionsValuesContainerType shape_functions_values =
        {
            {
                Hexahedra3D8<TPointType>::CalculateShapeFunctionsIntegrationPointsValues( GeometryData::IntegrationMethod::GI_GAUSS_1 ),
                Hexahedra3D8<TPointType>::CalculateShapeFunctionsIntegrationPointsValues( GeometryData::IntegrationMethod::GI_GAUSS_2 ),
                Hexahedra3D8<TPointType>::CalculateShapeFunctionsIntegrationPointsValues( GeometryData::IntegrationMethod::GI_GAUSS_3 ),
                Hexahedra3D8<TPointType>::CalculateShapeFunctionsIntegrationPointsValues( GeometryData::IntegrationMethod::GI_GAUSS_4 ),
                Hexahedra3D8<TPointType>::CalculateShapeFunctionsIntegrationPointsValues( GeometryData::IntegrationMethod::GI_GAUSS_5 ),
                Hexahedra3D8<TPointType>::CalculateShapeFunctionsIntegrationPointsValues( GeometryData::IntegrationMethod::GI_EXTENDED_GAUSS_1 ),
                Hexahedra3D8<TPointType>::CalculateShapeFunctionsIntegrationPointsValues( GeometryData::IntegrationMethod::GI_EXTENDED_GAUSS_2 ),
                Hexahedra3D8<TPointType>::CalculateShapeFunctionsIntegrationPointsValues( GeometryData::IntegrationMethod::GI_LOBATTO_1 )
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
                Hexahedra3D8<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients( GeometryData::IntegrationMethod::GI_GAUSS_1 ),
                Hexahedra3D8<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients( GeometryData::IntegrationMethod::GI_GAUSS_2 ),
                Hexahedra3D8<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients( GeometryData::IntegrationMethod::GI_GAUSS_3 ),
                Hexahedra3D8<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients( GeometryData::IntegrationMethod::GI_GAUSS_4 ),
                Hexahedra3D8<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients( GeometryData::IntegrationMethod::GI_GAUSS_5 ),
                Hexahedra3D8<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients( GeometryData::IntegrationMethod::GI_EXTENDED_GAUSS_1 ),
                Hexahedra3D8<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients( GeometryData::IntegrationMethod::GI_EXTENDED_GAUSS_2 ),
                Hexahedra3D8<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients( GeometryData::IntegrationMethod::GI_LOBATTO_1 )
            }
        };
        return shape_functions_local_gradients;
    }


    /**
     * Private Friends
     */

    template<class TOtherPointType> friend class Hexahedra3D8;


    /**
     * Un accessible methods
     */

};// Class Hexahedra3D8


/**
 * Input and output
 */

/**
 * input stream function
 */
template<class TPointType> inline std::istream& operator >> (
    std::istream& rIStream, Hexahedra3D8<TPointType>& rThis );

/**
 * output stream function
 */
template<class TPointType> inline std::ostream& operator << (
    std::ostream& rOStream, const Hexahedra3D8<TPointType>& rThis )
{
    rThis.PrintInfo( rOStream );
    rOStream << std::endl;
    rThis.PrintData( rOStream );

    return rOStream;
}


template<class TPointType> const
GeometryData Hexahedra3D8<TPointType>::msGeometryData(
    &msGeometryDimension,
    GeometryData::IntegrationMethod::GI_GAUSS_2,
    Hexahedra3D8<TPointType>::AllIntegrationPoints(),
    Hexahedra3D8<TPointType>::AllShapeFunctionsValues(),
    AllShapeFunctionsLocalGradients()
);

template<class TPointType> const
GeometryDimension Hexahedra3D8<TPointType>::msGeometryDimension(3, 3);

}// namespace Kratos.
