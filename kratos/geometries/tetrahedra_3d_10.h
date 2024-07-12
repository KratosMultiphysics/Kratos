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
//  Contributors:    Hoang Giang Bui
//                   Josep Maria Carbonell
//

#pragma once

// System includes
#include <numeric>

// External includes

// Project includes
#include "geometries/triangle_3d_6.h"
#include "geometries/tetrahedra_3d_4.h"
#include "utilities/geometry_utilities.h"
#include "utilities/integration_utilities.h"
#include "integration/tetrahedron_gauss_legendre_integration_points.h"

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
 * @class Tetrahedra3D10
 * @ingroup KratosCore
 * @brief A ten node tetrahedra geometry with quadratic shape functions
 * @details The node ordering corresponds with:
 *                     3
 *                   ,/|`\
 *                 ,/  |  `\
 *               ,7    '.   `9
 *             ,/       8     `\
 *          ,/          |       `\
 *         0--------6--'.--------2
 *          `\.         |      ,/
 *             `\.      |    ,5
 *                `4.   '. ,/
 *                   `\. |/
 *                      `1
 * @author Riccardo Rossi
 * @author Janosch Stascheit
 * @author Felix Nagel
 */
template<class TPointType>
class Tetrahedra3D10
    : public Geometry<TPointType>
{
public:
    ///@name Type Definitions
    ///@{

    /// Geometry as base class.
    using BaseType = Geometry<TPointType>;

    /// Type of edge and face geometries
    using EdgeType = Line3D3<TPointType>;
    using FaceType = Triangle3D6<TPointType>;

    /// Pointer definition of Tetrahedra3D10
    KRATOS_CLASS_POINTER_DEFINITION( Tetrahedra3D10 );

    /// Integration methods implemented in geometry.
    using IntegrationMethod = GeometryData::IntegrationMethod;

    /// A Vector of counted pointers to Geometries. Used for
    /// returning edges of the geometry.
    using GeometriesArrayType = typename BaseType::GeometriesArrayType;

    /// Redefinition of template parameter TPointType.
    using PointType = TPointType;

    /// Type used for indexing in geometry class. std::size_t used for indexing
    /// point or integration point access methods and also all other
    /// methods which need point or integration point index.
    using IndexType = typename BaseType::IndexType;

    /// This type used to return size or dimension in
    /// geometry. Dimension, WorkingDimension, PointsNumber, and
    /// ... return this type as their results.
    using SizeType = typename BaseType::SizeType;

    /// Array of counted pointers to point. This type used to hold
    /// geometry's points.
    using PointsArrayType = typename BaseType::PointsArrayType;

    /// This type used for representing an integration point in
    /// geometry. This integration point is a point with an
    /// additional weight component.
    using IntegrationPointType = typename BaseType::IntegrationPointType;

    /// A Vector of IntegrationPointType which used to hold
    /// integration points related to an integration
    /// method. IntegrationPoints functions used this type to return
    /// their results.
    using IntegrationPointsArrayType = typename BaseType::IntegrationPointsArrayType;

    /// A Vector of IntegrationPointsArrayType which used to hold
    /// integration points related to different integration method
    /// implemented in geometry.
    using IntegrationPointsContainerType = typename BaseType::IntegrationPointsContainerType;

    /// A third order tensor used as shape functions' values
    /// container.
    using ShapeFunctionsValuesContainerType = typename BaseType::ShapeFunctionsValuesContainerType;

    /// A fourth order tensor used as shape functions' local
    /// gradients container in geometry.
    using ShapeFunctionsLocalGradientsContainerType = typename BaseType::ShapeFunctionsLocalGradientsContainerType;

    /// A third order tensor to hold Jacobian matrices evaluated at
    /// integration points. Jacobian and InverseOfJacobian functions
    /// return this type as their result.
    using JacobiansType = typename BaseType::JacobiansType;

    /// A third order tensor to hold shape functions' local
    /// gradients. ShapefunctionsLocalGradients function return this
    /// type as its result.
    using ShapeFunctionsGradientsType = typename BaseType::ShapeFunctionsGradientsType;

    /// A third order tensor to hold shape functions' local second derivatives.
    /// ShapeFunctionsSecondDerivatives function return this type as its result.
    using ShapeFunctionsSecondDerivativesType = typename BaseType::ShapeFunctionsSecondDerivativesType;

    /// Type of the normal vector used for normal to edges in geometry.
    using NormalType = typename BaseType::NormalType;

    /// Type of coordinates array
    using CoordinatesArrayType = typename BaseType::CoordinatesArrayType;

    /// Type of Matrix
    using MatrixType = Matrix;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
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

    /**
     * @brief Constructor with points array
     * @param ThisPoints The points of the current geometry
     */
    Tetrahedra3D10( const PointsArrayType& ThisPoints )
        : BaseType( ThisPoints, &msGeometryData )
    {
        KRATOS_ERROR_IF( this->PointsNumber() != 10 ) << "Invalid points number. Expected 10, given " << this->PointsNumber() << std::endl;
    }

    /// Constructor with Geometry Id
    explicit Tetrahedra3D10(
        const IndexType GeometryId,
        const PointsArrayType& rThisPoints
    ) : BaseType(GeometryId, rThisPoints, &msGeometryData)
    {
        KRATOS_ERROR_IF( this->PointsNumber() != 10 ) << "Invalid points number. Expected 10, given " << this->PointsNumber() << std::endl;
    }

    /// Constructor with Geometry Name
    explicit Tetrahedra3D10(
        const std::string& rGeometryName,
        const PointsArrayType& rThisPoints
    ) : BaseType(rGeometryName, rThisPoints, &msGeometryData)
    {
        KRATOS_ERROR_IF(this->PointsNumber() != 10) << "Invalid points number. Expected 10, given " << this->PointsNumber() << std::endl;
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
    ~Tetrahedra3D10() override {}

    GeometryData::KratosGeometryFamily GetGeometryFamily() const override
    {
        return GeometryData::KratosGeometryFamily::Kratos_Tetrahedra;
    }

    GeometryData::KratosGeometryType GetGeometryType() const override
    {
        return GeometryData::KratosGeometryType::Kratos_Tetrahedra3D10;
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
        return typename BaseType::Pointer( new Tetrahedra3D10( rThisPoints ) );
    }

    /**
     * @brief Creates a new geometry pointer
     * @param NewGeometryId the ID of the new geometry
     * @param ThisPoints the nodes of the new geometry
     * @return Pointer to the new geometry
     */
    typename BaseType::Pointer Create(
        const IndexType NewGeometryId,
        PointsArrayType const& rThisPoints
        ) const override
    {
        return typename BaseType::Pointer( new Tetrahedra3D10( NewGeometryId, rThisPoints ) );
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
        auto p_geometry = typename BaseType::Pointer( new Tetrahedra3D10( rGeometry.Points() ) );
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
        auto p_geometry = typename BaseType::Pointer( new Tetrahedra3D10( NewGeometryId, rGeometry.Points() ) );
        p_geometry->SetData(rGeometry.GetData());
        return p_geometry;
    }

    ///@}
    ///@name Information
    ///@{

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
     */
    double DomainSize() const override
    {
        return Volume();
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
        // Using linear approximation for planar faces
        if (this->FacesArePlanar()) {
            return GeometryUtils::PointLocalCoordinatesPlanarFaceTetrahedra(*this, rResult, rPoint);
        } else {
            return BaseType::PointLocalCoordinates( rResult, rPoint );
        }
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

        if ( (rResult[0] >= (0.0 - Tolerance)) && (rResult[0] <= (1.0 + Tolerance)) ) {
            if ( (rResult[1] >= (0.0 - Tolerance)) && (rResult[1] <= (1.0 + Tolerance)) ) {
                if ( (rResult[2] >= (0.0 - Tolerance)) && (rResult[2] <= (1.0 + Tolerance)) ) {
                    if ((( 1.0 - ( rResult[0] + rResult[1] + rResult[2] ) ) >= (0.0 - Tolerance) ) && (( 1.0 - ( rResult[0] + rResult[1] + rResult[2] ) ) <= (1.0 + Tolerance) ) ) {
                        return true;
                    }
                }
            }
        }

        return false;
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
        return 6;
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
                                              this->pGetPoint( 1 ),
                                              this->pGetPoint( 4 ) ) ) );
        edges.push_back( EdgePointerType( new EdgeType(
                                              this->pGetPoint( 1 ),
                                              this->pGetPoint( 2 ),
                                              this->pGetPoint( 5 ) ) ) );
        edges.push_back( EdgePointerType( new EdgeType(
                                              this->pGetPoint( 2 ),
                                              this->pGetPoint( 0 ),
                                              this->pGetPoint( 6 ) ) ) );
        edges.push_back( EdgePointerType( new EdgeType(
                                              this->pGetPoint( 0 ),
                                              this->pGetPoint( 3 ),
                                              this->pGetPoint( 7 ) ) ) );
        edges.push_back( EdgePointerType( new EdgeType(
                                              this->pGetPoint( 1 ),
                                              this->pGetPoint( 3 ),
                                              this->pGetPoint( 8 ) ) ) );
        edges.push_back( EdgePointerType( new EdgeType(
                                              this->pGetPoint( 2 ),
                                              this->pGetPoint( 3 ),
                                              this->pGetPoint( 9 ) ) ) );

        return edges;
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
        return 4;
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

    /** This method calculates and returns the average edge length of the geometry
     *
     * @return double value with the average edge length
     *
     */
    double AverageEdgeLength() const override {
        const GeometriesArrayType edges = this->GenerateEdges();
        return std::accumulate(
            edges.begin(),
            edges.end(),
            0.0,
            [](double sum, const auto& rEdge) {return sum + rEdge.Length();}
        ) * 0.16666666666666666667;
    }

    Matrix& PointsLocalCoordinates( Matrix& rResult ) const override
    {
        if(rResult.size1()!= 10 || rResult.size2()!= 3)
            rResult.resize(10, 3, false);
        rResult(0,0)=0.0;
        rResult(0,1)=0.0;
        rResult(0,2)=0.0;
        rResult(1,0)=1.0;
        rResult(1,1)=0.0;
        rResult(1,2)=0.0;
        rResult(2,0)=0.0;
        rResult(2,1)=1.0;
        rResult(2,2)=0.0;
        rResult(3,0)=0.0;
        rResult(3,1)=0.0;
        rResult(3,2)=1.0;
        rResult(4,0)=0.5;
        rResult(4,1)=0.0;
        rResult(4,2)=0.0;
        rResult(5,0)=0.5;
        rResult(5,1)=0.5;
        rResult(5,2)=0.0;
        rResult(6,0)=0.0;
        rResult(6,1)=0.5;
        rResult(6,2)=0.0;
        rResult(7,0)=0.0;
        rResult(7,1)=0.0;
        rResult(7,2)=0.5;
        rResult(8,0)=0.5;
        rResult(8,1)=0.0;
        rResult(8,2)=0.5;
        rResult(9,0)=0.0;
        rResult(9,1)=0.5;
        rResult(9,2)=0.5;
        return rResult;
    }

    ///@}
    ///@name Shape Function
    ///@{

    Vector& ShapeFunctionsValues(Vector &rResult, const CoordinatesArrayType& rCoordinates) const override
    {
        ShapeFunctionsValuesImpl(rResult, rCoordinates);
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
            KRATOS_ERROR << "Wrong index of shape function!" << *this << std::endl;
        }

        return 0;
    }

    Matrix& ShapeFunctionsLocalGradients(Matrix& rResult,
            const CoordinatesArrayType& rPoint) const override
    {
        const double fourthCoord = 1.0 - (rPoint[0] + rPoint[1] + rPoint[2]);
        if (rResult.size1() != this->size() || rResult.size2() != this->LocalSpaceDimension())
            rResult.resize(this->size(), this->LocalSpaceDimension(), false);

        rResult(0, 0) = -(4.0 * fourthCoord - 1.0);
        rResult(0, 1) = -(4.0 * fourthCoord - 1.0);
        rResult(0, 2) = -(4.0 * fourthCoord - 1.0);
        rResult(1, 0) =  4.0 * rPoint[0] - 1.0;
        rResult(1, 1) =  0.0;
        rResult(1, 2) =  0.0;
        rResult(2, 0) =  0.0;
        rResult(2, 1) =  4.0 * rPoint[1] - 1.0;
        rResult(2, 2) =  0.0;
        rResult(3, 0) =  0.0;
        rResult(3, 1) =  0.0;
        rResult(3, 2) =  4.0 * rPoint[2] - 1.0;
        rResult(4, 0) = -4.0 * rPoint[0] + 4.0 * fourthCoord;
        rResult(4, 1) = -4.0 * rPoint[0];
        rResult(4, 2) = -4.0 * rPoint[0];
        rResult(5, 0) =  4.0 * rPoint[1];
        rResult(5, 1) =  4.0 * rPoint[0];
        rResult(5, 2) =  0.0;
        rResult(6, 0) = -4.0 * rPoint[1];
        rResult(6, 1) = -4.0 * rPoint[1] + 4.0 * fourthCoord;
        rResult(6, 2) = -4.0 * rPoint[1];
        rResult(7, 0) = -4.0 * rPoint[2];
        rResult(7, 1) = -4.0 * rPoint[2];
        rResult(7, 2) = -4.0 * rPoint[2] + 4.0 * fourthCoord;
        rResult(8, 0) =  4.0 * rPoint[2];
        rResult(8, 1) =  0.0;
        rResult(8, 2) =  4.0 * rPoint[0];
        rResult(9, 0) =  0.0;
        rResult(9, 1) =  4.0 * rPoint[2];
        rResult(9, 2) =  4.0 * rPoint[1];

        return rResult;
    }

    ShapeFunctionsSecondDerivativesType& ShapeFunctionsSecondDerivatives(
        ShapeFunctionsSecondDerivativesType& rResult,
        const CoordinatesArrayType& rPoint) const override
    {
        // Check and resize results container
        if (rResult.size() != this->PointsNumber()) {
            rResult.resize(this->PointsNumber(), false);
        }

        for (IndexType i = 0; i < this->PointsNumber(); ++i) {
            auto& r_DDN_i = rResult[i];
            if (r_DDN_i.size1() != 3 || r_DDN_i.size2() != 3) {
                r_DDN_i.resize(3,3, false);
            }
        }

        // Node 0
        rResult[0](0,0) = 4.0; rResult[0](0,1) = 4.0; rResult[0](0,2) = 4.0;
        rResult[0](1,0) = 4.0; rResult[0](1,1) = 4.0; rResult[0](1,2) = 4.0;
        rResult[0](2,0) = 4.0; rResult[0](2,1) = 4.0; rResult[0](2,2) = 4.0;

        // Node 1
        rResult[1](0,0) = 4.0; rResult[1](0,1) = 0.0; rResult[1](0,2) = 0.0;
        rResult[1](1,0) = 0.0; rResult[1](1,1) = 0.0; rResult[1](1,2) = 0.0;
        rResult[1](2,0) = 0.0; rResult[1](2,1) = 0.0; rResult[1](2,2) = 0.0;

        // Node 2
        rResult[2](0,0) = 0.0; rResult[2](0,1) = 0.0; rResult[2](0,2) = 0.0;
        rResult[2](1,0) = 0.0; rResult[2](1,1) = 4.0; rResult[2](1,2) = 0.0;
        rResult[2](2,0) = 0.0; rResult[2](2,1) = 0.0; rResult[2](2,2) = 0.0;

        // Node 3
        rResult[3](0,0) = 0.0; rResult[3](0,1) = 0.0; rResult[3](0,2) = 0.0;
        rResult[3](1,0) = 0.0; rResult[3](1,1) = 0.0; rResult[3](1,2) = 0.0;
        rResult[3](2,0) = 0.0; rResult[3](2,1) = 0.0; rResult[3](2,2) = 4.0;

        // Node 4
        rResult[4](0,0) = -8.0; rResult[4](0,1) = -4.0; rResult[4](0,2) = -4.0;
        rResult[4](1,0) = -4.0; rResult[4](1,1) = 0.0; rResult[4](1,2) = 0.0;
        rResult[4](2,0) = -4.0; rResult[4](2,1) = 0.0; rResult[4](2,2) = 0.0;

        // Node 5
        rResult[5](0,0) = 0.0; rResult[5](0,1) = 4.0; rResult[5](0,2) = 0.0;
        rResult[5](1,0) = 4.0; rResult[5](1,1) = 0.0; rResult[5](1,2) = 0.0;
        rResult[5](2,0) = 0.0; rResult[5](2,1) = 0.0; rResult[5](2,2) = 0.0;

        // Node 6
        rResult[6](0,0) = 0.0; rResult[6](0,1) = -4.0; rResult[6](0,2) = 0.0;
        rResult[6](1,0) = -4.0; rResult[6](1,1) = -8.0; rResult[6](1,2) = -4.0;
        rResult[6](2,0) = 0.0; rResult[6](2,1) = -4.0; rResult[6](2,2) = 0.0;

        // Node 7
        rResult[7](0,0) = 0.0; rResult[7](0,1) = 0.0; rResult[7](0,2) = -4.0;
        rResult[7](1,0) = 0.0; rResult[7](1,1) = 0.0; rResult[7](1,2) = -4.0;
        rResult[7](2,0) = -4.0; rResult[7](2,1) = -4.0; rResult[7](2,2) = -8.0;

        // Node 8
        rResult[8](0,0) = 0.0; rResult[8](0,1) = 0.0; rResult[8](0,2) = 4.0;
        rResult[8](1,0) = 0.0; rResult[8](1,1) = 0.0; rResult[8](1,2) = 0.0;
        rResult[8](2,0) = 4.0; rResult[8](2,1) = 0.0; rResult[8](2,2) = 0.0;

        // Node 9
        rResult[9](0,0) = 0.0; rResult[9](0,1) = 0.0; rResult[9](0,2) = 0.0;
        rResult[9](1,0) = 0.0; rResult[9](1,1) = 0.0; rResult[9](1,2) = 4.0;
        rResult[9](2,0) = 0.0; rResult[9](2,1) = 4.0; rResult[9](2,2) = 0.0;

        return rResult;
    }

    /** Tests the intersection of the geometry with
     * a 3D box defined by rLowPoint and rHighPoint.
     * The method is only implemented for simple tets
     * where the faces are planar.
     *
     * @param  rLowPoint  Lower point of the box to test the intersection
     * @param  rHighPoint Higher point of the box to test the intersection
     * @return            True if the geometry intersects the box, False in any other case.
     */
    bool HasIntersection(const Point& rLowPoint, const Point& rHighPoint) const override
    {
        // Check if the faces are planar
        if (this->FacesArePlanar()) {
            // TODO: Move implementation to GeometryUtils to avoid creating a new tetrahedra
            return Tetrahedra3D4<TPointType>(
                this->pGetPoint(0),
                this->pGetPoint(1),
                this->pGetPoint(2),
                this->pGetPoint(3)).HasIntersection(rLowPoint, rHighPoint);
        } else {
             KRATOS_ERROR << "\"HasIntersection\" is not implemented for non-planar 10 noded tetrahedra.";
        }
        return false;
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
        std::array<double, 4> distances;
        distances[0] = GeometryUtils::PointDistanceToTriangle3D(this->GetPoint(0), this->GetPoint(2), this->GetPoint(1), this->GetPoint(6), this->GetPoint(5), this->GetPoint(4), point);
        distances[1] = GeometryUtils::PointDistanceToTriangle3D(this->GetPoint(0), this->GetPoint(3), this->GetPoint(2), this->GetPoint(7), this->GetPoint(9), this->GetPoint(6), point);
        distances[2] = GeometryUtils::PointDistanceToTriangle3D(this->GetPoint(0), this->GetPoint(1), this->GetPoint(3), this->GetPoint(4), this->GetPoint(8), this->GetPoint(7), point);
        distances[3] = GeometryUtils::PointDistanceToTriangle3D(this->GetPoint(2), this->GetPoint(3), this->GetPoint(1), this->GetPoint(9), this->GetPoint(8), this->GetPoint(5), point);
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
        return "3 dimensional tetrahedra with ten nodes in 3D space";
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

    ///@}
    ///@name Private Life Cycle
    ///@{

    Tetrahedra3D10(): BaseType( PointsArrayType(), &msGeometryData ) {}

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief Returns vector of shape function values at local coordinate.
     *
     * For a definition of the shape functions see, e.g.,
     * P. Wriggers, Nonlinear Finite Element Methods, Springer, 2008, Sec. 4.1.
     */
    static void ShapeFunctionsValuesImpl(Vector &rResult, const CoordinatesArrayType& rCoordinates)
    {
        if (rResult.size() != 10)
            rResult.resize(10, false);
        const double fourthCoord = 1.0 - rCoordinates[0] - rCoordinates[1] - rCoordinates[2];
        rResult[0] = fourthCoord * (2.0 * fourthCoord - 1.0);
        rResult[1] = rCoordinates[0] * (2.0 * rCoordinates[0] - 1.0);
        rResult[2] = rCoordinates[1] * (2.0 * rCoordinates[1] - 1.0);
        rResult[3] = rCoordinates[2] * (2.0 * rCoordinates[2] - 1.0);
        rResult[4] = 4.0 * fourthCoord * rCoordinates[0];
        rResult[5] = 4.0 * rCoordinates[0] * rCoordinates[1];
        rResult[6] = 4.0 * rCoordinates[1] * fourthCoord;
        rResult[7] = 4.0 * rCoordinates[2] * fourthCoord;
        rResult[8] = 4.0 * rCoordinates[0] * rCoordinates[2];
        rResult[9] = 4.0 * rCoordinates[1] * rCoordinates[2];
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
        const std::size_t points_number = 10;
        IntegrationPointsContainerType all_integration_points = AllIntegrationPoints();
        IntegrationPointsArrayType integration_points = all_integration_points[static_cast<int>(ThisMethod)];
        Matrix shape_function_values(integration_points.size(), points_number);
        //loop over all integration points
        Vector N(points_number);
        for (std::size_t pnt = 0; pnt < integration_points.size(); ++pnt)
        {
            ShapeFunctionsValuesImpl(N, integration_points[pnt]);
            for (std::size_t i = 0; i < N.size(); ++i)
                shape_function_values(pnt, i) = N[i];
        }
        return shape_function_values;
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
                    GeometryData::IntegrationMethod::GI_GAUSS_1 ),
                Tetrahedra3D10<TPointType>::CalculateShapeFunctionsIntegrationPointsValues(
                    GeometryData::IntegrationMethod::GI_GAUSS_2 ),
                Tetrahedra3D10<TPointType>::CalculateShapeFunctionsIntegrationPointsValues(
                    GeometryData::IntegrationMethod::GI_GAUSS_3 ),
                Tetrahedra3D10<TPointType>::CalculateShapeFunctionsIntegrationPointsValues(
                    GeometryData::IntegrationMethod::GI_GAUSS_4 ),
                Tetrahedra3D10<TPointType>::CalculateShapeFunctionsIntegrationPointsValues(
                    GeometryData::IntegrationMethod::GI_GAUSS_5 )
            }
        };
        return shape_functions_values;
    }

    static const ShapeFunctionsLocalGradientsContainerType
    AllShapeFunctionsLocalGradients()
    {
        ShapeFunctionsLocalGradientsContainerType shape_functions_local_gradients =
        {
            {
                Tetrahedra3D10<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients(
                    GeometryData::IntegrationMethod::GI_GAUSS_1 ),
                Tetrahedra3D10<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients(
                    GeometryData::IntegrationMethod::GI_GAUSS_2 ),
                Tetrahedra3D10<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients(
                    GeometryData::IntegrationMethod::GI_GAUSS_3 ),
                Tetrahedra3D10<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients(
                    GeometryData::IntegrationMethod::GI_GAUSS_4 ),
                Tetrahedra3D10<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients(
                    GeometryData::IntegrationMethod::GI_GAUSS_5 )
            }
        };
        return shape_functions_local_gradients;
    }

    /**
     * @brief Checks if faces are planar. We iterate for all edges and check
     * that the sum of 0-2 and 2-1 segments is no bigger than 0-1.
     * @return bool faces are planar or not
     */
    bool FacesArePlanar() const
    {
        constexpr double tol = 1e-6;
        constexpr std::array<std::array<size_t, 3>, 6> edges{
            {{0, 1, 4}, {1, 2, 5}, {2, 0, 6}, {0, 3, 7}, {1, 3, 8}, {2, 3, 9}}};
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
    ///@name Private  Friends
    ///@{

    template<class TOtherPointType> friend class Tetrahedra3D10;

    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}
};// Class Tetrahedra3D10

///@}
///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
template<class TPointType> inline std::istream& operator >> (
    std::istream& rIStream, Tetrahedra3D10<TPointType>& rThis );

/// output stream function
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
    &msGeometryDimension,
    GeometryData::IntegrationMethod::GI_GAUSS_2,
    Tetrahedra3D10<TPointType>::AllIntegrationPoints(),
    Tetrahedra3D10<TPointType>::AllShapeFunctionsValues(),
    AllShapeFunctionsLocalGradients()
);

template<class TPointType> const
GeometryDimension Tetrahedra3D10<TPointType>::msGeometryDimension(3, 3);

///@}

}// namespace Kratos.
