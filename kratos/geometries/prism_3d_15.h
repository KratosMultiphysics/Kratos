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
//                   Vicente Mataix Ferrandiz
//

#pragma once

// System includes

// External includes

// Project includes
#include "geometries/triangle_3d_6.h"
#include "geometries/quadrilateral_3d_8.h"
#include "utilities/integration_utilities.h"
#include "integration/prism_gauss_legendre_integration_points.h"

namespace Kratos
{
/**
 * @class Prism3D15
 * @ingroup KratosCore
 * @brief A fifteen node prism geometry with quadratic shape functions
 * @details The node ordering corresponds with:
 *                 3
 *               ,/|`\
 *             12  |  14
 *           ,/    |    `\
 *          4------13-----5
 *          |      9      |
 *          |      |      |
 *          |      |      |
 *          |      |      |
 *          10     |      11
 *          |      0      |
 *          |    ,/ `\    |
 *          |  ,6     `8  |
 *          |,/         `\|
 *          1------7------2
 * @author Riccardo Rossi
 * @author Janosch Stascheit
 * @author Felix Nagel
 */
template<class TPointType>
class Prism3D15
    : public Geometry<TPointType>
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
    typedef Triangle3D6<TPointType> FaceType1;
    typedef Quadrilateral3D8<TPointType> FaceType2;

    /**
     * Pointer definition of Prism3D15
     */
    KRATOS_CLASS_POINTER_DEFINITION( Prism3D15 );

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
     * Type of the normal vector used for normal to edges in geometry.
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
    Prism3D15( typename PointType::Pointer pPoint1,
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
               typename PointType::Pointer pPoint13,
               typename PointType::Pointer pPoint14,
               typename PointType::Pointer pPoint15 )
        : BaseType( PointsArrayType(), &msGeometryData )
    {
        this->Points().reserve( 15 );
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
        this->Points().push_back( pPoint11 );
        this->Points().push_back( pPoint12 );
        this->Points().push_back( pPoint13 );
        this->Points().push_back( pPoint14 );
        this->Points().push_back( pPoint15 );
    }

    Prism3D15( const PointsArrayType& rThisPoints )
        : BaseType( rThisPoints, &msGeometryData )
    {
        KRATOS_ERROR_IF( this->PointsNumber() != 15 ) << "Invalid points number. Expected 15, given " << this->PointsNumber() << std::endl;
    }

    /// Constructor with Geometry Id
    explicit Prism3D15(
        const IndexType GeometryId,
        const PointsArrayType& rThisPoints
        ) : BaseType(GeometryId, rThisPoints, &msGeometryData)
    {
        KRATOS_ERROR_IF( this->PointsNumber() != 15 ) << "Invalid points number. Expected 15, given " << this->PointsNumber() << std::endl;
    }

    /// Constructor with Geometry Name
    explicit Prism3D15(
        const std::string& rGeometryName,
        const PointsArrayType& rThisPoints
        ) : BaseType(rGeometryName, rThisPoints, &msGeometryData)
    {
        KRATOS_ERROR_IF(this->PointsNumber() != 15) << "Invalid points number. Expected 15, given " << this->PointsNumber() << std::endl;
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
    Prism3D15( Prism3D15 const& rOther )
        : BaseType( rOther )
    {
    }

    /**
     * Copy constructor from a geometry with other point type.
     * Construct this geometry as a copy of given geometry which
     * has different type of points. The given goemetry's
     * TOtherPointType* must be implicitly convertible to this
     * geometry PointType.
     *
     * @note This copy constructor don't copy the points and new
     * geometry shares points with given source geometry. It's
     * obvious that any change to this new geometry's point affect
     * source geometry's points too.
     */
    template<class TOtherPointType> Prism3D15( Prism3D15<TOtherPointType> const& rOther )
        : BaseType( rOther )
    {
    }

    /// Destructor. Does nothing!!!
    ~Prism3D15() override {}

    /**
     * @brief Gets the geometry family.
     * @details This function returns the family type of the geometry. The geometry family categorizes the geometry into a broader classification, aiding in its identification and processing.
     * @return GeometryData::KratosGeometryFamily The geometry family.
     */
    GeometryData::KratosGeometryFamily GetGeometryFamily() const override
    {
        return GeometryData::KratosGeometryFamily::Kratos_Prism;
    }

    /**
     * @brief Gets the geometry type.
     * @details This function returns the specific type of the geometry. The geometry type provides a more detailed classification of the geometry.
     * @return GeometryData::KratosGeometryType The specific geometry type.
     */
    GeometryData::KratosGeometryType GetGeometryType() const override
    {
        return GeometryData::KratosGeometryType::Kratos_Prism3D15;
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
    Prism3D15& operator=( const Prism3D15& rOther )
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
    Prism3D15& operator=( Prism3D15<TOtherPointType> const & rOther )
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
        return typename BaseType::Pointer( new Prism3D15( NewGeometryId, rThisPoints ) );
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
        auto p_geometry = typename BaseType::Pointer( new Prism3D15( NewGeometryId, rGeometry.Points() ) );
        p_geometry->SetData(rGeometry.GetData());
        return p_geometry;
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
     */
    double DomainSize() const override
    {
        return this->Volume();
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

    double Area() const override
    {
        return std::abs( this->DeterminantOfJacobian( PointType() ) ) * 0.5;
    }

    double Length() const override
    {
        const double volume = Volume();
        return std::pow(volume, 1.0/3.0)/3.0;
    }


    /**
    * Returns a matrix of the local coordinates of all points
    * @param rResult a Matrix that will be overwritten by the results
    * @return the coordinates of all points of the current geometry
    */
    Matrix& PointsLocalCoordinates( Matrix& rResult ) const override
    {
        if ( rResult.size1() != 15 || rResult.size2() != 3 )
            rResult.resize( 15, 3, false );

        rResult( 0, 0 ) = 0.0; //corners of bottom base
        rResult( 0, 1 ) = 0.0;
        rResult( 0, 2 ) = 0.0;

        rResult( 1, 0 ) = 1.0;
        rResult( 1, 1 ) = 0.0;
        rResult( 1, 2 ) = 0.0;

        rResult( 2, 0 ) = 0.0;
        rResult( 2, 1 ) = 1.0;
        rResult( 2, 2 ) = 0.0;

        rResult( 3, 0 ) = 0.0; //corners of top base
        rResult( 3, 1 ) = 0.0;
        rResult( 3, 2 ) = 1.0;

        rResult( 4, 0 ) = 1.0;
        rResult( 4, 1 ) = 0.0;
        rResult( 4, 2 ) = 1.0;

        rResult( 5, 0 ) = 0.0;
        rResult( 5, 1 ) = 1.0;
        rResult( 5, 2 ) = 1.0;

        rResult( 6, 0 ) = 0.5; //mid point of bottom base
        rResult( 6, 1 ) = 0.0;
        rResult( 6, 2 ) = -1.0;

        rResult( 7, 0 ) = 0.5;
        rResult( 7, 1 ) = 0.5;
        rResult( 7, 2 ) = -1.0;

        rResult( 8, 0 ) = 0.0;
        rResult( 8, 1 ) = 0.5;
        rResult( 8, 2 ) = -1.0;

        rResult( 9, 0 ) = 0.0; //veritices of mid height
        rResult( 9, 1 ) = 0.0;
        rResult( 9, 2 ) = 0.5;

        rResult( 10, 0 ) = 1.0;
        rResult( 10, 1 ) = 0.0;
        rResult( 10, 2 ) = 0.5;

        rResult( 11, 0 ) = 0.0;
        rResult( 11, 1 ) = 1.0;
        rResult( 11, 2 ) = 0.5;

        rResult( 12, 0 ) = 0.5; //mid point of top base
        rResult( 12, 1 ) = 0.0;
        rResult( 12, 2 ) = 1.0;

        rResult( 13, 0 ) = 0.5;
        rResult( 13, 1 ) = 0.5;
        rResult( 13, 2 ) = 1.0;

        rResult( 14, 0 ) = 0.0;
        rResult( 14, 1 ) = 0.5;
        rResult( 14, 2 ) = 1.0;

        return rResult;
    }

    /**
     * Returns whether given arbitrary point is inside the Geometry and the respective
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

        if ( (rResult[0] >= (0.0 - Tolerance)) && (rResult[0] <= (1.0 + Tolerance)) )
            if ( (rResult[1] >= (0.0 - Tolerance)) && (rResult[1] <= (1.0 + Tolerance)) )
                if ( (rResult[2] >= (-1.0 - Tolerance)) && (rResult[2] <= (1.0 + Tolerance)) )
                    if ((( 1.0 - ( rResult[0] + rResult[1] ) ) >= (0.0 - Tolerance) ) && (( 1.0 - ( rResult[0] + rResult[1] ) ) <= (1.0 + Tolerance) ) )
                        return true;

        return false;
    }

    ///@}
    ///@name Edge
    ///@{

    /**
     * @brief This method gives you number of all edges of this geometry.
     * @details For example, for a hexahedron, this would be 12
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
        return 9;
    }

    /**
     * @brief This method gives you all edges of this geometry.
     * @details This method will gives you all the edges with one dimension less than this geometry.
     * For example a triangle would return three lines as its edges or a tetrahedral would return four triangle as its edges but won't return its six edge lines by this method.
     * @return GeometriesArrayType contains this geometry edges.
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
                                              this->pGetPoint( 6 ) ) ) );
        edges.push_back( EdgePointerType( new EdgeType(
                                              this->pGetPoint( 1 ),
                                              this->pGetPoint( 2 ),
                                              this->pGetPoint( 7 ) ) ) );
        edges.push_back( EdgePointerType( new EdgeType(
                                              this->pGetPoint( 2 ),
                                              this->pGetPoint( 0 ),
                                              this->pGetPoint( 8 ) ) ) );
        edges.push_back( EdgePointerType( new EdgeType(
                                              this->pGetPoint( 3 ),
                                              this->pGetPoint( 4 ),
                                              this->pGetPoint( 12 ) ) ) );
        edges.push_back( EdgePointerType( new EdgeType(
                                              this->pGetPoint( 4 ),
                                              this->pGetPoint( 5 ),
                                              this->pGetPoint( 13 ) ) ) );
        edges.push_back( EdgePointerType( new EdgeType(
                                              this->pGetPoint( 5 ),
                                              this->pGetPoint( 3 ),
                                              this->pGetPoint( 14 ) ) ) );
        edges.push_back( EdgePointerType( new EdgeType(
                                              this->pGetPoint( 0 ),
                                              this->pGetPoint( 3 ),
                                              this->pGetPoint( 9 ) ) ) );
        edges.push_back( EdgePointerType( new EdgeType(
                                              this->pGetPoint( 1 ),
                                              this->pGetPoint( 4 ),
                                              this->pGetPoint( 10 ) ) ) );
        edges.push_back( EdgePointerType( new EdgeType(
                                              this->pGetPoint( 2 ),
                                              this->pGetPoint( 5 ),
                                              this->pGetPoint( 11 ) ) ) );
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
        return 5;
    }

    /**
     * @brief Returns all faces of the current geometry.
     * @details This is only implemented for 3D geometries, since 2D geometries only have edges but no faces
     * @return GeometriesArrayType contains this geometry faces.
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
                                              this->pGetPoint( 2 ),
                                              this->pGetPoint( 1 ),
                                              this->pGetPoint( 8 ),
                                              this->pGetPoint( 7 ),
                                              this->pGetPoint( 6 ) ) ) );
        faces.push_back( FacePointerType( new FaceType1(
                                              this->pGetPoint( 3 ),
                                              this->pGetPoint( 4 ),
                                              this->pGetPoint( 5 ),
                                              this->pGetPoint( 12 ),
                                              this->pGetPoint( 13 ),
                                              this->pGetPoint( 14 ) ) ) );
        faces.push_back( FacePointerType( new FaceType2(
                                              this->pGetPoint( 0 ),
                                              this->pGetPoint( 1 ),
                                              this->pGetPoint( 4 ),
                                              this->pGetPoint( 3 ),
                                              this->pGetPoint( 6 ),
                                              this->pGetPoint( 10 ),
                                              this->pGetPoint( 12 ),
                                              this->pGetPoint( 9 ) ) ) );
        faces.push_back( FacePointerType( new FaceType2(
                                              this->pGetPoint( 2 ),
                                              this->pGetPoint( 0 ),
                                              this->pGetPoint( 3 ),
                                              this->pGetPoint( 5 ),
                                              this->pGetPoint( 8 ),
                                              this->pGetPoint( 9 ),
                                              this->pGetPoint( 14 ),
                                              this->pGetPoint( 11 ) ) ) );
        faces.push_back( FacePointerType( new FaceType2(
                                              this->pGetPoint( 1 ),
                                              this->pGetPoint( 2 ),
                                              this->pGetPoint( 5 ),
                                              this->pGetPoint( 4 ),
                                              this->pGetPoint( 7 ),
                                              this->pGetPoint( 11 ),
                                              this->pGetPoint( 13 ),
                                              this->pGetPoint( 10 ) ) ) );
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
     */
    double ShapeFunctionValue(
        IndexType ShapeFunctionIndex,
        const CoordinatesArrayType& rPoint
        ) const override
    {
        return CalculateShapeFunctionValue(ShapeFunctionIndex, rPoint);
    }

    using Geometry<TPointType>::ShapeFunctionsValues;

    /** This method gives gradient of all shape functions evaluated
     * in given point.
     * There is no calculation and it just give it from
     * shape functions values container if they are
     * existing. Otherwise it gives you error which this value is
     * not exist.
     * @param rResult the given Container that will be overwritten by the solution
     * @param rPoint the given local coordinates the gradients will be evaluated for
     * @return a matrix of gradients for each shape function
     */
    Matrix& ShapeFunctionsLocalGradients( Matrix& rResult, const CoordinatesArrayType& rPoint ) const override
    {
        return CalculateShapeFunctionsLocalGradients( rResult, rPoint );
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
    std::string Info() const override
    {
        return "3 dimensional prism with fiftheen nodes in 3D space";
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
        rOStream << "3 dimensional prism with fifthteen nodes in 3D space";
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

    Prism3D15(): BaseType( PointsArrayType(), &msGeometryData ) {}


    /**
     * Private Operations
     */

    /**
     * Calculates the value of a given shape function at a given point.
     *
     * @param ShapeFunctionIndex The number of the desired shape function
     * @param rPoint the given point in local coordinates at which the
     * value of the shape function is calculated
     *
     * @return the value of the shape function at the given point
     */
    static double CalculateShapeFunctionValue(
        const IndexType ShapeFunctionIndex,
        const CoordinatesArrayType& rPoint
        )
    {
        const double x = rPoint[0];
        const double y = rPoint[1];
        const double z = rPoint[2];
        switch ( ShapeFunctionIndex )
        {
            case 0 : return  (1.0/2.0)*(2*z - 2.0)*(2*z - 1)*(-2.0*x - 2.0*y + 1.0)*(-x - y + 1.0) ;
            case 1 : return  (1.0/2.0)*x*(2.0*x - 1.0)*(2*z - 2.0)*(2*z - 1) ;
            case 2 : return  (1.0/2.0)*y*(2.0*y - 1.0)*(2*z - 2.0)*(2*z - 1) ;
            case 3 : return  z*(2*z - 1)*(-2.0*x - 2.0*y + 1.0)*(-x - y + 1.0) ;
            case 4 : return  x*z*(2.0*x - 1.0)*(2*z - 1) ;
            case 5 : return  y*z*(2.0*y - 1.0)*(2*z - 1) ;
            case 6 : return  (1.0/2.0)*x*(2*z - 2.0)*(2*z - 1)*(-4.0*x - 4.0*y + 4.0) ;
            case 7 : return  2.0*x*y*(2*z - 2.0)*(2*z - 1) ;
            case 8 : return  2.0*y*(2*z - 2.0)*(2*z - 1)*(-x - y + 1.0) ;
            case 9 : return  (1.0 - std::pow(2*z - 1, 2))*(-x - y + 1.0) ;
            case 10 : return  x*(1.0 - std::pow(2*z - 1, 2)) ;
            case 11 : return  y*(1.0 - std::pow(2*z - 1, 2)) ;
            case 12 : return  x*z*(2*z - 1)*(-4.0*x - 4.0*y + 4.0) ;
            case 13 : return  4.0*x*y*z*(2*z - 1) ;
            case 14 : return  4.0*y*z*(2*z - 1)*(-x - y + 1.0) ;
            default:
                KRATOS_ERROR << "Wrong index of shape function!" << ShapeFunctionIndex << std::endl;
        }
    }

        /**
     * This method gives all non-zero shape functions values
     * evaluated at the rCoordinates provided
     * @return Vector of values of shape functions \f$ F_{i} \f$ where i is the shape function index (for NURBS it is the index of the local enumeration in the element).
     * @see ShapeFunctionValue
     * @see ShapeFunctionsLocalGradients
     * @see ShapeFunctionLocalGradient
    */
    Vector& ShapeFunctionsValues (
        Vector& rResult,
        const CoordinatesArrayType& rCoordinates
        ) const override
    {
        const double x = rCoordinates[0];
        const double y = rCoordinates[1];
        const double z = rCoordinates[2];

        rResult(  0  ) = (1.0/2.0)*(2*z - 2.0)*(2*z - 1)*(-2.0*x - 2.0*y + 1.0)*(-x - y + 1.0) ;
        rResult(  1  ) = (1.0/2.0)*x*(2.0*x - 1.0)*(2*z - 2.0)*(2*z - 1) ;
        rResult(  2  ) = (1.0/2.0)*y*(2.0*y - 1.0)*(2*z - 2.0)*(2*z - 1) ;
        rResult(  3  ) = z*(2*z - 1)*(-2.0*x - 2.0*y + 1.0)*(-x - y + 1.0) ;
        rResult(  4  ) = x*z*(2.0*x - 1.0)*(2*z - 1) ;
        rResult(  5  ) = y*z*(2.0*y - 1.0)*(2*z - 1) ;
        rResult(  6  ) = (1.0/2.0)*x*(2*z - 2.0)*(2*z - 1)*(-4.0*x - 4.0*y + 4.0) ;
        rResult(  7  ) = 2.0*x*y*(2*z - 2.0)*(2*z - 1) ;
        rResult(  8  ) = 2.0*y*(2*z - 2.0)*(2*z - 1)*(-x - y + 1.0) ;
        rResult(  9  ) = (1.0 - std::pow(2*z - 1, 2))*(-x - y + 1.0) ;
        rResult(  10  ) = x*(1.0 - std::pow(2*z - 1, 2)) ;
        rResult(  11  ) = y*(1.0 - std::pow(2*z - 1, 2)) ;
        rResult(  12  ) = x*z*(2*z - 1)*(-4.0*x - 4.0*y + 4.0) ;
        rResult(  13  ) = 4.0*x*y*z*(2*z - 1) ;
        rResult(  14  ) = 4.0*y*z*(2*z - 1)*(-x - y + 1.0) ;

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
    static Matrix& CalculateShapeFunctionsLocalGradients( Matrix& DN, const CoordinatesArrayType& rPoint )
    {
        const double x = rPoint[0];
        const double y = rPoint[1];
        const double z = rPoint[2];
        DN.resize(15,3,false);

        DN( 0 , 0 )= (1.0/2.0)*(2*z - 2.0)*(2*z - 1)*(4.0*x + 4.0*y - 3.0) ;
        DN( 0 , 1 )= (1.0/2.0)*(2*z - 2.0)*(2*z - 1)*(4.0*x + 4.0*y - 3.0) ;
        DN( 0 , 2 )= (4*z - 3.0)*(x + y - 1.0)*(2.0*x + 2.0*y - 1.0) ;
        DN( 1 , 0 )= (1.0/2.0)*(4.0*x - 1.0)*(2*z - 2.0)*(2*z - 1) ;
        DN( 1 , 1 )= 0 ;
        DN( 1 , 2 )= x*(2.0*x - 1.0)*(4*z - 3.0) ;
        DN( 2 , 0 )= 0 ;
        DN( 2 , 1 )= (1.0/2.0)*(4.0*y - 1.0)*(2*z - 2.0)*(2*z - 1) ;
        DN( 2 , 2 )= y*(2.0*y - 1.0)*(4*z - 3.0) ;
        DN( 3 , 0 )= z*(2*z - 1)*(4.0*x + 4.0*y - 3.0) ;
        DN( 3 , 1 )= z*(2*z - 1)*(4.0*x + 4.0*y - 3.0) ;
        DN( 3 , 2 )= (4*z - 1)*(x + y - 1.0)*(2.0*x + 2.0*y - 1.0) ;
        DN( 4 , 0 )= z*(4.0*x - 1.0)*(2*z - 1) ;
        DN( 4 , 1 )= 0 ;
        DN( 4 , 2 )= x*(2.0*x - 1.0)*(4*z - 1) ;
        DN( 5 , 0 )= 0 ;
        DN( 5 , 1 )= z*(4.0*y - 1.0)*(2*z - 1) ;
        DN( 5 , 2 )= y*(2.0*y - 1.0)*(4*z - 1) ;
        DN( 6 , 0 )= 2.0*(2*z - 2.0)*(2*z - 1)*(-2*x - y + 1) ;
        DN( 6 , 1 )= x*(-8.0*std::pow(z, 2) + 12.0*z - 4.0) ;
        DN( 6 , 2 )= 4.0*x*(3.0 - 4*z)*(x + y - 1) ;
        DN( 7 , 0 )= y*(8.0*std::pow(z, 2) - 12.0*z + 4.0) ;
        DN( 7 , 1 )= x*(8.0*std::pow(z, 2) - 12.0*z + 4.0) ;
        DN( 7 , 2 )= x*y*(16.0*z - 12.0) ;
        DN( 8 , 0 )= y*(-8.0*std::pow(z, 2) + 12.0*z - 4.0) ;
        DN( 8 , 1 )= -(2*z - 2.0)*(2.0*y*(2*z - 1) + (4.0*z - 2.0)*(x + y - 1.0)) ;
        DN( 8 , 2 )= 4.0*y*(3.0 - 4*z)*(x + y - 1.0) ;
        DN( 9 , 0 )= 4*z*(z - 1) ;
        DN( 9 , 1 )= 4*z*(z - 1) ;
        DN( 9 , 2 )= 4*(2*z - 1)*(x + y - 1.0) ;
        DN( 10 , 0 )= 4*z*(1 - z) ;
        DN( 10 , 1 )= 0 ;
        DN( 10 , 2 )= 4*x*(1 - 2*z) ;
        DN( 11 , 0 )= 0 ;
        DN( 11 , 1 )= 4*z*(1 - z) ;
        DN( 11 , 2 )= 4*y*(1 - 2*z) ;
        DN( 12 , 0 )= 4.0*z*(2*z - 1)*(-2*x - y + 1) ;
        DN( 12 , 1 )= x*z*(4.0 - 8.0*z) ;
        DN( 12 , 2 )= x*(4.0 - 16.0*z)*(x + y - 1) ;
        DN( 13 , 0 )= y*z*(8.0*z - 4.0) ;
        DN( 13 , 1 )= x*z*(8.0*z - 4.0) ;
        DN( 13 , 2 )= x*y*(16.0*z - 4.0) ;
        DN( 14 , 0 )= y*z*(4.0 - 8.0*z) ;
        DN( 14 , 1 )= 4.0*z*(2*z - 1)*(-x - 2*y + 1.0) ;
        DN( 14 , 2 )= y*(4.0 - 16.0*z)*(x + y - 1.0) ;

        return DN;
    }

    /**
     * Calculates the values of all shape function in all integration points.
     * Integration points are expected to be given in local coordinates
     *
     * @param ThisMethod the current integration method
     * @return the matrix of values of every shape function in each integration point
     *
     */
    static Matrix CalculateShapeFunctionsIntegrationPointsValues(typename BaseType::IntegrationMethod ThisMethod )
    {
        IntegrationPointsContainerType all_integration_points = AllIntegrationPoints();
        IntegrationPointsArrayType integration_points = all_integration_points[static_cast<int>(ThisMethod)];

        // Number of integration points
        const std::size_t integration_points_number = integration_points.size();

        //Setting up return matrix
        Matrix shape_function_values( integration_points_number, 15 );

        // Loop over all integration points
        double x, y, z;
        for ( std::size_t pnt = 0; pnt < integration_points_number; pnt++ ) {
            const auto& r_point = integration_points[pnt];
            x = r_point[0];
            y = r_point[1];
            z = r_point[2];

            shape_function_values( pnt,  0  ) = (1.0/2.0)*(2*z - 2.0)*(2*z - 1)*(-2.0*x - 2.0*y + 1.0)*(-x - y + 1.0) ;
            shape_function_values( pnt,  1  ) = (1.0/2.0)*x*(2.0*x - 1.0)*(2*z - 2.0)*(2*z - 1) ;
            shape_function_values( pnt,  2  ) = (1.0/2.0)*y*(2.0*y - 1.0)*(2*z - 2.0)*(2*z - 1) ;
            shape_function_values( pnt,  3  ) = z*(2*z - 1)*(-2.0*x - 2.0*y + 1.0)*(-x - y + 1.0) ;
            shape_function_values( pnt,  4  ) = x*z*(2.0*x - 1.0)*(2*z - 1) ;
            shape_function_values( pnt,  5  ) = y*z*(2.0*y - 1.0)*(2*z - 1) ;
            shape_function_values( pnt,  6  ) = (1.0/2.0)*x*(2*z - 2.0)*(2*z - 1)*(-4.0*x - 4.0*y + 4.0) ;
            shape_function_values( pnt,  7  ) = 2.0*x*y*(2*z - 2.0)*(2*z - 1) ;
            shape_function_values( pnt,  8  ) = 2.0*y*(2*z - 2.0)*(2*z - 1)*(-x - y + 1.0) ;
            shape_function_values( pnt,  9  ) = (1.0 - std::pow(2*z - 1, 2))*(-x - y + 1.0) ;
            shape_function_values( pnt,  10  ) = x*(1.0 - std::pow(2*z - 1, 2)) ;
            shape_function_values( pnt,  11  ) = y*(1.0 - std::pow(2*z - 1, 2)) ;
            shape_function_values( pnt,  12  ) = x*z*(2*z - 1)*(-4.0*x - 4.0*y + 4.0) ;
            shape_function_values( pnt,  13  ) = 4.0*x*y*z*(2*z - 1) ;
            shape_function_values( pnt,  14  ) = 4.0*y*z*(2*z - 1)*(-x - y + 1.0) ;
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
    static ShapeFunctionsGradientsType CalculateShapeFunctionsIntegrationPointsLocalGradients(typename BaseType::IntegrationMethod ThisMethod )
    {
        IntegrationPointsContainerType all_integration_points = AllIntegrationPoints();
        IntegrationPointsArrayType integration_points = all_integration_points[static_cast<int>(ThisMethod)];

        // Number of integration points
        const std::size_t integration_points_number = integration_points.size();
        ShapeFunctionsGradientsType d_shape_f_values( integration_points_number );

        // Initialising container
        Matrix result = ZeroMatrix( 15, 3 );

        // Loop over all integration points
        for ( std::size_t pnt = 0; pnt < integration_points_number; ++pnt ) {
            CalculateShapeFunctionsLocalGradients(result, integration_points[pnt]);
            d_shape_f_values[pnt] = result;
        }

        return d_shape_f_values;
    }

    static const IntegrationPointsContainerType AllIntegrationPoints()
    {
        IntegrationPointsContainerType integration_points =
        {
            {
                Quadrature < PrismGaussLegendreIntegrationPoints1,
                3, IntegrationPoint<3> >::GenerateIntegrationPoints(),
                Quadrature < PrismGaussLegendreIntegrationPoints2,
                3, IntegrationPoint<3> >::GenerateIntegrationPoints(),
                Quadrature < PrismGaussLegendreIntegrationPoints3,
                3, IntegrationPoint<3> >::GenerateIntegrationPoints(),
                Quadrature < PrismGaussLegendreIntegrationPoints4,
                3, IntegrationPoint<3> >::GenerateIntegrationPoints(),
                Quadrature < PrismGaussLegendreIntegrationPoints5,
                3, IntegrationPoint<3> >::GenerateIntegrationPoints(),
                Quadrature < PrismGaussLegendreIntegrationPointsExt1,
                3, IntegrationPoint<3> >::GenerateIntegrationPoints(),
                Quadrature < PrismGaussLegendreIntegrationPointsExt2,
                3, IntegrationPoint<3> >::GenerateIntegrationPoints(),
                Quadrature < PrismGaussLegendreIntegrationPointsExt3,
                3, IntegrationPoint<3> >::GenerateIntegrationPoints(),
                Quadrature < PrismGaussLegendreIntegrationPointsExt4,
                3, IntegrationPoint<3> >::GenerateIntegrationPoints(),
                Quadrature < PrismGaussLegendreIntegrationPointsExt5,
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
                Prism3D15<TPointType>::CalculateShapeFunctionsIntegrationPointsValues(
                    GeometryData::IntegrationMethod::GI_GAUSS_1 ),
                Prism3D15<TPointType>::CalculateShapeFunctionsIntegrationPointsValues(
                    GeometryData::IntegrationMethod::GI_GAUSS_2 ),
                Prism3D15<TPointType>::CalculateShapeFunctionsIntegrationPointsValues(
                    GeometryData::IntegrationMethod::GI_GAUSS_3 ),
                Prism3D15<TPointType>::CalculateShapeFunctionsIntegrationPointsValues(
                    GeometryData::IntegrationMethod::GI_GAUSS_4 ),
                Prism3D15<TPointType>::CalculateShapeFunctionsIntegrationPointsValues(
                    GeometryData::IntegrationMethod::GI_GAUSS_5 ),
                Prism3D15<TPointType>::CalculateShapeFunctionsIntegrationPointsValues(
                    GeometryData::IntegrationMethod::GI_EXTENDED_GAUSS_1 ),
                Prism3D15<TPointType>::CalculateShapeFunctionsIntegrationPointsValues(
                    GeometryData::IntegrationMethod::GI_EXTENDED_GAUSS_2 ),
                Prism3D15<TPointType>::CalculateShapeFunctionsIntegrationPointsValues(
                    GeometryData::IntegrationMethod::GI_EXTENDED_GAUSS_3 ),
                Prism3D15<TPointType>::CalculateShapeFunctionsIntegrationPointsValues(
                    GeometryData::IntegrationMethod::GI_EXTENDED_GAUSS_4 ),
                Prism3D15<TPointType>::CalculateShapeFunctionsIntegrationPointsValues(
                    GeometryData::IntegrationMethod::GI_EXTENDED_GAUSS_5 )
            }
        };
        return shape_functions_values;
    }

    static const ShapeFunctionsLocalGradientsContainerType AllShapeFunctionsLocalGradients()
    {
        ShapeFunctionsLocalGradientsContainerType shape_functions_local_gradients =
        {
            {
                Prism3D15<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients(
                    GeometryData::IntegrationMethod::GI_GAUSS_1 ),
                Prism3D15<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients(
                    GeometryData::IntegrationMethod::GI_GAUSS_2 ),
                Prism3D15<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients(
                    GeometryData::IntegrationMethod::GI_GAUSS_3 ),
                Prism3D15<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients(
                    GeometryData::IntegrationMethod::GI_GAUSS_4 ),
                Prism3D15<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients(
                    GeometryData::IntegrationMethod::GI_GAUSS_5 ),
                Prism3D15<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients(
                    GeometryData::IntegrationMethod::GI_EXTENDED_GAUSS_1 ),
                Prism3D15<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients(
                    GeometryData::IntegrationMethod::GI_EXTENDED_GAUSS_2 ),
                Prism3D15<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients(
                    GeometryData::IntegrationMethod::GI_EXTENDED_GAUSS_3 ),
                Prism3D15<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients(
                    GeometryData::IntegrationMethod::GI_EXTENDED_GAUSS_4 ),
                Prism3D15<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients(
                    GeometryData::IntegrationMethod::GI_EXTENDED_GAUSS_5 )
            }
        };
        return shape_functions_local_gradients;
    }

    /**
     * Private Friends
     */

    template<class TOtherPointType> friend class Prism3D15;


    /**
     * Un accessible methods
     */

};// Class Prism3D15


/**
 * Input and output
 */

/// input stream function
template<class TPointType> inline std::istream& operator >> (
    std::istream& rIStream, Prism3D15<TPointType>& rThis );

/// output stream function
template<class TPointType> inline std::ostream& operator << (
    std::ostream& rOStream, const Prism3D15<TPointType>& rThis )
{
    rThis.PrintInfo( rOStream );
    rOStream << std::endl;
    rThis.PrintData( rOStream );

    return rOStream;
}


template<class TPointType> const
GeometryData Prism3D15<TPointType>::msGeometryData(
    &msGeometryDimension,
    GeometryData::IntegrationMethod::GI_GAUSS_3,
    Prism3D15<TPointType>::AllIntegrationPoints(),
    Prism3D15<TPointType>::AllShapeFunctionsValues(),
    AllShapeFunctionsLocalGradients()
);

template<class TPointType> const
GeometryDimension Prism3D15<TPointType>::msGeometryDimension(3, 3);

}// namespace Kratos.
