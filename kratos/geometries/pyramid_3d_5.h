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
//

#if !defined (KRATOS_PYRAMID_3D_5_H_INCLUDED)
#define KRATOS_PYRAMID_3D_5_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "geometries/geometry.h"
#include "integration/pyramid_gauss_legendre_integration_points.h"


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
    ///@name Informations
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
        return 8;
    }

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
        Vector temp;
        this->DeterminantOfJacobian(temp, msGeometryData.DefaultIntegrationMethod());
        const IntegrationPointsArrayType& integration_points = this->IntegrationPoints(msGeometryData.DefaultIntegrationMethod());
        double vol = 0.00;

        for (std::size_t i=0; i<integration_points.size(); ++i) {
            vol += temp[i] * integration_points[i].Weight();
        }

        return vol;
    }


    /**
     * Returns whether given arbitrary point is inside the Geometry and the respective
     * local point for the given global point
     * @param rPoint The point to be checked if is inside o note in global coordinates
     * @param rResult The local coordinates of the point
     * @param Tolerance The  tolerance that will be considered to check if the point is inside or not
     * @ref (https://stackoverflow.com/questions/68641598/check-if-a-3d-point-is-in-a-square-based-pyramid-or-not)
     * @return True if the point is inside, false otherwise
     */
    bool IsInside(
        const CoordinatesArrayType& rPoint,
        CoordinatesArrayType& rResult,
        const double Tolerance = std::numeric_limits<double>::epsilon()
        ) const override
    {
        this->PointLocalCoordinates( rResult, rPoint );

        //Calculation of all the surface normals of a Pyramid ABCDE
        Matrix LocalPyramid;
        LocalPyramid = PointsLocalCoordinates(LocalPyramid); //Fix values

        //Surface Normal for ABE
        array_1d<double,3> AB;
        array_1d<double,3> AE;
        array_1d<double,3> nABE;

        //(B - A)
        AB[0] = LocalPyramid(1,0) - LocalPyramid(0,0);
        AB[1] = LocalPyramid(1,1) - LocalPyramid(0,1);
        AB[2] = LocalPyramid(1,2) - LocalPyramid(0,2);

        //(E - A)
        AE[0] = LocalPyramid(4,0) - LocalPyramid(0,0);
        AE[1] = LocalPyramid(4,1) - LocalPyramid(0,1);
        AE[2] = LocalPyramid(4,2) - LocalPyramid(0,2);

        MathUtils<double>::UnitCrossProduct(nABE, AB, AE);

        //Surface Normal for BCE
        array_1d<double,3> BC;
        array_1d<double,3> BE;
        array_1d<double,3> nBCE;

        //(C - B)
        BC[0] = LocalPyramid(2,0) - LocalPyramid(1,0);
        BC[1] = LocalPyramid(2,1) - LocalPyramid(1,1);
        BC[2] = LocalPyramid(2,2) - LocalPyramid(1,2);

        //(E - B)
        BE[0] = LocalPyramid(4,0) - LocalPyramid(1,0);
        BE[1] = LocalPyramid(4,1) - LocalPyramid(1,1);
        BE[2] = LocalPyramid(4,2) - LocalPyramid(1,2);

        MathUtils<double>::UnitCrossProduct(nBCE, BC, BE);

        //Surface Normal for CDE
        array_1d<double,3> CD;
        array_1d<double,3> CE;
        array_1d<double,3> nCDE;

        //(D - C)
        CD[0] = LocalPyramid(3,0) - LocalPyramid(2,0);
        CD[1] = LocalPyramid(3,1) - LocalPyramid(2,1);
        CD[2] = LocalPyramid(3,2) - LocalPyramid(2,2);

        //(E - C)
        CE[0] = LocalPyramid(4,0) - LocalPyramid(2,0);
        CE[1] = LocalPyramid(4,1) - LocalPyramid(2,1);
        CE[2] = LocalPyramid(4,2) - LocalPyramid(2,2);

        MathUtils<double>::UnitCrossProduct(nCDE, CD, CE);

        //Surface Normal for DAE
        array_1d<double,3> DA;
        array_1d<double,3> DE;
        array_1d<double,3> nDAE;

        //(A - D)
        DA[0] = LocalPyramid(0,0) - LocalPyramid(3,0);
        DA[1] = LocalPyramid(0,1) - LocalPyramid(3,1);
        DA[2] = LocalPyramid(0,2) - LocalPyramid(3,2);

        //(E - D)
        DE[0] = LocalPyramid(4,0) - LocalPyramid(3,0);
        DE[1] = LocalPyramid(4,1) - LocalPyramid(3,1);
        DE[2] = LocalPyramid(4,2) - LocalPyramid(3,2);

        MathUtils<double>::UnitCrossProduct(nDAE, DA, DE);

        //Surface Normal for ABCD , using AB and DA
        array_1d<double,3> nABCD;

        MathUtils<double>::UnitCrossProduct(nABCD, AB, DA);

        //Direction Vector from Point to the point in the Plane
        array_1d<double,3> PE; //As point E is been shared between all 4 planes ABE, BCE, CDE, DAE
        array_1d<double,3> PD; //For the base place ABCD

        //PE = (E - P)
        PE[0] = LocalPyramid(4,0) - rResult[0];
        PE[1] = LocalPyramid(4,1) - rResult[1];
        PE[2] = LocalPyramid(4,2) - rResult[2];

        //PD = (D - P)
        PD[0] = LocalPyramid(3,0) - rResult[0];
        PD[1] = LocalPyramid(3,1) - rResult[1];
        PD[2] = LocalPyramid(3,2) - rResult[2];

        //Dot products of direction vector with Planar normal
        return ((MathUtils<double>::Dot(PE, nABE) >= 0) && (MathUtils<double>::Dot(PE, nBCE) >= 0) && (MathUtils<double>::Dot(PE, nCDE) >= 0) && (MathUtils<double>::Dot(PE, nDAE) >= 0) && (MathUtils<double>::Dot(PD, nABCD) >= 0));
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
        if ( std::abs( rPointLocalCoordinates[0] ) <= (1.0 + Tolerance) )
        {
            if ( std::abs( rPointLocalCoordinates[1] ) <= (1.0 + Tolerance) )
            {
                if ( std::abs( rPointLocalCoordinates[2] ) <= (1.0 + Tolerance) )
                {
                    return 1;
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
        if(rResult.size() != 5)
            rResult.resize(5,false);

        rResult[0] = (0.125) * (1 - rCoordinates[0]) * (1 - rCoordinates[1]) * (1 - rCoordinates[2]);
        rResult[1] = (0.125) * (1 + rCoordinates[0]) * (1 - rCoordinates[1]) * (1 - rCoordinates[2]);
        rResult[2] = (0.125) * (1 + rCoordinates[0]) * (1 + rCoordinates[1]) * (1 - rCoordinates[2]);
        rResult[3] = (0.125) * (1 - rCoordinates[0]) * (1 + rCoordinates[1]) * (1 - rCoordinates[2]);
        rResult[4] = (0.5) * (1 + rCoordinates[2]);

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
            KRATOS_ERROR << "Wrong index of shape function!" << *this  << std::endl;
        }

        return 0;
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
        IntegrationPointsArrayType integration_points = all_integration_points[ThisMethod];
        //number of integration points
        const int integration_points_number = integration_points.size();
        //number of nodes in current geometry
        const int points_number = 5;
        //setting up return matrix
        Matrix shape_function_values( integration_points_number, points_number );
        //loop over all integration points

        for ( int pnt = 0; pnt < integration_points_number; pnt++ ) {
            shape_function_values( pnt, 0 ) = (0.125) * (1 - integration_points[pnt].X()) * (1 - integration_points[pnt].Y()) * (1 - integration_points[pnt].Z()) ;
            shape_function_values( pnt, 1 ) = (0.125) * (1 + integration_points[pnt].X()) * (1 - integration_points[pnt].Y()) * (1 - integration_points[pnt].Z()) ;
            shape_function_values( pnt, 2 ) = (0.125) * (1 + integration_points[pnt].X()) * (1 + integration_points[pnt].Y()) * (1 - integration_points[pnt].Z()) ;
            shape_function_values( pnt, 3 ) = (0.125) * (1 - integration_points[pnt].X()) * (1 + integration_points[pnt].Y()) * (1 - integration_points[pnt].Z()) ;
            shape_function_values( pnt, 4 ) = (0.5) * (1 + integration_points[pnt].Z()) ;
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
        IntegrationPointsArrayType integration_points = all_integration_points[ThisMethod]; //number of integration points

        const int integration_points_number = integration_points.size();
        ShapeFunctionsGradientsType d_shape_f_values( integration_points_number ); //initialising container

        //loop over all integration points
        for ( int pnt = 0; pnt < integration_points_number; pnt++ )
        {
            Matrix result = ZeroMatrix( 5, 3 );

            result( 0, 0 ) =  (-0.125) * ( 1 - integration_points[pnt].Y() ) * ( 1 - integration_points[pnt].Z() ) ;
            result( 0, 1 ) =  (-0.125) * ( 1 - integration_points[pnt].X() ) * ( 1 - integration_points[pnt].Z() ) ;
            result( 0, 2 ) =  (-0.125) * ( 1 - integration_points[pnt].X() ) * ( 1 - integration_points[pnt].Y() ) ;

            result( 1, 0 ) =  (+0.125) * ( 1 - integration_points[pnt].Y() ) * ( 1 - integration_points[pnt].Z() ) ;
            result( 1, 1 ) =  (-0.125) * ( 1 + integration_points[pnt].X() ) * ( 1 - integration_points[pnt].Z() ) ;
            result( 1, 2 ) =  (-0.125) * ( 1 + integration_points[pnt].X() ) * ( 1 - integration_points[pnt].Y() ) ;

            result( 2, 0 ) =  (+0.125) * ( 1 + integration_points[pnt].Y() ) * ( 1 - integration_points[pnt].Z() ) ;
            result( 2, 1 ) =  (+0.125) * ( 1 + integration_points[pnt].X() ) * ( 1 - integration_points[pnt].Z() ) ;
            result( 2, 2 ) =  (-0.125) * ( 1 + integration_points[pnt].X() ) * ( 1 + integration_points[pnt].Y() ) ;

            result( 3, 0 ) =  (-0.125) * ( 1 + integration_points[pnt].Y() ) * ( 1 - integration_points[pnt].Z() ) ;
            result( 3, 1 ) =  (+0.125) * ( 1 - integration_points[pnt].X() ) * ( 1 - integration_points[pnt].Z() ) ;
            result( 3, 2 ) =  (-0.125) * ( 1 - integration_points[pnt].X() ) * ( 1 + integration_points[pnt].Y() ) ;

            result( 4, 0 ) =   0.00 ;
            result( 4, 1 ) =   0.00 ;
            result( 4, 2 ) =  +0.50 ;

            d_shape_f_values[pnt] = result;
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
                Pyramid3D5<TPointType>::CalculateShapeFunctionsIntegrationPointsValues(GeometryData::GI_GAUSS_1),
                Pyramid3D5<TPointType>::CalculateShapeFunctionsIntegrationPointsValues(GeometryData::GI_GAUSS_2),
                Pyramid3D5<TPointType>::CalculateShapeFunctionsIntegrationPointsValues(GeometryData::GI_GAUSS_3),
                Pyramid3D5<TPointType>::CalculateShapeFunctionsIntegrationPointsValues(GeometryData::GI_GAUSS_4),
                Pyramid3D5<TPointType>::CalculateShapeFunctionsIntegrationPointsValues(GeometryData::GI_GAUSS_5)
            }
        };
        return shape_functions_values;
    }

    /**
     * TODO: TO BE VERIFIED
     */
    static const ShapeFunctionsLocalGradientsContainerType AllShapeFunctionsLocalGradients()
    {
        ShapeFunctionsLocalGradientsContainerType shape_functions_local_gradients =
        {
            {
                Pyramid3D5<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients(GeometryData::GI_GAUSS_1),
                Pyramid3D5<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients(GeometryData::GI_GAUSS_2),
                Pyramid3D5<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients(GeometryData::GI_GAUSS_3),
                Pyramid3D5<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients(GeometryData::GI_GAUSS_4),
                Pyramid3D5<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients(GeometryData::GI_GAUSS_5)
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
    GeometryData::GI_GAUSS_2,
    Pyramid3D5<TPointType>::AllIntegrationPoints(),
    Pyramid3D5<TPointType>::AllShapeFunctionsValues(),
    AllShapeFunctionsLocalGradients()
);

template<class TPointType> const
GeometryDimension Pyramid3D5<TPointType>::msGeometryDimension(
    3, 3, 3);

}  // namespace Kratos.

#endif // KRATOS_PYRAMID_3D_5_H_INCLUDED defined
