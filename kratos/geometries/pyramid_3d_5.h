//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher (https://github.com/philbucher)
//  Contributors:    Ashish Darekar

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
 *               ,/   | | \
 *             ,/    .' | `.
 *           ,/      |  '.  \
 *         ,/       .' w |   \
 *       ,/         |  ^ |    \
 *      3----------.'--|-2    `.
 *       `\        |   |  `\    \
 *         `\     .'   +----`\ - \ -> v
 *           `\   |    `\     `\  \
 *             `\.'      `\     `\`
 *                0----------------1
 *                          `\
 *                             u
 * @author Philipp Bucher
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

    //typedef typename::BaseType::Line3D2<TPointType> EdgeType;


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
    // Pyramid3D5(Pyramid3D5 const& rOther)
    // : BaseType(rOther)
    // {
    // }

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
    // template<class TOtherPointType> Pyramid3D5(Pyramid3D5<TOtherPointType> const& rOther)
    // : BaseType(rOther)
    // {
    // }

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
     * :TODO: TO BE TESTED
     */
    /** This method calculate and return the height of the pyramid.
     @return double value contain the height
     @ref (https://math.stackexchange.com/questions/2177006/how-to-define-a-plane-based-on-4-points)
     @ref (https://stackoverflow.com/questions/9605556/how-to-project-a-point-onto-a-plane-in-3d)
     */
    double HeightOfPyramid() const
    {
        //Step 1 : find the Unit Normal Vector of the Base using 4 points A,B,C,D
        const PointsArrayType& vertices = this->Points();

        array_1d<double,3> AB;
        array_1d<double,3> BC;
        array_1d<double,3> UnitNormal;

        AB[0] = vertices[1].X() - vertices[0].X();
        AB[1] = vertices[1].Y() - vertices[0].Y();
        AB[2] = vertices[1].Z() - vertices[0].Z();

        BC[0] = vertices[2].X() - vertices[1].X();
        BC[1] = vertices[2].Y() - vertices[1].Y();
        BC[2] = vertices[2].Z() - vertices[1].Z();

        MathUtils<double>::UnitCrossProduct(UnitNormal, AB, BC);

        //Step2 : Make a Vector from any Base point to the apex point
        array_1d<double,3> AE;

        AE[0] = vertices[4].X() - vertices[0].X();
        AE[1] = vertices[4].Y() - vertices[0].Y();
        AE[2] = vertices[4].Z() - vertices[0].Z();

        //Step3 : Take a dot product of this Vector with the Unit Normal Vector to get height
        double height = AE[0] * UnitNormal[0] + AE[1] * UnitNormal[1] + AE[2] * UnitNormal[2] ;

        return height;
    }


    /**
     * :TODO: TO BE TESTED
     */
    /**
     * Determinant of jacobians for given integration method.
     * This method calculate determinant of jacobian in all
     * integrations points of given integration method.
     *
     * @return Vector of double which is vector of determinants of
     * jacobians \f$ |J|_i \f$ where \f$ i=1,2,...,n \f$ is the
     * integration point index of given integration method.
     *
     * @see Jacobian
     * @see InverseOfJacobian
     */

    Vector& DeterminantOfJacobian( Vector& rResult,
                                   IntegrationMethod ThisMethod ) const override
    {
        double side1, side2, height;

        const PointsArrayType& vertices = this->Points();

        side1 = std::sqrt(std::pow((vertices[0].X() - vertices[1].X()), 2.0) +
            std::pow((vertices[0].Y() - vertices[1].Y()), 2.0) + std::pow((vertices[0].Z() - vertices[1].Z()), 2.0));

        side2 = std::sqrt(std::pow((vertices[2].X() - vertices[1].X()), 2.0) +
            std::pow((vertices[2].Y() - vertices[1].Y()), 2.0) + std::pow((vertices[2].Z() - vertices[1].Z()), 2.0));

        height = HeightOfPyramid(); // to get Height of pyramid

        IntegrationPointsContainerType all_integration_points = AllIntegrationPoints();
        IntegrationPointsArrayType integration_points = all_integration_points[ThisMethod];

        if( rResult.size() != this->IntegrationPointsNumber( ThisMethod ) )
            rResult.resize( this->IntegrationPointsNumber( ThisMethod ), false );

        for ( unsigned int pnt = 0; pnt < this->IntegrationPointsNumber( ThisMethod ); pnt++ )
        {
            rResult[pnt] = (0.03125) * side1 * side2 * height * (1 - integration_points[pnt].Z()) * (1 - integration_points[pnt].Z()) ;
        }
        return rResult;
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

        rResult( 0, 0 ) = +1.0;
        rResult( 0, 1 ) = -1.0;
        rResult( 0, 2 ) =  0.0;

        rResult( 1, 0 ) = +1.0;
        rResult( 1, 1 ) = +1.0;
        rResult( 1, 2 ) =  0.0;

        rResult( 2, 0 ) = -1.0;
        rResult( 2, 1 ) =  1.0;
        rResult( 2, 2 ) =  0.0;

        rResult( 3, 0 ) = -1.0;
        rResult( 3, 1 ) = -1.0;
        rResult( 3, 2 ) = 0.0;

        rResult( 4, 0 ) = 0.0;
        rResult( 4, 1 ) = 0.0;
        rResult( 4, 2 ) = 1.0;

        return rResult;
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
    double ShapeFunctionValue( IndexType ShapeFunctionIndex,
                                       const CoordinatesArrayType& rPoint ) const override
    {
        switch ( ShapeFunctionIndex )
        {
        case 0:
            return( (0.125) * (1 - rPoint[1]) * (1 + rPoint[0]) * (1 + rPoint[2]) );
        case 1:
            return( (0.125) * (1 + rPoint[1]) * (1 + rPoint[0]) * (1 + rPoint[2]) );
        case 2:
            return( (0.125) * (1 + rPoint[1]) * (1 - rPoint[0]) * (1 + rPoint[2]) );
        case 3:
            return( (0.125) * (1 - rPoint[1]) * (1 - rPoint[0]) * (1 + rPoint[2]) );
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
            shape_function_values( pnt, 0 ) = (0.125) * (1 - integration_points[pnt].Y()) * (1 + integration_points[pnt].X()) * (1 + integration_points[pnt].Z()) ;
            shape_function_values( pnt, 1 ) = (0.125) * (1 + integration_points[pnt].Y()) * (1 + integration_points[pnt].X()) * (1 + integration_points[pnt].Z()) ;
            shape_function_values( pnt, 2 ) = (0.125) * (1 + integration_points[pnt].Y()) * (1 - integration_points[pnt].X()) * (1 + integration_points[pnt].Z()) ;
            shape_function_values( pnt, 3 ) = (0.125) * (1 - integration_points[pnt].Y()) * (1 - integration_points[pnt].X()) * (1 + integration_points[pnt].Z()) ;
            shape_function_values( pnt, 4 ) = (0.5) * (1 + integration_points[pnt].Z()) ;
        }

        return shape_function_values;
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
        rResult( 0, 0 ) =  (+0.125) * ( 1 - rPoint[1] ) * ( 1 + rPoint[2] ) ;
        rResult( 0, 1 ) =  (-0.125) * ( 1 + rPoint[0] ) * ( 1 + rPoint[2] ) ;
        rResult( 0, 2 ) =  (+0.125) * ( 1 - rPoint[1] ) * ( 1 + rPoint[0] ) ;

        rResult( 1, 0 ) =  (+0.125) * ( 1 + rPoint[1] ) * ( 1 + rPoint[2] ) ;
        rResult( 1, 1 ) =  (+0.125) * ( 1 + rPoint[0] ) * ( 1 + rPoint[2] ) ;
        rResult( 1, 2 ) =  (+0.125) * ( 1 + rPoint[1] ) * ( 1 + rPoint[0] ) ;

        rResult( 2, 0 ) =  (-0.125) * ( 1 + rPoint[1] ) * ( 1 + rPoint[2] ) ;
        rResult( 2, 1 ) =  (+0.125) * ( 1 - rPoint[0] ) * ( 1 + rPoint[2] ) ;
        rResult( 2, 2 ) =  (+0.125) * ( 1 + rPoint[1] ) * ( 1 - rPoint[0] ) ;

        rResult( 3, 0 ) =  (-0.125) * ( 1 - rPoint[1] ) * ( 1 + rPoint[2] ) ;
        rResult( 3, 1 ) =  (-0.125) * ( 1 - rPoint[0] ) * ( 1 + rPoint[2] ) ;
        rResult( 3, 2 ) =  (+0.125) * ( 1 - rPoint[1] ) * ( 1 - rPoint[0] ) ;

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

            result( 0, 0 ) =  (+0.125) * ( 1 - integration_points[pnt].Y() ) * ( 1 + integration_points[pnt].Z() ) ;
            result( 0, 1 ) =  (-0.125) * ( 1 + integration_points[pnt].X() ) * ( 1 + integration_points[pnt].Z() ) ;
            result( 0, 2 ) =  (+0.125) * ( 1 - integration_points[pnt].Y() ) * ( 1 + integration_points[pnt].X() ) ;

            result( 1, 0 ) =  (+0.125) * ( 1 + integration_points[pnt].Y() ) * ( 1 + integration_points[pnt].Z() ) ;
            result( 1, 1 ) =  (+0.125) * ( 1 + integration_points[pnt].X() ) * ( 1 + integration_points[pnt].Z() ) ;
            result( 1, 2 ) =  (+0.125) * ( 1 + integration_points[pnt].Y() ) * ( 1 + integration_points[pnt].X() ) ;

            result( 2, 0 ) =  (-0.125) * ( 1 + integration_points[pnt].Y() ) * ( 1 + integration_points[pnt].Z() ) ;
            result( 2, 1 ) =  (+0.125) * ( 1 - integration_points[pnt].X() ) * ( 1 + integration_points[pnt].Z() ) ;
            result( 2, 2 ) =  (+0.125) * ( 1 + integration_points[pnt].Y() ) * ( 1 - integration_points[pnt].X() ) ;

            result( 3, 0 ) =  (-0.125) * ( 1 - integration_points[pnt].Y() ) * ( 1 + integration_points[pnt].Z() ) ;
            result( 3, 1 ) =  (-0.125) * ( 1 - integration_points[pnt].X() ) * ( 1 + integration_points[pnt].Z() ) ;
            result( 3, 2 ) =  (+0.125) * ( 1 - integration_points[pnt].Y() ) * ( 1 - integration_points[pnt].X() ) ;

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
        return "3 dimensional pyramid with four nodes in 3D space";
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
