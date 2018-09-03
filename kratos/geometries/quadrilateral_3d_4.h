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
//                   Bodhinanda Chandra
//

#if !defined(KRATOS_QUADRILATERAL_3D_4_H_INCLUDED )
#define  KRATOS_QUADRILATERAL_3D_4_H_INCLUDED

// System includes

// External includes

// Project includes
#include "geometries/line_3d_2.h"
#include "geometries/triangle_3d_3.h"
#include "integration/quadrilateral_gauss_legendre_integration_points.h"
#include "integration/quadrilateral_collocation_integration_points.h"

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
 * @class Quadrilateral3D4
 * @ingroup KratosCore
 * @brief A four node 3D quadrilateral geometry with bi-linear shape functions
 * @details While the shape functions are only defined in 2D it is possible to define an arbitrary orientation in space. Thus it can be used for defining surfaces on 3D elements.
 * The node ordering corresponds with: 
 *            v
 *            ^
 *            |
 *      3-----------2 
 *      |     |     |         
 *      |     |     |          
 *      |     +---- | --> u    
 *      |           |          
 *      |           |          
 *      0-----------1      
 * @author Riccardo Rossi
 * @author Janosch Stascheit
 * @author Felix Nagel
 */
template<class TPointType> class Quadrilateral3D4
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
    typedef Line3D2<TPointType> EdgeType;

    /**
     * Pointer definition of Quadrilateral3D4
     */
    KRATOS_CLASS_POINTER_DEFINITION( Quadrilateral3D4 );

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
    * A third order tensor to hold shape functions' local second derivatives.
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

//     Quadrilateral3D4( const PointType& FirstPoint,
//                       const PointType& SecondPoint,
//                       const PointType& ThirdPoint,
//                       const PointType& FourthPoint )
//         : BaseType( PointsArrayType(), &msGeometryData )
//     {
//         this->Points().push_back( typename PointType::Pointer( new PointType( FirstPoint ) ) );
//         this->Points().push_back( typename PointType::Pointer( new PointType( SecondPoint ) ) );
//         this->Points().push_back( typename PointType::Pointer( new PointType( ThirdPoint ) ) );
//         this->Points().push_back( typename PointType::Pointer( new PointType( FourthPoint ) ) );
//     }

    Quadrilateral3D4( typename PointType::Pointer pFirstPoint,
                      typename PointType::Pointer pSecondPoint,
                      typename PointType::Pointer pThirdPoint,
                      typename PointType::Pointer pFourthPoint )
        : BaseType( PointsArrayType(), &msGeometryData )
    {
        this->Points().push_back( pFirstPoint );
        this->Points().push_back( pSecondPoint );
        this->Points().push_back( pThirdPoint );
        this->Points().push_back( pFourthPoint );
    }

    Quadrilateral3D4( const PointsArrayType& ThisPoints )
        : BaseType( ThisPoints, &msGeometryData )
    {
        if ( this->PointsNumber() != 4 )
            KRATOS_ERROR << "Invalid points number. Expected 4, given " << this->PointsNumber() << std::endl;
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
    Quadrilateral3D4( Quadrilateral3D4 const& rOther )
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
    template<class TOtherPointType> Quadrilateral3D4( Quadrilateral3D4<TOtherPointType> const& rOther )
        : BaseType( rOther )
    {
    }

    /**
     * Destructor. Does nothing!!!
     */
    ~Quadrilateral3D4() override {}

    GeometryData::KratosGeometryFamily GetGeometryFamily() const override
    {
        return GeometryData::Kratos_Quadrilateral;
    }

    GeometryData::KratosGeometryType GetGeometryType() const override
    {
        return GeometryData::Kratos_Quadrilateral3D4;
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
    Quadrilateral3D4& operator=( const Quadrilateral3D4& rOther )
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
    Quadrilateral3D4& operator=( Quadrilateral3D4<TOtherPointType> const & rOther )
    {
        BaseType::operator=( rOther );
        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    typename BaseType::Pointer Create( PointsArrayType const& ThisPoints ) const override
    {
        return typename BaseType::Pointer( new Quadrilateral3D4( ThisPoints ) );
    }


    // Geometry< Point<3> >::Pointer  Clone() const override
    // {
    //     Geometry< Point<3> >::PointsArrayType NewPoints;

    //     //making a copy of the nodes TO POINTS (not Nodes!!!)
    //     for ( IndexType i = 0 ; i < this->size() ; i++ )
    //     {
    //         NewPoints.push_back(Kratos::make_shared< Point<3> >(( *this )[i]));
    //     }

    //     //creating a geometry with the new points
    //     Geometry< Point<3> >::Pointer p_clone( new Quadrilateral3D4< Point<3> >( NewPoints ) );

    //     return p_clone;
    // }

    /**
     * returns the local coordinates of all nodes of the current geometry
     * @param rResult a Matrix object that will be overwritten by the result
     * @return the local coordinates of all nodes
     */
    Matrix& PointsLocalCoordinates( Matrix& rResult ) const override
    {
        rResult.resize( 4, 2, false );
        noalias( rResult ) = ZeroMatrix( 4, 2 );
        rResult( 0, 0 ) = -1.0;
        rResult( 0, 1 ) = -1.0;
        rResult( 1, 0 ) =  1.0;
        rResult( 1, 1 ) = -1.0;
        rResult( 2, 0 ) =  1.0;
        rResult( 2, 1 ) =  1.0;
        rResult( 3, 0 ) = -1.0;
        rResult( 3, 1 ) =  1.0;
        return rResult;
    }

    //lumping factors for the calculation of the lumped mass matrix
    Vector& LumpingFactors( Vector& rResult ) const override
    {
        if(rResult.size() != 4) rResult.resize( 4, false );
        std::fill( rResult.begin(), rResult.end(), 1.00 / 4.00 );
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
     *
     * @return double value contains length or Characteristic
     * length
     * @see Area()
     * @see Volume()
     * @see DomainSize()
     */
    /**
     * :TODO: could be replaced by something more suitable
     * (comment by janosch)
     */
    double Length() const override
    {
        return std::sqrt( Area() );
    }

    /** This method calculates and returns area or surface area of
     * this geometry depending to it's dimension. For one dimensional
     * geometry it returns zero, for two dimensional it gives area
     * and for three dimensional geometries it gives surface area.
     *
     * @return double value contains area or surface
     * area.
     * @see Length()
     * @see Volume()
     * @see DomainSize()
     */
    /**
     * :TODO: could be replaced by something more suitable
     * (comment by janosch)
     */
    double Area() const override
    {
        // Finite element way
        Vector temp;
        DeterminantOfJacobian( temp, GeometryData::GI_GAUSS_3 );
        const IntegrationPointsArrayType& integration_points = this->IntegrationPoints( GeometryData::GI_GAUSS_3 );
        double area = 0.0;

        for ( unsigned int i = 0; i < integration_points.size(); i++ )
        {
           area += temp[i] * integration_points[i].Weight();
        }

        return area;

//         // 24/01/2014 - Massimo Petracca
//         // the following procedure calculates the area of a general
//         // quadrilateral (flat or warped) using the parametric representation
//         // of ruled hyperbolic paraboloid surface.
//         // the integration of the normal is then performed with a 2x2 gauss quadrature
//         // in the U-V domain [0,1].
//         // results explicitly written after symbolic calculation.
//
//         const TPointType& p1 = this->Points()[0];
//         const TPointType& p2 = this->Points()[1];
//         const TPointType& p3 = this->Points()[2];
//         const TPointType& p4 = this->Points()[3];
//
//         const double& p1x = p1.X();
//         const double& p1y = p1.Y();
//         const double& p1z = p1.Z();
//
//         const double& p2x = p2.X();
//         const double& p2y = p2.Y();
//         const double& p2z = p2.Z();
//
//         const double& p3x = p3.X();
//         const double& p3y = p3.Y();
//         const double& p3z = p3.Z();
//
//         const double& p4x = p4.X();
//         const double& p4y = p4.Y();
//         const double& p4z = p4.Z();
//
//         const double pos = 0.5 + 0.5 / std::sqrt(3.0);
//         const double w = 0.25;
//
//         const double C1  = pos*(p1z - p2z + p3z - p4z);
//         const double C2  = pos*(p1y - p2y + p3y - p4y);
//         const double C3  = pos*(p1x - p2x + p3x - p4x);
//         const double C4  = C1 - p1z + p2z;
//         const double C5  = C1 + p1z - p2z;
//         const double C6  = C2 + p1y - p2y;
//         const double C7  = C2 - p1y + p2y;
//         const double C8  = C3 - p1x + p2x;
//         const double C9  = C3 + p1x - p2x;
//         const double C10 = C1 - p1z + p4z;
//         const double C11 = C2 - p1y + p4y;
//         const double C12 = C3 - p1x + p4x;
//         const double C13 = C1 + p1z - p4z;
//         const double C14 = C2 + p1y - p4y;
//         const double C15 = C3 + p1x - p4x;
//
//         return w * (
//             std::sqrt( std::pow(C4*C11 - C7*C10, 2) + std::pow(C4*C12 - C8*C10, 2) + std::pow(C7*C12 - C8*C11, 2)) +
//             std::sqrt( std::pow(C5*C11 - C6*C10, 2) + std::pow(C5*C12 - C9*C10, 2) + std::pow(C6*C12 - C9*C11, 2)) +
//             std::sqrt( std::pow(C4*C14 - C7*C13, 2) + std::pow(C4*C15 - C8*C13, 2) + std::pow(C7*C15 - C8*C14, 2)) +
//             std::sqrt( std::pow(C5*C14 - C6*C13, 2) + std::pow(C5*C15 - C9*C13, 2) + std::pow(C6*C15 - C9*C14, 2))
//             );
    }

    /** This method calculates and returns length, area or volume of
     * this geometry depending to it's dimension. For one dimensional
     * geometry it returns its length, for two dimensional it gives area
     * and for three dimensional geometries it gives its volume.
     *
     * @return double value contains length, area or volume.
     * @see Length()
     * @see Area()
     * @see Volume()
     */
    /**
     * :TODO: could be replaced by something more suitable
     * (comment by janosch)
     */
    double DomainSize() const override
    {
        return Area();
    }


    double Volume() const override
    {
        return Area();
    }

    /**
     * Returns whether given arbitrary point is inside the Geometry and the respective
     * local point for the given global point
     * @param rPoint The point to be checked if is inside o note in global coordinates
     * @param rResult The local coordinates of the point
     * @param Tolerance The  tolerance that will be considered to check if the point is inside or not
     * @return True if the point is inside, false otherwise
     */
    virtual bool IsInside(
        const CoordinatesArrayType& rPoint,
        CoordinatesArrayType& rResult,
        const double Tolerance = std::numeric_limits<double>::epsilon()
        ) override
    {
        PointLocalCoordinatesImplementation( rResult, rPoint, true );

        if ( std::abs(rResult[0]) <= (1.0+Tolerance) )
        {
            if ( std::abs(rResult[1]) <= (1.0+Tolerance) )
            {
                return true;
            }
        }

        return false;
    }

    /**
     * Returns the local coordinates of a given arbitrary point
     * @param rResult The vector containing the local coordinates of the point
     * @param rPoint The point in global coordinates
     * @return The vector containing the local coordinates of the point
     */
    CoordinatesArrayType& PointLocalCoordinates(
        CoordinatesArrayType& rResult,
        const CoordinatesArrayType& rPoint
        ) override
    {
        return PointLocalCoordinatesImplementation(rResult, rPoint);
    }

    ///@}
    ///@name Jacobian
    ///@{
    /**
     * TODO: implemented but not yet tested
     */
    /**
     * Jacobians for given method.
     * This method calculates jacobians matrices in all
     * integrations points of given integration method.
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
    JacobiansType& Jacobian(
        JacobiansType& rResult,
        IntegrationMethod ThisMethod
        ) const override
    {
        // Getting derivatives of shape functions
        const ShapeFunctionsGradientsType& shape_functions_gradients =
        msGeometryData.ShapeFunctionsLocalGradients( ThisMethod );
        // Getting values of shape functions
        Matrix shape_functions_values =
        CalculateShapeFunctionsIntegrationPointsValues( ThisMethod );

        if ( rResult.size() != this->IntegrationPointsNumber( ThisMethod ) )
        {
            // KLUDGE: While there is a bug in ublas
            // vector resize, I have to put this beside resizing!!
            JacobiansType temp( this->IntegrationPointsNumber( ThisMethod ) );
            rResult.swap( temp );
        }

        // Loop over all integration points
        for ( unsigned int pnt = 0; pnt < this->IntegrationPointsNumber( ThisMethod ); pnt++ )
        {
            // Defining single jacobian matrix
            Matrix jacobian = ZeroMatrix( 3, 2 );

            // Loop over all nodes
            for ( unsigned int i = 0; i < this->PointsNumber(); i++ )
            {
                jacobian( 0, 0 ) +=
                    ( this->GetPoint( i ).X() ) * ( shape_functions_gradients[pnt]( i, 0 ) );
                jacobian( 0, 1 ) +=
                    ( this->GetPoint( i ).X() ) * ( shape_functions_gradients[pnt]( i, 1 ) );
                jacobian( 1, 0 ) +=
                    ( this->GetPoint( i ).Y() ) * ( shape_functions_gradients[pnt]( i, 0 ) );
                jacobian( 1, 1 ) +=
                    ( this->GetPoint( i ).Y() ) * ( shape_functions_gradients[pnt]( i, 1 ) );
                jacobian( 2, 0 ) +=
                    ( this->GetPoint( i ).Z() ) * ( shape_functions_gradients[pnt]( i, 0 ) );
                jacobian( 2, 1 ) +=
                    ( this->GetPoint( i ).Z() ) * ( shape_functions_gradients[pnt]( i, 1 ) );
            }

            rResult[pnt] = jacobian;
        } //end of loop over all integration points

        return rResult;
    }

    /**
     * Jacobians for given method.
     * This method calculates jacobians matrices in all
     * integrations points of given integration method.
     *
     * @param ThisMethod integration method which jacobians has to
     * be calculated in its integration points.
     *
     * @return JacobiansType a Vector of jacobian
     * matrices \f$ J_i \f$ where \f$ i=1,2,...,n \f$ is the integration
     * point index of given integration method.
     *
     * @param DeltaPosition Matrix with the nodes position increment which describes
     * the configuration where the jacobian has to be calculated.
     *
     * @see DeterminantOfJacobian
     * @see InverseOfJacobian
     */
    JacobiansType& Jacobian(
        JacobiansType& rResult,
        IntegrationMethod ThisMethod,
        Matrix & DeltaPosition
        ) const override
    {
        // Getting derivatives of shape functions
        const ShapeFunctionsGradientsType& shape_functions_gradients =
        msGeometryData.ShapeFunctionsLocalGradients( ThisMethod );
        // Getting values of shape functions
        Matrix shape_functions_values = CalculateShapeFunctionsIntegrationPointsValues( ThisMethod );

        if ( rResult.size() != this->IntegrationPointsNumber( ThisMethod ) )
        {
            // KLUDGE: While there is a bug in ublas
            // vector resize, I have to put this beside resizing!!
            JacobiansType temp( this->IntegrationPointsNumber( ThisMethod ) );
            rResult.swap( temp );
        }

        //loop over all integration points
        for ( unsigned int pnt = 0; pnt < this->IntegrationPointsNumber( ThisMethod ); pnt++ )
        {
            //defining single jacobian matrix
            Matrix jacobian = ZeroMatrix( 3, 2 );
            //loop over all nodes

            for ( unsigned int i = 0; i < this->PointsNumber(); i++ )
            {
                jacobian( 0, 0 ) +=
                    ( this->GetPoint( i ).X() - DeltaPosition(i,0) ) * ( shape_functions_gradients[pnt]( i, 0 ) );
                jacobian( 0, 1 ) +=
                    ( this->GetPoint( i ).X() - DeltaPosition(i,0) ) * ( shape_functions_gradients[pnt]( i, 1 ) );
                jacobian( 1, 0 ) +=
                    ( this->GetPoint( i ).Y() - DeltaPosition(i,1) ) * ( shape_functions_gradients[pnt]( i, 0 ) );
                jacobian( 1, 1 ) +=
                    ( this->GetPoint( i ).Y() - DeltaPosition(i,1) ) * ( shape_functions_gradients[pnt]( i, 1 ) );
                jacobian( 2, 0 ) +=
                    ( this->GetPoint( i ).Z() - DeltaPosition(i,2) ) * ( shape_functions_gradients[pnt]( i, 0 ) );
                jacobian( 2, 1 ) +=
                    ( this->GetPoint( i ).Z() - DeltaPosition(i,2) ) * ( shape_functions_gradients[pnt]( i, 1 ) );
            }

            rResult[pnt] = jacobian;
        }//end of loop over all integration points

        return rResult;
    }
    /**
     * TODO: implemented but not yet tested
     */
    /**
     * Jacobian in specific integration point of given integration
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
     *
     * @see DeterminantOfJacobian
     * @see InverseOfJacobian
     */
    Matrix& Jacobian(
        Matrix& rResult,
        IndexType IntegrationPointIndex,
        IntegrationMethod ThisMethod
        ) const override
    {
        // Setting up size of jacobian matrix
        if (rResult.size1() != 3 || rResult.size2() != 2 )
            rResult.resize( 3, 2, false );
        noalias(rResult) = ZeroMatrix(3, 2);
        // Derivatives of shape functions
        Matrix shape_functions_gradients = msGeometryData.ShapeFunctionLocalGradient(IntegrationPointIndex, ThisMethod );

        //Elements of jacobian matrix (e.g. J(1,1) = dX1/dXi1)
        //loop over all nodes
        for ( unsigned int i = 0; i < this->PointsNumber(); i++ )
        {
            rResult( 0, 0 ) +=
                ( this->GetPoint( i ).X() ) * ( shape_functions_gradients( i, 0 ) );
            rResult( 0, 1 ) +=
                ( this->GetPoint( i ).X() ) * ( shape_functions_gradients( i, 1 ) );
            rResult( 1, 0 ) +=
                ( this->GetPoint( i ).Y() ) * ( shape_functions_gradients( i, 0 ) );
            rResult( 1, 1 ) +=
                ( this->GetPoint( i ).Y() ) * ( shape_functions_gradients( i, 1 ) );
            rResult( 2, 0 ) +=
                ( this->GetPoint( i ).Z() ) * ( shape_functions_gradients( i, 0 ) );
            rResult( 2, 1 ) +=
                ( this->GetPoint( i ).Z() ) * ( shape_functions_gradients( i, 1 ) );
        }

        return rResult;
    }

    /**
     * TODO: implemented but not yet tested
     */
    /**
    * Jacobian in given point. This method calculate jacobian
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
        if (rResult.size1() != 3 || rResult.size2() != 2 )
            rResult.resize( 3, 2, false );
        noalias(rResult) = ZeroMatrix(3, 2);

        // Derivatives of shape functions
        Matrix shape_functions_gradients;
        shape_functions_gradients = ShapeFunctionsLocalGradients(shape_functions_gradients, rPoint );
        //Elements of jacobian matrix (e.g. J(1,1) = dX1/dXi1)

        //loop over all nodes
        for ( unsigned int i = 0; i < this->PointsNumber(); i++ )
        {
            rResult( 0, 0 ) += ( this->GetPoint( i ).X() ) * ( shape_functions_gradients( i, 0 ) );
            rResult( 0, 1 ) += ( this->GetPoint( i ).X() ) * ( shape_functions_gradients( i, 1 ) );
            rResult( 1, 0 ) += ( this->GetPoint( i ).Y() ) * ( shape_functions_gradients( i, 0 ) );
            rResult( 1, 1 ) += ( this->GetPoint( i ).Y() ) * ( shape_functions_gradients( i, 1 ) );
            rResult( 2, 0 ) += ( this->GetPoint( i ).Z() ) * ( shape_functions_gradients( i, 0 ) );
            rResult( 2, 1 ) += ( this->GetPoint( i ).Z() ) * ( shape_functions_gradients( i, 1 ) );
        }

        return rResult;
    }

    /**
     * TODO: implemented but not yet tested
     */
    /**
     * Determinant of jacobians for given integration method.
     * This method calculates determinant of jacobian in all
     * integrations points of given integration method.
     *
     * @return Vector of double which is vector of determinants of
     * jacobians \f$ |J|_i \f$ where \f$ i=1,2,...,n \f$ is the
     * integration point index of given integration method.
     *
     * @see Jacobian
     * @see InverseOfJacobian
     */
    Vector& DeterminantOfJacobian(
        Vector& rResult,
        IntegrationMethod ThisMethod
        ) const override
    {
        const unsigned int integration_points_number = msGeometryData.IntegrationPointsNumber( ThisMethod );
        if(rResult.size() != integration_points_number)
        {
            rResult.resize(integration_points_number,false);
        }

        JacobiansType jacobian;
        this->Jacobian( jacobian, ThisMethod);

        for ( unsigned int pnt = 0; pnt < integration_points_number; pnt++ )
        {
            const double det_j = std::pow(jacobian[pnt](0,1),2)*(std::pow(jacobian[pnt](1,0),2) + std::pow(jacobian[pnt](2,0),2)) + std::pow(jacobian[pnt](1,1)*jacobian[pnt](2,0) - jacobian[pnt](1,0)*jacobian[pnt](2,1),2) - 2.0*jacobian[pnt](0,0)*jacobian[pnt](0,1)*(jacobian[pnt](1,0)*jacobian[pnt](1,1) + jacobian[pnt](2,0)*jacobian[pnt](2,1)) + std::pow(jacobian[pnt](0,0),2)*(std::pow(jacobian[pnt](1,1),2) + std::pow(jacobian[pnt](2,1),2));

            if (det_j < 0.0) KRATOS_ERROR << "WARNING::NEGATIVE VALUE: NOT POSSIBLE TO EVALUATE THE JACOBIAN DETERMINANT" << std::endl;

            rResult[pnt] = std::sqrt(det_j);
        }

        return rResult;
    }

    /**
     * TODO: implemented but not yet tested
     */
    /**
     * Determinant of jacobian in specific integration point of
     * given integration method. This method calculate determinant
     * of jacobian in given integration point of given integration
     * method.
     *
     * @param IntegrationPointIndex index of integration point which jacobians has to
     * be calculated in it.
     *
     * @param IntegrationPointIndex index of integration point
     * which determinant of jacobians has to be calculated in it.
     *
     * @return Determinamt of jacobian matrix \f$ |J|_i \f$ where \f$
     * i \f$ is the given integration point index of given
     * integration method.
     *
     * @see Jacobian
     * @see InverseOfJacobian
     */
    double DeterminantOfJacobian(
        IndexType IntegrationPointIndex,
        IntegrationMethod ThisMethod
        ) const override
    {
        Matrix jacobian( 3, 2 );

        this->Jacobian( jacobian, IntegrationPointIndex, ThisMethod);

        const double det_j = std::pow(jacobian(0,1),2)*(std::pow(jacobian(1,0),2) + std::pow(jacobian(2,0),2)) + std::pow(jacobian(1,1)*jacobian(2,0) - jacobian(1,0)*jacobian(2,1),2) - 2.0*jacobian(0,0)*jacobian(0,1)*(jacobian(1,0)*jacobian(1,1) + jacobian(2,0)*jacobian(2,1)) + std::pow(jacobian(0,0),2)*(std::pow(jacobian(1,1),2) + std::pow(jacobian(2,1),2));

        if (det_j < 0.0) KRATOS_ERROR << "WARNING::NEGATIVE VALUE: NOT POSSIBLE TO EVALUATE THE JACOBIAN DETERMINANT" << std::endl;

        return std::sqrt(det_j);
    }

    /**
     * TODO: implemented but not yet tested
     */
    /**
     * Determinant of jacobian in given point.
     * This method calculate determinant of jacobian
     * matrix in given point.
     * @param rPoint point which determinant of jacobians has to
     * be calculated in it.
     *
     * @return Determinamt of jacobian matrix \f$ |J| \f$ in given
     * point.
     *
     * @see DeterminantOfJacobian
     * @see InverseOfJacobian
     *
     * KLUDGE: PointType needed for proper functionality
     * KLUDGE: works only with explicitly generated Matrix object
     */
    /**
     * :TODO: needs to be changed to Point<3> again. As PointType can
     * be a Node with unique ID or an IntegrationPoint or any arbitrary
     * point in space this needs to be reviewed
     * (comment by janosch)
     */
    double DeterminantOfJacobian( const CoordinatesArrayType& rPoint ) const override
    {
        Matrix jacobian( 3, 2 );

        this->Jacobian( jacobian, rPoint);

        const double det_j = std::pow(jacobian(0,1),2)*(std::pow(jacobian(1,0),2) + std::pow(jacobian(2,0),2)) + std::pow(jacobian(1,1)*jacobian(2,0) - jacobian(1,0)*jacobian(2,1),2) - 2.0*jacobian(0,0)*jacobian(0,1)*(jacobian(1,0)*jacobian(1,1) + jacobian(2,0)*jacobian(2,1)) + std::pow(jacobian(0,0),2)*(std::pow(jacobian(1,1),2) + std::pow(jacobian(2,1),2));

        if (det_j < 0.0) KRATOS_ERROR << "WARNING::NEGATIVE VALUE: NOT POSSIBLE TO EVALUATE THE JACOBIAN DETERMINANT" << std::endl;

        return std::sqrt(det_j);
    }

    /**
     * TODO: implemented but not yet tested
     */
    /**
     * Inverse of jacobians for given integration method.
     * This method calculates inverse of jacobians matrices
     * in all integrations points of
     * given integration method.
     *
     * @param ThisMethod integration method which inverse of jacobians has to
     * be calculated in its integration points.
     *
     * @return Inverse of jacobian
     * matrices \f$ J^{-1}_i \f$ where \f$ i=1,2,...,n \f$ is the integration
     * point index of given integration method.
     *
     * @see Jacobian
     * @see DeterminantOfJacobian
     *
     * KLUDGE: works only with explicitly generated Matrix object
     */
    JacobiansType& InverseOfJacobian( JacobiansType& rResult,
            IntegrationMethod ThisMethod ) const override
    {
        KRATOS_ERROR << "Quadrilateral3D4::InverseOfJacobian" << "Jacobian is not square" << std::endl;
        return rResult;
    }

    /**
     * TODO: implemented but not yet tested
     */
    /**
     * Inverse of jacobian in specific integration point of given integration
     * method. This method calculate Inverse of jacobian matrix in given
     * integration point of given integration method.
     *
     * @param IntegrationPointIndex index of integration point
     * which inverse of jacobians has to
     * be calculated in it.
     * @param ThisMethod integration method which inverse of jacobians has to
     * be calculated in its integration points.
     *
     * @return Inverse of jacobian matrix \f$ J^{-1}_i \f$ where \f$
     * i \f$ is the given integration point index of given
     * integration method.
     *
     * @see Jacobian
     * @see DeterminantOfJacobian
     *
     * KLUDGE: works only with explicitly generated Matrix object
     */
    Matrix& InverseOfJacobian( Matrix& rResult,
                                       IndexType IntegrationPointIndex,
                                       IntegrationMethod ThisMethod ) const override
    {
        KRATOS_ERROR << "Quadrilateral3D4::InverseOfJacobian" << "Jacobian is not square" << std::endl;
        return rResult;
    }

    /**
     * TODO: implemented but not yet tested
     */
    /**
     * Inverse of jacobian in given point.
     * This method calculates inverse of jacobian
     * matrix in given point.
     * @param rPoint point which inverse of jacobians has to
     * be calculated in it.
     * @return Inverse of jacobian matrix \f$ J^{-1} \f$ in given point.
     *
     * @see DeterminantOfJacobian
     * @see InverseOfJacobian
     *
     * KLUDGE: works only with explicitly generated Matrix object
     */
    Matrix& InverseOfJacobian( Matrix& rResult,
                                       const CoordinatesArrayType& rPoint ) const override
    {
        KRATOS_ERROR << "Quadrilateral3D4::InverseOfJacobian" << "Jacobian is not square" << std::endl;
        return rResult;
    }

    /** This method gives you number of all edges of this
    geometry. This method will gives you number of all the edges
    with one dimension less than this geometry. for example a
    triangle would return three or a tetrahedral would return
    four but won't return nine related to its six edge lines.

    @return SizeType containes number of this geometry edges.
    @see Edges()
    @see Edge()
    */
    SizeType EdgesNumber() const override
    {
        return 4;
    }


    /** FacesNumber
    @return SizeType containes number of this geometry edges/faces.
    */
    SizeType FacesNumber() const override
    {
        return EdgesNumber();
    }

    /** This method gives you all edges of this geometry. This
    method will gives you all the edges with one dimension less
    than this geometry. for example a triangle would return
    three lines as its edges or a tetrahedral would return four
    triangle as its edges but won't return its six edge
    lines by this method.

    @return GeometriesArrayType containes this geometry edges.
    @see EdgesNumber()
    @see Edge()
    */
    GeometriesArrayType Edges( void ) override
    {
        GeometriesArrayType edges = GeometriesArrayType();
        typedef typename Geometry<TPointType>::Pointer EdgePointerType;
        edges.push_back( EdgePointerType( new EdgeType( this->pGetPoint( 0 ), this->pGetPoint( 1 ) ) ) );
        edges.push_back( EdgePointerType( new EdgeType( this->pGetPoint( 1 ), this->pGetPoint( 2 ) ) ) );
        edges.push_back( EdgePointerType( new EdgeType( this->pGetPoint( 2 ), this->pGetPoint( 3 ) ) ) );
        edges.push_back( EdgePointerType( new EdgeType( this->pGetPoint( 3 ), this->pGetPoint( 0 ) ) ) );
        return edges;
    }

    //Connectivities of faces required
    void NumberNodesInFaces (DenseVector<unsigned int>& NumberNodesInFaces) const override
    {
        if(NumberNodesInFaces.size() != 4 )
            NumberNodesInFaces.resize(4,false);

        NumberNodesInFaces[0]=2;
        NumberNodesInFaces[1]=2;
        NumberNodesInFaces[2]=2;
	NumberNodesInFaces[3]=2;

    }

    void NodesInFaces (DenseMatrix<unsigned int>& NodesInFaces) const override
    {
        if(NodesInFaces.size1() != 3 || NodesInFaces.size2() != 4)
            NodesInFaces.resize(3,4,false);
        //face 1
        NodesInFaces(0,0)=0;//contrary node to the face
        NodesInFaces(1,0)=2;
        NodesInFaces(2,0)=3;
        //face 2
        NodesInFaces(0,1)=1;//contrary node to the face
        NodesInFaces(1,1)=3;
        NodesInFaces(2,1)=0;
        //face 3
        NodesInFaces(0,2)=2;//contrary node to the face
        NodesInFaces(1,2)=0;
        NodesInFaces(2,2)=1;
        //face 4
        NodesInFaces(0,3)=3;//contrary node to the face
        NodesInFaces(1,3)=1;
        NodesInFaces(2,3)=2;
    }

    /** This method checks if an axis-aliged bounding box (AABB)
    intersects the quadrilateral

    @return bool if the quadrilateral overlaps the box
    @param rLowPoint first corner of the box
    @param rHighPoint second corner of the box
    @see Triangle3D3::HasIntersection
    */
    bool HasIntersection( const Point& rLowPoint, const Point& rHighPoint ) override
    {
        Triangle3D3<PointType> triangle_0 (this->pGetPoint( 0 ),
                                           this->pGetPoint( 1 ),
                                           this->pGetPoint( 2 )
        );
        Triangle3D3<PointType> triangle_1 (this->pGetPoint( 2 ),
                                           this->pGetPoint( 3 ),
                                           this->pGetPoint( 0 )
        );

        if      ( triangle_0.HasIntersection(rLowPoint, rHighPoint) ) return true;
        else if ( triangle_1.HasIntersection(rLowPoint, rHighPoint) ) return true;
        else return false;
    }


    /**
     * Returns all faces of the current geometry.
     * This is only implemented for 3D geometries, since 2D geometries
     * only have edges but no faces
     * @see EdgesNumber
     * @see Edges
     * @see FacesNumber
    */
    GeometriesArrayType Faces( void ) override
    {
        return GeometriesArrayType();
    }

    ///@}
    ///@name Shape Function
    ///@{

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
        switch ( ShapeFunctionIndex )
        {
        case 0:
            return( 0.25*( 1.0 - rPoint[0] )*( 1.0 - rPoint[1] ) );
        case 1:
            return( 0.25*( 1.0 + rPoint[0] )*( 1.0 - rPoint[1] ) );
        case 2:
            return( 0.25*( 1.0 + rPoint[0] )*( 1.0 + rPoint[1] ) );
        case 3:
            return( 0.25*( 1.0 - rPoint[0] )*( 1.0 + rPoint[1] ) );
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
        if(rResult.size() != 4) rResult.resize(4,false);
        rResult[0] =  0.25*( 1.0 - rCoordinates[0] )*( 1.0 - rCoordinates[1] );
        rResult[1] =  0.25*( 1.0 + rCoordinates[0] )*( 1.0 - rCoordinates[1] );
        rResult[2] =  0.25*( 1.0 + rCoordinates[0] )*( 1.0 + rCoordinates[1] );
        rResult[3] =  0.25*( 1.0 - rCoordinates[0] )*( 1.0 + rCoordinates[1] );

        return rResult;
    }

    /**
     * Calculates the Gradients of the shape functions.
     * Calculates the gradients of the shape functions with regard to
     * the global coordinates in all
     * integration points (\f$ \frac{\partial N^i}{\partial X_j} \f$)
     *
     * @param rResult a container which takes the calculated gradients
     * @param ThisMethod the given IntegrationMethod
     *
     * @return the gradients of all shape functions with regard to the global coordinates
     * KLUDGE: method call only works with explicit JacobiansType rather than creating
     * JacobiansType within argument list
    */
    ShapeFunctionsGradientsType& ShapeFunctionsIntegrationPointsGradients(
        ShapeFunctionsGradientsType& rResult,
        IntegrationMethod ThisMethod ) const override
    {
        const unsigned int integration_points_number = msGeometryData.IntegrationPointsNumber( ThisMethod );

        if ( integration_points_number == 0 )
        {
            KRATOS_ERROR << "This integration method is not supported" << *this << std::endl;
        }

        if ( rResult.size() != integration_points_number )
        {
            // KLUDGE: While there is a bug in ublas
            // vector resize, I have to put this beside resizing!!
            ShapeFunctionsGradientsType temp( integration_points_number );
            rResult.swap( temp );
        }

        // Calculating the local gradients
        const ShapeFunctionsGradientsType& shape_functions_local_gradient = msGeometryData.ShapeFunctionsLocalGradients( ThisMethod );

        //getting the inverse jacobian matrices
        JacobiansType temp( integration_points_number );

        JacobiansType invJ = InverseOfJacobian( temp, ThisMethod );

        //loop over all integration points
        for ( unsigned int pnt = 0; pnt < integration_points_number; pnt++ )
        {
            rResult[pnt].resize( 4, 2, false );

            for ( int i = 0; i < 4; i++ )
            {
                for ( int j = 0; j < 2; j++ )
                {
                    rResult[pnt]( i, j ) =
                        ( shape_functions_local_gradient[pnt]( i, 0 ) * invJ[pnt]( j, 0 ) )
                        + ( shape_functions_local_gradient[pnt]( i, 1 ) * invJ[pnt]( j, 1 ) );
                }
            }
        }//end of loop over integration points

        return rResult;
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
        return "2 dimensional quadrilateral with four nodes in 3D space";
    }

    /**
     * Print information about this object.
     * @param rOStream Stream to print into it.
     * @see PrintData()
     * @see Info()
     */
    void PrintInfo( std::ostream& rOStream ) const override
    {
        rOStream << Info();
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
        BaseType::PrintData( rOStream );
        std::cout << std::endl;
        Matrix jacobian;
        Jacobian( jacobian, PointType() );
        rOStream << "    Jacobian in the origin\t : " << jacobian;
    }

    /**
     * Calculates the local gradients for all integration points for
     * given integration method
     */
    virtual ShapeFunctionsGradientsType ShapeFunctionsLocalGradients(
        IntegrationMethod ThisMethod )
    {
        const ShapeFunctionsGradientsType& shape_function_local_gradient
        = msGeometryData.ShapeFunctionsLocalGradients( ThisMethod );
        const int& integration_points_number
        = msGeometryData.IntegrationPointsNumber( ThisMethod );
        ShapeFunctionsGradientsType Result( integration_points_number );

        for ( int pnt = 0; pnt < integration_points_number; pnt++ )
        {
            Result[pnt] = shape_function_local_gradient[pnt];
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
        const ShapeFunctionsGradientsType& shape_function_local_gradient
        = msGeometryData.ShapeFunctionsLocalGradients( ThisMethod );
        const int integration_points_number
        = msGeometryData.IntegrationPointsNumber( ThisMethod );
        ShapeFunctionsGradientsType Result( integration_points_number );

        for ( int pnt = 0; pnt < integration_points_number; pnt++ )
        {
            Result[pnt] = shape_function_local_gradient[pnt];
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
    Matrix& ShapeFunctionsLocalGradients(
        Matrix& rResult,
        const CoordinatesArrayType& rPoint
        ) const override
    {
        rResult.resize( 4, 2, false );
        noalias( rResult ) = ZeroMatrix( 4, 2 );
        rResult( 0, 0 ) = -0.25 * ( 1.0 - rPoint[1] );
        rResult( 0, 1 ) = -0.25 * ( 1.0 - rPoint[0] );
        rResult( 1, 0 ) =  0.25 * ( 1.0 - rPoint[1] );
        rResult( 1, 1 ) = -0.25 * ( 1.0 + rPoint[0] );
        rResult( 2, 0 ) =  0.25 * ( 1.0 + rPoint[1] );
        rResult( 2, 1 ) =  0.25 * ( 1.0 + rPoint[0] );
        rResult( 3, 0 ) = -0.25 * ( 1.0 + rPoint[1] );
        rResult( 3, 1 ) =  0.25 * ( 1.0 - rPoint[0] );
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
    virtual Matrix& ShapeFunctionsGradients(
        Matrix& rResult,
        PointType& rPoint
        )
    {
        rResult.resize( 4, 2, false );
        rResult( 0, 0 ) = -0.25 * ( 1.0 - rPoint.Y() );
        rResult( 0, 1 ) = -0.25 * ( 1.0 - rPoint.X() );
        rResult( 1, 0 ) =  0.25 * ( 1.0 - rPoint.Y() );
        rResult( 1, 1 ) = -0.25 * ( 1.0 + rPoint.X() );
        rResult( 2, 0 ) =  0.25 * ( 1.0 + rPoint.Y() );
        rResult( 2, 1 ) =  0.25 * ( 1.0 + rPoint.X() );
        rResult( 3, 0 ) = -0.25 * ( 1.0 + rPoint.Y() );
        rResult( 3, 1 ) =  0.25 * ( 1.0 - rPoint.X() );
        return rResult;
    }

    /**
     * Returns the second order derivatives of all shape functions
     * in given arbitrary pointers
     * @param rResult a third order tensor which contains the second derivatives
     * @param rPoint the given point the second order derivatives are calculated in
     */
    ShapeFunctionsSecondDerivativesType& ShapeFunctionsSecondDerivatives(
        ShapeFunctionsSecondDerivativesType& rResult,
        const CoordinatesArrayType& rPoint
        ) const override
    {
        if ( rResult.size() != this->PointsNumber() )
        {
            // KLUDGE: While there is a bug in
            // ublas vector resize, I have to put this beside resizing!!
            ShapeFunctionsGradientsType temp( this->PointsNumber() );
            rResult.swap( temp );
        }

        rResult[0].resize( 2, 2, false );
        rResult[1].resize( 2, 2, false );
        rResult[2].resize( 2, 2, false );
        rResult[3].resize( 2, 2, false );

        rResult[0]( 0, 0 ) = 0.0;
        rResult[0]( 0, 1 ) = 0.25;
        rResult[0]( 1, 0 ) = 0.25;
        rResult[0]( 1, 1 ) = 0.0;
        rResult[1]( 0, 0 ) = 0.0;
        rResult[1]( 0, 1 ) = -0.25;
        rResult[1]( 1, 0 ) = -0.25;
        rResult[1]( 1, 1 ) = 0.0;
        rResult[2]( 0, 0 ) = 0.0;
        rResult[2]( 0, 1 ) = 0.25;
        rResult[2]( 1, 0 ) = 0.25;
        rResult[2]( 1, 1 ) = 0.0;
        rResult[3]( 0, 0 ) = 0.0;
        rResult[3]( 0, 1 ) = -0.25;
        rResult[3]( 1, 0 ) = -0.25;
        rResult[3]( 1, 1 ) = 0.0;
        return rResult;
    }

    /**
     * returns the third order derivatives of all shape functions
     * in given arbitrary pointers
     * @param rResult a fourth order tensor which contains the third derivatives
     * @param rPoint the given point the third order derivatives are calculated in
     */
    ShapeFunctionsThirdDerivativesType& ShapeFunctionsThirdDerivatives(
        ShapeFunctionsThirdDerivativesType& rResult,
        const CoordinatesArrayType& rPoint
        ) const override
    {
        if ( rResult.size() != this->PointsNumber() )
        {
            // KLUDGE: While there is a bug in
            // ublas vector resize, I have to put this beside resizing!!
            ShapeFunctionsThirdDerivativesType temp( this->PointsNumber() );
            rResult.swap( temp );
        }

        for ( IndexType i = 0; i < rResult.size(); i++ )
        {
            DenseVector<Matrix> temp( this->PointsNumber() );
            rResult[i].swap( temp );
        }

        rResult[0][0].resize( 2, 2, false );
        rResult[0][1].resize( 2, 2, false );
        rResult[1][0].resize( 2, 2, false );
        rResult[1][1].resize( 2, 2, false );
        rResult[2][0].resize( 2, 2, false );
        rResult[2][1].resize( 2, 2, false );
        rResult[3][0].resize( 2, 2, false );
        rResult[3][1].resize( 2, 2, false );

        for ( int i = 0; i < 4; i++ )
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

    ///@}
    ///@name Friends
    ///@{

    ///@}

protected:
    /**
     * There are no protected members in class Quadrilateral3D4
     */

private:
    ///@name Static Member Variables
    ///@{
    static const GeometryData msGeometryData;


    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save( Serializer& rSerializer ) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, PointsArrayType );
    }

    void load( Serializer& rSerializer ) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, PointsArrayType );
    }

    Quadrilateral3D4(): BaseType( PointsArrayType(), &msGeometryData ) {}

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
     * @brief Returns the local coordinates of a given arbitrary point
     * @param rResult The vector containing the local coordinates of the point
     * @param rPoint The point in global coordinates
     * @param IsInside The flag that checks if we are computing IsInside (is common for seach to have the nodes outside the geometry)
     * @return The vector containing the local coordinates of the point
     */
    CoordinatesArrayType& PointLocalCoordinatesImplementation(
        CoordinatesArrayType& rResult,
        const CoordinatesArrayType& rPoint,
        const bool IsInside = false
        )
    {
        BoundedMatrix<double,3,4> X;
        BoundedMatrix<double,3,2> DN;
        for(IndexType i=0; i<this->size();i++) {
            X(0, i) = this->GetPoint( i ).X();
            X(1, i) = this->GetPoint( i ).Y();
            X(2, i) = this->GetPoint( i ).Z();
        }

        static constexpr double MaxNormPointLocalCoordinates = 300.0;
        static constexpr std::size_t MaxIteratioNumberPointLocalCoordinates = 500;
        static constexpr double MaxTolerancePointLocalCoordinates = 1.0e-8;

        Matrix J = ZeroMatrix( 2, 2 );
        Matrix invJ = ZeroMatrix( 2, 2 );

        // Starting with xi = 0
        rResult = ZeroVector( 3 );
        array_1d<double, 2> DeltaXi( 2, 0.0 );
	const array_1d<double, 3> zero_array(3, 0.0);
        array_1d<double, 3> CurrentGlobalCoords;

        //Newton iteration:
        for ( IndexType k = 0; k < MaxIteratioNumberPointLocalCoordinates; k++ ) {
            noalias(CurrentGlobalCoords) = zero_array;
            this->GlobalCoordinates( CurrentGlobalCoords, rResult );

            noalias( CurrentGlobalCoords ) = rPoint - CurrentGlobalCoords;

            // Derivatives of shape functions
            Matrix shape_functions_gradients;
            shape_functions_gradients = ShapeFunctionsLocalGradients(shape_functions_gradients, rResult );
            noalias(DN) = prod(X,shape_functions_gradients);

            noalias(J) = prod(trans(DN),DN);
            const array_1d<double, 2> res = prod(trans(DN), CurrentGlobalCoords);

            // Deteminant of Jacobian
            const double det_j = J( 0, 0 ) * J( 1, 1 ) - J( 0, 1 ) * J( 1, 0 );

            // Filling matrix
            invJ( 0, 0 ) = ( J( 1, 1 ) ) / ( det_j );
            invJ( 1, 0 ) = -( J( 1, 0 ) ) / ( det_j );
            invJ( 0, 1 ) = -( J( 0, 1 ) ) / ( det_j );
            invJ( 1, 1 ) = ( J( 0, 0 ) ) / ( det_j );

            DeltaXi( 0 ) = invJ( 0, 0 ) * res[0] + invJ( 0, 1 ) * res[1];
            DeltaXi( 1 ) = invJ( 1, 0 ) * res[0] + invJ( 1, 1 ) * res[1];

            rResult[0] += DeltaXi[0];
            rResult[1] += DeltaXi[1];

            if ( norm_2( DeltaXi ) > MaxNormPointLocalCoordinates ) {
                KRATOS_WARNING_IF("Quadrilateral3D4", IsInside == false && k > 0) << "detJ =\t" << det_j << " DeltaX =\t" << DeltaXi << " stopping calculation. Iteration:\t" << k << std::endl;
                break;
            }

            if ( norm_2( DeltaXi ) < MaxTolerancePointLocalCoordinates )
                break;
        }

        return rResult;
    }

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
        IntegrationPointsArrayType integration_points =
            all_integration_points[ThisMethod];
        //number of integration points
        const int integration_points_number = integration_points.size();
        //number of nodes in current geometry
        const int points_number = 4;
        //setting up return matrix
        Matrix shape_function_values( integration_points_number, points_number );
        //loop over all integration points

        for ( int pnt = 0; pnt < integration_points_number; pnt++ )
        {
            shape_function_values( pnt, 0 ) =
                0.25 * ( 1.0 - integration_points[pnt].X() )
                * ( 1.0 - integration_points[pnt].Y() );
            shape_function_values( pnt, 1 ) =
                0.25 * ( 1.0 + integration_points[pnt].X() )
                * ( 1.0 - integration_points[pnt].Y() );
            shape_function_values( pnt, 2 ) =
                0.25 * ( 1.0 + integration_points[pnt].X() )
                * ( 1.0 + integration_points[pnt].Y() );
            shape_function_values( pnt, 3 ) =
                0.25 * ( 1.0 - integration_points[pnt].X() )
                * ( 1.0 + integration_points[pnt].Y() );
        }

        return shape_function_values;
    }

    /**
     * Calculates the local gradients of all shape functions
     * in all integration points.
     * Integration points are expected to be given in local coordinates
     *
     * @param ThisMethod the current integration method
     * @return the vector of the gradients of all shape functions
     * in each integration point
     */
    static ShapeFunctionsGradientsType CalculateShapeFunctionsIntegrationPointsLocalGradients( typename BaseType::IntegrationMethod ThisMethod )
    {
        IntegrationPointsContainerType all_integration_points = AllIntegrationPoints();
        IntegrationPointsArrayType integration_points = all_integration_points[ThisMethod];
        //number of integration points
        const int integration_points_number = integration_points.size();
        ShapeFunctionsGradientsType d_shape_f_values( integration_points_number );
        //initialising container
        //std::fill(d_shape_f_values.begin(), d_shape_f_values.end(), Matrix(4,2));
        //loop over all integration points

        for ( int pnt = 0; pnt < integration_points_number; pnt++ )
        {
            Matrix result( 4, 2 );
            result( 0, 0 ) = -0.25 * ( 1.0 - integration_points[pnt].Y() );
            result( 0, 1 ) = -0.25 * ( 1.0 - integration_points[pnt].X() );
            result( 1, 0 ) =  0.25 * ( 1.0 - integration_points[pnt].Y() );
            result( 1, 1 ) = -0.25 * ( 1.0 + integration_points[pnt].X() );
            result( 2, 0 ) =  0.25 * ( 1.0 + integration_points[pnt].Y() );
            result( 2, 1 ) =  0.25 * ( 1.0 + integration_points[pnt].X() );
            result( 3, 0 ) = -0.25 * ( 1.0 + integration_points[pnt].Y() );
            result( 3, 1 ) =  0.25 * ( 1.0 - integration_points[pnt].X() );
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
                Quadrature < QuadrilateralGaussLegendreIntegrationPoints1,
                2, IntegrationPoint<3> >::GenerateIntegrationPoints(),
                Quadrature < QuadrilateralGaussLegendreIntegrationPoints2,
                2, IntegrationPoint<3> >::GenerateIntegrationPoints(),
                Quadrature < QuadrilateralGaussLegendreIntegrationPoints3,
                2, IntegrationPoint<3> >::GenerateIntegrationPoints(),
                Quadrature < QuadrilateralGaussLegendreIntegrationPoints4,
                2, IntegrationPoint<3> >::GenerateIntegrationPoints(),
                Quadrature < QuadrilateralGaussLegendreIntegrationPoints5,
                2, IntegrationPoint<3> >::GenerateIntegrationPoints(),
                Quadrature < QuadrilateralCollocationIntegrationPoints1,
                2, IntegrationPoint<3> >::GenerateIntegrationPoints(),
                Quadrature < QuadrilateralCollocationIntegrationPoints2,
                2, IntegrationPoint<3> >::GenerateIntegrationPoints(),
                Quadrature < QuadrilateralCollocationIntegrationPoints3,
                2, IntegrationPoint<3> >::GenerateIntegrationPoints(),
                Quadrature < QuadrilateralCollocationIntegrationPoints4,
                2, IntegrationPoint<3> >::GenerateIntegrationPoints(),
                Quadrature < QuadrilateralCollocationIntegrationPoints5,
                2, IntegrationPoint<3> >::GenerateIntegrationPoints()
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
                Quadrilateral3D4<TPointType>::CalculateShapeFunctionsIntegrationPointsValues(
                    GeometryData::GI_GAUSS_1 ),
                Quadrilateral3D4<TPointType>::CalculateShapeFunctionsIntegrationPointsValues(
                    GeometryData::GI_GAUSS_2 ),
                Quadrilateral3D4<TPointType>::CalculateShapeFunctionsIntegrationPointsValues(
                    GeometryData::GI_GAUSS_3 ),
                Quadrilateral3D4<TPointType>::CalculateShapeFunctionsIntegrationPointsValues(
                    GeometryData::GI_GAUSS_4 ),
                Quadrilateral3D4<TPointType>::CalculateShapeFunctionsIntegrationPointsValues(
                    GeometryData::GI_GAUSS_5 ),
                Quadrilateral3D4<TPointType>::CalculateShapeFunctionsIntegrationPointsValues(
                    GeometryData::GI_EXTENDED_GAUSS_1 ),
                Quadrilateral3D4<TPointType>::CalculateShapeFunctionsIntegrationPointsValues(
                    GeometryData::GI_EXTENDED_GAUSS_2 ),
                Quadrilateral3D4<TPointType>::CalculateShapeFunctionsIntegrationPointsValues(
                    GeometryData::GI_EXTENDED_GAUSS_3 ),
                Quadrilateral3D4<TPointType>::CalculateShapeFunctionsIntegrationPointsValues(
                    GeometryData::GI_EXTENDED_GAUSS_4 ),
                Quadrilateral3D4<TPointType>::CalculateShapeFunctionsIntegrationPointsValues(
                    GeometryData::GI_EXTENDED_GAUSS_5 )
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
                Quadrilateral3D4<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients( GeometryData::GI_GAUSS_1 ),
                Quadrilateral3D4<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients( GeometryData::GI_GAUSS_2 ),
                Quadrilateral3D4<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients( GeometryData::GI_GAUSS_3 ),
                Quadrilateral3D4<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients( GeometryData::GI_GAUSS_4 ),
                Quadrilateral3D4<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients( GeometryData::GI_GAUSS_5 ),
                Quadrilateral3D4<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients( GeometryData::GI_EXTENDED_GAUSS_1 ),
                Quadrilateral3D4<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients( GeometryData::GI_EXTENDED_GAUSS_2 ),
                Quadrilateral3D4<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients( GeometryData::GI_EXTENDED_GAUSS_3 ),
                Quadrilateral3D4<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients( GeometryData::GI_EXTENDED_GAUSS_4 ),
                Quadrilateral3D4<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients( GeometryData::GI_EXTENDED_GAUSS_5 )
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

    template<class TOtherPointType> friend class Quadrilateral3D4;

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
    Quadrilateral3D4<TPointType>& rThis );
/**
 * output stream functions
 */
template<class TPointType> inline std::ostream& operator << (
    std::ostream& rOStream,
    const Quadrilateral3D4<TPointType>& rThis )
{
    rThis.PrintInfo( rOStream );
    rOStream << std::endl;
    rThis.PrintData( rOStream );
    return rOStream;
}

///@}

template<class TPointType> const
GeometryData Quadrilateral3D4<TPointType>::msGeometryData(
    2, 3, 2,
    GeometryData::GI_GAUSS_2,
    Quadrilateral3D4<TPointType>::AllIntegrationPoints(),
    Quadrilateral3D4<TPointType>::AllShapeFunctionsValues(),
    AllShapeFunctionsLocalGradients()
);
}// namespace Kratos.

#endif // KRATOS_QUADRILATERAL_3D_4_H_INCLUDED  defined
