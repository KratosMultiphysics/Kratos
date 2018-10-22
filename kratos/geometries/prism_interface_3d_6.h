//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Ignasi de Pouplana
//

#if !defined(KRATOS_PRISM_INTERFACE_3D_6_H_INCLUDED )
#define  KRATOS_PRISM_INTERFACE_3D_6_H_INCLUDED

// System includes

// External includes

// Project includes
#include "geometries/triangle_3d_3.h"
#include "geometries/quadrilateral_3d_4.h"
#include "integration/prism_gauss_lobatto_integration_points.h"

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
 * A six node prism interface geometry. The shape functions are the same as for the
 * prism_3d_6 element, but the jacobian is computed as in the triangle_3d_3 element to
 * to avoid having detJ = 0. Default integration method is Lobatto.
 */

template<class TPointType> class PrismInterface3D6
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
    typedef Triangle3D3<TPointType> FaceType1;
    typedef Quadrilateral3D4<TPointType> FaceType2;

    /**
     * Pointer definition of PrismInterface3D6
     */
    KRATOS_CLASS_POINTER_DEFINITION( PrismInterface3D6 );

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

//     PrismInterface3D6( const PointType& Point1, const PointType& Point2,
//               const PointType& Point3, const PointType& Point4,
//               const PointType& Point5, const PointType& Point6 )
//         : BaseType( PointsArrayType(), &msGeometryData )
//     {
//         this->Points().reserve( 6 );
//         this->Points().push_back( typename PointType::Pointer( new PointType( Point1 ) ) );
//         this->Points().push_back( typename PointType::Pointer( new PointType( Point2 ) ) );
//         this->Points().push_back( typename PointType::Pointer( new PointType( Point3 ) ) );
//         this->Points().push_back( typename PointType::Pointer( new PointType( Point4 ) ) );
//         this->Points().push_back( typename PointType::Pointer( new PointType( Point5 ) ) );
//         this->Points().push_back( typename PointType::Pointer( new PointType( Point6 ) ) );
//     }

    PrismInterface3D6( typename PointType::Pointer pPoint1,
              typename PointType::Pointer pPoint2,
              typename PointType::Pointer pPoint3,
              typename PointType::Pointer pPoint4,
              typename PointType::Pointer pPoint5,
              typename PointType::Pointer pPoint6 )
        : BaseType( PointsArrayType(), &msGeometryData )
    {
        this->Points().reserve( 6 );
        this->Points().push_back( pPoint1 );
        this->Points().push_back( pPoint2 );
        this->Points().push_back( pPoint3 );
        this->Points().push_back( pPoint4 );
        this->Points().push_back( pPoint5 );
        this->Points().push_back( pPoint6 );
    }

    PrismInterface3D6( const PointsArrayType& ThisPoints )
        : BaseType( ThisPoints, &msGeometryData )
    {
        if ( this->PointsNumber() != 6 )
            KRATOS_ERROR << "Invalid points number. Expected 6, given " << this->PointsNumber() << std::endl;
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
    PrismInterface3D6( PrismInterface3D6 const& rOther )
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
    template<class TOtherPointType> PrismInterface3D6( PrismInterface3D6<TOtherPointType> const& rOther )
        : BaseType( rOther )
    {
    }

    /**
     * Destructor. Does nothing!!!
     */
    virtual ~PrismInterface3D6() {}

    GeometryData::KratosGeometryFamily GetGeometryFamily() const override
    {
        return GeometryData::Kratos_Prism;
    }

    GeometryData::KratosGeometryType GetGeometryType() const override
    {
        return GeometryData::Kratos_Prism3D6;
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
    PrismInterface3D6& operator=( const PrismInterface3D6& rOther )
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
    PrismInterface3D6& operator=( PrismInterface3D6<TOtherPointType> const & rOther )
    {
        BaseType::operator=( rOther );
        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    typename BaseType::Pointer Create( PointsArrayType const& ThisPoints ) const override
    {
        return typename BaseType::Pointer( new PrismInterface3D6( ThisPoints ) );
    }

    // Geometry< Point<3> >::Pointer Clone() const override
    // {
    //     Geometry< Point<3> >::PointsArrayType NewPoints;

    //     //making a copy of the nodes TO POINTS (not Nodes!!!)
    //     for ( IndexType i = 0 ; i < this->size() ; i++ )
    //     {
    //             NewPoints.push_back(Kratos::make_shared< Point<3> >(( *this )[i]));
    //     }

    //     //creating a geometry with the new points
    //     Geometry< Point<3> >::Pointer p_clone( new PrismInterface3D6< Point<3> >( NewPoints ) );

    //     return p_clone;
    // }


    /**
     * returns the local coordinates of all nodes of the current geometry
     * @param rResult a Matrix object that will be overwritten by the result
     * @return the local coordinates of all nodes
     */
    Matrix& PointsLocalCoordinates( Matrix& rResult ) const override
    {
        if ( rResult.size1() != 6 || rResult.size2() != 3 )
            rResult.resize( 6, 3, false );

        rResult( 0, 0 ) = 0.0;
        rResult( 0, 1 ) = 0.0;
        rResult( 0, 2 ) = 0.0;

        rResult( 1, 0 ) = 1.0;
        rResult( 1, 1 ) = 0.0;
        rResult( 1, 2 ) = 0.0;

        rResult( 2, 0 ) = 0.0;
        rResult( 2, 1 ) = 1.0;
        rResult( 2, 2 ) = 0.0;

        rResult( 3, 0 ) = 0.0;
        rResult( 3, 1 ) = 0.0;
        rResult( 3, 2 ) = 1.0;

        rResult( 4, 0 ) = 1.0;
        rResult( 4, 1 ) = 0.0;
        rResult( 4, 2 ) = 1.0;

        rResult( 5, 0 ) = 0.0;
        rResult( 5, 1 ) = 1.0;
        rResult( 5, 2 ) = 1.0;

        return rResult;
    }

    //lumping factors for the calculation of the lumped mass matrix
    Vector& LumpingFactors( Vector& rResult ) const override
    {
        rResult.resize( 6, false );
        std::fill( rResult.begin(), rResult.end(), 1.00 / 6.00 );
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
        return std::sqrt(2.0*Area());
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
		array_1d<double, 3> p0 = 0.5 * (BaseType::GetPoint( 0 ) + BaseType::GetPoint( 3 ));
		array_1d<double, 3> p1 = 0.5 * (BaseType::GetPoint( 1 ) + BaseType::GetPoint( 4 ));
		array_1d<double, 3> p2 = 0.5 * (BaseType::GetPoint( 2 ) + BaseType::GetPoint( 5 ));

        //heron formula
        Vector side_a = ( p0 - p1 );
        double a = MathUtils<double>::Norm3( side_a );
        Vector side_b = ( p1 - p2 );
        double b = MathUtils<double>::Norm3( side_b );
        Vector side_c = ( p2 - p0 );
        double c = MathUtils<double>::Norm3( side_c );
        double s = ( a + b + c ) / 2;
        return( sqrt( s*( s - a )*( s - b )*( s - c ) ) );
    }

    double Volume() const override
    {
        return Area();
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
        ) override
    {
        this->PointLocalCoordinates( rResult, rPoint );

        if ( rResult[0] >= (0.0 - Tolerance) )
            if ( rResult[1] >= (0.0 - Tolerance) )
                if ( rResult[2] >= (0.0 - Tolerance) )
                    if( rResult[2] <= (1.0 + Tolerance) )
                        if ( (rResult[0] + rResult[1]) <= (1.0 + Tolerance) )
                            return true;

        return false;
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
        BoundedMatrix<double,3,3> X;
        BoundedMatrix<double,3,2> DN;
        for(unsigned int i=0; i<3;i++) {
            X(0,i ) = 0.5*(this->GetPoint( i ).X() + this->GetPoint( i+3 ).X());
            X(1,i ) = 0.5*(this->GetPoint( i ).Y() + this->GetPoint( i+3 ).Y());
            X(2,i ) = 0.5*(this->GetPoint( i ).Z() + this->GetPoint( i+3 ).Z());
        }

        double tol = 1.0e-8;
        int maxiter = 1000;

        Matrix J = ZeroMatrix( 2, 2 );
        Matrix invJ = ZeroMatrix( 2, 2 );

        //starting with xi = 0
        rResult = ZeroVector( 3 );
        Vector DeltaXi = ZeroVector( 2 );
        array_1d<double,3> CurrentGlobalCoords;


        //Newton iteration:
        for ( int k = 0; k < maxiter; k++ )
        {
            noalias(CurrentGlobalCoords) = ZeroVector( 3 );
            this->GlobalCoordinates( CurrentGlobalCoords, rResult );

            noalias( CurrentGlobalCoords ) = rPoint - CurrentGlobalCoords;


            //derivatives of shape functions
            Matrix shape_functions_gradients;
            shape_functions_gradients = ShapeFunctionsLocalGradients(shape_functions_gradients, rResult );
            noalias(DN) = prod(X,shape_functions_gradients);

            noalias(J) = prod(trans(DN),DN);
            Vector res = prod(trans(DN),CurrentGlobalCoords);

            //deteminant of Jacobian
            const double det_j = J( 0, 0 ) * J( 1, 1 ) - J( 0, 1 ) * J( 1, 0 );

            //filling matrix
            invJ( 0, 0 ) = ( J( 1, 1 ) ) / ( det_j );
            invJ( 1, 0 ) = -( J( 1, 0 ) ) / ( det_j );
            invJ( 0, 1 ) = -( J( 0, 1 ) ) / ( det_j );
            invJ( 1, 1 ) = ( J( 0, 0 ) ) / ( det_j );


            DeltaXi( 0 ) = invJ( 0, 0 ) * res[0] + invJ( 0, 1 ) * res[1];
            DeltaXi( 1 ) = invJ( 1, 0 ) * res[0] + invJ( 1, 1 ) * res[1];

            rResult[0] += DeltaXi[0];
            rResult[1] += DeltaXi[1];
            rResult[2] = 0.0;

            if ( k>0 && norm_2( DeltaXi ) > 30 )
            {
                KRATOS_ERROR << "Computation of local coordinates failed at iteration " << k<< std::endl;
            }

            if ( norm_2( DeltaXi ) < tol )
            {
                break;
            }
        }

        return( rResult );
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
    Matrix& Jacobian( Matrix& rResult,
                      IndexType IntegrationPointIndex,
                      IntegrationMethod ThisMethod ) const override
    {
		array_1d<double, 3> p0 = 0.5 * (BaseType::GetPoint( 0 ) + BaseType::GetPoint( 3 ));
		array_1d<double, 3> p1 = 0.5 * (BaseType::GetPoint( 1 ) + BaseType::GetPoint( 4 ));
		array_1d<double, 3> p2 = 0.5 * (BaseType::GetPoint( 2 ) + BaseType::GetPoint( 5 ));

		rResult.resize( 3, 2, false );
        rResult( 0, 0 ) = -( p0[0] ) + ( p1[0] ); //on the Gauss points (J is constant at each element)
        rResult( 1, 0 ) = -( p0[1] ) + ( p1[1] );
        rResult( 2, 0 ) = -( p0[2] ) + ( p1[2] );
        rResult( 0, 1 ) = -( p0[0] ) + ( p2[0] );
        rResult( 1, 1 ) = -( p0[1] ) + ( p2[1] );
        rResult( 2, 1 ) = -( p0[2] ) + ( p2[2] );
        return rResult;
    }

	/** Jacobian in specific integration point of given integration
    method. This method calculate jacobian matrix in given
    integration point of given integration method.

    @param IntegrationPointIndex index of integration point which jacobians has to
    be calculated in it.

    @param ThisMethod integration method which jacobians has to
    be calculated in its integration points.

    @param DeltaPosition Matrix with the nodes position increment which describes
    the configuration where the jacobian has to be calculated.

    @return Matrix<double> Jacobian matrix \f$ J_i \f$ where \f$
    i \f$ is the given integration point index of given
    integration method.

    @see DeterminantOfJacobian
    @see InverseOfJacobian
    */
    Matrix& Jacobian( Matrix& rResult,
                      IndexType IntegrationPointIndex,
                      IntegrationMethod ThisMethod,
                      Matrix& DeltaPosition ) const override
    {
        array_1d<double, 3> p0 = 0.5 * (BaseType::GetPoint( 0 ) + BaseType::GetPoint( 3 ));
		array_1d<double, 3> p1 = 0.5 * (BaseType::GetPoint( 1 ) + BaseType::GetPoint( 4 ));
		array_1d<double, 3> p2 = 0.5 * (BaseType::GetPoint( 2 ) + BaseType::GetPoint( 5 ));

		Matrix deltaPosMid(3, 3);
		for(unsigned int i = 0; i < 3; i++)
			for(unsigned int j = 0; j < 3; j++)
				deltaPosMid(i, j) = 0.5*( DeltaPosition(i, j) + DeltaPosition(i+3, j) );

        if(rResult.size1() != 3 || rResult.size2() != 2)
			rResult.resize(3, 2, false);

        rResult( 0, 0 ) = -( p0[0] - deltaPosMid(0,0) ) + ( p1[0] - deltaPosMid(1,0) ); //on the Gauss points (J is constant at each element)
        rResult( 1, 0 ) = -( p0[1] - deltaPosMid(0,1) ) + ( p1[1] - deltaPosMid(1,1) );
        rResult( 2, 0 ) = -( p0[2] - deltaPosMid(0,2) ) + ( p1[2] - deltaPosMid(1,2) );
        rResult( 0, 1 ) = -( p0[0] - deltaPosMid(0,0) ) + ( p2[0] - deltaPosMid(2,0) );
        rResult( 1, 1 ) = -( p0[1] - deltaPosMid(0,1) ) + ( p2[1] - deltaPosMid(2,1) );
        rResult( 2, 1 ) = -( p0[2] - deltaPosMid(0,2) ) + ( p2[2] - deltaPosMid(2,2) );

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
        array_1d<double, 3> p0 = 0.5 * (BaseType::GetPoint( 0 ) + BaseType::GetPoint( 3 ));
		array_1d<double, 3> p1 = 0.5 * (BaseType::GetPoint( 1 ) + BaseType::GetPoint( 4 ));
		array_1d<double, 3> p2 = 0.5 * (BaseType::GetPoint( 2 ) + BaseType::GetPoint( 5 ));

		rResult.resize( 3, 2 , false);
        rResult( 0, 0 ) = -( p0[0] ) + ( p1[0] ); //on the Gauss points (J is constant at each element)
        rResult( 1, 0 ) = -( p0[1] ) + ( p1[1] );
        rResult( 2, 0 ) = -( p0[2] ) + ( p1[2] );
        rResult( 0, 1 ) = -( p0[0] ) + ( p2[0] );
        rResult( 1, 1 ) = -( p0[1] ) + ( p2[1] );
        rResult( 2, 1 ) = -( p0[2] ) + ( p2[2] );
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
    Vector& DeterminantOfJacobian( Vector& rResult,
                                   IntegrationMethod ThisMethod ) const override
    {
		array_1d<double, 3> p0 = 0.5 * (BaseType::GetPoint( 0 ) + BaseType::GetPoint( 3 ));
		array_1d<double, 3> p1 = 0.5 * (BaseType::GetPoint( 1 ) + BaseType::GetPoint( 4 ));
		array_1d<double, 3> p2 = 0.5 * (BaseType::GetPoint( 2 ) + BaseType::GetPoint( 5 ));

		array_1d<double, 3> Tangent1( p1 - p0 );
		array_1d<double, 3> Tangent2( p2 - p0 );

		array_1d<double, 3> Normal;
		MathUtils<double>::CrossProduct( Normal, Tangent1, Tangent2);

        double detJ = norm_2(Normal);

        //for all integration points
        if ( rResult.size() != this->IntegrationPointsNumber( ThisMethod ) )
			rResult.resize( this->IntegrationPointsNumber( ThisMethod ), false );

        for ( unsigned int pnt = 0; pnt < this->IntegrationPointsNumber( ThisMethod ); pnt++ )
            rResult[pnt] = detJ;

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
    double DeterminantOfJacobian( IndexType IntegrationPointIndex,
                                  IntegrationMethod ThisMethod ) const override
    {
		array_1d<double, 3> p0 = 0.5 * (BaseType::GetPoint( 0 ) + BaseType::GetPoint( 3 ));
		array_1d<double, 3> p1 = 0.5 * (BaseType::GetPoint( 1 ) + BaseType::GetPoint( 4 ));
		array_1d<double, 3> p2 = 0.5 * (BaseType::GetPoint( 2 ) + BaseType::GetPoint( 5 ));

		array_1d<double, 3> Tangent1( p1 - p0 );
		array_1d<double, 3> Tangent2( p2 - p0 );

		array_1d<double, 3> Normal;
		MathUtils<double>::CrossProduct( Normal, Tangent1, Tangent2);

        return norm_2(Normal);
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
        array_1d<double, 3> p0 = 0.5 * (BaseType::GetPoint( 0 ) + BaseType::GetPoint( 3 ));
		array_1d<double, 3> p1 = 0.5 * (BaseType::GetPoint( 1 ) + BaseType::GetPoint( 4 ));
		array_1d<double, 3> p2 = 0.5 * (BaseType::GetPoint( 2 ) + BaseType::GetPoint( 5 ));

		array_1d<double, 3> Tangent1( p1 - p0 );
		array_1d<double, 3> Tangent2( p2 - p0 );

		array_1d<double, 3> Normal;
		MathUtils<double>::CrossProduct( Normal, Tangent1, Tangent2);

        return norm_2(Normal);
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
        KRATOS_ERROR << "Jacobian is not square" << std::endl;;
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
        KRATOS_ERROR << "Jacobian is not square" << std::endl;
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
        KRATOS_ERROR << "Jacobian is not square" << std::endl;
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
        return 9;
    }

    SizeType FacesNumber() const override
    {
        return 5;
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
        edges.push_back( EdgePointerType( new EdgeType(
                                              this->pGetPoint( 0 ),
                                              this->pGetPoint( 1 ) ) ) );
        edges.push_back( EdgePointerType( new EdgeType(
                                              this->pGetPoint( 1 ),
                                              this->pGetPoint( 2 ) ) ) );
        edges.push_back( EdgePointerType( new EdgeType(
                                              this->pGetPoint( 2 ),
                                              this->pGetPoint( 0 ) ) ) );
        edges.push_back( EdgePointerType( new EdgeType(
                                              this->pGetPoint( 3 ),
                                              this->pGetPoint( 4 ) ) ) );
        edges.push_back( EdgePointerType( new EdgeType(
                                              this->pGetPoint( 4 ),
                                              this->pGetPoint( 5 ) ) ) );
        edges.push_back( EdgePointerType( new EdgeType(
                                              this->pGetPoint( 5 ),
                                              this->pGetPoint( 3 ) ) ) );
        edges.push_back( EdgePointerType( new EdgeType(
                                              this->pGetPoint( 0 ),
                                              this->pGetPoint( 3 ) ) ) );
        edges.push_back( EdgePointerType( new EdgeType(
                                              this->pGetPoint( 1 ),
                                              this->pGetPoint( 4 ) ) ) );
        edges.push_back( EdgePointerType( new EdgeType(
                                              this->pGetPoint( 2 ),
                                              this->pGetPoint( 5 ) ) ) );
        return edges;
    }

    GeometriesArrayType Faces( void ) override
    {
        GeometriesArrayType faces = GeometriesArrayType();
        typedef typename Geometry<TPointType>::Pointer FacePointerType;
        faces.push_back( FacePointerType( new FaceType1(
                                              this->pGetPoint( 0 ),
                                              this->pGetPoint( 2 ),
                                              this->pGetPoint( 1 ) ) ) );
        faces.push_back( FacePointerType( new FaceType1(
                                              this->pGetPoint( 3 ),
                                              this->pGetPoint( 4 ),
                                              this->pGetPoint( 5 ) ) ) );
        faces.push_back( FacePointerType( new FaceType2(
                                              this->pGetPoint( 1 ),
                                              this->pGetPoint( 2 ),
                                              this->pGetPoint( 5 ),
                                              this->pGetPoint( 4 ) ) ) );
        faces.push_back( FacePointerType( new FaceType2(
                                              this->pGetPoint( 0 ),
                                              this->pGetPoint( 3 ),
                                              this->pGetPoint( 5 ),
                                              this->pGetPoint( 2 ) ) ) );
        faces.push_back( FacePointerType( new FaceType2(
                                              this->pGetPoint( 0 ),
                                              this->pGetPoint( 1 ),
                                              this->pGetPoint( 4 ),
                                              this->pGetPoint( 3 ) ) ) );
        return faces;
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
        switch ( ShapeFunctionIndex )
        {
        case 0:
            return( 1.0 -( rPoint[0] + rPoint[1] + rPoint[2] - ( rPoint[0]*rPoint[2] ) - ( rPoint[1]*rPoint[2] ) ) );
        case 1:
            return( rPoint[0] - ( rPoint[0]*rPoint[2] ) );
        case 2:
            return( rPoint[1] - ( rPoint[1]*rPoint[2] ) );
        case 3:
            return( rPoint[2] - ( rPoint[0]*rPoint[2] ) - ( rPoint[1]*rPoint[2] ) );
        case 4:
            return( rPoint[0]*rPoint[2] );
        case 5:
            return( rPoint[1]*rPoint[2] );
        default:
            KRATOS_ERROR << "Wrong index of shape function!" << *this  << std::endl;
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
      if(rResult.size() != 6) rResult.resize(6,false);
      rResult[0] =  1.0 -( rCoordinates[0] + rCoordinates[1] + rCoordinates[2] - ( rCoordinates[0]*rCoordinates[2] ) - ( rCoordinates[1]*rCoordinates[2] ) );
      rResult[1] =  rCoordinates[0] - ( rCoordinates[0]*rCoordinates[2] );
      rResult[2] =  rCoordinates[1] - ( rCoordinates[1]*rCoordinates[2] );
      rResult[3] =  rCoordinates[2] - ( rCoordinates[0]*rCoordinates[2] ) - ( rCoordinates[1]*rCoordinates[2] );
      rResult[4] =  rCoordinates[0]*rCoordinates[2];
      rResult[5] =  rCoordinates[1]*rCoordinates[2];

        return rResult;
    }

    /**
     * TODO: implemented but not yet tested
     */
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
        const unsigned int integration_points_number =
            msGeometryData.IntegrationPointsNumber( ThisMethod );

        if ( integration_points_number == 0 )
            KRATOS_ERROR << "This integration method is not supported" << *this << std::endl;

        //workaround by riccardo
        if ( rResult.size() != integration_points_number )
        {
            // KLUDGE: While there is a bug in ublas
            // vector resize, I have to put this beside resizing!!
            ShapeFunctionsGradientsType temp( integration_points_number );
            rResult.swap( temp );
        }

        //calculating the local gradients
        ShapeFunctionsGradientsType locG =
            CalculateShapeFunctionsIntegrationPointsLocalGradients( ThisMethod );

        //getting the inverse jacobian matrices
        JacobiansType temp( integration_points_number );

        JacobiansType invJ = InverseOfJacobian( temp, ThisMethod );

        //loop over all integration points
        for ( unsigned int pnt = 0; pnt < integration_points_number; pnt++ )
        {
            rResult[pnt].resize( 6, 3, false );

            for ( int i = 0; i < 6; i++ )
            {
                for ( int j = 0; j < 3; j++ )
                {
                    rResult[pnt]( i, j ) =
                        ( locG[pnt]( i, 0 ) * invJ[pnt]( j, 0 ) )
                        + ( locG[pnt]( i, 1 ) * invJ[pnt]( j, 1 ) )
                        + ( locG[pnt]( i, 2 ) * invJ[pnt]( j, 2 ) );
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
        return "3 dimensional interface Prism with six nodes in 3D space";
    }

    /**
     * Print information about this object.
     * @param rOStream Stream to print into it.
     * @see PrintData()
     * @see Info()
     */
    void PrintInfo( std::ostream& rOStream ) const override
    {
        rOStream << "3 dimensional interface Prism with six nodes in 3D space";
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
        rResult.resize( 6, 3, false );
        noalias( rResult ) = ZeroMatrix( 6, 3 );

        rResult( 0, 0 ) = -1.0 + rPoint[2];
        rResult( 0, 1 ) = -1.0 + rPoint[2];
        rResult( 0, 2 ) =  -1.0 + rPoint[0] + rPoint[1];

        rResult( 1, 0 ) =  1.0 - rPoint[2];
        rResult( 1, 1 ) =  0.0;
        rResult( 1, 2 ) =  -rPoint[0];

        rResult( 2, 0 ) =  0.0;
        rResult( 2, 1 ) =  1.0 - rPoint[2];
        rResult( 2, 2 ) =  -rPoint[1];

        rResult( 3, 0 ) = -rPoint[2];
        rResult( 3, 1 ) = -rPoint[2];
        rResult( 3, 2 ) =  1.0 - rPoint[0] - rPoint[1];

        rResult( 4, 0 ) =  rPoint[2];
        rResult( 4, 1 ) =  0.0;
        rResult( 4, 2 ) =  rPoint[0];

        rResult( 5, 0 ) =  0.0;
        rResult( 5, 1 ) =  rPoint[2];
        rResult( 5, 2 ) =  rPoint[1];

        return rResult;
    }


    /**
     * returns the second order derivatives of all shape functions
     * in given arbitrary pointers
     * @param rResult a third order tensor which contains the second derivatives
     * @param rPoint the given point the second order derivatives are calculated in
     */
    ShapeFunctionsSecondDerivativesType& ShapeFunctionsSecondDerivatives( ShapeFunctionsSecondDerivativesType& rResult,
                                                                          const CoordinatesArrayType& rPoint ) const override
    {
        if ( rResult.size() != this->PointsNumber() )
        {
            // KLUDGE: While there is a bug in
            // ublas vector resize, I have to put this beside resizing!!
            ShapeFunctionsGradientsType temp( this->PointsNumber() );
            rResult.swap( temp );
        }

        //TODO: this is not correct for a prism

        rResult[0].resize( 2, 2 , false);
        rResult[1].resize( 2, 2 , false);
        rResult[2].resize( 2, 2 , false);

        rResult[0]( 0, 0 ) = 0.0;
        rResult[0]( 0, 1 ) = 0.0;
        rResult[0]( 1, 0 ) = 0.0;
        rResult[0]( 1, 1 ) = 0.0;
        rResult[1]( 0, 0 ) = 0.0;
        rResult[1]( 0, 1 ) = 0.0;
        rResult[1]( 1, 0 ) = 0.0;
        rResult[1]( 1, 1 ) = 0.0;
        rResult[2]( 0, 0 ) = 0.0;
        rResult[2]( 0, 1 ) = 0.0;
        rResult[2]( 1, 0 ) = 0.0;
        rResult[2]( 1, 1 ) = 0.0;
        return rResult;
    }

    /**
     * returns the third order derivatives of all shape functions
     * in given arbitrary pointers
     * @param rResult a fourth order tensor which contains the third derivatives
     * @param rPoint the given point the third order derivatives are calculated in
     */
    ShapeFunctionsThirdDerivativesType& ShapeFunctionsThirdDerivatives( ShapeFunctionsThirdDerivativesType& rResult,
                                                                        const CoordinatesArrayType& rPoint ) const override
    {
        if ( rResult.size() != this->PointsNumber() )
        {
            // KLUDGE: While there is a bug in
            // ublas vector resize, I have to put this beside resizing!!
//                 ShapeFunctionsGradientsType
            ShapeFunctionsThirdDerivativesType temp( this->PointsNumber() );
            rResult.swap( temp );
        }

        for ( IndexType i = 0; i < rResult.size(); i++ )
        {
            DenseVector<Matrix> temp( this->PointsNumber() );
            rResult[i].swap( temp );
        }

        //TODO: this is not correct for a prism

        rResult[0][0].resize( 2, 2 , false);
        rResult[0][1].resize( 2, 2 , false);
        rResult[1][0].resize( 2, 2 , false);
        rResult[1][1].resize( 2, 2 , false);
        rResult[2][0].resize( 2, 2 , false);
        rResult[2][1].resize( 2, 2 , false);

        for ( int i = 0; i < 3; i++ )
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
     * There are no protected members in class PrismInterface3D6
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
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BaseType );
    }

    void load( Serializer& rSerializer ) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BaseType );
    }

    PrismInterface3D6(): BaseType( PointsArrayType(), &msGeometryData ) {}

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
        IntegrationPointsArrayType integration_points =
            all_integration_points[ThisMethod];
        //number of integration points
        const int integration_points_number = integration_points.size();
        //number of nodes in current geometry
        const int points_number = 6;
        //setting up return matrix
        Matrix shape_function_values( integration_points_number, points_number );
        //loop over all integration points

        for ( int pnt = 0; pnt < integration_points_number; pnt++ )
        {
            shape_function_values ( pnt, 0 ) =  ( 1.0
                                                - integration_points[pnt].X()
                                                - integration_points[pnt].Y()
                                                - integration_points[pnt].Z()
                                                + ( integration_points[pnt].X() * integration_points[pnt].Z() )
                                                + ( integration_points[pnt].Y() * integration_points[pnt].Z() ) );
            shape_function_values( pnt, 1 ) = integration_points[pnt].X()
                                              - ( integration_points[pnt].X() * integration_points[pnt].Z() );
            shape_function_values( pnt, 2 ) = integration_points[pnt].Y()
                                              - ( integration_points[pnt].Y() * integration_points[pnt].Z() );
            shape_function_values( pnt, 3 ) = integration_points[pnt].Z()
                                              - ( integration_points[pnt].X() * integration_points[pnt].Z() )
                                              - ( integration_points[pnt].Y() * integration_points[pnt].Z() );
            shape_function_values( pnt, 4 ) = ( integration_points[pnt].X() * integration_points[pnt].Z() );
            shape_function_values( pnt, 5 ) = ( integration_points[pnt].Y() * integration_points[pnt].Z() );
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
    static ShapeFunctionsGradientsType CalculateShapeFunctionsIntegrationPointsLocalGradients(
        typename BaseType::IntegrationMethod ThisMethod )
    {
        IntegrationPointsContainerType all_integration_points =
            AllIntegrationPoints();
        IntegrationPointsArrayType integration_points =
            all_integration_points[ThisMethod];
        //number of integration points
        const int integration_points_number = integration_points.size();
        ShapeFunctionsGradientsType d_shape_f_values( integration_points_number );
        //initialising container
        //std::fill(d_shape_f_values.begin(), d_shape_f_values.end(), Matrix(4,2));
        //loop over all integration points

        for ( int pnt = 0; pnt < integration_points_number; pnt++ )
        {
            Matrix result = ZeroMatrix( 6, 3 );
            result( 0, 0 ) = -1.0 + integration_points[pnt].Z();
            result( 0, 1 ) = -1.0 + integration_points[pnt].Z();
            result( 0, 2 ) = -1.0 + integration_points[pnt].X() + integration_points[pnt].Y();
            result( 1, 0 ) =  1.0 - integration_points[pnt].Z();
            result( 1, 1 ) =  0.0;
            result( 1, 2 ) =  -integration_points[pnt].X();
            result( 2, 0 ) =  0.0;
            result( 2, 1 ) =  1.0 - integration_points[pnt].Z();
            result( 2, 2 ) =  -integration_points[pnt].Y();
            result( 3, 0 ) =  -integration_points[pnt].Z();
            result( 3, 1 ) =  -integration_points[pnt].Z();
            result( 3, 2 ) =  1.0 - integration_points[pnt].X() - integration_points[pnt].Y();
            result( 4, 0 ) =  integration_points[pnt].Z();
            result( 4, 1 ) =  0.0;
            result( 4, 2 ) =  integration_points[pnt].X();
            result( 5, 0 ) =  0.0;
            result( 5, 1 ) =  integration_points[pnt].Z();
            result( 5, 2 ) =  integration_points[pnt].Y();
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
                Quadrature<PrismGaussLobattoIntegrationPoints1,
                3, IntegrationPoint<3> >::GenerateIntegrationPoints(),
                Quadrature<PrismGaussLobattoIntegrationPoints2,
                3, IntegrationPoint<3> >::GenerateIntegrationPoints(),
                IntegrationPointsArrayType(),
                IntegrationPointsArrayType()
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
                PrismInterface3D6<TPointType>::CalculateShapeFunctionsIntegrationPointsValues(
                    GeometryData::GI_GAUSS_1 ),
                PrismInterface3D6<TPointType>::CalculateShapeFunctionsIntegrationPointsValues(
                    GeometryData::GI_GAUSS_2 ),
                Matrix(),
                Matrix()
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
                PrismInterface3D6<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients( GeometryData::GI_GAUSS_1 ),
                PrismInterface3D6<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients( GeometryData::GI_GAUSS_2 ),
                ShapeFunctionsGradientsType(),
                ShapeFunctionsGradientsType(),
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

    template<class TOtherPointType> friend class PrismInterface3D6;

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
    PrismInterface3D6<TPointType>& rThis );
/**
 * output stream functions
 */
template<class TPointType> inline std::ostream& operator << (
    std::ostream& rOStream,
    const PrismInterface3D6<TPointType>& rThis )
{
    rThis.PrintInfo( rOStream );
    rOStream << std::endl;
    rThis.PrintData( rOStream );
    return rOStream;
}

///@}

template<class TPointType> const
GeometryData PrismInterface3D6<TPointType>::msGeometryData(
    3, 3, 2, GeometryData::GI_GAUSS_2,
    PrismInterface3D6<TPointType>::AllIntegrationPoints(),
    PrismInterface3D6<TPointType>::AllShapeFunctionsValues(),
    AllShapeFunctionsLocalGradients()
);
}// namespace Kratos.

#endif // KRATOS_PRISM_INTERFACE_3D_6_H_INCLUDED  defined
