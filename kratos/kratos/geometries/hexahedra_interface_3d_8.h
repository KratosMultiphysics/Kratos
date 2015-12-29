//
//   Project Name:        Kratos
//   Last Modified by:    $Author:   JMCarbonell $
//   Date:                $Date:   December 2015 $
//   Revision:            $Revision:         1.7 $
//
//

#if !defined(KRATOS_HEXAHEDRA_INTERFACE_3D_8_H_INCLUDED )
#define  KRATOS_HEXAHEDRA_INTERFACE_3D_8_H_INCLUDED


// System includes

// External includes

// Project includes
#include "geometries/quadrilateral_3d_4.h"
#include "integration/hexahedron_gauss_lobatto_integration_points.h"



namespace Kratos
{
/**
 * An eight node hexahedra geometry with linear shape functions
 */

template<class TPointType> class HexahedraInterface3D8 : public Geometry<TPointType>
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
     * Pointer definition of HexahedraInterface3D8
     */
    KRATOS_CLASS_POINTER_DEFINITION( HexahedraInterface3D8 );

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
    HexahedraInterface3D8( const PointType& Point1, const PointType& Point2,
                  const PointType& Point3, const PointType& Point4,
                  const PointType& Point5, const PointType& Point6,
                  const PointType& Point7, const PointType& Point8 )
        : BaseType( PointsArrayType(), &msGeometryData )
    { 
        array_1d< double , 3 > vx;
        vx.clear();
        noalias(vx) += - Point1 - Point5;
        noalias(vx) += - Point4 - Point8;
        noalias(vx) += Point2 + Point6;
        noalias(vx) += Point3 + Point7;

        vx *= 0.25;

        array_1d< double , 3 > vy;
        vy.clear();
        noalias(vy) += Point4 + Point3;
        noalias(vy) += Point7 + Point8;
        noalias(vy) += - Point1 - Point2;
        noalias(vy) += - Point6 - Point5;

        vy *= 0.25;

        array_1d< double , 3 > vz;
        vz.clear();
        noalias(vz) += Point6 + Point6;
        noalias(vz) += Point7 + Point8;
        noalias(vz) += - Point1 - Point2;
        noalias(vz) += - Point3 - Point4;

        vz *= 0.25;

        double lx = MathUtils<double>::Norm3(vx);
        double ly = MathUtils<double>::Norm3(vy);
        double lz = MathUtils<double>::Norm3(vz);
		if(lz < lx)
		{
			if(lz < ly)
			{
				// LZ
				this->Points().push_back( typename PointType::Pointer( new PointType( Point1 ) ) );
				this->Points().push_back( typename PointType::Pointer( new PointType( Point2 ) ) );
				this->Points().push_back( typename PointType::Pointer( new PointType( Point3 ) ) );
				this->Points().push_back( typename PointType::Pointer( new PointType( Point4 ) ) );
				this->Points().push_back( typename PointType::Pointer( new PointType( Point5 ) ) );
				this->Points().push_back( typename PointType::Pointer( new PointType( Point6 ) ) );
				this->Points().push_back( typename PointType::Pointer( new PointType( Point7 ) ) );
				this->Points().push_back( typename PointType::Pointer( new PointType( Point8 ) ) );
			}
			else
			{
				if(ly < lx)
				{
					// LY
					this->Points().push_back( typename PointType::Pointer( new PointType( Point2 ) ) );
					this->Points().push_back( typename PointType::Pointer( new PointType( Point1 ) ) );
					this->Points().push_back( typename PointType::Pointer( new PointType( Point5 ) ) );
					this->Points().push_back( typename PointType::Pointer( new PointType( Point6 ) ) );
					this->Points().push_back( typename PointType::Pointer( new PointType( Point3 ) ) );
					this->Points().push_back( typename PointType::Pointer( new PointType( Point4 ) ) );
					this->Points().push_back( typename PointType::Pointer( new PointType( Point8 ) ) );
					this->Points().push_back( typename PointType::Pointer( new PointType( Point7 ) ) );
				}
				else
				{
					// LX
					this->Points().push_back( typename PointType::Pointer( new PointType( Point1 ) ) );
					this->Points().push_back( typename PointType::Pointer( new PointType( Point4 ) ) );
					this->Points().push_back( typename PointType::Pointer( new PointType( Point8 ) ) );
					this->Points().push_back( typename PointType::Pointer( new PointType( Point5 ) ) );
					this->Points().push_back( typename PointType::Pointer( new PointType( Point2 ) ) );
					this->Points().push_back( typename PointType::Pointer( new PointType( Point3 ) ) );
					this->Points().push_back( typename PointType::Pointer( new PointType( Point7 ) ) );
					this->Points().push_back( typename PointType::Pointer( new PointType( Point6 ) ) );
				}
			}
		}
		else
		{
			if(lx < ly)
			{
				// LX
				this->Points().push_back( typename PointType::Pointer( new PointType( Point1 ) ) );
				this->Points().push_back( typename PointType::Pointer( new PointType( Point4 ) ) );
				this->Points().push_back( typename PointType::Pointer( new PointType( Point8 ) ) );
				this->Points().push_back( typename PointType::Pointer( new PointType( Point5 ) ) );
				this->Points().push_back( typename PointType::Pointer( new PointType( Point2 ) ) );
				this->Points().push_back( typename PointType::Pointer( new PointType( Point3 ) ) );
				this->Points().push_back( typename PointType::Pointer( new PointType( Point7 ) ) );
				this->Points().push_back( typename PointType::Pointer( new PointType( Point6 ) ) );
			}
			else
			{
				if(lz < ly)
				{
					// LZ
					this->Points().push_back( typename PointType::Pointer( new PointType( Point1 ) ) );
					this->Points().push_back( typename PointType::Pointer( new PointType( Point2 ) ) );
					this->Points().push_back( typename PointType::Pointer( new PointType( Point3 ) ) );
					this->Points().push_back( typename PointType::Pointer( new PointType( Point4 ) ) );
					this->Points().push_back( typename PointType::Pointer( new PointType( Point5 ) ) );
					this->Points().push_back( typename PointType::Pointer( new PointType( Point6 ) ) );
					this->Points().push_back( typename PointType::Pointer( new PointType( Point7 ) ) );
					this->Points().push_back( typename PointType::Pointer( new PointType( Point8 ) ) );
				}
				else
				{
					// LY
					this->Points().push_back( typename PointType::Pointer( new PointType( Point2 ) ) );
					this->Points().push_back( typename PointType::Pointer( new PointType( Point1 ) ) );
					this->Points().push_back( typename PointType::Pointer( new PointType( Point5 ) ) );
					this->Points().push_back( typename PointType::Pointer( new PointType( Point6 ) ) );
					this->Points().push_back( typename PointType::Pointer( new PointType( Point3 ) ) );
					this->Points().push_back( typename PointType::Pointer( new PointType( Point4 ) ) );
					this->Points().push_back( typename PointType::Pointer( new PointType( Point8 ) ) );
					this->Points().push_back( typename PointType::Pointer( new PointType( Point7 ) ) );
				}
			}
		}
	}

    HexahedraInterface3D8( typename PointType::Pointer pPoint1,
                  typename PointType::Pointer pPoint2,
                  typename PointType::Pointer pPoint3,
                  typename PointType::Pointer pPoint4,
                  typename PointType::Pointer pPoint5,
                  typename PointType::Pointer pPoint6,
                  typename PointType::Pointer pPoint7,
                  typename PointType::Pointer pPoint8 )
        : BaseType( PointsArrayType(), &msGeometryData )
    {
        array_1d< double , 3 > vx;
        vx.clear();
        noalias(vx) += - *pPoint1 - *pPoint5;
        noalias(vx) += - *pPoint4 - *pPoint8;
        noalias(vx) += *pPoint2 + *pPoint6;
        noalias(vx) += *pPoint3 + *pPoint7;

        vx *= 0.25;

        array_1d< double , 3 > vy;
        vy.clear();
        noalias(vy) += *pPoint4 + *pPoint3;
        noalias(vy) += *pPoint7 + *pPoint8;
        noalias(vy) += - *pPoint1 - *pPoint2;
        noalias(vy) += - *pPoint6 - *pPoint5;

        vy *= 0.25;

        array_1d< double , 3 > vz;
        vz.clear();
        noalias(vz) += *pPoint6 + *pPoint6;
        noalias(vz) += *pPoint7 + *pPoint8;
        noalias(vz) += - *pPoint1 - *pPoint2;
        noalias(vz) += - *pPoint3 - *pPoint4;

        vz *= 0.25;

        double lx = MathUtils<double>::Norm3(vx);
        double ly = MathUtils<double>::Norm3(vy);
        double lz = MathUtils<double>::Norm3(vz);
		if(lz < lx)
		{
			if(lz < ly)
			{
				// LZ
				this->Points().push_back( pPoint1 );
				this->Points().push_back( pPoint2 );
				this->Points().push_back( pPoint3 );
				this->Points().push_back( pPoint4 );
				this->Points().push_back( pPoint5 );
				this->Points().push_back( pPoint6 );
				this->Points().push_back( pPoint7 );
				this->Points().push_back( pPoint8 );
			}
			else
			{
				if(ly < lx)
				{
					// LY
					this->Points().push_back( pPoint2 );
					this->Points().push_back( pPoint1 );
					this->Points().push_back( pPoint5 );
					this->Points().push_back( pPoint6 );
					this->Points().push_back( pPoint3 );
					this->Points().push_back( pPoint4 );
					this->Points().push_back( pPoint8 );
					this->Points().push_back( pPoint7 );
				}
				else
				{
					// LX
					this->Points().push_back( pPoint1 );
					this->Points().push_back( pPoint4 );
					this->Points().push_back( pPoint8 );
					this->Points().push_back( pPoint5 );
					this->Points().push_back( pPoint2 );
					this->Points().push_back( pPoint3 );
					this->Points().push_back( pPoint7 );
					this->Points().push_back( pPoint6 );
				}
			}
		}
		else
		{
			if(lx < ly)
			{
				// LX
				this->Points().push_back( pPoint1 );
				this->Points().push_back( pPoint4 );
				this->Points().push_back( pPoint8 );
				this->Points().push_back( pPoint5 );
				this->Points().push_back( pPoint2 );
				this->Points().push_back( pPoint3 );
				this->Points().push_back( pPoint7 );
				this->Points().push_back( pPoint6 );
			}
			else
			{
				if(lz < ly)
				{
					// LZ
					this->Points().push_back( pPoint1 );
					this->Points().push_back( pPoint2 );
					this->Points().push_back( pPoint3 );
					this->Points().push_back( pPoint4 );
					this->Points().push_back( pPoint5 );
					this->Points().push_back( pPoint6 );
					this->Points().push_back( pPoint7 );
					this->Points().push_back( pPoint8 );
				}
				else
				{
					// LY
					this->Points().push_back( pPoint2 );
					this->Points().push_back( pPoint1 );
					this->Points().push_back( pPoint5 );
					this->Points().push_back( pPoint6 );
					this->Points().push_back( pPoint3 );
					this->Points().push_back( pPoint4 );
					this->Points().push_back( pPoint8 );
					this->Points().push_back( pPoint7 );
				}
			}
		}
    }

    HexahedraInterface3D8( const PointsArrayType& ThisPoints )
        : BaseType( PointsArrayType(), &msGeometryData )
    {
        if ( ThisPoints.size() != 8 )
            KRATOS_THROW_ERROR( std::invalid_argument,
                          "Invalid points number. Expected 8, given " ,
                          this->PointsNumber() );

		const typename PointType::Pointer& pPoint1 = ThisPoints(0);
		const typename PointType::Pointer& pPoint2 = ThisPoints(1);
		const typename PointType::Pointer& pPoint3 = ThisPoints(2);
		const typename PointType::Pointer& pPoint4 = ThisPoints(3);
		const typename PointType::Pointer& pPoint5 = ThisPoints(4);
		const typename PointType::Pointer& pPoint6 = ThisPoints(5);
		const typename PointType::Pointer& pPoint7 = ThisPoints(6);
		const typename PointType::Pointer& pPoint8 = ThisPoints(7);

        array_1d< double , 3 > vx;
        vx.clear();
        vx = - (*pPoint1) - (*pPoint5) - (*pPoint4) - (*pPoint8) + (*pPoint2) + (*pPoint6) + (*pPoint3) + (*pPoint7);
        vx *= 0.25;

        array_1d< double , 3 > vy;
        vy.clear();
        vy = (*pPoint4) + (*pPoint3) +(*pPoint7) + (*pPoint8) - (*pPoint1) - (*pPoint2) - (*pPoint6) - (*pPoint5);
        vy *= 0.25;

        Vector vz(3);
        vz.clear();
        vz = (*pPoint6) + (*pPoint5) + (*pPoint7) + (*pPoint8) - (*pPoint1) - (*pPoint2) - (*pPoint3) - (*pPoint4);
        vz *= 0.25;

        double lx = MathUtils<double>::Norm3(vx);
        double ly = MathUtils<double>::Norm3(vy);
        double lz = MathUtils<double>::Norm3(vz);
		if(lz < lx)
		{
			if(lz < ly)
			{
				// LZ
				this->Points().push_back(pPoint1);
				this->Points().push_back(pPoint2);
				this->Points().push_back(pPoint3);
				this->Points().push_back(pPoint4);
				this->Points().push_back(pPoint5);
				this->Points().push_back(pPoint6);
				this->Points().push_back(pPoint7);
				this->Points().push_back(pPoint8);
			}
			else
			{
				if(ly < lx)
				{
					// LY
					this->Points().push_back(pPoint2);
					this->Points().push_back(pPoint1);
					this->Points().push_back(pPoint5);
					this->Points().push_back(pPoint6);
					this->Points().push_back(pPoint3);
					this->Points().push_back(pPoint4);
					this->Points().push_back(pPoint8);
					this->Points().push_back(pPoint7);
				}
				else
				{
					// LX
					this->Points().push_back(pPoint1);
					this->Points().push_back(pPoint4);
					this->Points().push_back(pPoint8);
					this->Points().push_back(pPoint5);
					this->Points().push_back(pPoint2);
					this->Points().push_back(pPoint3);
					this->Points().push_back(pPoint7);
					this->Points().push_back(pPoint6);
				}
			}
		}
		else
		{
			if(lx < ly)
			{
				// LX
				this->Points().push_back(pPoint1);
				this->Points().push_back(pPoint4);
				this->Points().push_back(pPoint8);
				this->Points().push_back(pPoint5);
				this->Points().push_back(pPoint2);
				this->Points().push_back(pPoint3);
				this->Points().push_back(pPoint7);
				this->Points().push_back(pPoint6);
			}
			else
			{
				if(lz < ly)
				{
					// LZ
					this->Points().push_back(pPoint1);
					this->Points().push_back(pPoint2);
					this->Points().push_back(pPoint3);
					this->Points().push_back(pPoint4);
					this->Points().push_back(pPoint5);
					this->Points().push_back(pPoint6);
					this->Points().push_back(pPoint7);
					this->Points().push_back(pPoint8);
				}
				else
				{
					// LY
					this->Points().push_back(pPoint2);
					this->Points().push_back(pPoint1);
					this->Points().push_back(pPoint5);
					this->Points().push_back(pPoint6);
					this->Points().push_back(pPoint3);
					this->Points().push_back(pPoint4);
					this->Points().push_back(pPoint8);
					this->Points().push_back(pPoint7);
				}
			}
		}
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
    HexahedraInterface3D8( HexahedraInterface3D8 const& rOther )
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
    template<class TOtherPointType> HexahedraInterface3D8( HexahedraInterface3D8<TOtherPointType> const& rOther )
        : BaseType( rOther )
    {
    }

    /// Destructor. Does nothing!!!
    virtual ~HexahedraInterface3D8() {}

    GeometryData::KratosGeometryFamily GetGeometryFamily()
    {
        return GeometryData::Kratos_Hexahedra;
    }

    GeometryData::KratosGeometryType GetGeometryType()
    {
        return GeometryData::Kratos_Hexahedra3D8;
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
    HexahedraInterface3D8& operator=( const HexahedraInterface3D8& rOther )
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
    HexahedraInterface3D8& operator=( HexahedraInterface3D8<TOtherPointType> const & rOther )
    {
        BaseType::operator=( rOther );

        return *this;
    }


    /**
     * Operations
     */

    typename BaseType::Pointer Create( PointsArrayType const& ThisPoints ) const
    {
        return typename BaseType::Pointer( new HexahedraInterface3D8( ThisPoints ) );
    }

    virtual boost::shared_ptr< Geometry< Point<3> > > Clone() const
    {
        Geometry< Point<3> >::PointsArrayType NewPoints;
        //making a copy of the nodes TO POINTS (not Nodes!!!)

        for ( IndexType i = 0 ; i < this->Points().size() ; i++ )
            NewPoints.push_back( this->Points()[i] );

        //creating a geometry with the new points
        boost::shared_ptr< Geometry< Point<3> > >
        p_clone( new HexahedraInterface3D8< Point<3> >( NewPoints ) );

        p_clone->ClonePoints();

        return p_clone;
    }

    //lumping factors for the calculation of the lumped mass matrix
    virtual Vector& LumpingFactors( Vector& rResult ) const
    {
	if(rResult.size() != 8)
           rResult.resize( 8, false );
        std::fill( rResult.begin(), rResult.end(), 1.00 / 8.00 );
        return rResult;
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
    virtual double Length() const
    {
        return std::sqrt(Area());
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
    virtual double Area() const
    {
        //return fabs( DeterminantOfJacobian( PointType() ) ) * 0.5;
        
        array_1d<double, 3> p1 = 0.5 * (BaseType::GetPoint( 0 ) + BaseType::GetPoint( 4 ));
		array_1d<double, 3> p2 = 0.5 * (BaseType::GetPoint( 1 ) + BaseType::GetPoint( 5 ));
		array_1d<double, 3> p3 = 0.5 * (BaseType::GetPoint( 2 ) + BaseType::GetPoint( 6 ));
		array_1d<double, 3> p4 = 0.5 * (BaseType::GetPoint( 3 ) + BaseType::GetPoint( 7 ));
		
        //const TPointType& p1 = this->Points()[0];
		//const TPointType& p2 = this->Points()[1];
		//const TPointType& p3 = this->Points()[2];
		//const TPointType& p4 = this->Points()[3];

		double p1x = p1[0];
		double p1y = p1[1];
		double p1z = p1[2];

		double p2x = p2[0];
		double p2y = p2[1];
		double p2z = p2[2];

		double p3x = p3[0];
		double p3y = p3[1];
		double p3z = p3[2];

		double p4x = p4[0];
		double p4y = p4[1];
		double p4z = p4[2];

		double pos = 0.5 + 0.5 / std::sqrt(3.0);
		double w = 0.25;

		double C1  = pos*(p1z - p2z + p3z - p4z);
		double C2  = pos*(p1y - p2y + p3y - p4y);
		double C3  = pos*(p1x - p2x + p3x - p4x);
		double C4  = C1 - p1z + p2z;
		double C5  = C1 + p1z - p2z;
		double C6  = C2 + p1y - p2y;
		double C7  = C2 - p1y + p2y;
		double C8  = C3 - p1x + p2x;
		double C9  = C3 + p1x - p2x;
		double C10 = C1 - p1z + p4z;
		double C11 = C2 - p1y + p4y;
		double C12 = C3 - p1x + p4x;
		double C13 = C1 + p1z - p4z;
		double C14 = C2 + p1y - p4y;
		double C15 = C3 + p1x - p4x;

		return w * (
			std::sqrt( std::pow(C4*C11 - C7*C10, 2) + std::pow(C4*C12 - C8*C10, 2) + std::pow(C7*C12 - C8*C11, 2)) + 
			std::sqrt( std::pow(C5*C11 - C6*C10, 2) + std::pow(C5*C12 - C9*C10, 2) + std::pow(C6*C12 - C9*C11, 2)) + 
			std::sqrt( std::pow(C4*C14 - C7*C13, 2) + std::pow(C4*C15 - C8*C13, 2) + std::pow(C7*C15 - C8*C14, 2)) + 
			std::sqrt( std::pow(C5*C14 - C6*C13, 2) + std::pow(C5*C15 - C9*C13, 2) + std::pow(C6*C15 - C9*C14, 2))
			);
    }

    virtual double Volume() const
    {
		return Area();
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
    virtual double DomainSize() const
    {
        return Area();
    }

    /**
     * Returns a matrix of the local coordinates of all points
     * @param rResult a Matrix that will be overwritten by the results
     * @return the coordinates of all points of the current geometry
     */
    virtual Matrix& PointsLocalCoordinates( Matrix& rResult ) const
    {
        if ( rResult.size1() != 8 || rResult.size2() != 3 )
            rResult.resize( 8, 3 );

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
     * Returns whether given arbitrary point is inside the Geometry
     */
    virtual bool IsInside( const CoordinatesArrayType& rPoint, CoordinatesArrayType& rResult )
    {
        this->PointLocalCoordinates( rResult, rPoint );

        if ( fabs( rResult[0] ) < 1 + 1.0e-8 )
            if ( fabs( rResult[1] ) < 1 + 1.0e-8 )
                if ( fabs( rResult[2] ) < 1 + 1.0e-8 )
                    return true;

        return false;
    }


//             virtual void Bounding_Box(BoundingBox<TPointType, BaseType>& rResult) const
//             {
//                 //rResult.Geometry() = *(this);
//                 BaseType::Bounding_Box(rResult.LowPoint(), rResult.HighPoint());
//             }

    /**
     * Jacobian
     */

    /**
     * :TODO: TO BE TESTED
     */
    /**
     * Jacobians for given method.
     * This method calculates jacobians matrices in all integrations
     * points of given integration method.
     *
     * @param ThisMethod integration method which jacobians has to
     * be calculated in its integration points.
     * @return JacobiansType a Vector of jacobian
     * matrices \f$ J_i \f$ where \f$ i=1,2,...,n \f$ is the
     * integration point index of given integration method.
     *
     * @see DeterminantOfJacobian
     * @see InverseOfJacobian
     */
    virtual JacobiansType& Jacobian( JacobiansType& rResult,
                                     IntegrationMethod ThisMethod ) const
    {
		array_1d< array_1d<double, 3>, 4 > pmid;
		pmid[0] = 0.5 * (BaseType::GetPoint( 0 ) + BaseType::GetPoint( 4 ));
		pmid[1] = 0.5 * (BaseType::GetPoint( 1 ) + BaseType::GetPoint( 5 ));
		pmid[2] = 0.5 * (BaseType::GetPoint( 2 ) + BaseType::GetPoint( 6 ));
		pmid[3] = 0.5 * (BaseType::GetPoint( 3 ) + BaseType::GetPoint( 7 ));

        //getting derivatives of shape functions
        ShapeFunctionsGradientsType shape_functions_gradients =
            CalculateShapeFunctionsIntegrationPointsLocalGradients( ThisMethod );
        //getting values of shape functions
        Matrix shape_functions_values =
            CalculateShapeFunctionsIntegrationPointsValues( ThisMethod );
        //workaround by riccardo...

        if ( rResult.size() != this->IntegrationPointsNumber( ThisMethod ) )
        {
             //KLUDGE: While there is a bug in ublas
             //vector resize, I have to put this beside resizing!!
            JacobiansType temp( this->IntegrationPointsNumber( ThisMethod ) );
            rResult.swap( temp );
        }

        //loop over all integration points
        for ( unsigned int pnt = 0; pnt < this->IntegrationPointsNumber( ThisMethod )/2; pnt++ )
        {
            //defining single jacobian matrix
            Matrix jacobian(3,2,0.0);
            
            //loop over all nodes
            for ( unsigned int i = 0; i < pmid.size(); i++ )
            {
                jacobian( 0, 0 ) += ( pmid[i][0] ) * ( shape_functions_gradients[pnt]( i, 0 ) );
                jacobian( 0, 1 ) += ( pmid[i][0] ) * ( shape_functions_gradients[pnt]( i, 1 ) );
                jacobian( 1, 0 ) += ( pmid[i][1] ) * ( shape_functions_gradients[pnt]( i, 0 ) );
                jacobian( 1, 1 ) += ( pmid[i][1] ) * ( shape_functions_gradients[pnt]( i, 1 ) );
                jacobian( 2, 0 ) += ( pmid[i][2] ) * ( shape_functions_gradients[pnt]( i, 0 ) );
                jacobian( 2, 1 ) += ( pmid[i][2] ) * ( shape_functions_gradients[pnt]( i, 1 ) );
			}
			
			rResult[pnt] = jacobian;
			rResult[pnt+4] = jacobian;
        }//end of loop over all integration points
        return rResult;
    }

    /**
     * Jacobians for given  method.
     * This method calculates jacobians matrices in all integrations
     * points of given integration method.
     *
     * @param ThisMethod integration method which jacobians has to
     * be calculated in its integration points.
     *
     * @return JacobiansType a Vector of jacobian
     * matrices \f$ J_i \f$ where \f$ i=1,2,...,n \f$ is the
     * integration point index of given integration method.
     *
     * @param DeltaPosition Matrix with the nodes position increment which describes
     * the configuration where the jacobian has to be calculated.     
     *
     * @see DeterminantOfJacobian
     * @see InverseOfJacobian
     */
    virtual JacobiansType& Jacobian( JacobiansType& rResult,
                                     IntegrationMethod ThisMethod, 
									 Matrix & DeltaPosition ) const

    {
		array_1d< array_1d<double, 3>, 4 > pmid;
		pmid[0] = 0.5 * (BaseType::GetPoint( 0 ) + BaseType::GetPoint( 4 ));
		pmid[1] = 0.5 * (BaseType::GetPoint( 1 ) + BaseType::GetPoint( 5 ));
		pmid[2] = 0.5 * (BaseType::GetPoint( 2 ) + BaseType::GetPoint( 6 ));
		pmid[3] = 0.5 * (BaseType::GetPoint( 3 ) + BaseType::GetPoint( 7 ));
		
		Matrix deltaPosMid(3, 3);
		for(unsigned int i = 0; i < 3; i++)
			for(unsigned int j = 0; j < 3; j++)
				deltaPosMid(i, j) = 0.5*( DeltaPosition(i, j) + DeltaPosition(i+4, j) );
		
        //getting derivatives of shape functions
        ShapeFunctionsGradientsType shape_functions_gradients =
            CalculateShapeFunctionsIntegrationPointsLocalGradients( ThisMethod );

        //workaround by riccardo...
        if ( rResult.size() != this->IntegrationPointsNumber( ThisMethod ) )
        {
            // KLUDGE: While there is a bug in ublas
            // vector resize, I have to put this beside resizing!!
            JacobiansType temp( this->IntegrationPointsNumber( ThisMethod ) );
            rResult.swap( temp );
        }

        //loop over all integration points
        for ( unsigned int pnt = 0; pnt < this->IntegrationPointsNumber( ThisMethod )/2; pnt++ )
        {
            //defining single jacobian matrix
            Matrix jacobian(3,2,0.0);
            //loop over all nodes
            for ( unsigned int i = 0; i < pmid.size(); i++ )
            {
                jacobian( 0, 0 ) += ( pmid[i][0] - deltaPosMid(i,0) ) * ( shape_functions_gradients[pnt]( i, 0 ) );
                jacobian( 0, 1 ) += ( pmid[i][0] - deltaPosMid(i,0) ) * ( shape_functions_gradients[pnt]( i, 1 ) );
                jacobian( 1, 0 ) += ( pmid[i][1] - deltaPosMid(i,1) ) * ( shape_functions_gradients[pnt]( i, 0 ) );
                jacobian( 1, 1 ) += ( pmid[i][1] - deltaPosMid(i,1) ) * ( shape_functions_gradients[pnt]( i, 1 ) );
                jacobian( 2, 0 ) += ( pmid[i][2] - deltaPosMid(i,2) ) * ( shape_functions_gradients[pnt]( i, 0 ) );
                jacobian( 2, 1 ) += ( pmid[i][2] - deltaPosMid(i,2) ) * ( shape_functions_gradients[pnt]( i, 1 ) );
            }

            rResult[pnt] = jacobian;
			rResult[pnt+4] = jacobian;
        }//end of loop over all integration points

        return rResult;
    }

    /**
     * :TODO: TO BE TESTED
     */
    /** Jacobian in specific integration point of given integration
     *  method. This method calculate jacobian matrix in given
     * integration point of given integration method.
     *
     * @param IntegrationPointIndex index of integration point which jacobians has to
     * be calculated in it.
     *
     * @param ThisMethod integration method which jacobians has to
     * be calculated in its integration points.
     *
     * @return Matrix(double) Jacobian matrix \f$ J_i \f$ where \f$
     * i \f$ is the given integration point index of given
     * integration method.
     *
     * @see DeterminantOfJacobian
     * @see InverseOfJacobian
     */
    virtual Matrix& Jacobian( Matrix& rResult,
                              IndexType IntegrationPointIndex,
                              IntegrationMethod ThisMethod ) const
    {
		array_1d< array_1d<double, 3>, 4 > pmid;
		pmid[0] = 0.5 * (BaseType::GetPoint( 0 ) + BaseType::GetPoint( 4 ));
		pmid[1] = 0.5 * (BaseType::GetPoint( 1 ) + BaseType::GetPoint( 5 ));
		pmid[2] = 0.5 * (BaseType::GetPoint( 2 ) + BaseType::GetPoint( 6 ));
		pmid[3] = 0.5 * (BaseType::GetPoint( 3 ) + BaseType::GetPoint( 7 ));
		
        //setting up size of jacobian matrix
        rResult.resize( 3, 2 );
        //derivatives of shape functions
        ShapeFunctionsGradientsType shape_functions_gradients =
                CalculateShapeFunctionsIntegrationPointsLocalGradients( ThisMethod );
        Matrix& ShapeFunctionsGradientInIntegrationPoint = shape_functions_gradients( IntegrationPointIndex );

        //Elements of jacobian matrix (e.g. J(1,1) = dX1/dXi1)
        //loop over all nodes

        for ( unsigned int i = 0; i < pmid.size(); i++ )
        {
            rResult( 0, 0 ) += ( pmid[i][0] ) * ( ShapeFunctionsGradientInIntegrationPoint( i, 0 ) );
            rResult( 0, 1 ) += ( pmid[i][0] ) * ( ShapeFunctionsGradientInIntegrationPoint( i, 1 ) );
            rResult( 1, 0 ) += ( pmid[i][1] ) * ( ShapeFunctionsGradientInIntegrationPoint( i, 0 ) );
            rResult( 1, 1 ) += ( pmid[i][1] ) * ( ShapeFunctionsGradientInIntegrationPoint( i, 1 ) );
            rResult( 2, 0 ) += ( pmid[i][2] ) * ( ShapeFunctionsGradientInIntegrationPoint( i, 0 ) );
            rResult( 2, 1 ) += ( pmid[i][2] ) * ( ShapeFunctionsGradientInIntegrationPoint( i, 1 ) );
        }

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
    virtual Matrix& Jacobian( Matrix& rResult, 
                              IndexType IntegrationPointIndex, 
                              IntegrationMethod ThisMethod, 
                              Matrix& DeltaPosition ) const
    {
        array_1d< array_1d<double, 3>, 4 > pmid;
		pmid[0] = 0.5 * (BaseType::GetPoint( 0 ) + BaseType::GetPoint( 4 ));
		pmid[1] = 0.5 * (BaseType::GetPoint( 1 ) + BaseType::GetPoint( 5 ));
		pmid[2] = 0.5 * (BaseType::GetPoint( 2 ) + BaseType::GetPoint( 6 ));
		pmid[3] = 0.5 * (BaseType::GetPoint( 3 ) + BaseType::GetPoint( 7 ));
		
		Matrix deltaPosMid(4, 3);
		for(unsigned int i = 0; i < 4; i++)
			for(unsigned int j = 0; j < 3; j++)
				deltaPosMid(i, j) = 0.5*( DeltaPosition(i, j) + DeltaPosition(i+4, j) );
		
        //getting derivatives of shape functions
        ShapeFunctionsGradientsType shape_functions_gradients =
            CalculateShapeFunctionsIntegrationPointsLocalGradients( ThisMethod );

        if(rResult.size1() != 3 || rResult.size2() != 2)
			rResult.resize(3, 2, false);
		noalias(rResult) = ZeroMatrix(3, 2);

        //loop over all nodes
        for ( unsigned int i = 0; i < pmid.size(); i++ )
        {
            rResult( 0, 0 ) += ( pmid[i][0] - deltaPosMid(i,0) ) * ( shape_functions_gradients[IntegrationPointIndex]( i, 0 ) );
            rResult( 0, 1 ) += ( pmid[i][0] - deltaPosMid(i,0) ) * ( shape_functions_gradients[IntegrationPointIndex]( i, 1 ) );
            rResult( 1, 0 ) += ( pmid[i][1] - deltaPosMid(i,1) ) * ( shape_functions_gradients[IntegrationPointIndex]( i, 0 ) );
            rResult( 1, 1 ) += ( pmid[i][1] - deltaPosMid(i,1) ) * ( shape_functions_gradients[IntegrationPointIndex]( i, 1 ) );
            rResult( 2, 0 ) += ( pmid[i][2] - deltaPosMid(i,2) ) * ( shape_functions_gradients[IntegrationPointIndex]( i, 0 ) );
            rResult( 2, 1 ) += ( pmid[i][2] - deltaPosMid(i,2) ) * ( shape_functions_gradients[IntegrationPointIndex]( i, 1 ) );
        }

        return rResult;
    }
    
    /**
     * :TODO: TO BE TESTED
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
    virtual Matrix& Jacobian( Matrix& rResult, const CoordinatesArrayType& rPoint ) const
    {
		array_1d< array_1d<double, 3>, 4 > pmid;
		pmid[0] = 0.5 * (BaseType::GetPoint( 0 ) + BaseType::GetPoint( 4 ));
		pmid[1] = 0.5 * (BaseType::GetPoint( 1 ) + BaseType::GetPoint( 5 ));
		pmid[2] = 0.5 * (BaseType::GetPoint( 2 ) + BaseType::GetPoint( 6 ));
		pmid[3] = 0.5 * (BaseType::GetPoint( 3 ) + BaseType::GetPoint( 7 ));
		
        //setting up size of jacobian matrix
        rResult.resize( 3, 2 );
        //derivatives of shape functions
        Matrix shape_functions_gradients = ShapeFunctionsLocalGradients( rPoint );
        //Elements of jacobian matrix (e.g. J(1,1) = dX1/dXi1)
        //loop over all nodes

        for ( unsigned int i = 0; i < pmid.size(); i++ )
        {
            rResult( 0, 0 ) += ( pmid[i][0] ) * ( shape_functions_gradients( i, 0 ) );
            rResult( 0, 1 ) += ( pmid[i][0] ) * ( shape_functions_gradients( i, 1 ) );
            rResult( 1, 0 ) += ( pmid[i][1] ) * ( shape_functions_gradients( i, 0 ) );
            rResult( 1, 1 ) += ( pmid[i][1] ) * ( shape_functions_gradients( i, 1 ) );
            rResult( 2, 0 ) += ( pmid[i][2] ) * ( shape_functions_gradients( i, 0 ) );
            rResult( 2, 1 ) += ( pmid[i][2] ) * ( shape_functions_gradients( i, 1 ) );
        }

        return rResult;
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
    virtual Vector& DeterminantOfJacobian( Vector& rResult,
                                           IntegrationMethod ThisMethod ) const
    {
        KRATOS_THROW_ERROR( std::logic_error, "Quadrilateral3D4::DeterminantOfJacobian", "Jacobian is not square" );
        return rResult;
    }

    /**
     * :TODO: TO BE TESTED
     */
    /**
     * Determinant of jacobian in specific integration point of
     * given integration method. This method calculate determinant
     * of jacobian in given integration point of given integration
     * method.
     *
     * @param IntegrationPointIndex index of integration point which
     * jacobians has to be calculated in it.
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
    virtual double DeterminantOfJacobian( IndexType IntegrationPointIndex,
                                          IntegrationMethod ThisMethod ) const
    {
        KRATOS_THROW_ERROR( std::logic_error, "Quadrilateral3D4::DeterminantOfJacobian", "Jacobian is not square" );
        return 0.0;
    }

    /**
     * :TODO: TO BE TESTED
     */
    /**
     * Determinant of jacobian in given point.
     * This method calculate determinant of jacobian
     * matrix in given point.
     *
     * @param rPoint point which determinant of jacobians has to
     * be calculated in it.
     * @return Determinamt of jacobian matrix \f$ |J| \f$ in given
     * point.
     *
     * @see DeterminantOfJacobian
     * @see InverseOfJacobian
     *
     * KLUDGE: PointType needed for proper functionality
     * KLUDGE: works only with explicitly generated Matrix object
     */
    virtual double DeterminantOfJacobian( const CoordinatesArrayType& rPoint ) const
    {
        KRATOS_THROW_ERROR( std::logic_error, "Quadrilateral3D4::DeterminantOfJacobian", "Jacobian is not square" );
        return 0.0;
    }

    /**
     * :TODO: TO BE TESTED
     */
    /**
     * Inverse of jacobians for given integration method.
     * This method calculate inverse of jacobians matrices in all
     * integrations points of given integration method.
     *
     * @param ThisMethod integration method which inverse of jacobians has to
     * be calculated in its integration points.
     * @return Inverse of jacobian
     * matrices \f$ J^{-1}_i \f$ where \f$ i=1,2,...,n \f$ is the integration
     * point index of given integration method.
     *
     * @see Jacobian
     * @see DeterminantOfJacobian
     *
     * KLUDGE: works only with explicitly generated Matrix object
     */
    virtual JacobiansType& InverseOfJacobian( JacobiansType& rResult,
            IntegrationMethod ThisMethod ) const
    {
        KRATOS_THROW_ERROR( std::logic_error, "Quadrilateral3D4::DeterminantOfJacobian", "Jacobian is not square" );
        return rResult;
    }

    /**
     * :TODO: TO BE TESTED
     */
    /**
     * Inverse of jacobian in specific integration point of given integration
     * method.
     * This method calculate Inverse of jacobian matrix in given
     * integration point of given integration method.
     *
     * @param IntegrationPointIndex index of integration point which
     * inverse of jacobians has to be calculated in it.
     *
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
    virtual Matrix& InverseOfJacobian( Matrix& rResult,
                                       IndexType IntegrationPointIndex,
                                       IntegrationMethod ThisMethod ) const
    {
        KRATOS_THROW_ERROR( std::logic_error, "Quadrilateral3D4::DeterminantOfJacobian", "Jacobian is not square" );
        return rResult;
    }

    /**
     * :TODO: TO BE TESTED
     */
    /**
     * Inverse of jacobian in given point.
     * This method calculate inverse of jacobian matrix in given point.
     *
     * @param rPoint point which inverse of jacobians has to
     * be calculated in it.
     * @return Inverse of jacobian matrix \f$ J^{-1} \f$ in given point.
     *
     * @see DeterminantOfJacobian
     * @see InverseOfJacobian
     *
     * KLUDGE: works only with explicitly generated Matrix object
     */
    virtual Matrix& InverseOfJacobian( Matrix& rResult, const CoordinatesArrayType& rPoint ) const
    {
        KRATOS_THROW_ERROR( std::logic_error, "Quadrilateral3D4::DeterminantOfJacobian", "Jacobian is not square" );
        return rResult;
    }


    /** This method gives you number of all edges of this
    geometry.
    @return SizeType containes number of this geometry edges.
    @see Edges()
    @see Edge()
     */
    // will be used by refinement algorithm, thus uncommented. janosch.
    virtual SizeType EdgesNumber() const
    {
        return 12;
    }

    virtual SizeType FacesNumber() const
    {
        return 6;
    }

    /** This method gives you all edges of this geometry.

    @return GeometriesArrayType containes this geometry edges.
    @see EdgesNumber()
    @see Edge()
     */
    virtual GeometriesArrayType Edges( void )
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

    virtual GeometriesArrayType Faces( void )
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
    virtual double ShapeFunctionValue( IndexType ShapeFunctionIndex,
                                       const CoordinatesArrayType& rPoint ) const
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
        case 4:
            return( 0.25*( 1.0 - rPoint[0] )*( 1.0 - rPoint[1] ) );
        case 5:
            return( 0.25*( 1.0 + rPoint[0] )*( 1.0 - rPoint[1] ) );
        case 6:
            return( 0.25*( 1.0 + rPoint[0] )*( 1.0 + rPoint[1] ) );
        case 7:
            return( 0.25*( 1.0 - rPoint[0] )*( 1.0 + rPoint[1] ) );
        default:
            KRATOS_THROW_ERROR( std::logic_error,
                          "Wrong index of shape function!" , *this );
        }

        return 0;
    }


    /**
     * Calculates the gradients in terms of local coordinateds
     * of all shape functions in a given point.
     *
     * @param rPoint the current point at which the gradients are calculated
     * @return the gradients of all shape functions
     * \f$ \frac{\partial N^i}{\partial \xi_j} \f$
     */
    virtual Matrix& ShapeFunctionsLocalGradients( Matrix& rResult, const CoordinatesArrayType& rPoint ) const
    {
		KRATOS_THROW_ERROR( std::logic_error, "This integration method is not supported" , *this );
		return rResult;
    }


    /**
     * Calculates the Gradients of the shape functions.
     * Calculates the gradients of the shape functions with regard to the global
     * coordinates in all
     * integration points (\f$ \frac{\partial N^i}{\partial X_j} \f$)
     *
     * @param rResult a container which takes the calculated gradients
     * @param ThisMethod the given IntegrationMethod
     * @return the gradients of all shape functions with regard to the global coordinates
     *
     * KLUDGE: method call only works with explicit JacobiansType rather than creating
     * JacobiansType within argument list
     *
     * :TODO: TESTING!!!
     */
    virtual ShapeFunctionsGradientsType& ShapeFunctionsIntegrationPointsGradients(
        ShapeFunctionsGradientsType& rResult,
        IntegrationMethod ThisMethod ) const
    {
        const unsigned int integration_points_number =
            msGeometryData.IntegrationPointsNumber( ThisMethod );

        if ( integration_points_number == 0 )
            KRATOS_THROW_ERROR( std::logic_error,
                          "This integration method is not supported" , *this );

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
            rResult[pnt].resize( 4, 3 );

            for ( int i = 0; i < 4; i++ )
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


    virtual ShapeFunctionsGradientsType& ShapeFunctionsIntegrationPointsGradients(
            ShapeFunctionsGradientsType& rResult,
            Vector& determinants_of_jacobian,
            IntegrationMethod ThisMethod) const
    {
        const unsigned int integration_points_number = msGeometryData.IntegrationPointsNumber(ThisMethod);

        if ( integration_points_number == 0 )
            KRATOS_THROW_ERROR( std::logic_error,"This integration method is not supported" , *this );

        //workaround by riccardo
        if ( rResult.size() != integration_points_number )
        {
            // KLUDGE: While there is a bug in ublas
            // vector resize, I have to put this beside resizing!!
            ShapeFunctionsGradientsType temp( integration_points_number );
            rResult.swap( temp );
        }

        if ( determinants_of_jacobian.size() != integration_points_number)
            determinants_of_jacobian.resize(integration_points_number,false);

        //calculating the local gradients
        ShapeFunctionsGradientsType locG =
            CalculateShapeFunctionsIntegrationPointsLocalGradients( ThisMethod );

        JacobiansType J( integration_points_number );
        Jacobian(J,ThisMethod);
//        JacobiansType invJ = InverseOfJacobian( temp, ThisMethod );

        //loop over all integration points
        for ( unsigned int pnt = 0; pnt < integration_points_number; pnt++ )
        {
            //current jacobian
//            Matrix J = ZeroMatrix( 3, 3 );
//            this->Jacobian( J, pnt, ThisMethod );
//            double DetJ = DeterminantOfJacobian(pnt,ThisMethod);

            Matrix invJ = ZeroMatrix( 3, 3 );
            double DetJ;

            MathUtils<double>::InvertMatrix3( J[pnt], invJ, DetJ );

            determinants_of_jacobian[pnt] = DetJ;

            rResult[pnt].resize( 4, 3 );

            for ( int i = 0; i < 4; i++ )
            {
                for ( int j = 0; j < 3; j++ )
                {
                    rResult[pnt]( i, j ) =
//                        ( locG[pnt]( i, 0 ) * invJ[pnt]( j, 0 ) )
//                        + ( locG[pnt]( i, 1 ) * invJ[pnt]( j, 1 ) )
//                        + ( locG[pnt]( i, 2 ) * invJ[pnt]( j, 2 ) );
                    ( locG[pnt]( i, 0 ) * invJ( 0, j ) )
                    + ( locG[pnt]( i, 1 ) * invJ( 1, j ) )
                    + ( locG[pnt]( i, 2 ) * invJ( 2, j ) );
                }
            }
        }

        return rResult;
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
    virtual std::string Info() const
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
    virtual void PrintInfo( std::ostream& rOStream ) const
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
    virtual void PrintData( std::ostream& rOStream ) const
    {
        BaseType::PrintData( rOStream );
        std::cout << std::endl;
        Matrix jacobian;
        Jacobian( jacobian, PointType() );
        rOStream << "    Jacobian in the origin\t : " << jacobian;
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


    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    virtual void save( Serializer& rSerializer ) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, PointsArrayType );
    }

    virtual void load( Serializer& rSerializer )
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, PointsArrayType );
    }

    HexahedraInterface3D8(): BaseType( PointsArrayType(), &msGeometryData ) {}

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
        result( 0, 0 ) = -0.25 * ( 1.0 - rPoint[1] );
        result( 0, 1 ) = -0.25 * ( 1.0 - rPoint[0] );
        result( 0, 2 ) = 0.0;
        result( 1, 0 ) = 0.25 * ( 1.0 - rPoint[1] );
        result( 1, 1 ) = -0.25 * ( 1.0 + rPoint[0] );
        result( 1, 2 ) = 0.0;
        result( 2, 0 ) = 0.25 * ( 1.0 + rPoint[1] );
        result( 2, 1 ) = 0.25 * ( 1.0 + rPoint[0] );
        result( 2, 2 ) = 0.0;
        result( 3, 0 ) = -0.25 * ( 1.0 + rPoint[1] );
        result( 3, 1 ) = 0.25 * ( 1.0 - rPoint[0] );
        result( 3, 2 ) = 0.0;
        result( 4, 0 ) = -0.25 * ( 1.0 - rPoint[1] );
        result( 4, 1 ) = -0.25 * ( 1.0 - rPoint[0] );
        result( 4, 2 ) = 0.0;
        result( 5, 0 ) = 0.25 * ( 1.0 - rPoint[1] );
        result( 5, 1 ) = -0.25 * ( 1.0 + rPoint[0] );
        result( 5, 2 ) = 0.0;
        result( 6, 0 ) = 0.25 * ( 1.0 + rPoint[1] );
        result( 6, 1 ) = 0.25 * ( 1.0 + rPoint[0] );
        result( 6, 2 ) = 0.0;
        result( 7, 0 ) = -0.25 * ( 1.0 + rPoint[1] );
        result( 7, 1 ) = 0.25 * ( 1.0 - rPoint[0] );
        result( 7, 2 ) = 0.0;
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
        IntegrationPointsArrayType& integration_points = all_integration_points[ThisMethod];

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
            shape_function_values( pnt, 4 ) =
                0.25 * ( 1.0 - integration_points[pnt].X() )
                * ( 1.0 - integration_points[pnt].Y() );
            shape_function_values( pnt, 5 ) =
                0.25 * ( 1.0 + integration_points[pnt].X() )
                * ( 1.0 - integration_points[pnt].Y() );
            shape_function_values( pnt, 6 ) =
                0.25 * ( 1.0 + integration_points[pnt].X() )
                * ( 1.0 + integration_points[pnt].Y() );
            shape_function_values( pnt, 7 ) =
                0.25 * ( 1.0 - integration_points[pnt].X() )
                * ( 1.0 + integration_points[pnt].Y() );
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
        IntegrationPointsArrayType integration_points =
            all_integration_points[ThisMethod];
        //number of integration points
        const int integration_points_number = integration_points.size();
        ShapeFunctionsGradientsType d_shape_f_values( integration_points_number );
        //initialising container
        //loop over all integration points

        for ( int pnt = 0; pnt < integration_points_number; pnt++ )
        {
            Matrix& result = d_shape_f_values[pnt];
            result = ZeroMatrix( 8, 3 );
            result( 0, 0 ) = -0.25 * ( 1.0 - integration_points[pnt].Y() );
            result( 0, 1 ) = -0.25 * ( 1.0 - integration_points[pnt].X() );
            result( 0, 2 ) = 0.0;
            result( 1, 0 ) = 0.25 * ( 1.0 - integration_points[pnt].Y() );
            result( 1, 1 ) = -0.25 * ( 1.0 + integration_points[pnt].X() );
            result( 1, 2 ) = 0.0;
            result( 2, 0 ) = 0.25 * ( 1.0 + integration_points[pnt].Y() );
            result( 2, 1 ) = 0.25 * ( 1.0 + integration_points[pnt].X() );
            result( 2, 2 ) = 0.0;
            result( 3, 0 ) = -0.25 * ( 1.0 + integration_points[pnt].Y() );
            result( 3, 1 ) = 0.25 * ( 1.0 - integration_points[pnt].X() );
            result( 3, 2 ) = 0.0;
            result( 4, 0 ) = -0.25 * ( 1.0 - integration_points[pnt].Y() );
            result( 4, 1 ) = -0.25 * ( 1.0 - integration_points[pnt].X() );
            result( 4, 2 ) = 0.0;
            result( 5, 0 ) = 0.25 * ( 1.0 - integration_points[pnt].Y() );
            result( 5, 1 ) = -0.25 * ( 1.0 + integration_points[pnt].X() );
            result( 5, 2 ) = 0.0;
            result( 6, 0 ) = 0.25 * ( 1.0 + integration_points[pnt].Y() );
            result( 6, 1 ) = 0.25 * ( 1.0 + integration_points[pnt].X() );
            result( 6, 2 ) = 0.0;
            result( 7, 0 ) = -0.25 * ( 1.0 + integration_points[pnt].Y() );
            result( 7, 1 ) = 0.25 * ( 1.0 - integration_points[pnt].X() );
            result( 7, 2 ) = 0.0;
        }

        return d_shape_f_values;
    }

    static const IntegrationPointsContainerType AllIntegrationPoints()
    {
        IntegrationPointsContainerType integration_points =
        {
            {
                IntegrationPointsArrayType(),
                Quadrature < HexaedronGaussLobattoIntegrationPoints2,
                3, IntegrationPoint<3> >::GenerateIntegrationPoints(),
                IntegrationPointsArrayType(),
                IntegrationPointsArrayType()
            }
        };
        return integration_points;
    }

    static const ShapeFunctionsValuesContainerType AllShapeFunctionsValues()
    {
        ShapeFunctionsValuesContainerType shape_functions_values =
        {
            {
                Matrix(),
                HexahedraInterface3D8<TPointType>::CalculateShapeFunctionsIntegrationPointsValues(
                    GeometryData::GI_GAUSS_2 ),
                Matrix(),
                Matrix()
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
                ShapeFunctionsGradientsType(),
                HexahedraInterface3D8<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients(
                    GeometryData::GI_GAUSS_2 ),
                ShapeFunctionsGradientsType(),
                ShapeFunctionsGradientsType()
            }
        };
        return shape_functions_local_gradients;
    }


    /**
     * Private Friends
     */

    template<class TOtherPointType> friend class HexahedraInterface3D8;


    /**
     * Un accessible methods
     */

};// Class HexahedraInterface3D8


/**
 * Input and output
 */

/**
 * input stream function
 */
template<class TPointType> inline std::istream& operator >> (
    std::istream& rIStream, HexahedraInterface3D8<TPointType>& rThis );

/**
 * output stream function
 */
template<class TPointType> inline std::ostream& operator << (
    std::ostream& rOStream, const HexahedraInterface3D8<TPointType>& rThis )
{
    rThis.PrintInfo( rOStream );
    rOStream << std::endl;
    rThis.PrintData( rOStream );

    return rOStream;
}


template<class TPointType> const
GeometryData HexahedraInterface3D8<TPointType>::msGeometryData(
    3, 3, 3, GeometryData::GI_GAUSS_2,
    HexahedraInterface3D8<TPointType>::AllIntegrationPoints(),
    HexahedraInterface3D8<TPointType>::AllShapeFunctionsValues(),
    AllShapeFunctionsLocalGradients()
);

}// namespace Kratos.

#endif // KRATOS_HEXAHEDRA_INTERFACE_3D_8_H_INCLUDED  defined 

