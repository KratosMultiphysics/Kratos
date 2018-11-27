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
//

#if !defined(KRATOS_TRIANGLE_3D_3_H_INCLUDED )
#define  KRATOS_TRIANGLE_3D_3_H_INCLUDED

// System includes

// External includes

// Project includes
#include "geometries/plane_3d.h"
#include "geometries/line_3d_2.h"
#include "integration/triangle_gauss_legendre_integration_points.h"
#include "integration/triangle_collocation_integration_points.h"

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
 * @class Triangle3D3
 * @ingroup KratosCore
 * @brief A three node 3D triangle geometry with linear shape functions
 * @details While the shape functions are only defined in 2D it is possible to define an arbitrary orientation in space. Thus it can be used for defining surfaces on 3D elements.
 * The node ordering corresponds with: 
 *      v                                                              
 *      ^                                                               
 *      |                                                              
 *      2                                   
 *      |`\                   
 *      |  `\                   
 *      |    `\                 
 *      |      `\                
 *      |        `\                 
 *      0----------1 --> u  
 * @author Riccardo Rossi
 * @author Janosch Stascheit
 * @author Felix Nagel
 */
template<class TPointType> class Triangle3D3
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

    typedef Geometry<TPointType> GeometryType;

    /**
     * Type of edge geometry
     */
    typedef Line3D2<TPointType> EdgeType;

    /**
     * Pointer definition of Triangle3D3
     */
    KRATOS_CLASS_POINTER_DEFINITION( Triangle3D3 );

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

//     Triangle3D3( const PointType& FirstPoint,
//                  const PointType& SecondPoint,
//                  const PointType& ThirdPoint )
//         : BaseType( PointsArrayType(), &msGeometryData )
//     {
//         this->Points().push_back( typename PointType::Pointer( new PointType( FirstPoint ) ) );
//         this->Points().push_back( typename PointType::Pointer( new PointType( SecondPoint ) ) );
//         this->Points().push_back( typename PointType::Pointer( new PointType( ThirdPoint ) ) );
//     }

    Triangle3D3( typename PointType::Pointer pFirstPoint,
                 typename PointType::Pointer pSecondPoint,
                 typename PointType::Pointer pThirdPoint )
        : BaseType( PointsArrayType(), &msGeometryData )
    {
        this->Points().push_back( pFirstPoint );
        this->Points().push_back( pSecondPoint );
        this->Points().push_back( pThirdPoint );
    }

    explicit Triangle3D3( const PointsArrayType& ThisPoints )
        : BaseType( ThisPoints, &msGeometryData )
    {
        if ( this->PointsNumber() != 3 )
            KRATOS_ERROR << "Invalid points number. Expected 3, given " << this->PointsNumber() << std::endl;
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
    Triangle3D3( Triangle3D3 const& rOther )
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
    template<class TOtherPointType> explicit Triangle3D3( Triangle3D3<TOtherPointType> const& rOther )
        : BaseType( rOther )
    {
    }

    /**
     * Destructor. Does nothing!!!
     */
    ~Triangle3D3() override {}

    GeometryData::KratosGeometryFamily GetGeometryFamily() const override
    {
        return GeometryData::Kratos_Triangle;
    }

    GeometryData::KratosGeometryType GetGeometryType() const override
    {
        return GeometryData::Kratos_Triangle3D3;
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
    Triangle3D3& operator=( const Triangle3D3& rOther )
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
    Triangle3D3& operator=( Triangle3D3<TOtherPointType> const & rOther )
    {
        BaseType::operator=( rOther );
        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    typename BaseType::Pointer Create( PointsArrayType const& ThisPoints ) const override
    {
        return typename BaseType::Pointer( new Triangle3D3( ThisPoints ) );
    }


    // Kratos::shared_ptr< Geometry< Point<3> > > Clone() const override
    // {
    //     Geometry< Point<3> >::PointsArrayType NewPoints;

    //     //making a copy of the nodes TO POINTS (not Nodes!!!)
    //     for ( IndexType i = 0 ; i < this->size() ; i++ )
    //     {
    //             NewPoints.push_back(Kratos::make_shared< Point<3> >(( *this )[i]));
    //     }

    //     //creating a geometry with the new points
    //     Geometry< Point<3> >::Pointer p_clone( new Triangle3D3< Point<3> >( NewPoints ) );

    //     return p_clone;
    // }

    /**
     * returns the local coordinates of all nodes of the current geometry
     * @param rResult a Matrix object that will be overwritten by the result
     * @return the local coordinates of all nodes
     */
    Matrix& PointsLocalCoordinates( Matrix& rResult ) const override
    {
        rResult.resize( 3, 2 ,false);
        noalias( rResult ) = ZeroMatrix( 3, 2 );
        rResult( 0, 0 ) =  0.0;
        rResult( 0, 1 ) =  0.0;
        rResult( 1, 0 ) =  1.0;
        rResult( 1, 1 ) =  0.0;
        rResult( 2, 0 ) =  0.0;
        rResult( 2, 1 ) =  1.0;
        return rResult;
    }

    //lumping factors for the calculation of the lumped mass matrix
    Vector& LumpingFactors( Vector& rResult ) const override
    {
        rResult.resize( 3, false );
        std::fill( rResult.begin(), rResult.end(), 1.00 / 3.00 );
        return rResult;
    }

    ///@}
    ///@name Information
    ///@{

    /** Calculates and returns Length or charactereistic length of this geometry.
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
        return std::sqrt(2.0 * Area());
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
        const double a = MathUtils<double>::Norm3(this->GetPoint(0)-this->GetPoint(1));
        const double b = MathUtils<double>::Norm3(this->GetPoint(1)-this->GetPoint(2));
        const double c = MathUtils<double>::Norm3(this->GetPoint(2)-this->GetPoint(0));

        const double s = (a+b+c) / 2.0;

        return std::sqrt(s*(s-a)*(s-b)*(s-c));
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

    /// Class Interface

    /** This method calculates and returns the minimum edge
     * length of the geometry
     *
     * @return double value with the minimum edge length
     *
     * @see MaxEdgeLength()
     * @see AverageEdgeLength()
     */
    double MinEdgeLength() const override
    {
        const array_1d<double, 3> a = this->GetPoint(0) - this->GetPoint(1);
        const array_1d<double, 3> b = this->GetPoint(1) - this->GetPoint(2);
        const array_1d<double, 3> c = this->GetPoint(2) - this->GetPoint(0);

        const double sa = (a[0]*a[0])+(a[1]*a[1])+(a[2]*a[2]);
        const double sb = (b[0]*b[0])+(b[1]*b[1])+(b[2]*b[2]);
        const double sc = (c[0]*c[0])+(c[1]*c[1])+(c[2]*c[2]);

        return CalculateMinEdgeLength(sa, sb, sc);
    }

    /** This method calculates and returns the maximum edge
     * length of the geometry
     *
     * @return double value with the maximum edge length
     *
     * @see MinEdgeLength()
     * @see AverageEdgeLength()
     */
    double MaxEdgeLength() const override
    {
        const array_1d<double, 3> a = this->GetPoint(0) - this->GetPoint(1);
        const array_1d<double, 3> b = this->GetPoint(1) - this->GetPoint(2);
        const array_1d<double, 3> c = this->GetPoint(2) - this->GetPoint(0);

        const double sa = (a[0]*a[0])+(a[1]*a[1])+(a[2]*a[2]);
        const double sb = (b[0]*b[0])+(b[1]*b[1])+(b[2]*b[2]);
        const double sc = (c[0]*c[0])+(c[1]*c[1])+(c[2]*c[2]);

        return CalculateMaxEdgeLength(sa, sb, sc);
    }

    /** This method calculates and returns the average edge
     * length of the geometry
     *
     * @return double value with the average edge length
     *
     * @see MinEdgeLength()
     * @see MaxEdgeLength()
     */
    double AverageEdgeLength() const override
    {
        return CalculateAvgEdgeLength(
        MathUtils<double>::Norm3(this->GetPoint(0)-this->GetPoint(1)),
        MathUtils<double>::Norm3(this->GetPoint(1)-this->GetPoint(2)),
        MathUtils<double>::Norm3(this->GetPoint(2)-this->GetPoint(0))
        );
    }

    /** Calculates the circumradius of the geometry.
     * Calculates the circumradius of the geometry.
     *
     * @return Circumradius of the geometry.
     *
     * @see Inradius()
     */
    double Circumradius() const override
    {
        return CalculateCircumradius(
        MathUtils<double>::Norm3(this->GetPoint(0)-this->GetPoint(1)),
        MathUtils<double>::Norm3(this->GetPoint(1)-this->GetPoint(2)),
        MathUtils<double>::Norm3(this->GetPoint(2)-this->GetPoint(0))
        );
    }

    /** Calculates the inradius of the geometry.
     * Calculates the inradius of the geometry.
     *
     * @return Inradius of the geometry.
     *
     * @see Circumradius()
     */
    double Inradius() const override
    {
        return CalculateInradius(
        MathUtils<double>::Norm3(this->GetPoint(0)-this->GetPoint(1)),
        MathUtils<double>::Norm3(this->GetPoint(1)-this->GetPoint(2)),
        MathUtils<double>::Norm3(this->GetPoint(2)-this->GetPoint(0))
        );
    }


	bool AllSameSide(array_1d<double, 3> const& Distances)
    {
        constexpr double epsilon = std::numeric_limits<double>::epsilon();

        // put U0,U1,U2 into plane equation 1 to compute signed distances to the plane//
        double du0 = Distances[0];
        double du1 = Distances[1];
        double du2 = Distances[2];

        // coplanarity robustness check //
        if (std::abs(du0)<epsilon) du0 = 0.0;
        if (std::abs(du1)<epsilon) du1 = 0.0;
        if (std::abs(du2)<epsilon) du2 = 0.0;

        const double du0du1 = du0*du1;
        const double du0du2 = du0*du2;

        if (du0du1>0.00 && du0du2>0.00)// same sign on all of them + not equal 0 ? //
            return true;                   // no intersection occurs //

        return false;

	}

	int GetMajorAxis(array_1d<double, 3> const& V)
    {
        int index = static_cast<int>(std::abs(V[0]) < std::abs(V[1]));
        return (std::abs(V[index]) > std::abs(V[2])) ? index : 2;
	}

	bool HasIntersection(const GeometryType& ThisGeometry) override
	{
        // Based on code develop by Moller: http://fileadmin.cs.lth.se/cs/Personal/Tomas_Akenine-Moller/code/opttritri.txt
        // and the article "A Fast Triangle-Triangle Intersection Test", Journal of Graphics Tools, 2(2), 1997:
        // http://web.stanford.edu/class/cs277/resources/papers/Moller1997b.pdf

        Plane3D plane_1(this->GetPoint(0), this->GetPoint(1), this->GetPoint(2));
        array_1d<double, 3> distances_1;
        for (int i = 0; i < 3; i++)
            distances_1[i] = plane_1.CalculateSignedDistance(ThisGeometry[i]);
        if (AllSameSide(distances_1))
            return false;

        Plane3D plane_2(ThisGeometry[0], ThisGeometry[1], ThisGeometry[2]);
        array_1d<double, 3> distances_2;
        for (int i = 0; i < 3; i++)
            distances_2[i] = plane_2.CalculateSignedDistance(this->GetPoint(i));
        if (AllSameSide(distances_2))
            return false;

        // compute direction of intersection line //
        array_1d<double, 3> intersection_direction;
        MathUtils<double>::CrossProduct(intersection_direction, plane_1.GetNormal(), plane_2.GetNormal());

        int index = GetMajorAxis(intersection_direction);

        // this is the simplified projection onto L//
        double vp0 = this->GetPoint(0)[index];
        double vp1 = this->GetPoint(1)[index];
        double vp2 = this->GetPoint(2)[index];

        double up0 = ThisGeometry[0][index];
        double up1 = ThisGeometry[1][index];
        double up2 = ThisGeometry[2][index];


        // compute interval for triangle 1 //
        double a, b, c, x0, x1;
        if (ComputeIntervals(vp0, vp1, vp2, distances_2[0], distances_2[1], distances_2[2], a, b, c, x0, x1) == true)
        {
            return CoplanarIntersectionCheck(plane_1.GetNormal(), ThisGeometry);
        }

        // compute interval for triangle 2 //
        double d, e, f, y0, y1;
        if (ComputeIntervals(up0, up1, up2, distances_1[0], distances_1[1], distances_1[2], d, e, f, y0, y1) == true)
        {
            return CoplanarIntersectionCheck(plane_1.GetNormal(), ThisGeometry);
        }


        double xx, yy, xxyy, tmp;
        xx = x0*x1;
        yy = y0*y1;
        xxyy = xx*yy;

        array_1d<double, 2> isect1, isect2;

        tmp = a*xxyy;
        isect1[0] = tmp + b*x1*yy;
        isect1[1] = tmp + c*x0*yy;

        tmp = d*xxyy;
        isect2[0] = tmp + e*xx*y1;
        isect2[1] = tmp + f*xx*y0;

        std::sort(isect1.begin(), isect1.end());
        std::sort(isect2.begin(), isect2.end());

        if (isect1[1]<isect2[0] || isect2[1]<isect1[0]) return false;
        return true;

	}

    /**
     * Check if an axis-aliged bounding box (AABB) intersects a triangle
     *
     * Based on code develop by Moller: http://fileadmin.cs.lth.se/cs/personal/tomas_akenine-moller/code/tribox3.txt
     * and the article "A Fast Triangle-Triangle Intersection Test", SIGGRAPH '05 ACM, Art.8, 2005:
     * http://fileadmin.cs.lth.se/cs/personal/tomas_akenine-moller/code/tribox_tam.pdf
     *
     * @return bool if the triangle overlaps a box
     * @param rLowPoint first corner of the box
     * @param rHighPoint second corner of the box
     */
    bool HasIntersection( const Point& rLowPoint, const Point& rHighPoint) override
    {
        Point box_center;
        Point box_half_size;

        box_center[0] = 0.5 * (rLowPoint[0] + rHighPoint[0]);
        box_center[1] = 0.5 * (rLowPoint[1] + rHighPoint[1]);
        box_center[2] = 0.5 * (rLowPoint[2] + rHighPoint[2]);

        box_half_size[0] = 0.5 * std::abs(rHighPoint[0] - rLowPoint[0]);
        box_half_size[1] = 0.5 * std::abs(rHighPoint[1] - rLowPoint[1]);
        box_half_size[2] = 0.5 * std::abs(rHighPoint[2] - rLowPoint[2]);

        return TriBoxOverlap(box_center, box_half_size);
    }

    /// Quality functions

    /** Calculates the inradius to circumradius quality metric.
     * Calculates the inradius to circumradius quality metric.
     * This metric is bounded by the interval (0,1) being:
     *  1 -> Optimal value
     *  0 -> Worst value
     *
     * @formulae $$ \frac{r_K}{\rho_K} $$
     *
     * @return The inradius to circumradius quality metric.
     */
    double InradiusToCircumradiusQuality() const override
    {
        constexpr double normFactor = 1.0;

        const double a = MathUtils<double>::Norm3(this->GetPoint(0)-this->GetPoint(1));
        const double b = MathUtils<double>::Norm3(this->GetPoint(1)-this->GetPoint(2));
        const double c = MathUtils<double>::Norm3(this->GetPoint(2)-this->GetPoint(0));

        return normFactor * CalculateInradius(a,b,c) / CalculateCircumradius(a,b,c);
    };

    /** Calculates the inradius to longest edge quality metric.
     * Calculates the inradius to longest edge quality metric.
     * This metric is bounded by the interval (0,1) being:
     *  1 -> Optimal value
     *  0 -> Worst value
     *
     * @formulae $$ \frac{r_K}{h_{max}} $$
     *
     * @return The inradius to longest edge quality metric.
     */
    double InradiusToLongestEdgeQuality() const override
    {
        constexpr double normFactor = 1.0; // TODO: This normalization coeficient is not correct.

        const array_1d<double, 3> a = this->GetPoint(0) - this->GetPoint(1);
        const array_1d<double, 3> b = this->GetPoint(1) - this->GetPoint(2);
        const array_1d<double, 3> c = this->GetPoint(2) - this->GetPoint(0);

        const double sa = (a[0]*a[0])+(a[1]*a[1])+(a[2]*a[2]);
        const double sb = (b[0]*b[0])+(b[1]*b[1])+(b[2]*b[2]);
        const double sc = (c[0]*c[0])+(c[1]*c[1])+(c[2]*c[2]);

        return normFactor * CalculateInradius(std::sqrt(sa),std::sqrt(sb),std::sqrt(sc)) / CalculateMaxEdgeLength(sa,sb,sc);
    }

    /** Calculates the area to edge length quality metric.
     * Calculates the area to edge length quality metric using the
     * sum of the distances squares.
     * This metric is bounded by the interval (0,1) being:
     *  1 -> Optimal value
     *  0 -> Worst value
     *
     * @formulae $$ \frac{h_{min}}{h_{max}} $$
     *
     * @return The Inradius to Circumradius Quality metric.
     */
    double AreaToEdgeLengthRatio() const override
    {
        constexpr double normFactor = 1.0;

        const array_1d<double, 3> a = this->GetPoint(0) - this->GetPoint(1);
        const array_1d<double, 3> b = this->GetPoint(1) - this->GetPoint(2);
        const array_1d<double, 3> c = this->GetPoint(2) - this->GetPoint(0);

        const double sa = (a[0]*a[0])+(a[1]*a[1])+(a[2]*a[2]);
        const double sb = (b[0]*b[0])+(b[1]*b[1])+(b[2]*b[2]);
        const double sc = (c[0]*c[0])+(c[1]*c[1])+(c[2]*c[2]);

        return normFactor * Area() / (sa+sb+sc);
    }

    /** Calculates the shortest altitude to edge length quality metric.
     * Calculates the shortest altitude to edge length quality metric.
     * This metric is bounded by the interval (0,1) being:
     *  1 -> Optimal value
     *  0 -> Worst value
     * -1 -> Optimal value with inverted volume.
     *
     * @formulae $$ \frac{h_{min}}{h_{max}} $$
     *
     * @return The shortest altitude to edge length quality metric.
     */
    double ShortestAltitudeToEdgeLengthRatio() const override {
      constexpr double normFactor = 1.0;

      const array_1d<double, 3> a = this->GetPoint(0) - this->GetPoint(1);
      const array_1d<double, 3> b = this->GetPoint(1) - this->GetPoint(2);
      const array_1d<double, 3> c = this->GetPoint(2) - this->GetPoint(0);

      const double sa = (a[0]*a[0])+(a[1]*a[1])+(a[2]*a[2]);
      const double sb = (b[0]*b[0])+(b[1]*b[1])+(b[2]*b[2]);
      const double sc = (c[0]*c[0])+(c[1]*c[1])+(c[2]*c[2]);

      // Shortest altitude is the one intersecting the largest base.
      double base = CalculateMaxEdgeLength(sa,sb,sc);

      return normFactor * (Area() * 2 / base ) / std::sqrt(sa+sb+sc);
    }

    /** Calculates the area to edge length quality metric.
     * Calculates the area to edge length quality metric using the
     * square of the sum of the distances.
     * This metric is bounded by the interval (0,1) being:
     *  1 -> Optimal value
     *  0 -> Worst value
     *
     * @formulae $$ \frac{h_{min}}{h_{max}} $$
     *
     * @return The Inradius to Circumradius Quality metric.
     */
    virtual double AreaToEdgeLengthSquareRatio() const {
      constexpr double normFactor = 1.0;

      const double a = MathUtils<double>::Norm3(this->GetPoint(0)-this->GetPoint(1));
      const double b = MathUtils<double>::Norm3(this->GetPoint(1)-this->GetPoint(2));
      const double c = MathUtils<double>::Norm3(this->GetPoint(2)-this->GetPoint(0));

      return normFactor * Area() / std::pow(a+b+c, 2);
    }

    /** Calculates the shortest altitude to longest edge.
     * Calculates the shortest altitude to longest edge.
     * This metric is bounded by the interval (0,1) being:
     *  1 -> Optimal value
     *  0 -> Worst value
     * -1 -> Optimal value with inverted volume.
     *
     * @formulae $$ \frac{h_{min}}{h_{max}} $$
     *
     * @return The shortest altitude to edge length quality metric.
     */
    virtual double ShortestAltitudeToLongestEdge() const
    {
        constexpr double normFactor = 1.0;

        const array_1d<double, 3> a = this->GetPoint(0) - this->GetPoint(1);
        const array_1d<double, 3> b = this->GetPoint(1) - this->GetPoint(2);
        const array_1d<double, 3> c = this->GetPoint(2) - this->GetPoint(0);

        const double sa = (a[0]*a[0])+(a[1]*a[1])+(a[2]*a[2]);
        const double sb = (b[0]*b[0])+(b[1]*b[1])+(b[2]*b[2]);
        const double sc = (c[0]*c[0])+(c[1]*c[1])+(c[2]*c[2]);

        // Shortest altitude is the one intersecting the largest base (or edge).
        const double base = CalculateMaxEdgeLength(sa,sb,sc);

        return normFactor * (Area() * 2 / base ) / base;
    }

    /**
     * It computes the area normal of the geometry
     * @param rPointLocalCoordinates Local coordinates of the point
     * in where the area normal is to be computed
     * @return The area normal in the given point
     */
    array_1d<double, 3> AreaNormal(const CoordinatesArrayType& rPointLocalCoordinates) const override
    {
        const array_1d<double, 3> tangent_xi  = this->GetPoint(1) - this->GetPoint(0);
        const array_1d<double, 3> tangent_eta = this->GetPoint(2) - this->GetPoint(0);

        array_1d<double, 3> normal;
        MathUtils<double>::CrossProduct(normal, tangent_xi, tangent_eta);

        return 0.5 * normal;
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
        PointLocalCoordinatesImplementation( rResult, rPoint, true );

        if ( (rResult[0] >= (0.0-Tolerance)) && (rResult[0] <= (1.0+Tolerance)) )
        {
            if ( (rResult[1] >= (0.0-Tolerance)) && (rResult[1] <= (1.0+Tolerance)) )
            {
                if ( (rResult[0] + rResult[1]) <= (1.0+Tolerance) )
                {
                    return true;
                }
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
    JacobiansType& Jacobian( JacobiansType& rResult,
                                     IntegrationMethod ThisMethod ) const override
    {
        Matrix jacobian( 3, 2 );
        jacobian( 0, 0 ) = -( BaseType::GetPoint( 0 ).X() ) + ( BaseType::GetPoint( 1 ).X() ); //on the Gauss points (J is constant at each element)
        jacobian( 1, 0 ) = -( BaseType::GetPoint( 0 ).Y() ) + ( BaseType::GetPoint( 1 ).Y() );
        jacobian( 2, 0 ) = -( BaseType::GetPoint( 0 ).Z() ) + ( BaseType::GetPoint( 1 ).Z() );
        jacobian( 0, 1 ) = -( BaseType::GetPoint( 0 ).X() ) + ( BaseType::GetPoint( 2 ).X() );
        jacobian( 1, 1 ) = -( BaseType::GetPoint( 0 ).Y() ) + ( BaseType::GetPoint( 2 ).Y() );
        jacobian( 2, 1 ) = -( BaseType::GetPoint( 0 ).Z() ) + ( BaseType::GetPoint( 2 ).Z() );

        if ( rResult.size() != this->IntegrationPointsNumber( ThisMethod ) )
        {
            // KLUDGE: While there is a bug in ublas vector resize, I have to put this beside resizing!!
            JacobiansType temp( this->IntegrationPointsNumber( ThisMethod ) );
            rResult.swap( temp );
        }

        std::fill( rResult.begin(), rResult.end(), jacobian );

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
    JacobiansType& Jacobian( JacobiansType& rResult,
                                     IntegrationMethod ThisMethod,
				     Matrix & DeltaPosition ) const override
    {
        Matrix jacobian( 3, 2 );
        jacobian( 0, 0 ) = -( BaseType::GetPoint( 0 ).X() - DeltaPosition(0,0) ) + ( BaseType::GetPoint( 1 ).X() - DeltaPosition(1,0) ); //on the Gauss points (J is constant at each element)
        jacobian( 1, 0 ) = -( BaseType::GetPoint( 0 ).Y() - DeltaPosition(0,1) ) + ( BaseType::GetPoint( 1 ).Y() - DeltaPosition(1,1) );
        jacobian( 2, 0 ) = -( BaseType::GetPoint( 0 ).Z() - DeltaPosition(0,2) ) + ( BaseType::GetPoint( 1 ).Z() - DeltaPosition(1,2) );
        jacobian( 0, 1 ) = -( BaseType::GetPoint( 0 ).X() - DeltaPosition(0,0) ) + ( BaseType::GetPoint( 2 ).X() - DeltaPosition(2,0) );
        jacobian( 1, 1 ) = -( BaseType::GetPoint( 0 ).Y() - DeltaPosition(0,1) ) + ( BaseType::GetPoint( 2 ).Y() - DeltaPosition(2,1) );
        jacobian( 2, 1 ) = -( BaseType::GetPoint( 0 ).Z() - DeltaPosition(0,2) ) + ( BaseType::GetPoint( 2 ).Z() - DeltaPosition(2,2) );

        if ( rResult.size() != this->IntegrationPointsNumber( ThisMethod ) )
        {
            // KLUDGE: While there is a bug in ublas vector resize, I have to put this beside resizing!!
            JacobiansType temp( this->IntegrationPointsNumber( ThisMethod ) );
            rResult.swap( temp );
        }

        std::fill( rResult.begin(), rResult.end(), jacobian );

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
    Matrix& Jacobian( Matrix& rResult,
                              IndexType IntegrationPointIndex,
                              IntegrationMethod ThisMethod ) const override
    {
        rResult.resize( 3, 2,false );
        rResult( 0, 0 ) = -( BaseType::GetPoint( 0 ).X() ) + ( BaseType::GetPoint( 1 ).X() ); //on the Gauss points (J is constant at each element)
        rResult( 1, 0 ) = -( BaseType::GetPoint( 0 ).Y() ) + ( BaseType::GetPoint( 1 ).Y() );
        rResult( 2, 0 ) = -( BaseType::GetPoint( 0 ).Z() ) + ( BaseType::GetPoint( 1 ).Z() );
        rResult( 0, 1 ) = -( BaseType::GetPoint( 0 ).X() ) + ( BaseType::GetPoint( 2 ).X() );
        rResult( 1, 1 ) = -( BaseType::GetPoint( 0 ).Y() ) + ( BaseType::GetPoint( 2 ).Y() );
        rResult( 2, 1 ) = -( BaseType::GetPoint( 0 ).Z() ) + ( BaseType::GetPoint( 2 ).Z() );
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
        rResult.resize( 3, 2 ,false);
        rResult( 0, 0 ) = -( BaseType::GetPoint( 0 ).X() ) + ( BaseType::GetPoint( 1 ).X() );
        rResult( 1, 0 ) = -( BaseType::GetPoint( 0 ).Y() ) + ( BaseType::GetPoint( 1 ).Y() );
        rResult( 2, 0 ) = -( BaseType::GetPoint( 0 ).Z() ) + ( BaseType::GetPoint( 1 ).Z() );
        rResult( 0, 1 ) = -( BaseType::GetPoint( 0 ).X() ) + ( BaseType::GetPoint( 2 ).X() );
        rResult( 1, 1 ) = -( BaseType::GetPoint( 0 ).Y() ) + ( BaseType::GetPoint( 2 ).Y() );
        rResult( 2, 1 ) = -( BaseType::GetPoint( 0 ).Z() ) + ( BaseType::GetPoint( 2 ).Z() );
        return rResult;
    }

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
    Vector& DeterminantOfJacobian( Vector& rResult, IntegrationMethod ThisMethod ) const override
    {
        const unsigned int integration_points_number = msGeometryData.IntegrationPointsNumber( ThisMethod );
        if(rResult.size() != integration_points_number)
        {
            rResult.resize(integration_points_number,false);
        }

        const double detJ = 2.0*(this->Area());

        for ( unsigned int pnt = 0; pnt < integration_points_number; pnt++ )
        {
            rResult[pnt] = detJ;
        }
        return rResult;
    }

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
    double DeterminantOfJacobian( IndexType IntegrationPoint,
                                          IntegrationMethod ThisMethod ) const override
    {
        return 2.0*(this->Area());
    }

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
     */
    double DeterminantOfJacobian( const CoordinatesArrayType& rPoint ) const override
    {
        return 2.0*(this->Area());
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
        KRATOS_ERROR << "Triangle3D::InverseOfJacobian" << "Jacobian is not square" << std::endl;
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
        KRATOS_ERROR << "Triangle3D::InverseOfJacobian" << "Jacobian is not square" << std::endl;
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
        KRATOS_ERROR << "Triangle3D::InverseOfJacobian" << "Jacobian is not square" << std::endl;
        return rResult;
    }

    /** EdgesNumber
    @return SizeType containes number of this geometry edges.
    */
    SizeType EdgesNumber() const override
    {
        return 3;
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

        edges.push_back( Kratos::make_shared<EdgeType>( this->pGetPoint( 0 ), this->pGetPoint( 1 ) ) );
        edges.push_back( Kratos::make_shared<EdgeType>( this->pGetPoint( 1 ), this->pGetPoint( 2 ) ) );
        edges.push_back( Kratos::make_shared<EdgeType>( this->pGetPoint( 2 ), this->pGetPoint( 0 ) ) );
        return edges;
    }


    //Connectivities of faces required
    void NumberNodesInFaces (DenseVector<unsigned int>& NumberNodesInFaces) const override
    {
        if(NumberNodesInFaces.size() != 3 )
            NumberNodesInFaces.resize(3,false);
        // Linear Triangles have elements of 2 nodes as faces
        NumberNodesInFaces[0]=2;
        NumberNodesInFaces[1]=2;
        NumberNodesInFaces[2]=2;

    }

    void NodesInFaces (DenseMatrix<unsigned int>& NodesInFaces) const override
    {
        // faces in columns
        if(NodesInFaces.size1() != 3 || NodesInFaces.size2() != 3)
            NodesInFaces.resize(3,3,false);

        //face 1
        NodesInFaces(0,0)=0;//contrary node to the face
        NodesInFaces(1,0)=1;
        NodesInFaces(2,0)=2;
        //face 2
        NodesInFaces(0,1)=1;//contrary node to the face
        NodesInFaces(1,1)=2;
        NodesInFaces(2,1)=0;
        //face 3
        NodesInFaces(0,2)=2;//contrary node to the face
        NodesInFaces(1,2)=0;
        NodesInFaces(2,2)=1;

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
            return( 1.0 -rPoint[0] - rPoint[1] );
        case 1:
            return( rPoint[0] );
        case 2:
            return( rPoint[1] );
        default:
            KRATOS_ERROR << "Wrong index of shape function!" << *this << std::endl;
        }

        return 0;
    }

    /** This method gives all non-zero shape functions values
     * evaluated at the rCoordinates provided
     * \note There is no control if the return vector is empty or not!
     * @return Vector of values of shape functions \f$ F_{i} \f$
     * where i is the shape function index (for NURBS it is the index
     * of the local enumeration in the element).
     *
     * @see ShapeFunctionValue
     * @see ShapeFunctionsLocalGradients
     * @see ShapeFunctionLocalGradient
     */

    Vector& ShapeFunctionsValues (Vector &rResult, const CoordinatesArrayType& rCoordinates) const override
    {
        if(rResult.size() != 3)
        {
            rResult.resize(3,false);
        }

        rResult[0] =  1.0 -rCoordinates[0] - rCoordinates[1];
        rResult[1] =  rCoordinates[0] ;
        rResult[2] =  rCoordinates[1] ;

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
            rResult[pnt].resize( 3, 2,false );

            for ( int i = 0; i < 3; i++ )
            {
                for ( int j = 0; j < 2; j++ )
                {
                    rResult[pnt]( i, j ) =
                        ( locG[pnt]( i, 0 ) * invJ[pnt]( j, 0 ) )
                        + ( locG[pnt]( i, 1 ) * invJ[pnt]( j, 1 ) );
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
        return "2 dimensional triangle with three nodes in 3D space";
    }

    /**
     * Print information about this object.
     * @param rOStream Stream to print into it.
     * @see PrintData()
     * @see Info()
     */
    void PrintInfo( std::ostream& rOStream ) const override
    {
        rOStream << "2 dimensional triangle with three nodes in 3D space";
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
        ShapeFunctionsGradientsType localGradients
        = CalculateShapeFunctionsIntegrationPointsLocalGradients( ThisMethod );
        const int integration_points_number
        = msGeometryData.IntegrationPointsNumber( ThisMethod );
        ShapeFunctionsGradientsType Result( integration_points_number );

        for ( int pnt = 0; pnt < integration_points_number; pnt++ )
        {
            Result[pnt] = localGradients[pnt];
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
        ShapeFunctionsGradientsType localGradients
        = CalculateShapeFunctionsIntegrationPointsLocalGradients( ThisMethod );
        const int integration_points_number
        = msGeometryData.IntegrationPointsNumber( ThisMethod );
        ShapeFunctionsGradientsType Result( integration_points_number );

        for ( int pnt = 0; pnt < integration_points_number; pnt++ )
        {
            Result[pnt] = localGradients[pnt];
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
    Matrix& ShapeFunctionsLocalGradients( Matrix& rResult,
            const CoordinatesArrayType& rPoint ) const override
    {
        rResult.resize( 3, 2 ,false);
        noalias( rResult ) = ZeroMatrix( 3, 2 );
        rResult( 0, 0 ) = -1.0;
        rResult( 0, 1 ) = -1.0;
        rResult( 1, 0 ) =  1.0;
        rResult( 1, 1 ) =  0.0;
        rResult( 2, 0 ) =  0.0;
        rResult( 2, 1 ) =  1.0;
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
    virtual Matrix& ShapeFunctionsGradients( Matrix& rResult, PointType& rPoint )
    {
        rResult.resize( 3, 2,false );
        noalias( rResult ) = ZeroMatrix( 3, 2 );
        rResult( 0, 0 ) = -1.0;
        rResult( 0, 1 ) = -1.0;
        rResult( 1, 0 ) =  1.0;
        rResult( 1, 1 ) =  0.0;
        rResult( 2, 0 ) =  0.0;
        rResult( 2, 1 ) =  1.0;
        return rResult;
    }

    /**
     * returns the second order derivatives of all shape functions
     * in given arbitrary pointers
     * @param rResult a third order tensor which contains the second derivatives
     * @param rPoint the given point the second order derivatives are calculated in
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

        rResult[0].resize( 2, 2 ,false);

        rResult[1].resize( 2, 2 ,false);
        rResult[2].resize( 2, 2 ,false);
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
    ShapeFunctionsThirdDerivativesType& ShapeFunctionsThirdDerivatives( ShapeFunctionsThirdDerivativesType& rResult, const CoordinatesArrayType& rPoint ) const override
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

        rResult[0][0].resize( 2, 2 ,false);

        rResult[0][1].resize( 2, 2 ,false);
        rResult[1][0].resize( 2, 2 ,false);
        rResult[1][1].resize( 2, 2 ,false);
        rResult[2][0].resize( 2, 2 ,false);
        rResult[2][1].resize( 2, 2 ,false);

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
     * There are no protected members in class Triangle3D3
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

    Triangle3D3(): BaseType( PointsArrayType(), &msGeometryData ) {}

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
     * Returns the local coordinates of a given arbitrary point
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
        BoundedMatrix<double, 3, 3> X;
        BoundedMatrix<double, 3, 2> DN;
        for(unsigned int i=0; i<this->size();i++)
        {
            X(0,i ) = this->GetPoint( i ).X();
            X(1,i ) = this->GetPoint( i ).Y();
            X(2,i ) = this->GetPoint( i ).Z();
        }

        static constexpr double MaxNormPointLocalCoordinates = 30.0;
        static constexpr std::size_t MaxIteratioNumberPointLocalCoordinates = 1000;
        static constexpr double MaxTolerancePointLocalCoordinates = 1.0e-8;

        BoundedMatrix<double, 2, 2> J = ZeroMatrix( 2, 2 );
        BoundedMatrix<double, 2, 2> invJ = ZeroMatrix( 2, 2 );

        //starting with xi = 0
        noalias(rResult) = ZeroVector( 3 );
        Vector delta_xi = ZeroVector( 2 );
        array_1d<double,3> current_global_coords;

        //Newton iteration:
        for ( std::size_t k = 0; k < MaxIteratioNumberPointLocalCoordinates; k++ )
        {
            noalias(current_global_coords) = ZeroVector( 3 );
            this->GlobalCoordinates( current_global_coords, rResult );

            noalias( current_global_coords ) = rPoint - current_global_coords;

            //derivatives of shape functions
            Matrix shape_functions_gradients;
            shape_functions_gradients = ShapeFunctionsLocalGradients(shape_functions_gradients, rResult );
            noalias(DN) = prod(X,shape_functions_gradients);

            noalias(J) = prod(trans(DN),DN);
            Vector res = prod(trans(DN),current_global_coords);

            //deteminant of Jacobian
            const double det_j = J( 0, 0 ) * J( 1, 1 ) - J( 0, 1 ) * J( 1, 0 );

            //filling matrix
            invJ( 0, 0 ) = ( J( 1, 1 ) ) / ( det_j );
            invJ( 1, 0 ) = -( J( 1, 0 ) ) / ( det_j );
            invJ( 0, 1 ) = -( J( 0, 1 ) ) / ( det_j );
            invJ( 1, 1 ) = ( J( 0, 0 ) ) / ( det_j );

            delta_xi( 0 ) = invJ( 0, 0 ) * res[0] + invJ( 0, 1 ) * res[1];
            delta_xi( 1 ) = invJ( 1, 0 ) * res[0] + invJ( 1, 1 ) * res[1];

            rResult[0] += delta_xi[0];
            rResult[1] += delta_xi[1];
            rResult[2] = 0.0;

            if ( k>0 && norm_2( delta_xi ) > MaxNormPointLocalCoordinates ) {
                KRATOS_WARNING_IF("Triangle3D3", IsInside == false) << "detJ =\t" << det_j << " DeltaX =\t" << delta_xi << " stopping calculation. Iteration:\t" << k << std::endl;
                break;
            }

            if ( norm_2( delta_xi ) < MaxTolerancePointLocalCoordinates ) {
                break;
            }
        }

        return( rResult );
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
        const int points_number = 3;
        //setting up return matrix
        Matrix shape_function_values( integration_points_number, points_number );
        //loop over all integration points

        for ( int pnt = 0; pnt < integration_points_number; pnt++ )
        {
            row( shape_function_values, pnt )[0] = 1.0
                                                   - integration_points[pnt].X()
                                                   - integration_points[pnt].Y();
            row( shape_function_values, pnt )[1] = integration_points[pnt].X();
            row( shape_function_values, pnt )[2] = integration_points[pnt].Y();
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
        //std::fill(d_shape_f_values.begin(), d_shape_f_values.end(), Matrix(4,2));
        //loop over all integration points

        for ( int pnt = 0; pnt < integration_points_number; pnt++ )
        {
            Matrix result( 3, 2 );
            result( 0, 0 ) = -1.0;
            result( 0, 1 ) = -1.0;
            result( 1, 0 ) =  1.0;
            result( 1, 1 ) =  0.0;
            result( 2, 0 ) =  0.0;
            result( 2, 1 ) =  1.0;
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
                Quadrature<TriangleGaussLegendreIntegrationPoints1, 2, IntegrationPoint<3> >::GenerateIntegrationPoints(),
                Quadrature<TriangleGaussLegendreIntegrationPoints2, 2, IntegrationPoint<3> >::GenerateIntegrationPoints(),
                Quadrature<TriangleGaussLegendreIntegrationPoints3, 2, IntegrationPoint<3> >::GenerateIntegrationPoints(),
                Quadrature<TriangleGaussLegendreIntegrationPoints4, 2, IntegrationPoint<3> >::GenerateIntegrationPoints(),
                Quadrature<TriangleGaussLegendreIntegrationPoints5, 2, IntegrationPoint<3> >::GenerateIntegrationPoints(),
                Quadrature<TriangleCollocationIntegrationPoints1, 2, IntegrationPoint<3> >::GenerateIntegrationPoints(),
                Quadrature<TriangleCollocationIntegrationPoints2, 2, IntegrationPoint<3> >::GenerateIntegrationPoints(),
                Quadrature<TriangleCollocationIntegrationPoints3, 2, IntegrationPoint<3> >::GenerateIntegrationPoints(),
                Quadrature<TriangleCollocationIntegrationPoints4, 2, IntegrationPoint<3> >::GenerateIntegrationPoints(),
                Quadrature<TriangleCollocationIntegrationPoints5, 2, IntegrationPoint<3> >::GenerateIntegrationPoints()
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
                Triangle3D3<TPointType>::CalculateShapeFunctionsIntegrationPointsValues(
                    GeometryData::GI_GAUSS_1 ),
                Triangle3D3<TPointType>::CalculateShapeFunctionsIntegrationPointsValues(
                    GeometryData::GI_GAUSS_2 ),
                Triangle3D3<TPointType>::CalculateShapeFunctionsIntegrationPointsValues(
                    GeometryData::GI_GAUSS_3 ),
                Triangle3D3<TPointType>::CalculateShapeFunctionsIntegrationPointsValues(
                    GeometryData::GI_GAUSS_4 ),
                Triangle3D3<TPointType>::CalculateShapeFunctionsIntegrationPointsValues(
                    GeometryData::GI_GAUSS_5 ),
                 Triangle3D3<TPointType>::CalculateShapeFunctionsIntegrationPointsValues(
                     GeometryData::GI_EXTENDED_GAUSS_1 ),
                 Triangle3D3<TPointType>::CalculateShapeFunctionsIntegrationPointsValues(
                     GeometryData::GI_EXTENDED_GAUSS_2 ),
                 Triangle3D3<TPointType>::CalculateShapeFunctionsIntegrationPointsValues(
                     GeometryData::GI_EXTENDED_GAUSS_3 ),
                 Triangle3D3<TPointType>::CalculateShapeFunctionsIntegrationPointsValues(
                     GeometryData::GI_EXTENDED_GAUSS_4 ),
                 Triangle3D3<TPointType>::CalculateShapeFunctionsIntegrationPointsValues(
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
                Triangle3D3<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients( GeometryData::GI_GAUSS_1 ),
                Triangle3D3<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients( GeometryData::GI_GAUSS_2 ),
                Triangle3D3<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients( GeometryData::GI_GAUSS_3 ),
                Triangle3D3<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients( GeometryData::GI_GAUSS_4 ),
                Triangle3D3<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients( GeometryData::GI_GAUSS_5 ),
                Triangle3D3<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients( GeometryData::GI_EXTENDED_GAUSS_1 ),
                Triangle3D3<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients( GeometryData::GI_EXTENDED_GAUSS_2 ),
                Triangle3D3<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients( GeometryData::GI_EXTENDED_GAUSS_3 ),
                Triangle3D3<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients( GeometryData::GI_EXTENDED_GAUSS_4 ),
                Triangle3D3<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients( GeometryData::GI_EXTENDED_GAUSS_5 )
            }
        };
        return shape_functions_local_gradients;
    }

    /** Implements the calculus of the minimum edge length
     * Implements the calculus of the minimum edge length given the length of the geometry edges.
     *
     * @param  a Squared Length of the edge a
     * @param  b Squared Length of the edge b
     * @param  c Squared Length of the edge c
     *
     * @return   The minimum edge length of the geometry with edges a,b,c
     */
    inline double CalculateMinEdgeLength(const double a, const double b, const double c) const {
      return std::sqrt(std::min({a, b, c}));
    }

    /** Implements the calculus of the maximum edge length
     * Implements the calculus of the maximum edge length given the length of the geometry edges.
     *
     * @param  a Squared Length of the edge a
     * @param  b Squared Length of the edge b
     * @param  c Squared Length of the edge c
     *
     * @return   The maximum edge length of the geometry with edges a,b,c
     */
    inline double CalculateMaxEdgeLength(const double a, const double b, const double c) const {
      return std::sqrt(std::max({a, b, c}));
    }

    /** Implements the calculus of the average edge length
     * Implements the calculus of the average edge length given the length of the geometry edges.
     *
     * @param  a Length of the edge a
     * @param  b Length of the edge b
     * @param  c Length of the edge c
     *
     * @return   The average edge length of the geometry with edges a,b,c
     */
    inline double CalculateAvgEdgeLength(const double a, const double b, const double c) const {
      constexpr double onethird = 1.0 / 3.0;
      return (a+b+c) * onethird;
    }

    /** Implements the calculus of the circumradius
     * Implements the calculus of the circumradius given the length of the geometry edges.
     *
     * @param  a Length of the edge a
     * @param  b Length of the edge b
     * @param  c Length of the edge c
     *
     * @return   The circumradius of the geometry with edges a,b,c
     */
    inline double CalculateCircumradius(const double a, const double b, const double c) const {
      return (a*b*c) / std::sqrt((a+b+c) * (b+c-a) * (c+a-b) * (a+b-c));
    }

    /** Implements the calculus of the inradius
     * Implements the calculus of the inradius given the length of the geometry edges.
     *
     * @param  a Length of the edge a
     * @param  b Length of the edge b
     * @param  c Length of the edge c
     *
     * @return   The inradius of the geometry with edges a,b,c
     */
    inline double CalculateInradius(const double a, const double b, const double c) const {
      return 0.5 * std::sqrt((b+c-a) * (c+a-b) * (a+b-c) / (a+b+c));
    }

	bool ComputeIntervals(double& VV0,
		double& VV1,
		double& VV2,
		double& D0,
		double& D1,
		double& D2,
		double& A,
		double& B,
		double& C,
		double& X0,
		double& X1
		)
	{
		double D0D1 = D0*D1;
		double D0D2 = D0*D2;

		if (D0D1>0.00)
		{
			// here we know that D0D2<=0.0 //
			// that is D0, D1 are on the same side, D2 on the other or on the plane //
			A = VV2;
			B = (VV0 - VV2)*D2;
			C = (VV1 - VV2)*D2;
			X0 = D2 - D0;
			X1 = D2 - D1;
		}
		else if (D0D2>0.00)
		{
			// here we know that d0d1<=0.0 //
			A = VV1;
			B = (VV0 - VV1)*D1;
			C = (VV2 - VV1)*D1;
			X0 = D1 - D0;
			X1 = D1 - D2;
		}
		else if (D1*D2>0.00 || D0 != 0.00)
		{
			// here we know that d0d1<=0.0 or that D0!=0.0 //
			A = VV0;
			B = (VV1 - VV0)*D0;
			C = (VV2 - VV0)*D0;
			X0 = D0 - D1;
			X1 = D0 - D2;
		}
		else if (D1 != 0.00)
		{
			A = VV1;
			B = (VV0 - VV1)*D1;
			C = (VV2 - VV1)*D1;
			X0 = D1 - D0;
			X1 = D1 - D2;
		}
		else if (D2 != 0.00)
		{
			A = VV2;
			B = (VV0 - VV2)*D2;
			C = (VV1 - VV2)*D2;
			X0 = D2 - D0;
			X1 = D2 - D1;
		}
		else
		{
			///Triangles are coplanar
			return true;
		}

		return false;

	}

	bool CoplanarIntersectionCheck(const array_1d<double, 3>& N,
		const GeometryType& OtherTriangle)
	{
		array_1d<double, 3 > A;
		short i0, i1;

		// first project onto an axis-aligned plane, that maximizes the area //
		// of the triangles, compute indices: i0,i1. //
		A[0] = std::abs(N[0]);
		A[1] = std::abs(N[1]);
		A[2] = std::abs(N[2]);
		if (A[0]>A[1])
		{
			if (A[0]>A[2])
			{
				i0 = 1;      // A[0] is greatest //
				i1 = 2;
			}
			else
			{
				i0 = 0;      // A[2] is greatest //
				i1 = 1;
			}
		}
		else   // A[0]<=A[1] //
		{
			if (A[2]>A[1])
			{
				i0 = 0;      // A[2] is greatest //
				i1 = 1;
			}
			else
			{
				i0 = 0;      // A[1] is greatest //
				i1 = 2;
			}
		}

		// test all edges of triangle 1 against the edges of triangle 2 //
		if (EdgeToTriangleEdgesCheck(i0, i1, this->GetPoint(0), this->GetPoint(1), OtherTriangle[0], OtherTriangle[1], OtherTriangle[2]) == true) return true;

		if (EdgeToTriangleEdgesCheck(i0, i1, this->GetPoint(1), this->GetPoint(2), OtherTriangle[0], OtherTriangle[1], OtherTriangle[2]) == true) return true;

		if (EdgeToTriangleEdgesCheck(i0, i1, this->GetPoint(2), this->GetPoint(0), OtherTriangle[0], OtherTriangle[1], OtherTriangle[2]) == true) return true;

		// finally, test if tri1 is totally contained in tri2 or vice versa //
		array_1d<double, 3> local_coordinates;
		// TODO: I should add the const to the is inside method in all geometries. Pooyan.
		if (const_cast<GeometryType&>(OtherTriangle).IsInside(this->GetPoint(0), local_coordinates) == true) return true;
		if (IsInside(OtherTriangle[0], local_coordinates) == true) return true;

		return false;
	}

	bool EdgeToTriangleEdgesCheck(const short& i0,
		const short& i1,
		const Point& V0,
		const Point& V1,
		const Point&U0,
		const Point&U1,
		const Point&U2)
	{

		double Ax, Ay, Bx, By, Cx, Cy, e, d, f;
		Ax = V1[i0] - V0[i0];
		Ay = V1[i1] - V0[i1];
		// test edge U0,U1 against V0,V1 //

		//std::cout<< "Proof One B " << std::endl;
		if (EdgeToEdgeIntersectionCheck(Ax, Ay, Bx, By, Cx, Cy, e, d, f, i0, i1, V0, U0, U1) == true) return true;
		// test edge U1,U2 against V0,V1 //
		//std::cout<< "Proof Two B " << std::endl;
		if (EdgeToEdgeIntersectionCheck(Ax, Ay, Bx, By, Cx, Cy, e, d, f, i0, i1, V0, U1, U2) == true) return true;
		// test edge U2,U1 against V0,V1 //
		if (EdgeToEdgeIntersectionCheck(Ax, Ay, Bx, By, Cx, Cy, e, d, f, i0, i1, V0, U2, U0) == true) return true;

		return false;
	}

	//   this edge to edge test is based on Franlin Antonio's gem:
	//   "Faster Line Segment Intersection", in Graphics Gems III,
	//   pp. 199-202
	bool EdgeToEdgeIntersectionCheck(double& Ax,
		double& Ay,
		double& Bx,
		double& By,
		double& Cx,
		double& Cy,
		double& e,
		double& d,
		double& f,
		const short& i0,
		const short& i1,
		const Point& V0,
		const Point& U0,
		const Point& U1)
	{
		Bx = U0[i0] - U1[i0];
		By = U0[i1] - U1[i1];
		Cx = V0[i0] - U0[i0];
		Cy = V0[i1] - U0[i1];
		f = Ay*Bx - Ax*By;
		d = By*Cx - Bx*Cy;

		if (std::abs(f)<1E-10) f = 0.00;
		if (std::abs(d)<1E-10) d = 0.00;


		if ((f>0.00 && d >= 0.00 && d <= f) || (f<0.00 && d <= 0.00 && d >= f))
		{
			e = Ax*Cy - Ay*Cx;

			if (f>0.00)
			{
				if (e >= 0.00 && e <= f) return true;
			}
			else
			{
				if (e <= 0.00 && e >= f) return true;
			}
		}
		return false;
	}


    /**
     * @see HasIntersection
     * use separating axis theorem to test overlap between triangle and box
     * need to test for overlap in these directions:
     * 1) the {x,y,(z)}-directions
     * 2) normal of the triangle
     * 3) crossproduct (edge from tri, {x,y,z}-direction) gives 3x3=9 more tests
     */
    inline bool TriBoxOverlap(Point& rBoxCenter, Point& rBoxHalfSize)
    {
        double abs_ex, abs_ey, abs_ez, distance;
        array_1d<double,3 > vert0, vert1, vert2;
        array_1d<double,3 > edge0, edge1, edge2, normal;
        std::pair<double, double> min_max;

        // move everything so that the boxcenter is in (0,0,0)
        noalias(vert0) = this->GetPoint(0) - rBoxCenter;
        noalias(vert1) = this->GetPoint(1) - rBoxCenter;
        noalias(vert2) = this->GetPoint(2) - rBoxCenter;

        // compute triangle edges
        noalias(edge0) = vert1 - vert0;
        noalias(edge1) = vert2 - vert1;
        noalias(edge2) = vert0 - vert2;

        // Bullet 3:
        // test the 12 tests first (this was faster)
        abs_ex = std::abs(edge0[0]);
        abs_ey = std::abs(edge0[1]);
        abs_ez = std::abs(edge0[2]);
        if (AxisTestX(edge0[1],edge0[2],abs_ey,abs_ez,vert0,vert2,rBoxHalfSize)) return false;
        if (AxisTestY(edge0[0],edge0[2],abs_ex,abs_ez,vert0,vert2,rBoxHalfSize)) return false;
        if (AxisTestZ(edge0[0],edge0[1],abs_ex,abs_ey,vert0,vert2,rBoxHalfSize)) return false;

        abs_ex = std::abs(edge1[0]);
        abs_ey = std::abs(edge1[1]);
        abs_ez = std::abs(edge1[2]);
        if (AxisTestX(edge1[1],edge1[2],abs_ey,abs_ez,vert1,vert0,rBoxHalfSize)) return false;
        if (AxisTestY(edge1[0],edge1[2],abs_ex,abs_ez,vert1,vert0,rBoxHalfSize)) return false;
        if (AxisTestZ(edge1[0],edge1[1],abs_ex,abs_ey,vert1,vert0,rBoxHalfSize)) return false;

        abs_ex = std::abs(edge2[0]);
        abs_ey = std::abs(edge2[1]);
        abs_ez = std::abs(edge2[2]);
        if (AxisTestX(edge2[1],edge2[2],abs_ey,abs_ez,vert2,vert1,rBoxHalfSize)) return false;
        if (AxisTestY(edge2[0],edge2[2],abs_ex,abs_ez,vert2,vert1,rBoxHalfSize)) return false;
        if (AxisTestZ(edge2[0],edge2[1],abs_ex,abs_ey,vert2,vert1,rBoxHalfSize)) return false;

        // Bullet 1:
        //  first test overlap in the {x,y,z}-directions
        //  find min, max of the triangle for each direction, and test for
        //  overlap in that direction -- this is equivalent to testing a minimal
        //  AABB around the triangle against the AABB

        // test in X-direction
        min_max = std::minmax({vert0[0], vert1[0], vert2[0]});
        if(min_max.first>rBoxHalfSize[0] || min_max.second<-rBoxHalfSize[0]) return false;

        // test in Y-direction
        min_max = std::minmax({vert0[0], vert1[0], vert2[0]});
        if(min_max.first>rBoxHalfSize[1] || min_max.second<-rBoxHalfSize[1]) return false;

        // test in Z-direction
        min_max = std::minmax({vert0[0], vert1[0], vert2[0]});
        if(min_max.first>rBoxHalfSize[2] || min_max.second<-rBoxHalfSize[2]) return false;

        // Bullet 2:
        //  test if the box intersects the plane of the triangle
        //  compute plane equation of triangle: normal*x+distance=0
        MathUtils<double>::CrossProduct(normal, edge0, edge1);
        distance = -inner_prod(normal, vert0);
        if(!PlaneBoxOverlap(normal, distance, rBoxHalfSize)) return false;

        return true;  // box and triangle overlaps
    }

    /**
     * Check if a plane intersects a box
     * @see TriBoxOverlap
     *
     * @return bool intersection flagg
     * @param rNormal the plane normal
     * @param rDist   distance to origin
     * @param rMaxBox box corner from the origin
     *
     * plane equation: rNormal*x+rDist=0
     */
    bool PlaneBoxOverlap(const array_1d<double,3>& rNormal, const double& rDist, const array_1d<double,3>& rMaxBox)
    {
        array_1d<double,3> vmin, vmax;
        for(int q = 0; q < 3; q++)
        {
            if(rNormal[q] > 0.00)
            {
                vmin[q] = -rMaxBox[q];
                vmax[q] =  rMaxBox[q];
            }
            else
            {
                vmin[q] =  rMaxBox[q];
                vmax[q] = -rMaxBox[q];
            }
        }
        if(inner_prod(rNormal, vmin) + rDist >  0.00) return false;
        if(inner_prod(rNormal, vmax) + rDist >= 0.00) return true;

        return false;
    }

    /** AxisTestX
     * This method returns true if there is a separating axis
     *
     * @param rEdgeY, rEdgeZ: i-edge corrdinates
     * @param rAbsEdgeY, rAbsEdgeZ: i-edge abs coordinates
     * @param rVertA: i   vertex
     * @param rVertB: i+1 vertex (omitted, proj_a = proj_b)
     * @param rVertC: i+2 vertex
     * @param rBoxHalfSize
     */
    bool AxisTestX(double& rEdgeY, double& rEdgeZ,
                   double& rAbsEdgeY, double& rAbsEdgeZ,
                   array_1d<double,3>& rVertA,
                   array_1d<double,3>& rVertC,
                   Point& rBoxHalfSize)
    {
        double proj_a, proj_c, rad;
        proj_a = rEdgeY*rVertA[2] - rEdgeZ*rVertA[1];
        proj_c = rEdgeY*rVertC[2] - rEdgeZ*rVertC[1];
        std::pair<double, double> min_max = std::minmax(proj_a, proj_c);

        rad = rAbsEdgeZ*rBoxHalfSize[1] + rAbsEdgeY*rBoxHalfSize[2];

        if(min_max.first>rad || min_max.second<-rad) return true;
        else return false;
    }

    /** AxisTestY
     * This method returns true if there is a separating axis
     *
     * @param rEdgeX, rEdgeZ: i-edge corrdinates
     * @param rAbsEdgeX, rAbsEdgeZ: i-edge fabs coordinates
     * @param rVertA: i   vertex
     * @param rVertB: i+1 vertex (omitted, proj_a = proj_b)
     * @param rVertC: i+2 vertex
     * @param rBoxHalfSize
     */
    bool AxisTestY(double& rEdgeX, double& rEdgeZ,
                   double& rAbsEdgeX, double& rAbsEdgeZ,
                   array_1d<double,3>& rVertA,
                   array_1d<double,3>& rVertC,
                   Point& rBoxHalfSize)
    {
        double proj_a, proj_c, rad;
        proj_a = rEdgeZ*rVertA[0] - rEdgeX*rVertA[2];
        proj_c = rEdgeZ*rVertC[0] - rEdgeX*rVertC[2];
        std::pair<double, double> min_max = std::minmax(proj_a, proj_c);

        rad = rAbsEdgeZ*rBoxHalfSize[0] + rAbsEdgeX*rBoxHalfSize[2];

        if(min_max.first>rad || min_max.second<-rad) return true;
        else return false;
    }

    /** AxisTestZ
     * This method returns true if there is a separating axis
     *
     * @param rEdgeX, rEdgeY: i-edge corrdinates
     * @param rAbsEdgeX, rAbsEdgeY: i-edge fabs coordinates
     * @param rVertA: i   vertex
     * @param rVertB: i+1 vertex (omitted, proj_a = proj_b)
     * @param rVertC: i+2 vertex
     * @param rBoxHalfSize
     */
    bool AxisTestZ(double& rEdgeX, double& rEdgeY,
                   double& rAbsEdgeX, double& rAbsEdgeY,
                   array_1d<double,3>& rVertA,
                   array_1d<double,3>& rVertC,
                   Point& rBoxHalfSize)
    {
        double proj_a, proj_c, rad;
        proj_a = rEdgeX*rVertA[1] - rEdgeY*rVertA[0];
        proj_c = rEdgeX*rVertC[1] - rEdgeY*rVertC[0];
        std::pair<double, double> min_max = std::minmax(proj_a, proj_c);

        rad = rAbsEdgeY*rBoxHalfSize[0] + rAbsEdgeX*rBoxHalfSize[1];

        if(min_max.first>rad || min_max.second<-rad) return true;
        else return false;
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

    template<class TOtherPointType> friend class Triangle3D3;

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
    Triangle3D3<TPointType>& rThis );
/**
 * output stream functions
 */
template<class TPointType> inline std::ostream& operator << (
    std::ostream& rOStream,
    const Triangle3D3<TPointType>& rThis )
{
    rThis.PrintInfo( rOStream );
    rOStream << std::endl;
    rThis.PrintData( rOStream );
    return rOStream;
}

///@}

template<class TPointType> const
GeometryData Triangle3D3<TPointType>::msGeometryData(
    2, 3, 2,
    GeometryData::GI_GAUSS_1,
    Triangle3D3<TPointType>::AllIntegrationPoints(),
    Triangle3D3<TPointType>::AllShapeFunctionsValues(),
    AllShapeFunctionsLocalGradients()
);
}// namespace Kratos.

#endif // KRATOS_QUADRILATERAL_3D_4_H_INCLUDED  defined
