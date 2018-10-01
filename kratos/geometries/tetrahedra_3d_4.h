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

#if !defined(KRATOS_TETRAHEDRA_3D_4_H_INCLUDED )
#define  KRATOS_TETRAHEDRA_3D_4_H_INCLUDED

// System includes

// External includes

// Project includes
#include "geometries/triangle_3d_3.h"
#include "integration/tetrahedron_gauss_legendre_integration_points.h"

#include "geometries/plane.h"

namespace Kratos
{
/**
 * @class Tetrahedra3D4
 * @ingroup KratosCore
 * @brief A four node tetrahedra geometry with linear shape functions
 * @details The node ordering corresponds with:       
 *                             v
 *                            .
 *                          ,/
 *                         /
 *                      2                                                            
 *                    ,/|`\                                                       
 *                  ,/  |  `\                         
 *                ,/    '.   `\                      
 *              ,/       |     `\                 
 *            ,/         |       `\                
 *           0-----------'.--------1 --> u        
 *            `\.         |      ,/               
 *               `\.      |    ,/                     
 *                  `\.   '. ,/                      
 *                     `\. |/                                
 *                        `3                                   
 *                           `\.
 *                              ` w       
 * @author Riccardo Rossi
 * @author Janosch Stascheit
 * @author Felix Nagel
 */
template<class TPointType> class Tetrahedra3D4 : public Geometry<TPointType>
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
    typedef Triangle3D3<TPointType> FaceType;

    /**
     * Pointer definition of Tetrahedra3D4
     */
    KRATOS_CLASS_POINTER_DEFINITION(Tetrahedra3D4);

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

//     Tetrahedra3D4( const PointType& Point1, const PointType& Point2,
//                    const PointType& Point3, const PointType& Point4 )
//         : BaseType(PointsArrayType(), &msGeometryData)
//     {
//         this->Points().reserve(4);
//         this->Points().push_back(typename PointType::Pointer(new PointType(Point1)));
//         this->Points().push_back(typename PointType::Pointer(new PointType(Point2)));
//         this->Points().push_back(typename PointType::Pointer(new PointType(Point3)));
//         this->Points().push_back(typename PointType::Pointer(new PointType(Point4)));
//     }

    Tetrahedra3D4( typename PointType::Pointer pPoint1,
                   typename PointType::Pointer pPoint2,
                   typename PointType::Pointer pPoint3,
                   typename PointType::Pointer pPoint4 )
        : BaseType(PointsArrayType(), &msGeometryData)
    {
        this->Points().reserve(4);
        this->Points().push_back(pPoint1);
        this->Points().push_back(pPoint2);
        this->Points().push_back(pPoint3);
        this->Points().push_back(pPoint4);
    }

    Tetrahedra3D4( const PointsArrayType& ThisPoints)
        : BaseType(ThisPoints, &msGeometryData)
    {
        if( this->PointsNumber() != 4)
            KRATOS_ERROR << "Invalid points number. Expected 4, given " << this->PointsNumber() << std::endl;
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
    Tetrahedra3D4(Tetrahedra3D4 const& rOther)
        : BaseType(rOther)
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
    template<class TOtherPointType> Tetrahedra3D4(Tetrahedra3D4<TOtherPointType> const& rOther)
        : BaseType(rOther)
    {
    }

    /// Destructor. Does nothing!!!
    ~Tetrahedra3D4() override {}

    GeometryData::KratosGeometryFamily GetGeometryFamily() const override
    {
        return GeometryData::Kratos_Tetrahedra;
    }
    GeometryData::KratosGeometryType GetGeometryType() const override
    {
        return GeometryData::Kratos_Tetrahedra3D4;
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
    Tetrahedra3D4& operator=(const Tetrahedra3D4& rOther)
    {
        BaseType::operator=(rOther);
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
    Tetrahedra3D4& operator=(Tetrahedra3D4<TOtherPointType> const & rOther)
    {
        BaseType::operator=(rOther);

        return *this;
    }


    /**
     * Operations
     */

    typename BaseType::Pointer Create(PointsArrayType const& ThisPoints) const override
    {
        return typename BaseType::Pointer(new Tetrahedra3D4(ThisPoints));
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
    //     Geometry< Point<3> >::Pointer p_clone( new Tetrahedra3D4< Point<3> >( NewPoints ) );

    //     return p_clone;
    // }

    //lumping factors for the calculation of the lumped mass matrix
    Vector& LumpingFactors(Vector& rResult) const override
    {
        if(rResult.size() != 4)
            rResult.resize(4, false);
        std::fill(rResult.begin(), rResult.end(), 1.00 / 4.00);
        return rResult;
    }

    /**
     * Informations
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
        constexpr double factor = 2.0396489026555;                              // (12/sqrt(2)) ^ 1/3);
        return  factor * std::pow(std::fabs(Volume()), 0.33333333333333);            // sqrt(fabs( DeterminantOfJacobian(PointType())));
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
        return this->Volume();
    }

    /** This method calculates and returns the volume of this geometry.
     * This method calculates and returns the volume of this geometry.
     *
     * This method uses the V = (A x B) * C / 6 formula.
     *
     * @return double value contains length, area or volume.
     *
     * @see Length()
     * @see Area()
     * @see Volume()
     *
     * :TODO: might be necessary to reimplement
     */
    double Volume() const override {
        //closed formula for the linear triangle
        const double onesixth = 1.0/6.0;

        const CoordinatesArrayType& rP0 = this->Points()[0].Coordinates();
        const CoordinatesArrayType& rP1 = this->Points()[1].Coordinates();
        const CoordinatesArrayType& rP2 = this->Points()[2].Coordinates();
        const CoordinatesArrayType& rP3 = this->Points()[3].Coordinates();

        double x10 = rP1[0] - rP0[0];
        double y10 = rP1[1] - rP0[1];
        double z10 = rP1[2] - rP0[2];

        double x20 = rP2[0] - rP0[0];
        double y20 = rP2[1] - rP0[1];
        double z20 = rP2[2] - rP0[2];

        double x30 = rP3[0] - rP0[0];
        double y30 = rP3[1] - rP0[1];
        double z30 = rP3[2] - rP0[2];

        double detJ = x10 * y20 * z30 - x10 * y30 * z20 + y10 * z20 * x30 - y10 * x20 * z30 + z10 * x20 * y30 - z10 * y20 * x30;
        return  detJ*onesixth;
    }

    double DomainSize() const override {
      return Volume();
    }

    /** This method calculates and returns the minimum edge length of the geometry.
     *
     * @return double value with the minimum edge length
     *
     * @see MaxEdgeLength()
     * @see AverageEdgeLength()
     */
    double MinEdgeLength() const override {
      auto a = this->GetPoint(0) - this->GetPoint(1);
      auto b = this->GetPoint(1) - this->GetPoint(2);
      auto c = this->GetPoint(2) - this->GetPoint(0);
      auto d = this->GetPoint(3) - this->GetPoint(0);
      auto e = this->GetPoint(3) - this->GetPoint(1);
      auto f = this->GetPoint(3) - this->GetPoint(2);

      double sa = (a[0]*a[0])+(a[1]*a[1])+(a[2]*a[2]);
      double sb = (b[0]*b[0])+(b[1]*b[1])+(b[2]*b[2]);
      double sc = (c[0]*c[0])+(c[1]*c[1])+(c[2]*c[2]);
      double sd = (d[0]*d[0])+(d[1]*d[1])+(d[2]*d[2]);
      double se = (e[0]*e[0])+(e[1]*e[1])+(e[2]*e[2]);
      double sf = (f[0]*f[0])+(f[1]*f[1])+(f[2]*f[2]);

      return CalculateMinEdgeLength(sa, sb, sc, sd, se, sf);
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
      auto a = this->GetPoint(0) - this->GetPoint(1);
      auto b = this->GetPoint(1) - this->GetPoint(2);
      auto c = this->GetPoint(2) - this->GetPoint(0);
      auto d = this->GetPoint(3) - this->GetPoint(0);
      auto e = this->GetPoint(3) - this->GetPoint(1);
      auto f = this->GetPoint(3) - this->GetPoint(2);

      double sa = (a[0]*a[0])+(a[1]*a[1])+(a[2]*a[2]);
      double sb = (b[0]*b[0])+(b[1]*b[1])+(b[2]*b[2]);
      double sc = (c[0]*c[0])+(c[1]*c[1])+(c[2]*c[2]);
      double sd = (d[0]*d[0])+(d[1]*d[1])+(d[2]*d[2]);
      double se = (e[0]*e[0])+(e[1]*e[1])+(e[2]*e[2]);
      double sf = (f[0]*f[0])+(f[1]*f[1])+(f[2]*f[2]);

      return CalculateMaxEdgeLength(sa, sb, sc, sd, se, sf);
    }

    /** This method calculates and returns the average edge length of the geometry
     *
     * @return double value with the average edge length
     *
     * @see MinEdgeLength()
     * @see MaxEdgeLength()
     */
    double AverageEdgeLength() const override {
      return CalculateAvgEdgeLength(
        MathUtils<double>::Norm3(this->GetPoint(0)-this->GetPoint(1)),
        MathUtils<double>::Norm3(this->GetPoint(1)-this->GetPoint(2)),
        MathUtils<double>::Norm3(this->GetPoint(2)-this->GetPoint(0)),
        MathUtils<double>::Norm3(this->GetPoint(3)-this->GetPoint(0)),
        MathUtils<double>::Norm3(this->GetPoint(3)-this->GetPoint(1)),
        MathUtils<double>::Norm3(this->GetPoint(3)-this->GetPoint(2))
      );
    }

    /** This method calculates the circumradius of the geometry
     *
     * @return The circumradius of the geometry
     */
    double Circumradius() const override {
      //    | s10 y10 z10 |           | s10 x10 z10 |           | s10 x10 y10 |           | x10 y10 z10 |
      // MX | s20 y20 z20 |        MY | s20 x20 z20 |        MZ | s20 x20 y20 |        AL | x20 y20 z20 |
      //    | s30 y30 z30 |           | s30 x30 z30 |           | s30 x30 y30 |           | x30 y30 z30 |

      const CoordinatesArrayType& rP0 = this->Points()[0].Coordinates();
      const CoordinatesArrayType& rP1 = this->Points()[1].Coordinates();
      const CoordinatesArrayType& rP2 = this->Points()[2].Coordinates();
      const CoordinatesArrayType& rP3 = this->Points()[3].Coordinates();

      double aDot = rP0[0] * rP0[0] + rP0[1] * rP0[1] + rP0[2] * rP0[2];
      double bDot = rP1[0] * rP1[0] + rP1[1] * rP1[1] + rP1[2] * rP1[2];
      double cDot = rP2[0] * rP2[0] + rP2[1] * rP2[1] + rP2[2] * rP2[2];
      double dDot = rP3[0] * rP3[0] + rP3[1] * rP3[1] + rP3[2] * rP3[2];

      // Build the simplified matrices
      double s10 = aDot - dDot;
      double x10 = rP0[0] - rP3[0];
      double y10 = rP0[1] - rP3[1];
      double z10 = rP0[2] - rP3[2];

      double s20 = bDot - dDot;
      double x20 = rP1[0] - rP3[0];
      double y20 = rP1[1] - rP3[1];
      double z20 = rP1[2] - rP3[2];

      double s30 = cDot - dDot;
      double x30 = rP2[0] - rP3[0];
      double y30 = rP2[1] - rP3[1];
      double z30 = rP2[2] - rP3[2];

      double detJMX = s10 * y20 * z30 + y10 * z20 * s30 + z10 * s20 * y30 - s30 * y20 * z10 - y30 * z20 * s10 - z30 * s20 * y10;
      double detJMY = s10 * x20 * z30 + x10 * z20 * s30 + z10 * s20 * x30 - s30 * x20 * z10 - x30 * z20 * s10 - z30 * s20 * x10;
      double detJMZ = s10 * x20 * y30 + x10 * y20 * s30 + y10 * s20 * x30 - s30 * x20 * y10 - x30 * y20 * s10 - y30 * s20 * x10;
      double detJAL = x10 * y20 * z30 + y10 * z20 * x30 + z10 * x20 * y30 - x30 * y20 * z10 - y30 * z20 * x10 - z30 * x20 * y10;

      return std::sqrt( (detJMX*detJMX + detJMY*detJMY + detJMZ*detJMZ)) / (2*std::abs(detJAL));
    }

    /** This method calculates the inradius of the
     * geometry
     *
     * @return The inradius of the geometry
     */
    double Inradius() const override {
      //    | x10 y10 z10 |
      // AL | x20 y20 z20 |
      //    | x30 y30 z30 |

      const CoordinatesArrayType& rP0 = this->Points()[0].Coordinates();
      const CoordinatesArrayType& rP1 = this->Points()[1].Coordinates();
      const CoordinatesArrayType& rP2 = this->Points()[2].Coordinates();
      const CoordinatesArrayType& rP3 = this->Points()[3].Coordinates();

      array_1d<double, 3> c012, c013, c023, c123;

      MathUtils<double>::CrossProduct(c012, rP2-rP0, rP1-rP0);
      MathUtils<double>::CrossProduct(c013, rP3-rP0, rP1-rP0);
      MathUtils<double>::CrossProduct(c023, rP3-rP0, rP2-rP0);
      MathUtils<double>::CrossProduct(c123, rP3-rP1, rP2-rP1);

      double n012 = std::sqrt(c012[0]*c012[0] + c012[1]*c012[1] + c012[2]*c012[2]);
      double n013 = std::sqrt(c013[0]*c013[0] + c013[1]*c013[1] + c013[2]*c013[2]);
      double n023 = std::sqrt(c023[0]*c023[0] + c023[1]*c023[1] + c023[2]*c023[2]);
      double n123 = std::sqrt(c123[0]*c123[0] + c123[1]*c123[1] + c123[2]*c123[2]);

      // Build the simplified matrices
      double x10 = rP0[0] - rP3[0];
      double y10 = rP0[1] - rP3[1];
      double z10 = rP0[2] - rP3[2];

      double x20 = rP1[0] - rP3[0];
      double y20 = rP1[1] - rP3[1];
      double z20 = rP1[2] - rP3[2];

      double x30 = rP2[0] - rP3[0];
      double y30 = rP2[1] - rP3[1];
      double z30 = rP2[2] - rP3[2];

      double detJAL = x10 * y20 * z30 + y10 * z20 * x30 + z10 * x20 * y30 - x30 * y20 * z10 - y30 * z20 * x10 - z30 * x20 * y10;

      return std::abs(detJAL) / (n012 + n013 + n023 + n123);
    }

    /// Quality functions

    /** Calculates the inradius to circumradius quality metric.
     * Calculates the inradius to circumradius quality metric.
     * This metric is bounded by the interval (0,1) being:
     *  1 -> Optimal value
     *  0 -> Worst value
     *
     * \f$ \frac{r}{\rho} \f$
     *
     * @return The inradius to circumradius quality metric.
     */
    double InradiusToCircumradiusQuality() const override {
      constexpr double normFactor = 3.0;

      return normFactor * Inradius() / Circumradius();
    };

    /** Calculates the inradius to longest edge quality metric.
     * Calculates the inradius to longest edge quality metric.
     * This metric is bounded by the interval (0,1) being:
     *  1 -> Optimal value
     *  0 -> Worst value
     *
     * \f$ \frac{r}{L} \f$
     *
     * @return The inradius to longest edge quality metric.
     */
    double InradiusToLongestEdgeQuality() const override {
      constexpr double normFactor = 4.89897982161;

      auto a = this->GetPoint(0) - this->GetPoint(1);
      auto b = this->GetPoint(1) - this->GetPoint(2);
      auto c = this->GetPoint(2) - this->GetPoint(0);
      auto d = this->GetPoint(3) - this->GetPoint(0);
      auto e = this->GetPoint(3) - this->GetPoint(1);
      auto f = this->GetPoint(3) - this->GetPoint(2);

      double sa = (a[0]*a[0])+(a[1]*a[1])+(a[2]*a[2]);
      double sb = (b[0]*b[0])+(b[1]*b[1])+(b[2]*b[2]);
      double sc = (c[0]*c[0])+(c[1]*c[1])+(c[2]*c[2]);
      double sd = (d[0]*d[0])+(d[1]*d[1])+(d[2]*d[2]);
      double se = (e[0]*e[0])+(e[1]*e[1])+(e[2]*e[2]);
      double sf = (f[0]*f[0])+(f[1]*f[1])+(f[2]*f[2]);

      return normFactor * Inradius() / CalculateMaxEdgeLength(sa,sb,sc,sd,se,sf);
    }

    /** Calculates the shortest to longest edge quality metric.
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
      auto a = this->GetPoint(0) - this->GetPoint(1);
      auto b = this->GetPoint(1) - this->GetPoint(2);
      auto c = this->GetPoint(2) - this->GetPoint(0);
      auto d = this->GetPoint(3) - this->GetPoint(0);
      auto e = this->GetPoint(3) - this->GetPoint(1);
      auto f = this->GetPoint(3) - this->GetPoint(2);

      double sa = (a[0]*a[0])+(a[1]*a[1])+(a[2]*a[2]);
      double sb = (b[0]*b[0])+(b[1]*b[1])+(b[2]*b[2]);
      double sc = (c[0]*c[0])+(c[1]*c[1])+(c[2]*c[2]);
      double sd = (d[0]*d[0])+(d[1]*d[1])+(d[2]*d[2]);
      double se = (e[0]*e[0])+(e[1]*e[1])+(e[2]*e[2]);
      double sf = (f[0]*f[0])+(f[1]*f[1])+(f[2]*f[2]);

      return CalculateMinEdgeLength(sa,sb,sc,sd,se,sf) / CalculateMaxEdgeLength(sa,sb,sc,sd,se,sf);
    }

    /** Calculates the Regularity quality metric.
     * Calculates the Regularity quality metric.
     *  1 -> Optimal value
     *  0 -> Worst value
     *
     * \f$ \frac{4r}{H} \f$
     *
     * @return regularity quality.
     */
    double RegularityQuality() const override {
      KRATOS_ERROR << "Method 'RegularityQuality' is not yet implemented for Tetrahedra3D4" << std::endl;
      return 0.0;
    }

    /** Calculates the volume to surface area quality metric.
     * Calculates the volume to surface area quality metric.
     *  1 -> Optimal value
     *  0 -> Worst value
     *
     * \f$ \frac{V^4}{(\sum{A_{i}^{2}})^{3}} \f$
     *
     * @return volume to surface quality.
     */
    double VolumeToSurfaceAreaQuality() const override {
      KRATOS_ERROR << "Method 'VolumeToSurfaceAreaQuality' is not yet implemented for Tetrahedra3D4" << std::endl;
      return 0.0;
    }

    /** Calculates the Volume to edge length quaility metric.
     * Calculates the Volume to edge length quaility metric.
     *  1 -> Optimal value
     *  0 -> Worst value
     *
     * \f$ \frac{V^{2/3}}{\sum{l_{i}^{2}}} \f$
     *
     * @return Volume to edge length quality.
     */
    double VolumeToEdgeLengthQuality() const override {
      constexpr double normFactor = 12.0;

      auto a = this->GetPoint(0) - this->GetPoint(1);
      auto b = this->GetPoint(1) - this->GetPoint(2);
      auto c = this->GetPoint(2) - this->GetPoint(0);
      auto d = this->GetPoint(3) - this->GetPoint(0);
      auto e = this->GetPoint(3) - this->GetPoint(1);
      auto f = this->GetPoint(3) - this->GetPoint(2);

      double sa = (a[0]*a[0]) + (a[1]*a[1]) + (a[2]*a[2]);
      double sb = (b[0]*b[0]) + (b[1]*b[1]) + (b[2]*b[2]);
      double sc = (c[0]*c[0]) + (c[1]*c[1]) + (c[2]*c[2]);
      double sd = (d[0]*d[0]) + (d[1]*d[1]) + (d[2]*d[2]);
      double se = (e[0]*e[0]) + (e[1]*e[1]) + (e[2]*e[2]);
      double sf = (f[0]*f[0]) + (f[1]*f[1]) + (f[2]*f[2]);

      double vol = Volume();

      return std::abs(normFactor * std::pow(9 * vol * vol, 1.0 / 3.0) / (sa + sb + sc + sd + se + sf)) * (vol < 0 ? -1 : 1);
    }

    /** Calculates the volume to average edge lenght quality metric.
     * Calculates the volume to average edge lenght quality metric.
     *  1 -> Optimal value
     *  0 -> Worst value
     *
     * \f$ \frac{V}{\frac{1}{6}\sum{l_i}} \f$
     *
     * @return [description]
     */
    double VolumeToAverageEdgeLength() const override {
      constexpr double normFactor = 6.0 * 1.41421356237309504880;

      return normFactor * Volume() / std::pow(AverageEdgeLength(), 3);
    }

    /** Calculates the volume to average edge length quality metric.
     * Calculates the volume to average edge length quality metric.
     * The average edge lenght is calculated using the RMS.
     * This metric is bounded by the interval (0,1) being:
     *  1 -> Optimal value
     *  0 -> Worst value
     *
     * \f$ \frac{V}{\sqrt{\frac{1}{6}\sum{A_{i}^{2}}}} \f$
     *
     * @return [description]
     */
    double VolumeToRMSEdgeLength() const override {
      constexpr double normFactor = 6.0 * 1.41421356237309504880;

      auto a = this->GetPoint(0) - this->GetPoint(1);
      auto b = this->GetPoint(1) - this->GetPoint(2);
      auto c = this->GetPoint(2) - this->GetPoint(0);
      auto d = this->GetPoint(3) - this->GetPoint(0);
      auto e = this->GetPoint(3) - this->GetPoint(1);
      auto f = this->GetPoint(3) - this->GetPoint(2);

      double sa = (a[0]*a[0])+(a[1]*a[1])+(a[2]*a[2]);
      double sb = (b[0]*b[0])+(b[1]*b[1])+(b[2]*b[2]);
      double sc = (c[0]*c[0])+(c[1]*c[1])+(c[2]*c[2]);
      double sd = (d[0]*d[0])+(d[1]*d[1])+(d[2]*d[2]);
      double se = (e[0]*e[0])+(e[1]*e[1])+(e[2]*e[2]);
      double sf = (f[0]*f[0])+(f[1]*f[1])+(f[2]*f[2]);

      return normFactor * Volume() / std::pow(std::sqrt(1.0/6.0 * (sa + sb + sc + sd + se + sf)), 3.0);
    }

    /**
    * Returns a matrix of the local coordinates of all points
    * @param rResult a Matrix that will be overwritten by the results
    * @return the coordinates of all points of the current geometry
    */
    Matrix& PointsLocalCoordinates( Matrix& rResult ) const override
    {
        if(rResult.size1()!= 4 || rResult.size2()!= 3)
            rResult.resize(4, 3, false);

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

        return rResult;
    }

    /**
     * @brief Returns the local coordinates of a given arbitrary point
     * @details Based on https://www.colorado.edu/engineering/CAS/courses.d/AFEM.d/AFEM.Ch09.d/AFEM.Ch09.pdf. Section 9.1.6
     * @param rResult The vector containing the local coordinates of the point
     * @param rPoint The point in global coordinates
     * @return The vector containing the local coordinates of the point
     */
    CoordinatesArrayType& PointLocalCoordinates(
            CoordinatesArrayType& rResult,
            const CoordinatesArrayType& rPoint
            ) override
    {
        // Compute RHS
        array_1d<double,4> X;
        X[0] = 1.0;
        X[1] = rPoint[0];
        X[2] = rPoint[1];
        X[3] = rPoint[2];

        // Auxiliar coordinates
        const double x1 = this->GetPoint( 0 ).X();
        const double x2 = this->GetPoint( 1 ).X();
        const double x3 = this->GetPoint( 2 ).X();
        const double x4 = this->GetPoint( 3 ).X();
        const double y1 = this->GetPoint( 0 ).Y();
        const double y2 = this->GetPoint( 1 ).Y();
        const double y3 = this->GetPoint( 2 ).Y();
        const double y4 = this->GetPoint( 3 ).Y();
        const double z1 = this->GetPoint( 0 ).Z();
        const double z2 = this->GetPoint( 1 ).Z();
        const double z3 = this->GetPoint( 2 ).Z();
        const double z4 = this->GetPoint( 3 ).Z();

        // Auxiliar diff
        const double x12 = x1 - x2;
        const double x13 = x1 - x3;
        const double x14 = x1 - x4;
        const double x21 = x2 - x1;
        const double x24 = x2 - x4;
        const double x31 = x3 - x1;
        const double x32 = x3 - x2;
        const double x34 = x3 - x4;
        const double x42 = x4 - x2;
        const double x43 = x4 - x3;
        const double y12 = y1 - y2;
        const double y13 = y1 - y3;
        const double y14 = y1 - y4;
        const double y21 = y2 - y1;
        const double y24 = y2 - y4;
        const double y31 = y3 - y1;
        const double y32 = y3 - y2;
        const double y34 = y3 - y4;
        const double y42 = y4 - y2;
        const double y43 = y4 - y3;
        const double z12 = z1 - z2;
        const double z13 = z1 - z3;
        const double z14 = z1 - z4;
        const double z21 = z2 - z1;
        const double z24 = z2 - z4;
        const double z31 = z3 - z1;
        const double z32 = z3 - z2;
        const double z34 = z3 - z4;
        const double z42 = z4 - z2;
        const double z43 = z4 - z3;

        // Compute LHS
        BoundedMatrix<double, 4,4> invJ;
        const double aux_volume = 1.0/(6.0*this->Volume());
        invJ(0,0) = aux_volume * (x2*(y3*z4-y4*z3)+x3*(y4*z2-y2*z4)+x4*(y2*z3-y3*z2));
        invJ(1,0) = aux_volume * (x1*(y4*z3-y3*z4)+x3*(y1*z4-y4*z1)+x4*(y3*z1-y1*z3));
        invJ(2,0) = aux_volume * (x1*(y2*z4-y4*z2)+x2*(y4*z1-y1*z4)+x4*(y1*z2-y2*z1));
        invJ(3,0) = aux_volume * (x1*(y3*z2-y2*z3)+x2*(y1*z3-y3*z1)+x3*(y2*z1-y1*z2));
        invJ(0,1) = aux_volume * (y42*z32 - y32*z42);
        invJ(1,1) = aux_volume * (y31*z43 - y34*z13);
        invJ(2,1) = aux_volume * (y24*z14 - y14*z24);
        invJ(3,1) = aux_volume * (y13*z21 - y12*z31);
        invJ(0,2) = aux_volume * (x32*z42 - x42*z32);
        invJ(1,2) = aux_volume * (x43*z31 - x13*z34);
        invJ(2,2) = aux_volume * (x14*z24 - x24*z14);
        invJ(3,2) = aux_volume * (x21*z13 - x31*z12);
        invJ(0,3) = aux_volume * (x42*y32 - x32*y42);
        invJ(1,3) = aux_volume * (x31*y43 - x34*y13);
        invJ(2,3) = aux_volume * (x24*y14 - x14*y24);
        invJ(3,3) = aux_volume * (x13*y21 - x12*y31);

        const array_1d<double,4> result = prod(invJ, X);

        if (rResult.size() != 3)
            rResult.resize(3, false);

        rResult[0] = result[1];
        rResult[1] = result[2];
        rResult[2] = result[3];

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
        ) override
    {
        this->PointLocalCoordinates( rResult, rPoint );

        if( rResult[0] >= 0.0-Tolerance )
        {
            if( rResult[1] >= 0.0-Tolerance )
            {
                if( rResult[2] >= 0.0-Tolerance )
                {
                    if( (rResult[0] + rResult[1] + rResult[2]) <= (1.0+Tolerance))
                    {
                        return true;
                    }
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
        return 6;
    }

    SizeType FacesNumber() const override
    {
        return 4;
    }

    /** This method gives you all edges of this geometry.

    @return GeometriesArrayType containes this geometry edges.
    @see EdgesNumber()
    @see Edge()
     */
    GeometriesArrayType Edges(void) override
    {
        GeometriesArrayType edges = GeometriesArrayType();
        typedef typename Geometry<TPointType>::Pointer EdgePointerType;
        edges.push_back( EdgePointerType(new EdgeType(
                                             this->pGetPoint(0),
                                             this->pGetPoint(1))) );
        edges.push_back( EdgePointerType(new EdgeType(
                                             this->pGetPoint(1),
                                             this->pGetPoint(2))) );
        edges.push_back( EdgePointerType(new EdgeType(
                                             this->pGetPoint(2),
                                             this->pGetPoint(0))) );
        edges.push_back( EdgePointerType(new EdgeType(
                                             this->pGetPoint(0),
                                             this->pGetPoint(3))) );
        edges.push_back( EdgePointerType(new EdgeType(
                                             this->pGetPoint(1),
                                             this->pGetPoint(3))) );
        edges.push_back( EdgePointerType(new EdgeType(
                                             this->pGetPoint(2),
                                             this->pGetPoint(3))) );

        return edges;
    }

    GeometriesArrayType Faces(void) override
    {
        GeometriesArrayType faces = GeometriesArrayType();
        typedef typename Geometry<TPointType>::Pointer FacePointerType;
        faces.push_back( FacePointerType(new FaceType(
                                             this->pGetPoint(0),
                                             this->pGetPoint(2),
                                             this->pGetPoint(1))) );
        faces.push_back( FacePointerType(new FaceType(
                                             this->pGetPoint(0),
                                             this->pGetPoint(3),
                                             this->pGetPoint(2))) );
        faces.push_back( FacePointerType(new FaceType(
                                             this->pGetPoint(0),
                                             this->pGetPoint(1),
                                             this->pGetPoint(3))) );
        faces.push_back( FacePointerType(new FaceType(
                                             this->pGetPoint(2),
                                             this->pGetPoint(3),
                                             this->pGetPoint(1))) );
        return faces;
    }

    //Connectivities of faces required
    void NumberNodesInFaces (DenseVector<unsigned int>& NumberNodesInFaces) const override
    {
        NumberNodesInFaces.resize(4, false);
        // Linear Tetrahedra have elements of 3 nodes as faces
        NumberNodesInFaces[0]=3;
        NumberNodesInFaces[1]=3;
        NumberNodesInFaces[2]=3;
        NumberNodesInFaces[3]=3;

    }


    void NodesInFaces (DenseMatrix<unsigned int>& NodesInFaces) const override
    {
        // faces in columns
      if(NodesInFaces.size1() != 4 || NodesInFaces.size2() != 4)
        NodesInFaces.resize(4, 4, false);

      //face 1
      NodesInFaces(0,0)=0;//contrary node to the face
      NodesInFaces(1,0)=1;
      NodesInFaces(2,0)=2;
      NodesInFaces(3,0)=3;
      //face 2
      NodesInFaces(0,1)=1;//contrary node to the face
      NodesInFaces(1,1)=2;
      NodesInFaces(2,1)=0;
      NodesInFaces(3,1)=3;
      //face 3
      NodesInFaces(0,2)=2;//contrary node to the face
      NodesInFaces(1,2)=0;
      NodesInFaces(2,2)=1;
      NodesInFaces(3,2)=3;
      //face 4
      NodesInFaces(0,3)=3;//contrary node to the face
      NodesInFaces(1,3)=0;
      NodesInFaces(2,3)=2;
      NodesInFaces(3,3)=1;

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
                                       const CoordinatesArrayType& rPoint) const override
    {
        switch( ShapeFunctionIndex )
        {
        case 0:
            return( 1.0-(rPoint[0]+rPoint[1]+rPoint[2]) );
        case 1:
            return( rPoint[0] );
        case 2:
            return( rPoint[1] );
        case 3:
            return( rPoint[2] );
        default:
            KRATOS_ERROR << "Wrong index of shape function!" << *this << std::endl;
        }
        return 0;
    }

    /** This method gives all non-zero shape functions values
    evaluated at the rCoordinates provided

    \note There is no control if the return vector is empty or not!

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
      rResult[0] =  1.0-(rCoordinates[0]+rCoordinates[1]+rCoordinates[2]);
      rResult[1] =  rCoordinates[0] ;
      rResult[2] =  rCoordinates[1] ;
      rResult[3] =  rCoordinates[2] ;

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
        if(rResult.size1() != this->PointsNumber() || rResult.size2() != this->LocalSpaceDimension())
            rResult.resize(this->PointsNumber(),this->LocalSpaceDimension(),false);

        CalculateShapeFunctionsLocalGradients(rResult, rPoint);

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
    ShapeFunctionsGradientsType& ShapeFunctionsIntegrationPointsGradients(
        ShapeFunctionsGradientsType& rResult,
        IntegrationMethod ThisMethod) const override
    {
        const unsigned int integration_points_number =
            msGeometryData.IntegrationPointsNumber(ThisMethod);
        if(integration_points_number == 0)
            KRATOS_ERROR << "This integration method is not supported" << *this << std::endl;

        BoundedMatrix<double,4,3> DN_DX;
        const double x10 = this->Points()[1].X() - this->Points()[0].X();
        const double y10 = this->Points()[1].Y() - this->Points()[0].Y();
        const double z10 = this->Points()[1].Z() - this->Points()[0].Z();

        const double x20 = this->Points()[2].X() - this->Points()[0].X();
        const double y20 = this->Points()[2].Y() - this->Points()[0].Y();
        const double z20 = this->Points()[2].Z() - this->Points()[0].Z();

        const double x30 = this->Points()[3].X() - this->Points()[0].X();
        const double y30 = this->Points()[3].Y() - this->Points()[0].Y();
        const double z30 = this->Points()[3].Z() - this->Points()[0].Z();

        const double detJ = x10 * y20 * z30 - x10 * y30 * z20 + y10 * z20 * x30 - y10 * x20 * z30 + z10 * x20 * y30 - z10 * y20 * x30;

        DN_DX(0,0) = -y20 * z30 + y30 * z20 + y10 * z30 - z10 * y30 - y10 * z20 + z10 * y20;
        DN_DX(0,1) = -z20 * x30 + x20 * z30 - x10 * z30 + z10 * x30 + x10 * z20 - z10 * x20;
        DN_DX(0,2) = -x20 * y30 + y20 * x30 + x10 * y30 - y10 * x30 - x10 * y20 + y10 * x20;
        DN_DX(1,0) = y20 * z30 - y30 * z20;
        DN_DX(1,1) = z20 * x30 - x20 * z30;
        DN_DX(1,2) = x20 * y30 - y20 * x30;
        DN_DX(2,0) = -y10 * z30 + z10 * y30;
        DN_DX(2,1) = x10 * z30 - z10 * x30;
        DN_DX(2,2) = -x10 * y30 + y10 * x30;
        DN_DX(3,0) = y10 * z20 - z10 * y20;
        DN_DX(3,1) = -x10 * z20 + z10 * x20;
        DN_DX(3,2) = x10 * y20 - y10 * x20;

        DN_DX /= detJ;
        if(rResult.size() != integration_points_number)
        {
            rResult.resize(integration_points_number,false);
        }

        for(unsigned int i=0; i<integration_points_number; i++)
                rResult[i] = DN_DX;


        return rResult;
    }


    ShapeFunctionsGradientsType& ShapeFunctionsIntegrationPointsGradients(
        ShapeFunctionsGradientsType& rResult
        , Vector& determinants_of_jacobian
        , IntegrationMethod ThisMethod) const override
    {
        const unsigned int integration_points_number =
            msGeometryData.IntegrationPointsNumber(ThisMethod);
        if(integration_points_number == 0)
            KRATOS_ERROR << "This integration method is not supported" << *this << std::endl;

        BoundedMatrix<double,4,3> DN_DX;
        const double x10 = this->Points()[1].X() - this->Points()[0].X();
        const double y10 = this->Points()[1].Y() - this->Points()[0].Y();
        const double z10 = this->Points()[1].Z() - this->Points()[0].Z();

        const double x20 = this->Points()[2].X() - this->Points()[0].X();
        const double y20 = this->Points()[2].Y() - this->Points()[0].Y();
        const double z20 = this->Points()[2].Z() - this->Points()[0].Z();

        const double x30 = this->Points()[3].X() - this->Points()[0].X();
        const double y30 = this->Points()[3].Y() - this->Points()[0].Y();
        const double z30 = this->Points()[3].Z() - this->Points()[0].Z();

        const double detJ = x10 * y20 * z30 - x10 * y30 * z20 + y10 * z20 * x30 - y10 * x20 * z30 + z10 * x20 * y30 - z10 * y20 * x30;

        DN_DX(0,0) = -y20 * z30 + y30 * z20 + y10 * z30 - z10 * y30 - y10 * z20 + z10 * y20;
        DN_DX(0,1) = -z20 * x30 + x20 * z30 - x10 * z30 + z10 * x30 + x10 * z20 - z10 * x20;
        DN_DX(0,2) = -x20 * y30 + y20 * x30 + x10 * y30 - y10 * x30 - x10 * y20 + y10 * x20;
        DN_DX(1,0) = y20 * z30 - y30 * z20;
        DN_DX(1,1) = z20 * x30 - x20 * z30;
        DN_DX(1,2) = x20 * y30 - y20 * x30;
        DN_DX(2,0) = -y10 * z30 + z10 * y30;
        DN_DX(2,1) = x10 * z30 - z10 * x30;
        DN_DX(2,2) = -x10 * y30 + y10 * x30;
        DN_DX(3,0) = y10 * z20 - z10 * y20;
        DN_DX(3,1) = -x10 * z20 + z10 * x20;
        DN_DX(3,2) = x10 * y20 - y10 * x20;

        DN_DX /= detJ;

        if(determinants_of_jacobian.size() != integration_points_number )
            determinants_of_jacobian.resize(integration_points_number,false);

        for(unsigned int i=0; i<integration_points_number; i++)
            determinants_of_jacobian[i] = detJ;

        // Volume = detJ*0.1666666666666666666667;

        // Workaround by riccardo
        if(rResult.size() != integration_points_number)
        {
            rResult.resize(integration_points_number,false);
        }
        for(unsigned int i=0; i<integration_points_number; i++)
                rResult[i] = DN_DX;


        return rResult;
    }


    /// detect if two tetrahedra are intersected
    bool HasIntersection( const BaseType& rThisGeometry) override
    {

        array_1d<Plane, 4>  plane;
        std::vector<BaseType> Intersection;

        //const BaseType& geom_1 = *this;
        const BaseType& geom_2 = rThisGeometry;

        GetPlanes(plane);
        Intersection.push_back(geom_2);
        for (unsigned int i = 0; i < 4; ++i)
        {
            std::vector<BaseType> inside;
            for (unsigned int j = 0; j < Intersection.size(); ++j)
            {
                SplitAndDecompose(Intersection[j], plane[i], inside);
            }
            Intersection = inside;
        }

        return bool (Intersection.size() > 0);
    }


    bool HasIntersection(const Point& rLowPoint, const Point& rHighPoint) override
    {
        using Triangle3D3Type = Triangle3D3<TPointType>;
        // Check if faces have intersection
        if(Triangle3D3Type(this->pGetPoint(0),this->pGetPoint(2), this->pGetPoint(1)).HasIntersection(rLowPoint, rHighPoint))
            return true;
        if(Triangle3D3Type(this->pGetPoint(0),this->pGetPoint(3), this->pGetPoint(2)).HasIntersection(rLowPoint, rHighPoint))
            return true;
        if(Triangle3D3Type(this->pGetPoint(0),this->pGetPoint(1), this->pGetPoint(3)).HasIntersection(rLowPoint, rHighPoint))
            return true;
        if(Triangle3D3Type(this->pGetPoint(2),this->pGetPoint(3), this->pGetPoint(1)).HasIntersection(rLowPoint, rHighPoint))
            return true;

        CoordinatesArrayType local_coordinates;
        // if there are no faces intersecting the box then or the box is inside the tetrahedron or it does not have intersection
        if(IsInside(rLowPoint,local_coordinates))
            return true;

        return false;
    }


    void SplitAndDecompose(
        const BaseType& tetra, Plane& plane,
        std::vector<BaseType>& inside)
    {

        // Determine on which side of the plane the points of the tetrahedron lie.
        int i = 0;
        int positive = 0, negative = 0, zero = 0;
        array_1d<double,4> C;
        array_1d<int,4>  pos, neg, zer;

        noalias(C)   = ZeroVector(4);

        pos[0] = 0;
        neg[0] = 0;
        zer[0] = 0;
        pos[1] = 0;
        neg[1] = 0;
        zer[1] = 0;
        pos[2] = 0;
        neg[2] = 0;
        zer[2] = 0;
        pos[3] = 0;
        neg[3] = 0;
        zer[3] = 0;

        for (i = 0; i < 4; ++i)
        {
            const array_1d<double,3>& p = tetra[i].Coordinates();
            C[i] = plane.DistanceTo(p);
            if (C[i] > 0.00)
                pos[positive++] = i;
            else if (C[i] < 0.00)
                neg[negative++] = i;
            else
                zer[zero++] = i;
        }


        // For a split to occur, one of the c_i must be positive and one must
        // be negative.
        if (negative == 0)
        {
            // Tetrahedron is completely on the positive side of plane, full clip.
            return;
        }

        if (positive == 0)
        {
            // Tetrahedron is completely on the negative side of plane.
            inside.push_back(tetra);
            return;
        }

        double w0, w1, invCDiff;
        array_1d<array_1d<double, 3 >, 4> intp;
        array_1d<array_1d<double, 3 >, 4> V;

        if (positive == 3)
        {
            // +++-
            for (i = 0; i < positive; ++i)
            {
                invCDiff = 1.00/(C[pos[i]] - C[neg[0]]);
                w0        = -C[neg[0]]*invCDiff;
                w1        = +C[pos[i]]*invCDiff;
                V[pos[i]] = w0*tetra[pos[i]].Coordinates() +  w1*tetra[neg[0]].Coordinates();
            }
            inside.push_back(tetra);
        }

        else if (positive == 2)
        {
            if (negative == 2)
            {
                // ++--
                for (i = 0; i < positive; ++i)
                {
                    invCDiff = (1.00)/(C[pos[i]] - C[neg[0]]);
                    w0 = -C[neg[0]]*invCDiff;
                    w1 = +C[pos[i]]*invCDiff;
                    intp[i] = w0*tetra[pos[i]].Coordinates() + w1*tetra[neg[0]].Coordinates();
                }
                for (i = 0; i < negative; ++i)
                {
                    invCDiff = (1.00)/(C[pos[i]] - C[neg[1]]);
                    w0 = -C[neg[1]]*invCDiff;
                    w1 = +C[pos[i]]*invCDiff;
                    intp[i+2] = w0*tetra[pos[i]].Coordinates() + w1*tetra[neg[1]].Coordinates();
                }

                V[pos[0]] = intp[2];
                V[pos[1]] = intp[1];
                inside.push_back(tetra);
            }
            else
            {
                // ++-0
                for (i = 0; i < positive; ++i)
                {
                    invCDiff = (1.00)/(C[pos[i]] - C[neg[0]]);
                    w0 = -C[neg[0]]*invCDiff;
                    w1 = +C[pos[i]]*invCDiff;
                    V[pos[i]] = w0*tetra[pos[i]].Coordinates() + w1*tetra[neg[0]].Coordinates();
                }
                inside.push_back(tetra);
            }
        }


        else if (positive == 1)
        {
            if (negative == 3)
            {
                // +---
                for (i = 0; i < negative; ++i)
                {
                    invCDiff = (1.00)/(C[pos[0]] - C[neg[i]]);
                    w0       = -C[neg[i]]*invCDiff;
                    w1       = +C[pos[0]]*invCDiff;
                    intp[i]  = w0*tetra[pos[0]] + w1*tetra[neg[i]];
                }

                V[pos[0]] = intp[0];
                inside.push_back(tetra);
            }
            else if (negative == 2)
            {
                // +--0
                for (i = 0; i < negative; ++i)
                {
                    invCDiff = (1.00)/(C[pos[0]] - C[neg[i]]);
                    w0 = -C[neg[i]]*invCDiff;
                    w1 = +C[pos[0]]*invCDiff;
                    intp[i] = w0*tetra[pos[0]] + w1*tetra[neg[i]];
                }

                V[pos[0]] = intp[0];
                inside.push_back(tetra);
            }
            else
            {
                // +-00
                invCDiff        = (1.00)/(C[pos[0]] - C[neg[0]]);
                w0              = -C[neg[0]]*invCDiff;
                w1              = +C[pos[0]]*invCDiff;
                V[pos[0]]       = w0*tetra[pos[0]] + w1*tetra[neg[0]];
                inside.push_back(tetra);
            }
        }
    }


    void GetPlanes(array_1d<Plane, 4>& plane)
    {
        const BaseType& geom_1 = *this;
        array_1d<double, 3> edge10 = geom_1[1].Coordinates() - geom_1[0].Coordinates();
        array_1d<double, 3> edge20 = geom_1[2].Coordinates() - geom_1[0].Coordinates();
        array_1d<double, 3> edge30 = geom_1[3].Coordinates() - geom_1[0].Coordinates();
        array_1d<double, 3> edge21 = geom_1[2].Coordinates() - geom_1[1].Coordinates();
        array_1d<double, 3> edge31 = geom_1[3].Coordinates() - geom_1[1].Coordinates();

        MathUtils<double>::UnitCrossProduct(plane[0].mNormal, edge10, edge20);  // <v0,v2,v1>
        MathUtils<double>::UnitCrossProduct(plane[1].mNormal, edge30, edge10);  // <v0,v1,v3>
        MathUtils<double>::UnitCrossProduct(plane[2].mNormal, edge20, edge30);  // <v0,v3,v2>
        MathUtils<double>::UnitCrossProduct(plane[3].mNormal, edge31, edge21);  // <v1,v2,v3>

        double det = inner_prod(edge10, plane[3].mNormal);
        if (det < 0.00)
        {
            // The normals are inner pointing, reverse their directions.
            for (int i = 0; i < 4; ++i)
            {
                plane[i].mNormal = -plane[i].mNormal;
            }
        }

        for (int i = 0; i < 4; ++i)
        {
            plane[i].mConstant = inner_prod(geom_1[i].Coordinates(), plane[i].mNormal);
        }

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
        return "3 dimensional tetrahedra with four nodes in 3D space";
    }

    /**
     * Print information about this object.
     *
     * @param rOStream Stream to print into it.
     * @see PrintData()
     * @see Info()
     */
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "3 dimensional tetrahedra with four nodes in 3D space";
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
    void PrintData(std::ostream& rOStream) const override
    {
        BaseType::PrintData(rOStream);
        std::cout << std::endl;
        rOStream << "    in Tetrahedra3D4 PrintData\t : " << std::endl;
        Matrix jacobian(3,3);
        this->Jacobian(jacobian, PointType());
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

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseType );
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseType );
    }

    // serialization needs the default constructor
    Tetrahedra3D4(): BaseType(PointsArrayType(), &msGeometryData) {}

    /**
     * Private Operations
     */

    /**
     * Calculates the gradients in terms of local coordinateds
     * of all shape functions in a given point.
     *
     * @param rPoint the current point at which the gradients are calculated
     * @return the gradients of all shape functions
     * \f$ \frac{\partial N^i}{\partial \xi_j} \f$
     */

    static Matrix& CalculateShapeFunctionsLocalGradients( Matrix& rResult, const CoordinatesArrayType& rPoint )
    {
        rResult(0,0) = -1.0;
        rResult(0,1) = -1.0;
        rResult(0,2) = -1.0;
        rResult(1,0) =  1.0;
        rResult(1,1) =  0.0;
        rResult(1,2) =  0.0;
        rResult(2,0) =  0.0;
        rResult(2,1) =  1.0;
        rResult(2,2) =  0.0;
        rResult(3,0) =  0.0;
        rResult(3,1) =  0.0;
        rResult(3,2) =  1.0;
        return rResult;
    }

    /**
     * Calculates the values of all shape function in all integration points.
     * Integration points are expected to be given in local coordinates
     *
     * @param ThisMethod the current integration method
     * @return the matrix of values of every shape function in each integration point
     *
     */
    static Matrix CalculateShapeFunctionsIntegrationPointsValues(
        typename BaseType::IntegrationMethod ThisMethod)
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
        Matrix shape_function_values(integration_points_number, points_number );
        //loop over all integration points
        for(int pnt = 0; pnt < integration_points_number; pnt++)
        {
            shape_function_values(pnt,0) = ( 1.0
                                             -integration_points[pnt].X()
                                             -integration_points[pnt].Y()
                                             -integration_points[pnt].Z() );
            shape_function_values(pnt,1) = integration_points[pnt].X();
            shape_function_values(pnt,2) = integration_points[pnt].Y();
            shape_function_values(pnt,3) = integration_points[pnt].Z();
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
        typename BaseType::IntegrationMethod ThisMethod)
    {
        IntegrationPointsContainerType all_integration_points =
            AllIntegrationPoints();
        IntegrationPointsArrayType integration_points =
            all_integration_points[ThisMethod];
        //number of integration points
        const int integration_points_number = integration_points.size();
        ShapeFunctionsGradientsType d_shape_f_values(integration_points_number);
        //initialising container
        //loop over all integration points
        for( int pnt=0; pnt<integration_points_number; pnt++ )
        {
            Matrix result = ZeroMatrix(4,3);
            result(0,0) = -1.0;
            result(0,1) = -1.0;
            result(0,2) = -1.0;
            result(1,0) =  1.0;
            result(1,1) =  0.0;
            result(1,2) =  0.0;
            result(2,0) =  0.0;
            result(2,1) =  1.0;
            result(2,2) =  0.0;
            result(3,0) =  0.0;
            result(3,1) =  0.0;
            result(3,2) =  1.0;
            d_shape_f_values[pnt] = result;
        }
        return d_shape_f_values;
    }

    static const IntegrationPointsContainerType AllIntegrationPoints()
    {
        IntegrationPointsContainerType integration_points =
        {
            {
                Quadrature<TetrahedronGaussLegendreIntegrationPoints1,
                3, IntegrationPoint<3> >::GenerateIntegrationPoints(),
                Quadrature<TetrahedronGaussLegendreIntegrationPoints2,
                3, IntegrationPoint<3> >::GenerateIntegrationPoints(),
                Quadrature<TetrahedronGaussLegendreIntegrationPoints3,
                3, IntegrationPoint<3> >::GenerateIntegrationPoints(),
                Quadrature<TetrahedronGaussLegendreIntegrationPoints4,
                3, IntegrationPoint<3> >::GenerateIntegrationPoints(),
                Quadrature<TetrahedronGaussLegendreIntegrationPoints5,
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
                Tetrahedra3D4<TPointType>::CalculateShapeFunctionsIntegrationPointsValues(
                    GeometryData::GI_GAUSS_1),
                Tetrahedra3D4<TPointType>::CalculateShapeFunctionsIntegrationPointsValues(
                    GeometryData::GI_GAUSS_2),
                Tetrahedra3D4<TPointType>::CalculateShapeFunctionsIntegrationPointsValues(
                    GeometryData::GI_GAUSS_3),
                Tetrahedra3D4<TPointType>::CalculateShapeFunctionsIntegrationPointsValues(
                    GeometryData::GI_GAUSS_4),
                Tetrahedra3D4<TPointType>::CalculateShapeFunctionsIntegrationPointsValues(
                    GeometryData::GI_GAUSS_5)
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
                Tetrahedra3D4<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients(
                    GeometryData::GI_GAUSS_1),
                Tetrahedra3D4<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients(
                    GeometryData::GI_GAUSS_2),
                Tetrahedra3D4<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients(
                    GeometryData::GI_GAUSS_3),
                Tetrahedra3D4<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients(
                    GeometryData::GI_GAUSS_4),
                Tetrahedra3D4<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients(
                    GeometryData::GI_GAUSS_5)
            }
        };
        return shape_functions_local_gradients;
    }

    /** Implements the calculus of the minimum edge length
     * Implements the calculus of the minimum edge length given the length of the geometry edges.
     *
     * @param  a Length of the edge a
     * @param  b Length of the edge b
     * @param  c Length of the edge c
     * @param  d Length of the edge d
     * @param  e Length of the edge e
     * @param  f Length of the edge f
     *
     * @return   The minimum edge length of the geometry with edges a,b,c,d,e,f
     */
    inline double CalculateMinEdgeLength(double sa, double sb, double sc, double sd, double se, double sf) const {
      return std::sqrt(std::min({sa, sb, sc, sd, se, sf}));
    }

    /** Implements the calculus of the maximum edge length
     * Implements the calculus of the maximum edge length given the length of the geometry edges.
     *
     * @param  a Length of the edge a
     * @param  b Length of the edge b
     * @param  c Length of the edge c
     * @param  d Length of the edge d
     * @param  e Length of the edge e
     * @param  f Length of the edge f
     *
     * @return   The maximum edge length of the geometry with edges a,b,c,d,e,f
     */
    inline double CalculateMaxEdgeLength(double sa, double sb, double sc, double sd, double se, double sf) const {
      return std::sqrt(std::max({sa, sb, sc, sd, se, sf}));
    }

    /** Implements the calculus of the average edge length
     * Implements the calculus of the average edge length given the length of the geometry edges.
     *
     * @param  a Length of the edge a
     * @param  b Length of the edge b
     * @param  c Length of the edge c
     * @param  d Length of the edge d
     * @param  e Length of the edge e
     * @param  f Length of the edge f
     *
     * @return   The average edge length of the geometry with edges a,b,c,d,e,f
     */
    inline double CalculateAvgEdgeLength(double a, double b, double c, double d, double e, double f) const {
      return (a + b + c + d + e + f) * 1.0/6.0;
    }

    /** Implements the calculus of the circumradius
     * Implements the calculus of the circumradius given the length of the geometry edges.
     *
     * @param  a Length of the edge a
     * @param  b Length of the edge b
     * @param  c Length of the edge c
     * @param  d Length of the edge d
     * @param  e Length of the edge e
     * @param  f Length of the edge f
     *
     * @return   The circumradius of the geometry with edges a,b,c
     */
    inline double CalculateCircumradius(double a, double b, double c, double d, double e, double f) const {
      KRATOS_ERROR << "Circumradius function hasn't been implemented yet" << std::endl;
      return 0.0;
    }

    /** Implements the calculus of the inradius
     * Implements the calculus of the inradius given the length of the geometry edges.
     *
     * @param  a Length of the edge a
     * @param  b Length of the edge b
     * @param  c Length of the edge c
     * @param  d Length of the edge d
     * @param  e Length of the edge e
     * @param  f Length of the edge f
     *
     * @return   The inradius of the geometry with edges a,b,c
     */
    inline double CalculateInradius(double a, double b, double c, double d, double e, double f) const {
      KRATOS_ERROR << "Inradius function hasn't been implemented yet" << std::endl;
      return 0.0;
    }


    /**
     * Private Friends
     */

    template<class TOtherPointType> friend class Tetrahedra3D4;


    /**
     * Un accessible methods
     */


};// Class Tetrahedra3D4


/**
 * Input and output
 */

/**
 * input stream function
 */
template<class TPointType> inline std::istream& operator >> (
    std::istream& rIStream, Tetrahedra3D4<TPointType>& rThis);

/**
 * output stream function
 */
template<class TPointType> inline std::ostream& operator << (
    std::ostream& rOStream, const Tetrahedra3D4<TPointType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}


template<class TPointType> const
GeometryData Tetrahedra3D4<TPointType>::msGeometryData(
    3, 3, 3, GeometryData::GI_GAUSS_1,
    Tetrahedra3D4<TPointType>::AllIntegrationPoints(),
    Tetrahedra3D4<TPointType>::AllShapeFunctionsValues(),
    AllShapeFunctionsLocalGradients()
);

}// namespace Kratos.

#endif // KRATOS_TETRAHEDRA_3D_4_H_INCLUDED  defined
