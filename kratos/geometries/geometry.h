//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//                   Riccardo Rossi
//                   Janosch Stascheit
//                   Felix Nagel
//  contributors:    Hoang Giang Bui
//                   Josep Maria Carbonell
//                   Carlos Roig
//

#if !defined(KRATOS_GEOMETRY_H_INCLUDED )
#define  KRATOS_GEOMETRY_H_INCLUDED

// System includes

// External includes

// Project includes
#include "geometries/geometry_data.h"
#include "geometries/point.h"
#include "containers/pointer_vector.h"

#include "utilities/math_utils.h"
#include "input_output/logger.h"

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

///Geometry base class.
/** As a base class Geometry has all the common
 * interface of Kratos' geometries. Also it contains array of
 * pointers to its points, reference to shape functions values in
 * all integrations points and also local gradients of shape
 * functions evaluated in all integrations points.
 *
 * Geometry is a template class with just one template parameter:
 * - TPointType which reperesent the type of the point this geometry
 * type contain and build on.
 *
 * @see Point
 * @see Node
 * @see Formulation
 * @see GeometryAndFormulationElement
 */
template<class TPointType>
class Geometry : public PointerVector<TPointType, 
                                      typename TPointType::Pointer, 
                                      std::vector<typename TPointType::Pointer> 
                                      >
{
public:
    ///@}
    ///@name Type Definitions
    ///@{

    /// This Geometry type.
    typedef Geometry<TPointType> GeometryType;

    /// Pointer definition of Geometry
    KRATOS_CLASS_POINTER_DEFINITION( Geometry );

    /** Different criteria to evaluate the quality of a geometry.
     * Different criteria to evaluate the quality of a geometry.
     */
    enum class QualityCriteria {
      INRADIUS_TO_CIRCUMRADIUS,
      AREA_TO_LENGTH,
      SHORTEST_ALTITUDE_TO_LENGTH,
      INRADIUS_TO_LONGEST_EDGE,
      SHORTEST_TO_LONGEST_EDGE,
      REGULARITY,
      VOLUME_TO_SURFACE_AREA,
      VOLUME_TO_EDGE_LENGTH,
      VOLUME_TO_AVERAGE_EDGE_LENGTH,
      VOLUME_TO_RMS_EDGE_LENGTH
    };


    /** Base type for geometry.
    */
    typedef PointerVector<TPointType, typename TPointType::Pointer, std::vector<typename TPointType::Pointer> > BaseType;


    /** The bounding box */
    /*typedef BoundingBox<TPointType, GeometryType>  BoundingBoxType; */

    /** Array of counted pointers to point. This type used to hold
    geometry's points.
    */
    typedef PointerVector<TPointType, typename TPointType::Pointer, std::vector<typename TPointType::Pointer> > PointsArrayType;

    /** Integration methods implemented in geometry.
    */
    typedef GeometryData::IntegrationMethod IntegrationMethod;

    /** A Vector of counted pointers to Geometries. Used for
    returning edges of the geometry.
     */
    typedef PointerVector<GeometryType> GeometriesArrayType;

    /** Redefinition of geometry template parameter TPointType as this geometry point type.
     */
    typedef TPointType PointType;

    /** Type used for indexing in geometry class.std::size_t used for indexing
    point or integration point access methods and also all other
    methods which need point or integration point index.
    */
    typedef std::size_t IndexType;


    /** This typed used to return size or dimension in
    geometry. Dimension, WorkingDimension, PointsNumber and
    ... return this type as their results.
    */
    typedef std::size_t SizeType;

    typedef typename PointType::CoordinatesArrayType CoordinatesArrayType;


    /** This type used for representing an integration point in
    geometry. This integration point is a point with an
    additional weight component.
    */
    typedef IntegrationPoint<3> IntegrationPointType;

    /** A Vector of IntegrationPointType which used to hold
    integration points related to an integration
    method. IntegrationPoints functions used this type to return
    their results.
    */
    typedef std::vector<IntegrationPointType> IntegrationPointsArrayType;

    /** A Vector of IntegrationPointsArrayType which used to hold
    integration points related to different integration method
    implemented in geometry.
    */
    typedef std::array<IntegrationPointsArrayType, GeometryData::NumberOfIntegrationMethods> IntegrationPointsContainerType;

    /** A third order tensor used as shape functions' values
    continer.
    */
    typedef std::array<Matrix, GeometryData::NumberOfIntegrationMethods> ShapeFunctionsValuesContainerType;

    /** A fourth order tensor used as shape functions' local
    gradients container in geometry.
    */
    typedef GeometryData::ShapeFunctionsLocalGradientsContainerType ShapeFunctionsLocalGradientsContainerType;

    /** A third order tensor to hold jacobian matrices evaluated at
    integration points. Jacobian and InverseOfJacobian functions
    return this type as their result.
    */
    typedef DenseVector<Matrix > JacobiansType;

    /** A third order tensor to hold shape functions'  gradients.
    ShapefunctionsGradients function return this
    type as its result.
    */
    typedef GeometryData::ShapeFunctionsGradientsType ShapeFunctionsGradientsType;

    /** A third order tensor to hold shape functions' local second derivatives.
    ShapefunctionsLocalGradients function return this
    type as its result.
    */
    typedef GeometryData::ShapeFunctionsSecondDerivativesType ShapeFunctionsSecondDerivativesType;

    /** A fourth order tensor to hold shape functions' local third order derivatives
     */
    typedef GeometryData::ShapeFunctionsThirdDerivativesType ShapeFunctionsThirdDerivativesType;

    /** Type of the normal vector used for normal to edges in geomety.
     */
    typedef DenseVector<double> NormalType;


    typedef typename BaseType::iterator              iterator;
    typedef typename BaseType::const_iterator          const_iterator;
    typedef typename BaseType::reverse_iterator        reverse_iterator;
    typedef typename BaseType::const_reverse_iterator  const_reverse_iterator;

    typedef typename BaseType::ptr_iterator ptr_iterator;
    typedef typename BaseType::ptr_const_iterator ptr_const_iterator;
    typedef typename BaseType::ptr_reverse_iterator ptr_reverse_iterator;
    typedef typename BaseType::ptr_const_reverse_iterator ptr_const_reverse_iterator;
    ///@}
    ///@name Life Cycle
    ///@{

    Geometry() : mpGeometryData(&GeometryDataInstance())
    {

    }

    /** Complete argument constructor. This constructor gives a
    complete set of arguments to pass all the initial value of
    all the member variables of geometry class. Also it has
    default value for integration variables to make it usefull
    in the case of constructing new geometry without mapping and
    integrating properties.

    @param ThisPoints Vector of pointers to points which this
    geometry constructing on them. Points must have dimension
    equal or greater than working space dimension though there
    is no control on it.

    @param ThisDefaultMethod Default integration method. Its
    default value is gaussian integration with orden one which
    make no deference while in this condition there is no shape
    function database exist and integrating is not possible
    including by default method.

    @param ThisIntegrationPoints All the integration points in
    all methods. This is a Vector of IntegrationPointsArrayType
    and It must have at least four component correspounding to
    four integration method defined now. If there is some
    geometry which don't have all this method implemented
    related points Vector must exist but with zero size. For
    example if a geometry don't have gaussian orden one
    ThisIntegrationPoints[GI_GAUSS_1] must be an empty
    IntegrationPointsArrayType.

    @param ThisShapeFunctionsValues Values of all the shape
    functions evaluated in all integrations points of all
    integration methods. It's a three dimensional array \f$
    F_{ijk} \f$ where i = GI_GAUSS_1,..., GI_GAUSS_4 and j is
    the integration point index and k is the shape function
    index. In the other word component \f$ f_{ijk} \f$ is the
    value of the shape function related to node k evaluated in
    integration point j of i integration method point set. Again
    if there is some integration method unsupported an empty
    Matrix must assigned to related place. For example if a
    geometry don't have gaussian orden four
    ThisShapeFunctionsValues[GI_GAUSS_4] must be an empty
    Matrix.

    @param ThisShapeFunctionsLocalGradients Values of local
    gradients respected to all local coordinates of all the
    shape functions evaluated in all integrations points of all
    integration methods. It's a four dimensional array \f$
    F_{ijkh} \f$ where i = GI_GAUSS_1,..., GI_GAUSS_4 and j is
    the integration point index and k is the shape function
    index and h is local coordinate index. In the other word
    component \f$ f_{ijkh} \f$ is the value of h'th component of
    local gradient of the shape function related to node k
    evaluated in integration point j of i integration method
    point set. Again if there is some integration method
    unsupported an empty ShapeFunctionsGradientsType must
    assigned to related place. For example if a geometry don't
    have gaussian orden two ThisShapeFunctionsValues[GI_GAUSS_2]
    must be an empty ShapeFunctionsGradientsType.
    */
    Geometry(const PointsArrayType &ThisPoints,
             GeometryData const *pThisGeometryData = &GeometryDataInstance())
        : BaseType(ThisPoints), mpGeometryData(pThisGeometryData)
    {
    }

//       Geometry(const PointsArrayType& ThisPoints,
//         GeometryData const& ThisGeometryData= msEmptyGeometryData)
//  : BaseType(ThisPoints)
//  , mGeometryData(ThisGeometryData)
//  {
//  }

    /** Copy constructor.
    Construct this geometry as a copy of given geometry.

    @note This copy constructor don't copy the points and new
    geometry shares points with given source geometry. It's
    obvious that any change to this new geometry's point affect
    source geometry's points too.
    */
    Geometry( const Geometry& rOther )
        : BaseType( rOther )
        , mpGeometryData( rOther.mpGeometryData )
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
    template<class TOtherPointType> Geometry( Geometry<TOtherPointType> const & rOther )
        : BaseType( rOther.begin(), rOther.end() )
        , mpGeometryData( rOther.mpGeometryData )
    {
    }

    /// Destructor. Do nothing!!!
    ~Geometry() override {}

    virtual GeometryData::KratosGeometryFamily GetGeometryFamily() const
    {
        return GeometryData::Kratos_generic_family;
    }

    virtual GeometryData::KratosGeometryType GetGeometryType() const
    {
        return GeometryData::Kratos_generic_type;
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
    Geometry& operator=( const Geometry& rOther )
    {
        BaseType::operator=( rOther );
        mpGeometryData = rOther.mpGeometryData;

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
    Geometry& operator=( Geometry<TOtherPointType> const & rOther )
    {
        this->clear();

        for ( typename Geometry<TOtherPointType>::ptr_const_iterator i = rOther.ptr_begin() ; i != rOther.ptr_end() ; ++i )
            push_back( typename PointType::Pointer( new PointType( **i ) ) );

        mpGeometryData = rOther.mpGeometryData;

        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    virtual Pointer Create( PointsArrayType const& ThisPoints ) const
    {
        return Pointer( new Geometry( ThisPoints, mpGeometryData ) );
    }

    /** This methods will create a duplicate of all its points and
    substitute them with its points. */
    void ClonePoints()
    {
        for ( ptr_iterator i = this->ptr_begin() ; i != this->ptr_end() ; i++ )
            *i = typename PointType::Pointer( new PointType( **i ) );
    }

    // virtual Kratos::shared_ptr< Geometry< Point > > Clone() const
    // {
    //     Geometry< Point >::PointsArrayType NewPoints;

    //     //making a copy of the nodes TO POINTS (not Nodes!!!)

    //     for ( IndexType i = 0 ; i < this->size() ; i++ )
    //     {
    //         NewPoints.push_back(Kratos::make_shared< Point >((*this)[i]));
    //     }

    //     //NewPoints[i] = typename Point::Pointer(new Point(*mPoints[i]));

    //     //creating a geometry with the new points
    //     Geometry< Point >::Pointer p_clone( new Geometry< Point >( NewPoints ) );

    //     p_clone->ClonePoints();

    //     return p_clone;
    // }

    //lumping factors for the calculation of the lumped mass matrix
    virtual Vector& LumpingFactors( Vector& rResult )  const
    {
        KRATOS_ERROR << "Called the virtual function for LumpingFactors " << *this << std::endl;
        return rResult;
    }

    ///@}
    ///@name Informations
    ///@{

    /** Dimension of the geometry for example a triangle2d is a 2
    dimensional shape

    @return SizeType, dimension of this geometry.
    @see WorkingSpaceDimension()
    @see LocalSpaceDimension()
    */
    inline SizeType Dimension() const
    {
        return mpGeometryData->Dimension();
    }

    /** Working space dimension. for example a triangle is a 2
    dimensional shape but can be used in 3 dimensional space.

    @return SizeType, working space dimension of this geometry.
    @see Dimension()
    @see LocalSpaceDimension()
    */
    inline SizeType WorkingSpaceDimension() const
    {
        return mpGeometryData->WorkingSpaceDimension();
    }

    /** Local space dimension. for example a triangle is a 2
    dimensional shape but can have 3 dimensional area
    coordinates l1, l2, l3.

    @return SizeType, local space dimension of this geometry.
    @see Dimension()
    @see WorkingSpaceDimension()
    */
    inline SizeType LocalSpaceDimension() const
    {
        return mpGeometryData->LocalSpaceDimension();
    }

    /** Returns number of the points which this geometry has.
     *
     * @return SizeType, number of the points in this geometry.
     */
    SizeType PointsNumber() const {
      return this->size();
    }

    /** This method calculate and return Length or charactereistic
     * length of this geometry depending to it's dimension. For one
     * dimensional geometry for example Line it returns length of it
     * and for the other geometries it gives Characteristic length
     * otherwise.
     *
     * @return double value contains length or Characteristic
     * length
     *
     * @see Area()
     * @see Volume()
     * @see DomainSize()
     */
    virtual double Length() const {
      KRATOS_ERROR << "Calling base class 'Length' method instead of derived class one. Please check the definition of derived class. " << *this << std::endl;
      return 0.0;
    }

    /** This method calculate and return area or surface area of
     * this geometry depending to it's dimension. For one dimensional
     * geometry it returns length, for two dimensional it gives area
     * and for three dimensional geometries it gives surface area.
     *
     * @return double value contains area or surface
     * area.
     *
     * @see Length()
     * @see Volume()
     * @see DomainSize()
     */
    virtual double Area() const {
      KRATOS_ERROR << "Calling base class 'Area' method instead of derived class one. Please check the definition of derived class. " << *this << std::endl;
      return 0.0;
    }

    /** This method calculate and return volume of this
     * geometry. For one and two dimensional geometry it returns
     * zero and for three dimensional it gives volume of geometry.
     *
     * @return double value contains volume.
     *
     * @see Length()
     * @see Area()
     * @see DomainSize()
     */
    virtual double Volume() const {
      KRATOS_ERROR << "Calling base class 'Volume' method instead of derived class one. Please check the definition of derived class. " << *this << std::endl;
      return 0.0;
    }

    /** This method calculate and return length, area or volume of
     * this geometry depending to it's dimension. For one dimensional
     * geometry it returns its length, for two dimensional it gives area
     * and for three dimensional geometries it gives its volume.
     *
     * @return double value contains length, area or volume.
     *
     * @see Length()
     * @see Area()
     * @see Volume()
     */
    virtual double DomainSize() const {
      KRATOS_ERROR << "Calling base class 'DomainSize' method instead of derived class one. Please check the definition of derived class. " << *this << std::endl;
      return 0.0;
    }

    /** This method calculates and returns the minimum edge.
     * length of the geometry.
     *
     * @return double value with the minimum edge length.
     *
     * @see MaxEdgeLength()
     * @see AverageEdgeLength()
     */
    virtual double MinEdgeLength() const {
      KRATOS_ERROR << "Calling base class 'MinEdgeLength' method instead of derived class one. Please check the definition of derived class. " << *this << std::endl;
      return 0.0;
    }

    /** This method calculates and returns the maximum edge.
     * length of the geometry.
     *
     * @return double value with the maximum edge length.
     *
     * @see MinEdgeLength()
     * @see AverageEdgeLength()
     */
    virtual double MaxEdgeLength() const {
      KRATOS_ERROR << "Calling base class 'MaxEdgeLength' method instead of derived class one. Please check the definition of derived class. " << *this << std::endl;
      return 0.0;
    }

    /** This method calculates and returns the average edge.
     * length of the geometry.
     *
     * @return double value with the average edge length
     *
     * @see MinEdgeLength()
     * @see MaxEdgeLength()
     */
    virtual double AverageEdgeLength() const {
      KRATOS_ERROR << "Calling base class 'AverageEdgeLength' method instead of derived class one. Please check the definition of derived class. " << *this << std::endl;
      return 0.0;
    }

    /** Calculates the circumradius of the geometry.
     * Calculates the circumradius of the geometry.
     *
     * @return Circumradius of the geometry.
     *
     * @see Inradius()
     */
    virtual double Circumradius() const {
      KRATOS_ERROR << "Calling base class 'Circumradius' method instead of derived class one. Please check the definition of derived class. " << *this << std::endl;
      return 0.0;
    }

    /** Calculates the inradius of the geometry.
     * Calculates the inradius of the geometry.
     *
     * @return Inradius of the geometry.
     *
     * @see Circumradius()
     */
    virtual double Inradius() const {
      KRATOS_ERROR << "Calling base class 'Inradius' method instead of derived class one. Please check the definition of derived class. " << *this << std::endl;
      return 0.0;
    }

    /** Test the intersection with another geometry
     *
     * Test if this geometry intersects with other geometry
     *
     * @param  ThisGeometry Geometry to intersect with
     * @return              True if the geometries intersect, False in any other case.
     */
    virtual bool HasIntersection(const GeometryType& ThisGeometry) {
      KRATOS_ERROR << "Calling base class 'HasIntersection' method instead of derived class one. Please check the definition of derived class. " << *this << std::endl;
      return false;
    }

    /** Test intersection of the geometry with a box
     *
     * Tests the intersection of the geometry with
     * a 3D box defined by rLowPoint and rHighPoint
     *
     * @param  rLowPoint  Lower point of the box to test the intersection
     * @param  rHighPoint Higher point of the box to test the intersection
     * @return            True if the geometry intersects the box, False in any other case.
     */
    virtual bool HasIntersection(const Point& rLowPoint, const Point& rHighPoint) {
      KRATOS_ERROR << "Calling base class 'HasIntersection' method instead of derived class one. Please check the definition of derived class. " << *this << std::endl;
      return false;
    }

    // virtual void BoundingBox(BoundingBoxType& rResult) const
    // {
    //
    //   Bounding_Box(rResult.LowPoint(), rResult.HighPoint());
    // }

    /** Calculates the boundingbox of the geometry.
     * Calculates the boundingbox of the geometry.
     *
     * @param rLowPoint  Lower point of the boundingbox.
     * @param rHighPoint Higher point of the boundingbox.
     */
    virtual void BoundingBox( TPointType& rLowPoint, TPointType& rHighPoint ) const
    {
        rHighPoint         = this->GetPoint( 0 );
        rLowPoint          = this->GetPoint( 0 );
        const SizeType dim = WorkingSpaceDimension();

        for ( unsigned int point = 0; point < PointsNumber(); point++ )
        {
            for ( unsigned int i = 0; i < dim; i++ )
            {
                rHighPoint[i] = ( rHighPoint[i] < this->GetPoint( point )[i] ) ? this->GetPoint( point )[i] : rHighPoint[i];
                rLowPoint[i]  = ( rLowPoint[i]  > this->GetPoint( point )[i] ) ? this->GetPoint( point )[i] : rLowPoint[i];
            }
        }
    }

    /** Calculates center of this geometry by a simple averaging algorithm.
    Each center point component calculated using:
    \f[
    c_i = \sum_j^n(x_i^j) / n
    \f]

    where \f$ c_i \f$ is component i of center point and \f$
    X_i^j \f$ is component i of j'th point of geometry and n is
    number of the points in this geometry.

    @return PointType which is the calculated center of this geometry.
    */
    virtual Point Center() const
    {
        const SizeType points_number = this->size();

        if ( points_number == 0 )
        {
            KRATOS_ERROR << "can not compute the ceneter of a geometry of zero points" << std::endl;
            // return PointType();
        }

        Point result = ( *this )[0];

        for ( IndexType i = 1 ; i < points_number ; i++ )
        {
            result.Coordinates() += ( *this )[i];
        }

        const double temp = 1.0 / double( points_number );

        result.Coordinates() *= temp;

        return result;
    }

    /**
     * It computes the unit normal of the geometry, if possible
     * @return The normal of the geometry
     */
    virtual array_1d<double, 3> AreaNormal(const CoordinatesArrayType& rPointLocalCoordinates) const
    {
        const unsigned int local_space_dimension = this->LocalSpaceDimension();
        const unsigned int dimension = this->WorkingSpaceDimension();

        if (dimension == local_space_dimension)
        {
            KRATOS_ERROR << "Remember the normal can be computed just in geometries with a local dimension: "<< this->LocalSpaceDimension() << "smaller than the spatial dimension: " << this->WorkingSpaceDimension() << std::endl;
        }

        // We define the normal and tangents
        array_1d<double,3> tangent_xi(3, 0.0);
        array_1d<double,3> tangent_eta(3, 0.0);

        Matrix j_node = ZeroMatrix( dimension, local_space_dimension );
        this->Jacobian( j_node, rPointLocalCoordinates);

        // Using the Jacobian tangent directions
        if (dimension == 2)
        {
            tangent_eta[2] = 1.0;
            for (unsigned int i_dim = 0; i_dim < dimension; i_dim++)
            {
                tangent_xi[i_dim]  = j_node(i_dim, 0);
            }
        }
        else
        {
            for (unsigned int i_dim = 0; i_dim < dimension; i_dim++)
            {
                tangent_xi[i_dim]  = j_node(i_dim, 0);
                tangent_eta[i_dim] = j_node(i_dim, 1);
            }
        }

        array_1d<double, 3> normal;
        MathUtils<double>::CrossProduct(normal, tangent_xi, tangent_eta);
        return normal;
    }

    /**
     * It computes the unit normal of the geometry
     * @param rPointLocalCoordinates Refernce to the local coordinates of the
     * point in where the unit normal is to be computed
     * @return The unit normal in the given point
     */
    virtual array_1d<double, 3> UnitNormal(const CoordinatesArrayType& rPointLocalCoordinates) const
    {
        array_1d<double, 3> normal = AreaNormal(rPointLocalCoordinates);
        const double norm_normal = norm_2(normal);
        if (norm_normal > std::numeric_limits<double>::epsilon()) normal /= norm_normal;
        else KRATOS_ERROR << "ERROR: The normal norm is zero or almost zero. Norm. normal: " << norm_normal << std::endl;
        return normal;
    }

    /** Calculates the quality of the geometry according to a given criteria.
     *
     * Calculates the quality of the geometry according to a given criteria. In General
     * The quality of the result is normalized being 1.0 for best quality, 0.0 for degenerated elements and -1.0 for
     * inverted elements.
     *
     * Different crtieria can be used to stablish the quality of the geometry.
     *
     * @return double value contains quality of the geometry
     *
     * @see QualityCriteria
     * @see QualityAspectRatio
     * @see QualityAverageEdgeLenght
     */
     double Quality(const QualityCriteria qualityCriteria) const {
       double quality = 0.0f;

       if(qualityCriteria == QualityCriteria::INRADIUS_TO_CIRCUMRADIUS) {
         quality = InradiusToCircumradiusQuality();
       } else if(qualityCriteria == QualityCriteria::AREA_TO_LENGTH) {
         quality = AreaToEdgeLengthRatio();
       } else if(qualityCriteria == QualityCriteria::SHORTEST_ALTITUDE_TO_LENGTH) {
         quality = ShortestAltitudeToEdgeLengthRatio();
       } else if(qualityCriteria == QualityCriteria::INRADIUS_TO_LONGEST_EDGE) {
         quality = InradiusToLongestEdgeQuality();
       } else if(qualityCriteria == QualityCriteria::SHORTEST_TO_LONGEST_EDGE) {
         quality = ShortestToLongestEdgeQuality();
       } else if(qualityCriteria == QualityCriteria::REGULARITY) {
         quality = RegularityQuality();
       } else if(qualityCriteria == QualityCriteria::VOLUME_TO_SURFACE_AREA) {
         quality = VolumeToSurfaceAreaQuality();
       } else if(qualityCriteria == QualityCriteria::VOLUME_TO_EDGE_LENGTH) {
         quality = VolumeToEdgeLengthQuality();
       } else if(qualityCriteria == QualityCriteria::VOLUME_TO_AVERAGE_EDGE_LENGTH) {
         quality = VolumeToAverageEdgeLength();
       } else if(qualityCriteria == QualityCriteria::VOLUME_TO_RMS_EDGE_LENGTH) {
         quality = VolumeToRMSEdgeLength();
       }

       return quality;
     }

    ///@}
    ///@name Access
    ///@{

    /** A constant access method to the Vector of the points stored in
    this geometry.

    @return A constant reference to PointsArrayType contains
    pointers to the points.
    */
    const PointsArrayType& Points() const
    {
        return *this;
    }

    /** An access method to the Vector of the points stored in
    this geometry.

    @return A reference to PointsArrayType contains pointers to
    the points.
    */
    PointsArrayType& Points()
    {
        return *this;
    }

    /** A constant access method to the i'th points stored in
    this geometry.

    @return A constant counted pointer to i'th point of
    geometry.
    */
    const typename TPointType::Pointer pGetPoint( const int Index ) const
    {
        KRATOS_TRY
        return ( *this )( Index );
        KRATOS_CATCH( *this )
    }

    /** An access method to the i'th points stored in
    this geometry.

    @return A counted pointer to i'th point of
    geometry.
    */
    typename TPointType::Pointer pGetPoint( const int Index )
    {
        KRATOS_TRY
        return ( *this )( Index );
        KRATOS_CATCH( *this );
    }

    /** A constant access method to the i'th points stored in
    this geometry.

    @return A constant counted pointer to i'th point of
    geometry.
    */
    TPointType const& GetPoint( const int Index ) const
    {
        KRATOS_TRY
        return ( *this )[Index];
        KRATOS_CATCH( *this );
    }


    /** An access method to the i'th points stored in
    this geometry.

    @return A counted pointer to i'th point of
    geometry.
    */
    TPointType& GetPoint( const int Index )
    {
        KRATOS_TRY
        return ( *this )[Index];
        KRATOS_CATCH( *this );
    }

    /**
     * Returns a matrix of the local coordinates of all points
     * @param rResult a Matrix that will be overwritten by the results
     * @return the coordinates of all points of the current geometry
     */
    virtual Matrix& PointsLocalCoordinates( Matrix& rResult ) const
    {
        KRATOS_ERROR << "Calling base class 'PointsLocalCoordinates' method instead of derived class one. Please check the definition of derived class. " << *this << std::endl;
        return rResult;
    }

    /**
     * @brief Returns the local coordinates of a given arbitrary point
     * @param rResult The vector containing the local coordinates of the point
     * @param rPoint The point in global coordinates
     * @return The vector containing the local coordinates of the point
     */
    virtual CoordinatesArrayType& PointLocalCoordinates(
            CoordinatesArrayType& rResult,
            const CoordinatesArrayType& rPoint
            ) const
    {
        KRATOS_ERROR_IF(WorkingSpaceDimension() != LocalSpaceDimension()) << "ERROR:: Attention, the Point Local Coordinates must be specialized for the current geometry" << std::endl;

        Matrix J = ZeroMatrix( WorkingSpaceDimension(), LocalSpaceDimension() );

        rResult.clear();

        Vector DeltaXi = ZeroVector( LocalSpaceDimension() );

        CoordinatesArrayType CurrentGlobalCoords( ZeroVector( 3 ) );

        static constexpr double MaxNormPointLocalCoordinates = 30.0;
        static constexpr std::size_t MaxIteratioNumberPointLocalCoordinates = 1000;
        static constexpr double MaxTolerancePointLocalCoordinates = 1.0e-8;

        //Newton iteration:
        for(std::size_t k = 0; k < MaxIteratioNumberPointLocalCoordinates; k++) {
            CurrentGlobalCoords.clear();
            DeltaXi.clear();

            GlobalCoordinates( CurrentGlobalCoords, rResult );
            noalias( CurrentGlobalCoords ) = rPoint - CurrentGlobalCoords;
            InverseOfJacobian( J, rResult );
            for(unsigned int i = 0; i < WorkingSpaceDimension(); i++) {
                for(unsigned int j = 0; j < WorkingSpaceDimension(); j++) {
                    DeltaXi[i] += J(i,j)*CurrentGlobalCoords[j];
                }
                rResult[i] += DeltaXi[i];
            }

            const double norm2DXi = norm_2(DeltaXi);

            if(norm2DXi > MaxNormPointLocalCoordinates) {
                KRATOS_WARNING("Geometry") << "Computation of local coordinates failed at iteration " << k << std::endl;
                break;
            }

            if(norm2DXi < MaxTolerancePointLocalCoordinates) {
                break;
            }
        }

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
    virtual bool IsInside(
        const CoordinatesArrayType& rPoint,
        CoordinatesArrayType& rResult,
        const double Tolerance = std::numeric_limits<double>::epsilon()
        )
    {
        KRATOS_ERROR << "Calling base class IsInside method instead of derived class one. Please check the definition of derived class. " << *this << std::endl;
        return false;
    }

    ///@}
    ///@name Inquiry
    ///@{

    /** This method confirm you if this geometry has a specific
    integration method or not. This method will be usefull to
    control the geometry before intagrating using a specific
    method. In Geometry class this method controls if the
    integration points vector respecting to this method is empty
    or not.

    @return bool true if this integration method exist and false if this
    method is not imeplemented for this geometry.
    */
    bool HasIntegrationMethod( IntegrationMethod ThisMethod ) const
    {
        return ( mpGeometryData->HasIntegrationMethod( ThisMethod ) );
    }

    /**
    * @return default integration method
    */

    IntegrationMethod GetDefaultIntegrationMethod() const
    {
        return mpGeometryData->DefaultIntegrationMethod();
    }

    /** This method is to know if this geometry is symmetric or
    not.

    @todo Making some method related to symmetry axis and more...

    @return bool true if this geometry is symmetric and false if
    it's not.
    */
    virtual bool IsSymmetric() const
    {
        return false;
    }

    ///@}
    ///@name Edge
    ///@{

    /** This method gives you number of all edges of this
     * geometry. For example, for a hexahedron, this would be
     * 12
     *
     * @return SizeType containes number of this geometry edges.
     * @see EdgesNumber()
     * @see Edges()
     * @see FacesNumber()
     * @see Faces()
     */
    // will be used by refinement algorithm, thus uncommented. janosch.
    virtual SizeType EdgesNumber() const
    {
        KRATOS_ERROR << "Calling base class EdgesNumber method instead of derived class one. Please check the definition of derived class. " << *this << std::endl;

        return SizeType();
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
    // will be used by refinement algorithm, thus uncommented. janosch.
    virtual GeometriesArrayType Edges( void )
    {
        KRATOS_ERROR << "Calling base class Edges method instead of derived class one. Please check the definition of derived class. " << *this << std::endl;

        return GeometriesArrayType();
    }

    /** This method gives you an edge of this geometry which holds
    given points. This method will gives you an edge with
    dimension related to given points number. for example a
    tetrahedral would return a triangle for given three points or
    return an edge line for given two nodes by this method.

    @return Geometry::Pointer to this geometry specific edge.
    @see EdgesNumber()
    @see Edges()
    */
    // Commented for possible change in Edge interface of geometry. Pooyan.
//       virtual Pointer Edge(const PointsArrayType& EdgePoints)
//  {
//    KRATOS_ERROR << "Calling base class Edge method instead of derived class one. Please check the definition of derived class." << *this << std::endl;
//
//  }

    /**
     * Returns the number of faces of the current geometry.
     * This is only implemented for 3D geometries, since 2D geometries
     * only have edges but no faces
     * @see EdgesNumber
     * @see Edges
     * @see Faces
     */
    virtual SizeType FacesNumber() const
    {
        KRATOS_ERROR << "Calling base class FacesNumber method instead of derived class one. Please check the definition of derived class. " << *this << std::endl;

        return SizeType();
    }

    /**
     * Returns all faces of the current geometry.
     * This is only implemented for 3D geometries, since 2D geometries
     * only have edges but no faces
     * @see EdgesNumber
     * @see Edges
     * @see FacesNumber
     */
    virtual GeometriesArrayType Faces( void )
    {
        KRATOS_ERROR << "Calling base class Faces method instead of derived class one. Please check the definition of derived class. " << *this << std::endl;

        return GeometriesArrayType();
    }

    //Connectivities of faces required
    virtual void NumberNodesInFaces (DenseVector<unsigned int>& rNumberNodesInFaces) const
    {
        KRATOS_ERROR << "Calling base class NumberNodesInFaces method instead of derived class one. Please check the definition of derived class. " << *this << std::endl;
    }

    virtual void NodesInFaces (DenseMatrix<unsigned int>& rNodesInFaces) const
    {
        KRATOS_ERROR << "Calling base class NodesInFaces method instead of derived class one. Please check the definition of derived class. " << *this << std::endl;
    }


    /** This method gives you an edge of this geometry related to
    given index. The numbering order of each geometries edges is
    depended to type of that geometry.

    @return Geometry::Pointer to this geometry specific edge.
    @see EdgesNumber()
    @see Edges()
    */
    // Commented for possible change in Edge interface of geometry. Pooyan.
//       virtual Pointer Edge(IndexType EdgeIndex)
//  {
//    KRATOS_ERROR << "Calling base class Edge method instead of derived class one. Please check the definition of derived class." << *this << std::endl;

//  }

    /** This method gives you normal edge of this geometry which holds
    given points.

    @return NormalType which is normal to this geometry specific edge.
    @see Edge()
    */
    // Commented for possible change in Edge interface of geometry. Pooyan.
//       virtual NormalType NormalEdge(const PointsArrayType& EdgePoints)
//  {
//    KRATOS_ERROR << "Calling base class NormalEdge method instead of derived class one. Please check the definition of derived class." << *this << std::endl

//    return NormalType();
//  }

    /** This method gives you normal to edge of this geometry related to
    given index. The numbering order of each geometries edges is
    depended to type of that geometry.

    @return NormalType which is normal to this geometry specific edge.
    @see Edge()
    */
    // Commented for possible change in Edge interface of geometry. Pooyan.
//       virtual NormalType NormalEdge(IndexType EdgeIndex)
//  {
//    KRATOS_ERROR << "Calling base class NormalEdge method instead of derived class one. Please check the definition of derived class." << *this << std::endl;

//    return NormalType();
//  }

    ///@}
    ///@name Integration Points
    ///@{

    /** Number of integtation points for default integration
    method. This method just call IntegrationPointsNumber(enum
    IntegrationMethod ThisMethod) with default integration
    method.

    @return SizeType which is the number of integration points
    for default integrating method.
    */
    SizeType IntegrationPointsNumber() const
    {
        return mpGeometryData->IntegrationPoints().size();
    }

    /** Number of integtation points for given integration
    method. This method use integration points data base to
    obtain size of the integration points Vector respected to
    given method.

    @return SizeType which is the number of integration points
    for given integrating method.
    */
    SizeType IntegrationPointsNumber( IntegrationMethod ThisMethod ) const
    {
        return mpGeometryData->IntegrationPointsNumber( ThisMethod );
    }


    /** Integtation points for default integration
    method. This method just call IntegrationPoints(enum
    IntegrationMethod ThisMethod) with default integration
    method.

    @return const IntegrationPointsArrayType which is Vector of integration points
    for default integrating method.
    */
    const IntegrationPointsArrayType& IntegrationPoints() const
    {
        return mpGeometryData->IntegrationPoints();
    }

    /** Integtation points for given integration
    method. This method use integration points data base to
    obtain integration points Vector respected to
    given method.

    @return const IntegrationPointsArrayType which is Vector of integration points
    for default integrating method.
    */
    const IntegrationPointsArrayType& IntegrationPoints( IntegrationMethod ThisMethod ) const
    {
        return mpGeometryData->IntegrationPoints( ThisMethod );
    }

    ///@}
    ///@name Jacobian
    ///@{

    /** This method provides the global coordinates corresponding to the local coordinates provided
     * @param rResult The array containing the global coordinates corresponding to the local coordinates provided
     * @param LocalCoordinates The local coordinates provided
     * @return An array containing the global coordinates corresponding to the local coordinates provides
     * @see PointLocalCoordinates
     */
    virtual CoordinatesArrayType& GlobalCoordinates(
        CoordinatesArrayType& rResult,
        CoordinatesArrayType const& LocalCoordinates
        ) const
    {
        noalias( rResult ) = ZeroVector( 3 );

        Vector N( this->size() );
        ShapeFunctionsValues( N, LocalCoordinates );

        for ( IndexType i = 0 ; i < this->size() ; i++ )
            noalias( rResult ) += N[i] * (*this)[i];

        return rResult;
    }

    /** This method provides the global coordinates corresponding to the local coordinates provided, considering additionally a certain increment in the coordinates
     * @param rResult The array containing the global coordinates corresponding to the local coordinates provided
     * @param LocalCoordinates The local coordinates provided
     * @param DeltaPosition The increment of position considered
     * @return An array containing the global coordinates corresponding to the local coordinates provides
     * @see PointLocalCoordinates
     */
    virtual CoordinatesArrayType& GlobalCoordinates(
        CoordinatesArrayType& rResult,
        CoordinatesArrayType const& LocalCoordinates,
        Matrix& DeltaPosition
        ) const
    {
        constexpr std::size_t dimension = 3;
        noalias( rResult ) = ZeroVector( 3 );
        if (DeltaPosition.size2() != 3)
            DeltaPosition.resize(DeltaPosition.size1(), dimension,false);

        Vector N( this->size() );
        ShapeFunctionsValues( N, LocalCoordinates );

        for ( IndexType i = 0 ; i < this->size() ; i++ )
            noalias( rResult ) += N[i] * ((*this)[i] + row(DeltaPosition, i));

        return rResult;
    }

    /** Jacobians for default integration method. This method just
    call Jacobian(enum IntegrationMethod ThisMethod) with
    default integration method.

    @return JacobiansType a Vector of jacobian
    matrices \f$ J_i \f$ where \f$ i=1,2,...,n \f$ is the integration
    point index of default integration method.

    @see DeterminantOfJacobian
    @see InverseOfJacobian
    */
    JacobiansType& Jacobian( JacobiansType& rResult ) const
    {
        Jacobian( rResult, mpGeometryData->DefaultIntegrationMethod() );
        return rResult;
    }

    /** Jacobians for given  method. This method
    calculate jacobians matrices in all integrations points of
    given integration method.

    @param ThisMethod integration method which jacobians has to
    be calculated in its integration points.

    @return JacobiansType a Vector of jacobian
    matrices \f$ J_i \f$ where \f$ i=1,2,...,n \f$ is the integration
    point index of given integration method.

    @see DeterminantOfJacobian
    @see InverseOfJacobian
    */
    virtual JacobiansType& Jacobian( JacobiansType& rResult,
                                     IntegrationMethod ThisMethod ) const
    {
        if( rResult.size() != this->IntegrationPointsNumber( ThisMethod ) )
            rResult.resize( this->IntegrationPointsNumber( ThisMethod ), false );

        for ( unsigned int pnt = 0; pnt < this->IntegrationPointsNumber( ThisMethod ); pnt++ )
        {
            this->Jacobian( rResult[pnt], pnt, ThisMethod);
        }

        return rResult;
    }

    /** Jacobians for given  method. This method
    calculate jacobians matrices in all integrations points of
    given integration method.

    @param ThisMethod integration method which jacobians has to
    be calculated in its integration points.

    @return JacobiansType a Vector of jacobian
    matrices \f$ J_i \f$ where \f$ i=1,2,...,n \f$ is the integration
    point index of given integration method.

    @param DeltaPosition Matrix with the nodes position increment which describes
    the configuration where the jacobian has to be calculated.

    @see DeterminantOfJacobian
    @see InverseOfJacobian
    */
    virtual JacobiansType& Jacobian( JacobiansType& rResult, IntegrationMethod ThisMethod, Matrix & DeltaPosition ) const
    {
        if( rResult.size() != this->IntegrationPointsNumber( ThisMethod ) )
            rResult.resize( this->IntegrationPointsNumber( ThisMethod ), false );

        for ( unsigned int pnt = 0; pnt < this->IntegrationPointsNumber( ThisMethod ); pnt++ )
        {
            this->Jacobian( rResult[pnt], pnt, ThisMethod, DeltaPosition);
        }
        return rResult;
    }

    /** Jacobian in specific integration point of default integration method. This method just
    call Jacobian(IndexType IntegrationPointIndex, enum IntegrationMethod ThisMethod) with
    default integration method.

    @param IntegrationPointIndex index of integration point which jacobians has to
    be calculated in it.

    @return Matrix<double> Jacobian matrix \f$ J_i \f$ where \f$
    i \f$ is the given integration point index of default
    integration method.

    @see DeterminantOfJacobian
    @see InverseOfJacobian
    */
    Matrix& Jacobian( Matrix& rResult, IndexType IntegrationPointIndex ) const
    {
        Jacobian( rResult, IntegrationPointIndex, mpGeometryData->DefaultIntegrationMethod() );
        return rResult;
    }

    /** Jacobian in specific integration point of given integration
    method. This method calculate jacobian matrix in given
    integration point of given integration method.

    @param IntegrationPointIndex index of integration point which jacobians has to
    be calculated in it.

    @param ThisMethod integration method which jacobians has to
    be calculated in its integration points.

    @return Matrix<double> Jacobian matrix \f$ J_i \f$ where \f$
    i \f$ is the given integration point index of given
    integration method.

    @see DeterminantOfJacobian
    @see InverseOfJacobian
    */
    virtual Matrix& Jacobian( Matrix& rResult, IndexType IntegrationPointIndex, IntegrationMethod ThisMethod ) const
    {
        if(rResult.size1() != this->WorkingSpaceDimension() || rResult.size2() != this->LocalSpaceDimension())
            rResult.resize( this->WorkingSpaceDimension(), this->LocalSpaceDimension(), false );

        const Matrix& ShapeFunctionsGradientInIntegrationPoint = ShapeFunctionsLocalGradients( ThisMethod )[ IntegrationPointIndex ];

        rResult.clear();
        for ( unsigned int i = 0; i < this->PointsNumber(); i++ )
        {
            for(unsigned int k=0; k<this->WorkingSpaceDimension(); k++)
            {
                for(unsigned int m=0; m<this->LocalSpaceDimension(); m++)
                {
                    rResult(k,m) += (( *this )[i]).Coordinates()[k]*ShapeFunctionsGradientInIntegrationPoint(i,m);
                }
            }
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
    virtual Matrix& Jacobian( Matrix& rResult, IndexType IntegrationPointIndex, IntegrationMethod ThisMethod, Matrix& DeltaPosition ) const
    {
        if(rResult.size1() != this->WorkingSpaceDimension() || rResult.size2() != this->LocalSpaceDimension())
            rResult.resize( this->WorkingSpaceDimension(), this->LocalSpaceDimension(), false );

        const Matrix& ShapeFunctionsGradientInIntegrationPoint = ShapeFunctionsLocalGradients( ThisMethod )[ IntegrationPointIndex ];

        rResult.clear();
        for ( unsigned int i = 0; i < this->PointsNumber(); i++ )
        {
            for(unsigned int k=0; k<this->WorkingSpaceDimension(); k++)
            {
                for(unsigned int m=0; m<this->LocalSpaceDimension(); m++)
                {
                    rResult(k,m) += ( (( *this )[i]).Coordinates()[k]  - DeltaPosition(i,k)  )*ShapeFunctionsGradientInIntegrationPoint(i,m);
                }
            }
        }

        return rResult;
    }

    /** Jacobian in given point. This method calculate jacobian
    matrix in given point.

    @param rCoordinates point which jacobians has to
    be calculated in it.

    @return Matrix of double which is jacobian matrix \f$ J \f$ in given point.

    @see DeterminantOfJacobian
    @see InverseOfJacobian
    */
    virtual Matrix& Jacobian( Matrix& rResult, const CoordinatesArrayType& rCoordinates ) const
    {
        if(rResult.size1() != this->WorkingSpaceDimension() || rResult.size2() != this->LocalSpaceDimension())
            rResult.resize( this->WorkingSpaceDimension(), this->LocalSpaceDimension(), false );

        Matrix shape_functions_gradients(this->PointsNumber(), this->LocalSpaceDimension());
        ShapeFunctionsLocalGradients( shape_functions_gradients, rCoordinates );

        rResult.clear();
        for ( unsigned int i = 0; i < this->PointsNumber(); i++ )
        {
            for(unsigned int k=0; k<this->WorkingSpaceDimension(); k++)
            {
                for(unsigned int m=0; m<this->LocalSpaceDimension(); m++)
                {
                    rResult(k,m) += (( *this )[i]).Coordinates()[k]*shape_functions_gradients(i,m);
                }
            }
        }

        return rResult;
    }

    /** Jacobian in given point. This method calculate jacobian
    matrix in given point.

    @param rCoordinates point which jacobians has to
    be calculated in it.

    @param DeltaPosition Matrix with the nodes position increment which describes
    the configuration where the jacobian has to be calculated.

    @return Matrix of double which is jacobian matrix \f$ J \f$ in given point.

    @see DeterminantOfJacobian
    @see InverseOfJacobian
    */

    virtual Matrix& Jacobian( Matrix& rResult, const CoordinatesArrayType& rCoordinates, Matrix& DeltaPosition ) const
    {
        if(rResult.size1() != this->WorkingSpaceDimension() || rResult.size2() != this->LocalSpaceDimension())
            rResult.resize( this->WorkingSpaceDimension(), this->LocalSpaceDimension(), false );

        Matrix shape_functions_gradients(this->PointsNumber(), this->LocalSpaceDimension());
        ShapeFunctionsLocalGradients( shape_functions_gradients, rCoordinates );

        rResult.clear();
        for ( unsigned int i = 0; i < this->PointsNumber(); i++ )
        {
            for(unsigned int k=0; k<this->WorkingSpaceDimension(); k++)
            {
                for(unsigned int m=0; m<this->LocalSpaceDimension(); m++)
                {
                    rResult(k,m) += ( (( *this )[i]).Coordinates()[k] - DeltaPosition(i,k))*shape_functions_gradients(i,m);
                }
            }
        }

        return rResult;
    }

    /** Determinant of jacobians for default integration method. This method just
    call DeterminantOfJacobian(enum IntegrationMethod ThisMethod) with
    default integration method.

    @return Vector of double which is vector of determinants of
    jacobians \f$ |J|_i \f$ where \f$ i=1,2,...,n \f$ is the
    integration point index of default integration method.

    @see Jacobian
    @see InverseOfJacobian
    */
    Vector& DeterminantOfJacobian( Vector& rResult ) const
    {
        DeterminantOfJacobian( rResult, mpGeometryData->DefaultIntegrationMethod() );
        return rResult;
    }

    /** Determinant of jacobians for given integration method. This
    method calculate determinant of jacobian in all
    integrations points of given integration method.

    @return Vector of double which is vector of determinants of
    jacobians \f$ |J|_i \f$ where \f$ i=1,2,...,n \f$ is the
    integration point index of given integration method.

    @see Jacobian
    @see InverseOfJacobian
    */
    virtual Vector& DeterminantOfJacobian( Vector& rResult, IntegrationMethod ThisMethod ) const
    {
        if( rResult.size() != this->IntegrationPointsNumber( ThisMethod ) )
            rResult.resize( this->IntegrationPointsNumber( ThisMethod ), false );

        Matrix J;
        for ( unsigned int pnt = 0; pnt < this->IntegrationPointsNumber( ThisMethod ); pnt++ )
        {
            this->Jacobian( J, pnt, ThisMethod);
            rResult[pnt] = MathUtils<double>::GeneralizedDet(J);
        }
        return rResult;
    }

    /** Determinant of jacobian in specific integration point of
    default integration method. This method just call
    DeterminantOfJacobian(IndexType IntegrationPointIndex, enum
    IntegrationMethod ThisMethod) with default integration
    method.

    @param IntegrationPointIndex index of integration point
    which determinant jacobians has to be calculated in it.

    @return Determinamt of jacobian matrix \f$ |J|_i \f$ where \f$
    i \f$ is the given integration point index of default
    integration method.

    @see Jacobian
    @see InverseOfJacobian
    */
    double DeterminantOfJacobian( IndexType IntegrationPointIndex ) const
    {
        return DeterminantOfJacobian( IntegrationPointIndex, mpGeometryData->DefaultIntegrationMethod() );
    }

    /** Determinant of jacobian in specific integration point of
    given integration method. This method calculate determinant
    of jacobian in given integration point of given integration
    method.

    @param IntegrationPointIndex index of integration point which jacobians has to
    be calculated in it.

    @param IntegrationPointIndex index of integration point
    which determinant of jacobians has to be calculated in it.

    @return Determinamt of jacobian matrix \f$ |J|_i \f$ where \f$
    i \f$ is the given integration point index of given
    integration method.

    @see Jacobian
    @see InverseOfJacobian
    */
    virtual double DeterminantOfJacobian( IndexType IntegrationPointIndex, IntegrationMethod ThisMethod ) const
    {
        Matrix J;
        this->Jacobian( J, IntegrationPointIndex, ThisMethod);
        return MathUtils<double>::GeneralizedDet(J);
    }


    /** Determinant of jacobian in given point. This method calculate determinant of jacobian
    matrix in given point.

    @param rPoint point which determinant of jacobians has to
    be calculated in it.

    @return Determinamt of jacobian matrix \f$ |J| \f$ in given
    point.

    @see DeterminantOfJacobian
    @see InverseOfJacobian
    */
    virtual double DeterminantOfJacobian( const CoordinatesArrayType& rPoint ) const
    {
        Matrix J;
        this->Jacobian( J, rPoint);
        return MathUtils<double>::GeneralizedDet(J);
    }


    /** Inverse of jacobians for default integration method. This method just
    call InverseOfJacobian(enum IntegrationMethod ThisMethod) with
    default integration method.

    @return Inverse of jacobian
    matrices \f$ J_i^{-1} \f$ where \f$ i=1,2,...,n \f$ is the integration
    point index of default integration method.

    @see Jacobian
    @see DeterminantOfJacobian
    */
    JacobiansType& InverseOfJacobian( JacobiansType& rResult ) const
    {
        InverseOfJacobian( rResult, mpGeometryData->DefaultIntegrationMethod() );
        return rResult;
    }

    /** Inverse of jacobians for given integration method. This method
    calculate inverse of jacobians matrices in all integrations points of
    given integration method.

    @param ThisMethod integration method which inverse of jacobians has to
    be calculated in its integration points.

    @return Inverse of jacobian
    matrices \f$ J^{-1}_i \f$ where \f$ i=1,2,...,n \f$ is the integration
    point index of given integration method.

    @see Jacobian
    @see DeterminantOfJacobian
    */
    virtual JacobiansType& InverseOfJacobian( JacobiansType& rResult, IntegrationMethod ThisMethod ) const
    {
        Jacobian(rResult, ThisMethod); //this will be overwritten

        double detJ;
        Matrix Jinv(this->WorkingSpaceDimension(), this->WorkingSpaceDimension());
        for ( unsigned int pnt = 0; pnt < this->IntegrationPointsNumber( ThisMethod ); pnt++ )
        {
            MathUtils<double>::GeneralizedInvertMatrix(rResult[pnt], Jinv, detJ);
            noalias(rResult[pnt]) = Jinv;
        }
        return rResult;
    }

    /** Inverse of jacobian in specific integration point of default integration method. This method just
    call InverseOfJacobian(IndexType IntegrationPointIndex, enum IntegrationMethod ThisMethod) with
    default integration method.

    @param IntegrationPointIndex index of integration point which inverse of jacobians has to
    be calculated in it.

    @return Inverse of jacobian matrix \f$ J^{-1}_i \f$ where \f$
    i \f$ is the given integration point index of default
    integration method.

    @see Jacobian
    @see DeterminantOfJacobian
    */
    Matrix& InverseOfJacobian( Matrix& rResult, IndexType IntegrationPointIndex ) const
    {
        InverseOfJacobian( rResult, IntegrationPointIndex, mpGeometryData->DefaultIntegrationMethod() );
        return rResult;
    }

    /** Inverse of jacobian in specific integration point of given integration
    method. This method calculate Inverse of jacobian matrix in given
    integration point of given integration method.

    @param IntegrationPointIndex index of integration point which inverse of jacobians has to
    be calculated in it.

    @param ThisMethod integration method which inverse of jacobians has to
    be calculated in its integration points.

    @return Inverse of jacobian matrix \f$ J^{-1}_i \f$ where \f$
    i \f$ is the given integration point index of given
    integration method.

    @see Jacobian
    @see DeterminantOfJacobian
    */
    virtual Matrix& InverseOfJacobian( Matrix& rResult, IndexType IntegrationPointIndex, IntegrationMethod ThisMethod ) const
    {
        Jacobian(rResult,IntegrationPointIndex, ThisMethod); //this will be overwritten

        double detJ;
        Matrix Jinv(this->WorkingSpaceDimension(), this->WorkingSpaceDimension());

        MathUtils<double>::GeneralizedInvertMatrix(rResult, Jinv, detJ);
        noalias(rResult) = Jinv;

        return rResult;
    }

    /** Inverse of jacobian in given point. This method calculate inverse of jacobian
    matrix in given point.

    @param rPoint point which inverse of jacobians has to
    be calculated in it.

    @return Inverse of jacobian matrix \f$ J^{-1} \f$ in given point.

    @see DeterminantOfJacobian
    @see InverseOfJacobian
    */
    virtual Matrix& InverseOfJacobian( Matrix& rResult, const CoordinatesArrayType& rCoordinates ) const
    {
        Jacobian(rResult,rCoordinates); //this will be overwritten

        double detJ;
        Matrix Jinv(this->WorkingSpaceDimension(), this->WorkingSpaceDimension());

        MathUtils<double>::GeneralizedInvertMatrix(rResult, Jinv, detJ);
        noalias(rResult) = Jinv;

        return rResult;
    }



    ///@}
    ///@name Shape Function
    ///@{

    /** This method gives all shape functions values evaluated in all
    integration points of default integration method. It just
    call ShapeFunctionsValues(enum IntegrationMethod ThisMethod)
    with default integration method.There is no calculation and
    it just give it from shape functions values container.

    \note There is no control if the return matrix is empty or not!

    @return Matrix of values of shape functions \f$ F_{ij} \f$
    where i is the integration point index and j is the shape
    function index. In other word component \f$ f_{ij} \f$ is value
    of the shape function corresponding to node j evaluated in
    integration point i of default integration method.

    @see ShapeFunctionValue
    @see ShapeFunctionsLocalGradients
    @see ShapeFunctionLocalGradient
    */
    const Matrix& ShapeFunctionsValues() const
    {
        return mpGeometryData->ShapeFunctionsValues();
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

    virtual Vector& ShapeFunctionsValues (Vector &rResult, const CoordinatesArrayType& rCoordinates) const
    {
        KRATOS_ERROR << "Calling base class ShapeFunctionsValues method instead of derived class one. Please check the definition of derived class. " << *this << std::endl;
        return rResult;
    }

    /** This method gives all shape functions values evaluated in all
    integration points of given integration method. There is no
    calculation and it just give it from shape functions values
    container.

    \note There is no control if the return matrix is empty or not!

    @param ThisMethod integration method which shape functions
    evaluated in its integration points.

    @return Matrix of values of shape functions \f$ F_{ij} \f$
    where i is the integration point index and j is the shape
    function index. In other word component \f$ f_{ij} \f$ is value
    of the shape function corresponding to node j evaluated in
    integration point i of given integration method.

    @see ShapeFunctionValue
    @see ShapeFunctionsLocalGradients
    @see ShapeFunctionLocalGradient
    */
    const Matrix& ShapeFunctionsValues( IntegrationMethod ThisMethod )  const
    {
        return mpGeometryData->ShapeFunctionsValues( ThisMethod );
    }

    /** This method gives value of given shape function evaluated in
    given integration point of default integration method. It just
    call ShapeFunctionValue(IndexType IntegrationPointIndex,
    IndexType ShapeFunctionIndex, enum IntegrationMethod
    ThisMethod) with default integration method. There is no
    calculation and it just give it from shape functions values
    container if they are existing. Otherwise it gives you error
    which this value is not exist.

    @param IntegrationPointIndex index of integration point
    which shape functions evaluated in it.

    @param ShapeFunctionIndex index of node which correspounding
    shape function evaluated in given integration point.

    @return Value of given shape function in given integration
    point of default integration method.

    @see ShapeFunctionsValues
    @see ShapeFunctionsLocalGradients
    @see ShapeFunctionLocalGradient
    */
    double ShapeFunctionValue( IndexType IntegrationPointIndex, IndexType ShapeFunctionIndex ) const
    {
        return mpGeometryData->ShapeFunctionValue( IntegrationPointIndex, ShapeFunctionIndex );
    }

    /** This method gives value of given shape function evaluated in given
    integration point of given integration method. There is no
    calculation and it just give it from shape functions values
    container if they are existing. Otherwise it gives you error
    which this value is not exist.

    @param IntegrationPointIndex index of integration point
    which shape functions evaluated in it.

    @param ShapeFunctionIndex index of node which correspounding
    shape function evaluated in given integration point.

    @param ThisMethod integration method which shape function
    evaluated in its integration point.

    @return Value of given shape function in given integration
    point of given integration method.

    @see ShapeFunctionsValues
    @see ShapeFunctionsLocalGradients
    @see ShapeFunctionLocalGradient
    */
    double ShapeFunctionValue( IndexType IntegrationPointIndex, IndexType ShapeFunctionIndex, IntegrationMethod ThisMethod ) const
    {
        return mpGeometryData->ShapeFunctionValue( IntegrationPointIndex, ShapeFunctionIndex, ThisMethod );
    }

    /** This method gives value of given shape function evaluated in given
    point.

    @param rPoint Point of evaluation of the shape
    function. This point must be in local coordinate.

    @param ShapeFunctionIndex index of node which correspounding
    shape function evaluated in given integration point.

    @return Value of given shape function in given point.

    @see ShapeFunctionsValues
    @see ShapeFunctionsLocalGradients
    @see ShapeFunctionLocalGradient
    */
    virtual double ShapeFunctionValue( IndexType ShapeFunctionIndex, const CoordinatesArrayType& rCoordinates ) const
    {
        KRATOS_ERROR << "Calling base class ShapeFunctionValue method instead of derived class one. Please check the definition of derived class. " << *this << std::endl;

        return 0;
    }

    /** This method gives all shape functions gradients evaluated in all
    integration points of default integration method. It just
    call ShapeFunctionsLocalGradients(enum IntegrationMethod ThisMethod)
    with default integration method. There is no calculation and
    it just give it from shape functions values container.

    \note There is no control if there is any gradient calculated or not!

    @return shape functions' gradients \f$ F_{ijk} \f$ where i
    is the integration point index and j is the shape function
    index and k is local coordinate index. In other word
    component \f$ f_{ijk} \f$ is k'th component of gradient of
    the shape function corresponding to node j evaluated in
    integration point i of default integration method.

    @see ShapeFunctionsValues
    @see ShapeFunctionValue
    @see ShapeFunctionLocalGradient
    */

    const ShapeFunctionsGradientsType& ShapeFunctionsLocalGradients() const
    {
        return mpGeometryData->ShapeFunctionsLocalGradients();
    }

    /** This method gives all shape functions gradients evaluated in
    all integration points of given integration method. There is
    no calculation and it just give it from shape functions
    values container.

    \note There is no control if there is any gradient calculated or not!

    @param ThisMethod integration method which shape functions
    gradients evaluated in its integration points.

    @return shape functions' gradients \f$ F_{ijk} \f$ where i
    is the integration point index and j is the shape function
    index and k is local coordinate index. In other word
    component \f$ f_{ijk} \f$ is k'th component of gradient of
    the shape function corresponding to node j evaluated in
    integration point i of given integration method.

    @see ShapeFunctionsValues
    @see ShapeFunctionValue
    @see ShapeFunctionLocalGradient
    */
    const ShapeFunctionsGradientsType& ShapeFunctionsLocalGradients( IntegrationMethod ThisMethod ) const
    {
        return mpGeometryData->ShapeFunctionsLocalGradients( ThisMethod );
    }

    /** This method gives gradient of given shape function evaluated in
    given integration point of default integration method. It just
    call ShapeFunctionLocalGradient(IndexType IntegrationPointIndex,
    IndexType ShapeFunctionIndex, enum IntegrationMethod
    ThisMethod) with default integration method. There is no
    calculation and it just give it from shape functions values
    container if they are existing. Otherwise it gives you error
    which this value is not exist.

    @param IntegrationPointIndex index of integration point
    which shape function gradient evaluated in it.

    @param ShapeFunctionIndex index of node which correspounding
    shape function gradient evaluated in given integration point.

    @return Gradient of given shape function in given integration
    point of default integration method.

    @see ShapeFunctionsValues
    @see ShapeFunctionValue
    @see ShapeFunctionsLocalGradients
    */
    const Matrix& ShapeFunctionLocalGradient( IndexType IntegrationPointIndex )  const
    {
        return mpGeometryData->ShapeFunctionLocalGradient( IntegrationPointIndex );
    }

    /** This method gives gradient of given shape function evaluated
    in given integration point of given integration
    method. There is no calculation and it just give it from
    shape functions values container if they are
    existing. Otherwise it gives you error which this value is
    not exist.

    @param IntegrationPointIndex index of integration point
    which shape function gradient evaluated in it.

    @param ShapeFunctionIndex index of node which correspounding
    shape function gradient evaluated in given integration point.

    @param ThisMethod integration method which shape function gradient
    evaluated in its integration points.

    @return Gradient of given shape function in given integration
    point of given integration method.

    @see ShapeFunctionsValues
    @see ShapeFunctionValue
    @see ShapeFunctionsLocalGradients
    */
    const Matrix& ShapeFunctionLocalGradient(IndexType IntegrationPointIndex , IntegrationMethod ThisMethod)  const
    {
        return mpGeometryData->ShapeFunctionLocalGradient(IntegrationPointIndex, ThisMethod);
    }

    const Matrix& ShapeFunctionLocalGradient(IndexType IntegrationPointIndex, IndexType ShapeFunctionIndex, IntegrationMethod ThisMethod)  const
    {
        return mpGeometryData->ShapeFunctionLocalGradient(IntegrationPointIndex, ShapeFunctionIndex, ThisMethod);
    }


    /** This method gives gradient of all shape functions evaluated
     * in given point.
     * There is no calculation and it just give it from
     * shape functions values container if they are
     * existing. Otherwise it gives you error which this value is
     * not exist.
     *
     * @param rResult the given Container that will be overwritten by the solution
     * @param rPoint the given local coordinates the gradients will be evaluated for
     * @return a matrix of gradients for each shape function
     */
    virtual Matrix& ShapeFunctionsLocalGradients( Matrix& rResult, const CoordinatesArrayType& rPoint ) const
    {
        KRATOS_ERROR << "Calling base class ShapeFunctionsLocalGradients method instead of derived class one. Please check the definition of derived class. " << *this << std::endl;
        return rResult;
    }

    /** This method gives second order derivatives of all shape
     * functions evaluated in given point.
     *
     * @param rResult the given container will be overwritten by the results
     * @param rPoint the given local coordinates the derivatives will be evaluated for.
     * @return a third order tensor containing the second order derivatives for each shape function
     */
    virtual ShapeFunctionsSecondDerivativesType& ShapeFunctionsSecondDerivatives( ShapeFunctionsSecondDerivativesType& rResult, const CoordinatesArrayType& rPoint ) const
    {
        KRATOS_ERROR << "Calling base class ShapeFunctionsSecondDerivatives method instead of derived class one. Please check the definition of derived class. " << *this << std::endl;
        return rResult;
    }

    /** This method gives third order derivatives of all shape
     * functions evaluated in given point.
     *
     * @param rResult the given container will be overwritten by the results
     * @param rPoint the given local coordinates the derivatives will be evaluated for.
     * @return a fourth order tensor containing the second order derivatives for each shape function
     */
    virtual ShapeFunctionsThirdDerivativesType& ShapeFunctionsThirdDerivatives( ShapeFunctionsThirdDerivativesType& rResult, const CoordinatesArrayType& rPoint ) const
    {
        KRATOS_ERROR << "Calling base class ShapeFunctionsThirdDerivatives method instead of derived class one. Please check the definition of derived class. " << *this << std::endl;
        return rResult;
    }


    ShapeFunctionsGradientsType& ShapeFunctionsIntegrationPointsGradients( ShapeFunctionsGradientsType& rResult ) const
    {
        ShapeFunctionsIntegrationPointsGradients( rResult, mpGeometryData->DefaultIntegrationMethod() );
        return rResult;
    }

    virtual ShapeFunctionsGradientsType& ShapeFunctionsIntegrationPointsGradients(
        ShapeFunctionsGradientsType& rResult,
        IntegrationMethod ThisMethod ) const
    {
        const unsigned int integration_points_number = this->IntegrationPointsNumber( ThisMethod );

        if ( integration_points_number == 0 )
            KRATOS_ERROR << "This integration method is not supported" << *this << std::endl;

        if ( rResult.size() != integration_points_number )
            rResult.resize(  this->IntegrationPointsNumber( ThisMethod ), false  );

        //calculating the local gradients
        const ShapeFunctionsGradientsType& DN_De = ShapeFunctionsLocalGradients( ThisMethod );

        //loop over all integration points
        Matrix J(this->WorkingSpaceDimension(),this->LocalSpaceDimension()),Jinv(this->WorkingSpaceDimension(),this->LocalSpaceDimension());
        double DetJ;
        for ( unsigned int pnt = 0; pnt < integration_points_number; pnt++ )
        {
            if(rResult[pnt].size1() != this->WorkingSpaceDimension() ||  rResult[pnt].size2() != this->LocalSpaceDimension())
                rResult[pnt].resize( (*this).size(), this->LocalSpaceDimension(), false );
            this->Jacobian(J,pnt, ThisMethod);
            MathUtils<double>::GeneralizedInvertMatrix( J, Jinv, DetJ );
            noalias(rResult[pnt]) =  prod( DN_De[pnt], Jinv );
        }

        return rResult;
    }

    virtual ShapeFunctionsGradientsType& ShapeFunctionsIntegrationPointsGradients( ShapeFunctionsGradientsType& rResult, Vector& determinants_of_jacobian, IntegrationMethod ThisMethod ) const
    {
        const unsigned int integration_points_number = this->IntegrationPointsNumber( ThisMethod );

        if ( integration_points_number == 0 )
            KRATOS_ERROR << "This integration method is not supported " << *this << std::endl;

        if ( rResult.size() != integration_points_number )
            rResult.resize(  this->IntegrationPointsNumber( ThisMethod ), false  );
        if ( determinants_of_jacobian.size() != integration_points_number )
            determinants_of_jacobian.resize(  this->IntegrationPointsNumber( ThisMethod ), false  );

        //calculating the local gradients
        const ShapeFunctionsGradientsType& DN_De = ShapeFunctionsLocalGradients( ThisMethod );

        //loop over all integration points
        Matrix J(this->WorkingSpaceDimension(),this->LocalSpaceDimension());
        Matrix Jinv(this->WorkingSpaceDimension(),this->LocalSpaceDimension());
        double DetJ;
        for ( unsigned int pnt = 0; pnt < integration_points_number; pnt++ )
        {
            if(rResult[pnt].size1() != this->WorkingSpaceDimension() ||  rResult[pnt].size2() != this->LocalSpaceDimension())
                rResult[pnt].resize( (*this).size(), this->LocalSpaceDimension(), false );
            this->Jacobian(J,pnt, ThisMethod);
            MathUtils<double>::GeneralizedInvertMatrix( J, Jinv, DetJ );
            noalias(rResult[pnt]) =  prod( DN_De[pnt], Jinv );
            determinants_of_jacobian[pnt] = DetJ;
        }

        return rResult;
    }

    virtual ShapeFunctionsGradientsType& ShapeFunctionsIntegrationPointsGradients( ShapeFunctionsGradientsType& rResult, Vector& determinants_of_jacobian, IntegrationMethod ThisMethod, Matrix& ShapeFunctionsIntegrationPointsValues ) const
    {

        ShapeFunctionsIntegrationPointsGradients(rResult,determinants_of_jacobian,ThisMethod);
        ShapeFunctionsIntegrationPointsValues = ShapeFunctionsValues(ThisMethod);
        return rResult;
    }

    ///@}
    ///@name Input and output
    ///@{

    /** Returns geometry information as a string.
     * Returns geometry information as a string.
     *
     * @return String contains information about this geometry.
     *
     * @see Name()
     */
    std::string Info() const override {
      std::stringstream buffer;
      buffer << Dimension() << " dimensional geometry in " << WorkingSpaceDimension() << "D space";

      return buffer.str();
    }

    /** Returns the name of the geometry as a string.
     * Returns the name of the geometry as a string.
     *
     * Note: compiler's RVO should optimize this code automatically.
     *
     * @return String with the name of the geometry.
     *
     * @see Info()
     */
    virtual std::string Name() const {
      std::string geometryName = "BaseGeometry";
      KRATOS_ERROR << "Base geometry does not have a name." << std::endl;
      return geometryName;
    }

    /** Prints information about this object.
     * Prints information about this object.
     *
     * @param rOStream Output Stream.
     *
     * @see PrintName()
     * @see PrintData()
     */
    void PrintInfo(std::ostream& rOStream) const override {
      rOStream << Dimension()  << " dimensional geometry in " << WorkingSpaceDimension() << "D space";
    }

    /** Prints the name of the geometry.
     * Prints the name of the geometry.
     *
     * @param rOStream Output Stream.
     *
     * @see PrintInfo()
     * @see PrintData()
     */
    virtual void PrintName(std::ostream& rOstream) const {
      rOstream << Name() << std::endl;
    }

    /** Print geometry's data into given stream.
     * Prints it's points by the order they stored in the
     * geometry and then center point of geometry.
     *
     * @param rOStream Output Stream.
     *
     * @see PrintInfo()
     * @see PrintName()
     */
    void PrintData( std::ostream& rOStream ) const override {
      if(mpGeometryData) {
        mpGeometryData->PrintData( rOStream );
      }

      rOStream << std::endl;
      rOStream << std::endl;

      for (unsigned int i = 0; i < this->size(); ++i) {
        rOStream << "\tPoint " << i + 1 << "\t : ";
        (*this)[i].PrintData(rOStream);
        rOStream << std::endl;
      }

      rOStream << "\tCenter\t : ";

      Center().PrintData( rOStream );

      rOStream << std::endl;
      rOStream << std::endl;
      // rOStream << "\tLength\t : " << Length() << std::endl;
      // rOStream << "\tArea\t : " << Area() << std::endl;

      // Charlie: Volume is not defined by every geometry (2D geometries),
      // which can cause this call to generate a KRATOS_ERROR while trying
      // to call the base class Volume() method.

      // rOStream << "\tVolume\t : " << Volume() << std::endl;

      // Charlie: Can this be deleted?

      // for(unsigned int i = 0 ; i < mPoints.size() ; ++i) {
      //   rOStream << "    Point " << i+1 << "\t            : ";
      //   mPoints[i].PrintData(rOStream);
      //   rOStream << std::endl;
      // }
      //
      // rOStream << "    Center\t            : ";
      // Center().PrintData(rOStream);
      // rOStream << std::endl;
      // rOStream << std::endl;
      // rOStream << "    Length                  : " << Length() << std::endl;
      // rOStream << "    Area                    : " << Area() << std::endl;
      // rOStream << "    Volume                  : " << Volume();
    }


    ///@}
    ///@name Friends
    ///@{


    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{


    ///@}
    ///@name Protected member Variables
    ///@{


    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{

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
    virtual double InradiusToCircumradiusQuality() const {
      KRATOS_ERROR << "Calling base class 'InradiusToCircumradiusQuality' method instead of derived class one. Please check the definition of derived class. " << *this << std::endl;
      return 0.0;
    }

    /** Calculates the minimum to maximum edge length quality metric.
     * Calculates the minimum to maximum edge length quality metric.
     * This metric is bounded by the interval (0,1) being:
     *  1 -> Optimal value
     *  0 -> Worst value
     *
     * @formulae $$ \frac{h_{min}}{h_{max}} $$
     *
     * @return The Inradius to Circumradius Quality metric.
     */
    virtual double AreaToEdgeLengthRatio() const {
      KRATOS_ERROR << "Calling base class 'AreaToEdgeLengthRatio' method instead of derived class one. Please check the definition of derived class. " << *this << std::endl;
      return 0.0;
    }

    /** Calculates the shortest altitude to edge length quality metric.
     * Calculates the shortest altitude to edge length quality metric.
     * This metric is bounded by the interval (0,1) being:
     *  1 -> Optimal value
     *  0 -> Worst value
     *
     * @formulae $$ \frac{h_{min}}{h_{max}} $$
     *
     * @return The shortest altitude to edge length quality metric.
     */
    virtual double ShortestAltitudeToEdgeLengthRatio() const {
      KRATOS_ERROR << "Calling base class 'ShortestAltitudeToEdgeLengthRatio' method instead of derived class one. Please check the definition of derived class. " << *this << std::endl;
      return 0.0;
    }

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
    virtual double InradiusToLongestEdgeQuality() const {
      KRATOS_ERROR << "Calling base class 'InradiusToLongestEdgeQuality' method instead of derived class one. Please check the definition of derived class. " << *this << std::endl;
      return 0.0;
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
    virtual double ShortestToLongestEdgeQuality() const {
      KRATOS_ERROR << "Calling base class 'ShortestToLongestEdgeQuality' method instead of derived class one. Please check the definition of derived class. " << *this << std::endl;
      return 0.0;
    }

    /** Calculates the Regularity quality metric.
     * Calculates the Regularity quality metric.
     * This metric is bounded by the interval (-1,1) being:
     *  1 -> Optimal value
     *  0 -> Worst value
     *  -1 -> Negative volume
     *
     * \f$ \frac{4r}{H} \f$
     *
     * @return regularity quality.
     */
    virtual double RegularityQuality() const {
      KRATOS_ERROR << "Calling base class 'RegularityQuality' method instead of derived class one. Please check the definition of derived class. " << *this << std::endl;
      return 0.0;
    }

    /** Calculates the volume to surface area quality metric.
     * Calculates the volume to surface area quality metric.
     * This metric is bounded by the interval (-1,1) being:
     *   1 -> Optimal value
     *   0 -> Worst value
     *  -1 -> Negative volume
     *
     * \f$ \frac{V^4}{(\sum{A_{i}^{2}})^{3}} \f$
     *
     * @return volume to surface quality.
     */
    virtual double VolumeToSurfaceAreaQuality() const {
      KRATOS_ERROR << "Calling base class 'VolumeToSurfaceAreaQuality' method instead of derived class one. Please check the definition of derived class. " << *this << std::endl;
      return 0.0;
    }

    /** Calculates the Volume to edge length quaility metric.
     * Calculates the Volume to edge length quaility metric.
     * This metric is bounded by the interval (-1,1) being:
     *  1 -> Optimal value
     *  0 -> Worst value
     *  -1 -> Negative volume
     *
     * \f$ \frac{V^{2/3}}{\sum{l_{i}^{2}}} \f$
     *
     * @return Volume to edge length quality.
     */
    virtual double VolumeToEdgeLengthQuality() const {
      KRATOS_ERROR << "Calling base class 'VolumeToEdgeLengthQuality' method instead of derived class one. Please check the definition of derived class. " << *this << std::endl;
      return 0.0;
    }

    /** Calculates the volume to average edge lenght quality metric.
     * Calculates the volume to average edge lenght quality metric.
     * This metric is bounded by the interval (-1,1) being:
     *  1 -> Optimal value
     *  0 -> Worst value
     *  -1 -> Negative volume
     *
     * \f$ \frac{V}{\frac{1}{6}\sum{l_i}} \f$
     *
     * @return [description]
     */
    virtual double VolumeToAverageEdgeLength() const {
      KRATOS_ERROR << "Calling base class 'VolumeToAverageEdgeLength' method instead of derived class one. Please check the definition of derived class. " << *this << std::endl;
      return 0.0;
    }

    /** Calculates the volume to average edge length quality metric.
     * Calculates the volume to average edge length quality metric.
     * The average edge lenght is calculated using the RMS.
     * This metric is bounded by the interval (-1,1) being:
     *  1 -> Optimal value
     *  0 -> Worst value
     *  -1 -> Negative volume
     *
     * \f$ \frac{V}{\sqrt{\frac{1}{6}\sum{A_{i}^{2}}}} \f$
     *
     * @return [description]
     */
    virtual double VolumeToRMSEdgeLength() const {
      KRATOS_ERROR << "Calling base class 'VolumeToRMSEdgeLength' method instead of derived class one. Please check the definition of derived class. " << *this << std::endl;
      return 0.0;
    }

    ///@}
    ///@name Protected  Access
    ///@{


    ///@}
    ///@name Protected Inquiry
    ///@{


    ///@}
    ///@name Protected LifeCycle
    ///@{

    /** Protected Constructor.
    Avoids object to be created Except for derived classes
    */


    ///@}



private:
    ///@name Static Member Variables
    ///@{

    // static const GeometryData msEmptyGeometryData;

    ///@}
    ///@name Member Variables
    ///@{

    GeometryData const* mpGeometryData;


    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save( Serializer& rSerializer ) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BaseType );
//                 rSerializer.save( "Geometry Data", mpGeometryData );
    }

    void load( Serializer& rSerializer ) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BaseType );
        //rSerializer.load( "Geometry Data", const_cast<GeometryData*>( mpGeometryData ) );
    }


    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    static const GeometryData& GeometryDataInstance()
    {
        IntegrationPointsContainerType integration_points = {};
        ShapeFunctionsValuesContainerType shape_functions_values = {};
        ShapeFunctionsLocalGradientsContainerType shape_functions_local_gradients = {};
        static GeometryData s_geometry_data(3,
                            3,
                            3,
                            GeometryData::GI_GAUSS_1,
                            integration_points,
                            shape_functions_values,
                            shape_functions_local_gradients);

        return s_geometry_data;
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

    template<class TOtherPointType> friend class Geometry;

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


/// input stream function
template<class TPointType>
inline std::istream& operator >> ( std::istream& rIStream,
                                   Geometry<TPointType>& rThis );

/// output stream function
template<class TPointType>
inline std::ostream& operator << ( std::ostream& rOStream,
                                   const Geometry<TPointType>& rThis )
{
    rThis.PrintInfo( rOStream );
    rOStream << std::endl;
    rThis.PrintData( rOStream );

    return rOStream;
}

///@}


}  // namespace Kratos.

#endif // KRATOS_GEOMETRY_H_INCLUDED  defined
