//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//                   Riccardo Rossi
//                   Janosch Stascheit
//                   Felix Nagel
//  contributors:    Hoang Giang Bui
//                   Josep Maria Carbonell
//                   Carlos Roig
//                   Vicente Mataix Ferrandiz
//

#pragma once

// System includes
#include <typeinfo>

// External includes

// Project includes
#include "geometries/geometry_data.h"
#include "geometries/point.h"
#include "containers/pointer_vector.h"
#include "containers/data_value_container.h"
#include "utilities/math_utils.h"
#include "input_output/logger.h"
#include "integration/integration_info.h"

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
class Geometry
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
      VOLUME_TO_RMS_EDGE_LENGTH,
      MIN_DIHEDRAL_ANGLE,
      MAX_DIHEDRAL_ANGLE,
      MIN_SOLID_ANGLE
    };

    /**
     * @brief This defines the different methods to compute the lumping methods
     * @details The three methods available are:
     *      - The row sum method
     *      - Diagonal scaling
     *      - Evaluation of M using a quadrature involving only the nodal points and thus automatically yielding a diagonal matrix for standard element shape function
     */
    enum class LumpingMethods {
        ROW_SUM,
        DIAGONAL_SCALING,
        QUADRATURE_ON_NODES
    };

    /** Array of counted pointers to point. This type used to hold
    geometry's points.
    */
    typedef PointerVector<TPointType> PointsArrayType;

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
    typedef std::array<IntegrationPointsArrayType, static_cast<int>(GeometryData::IntegrationMethod::NumberOfIntegrationMethods)> IntegrationPointsContainerType;

    /** A third order tensor used as shape functions' values
    continer.
    */
    typedef std::array<Matrix, static_cast<int>(GeometryData::IntegrationMethod::NumberOfIntegrationMethods)> ShapeFunctionsValuesContainerType;

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

    /// data type stores in this container.
    typedef typename PointType::Pointer PointPointerType;
    typedef const PointPointerType ConstPointPointerType;
    typedef TPointType& PointReferenceType;
    typedef const TPointType& ConstPointReferenceType;
    typedef std::vector<PointPointerType> PointPointerContainerType;

    /// PointsArrayType typedefs
    typedef typename PointsArrayType::iterator iterator;
    typedef typename PointsArrayType::const_iterator const_iterator;

    typedef typename PointsArrayType::ptr_iterator ptr_iterator;
    typedef typename PointsArrayType::ptr_const_iterator ptr_const_iterator;
    typedef typename PointsArrayType::difference_type difference_type;

    static constexpr IndexType BACKGROUND_GEOMETRY_INDEX = std::numeric_limits<IndexType>::max();

    ///@}
    ///@name Life Cycle
    ///@{

    /// Standard Constructor. Generates self assigned id.
    Geometry()
        : mId(GenerateSelfAssignedId())
        , mpGeometryData(&GeometryDataInstance())
    {
    }

    /// Standard Constructor with a geometry Id
    Geometry(IndexType GeometryId)
        : mpGeometryData(&GeometryDataInstance())
    {
        SetId(GeometryId);
    }

    /// Standard Constructor with a Name
    Geometry(const std::string& GeometryName)
        : mId(GenerateId(GeometryName))
        , mpGeometryData(&GeometryDataInstance())
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
    Geometry(
        const PointsArrayType &ThisPoints,
        GeometryData const *pThisGeometryData = &GeometryDataInstance())
        : mId(GenerateSelfAssignedId())
        , mpGeometryData(pThisGeometryData)
        , mPoints(ThisPoints)
    {
    }

    Geometry(
        IndexType GeometryId,
        const PointsArrayType& ThisPoints,
        GeometryData const* pThisGeometryData = &GeometryDataInstance())
        : mpGeometryData(pThisGeometryData)
        , mPoints(ThisPoints)
    {
        SetId(GeometryId);
    }

    Geometry(
        const std::string& GeometryName,
        const PointsArrayType& ThisPoints,
        GeometryData const* pThisGeometryData = &GeometryDataInstance())
        : mId(GenerateId(GeometryName))
        , mpGeometryData(pThisGeometryData)
        , mPoints(ThisPoints)
    {
    }

    /**
    * @brief Copy constructor
    *
    * @note Does not copy the points but shares same points with
    *       the original geometry. Any change to the points of the
    *       copied geometry affect point of original geometry, too.
    * @note Copied geometry shares the same Id as the
    *       original geometry.
    */
    Geometry( const Geometry& rOther )
        : mId(rOther.mId),
          mpGeometryData(rOther.mpGeometryData),
          mPoints(rOther.mPoints),
          mData(rOther.mData)
    {
    }

    /**
    * @brief Copy constructor with TOtherPointType
    *
    *        Copies geometry with a different type of points.
    *        TOtherPointType* must be implicity convertible
    *        to TPointType of the original geometry.
    *
    * @note Does not copy the points but shares same points with
    *       the original geometry. Any change to the points of the
    *       copied geometry affect point of original geometry, too.
    * @note Copied geometry shares the same Id as the
    *       original geometry.
    */
    template<class TOtherPointType>
    Geometry( Geometry<TOtherPointType> const & rOther )
        : mId(rOther.mId),
          mpGeometryData(rOther.mpGeometryData),
          mData(rOther.mData)
    {
        mPoints = new PointsArrayType(rOther.begin(), rOther.end());
    }

    /// Destructor. Do nothing!!!
    virtual ~Geometry() {}

    virtual GeometryData::KratosGeometryFamily GetGeometryFamily() const
    {
        return GeometryData::KratosGeometryFamily::Kratos_generic_family;
    }

    virtual GeometryData::KratosGeometryType GetGeometryType() const
    {
        return GeometryData::KratosGeometryType::Kratos_generic_type;
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
        mpGeometryData = rOther.mpGeometryData;
        mPoints = rOther.mPoints;
        mData = rOther.mData;

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

     operator PointsArrayType&()
    {
        return mPoints;
    }

    ///@}
    ///@name PointerVector Operators
    ///@{

    TPointType& operator[](const SizeType& i)
    {
        return mPoints[i];
    }

    TPointType const& operator[](const SizeType& i) const
    {
        return mPoints[i];
    }

    PointPointerType& operator()(const SizeType& i)
    {
        return mPoints(i);
    }

    ConstPointPointerType& operator()(const SizeType& i) const
    {
        return mPoints(i);
    }

    ///@}
    ///@name PointerVector Operations
    ///@{

    iterator                   begin()
    {
        return iterator(mPoints.begin());
    }
    const_iterator             begin() const
    {
        return const_iterator(mPoints.begin());
    }
    iterator                   end()
    {
        return iterator(mPoints.end());
    }
    const_iterator             end() const
    {
        return const_iterator(mPoints.end());
    }
    ptr_iterator               ptr_begin()
    {
        return mPoints.ptr_begin();
    }
    ptr_const_iterator         ptr_begin() const
    {
        return mPoints.ptr_begin();
    }
    ptr_iterator               ptr_end()
    {
        return mPoints.ptr_end();
    }
    ptr_const_iterator         ptr_end() const
    {
        return mPoints.ptr_end();
    }
    PointReferenceType        front()       /* nothrow */
    {
        assert(!empty());
        return mPoints.front();
    }
    ConstPointReferenceType  front() const /* nothrow */
    {
        assert(!empty());
        return mPoints.front();
    }
    PointReferenceType        back()        /* nothrow */
    {
        assert(!empty());
        return mPoints.back();
    }
    ConstPointReferenceType  back() const  /* nothrow */
    {
        assert(!empty());
        return mPoints.back();
    }

    SizeType size() const
    {
        return mPoints.size();
    }

    /**
    * @detail Returns the number of the points/ nodes
    *         belonging to this geometry.
    * @return Number of points/ nodes.
    */
    SizeType PointsNumber() const {
        return this->size();
    }

    /// Returns number of points per direction.
    virtual SizeType PointsNumberInDirection(IndexType LocalDirectionIndex) const
    {
        KRATOS_ERROR << "Trying to access PointsNumberInDirection from geometry base class." << std::endl;
    }

    SizeType max_size() const
    {
        return mPoints.max_size();
    }

    void swap(GeometryType& rOther)
    {
        mPoints.swap(rOther.mPoints);
    }

    void push_back(PointPointerType x)
    {
        mPoints.push_back(x);
    }

    void clear()
    {
        mPoints.clear();
    }

    void reserve(int dim)
    {
        mPoints.reserve(dim);
    }

    int capacity()
    {
        return mPoints.capacity();
    }

    /////@}
    /////@name Access
    /////@{

    ///** Gives a reference to underly normal container. */
    PointPointerContainerType& GetContainer()
    {
        return mPoints.GetContainer();
    }

    /** Gives a constant reference to underly normal container. */
    const PointPointerContainerType& GetContainer() const
    {
        return mPoints.GetContainer();
    }

    ///@}
    ///@name Data Container
    ///@{

    /**
     * Access Data:
     */
    DataValueContainer& GetData()
    {
      return mData;
    }

    DataValueContainer const& GetData() const
    {
      return mData;
    }

    void SetData(DataValueContainer const& rThisData)
    {
      mData = rThisData;
    }

    /**
     * Check if the Data exists with Has(..) methods:
     */
    template<class TDataType> bool Has(const Variable<TDataType>& rThisVariable) const
    {
        return mData.Has(rThisVariable);
    }

    /**
     * Set Data with SetValue and the Variable to set:
     */
    template<class TVariableType> void SetValue(
        const TVariableType& rThisVariable,
        typename TVariableType::Type const& rValue)
    {
        mData.SetValue(rThisVariable, rValue);
    }

    /**
     * Get Data with GetValue and the Variable to get:
     */
    template<class TVariableType> typename TVariableType::Type& GetValue(
        const TVariableType& rThisVariable)
    {
        return mData.GetValue(rThisVariable);
    }

    template<class TVariableType> typename TVariableType::Type const& GetValue(
        const TVariableType& rThisVariable) const
    {
        return mData.GetValue(rThisVariable);
    }

    ///@}
    ///@name Dynamic access to internals
    ///@{

    /* Assigns a value to the geometry,
     * according to a variable.
     * Allows dynamic interfaces with each respective geometry.
     */

    /// Assign with bool
    virtual void Assign(
        const Variable<bool>& rVariable,
        const bool Input) {}

    /// Assign with int
    virtual void Assign(
        const Variable<int>& rVariable,
        const int Input) {}

    /// Assign with double
    virtual void Assign(
        const Variable<double>& rVariable,
        const double Input) {}

    /// Assign with array_1d<double, 2>
    virtual void Assign(
        const Variable<array_1d<double, 2>>& rVariable,
        const array_1d<double, 2>& rInput) {}

    /// Assign with array_1d<double, 3>
    virtual void Assign(
        const Variable<array_1d<double, 3>>& rVariable,
        const array_1d<double, 3>& rInput) {}

    /// Assign with array_1d<double, 6>
    virtual void Assign(
        const Variable<array_1d<double, 6>>& rVariable,
        const array_1d<double, 6>& rInput) {}

    /// Assign with Vector
    virtual void Assign(
        const Variable<Vector>& rVariable,
        const Vector& rInput) {}

    /// Assign with Matrix
    virtual void Assign(
        const Variable<Matrix>& rVariable,
        const Matrix& rInput) {}

    /* Calculate either provides, gets or calculates a certain value,
     * according to a variable.
     */

    /// Calculate with bool
    virtual void Calculate(
        const Variable<bool>& rVariable,
        bool& rOutput) const {}

    /// Calculate with int
    virtual void Calculate(
        const Variable<int>& rVariable,
        int& rOutput) const {}

    /// Calculate with double
    virtual void Calculate(
        const Variable<double>& rVariable,
        double& rOutput) const {}

    /// Calculate with array_1d<double, 2>
    virtual void Calculate(
        const Variable<array_1d<double, 2>>& rVariable,
        array_1d<double, 2>& rOutput) const {}

    /// Calculate with array_1d<double, 3>
    virtual void Calculate(
        const Variable<array_1d<double, 3>>& rVariable,
        array_1d<double, 3>& rOutput) const {}

    /// Calculate with array_1d<double, 6>
    virtual void Calculate(
        const Variable<array_1d<double, 6>>& rVariable,
        array_1d<double, 6>& rOutput) const {}

    /// Calculate with Vector
    virtual void Calculate(
        const Variable<Vector>& rVariable,
        Vector& rOutput) const {}

    /// Calculate with Matrix
    virtual void Calculate(
        const Variable<Matrix>& rVariable,
        Matrix& rOutput) const {}

    ///@}
    ///@name Inquiry
    ///@{

    /**
     * @brief Checks if two GeometryType have the same type
     * @return True if the objects are the same type, false otherwise
     */
    inline static bool HasSameType(
        const GeometryType& rLHS,
        const GeometryType& rRHS)
    {
        return (typeid(rLHS) == typeid(rRHS));
    }

    /**
     * @brief Checks if two GeometryType have the same type (pointer version)
     * @return True if the objects are the same type, false otherwise
     */
    inline static bool HasSameType(
        const GeometryType * rLHS,
        const GeometryType* rRHS)
    {
        return GeometryType::HasSameType(*rLHS, *rRHS);
    }

    /**
     * @brief Checks if two GeometryType have the same geometry type
     * @return True if the geometries are the same type, false otherwise
     */
    inline static bool HasSameGeometryType(const GeometryType& rLHS, const GeometryType& rRHS) {
        return (rLHS.GetGeometryType() == rRHS.GetGeometryType());
    }

    /**
     * @brief Checks if two GeometryType have the same geometry type (pointer version)
     * @return True if the geometries are the same type, false otherwise
     */
    inline static bool HasSameGeometryType(
        const GeometryType* rLHS,
        const GeometryType* rRHS)
    {
        return GeometryType::HasSameGeometryType(*rLHS, *rRHS);
    }

    /**
     * @brief Checks if two GeometryType are the same
     * @return True if the object is the same, false otherwise
     */
    inline static bool IsSame(
        const GeometryType& rLHS,
        const GeometryType& rRHS)
    {
        return GeometryType::HasSameType(rLHS, rRHS) && GeometryType::HasSameGeometryType(rLHS, rRHS);
    }

    /**
     * @brief Checks if two GeometryType are the same (pointer version)
     * @return True if the object is the same, false otherwise
     */
    inline static bool IsSame(
        const GeometryType* rLHS,
        const GeometryType* rRHS)
    {
        return GeometryType::HasSameType(*rLHS, *rRHS) && GeometryType::HasSameGeometryType(*rLHS, *rRHS);
    }

    bool empty() const
    {
        return mPoints.empty();
    }

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Creates a new geometry pointer
     * @param rThisPoints the nodes of the new geometry
     * @return Pointer to the new geometry
     */
    virtual Pointer Create(
        PointsArrayType const& rThisPoints
    ) const
    {
        // Create geometry
        auto p_geom = this->Create(0, rThisPoints);

        // Generate Id
        IndexType id = reinterpret_cast<IndexType>(p_geom.get());

        // Sets second bit to zero.
        p_geom->SetIdSelfAssigned(id);

        // Sets first bit to zero.
        p_geom->SetIdNotGeneratedFromString(id);

        // Sets Id
        p_geom->SetIdWithoutCheck(id);

        return p_geom;
    }

    /**
     * @brief Creates a new geometry pointer
     * @param NewGeometryId the ID of the new geometry
     * @param rThisPoints the nodes of the new geometry
     * @return Pointer to the new geometry
     */
    virtual Pointer Create(
        const IndexType NewGeometryId,
        PointsArrayType const& rThisPoints
    ) const
    {
        return Pointer( new Geometry( NewGeometryId, rThisPoints, mpGeometryData));
    }

    /**
     * @brief Creates a new geometry pointer
     * @param rNewGeometryName the name of the new geometry
     * @param rThisPoints the nodes of the new geometry
     * @return Pointer to the new geometry
     */
    Pointer Create(
        const std::string& rNewGeometryName,
        PointsArrayType const& rThisPoints
        ) const
    {
        auto p_geom = this->Create(0, rThisPoints);
        p_geom->SetId(rNewGeometryName);
        return p_geom;
    }

    /**
     * @brief Creates a new geometry pointer
     * @param rGeometry Reference to an existing geometry
     * @return Pointer to the new geometry
     */
    virtual Pointer Create(
        const GeometryType& rGeometry
    ) const
    {
        // Create geometry
        auto p_geom = this->Create(0, rGeometry);

        // Generate Id
        IndexType id = reinterpret_cast<IndexType>(p_geom.get());

        // Sets second bit to zero.
        p_geom->SetIdSelfAssigned(id);

        // Sets first bit to zero.
        p_geom->SetIdNotGeneratedFromString(id);

        // Sets Id
        p_geom->SetIdWithoutCheck(id);

        return p_geom;
    }

    /**
     * @brief Creates a new geometry pointer
     * @param NewGeometryId the ID of the new geometry
     * @param rGeometry Reference to an existing geometry
     * @return Pointer to the new geometry
     */
    virtual Pointer Create(
        const IndexType NewGeometryId,
        const GeometryType& rGeometry
    ) const
    {
        auto p_geometry = Pointer( new Geometry( NewGeometryId, rGeometry.Points(), mpGeometryData));
        p_geometry->SetData(rGeometry.GetData());
        return p_geometry;
    }

    /**
     * @brief Creates a new geometry pointer
     * @param rNewGeometryName the name of the new geometry
     * @param rGeometry Reference to an existing geometry
     * @return Pointer to the new geometry
     */
    Pointer Create(
        const std::string& rNewGeometryName,
        const GeometryType& rGeometry
        ) const
    {
        auto p_geom = this->Create(0, rGeometry);
        p_geom->SetId(rNewGeometryName);
        return p_geom;
    }

    /** This methods will create a duplicate of all its points and
    substitute them with its points. */
    void ClonePoints()
    {
        for ( ptr_iterator i = this->ptr_begin() ; i != this->ptr_end() ; i++ )
            *i = typename PointType::Pointer( new PointType( **i ) );
    }

    ///@}
    ///@name Geometry Data and Geometry Shape Function Container
    ///@{

    /**
    * @brief GeometryData contains all information about dimensions
    *        and has a set of precomputed values for integration points
    *        and shape functions, including derivatives.
    * @return the geometry data of a certain geometry class.
    */
    GeometryData const& GetGeometryData() const
    {
        return *mpGeometryData;
    }

    /* @brief SetGeometryShapeFunctionContainer updates the GeometryShapeFunctionContainer within
     *        the GeometryData. This function works only for geometries with a non-const GeometryData.
     *        E.g. QuadraturePointGeometries.
     */
    virtual void SetGeometryShapeFunctionContainer(
        const GeometryShapeFunctionContainer<GeometryData::IntegrationMethod>&  rGeometryShapeFunctionContainer)
    {
        KRATOS_ERROR <<
            "Calling SetGeometryShapeFunctionContainer from base geometry class."
            << std::endl;
    }

    ///@}
    ///@name Id
    ///@{

    /// Id of this Geometry
    IndexType const& Id() const
    {
        return mId;
    }

    /// Returns if id was generated from a geometry name
    bool IsIdGeneratedFromString()
    {
        return IsIdGeneratedFromString(mId);
    }

    /// Returns if id was generated by itself
    bool IsIdSelfAssigned()
    {
        return IsIdSelfAssigned(mId);
    }

    /// Sets Id of this Geometry
    void SetId(const IndexType Id)
    {
        // The first bit of the Id is used to detect if Id
        // is int or hash of name. Second bit defines if Id
        // is self assigned or not.
        KRATOS_ERROR_IF(IsIdGeneratedFromString(Id)
            || IsIdSelfAssigned(Id))
            << "Id: " << Id << " out of range. The Id must me lower than 2^62 = 4.61e+18. "
            << "Geometry being recognized as generated from string: " << IsIdGeneratedFromString(Id)
            << ", self assigned: " << IsIdSelfAssigned(Id) << "."
            << std::endl;

        mId = Id;
    }

    /// Sets Id with the use of the name of this geometry
    void SetId(const std::string& rName)
    {
        mId = GenerateId(rName);
    }

    /// Gets the corresponding hash-Id to a string name
    static inline IndexType GenerateId(const std::string& rName)
    {
        // Create id hash from provided name.
        std::hash<std::string> string_hash_generator;
        auto id = string_hash_generator(rName);

        // Sets first bit to one.
        SetIdGeneratedFromString(id);

        // Sets second bit to zero.
        SetIdNotSelfAssigned(id);

        return id;
    }

    ///@}
    ///@name Parent
    ///@{

    /**
    * @brief Some geometries require relations to other geometries. This is the
    *        case for e.g. quadrature points. To reach the parent geometry
    *        this function can be used.
    * @return Parent geometry of this geometry object.
    */
    virtual GeometryType& GetGeometryParent(IndexType Index) const
    {
        KRATOS_ERROR <<
            "Calling GetGeometryParent from base geometry class."
            << std::endl;
    }

    /**
    * @brief Some geometries require relations to other geometries. This is the
    *        case for e.g. quadrature points. To set or change the parent geometry
    *        this function can be used.
    * @param Parent geometry of this geometry object.
    */
    virtual void SetGeometryParent(GeometryType* pGeometryParent)
    {
        KRATOS_ERROR <<
            "Calling SetGeometryParent from base geometry class."
            << std::endl;
    }

    ///@}
    ///@name Geometry part functions
    ///@{

    /**
    * @brief Used for composite geometries. It returns the
    *        the geometry part, corresponding to the Index.
    * @param Index of the geometry part. This index can be used differently
    *        within the derived classes.
    * @return reference to corresponding geometry.
     */
    virtual GeometryType& GetGeometryPart(const IndexType Index)
    {
        return *pGetGeometryPart(Index);
    }

    /**
    * @brief Used for composite geometries. It returns the
    *        the geometry part, corresponding to the Index.
    * @param Index of the geometry part. This index can be used differently
    *        within the derived classes.
    * @return const reference to corresponding geometry.
    */
    virtual const GeometryType& GetGeometryPart(const IndexType Index) const
    {
        return *pGetGeometryPart(Index);
    }

    /**
    * @brief Used for composite geometries. It returns the pointer
    *        of a geometry part, corresponding to the Index.
    * @param Index of the geometry part. This index can be used differently
    *        within the derived classes.
    * @return pointer to corresponding geometry.
    */
    virtual typename GeometryType::Pointer pGetGeometryPart(const IndexType Index)
    {
        KRATOS_ERROR << "Calling base class 'pGetGeometryPart' method instead of derived function."
            << " Please check the definition in the derived class. " << *this << std::endl;
    }

    /**
    * @brief Used for composite geometries. It returns the const pointer
    *        of a geometry part, corresponding to the Index.
    * @details This index is dependent on the derived implementation.
    * @param Index of the geometry part. This index can be used differently
    *        within the derived classes.
    * @return const pointer to corresponding geometry.
    */
    virtual const typename GeometryType::Pointer pGetGeometryPart(const IndexType Index) const
    {
        KRATOS_ERROR << "Calling base class 'pGetGeometryPart' method instead of derived function."
            << " Please check the definition in the derived class. " << *this << std::endl;
    }

    /**
     * @brief Allows to exchange certain geometries.
     * @param Index of the geometry part. 0->Master; 1->Slave
     * @param pGeometry The new geometry to add
     */
    virtual void SetGeometryPart(
        const IndexType Index,
        GeometryType::Pointer pGeometry
        )
    {
        KRATOS_ERROR << "Calling base class 'SetGeometryPart' method instead of derived function."
            << " Please check the definition in the derived class. " << *this << std::endl;
    }

    /**
     * @brief Allows to enhance the coupling geometry, with another geometry.
     * @param pGeometry The new geometry to add
     */
    virtual IndexType AddGeometryPart(GeometryType::Pointer pGeometry)
    {
        KRATOS_ERROR << "Calling base class 'AddGeometryPart' method instead of derived function."
            << " Please check the definition in the derived class. " << *this << std::endl;
    }

    /**
     * @brief Removes a geometry part
     * @param pGeometry The new geometry to remove
     */
    virtual void RemoveGeometryPart(GeometryType::Pointer pGeometry)
    {
        KRATOS_ERROR << "Calling base class 'RemoveGeometryPart' method instead of derived function."
            << " Please check the definition in the derived class. " << *this << std::endl;
    }

    /**
     * @brief Removes a geometry part
     * @param Index of the geometry part.
     */
    virtual void RemoveGeometryPart(const IndexType Index)
    {
        KRATOS_ERROR << "Calling base class 'RemoveGeometryPart' method instead of derived function."
            << " Please check the definition in the derived class. " << *this << std::endl;
    }

    /**
    * @brief Use to check if certain Indexed object is
    *        within the geometry parts of this geometry.
    * @param Index of the geometry part. This index can be used differently
    *        within the derived classes.
    * @return true if has geometry part
    */
    virtual bool HasGeometryPart(const IndexType Index) const
    {
        KRATOS_ERROR << "Calling base class 'HasGeometryPart' method instead of derived function."
            << " Please check the definition in the derived class. " << *this << std::endl;
    }

    /**
    * @return the number of geometry parts that this geometry contains.
    */
    virtual SizeType NumberOfGeometryParts() const
    {
        return 0;
    }

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Lumping factors for the calculation of the lumped mass matrix
     * @param rResult Vector containing the lumping factors
     * @param LumpingMethod The lumping method considered. The three methods available are:
     *      - The row sum method
     *      - Diagonal scaling
     *      - Evaluation of M using a quadrature involving only the nodal points and thus automatically yielding a diagonal matrix for standard element shape function
     */
    virtual Vector& LumpingFactors(
        Vector& rResult,
        const LumpingMethods LumpingMethod = LumpingMethods::ROW_SUM
        )  const
    {
        const SizeType number_of_nodes = this->size();
        const SizeType local_space_dimension = this->LocalSpaceDimension();

        // Clear lumping factors
        if (rResult.size() != number_of_nodes)
            rResult.resize(number_of_nodes, false);
        noalias(rResult) = ZeroVector(number_of_nodes);

        if (LumpingMethod == LumpingMethods::ROW_SUM) {
            const IntegrationMethod integration_method = GetDefaultIntegrationMethod();
            const GeometryType::IntegrationPointsArrayType& r_integrations_points = this->IntegrationPoints( integration_method );
            const Matrix& r_Ncontainer = this->ShapeFunctionsValues(integration_method);

            // Vector fo jacobians
            Vector detJ_vector(r_integrations_points.size());
            DeterminantOfJacobian(detJ_vector, integration_method);

            // Iterate over the integration points
            double domain_size = 0.0;
            for ( IndexType point_number = 0; point_number < r_integrations_points.size(); ++point_number ) {
                const double integration_weight = r_integrations_points[point_number].Weight() * detJ_vector[point_number];
                const Vector& rN = row(r_Ncontainer,point_number);

                // Computing domain size
                domain_size += integration_weight;

                for ( IndexType i = 0; i < number_of_nodes; ++i ) {
                    rResult[i] += rN[i] * integration_weight;
                }
            }

            // Divide by the domain size
            for ( IndexType i = 0; i < number_of_nodes; ++i ) {
                rResult[i] /= domain_size;
            }
        } else if (LumpingMethod == LumpingMethods::DIAGONAL_SCALING) {
            IntegrationMethod integration_method = GetDefaultIntegrationMethod();
            int j = std::min(static_cast<int>(integration_method) + 1, 4);
            integration_method = static_cast<IntegrationMethod>(j);
            const GeometryType::IntegrationPointsArrayType& r_integrations_points = this->IntegrationPoints( integration_method );
            const Matrix& r_Ncontainer = this->ShapeFunctionsValues(integration_method);

            // Vector fo jacobians
            Vector detJ_vector(r_integrations_points.size());
            DeterminantOfJacobian(detJ_vector, integration_method);

            // Iterate over the integration points
            for ( IndexType point_number = 0; point_number < r_integrations_points.size(); ++point_number ) {
                const double detJ = detJ_vector[point_number];
                const double integration_weight = r_integrations_points[point_number].Weight() * detJ;
                const Vector& rN = row(r_Ncontainer,point_number);

                for ( IndexType i = 0; i < number_of_nodes; ++i ) {
                    rResult[i] += std::pow(rN[i], 2) * integration_weight;
                }
            }

            // Computing diagonal scaling coefficient
            double total_value = 0.0;
            for ( IndexType i = 0; i < number_of_nodes; ++i ) {
                total_value += rResult[i];
            }
            for ( IndexType i = 0; i < number_of_nodes; ++i ) {
                rResult[i] /= total_value;
            }
        } else if (LumpingMethod == LumpingMethods::QUADRATURE_ON_NODES) {
            // Divide by the domain size
            const double domain_size = DomainSize();

            // Getting local coordinates
            Matrix local_coordinates(number_of_nodes, local_space_dimension);
            PointsLocalCoordinates(local_coordinates);
            Point local_point(ZeroVector(3));
            array_1d<double, 3>& r_local_coordinates = local_point.Coordinates();

            // Iterate over integration points
            const GeometryType::IntegrationPointsArrayType& r_integrations_points = this->IntegrationPoints( GeometryData::IntegrationMethod::GI_GAUSS_1 ); // First order
            const double weight = r_integrations_points[0].Weight()/static_cast<double>(number_of_nodes);
            for ( IndexType point_number = 0; point_number < number_of_nodes; ++point_number ) {
                for ( IndexType dim = 0; dim < local_space_dimension; ++dim ) {
                    r_local_coordinates[dim] = local_coordinates(point_number, dim);
                }
                const double detJ = DeterminantOfJacobian(local_point);
                rResult[point_number] = weight * detJ/domain_size;
            }
        }

        return rResult;
    }

    ///@}
    ///@name Informations
    ///@{

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

    ///@}
    ///@name Mathematical Informations
    ///@{

    /// Return polynomial degree of the geometry in a certain direction
    virtual SizeType PolynomialDegree(IndexType LocalDirectionIndex) const
    {
        KRATOS_ERROR << "Trying to access PolynomialDegree from geometry base class." << std::endl;
    }

    ///@}
    ///@name Geometrical Informations
    ///@{

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

    /**
     * @brief This method calculate and return area or surface area of this geometry depending to it's dimension.
     * @details For one dimensional geometry it returns length, for two dimensional it gives area and for three dimensional geometries it gives surface area.
     * @return double value contains area or surface area.
     * @see Length()
     * @see Volume()
     * @see DomainSize()
     */
    virtual double Area() const {
        KRATOS_ERROR << "Calling base class 'Area' method instead of derived class one. Please check the definition of derived class. " << *this << std::endl;
        return 0.0;
    }

    /**
     * @brief This method calculate and return volume of this geometry.
     * @details For one and two dimensional geometry it returns zero and for three dimensional it gives volume of geometry.
     * @return double value contains volume.
     * @see Length()
     * @see Area()
     * @see DomainSize()
     */
    virtual double Volume() const {
        KRATOS_ERROR << "Calling base class 'Volume' method instead of derived class one. Please check the definition of derived class. " << *this << std::endl;
        return 0.0;
    }

    /**
     * @brief This method calculate and return length, area or volume of this geometry depending to it's dimension.
     * @details For one dimensional geometry it returns its length, for two dimensional it gives area and for three dimensional geometries it gives its volume.
     * @return double value contains length, area or volume.
     * @see Length()
     * @see Area()
     * @see Volume()
     */
    virtual double DomainSize() const {
        const SizeType local_dimension = this->LocalSpaceDimension();
        if (local_dimension == 1) { // 1D geometry
            return this->Length();
        } else if (local_dimension == 2) { // 2D geometry
            return this->Area();
        } else { // 3D geometry
            return this->Volume();
        }
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
    virtual bool HasIntersection(const GeometryType& ThisGeometry) const {
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
    virtual bool HasIntersection(const Point& rLowPoint, const Point& rHighPoint) const {
      KRATOS_ERROR << "Calling base class 'HasIntersection' method instead of derived class one. Please check the definition of derived class. " << *this << std::endl;
      return false;
    }

    // virtual void BoundingBox(BoundingBoxType& rResult) const
    // {
    //
    //   Bounding_Box(rResult.LowPoint(), rResult.HighPoint());
    // }

    /**
     * @brief Calculates the boundingbox of the geometry.
     * @details Corresponds with the highest and lowest point in space
     * @param rLowPoint  Lower point of the boundingbox.
     * @param rHighPoint Higher point of the boundingbox.
     */
    virtual void BoundingBox(
        TPointType& rLowPoint,
        TPointType& rHighPoint
        ) const
    {
        rHighPoint = this->GetPoint( 0 );
        rLowPoint  = this->GetPoint( 0 );
        const SizeType dim = WorkingSpaceDimension();

        for ( IndexType point = 1; point < PointsNumber(); ++point ) { //The first node is already assigned, so we can start from 1
            const auto& r_point = this->GetPoint( point );
            for ( IndexType i = 0; i < dim; ++i ) {
                rHighPoint[i] = ( rHighPoint[i] < r_point[i] ) ? r_point[i] : rHighPoint[i];
                rLowPoint[i]  = ( rLowPoint[i]  > r_point[i] ) ? r_point[i] : rLowPoint[i];
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
     * @brief It returns a vector that is normal to its corresponding geometry in the given local point
     * @param rPointLocalCoordinates Reference to the local coordinates of the point in where the normal is to be computed
     * @return The normal in the given point
     */
    virtual array_1d<double, 3> Normal(const CoordinatesArrayType& rPointLocalCoordinates) const
    {
        const SizeType local_space_dimension = this->LocalSpaceDimension();
        const SizeType dimension = this->WorkingSpaceDimension();

        KRATOS_ERROR_IF(dimension == local_space_dimension) << "Remember the normal can be computed just in geometries with a local dimension: "<< this->LocalSpaceDimension() << "smaller than the spatial dimension: " << this->WorkingSpaceDimension() << std::endl;

        // We define the normal and tangents
        array_1d<double,3> tangent_xi = ZeroVector(3);
        array_1d<double,3> tangent_eta = ZeroVector(3);

        Matrix j_node = ZeroMatrix( dimension, local_space_dimension );
        this->Jacobian( j_node, rPointLocalCoordinates);

        // Using the Jacobian tangent directions
        if (dimension == 2) {
            tangent_eta[2] = 1.0;
            for (unsigned int i_dim = 0; i_dim < dimension; i_dim++) {
                tangent_xi[i_dim]  = j_node(i_dim, 0);
            }
        } else {
            for (unsigned int i_dim = 0; i_dim < dimension; i_dim++) {
                tangent_xi[i_dim]  = j_node(i_dim, 0);
                tangent_eta[i_dim] = j_node(i_dim, 1);
            }
        }

        array_1d<double, 3> normal;
        MathUtils<double>::CrossProduct(normal, tangent_xi, tangent_eta);
        return normal;
    }

    /**
     * @brief It returns the vector, which is normal to its corresponding
     *        geometry in the given integration point for the default
     *        integration method.
     * @param IntegrationPointIndex index in internal integration point list
     * @return The normal in the given integration point
     */
    virtual array_1d<double, 3> Normal(
        IndexType IntegrationPointIndex) const
    {
        return Normal(IntegrationPointIndex, mpGeometryData->DefaultIntegrationMethod());
    }

    /**
     * @brief It returns the vector, which is normal to its corresponding
     *        geometry in the given integration point.
     * @param IntegrationPointIndex index in internal integration point list
     * @param ThisMethod the integration point is dependent on the used
     *        integration method
     * @return The normal in the given integration point
     */
    virtual array_1d<double, 3> Normal(
        IndexType IntegrationPointIndex,
        IntegrationMethod ThisMethod) const
    {
        const SizeType local_space_dimension = this->LocalSpaceDimension();
        const SizeType dimension = this->WorkingSpaceDimension();

        KRATOS_DEBUG_ERROR_IF(dimension == local_space_dimension)
            << "Remember the normal can be computed just in geometries with a local dimension: "
            << this->LocalSpaceDimension() << "smaller than the spatial dimension: "
            << this->WorkingSpaceDimension() << std::endl;

        // We define the normal and tangents
        array_1d<double, 3> tangent_xi = ZeroVector(3);
        array_1d<double, 3> tangent_eta = ZeroVector(3);

        Matrix j_node = ZeroMatrix(dimension, local_space_dimension);
        this->Jacobian(j_node, IntegrationPointIndex, ThisMethod);

        // Using the Jacobian tangent directions
        if (dimension == 2) {
            tangent_eta[2] = 1.0;
            for (IndexType i_dim = 0; i_dim < dimension; i_dim++) {
                tangent_xi[i_dim] = j_node(i_dim, 0);
            }
        }
        else {
            for (IndexType i_dim = 0; i_dim < dimension; i_dim++) {
                tangent_xi[i_dim] = j_node(i_dim, 0);
                tangent_eta[i_dim] = j_node(i_dim, 1);
            }
        }

        array_1d<double, 3> normal;
        MathUtils<double>::CrossProduct(normal, tangent_xi, tangent_eta);
        return normal;
    }

    /**
     * @brief It computes the unit normal of the geometry in the given local point
     * @param rPointLocalCoordinates Refernce to the local coordinates of the point in where the unit normal is to be computed
     * @return The unit normal in the given point
     */
    virtual array_1d<double, 3> UnitNormal(
        const CoordinatesArrayType& rPointLocalCoordinates) const
    {
        array_1d<double, 3> normal = Normal(rPointLocalCoordinates);
        const double norm_normal = norm_2(normal);
        if (norm_normal > std::numeric_limits<double>::epsilon()) normal /= norm_normal;
        else KRATOS_ERROR << "ERROR: The normal norm is zero or almost zero. Norm. normal: " << norm_normal << std::endl;
        return normal;
    }

    /**
     * @brief It returns the normalized normal vector
     *        in the given integration point.
     * @param IntegrationPointIndex index in internal integration point list
     * @param ThisMethod the integration point is dependent on the used
     *        integration method
     * @return The normal in the given integration point
     */
    virtual array_1d<double, 3> UnitNormal(
        IndexType IntegrationPointIndex) const
    {
        return UnitNormal(IntegrationPointIndex, mpGeometryData->DefaultIntegrationMethod());
    }

    /**
     * @brief It returns the normalized normal vector
     *        in the given integration point.
     * @param IntegrationPointIndex index in internal integration point list
     * @param ThisMethod the integration point is dependent on the used
     *        integration method
     * @return The normal in the given integration point
     */
    virtual array_1d<double, 3> UnitNormal(
        IndexType IntegrationPointIndex,
        IntegrationMethod ThisMethod) const
    {
        array_1d<double, 3> normal_vector = Normal(IntegrationPointIndex, ThisMethod);
        const double norm_normal = norm_2(normal_vector);
        if (norm_normal > std::numeric_limits<double>::epsilon())
            normal_vector /= norm_normal;
        else
            KRATOS_ERROR
            << "ERROR: The normal norm is zero or almost zero: "
            << norm_normal << std::endl;
        return normal_vector;
    }

    ///@}
    ///@name Quality
    ///@{

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
       } else if(qualityCriteria == QualityCriteria::MIN_DIHEDRAL_ANGLE) {
         quality = MinDihedralAngle();
       } else if (qualityCriteria == QualityCriteria::MAX_DIHEDRAL_ANGLE) {
         quality = MaxDihedralAngle();
       } else if(qualityCriteria == QualityCriteria::MIN_SOLID_ANGLE) {
         quality = MinSolidAngle();
       }

       return quality;
     }

    /** Calculates the dihedral angles of the geometry.
     * Calculates the dihedral angles of the geometry.
     *
     * @return a vector of dihedral angles of the geometry..
     */
    virtual inline void ComputeDihedralAngles(Vector& rDihedralAngles )  const
    {
        KRATOS_ERROR << "Called the virtual function for ComputeDihedralAngles " << *this << std::endl;
    }

    /** Calculates the solid angles of the geometry.
     * Calculates the solid angles of the geometry.
     *
     * @return a vector of dihedral angles of the geometry..
     */
    virtual inline void ComputeSolidAngles(Vector& rSolidAngles )  const
    {
        KRATOS_ERROR << "Called the virtual function for ComputeDihedralAngles " << *this << std::endl;
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
        return mPoints;
    }

    /** An access method to the Vector of the points stored in
    this geometry.

    @return A reference to PointsArrayType contains pointers to
    the points.
    */
    PointsArrayType& Points()
    {
        return mPoints;
    }

    /** A constant access method to the i'th points stored in
    this geometry.

    @return A constant counted pointer to i'th point of
    geometry.
    */
    const typename TPointType::Pointer pGetPoint( const int Index ) const
    {
        KRATOS_TRY
        return mPoints( Index );
        KRATOS_CATCH(mPoints)
    }

    /** An access method to the i'th points stored in
    this geometry.

    @return A counted pointer to i'th point of
    geometry.
    */
    typename TPointType::Pointer pGetPoint( const int Index )
    {
        KRATOS_TRY
        return mPoints( Index );
        KRATOS_CATCH(mPoints);
    }

    /** A constant access method to the i'th points stored in
    this geometry.

    @return A constant counted pointer to i'th point of
    geometry.
    */
    TPointType const& GetPoint( const int Index ) const
    {
        KRATOS_TRY
        return mPoints[Index];
        KRATOS_CATCH( mPoints);
    }


    /** An access method to the i'th points stored in
    this geometry.

    @return A counted pointer to i'th point of
    geometry.
    */
    TPointType& GetPoint( const int Index )
    {
        KRATOS_TRY
        return mPoints[Index];
        KRATOS_CATCH(mPoints);
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

    ///@}
    ///@name IsInside
    ///@{

    /**
    * @brief Checks if given point in global space coordinates
    *        is inside the geometry boundaries. This function
    *        computes the local coordinates and checks then if
    *        this point lays within the boundaries.
    * @param rPointGlobalCoordinates the global coordinates of the
    *        external point.
    * @param rResult the local coordinates of the point.
    * @param Tolerance the tolerance to the boundary.
    * @return true if the point is inside, false otherwise
    */
    virtual bool IsInside(
        const CoordinatesArrayType& rPointGlobalCoordinates,
        CoordinatesArrayType& rResult,
        const double Tolerance = std::numeric_limits<double>::epsilon()
        ) const
    {
        PointLocalCoordinates(
            rResult,
            rPointGlobalCoordinates);

        if (IsInsideLocalSpace(rResult, Tolerance) == 0) {
            return false;
        }
        return true;
    }

    /**
    * @brief Checks if given point in local space coordinates of this geometry
    *        is inside the geometry boundaries.
    * @param rPointLocalCoordinates the point on the geometry,
    *        which shall be checked if it lays within
    *        the boundaries.
    * @param Tolerance the tolerance to the boundary.
    * @return -1 -> failed
    *          0 -> outside
    *          1 -> inside
    *          2 -> on the boundary
    */
    virtual int IsInsideLocalSpace(
        const CoordinatesArrayType& rPointLocalCoordinates,
        const double Tolerance = std::numeric_limits<double>::epsilon()
    ) const
    {
        KRATOS_ERROR << "Calling IsInsideLocalSpace from base class."
            << " Please check the definition of derived class. "
            << *this << std::endl;
        return 0;
    }

    ///@}
    ///@name Spans
    ///@{

    /* @brief Provides spans in local paramater coordinates of the geometry
     *        according to its direction from LocalDirectionIndex.
     *        For NurbsSurface this is equivalent to the knot vector per direction,
     *        whereby NurbsCurve also provide its knot vector.
     *        Linear geometries shall provide the delimiters, which might be -1 and 1
     *        in standard cases.
     *        Qudartic geometries, may provide additionally the middle point.
     *
     * @param resulting vector of span intervals.
     * @param LocalDirectionIndex of chosen direction, for curves always 0.
     */
    virtual void SpansLocalSpace(
        std::vector<double>& rSpans,
        IndexType LocalDirectionIndex = 0) const
    {
        KRATOS_ERROR <<
            "Calling SpansLocalSpace of geometry base class. Please check derived definitions. "
            << *this << std::endl;
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

    /// Provides the default integration per geometry.
    virtual IntegrationInfo GetDefaultIntegrationInfo() const
    {
        return IntegrationInfo(LocalSpaceDimension(), GetDefaultIntegrationMethod());
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
    ///@name Boundaries
    ///@{

    /**
     * @brief This method gives you all boundaries entities of this geometry.
     * @details This method will gives you all the boundaries entities
     * @return GeometriesArrayType containes this geometry boundaries entities.
     * @see GeneratePoints()
     * @see GenerateEdges()
     * @see GenerateFaces()
     */
    virtual GeometriesArrayType GenerateBoundariesEntities() const
    {
        const SizeType dimension = this->LocalSpaceDimension();
        if (dimension == 3) {
            return this->GenerateFaces();
        } else if (dimension == 2) {
            return this->GenerateEdges();
        } else { // Let's assume is one
            return this->GeneratePoints();
        }
    }

    ///@}
    ///@name Points
    ///@{

    /**
     * @brief This method gives you all points of this geometry.
     * @details This method will gives you all the points
     * @return GeometriesArrayType containes this geometry points.
     * @see Points()
     */
    virtual GeometriesArrayType GeneratePoints() const
    {
        GeometriesArrayType points;

        const auto& p_points = this->Points();
        for (IndexType i_point = 0; i_point < p_points.size(); ++i_point) {
            PointsArrayType point_array;
            point_array.push_back(p_points(i_point));
            auto p_point_geometry = Kratos::make_shared<Geometry<TPointType>>(point_array);
            points.push_back(p_point_geometry);
        }

        return points;
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
    virtual SizeType EdgesNumber() const
    {
        KRATOS_ERROR << "Calling base class EdgesNumber method instead of derived class one. Please check the definition of derived class. " << *this << std::endl;
    }

    /**
     * @brief This method gives you all edges of this geometry.
     * @details This method will gives you all the edges with one dimension less than this geometry.
     * For example a triangle would return three lines as its edges or a tetrahedral would return four triangle as its edges but won't return its six edge lines by this method.
     * @return GeometriesArrayType containes this geometry edges.
     * @deprecated This is legacy version, move to GenerateFaces
     * @see EdgesNumber()
     * @see Edge()
     */
    KRATOS_DEPRECATED_MESSAGE("This is legacy version (use GenerateEdges instead)") virtual GeometriesArrayType Edges( void )
    {
        return this->GenerateEdges();
    }

    /**
     * @brief This method gives you all edges of this geometry.
     * @details This method will gives you all the edges with one dimension less than this geometry.
     * For example a triangle would return three lines as its edges or a tetrahedral would return four triangle as its edges but won't return its six edge lines by this method.
     * @return GeometriesArrayType containes this geometry edges.
     * @see EdgesNumber()
     * @see Edge()
     */
    virtual GeometriesArrayType GenerateEdges() const
    {
        KRATOS_ERROR << "Calling base class Edges method instead of derived class one. Please check the definition of derived class. " << *this << std::endl;
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
    // Commented for possible change in Edge interface of geometry. Pooyan. // NOTE: We should rethink this because is aligned with the current PR
//       virtual Pointer Edge(const PointsArrayType& EdgePoints)
//  {
//    KRATOS_ERROR << "Calling base class Edge method instead of derived class one. Please check the definition of derived class." << *this << std::endl;
//
//  }

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
    virtual SizeType FacesNumber() const
    {
        KRATOS_ERROR << "Calling base class FacesNumber method instead of derived class one. Please check the definition of derived class. " << *this << std::endl;
    }

    /**
     * @brief Returns all faces of the current geometry.
     * @details This is only implemented for 3D geometries, since 2D geometries only have edges but no faces
     * @return GeometriesArrayType containes this geometry faces.
     * @deprecated This is legacy version, move to GenerateFaces
     * @see EdgesNumber
     * @see Edges
     * @see FacesNumber
     */
    KRATOS_DEPRECATED_MESSAGE("This is legacy version (use GenerateFaces instead)") virtual GeometriesArrayType Faces( void )
    {
        const SizeType dimension = this->LocalSpaceDimension();
        if (dimension == 3) {
            return this->GenerateFaces();
        } else {
            return this->GenerateEdges();
        }
    }

    /**
     * @brief Returns all faces of the current geometry.
     * @details This is only implemented for 3D geometries, since 2D geometries only have edges but no faces
     * @return GeometriesArrayType containes this geometry faces.
     * @see EdgesNumber
     * @see GenerateEdges
     * @see FacesNumber
     */
    virtual GeometriesArrayType GenerateFaces() const
    {
        KRATOS_ERROR << "Calling base class GenerateFaces method instead of derived class one. Please check the definition of derived class. " << *this << std::endl;
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
    // Commented for possible change in Edge interface of geometry. Pooyan. // NOTE: We should rethink this because is aligned with the current PR
//       virtual Pointer Edge(IndexType EdgeIndex)
//  {
//    KRATOS_ERROR << "Calling base class Edge method instead of derived class one. Please check the definition of derived class." << *this << std::endl;

//  }

    /** This method gives you normal edge of this geometry which holds
    given points.

    @return NormalType which is normal to this geometry specific edge.
    @see Edge()
    */
    // Commented for possible change in Edge interface of geometry. Pooyan. // NOTE: We should rethink this because is aligned with the current PR
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
    // Commented for possible change in Edge interface of geometry. Pooyan. // NOTE: We should rethink this because is aligned with the current PR
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

    /* Creates integration points according to its quadrature rule.
     * @return integration points.
     */
    virtual void CreateIntegrationPoints(
        IntegrationPointsArrayType& rIntegrationPoints,
        IntegrationInfo& rIntegrationInfo) const
    {
        IntegrationMethod integration_method = rIntegrationInfo.GetIntegrationMethod(0);
        for (IndexType i = 1; i < LocalSpaceDimension(); ++i) {
            KRATOS_ERROR_IF(integration_method != rIntegrationInfo.GetIntegrationMethod(i))
                << "Default creation of integration points only valid if integration method is not varying per direction." << std::endl;
        }
        rIntegrationPoints = IntegrationPoints(integration_method);
    }

    ///@}
    ///@name Quadrature Point Geometries
    ///@{

    /* @brief This method creates a list of quadrature point geometries
     *        from a list of integration points.
     *
     * @param rResultGeometries list of quadrature point geometries.
     * @param rIntegrationPoints list of integration points.
     * @param NumberOfShapeFunctionDerivatives the number of evaluated
     *        derivatives of shape functions at the quadrature point geometries.
     *
     * @see quadrature_point_geometry.h
     */
    virtual void CreateQuadraturePointGeometries(
        GeometriesArrayType& rResultGeometries,
        IndexType NumberOfShapeFunctionDerivatives,
        const IntegrationPointsArrayType& rIntegrationPoints,
        IntegrationInfo& rIntegrationInfo)
    {
        KRATOS_ERROR << "Calling CreateQuadraturePointGeometries from geometry base class."
            << " Please check the definition of derived class. "
            << *this << std::endl;
    }

    /* @brief This method creates a list of quadrature point geometries
     *        from a list of integration points. It creates the list of
     *        integration points byitself.
     *
     * @param rResultGeometries list of quadrature point geometries.
     * @param NumberOfShapeFunctionDerivatives the number of evaluated
     *        derivatives of shape functions at the quadrature point geometries.
     *
     * @see quadrature_point_geometry.h
     */
    virtual void CreateQuadraturePointGeometries(
        GeometriesArrayType& rResultGeometries,
        IndexType NumberOfShapeFunctionDerivatives,
        IntegrationInfo& rIntegrationInfo)
    {
        IntegrationPointsArrayType IntegrationPoints;
        CreateIntegrationPoints(IntegrationPoints, rIntegrationInfo);

        this->CreateQuadraturePointGeometries(
            rResultGeometries,
            NumberOfShapeFunctionDerivatives,
            IntegrationPoints,
            rIntegrationInfo);
    }

    ///@}
    ///@name Operation within Global Space
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

    /** This method provides the global coordinates to
    *   the corresponding integration point
    * @param rResult The global coordinates
    * @param IntegrationPointIndex The index of the integration point
    * @return the global coordinates
    */
    void GlobalCoordinates(
        CoordinatesArrayType& rResult,
        IndexType IntegrationPointIndex
        ) const
    {
        this->GlobalCoordinates(rResult, IntegrationPointIndex, GetDefaultIntegrationMethod());
    }

    /**
    * @brief This method provides the global coordinates to the corresponding integration point
    * @param rResult The global coordinates
    * @param IntegrationPointIndex The index of the integration point
    * @param ThisMethod The integration method
    * @return The global coordinates
    */
    void GlobalCoordinates(
        CoordinatesArrayType& rResult,
        IndexType IntegrationPointIndex,
        const IntegrationMethod ThisMethod
        ) const
    {
        noalias(rResult) = ZeroVector(3);

        const Matrix& N = this->ShapeFunctionsValues(ThisMethod);

        for (IndexType i = 0; i < this->size(); i++)
            noalias(rResult) += N(IntegrationPointIndex, i) * (*this)[i];
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

    /**
    * @brief This method maps from dimension space to working space and computes the
    *        number of derivatives at the dimension parameter.
    * @param rGlobalSpaceDerivatives The derivative in global space.
    * @param rLocalCoordinates the local coordinates
    * @param rDerivativeOrder of computed derivatives
    * @return std::vector<array_1d<double, 3>> with the coordinates in working space
    *         The list is structured as following:
    *           [0] - global coordinates
    *           [1 - loc_space_dim] - base vectors (du, dv, dw)
    *           [...] - second order vectors:
    *                       1D: du^2
    *                       2D: du^2, dudv, dv^2
    *                       3D: du^2, dudv, dudw, dv^2, dvdw, dw^2
    *           [...] - third order vectors:
    *                       1D: du^3
    *                       2D: du^3, du^2dv, dudv^2, dv^3
    */
    virtual void GlobalSpaceDerivatives(
        std::vector<CoordinatesArrayType>& rGlobalSpaceDerivatives,
        const CoordinatesArrayType& rLocalCoordinates,
        const SizeType DerivativeOrder) const
    {
        if (DerivativeOrder == 0)
        {
            if (rGlobalSpaceDerivatives.size() != 1)
                rGlobalSpaceDerivatives.resize(1);

            this->GlobalCoordinates(
                rGlobalSpaceDerivatives[0],
                rLocalCoordinates);
        }
        else if (DerivativeOrder == 1)
        {
            const double local_space_dimension = LocalSpaceDimension();
            const SizeType points_number = this->size();

            if (rGlobalSpaceDerivatives.size() != 1 + local_space_dimension)
                rGlobalSpaceDerivatives.resize(1 + local_space_dimension);

            this->GlobalCoordinates(
                rGlobalSpaceDerivatives[0],
                rLocalCoordinates);

            Matrix shape_functions_gradients(points_number, local_space_dimension);
            this->ShapeFunctionsLocalGradients(shape_functions_gradients, rLocalCoordinates);

            for (IndexType i = 0; i < points_number; ++i) {
                const array_1d<double, 3>& r_coordinates = (*this)[i].Coordinates();
                for (IndexType k = 0; k < WorkingSpaceDimension(); ++k) {
                    const double value = r_coordinates[k];
                    for (IndexType m = 0; m < local_space_dimension; ++m) {
                        rGlobalSpaceDerivatives[m + 1][k] += value * shape_functions_gradients(i, m);
                    }
                }
            }

            return;
        }
        else
        {
            KRATOS_ERROR << "Calling GlobalDerivatives within geometry.h."
                << " Please check the definition within derived class. "
                << *this << std::endl;
        }
    }

    /**
    * @brief This method maps from dimension space to working space and computes the
    *        number of derivatives at the dimension parameter.
    * @param IntegrationPointIndex the coordinates of a certain integration point.
    * @param rDerivativeOrder of computed derivatives
    * @return std::vector<array_1d<double, 3>> with the coordinates in working space
    *         The list is structured as following:
    *           [0] - global coordinates
    *           [1 - loc_space_dim] - base vectors
    *           [...] - higher order vectors (2D: du^2, dudv, dv^2)
    */
    virtual void GlobalSpaceDerivatives(
        std::vector<CoordinatesArrayType>& rGlobalSpaceDerivatives,
        IndexType IntegrationPointIndex,
        const SizeType DerivativeOrder) const
    {
        if (DerivativeOrder == 0)
        {
            if (rGlobalSpaceDerivatives.size() != 1)
                rGlobalSpaceDerivatives.resize(1);

            GlobalCoordinates(
                rGlobalSpaceDerivatives[0],
                IntegrationPointIndex);
        }
        else if (DerivativeOrder == 1)
        {
            const double local_space_dimension = LocalSpaceDimension();
            const SizeType points_number = this->size();

            if (rGlobalSpaceDerivatives.size() != 1 + local_space_dimension)
                rGlobalSpaceDerivatives.resize(1 + local_space_dimension);

            this->GlobalCoordinates(
                rGlobalSpaceDerivatives[0],
                IntegrationPointIndex);

            for (IndexType k = 0; k < local_space_dimension; ++k)
            {
                rGlobalSpaceDerivatives[1 + k] = ZeroVector(3);
            }

            const Matrix& r_shape_functions_gradient_in_integration_point = this->ShapeFunctionLocalGradient(IntegrationPointIndex);

            for (IndexType i = 0; i < points_number; ++i) {
                const array_1d<double, 3>& r_coordinates = (*this)[i].Coordinates();
                for (IndexType k = 0; k < WorkingSpaceDimension(); ++k) {
                    const double value = r_coordinates[k];
                    for (IndexType m = 0; m < local_space_dimension; ++m) {
                        rGlobalSpaceDerivatives[m + 1][k] += value * r_shape_functions_gradient_in_integration_point(i, m);
                    }
                }
            }
        }
        else
        {
            KRATOS_ERROR << "Calling GlobalDerivatives within geometry.h."
                << " Please check the definition within derived class. "
                << *this << std::endl;
        }
    }

    ///@}
    ///@name Spatial Operations
    ///@{

    /**
    * @brief Projects a certain point on the geometry, or finds
    *        the closest point, depending on the provided
    *        initial guess. The external point does not necessary
    *        lay on the geometry.
    *        It shall deal as the interface to the mathematical
    *        projection function e.g. the Newton-Raphson.
    *        Thus, the breaking criteria does not necessarily mean
    *        that it found a point on the surface, if it is really
    *        the closest if or not. It shows only if the breaking
    *        criteria, defined by the tolerance is reached.
    *
    *        This function requires an initial guess, provided by
    *        rProjectedPointLocalCoordinates.
    *        This function can be a very costly operation.
    *
    * @param rPointGlobalCoordinates the point to which the
    *        projection has to be found.
    * @param rProjectedPointGlobalCoordinates the location of the
    *        projection in global coordinates.
    * @param rProjectedPointLocalCoordinates the location of the
    *        projection in local coordinates.
    *        The variable is as initial guess!
    * @param Tolerance accepted of orthogonal error to projection.
    * @return It is chosen to take an int as output parameter to
    *         keep more possibilities within the interface.
    *         0 -> failed
    *         1 -> converged
    */
    KRATOS_DEPRECATED_MESSAGE("This method is deprecated. Use either \'ProjectionPointLocalToLocalSpace\' or \'ProjectionPointGlobalToLocalSpace\' instead.")
    virtual int ProjectionPoint(
        const CoordinatesArrayType& rPointGlobalCoordinates,
        CoordinatesArrayType& rProjectedPointGlobalCoordinates,
        CoordinatesArrayType& rProjectedPointLocalCoordinates,
        const double Tolerance = std::numeric_limits<double>::epsilon()
    ) const
    {
        KRATOS_ERROR << "Calling ProjectionPoint within geometry base class."
            << " Please check the definition within derived class. "
            << *this << std::endl;
    }

    /**
     * @brief Projects a point onto the geometry
     * Projects a certain point on the geometry, or finds the closest point, depending on the provided initial guess.
     * The external point does not necessary lay on the geometry.
     * It shall deal as the interface to the mathematical projection function e.g. the Newton-Raphson.
     * Thus, the breaking criteria does not necessarily mean that it found a point on the surface, if it is really
     * the closest if or not.
     * It shows only if the breaking criteria, defined by the tolerance is reached.
     * This function requires an initial guess, provided by rProjectionPointLocalCoordinates.
     * This function can be a very costly operation.
     * @param rPointLocalCoordinates Local coordinates of the point to be projected
     * @param rProjectionPointLocalCoordinates Projection point local coordinates. This should be initialized with the initial guess
     * @param Tolerance Accepted orthogonal error
     * @return int 0 -> failed
     *             1 -> converged
     */
    virtual int ProjectionPointLocalToLocalSpace(
        const CoordinatesArrayType& rPointLocalCoordinates,
        CoordinatesArrayType& rProjectionPointLocalCoordinates,
        const double Tolerance = std::numeric_limits<double>::epsilon()
    ) const
    {
        KRATOS_ERROR << "Calling ProjectionPointLocalToLocalSpace within geometry base class."
            << " Please check the definition within derived class. "
            << *this << std::endl;
    }

    /**
     * @brief Projects a point onto the geometry
     * Projects a certain point on the geometry, or finds the closest point, depending on the provided initial guess.
     * The external point does not necessary lay on the geometry.
     * It shall deal as the interface to the mathematical projection function e.g. the Newton-Raphson.
     * Thus, the breaking criteria does not necessarily mean that it found a point on the surface, if it is really
     * the closest if or not.
     * It shows only if the breaking criteria, defined by the tolerance is reached.
     * This function requires an initial guess, provided by rProjectionPointLocalCoordinates.
     * This function can be a very costly operation.
     * @param rPointLocalCoordinates Global coordinates of the point to be projected
     * @param rProjectionPointLocalCoordinates Projection point local coordinates. This should be initialized with the initial guess
     * @param Tolerance Accepted orthogonal error
     * @return int 0 -> failed
     *             1 -> converged
     */
    virtual int ProjectionPointGlobalToLocalSpace(
        const CoordinatesArrayType& rPointGlobalCoordinates,
        CoordinatesArrayType& rProjectionPointLocalCoordinates,
        const double Tolerance = std::numeric_limits<double>::epsilon()
    ) const
    {
        KRATOS_ERROR << "Calling ProjectionPointGlobalToLocalSpace within geometry base class."
            << " Please check the definition within derived class. "
            << *this << std::endl;
    }

    /**
    * @brief Returns all coordinates of the closest point on
    *        the geometry given to an arbitrary point in global coordinates.
    *        The basic concept is to first do a projection towards
    *        this geometry and second checking if the projection
    *        was successfull or if no point on the geometry was found.
    * @param rPointGlobalCoordinates the point to which the
    *        closest point has to be found.
    * @param rClosestPointGlobalCoordinates the location of the
    *        closest point in global coordinates.
    * @param rClosestPointLocalCoordinates the location of the
    *        closest point in local coordinates.
    *        IMPORTANT: The variable can also be used as initial guess.
    * @param Tolerance accepted orthogonal error.
    * @return -1 -> failed
    *          0 -> outside
    *          1 -> inside
    *          2 -> on the boundary
    */
    KRATOS_DEPRECATED_MESSAGE("This method is deprecated. Use either \'ClosestPointLocalToLocalSpace\' or \'ClosestPointGlobalToLocalSpace\' instead. Please note that \'rClosestPointGlobalCoordinates\' returns unmodified original value.")
    virtual int ClosestPoint(
        const CoordinatesArrayType& rPointGlobalCoordinates,
        CoordinatesArrayType& rClosestPointGlobalCoordinates,
        CoordinatesArrayType& rClosestPointLocalCoordinates,
        const double Tolerance = std::numeric_limits<double>::epsilon()
    ) const
    {
        return ClosestPointGlobalToLocalSpace(rPointGlobalCoordinates, rClosestPointLocalCoordinates, Tolerance);
    }

    /**
    * @brief Returns global coordinates of the closest point on
    *        the geometry given to an arbitrary point in global coordinates.
    *        The basic concept is to first do a projection towards
    *        this geometry and second checking if the projection
    *        was successfull or if no point on the geometry was found.
    * @param rPointGlobalCoordinates the point to which the
    *        closest point has to be found.
    * @param rClosestPointGlobalCoordinates the location of the
    *        closest point in global coordinates.
    *
    *        WARNING: This function does not provide the possibility
    *                 to use an initial guess!!
    *
    * @param Tolerance accepted orthogonal error.
    * @return -1 -> failed
    *          0 -> outside
    *          1 -> inside
    *          2 -> on the boundary
    */
    KRATOS_DEPRECATED_MESSAGE("This method is deprecated. Use either \'ClosestPointLocalToLocalSpace\' or \'ClosestPointGlobalToLocalSpace\' instead.")
    virtual int ClosestPoint(
        const CoordinatesArrayType& rPointGlobalCoordinates,
        CoordinatesArrayType& rClosestPointGlobalCoordinates,
        const double Tolerance = std::numeric_limits<double>::epsilon()
    ) const
    {
        CoordinatesArrayType local_coordinates(ZeroVector(3));
        const int result = ClosestPointGlobalToLocalSpace(rPointGlobalCoordinates, local_coordinates, Tolerance);

        if (result == 1) {
            this->GlobalCoordinates(rClosestPointGlobalCoordinates, local_coordinates);
        }

        return result;
    }

    /**
    * @brief Returns local coordinates of the closest point on
    *        the geometry given to an arbitrary point in global coordinates.
    *        The basic concept is to first do a projection towards
    *        this geometry and second checking if the projection
    *        was successfull or if no point on the geometry was found.
    * @param rPointGlobalCoordinates the point to which the
    *        closest point has to be found.
    * @param rClosestPointLocalCoordinates the location of the
    *        closest point in local coordinates.
    *
    *        IMPORTANT: The rClosestPointLocalCoordinates can
    *                   also be used as initial guess.
    *
    * @param Tolerance accepted orthogonal error.
    * @return -1 -> failed
    *          0 -> outside
    *          1 -> inside
    *          2 -> on the boundary
    */
    KRATOS_DEPRECATED_MESSAGE("This method is deprecated. Use either \'ClosestPointLocalToLocalSpace\' or \'ClosestPointGlobalToLocalSpace\' instead.")
    virtual int ClosestPointLocalCoordinates(
        const CoordinatesArrayType& rPointGlobalCoordinates,
        CoordinatesArrayType& rClosestPointLocalCoordinates,
        const double Tolerance = std::numeric_limits<double>::epsilon()
    ) const
    {
        return ClosestPointGlobalToLocalSpace(rPointGlobalCoordinates, rClosestPointLocalCoordinates, Tolerance);
    }

    /**
     * @brief Calculates the closes point projection
     * This method calculates the closest point projection of a point in local space coordinates
     * @param rPointLocalCoordinates Input local coordinates
     * @param rClosestPointLocalCoordinates Closest point local coordinates. This should be initialized with the initial guess
     * @param Tolerance Accepted orthogonal error
     * @return int -1 -> failed
     *             0 -> outside
     *             1 -> inside
     *             2 -> on the boundary
     */
    virtual int ClosestPointLocalToLocalSpace(
        const CoordinatesArrayType& rPointLocalCoordinates,
        CoordinatesArrayType& rClosestPointLocalCoordinates,
        const double Tolerance = std::numeric_limits<double>::epsilon()
    ) const
    {
        // 1. Make projection on geometry
        const int projection_result = ProjectionPointLocalToLocalSpace(
            rPointLocalCoordinates,
            rClosestPointLocalCoordinates,
            Tolerance);

        if (projection_result == 1) {
            // 2. If projection converged check if solution lays
            // within the boundaries of this geometry
            // Returns either 0, 1 or 2
            // Or -1 if IsInsideLocalSpace failed
            return IsInsideLocalSpace(
                rClosestPointLocalCoordinates,
                Tolerance);
        } else {
            // Projection failed
            return -1;
        }
    }

    /**
     * @brief Calculates the closes point projection
     * This method calculates the closest point projection of a point in global space coordinates
     * @param rPointLocalCoordinates Input global coordinates
     * @param rClosestPointLocalCoordinates Closest point local coordinates. This should be initialized with the initial guess
     * @param Tolerance Accepted orthogonal error
     * @return int -1 -> failed
     *             0 -> outside
     *             1 -> inside
     *             2 -> on the boundary
     */
    virtual int ClosestPointGlobalToLocalSpace(
        const CoordinatesArrayType& rPointGlobalCoordinates,
        CoordinatesArrayType& rClosestPointLocalCoordinates,
        const double Tolerance = std::numeric_limits<double>::epsilon()
    ) const
    {
        // 1. Make projection on geometry
        const int projection_result = ProjectionPointGlobalToLocalSpace(
            rPointGlobalCoordinates,
            rClosestPointLocalCoordinates,
            Tolerance);

        if (projection_result == 1) {
            // 2. If projection converged check if solution lays
            // within the boundaries of this geometry
            // Returns either 0, 1 or 2
            // Or -1 if IsInsideLocalSpace failed
            return IsInsideLocalSpace(
                rClosestPointLocalCoordinates,
                Tolerance);
        } else {
            // Projection failed
            return -1;
        }
    }

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
    virtual double CalculateDistance(
        const CoordinatesArrayType& rPointGlobalCoordinates,
        const double Tolerance = std::numeric_limits<double>::epsilon()
    ) const
    {
        CoordinatesArrayType local_coordinates(ZeroVector(3));
        if (ClosestPointGlobalToLocalSpace(rPointGlobalCoordinates, local_coordinates, Tolerance) < 1) {
            // If projection fails, double::max will be returned
            return std::numeric_limits<double>::max();
        }

        // Global coordinates of projected point
        CoordinatesArrayType global_coordinates(ZeroVector(3));
        this->GlobalCoordinates(global_coordinates, local_coordinates);

        // Distance to projected point
        return norm_2(rPointGlobalCoordinates - global_coordinates);
    }

    ///@}
    ///@name Jacobian
    ///@{

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

        for ( unsigned int pnt = 0; pnt < this->IntegrationPointsNumber( ThisMethod ); pnt++ ) {
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

        for ( unsigned int pnt = 0; pnt < this->IntegrationPointsNumber( ThisMethod ); pnt++ ) {
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
        const SizeType working_space_dimension = this->WorkingSpaceDimension();
        const SizeType local_space_dimension = this->LocalSpaceDimension();
        if(rResult.size1() != working_space_dimension || rResult.size2() != local_space_dimension)
            rResult.resize( working_space_dimension, local_space_dimension, false );

        const Matrix& r_shape_functions_gradient_in_integration_point = ShapeFunctionsLocalGradients( ThisMethod )[ IntegrationPointIndex ];

        rResult.clear();
        const SizeType points_number = this->PointsNumber();
        for (IndexType i = 0; i < points_number; ++i ) {
            const array_1d<double, 3>& r_coordinates = (*this)[i].Coordinates();
            for(IndexType k = 0; k< working_space_dimension; ++k) {
                const double value = r_coordinates[k];
                for(IndexType m = 0; m < local_space_dimension; ++m) {
                    rResult(k,m) += value * r_shape_functions_gradient_in_integration_point(i,m);
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

    @param rDeltaPosition Matrix with the nodes position increment which describes
    the configuration where the jacobian has to be calculated.

    @return Matrix<double> Jacobian matrix \f$ J_i \f$ where \f$
    i \f$ is the given integration point index of given
    integration method.

    @see DeterminantOfJacobian
    @see InverseOfJacobian
    */
    virtual Matrix& Jacobian( Matrix& rResult, IndexType IntegrationPointIndex, IntegrationMethod ThisMethod, const Matrix& rDeltaPosition ) const
    {
        const SizeType working_space_dimension = this->WorkingSpaceDimension();
        const SizeType local_space_dimension = this->LocalSpaceDimension();
        if(rResult.size1() != working_space_dimension || rResult.size2() != local_space_dimension)
            rResult.resize( working_space_dimension, local_space_dimension, false );

        const Matrix& r_shape_functions_gradient_in_integration_point = ShapeFunctionsLocalGradients( ThisMethod )[ IntegrationPointIndex ];

        rResult.clear();
        const SizeType points_number = this->PointsNumber();
        for (IndexType i = 0; i < points_number; ++i ) {
            const array_1d<double, 3>& r_coordinates = (*this)[i].Coordinates();
            for(IndexType k = 0; k< working_space_dimension; ++k) {
                const double value = r_coordinates[k] - rDeltaPosition(i,k);
                for(IndexType m = 0; m < local_space_dimension; ++m) {
                    rResult(k,m) += value * r_shape_functions_gradient_in_integration_point(i,m);
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
        const SizeType working_space_dimension = this->WorkingSpaceDimension();
        const SizeType local_space_dimension = this->LocalSpaceDimension();
        const SizeType points_number = this->PointsNumber();
        if(rResult.size1() != working_space_dimension || rResult.size2() != local_space_dimension)
            rResult.resize( working_space_dimension, local_space_dimension, false );

        Matrix shape_functions_gradients(points_number, local_space_dimension);
        ShapeFunctionsLocalGradients( shape_functions_gradients, rCoordinates );

        rResult.clear();
        for (IndexType i = 0; i < points_number; ++i ) {
            const array_1d<double, 3>& r_coordinates = (*this)[i].Coordinates();
            for(IndexType k = 0; k< working_space_dimension; ++k) {
                const double value = r_coordinates[k];
                for(IndexType m = 0; m < local_space_dimension; ++m) {
                    rResult(k,m) += value * shape_functions_gradients(i,m);
                }
            }
        }

        return rResult;
    }

    /** Jacobian in given point. This method calculate jacobian
    matrix in given point.

    @param rCoordinates point which jacobians has to
    be calculated in it.

    @param rDeltaPosition Matrix with the nodes position increment which describes
    the configuration where the jacobian has to be calculated.

    @return Matrix of double which is jacobian matrix \f$ J \f$ in given point.

    @see DeterminantOfJacobian
    @see InverseOfJacobian
    */

    virtual Matrix& Jacobian( Matrix& rResult, const CoordinatesArrayType& rCoordinates, Matrix& rDeltaPosition ) const
    {
        const SizeType working_space_dimension = this->WorkingSpaceDimension();
        const SizeType local_space_dimension = this->LocalSpaceDimension();
        const SizeType points_number = this->PointsNumber();
        if(rResult.size1() != working_space_dimension || rResult.size2() != local_space_dimension)
            rResult.resize( working_space_dimension, local_space_dimension, false );

        Matrix shape_functions_gradients(points_number, local_space_dimension);
        ShapeFunctionsLocalGradients( shape_functions_gradients, rCoordinates );

        rResult.clear();
        for (IndexType i = 0; i < points_number; ++i ) {
            const array_1d<double, 3>& r_coordinates = (*this)[i].Coordinates();
            for(IndexType k = 0; k< working_space_dimension; ++k) {
                const double value = r_coordinates[k] - rDeltaPosition(i,k);
                for(IndexType m = 0; m < local_space_dimension; ++m) {
                    rResult(k,m) += value * shape_functions_gradients(i,m);
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

        Matrix J( this->WorkingSpaceDimension(), this->LocalSpaceDimension());
        for ( unsigned int pnt = 0; pnt < this->IntegrationPointsNumber( ThisMethod ); pnt++ ) {
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
        Matrix J( this->WorkingSpaceDimension(), this->LocalSpaceDimension());
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
        Matrix J( this->WorkingSpaceDimension(), this->LocalSpaceDimension());
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
        Matrix Jinv(this->LocalSpaceDimension(), this->WorkingSpaceDimension());
        for ( unsigned int pnt = 0; pnt < this->IntegrationPointsNumber( ThisMethod ); pnt++ ) {
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

    /*
    * @brief access to the shape function derivatives.
    * @param DerivativeOrderIndex defines the wanted order of the derivative
    *        0 is NOT accessible
    * @param IntegrationPointIndex the corresponding contorl point of this geometry
    * @return the shape function derivative matrix.
    *         The matrix is structured: (derivative dN_de / dN_du , the corresponding node)
    */
    const Matrix& ShapeFunctionDerivatives(
        IndexType DerivativeOrderIndex,
        IndexType IntegrationPointIndex,
        IntegrationMethod ThisMethod) const
    {
        return mpGeometryData->ShapeFunctionDerivatives(
            DerivativeOrderIndex, IntegrationPointIndex, ThisMethod);
    }

    /*
    * @brief access to the shape function derivatives.
    * @param DerivativeOrderIndex defines the wanted order of the derivative
    *        0 is NOT accessible
    * @param IntegrationPointIndex the corresponding contorl point of this geometry
    * @return the shape function derivative matrix.
    *         The matrix is structured: (derivative dN_de / dN_du , the corresponding node)
    */
    const Matrix& ShapeFunctionDerivatives(
        IndexType DerivativeOrderIndex,
        IndexType IntegrationPointIndex) const
    {
        return mpGeometryData->ShapeFunctionDerivatives(
            DerivativeOrderIndex, IntegrationPointIndex, GetDefaultIntegrationMethod());
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


    void ShapeFunctionsIntegrationPointsGradients( ShapeFunctionsGradientsType& rResult ) const
    {
        ShapeFunctionsIntegrationPointsGradients( rResult, mpGeometryData->DefaultIntegrationMethod() );
    }

    virtual void ShapeFunctionsIntegrationPointsGradients(
        ShapeFunctionsGradientsType& rResult,
        IntegrationMethod ThisMethod ) const
    {
        //check that current geometry can perform this operation
        KRATOS_ERROR_IF_NOT(this->WorkingSpaceDimension() == this->LocalSpaceDimension())
            << "\'ShapeFunctionsIntegrationPointsGradients\' is not defined for current geometry type as gradients are only defined in the local space." << std::endl;

        const unsigned int integration_points_number = this->IntegrationPointsNumber( ThisMethod );

        if ( integration_points_number == 0 )
            KRATOS_ERROR << "This integration method is not supported" << *this << std::endl;

        if ( rResult.size() != integration_points_number )
            rResult.resize(  this->IntegrationPointsNumber( ThisMethod ), false  );

        //calculating the local gradients
        const ShapeFunctionsGradientsType& DN_De = ShapeFunctionsLocalGradients( ThisMethod );

        //loop over all integration points
        Matrix Jinv(this->LocalSpaceDimension(), this->WorkingSpaceDimension());
        for ( unsigned int pnt = 0; pnt < integration_points_number; pnt++ ) {
            if (rResult[pnt].size1() != (*this).size() || rResult[pnt].size2() != this->LocalSpaceDimension())
                rResult[pnt].resize( (*this).size(), this->LocalSpaceDimension(), false );
            this->InverseOfJacobian(Jinv,pnt, ThisMethod);
            noalias(rResult[pnt]) =  prod( DN_De[pnt], Jinv );
        }
    }

    virtual void ShapeFunctionsIntegrationPointsGradients(
        ShapeFunctionsGradientsType& rResult,
        Vector& rDeterminantsOfJacobian,
        IntegrationMethod ThisMethod ) const
    {
        //check that current geometry can perform this operation
        KRATOS_ERROR_IF_NOT(this->WorkingSpaceDimension() == this->LocalSpaceDimension())
            << "\'ShapeFunctionsIntegrationPointsGradients\' is not defined for current geometry type as gradients are only defined in the local space." << std::endl;

        const unsigned int integration_points_number = this->IntegrationPointsNumber( ThisMethod );

        if ( integration_points_number == 0 )
            KRATOS_ERROR << "This integration method is not supported " << *this << std::endl;

        if ( rResult.size() != integration_points_number )
            rResult.resize(  this->IntegrationPointsNumber( ThisMethod ), false  );
        if (rDeterminantsOfJacobian.size() != integration_points_number)
            rDeterminantsOfJacobian.resize(this->IntegrationPointsNumber(ThisMethod), false);

        //calculating the local gradients
        const ShapeFunctionsGradientsType& DN_De = ShapeFunctionsLocalGradients( ThisMethod );

        //loop over all integration points
        Matrix J(this->WorkingSpaceDimension(), this->LocalSpaceDimension());
        Matrix Jinv(this->LocalSpaceDimension(), this->WorkingSpaceDimension());
        double DetJ;
        for ( unsigned int pnt = 0; pnt < integration_points_number; pnt++ ) {
            if (rResult[pnt].size1() != (*this).size() || rResult[pnt].size2() != this->LocalSpaceDimension())
                rResult[pnt].resize( (*this).size(), this->LocalSpaceDimension(), false );
            this->Jacobian(J,pnt, ThisMethod);
            MathUtils<double>::GeneralizedInvertMatrix(J, Jinv, DetJ);
            noalias(rResult[pnt]) =  prod( DN_De[pnt], Jinv );
            rDeterminantsOfJacobian[pnt] = DetJ;
        }
    }

    KRATOS_DEPRECATED_MESSAGE("This is signature of \'ShapeFunctionsIntegrationPointsGradients\' is legacy (use any of the alternatives without shape functions calculation).")
    virtual void ShapeFunctionsIntegrationPointsGradients(
        ShapeFunctionsGradientsType& rResult,
        Vector& rDeterminantsOfJacobian,
        IntegrationMethod ThisMethod,
        Matrix& ShapeFunctionsIntegrationPointsValues) const
    {

        ShapeFunctionsIntegrationPointsGradients(rResult, rDeterminantsOfJacobian, ThisMethod);
        ShapeFunctionsIntegrationPointsValues = ShapeFunctionsValues(ThisMethod);
    }

    virtual int Check() const
    {
        KRATOS_TRY

        return 0;

        KRATOS_CATCH("")
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Return geometry information as a string.
    virtual std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Geometry # "
            << std::to_string(mId) << ": "
            << LocalSpaceDimension() << "-dimensional geometry in "
            << WorkingSpaceDimension() << "D space";

        return buffer.str();
    }

    /// Returns name.
    virtual std::string Name() const {
        std::string geometryName = "BaseGeometry";
        KRATOS_ERROR << "Base geometry does not have a name." << std::endl;
        return geometryName;
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info();
    }

    /// Print name.
    virtual void PrintName(std::ostream& rOstream) const {
        rOstream << Name() << std::endl;
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
        if (mpGeometryData) {
            mpGeometryData->PrintData(rOStream);
        }

        rOStream << std::endl;
        rOStream << std::endl;

        for (unsigned int i = 0; i < this->size(); ++i) {
            rOStream << "\tPoint " << i + 1 << "\t : ";
            if (mPoints(i) != nullptr) {
                mPoints[i].PrintData(rOStream);
            } else {
                rOStream << "point is empty (nullptr)." << std::endl;
            }
            rOStream << std::endl;
        }

        if (AllPointsAreValid()) {
            rOStream << "\tCenter\t : ";
            Center().PrintData(rOStream);
        }

        rOStream << std::endl;
        rOStream << std::endl;
    }


    ///@}
    ///@name Friends
    ///@{

    ///@}

protected:
    ///@name Geometry Data
    ///@{

    /**
    * @brief updates the pointer to GeometryData of the
    *        respective geometry.
    * @param pGeometryData pointer to const GeometryData.
    */
    void SetGeometryData(GeometryData const* pGeometryData)
    {
        mpGeometryData = pGeometryData;
    }

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

    /** Calculates the min dihedral angle quality metric.
     * Calculates the min dihedral angle quality metric.
     * The min dihedral angle is min angle between two faces of the element
     * In radians
     * @return [description]
     */
    virtual double MinDihedralAngle() const {
        KRATOS_ERROR << "Calling base class 'MinDihedralAngle' method instead of derived class one. Please check the definition of derived class. " << *this << std::endl;
        return 0.0;
    }

    /** Calculates the max dihedral angle quality metric.
     * Calculates the max dihedral angle quality metric.
     * The max dihedral angle is max angle between two faces of the element
     * In radians
     * @return [description]
     */
    virtual double MaxDihedralAngle() const {
        KRATOS_ERROR << "Calling base class 'MaxDihedralAngle' method instead of derived class one. Please check the definition of derived class. " << *this << std::endl;
        return 0.0;
    }

    /** Calculates the min solid angle quality metric.
     * Calculates the min solid angle quality metric.
     * The min solid angle  [stereoradians] is the lowest solid angle "seen" from any of the 4 nodes of the geometry. Valid only for 3d elems!
     * In stereo radians
     * @return [description]
     */
    virtual double MinSolidAngle() const {
        KRATOS_ERROR << "Calling base class 'MinSolidAngle' method instead of derived class one. Please check the definition of derived class. " << *this << std::endl;
        return 0.0;
    }

    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    /**
     * @brief Checks if the geometry points are valid
     * Checks if the geometry points are valid from the pointer value
     * Points are not valid when the pointer value is null
     * @return true All points are valid
     * @return false At least one point has nullptr value
     */
    bool AllPointsAreValid() const
    {
        return std::none_of(mPoints.ptr_begin(), mPoints.ptr_end(), [](const auto& pPoint){return pPoint == nullptr;});
    }

    ///@}
    ///@name Protected LifeCycle
    ///@{

    /** Protected Constructor.
    Avoids object to be created Except for derived classes
    */

    ///@}
private:
    ///@name Member Variables
    ///@{

    IndexType mId;

    GeometryData const* mpGeometryData;

    static const GeometryDimension msGeometryDimension;

    PointsArrayType mPoints;

    DataValueContainer mData;


    ///@}
    ///@name Id Bit-Change Operations
    ///@{

    /// Gets the corresponding self assigned id from pointer
    IndexType GenerateSelfAssignedId() const
    {
        // Create id hash from provided name.
        IndexType id = reinterpret_cast<IndexType>(this);

        // Sets second bit to zero.
        SetIdSelfAssigned(id);

        // Sets first bit to zero.
        SetIdNotGeneratedFromString(id);

        return id;
    }

    /// Checks first bit in Id. 0 -> id; 1 -> name/ string
    static inline bool IsIdGeneratedFromString(IndexType Id)
    {
        return Id & (IndexType(1) << (sizeof(IndexType) * 8 - 1));
    }

    /// Sets first bit in Id to 1 -> name/ string
    static inline void SetIdGeneratedFromString(IndexType& Id)
    {
        Id |= (IndexType(1) << (sizeof(IndexType) * 8 - 1));
    }

    /// Sets first bit in Id to 0 -> no name/ string
    static inline void SetIdNotGeneratedFromString(IndexType& Id)
    {
        Id &= ~(IndexType(1) << (sizeof(IndexType) * 8 - 1));
    }

    /// Checks second bit in Id. 0 -> defined id; 1 -> self assigned
    static inline bool IsIdSelfAssigned(IndexType Id)
    {
        return Id & (IndexType(1) << (sizeof(IndexType) * 8 - 2));
    }

    /// Sets second bit in Id to 1 -> self assigned
    static inline void SetIdSelfAssigned(IndexType& Id)
    {
        Id |= (IndexType(1) << (sizeof(IndexType) * 8 - 2));
    }

    /// Sets second bit in Id to 0 -> not self assigned
    static inline void SetIdNotSelfAssigned(IndexType& Id)
    {
        Id &= ~(IndexType(1) << (sizeof(IndexType) * 8 - 2));
    }

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    virtual void save( Serializer& rSerializer ) const
    {
        rSerializer.save("Id", mId);
        rSerializer.save( "Points", mPoints);
        rSerializer.save("Data", mData);
    }

    virtual void load( Serializer& rSerializer )
    {
        rSerializer.load("Id", mId);
        rSerializer.load( "Points", mPoints );
        rSerializer.load("Data", mData);
   }

    ///@}
    ///@name Private Operations
    ///@{

    /// Sets Id of this Geometry (avoids checks, can be used only as private)
    void SetIdWithoutCheck(const IndexType Id)
    {
        mId = Id;
    }

    static const GeometryData& GeometryDataInstance()
    {
        IntegrationPointsContainerType integration_points = {};
        ShapeFunctionsValuesContainerType shape_functions_values = {};
        ShapeFunctionsLocalGradientsContainerType shape_functions_local_gradients = {};
        static GeometryData s_geometry_data(
                            &msGeometryDimension,
                            GeometryData::IntegrationMethod::GI_GAUSS_1,
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

template<class TPointType>
const GeometryDimension Geometry<TPointType>::msGeometryDimension(3, 3);

}  // namespace Kratos.
