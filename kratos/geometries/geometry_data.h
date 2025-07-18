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
//
//

#pragma once

// System includes

// External includes

// Project includes
#include "includes/ublas_interface.h"
#include "integration/integration_point.h"
#include "geometries/geometry_dimension.h"
#include "geometries/geometry_shape_function_container.h"

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

/** GeometryData base class.As a base class GeometryData has all the common
    interface of Kratos' geometries. Also it contains array of
    pointers to its points, reference to shape functions values in
    all integrations points and also local gradients of shape
    functions evaluated in all integrations points.

    @see Geometry
    @see Element
*/

class GeometryData
{
public:
    ///@name  Enum's
    ///@{

    /** Integration methods implemented in geometry. Each geometry
    supposed to have different integration method for
    integrating. These enums are used to refer to each
    different integration methods:

    - GI_GAUSS_1 gaussian integration with order 1.
    - GI_GAUSS_2 gaussian integration with order 2.
    - GI_GAUSS_3 gaussian integration with order 3.
    - GI_GAUSS_4 gaussian integration with order 4.
    - GI_GAUSS_5 gaussian integration with order 5.
    */
    enum class IntegrationMethod {
        GI_GAUSS_1,
        GI_GAUSS_2,
        GI_GAUSS_3,
        GI_GAUSS_4,
        GI_GAUSS_5,
        GI_LOBATTO_1,
        NumberOfIntegrationMethods // Note that this entry needs to be always the last to be used as integration methods counter
    };

    enum class KratosGeometryFamily
    {
        Kratos_NoElement,
        Kratos_Point,
        Kratos_Linear,
        Kratos_Triangle,
        Kratos_Quadrilateral,
        Kratos_Tetrahedra,
        Kratos_Hexahedra,
        Kratos_Prism,
        Kratos_Pyramid,
        Kratos_Nurbs,
        Kratos_Brep,
        Kratos_Quadrature_Geometry,
        Kratos_Composite,
        Kratos_generic_family,
        NumberOfGeometryFamilies // Note that this entry needs to be always the last to be used as geometry families counter
    };

    enum class KratosGeometryType
    {
        Kratos_generic_type,
        Kratos_Hexahedra3D20,
        Kratos_Hexahedra3D27,
        Kratos_Hexahedra3D8,
        Kratos_Prism3D15,
        Kratos_Prism3D6,
        Kratos_Pyramid3D13,
        Kratos_Pyramid3D5,
        Kratos_Quadrilateral2D4,
        Kratos_Quadrilateral2D8,
        Kratos_Quadrilateral2D9,
        Kratos_Quadrilateral3D4,
        Kratos_Quadrilateral3D8,
        Kratos_Quadrilateral3D9,
        Kratos_Tetrahedra3D10,
        Kratos_Tetrahedra3D4,
        Kratos_Triangle2D3,
        Kratos_Triangle2D6,
        Kratos_Triangle2D10,
        Kratos_Triangle2D15,
        Kratos_Triangle3D3,
        Kratos_Triangle3D6,
        Kratos_Line2D2,
        Kratos_Line2D3,
        Kratos_Line2D4,
        Kratos_Line2D5,
        Kratos_Line3D2,
        Kratos_Line3D3,
        Kratos_Point2D,
        Kratos_Point3D,
        Kratos_Sphere3D1,
        Kratos_Nurbs_Curve,
        Kratos_Nurbs_Surface,
        Kratos_Nurbs_Volume,
        Kratos_Nurbs_Curve_On_Surface,
        Kratos_Surface_In_Nurbs_Volume,
        Kratos_Brep_Curve,
        Kratos_Brep_Surface,
        Kratos_Brep_Curve_On_Surface,
        Kratos_Quadrature_Point_Geometry,
        Kratos_Coupling_Geometry,
        Kratos_Quadrature_Point_Curve_On_Surface_Geometry,
        Kratos_Quadrature_Point_Surface_In_Volume_Geometry,
        NumberOfGeometryTypes // Note that this entry needs to be always the last to be used as geometry types counter
    };

    /**
     * @enum KratosGeometryOrderType
     * @brief Defines the order of geometries in Kratos.
     * @details This enumeration lists the different polynomial orders that can be used for defining the geometry of elements in Kratos. The geometry order impacts the accuracy and computational cost of simulations.
     * @var KratosGeometryOrderType::Kratos_Zero_Order: Defines a zero-order geometry. This is 0D geometry, aka discrete points, and is the simplest form of geometry.
     * @var KratosGeometryOrderType::Kratos_Linear_Order: Defines a linear geometry. This is the simplest form, with minimal computational cost but lower accuracy.
     * @var KratosGeometryOrderType::Kratos_Quadratic_Order: Defines a quadratic geometry. Offers a balance between computational cost and accuracy.
     * @var KratosGeometryOrderType::Kratos_Cubic_Order: Defines a cubic geometry. Higher computational cost than quadratic, but offers improved accuracy.
     * @var KratosGeometryOrderType::Kratos_Quartic_Order: Defines a quartic geometry. Further increases computational cost and accuracy.
     * @var KratosGeometryOrderType::Kratos_Quintic_Order: Defines a quintic geometry. Represents the highest polynomial order here, with the highest accuracy and computational cost.
     * @var KratosGeometryOrderType::Kratos_Unknown_Order: Used when the geometry order is unknown or not specified. Allows for flexibility in geometry definition.
    */
    enum KratosGeometryOrderType {
        Kratos_Zero_Order,
        Kratos_Linear_Order,
        Kratos_Quadratic_Order,
        Kratos_Cubic_Order,
        Kratos_Quartic_Order,
        Kratos_Quintic_Order,
        Kratos_Unknown_Order
    };

    ///@}
    ///@name Type Definitions
    ///@{

    /// Pointer definition of GeometryData
    KRATOS_CLASS_POINTER_DEFINITION( GeometryData );

    /** Type used for indexing in geometry data class. Unsigned int used for indexing
    point or integration point access methods and also all other
    methods which need point or integration point index.
    */
    using IndexType = std::size_t;

    /** This typed used to return size or dimension in
    geometry data. Dimension, WorkingSpaceDimension, PointsNumber and
    ... return this type as their results.
    */
    using SizeType = std::size_t;

    /** This type used for representing an integration point in
    geometry data. This integration point is a point with an
    additional weight component.
    */
    using IntegrationPointType = IntegrationPoint<3>;

    /** A Vector of IntegrationPointType which used to hold
    integration points related to an integration
    method. IntegrationPoints functions used this type to return
    their results.
    */
    using IntegrationPointsArrayType = GeometryShapeFunctionContainer<IntegrationMethod>::IntegrationPointsArrayType;

    /** A Vector of IntegrationPointsArrayType which used to hold
    integration points related to different integration method
    implemented in geometry.
    */
    using IntegrationPointsContainerType = GeometryShapeFunctionContainer<IntegrationMethod>::IntegrationPointsContainerType;

    /** A third order tensor used as shape functions' values
    container.
    */
    using ShapeFunctionsValuesContainerType = GeometryShapeFunctionContainer<IntegrationMethod>::ShapeFunctionsValuesContainerType;

    /** A fourth order tensor used as shape functions' local
    gradients container in geometry data.
    */
    using ShapeFunctionsLocalGradientsContainerType = GeometryShapeFunctionContainer<IntegrationMethod>::ShapeFunctionsLocalGradientsContainerType;

    /** A third order tensor to hold shape functions'
    gradients. ShapefunctionsLocalGradients function return this
    type as its result.
    */
    using ShapeFunctionsGradientsType = DenseVector<Matrix>;

    using ShapeFunctionsSecondDerivativesType = DenseVector<Matrix>;

    /**
     * fourth order tensor to hold the third order derivatives of the
     * shape functions
     */
    using ShapeFunctionsThirdDerivativesType = DenseVector<DenseVector<Matrix>>;

    ///@}
    ///@name Life Cycle
    ///@{

    /** Complete argument constructor. This constructor gives a
    complete set of arguments to pass all the initial value of
    all the member variables of geometry class. Also it has
    default value for integration variables to make it useful
    in the case of constructing new geometry without mapping and
    integrating properties.

    @param ThisDimension Dimension of this geometry.

    @param ThisWorkingSpaceDimension Working space dimension. for
    example a triangle 3d is a 2 dimensional shape but can be used
    in 3 dimensional space.

    @param ThisLocalSpaceDimension Local space dimension.
    for example a triangle is a 2 dimensional shape but
    can have 3 dimensional area coordinates l1, l2, l3.

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
    have gaussian order "2" ThisShapeFunctionsValues[GI_GAUSS_2]
    must be an empty ShapeFunctionsGradientsType.
    */
    GeometryData(GeometryDimension const *pThisGeometryDimension,
        IntegrationMethod ThisDefaultMethod,
        const IntegrationPointsContainerType& ThisIntegrationPoints,
        const ShapeFunctionsValuesContainerType& ThisShapeFunctionsValues,
        const ShapeFunctionsLocalGradientsContainerType& ThisShapeFunctionsLocalGradients)
        : mpGeometryDimension(pThisGeometryDimension)
        , mGeometryShapeFunctionContainer(
                ThisDefaultMethod,
                ThisIntegrationPoints,
                ThisShapeFunctionsValues,
                ThisShapeFunctionsLocalGradients)
    {
    }

    /*
    * Constructor which has a precomputed shape function container.
    * @param pThisGeometryDimension pointer to the dimensional data
    * @param ThisGeometryShapeFunctionContainer including the evaluated
    *        values for the shape functions, it's derivatives and the
    *        integration points.
    */
    GeometryData(GeometryDimension const *pThisGeometryDimension,
        const GeometryShapeFunctionContainer<IntegrationMethod>& ThisGeometryShapeFunctionContainer)
        : mpGeometryDimension(pThisGeometryDimension)
        , mGeometryShapeFunctionContainer(
                ThisGeometryShapeFunctionContainer)
    {
    }

    /// Copy constructor.
    GeometryData( const GeometryData& rOther )
        : mpGeometryDimension( rOther.mpGeometryDimension)
        , mGeometryShapeFunctionContainer( rOther.mGeometryShapeFunctionContainer)
    {
    }

    /// Destructor.
    virtual ~GeometryData() = default;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    GeometryData& operator=( const GeometryData& rOther )
    {
        mpGeometryDimension = rOther.mpGeometryDimension;
        mGeometryShapeFunctionContainer = rOther.mGeometryShapeFunctionContainer;

        return *this;
    }

    ///@}
    ///@name GeometryDimension
    ///@{

    void SetGeometryDimension(GeometryDimension const* pGeometryDimension)
    {
        mpGeometryDimension = pGeometryDimension;
    }

    ///@}
    ///@name Geometry Shape Function Container
    ///@{

    /// SetGeometryShapeFunctionContainer updates the GeometryShapeFunctionContainer.
    void SetGeometryShapeFunctionContainer(
        const GeometryShapeFunctionContainer<IntegrationMethod>& rGeometryShapeFunctionContainer)
    {
        mGeometryShapeFunctionContainer = rGeometryShapeFunctionContainer;
    }

    /// Returns the GeometryShapeFunctionContainer.
    const GeometryShapeFunctionContainer<IntegrationMethod>& GetGeometryShapeFunctionContainer() const
    {
        return mGeometryShapeFunctionContainer;
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
    SizeType WorkingSpaceDimension() const
    {
        return mpGeometryDimension->WorkingSpaceDimension();
    }

    /** Local space dimension. for example a triangle is a 2
    dimensional shape but can have 3 dimensional area
    coordinates l1, l2, l3.

    @return SizeType, local space dimension of this geometry.
    @see Dimension()
    @see WorkingSpaceDimension()
    */
    SizeType LocalSpaceDimension() const
    {
        return mpGeometryDimension->LocalSpaceDimension();
    }


    ///@}
    ///@name Inquiry
    ///@{

    /** This method confirm you if this geometry has a specific
    integration method or not. This method will be useful to
    control the geometry before integrating using a specific
    method. In GeometryData class this method controls if the
    integration points vector respecting to this method is empty
    or not.

    @return bool true if this integration method exist and false if this
    method is not imeplemented for this geometry.
    */
    bool HasIntegrationMethod( IntegrationMethod ThisMethod ) const
    {
        return mGeometryShapeFunctionContainer.HasIntegrationMethod(ThisMethod);
    }

    ///@}
    ///@name Integration
    ///@{

    /** Number of integration points for default integration
    method. This method just call IntegrationPointsNumber(enum
    IntegrationMethod ThisMethod) with default integration
    method.

    @return SizeType which is the number of integration points
    for default integrating method.
    */

    IntegrationMethod DefaultIntegrationMethod() const
    {
        return mGeometryShapeFunctionContainer.DefaultIntegrationMethod();
    }

    SizeType IntegrationPointsNumber() const
    {
        return mGeometryShapeFunctionContainer.IntegrationPointsNumber();
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
        return mGeometryShapeFunctionContainer.IntegrationPointsNumber(ThisMethod);
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
        return mGeometryShapeFunctionContainer.IntegrationPoints();
    }

    /** Integtation points for given integration
    method. This method use integration points data base to
    obtain integration points Vector respected to
    given method.

    @return const IntegrationPointsArrayType which is Vector of integration points
    for default integrating method.
    */
    const IntegrationPointsArrayType& IntegrationPoints(  IntegrationMethod ThisMethod ) const
    {
        return mGeometryShapeFunctionContainer.IntegrationPoints(ThisMethod);
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
        return mGeometryShapeFunctionContainer.ShapeFunctionsValues();
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
    const Matrix& ShapeFunctionsValues(  IntegrationMethod ThisMethod ) const
    {
        return mGeometryShapeFunctionContainer.ShapeFunctionsValues(
            ThisMethod);
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
        return mGeometryShapeFunctionContainer.ShapeFunctionValue(
            IntegrationPointIndex,
            ShapeFunctionIndex);
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
    double ShapeFunctionValue(
        IndexType IntegrationPointIndex,
        IndexType ShapeFunctionIndex,
        IntegrationMethod ThisMethod ) const
    {
        return mGeometryShapeFunctionContainer.ShapeFunctionValue(
            IntegrationPointIndex, ShapeFunctionIndex, ThisMethod );
    }

    /** This method gives all shape functions gradients evaluated in all
    integration points of default integration method. It just
    call ShapeFunctionsLocalGradients( IntegrationMethod ThisMethod)
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
        return mGeometryShapeFunctionContainer.ShapeFunctionsLocalGradients();
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
        return mGeometryShapeFunctionContainer.ShapeFunctionsLocalGradients(ThisMethod);
    }

    /** This method gives gradient of given shape function evaluated in
    given integration point of default integration method. It just
    call ShapeFunctionLocalGradient(IndexType IntegrationPointIndex,
    IndexType ShapeFunctionIndex,  IntegrationMethod
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
    const Matrix& ShapeFunctionLocalGradient( IndexType IntegrationPointIndex ) const
    {
        return mGeometryShapeFunctionContainer.ShapeFunctionLocalGradient(IntegrationPointIndex);
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
    const Matrix& ShapeFunctionLocalGradient( IndexType IntegrationPointIndex,  IntegrationMethod ThisMethod ) const
    {
        return mGeometryShapeFunctionContainer.ShapeFunctionLocalGradient(IntegrationPointIndex, ThisMethod);
    }

    const Matrix& ShapeFunctionLocalGradient( IndexType IntegrationPointIndex, IndexType ShapeFunctionIndex,  IntegrationMethod ThisMethod ) const
    {
        return mGeometryShapeFunctionContainer.ShapeFunctionLocalGradient(IntegrationPointIndex, ThisMethod);
    }

    /*
    * @brief access to the shape function derivatives.
    * @param DerivativeOrderIndex defines the wanted order of the derivative
    * @param IntegrationPointIndex the corresponding control point of this geometry
    * @return the shape function or derivative value related to the input parameters
    *         the matrix is structured: (derivative dN_de / dN_du , the corresponding node)
    */
    const Matrix& ShapeFunctionDerivatives(
        IndexType DerivativeOrderIndex,
        IndexType IntegrationPointIndex,
        IntegrationMethod ThisMethod) const
    {
        return mGeometryShapeFunctionContainer.ShapeFunctionDerivatives(
            DerivativeOrderIndex, IntegrationPointIndex, ThisMethod);
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "geometry data";
    }

    /// Print information about this object.
    virtual void PrintInfo( std::ostream& rOStream ) const
    {
        rOStream << "geometry data";
    }

    /// Print object's data.
    virtual void PrintData( std::ostream& rOStream ) const
    {
        rOStream << "    Working space dimension : " << mpGeometryDimension->WorkingSpaceDimension() << std::endl;
        rOStream << "    Local space dimension   : " << mpGeometryDimension->LocalSpaceDimension();
    }

    ///@}
private:
    ///@name Member Variables
    ///@{

    GeometryDimension const* mpGeometryDimension;

    GeometryShapeFunctionContainer<IntegrationMethod> mGeometryShapeFunctionContainer;

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    virtual void save( Serializer& rSerializer ) const
    {
        rSerializer.save("GeometryDimension", mpGeometryDimension);
        rSerializer.save("GeometryShapeFunctionContainer", mGeometryShapeFunctionContainer);
    }

    virtual void load( Serializer& rSerializer )
    {
        rSerializer.load("GeometryDimension", const_cast<GeometryDimension*>(mpGeometryDimension));
        rSerializer.load("GeometryShapeFunctionContainer", mGeometryShapeFunctionContainer);
    }

    // Private default constructor for serialization
    GeometryData()
    {
    }

    ///@}

}; // Class GeometryData

///@}
///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >> ( std::istream& rIStream,
                                   GeometryData& rThis );

/// output stream function
inline std::ostream& operator << ( std::ostream& rOStream,
                                   const GeometryData& rThis )
{
    rThis.PrintInfo( rOStream );
    rOStream << std::endl;
    rThis.PrintData( rOStream );

    return rOStream;
}

///@}

}  // namespace Kratos.
