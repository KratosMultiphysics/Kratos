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
//
//

#if !defined(KRATOS_GEOMETRY_DATA_H_INCLUDED )
#define  KRATOS_GEOMETRY_DATA_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/ublas_interface.h"
#include "integration/integration_point.h"


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
    integrating. These enums are used to refere to each
    different integration methods:

    - GI_GAUSS_1 gaussian integration with order 1.
    - GI_GAUSS_2 gaussian integration with order 2.
    - GI_GAUSS_3 gaussian integration with order 3.
    - GI_GAUSS_4 gaussian integration with order 4.
    - GI_GAUSS_5 gaussian integration with order 5.
    */
    enum IntegrationMethod {GI_GAUSS_1, GI_GAUSS_2, GI_GAUSS_3, GI_GAUSS_4, GI_GAUSS_5, GI_EXTENDED_GAUSS_1, GI_EXTENDED_GAUSS_2, GI_EXTENDED_GAUSS_3, GI_EXTENDED_GAUSS_4, GI_EXTENDED_GAUSS_5, NumberOfIntegrationMethods };

    enum KratosGeometryFamily
    {
        Kratos_NoElement,
        Kratos_Point,
        Kratos_Linear,
        Kratos_Triangle,
        Kratos_Quadrilateral,
        Kratos_Tetrahedra,
        Kratos_Hexahedra,
        Kratos_Prism,
        Kratos_generic_family
    };

    enum KratosGeometryType
    {
        Kratos_generic_type,
        Kratos_Hexahedra3D20,
        Kratos_Hexahedra3D27,
        Kratos_Hexahedra3D8,
        Kratos_Prism3D15,
        Kratos_Prism3D6,
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
        Kratos_Triangle3D3,
        Kratos_Triangle3D6,
        Kratos_Line2D2,
        Kratos_Line2D3,
        Kratos_Line3D2,
        Kratos_Line3D3,
        Kratos_Point2D,
        Kratos_Point3D,
        Kratos_Sphere3D1
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
    typedef std::size_t IndexType;

    /** This typed used to return size or dimension in
    geometry data. Dimension, WorkingSpaceDimension, PointsNumber and
    ... return this type as their results.
    */
    typedef std::size_t SizeType;

    /** This type used for representing an integration point in
    geometry data. This integration point is a point with an
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
    typedef std::array<IntegrationPointsArrayType, NumberOfIntegrationMethods> IntegrationPointsContainerType;

    /** A third order tensor used as shape functions' values
    continer.
    */
    typedef std::array<Matrix, NumberOfIntegrationMethods> ShapeFunctionsValuesContainerType;

    /** A fourth order tensor used as shape functions' local
    gradients container in geometry data.
    */
    typedef std::array<DenseVector<Matrix>, NumberOfIntegrationMethods> ShapeFunctionsLocalGradientsContainerType;

    /** A third order tensor to hold shape functions'
    gradients. ShapefunctionsLocalGradients function return this
    type as its result.
    */
    typedef DenseVector<Matrix> ShapeFunctionsGradientsType;

    typedef DenseVector<Matrix> ShapeFunctionsSecondDerivativesType;

    /**
     * fourth order tensor to hold the third order derivatives of the
     * shape functions
     */
    typedef DenseVector<DenseVector<Matrix> > ShapeFunctionsThirdDerivativesType;

    ///@}
    ///@name Life Cycle
    ///@{

    /** Complete argument constructor. This constructor gives a
    complete set of arguments to pass all the initial value of
    all the member variables of geometry class. Also it has
    default value for integration variables to make it usefull
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
    GeometryData( SizeType ThisDimension,
                  SizeType ThisWorkingSpaceDimension,
                  SizeType ThisLocalSpaceDimension,
                  enum IntegrationMethod ThisDefaultMethod,
                  const IntegrationPointsContainerType& ThisIntegrationPoints,
                  const ShapeFunctionsValuesContainerType& ThisShapeFunctionsValues,
                  const ShapeFunctionsLocalGradientsContainerType& ThisShapeFunctionsLocalGradients )
        : mDimension( ThisDimension )
        , mWorkingSpaceDimension( ThisWorkingSpaceDimension )
        , mLocalSpaceDimension( ThisLocalSpaceDimension )
        , mDefaultMethod( ThisDefaultMethod )
        , mIntegrationPoints( ThisIntegrationPoints )
        , mShapeFunctionsValues( ThisShapeFunctionsValues )
        , mShapeFunctionsLocalGradients( ThisShapeFunctionsLocalGradients )
    {
    }

    /** Copy constructor.
    Construct this geometry data as a copy of given geometry data.
    */
    GeometryData( const GeometryData& rOther )
        : mDimension( rOther.mDimension )
        , mWorkingSpaceDimension( rOther.mWorkingSpaceDimension )
        , mLocalSpaceDimension( rOther.mLocalSpaceDimension )
        , mDefaultMethod( rOther.mDefaultMethod )
        , mIntegrationPoints( rOther.mIntegrationPoints )
        , mShapeFunctionsValues( rOther.mShapeFunctionsValues )
        , mShapeFunctionsLocalGradients( rOther.mShapeFunctionsLocalGradients )
    {
    }



    /// Destructor. Do nothing!!!
    virtual ~GeometryData() {}


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
    GeometryData& operator=( const GeometryData& rOther )
    {
        mDimension = rOther.mDimension;
        mWorkingSpaceDimension = rOther.mWorkingSpaceDimension;
        mLocalSpaceDimension = rOther.mLocalSpaceDimension;
        mDefaultMethod = rOther.mDefaultMethod;
        mIntegrationPoints = rOther.mIntegrationPoints;
        mShapeFunctionsValues = rOther.mShapeFunctionsValues;
        mShapeFunctionsLocalGradients = rOther.mShapeFunctionsLocalGradients;

        return *this;
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
    SizeType Dimension() const
    {
        return mDimension;
    }

    /** Working space dimension. for example a triangle is a 2
    dimensional shape but can be used in 3 dimensional space.

    @return SizeType, working space dimension of this geometry.
    @see Dimension()
    @see LocalSpaceDimension()
    */
    SizeType WorkingSpaceDimension() const
    {
        return mWorkingSpaceDimension;
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
        return mLocalSpaceDimension;
    }


    ///@}
    ///@name Inquiry
    ///@{

    /** This method confirm you if this geometry has a specific
    integration method or not. This method will be usefull to
    control the geometry before intagrating using a specific
    method. In GeometryData class this method controls if the
    integration points vector respecting to this method is empty
    or not.

    @return bool true if this integration method exist and false if this
    method is not imeplemented for this geometry.
    */
    bool HasIntegrationMethod( enum IntegrationMethod ThisMethod ) const
    {
        return ( !mIntegrationPoints[ThisMethod].empty() );
    }

    ///@}
    ///@name Integration
    ///@{

    /** Number of integtation points for default integration
    method. This method just call IntegrationPointsNumber(enum
    IntegrationMethod ThisMethod) with default integration
    method.

    @return SizeType which is the number of integration points
    for default integrating method.
    */

    IntegrationMethod DefaultIntegrationMethod() const
    {
        return mDefaultMethod;
    }

    SizeType IntegrationPointsNumber() const
    {
        return mIntegrationPoints[mDefaultMethod].size();
    }

    /** Number of integtation points for given integration
    method. This method use integration points data base to
    obtain size of the integration points Vector respected to
    given method.

    @return SizeType which is the number of integration points
    for given integrating method.
    */
    SizeType IntegrationPointsNumber( enum IntegrationMethod ThisMethod ) const
    {
        return mIntegrationPoints[ThisMethod].size();
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
        return mIntegrationPoints[mDefaultMethod];
    }

    /** Integtation points for given integration
    method. This method use integration points data base to
    obtain integration points Vector respected to
    given method.

    @return const IntegrationPointsArrayType which is Vector of integration points
    for default integrating method.
    */
    const IntegrationPointsArrayType& IntegrationPoints( enum IntegrationMethod ThisMethod ) const
    {
        return mIntegrationPoints[ThisMethod];
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
        return mShapeFunctionsValues[mDefaultMethod];
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
    const Matrix& ShapeFunctionsValues( enum IntegrationMethod ThisMethod ) const
    {
        return mShapeFunctionsValues[ThisMethod];
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
        if ( mShapeFunctionsValues[mDefaultMethod].size1() <= IntegrationPointIndex )
            KRATOS_ERROR << "No existing integration point" << std::endl;

        if ( mShapeFunctionsValues[mDefaultMethod].size2() <= ShapeFunctionIndex )
            KRATOS_ERROR << "No existing shape function value" << std::endl;

        return mShapeFunctionsValues[mDefaultMethod]( IntegrationPointIndex, ShapeFunctionIndex );
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
    double ShapeFunctionValue( IndexType IntegrationPointIndex, IndexType ShapeFunctionIndex, enum IntegrationMethod ThisMethod ) const
    {
        if ( mShapeFunctionsValues[ThisMethod].size1() <= IntegrationPointIndex )
            KRATOS_ERROR << "No existing integration point" << std::endl;

        if ( mShapeFunctionsValues[ThisMethod].size2() <= ShapeFunctionIndex )
            KRATOS_ERROR << "No existing shape function value" << std::endl;

        return mShapeFunctionsValues[ThisMethod]( IntegrationPointIndex, ShapeFunctionIndex );
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
        return mShapeFunctionsLocalGradients[mDefaultMethod];
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
    const ShapeFunctionsGradientsType& ShapeFunctionsLocalGradients( enum IntegrationMethod ThisMethod ) const
    {
        return mShapeFunctionsLocalGradients[ThisMethod];
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
    const Matrix& ShapeFunctionLocalGradient( IndexType IntegrationPointIndex ) const
    {
        if ( mShapeFunctionsLocalGradients[mDefaultMethod].size() <= IntegrationPointIndex )
            KRATOS_ERROR << "No existing integration point" << std::endl;

        return mShapeFunctionsLocalGradients[mDefaultMethod][IntegrationPointIndex];
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
    const Matrix& ShapeFunctionLocalGradient( IndexType IntegrationPointIndex, enum IntegrationMethod ThisMethod ) const
    {
        if ( mShapeFunctionsLocalGradients[ThisMethod].size() <= IntegrationPointIndex )
        {
            KRATOS_ERROR << "No existing integration point" << std::endl;
        }

        return mShapeFunctionsLocalGradients[ThisMethod][IntegrationPointIndex];
    }

    const Matrix& ShapeFunctionLocalGradient( IndexType IntegrationPointIndex, IndexType ShapeFunctionIndex, enum IntegrationMethod ThisMethod ) const
    {
        if ( mShapeFunctionsLocalGradients[ThisMethod].size() <= IntegrationPointIndex )
            KRATOS_ERROR << "No existing integration point" << std::endl;

        return mShapeFunctionsLocalGradients[ThisMethod][IntegrationPointIndex];
    }


    ///@}
    ///@name Input and output
    ///@{

    /** Turn back information as a string.

    @return String contains information about this geometry.
    @see PrintData()
    @see PrintInfo()
    */
    virtual std::string Info() const
    {
        return "geometry data";
    }

    /** Print information about this object.

    @param rOStream Stream to print into it.
    @see PrintData()
    @see Info()
    */
    virtual void PrintInfo( std::ostream& rOStream ) const
    {
        rOStream << "geometry data";
    }

    /** Print geometry's data into given stream. Prints it's points
    by the order they stored in the geometry and then center
    point of geometry.

    @param rOStream Stream to print into it.
    @see PrintInfo()
    @see Info()
    */
    virtual void PrintData( std::ostream& rOStream ) const
    {
        rOStream << "    Dimension               : " << mDimension << std::endl;
        rOStream << "    working space dimension : " << mWorkingSpaceDimension << std::endl;
        rOStream << "    Local space dimension   : " << mLocalSpaceDimension;
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


    ///@}
    ///@name Protected  Access
    ///@{


    ///@}
    ///@name Protected Inquiry
    ///@{


    ///@}
    ///@name Protected LifeCycle
    ///@{


    ///@}

private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{

    SizeType mDimension;

    SizeType mWorkingSpaceDimension;

    SizeType mLocalSpaceDimension;

    enum IntegrationMethod mDefaultMethod;

    IntegrationPointsContainerType mIntegrationPoints;

    ShapeFunctionsValuesContainerType mShapeFunctionsValues;

    ShapeFunctionsLocalGradientsContainerType mShapeFunctionsLocalGradients;

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    virtual void save( Serializer& rSerializer ) const
    {
        rSerializer.save( "Dimension", mDimension );
        rSerializer.save( "Working Space Dimension", mWorkingSpaceDimension );
        rSerializer.save( "Local Space Dimension", mLocalSpaceDimension );
//    rSerializer.save("Default Method", mDefaultMethod);
//    rSerializer.save("Integration Points", mIntegrationPoints);
//    rSerializer.save("Shape Functions Values", mShapeFunctionsValues);
//    rSerializer.save("Shape Functions Local Gradients", mShapeFunctionsLocalGradients);
    }

    virtual void load( Serializer& rSerializer )
    {
        rSerializer.load( "Dimension", mDimension );
        rSerializer.load( "Working Space Dimension", mWorkingSpaceDimension );
        rSerializer.load( "Local Space Dimension", mLocalSpaceDimension );
//    rSerializer.load("Default Method", mDefaultMethod);
//    rSerializer.load("Integration Points", mIntegrationPoints);
//    rSerializer.load("Shape Functions Values", mShapeFunctionsValues);
//    rSerializer.load("Shape Functions Local Gradients", mShapeFunctionsLocalGradients);
    }

    // Private default constructor for serialization
    GeometryData()
    {
    }


    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{


    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Private Friends
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{


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

#endif // KRATOS_GEOMETRY_DATA_H_INCLUDED  defined


