//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Tobias Teschemacher
//                   Pooyan Dadvand
//
//

#if !defined(KRATOS_GEOMETRY_SHAPE_FUNCTION_CONTAINER_H_INCLUDED )
#define  KRATOS_GEOMETRY_SHAPE_FUNCTION_CONTAINER_H_INCLUDED

// System includes

// External includes

// Project includes

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

/** GeometryShapeFunctionContainer.
    Holds the evaluated values for the shape functions and its derivatives.
    Needs to be templated by the IntegrationMethod to import the
    IntegrationMethod enum which is defined in the GeometryData.
*/
template<typename TIntegrationMethodType>
class GeometryShapeFunctionContainer
{
public:

    ///@name Type Definitions
    ///@{

    /// Pointer definition of GeometryShapeFunctionContainer
    KRATOS_CLASS_POINTER_DEFINITION( GeometryShapeFunctionContainer );

    using IntegrationMethod = TIntegrationMethodType;

    typedef std::size_t IndexType;
    typedef std::size_t SizeType;

    typedef IntegrationPoint<3> IntegrationPointType;
    typedef std::vector<IntegrationPointType> IntegrationPointsArrayType;
    typedef std::array<IntegrationPointsArrayType, IntegrationMethod::NumberOfIntegrationMethods> IntegrationPointsContainerType;

    // Shape Function Values
    typedef Matrix ShapeFunctionsValuesType;
    typedef std::array<ShapeFunctionsValuesType, IntegrationMethod::NumberOfIntegrationMethods> ShapeFunctionsValuesContainerType;

    // Shape Function Gradients Values
    typedef DenseVector<Matrix> ShapeFunctionsGradientsType;
    typedef std::array<DenseVector<Matrix>, IntegrationMethod::NumberOfIntegrationMethods> ShapeFunctionsLocalGradientsContainerType;

    // Shape Function Derivatives Values
    typedef DenseVector<Matrix>
        ShapeFunctionsDerivativesType;
    typedef DenseVector<ShapeFunctionsDerivativesType>
        ShapeFunctionsDerivativesIntegrationPointsArrayType;
    typedef std::array<ShapeFunctionsDerivativesIntegrationPointsArrayType, IntegrationMethod::NumberOfIntegrationMethods>
        ShapeFunctionsDervativesContainerType;

    ///@}
    ///@name Life Cycle
    ///@{

    GeometryShapeFunctionContainer(
                  IntegrationMethod ThisDefaultMethod,
                  const IntegrationPointsContainerType& ThisIntegrationPoints,
                  const ShapeFunctionsValuesContainerType& ThisShapeFunctionsValues,
                  const ShapeFunctionsLocalGradientsContainerType& ThisShapeFunctionsLocalGradients )
        : mShapeFunctionsValues( ThisShapeFunctionsValues )
    {
        for (IndexType i = 0; i < ThisShapeFunctionsLocalGradients.size(); ++i)
        {
            ShapeFunctionsDerivativesIntegrationPointsArrayType integration_points(ThisShapeFunctionsLocalGradients[i].size());

            for (IndexType j = 0; j < IntegrationMethod::NumberOfIntegrationMethods; ++j)
            {
                ShapeFunctionsDerivativesType gradients(1);
                gradients[0] = ThisShapeFunctionsLocalGradients[i][j];

                integration_points[j] = gradients;
            }

            mShapeFunctionsDervativesContainer[i] = integration_points;
        }
    }

    GeometryShapeFunctionContainer(
        IntegrationMethod ThisDefaultMethod,
        const IntegrationPointsContainerType& ThisIntegrationPoints,
        const ShapeFunctionsValuesContainerType& ThisShapeFunctionsValues,
        const ShapeFunctionsDervativesContainerType& ThisShapeFunctionsDervativesContainer)
        : mShapeFunctionsValues(ThisShapeFunctionsValues)
        , mShapeFunctionsDervativesContainer(ThisShapeFunctionsDervativesContainer)
    {
    }

    GeometryShapeFunctionContainer(
        IntegrationMethod ThisIntegrationMethod,
        const IntegrationPointsArrayType& ThisIntegrationPoints,
        const ShapeFunctionsValuesType& ThisShapeFunctionsValues,
        const ShapeFunctionsDerivativesIntegrationPointsArrayType& ThisShapeFunctionsDerivativesIntegrationPointsArray)
        : mDefaultMethod(ThisIntegrationMethod)
    {
        mIntegrationPoints[ThisIntegrationMethod] = ThisIntegrationPoints;
        mShapeFunctionsValues[ThisIntegrationMethod] = ThisShapeFunctionsValues;
        mShapeFunctionsDervativesContainer[ThisIntegrationMethod] = ThisShapeFunctionsDerivativesIntegrationPointsArray;
    }

    GeometryShapeFunctionContainer()
        : mIntegrationPoints({})
        , mShapeFunctionsValues({})
        , mShapeFunctionsDervativesContainer({})
    {
    }

    /** Copy constructor.
    Construct this geometry shape function container as a copy of given geometry data.
    */
    GeometryShapeFunctionContainer( const GeometryShapeFunctionContainer& rOther )
        : mShapeFunctionsValues( rOther.mShapeFunctionsValues )
        , mShapeFunctionsDervativesContainer( rOther.mShapeFunctionsDervativesContainer)
    {
    }

    /// Destructor. Do nothing!!!
    virtual ~GeometryShapeFunctionContainer() {}

    ///@}
    ///@name Operators
    ///@{

    GeometryShapeFunctionContainer& operator=( 
        const GeometryShapeFunctionContainer& rOther )
    {
        mShapeFunctionsValues = rOther.mShapeFunctionsValues;
        mShapeFunctionsDervativesContainer = rOther.mShapeFunctionsDervativesContainer;

        return *this;
    }

    ///@}
    ///@name Integration
    ///@{

    inline bool HasIntegrationMethod(IntegrationMethod ThisMethod) const
    {
        return (!mIntegrationPoints[ThisMethod].empty());
    }

    inline IntegrationMethod DefaultIntegrationMethod() const
    {
        return mDefaultMethod;
    }

    inline SizeType IntegrationPointsNumber() const
    {
        return mIntegrationPoints[mDefaultMethod].size();
    }

    inline SizeType IntegrationPointsNumber(IntegrationMethod ThisMethod) const
    {
        return mIntegrationPoints[ThisMethod].size();
    }

    inline const IntegrationPointsArrayType& IntegrationPoints() const
    {
        return mIntegrationPoints[mDefaultMethod];
    }

    const IntegrationPointsArrayType& IntegrationPoints(IntegrationMethod ThisMethod) const
    {
        return mIntegrationPoints[ThisMethod];
    }

    ///@}
    ///@name Shape Function
    ///@{

    const Matrix& ShapeFunctionsValues() const
    {
        return ShapeFunctionsValues(mDefaultMethod);
    }

    const Matrix& ShapeFunctionsValues( 
        IntegrationMethod ThisMethod ) const
    {
        return mShapeFunctionsValues[ThisMethod];
    }

    double ShapeFunctionValue(IndexType IntegrationPointIndex, IndexType ShapeFunctionIndex) const
    {
        return ShapeFunctionValue(
            IntegrationPointIndex,
            ShapeFunctionIndex,
            mDefaultMethod);
    }

    double ShapeFunctionValue(
        IndexType IntegrationPointIndex,
        IndexType ShapeFunctionIndex,
        IntegrationMethod ThisMethod ) const
    {
        KRATOS_DEBUG_ERROR_IF(mShapeFunctionsValues[ThisMethod].size1() <= IntegrationPointIndex )
            << "Not existing integration point with index: " << IntegrationPointIndex << std::endl;

        KRATOS_DEBUG_ERROR_IF(mShapeFunctionsValues[ThisMethod].size2() <= ShapeFunctionIndex )
            << "No existing shape function value" << std::endl;

        return mShapeFunctionsValues[ThisMethod]( IntegrationPointIndex, ShapeFunctionIndex );
    }

    const ShapeFunctionsGradientsType& ShapeFunctionsLocalGradients() const
    {
        return ShapeFunctionsLocalGradients(mDefaultMethod);
    }

    const ShapeFunctionsGradientsType& ShapeFunctionsLocalGradients(
        IntegrationMethod ThisMethod ) const
    {
        KRATOS_ERROR << "const ShapeFunctionsGradientsType& ShapeFunctionsLocalGradients(IntegrationMethod ThisMethod ) const does not exist anymore." <<
            "Shape functions can only be accessed directly by the integration point." << std::endl;
        //return mShapeFunctionsLocalGradients[ThisMethod];
    }

    const Matrix& ShapeFunctionLocalGradient(IndexType IntegrationPointIndex) const
    {
        return ShapeFunctionLocalGradient(IntegrationPointIndex, mDefaultMethod);
    }

    const Matrix& ShapeFunctionLocalGradient(
        IndexType IntegrationPointIndex,
        IntegrationMethod ThisMethod ) const
    {
        KRATOS_DEBUG_ERROR_IF(mShapeFunctionsDervativesContainer[ThisMethod].size() <= IntegrationPointIndex )
            << "Not existing integration point with index: " << IntegrationPointIndex << std::endl;

        return mShapeFunctionsDervativesContainer[ThisMethod][IntegrationPointIndex][0];
    }

    const ShapeFunctionsDerivativesType& ShapeFunctionDerivatives(
        IndexType IntegrationPointIndex,
        IndexType DerivativeOrder,
        IntegrationMethod ThisMethod) const
    {
        KRATOS_DEBUG_ERROR_IF(mShapeFunctionsDervativesContainer[ThisMethod].size() < IntegrationPointIndex)
            << "Not enough integration points for integration point index: " << IntegrationPointIndex << std::endl;
        KRATOS_DEBUG_ERROR_IF(mShapeFunctionsDervativesContainer[ThisMethod][IntegrationPointIndex].size() < DerivativeOrder)
            << "Not enough derivatives for derivative order: " << DerivativeOrder << std::endl;

        return mShapeFunctionsDervativesContainer[ThisMethod][IntegrationPointIndex][DerivativeOrder];
    }

    const Matrix& ShapeFunctionDerivativesAll(
        IndexType IntegrationPointIndex,
        IntegrationMethod ThisMethod) const
    {
        KRATOS_DEBUG_ERROR_IF(mShapeFunctionsDervativesContainer[ThisMethod].size() < IntegrationPointIndex)
            << "Not enough integration points for integration point index: " << IntegrationPointIndex << std::endl;

        return mShapeFunctionsDervativesContainer[ThisMethod][IntegrationPointIndex];
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
        rOStream << "shape function container";
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
    }

    ///@}

protected:

private:
    ///@name Member Variables
    ///@{

    IntegrationMethod mDefaultMethod;

    IntegrationPointsContainerType mIntegrationPoints;

    ShapeFunctionsValuesContainerType mShapeFunctionsValues;

    //ShapeFunctionsLocalGradientsContainerType mShapeFunctionsLocalGradients;

    ShapeFunctionsDervativesContainerType mShapeFunctionsDervativesContainer;

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    virtual void save( Serializer& rSerializer ) const
    {
        //rSerializer.save("DefaultMethod", mDefaultMethod);
        //rSerializer.save("IntegrationPoints", mIntegrationPoints);
        rSerializer.save("ShapeFunctionsValues", mShapeFunctionsValues);
        //rSerializer.save("ShapeFunctionsLocalGradients", mShapeFunctionsLocalGradients);
        rSerializer.save("ShapeFunctionsDervativesContainer", mShapeFunctionsDervativesContainer);
    }

    virtual void load( Serializer& rSerializer )
    {
        //rSerializer.load("DefaultMethod", mDefaultMethod);
        //rSerializer.load("IntegrationPoints", mIntegrationPoints);
        rSerializer.load("ShapeFunctionsValues", mShapeFunctionsValues);
        //rSerializer.load("ShapeFunctionsLocalGradients", mShapeFunctionsLocalGradients);
        rSerializer.load("ShapeFunctionsDervativesContainer", mShapeFunctionsDervativesContainer);
    }

    ///@}

}; // Class GeometryShapeFunctionContainer

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template<typename TIntegrationMethodType>
inline std::istream& operator >> ( std::istream& rIStream,
                                   GeometryShapeFunctionContainer<TIntegrationMethodType>& rThis );

/// output stream function
template<typename TIntegrationMethodType>
inline std::ostream& operator << ( std::ostream& rOStream,
                                   const GeometryShapeFunctionContainer<TIntegrationMethodType>& rThis )
{
    rThis.PrintInfo( rOStream );
    rOStream << std::endl;
    rThis.PrintData( rOStream );

    return rOStream;
}

///@}


}  // namespace Kratos.

#endif // KRATOS_GEOMETRY_SHAPE_FUNCTION_CONTAINER_H_INCLUDED  defined


