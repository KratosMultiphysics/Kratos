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

    /// Defining enum of integration method
    using IntegrationMethod = TIntegrationMethodType;

    /// Size types
    typedef std::size_t IndexType;
    typedef std::size_t SizeType;

    /// Integration points
    typedef IntegrationPoint<3> IntegrationPointType;
    typedef std::vector<IntegrationPointType> IntegrationPointsArrayType;
    typedef std::array<IntegrationPointsArrayType, IntegrationMethod::NumberOfIntegrationMethods> IntegrationPointsContainerType;

    /// Shape functions
    typedef std::array<Matrix, IntegrationMethod::NumberOfIntegrationMethods> ShapeFunctionsValuesContainerType;

    /// First derivatives/ gradients
    typedef DenseVector<Matrix> ShapeFunctionsGradientsType;
    typedef std::array<DenseVector<Matrix>, IntegrationMethod::NumberOfIntegrationMethods> ShapeFunctionsLocalGradientsContainerType;

    /// Higher order derivatives
    typedef DenseVector<Matrix>
        ShapeFunctionsDerivativesType;
    typedef DenseVector<ShapeFunctionsDerivativesType >
        ShapeFunctionsDerivativesIntegrationPointArrayType;
    typedef std::array<ShapeFunctionsDerivativesIntegrationPointArrayType, IntegrationMethod::NumberOfIntegrationMethods>
        ShapeFunctionsDerivativesContainerType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor for single integration point having the full containers
    GeometryShapeFunctionContainer(
        IntegrationMethod ThisDefaultMethod,
        const IntegrationPointsContainerType& ThisIntegrationPoints,
        const ShapeFunctionsValuesContainerType& ThisShapeFunctionsValues,
        const ShapeFunctionsLocalGradientsContainerType& ThisShapeFunctionsLocalGradients )
        : mDefaultMethod(ThisDefaultMethod)
        , mIntegrationPoints(ThisIntegrationPoints)
        , mShapeFunctionsValues( ThisShapeFunctionsValues )
        , mShapeFunctionsLocalGradients( ThisShapeFunctionsLocalGradients )
    {
    }

    /// Constructor ONLY for single integration point with first derivatives
    GeometryShapeFunctionContainer(
        IntegrationMethod ThisDefaultMethod,
        const IntegrationPointType& ThisIntegrationPoint,
        const Matrix& ThisShapeFunctionsValues,
        const Matrix& ThisShapeFunctionsGradients)
        : mDefaultMethod(ThisDefaultMethod)
    {
        IntegrationPointsArrayType ips(1);
        ips[0] = ThisIntegrationPoint;
        mIntegrationPoints[ThisDefaultMethod] = ips;

        mShapeFunctionsValues[ThisDefaultMethod] = ThisShapeFunctionsValues;

        ShapeFunctionsGradientsType DN_De_array(1);
        DN_De_array[0] = ThisShapeFunctionsGradients;
        mShapeFunctionsLocalGradients[ThisDefaultMethod] = DN_De_array;
    }

    /// Constructor ONLY for single integration point with multiple derivatives
    GeometryShapeFunctionContainer(
        IntegrationMethod ThisDefaultMethod,
        const IntegrationPointType& ThisIntegrationPoint,
        const Matrix& ThisShapeFunctionsValues,
        const DenseVector<Matrix>& ThisShapeFunctionsDerivatives)
        : mDefaultMethod(ThisDefaultMethod)
    {
        IntegrationPointsArrayType ips(1);
        ips[0] = ThisIntegrationPoint;
        mIntegrationPoints[ThisDefaultMethod] = ips;

        mShapeFunctionsValues[ThisDefaultMethod] = ThisShapeFunctionsValues;

        if (ThisShapeFunctionsDerivatives.size() > 0)
        {
            ShapeFunctionsGradientsType DN_De_array(1);
            DN_De_array[0] = ThisShapeFunctionsDerivatives[0];
            mShapeFunctionsLocalGradients[ThisDefaultMethod] = DN_De_array;
        }
        if (ThisShapeFunctionsDerivatives.size() > 1)
        {
            ShapeFunctionsDerivativesIntegrationPointArrayType derivatives_array(ThisShapeFunctionsDerivatives.size() - 1);
            for (int i = 1; i < ThisShapeFunctionsDerivatives.size(); ++i)
            {
                ShapeFunctionsDerivativesType DN_De_i_array(1);
                DN_De_i_array[0] = ThisShapeFunctionsDerivatives[i];
                derivatives_array[i - 1] = DN_De_i_array;
            }
            mShapeFunctionsDerivatives[ThisDefaultMethod] = derivatives_array;
        }
    }

    /// Default Constructor
    GeometryShapeFunctionContainer()
        : mIntegrationPoints({})
        , mShapeFunctionsValues({})
        , mShapeFunctionsLocalGradients({})
    {
    }

    /*
    * Copy constructor.
    * Construct this geometry shape function container as a copy of given geometry data.
    */
    GeometryShapeFunctionContainer( const GeometryShapeFunctionContainer& rOther )
        : mDefaultMethod(rOther.mDefaultMethod)
        , mIntegrationPoints(rOther.mIntegrationPoints)
        , mShapeFunctionsValues( rOther.mShapeFunctionsValues )
        , mShapeFunctionsLocalGradients( rOther.mShapeFunctionsLocalGradients )
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
        mDefaultMethod = rOther.mDefaultMethod;
        mIntegrationPoints = rOther.mIntegrationPoints;
        mShapeFunctionsValues = rOther.mShapeFunctionsValues;
        mShapeFunctionsLocalGradients = rOther.mShapeFunctionsLocalGradients;

        return *this;
    }

    ///@}
    ///@name Integration Method
    ///@{

    IntegrationMethod DefaultIntegrationMethod() const
    {
        return mDefaultMethod;
    }

    bool HasIntegrationMethod(IntegrationMethod ThisMethod) const
    {
        return (!mIntegrationPoints[ThisMethod].empty());
    }

    ///@}
    ///@name Integration Points
    ///@{

    SizeType IntegrationPointsNumber() const
    {
        return mIntegrationPoints[mDefaultMethod].size();
    }

    SizeType IntegrationPointsNumber(IntegrationMethod ThisMethod) const
    {
        return mIntegrationPoints[ThisMethod].size();
    }

    const IntegrationPointsArrayType& IntegrationPoints() const
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
        return mShapeFunctionsValues[mDefaultMethod];
    }

    const Matrix& ShapeFunctionsValues( 
        IntegrationMethod ThisMethod ) const
    {
        return mShapeFunctionsValues[ThisMethod];
    }

    double ShapeFunctionValue(
        IndexType IntegrationPointIndex,
        IndexType ShapeFunctionIndex,
        IntegrationMethod ThisMethod ) const
    {
        KRATOS_DEBUG_ERROR_IF(mShapeFunctionsValues[ThisMethod].size1() <= IntegrationPointIndex )
            << "No existing integration point" << std::endl;

        KRATOS_DEBUG_ERROR_IF(mShapeFunctionsValues[ThisMethod].size2() <= ShapeFunctionIndex )
            << "No existing shape function value" << std::endl;

        return mShapeFunctionsValues[ThisMethod]( IntegrationPointIndex, ShapeFunctionIndex );
    }

    double ShapeFunctionValue(IndexType IntegrationPointIndex, IndexType ShapeFunctionIndex) const
    {
        return mShapeFunctionsValues[mDefaultMethod](IntegrationPointIndex, ShapeFunctionIndex);
    }

    const ShapeFunctionsGradientsType& ShapeFunctionsLocalGradients() const
    {
        return mShapeFunctionsLocalGradients[mDefaultMethod];
    }

    const ShapeFunctionsGradientsType& ShapeFunctionsLocalGradients(
        IntegrationMethod ThisMethod ) const
    {
        return mShapeFunctionsLocalGradients[ThisMethod];
    }

    const Matrix& ShapeFunctionLocalGradient(
        IndexType IntegrationPointIndex) const
    {
        return ShapeFunctionLocalGradient(IntegrationPointIndex, mDefaultMethod);
    }

    const Matrix& ShapeFunctionLocalGradient(
        IndexType IntegrationPointIndex,
        IntegrationMethod ThisMethod ) const
    {
        KRATOS_DEBUG_ERROR_IF(mShapeFunctionsLocalGradients[ThisMethod].size() <= IntegrationPointIndex )
            << "No existing integration point" << std::endl;

        return mShapeFunctionsLocalGradients[ThisMethod][IntegrationPointIndex];
    }

    const Matrix& ShapeFunctionLocalGradient(
        IndexType IntegrationPointIndex,
        IndexType ShapeFunctionIndex,
        IntegrationMethod ThisMethod ) const
    {
        KRATOS_DEBUG_ERROR_IF(mShapeFunctionsLocalGradients[ThisMethod].size() <= IntegrationPointIndex )
            << "No existing integration point" << std::endl;

        return mShapeFunctionsLocalGradients[ThisMethod][IntegrationPointIndex];
    }

    /*
    * @brief access to the shape function derivatives.
    * @param DerivativeOrderIndex defines the wanted order of the derivative
    * @param IntegrationPointIndex the corresponding contorl point of this geometry
    * @return the shape function or derivative value related to the input parameters
    *         the matrix is structured: (derivative dN_de / dN_du , the corresponding node)
    */
    const Matrix& ShapeFunctionDerivatives(
        IndexType IntegrationPointIndex,
        IndexType DerivativeOrderIndex,
        IntegrationMethod ThisMethod) const
    {
        if (DerivativeOrderIndex > 1)
        {
            return mShapeFunctionsDerivatives[ThisMethod][IntegrationPointIndex][DerivativeOrderIndex - 2];
        }
        if (DerivativeOrderIndex == 1)
        {
            return mShapeFunctionsLocalGradients[ThisMethod][IntegrationPointIndex];
        }

        /* Shape function values are stored within a Matrix, however, only one row
           should be provided here. Thus, currently it is not possible to provide the 
           needed source to this object.*/
        KRATOS_ERROR_IF(DerivativeOrderIndex == 0)
            << "Shape functions cannot be accessed through ShapeFunctionDerivatives()" << std::endl;
    }

    /*
    * @brief access each item separateley.
    * @param DerivativeOrderIndex defines the wanted order of the derivative
    * @param DerivativeOrderRowIndex within each derivative the entries can
    *        be accessed differently.
    * @return the shape function or derivative value related to the input parameters.
    */
    double& ShapeFunctionDerivativeValue(
        IndexType IntegrationPointIndex,
        IndexType DerivativeOrderIndex,
        IndexType DerivativeOrderRowIndex,
        IndexType ShapeFunctionIndex,
        IntegrationMethod ThisMethod)
    {
        if (DerivativeOrderIndex == 0)
            return mShapeFunctionsValues[ThisMethod](IntegrationPointIndex, ShapeFunctionIndex);
        if (DerivativeOrderIndex == 1)
            return mShapeFunctionsLocalGradients[ThisMethod][IntegrationPointIndex](DerivativeOrderRowIndex, ShapeFunctionIndex);
        if (DerivativeOrderIndex > 1)
            return mShapeFunctionsDerivatives[ThisMethod][IntegrationPointIndex][DerivativeOrderIndex - 2](DerivativeOrderRowIndex, ShapeFunctionIndex);
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "shape function container";
    }

    /// Print information about this object.
    virtual void PrintInfo( std::ostream& rOStream ) const
    {
        rOStream << "shape function container";
    }

    /// Print object's data.
    virtual void PrintData( std::ostream& rOStream ) const
    {
    }

    ///@}

private:
    ///@name Member Variables
    ///@{

    IntegrationMethod mDefaultMethod;

    IntegrationPointsContainerType mIntegrationPoints;

    ShapeFunctionsValuesContainerType mShapeFunctionsValues;

    ShapeFunctionsLocalGradientsContainerType mShapeFunctionsLocalGradients;

    ShapeFunctionsDerivativesContainerType mShapeFunctionsDerivatives;

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    virtual void save( Serializer& rSerializer ) const
    {
        rSerializer.save("IntegrationPoints", mIntegrationPoints);
        rSerializer.save("ShapeFunctionsValues", mShapeFunctionsValues);
        rSerializer.save("ShapeFunctionsLocalGradients", mShapeFunctionsLocalGradients);
        rSerializer.save("ShapeFunctionsDerivatives", mShapeFunctionsDerivatives);
    }

    virtual void load( Serializer& rSerializer )
    {
        rSerializer.load("IntegrationPoints", mIntegrationPoints);
        rSerializer.load("ShapeFunctionsValues", mShapeFunctionsValues);
        rSerializer.load("ShapeFunctionsLocalGradients", mShapeFunctionsLocalGradients);
        rSerializer.load("ShapeFunctionsDerivatives", mShapeFunctionsDerivatives);
    }

    ///@}

}; // Class GeometryShapeFunctionContainer

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


