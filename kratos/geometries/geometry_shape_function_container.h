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

    typedef std::array<Matrix, IntegrationMethod::NumberOfIntegrationMethods> ShapeFunctionsValuesContainerType;
    typedef std::array<DenseVector<Matrix>, IntegrationMethod::NumberOfIntegrationMethods> ShapeFunctionsLocalGradientsContainerType;

    typedef DenseVector<Matrix> ShapeFunctionsGradientsType;
    typedef DenseVector<Matrix> ShapeFunctionsSecondDerivativesType;
    typedef DenseVector<DenseVector<Matrix> > ShapeFunctionsThirdDerivativesType;

    ///@}
    ///@name Life Cycle
    ///@{

    GeometryShapeFunctionContainer( 
                  const ShapeFunctionsValuesContainerType& ThisShapeFunctionsValues,
                  const ShapeFunctionsLocalGradientsContainerType& ThisShapeFunctionsLocalGradients )
        : mShapeFunctionsValues( ThisShapeFunctionsValues )
        , mShapeFunctionsLocalGradients( ThisShapeFunctionsLocalGradients )
    {
    }

    GeometryShapeFunctionContainer()
        : mShapeFunctionsValues({})
        , mShapeFunctionsLocalGradients({})
    {
    }

    /** Copy constructor.
    Construct this geometry shape function container as a copy of given geometry data.
    */
    GeometryShapeFunctionContainer( const GeometryShapeFunctionContainer& rOther )
        : mShapeFunctionsValues( rOther.mShapeFunctionsValues )
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
        mShapeFunctionsValues = rOther.mShapeFunctionsValues;
        mShapeFunctionsLocalGradients = rOther.mShapeFunctionsLocalGradients;

        return *this;
    }

    ///@}
    ///@name Shape Function
    ///@{

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

    const ShapeFunctionsGradientsType& ShapeFunctionsLocalGradients(
        IntegrationMethod ThisMethod ) const
    {
        return mShapeFunctionsLocalGradients[ThisMethod];
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

    ShapeFunctionsValuesContainerType mShapeFunctionsValues;

    ShapeFunctionsLocalGradientsContainerType mShapeFunctionsLocalGradients;

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    virtual void save( Serializer& rSerializer ) const
    {
    }

    virtual void load( Serializer& rSerializer )
    {

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


