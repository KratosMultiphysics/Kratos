//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher
//

#if !defined(KRATOS_PROPERTIES_ACCESSOR_H)
#define  KRATOS_PROPERTIES_ACCESSOR_H

// System includes
#include <map>

// External includes

// Project includes
#include "includes/define.h"
#include "containers/data_value_container.h"
#include "geometries/geometry.h"
#include "includes/node.h"
#include "includes/process_info.h"
#include "includes/properties.h"


namespace Kratos
{
///@addtogroup KratosCore
///@{

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

/// Short class definition.
/** Detail class definition.
*/
class PropertiesAccessor
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of PropertiesAccessor
    KRATOS_CLASS_POINTER_DEFINITION(PropertiesAccessor);

    using NodeType = Node<3>;
    using GeometryType = Geometry<NodeType>;

    using IndexType = std::size_t;
    using SizeType = std::size_t;

    using DoubleVariableType = Variable<double>;
    using ComponentVariableType = VariableComponent< VectorComponentAdaptor<array_1d<double, 3> > >;
    using Array3VariableType = Variable<array_1d<double, 3>>;

    using DoubleFunction = std::function<double(const Properties*,
                                                const GeometryType&,
                                                const DataValueContainer&,
                                                const ProcessInfo&,
                                                const Vector&,
                                                const int)>;

    using DoubleConfiguration = std::map<DoubleVariableType, DoubleFunction>;

    // std::function<void(int)>

    // using std::unordered_map<DoubleVariableType>

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    PropertiesAccessor(Properties* pParentProperties)
        : mpParentProperties(pParentProperties) {}

    /// Destructor.
    virtual ~PropertiesAccessor() {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    template<class TVariableType>
    typename TVariableType::Type& GetValue(const TVariableType& rV,
                                           const GeometryType& rpGeometry,
                                           const DataValueContainer& rDataContainer,
                                           const ProcessInfo& rCurrentProcessInfo,
                                           const Vector& rShapeFunctionValues,
                                           const int GaussPointIndex=0)
    {
        mpParentProperties->GetValue(rV);
    }


    double GetValue(const DoubleVariableType& rV,
                     const GeometryType& rGeometry,
                     const DataValueContainer& rDataContainer,
                     const ProcessInfo& rCurrentProcessInfo,
                     const Vector& rShapeFunctionValues,
                     const int GaussPointIndex=0)
    {
        if(mConfigurationsDouble.find(rV) != mConfigurationsDouble.end())
            return mConfigurationsDouble[rV](mpParentProperties,
                                             rGeometry,
                                             rDataContainer,
                                             rCurrentProcessInfo,
                                             rShapeFunctionValues,
                                             GaussPointIndex);
        else
            return mpParentProperties->GetValue(rV);
    }

    void Configure(const DoubleVariableType& rV, std::function<double(
                                                const Properties*,
                                                const GeometryType&,
                                                const DataValueContainer&,
                                                const ProcessInfo&,
                                                const Vector&,
                                                const int)> DoubleFunction)
    {
        mConfigurationsDouble[rV] = DoubleFunction;
    }



    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "PropertiesAccessor" ;
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const {rOStream << "PropertiesAccessor";}

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const {}


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

    Properties* mpParentProperties;
    DoubleConfiguration mConfigurationsDouble;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const
    {
        // rSerializer.save("ParentProps", mpParentProperties);
        // TODO check if this works... (=> circular dependencies)
        // If it is not working I could make a "SetAccessorProps"
        // that is called in the load fct of the Properties and passes a "this"
    }

    void load(Serializer& rSerializer)
    {
        // rSerializer.load("ParentProps", mpParentProperties);
    }

    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    // /// Assignment operator.
    // PropertiesAccessor& operator=(PropertiesAccessor const& rOther){}

    // /// Copy constructor.
    // PropertiesAccessor(PropertiesAccessor const& rOther){}


    ///@}

}; // Class PropertiesAccessor

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


// /// input stream function
// inline std::istream& operator >> (std::istream& rIStream,
//                 PropertiesAccessor& rThis){}

// /// output stream function
// inline std::ostream& operator << (std::ostream& rOStream,
//                 const PropertiesAccessor& rThis)
// {
//     rThis.PrintInfo(rOStream);
//     rOStream << std::endl;
//     rThis.PrintData(rOStream);

//     return rOStream;
// }
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_PROPERTIES_ACCESSOR_H  defined
