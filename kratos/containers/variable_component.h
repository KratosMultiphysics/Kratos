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
//
//


#if !defined(KRATOS_VARIABLE_COMPONENT_H_INCLUDED )
#define  KRATOS_VARIABLE_COMPONENT_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "containers/variable_data.h"

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
 * @class VariableComponent
 * @brief Provide information for store or retrive a component of a variable in data container.
 * @details Provide information for store or retrive a component of a variable in data container. This class also provide a method to extract its component value from the source variable value in container.
 * @tparam TAdaptorType The adaptor variable type
 * @ingroup KratosCore
 * @author Pooyan Dadvand
 */
template<class TAdaptorType>
class VariableComponent : public VariableData
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of VariableComponent
    KRATOS_CLASS_POINTER_DEFINITION(VariableComponent);

    /// Type of this component.
    typedef typename TAdaptorType::Type Type;

    /// Data value type of this component.
    typedef typename TAdaptorType::Type DataType;

    /// Source value type of this component.
    typedef typename TAdaptorType::SourceType SourceType;

    /// Source variable type of this component.
    typedef Variable<typename TAdaptorType::SourceType> SourceVariableType;

    /// Base class type definition.
    typedef VariableData BaseType;

    /// Adaptor type.
    typedef TAdaptorType AdaptorType;

    typedef VariableComponent<TAdaptorType> VariableComponentType;

    ///@}
    ///@name Life Cycle
    ///@{

    VariableComponent(
        const std::string& rComponentName,
        const std::string& rSourceName,
        int ComponentIndex,
        const AdaptorType& rNewAdaptor,
        const VariableComponentType* pTimeDerivativeVariable = nullptr
        )
        : BaseType(rComponentName, sizeof(DataType),&rNewAdaptor.GetSourceVariable(), rNewAdaptor.GetComponentIndex()), mpSourceVariable(&rNewAdaptor.GetSourceVariable()),
          mpTimeDerivativeVariable(pTimeDerivativeVariable)
    {
        SetKey(GenerateKey(rSourceName, sizeof(DataType), true, ComponentIndex));
    }

    /// Copy constructor.
    VariableComponent(const VariableComponent& rOther)
        : BaseType(rOther), mpSourceVariable(rOther.mpSourceVariable) {}

    /// Destructor.
    ~VariableComponent() override {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    ///@}
    ///@name Access
    ///@{

    /**
     * @brief This method returns the time derivative component variable
     * @return The reference of the time derivative component variable (if any)
     */
    const VariableComponentType& GetTimeDerivative() const
    {
        KRATOS_DEBUG_ERROR_IF(mpTimeDerivativeVariable == nullptr) << "Time derivative for Variable \"" << Name() << "\" was not assigned" << std::endl;
        return *mpTimeDerivativeVariable;
    }

    const SourceVariableType& GetSourceVariable() const
    {
        return *mpSourceVariable;
    }

    DataType& GetValue(SourceType& SourceValue) const
    {
        return GetValueByIndex(SourceValue,GetComponentIndex());
    }

    const DataType& GetValue(const SourceType& SourceValue) const
    {
        return GetValueByIndex(SourceValue,GetComponentIndex());
    }

    static VariableComponent const& StaticObject()
    {
        static const VariableComponent<TAdaptorType> static_object("NONE", "NONE", 0, TAdaptorType::StaticObject());
        return static_object;
    }

    void Print(const void* pSource, std::ostream& rOStream) const override
    {
        rOStream << Name() << " component of " <<  GetSourceVariable().Name() << " variable : " <<  *static_cast<const DataType* >(pSource) ;
    }

    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << Name() << " component of " <<  GetSourceVariable().Name() << " variable";
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Name() << " component of " <<  GetSourceVariable().Name() << " variable";
    }

    /// Print object's data.
//       virtual void PrintData(std::ostream& rOStream) const;


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

    /// Assignment operator.
    VariableComponent& operator=(const VariableComponent& rOther)
    {
        BaseType::operator=(rOther);
        mpSourceVariable = rOther.mpSourceVariable;
    }

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

    static const VariableComponent  msStaticObject;

    ///@}
    ///@name Member Variables
    ///@{

    const SourceVariableType* mpSourceVariable;

    const VariableComponentType* mpTimeDerivativeVariable = nullptr; /// Definition of the pointer to the variable for the time derivative

    ///@}
    ///@name Serialization
    ///@{

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    /// This is the default function for getting a value by index using operator[]
    /** It is templated so one can create specialized version of this for types without operator[]
    **/
    template<typename TValueType>
    DataType& GetValueByIndex(TValueType& rValue, std::size_t index) const
    {
        return rValue[index];
    }

    template<typename TValueType>
    const DataType& GetValueByIndex(const TValueType& rValue, std::size_t index) const
    {
        return rValue[index];
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

    /// Default constructor.
    VariableComponent() {}

    ///@}

}; // Class VariableComponent

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template<class TAdaptorType>
inline std::istream& operator >> (std::istream& IStream,
                                  VariableComponent<TAdaptorType>& rThis);

/// output stream function
template<class TAdaptorType>
inline std::ostream& operator << (std::ostream& OStream,
                                  const VariableComponent<TAdaptorType>& rThis)
{
    rThis.PrintInfo(OStream);
    rThis.PrintData(OStream);

    return OStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_FILENAME_H_INCLUDED  defined


