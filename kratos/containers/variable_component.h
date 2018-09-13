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
#include "variable_data.h"


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

/// Provide information for store or retrive a component of a variable in data container.
/** Provide information for store or retrive a component of a
    variable in data container. This class also provide a method to
    extract its component value from the source variable value in
    container.
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

    ///@}
    ///@name Life Cycle
    ///@{

    VariableComponent(const std::string& ComponentName, const std::string& SourceName, int ComponentIndex, const AdaptorType& NewAdaptor)
        : BaseType(ComponentName, sizeof(DataType), true, NewAdaptor.GetComponentIndex()), mAdaptor(NewAdaptor)
    {
        SetKey(GenerateKey(SourceName, sizeof(DataType), true,  ComponentIndex));
    }

    /// Copy constructor.
    VariableComponent(const VariableComponent& rOther)
        : BaseType(rOther), mAdaptor(rOther.mAdaptor) {}

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

    const SourceVariableType& GetSourceVariable() const
    {
        return mAdaptor.GetSourceVariable();
    }

    const AdaptorType& GetAdaptor() const
    {
        return mAdaptor;
    }

    DataType& GetValue(SourceType& SourceValue) const
    {
        return mAdaptor.GetValue(SourceValue);
    }

    const DataType& GetValue(const SourceType& SourceValue) const
    {
        return mAdaptor.GetValue(SourceValue);
    }

    static VariableComponent const& StaticObject()
    {
        return msStaticObject;
    }

    void Print(const void* pSource, std::ostream& rOStream) const override
    {
        rOStream << Name() << " component of " <<  mAdaptor.GetSourceVariable().Name() << " variable : " <<  *static_cast<const DataType* >(pSource) ;
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
        buffer << Name() << " component of " <<  mAdaptor.GetSourceVariable().Name() << " variable";
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Name() << " component of " <<  mAdaptor.GetSourceVariable().Name() << " variable";
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
        mAdaptor = rOther.mAdaptor;
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

    TAdaptorType mAdaptor;

    ///@}
    ///@name Serialization
    ///@{

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
    ///@name Un accessible methods
    ///@{

    /// Default constructor.
    VariableComponent() {}

    ///@}

}; // Class VariableComponent

///@}

template<class TAdaptorType>
const VariableComponent<TAdaptorType> VariableComponent<TAdaptorType>::msStaticObject("NONE", "NONE", 0, TAdaptorType::StaticObject());

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


