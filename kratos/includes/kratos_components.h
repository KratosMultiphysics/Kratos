//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//                   Riccardo Rossi
//

#if !defined(KRATOS_KRATOS_COMPONENTS_H_INCLUDED )
#define  KRATOS_KRATOS_COMPONENTS_H_INCLUDED

// System includes
#include <string>
#include <iostream>
#include <map>

// External includes

// Project includes
#include "includes/define.h"
#include "containers/variable_component.h"
#include "containers/vector_component_adaptor.h"
#include "containers/flags.h"
#include "utilities/quaternion.h"

namespace Kratos
{

/// KratosComponents class encapsulates a lookup table for a family of classes in a generic way.
/** KratosComponents class encapsulates a lookup table for a family of classes in a generic way.
 * Prototypes must be added to this table by unique names to be accessible by IO.
 * These names can be created automatically using C++ RTTI or given manually for each component.
 * In this design the manual approach is chosen, so shorter and more clear names can be given
 * to each component and also there is a flexibility to give different names to different
 * states of an object and create them via different prototypes.
 * For example having TriangularThermal and  both
*/

template<class TComponentType>
class KRATOS_API(KRATOS_CORE) KratosComponents
{

public:

    ///@name Type Definitions
    ///@{

    /// Pointer definition of KratosComponents
    KRATOS_CLASS_POINTER_DEFINITION(KratosComponents);

    typedef std::map<std::string, const TComponentType* > ComponentsContainerType;
    typedef typename ComponentsContainerType::value_type ValueType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.

    KratosComponents() {}

    /// Destructor.

    virtual ~KratosComponents() {}

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    static void Add(std::string const& Name, TComponentType const& ThisComponent)
    {
        msComponents.insert(ValueType(Name , &ThisComponent));
    }

    static TComponentType const& Get(std::string const& Name)
    {
        auto it_comp =  msComponents.find(Name);
        if (it_comp != msComponents.end()) {
            return *(it_comp->second);
        }

        std::stringstream err_msg;
        err_msg << "The component \"" << Name << "\" is not registered!\n"
                << "Maybe you need to import the application where it is defined?\n"
                << "The following components of this type are registered:" << std::endl;
        PrintData(err_msg);
        KRATOS_ERROR << err_msg.str() << std::endl;
    }

    static ComponentsContainerType & GetComponents()
    {
        return msComponents;
    }

    static ComponentsContainerType * pGetComponents()
    {
        return &msComponents;
    }

    static void Register()
    {

    }

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{

    static bool Has(std::string const& Name)
    {
        return (msComponents.find(Name) != msComponents.end());
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "Kratos components";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "Kratos components";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
        for (const auto& r_comp : msComponents) {
            rOStream << "    " << r_comp.first << std::endl;
        }
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

    static ComponentsContainerType msComponents;

    ///@}
    ///@name Member Variables
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

    /// Assignment operator.
    KratosComponents& operator=(KratosComponents const& rOther);

    /// Copy constructor.
    KratosComponents(KratosComponents const& rOther);

    ///@}

}; // Class KratosComponents


///@}
template<>
class KRATOS_API(KRATOS_CORE) KratosComponents<VariableData>
{
public:

    ///@name Type Definitions
    ///@{

    /// Pointer definition of KratosComponents
    KRATOS_CLASS_POINTER_DEFINITION(KratosComponents);

    typedef std::map<std::string, VariableData* > ComponentsContainerType;
    typedef ComponentsContainerType::value_type ValueType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    KratosComponents() {}

    /// Destructor.
    virtual ~KratosComponents() {}

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    static void Add(std::string const& Name, VariableData& ThisComponent)
    {
        msComponents.insert(ValueType(Name ,&ThisComponent));
    }

    static std::size_t Size()
    {
        return msComponents.size();
    }

    static VariableData & Get(std::string const& Name)
    {
        return *(msComponents.find(Name)->second);
    }

    static VariableData* pGet(std::string const& Name)
    {
        return (msComponents.find(Name)->second);
    }

    static ComponentsContainerType & GetComponents()
    {
        return msComponents;
    }

    static ComponentsContainerType * pGetComponents()
    {
        return &msComponents;
    }

    static void Register()
    {

    }

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{

    static bool Has(std::string const& Name)
    {
        return (msComponents.find(Name) != msComponents.end());
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "Kratos components";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "Kratos components";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
        for (const auto& r_comp : msComponents) {
            rOStream << "    " << r_comp.second << std::endl;
        }
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

    static ComponentsContainerType msComponents;

    ///@}
    ///@name Member Variables
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

    /// Assignment operator.
    KratosComponents& operator=(KratosComponents const& rOther);

    /// Copy constructor.
    KratosComponents(KratosComponents const& rOther);

    ///@}

}; // Class KratosComponents


template class KRATOS_API(KRATOS_CORE) KratosComponents<Variable<bool> >;
template class KRATOS_API(KRATOS_CORE) KratosComponents<Variable<int> >;
template class KRATOS_API(KRATOS_CORE) KratosComponents<Variable<unsigned int> >;
template class KRATOS_API(KRATOS_CORE) KratosComponents<Variable<double> >;
template class KRATOS_API(KRATOS_CORE) KratosComponents<Variable<array_1d<double, 3> > >;
template class KRATOS_API(KRATOS_CORE) KratosComponents<Variable<array_1d<double, 4> > >;
template class KRATOS_API(KRATOS_CORE) KratosComponents<Variable<array_1d<double, 6> > >;
template class KRATOS_API(KRATOS_CORE) KratosComponents<Variable<array_1d<double, 9> > >;
template class KRATOS_API(KRATOS_CORE) KratosComponents<Variable<Quaternion<double> > >;
template class KRATOS_API(KRATOS_CORE) KratosComponents<Variable<Vector> >;
template class KRATOS_API(KRATOS_CORE) KratosComponents<Variable<Matrix> >;
template class KRATOS_API(KRATOS_CORE) KratosComponents<Variable<std::string> >;
template class KRATOS_API(KRATOS_CORE) KratosComponents<VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > > >;
template class KRATOS_API(KRATOS_CORE) KratosComponents<VariableComponent<VectorComponentAdaptor<array_1d<double, 4> > > >;
template class KRATOS_API(KRATOS_CORE) KratosComponents<VariableComponent<VectorComponentAdaptor<array_1d<double, 6> > > >;
template class KRATOS_API(KRATOS_CORE) KratosComponents<VariableComponent<VectorComponentAdaptor<array_1d<double, 9> > > >;
template class KRATOS_API(KRATOS_CORE) KratosComponents<Variable<Flags> >;
template class KRATOS_API(KRATOS_CORE) KratosComponents<Flags>;

#ifdef KratosCore_EXPORTS
template<class TComponentType>
typename KratosComponents<TComponentType>::ComponentsContainerType KratosComponents<TComponentType>::msComponents;
#endif

///@name Input and output
///@{

/// output stream function
template<class TComponentType>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const KratosComponents<TComponentType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

void KRATOS_API(KRATOS_CORE) AddKratosComponent(std::string const& Name, Variable<bool> const& ThisComponent);
void KRATOS_API(KRATOS_CORE) AddKratosComponent(std::string const& Name, Variable<int> const& ThisComponent);
void KRATOS_API(KRATOS_CORE) AddKratosComponent(std::string const& Name, Variable<unsigned int> const& ThisComponent);
void KRATOS_API(KRATOS_CORE) AddKratosComponent(std::string const& Name, Variable<double> const& ThisComponent);
void KRATOS_API(KRATOS_CORE) AddKratosComponent(std::string const& Name, Variable<array_1d<double, 3> > const& ThisComponent);
void KRATOS_API(KRATOS_CORE) AddKratosComponent(std::string const& Name, Variable<array_1d<double, 4> > const& ThisComponent);
void KRATOS_API(KRATOS_CORE) AddKratosComponent(std::string const& Name, Variable<array_1d<double, 6> > const& ThisComponent);
void KRATOS_API(KRATOS_CORE) AddKratosComponent(std::string const& Name, Variable<array_1d<double, 9> > const& ThisComponent);
void KRATOS_API(KRATOS_CORE) AddKratosComponent(std::string const& Name, Variable<Quaternion<double> > const& ThisComponent);
void KRATOS_API(KRATOS_CORE) AddKratosComponent(std::string const& Name, Variable<Vector> const& ThisComponent);
void KRATOS_API(KRATOS_CORE) AddKratosComponent(std::string const& Name, Variable<Matrix> const& ThisComponent);
void KRATOS_API(KRATOS_CORE) AddKratosComponent(std::string const& Name, Variable<std::string> const& ThisComponent);
void KRATOS_API(KRATOS_CORE) AddKratosComponent(std::string const& Name, VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > > const& ThisComponent);
void KRATOS_API(KRATOS_CORE) AddKratosComponent(std::string const& Name, VariableComponent<VectorComponentAdaptor<array_1d<double, 4> > > const& ThisComponent);
void KRATOS_API(KRATOS_CORE) AddKratosComponent(std::string const& Name, VariableComponent<VectorComponentAdaptor<array_1d<double, 6> > > const& ThisComponent);
void KRATOS_API(KRATOS_CORE) AddKratosComponent(std::string const& Name, VariableComponent<VectorComponentAdaptor<array_1d<double, 9> > > const& ThisComponent);
void KRATOS_API(KRATOS_CORE) AddKratosComponent(std::string const& Name, Flags const& ThisComponent);
void KRATOS_API(KRATOS_CORE) AddKratosComponent(std::string const& Name, Variable<Flags> const& ThisComponent);

template<class TComponentType> void AddKratosComponent(std::string const& Name, TComponentType const& ThisComponent)
{
}

}  // namespace Kratos.

#endif // KRATOS_KRATOS_COMPONENTS_H_INCLUDED  defined
