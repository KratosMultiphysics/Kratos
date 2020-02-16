//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                     Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#if !defined(KRATOS_VARIABLES_TIME_DERIVATIVES_H_INCLUDED )
#define  KRATOS_VARIABLES_TIME_DERIVATIVES_H_INCLUDED

// System includes
#include <map>

// External includes

// Project includes
#include "includes/define.h"
#include "containers/variable_component.h"
#include "containers/vector_component_adaptor.h"

namespace Kratos
{

/**
 * @class VariablesTimeDerivatives
 * @ingroup KratosCore
 * @brief This container defines the time derivatives of the variables registered
 * @details It stores in a similar way to the KratosComponents but applied to derivatives. The map uses the key of the variable
 * @author Vicente Mataix Ferrandiz
 * @tparam TComponentType The component type (type of variable)
 */
template<class TComponentType>
class KRATOS_API(KRATOS_CORE) VariablesTimeDerivatives
{
public:

    ///@name Type Definitions
    ///@{

    /// Pointer definition of VariablesTimeDerivatives
    KRATOS_CLASS_POINTER_DEFINITION(VariablesTimeDerivatives);

    /// Definition of the database for the derivatives
    typedef std::map<std::size_t, const TComponentType* > DerivativesDatabaseType;
    typedef typename DerivativesDatabaseType::value_type ValueType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    VariablesTimeDerivatives() {}

    /// Destructor.
    virtual ~VariablesTimeDerivatives() {}

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    static void Add(TComponentType const& rVariable, TComponentType const& rDerivativeVariable)
    {
        // check if a different object was already registered with this name, since this is undefined behavior
        const std::size_t key = rVariable.Key();
        const auto it_der = msVariablesTimeDerivatives.find(key);
        KRATOS_ERROR_IF(it_der != msVariablesTimeDerivatives.end() && typeid(*(it_der->second)) != typeid(rDerivativeVariable)) << "An object of different type was already registered with name \"" << rVariable.Key() << "\"!" << std::endl;
        msVariablesTimeDerivatives.insert(ValueType(rVariable.Key(), &rDerivativeVariable));
    }

    static void Remove(TComponentType const& rVariable)
    {
        const std::size_t num_erased = msVariablesTimeDerivatives.erase(rVariable.Key());
        KRATOS_ERROR_IF(num_erased == 0) << "Trying to remove inexistent component \"" << rVariable.Key() << "\"." << std::endl;
    }

    static TComponentType const& GetFirstDerivative(TComponentType const& rVariable)
    {
        const auto it_der = msVariablesTimeDerivatives.find(rVariable.Key());
        KRATOS_DEBUG_ERROR_IF(it_der == msVariablesTimeDerivatives.end()) << GetMessageUnregisteredDerivative(rVariable) << std::endl;
        return *(it_der->second);
    }

    static DerivativesDatabaseType & GetVariableTimeDerivatives()
    {
        return msVariablesTimeDerivatives;
    }

    static DerivativesDatabaseType * pGetVariableTimeDerivatives()
    {
        return &msVariablesTimeDerivatives;
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

    static bool Has(TComponentType const& rVariable)
    {
        return (msVariablesTimeDerivatives.find(rVariable.Key()) != msVariablesTimeDerivatives.end());
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "Variables time derivatives";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "Variables time derivatives";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
        for (const auto& r_comp : msVariablesTimeDerivatives) {
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

    static DerivativesDatabaseType msVariablesTimeDerivatives;

    ///@}
    ///@name Member Variables
    ///@{


    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    static std::string GetMessageUnregisteredDerivative(TComponentType const& rVariable)
    {
        std::stringstream msg;
        msg << "The derivative for \"" << rVariable.Key() << "\" is not registered!\nMaybe you need to import the application where it is defined?\nThe following components of this type are registered:" << std::endl;
        VariablesTimeDerivatives instance; // creating an instance for using "PrintData"
        instance.PrintData(msg);
        return msg.str();
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

    /// Assignment operator.
    VariablesTimeDerivatives& operator=(VariablesTimeDerivatives const& rOther);

    /// Copy constructor.
    VariablesTimeDerivatives(VariablesTimeDerivatives const& rOther);

    ///@}

}; // Class VariablesTimeDerivatives


///@}
template<>
class KRATOS_API(KRATOS_CORE) VariablesTimeDerivatives<VariableData>
{
public:

    ///@name Type Definitions
    ///@{

    /// Pointer definition of VariablesTimeDerivatives
    KRATOS_CLASS_POINTER_DEFINITION(VariablesTimeDerivatives);

    /// Definition of the database for the derivatives
    typedef std::map<std::size_t, const VariableData* > DerivativesDatabaseType;
    typedef DerivativesDatabaseType::value_type ValueType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    VariablesTimeDerivatives() {}

    /// Destructor.
    virtual ~VariablesTimeDerivatives() {}

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    static void Add(VariableData const& rVariable, VariableData& rDerivativeVariable)
    {
        msVariablesTimeDerivatives.insert(ValueType(rVariable.Key(), &rDerivativeVariable));
    }

    static void Remove(VariableData const& rVariable)
    {
        std::size_t num_erased = msVariablesTimeDerivatives.erase(rVariable.Key());
        KRATOS_ERROR_IF(num_erased == 0) << "Trying to remove inexistent derivative \"" << rVariable.Key() << "\"." << std::endl;
    }

    static std::size_t Size()
    {
        return msVariablesTimeDerivatives.size();
    }

    static const VariableData & GetFirstDerivative(VariableData const& rVariable)
    {
        const auto it_der = msVariablesTimeDerivatives.find(rVariable.Key());
        KRATOS_DEBUG_ERROR_IF(it_der == msVariablesTimeDerivatives.end()) << GetMessageUnregisteredVariable(rVariable) << std::endl;
        return *(it_der->second);
    }

    static const VariableData* pGetFirstDerivative(VariableData const& rVariable)
    {
        const auto it_der = msVariablesTimeDerivatives.find(rVariable.Key());
        KRATOS_DEBUG_ERROR_IF(it_der == msVariablesTimeDerivatives.end()) << GetMessageUnregisteredVariable(rVariable) << std::endl;
        return it_der->second;
    }

    static const VariableData& GetSecondDerivative(VariableData const& rVariable)
    {
        const auto it_der = msVariablesTimeDerivatives.find(rVariable.Key());
        KRATOS_DEBUG_ERROR_IF(it_der == msVariablesTimeDerivatives.end()) << GetMessageUnregisteredVariable(rVariable) << std::endl;
        return GetFirstDerivative(*(it_der->second));
    }

    static const VariableData* pGetSecondDerivative(VariableData const& rVariable)
    {
        const auto it_der = msVariablesTimeDerivatives.find(rVariable.Key());
        KRATOS_DEBUG_ERROR_IF(it_der == msVariablesTimeDerivatives.end()) << GetMessageUnregisteredVariable(rVariable) << std::endl;
        return pGetFirstDerivative(*(it_der->second));
    }

    static DerivativesDatabaseType & GetVariableTimeDerivatives()
    {
        return msVariablesTimeDerivatives;
    }

    static DerivativesDatabaseType * pGetVariableTimeDerivatives()
    {
        return &msVariablesTimeDerivatives;
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

    static bool Has(VariableData const& rVariable)
    {
        return (msVariablesTimeDerivatives.find(rVariable.Key()) != msVariablesTimeDerivatives.end());
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "Variables time derivatives <VariableData>";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "Variables time derivatives <VariableData>";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
        for (const auto& r_comp : msVariablesTimeDerivatives) {
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

    static DerivativesDatabaseType msVariablesTimeDerivatives;

    ///@}
    ///@name Member Variables
    ///@{


    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    static std::string GetMessageUnregisteredVariable(VariableData const& rVariable)
    {
        std::stringstream msg;
        msg << "The variable \"" << rVariable.Key() << "\" is not registered!\nMaybe you need to import the application where it is defined?\nThe following variables are registered:" << std::endl;
        VariablesTimeDerivatives instance; // creating an instance for using "PrintData"
        instance.PrintData(msg);
        return msg.str();
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

    /// Assignment operator.
    VariablesTimeDerivatives& operator=(VariablesTimeDerivatives const& rOther);

    /// Copy constructor.
    VariablesTimeDerivatives(VariablesTimeDerivatives const& rOther);

    ///@}

}; // Class VariablesTimeDerivatives

template<class TComponentType>
typename VariablesTimeDerivatives<TComponentType>::DerivativesDatabaseType VariablesTimeDerivatives<TComponentType>::msVariablesTimeDerivatives;

KRATOS_API_EXTERN template class KRATOS_API(KRATOS_CORE) VariablesTimeDerivatives<Variable<double>>;
KRATOS_API_EXTERN template class KRATOS_API(KRATOS_CORE) VariablesTimeDerivatives<Variable<array_1d<double, 3>>>;
KRATOS_API_EXTERN template class KRATOS_API(KRATOS_CORE) VariablesTimeDerivatives<Variable<array_1d<double, 4>>>;
KRATOS_API_EXTERN template class KRATOS_API(KRATOS_CORE) VariablesTimeDerivatives<Variable<array_1d<double, 6>>>;
KRATOS_API_EXTERN template class KRATOS_API(KRATOS_CORE) VariablesTimeDerivatives<Variable<array_1d<double, 9>>>;
KRATOS_API_EXTERN template class KRATOS_API(KRATOS_CORE) VariablesTimeDerivatives<Variable<Vector>>;
KRATOS_API_EXTERN template class KRATOS_API(KRATOS_CORE) VariablesTimeDerivatives<Variable<Matrix>>;
KRATOS_API_EXTERN template class KRATOS_API(KRATOS_CORE) VariablesTimeDerivatives<VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>>>;
KRATOS_API_EXTERN template class KRATOS_API(KRATOS_CORE) VariablesTimeDerivatives<VariableComponent<VectorComponentAdaptor<array_1d<double, 4>>>>;
KRATOS_API_EXTERN template class KRATOS_API(KRATOS_CORE) VariablesTimeDerivatives<VariableComponent<VectorComponentAdaptor<array_1d<double, 6>>>>;
KRATOS_API_EXTERN template class KRATOS_API(KRATOS_CORE) VariablesTimeDerivatives<VariableComponent<VectorComponentAdaptor<array_1d<double, 9>>>>;

///@name Input and output
///@{

/// output stream function
template<class TComponentType>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const VariablesTimeDerivatives<TComponentType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

void KRATOS_API(KRATOS_CORE) AddVariableTimeDerivative(Variable<double> const& rVariable, Variable<double> const& rDerivativeVariable);
void KRATOS_API(KRATOS_CORE) AddVariableTimeDerivative(Variable<array_1d<double, 3>> const& rVariable, Variable<array_1d<double, 3>> const& rDerivativeVariable);
void KRATOS_API(KRATOS_CORE) AddVariableTimeDerivative(Variable<array_1d<double, 4>> const& rVariable, Variable<array_1d<double, 4>> const& rDerivativeVariable);
void KRATOS_API(KRATOS_CORE) AddVariableTimeDerivative(Variable<array_1d<double, 6>> const& rVariable, Variable<array_1d<double, 6>> const& rDerivativeVariable);
void KRATOS_API(KRATOS_CORE) AddVariableTimeDerivative(Variable<array_1d<double, 9>> const& rVariable, Variable<array_1d<double, 9>> const& rDerivativeVariable);
void KRATOS_API(KRATOS_CORE) AddVariableTimeDerivative(Variable<Vector> const& rVariable, Variable<Vector> const& rDerivativeVariable);
void KRATOS_API(KRATOS_CORE) AddVariableTimeDerivative(Variable<Matrix> const& rVariable, Variable<Matrix> const& rDerivativeVariable);
void KRATOS_API(KRATOS_CORE) AddVariableTimeDerivative(VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>> const& rVariable, VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>> const& rDerivativeVariable);
void KRATOS_API(KRATOS_CORE) AddVariableTimeDerivative(VariableComponent<VectorComponentAdaptor<array_1d<double, 4>>> const& rVariable, VariableComponent<VectorComponentAdaptor<array_1d<double, 4>>> const& rDerivativeVariable);
void KRATOS_API(KRATOS_CORE) AddVariableTimeDerivative(VariableComponent<VectorComponentAdaptor<array_1d<double, 6>>> const& rVariable, VariableComponent<VectorComponentAdaptor<array_1d<double, 6>>> const& rDerivativeVariable);
void KRATOS_API(KRATOS_CORE) AddVariableTimeDerivative(VariableComponent<VectorComponentAdaptor<array_1d<double, 9>>> const& rVariable, VariableComponent<VectorComponentAdaptor<array_1d<double, 9>>> const& rDerivativeVariable);

template<class TComponentType> void AddVariableTimeDerivative(TComponentType const& rVariable, TComponentType const& rDerivativeVariable)
{
}

}  // namespace Kratos.

#endif // KRATOS_VARIABLES_TIME_DERIVATIVES_H_INCLUDED  defined
