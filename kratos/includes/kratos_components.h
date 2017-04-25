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
#include <boost/ref.hpp>

// Project includes
#include "includes/define.h"
#include "containers/variable_data.h"


namespace Kratos
{

/// Short class definition.
/** Detail class definition.
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
        KratosComponents<TComponentType>::GetComponentsInstance().insert(typename ComponentsContainerType::value_type(Name , &ThisComponent));
    }

    static TComponentType const& Get(std::string const& Name)
    {
        typename KratosComponents<TComponentType>::ComponentsContainerType::iterator i =  GetComponentsInstance().find(Name);
        if(i == KratosComponents<TComponentType>::GetComponentsInstance().end())
          KRATOS_THROW_ERROR(std::invalid_argument, "The component is not registered!", Name);
        return *(i->second);
    }

    static ComponentsContainerType & GetComponents()
    {
        return KratosComponents<TComponentType>::GetComponentsInstance();
    }

    static ComponentsContainerType * pGetComponents()
    {
        return &(KratosComponents<TComponentType>::GetComponentsInstance());
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
        return (GetComponentsInstance().find(Name) != GetComponentsInstance().end());
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
        for(typename ComponentsContainerType::const_iterator i = GetComponentsInstance().begin() ; i != GetComponentsInstance().end() ; ++i)
            rOStream << "    " << i->first << std::endl;
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


    ///@}
    ///@name Private Operators
    ///@{

    static ComponentsContainerType& GetComponentsInstance();

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

#define REGISTER_COMPONENT(component_type) \
    template<>  KratosComponents<component_type >::ComponentsContainerType& KratosComponents<component_type >::GetComponentsInstance() \
    { \
            static KratosComponents<component_type >::ComponentsContainerType instance; \
            return instance; \
    }
    

///@name Input and output
///@{

/// input stream function
//   template<class TComponentType>
//   inline std::istream& operator >> (std::istream& rIStream,
// 				    KratosComponents<TComponentType>& rThis);



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


template<class TComponentType> void KRATOS_API(KRATOS_CORE) AddKratosComponent(std::string const& Name, TComponentType const& ThisComponent)
{
    KratosComponents<TComponentType>::Add(Name, ThisComponent);
}



}  // namespace Kratos.


#endif // KRATOS_KRATOS_COMPONENTS_H_INCLUDED  defined
