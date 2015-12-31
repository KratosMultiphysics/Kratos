// Kratos Multi-Physics
//
// Copyright (c) 2016 Pooyan Dadvand, Riccardo Rossi, CIMNE (International Center for Numerical Methods in Engineering)
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
//
// 	-	Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
// 	-	Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer
// 		in the documentation and/or other materials provided with the distribution.
// 	-	All advertising materials mentioning features or use of this software must display the following acknowledgement:
// 			This product includes Kratos Multi-Physics technology.
// 	-	Neither the name of the CIMNE nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS ''AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
// THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// HOLDERS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED ANDON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF
// THE USE OF THISSOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


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
#include "containers/variable_component.h"
#include "containers/vector_component_adaptor.h"
#include "containers/flags.h"

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

//       KratosComponents(std::string const& Name, TComponentType const& ThisComponent)
//       {
//  	msComponents.insert(typename ComponentsContainerType::value_type(Name ,boost::cref(ThisComponent)));
//       }

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
        msComponents.insert(typename ComponentsContainerType::value_type(Name , &ThisComponent));
    }


//     static void Add(std::string const& Name, TComponentType const& ThisComponent, ComponentsContainerType& ThisComponents)
//     {
    ////ThisComponents.insert(typename ComponentsContainerType::value_type(Name ,boost::cref(ThisComponent)));
    //msComponents.insert(typename ComponentsContainerType::value_type(Name ,boost::cref(ThisComponent)));
//     }


    static TComponentType const& Get(std::string const& Name)
    {
        typename ComponentsContainerType::iterator i =  msComponents.find(Name);
        /* 	if(i == msComponents.end()) */
        /* 	   KRATOS_THROW_ERROR(std::invalid_argument, "The component is not registered!", Name); */
        return *(i->second);
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
        for(typename ComponentsContainerType::const_iterator i = msComponents.begin() ; i != msComponents.end() ; ++i)
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

//       KratosComponents(std::string const& Name, TComponentType const& ThisComponent)
//       {
//  	msComponents.insert(typename ComponentsContainerType::value_type(Name ,boost::cref(ThisComponent)));
//       }

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
        msComponents.insert(ComponentsContainerType::value_type(Name ,&ThisComponent));
    }

    static std::size_t Size()
    {
        return msComponents.size();
    }

//     static void Add(std::string const& Name, VariableData& ThisComponent, ComponentsContainerType& ThisComponents)
//     {
    ////ThisComponents.insert(typename ComponentsContainerType::value_type(Name ,boost::cref(ThisComponent)));
    //msComponents.insert(typename ComponentsContainerType::value_type(Name ,boost::ref(ThisComponent)));
//     }

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
        for(ComponentsContainerType::const_iterator i = msComponents.begin() ; i != msComponents.end() ; ++i)
            rOStream << "    " << *(i->second) << std::endl;
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
template class KRATOS_API(KRATOS_CORE) KratosComponents<Variable<Vector> >;
template class KRATOS_API(KRATOS_CORE) KratosComponents<Variable<Matrix> >;
template class KRATOS_API(KRATOS_CORE) KratosComponents<Variable<std::string> >;
template class KRATOS_API(KRATOS_CORE) KratosComponents<VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > > >;
template class KRATOS_API(KRATOS_CORE) KratosComponents<Variable<Flags> >;
template class KRATOS_API(KRATOS_CORE) KratosComponents<Flags>;

#ifdef KratosCore_EXPORTS
template<class TComponentType>
typename KratosComponents<TComponentType>::ComponentsContainerType KratosComponents<TComponentType>::msComponents;
#endif

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

void KRATOS_API(KRATOS_CORE) AddKratosComponent(std::string const& Name, Variable<bool> const& ThisComponent);
void KRATOS_API(KRATOS_CORE) AddKratosComponent(std::string const& Name, Variable<int> const& ThisComponent);
void KRATOS_API(KRATOS_CORE) AddKratosComponent(std::string const& Name, Variable<unsigned int> const& ThisComponent);
void KRATOS_API(KRATOS_CORE) AddKratosComponent(std::string const& Name, Variable<double> const& ThisComponent);
void KRATOS_API(KRATOS_CORE) AddKratosComponent(std::string const& Name, Variable<array_1d<double, 3> > const& ThisComponent);
void KRATOS_API(KRATOS_CORE) AddKratosComponent(std::string const& Name, Variable<Vector> const& ThisComponent);
void KRATOS_API(KRATOS_CORE) AddKratosComponent(std::string const& Name, Variable<Matrix> const& ThisComponent);
void KRATOS_API(KRATOS_CORE) AddKratosComponent(std::string const& Name, Variable<std::string> const& ThisComponent);
void KRATOS_API(KRATOS_CORE) AddKratosComponent(std::string const& Name, VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > > const& ThisComponent);
void KRATOS_API(KRATOS_CORE) AddKratosComponent(std::string const& Name, Flags const& ThisComponent);
void KRATOS_API(KRATOS_CORE) AddKratosComponent(std::string const& Name, Variable<Flags> const& ThisComponent);

template<class TComponentType> void AddKratosComponent(std::string const& Name, TComponentType const& ThisComponent)
{
}

}  // namespace Kratos.

#endif // KRATOS_KRATOS_COMPONENTS_H_INCLUDED  defined 





