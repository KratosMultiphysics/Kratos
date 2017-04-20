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
#include "containers/variable_component.h"
#include "containers/vector_component_adaptor.h"
#include "containers/flags.h"
#include "utilities/quaternion.h"
#include "geometries/point.h"

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
//  	GetComponentsInstance().insert(typename ComponentsContainerType::value_type(Name ,boost::cref(ThisComponent)));
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
        KratosComponents<TComponentType>::GetComponentsInstance().insert(typename ComponentsContainerType::value_type(Name , &ThisComponent));
    }


//     static void Add(std::string const& Name, TComponentType const& ThisComponent, ComponentsContainerType& ThisComponents)
//     {
    ////ThisComponents.insert(typename ComponentsContainerType::value_type(Name ,boost::cref(ThisComponent)));
    //GetComponentsInstance().insert(typename ComponentsContainerType::value_type(Name ,boost::cref(ThisComponent)));
//     }


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

    static ComponentsContainerType& GetComponentsInstance()
	  {
		  static ComponentsContainerType instance;
		  return instance;
	  }



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
//  	GetComponentsInstance().insert(typename ComponentsContainerType::value_type(Name ,boost::cref(ThisComponent)));
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
        GetComponentsInstance().insert(ComponentsContainerType::value_type(Name ,&ThisComponent));
    }

    static std::size_t Size()
    {
        return GetComponentsInstance().size();
    }

//     static void Add(std::string const& Name, VariableData& ThisComponent, ComponentsContainerType& ThisComponents)
//     {
    ////ThisComponents.insert(typename ComponentsContainerType::value_type(Name ,boost::cref(ThisComponent)));
    //GetComponentsInstance().insert(typename ComponentsContainerType::value_type(Name ,boost::ref(ThisComponent)));
//     }

    static VariableData & Get(std::string const& Name)
    {
        return *(GetComponentsInstance().find(Name)->second);
    }

    static VariableData* pGet(std::string const& Name)
    {
        return (GetComponentsInstance().find(Name)->second);
    }

    static ComponentsContainerType & GetComponents()
    {
        return GetComponentsInstance();
    }

    static ComponentsContainerType * pGetComponents()
    {
        return &GetComponentsInstance();
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
        for(ComponentsContainerType::const_iterator i = GetComponentsInstance().begin() ; i != GetComponentsInstance().end() ; ++i)
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

    ///@}
    ///@name Member Variables
    ///@{


    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{


        static ComponentsContainerType& GetComponentsInstance()
    	  {
    		  static ComponentsContainerType instance;
    		  return instance;
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
    KratosComponents& operator=(KratosComponents const& rOther);

    /// Copy constructor.
    KratosComponents(KratosComponents const& rOther);

    ///@}

}; // Class KratosComponents




#ifdef KratosCore_EXPORTS
    template class KRATOS_API(KRATOS_CORE) KratosComponents<Variable<bool> >;
    template class KRATOS_API(KRATOS_CORE) KratosComponents<Variable<int> >;
    template class KRATOS_API(KRATOS_CORE) KratosComponents<Variable<unsigned int> >;
    template class KRATOS_API(KRATOS_CORE) KratosComponents<Variable<double> >;
    template class KRATOS_API(KRATOS_CORE) KratosComponents<Variable<array_1d<double, 3> > >;
    template class KRATOS_API(KRATOS_CORE) KratosComponents<Variable<Quaternion<double> > >;
    template class KRATOS_API(KRATOS_CORE) KratosComponents<Variable<Vector> >;
    template class KRATOS_API(KRATOS_CORE) KratosComponents<Variable<Matrix> >;
    template class KRATOS_API(KRATOS_CORE) KratosComponents<Variable<std::string> >;
    template class KRATOS_API(KRATOS_CORE) KratosComponents<VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > > >;
    template class KRATOS_API(KRATOS_CORE) KratosComponents<Variable<Flags> >;
    template class KRATOS_API(KRATOS_CORE) KratosComponents<Flags>;
    template class KRATOS_API(KRATOS_CORE) KratosComponents<Point<3,double> >;
#else
//     extern template< class TComponentType> typename KratosComponents<TComponentType>::ComponentsContainerType KratosComponents<TComponentType>::GetComponentsInstance();

    // extern template class  KratosComponents<Variable<bool> >;
    // extern template class  KratosComponents<Variable<int> >;
    // extern template class  KratosComponents<Variable<unsigned int> >;
    // extern template class  KratosComponents<Variable<double> >;
    // extern template class  KratosComponents<Variable<array_1d<double, 3> > >;
    // extern template class  KratosComponents<Variable<Quaternion<double> > >;
    // extern template class  KratosComponents<Variable<Vector> >;
    // extern template class  KratosComponents<Variable<Matrix> >;
    // extern template class  KratosComponents<Variable<std::string> >;
    // extern template class  KratosComponents<VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > > >;
    // extern template class  KratosComponents<Variable<Flags> >;
    // extern template class  KratosComponents<Flags>;
    // extern template class  KratosComponents<Point<3,double> >;
    //
    // class Element;
    // extern template class  KratosComponents<Element>;
    //
    // class Condition;
    // extern template class  KratosComponents<Condition>;
    //
    // class ConstitutiveLaw;
    // extern template class  KratosComponents<ConstitutiveLaw>;
    //
    // class PeriodicCondition;
    // extern template class  KratosComponents<PeriodicCondition>;

//
//     template<std::size_t TDimension, class TDataType = double> class Point;
//     extern template class  KratosComponents<Point<2,double> >;
//     extern template class  KratosComponents<Point<3,double> >;

    //extern template typename KratosComponents<Variable<bool> >::ComponentsContainerType KratosComponents<Variable<bool> >::GetComponentsInstance();
//     extern KratosComponents<Variable<int> >::GetComponentsInstance();
//     extern KratosComponents<Variable<unsigned int> >::GetComponentsInstance();
//     extern KratosComponents<Variable<double> >::GetComponentsInstance();
//     extern KratosComponents<Variable<array_1d<double, 3> > >::GetComponentsInstance();
//     extern KratosComponents<Variable<Quaternion<double> > >::GetComponentsInstance();
//     extern KratosComponents<Variable<Vector> >::GetComponentsInstance();
//     extern KratosComponents<Variable<Matrix> >::GetComponentsInstance();
//     extern KratosComponents<Variable<std::string> >::GetComponentsInstance();
//     extern KratosComponents<VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > > >::GetComponentsInstance();
//     extern KratosComponents<Variable<Flags> >::GetComponentsInstance();
//     extern KratosComponents<Flags>::GetComponentsInstance();
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
void KRATOS_API(KRATOS_CORE) AddKratosComponent(std::string const& Name, Variable<Quaternion<double> > const& ThisComponent);
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
