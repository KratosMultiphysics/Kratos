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
//












// System includes

// External includes


// Project includes
#include "includes/kratos_components.h"
#include "includes/element.h"
#include "includes/condition.h"
#include "includes/constitutive_law.h"
#include "utilities/quaternion.h"

namespace Kratos
{

void AddKratosComponent(std::string const& Name, Variable<bool> const& ThisComponent)
{
    KratosComponents<Variable<bool> >::Add(Name, ThisComponent);
}

void AddKratosComponent(std::string const& Name, Variable<int> const& ThisComponent)
{
    KratosComponents<Variable<int> >::Add(Name, ThisComponent);
}

void AddKratosComponent(std::string const& Name, Variable<unsigned int> const& ThisComponent)
{
    KratosComponents<Variable<unsigned int> >::Add(Name, ThisComponent);
}

void AddKratosComponent(std::string const& Name, Variable<double> const& ThisComponent)
{
    KratosComponents<Variable<double> >::Add(Name, ThisComponent);
}

void AddKratosComponent(std::string const& Name, Variable<array_1d<double, 3> > const& ThisComponent)
{
    KratosComponents<Variable<array_1d<double, 3> > >::Add(Name, ThisComponent);
}

void AddKratosComponent(std::string const& Name, Variable<Quaternion<double> > const& ThisComponent)
{
    KratosComponents<Variable<Quaternion<double> > >::Add(Name, ThisComponent);
}

void AddKratosComponent(std::string const& Name, Variable<Vector> const& ThisComponent)
{
    KratosComponents<Variable<Vector> >::Add(Name, ThisComponent);
}

void AddKratosComponent(std::string const& Name, Variable<Matrix> const& ThisComponent)
{
    KratosComponents<Variable<Matrix> >::Add(Name, ThisComponent);
}

void AddKratosComponent(std::string const& Name, Variable<std::string> const& ThisComponent)
{
    KratosComponents<Variable<std::string> >::Add(Name, ThisComponent);
}

void AddKratosComponent(std::string const& Name, VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > > const& ThisComponent)
{
    KratosComponents<VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > > >::Add(Name, ThisComponent);
}

void AddKratosComponent(std::string const& Name, Variable<Flags> const& ThisComponent)
{
    KratosComponents<Variable<Flags> >::Add(Name, ThisComponent);
}

void AddKratosComponent(std::string const& Name, Element const& ThisComponent)
{
    KratosComponents<Element>::Add(Name, ThisComponent);
}

void AddKratosComponent(std::string const& Name, Condition const& ThisComponent)
{
    KratosComponents<Condition>::Add(Name, ThisComponent);
}

void AddKratosComponent(std::string const& Name, ConstitutiveLaw const& ThisComponent)
{
    KratosComponents<ConstitutiveLaw>::Add(Name, ThisComponent);
}

// Specialize array of compenents for VariableData
KratosComponents<VariableData>::ComponentsContainerType KratosComponents<VariableData>::msComponents;

// Explicit instantiation definition
template class KratosComponents<Variable<bool> >;
template class KratosComponents<Variable<int> >;
template class KratosComponents<Variable<unsigned int> >;
template class KratosComponents<Variable<double> >;
template class KratosComponents<Variable<array_1d<double, 3> > >;
template class KratosComponents<Variable<Quaternion<double> > >;
template class KratosComponents<Variable<Vector> >;
template class KratosComponents<Variable<Matrix> >;
template class KratosComponents<Variable<std::string> >;
template class KratosComponents<VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > > >;
template class KratosComponents<Variable<Flags> >;
template class KratosComponents<Flags>;

}  // namespace Kratos.



