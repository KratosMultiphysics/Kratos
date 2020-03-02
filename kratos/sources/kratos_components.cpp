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
#include "geometries/register_kratos_components_for_geometry.h"
#include "includes/element.h"
#include "includes/condition.h"
#include "includes/constitutive_law.h"
#include "includes/master_slave_constraint.h"

/* Utilities */
#include "utilities/quaternion.h"

/* Factories */
#include "factories/linear_solver_factory.h"
#include "factories/preconditioner_factory.h"

namespace Kratos {

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

void AddKratosComponent(std::string const& Name, Variable<array_1d<double, 4> > const& ThisComponent)
{
    KratosComponents<Variable<array_1d<double, 4> > >::Add(Name, ThisComponent);
}

void AddKratosComponent(std::string const& Name, Variable<array_1d<double, 6> > const& ThisComponent)
{
    KratosComponents<Variable<array_1d<double, 6> > >::Add(Name, ThisComponent);
}

void AddKratosComponent(std::string const& Name, Variable<array_1d<double, 9> > const& ThisComponent)
{
    KratosComponents<Variable<array_1d<double, 9> > >::Add(Name, ThisComponent);
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

void AddKratosComponent(std::string const& Name, VariableComponent<VectorComponentAdaptor<array_1d<double, 4> > > const& ThisComponent)
{
    KratosComponents<VariableComponent<VectorComponentAdaptor<array_1d<double, 4> > > >::Add(Name, ThisComponent);
}

void AddKratosComponent(std::string const& Name, VariableComponent<VectorComponentAdaptor<array_1d<double, 6> > > const& ThisComponent)
{
    KratosComponents<VariableComponent<VectorComponentAdaptor<array_1d<double, 6> > > >::Add(Name, ThisComponent);
}

void AddKratosComponent(std::string const& Name, VariableComponent<VectorComponentAdaptor<array_1d<double, 9> > > const& ThisComponent)
{
    KratosComponents<VariableComponent<VectorComponentAdaptor<array_1d<double, 9> > > >::Add(Name, ThisComponent);
}

void AddKratosComponent(std::string const& Name, Variable<Flags> const& ThisComponent)
{
    KratosComponents<Variable<Flags> >::Add(Name, ThisComponent);
}

void AddKratosComponent(std::string const& Name, Geometry<Node<3>> const& ThisComponent)
{
    KratosComponents<Geometry<Node<3>>>::Add(Name, ThisComponent);
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

void AddKratosComponent(std::string const& Name, Variable<ConstitutiveLaw::Pointer> const& ThisComponent)
{
    KratosComponents<Variable<ConstitutiveLaw::Pointer> >::Add(Name, ThisComponent);
}

template class KratosComponents<Variable<bool> >;
template class KratosComponents<Variable<int> >;
template class KratosComponents<Variable<unsigned int> >;
template class KratosComponents<Variable<double> >;
template class KratosComponents<Variable<array_1d<double, 3> > >;
template class KratosComponents<Variable<array_1d<double, 4> > >;
template class KratosComponents<Variable<array_1d<double, 6> > >;
template class KratosComponents<Variable<array_1d<double, 9> > >;
template class KratosComponents<Variable<Quaternion<double> > >;
template class KratosComponents<Variable<Vector> >;
template class KratosComponents<Variable<Matrix> >;
template class KratosComponents<Variable<std::string> >;
template class KratosComponents<VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > > >;
template class KratosComponents<VariableComponent<VectorComponentAdaptor<array_1d<double, 4> > > >;
template class KratosComponents<VariableComponent<VectorComponentAdaptor<array_1d<double, 6> > > >;
template class KratosComponents<VariableComponent<VectorComponentAdaptor<array_1d<double, 9> > > >;
template class KratosComponents<Variable<Flags> >;
template class KratosComponents<Flags>;
template class KratosComponents<DataCommunicator>;

template class KratosComponents<Geometry<Node<3>>>;
template class KratosComponents<Element>;
template class KratosComponents<Condition>;
template class KratosComponents<ConstitutiveLaw>;
template class KratosComponents<Variable<ConstitutiveLaw::Pointer>>;
template class KratosComponents<MasterSlaveConstraint>;

using RealSparseSpace = UblasSpace<double, boost::numeric::ublas::compressed_matrix<double>, boost::numeric::ublas::vector<double>>;
using RealDenseSpace = UblasSpace<double, DenseMatrix<double>, DenseVector<double>>;
using ComplexSparseSpace = UblasSpace<std::complex<double>, boost::numeric::ublas::compressed_matrix<std::complex<double>>, boost::numeric::ublas::vector<std::complex<double>>>;
using ComplexDenseSpace = UblasSpace<std::complex<double>, DenseMatrix<std::complex<double>>, DenseVector<std::complex<double>>>;

template class KratosComponents<LinearSolverFactory<RealSparseSpace, RealDenseSpace>>;
template class KratosComponents<LinearSolverFactory<ComplexSparseSpace, ComplexDenseSpace>>;
template class KratosComponents<PreconditionerFactory<RealSparseSpace, RealDenseSpace>>;

// Specialize array of compenents for VariableData
KratosComponents<VariableData>::ComponentsContainerType KratosComponents<VariableData>::msComponents;

}  // namespace Kratos.



