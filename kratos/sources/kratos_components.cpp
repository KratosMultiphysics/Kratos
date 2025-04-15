//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
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
#include "modeler/modeler.h"

/* Utilities */
#include "utilities/quaternion.h"

/* Factories */
#include "factories/register_factories.h"
#include "factories/linear_solver_factory.h"
#include "factories/preconditioner_factory.h"

namespace Kratos {

void AddKratosComponent(const std::string& rName, const Variable<bool>& rComponent)
{
    KratosComponents<Variable<bool>>::Add(rName, rComponent);
}

void AddKratosComponent(const std::string& rName, const Variable<int>& rComponent)
{
    KratosComponents<Variable<int>>::Add(rName, rComponent);
}

void AddKratosComponent(const std::string& rName, const Variable<unsigned int>& rComponent)
{
    KratosComponents<Variable<unsigned int>>::Add(rName, rComponent);
}

void AddKratosComponent(const std::string& rName, const Variable<double>& rComponent)
{
    KratosComponents<Variable<double>>::Add(rName, rComponent);
}

void AddKratosComponent(const std::string& rName, const Variable<array_1d<double, 3>>& rComponent)
{
    KratosComponents<Variable<array_1d<double, 3>>>::Add(rName, rComponent);
}

void AddKratosComponent(const std::string& rName, const Variable<array_1d<double, 4>>& rComponent)
{
    KratosComponents<Variable<array_1d<double, 4>>>::Add(rName, rComponent);
}

void AddKratosComponent(const std::string& rName, const Variable<array_1d<double, 6>>& rComponent)
{
    KratosComponents<Variable<array_1d<double, 6>>>::Add(rName, rComponent);
}

void AddKratosComponent(const std::string& rName, const Variable<array_1d<double, 9>>& rComponent)
{
    KratosComponents<Variable<array_1d<double, 9>>>::Add(rName, rComponent);
}

void AddKratosComponent(const std::string& rName, const Variable<Quaternion<double>>& rComponent)
{
    KratosComponents<Variable<Quaternion<double>>>::Add(rName, rComponent);
}

void AddKratosComponent(const std::string& rName, const Variable<Vector>& rComponent)
{
    KratosComponents<Variable<Vector>>::Add(rName, rComponent);
}

void AddKratosComponent(const std::string& rName, const Variable<Matrix>& rComponent)
{
    KratosComponents<Variable<Matrix>>::Add(rName, rComponent);
}

void AddKratosComponent(const std::string& rName, const Variable<std::string>& rComponent)
{
    KratosComponents<Variable<std::string>>::Add(rName, rComponent);
}

void AddKratosComponent(const std::string& rName, const Variable<Flags>& rComponent)
{
    KratosComponents<Variable<Flags>>::Add(rName, rComponent);
}

void AddKratosComponent(const std::string& rName, const Geometry<Node>& rComponent)
{
    KratosComponents<Geometry<Node>>::Add(rName, rComponent);
}

void AddKratosComponent(const std::string& rName, const Element& rComponent)
{
    KratosComponents<Element>::Add(rName, rComponent);
}

void AddKratosComponent(const std::string& rName, const Condition& rComponent)
{
    KratosComponents<Condition>::Add(rName, rComponent);
}

void AddKratosComponent(const std::string& rName, const Modeler& rComponent)
{
    KratosComponents<Modeler>::Add(rName, rComponent);
}

void AddKratosComponent(const std::string& rName, const ConstitutiveLaw& rComponent)
{
    KratosComponents<ConstitutiveLaw>::Add(rName, rComponent);
}

void AddKratosComponent(const std::string& rName, const Variable<ConstitutiveLaw::Pointer>& rComponent)
{
    KratosComponents<Variable<ConstitutiveLaw::Pointer>>::Add(rName, rComponent);
}

template class KratosComponents<Variable<bool>>;
template class KratosComponents<Variable<int>>;
template class KratosComponents<Variable<unsigned int>>;
template class KratosComponents<Variable<double>>;
template class KratosComponents<Variable<array_1d<double, 3>>>;
template class KratosComponents<Variable<array_1d<double, 4>>>;
template class KratosComponents<Variable<array_1d<double, 6>>>;
template class KratosComponents<Variable<array_1d<double, 9>>>;
template class KratosComponents<Variable<Quaternion<double>>>;
template class KratosComponents<Variable<Vector>>;
template class KratosComponents<Variable<Matrix>>;
template class KratosComponents<Variable<std::string>>;
template class KratosComponents<Variable<Flags>>;
template class KratosComponents<Flags>;
template class KratosComponents<DataCommunicator>;

template class KratosComponents<Geometry<Node>>;
template class KratosComponents<Element>;
template class KratosComponents<Condition>;
template class KratosComponents<ConstitutiveLaw>;
template class KratosComponents<Variable<ConstitutiveLaw::Pointer>>;
template class KratosComponents<MasterSlaveConstraint>;
template class KratosComponents<Modeler>;

using RealSparseSpace = UblasSpace<double, boost::numeric::ublas::compressed_matrix<double>, boost::numeric::ublas::vector<double>>;
using SinglePrecisionRealSparseSpace = UblasSpace<float, boost::numeric::ublas::compressed_matrix<float>, boost::numeric::ublas::vector<float>>;
using RealDenseSpace = UblasSpace<double, DenseMatrix<double>, DenseVector<double>>;
using ComplexSparseSpace = UblasSpace<std::complex<double>, boost::numeric::ublas::compressed_matrix<std::complex<double>>, boost::numeric::ublas::vector<std::complex<double>>>;
using ComplexDenseSpace = UblasSpace<std::complex<double>, DenseMatrix<std::complex<double>>, DenseVector<std::complex<double>>>;

template class KratosComponents<LinearSolverFactory<RealSparseSpace, RealDenseSpace>>;
template class KratosComponents<LinearSolverFactory<SinglePrecisionRealSparseSpace, RealDenseSpace>>;
template class KratosComponents<LinearSolverFactory<ComplexSparseSpace, ComplexDenseSpace>>;
template class KratosComponents<PreconditionerFactory<RealSparseSpace, RealDenseSpace>>;
template class KratosComponents<ExplicitBuilder<RealSparseSpace, RealDenseSpace>>;

// Specialize array of components for VariableData
KratosComponents<VariableData>::ComponentsContainerType KratosComponents<VariableData>::msComponents;

}  // namespace Kratos.



