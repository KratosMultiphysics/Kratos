//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "includes/variables_derivatives.h"

namespace Kratos {

void AddVariableTimeDerivative(Variable<double> const& rVariable, Variable<double> const& rDerivativeVariable)
{
    VariablesDerivatives<Variable<double>>::AddTimeDerivative(rVariable, rDerivativeVariable);
}

void AddVariableTimeDerivative(Variable<array_1d<double, 3>> const& rVariable, Variable<array_1d<double, 3>> const& rDerivativeVariable)
{
    VariablesDerivatives<Variable<array_1d<double, 3>>>::AddTimeDerivative(rVariable, rDerivativeVariable);
}

void AddVariableTimeDerivative(Variable<array_1d<double, 4>> const& rVariable, Variable<array_1d<double, 4>> const& rDerivativeVariable)
{
    VariablesDerivatives<Variable<array_1d<double, 4>>>::AddTimeDerivative(rVariable, rDerivativeVariable);
}

void AddVariableTimeDerivative(Variable<array_1d<double, 6>> const& rVariable, Variable<array_1d<double, 6>> const& rDerivativeVariable)
{
    VariablesDerivatives<Variable<array_1d<double, 6>>>::AddTimeDerivative(rVariable, rDerivativeVariable);
}

void AddVariableTimeDerivative(Variable<array_1d<double, 9>> const& rVariable, Variable<array_1d<double, 9>> const& rDerivativeVariable)
{
    VariablesDerivatives<Variable<array_1d<double, 9>>>::AddTimeDerivative(rVariable, rDerivativeVariable);
}

void AddVariableTimeDerivative(Variable<Vector> const& rVariable, Variable<Vector> const& rDerivativeVariable)
{
    VariablesDerivatives<Variable<Vector>>::AddTimeDerivative(rVariable, rDerivativeVariable);
}

void AddVariableTimeDerivative(Variable<Matrix> const& rVariable, Variable<Matrix> const& rDerivativeVariable)
{
    VariablesDerivatives<Variable<Matrix>>::AddTimeDerivative(rVariable, rDerivativeVariable);
}

void AddVariableTimeDerivative(VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>> const& rVariable, VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>> const& rDerivativeVariable)
{
    VariablesDerivatives<VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>>>::AddTimeDerivative(rVariable, rDerivativeVariable);
}

void AddVariableTimeDerivative(VariableComponent<VectorComponentAdaptor<array_1d<double, 4>>> const& rVariable, VariableComponent<VectorComponentAdaptor<array_1d<double, 4>>> const& rDerivativeVariable)
{
    VariablesDerivatives<VariableComponent<VectorComponentAdaptor<array_1d<double, 4>>>>::AddTimeDerivative(rVariable, rDerivativeVariable);
}

void AddVariableTimeDerivative(VariableComponent<VectorComponentAdaptor<array_1d<double, 6>>> const& rVariable, VariableComponent<VectorComponentAdaptor<array_1d<double, 6>>> const& rDerivativeVariable)
{
    VariablesDerivatives<VariableComponent<VectorComponentAdaptor<array_1d<double, 6>>>>::AddTimeDerivative(rVariable, rDerivativeVariable);
}

void AddVariableTimeDerivative(VariableComponent<VectorComponentAdaptor<array_1d<double, 9>>> const& rVariable, VariableComponent<VectorComponentAdaptor<array_1d<double, 9>>> const& rDerivativeVariable)
{
    VariablesDerivatives<VariableComponent<VectorComponentAdaptor<array_1d<double, 9>>>>::AddTimeDerivative(rVariable, rDerivativeVariable);
}

template class VariablesDerivatives<Variable<double>>;
template class VariablesDerivatives<Variable<array_1d<double, 3>>>;
template class VariablesDerivatives<Variable<array_1d<double, 4>>>;
template class VariablesDerivatives<Variable<array_1d<double, 6>>>;
template class VariablesDerivatives<Variable<array_1d<double, 9>>>;
template class VariablesDerivatives<Variable<Vector>>;
template class VariablesDerivatives<Variable<Matrix>>;
template class VariablesDerivatives<VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>>>;
template class VariablesDerivatives<VariableComponent<VectorComponentAdaptor<array_1d<double, 4>>>>;
template class VariablesDerivatives<VariableComponent<VectorComponentAdaptor<array_1d<double, 6>>>>;
template class VariablesDerivatives<VariableComponent<VectorComponentAdaptor<array_1d<double, 9>>>>;

// Specialize array of compenents for VariableData
VariablesDerivatives<VariableData>::DerivativesDatabaseType VariablesDerivatives<VariableData>::msVariablesTimeDerivatives;

}  // namespace Kratos.



