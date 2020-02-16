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
#include "includes/variables_time_derivatives.h"

namespace Kratos {

void AddVariableTimeDerivative(Variable<double> const& rVariable, Variable<double> const& rDerivativeVariable)
{
    VariablesTimeDerivatives<Variable<double>>::Add(rVariable, rDerivativeVariable);
}

void AddVariableTimeDerivative(Variable<array_1d<double, 3>> const& rVariable, Variable<array_1d<double, 3>> const& rDerivativeVariable)
{
    VariablesTimeDerivatives<Variable<array_1d<double, 3>>>::Add(rVariable, rDerivativeVariable);
}

void AddVariableTimeDerivative(Variable<array_1d<double, 4>> const& rVariable, Variable<array_1d<double, 4>> const& rDerivativeVariable)
{
    VariablesTimeDerivatives<Variable<array_1d<double, 4>>>::Add(rVariable, rDerivativeVariable);
}

void AddVariableTimeDerivative(Variable<array_1d<double, 6>> const& rVariable, Variable<array_1d<double, 6>> const& rDerivativeVariable)
{
    VariablesTimeDerivatives<Variable<array_1d<double, 6>>>::Add(rVariable, rDerivativeVariable);
}

void AddVariableTimeDerivative(Variable<array_1d<double, 9>> const& rVariable, Variable<array_1d<double, 9>> const& rDerivativeVariable)
{
    VariablesTimeDerivatives<Variable<array_1d<double, 9>>>::Add(rVariable, rDerivativeVariable);
}

void AddVariableTimeDerivative(Variable<Vector> const& rVariable, Variable<Vector> const& rDerivativeVariable)
{
    VariablesTimeDerivatives<Variable<Vector>>::Add(rVariable, rDerivativeVariable);
}

void AddVariableTimeDerivative(Variable<Matrix> const& rVariable, Variable<Matrix> const& rDerivativeVariable)
{
    VariablesTimeDerivatives<Variable<Matrix>>::Add(rVariable, rDerivativeVariable);
}

void AddVariableTimeDerivative(VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>> const& rVariable, VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>> const& rDerivativeVariable)
{
    VariablesTimeDerivatives<VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>>>::Add(rVariable, rDerivativeVariable);
}

void AddVariableTimeDerivative(VariableComponent<VectorComponentAdaptor<array_1d<double, 4>>> const& rVariable, VariableComponent<VectorComponentAdaptor<array_1d<double, 4>>> const& rDerivativeVariable)
{
    VariablesTimeDerivatives<VariableComponent<VectorComponentAdaptor<array_1d<double, 4>>>>::Add(rVariable, rDerivativeVariable);
}

void AddVariableTimeDerivative(VariableComponent<VectorComponentAdaptor<array_1d<double, 6>>> const& rVariable, VariableComponent<VectorComponentAdaptor<array_1d<double, 6>>> const& rDerivativeVariable)
{
    VariablesTimeDerivatives<VariableComponent<VectorComponentAdaptor<array_1d<double, 6>>>>::Add(rVariable, rDerivativeVariable);
}

void AddVariableTimeDerivative(VariableComponent<VectorComponentAdaptor<array_1d<double, 9>>> const& rVariable, VariableComponent<VectorComponentAdaptor<array_1d<double, 9>>> const& rDerivativeVariable)
{
    VariablesTimeDerivatives<VariableComponent<VectorComponentAdaptor<array_1d<double, 9>>>>::Add(rVariable, rDerivativeVariable);
}

template class VariablesTimeDerivatives<Variable<double>>;
template class VariablesTimeDerivatives<Variable<array_1d<double, 3>>>;
template class VariablesTimeDerivatives<Variable<array_1d<double, 4>>>;
template class VariablesTimeDerivatives<Variable<array_1d<double, 6>>>;
template class VariablesTimeDerivatives<Variable<array_1d<double, 9>>>;
template class VariablesTimeDerivatives<Variable<Vector>>;
template class VariablesTimeDerivatives<Variable<Matrix>>;
template class VariablesTimeDerivatives<VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>>>;
template class VariablesTimeDerivatives<VariableComponent<VectorComponentAdaptor<array_1d<double, 4>>>>;
template class VariablesTimeDerivatives<VariableComponent<VectorComponentAdaptor<array_1d<double, 6>>>>;
template class VariablesTimeDerivatives<VariableComponent<VectorComponentAdaptor<array_1d<double, 9>>>>;

// Specialize array of compenents for VariableData
VariablesTimeDerivatives<VariableData>::DerivativesDatabaseType VariablesTimeDerivatives<VariableData>::msVariablesTimeDerivatives;

}  // namespace Kratos.



