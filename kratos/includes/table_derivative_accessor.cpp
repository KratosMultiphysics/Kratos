//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//

// System includes

// External includes

// Project includes
#include "includes/table_derivative_accessor.h"
#include "includes/properties.h"

namespace Kratos
{

/***********************************************************************************/
/***********************************************************************************/

double TableDerivativeAccessor::GetValue(
    const Variable<double>& rVariable,
    const Properties& rProperties,
    const GeometryType& rGeometry,
    const Vector& rShapeFunctionVector,
    const ProcessInfo& rProcessInfo
    ) const
{
    return GetValueFromTable(*mpInputVariable, rVariable, rProperties, rGeometry, rShapeFunctionVector, rProcessInfo);
}

/***********************************************************************************/
/***********************************************************************************/

double TableDerivativeAccessor::GetValueFromTable(
        const Variable<double> &rIndependentVariable,
        const Variable<double> &rDependentVariable,
        const Properties& rProperties,
        const GeometryType& rGeometry,
        const Vector& rShapeFunctionVector,
        const ProcessInfo& rProcessInfo) const
{
    double independent_at_gauss = 0.0;
    if (mInputVariableType == Globals::DataLocation::NodeHistorical) {
        for (SizeType i = 0; i < rShapeFunctionVector.size(); ++i) {
            KRATOS_DEBUG_ERROR_IF_NOT(rGeometry[i].SolutionStepsDataHas(rIndependentVariable)) << "The Variable " << rIndependentVariable.Name() << " is not available at the nodes of the Geometry to retrieve Table values." << std::endl;
            const double nodal_value = rGeometry[i].FastGetSolutionStepValue(rIndependentVariable);
            independent_at_gauss += nodal_value * rShapeFunctionVector[i];
        }
    } else if (mInputVariableType == Globals::DataLocation::NodeNonHistorical) {
        for (SizeType i = 0; i < rShapeFunctionVector.size(); ++i) {
            KRATOS_DEBUG_ERROR_IF_NOT(rGeometry[i].Has(rIndependentVariable)) << "The Variable " << rIndependentVariable.Name() << " is not available at the nodes of the Geometry to retrieve Table values." << std::endl;
            const double nodal_value = rGeometry[i].GetValue(rIndependentVariable);
            independent_at_gauss += nodal_value * rShapeFunctionVector[i];
        }
    } else if (mInputVariableType == Globals::DataLocation::Element) {
        KRATOS_DEBUG_ERROR_IF_NOT(rGeometry.Has(rIndependentVariable)) << "The Variable " << rIndependentVariable.Name() << " is not available at the Geometry to retrieve Table values." << std::endl;
        independent_at_gauss = rGeometry.GetValue(rIndependentVariable);
    } else {
        KRATOS_ERROR << "The table_input_variable_type is incorrect or not supported. Types available are : nodal_historical, nodal_non_historical and elemental_non_historical" << std::endl;
    }

    // Retrieve the dependent variable from the table
    const auto& r_table = rProperties.GetTable(rIndependentVariable, rDependentVariable);
    return r_table.GetDerivative(independent_at_gauss);
}

/***********************************************************************************/
/***********************************************************************************/

void TableDerivativeAccessor::save(Serializer& rSerializer) const
{
    rSerializer.save("InputVariable", mpInputVariable->Name());
    // // we must do the int cast to be able to compile
    rSerializer.save("InputVariableType", static_cast<int>(mInputVariableType));
}
void TableDerivativeAccessor::load(Serializer& rSerializer)
{
    std::string variable_name;
    rSerializer.load("InputVariable", variable_name);
    mpInputVariable = static_cast<Variable<double> *>(KratosComponents<VariableData>::pGet(variable_name));

    // // we must do the int cast to be able to compile
    int enum_value;
    rSerializer.load("InputVariableType", enum_value);
    mInputVariableType = static_cast<Globals::DataLocation>(enum_value);
}

/***********************************************************************************/
/***********************************************************************************/

Accessor::UniquePointer TableDerivativeAccessor::Clone() const
{
    return Kratos::make_unique<TableDerivativeAccessor>(*this);
}

} // namespace Kratos
