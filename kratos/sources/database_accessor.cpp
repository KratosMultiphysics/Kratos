//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Alejandro Cornejo
//                   Riccardo Rossi
//                   Carlos Roig
//
//

// System includes

// External includes

// Project includes
#include "includes/database_accessor.h"
#include "includes/properties.h"

namespace Kratos
{

/***********************************************************************************/
/***********************************************************************************/

double DatabaseAccessor::GetValue(
    const Variable<double>& rVariable,
    const Properties& rProperties,
    const GeometryType& rGeometry,
    const Vector& rShapeFunctionVector,
    const ProcessInfo& rProcessInfo
    ) const
{
    return GetValueFromDatabase(rVariable, rProperties, rGeometry, rShapeFunctionVector, rProcessInfo);
}

/***********************************************************************************/
/***********************************************************************************/

double DatabaseAccessor::GetValueFromDatabase(
        const Variable<double> &rVariable,
        const Properties& rProperties,
        const GeometryType& rGeometry,
        const Vector& rShapeFunctionVector,
        const ProcessInfo& rProcessInfo) const
{
    double value = 0.0;
    if (mInputVariableType == Globals::DataLocation::NodeHistorical) {
        for (SizeType i = 0; i < rShapeFunctionVector.size(); ++i) {
            KRATOS_DEBUG_ERROR_IF_NOT(rGeometry[i].SolutionStepsDataHas(rVariable)) << "The Variable " << rVariable.Name() << " is not available at the nodes of the Geometry to retrieve Table values." << std::endl;
            const double nodal_value = rGeometry[i].FastGetSolutionStepValue(rVariable);
            value += nodal_value * rShapeFunctionVector[i];
        }
    } else if (mInputVariableType == Globals::DataLocation::NodeNonHistorical) {
        for (SizeType i = 0; i < rShapeFunctionVector.size(); ++i) {
            KRATOS_DEBUG_ERROR_IF_NOT(rGeometry[i].Has(rVariable)) << "The Variable " << rVariable.Name() << " is not available at the nodes of the Geometry to retrieve Table values." << std::endl;
            const double nodal_value = rGeometry[i].GetValue(rVariable);
            value += nodal_value * rShapeFunctionVector[i];
        }
    } else if (mInputVariableType == Globals::DataLocation::Element) {
        KRATOS_DEBUG_ERROR_IF_NOT(rGeometry.Has(rVariable)) << "The Variable " << rVariable.Name() << " is not available at the Geometry to retrieve Table values." << std::endl;
        value = rGeometry.GetValue(rVariable);
    } else if (mInputVariableType == Globals::DataLocation::ProcessInfo) {
        KRATOS_DEBUG_ERROR_IF_NOT(rProcessInfo.Has(rVariable)) << "The Variable " << rVariable.Name() << " is not available at the rProcessInfo to retrieve Table values." << std::endl;
        value = rProcessInfo.GetValue(rVariable);
    } else {
        KRATOS_ERROR << "The table_input_variable_type is incorrect or not supported. Types available are : nodal_historical, nodal_non_historical and elemental_non_historical" << std::endl;
    }

    // Retrieve the dependent variable from the table
    return value;
}

/***********************************************************************************/
/***********************************************************************************/

void DatabaseAccessor::save(Serializer& rSerializer) const
{
    // we must do the int cast to be able to compile
    rSerializer.save("InputVariableType", static_cast<int>(mInputVariableType));
}
void DatabaseAccessor::load(Serializer& rSerializer)
{
    // we must do the int cast to be able to compile
    int enum_value;
    rSerializer.load("InputVariableType", enum_value);
    mInputVariableType = static_cast<Globals::DataLocation>(enum_value);
}

/***********************************************************************************/
/***********************************************************************************/

Accessor::UniquePointer DatabaseAccessor::Clone() const
{
    return Kratos::make_unique<DatabaseAccessor>(*this);
}

} // namespace Kratos
