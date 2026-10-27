//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//

#pragma once

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/table.h"
#include "includes/variables.h"

namespace Kratos
{

/**
 * @brief Historical nodal-field assignment helpers shared by the Dam nodal
 * field producers (nodal reference temperature, input-table and spatially
 * varying Young's modulus, random fields and chemo-mechanical aging).
 * @details All these producers write the historical value of a nodal variable
 * on mesh 0, either as a uniform constant or as a per-node value interpolated
 * from a model-part table using the node Id as abscissa. The loops live here so
 * that the producers keep only their parameter handling and their own
 * constitutive/aging law.
 */
class DamNodalFieldUtilities
{
public:

    using VariableType = Variable<double>;
    using TableType    = Table<double,double>;

    /**
     * @brief Writes the same historical value to every node of mesh 0.
     */
    static void SetUniformHistoricalValue(
        ModelPart& rModelPart,
        const VariableType& rVariable,
        const double Value)
    {
        auto& r_nodes = rModelPart.GetMesh(0).Nodes();
        const int number_of_nodes = static_cast<int>(r_nodes.size());

        #pragma omp parallel for
        for (int i = 0; i < number_of_nodes; ++i) {
            (r_nodes.begin() + i)->FastGetSolutionStepValue(rVariable) = Value;
        }
    }

    /**
     * @brief Writes a different historical value to every node of mesh 0,
     * interpolated from the table with the node Id as abscissa.
     */
    static void SetHistoricalValuesFromTable(
        ModelPart& rModelPart,
        const VariableType& rVariable,
        const TableType& rTable)
    {
        auto& r_nodes = rModelPart.GetMesh(0).Nodes();
        const int number_of_nodes = static_cast<int>(r_nodes.size());

        #pragma omp parallel for
        for (int i = 0; i < number_of_nodes; ++i) {
            auto it_node = r_nodes.begin() + i;
            it_node->FastGetSolutionStepValue(rVariable) = rTable.GetValue(it_node->Id());
        }
    }

    /**
     * @brief Writes a uniform value or a per-node table value, according to the
     * caller's input-file selection.
     */
    static void SetHistoricalValuesFromTableOrConstant(
        ModelPart& rModelPart,
        const VariableType& rVariable,
        const double Value,
        const TableType& rTable,
        const bool UseTable)
    {
        if (UseTable) {
            SetHistoricalValuesFromTable(rModelPart, rVariable, rTable);
        } else {
            SetUniformHistoricalValue(rModelPart, rVariable, Value);
        }
    }

}; // class DamNodalFieldUtilities

} // namespace Kratos
