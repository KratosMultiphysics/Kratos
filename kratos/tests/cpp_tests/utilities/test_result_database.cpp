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
//

// System includes

// External includes

// Project includes
#include "containers/model.h"
#include "testing/testing.h"
#include "utilities/cpp_tests_utilities.h"
#include "utilities/result_dabatase.h"

namespace Kratos {
namespace Testing {

KRATOS_TEST_CASE_IN_SUITE(ResultDatabase, KratosCoreFastSuite)
{
    // Create model part
    Model current_model;

    ModelPart& r_model_part = current_model.CreateModelPart("Main");
    r_model_part.SetBufferSize(2);

    r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
    r_model_part.AddNodalSolutionStepVariable(TEMPERATURE);

    auto& r_process_info = r_model_part.GetProcessInfo();
    r_process_info.SetValue(DOMAIN_SIZE, 2);
    r_process_info.SetValue(NL_ITERATION_NUMBER, 1);

    CppTestsUtilities::Create2DGeometry(r_model_part);

    // Fill values
    const auto& r_nodes_array = r_model_part.Nodes();
    const auto it_node_begin = r_nodes_array.begin();

    for (int i = 0; i < 2; ++i) {
        const double time = i * 1.0;
        r_process_info.SetValue(STEP, i + 1);
        r_model_part.CloneTimeStep(time);
        Vector aux_vector = ZeroVector(3);
        for (auto& r_node : r_nodes_array) {
            aux_vector[0] = time * r_node.X()-1.0;
            aux_vector[1] = time * r_node.Y()+2.0;
            r_node.FastGetSolutionStepValue(DISPLACEMENT) = aux_vector;
            r_node.FastGetSolutionStepValue(TEMPERATURE) = time * (std::pow(r_node.X(),2)+r_node.Y());
        }
    }

    // Fill database
    ResultDatabase database_nodes;
    std::vector<IndexType> nodal_variables_ids(4);
    std::vector<IndexType> nodal_values_sizes(4, 1);

    const std::vector<Variable<double>*> variable_vector = {&DISPLACEMENT_X,&DISPLACEMENT_Y,&DISPLACEMENT_Z,&TEMPERATURE};

    for (std::size_t i = 0; i < 4; ++i) {
        nodal_variables_ids[i] = variable_vector[i]->Key();
    }

    database_nodes.Initialize(nodal_variables_ids, nodal_values_sizes, r_nodes_array.size());

    Vector time = ZeroVector(2);
    time[1] = 1.0;
    database_nodes.SetCommonColumn(time);


    for (auto& p_var_double : variable_vector) {
        auto& r_var_database = database_nodes.GetVariableData(*p_var_double);
        for (int i = 0; i < static_cast<int>(r_nodes_array.size()); ++i) {
            auto it_node = it_node_begin + i;
            Vector values_vector = ZeroVector(2);
            for (int j = 0; j < 2; ++j) {
                values_vector[j] = it_node->FastGetSolutionStepValue(*p_var_double);
            }
            r_var_database.SetValues(time, values_vector, i);
        }
    }

    // Check database
    for (int i = 0; i < 2; ++i) {
        const double time = i * 1.0;

        for (auto& p_var_double : variable_vector) {
            const auto& r_var_database = database_nodes.GetVariableData(*p_var_double);

            for (int i = 0; i < static_cast<int>(r_nodes_array.size()); ++i) {
                auto it_node = it_node_begin + i;

                const double result = it_node->FastGetSolutionStepValue(*p_var_double);
                const double reference = r_var_database.GetValue(i, time);

                KRATOS_CHECK_NEAR(result, reference, 1.0e-12);
            }
        }
    }

}

}   // namespace Testing
}  // namespace Kratos.
