// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Martin Fusseder, https://github.com/MFusseder
//

// System includes
#include <iostream>
#include <fstream>
#include <ostream>

// External includes
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

// Project includes
#include "adjoint_nodal_displacement_measurement_residual_response_function.h"
#include "processes/generic_find_elements_neighbours_process.h"
// #include "kratos_parameters.h"

namespace Kratos
{

    AdjointNodalDisplacementMeasurementResidualResponseFunction::AdjointNodalDisplacementMeasurementResidualResponseFunction(ModelPart &rModelPart, Parameters ResponseSettings)
        : AdjointStructuralResponseFunction(rModelPart, ResponseSettings)
    {
        mResponsePartName = ResponseSettings["response_part_name"].GetString();
        mResponseDirection = ResponseSettings["direction"].GetVector();
        mTracedDofLabel = ResponseSettings["traced_dof"].GetString();

        std::string mMeasurementFileName = "";
        int mMeasurementLoadCaseToUse = -1;

        try
        {
            mMeasurementFileName = ResponseSettings["measurement_file_name"].GetString();
            mMeasurementLoadCaseToUse = ResponseSettings["measurement_load_case_to_use"].GetInt();
        }
        catch (const std::exception &e)
        {
            KRATOS_ERROR_IF_NOT(false) << "\n\nAdjointNodalDisplacementMeasurementResidualResponseFunction: Properties specific to AdjointNodalDisplacementMeasurementResidualResponseFunction are not completely specified. Please provide 'measurement_file_name' (as str) and 'measurement_load_case_to_use' (as int)." << ResponseSettings << std::endl;
        }

        std::ifstream t(mMeasurementFileName);
        KRATOS_ERROR_IF_NOT(t.good() == 1) << "\n\nAdjointNodalDisplacementMeasurementResidualResponseFunction: The file " << mMeasurementFileName << " which contains the measurement data was not found. Please use the parameter option 'measurement_file_name' to specify the path.\n\n"
                                           << std::endl;

        std::stringstream buffer;
        buffer << t.rdbuf();
        auto json_string = buffer.str();
        std::cout << "::[AdjointNodalDisplacementMeasurementResidualResponseFunction]:: Start reading measurement json ..." << std::endl;

        mMeasurementData = Parameters(json_string)["load_cases"][mMeasurementLoadCaseToUse];

        if (norm_2(mResponseDirection) > 1.0e-7)
        {
            mResponseDirection /= norm_2(mResponseDirection);
        }
        else
        {
            KRATOS_ERROR << "AdjointNodalDisplacementMeasurementResidualResponseFunction: 'response_direction' must not have a norm of 1.0." << std::endl;
        }

        // Check if variable for traced dof is valid
        KRATOS_ERROR_IF_NOT(KratosComponents<ArrayVariableType>::Has(mTracedDofLabel))
            << "AdjointNodalDisplacementMeasurementResidualResponseFunction: Specified traced DOF is not available. Specified DOF: " << mTracedDofLabel << std::endl;

        // Check if variable for traced adjoint dof is valid
        KRATOS_ERROR_IF_NOT(KratosComponents<ArrayVariableType>::Has(std::string("ADJOINT_") + mTracedDofLabel))
            << "AdjointNodalDisplacementMeasurementResidualResponseFunction: Specified traced adjoint DOF is not available." << mTracedDofLabel << std::endl;

        ModelPart &response_part = rModelPart.GetSubModelPart(mResponsePartName);
        const ArrayVariableType &r_traced_dof = KratosComponents<ArrayVariableType>::Get(mTracedDofLabel);
        for (auto &node_i : response_part.Nodes())
        {
            KRATOS_ERROR_IF_NOT(node_i.SolutionStepsDataHas(r_traced_dof))
                << "AdjointNodalDisplacementMeasurementResidualResponseFunction: Specified DOF is not available at traced node." << std::endl;
        }

        // std::cout << "::[AdjointNodalDisplacementMeasurementResidualResponseFunction]:: Start computing neighboring list..." << std::endl;
        this->ComputeNeighboringElementNodeMap();
        std::cout << "::[AdjointNodalDisplacementMeasurementResidualResponseFunction]:: Finished computing neighboring list." << std::endl;
    }

    AdjointNodalDisplacementMeasurementResidualResponseFunction::~AdjointNodalDisplacementMeasurementResidualResponseFunction() {}

    void AdjointNodalDisplacementMeasurementResidualResponseFunction::CalculateGradient(const Element &rAdjointElement,
                                                                                        const Matrix &rResidualGradient,
                                                                                        Vector &rResponseGradient,
                                                                                        const ProcessInfo &rProcessInfo)
    {
        KRATOS_TRY;

        // calculate the gradient of the objective function with respect to displacement

        if (rResponseGradient.size() != rResidualGradient.size1())
            rResponseGradient.resize(rResidualGradient.size1(), false);

        rResponseGradient.clear();
        const Variable<double> *adjoint_solution_variable = &KratosComponents<Variable<double>>::Get("ADJOINT_" + mTracedDofLabel + "_X");
        DofsVectorType dofs_of_element;

        auto it_map = mElementNodeMap.find(rAdjointElement.Id());
        if (it_map == mElementNodeMap.end())
        {
            return;
        }

        double measurement_value;
        Vector measurement_normal;
        Vector simulated_displacement;
        double simulated_displacement_projected_on_measurement;

        const ModelPart &response_part = mrModelPart.GetSubModelPart(mResponsePartName);
        const ArrayVariableType &r_traced_dof = KratosComponents<ArrayVariableType>::Get(mTracedDofLabel);

        rAdjointElement.GetDofList(dofs_of_element, rProcessInfo);
        for (auto const &node_id : it_map->second) // iterate over the node ids of all nodes that are part of the element
        {
            for (IndexType i = 0; i < dofs_of_element.size(); ++i) // iterate over all available dof in the element
            {

                if (dofs_of_element[i]->Id() == node_id &&
                    dofs_of_element[i]->GetVariable() == *adjoint_solution_variable)
                {
                    for (auto &sensor_data : mMeasurementData["sensors_infos"])
                    {

                        const long unsigned int sensor_node_id = sensor_data["mesh_node_id"].GetInt();

                        if (sensor_node_id == node_id)
                        {
                            measurement_value = sensor_data["measured_value"].GetDouble();
                            measurement_normal = sensor_data["measurement_direction_normal"].GetVector();

                            const ArrayVariableType &measurement_type = KratosComponents<ArrayVariableType>::Get(sensor_data["type_of_sensor"].GetString());
                            simulated_displacement = response_part.GetNode(node_id).FastGetSolutionStepValue(measurement_type);

                            simulated_displacement_projected_on_measurement = inner_prod(simulated_displacement, measurement_normal);

                            rResponseGradient[i + 0] = measurement_normal[0] * (measurement_value - simulated_displacement_projected_on_measurement);
                            rResponseGradient[i + 1] = measurement_normal[1] * (measurement_value - simulated_displacement_projected_on_measurement);
                            rResponseGradient[i + 2] = measurement_normal[2] * (measurement_value - simulated_displacement_projected_on_measurement);

                            break;
                        }
                    }
                }
            }
        }

        KRATOS_CATCH("");
    }

    void AdjointNodalDisplacementMeasurementResidualResponseFunction::CalculateFirstDerivativesGradient(
        const Element &rAdjointElement,
        const Matrix &rResidualGradient,
        Vector &rResponseGradient,
        const ProcessInfo &rProcessInfo)
    {
        KRATOS_TRY;
        rResponseGradient = ZeroVector(rResidualGradient.size1());
        KRATOS_CATCH("");
    }

    void AdjointNodalDisplacementMeasurementResidualResponseFunction::CalculateFirstDerivativesGradient(
        const Condition &rAdjointCondition,
        const Matrix &rResidualGradient,
        Vector &rResponseGradient,
        const ProcessInfo &rProcessInfo)
    {
        KRATOS_TRY;
        rResponseGradient = ZeroVector(rResidualGradient.size1());
        KRATOS_CATCH("");
    }

    void AdjointNodalDisplacementMeasurementResidualResponseFunction::CalculateSecondDerivativesGradient(
        const Element &rAdjointElement,
        const Matrix &rResidualGradient,
        Vector &rResponseGradient,
        const ProcessInfo &rProcessInfo)
    {
        KRATOS_TRY;
        rResponseGradient = ZeroVector(rResidualGradient.size1());
        KRATOS_CATCH("");
    }

    void AdjointNodalDisplacementMeasurementResidualResponseFunction::CalculateSecondDerivativesGradient(
        const Condition &rAdjointCondition,
        const Matrix &rResidualGradient,
        Vector &rResponseGradient,
        const ProcessInfo &rProcessInfo)
    {
        KRATOS_TRY;
        rResponseGradient = ZeroVector(rResidualGradient.size1());
        KRATOS_CATCH("");
    }

    void AdjointNodalDisplacementMeasurementResidualResponseFunction::CalculatePartialSensitivity(Element &rAdjointElement,
                                                                                                  const Variable<double> &rVariable,
                                                                                                  const Matrix &rSensitivityMatrix,
                                                                                                  Vector &rSensitivityGradient,
                                                                                                  const ProcessInfo &rProcessInfo)
    {
        KRATOS_TRY;
        rSensitivityGradient = ZeroVector(rSensitivityMatrix.size1());
        KRATOS_CATCH("");
    }

    void AdjointNodalDisplacementMeasurementResidualResponseFunction::CalculatePartialSensitivity(Condition &rAdjointCondition,
                                                                                                  const Variable<double> &rVariable,
                                                                                                  const Matrix &rSensitivityMatrix,
                                                                                                  Vector &rSensitivityGradient,
                                                                                                  const ProcessInfo &rProcessInfo)
    {
        KRATOS_TRY;
        rSensitivityGradient = ZeroVector(rSensitivityMatrix.size1());
        KRATOS_CATCH("");
    }

    void AdjointNodalDisplacementMeasurementResidualResponseFunction::CalculatePartialSensitivity(Element &rAdjointElement,
                                                                                                  const Variable<array_1d<double, 3>> &rVariable,
                                                                                                  const Matrix &rSensitivityMatrix,
                                                                                                  Vector &rSensitivityGradient,
                                                                                                  const ProcessInfo &rProcessInfo)
    {
        KRATOS_TRY;
        rSensitivityGradient = ZeroVector(rSensitivityMatrix.size1());
        KRATOS_CATCH("");
    }

    void AdjointNodalDisplacementMeasurementResidualResponseFunction::CalculatePartialSensitivity(Condition &rAdjointCondition,
                                                                                                  const Variable<array_1d<double, 3>> &rVariable,
                                                                                                  const Matrix &rSensitivityMatrix,
                                                                                                  Vector &rSensitivityGradient,
                                                                                                  const ProcessInfo &rProcessInfo)
    {
        KRATOS_TRY;
        rSensitivityGradient = ZeroVector(rSensitivityMatrix.size1());
        KRATOS_CATCH("");
    }

    double AdjointNodalDisplacementMeasurementResidualResponseFunction::CalculateValue(ModelPart &rModelPart)
    {
        KRATOS_TRY;

        const ArrayVariableType &r_traced_dof = KratosComponents<ArrayVariableType>::Get(mTracedDofLabel);

        double response_value = 0.0;
        double measurement_value;
        Vector measurement_normal;
        int sensor_node_id;
        Vector simulated_response_vector;
        double projected_simulated_response;
        double measurement_difference;
        ModelPart &response_part = rModelPart.GetSubModelPart(mResponsePartName);

        for (auto &sensor_data : mMeasurementData["sensors_infos"])
        {
            measurement_value = sensor_data["measured_value"].GetDouble();
            measurement_normal = sensor_data["measurement_direction_normal"].GetVector();

            sensor_node_id = sensor_data["mesh_node_id"].GetInt();

            const ArrayVariableType &measurement_type = KratosComponents<ArrayVariableType>::Get(sensor_data["type_of_sensor"].GetString());
            simulated_response_vector = response_part.GetNode(sensor_node_id).FastGetSolutionStepValue(measurement_type);

            projected_simulated_response = inner_prod(measurement_normal, simulated_response_vector);

            measurement_difference = measurement_value - projected_simulated_response;
            response_value += measurement_difference * measurement_difference;
        }
        response_value *= 0.5;

        return response_value;

        KRATOS_CATCH("");
    }

    /// Find one element which is bounded by one of the traced nodes. The elements are needed for assembling the adjoint load.
    void AdjointNodalDisplacementMeasurementResidualResponseFunction::ComputeNeighboringElementNodeMap()
    {
        KRATOS_TRY;

        ModelPart &response_part = mrModelPart.GetSubModelPart(mResponsePartName);
        GenericFindElementalNeighboursProcess neighbour_elements_finder(mrModelPart);
        neighbour_elements_finder.Execute();

        for (auto &node_i : response_part.Nodes())
        {
            auto const &r_neighbours = node_i.GetValue(NEIGHBOUR_ELEMENTS);
            KRATOS_ERROR_IF(r_neighbours.size() == 0) << "AdjointNodalDisplacementMeasurementResidualResponseFunction: Node " << node_i.Id() << " has no neighbouring element" << std::endl;
            // take the first element since only one neighbour element is required
            auto it_map = mElementNodeMap.find(r_neighbours[0].Id());
            if (it_map == mElementNodeMap.end())
            {
                std::vector<IndexType> node_ids = {node_i.Id()};
                mElementNodeMap[r_neighbours[0].Id()] = node_ids;
            }
            else
            {
                (it_map->second).push_back(node_i.Id());
            }
        }

        KRATOS_CATCH("");
    }

} // namespace Kratos.
