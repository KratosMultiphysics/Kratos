// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Jonas Hillenbrand, M.Sc.
//

// System includes

// External includes

// Project includes
#include "adjoint_least_squares_response_function.h"
#include "processes/find_elements_neighbours_process.h"

namespace Kratos
{
    AdjointLeastSquaresResponseFunction::AdjointLeastSquaresResponseFunction(ModelPart& rModelPart, Parameters ResponseSettings, Parameters MeasurementData)
    : AdjointStructuralResponseFunction(rModelPart, ResponseSettings)
    { 
        mMeasurementData = MeasurementData;
        mResponsePartName = ResponseSettings["response_part_name"].GetString();

        // Assuming constant measurement time step
        mMeasurementDeltaTime = MeasurementData["TIME"][1].GetDouble() - MeasurementData["TIME"][0].GetDouble();
        mMeasurementTimeSteps = (MeasurementData["TIME"].GetVector()).size();

        // Check if measurement data for all nodes in response model part is available
        ModelPart& response_part = rModelPart.GetSubModelPart(mResponsePartName);
        for(auto& node_i : response_part.Nodes()){
            KRATOS_ERROR_IF_NOT( mMeasurementData.Has("NODE_" + std::to_string(node_i.Id())) )
                << "AdjointNodalDisplacementResponseFunction: Node" << node_i.Id() << " has no measurement data." << std::endl;
        }

        this->ComputeNeighboringElementNodeMap();
    }

    AdjointLeastSquaresResponseFunction::~AdjointLeastSquaresResponseFunction(){}

    void AdjointLeastSquaresResponseFunction::CalculateGradient(const Element& rAdjointElement,
                                   const Matrix& rResidualGradient,
                                   Vector& rResponseGradient,
                                   const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;

        if (rResponseGradient.size() != rResidualGradient.size1())
            rResponseGradient.resize(rResidualGradient.size1(), false);

        rResponseGradient.clear();
        const Variable<double>& r_x_dof = KratosComponents<Variable<double>>::Get("DISPLACEMENT_X");
        const Variable<double>& r_y_dof = KratosComponents<Variable<double>>::Get("DISPLACEMENT_Y");
        const Variable<double>& r_z_dof = KratosComponents<Variable<double>>::Get("DISPLACEMENT_Z");
        
        const Variable<double>& r_time = KratosComponents<Variable<double>>::Get("TIME");
        double simulation_time = rProcessInfo.GetValue(r_time);

        Vector u = ZeroVector(3);
        Vector u_m = ZeroVector(3);

        const Variable<double>& r_adjoint_x_dof = KratosComponents<Variable<double>>::Get("ADJOINT_DISPLACEMENT_X");
        DofsVectorType dofs_of_element;
        auto it_map = mElementNodeMap.find(rAdjointElement.Id());
        if (it_map != mElementNodeMap.end()) {
            rAdjointElement.GetDofList(dofs_of_element, rProcessInfo);
            for(auto const& node_id: it_map->second) {
                for(IndexType i = 0; i < dofs_of_element.size(); ++i) {
                    if (dofs_of_element[i]->Id() == node_id &&
                        dofs_of_element[i]->GetVariable() == r_adjoint_x_dof) { 
                        u[0] = mrModelPart.GetNode(node_id).FastGetSolutionStepValue(r_x_dof, 0);
                        u[1] = mrModelPart.GetNode(node_id).FastGetSolutionStepValue(r_y_dof, 0);
                        u[2] = mrModelPart.GetNode(node_id).FastGetSolutionStepValue(r_z_dof, 0);
                        this->InterpolateNodalMeasurementData(u_m, node_id, simulation_time); 
                        rResponseGradient[i] = u[0] - u_m[0];
                        rResponseGradient[i+1] = u[1] - u_m[1];
                        rResponseGradient[i+2] = u[2] - u_m[2];
                        break;
                    }
                }
            }
        }

        KRATOS_CATCH("");
    }

    void AdjointLeastSquaresResponseFunction::CalculateFirstDerivativesGradient(
        const Element& rAdjointElement,
        const Matrix& rResidualGradient,
        Vector& rResponseGradient,
        const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;

        rResponseGradient = ZeroVector(rResidualGradient.size1());
        
        KRATOS_CATCH("");
    }

    void AdjointLeastSquaresResponseFunction::CalculateFirstDerivativesGradient(
        const Condition& rAdjointCondition,
        const Matrix& rResidualGradient,
        Vector& rResponseGradient,
        const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;
        rResponseGradient = ZeroVector(rResidualGradient.size1());
        KRATOS_CATCH("");
    }

    void AdjointLeastSquaresResponseFunction::CalculateSecondDerivativesGradient(
        const Element& rAdjointElement,
        const Matrix& rResidualGradient,
        Vector& rResponseGradient,
        const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;

        rResponseGradient = ZeroVector(rResidualGradient.size1());
        
        KRATOS_CATCH("");
    }

    void AdjointLeastSquaresResponseFunction::CalculateSecondDerivativesGradient(
        const Condition& rAdjointCondition,
        const Matrix& rResidualGradient,
        Vector& rResponseGradient,
        const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;
        rResponseGradient = ZeroVector(rResidualGradient.size1());
        KRATOS_CATCH("");
    }

    void AdjointLeastSquaresResponseFunction::CalculatePartialSensitivity(Element& rAdjointElement,
                                             const Variable<double>& rVariable,
                                             const Matrix& rSensitivityMatrix,
                                             Vector& rSensitivityGradient,
                                             const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;
        rSensitivityGradient = ZeroVector(rSensitivityMatrix.size1());
        KRATOS_CATCH("");
    }

    void AdjointLeastSquaresResponseFunction::CalculatePartialSensitivity(Condition& rAdjointCondition,
                                             const Variable<double>& rVariable,
                                             const Matrix& rSensitivityMatrix,
                                             Vector& rSensitivityGradient,
                                             const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;
        rSensitivityGradient = ZeroVector(rSensitivityMatrix.size1());
        KRATOS_CATCH("");
    }

    void AdjointLeastSquaresResponseFunction::CalculatePartialSensitivity(Element& rAdjointElement,
                                             const Variable<array_1d<double, 3>>& rVariable,
                                             const Matrix& rSensitivityMatrix,
                                             Vector& rSensitivityGradient,
                                             const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;
        rSensitivityGradient = ZeroVector(rSensitivityMatrix.size1());
        KRATOS_CATCH("");
    }

    void AdjointLeastSquaresResponseFunction::CalculatePartialSensitivity(Condition& rAdjointCondition,
                                             const Variable<array_1d<double, 3>>& rVariable,
                                             const Matrix& rSensitivityMatrix,
                                             Vector& rSensitivityGradient,
                                             const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;
        rSensitivityGradient = ZeroVector(rSensitivityMatrix.size1());
        KRATOS_CATCH("");
    }

    double AdjointLeastSquaresResponseFunction::CalculateValue(ModelPart& rModelPart)
    {
        KRATOS_TRY;
        
        const Variable<double>& r_time = KratosComponents<Variable<double>>::Get("TIME");
        double simulation_time = rModelPart.GetProcessInfo().GetValue(r_time);

        const Variable<double>& r_x_dof = KratosComponents<Variable<double>>::Get("DISPLACEMENT_X");
        const Variable<double>& r_y_dof = KratosComponents<Variable<double>>::Get("DISPLACEMENT_Y");
        const Variable<double>& r_z_dof = KratosComponents<Variable<double>>::Get("DISPLACEMENT_Z");
        ModelPart& response_part = rModelPart.GetSubModelPart(mResponsePartName);

        double value = 0;
        Vector u = ZeroVector(3);
        Vector u_m = ZeroVector(3);
        for(auto& node_i : response_part.Nodes()){
            u[0] = rModelPart.GetNode(node_i.Id()).FastGetSolutionStepValue(r_x_dof , 0);
            u[1] = rModelPart.GetNode(node_i.Id()).FastGetSolutionStepValue(r_y_dof , 0);
            u[2] = rModelPart.GetNode(node_i.Id()).FastGetSolutionStepValue(r_z_dof , 0);
            this->InterpolateNodalMeasurementData(u_m, node_i.Id(), simulation_time);
            value += inner_prod((u - u_m), (u - u_m));
        }
        return 0.5 * value;
        
        KRATOS_CATCH("");
    }

    /// Find one element which is bounded by the traced node. The element is needed for assembling the adjoint load.
    void AdjointLeastSquaresResponseFunction::ComputeNeighboringElementNodeMap()
    {
        KRATOS_TRY;

        ModelPart& response_part = mrModelPart.GetSubModelPart(mResponsePartName);
        FindElementalNeighboursProcess neighbour_elements_finder(mrModelPart, 10, 10);
        neighbour_elements_finder.Execute();

        for(auto& node_i : response_part.Nodes()) {
            auto const& r_neighbours = node_i.GetValue(NEIGHBOUR_ELEMENTS);
            KRATOS_ERROR_IF(r_neighbours.size() == 0) << "AdjointNodalDisplacementResponseFunction: Node " << node_i.Id() << " has no neighbouring element" << std::endl;
            // take the first element since only one neighbour element is required
            auto it_map = mElementNodeMap.find(r_neighbours[0].Id());
            if (it_map == mElementNodeMap.end()) {
                std::vector<IndexType> node_ids = {node_i.Id()};
                mElementNodeMap[r_neighbours[0].Id()] = node_ids;
            }
            else {
                (it_map->second).push_back(node_i.Id());
            }
        }

        KRATOS_CATCH("");
    }

    void AdjointLeastSquaresResponseFunction::InterpolateNodalMeasurementData(Vector& rNodalMeasuredDisplacement,
                                                                              const IndexType& rNodeId,
                                                                              const double& rSimulationTime)
    {
        double step = rSimulationTime / mMeasurementDeltaTime;
        IndexType lower_step = int(step);
        double lower_weight = int(step) + 1 - step;
        if (lower_step > 0 && lower_step < mMeasurementTimeSteps){
            IndexType upper_step = int(step) + 1;
            double upper_weight = step - int(step);

            rNodalMeasuredDisplacement[0] = lower_weight * mMeasurementData["NODE_"+ std::to_string(rNodeId)]["DISPLACEMENT_X"][lower_step-1].GetDouble() +
                                            upper_weight * mMeasurementData["NODE_"+ std::to_string(rNodeId)]["DISPLACEMENT_X"][upper_step-1].GetDouble();
            rNodalMeasuredDisplacement[1] = lower_weight * mMeasurementData["NODE_"+ std::to_string(rNodeId)]["DISPLACEMENT_Y"][lower_step-1].GetDouble() +
                                            upper_weight * mMeasurementData["NODE_"+ std::to_string(rNodeId)]["DISPLACEMENT_Y"][upper_step-1].GetDouble();
            rNodalMeasuredDisplacement[2] = lower_weight * mMeasurementData["NODE_"+ std::to_string(rNodeId)]["DISPLACEMENT_Z"][lower_step-1].GetDouble() +
                                            upper_weight * mMeasurementData["NODE_"+ std::to_string(rNodeId)]["DISPLACEMENT_Z"][upper_step-1].GetDouble();
        }
        else if (lower_step == 0){
            rNodalMeasuredDisplacement[0] = mMeasurementData["NODE_"+ std::to_string(rNodeId)]["DISPLACEMENT_X"][lower_step].GetDouble();
            rNodalMeasuredDisplacement[1] = mMeasurementData["NODE_"+ std::to_string(rNodeId)]["DISPLACEMENT_Y"][lower_step].GetDouble();
            rNodalMeasuredDisplacement[2] = mMeasurementData["NODE_"+ std::to_string(rNodeId)]["DISPLACEMENT_Z"][lower_step].GetDouble();
        }
        else{
            rNodalMeasuredDisplacement[0] = mMeasurementData["NODE_"+ std::to_string(rNodeId)]["DISPLACEMENT_X"][lower_step-1].GetDouble();
            rNodalMeasuredDisplacement[1] = mMeasurementData["NODE_"+ std::to_string(rNodeId)]["DISPLACEMENT_Y"][lower_step-1].GetDouble();
            rNodalMeasuredDisplacement[2] = mMeasurementData["NODE_"+ std::to_string(rNodeId)]["DISPLACEMENT_Z"][lower_step-1].GetDouble();
        }

    }

    
} // namespace Kratos.


