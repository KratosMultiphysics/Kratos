// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Martin Fusseder, https://github.com/MFusseder
//

// System includes

// External includes

// Project includes
#include "adjoint_nodal_displacement_aggregation_response_function.h"
#include "processes/find_elements_neighbours_process.h"

namespace Kratos
{
    AdjointNodalDisplacementAggregationResponseFunction::AdjointNodalDisplacementAggregationResponseFunction(ModelPart& rModelPart, Parameters ResponseSettings)
    : AdjointStructuralResponseFunction(rModelPart, ResponseSettings)
    {
        mResponsePartName = ResponseSettings["response_part_name"].GetString();
        mResponseDirection = ResponseSettings["direction"].GetVector();
        mTracedDofLabel = ResponseSettings["traced_dof"].GetString();
        mObservedResponseVariableName = ResponseSettings["observed_response_variable"].GetString();
        mRho = ResponseSettings["rho"].GetDouble();
        mDeformationScalingFactor = ResponseSettings["deformation_scaling_factor"].GetDouble();

        if ( norm_2( mResponseDirection ) > 1.0e-7 ) {
            mResponseDirection /= norm_2( mResponseDirection );
        } else {
            KRATOS_ERROR << "AdjointNodalDisplacementAggregationResponseFunction: 'response_direction' must not have a norm of 0.0." << std::endl;
        }

        // Check if variable for traced dof is valid
        KRATOS_ERROR_IF_NOT( KratosComponents<ArrayVariableType>::Has(mTracedDofLabel) )
            << "AdjointNodalDisplacementAggregationResponseFunction: Specified traced DOF is not available. Specified DOF: " << mTracedDofLabel << std::endl;

        // Check if variable for traced adjoint dof is valid
        KRATOS_ERROR_IF_NOT( KratosComponents<ArrayVariableType>::Has(std::string("ADJOINT_") + mTracedDofLabel) )
            << "AdjointNodalDisplacementAggregationResponseFunction: Specified traced adjoint DOF is not available." << mTracedDofLabel << std::endl;

        // Check if variable for traced response variable is valid
        KRATOS_ERROR_IF_NOT( KratosComponents<ArrayVariableType>::Has(mObservedResponseVariableName) )
            << "AdjointNodalDisplacementAggregationResponseFunction: Specified traced response variable is not available. Specified: " << mObservedResponseVariableName << std::endl;

        ModelPart& response_part = rModelPart.GetSubModelPart(mResponsePartName);
        const ArrayVariableType& r_traced_dof = KratosComponents<ArrayVariableType>::Get(mTracedDofLabel);
        for(auto& node_i : response_part.Nodes()){
            KRATOS_ERROR_IF_NOT( node_i.SolutionStepsDataHas(r_traced_dof) )
                << "AdjointNodalDisplacementAggregationResponseFunction: Specified DOF is not available at traced node." << std::endl;
        }

        // Initialize
        for(auto& node_i : response_part.Nodes()) {
            mKSPrefactors[node_i.Id()] = 0.0;
        }

        this->ComputeNeighboringElementNodeMap();
    }

    AdjointNodalDisplacementAggregationResponseFunction::~AdjointNodalDisplacementAggregationResponseFunction(){}

    void AdjointNodalDisplacementAggregationResponseFunction::CalculateGradient(const Element& rAdjointElement,
                                   const Matrix& rResidualGradient,
                                   Vector& rResponseGradient,
                                   const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;

        if (rResponseGradient.size() != rResidualGradient.size1())
            rResponseGradient.resize(rResidualGradient.size1(), false);

        rResponseGradient.clear();
        const Variable<double>* adjoint_solution_variable = &KratosComponents<Variable<double>>::Get("ADJOINT_" + mTracedDofLabel + "_X");
        DofsVectorType dofs_of_element;

        auto it_map = mElementNodeMap.find(rAdjointElement.Id());
        if (it_map != mElementNodeMap.end()) {
            rAdjointElement.GetDofList(dofs_of_element, rProcessInfo);
            for(auto const& node_id: it_map->second) {
                for(IndexType i = 0; i < dofs_of_element.size(); ++i) {
                    if (dofs_of_element[i]->Id() == node_id &&
                        dofs_of_element[i]->GetVariable() == *adjoint_solution_variable) {
                        rResponseGradient[i]   = -1 * mResponseDirection[0] * mKSPrefactors[node_id] / (mDeformationScalingFactor * mSumPrefactors);
                        rResponseGradient[i+1] = -1 * mResponseDirection[1] * mKSPrefactors[node_id] / (mDeformationScalingFactor * mSumPrefactors);
                        rResponseGradient[i+2] = -1 * mResponseDirection[2] * mKSPrefactors[node_id] / (mDeformationScalingFactor * mSumPrefactors);
                        break;
                    }
                }
            }
        }

        KRATOS_CATCH("");
    }

    void AdjointNodalDisplacementAggregationResponseFunction::CalculateFirstDerivativesGradient(
        const Element& rAdjointElement,
        const Matrix& rResidualGradient,
        Vector& rResponseGradient,
        const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;
        rResponseGradient = ZeroVector(rResidualGradient.size1());
        KRATOS_CATCH("");
    }

    void AdjointNodalDisplacementAggregationResponseFunction::CalculateFirstDerivativesGradient(
        const Condition& rAdjointCondition,
        const Matrix& rResidualGradient,
        Vector& rResponseGradient,
        const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;
        rResponseGradient = ZeroVector(rResidualGradient.size1());
        KRATOS_CATCH("");
    }

    void AdjointNodalDisplacementAggregationResponseFunction::CalculateSecondDerivativesGradient(
        const Element& rAdjointElement,
        const Matrix& rResidualGradient,
        Vector& rResponseGradient,
        const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;
        rResponseGradient = ZeroVector(rResidualGradient.size1());
        KRATOS_CATCH("");
    }

    void AdjointNodalDisplacementAggregationResponseFunction::CalculateSecondDerivativesGradient(
        const Condition& rAdjointCondition,
        const Matrix& rResidualGradient,
        Vector& rResponseGradient,
        const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;
        rResponseGradient = ZeroVector(rResidualGradient.size1());
        KRATOS_CATCH("");
    }

    void AdjointNodalDisplacementAggregationResponseFunction::CalculatePartialSensitivity(Element& rAdjointElement,
                                             const Variable<double>& rVariable,
                                             const Matrix& rSensitivityMatrix,
                                             Vector& rSensitivityGradient,
                                             const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;
        rSensitivityGradient = ZeroVector(rSensitivityMatrix.size1());
        KRATOS_CATCH("");
    }

    void AdjointNodalDisplacementAggregationResponseFunction::CalculatePartialSensitivity(Condition& rAdjointCondition,
                                             const Variable<double>& rVariable,
                                             const Matrix& rSensitivityMatrix,
                                             Vector& rSensitivityGradient,
                                             const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;
        rSensitivityGradient = ZeroVector(rSensitivityMatrix.size1());
        KRATOS_CATCH("");
    }

    void AdjointNodalDisplacementAggregationResponseFunction::CalculatePartialSensitivity(Element& rAdjointElement,
                                             const Variable<array_1d<double, 3>>& rVariable,
                                             const Matrix& rSensitivityMatrix,
                                             Vector& rSensitivityGradient,
                                             const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;
        rSensitivityGradient = ZeroVector(rSensitivityMatrix.size1());
        KRATOS_CATCH("");
    }

    void AdjointNodalDisplacementAggregationResponseFunction::CalculatePartialSensitivity(Condition& rAdjointCondition,
                                             const Variable<array_1d<double, 3>>& rVariable,
                                             const Matrix& rSensitivityMatrix,
                                             Vector& rSensitivityGradient,
                                             const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;
        rSensitivityGradient = ZeroVector(rSensitivityMatrix.size1());
        KRATOS_CATCH("");
    }

    double AdjointNodalDisplacementAggregationResponseFunction::CalculateValue(ModelPart& rModelPart)
    {
        KRATOS_TRY;

        ModelPart& response_part = rModelPart.GetSubModelPart(mResponsePartName);

        // Find max value of deformation stiffness in primal model part
        std::vector<double> deformation_stiffnesses;
        deformation_stiffnesses.reserve(response_part.NumberOfNodes());

        double max_deformation_stiffness = 0.0;
        IndexType node_id_at_max = 0;

        const ArrayVariableType& r_traced_response_var = KratosComponents<ArrayVariableType>::Get(mObservedResponseVariableName);
        for(auto& node_i : response_part.Nodes()){
            KRATOS_ERROR_IF_NOT( node_i.Has(r_traced_response_var) )
                << "AdjointNodalDisplacementAggregationResponseFunction: Specified traced response variable is not available at traced node." << std::endl;
        }

        // Determine max value
        for(auto& node_i : response_part.Nodes()) {
            double deformation_stiffness = inner_prod(mResponseDirection, node_i.GetValue(r_traced_response_var));
            deformation_stiffnesses.push_back(deformation_stiffness);

            if(deformation_stiffness > max_deformation_stiffness) {
                max_deformation_stiffness = deformation_stiffness;
                node_id_at_max = node_i.Id();
            }
        }

        KRATOS_INFO_IF("  AdjointNodalDisplacementAggregationResponseFunction::CalculateValue", mEchoLevel > 0) << "Id of element with max deformation stiffness = " << node_id_at_max << std::endl;
        KRATOS_INFO_IF("  AdjointNodalDisplacementAggregationResponseFunction::CalculateValue", mEchoLevel > 0) << "Max  deformation stiffness value = " << max_deformation_stiffness << std::endl;

        /*if(std::abs(mDeformationScalingFactor) < 1e-13)
            mDeformationScalingFactor = max_deformation_stiffness;*/

        // Aggegration according Kreiselmeier-Steinhauser Function
        int node_index = 0;
        mSumPrefactors = 0.0;

        for(auto& node_i : response_part.Nodes()){
            double KS_value_contribution = std::exp(mRho * deformation_stiffnesses[node_index] / mDeformationScalingFactor);
            mKSPrefactors[node_i.Id()] = KS_value_contribution;
            mSumPrefactors += KS_value_contribution;
            node_index++;
        }
        mArePrefactorsInitialized = true;

        double KS_value = 1/mRho * std::log(mSumPrefactors);
        return KS_value;

        KRATOS_CATCH("");
    }

    /// Find one element which is bounded by one of the traced nodes. The elements are needed for assembling the adjoint load.
    void AdjointNodalDisplacementAggregationResponseFunction::ComputeNeighboringElementNodeMap()
    {
        KRATOS_TRY;

        ModelPart& response_part = mrModelPart.GetSubModelPart(mResponsePartName);
        FindElementalNeighboursProcess neighbour_elements_finder(mrModelPart, 10, 10);
        neighbour_elements_finder.Execute();

        for(auto& node_i : response_part.Nodes()) {
            auto const& r_neighbours = node_i.GetValue(NEIGHBOUR_ELEMENTS);
            KRATOS_ERROR_IF(r_neighbours.size() == 0) << "AdjointNodalDisplacementAggregationResponseFunction: Node " << node_i.Id() << " has no neighbouring element" << std::endl;
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

} // namespace Kratos.


