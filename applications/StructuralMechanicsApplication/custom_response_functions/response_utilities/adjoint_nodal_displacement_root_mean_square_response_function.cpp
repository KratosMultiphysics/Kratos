// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    
//

// System includes

// External includes

// Project includes
#include "adjoint_nodal_displacement_root_mean_square_response_function.h"
#include "processes/find_elements_neighbours_process.h"

namespace Kratos
{
    AdjointNodalDisplacementRootMeanSquareResponseFunction::AdjointNodalDisplacementRootMeanSquareResponseFunction(ModelPart& rModelPart, Parameters ResponseSettings)
    : AdjointStructuralResponseFunction(rModelPart, ResponseSettings)
    { 
        mTimeDomain = ResponseSettings["time_domain"].GetDouble();
        mResponsePartName = ResponseSettings["response_part_name"].GetString();
        mTracedDofLabel = ResponseSettings["traced_dof"].GetString();

        // Check if variable for traced dof is valid
        KRATOS_ERROR_IF_NOT( KratosComponents<Variable<double>>::Has(mTracedDofLabel) )
            << "AdjointNodalDisplacementResponseFunction: Specified traced DOF is not available. Specified DOF: " << mTracedDofLabel << std::endl;

        // Check if variable for traced adjoint dof is valid
        KRATOS_ERROR_IF_NOT( KratosComponents<Variable<double>>::Has(std::string("ADJOINT_") + mTracedDofLabel) )
            << "AdjointNodalDisplacementResponseFunction: Specified traced adjoint DOF is not available." << mTracedDofLabel << std::endl;

        ModelPart& response_part = rModelPart.GetSubModelPart(mResponsePartName);
        const Variable<double>* r_traced_dof = &KratosComponents<Variable<double>>::Get(mTracedDofLabel);
        for(auto& node_i : response_part.Nodes()){
            KRATOS_ERROR_IF_NOT( node_i.SolutionStepsDataHas(*r_traced_dof) )
                << "AdjointNodalDisplacementResponseFunction: Specified DOF is not available at traced node." << std::endl;
        }

        this->ComputeNeighboringElementNodeMap();
    }

    AdjointNodalDisplacementRootMeanSquareResponseFunction::~AdjointNodalDisplacementRootMeanSquareResponseFunction(){}

    void AdjointNodalDisplacementRootMeanSquareResponseFunction::CalculateGradient(const Element& rAdjointElement,
                                   const Matrix& rResidualGradient,
                                   Vector& rResponseGradient,
                                   const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;

        if (rResponseGradient.size() != rResidualGradient.size1())
            rResponseGradient.resize(rResidualGradient.size1(), false);

        rResponseGradient.clear();
        const Variable<double>* p_traced_dof = &KratosComponents<Variable<double>>::Get(mTracedDofLabel);
        const Variable<double>* p_traced_adjoint_dof = &KratosComponents<Variable<double>>::Get("ADJOINT_" + mTracedDofLabel);
        DofsVectorType dofs_of_element;
        
        auto it_map = mElementNodeMap.find(rAdjointElement.Id());
        if (it_map != mElementNodeMap.end()) {
            rAdjointElement.GetDofList(dofs_of_element, rProcessInfo);
            for(auto const& node_id: it_map->second) {
                for(IndexType i = 0; i < dofs_of_element.size(); ++i) {
                    if (dofs_of_element[i]->Id() == node_id &&
                        dofs_of_element[i]->GetVariable() == *p_traced_adjoint_dof) {
                        rResponseGradient[i]   = 2 / mTimeDomain * mrModelPart.GetNode(node_id).FastGetSolutionStepValue(*p_traced_dof, 0);
                        break;
                    }
                }
            }
        }

        KRATOS_CATCH("");
    }

    void AdjointNodalDisplacementRootMeanSquareResponseFunction::CalculateFirstDerivativesGradient(
        const Element& rAdjointElement,
        const Matrix& rResidualGradient,
        Vector& rResponseGradient,
        const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;
        rResponseGradient = ZeroVector(rResidualGradient.size1());
        KRATOS_CATCH("");
    }

    void AdjointNodalDisplacementRootMeanSquareResponseFunction::CalculateFirstDerivativesGradient(
        const Condition& rAdjointCondition,
        const Matrix& rResidualGradient,
        Vector& rResponseGradient,
        const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;
        rResponseGradient = ZeroVector(rResidualGradient.size1());
        KRATOS_CATCH("");
    }

    void AdjointNodalDisplacementRootMeanSquareResponseFunction::CalculateSecondDerivativesGradient(
        const Element& rAdjointElement,
        const Matrix& rResidualGradient,
        Vector& rResponseGradient,
        const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;
        rResponseGradient = ZeroVector(rResidualGradient.size1());
        KRATOS_CATCH("");
    }

    void AdjointNodalDisplacementRootMeanSquareResponseFunction::CalculateSecondDerivativesGradient(
        const Condition& rAdjointCondition,
        const Matrix& rResidualGradient,
        Vector& rResponseGradient,
        const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;
        rResponseGradient = ZeroVector(rResidualGradient.size1());
        KRATOS_CATCH("");
    }

    void AdjointNodalDisplacementRootMeanSquareResponseFunction::CalculatePartialSensitivity(Element& rAdjointElement,
                                             const Variable<double>& rVariable,
                                             const Matrix& rSensitivityMatrix,
                                             Vector& rSensitivityGradient,
                                             const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;
        rSensitivityGradient = ZeroVector(rSensitivityMatrix.size1());
        KRATOS_CATCH("");
    }

    void AdjointNodalDisplacementRootMeanSquareResponseFunction::CalculatePartialSensitivity(Condition& rAdjointCondition,
                                             const Variable<double>& rVariable,
                                             const Matrix& rSensitivityMatrix,
                                             Vector& rSensitivityGradient,
                                             const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;
        rSensitivityGradient = ZeroVector(rSensitivityMatrix.size1());
        KRATOS_CATCH("");
    }

    void AdjointNodalDisplacementRootMeanSquareResponseFunction::CalculatePartialSensitivity(Element& rAdjointElement,
                                             const Variable<array_1d<double, 3>>& rVariable,
                                             const Matrix& rSensitivityMatrix,
                                             Vector& rSensitivityGradient,
                                             const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;
        rSensitivityGradient = ZeroVector(rSensitivityMatrix.size1());
        KRATOS_CATCH("");
    }

    void AdjointNodalDisplacementRootMeanSquareResponseFunction::CalculatePartialSensitivity(Condition& rAdjointCondition,
                                             const Variable<array_1d<double, 3>>& rVariable,
                                             const Matrix& rSensitivityMatrix,
                                             Vector& rSensitivityGradient,
                                             const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;
        rSensitivityGradient = ZeroVector(rSensitivityMatrix.size1());
        KRATOS_CATCH("");
    }

    double AdjointNodalDisplacementRootMeanSquareResponseFunction::CalculateValue(ModelPart& rModelPart)
    {
        KRATOS_TRY;

        const Variable<double>* p_traced_dof = &KratosComponents<Variable<double>>::Get(mTracedDofLabel);
        ModelPart& response_part = rModelPart.GetSubModelPart(mResponsePartName);

        double value = 0;
        double x_i;
        for(auto& node_i : response_part.Nodes()){
            x_i = 1/ mTimeDomain * rModelPart.GetNode(node_i.Id()).FastGetSolutionStepValue(*p_traced_dof , 0);
            value += x_i * x_i;
        }

        return value;

        KRATOS_CATCH("");
    }

    /// Find one element which is bounded by the traced node. The element is needed for assembling the adjoint load.
    void AdjointNodalDisplacementRootMeanSquareResponseFunction::ComputeNeighboringElementNodeMap()
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

    
} // namespace Kratos.


