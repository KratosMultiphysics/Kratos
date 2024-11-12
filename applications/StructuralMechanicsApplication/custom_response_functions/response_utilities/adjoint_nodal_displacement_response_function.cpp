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

// External includes

// Project includes
#include "adjoint_nodal_displacement_response_function.h"
#include "processes/generic_find_elements_neighbours_process.h"

namespace Kratos
{
    AdjointNodalDisplacementResponseFunction::AdjointNodalDisplacementResponseFunction(ModelPart& rModelPart, Parameters ResponseSettings)
    : AdjointStructuralResponseFunction(rModelPart, ResponseSettings)
    {
        mResponsePartName = ResponseSettings["response_part_name"].GetString();
        mResponseDirection = ResponseSettings["direction"].GetVector();
        mTracedDofLabel = ResponseSettings["traced_dof"].GetString();

        if ( norm_2( mResponseDirection ) > 1.0e-7 ) {
            mResponseDirection /= norm_2( mResponseDirection );
        } else {
            KRATOS_ERROR << "AdjointNodalDisplacementResponseFunction: 'response_direction' must not have a norm of 0.0." << std::endl;
        }

        // Check if variable for traced dof is valid
        KRATOS_ERROR_IF_NOT( KratosComponents<ArrayVariableType>::Has(mTracedDofLabel) )
            << "AdjointNodalDisplacementResponseFunction: Specified traced DOF is not available. Specified DOF: " << mTracedDofLabel << std::endl;

        // Check if variable for traced adjoint dof is valid
        KRATOS_ERROR_IF_NOT( KratosComponents<ArrayVariableType>::Has(std::string("ADJOINT_") + mTracedDofLabel) )
            << "AdjointNodalDisplacementResponseFunction: Specified traced adjoint DOF is not available." << mTracedDofLabel << std::endl;

        ModelPart& response_part = rModelPart.GetSubModelPart(mResponsePartName);
        const ArrayVariableType& r_traced_dof = KratosComponents<ArrayVariableType>::Get(mTracedDofLabel);
        for(auto& node_i : response_part.Nodes()){
            KRATOS_ERROR_IF_NOT( node_i.SolutionStepsDataHas(r_traced_dof) )
                << "AdjointNodalDisplacementResponseFunction: Specified DOF is not available at traced node." << std::endl;
        }

        this->ComputeNeighboringElementNodeMap();
    }

    AdjointNodalDisplacementResponseFunction::~AdjointNodalDisplacementResponseFunction(){}

    void AdjointNodalDisplacementResponseFunction::CalculateGradient(const Element& rAdjointElement,
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
                        rResponseGradient[i]   = 1 * mResponseDirection[0];
                        rResponseGradient[i+1] = 1 * mResponseDirection[1];
                        rResponseGradient[i+2] = 1 * mResponseDirection[2];
                        break;
                    }
                }
            }
        }

        KRATOS_CATCH("");
    }

    void AdjointNodalDisplacementResponseFunction::CalculateFirstDerivativesGradient(
        const Element& rAdjointElement,
        const Matrix& rResidualGradient,
        Vector& rResponseGradient,
        const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;
        rResponseGradient = ZeroVector(rResidualGradient.size1());
        KRATOS_CATCH("");
    }

    void AdjointNodalDisplacementResponseFunction::CalculateFirstDerivativesGradient(
        const Condition& rAdjointCondition,
        const Matrix& rResidualGradient,
        Vector& rResponseGradient,
        const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;
        rResponseGradient = ZeroVector(rResidualGradient.size1());
        KRATOS_CATCH("");
    }

    void AdjointNodalDisplacementResponseFunction::CalculateSecondDerivativesGradient(
        const Element& rAdjointElement,
        const Matrix& rResidualGradient,
        Vector& rResponseGradient,
        const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;
        rResponseGradient = ZeroVector(rResidualGradient.size1());
        KRATOS_CATCH("");
    }

    void AdjointNodalDisplacementResponseFunction::CalculateSecondDerivativesGradient(
        const Condition& rAdjointCondition,
        const Matrix& rResidualGradient,
        Vector& rResponseGradient,
        const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;
        rResponseGradient = ZeroVector(rResidualGradient.size1());
        KRATOS_CATCH("");
    }

    void AdjointNodalDisplacementResponseFunction::CalculatePartialSensitivity(Element& rAdjointElement,
                                             const Variable<double>& rVariable,
                                             const Matrix& rSensitivityMatrix,
                                             Vector& rSensitivityGradient,
                                             const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;
        rSensitivityGradient = ZeroVector(rSensitivityMatrix.size1());
        KRATOS_CATCH("");
    }

    void AdjointNodalDisplacementResponseFunction::CalculatePartialSensitivity(Condition& rAdjointCondition,
                                             const Variable<double>& rVariable,
                                             const Matrix& rSensitivityMatrix,
                                             Vector& rSensitivityGradient,
                                             const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;
        rSensitivityGradient = ZeroVector(rSensitivityMatrix.size1());
        KRATOS_CATCH("");
    }

    void AdjointNodalDisplacementResponseFunction::CalculatePartialSensitivity(Element& rAdjointElement,
                                             const Variable<array_1d<double, 3>>& rVariable,
                                             const Matrix& rSensitivityMatrix,
                                             Vector& rSensitivityGradient,
                                             const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;
        rSensitivityGradient = ZeroVector(rSensitivityMatrix.size1());
        KRATOS_CATCH("");
    }

    void AdjointNodalDisplacementResponseFunction::CalculatePartialSensitivity(Condition& rAdjointCondition,
                                             const Variable<array_1d<double, 3>>& rVariable,
                                             const Matrix& rSensitivityMatrix,
                                             Vector& rSensitivityGradient,
                                             const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;
        rSensitivityGradient = ZeroVector(rSensitivityMatrix.size1());
        KRATOS_CATCH("");
    }

    double AdjointNodalDisplacementResponseFunction::CalculateValue(ModelPart& rModelPart)
    {
        KRATOS_TRY;

        const ArrayVariableType& r_traced_dof =
            KratosComponents<ArrayVariableType>::Get(mTracedDofLabel);

        double response_value = 0.0;
        ModelPart& response_part = rModelPart.GetSubModelPart(mResponsePartName);
        for(auto& node_i : response_part.Nodes()){
            // project displacement vector in the traced direction. As mResponseDirection is a normalized vector
            // the result of the inner product is already the displacement value in the traced direction.
            response_value += inner_prod(mResponseDirection, node_i.FastGetSolutionStepValue(r_traced_dof, 0));
        }

        return response_value;

        KRATOS_CATCH("");
    }

    /// Find one element which is bounded by one of the traced nodes. The elements are needed for assembling the adjoint load.
    void AdjointNodalDisplacementResponseFunction::ComputeNeighboringElementNodeMap()
    {
        KRATOS_TRY;

        ModelPart& response_part = mrModelPart.GetSubModelPart(mResponsePartName);
        GenericFindElementalNeighboursProcess neighbour_elements_finder(mrModelPart);
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


