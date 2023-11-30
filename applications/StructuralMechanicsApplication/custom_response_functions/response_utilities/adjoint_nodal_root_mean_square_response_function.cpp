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
#include "adjoint_nodal_root_mean_square_response_function.h"
#include "processes/find_elements_neighbours_process.h"

namespace Kratos
{
    AdjointNodalRootMeanSquareResponseFunction::AdjointNodalRootMeanSquareResponseFunction(ModelPart& rModelPart, Parameters ResponseSettings)
    : AdjointStructuralResponseFunction(rModelPart, ResponseSettings)
    { 
        mTimeDomain = ResponseSettings["time_domain"].GetDouble();
        mResponsePartName = ResponseSettings["response_part_name"].GetString();
        mTracedDofLabel = ResponseSettings["traced_dof"].GetString();
        mTracedDofDirection = ResponseSettings["direction"].GetString();
        mTracedDofTimeDerivativeOrder = ResponseSettings["time_derivative_order"].GetInt();

        // Check if variable for traced dof is valid
        KRATOS_ERROR_IF_NOT( KratosComponents<ArrayVariableType>::Has(mTracedDofLabel) && (mTracedDofLabel == "DISPLACEMENT" || mTracedDofLabel == "ROTATION") )
            << "AdjointNodalRootMeanSquareResponseFunction: Specified traced DOF is not available. Specified DOF: " << mTracedDofLabel << std::endl;

        // Check if time derivative order for traced dof is valid
        KRATOS_ERROR_IF_NOT((mTracedDofTimeDerivativeOrder == 0 || mTracedDofTimeDerivativeOrder == 1 || mTracedDofTimeDerivativeOrder == 2))
            << "AdjointNodalRootMeanSquareResponseFunction: Specified time derivative order is not available. Use 0, 1 or 2 as derivative order." << std::endl;

        if (mTracedDofLabel == "DISPLACEMENT"){
            switch (mTracedDofTimeDerivativeOrder) {
                case 0:
                    mTracedDynDofLabel = "DISPLACEMENT";
                    break;
                case 1:
                    mTracedDynDofLabel = "VELOCITY";
                    break;
                case 2:
                    mTracedDynDofLabel = "ACCELERATION";
            }
        }
        else if (mTracedDofLabel == "ROTATION"){
            switch (mTracedDofTimeDerivativeOrder) {
                case 0:
                    mTracedDynDofLabel = "ROTATION";
                    break;
                case 1:
                    mTracedDynDofLabel = "ANGULAR_VELOCITY";
                    break;
                case 2:
                    mTracedDynDofLabel = "ANGULAR_ACCELERATION";
            }
        }

        KRATOS_WATCH(mTracedDofLabel)
        KRATOS_WATCH(mTracedDofDirection)
        KRATOS_WATCH(mTracedDofTimeDerivativeOrder)
        KRATOS_WATCH(mTracedDynDofLabel)

        // Check if direction for traced dof is valid
        KRATOS_ERROR_IF_NOT((mTracedDofDirection == "X" || mTracedDofDirection == "Y" || mTracedDofDirection == "Z"))
            << "AdjointNodalDisplacementResponseFunction: Specified direction is not available. Use X, Y or Z as direction." << std::endl;

        ModelPart& response_part = rModelPart.GetSubModelPart(mResponsePartName);
        const ArrayVariableType& r_traced_dof = KratosComponents<ArrayVariableType>::Get(mTracedDofLabel);
        for(auto& node_i : response_part.Nodes()){
            KRATOS_ERROR_IF_NOT( node_i.SolutionStepsDataHas(r_traced_dof) )
                << "AdjointNodalDisplacementResponseFunction: Specified DOF is not available at node" << node_i.Id() << "." << std::endl;
        }

        this->ComputeNeighboringElementNodeMap();
    }

    AdjointNodalRootMeanSquareResponseFunction::~AdjointNodalRootMeanSquareResponseFunction(){}

    void AdjointNodalRootMeanSquareResponseFunction::CalculateGradient(const Element& rAdjointElement,
                                   const Matrix& rResidualGradient,
                                   Vector& rResponseGradient,
                                   const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;

        if (mTracedDofTimeDerivativeOrder == 0) {
            if (rResponseGradient.size() != rResidualGradient.size1())
                rResponseGradient.resize(rResidualGradient.size1(), false);

            rResponseGradient.clear();
            const Variable<double>& r_traced_dyn_dof = KratosComponents<Variable<double>>::Get(mTracedDynDofLabel + "_" + mTracedDofDirection);
            const Variable<double>& r_traced_adjoint_dof = KratosComponents<Variable<double>>::Get("ADJOINT_" + mTracedDofLabel + "_" + mTracedDofDirection);
            DofsVectorType dofs_of_element;
            auto it_map = mElementNodeMap.find(rAdjointElement.Id());
            if (it_map != mElementNodeMap.end()) {
                rAdjointElement.GetDofList(dofs_of_element, rProcessInfo);
                for(auto const& node_id: it_map->second) {
                    for(IndexType i = 0; i < dofs_of_element.size(); ++i) {
                        if (dofs_of_element[i]->Id() == node_id &&
                            dofs_of_element[i]->GetVariable() == r_traced_adjoint_dof) {
                            rResponseGradient[i]   = 2 / mTimeDomain * mrModelPart.GetNode(node_id).FastGetSolutionStepValue(r_traced_dyn_dof, 0);
                            break;
                        }
                    }
                }
            }
        }
        else {
            rResponseGradient = ZeroVector(rResidualGradient.size1());
        }
        
        KRATOS_CATCH("");
    }

    void AdjointNodalRootMeanSquareResponseFunction::CalculateFirstDerivativesGradient(
        const Element& rAdjointElement,
        const Matrix& rResidualGradient,
        Vector& rResponseGradient,
        const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;

        if (mTracedDofTimeDerivativeOrder == 1) {
            if (rResponseGradient.size() != rResidualGradient.size1())
                rResponseGradient.resize(rResidualGradient.size1(), false);

            rResponseGradient.clear();
            const Variable<double>& r_traced_dyn_dof = KratosComponents<Variable<double>>::Get(mTracedDynDofLabel + "_" + mTracedDofDirection);
            const Variable<double>& r_traced_adjoint_dof = KratosComponents<Variable<double>>::Get("ADJOINT_" + mTracedDofLabel + "_" + mTracedDofDirection);
            DofsVectorType dofs_of_element;
            auto it_map = mElementNodeMap.find(rAdjointElement.Id());
            if (it_map != mElementNodeMap.end()) {
                rAdjointElement.GetDofList(dofs_of_element, rProcessInfo);
                for(auto const& node_id: it_map->second) {
                    for(IndexType i = 0; i < dofs_of_element.size(); ++i) {
                        if (dofs_of_element[i]->Id() == node_id &&
                            dofs_of_element[i]->GetVariable() == r_traced_adjoint_dof) {
                            rResponseGradient[i]   = 2 / mTimeDomain * mrModelPart.GetNode(node_id).FastGetSolutionStepValue(r_traced_dyn_dof, 0);
                            KRATOS_WATCH(rAdjointElement.Id())
                            KRATOS_WATCH(mTimeDomain)
                            KRATOS_WATCH(mTracedDynDofLabel + "_" + mTracedDofDirection)
                            KRATOS_WATCH(mrModelPart.GetNode(node_id).FastGetSolutionStepValue(r_traced_dyn_dof, 0))
                            break;
                        }
                    }
                }
            }
            KRATOS_WATCH(rResponseGradient)
        }
        else {
            rResponseGradient = ZeroVector(rResidualGradient.size1());
        }
        
        KRATOS_CATCH("");
    }

    void AdjointNodalRootMeanSquareResponseFunction::CalculateFirstDerivativesGradient(
        const Condition& rAdjointCondition,
        const Matrix& rResidualGradient,
        Vector& rResponseGradient,
        const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;
        rResponseGradient = ZeroVector(rResidualGradient.size1());
        KRATOS_CATCH("");
    }

    void AdjointNodalRootMeanSquareResponseFunction::CalculateSecondDerivativesGradient(
        const Element& rAdjointElement,
        const Matrix& rResidualGradient,
        Vector& rResponseGradient,
        const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;

        if (mTracedDofTimeDerivativeOrder == 2) {
            if (rResponseGradient.size() != rResidualGradient.size1())
                rResponseGradient.resize(rResidualGradient.size1(), false);

            rResponseGradient.clear();
            const Variable<double>& r_traced_dyn_dof = KratosComponents<Variable<double>>::Get(mTracedDynDofLabel + "_" + mTracedDofDirection);
            const Variable<double>& r_traced_adjoint_dof = KratosComponents<Variable<double>>::Get("ADJOINT_" + mTracedDofLabel + "_" + mTracedDofDirection);
            DofsVectorType dofs_of_element;
            auto it_map = mElementNodeMap.find(rAdjointElement.Id());
            if (it_map != mElementNodeMap.end()) {
                rAdjointElement.GetDofList(dofs_of_element, rProcessInfo);
                for(auto const& node_id: it_map->second) {
                    for(IndexType i = 0; i < dofs_of_element.size(); ++i) {
                        if (dofs_of_element[i]->Id() == node_id &&
                            dofs_of_element[i]->GetVariable() == r_traced_adjoint_dof) {
                            rResponseGradient[i]   = 2 / mTimeDomain * mrModelPart.GetNode(node_id).FastGetSolutionStepValue(r_traced_dyn_dof, 0);
                            break;
                        }
                    }
                }
            }
        }
        else {
            rResponseGradient = ZeroVector(rResidualGradient.size1());
        }
        
        KRATOS_CATCH("");
    }

    void AdjointNodalRootMeanSquareResponseFunction::CalculateSecondDerivativesGradient(
        const Condition& rAdjointCondition,
        const Matrix& rResidualGradient,
        Vector& rResponseGradient,
        const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;
        rResponseGradient = ZeroVector(rResidualGradient.size1());
        KRATOS_CATCH("");
    }

    void AdjointNodalRootMeanSquareResponseFunction::CalculatePartialSensitivity(Element& rAdjointElement,
                                             const Variable<double>& rVariable,
                                             const Matrix& rSensitivityMatrix,
                                             Vector& rSensitivityGradient,
                                             const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;
        rSensitivityGradient = ZeroVector(rSensitivityMatrix.size1());
        KRATOS_CATCH("");
    }

    void AdjointNodalRootMeanSquareResponseFunction::CalculatePartialSensitivity(Condition& rAdjointCondition,
                                             const Variable<double>& rVariable,
                                             const Matrix& rSensitivityMatrix,
                                             Vector& rSensitivityGradient,
                                             const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;
        rSensitivityGradient = ZeroVector(rSensitivityMatrix.size1());
        KRATOS_CATCH("");
    }

    void AdjointNodalRootMeanSquareResponseFunction::CalculatePartialSensitivity(Element& rAdjointElement,
                                             const Variable<array_1d<double, 3>>& rVariable,
                                             const Matrix& rSensitivityMatrix,
                                             Vector& rSensitivityGradient,
                                             const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;
        rSensitivityGradient = ZeroVector(rSensitivityMatrix.size1());
        KRATOS_CATCH("");
    }

    void AdjointNodalRootMeanSquareResponseFunction::CalculatePartialSensitivity(Condition& rAdjointCondition,
                                             const Variable<array_1d<double, 3>>& rVariable,
                                             const Matrix& rSensitivityMatrix,
                                             Vector& rSensitivityGradient,
                                             const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;
        rSensitivityGradient = ZeroVector(rSensitivityMatrix.size1());
        KRATOS_CATCH("");
    }

    double AdjointNodalRootMeanSquareResponseFunction::CalculateValue(ModelPart& rModelPart)
    {
        KRATOS_TRY;

        const Variable<double>& r_traced_dyn_dof = KratosComponents<Variable<double>>::Get(mTracedDynDofLabel + "_" + mTracedDofDirection);
        ModelPart& response_part = rModelPart.GetSubModelPart(mResponsePartName);

        double value = 0;
        double x_i;
        for(auto& node_i : response_part.Nodes()){
            x_i = rModelPart.GetNode(node_i.Id()).FastGetSolutionStepValue(r_traced_dyn_dof , 0);
            value += x_i * x_i * 1/ mTimeDomain;
            KRATOS_WATCH(node_i.Id())
            KRATOS_WATCH(mTracedDynDofLabel + "_" + mTracedDofDirection)
            KRATOS_WATCH(x_i)
        }

        return value;

        KRATOS_CATCH("");
    }

    /// Find one element which is bounded by the traced node. The element is needed for assembling the adjoint load.
    void AdjointNodalRootMeanSquareResponseFunction::ComputeNeighboringElementNodeMap()
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


