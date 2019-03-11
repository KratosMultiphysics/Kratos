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
#include "adjoint_nodal_reaction_response_function.h"
#include "utilities/variable_utils.h"
#include "processes/find_nodal_neighbours_process.h"
#include "processes/find_conditions_neighbours_process.h"

namespace Kratos
{
    AdjointNodalReactionResponseFunction::AdjointNodalReactionResponseFunction(ModelPart& rModelPart, Parameters ResponseSettings)
    : AdjointStructuralResponseFunction(rModelPart, ResponseSettings)
    {
        // Get id of node where a reaction should be traced
        const int id_traced_node = ResponseSettings["traced_node_id"].GetInt();

        // Get the corresponding dof to the reaction which should be traced
        // by this response function e.g. REACTION_X, REACTION_MOMENT_X,...
        mTracedReactionLabel = ResponseSettings["traced_reaction"].GetString();
        // Get the corresponding dof to the displacement e.g. REACTION_X --> DISPLACEMENT_X
        mTracedDisplacementDofLabel = this->GetCorrespondingDisplacementDofLabel(mTracedReactionLabel);

        if(ResponseSettings.Has("adjust_influence_function"))
            mAdjustInfluenceFunction = ResponseSettings["adjust_influence_function"].GetBool();

        // Get pointer to traced node
        mpTracedNode = rModelPart.pGetNode(id_traced_node);

        // Check the given variables after reading
        this->PerformResponseVariablesCheck();

        // Find neighbour elements and conditions because they are needed to construct the partial derivatives
        FindNodalNeighboursProcess neigbhor_elements_finder = FindNodalNeighboursProcess(mrModelPart, 10, 10);
        FindConditionsNeighboursProcess neigbhor_conditions_finder = FindConditionsNeighboursProcess(mrModelPart, 10, 10);
        neigbhor_elements_finder.Execute();
        neigbhor_conditions_finder.Execute();
        mpNeighborElements = mpTracedNode->GetValue(NEIGHBOUR_ELEMENTS);
        mpNeighborConditions = mpTracedNode->GetValue(NEIGHBOUR_CONDITIONS);
    }

    AdjointNodalReactionResponseFunction::~AdjointNodalReactionResponseFunction(){}

    void AdjointNodalReactionResponseFunction::InitializeSolutionStep()
    {
        KRATOS_TRY;

        const VariableComponentType& r_traced_adjoint_dof =
            KratosComponents<VariableComponentType>::Get(std::string("ADJOINT_") + mTracedDisplacementDofLabel);

        // check if the given node is really fixed in the direction of the given dof.
        if((mpTracedNode->GetDof(r_traced_adjoint_dof)).IsFree())
            KRATOS_ERROR << "Chosen reaction '" << mTracedReactionLabel << "' is not fixed on node #" << mpTracedNode->Id() << "!" << std::endl;

        KRATOS_CATCH("");
    }

    void AdjointNodalReactionResponseFunction::FinalizeSolutionStep()
    {
        KRATOS_TRY;

        if(mAdjustInfluenceFunction)
        {
            const VariableComponentType& r_traced_adjoint_dof =
                KratosComponents<VariableComponentType>::Get(std::string("ADJOINT_") + mTracedDisplacementDofLabel);
            mpTracedNode->FastGetSolutionStepValue(r_traced_adjoint_dof) = -1.0;

            mPerformedModificationOfAdjointDisplacements = true;
        }

        KRATOS_CATCH("");
    }

    void AdjointNodalReactionResponseFunction::CalculateGradient(const Element& rAdjointElement,
                                   const Matrix& rResidualGradient,
                                   Vector& rResponseGradient,
                                   const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;

        if (rResponseGradient.size() != rResidualGradient.size1())
            rResponseGradient.resize(rResidualGradient.size1(), false);
        rResponseGradient.clear();

        for(IndexType i = 0; i < mpNeighborElements.size(); ++i)
        {
            Kratos::Element& ng_elem_i = mpNeighborElements[i];

            if( rAdjointElement.Id() == ng_elem_i.Id() )
            {
                Vector filter_vector;
                ProcessInfo process_info = rProcessInfo;
                this->ConstructFilterVector(ng_elem_i, filter_vector, process_info);

                Matrix left_hand_side;
                ng_elem_i.CalculateLeftHandSide(left_hand_side, process_info);

                noalias(rResponseGradient) = prod(left_hand_side, filter_vector);
            }
        }

        KRATOS_CATCH("");
    }

    void AdjointNodalReactionResponseFunction::CalculateFirstDerivativesGradient(
        const Element& rAdjointElement,
        const Matrix& rResidualGradient,
        Vector& rResponseGradient,
        const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;
        rResponseGradient = ZeroVector(rResidualGradient.size1());
        KRATOS_CATCH("");
    }

    void AdjointNodalReactionResponseFunction::CalculateFirstDerivativesGradient(
        const Condition& rAdjointCondition,
        const Matrix& rResidualGradient,
        Vector& rResponseGradient,
        const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;
        rResponseGradient = ZeroVector(rResidualGradient.size1());
        KRATOS_CATCH("");
    }

    void AdjointNodalReactionResponseFunction::CalculateSecondDerivativesGradient(
        const Element& rAdjointElement,
        const Matrix& rResidualGradient,
        Vector& rResponseGradient,
        const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;
        rResponseGradient = ZeroVector(rResidualGradient.size1());
        KRATOS_CATCH("");
    }

    void AdjointNodalReactionResponseFunction::CalculateSecondDerivativesGradient(
        const Condition& rAdjointCondition,
        const Matrix& rResidualGradient,
        Vector& rResponseGradient,
        const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;
        rResponseGradient = ZeroVector(rResidualGradient.size1());
        KRATOS_CATCH("");
    }

    void AdjointNodalReactionResponseFunction::CalculatePartialSensitivity(Element& rAdjointElement,
                                             const Variable<double>& rVariable,
                                             const Matrix& rSensitivityMatrix,
                                             Vector& rSensitivityGradient,
                                             const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;

        rSensitivityGradient = ZeroVector(rSensitivityMatrix.size1());

        if(!mPerformedModificationOfAdjointDisplacements)
            this->CalculateContributionToPartialSensitivity(rAdjointElement,
                                rSensitivityMatrix, rSensitivityGradient, rProcessInfo);

        KRATOS_CATCH("");
    }

    void AdjointNodalReactionResponseFunction::CalculatePartialSensitivity(Condition& rAdjointCondition,
                                             const Variable<double>& rVariable,
                                             const Matrix& rSensitivityMatrix,
                                             Vector& rSensitivityGradient,
                                             const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;

        rSensitivityGradient = ZeroVector(rSensitivityMatrix.size1());

        if(!mPerformedModificationOfAdjointDisplacements)
            this->CalculateContributionToPartialSensitivity(rAdjointCondition,
                                rSensitivityMatrix, rSensitivityGradient, rProcessInfo);

        KRATOS_CATCH("");
    }

    void AdjointNodalReactionResponseFunction::CalculatePartialSensitivity(Element& rAdjointElement,
                                             const Variable<array_1d<double, 3>>& rVariable,
                                             const Matrix& rSensitivityMatrix,
                                             Vector& rSensitivityGradient,
                                             const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;

        rSensitivityGradient = ZeroVector(rSensitivityMatrix.size1());

        if(!mPerformedModificationOfAdjointDisplacements)
            this->CalculateContributionToPartialSensitivity(rAdjointElement,
                                rSensitivityMatrix, rSensitivityGradient, rProcessInfo);

        KRATOS_CATCH("");
    }

    void AdjointNodalReactionResponseFunction::CalculatePartialSensitivity(Condition& rAdjointCondition,
                                             const Variable<array_1d<double, 3>>& rVariable,
                                             const Matrix& rSensitivityMatrix,
                                             Vector& rSensitivityGradient,
                                             const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;

        rSensitivityGradient = ZeroVector(rSensitivityMatrix.size1());

        if(!mPerformedModificationOfAdjointDisplacements)
            this->CalculateContributionToPartialSensitivity(rAdjointCondition,
                                rSensitivityMatrix, rSensitivityGradient, rProcessInfo);

        KRATOS_CATCH("");
    }

    double AdjointNodalReactionResponseFunction::CalculateValue(ModelPart& rModelPart)
    {
        KRATOS_TRY;

        const VariableComponentType& r_traced_dof =
            KratosComponents<VariableComponentType>::Get(mTracedReactionLabel);

        return rModelPart.GetNode(mpTracedNode->Id()).FastGetSolutionStepValue(r_traced_dof, 0);

        KRATOS_CATCH("");
    }

    void AdjointNodalReactionResponseFunction::CalculateContributionToPartialSensitivity(Element& rAdjointElement,
                                             const Matrix& rSensitivityMatrix,
                                             Vector& rSensitivityGradient,
                                             const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;

        for(IndexType i = 0; i < mpNeighborElements.size(); ++i)
        {
            Kratos::Element& ng_elem_i = mpNeighborElements[i];

            if( rAdjointElement.Id() == ng_elem_i.Id() )
            {
                Vector filter_vector;
                ProcessInfo process_info = rProcessInfo;
                this->ConstructFilterVector(rAdjointElement, filter_vector, process_info);
                noalias(rSensitivityGradient) = prod(rSensitivityMatrix, filter_vector) ;
            }
        }

        KRATOS_CATCH("");
    }

    void AdjointNodalReactionResponseFunction::CalculateContributionToPartialSensitivity(Condition& rAdjointCondition,
                                             const Matrix& rSensitivityMatrix,
                                             Vector& rSensitivityGradient,
                                             const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;

        for(IndexType i = 0; i < mpNeighborConditions.size(); ++i)
        {
            Kratos::Condition& ng_cond_i = mpNeighborConditions[i];

            if( rAdjointCondition.Id() == ng_cond_i.Id() )
            {
                Vector filter_vector;
                ProcessInfo process_info = rProcessInfo;
                this->ConstructFilterVector(rAdjointCondition, filter_vector, process_info);
                noalias(rSensitivityGradient) = prod(rSensitivityMatrix, filter_vector);
            }
        }

        KRATOS_CATCH("");
    }


    std::string AdjointNodalReactionResponseFunction::GetCorrespondingDisplacementDofLabel(std::string& rReactionDofLabel) const
    {
        KRATOS_TRY;

        if (rReactionDofLabel == "REACTION_X")
            return "DISPLACEMENT_X";
        else if (rReactionDofLabel == "REACTION_Y")
            return "DISPLACEMENT_Y";
        else if (rReactionDofLabel == "REACTION_Z")
            return "DISPLACEMENT_Z";
        else if (rReactionDofLabel == "REACTION_MOMENT_X")
            return "ROTATION_X";
        else if (rReactionDofLabel == "REACTION_MOMENT_Y")
            return "ROTATION_Y";
        else if (rReactionDofLabel == "REACTION_MOMENT_Z")
            return "ROTATION_Z";
        else
            KRATOS_ERROR << "Invalid reaction dof label!" << std::endl;

        KRATOS_CATCH("");
    }

    void AdjointNodalReactionResponseFunction::PerformResponseVariablesCheck()
    {
        KRATOS_TRY;

        // Check if variable for traced reaction is valid
        if( !( KratosComponents< VariableComponent< VectorComponentAdaptor<array_1d<double, 3> > > >::Has(mTracedReactionLabel)) )
            KRATOS_ERROR << "Specified traced DOF is not available. Specified DOF: " << mTracedReactionLabel << std::endl;
        else
        {
            const VariableComponentType& r_traced_dof =
                KratosComponents<VariableComponentType>::Get(mTracedReactionLabel);
            KRATOS_ERROR_IF_NOT( mpTracedNode->SolutionStepsDataHas(r_traced_dof) )
                << "Specified DOF is not available at traced node." << std::endl;
        }

        // Check if variable for traced dof is valid
        if( !( KratosComponents< VariableComponent< VectorComponentAdaptor<array_1d<double, 3> > > >::Has(mTracedDisplacementDofLabel)) )
            KRATOS_ERROR << "Specified traced DOF is not available. Specified DOF: " << mTracedDisplacementDofLabel << std::endl;
        else
        {
            const VariableComponentType& r_traced_dof =
                KratosComponents<VariableComponentType>::Get(mTracedDisplacementDofLabel);
            KRATOS_ERROR_IF_NOT( mpTracedNode->SolutionStepsDataHas(r_traced_dof) )
                << "Specified DOF is not available at traced node." << std::endl;
        }

        // Check if variable for traced adjoint dof is valid
        if( !(KratosComponents< VariableComponent< VectorComponentAdaptor<array_1d<double, 3> > > >::Has(std::string("ADJOINT_") + mTracedDisplacementDofLabel)) )
        {
            KRATOS_ERROR << "Specified traced adjoint DOF is not available." << std::endl;
        }

        KRATOS_CATCH("");
    }

} // namespace Kratos.


