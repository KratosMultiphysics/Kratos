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

        // Get label of the reaction which should be traced
        // by this response function e.g. REACTION_X, REACTION_MOMENT_X,...
        mTracedReactionLabel = ResponseSettings["traced_reaction"].GetString();
        // Get the label of the correspondig displacement e.g. REACTION_X --> DISPLACEMENT_X
        mTracedDisplacementLabel = this->GetCorrespondingDisplacementLabel(mTracedReactionLabel);

        if(ResponseSettings.Has("adjust_adjoint_displacement"))
            mAdjustAdjointDisplacement = ResponseSettings["adjust_adjoint_displacement"].GetBool();

        // Get pointer to traced node
        mpTracedNode = rModelPart.pGetNode(id_traced_node);

        // Check the given variables after reading
        this->PerformResponseVariablesCheck();

        // Find neighbour elements and conditions because they are needed to construct the partial derivatives
        FindNodalNeighboursProcess neigbhour_elements_finder(mrModelPart);
        FindConditionsNeighboursProcess neigbhour_conditions_finder = FindConditionsNeighboursProcess(mrModelPart);
        neigbhour_elements_finder.Execute();
        neigbhour_conditions_finder.Execute();
        mpNeighbourElements = mpTracedNode->GetValue(NEIGHBOUR_ELEMENTS);
        mpNeighbourConditions = mpTracedNode->GetValue(NEIGHBOUR_CONDITIONS);
    }

    AdjointNodalReactionResponseFunction::~AdjointNodalReactionResponseFunction(){}

    void AdjointNodalReactionResponseFunction::InitializeSolutionStep()
    {
        KRATOS_TRY;

        const Variable<double>& r_corresponding_adjoint_dof =
            KratosComponents<Variable<double>>::Get(std::string("ADJOINT_") + mTracedDisplacementLabel);

        // check if the given node is really fixed in the direction of the traced reaction.
        if((mpTracedNode->GetDof(r_corresponding_adjoint_dof)).IsFree())
            KRATOS_ERROR << "The dof '" << r_corresponding_adjoint_dof << "' corresponding to the chosen reaction '"<< mTracedReactionLabel << "' is not fixed on node #" << mpTracedNode->Id() << "!" << std::endl;

        KRATOS_CATCH("");
    }

    void AdjointNodalReactionResponseFunction::FinalizeSolutionStep()
    {
        KRATOS_TRY;

        /*
        With this adjustment the adjoint displacements will be equivalent to the discrete influence function.
        The modification can be justified by the partial derivative of the response w.r.t. to the design variable.
        Please note: adjusting the adjoint displacement replaces the upcoming computation of the partial derivative of the response
        w.r.t. to the design variable.
        */
        if(mAdjustAdjointDisplacement)
        {
            const Variable<double>& r_corresponding_adjoint_dof =
                KratosComponents<Variable<double>>::Get(std::string("ADJOINT_") + mTracedDisplacementLabel);
            mpTracedNode->FastGetSolutionStepValue(r_corresponding_adjoint_dof) = -1.0;
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

        for(IndexType i = 0; i < mpNeighbourElements.size(); ++i)
        {
            Kratos::Element& ng_elem_i = mpNeighbourElements[i];

            if( rAdjointElement.Id() == ng_elem_i.Id() )
            {
                Matrix left_hand_side;
                ng_elem_i.CalculateLeftHandSide(left_hand_side, rProcessInfo);
                auto dof_index = this->GetDofIndex(ng_elem_i, rProcessInfo);
                rResponseGradient = -1.0 * (this->GetColumnCopy(left_hand_side, dof_index));
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

        if(!mAdjustAdjointDisplacement)
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

        if(!mAdjustAdjointDisplacement)
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

        if(!mAdjustAdjointDisplacement)
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

        if(!mAdjustAdjointDisplacement)
            this->CalculateContributionToPartialSensitivity(rAdjointCondition,
                                rSensitivityMatrix, rSensitivityGradient, rProcessInfo);

        KRATOS_CATCH("");
    }

    double AdjointNodalReactionResponseFunction::CalculateValue(ModelPart& rModelPart)
    {
        KRATOS_TRY;

        const Variable<double>& r_traced_reaction =
            KratosComponents<Variable<double>>::Get(mTracedReactionLabel);

        return rModelPart.GetNode(mpTracedNode->Id()).FastGetSolutionStepValue(r_traced_reaction, 0);

        KRATOS_CATCH("");
    }

    void AdjointNodalReactionResponseFunction::CalculateContributionToPartialSensitivity(Element& rAdjointElement,
                                             const Matrix& rSensitivityMatrix,
                                             Vector& rSensitivityGradient,
                                             const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;

        for(IndexType i = 0; i < mpNeighbourElements.size(); ++i)
        {
            Kratos::Element& ng_elem_i = mpNeighbourElements[i];

            if( rAdjointElement.Id() == ng_elem_i.Id() )
            {
                auto dof_index = this->GetDofIndex(ng_elem_i, rProcessInfo);
                rSensitivityGradient = -1.0 * (this->GetColumnCopy(rSensitivityMatrix, dof_index));
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

        for(IndexType i = 0; i < mpNeighbourConditions.size(); ++i)
        {
            Kratos::Condition& ng_cond_i = mpNeighbourConditions[i];

            if( rAdjointCondition.Id() == ng_cond_i.Id() )
            {
                auto dof_index = this->GetDofIndex(ng_cond_i, rProcessInfo);
                rSensitivityGradient = -1.0 * (this->GetColumnCopy(rSensitivityMatrix, dof_index));
            }
        }

        KRATOS_CATCH("");
    }

    Vector AdjointNodalReactionResponseFunction::GetColumnCopy(const Matrix& rMatrix, size_t ColumnIndex)
    {
        KRATOS_TRY;

        auto num_rows = rMatrix.size1();
        auto num_columns = rMatrix.size2();
        KRATOS_ERROR_IF(ColumnIndex > num_columns) << "The provided column index is not valid!" << std::endl;
        Vector column = ZeroVector(num_rows);
        for(IndexType i = 0; i < num_rows; ++i)
            column[i] = rMatrix(i, ColumnIndex);
        return column;

        KRATOS_CATCH("");
    }

    Vector AdjointNodalReactionResponseFunction::GetRowCopy(const Matrix& rMatrix, size_t RowIndex)
    {
        KRATOS_TRY;

        auto num_rows = rMatrix.size1();
        auto num_columns = rMatrix.size2();
        KRATOS_ERROR_IF(RowIndex > num_rows) << "The provided row index is not valid!" << std::endl;
        Vector row = ZeroVector(num_columns);
        for(IndexType i = 0; i < num_columns; ++i)
            row[i] = rMatrix(RowIndex, i);
        return row;

        KRATOS_CATCH("");
    }


    std::string AdjointNodalReactionResponseFunction::GetCorrespondingDisplacementLabel(std::string& rReactionLabel) const
    {
        KRATOS_TRY;

        std::map<std::string, std::string> reaction_displacement_table;
        reaction_displacement_table["REACTION_X"] = "DISPLACEMENT_X";
        reaction_displacement_table["REACTION_Y"] = "DISPLACEMENT_Y";
        reaction_displacement_table["REACTION_Z"] = "DISPLACEMENT_Z";
        reaction_displacement_table["REACTION_MOMENT_X"] = "ROTATION_X";
        reaction_displacement_table["REACTION_MOMENT_Y"] = "ROTATION_Y";
        reaction_displacement_table["REACTION_MOMENT_Z"] = "ROTATION_Z";

        auto it_table = reaction_displacement_table.find(rReactionLabel);
        KRATOS_ERROR_IF(it_table == reaction_displacement_table.end()) << "Given reaction label is not valid!" << std::endl;

        return it_table->second;

        KRATOS_CATCH("");
    }

    void AdjointNodalReactionResponseFunction::PerformResponseVariablesCheck()
    {
        KRATOS_TRY;

        // Check if variable for traced reaction is valid
        if( !( KratosComponents< Variable<double> >::Has(mTracedReactionLabel)) )
            KRATOS_ERROR << mTracedReactionLabel << " is not available!" << std::endl;
        else
        {
            const Variable<double>& r_traced_reaction =
                KratosComponents<Variable<double>>::Get(mTracedReactionLabel);
            KRATOS_ERROR_IF_NOT( mpTracedNode->SolutionStepsDataHas(r_traced_reaction) )
                << mTracedReactionLabel << " is not available at traced node." << std::endl;
        }

        // Check if the displacement variable corresponding to the chosen reaction is valid
        if( !( KratosComponents< Variable<double> >::Has(mTracedDisplacementLabel)) )
            KRATOS_ERROR << mTracedDisplacementLabel << " is not available!" << std::endl;
        else
        {
            const Variable<double>& r_corresponding_dof =
                KratosComponents<Variable<double>>::Get(mTracedDisplacementLabel);
            KRATOS_ERROR_IF_NOT( mpTracedNode->SolutionStepsDataHas(r_corresponding_dof) )
                << mTracedDisplacementLabel << " is not available at traced node." << std::endl;
        }

        // Check if the adjoint displacement variable corresponding to the chosen reaction is valid
        if(!(KratosComponents< Variable<double> >::Has(std::string("ADJOINT_") + mTracedDisplacementLabel)))
        {
            KRATOS_ERROR << (std::string("ADJOINT_") + mTracedDisplacementLabel) << " is not available!" << std::endl;
        }

        KRATOS_CATCH("");
    }

} // namespace Kratos.


