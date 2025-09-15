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

#pragma once


// System includes

// External includes

// Project includes
#include "adjoint_structural_response_function.h"


namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/** \brief AdjointNodalReactionResponseFunction
*
* This is a response function which traces a chosen reaction force of a single
* node as response. It is designed to be used in adjoint sensitivity analysis.
*/
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) AdjointNodalReactionResponseFunction : public AdjointStructuralResponseFunction
{
public:
    ///@name Type Definitions
    ///@{

    typedef Element::DofsVectorType DofsVectorType;
    typedef Node::Pointer PointTypePointer;
    typedef matrix_column< Matrix > MatrixColumnType;
    typedef matrix_row< Matrix > MatrixRowType;

    ///@}
    ///@name Pointer Definitions
    /// Pointer definition of AdjointNodalReactionResponseFunction
    KRATOS_CLASS_POINTER_DEFINITION(AdjointNodalReactionResponseFunction);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    AdjointNodalReactionResponseFunction(ModelPart& rModelPart, Parameters ResponseSettings);

    /// Destructor.
    ~AdjointNodalReactionResponseFunction() override;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    void InitializeSolutionStep() override;

    void FinalizeSolutionStep() override;

    using AdjointStructuralResponseFunction::CalculateGradient;

    void CalculateGradient(const Element& rAdjointElement,
                                   const Matrix& rResidualGradient,
                                   Vector& rResponseGradient,
                                   const ProcessInfo& rProcessInfo) override;

    void CalculateFirstDerivativesGradient(const Element& rAdjointElement,
                                           const Matrix& rResidualGradient,
                                           Vector& rResponseGradient,
                                           const ProcessInfo& rProcessInfo) override;

    void CalculateFirstDerivativesGradient(const Condition& rAdjointCondition,
                                           const Matrix& rResidualGradient,
                                           Vector& rResponseGradient,
                                           const ProcessInfo& rProcessInfo) override;

    void CalculateSecondDerivativesGradient(const Element& rAdjointElement,
                                            const Matrix& rResidualGradient,
                                            Vector& rResponseGradient,
                                            const ProcessInfo& rProcessInfo) override;

    void CalculateSecondDerivativesGradient(const Condition& rAdjointCondition,
                                            const Matrix& rResidualGradient,
                                            Vector& rResponseGradient,
                                            const ProcessInfo& rProcessInfo) override;

    void CalculatePartialSensitivity(Element& rAdjointElement,
                                             const Variable<double>& rVariable,
                                             const Matrix& rSensitivityMatrix,
                                             Vector& rSensitivityGradient,
                                             const ProcessInfo& rProcessInfo) override;

    void CalculatePartialSensitivity(Condition& rAdjointCondition,
                                             const Variable<double>& rVariable,
                                             const Matrix& rSensitivityMatrix,
                                             Vector& rSensitivityGradient,
                                             const ProcessInfo& rProcessInfo) override;

    void CalculatePartialSensitivity(Element& rAdjointElement,
                                             const Variable<array_1d<double, 3>>& rVariable,
                                             const Matrix& rSensitivityMatrix,
                                             Vector& rSensitivityGradient,
                                             const ProcessInfo& rProcessInfo) override;

    void CalculatePartialSensitivity(Condition& rAdjointCondition,
                                             const Variable<array_1d<double, 3>>& rVariable,
                                             const Matrix& rSensitivityMatrix,
                                             Vector& rSensitivityGradient,
                                             const ProcessInfo& rProcessInfo) override;


    double CalculateValue(ModelPart& rModelPart) override;

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    ///@}
    ///@name Friends
    ///@{

    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{

    ///@}

private:

    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    std::string mTracedDisplacementLabel;
    std::string mTracedReactionLabel;
    PointTypePointer  mpTracedNode;
    GlobalPointersVector<Element> mpNeighbourElements;
    GlobalPointersVector<Condition> mpNeighbourConditions;
    bool mAdjustAdjointDisplacement = false;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{


    void CalculateContributionToPartialSensitivity(Element& rAdjointElement,
                                             const Matrix& rSensitivityMatrix,
                                             Vector& rSensitivityGradient,
                                             const ProcessInfo& rProcessInfo);

    void CalculateContributionToPartialSensitivity(Condition& rAdjointCondition,
                                             const Matrix& rSensitivityMatrix,
                                             Vector& rSensitivityGradient,
                                             const ProcessInfo& rProcessInfo);


    template <typename TObjectType>
    size_t GetDofIndex(TObjectType& rAdjointObject, const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;

        const Variable<double>& r_corresponding_adjoint_dof =
        KratosComponents<Variable<double>>::Get(std::string("ADJOINT_") + mTracedDisplacementLabel);

        DofsVectorType dof_list;
        rAdjointObject.GetDofList(dof_list, rProcessInfo);

        size_t dof_index = 0;
        for(IndexType dof_it = 0; dof_it < dof_list.size(); ++dof_it)
        {
            if (dof_list[dof_it]->Id() == mpTracedNode->Id() &&
                dof_list[dof_it]->GetVariable() == r_corresponding_adjoint_dof)
            {
                dof_index = dof_it;
                break;
            }
        }
        return dof_index;

        KRATOS_CATCH("");
    }

    Vector GetColumnCopy(const Matrix& rMatrix, size_t ColumnIndex);

    Vector GetRowCopy(const Matrix& rMatrix, size_t RowIndex);

    std::string GetCorrespondingDisplacementLabel(std::string& rReactionLabel) const;

    void PerformResponseVariablesCheck();

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}

}; // Class AdjointNodalReactionResponseFunction

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

} // namespace Kratos.
