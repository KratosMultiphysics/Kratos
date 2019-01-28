// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Baumgaertner Daniel, https://github.com/dbaumgaertner
//

// System includes

// External includes

// Project includes
#include "adjoint_external_response_function.h"

namespace Kratos
{
    AdjointExternalResponseFunction::AdjointExternalResponseFunction(ModelPart& rModelPart, Parameters ResponseSettings)
    : AdjointStructuralResponseFunction(rModelPart, ResponseSettings)
    {
        for(auto& elem_i : mrModelPart.Elements())
        {
            DofsVectorType dofs_of_element;
            elem_i.GetDofList(dofs_of_element, rModelPart.GetProcessInfo());
            unsigned int ndof_per_node = dofs_of_element.size() / elem_i.GetGeometry().size();
            if(ndof_per_node != mRequiredNumDofPerNode)
                KRATOS_ERROR << "AdjointExternalResponseFunction::CalculateGradient --> function not implemented for other than 3 dofs per node!" << std::endl;
        }
    }

    AdjointExternalResponseFunction::~AdjointExternalResponseFunction(){}

    void AdjointExternalResponseFunction::InitializeSolutionStep()
    {
        KRATOS_TRY;

        AdjointStructuralResponseFunction::InitializeSolutionStep();

        for(auto& node_i : mrModelPart.Nodes())
        {
            mIsNodeAlreadyConsideredInCalculateGradient[node_i.Id()] = false;
            mIsNodeAlreadyConsideredInCalculatePartialSensitivity[node_i.Id()] = false;
        }

        KRATOS_CATCH("");
    }

    void AdjointExternalResponseFunction::CalculateGradient(const Element& rAdjointElement,
                                   const Matrix& rResidualGradient,
                                   Vector& rResponseGradient,
                                   const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;

        if (rResponseGradient.size() != rResidualGradient.size1())
            rResponseGradient.resize(rResidualGradient.size1(), false);

        rResponseGradient.clear();

        unsigned int node_counter = 0;

        for(auto& node_i : rAdjointElement.GetGeometry())
        {
            if(mIsNodeAlreadyConsideredInCalculateGradient[node_i.Id()] == false)
            {
                array_1d<double,3>& local_derivative = node_i.FastGetSolutionStepValue(DFDU);

                // the Minus sign has to be included since the ResidualBasedLinearStrategy is to solve the adjoint system in the StructuralMechanicsApp
                rResponseGradient[node_counter*mRequiredNumDofPerNode + 0] = -local_derivative[0];
                rResponseGradient[node_counter*mRequiredNumDofPerNode + 1] = -local_derivative[1];
                rResponseGradient[node_counter*mRequiredNumDofPerNode + 2] = -local_derivative[2];

                mIsNodeAlreadyConsideredInCalculateGradient[node_i.Id()] = true;
            }
            node_counter++;
        }

        KRATOS_CATCH("");
    }

    void AdjointExternalResponseFunction::CalculatePartialSensitivity(Element& rAdjointElement,
                                             const Variable<double>& rVariable,
                                             const Matrix& rSensitivityMatrix,
                                             Vector& rSensitivityGradient,
                                             const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY

        if (rSensitivityGradient.size() != 0)
            rSensitivityGradient.resize(0, false);

        KRATOS_CATCH("")
    }

    void AdjointExternalResponseFunction::CalculatePartialSensitivity(Condition& rAdjointCondition,
                                             const Variable<double>& rVariable,
                                             const Matrix& rSensitivityMatrix,
                                             Vector& rSensitivityGradient,
                                             const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;

        if (rSensitivityGradient.size() != 0)
            rSensitivityGradient.resize(0, false);

        KRATOS_CATCH("");
    }

    void AdjointExternalResponseFunction::CalculatePartialSensitivity(Element& rAdjointElement,
                                             const Variable<array_1d<double, 3>>& rVariable,
                                             const Matrix& rSensitivityMatrix,
                                             Vector& rSensitivityGradient,
                                             const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY

        if(rVariable != SHAPE)
            KRATOS_ERROR << "AdjointExternalResponseFunction::CalculatePartialSensitivity --> Not implemented for other design variables than SHAPE!" << std::endl;

        if (rSensitivityGradient.size() != rSensitivityMatrix.size2())
            rSensitivityGradient.resize(rSensitivityMatrix.size2(), false);

        rSensitivityGradient.clear();

        unsigned int node_counter = 0;

        for(auto& node_i : rAdjointElement.GetGeometry())
        {
            if(mIsNodeAlreadyConsideredInCalculatePartialSensitivity[node_i.Id()] == false)
            {
                array_1d<double,3>& local_derivative = node_i.FastGetSolutionStepValue(DFDX);

                rSensitivityGradient[node_counter*3 + 0] = local_derivative[0];
                rSensitivityGradient[node_counter*3 + 1] = local_derivative[1];
                rSensitivityGradient[node_counter*3 + 2] = local_derivative[2];

                mIsNodeAlreadyConsideredInCalculatePartialSensitivity[node_i.Id()] = true;
            }
            node_counter++;
        }

        KRATOS_CATCH("")
    }

    void AdjointExternalResponseFunction::CalculatePartialSensitivity(Condition& rAdjointCondition,
                                             const Variable<array_1d<double, 3>>& rVariable,
                                             const Matrix& rSensitivityMatrix,
                                             Vector& rSensitivityGradient,
                                             const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;

        if (rSensitivityGradient.size() != 0)
            rSensitivityGradient.resize(0, false);

        KRATOS_CATCH("");
    }

    double AdjointExternalResponseFunction::CalculateValue(ModelPart& rModelPart)
    {
        KRATOS_TRY;

        return rModelPart.GetProcessInfo()[RESPONSE_VALUE];

        KRATOS_CATCH("");
    }

} // namespace Kratos.


