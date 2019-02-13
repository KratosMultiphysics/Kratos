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
#include "adjoint_nodal_displacement_response_function.h"

namespace Kratos
{
    AdjointNodalDisplacementResponseFunction::AdjointNodalDisplacementResponseFunction(ModelPart& rModelPart, Parameters ResponseSettings)
    : AdjointStructuralResponseFunction(rModelPart, ResponseSettings)
    {
        // Get id of node where a displacement should be traced
        const int id_traced_node = ResponseSettings["traced_node_id"].GetInt();

        // Get the corresponding dof to the displacement which should be traced
        // by this response function e.g. DISPLACEMENT_X, ROTATION_X,...
        mTracedDofLabel = ResponseSettings["traced_dof"].GetString();

        // Get pointer to traced node
        mpTracedNode = rModelPart.pGetNode(id_traced_node);

        // Check if variable for traced dof is valid
        if( !( KratosComponents< VariableComponent< VectorComponentAdaptor<array_1d<double, 3> > > >::Has(mTracedDofLabel)) )
            KRATOS_ERROR << "Specified traced DOF is not available. Specified DOF: " << mTracedDofLabel << std::endl;
        else
        {
            const VariableComponentType& r_traced_dof =
                KratosComponents<VariableComponentType>::Get(mTracedDofLabel);
            KRATOS_ERROR_IF_NOT( mpTracedNode->SolutionStepsDataHas(r_traced_dof) )
                << "Specified DOF is not available at traced node." << std::endl;
        }

        // Check if variable for traced adjoint dof is valid
        if( !(KratosComponents< VariableComponent< VectorComponentAdaptor<array_1d<double, 3> > > >::Has(std::string("ADJOINT_") + mTracedDofLabel)) )
        {
            KRATOS_ERROR << "Specified traced adjoint DOF is not available." << std::endl;
        }

        this->GetNeighboringElementPointer();
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

        if( rAdjointElement.Id() == mpNeighboringElement->Id() )
        {
            DofsVectorType dofs_of_element;
            ProcessInfo process_info = rProcessInfo;
            mpNeighboringElement->GetDofList(dofs_of_element, process_info);

            const VariableComponentType& r_traced_adjoint_dof =
                KratosComponents<VariableComponentType>::Get(std::string("ADJOINT_") + mTracedDofLabel);

            for(IndexType i = 0; i < dofs_of_element.size(); ++i)
            {
                if (dofs_of_element[i]->Id() == mpTracedNode->Id() &&
                    dofs_of_element[i]->GetVariable() == r_traced_adjoint_dof)
                {
                    rResponseGradient[i] = -1;
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

        const VariableComponentType& r_traced_dof =
            KratosComponents<VariableComponentType>::Get(mTracedDofLabel);

        return mpTracedNode->FastGetSolutionStepValue(r_traced_dof, 0);

        KRATOS_CATCH("");
    }

    /// Find one element which is bounded by the traced node. The element is needed for assembling the adjoint load.
    void AdjointNodalDisplacementResponseFunction::GetNeighboringElementPointer()
    {
        KRATOS_TRY;

        for (auto elem_it = mrModelPart.Elements().ptr_begin(); elem_it != mrModelPart.Elements().ptr_end(); ++elem_it)
        {
            const SizeType number_of_nodes = (*elem_it)->GetGeometry().PointsNumber();
            for(IndexType i = 0; i < number_of_nodes; ++i)
            {
                if((*elem_it)->GetGeometry()[i].Id() == mpTracedNode->Id())
                {
                    mpNeighboringElement = (*elem_it);
                    return;
                }
            }
        }
        KRATOS_ERROR << "No neighboring element is available for the traced node." << std::endl;

        KRATOS_CATCH("");
    }

} // namespace Kratos.


