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
        ModelPart& r_model_part = this->GetModelPart();

        // This response function currently only works in 3D!
        ProcessInfo& r_current_process_info = r_model_part.GetProcessInfo();
        const int domain_size = r_current_process_info[DOMAIN_SIZE];
        KRATOS_ERROR_IF(domain_size != 3) << "Invalid DOMAIN_SIZE: " << domain_size << std::endl;

        // Get id of node where a displacement should be traced
        const int id_traced_node = ResponseSettings["traced_node_id"].GetInt();

        // Get the corresponding dof to the displacement which should be traced
        // by this response function e.g. DISPLACEMENT_X, ROTATION_X,...
        mTracedDofLabel = ResponseSettings["traced_dof"].GetString();

        // Get pointer to traced node
        mpTracedNode = r_model_part.pGetNode(id_traced_node);

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
        else
        {
            const VariableComponentType& r_traced_adjoint_dof =
                KratosComponents<VariableComponentType>::Get(std::string("ADJOINT_") + mTracedDofLabel);
            KRATOS_ERROR_IF_NOT( mpTracedNode->SolutionStepsDataHas(r_traced_adjoint_dof) )
                << "Specified adjoint DOF is not available at traced node." << std::endl;
        }

        this->GetNeighboringElementPointer();
    }

    AdjointNodalDisplacementResponseFunction::~AdjointNodalDisplacementResponseFunction(){}

    /// Find one element which is bounded by the traced node. The element is needed for assembling the adjoint load.
    void AdjointNodalDisplacementResponseFunction::GetNeighboringElementPointer()
    {
        KRATOS_TRY;

        ModelPart& r_model_part = this->GetModelPart();
        bool neighboring_element_found = false;

        for (auto& elem_i : r_model_part.Elements())
        {
            const SizeType number_of_nodes = elem_i.GetGeometry().size();
            for(IndexType i = 0; i < number_of_nodes; ++i)
            {
                if(elem_i.GetGeometry()[i].Id() == mpTracedNode->Id())
                {
                    mpNeighboringElement = r_model_part.pGetElement(elem_i.Id());
                    neighboring_element_found = true;
                    break;
                }
            }
            if(neighboring_element_found) { break; }
        }
        KRATOS_ERROR_IF_NOT(neighboring_element_found)
             << "No neighboring element is available for the traced node." << std::endl;

        KRATOS_CATCH("");

    }

    double AdjointNodalDisplacementResponseFunction::CalculateValue(ModelPart& rModelPart)
    {
        KRATOS_TRY;

        const VariableComponentType& r_traced_dof =
            KratosComponents<VariableComponentType>::Get(mTracedDofLabel);

        return rModelPart.Nodes()[(mpTracedNode->Id())].FastGetSolutionStepValue(r_traced_dof, 0);

        KRATOS_CATCH("");
    }

    void AdjointNodalDisplacementResponseFunction::CalculateGradient(const Element& rAdjointElem, const Matrix& rAdjointMatrix,
                                   Vector& rResponseGradient,
                                   ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;

        if (rResponseGradient.size() != rAdjointMatrix.size1())
            rResponseGradient.resize(rAdjointMatrix.size1(), false);

        rResponseGradient.clear();

        if( rAdjointElem.Id() == mpNeighboringElement->Id() )
        {
            DofsVectorType dofs_of_lement;
            mpNeighboringElement->GetDofList(dofs_of_lement,rProcessInfo);

            const VariableComponentType& r_traced_adjoint_dof =
                KratosComponents<VariableComponentType>::Get(std::string("ADJOINT_") + mTracedDofLabel);

            for(IndexType i = 0; i < dofs_of_lement.size(); ++i)
            {
                if(mpTracedNode->pGetDof(r_traced_adjoint_dof) == dofs_of_lement[i])
                {
                    rResponseGradient[i] = -1;
                }
            }
        }

        KRATOS_CATCH("");
    }

    void AdjointNodalDisplacementResponseFunction::CalculateSensitivityGradient(Element& rAdjointElem,
                                      const Variable<double>& rVariable,
                                      const Matrix& rDerivativesMatrix,
                                      Vector& rResponseGradient,
                                      ProcessInfo& rProcessInfo)
    {
          KRATOS_TRY

        if (rResponseGradient.size() != rDerivativesMatrix.size1())
            rResponseGradient.resize(rDerivativesMatrix.size1(), false);

        rResponseGradient.clear();


        KRATOS_CATCH("")
    }

    void AdjointNodalDisplacementResponseFunction::CalculateSensitivityGradient(Condition& rAdjointCondition,
                                     const Variable<double>& rVariable,
                                     const Matrix& rDerivativesMatrix,
                                     Vector& rResponseGradient,
                                     ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;

        if (rResponseGradient.size() != rDerivativesMatrix.size1())
            rResponseGradient.resize(rDerivativesMatrix.size1(), false);

        rResponseGradient.clear();

        KRATOS_CATCH("");
    }

    void AdjointNodalDisplacementResponseFunction::CalculateSensitivityGradient(Element& rAdjointElem,
                                      const Variable<array_1d<double,3>>& rVariable,
                                      const Matrix& rDerivativesMatrix,
                                      Vector& rResponseGradient,
                                      ProcessInfo& rProcessInfo)
    {
          KRATOS_TRY

        if (rResponseGradient.size() != rDerivativesMatrix.size1())
            rResponseGradient.resize(rDerivativesMatrix.size1(), false);

        rResponseGradient.clear();

          KRATOS_CATCH("")
    }

    void AdjointNodalDisplacementResponseFunction::CalculateSensitivityGradient(Condition& rAdjointCondition,
                                      const Variable<array_1d<double,3>>& rVariable,
                                      const Matrix& rDerivativesMatrix,
                                      Vector& rResponseGradient,
                                      ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;

        if (rResponseGradient.size() != rDerivativesMatrix.size1())
            rResponseGradient.resize(rDerivativesMatrix.size1(), false);

        rResponseGradient.clear();

        KRATOS_CATCH("");
    }

} // namespace Kratos.


