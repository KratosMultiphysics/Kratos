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

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------

// ------------------------------------------------------------------------------
// External includes
// ------------------------------------------------------------------------------

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "adjoint_nodal_displacement_response_function.h"


namespace Kratos
{
    AdjointNodalDisplacementResponseFunction::AdjointNodalDisplacementResponseFunction(ModelPart& rModelPart, Parameters ResponseSettings)
    : AdjointStructuralResponseFunction(rModelPart, ResponseSettings)
    {
        ModelPart& r_model_part = this->GetModelPart();
        
        // This response function currently only works in 3D!
        ProcessInfo& r_current_process_info = r_model_part.GetProcessInfo();
        const unsigned int domain_size =
            static_cast<unsigned int>(r_current_process_info[DOMAIN_SIZE]);
        KRATOS_ERROR_IF(domain_size != 3) << "Invalid DOMAIN_SIZE: " << domain_size << std::endl;

        // Get id of node where a displacement should be traced
        mIdOfTracedNode = ResponseSettings["traced_node_id"].GetInt();

        // Get the corresponding dof to the displacement which should be traced
        // by this response function e.g. DISPLACEMENT_X, ROTATION_X,...
        mTracedDofLabel = ResponseSettings["traced_dof"].GetString();

        // Get pointer to traced node
        mpTracedNode = r_model_part.pGetNode(mIdOfTracedNode);

        // Check if variable for traced dof is valid
        if( !( KratosComponents< VariableComponent< VectorComponentAdaptor<array_1d<double, 3> > > >::Has(mTracedDofLabel)) )
            KRATOS_ERROR << "Specified traced DOF is not available. Specified DOF: " << mTracedDofLabel << std::endl;
        else
        {
            const VariableComponentType& rTRACED_DOF =
                KratosComponents<VariableComponentType>::Get(mTracedDofLabel);
            KRATOS_ERROR_IF_NOT( mpTracedNode->SolutionStepsDataHas(rTRACED_DOF) )
                << "Specified DOF is not available at traced node." << std::endl;
        }

        // Check if variable for traced adjoint dof is valid
        if( !(KratosComponents< VariableComponent< VectorComponentAdaptor<array_1d<double, 3> > > >::Has(std::string("ADJOINT_") + mTracedDofLabel)) )
        {
            KRATOS_ERROR << "Specified traced adjoint DOF is not available." << std::endl;
        }
        else
        {
            const VariableComponentType& rTRACED_ADJOINT_DOF =
                KratosComponents<VariableComponentType>::Get(std::string("ADJOINT_") + mTracedDofLabel);
            KRATOS_ERROR_IF_NOT( mpTracedNode->SolutionStepsDataHas(rTRACED_ADJOINT_DOF) )
                << "Specified adjoint DOF is not available at traced node." << std::endl;
        }

        mDisplacementValue = 0.0;
        this->GetNeighboringElementPointer();
    }

    // ==============================================================================
    AdjointNodalDisplacementResponseFunction::~AdjointNodalDisplacementResponseFunction(){}

    // ==============================================================================
    void AdjointNodalDisplacementResponseFunction::Initialize() 
    {
        KRATOS_TRY;

        BaseType::Initialize();

        KRATOS_CATCH("");
    }

    // ==============================================================================
    void AdjointNodalDisplacementResponseFunction::GetNeighboringElementPointer()
    {
        KRATOS_TRY;

        ModelPart& r_model_part = this->GetModelPart();
        bool neighboring_element_found = false;

        for (ModelPart::ElementIterator it = r_model_part.ElementsBegin(); it != r_model_part.ElementsEnd(); ++it)
        {
            const unsigned int number_of_nodes = it->GetGeometry().size();
            for(unsigned int i = 0; i < number_of_nodes; ++i)
            {
                int current_node_id = it->GetGeometry()[i].Id();
                if(current_node_id == mIdOfTracedNode)
                {
                    mpNeighboringElement = r_model_part.pGetElement(it->Id());
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

    // ==============================================================================
    double AdjointNodalDisplacementResponseFunction::CalculateValue(ModelPart& rModelPart)
    {
        KRATOS_TRY;

        const VariableComponentType& rTRACED_DOF =
            KratosComponents<VariableComponentType>::Get(mTracedDofLabel);

        mDisplacementValue = rModelPart.Nodes()[(mpTracedNode->Id())].FastGetSolutionStepValue(rTRACED_DOF, 0);

        return mDisplacementValue;

        KRATOS_CATCH("");
    }

    // ==============================================================================
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

            const VariableComponentType& rTRACED_ADJOINT_DOF =
                KratosComponents<VariableComponentType>::Get(std::string("ADJOINT_") + mTracedDofLabel);

            for(unsigned int i = 0; i < dofs_of_lement.size(); ++i)
            {
                if(mpTracedNode->pGetDof(rTRACED_ADJOINT_DOF) == dofs_of_lement[i])
                {
                    rResponseGradient[i] = -1;
                }
            }
        }

        KRATOS_CATCH("");
    }

    // ==============================================================================
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

    // ==============================================================================
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

    // ==============================================================================
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

    // ==============================================================================
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


