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

#ifndef ADJOINT_NODAL_DISPLACEMENT_RESPONSE_FUNCTION_H
#define ADJOINT_NODAL_DISPLACEMENT_RESPONSE_FUNCTION_H

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------

// ------------------------------------------------------------------------------
// External includes
// ------------------------------------------------------------------------------

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "includes/define.h"
#include "adjoint_structural_response_function.h"


// ==============================================================================

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

/// Short class definition.
/** Detail class definition.

 */

//template<class TDenseSpace>

class AdjointNodalDisplacementResponseFunction : public AdjointStructuralResponseFunction
{
public:
    ///@name Type Definitions
    ///@{

    typedef AdjointStructuralResponseFunction BaseType;
    //typedef array_1d<double, 3> array_3d;
    typedef Element::DofsVectorType DofsVectorType;
    typedef Node<3>::Pointer PointTypePointer; //try to ensure that this response function also works for 2D
    typedef VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>> VariableComponentType;


    /// Pointer definition of AdjointNodalDisplacementResponseFunction
    KRATOS_CLASS_POINTER_DEFINITION(AdjointNodalDisplacementResponseFunction);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    AdjointNodalDisplacementResponseFunction(ModelPart& model_part, Parameters& responseSettings)
    : AdjointStructuralResponseFunction(model_part, responseSettings)
    {
        ModelPart& r_model_part = this->GetModelPart();

        // Get id of node where a displacement should be traced
        mIdOfTracedNode = responseSettings["traced_node"].GetInt();

        // Get the corresponding dof to the displacement which should be traced
        // By this response function e.g. DISPLACEMENT_X, ROTATION_X,...
        mTracedDofLabel = responseSettings["traced_dof"].GetString();


        // Get pointer to traced node
        mpTracedNode = r_model_part.pGetNode(mIdOfTracedNode);

        // Check if variable for traced dof is valid
        if( !( KratosComponents< VariableComponent< VectorComponentAdaptor<array_1d<double, 3> > > >::Has(mTracedDofLabel)) )
            KRATOS_ERROR << "Specified traced DOF is not availible. Specified DOF: " << mTracedDofLabel << std::endl;
        else
        {
            const VariableComponentType& rTRACED_DOF =
                KratosComponents<VariableComponentType>::Get(mTracedDofLabel);
            KRATOS_ERROR_IF_NOT( mpTracedNode->SolutionStepsDataHas(rTRACED_DOF) )
                << "Specified DOF is not availible at traced node." << std::endl;
        }

        // Check if variable for traced adjoint dof is valid
        if( !(KratosComponents< VariableComponent< VectorComponentAdaptor<array_1d<double, 3> > > >::Has(std::string("ADJOINT_") + mTracedDofLabel)) )
        {
            KRATOS_ERROR << "Specified traced adjoint DOF is not availible." << std::endl;
        }
        else
        {
            const VariableComponentType& rTRACED_ADJOINT_DOF =
                KratosComponents<VariableComponentType>::Get(std::string("ADJOINT_") + mTracedDofLabel);
            KRATOS_ERROR_IF_NOT( mpTracedNode->SolutionStepsDataHas(rTRACED_ADJOINT_DOF) )
                << "Specified adjoint DOF is not availible at traced node." << std::endl;
        }

        mDisplacementValue = 0.0;
        this->GetNeighboringElementPointer();
    }

    /// Destructor.
    virtual ~AdjointNodalDisplacementResponseFunction()
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    void Initialize() override
    {
        KRATOS_TRY;

        BaseType::Initialize();

        KRATOS_CATCH("");
    }

    // ==============================================================================
    void GetNeighboringElementPointer()
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
             << "No neighboring element is availible for the traced node." << std::endl;

        KRATOS_CATCH("");

    }

    // ==============================================================================
    double CalculateValue(ModelPart& rModelPart) override  
    {
        KRATOS_TRY;

        const VariableComponentType& rTRACED_DOF =
            KratosComponents<VariableComponentType>::Get(mTracedDofLabel);

        mDisplacementValue = rModelPart.Nodes()[(mpTracedNode->Id())].FastGetSolutionStepValue(rTRACED_DOF, 0);

        return mDisplacementValue;

        KRATOS_CATCH("");
    }

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
    // ==============================================================================
    void CalculateGradient(const Element& rAdjointElem, const Matrix& rAdjointMatrix,
                                   Vector& rResponseGradient,
                                   ProcessInfo& rProcessInfo) override
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

    // ==============================================================================
    void CalculateSensitivityGradient(Element& rAdjointElem,
                                      const Variable<double>& rVariable,
                                      const Matrix& rDerivativesMatrix,
                                      Vector& rResponseGradient,
                                      ProcessInfo& rProcessInfo) override
    {
          KRATOS_TRY

        if (rResponseGradient.size() != rDerivativesMatrix.size1())
            rResponseGradient.resize(rDerivativesMatrix.size1(), false);

        rResponseGradient.clear();


        KRATOS_CATCH("")
    }

    // ==============================================================================
    void CalculateSensitivityGradient(Condition& rAdjointCondition,
                                     const Variable<double>& rVariable,
                                     const Matrix& rDerivativesMatrix,
                                     Vector& rResponseGradient,
                                     ProcessInfo& rProcessInfo) override
    {
        KRATOS_TRY;

        if (rResponseGradient.size() != rDerivativesMatrix.size1())
            rResponseGradient.resize(rDerivativesMatrix.size1(), false);

        rResponseGradient.clear();

        KRATOS_CATCH("");
    }

    // ==============================================================================
    void CalculateSensitivityGradient(Element& rAdjointElem,
                                      const Variable<array_1d<double,3>>& rVariable,
                                      const Matrix& rDerivativesMatrix,
                                      Vector& rResponseGradient,
                                      ProcessInfo& rProcessInfo) override
    {
          KRATOS_TRY

        if (rResponseGradient.size() != rDerivativesMatrix.size1())
            rResponseGradient.resize(rDerivativesMatrix.size1(), false);

        rResponseGradient.clear();

          KRATOS_CATCH("")
    }

    // ==============================================================================
    void CalculateSensitivityGradient(Condition& rAdjointCondition,
                                      const Variable<array_1d<double,3>>& rVariable,
                                      const Matrix& rDerivativesMatrix,
                                      Vector& rResponseGradient,
                                      ProcessInfo& rProcessInfo) override
    {
        KRATOS_TRY;

        if (rResponseGradient.size() != rDerivativesMatrix.size1())
            rResponseGradient.resize(rDerivativesMatrix.size1(), false);

        rResponseGradient.clear();

        KRATOS_CATCH("");
    }



    // ==============================================================================

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
    double mDisplacementValue;
    int mIdOfTracedNode;
    std::string mTracedDofLabel;
    PointTypePointer  mpTracedNode;
    Element::Pointer mpNeighboringElement;


    //std::string m_traced_displacement_type;
    //std::string m_traced_direction;

    ///@}
///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    //      AdjointNodalDisplacementResponseFunction& operator=(SAdjointNodalDisplacementResponseFunction const& rOther);

    /// Copy constructor.
    //      AdjointNodalDisplacementResponseFunction(AdjointNodalDisplacementResponseFunction const& rOther);

    ///@}

}; // Class AdjointNodalDisplacementResponseFunction

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

} // namespace Kratos.

#endif // ADJOINT_NODAL_DISPLACEMENT_RESPONSE_FUNCTION_H
