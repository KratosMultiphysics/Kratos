// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Kevin Braun, https://github.com/MFusseder
//

#ifndef DIRECT_SENSITIVITY_NODAL_DISPLACEMENT_RESPONSE_FUNCTION_H
#define DIRECT_SENSITIVITY_NODAL_DISPLACEMENT_RESPONSE_FUNCTION_H

// System includes

// External includes

// Project includes
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"
#include "structural_mechanics_application_variables.h"
#include "direct_sensitivity_response_function.h"


namespace Kratos
{
///@addtogroup StructuralMechanicsApplication
///@{

///@name Kratos Classes
///@{

/** \brief DirectSensitivityNodalDisplacementResponseFunction
*
* This is the response class for nodal displacements used in the direct sensitivity analysis.
*/
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) DirectSensitivityNodalDisplacementResponseFunction : public DirectSensitivityResponseFunction
{
public:
    ///@name Type Definitions
    ///@{
    
    typedef DirectSensitivityResponseFunction BaseType;

    KRATOS_CLASS_POINTER_DEFINITION(DirectSensitivityNodalDisplacementResponseFunction);

    typedef std::size_t IndexType;

    typedef std::size_t SizeType;

    typedef Element::DofsVectorType DofsVectorType;
    typedef Node<3>::Pointer PointTypePointer;
    typedef VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>> VariableComponentType;

    ///@}
    ///@name Pointer Definitions

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    DirectSensitivityNodalDisplacementResponseFunction(ModelPart& rModelPart, Parameters ResponseSettings, std::string& ResponseVariableName);

    /// Destructor.
     ~DirectSensitivityNodalDisplacementResponseFunction();
    
    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    void Initialize() override;    

    void CalculateGradient(Node<3>& rNode,                            
                            Variable<array_1d<double, 3>> const& rStressVariable,
                            std::vector<array_1d<double, 3>>& rOutput, 
                            const ProcessInfo& rProcessInfo) override;

    void CalculatePartialSensitivity(Node<3>& rNode, 
                            DirectSensitivityVariable& rDesignVariable,
                            Variable<array_1d<double, 3>> const& rStressVariable, 
                            array_1d<double, 3>& rOutput, 
                            const ProcessInfo& rProcessInfo) override;

    std::string GetEvaluationFlag() override;
    
    ///@}

protected:
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    ///@}

private:
    ///@name private Member Variables
    ///@{
    
    
    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{
    
    

    ///@}
};

///@} // Kratos Classes

///@} //Structural Mechanics Application group

} /* namespace Kratos.*/

#endif /* DIRECT_SENSITIVITY_NODAL_DISPLACEMENT_RESPONSE_FUNCTION_H defined */
