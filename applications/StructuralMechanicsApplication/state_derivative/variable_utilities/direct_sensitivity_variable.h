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

#ifndef DIRECT_SENSITIVITY_VARIABLE_H
#define DIRECT_SENSITIVITY_VARIABLE_H
// System includes

// External includes

// Project includes
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"
#include "structural_mechanics_application_variables.h"
#include "includes/define.h"
#include "includes/element.h"
#include "includes/condition.h"
#include "includes/process_info.h"




namespace Kratos
{
///@addtogroup StructuralMechanicsApplication
///@{

///@name Kratos Classes
///@{

/** \brief DirectSensitivityVariable
*
* This is the response base class for responses in structural mechanics.
* It is designed to be used in adjoint sensitivity analysis.
*/
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) DirectSensitivityVariable
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(DirectSensitivityVariable);

    typedef VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>> VariableComponentType;

    typedef Variable<array_1d<double, 3>> VariableWithComponentsType;

    typedef ModelPart::NodesContainerType NodesContainerType;

    typedef ModelPart::ElementsContainerType ElementsContainerType;

    typedef ModelPart::ConditionsContainerType ConditionsContainerType;

    typedef Vector VectorType;

    typedef std::size_t IndexType;

    typedef std::size_t SizeType;

    

    ///@}
    ///@name Pointer Definitions

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    DirectSensitivityVariable(ModelPart& rModelPart, Parameters ResponseSettings);

    /// Destructor.
    virtual ~DirectSensitivityVariable();

    ///@}
    ///@name Operators
    ///@{
        
    ///@}
    ///@name Operations
    ///@{

    
    virtual void Initialize();

    virtual void InitializeSolutionStep(){};

    virtual void FinalizeSolutionStep(){};
       
    virtual void CalculatePseudoLoadVector(Element& rDirectElement, const Matrix& rRHS, Vector& rPseudoLoadVector, const ProcessInfo& rProcessInfo);

    virtual void CalculatePseudoLoadVector(Condition& rDirectCondition, const Matrix& rLHS, Vector& rPseudoLoadVector, const ProcessInfo& rProcessInfo);

    virtual void PerturbDesignVariable(Element& rDirectElement);
    
    virtual void UnperturbDesignVariable(Element& rDirectElement);        
       
    virtual std::string GetDesignVariableName();

    virtual std::string GetDesignVariableType();

    virtual unsigned int GetTracedElementId();

    virtual double GetPerturbationSize();

    

    ///@}

protected:
    ///@name protected member Variables
    ///@{

    ModelPart& mrModelPart;
    double mDelta;
    std::string mDesignVariableName;
    std::string mVariableType;
    unsigned int mTracedElement;

     
    ///@}
    ///@name protected Operators
    ///@{

    ///@}
    ///@name protected Operations
    ///@{

    ///@} 

private:
    ///@name private Member Variables
    ///@{

    /*std::string mSensitivityModelPartName;
    unsigned int mGradientMode;
    */

    ///@}
    ///@name private Operators
    ///@{

    ///@}
    ///@name private Operations
    ///@{
    virtual Variable<double> ReadScalarDesignVariables(std::string const& rVariableName);

    virtual Variable<array_1d<double,3>> ReadVectorDesignVariables(std::string const& rVariableName);    
    ///@} 
        
};

///@} // Kratos Classes

///@} //Structural Mechanics Application group

} /* namespace Kratos.*/

#endif /* DIRECT_SENSITIVITY_VARIABLE_H defined */
