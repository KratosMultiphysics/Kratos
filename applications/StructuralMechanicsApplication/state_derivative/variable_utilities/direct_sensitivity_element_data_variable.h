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

#ifndef DIRECT_SENSITIVITY_ELEMENT_DATA_VARIABLE_H
#define DIRECT_SENSITIVITY_ELEMENT_DATA_VARIABLE_H

// System includes

// External includes

// Project includes
#include "direct_sensitivity_variable.h"


namespace Kratos
{

class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) DirectSensitivityElementDataVariable : public DirectSensitivityVariable
{
public:
    ///@name Type Definitions
    ///@{

    typedef DirectSensitivityVariable BaseType;

    ///@}
    ///@name Pointer Definitions
    /// Pointer definition of DirectSensitivityElementDataVariable
    KRATOS_CLASS_POINTER_DEFINITION(DirectSensitivityElementDataVariable);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    DirectSensitivityElementDataVariable(ModelPart& rModelPart, Parameters VariableSettings);

    /// Destructor.
    ~DirectSensitivityElementDataVariable();
    
    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    void Initialize() override;

    void CalculatePseudoLoadVector(Element& rDirectElem, const Matrix& rRHS, Vector& rPseudoLoadVector, const ProcessInfo& rProcessInfo) override;

    void CalculatePseudoLoadVector(Condition& rDirectCondition, const Matrix& rLHS, Vector& rPseudoLoadVector, const ProcessInfo& rProcessInfo) override;    

    void PerturbDesignVariable(Element& rElement) override;
    
    void UnperturbDesignVariable(Element& rElement) override;
      
    unsigned int GetTracedElementId() override;

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
    ///@name protected member Variables
    ///@{
    ///@}

    ///@name protected Operators
    ///@{

    ///@}
    ///@name protected Operations
    ///@{

    ///@}
    
private:
    ///@name private member Variables
    ///@{
        unsigned int mIdTracedElement;
        Element mTracedElement;   
        Element::Pointer mpTracedElement;        
        ProcessInfo mProcessInfo;

    ///@}

    ///@name private Operators
    ///@{

    ///@}
    ///@name private Operations
    ///@{
        Variable<double> ReadScalarDesignVariables(std::string const& rVariableName) override;

        Variable<array_1d<double,3>> ReadVectorDesignVariables(std::string const& rVariableName) override;        
        
    ///@}

}; // Class DirectSensitivityElementDataVariable

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

} // namespace Kratos.

#endif // DIRECT_SENSITIVITY_ELEMENT_DATA_VARIABLE_H_INCLUDED
