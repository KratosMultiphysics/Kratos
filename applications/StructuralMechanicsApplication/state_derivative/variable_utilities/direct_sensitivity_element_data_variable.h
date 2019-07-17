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
    ///@addtogroup StructuralMechanicsApplication
    ///@{

    ///@name Kratos Classes
    ///@{

    /** \brief DirectSensitivityElementDataVariable
    *
    * This variable class describes parameters concerning parameters of elements (e.g. YOUNG_MODULUS).
    * It is designed to be used in direct sensitivity analysis.
    */

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
         
    std::vector<unsigned int> GetTracedElementId() override;

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
        std::vector<unsigned int> mTracedElementIdVector;        
    ///@}

    ///@name private Operators
    ///@{

    ///@}
    ///@name private Operations
    ///@{
                    
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
