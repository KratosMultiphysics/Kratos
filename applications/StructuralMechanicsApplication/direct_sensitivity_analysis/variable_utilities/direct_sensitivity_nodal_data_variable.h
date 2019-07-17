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

#ifndef DIRECT_SENSITIVITY_NODAL_DATA_VARIABLE_H
#define DIRECT_SENSITIVITY_NODAL_DATA_VARIABLE_H

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

    /** \brief DirectSensitivityNodalDataVariable
    *
    * This variable class describes parameters concerning the coordinates of nodes.
    * It is designed to be used in direct sensitivity analysis.
    */

class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) DirectSensitivityNodalDataVariable : public DirectSensitivityVariable
{
public:
    ///@name Type Definitions
    ///@{

    typedef DirectSensitivityVariable BaseType;

    ///@}
    ///@name Pointer Definitions
    /// Pointer definition of DirectSensitivityNodalDataVariable
    KRATOS_CLASS_POINTER_DEFINITION(DirectSensitivityNodalDataVariable);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    DirectSensitivityNodalDataVariable(ModelPart& rModelPart, Parameters VariableSettings);

    /// Destructor.
    ~DirectSensitivityNodalDataVariable();
    
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

    void ExtractDataFromDerivativeMatrix(Element& rDirectElement,
                                        Matrix& rExtractedDerivativeMatrix,
                                        const Matrix& rDerivativeMatrix) override;

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
    
    std::vector<unsigned int> mNeighboringElementIdVector;
    std::vector<unsigned int> mNeighboringConditionIdVector;
    std::string mCoordinateDirection;

    ///@}

    ///@name private Operators
    ///@{

    ///@}
    ///@name private Operations
    ///@{    

    void GetNeighboringElementIdVector(std::vector<unsigned int>& rTracedElementIdVector);

    void GetNeighboringConditionIdVector(std::vector<unsigned int>& rNeighboringConditionIdVector);

    unsigned int GetCoordinateDirection();                
    ///@}

}; // Class DirectSensitivityNodalDataVariable

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

} // namespace Kratos.

#endif // DIRECT_SENSITIVITY_NODAL_DATA_VARIABLE_H_INCLUDED
