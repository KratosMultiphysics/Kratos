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

#if !defined(SHELL_THIN_ADJOINT_ELEMENT_3D3N_H_INCLUDED )
#define  SHELL_THIN_ADJOINT_ELEMENT_3D3N_H_INCLUDED


// System includes

// External includes

// Project includes
#include "includes/element.h"
#include "custom_utilities/shell_cross_section.hpp"
#include "custom_utilities/shellt3_local_coordinate_system.hpp"
#include "utilities/quaternion.h"
#include "custom_elements/shell_thin_element_3D3N.hpp"

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

/** \brief ShellThinAdjointElement3D3N
 *
 * This element is inherited from ShellThinElement3D3N. It is the corresponding 
 * element to it and is used for solving the adjoint problem and for computing sensitivities.
 */
class ShellThinAdjointElement3D3N : public ShellThinElement3D3N 
{
public:

    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(ShellThinAdjointElement3D3N);

     typedef Element::PropertiesType PropertiesType;

     typedef Element::DofsArrayType DofsArrayType;

    ///@}

    ///@name Classes
    ///@{
    ///@}

    ///@name Life Cycle
    ///@{

    ShellThinAdjointElement3D3N(IndexType NewId,
                         GeometryType::Pointer pGeometry,
                         bool NLGeom = false);

    ShellThinAdjointElement3D3N(IndexType NewId,
                         GeometryType::Pointer pGeometry,
                         PropertiesType::Pointer pProperties,
                         bool NLGeom = false);

    ShellThinAdjointElement3D3N(IndexType NewId,
                         GeometryType::Pointer pGeometry,
                         PropertiesType::Pointer pProperties,
                         CoordinateTransformationBasePointerType pCoordinateTransformation);

    ~ShellThinAdjointElement3D3N() override;

    ///@}

    ///@name Operations
    ///@{

    // Basic
    
    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override;

    void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo) override;

    void GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& CurrentProcessInfo) override;

    int Check(const ProcessInfo& rCurrentProcessInfo) override;

    void GetValuesVector(Vector& values, int Step = 0) override;

    double GetDisturbanceMeasureCorrectionFactor(const Variable<double>& rVariable);

    double GetDisturbanceMeasureCorrectionFactor(const Variable<array_1d<double,3>>& rDesignVariable);

    void CalculateSensitivityMatrix(const Variable<double>& rDesignVariable, Matrix& rOutput, 
                                            const ProcessInfo& rCurrentProcessInfo) override;
    
    void CalculateSensitivityMatrix(const Variable<array_1d<double,3>>& rDesignVariable, Matrix& rOutput, 
                                            const ProcessInfo& rCurrentProcessInfo) override;

    void Calculate(const Variable<Vector >& rVariable, Vector& rOutput,
                           const ProcessInfo& rCurrentProcessInfo) override;                                         

    void Calculate(const Variable<Matrix >& rVariable, Matrix& rOutput,
                           const ProcessInfo& rCurrentProcessInfo) override;   

    void CalculateStressDisplacementDerivative(const Variable<Vector>& rStressVariable,
                                    Matrix& rOutput, const ProcessInfo& rCurrentProcessInfo);    

    void CalculateStressDesignVariableDerivative(const Variable<double>& rDesignVariable, const Variable<Vector>& rStressVariable,
                                        Matrix& rOutput, const ProcessInfo& rCurrentProcessInfo);
    
    void CalculateStressDesignVariableDerivative(const Variable<array_1d<double,3>>& rDesignVariable, 
                                            const Variable<Vector>& rStressVariable,
                                             Matrix& rOutput, const ProcessInfo& rCurrentProcessInfo);                       

    void CalculateLeftHandSide( MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo) override;    
                      
    // Results calculation on integration points

    void CalculateOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rOutput,
                                        const ProcessInfo& rCurrentProcessInfo) override;

    void GetValueOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo) override;

    ///@}

    ///@name Public specialized Access - Temporary
    ///@{
    ///@}

protected:

    ///@name Protected Lyfe Cycle
    ///@{

    /**
     * Protected empty constructor
     */
    ShellThinAdjointElement3D3N() : ShellThinElement3D3N()
    {
    }

    ///@}

private:

    ///@name Private Classes
    ///@{
    ///@}

    ///@name Private Operations
    ///@{
    ///@}

    ///@name Static Member Variables
    ///@{
    ///@}

    ///@name Member Variables
    ///@{
    ///@}

    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override;

    void load(Serializer& rSerializer) override;

    ///@}

    ///@name Private  Access
    ///@{
    ///@}

    ///@name Private Inquiry
    ///@{
    ///@}

    ///@name Un accessible methods
    ///@{
    ///@}

};

}
#endif // SHELL_THIN_ADJOINT_ELEMENT_3D3N_H_INCLUDED
