// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Armin Geiser, https://github.com/armingeiser
//

#if !defined(ADJOINT_FINITE_DIFFERENCE_SHELL_ELEMENT_H_INCLUDED )
#define  ADJOINT_FINITE_DIFFERENCE_SHELL_ELEMENT_H_INCLUDED


// System includes

// External includes

// Project includes
#include "adjoint_finite_difference_base_element.h"
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

/** \brief AdjointFiniteDifferencingShellElement
 *
 * This element is inherited from AdjointFiniteDifferencingBaseElement.
 * It overwrites some functions necessary to do proper finite differencing with shell elements
 */
class AdjointFiniteDifferencingShellElement : public AdjointFiniteDifferencingBaseElement
{
public:

    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(AdjointFiniteDifferencingShellElement);

     typedef Element::PropertiesType PropertiesType;

     typedef Element::DofsArrayType DofsArrayType;

    ///@}

    ///@name Classes
    ///@{
    ///@}

    ///@name Life Cycle
    ///@{

    AdjointFiniteDifferencingShellElement(IndexType NewId,
                         GeometryType::Pointer pGeometry);

    AdjointFiniteDifferencingShellElement(IndexType NewId,
                         GeometryType::Pointer pGeometry,
                         PropertiesType::Pointer pProperties,
                         Element::Pointer pPrimalElement);

    ~AdjointFiniteDifferencingShellElement() override;

    ///@}

    ///@name Operations
    ///@{

    // Basic

    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties,
                Element::Pointer pPrimalElement) const;

    // TODO Element::Pointer Clone (IndexType NewId, NodesArrayType const& ThisNodes) const override;

    void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo) override;

    void GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& CurrentProcessInfo) override;

    void GetValuesVector(Vector& values, int Step = 0) override;

    void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
                                       ProcessInfo& rCurrentProcessInfo) override
    {
        Vector dummy;
        mpPrimalElement->CalculateLocalSystem(rLeftHandSideMatrix, dummy, rCurrentProcessInfo);
        // TODO HACK necessary because shell doe not LHS...!!
        // mpPrimalElement->CalculateLeftHandSide(rLeftHandSideMatrix,
        //                                       rCurrentProcessInfo);
    }

    // TODO add functions from element.h line 641 - 710

    void Calculate(const Variable<Vector >& rVariable, Vector& rOutput,
                           const ProcessInfo& rCurrentProcessInfo) override;

    void Calculate(const Variable<Matrix >& rVariable, Matrix& rOutput,
                           const ProcessInfo& rCurrentProcessInfo) override;

    // TODO evaluate if other Calculate functions are necessary

    // TODO add functions from element.h line 744 - 882
    // Results calculation on integration points

    void CalculateOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rOutput,
                                        const ProcessInfo& rCurrentProcessInfo) override;

    void GetValueOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo) override;


    int Check(const ProcessInfo& rCurrentProcessInfo) override;

    // Sensitivity functions

    void CalculateSensitivityMatrix(const Variable<double>& rDesignVariable, Matrix& rOutput,
                                            const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateSensitivityMatrix(const Variable<array_1d<double,3>>& rDesignVariable, Matrix& rOutput,
                                            const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateStressDisplacementDerivative(const Variable<Vector>& rStressVariable,
                                    Matrix& rOutput, const ProcessInfo& rCurrentProcessInfo);

    void CalculateStressDesignVariableDerivative(const Variable<double>& rDesignVariable, const Variable<Vector>& rStressVariable,
                                        Matrix& rOutput, const ProcessInfo& rCurrentProcessInfo);

    void CalculateStressDesignVariableDerivative(const Variable<array_1d<double,3>>& rDesignVariable,
                                            const Variable<Vector>& rStressVariable,
                                             Matrix& rOutput, const ProcessInfo& rCurrentProcessInfo);

    ///@}

    ///@name Public specialized Access - Temporary
    ///@{
    ///@}

protected:

    ///@name Protected Lyfe Cycle
    ///@{

    /**
     * Protected empty constructor TODO needed?
     */
    AdjointFiniteDifferencingShellElement() : AdjointFiniteDifferencingBaseElement()
    {
    }

    ///@}

private:

    ///@name Private Classes
    ///@{
    ///@}

    ///@name Private Operations
    ///@{

    double GetDisturbanceMeasureCorrectionFactor(const Variable<double>& rVariable);

    double GetDisturbanceMeasureCorrectionFactor(const Variable<array_1d<double,3>>& rDesignVariable);

    ///@}

    ///@name Static Member Variables
    ///@{
    ///@}

    ///@name Member Variables
    ///@{

    /**
     * pointer to the primal element
     */
    ShellThinElement3D3N::Pointer mpPrimalShellElement;

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
#endif // ADJOINT_FINITE_DIFFERENCE_SHELL_ELEMENT_H_INCLUDED
