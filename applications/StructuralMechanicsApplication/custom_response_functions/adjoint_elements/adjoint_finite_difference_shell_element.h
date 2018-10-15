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
    ///@}

    ///@name Classes
    ///@{
    ///@}

    ///@name Life Cycle
    ///@{
    AdjointFiniteDifferencingShellElement() : AdjointFiniteDifferencingBaseElement()
    {
    }

    AdjointFiniteDifferencingShellElement(Element::Pointer pPrimalElement);

    ~AdjointFiniteDifferencingShellElement() override;

    ///@}

    ///@name Operations
    ///@{

    // Basic

    int Check(const ProcessInfo& rCurrentProcessInfo) override;

    ///@}

    ///@name Public specialized Access - Temporary
    ///@{
    ///@}

protected:

    ///@name Protected Lyfe Cycle
    ///@{

    ///@}

private:

    ///@name Private Classes
    ///@{
    ///@}

    ///@name Private Operations
    ///@{

    void CheckVariables();
    void CheckDofs();
    void CheckProperties(const ProcessInfo& rCurrentProcessInfo);
    void CheckSpecificProperties();

    double GetPerturbationSizeModificationFactor(const Variable<array_1d<double,3>>& rDesignVariable) override;

    ///@}
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
#endif // ADJOINT_FINITE_DIFFERENCE_SHELL_ELEMENT_H_INCLUDED
