// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:     BSD License
//           license: structural_mechanics_application/license.txt
//
//  Main authors:    Armin Geiser, https://github.com/armingeiser
//


#if !defined(ADJOINT_FINITE_DIFFERENCE_SMALL_DISPLACEMENT_ELEMENT_H_INCLUDED )
#define  ADJOINT_FINITE_DIFFERENCE_SMALL_DISPLACEMENT_ELEMENT_H_INCLUDED

#include "adjoint_finite_difference_base_element.h"

namespace Kratos
{


/** \brief AdjointFiniteDifferenceSmallDisplacementElement
 *
 * This element is a wrapper for a primal small displacement element. It is responsible to deliver local stresses and
 * the stress displacement derivative. It is designed to be used in adjoint
 * sensitivity analysis.
 */
class AdjointFiniteDifferenceSmallDisplacementElement : public AdjointFiniteDifferencingBaseElement
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(AdjointFiniteDifferenceSmallDisplacementElement);

    AdjointFiniteDifferenceSmallDisplacementElement(): AdjointFiniteDifferencingBaseElement()
    {
    }

    AdjointFiniteDifferenceSmallDisplacementElement(Element::Pointer pPrimalElement);

    ~AdjointFiniteDifferenceSmallDisplacementElement() override;

    int Check(const ProcessInfo& rCurrentProcessInfo) override;

private:

    friend class Serializer;
    void save(Serializer& rSerializer) const override;
    void load(Serializer& rSerializer) override;
};


}

#endif
