// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:     BSD License
//           license: structural_mechanics_application/license.txt
//
//  Main authors:    Martin Fusseder, https://github.com/MFusseder
//


#include "adjoint_finite_difference_cable_element_3D2N.h"
#include "structural_mechanics_application_variables.h"
#include "custom_response_functions/response_utilities/stress_response_definitions.h"
#include "custom_elements/cable_element_3D2N.hpp"
#include "custom_utilities/structural_mechanics_element_utilities.h"


namespace Kratos
{

template <class TPrimalElement>
void AdjointFiniteDifferenceCableElement<TPrimalElement>::InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // to unsure that flag 'mIsCompressed' is set correctly in primal cable element
    Vector RHS;
    this->pGetPrimalElement()->CalculateRightHandSide(RHS, rCurrentProcessInfo);
    BaseType::InitializeSolutionStep(rCurrentProcessInfo);

    KRATOS_CATCH("")
}

template <class TPrimalElement>
void AdjointFiniteDifferenceCableElement<TPrimalElement>::CalculateStressDisplacementDerivative(const Variable<Vector>& rStressVariable,
                                    Matrix& rOutput, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    if(rStressVariable == STRESS_ON_GP)
    {

        BaseType::CalculateStressDisplacementDerivative(rStressVariable, rOutput, rCurrentProcessInfo);

        /*const double l_0 = StructuralMechanicsElementUtilities::CalculateReferenceLength3D2N(*this);
        const double l = StructuralMechanicsElementUtilities::CalculateCurrentLength3D2N(*this);

        if (l > l_0) {
            BaseType::CalculateStressDisplacementDerivative(rStressVariable, rOutput, rCurrentProcessInfo);
        }
        else {
            const SizeType num_nodes = this->mpPrimalElement->GetGeometry().PointsNumber();
            const SizeType dimension = this->mpPrimalElement->GetGeometry().WorkingSpaceDimension();
            const SizeType num_dofs = num_nodes * dimension;
            const SizeType num_GP = (this->mpPrimalElement->GetGeometry().IntegrationPoints()).size();
            if ( (rOutput.size1() != num_dofs) || (rOutput.size2() != num_GP ) ) {
                rOutput.resize(num_dofs, num_GP);
            }
            noalias(rOutput) = ZeroMatrix(num_dofs, num_GP);
        }*/
    }
    else
        KRATOS_ERROR << "Stress displacement derivative only available for Gauss-points quantities!" << std::endl;

    KRATOS_CATCH("")
}

template <class TPrimalElement>
void AdjointFiniteDifferenceCableElement<TPrimalElement>::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseType);
}

template <class TPrimalElement>
void AdjointFiniteDifferenceCableElement<TPrimalElement>::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseType);
}

template class AdjointFiniteDifferenceCableElement<CableElement3D2N>;

} // namespace Kratos.


