//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

#if !defined(KRATOS_VMS_K_EPSILON_ADJOINT_ELEMENT_DATA_H_INCLUDED)
#define KRATOS_VMS_K_EPSILON_ADJOINT_ELEMENT_DATA_H_INCLUDED

// System includes

// Project includes

// Application includes
#include "custom_elements/convection_diffusion_reaction_residual_based_flux_corrected_adjoint_utilities.h"
#include "custom_elements/data_containers/k_epsilon/epsilon_adjoint_element_data.h"
#include "custom_elements/data_containers/k_epsilon/k_adjoint_element_data.h"
// #include "custom_elements/data_containers/vms_adjoint_element_data.h"

namespace Kratos
{
template <unsigned int TDim, unsigned int TNumNodes>
class VMSRFCKEpsilonAdjintElementData
{
    template <class TStateDerivativesType>
    using stabilization_state_derivatives =
        ConvectionDiffusionReactionResidualBasedFluxCorrectedAdjointUtilities::StabilizationStateDerivatives<TDim, TNumNodes, TStateDerivativesType>;

    constexpr static unsigned int TEquation1Offset = TDim + 1;
    constexpr static unsigned int TEquation2Offset = TDim + 2;

public:
    class TurbulenceEquation1
        : public stabilization_state_derivatives<KEpsilonAdjointElementData::KAdjointStateDerivatives<TDim, TNumNodes>>
    {
    public:
        using BaseType =
            stabilization_state_derivatives<KEpsilonAdjointElementData::KAdjointStateDerivatives<TDim, TNumNodes>>;

        class FirstDerivatives
        {
        public:
            using Data = typename BaseType::FirstDerivativesData;
            using VelocityPressure = typename BaseType::FirstDerivatives<
                typename KEpsilonAdjointElementData::KAdjointStateDerivatives<TDim, TNumNodes>::VelocityDerivatives, TEquation1Offset, 0>;
            using TurbulenceDof1 = typename BaseType::FirstDerivatives<
                typename KEpsilonAdjointElementData::KAdjointStateDerivatives<TDim, TNumNodes>::KDerivatives, TEquation1Offset, TEquation1Offset>;
            using TurbulenceDof2 = typename BaseType::FirstDerivatives<
                typename KEpsilonAdjointElementData::KAdjointStateDerivatives<TDim, TNumNodes>::EpsilonDerivatives, TEquation1Offset, TEquation2Offset>;
        };

        using SecondDerivatives = typename BaseType::SecondDerivatives<TEquation1Offset>;
    };

    class TurbulenceEquation2
        : public stabilization_state_derivatives<KEpsilonAdjointElementData::EpsilonAdjointStateDerivatives<TDim, TNumNodes>>
    {
    public:
        using BaseType =
            stabilization_state_derivatives<KEpsilonAdjointElementData::EpsilonAdjointStateDerivatives<TDim, TNumNodes>>;

        class FirstDerivatives
        {
        public:
            using Data = typename BaseType::FirstDerivativesData;;

            using VelocityPressure = typename BaseType::FirstDerivatives<
                typename KEpsilonAdjointElementData::EpsilonAdjointStateDerivatives<TDim, TNumNodes>::VelocityDerivatives, TEquation2Offset, 0>;
            using TurbulenceDof1 = typename BaseType::FirstDerivatives<
                typename KEpsilonAdjointElementData::EpsilonAdjointStateDerivatives<TDim, TNumNodes>::KDerivatives, TEquation2Offset, TEquation1Offset>;
            using TurbulenceDof2 = typename BaseType::FirstDerivatives<
                typename KEpsilonAdjointElementData::EpsilonAdjointStateDerivatives<TDim, TNumNodes>::EpsilonDerivatives, TEquation2Offset, TEquation2Offset>;
        };

        using SecondDerivatives = typename BaseType::SecondDerivatives<TEquation2Offset>;
    };
};
} // namespace Kratos

#endif // KRATOS_VMS_K_EPSILON_ADJOINT_ELEMENT_DATA_H_INCLUDED