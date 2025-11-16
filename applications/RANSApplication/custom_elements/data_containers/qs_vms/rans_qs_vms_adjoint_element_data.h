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

#if !defined(KRATOS_RANS_QS_VMS_ADJOINT_ELEMENT_DATA_H)
#define KRATOS_RANS_QS_VMS_ADJOINT_ELEMENT_DATA_H

// System includes

// External includes

// Project includes

// Application includes
#include "custom_elements/data_containers/fluid_adjoint_derivatives.h"
#include "custom_elements/data_containers/qs_vms/qs_vms_derivative_utilities.h"
#include "custom_elements/data_containers/qs_vms/qs_vms_residual_derivatives.h"
#include "custom_elements/data_containers/qs_vms/rans_qs_vms_derivative_utilities.h"

namespace Kratos
{
///@name Kratos Classes
///@{

template <
    unsigned int TDim,
    unsigned int TNumNodes,
    template<unsigned int, unsigned int> class TVelocityDerivative,
    template<unsigned int, unsigned int> class TShapeDerivative,
    class TTurbulenceModelElementData1,
    class TTurbulenceModelElementData2>
class RansQSVMSAdjointElementData
{
private:
    ///@name Private Type Definitions
    ///@{

    using IndexType = std::size_t;

    using TResidualsDerivatives = QSVMSResidualDerivatives<TDim, TNumNodes>;

    using Data = typename TResidualsDerivatives::Data;

    using ResidualsContributions = typename TResidualsDerivatives::ResidualsContributions;

    using PressureDerivativeContributions = typename TResidualsDerivatives::template VariableDerivatives<typename QSVMSDerivativeUtilities<TDim>::template PressureDerivative<TNumNodes>>;

    template<unsigned int TDirectionIndex>
    using VelocityDerivativeContributions = typename TResidualsDerivatives::template VariableDerivatives<TVelocityDerivative<TNumNodes, TDirectionIndex>>;

    template<unsigned int TDirectionIndex>
    using ShapeDerivatives = typename TResidualsDerivatives::template VariableDerivatives<TShapeDerivative<TNumNodes, TDirectionIndex>>;

    template<unsigned int TDirectionIndex>
    using AccelerationDerivativeContributions = typename TResidualsDerivatives::template SecondDerivatives<TDirectionIndex>;

    using TurbulenceModelVariable1DerivativeContributions = typename TResidualsDerivatives::template VariableDerivatives<typename RansQSVMSDerivativeUtilities<TDim>::template TurbulenceVariableDerivative<TNumNodes, TTurbulenceModelElementData1>>;

    using TurbulenceModelVariable2DerivativeContributions = typename TResidualsDerivatives::template VariableDerivatives<typename RansQSVMSDerivativeUtilities<TDim>::template TurbulenceVariableDerivative<TNumNodes, TTurbulenceModelElementData2>>;

    static constexpr IndexType ElementDataContainerIndex = 0;

    static constexpr IndexType ResidualColumnOffset = 0;

    ///@}

public:
    ///@name Template Type Definitions
    ///@{

    using Residual = CalculationDataContainers<
                        std::tuple<
                            Data>,
                        std::tuple<
                            SubAssembly<ResidualsContributions, ElementDataContainerIndex, 0, ResidualColumnOffset>>
                        >;

    using ResidualStateVariableFirstDerivatives = typename std::conditional<
                                                        TDim == 2,
                                                        CalculationDataContainers<
                                                            std::tuple<
                                                                Data>,
                                                            std::tuple<
                                                                SubAssembly<VelocityDerivativeContributions<0>,              ElementDataContainerIndex, 0, ResidualColumnOffset>,
                                                                SubAssembly<VelocityDerivativeContributions<1>,              ElementDataContainerIndex, 1, ResidualColumnOffset>,
                                                                SubAssembly<PressureDerivativeContributions,                 ElementDataContainerIndex, 2, ResidualColumnOffset>,
                                                                SubAssembly<TurbulenceModelVariable1DerivativeContributions, ElementDataContainerIndex, 3, ResidualColumnOffset>,
                                                                SubAssembly<TurbulenceModelVariable2DerivativeContributions, ElementDataContainerIndex, 4, ResidualColumnOffset>>
                                                            >,
                                                        CalculationDataContainers<
                                                            std::tuple<
                                                                Data>,
                                                            std::tuple<
                                                                SubAssembly<VelocityDerivativeContributions<0>,              ElementDataContainerIndex, 0, ResidualColumnOffset>,
                                                                SubAssembly<VelocityDerivativeContributions<1>,              ElementDataContainerIndex, 1, ResidualColumnOffset>,
                                                                SubAssembly<VelocityDerivativeContributions<2>,              ElementDataContainerIndex, 2, ResidualColumnOffset>,
                                                                SubAssembly<PressureDerivativeContributions,                 ElementDataContainerIndex, 3, ResidualColumnOffset>,
                                                                SubAssembly<TurbulenceModelVariable1DerivativeContributions, ElementDataContainerIndex, 4, ResidualColumnOffset>,
                                                                SubAssembly<TurbulenceModelVariable2DerivativeContributions, ElementDataContainerIndex, 5, ResidualColumnOffset>>
                                                            >
                                                        >::type;

    using ResidualStateVariableSecondDerivatives = typename std::conditional<
                                                        TDim == 2,
                                                        CalculationDataContainers<
                                                            std::tuple<
                                                                Data>,
                                                            std::tuple<
                                                                SubAssembly<AccelerationDerivativeContributions<0>, ElementDataContainerIndex, 0, ResidualColumnOffset>,
                                                                SubAssembly<AccelerationDerivativeContributions<1>, ElementDataContainerIndex, 1, ResidualColumnOffset>,
                                                                SubAssembly<ZeroDerivatives<TNumNodes, 5>,          ElementDataContainerIndex, 2, ResidualColumnOffset>,
                                                                SubAssembly<ZeroDerivatives<TNumNodes, 5>,          ElementDataContainerIndex, 3, ResidualColumnOffset>,
                                                                SubAssembly<ZeroDerivatives<TNumNodes, 5>,          ElementDataContainerIndex, 4, ResidualColumnOffset>>
                                                            >,
                                                        CalculationDataContainers<
                                                            std::tuple<
                                                                Data>,
                                                            std::tuple<
                                                                SubAssembly<AccelerationDerivativeContributions<0>, ElementDataContainerIndex, 0, ResidualColumnOffset>,
                                                                SubAssembly<AccelerationDerivativeContributions<1>, ElementDataContainerIndex, 1, ResidualColumnOffset>,
                                                                SubAssembly<AccelerationDerivativeContributions<2>, ElementDataContainerIndex, 2, ResidualColumnOffset>,
                                                                SubAssembly<ZeroDerivatives<TNumNodes, 6>,          ElementDataContainerIndex, 3, ResidualColumnOffset>,
                                                                SubAssembly<ZeroDerivatives<TNumNodes, 6>,          ElementDataContainerIndex, 4, ResidualColumnOffset>,
                                                                SubAssembly<ZeroDerivatives<TNumNodes, 6>,          ElementDataContainerIndex, 5, ResidualColumnOffset>>
                                                            >
                                                        >::type;

    using ResidualShapeDerivatives = typename std::conditional<
                                            TDim == 2,
                                            CalculationDataContainers<
                                                std::tuple<
                                                    Data>,
                                                std::tuple<
                                                    SubAssembly<ShapeDerivatives<0>, ElementDataContainerIndex, 0, ResidualColumnOffset>,
                                                    SubAssembly<ShapeDerivatives<1>, ElementDataContainerIndex, 1, ResidualColumnOffset>>
                                                >,
                                            CalculationDataContainers<
                                                std::tuple<
                                                    Data>,
                                                std::tuple<
                                                    SubAssembly<ShapeDerivatives<0>, ElementDataContainerIndex, 0, ResidualColumnOffset>,
                                                    SubAssembly<ShapeDerivatives<1>, ElementDataContainerIndex, 1, ResidualColumnOffset>,
                                                    SubAssembly<ShapeDerivatives<2>, ElementDataContainerIndex, 2, ResidualColumnOffset>>
                                                >
                                            >::type;

    ///@}
    ///@name Static Operations
    ///@{

    static void Check(
        const Element& rElement,
        const ProcessInfo& rProcessInfo);

    static std::vector<const Variable<double>*> GetDofVariablesList();

    static GeometryData::IntegrationMethod GetIntegrationMethod();

    ///@}
};
} // namespace Kratos
#endif // KRATOS_RANS_QS_VMS_ADJOINT_ELEMENT_DATA_H