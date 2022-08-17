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

#if !defined(KRATOS_QS_VMS_ADJOINT_ELEMENT_DATA_H)
#define KRATOS_QS_VMS_ADJOINT_ELEMENT_DATA_H

// System includes

// External includes

// Project includes

// Application includes
#include "custom_elements/data_containers/derivatives.h"
#include "custom_elements/data_containers/qs_vms/qs_vms_derivative_utilities.h"
#include "custom_elements/data_containers/qs_vms/qs_vms_residual_derivatives.h"

namespace Kratos
{
///@name Kratos Classes
///@{

template <unsigned int TDim, unsigned int TNumNodes>
class QSVMSAdjointElementData
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
    using VelocityDerivativeContributions = typename TResidualsDerivatives::template VariableDerivatives<typename QSVMSDerivativeUtilities<TDim>::template VelocityDerivative<TNumNodes, TDirectionIndex>>;

    template<unsigned int TDirectionIndex>
    using ShapeDerivatives = typename TResidualsDerivatives::template VariableDerivatives<typename QSVMSDerivativeUtilities<TDim>::template ShapeDerivative<TNumNodes, TDirectionIndex>>;

    template<unsigned int TDirectionIndex>
    using AccelerationDerivativeContributions = typename TResidualsDerivatives::template SecondDerivatives<TDirectionIndex>;

    static constexpr IndexType ElementDataContainerIndex = 0;

    static constexpr IndexType ResidualColumnOffset = 0;

    ///@}

public:
    ///@name Template Type Definitions
    ///@{

    static GeometryData::IntegrationMethod GetIntegrationMethod();

    using Residual = CalculationDataContainers<
                        std::tuple<
                            Data>,
                        std::tuple<
                            SubAssembly<ResidualsContributions, ElementDataContainerIndex, 0, ResidualColumnOffset>>
                        >;

    using ResidualStateVariableFirstDerivatives = CalculationDataContainers<
                                                    std::tuple<
                                                        Data>,
                                                    std::tuple<
                                                        SubAssembly<VelocityDerivativeContributions<0>, ElementDataContainerIndex, 0, ResidualColumnOffset>,
                                                        SubAssembly<VelocityDerivativeContributions<1>, ElementDataContainerIndex, 1, ResidualColumnOffset>,
                                                        SubAssembly<PressureDerivativeContributions   , ElementDataContainerIndex, 2, ResidualColumnOffset>>
                                                    >;

    using ResidualStateVariableSecondDerivatives = CalculationDataContainers<
                                                    std::tuple<
                                                        Data>,
                                                    std::tuple<
                                                        SubAssembly<AccelerationDerivativeContributions<0>, ElementDataContainerIndex, 0, ResidualColumnOffset>,
                                                        SubAssembly<AccelerationDerivativeContributions<1>, ElementDataContainerIndex, 1, ResidualColumnOffset>,
                                                        SubAssembly<ZeroDerivatives<3, 3>                 , ElementDataContainerIndex, 2, ResidualColumnOffset>>
                                                    >;

    using ResidualShapeDerivatives = CalculationDataContainers<
                                        std::tuple<
                                            Data>,
                                        std::tuple<
                                            SubAssembly<ShapeDerivatives<0>, ElementDataContainerIndex, 0, ResidualColumnOffset>,
                                            SubAssembly<ShapeDerivatives<1>, ElementDataContainerIndex, 1, ResidualColumnOffset>>
                                        >;

    ///@}
    ///@name Static Operations
    ///@{

    static void Check(
        const Element& rElement,
        const ProcessInfo& rProcessInfo);

    static std::vector<const Variable<double>*> GetDofVariablesList();

    ///@}
};
} // namespace Kratos
#endif // KRATOS_QS_VMS_ADJOINT_ELEMENT_DATA_H