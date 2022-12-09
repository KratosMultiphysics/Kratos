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

#pragma once

// System includes

// External includes

// Project includes

// Application includes
#include "custom_elements/data_containers/fluid_adjoint_derivatives.h"
#include "custom_elements/data_containers/qs_vms/qs_vms_derivative_utilities.h"
#include "custom_elements/data_containers/qs_vms/qs_vms_residual_derivatives.h"

namespace Kratos
{
///@name Kratos Classes
///@{

/**
 * @brief Traits class to hold QSVMS adjoint element data
 *
 * This class provides following type information about element data
 *      EquationAuxiliaries: Needs to provide following auxiliary static methods
 *                              static void Check(Element&, const ProcessInfo&) method -> Checks whether required data is there
 *                              static GeometryData::IntegrationMethod GetIntegrationMethod() -> Provides the integration method
 *
 *      Residual: Provides CalculationContainerTraits with the data containers and calculation containers
 *                to compute the residuals.
 *
 *      ResidualStateVariableFirstDerivatives: Provides CalculationContainerTraits with data containers and calculation
 *                                             containers to compute the residual state first derivatives.
 *
 *      ResidualStateVariableSecondDerivatives: Provides CalculationContainerTraits with data containers and calculation
 *                                              containers to compute the residual state second derivatives.
 *
 *      ResidualShapeDerivatives: Provides CalculationContainerTraits with data containers and calculation
 *                                containers to compute the residual shape derivatives.
 *
 * @tparam TDim             Dimensionality of the element.
 * @tparam TNumNodes        Number of nodes in the element.
 */
template <unsigned int TDim, unsigned int TNumNodes>
class QSVMSAdjointElementData
{
private:
    ///@name Private Type Definitions
    ///@{

    using IndexType = std::size_t;

    using TResidualsDerivatives = QSVMSResidualDerivatives<TDim, TNumNodes>;

    using Data = typename TResidualsDerivatives::QSVMSResidualData;

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

    /**
     * @brief This holds static helper functions such as Check and GetIntegrationMethod for the equations.
     */
    using EquationAuxiliaries = TResidualsDerivatives;

    /**
     * @brief This holds the container traits for primal residual computation which is required for analytical sensitivity computation.
     */
    using Residual = CalculationContainerTraits<
                        std::tuple<
                            Data>,
                        std::tuple<
                            SubAssembly<ResidualsContributions, ElementDataContainerIndex, 0, ResidualColumnOffset>>
                        >;
    /**
     * @brief This holds the container traits for analytical first derivative computations
     */
    using ResidualStateVariableFirstDerivatives = std::conditional_t<
                                                        TDim == 2,
                                                        CalculationContainerTraits<
                                                            std::tuple<
                                                                Data>,
                                                            std::tuple<
                                                                SubAssembly<VelocityDerivativeContributions<0>, ElementDataContainerIndex, 0, ResidualColumnOffset>,
                                                                SubAssembly<VelocityDerivativeContributions<1>, ElementDataContainerIndex, 1, ResidualColumnOffset>,
                                                                SubAssembly<PressureDerivativeContributions,    ElementDataContainerIndex, 2, ResidualColumnOffset>>
                                                            >,
                                                        CalculationContainerTraits<
                                                            std::tuple<
                                                                Data>,
                                                            std::tuple<
                                                                SubAssembly<VelocityDerivativeContributions<0>, ElementDataContainerIndex, 0, ResidualColumnOffset>,
                                                                SubAssembly<VelocityDerivativeContributions<1>, ElementDataContainerIndex, 1, ResidualColumnOffset>,
                                                                SubAssembly<VelocityDerivativeContributions<2>, ElementDataContainerIndex, 2, ResidualColumnOffset>,
                                                                SubAssembly<PressureDerivativeContributions,    ElementDataContainerIndex, 3, ResidualColumnOffset>>
                                                            >
                                                        >;
    /**
     * @brief This holds the traits for the analytical second derivative computations
     */
    using ResidualStateVariableSecondDerivatives = std::conditional_t<
                                                        TDim == 2,
                                                        CalculationContainerTraits<
                                                            std::tuple<
                                                                Data>,
                                                            std::tuple<
                                                                SubAssembly<AccelerationDerivativeContributions<0>, ElementDataContainerIndex, 0, ResidualColumnOffset>,
                                                                SubAssembly<AccelerationDerivativeContributions<1>, ElementDataContainerIndex, 1, ResidualColumnOffset>,
                                                                SubAssembly<ZeroDerivatives<TNumNodes, 3>,          ElementDataContainerIndex, 2, ResidualColumnOffset>>
                                                            >,
                                                        CalculationContainerTraits<
                                                            std::tuple<
                                                                Data>,
                                                            std::tuple<
                                                                SubAssembly<AccelerationDerivativeContributions<0>, ElementDataContainerIndex, 0, ResidualColumnOffset>,
                                                                SubAssembly<AccelerationDerivativeContributions<1>, ElementDataContainerIndex, 1, ResidualColumnOffset>,
                                                                SubAssembly<AccelerationDerivativeContributions<2>, ElementDataContainerIndex, 2, ResidualColumnOffset>,
                                                                SubAssembly<ZeroDerivatives<TNumNodes, 4>,          ElementDataContainerIndex, 3, ResidualColumnOffset>>
                                                            >
                                                        >;

    /**
     * @brief This holds the traits for the analytical shape derivative computations
     */
    using ResidualShapeDerivatives = std::conditional_t<
                                            TDim == 2,
                                            CalculationContainerTraits<
                                                std::tuple<
                                                    Data>,
                                                std::tuple<
                                                    SubAssembly<ShapeDerivatives<0>, ElementDataContainerIndex, 0, ResidualColumnOffset>,
                                                    SubAssembly<ShapeDerivatives<1>, ElementDataContainerIndex, 1, ResidualColumnOffset>>
                                                >,
                                            CalculationContainerTraits<
                                                std::tuple<
                                                    Data>,
                                                std::tuple<
                                                    SubAssembly<ShapeDerivatives<0>, ElementDataContainerIndex, 0, ResidualColumnOffset>,
                                                    SubAssembly<ShapeDerivatives<1>, ElementDataContainerIndex, 1, ResidualColumnOffset>,
                                                    SubAssembly<ShapeDerivatives<2>, ElementDataContainerIndex, 2, ResidualColumnOffset>>
                                                >
                                            >;

    ///@}
};
} // namespace Kratos