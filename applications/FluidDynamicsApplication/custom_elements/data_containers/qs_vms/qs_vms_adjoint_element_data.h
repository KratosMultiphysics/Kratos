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
#include "includes/node.h"
#include "includes/process_info.h"
#include "geometries/geometry.h"
#include "geometries/geometry_data.h"

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
public:
    ///@name Local Type Definitions
    ///@{

    using IndexType = std::size_t;

    using NodeType = Node<3>;

    using GeometryType = Geometry<NodeType>;

    using TResidualsDerivatives = QSVMSResidualDerivatives<TDim, TNumNodes>;

    ///@}
    ///@name Template Type Definitions
    ///@{

    using Data = typename TResidualsDerivatives::Data;

    using ResidualsContributions = typename TResidualsDerivatives::ResidualsContributions;

    template<unsigned int TDirectionIndex>
    using VelocityDerivativeContributions = typename TResidualsDerivatives::template VariableDerivatives<typename QSVMSDerivativeUtilities<TDim>::template VelocityDerivative<TNumNodes>, TDirectionIndex>;

    // using VelocityDerivativeContributions = typename TResidualsDerivatives::template VariableDerivatives<typename QSVMSDerivativeUtilities<TDim>::template VelocityDerivative<TNumNodes>>;

    using PressureDerivativeContributions = typename TResidualsDerivatives::template VariableDerivatives<typename QSVMSDerivativeUtilities<TDim>::template PressureDerivative<TNumNodes>>;

    template<unsigned int TDirectionIndex>
    using AccelerationDerivativeContributions = typename TResidualsDerivatives::template SecondDerivatives<TDirectionIndex>;

    using ShapeDerivatives = typename TResidualsDerivatives::template VariableDerivatives<typename QSVMSDerivativeUtilities<TDim>::template ShapeDerivative<TNumNodes>>;


    using ResidualStateVariableFirstDerivatives = FirstDerivatives<
                                                    Data,
                                                    SubAssembly<0, 0, VelocityDerivativeContributions<0>>,
                                                    SubAssembly<1, 0, VelocityDerivativeContributions<1>>,
                                                    SubAssembly<2, 0, PressureDerivativeContributions>
                                                    >;

    using ResidualStateVariableSecondDerivatives = SecondDerivatives<
                                                    Data,
                                                    SubAssembly<0, 0, AccelerationDerivativeContributions<0>>,
                                                    SubAssembly<1, 0, AccelerationDerivativeContributions<1>>,
                                                    SubAssembly<2, 0, ZeroDerivatives>
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