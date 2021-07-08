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

#if !defined(KRATOS_K_EPSILON_QS_VMS_RFC_ADJOINT_ELEMENT_DATA_H)
#define KRATOS_K_EPSILON_QS_VMS_RFC_ADJOINT_ELEMENT_DATA_H

// System includes

// External includes

// Project includes
#include "includes/node.h"
#include "includes/process_info.h"
#include "geometries/geometry.h"
#include "geometries/geometry_data.h"

// Application includes
#include "custom_elements/convection_diffusion_reaction_residual_based_flux_corrected_derivatives.h"
#include "custom_elements/data_containers/k_epsilon/epsilon_element_data_derivatives.h"
#include "custom_elements/data_containers/k_epsilon/k_element_data_derivatives.h"
#include "custom_elements/data_containers/qs_vms/qs_vms_derivative_utilities.h"
#include "custom_elements/data_containers/qs_vms/qs_vms_residual_derivatives.h"
#include "custom_elements/data_containers/qs_vms/rans_qs_vms_derivative_utilities.h"
#include "custom_elements/data_containers/k_epsilon/k_element_data.h"
#include "custom_elements/data_containers/k_epsilon/epsilon_element_data.h"

namespace Kratos
{
///@name Kratos Classes
///@{

namespace KEpsilonElementData
{

template <unsigned int TDim, unsigned int TNumNodes>
class QSVMSRFCAdjointElementData
{
public:
    ///@name Type Definitions
    ///@{

    using IndexType = std::size_t;

    using NodeType = Node<3>;

    using GeometryType = Geometry<NodeType>;

    ///@}
    ///@name Static Operations
    ///@{

    static void Check(
        const Element& rElement,
        const ProcessInfo& rProcessInfo);

    static std::vector<const Variable<double>*> GetDofVariablesList();

    static std::vector<const Variable<double>*> GetPrimalSecondDerivativeVariablesList();

    ///@}
    ///@name Classes
    ///@{

    class Fluid
    {
    public:
        ///@name Type Definitions
        ///@{

        using TResidualsDerivatives = QSVMSResidualDerivatives<TDim, TNumNodes>;

        ///@}
        ///@name Classes
        ///@{

        class Primal
        {
        public:
            ///@name Type Definitions
            ///@{

            using Data = typename TResidualsDerivatives::Data;

            using ResidualsContributions = typename TResidualsDerivatives::ResidualsContributions;

            ///@}
        };

        class StateDerivatives
        {
        public:
            ///@name Type Definitions
            ///@{

            class SecondDerivatives
            {
            public:
                ///@name Type Definitions
                ///@{

                using Data = typename TResidualsDerivatives::Data;

                using Acceleration = typename TResidualsDerivatives::SecondDerivatives;

                ///@}

            };

            ///@}
            ///@name Classes
            ///@{

            class FirstDerivatives
            {
            public:
                ///@name Type Definitions
                ///@{

                using Data = typename TResidualsDerivatives::Data;

                using Velocity = typename TResidualsDerivatives::template VariableDerivatives<typename QSVMSDerivativeUtilities<TDim>::template VelocityDerivative<TNumNodes>>;

                using Pressure = typename TResidualsDerivatives::template VariableDerivatives<typename QSVMSDerivativeUtilities<TDim>::template PressureDerivative<TNumNodes>>;

                using TurbulenceModelVariable1 = typename TResidualsDerivatives::template VariableDerivatives<typename RansQSVMSDerivativeUtilities<TDim>::template TurbulenceVariableDerivative<TNumNodes, KEpsilonElementData::KElementData<TDim>>>;

                using TurbulenceModelVariable2 = typename TResidualsDerivatives::template VariableDerivatives<typename RansQSVMSDerivativeUtilities<TDim>::template TurbulenceVariableDerivative<TNumNodes, KEpsilonElementData::EpsilonElementData<TDim>>>;

                ///@}
            };

            ///@}
        };

        class SensitivityDerivatives
        {
        public:
            ///@name Type Definitions
            ///@{

            using Data = typename TResidualsDerivatives::Data;

            using Shape = typename TResidualsDerivatives::template VariableDerivatives<typename QSVMSDerivativeUtilities<TDim>::template ShapeDerivative<TNumNodes>>;

            ///@}
        };

    ///@}
    };

    class TurbulenceModelEquation1
    {
    public:
        ///@name Classes
        ///@{

        using EquationDataType = KEpsilonElementData::KElementDataDerivatives<TDim, TNumNodes>;

        using TResidualsDerivatives = ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<TDim, TNumNodes, typename EquationDataType::Data>;

        class Primal
        {
        public:
            ///@name Type Definitions
            ///@{

            using Data = typename TResidualsDerivatives::Data;

            using ResidualsContributions = typename TResidualsDerivatives::ResidualsContributions;

            ///@}
        };

        class StateDerivatives
        {
        public:
            ///@name Classes
            ///@{

            class FirstDerivatives
            {
            public:
                ///@name Type Definitions
                ///@{

                using Data = typename TResidualsDerivatives::Data;

                using Velocity = typename TResidualsDerivatives::template VariableDerivatives<typename EquationDataType::UDerivative>;

                using TurbulenceModelVariable1 = typename TResidualsDerivatives::template VariableDerivatives<typename EquationDataType::KDerivative>;

                using TurbulenceModelVariable2 = typename TResidualsDerivatives::template VariableDerivatives<typename EquationDataType::EpsilonDerivative>;

                ///@}
            };

            class SecondDerivatives
            {
            public:
                ///@name Type Definitions
                ///@{

                using Data = typename TResidualsDerivatives::Data;

                using TurbulenceModelVariableRate1 = typename TResidualsDerivatives::SecondDerivatives;

                ///@}

            };

            ///@}
        };

        class SensitivityDerivatives
        {
        public:
            ///@name Type Definitions
            ///@{

            using Data = typename TResidualsDerivatives::Data;

            using Shape = typename TResidualsDerivatives::template VariableDerivatives<typename EquationDataType::ShapeDerivative>;

            ///@}
        };

        ///@}
    };

    class TurbulenceModelEquation2
    {
    public:
        ///@name Classes
        ///@{

        using EquationDataType = KEpsilonElementData::EpsilonElementDataDerivatives<TDim, TNumNodes>;

        using TResidualsDerivatives = ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<TDim, TNumNodes, typename EquationDataType::Data>;

        class Primal
        {
        public:
            ///@name Type Definitions
            ///@{

            using Data = typename TResidualsDerivatives::Data;

            using ResidualsContributions = typename TResidualsDerivatives::ResidualsContributions;

            ///@}
        };

        class StateDerivatives
        {
        public:
            ///@name Classes
            ///@{

            class FirstDerivatives
            {
            public:
                ///@name Type Definitions
                ///@{

                using Data = typename TResidualsDerivatives::Data;

                using Velocity = typename TResidualsDerivatives::template VariableDerivatives<typename EquationDataType::UDerivative>;

                using TurbulenceModelVariable1 = typename TResidualsDerivatives::template VariableDerivatives<typename EquationDataType::KDerivative>;

                using TurbulenceModelVariable2 = typename TResidualsDerivatives::template VariableDerivatives<typename EquationDataType::EpsilonDerivative>;

                ///@}
            };

            class SecondDerivatives
            {
            public:
                ///@name Type Definitions
                ///@{

                using Data = typename TResidualsDerivatives::Data;

                using TurbulenceModelVariableRate2 = typename TResidualsDerivatives::SecondDerivatives;

                ///@}
            };

            ///@}
        };

        class SensitivityDerivatives
        {
        public:
            ///@name Type Definitions
            ///@{

            using Data = typename TResidualsDerivatives::Data;

            using Shape = typename TResidualsDerivatives::template VariableDerivatives<typename EquationDataType::ShapeDerivative>;

            ///@}
        };

        ///@}
    };

    ///@}
};
} // namespace KEpsilonElementData
} // namespace Kratos

#endif // KRATOS_K_EPSILON_QS_VMS_RFC_ADJOINT_ELEMENT_DATA_H