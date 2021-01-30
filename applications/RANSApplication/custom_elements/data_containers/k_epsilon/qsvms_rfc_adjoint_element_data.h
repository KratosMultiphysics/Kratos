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
#include "custom_elements/data_containers/qs_vms/qs_vms_derivative_utilities.h"
#include "custom_elements/data_containers/qs_vms/qs_vms_residual_derivatives.h"
#include "custom_elements/data_containers/k_epsilon/k_element_data_derivatives.h"
#include "custom_elements/data_containers/k_epsilon/epsilon_element_data_derivatives.h"
#include "custom_elements/convection_diffusion_reaction_residual_based_flux_corrected_derivatives.h"

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

    static int Check(
        const GeometryType& rGeometry,
        const ProcessInfo& rProcessInfo);

    static GeometryData::IntegrationMethod GetIntegrationMethod();

    static std::vector<const Variable<double>*> GetDofVariablesList();

    ///@}
    ///@name Classes
    ///@{

    class Fluid
    {
    public:
        ///@name Type Definitions
        ///@{

        using TResidualDerivatives = QSVMSResidualDerivatives<TDim, TNumNodes>;

        ///@}
        ///@name Classes
        ///@{

        class Primal
        {
        public:
            ///@name Type Definitions
            ///@{

            using Data = typename TResidualDerivatives::Data;

            using ResidualContributions = typename TResidualDerivatives::ResidualContributions;

            ///@}
        };

        class StateDerivatives
        {
        public:
            ///@name Type Definitions
            ///@{

            using SecondDerivatives = typename TResidualDerivatives::SecondDerivatives;

            ///@}
            ///@name Classes
            ///@{

            class FirstDerivatives
            {
            public:
                ///@name Type Definitions
                ///@{

                using Data = typename TResidualDerivatives::Data;

                using Velocity = typename TResidualDerivatives::template VariableDerivatives<typename QSVMSDerivativeUtilities<TDim>::template VelocityDerivative<TNumNodes>, 0>;

                using Pressure = typename TResidualDerivatives::template VariableDerivatives<typename QSVMSDerivativeUtilities<TDim>::template PressureDerivative<TNumNodes>, 0>;

                using TurbulenceModelVariable1 = typename TResidualDerivatives::template VariableDerivatives<typename QSVMSDerivativeUtilities<TDim>::template PressureDerivative<TNumNodes>, 0>;

                using TurbulenceModelVariable2 = typename TResidualDerivatives::template VariableDerivatives<typename QSVMSDerivativeUtilities<TDim>::template PressureDerivative<TNumNodes>, 0>;

                ///@}
            };

            ///@}
        };

        class SensitivityDerivatives
        {
        public:
            ///@name Type Definitions
            ///@{

            using Data = typename TResidualDerivatives::Data;

            using Shape = typename TResidualDerivatives::template VariableDerivatives<typename QSVMSDerivativeUtilities<TDim>::template ShapeDerivative<TNumNodes>, 0>;

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
            ///@name Type Definitions
            ///@{

            using SecondDerivatives = typename TResidualsDerivatives::SecondDerivatives;

            ///@}
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
            ///@name Type Definitions
            ///@{

            using SecondDerivatives = typename TResidualsDerivatives::SecondDerivatives;

            ///@}
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