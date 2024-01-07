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

#if !defined(KRATOS_STABILIZATION_VALIDATION_DIFFUSION_RFC_ADJOINT_ELEMENT_DATA_H)
#define KRATOS_STABILIZATION_VALIDATION_DIFFUSION_RFC_ADJOINT_ELEMENT_DATA_H

// System includes

// External includes

// Project includes
#include "includes/node.h"
#include "includes/process_info.h"
#include "geometries/geometry.h"
#include "geometries/geometry_data.h"

// Application includes
#include "custom_elements/convection_diffusion_reaction_residual_based_flux_corrected_derivatives.h"
#include "custom_elements/data_containers/stabilization_validation/diffusion_element_data.h"
#include "custom_elements/data_containers/stabilization_validation/diffusion_element_data_derivatives.h"

namespace Kratos
{
///@name Kratos Classes
///@{

namespace StabilizationValidationElementData
{

class DiffusionRFCAdjointElementData
{
public:
    ///@name Type Definitions
    ///@{

    using IndexType = std::size_t;

    using NodeType = Node<3>;

    using GeometryType = Geometry<NodeType>;

    static constexpr IndexType TDim = 2;

    static constexpr IndexType TNumNodes = 3;

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

     class ScalarEquation
    {
    public:
        ///@name Classes
        ///@{

        using EquationDataType = StabilizationValidationElementData::DiffusionElementDataDerivatives;

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

                using ScalarVariable = typename TResidualsDerivatives::template VariableDerivatives<typename EquationDataType::PhiDerivative>;

                ///@}
            };

            class SecondDerivatives
            {
            public:
                ///@name Type Definitions
                ///@{

                using Data = typename TResidualsDerivatives::Data;

                using ScalarVariableRate = typename TResidualsDerivatives::SecondDerivatives;

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
} // namespace StabilizationValidationElementData
} // namespace Kratos

#endif // KRATOS_STABILIZATION_VALIDATION_DIFFUSION_RFC_ADJOINT_ELEMENT_DATA_H