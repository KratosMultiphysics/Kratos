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

    using TResidualDerivatives = QSVMSResidualDerivatives<TDim, TNumNodes>;

    static constexpr IndexType TBlockSize = TResidualDerivatives::TBlockSize;

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

            using U = typename TResidualDerivatives::template VariableDerivatives<typename QSVMSDerivativeUtilities<TDim>::template VelocityDerivative<TNumNodes>, 0>;

            using P = typename TResidualDerivatives::template VariableDerivatives<typename QSVMSDerivativeUtilities<TDim>::template PressureDerivative<TNumNodes>, 0>;

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
} // namespace KEpsilonElementData
} // namespace Kratos

#endif // KRATOS_K_EPSILON_QS_VMS_RFC_ADJOINT_ELEMENT_DATA_H