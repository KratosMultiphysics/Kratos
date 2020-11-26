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
#include "custom_elements/data_containers/qs_vms_adjoint_data/qs_vms_derivative_utilities.h"
#include "custom_elements/data_containers/qs_vms_adjoint_data/qs_vms_residual_derivatives.h"

namespace Kratos
{
///@name Kratos Classes
///@{

template <unsigned int TDim, unsigned int TNumNodes>
class QSVMSAdjointElementData
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
        const ProcessInfo& rProcessInfo)
    {
        // return TResidualDerivatives::Check(rGeometry, rProcessInfo);
        return 0;
    }

    static GeometryData::IntegrationMethod GetIntegrationMethod()
    {
        return TResidualDerivatives::GetIntegrationMethod();
    }

    ///@}
    ///@name Classes
    ///@{

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

            using Data = typename TResidualDerivatives::Data;

            using Velocity = typename TResidualDerivatives::VariableDerivatives<typename QSVMSDerivativeUtilities<TDim>::VelocityDerivative<TNumNodes>, 0>;

            using Pressure = typename TResidualDerivatives::VariableDerivatives<typename QSVMSDerivativeUtilities<TDim>::PressureDerivative<TNumNodes>, 0>;

            ///@}
        };

        class SecondDerivatives
        {
        public:
            ///@name Type Definitions
            ///@{

            // using Acceleration = nullptr;

            ///@}
        };

        ///@}
    };

    class SensitivityDerivatives
    {
    public:
        ///@name Classes
        ///@{

        class ShapeSensitivities
        {
        public:
            ///@name Type Definitions
            ///@{

            constexpr static IndexType TDerivativesDimension = TDim;

            using Data = typename TResidualDerivatives::Data;

            using Sensitivity = typename TResidualDerivatives::VariableDerivatives<typename QSVMSDerivativeUtilities<TDim>::ShapeDerivative<TNumNodes>, 0>;

            ///@}
        };

        ///@}
    };

    ///@}
};
} // namespace Kratos
#endif // KRATOS_QS_VMS_ADJOINT_ELEMENT_DATA_H