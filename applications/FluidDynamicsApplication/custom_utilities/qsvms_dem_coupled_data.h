//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Joaquin Gonzalez-Usua
//

#if !defined(KRATOS_QSVMSDEMCOUPLED_DATA_H)
#define KRATOS_QSVMSDEMCOUPLED_DATA_H

// System includes

// External includes

// Project includes
#include "geometries/geometry.h"
#include "utilities/element_size_calculator.h"
#include "includes/checks.h"

// Application includes
#include "fluid_dynamics_application_variables.h"
#include "custom_utilities/fluid_element_data.h"
#include "custom_utilities/qsvms_data.h"

namespace Kratos {

///@addtogroup FluidDynamicsApplication
///@{

///@name Kratos classes
///@{

template< size_t TDim, size_t TNumNodes, bool TElementIntegratesInTime = false>
class QSVMSDEMCoupledData : public QSVMSData<TDim, TNumNodes, TElementIntegratesInTime>
{
public:

///@name Type Definitions
///@{


constexpr static unsigned int Dim = TDim;
constexpr static unsigned int NumNodes = TNumNodes;
constexpr static unsigned int BlockSize = Dim + 2;
constexpr static unsigned int DofSize = NumNodes * BlockSize;
constexpr static bool ElementIntegratesInTime = TElementIntegratesInTime;

using NodalScalarData = typename FluidElementData<TDim,TNumNodes, false>::NodalScalarData;
using NodalVectorData = typename FluidElementData<TDim,TNumNodes, false>::NodalVectorData;
using NodalTensorData = typename FluidElementData<TDim,TNumNodes, false>::NodalTensorData;
using ShapeFunctionsSecondDerivativesType = DenseVector<Matrix>;
using MatrixRowType = MatrixRow<Matrix>;
using ShapeDerivativesType = BoundedMatrix<double,TNumNodes,TDim>;

///@}
///@name Public Members
///@{

NodalScalarData FluidFraction;
NodalScalarData FluidFractionRate;
NodalScalarData MassSource;

NodalVectorData FluidFractionGradient;
NodalVectorData Acceleration;
NodalVectorData BodyForce;

NodalTensorData Permeability;

ShapeFunctionsSecondDerivativesType DDN_DDX;

double ElementSize;

///@}
///@name Public Operations
///@{

void Initialize(
    const Element& rElement,
    const ProcessInfo& rProcessInfo) override
{
    // Base class Initialize manages constitutive law parameters
    QSVMSData<TDim, TNumNodes, TElementIntegratesInTime>::Initialize(rElement,rProcessInfo);
    const auto& r_geometry = rElement.GetGeometry();
    this->FillFromHistoricalNodalData(FluidFraction, FLUID_FRACTION, r_geometry);
    this->FillFromHistoricalNodalData(FluidFractionRate, FLUID_FRACTION_RATE, r_geometry);
    this->FillFromHistoricalNodalData(FluidFractionGradient, FLUID_FRACTION_GRADIENT, r_geometry);
    this->FillFromHistoricalNodalData(Permeability, PERMEABILITY, r_geometry);
    this->FillFromHistoricalNodalData(MassSource, MASS_SOURCE, r_geometry);
    this->FillFromHistoricalNodalData(Acceleration, ACCELERATION, r_geometry);
    this->FillFromHistoricalNodalData(BodyForce,BODY_FORCE,r_geometry);

    ElementSize = ElementSizeCalculator<TDim,TNumNodes>::AverageElementSize(r_geometry);
}

void UpdateSecondDerivativesValues(const ShapeFunctionsSecondDerivativesType& rDDN_DDX)
{
    ShapeFunctionsSecondDerivativesType& ddn_ddx = this->DDN_DDX;
    ddn_ddx = rDDN_DDX;
}

///@}
};

///@}

///@}

}

#endif