//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//

#if !defined(KRATOS_SYMBOLIC_STOKES_DATA_H)
#define KRATOS_SYMBOLIC_STOKES_DATA_H

// External includes

// System includes

// Project includes
#include "includes/constitutive_law.h"

// Application includes
#include "fluid_dynamics_application_variables.h"
#include "custom_elements/data_containers/fluid_element_data.h"
#include "utilities/element_size_calculator.h"

namespace Kratos {

///@addtogroup FluidDynamicsApplication
///@{

///@name Kratos classes
///@{

template< size_t TDim, size_t TNumNodes >
class SymbolicStokesData : public FluidElementData<TDim,TNumNodes, true>
{
public:

///@name Type Definitions
///@{

using NodalScalarData = typename FluidElementData<TDim, TNumNodes, true>::NodalScalarData;
using NodalVectorData = typename FluidElementData<TDim, TNumNodes, true>::NodalVectorData;
using ShapeFunctionsType = typename FluidElementData<TDim, TNumNodes, true>::ShapeFunctionsType;
using MatrixRowType = typename FluidElementData<TDim, TNumNodes, true>::MatrixRowType;

static constexpr std::size_t BlockSize = TDim + 1;

///@}
///@name Public Members
///@{

NodalVectorData Velocity;
NodalVectorData Velocity_OldStep1;
NodalVectorData Velocity_OldStep2;
NodalVectorData BodyForce;

NodalScalarData Pressure;

double Density;
double DynamicViscosity;
double DeltaTime;      // Time increment
double DynamicTau;     // Dynamic tau considered in ASGS stabilization coefficients

double bdf0;
double bdf1;
double bdf2;

// Auxiliary containers for the symbolically-generated matrices
BoundedMatrix<double,TNumNodes*(TDim+1),TNumNodes*(TDim+1)> lhs;
array_1d<double,TNumNodes*(TDim+1)> rhs;

double ElementSize;

///@}
///@name Public Operations
///@{

void Initialize(
    const Element& rElement,
    const ProcessInfo& rProcessInfo) override
{
    // Base class Initialize manages constitutive law parameters
    FluidElementData<TDim, TNumNodes, true>::Initialize(rElement, rProcessInfo);

    const Geometry< Node >& r_geometry = rElement.GetGeometry();
    const Properties& r_properties = rElement.GetProperties();
    this->FillFromHistoricalNodalData(Velocity, VELOCITY, r_geometry);
    this->FillFromHistoricalNodalData(Velocity_OldStep1, VELOCITY, r_geometry, 1);
    this->FillFromHistoricalNodalData(Velocity_OldStep2, VELOCITY, r_geometry, 2);
    this->FillFromHistoricalNodalData(BodyForce, BODY_FORCE, r_geometry);
    this->FillFromHistoricalNodalData(Pressure, PRESSURE, r_geometry);
    this->FillFromProperties(Density, DENSITY, r_properties);
    this->FillFromProperties(DynamicViscosity, DYNAMIC_VISCOSITY, r_properties);
    this->FillFromProcessInfo(DeltaTime, DELTA_TIME, rProcessInfo);
    this->FillFromProcessInfo(DynamicTau, DYNAMIC_TAU, rProcessInfo);

    const Vector& BDFVector = rProcessInfo[BDF_COEFFICIENTS];
    bdf0 = BDFVector[0];
    bdf1 = BDFVector[1];
    bdf2 = BDFVector[2];

    ElementSize = ElementSizeCalculator<TDim,TNumNodes>::MinimumElementSize(r_geometry);

    noalias(lhs) = ZeroMatrix(TNumNodes*(TDim+1),TNumNodes*(TDim+1));
    noalias(rhs) = ZeroVector(TNumNodes*(TDim+1));
}

void UpdateGeometryValues(
    const unsigned int IntegrationPointIndex,
    double NewWeight,
    const MatrixRowType& rN,
    const BoundedMatrix<double, TNumNodes, TDim>& rDN_DX) override
{
    FluidElementData<TDim,TNumNodes, true>::UpdateGeometryValues(IntegrationPointIndex, NewWeight, rN, rDN_DX);
}

static int Check(const Element& rElement, const ProcessInfo& rProcessInfo)
{
    const Geometry< Node >& r_geometry = rElement.GetGeometry();

    for (unsigned int i = 0; i < TNumNodes; i++)
    {
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VELOCITY, r_geometry[i]);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(BODY_FORCE, r_geometry[i]);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(PRESSURE, r_geometry[i]);
    }

    return 0;
}

///@}

};

///@}

///@}

}

#endif // KRATOS_SYMBOLIC_STOKES_DATA_H
