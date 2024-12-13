//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela
//

#if !defined(KRATOS_QSVMS_DATA_H)
#define KRATOS_QSVMS_DATA_H

#include "fluid_dynamics_application_variables.h"
#include "data_containers/fluid_element_data.h"
#include "utilities/element_size_calculator.h"

namespace Kratos {

///@addtogroup FluidDynamicsApplication
///@{

///@name Kratos classes
///@{

template< size_t TDim, size_t TNumNodes,bool TElementIntegratesInTime = false >
class QSVMSData : public FluidElementData<TDim,TNumNodes, TElementIntegratesInTime>
{
public:

///@name Type Definitions
///@{

using NodalScalarData = typename FluidElementData<TDim,TNumNodes, false>::NodalScalarData;
using NodalVectorData = typename FluidElementData<TDim,TNumNodes, false>::NodalVectorData;

static constexpr std::size_t BlockSize = TDim + 1;

///@}
///@name Public Members
///@{

NodalVectorData Velocity;
NodalVectorData MeshVelocity;
NodalVectorData BodyForce;
NodalVectorData MomentumProjection;

NodalScalarData Pressure;
NodalScalarData MassProjection;

double Density;
double DynamicViscosity;
double CSmagorinsky;
double DeltaTime;
double DynamicTau;
int UseOSS;

double ElementSize;

/// Auxiliary container for the local matrix at the integration point (stored to save reallocation at each point)
BoundedMatrix<double,TNumNodes*(TDim+1),TNumNodes*(TDim+1)> LHS;

///@}
///@name Public Operations
///@{

void Initialize(const Element& rElement, const ProcessInfo& rProcessInfo) override
{
    // Base class Initialize manages constitutive law parameters
    FluidElementData<TDim,TNumNodes, TElementIntegratesInTime>::Initialize(rElement,rProcessInfo);

    const Geometry< Node >& r_geometry = rElement.GetGeometry();
    const Properties& r_properties = rElement.GetProperties();
    this->FillFromHistoricalNodalData(Velocity,VELOCITY,r_geometry);
    this->FillFromHistoricalNodalData(MeshVelocity,MESH_VELOCITY,r_geometry);
    this->FillFromHistoricalNodalData(BodyForce,BODY_FORCE,r_geometry);
    this->FillFromHistoricalNodalData(MomentumProjection,ADVPROJ,r_geometry);
    this->FillFromHistoricalNodalData(Pressure,PRESSURE,r_geometry);
    this->FillFromHistoricalNodalData(MassProjection,DIVPROJ,r_geometry);
    this->FillFromProperties(Density,DENSITY,r_properties);
    this->FillFromProperties(DynamicViscosity,DYNAMIC_VISCOSITY,r_properties); //TODO: remove once we have a Smagorinky constitutive law
    this->FillFromElementData(CSmagorinsky,C_SMAGORINSKY,rElement); //TODO: remove once we have a Smagorinky constitutive law
    this->FillFromProcessInfo(DeltaTime,DELTA_TIME,rProcessInfo);
    this->FillFromProcessInfo(DynamicTau,DYNAMIC_TAU,rProcessInfo);
    this->FillFromProcessInfo(UseOSS,OSS_SWITCH,rProcessInfo);

    ElementSize = ElementSizeCalculator<TDim,TNumNodes>::MinimumElementSize(r_geometry);
}

static int Check(const Element& rElement, const ProcessInfo& rProcessInfo)
{
    const Geometry< Node >& r_geometry = rElement.GetGeometry();

    for (unsigned int i = 0; i < TNumNodes; i++)
    {
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VELOCITY,r_geometry[i]);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(MESH_VELOCITY,r_geometry[i]);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(BODY_FORCE,r_geometry[i]);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(ADVPROJ,r_geometry[i]);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(PRESSURE,r_geometry[i]);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DIVPROJ,r_geometry[i]);
    }

    return 0;
}

///@}

};

///@}

///@}

}

#endif
