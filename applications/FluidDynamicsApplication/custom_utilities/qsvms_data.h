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
#include "custom_utilities/fluid_element_data.h"

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

///@}
///@name Public Members
///@{

NodalVectorData Velocity;
NodalVectorData MeshVelocity;
NodalVectorData BodyForce;
NodalVectorData MomentumProjection;

NodalScalarData Pressure;
NodalScalarData Density;
NodalScalarData DynamicViscosity;
NodalScalarData MassProjection;

double CSmagorinsky;
double DeltaTime;
double DynamicTau;
int UseOSS;

///@}
///@name Public Operations
///@{

void Initialize(const Element& rElement, const ProcessInfo& rProcessInfo) override
{
    const Geometry< Node<3> >& r_geometry = rElement.GetGeometry();
    this->FillFromNodalData(Velocity,VELOCITY,r_geometry);
    this->FillFromNodalData(MeshVelocity,MESH_VELOCITY,r_geometry);
    this->FillFromNodalData(BodyForce,BODY_FORCE,r_geometry);
    this->FillFromNodalData(MomentumProjection,ADVPROJ,r_geometry);
    this->FillFromNodalData(Pressure,PRESSURE,r_geometry);
    this->FillFromNodalData(Density,DENSITY,r_geometry);
    this->FillFromNodalData(DynamicViscosity,DYNAMIC_VISCOSITY,r_geometry);
    this->FillFromNodalData(MassProjection,DIVPROJ,r_geometry);
    this->FillFromElementData(CSmagorinsky,C_SMAGORINSKY,rElement);
    this->FillFromProcessInfo(DeltaTime,DELTA_TIME,rProcessInfo);
    this->FillFromProcessInfo(DynamicTau,DYNAMIC_TAU,rProcessInfo);
    this->FillFromProcessInfo(UseOSS,OSS_SWITCH,rProcessInfo);
}

static int Check(const Element& rElement, const ProcessInfo& rProcessInfo)
{
    const Geometry< Node<3> >& r_geometry = rElement.GetGeometry();
    
    KRATOS_CHECK_VARIABLE_KEY(VELOCITY);
    KRATOS_CHECK_VARIABLE_KEY(MESH_VELOCITY);
    KRATOS_CHECK_VARIABLE_KEY(BODY_FORCE);
    KRATOS_CHECK_VARIABLE_KEY(ADVPROJ);
    KRATOS_CHECK_VARIABLE_KEY(PRESSURE);
    KRATOS_CHECK_VARIABLE_KEY(DENSITY);
    KRATOS_CHECK_VARIABLE_KEY(DYNAMIC_VISCOSITY);
    KRATOS_CHECK_VARIABLE_KEY(DIVPROJ);

    for (unsigned int i = 0; i < TNumNodes; i++)
    {
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VELOCITY,r_geometry[i]);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(MESH_VELOCITY,r_geometry[i]);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(BODY_FORCE,r_geometry[i]);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(ADVPROJ,r_geometry[i]);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(PRESSURE,r_geometry[i]);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DENSITY,r_geometry[i]);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DYNAMIC_VISCOSITY,r_geometry[i]);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DIVPROJ,r_geometry[i]);
    }

    KRATOS_CHECK_VARIABLE_KEY(C_SMAGORINSKY);
    KRATOS_CHECK_VARIABLE_KEY(DELTA_TIME);
    KRATOS_CHECK_VARIABLE_KEY(DYNAMIC_TAU);
    KRATOS_CHECK_VARIABLE_KEY(OSS_SWITCH);

    return 0;
}

///@}

};

///@}

///@}

}

#endif