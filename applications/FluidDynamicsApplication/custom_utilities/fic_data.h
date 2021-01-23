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


#if !defined(KRATOS_FIC_DATA_H)
#define KRATOS_FIC_DATA_H

#include "fluid_dynamics_application_variables.h"
#include "custom_utilities/fluid_element_data.h"

namespace Kratos {

///@addtogroup FluidDynamicsApplication
///@{

///@name Kratos classes
///@{

template< size_t TDim, size_t TNumNodes,bool TElementIntegratesInTime = false >
class FICData : public FluidElementData<TDim,TNumNodes, TElementIntegratesInTime>
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

NodalScalarData Pressure;

double Density;
double DeltaTime;      // Time increment
double FICBeta;        // FIC Beta parameter
double DynamicTau;     // Time stabilization coefficient for incompressibility stabilization

/// Auxiliary container for the local matrix at the integration point (stored to save reallocation at each point)
BoundedMatrix<double,TNumNodes*(TDim+1),TNumNodes*(TDim+1)> LHS;

///@}
///@name Public Operations
///@{

void Initialize(const Element& rElement, const ProcessInfo& rProcessInfo) override
{
    // Base class Initialize manages constitutive law parameters
    FluidElementData<TDim,TNumNodes, TElementIntegratesInTime>::Initialize(rElement,rProcessInfo);

    const Geometry< Node<3> >& r_geometry = rElement.GetGeometry();
    const Properties& r_properties = rElement.GetProperties();
    this->FillFromHistoricalNodalData(Velocity,VELOCITY,r_geometry);
    this->FillFromHistoricalNodalData(MeshVelocity,MESH_VELOCITY,r_geometry);
    this->FillFromHistoricalNodalData(BodyForce,BODY_FORCE,r_geometry);
    this->FillFromHistoricalNodalData(Pressure,PRESSURE,r_geometry);
    this->FillFromProperties(Density,DENSITY,r_properties);
    this->FillFromProcessInfo(DeltaTime,DELTA_TIME,rProcessInfo);
    this->FillFromProcessInfo(FICBeta,FIC_BETA,rProcessInfo);
    this->FillFromProcessInfo(DynamicTau,DYNAMIC_TAU,rProcessInfo);
}

static int Check(const Element& rElement, const ProcessInfo& rProcessInfo)
{
    const Geometry< Node<3> >& r_geometry = rElement.GetGeometry();

    for (unsigned int i = 0; i < TNumNodes; i++)
    {
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VELOCITY,r_geometry[i]);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(MESH_VELOCITY,r_geometry[i]);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(BODY_FORCE,r_geometry[i]);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(PRESSURE,r_geometry[i]);
    }

    return 0;
}

///@}

};

///@}

///@}

}

#endif