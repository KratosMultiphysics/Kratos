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

#if !defined(KRATOS_TIME_INTEGRATED_FIC_DATA_H)
#define KRATOS_TIME_INTEGRATED_FIC_DATA_H

#include "fluid_dynamics_application_variables.h"
#include "custom_utilities/fic_data.h"

namespace Kratos {

///@addtogroup FluidDynamicsApplication
///@{

///@name Kratos classes
///@{

template< size_t TDim, size_t TNumNodes >
class TimeIntegratedFICData : public FICData<TDim,TNumNodes,true>
{
public:

///@name Type Definitions
///@{

using NodalScalarData = typename FluidElementData<TDim,TNumNodes, false>::NodalScalarData;
using NodalVectorData = typename FluidElementData<TDim,TNumNodes, false>::NodalVectorData;

///@}
///@name Public Members
///@{

NodalVectorData Velocity_OldStep1;
NodalVectorData Velocity_OldStep2;

double bdf0;
double bdf1;
double bdf2;

///@}
///@name Public Operations
///@{

void Initialize(const Element& rElement, const ProcessInfo& rProcessInfo) override
{
    FICData<TDim,TNumNodes,true>::Initialize(rElement,rProcessInfo);
    const Geometry< Node<3> >& r_geometry = rElement.GetGeometry();
    this->FillFromHistoricalNodalData(Velocity_OldStep1,VELOCITY,r_geometry,1);
    this->FillFromHistoricalNodalData(Velocity_OldStep2,VELOCITY,r_geometry,2);

    const Vector& BDFVector = rProcessInfo[BDF_COEFFICIENTS];
    bdf0 = BDFVector[0];
    bdf1 = BDFVector[1];
    bdf2 = BDFVector[2];
}


static int Check(const Element& rElement, const ProcessInfo& rProcessInfo)
{
    FICData<TDim,TNumNodes,true>::Check(rElement,rProcessInfo);
    KRATOS_CHECK_VARIABLE_KEY(BDF_COEFFICIENTS);

    return 0;
}

///@}

};

///@}

///@}

}

#endif