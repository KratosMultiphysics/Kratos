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

#if !defined(KRATOS_EMBEDDED_DISCONTINUOUS_DATA_H)
#define KRATOS_EMBEDDED_DISCONTINUOUS_DATA_H

#include "fluid_dynamics_application_variables.h"
#include "custom_utilities/fluid_element_data.h"

namespace Kratos {

///@addtogroup FluidDynamicsApplication
///@{

///@name Kratos classes
///@{

template< class TFluidData >
class EmbeddedDiscontinuousData : public TFluidData
{
public:

///@name Type Definitions
///@{

using NodalScalarData = typename TFluidData::NodalScalarData;
using NodalVectorData = typename TFluidData::NodalVectorData;

typedef GeometryData::ShapeFunctionsGradientsType ShapeFunctionsGradientsType;
typedef std::vector< Vector > InterfaceNormalsType;

///@}
///@name Public Members
///@{

NodalScalarData ElementalDistances;

Matrix PositiveSideN;
Matrix NegativeSideN;
ShapeFunctionsGradientsType PositiveSideDNDX;
ShapeFunctionsGradientsType NegativeSideDNDX;
Vector PositiveSideWeights;
Vector NegativeSideWeights;

Matrix PositiveInterfaceN;
Matrix NegativeInterfaceN;
ShapeFunctionsGradientsType PositiveInterfaceDNDX;
ShapeFunctionsGradientsType NegativeInterfaceDNDX;
Vector PositiveInterfaceWeights;
Vector NegativeInterfaceWeights;
InterfaceNormalsType PositiveInterfaceUnitNormals;
InterfaceNormalsType NegativeInterfaceUnitNormals;

std::vector< size_t > PositiveIndices;
std::vector< size_t > NegativeIndices;

size_t NumPositiveNodes;
size_t NumNegativeNodes;

///@}
///@name Public Operations
///@{

void Initialize(
    const Element& rElement,
    const ProcessInfo& rProcessInfo) override
{
    TFluidData::Initialize(rElement, rProcessInfo);
    this->FillFromElementData(ElementalDistances, ELEMENTAL_DISTANCES, rElement);

    NumPositiveNodes = 0;
    NumNegativeNodes = 0;
}

static int Check(
    const Element& rElement,
    const ProcessInfo& rProcessInfo)
{
    KRATOS_CHECK_VARIABLE_KEY(ELEMENTAL_DISTANCES);

    int out = TFluidData::Check(rElement,rProcessInfo);
    return out;
}

bool IsCut()
{
    return (NumPositiveNodes > 0) && (NumNegativeNodes > 0);
}

///@}

};

///@}

///@}

}

#endif