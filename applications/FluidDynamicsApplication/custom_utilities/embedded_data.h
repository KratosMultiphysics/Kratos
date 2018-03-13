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

#if !defined(KRATOS_EMBEDDED_DATA_H)
#define KRATOS_EMBEDDED_DATA_H

#include "fluid_dynamics_application_variables.h"
#include "custom_utilities/fluid_element_data.h"

namespace Kratos {

///@addtogroup FluidDynamicsApplication
///@{

///@name Kratos classes
///@{

template< class TFluidData >
class EmbeddedData : public TFluidData
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

NodalScalarData Distance;

Matrix PositiveSideN;
ShapeFunctionsGradientsType PositiveSideDNDX;
Vector PositiveSideWeights;

Matrix PositiveInterfaceN;
ShapeFunctionsGradientsType PositiveInterfaceDNDX;
Vector PositiveInterfaceWeights;
InterfaceNormalsType PositiveInterfaceUnitNormals;

std::vector< size_t > PositiveIndices;
std::vector< size_t > NegativeIndices;

size_t NumPositiveNodes;
size_t NumNegativeNodes;

///@}
///@name Public Operations
///@{

void Initialize(
    const Element& rElement, const ProcessInfo& rProcessInfo) override {
    TFluidData::Initialize(rElement, rProcessInfo);
    const Geometry<Node<3> >& r_geometry = rElement.GetGeometry();
    this->FillFromNodalData(Distance, DISTANCE, r_geometry);

    NumPositiveNodes = 0;
    NumNegativeNodes = 0;
}

static int Check(const Element& rElement, const ProcessInfo& rProcessInfo)
{
    const Geometry< Node<3> >& r_geometry = rElement.GetGeometry();
    
    KRATOS_CHECK_VARIABLE_KEY(DISTANCE);

    for (unsigned int i = 0; i < TFluidData::NumNodes; i++) {
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DISTANCE,r_geometry[i]);
    }

    int out = TFluidData::Check(rElement,rProcessInfo);
    return out;
}

bool IsCut() {
    return (NumPositiveNodes > 0) && (NumNegativeNodes > 0);
}

///@}

};

///@}

///@}

}

#endif