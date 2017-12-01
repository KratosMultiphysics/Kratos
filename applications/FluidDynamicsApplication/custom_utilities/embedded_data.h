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

#include "modified_shape_functions/modified_shape_functions.h"
#include "modified_shape_functions/triangle_2d_3_modified_shape_functions.h"
#include "modified_shape_functions/tetrahedra_3d_4_modified_shape_functions.h"

namespace Kratos {

///@addtogroup FluidDynamicsApplication
///@{

///@name Kratos classes
///@{

namespace Internals {

template <size_t TDim, size_t TNumNodes>
ModifiedShapeFunctions::Pointer GetShapeFunctionCalculator(
    const Element& rElement, const Vector& rDistance);

template <>
ModifiedShapeFunctions::Pointer GetShapeFunctionCalculator<2, 3>(
    const Element& rElement, const Vector& rDistance) {
    return ModifiedShapeFunctions::Pointer(new Triangle2D3ModifiedShapeFunctions(rElement.pGetGeometry(),rDistance));
}

template <>
ModifiedShapeFunctions::Pointer GetShapeFunctionCalculator<3, 4>(
    const Element& rElement, const Vector& rDistance) {
    return ModifiedShapeFunctions::Pointer(new Tetrahedra3D4ModifiedShapeFunctions(rElement.pGetGeometry(),rDistance));
}
}

template< class TFluidData >
class EmbeddedData : public TFluidData
{
public:

///@name Type Definitions
///@{

using NodalScalarData = typename TFluidData::NodalScalarData;
using NodalVectorData = typename TFluidData::NodalVectorData;

///@}
///@name Public Members
///@{

NodalScalarData Distance;

Matrix PositiveSideN;
std::vector< Matrix > PositiveSideDNDX;
Vector PositiveSideWeights;

Matrix PositiveInterfaceN;
std::vector< Matrix > PositiveInterfaceDNDX;
Vector PositiveInterfaceWeights;
std::vector< Vector > PositiveInterfaceUnitNormals;

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

    // Auxiliary distance vector for the element subdivision utility
    Vector distances = this->Distance;

    ModifiedShapeFunctions::Pointer p_calculator =
        Internals::GetShapeFunctionCalculator<TFluidData::Dim,
            TFluidData::NumNodes>(rElement, distances);

    // Fluid side
    p_calculator->ComputePositiveSideShapeFunctionsAndGradientsValues(
        PositiveSideN, PositiveSideDNDX, PositiveSideWeights,
        GeometryData::GI_GAUSS_2);

    // Fluid side interface
    p_calculator->ComputeInterfacePositiveSideShapeFunctionsAndGradientsValues(
        PositiveInterfaceN, PositiveInterfaceDNDX, PositiveInterfaceWeights,
        GeometryData::GI_GAUSS_2);

    // Fluid side interface normals
    p_calculator->ComputePositiveSideInterfaceUnitNormals(
        PositiveInterfaceUnitNormals, GeometryData::GI_GAUSS_2);
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

///@}

};

///@}

///@}

}

#endif