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

#include "fluid_element_data.h"

namespace Kratos
{
template <unsigned int TDim, unsigned int TNumNodes>
FluidElementData<TDim, TNumNodes>::FluidElementData(Geometry<Node<3>>& rGeom)
    : Pressure(TDim, 0.0), Density(TDim, 0.0), Viscosity(TDim, 0.0)
{
    for (unsigned int i = 0; i < TNumNodes; i++)
    {
        Node<3>& rNode = rGeom[i];
        const array_1d<double, 3>& rNodalVel = rNode.FastGetSolutionStepValue(VELOCITY);
        const array_1d<double, 3>& rNodalMeshVel =
            rNode.FastGetSolutionStepValue(MESH_VELOCITY);

        for (unsigned int d = 0; d < 3; d++)
        {
            Velocity(i, d) = rNodalVel[d];
            MeshVelocity(i, d) = rNodalMeshVel[d];
        }

        Pressure[i] = rNode.FastGetSolutionStepValue(PRESSURE);
        Density[i] = rNode.FastGetSolutionStepValue(DENSITY);
        Viscosity[i] = rNode.FastGetSolutionStepValue(VISCOSITY);
    }
}

template< unsigned int TDim, unsigned int TNumNodes >
constexpr unsigned int FluidElementData<TDim,TNumNodes>::Dim;

template< unsigned int TDim, unsigned int TNumNodes >
constexpr unsigned int FluidElementData<TDim,TNumNodes>::NumNodes;

template< unsigned int TDim, unsigned int TNumNodes >
constexpr unsigned int FluidElementData<TDim,TNumNodes>::BlockSize;

template< unsigned int TDim, unsigned int TNumNodes >
constexpr unsigned int FluidElementData<TDim,TNumNodes>::LocalSize;

template class FluidElementData<2, 3>; // Trianlge2D3
template class FluidElementData<3, 4>; // Tetrahedra3D4
}