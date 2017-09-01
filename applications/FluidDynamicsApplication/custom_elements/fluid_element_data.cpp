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
#include "includes/cfd_variables.h"
#include "includes/checks.h"

namespace Kratos
{
template <unsigned int TDim, unsigned int TNumNodes>
FluidElementData<TDim, TNumNodes>::FluidElementData(Geometry<Node<3>>& rGeometry)
    : Pressure(TNumNodes, 0.0), Density(TNumNodes, 0.0), Viscosity(TNumNodes, 0.0), MassProjection(TNumNodes, 0.0)
{
    for (unsigned int i = 0; i < TNumNodes; i++)
    {
        Node<3>& rNode = rGeometry[i];
        const array_1d<double, 3>& rNodalVel = rNode.FastGetSolutionStepValue(VELOCITY);
        const array_1d<double, 3>& rNodalMeshVel = rNode.FastGetSolutionStepValue(MESH_VELOCITY);
        const array_1d<double, 3>& rNodalBodyForce = rNode.FastGetSolutionStepValue(BODY_FORCE);
        const array_1d<double, 3>& rNodalMomentumProjection = rNode.FastGetSolutionStepValue(ADVPROJ);

        for (unsigned int d = 0; d < TDim; d++)
        {
            Velocity(i, d) = rNodalVel[d];
            MeshVelocity(i, d) = rNodalMeshVel[d];
            BodyForce(i,d) = rNodalBodyForce[d];
            MomentumProjection(i,d) = rNodalMomentumProjection[d];
        }

        Pressure[i] = rNode.FastGetSolutionStepValue(PRESSURE);
        Density[i] = rNode.FastGetSolutionStepValue(DENSITY);
        Viscosity[i] = rNode.FastGetSolutionStepValue(VISCOSITY);
        MassProjection[i] = rNode.FastGetSolutionStepValue(DIVPROJ);
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
int FluidElementData<TDim, TNumNodes>::Check(Element& rElement)
{
    KRATOS_TRY;

    KRATOS_CHECK_VARIABLE_KEY(VELOCITY);
    KRATOS_CHECK_VARIABLE_KEY(MESH_VELOCITY);
    KRATOS_CHECK_VARIABLE_KEY(BODY_FORCE);
    KRATOS_CHECK_VARIABLE_KEY(ADVPROJ);
    KRATOS_CHECK_VARIABLE_KEY(PRESSURE);
    KRATOS_CHECK_VARIABLE_KEY(DENSITY);
    KRATOS_CHECK_VARIABLE_KEY(VISCOSITY);
    KRATOS_CHECK_VARIABLE_KEY(DIVPROJ);

    Geometry< Node<3> >& rGeometry = rElement.GetGeometry();

    for (unsigned int i = 0; i < NumNodes; i++)
    {
        Node<3>& rNode = rGeometry[i];
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VELOCITY,rNode);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(MESH_VELOCITY,rNode);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(BODY_FORCE,rNode);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(ADVPROJ,rNode);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(PRESSURE,rNode);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DENSITY,rNode);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VISCOSITY,rNode);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DIVPROJ,rNode);
    }

    return 0;
    
    KRATOS_CATCH("");
}

///////////////////////////////////////////////////////////////////////////////

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