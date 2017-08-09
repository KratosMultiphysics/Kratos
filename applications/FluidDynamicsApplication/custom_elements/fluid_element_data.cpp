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

    CheckVariableKey(VELOCITY);
    CheckVariableKey(MESH_VELOCITY);
    CheckVariableKey(BODY_FORCE);
    CheckVariableKey(ADVPROJ);
    CheckVariableKey(PRESSURE);
    CheckVariableKey(DENSITY);
    CheckVariableKey(VISCOSITY);
    CheckVariableKey(DIVPROJ);

    Geometry< Node<3> >& rGeometry = rElement.GetGeometry();

    for (unsigned int i = 0; i < NumNodes; i++)
    {
        Node<3>& rNode = rGeometry[i];
        CheckVariableInNodalData(VELOCITY,rNode);
        CheckVariableInNodalData(MESH_VELOCITY,rNode);
        CheckVariableInNodalData(BODY_FORCE,rNode);
        CheckVariableInNodalData(ADVPROJ,rNode);
        CheckVariableInNodalData(PRESSURE,rNode);
        CheckVariableInNodalData(DENSITY,rNode);
        CheckVariableInNodalData(VISCOSITY,rNode);
        CheckVariableInNodalData(DIVPROJ,rNode);

        CheckDofInNode(VELOCITY_X,rNode);
        CheckDofInNode(VELOCITY_Y,rNode);
        if (Dim == 3)
            CheckDofInNode(VELOCITY_Z,rNode);
        CheckDofInNode(PRESSURE,rNode);
    }

    return 0;
    
    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
void FluidElementData<TDim, TNumNodes>::CheckVariableKey(const Kratos::VariableData& rVar)
{
    KRATOS_ERROR_IF( rVar.Key() == 0 ) << rVar.Name() << " Key is 0." << std::endl <<
    "Check that Kratos variables have been correctly registered and all required applications have been imported." << std::endl;
}

template <unsigned int TDim, unsigned int TNumNodes>
void FluidElementData<TDim, TNumNodes>::CheckVariableInNodalData(const Kratos::Variable<double>& rVar, Node<3>& rNode)
{
    KRATOS_ERROR_IF_NOT( rNode.SolutionStepsDataHas(rVar) ) << 
    "Missing " << rVar.Name() << " variable in solution step data for node " << rNode.Id() << "." << std::endl;
}

template <unsigned int TDim, unsigned int TNumNodes>
void FluidElementData<TDim, TNumNodes>::CheckVariableInNodalData(const Kratos::Variable< array_1d<double,3> >& rVar, Node<3>& rNode)
{
    KRATOS_ERROR_IF_NOT( rNode.SolutionStepsDataHas(rVar) ) << 
    "Missing " << rVar.Name() << " variable in solution step data for node " << rNode.Id() << "." << std::endl;
}

template <unsigned int TDim, unsigned int TNumNodes>
void FluidElementData<TDim, TNumNodes>::CheckDofInNode(const Kratos::VariableData& rVar, Node<3>& rNode)
{
    KRATOS_ERROR_IF_NOT( rNode.HasDofFor(rVar) ) << 
    "Missing Degree of Freedom for " << rVar.Name() << " in node " << rNode.Id() << "." << std::endl;
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