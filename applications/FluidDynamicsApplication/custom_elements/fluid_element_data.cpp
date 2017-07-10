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

template< unsigned int TNumNodes >
FluidElementData<TNumNodes>::FluidElementData(Geometry< Node<3> > &rGeom):
  mNodalData(size, TNumNodes)
{
    for (unsigned int i = 0; i < TNumNodes; i++)
    {
        Node<3> &rNode = rGeom[i];
        const array_1d<double,3> &rVel = rNode.FastGetSolutionStepValue(VELOCITY);
        mNodalData(velocity_x, i) = rVel[0];
        mNodalData(velocity_y, i) = rVel[1];
        mNodalData(velocity_z, i) = rVel[2];

        mNodalData(pressure, i)  = rNode.FastGetSolutionStepValue(PRESSURE);
        mNodalData(density, i)   = rNode.FastGetSolutionStepValue(DENSITY);
        mNodalData(viscosity, i) = rNode.FastGetSolutionStepValue(VISCOSITY);

        const array_1d<double,3> &rMeshVel = rNode.FastGetSolutionStepValue(MESH_VELOCITY);
        mNodalData(mesh_velocity_x, i) = rMeshVel[0];
        mNodalData(mesh_velocity_y, i) = rMeshVel[1];
        mNodalData(mesh_velocity_z, i) = rMeshVel[2];
    }
}

template< unsigned int TNumNodes >
FluidElementData<TNumNodes>::~FluidElementData()
{

}

template< unsigned int TNumNodes >
double FluidElementData<TNumNodes>::GetValue(NodalData Value, unsigned int node)
{
    return mNodalData(Value,node);
}

template< unsigned int TNumNodes >
void FluidElementData<TNumNodes>::Evaluate(const Kratos::Vector& rShapeFunctions, NodalData Value, double &rOutput)
{
    rOutput = rShapeFunctions[0] * mNodalData(Value,0);
    for (unsigned int i = 1; i < TNumNodes; i++)
    {
        rOutput += rShapeFunctions[i] * mNodalData(Value,i);
    }
}

template< unsigned int TNumNodes >
void FluidElementData<TNumNodes>::EvaluateGradient(const Kratos::Matrix& rShapeFunctionGradients, NodalData Value, array_1d<double,3> &rGradient)
{
    rGradient[0] = 0.0;
    rGradient[1] = 0.0;
    rGradient[2] = 0.0;

    for (unsigned int i = 0; i < TNumNodes; i++)
    {
        for (unsigned int d = 0; d < rShapeFunctionGradients.size2(); d++)
        {
            rGradient[d] += rShapeFunctionGradients(i,d) * mNodalData(Value,i);
        }
    }
}

template class FluidElementData<3>; // Trianlge3D3
template class FluidElementData<4>; // Tetrahedra3D4 

}