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
FluidElementData<TNumNodes>::FluidElementData(Geometry< Node<3> > &rGeom)
{
    for (unsigned int i = 0; i < TNumNodes; i++)
    {
        Node<3> &rNode = rGeom[i];
        const array_1d<double,3> &rNodalVel = rNode.FastGetSolutionStepValue(VELOCITY);
        const array_1d<double,3> &rNodalMeshVel = rNode.FastGetSolutionStepValue(MESH_VELOCITY);
        auto & rVel = mVectorData[Velocity];
        auto & rMeshVel = mVectorData[MeshVelocity];

        for (unsigned int d = 0; d < 3; d++)
        {
            rVel(i,d) = rNodalVel[d];
            rMeshVel(i,d) = rNodalMeshVel[d];
        }

        mScalarData[Pressure ](i) = rNode.FastGetSolutionStepValue(PRESSURE);
        mScalarData[Density  ](i) = rNode.FastGetSolutionStepValue(DENSITY);
        mScalarData[Viscosity](i) = rNode.FastGetSolutionStepValue(VISCOSITY);
    }
}

template< unsigned int TNumNodes >
FluidElementData<TNumNodes>::~FluidElementData()
{

}


template< unsigned int TNumNodes >
const array_1d<double,TNumNodes>& FluidElementData<TNumNodes>::GetNodalValues(ScalarValue Value)
{
    return mScalarData[Value];
}

template< unsigned int TNumNodes >
const boost::numeric::ublas::bounded_matrix<double, 3, TNumNodes >& FluidElementData<TNumNodes>::GetNodalValues(VectorValue Value)
{
    return mVectorData[Value];
}

template< unsigned int TNumNodes >
void FluidElementData<TNumNodes>::Evaluate(const Kratos::Vector& rShapeFunctions, ScalarValue Value, double &rOutput)
{
    const array_1d<double,TNumNodes> &rNodalData = mScalarData[Value];
    rOutput = rShapeFunctions[0] * rNodalData[0];
    for (unsigned int i = 1; i < TNumNodes; i++)
    {
        rOutput += rShapeFunctions[i] * rNodalData[Value];
    }
}

template< unsigned int TNumNodes >
void FluidElementData<TNumNodes>::Evaluate(const Kratos::Vector& rShapeFunctions, VectorValue Value, array_1d<double,3>& rOutput)
{
    const boost::numeric::ublas::bounded_matrix<double,3,TNumNodes> &rNodalData = mVectorData[Value];

    rOutput[0] = 0.0;
    rOutput[1] = 0.0;
    rOutput[2] = 0.0;

    for (unsigned int i = 0; i < TNumNodes; i++)
    {
        const double N = rShapeFunctions[i];
        for (unsigned int d = 0; d < 3; d++)
        {
            rOutput[d] += N * rNodalData(d,i);
        }
    }
}

template< unsigned int TNumNodes >
void FluidElementData<TNumNodes>::EvaluateGradient(const Kratos::Matrix& rShapeFunctionGradients, ScalarValue Value, array_1d<double,3> &rGradient)
{
    const array_1d<double,TNumNodes> &rNodalData = mScalarData[Value];
    rGradient[0] = 0.0;
    rGradient[1] = 0.0;
    rGradient[2] = 0.0;

    for (unsigned int i = 0; i < TNumNodes; i++)
    {
        for (unsigned int d = 0; d < 3; d++)
        {
            rGradient[d] += rShapeFunctionGradients(i,d) * rNodalData[i];
        }
    }
}

template< unsigned int TNumNodes >
void FluidElementData<TNumNodes>::EvaluateGradient(const Kratos::Matrix& rShapeFunctionGradients, VectorValue Value, boost::numeric::ublas::bounded_matrix<double, 3,3> &rGradient)
{
    const boost::numeric::ublas::bounded_matrix<double,3,TNumNodes> &rNodalData = mVectorData[Value];

    noalias(rGradient) = ZeroMatrix(3,3);
    
    for (unsigned int i = 0; i < TNumNodes; i++)
    {
        for (unsigned int j = 0; j < 3; j++)
        {
            for (unsigned int k = 0; k < 3; k++)
                rGradient(j,k) += rShapeFunctionGradients(i,k) * rNodalData(j,i);
        }
    }
}

template class FluidElementData<3>; // Trianlge3D3
template class FluidElementData<4>; // Tetrahedra3D4 

}