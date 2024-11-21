//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Uxue Chasco
//


#if !defined(KRATOS_VECTORIAL_CONVECTION_FRACTIONAL_ELEMENT_DATA_H)
#define KRATOS_VECTORIAL_CONVECTION_FRACTIONAL_ELEMENT_DATA_H

#include "fluid_dynamics_application_variables.h"
#include "custom_utilities/fluid_element_data.h"
#include "utilities/element_size_calculator.h"
#include "custom_utilities/fluid_element_utilities.h"

namespace Kratos {

///@addtogroup FluidDynamicsApplication
///@{

///@name Kratos classes
///@{

template< size_t TDim, size_t TNumNodes >
class VectorialConvectionFractionalElementData : public FluidElementData<TDim,TNumNodes, true>
{
public:

///@name Type Definitions
///@{

using NodalScalarData = typename FluidElementData<TDim,TNumNodes, true>::NodalScalarData;
using NodalVectorData = typename FluidElementData<TDim,TNumNodes, true>::NodalVectorData;
using ShapeFunctionsType = typename FluidElementData<TDim, TNumNodes, true>::ShapeFunctionsType;
using ShapeDerivativesType = typename FluidElementData<TDim, TNumNodes, true>::ShapeDerivativesType;
using MatrixRowType = typename FluidElementData<TDim, TNumNodes, true>::MatrixRowType;

///@}
///@name Public Members
///@{


NodalVectorData Velocity_Fractional;
NodalVectorData Velocity_OldStep1;
NodalVectorData Velocity_OldStep2;
NodalVectorData Velocity_OldStep3;

NodalVectorData MeshVelocity;

double DynamicViscosity;
double DeltaTime;		   // Time increment
double DynamicTau;         // Dynamic tau considered in ASGS stabilization coefficients
double bdf0;
double bdf1;
double bdf2;

// Auxiliary containers for the symbolically-generated matrices
BoundedMatrix<double,TNumNodes*(TDim),TNumNodes*(TDim)> lhs;
array_1d<double,TNumNodes*(TDim)> rhs;

double ElementSize;



///@}
///@name Public Operations
///@{

void Initialize(const Element& rElement, const ProcessInfo& rProcessInfo) override
{


    const Geometry< Node >& r_geometry = rElement.GetGeometry();
    const Properties& r_properties = rElement.GetProperties();
    this->FillFromHistoricalNodalData(Velocity_OldStep1,VELOCITY,r_geometry,1);
    this->FillFromHistoricalNodalData(Velocity_OldStep2,VELOCITY,r_geometry,2);
    this->FillFromHistoricalNodalData(Velocity_OldStep3, VELOCITY, r_geometry, 3);
    this->FillFromHistoricalNodalData(MeshVelocity,MESH_VELOCITY,r_geometry);
    this->FillFromHistoricalNodalData(Velocity_Fractional, FRACTIONAL_VELOCITY, r_geometry, 0);
    this->FillFromProcessInfo(DeltaTime,DELTA_TIME,rProcessInfo);
    this->FillFromProcessInfo(DynamicTau,DYNAMIC_TAU,rProcessInfo);
    const Vector& BDFVector = rProcessInfo[BDF_COEFFICIENTS];

    bdf0 = BDFVector[0];
    bdf1 = BDFVector[1];
    bdf2 = BDFVector[2];
    KRATOS_WATCH(BDFVector)

    noalias(lhs) = ZeroMatrix(TNumNodes*(TDim),TNumNodes*(TDim));
    noalias(rhs) = ZeroVector(TNumNodes*(TDim));
}

void UpdateGeometryValues(
    unsigned int IntegrationPointIndex,
    double NewWeight,
    const MatrixRowType& rN,
    const BoundedMatrix<double, TNumNodes, TDim>& rDN_DX) override
{
    FluidElementData<TDim,TNumNodes, true>::UpdateGeometryValues(IntegrationPointIndex, NewWeight,rN,rDN_DX);
    ElementSize = ElementSizeCalculator<TDim,TNumNodes>::GradientsElementSize(rDN_DX);
}


static int Check(const Element& rElement, const ProcessInfo& rProcessInfo)
{
    const Geometry< Node >& r_geometry = rElement.GetGeometry();

    for (unsigned int i = 0; i < TNumNodes; i++)
    {

        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(FRACTIONAL_VELOCITY, r_geometry[i]);

    }

    return 0;
}

///@}





};

///@}

///@}

}

#endif