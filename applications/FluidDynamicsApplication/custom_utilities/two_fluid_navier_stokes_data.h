//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Daniel Diez, Jordi Cotela
//


#if !defined(KRATOS_TWO_FLUID_NAVIER_STOKES_DATA_H)
#define KRATOS_TWO_FLUID_NAVIER_STOKES_DATA_H

#include "includes/constitutive_law.h"

#include "fluid_dynamics_application_variables.h"
#include "custom_utilities/fluid_element_data.h"
#include "custom_utilities/element_size_calculator.h"

namespace Kratos {

///@addtogroup FluidDynamicsApplication
///@{

///@name Kratos classes
///@{

template< size_t TDim, size_t TNumNodes >
class TwoFluidNavierStokesData : public FluidElementData<TDim,TNumNodes, true>
{
public:

///@name Type Definitions
///@{

using NodalScalarData = typename FluidElementData<TDim,TNumNodes, true>::NodalScalarData;
using NodalVectorData = typename FluidElementData<TDim,TNumNodes, true>::NodalVectorData;

///@}
///@name Public Members
///@{

NodalVectorData Velocity;
NodalVectorData Velocity_OldStep1;
NodalVectorData Velocity_OldStep2;
NodalVectorData MeshVelocity;
NodalVectorData BodyForce;

NodalScalarData Pressure;
NodalScalarData Pressure_OldStep1; //BORRAR
NodalScalarData Pressure_OldStep2;  //BORRAR
NodalScalarData Distance;

double Density;
double DynamicViscosity;
double DeltaTime;      // Time increment
double DynamicTau;     // Dynamic tau considered in ASGS stabilization coefficients

double bdf0;
double bdf1;
double bdf2;

// Auxiliary containers for the symbolically-generated matrices
boost::numeric::ublas::bounded_matrix<double,TNumNodes*(TDim+1),TNumNodes*(TDim+1)> lhs;
array_1d<double,TNumNodes*(TDim+1)> rhs;
boost::numeric::ublas::bounded_matrix<double, TNumNodes*(TDim + 1), TNumNodes> V;
boost::numeric::ublas::bounded_matrix<double, TNumNodes, TNumNodes*(TDim + 1)> H;
boost::numeric::ublas::bounded_matrix<double, TNumNodes, TNumNodes> Kee;
array_1d<double, TNumNodes> rhs_ee;

double ElementSize;

ShapeFunctionsType Nenr;
ShapeDerivativesType DN_DXenr;

size_t NumPositiveNodes;
size_t NumNegativeNodes;
unsigned int NumberOfDivisions;
array_1d<double, NumNodes> PartitionsVolumes;
array_1d<double, NumNodes> PartitionsSigns; //ATTENTION: this shall be initialized of size 6 ---> Mentira

///@}
///@name Public Operations
///@{

void Initialize(const Element& rElement, const ProcessInfo& rProcessInfo) override
{
    // Base class Initialize manages constitutive law parameters
    FluidElementData<TDim,TNumNodes, true>::Initialize(rElement,rProcessInfo);
    
    const Geometry< Node<3> >& r_geometry = rElement.GetGeometry();
    const Properties& r_properties = rElement.GetProperties();
    this->FillFromNodalData(Velocity,VELOCITY,r_geometry);
    this->FillFromHistoricalNodalData(Velocity_OldStep1,VELOCITY,r_geometry,1);
    this->FillFromHistoricalNodalData(Velocity_OldStep2,VELOCITY,r_geometry,2);
	this->FillFromNodalData(Distance, DISTANCE, r_geometry);
    this->FillFromNodalData(MeshVelocity,MESH_VELOCITY,r_geometry);
    this->FillFromNodalData(BodyForce,BODY_FORCE,r_geometry);
    this->FillFromNodalData(Pressure,PRESSURE,r_geometry);
    this->FillFromHistoricalNodalData(Pressure_OldStep1,PRESSURE,r_geometry,1); //BORRAR
    this->FillFromHistoricalNodalData(Pressure_OldStep2,PRESSURE,r_geometry,2); //BORRAR
    this->FillFromProperties(Density,DENSITY,r_properties);
    this->FillFromProperties(DynamicViscosity,DYNAMIC_VISCOSITY,r_properties);
    this->FillFromProcessInfo(DeltaTime,DELTA_TIME,rProcessInfo);
    this->FillFromProcessInfo(DynamicTau,DYNAMIC_TAU,rProcessInfo);

    const Vector& BDFVector = rProcessInfo[BDF_COEFFICIENTS];
    bdf0 = BDFVector[0];
    bdf1 = BDFVector[1];
    bdf2 = BDFVector[2];

    noalias(lhs) = ZeroMatrix(TNumNodes*(TDim+1),TNumNodes*(TDim+1));
    noalias(rhs) = ZeroVector(TNumNodes*(TDim+1));
    noalias(V) = ZeroMatrix(TNumNodes*(TDim + 1), TNumNodes);
    noalias(H) = ZeroMatrix(TNumNodes, TNumNodes*(TDim + 1));
    noalias(Kee) = ZeroMatrix(TNumNodes, TNumNodes);
    noalias(rhs_ee) = ZeroVector(TNumNodes);

	NumPositiveNodes = 0;
	NumNegativeNodes = 0;
	NumberOfDivisions = 1;
}

void UpdateGeometryValues(
    double NewWeight,
    const boost::numeric::ublas::matrix_row<Kratos::Matrix> rN,
    const boost::numeric::ublas::bounded_matrix<double, TNumNodes, TDim>& rDN_DX) override
{
    FluidElementData<TDim,TNumNodes, true>::UpdateGeometryValues(NewWeight,rN,rDN_DX);
    ElementSize = ElementSizeCalculator<TDim,TNumNodes>::GradientsElementSize(rDN_DX);
}

void UpdateGeometryValues(
	double NewWeight,
	const boost::numeric::ublas::matrix_row<Kratos::Matrix> rN,
	const boost::numeric::ublas::bounded_matrix<double, TNumNodes, TDim>& rDN_DX,
	const boost::numeric::ublas::matrix_row<Kratos::Matrix> rNenr,
	const boost::numeric::ublas::bounded_matrix<double, TNumNodes, TDim>& rDN_DXenr)
{
	FluidElementData<TDim, TNumNodes, true>::UpdateGeometryValues(NewWeight, rN, rDN_DX);
	ElementSize = ElementSizeCalculator<TDim, TNumNodes>::GradientsElementSize(rDN_DX);
	noalias(this->Nenr) = rNenr;
	noalias(this->DN_DXenr) = rDN_DXenr;
}

static int Check(const Element& rElement, const ProcessInfo& rProcessInfo)
{
    const Geometry< Node<3> >& r_geometry = rElement.GetGeometry();
    
    KRATOS_CHECK_VARIABLE_KEY(VELOCITY);
	KRATOS_CHECK_VARIABLE_KEY(DISTANCE);
    KRATOS_CHECK_VARIABLE_KEY(MESH_VELOCITY);
    KRATOS_CHECK_VARIABLE_KEY(BODY_FORCE);
    KRATOS_CHECK_VARIABLE_KEY(PRESSURE);

    for (unsigned int i = 0; i < TNumNodes; i++)
    {
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VELOCITY,r_geometry[i]);
		KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DISTANCE, r_geometry[i]);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(MESH_VELOCITY,r_geometry[i]);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(BODY_FORCE,r_geometry[i]);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(PRESSURE,r_geometry[i]);
    }

    KRATOS_CHECK_VARIABLE_KEY(DENSITY);
    KRATOS_CHECK_VARIABLE_KEY(DYNAMIC_VISCOSITY);
    KRATOS_CHECK_VARIABLE_KEY(SOUND_VELOCITY);
    KRATOS_CHECK_VARIABLE_KEY(DELTA_TIME);
    KRATOS_CHECK_VARIABLE_KEY(DYNAMIC_TAU);
    KRATOS_CHECK_VARIABLE_KEY(BDF_COEFFICIENTS);

    return 0;
}

bool IsCut() {
	return (NumPositiveNodes > 0) && (NumNegativeNodes > 0);
}

bool IsAir() {
	return (NumPositiveNodes == TNumNodes);
}

void CalculateAirMaterialResponse() {
	unsigned int strain_size = 6;

	if (Dim == 2) strain_size = 3;

	if(C.size1() != strain_size)
		C.resize(strain_size,strain_size,false);
	if(ShearStress.size() != strain_size)
		ShearStress.resize(strain_size,false);

	ComputeStrain(strain_size);

	const double nu = DynamicViscosity;
	const double c1 = 2.0*nu;
	const double c2 = nu;

	//here we shall call the constitutive law
	C.clear();
	if (Dim == 2) {
		C(0, 0) = 2.0*nu;
		C(1, 1) = 2.0*nu;
		C(2, 2) = nu;

		ShearStress[0] = StrainRate[0];
		ShearStress[1] = StrainRate[1];
		ShearStress[2] = StrainRate[2];

	}

	else
	{
		C(0, 0) = 2.0*nu;
		C(1, 1) = 2.0*nu;
		C(2, 2) = 2.0*nu;
		C(3, 3) = nu;
		C(4, 4) = nu;
		C(5, 5) = nu;
		ShearStress[0] = c1*StrainRate[0];
		ShearStress[1] = c1*StrainRate[1];
		ShearStress[2] = c1*StrainRate[2];
		ShearStress[3] = c2*StrainRate[3];
		ShearStress[4] = c2*StrainRate[4];
		ShearStress[5] = c2*StrainRate[5];
	}
}

void ComputeStrain(const unsigned int& strain_size)
{
    const bounded_matrix<double, NumNodes, Dim>& v = Velocity;
    const bounded_matrix<double, NumNodes, Dim>& DN = DN_DX;
    
    // Compute strain (B*v)
    // 3D strain computation
    if (strain_size == 6)
    {
		StrainRate[0] = DN(0,0)*v(0,0) + DN(1,0)*v(1,0) + DN(2,0)*v(2,0) + DN(3,0)*v(3,0);
		StrainRate[1] = DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1) + DN(3,1)*v(3,1);
		StrainRate[2] = DN(0,2)*v(0,2) + DN(1,2)*v(1,2) + DN(2,2)*v(2,2) + DN(3,2)*v(3,2);
		StrainRate[3] = DN(0,0)*v(0,1) + DN(0,1)*v(0,0) + DN(1,0)*v(1,1) + DN(1,1)*v(1,0) + DN(2,0)*v(2,1) + DN(2,1)*v(2,0) + DN(3,0)*v(3,1) + DN(3,1)*v(3,0);
		StrainRate[4] = DN(0,1)*v(0,2) + DN(0,2)*v(0,1) + DN(1,1)*v(1,2) + DN(1,2)*v(1,1) + DN(2,1)*v(2,2) + DN(2,2)*v(2,1) + DN(3,1)*v(3,2) + DN(3,2)*v(3,1);
		StrainRate[5] = DN(0,0)*v(0,2) + DN(0,2)*v(0,0) + DN(1,0)*v(1,2) + DN(1,2)*v(1,0) + DN(2,0)*v(2,2) + DN(2,2)*v(2,0) + DN(3,0)*v(3,2) + DN(3,2)*v(3,0);
    }
    // 2D strain computation
    else if (strain_size == 3)
    {                
		StrainRate[0] = DN(0,0)*v(0,0) + DN(1,0)*v(1,0) + DN(2,0)*v(2,0);
		StrainRate[1] = DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1);
		StrainRate[2] = DN(0,1)*v(0,0) + DN(0,0)*v(0,1) + DN(1,1)*v(1,0) + DN(1,0)*v(1,1) + DN(2,1)*v(2,0) + DN(2,0)*v(2,1);
    }
}
///@}

};

///@}

///@}

}

#endif