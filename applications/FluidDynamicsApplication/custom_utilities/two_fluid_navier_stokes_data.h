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
#include "custom_utilities/fluid_element_utilities.h"

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
using ShapeFunctionsType = typename FluidElementData<TDim, TNumNodes, true>::ShapeFunctionsType;
using ShapeDerivativesType = typename FluidElementData<TDim, TNumNodes, true>::ShapeDerivativesType;
using MatrixRowType = typename FluidElementData<TDim, TNumNodes, true>::MatrixRowType;
typedef Geometry<Node<3>> GeometryType;
typedef GeometryType::ShapeFunctionsGradientsType ShapeFunctionsGradientsType;

///@}
///@name Public Members
///@{

NodalVectorData Velocity;
NodalVectorData Velocity_OldStep1;
NodalVectorData Velocity_OldStep2;
NodalVectorData MeshVelocity;
NodalVectorData BodyForce;

NodalScalarData Pressure;
NodalScalarData Distance;
NodalScalarData NodalDensity;
NodalScalarData NodalDynamicViscosity;

double Density;
double DynamicViscosity;
double DeltaTime;		   // Time increment
double DynamicTau;         // Dynamic tau considered in ASGS stabilization coefficients
double SmagorinskyConstant;

double bdf0;
double bdf1;
double bdf2;

// Auxiliary containers for the symbolically-generated matrices
BoundedMatrix<double,TNumNodes*(TDim+1),TNumNodes*(TDim+1)> lhs;
array_1d<double,TNumNodes*(TDim+1)> rhs;
BoundedMatrix<double, TNumNodes*(TDim + 1), TNumNodes> V;
BoundedMatrix<double, TNumNodes, TNumNodes*(TDim + 1)> H;
BoundedMatrix<double, TNumNodes, TNumNodes> Kee;
array_1d<double, TNumNodes> rhs_ee;

double ElementSize;

Matrix N_pos_side;
Matrix N_neg_side;
ShapeFunctionsGradientsType DN_DX_pos_side;
ShapeFunctionsGradientsType DN_DX_neg_side;

BoundedMatrix<double,TNumNodes,TNumNodes> Enr_Pos_Interp;
BoundedMatrix<double,TNumNodes,TNumNodes> Enr_Neg_Interp;

Vector w_gauss_pos_side;
Vector w_gauss_neg_side;

ShapeFunctionsType Nenr;
ShapeDerivativesType DN_DXenr;

size_t NumPositiveNodes;
size_t NumNegativeNodes;
unsigned int NumberOfDivisions;


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
    this->FillFromNodalData(NodalDensity, DENSITY, r_geometry);
    this->FillFromNodalData(NodalDynamicViscosity, DYNAMIC_VISCOSITY, r_geometry);
    this->FillFromProperties(SmagorinskyConstant, C_SMAGORINSKY, r_properties);
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

    for (unsigned int i = 0; i < TNumNodes; i++){
        if(Distance[i] > 0)
            NumPositiveNodes++;
        else
            NumNegativeNodes++;
    }
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

void UpdateGeometryValues(
    unsigned int IntegrationPointIndex,
	double NewWeight,
	const MatrixRowType& rN,
	const BoundedMatrix<double, TNumNodes, TDim>& rDN_DX,
	const MatrixRowType& rNenr,
	const BoundedMatrix<double, TNumNodes, TDim>& rDN_DXenr)
{
	FluidElementData<TDim, TNumNodes, true>::UpdateGeometryValues(IntegrationPointIndex, NewWeight, rN, rDN_DX);
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
	const unsigned int strain_size = 3 * (TDim - 1);

	if(this->C.size1() != strain_size)
		this->C.resize(strain_size,strain_size,false);
	if(this->ShearStress.size() != strain_size)
		this->ShearStress.resize(strain_size,false);

    ComputeStrain();

    CalculateEffectiveViscosityAtGaussPoint();

	const double mu = this->EffectiveViscosity;
	const double c1 = 2.0*mu;
	const double c2 = mu;

	this->C.clear();
    BoundedMatrix<double, strain_size, strain_size> c_mat = this->C;
    Vector& stress = this->ShearStress;
    Vector& strain = this->StrainRate;

    FluidElementUtilities<TNumNodes>::GetNewtonianConstitutiveMatrix(mu, c_mat);
    this->C = c_mat;

	if (TDim == 2)
	{
        const double trace = strain[0] + strain[1];
        const double volumetric_part = trace/3.0; // Note: this should be small for an incompressible fluid (it is basically the incompressibility error)

		stress[0] = c1 * strain[0] - volumetric_part;
		stress[1] = c1 * strain[1] - volumetric_part;
		stress[2] = c2 * strain[2];
	}

	else if (TDim == 3)
	{
        const double trace = strain[0] + strain[1] + strain[2];
        const double volumetric_part = trace/3.0; // Note: this should be small for an incompressible fluid (it is basically the incompressibility error)

		stress[0] = c1*(strain[0] - volumetric_part);
		stress[1] = c1*(strain[1] - volumetric_part);
		stress[2] = c1*(strain[2] - volumetric_part);
		stress[3] = c2*strain[3];
		stress[4] = c2*strain[4];
		stress[5] = c2*strain[5];
	}
}

void ComputeStrain()
{
    const BoundedMatrix<double, TNumNodes, TDim>& v = Velocity;
    const BoundedMatrix<double, TNumNodes, TDim>& DN = this->DN_DX;

    // Compute strain (B*v)
    // 3D strain computation
    if (TDim == 3)
    {
		this->StrainRate[0] = DN(0,0)*v(0,0) + DN(1,0)*v(1,0) + DN(2,0)*v(2,0) + DN(3,0)*v(3,0);
		this->StrainRate[1] = DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1) + DN(3,1)*v(3,1);
		this->StrainRate[2] = DN(0,2)*v(0,2) + DN(1,2)*v(1,2) + DN(2,2)*v(2,2) + DN(3,2)*v(3,2);
		this->StrainRate[3] = DN(0,0)*v(0,1) + DN(0,1)*v(0,0) + DN(1,0)*v(1,1) + DN(1,1)*v(1,0) + DN(2,0)*v(2,1) + DN(2,1)*v(2,0) + DN(3,0)*v(3,1) + DN(3,1)*v(3,0);
		this->StrainRate[4] = DN(0,1)*v(0,2) + DN(0,2)*v(0,1) + DN(1,1)*v(1,2) + DN(1,2)*v(1,1) + DN(2,1)*v(2,2) + DN(2,2)*v(2,1) + DN(3,1)*v(3,2) + DN(3,2)*v(3,1);
		this->StrainRate[5] = DN(0,0)*v(0,2) + DN(0,2)*v(0,0) + DN(1,0)*v(1,2) + DN(1,2)*v(1,0) + DN(2,0)*v(2,2) + DN(2,2)*v(2,0) + DN(3,0)*v(3,2) + DN(3,2)*v(3,0);
    }
    // 2D strain computation
    else if (TDim == 2)
    {
		this->StrainRate[0] = DN(0,0)*v(0,0) + DN(1,0)*v(1,0) + DN(2,0)*v(2,0);
		this->StrainRate[1] = DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1);
		this->StrainRate[2] = DN(0,1)*v(0,0) + DN(0,0)*v(0,1) + DN(1,1)*v(1,0) + DN(1,0)*v(1,1) + DN(2,1)*v(2,0) + DN(2,0)*v(2,1);
    }
}

double ComputeStrainNorm()
{
    double strain_rate_norm;
    Vector& S = this->StrainRate;
    if (TDim == 3)
    {
        strain_rate_norm = std::sqrt(2.*S[0] * S[0] + 2.*S[1] * S[1] + 2.*S[2] * S[2] +
            S[3] * S[3] + S[4] * S[4] + S[5] * S[5]);
    }

    else if (TDim == 2)
    {
        strain_rate_norm = std::sqrt(2.*S[0] * S[0] + 2.*S[1] * S[1] + S[2] * S[2]);
    }
    return strain_rate_norm;
}

void CalculateDensityAtGaussPoint()
{
    double dist = 0.0;
    for (unsigned int i = 0; i < TNumNodes; i++)
        dist += this->N[i] * Distance[i];

    int navg = 0;
    double density = 0.0;
    for (unsigned int i = 0; i < TNumNodes; i++)
    {
        if (dist * Distance[i] > 0.0)
        {
            navg += 1;
            density += NodalDensity[i];
        }
    }

    Density = density / navg;
}

void CalculateEffectiveViscosityAtGaussPoint()
{
    double dist = 0.0;
    for (unsigned int i = 0; i < TNumNodes; i++)
        dist += this->N[i] * Distance[i];

    int navg = 0;
    double dynamic_viscosity = 0.0;
    for (unsigned int i = 0; i < TNumNodes; i++)
    {
        if (dist * Distance[i] > 0.0)
        {
            navg += 1;
            dynamic_viscosity += NodalDynamicViscosity[i];
        }
    }
    DynamicViscosity = dynamic_viscosity / navg;


    if (SmagorinskyConstant > 0.0)
    {
        const double strain_rate_norm = ComputeStrainNorm();

        double length_scale = SmagorinskyConstant*ElementSize;
        length_scale *= length_scale; // square
        this->EffectiveViscosity = DynamicViscosity + 2.0*length_scale*strain_rate_norm;
    }
    else this->EffectiveViscosity = DynamicViscosity;
}
///@}





};

///@}

///@}

}

#endif