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


#if !defined(KRATOS_TWO_FLUID_NAVIER_STOKES_DATA_ALPHA_METHOD_H)
#define KRATOS_TWO_FLUID_NAVIER_STOKES_DATA_ALPHA_METHOD_H

#include "includes/constitutive_law.h"

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
class TwoFluidNavierStokesAlphaMethodData : public FluidElementData<TDim,TNumNodes, true>
{
public:

///@name Type Definitions
///@{

using NodalScalarData = typename FluidElementData<TDim,TNumNodes, true>::NodalScalarData;
using NodalVectorData = typename FluidElementData<TDim,TNumNodes, true>::NodalVectorData;
using ShapeFunctionsType = typename FluidElementData<TDim, TNumNodes, true>::ShapeFunctionsType;
using ShapeDerivativesType = typename FluidElementData<TDim, TNumNodes, true>::ShapeDerivativesType;
using MatrixRowType = typename FluidElementData<TDim, TNumNodes, true>::MatrixRowType;
typedef Geometry<Node> GeometryType;
typedef GeometryType::ShapeFunctionsGradientsType ShapeFunctionsGradientsType;

///@}
///@name Public Members
///@{

NodalVectorData Velocity;
NodalVectorData Velocity_OldStep1;
NodalScalarData Pressure;
NodalVectorData AccelerationAlphaMethod;

NodalVectorData MeshVelocity;
NodalVectorData MeshVelocityOldStep;

NodalVectorData BodyForce;
NodalVectorData BodyForce_OldStep1;

NodalScalarData Distance;

NodalScalarData NodalDensity;
NodalScalarData NodalDensityOldStep;
NodalScalarData NodalDynamicViscosity;
NodalScalarData NodalDynamicViscosityOldStep;

Vector ShearStressOldStep;

double Density;
double DynamicViscosity;
double DeltaTime;		   // Time increment
double DynamicTau;         // Dynamic tau considered in ASGS stabilization coefficients
double VolumeErrorRate;    // Mass loss time rate (m^3/s) to be used as source term in the mass conservation equation
double MaxSpectralRadius;
double ArtificialDynamicViscosity;

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

    const Geometry< Node >& r_geometry = rElement.GetGeometry();

    this->FillFromHistoricalNodalData(Velocity,VELOCITY,r_geometry);
    this->FillFromHistoricalNodalData(Velocity_OldStep1,VELOCITY,r_geometry,1);
    this->FillFromHistoricalNodalData(Pressure,PRESSURE,r_geometry);

    this->FillFromHistoricalNodalData(Distance, DISTANCE, r_geometry);

    this->FillFromHistoricalNodalData(MeshVelocity,MESH_VELOCITY,r_geometry);
    this->FillFromHistoricalNodalData(MeshVelocityOldStep,MESH_VELOCITY,r_geometry,1);

    this->FillFromHistoricalNodalData(BodyForce,BODY_FORCE,r_geometry);
    this->FillFromHistoricalNodalData(BodyForce_OldStep1,BODY_FORCE,r_geometry,1);

    this->FillFromHistoricalNodalData(NodalDensity, DENSITY, r_geometry);
    this->FillFromHistoricalNodalData(NodalDensityOldStep, DENSITY, r_geometry, 1);
    this->FillFromHistoricalNodalData(NodalDynamicViscosity, DYNAMIC_VISCOSITY, r_geometry);
    this->FillFromHistoricalNodalData(NodalDynamicViscosityOldStep, DYNAMIC_VISCOSITY, r_geometry, 1);
    this->FillFromNonHistoricalNodalData(AccelerationAlphaMethod,ACCELERATION,r_geometry);
    this->FillFromProcessInfo(DeltaTime,DELTA_TIME,rProcessInfo);
    this->FillFromProcessInfo(DynamicTau,DYNAMIC_TAU,rProcessInfo);
    this->FillFromProcessInfo(MaxSpectralRadius,SPECTRAL_RADIUS_LIMIT,rProcessInfo);


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

    ArtificialDynamicViscosity = r_geometry.Has(ARTIFICIAL_DYNAMIC_VISCOSITY) ? r_geometry.GetValue(ARTIFICIAL_DYNAMIC_VISCOSITY) : 0.0;

    // In here we calculate the volume error temporary ratio (note that the input value is a relative measure of the volume loss)
    // Also note that we do consider time varying time step but a constant theta (we incur in a small error when switching from BE to CN)
    // Note as well that there is a minus sign (this comes from the divergence sign)
    if (IsCut()) {
        // Get the previous time increment. Note that we check its value in case the previous ProcessInfo is empty (e.g. first step)
        double previous_dt = rProcessInfo.GetPreviousTimeStepInfo()[DELTA_TIME];
        if (previous_dt < 1.0e-12) {
            previous_dt = rProcessInfo[DELTA_TIME];
        }
        // Get the absolute volume error from the ProcessInfo and calculate the time rate
        this->FillFromProcessInfo(VolumeErrorRate,VOLUME_ERROR,rProcessInfo);
        VolumeErrorRate /= -previous_dt;
    } else {
        VolumeErrorRate = 0.0;
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
    CalculateDensityAtGaussPoint();
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
    CalculateDensityAtGaussPoint();
}

static int Check(const Element& rElement, const ProcessInfo& rProcessInfo)
{
    const Geometry< Node >& r_geometry = rElement.GetGeometry();

    for (unsigned int i = 0; i < TNumNodes; i++)
    {
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VELOCITY,r_geometry[i]);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DISTANCE, r_geometry[i]);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(MESH_VELOCITY,r_geometry[i]);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(BODY_FORCE,r_geometry[i]);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(PRESSURE,r_geometry[i]);
    }

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

    if constexpr (TDim == 2)
    {
        const double trace = strain[0] + strain[1];
        const double volumetric_part = trace/2.0; // Note: this should be small for an incompressible fluid (it is basically the incompressibility error)

        stress[0] = c1 * (strain[0] - volumetric_part);
        stress[1] = c1 * (strain[1] - volumetric_part);
        stress[2] = c2 * strain[2];
    }

    else if constexpr (TDim == 3)
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
    const double rho_inf=this->MaxSpectralRadius;
    const double alpha_f= 1/(rho_inf+1);
    const BoundedMatrix<double, TNumNodes, TDim>& v = Velocity_OldStep1+alpha_f*(Velocity-Velocity_OldStep1);
    const BoundedMatrix<double, TNumNodes, TDim>& DN = this->DN_DX;

    // Compute strain (B*v)
    // 3D strain computation
    if constexpr (TDim == 3)
    {
        this->StrainRate[0] = DN(0,0)*v(0,0) + DN(1,0)*v(1,0) + DN(2,0)*v(2,0) + DN(3,0)*v(3,0);
        this->StrainRate[1] = DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1) + DN(3,1)*v(3,1);
        this->StrainRate[2] = DN(0,2)*v(0,2) + DN(1,2)*v(1,2) + DN(2,2)*v(2,2) + DN(3,2)*v(3,2);
        this->StrainRate[3] = DN(0,0)*v(0,1) + DN(0,1)*v(0,0) + DN(1,0)*v(1,1) + DN(1,1)*v(1,0) + DN(2,0)*v(2,1) + DN(2,1)*v(2,0) + DN(3,0)*v(3,1) + DN(3,1)*v(3,0);
        this->StrainRate[4] = DN(0,1)*v(0,2) + DN(0,2)*v(0,1) + DN(1,1)*v(1,2) + DN(1,2)*v(1,1) + DN(2,1)*v(2,2) + DN(2,2)*v(2,1) + DN(3,1)*v(3,2) + DN(3,2)*v(3,1);
        this->StrainRate[5] = DN(0,0)*v(0,2) + DN(0,2)*v(0,0) + DN(1,0)*v(1,2) + DN(1,2)*v(1,0) + DN(2,0)*v(2,2) + DN(2,2)*v(2,0) + DN(3,0)*v(3,2) + DN(3,2)*v(3,0);
    }
    // 2D strain computation
    else if constexpr (TDim == 2)
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
    if constexpr (TDim == 3)
    {
        strain_rate_norm = std::sqrt(2.*S[0] * S[0] + 2.*S[1] * S[1] + 2.*S[2] * S[2] +
            S[3] * S[3] + S[4] * S[4] + S[5] * S[5]);
    }

    else if constexpr (TDim == 2)
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
    this->EffectiveViscosity = DynamicViscosity + ArtificialDynamicViscosity;
}

void ComputeDarcyTerm()
{
    //TODO: We need to implement this in order to do the explicit template instantiation in the base TwoFluidNavierStokesElement
    //TODO: Properly implement it (with the required member variables) once we add the Darcy contribution to the Alpha method element
}

///@}

};

///@}

///@}

}

#endif