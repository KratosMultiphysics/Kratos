//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Joaquin Gonzalez-Usua
//

// Project includes
#include "includes/cfd_variables.h"
#include "includes/dem_variables.h"
#include "includes/checks.h"

// Aplication includes
#include "d_vms_dem_coupled.h"
#include "custom_utilities/fluid_element_utilities.h"
#include "fluid_dynamics_application_variables.h"

namespace Kratos
{

///////////////////////////////////////////////////////////////////////////////////////////////////
// Life cycle

template< class TElementData >
DVMSDEMCoupled<TElementData>::DVMSDEMCoupled(IndexType NewId):
    DVMS<TElementData>(NewId)
{}

template< class TElementData >
DVMSDEMCoupled<TElementData>::DVMSDEMCoupled(IndexType NewId, const NodesArrayType& ThisNodes):
    DVMS<TElementData>(NewId,ThisNodes)
{}


template< class TElementData >
DVMSDEMCoupled<TElementData>::DVMSDEMCoupled(IndexType NewId, GeometryType::Pointer pGeometry):
    DVMS<TElementData>(NewId,pGeometry)
{}


template< class TElementData >
DVMSDEMCoupled<TElementData>::DVMSDEMCoupled(IndexType NewId, GeometryType::Pointer pGeometry, Properties::Pointer pProperties):
    DVMS<TElementData>(NewId,pGeometry,pProperties)
{}


template< class TElementData >
DVMSDEMCoupled<TElementData>::~DVMSDEMCoupled()
{}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Public Operations

template< class TElementData >
Element::Pointer DVMSDEMCoupled<TElementData>::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<DVMSDEMCoupled>(NewId, this->GetGeometry().Create(ThisNodes), pProperties);
}


template< class TElementData >
Element::Pointer DVMSDEMCoupled<TElementData>::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<DVMSDEMCoupled>(NewId, pGeom, pProperties);
}

template <class TElementData>
void DVMSDEMCoupled<TElementData>::Calculate(
    const Variable<double>& rVariable,
    double& rOutput,
    const ProcessInfo& rCurrentProcessInfo) {}

template <class TElementData>
void DVMSDEMCoupled<TElementData>::Calculate(
    const Variable<array_1d<double, 3>>& rVariable,
    array_1d<double, 3>& rOutput,
    const ProcessInfo& rCurrentProcessInfo) {
    // Lumped projection terms
    if (rVariable == ADVPROJ) {
        this->CalculateProjections(rCurrentProcessInfo);
    }
}

template <class TElementData>
void DVMSDEMCoupled<TElementData>::Calculate(
    const Variable<Vector>& rVariable,
    Vector& rOutput,
    const ProcessInfo& rCurrentProcessInfo) {}

template <class TElementData>
void DVMSDEMCoupled<TElementData>::Calculate(
    const Variable<Matrix>& rVariable,
    Matrix& rOutput,
    const ProcessInfo& rCurrentProcessInfo) {}

template <class TElementData>
void DVMSDEMCoupled<TElementData>::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    // Base class does things with constitutive law here.
    DVMS<TElementData>::Initialize(rCurrentProcessInfo);
}

template <class TElementData>
void DVMSDEMCoupled<TElementData>::FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
    DVMS<TElementData>::FinalizeSolutionStep(rCurrentProcessInfo);
}

template <class TElementData>
void DVMSDEMCoupled<TElementData>::InitializeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo)
{
    DVMS<TElementData>::InitializeNonLinearIteration(rCurrentProcessInfo);
}
///////////////////////////////////////////////////////////////////////////////////////////////////
// Inquiry

template< class TElementData >
int DVMSDEMCoupled<TElementData>::Check(const ProcessInfo &rCurrentProcessInfo) const
{
    int out = DVMS<TElementData>::Check(rCurrentProcessInfo);
    KRATOS_ERROR_IF_NOT(out == 0)
        << "Error in base class Check for Element " << this->Info() << std::endl
        << "Error code is " << out << std::endl;

    // Output variables (for Calculate() functions)
    KRATOS_CHECK_VARIABLE_KEY(SUBSCALE_VELOCITY);
    KRATOS_CHECK_VARIABLE_KEY(SUBSCALE_PRESSURE);

    return out;
}

///////////////////////////////////////////////////////////////////////////////////////////////////

template< class TElementData >
void DVMSDEMCoupled<TElementData>::GetValueOnIntegrationPoints(
    Variable<array_1d<double, 3 > > const& rVariable,
    std::vector<array_1d<double, 3 > >& rValues,
    ProcessInfo const& rCurrentProcessInfo)
{
    DVMS<TElementData>::GetValueOnIntegrationPoints(rVariable,rValues,rCurrentProcessInfo);
}

template< class TElementData >
void DVMSDEMCoupled<TElementData>::GetValueOnIntegrationPoints(
    Variable<double> const& rVariable,
    std::vector<double>& rValues,
    ProcessInfo const& rCurrentProcessInfo)
{
    DVMS<TElementData>::GetValueOnIntegrationPoints(rVariable,rValues,rCurrentProcessInfo);
}

template <class TElementData>
void DVMSDEMCoupled<TElementData>::GetValueOnIntegrationPoints(
    Variable<array_1d<double, 6>> const& rVariable,
    std::vector<array_1d<double, 6>>& rValues,
    ProcessInfo const& rCurrentProcessInfo)
{
    DVMS<TElementData>::GetValueOnIntegrationPoints(rVariable,rValues,rCurrentProcessInfo);
}

template <class TElementData>
void DVMSDEMCoupled<TElementData>::GetValueOnIntegrationPoints(
    Variable<Vector> const& rVariable,
    std::vector<Vector>& rValues,
    ProcessInfo const& rCurrentProcessInfo)
{
    DVMS<TElementData>::GetValueOnIntegrationPoints(rVariable,rValues,rCurrentProcessInfo);
}

template <class TElementData>
void DVMSDEMCoupled<TElementData>::GetValueOnIntegrationPoints(
    Variable<Matrix> const& rVariable,
    std::vector<Matrix>& rValues,
    ProcessInfo const& rCurrentProcessInfo)
{
    DVMS<TElementData>::GetValueOnIntegrationPoints(rVariable,rValues,rCurrentProcessInfo);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Input and output

template< class TElementData >
std::string DVMSDEMCoupled<TElementData>::Info() const
{
    std::stringstream buffer;
    buffer << "DVMSDEMCoupled #" << this->Id();
    return buffer.str();
}


template< class TElementData >
void DVMSDEMCoupled<TElementData>::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "DVMSDEMCoupled" << Dim << "D";
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Protected functions

///////////////////////////////////////////////////////////////////////////////////////////////////
// Evaluation of system terms on Gauss Points

template< class TElementData >
void DVMSDEMCoupled<TElementData>::AddVelocitySystem(
    TElementData& rData,
    MatrixType &rLocalLHS,
    VectorType &rLocalRHS)
{
    auto& LHS = rData.LHS;
    LHS.clear();

    const double density = this->GetAtCoordinate(rData.Density,rData.N);
    const array_1d<double,3> body_force = density * this->GetAtCoordinate(rData.BodyForce,rData.N);

    const array_1d<double,3> convective_velocity = this->FullConvectiveVelocity(rData);

    double tau_one;
    double tau_two;
    double tau_p;
    this->CalculateStabilizationParameters(rData,convective_velocity,tau_one,tau_two,tau_p);

    const double dt = rData.DeltaTime;
    DenseVector< array_1d<double,Dim> > mOldSubscaleVelocity = DVMS<TElementData>::mOldSubscaleVelocity;
    // small scale velocity contributions (subscale tracking)
    array_1d<double,Dim> OldUssTerm = (density/dt) * mOldSubscaleVelocity[rData.IntegrationPointIndex]; // rho * u_ss^{n-1}/dt

    // Old mass residual for dynamic pressure subscale: -Div(u^n)
    double OldResidual = 0.0;
    for (unsigned int a = 0; a < NumNodes; a++)
    {
        const array_1d<double,3>& rOldVel = this->GetGeometry()[a].FastGetSolutionStepValue(VELOCITY,1);
        double OldDivProj = this->GetGeometry()[a].FastGetSolutionStepValue(DIVPROJ,1);
        for (unsigned int d = 0; d < Dim; d++)
            OldResidual -= rData.DN_DX(a,d)*rOldVel[d] + rData.N[a]*OldDivProj;
    }

    Vector AGradN;
    this->ConvectionOperator(AGradN,convective_velocity,rData.DN_DX);

    // These two should be zero unless we are using OSS
    const array_1d<double,3> MomentumProj = this->GetAtCoordinate(rData.MomentumProjection,rData.N);
    const double MassProj = this->GetAtCoordinate(rData.MassProjection,rData.N);
    const double fluid_fraction = this->GetAtCoordinate(rData.FluidFraction, rData.N);
    const double fluid_fraction_rate = this->GetAtCoordinate(rData.FluidFractionRate, rData.N);
    const double mass_source = this->GetAtCoordinate(rData.MassSource, rData.N);
    array_1d<double, 3> fluid_fraction_gradient = this->GetAtCoordinate(rData.FluidFractionGradient, rData.N);
    // Multiplying convective operator by density to have correct units
    AGradN *= density;

    // Note: Dof order is (u,v,[w,]p) for each node
    for (unsigned int i = 0; i < NumNodes; i++) {

        unsigned int row = i*BlockSize;

        // LHS terms
        for (unsigned int j = 0; j < NumNodes; j++) {
            unsigned int col = j*BlockSize;

            // Some terms are the same for all velocity components, calculate them once for each i,j

            // Skew-symmetric convective term 1/2( v*grad(u)*u - grad(v) uu )
            //double K = 0.5*(rN[i]*AGradN[j] - AGradN[i]*rN[j]);
            double K = rData.N[i]*AGradN[j];

            // Stabilization: u*grad(v) * TauOne * u*grad(u) - vh * TauOne/Dt u*grad(u)
            // The last term comes from vh*d(u_ss)/dt
            K += (AGradN[i] - density*rData.N[i]/dt)*tau_one*(AGradN[j]);
            K *= rData.Weight;

            // q-p stabilization block (reset result)
            double G = 0;

            for (unsigned int d = 0; d < Dim; d++) {
                LHS(row+d,col+d) += K;

                /* v * Grad(p) block */
                // Stabilization: (a * Grad(v)) * TauOne * Grad(p)
                double GAlphaA = tau_one * AGradN[i] * fluid_fraction * rData.DN_DX(j,d);

                double AG = tau_one * AGradN[i] * rData.DN_DX(j,d);

                // From vh*d(u_ss)/dt: vh * TauOne/Dt * Grad(p)
                double VP = tau_one * density*rData.N[i]/dt * rData.DN_DX(j,d);

                // Galerkin pressure term: Div(v) * p
                double P = rData.DN_DX(i,d) * rData.N[j];

                double QAlpha = fluid_fraction * rData.DN_DX(i,d) * rData.N[j];

                double U = fluid_fraction_gradient[d] * rData.N[j] * rData.N[i];

                // Write v * Grad(p) component
                LHS(row+d,col+Dim) += rData.Weight * (AG - VP - P);
                // Use symmetry to write the q * Div(u) component
                LHS(col+Dim,row+d) += rData.Weight * (GAlphaA + U + QAlpha);

                /* q-p stabilization block */
                // Stabilization: Grad(q) * TauOne * Grad(p)
                G += fluid_fraction * rData.DN_DX(i,d) * rData.DN_DX(j,d);

                /* v-u block */
                // Stabilization: Div(v) * TauTwo * Div(u)
                for (unsigned int e = 0; e < Dim; e++){
                    double DAlphaD = fluid_fraction * (tau_two + tau_p)*rData.DN_DX(i,d)*rData.DN_DX(j,e);
                    double DU = fluid_fraction_gradient[e]*(tau_two + tau_p)*rData.DN_DX(i,d)*rData.N[j];
                    LHS(row+d,col+e) += rData.Weight * (DAlphaD + DU);
                }
            }

            // Write q-p term
            LHS(row+Dim,col+Dim) += rData.Weight * tau_one * G;
        }

        // RHS terms
        double QAlphaF = 0.0;
        for (unsigned int d = 0; d < Dim; ++d)
        {
            // v*BodyForce + v * du_ss/dt
            double VF = rData.N[i] * (body_force[d] + OldUssTerm[d]);

            // ( a * Grad(v) ) * TauOne * (Density * BodyForce - Projection)
            // vh * TauOne/Dt * f (from vh*d(uss)/dt
            double VI = tau_one * density * rData.N[i] / dt * (body_force[d] - MomentumProj[d] + OldUssTerm[d]);
            double AF = tau_one * AGradN[i] * (body_force[d] - MomentumProj[d] + OldUssTerm[d]);

            // OSS pressure subscale projection
            double DPhi = rData.DN_DX(i,d) * (tau_two + tau_p ) * (mass_source - fluid_fraction_rate - MassProj);

            // Dynamic term in pressure subscale div(vh) * h^2/(c1*dt) * (-div(uh^n) )
            rLocalRHS[row+d] -= rData.Weight * rData.DN_DX(i,d) * tau_p * OldResidual;

            rLocalRHS[row+d] += rData.Weight * (VF - VI + AF - DPhi);

            // Grad(q) * TauOne * (Density * BodyForce - Projection)
            QAlphaF += tau_one * rData.DN_DX(i, d) * fluid_fraction * (body_force[d] - MomentumProj[d] + OldUssTerm[d]);
        }
        double Q = rData.N[i] * (mass_source - fluid_fraction_rate);
        rLocalRHS[row+Dim] += rData.Weight * (QAlphaF + Q); // Grad(q) * TauOne * (Density * BodyForce)
    }

    // Write (the linearized part of the) local contribution into residual form (A*dx = b - A*x)
    array_1d<double,LocalSize> values;
    this->GetCurrentValuesVector(rData,values);
    noalias(rLocalRHS) -= prod(LHS, values);

    /* Viscous contribution (with symmetric gradient 2*nu*{E(u) - 1/3 Tr(E)} )
     * For a generic (potentially non-linear) constitutive law, one cannot assume that RHS = F - LHS*current_values.
     * Because of this, the AddViscousTerm function manages both the LHS and the RHS.
     */
    this->AddViscousTerm(rData,LHS,rLocalRHS);

    noalias(rLocalLHS) += LHS;
}

template<class TElementData>
void DVMSDEMCoupled<TElementData>::MassProjTerm(
    const TElementData& rData,
    double &rMassRHS) const
{
        const auto velocities = rData.Velocity;

        const double fluid_fraction = this->GetAtCoordinate(rData.FluidFraction, rData.N);
        const auto fluid_fraction_gradient = this->GetAtCoordinate(rData.FluidFractionGradient, rData.N);
        const double mass_source = this->GetAtCoordinate(rData.MassSource, rData.N);
        const double fluid_fraction_rate = this->GetAtCoordinate(rData.FluidFractionRate, rData.N);

        // Compute this node's contribution to the residual (evaluated at integration point)
        for (unsigned int i = 0; i < NumNodes; i++) {
            for (unsigned int d = 0; d < Dim; ++d)
            {
                rMassRHS -= (fluid_fraction * rData.DN_DX(i, d) * velocities(i, d)) + fluid_fraction_gradient[d] * rData.N[i] * velocities(i, d);
            }
        }
        rMassRHS += mass_source - fluid_fraction_rate;
}

///////////////////////////////////////////////////////////////////////////////////////////////////

template< class TElementData >
void DVMSDEMCoupled<TElementData>::AddMassStabilization(
    TElementData& rData,
    MatrixType &rMassMatrix)
{
    const double density = this->GetAtCoordinate(rData.Density,rData.N);
    const array_1d<double,3> convective_velocity = this->FullConvectiveVelocity(rData);

    double tau_one;
    double tau_two;
    double tau_p;
    this->CalculateStabilizationParameters(rData,convective_velocity,tau_one,tau_two,tau_p);

    const double dt = rData.DeltaTime;

    Vector AGradN;
    this->ConvectionOperator(AGradN,convective_velocity,rData.DN_DX);

    // Multiplying convective operator by density to have correct units
    AGradN *= density;

    const double fluid_fraction = this->GetAtCoordinate(rData.FluidFraction, rData.N);

    double W = rData.Weight * tau_one * density; // This density is for the dynamic term in the residual (rho*Du/Dt)

    // Note: Dof order is (u,v,[w,]p) for each node
    for (unsigned int i = 0; i < NumNodes; i++) {
        unsigned int row = i*BlockSize;

        for (unsigned int j = 0; j < NumNodes; j++) {
            unsigned int col = j*BlockSize;

            // u*grad(v) * TauOne * du/dt
            // v * TauOne/dt * du/dt (from v*d(uss)/dt)
            double K = W * (AGradN[i] - rData.N[i]/dt) * rData.N[j];

            for (unsigned int d = 0; d < Dim; d++)
            {
                rMassMatrix(row+d,col+d) += K;
                // grad(q) * TauOne * du/dt
                double UGAlpha = W * fluid_fraction * rData.DN_DX(i,d) * rData.N[j];
                rMassMatrix(row+Dim,col+d) += UGAlpha;
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Private functions
///////////////////////////////////////////////////////////////////////////////////////////////////

// Implementation details

// serializer

template< class TElementData >
void DVMSDEMCoupled<TElementData>::save(Serializer& rSerializer) const
{
    typedef DVMS<TElementData> BaseElement;
    DenseVector< array_1d<double,Dim> > mOldSubscaleVelocity = DVMS<TElementData>::mOldSubscaleVelocity;
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseElement );
    rSerializer.save("mOldSubscaleVelocity",mOldSubscaleVelocity);
}


template< class TElementData >
void DVMSDEMCoupled<TElementData>::load(Serializer& rSerializer)
{
    typedef DVMS<TElementData> BaseElement;
    DenseVector< array_1d<double,Dim> > mOldSubscaleVelocity = DVMS<TElementData>::mOldSubscaleVelocity;
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseElement);
    rSerializer.load("mOldSubscaleVelocity",mOldSubscaleVelocity);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Class template instantiation

template class DVMSDEMCoupled< QSVMSDEMCoupledData<2,3> >;
template class DVMSDEMCoupled< QSVMSDEMCoupledData<3,4> >;

} // namespace Kratos
