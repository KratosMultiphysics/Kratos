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
#include "utilities/math_utils.h"

// Aplication includes
#include "d_vms_dem_coupled.h"
#include "custom_utilities/qsvms_dem_coupled_data.h"
#include "custom_utilities/fluid_element_utilities.h"
#include "fluid_dynamics_application_variables.h"

namespace Kratos
{

///////////////////////////////////////////////////////////////////////////////////////////////////
// Life cycle

template< class TElementData >
DVMSDEMCoupled<TElementData>::DVMSDEMCoupled(IndexType NewId):
    DVMS<TElementData>(NewId),
    mPredictedSubscaleVelocity(),
    mOldSubscaleVelocity()
{}

template< class TElementData >
DVMSDEMCoupled<TElementData>::DVMSDEMCoupled(IndexType NewId, const NodesArrayType& ThisNodes):
    DVMS<TElementData>(NewId,ThisNodes),
    mPredictedSubscaleVelocity(),
    mOldSubscaleVelocity()
{}


template< class TElementData >
DVMSDEMCoupled<TElementData>::DVMSDEMCoupled(IndexType NewId, GeometryType::Pointer pGeometry):
    DVMS<TElementData>(NewId,pGeometry),
    mPredictedSubscaleVelocity(),
    mOldSubscaleVelocity()
{}


template< class TElementData >
DVMSDEMCoupled<TElementData>::DVMSDEMCoupled(IndexType NewId, GeometryType::Pointer pGeometry, Properties::Pointer pProperties):
    DVMS<TElementData>(NewId,pGeometry,pProperties),
    mPredictedSubscaleVelocity(),
    mOldSubscaleVelocity()
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

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
// Input and output

template <class TElementData>
void DVMSDEMCoupled<TElementData>::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    // Base class does things with constitutive law here.
    DVMS<TElementData>::Initialize(rCurrentProcessInfo);

    const unsigned int number_of_gauss_points = this->GetGeometry().IntegrationPointsNumber(this->GetIntegrationMethod());

    // The prediction is updated before each non-linear iteration:
    // It is not stored in a restart and can be safely initialized.
    mPredictedSubscaleVelocity.resize(number_of_gauss_points);
    for (unsigned int g = 0; g < number_of_gauss_points; g++)
        mPredictedSubscaleVelocity[g] = ZeroVector(Dim);

    // The old velocity may be already defined (if restarting)
    // and we want to keep the loaded values in that case.
    if (mOldSubscaleVelocity.size() != number_of_gauss_points)
    {
        mOldSubscaleVelocity.resize(number_of_gauss_points);
        for (unsigned int g = 0; g < number_of_gauss_points; g++)
            mOldSubscaleVelocity[g] = ZeroVector(Dim);
    }

}

template <class TElementData>
void DVMSDEMCoupled<TElementData>::FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
    // Get Shape function data
    Vector gauss_weights;
    Matrix shape_functions;
    ShapeFunctionDerivativesArrayType shape_function_derivatives;
    this->CalculateGeometryData(gauss_weights,shape_functions,shape_function_derivatives);
    const unsigned int number_of_integration_points = gauss_weights.size();

    TElementData data;
    data.Initialize(*this,rCurrentProcessInfo);

    for (unsigned int g = 0; g < number_of_integration_points; g++) {
        this->UpdateIntegrationPointData(data, g, gauss_weights[g],row(shape_functions,g),shape_function_derivatives[g]);

        // Not doing the update "in place" because SubscaleVelocity uses mOldSubscaleVelocity
        array_1d<double,3> UpdatedValue = ZeroVector(3);
        this->SubscaleVelocity(data,UpdatedValue);
        array_1d<double,Dim>& r_value = mOldSubscaleVelocity[g];
        for (unsigned int d = 0; d < Dim; d++) {
            r_value[d] = UpdatedValue[d];
        }
    }
}

template <class TElementData>
void DVMSDEMCoupled<TElementData>::InitializeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo)
{
    // Get Shape function data
    Vector gauss_weights;
    Matrix shape_functions;
    ShapeFunctionDerivativesArrayType shape_function_derivatives;
    this->CalculateGeometryData(gauss_weights,shape_functions,shape_function_derivatives);
    const unsigned int number_of_integration_points = gauss_weights.size();

    TElementData data;
    data.Initialize(*this,rCurrentProcessInfo);

    for (unsigned int g = 0; g < number_of_integration_points; g++) {
        this->UpdateIntegrationPointData(data, g, gauss_weights[g],row(shape_functions,g),shape_function_derivatives[g]);

        this->UpdateSubscaleVelocityPrediction(data);
    }
}

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
void DVMSDEMCoupled<TElementData>::AlgebraicMomentumResidual(
    const TElementData& rData,
    const array_1d<double,3> &rConvectionVelocity,
    array_1d<double,3>& rResidual) const
{
    const GeometryType rGeom = this->GetGeometry();

    Vector convection; // u * grad(N)
    this->ConvectionOperator(convection,rConvectionVelocity,rData.DN_DX);

    const double density = this->GetAtCoordinate(rData.Density,rData.N);
    const double viscosity = this->GetAtCoordinate(rData.DynamicViscosity, rData.N);
    BoundedMatrix<double,Dim,Dim> permeability = this->GetAtCoordinate(rData.Permeability, rData.N);
    BoundedMatrix<double,Dim,Dim> sigma = ZeroMatrix(Dim, Dim);
    const auto& r_body_forces = rData.BodyForce;
    const auto& r_velocities = rData.Velocity;
    const auto& r_pressures = rData.Pressure;

    double det_permeability = MathUtils<double>::Det(permeability);
    MathUtils<double>::InvertMatrix(permeability, sigma, det_permeability, -1.0);

    sigma *= viscosity;

    for (unsigned int i = 0; i < NumNodes; i++) {
        const array_1d<double,3>& r_acceleration = rGeom[i].FastGetSolutionStepValue(ACCELERATION);
        Vector sigma_U = ZeroVector(Dim);
        for (unsigned int d = 0; d < Dim; d++) {
            for (unsigned int e = 0; e < Dim; e++){
                sigma_U[d] += sigma(d,e) * rData.N[i] * r_velocities(i,e);
            }
            rResidual[d] += density * ( rData.N[i]*(r_body_forces(i,d) - r_acceleration[d]) - convection[i]*r_velocities(i,d)) - rData.DN_DX(i,d)*r_pressures[i] - sigma_U[d];
        }
    }
}

template< class TElementData >
void DVMSDEMCoupled<TElementData>::MomentumProjTerm(
    const TElementData& rData,
    const array_1d<double,3>& rConvectionVelocity,
    array_1d<double,3> &rMomentumRHS) const
{
    Vector AGradN;
    this->ConvectionOperator(AGradN,rConvectionVelocity,rData.DN_DX);

    const double density = this->GetAtCoordinate(rData.Density,rData.N);
    const double viscosity = this->GetAtCoordinate(rData.DynamicViscosity, rData.N);
    BoundedMatrix<double,Dim,Dim> permeability = this->GetAtCoordinate(rData.Permeability, rData.N);
    BoundedMatrix<double,Dim,Dim> sigma = ZeroMatrix(Dim, Dim);

    double det_permeability = MathUtils<double>::Det(permeability);
    MathUtils<double>::InvertMatrix(permeability, sigma, det_permeability, -1.0);

    sigma *= viscosity;

    for (unsigned int i = 0; i < NumNodes; i++) {
        array_1d<double,Dim> sigma_u = ZeroVector(Dim);
        for (unsigned int d = 0; d < Dim; d++) {
            for (unsigned int e = 0; e < Dim; e++){
                sigma_u[d] += sigma(d,e) * rData.N[i] * rData.Velocity(i,e);
            }
            rMomentumRHS[d] += density * ( rData.N[i]*(rData.BodyForce(i,d) /*- rAcc[d]*/) - AGradN[i]*rData.Velocity(i,d)) - rData.DN_DX(i,d)*rData.Pressure[i] - sigma_u[d];
        }
    }
}

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

    BoundedMatrix<double,Dim,Dim> tau_one = ZeroMatrix(Dim, Dim);
    double tau_two;
    this->CalculateStabilizationParameters(rData,convective_velocity,tau_one,tau_two);

    const double dt = rData.DeltaTime;

    // small scale velocity contributions (subscale tracking)
    array_1d<double,Dim> OldUssTerm = (density/dt) * mOldSubscaleVelocity[rData.IntegrationPointIndex]; // rho * u_ss^{n-1}/dt

    Vector AGradN;
    this->ConvectionOperator(AGradN,convective_velocity,rData.DN_DX);

    // These two should be zero unless we are using OSS
    const array_1d<double,3> MomentumProj = this->GetAtCoordinate(rData.MomentumProjection,rData.N);
    const double MassProj = this->GetAtCoordinate(rData.MassProjection,rData.N);

    double viscosity = this->GetAtCoordinate(rData.DynamicViscosity, rData.N);
    const double fluid_fraction = this->GetAtCoordinate(rData.FluidFraction, rData.N);
    const double fluid_fraction_rate = this->GetAtCoordinate(rData.FluidFractionRate, rData.N);
    const double mass_source = this->GetAtCoordinate(rData.MassSource, rData.N);

    BoundedMatrix<double,Dim,Dim> permeability = this->GetAtCoordinate(rData.Permeability, rData.N);
    BoundedMatrix<double,Dim,Dim> sigma = ZeroMatrix(Dim, Dim);
    array_1d<double, 3> fluid_fraction_gradient = this->GetAtCoordinate(rData.FluidFractionGradient, rData.N);

    double det_permeability = MathUtils<double>::Det(permeability);
    MathUtils<double>::InvertMatrix(permeability, sigma, det_permeability, -1.0);

    sigma *= viscosity;
    // Multiplying convective operator by density to have correct units
    AGradN *= density;
    // Multiplying convective operator by density to have correct units

    // Note: Dof order is (u,v,[w,]p) for each node
    for (unsigned int i = 0; i < NumNodes; i++) {

        unsigned int row = i*BlockSize;

        double divergence_convective = 0.0;
        for (unsigned int d = 0; d < Dim; d++) {
            divergence_convective += (rData.Velocity(i,d) - rData.MeshVelocity(i,d)) * rData.DN_DX(i,d);
        }

        // LHS terms
        for (unsigned int j = 0; j < NumNodes; j++) {
            unsigned int col = j*BlockSize;

            // Some terms are the same for all velocity components, calculate them once for each i,j

            // Skew-symmetric convective term 1/2( v*grad(u)*u - grad(v) uu )
            //double K = 0.5*(rN[i]*AGradN[j] - AGradN[i]*rN[j]);
            double V = rData.N[i]*AGradN[j];

            // q-p stabilization block (reset result)
            double G = 0;
            for (unsigned int d = 0; d < Dim; d++) {
                // Stabilization: u*grad(v) * TauOne * u*grad(u) - vh * TauOne/Dt u*grad(u)
                // The last term comes from vh*d(u_ss)/dt
                double AA = AGradN[i] * tau_one(d,d) * AGradN[j];

                double A = tau_one(d,d) * density * rData.N[i]/dt * AGradN[j];

                LHS(row+d,col+d) += rData.Weight * (V + AA - A);
                // Galerkin pressure term: Div(v) * p
                double P = rData.DN_DX(i,d) * rData.N[j];

                double QD = fluid_fraction * rData.DN_DX(j,d) * rData.N[i];

                double U = fluid_fraction_gradient[d] * rData.N[j] * rData.N[i];

                /* q-p stabilization block */
                // Stabilization: Grad(q) * TauOne * Grad(p)
                G += tau_one(d,d) * fluid_fraction * rData.DN_DX(i,d) * rData.DN_DX(j,d);
                /* v-u block */
                // Stabilization: Div(v) * TauTwo * Div(u)
                double GAlphaR = 0.0;
                double RSigmaG = 0.0;
                double GAlphaA = tau_one(d,d) * AGradN[j] * fluid_fraction * rData.DN_DX(i,d);
                // Stabilization: (a * Grad(v)) * TauOne * Grad(p)
                double AG = tau_one(d,d) * AGradN[i] * rData.DN_DX(j,d);
                // From vh*d(u_ss)/dt: vh * TauOne/Dt * Grad(p)
                double VP = tau_one(d,d) * density * rData.N[i] / dt * rData.DN_DX(j,d);
                for (unsigned int e = 0; e < Dim; e++){
                    double RSigma = rData.N[i] * sigma(d,e) * rData.N[j];
                    double VSigma = tau_one(d,d) * density/dt * rData.N[i] * rData.N[j] * sigma(d,e);
                    double ASigma = tau_one(d,d) * AGradN[i] * sigma(d,e) * rData.N[j];
                    double RRSigma = tau_one(d,d) * sigma(d,e) * rData.N[i] * sigma(e,d) * rData.N[j];
                    double RSigmaA = tau_one(d,d) * sigma(d,e) * rData.N[i] * AGradN[j];
                    double DD = fluid_fraction * tau_two * rData.DN_DX(i,d) * rData.DN_DX(j,e);
                    double DU = fluid_fraction_gradient[e] * tau_two * rData.DN_DX(i,d) * rData.N[j];
                    GAlphaR += tau_one(d,d) * fluid_fraction * rData.DN_DX(i,d) * sigma(d,e) * rData.N[j];
                    RSigmaG += tau_one(d,d) * sigma(d,e) * rData.N[i] * rData.DN_DX(j,e);
                    LHS(row+d,col+e) += rData.Weight * (DD + DU + RSigma - VSigma + ASigma - RRSigma - RSigmaA);
                }

                LHS(row+Dim,col+d) += rData.Weight * (GAlphaA + U + QD + GAlphaR);
                LHS(row+d,col+Dim) += rData.Weight * (AG - VP - P - RSigmaG);

            }

            // Write q-p term
            LHS(row+Dim,col+Dim) += rData.Weight * G;
        }

        // RHS terms
        double QAlphaF = 0.0;
        for (unsigned int d = 0; d < Dim; ++d)
        {
            // v*BodyForce + v * du_ss/dt
            double VF = rData.N[i] * (body_force[d] + OldUssTerm[d]);
            // ( a * Grad(v) ) * TauOne * (Density * BodyForce - Projection)
            // vh * TauOne/Dt * f (from vh*d(uss)/dt
            double VI = tau_one(d,d) * density * rData.N[i] / dt * (body_force[d] - MomentumProj[d] + OldUssTerm[d]);
            double AF = tau_one(d,d) * AGradN[i] * (body_force[d] - MomentumProj[d] + OldUssTerm[d]);
            double RSigmaF = 0.0;
            for (unsigned int e = 0; e < Dim; ++e){
                RSigmaF += tau_one(d,d) * sigma(d,e) * rData.N[i] * (body_force[e] - MomentumProj[e] + OldUssTerm[e]);
            }
            // Grad(q) * TauOne * (Density * BodyForce - Projection)
            QAlphaF += tau_one(d,d) * rData.DN_DX(i,d) * fluid_fraction * (body_force[d] - MomentumProj[d] + OldUssTerm[d]);
            // OSS pressure subscale projection
            double DPhi = rData.DN_DX(i,d) * tau_two * (mass_source - fluid_fraction_rate - MassProj);
            rLocalRHS[row+d] += rData.Weight * (VF - VI + AF + DPhi - RSigmaF);
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

    BoundedMatrix<double,Dim,Dim> tau_one = ZeroMatrix(Dim, Dim);
    double tau_two;
    this->CalculateStabilizationParameters(rData,convective_velocity,tau_one,tau_two);

    const double dt = rData.DeltaTime;

    Vector AGradN;
    this->ConvectionOperator(AGradN,convective_velocity,rData.DN_DX);

    // Multiplying convective operator by density to have correct units
    AGradN *= density;

    const double fluid_fraction = this->GetAtCoordinate(rData.FluidFraction, rData.N);
    double viscosity = this->GetAtCoordinate(rData.DynamicViscosity, rData.N);
    BoundedMatrix<double,Dim,Dim> permeability = this->GetAtCoordinate(rData.Permeability, rData.N);
    BoundedMatrix<double,Dim,Dim> sigma = ZeroMatrix(Dim, Dim);

    double det_permeability = MathUtils<double>::Det(permeability);
    MathUtils<double>::InvertMatrix(permeability, sigma, det_permeability, -1.0);

    double W = rData.Weight * density; // This density is for the dynamic term in the residual (rho*Du/Dt)
    sigma *= viscosity;

    // Note: Dof order is (u,v,[w,]p) for each node
    for (unsigned int i = 0; i < NumNodes; i++) {
        unsigned int row = i*BlockSize;

        double divergence_convective = 0.0;
        for (unsigned int d = 0; d < Dim; d++) {
            divergence_convective += (rData.Velocity(i,d) - rData.MeshVelocity(i,d)) * rData.DN_DX(i,d);
        }

        for (unsigned int j = 0; j < NumNodes; j++) {
            unsigned int col = j*BlockSize;

            for (unsigned int d = 0; d < Dim; d++)
            {
                // grad(q) * TauOne * du/dt
                double UGAlpha = tau_one(d,d) * fluid_fraction * rData.DN_DX(i,d) * rData.N[j];
                // u*grad(v) * TauOne * du/dt
                // v * TauOne/dt * du/dt (from v*d(uss)/dt)
                double AU = tau_one(d,d) * AGradN[i] * rData.N[j];
                double IU = tau_one(d,d) * rData.N[i]/dt * rData.N[j];
                for (unsigned int e = 0; e < Dim; ++e){
                    double RSigmaU = tau_one(d,d) * sigma(d,e) * rData.N[i] * rData.N[j];
                    rMassMatrix(row+d, col+e) -= W * RSigmaU;
                }
                rMassMatrix(row+d,col+d) += W * (AU - IU);
                rMassMatrix(row+Dim,col+d) += W * UGAlpha;
            }
        }
    }
}

template< class TElementData >
void DVMSDEMCoupled<TElementData>::CalculateStabilizationParameters(
    const TElementData& rData,
    const array_1d<double,3> &Velocity,
    BoundedMatrix<double,Dim,Dim> &TauOne,
    double &TauTwo) const
{
    double tau_one;
    double inv_tau;
    const double h = rData.ElementSize;
    const double density = this->GetAtCoordinate(rData.Density,rData.N);
    const double viscosity = this->GetAtCoordinate(rData.EffectiveViscosity,rData.N);
    constexpr double c1 = DVMS<TElementData>::mTauC1;
    constexpr double c2 = DVMS<TElementData>::mTauC2;
    BoundedMatrix<double,Dim,Dim> permeability = this->GetAtCoordinate(rData.Permeability, rData.N);
    BoundedMatrix<double,Dim,Dim> sigma = ZeroMatrix(Dim, Dim);
    BoundedMatrix<double,Dim,Dim> I = IdentityMatrix(Dim, Dim);

    double det_permeability = MathUtils<double>::Det(permeability);
    MathUtils<double>::InvertMatrix(permeability, sigma, det_permeability, -1.0);

    double velocity_norm = 0.0;
    double sigma_term = 0.0;
    for (unsigned int d = 0; d < Dim; d++){
        velocity_norm += Velocity[d] * Velocity[d];
        for (unsigned int e = d; e < Dim; e++){
            sigma_term += std::pow(sigma(d,e),2);
        }
    }

    velocity_norm = std::sqrt(velocity_norm);

    inv_tau = (c1 * viscosity / (h * h) + density * ( 1.0 / rData.DeltaTime + c2 * velocity_norm / h ) + viscosity * std::sqrt(sigma_term));

    tau_one = 1 / inv_tau;
    TauOne = tau_one * I;
    TauTwo = viscosity + density * c2 * velocity_norm * h / c1;
}

template< class TElementData >
void DVMSDEMCoupled<TElementData>::SubscaleVelocity(
    const TElementData& rData,
    array_1d<double,3>& rVelocitySubscale) const
{
    const double density = this->GetAtCoordinate(rData.Density,rData.N);
    array_1d<double,3> convective_velocity = this->FullConvectiveVelocity(rData);

    BoundedMatrix<double,Dim,Dim> tau_one = ZeroMatrix(Dim,Dim);
    double tau_two;
    this->CalculateStabilizationParameters(rData,convective_velocity,tau_one,tau_two);

    const double dt = rData.DeltaTime;

    array_1d<double,3> residual = ZeroVector(3);

    if (!rData.UseOSS) {
        this->AlgebraicMomentumResidual(rData,convective_velocity,residual);
    }
    else {
        this->OrthogonalMomentumResidual(rData,convective_velocity,residual);
    }

    // Note: residual is always of size 3, but stored subscale is of size Dim
    const auto& rOldSubscaleVelocity = mOldSubscaleVelocity[rData.IntegrationPointIndex];
    for (unsigned int d = 0; d < Dim; d++) {
        rVelocitySubscale[d] = tau_one(d,d)*(residual[d] + (density/dt)*rOldSubscaleVelocity[d]);
    }
}

template< class TElementData >
void DVMSDEMCoupled<TElementData>::SubscalePressure(
    const TElementData& rData,
    double &rPressureSubscale) const
{
    array_1d<double,3> convective_velocity = this->FullConvectiveVelocity(rData);

    BoundedMatrix<double,Dim,Dim> tau_one = ZeroMatrix(Dim,Dim);
    double tau_two;
    this->CalculateStabilizationParameters(rData,convective_velocity,tau_one,tau_two);

    double residual = 0.0;

    if (!rData.UseOSS)
        this->AlgebraicMassResidual(rData,residual);
    else
        this->OrthogonalMassResidual(rData,residual);

    rPressureSubscale = tau_two*residual;
}

template< class TElementData >
array_1d<double,3> DVMSDEMCoupled<TElementData>::FullConvectiveVelocity(
    const TElementData& rData) const
{
    array_1d<double,3> convective_velocity = this->GetAtCoordinate(rData.Velocity,rData.N) - this->GetAtCoordinate(rData.MeshVelocity,rData.N);
    // Adding subscale term componentwise because return type is of size 3, but subscale is of size Dim
    const array_1d<double,Dim>& r_predicted_subscale = mPredictedSubscaleVelocity[rData.IntegrationPointIndex];

    for (unsigned int d = 0; d < Dim; d++) {
        convective_velocity[d] += r_predicted_subscale[d];
    }

    return convective_velocity;
}

template< class TElementData >
void DVMSDEMCoupled<TElementData>::UpdateSubscaleVelocityPrediction(
    const TElementData& rData)
{
    const double density = this->GetAtCoordinate(rData.Density,rData.N);
    const double viscosity = this->GetAtCoordinate(rData.EffectiveViscosity,rData.N);
    array_1d<double,3> resolved_convection_velocity = this->GetAtCoordinate(rData.Velocity,rData.N) - this->GetAtCoordinate(rData.MeshVelocity,rData.N);

    const double dt = rData.DeltaTime;
    const double h = rData.ElementSize;

    // Elemental large-scale velocity gradient
    BoundedMatrix<double,Dim,Dim> resolved_velocity_gradient = ZeroMatrix(Dim,Dim);

    const auto& r_resolved_velocities = rData.Velocity;
    for (unsigned int i = 0; i < NumNodes; i++) {
        for (unsigned int m = 0; m < Dim; m++) {
            for (unsigned int n = 0; n < Dim; n++) {
                resolved_velocity_gradient(m,n) += rData.DN_DX(i,n) * r_resolved_velocities(i,m);
            }
        }
    }

    const array_1d<double,Dim>& old_subscale_velocity = mOldSubscaleVelocity[rData.IntegrationPointIndex];

    // Part of the residual that does not depend on the subscale
    array_1d<double,3> static_residual = ZeroVector(3);

    // Note I'm only using large scale convection here, small-scale convection is re-evaluated at each iteration.
    if (!rData.UseOSS)
        this->AlgebraicMomentumResidual(rData,resolved_convection_velocity,static_residual);
    else
        this->OrthogonalMomentumResidual(rData,resolved_convection_velocity,static_residual);

    // Add the time discretization term to obtain the part of the residual that does not change during iteration
    for (unsigned int d = 0; d < Dim; d++)
        static_residual[d] += density/dt * old_subscale_velocity[d];

    constexpr double c1 = DVMS<TElementData>::mTauC1;
    constexpr double c2 = DVMS<TElementData>::mTauC2;
    constexpr static double subscale_prediction_velocity_tolerance = DVMS<TElementData>::mSubscalePredictionVelocityTolerance;
    constexpr static double subscale_prediction_maximum_iterations = DVMS<TElementData>::mSubscalePredictionMaxIterations;
    constexpr static double subscale_prediction_residual_tolerance = DVMS<TElementData>::mSubscalePredictionResidualTolerance;

    // Newton-Raphson iterations for the subscale
    unsigned int iter = 0;
    bool converged = false;
    double subscale_velocity_error;

    BoundedMatrix<double,Dim,Dim> J = ZeroMatrix(Dim,Dim);
    array_1d<double,Dim> rhs = ZeroVector(Dim);
    array_1d<double,Dim> u = mPredictedSubscaleVelocity[rData.IntegrationPointIndex]; // Use last result as initial guess
    array_1d<double,Dim> du = ZeroVector(Dim);

    BoundedMatrix<double,Dim,Dim> sigma = ZeroMatrix(Dim, Dim);
    BoundedMatrix<double,Dim,Dim> I = IdentityMatrix(Dim, Dim);
    BoundedMatrix<double,Dim,Dim> permeability = this->GetAtCoordinate(rData.Permeability, rData.N);

    double det_permeability = MathUtils<double>::Det(permeability);
    MathUtils<double>::InvertMatrix(permeability, sigma, det_permeability, -1.0);
    while ( (!converged) && (iter++ < subscale_prediction_maximum_iterations) ) {

        double sigma_term = 0.0;

        // Calculate new Tau
        double convection_velocity_norm = 0.0;
        for (unsigned int d = 0; d < Dim; d++) {
            double v_d = resolved_convection_velocity[d] + u[d];
            convection_velocity_norm += v_d*v_d;
            for (unsigned int e = d; e < Dim; e++){
                sigma_term += std::pow(sigma(d,e),2);
            }
        }
        convection_velocity_norm = std::sqrt(convection_velocity_norm);

        BoundedMatrix<double,Dim,Dim> inv_tau = (c1 * viscosity / (h * h) + density * ( 1.0 / dt + c2 * convection_velocity_norm / h ) + viscosity * std::sqrt(sigma_term)) * I;

        // Newton-Raphson LHS
        noalias(J) = density * resolved_velocity_gradient;
        for (unsigned int d = 0; d < Dim; d++)
            J(d,d) += inv_tau(d,d);

        // Newton-Raphson RHS
        for (unsigned int d = 0; d < Dim; d++)
            rhs[d] = static_residual[d];
        noalias(rhs) -= prod(J,u);

        double residual_norm = rhs[0]*rhs[0];
        for (unsigned int d = 1; d < Dim; d++)
            residual_norm += rhs[d]*rhs[d];

        FluidElementUtilities<NumNodes>::DenseSystemSolve(J,rhs,du);

        // Update
        noalias(u) += du;

        // Convergence check
        subscale_velocity_error = du[0]*du[0];
        double subscale_velocity_norm = u[0]*u[0];
        for(unsigned int d = 1; d < Dim; ++d) {
            subscale_velocity_error += du[d]*du[d];
            subscale_velocity_norm += u[d]*u[d];
        }

        if (subscale_velocity_norm > subscale_prediction_velocity_tolerance)
            subscale_velocity_error /= subscale_velocity_norm;

        if ( (subscale_velocity_error <= subscale_prediction_velocity_tolerance) ||
            (residual_norm <= subscale_prediction_residual_tolerance) ) // If RHS is zero, dU is zero too, converged.
        {
            converged = true; // If RHS is zero, dU is zero too, converged.
        }
    }

    // Store new subscale values or discard the calculation
    // If not converged, we will not use the subscale in the convective term.
    noalias(mPredictedSubscaleVelocity[rData.IntegrationPointIndex]) = converged ? u : ZeroVector(Dim);
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
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseElement );
    rSerializer.save("mOldSubscaleVelocity",mOldSubscaleVelocity);
}


template< class TElementData >
void DVMSDEMCoupled<TElementData>::load(Serializer& rSerializer)
{
    typedef DVMS<TElementData> BaseElement;
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseElement);
    rSerializer.load("mOldSubscaleVelocity",mOldSubscaleVelocity);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Class template instantiation

template class DVMSDEMCoupled< QSVMSDEMCoupledData<2,3> >;
template class DVMSDEMCoupled< QSVMSDEMCoupledData<3,4> >;

template class DVMSDEMCoupled< QSVMSDEMCoupledData<2,4> >;
template class DVMSDEMCoupled< QSVMSDEMCoupledData<3,8> >;

} // namespace Kratos
