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
#include "alternative_d_vms_dem_coupled.h"
#include "custom_utilities/qsvms_dem_coupled_data.h"
#include "custom_utilities/fluid_element_utilities.h"
#include "fluid_dynamics_application_variables.h"

namespace Kratos
{

///////////////////////////////////////////////////////////////////////////////////////////////////
// Life cycle

template< class TElementData >
AlternativeDVMSDEMCoupled<TElementData>::AlternativeDVMSDEMCoupled(IndexType NewId):
    DVMS<TElementData>(NewId),
    mPredictedSubscaleVelocity(),
    mOldSubscaleVelocity(),
    mPreviousVelocity()
{}

template< class TElementData >
AlternativeDVMSDEMCoupled<TElementData>::AlternativeDVMSDEMCoupled(IndexType NewId, const NodesArrayType& ThisNodes):
    DVMS<TElementData>(NewId,ThisNodes),
    mPredictedSubscaleVelocity(),
    mOldSubscaleVelocity(),
    mPreviousVelocity()
{}


template< class TElementData >
AlternativeDVMSDEMCoupled<TElementData>::AlternativeDVMSDEMCoupled(IndexType NewId, GeometryType::Pointer pGeometry):
    DVMS<TElementData>(NewId,pGeometry),
    mPredictedSubscaleVelocity(),
    mOldSubscaleVelocity(),
    mPreviousVelocity()
{}


template< class TElementData >
AlternativeDVMSDEMCoupled<TElementData>::AlternativeDVMSDEMCoupled(IndexType NewId, GeometryType::Pointer pGeometry, Properties::Pointer pProperties):
    DVMS<TElementData>(NewId,pGeometry,pProperties),
    mPredictedSubscaleVelocity(),
    mOldSubscaleVelocity(),
    mPreviousVelocity()
{}


template< class TElementData >
AlternativeDVMSDEMCoupled<TElementData>::~AlternativeDVMSDEMCoupled()
{}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Public Operations

template< class TElementData >
Element::Pointer AlternativeDVMSDEMCoupled<TElementData>::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<AlternativeDVMSDEMCoupled>(NewId, this->GetGeometry().Create(ThisNodes), pProperties);
}


template< class TElementData >
Element::Pointer AlternativeDVMSDEMCoupled<TElementData>::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<AlternativeDVMSDEMCoupled>(NewId, pGeom, pProperties);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
// Input and output

template <class TElementData>
void AlternativeDVMSDEMCoupled<TElementData>::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    // Base class does things with constitutive law here.
    DVMS<TElementData>::Initialize(rCurrentProcessInfo);

    const unsigned int number_of_gauss_points = this->GetGeometry().IntegrationPointsNumber(this->GetIntegrationMethod());

    // The prediction is updated before each non-linear iteration:
    // It is not stored in a restart and can be safely initialized.
    mPreviousVelocity.resize(number_of_gauss_points);

    // The old velocity may be already defined (if restarting)
    // and we want to keep the loaded values in that case.
    if (mPreviousVelocity.size() != number_of_gauss_points)
    {
        mPreviousVelocity.resize(number_of_gauss_points);
        for (unsigned int g = 0; g < number_of_gauss_points; g++)
            mPreviousVelocity[g] = ZeroVector(Dim);
    }

    mPredictedSubscaleVelocity.resize(number_of_gauss_points);
    for (unsigned int g = 0; g < number_of_gauss_points; g++)
    // The old velocity may be already defined (if restarting)
    // and we want to keep the loaded values in that case.
    if (mOldSubscaleVelocity.size() != number_of_gauss_points)
    {
        mOldSubscaleVelocity.resize(number_of_gauss_points);
        for (unsigned int g = 0; g < number_of_gauss_points; g++)
            mOldSubscaleVelocity[g] = ZeroVector(Dim);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////

template <class TElementData>
void AlternativeDVMSDEMCoupled<TElementData>::FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
    // Get Shape function data
    Vector gauss_weights;
    Matrix shape_functions;
    ShapeFunctionDerivativesArrayType shape_function_derivatives;
    this->CalculateGeometryData(gauss_weights,shape_functions,shape_function_derivatives);
    const unsigned int number_of_integration_points = gauss_weights.size();

    TElementData data;
    data.Initialize(*this,rCurrentProcessInfo);
    array_1d<double,3> UpdatedValue;
    for (unsigned int g = 0; g < number_of_integration_points; g++) {
        this->UpdateIntegrationPointData(data, g, gauss_weights[g],row(shape_functions,g),shape_function_derivatives[g]);

        // Not doing the update "in place" because SubscaleVelocity uses mOldSubscaleVelocity
        UpdatedValue = ZeroVector(3);
        this->SubscaleVelocity(data,UpdatedValue);
        array_1d<double,Dim>& r_value = mOldSubscaleVelocity[g];
        for (unsigned int d = 0; d < Dim; d++) {
            r_value[d] = UpdatedValue[d];
        }
    }
}

template <class TElementData>
void AlternativeDVMSDEMCoupled<TElementData>::InitializeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo)
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

template <class TElementData>
void AlternativeDVMSDEMCoupled<TElementData>::FinalizeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo)
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

        this->UpdateSubscaleVelocity(data);
    }
}

template< class TElementData >
std::string AlternativeDVMSDEMCoupled<TElementData>::Info() const
{
    std::stringstream buffer;
    buffer << "AlternativeDVMSDEMCoupled #" << this->Id();
    return buffer.str();
}


template< class TElementData >
void AlternativeDVMSDEMCoupled<TElementData>::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "AlternativeDVMSDEMCoupled" << Dim << "D";
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Protected functions

///////////////////////////////////////////////////////////////////////////////////////////////////
// Evaluation of system terms on Gauss Points

template< class TElementData >
void AlternativeDVMSDEMCoupled<TElementData>::AlgebraicMomentumResidual(
    const TElementData& rData,
    const array_1d<double,3> &rConvectionVelocity,
    array_1d<double,3>& rResidual) const
{
    const GeometryType rGeom = this->GetGeometry();

    Vector convection; // u * grad(N)
    this->ConvectionOperator(convection,rConvectionVelocity,rData.DN_DX);

    const double density = this->GetAtCoordinate(rData.Density,rData.N);
    const double viscosity = this->GetAtCoordinate(rData.DynamicViscosity, rData.N);
    const double fluid_fraction = this->GetAtCoordinate(rData.FluidFraction, rData.N);
    BoundedMatrix<double,Dim,Dim> permeability = this->GetAtCoordinate(rData.Permeability, rData.N);
    BoundedMatrix<double,Dim,Dim> sigma = ZeroMatrix(Dim, Dim);

    const auto& r_body_forces = rData.BodyForce;
    const auto& r_velocities = rData.Velocity;
    const auto& r_pressures = rData.Pressure;
    const auto& r_fluid_fraction_gradient = this->GetAtCoordinate(rData.FluidFractionGradient, rData.N);

    double det_permeability = MathUtils<double>::Det(permeability);
    MathUtils<double>::InvertMatrix(permeability, sigma, det_permeability, -1.0);

    sigma *= viscosity;
    Vector sigma_U, grad_alpha_sym_grad_u;
    BoundedMatrix<double,Dim,Dim> sym_gradient_u;
    for (unsigned int i = 0; i < NumNodes; i++) {
        const array_1d<double,3>& r_acceleration = rGeom[i].FastGetSolutionStepValue(ACCELERATION);
        sigma_U = ZeroVector(Dim);
        sym_gradient_u = ZeroMatrix(Dim, Dim);
        grad_alpha_sym_grad_u = ZeroVector(Dim);
        for (unsigned int d = 0; d < Dim; d++) {
            double div_u = 0.0;
            for (unsigned int e = 0; e < Dim; e++){
                sigma_U[d] += sigma(d,e) * rData.N[i] * r_velocities(i,e);
                sym_gradient_u(d,e) += 1.0 / 2.0 * (rData.DN_DX(i,d) * r_velocities(i,e) + rData.DN_DX(i,e) * r_velocities(i,d));
                grad_alpha_sym_grad_u[d] += r_fluid_fraction_gradient[d] * sym_gradient_u(d,e);
                div_u += rData.DN_DX(i,e) * r_velocities(i,e);
            }
            rResidual[d] += density * (rData.N[i] * r_body_forces(i,d) - fluid_fraction * rData.N[i] * r_acceleration[d] - fluid_fraction * convection[i] * r_velocities(i,d)) + 2 * grad_alpha_sym_grad_u[d] * viscosity - 2.0 / 3.0 * viscosity * r_fluid_fraction_gradient[d] * div_u - fluid_fraction * rData.DN_DX(i,d) * r_pressures[i] - sigma_U[d];
        }
    }
}

template< class TElementData >
void AlternativeDVMSDEMCoupled<TElementData>::MomentumProjTerm(
    const TElementData& rData,
    const array_1d<double,3>& rConvectionVelocity,
    array_1d<double,3> &rMomentumRHS) const
{
    Vector AGradN;
    this->ConvectionOperator(AGradN,rConvectionVelocity,rData.DN_DX);

    const double density = this->GetAtCoordinate(rData.Density,rData.N);
    const double viscosity = this->GetAtCoordinate(rData.DynamicViscosity, rData.N);
    const double fluid_fraction = this->GetAtCoordinate(rData.FluidFraction, rData.N);
    BoundedMatrix<double,Dim,Dim> permeability = this->GetAtCoordinate(rData.Permeability, rData.N);
    BoundedMatrix<double,Dim,Dim> sigma = ZeroMatrix(Dim, Dim);

    sigma *= viscosity;

    const auto& fluid_fraction_gradient = this->GetAtCoordinate(rData.FluidFractionGradient, rData.N);

    double det_permeability = MathUtils<double>::Det(permeability);
    MathUtils<double>::InvertMatrix(permeability, sigma, det_permeability, -1.0);

    sigma *= viscosity;
    Vector grad_alpha_sym_grad_u;
    array_1d<double,Dim> sigma_u;
    for (unsigned int i = 0; i < NumNodes; i++) {
        sigma_u = ZeroVector(Dim);
        BoundedMatrix<double,Dim,Dim> sym_gradient_u = ZeroMatrix(Dim, Dim);
        grad_alpha_sym_grad_u = ZeroVector(Dim);
        for (unsigned int d = 0; d < Dim; d++) {
            double div_u = 0.0;
            for (unsigned int e = 0; e < Dim; e++){
                sigma_u[d] += sigma(d,e) * rData.N[i] * rData.Velocity(i,e);
                sym_gradient_u(d,e) += 1.0 / 2.0 * (rData.DN_DX(i,d) * rData.Velocity(i,e) + rData.DN_DX(i,e) * rData.Velocity(i,d));
                grad_alpha_sym_grad_u[d] += fluid_fraction_gradient[d] * sym_gradient_u(d,e);
                div_u += rData.DN_DX(i,e) * rData.Velocity(i,e);
            }
            rMomentumRHS[d] += density * ( rData.N[i] * (rData.BodyForce(i,d)) /*- fluid_fraction *  rData.N[i] * rAcc[d]*/
                                - fluid_fraction * AGradN[i] * rData.Velocity(i,d)) + 2 * grad_alpha_sym_grad_u[d] * viscosity - 2.0 / 3.0 * viscosity * fluid_fraction_gradient[d] * div_u - fluid_fraction * rData.DN_DX(i,d)*rData.Pressure[i] - sigma_u[d];
        }
    }
}

template< class TElementData >
void AlternativeDVMSDEMCoupled<TElementData>::AddVelocitySystem(
    TElementData& rData,
    MatrixType &rLocalLHS,
    VectorType &rLocalRHS)
{
    auto& LHS = rData.LHS;
    LHS.clear();

    const double density = this->GetAtCoordinate(rData.Density,rData.N);
    const array_1d<double,3> body_force = density * this->GetAtCoordinate(rData.BodyForce,rData.N);

    const array_1d<double,3> convective_velocity = this->FullConvectiveVelocity(rData);
    array_1d<double,3> velocity = this->GetAtCoordinate(rData.Velocity,rData.N);

    array_1d<double,Dim>& r_prev_velocity = mPreviousVelocity[rData.IntegrationPointIndex];
    for (unsigned int n = 0; n < Dim; n++){
        r_prev_velocity[n] = velocity[n];
    }

    BoundedMatrix<double,Dim,Dim> tau_one = ZeroMatrix(Dim, Dim);
    double tau_two;
    this->CalculateStabilizationParameters(rData,convective_velocity,tau_one,tau_two);

    const double dt = rData.DeltaTime;
    const double fluid_fraction = this->GetAtCoordinate(rData.FluidFraction, rData.N);
    // small scale velocity contributions (subscale tracking)
    array_1d<double,Dim> OldUssTerm = fluid_fraction * (density/dt) * mOldSubscaleVelocity[rData.IntegrationPointIndex]; // rho * u_ss^{n-1}/dt

    Vector AGradN;
    this->ConvectionOperator(AGradN,convective_velocity,rData.DN_DX);

    // These two should be zero unless we are using OSS
    const array_1d<double,3> MomentumProj = this->GetAtCoordinate(rData.MomentumProjection,rData.N);
    const double MassProj = this->GetAtCoordinate(rData.MassProjection,rData.N);

    // Multiplying convective operator by density to have correct units
    AGradN *= density;

    double viscosity = this->GetAtCoordinate(rData.DynamicViscosity, rData.N);
    double kin_viscosity = viscosity / density;
    const double fluid_fraction_rate = this->GetAtCoordinate(rData.FluidFractionRate, rData.N);
    const double mass_source = this->GetAtCoordinate(rData.MassSource, rData.N);
    BoundedMatrix<double,Dim,Dim> permeability = this->GetAtCoordinate(rData.Permeability, rData.N);
    BoundedMatrix<double,Dim,Dim> sigma = ZeroMatrix(Dim, Dim);
    array_1d<double, 3> fluid_fraction_gradient = this->GetAtCoordinate(rData.FluidFractionGradient, rData.N);

    double det_permeability = MathUtils<double>::Det(permeability);
    MathUtils<double>::InvertMatrix(permeability, sigma, det_permeability, -1.0);

    sigma *= viscosity;
    // Multiplying convective operator by density to have correct units
    // Note: Dof order is (u,v,[w,]p) for each node
    for (unsigned int i = 0; i < NumNodes; i++) {

        unsigned int row = i*BlockSize;

        // LHS terms
        for (unsigned int j = 0; j < NumNodes; j++) {
            unsigned int col = j*BlockSize;

            // Some terms are the same for all velocity components, calculate them once for each i,j

            // Skew-symmetric convective term 1/2( v*grad(u)*u - grad(v) uu )
            //double K = 0.5*(rN[i]*AGradN[j] - AGradN[i]*rN[j]);
            double V = fluid_fraction * rData.N[i] * AGradN[j];

            // q-p stabilization block (reset result)
            double G = 0.0;
            for (unsigned int d = 0; d < Dim; d++) {
                // Stabilization: u*grad(v) * TauOne * u*grad(u) - vh * TauOne/Dt u*grad(u)
                // The last term comes from vh*d(u_ss)/dt
                double AA = tau_one(d,d) * AGradN[i] * std::pow(fluid_fraction, 2) * AGradN[j];

                double A = tau_one(d,d) * density * std::pow(fluid_fraction, 2) * rData.N[i]/dt * AGradN[j];

                double GBetaDiag = tau_one(d,d) * fluid_fraction * kin_viscosity * rData.N[i] / dt * rData.DN_DX(j,d) * fluid_fraction_gradient[d];
                double AGBetaDiag = tau_one(d,d) * fluid_fraction * kin_viscosity * AGradN[i] * (rData.DN_DX(j,d) * fluid_fraction_gradient[d]);
                double GBetaADiag = 2.0 / 3.0 * kin_viscosity * fluid_fraction * tau_one(d,d) * AGradN[j] * fluid_fraction_gradient[d] * rData.DN_DX(i,d);

                LHS(row+d,col+d) += rData.Weight * (V + AA - A + GBetaDiag - AGBetaDiag + GBetaADiag);
                // Galerkin pressure term: Div(v) * p
                double P = fluid_fraction * rData.DN_DX(i,d) * rData.N[j];

                double QD = fluid_fraction * rData.DN_DX(j,d) * rData.N[i];

                double GP = fluid_fraction_gradient[d] * rData.N[j] * rData.N[i];

                double GD = fluid_fraction_gradient[d] * rData.N[j] * rData.N[i];

                /* q-p stabilization block */
                // Stabilization: Grad(q) * TauOne * Grad(p)
                G += tau_one(d,d) * fluid_fraction * fluid_fraction * rData.DN_DX(i,d) * rData.DN_DX(j,d);
                /* v-u block */
                // Stabilization: Div(v) * TauTwo * Div(u)
                double GR = 0.0;
                double RSigmaG = 0.0;
                double GG = 0.0;
                double DG = 0.0;
                double GGBeta = 0.0;
                double GDBeta = 0.0;
                double GA = tau_one(d,d) * fluid_fraction * fluid_fraction * rData.DN_DX(i,d) * AGradN[j];
                // Stabilization: (a * Grad(v)) * TauOne * Grad(p)
                double AG = fluid_fraction * fluid_fraction * tau_one(d,d) * AGradN[i] * rData.DN_DX(j,d);
                // From vh*d(u_ss)/dt: vh * TauOne/Dt * Grad(p)
                double VP = tau_one(d,d) * density * rData.N[i] / dt * rData.DN_DX(j,d) * fluid_fraction * fluid_fraction;
                for (unsigned int e = 0; e < Dim; e++){
                    double GBetaGDiag = tau_one(d,d) * std::pow(kin_viscosity, 2) * fluid_fraction_gradient[d] * rData.DN_DX(i,d) * (rData.DN_DX(j,e) * fluid_fraction_gradient[e]);
                    GGBeta += tau_one(d,d) * fluid_fraction * kin_viscosity * (rData.DN_DX(i,e) * rData.DN_DX(j,e) * fluid_fraction_gradient[d]);
                    GGBeta += tau_one(d,d) * fluid_fraction * kin_viscosity * (rData.DN_DX(i,d) * fluid_fraction_gradient[e] * rData.DN_DX(j,e));
                    GDBeta += 2.0 / 3.0 * tau_one(d,d) * fluid_fraction * kin_viscosity * fluid_fraction_gradient[e] * rData.DN_DX(i,e) * rData.DN_DX(j,d);
                    double RSigma = rData.N[i] * sigma(d,e) * rData.N[j];
                    double VSigma = tau_one(d,d) * fluid_fraction * density/dt * rData.N[i] * rData.N[j] * sigma(d,e);
                    double ASigma = tau_one(d,d) * AGradN[i] * sigma(d,e) * rData.N[j] * fluid_fraction;
                    double ADBeta = 2.0 / 3.0 * tau_one(d,d) * fluid_fraction * kin_viscosity * AGradN[i] * fluid_fraction_gradient[d] * rData.DN_DX(j,e);
                    double RRSigma = tau_one(d,d) * sigma(d,e) * rData.N[i] * sigma(e,d) * rData.N[j];
                    double RSigmaA = tau_one(d,d) * fluid_fraction * sigma(d,e) * rData.N[i] * AGradN[j];
                    double DD = tau_two * fluid_fraction * fluid_fraction * rData.DN_DX(i,d) * rData.DN_DX(j,e);
                    double DU = tau_two * fluid_fraction * fluid_fraction_gradient[e] * rData.DN_DX(i,d) * rData.N[j];
                    double DBeta = 2.0 / 3.0 * tau_one(d,d) * fluid_fraction * kin_viscosity * rData.N[i] / dt * fluid_fraction_gradient[d] * rData.DN_DX(j,e);
                    double DBetaA = 2.0 /3.0 * kin_viscosity * fluid_fraction * tau_one(d,d) * AGradN[j] * fluid_fraction_gradient[e] * rData.DN_DX(i,d);
                    double GU = tau_two * rData.N[j] * rData.N[i] * fluid_fraction_gradient[d] * fluid_fraction_gradient[e];
                    double GD = tau_two * fluid_fraction * rData.N[i] * fluid_fraction_gradient[d] * rData.DN_DX(j,e);
                    double GBeta = tau_one(d,d) * fluid_fraction * kin_viscosity * rData.N[i] / dt * (rData.DN_DX(j,d) * fluid_fraction_gradient[e]);
                    double AGBeta = tau_one(d,d) * fluid_fraction * kin_viscosity * AGradN[i] * (rData.DN_DX(j,d) * fluid_fraction_gradient[e]);
                    double GBetaA = tau_one(d,d) * fluid_fraction * kin_viscosity * AGradN[j] * rData.DN_DX(i,e) * fluid_fraction_gradient[d];
                    double GBetaG = 0.0;
                    double GBetaD = 0.0;
                    double GBetaSigma = 0.0;
                    double DBetaG = 0.0;
                    double DBetaD = 0.0;
                    double DBetaSigma = 0.0;
                    double RGBeta = 0.0;
                    double RDBeta = 0.0;
                    GG += tau_one(d,d) * fluid_fraction * kin_viscosity * (fluid_fraction_gradient[d] * rData.DN_DX(i,e) * rData.DN_DX(j,e) + fluid_fraction_gradient[e] * rData.DN_DX(i,e) * rData.DN_DX(j,d));
                    DG += 2.0 / 3.0 * tau_one(d,d) * fluid_fraction * kin_viscosity * rData.DN_DX(i,d) * fluid_fraction_gradient[e] * rData.DN_DX(j,e);
                    for (unsigned int f = 0; f < Dim; f++){
                        GBetaG += tau_one(d,d) * std::pow(kin_viscosity, 2) * fluid_fraction_gradient[d] * (rData.DN_DX(i,f) * rData.DN_DX(j,f) * fluid_fraction_gradient[e]);
                        GBetaG += tau_one(d,d) * std::pow(kin_viscosity, 2) * fluid_fraction_gradient[d] * (rData.DN_DX(j,f) * rData.DN_DX(i,e) * fluid_fraction_gradient[f]);
                        GBetaG += tau_one(d,d) * std::pow(kin_viscosity, 2) * fluid_fraction_gradient[f] * rData.DN_DX(i,f) * (rData.DN_DX(j,d) * fluid_fraction_gradient[e]);
                        GBetaD += 2.0 * (2.0 / 3.0 * tau_one(d,d) * std::pow(kin_viscosity, 2) * (fluid_fraction_gradient[f] * rData.DN_DX(i,f)) * fluid_fraction_gradient[d] * rData.DN_DX(j,e));
                        DBetaG += 4.0 / 3.0 * tau_one(d,d) * std::pow(kin_viscosity, 2) * rData.DN_DX(i,d) * fluid_fraction_gradient[f] * rData.DN_DX(j,f) * fluid_fraction_gradient[e];
                        DBetaD += 4.0 / 9.0 * tau_one(d,d) * std::pow(kin_viscosity, 2) * fluid_fraction_gradient[f] * fluid_fraction_gradient[f] * rData.DN_DX(i,d) * rData.DN_DX(j,e);
                        DBetaSigma += 2.0 / 3.0 * tau_one(d,d) * kin_viscosity * rData.DN_DX(i,d) * rData.N[j] * fluid_fraction_gradient[f] * sigma(f,e);
                        RDBeta += 2.0 / 3.0 * tau_one(d,d) * kin_viscosity * rData.N[i] * sigma(d,f) * fluid_fraction_gradient[f] * rData.DN_DX(j,e);
                        for (unsigned int g = 0; g < Dim; g++){
                            GBetaSigma += tau_one(d,d) * kin_viscosity * (fluid_fraction_gradient[g] * rData.DN_DX(i,f) * sigma(f,e) + fluid_fraction_gradient[g] * rData.DN_DX(i,g) * sigma(d,f));
                            for(unsigned int h = 0; h < Dim; h++){
                                RGBeta += tau_one(d,d) * kin_viscosity * (rData.N[d] * sigma(d,f) * rData.DN_DX(i,f) * fluid_fraction_gradient[e] + rData.N[d] * fluid_fraction_gradient[h] * rData.DN_DX(i,h) * sigma(d,e));
                                }
                        }
                    }

                    GR += tau_one(d,d) * fluid_fraction * rData.DN_DX(i,d) * sigma(d,e) * rData.N[j];
                    RSigmaG += tau_one(d,d) * sigma(d,e) * rData.N[i] * rData.DN_DX(j,e);
                    LHS(row+d,col+e) += rData.Weight * (GBeta - DBeta + DD + DU + GU + GD + GBetaA - GBetaG + GBetaD + GBetaSigma - DBetaA + DBetaG - DBetaD - AGBeta + ADBeta - DBetaSigma + RGBeta - RDBeta  + RSigma - VSigma + ASigma - RRSigma - RSigmaA);
                    if (d == e){
                        LHS(row+d,col+e) -= rData.Weight * (GBetaGDiag);
                    }
                }

                LHS(row+Dim,col+d) += rData.Weight * (GA + QD + GR - GGBeta + GDBeta + GD);

                LHS(row+d,col+Dim) += rData.Weight * (AG - VP - P - GP - RSigmaG + GG - DG);

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
            double VI = tau_one(d,d) * density * fluid_fraction * rData.N[i] / dt * (body_force[d] - MomentumProj[d] + OldUssTerm[d]);
            double AF = tau_one(d,d) * fluid_fraction * AGradN[i] * (body_force[d] - MomentumProj[d] + OldUssTerm[d]);
            double RSigmaF = 0.0;
            double GBetaF = 0.0;
            double DBetaF = 0.0;
            for (unsigned int e = 0; e < Dim; ++e){
                RSigmaF += tau_one(d,d) * sigma(d,e) * rData.N[i] * (body_force[e] - MomentumProj[e] + OldUssTerm[e]);
                DBetaF += 2.0 / 3.0 * tau_one(d,d) * kin_viscosity * rData.DN_DX(i,d) * fluid_fraction_gradient[e] * (body_force[e] - MomentumProj[e] + OldUssTerm[e]);
                GBetaF += (tau_one(d,d) * kin_viscosity * (fluid_fraction_gradient[d] * rData.DN_DX(i,e)) * (body_force[e] - MomentumProj[e] + OldUssTerm[e]));
                // GBetaF += tau_one(d,d) * kin_viscosity * (fluid_fraction_gradient[d] * rData.DN_DX(i,d)) * (body_force[e] - MomentumProj[e] + OldUssTerm[e]);
                for (unsigned int f = 0; f < Dim; ++f){
                    GBetaF += tau_one(d,d) * kin_viscosity * (fluid_fraction_gradient[f] * rData.DN_DX(i,f)) * (body_force[e] - MomentumProj[e] + OldUssTerm[e]);
                }
            }
            // Grad(q) * TauOne * (Density * BodyForce - Projection)
            QAlphaF += tau_one(d,d) * rData.DN_DX(i,d) * fluid_fraction * (body_force[d] - MomentumProj[d] + OldUssTerm[d]);
            double VPhi = tau_two * rData.N[i] * fluid_fraction_gradient[d] * (mass_source - fluid_fraction_rate - MassProj);
            // OSS pressure subscale projection
            double DPhi = rData.DN_DX(i,d) * tau_two * fluid_fraction * (mass_source - fluid_fraction_rate - MassProj);
            rLocalRHS[row+d] += rData.Weight * (VF - VI + AF + DPhi - RSigmaF + GBetaF - DBetaF + VPhi);
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

template< class TElementData >
void AlternativeDVMSDEMCoupled<TElementData>::AddMassLHS(
    TElementData& rData,
    MatrixType &rMassMatrix)
{
    const double density = this->GetAtCoordinate(rData.Density,rData.N);
    const double fluid_fraction = this->GetAtCoordinate(rData.FluidFraction, rData.N);
    // Note: Dof order is (u,v,[w,]p) for each node
    for (unsigned int i = 0; i < NumNodes; i++)
    {
        unsigned int row = i*BlockSize;
        for (unsigned int j = 0; j < NumNodes; j++)
        {
            unsigned int col = j*BlockSize;
            const double Mij = rData.Weight * density * fluid_fraction * rData.N[i] * rData.N[j];
            for (unsigned int d = 0; d < Dim; d++)
                rMassMatrix(row+d,col+d) += Mij;
        }
    }

    /* Note on OSS and full projection: Riccardo says that adding the terms provided by
     * AddMassStabilization (and incluiding their corresponding terms in the projeciton)
     * could help reduce the non-linearity of the coupling between projection and u,p
     * However, leaving them on gives a lot of trouble whith the Bossak scheme:
     * think that we solve F - (1-alpha)*M*u^(n+1) - alpha*M*u^(n) - K(u^(n+1)) = 0
     * so the projection of the dynamic terms should be Pi( (1-alpha)*u^(n+1) - alpha*u^(n) )
     */
    if ( !rData.UseOSS)
        this->AddMassStabilization(rData,rMassMatrix);
}

template<class TElementData>
void AlternativeDVMSDEMCoupled<TElementData>::MassProjTerm(
    const TElementData& rData,
    double &rMassRHS) const
{
        const auto velocities = rData.Velocity;

        const double fluid_fraction = this->GetAtCoordinate(rData.FluidFraction, rData.N);
        const double mass_source = this->GetAtCoordinate(rData.MassSource, rData.N);
        const double fluid_fraction_rate = this->GetAtCoordinate(rData.FluidFractionRate, rData.N);
        const auto fluid_fraction_gradient = this->GetAtCoordinate(rData.FluidFractionGradient, rData.N);

        // Compute this node's contribution to the residual (evaluated at integration point)
        for (unsigned int i = 0; i < NumNodes; i++) {
            for (unsigned int d = 0; d < Dim; ++d)
            {
                rMassRHS -= fluid_fraction * rData.DN_DX(i, d) * velocities(i, d) + fluid_fraction_gradient[d] * rData.N[i] * velocities(i,d);
            }
            rMassRHS += mass_source - fluid_fraction_rate;
        }
}

///////////////////////////////////////////////////////////////////////////////////////////////////

template< class TElementData >
void AlternativeDVMSDEMCoupled<TElementData>::AddMassStabilization(
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
    array_1d<double, 3> fluid_fraction_gradient = this->GetAtCoordinate(rData.FluidFractionGradient, rData.N);
    double viscosity = this->GetAtCoordinate(rData.DynamicViscosity, rData.N);
    double kin_viscosity = viscosity / density;
    BoundedMatrix<double,Dim,Dim> permeability = this->GetAtCoordinate(rData.Permeability, rData.N);
    BoundedMatrix<double,Dim,Dim> sigma = ZeroMatrix(Dim, Dim);

    double det_permeability = MathUtils<double>::Det(permeability);
    MathUtils<double>::InvertMatrix(permeability, sigma, det_permeability, -1.0);

    double W = rData.Weight * density; // This density is for the dynamic term in the residual (rho*Du/Dt)
    sigma *= viscosity;

    // Note: Dof order is (u,v,[w,]p) for each node
    for (unsigned int i = 0; i < NumNodes; i++) {
        unsigned int row = i*BlockSize;

        for (unsigned int j = 0; j < NumNodes; j++) {
            unsigned int col = j*BlockSize;

            for (unsigned int d = 0; d < Dim; d++)
            {
                // grad(q) * TauOne * du/dt
                double UGAlpha = tau_one(d,d) * fluid_fraction * fluid_fraction * rData.DN_DX(i,d) * rData.N[j];
                // u*grad(v) * TauOne * du/dt
                // v * TauOne/dt * du/dt (from v*d(uss)/dt)
                double AU = tau_one(d,d) * fluid_fraction * fluid_fraction * AGradN[i] * rData.N[j];
                double IU = tau_one(d,d) * std::pow(fluid_fraction, 2) * rData.N[i]/dt * rData.N[j];
                double GBetaUDiag = tau_one(d,d) * fluid_fraction * kin_viscosity * (fluid_fraction_gradient[d] * rData.DN_DX(i,d)) * rData.N[j];
                for (unsigned int e = 0; e < Dim; ++e){
                    double RSigmaU = tau_one(d,d) * sigma(d,e) * rData.N[i] * rData.N[j];
                    double DBetaU = 2.0 / 3.0 * tau_one(d,d) * fluid_fraction * kin_viscosity * rData.DN_DX(i,d) * rData.N[j] * fluid_fraction_gradient[e];
                    double GBetaU = tau_one(d,d) * fluid_fraction * kin_viscosity * (fluid_fraction_gradient[d] * rData.DN_DX(i,e)) * rData.N[j];

                    rMassMatrix(row+d, col+e) += W * (GBetaU - RSigmaU - DBetaU);
                }
                rMassMatrix(row+d,col+d) += W * (AU - IU + GBetaUDiag);
                rMassMatrix(row+Dim,col+d) += W * UGAlpha;
            }
        }
    }
}

template <class TElementData>
void AlternativeDVMSDEMCoupled<TElementData>::AddViscousTerm(
    const TElementData& rData,
    BoundedMatrix<double,LocalSize,LocalSize>& rLHS,
    VectorType& rRHS) {

    const double fluid_fraction = this->GetAtCoordinate(rData.FluidFraction, rData.N);
    BoundedMatrix<double,StrainSize,LocalSize> strain_matrix = ZeroMatrix(StrainSize,LocalSize);
    FluidElementUtilities<NumNodes>::GetStrainMatrix(rData.DN_DX,strain_matrix);

    const auto& constitutive_matrix = rData.C;
    BoundedMatrix<double,StrainSize,LocalSize> shear_stress_matrix = prod(constitutive_matrix,strain_matrix);

    // Multiply times integration point weight (I do this here to avoid a temporal in LHS += weight * Bt * C * B)
    strain_matrix *= rData.Weight;

    noalias(rLHS) += prod(trans(strain_matrix), fluid_fraction * shear_stress_matrix);
    noalias(rRHS) -= prod(trans(strain_matrix), fluid_fraction * rData.ShearStress);
}

template< class TElementData >
void AlternativeDVMSDEMCoupled<TElementData>::CalculateStabilizationParameters(
    const TElementData& rData,
    const array_1d<double,3> &Velocity,
    BoundedMatrix<double,Dim,Dim> &TauOne,
    double &TauTwo) const
{
    double inv_tau;
    double inv_tau_NS;
    const double h = rData.ElementSize;
    const double density = this->GetAtCoordinate(rData.Density,rData.N);
    const double viscosity = this->GetAtCoordinate(rData.EffectiveViscosity,rData.N);
    double fluid_fraction = this->GetAtCoordinate(rData.FluidFraction, rData.N);
    constexpr double c1 = DVMS<TElementData>::mTauC1;
    constexpr double c2 = DVMS<TElementData>::mTauC2;
    BoundedMatrix<double,Dim,Dim> permeability = this->GetAtCoordinate(rData.Permeability, rData.N);
    BoundedMatrix<double,Dim,Dim> sigma = ZeroMatrix(Dim, Dim);
    BoundedMatrix<double,Dim,Dim> I = IdentityMatrix(Dim, Dim);
    array_1d<double, 3> fluid_fraction_gradient = this->GetAtCoordinate(rData.FluidFractionGradient, rData.N);

    double det_permeability = MathUtils<double>::Det(permeability);
    MathUtils<double>::InvertMatrix(permeability, sigma, det_permeability, -1.0);

    double velocity_modulus = 0.0;
    double fluid_fraction_gradient_modulus = 0.0;
    double sigma_term = 0.0;
    for (unsigned int d = 0; d < Dim; d++){
        velocity_modulus += Velocity[d] * Velocity[d];
        fluid_fraction_gradient_modulus += std::pow(fluid_fraction_gradient[d],2);
        for (unsigned int e = d; e < Dim; e++){
            sigma_term += std::pow(sigma(d,e),2);
        }
    }

    double velocity_norm = std::sqrt(velocity_modulus);
    double fluid_fraction_gradient_norm = std::sqrt(fluid_fraction_gradient_modulus);
    double c_alpha = fluid_fraction + h / c1 * fluid_fraction_gradient_norm;
    inv_tau = (c1 * viscosity / (h*h) + density * (c2 * velocity_norm / h ) ) * c_alpha + density * fluid_fraction / rData.DeltaTime + std::sqrt(sigma_term);
    inv_tau_NS = c1 * viscosity / (h*h) + density * (c2 * velocity_norm / h ) + std::sqrt(sigma_term);
    double tau_one = 1 / inv_tau;
    double tau_one_NS = 1 / inv_tau_NS;

    TauOne = tau_one * I;
    TauTwo = h * h / (c1 * fluid_fraction * tau_one_NS);
}

template< class TElementData >
void AlternativeDVMSDEMCoupled<TElementData>::SubscaleVelocity(
    const TElementData& rData,
    array_1d<double,3>& rVelocitySubscale) const
{
    const double density = this->GetAtCoordinate(rData.Density,rData.N);
    array_1d<double,3> convective_velocity = this->FullConvectiveVelocity(rData);
    BoundedMatrix<double,Dim,Dim> tau_one = ZeroMatrix(Dim,Dim);
    double fluid_fraction = this->GetAtCoordinate(rData.FluidFraction, rData.N);
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
        rVelocitySubscale[d] = tau_one(d,d)*(residual[d] + fluid_fraction * (density/dt)*rOldSubscaleVelocity[d]);
    }
}

template< class TElementData >
void AlternativeDVMSDEMCoupled<TElementData>::SubscalePressure(
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
array_1d<double,3> AlternativeDVMSDEMCoupled<TElementData>::FullConvectiveVelocity(
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
void AlternativeDVMSDEMCoupled<TElementData>::UpdateSubscaleVelocity(
    const TElementData& rData)
{
    const double density = this->GetAtCoordinate(rData.Density,rData.N);
    double fluid_fraction = this->GetAtCoordinate(rData.FluidFraction,rData.N);
    array_1d<double,3> resolved_convection_velocity = this->FullConvectiveVelocity(rData);
    const double dt = rData.DeltaTime;
    array_1d<double,3> predicted_subscale_velocity = ZeroVector(3);

    const array_1d<double,Dim>& r_old_subscale_velocity = mOldSubscaleVelocity[rData.IntegrationPointIndex];

    array_1d<double, 3> previous_subscale_velocity = ZeroVector(3);

    array_1d<double, 3> previous_velocity = mPreviousVelocity[rData.IntegrationPointIndex];

    for (size_t i = 0; i < NumNodes; i++) {
        array_1d<double, 3> subscale_velocity_on_previous_iteration = mPredictedSubscaleVelocity[rData.IntegrationPointIndex];
        for (size_t d = 0; d < Dim; d++) {
            previous_subscale_velocity[d] += rData.N[i] * subscale_velocity_on_previous_iteration[d];
        }
    }

    // Part of the residual that does not depend on the subscale
    array_1d<double,3> static_residual = ZeroVector(3);

    // Note I'm only using large scale convection here, small-scale convection is re-evaluated at each iteration.
    if (!rData.UseOSS)
        this->AlgebraicMomentumResidual(rData,resolved_convection_velocity,static_residual);
    else
        this->OrthogonalMomentumResidual(rData,resolved_convection_velocity,static_residual);

    array_1d<double,3> v_d = ZeroVector(Dim);
    for (unsigned int d = 0; d < Dim; d++)
    {
        v_d[d] = previous_velocity[d] + previous_subscale_velocity[d];
    }
    BoundedMatrix<double,Dim,Dim> tau_one = ZeroMatrix(Dim, Dim);
    double tau_two;

    this->CalculateStabilizationParameters(rData,v_d,tau_one,tau_two);

    for (unsigned int d = 0; d < Dim; d++)
    {
        predicted_subscale_velocity[d] = tau_one(d,d) * (static_residual[d] + fluid_fraction * (density/dt)*r_old_subscale_velocity[d]);
    }
    noalias(mPredictedSubscaleVelocity[rData.IntegrationPointIndex]) = predicted_subscale_velocity;

}

template< class TElementData >
void AlternativeDVMSDEMCoupled<TElementData>::UpdateSubscaleVelocityPrediction(
    const TElementData& rData)
{
    const double density = this->GetAtCoordinate(rData.Density,rData.N);
    const double viscosity = this->GetAtCoordinate(rData.EffectiveViscosity,rData.N);
    double fluid_fraction = this->GetAtCoordinate(rData.FluidFraction,rData.N);
    array_1d<double,3> resolved_convection_velocity = this->GetAtCoordinate(rData.Velocity,rData.N) - this->GetAtCoordinate(rData.MeshVelocity,rData.N);
    const double dt = rData.DeltaTime;
    const double h = rData.ElementSize;

    // Elemental large-scale velocity gradient
    BoundedMatrix<double,Dim,Dim> resolved_velocity_gradient = ZeroMatrix(Dim,Dim);

    BoundedMatrix<double,Dim,Dim> sigma = ZeroMatrix(Dim, Dim);
    BoundedMatrix<double,Dim,Dim> permeability = this->GetAtCoordinate(rData.Permeability, rData.N);
    double det_permeability = MathUtils<double>::Det(permeability);
    MathUtils<double>::InvertMatrix(permeability, sigma, det_permeability, -1.0);
    array_1d<double, Dim> fluid_fraction_gradient = this->GetAtCoordinate(rData.FluidFractionGradient, rData.N);

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

    double sigma_term = 0.0;
    double fluid_fraction_gradient_norm_squared = 0.0;

    // Add the time discretization term to obtain the part of the residual that does not change during iteration
    for (unsigned int d = 0; d < Dim; d++){
        fluid_fraction_gradient_norm_squared += std::pow(fluid_fraction_gradient[d],2);
        static_residual[d] += fluid_fraction * density/dt * old_subscale_velocity[d];
        for (unsigned int e = d; e < Dim; e++){
                sigma_term += std::pow(sigma(d,e),2);
        }
    }

    const double fluid_fraction_gradient_norm = std::sqrt(fluid_fraction_gradient_norm_squared);

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
    array_1d<double,Dim>& u = mPredictedSubscaleVelocity[rData.IntegrationPointIndex]; // Use last result as initial guess
    array_1d<double,Dim> du = ZeroVector(Dim);

    while ( (!converged) && (iter++ < subscale_prediction_maximum_iterations) ) {


        // Calculate new Tau
        double convection_velocity_norm_squared = 0.0;
        Vector v_d = ZeroVector(Dim);
        for (unsigned int d = 0; d < Dim; d++)
        {
            v_d[d] = resolved_convection_velocity[d] + u[d];
            convection_velocity_norm_squared += v_d[d] * v_d[d];
        }
        double convection_velocity_norm = std::sqrt(convection_velocity_norm_squared);

        double c_alpha = fluid_fraction + h / c1 * fluid_fraction_gradient_norm;
        double inv_tau_t = (c1 * viscosity / (h * h) + density * (c2 * convection_velocity_norm / h ) + std::sqrt(sigma_term)) * c_alpha + density * fluid_fraction / dt;
        double inv_tau = (c1 * viscosity / (h * h) + density * (c2 * convection_velocity_norm / h ) + std::sqrt(sigma_term)) * c_alpha;
        double tau_one = 1 / inv_tau;

        Vector grad_tau_one = - 9 * c2 * h * v_d / (convection_velocity_norm * c_alpha * std::pow(3 * c2 * convection_velocity_norm + 4 * c1 * std::pow(h, 3) * viscosity + 3 * h * std::sqrt(sigma_term), 2));
        //double c_alpha = 1.0;
        Vector grad_tau_one_t = (dt * grad_tau_one * (tau_one * density * fluid_fraction + dt) - dt * tau_one * (grad_tau_one * density * fluid_fraction)) / std::pow(tau_one * density * fluid_fraction + dt, 2);

        noalias(J) = density * fluid_fraction * resolved_velocity_gradient;
        for (unsigned int d = 0; d < Dim; d++)
        {
            J(d,d) += inv_tau_t;
            for (unsigned int e = 0; e < Dim; e++){
                J(d,e) -= grad_tau_one_t[d] * u[e] * std::pow(inv_tau_t, 2);
            }
        }
        // Newton-Raphson RHSs
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
void AlternativeDVMSDEMCoupled<TElementData>::save(Serializer& rSerializer) const
{
    typedef DVMS<TElementData> BaseElement;
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseElement );
    rSerializer.save("mOldSubscaleVelocity",mOldSubscaleVelocity);
}


template< class TElementData >
void AlternativeDVMSDEMCoupled<TElementData>::load(Serializer& rSerializer)
{
    typedef DVMS<TElementData> BaseElement;
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseElement);
    rSerializer.load("mOldSubscaleVelocity",mOldSubscaleVelocity);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Class template instantiation

template class AlternativeDVMSDEMCoupled< QSVMSDEMCoupledData<2,3> >;
template class AlternativeDVMSDEMCoupled< QSVMSDEMCoupledData<3,4> >;

template class AlternativeDVMSDEMCoupled< QSVMSDEMCoupledData<2,4> >;
template class AlternativeDVMSDEMCoupled< QSVMSDEMCoupledData<3,8> >;

} // namespace Kratos
