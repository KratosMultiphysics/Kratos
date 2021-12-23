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
#include "alternative_qs_vms_dem_coupled.h"
#include "custom_utilities/qsvms_dem_coupled_data.h"
#include "custom_utilities/fluid_element_utilities.h"
#include "fluid_dynamics_application_variables.h"

namespace Kratos
{

//////////////////////////Life cycle

template< class TElementData >
AlternativeQSVMSDEMCoupled<TElementData>::AlternativeQSVMSDEMCoupled(IndexType NewId):
    QSVMS<TElementData>(NewId)
{}

template< class TElementData >
AlternativeQSVMSDEMCoupled<TElementData>::AlternativeQSVMSDEMCoupled(IndexType NewId, const NodesArrayType& ThisNodes):
    QSVMS<TElementData>(NewId,ThisNodes)
{}


template< class TElementData >
AlternativeQSVMSDEMCoupled<TElementData>::AlternativeQSVMSDEMCoupled(IndexType NewId, GeometryType::Pointer pGeometry):
    QSVMS<TElementData>(NewId,pGeometry)
{}


template< class TElementData >
AlternativeQSVMSDEMCoupled<TElementData>::AlternativeQSVMSDEMCoupled(IndexType NewId, GeometryType::Pointer pGeometry, Properties::Pointer pProperties):
    QSVMS<TElementData>(NewId,pGeometry,pProperties)
{}

///////////Destructor

template< class TElementData >
AlternativeQSVMSDEMCoupled<TElementData>::~AlternativeQSVMSDEMCoupled()
{}

template< class TElementData >
Element::Pointer AlternativeQSVMSDEMCoupled<TElementData>::Create(IndexType NewId,NodesArrayType const& ThisNodes,Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<AlternativeQSVMSDEMCoupled>(NewId, this->GetGeometry().Create(ThisNodes), pProperties);
}

template< class TElementData >
Element::Pointer AlternativeQSVMSDEMCoupled<TElementData>::Create(IndexType NewId,GeometryType::Pointer pGeom,Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<AlternativeQSVMSDEMCoupled>(NewId, pGeom, pProperties);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Inquiry

template< class TElementData >
int AlternativeQSVMSDEMCoupled<TElementData>::Check(const ProcessInfo &rCurrentProcessInfo) const
{
    int out = QSVMS<TElementData>::Check(rCurrentProcessInfo);
    KRATOS_ERROR_IF_NOT(out == 0)
        << "Error in base class Check for Element " << this->Info() << std::endl
        << "Error code is " << out << std::endl;

    const auto &r_geom = this->GetGeometry();
    for (unsigned int i = 0; i < NumNodes; ++i)
    {
        const auto& rNode = r_geom[i];
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(ACCELERATION,rNode);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(NODAL_AREA,rNode);
    }

    return out;
}
///////////////////////////////////////////////////////////////////////////////////////////////////

template< class TElementData >
std::string AlternativeQSVMSDEMCoupled<TElementData>::Info() const
{
    std::stringstream buffer;
    buffer << "AlternativeQSVMSDEMCoupled #" << this->Id();
    return buffer.str();
}


template< class TElementData >
void AlternativeQSVMSDEMCoupled<TElementData>::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "AlternativeQSVMSDEMCoupled" << Dim << "D";
}

template< class TElementData >
void AlternativeQSVMSDEMCoupled<TElementData>::AlgebraicMomentumResidual(
    const TElementData& rData,
    const array_1d<double,3> &rConvectionVelocity,
    array_1d<double,3>& rResidual) const
{
    const auto& rGeom = this->GetGeometry();

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
    const auto& fluid_fraction_gradient = this->GetAtCoordinate(rData.FluidFractionGradient, rData.N);

    double det_permeability = MathUtils<double>::Det(permeability);
    MathUtils<double>::InvertMatrix(permeability, sigma, det_permeability, -1.0);
    sigma *= viscosity;
    Vector grad_alpha_sym_grad_u, sigma_U;
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
                    sym_gradient_u(d,e) += 1.0/2.0 * (rData.DN_DX(i,d) * r_velocities(i,e) + rData.DN_DX(i,e) * r_velocities(i,d));
                    grad_alpha_sym_grad_u[d] += fluid_fraction_gradient[d] * sym_gradient_u(d,e);
                    div_u += rData.DN_DX(i,e) * r_velocities(i,e);
                }
                rResidual[d] += density * (rData.N[i] * r_body_forces(i,d) - fluid_fraction * rData.N[i] * r_acceleration[d] - fluid_fraction * convection[i] * r_velocities(i,d)) + 2 * grad_alpha_sym_grad_u[d] * viscosity - 2.0 / 3.0 * viscosity * fluid_fraction_gradient[d] * div_u - fluid_fraction * rData.DN_DX(i,d) * r_pressures[i] - sigma_U[d];
            }
        }
}

template< class TElementData >
void AlternativeQSVMSDEMCoupled<TElementData>::MomentumProjTerm(
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
        sigma_U = ZeroVector(Dim);
        sym_gradient_u = ZeroMatrix(Dim, Dim);
        grad_alpha_sym_grad_u = ZeroVector(Dim);
        for (unsigned int d = 0; d < Dim; d++) {
            double div_u = 0.0;
            for (unsigned int e = 0; e < Dim; e++){
                sigma_U[d] += sigma(d,e) * rData.N[i] * r_velocities(i,e);
                sym_gradient_u(d,e) += 1.0/2.0 * (rData.DN_DX(i,d) * r_velocities(i,e) + rData.DN_DX(i,e) * r_velocities(i,d));
                grad_alpha_sym_grad_u[d] += r_fluid_fraction_gradient[d] * sym_gradient_u(d,e);
                div_u += rData.DN_DX(i,e) * r_velocities(i,e);
            }
            rMomentumRHS[d] += density * (rData.N[i] * r_body_forces(i,d) /*- fluid_fraction * rData.N[i] * r_acceleration[d]*/ - fluid_fraction * AGradN[i] * rData.Velocity(i,d)) + 2 * grad_alpha_sym_grad_u[d] * viscosity - 2.0/3.0 * viscosity * r_fluid_fraction_gradient[d] * div_u - fluid_fraction * rData.DN_DX(i,d) * r_pressures[i] - sigma_U[d];
        }
    }
}

template<class TElementData>
void AlternativeQSVMSDEMCoupled<TElementData>::AddMassStabilization(
    TElementData& rData,
    MatrixType& rMassMatrix)
{
    const double density = this->GetAtCoordinate(rData.Density,rData.N);
    const array_1d<double, 3> convective_velocity=
        this->GetAtCoordinate(rData.Velocity, rData.N) -
        this->GetAtCoordinate(rData.MeshVelocity, rData.N);

    BoundedMatrix<double,Dim,Dim> tau_one = ZeroMatrix(Dim, Dim);
    double tau_two;
    this->CalculateTau(rData,convective_velocity,tau_one,tau_two);

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

    double W = rData.Weight * density; // This density is for the dynamic term in the residual (rho*Du)
    sigma *= viscosity;

    // Note: Dof order is (u,v,[w,]p) for each node
    for (unsigned int i = 0; i < NumNodes; i++) {
        unsigned int row = i*BlockSize;

        for (unsigned int j = 0; j < NumNodes; j++) {
            unsigned int col = j*BlockSize;

            for (unsigned int d = 0; d < Dim; d++)
            {
                // grad(q) * TauOne * du
                double UGAlpha = tau_one(d,d) * fluid_fraction * fluid_fraction * rData.DN_DX(i,d) * rData.N[j];
                // u*grad(v) * TauOne * du
                double AU = tau_one(d,d) * fluid_fraction * fluid_fraction * AGradN[i] * rData.N[j];
                double GBetaUDiag = tau_one(d,d) * fluid_fraction * kin_viscosity * (fluid_fraction_gradient[d] * rData.DN_DX(i,d)) * rData.N[j];
                for (unsigned int e = 0; e < Dim; ++e){
                    double RSigmaU = tau_one(d,d) * sigma(d,e) * rData.N[i] * rData.N[j];
                    double DBetaU = 2.0 / 3.0 * tau_one(d,d) * fluid_fraction * kin_viscosity * rData.DN_DX(i,d) * rData.N[j] * fluid_fraction_gradient[e];
                    double GBetaU = tau_one(d,d) * fluid_fraction * kin_viscosity * (fluid_fraction_gradient[d] * rData.DN_DX(i,e)) * rData.N[j];

                    rMassMatrix(row+d, col+e) += W * (GBetaU - RSigmaU - DBetaU);
                }
                rMassMatrix(row+d,col+d) += W * (AU + GBetaUDiag);
                rMassMatrix(row+Dim,col+d) += W * UGAlpha;
            }
        }
    }
}

template <class TElementData>
void AlternativeQSVMSDEMCoupled<TElementData>::AddViscousTerm(
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

// Add a the contribution from a single integration point to the velocity contribution
template< class TElementData >
void AlternativeQSVMSDEMCoupled<TElementData>::AddVelocitySystem(
    TElementData& rData,
    MatrixType& rLocalLHS,
    VectorType& rLocalRHS)
{
    auto& LHS = rData.LHS;
    LHS.clear();

    const double density = this->GetAtCoordinate(rData.Density,rData.N);
    array_1d<double,3> body_force = density * this->GetAtCoordinate(rData.BodyForce,rData.N);

    const array_1d<double, 3> convective_velocity =
        this->GetAtCoordinate(rData.Velocity, rData.N) -
        this->GetAtCoordinate(rData.MeshVelocity, rData.N);

    BoundedMatrix<double,Dim,Dim> tau_one = ZeroMatrix(Dim, Dim);
    double tau_two;
    this->CalculateTau(rData,convective_velocity,tau_one,tau_two);

    Vector AGradN;
    this->ConvectionOperator(AGradN,convective_velocity,rData.DN_DX);

    // These two should be zero unless we are using OSS
    const array_1d<double,3> MomentumProj = this->GetAtCoordinate(rData.MomentumProjection,rData.N);
    const double MassProj = this->GetAtCoordinate(rData.MassProjection,rData.N);

    double viscosity = this->GetAtCoordinate(rData.DynamicViscosity, rData.N);
    double kin_viscosity = viscosity / density;
    const double fluid_fraction = this->GetAtCoordinate(rData.FluidFraction, rData.N);
    const double fluid_fraction_rate = this->GetAtCoordinate(rData.FluidFractionRate, rData.N);
    const double mass_source = this->GetAtCoordinate(rData.MassSource, rData.N);

    BoundedMatrix<double,Dim,Dim> permeability = this->GetAtCoordinate(rData.Permeability, rData.N);
    BoundedMatrix<double,Dim,Dim> sigma = ZeroMatrix(Dim, Dim);
    array_1d<double, 3> fluid_fraction_gradient = this->GetAtCoordinate(rData.FluidFractionGradient, rData.N);

    double det_permeability = MathUtils<double>::Det(permeability);
    MathUtils<double>::InvertMatrix(permeability, sigma, det_permeability, -1.0);

    sigma *= viscosity;
    body_force *= density; // Force per unit of volume
    AGradN *= density; // Convective term is always multiplied by density


    // Multiplying convective operator by density to have correct units
    // Note: Dof order is (u,v,[w,]p) for each node
    for (unsigned int i = 0; i < NumNodes; i++) {

        unsigned int row = i*BlockSize;

        // LHS terms
        for (unsigned int j = 0; j < NumNodes; j++) {
            unsigned int col = j*BlockSize;

            // Some terms are the same for all velocity components, calculate them once for each i,j

            // Skew-symmetric convective term 1/2( v*grad(u)*u - grad(v) uu )
            double V = fluid_fraction * rData.N[i] * AGradN[j];

            // q-p stabilization block (reset result)
            double G = 0.0;
            for (unsigned int d = 0; d < Dim; d++) {
                // Stabilization: u*grad(v) * TauOne * u*grad(u)
                // The last term comes from vh*d(u_ss)
                double AA = tau_one(d,d) * AGradN[i] * std::pow(fluid_fraction, 2) * AGradN[j];
                double AGBetaDiag = tau_one(d,d) * fluid_fraction * kin_viscosity * AGradN[i] * (rData.DN_DX(j,d) * fluid_fraction_gradient[d]);
                double GBetaADiag = 2.0 / 3.0 * kin_viscosity * fluid_fraction * tau_one(d,d) * AGradN[j] * fluid_fraction_gradient[d] * rData.DN_DX(i,d);

                LHS(row+d,col+d) += rData.Weight * (V + AA - AGBetaDiag + GBetaADiag);
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
                for (unsigned int e = 0; e < Dim; e++){
                    double GBetaGDiag = tau_one(d,d) * std::pow(kin_viscosity, 2) * fluid_fraction_gradient[d] * rData.DN_DX(i,d) * (rData.DN_DX(j,e) * fluid_fraction_gradient[e]);
                    GGBeta += tau_one(d,d) * fluid_fraction * kin_viscosity * (rData.DN_DX(i,e) * rData.DN_DX(j,e) * fluid_fraction_gradient[d]);
                    GGBeta += tau_one(d,d) * fluid_fraction * kin_viscosity * (rData.DN_DX(i,d) * fluid_fraction_gradient[e] * rData.DN_DX(j,e));
                    GDBeta += 2.0 / 3.0 * tau_one(d,d) * fluid_fraction * kin_viscosity * fluid_fraction_gradient[e] * rData.DN_DX(i,e) * rData.DN_DX(j,d);
                    double RSigma = rData.N[i] * sigma(d,e) * rData.N[j];
                    double ASigma = tau_one(d,d) * AGradN[i] * sigma(d,e) * rData.N[j] * fluid_fraction;
                    double ADBeta = 2.0 / 3.0 * tau_one(d,d) * fluid_fraction * kin_viscosity * AGradN[i] * fluid_fraction_gradient[d] * rData.DN_DX(j,e);
                    double RRSigma = tau_one(d,d) * sigma(d,e) * rData.N[i] * sigma(e,d) * rData.N[j];
                    double RSigmaA = tau_one(d,d) * fluid_fraction * sigma(d,e) * rData.N[i] * AGradN[j];
                    double DD = tau_two * fluid_fraction * fluid_fraction * rData.DN_DX(i,d) * rData.DN_DX(j,e);
                    double DU = tau_two * fluid_fraction * fluid_fraction_gradient[e] * rData.DN_DX(i,d) * rData.N[j];
                    double DBetaA = 2.0 / 3.0 * kin_viscosity * fluid_fraction * tau_one(d,d) * AGradN[j] * fluid_fraction_gradient[e] * rData.DN_DX(i,d);
                    double GU = tau_two * rData.N[j] * rData.N[i] * fluid_fraction_gradient[d] * fluid_fraction_gradient[e];
                    double GD = tau_two * fluid_fraction * rData.N[i] * fluid_fraction_gradient[d] * rData.DN_DX(j,e);
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
                    LHS(row+d,col+e) += rData.Weight * (DD + DU + GU + GD + GBetaA - GBetaG + GBetaD + GBetaSigma - DBetaA + DBetaG - DBetaD - AGBeta + ADBeta - DBetaSigma + RGBeta - RDBeta  + RSigma + ASigma - RRSigma - RSigmaA);
                    if (d == e){
                        LHS(row+d,col+e) -= rData.Weight * (GBetaGDiag);
                    }
                }

                LHS(row+Dim,col+d) += rData.Weight * (GA + QD + GR - GGBeta + GDBeta + GD);

                LHS(row+d,col+Dim) += rData.Weight * (AG - P - GP - RSigmaG + GG - DG);

            }

            // Write q-p term
            LHS(row+Dim,col+Dim) += rData.Weight * G;
        }

        // RHS terms
        double QAlphaF = 0.0;
        for (unsigned int d = 0; d < Dim; ++d)
        {
            // v*BodyForce
            double VF = rData.N[i] * (body_force[d]);
            // ( a * Grad(v) ) * TauOne * (Density * BodyForce - Projection)
            double AF = tau_one(d,d) * fluid_fraction * AGradN[i] * (body_force[d] - MomentumProj[d]);
            double RSigmaF = 0.0;
            double GBetaF = 0.0;
            double DBetaF = 0.0;
            for (unsigned int e = 0; e < Dim; ++e){
                RSigmaF += tau_one(d,d) * sigma(d,e) * rData.N[i] * (body_force[e] - MomentumProj[e]);
                DBetaF += 2.0 / 3.0 * tau_one(d,d) * kin_viscosity * rData.DN_DX(i,d) * fluid_fraction_gradient[e] * (body_force[e] - MomentumProj[e]);
                GBetaF += tau_one(d,d) * kin_viscosity * (fluid_fraction_gradient[d] * rData.DN_DX(i,e)) * (body_force[e] - MomentumProj[e]);
                for (unsigned int f = 0; f < Dim; ++f){
                    GBetaF += tau_one(d,d) * kin_viscosity * (fluid_fraction_gradient[f] * rData.DN_DX(i,f)) * (body_force[e] - MomentumProj[e]);
                }
            }
            // Grad(q) * TauOne * (Density * BodyForce - Projection)
            QAlphaF += tau_one(d,d) * rData.DN_DX(i,d) * fluid_fraction * (body_force[d] - MomentumProj[d]);
            double VPhi = tau_two * rData.N[i] * fluid_fraction_gradient[d] * (mass_source - fluid_fraction_rate - MassProj);
            // OSS pressure subscale projection
            double DPhi = rData.DN_DX(i,d) * tau_two * fluid_fraction * (mass_source - fluid_fraction_rate - MassProj);
            rLocalRHS[row+d] += rData.Weight * (VF + AF + DPhi - RSigmaF + GBetaF - DBetaF + VPhi);
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
void AlternativeQSVMSDEMCoupled<TElementData>::AddMassLHS(
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
    if (!rData.UseOSS)
        this->AddMassStabilization(rData,rMassMatrix);
}

template<class TElementData>
void AlternativeQSVMSDEMCoupled<TElementData>::MassProjTerm(
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
                rMassRHS -= fluid_fraction * rData.DN_DX(i, d) * velocities(i, d) + fluid_fraction_gradient[d] * rData.N[i] * velocities(i,d);
            }
        }
        rMassRHS += mass_source - fluid_fraction_rate;
}

template< class TElementData >
void AlternativeQSVMSDEMCoupled<TElementData>::CalculateTau(
    const TElementData& rData,
    const array_1d<double,3> &Velocity,
    BoundedMatrix<double,Dim,Dim> &TauOne,
    double &TauTwo) const
{
    constexpr double c1 = 8.0;
    constexpr double c2 = 2.0;
    double inv_tau;
    double inv_tau_NS;
    const double h = rData.ElementSize;
    const double density = this->GetAtCoordinate(rData.Density,rData.N);
    const double viscosity = this->GetAtCoordinate(rData.EffectiveViscosity,rData.N);
    double fluid_fraction = this->GetAtCoordinate(rData.FluidFraction, rData.N);
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
    inv_tau = (c1 * viscosity / (h*h) + density * (c2 * velocity_norm / h )) * c_alpha + std::sqrt(sigma_term);
    inv_tau_NS = c1 * viscosity / (h*h) + density * (c2 * velocity_norm / h ) + std::sqrt(sigma_term);

    double tau_one = 1 / inv_tau;
    double tau_one_NS = 1 / inv_tau_NS;
    TauOne = tau_one * I;
    TauTwo = h * h / (c1 * fluid_fraction * tau_one_NS);
}

template< class TElementData >
void AlternativeQSVMSDEMCoupled<TElementData>::SubscaleVelocity(
    const TElementData& rData,
    array_1d<double,3> &rVelocitySubscale) const
{
    BoundedMatrix<double,Dim,Dim> tau_one = ZeroMatrix(Dim, Dim);
    double tau_two;
    array_1d<double,3> convective_velocity = this->GetAtCoordinate(rData.Velocity,rData.N) - this->GetAtCoordinate(rData.MeshVelocity,rData.N);
    this->CalculateTau(rData,convective_velocity,tau_one,tau_two);

    array_1d<double,3> Residual = ZeroVector(3);

    if (!rData.UseOSS)
        this->AlgebraicMomentumResidual(rData,convective_velocity,Residual);
    else
        this->OrthogonalMomentumResidual(rData,convective_velocity,Residual);

    for (unsigned int d = 0; d < Dim; ++d)
        rVelocitySubscale[d] = tau_one(d,d) * Residual[d];
}

template< class TElementData >
void AlternativeQSVMSDEMCoupled<TElementData>::SubscalePressure(
        const TElementData& rData,
        double &rPressureSubscale) const
{
    BoundedMatrix<double,Dim,Dim> tau_one = ZeroMatrix(Dim, Dim);
    double tau_two;
    array_1d<double, 3> convective_velocity =
        this->GetAtCoordinate(rData.Velocity, rData.N) -
        this->GetAtCoordinate(rData.MeshVelocity, rData.N);
    this->CalculateTau(rData, convective_velocity, tau_one, tau_two);

    double Residual = 0.0;

    if (!rData.UseOSS)
        this->AlgebraicMassResidual(rData,Residual);
    else
        this->OrthogonalMassResidual(rData,Residual);

    rPressureSubscale = tau_two*Residual;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Private functions
///////////////////////////////////////////////////////////////////////////////////////////////////

// serializer

template< class TElementData >
void AlternativeQSVMSDEMCoupled<TElementData>::save(Serializer& rSerializer) const
{
    typedef QSVMS<TElementData> BaseElement;
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseElement );
}


template< class TElementData >
void AlternativeQSVMSDEMCoupled<TElementData>::load(Serializer& rSerializer)
{
    typedef QSVMS<TElementData> BaseElement;
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseElement);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Class template instantiation

template class AlternativeQSVMSDEMCoupled<QSVMSDEMCoupledData<2,3> >;
template class AlternativeQSVMSDEMCoupled<QSVMSDEMCoupledData<3,4> >;

template class AlternativeQSVMSDEMCoupled< QSVMSDEMCoupledData<2,4> >;
template class AlternativeQSVMSDEMCoupled< QSVMSDEMCoupledData<3,8> >;
}

