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
#include "custom_utilities/qsvms_dem_coupled_data.h"
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

///////////////////////////////////////////////////////////////////////////////////////////////////
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

template <class TElementData>
Matrix DVMSDEMCoupled<TElementData>::GetAtCoordinate(
    const typename TElementData::NodalTensorData& rValues,
    const typename TElementData::ShapeFunctionsType& rN) const
{
    Matrix result = ZeroMatrix(3);

    for (size_t i = 0; i < NumNodes; i++) {
        for (size_t j = 0; j < Dim; j++) {
            for (size_t k = 0; k < Dim; k++) {
                result(j, k) += rN[i] * rValues[i](j, k);
            }
        }
    }

    return result;
}

template <class TElementData>
double DVMSDEMCoupled<TElementData>::GetAtCoordinate(
    const typename TElementData::NodalScalarData& rValues,
    const typename TElementData::ShapeFunctionsType& rN) const
{
    return DVMS<TElementData>::GetAtCoordinate(rValues, rN);
}

template <class TElementData>
array_1d<double, 3> DVMSDEMCoupled<TElementData>::GetAtCoordinate(
    const typename TElementData::NodalVectorData& rValues,
    const typename TElementData::ShapeFunctionsType& rN) const
{
    return DVMS<TElementData>::GetAtCoordinate(rValues, rN);
}

template <class TElementData>
double DVMSDEMCoupled<TElementData>::GetAtCoordinate(
    const double Value,
    const typename TElementData::ShapeFunctionsType& rN) const
{
    return DVMS<TElementData>::GetAtCoordinate(Value, rN);
}

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
    Matrix permeability = this->GetAtCoordinate(rData.Permeability, rData.N);
    const auto& r_body_forces = rData.BodyForce;
    const auto& r_velocities = rData.Velocity;
    const auto& r_pressures = rData.Pressure;

    for (unsigned int i = 0; i < NumNodes; i++) {
        const array_1d<double,3>& r_acceleration = rGeom[i].FastGetSolutionStepValue(ACCELERATION);
        for (unsigned int d = 0; d < Dim; d++) {
            Vector sigma = ZeroVector(Dim);
            for (unsigned int e = 0; e < Dim; e++){
                sigma[d] += (viscosity/permeability(d,e)) * rData.N[i] * r_velocities(i,e);
            }
            rResidual[d] += density * ( rData.N[i]*(r_body_forces(i,d) - r_acceleration[d]) - convection[i]*r_velocities(i,d)) - rData.DN_DX(i,d)*r_pressures[i] - sigma[d];
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
    Matrix permeability = this->GetAtCoordinate(rData.Permeability, rData.N);

    for (unsigned int i = 0; i < NumNodes; i++) {
        for (unsigned int d = 0; d < Dim; d++) {
            Vector sigma = ZeroVector(Dim);
            for (unsigned int e = 0; e < Dim; e++){
                sigma[d] += (viscosity/permeability(d,e)) * rData.N[i] * rData.Velocity(i,e);
            }
            rMomentumRHS[d] += density * ( rData.N[i]*(rData.BodyForce(i,d) /*- rAcc[d]*/) - AGradN[i]*rData.Velocity(i,d)) - rData.DN_DX(i,d)*rData.Pressure[i] - sigma[d];
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

    double tau_one;
    double tau_two;
    this->CalculateStabilizationParameters(rData,convective_velocity,tau_one,tau_two);

    const double dt = rData.DeltaTime;
    DenseVector< array_1d<double,Dim> > mOldSubscaleVelocity = DVMS<TElementData>::mOldSubscaleVelocity;
    // small scale velocity contributions (subscale tracking)
    array_1d<double,Dim> OldUssTerm = (density/dt) * mOldSubscaleVelocity[rData.IntegrationPointIndex]; // rho * u_ss^{n-1}/dt

    Vector AGradN;
    this->ConvectionOperator(AGradN,convective_velocity,rData.DN_DX);

    // These two should be zero unless we are using OSS
    double viscosity = this->GetAtCoordinate(rData.DynamicViscosity, rData.N);
    const array_1d<double,3> MomentumProj = this->GetAtCoordinate(rData.MomentumProjection,rData.N);
    const double MassProj = this->GetAtCoordinate(rData.MassProjection,rData.N);
    const double fluid_fraction = this->GetAtCoordinate(rData.FluidFraction, rData.N);
    const double fluid_fraction_rate = this->GetAtCoordinate(rData.FluidFractionRate, rData.N);
    const double mass_source = this->GetAtCoordinate(rData.MassSource, rData.N);
    Matrix permeability = this->GetAtCoordinate(rData.Permeability, rData.N);
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

                /* q-p stabilization block */
                // Stabilization: Grad(q) * TauOne * Grad(p)
                G += fluid_fraction * rData.DN_DX(i,d) * rData.DN_DX(j,d);

                /* v-u block */
                // Stabilization: Div(v) * TauTwo * Div(u)
                double GAlphaR = 0.0;
                double RSigmaG = 0.0;
                for (unsigned int e = 0; e < Dim; e++){
                    double RSigma = rData.N[i] * (viscosity / permeability(d,e)) * rData.N[j];
                    double VSigma = density/dt * rData.N[i] * rData.N[j] * (viscosity / permeability(d,e));
                    double ASigma = tau_one * AGradN[i] * (viscosity / permeability(d,e)) * rData.N[j];
                    double RRSigma = tau_one * (viscosity / permeability(d,e)) * rData.N[i] * (viscosity / permeability(d,e)) * rData.N[j];
                    double RSigmaA = tau_one * (viscosity / permeability(d,e)) * rData.N[i] * AGradN[j];
                    double DAlphaD = fluid_fraction * tau_two *rData.DN_DX(i,d)*rData.DN_DX(j,e);
                    double DU = fluid_fraction_gradient[e] * tau_two * rData.DN_DX(i,d)*rData.N[j];
                    GAlphaR += tau_one * fluid_fraction * rData.DN_DX(i,d) * (viscosity / permeability(d,e)) * rData.N[j];
                    RSigmaG += tau_one * (viscosity / permeability(d,e)) * rData.N[i] * rData.DN_DX(j,d);

                    LHS(row+d,col+e) += rData.Weight * (DAlphaD + DU + RSigma + VSigma + ASigma + RRSigma + RSigmaA);
                }

                LHS(col+Dim,row+d) += rData.Weight * (GAlphaA + U + QAlpha + GAlphaR);
                LHS(row+d,col+Dim) += rData.Weight * (AG - VP - P + RSigmaG);

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
            double RSigmaF =0.0;
            for (unsigned int e = 0; e < Dim; ++e){
                RSigmaF += tau_one * (viscosity/permeability(d,e)) * rData.N[i] * (body_force[e] - MomentumProj[e] + OldUssTerm[e]);
            }
            // ( a * Grad(v) ) * TauOne * (Density * BodyForce - Projection)
            // vh * TauOne/Dt * f (from vh*d(uss)/dt
            double VI = tau_one * density * rData.N[i] / dt * (body_force[d] - MomentumProj[d] + OldUssTerm[d]);
            double AF = tau_one * AGradN[i] * (body_force[d] - MomentumProj[d] + OldUssTerm[d]);

            // OSS pressure subscale projection
            double DPhi = rData.DN_DX(i,d) * tau_two * (mass_source - fluid_fraction_rate - MassProj);

            rLocalRHS[row+d] += rData.Weight * (VF - VI + AF - DPhi + RSigmaF);

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
    this->CalculateStabilizationParameters(rData,convective_velocity,tau_one,tau_two);

    const double dt = rData.DeltaTime;

    Vector AGradN;
    this->ConvectionOperator(AGradN,convective_velocity,rData.DN_DX);

    // Multiplying convective operator by density to have correct units
    AGradN *= density;

    const double fluid_fraction = this->GetAtCoordinate(rData.FluidFraction, rData.N);
    double viscosity = this->GetAtCoordinate(rData.DynamicViscosity, rData.N);
    Matrix permeability = this->GetAtCoordinate(rData.Permeability, rData.N);
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
                for (unsigned int e = 0; e < Dim; ++e){
                    double RSigmaU = (viscosity / permeability(d,e)) * rData.N[i] * AGradN[j];
                    rMassMatrix(row+d, col+e) += W * RSigmaU;
                }

            }
        }
    }
}

template< class TElementData >
void DVMSDEMCoupled<TElementData>::CalculateStabilizationParameters(
    const TElementData& rData,
    const array_1d<double,3> &Velocity,
    double &TauOne,
    double &TauTwo) const
{
    const double h = rData.ElementSize;
    const double density = this->GetAtCoordinate(rData.Density,rData.N);
    const double viscosity = this->GetAtCoordinate(rData.EffectiveViscosity,rData.N);
    constexpr double mTauC1 = DVMS<TElementData>::mTauC1;
    constexpr double mTauC2 = DVMS<TElementData>::mTauC2;

    double velocity_norm = Velocity[0]*Velocity[0];
    for (unsigned int d = 1; d < Dim; d++)
        velocity_norm += Velocity[d]*Velocity[d];
    velocity_norm = std::sqrt(velocity_norm);

    double inv_tau = mTauC1 * viscosity / (h*h) + density * ( 1.0/rData.DeltaTime + mTauC2 * velocity_norm / h );
    TauOne = 1.0/inv_tau;
    TauTwo = viscosity + density * mTauC2 * velocity_norm * h / mTauC1;

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

template class DVMSDEMCoupled< QSVMSDEMCoupledData<2,4> >;
template class DVMSDEMCoupled< QSVMSDEMCoupledData<3,8> >;

} // namespace Kratos
