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
#include "qs_vms_dem_coupled.h"
#include "custom_utilities/qsvms_dem_coupled_data.h"
#include "custom_utilities/fluid_element_utilities.h"
#include "fluid_dynamics_application_variables.h"

namespace Kratos
{

//////////////////////////Life cycle

template< class TElementData >
QSVMSDEMCoupled<TElementData>::QSVMSDEMCoupled(IndexType NewId):
    QSVMS<TElementData>(NewId)
{}

template< class TElementData >
QSVMSDEMCoupled<TElementData>::QSVMSDEMCoupled(IndexType NewId, const NodesArrayType& ThisNodes):
    QSVMS<TElementData>(NewId,ThisNodes)
{}


template< class TElementData >
QSVMSDEMCoupled<TElementData>::QSVMSDEMCoupled(IndexType NewId, GeometryType::Pointer pGeometry):
    QSVMS<TElementData>(NewId,pGeometry)
{}


template< class TElementData >
QSVMSDEMCoupled<TElementData>::QSVMSDEMCoupled(IndexType NewId, GeometryType::Pointer pGeometry, Properties::Pointer pProperties):
    QSVMS<TElementData>(NewId,pGeometry,pProperties)
{}

///////////Destructor

template< class TElementData >
QSVMSDEMCoupled<TElementData>::~QSVMSDEMCoupled()
{}

template< class TElementData >
Element::Pointer QSVMSDEMCoupled<TElementData>::Create(IndexType NewId,NodesArrayType const& ThisNodes,Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<QSVMSDEMCoupled>(NewId, this->GetGeometry().Create(ThisNodes), pProperties);
}

template< class TElementData >
Element::Pointer QSVMSDEMCoupled<TElementData>::Create(IndexType NewId,GeometryType::Pointer pGeom,Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<QSVMSDEMCoupled>(NewId, pGeom, pProperties);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Inquiry

template< class TElementData >
int QSVMSDEMCoupled<TElementData>::Check(const ProcessInfo &rCurrentProcessInfo) const
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

template <class TElementData>
void QSVMSDEMCoupled<TElementData>::Calculate(
    const Variable<double>& rVariable,
    double& rOutput,
    const ProcessInfo& rCurrentProcessInfo)
{
    QSVMS<TElementData>::Calculate(rVariable, rOutput, rCurrentProcessInfo);
}

template <class TElementData>
void QSVMSDEMCoupled<TElementData>::Calculate(
    const Variable<array_1d<double, 3>>& rVariable,
    array_1d<double, 3>& rOutput,
    const ProcessInfo& rCurrentProcessInfo)
{
    QSVMS<TElementData>::Calculate(rVariable, rOutput, rCurrentProcessInfo);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @see QSVMSDEMCoupled::EquationIdVector
 **/
template < class TElementData >
void QSVMSDEMCoupled<TElementData>::EquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo) const
{
    QSVMS<TElementData>::EquationIdVector(rResult, rCurrentProcessInfo);
}

template< class TElementData >
std::string QSVMSDEMCoupled<TElementData>::Info() const
{
    std::stringstream buffer;
    buffer << "QSVMSDEMCoupled #" << this->Id();
    return buffer.str();
}


template< class TElementData >
void QSVMSDEMCoupled<TElementData>::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "QSVMSDEMCoupled" << Dim << "D";
}

template<class TElementData>
void QSVMSDEMCoupled<TElementData>::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
        TElementData data;
        data.Initialize(*this, rCurrentProcessInfo);

        // Calculate this element RHS contribution
        QSVMS<TElementData>::CalculateRightHandSide(rRightHandSideVector, rCurrentProcessInfo);
}

template< class TElementData >
void QSVMSDEMCoupled<TElementData>::AlgebraicMomentumResidual(
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
        const array_1d<double,Dim>& r_acceleration = rGeom[i].FastGetSolutionStepValue(ACCELERATION);
        array_1d<double,Dim> sigma_U = ZeroVector(Dim);
        for (unsigned int d = 0; d < Dim; d++) {
            for (unsigned int e = 0; e < Dim; e++){
                sigma_U[d] += sigma(d,e) * rData.N[i] * r_velocities(i,e);
            }
            rResidual[d] += density * ( rData.N[i]*(r_body_forces(i,d) - r_acceleration[d]) - convection[i]*r_velocities(i,d))
                            - rData.DN_DX(i,d)*r_pressures[i] - sigma_U[d];
        }
    }
}

template< class TElementData >
void QSVMSDEMCoupled<TElementData>::MomentumProjTerm(
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
        Vector sigma_U = ZeroVector(Dim);
        for (unsigned int d = 0; d < Dim; d++) {
            for (unsigned int e = 0; e < Dim; e++){
                sigma_U[d] += sigma(d,e) * rData.N[i] * rData.Velocity(i,e);
            }
            rMomentumRHS[d] += density * ( rData.N[i]*(rData.BodyForce(i,d) /*- rAcc[d]*/) - AGradN[i]*rData.Velocity(i,d))
                                - rData.DN_DX(i,d)*rData.Pressure[i] - sigma_U[d];
        }
    }
}

template<class TElementData>
void QSVMSDEMCoupled<TElementData>::AddMassStabilization(
    TElementData& rData,
    MatrixType& rMassMatrix)
{

    const double density = this->GetAtCoordinate(rData.Density, rData.N);

    BoundedMatrix<double,Dim,Dim> tau_one = ZeroMatrix(Dim, Dim);
    double tau_two;
    const array_1d<double, 3> convective_velocity=
        this->GetAtCoordinate(rData.Velocity, rData.N) -
        this->GetAtCoordinate(rData.MeshVelocity, rData.N);

    this->CalculateTau(rData, convective_velocity, tau_one, tau_two);

    const double weight = rData.Weight * density; // This density is for the dynamic term in the residual (rho*Du/Dt)
        // If we want to use more than one Gauss point to integrate the convective term, this has to be evaluated once per integration point

    Vector AGradN;
    this->ConvectionOperator(AGradN, convective_velocity, rData.DN_DX); // Get a * grad(Ni)

    AGradN *= density;

    const double fluid_fraction = this->GetAtCoordinate(rData.FluidFraction, rData.N);
    double viscosity = this->GetAtCoordinate(rData.DynamicViscosity, rData.N);
    BoundedMatrix<double,Dim,Dim> permeability = this->GetAtCoordinate(rData.Permeability, rData.N);
    BoundedMatrix<double,Dim,Dim> sigma = ZeroMatrix(Dim, Dim);

    double det_permeability = MathUtils<double>::Det(permeability);
    MathUtils<double>::InvertMatrix(permeability, sigma, det_permeability, -1.0);

    sigma *= viscosity;

    // Note: Dof order is (vx,vy,[vz,]p) for each node
    for (unsigned int i = 0; i < NumNodes; ++i)
    {
        unsigned int row = i*BlockSize;
        // Loop over columns
        for (unsigned int j = 0; j < NumNodes; ++j)
        {
            unsigned int col = j*BlockSize;

            for (unsigned int d = 0; d < Dim; ++d) // iterate over dimensions for velocity Dofs in this node combination
            {
                double K = weight * tau_one(d,d) * AGradN[i] * rData.N[j];
                for (unsigned int e = 0; e < Dim; ++e){
                    double RSigmaU = tau_one(d,d) * sigma(d,e) * rData.N[i] * AGradN[j];
                    rMassMatrix(row+d, col+e) += weight * RSigmaU;
                }
                rMassMatrix(row+d, col+d) += K;
                rMassMatrix(row+Dim,col+d) += weight * tau_one(d,d) * (fluid_fraction * rData.DN_DX(i,d) * rData.N[j]);
            }
        }
    }
}

template<class TElementData>
void QSVMSDEMCoupled<TElementData>::AddMassRHS(
    VectorType& rRightHandSideVector,
    TElementData& rData)
{
        double fluid_fraction_rate = 0.0;
        double mass_source = 0.0;
        mass_source = this->GetAtCoordinate(rData.MassSource, rData.N);
        fluid_fraction_rate = this->GetAtCoordinate(rData.FluidFractionRate, rData.N);

        // Add the results to the pressure components (Local Dofs are vx, vy, [vz,] p for each node)
        int LocalIndex = Dim;
        for (unsigned int i = 0; i < NumNodes; ++i){
            for (unsigned int d = 0; d < Dim;++d)
                rRightHandSideVector[LocalIndex] -= rData.Weight * rData.N[i] * (fluid_fraction_rate - mass_source);
            LocalIndex += Dim + 1;
        }

}

// Add a the contribution from a single integration point to the velocity contribution
template< class TElementData >
void QSVMSDEMCoupled<TElementData>::AddVelocitySystem(
    TElementData& rData,
    MatrixType& rLocalLHS,
    VectorType& rLocalRHS)
{
    auto& LHS = rData.LHS;
    LHS.clear();

    // Interpolate nodal data on the integration point
    const double density = this->GetAtCoordinate(rData.Density, rData.N);
    const array_1d<double, 3> body_force = density * this->GetAtCoordinate(rData.BodyForce,rData.N);
    const array_1d<double,3> momentum_projection = this->GetAtCoordinate(rData.MomentumProjection, rData.N);
    double mass_projection = this->GetAtCoordinate(rData.MassProjection, rData.N);

    BoundedMatrix<double,Dim,Dim> tau_one = ZeroMatrix(Dim, Dim);
    double tau_two;
    const array_1d<double, 3> convective_velocity =
        this->GetAtCoordinate(rData.Velocity, rData.N) -
        this->GetAtCoordinate(rData.MeshVelocity, rData.N);

    this->CalculateTau(rData, convective_velocity, tau_one, tau_two);

    Vector AGradN;
    this->ConvectionOperator(AGradN, convective_velocity, rData.DN_DX);

    // Multiplying some quantities by density to have correct units
    AGradN *= density; // Convective term is always multiplied by density

    double viscosity = this->GetAtCoordinate(rData.DynamicViscosity, rData.N);
    const double fluid_fraction = this->GetAtCoordinate(rData.FluidFraction, rData.N);
    const double fluid_fraction_rate = this->GetAtCoordinate(rData.FluidFractionRate, rData.N);
    const double mass_source = this->GetAtCoordinate(rData.MassSource, rData.N);
    BoundedMatrix<double,Dim,Dim> permeability = this->GetAtCoordinate(rData.Permeability, rData.N);
    array_1d<double, 3> fluid_fraction_gradient = this->GetAtCoordinate(rData.FluidFractionGradient, rData.N);
    BoundedMatrix<double,Dim,Dim> sigma = ZeroMatrix(Dim, Dim);

    double det_permeability = MathUtils<double>::Det(permeability);
    MathUtils<double>::InvertMatrix(permeability, sigma, det_permeability, -1.0);

    sigma *= viscosity;

    // Temporary containers
    double V, P, U, QAlpha, DAlphaD, DU, RSigma, ASigma, RRSigma, RSigmaA;

    // Note: Dof order is (u,v,[w,]p) for each node
    for (unsigned int i = 0; i < NumNodes; i++)
    {

        unsigned int row = i*BlockSize;

        // LHS terms
        for (unsigned int j = 0; j < NumNodes; j++)
        {
            unsigned int col = j*BlockSize;

            // Some terms are the same for all velocity components, calculate them once for each i,j
            V = rData.Weight * rData.N[i] * AGradN[j];

            // q-p stabilization block (initialize result)
            double G = 0;
            for (unsigned int d = 0; d < Dim; d++)
            {

                // Stabilization: (a * Grad(v)) * tau_one * Grad(p)
                P = rData.DN_DX(i,d) * rData.N[j]; // Div(v) * p
                U = fluid_fraction_gradient[d] * rData.N[j] * rData.N[i];
                QAlpha = fluid_fraction * rData.DN_DX(j,d) * rData.N[i];

                double GAlphaR = 0.0;
                double RSigmaG = 0.0;
                double GAlphaA = tau_one(d,d) * AGradN[j] * (fluid_fraction * rData.DN_DX(i,d));
                double AG = tau_one(d,d) * AGradN[i] * rData.DN_DX(j,d);
                G += tau_one(d,d) * fluid_fraction * rData.DN_DX(i,d) * rData.DN_DX(j,d);
                double AA = rData.Weight * tau_one(d,d) * AGradN[j] * AGradN[i]; // Stabilization: u*grad(v) * tau_one * u*grad(u);

                for (unsigned int e = 0; e < Dim; e++){ // Stabilization: Div(v) * tau_two * Div(u)
                    RSigma = rData.N[i] * sigma(d,e) * rData.N[j];
                    ASigma = tau_one(d,d) * AGradN[i] * sigma(d,e) * rData.N[j];
                    RRSigma = tau_one(d,d) * sigma(d,e) * rData.N[i] * sigma(e,d) * rData.N[j];
                    RSigmaA = tau_one(d,d) * sigma(d,e) * rData.N[i] * AGradN[j];
                    DAlphaD = tau_two * fluid_fraction * rData.DN_DX(i,d) * rData.DN_DX(j,e);
                    DU = tau_two * rData.DN_DX(i,d) * fluid_fraction_gradient[e] * rData.N[j];
                    GAlphaR += tau_one(d,d) * fluid_fraction * rData.DN_DX(i,d) * sigma(d,e) * rData.N[j];
                    RSigmaG += tau_one(d,d) * sigma(d,e) * rData.N[i] * rData.DN_DX(j,d);
                    LHS(row+d,col+e) += rData.Weight * (DAlphaD + DU + RSigma + ASigma + RRSigma + RSigmaA);
                }

                LHS(row+d,col+d) += V + AA;
                LHS(row+Dim,col+d) += rData.Weight * (GAlphaA + U + QAlpha + GAlphaR);
                LHS(row+d,col+Dim) += rData.Weight * (AG - P + RSigmaG);

            }
            // Write q-p term
            LHS(row+Dim,col+Dim) += rData.Weight * G;

        }

        // RHS terms
        double QAlphaFplusGAlphaF = 0.0;
        for (unsigned int d = 0; d < Dim; ++d)
        {
            rLocalRHS[row+d] += rData.Weight * rData.N[i] * body_force[d]; // v*BodyForce
            rLocalRHS[row+d] += rData.Weight * tau_one(d,d) * AGradN[i] * (body_force[d] - momentum_projection[d]); // A_F: ( a * Grad(v) ) * tau_one * (Density * BodyForce)
            rLocalRHS[row+d] -= rData.Weight * tau_two * rData.DN_DX(i,d) * (mass_projection + mass_source - fluid_fraction_rate);
            QAlphaFplusGAlphaF += tau_one(d,d) * (body_force[d] - momentum_projection[d]) * fluid_fraction * rData.DN_DX(i,d);
        }
        rLocalRHS[row+Dim] += rData.Weight * (QAlphaFplusGAlphaF + rData.N[i] * (mass_source - fluid_fraction_rate));
    }

    // Write (the linearized part of the) local contribution into residual form (A*dx = b - A*x)
    array_1d<double,LocalSize> values;
    this->GetCurrentValuesVector(rData, values);
    noalias(rLocalRHS) -= prod(LHS, values);
    /* Viscous contribution (with symmetric gradient 2*nu*{E(u) - 1/3 Tr(E)} )
     * For a generic (potentially non-linear) constitutive law, one cannot assume that RHS = F - LHS*current_values.
     * Because of this, the AddViscousTerm function manages both the LHS and the RHS.
     */
    this->AddViscousTerm(rData, LHS, rLocalRHS);

    noalias(rLocalLHS) += LHS;

}

template<class TElementData>
void QSVMSDEMCoupled<TElementData>::MassProjTerm(
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

template< class TElementData >
void QSVMSDEMCoupled<TElementData>::CalculateTau(
    const TElementData& rData,
    const array_1d<double,3> &Velocity,
    BoundedMatrix<double,Dim,Dim> &TauOne,
    double &TauTwo) const
{
    constexpr double c1 = 8.0;
    constexpr double c2 = 2.0;

    const double h = rData.ElementSize;
    const double density = this->GetAtCoordinate(rData.Density,rData.N);
    const double viscosity = this->GetAtCoordinate(rData.EffectiveViscosity,rData.N);
    BoundedMatrix<double,Dim,Dim> permeability = this->GetAtCoordinate(rData.Permeability, rData.N);
    BoundedMatrix<double,Dim,Dim> sigma = ZeroMatrix(Dim, Dim);
    BoundedMatrix<double,Dim,Dim> non_diag_tau_one = ZeroMatrix(Dim, Dim);
    BoundedMatrix<double,Dim,Dim> inv_tau = ZeroMatrix(Dim, Dim);
    BoundedMatrix<double,Dim,Dim> I = IdentityMatrix(Dim, Dim);
    BoundedMatrix<double,Dim,Dim> eigen_values_matrix, eigen_vectors_matrix;
    BoundedMatrix<double,Dim,Dim> inv_eigen_matrix = ZeroMatrix(Dim, Dim);

    double det_permeability = MathUtils<double>::Det(permeability);
    MathUtils<double>::InvertMatrix(permeability, sigma, det_permeability, -1.0);

    double velocity_norm = Velocity[0]*Velocity[0];
    for (unsigned int d = 1; d < Dim; d++)
        velocity_norm += Velocity[d]*Velocity[d];
    velocity_norm = std::sqrt(velocity_norm);

    inv_tau = (c1 * viscosity / (h*h) + density * ( rData.DynamicTau/rData.DeltaTime + c2 * velocity_norm / h )) * I + viscosity * sigma;

    double det_inv_tau = MathUtils<double>::Det(inv_tau);
    MathUtils<double>::InvertMatrix(inv_tau, non_diag_tau_one, det_inv_tau, -1.0);

    MathUtils<double>::GaussSeidelEigenSystem<BoundedMatrix<double,Dim,Dim>, BoundedMatrix<double,Dim,Dim>>(non_diag_tau_one, eigen_vectors_matrix, eigen_values_matrix);

    double det_eigen_vectors_matrix = MathUtils<double>::Det(eigen_vectors_matrix);
    MathUtils<double>::InvertMatrix(eigen_vectors_matrix, inv_eigen_matrix, det_eigen_vectors_matrix, -1.0);

    BoundedMatrix<double,Dim,Dim> inv_PTau = prod(inv_eigen_matrix, non_diag_tau_one);
    TauOne = prod(inv_PTau, eigen_vectors_matrix);
    TauTwo = viscosity + c2 * density * velocity_norm * h / c1;
}

template< class TElementData >
void QSVMSDEMCoupled<TElementData>::SubscaleVelocity(
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
void QSVMSDEMCoupled<TElementData>::SubscalePressure(
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
void QSVMSDEMCoupled<TElementData>::save(Serializer& rSerializer) const
{
    typedef QSVMS<TElementData> BaseElement;
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseElement );
}


template< class TElementData >
void QSVMSDEMCoupled<TElementData>::load(Serializer& rSerializer)
{
    typedef QSVMS<TElementData> BaseElement;
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseElement);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Class template instantiation

template class QSVMSDEMCoupled<QSVMSDEMCoupledData<2,3> >;
template class QSVMSDEMCoupled<QSVMSDEMCoupledData<3,4> >;

template class QSVMSDEMCoupled< QSVMSDEMCoupledData<2,4> >;
template class QSVMSDEMCoupled< QSVMSDEMCoupledData<3,8> >;
}

