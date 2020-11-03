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
    ProcessInfo& rCurrentProcessInfo)
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
    ProcessInfo& rCurrentProcessInfo)
{
        TElementData data;
        data.Initialize(*this, rCurrentProcessInfo);

        // Calculate this element RHS contribution
        QSVMS<TElementData>::CalculateRightHandSide(rRightHandSideVector, rCurrentProcessInfo);
        this->AddMassRHS(rRightHandSideVector, data);
}

template<class TElementData>
void QSVMSDEMCoupled<TElementData>::AddMassStabilization(
    TElementData& rData,
    MatrixType& rMassMatrix)
{

        const double density = this->GetAtCoordinate(rData.Density, rData.N);

        double tau_one;
        double tau_two;
        const array_1d<double, 3> convective_velocity=
            this->GetAtCoordinate(rData.Velocity, rData.N) -
            this->GetAtCoordinate(rData.MeshVelocity, rData.N);

        this->CalculateTau(rData, convective_velocity, tau_one, tau_two);

        double K; // Temporary results
        const double weight = rData.Weight * tau_one * density; // This density is for the dynamic term in the residual (rho*Du/Dt)
        // If we want to use more than one Gauss point to integrate the convective term, this has to be evaluated once per integration point

        Vector AGradN;
        this->ConvectionOperator(AGradN, convective_velocity, rData.DN_DX); // Get a * grad(Ni)

        AGradN *= density;

        const double fluid_fraction = this->GetAtCoordinate(rData.FluidFraction, rData.N);
        array_1d<double, 3> fluid_fraction_gradient = this->GetAtCoordinate(rData.FluidFractionGradient, rData.N);
        // Note: Dof order is (vx,vy,[vz,]p) for each node
        for (unsigned int i = 0; i < NumNodes; ++i)
        {
            unsigned int row = i*BlockSize;
            // Loop over columns
            for (unsigned int j = 0; j < NumNodes; ++j)
            {
                unsigned int col = j*BlockSize;
                K = weight * AGradN[i] * rData.N[j];

                for (unsigned int d = 0; d < Dim; ++d) // iterate over dimensions for velocity Dofs in this node combination
                {
                    rMassMatrix(row+d, col+d) += K;
                    rMassMatrix(row+Dim,col+d) += weight * fluid_fraction * rData.DN_DX(i,d) * rData.N[j];
                    rMassMatrix(row+Dim,col+d) += weight * fluid_fraction_gradient[d] * rData.N[i] * rData.N[j]; // Delta(u) * TauOne * alpha * Grad(q)
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
        fluid_fraction_rate = this->GetAtCoordinate(rData.FluidFractionRate, rData.N);
        const auto& r_geom = this->GetGeometry();

        for (unsigned int i = 0; i < NumNodes; ++i)
        {
            mass_source += rData.N[i] * r_geom[i].FastGetSolutionStepValue(MASS_SOURCE);
        }
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

    double tau_one;
    double tau_two;
    const array_1d<double, 3> convective_velocity =
        this->GetAtCoordinate(rData.Velocity, rData.N) -
        this->GetAtCoordinate(rData.MeshVelocity, rData.N);

    this->CalculateTau(rData, convective_velocity, tau_one, tau_two);

    Vector AGradN;
    this->ConvectionOperator(AGradN, convective_velocity, rData.DN_DX);

    // Multiplying some quantities by density to have correct units
    AGradN *= density; // Convective term is always multiplied by density

    const double fluid_fraction = this->GetAtCoordinate(rData.FluidFraction, rData.N);

    array_1d<double, 3> fluid_fraction_gradient = this->GetAtCoordinate(rData.FluidFractionGradient, rData.N);

    // Temporary containers
    double V, AA, P, GAlpha, AG, U, Q, DD, UD;

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
            AA = rData.Weight * AGradN[j] * tau_one * (AGradN[i]); // Stabilization: u*grad(v) * tau_one * u*grad(u)

            // q-p stabilization block (initialize result)
            double G = 0;
            for (unsigned int d = 0; d < Dim; d++)
            {
                LHS(row+d,col+d) += V + AA;

                // Stabilization: (a * Grad(v)) * tau_one * Grad(p)
                P = rData.DN_DX(i,d) * rData.N[j]; // Div(v) * p
                GAlpha = tau_one * AGradN[j] * (fluid_fraction * rData.DN_DX(i,d));
                AG = tau_one * AGradN[i] * rData.DN_DX(j,d);
                U = fluid_fraction_gradient[d] * rData.N[j] * rData.N[i];
                Q = fluid_fraction * rData.DN_DX(j,d) * rData.N[i];

                LHS(row+d,col+Dim) += rData.Weight * (AG - P);
                LHS(row+Dim,col+d) += rData.Weight * (GAlpha + U + Q);

                G += tau_one * fluid_fraction * rData.DN_DX(j,d) * rData.DN_DX(i,d);

                for (unsigned int e = 0; e < Dim; e++){ // Stabilization: Div(v) * tau_two * Div(u)
                    DD = tau_two * (rData.DN_DX(i,d) * fluid_fraction * rData.DN_DX(j,e));
                    UD = tau_two * rData.DN_DX(i,d) * fluid_fraction_gradient[e] * rData.N[j];
                    LHS(row+d,col+e) += rData.Weight * (DD + UD);
                }
            }
        // Write q-p term
        LHS(row+Dim,col+Dim) += rData.Weight * G;

        }

        // RHS terms
        double QF = 0.0;
        for (unsigned int d = 0; d < Dim; ++d)
        {
            rLocalRHS[row+d] += rData.Weight * rData.N[i] * body_force[d]; // v*BodyForce
            rLocalRHS[row+d] += rData.Weight * tau_one * AGradN[i] * (body_force[d] - momentum_projection[d]); // ( a * Grad(v) ) * tau_one * (Density * BodyForce)
            rLocalRHS[row+d] -= rData.Weight * tau_two * rData.DN_DX(i,d) * (mass_projection);
            QF += tau_one * (body_force[d] - momentum_projection[d]) * (fluid_fraction * rData.DN_DX(i,d));

        }
        rLocalRHS[row+Dim] += rData.Weight * (QF);
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

        const auto fluid_fraction_gradient = rData.FluidFractionGradient;

        const double fluid_fraction_rate = this->GetAtCoordinate(rData.FluidFractionRate, rData.N);
        // Compute this node's contribution to the residual (evaluated at integration point)
        for (unsigned int i = 0; i < NumNodes; i++) {
            for (unsigned int d = 0; d < Dim; ++d)
            {
                rMassRHS -= (fluid_fraction * rData.DN_DX(i, d) * velocities(i, d)) + fluid_fraction_gradient(i,d) * rData.N[i] * velocities(i, d);
            }
        }
        rMassRHS -= fluid_fraction_rate;
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

