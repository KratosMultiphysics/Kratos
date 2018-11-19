//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela
//                   Ignasi de Pouplana
//

#include "fic.h"
#include "includes/cfd_variables.h"
#include "includes/checks.h"

#include "custom_utilities/fic_data.h"
#include "custom_utilities/time_integrated_fic_data.h"
#include "custom_utilities/fluid_element_utilities.h"
#include "custom_utilities/element_size_calculator.h"
#include "custom_utilities/fluid_element_time_integration_detail.h"

namespace Kratos
{

///////////////////////////////////////////////////////////////////////////////////////////////////
// Life cycle

template< class TElementData >
FIC<TElementData>::FIC(IndexType NewId):
    FluidElement<TElementData>(NewId)
{}

template< class TElementData >
FIC<TElementData>::FIC(IndexType NewId, const NodesArrayType& ThisNodes):
    FluidElement<TElementData>(NewId,ThisNodes)
{}


template< class TElementData >
FIC<TElementData>::FIC(IndexType NewId, GeometryType::Pointer pGeometry):
    FluidElement<TElementData>(NewId,pGeometry)
{}


template< class TElementData >
FIC<TElementData>::FIC(IndexType NewId, GeometryType::Pointer pGeometry, Properties::Pointer pProperties):
    FluidElement<TElementData>(NewId,pGeometry,pProperties)
{}


template< class TElementData >
FIC<TElementData>::~FIC()
{}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Public Operations

template< class TElementData >
Element::Pointer FIC<TElementData>::Create(IndexType NewId,NodesArrayType const& ThisNodes,Properties::Pointer pProperties) const
{
    return Kratos::make_shared<FIC>(NewId, this->GetGeometry().Create(ThisNodes), pProperties);
}


template< class TElementData >
Element::Pointer FIC<TElementData>::Create(IndexType NewId,GeometryType::Pointer pGeom,Properties::Pointer pProperties) const
{
    return Kratos::make_shared<FIC>(NewId, pGeom, pProperties);
}

template <class TElementData>
void FIC<TElementData>::Calculate(const Variable<double>& rVariable,
    double& rOutput, const ProcessInfo& rCurrentProcessInfo) {}

template <class TElementData>
void FIC<TElementData>::Calculate(
    const Variable<array_1d<double, 3>>& rVariable,
    array_1d<double, 3>& rOutput, const ProcessInfo& rCurrentProcessInfo) {}

template <class TElementData>
void FIC<TElementData>::Calculate(const Variable<Vector>& rVariable,
    Vector& rOutput, const ProcessInfo& rCurrentProcessInfo) {}

template <class TElementData>
void FIC<TElementData>::Calculate(const Variable<Matrix>& rVariable,
    Matrix& rOutput, const ProcessInfo& rCurrentProcessInfo) {}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Inquiry

template< class TElementData >
int FIC<TElementData>::Check(const ProcessInfo &rCurrentProcessInfo)
{
    int out = FluidElement<TElementData>::Check(rCurrentProcessInfo);
    KRATOS_ERROR_IF_NOT(out == 0)
        << "Error in base class Check for Element " << this->Info() << std::endl
        << "Error code is " << out << std::endl;

    // Extra variables
    KRATOS_CHECK_VARIABLE_KEY(ACCELERATION);

    for(unsigned int i=0; i<NumNodes; ++i)
    {
        Node<3>& rNode = this->GetGeometry()[i];
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(ACCELERATION,rNode);
    }

    return out;
}

///////////////////////////////////////////////////////////////////////////////////////////////////

template< class TElementData >
void FIC<TElementData>::GetValueOnIntegrationPoints(
    Variable<array_1d<double, 3 > > const& rVariable,
    std::vector<array_1d<double, 3 > >& rValues,
    ProcessInfo const& rCurrentProcessInfo)
{
    FluidElement<TElementData>::GetValueOnIntegrationPoints(rVariable,rValues,rCurrentProcessInfo);
}

template< class TElementData >
void FIC<TElementData>::GetValueOnIntegrationPoints(
    Variable<double> const& rVariable,
    std::vector<double>& rValues,
    ProcessInfo const& rCurrentProcessInfo)
{
    FluidElement<TElementData>::GetValueOnIntegrationPoints(rVariable,rValues,rCurrentProcessInfo);
}

template <class TElementData>
void FIC<TElementData>::GetValueOnIntegrationPoints(
    Variable<array_1d<double, 6>> const& rVariable,
    std::vector<array_1d<double, 6>>& rValues,
    ProcessInfo const& rCurrentProcessInfo)
{
    FluidElement<TElementData>::GetValueOnIntegrationPoints(rVariable,rValues,rCurrentProcessInfo);
}

template <class TElementData>
void FIC<TElementData>::GetValueOnIntegrationPoints(
    Variable<Vector> const& rVariable,
    std::vector<Vector>& rValues,
    ProcessInfo const& rCurrentProcessInfo)
{
    FluidElement<TElementData>::GetValueOnIntegrationPoints(rVariable,rValues,rCurrentProcessInfo);
}

template <class TElementData>
void FIC<TElementData>::GetValueOnIntegrationPoints(
    Variable<Matrix> const& rVariable,
    std::vector<Matrix>& rValues,
    ProcessInfo const& rCurrentProcessInfo)
{
    FluidElement<TElementData>::GetValueOnIntegrationPoints(rVariable,rValues,rCurrentProcessInfo);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Input and output

template< class TElementData >
std::string FIC<TElementData>::Info() const
{
    std::stringstream buffer;
    buffer << "FIC #" << this->Id();
    return buffer.str();
}


template< class TElementData >
void FIC<TElementData>::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "FIC" << Dim << "D";
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Protected functions

template< class TElementData >
void FIC<TElementData>::AlgebraicMomentumResidual(
    const TElementData& rData,
    const Vector& rConvection,
    array_1d<double,3>& rResidual) const
{
    const GeometryType rGeom = this->GetGeometry();

    const double density = this->GetAtCoordinate(rData.Density,rData.N);
    const auto& r_body_forces = rData.BodyForce;
    const auto& r_velocities = rData.Velocity;
    const auto& r_pressures = rData.Pressure;

    for (unsigned int i = 0; i < NumNodes; i++) {
        const array_1d<double,3>& r_acceleration = rGeom[i].FastGetSolutionStepValue(ACCELERATION);
        for (unsigned int d = 0; d < Dim; d++) {
            rResidual[d] += density * ( rData.N[i]*(r_body_forces(i,d) - r_acceleration[d]) - rConvection[i]*r_velocities(i,d)) - rData.DN_DX(i,d)*r_pressures[i];
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Evaluation of system terms on Gauss Points

template <class TElementData>
void FIC<TElementData>::AddTimeIntegratedSystem(
    TElementData& rData, MatrixType& rLHS, VectorType& rRHS) {

    // Call specialized implementation (it is on a helper class to avoid partial template specialization problems)
    Internals::FluidElementTimeIntegrationDetail<TElementData,
        TElementData::ElementManagesTimeIntegration>::AddTimeIntegratedSystem(this, rData,
        rLHS, rRHS);
}

template <class TElementData>
void FIC<TElementData>::AddTimeIntegratedLHS(
    TElementData& rData, MatrixType& rLHS) {
        KRATOS_ERROR << "AddTimeIntegratedLHS is not implemented." << std::endl;
    }

template <class TElementData>
void FIC<TElementData>::AddTimeIntegratedRHS(
    TElementData& rData, VectorType& rRHS) {
        KRATOS_ERROR << "AddTimeIntegratedRHS is not implemented." << std::endl;
    }

template< class TElementData >
void FIC<TElementData>::AddVelocitySystem(
    TElementData& rData,
    MatrixType &rLocalLHS,
    VectorType &rLocalRHS)
{
    auto& LHS = rData.LHS;
    LHS.clear();

    // Interpolate nodal data on the integration point
    const double density = this->GetAtCoordinate(rData.Density,rData.N);
    array_1d<double,3> body_force = this->GetAtCoordinate(rData.BodyForce,rData.N);

    double TauIncompr;
    double TauMomentum;
    array_1d<double,3> convective_velocity = this->GetAtCoordinate(rData.Velocity,rData.N) - this->GetAtCoordinate(rData.MeshVelocity,rData.N);

    array_1d<double,3> TauGrad = ZeroVector(3);
    this->CalculateTau(rData,convective_velocity,TauIncompr,TauMomentum,TauGrad);

    Vector AGradN;
    this->ConvectionOperator(AGradN,convective_velocity,rData.DN_DX);

    // Residual (used by FIC shock-capturing term)
    array_1d<double,3> MomRes = ZeroVector(3);
    this->AlgebraicMomentumResidual(rData,AGradN,MomRes);

    // Multiplying some quantities by density to have correct units
    body_force *= density; // Force per unit of volume
    AGradN *= density; // Convective term is always multiplied by density

    // Temporary containers
    double K,G,PDivV;

    // Note: Dof order is (u,v,[w,]p) for each node
    for (unsigned int i = 0; i < NumNodes; i++)
    {
        unsigned int row = i*BlockSize;

        // LHS terms
        for (unsigned int j = 0; j < NumNodes; j++)
        {
            unsigned int col = j*BlockSize;

            // Some terms are the same for all velocity components, calculate them once for each i,j
            //K = 0.5*(rData.N[i]*AGradN[j] - AGradN[i]*rData.N[j]); // Skew-symmetric convective term 1/2( v*grad(u)*u - grad(v) uu )
            K = rData.N[i]*AGradN[j];
            K += AGradN[i]*TauMomentum*(AGradN[j]); // Stabilization: u*grad(v) * TauOne * u*grad(u)
            K *= rData.Weight;

            // q-p stabilization block (reset result)
            double laplacian = 0;

            // The following lines implement the viscous term as a Laplacian
            //for (unsigned int d = 0; d < Dim; d++)
            //    K += rData.Weight * density * rData.EffectiveViscosity * rData.DN_DX(i, d) * rData.DN_DX(j, d);

            for (unsigned int d = 0; d < Dim; d++)
            {
                //K += rData.Weight * density * rData.EffectiveViscosity * rData.DN_DX(i, d) * rData.DN_DX(j, d);
                LHS(row+d,col+d) += K;

                // v * Grad(p) block
                G = AGradN[i] * rData.DN_DX(j,d); // Stabilization: (a * Grad(v)) * TauOne * Grad(p)
                PDivV = rData.DN_DX(i,d) * rData.N[j]; // Div(v) * p

                // Write v * Grad(p) component
                LHS(row+d,col+Dim) += rData.Weight * (TauMomentum*G - PDivV);
                // Use symmetry to write the q * Div(u) component
                LHS(col+Dim,row+d) += rData.Weight * (TauIncompr*G + PDivV);

                // q-p stabilization block
                laplacian += rData.DN_DX(i,d) * rData.DN_DX(j,d); // Stabilization: Grad(q) * TauOne * Grad(p)
            }

            // Write q-p term
            LHS(row+Dim,col+Dim) += rData.Weight*TauIncompr*laplacian;

            // FIC shock capturing term
            for (unsigned int d = 0; d < Dim; d++)
                LHS(row+d,col+d) += rData.Weight * density * std::fabs(TauGrad[d]*MomRes[d]) * laplacian;
        }

        // RHS terms
        double forcing = 0.0;
        for (unsigned int d = 0; d < Dim; ++d)
        {
            rLocalRHS[row+d] += rData.Weight * rData.N[i] * body_force[d]; // v*body_force
            rLocalRHS[row+d] += rData.Weight * TauMomentum * AGradN[i] * body_force[d]; // ( a * Grad(v) ) * TauOne * (Density * BodyForce)
            forcing += rData.DN_DX(i, d) * body_force[d];
        }
        rLocalRHS[row + Dim] += rData.Weight * TauIncompr * forcing; // Grad(q) * TauOne * (density * body_force)
    }

    // Write (the linearized part of the) local contribution into residual form (A*dx = b - A*x)
    array_1d<double,LocalSize> values;
    this->GetCurrentValuesVector(rData,values);
    noalias(rLocalRHS) -= prod(LHS, values);

    /* Viscous contribution (symmetric gradient E(u) - 1/3 Tr(E) )
     * For a generic (potentially non-linear) constitutive law, one cannot assume that RHS = F - LHS*current_values.
     * Because of this, the AddViscousTerm function manages both the LHS and the RHS.
     */
    this->AddViscousTerm(rData,LHS,rLocalRHS);

    noalias(rLocalLHS) += LHS;
}

///////////////////////////////////////////////////////////////////////////////////////////////////

template< class TElementData >
void FIC<TElementData>::AddMassLHS(
    TElementData& rData,
    MatrixType &rMassMatrix)
{
    const double density = this->GetAtCoordinate(rData.Density,rData.N);

    // Note: Dof order is (u,v,[w,]p) for each node
    for (unsigned int i = 0; i < NumNodes; i++)
    {
        unsigned int row = i*BlockSize;
        for (unsigned int j = 0; j < NumNodes; j++)
        {
            unsigned int col = j*BlockSize;
            const double Mij = rData.Weight * density * rData.N[i] * rData.N[j];
            for (unsigned int d = 0; d < Dim; d++)
                rMassMatrix(row+d,col+d) += Mij;
        }
    }

    /* NOTE: in FIC we have momentum stabilization terms on the mass matrix irregardless of OSS
        */
    this->AddMassStabilization(rData,rMassMatrix);

}

///////////////////////////////////////////////////////////////////////////////////////////////////

template< class TElementData >
void FIC<TElementData>::AddMassStabilization(
    TElementData& rData,
    MatrixType &rMassMatrix)
{
    const double density = this->GetAtCoordinate(rData.Density,rData.N);

    double TauIncompr;
    double TauMomentum;
    array_1d<double,3> convective_velocity = this->GetAtCoordinate(rData.Velocity,rData.N) - this->GetAtCoordinate(rData.MeshVelocity,rData.N);

    array_1d<double,3> TauGrad = ZeroVector(3);
    this->CalculateTau(rData,convective_velocity,TauIncompr,TauMomentum,TauGrad);

    Vector AGradN;
    this->ConvectionOperator(AGradN,convective_velocity,rData.DN_DX);

    // Multiplying some quantities by density to have correct units
    AGradN *= density; // Convective term is always multiplied by density

    // Temporary container
    double K;
    double weight = rData.Weight * density; // This density is for the dynamic term in the residual (rho*Du/Dt)

    // Note: Dof order is (u,v,[w,]p) for each node
    for (unsigned int i = 0; i < NumNodes; i++)
    {
        unsigned int row = i*BlockSize;

        for (unsigned int j = 0; j < NumNodes; j++)
        {
            unsigned int col = j*BlockSize;

            K = TauMomentum * weight * AGradN[i] * rData.N[j];

            for (unsigned int d = 0; d < Dim; d++)
            {
                rMassMatrix(row+d,col+d) += K;
                rMassMatrix(row+Dim,col+d) += TauIncompr*weight*rData.DN_DX(i,d)*rData.N[j];
            }
        }
    }
}

template <class TElementData>
void FIC<TElementData>::AddBoundaryTraction(TElementData& rData,
    const Vector& rUnitNormal, MatrixType& rLHS, VectorType& rRHS) {

    BoundedMatrix<double,StrainSize,LocalSize> strain_matrix = ZeroMatrix(StrainSize,LocalSize);
    FluidElementUtilities<NumNodes>::GetStrainMatrix(rData.DN_DX,strain_matrix);

    const auto& constitutive_matrix = rData.C;

    BoundedMatrix<double,StrainSize,LocalSize> shear_stress_matrix = prod(constitutive_matrix,strain_matrix);

    BoundedMatrix<double,Dim,StrainSize> normal_projection = ZeroMatrix(Dim,StrainSize);
    FluidElementUtilities<NumNodes>::VoigtTransformForProduct(rUnitNormal,normal_projection);

    // Contribution to boundary stress from 2*mu*symmetric_gradient(velocity)*n
    BoundedMatrix<double,Dim,LocalSize> normal_stress_operator = prod(normal_projection,shear_stress_matrix);

    // Contribution to boundary stress from p*n
    for (unsigned int i = 0; i < NumNodes; i++) {
        const double ni = rData.N[i];
        for (unsigned int d = 0; d < Dim; d++) {
            const std::size_t pressure_column = i*BlockSize + Dim;
            normal_stress_operator(d,pressure_column) = -rUnitNormal[d]*ni;
        }
    }

    // RHS: stress computed using current solution
    array_1d<double,Dim> shear_stress = prod(normal_projection,rData.ShearStress);
    const double p_gauss = this->GetAtCoordinate(rData.Pressure,rData.N);

    // Add -Ni*normal_stress_operator to the LHS, Ni*current_stress to the RHS
    for (unsigned int i = 0; i < NumNodes; i++) {
        const double wni = rData.Weight*rData.N[i];
        for (unsigned int d = 0; d < Dim; d++) {
            const unsigned int row = i*BlockSize + d;
            for (unsigned int col = 0; col < LocalSize; col++) {
                rLHS(row,col) += wni*normal_stress_operator(d,col);
            }
            rRHS[row] -= wni*(shear_stress[d]-p_gauss*rUnitNormal[d]);
        }
    }
}

template <class TElementData>
void FIC<TElementData>::AddViscousTerm(
    const TElementData& rData,
    BoundedMatrix<double,LocalSize,LocalSize>& rLHS,
    VectorType& rRHS) {

    BoundedMatrix<double,StrainSize,LocalSize> strain_matrix = ZeroMatrix(StrainSize,LocalSize);
    FluidElementUtilities<NumNodes>::GetStrainMatrix(rData.DN_DX,strain_matrix);

    const auto& constitutive_matrix = rData.C;
    BoundedMatrix<double,StrainSize,LocalSize> shear_stress_matrix = prod(constitutive_matrix,strain_matrix);

    // Multiply times integration point weight (I do this here to avoid a temporal in LHS += weight * Bt * C * B)
    strain_matrix *= rData.Weight;

    noalias(rLHS) += prod(trans(strain_matrix),shear_stress_matrix);
    noalias(rRHS) -= prod(trans(strain_matrix),rData.ShearStress);
}

///////////////////////////////////////////////////////////////////////////////////////////////////

template< class TElementData >
void FIC<TElementData>::CalculateTau(
    const TElementData& rData,
    const array_1d<double,3> &Velocity,
    double &TauIncompr,
    double &TauMomentum,
    array_1d<double,3> &TauGrad) const
{
    const Geometry< Node<3> >& r_geometry = this->GetGeometry();

    constexpr double c1 = 8.0;
    constexpr double c2 = 2.0;

    const double Beta = rData.FICBeta;
    const double Nobeta = 1.0-Beta;

    double Havg = ElementSizeCalculator<Dim,NumNodes>::AverageElementSize(r_geometry);

    double velocity_norm = Velocity[0]*Velocity[0];
    for (unsigned int d = 1; d < Dim; d++)
        velocity_norm += Velocity[d]*Velocity[d];
    velocity_norm = std::sqrt(velocity_norm);

    double Hvel = Havg;
    if (velocity_norm > 1.0e-6)
    {
        Hvel = ElementSizeCalculator<Dim,NumNodes>::ProjectedElementSize(r_geometry,Velocity);
    }

    const double density = this->GetAtCoordinate(rData.Density,rData.N);
    const double viscosity = this->GetAtCoordinate(rData.EffectiveViscosity,rData.N);

    double InvTau = c1 * viscosity / (Havg*Havg) + density * c2 * velocity_norm / Havg;
    TauIncompr = 1.0/InvTau;
    TauMomentum = (Hvel / (density * c2 * velocity_norm) );

    // TAU limiter for momentum equation: tau = min{ h/2u, dt }
    double TimeTerm = rData.DeltaTime/density;
    if (TauMomentum > TimeTerm)
    {
        TauMomentum = TimeTerm;
    }

    TauMomentum *= Beta;

    // Coefficients for FIC shock-capturing term
    this->CalculateTauGrad(rData,TauGrad);
    TauGrad /= density;
    for (unsigned int d = 0; d < Dim; d++)
        if (TauGrad[d] > Havg*TimeTerm)
            TauGrad[d] = Havg*TimeTerm;

    TauGrad *= Nobeta;
}

template< class TElementData >
void FIC<TElementData>::CalculateTauGrad(
    const TElementData& rData,
    array_1d<double,3> &TauGrad) const
{
    // Small constant to prevent division by zero
    const double Small = 1.0e-12;

    // Evaluate velocity gradient
    const auto& r_velocities = rData.Velocity;
    BoundedMatrix<double,3,3> Gradient = ZeroMatrix(3,3);
    for (unsigned int n = 0; n < NumNodes; n++)
    {
        for (unsigned int i = 0; i < Dim; i++)
        {
            for (unsigned int j = 0; j < Dim; j++)
            {
                Gradient(i,j) += rData.DN_DX(n,j)*r_velocities(n,i);
            }
        }
    }

    // Calculate characteristic lenghts on the gradient directions and gradient norms
    const Geometry< Node<3> >& r_geometry = this->GetGeometry();
    array_1d<double,3> Hg = ZeroVector(3);
    array_1d<double,3> GradNorm;
    for (unsigned int d = 0; d < Dim; d++)
    {
        array_1d<double,3> Gi;
        Gi[0] = Gradient(d,0);
        Gi[1] = Gradient(d,1);
        Gi[2] = Gradient(d,2);

        Hg[d] = ElementSizeCalculator<Dim,NumNodes>::ProjectedElementSize(r_geometry,Gi);

        GradNorm[d] = std::sqrt(Gi[0]*Gi[0]+Gi[1]*Gi[1]+Gi[2]*Gi[2]);

        // Coefficients for shock capturing term (remember to multiply by momentum residual!)
        TauGrad[d] = Hg[d] / (2.*GradNorm[d] + Small);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Private functions
///////////////////////////////////////////////////////////////////////////////////////////////////

// serializer

template< class TElementData >
void FIC<TElementData>::save(Serializer& rSerializer) const
{
    typedef FluidElement<TElementData> BaseElement;
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseElement );
}


template< class TElementData >
void FIC<TElementData>::load(Serializer& rSerializer)
{
    typedef FluidElement<TElementData> BaseElement;
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseElement);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Class template instantiation

template class FIC< FICData<2,3> >;
template class FIC< FICData<3,4> >;

template class FIC< FICData<2,4> >;
template class FIC< FICData<3,8> >;

template class FIC< TimeIntegratedFICData<2,3> >;
template class FIC< TimeIntegratedFICData<3,4> >;

} // namespace Kratos
