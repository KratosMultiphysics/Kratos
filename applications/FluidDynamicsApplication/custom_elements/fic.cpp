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
//

#include "fic.h"
#include "includes/cfd_variables.h"
#include "includes/checks.h"

#include "custom_utilities/qsvms_data.h"
#include "custom_utilities/time_integrated_qsvms_data.h"
#include "custom_utilities/fluid_element_utilities.h"
#include "custom_utilities/element_size_calculator.h"

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
    array_1d<double, 3>& rOutput, const ProcessInfo& rCurrentProcessInfo) {
    // Lumped projection terms
    if (rVariable == ADVPROJ) {
        this->CalculateProjections(rCurrentProcessInfo);
    }
}

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
    KRATOS_CHECK_VARIABLE_KEY(NODAL_AREA);

    for(unsigned int i=0; i<NumNodes; ++i)
    {
        Node<3>& rNode = this->GetGeometry()[i];
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(ACCELERATION,rNode);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(NODAL_AREA,rNode);
    }

    return out;
}

///////////////////////////////////////////////////////////////////////////////////////////////////

template< class TElementData >
void FIC<TElementData>::GetValueOnIntegrationPoints(Variable<array_1d<double, 3 > > const& rVariable,
                                            std::vector<array_1d<double, 3 > >& rValues,
                                            ProcessInfo const& rCurrentProcessInfo)
{
    if (rVariable == VORTICITY)
    {
        // Get Shape function data
        Vector GaussWeights;
        Matrix ShapeFunctions;
        ShapeFunctionDerivativesArrayType ShapeDerivatives;
        this->CalculateGeometryData(GaussWeights,ShapeFunctions,ShapeDerivatives);
        const unsigned int NumGauss = GaussWeights.size();

        rValues.resize(NumGauss);

        for (unsigned int g = 0; g < NumGauss; g++)
        {
            this->IntegrationPointVorticity(ShapeDerivatives[g],rValues[g]);
        }
    }
}


template< class TElementData >
void FIC<TElementData>::GetValueOnIntegrationPoints(Variable<double> const& rVariable,
                                            std::vector<double>& rValues,
                                            ProcessInfo const& rCurrentProcessInfo)
{
    FluidElement<TElementData>::GetValueOnIntegrationPoints(rVariable,rValues,rCurrentProcessInfo);
}

template <class TElementData>
void FIC<TElementData>::GetValueOnIntegrationPoints(Variable<array_1d<double, 6>> const& rVariable,
                                                    std::vector<array_1d<double, 6>>& rValues,
                                                    ProcessInfo const& rCurrentProcessInfo)
{
    FluidElement<TElementData>::GetValueOnIntegrationPoints(rVariable,rValues,rCurrentProcessInfo);
}

template <class TElementData>
void FIC<TElementData>::GetValueOnIntegrationPoints(Variable<Vector> const& rVariable,
                                                    std::vector<Vector>& rValues,
                                                    ProcessInfo const& rCurrentProcessInfo)
{
    FluidElement<TElementData>::GetValueOnIntegrationPoints(rVariable,rValues,rCurrentProcessInfo);
}

template <class TElementData>
void FIC<TElementData>::GetValueOnIntegrationPoints(Variable<Matrix> const& rVariable,
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

///////////////////////////////////////////////////////////////////////////////////////////////////
// Evaluation of system terms on Gauss Points

template <class TElementData>
void FIC<TElementData>::AddTimeIntegratedSystem(
    TElementData& rData, MatrixType& rLHS, VectorType& rRHS) {

    // Call specialized implementation (it is on a helper class to avoid partial template specialization problems)
    Internals::SpecializedAddTimeIntegratedSystem<TElementData,
        TElementData::ElementManagesTimeIntegration>::AddSystem(this, rData,
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
    double density = rData.Density;
    array_1d<double,3> body_force = this->Interpolate(rData.BodyForce,rData.N);
    array_1d<double,3> momentum_projection = this->Interpolate(rData.MomentumProjection,rData.N);

    array_1d<double,3> convective_velocity = 
        this->Interpolate(rData.Velocity,rData.N) - 
        this->Interpolate(rData.MeshVelocity,rData.N);

    double TauIncompr;
    double TauMomentum;
    array_1d<double,3> TauGrad(3,0.0);
    this->CalculateTau(rData,convective_velocity,TauIncompr,TauMomentum,TauGrad);

//TODO: seguir

    Vector AGradN;
    this->ConvectionOperator(AGradN,convective_velocity,rDN_DX);

    // Residual (used by FIC shock-capturing term)
    array_1d<double,3> MomRes(3,0.0);
    this->ASGSMomentumResidual(GaussIndex,rN,rDN_DX,MomRes);

    // Multiplying some quantities by density to have correct units
    Viscosity *= density; // Dynamic viscosity
    body_force *= density; // Force per unit of volume
    AGradN *= density; // Convective term is always multiplied by density

    // Auxiliary variables for matrix looping
    const unsigned int NumNodes = rN.size();
    const unsigned int BlockSize = TDim+1;
    unsigned int Row = 0;
    unsigned int Col = 0;

    // Temporary containers
    double K,L,G,PDivV,qF;


    // Note: Dof order is (u,v,[w,]p) for each node
    for (unsigned int i = 0; i < NumNodes; i++)
    {
        Row = i*BlockSize;

        // LHS terms
        for (unsigned int j = 0; j < NumNodes; j++)
        {
            Col = j*BlockSize;

            // Some terms are the same for all velocity components, calculate them once for each i,j
            K = 0.5*(rN[i]*AGradN[j] - AGradN[i]*rN[j]); // Skew-symmetric convective term 1/2( v*grad(u)*u - grad(v) uu )
            K += AGradN[i]*TauMomentum*(AGradN[j]); // Stabilization: u*grad(v) * TauOne * u*grad(u)
            K *= GaussWeight;

            // q-p stabilization block (reset result)
            L = 0;

            // The following lines implement the viscous term as a Laplacian
            //for (unsigned int d = 0; d < TDim; d++)
            //    K += GaussWeight * density * Viscosity * rDN_DX(i, d) * rDN_DX(j, d);

            for (unsigned int d = 0; d < TDim; d++)
            {
                //K += GaussWeight * density * Viscosity * rDN_DX(i, d) * rDN_DX(j, d);
                rLHS(Row+d,Col+d) += K;

                // v * Grad(p) block
                G = AGradN[i] * rDN_DX(j,d); // Stabilization: (a * Grad(v)) * TauOne * Grad(p)
                PDivV = rDN_DX(i,d) * rN[j]; // Div(v) * p

                // Write v * Grad(p) component
                rLHS(Row+d,Col+TDim) += GaussWeight * (TauMomentum*G - PDivV);
                // Use symmetry to write the q * Div(u) component
                rLHS(Col+TDim,Row+d) += GaussWeight * (TauIncompr*G + PDivV);

                // q-p stabilization block
                L += rDN_DX(i,d) * rDN_DX(j,d); // Stabilization: Grad(q) * TauOne * Grad(p)
            }

            // Write q-p term
            rLHS(Row+TDim,Col+TDim) += GaussWeight*TauIncompr*L;

            // FIC shock capturing term
            for (unsigned int d = 0; d < TDim; d++)
                rLHS(Row+d,Col+d) += GaussWeight * density * std::fabs(TauGrad[d]*MomRes[d]) * L;
        }

        // RHS terms
        qF = 0.0;
        for (unsigned int d = 0; d < TDim; ++d)
        {
            rRHS[Row+d] += GaussWeight * rN[i] * body_force[d]; // v*body_force
            rRHS[Row+d] += GaussWeight * TauMomentum * AGradN[i] * body_force[d]; // ( a * Grad(v) ) * TauOne * (Density * BodyForce)
            qF += rDN_DX(i, d) * (body_force[d] - momentum_projection[d]);
        }
        rRHS[Row + TDim] += GaussWeight * TauIncompr * qF; // Grad(q) * TauOne * (density * body_force)
    }

    // Viscous contribution (symmetric gradient E(u) - 1/3 Tr(E) )
    // This could potentially be optimized, as it can be integrated exactly using one less integration order when compared to previous terms.
    this->AddViscousTerm(Viscosity,GaussWeight,rDN_DX,rLHS);

}

///////////////////////////////////////////////////////////////////////////////////////////////////

template< class TElementData >
void FIC<TElementData>::AddMassLHS(
    TElementData& rData,
    MatrixType &rMassMatrix)
{
    double density = rData.Density;

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
    double density = rData.Density;

    array_1d<double,3> convective_velocity = 
        this->Interpolate(rData.Velocity,rData.N) - 
        this->Interpolate(rData.MeshVelocity,rData.N);

    double TauIncompr;
    double TauMomentum;
    array_1d<double,3> TauGrad(3,0.0);
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
                if ( rData.UseOSS != 1.0 )
                    rMassMatrix(row+Dim,col+d) += TauIncompr*weight*rData.DN_DX(i,d)*rData.N[j];
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////

template <class TElementData>
void FIC<TElementData>::AddBoundaryIntegral(TElementData& rData,
    const Vector& rUnitNormal, MatrixType& rLHS, VectorType& rRHS) {

}

///////////////////////////////////////////////////////////////////////////////////////////////////

template< class TElementData >
void FIC<TElementData>::CalculateTau(
    const TElementData& rData,
    const array_1d<double,3> &Velocity,
    double &TauIncompr,
    double &TauMomentum,
    array_1d<double,3> &TauGrad)
{
    GeometryType& rGeom = this->GetGeometry();

    constexpr double c1 = 8.0;
    constexpr double c2 = 2.0;

    const double Beta = rData.FICBeta;
    const double Nobeta = 1.0-Beta;

    double Havg = ElementSizeCalculator<Dim,NumNodes>::AverageElementSize(rGeom);

    double velocity_norm = Velocity[0]*Velocity[0];
    for (unsigned int d = 1; d < Dim; d++)
        velocity_norm += Velocity[d]*Velocity[d];
    velocity_norm = std::sqrt(velocity_norm);

    double Hvel = Havg;
    if (velocity_norm > 1.0e-6)
    {
        Hvel = ElementSizeCalculator<Dim,NumNodes>::ProjectedElementSize(rGeom,Velocity);
    }

    double InvTau = rData.Density * ( c1 * rData.EffectiveViscosity / (Havg*Havg) + c2 * velocity_norm / Havg );
    TauIncompr = 1.0/InvTau;
    TauMomentum = (Hvel / (rData.Density * c2 * velocity_norm) );

    // TAU limiter for momentum equation: tau = min{ h/2u, dt }
    double TimeTerm = rData.DeltaTime/rData.Density;
    if (TauMomentum > TimeTerm)
    {
        TauMomentum = TimeTerm;
    }

    TauMomentum *= Beta;

    // Coefficients for FIC shock-capturing term
    this->CalculateTauGrad(TauGrad);
    TauGrad /= rData.Density;
    for (unsigned int d = 0; d < Dim; d++)
        if (TauGrad[d] > Havg*TimeTerm)
            TauGrad[d] = Havg*TimeTerm;

    TauGrad *= Nobeta;
}

template< class TElementData >
void FIC<TElementData>::CalculateTauGrad(array_1d<double,3> &TauGrad)
{
    // Small constant to prevent division by zero
    const double Small = 1.0e-12;

    // Evaluate velocity gradient
    GeometryType& rGeom = this->GetGeometry();

    ShapeFunctionDerivativesArrayType DN_DX;
    rGeom.ShapeFunctionsIntegrationPointsGradients(DN_DX,GeometryData::GI_GAUSS_1);
    ShapeFunctionDerivativesType& rDN_DX = DN_DX[0];

    boost::numeric::ublas::bounded_matrix<double,3,3> Gradient = ZeroMatrix(3,3);
    for (unsigned int n = 0; n < NumNodes; n++)
    {
        const array_1d<double,3>& rU = rGeom[n].FastGetSolutionStepValue(VELOCITY);
        for (unsigned int i = 0; i < Dim; i++)
        {
            for (unsigned int j = 0; j < Dim; j++)
            {
                Gradient(i,j) += rDN_DX(n,j)*rU[i];
            }
        }
    }

    // Calculate characteristic lenghts on the gradient directions and gradient norms
    array_1d<double,3> Hg(3,0.0);
    array_1d<double,3> GradNorm(3,0.0);
    for (unsigned int d = 0; d < Dim; d++)
    {
        array_1d<double,3> Gi(3,0.0);
        Gi[0] = Gradient(d,0);
        Gi[1] = Gradient(d,1);
        Gi[2] = Gradient(d,2);

        Hg[d] = ElementSizeCalculator<Dim,NumNodes>::ProjectedElementSize(rGeom,Gi);

        GradNorm[d] = std::sqrt(Gi[0]*Gi[0]+Gi[1]*Gi[1]+Gi[2]*Gi[2]);

        // Coefficients for shock capturing term (remember to multiply by momentum residual!)
        TauGrad[d] = Hg[d] / (2.*GradNorm[d] + Small);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////

template< class TElementData >
void FIC<TElementData>::CalculateProjections(const ProcessInfo &rCurrentProcessInfo)
{
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
// Internals
///////////////////////////////////////////////////////////////////////////////////////////////////
namespace Internals {

///////////////////////////////////////////////////////////////////////////////////////////////////
// For Dim == 2
///////////////////////////////////////////////////////////////////////////////////////////////////

template <>
void AddViscousTerm<2>(double DynamicViscosity,
                       double GaussWeight,
                       const Kratos::Matrix& rDN_DX,
                       Kratos::Matrix& rLHS)
{
    double weight = GaussWeight * DynamicViscosity;

    constexpr double four_thirds = 4.0 / 3.0;
    constexpr double minus_two_thirds = -2.0 / 3.0;

    const unsigned int num_nodes = rDN_DX.size1();
    const unsigned int block_size = 3;

    for (unsigned int a = 0; a < num_nodes; ++a)
    {
        unsigned int row = a*block_size;
        for (unsigned int b = 0; b < num_nodes; ++b)
        {
            unsigned int col = b*block_size;

            // First row
            rLHS(row,col) += weight * ( four_thirds * rDN_DX(a,0) * rDN_DX(b,0) + rDN_DX(a,1) * rDN_DX(b,1) );
            rLHS(row,col+1) += weight * ( minus_two_thirds * rDN_DX(a,0) * rDN_DX(b,1) + rDN_DX(a,1) * rDN_DX(b,0) );

            // Second row
            rLHS(row+1,col) += weight * ( minus_two_thirds * rDN_DX(a,1) * rDN_DX(b,0) + rDN_DX(a,0) * rDN_DX(b,1) );
            rLHS(row+1,col+1) += weight * ( four_thirds * rDN_DX(a,1) * rDN_DX(b,1) + rDN_DX(a,0) * rDN_DX(b,0) );
        }
    }
}


///////////////////////////////////////////////////////////////////////////////////////////////////
// For Dim == 3
///////////////////////////////////////////////////////////////////////////////////////////////////

template <>
void AddViscousTerm<3>(double DynamicViscosity,
                       double GaussWeight,
                       const Kratos::Matrix& rDN_DX,
                       Kratos::Matrix& rLHS)
{
    double weight = GaussWeight * DynamicViscosity;

    constexpr double one_third = 1.0 / 3.0;
    constexpr double minus_two_thirds = -2.0 / 3.0;

    const unsigned int num_nodes = rDN_DX.size1();
    const unsigned int block_size = 4;

    unsigned int row(0),col(0);

    for (unsigned int i = 0; i < num_nodes; ++i)
    {
        row = i*block_size;
        for (unsigned int j = 0; j < num_nodes; ++j)
        {
            col = j*block_size;
            // (dN_i/dx_k dN_j/dx_k)
            const double diag =  rDN_DX(i,0) * rDN_DX(j,0) + rDN_DX(i,1) * rDN_DX(j,1) + rDN_DX(i,2) * rDN_DX(j,2);

            // First row
            rLHS(row,col) += weight * ( one_third * rDN_DX(i,0) * rDN_DX(j,0) + diag );
            rLHS(row,col+1) += weight * ( minus_two_thirds * rDN_DX(i,0) * rDN_DX(j,1) + rDN_DX(i,1) * rDN_DX(j,0) );
            rLHS(row,col+2) += weight * ( minus_two_thirds * rDN_DX(i,0) * rDN_DX(j,2) + rDN_DX(i,2) * rDN_DX(j,0) );

            // Second row
            rLHS(row+1,col) += weight * ( minus_two_thirds * rDN_DX(i,1) * rDN_DX(j,0) + rDN_DX(i,0) * rDN_DX(j,1) );
            rLHS(row+1,col+1) += weight * ( one_third * rDN_DX(i,1) * rDN_DX(j,1) + diag );
            rLHS(row+1,col+2) += weight * ( minus_two_thirds * rDN_DX(i,1) * rDN_DX(j,2) + rDN_DX(i,2) * rDN_DX(j,1) );

            // Third row
            rLHS(row+2,col) += weight * ( minus_two_thirds * rDN_DX(i,2) * rDN_DX(j,0) + rDN_DX(i,0) * rDN_DX(j,2) );
            rLHS(row+2,col+1) += weight * ( minus_two_thirds * rDN_DX(i,2) * rDN_DX(j,1) + rDN_DX(i,1) * rDN_DX(j,2) );
            rLHS(row+2,col+2) += weight * ( one_third * rDN_DX(i,2) * rDN_DX(j,2) + diag );
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// For Standard data: Time integration is not available
///////////////////////////////////////////////////////////////////////////////////////////////////

template <class TElementData>
void SpecializedAddTimeIntegratedSystem<TElementData, false>::AddSystem(
    FIC<TElementData>* pElement, TElementData& rData, Matrix& rLHS,
    Vector& rRHS) {
    KRATOS_TRY;
    KRATOS_ERROR << "Trying to use time-integrated element functions with a "
                    "data type that does not know previous time step data"
                 << std::endl;
    KRATOS_CATCH("");
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Specialized time integration
///////////////////////////////////////////////////////////////////////////////////////////////////

template <class TElementData>
void SpecializedAddTimeIntegratedSystem<TElementData, true>::AddSystem(
    FIC<TElementData>* pElement, TElementData& rData, Matrix& rLHS,
    Vector& rRHS) {
        Matrix mass_matrix = ZeroMatrix(rLHS.size1(),rLHS.size2());
        Matrix velocity_lhs = ZeroMatrix(rLHS.size1(),rLHS.size2());

        pElement->AddVelocitySystem(rData,velocity_lhs,rRHS);
        pElement->AddMassLHS(rData,mass_matrix);

        noalias(rLHS) += rData.bdf0*mass_matrix + velocity_lhs;
        
        Vector values = ZeroVector(rRHS.size());
        Vector acceleration = ZeroVector(rRHS.size());

        int LocalIndex = 0;
        const auto& r_velocities = rData.Velocity;
        const auto& r_velocities_step1 = rData.Velocity_OldStep1;
        const auto& r_velocities_step2 = rData.Velocity_OldStep2;
        const auto& r_pressures = rData.Pressure;

        for (unsigned int i = 0; i < TElementData::NumNodes; ++i) {
            for (unsigned int d = 0; d < TElementData::Dim; ++d)  {
                values[LocalIndex] = r_velocities(i,d);
                // Velocity Dofs
                acceleration[LocalIndex] = rData.bdf0*r_velocities(i,d);
                acceleration[LocalIndex] += rData.bdf1*r_velocities_step1(i,d);
                acceleration[LocalIndex] += rData.bdf2*r_velocities_step2(i,d);
                ++LocalIndex;
            }
            values[LocalIndex] = r_pressures[i];
            ++LocalIndex;
        }

        noalias(rRHS) -= prod(velocity_lhs,values); //TODO: should this be here? (Ignasi)
        noalias(rRHS) -= prod(mass_matrix,acceleration);
}

} // namespace Internals

///////////////////////////////////////////////////////////////////////////////////////////////////
// Class template instantiation

template class FIC< FICData<2,3> >;
template class FIC< FICData<3,4> >;

template class FIC< FICData<2,4> >;
template class FIC< FICData<3,8> >;

} // namespace Kratos
