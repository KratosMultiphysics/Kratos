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

#include "qs_vms.h"
#include "includes/cfd_variables.h"
#include "includes/checks.h"

#include "custom_utilities/qsvms_data.h"
#include "custom_utilities/time_integrated_qsvms_data.h"
#include "custom_utilities/qsvms_dem_coupled_data.h"
#include "custom_utilities/fluid_element_utilities.h"
#include "custom_utilities/fluid_element_time_integration_detail.h"

namespace Kratos
{

///////////////////////////////////////////////////////////////////////////////////////////////////
// Life cycle

template< class TElementData >
QSVMS<TElementData>::QSVMS(IndexType NewId):
    FluidElement<TElementData>(NewId)
{}

template< class TElementData >
QSVMS<TElementData>::QSVMS(IndexType NewId, const NodesArrayType& ThisNodes):
    FluidElement<TElementData>(NewId,ThisNodes)
{}


template< class TElementData >
QSVMS<TElementData>::QSVMS(IndexType NewId, GeometryType::Pointer pGeometry):
    FluidElement<TElementData>(NewId,pGeometry)
{}


template< class TElementData >
QSVMS<TElementData>::QSVMS(IndexType NewId, GeometryType::Pointer pGeometry, Properties::Pointer pProperties):
    FluidElement<TElementData>(NewId,pGeometry,pProperties)
{}


template< class TElementData >
QSVMS<TElementData>::~QSVMS()
{}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Public Operations

template< class TElementData >
Element::Pointer QSVMS<TElementData>::Create(IndexType NewId,NodesArrayType const& ThisNodes,Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<QSVMS>(NewId, this->GetGeometry().Create(ThisNodes), pProperties);
}


template< class TElementData >
Element::Pointer QSVMS<TElementData>::Create(IndexType NewId,GeometryType::Pointer pGeom,Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<QSVMS>(NewId, pGeom, pProperties);
}

template <class TElementData>
void QSVMS<TElementData>::Calculate(
    const Variable<double>& rVariable,
    double& rOutput,
    const ProcessInfo& rCurrentProcessInfo)
{
    BaseType::Calculate(rVariable, rOutput, rCurrentProcessInfo);
}

template <class TElementData>
void QSVMS<TElementData>::Calculate(
    const Variable<array_1d<double, 3>>& rVariable,
    array_1d<double, 3>& rOutput, const ProcessInfo& rCurrentProcessInfo) {
    // Lumped projection terms
    if (rVariable == ADVPROJ) {
        this->CalculateProjections(rCurrentProcessInfo);
    } else {
        BaseType::Calculate(rVariable, rOutput, rCurrentProcessInfo);
    }
}

template <class TElementData>
void QSVMS<TElementData>::Calculate(
    const Variable<Vector>& rVariable,
    Vector& rOutput,
    const ProcessInfo& rCurrentProcessInfo)
{
    BaseType::Calculate(rVariable, rOutput, rCurrentProcessInfo);
}

template <class TElementData>
void QSVMS<TElementData>::Calculate(
    const Variable<Matrix>& rVariable,
    Matrix& rOutput,
    const ProcessInfo& rCurrentProcessInfo)
{
    BaseType::Calculate(rVariable, rOutput, rCurrentProcessInfo);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Inquiry

template< class TElementData >
int QSVMS<TElementData>::Check(const ProcessInfo &rCurrentProcessInfo) const
{
    int out = FluidElement<TElementData>::Check(rCurrentProcessInfo);
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
void QSVMS<TElementData>::CalculateOnIntegrationPoints(
    Variable<array_1d<double, 3 > > const& rVariable,
    std::vector<array_1d<double, 3 > >& rValues,
    ProcessInfo const& rCurrentProcessInfo)
{
    if (rVariable == SUBSCALE_VELOCITY) {
        // Get Shape function data
        Vector GaussWeights;
        Matrix ShapeFunctions;
        ShapeFunctionDerivativesArrayType ShapeDerivatives;
        this->CalculateGeometryData(GaussWeights,ShapeFunctions,ShapeDerivatives);
        const unsigned int NumGauss = GaussWeights.size();

        rValues.resize(NumGauss);

        TElementData data;
        data.Initialize(*this, rCurrentProcessInfo);

        for (unsigned int g = 0; g < NumGauss; g++)
        {
            this->UpdateIntegrationPointData(data, g, GaussWeights[g], row(ShapeFunctions, g), ShapeDerivatives[g]);

            this->SubscaleVelocity(data, rValues[g]);
        }
    }
    else {
        FluidElement<TElementData>::CalculateOnIntegrationPoints(rVariable,rValues,rCurrentProcessInfo);
    }
}


template< class TElementData >
void QSVMS<TElementData>::CalculateOnIntegrationPoints(
    Variable<double> const& rVariable,
    std::vector<double>& rValues,
    ProcessInfo const& rCurrentProcessInfo)
{
    if (rVariable == SUBSCALE_PRESSURE) {
        // Get Shape function data
        Vector GaussWeights;
        Matrix ShapeFunctions;
        ShapeFunctionDerivativesArrayType ShapeDerivatives;
        this->CalculateGeometryData(GaussWeights,ShapeFunctions,ShapeDerivatives);
        const unsigned int NumGauss = GaussWeights.size();

        rValues.resize(NumGauss);

        TElementData data;
        data.Initialize(*this, rCurrentProcessInfo);

        for (unsigned int g = 0; g < NumGauss; g++)
        {
            this->UpdateIntegrationPointData(data, g, GaussWeights[g], row(ShapeFunctions, g), ShapeDerivatives[g]);

            this->SubscalePressure(data,rValues[g]);
        }

    }
    else {
        FluidElement<TElementData>::CalculateOnIntegrationPoints(rVariable,rValues,rCurrentProcessInfo);
    }
}

template <class TElementData>
void QSVMS<TElementData>::CalculateOnIntegrationPoints(
    Variable<array_1d<double, 6>> const& rVariable,
    std::vector<array_1d<double, 6>>& rValues,
    ProcessInfo const& rCurrentProcessInfo)
{
    FluidElement<TElementData>::CalculateOnIntegrationPoints(rVariable,rValues,rCurrentProcessInfo);
}

template <class TElementData>
void QSVMS<TElementData>::CalculateOnIntegrationPoints(
    Variable<Vector> const& rVariable,
    std::vector<Vector>& rValues,
    ProcessInfo const& rCurrentProcessInfo)
{
    FluidElement<TElementData>::CalculateOnIntegrationPoints(rVariable,rValues,rCurrentProcessInfo);
}

template <class TElementData>
void QSVMS<TElementData>::CalculateOnIntegrationPoints(
    Variable<Matrix> const& rVariable,
    std::vector<Matrix>& rValues,
    ProcessInfo const& rCurrentProcessInfo)
{
    FluidElement<TElementData>::CalculateOnIntegrationPoints(rVariable,rValues,rCurrentProcessInfo);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Input and output

template<class TElementData>
const Parameters QSVMS<TElementData>::GetSpecifications() const
{
    const Parameters specifications = Parameters(R"({
        "time_integration"           : ["implicit"],
        "framework"                  : "ale",
        "symmetric_lhs"              : false,
        "positive_definite_lhs"      : true,
        "output"                     : {
            "gauss_point"            : ["SUBSCALE_VELOCITY","SUBSCALE_PRESSURE","VORTICITY","Q_VALUE","VORTICITY_MAGNITUDE"],
            "nodal_historical"       : ["VELOCITY","PRESSURE"],
            "nodal_non_historical"   : [],
            "entity"                 : ["ADVPROJ"]
        },
        "required_variables"         : ["VELOCITY","ACCELERATION","MESH_VELOCITY","PRESSURE","IS_STRUCTURE","DISPLACEMENT","BODY_FORCE","NODAL_AREA","NODAL_H","ADVPROJ","DIVPROJ","REACTION","REACTION_WATER_PRESSURE","EXTERNAL_PRESSURE","NORMAL","Y_WALL","Q_VALUE"]
        "required_dofs"              : [],
        "flags_used"                 : [],
        "compatible_geometries"      : ["Triangle2D3","Quadrilateral2D4","Tetrahedra3D4","Hexahedra3D8"],
        "element_integrates_in_time" : false,
        "required_polynomial_degree_of_geometry" : 1,
        "documentation"   : "This implements a Navier-Stokes element with quasi-static Variational MultiScales (VMS) stabilization."
    })");

    if (Dim == 2) {
        std::vector<std::string> dofs_2d({"VELOCITY_X","VELOCITY_Y","PRESSURE"});
        specifications["required_dofs"].SetStringArray(dofs_2d);
    } else {
        std::vector<std::string> dofs_3d({"VELOCITY_X","VELOCITY_Y","VELOCITY_Z","PRESSURE"});
        specifications["required_dofs"].SetStringArray(dofs_3d);
    }

    return specifications;
}

template< class TElementData >
std::string QSVMS<TElementData>::Info() const
{
    std::stringstream buffer;
    buffer << "QSVMS #" << this->Id();
    return buffer.str();
}


template< class TElementData >
void QSVMS<TElementData>::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "QSVMS" << Dim << "D";
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Protected functions

template< class TElementData >
void QSVMS<TElementData>::AlgebraicMomentumResidual(
    const TElementData& rData,
    const array_1d<double,3> &rConvectionVelocity,
    array_1d<double,3>& rResidual) const
{
    const GeometryType rGeom = this->GetGeometry();

    Vector convection; // u * grad(N)
    this->ConvectionOperator(convection,rConvectionVelocity,rData.DN_DX);

    const double density = this->GetAtCoordinate(rData.Density,rData.N);
    const auto& r_body_forces = rData.BodyForce;
    const auto& r_velocities = rData.Velocity;
    const auto& r_pressures = rData.Pressure;

    for (unsigned int i = 0; i < NumNodes; i++) {
        const array_1d<double,3>& r_acceleration = rGeom[i].FastGetSolutionStepValue(ACCELERATION);
        for (unsigned int d = 0; d < Dim; d++) {
            rResidual[d] += density * ( rData.N[i]*(r_body_forces(i,d) - r_acceleration[d]) - convection[i]*r_velocities(i,d)) - rData.DN_DX(i,d)*r_pressures[i];
        }
    }
}

template< class TElementData >
void QSVMS<TElementData>::AlgebraicMassResidual(
    const TElementData& rData,
    double &rMomentumRes) const
{
    this->MassProjTerm(rData,rMomentumRes);
}

template< class TElementData >
void QSVMS<TElementData>::OrthogonalMomentumResidual(
    const TElementData& rData,
    const array_1d<double,3> &rConvectionVelocity,
    array_1d<double,3>& rResidual) const
{
    this->MomentumProjTerm(rData,rConvectionVelocity,rResidual);

    const array_1d<double,3> momentum_projection = this->GetAtCoordinate(rData.MomentumProjection,rData.N);
    noalias(rResidual) -= momentum_projection;
}

template< class TElementData >
void QSVMS<TElementData>::OrthogonalMassResidual(
    const TElementData& rData,
    double &rMassRes) const
{
    this->MassProjTerm(rData,rMassRes);
    double mass_projection = this->GetAtCoordinate(rData.MassProjection,rData.N);
    rMassRes -= mass_projection;
}


template< class TElementData >
void QSVMS<TElementData>::MomentumProjTerm(
    const TElementData& rData,
    const array_1d<double,3>& rConvectionVelocity,
    array_1d<double,3> &rMomentumRHS) const
{
    Vector AGradN;
    this->ConvectionOperator(AGradN,rConvectionVelocity,rData.DN_DX);

    const double density = this->GetAtCoordinate(rData.Density,rData.N);

    for (unsigned int i = 0; i < NumNodes; i++) {
        for (unsigned int d = 0; d < Dim; d++) {
            rMomentumRHS[d] += density * ( rData.N[i]*(rData.BodyForce(i,d) /*- rAcc[d]*/) - AGradN[i]*rData.Velocity(i,d)) - rData.DN_DX(i,d)*rData.Pressure[i];
        }
    }
}

template< class TElementData >
void QSVMS<TElementData>::MassProjTerm(
    const TElementData& rData,
    double &rMassRHS) const
{
    for (unsigned int i = 0; i < NumNodes; i++) {
        for (unsigned int d = 0; d < Dim; d++)
            rMassRHS -= rData.DN_DX(i,d)*rData.Velocity(i,d);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Evaluation of system terms on Gauss Points

template <class TElementData>
void QSVMS<TElementData>::AddTimeIntegratedSystem(
    TElementData& rData, MatrixType& rLHS, VectorType& rRHS) {

    // Call specialized implementation (it is on a helper class to avoid partial template specialization problems)
    Internals::FluidElementTimeIntegrationDetail<TElementData,
        TElementData::ElementManagesTimeIntegration>::AddTimeIntegratedSystem(this, rData,
        rLHS, rRHS);
}

template <class TElementData>
void QSVMS<TElementData>::AddTimeIntegratedLHS(
    TElementData& rData, MatrixType& rLHS) {

    KRATOS_ERROR << "AddTimeIntegratedLHS is not implemented." << std::endl;
}

template <class TElementData>
void QSVMS<TElementData>::AddTimeIntegratedRHS(
    TElementData& rData, VectorType& rRHS) {

    KRATOS_ERROR << "AddTimeIntegratedRHS is not implemented." << std::endl;
}

template< class TElementData >
void QSVMS<TElementData>::AddVelocitySystem(
    TElementData& rData,
    MatrixType &rLocalLHS,
    VectorType &rLocalRHS)
{
    auto& LHS = rData.LHS;
    LHS.clear();

    // Interpolate nodal data on the integration point
    const double density = this->GetAtCoordinate(rData.Density,rData.N);
    array_1d<double,3> body_force = this->GetAtCoordinate(rData.BodyForce,rData.N);
    array_1d<double,3> momentum_projection = this->GetAtCoordinate(rData.MomentumProjection,rData.N);
    double mass_projection = this->GetAtCoordinate(rData.MassProjection,rData.N);

    double tau_one;
    double tau_two;
    array_1d<double, 3> convective_velocity =
        this->GetAtCoordinate(rData.Velocity, rData.N) -
        this->GetAtCoordinate(rData.MeshVelocity, rData.N);

    this->CalculateTau(rData,convective_velocity,tau_one,tau_two);

    Vector AGradN;
    this->ConvectionOperator(AGradN,convective_velocity,rData.DN_DX);

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
            //K = 0.5*(rN[i]*AGradN[j] - AGradN[i]*rN[j]); // Skew-symmetric convective term 1/2( v*grad(u)*u - grad(v) uu )
            K = rData.N[i]*AGradN[j];
            K += AGradN[i]*tau_one*(AGradN[j]); // Stabilization: u*grad(v) * tau_one * u*grad(u)
            K *= rData.Weight;

            // q-p stabilization block (initialize result)
            double laplacian = 0;

            // The following lines implement the viscous term as a Laplacian
            //for (unsigned int d = 0; d < Dim; d++)
            //    K += GaussWeight * Density * Viscosity * rDN_DX(i, d) * rDN_DX(j, d);

            for (unsigned int d = 0; d < Dim; d++)
            {
                //K += GaussWeight * Density * Viscosity * rDN_DX(i, d) * rDN_DX(j, d);
                LHS(row+d,col+d) += K;

                // v * Grad(p) block
                G = tau_one * AGradN[i] * rData.DN_DX(j,d); // Stabilization: (a * Grad(v)) * tau_one * Grad(p)
                PDivV = rData.DN_DX(i,d) * rData.N[j]; // Div(v) * p

                // Write v * Grad(p) component
                LHS(row+d,col+Dim) += rData.Weight * (G - PDivV);
                // Use symmetry to write the q * Div(u) component
                LHS(col+Dim,row+d) += rData.Weight * (G + PDivV);

                // q-p stabilization block
                laplacian += rData.DN_DX(i,d) * rData.DN_DX(j,d); // Stabilization: Grad(q) * tau_one * Grad(p)

                for (unsigned int e = 0; e < Dim; e++) // Stabilization: Div(v) * tau_two * Div(u)
                    LHS(row+d,col+e) += rData.Weight*tau_two*rData.DN_DX(i,d)*rData.DN_DX(j,e);
            }

            // Write q-p term
            LHS(row+Dim,col+Dim) += rData.Weight*tau_one*laplacian;
        }

        // RHS terms
        double forcing = 0.0;
        for (unsigned int d = 0; d < Dim; ++d)
        {
            rLocalRHS[row+d] += rData.Weight * rData.N[i] * body_force[d]; // v*BodyForce
            rLocalRHS[row+d] += rData.Weight * tau_one * AGradN[i] * ( body_force[d] - momentum_projection[d]); // ( a * Grad(v) ) * tau_one * (Density * BodyForce)
            rLocalRHS[row+d] -= rData.Weight * tau_two * rData.DN_DX(i,d) * mass_projection;
            forcing += rData.DN_DX(i, d) * (body_force[d] - momentum_projection[d]);
        }
        rLocalRHS[row + Dim] += rData.Weight * tau_one * forcing; // Grad(q) * tau_one * (Density * BodyForce)
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

///////////////////////////////////////////////////////////////////////////////////////////////////

template< class TElementData >
void QSVMS<TElementData>::AddMassLHS(
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

    /* Note on OSS and full projection: Riccardo says that adding the terms provided by
     * AddMassStabilization (and incluiding their corresponding terms in the projeciton)
     * could help reduce the non-linearity of the coupling between projection and u,p
     * However, leaving them on gives a lot of trouble whith the Bossak scheme:
     * think that we solve F - (1-alpha)*M*u^(n+1) - alpha*M*u^(n) - K(u^(n+1)) = 0
     * so the projection of the dynamic terms should be Pi( (1-alpha)*u^(n+1) - alpha*u^(n) )
     */
    if ( rData.UseOSS != 1.0 )
        this->AddMassStabilization(rData,rMassMatrix);
}

///////////////////////////////////////////////////////////////////////////////////////////////////

template< class TElementData >
void QSVMS<TElementData>::AddMassStabilization(
    TElementData& rData,
    MatrixType &rMassMatrix)
{
    const double density = this->GetAtCoordinate(rData.Density,rData.N);

    double tau_one;
    double tau_two;
    array_1d<double,3> convective_velocity = this->GetAtCoordinate(rData.Velocity,rData.N) - this->GetAtCoordinate(rData.MeshVelocity,rData.N);
    this->CalculateTau(rData,convective_velocity,tau_one,tau_two);

    Vector AGradN;
    this->ConvectionOperator(AGradN,convective_velocity,rData.DN_DX);

    // Multiplying some quantities by density to have correct units
    AGradN *= density; // Convective term is always multiplied by density

    // Temporary container
    double K;
    double weight = rData.Weight * tau_one * density; // This density is for the dynamic term in the residual (rho*Du/Dt)

    // Note: Dof order is (u,v,[w,]p) for each node
    for (unsigned int i = 0; i < NumNodes; i++)
    {
        unsigned int row = i*BlockSize;

        for (unsigned int j = 0; j < NumNodes; j++)
        {
            unsigned int col = j*BlockSize;

            K = weight * AGradN[i] * rData.N[j];

            for (unsigned int d = 0; d < Dim; d++)
            {
                rMassMatrix(row+d,col+d) += K;
                rMassMatrix(row+Dim,col+d) += weight*rData.DN_DX(i,d)*rData.N[j];
            }
        }
    }
}

template <class TElementData>
void QSVMS<TElementData>::AddBoundaryTraction(TElementData& rData,
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
                rLHS(row,col) -= wni*normal_stress_operator(d,col);
            }
            rRHS[row] += wni*(shear_stress[d]-p_gauss*rUnitNormal[d]);
        }
    }
}

template <class TElementData>
void QSVMS<TElementData>::AddViscousTerm(
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

template <class TElementData>
double QSVMS<TElementData>::EffectiveViscosity(
    TElementData& rData, double ElementSize) {

    double c_s = rData.CSmagorinsky;
    double viscosity = rData.DynamicViscosity; //this->GetAtCoordinate(rData.DynamicViscosity, rData.N);

    if (c_s != 0.0) {
        const double density = this->GetAtCoordinate(rData.Density,rData.N);
        const auto& r_velocities = rData.Velocity;
        const auto& r_dndx = rData.DN_DX;

        // Calculate Symetric gradient
        MatrixType strain_rate = ZeroMatrix(Dim, Dim);
        for (unsigned int n = 0; n < NumNodes; ++n) {
            for (unsigned int i = 0; i < Dim; ++i)
                for (unsigned int j = 0; j < Dim; ++j)
                    strain_rate(i, j) +=
                        0.5 * (r_dndx(n, j) * r_velocities(n, i) +
                                  r_dndx(n, i) * r_velocities(n, j));
        }

        // Norm of symetric gradient
        double strain_rate_norm = 0.0;
        for (unsigned int i = 0; i < Dim; ++i)
            for (unsigned int j = 0; j < Dim; ++j)
                strain_rate_norm += strain_rate(i, j) * strain_rate(i, j);
        strain_rate_norm = sqrt(2.0 * strain_rate_norm);

        // Nu_sgs = (c_s * Delta)^2 * (2*Sij*Sij)^(1/2)
        viscosity +=
            density * c_s * c_s * ElementSize * ElementSize * strain_rate_norm;
    }

    return viscosity;
}

///////////////////////////////////////////////////////////////////////////////////////////////////

template< class TElementData >
void QSVMS<TElementData>::CalculateTau(
    const TElementData& rData,
    const array_1d<double,3> &Velocity,
    double &TauOne,
    double &TauTwo) const
{
    constexpr double c1 = 8.0;
    constexpr double c2 = 2.0;

    const double h = rData.ElementSize;
    const double density = this->GetAtCoordinate(rData.Density,rData.N);
    const double viscosity = this->GetAtCoordinate(rData.EffectiveViscosity,rData.N);

    double velocity_norm = Velocity[0]*Velocity[0];
    for (unsigned int d = 1; d < Dim; d++)
        velocity_norm += Velocity[d]*Velocity[d];
    velocity_norm = std::sqrt(velocity_norm);

    double inv_tau = c1 * viscosity / (h*h) + density * ( rData.DynamicTau/rData.DeltaTime + c2 * velocity_norm / h );
    TauOne = 1.0/inv_tau;
    TauTwo = viscosity + c2 * density * velocity_norm * h / c1;
}


///////////////////////////////////////////////////////////////////////////////////////////////////

template< class TElementData >
void QSVMS<TElementData>::CalculateProjections(const ProcessInfo &rCurrentProcessInfo)
{
    // Get Shape function data
    Vector GaussWeights;
    Matrix ShapeFunctions;
    ShapeFunctionDerivativesArrayType ShapeDerivatives;
    this->CalculateGeometryData(GaussWeights,ShapeFunctions,ShapeDerivatives);
    const unsigned int NumGauss = GaussWeights.size();

    array_1d<double,NumNodes*Dim> momentum_rhs = ZeroVector(NumNodes*Dim);
    VectorType MassRHS = ZeroVector(NumNodes);
    VectorType NodalArea = ZeroVector(NumNodes);

    TElementData data;
    data.Initialize(*this, rCurrentProcessInfo);

    for (unsigned int g = 0; g < NumGauss; g++)
    {
        this->UpdateIntegrationPointData(data, g, GaussWeights[g], row(ShapeFunctions, g), ShapeDerivatives[g]);

        array_1d<double, 3> MomentumRes = ZeroVector(3);
        double MassRes = 0.0;

        array_1d<double,3> convective_velocity = this->GetAtCoordinate(data.Velocity,data.N) - this->GetAtCoordinate(data.MeshVelocity,data.N);

        this->MomentumProjTerm(data, convective_velocity, MomentumRes);
        this->MassProjTerm(data,MassRes);

        for (unsigned int i = 0; i < NumNodes; i++)
        {
            double W = data.Weight*data.N[i];
            unsigned int row = i*Dim;
            for (unsigned int d = 0; d < Dim; d++)
                momentum_rhs[row+d] += W*MomentumRes[d];
            NodalArea[i] += W;
            MassRHS[i] += W*MassRes;
        }
    }

    // Add carefully to nodal variables to avoid OpenMP race condition
    GeometryType& r_geometry = this->GetGeometry();
    for (SizeType i = 0; i < NumNodes; ++i)
    {
        r_geometry[i].SetLock(); // So it is safe to write in the node in OpenMP
        array_1d<double,3>& rMomValue = r_geometry[i].FastGetSolutionStepValue(ADVPROJ);
        unsigned int row = i*Dim;
        for (unsigned int d = 0; d < Dim; d++)
            rMomValue[d] += momentum_rhs[row + d];
        r_geometry[i].FastGetSolutionStepValue(DIVPROJ) += MassRHS[i];
        r_geometry[i].FastGetSolutionStepValue(NODAL_AREA) += NodalArea[i];
        r_geometry[i].UnSetLock(); // Free the node for other threads
    }
}

template< class TElementData >
void QSVMS<TElementData>::SubscaleVelocity(
    const TElementData& rData,
    array_1d<double,3> &rVelocitySubscale) const
{
    double tau_one;
    double tau_two;
    array_1d<double,3> convective_velocity = this->GetAtCoordinate(rData.Velocity,rData.N) - this->GetAtCoordinate(rData.MeshVelocity,rData.N);
    this->CalculateTau(rData,convective_velocity,tau_one,tau_two);

    array_1d<double,3> Residual = ZeroVector(3);

    if (rData.UseOSS != 1.0)
        this->AlgebraicMomentumResidual(rData,convective_velocity,Residual);
    else
        this->OrthogonalMomentumResidual(rData,convective_velocity,Residual);

    rVelocitySubscale = tau_one*Residual;
}

template< class TElementData >
void QSVMS<TElementData>::SubscalePressure(
        const TElementData& rData,
        double &rPressureSubscale) const
{
    double tau_one;
    double tau_two;
    array_1d<double, 3> convective_velocity =
        this->GetAtCoordinate(rData.Velocity, rData.N) -
        this->GetAtCoordinate(rData.MeshVelocity, rData.N);
    this->CalculateTau(rData, convective_velocity, tau_one, tau_two);

    double Residual = 0.0;

    if (rData.UseOSS != 1.0)
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
void QSVMS<TElementData>::save(Serializer& rSerializer) const
{
    typedef FluidElement<TElementData> BaseElement;
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseElement );
}


template< class TElementData >
void QSVMS<TElementData>::load(Serializer& rSerializer)
{
    typedef FluidElement<TElementData> BaseElement;
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseElement);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Class template instantiation

template class QSVMS< QSVMSData<2,3> >;
template class QSVMS< QSVMSData<3,4> >;

template class QSVMS< QSVMSData<2,4> >;
template class QSVMS< QSVMSData<3,8> >;

template class QSVMS< TimeIntegratedQSVMSData<2,3> >;
template class QSVMS< TimeIntegratedQSVMSData<3,4> >;

template class QSVMS< QSVMSDEMCoupledData<2,3> >;
template class QSVMS< QSVMSDEMCoupledData<2,6> >;
template class QSVMS< QSVMSDEMCoupledData<3,4> >;

template class QSVMS< QSVMSDEMCoupledData<2,4> >;
template class QSVMS< QSVMSDEMCoupledData<2,9> >;
template class QSVMS< QSVMSDEMCoupledData<3,8> >;
template class QSVMS< QSVMSDEMCoupledData<3,27> >;

} // namespace Kratos
