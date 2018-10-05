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

#include "d_vms.h"
#include "includes/cfd_variables.h"
#include "includes/checks.h"

#include "custom_utilities/qsvms_data.h"
//#include "custom_utilities/time_integrated_qsvms_data.h"
#include "custom_utilities/fluid_element_utilities.h"

namespace Kratos
{

///////////////////////////////////////////////////////////////////////////////////////////////////
// Life cycle

template< class TElementData >
DVMS<TElementData>::DVMS(IndexType NewId):
    QSVMS<TElementData>(NewId),
    mPredictedSubscaleVelocity(),
    mOldSubscaleVelocity()
{}

template< class TElementData >
DVMS<TElementData>::DVMS(IndexType NewId, const NodesArrayType& ThisNodes):
    QSVMS<TElementData>(NewId,ThisNodes),
    mPredictedSubscaleVelocity(),
    mOldSubscaleVelocity()
{}


template< class TElementData >
DVMS<TElementData>::DVMS(IndexType NewId, GeometryType::Pointer pGeometry):
    QSVMS<TElementData>(NewId,pGeometry),
    mPredictedSubscaleVelocity(),
    mOldSubscaleVelocity()
{}


template< class TElementData >
DVMS<TElementData>::DVMS(IndexType NewId, GeometryType::Pointer pGeometry, Properties::Pointer pProperties):
    QSVMS<TElementData>(NewId,pGeometry,pProperties),
    mPredictedSubscaleVelocity(),
    mOldSubscaleVelocity()
{}


template< class TElementData >
DVMS<TElementData>::~DVMS()
{}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Public Operations

template< class TElementData >
Element::Pointer DVMS<TElementData>::Create(IndexType NewId,NodesArrayType const& ThisNodes,Properties::Pointer pProperties) const
{
    return Kratos::make_shared<DVMS>(NewId, this->GetGeometry().Create(ThisNodes), pProperties);
}


template< class TElementData >
Element::Pointer DVMS<TElementData>::Create(IndexType NewId,GeometryType::Pointer pGeom,Properties::Pointer pProperties) const
{
    return Kratos::make_shared<DVMS>(NewId, pGeom, pProperties);
}

template <class TElementData>
void DVMS<TElementData>::Calculate(const Variable<double>& rVariable,
    double& rOutput, const ProcessInfo& rCurrentProcessInfo) {}

template <class TElementData>
void DVMS<TElementData>::Calculate(
    const Variable<array_1d<double, 3>>& rVariable,
    array_1d<double, 3>& rOutput, const ProcessInfo& rCurrentProcessInfo) {
    // Lumped projection terms
    if (rVariable == ADVPROJ) {
        this->CalculateProjections(rCurrentProcessInfo);
    }
}

template <class TElementData>
void DVMS<TElementData>::Calculate(const Variable<Vector>& rVariable,
    Vector& rOutput, const ProcessInfo& rCurrentProcessInfo) {}

template <class TElementData>
void DVMS<TElementData>::Calculate(const Variable<Matrix>& rVariable,
    Matrix& rOutput, const ProcessInfo& rCurrentProcessInfo) {}

template <class TElementData>
void DVMS<TElementData>::Initialize()
{
    // Base class does things with constitutive law here.
    QSVMS<TElementData>::Initialize();

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

    #ifdef KRATOS_D_VMS_SUBSCALE_ERROR_INSTRUMENTATION
    mSubscaleIterationError.resize(number_of_gauss_points);
    mSubscaleIterationCount.resize(number_of_gauss_points);
    #endif
}

template <class TElementData>
void DVMS<TElementData>::FinalizeSolutionStep(ProcessInfo &rCurrentProcessInfo)
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
        data.UpdateGeometryValues(g, gauss_weights[g],row(shape_functions,g),shape_function_derivatives[g]);

        this->CalculateMaterialResponse(data);

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
void DVMS<TElementData>::InitializeNonLinearIteration(ProcessInfo &rCurrentProcessInfo)
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
        data.UpdateGeometryValues(g, gauss_weights[g],row(shape_functions,g),shape_function_derivatives[g]);
        this->CalculateMaterialResponse(data);

        this->UpdateSubscaleVelocityPrediction(data);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Inquiry

template< class TElementData >
int DVMS<TElementData>::Check(const ProcessInfo &rCurrentProcessInfo)
{
    int out = QSVMS<TElementData>::Check(rCurrentProcessInfo);
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
void DVMS<TElementData>::CalculateOnIntegrationPoints(
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
        const unsigned int number_of_integration_points = GaussWeights.size();

        rValues.resize(number_of_integration_points);

        // fill only if the subscale has already been initialized
        if (mPredictedSubscaleVelocity.size() > 0) {

            TElementData data;
            data.Initialize(*this, rCurrentProcessInfo);

            for (unsigned int g = 0; g < number_of_integration_points; g++) {

                data.UpdateGeometryValues(g, GaussWeights[g], row(ShapeFunctions, g), ShapeDerivatives[g]);
                this->CalculateMaterialResponse(data);

                this->SubscaleVelocity(data, rValues[g]);
            }
        }
        else {
            for (unsigned int g = 0; g < number_of_integration_points; g++) {
                rValues[g].clear();
            }
        }
    }
    else {
        QSVMS<TElementData>::CalculateOnIntegrationPoints(rVariable,rValues,rCurrentProcessInfo);
    }
}


template< class TElementData >
void DVMS<TElementData>::CalculateOnIntegrationPoints(
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
        const unsigned int number_of_integration_points = GaussWeights.size();

        rValues.resize(number_of_integration_points);

        // fill only if the subscale has already been initialized
        if (mPredictedSubscaleVelocity.size() > 0) {

            TElementData data;
            data.Initialize(*this, rCurrentProcessInfo);

            for (unsigned int g = 0; g < number_of_integration_points; g++) {

                data.UpdateGeometryValues(g, GaussWeights[g], row(ShapeFunctions, g), ShapeDerivatives[g]);
                this->CalculateMaterialResponse(data);

                this->SubscalePressure(data, rValues[g]);
            }
        }
        else {
            for (unsigned int g = 0; g < number_of_integration_points; g++) {
                rValues[g] = 0.0;
            }
        }
    }
    #ifdef KRATOS_D_VMS_SUBSCALE_ERROR_INSTRUMENTATION
    else if (rVariable == RESIDUAL_NORM) {
        if (mSubscaleIterationError.size() != 0) {
            rValues.resize(mSubscaleIterationError.size());
            for (unsigned int g = 0; g < mSubscaleIterationError.size(); g++)
                rValues[g] = mSubscaleIterationError[g];
        }
        else {
            rValues.resize(this->GetGeometry().IntegrationPointsNumber(this->GetIntegrationMethod()));
            for (unsigned int g = 0; g < rValues.size(); g++)
                rValues[g] = -1.0;
        }
    }
    else if (rVariable == FLAG_VARIABLE) {
        if (mSubscaleIterationCount.size() != 0) {
            rValues.resize(mSubscaleIterationCount.size());
            for (unsigned int g = 0; g < mSubscaleIterationCount.size(); g++)
                rValues[g] = mSubscaleIterationCount[g];
        }
        else {
            rValues.resize(this->GetGeometry().IntegrationPointsNumber(this->GetIntegrationMethod()));
            for (unsigned int g = 0; g < rValues.size(); g++)
                rValues[g] = -1.0;
        }
    }
    #endif
    else {
        QSVMS<TElementData>::CalculateOnIntegrationPoints(rVariable,rValues,rCurrentProcessInfo);
    }
}

template <class TElementData>
void DVMS<TElementData>::CalculateOnIntegrationPoints(
    Variable<array_1d<double, 6>> const& rVariable,
    std::vector<array_1d<double, 6>>& rValues,
    ProcessInfo const& rCurrentProcessInfo)
{
    QSVMS<TElementData>::CalculateOnIntegrationPoints(rVariable,rValues,rCurrentProcessInfo);
}

template <class TElementData>
void DVMS<TElementData>::CalculateOnIntegrationPoints(
    Variable<Vector> const& rVariable,
    std::vector<Vector>& rValues,
    ProcessInfo const& rCurrentProcessInfo)
{
    QSVMS<TElementData>::CalculateOnIntegrationPoints(rVariable,rValues,rCurrentProcessInfo);
}

template <class TElementData>
void DVMS<TElementData>::CalculateOnIntegrationPoints(
    Variable<Matrix> const& rVariable,
    std::vector<Matrix>& rValues,
    ProcessInfo const& rCurrentProcessInfo)
{
    QSVMS<TElementData>::CalculateOnIntegrationPoints(rVariable,rValues,rCurrentProcessInfo);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Input and output

template< class TElementData >
std::string DVMS<TElementData>::Info() const
{
    std::stringstream buffer;
    buffer << "DVMS #" << this->Id();
    return buffer.str();
}


template< class TElementData >
void DVMS<TElementData>::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "DVMS" << Dim << "D";
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Protected functions

///////////////////////////////////////////////////////////////////////////////////////////////////
// Evaluation of system terms on Gauss Points

template< class TElementData >
void DVMS<TElementData>::AddVelocitySystem(
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
            double L = 0;

            for (unsigned int d = 0; d < Dim; d++) {
                LHS(row+d,col+d) += K;

                /* v * Grad(p) block */
                // Stabilization: (a * Grad(v)) * TauOne * Grad(p)
                double G = tau_one * AGradN[i] * rData.DN_DX(j,d);

                // From vh*d(u_ss)/dt: vh * TauOne/Dt * Grad(p)
                double Sg = tau_one * density*rData.N[i]/dt * rData.DN_DX(j,d);

                // Galerkin pressure term: Div(v) * p
                double PDivV = rData.DN_DX(i,d) * rData.N[j];

                // Write v * Grad(p) component
                LHS(row+d,col+Dim) += rData.Weight * (G - Sg - PDivV);
                // Use symmetry to write the q * Div(u) component
                LHS(col+Dim,row+d) += rData.Weight * (G + PDivV);

                /* q-p stabilization block */
                // Stabilization: Grad(q) * TauOne * Grad(p)
                L += rData.DN_DX(i,d) * rData.DN_DX(j,d);

                /* v-u block */
                // Stabilization: Div(v) * TauTwo * Div(u)
                for (unsigned int e = 0; e < Dim; e++)
                    LHS(row+d,col+e) += rData.Weight*(tau_two + tau_p)*rData.DN_DX(i,d)*rData.DN_DX(j,e);
            }

            // Write q-p term
            LHS(row+Dim,col+Dim) += rData.Weight*tau_one*L;
        }

        // RHS terms
        double qF = 0.0;
        for (unsigned int d = 0; d < Dim; ++d)
        {
            // v*BodyForce + v * du_ss/dt
            rLocalRHS[row+d] += rData.Weight * rData.N[i] * (body_force[d] + OldUssTerm[d]);

            // ( a * Grad(v) ) * TauOne * (Density * BodyForce - Projection)
            // vh * TauOne/Dt * f (from vh*d(uss)/dt
            rLocalRHS[row+d] += rData.Weight * tau_one * (AGradN[i] - density*rData.N[i]/dt ) * ( body_force[d] - MomentumProj[d] + OldUssTerm[d] );

            // OSS pressure subscale projection
            rLocalRHS[row+d] -= rData.Weight * rData.DN_DX(i,d) * (tau_two + tau_p ) * MassProj;

            // Dynamic term in pressure subscale div(vh) * h^2/(c1*dt) * (-div(uh^n) )
            rLocalRHS[row+d] -= rData.Weight * rData.DN_DX(i,d) * tau_p * OldResidual;

            // Grad(q) * TauOne * (Density * BodyForce - Projection)
            qF += rData.DN_DX(i, d) * (body_force[d] - MomentumProj[d] + OldUssTerm[d]);
        }
        rLocalRHS[row+Dim] += rData.Weight * tau_one * qF; // Grad(q) * TauOne * (Density * BodyForce)
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
void DVMS<TElementData>::AddMassLHS(
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
void DVMS<TElementData>::AddMassStabilization(
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
                rMassMatrix(row+Dim,col+d) += W*rData.DN_DX(i,d)*rData.N[j];
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////

template< class TElementData >
void DVMS<TElementData>::CalculateProjections(const ProcessInfo &rCurrentProcessInfo)
{
    // Get Shape function data
    Vector GaussWeights;
    Matrix ShapeFunctions;
    ShapeFunctionDerivativesArrayType ShapeDerivatives;
    this->CalculateGeometryData(GaussWeights,ShapeFunctions,ShapeDerivatives);
    const unsigned int NumGauss = GaussWeights.size();

    VectorType MomentumRHS = ZeroVector(NumNodes * Dim);
    VectorType MassRHS = ZeroVector(NumNodes);
    VectorType NodalArea = ZeroVector(NumNodes);

    TElementData data;
    data.Initialize(*this, rCurrentProcessInfo);

    for (unsigned int g = 0; g < NumGauss; g++)
    {
        data.UpdateGeometryValues(g, GaussWeights[g], row(ShapeFunctions, g), ShapeDerivatives[g]);
        this->CalculateMaterialResponse(data);

        array_1d<double, 3> MomentumRes = ZeroVector(3);
        double MassRes = 0.0;

        array_1d<double, 3> convection_velocity = this->FullConvectiveVelocity(data);

        this->MomentumProjTerm(data, convection_velocity, MomentumRes);
        this->MassProjTerm(data, MassRes);

        for (unsigned int i = 0; i < NumNodes; i++)
        {
            double W = data.Weight*data.N[i];
            unsigned int Row = i*Dim;
            for (unsigned int d = 0; d < Dim; d++)
                MomentumRHS[Row+d] += W*MomentumRes[d];
            NodalArea[i] += W;
            MassRHS[i] += W*MassRes;
        }
    }

    // Add carefully to nodal variables to avoid OpenMP race condition
    GeometryType& r_geometry = this->GetGeometry();
    unsigned int Row = 0;
    for (SizeType i = 0; i < NumNodes; ++i)
    {
        r_geometry[i].SetLock(); // So it is safe to write in the node in OpenMP
        array_1d<double,3>& rMomValue = r_geometry[i].FastGetSolutionStepValue(ADVPROJ);
        for (unsigned int d = 0; d < Dim; ++d)
            rMomValue[d] += MomentumRHS[Row++];
        r_geometry[i].FastGetSolutionStepValue(DIVPROJ) += MassRHS[i];
        r_geometry[i].FastGetSolutionStepValue(NODAL_AREA) += NodalArea[i];
        r_geometry[i].UnSetLock(); // Free the node for other threads
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Private functions
///////////////////////////////////////////////////////////////////////////////////////////////////


// Implementation details

template< class TElementData >
void DVMS<TElementData>::CalculateStabilizationParameters(
    const TElementData& rData,
    const array_1d<double,3> &Velocity,
    double &TauOne,
    double &TauTwo,
    double &TauP) const
{
    const double h = rData.ElementSize;
    const double density = this->GetAtCoordinate(rData.Density,rData.N);
    const double viscosity = this->GetAtCoordinate(rData.EffectiveViscosity,rData.N);

    double velocity_norm = Velocity[0]*Velocity[0];
    for (unsigned int d = 1; d < Dim; d++)
        velocity_norm += Velocity[d]*Velocity[d];
    velocity_norm = std::sqrt(velocity_norm);

    double inv_tau = mTauC1 * viscosity / (h*h) + density * ( 1.0/rData.DeltaTime + mTauC2 * velocity_norm / h );
    TauOne = 1.0/inv_tau;
    TauTwo = viscosity + density * mTauC2 * velocity_norm * h / mTauC1;

    // Auxiliary coefficient StaticTauOne*TauTwo/Dt that appears on the pressure subscale model
    TauP = density * h*h / (mTauC1*rData.DeltaTime);
}

template< class TElementData >
void DVMS<TElementData>::SubscaleVelocity(
    const TElementData& rData,
    array_1d<double,3>& rVelocitySubscale) const
{
    const double density = this->GetAtCoordinate(rData.Density,rData.N);
    array_1d<double,3> convective_velocity = this->FullConvectiveVelocity(rData);

    double tau_one;
    double tau_two;
    double tau_p;
    this->CalculateStabilizationParameters(rData,convective_velocity,tau_one,tau_two,tau_p);

    const double dt = rData.DeltaTime;

    array_1d<double,3> residual = ZeroVector(3);

    if (rData.UseOSS != 1.0) {
        this->AlgebraicMomentumResidual(rData,convective_velocity,residual);
    }
    else {
        this->OrthogonalMomentumResidual(rData,convective_velocity,residual);
    }

    // Note: residual is always of size 3, but stored subscale is of size Dim
    const auto& rOldSubscaleVelocity = mOldSubscaleVelocity[rData.IntegrationPointIndex];
    for (unsigned int d = 0; d < Dim; d++) {
        rVelocitySubscale[d] = tau_one*(residual[d] + (density/dt)*rOldSubscaleVelocity[d]);
    }
}

template< class TElementData >
void DVMS<TElementData>::SubscalePressure(
    const TElementData& rData,
    double &rPressureSubscale) const
{
    array_1d<double,3> convective_velocity = this->FullConvectiveVelocity(rData);

    double tau_one;
    double tau_two;
    double tau_p;
    this->CalculateStabilizationParameters(rData,convective_velocity,tau_one,tau_two,tau_p);

    // Old mass residual for dynamic pressure subscale
    // Note: Residual is defined as -Div(u) [- Projection (if OSS)]
    double old_residual = 0.0;
    const Geometry<Node<3>>& r_geometry = this->GetGeometry();
    for (unsigned int a = 0; a < NumNodes; a++) {
        const array_1d<double,3>& r_old_velocity = r_geometry[a].FastGetSolutionStepValue(VELOCITY,1);
        double old_divergence_projection = r_geometry[a].FastGetSolutionStepValue(DIVPROJ,1);
        for (unsigned int d = 0; d < Dim; d++) {
            old_residual -= rData.DN_DX(a,d)*r_old_velocity[d] + rData.N[a]*old_divergence_projection;
        }
    }

    double residual = 0.0;

    if (rData.UseOSS != 1.0)
        this->AlgebraicMassResidual(rData,residual);
    else
        this->OrthogonalMassResidual(rData,residual);

    rPressureSubscale = (tau_two+tau_p)*residual - tau_p*old_residual;
}

template< class TElementData >
array_1d<double,3> DVMS<TElementData>::FullConvectiveVelocity(
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
void DVMS<TElementData>::UpdateSubscaleVelocityPrediction(
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
    if (rData.UseOSS != 1.0)
        this->AlgebraicMomentumResidual(rData,resolved_convection_velocity,static_residual);
    else
        this->OrthogonalMomentumResidual(rData,resolved_convection_velocity,static_residual);

    // Add the time discretization term to obtain the part of the residual that does not change during iteration
    for (unsigned int d = 0; d < Dim; d++)
        static_residual[d] += density/dt * old_subscale_velocity[d];

    // Newton-Raphson iterations for the subscale
    unsigned int iter = 0;
    bool converged = false;
    double subscale_velocity_error = 2.0 * mSubscalePredictionVelocityTolerance;

    BoundedMatrix<double,Dim,Dim> J = ZeroMatrix(Dim,Dim);
    array_1d<double,Dim> rhs = ZeroVector(Dim);
    array_1d<double,Dim> u = mPredictedSubscaleVelocity[rData.IntegrationPointIndex]; // Use last result as initial guess
    array_1d<double,Dim> du = ZeroVector(Dim);

    while ( (!converged) && (iter++ < mSubscalePredictionMaxIterations) ) {

        // Calculate new Tau
        double convection_velocity_norm = 0.0;
        for (unsigned int d = 0; d < Dim; d++) {
            double v_d = resolved_convection_velocity[d] + u[d];
            convection_velocity_norm += v_d*v_d;
        }
        convection_velocity_norm = sqrt(convection_velocity_norm);
        double inv_tau = mTauC1*viscosity/(h*h) + density * ( 1.0/dt + mTauC2*convection_velocity_norm/h );

        // Newton-Raphson LHS
        noalias(J) = density * resolved_velocity_gradient;
        for (unsigned int d = 0; d < Dim; d++)
            J(d,d) += inv_tau;

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

        if (subscale_velocity_norm > mSubscalePredictionVelocityTolerance)
            subscale_velocity_error /= subscale_velocity_norm;

        if ( (subscale_velocity_error <= mSubscalePredictionVelocityTolerance) ||
            (residual_norm <= mSubscalePredictionResidualTolerance) ) // If RHS is zero, dU is zero too, converged.
        {
            converged = true; // If RHS is zero, dU is zero too, converged.
        }
    }

    #ifdef KRATOS_D_VMS_SUBSCALE_ERROR_INSTRUMENTATION
    mSubscaleIterationError[rData.IntegrationPointIndex] = subscale_velocity_error;
    mSubscaleIterationCount[rData.IntegrationPointIndex] = iter;
    #endif

    // Store new subscale values or discard the calculation
    // If not converged, we will not use the subscale in the convective term.
    noalias(mPredictedSubscaleVelocity[rData.IntegrationPointIndex]) = converged ? u : ZeroVector(Dim);
}

// serializer

template< class TElementData >
void DVMS<TElementData>::save(Serializer& rSerializer) const
{
    typedef QSVMS<TElementData> BaseElement;
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseElement );
    rSerializer.save("mOldSubscaleVelocity",mOldSubscaleVelocity);
}


template< class TElementData >
void DVMS<TElementData>::load(Serializer& rSerializer)
{
    typedef QSVMS<TElementData> BaseElement;
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseElement);
    rSerializer.load("mOldSubscaleVelocity",mOldSubscaleVelocity);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Class template instantiation

template class DVMS< QSVMSData<2,3> >;
template class DVMS< QSVMSData<3,4> >;

} // namespace Kratos
