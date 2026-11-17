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
#include <cmath>
#include "includes/cfd_variables.h"
#include "includes/checks.h"

#include "data_containers/qs_vms/qs_vms_data.h"
#include "data_containers/qs_vms_dem_coupled/qs_vms_dem_coupled_data.h"
#include "data_containers/time_integrated_qs_vms/time_integrated_qs_vms_data.h"
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
    return Kratos::make_intrusive<DVMS>(NewId, this->GetGeometry().Create(ThisNodes), pProperties);
}


template< class TElementData >
Element::Pointer DVMS<TElementData>::Create(IndexType NewId,GeometryType::Pointer pGeom,Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<DVMS>(NewId, pGeom, pProperties);
}

template <class TElementData>
void DVMS<TElementData>::Calculate(
    const Variable<double>& rVariable,
    double& rOutput,
    const ProcessInfo& rCurrentProcessInfo)
{
    BaseType::Calculate(rVariable, rOutput, rCurrentProcessInfo);
}

template <class TElementData>
void DVMS<TElementData>::Calculate(
    const Variable<array_1d<double, 3>>& rVariable,
    array_1d<double, 3>& rOutput,
    const ProcessInfo& rCurrentProcessInfo)
{
    // Lumped projection terms
    if (rVariable == ADVPROJ) {
        this->CalculateProjections(rCurrentProcessInfo);
    } else {
        BaseType::Calculate(rVariable, rOutput, rCurrentProcessInfo);
    }
}

template <class TElementData>
void DVMS<TElementData>::Calculate(
    const Variable<Vector>& rVariable,
    Vector& rOutput,
    const ProcessInfo& rCurrentProcessInfo)
{
    BaseType::Calculate(rVariable, rOutput, rCurrentProcessInfo);
}

template <class TElementData>
void DVMS<TElementData>::Calculate(
    const Variable<Matrix>& rVariable,
    Matrix& rOutput,
    const ProcessInfo& rCurrentProcessInfo)
{
    BaseType::Calculate(rVariable, rOutput, rCurrentProcessInfo);
}

template <class TElementData>
void DVMS<TElementData>::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    // Base class does things with constitutive law here.
    QSVMS<TElementData>::Initialize(rCurrentProcessInfo);

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

    if constexpr (TElementData::ElementManagesTimeIntegration) {
        if (mOldContinuitySource.size() != number_of_gauss_points) {
            mOldContinuitySource = ZeroVector(number_of_gauss_points);
        }
    }

    #ifdef KRATOS_D_VMS_SUBSCALE_ERROR_INSTRUMENTATION
    mSubscaleIterationError.resize(number_of_gauss_points);
    mSubscaleIterationCount.resize(number_of_gauss_points);
    #endif
}

template <class TElementData>
void DVMS<TElementData>::FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
    if constexpr (TElementData::ElementManagesTimeIntegration) {
        KRATOS_ERROR_IF_NOT(rCurrentProcessInfo.Has(STEP)) << "DVMS requires STEP to commit subscale history." << std::endl;
        if (mLastCommittedStep == rCurrentProcessInfo[STEP]) {
            return;
        }
    }
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
        if constexpr (TElementData::ElementManagesTimeIntegration) {
            mOldContinuitySource[g] = this->Has(CONTINUITY_SOURCE) ? this->GetValue(CONTINUITY_SOURCE) : 0.0;
        }
    }
    if constexpr (TElementData::ElementManagesTimeIntegration) {
        mLastCommittedStep = rCurrentProcessInfo[STEP];
    }
}


template <class TElementData>
void DVMS<TElementData>::InitializeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo)
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

///////////////////////////////////////////////////////////////////////////////////////////////////
// Inquiry

template< class TElementData >
int DVMS<TElementData>::Check(const ProcessInfo &rCurrentProcessInfo) const
{
    if constexpr (TElementData::ElementManagesTimeIntegration) {
        KRATOS_ERROR_IF_NOT(rCurrentProcessInfo.Has(DELTA_TIME) && std::isfinite(rCurrentProcessInfo[DELTA_TIME]) && rCurrentProcessInfo[DELTA_TIME] > 0.0)
            << "Time-integrated DVMS requires a positive DELTA_TIME." << std::endl;
        KRATOS_ERROR_IF_NOT(rCurrentProcessInfo.Has(BDF_COEFFICIENTS))
            << "Time-integrated DVMS requires BDF_COEFFICIENTS." << std::endl;
        const auto& r_bdf = rCurrentProcessInfo[BDF_COEFFICIENTS];
        KRATOS_ERROR_IF(r_bdf.size() != 3) << "Time-integrated DVMS requires three BDF coefficients." << std::endl;
        for (const auto coefficient : r_bdf) {
            KRATOS_ERROR_IF_NOT(std::isfinite(coefficient)) << "BDF coefficients must be finite." << std::endl;
        }
        for (const auto& r_node : this->GetGeometry()) {
            KRATOS_ERROR_IF(r_node.GetBufferSize() < 3) << "Time-integrated DVMS requires a buffer of at least three steps." << std::endl;
        }
        KRATOS_ERROR_IF_NOT(rCurrentProcessInfo.Has(STEP)) << "Time-integrated DVMS requires STEP." << std::endl;
    }
    int out = QSVMS<TElementData>::Check(rCurrentProcessInfo);
    KRATOS_ERROR_IF_NOT(out == 0)
        << "Error in base class Check for Element " << this->Info() << std::endl
        << "Error code is " << out << std::endl;

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
        Vector gauss_weights;
        Matrix shape_functions;
        ShapeFunctionDerivativesArrayType shape_function_derivatives;
        this->CalculateGeometryData(gauss_weights,shape_functions,shape_function_derivatives);
        const unsigned int number_of_integration_points = gauss_weights.size();

        rValues.resize(number_of_integration_points);

        // fill only if the subscale has already been initialized
        if (mPredictedSubscaleVelocity.size() > 0) {

            TElementData data;
            data.Initialize(*this, rCurrentProcessInfo);

            for (unsigned int g = 0; g < number_of_integration_points; g++) {
                this->UpdateIntegrationPointData(
                    data, g, gauss_weights[g],
                    row(shape_functions,g),shape_function_derivatives[g]);

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
        Vector gauss_weights;
        Matrix shape_functions;
        ShapeFunctionDerivativesArrayType shape_function_derivatives;
        this->CalculateGeometryData(gauss_weights,shape_functions,shape_function_derivatives);
        const unsigned int number_of_integration_points = gauss_weights.size();

        rValues.resize(number_of_integration_points);

        // fill only if the subscale has already been initialized
        if (mPredictedSubscaleVelocity.size() > 0) {

            TElementData data;
            data.Initialize(*this, rCurrentProcessInfo);

            for (unsigned int g = 0; g < number_of_integration_points; g++) {
                this->UpdateIntegrationPointData(
                    data, g, gauss_weights[g],
                    row(shape_functions,g),shape_function_derivatives[g]);
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


template<class TElementData>
const Parameters DVMS<TElementData>::GetSpecifications() const
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
        "compatible_geometries"      : ["Triangle2D3","Tetrahedra3D4"],
        "element_integrates_in_time" : false,
        "required_polynomial_degree_of_geometry" : 1,
        "documentation"   : "This implements a Navier-Stokes element with dynamic Variational MultiScales (VMS) stabilization."
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

template<class TElementData>
Vector DVMS<TElementData>::CalculateContinuitySourceRHS(const ProcessInfo& rCurrentProcessInfo)
{
    Vector weights;
    Matrix shape_functions;
    ShapeFunctionDerivativesArrayType derivatives;
    this->CalculateGeometryData(weights, shape_functions, derivatives);
    TElementData data;
    data.Initialize(*this, rCurrentProcessInfo);
    Vector result = ZeroVector(LocalSize);
    if constexpr (TElementData::ElementManagesTimeIntegration) {
        KRATOS_ERROR_IF(mOldContinuitySource.size() != weights.size())
            << "DVMS source history has not been initialized." << std::endl;
        const double source = this->Has(CONTINUITY_SOURCE) ? this->GetValue(CONTINUITY_SOURCE) : 0.0;
        for (std::size_t g = 0; g < weights.size(); ++g) {
            this->UpdateIntegrationPointData(data, g, weights[g], row(shape_functions, g), derivatives[g]);
            const auto velocity = this->FullConvectiveVelocity(data);
            double tau_one, tau_two, tau_p;
            this->CalculateStabilizationParameters(data, velocity, tau_one, tau_two, tau_p);
            for (unsigned int i = 0; i < NumNodes; ++i) {
                for (unsigned int d = 0; d < Dim; ++d) {
                    result[i * BlockSize + d] += data.Weight * data.DN_DX(i, d) *
                        ((tau_two + tau_p) * source - tau_p * mOldContinuitySource[g]);
                }
                result[i * BlockSize + Dim] += data.Weight * data.N[i] * source;
            }
        }
    }
    return result;
}

template<class TElementData>
Matrix DVMS<TElementData>::GetAcceptedSubscaleState() const
{
    Vector weights;
    Matrix shape_functions;
    ShapeFunctionDerivativesArrayType derivatives;
    this->CalculateGeometryData(weights, shape_functions, derivatives);
    KRATOS_ERROR_IF(mOldSubscaleVelocity.size() != weights.size() || mOldContinuitySource.size() != weights.size())
        << "DVMS accepted subscale state has not been initialized." << std::endl;
    Matrix state = ZeroMatrix(weights.size(), 6);
    for (std::size_t g = 0; g < weights.size(); ++g) {
        for (unsigned int i = 0; i < NumNodes; ++i) {
            state(g, 0) += shape_functions(g, i) * this->GetGeometry()[i].X();
            state(g, 1) += shape_functions(g, i) * this->GetGeometry()[i].Y();
        }
        state(g, 2) = weights[g];
        state(g, 3) = mOldSubscaleVelocity[g][0];
        state(g, 4) = mOldSubscaleVelocity[g][1];
        state(g, 5) = mOldContinuitySource[g];
    }
    return state;
}

template<class TElementData>
void DVMS<TElementData>::SetAcceptedSubscaleState(const Matrix& rState, int LastCommittedStep)
{
    const auto count = this->GetGeometry().IntegrationPointsNumber(this->GetIntegrationMethod());
    KRATOS_ERROR_IF(rState.size1() != count || rState.size2() != 3)
        << "DVMS subscale state must contain two velocities and one source per integration point." << std::endl;
    mOldSubscaleVelocity.resize(count);
    mPredictedSubscaleVelocity.resize(count);
    mOldContinuitySource.resize(count, false);
    for (std::size_t g = 0; g < count; ++g) {
        mOldSubscaleVelocity[g].clear();
        mPredictedSubscaleVelocity[g].clear();
        for (unsigned int d = 0; d < Dim; ++d) {
            if (d < 2) mOldSubscaleVelocity[g][d] = rState(g, d);
        }
        mOldContinuitySource[g] = rState(g, 2);
    }
    mLastCommittedStep = LastCommittedStep;
}

template<class TElementData>
void DVMS<TElementData>::AddTimeIntegratedLHS(TElementData& rData, MatrixType& rLHS)
{
    Vector rhs = ZeroVector(LocalSize);
    this->AddTimeIntegratedSystem(rData, rLHS, rhs);
}

template<class TElementData>
void DVMS<TElementData>::AddTimeIntegratedRHS(TElementData& rData, VectorType& rRHS)
{
    Matrix lhs = ZeroMatrix(LocalSize, LocalSize);
    this->AddTimeIntegratedSystem(rData, lhs, rRHS);
}

template<class TElementData>
void DVMS<TElementData>::AlgebraicMomentumResidual(
    const TElementData& rData,
    const array_1d<double,3>& rConvectionVelocity,
    array_1d<double,3>& rResidual) const
{
    BaseType::AlgebraicMomentumResidual(rData, rConvectionVelocity, rResidual);
    if constexpr (TElementData::ElementManagesTimeIntegration) {
        const double density = this->GetAtCoordinate(rData.Density, rData.N);
        for (unsigned int i = 0; i < NumNodes; ++i) {
            const auto& r_acceleration = this->GetGeometry()[i].FastGetSolutionStepValue(ACCELERATION);
            for (unsigned int d = 0; d < Dim; ++d) {
                const double bdf_acceleration = rData.bdf0 * rData.Velocity(i,d)
                    + rData.bdf1 * rData.Velocity_OldStep1(i,d)
                    + rData.bdf2 * rData.Velocity_OldStep2(i,d);
                rResidual[d] += density * rData.N[i] * (r_acceleration[d] - bdf_acceleration);
            }
        }
    }
}

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

    const double continuity_source = TElementData::ElementManagesTimeIntegration && this->Has(CONTINUITY_SOURCE)
        ? this->GetValue(CONTINUITY_SOURCE) : 0.0;
    const double old_source = TElementData::ElementManagesTimeIntegration
        ? mOldContinuitySource[rData.IntegrationPointIndex] : 0.0;

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
            rLocalRHS[row+d] += rData.Weight * rData.DN_DX(i,d) *
                ((tau_two + tau_p) * continuity_source - tau_p * old_source);

            // Grad(q) * TauOne * (Density * BodyForce - Projection)
            qF += rData.DN_DX(i, d) * (body_force[d] - MomentumProj[d] + OldUssTerm[d]);
        }
        rLocalRHS[row+Dim] += rData.Weight * (tau_one * qF + rData.N[i] * continuity_source);
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
     * AddMassStabilization (and including their corresponding terms in the projection)
     * could help reduce the non-linearity of the coupling between projection and u,p
     * However, leaving them on gives a lot of trouble with the Bossak scheme:
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
    Vector gauss_weights;
    Matrix shape_functions;
    ShapeFunctionDerivativesArrayType shape_function_derivatives;
    this->CalculateGeometryData(gauss_weights,shape_functions,shape_function_derivatives);
    const unsigned int NumGauss = gauss_weights.size();

    VectorType MomentumRHS = ZeroVector(NumNodes * Dim);
    VectorType MassRHS = ZeroVector(NumNodes);
    VectorType NodalArea = ZeroVector(NumNodes);

    TElementData data;
    data.Initialize(*this, rCurrentProcessInfo);

    for (unsigned int g = 0; g < NumGauss; g++)
    {
        this->UpdateIntegrationPointData(
            data, g, gauss_weights[g],
            row(shape_functions,g),shape_function_derivatives[g]);

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
    const Geometry<Node>& r_geometry = this->GetGeometry();
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

    if constexpr (TElementData::ElementManagesTimeIntegration) {
        residual += this->Has(CONTINUITY_SOURCE) ? this->GetValue(CONTINUITY_SOURCE) : 0.0;
        old_residual += mOldContinuitySource[rData.IntegrationPointIndex];
    }
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
    if constexpr (TElementData::ElementManagesTimeIntegration) {
        rSerializer.save("LastCommittedStep", mLastCommittedStep);
        rSerializer.save("OldContinuitySource", mOldContinuitySource);
    }
}


template< class TElementData >
void DVMS<TElementData>::load(Serializer& rSerializer)
{
    typedef QSVMS<TElementData> BaseElement;
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseElement);
    rSerializer.load("mOldSubscaleVelocity",mOldSubscaleVelocity);
    if constexpr (TElementData::ElementManagesTimeIntegration) {
        rSerializer.load("LastCommittedStep", mLastCommittedStep);
        rSerializer.load("OldContinuitySource", mOldContinuitySource);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Class template instantiation

template class DVMS< QSVMSData<2,3> >;
template class DVMS< QSVMSData<2,4> >;
template class DVMS< TimeIntegratedQSVMSData<2,4> >;
template class DVMS< QSVMSData<3,4> >;

template class DVMS< QSVMSDEMCoupledData<2,3> >;
template class DVMS< QSVMSDEMCoupledData<3,4> >;
template class DVMS< QSVMSDEMCoupledData<3,10> >;

template class DVMS< QSVMSDEMCoupledData<2,4> >;
template class DVMS< QSVMSDEMCoupledData<2,6> >;
template class DVMS< QSVMSDEMCoupledData<2,9> >;
template class DVMS< QSVMSDEMCoupledData<3,8> >;
template class DVMS< QSVMSDEMCoupledData<3,27> >;

} // namespace Kratos
