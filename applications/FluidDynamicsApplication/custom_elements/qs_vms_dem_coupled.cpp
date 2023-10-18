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

template< class TElementData >
void QSVMSDEMCoupled<TElementData>::CalculateOnIntegrationPoints(
    const Variable<double>& rVariable,
    std::vector<double>& rOutput,
    const ProcessInfo& rCurrentProcessInfo)
{
    const GeometryType::IntegrationPointsArrayType integration_points = this->GetGeometry().IntegrationPoints(this->GetIntegrationMethod());
    const SizeType number_of_integration_points = integration_points.size();

    Vector gauss_weights;
    Matrix shape_functions;
    ShapeFunctionDerivativesArrayType shape_derivatives;
    this->CalculateGeometryData(
        gauss_weights, shape_functions, shape_derivatives);

    if (rOutput.size() != number_of_integration_points)
        rOutput.resize(number_of_integration_points);
    TElementData data;
    data.Initialize(*this, rCurrentProcessInfo);
    for (unsigned int g = 0; g < number_of_integration_points; g++ ) {
        data.UpdateGeometryValues(g, gauss_weights[g], row(shape_functions, g),shape_derivatives[g]);
        if (rVariable == PRESSURE) {
            const auto& r_pressure = data.Pressure;
            double value = this->GetAtCoordinate(r_pressure,data.N);
            rOutput[g] = value;
        }
    }
}

template< class TElementData >
void QSVMSDEMCoupled<TElementData>::CalculateOnIntegrationPoints(
    const Variable<Matrix>& rVariable,
    std::vector<Matrix>& rValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    const GeometryType::IntegrationPointsArrayType integration_points = this->GetGeometry().IntegrationPoints(this->GetIntegrationMethod());
    const SizeType number_of_integration_points = integration_points.size();

    Vector gauss_weights;
    Matrix shape_functions;
    ShapeFunctionDerivativesArrayType shape_derivatives;
    this->CalculateGeometryData(
        gauss_weights, shape_functions, shape_derivatives);

    if (rValues.size() != number_of_integration_points)
        rValues.resize(number_of_integration_points);

    TElementData data;
    data.Initialize(*this, rCurrentProcessInfo);

    for (unsigned int g = 0; g < number_of_integration_points; g++ ) {
        data.UpdateGeometryValues(g, gauss_weights[g], row(shape_functions, g),shape_derivatives[g]);
        Matrix value = ZeroMatrix(Dim, Dim);
        if (rVariable == VELOCITY_GRADIENT) {
            const auto& r_velocity = data.Velocity;
            for (unsigned int i = 0; i < NumNodes; i++) {
                for (unsigned int d = 0; d < Dim; d++) {
                    for (unsigned int e = 0; e < Dim; e++)
                        value(d,e) += data.DN_DX(i,d) * r_velocity(i,e);
                }
            }
        }
        rValues[g] = value;
    }
}

template< class TElementData >
void QSVMSDEMCoupled<TElementData>::CalculateOnIntegrationPoints(
    const Variable<array_1d<double,3>>& rVariable,
    std::vector<array_1d<double,3>>& rOutput,
    const ProcessInfo& rCurrentProcessInfo)
{
    const GeometryType::IntegrationPointsArrayType integration_points = this->GetGeometry().IntegrationPoints(this->GetIntegrationMethod());
    const SizeType number_of_integration_points = integration_points.size();

    Vector gauss_weights;
    Matrix shape_functions;
    ShapeFunctionDerivativesArrayType shape_derivatives;
    this->CalculateGeometryData(
        gauss_weights, shape_functions, shape_derivatives);

    if (rOutput.size() != number_of_integration_points)
        rOutput.resize(number_of_integration_points);

    TElementData data;
    data.Initialize(*this, rCurrentProcessInfo);


    for (unsigned int g = 0; g < number_of_integration_points; g++ ) {
        data.UpdateGeometryValues(g, gauss_weights[g], row(shape_functions, g),shape_derivatives[g]);
        array_1d<double,3> value(3,0.0);
        if (rVariable == VELOCITY) {
            const auto& r_velocity = data.Velocity;
            value = this->GetAtCoordinate(r_velocity,data.N);
        }
        if (rVariable == BODY_FORCE) {
            value = this->GetAtCoordinate(data.BodyForce,data.N);
        }
        if (rVariable == PRESSURE_GRADIENT){
            const auto& r_pressure = data.Pressure;
            for (unsigned int i = 0; i < NumNodes; i++) {
                for (unsigned int d = 0; d < Dim; d++) {
                    value[d] += r_pressure[i] * data.DN_DX(i,d);
                }
            }
        }
        rOutput[g] = value;
    }
}

template <class TElementData>
void QSVMSDEMCoupled<TElementData>::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    // Base class does things with constitutive law here.
    QSVMS<TElementData>::Initialize(rCurrentProcessInfo);

    if(Dim == 2){
        if (NumNodes == 9 || NumNodes == 6 || NumNodes == 4)
            mInterpolationOrder = 2;
    }
    else if(Dim == 3){
        if (NumNodes == 10 || NumNodes == 27)
            mInterpolationOrder = 2;
    }

    const unsigned int number_of_gauss_points = this->GetGeometry().IntegrationPointsNumber(this->GetIntegrationMethod());

    // The old velocity may be already defined (if restarting)
    // and we want to keep the loaded values in that case.
    if (mPreviousVelocity.size() != number_of_gauss_points)
    {
        mPreviousVelocity.resize(number_of_gauss_points);
        for (unsigned int g = 0; g < number_of_gauss_points; g++)
            mPreviousVelocity[g] = ZeroVector(Dim);
    }

    if (mPredictedSubscaleVelocity.size() != number_of_gauss_points)
    {
        mPredictedSubscaleVelocity.resize(number_of_gauss_points);
        for (unsigned int g = 0; g < number_of_gauss_points; g++)
            mPredictedSubscaleVelocity[g] = ZeroVector(Dim);
    }
    // The prediction is updated before each non-linear iteration:
    // It is not stored in a restart and can be safely initialized.

    // The old velocity may be already defined (if restarting)
    // and we want to keep the loaded values in that case.
    if (mViscousResistanceTensor.size() != number_of_gauss_points)
    {
        mViscousResistanceTensor.resize(number_of_gauss_points);
        for (unsigned int g = 0; g < number_of_gauss_points; g++)
            mViscousResistanceTensor[g] = ZeroMatrix(Dim,Dim);
    }
}



template< class TElementData >
GeometryData::IntegrationMethod QSVMSDEMCoupled<TElementData>::GetIntegrationMethod() const
{
    if(mInterpolationOrder == 1)
        return GeometryData::IntegrationMethod::GI_GAUSS_2;
    else
        return GeometryData::IntegrationMethod::GI_GAUSS_3;
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

template <class TElementData>
void QSVMSDEMCoupled<TElementData>::InitializeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo)
{
    // Get Shape function data
    Vector gauss_weights;
    Matrix shape_functions;
    ShapeFunctionDerivativesArrayType shape_function_derivatives;
    DenseVector<DenseVector<Matrix>> shape_function_second_derivatives;
    this->CalculateGeometryData(gauss_weights,shape_functions,shape_function_derivatives);
    const unsigned int number_of_integration_points = gauss_weights.size();
    GeometryUtils::ShapeFunctionsSecondDerivativesTransformOnAllIntegrationPoints(
            shape_function_second_derivatives,this->GetGeometry(),this->GetIntegrationMethod());
    TElementData data;
    data.Initialize(*this,rCurrentProcessInfo);

    for (unsigned int g = 0; g < number_of_integration_points; g++) {
        this->UpdateIntegrationPointDataSecondDerivatives(data, g, gauss_weights[g],row(shape_functions,g),shape_function_derivatives[g],shape_function_second_derivatives[g]);

        this->CalculateResistanceTensor(data);
        //this->UpdateSubscaleVelocity(data);
    }
}

template <class TElementData>
void QSVMSDEMCoupled<TElementData>::FinalizeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo)
{
    //Get Shape function data
    Vector gauss_weights;
    Matrix shape_functions;
    ShapeFunctionDerivativesArrayType shape_function_derivatives;
    DenseVector<DenseVector<Matrix>> shape_function_second_derivatives;
    this->CalculateGeometryData(gauss_weights,shape_functions,shape_function_derivatives);
    const unsigned int number_of_integration_points = gauss_weights.size();
    GeometryUtils::ShapeFunctionsSecondDerivativesTransformOnAllIntegrationPoints(
            shape_function_second_derivatives,this->GetGeometry(),this->GetIntegrationMethod());
    TElementData data;
    data.Initialize(*this,rCurrentProcessInfo);

    for (unsigned int g = 0; g < number_of_integration_points; g++) {
        this->UpdateIntegrationPointDataSecondDerivatives(data, g, gauss_weights[g],row(shape_functions,g),shape_function_derivatives[g],shape_function_second_derivatives[g]);

        this->UpdateSubscaleVelocity(data);
    }
}

template <class TElementData>
void QSVMSDEMCoupled<TElementData>::UpdateIntegrationPointDataSecondDerivatives(
    TElementData& rData,
    unsigned int IntegrationPointIndex,
    double Weight,
    const typename TElementData::MatrixRowType& rN,
    const typename TElementData::ShapeDerivativesType& rDN_DX,
    const typename TElementData::ShapeFunctionsSecondDerivativesType& rDDN_DDX) const
{
    this->UpdateIntegrationPointData(rData, IntegrationPointIndex, Weight, rN, rDN_DX);
    rData.UpdateSecondDerivativesValues(rDDN_DDX);
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
    const auto& r_body_forces = rData.BodyForce;
    const auto& r_velocities = rData.Velocity;
    const auto& r_pressures = rData.Pressure;

    MatrixType sigma = mViscousResistanceTensor[rData.IntegrationPointIndex];
    Vector sigma_U, grad_div_u, div_sym_grad_u;

    for (unsigned int i = 0; i < NumNodes; i++) {
        const array_1d<double,Dim>& r_acceleration = rGeom[i].FastGetSolutionStepValue(ACCELERATION);
        array_1d<double,Dim> sigma_U = ZeroVector(Dim);
        grad_div_u = ZeroVector(Dim);
        div_sym_grad_u = ZeroVector(Dim);
        for (unsigned int d = 0; d < Dim; d++) {
            for (unsigned int e = 0; e < Dim; e++){
                sigma_U[d] += sigma(d,e) * rData.N[i] * r_velocities(i,e);
                grad_div_u[d] += rData.DDN_DDX[i](d,e) *  r_velocities(i,e);
                if (d == e)
                    div_sym_grad_u[d] += rData.DDN_DDX[i](e,e) * r_velocities(i,d);
                else
                    div_sym_grad_u[d] += 1.0/2.0 * (rData.DDN_DDX[i](e,d) * r_velocities(i,e) + rData.DDN_DDX[i](e,e) * r_velocities(i,d));
            }
            rResidual[d] += density * ( rData.N[i]*(r_body_forces(i,d) - r_acceleration[d]) - convection[i]*r_velocities(i,d)) + 2.0 * viscosity * div_sym_grad_u[d] - 2.0/3.0 * viscosity * grad_div_u[d] - rData.DN_DX(i,d) * r_pressures[i] - sigma_U[d];
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

    MatrixType sigma = mViscousResistanceTensor[rData.IntegrationPointIndex];
    const auto body_force = this->GetAtCoordinate(rData.BodyForce, rData.N);

    const auto r_velocities = rData.Velocity;
    const auto r_pressures = rData.Pressure;
    const auto velocity = this->GetAtCoordinate(rData.Velocity,rData.N);
    Vector grad_div_u, div_sym_grad_u;

    for (unsigned int i = 0; i < NumNodes; i++) {
        array_1d<double,Dim> sigma_U = ZeroVector(Dim);
        grad_div_u = ZeroVector(Dim);
        div_sym_grad_u = ZeroVector(Dim);
        for (unsigned int d = 0; d < Dim; d++) {
            for (unsigned int e = 0; e < Dim; e++){
                sigma_U[d] += sigma(d,e) * rData.N[i] * r_velocities(i,e);
                grad_div_u[d] += rData.DDN_DDX[i](d,e) *  r_velocities(i,e);
                if (d == e)
                    div_sym_grad_u[d] += rData.DDN_DDX[i](e,e) * r_velocities(i,d);
                else
                    div_sym_grad_u[d] += 1.0/2.0 * (rData.DDN_DDX[i](e,d) * r_velocities(i,e) + rData.DDN_DDX[i](e,e) * r_velocities(i,d));
            }
            rMomentumRHS[d] += density * (- AGradN[i]*r_velocities(i,d)) + 2.0 * viscosity * div_sym_grad_u[d] - 2.0/3.0 * viscosity * grad_div_u[d] - rData.DN_DX(i,d) * r_pressures[i] - sigma_U[d];
        }
    }
    for (unsigned int d = 0; d < Dim; d++)
        rMomentumRHS[d] += density * body_force[d];
}

template<class TElementData>
void QSVMSDEMCoupled<TElementData>::AddReactionStabilization(
    TElementData& rData,
    BoundedMatrix<double,NumNodes*(Dim+1),NumNodes*(Dim+1)>& rLHS,
    VectorType& rLocalRHS)
{

    const double density = this->GetAtCoordinate(rData.Density, rData.N);

    BoundedMatrix<double,Dim,Dim> tau_one = ZeroMatrix(Dim, Dim);
    double tau_two;
    const array_1d<double, 3> convective_velocity =
        this->GetAtCoordinate(rData.Velocity, rData.N) -
        this->GetAtCoordinate(rData.MeshVelocity, rData.N);

    this->CalculateTau(rData, convective_velocity, tau_one, tau_two);

    array_1d<double,3> body_force = density * this->GetAtCoordinate(rData.BodyForce, rData.N);

    Vector AGradN;
    this->ConvectionOperator(AGradN, convective_velocity, rData.DN_DX); // Get a * grad(Ni)

    AGradN *= density;

    const double fluid_fraction = this->GetAtCoordinate(rData.FluidFraction, rData.N);
    double viscosity = this->GetAtCoordinate(rData.DynamicViscosity, rData.N);

    MatrixType sigma = mViscousResistanceTensor[rData.IntegrationPointIndex];

    // Note: Dof order is (vx,vy,[vz,]p) for each node
    for (unsigned int i = 0; i < NumNodes; i++)
    {
        unsigned int row = i*BlockSize;
        // Loop over columns
        for (unsigned int j = 0; j < NumNodes; j++)
        {
            unsigned int col = j*BlockSize;

            for (unsigned int d = 0; d < Dim; d++) // iterate over dimensions for velocity Dofs in this node combination
            {
                double RSigmaG = 0.0;
                double GAlphaR = 0.0;
                for (unsigned int e = 0; e < Dim; e++){
                    double ASigma = tau_one(d,d) * AGradN[i] * sigma(d,e) * rData.N[j];
                    double RRSigma = 0.0;
                    double RSigmaA = tau_one(d,d) * sigma(d,e) * rData.N[i] * AGradN[j];
                    double LSigma_1 = 0.0;
                    double LSigma_2 = 0.0;
                    double CSigma = 0.0;
                    double RSigmaL_1 = 0.0;
                    double RSigmaL_2 = 0.0;
                    double RSigmaC = 0.0;
                    GAlphaR += tau_one(d,d) * fluid_fraction * rData.DN_DX(i,e) * sigma(e,d) * rData.N[j];
                    RSigmaG += tau_one(d,d) * sigma(d,e) * rData.N[i] * rData.DN_DX(j,e);
                    for (unsigned int f = 0; f < Dim; f++){
                        LSigma_1 += tau_one(d,d) * viscosity * rData.DDN_DDX[i](f,f) * sigma(d,e) * rData.N[j];
                        LSigma_2 += tau_one(d,d) * viscosity * rData.DDN_DDX[i](d,f) * sigma(f,e) * rData.N[j];
                        CSigma += 2.0/3.0 * tau_one(d,d) * viscosity * rData.DDN_DDX[i](f,d) * sigma(f,e) * rData.N[j];
                        RRSigma += tau_one(d,d) * sigma(d,f) * rData.N[i] * sigma(f,e) * rData.N[j];
                        RSigmaL_1 += tau_one(d,d) * viscosity * rData.N[i] * sigma(d,e) * rData.DDN_DDX[j](f,f);
                        RSigmaL_2 += tau_one(d,d) * viscosity * rData.N[i] * sigma(d,f) * rData.DDN_DDX[j](e,f);
                        RSigmaC += 2.0/3.0 * tau_one(d,d) * viscosity * sigma(d,f) * rData.N[i] * rData.DDN_DDX[j](f,e);
                    }
                    double LSigma = LSigma_1 + LSigma_2;
                    double RSigmaL = RSigmaL_1 + RSigmaL_2;
                    rLHS(row+d,col+e) += rData.Weight * (ASigma - RRSigma - RSigmaA - CSigma + LSigma - RSigmaL + RSigmaC);
                }
                rLHS(row+Dim,col+d) += rData.Weight * (GAlphaR);
                rLHS(row+d,col+Dim) += rData.Weight * (-RSigmaG);
            }
        }

        // RHS terms
        for (unsigned int d = 0; d < Dim; ++d)
        {
            double RSigmaF = 0.0;
            for (unsigned int e = 0; e < Dim; ++e){
                RSigmaF += tau_one(d,d) * rData.N[i] * sigma(d,e) * body_force[e]; /*- momentum_projection[e]*/ //momentum_projection 0 because is ASGS
            }
            rLocalRHS[row+d] += rData.Weight * (- RSigmaF);
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

    MatrixType sigma = mViscousResistanceTensor[rData.IntegrationPointIndex];

    // Note: Dof order is (vx,vy,[vz,]p) for each node
    for (unsigned int i = 0; i < NumNodes; i++)
    {
        unsigned int row = i*BlockSize;
        // Loop over columns
        for (unsigned int j = 0; j < NumNodes; j++)
        {
            unsigned int col = j*BlockSize;

            for (unsigned int d = 0; d < Dim; d++) // iterate over dimensions for velocity Dofs in this node combination
            {
                double UGAlpha = tau_one(d,d) * fluid_fraction * rData.DN_DX(i,d) * rData.N[j];
                double AU = tau_one(d,d) * AGradN[i] * rData.N[j];
                for (unsigned int e = 0; e < Dim; e++){
                    double LI = tau_one(d,d) * viscosity * rData.N[j] * rData.DDN_DDX[i](d,e);
                    double CI = 2.0/3.0 * tau_one(d,d) * viscosity * rData.DDN_DDX[i](d,e) * rData.N[j];
                    double RSigmaU = tau_one(d,d) * sigma(d,e) * rData.N[i] * rData.N[j];
                    for (unsigned int f = 0; f < Dim; f++)
                        if (d == e)
                            LI += tau_one(d,d) * viscosity * rData.DDN_DDX[i](f,f) * rData.N[j];
                    rMassMatrix(row+d, col+e) += weight * (LI - CI - RSigmaU);
                }
                rMassMatrix(row+d, col+d) += weight * AU;
                rMassMatrix(row+Dim,col+d) += weight * UGAlpha;
            }
        }
    }
}

// Add a the contribution from a single integration point to the velocity contribution
template< class TElementData >
void QSVMSDEMCoupled<TElementData>::AddVelocitySystem(
    TElementData& rData,
    MatrixType& rLocalLHS,
    VectorType& rLocalRHS)
{
    BoundedMatrix<double,NumNodes*(Dim+1),NumNodes*(Dim+1)>& LHS = rData.LHS;
    LHS.clear();

    // Interpolate nodal data on the integration point
    const double density = this->GetAtCoordinate(rData.Density, rData.N);
    array_1d<double,3> body_force = density * this->GetAtCoordinate(rData.BodyForce,rData.N); // Force per unit of volume

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
    array_1d<double,3> fluid_fraction_gradient = this->GetAtCoordinate(rData.FluidFractionGradient, rData.N);

    MatrixType sigma = mViscousResistanceTensor[rData.IntegrationPointIndex];

    // Temporary containers
    //double U;

    // Note: Dof order is (u,v,[w,]p) for each node
    for (unsigned int i = 0; i < NumNodes; i++)
    {

        unsigned int row = i*BlockSize;

        // LHS terms
        for (unsigned int j = 0; j < NumNodes; j++)
        {
            unsigned int col = j*BlockSize;

            // Some terms are the same for all velocity components, calculate them once for each i,j
            double V = rData.N[i] * AGradN[j];

            // q-p stabilization block (initialize result)
            double G = 0;
            for (unsigned int d = 0; d < Dim; d++)
            {

                double AA = tau_one(d,d) * AGradN[i] * AGradN[j]; // Stabilization: u*grad(v) * tau_one * u*grad(u);
                // Stabilization: (a * Grad(v)) * tau_one * Grad(p)
                double P = rData.DN_DX(i,d) * rData.N[j]; // Div(v) * p
                double U = fluid_fraction_gradient[d] * rData.N[j] * rData.N[i];
                double QD = fluid_fraction * rData.N[i] * rData.DN_DX(j,d);

                double GAlphaA = tau_one(d,d) * fluid_fraction * rData.DN_DX(i,d) * AGradN[j];
                double GAlphaL_1 = 0.0;
                double GAlphaL_2 = 0.0;
                double GAlphaC = 0.0;
                double AG = tau_one(d,d) * AGradN[i] * rData.DN_DX(j,d);
                double LG_1 = 0.0;
                double LG_2 = 0.0;
                double CG = 0.0;
                G += tau_one(d,d) * fluid_fraction * rData.DN_DX(i,d) * rData.DN_DX(j,d);
                for (unsigned int e = 0; e < Dim; e++){ // Stabilization: Div(v) * tau_two * Div(u)
                    double DnuD = 2.0/3.0 * viscosity * rData.DN_DX(i,d) * rData.DN_DX(j,e);
                    double GS = viscosity * rData.DN_DX(i,e) * rData.DN_DX(j,d);
                    double DU = tau_two * rData.DN_DX(i,d) * fluid_fraction_gradient[e] * rData.N[j];
                    double DD = tau_two * fluid_fraction * rData.DN_DX(i,d) * rData.DN_DX(j,e);
                    double RSigma = rData.N[i] * sigma(d,e) * rData.N[j];
                    double AL = tau_one(d,d) * viscosity * AGradN[i] * rData.DDN_DDX[j](d,e);
                    double AC = 2.0 / 3.0 * tau_one(d,d) * viscosity * AGradN[i] * rData.DDN_DDX[j](d,e);
                    double LA = tau_one(d,d) * viscosity * rData.DDN_DDX[i](d,e) * AGradN[j];
                    double LL_diag_1 = 0.0;
                    double LL_diag_2 = 0.0;
                    double LL_2 = 0.0;
                    double LL_3 = 0.0;
                    double LL_4 = 0.0;
                    double LC_1 = 0.0;
                    double LC_2 = 0.0;
                    double CA = 2.0/3.0 * tau_one(d,d) * viscosity * rData.DDN_DDX[i](d,e) * AGradN[j];
                    double CL_1 = 0.0;
                    double CL_2 = 0.0;
                    double CC = 0.0;
                    GAlphaL_1 += tau_one(d,d) * viscosity * fluid_fraction * rData.DN_DX(i,d) * rData.DDN_DDX[j](e,e);
                    GAlphaL_2 += tau_one(d,d) * viscosity * fluid_fraction * rData.DN_DX(i,e) * rData.DDN_DDX[j](d,e);
                    GAlphaC += 2.0/3.0 * tau_one(d,d) * viscosity * fluid_fraction * rData.DN_DX(i,e) * rData.DDN_DDX[j](e,d);
                    LG_1 += tau_one(d,d) * viscosity * rData.DDN_DDX[i](e,e) * rData.DN_DX(j,d);
                    LG_2 += tau_one(d,d) * viscosity * rData.DDN_DDX[i](d,e) * rData.DN_DX(j,e);
                    CG += 2.0/3.0 * tau_one(d,d) * viscosity * rData.DDN_DDX[i](d,e) * rData.DN_DX(j,e);
                    for (unsigned int f = 0; f < Dim; f++){
                        if (d == e){
                            GS += viscosity * rData.DN_DX(i,f) * rData.DN_DX(j,f);
                            AL += tau_one(d,d) * viscosity * AGradN[i] * rData.DDN_DDX[j](f,f);
                            LA += tau_one(d,d) * viscosity * rData.DDN_DDX[i](f,f) * AGradN[j];
                            LL_diag_1 += tau_one(d,d) * viscosity * viscosity * rData.DDN_DDX[i](f,f);
                            LL_diag_2 += rData.DDN_DDX[j](f,f);
                        }
                        LL_2 += tau_one(d,d) * viscosity * viscosity * rData.DDN_DDX[i](f,f) * rData.DDN_DDX[j](d,e);
                        LL_3 += tau_one(d,d) * viscosity * viscosity * rData.DDN_DDX[i](d,e) * rData.DDN_DDX[j](f,f);
                        LL_4 += tau_one(d,d) * viscosity * viscosity * rData.DDN_DDX[i](d,f) * rData.DDN_DDX[j](e,f);
                        LC_1 += 2.0/3.0 * tau_one(d,d) * viscosity * viscosity * rData.DDN_DDX[i](f,f) * rData.DDN_DDX[j](d,e);
                        LC_2 += 2.0/3.0 * tau_one(d,d) * viscosity * viscosity * rData.DDN_DDX[i](d,f) * rData.DDN_DDX[j](f,e);
                        CL_2 += 2.0/3.0 * tau_one(d,d) * viscosity * viscosity * rData.DDN_DDX[i](d,f) * rData.DDN_DDX[j](f,e);
                        CL_1 += 2.0/3.0 * tau_one(d,d) * viscosity * viscosity * rData.DDN_DDX[i](d,e) * rData.DDN_DDX[j](f,f);
                        CC += 4.0/9.0 * tau_one(d,d) * viscosity * viscosity * rData.DDN_DDX[i](f,d) * rData.DDN_DDX[j](f,e);
                    }
                    double LL = (LL_diag_1 * LL_diag_2) + LL_2 + LL_3 + LL_4;
                    double LC = LC_1 + LC_2;
                    double CL = CL_1 + CL_2;
                    LHS(row+d,col+e) += rData.Weight * (GS - DnuD - AL + AC + LA - LL + LC - CA + CL - CC + DD + DU + RSigma);
                }
                double GAlphaL = GAlphaL_1 + GAlphaL_2;
                double LG = LG_1 + LG_2;

                LHS(row+d,col+d) += rData.Weight * (V + AA);
                LHS(row+Dim,col+d) += rData.Weight * (GAlphaA + U + QD - GAlphaL + GAlphaC);
                LHS(row+d,col+Dim) += rData.Weight * (AG - P + LG - CG);

            }
            // Write q-p term
            LHS(row+Dim,col+Dim) += rData.Weight * G;

        }

        // RHS terms
        double QAlphaF = 0.0;
        double Q = rData.N[i] * (mass_source - fluid_fraction_rate);
        for (unsigned int d = 0; d < Dim; ++d)
        {
            double VF = rData.N[i] * body_force[d];
            double AF = tau_one(d,d) * AGradN[i] * (body_force[d] - momentum_projection[d]);
            double DPhi = tau_two * rData.DN_DX(i,d) * (mass_source - fluid_fraction_rate - mass_projection);
            double LF_1 = 0.0;
            double LF_2 = 0.0;
            double CF = 0.0;
            QAlphaF += tau_one(d,d) * fluid_fraction * rData.DN_DX(i,d) * (body_force[d] - momentum_projection[d]);
            for (unsigned int e = 0; e < Dim; ++e){
                LF_1 += tau_one(d,d) * viscosity * rData.DDN_DDX[i](e,e) * (body_force[d] - momentum_projection[d]);
                LF_2 += tau_one(d,d) * viscosity * rData.DDN_DDX[i](d,e) * (body_force[e] - momentum_projection[e]);
                CF += 2.0/3.0 * tau_one(d,d) * viscosity * rData.DDN_DDX[i](d,e) * (body_force[e] - momentum_projection[e]);
            }
            double LF = LF_1 + LF_2;
            rLocalRHS[row+d] += rData.Weight * (VF + AF + LF - CF + DPhi);
        }
        rLocalRHS[row+Dim] += rData.Weight * (QAlphaF + Q);
    }
    // Adding reactive terms to the stabilization
    if(!rData.UseOSS)
        this->AddReactionStabilization(rData,LHS,rLocalRHS);

    // Write (the linearized part of the) local contribution into residual form (A*dx = b - A*x)
    array_1d<double,LocalSize> values;
    this->GetCurrentValuesVector(rData, values);
    noalias(rLocalRHS) -= prod(LHS, values);
    /* Viscous contribution (with symmetric gradient 2*nu*{E(u) - 1/3 Tr(E)} )
     * For a generic (potentially non-linear) constitutive law, one cannot assume that RHS = F - LHS*current_values.
     * Because of this, the AddViscousTerm function manages both the LHS and the RHS.
     */
    //this->AddViscousTerm(rData, LHS, rLocalRHS);

    noalias(rLocalLHS) += LHS;

}

template <class TElementData>
void QSVMSDEMCoupled<TElementData>::CalculateMassMatrix(MatrixType& rMassMatrix,
                                                                   const ProcessInfo& rCurrentProcessInfo)
{
    // Resize and intialize output
    if (rMassMatrix.size1() != LocalSize)
        rMassMatrix.resize(LocalSize, LocalSize, false);

    noalias(rMassMatrix) = ZeroMatrix(LocalSize, LocalSize);

    if (!TElementData::ElementManagesTimeIntegration) {
        // Get Shape function data
        Vector gauss_weights;
        Matrix shape_functions;
        ShapeFunctionDerivativesArrayType shape_derivatives;
        DenseVector<DenseVector<Matrix>> shape_function_second_derivatives;
        this->CalculateGeometryData(
            gauss_weights, shape_functions, shape_derivatives);
        GeometryUtils::ShapeFunctionsSecondDerivativesTransformOnAllIntegrationPoints(
            shape_function_second_derivatives,this->GetGeometry(),this->GetIntegrationMethod());
        const unsigned int number_of_gauss_points = gauss_weights.size();

        TElementData data;
        data.Initialize(*this, rCurrentProcessInfo);

        // Iterate over integration points to evaluate local contribution
        for (unsigned int g = 0; g < number_of_gauss_points; g++) {
            this->UpdateIntegrationPointDataSecondDerivatives(
                data, g, gauss_weights[g],
                row(shape_functions, g),shape_derivatives[g],shape_function_second_derivatives[g]);

            this->AddMassLHS(data, rMassMatrix);
        }
    }
}

template <class TElementData>
void QSVMSDEMCoupled<TElementData>::CalculateLocalVelocityContribution(MatrixType& rDampMatrix,
                                                                                  VectorType& rRightHandSideVector,
                                                                                  const ProcessInfo& rCurrentProcessInfo)
{
    // Resize and intialize output
    if( rDampMatrix.size1() != LocalSize )
        rDampMatrix.resize(LocalSize,LocalSize,false);

    if( rRightHandSideVector.size() != LocalSize )
        rRightHandSideVector.resize(LocalSize,false);

    noalias(rDampMatrix) = ZeroMatrix(LocalSize,LocalSize);
    noalias(rRightHandSideVector) = ZeroVector(LocalSize);

    if (!TElementData::ElementManagesTimeIntegration) {
        // Get Shape function data
        Vector gauss_weights;
        Matrix shape_functions;
        ShapeFunctionDerivativesArrayType shape_derivatives;
        DenseVector<DenseVector<Matrix>> shape_function_second_derivatives;
        this->CalculateGeometryData(
            gauss_weights, shape_functions, shape_derivatives);
        const unsigned int number_of_gauss_points = gauss_weights.size();

        GeometryUtils::ShapeFunctionsSecondDerivativesTransformOnAllIntegrationPoints(
            shape_function_second_derivatives,this->GetGeometry(),this->GetIntegrationMethod());

        TElementData data;
        data.Initialize(*this, rCurrentProcessInfo);

        // Iterate over integration points to evaluate local contribution
        for (unsigned int g = 0; g < number_of_gauss_points; g++) {
            const auto& r_dndx = shape_derivatives[g];
            this->UpdateIntegrationPointDataSecondDerivatives(
                data, g, gauss_weights[g],
                row(shape_functions, g),r_dndx, shape_function_second_derivatives[g]);
            this->AddVelocitySystem(data, rDampMatrix, rRightHandSideVector);
        }
    }
}

template < class TElementData >
void QSVMSDEMCoupled<TElementData>::CalculateResistanceTensor(
    const TElementData& rData)
{
    BoundedMatrix<double,Dim,Dim>& rsigma = mViscousResistanceTensor[rData.IntegrationPointIndex];

    BoundedMatrix<double,Dim,Dim> permeability = this->GetAtCoordinate(rData.Permeability, rData.N);

    rsigma = permeability;
}

template< class TElementData >
void QSVMSDEMCoupled<TElementData>::AddMassLHS(
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
    if (!rData.UseOSS)
        this->AddMassStabilization(rData,rMassMatrix);
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
                rMassRHS -= (fluid_fraction * rData.DN_DX(i, d) * velocities(i, d) + fluid_fraction_gradient[d] * rData.N[i] * velocities(i, d));
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
    const int p = mInterpolationOrder;
    const double h = rData.ElementSize;
    const double density = this->GetAtCoordinate(rData.Density,rData.N);
    const double viscosity = this->GetAtCoordinate(rData.EffectiveViscosity,rData.N);
    const double fluid_fraction = this->GetAtCoordinate(rData.FluidFraction, rData.N);
    MatrixType sigma = ZeroMatrix(Dim+1, Dim+1);
    BoundedMatrix<double,Dim,Dim> I = IdentityMatrix(Dim, Dim);
    array_1d<double,3> fluid_fraction_gradient = this->GetAtCoordinate(rData.FluidFractionGradient, rData.N);

    double velocity_modulus = 0.0;
    double fluid_fraction_gradient_modulus = 0.0;

    for (unsigned int d = 0; d < Dim; d++){
        velocity_modulus += Velocity[d] * Velocity[d];
        fluid_fraction_gradient_modulus += std::pow(fluid_fraction_gradient[d],2);
        sigma(d,d) = mViscousResistanceTensor[rData.IntegrationPointIndex](d,d);
    }

    double velocity_norm = std::sqrt(velocity_modulus);
    double fluid_fraction_gradient_norm = std::sqrt(fluid_fraction_gradient_modulus);

    double c_alpha = 1.0 + h / c1 * fluid_fraction_gradient_norm;

    double inv_tau_NS = c1 * viscosity / std::pow(h/(p*p),2.0) + density * (c2 * velocity_norm / (h/p) );
    double tau_one_NS = 1.0 / inv_tau_NS;

    double inv_tau = inv_tau_NS * c_alpha;
    //inv_tau = c1 * viscosity / (h*h) + density * (c2 * velocity_norm / h );

    double tau_one = 1 / (inv_tau + sigma(0,0));
    TauOne = tau_one * I;
    TauTwo = std::pow(h/p,2.0) / (c1 * fluid_fraction * tau_one_NS);
}

template< class TElementData >
void QSVMSDEMCoupled<TElementData>::CalculateProjections(const ProcessInfo &rCurrentProcessInfo)
{
    // Get Shape function data

    DenseVector<DenseVector<Matrix>> ShapeFunctionSecondDerivatives;
    Vector GaussWeights;
    Matrix ShapeFunctions;
    ShapeFunctionDerivativesArrayType ShapeDerivatives;
    this->CalculateGeometryData(GaussWeights,ShapeFunctions,ShapeDerivatives);
    const long unsigned int NumGauss = GaussWeights.size();

    GeometryType& r_geometry = this->GetGeometry();
    GeometryUtils::ShapeFunctionsSecondDerivativesTransformOnAllIntegrationPoints(
            ShapeFunctionSecondDerivatives,this->GetGeometry(),this->GetIntegrationMethod());

    TElementData data;
    data.Initialize(*this, rCurrentProcessInfo);

    array_1d<double,NumNodes*Dim> momentum_rhs = ZeroVector(NumNodes*Dim);
    VectorType MassRHS = ZeroVector(NumNodes);
    VectorType NodalArea = ZeroVector(NumNodes);

    for (unsigned int g = 0; g < NumGauss; g++)
    {
        this->UpdateIntegrationPointDataSecondDerivatives(data, g, GaussWeights[g], row(ShapeFunctions, g), ShapeDerivatives[g],ShapeFunctionSecondDerivatives[g]);

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
            MassRHS[i] += W*MassRes;
            NodalArea[i] += data.Weight*data.N[i];
        }
    }

    // Add carefully to nodal variables to avoid OpenMP race condition
    for (SizeType i = 0; i < NumNodes; ++i)
    {
        r_geometry[i].SetLock(); // So it is safe to write in the node in OpenMP
        array_1d<double,3>& rMomValue = r_geometry[i].FastGetSolutionStepValue(ADVPROJ);
        unsigned int row = i*Dim;
        for (unsigned int d = 0; d < Dim; d++)
            rMomValue[d] += momentum_rhs[row+d];
        r_geometry[i].FastGetSolutionStepValue(NODAL_AREA) += NodalArea[i];
        r_geometry[i].FastGetSolutionStepValue(DIVPROJ) += MassRHS[i];
        r_geometry[i].UnSetLock(); // Free the node for other threads
    }
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

template< class TElementData >
void QSVMSDEMCoupled<TElementData>::UpdateSubscaleVelocity(
    const TElementData& rData)
{
    array_1d<double,Dim> predicted_subscale_velocity;

    array_1d<double,Dim> previous_velocity = mPreviousVelocity[rData.IntegrationPointIndex];

    //for (size_t i = 0; i < NumNodes; i++) {
    array_1d<double,Dim> subscale_velocity_on_previous_iteration = mPredictedSubscaleVelocity[rData.IntegrationPointIndex];

    array_1d<double,3> v_d = ZeroVector(3);

    for (unsigned int d = 0; d < Dim; d++)
        v_d[d] = previous_velocity[d] + subscale_velocity_on_previous_iteration[d];

    // Part of the residual that does not depend on the subscale
    array_1d<double,3> static_residual = ZeroVector(3);

    if (!rData.UseOSS)
        this->AlgebraicMomentumResidual(rData,v_d,static_residual);
    else
        this->OrthogonalMomentumResidual(rData,v_d,static_residual);


    BoundedMatrix<double,Dim,Dim> tau_one = ZeroMatrix(Dim, Dim);
    double tau_two;

    this->CalculateTau(rData,v_d,tau_one,tau_two);

    for (unsigned int d = 0; d < Dim; d++)
        predicted_subscale_velocity[d] = tau_one(d,d) * static_residual[d];

    noalias(mPredictedSubscaleVelocity[rData.IntegrationPointIndex]) = predicted_subscale_velocity;

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
template class QSVMSDEMCoupled< QSVMSDEMCoupledData<2,6> >;
template class QSVMSDEMCoupled< QSVMSDEMCoupledData<2,9> >;

template class QSVMSDEMCoupled< QSVMSDEMCoupledData<3,8> >;
template class QSVMSDEMCoupled< QSVMSDEMCoupledData<3,27> >;
}