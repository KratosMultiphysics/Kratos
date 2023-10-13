// KRATOS ___ ___  _  ___   __   ___ ___ ___ ___
//       / __/ _ \| \| \ \ / /__|   \_ _| __| __|
//      | (_| (_) | .` |\ V /___| |) | || _|| _|
//       \___\___/|_|\_| \_/    |___/___|_| |_|  APPLICATION
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela
//                   Riccardo Rossi
//                   Ruben Zorrilla
//


#include "thermal_face.h"
#include "utilities/integration_utilities.h"
#include "includes/convection_diffusion_settings.h"

namespace Kratos
{

// Public Life Cycle //////////////////////////////////////////////////////////

ThermalFace::ThermalFace(IndexType NewId, Geometry< Node >::Pointer pGeometry):
    Condition(NewId,pGeometry)
{
}

ThermalFace::ThermalFace(
    IndexType NewId,
    Geometry< Node >::Pointer pGeometry,
    Properties::Pointer pProperties):
    Condition(NewId,pGeometry,pProperties)
{
}

ThermalFace::~ThermalFace()
{
}

// Public Operations //////////////////////////////////////////////////////////

Condition::Pointer ThermalFace::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<ThermalFace>(NewId, GetGeometry().Create(ThisNodes), pProperties);
}

Condition::Pointer ThermalFace::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<ThermalFace>(NewId, pGeom, pProperties);
}

void ThermalFace::EquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    // Get unknown variable from convection diffusion settings
    const ProcessInfo& r_process_info = rCurrentProcessInfo; // To ensure that Gets are threadsafe
    ConvectionDiffusionSettings &r_conv_diff_settings = *(r_process_info[CONVECTION_DIFFUSION_SETTINGS]);
    const Variable<double> &r_unknown_var = r_conv_diff_settings.GetUnknownVariable();

    // Resize the equation ids. vector
    auto &r_geometry = this->GetGeometry();
    const unsigned int n_nodes = r_geometry.PointsNumber();
    if (rResult.size() != n_nodes) {
        rResult.resize(n_nodes, false);
    }

    // Fill the equation ids. vector from the condition DOFs
    for (unsigned int i = 0; i < n_nodes; ++i){
        rResult[i] = r_geometry[i].GetDof(r_unknown_var).EquationId();
    }

    KRATOS_CATCH("")
}

void ThermalFace::GetDofList(
    DofsVectorType& rConditionalDofList,
    const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    // Get unknown variable from convection diffusion settings
    const ProcessInfo &r_process_info = rCurrentProcessInfo; // To ensure that Gets are threadsafe
    ConvectionDiffusionSettings &r_conv_diff_settings = *(r_process_info[CONVECTION_DIFFUSION_SETTINGS]);
    const Variable<double> &r_unknown_var = r_conv_diff_settings.GetUnknownVariable();

    // Resize the DOFs vector
    auto &r_geometry = this->GetGeometry();
    const unsigned int n_nodes = r_geometry.PointsNumber();
    if (rConditionalDofList.size() != n_nodes){
        rConditionalDofList.resize(n_nodes);
    }

    // Fill the DOFs vector from the condition nodes
    for (unsigned int i = 0; i < n_nodes; ++i){
        rConditionalDofList[i] = r_geometry[i].pGetDof(r_unknown_var);
    }

    KRATOS_CATCH("")
}

void ThermalFace::SetIntegrationWeight(
    const IndexType IntegrationPointIndex,
    const typename GeometryType::IntegrationPointsArrayType &rIntegrationPoints,
    const Vector &rJacobianDeterminantsVector,
    ConditionDataStruct &rData)
{
    rData.Weight = rJacobianDeterminantsVector[IntegrationPointIndex] * rIntegrationPoints[IntegrationPointIndex].Weight();
}

void ThermalFace::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    this->CalculateLeftHandSide(rLeftHandSideMatrix, rCurrentProcessInfo);
    this->CalculateRightHandSide(rRightHandSideVector, rCurrentProcessInfo);

    KRATOS_CATCH("")
}

void ThermalFace::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Check (and resize) LHS matrix
    auto &r_geometry = this->GetGeometry();
    const unsigned int n_nodes = r_geometry.PointsNumber();
    if (rLeftHandSideMatrix.size1() != n_nodes || rLeftHandSideMatrix.size2() != n_nodes) {
        rLeftHandSideMatrix.resize(n_nodes, n_nodes, false);
    }

    // Set LHS to zero
    noalias(rLeftHandSideMatrix) = ZeroMatrix(n_nodes,n_nodes);

    // Declare Gauss pt. data container
    ConditionDataStruct gauss_pt_data;
    this->FillConditionDataStructure(rCurrentProcessInfo, gauss_pt_data);

    // Get geometry data
    const auto &r_gauss_pts = r_geometry.IntegrationPoints(this->GetIntegrationMethod());
    const unsigned int n_gauss = r_gauss_pts.size();
    Vector gauss_pts_J_det = ZeroVector(n_gauss);
    r_geometry.DeterminantOfJacobian(gauss_pts_J_det, this->GetIntegrationMethod());
    const MatrixType N_container = r_geometry.ShapeFunctionsValues(this->GetIntegrationMethod());

    // Gauss pts. loop
    for (unsigned int g = 0; g < n_gauss; g++) {
        gauss_pt_data.N = row(N_container, g);
        SetIntegrationWeight(g, r_gauss_pts, gauss_pts_J_det, gauss_pt_data);
        this->AddIntegrationPointLHSContribution(rLeftHandSideMatrix, gauss_pt_data);
    }

    KRATOS_CATCH("")
}

void ThermalFace::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Check (and resize) RHS vector
    auto &r_geometry = this->GetGeometry();
    const unsigned int n_nodes = r_geometry.PointsNumber();
    if (rRightHandSideVector.size() != n_nodes) {
        rRightHandSideVector.resize(n_nodes,false);
    }

    // Initialize RHS vector
    noalias(rRightHandSideVector) = ZeroVector(n_nodes);

    // Declare Gauss pt. data container
    ConditionDataStruct gauss_pt_data;
    this->FillConditionDataStructure(rCurrentProcessInfo, gauss_pt_data);

    // Get geometry data
    const auto &r_gauss_pts = r_geometry.IntegrationPoints(this->GetIntegrationMethod());
    const unsigned int n_gauss = r_gauss_pts.size();
    Vector gauss_pts_J_det = ZeroVector(n_gauss);
    r_geometry.DeterminantOfJacobian(gauss_pts_J_det, this->GetIntegrationMethod());
    const MatrixType N_container = r_geometry.ShapeFunctionsValues(this->GetIntegrationMethod());

    // Gauss pts. loop
    for (unsigned int g = 0; g < n_gauss; g++) {
        gauss_pt_data.N = row(N_container, g);
        SetIntegrationWeight(g, r_gauss_pts, gauss_pts_J_det, gauss_pt_data);
        this->AddIntegrationPointRHSContribution(rRightHandSideVector, gauss_pt_data);
    }

    KRATOS_CATCH("")
}


inline GeometryData::IntegrationMethod ThermalFace::GetIntegrationMethod() const
{
    return IntegrationUtilities::GetIntegrationMethodForExactMassMatrixEvaluation(this->GetGeometry());
}

void ThermalFace::CalculateOnIntegrationPoints(
    const Variable<array_1d<double,3> > &rVariable,
    std::vector<array_1d<double,3> > &rValues,
    const ProcessInfo &rCurrentProcessInfo)
{
    const unsigned int n_gauss = this->GetGeometry().IntegrationPointsNumber(this->GetIntegrationMethod());
    rValues.resize(n_gauss);
    if (rVariable == NORMAL) {
        const auto &r_geometry = this->GetGeometry();
        const auto &r_gauss_pts = r_geometry.IntegrationPoints(this->GetIntegrationMethod());
        for (unsigned int g = 0; g < n_gauss; ++g) {
            rValues[g] = r_geometry.Normal(r_gauss_pts[g].Coordinates());
        }
    } else {
        /* The cast is done to avoid modification of the element's data. Data modification
         * would happen if rVariable is not stored now (would initialize a pointer to &rVariable
         * with associated value of 0.0). This is catastrophic if the variable referenced
         * goes out of scope.
         */
        const ThermalFace* const_this = static_cast< const ThermalFace* >(this);
        rValues[0] = const_this->GetValue(rVariable);
        // Copy the values to the different gauss points
        for (unsigned int g = 1; g < n_gauss; ++g) {
            noalias(rValues[g]) = rValues[0];
        }
    }
}

void ThermalFace::CalculateOnIntegrationPoints(
    const Variable<double>& rVariable,
    std::vector<double>& rValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    const unsigned int n_gauss = this->GetGeometry().IntegrationPointsNumber(this->GetIntegrationMethod());
    rValues.resize(n_gauss);
    const ThermalFace* const_this = static_cast< const ThermalFace* >(this);
    rValues[0] = const_this->GetValue(rVariable);
    for (unsigned int g = 1; g < n_gauss; ++g) {
        rValues[g] = rValues[0];
    }
}

void ThermalFace::CalculateOnIntegrationPoints(
    const Variable<array_1d<double, 6 > >& rVariable,
    std::vector<array_1d<double, 6 > >& rValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    const unsigned int n_gauss = this->GetGeometry().IntegrationPointsNumber(this->GetIntegrationMethod());
    rValues.resize(n_gauss);
    const ThermalFace* const_this = static_cast< const ThermalFace* >(this);
    rValues[0] = const_this->GetValue(rVariable);
    for (unsigned int g = 1; g < n_gauss; ++g) {
        noalias(rValues[g]) = rValues[0];
    }
}

void ThermalFace::CalculateOnIntegrationPoints(
    const Variable<Vector>& rVariable,
    std::vector<Vector>& rValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    const unsigned int n_gauss = this->GetGeometry().IntegrationPointsNumber(this->GetIntegrationMethod());
    rValues.resize(n_gauss);
    const ThermalFace* const_this = static_cast< const ThermalFace* >(this);
    rValues[0] = const_this->GetValue(rVariable);
    for (unsigned int g = 1; g < n_gauss; ++g) {
        noalias(rValues[g]) = rValues[0];
    }
}

void ThermalFace::CalculateOnIntegrationPoints(
    const Variable<Matrix>& rVariable,
    std::vector<Matrix>& rValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    const unsigned int n_gauss = this->GetGeometry().IntegrationPointsNumber(this->GetIntegrationMethod());
    rValues.resize(n_gauss);
    const ThermalFace* const_this = static_cast< const ThermalFace* >(this);
    rValues[0] = const_this->GetValue(rVariable);
    for (unsigned int g = 1; g < n_gauss; ++g) {
        noalias(rValues[g]) = rValues[0];
    }
}

// Input and Output ///////////////////////////////////////////////////////////

std::string ThermalFace::Info() const
{
    std::stringstream buffer;
    buffer << "ThermalFace #" << Id();
    return buffer.str();
}

void ThermalFace::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "ThermalFace #" << Id();
}

void ThermalFace::PrintData(std::ostream& rOStream) const
{
    rOStream << "ThermalFace #" << Id() << std::endl;
    this->GetGeometry().PrintData(rOStream);
}


// Finite element functions ///////////////////////////////////////////////////

void ThermalFace::AddIntegrationPointRHSContribution(
    VectorType& rRightHandSideVector,
    const ConditionDataStruct &rData)
{
    const double gauss_pt_unknown = rData.GaussPointUnknown();
    const double gauss_pt_flux = rData.GaussPointFaceHeatFlux();
    const double aux_rad_rhs = rData.Emissivity * StefanBoltzmann * (std::pow(gauss_pt_unknown, 4) - pow(rData.AmbientTemperature, 4));
    const double aux_conv_rhs = rData.ConvectionCoefficient * (gauss_pt_unknown - rData.AmbientTemperature);
    for (unsigned int i = 0; i < (this->GetGeometry()).PointsNumber(); ++i) {
        // Add external face heat flux contribution
        rRightHandSideVector[i] += rData.N(i) * gauss_pt_flux * rData.Weight;
        // Add ambient radiation contribution
        rRightHandSideVector[i] -= rData.N(i) * aux_rad_rhs * rData.Weight;
        // Add ambient temperature convection contribution
        rRightHandSideVector[i] -= rData.N(i) * aux_conv_rhs * rData.Weight;
    }
}

void ThermalFace::AddIntegrationPointLHSContribution(
    MatrixType& rLeftHandSideMatrix,
    const ConditionDataStruct &rData)
{
    const unsigned int n_nodes = (this->GetGeometry()).PointsNumber();
    const double aux_rad_lhs = rData.Emissivity * StefanBoltzmann * 4.0 * std::pow(rData.GaussPointUnknown(), 3);
    for (unsigned int i = 0; i < n_nodes; ++i) {
        for (unsigned int j = 0; j < n_nodes; ++j) {
            // Add ambient radiation contribution
            rLeftHandSideMatrix(i,j) += rData.N(i) * aux_rad_lhs * rData.N(j) * rData.Weight;
            // Ambient temperature convection contribution
            rLeftHandSideMatrix(i,j) += rData.N(i) * rData.ConvectionCoefficient * rData.N(j) * rData.Weight;
        }
    }
}

void ThermalFace::FillConditionDataStructure(
    const ProcessInfo &rCurrentProcessInfo,
    ConditionDataStruct &rData)
{
    // Set user-defined variables values in data container
    ConvectionDiffusionSettings &r_conv_diff_settings = *(rCurrentProcessInfo[CONVECTION_DIFFUSION_SETTINGS]);
    const Variable<double> &r_unknown_var = r_conv_diff_settings.GetUnknownVariable();
    const Variable<double> &r_flux_var = r_conv_diff_settings.GetSurfaceSourceVariable();
    const auto &r_geometry = this->GetGeometry();
    const unsigned int n_nodes = r_geometry.PointsNumber();
    rData.UnknownValues.resize(n_nodes, false);
    rData.FaceHeatFluxValues.resize(n_nodes, false);

    for (unsigned int i = 0; i < n_nodes; ++i) {
        rData.UnknownValues[i] = r_geometry[i].FastGetSolutionStepValue(r_unknown_var);
        rData.FaceHeatFluxValues[i] = r_geometry[i].FastGetSolutionStepValue(r_flux_var);
    }

    // Check (and resize) and fill data container arrays
    if (rData.UnknownValues.size() != n_nodes) {
        rData.UnknownValues.resize(n_nodes, false);
    }
    if (rData.FaceHeatFluxValues.size() != n_nodes) {
        rData.FaceHeatFluxValues.resize(n_nodes, false);
    }

    // Fill data container values from properties
    // Const reference is required to have thread-safe access
    const auto &r_prop = this->GetProperties();
    rData.Emissivity = r_prop[EMISSIVITY];
    rData.AmbientTemperature = r_prop[AMBIENT_TEMPERATURE];
    rData.ConvectionCoefficient = r_prop[CONVECTION_COEFFICIENT];
}

// Serialization //////////////////////////////////////////////////////////////

ThermalFace::ThermalFace():
    Condition()
{
}

void ThermalFace::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Condition);
}

void ThermalFace::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Condition);
}

}
