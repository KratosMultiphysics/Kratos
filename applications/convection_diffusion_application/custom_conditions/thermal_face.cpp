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


#include "flux_condition.h"
#include "includes/convection_diffusion_settings.h"
#include "includes/variables.h"
#include "utilities/integration_utilities.h"

namespace Kratos
{

// Public Life Cycle //////////////////////////////////////////////////////////

template<unsigned int TDim, unsigned int TNodesNumber >
ThermalFace<TDim, TNodesNumber>::ThermalFace(IndexType NewId, Geometry< Node<3> >::Pointer pGeometry):
    Condition(NewId,pGeometry)
{
}

template<unsigned int TDim, unsigned int TNodesNumber >
ThermalFace<TDim, TNodesNumber>::ThermalFace(
    IndexType NewId,
    Geometry< Node<3> >::Pointer pGeometry,
    Properties::Pointer pProperties):
    Condition(NewId,pGeometry,pProperties)
{
}

template<unsigned int TDim, unsigned int TNodesNumber >
ThermalFace<TDim, TNodesNumber>::~ThermalFace()
{
}

// Public Operations //////////////////////////////////////////////////////////

template<unsigned int TDim, unsigned int TNodesNumber >
Condition::Pointer ThermalFace<TDim, TNodesNumber>::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    Properties::Pointer pProperties) const
{
    return Kratos::make_shared<ThermalFace<TDim,TNodesNumber>>(NewId, GetGeometry().Create(ThisNodes), pProperties);
}

template<unsigned int TDim, unsigned int TNodesNumber>
Condition::Pointer ThermalFace<TDim, TNodesNumber>::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    Properties::Pointer pProperties) const
{
    return Kratos::make_shared<ThermalFace<TDim,TNodesNumber>>(NewId, pGeom, pProperties);
}

template<unsigned int TDim, unsigned int TNodesNumber>
void ThermalFace<TDim, TNodesNumber>::EquationIdVector(
    EquationIdVectorType& rResult,
    ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Get unknown variable from convection diffusion settings
    const ProcessInfo& r_process_info = rCurrentProcessInfo; // To ensure that Gets are threadsafe
    ConvectionDiffusionSettings &r_conv_diff_settings = *(r_process_info[CONVECTION_DIFFUSION_SETTINGS]);
    const Variable<double> &r_unknown_var = r_conv_diff_settings.GetUnknownVariable();

    // Resize the equation ids. vector
    if (rResult.size() != TNodesNumber){
        rResult.resize(TNodesNumber, false);
    }

    // Fill the equation ids. vector from the condition DOFs
    auto &r_geometry = this->GetGeometry();
    for (unsigned int i = 0; i<TNodesNumber; ++i){
        rResult[i] = r_geometry[i].GetDof(r_unknown_var).EquationId();
    }

    KRATOS_CATCH("")
}

template<unsigned int TDim, unsigned int TNodesNumber>
void ThermalFace<TDim, TNodesNumber>::GetDofList(
    DofsVectorType& rConditionalDofList,
    ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Get unknown variable from convection diffusion settings
    const ProcessInfo &r_process_info = rCurrentProcessInfo; // To ensure that Gets are threadsafe
    ConvectionDiffusionSettings &r_conv_diff_settings = *(r_process_info[CONVECTION_DIFFUSION_SETTINGS]);
    const Variable<double> &r_unknown_var = r_conv_diff_settings.GetUnknownVariable();

    // Resize the DOFs vector
    if (rConditionalDofList.size() != TNodesNumber){
        rConditionalDofList.resize(TNodesNumber);
    }

    // Fill the DOFs vector from the condition nodes
    auto &r_geometry = this->GetGeometry();
    for (unsigned int i = 0; i < TNodesNumber; ++i){
        rConditionalDofList[i] = r_geometry[i].pGetDof(r_unknown_var);
    }

    KRATOS_CATCH("")
}

template<unsigned int TDim, unsigned int TNodesNumber>
void ThermalFace<TDim, TNodesNumber>::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    this->CalculateLeftHandSide(rLeftHandSideMatrix, rCurrentProcessInfo);
    this->CalculateRightHandSide(rRightHandSideVector, rCurrentProcessInfo);

    KRATOS_CATCH("")
}

template<unsigned int TDim, unsigned int TNodesNumber>
void ThermalFace<TDim, TNodesNumber>::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Check (and resize) LHS matrix
    if (rLeftHandSideMatrix.size1() != TNodesNumber || rLeftHandSideMatrix.size2() == TNodesNumber) {
        rLeftHandSideMatrix.resize(TNodesNumber, TNodesNumber, false);
    }

    // Set LHS to zero
    noalias(rLeftHandSideMatrix) = ZeroMatrix(TNodesNumber,TNodesNumber);

    // Declare Gauss pt. data container
    ConditionDataStruct gauss_pt_data;

    // Get geometry data
    auto &r_geometry = this->GetGeometry();
    const auto &r_gauss_pts = r_geometry.IntegrationPoints(this->GetIntegrationMethod());
    const unsigned int n_gauss = r_gauss_pts.size();
    Vector gauss_pts_J_det = ZeroVector(n_gauss);
    r_geometry.DeterminantOfJacobian(gauss_pts_J_det, this->GetIntegrationMethod());
    const MatrixType N_container = r_geometry.ShapeFunctionsValues(this->GetIntegrationMethod());

    // Set user-defined variables values in data container
    const ProcessInfo& r_process_info = rCurrentProcessInfo; // to ensure that Gets are threadsafe
    ConvectionDiffusionSettings &r_conv_diff_settings = *(r_process_info[CONVECTION_DIFFUSION_SETTINGS]);
    const Variable<double> &r_unknown_var = r_conv_diff_settings.GetUnknownVariable();
    const Variable<double> &r_flux_var = r_conv_diff_settings.GetSurfaceSourceVariable();

    for (unsigned int i = 0; i < TNodesNumber; ++i) {
        gauss_pt_data.UnknownValues[i] = r_geometry[i].FastGetSolutionStepValue(r_unknown_var);
        gauss_pt_data.FaceHeatFluxValues[i] = r_geometry[i].FastGetSolutionStepValue(r_flux_var);
    }

    // Fill data container values from properties
    gauss_pt_data.Emissivity = (this->GetProperties())[EMISSIVITY];
    gauss_pt_data.AmbientTemperature = (this->GetProperties())[AMBIENT_TEMPERATURE];
    gauss_pt_data.ConvectionCoefficient = (this->GetProperties())[CONVECTION_COEFFICIENT];

    // Gauss pts. loop
    for (unsigned int g = 0; g < n_gauss; g++) {
        gauss_pt_data.N = row(N_container, g);
        gauss_pt_data.Weight = gauss_pts_J_det[g] * r_gauss_pts[g].Weight();
        this->AddIntegrationPointLHSContribution(rLeftHandSideMatrix, gauss_pt_data);
    }

    KRATOS_CATCH("")
}

template<unsigned int TDim, unsigned int TNodesNumber>
void ThermalFace<TDim, TNodesNumber>::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Check (and resize) RHS vecto
    if (rRightHandSideVector.size() != TNodesNumber) {
        rRightHandSideVector.resize(TNodesNumber,false);
    }

    // Initialize RHS vector
    noalias(rRightHandSideVector) = ZeroVector(TNodesNumber);

    // Declare Gauss pt. data container
    ConditionDataStruct gauss_pt_data;

    // Get geometry data
    auto &r_geometry = this->GetGeometry();
    const auto &r_gauss_pts = r_geometry.IntegrationPoints(this->GetIntegrationMethod());
    const unsigned int n_gauss = r_gauss_pts.size();
    Vector gauss_pts_J_det = ZeroVector(n_gauss);
    r_geometry.DeterminantOfJacobian(gauss_pts_J_det, this->GetIntegrationMethod());
    const MatrixType N_container = r_geometry.ShapeFunctionsValues(this->GetIntegrationMethod());

    // Set user-defined variables values in data container
    const ProcessInfo& r_process_info = rCurrentProcessInfo; // to ensure that Gets are threadsafe
    ConvectionDiffusionSettings &r_conv_diff_settings = *(r_process_info[CONVECTION_DIFFUSION_SETTINGS]);
    const Variable<double> &r_unknown_var = r_conv_diff_settings.GetUnknownVariable();
    const Variable<double> &r_flux_var = r_conv_diff_settings.GetSurfaceSourceVariable();

    for (unsigned int i = 0; i < TNodesNumber; ++i) {
        gauss_pt_data.UnknownValues[i] = r_geometry[i].FastGetSolutionStepValue(r_unknown_var);
        gauss_pt_data.FaceHeatFluxValues[i] = r_geometry[i].FastGetSolutionStepValue(r_flux_var);
    }

    // Fill data container values from properties
    gauss_pt_data.Emissivity = (this->GetProperties())[EMISSIVITY];
    gauss_pt_data.AmbientTemperature = (this->GetProperties())[AMBIENT_TEMPERATURE];
    gauss_pt_data.ConvectionCoefficient = (this->GetProperties())[CONVECTION_COEFFICIENT];

    // Gauss pts. loop
    for (unsigned int g = 0; g < n_gauss; g++) {
        gauss_pt_data.N = row(N_container, g);
        gauss_pt_data.Weight = gauss_pts_J_det[g] * r_gauss_pts[g].Weight();
        this->AddIntegrationPointRHSContribution(rRightHandSideVector, gauss_pt_data);
    }

    KRATOS_CATCH("")
}


template<unsigned int TDim, unsigned int TNodesNumber>
inline GeometryData::IntegrationMethod ThermalFace<TDim, TNodesNumber>::GetIntegrationMethod()
{
    return IntegrationUtilities::GetIntegrationMethodForExactMassMatrixEvaluation(this->GetGeometry());
}

template<unsigned int TDim, unsigned int TNodesNumber>
void ThermalFace<TDim, TNodesNumber>::GetValueOnIntegrationPoints(
    const Variable<array_1d<double,3> > &rVariable,
    std::vector<array_1d<double,3> > &rValues,
    const ProcessInfo &rCurrentProcessInfo)
{
    const unsigned int NumGauss = this->GetGeometry().IntegrationPointsNumber(this->GetIntegrationMethod());
    rValues.resize(NumGauss);
    if (rVariable == NORMAL)
    {
        this->CalculateNormal(rValues[0]);
    }
    else
    {
        /* The cast is done to avoid modification of the element's data. Data modification
         * would happen if rVariable is not stored now (would initialize a pointer to &rVariable
         * with associated value of 0.0). This is catastrophic if the variable referenced
         * goes out of scope.
         */
        const ThermalFace* const_this = static_cast< const ThermalFace* >(this);
        rValues[0] = const_this->GetValue(rVariable);
    }

    // Copy the values to the different gauss points
    for (unsigned int g = 1; g < NumGauss; g++)
    {
        noalias(rValues[g]) = rValues[0];
    }
}

template<unsigned int TDim, unsigned int TNodesNumber>
void ThermalFace<TDim, TNodesNumber>::GetValueOnIntegrationPoints(
    const Variable<double>& rVariable,
    std::vector<double>& rValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    const unsigned int NumGauss = this->GetGeometry().IntegrationPointsNumber(this->GetIntegrationMethod());
    rValues.resize(NumGauss);
    const ThermalFace* const_this = static_cast< const ThermalFace* >(this);
    rValues[0] = const_this->GetValue(rVariable);
    for (unsigned int g = 1; g < NumGauss; g++)
    {
        rValues[g] = rValues[0];
    }
}

template<unsigned int TDim, unsigned int TNodesNumber>
void ThermalFace<TDim, TNodesNumber>::GetValueOnIntegrationPoints(
    const Variable<array_1d<double, 6 > >& rVariable,
    std::vector<array_1d<double, 6 > >& rValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    const unsigned int NumGauss = this->GetGeometry().IntegrationPointsNumber(this->GetIntegrationMethod());
    rValues.resize(NumGauss);
    const ThermalFace* const_this = static_cast< const ThermalFace* >(this);
    rValues[0] = const_this->GetValue(rVariable);
    for (unsigned int g = 1; g < NumGauss; g++)
    {
        noalias(rValues[g]) = rValues[0];
    }
}

template<unsigned int TDim, unsigned int TNodesNumber>
void ThermalFace<TDim, TNodesNumber>::GetValueOnIntegrationPoints(
    const Variable<Vector>& rVariable,
    std::vector<Vector>& rValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    const unsigned int NumGauss = this->GetGeometry().IntegrationPointsNumber(this->GetIntegrationMethod());
    rValues.resize(NumGauss);
    const ThermalFace* const_this = static_cast< const ThermalFace* >(this);
    rValues[0] = const_this->GetValue(rVariable);
    for (unsigned int g = 1; g < NumGauss; g++)
    {
        noalias(rValues[g]) = rValues[0];
    }
}

template<unsigned int TDim, unsigned int TNodesNumber>
void ThermalFace<TDim, TNodesNumber>::GetValueOnIntegrationPoints(
    const Variable<Matrix>& rVariable,
    std::vector<Matrix>& rValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    const unsigned int NumGauss = this->GetGeometry().IntegrationPointsNumber(this->GetIntegrationMethod());
    rValues.resize(NumGauss);
    const ThermalFace* const_this = static_cast< const ThermalFace* >(this);
    rValues[0] = const_this->GetValue(rVariable);
    for (unsigned int g = 1; g < NumGauss; g++)
    {
        noalias(rValues[g]) = rValues[0];
    }
}

// Input and Output ///////////////////////////////////////////////////////////

template<unsigned int TDim, unsigned int TNodesNumber>
std::string ThermalFace<TDim, TNodesNumber>::Info() const
{
    std::stringstream buffer;
    buffer << "ThermalFace #" << Id();
    return buffer.str();
}

template<unsigned int TDim, unsigned int TNodesNumber>
void ThermalFace<TDim, TNodesNumber>::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "ThermalFace #" << Id();
}

template<unsigned int TDim, unsigned int TNodesNumber>
void ThermalFace<TDim, TNodesNumber>::PrintData(std::ostream& rOStream) const
{
    rOStream << "ThermalFace #" << Id() << std::endl;
    this->GetGeometry().PrintData(rOStream);
}


// Finite element functions ///////////////////////////////////////////////////

template<unsigned int TDim, unsigned int TNodesNumber>
void ThermalFace<TDim, TNodesNumber>::AddIntegrationPointRHSContribution(
    VectorType& rRightHandSideVector,
    const ConditionDataStruct &rData)
{
    const double stefen_boltzmann = 5.67e-8;
    const double gauss_pt_unknown = rData.GaussPointUnknown();
    const double gauss_pt_flux = rData.GaussPointFaceHeatFlux();
    const double aux_rad_rhs = rData.Emissivity * stefen_boltzmann * (std::pow(gauss_pt_unknown, 4) - pow(rData.AmbientTemperature, 4));
    const double aux_conv_rhs = rData.ConvectionCoefficient * (gauss_pt_unknown - rData.AmbientTemperature);
    for (unsigned int i = 0; i < TNodesNumber; ++i) {
        // Add external face heat flux contribution
        rRightHandSideVector[i] += rData.N(i) * gauss_pt_flux * rData.Weight;
        // Add ambient radiation contribution
        rRightHandSideVector[i] -= rData.N(i) * aux_rad_rhs * rData.Weight;
        // Add ambient temperature convection contribution
        rRightHandSideVector[i] -= rData.N(i) * aux_conv_rhs * rData.Weight;
    }
}

template<unsigned int TDim, unsigned int TNodesNumber>
void ThermalFace<TDim, TNodesNumber>::AddIntegrationPointLHSContribution(
    MatrixType& rLeftHandSideMatrix,
    const ConditionDataStruct &rData)
{
    const double stefen_boltzmann = 5.67e-8;
    const double aux_rad_lhs = rData.Emissivity * stefen_boltzmann * 4.0 * std::pow(rData.GaussPointUnknown(), 3);
    for (unsigned int i = 0; i < TNodesNumber; ++i) {
        for (unsigned int j = 0; j < TNodesNumber; ++j) {
            // Add ambient radiation contribution
            rLeftHandSideMatrix(i,j) += rData.N(i) * aux_rad_lhs * rData.N(j) * rData.Weight;
            // Ambient temperature convection contribution
            rLeftHandSideMatrix(i,j) += rData.N(i) * rData.ConvectionCoefficient * rData.N(j) * rData.Weight;
        }
    }
}

template <>
void ThermalFace<2,2>::CalculateNormal(array_1d<double,3> &rNormal)
{
    const auto &r_geometry = this->GetGeometry();
    rNormal[0] = r_geometry[1].Y() - r_geometry[0].Y();
    rNormal[1] = r_geometry[0].X() - r_geometry[1].X();
    rNormal[2] = 0.0;
}


template <>
void ThermalFace<3,3>::CalculateNormal(array_1d<double,3> &rNormal)
{
    const auto &r_geometry = this->GetGeometry();

    const array_1d<double,3> v1 = r_geometry[1] - r_geometry[0];
    const array_1d<double,3> v2 = r_geometry[2] - r_geometry[0];

    MathUtils<double>::CrossProduct(rNormal, v1, v2);
    rNormal *= 0.5;
}

template <>
void ThermalFace<3,4>::CalculateNormal(array_1d<double,3> &rNormal)
{
    KRATOS_ERROR << "This function is not yet implemented" << std::endl;

}

// Serialization //////////////////////////////////////////////////////////////

template<unsigned int TDim, unsigned int TNodesNumber>
ThermalFace<TDim, TNodesNumber>::ThermalFace():
    Condition()
{
}

template<unsigned int TDim, unsigned int TNodesNumber>
void ThermalFace<TDim, TNodesNumber>::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Condition);
}

template<unsigned int TDim, unsigned int TNodesNumber>
void ThermalFace<TDim, TNodesNumber>::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Condition);
}

template class ThermalFace<2,2>;
template class ThermalFace<3,3>;
template class ThermalFace<3,4>;

}
