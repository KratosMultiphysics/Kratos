// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:     BSD License
//           license: structural_mechanics_application/license.txt
//
//  Main authors: Klaus B. Sautter
//
//
//

// System includes

// External includes

// Project includes
#include "custom_elements/truss_element_linear_3D2N.hpp"
#include "includes/define.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos {
TrussElementLinear3D2N::TrussElementLinear3D2N(IndexType NewId,
                                               GeometryType::Pointer pGeometry)
    : TrussElement3D2N(NewId, pGeometry) {}

TrussElementLinear3D2N::TrussElementLinear3D2N(
    IndexType NewId, GeometryType::Pointer pGeometry,
    PropertiesType::Pointer pProperties)
    : TrussElement3D2N(NewId, pGeometry, pProperties) {}

Element::Pointer
TrussElementLinear3D2N::Create(IndexType NewId,
                               NodesArrayType const &rThisNodes,
                               PropertiesType::Pointer pProperties) const {
  const GeometryType &rGeom = this->GetGeometry();
  return Kratos::make_shared<TrussElementLinear3D2N>(
      NewId, rGeom.Create(rThisNodes), pProperties);
}

TrussElementLinear3D2N::~TrussElementLinear3D2N() {}

BoundedMatrix<double, TrussElement3D2N::msLocalSize,
               TrussElement3D2N::msLocalSize>
TrussElementLinear3D2N::CreateElementStiffnessMatrix(
    ProcessInfo &rCurrentProcessInfo) {

  KRATOS_TRY
  BoundedMatrix<double, msLocalSize, msLocalSize> LocalStiffnessMatrix =
      ZeroMatrix(msLocalSize, msLocalSize);
  this->CalculateElasticStiffnessMatrix(LocalStiffnessMatrix,
                                        rCurrentProcessInfo);
  
  // test solution -> do this directly in CalculateElasticStiffnessMatrix!!!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  if(this->CheckIfIsPlasticRegime()) {
    const double youngs_modulus = this->GetProperties()[YOUNG_MODULUS];
    const double plastic_tangent_modulus = this->ReturnElastoPlasticTangentModulus();
    for (SizeType i=0;i<msLocalSize;++i){
      for (SizeType j=0;j<msLocalSize;++j){
        LocalStiffnessMatrix(i,j) = LocalStiffnessMatrix(i,j)*(plastic_tangent_modulus/youngs_modulus);
      }
    }
  }
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    

  return LocalStiffnessMatrix;
  KRATOS_CATCH("")
}

void TrussElementLinear3D2N::CalculateLocalSystem(
    MatrixType &rLeftHandSideMatrix, VectorType &rRightHandSideVector,
    ProcessInfo &rCurrentProcessInfo) {

  KRATOS_TRY;
  BoundedVector<double, msLocalSize> internal_forces = ZeroVector(msLocalSize);
  this->UpdateInternalForces(internal_forces);

  // resizing the matrices + create memory for LHS
  rLeftHandSideMatrix = ZeroMatrix(msLocalSize, msLocalSize);
  // creating LHS
  rLeftHandSideMatrix = this->CreateElementStiffnessMatrix(rCurrentProcessInfo);

  //Vector nodal_deformation = ZeroVector(msLocalSize);
  //this->GetValuesVector(nodal_deformation, 0);

  rRightHandSideVector = ZeroVector(msLocalSize);
  noalias(rRightHandSideVector) -= internal_forces;
  noalias(rRightHandSideVector) += this->CalculateBodyForces();
  this->AddPrestressLinear(rRightHandSideVector);

  KRATOS_CATCH("")
}

void TrussElementLinear3D2N::AddPrestressLinear(
    VectorType &rRightHandSideVector) {
  KRATOS_TRY;
  BoundedMatrix<double, msLocalSize, msLocalSize> transformation_matrix =
      ZeroMatrix(msLocalSize, msLocalSize);
  this->CreateTransformationMatrix(transformation_matrix);
  double prestress = 0.00;
  if (this->GetProperties().Has(TRUSS_PRESTRESS_PK2)) {
    prestress = this->GetProperties()[TRUSS_PRESTRESS_PK2];
  }
  const double A = this->GetProperties()[CROSS_AREA];
  const double N = prestress * A;

  // internal force vectors
  BoundedVector<double, msLocalSize> f_local = ZeroVector(msLocalSize);
  f_local[0] = -1.00 * N;
  f_local[3] = 1.00 * N;
  rRightHandSideVector -= prod(transformation_matrix, f_local);
  KRATOS_CATCH("")
}

void TrussElementLinear3D2N::CalculateRightHandSide(
    VectorType &rRightHandSideVector, ProcessInfo &rCurrentProcessInfo) {

  KRATOS_TRY

  std::cout << "CalculateRightHandSide" <<std::endl;

  Matrix left_hand_side_matrix = ZeroMatrix(msLocalSize, msLocalSize);
  this->CalculateLeftHandSide(left_hand_side_matrix, rCurrentProcessInfo);
  Vector nodal_deformation = ZeroVector(msLocalSize);
  this->GetValuesVector(nodal_deformation);
  rRightHandSideVector = ZeroVector(msLocalSize);
  noalias(rRightHandSideVector) -=
      prod(left_hand_side_matrix, nodal_deformation);
  this->AddPrestressLinear(rRightHandSideVector);

  // add bodyforces
  noalias(rRightHandSideVector) += this->CalculateBodyForces();
  KRATOS_CATCH("")
}

void TrussElementLinear3D2N::CalculateLeftHandSide(
    MatrixType &rLeftHandSideMatrix, ProcessInfo &rCurrentProcessInfo) {

  KRATOS_TRY

  // resizing the matrices + create memory for LHS
  rLeftHandSideMatrix = ZeroMatrix(msLocalSize, msLocalSize);
  // creating LHS
  rLeftHandSideMatrix = this->CreateElementStiffnessMatrix(rCurrentProcessInfo);
  KRATOS_CATCH("")
}

void TrussElementLinear3D2N::CalculateOnIntegrationPoints(
    const Variable<array_1d<double, 3>> &rVariable,
    std::vector<array_1d<double, 3>> &rOutput,
    const ProcessInfo &rCurrentProcessInfo) {

  const GeometryType::IntegrationPointsArrayType &integration_points =
      GetGeometry().IntegrationPoints();
  if (rOutput.size() != integration_points.size()) {
    rOutput.resize(integration_points.size());
  }

  if (rVariable == FORCE) {
    BoundedVector<double, msDimension> truss_forces = ZeroVector(msDimension);
    truss_forces[2] = 0.00;
    truss_forces[1] = 0.00;
    const double A = this->GetProperties()[CROSS_AREA];

    double prestress = 0.00;
    if (this->GetProperties().Has(TRUSS_PRESTRESS_PK2)) {
      prestress = this->GetProperties()[TRUSS_PRESTRESS_PK2];
    }

    Matrix left_hand_side_matrix = ZeroMatrix(msLocalSize, msLocalSize);
    ProcessInfo
        dummy_info; // CalculateLeftHandSide does not take const ProcessInfo
    this->CalculateLeftHandSide(left_hand_side_matrix, dummy_info);
    Vector nodal_deformation = ZeroVector(msLocalSize);
    this->GetValuesVector(nodal_deformation);
    BoundedMatrix<double, msLocalSize, msLocalSize> transformation_matrix =
        ZeroMatrix(msLocalSize, msLocalSize);
    this->CreateTransformationMatrix(transformation_matrix);
    Vector f_int = prod(left_hand_side_matrix, nodal_deformation);
    f_int = prod(Matrix(trans(transformation_matrix)), f_int);
    truss_forces[0] = f_int[3] + prestress * A;

    rOutput[0] = truss_forces;
  }
}

void TrussElementLinear3D2N::WriteTransformationCoordinates(
    BoundedVector<double, TrussElement3D2N::msLocalSize>
        &rReferenceCoordinates) {
  KRATOS_TRY;
  rReferenceCoordinates = ZeroVector(msLocalSize);
  rReferenceCoordinates[0] = this->GetGeometry()[0].X0();
  rReferenceCoordinates[1] = this->GetGeometry()[0].Y0();
  rReferenceCoordinates[2] = this->GetGeometry()[0].Z0();
  rReferenceCoordinates[3] = this->GetGeometry()[1].X0();
  rReferenceCoordinates[4] = this->GetGeometry()[1].Y0();
  rReferenceCoordinates[5] = this->GetGeometry()[1].Z0();
  KRATOS_CATCH("");
}

double TrussElementLinear3D2N::CalculateLinearStrain()  {
  KRATOS_TRY;
  const double length_0 = this->CalculateReferenceLength();
  const double length = this->CalculateCurrentLength();
  const double e = ((length*length-length_0*length_0)/(length+length_0))/length_0;
  return e;
  KRATOS_CATCH("");
}

double TrussElementLinear3D2N::ReturnElastoPlasticTangentModulus() const {
  KRATOS_TRY;
  double hardening_modulus = 0.00;
  if (this->GetProperties().Has(DP_K)) {
    hardening_modulus = this->GetProperties()[DP_K];
  }
  const double youngs_modulus = this->GetProperties()[YOUNG_MODULUS];
  const double tangent_modulus = (hardening_modulus*youngs_modulus)/(hardening_modulus+youngs_modulus);
  return tangent_modulus;
  KRATOS_CATCH("");
}

double TrussElementLinear3D2N::TrialStateStress() {
  KRATOS_TRY;
  const double elastic_trial_strain = this->CalculateLinearStrain()-this->plastic_strain;
  const double youngs_modulus = this->GetProperties()[YOUNG_MODULUS];
  KRATOS_WATCH(elastic_trial_strain);
  KRATOS_WATCH(this->plastic_strain);
  return (youngs_modulus*elastic_trial_strain);
  KRATOS_CATCH("");
}

double TrussElementLinear3D2N::TrialYieldFunction() {
  KRATOS_TRY;
  double yield_stress_limit = 0.00;
  if (this->GetProperties().Has(INFINITY_YIELD_STRESS)) {
    yield_stress_limit = this->GetProperties()[INFINITY_YIELD_STRESS];
  }
  double hardening_modulus = 0.00;
  if (this->GetProperties().Has(DP_K)) {
    hardening_modulus = this->GetProperties()[DP_K];
  }

  const double trial_stress = this->TrialStateStress();
  double trial_yield_function =  std::abs(trial_stress);
  trial_yield_function -= yield_stress_limit + (hardening_modulus*this->plastic_alpha);
  std::cout << trial_yield_function << '=' << std::abs(trial_stress) << '-' << yield_stress_limit << '+'  << hardening_modulus <<'*'<<this->plastic_alpha << std::endl;

  return trial_yield_function;
  KRATOS_CATCH("");  
}

bool TrussElementLinear3D2N::CheckIfIsPlasticRegime() {
  KRATOS_TRY;
  const double numerical_limit = std::numeric_limits<double>::epsilon();
  bool is_in_plastic_regime = false;
  const double trial_yield_function = this->TrialYieldFunction();
  if (trial_yield_function > 0.00) is_in_plastic_regime = true;

  if (std::abs(this->CalculateLinearStrain())<numerical_limit) is_in_plastic_regime=false;

  return is_in_plastic_regime;
  KRATOS_CATCH("");
}

void TrussElementLinear3D2N::UpdateInternalForces(BoundedVector<double,msLocalSize>& rinternalForces)
{
  KRATOS_TRY;
  ProcessInfo temp_process_information;
  rinternalForces = ZeroVector(msLocalSize);
  const double trial_stress = this->TrialStateStress();
  double current_stress = trial_stress;

  KRATOS_WATCH(this->CheckIfIsPlasticRegime());

  if (this->CheckIfIsPlasticRegime()) {

    double hardening_modulus = 0.00;
    if (this->GetProperties().Has(DP_K)) {
      hardening_modulus = this->GetProperties()[DP_K];
    }
    const double youngs_modulus = this->GetProperties()[YOUNG_MODULUS];  

    const double trial_yield_function = this->TrialYieldFunction();
    

    const double delta_gamma = trial_yield_function/(youngs_modulus+hardening_modulus);
    current_stress = 1.00 - ((delta_gamma*youngs_modulus)/std::abs(trial_stress));
    current_stress = current_stress * trial_stress;

    this->plastic_strain = this->plastic_strain + (delta_gamma*MathUtils<double>::Sign(trial_stress));
    this->plastic_alpha = this->plastic_alpha + delta_gamma;

    this->test_is_plas = 0.05; //test!!!

  }

  else {

    //MatrixType left_handside_matrix = ZeroMatrix(msLocalSize, msLocalSize);
    // creating LHS
    //left_handside_matrix = this->CreateElementStiffnessMatrix(temp_process_information);
    //Vector nodal_deformation = ZeroVector(msLocalSize);
    //this->GetValuesVector(nodal_deformation, 0);

    //rinternalForces = prod(left_handside_matrix, nodal_deformation);
    this->test_is_plas = 0.0; //test!!!
  }

  this->test_stress_total=current_stress;
  const double area = this->GetProperties()[CROSS_AREA];
  const double truss_axial_force = area*current_stress;
  rinternalForces[0] = -1.00 * truss_axial_force;
  rinternalForces[3] = 1.00 * truss_axial_force;


  //////////////////////////////////////
  ///////// ROTATATE TO GLOBAL CS !!!!!
  ////////////////////////////////////
  KRATOS_CATCH("");
}


void TrussElementLinear3D2N::CalculateOnIntegrationPoints(
    const Variable<double> &rVariable, std::vector<double> &rOutput,
    const ProcessInfo &rCurrentProcessInfo) {
  KRATOS_TRY;

  TrussElement3D2N::CalculateOnIntegrationPoints(rVariable,rOutput,rCurrentProcessInfo);

  const GeometryType::IntegrationPointsArrayType &integration_points =
      GetGeometry().IntegrationPoints();

  if (rOutput.size() != integration_points.size()) {
    rOutput.resize(integration_points.size());
  }

  //test!!!
  if (rVariable == DP_K)rOutput[0] = this->plastic_strain;
  if (rVariable == INFINITY_YIELD_STRESS)rOutput[0] = this->test_is_plas;
  if (rVariable == VON_MISES_STRESS_MIDDLE_SURFACE)rOutput[0] = this->test_stress_total;
  if (rVariable == LAMBDA_MAX)rOutput[0] = this->CalculateLinearStrain();
   

  KRATOS_CATCH("")
  }


void TrussElementLinear3D2N::save(Serializer &rSerializer) const {
  KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, TrussElement3D2N);
}
void TrussElementLinear3D2N::load(Serializer &rSerializer) {
  KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, TrussElement3D2N);
}

} // namespace Kratos.
