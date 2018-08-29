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

Element::Pointer
TrussElementLinear3D2N::Create(IndexType NewId,
                            GeometryType::Pointer pGeom,
                               PropertiesType::Pointer pProperties) const {
  return Kratos::make_shared<TrussElementLinear3D2N>(
      NewId, pGeom, pProperties);
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
  rRightHandSideVector = ZeroVector(msLocalSize);


  BoundedVector<double,msLocalSize> internal_forces =
    this->GetConstitutiveLawTrialResponse(rCurrentProcessInfo,false);

  BoundedMatrix<double, msLocalSize, msLocalSize> transformation_matrix =
      ZeroMatrix(msLocalSize, msLocalSize);
  this->CreateTransformationMatrix(transformation_matrix);

  internal_forces = prod(transformation_matrix, internal_forces);



  rRightHandSideVector -= internal_forces;
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

    array_1d<double, msDimension> temp_internal_stresses = ZeroVector(msDimension);
    ProcessInfo temp_process_information;
    ConstitutiveLaw::Parameters Values(this->GetGeometry(),this->GetProperties(),temp_process_information);

    Vector temp_strain = ZeroVector(1);
    temp_strain[0] = this->CalculateLinearStrain();
    Values.SetStrainVector(temp_strain);
    this->mpConstitutiveLaw->CalculateValue(Values,FORCE,temp_internal_stresses);
    truss_forces[0] = (temp_internal_stresses[0] + prestress) * A;

    rOutput[0] = truss_forces;
  }
}


void TrussElementLinear3D2N::CalculateOnIntegrationPoints(
    const Variable<Vector> &rVariable, std::vector<Vector> &rOutput,
    const ProcessInfo &rCurrentProcessInfo) {
  KRATOS_TRY
  const GeometryType::IntegrationPointsArrayType &integration_points =
      GetGeometry().IntegrationPoints();
  if (rOutput.size() != integration_points.size()) {
    rOutput.resize(integration_points.size());
  }
  if (rVariable == STRAIN) {
    Vector Strain = ZeroVector(msDimension);
    Strain[0] = this->CalculateLinearStrain();
    Strain[1] = 0.00;
    Strain[2] = 0.00;
    rOutput[0] = Strain;
  }
  KRATOS_CATCH("")
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

  Vector current_disp = ZeroVector(msLocalSize);
  this->GetValuesVector(current_disp);
  BoundedMatrix<double, msLocalSize, msLocalSize> transformation_matrix =
      ZeroMatrix(msLocalSize, msLocalSize);
  this->CreateTransformationMatrix(transformation_matrix);

  current_disp = prod(Matrix(trans(transformation_matrix)),current_disp);
  const double length_0 = this->CalculateReferenceLength();
  const double e = (current_disp[3]-current_disp[0])/length_0;

  return e;
  KRATOS_CATCH("");
}


void TrussElementLinear3D2N::UpdateInternalForces(BoundedVector<double,msLocalSize>& rInternalForces)
{
  KRATOS_TRY;

  Vector temp_internal_stresses = ZeroVector(msLocalSize);
  ProcessInfo temp_process_information;
  ConstitutiveLaw::Parameters Values(this->GetGeometry(),this->GetProperties(),temp_process_information);
  Vector temp_strain = ZeroVector(1);
  temp_strain[0] = this->CalculateLinearStrain();
  Values.SetStrainVector(temp_strain);
  this->mpConstitutiveLaw->CalculateValue(Values,NORMAL_STRESS,temp_internal_stresses);

  rInternalForces = temp_internal_stresses*this->GetProperties()[CROSS_AREA];


  BoundedMatrix<double, msLocalSize, msLocalSize> transformation_matrix =
      ZeroMatrix(msLocalSize, msLocalSize);
  this->CreateTransformationMatrix(transformation_matrix);

  rInternalForces = prod(transformation_matrix, rInternalForces);

  KRATOS_CATCH("");
}


void TrussElementLinear3D2N::InitializeNonLinearIteration(ProcessInfo& rCurrentProcessInfo)
{
  KRATOS_TRY;
  this->GetConstitutiveLawTrialResponse(rCurrentProcessInfo,true);
  KRATOS_CATCH("");
}


BoundedVector<double,TrussElementLinear3D2N::msLocalSize>
  TrussElementLinear3D2N::GetConstitutiveLawTrialResponse(
   ProcessInfo& rCurrentProcessInfo,const bool& rSaveInternalVariables)
{
    KRATOS_TRY;
    Vector strain_vector = ZeroVector(this->mpConstitutiveLaw->GetStrainSize());
    Vector stress_vector = ZeroVector(this->mpConstitutiveLaw->GetStrainSize());
    strain_vector[0] = this->CalculateLinearStrain();

    Matrix temp_matrix;
    Vector temp_vector;

    this->mpConstitutiveLaw->CalculateMaterialResponse(strain_vector,
    temp_matrix,stress_vector,temp_matrix,rCurrentProcessInfo,this->GetProperties(),
    this->GetGeometry(),temp_vector,true,true,rSaveInternalVariables);

    BoundedVector<double,msLocalSize> internal_forces = ZeroVector(msLocalSize);
    internal_forces[0] = -1.0 * this->GetProperties()[CROSS_AREA] * stress_vector[0];
    internal_forces[3] = +1.0 * this->GetProperties()[CROSS_AREA] * stress_vector[0];

    return internal_forces;
    KRATOS_CATCH("");
}


void TrussElementLinear3D2N::save(Serializer &rSerializer) const {
  KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, TrussElement3D2N);
  rSerializer.save("mConstitutiveLaw", mpConstitutiveLaw);
}
void TrussElementLinear3D2N::load(Serializer &rSerializer) {
  KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, TrussElement3D2N);
  rSerializer.load("mConstitutiveLaw", mpConstitutiveLaw);
}

} // namespace Kratos.
