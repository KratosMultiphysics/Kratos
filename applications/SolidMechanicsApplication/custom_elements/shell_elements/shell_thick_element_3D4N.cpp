//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:       Massimo Petracca $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:           September 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

#include "shell_thick_element_3D4N.hpp"
#include "custom_utilities/shellq4_corotational_coordinate_transformation.hpp"
#include "solid_mechanics_application.h"

#include "geometries/quadrilateral_3d_4.h"

#include <iomanip>
#include <string>

namespace Kratos {

// =====================================================================================
//
// Utilities
//
// =====================================================================================

namespace Utilities {

template <class TVec> inline void ShapeFunc(double xi, double eta, TVec &N) {
  N(0) = 0.25 * (1.0 - xi) * (1.0 - eta); // node 1
  N(1) = 0.25 * (1.0 + xi) * (1.0 - eta); // node 2
  N(2) = 0.25 * (1.0 + xi) * (1.0 + eta); // node 3
  N(3) = 0.25 * (1.0 - xi) * (1.0 + eta); // node 4
}

template <class TVec>
inline void ShapeFuncSerendipity(double xi, double eta, TVec &N) {
  N(0) = 0.5 * (1.0 - xi * xi) * (1.0 - eta);  // node 5
  N(1) = 0.5 * (1.0 + xi) * (1.0 - eta * eta); // node 6
  N(2) = 0.5 * (1.0 - xi * xi) * (1.0 + eta);  // node 7
  N(3) = 0.5 * (1.0 - xi) * (1.0 - eta * eta); // node 8
}

template <class TMat>
inline void ShapeFunc_NaturalDerivatives(double xi, double eta, TMat &dN) {
  dN(0, 0) = -(1.0 - eta) * 0.25;
  dN(1, 0) = (1.0 - eta) * 0.25;
  dN(2, 0) = (1.0 + eta) * 0.25;
  dN(3, 0) = -(1.0 + eta) * 0.25;

  dN(0, 1) = -(1.0 - xi) * 0.25;
  dN(1, 1) = -(1.0 + xi) * 0.25;
  dN(2, 1) = (1.0 + xi) * 0.25;
  dN(3, 1) = (1.0 - xi) * 0.25;
}

template <class TMat>
inline void ShapeFuncSerendipity_NaturalDerivatives(double xi, double eta,
                                                    TMat &dN) {
  dN(0, 0) = -xi * (1.0 - eta);
  dN(1, 0) = 0.5 * (1.0 - eta * eta);
  dN(2, 0) = -xi * (1.0 - eta);
  dN(3, 0) = -0.5 * (1.0 - eta * eta);

  dN(0, 1) = -0.5 * (1.0 - xi * xi);
  dN(1, 1) = -eta * (1.0 + xi);
  dN(2, 1) = 0.5 * (1.0 - xi * xi);
  dN(3, 1) = -eta * (1.0 + xi);
}
}

// =====================================================================================
//
// Class JacobianOperator
//
// =====================================================================================

ShellThickElement3D4N::JacobianOperator::JacobianOperator()
    : mJac(2, 2, 0.0), mInv(2, 2, 0.0), mXYDeriv(4, 2, 0.0), mDet(0.0) {}

void ShellThickElement3D4N::JacobianOperator::Calculate(
    const ShellQ4_LocalCoordinateSystem &CS, const Matrix &dN) {
  mJac(0, 0) = dN(0, 0) * CS.X1() + dN(1, 0) * CS.X2() + dN(2, 0) * CS.X3() +
               dN(3, 0) * CS.X4();
  mJac(0, 1) = dN(0, 0) * CS.Y1() + dN(1, 0) * CS.Y2() + dN(2, 0) * CS.Y3() +
               dN(3, 0) * CS.Y4();
  mJac(1, 0) = dN(0, 1) * CS.X1() + dN(1, 1) * CS.X2() + dN(2, 1) * CS.X3() +
               dN(3, 1) * CS.X4();
  mJac(1, 1) = dN(0, 1) * CS.Y1() + dN(1, 1) * CS.Y2() + dN(2, 1) * CS.Y3() +
               dN(3, 1) * CS.Y4();

  mDet = mJac(0, 0) * mJac(1, 1) - mJac(1, 0) * mJac(0, 1);
  double mult = 1.0 / mDet;

  mInv(0, 0) = mJac(1, 1) * mult;
  mInv(0, 1) = -mJac(0, 1) * mult;
  mInv(1, 0) = -mJac(1, 0) * mult;
  mInv(1, 1) = mJac(0, 0) * mult;

  noalias(mXYDeriv) = prod(dN, trans(mInv));
}

// =====================================================================================
//
// Class MITC4Params
//
// =====================================================================================

ShellThickElement3D4N::MITC4Params::MITC4Params(
    const ShellQ4_LocalCoordinateSystem &LCS)
    : Transformation(2, 2), ShearStrains(4, 24, 0.0) {

  double x21 = LCS.X2() - LCS.X1();
  double y21 = LCS.Y2() - LCS.Y1();
  double x34 = LCS.X3() - LCS.X4();
  double y34 = LCS.Y3() - LCS.Y4();
  double x41 = LCS.X4() - LCS.X1();
  double y41 = LCS.Y4() - LCS.Y1();
  double x32 = LCS.X3() - LCS.X2();
  double y32 = LCS.Y3() - LCS.Y2();

  Ax = -LCS.X1() + LCS.X2() + LCS.X3() - LCS.X4();
  Bx = LCS.X1() - LCS.X2() + LCS.X3() - LCS.X4();
  Cx = -LCS.X1() - LCS.X2() + LCS.X3() + LCS.X4();
  Ay = -LCS.Y1() + LCS.Y2() + LCS.Y3() - LCS.Y4();
  By = LCS.Y1() - LCS.Y2() + LCS.Y3() - LCS.Y4();
  Cy = -LCS.Y1() - LCS.Y2() + LCS.Y3() + LCS.Y4();

  double Alpha = std::atan(Ay / Ax);
  double Beta = Globals::Pi * 0.5 - std::atan(Cx / Cy);

  Transformation(0, 0) = std::sin(Beta);
  Transformation(0, 1) = -std::sin(Alpha);
  Transformation(1, 0) = -std::cos(Beta);
  Transformation(1, 1) = std::cos(Alpha);

  ShearStrains(0, 2) = -0.5;
  ShearStrains(0, 3) = -y41 * 0.25;
  ShearStrains(0, 4) = x41 * 0.25;

  ShearStrains(0, 20) = 0.5;
  ShearStrains(0, 21) = -y41 * 0.25;
  ShearStrains(0, 22) = x41 * 0.25;

  ShearStrains(1, 2) = -0.5;
  ShearStrains(1, 3) = -y21 * 0.25;
  ShearStrains(1, 4) = x21 * 0.25;

  ShearStrains(1, 8) = 0.5;
  ShearStrains(1, 9) = -y21 * 0.25;
  ShearStrains(1, 10) = x21 * 0.25;

  ShearStrains(2, 8) = -0.5;
  ShearStrains(2, 9) = -y32 * 0.25;
  ShearStrains(2, 10) = x32 * 0.25;

  ShearStrains(2, 14) = 0.5;
  ShearStrains(2, 15) = -y32 * 0.25;
  ShearStrains(2, 16) = x32 * 0.25;

  ShearStrains(3, 14) = 0.5;
  ShearStrains(3, 15) = -y34 * 0.25;
  ShearStrains(3, 16) = x34 * 0.25;

  ShearStrains(3, 20) = -0.5;
  ShearStrains(3, 21) = -y34 * 0.25;
  ShearStrains(3, 22) = x34 * 0.25;
}

// =====================================================================================
//
// Class EASOperatorStorage
//
// =====================================================================================

ShellThickElement3D4N::EASOperatorStorage::EASOperatorStorage()
    : mInitialized(false) {}

void ShellThickElement3D4N::EASOperatorStorage::Initialize(
    const GeometryType &geom) {
  if (!mInitialized) {
    noalias(alpha) = ZeroVector(5);
    noalias(alpha_converged) = ZeroVector(5);

    for (size_t i = 0; i < 4; i++) {
      size_t ii = i * 6;
      const array_1d<double, 3> &initialDispl =
          geom[i].FastGetSolutionStepValue(DISPLACEMENT);
      const array_1d<double, 3> &initialRot =
          geom[i].FastGetSolutionStepValue(ROTATION);

      displ(ii) = initialDispl(0);
      displ(ii + 1) = initialDispl(1);
      displ(ii + 2) = initialDispl(2);
      displ_converged(ii) = initialDispl(0);
      displ_converged(ii + 1) = initialDispl(1);
      displ_converged(ii + 2) = initialDispl(2);

      displ(ii + 3) = initialRot(0);
      displ(ii + 4) = initialRot(1);
      displ(ii + 5) = initialRot(2);
      displ_converged(ii + 3) = initialRot(0);
      displ_converged(ii + 4) = initialRot(1);
      displ_converged(ii + 5) = initialRot(2);
    }

    mInitialized = true;
  }
}

void ShellThickElement3D4N::EASOperatorStorage::InitializeSolutionStep(
    ProcessInfo &CurrentProcessInfo) {
  displ = displ_converged;
  alpha = alpha_converged;
}

void ShellThickElement3D4N::EASOperatorStorage::FinalizeSolutionStep(
    ProcessInfo &CurrentProcessInfo) {
  displ_converged = displ;
  alpha_converged = alpha;
}

void ShellThickElement3D4N::EASOperatorStorage::FinalizeNonLinearIteration(
    const Vector &displacementVector, ProcessInfo &CurrentProcessInfo) {
  Vector incrementalDispl(24);
  noalias(incrementalDispl) = displacementVector - displ;
  noalias(displ) = displacementVector;

  array_1d<double, 5> temp;
  noalias(temp) = prod(L, incrementalDispl);
  noalias(temp) -= residual;
  noalias(alpha) -= prod(Hinv, temp);
}

void ShellThickElement3D4N::EASOperatorStorage::save(
    Serializer &rSerializer) const {
  rSerializer.save("A0", alpha);
  rSerializer.save("A1", alpha_converged);
  rSerializer.save("U0", displ);
  rSerializer.save("U1", displ_converged);
  rSerializer.save("init", mInitialized);
}

void ShellThickElement3D4N::EASOperatorStorage::load(Serializer &rSerializer) {
  rSerializer.load("A0", alpha);
  rSerializer.load("A1", alpha_converged);
  rSerializer.load("U0", displ);
  rSerializer.load("U1", displ_converged);
  rSerializer.load("init", mInitialized);
}

// =====================================================================================
//
// Class EASOperator
//
// =====================================================================================

ShellThickElement3D4N::EASOperator::EASOperator(
    const ShellQ4_LocalCoordinateSystem &LCS, EASOperatorStorage &storage)
    : mF0inv(3, 3), mEnhancedStrains(3), mG(3, 5) {

  // compute the jacobian at the element center

  double xi(0.0);
  double eta(0.0);

  Matrix dN(4, 2);
  Utilities::ShapeFunc_NaturalDerivatives(xi, eta, dN);

  Matrix Jac0(2, 2);
  Jac0(0, 0) = dN(0, 0) * LCS.X1() + dN(1, 0) * LCS.X2() + dN(2, 0) * LCS.X3() +
               dN(3, 0) * LCS.X4();
  Jac0(0, 1) = dN(0, 0) * LCS.Y1() + dN(1, 0) * LCS.Y2() + dN(2, 0) * LCS.Y3() +
               dN(3, 0) * LCS.Y4();
  Jac0(1, 0) = dN(0, 1) * LCS.X1() + dN(1, 1) * LCS.X2() + dN(2, 1) * LCS.X3() +
               dN(3, 1) * LCS.X4();
  Jac0(1, 1) = dN(0, 1) * LCS.Y1() + dN(1, 1) * LCS.Y2() + dN(2, 1) * LCS.Y3() +
               dN(3, 1) * LCS.Y4();

  // save the jacobian determinant at center

  mJ0 = Jac0(0, 0) * Jac0(1, 1) - Jac0(1, 0) * Jac0(0, 1);

  // compute the transformation matrix used in the implementation of the EAS
  // method
  // which operates in the natural coordinate system

  double j11 = Jac0(0, 0);
  double j22 = Jac0(1, 1);
  double j12 = Jac0(0, 1);
  double j21 = Jac0(1, 0);

  Matrix F0(3, 3);
  F0(0, 0) = j11 * j11;
  F0(0, 1) = j21 * j12;
  F0(0, 2) = 2.0 * j11 * j12;
  F0(1, 0) = j12 * j21;
  F0(1, 1) = j22 * j22;
  F0(1, 2) = 2.0 * j21 * j22;
  F0(2, 0) = j11 * j21;
  F0(2, 1) = j12 * j22;
  F0(2, 2) = j11 * j22 + j12 * j21;

  double dummyDet;
  MathUtils<double>::InvertMatrix3(F0, mF0inv, dummyDet);

  // initialize these data to zero because they will
  // be integrated during the gauss loop
  storage.L.clear();
  storage.Hinv.clear();
  storage.residual.clear();
}

void ShellThickElement3D4N::EASOperator::GaussPointComputation_Step1(
    double xi, double eta, const JacobianOperator &jac,
    Vector &generalizedStrains, EASOperatorStorage &storage) {
  // construct the interpolation matrix in natural coordinate system
  Matrix E(3, 5, 0.0);
  E(0, 0) = xi;
  E(1, 1) = eta;
  E(2, 2) = xi;
  E(2, 3) = eta;
  E(0, 4) = xi * eta;
  E(1, 4) = -xi * eta;
  E(2, 4) = xi * xi - eta * eta;

  // construct the interpolation matrix in local coordinate system
  noalias(mG) = mJ0 / jac.Determinant() * prod(mF0inv, E);

  // assuming that the input generalized strains has been already calculated on
  // this
  // integration point, we just need to add the contribution of the enhanced
  // strains
  noalias(mEnhancedStrains) = prod(mG, storage.alpha);

  generalizedStrains(0) += mEnhancedStrains(0); // e.xx
  generalizedStrains(1) += mEnhancedStrains(1); // e.yy
  generalizedStrains(2) += mEnhancedStrains(2); // e.xy
}

void ShellThickElement3D4N::EASOperator::GaussPointComputation_Step2(
    const Matrix &D, const Matrix &B, const Vector &S,
    EASOperatorStorage &storage) {
  Matrix GTC(5, 3);
  noalias(GTC) = prod(trans(mG), project(D, range(0, 3), range(0, 3)));
  noalias(storage.Hinv) += prod(GTC, mG);
  noalias(storage.residual) -= prod(trans(mG), project(S, range(0, 3)));

  // compute L: [G'*Dmm, G'*Dmb, G'*Dms]*[Bm; Bm; Bs]
  int num_stress =
      D.size2(); // it can be 6 for thin shells or 8 for thick  shells
  Matrix GTD(5, num_stress, 0.0);
  noalias(project(GTD, range::all(), range(0, 3))) = GTC;
  noalias(project(GTD, range::all(), range(3, 6))) =
      prod(trans(mG), project(D, range(0, 3), range(3, 6)));
  if (num_stress == 8)
    noalias(project(GTD, range::all(), range(6, 8))) =
        prod(trans(mG), project(D, range(0, 3), range(6, 8)));
  noalias(storage.L) += prod(GTD, B);
}

void ShellThickElement3D4N::EASOperator::ComputeModfiedTangentAndResidual(
    Matrix &rLeftHandSideMatrix, Vector &rRightHandSideVector,
    EASOperatorStorage &storage) {
  // invert H
  Matrix Hcopy(storage.Hinv);
  permutation_matrix<Matrix::size_type> pm(5);
  lu_factorize(Hcopy, pm);
  noalias(storage.Hinv) = IdentityMatrix(5);
  lu_substitute(Hcopy, pm, storage.Hinv);

  // compute L' * H^-1
  Matrix LTHinv(24, 5);
  noalias(LTHinv) = prod(trans(storage.L), storage.Hinv);

  // modify the stiffness matrix and the residual vector
  // for the static condensation of the enhanced strain parameters
  noalias(rRightHandSideVector) -=
      prod(LTHinv, storage.residual); // R_mod = R - L' * H^-1 * residual
  noalias(rLeftHandSideMatrix) -=
      prod(LTHinv, storage.L); // K_mod = K - L' * H^-1 * L
}

// =====================================================================================
//
// Class ShellThickElement3D4N
//
// =====================================================================================

ShellThickElement3D4N::ShellThickElement3D4N(IndexType NewId,
                                             GeometryType::Pointer pGeometry,
                                             bool NLGeom)
    : Element(NewId, pGeometry),
      mpCoordinateTransformation(
          NLGeom ? new ShellQ4_CorotationalCoordinateTransformation(pGeometry)
                 : new ShellQ4_CoordinateTransformation(pGeometry)) {}

ShellThickElement3D4N::ShellThickElement3D4N(
    IndexType NewId, GeometryType::Pointer pGeometry,
    PropertiesType::Pointer pProperties, bool NLGeom)
    : Element(NewId, pGeometry, pProperties),
      mpCoordinateTransformation(
          NLGeom ? new ShellQ4_CorotationalCoordinateTransformation(pGeometry)
                 : new ShellQ4_CoordinateTransformation(pGeometry)) {}

ShellThickElement3D4N::ShellThickElement3D4N(
    IndexType NewId, GeometryType::Pointer pGeometry,
    PropertiesType::Pointer pProperties,
    CoordinateTransformationBasePointerType pCoordinateTransformation)
    : Element(NewId, pGeometry, pProperties),
      mpCoordinateTransformation(pCoordinateTransformation) {}

ShellThickElement3D4N::~ShellThickElement3D4N() {}

Element::Pointer
ShellThickElement3D4N::Create(IndexType NewId, NodesArrayType const &ThisNodes,
                              PropertiesType::Pointer pProperties) const {
  GeometryType::Pointer newGeom(GetGeometry().Create(ThisNodes));
  return Kratos::make_shared<ShellThickElement3D4N>(NewId, newGeom, pProperties,mpCoordinateTransformation->Create(newGeom));
}

void ShellThickElement3D4N::Initialize() {
  KRATOS_TRY

  const GeometryType &geom = GetGeometry();
  const PropertiesType &props = GetProperties();

  if (geom.PointsNumber() != 4)
    KRATOS_THROW_ERROR(
        std::logic_error,
        "ShellThickElement3D4N Element needs a geometry with 4 nodes",
        geom.PointsNumber());

  const GeometryType::IntegrationPointsArrayType &integrationPoints =
      geom.IntegrationPoints(GetIntegrationMethod());
  if (integrationPoints.size() != 4)
    KRATOS_THROW_ERROR(
        std::logic_error,
        "ShellThickElement3D4N Element needs a full integration scheme",
        integrationPoints.size());

  if (mSections.size() != 4) {
    const Matrix &shapeFunctionsValues =
        geom.ShapeFunctionsValues(GetIntegrationMethod());

    ShellCrossSection::Pointer theSection;
    if (props.Has(SHELL_CROSS_SECTION)) {
      theSection = props[SHELL_CROSS_SECTION];
    } else {
      theSection = Kratos::make_shared<ShellCrossSection>();
      theSection->BeginStack();
      theSection->AddPly(props[THICKNESS], 0.0, 5, this->pGetProperties());
      theSection->EndStack();
    }

    mSections.clear();
    for (int i = 0; i < 4; i++) {
      ShellCrossSection::Pointer sectionClone = theSection->Clone();
      sectionClone->SetSectionBehavior(ShellCrossSection::Thick);
      sectionClone->InitializeCrossSection(props, geom,
                                           row(shapeFunctionsValues, i));
      mSections.push_back(sectionClone);
    }
  }

  mpCoordinateTransformation->Initialize();

  this->SetupOrientationAngles();

  mEASStorage.Initialize(geom);

  KRATOS_CATCH("")
}

void ShellThickElement3D4N::ResetConstitutiveLaw() {
  KRATOS_TRY

  const GeometryType &geom = GetGeometry();
  const Matrix &shapeFunctionsValues =
      geom.ShapeFunctionsValues(GetIntegrationMethod());

  const Properties &props = GetProperties();
  for (size_t i = 0; i < mSections.size(); i++)
    mSections[i]->ResetCrossSection(props, geom, row(shapeFunctionsValues, i));

  KRATOS_CATCH("")
}

void ShellThickElement3D4N::EquationIdVector(EquationIdVectorType &rResult,
                                             ProcessInfo &rCurrentProcessInfo) {
  if (rResult.size() != 24)
    rResult.resize(24, false);

  GeometryType &geom = this->GetGeometry();

  for (int i = 0; i < 4; i++) {
    int index = i * 6;
    NodeType &iNode = geom[i];

    rResult[index] = iNode.GetDof(DISPLACEMENT_X).EquationId();
    rResult[index + 1] = iNode.GetDof(DISPLACEMENT_Y).EquationId();
    rResult[index + 2] = iNode.GetDof(DISPLACEMENT_Z).EquationId();

    rResult[index + 3] = iNode.GetDof(ROTATION_X).EquationId();
    rResult[index + 4] = iNode.GetDof(ROTATION_Y).EquationId();
    rResult[index + 5] = iNode.GetDof(ROTATION_Z).EquationId();
  }
}

void ShellThickElement3D4N::GetDofList(DofsVectorType &ElementalDofList,
                                       ProcessInfo &CurrentProcessInfo) {
  ElementalDofList.resize(0);
  ElementalDofList.reserve(24);

  GeometryType &geom = this->GetGeometry();

  for (int i = 0; i < 4; i++) {
    NodeType &iNode = geom[i];

    ElementalDofList.push_back(iNode.pGetDof(DISPLACEMENT_X));
    ElementalDofList.push_back(iNode.pGetDof(DISPLACEMENT_Y));
    ElementalDofList.push_back(iNode.pGetDof(DISPLACEMENT_Z));

    ElementalDofList.push_back(iNode.pGetDof(ROTATION_X));
    ElementalDofList.push_back(iNode.pGetDof(ROTATION_Y));
    ElementalDofList.push_back(iNode.pGetDof(ROTATION_Z));
  }
}

int ShellThickElement3D4N::Check(const ProcessInfo &rCurrentProcessInfo) {
  KRATOS_TRY

  GeometryType &geom = GetGeometry();

  // verify that the variables are correctly initialized
  if (DISPLACEMENT.Key() == 0)
    KRATOS_THROW_ERROR(std::invalid_argument, "DISPLACEMENT has Key zero! "
                                              "(check if the application is "
                                              "correctly registered",
                       "");
  if (ROTATION.Key() == 0)
    KRATOS_THROW_ERROR(std::invalid_argument, "ROTATION has Key zero! (check "
                                              "if the application is correctly "
                                              "registered",
                       "");
  if (VELOCITY.Key() == 0)
    KRATOS_THROW_ERROR(std::invalid_argument, "VELOCITY has Key zero! (check "
                                              "if the application is correctly "
                                              "registered",
                       "");
  if (ACCELERATION.Key() == 0)
    KRATOS_THROW_ERROR(std::invalid_argument, "ACCELERATION has Key zero! "
                                              "(check if the application is "
                                              "correctly registered",
                       "");
  if (DENSITY.Key() == 0)
    KRATOS_THROW_ERROR(std::invalid_argument, "DENSITY has Key zero! (check if "
                                              "the application is correctly "
                                              "registered",
                       "");
  if (SHELL_CROSS_SECTION.Key() == 0)
    KRATOS_THROW_ERROR(std::invalid_argument, "SHELL_CROSS_SECTION has Key "
                                              "zero! (check if the application "
                                              "is correctly registered",
                       "");
  if (THICKNESS.Key() == 0)
    KRATOS_THROW_ERROR(std::invalid_argument, "THICKNESS has Key zero! (check "
                                              "if the application is correctly "
                                              "registered",
                       "");
  if (CONSTITUTIVE_LAW.Key() == 0)
    KRATOS_THROW_ERROR(std::invalid_argument, "CONSTITUTIVE_LAW has Key zero! "
                                              "(check if the application is "
                                              "correctly registered",
                       "");

  // verify that the dofs exist
  for (unsigned int i = 0; i < geom.size(); i++) {
    if (geom[i].SolutionStepsDataHas(DISPLACEMENT) == false)
      KRATOS_THROW_ERROR(std::invalid_argument,
                         "missing variable DISPLACEMENT on node ",
                         geom[i].Id());
    if (geom[i].HasDofFor(DISPLACEMENT_X) == false ||
        geom[i].HasDofFor(DISPLACEMENT_Y) == false ||
        geom[i].HasDofFor(DISPLACEMENT_Z) == false)
      KRATOS_THROW_ERROR(
          std::invalid_argument,
          "missing one of the dofs for the variable DISPLACEMENT on node ",
          GetGeometry()[i].Id());
    if (geom[i].SolutionStepsDataHas(ROTATION) == false)
      KRATOS_THROW_ERROR(std::invalid_argument,
                         "missing variable ROTATION on node ", geom[i].Id());
    if (geom[i].HasDofFor(ROTATION_X) == false ||
        geom[i].HasDofFor(ROTATION_Y) == false ||
        geom[i].HasDofFor(ROTATION_Z) == false)
      KRATOS_THROW_ERROR(
          std::invalid_argument,
          "missing one of the dofs for the variable ROTATION on node ",
          geom[i].Id());

    if (geom[i].GetBufferSize() < 2)
      KRATOS_THROW_ERROR(std::logic_error,
                         "This Element needs at least a buffer size = 2", "");
  }

  // check properties
  if (this->pGetProperties() == NULL)
    KRATOS_THROW_ERROR(std::logic_error, "Properties not provided for element ",
                       this->Id());

  const PropertiesType &props = this->GetProperties();

  if (props.Has(
          SHELL_CROSS_SECTION)) // if the user specified a cross section ...
  {
    const ShellCrossSection::Pointer &section = props[SHELL_CROSS_SECTION];
    if (section == NULL)
      KRATOS_THROW_ERROR(std::logic_error,
                         "SHELL_CROSS_SECTION not provided for element ",
                         this->Id());

    section->Check(props, geom, rCurrentProcessInfo);
  } else // ... allow the automatic creation of a homogeneous section from a
         // material and a thickness
  {
    if (!props.Has(CONSTITUTIVE_LAW))
      KRATOS_THROW_ERROR(std::logic_error,
                         "CONSTITUTIVE_LAW not provided for element ",
                         this->Id());
    const ConstitutiveLaw::Pointer &claw = props[CONSTITUTIVE_LAW];
    if (claw == NULL)
      KRATOS_THROW_ERROR(std::logic_error,
                         "CONSTITUTIVE_LAW not provided for element ",
                         this->Id());

    if (!props.Has(THICKNESS))
      KRATOS_THROW_ERROR(std::logic_error,
                         "THICKNESS not provided for element ", this->Id());
    if (props[THICKNESS] <= 0.0)
      KRATOS_THROW_ERROR(std::logic_error,
                         "wrong THICKNESS value provided for element ",
                         this->Id());

    ShellCrossSection::Pointer dummySection = Kratos::make_shared<ShellCrossSection>();
    dummySection->BeginStack();
    dummySection->AddPly(props[THICKNESS], 0.0, 5, this->pGetProperties());
    dummySection->EndStack();
    dummySection->SetSectionBehavior(ShellCrossSection::Thick);
    dummySection->Check(props, geom, rCurrentProcessInfo);
  }

  return 0;

  KRATOS_CATCH("")
}

void ShellThickElement3D4N::CleanMemory() {}

void ShellThickElement3D4N::GetValuesVector(Vector &values, int Step) {
  if (values.size() != 24)
    values.resize(24, false);

  const GeometryType &geom = GetGeometry();

  for (int i = 0; i < 4; i++) {
    const NodeType &iNode = geom[i];
    const array_1d<double, 3> &disp =
        iNode.FastGetSolutionStepValue(DISPLACEMENT, Step);
    const array_1d<double, 3> &rot =
        iNode.FastGetSolutionStepValue(ROTATION, Step);

    int index = i * 6;
    values[index] = disp[0];
    values[index + 1] = disp[1];
    values[index + 2] = disp[2];

    values[index + 3] = rot[0];
    values[index + 4] = rot[1];
    values[index + 5] = rot[2];
  }
}

void ShellThickElement3D4N::GetFirstDerivativesVector(Vector &values,
                                                      int Step) {
  if (values.size() != 24)
    values.resize(24, false);

  const GeometryType &geom = GetGeometry();

  for (int i = 0; i < 4; i++) {
    const NodeType &iNode = geom[i];
    const array_1d<double, 3> &vel =
        iNode.FastGetSolutionStepValue(VELOCITY, Step);

    int index = i * 6;
    values[index] = vel[0];
    values[index + 1] = vel[1];
    values[index + 2] = vel[2];
    values[index + 3] = 0.0;
    values[index + 4] = 0.0;
    values[index + 5] = 0.0;
  }
}

void ShellThickElement3D4N::GetSecondDerivativesVector(Vector &values,
                                                       int Step) {
  if (values.size() != 24)
    values.resize(24, false);

  const GeometryType &geom = GetGeometry();

  for (int i = 0; i < 4; i++) {
    const NodeType &iNode = geom[i];
    const array_1d<double, 3> &acc =
        iNode.FastGetSolutionStepValue(ACCELERATION, Step);

    int index = i * 6;
    values[index] = acc[0];
    values[index + 1] = acc[1];
    values[index + 2] = acc[2];
    values[index + 3] = 0.0;
    values[index + 4] = 0.0;
    values[index + 5] = 0.0;
  }
}

void ShellThickElement3D4N::InitializeNonLinearIteration(
    ProcessInfo &CurrentProcessInfo) {
  mpCoordinateTransformation->InitializeNonLinearIteration(CurrentProcessInfo);

  const GeometryType &geom = this->GetGeometry();
  const Matrix &shapeFunctionsValues =
      geom.ShapeFunctionsValues(GetIntegrationMethod());
  for (int i = 0; i < 4; i++)
    mSections[i]->InitializeNonLinearIteration(GetProperties(), geom,
                                               row(shapeFunctionsValues, i),
                                               CurrentProcessInfo);
}

void ShellThickElement3D4N::FinalizeNonLinearIteration(
    ProcessInfo &CurrentProcessInfo) {
  mpCoordinateTransformation->FinalizeNonLinearIteration(CurrentProcessInfo);

  ShellQ4_LocalCoordinateSystem LCS(
      mpCoordinateTransformation->CreateLocalCoordinateSystem());
  Vector globalDisplacementVector(24);
  GetValuesVector(globalDisplacementVector);
  Vector localDisplacementVector(
      mpCoordinateTransformation->CalculateLocalDisplacements(
          LCS, globalDisplacementVector));

  mEASStorage.FinalizeNonLinearIteration(localDisplacementVector,
                                         CurrentProcessInfo);

  const GeometryType &geom = this->GetGeometry();
  const Matrix &shapeFunctionsValues =
      geom.ShapeFunctionsValues(GetIntegrationMethod());
  for (int i = 0; i < 4; i++)
    mSections[i]->FinalizeNonLinearIteration(GetProperties(), geom,
                                             row(shapeFunctionsValues, i),
                                             CurrentProcessInfo);
}

void ShellThickElement3D4N::InitializeSolutionStep(
    ProcessInfo &CurrentProcessInfo) {
  const PropertiesType &props = GetProperties();
  const GeometryType &geom = GetGeometry();
  const Matrix &shapeFunctionsValues =
      geom.ShapeFunctionsValues(GetIntegrationMethod());

  for (int i = 0; i < 4; i++)
    mSections[i]->InitializeSolutionStep(
        props, geom, row(shapeFunctionsValues, i), CurrentProcessInfo);

  mpCoordinateTransformation->InitializeSolutionStep(CurrentProcessInfo);

  mEASStorage.InitializeSolutionStep(CurrentProcessInfo);
}

void ShellThickElement3D4N::FinalizeSolutionStep(
    ProcessInfo &CurrentProcessInfo) {
  const PropertiesType &props = GetProperties();
  const GeometryType &geom = GetGeometry();
  const Matrix &shapeFunctionsValues =
      geom.ShapeFunctionsValues(GetIntegrationMethod());

  for (int i = 0; i < 4; i++)
    mSections[i]->FinalizeSolutionStep(
        props, geom, row(shapeFunctionsValues, i), CurrentProcessInfo);

  mpCoordinateTransformation->FinalizeSolutionStep(CurrentProcessInfo);

  mEASStorage.FinalizeSolutionStep(CurrentProcessInfo);
}

void ShellThickElement3D4N::CalculateMassMatrix(
    MatrixType &rMassMatrix, ProcessInfo &rCurrentProcessInfo) {
  if ((rMassMatrix.size1() != 24) || (rMassMatrix.size2() != 24))
    rMassMatrix.resize(24, 24, false);
  noalias(rMassMatrix) = ZeroMatrix(24, 24);

  // Compute the local coordinate system.

  ShellQ4_LocalCoordinateSystem referenceCoordinateSystem(
      mpCoordinateTransformation->CreateReferenceCoordinateSystem());

  // lumped area

  double lump_area = referenceCoordinateSystem.Area() / 4.0;

  // Calculate avarage mass per unit area
  double av_mass_per_unit_area = 0.0;
  for (size_t i = 0; i < 4; i++)
    av_mass_per_unit_area += mSections[i]->CalculateMassPerUnitArea();
  av_mass_per_unit_area /= 4.0;

  // Gauss Loop

  for (size_t i = 0; i < 4; i++) {
    size_t index = i * 6;

    double nodal_mass = av_mass_per_unit_area * lump_area;

    // translational mass
    rMassMatrix(index, index) = nodal_mass;
    rMassMatrix(index + 1, index + 1) = nodal_mass;
    rMassMatrix(index + 2, index + 2) = nodal_mass;

    // rotational mass - neglected for the moment...
  }
}

void ShellThickElement3D4N::CalculateDampingMatrix(
    MatrixType &rDampingMatrix, ProcessInfo &rCurrentProcessInfo) {
  if ((rDampingMatrix.size1() != 24) || (rDampingMatrix.size2() != 24))
    rDampingMatrix.resize(24, 24, false);

  noalias(rDampingMatrix) = ZeroMatrix(24, 24);
}

void ShellThickElement3D4N::CalculateLocalSystem(
    MatrixType &rLeftHandSideMatrix, VectorType &rRightHandSideVector,
    ProcessInfo &rCurrentProcessInfo) {
  CalculateAll(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo,
               true, true);
}

void ShellThickElement3D4N::CalculateRightHandSide(
    VectorType &rRightHandSideVector, ProcessInfo &rCurrentProcessInfo) {
  Matrix dummy;
  CalculateAll(dummy, rRightHandSideVector, rCurrentProcessInfo, true, true);
}

// =====================================================================================
//
// Class ShellThickElement3D4N - Results on Gauss Points
//
// =====================================================================================

void ShellThickElement3D4N::GetValueOnIntegrationPoints(
    const Variable<double> &rVariable, std::vector<double> &rValues,
    const ProcessInfo &rCurrentProcessInfo) {
  SizeType size = GetGeometry().size();
  if (rValues.size() != size)
    rValues.resize(size);

  std::vector<double> temp(size);

  for (SizeType i = 0; i < size; i++)
    mSections[i]->GetValue(rVariable, temp[i]);

  const Matrix &shapeFunctions = GetGeometry().ShapeFunctionsValues();
  Vector N(size);

  for (SizeType i = 0; i < size; i++) {
    noalias(N) = row(shapeFunctions, i);
    double &ival = rValues[i];
    ival = 0.0;
    for (SizeType j = 0; j < size; j++) {
      ival += N(j) * temp[j];
    }
  }
}

void ShellThickElement3D4N::GetValueOnIntegrationPoints(
    const Variable<Vector> &rVariable, std::vector<Vector> &rValues,
    const ProcessInfo &rCurrentProcessInfo) {}

void ShellThickElement3D4N::GetValueOnIntegrationPoints(
    const Variable<Matrix> &rVariable, std::vector<Matrix> &rValues,
    const ProcessInfo &rCurrentProcessInfo) {
  if (TryGetValueOnIntegrationPoints_GeneralizedStrainsOrStresses(
          rVariable, rValues, rCurrentProcessInfo))
    return;
}

void ShellThickElement3D4N::GetValueOnIntegrationPoints(
    const Variable<array_1d<double, 3>> &rVariable,
    std::vector<array_1d<double, 3>> &rValues,
    const ProcessInfo &rCurrentProcessInfo) {
  if (TryGetValueOnIntegrationPoints_MaterialOrientation(rVariable, rValues,
                                                         rCurrentProcessInfo))
    return;
}

void ShellThickElement3D4N::GetValueOnIntegrationPoints(
    const Variable<array_1d<double, 6>> &rVariable,
    std::vector<array_1d<double, 6>> &rValues,
    const ProcessInfo &rCurrentProcessInfo) {}

void ShellThickElement3D4N::SetCrossSectionsOnIntegrationPoints(
    std::vector<ShellCrossSection::Pointer> &crossSections) {
  KRATOS_TRY
  if (crossSections.size() != 4)
    KRATOS_THROW_ERROR(std::logic_error,
                       "Cannot set a number of cross section different from 4",
                       "");
  mSections.clear();
  for (SizeType i = 0; i < crossSections.size(); i++)
    mSections.push_back(crossSections[i]);
  this->SetupOrientationAngles();
  KRATOS_CATCH("")
}

// =====================================================================================
//
// Class ShellThickElement3D4N - Private methods
//
// =====================================================================================

void ShellThickElement3D4N::DecimalCorrection(Vector &a) {
  double norm = norm_2(a);
  double tolerance = std::max(norm * 1.0E-12, 1.0E-12);
  for (SizeType i = 0; i < a.size(); i++)
    if (std::abs(a(i)) < tolerance)
      a(i) = 0.0;
}

void ShellThickElement3D4N::SetupOrientationAngles() {
  ShellQ4_LocalCoordinateSystem lcs(
      mpCoordinateTransformation->CreateReferenceCoordinateSystem());

  Vector3Type normal;
  noalias(normal) = lcs.Vz();

  Vector3Type dZ;
  dZ(0) = 0.0;
  dZ(1) = 0.0;
  dZ(2) = 1.0; // for the moment let's take this. But the user can specify its
               // own triad! TODO

  Vector3Type dirX;
  MathUtils<double>::CrossProduct(dirX, dZ, normal);

  // try to normalize the x vector. if it is near zero it means that we need
  // to choose a default one.
  double dirX_norm = dirX(0) * dirX(0) + dirX(1) * dirX(1) + dirX(2) * dirX(2);
  if (dirX_norm < 1.0E-12) {
    dirX(0) = 1.0;
    dirX(1) = 0.0;
    dirX(2) = 0.0;
  } else if (dirX_norm != 1.0) {
    dirX_norm = std::sqrt(dirX_norm);
    dirX /= dirX_norm;
  }

  Vector3Type elem_dirX = lcs.Vx();

  // now calculate the angle between the element x direction and the material x
  // direction.
  Vector3Type &a = elem_dirX;
  Vector3Type &b = dirX;
  double a_dot_b = a(0) * b(0) + a(1) * b(1) + a(2) * b(2);
  if (a_dot_b < -1.0)
    a_dot_b = -1.0;
  if (a_dot_b > 1.0)
    a_dot_b = 1.0;
  double angle = std::acos(a_dot_b);

  // if they are not counter-clock-wise, let's change the sign of the angle
  if (angle != 0.0) {
    const MatrixType &R = lcs.Orientation();
    if (dirX(0) * R(1, 0) + dirX(1) * R(1, 1) + dirX(2) * R(1, 2) < 0.0)
      angle = -angle;
  }

  for (CrossSectionContainerType::iterator it = mSections.begin();
       it != mSections.end(); ++it)
    (*it)->SetOrientationAngle(angle);
}

void ShellThickElement3D4N::CalculateBMatrix(double xi, double eta,
                                             const JacobianOperator &Jac,
                                             const MITC4Params &mitc_params,
                                             const Vector &N, Matrix &B,
                                             Vector &Bdrill) {
  // some data

  const Matrix &dNxy = Jac.XYDerivatives();

  // membrane
  // ************************************************************************************************

  B(0, 0) = dNxy(0, 0);
  B(0, 6) = dNxy(1, 0);
  B(0, 12) = dNxy(2, 0);
  B(0, 18) = dNxy(3, 0);
  B(1, 1) = dNxy(0, 1);
  B(1, 7) = dNxy(1, 1);
  B(1, 13) = dNxy(2, 1);
  B(1, 19) = dNxy(3, 1);
  B(2, 0) = dNxy(0, 1);
  B(2, 6) = dNxy(1, 1);
  B(2, 12) = dNxy(2, 1);
  B(2, 18) = dNxy(3, 1);
  B(2, 1) = dNxy(0, 0);
  B(2, 7) = dNxy(1, 0);
  B(2, 13) = dNxy(2, 0);
  B(2, 19) = dNxy(3, 0);

  // bending
  // *************************************************************************************************

  B(3, 4) = dNxy(0, 0);
  B(3, 10) = dNxy(1, 0);
  B(3, 16) = dNxy(2, 0);
  B(3, 22) = dNxy(3, 0);
  B(4, 3) = -dNxy(0, 1);
  B(4, 9) = -dNxy(1, 1);
  B(4, 15) = -dNxy(2, 1);
  B(4, 21) = -dNxy(3, 1);
  B(5, 3) = -dNxy(0, 0);
  B(5, 9) = -dNxy(1, 0);
  B(5, 15) = -dNxy(2, 0);
  B(5, 21) = -dNxy(3, 0);
  B(5, 4) = dNxy(0, 1);
  B(5, 10) = dNxy(1, 1);
  B(5, 16) = dNxy(2, 1);
  B(5, 22) = dNxy(3, 1);

  // drilling
  // ************************************************************************************************

  Bdrill(0) = -0.5 * dNxy(0, 1);
  Bdrill(1) = 0.5 * dNxy(0, 0);
  Bdrill(5) = -N(0);
  Bdrill(6) = -0.5 * dNxy(1, 1);
  Bdrill(7) = 0.5 * dNxy(1, 0);
  Bdrill(11) = -N(1);
  Bdrill(12) = -0.5 * dNxy(2, 1);
  Bdrill(13) = 0.5 * dNxy(2, 0);
  Bdrill(17) = -N(2);
  Bdrill(18) = -0.5 * dNxy(3, 1);
  Bdrill(19) = 0.5 * dNxy(3, 0);
  Bdrill(23) = -N(3);

  // shear
  // ***************************************************************************************************

  // MITC modified shape functions
  Matrix MITCShapeFunctions(2, 4, 0.0);
  MITCShapeFunctions(1, 0) = 1.0 - xi;
  MITCShapeFunctions(0, 1) = 1.0 - eta;
  MITCShapeFunctions(1, 2) = 1.0 + xi;
  MITCShapeFunctions(0, 3) = 1.0 + eta;

  // strain displacement matrix in natural coordinate system.
  // interpolate the shear strains given in MITC4Params
  // using the modified shape function
  Matrix BN = prod(MITCShapeFunctions, mitc_params.ShearStrains);

  // modify the shear strain intensity in the tying points
  // to match the values that would be obtained using standard
  // interpolations
  double Temp1, Temp2, Temp3;
  Temp1 = mitc_params.Cx + xi * mitc_params.Bx;
  Temp3 = mitc_params.Cy + xi * mitc_params.By;
  Temp1 = Temp1 * Temp1 + Temp3 * Temp3;
  Temp1 = std::sqrt(Temp1) / (8.0 * Jac.Determinant());
  Temp2 = mitc_params.Ax + eta * mitc_params.Bx;
  Temp3 = mitc_params.Ay + eta * mitc_params.By;
  Temp2 = Temp2 * Temp2 + Temp3 * Temp3;
  Temp2 = std::sqrt(Temp2) / (8.0 * Jac.Determinant());

  row(BN, 0) *= Temp1;
  row(BN, 1) *= Temp2;

  // transform the strain-displacement matrix from natural
  // to local coordinate system taking into account the element distorsion

  project(B, slice(7, -1, 2), slice::all()) =
      prod(mitc_params.Transformation, BN);

  // Explanation of the 'slice':
  // The MITC4Params class and the first part of this method were coded
  // assuming the following notations for strain vector : [e.xx, e.yy, e.zz, 2
  // e.xy, 2 e.XZ, 2 e.YZ].
  // Since in Kratos the strain vector is assumed to be : [e.xx, e.yy, e.zz, 2
  // e.xy, 2 e.YZ, 2 e.XZ],
  // I had to "virtually" swap the 2 rows of the matrix B before the assignment
  // from the product on the Right-Hand-Side.
  // In this way the matrix B is consistent with the notation used in Kratos
  // without recoding all the MITC stuff.
}

void ShellThickElement3D4N::CalculateAll(MatrixType &rLeftHandSideMatrix,
                                         VectorType &rRightHandSideVector,
                                         ProcessInfo &rCurrentProcessInfo,
                                         const bool LHSrequired,
                                         const bool RHSrequired) {
  // Resize the Left Hand Side if necessary,
  // and initialize it to Zero

  if ((rLeftHandSideMatrix.size1() != 24) ||
      (rLeftHandSideMatrix.size2() != 24))
    rLeftHandSideMatrix.resize(24, 24, false);
  noalias(rLeftHandSideMatrix) = ZeroMatrix(24, 24);

  // Resize the Right Hand Side if necessary,
  // and initialize it to Zero

  if (rRightHandSideVector.size() != 24)
    rRightHandSideVector.resize(24, false);
  noalias(rRightHandSideVector) = ZeroVector(24);

  // Get some references.

  PropertiesType &props = GetProperties();
  GeometryType &geom = GetGeometry();
  const Matrix &shapeFunctions = geom.ShapeFunctionsValues();
  Vector iN(shapeFunctions.size2());

  // Compute the local coordinate system.

  ShellQ4_LocalCoordinateSystem localCoordinateSystem(
      mpCoordinateTransformation->CreateLocalCoordinateSystem());

  ShellQ4_LocalCoordinateSystem referenceCoordinateSystem(
      mpCoordinateTransformation->CreateReferenceCoordinateSystem());

  // Prepare all the parameters needed for the MITC formulation.
  // This is to be done here outside the Gauss Loop.

  MITC4Params shearParameters(referenceCoordinateSystem);

  // Instantiate the Jacobian Operator.
  // This will store:
  // the jacobian matrix, its inverse, its determinant
  // and the derivatives of the shape functions in the local
  // coordinate system

  JacobianOperator jacOp;
  array_1d<double, 4> dArea;

  // Instantiate all strain-displacement matrices.

  Matrix B(8, 24, 0.0);
  Vector Bdrilling(24, 0.0);

  // Instantiate all section tangent matrices.

  Matrix D(8, 8, 0.0);
  double Ddrilling(0.0);

  // Instantiate strain and stress-resultant vectors

  Vector generalizedStrains(8);
  Vector generalizedStresses(8);

  // Get the current displacements in global coordinate system

  Vector globalDisplacements(24);
  GetValuesVector(globalDisplacements, 0);

  // Get the current displacements in local coordinate system

  Vector localDisplacements(
      mpCoordinateTransformation->CalculateLocalDisplacements(
          localCoordinateSystem, globalDisplacements));

  // Instantiate the EAS Operator.
  // This will apply the Enhanced Assumed Strain Method for the calculation
  // of the membrane contribution.

  EASOperator EASOp(referenceCoordinateSystem, mEASStorage);

  // auxiliary data

  Matrix BTD(24, 8); // auxiliary matrix to store the product B'*D

  // Initialize parameters for the cross section calculation

  ShellCrossSection::Parameters parameters(geom, props, rCurrentProcessInfo);
  parameters.SetGeneralizedStrainVector(generalizedStrains);
  parameters.SetGeneralizedStressVector(generalizedStresses);
  parameters.SetConstitutiveMatrix(D);
  Flags &options = parameters.GetOptions();
  options.Set(ConstitutiveLaw::COMPUTE_STRESS, RHSrequired);
  options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, LHSrequired);

  // Gauss Loop.
  for (int i = 0; i < 4; i++) {

    // get a reference of the current integration point and shape functions

    const GeometryType::IntegrationPointType &ip = geom.IntegrationPoints()[i];

    noalias(iN) = row(shapeFunctions, i);

    // Compute Jacobian, Inverse of Jacobian, Determinant of Jacobian
    // and Shape functions derivatives in the local coordinate system

    jacOp.Calculate(referenceCoordinateSystem,
                    geom.ShapeFunctionLocalGradient(i));

    // compute the 'area' of the current integration point

    double dA = ip.Weight() * jacOp.Determinant();
    dArea[i] = dA;

    // Compute all strain-displacement matrices

    CalculateBMatrix(ip.X(), ip.Y(), jacOp, shearParameters, iN, B, Bdrilling);

    // Calculate strain vectors in local coordinate system

    noalias(generalizedStrains) = prod(B, localDisplacements);
    double drillingStrain = inner_prod(Bdrilling, localDisplacements);

    // Apply the EAS method to modify the membrane part of the strains computed
    // above.

    EASOp.GaussPointComputation_Step1(ip.X(), ip.Y(), jacOp, generalizedStrains,
                                      mEASStorage);

    // Calculate the response of the Cross Section

    ShellCrossSection::Pointer &section = mSections[i];

    parameters.SetShapeFunctionsValues(iN);
    parameters.SetShapeFunctionsDerivatives(jacOp.XYDerivatives());
    section->CalculateSectionResponse(parameters,
                                      ConstitutiveLaw::StressMeasure_Kirchhoff);
    Ddrilling = section->GetDrillingStiffness();

    // multiply the section tangent matrices and stress resultants by 'dA'

    D *= dA;
    Ddrilling *= dA;
    generalizedStresses *= dA;
    double drillingStress =
        Ddrilling * drillingStrain; // already multiplied by 'dA'

    // Add all contributions to the Stiffness Matrix

    noalias(BTD) = prod(trans(B), D);
    noalias(rLeftHandSideMatrix) += prod(BTD, B);
    noalias(rLeftHandSideMatrix) +=
        outer_prod(Ddrilling * Bdrilling, Bdrilling);

    // Add all contributions to the residual vector

    noalias(rRightHandSideVector) -= prod(trans(B), generalizedStresses);
    noalias(rRightHandSideVector) -= Bdrilling * drillingStress;

    // Continue the calculation of the EAS method now that the contitutive
    // response
    // has been computed

    EASOp.GaussPointComputation_Step2(D, B, generalizedStresses, mEASStorage);
  }

  // Now that the gauss integration is over, let the EAS operator modify
  // the local stiffness matrix and residual vector.
  // It will perform a static condensation to remove the enhanced strain
  // parameters
  // at the element level

  EASOp.ComputeModfiedTangentAndResidual(rLeftHandSideMatrix,
                                         rRightHandSideVector, mEASStorage);

  // Let the CoordinateTransformation finalize the calculation.
  // This will handle the transformation of the local matrices/vectors to
  // the global coordinate system.

  mpCoordinateTransformation->FinalizeCalculations(
      localCoordinateSystem, globalDisplacements, localDisplacements,
      rLeftHandSideMatrix, rRightHandSideVector, RHSrequired, LHSrequired);

  // Add body forces contributions. This doesn't depend on the coordinate system
  AddBodyForces(dArea, rRightHandSideVector);
}

void ShellThickElement3D4N::AddBodyForces(const array_1d<double, 4> &dA,
                                          VectorType &rRightHandSideVector) {
  const GeometryType &geom = GetGeometry();

  // Get shape functions
  const Matrix &N = geom.ShapeFunctionsValues();

  // auxiliary
  array_1d<double, 3> bf;

  // gauss loop to integrate the external force vector
  for (unsigned int igauss = 0; igauss < 4; igauss++) {
    // get mass per unit area
    double mass_per_unit_area = mSections[igauss]->CalculateMassPerUnitArea();

    // interpolate nodal volume accelerations to this gauss point
    // and obtain the body force vector
    bf.clear();
    for (unsigned int inode = 0; inode < 4; inode++) {
      if (geom[inode].SolutionStepsDataHas(VOLUME_ACCELERATION)) // temporary,
                                                                 // will be
                                                                 // checked once
                                                                 // at the
                                                                 // beginning
                                                                 // only
        bf += N(igauss, inode) *
              geom[inode].FastGetSolutionStepValue(VOLUME_ACCELERATION);
    }
    bf *= (mass_per_unit_area * dA[igauss]);

    // add it to the RHS vector
    for (unsigned int inode = 0; inode < 4; inode++) {
      unsigned int index = inode * 6;
      double iN = N(igauss, inode);
      rRightHandSideVector[index + 0] += iN * bf[0];
      rRightHandSideVector[index + 1] += iN * bf[1];
      rRightHandSideVector[index + 2] += iN * bf[2];
    }
  }
}

bool ShellThickElement3D4N::TryGetValueOnIntegrationPoints_MaterialOrientation(
    const Variable<array_1d<double, 3>> &rVariable,
    std::vector<array_1d<double, 3>> &rValues,
    const ProcessInfo &rCurrentProcessInfo) {
  // Check the required output

  int ijob = 0;
  if (rVariable == MATERIAL_ORIENTATION_DX)
    ijob = 1;
  else if (rVariable == MATERIAL_ORIENTATION_DY)
    ijob = 2;
  else if (rVariable == MATERIAL_ORIENTATION_DZ)
    ijob = 3;

  // quick return

  if (ijob == 0)
    return false;

  // resize output

  size_t size = 4;
  if (rValues.size() != size)
    rValues.resize(size);

  // Get some references.

  // GeometryType & geom = GetGeometry();

  // Compute the local coordinate system.

  ShellQ4_LocalCoordinateSystem localCoordinateSystem(
      mpCoordinateTransformation->CreateLocalCoordinateSystem());

  Vector3Type eZ = localCoordinateSystem.Vz();

  // Gauss Loop

  if (ijob == 1) {
    Vector3Type eX = localCoordinateSystem.Vx();
    for (int i = 0; i < 4; i++) {
      QuaternionType q = QuaternionType::FromAxisAngle(
          eZ(0), eZ(1), eZ(2), mSections[i]->GetOrientationAngle());
      q.RotateVector3(eX, rValues[i]);
    }
  } else if (ijob == 2) {
    Vector3Type eY = localCoordinateSystem.Vy();
    for (int i = 0; i < 4; i++) {
      QuaternionType q = QuaternionType::FromAxisAngle(
          eZ(0), eZ(1), eZ(2), mSections[i]->GetOrientationAngle());
      q.RotateVector3(eY, rValues[i]);
    }
  } else if (ijob == 3) {
    for (int i = 0; i < 4; i++) {
      noalias(rValues[i]) = eZ;
    }
  }

  return true;
}

bool ShellThickElement3D4N::
    TryGetValueOnIntegrationPoints_GeneralizedStrainsOrStresses(
        const Variable<Matrix> &rVariable, std::vector<Matrix> &rValues,
        const ProcessInfo &rCurrentProcessInfo) {
  // Check the required output

  int ijob = 0;
  bool bGlobal = false;
  if (rVariable == SHELL_STRAIN) {
    ijob = 1;
  } else if (rVariable == SHELL_STRAIN_GLOBAL) {
    ijob = 1;
    bGlobal = true;
  } else if (rVariable == SHELL_CURVATURE) {
    ijob = 2;
  } else if (rVariable == SHELL_CURVATURE_GLOBAL) {
    ijob = 2;
    bGlobal = true;
  } else if (rVariable == SHELL_FORCE) {
    ijob = 3;
  } else if (rVariable == SHELL_FORCE_GLOBAL) {
    ijob = 3;
    bGlobal = true;
  } else if (rVariable == SHELL_MOMENT) {
    ijob = 4;
  } else if (rVariable == SHELL_MOMENT_GLOBAL) {
    ijob = 4;
    bGlobal = true;
  }

  // quick return

  if (ijob == 0)
    return false;

  // resize output

  size_t size = 4;
  if (rValues.size() != size)
    rValues.resize(size);

  // Get some references.

  PropertiesType &props = GetProperties();
  GeometryType &geom = GetGeometry();
  const Matrix &shapeFunctions = geom.ShapeFunctionsValues();
  Vector iN(shapeFunctions.size2());

  // Compute the local coordinate system.

  ShellQ4_LocalCoordinateSystem localCoordinateSystem(
      mpCoordinateTransformation->CreateLocalCoordinateSystem());

  ShellQ4_LocalCoordinateSystem referenceCoordinateSystem(
      mpCoordinateTransformation->CreateReferenceCoordinateSystem());

  // Prepare all the parameters needed for the MITC formulation.
  // This is to be done here outside the Gauss Loop.

  MITC4Params shearParameters(referenceCoordinateSystem);

  // Instantiate the Jacobian Operator.
  // This will store:
  // the jacobian matrix, its inverse, its determinant
  // and the derivatives of the shape functions in the local
  // coordinate system

  JacobianOperator jacOp;

  // Instantiate all strain-displacement matrices.

  Matrix B(8, 24, 0.0);
  Vector Bdrilling(24, 0.0);

  // Instantiate all section tangent matrices.

  Matrix D(8, 8, 0.0);

  // Instantiate strain and stress-resultant vectors

  Vector generalizedStrains(8);
  Vector generalizedStresses(8);

  // Get the current displacements in global coordinate system

  Vector globalDisplacements(24);
  GetValuesVector(globalDisplacements, 0);

  // Get the current displacements in local coordinate system

  Vector localDisplacements(
      mpCoordinateTransformation->CalculateLocalDisplacements(
          localCoordinateSystem, globalDisplacements));

  // Instantiate the EAS Operator.
  // This will apply the Enhanced Assumed Strain Method for the calculation
  // of the membrane contribution.

  EASOperator EASOp(referenceCoordinateSystem, mEASStorage);

  // Just to store the rotation matrix for visualization purposes
  Matrix R(8, 8);
  Matrix aux33(3, 3);

  // Initialize parameters for the cross section calculation

  ShellCrossSection::Parameters parameters(geom, props, rCurrentProcessInfo);
  parameters.SetGeneralizedStrainVector(generalizedStrains);
  parameters.SetGeneralizedStressVector(generalizedStresses);
  parameters.SetConstitutiveMatrix(D);
  Flags &options = parameters.GetOptions();
  options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
  options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);

  // Gauss Loop

  for (unsigned int i = 0; i < size; i++) {

    // get a reference of the current integration point and shape functions

    const GeometryType::IntegrationPointType &ip = geom.IntegrationPoints()[i];

    noalias(iN) = row(shapeFunctions, i);

    // Compute Jacobian, Inverse of Jacobian, Determinant of Jacobian
    // and Shape functions derivatives in the local coordinate system

    jacOp.Calculate(referenceCoordinateSystem,
                    geom.ShapeFunctionLocalGradient(i));

    // Compute all strain-displacement matrices

    CalculateBMatrix(ip.X(), ip.Y(), jacOp, shearParameters, iN, B, Bdrilling);

    // Calculate strain vectors in local coordinate system

    noalias(generalizedStrains) = prod(B, localDisplacements);

    // Apply the EAS method to modify the membrane part of the strains computed
    // above.

    EASOp.GaussPointComputation_Step1(ip.X(), ip.Y(), jacOp, generalizedStrains,
                                      mEASStorage);

    // Calculate the response of the Cross Section

    ShellCrossSection::Pointer &section = mSections[i];
    if (ijob > 2) {
      parameters.SetShapeFunctionsValues(iN);
      parameters.SetShapeFunctionsDerivatives(jacOp.XYDerivatives());
      section->CalculateSectionResponse(
          parameters, ConstitutiveLaw::StressMeasure_Kirchhoff);
    }

    // save the results

    DecimalCorrection(generalizedStrains);
    DecimalCorrection(generalizedStresses);

    // now the results are in the element coordinate system
    // if necessary, rotate the results in the section (local) coordinate system
    if (section->GetOrientationAngle() != 0.0 && !bGlobal) {
      if (ijob > 2) {
        section->GetRotationMatrixForGeneralizedStresses(
            -(section->GetOrientationAngle()), R);
        generalizedStresses = prod(R, generalizedStresses);
      } else {
        section->GetRotationMatrixForGeneralizedStrains(
            -(section->GetOrientationAngle()), R);
        generalizedStrains = prod(R, generalizedStrains);
      }
    }

    Matrix &iValue = rValues[i];
    if (iValue.size1() != 3 || iValue.size2() != 3)
      iValue.resize(3, 3, false);

    if (ijob == 1) // strains
    {
      iValue(0, 0) = generalizedStrains(0);
      iValue(1, 1) = generalizedStrains(1);
      iValue(2, 2) = 0.0;
      iValue(0, 1) = iValue(1, 0) = 0.5 * generalizedStrains(2);
      iValue(0, 2) = iValue(2, 0) = 0.5 * generalizedStrains(7);
      iValue(1, 2) = iValue(2, 1) = 0.5 * generalizedStrains(6);
    } else if (ijob == 2) // curvatures
    {
      iValue(0, 0) = generalizedStrains(3);
      iValue(1, 1) = generalizedStrains(4);
      iValue(2, 2) = 0.0;
      iValue(0, 1) = iValue(1, 0) = 0.5 * generalizedStrains(5);
      iValue(0, 2) = iValue(2, 0) = 0.0;
      iValue(1, 2) = iValue(2, 1) = 0.0;
    } else if (ijob == 3) // forces
    {
      iValue(0, 0) = generalizedStresses(0);
      iValue(1, 1) = generalizedStresses(1);
      iValue(2, 2) = 0.0;
      iValue(0, 1) = iValue(1, 0) = generalizedStresses(2);
      iValue(0, 2) = iValue(2, 0) = generalizedStresses(7);
      iValue(1, 2) = iValue(2, 1) = generalizedStresses(6);
    } else if (ijob == 4) // moments
    {
      iValue(0, 0) = generalizedStresses(3);
      iValue(1, 1) = generalizedStresses(4);
      iValue(2, 2) = 0.0;
      iValue(0, 1) = iValue(1, 0) = generalizedStresses(5);
      iValue(0, 2) = iValue(2, 0) = 0.0;
      iValue(1, 2) = iValue(2, 1) = 0.0;
    }

    // if requested, rotate the results in the global coordinate system
    if (bGlobal) {
      const Matrix &RG = localCoordinateSystem.Orientation();
      noalias(aux33) = prod(trans(RG), iValue);
      noalias(iValue) = prod(aux33, RG);
    }

  } // Gauss Loop

  return true;
}

// =====================================================================================
//
// Class ShellThickElement3D4N - Serialization
//
// =====================================================================================

void ShellThickElement3D4N::save(Serializer &rSerializer) const {
  KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);
  rSerializer.save("CTr", mpCoordinateTransformation);
  rSerializer.save("Sec", mSections);
  rSerializer.save("EAS", mEASStorage);
}

void ShellThickElement3D4N::load(Serializer &rSerializer) {
  KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);
  rSerializer.load("CTr", mpCoordinateTransformation);
  rSerializer.load("Sec", mSections);
  rSerializer.load("EAS", mEASStorage);
}
}
