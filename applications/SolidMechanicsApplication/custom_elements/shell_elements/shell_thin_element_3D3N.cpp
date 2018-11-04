//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:       Massimo Petracca $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:           September 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

#include "shell_thin_element_3D3N.hpp"
#include "custom_utilities/shellt3_corotational_coordinate_transformation.hpp"
#include "solid_mechanics_application.h"

#include "geometries/triangle_3d_3.h"

#include <string>
#include <iomanip>

#define OPT_NUM_NODES 3
#define OPT_STRAIN_SIZE 6
#define OPT_NUM_DOFS 18

//----------------------------------------
// preprocessors for the integration
// method used by this element.

//#define OPT_1_POINT_INTEGRATION

#ifdef OPT_1_POINT_INTEGRATION
#define OPT_INTEGRATION_METHOD Kratos::GeometryData::GI_GAUSS_1
#define OPT_NUM_GP 1
#else
#define OPT_INTEGRATION_METHOD Kratos::GeometryData::GI_GAUSS_2
#define OPT_NUM_GP 3
#endif // OPT_1_POINT_INTEGRATION

//----------------------------------------
// preprocessors to handle the output
// in case of 3 integration points

//#define OPT_USES_INTERIOR_GAUSS_POINTS

#ifdef OPT_1_POINT_INTEGRATION
#define OPT_INTERPOLATE_RESULTS_TO_STANDARD_GAUSS_POINTS(X)
#else
#ifdef OPT_USES_INTERIOR_GAUSS_POINTS
#define OPT_INTERPOLATE_RESULTS_TO_STANDARD_GAUSS_POINTS(X)
#else
#define OPT_INTERPOLATE_RESULTS_TO_STANDARD_GAUSS_POINTS(X) Utilities::InterpToStandardGaussPoints(X)
#endif // OPT_USES_INTERIOR_GAUSS_POINTS
#endif // OPT_1_POINT_INTEGRATION

//#define OPT_AVARAGE_RESULTS

namespace Kratos
{

namespace Utilities
{
inline void InterpToStandardGaussPoints(double& v1, double& v2, double& v3)
{
  double vg1 = v1;
  double vg2 = v2;
  double vg3 = v3;
#ifdef OPT_AVARAGE_RESULTS
  v1 = (vg1+vg2+vg3)/3.0;
  v2 = (vg1+vg2+vg3)/3.0;
  v3 = (vg1+vg2+vg3)/3.0;
#else
  v1 = (2.0*vg1)/3.0 - vg2/3.0       + (2.0*vg3)/3.0;
  v2 = (2.0*vg1)/3.0 + (2.0*vg2)/3.0 - vg3/3.0;
  v3 = (2.0*vg2)/3.0 - vg1/3.0       + (2.0*vg3)/3.0;
#endif // OPT_AVARAGE_RESULTS
}

inline void InterpToStandardGaussPoints(std::vector< double >& v)
{
  if(v.size() != 3) return;
  InterpToStandardGaussPoints(v[0], v[1], v[2]);
}

inline void InterpToStandardGaussPoints(std::vector< array_1d<double,3> >& v)
{
  if(v.size() != 3) return;
  for(size_t i = 0; i < 3; i++)
    InterpToStandardGaussPoints(v[0][i], v[1][i], v[2][i]);
}

inline void InterpToStandardGaussPoints(std::vector< array_1d<double,6> >& v)
{
  if(v.size() != 3) return;
  for(size_t i = 0; i < 6; i++)
    InterpToStandardGaussPoints(v[0][i], v[1][i], v[2][i]);
}

inline void InterpToStandardGaussPoints(std::vector< Vector >& v)
{
  if(v.size() != 3) return;
  size_t ncomp = v[0].size();
  for(int i = 1; i < 3; i++)
    if(v[i].size() != ncomp)
      return;
  for(size_t i = 0; i < ncomp; i++)
    InterpToStandardGaussPoints(v[0][i], v[1][i], v[2][i]);
}

inline void InterpToStandardGaussPoints(std::vector< Matrix >& v)
{
  if(v.size() != 3) return;
  size_t nrows = v[0].size1();
  size_t ncols = v[0].size2();
  for(int i = 1; i < 3; i++)
    if(v[i].size1() != nrows || v[i].size2() != ncols)
      return;
  for(size_t i = 0; i < nrows; i++)
    for(size_t j = 0; j < ncols; j++)
      InterpToStandardGaussPoints(v[0](i,j), v[1](i,j), v[2](i,j));
}

}

// =====================================================================================
//
// CalculationData
//
// =====================================================================================

ShellThinElement3D3N::CalculationData::CalculationData(const CoordinateTransformationBasePointerType& pCoordinateTransformation,
                                                       const ProcessInfo& rCurrentProcessInfo)
    : LCS0( pCoordinateTransformation->CreateReferenceCoordinateSystem() )
    , LCS( pCoordinateTransformation->CreateLocalCoordinateSystem() )
    , CurrentProcessInfo(rCurrentProcessInfo)

{
}

// =====================================================================================
//
// Class ShellThinElement3D3N
//
// =====================================================================================

ShellThinElement3D3N::ShellThinElement3D3N(IndexType NewId,
                                           GeometryType::Pointer pGeometry,
                                           bool NLGeom)
    : Element(NewId, pGeometry)
    , mpCoordinateTransformation( NLGeom ?
                                  new ShellT3_CorotationalCoordinateTransformation(pGeometry) :
                                  new ShellT3_CoordinateTransformation(pGeometry))
{
  mThisIntegrationMethod = OPT_INTEGRATION_METHOD;
}

ShellThinElement3D3N::ShellThinElement3D3N(IndexType NewId,
                                           GeometryType::Pointer pGeometry,
                                           PropertiesType::Pointer pProperties,
                                           bool NLGeom)
    : Element(NewId, pGeometry, pProperties)
    , mpCoordinateTransformation( NLGeom ?
                                  new ShellT3_CorotationalCoordinateTransformation(pGeometry) :
                                  new ShellT3_CoordinateTransformation(pGeometry))
{
  mThisIntegrationMethod = OPT_INTEGRATION_METHOD;
}

ShellThinElement3D3N::ShellThinElement3D3N(IndexType NewId,
                                           GeometryType::Pointer pGeometry,
                                           PropertiesType::Pointer pProperties,
                                           CoordinateTransformationBasePointerType pCoordinateTransformation)
    : Element(NewId, pGeometry, pProperties)
    , mpCoordinateTransformation(pCoordinateTransformation)
{
  mThisIntegrationMethod = OPT_INTEGRATION_METHOD;
}

ShellThinElement3D3N::~ShellThinElement3D3N()
{
}

Element::Pointer ShellThinElement3D3N::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
{
  GeometryType::Pointer newGeom( GetGeometry().Create(ThisNodes) );
  return Kratos::make_shared< ShellThinElement3D3N >(NewId, newGeom, pProperties, mpCoordinateTransformation->Create(newGeom));
}

ShellThinElement3D3N::IntegrationMethod ShellThinElement3D3N::GetIntegrationMethod() const
{
  return mThisIntegrationMethod;
}

void ShellThinElement3D3N::Initialize()
{
  KRATOS_TRY

      const GeometryType & geom = GetGeometry();
  const PropertiesType & props = GetProperties();

  if(geom.PointsNumber() != OPT_NUM_NODES)
    KRATOS_THROW_ERROR(std::logic_error, "ShellThinElement3D3N Element - Wrong number of nodes", geom.PointsNumber());

  const GeometryType::IntegrationPointsArrayType & integrationPoints = geom.IntegrationPoints(GetIntegrationMethod());
  if(integrationPoints.size() != OPT_NUM_GP)
    KRATOS_THROW_ERROR(std::logic_error, "ShellThinElement3D3N Element - Wrong integration scheme", integrationPoints.size());

  if(mSections.size() != OPT_NUM_GP)
  {
    const Matrix & shapeFunctionsValues = geom.ShapeFunctionsValues(GetIntegrationMethod());

    ShellCrossSection::Pointer theSection;
    if(props.Has(SHELL_CROSS_SECTION))
    {
      theSection = props[SHELL_CROSS_SECTION];
    }
    else
    {
      theSection = Kratos::make_shared<ShellCrossSection>();
      theSection->BeginStack();
      theSection->AddPly(props[THICKNESS], 0.0, 5, this->pGetProperties());
      theSection->EndStack();
    }

    mSections.clear();
    for(int i = 0; i < OPT_NUM_GP; i++)
    {
      ShellCrossSection::Pointer sectionClone = theSection->Clone();
      sectionClone->SetSectionBehavior(ShellCrossSection::Thin);
      sectionClone->InitializeCrossSection(props, geom, row( shapeFunctionsValues, i ));
      mSections.push_back(sectionClone);
    }
  }

  mpCoordinateTransformation->Initialize();

  this->SetupOrientationAngles();

  KRATOS_CATCH("")
      }

void ShellThinElement3D3N::ResetConstitutiveLaw()
{
  KRATOS_TRY

      const GeometryType & geom = GetGeometry();
  const Matrix & shapeFunctionsValues = geom.ShapeFunctionsValues(GetIntegrationMethod());

  const Properties& props = GetProperties();
  for(SizeType i = 0; i < mSections.size(); i++)
    mSections[i]->ResetCrossSection(props, geom, row(shapeFunctionsValues, i));

  KRATOS_CATCH("")
      }

void ShellThinElement3D3N::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo)
{
  if(rResult.size() != OPT_NUM_DOFS)
    rResult.resize(OPT_NUM_DOFS, false);

  GeometryType & geom = this->GetGeometry();

  for(SizeType i = 0; i < geom.size(); i++)
  {
    int index = i * 6;
    NodeType & iNode = geom[i];

    rResult[index]     = iNode.GetDof(DISPLACEMENT_X).EquationId();
    rResult[index + 1] = iNode.GetDof(DISPLACEMENT_Y).EquationId();
    rResult[index + 2] = iNode.GetDof(DISPLACEMENT_Z).EquationId();

    rResult[index + 3] = iNode.GetDof(ROTATION_X).EquationId();
    rResult[index + 4] = iNode.GetDof(ROTATION_Y).EquationId();
    rResult[index + 5] = iNode.GetDof(ROTATION_Z).EquationId();
  }
}

void ShellThinElement3D3N::GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& CurrentProcessInfo)
{
  ElementalDofList.resize(0);
  ElementalDofList.reserve(OPT_NUM_DOFS);

  GeometryType & geom = this->GetGeometry();

  for (SizeType i = 0; i < geom.size(); i++)
  {
    NodeType & iNode = geom[i];

    ElementalDofList.push_back(iNode.pGetDof(DISPLACEMENT_X));
    ElementalDofList.push_back(iNode.pGetDof(DISPLACEMENT_Y));
    ElementalDofList.push_back(iNode.pGetDof(DISPLACEMENT_Z));

    ElementalDofList.push_back(iNode.pGetDof(ROTATION_X));
    ElementalDofList.push_back(iNode.pGetDof(ROTATION_Y));
    ElementalDofList.push_back(iNode.pGetDof(ROTATION_Z));
  }
}

int ShellThinElement3D3N::Check(const ProcessInfo& rCurrentProcessInfo)
{
  KRATOS_TRY

      GeometryType& geom = GetGeometry();

  // verify that the variables are correctly initialized
  if(DISPLACEMENT.Key() == 0)
    KRATOS_THROW_ERROR(std::invalid_argument,"DISPLACEMENT has Key zero! (check if the application is correctly registered","");
  if(ROTATION.Key() == 0)
    KRATOS_THROW_ERROR(std::invalid_argument,"ROTATION has Key zero! (check if the application is correctly registered","");
  if(VELOCITY.Key() == 0)
    KRATOS_THROW_ERROR(std::invalid_argument,"VELOCITY has Key zero! (check if the application is correctly registered","");
  if(ACCELERATION.Key() == 0)
    KRATOS_THROW_ERROR(std::invalid_argument,"ACCELERATION has Key zero! (check if the application is correctly registered","");
  if(DENSITY.Key() == 0)
    KRATOS_THROW_ERROR(std::invalid_argument,"DENSITY has Key zero! (check if the application is correctly registered","");
  if(SHELL_CROSS_SECTION.Key() == 0)
    KRATOS_THROW_ERROR(std::invalid_argument,"SHELL_CROSS_SECTION has Key zero! (check if the application is correctly registered","");
  if(THICKNESS.Key() == 0)
    KRATOS_THROW_ERROR(std::invalid_argument,"THICKNESS has Key zero! (check if the application is correctly registered","");
  if(CONSTITUTIVE_LAW.Key() == 0)
    KRATOS_THROW_ERROR(std::invalid_argument,"CONSTITUTIVE_LAW has Key zero! (check if the application is correctly registered","");

  // verify that the dofs exist
  for(unsigned int i=0; i<geom.size(); i++)
  {
    if(geom[i].SolutionStepsDataHas(DISPLACEMENT) == false)
      KRATOS_THROW_ERROR(std::invalid_argument,"missing variable DISPLACEMENT on node ",geom[i].Id());
    if(geom[i].HasDofFor(DISPLACEMENT_X) == false || geom[i].HasDofFor(DISPLACEMENT_Y) == false || geom[i].HasDofFor(DISPLACEMENT_Z) == false)
      KRATOS_THROW_ERROR(std::invalid_argument,"missing one of the dofs for the variable DISPLACEMENT on node ",GetGeometry()[i].Id());
    if(geom[i].SolutionStepsDataHas(ROTATION) == false)
      KRATOS_THROW_ERROR(std::invalid_argument,"missing variable ROTATION on node ",geom[i].Id());
    if(geom[i].HasDofFor(ROTATION_X) == false || geom[i].HasDofFor(ROTATION_Y) == false || geom[i].HasDofFor(ROTATION_Z) == false)
      KRATOS_THROW_ERROR(std::invalid_argument,"missing one of the dofs for the variable ROTATION on node ",geom[i].Id());

    if(geom[i].GetBufferSize() < 2)
      KRATOS_THROW_ERROR(std::logic_error, "This Element needs at least a buffer size = 2", "");
  }

  // check properties
  if(this->pGetProperties() == NULL)
    KRATOS_THROW_ERROR(std::logic_error, "Properties not provided for element ", this->Id());

  const PropertiesType & props = this->GetProperties();

  if(props.Has(SHELL_CROSS_SECTION)) // if the user specified a cross section ...
  {
    const ShellCrossSection::Pointer & section = props[SHELL_CROSS_SECTION];
    if(section == NULL)
      KRATOS_THROW_ERROR(std::logic_error,"SHELL_CROSS_SECTION not provided for element ",this->Id());

    section->Check(props, geom, rCurrentProcessInfo);
  }
  else // ... allow the automatic creation of a homogeneous section from a material and a thickness
  {
    if(!props.Has(CONSTITUTIVE_LAW))
      KRATOS_THROW_ERROR(std::logic_error, "CONSTITUTIVE_LAW not provided for element ", this->Id());
    const ConstitutiveLaw::Pointer& claw = props[CONSTITUTIVE_LAW];
    if(claw == NULL)
      KRATOS_THROW_ERROR(std::logic_error, "CONSTITUTIVE_LAW not provided for element ", this->Id());

    if(!props.Has(THICKNESS))
      KRATOS_THROW_ERROR(std::logic_error, "THICKNESS not provided for element ", this->Id());
    if(props[THICKNESS] <= 0.0)
      KRATOS_THROW_ERROR(std::logic_error, "wrong THICKNESS value provided for element ", this->Id());

    ShellCrossSection::Pointer dummySection = Kratos::make_shared<ShellCrossSection>();
    dummySection->BeginStack();
    dummySection->AddPly(props[THICKNESS], 0.0, 5, this->pGetProperties());
    dummySection->EndStack();
    dummySection->SetSectionBehavior(ShellCrossSection::Thin);
    dummySection->Check(props, geom, rCurrentProcessInfo);
  }

  return 0;

  KRATOS_CATCH("")
      }

void ShellThinElement3D3N::CleanMemory()
{
}

void ShellThinElement3D3N::GetValuesVector(Vector& values, int Step)
{
  if(values.size() != OPT_NUM_DOFS)
    values.resize(OPT_NUM_DOFS, false);

  const GeometryType & geom = GetGeometry();

  for (SizeType i = 0; i < geom.size(); i++)
  {
    const NodeType & iNode = geom[i];
    const array_1d<double,3>& disp = iNode.FastGetSolutionStepValue(DISPLACEMENT, Step);
    const array_1d<double,3>& rot = iNode.FastGetSolutionStepValue(ROTATION, Step);

    int index = i*6;
    values[index]     = disp[0];
    values[index + 1] = disp[1];
    values[index + 2] = disp[2];

    values[index + 3] = rot[0];
    values[index + 4] = rot[1];
    values[index + 5] = rot[2];
  }
}

void ShellThinElement3D3N::GetFirstDerivativesVector(Vector& values, int Step)
{
  if(values.size() != OPT_NUM_DOFS)
    values.resize(OPT_NUM_DOFS,false);

  const GeometryType & geom = GetGeometry();

  for (SizeType i = 0; i < geom.size(); i++)
  {
    const NodeType & iNode = geom[i];
    const array_1d<double,3>& vel = iNode.FastGetSolutionStepValue(VELOCITY, Step);

    int index = i * 6;
    values[index]        = vel[0];
    values[index + 1]    = vel[1];
    values[index + 2]    = vel[2];
    values[index + 3]    = 0.0;
    values[index + 4]    = 0.0;
    values[index + 5]    = 0.0;
  }
}

void ShellThinElement3D3N::GetSecondDerivativesVector(Vector& values, int Step)
{
  if(values.size() != OPT_NUM_DOFS)
    values.resize(OPT_NUM_DOFS,false);

  const GeometryType & geom = GetGeometry();

  for (SizeType i = 0; i < geom.size(); i++)
  {
    const NodeType & iNode = geom[i];
    const array_1d<double,3>& acc = iNode.FastGetSolutionStepValue(ACCELERATION, Step);

    int index = i * 6;
    values[index]        = acc[0];
    values[index + 1]    = acc[1];
    values[index + 2]    = acc[2];
    values[index + 3]    = 0.0;
    values[index + 4]    = 0.0;
    values[index + 5]    = 0.0;
  }
}

void ShellThinElement3D3N::InitializeNonLinearIteration(ProcessInfo& CurrentProcessInfo)
{
  mpCoordinateTransformation->InitializeNonLinearIteration(CurrentProcessInfo);

  const GeometryType & geom = this->GetGeometry();
  const Matrix & shapeFunctionsValues = geom.ShapeFunctionsValues(GetIntegrationMethod());
  for(SizeType i = 0; i < mSections.size(); i++)
    mSections[i]->InitializeNonLinearIteration(GetProperties(), geom, row(shapeFunctionsValues, i), CurrentProcessInfo);
}

void ShellThinElement3D3N::FinalizeNonLinearIteration(ProcessInfo& CurrentProcessInfo)
{
  mpCoordinateTransformation->FinalizeNonLinearIteration(CurrentProcessInfo);

  const GeometryType & geom = this->GetGeometry();
  const Matrix & shapeFunctionsValues = geom.ShapeFunctionsValues(GetIntegrationMethod());
  for(SizeType i = 0; i < mSections.size(); i++)
    mSections[i]->FinalizeNonLinearIteration(GetProperties(), geom, row(shapeFunctionsValues, i), CurrentProcessInfo);
}

void ShellThinElement3D3N::InitializeSolutionStep(ProcessInfo& CurrentProcessInfo)
{
  const PropertiesType& props = GetProperties();
  const GeometryType & geom = GetGeometry();
  const Matrix & shapeFunctionsValues = geom.ShapeFunctionsValues(GetIntegrationMethod());

  for(SizeType i = 0; i < mSections.size(); i++)
    mSections[i]->InitializeSolutionStep(props, geom, row(shapeFunctionsValues, i), CurrentProcessInfo);

  mpCoordinateTransformation->InitializeSolutionStep(CurrentProcessInfo);
}

void ShellThinElement3D3N::FinalizeSolutionStep(ProcessInfo& CurrentProcessInfo)
{
  const PropertiesType& props = GetProperties();
  const GeometryType& geom = GetGeometry();
  const Matrix & shapeFunctionsValues = geom.ShapeFunctionsValues(GetIntegrationMethod());

  for(SizeType i = 0; i < mSections.size(); i++)
    mSections[i]->FinalizeSolutionStep(props, geom, row(shapeFunctionsValues, i), CurrentProcessInfo);

  mpCoordinateTransformation->FinalizeSolutionStep(CurrentProcessInfo);
}

void ShellThinElement3D3N::CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo)
{
  if((rMassMatrix.size1() != OPT_NUM_DOFS) || (rMassMatrix.size2() != OPT_NUM_DOFS))
    rMassMatrix.resize(OPT_NUM_DOFS, OPT_NUM_DOFS, false);
  noalias(rMassMatrix) = ZeroMatrix(OPT_NUM_DOFS, OPT_NUM_DOFS);

  // Compute the local coordinate system.

  ShellT3_LocalCoordinateSystem referenceCoordinateSystem(
      mpCoordinateTransformation->CreateReferenceCoordinateSystem() );

  // lumped area

  double lump_area = referenceCoordinateSystem.Area() / 3.0;

  // Calculate avarage mass per unit area
  double av_mass_per_unit_area = 0.0;
  for(size_t i = 0; i < OPT_NUM_GP; i++)
    av_mass_per_unit_area += mSections[i]->CalculateMassPerUnitArea();
  av_mass_per_unit_area /= double(OPT_NUM_GP);

  // loop on nodes
  for(size_t i = 0; i < 3; i++)
  {
    size_t index = i * 6;

    double nodal_mass = av_mass_per_unit_area * lump_area;

    // translational mass
    rMassMatrix(index, index)            = nodal_mass;
    rMassMatrix(index + 1, index + 1)    = nodal_mass;
    rMassMatrix(index + 2, index + 2)    = nodal_mass;

    // rotational mass - neglected for the moment...
  }
}

void ShellThinElement3D3N::CalculateDampingMatrix(MatrixType& rDampingMatrix, ProcessInfo& rCurrentProcessInfo)
{
  if((rDampingMatrix.size1() != OPT_NUM_DOFS) || (rDampingMatrix.size2() != OPT_NUM_DOFS))
    rDampingMatrix.resize(OPT_NUM_DOFS, OPT_NUM_DOFS, false);

  noalias( rDampingMatrix ) = ZeroMatrix(OPT_NUM_DOFS, OPT_NUM_DOFS);
}

void ShellThinElement3D3N::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                                                VectorType& rRightHandSideVector,
                                                ProcessInfo& rCurrentProcessInfo)
{
  CalculateAll(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo, true, true);
}

void ShellThinElement3D3N::CalculateRightHandSide(VectorType& rRightHandSideVector,
                                                  ProcessInfo& rCurrentProcessInfo)
{
  Matrix dummy;
  CalculateAll(dummy, rRightHandSideVector, rCurrentProcessInfo, true, true);
}

// =====================================================================================
//
// Class ShellThinElement3D3N - Results on Gauss Points
//
// =====================================================================================

void ShellThinElement3D3N::GetValueOnIntegrationPoints(const Variable<double>& rVariable,
                                                       std::vector<double>& rValues,
                                                       const ProcessInfo& rCurrentProcessInfo)
{
  if(rValues.size() != OPT_NUM_GP)
    rValues.resize(OPT_NUM_GP);

  for(int i = 0; i < OPT_NUM_GP; i++)
    mSections[i]->GetValue(rVariable, rValues[i]);

  OPT_INTERPOLATE_RESULTS_TO_STANDARD_GAUSS_POINTS(rValues);
}

void ShellThinElement3D3N::GetValueOnIntegrationPoints(const Variable<Vector>& rVariable,
                                                       std::vector<Vector>& rValues,
                                                       const ProcessInfo& rCurrentProcessInfo)
{
}

void ShellThinElement3D3N::GetValueOnIntegrationPoints(const Variable<Matrix>& rVariable,
                                                       std::vector<Matrix>& rValues,
                                                       const ProcessInfo& rCurrentProcessInfo)
{
  if(TryGetValueOnIntegrationPoints_GeneralizedStrainsOrStresses(rVariable, rValues, rCurrentProcessInfo)) return;
}

void ShellThinElement3D3N::GetValueOnIntegrationPoints(const Variable<array_1d<double,3> >& rVariable,
                                                       std::vector<array_1d<double,3> >& rValues,
                                                       const ProcessInfo& rCurrentProcessInfo)
{
  if(TryGetValueOnIntegrationPoints_MaterialOrientation(rVariable, rValues, rCurrentProcessInfo)) return;
}

void ShellThinElement3D3N::GetValueOnIntegrationPoints(const Variable<array_1d<double,6> >& rVariable,
                                                       std::vector<array_1d<double,6> >& rValues,
                                                       const ProcessInfo& rCurrentProcessInfo)
{
}

void ShellThinElement3D3N::SetCrossSectionsOnIntegrationPoints(std::vector< ShellCrossSection::Pointer >& crossSections)
{
  KRATOS_TRY
      if(crossSections.size() != OPT_NUM_GP)
        KRATOS_THROW_ERROR(std::logic_error, "The number of cross section is wrong", crossSections.size());
  mSections.clear();
  for(SizeType i = 0; i <crossSections.size(); i++)
    mSections.push_back(crossSections[i]);
  this->SetupOrientationAngles();
  KRATOS_CATCH("")
      }

// =====================================================================================
//
// Class ShellThinElement3D3N - Private methods
//
// =====================================================================================

void ShellThinElement3D3N::DecimalCorrection(Vector& a)
{
  double norm = norm_2(a);
  double tolerance = std::max(norm * 1.0E-12, 1.0E-12);
  for(SizeType i = 0; i < a.size(); i++)
    if(std::abs(a(i)) < tolerance)
      a(i) = 0.0;
}

void ShellThinElement3D3N::SetupOrientationAngles()
{
  ShellT3_LocalCoordinateSystem lcs( mpCoordinateTransformation->CreateReferenceCoordinateSystem() );

  Vector3Type normal;
  noalias( normal ) = lcs.Vz();

  Vector3Type dZ;
  dZ(0) = 0.0;
  dZ(1) = 0.0;
  dZ(2) = 1.0; // for the moment let's take this. But the user can specify its own triad! TODO

  Vector3Type dirX;
  MathUtils<double>::CrossProduct(dirX,   dZ, normal);

  // try to normalize the x vector. if it is near zero it means that we need
  // to choose a default one.
  double dirX_norm = dirX(0)*dirX(0) + dirX(1)*dirX(1) + dirX(2)*dirX(2);
  if(dirX_norm < 1.0E-12) {
    dirX(0) = 1.0;
    dirX(1) = 0.0;
    dirX(2) = 0.0;
  }
  else if(dirX_norm != 1.0) {
    dirX_norm = std::sqrt(dirX_norm);
    dirX /= dirX_norm;
  }

  Vector3Type elem_dirX = lcs.Vx();

  // now calculate the angle between the element x direction and the material x direction.
  Vector3Type& a = elem_dirX;
  Vector3Type& b = dirX;
  double a_dot_b = a(0)*b(0) + a(1)*b(1) + a(2)*b(2);
  if(a_dot_b < -1.0) a_dot_b = -1.0;
  if(a_dot_b >  1.0) a_dot_b =  1.0;
  double angle = std::acos( a_dot_b );

  // if they are not counter-clock-wise, let's change the sign of the angle
  if(angle != 0.0) {
    const MatrixType& R = lcs.Orientation();
    if( dirX(0)*R(1, 0) + dirX(1)*R(1, 1) + dirX(2)*R(1, 2) < 0.0 )
      angle = -angle;
  }

  for(CrossSectionContainerType::iterator it = mSections.begin(); it != mSections.end(); ++it)
    (*it)->SetOrientationAngle(angle);
}

void ShellThinElement3D3N::InitializeCalculationData(CalculationData& data)
{
  //-------------------------------------
  // Computation of all stuff that remain
  // constant throughout the calculations

  //-------------------------------------
  // geometry data

  const double x12 = data.LCS0.X1() - data.LCS0.X2();
  const double x23 = data.LCS0.X2() - data.LCS0.X3();
  const double x31 = data.LCS0.X3() - data.LCS0.X1();
  const double x21 = -x12;
  const double x32 = -x23;
  const double x13 = -x31;

  const double y12 = data.LCS0.Y1() - data.LCS0.Y2();
  const double y23 = data.LCS0.Y2() - data.LCS0.Y3();
  const double y31 = data.LCS0.Y3() - data.LCS0.Y1();
  const double y21 = -y12;
  const double y32 = -y23;
  const double y13 = -y31;

  const double A  = 0.5*(y21*x13 - x21*y13);
  const double A2 = 2.0*A;
  const double A4 = 4.0*A;
  const double AA4 = A * A4;

  const double LL21 = x21*x21 + y21*y21;
  const double LL32 = x32*x32 + y32*y32;
  const double LL13 = x13*x13 + y13*y13;

  // Note: here we compute the avarage thickness,
  // since L is constant over the element.
  // Now it is not necessary to compute the avarage
  // because the current implementation of the cross section
  // doesn't have a variable thickness
  // (for example as a function of the spatial coordinates...).
  // This is just a place-holder for future
  // implementation of a variable thickness

  double h = 0.0;
  for(unsigned int i = 0; i < mSections.size(); i++)
    h += mSections[i]->GetThickness();
  h /= (double)mSections.size();

  data.hMean = h;
  data.TotalArea = A;
  data.TotalVolume = A * h;

  // this is the integration weight
  // used during the gauss loop.
  // it is dArea because it will
  // multiply section stress resultants
  // and section constitutive matrices
  // that already take into accout the
  // thickness

  data.dA = A / (double)OPT_NUM_GP;

  // crete the integration point locations
  if(data.gpLocations.size() != 0) data.gpLocations.clear();
  data.gpLocations.resize( OPT_NUM_GP );
#ifdef OPT_1_POINT_INTEGRATION
  array_1d<double,3>& gp0 = data.gpLocations[0];
  gp0[0] = 1.0/3.0; gp0[1] = 1.0/3.0; gp0[2] = 1.0/3.0;
#else
  array_1d<double,3>& gp0 = data.gpLocations[0];
  array_1d<double,3>& gp1 = data.gpLocations[1];
  array_1d<double,3>& gp2 = data.gpLocations[2];
#ifdef OPT_USES_INTERIOR_GAUSS_POINTS
  gp0[0] = 1.0/6.0; gp0[1] = 1.0/6.0; gp0[2] = 2.0/3.0;
  gp1[0] = 2.0/3.0; gp1[1] = 1.0/6.0; gp1[2] = 1.0/6.0;
  gp2[0] = 1.0/6.0; gp2[1] = 2.0/3.0; gp2[2] = 1.0/6.0;
#else
  gp0[0] = 0.5; gp0[1] = 0.5; gp0[2] = 0.0;
  gp1[0] = 0.0; gp1[1] = 0.5; gp1[2] = 0.5;
  gp2[0] = 0.5; gp2[1] = 0.0; gp2[2] = 0.5;
#endif // OPT_USES_INTERIOR_GAUSS_POINTS

#endif // OPT_1_POINT_INTEGRATION

  // cartesian derivatives
  data.dNxy.resize(3, 2, false);
  data.dNxy(0, 0) = (y13 - y12)/A2; data.dNxy(0, 1) = (x12 - x13)/A2;
  data.dNxy(1, 0) =        -y13/A2; data.dNxy(1, 1) =         x13/A2;
  data.dNxy(2, 0) =         y12/A2; data.dNxy(2, 1) =        -x12/A2;

  //-------------------------------------
  // template parameters

  const double b1    =  1.0;
  const double b2    =  2.0;
  const double b3    =  1.0;
  const double b4    =  0.0;
  const double b5    =  1.0;
  const double b6    = -1.0;
  const double b7    = -1.0;
  const double b8    = -1.0;
  const double b9    = -2.0;

  const double alpha   =  1.5;
  const double alpha_6 =  alpha / 6.0;

  //--------------------------------------
  // calculate L - Lumping matrix
  // for the construction of the basic
  // stiffness

  double L_mult = 0.5/A;

  data.L.resize(3, 9, false);

  data.L(0, 0) = L_mult * y23;
  data.L(1, 0) = 0.00;
  data.L(2, 0) = L_mult * x32;
  data.L(0, 1) = 0.00;
  data.L(1, 1) = L_mult * x32;
  data.L(2, 1) = L_mult * y23;
  data.L(0, 2) = L_mult * y23*(y13-y21)*alpha_6;
  data.L(1, 2) = L_mult * x32*(x31-x12)*alpha_6;
  data.L(2, 2) = L_mult * 2.00*(x31*y13-x12*y21)*alpha_6;
  data.L(0, 3) = L_mult * y31;
  data.L(1, 3) = 0.00;
  data.L(2, 3) = L_mult * x13;
  data.L(0, 4) = 0.00;
  data.L(1, 4) = L_mult * x13;
  data.L(2, 4) = L_mult * y31;
  data.L(0, 5) = L_mult * y31*(y21-y32)*alpha_6;
  data.L(1, 5) = L_mult * x13*(x12-x23)*alpha_6;
  data.L(2, 5) = L_mult * 2.00*(x12*y21-x23*y32)*alpha_6;
  data.L(0, 6) = L_mult * y12;
  data.L(1, 6) = 0.00;
  data.L(2, 6) = L_mult * x21;
  data.L(0, 7) = 0.00;
  data.L(1, 7) = L_mult * x21;
  data.L(2, 7) = L_mult * y12;
  data.L(0, 8) = L_mult * y12*(y32-y13)*alpha_6;
  data.L(1, 8) = L_mult * x21*(x23-x31)*alpha_6;
  data.L(2, 8) = L_mult * 2.00*(x23*y32-x31*y13)*alpha_6;

  //--------------------------------------
  // calculate Q1,Q2,Q3 - matrices
  // for the construction of the
  // higher order stiffness

  data.Q1.resize(3, 3, false);

  data.Q1(0,0) = b1*A2/(LL21*3.00);
  data.Q1(0,1) = b2*A2/(LL21*3.00);
  data.Q1(0,2) = b3*A2/(LL21*3.00);
  data.Q1(1,0) = b4*A2/(LL32*3.00);
  data.Q1(1,1) = b5*A2/(LL32*3.00);
  data.Q1(1,2) = b6*A2/(LL32*3.00);
  data.Q1(2,0) = b7*A2/(LL13*3.00);
  data.Q1(2,1) = b8*A2/(LL13*3.00);
  data.Q1(2,2) = b9*A2/(LL13*3.00);

  data.Q2.resize(3, 3, false);

  data.Q2(0,0) = b9*A2/(LL21*3.00);
  data.Q2(0,1) = b7*A2/(LL21*3.00);
  data.Q2(0,2) = b8*A2/(LL21*3.00);
  data.Q2(1,0) = b3*A2/(LL32*3.00);
  data.Q2(1,1) = b1*A2/(LL32*3.00);
  data.Q2(1,2) = b2*A2/(LL32*3.00);
  data.Q2(2,0) = b6*A2/(LL13*3.00);
  data.Q2(2,1) = b4*A2/(LL13*3.00);
  data.Q2(2,2) = b5*A2/(LL13*3.00);

  data.Q3.resize(3, 3, false);

  data.Q3(0,0) = b5*A2/(LL21*3.00);
  data.Q3(0,1) = b6*A2/(LL21*3.00);
  data.Q3(0,2) = b4*A2/(LL21*3.00);
  data.Q3(1,0) = b8*A2/(LL32*3.00);
  data.Q3(1,1) = b9*A2/(LL32*3.00);
  data.Q3(1,2) = b7*A2/(LL32*3.00);
  data.Q3(2,0) = b2*A2/(LL13*3.00);
  data.Q3(2,1) = b3*A2/(LL13*3.00);
  data.Q3(2,2) = b1*A2/(LL13*3.00);

  //--------------------------------------
  // calculate Te, TTu -
  // transformation matrices
  // for the construction of the
  // higher order stiffness

  data.Te.resize(3, 3, false);

  data.Te(0,0) = 1.0/AA4 * y23*y13*LL21;
  data.Te(0,1) = 1.0/AA4 * y31*y21*LL32;
  data.Te(0,2) = 1.0/AA4 * y12*y32*LL13;
  data.Te(1,0) = 1.0/AA4 * x23*x13*LL21;
  data.Te(1,1) = 1.0/AA4 * x31*x21*LL32;
  data.Te(1,2) = 1.0/AA4 * x12*x32*LL13;
  data.Te(2,0) = 1.0/AA4 * (y23*x31+x32*y13)*LL21;
  data.Te(2,1) = 1.0/AA4 * (y31*x12+x13*y21)*LL32;
  data.Te(2,2) = 1.0/AA4 * (y12*x23+x21*y32)*LL13;

  data.TTu.resize(3, 9, false);

  for(unsigned int i=0; i<3; i++)
  {
    data.TTu(i, 0) = 1.0/A4 * x32;
    data.TTu(i, 1) = 1.0/A4 * y32;
    data.TTu(i, 2) = 0.0;
    data.TTu(i, 3) = 1.0/A4 * x13;
    data.TTu(i, 4) = 1.0/A4 * y13;
    data.TTu(i, 5) = 0.0;
    data.TTu(i, 6) = 1.0/A4 * x21;
    data.TTu(i, 7) = 1.0/A4 * y21;
    data.TTu(i, 8) = 0.0;
  }
  data.TTu(0, 2) = 1.0;
  data.TTu(1, 5) = 1.0;
  data.TTu(2, 8) = 1.0;

  //--------------------------------------
  // calculate the displacement vector
  // in global and local coordinate systems

  data.globalDisplacements.resize(OPT_NUM_DOFS, false);
  GetValuesVector( data.globalDisplacements );

  data.localDisplacements =
      mpCoordinateTransformation->CalculateLocalDisplacements(
          data.LCS, data.globalDisplacements);

  //--------------------------------------
  // Finally allocate all auxiliary
  // matrices to be used later on
  // during the element integration.
  // Just to avoid re-allocations

  data.B.resize(OPT_STRAIN_SIZE, OPT_NUM_DOFS, false);
  data.D.resize(OPT_STRAIN_SIZE, OPT_STRAIN_SIZE, false);
  data.BTD.resize(OPT_NUM_DOFS, OPT_STRAIN_SIZE, false);

  data.generalizedStrains.resize(OPT_STRAIN_SIZE, false);
  data.generalizedStresses.resize(OPT_STRAIN_SIZE, false);

  data.N.resize(3, false);

  data.Q.resize(3, 3, false);
  data.Qh.resize(3, 9, false);
  data.TeQ.resize(3, 3, false);

  data.H1.resize(9, false);
  data.H2.resize(9, false);
  data.H3.resize(9, false);
  data.H4.resize(9, false);
  data.Bb.resize(3, 9, false);

  //--------------------------------------
  // Initialize the section parameters

  data.SectionParameters.SetElementGeometry( GetGeometry() );
  data.SectionParameters.SetMaterialProperties( GetProperties() );
  data.SectionParameters.SetProcessInfo( data.CurrentProcessInfo );

  data.SectionParameters.SetGeneralizedStrainVector( data.generalizedStrains );
  data.SectionParameters.SetGeneralizedStressVector( data.generalizedStresses );
  data.SectionParameters.SetConstitutiveMatrix( data.D );

  data.SectionParameters.SetShapeFunctionsDerivatives( data.dNxy );

  Flags& options = data.SectionParameters.GetOptions();
  options.Set(ConstitutiveLaw::COMPUTE_STRESS, data.CalculateRHS);
  options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, data.CalculateLHS);
}

void ShellThinElement3D3N::CalculateBMatrix(CalculationData& data)
{
  //---------------------------------------------
  // geom data
  array_1d<double, 3>& gpLoc = data.gpLocations[data.gpIndex];
  double loc1 = gpLoc[0];
  double loc2 = gpLoc[1];
  double loc3 = gpLoc[2];

  const double x12 = data.LCS0.X1() - data.LCS0.X2();
  const double x23 = data.LCS0.X2() - data.LCS0.X3();
  const double x31 = data.LCS0.X3() - data.LCS0.X1();
  const double y12 = data.LCS0.Y1() - data.LCS0.Y2();
  const double y23 = data.LCS0.Y2() - data.LCS0.Y3();
  const double y31 = data.LCS0.Y3() - data.LCS0.Y1();

  //---------------------------------------------
  // set to zero the total B matrix
  data.B.clear();

  //---------------------------------------------
  // membrane basic part L
  // already computed.it is constant over the
  // element

  //---------------------------------------------
  // membrane higher order part Qh
  noalias( data.Q )  = loc1 * data.Q1;
  noalias( data.Q ) += loc2 * data.Q2;
  noalias( data.Q ) += loc3 * data.Q3;
  noalias( data.TeQ ) = prod( data.Te, data.Q );
  noalias( data.Qh  ) = /*1.5**/std::sqrt(data.beta0) * prod( data.TeQ, data.TTu );

  //---------------------------------------------
  // compute the bending B matrix (DKT)

  const double LL12 = x12*x12 + y12*y12;
  const double LL23 = x23*x23 + y23*y23;
  const double LL31 = x31*x31 + y31*y31;

  const double p4 = -6.0*x23/(LL23);
  const double p5 = -6.0*x31/(LL31);
  const double p6 = -6.0*x12/(LL12);

  const double t4 = -6.0*y23/(LL23);
  const double t5 = -6.0*y31/(LL31);
  const double t6 = -6.0*y12/(LL12);

  const double q4 = 3.0*x23*y23/(LL23);
  const double q5 = 3.0*x31*y31/(LL31);
  const double q6 = 3.0*x12*y12/(LL12);

  const double r4 = 3.0*y23*y23/(LL23);
  const double r5 = 3.0*y31*y31/(LL31);
  const double r6 = 3.0*y12*y12/(LL12);

  const double area = 0.5*(y12*x31 - x12*y31);

  data.H1[0] = p6*(1.0-2.0*loc2) + (p5-p6)*loc3;
  data.H1[1] = q6*(1.0-2.0*loc2) - (q5+q6)*loc3;
  data.H1[2]= -4.0+6.0*(loc2+loc3) + r6*(1.0-2*loc2) - loc3*(r5+r6);
  data.H1[3] = -p6*(1.0-2.0*loc2) + (p4+p6)*loc3;
  data.H1[4] = q6*(1.0-2.0*loc2) - (q6-q4)*loc3;
  data.H1[5]= -2.0+6.0*loc2  + r6*(1-2.0*loc2) + loc3*(r4-r6);
  data.H1[6] = -(p4+p5)*loc3;
  data.H1[7] = (q4-q5)*loc3;
  data.H1[8]= -loc3*(r5-r4);

  data.H2[0] = t6*(1.0-2.0*loc2) + (t5-t6)*loc3;
  data.H2[1] = 1.0 + r6*(1.0-2.0*loc2) - loc3*(r5+r6);
  data.H2[2]= -q6*(1.0-2.0*loc2)+loc3*(q5+q6);
  data.H2[3] = -t6*(1.0-2.0*loc2) + (t4+t6)*loc3;
  data.H2[4] = -1.0 + r6*(1-2*loc2) + loc3*(r4-r6);
  data.H2[5]= -q6*(1.0-2.0*loc2)-loc3*(q4-q6);
  data.H2[6] = -(t5+t4)*loc3;
  data.H2[7] = (r4-r5)*loc3;
  data.H2[8]= -loc3*(q4-q5);

  data.H3[0] = -p5*(1.0-2.0*loc3) - (p6-p5)*loc2;
  data.H3[1] = q5*(1.0-2.0*loc3) - (q5+q6)*loc2;
  data.H3[2]= -4.0 + 6.0*(loc2+loc3) + r5*(1.0-2.0*loc3) - loc2*(r5+r6);
  data.H3[3] = loc2*(p4+p6);
  data.H3[4] = loc2*(q4-q6);
  data.H3[5]= -loc2*(r6-r4);
  data.H3[6] = p5*(1.0-2.0*loc3) - (p4+p5)*loc2;
  data.H3[7] = q5*(1.0-2.0*loc3)+loc2*(q4-q5);
  data.H3[8]= -2.0+6.0*loc3+r5*(1.0-2.0*loc3)+loc2*(r4-r5);

  data.H4[0] = -t5*(1.0-2.0*loc3) - (t6-t5)*loc2;
  data.H4[1] = 1.0+r5*(1.0-2.0*loc3)-loc2*(r5+r6);
  data.H4[2]= -q5*(1.0-2.0*loc3)+loc2*(q5+q6);
  data.H4[3] = (t4+t6)*loc2;
  data.H4[4] = (r4-r6)*loc2;
  data.H4[5]= -loc2*(q4-q6);
  data.H4[6] = t5*(1.0-2.0*loc3) - (t5+t4)*loc2;
  data.H4[7] = -1.0+r5*(1.0-2.0*loc3)+loc2*(r4-r5);
  data.H4[8]= -q5*(1.0-2.0*loc3)-loc2*(q4-q5);

  double temp = 0.5 / area;
  for(int i =0; i<9; i++)
  {
    data.Bb(0, i) = temp * ( y31*data.H1[i] + y12*data.H3[i]);
    data.Bb(1, i) = temp * (-x31*data.H2[i] - x12*data.H4[i]);
    data.Bb(2, i) = temp * (-x31*data.H1[i] - x12*data.H3[i] + y31*data.H2[i] + y12*data.H4[i] );
  }

  //---------------------------------------------
  // assemble the 2 mamebrane contributions
  // and the bending contribution
  // into the combined B matrix
  for(int nodeid = 0; nodeid < 3; nodeid++)
  {
    int i = nodeid * 3;
    int j = nodeid * 6;

    // membrane map: [0,1,5] <- [0,1,2]

    data.B(0, j  ) = data.L(0, i  ) + data.Qh(0, i  );
    data.B(0, j+1) = data.L(0, i+1) + data.Qh(0, i+1);
    data.B(0, j+5) = data.L(0, i+2) + data.Qh(0, i+2);

    data.B(1, j  ) = data.L(1, i  ) + data.Qh(1, i  );
    data.B(1, j+1) = data.L(1, i+1) + data.Qh(1, i+1);
    data.B(1, j+5) = data.L(1, i+2) + data.Qh(1, i+2);

    data.B(2, j  ) = data.L(2, i  ) + data.Qh(2, i  );
    data.B(2, j+1) = data.L(2, i+1) + data.Qh(2, i+1);
    data.B(2, j+5) = data.L(2, i+2) + data.Qh(2, i+2);

    // bending map: [2,3,4] <- [0,1,2]

    data.B(3, j+2) = data.Bb(0, i  );
    data.B(3, j+3) = data.Bb(0, i+1);
    data.B(3, j+4) = data.Bb(0, i+2);

    data.B(4, j+2) = data.Bb(1, i  );
    data.B(4, j+3) = data.Bb(1, i+1);
    data.B(4, j+4) = data.Bb(1, i+2);

    data.B(5, j+2) = data.Bb(2, i  );
    data.B(5, j+3) = data.Bb(2, i+1);
    data.B(5, j+4) = data.Bb(2, i+2);
  }
}

void ShellThinElement3D3N::CalculateBeta0(CalculationData& data)
{
  data.beta0 = 1.0; // to be changed!
}

void ShellThinElement3D3N::CalculateSectionResponse(CalculationData& data)
{
#ifdef OPT_USES_INTERIOR_GAUSS_POINTS
  const Matrix & shapeFunctions = GetGeometry().ShapeFunctionsValues(mThisIntegrationMethod);
  for(int nodeid = 0; nodeid < OPT_NUM_NODES; nodeid++)
    data.N(nodeid) = shapeFunctions(data.gpIndex, nodeid);
#else
  const array_1d<double,3>& loc = data.gpLocations[data.gpIndex];
  data.N(0) = 1.0 - loc[1] - loc[2];
  data.N(1) = loc[1];
  data.N(2) = loc[2];
#endif // !OPT_USES_INTERIOR_GAUSS_POINTS

  ShellCrossSection::Pointer& section = mSections[data.gpIndex];
  data.SectionParameters.SetShapeFunctionsValues( data.N );
  section->CalculateSectionResponse( data.SectionParameters, ConstitutiveLaw::StressMeasure_Kirchhoff );
}

void ShellThinElement3D3N::CalculateGaussPointContribution(CalculationData& data, MatrixType& LHS, VectorType& RHS)
{
  // calculate beta0
  CalculateBeta0( data );

  // calculate the total strain displ. matrix
  CalculateBMatrix( data );

  // compute generalized strains
  noalias( data.generalizedStrains ) = prod( data.B, data.localDisplacements );

  // calculate section response
  CalculateSectionResponse( data );

  // multiply the section tangent matrices and stress resultants by 'dA'
  data.D *= data.dA;
  //******
  Vector3Type& iSig = data.Sig[data.gpIndex];
  iSig(0) = data.generalizedStresses(0);
  iSig(1) = data.generalizedStresses(1);
  iSig(2) = data.generalizedStresses(2);
  //******
  data.generalizedStresses *= data.dA;

  // Add all contributions to the Stiffness Matrix
  noalias( data.BTD ) = prod( trans( data.B ), data.D );
  noalias( LHS ) += prod( data.BTD, data.B );

  // Add all contributions to the residual vector
  noalias( RHS ) -= prod( trans( data.B ), data.generalizedStresses );
}

void ShellThinElement3D3N::ApplyCorrectionToRHS(CalculationData& data, VectorType& RHS)
{
  Vector3Type meanS;
  meanS.clear();
  for(int i = 0; i < 3; i++)
    noalias( meanS ) += data.Sig[i];
  meanS /= 3.0;

  for(int i = 0; i < 3; i++)
  {
    int i1 = i;
    int i2 = i == 2 ? 0 : i+1;

    const Vector3Type& p1 = data.LCS0.Nodes()[i1];
    const Vector3Type& p2 = data.LCS0.Nodes()[i2];
    /*const Vector3Type& s1 = data.Sig[i1];
      const Vector3Type& s2 = data.Sig[i2];*/
    const Vector3Type& s1 = meanS;
    const Vector3Type& s2 = meanS;
    /*const Vector3Type& s1 = data.Sig[i];
      const Vector3Type& s2 = data.Sig[i];*/

    Vector3Type t = p2 - p1;
    Vector3Type z;
    z(0) = 0.0; z(1) = 0.0; z(2) = 1.0;
    Vector3Type n;
    MathUtils<double>::CrossProduct(n,  t,z);
    n /= MathUtils<double>::Norm3(n);

    double sx, sy;
    sx = s1(0)*n(0) + s1(2)*n(1);
    sy = s1(2)*n(0) + s1(1)*n(1);
    double q1 = std::sqrt(sx*sx + sy*sy);
    sx = s2(0)*n(0) + s2(2)*n(1);
    sy = s2(2)*n(0) + s2(1)*n(1);
    double q2 = std::sqrt(sx*sx + sy*sy);

    double q0 = (q1+q2)/2.0;

    double L = std::sqrt(t(0)*t(0) + t(1)*t(1));

    double m = 1.0/8.0*L*L*q0;

    RHS(i1*6 + 5) -= m;
    RHS(i2*6 + 5) += m;
  }
}

void ShellThinElement3D3N::AddBodyForces(CalculationData& data, VectorType& rRightHandSideVector)
{
  const GeometryType& geom = GetGeometry();

  // Get shape functions
#ifdef OPT_USES_INTERIOR_GAUSS_POINTS
  const Matrix & N = GetGeometry().ShapeFunctionsValues(mThisIntegrationMethod);
#else
  Matrix N(3,3);
  for(unsigned int igauss = 0; igauss < OPT_NUM_GP; igauss++)
  {
    const array_1d<double,3>& loc = data.gpLocations[igauss];
    N(igauss,0) = 1.0 - loc[1] - loc[2];
    N(igauss,1) = loc[1];
    N(igauss,2) = loc[2];
  }
#endif // !OPT_USES_INTERIOR_GAUSS_POINTS

  // auxiliary
  array_1d<double, 3> bf;

  // gauss loop to integrate the external force vector
  for(unsigned int igauss = 0; igauss < OPT_NUM_GP; igauss++)
  {
    // get mass per unit area
    double mass_per_unit_area = mSections[igauss]->CalculateMassPerUnitArea();

    // interpolate nodal volume accelerations to this gauss point
    // and obtain the body force vector
    bf.clear();
    for(unsigned int inode = 0; inode < 3; inode++)
    {
      if( geom[inode].SolutionStepsDataHas(VOLUME_ACCELERATION) ) //temporary, will be checked once at the beginning only
        bf += N(igauss,inode) * geom[inode].FastGetSolutionStepValue(VOLUME_ACCELERATION);
    }
    bf *= (mass_per_unit_area * data.dA);

    // add it to the RHS vector
    for(unsigned int inode = 0; inode < 3; inode++)
    {
      unsigned int index = inode*6;
      double iN = N(igauss,inode);
      rRightHandSideVector[index + 0] += iN * bf[0];
      rRightHandSideVector[index + 1] += iN * bf[1];
      rRightHandSideVector[index + 2] += iN * bf[2];
    }
  }
}

void ShellThinElement3D3N::CalculateAll(MatrixType& rLeftHandSideMatrix,
                                        VectorType& rRightHandSideVector,
                                        ProcessInfo& rCurrentProcessInfo,
                                        const bool LHSrequired,
                                        const bool RHSrequired)
{
  // Resize the Left Hand Side if necessary,
  // and initialize it to Zero

  if((rLeftHandSideMatrix.size1() != OPT_NUM_DOFS) || (rLeftHandSideMatrix.size2() != OPT_NUM_DOFS))
    rLeftHandSideMatrix.resize(OPT_NUM_DOFS, OPT_NUM_DOFS, false);
  noalias(rLeftHandSideMatrix) = ZeroMatrix(OPT_NUM_DOFS, OPT_NUM_DOFS);

  // Resize the Right Hand Side if necessary,
  // and initialize it to Zero

  if(rRightHandSideVector.size() != OPT_NUM_DOFS)
    rRightHandSideVector.resize(OPT_NUM_DOFS, false);
  noalias(rRightHandSideVector) = ZeroVector(OPT_NUM_DOFS);

  // Initialize common calculation variables

  CalculationData data(mpCoordinateTransformation, rCurrentProcessInfo);
  data.CalculateLHS = LHSrequired;
  data.CalculateRHS = RHSrequired;
  InitializeCalculationData(data);

  // Gauss Loop.

  for(size_t i = 0; i < OPT_NUM_GP; i++)
  {
    data.gpIndex = i;
    CalculateGaussPointContribution(data, rLeftHandSideMatrix, rRightHandSideVector);
  }

  //ApplyCorrectionToRHS(data, rRightHandSideVector);

  // Let the CoordinateTransformation finalize the calculation.
  // This will handle the transformation of the local matrices/vectors to
  // the global coordinate system.

  mpCoordinateTransformation->FinalizeCalculations(data.LCS,
                                                   data.globalDisplacements,
                                                   data.localDisplacements,
                                                   rLeftHandSideMatrix,
                                                   rRightHandSideVector,
                                                   RHSrequired,
                                                   LHSrequired);

  // Add body forces contributions. This doesn't depend on the coordinate system

  AddBodyForces(data, rRightHandSideVector);
}

bool ShellThinElement3D3N::TryGetValueOnIntegrationPoints_MaterialOrientation(const Variable<array_1d<double,3> >& rVariable,
                                                                              std::vector<array_1d<double,3> >& rValues,
                                                                              const ProcessInfo& rCurrentProcessInfo)
{
  // Check the required output

  int ijob = 0;
  if(rVariable == MATERIAL_ORIENTATION_DX)
    ijob = 1;
  else if(rVariable == MATERIAL_ORIENTATION_DY)
    ijob = 2;
  else if(rVariable == MATERIAL_ORIENTATION_DZ)
    ijob = 3;

  // quick return

  if(ijob == 0) return false;

  // resize output

  if(rValues.size() != OPT_NUM_GP)
    rValues.resize(OPT_NUM_GP);

  // Compute the local coordinate system.

  ShellT3_LocalCoordinateSystem localCoordinateSystem(
      mpCoordinateTransformation->CreateLocalCoordinateSystem() );

  Vector3Type eZ = localCoordinateSystem.Vz();

  // Gauss Loop

  if(ijob == 1)
  {
    Vector3Type eX = localCoordinateSystem.Vx();
    for(int i = 0; i < OPT_NUM_GP; i++)
    {
      QuaternionType q = QuaternionType::FromAxisAngle(eZ(0), eZ(1), eZ(2), mSections[i]->GetOrientationAngle());
      q.RotateVector3(eX, rValues[i]);
    }
  }
  else if(ijob == 2)
  {
    Vector3Type eY = localCoordinateSystem.Vy();
    for(int i = 0; i < OPT_NUM_GP; i++)
    {
      QuaternionType q = QuaternionType::FromAxisAngle(eZ(0), eZ(1), eZ(2), mSections[i]->GetOrientationAngle());
      q.RotateVector3(eY, rValues[i]);
    }
  }
  else if(ijob == 3)
  {
    for(int i = 0; i < OPT_NUM_GP; i++)
    {
      noalias( rValues[i] ) = eZ;
    }
  }

  return true;
}

bool ShellThinElement3D3N::TryGetValueOnIntegrationPoints_GeneralizedStrainsOrStresses(const Variable<Matrix>& rVariable,
                                                                                       std::vector<Matrix>& rValues,
                                                                                       const ProcessInfo& rCurrentProcessInfo)
{
  // Check the required output

  int ijob = 0;
  bool bGlobal = false;
  if(rVariable == SHELL_STRAIN) {
    ijob = 1;
  }
  else if(rVariable == SHELL_STRAIN_GLOBAL) {
    ijob = 1;
    bGlobal = true;
  }
  else if(rVariable == SHELL_CURVATURE) {
    ijob = 2;
  }
  else if(rVariable == SHELL_CURVATURE_GLOBAL) {
    ijob = 2;
    bGlobal = true;
  }
  else if(rVariable == SHELL_FORCE) {
    ijob = 3;
  }
  else if(rVariable == SHELL_FORCE_GLOBAL) {
    ijob = 3;
    bGlobal = true;
  }
  else if(rVariable == SHELL_MOMENT) {
    ijob = 4;
  }
  else if(rVariable == SHELL_MOMENT_GLOBAL) {
    ijob = 4;
    bGlobal = true;
  }

  // quick return

  if(ijob == 0) return false;

  // resize output

  if(rValues.size() != OPT_NUM_GP)
    rValues.resize(OPT_NUM_GP);

  // Just to store the rotation matrix for visualization purposes

  Matrix R(OPT_STRAIN_SIZE, OPT_STRAIN_SIZE);
  Matrix aux33(3, 3);

  // Initialize common calculation variables

  CalculationData data(mpCoordinateTransformation, rCurrentProcessInfo);
  data.CalculateLHS = false;
  data.CalculateRHS = true;
  InitializeCalculationData(data);

  // Gauss Loop.

  for(size_t i = 0; i < OPT_NUM_GP; i++)
  {
    // set the current integration point index
    data.gpIndex = i;
    ShellCrossSection::Pointer& section = mSections[i];

    // calculate beta0
    CalculateBeta0( data );

    // calculate the total strain displ. matrix
    CalculateBMatrix( data );

    // compute generalized strains
    noalias( data.generalizedStrains ) = prod( data.B, data.localDisplacements );

    // calculate section response
    if(ijob > 2) CalculateSectionResponse( data );

    // adjust output
    DecimalCorrection( data.generalizedStrains );
    DecimalCorrection( data.generalizedStresses );

    // store the results, but first rotate them back to the section coordinate system.
    // we want to visualize the results in that system not in the element one!
    if(section->GetOrientationAngle() != 0.0 && !bGlobal)
    {
      if (ijob > 2)
      {
        section->GetRotationMatrixForGeneralizedStresses( -(section->GetOrientationAngle()), R );
        data.generalizedStresses = prod( R, data.generalizedStresses );
      }
      else
      {
        section->GetRotationMatrixForGeneralizedStrains( -(section->GetOrientationAngle()), R );
        data.generalizedStrains = prod( R, data.generalizedStrains );
      }
    }

    // save results
    Matrix & iValue = rValues[i];
    if(iValue.size1() != 3 || iValue.size2() != 3)
      iValue.resize(3, 3, false);

    if(ijob == 1) // strains
    {
      iValue(0, 0) = data.generalizedStrains(0);
      iValue(1, 1) = data.generalizedStrains(1);
      iValue(2, 2) = 0.0;
      iValue(0, 1) = iValue(1, 0) = 0.5 * data.generalizedStrains(2);
      iValue(0, 2) = iValue(2, 0) = 0.0;
      iValue(1, 2) = iValue(2, 1) = 0.0;
    }
    else if(ijob == 2) // curvatures
    {
      iValue(0, 0) = data.generalizedStrains(3);
      iValue(1, 1) = data.generalizedStrains(4);
      iValue(2, 2) = 0.0;
      iValue(0, 1) = iValue(1, 0) = 0.5 * data.generalizedStrains(5);
      iValue(0, 2) = iValue(2, 0) = 0.0;
      iValue(1, 2) = iValue(2, 1) = 0.0;
    }
    else if(ijob == 3) // forces
    {
      iValue(0, 0) = data.generalizedStresses(0);
      iValue(1, 1) = data.generalizedStresses(1);
      iValue(2, 2) = 0.0;
      iValue(0, 1) = iValue(1, 0) = data.generalizedStresses(2);
      iValue(0, 2) = iValue(2, 0) = 0.0;
      iValue(1, 2) = iValue(2, 1) = 0.0;
    }
    else if(ijob == 4) // moments
    {
      iValue(0, 0) = data.generalizedStresses(3);
      iValue(1, 1) = data.generalizedStresses(4);
      iValue(2, 2) = 0.0;
      iValue(0, 1) = iValue(1, 0) = data.generalizedStresses(5);
      iValue(0, 2) = iValue(2, 0) = 0.0;
      iValue(1, 2) = iValue(2, 1) = 0.0;
    }

    // if requested, rotate the results in the global coordinate system
    if(bGlobal)
    {
      const Matrix& RG = data.LCS.Orientation();
      noalias( aux33 ) = prod( trans( RG ), iValue );
      noalias( iValue ) = prod( aux33, RG );
    }
  }

  OPT_INTERPOLATE_RESULTS_TO_STANDARD_GAUSS_POINTS(rValues);

  return true;
}

// =====================================================================================
//
// Class ShellThinElement3D3N - Serialization
//
// =====================================================================================

void ShellThinElement3D3N::save(Serializer& rSerializer) const
{
  KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer,  Element );
  rSerializer.save("CTr", mpCoordinateTransformation);
  rSerializer.save("Sec", mSections);
  rSerializer.save("IntM", (int)mThisIntegrationMethod);
}

void ShellThinElement3D3N::load(Serializer& rSerializer)
{
  KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer,  Element );
  rSerializer.load("CTr", mpCoordinateTransformation);
  rSerializer.load("Sec", mSections);
  int temp;
  rSerializer.load("IntM", temp);
  mThisIntegrationMethod = (IntegrationMethod)temp;
}

}
