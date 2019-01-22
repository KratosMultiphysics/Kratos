//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher
//

// System includes

// External includes

// Project includes
#include "structural_mechanics_element_utilities.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos {
namespace StructuralMechanicsElementUtilities {

bool ComputeLumpedMassMatrix(
    const Properties& rProperites,
    const ProcessInfo& rCurrentProcessInfo)
{
    // giving the globally defined setting (through ProcessInfo) priority
    // over the locally defined one (through Properties)
    // this is needed for the explicit solver, which requires the lumped
    // mass matrix and specifies it's computation through the ProcessInfo
    if (rCurrentProcessInfo.Has(COMPUTE_LUMPED_MASS_MATRIX)) {
        return rCurrentProcessInfo[COMPUTE_LUMPED_MASS_MATRIX];
    }
    else if (rProperites.Has(COMPUTE_LUMPED_MASS_MATRIX)) {
        return rProperites[COMPUTE_LUMPED_MASS_MATRIX];
    }

    // the default for all elements in StructuralMechanics is
    // to use the consistent-mass-matrix, hence returning false here
    return false;
}

bool HasRayleighDamping(
    const Properties& rProperites,
    const ProcessInfo& rCurrentProcessInfo)
{
    return (GetRayleighAlpha(rProperites, rCurrentProcessInfo) > 0.0 ||
            GetRayleighBeta(rProperites, rCurrentProcessInfo) > 0.0);
}

double GetRayleighAlpha(
    const Properties& rProperites,
    const ProcessInfo& rCurrentProcessInfo)
{
    // giving the locally defined setting (through Properties) priority
    // over the globally defined one (through ProcessInfo)
    if (rProperites.Has(RAYLEIGH_ALPHA)) {
        return rProperites[RAYLEIGH_ALPHA];
    }
    else if (rCurrentProcessInfo.Has(RAYLEIGH_ALPHA)) {
        return rCurrentProcessInfo[RAYLEIGH_ALPHA];
    }

    return 0.0;
}

double GetRayleighBeta(
    const Properties& rProperites,
    const ProcessInfo& rCurrentProcessInfo)
{
    // giving the locally defined setting (through Properties) priority
    // over the globally defined one (through ProcessInfo)
    if (rProperites.Has(RAYLEIGH_BETA)) {
        return rProperites[RAYLEIGH_BETA];
    }
    else if (rCurrentProcessInfo.Has(RAYLEIGH_BETA)) {
        return rCurrentProcessInfo[RAYLEIGH_BETA];
    }

    return 0.0;
}

void CalculateRayleighDampingMatrix(
    Element& rElement,
    Element::MatrixType& rDampingMatrix,
    /*const*/ProcessInfo& rCurrentProcessInfo,
    const std::size_t MatrixSize)
{
    KRATOS_TRY;

    // 1.-Resizing if needed
    if (rDampingMatrix.size1() != MatrixSize || rDampingMatrix.size2() != MatrixSize) {
        rDampingMatrix.resize( MatrixSize, MatrixSize, false );
    }
    noalias(rDampingMatrix) = ZeroMatrix(MatrixSize, MatrixSize);

    if (!HasRayleighDamping(rElement.GetProperties(), rCurrentProcessInfo)) {
        // Do nothing if no rayleigh-damping is specified
        return;
    }

    // 2.-Calculate StiffnessMatrix:
    Element::MatrixType stiffness_matrix;
    rElement.CalculateLeftHandSide(stiffness_matrix, rCurrentProcessInfo);

    // 3.-Calculate MassMatrix:
    Element::MatrixType mass_matrix;
    rElement.CalculateMassMatrix(mass_matrix, rCurrentProcessInfo);

    // 4.-Compose the Damping Matrix:
    // Rayleigh Damping Matrix: alpha*M + beta*K
    noalias(rDampingMatrix) += GetRayleighAlpha(rElement.GetProperties(), rCurrentProcessInfo) * mass_matrix;
    noalias(rDampingMatrix) += GetRayleighBeta(rElement.GetProperties(), rCurrentProcessInfo)  * stiffness_matrix;

    KRATOS_CATCH( "CalculateRayleighDampingMatrix" )
}


double CalculateCurrentLength2D2N(const Element& rElement) {
  KRATOS_TRY;
  const double numerical_limit = std::numeric_limits<double>::epsilon();
  const double du =
      rElement.GetGeometry()[1].FastGetSolutionStepValue(DISPLACEMENT_X) -
      rElement.GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT_X);
  const double dv =
      rElement.GetGeometry()[1].FastGetSolutionStepValue(DISPLACEMENT_Y) -
      rElement.GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT_Y);

  const double dx = rElement.GetGeometry()[1].X0() - rElement.GetGeometry()[0].X0();
  const double dy = rElement.GetGeometry()[1].Y0() - rElement.GetGeometry()[0].Y0();

  const double l = std::sqrt((du + dx) * (du + dx) + (dv + dy) * (dv + dy));

  KRATOS_ERROR_IF(l < numerical_limit)
   << "Current length of element " << rElement.Id() << " ~0" << std::endl;
  return l;
  KRATOS_CATCH("")
}

double CalculateReferenceLength2D2N(const Element& rElement) {
  KRATOS_TRY;
  const double numerical_limit = std::numeric_limits<double>::epsilon();
  const double dx = rElement.GetGeometry()[1].X0() - rElement.GetGeometry()[0].X0();
  const double dy = rElement.GetGeometry()[1].Y0() - rElement.GetGeometry()[0].Y0();
  const double L = std::sqrt((dx * dx) + (dy * dy));

  KRATOS_ERROR_IF(L < numerical_limit)
   << "Reference length of element " << rElement.Id() << " ~0" << std::endl;
  return L;
  KRATOS_CATCH("")
}


double CalculateReferenceLength3D2N(const Element& rElement) {

  KRATOS_TRY;
  const double numerical_limit = std::numeric_limits<double>::epsilon();
  const double dx = rElement.GetGeometry()[1].X0() - rElement.GetGeometry()[0].X0();
  const double dy = rElement.GetGeometry()[1].Y0() - rElement.GetGeometry()[0].Y0();
  const double dz = rElement.GetGeometry()[1].Z0() - rElement.GetGeometry()[0].Z0();
  const double L = std::sqrt(dx * dx + dy * dy + dz * dz);

  KRATOS_ERROR_IF(L<=numerical_limit)
   << "Reference length of element " << rElement.Id() << " ~0" << std::endl;
  return L;
  KRATOS_CATCH("")
}

double CalculateCurrentLength3D2N(const Element& rElement) {

  KRATOS_TRY;
  const double numerical_limit = std::numeric_limits<double>::epsilon();
  const double du =
      rElement.GetGeometry()[1].FastGetSolutionStepValue(DISPLACEMENT_X) -
      rElement.GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT_X);
  const double dv =
      rElement.GetGeometry()[1].FastGetSolutionStepValue(DISPLACEMENT_Y) -
      rElement.GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT_Y);
  const double dw =
      rElement.GetGeometry()[1].FastGetSolutionStepValue(DISPLACEMENT_Z) -
      rElement.GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT_Z);
  const double dx = rElement.GetGeometry()[1].X0() - rElement.GetGeometry()[0].X0();
  const double dy = rElement.GetGeometry()[1].Y0() - rElement.GetGeometry()[0].Y0();
  const double dz = rElement.GetGeometry()[1].Z0() - rElement.GetGeometry()[0].Z0();
  const double l = std::sqrt((du + dx) * (du + dx) + (dv + dy) * (dv + dy) +
                             (dw + dz) * (dw + dz));

  KRATOS_ERROR_IF(l<=numerical_limit)
   << "Current length of element " << rElement.Id() << " ~0" << std::endl;
  return l;
  KRATOS_CATCH("")
}


} // namespace StructuralMechanicsElementUtilities.
}  // namespace Kratos.


