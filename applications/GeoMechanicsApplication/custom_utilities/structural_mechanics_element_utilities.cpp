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
//                   Vahid Galavi

// System includes

// External includes

// Project includes
#include "structural_mechanics_element_utilities.h"
#include "structural_mechanics_application_variables.h"
#include "utilities/math_utils.h"

namespace Kratos {

bool GeoStructuralMechanicsElementUtilities::ComputeLumpedMassMatrix(
    const Properties& rProperites,
    const ProcessInfo& rCurrentProcessInfo)
{
    // Giving the globally defined setting (through ProcessInfo) priority
    // over the locally defined one (through Properties)
    // This is needed for the explicit solver, which requires the lumped
    // mass matrix and specifies it's computation through the ProcessInfo
    if (rCurrentProcessInfo.Has(COMPUTE_LUMPED_MASS_MATRIX)) {
        return rCurrentProcessInfo[COMPUTE_LUMPED_MASS_MATRIX];
    } else if (rProperites.Has(COMPUTE_LUMPED_MASS_MATRIX)) {
        return rProperites[COMPUTE_LUMPED_MASS_MATRIX];
    }

    // The default for all elements in StructuralMechanics is
    // to use the consistent-mass-matrix, hence returning false here
    return false;
}

/***********************************************************************************/
/***********************************************************************************/

bool GeoStructuralMechanicsElementUtilities::HasRayleighDamping(
    const Properties& rProperites,
    const ProcessInfo& rCurrentProcessInfo)
{
    return (std::abs(GetRayleighAlpha(rProperites, rCurrentProcessInfo)) > 0.0 ||
            std::abs(GetRayleighBeta(rProperites, rCurrentProcessInfo)) > 0.0);
}

/***********************************************************************************/
/***********************************************************************************/

double GeoStructuralMechanicsElementUtilities::GetRayleighAlpha(
    const Properties& rProperites,
    const ProcessInfo& rCurrentProcessInfo)
{
    // giving the locally defined setting (through Properties) priority
    // over the globally defined one (through ProcessInfo)
    if (rProperites.Has(RAYLEIGH_ALPHA)) {
        return rProperites[RAYLEIGH_ALPHA];
    } else if (rCurrentProcessInfo.Has(RAYLEIGH_ALPHA)) {
        return rCurrentProcessInfo[RAYLEIGH_ALPHA];
    }

    return 0.0;
}

/***********************************************************************************/
/***********************************************************************************/

double GeoStructuralMechanicsElementUtilities::GetRayleighBeta(
    const Properties& rProperites,
    const ProcessInfo& rCurrentProcessInfo)
{
    // Giving the locally defined setting (through Properties) priority
    // over the globally defined one (through ProcessInfo)
    if (rProperites.Has(RAYLEIGH_BETA)) {
        return rProperites[RAYLEIGH_BETA];
    } else if (rCurrentProcessInfo.Has(RAYLEIGH_BETA)) {
        return rCurrentProcessInfo[RAYLEIGH_BETA];
    }

    return 0.0;
}

/***********************************************************************************/
/***********************************************************************************/

double GeoStructuralMechanicsElementUtilities::GetDensityForMassMatrixComputation(const Element& rElement)
{
    // Getting the properties of the element
    const auto& r_prop = rElement.GetProperties();

    // Getting the density of the element
    const double density = r_prop[DENSITY];

    // Getting the mass factor
    if (rElement.Has(MASS_FACTOR)) {
        return rElement.GetValue(MASS_FACTOR) * density;
    } else if (r_prop.Has(MASS_FACTOR)) {
        return r_prop.GetValue(MASS_FACTOR) * density;
    } else {
        return density;
    }
}

/***********************************************************************************/
/***********************************************************************************/

void GeoStructuralMechanicsElementUtilities::CalculateRayleighDampingMatrix(
    Element& rElement,
    Element::MatrixType& rDampingMatrix,
    const ProcessInfo& rCurrentProcessInfo,
    const std::size_t MatrixSize)
{
    KRATOS_TRY;
    // Rayleigh Damping Matrix: alpha*M + beta*K

    // 1.-Resizing if needed
    if (rDampingMatrix.size1() != MatrixSize || rDampingMatrix.size2() != MatrixSize) {
        rDampingMatrix.resize(MatrixSize, MatrixSize, false);
    }
    noalias(rDampingMatrix) = ZeroMatrix(MatrixSize, MatrixSize);

    // 2.-Calculate StiffnessMatrix (if needed):
    const double beta = GetRayleighBeta(rElement.GetProperties(), rCurrentProcessInfo);
    if (std::abs(beta) > 0.0) {
        Element::MatrixType stiffness_matrix;
        rElement.CalculateLeftHandSide(stiffness_matrix, rCurrentProcessInfo);
        noalias(rDampingMatrix) += beta  * stiffness_matrix;
    }

    // 3.-Calculate MassMatrix (if needed):
    const double alpha = GetRayleighAlpha(rElement.GetProperties(), rCurrentProcessInfo);
    if (std::abs(alpha) > 0.0) {
        Element::MatrixType mass_matrix;
        rElement.CalculateMassMatrix(mass_matrix, rCurrentProcessInfo);
        noalias(rDampingMatrix) += alpha * mass_matrix;
    }

    KRATOS_CATCH("CalculateRayleighDampingMatrix")
}

/***********************************************************************************/
/***********************************************************************************/

double GeoStructuralMechanicsElementUtilities::CalculateReferenceLength2D2N(const Element& rElement)
{
    KRATOS_TRY;

    const array_1d<double, 3> delta_pos =
        rElement.GetGeometry()[1].GetInitialPosition().Coordinates() -
        rElement.GetGeometry()[0].GetInitialPosition().Coordinates();

    return std::sqrt((delta_pos[0] * delta_pos[0]) +
                     (delta_pos[1] * delta_pos[1]));

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

double GeoStructuralMechanicsElementUtilities::CalculateCurrentLength2D2N(const Element& rElement)
{
    KRATOS_TRY;

    const array_1d<double, 3> delta_pos =
        rElement.GetGeometry()[1].GetInitialPosition().Coordinates() -
        rElement.GetGeometry()[0].GetInitialPosition().Coordinates() +
        rElement.GetGeometry()[1].FastGetSolutionStepValue(DISPLACEMENT) -
        rElement.GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT);

    const double l = std::sqrt((delta_pos[0] * delta_pos[0]) +
                               (delta_pos[1] * delta_pos[1]));

    KRATOS_ERROR_IF(l <= std::numeric_limits<double>::epsilon())
            << "Element #" << rElement.Id() << " has a current length of zero!" << std::endl;

    return l;

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

double GeoStructuralMechanicsElementUtilities::CalculateReferenceLength3D2N(const Element& rElement)
{
    KRATOS_TRY;

    const array_1d<double, 3> delta_pos =
        rElement.GetGeometry()[1].GetInitialPosition().Coordinates() -
        rElement.GetGeometry()[0].GetInitialPosition().Coordinates();

    return MathUtils<double>::Norm3(delta_pos);

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

double GeoStructuralMechanicsElementUtilities::CalculateCurrentLength3D2N(const Element& rElement)
{
    KRATOS_TRY;

    const array_1d<double, 3> delta_pos =
        rElement.GetGeometry()[1].GetInitialPosition().Coordinates() -
        rElement.GetGeometry()[0].GetInitialPosition().Coordinates() +
        rElement.GetGeometry()[1].FastGetSolutionStepValue(DISPLACEMENT) -
        rElement.GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT);

    const double l = MathUtils<double>::Norm3(delta_pos);

    KRATOS_ERROR_IF(l <= std::numeric_limits<double>::epsilon())
            << "Element #" << rElement.Id() << " has a current length of zero!" << std::endl;

    return l;

    KRATOS_CATCH("")
}

}  // namespace Kratos.


