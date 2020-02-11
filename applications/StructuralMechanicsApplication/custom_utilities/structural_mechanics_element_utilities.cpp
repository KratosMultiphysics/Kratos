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
//                   Vicente Mataix Ferrandiz
//                   Riccardo Rossi
//                   Ruben Zorrilla
//

// System includes

// External includes

// Project includes
#include "includes/checks.h"
#include "structural_mechanics_element_utilities.h"
#include "structural_mechanics_application_variables.h"
#include "utilities/math_utils.h"

namespace Kratos {
namespace StructuralMechanicsElementUtilities {

int SolidElementCheck(
    const Element& rElement,
    const ProcessInfo& rCurrentProcessInfo,
    std::vector<ConstitutiveLaw::Pointer>& rConstitutiveLaws
    )
{
    const auto& r_geometry = rElement.GetGeometry();
    const auto& r_properties = rElement.GetProperties();
    const SizeType number_of_nodes = r_geometry.size();
    const SizeType dimension = r_geometry.WorkingSpaceDimension();

    // Verify that the variables are correctly initialized
    KRATOS_CHECK_VARIABLE_KEY(DISPLACEMENT)
    KRATOS_CHECK_VARIABLE_KEY(VELOCITY)
    KRATOS_CHECK_VARIABLE_KEY(ACCELERATION)
    KRATOS_CHECK_VARIABLE_KEY(DENSITY)
    KRATOS_CHECK_VARIABLE_KEY(VOLUME_ACCELERATION)
    KRATOS_CHECK_VARIABLE_KEY(THICKNESS)

    // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
    for ( IndexType i = 0; i < number_of_nodes; i++ ) {
        const NodeType &rnode = r_geometry[i];
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DISPLACEMENT,rnode)

        KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_X, rnode)
        KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_Y, rnode)
        KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_Z, rnode)
    }

    // Verify that the constitutive law exists
    KRATOS_ERROR_IF_NOT(r_properties.Has( CONSTITUTIVE_LAW )) << "Constitutive law not provided for property " << r_properties.Id() << std::endl;

    // Verify that the constitutive law has the correct dimension
    const SizeType strain_size = r_properties.GetValue( CONSTITUTIVE_LAW )->GetStrainSize();
    if ( dimension == 2 ) {
        KRATOS_ERROR_IF( strain_size < 3 || strain_size > 4) << "Wrong constitutive law used. This is a 2D element! expected strain size is 3 or 4 (el id = ) " << rElement.Id() << std::endl;
    } else {
        KRATOS_ERROR_IF_NOT(strain_size == 6) << "Wrong constitutive law used. This is a 3D element! expected strain size is 6 (el id = ) "<<  rElement.Id() << std::endl;
    }
    
    // Check constitutive law
    if ( rConstitutiveLaws.size() > 0 ) {
        return rConstitutiveLaws[0]->Check( r_properties, r_geometry, rCurrentProcessInfo );
    }

    return 0;
}

/***********************************************************************************/
/***********************************************************************************/

array_1d<double, 3> GetBodyForce(
    const Element& rElement,
    const GeometryType::IntegrationPointsArrayType& rIntegrationPoints,
    const IndexType PointNumber
    )
{
    array_1d<double, 3> body_force;
    for (IndexType i = 0; i < 3; ++i)
        body_force[i] = 0.0;

    const auto& r_properties = rElement.GetProperties();
    double density = 0.0;
    if (r_properties.Has( DENSITY ))
        density = r_properties[DENSITY];

    if (r_properties.Has( VOLUME_ACCELERATION ))
        noalias(body_force) += density * r_properties[VOLUME_ACCELERATION];

    const auto& r_geometry = rElement.GetGeometry();
    if( r_geometry[0].SolutionStepsDataHas(VOLUME_ACCELERATION) ) {
        Vector N(r_geometry.size());
        N = r_geometry.ShapeFunctionsValues(N, rIntegrationPoints[PointNumber].Coordinates());
        for (IndexType i_node = 0; i_node < r_geometry.size(); ++i_node)
            noalias(body_force) += N[i_node] * density * r_geometry[i_node].FastGetSolutionStepValue(VOLUME_ACCELERATION);
    }

    return body_force;
}

/***********************************************************************************/
/***********************************************************************************/

bool ComputeLumpedMassMatrix(
    const Properties& rProperties,
    const ProcessInfo& rCurrentProcessInfo)
{
    // Giving the globally defined setting (through ProcessInfo) priority
    // over the locally defined one (through Properties)
    // This is needed for the explicit solver, which requires the lumped
    // mass matrix and specifies it's computation through the ProcessInfo
    if (rCurrentProcessInfo.Has(COMPUTE_LUMPED_MASS_MATRIX)) {
        return rCurrentProcessInfo[COMPUTE_LUMPED_MASS_MATRIX];
    } else if (rProperties.Has(COMPUTE_LUMPED_MASS_MATRIX)) {
        return rProperties[COMPUTE_LUMPED_MASS_MATRIX];
    }

    // The default for all elements in StructuralMechanics is
    // to use the consistent-mass-matrix, hence returning false here
    return false;
}

/***********************************************************************************/
/***********************************************************************************/

bool HasRayleighDamping(
    const Properties& rProperties,
    const ProcessInfo& rCurrentProcessInfo)
{
    return (std::abs(GetRayleighAlpha(rProperties, rCurrentProcessInfo)) > 0.0 ||
            std::abs(GetRayleighBeta(rProperties, rCurrentProcessInfo)) > 0.0);
}

/***********************************************************************************/
/***********************************************************************************/

double GetRayleighAlpha(
    const Properties& rProperties,
    const ProcessInfo& rCurrentProcessInfo)
{
    // giving the locally defined setting (through Properties) priority
    // over the globally defined one (through ProcessInfo)
    if (rProperties.Has(RAYLEIGH_ALPHA)) {
        return rProperties[RAYLEIGH_ALPHA];
    } else if (rCurrentProcessInfo.Has(RAYLEIGH_ALPHA)) {
        return rCurrentProcessInfo[RAYLEIGH_ALPHA];
    }

    return 0.0;
}

/***********************************************************************************/
/***********************************************************************************/

double GetRayleighBeta(
    const Properties& rProperties,
    const ProcessInfo& rCurrentProcessInfo)
{
    // Giving the locally defined setting (through Properties) priority
    // over the globally defined one (through ProcessInfo)
    if (rProperties.Has(RAYLEIGH_BETA)) {
        return rProperties[RAYLEIGH_BETA];
    } else if (rCurrentProcessInfo.Has(RAYLEIGH_BETA)) {
        return rCurrentProcessInfo[RAYLEIGH_BETA];
    }

    return 0.0;
}

/***********************************************************************************/
/***********************************************************************************/

double GetDensityForMassMatrixComputation(const Element& rElement)
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

void CalculateRayleighDampingMatrix(
    Element& rElement,
    Element::MatrixType& rDampingMatrix,
    /*const*/ProcessInfo& rCurrentProcessInfo,
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

double CalculateReferenceLength2D2N(const Element& rElement)
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

double CalculateCurrentLength2D2N(const Element& rElement)
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

double CalculateReferenceLength3D2N(const Element& rElement)
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

double CalculateCurrentLength3D2N(const Element& rElement)
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

} // namespace StructuralMechanicsElementUtilities.
}  // namespace Kratos.


