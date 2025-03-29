// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Philipp Bucher
//                   Vicente Mataix Ferrandiz
//                   Riccardo Rossi
//                   Ruben Zorrilla
//                   Alejandro Cornejo
//

// System includes

// External includes

// Project includes
#include "includes/checks.h"
#include "structural_mechanics_element_utilities.h"
#include "structural_mechanics_application_variables.h"
#include "custom_utilities/constitutive_law_utilities.h"
#include "utilities/math_utils.h"

namespace Kratos {
namespace StructuralMechanicsElementUtilities {

int SolidElementCheck(
    const Element& rElement,
    const ProcessInfo& rCurrentProcessInfo,
    const std::vector<ConstitutiveLaw::Pointer>& rConstitutiveLaws
    )
{
    const auto& r_geometry = rElement.GetGeometry();
    const auto& r_properties = rElement.GetProperties();
    const SizeType number_of_nodes = r_geometry.size();
    const SizeType dimension = r_geometry.WorkingSpaceDimension();

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
    const ProcessInfo& rCurrentProcessInfo,
    const std::size_t MatrixSize)
{
    KRATOS_TRY;
    // Rayleigh Damping Matrix: alpha*M + beta*K

    const double alpha = GetRayleighAlpha(rElement.GetProperties(), rCurrentProcessInfo);
    const double beta = GetRayleighBeta(rElement.GetProperties(), rCurrentProcessInfo);

    if (std::abs(alpha) < 1E-12 && std::abs(beta) < 1E-12) {
        // no damping specified, only setting the matrix to zero
        if (rDampingMatrix.size1() != MatrixSize || rDampingMatrix.size2() != MatrixSize) {
            rDampingMatrix.resize(MatrixSize, MatrixSize, false);
        }
        noalias(rDampingMatrix) = ZeroMatrix(MatrixSize, MatrixSize);
    } else if (std::abs(alpha) > 1E-12 && std::abs(beta) < 1E-12) {
        // damping only required with the mass matrix
        rElement.CalculateMassMatrix(rDampingMatrix, rCurrentProcessInfo); // pass damping matrix to avoid creating a temporary
        rDampingMatrix *= alpha;
    } else if (std::abs(alpha) < 1E-12 && std::abs(beta) > 1E-12) {
        // damping only required with the stiffness matrix
        rElement.CalculateLeftHandSide(rDampingMatrix, rCurrentProcessInfo); // pass damping matrix to avoid creating a temporary
        rDampingMatrix *= beta;
    } else {
        // damping with both mass matrix and stiffness matrix required
        rElement.CalculateLeftHandSide(rDampingMatrix, rCurrentProcessInfo); // pass damping matrix to avoid creating a temporary
        rDampingMatrix *= beta;

        Matrix mass_matrix;
        rElement.CalculateMassMatrix(mass_matrix, rCurrentProcessInfo);
        noalias(rDampingMatrix) += alpha  * mass_matrix;
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

/***********************************************************************************/
/***********************************************************************************/

void InitialCheckLocalAxes(
    const array_1d<double, 3>& rv1,
    const array_1d<double, 3>& rv2,
    const array_1d<double, 3>& rv3,
    const double Tolerance
    )
{
    if (MathUtils<double>::Norm3(rv1) > 1.0 + Tolerance ||
        MathUtils<double>::Norm3(rv2) > 1.0 + Tolerance ||
        MathUtils<double>::Norm3(rv3) > 1.0 + Tolerance) {
            KRATOS_ERROR << "The norm of one of the LOCAL_AXIS is greater than 1.0!" << std::endl;
    }
}

/***********************************************************************************/
/***********************************************************************************/

void BuildRotationMatrix(
    BoundedMatrix<double, 3, 3>& rRotationMatrix,
    const array_1d<double, 3>& rv1,
    const array_1d<double, 3>& rv2,
    const array_1d<double, 3>& rv3
    )
{
    rRotationMatrix(0, 0) = rv1[0]; rRotationMatrix(0, 1) = rv1[1]; rRotationMatrix(0, 2) = rv1[2];
    rRotationMatrix(1, 0) = rv2[0]; rRotationMatrix(1, 1) = rv2[1]; rRotationMatrix(1, 2) = rv2[2];
    rRotationMatrix(2, 0) = rv3[0]; rRotationMatrix(2, 1) = rv3[1]; rRotationMatrix(2, 2) = rv3[2];
}

/***********************************************************************************/
/***********************************************************************************/

void BuildRotationMatrixForBeam(
    BoundedMatrix<double, 3, 3>& rRotationMatrix,
    const double AlphaAngle
)
{
    rRotationMatrix.clear();
    const double s = std::sin(AlphaAngle);
    const double c = std::cos(AlphaAngle);
    rRotationMatrix(0, 0) = c;
    rRotationMatrix(0, 1) = -s;
    rRotationMatrix(1, 0) = s;
    rRotationMatrix(1, 1) = c;
    rRotationMatrix(2, 2) = 1.0;
}

/***********************************************************************************/
/***********************************************************************************/

void BuildRotationMatrixForTruss(
    BoundedMatrix<double, 2, 2>& rRotationMatrix,
    const double AlphaAngle
)
{
    rRotationMatrix.clear();
    const double s = std::sin(AlphaAngle);
    const double c = std::cos(AlphaAngle);
    rRotationMatrix(0, 0) = c;
    rRotationMatrix(0, 1) = -s;
    rRotationMatrix(1, 0) = s;
    rRotationMatrix(1, 1) = c;
}

/***********************************************************************************/
/***********************************************************************************/

void BuildElementSizeRotationMatrixFor2D2NTruss(
    const BoundedMatrix<double, 2, 2>& rRotationMatrix,
    BoundedMatrix<double, 4, 4>& rElementSizeRotationMatrix
)
{
    rElementSizeRotationMatrix.clear();

    rElementSizeRotationMatrix(0, 0) = rRotationMatrix(0, 0);
    rElementSizeRotationMatrix(0, 1) = rRotationMatrix(0, 1);
    rElementSizeRotationMatrix(1, 0) = rRotationMatrix(1, 0);
    rElementSizeRotationMatrix(1, 1) = rRotationMatrix(1, 1);


    rElementSizeRotationMatrix(2, 2) = rRotationMatrix(0, 0);
    rElementSizeRotationMatrix(2, 3) = rRotationMatrix(0, 1);
    rElementSizeRotationMatrix(3, 2) = rRotationMatrix(1, 0);
    rElementSizeRotationMatrix(3, 3) = rRotationMatrix(1, 1);
}

/***********************************************************************************/
/***********************************************************************************/

void BuildElementSizeRotationMatrixFor2D3NTruss(
    const BoundedMatrix<double, 2, 2>& rRotationMatrix,
    BoundedMatrix<double, 6, 6>& rElementSizeRotationMatrix
)
{
    rElementSizeRotationMatrix.clear();

    rElementSizeRotationMatrix(0, 0) = rRotationMatrix(0, 0);
    rElementSizeRotationMatrix(0, 1) = rRotationMatrix(0, 1);
    rElementSizeRotationMatrix(1, 0) = rRotationMatrix(1, 0);
    rElementSizeRotationMatrix(1, 1) = rRotationMatrix(1, 1);

    rElementSizeRotationMatrix(2, 2) = rRotationMatrix(0, 0);
    rElementSizeRotationMatrix(2, 3) = rRotationMatrix(0, 1);
    rElementSizeRotationMatrix(3, 2) = rRotationMatrix(1, 0);
    rElementSizeRotationMatrix(3, 3) = rRotationMatrix(1, 1);

    rElementSizeRotationMatrix(4, 4) = rRotationMatrix(0, 0);
    rElementSizeRotationMatrix(4, 5) = rRotationMatrix(0, 1);
    rElementSizeRotationMatrix(5, 4) = rRotationMatrix(1, 0);
    rElementSizeRotationMatrix(5, 5) = rRotationMatrix(1, 1);
}

/***********************************************************************************/
/***********************************************************************************/

void BuildElementSizeRotationMatrixFor3D2NTruss(
    const BoundedMatrix<double, 3, 3>& rRotationMatrix,
    BoundedMatrix<double, 6, 6>& rElementSizeRotationMatrix)
{
    rElementSizeRotationMatrix.clear();

    rElementSizeRotationMatrix(0, 0) = rRotationMatrix(0, 0);
    rElementSizeRotationMatrix(0, 1) = rRotationMatrix(0, 1);
    rElementSizeRotationMatrix(0, 2) = rRotationMatrix(0, 2);

    rElementSizeRotationMatrix(1, 0) = rRotationMatrix(1, 0);
    rElementSizeRotationMatrix(1, 1) = rRotationMatrix(1, 1);
    rElementSizeRotationMatrix(1, 2) = rRotationMatrix(1, 2);

    rElementSizeRotationMatrix(2, 0) = rRotationMatrix(2, 0);
    rElementSizeRotationMatrix(2, 1) = rRotationMatrix(2, 1);
    rElementSizeRotationMatrix(2, 2) = rRotationMatrix(2, 2);

    rElementSizeRotationMatrix(3, 3) = rRotationMatrix(0, 0);
    rElementSizeRotationMatrix(3, 4) = rRotationMatrix(0, 1);
    rElementSizeRotationMatrix(3, 5) = rRotationMatrix(0, 2);

    rElementSizeRotationMatrix(4, 3) = rRotationMatrix(1, 0);
    rElementSizeRotationMatrix(4, 4) = rRotationMatrix(1, 1);
    rElementSizeRotationMatrix(4, 5) = rRotationMatrix(1, 2);

    rElementSizeRotationMatrix(5, 3) = rRotationMatrix(2, 0);
    rElementSizeRotationMatrix(5, 4) = rRotationMatrix(2, 1);
    rElementSizeRotationMatrix(5, 5) = rRotationMatrix(2, 2);
}

/***********************************************************************************/
/***********************************************************************************/

void BuildElementSizeRotationMatrixFor3D3NTruss(
    const BoundedMatrix<double, 3, 3>& rRotationMatrix,
    BoundedMatrix<double, 9, 9>& rElementSizeRotationMatrix)
{
    rElementSizeRotationMatrix.clear();

    rElementSizeRotationMatrix(0, 0) = rRotationMatrix(0, 0);
    rElementSizeRotationMatrix(0, 1) = rRotationMatrix(0, 1);
    rElementSizeRotationMatrix(0, 2) = rRotationMatrix(0, 2);

    rElementSizeRotationMatrix(1, 0) = rRotationMatrix(1, 0);
    rElementSizeRotationMatrix(1, 1) = rRotationMatrix(1, 1);
    rElementSizeRotationMatrix(1, 2) = rRotationMatrix(1, 2);

    rElementSizeRotationMatrix(2, 0) = rRotationMatrix(2, 0);
    rElementSizeRotationMatrix(2, 1) = rRotationMatrix(2, 1);
    rElementSizeRotationMatrix(2, 2) = rRotationMatrix(2, 2);

    rElementSizeRotationMatrix(3, 3) = rRotationMatrix(0, 0);
    rElementSizeRotationMatrix(3, 4) = rRotationMatrix(0, 1);
    rElementSizeRotationMatrix(3, 5) = rRotationMatrix(0, 2);

    rElementSizeRotationMatrix(4, 3) = rRotationMatrix(1, 0);
    rElementSizeRotationMatrix(4, 4) = rRotationMatrix(1, 1);
    rElementSizeRotationMatrix(4, 5) = rRotationMatrix(1, 2);

    rElementSizeRotationMatrix(5, 3) = rRotationMatrix(2, 0);
    rElementSizeRotationMatrix(5, 4) = rRotationMatrix(2, 1);
    rElementSizeRotationMatrix(5, 5) = rRotationMatrix(2, 2);

    rElementSizeRotationMatrix(6, 6) = rRotationMatrix(0, 0);
    rElementSizeRotationMatrix(6, 7) = rRotationMatrix(0, 1);
    rElementSizeRotationMatrix(6, 8) = rRotationMatrix(0, 2);

    rElementSizeRotationMatrix(7, 6) = rRotationMatrix(1, 0);
    rElementSizeRotationMatrix(7, 7) = rRotationMatrix(1, 1);
    rElementSizeRotationMatrix(7, 8) = rRotationMatrix(1, 2);

    rElementSizeRotationMatrix(8, 6) = rRotationMatrix(2, 0);
    rElementSizeRotationMatrix(8, 7) = rRotationMatrix(2, 1);
    rElementSizeRotationMatrix(8, 8) = rRotationMatrix(2, 2);
}

/***********************************************************************************/
/***********************************************************************************/

double GetReferenceRotationAngle2D2NBeam(const GeometryType& rGeometry)
{
    const auto &r_node_1 = rGeometry[0];
    const auto &r_node_2 = rGeometry[1];

    const double delta_x = r_node_2.X0() - r_node_1.X0();
    const double delta_y = r_node_2.Y0() - r_node_1.Y0();

    return std::atan2(delta_y, delta_x);
}

/***********************************************************************************/
/***********************************************************************************/

double GetReferenceRotationAngle2D3NBeam(const GeometryType& rGeometry)
{
    const auto &r_node_1 = rGeometry[0];
    const auto &r_node_2 = rGeometry[2];

    const double delta_x = r_node_2.X0() - r_node_1.X0();
    const double delta_y = r_node_2.Y0() - r_node_1.Y0();

    return std::atan2(delta_y, delta_x);
}

/***********************************************************************************/
/***********************************************************************************/

void BuildElementSizeRotationMatrixFor2D2NBeam(
    const BoundedMatrix<double, 3, 3>& rRotationMatrix,
    BoundedMatrix<double, 6, 6>& rElementSizeRotationMatrix
    )
{
    rElementSizeRotationMatrix.clear();

    rElementSizeRotationMatrix(0, 0) = rRotationMatrix(0, 0);
    rElementSizeRotationMatrix(0, 1) = rRotationMatrix(0, 1);
    rElementSizeRotationMatrix(0, 2) = rRotationMatrix(0, 2);

    rElementSizeRotationMatrix(1, 0) = rRotationMatrix(1, 0);
    rElementSizeRotationMatrix(1, 1) = rRotationMatrix(1, 1);
    rElementSizeRotationMatrix(1, 2) = rRotationMatrix(1, 2);

    rElementSizeRotationMatrix(2, 0) = rRotationMatrix(2, 0);
    rElementSizeRotationMatrix(2, 1) = rRotationMatrix(2, 1);
    rElementSizeRotationMatrix(2, 2) = rRotationMatrix(2, 2);

    rElementSizeRotationMatrix(3, 3) = rRotationMatrix(0, 0);
    rElementSizeRotationMatrix(3, 4) = rRotationMatrix(0, 1);
    rElementSizeRotationMatrix(3, 5) = rRotationMatrix(0, 2);

    rElementSizeRotationMatrix(4, 3) = rRotationMatrix(1, 0);
    rElementSizeRotationMatrix(4, 4) = rRotationMatrix(1, 1);
    rElementSizeRotationMatrix(4, 5) = rRotationMatrix(1, 2);

    rElementSizeRotationMatrix(5, 3) = rRotationMatrix(2, 0);
    rElementSizeRotationMatrix(5, 4) = rRotationMatrix(2, 1);
    rElementSizeRotationMatrix(5, 5) = rRotationMatrix(2, 2);
}

/***********************************************************************************/
/***********************************************************************************/

void BuildElementSizeRotationMatrixFor2D3NBeam(
    const BoundedMatrix<double, 3, 3>& rRotationMatrix,
    BoundedMatrix<double, 9, 9>& rElementSizeRotationMatrix
    )
{
    rElementSizeRotationMatrix.clear();

    rElementSizeRotationMatrix(0, 0) = rRotationMatrix(0, 0);
    rElementSizeRotationMatrix(0, 1) = rRotationMatrix(0, 1);
    rElementSizeRotationMatrix(0, 2) = rRotationMatrix(0, 2);

    rElementSizeRotationMatrix(1, 0) = rRotationMatrix(1, 0);
    rElementSizeRotationMatrix(1, 1) = rRotationMatrix(1, 1);
    rElementSizeRotationMatrix(1, 2) = rRotationMatrix(1, 2);

    rElementSizeRotationMatrix(2, 0) = rRotationMatrix(2, 0);
    rElementSizeRotationMatrix(2, 1) = rRotationMatrix(2, 1);
    rElementSizeRotationMatrix(2, 2) = rRotationMatrix(2, 2);

    rElementSizeRotationMatrix(3, 3) = rRotationMatrix(0, 0);
    rElementSizeRotationMatrix(3, 4) = rRotationMatrix(0, 1);
    rElementSizeRotationMatrix(3, 5) = rRotationMatrix(0, 2);

    rElementSizeRotationMatrix(4, 3) = rRotationMatrix(1, 0);
    rElementSizeRotationMatrix(4, 4) = rRotationMatrix(1, 1);
    rElementSizeRotationMatrix(4, 5) = rRotationMatrix(1, 2);

    rElementSizeRotationMatrix(5, 3) = rRotationMatrix(2, 0);
    rElementSizeRotationMatrix(5, 4) = rRotationMatrix(2, 1);
    rElementSizeRotationMatrix(5, 5) = rRotationMatrix(2, 2);

    rElementSizeRotationMatrix(6, 6) = rRotationMatrix(0, 0);
    rElementSizeRotationMatrix(6, 7) = rRotationMatrix(0, 1);
    rElementSizeRotationMatrix(6, 8) = rRotationMatrix(0, 2);

    rElementSizeRotationMatrix(7, 6) = rRotationMatrix(1, 0);
    rElementSizeRotationMatrix(7, 7) = rRotationMatrix(1, 1);
    rElementSizeRotationMatrix(7, 8) = rRotationMatrix(1, 2);

    rElementSizeRotationMatrix(8, 6) = rRotationMatrix(2, 0);
    rElementSizeRotationMatrix(8, 7) = rRotationMatrix(2, 1);
    rElementSizeRotationMatrix(8, 8) = rRotationMatrix(2, 2);
}

/***********************************************************************************/
/***********************************************************************************/

double CalculatePhi(const Properties& rProperties, const double L)
{
    const double E   = rProperties[YOUNG_MODULUS];
    const double I   = rProperties[I33];
    const double A_s = rProperties[AREA_EFFECTIVE_Y];
    const double G   = ConstitutiveLawUtilities<3>::CalculateShearModulus(rProperties);

    if (A_s == 0.0) // If effective area is null -> Euler Bernouilli case
        return 0.0;
    else
        return 12.0 * E * I / (G * A_s * std::pow(L, 2));
}

/***********************************************************************************/
/***********************************************************************************/

double CalculatePhiY(const Properties& rProperties, const double L)
{
    const double E   = rProperties[YOUNG_MODULUS];
    const double I   = rProperties[I22];
    const double A_s = rProperties[AREA_EFFECTIVE_Z];
    const double G   = ConstitutiveLawUtilities<3>::CalculateShearModulus(rProperties);

    if (A_s == 0.0) // If effective area is null -> Euler Bernouilli case
        return 0.0;
    else
        return 12.0 * E * I / (G * A_s * std::pow(L, 2));
}

/***********************************************************************************/
/***********************************************************************************/

void InitializeConstitutiveLawValuesForStressCalculation(ConstitutiveLaw::Parameters& rValues,
    Vector& rStrainVector, Vector& rStressVector, Matrix& rConstitutiveMatrix)
{
    InitializeConstitutiveLawValuesForStressCalculation(rValues, rStrainVector, rStressVector);
    rValues.SetConstitutiveMatrix(rConstitutiveMatrix);
}

/***********************************************************************************/
/***********************************************************************************/

void InitializeConstitutiveLawValuesForStressCalculation(ConstitutiveLaw::Parameters& rValues,
    Vector& rStrainVector, Vector& rStressVector)
{
    rValues.GetOptions().Set(ConstitutiveLaw::COMPUTE_STRESS             , true);
    rValues.GetOptions().Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);

    rStrainVector.clear();
    rValues.SetStrainVector(rStrainVector);
    rValues.SetStressVector(rStressVector);
}

/***********************************************************************************/
/***********************************************************************************/

BoundedMatrix<double, 3, 3> GetFrenetSerretMatrix3D(const GeometryType& rGeometry)
{
    BoundedMatrix<double, 3, 3> T;
    T.clear(); // global to local

    array_1d<double, 3> t;
    array_1d<double, 3> n;
    array_1d<double, 3> m;

    // t is the axis of the truss
    noalias(t) = rGeometry[1].GetInitialPosition() - rGeometry[0].GetInitialPosition();
    t /= norm_2(t);

    n.clear();

    if (rGeometry.Has(LOCAL_AXIS_2)) {
        noalias(n) = rGeometry.GetValue(LOCAL_AXIS_2);
    } else {
        // Default
        n[1] = 1.0;
    }

    if (norm_2(t-n) <= 1.0e-8) { // colineal, hence we use another aux vector
        n.clear();
        n[2] = 1.0;
    }

    // Gram-Schmidt ortogonalization
    n = n - inner_prod(t, n) / inner_prod(t, t) * t;
    n /= norm_2(n);

    noalias(m) = MathUtils<double>::CrossProduct(t, n);

    T(0, 0) = t[0];
    T(0, 1) = t[1];
    T(0, 2) = t[2];

    T(1, 0) = n[0];
    T(1, 1) = n[1];
    T(1, 2) = n[2];

    T(2, 0) = m[0];
    T(2, 1) = m[1];
    T(2, 2) = m[2];

    return T;
}

/***********************************************************************************/
/***********************************************************************************/

} // namespace StructuralMechanicsElementUtilities.
}  // namespace Kratos.


