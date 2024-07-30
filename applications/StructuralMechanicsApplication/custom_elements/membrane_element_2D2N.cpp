// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Ruben Zorrilla
//

// System includes

// External includes

// Project includes
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h"
#include "utilities/atomic_utilities.h"

// Application includes
#include "custom_elements/membrane_element_2D2N.h"
#include "structural_mechanics_application_variables.h"
#include "custom_utilities/structural_mechanics_element_utilities.h"
#include "custom_utilities/constitutive_law_utilities.h"

namespace Kratos
{

/***********************************************************************************/
/***********************************************************************************/

Element::Pointer MembraneElement2D2N::Create(
    IndexType NewId,
    NodesArrayType const &ThisNodes,
    PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<MembraneElement2D2N>(NewId, GetGeometry().Create(ThisNodes), pProperties);
}

/***********************************************************************************/
/***********************************************************************************/

Element::Pointer MembraneElement2D2N::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<MembraneElement2D2N>(NewId, pGeom, pProperties);
}

/***********************************************************************************/
/***********************************************************************************/

Element::Pointer MembraneElement2D2N::Clone(
    IndexType NewId,
    NodesArrayType const& rThisNodes) const
{
    KRATOS_TRY

    MembraneElement2D2N::Pointer p_new_elem = Kratos::make_intrusive<MembraneElement2D2N>(NewId, GetGeometry().Create(rThisNodes), pGetProperties());
    p_new_elem->SetData(this->GetData());
    p_new_elem->Set(Flags(*this));

    return p_new_elem;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void MembraneElement2D2N::EquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY;

    if (rResult.size() != mLocalSize) {
        rResult.resize(mLocalSize, false);
    }

    const auto& r_geom = GetGeometry();
    rResult[0] = r_geom[0].GetDof(DISPLACEMENT_X).EquationId();
    rResult[1] = r_geom[0].GetDof(DISPLACEMENT_Y).EquationId();
    rResult[2] = r_geom[1].GetDof(DISPLACEMENT_X).EquationId();
    rResult[3] = r_geom[1].GetDof(DISPLACEMENT_Y).EquationId();

    KRATOS_CATCH("")
};

/***********************************************************************************/
/***********************************************************************************/

void MembraneElement2D2N::GetDofList(
    DofsVectorType& rElementalDofList,
    const ProcessInfo& rCurrentProcessInfo
    ) const
{
    KRATOS_TRY;

    if (rElementalDofList.size() != mLocalSize) {
        rElementalDofList.resize(mLocalSize);
    }

    const auto& r_geom = GetGeometry();
    rElementalDofList[0] = r_geom[0].pGetDof(DISPLACEMENT_X);
    rElementalDofList[1] = r_geom[0].pGetDof(DISPLACEMENT_Y);
    rElementalDofList[2] = r_geom[1].pGetDof(DISPLACEMENT_X);
    rElementalDofList[3] = r_geom[1].pGetDof(DISPLACEMENT_Y);

    KRATOS_CATCH("")
};

/***********************************************************************************/
/***********************************************************************************/

void MembraneElement2D2N::GetFirstDerivativesVector(
    Vector& rValues,
    int Step) const
{
    KRATOS_TRY;

    if (rValues.size() != mLocalSize) {
        rValues.resize(mLocalSize, false);
    }

    const auto& r_geom = GetGeometry();
    rValues[0] = r_geom[0].FastGetSolutionStepValue(VELOCITY_X, Step);
    rValues[1] = r_geom[0].FastGetSolutionStepValue(VELOCITY_Y, Step);
    rValues[2] = r_geom[1].FastGetSolutionStepValue(VELOCITY_X, Step);
    rValues[3] = r_geom[1].FastGetSolutionStepValue(VELOCITY_Y, Step);

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void MembraneElement2D2N::GetSecondDerivativesVector(
    Vector& rValues,
    int Step) const
{
    KRATOS_TRY;

    if (rValues.size() != mLocalSize) {
        rValues.resize(mLocalSize, false);
    }

    const auto& r_geom = GetGeometry();
    rValues[0] = r_geom[0].FastGetSolutionStepValue(ACCELERATION_X, Step);
    rValues[1] = r_geom[0].FastGetSolutionStepValue(ACCELERATION_Y, Step);
    rValues[2] = r_geom[1].FastGetSolutionStepValue(ACCELERATION_X, Step);
    rValues[3] = r_geom[1].FastGetSolutionStepValue(ACCELERATION_Y, Step);

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void MembraneElement2D2N::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    // Get element data
    const auto& r_geom = GetGeometry();
    const auto& r_prop = GetProperties();

    // Get membrane data
    const double h = r_prop.GetValue(THICKNESS);
    const double young = r_prop.GetValue(YOUNG_MODULUS);

    // Calculate element lengths
    const double l = StructuralMechanicsElementUtilities::CalculateCurrentLength2D2N(*this);
    const double l0 = StructuralMechanicsElementUtilities::CalculateReferenceLength2D2N(*this);

    // Calculate material response
    const double l0_pow = std::pow(l0, 2);
    const double E = 0.5 * (std::pow(l, 2) - l0_pow) / l0_pow;
    const double S0 = GetMembranePrestress();
    const double S = young * E + S0;
    const bool is_compressed = S < 0.0 && std::abs(E) > 1.0e-12 ? true : false;

    // Check RHS size
    if (rRightHandSideVector.size() != mLocalSize) {
        rRightHandSideVector.resize(mLocalSize, false);
    }

    // Check LHS size
    if (rLeftHandSideMatrix.size1() != mLocalSize || rLeftHandSideMatrix.size2() != mLocalSize) {
        rLeftHandSideMatrix.resize(mLocalSize, mLocalSize, false);
    }

    // Set RHS and LHS to zero
    rLeftHandSideMatrix.clear();
    rRightHandSideVector.clear();

    // Calculate required kinematics
    array_1d<double, mLocalSize> B;
    const auto& r_x_0 = r_geom[0].Coordinates();
    const auto& r_x_1 = r_geom[1].Coordinates();
    const double dx = r_x_1[0] - r_x_0[0];
    const double dy = r_x_1[1] - r_x_0[1];
    B[0] = -dx / l0_pow;
    B[1] = -dy / l0_pow;
    B[2] = dx / l0_pow;
    B[3] = dy / l0_pow;

    // Calculate the internal force contribution
    array_1d<double, mLocalSize> f_int;
    if (is_compressed) {
        f_int.clear();
    } else {
        const double aux_f_int = h * l0 * S;
        noalias(f_int) = aux_f_int * B;
    }

    // Calculate the body force contribution
    array_1d<double, 3> f_ext;
    f_ext = StructuralMechanicsElementUtilities::GetBodyForce(*this, r_geom.IntegrationPoints(GeometryData::IntegrationMethod::GI_GAUSS_1), 0);
    f_ext *= h * l0 * 0.5; // Half the force to each node

    // Assemble the right hand side contributions
    rRightHandSideVector[0] += f_ext[0] - f_int[0];
    rRightHandSideVector[1] += f_ext[1] - f_int[1];
    rRightHandSideVector[2] += f_ext[0] - f_int[2];
    rRightHandSideVector[3] += f_ext[1] - f_int[3];

    // Calculate material stiffness matrix
    BoundedMatrix<double, mLocalSize, mLocalSize> K_mat;
    if (is_compressed) {
        K_mat.clear();
    } else {
        const double aux_K_mat = h * l0 * young;
        noalias(K_mat) = aux_K_mat * outer_prod(B, B);
    }

    // Calculate geometric stiffness matrix
    BoundedMatrix<double, mLocalSize, mLocalSize> K_geo;
    if (is_compressed) {
        K_geo.clear();
    } else {
        const double aux_K_geo = h * S / l0;

        K_geo.clear();
        for (IndexType i = 0; i < mNumNodes; ++i) {
            for (IndexType j = 0; j < mNumNodes; ++j) {
                for (IndexType d = 0; d < mDim; ++d) {
                    K_geo(i*mDim + d, j*mDim + d) = (i == j) ? aux_K_geo : - aux_K_geo;
                }
            }
        }
    }

    // Assemble the left hand side contributions
    noalias(rLeftHandSideMatrix) += K_mat + K_geo;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void MembraneElement2D2N::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    // Get element data
    const auto& r_geom = GetGeometry();
    const auto& r_prop = GetProperties();

    // Get membrane data
    const double h = r_prop.GetValue(THICKNESS);
    const double young = r_prop.GetValue(YOUNG_MODULUS);

    // Calculate element lengths
    const double l = StructuralMechanicsElementUtilities::CalculateCurrentLength2D2N(*this);
    const double l0 = StructuralMechanicsElementUtilities::CalculateReferenceLength2D2N(*this);

    // Calculate material response
    const double l0_pow = std::pow(l0, 2);
    const double E = 0.5 * (std::pow(l, 2) - l0_pow) / l0_pow;
    const double S0 = GetMembranePrestress();
    const double S = young * E + S0;
    const bool is_compressed = S < 0.0 && std::abs(E) > 1.0e-12 ? true : false;

    // Check LHS size
    if (rLeftHandSideMatrix.size1() != mLocalSize || rLeftHandSideMatrix.size2() != mLocalSize) {
        rLeftHandSideMatrix.resize(mLocalSize, mLocalSize, false);
    }

    // Set LHS to zero
    rLeftHandSideMatrix.clear();

    // Calculate required kinematics
    array_1d<double, mLocalSize> B;
    const auto& r_x_0 = r_geom[0].Coordinates();
    const auto& r_x_1 = r_geom[1].Coordinates();
    const double dx = r_x_1[0] - r_x_0[0];
    const double dy = r_x_1[1] - r_x_0[1];
    B[0] = -dx / l0_pow;
    B[1] = -dy / l0_pow;
    B[2] = dx / l0_pow;
    B[3] = dy / l0_pow;

    // Calculate material stiffness matrix
    BoundedMatrix<double, mLocalSize, mLocalSize> K_mat;
    if (is_compressed) {
        K_mat.clear();
    } else {
        const double aux_K_mat = h * l0 * young;
        noalias(K_mat) = aux_K_mat * outer_prod(B, B);
    }

    // Calculate geometric stiffness matrix
    BoundedMatrix<double, mLocalSize, mLocalSize> K_geo;
    if (is_compressed) {
        K_geo.clear();
    } else {
        const double aux_K_geo = h * S / l0;

        K_geo.clear();
        for (IndexType i = 0; i < mNumNodes; ++i) {
            for (IndexType j = 0; j < mNumNodes; ++j) {
                for (IndexType d = 0; d < mDim; ++d) {
                    K_geo(i*mDim + d, j*mDim + d) = (i == j) ? aux_K_geo : - aux_K_geo;
                }
            }
        }
    }

    // Assemble the left hand side contributions
    noalias(rLeftHandSideMatrix) += K_mat + K_geo;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void MembraneElement2D2N::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    // Get element data
    const auto& r_geom = GetGeometry();
    const auto& r_prop = GetProperties();

    // Get membrane data
    const double h = r_prop.GetValue(THICKNESS);
    const double young = r_prop.GetValue(YOUNG_MODULUS);

    // Calculate element lengths
    const double l = StructuralMechanicsElementUtilities::CalculateCurrentLength2D2N(*this);
    const double l0 = StructuralMechanicsElementUtilities::CalculateReferenceLength2D2N(*this);

    // Calculate material response
    const double l0_pow = std::pow(l0, 2);
    const double E = 0.5 * (std::pow(l, 2) - l0_pow) / l0_pow;
    const double S0 = GetMembranePrestress();
    const double S = young * E + S0;
    const bool is_compressed = (S < 0.0 && std::abs(E) > 1.0e-12) ? true : false;

    // Check RHS size
    if (rRightHandSideVector.size() != mLocalSize) {
        rRightHandSideVector.resize(mLocalSize, false);
    }

    // Set RHS to zero
    rRightHandSideVector.clear();

    // Calculate required kinematics
    array_1d<double, mLocalSize> B;
    const auto& r_x_0 = r_geom[0].Coordinates();
    const auto& r_x_1 = r_geom[1].Coordinates();
    const double dx = r_x_1[0] - r_x_0[0];
    const double dy = r_x_1[1] - r_x_0[1];
    B[0] = -dx / l0_pow;
    B[1] = -dy / l0_pow;
    B[2] = dx / l0_pow;
    B[3] = dy / l0_pow;

    // Calculate the internal force contributions
    array_1d<double, mLocalSize> f_int;
    if (is_compressed) {
        f_int.clear();
    } else {
        const double aux_f_int = h * l0 * S;
        noalias(f_int) = aux_f_int * B;
    }

    // Calculate the body force contribution
    array_1d<double, 3> f_ext;
    f_ext = StructuralMechanicsElementUtilities::GetBodyForce(*this, r_geom.IntegrationPoints(GeometryData::IntegrationMethod::GI_GAUSS_1), 0);
    f_ext *= h * l0 * 0.5; // Half the force to each node

    // Assemble the right hand side contributions
    rRightHandSideVector[0] += f_ext[0] - f_int[0];
    rRightHandSideVector[1] += f_ext[1] - f_int[1];
    rRightHandSideVector[2] += f_ext[0] - f_int[2];
    rRightHandSideVector[3] += f_ext[1] - f_int[3];

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void MembraneElement2D2N::CalculateMassMatrix(
    MatrixType& rMassMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    // Get element data
    const auto &r_prop = GetProperties();

    // Get membrane data
    const double h = r_prop.GetValue(THICKNESS);
    const double rho = r_prop.GetValue(DENSITY);

    // Calculate element lengths
    const double l0 = StructuralMechanicsElementUtilities::CalculateReferenceLength2D2N(*this);

    // Check the mass matrix size
    if (rMassMatrix.size1() != mLocalSize || rMassMatrix.size2() != mLocalSize) {
        rMassMatrix.resize(mLocalSize, mLocalSize, false);
    }

    // Calculate the mass matrix contributions
    rMassMatrix.clear();
    const double aux_M = rho * h * l0 / 6.0;
    for (IndexType i = 0; i < mNumNodes; ++i) {
        for (IndexType j = 0; j < mNumNodes; ++j) {
            rMassMatrix(i*mDim, j*mDim) = (i == j) ? 2.0 * aux_M : aux_M;
            rMassMatrix(i*mDim + 1, j*mDim + 1) = (i == j) ? 2.0 * aux_M : aux_M;
        }
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void MembraneElement2D2N::CalculateDampingMatrix(
    MatrixType& rDampingMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    StructuralMechanicsElementUtilities::CalculateRayleighDampingMatrix(
        *this,
        rDampingMatrix,
        rCurrentProcessInfo,
        2.0 * GetGeometry().PointsNumber());

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

int  MembraneElement2D2N::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY;

    int check = Element::Check(rCurrentProcessInfo);

    // Check that the required material data is present
    const auto& r_prop = GetProperties();
    KRATOS_ERROR_IF_NOT(r_prop.Has(THICKNESS)) << "THICKNESS value needs to be provided by properties." << std::endl;
    KRATOS_ERROR_IF_NOT(r_prop.Has(YOUNG_MODULUS)) << "YOUNG_MODULUS value needs to be provided by properties." << std::endl;
    if (r_prop.Has(PRESTRESS_VECTOR)) {
        const auto& r_prestress = r_prop.GetValue(PRESTRESS_VECTOR);
        KRATOS_ERROR_IF(r_prestress[0] < 0.0) << "Prestress needs to be positive." << std::endl;
        for (IndexType i = 1; i < r_prestress.size(); ++i) {
            KRATOS_ERROR_IF(r_prestress[i] > 1.0e-12) << "Only first component of PRESTRESS_VECTOR can be considered." << std::endl;
        }
    }

    return check;

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

const Parameters MembraneElement2D2N::GetSpecifications() const
{
    const Parameters specifications = Parameters(R"({
        "time_integration"           : ["implicit"],
        "framework"                  : "lagrangian",
        "symmetric_lhs"              : true,
        "positive_definite_lhs"      : true,
        "output"                     : {
            "gauss_point"            : [],
            "nodal_historical"       : ["DISPLACEMENT","VELOCITY","ACCELERATION"],
            "nodal_non_historical"   : [],
            "entity"                 : []
        },
        "required_variables"         : ["DISPLACEMENT","VELOCITY","ACCELERATION","POSITIVE_FACE_PRESSURE","NEGATIVE_FACE_PRESSURE"],
        "required_dofs"              : ["DISPLACEMENT_Y"],
        "flags_used"                 : [],
        "compatible_geometries"      : ["Line2D2N"],
        "element_integrates_in_time" : false,
        "compatible_constitutive_laws": {
            "type"        : [],
            "dimension"   : [],
            "strain_size" : []
        },
        "required_polynomial_degree_of_geometry" : -1,
        "documentation"   : "This element implements a simplified non-linear two-dimensional membrane model."
    })");

    return specifications;
}

/***********************************************************************************/
/***********************************************************************************/

double MembraneElement2D2N::GetMembranePrestress() const
{
    const auto& r_prop = GetProperties();
    if (r_prop.Has(PRESTRESS_VECTOR)) {
        return r_prop.GetValue(PRESTRESS_VECTOR)[0];
    } else {
        return 0.0;
    }
}

/***********************************************************************************/
/***********************************************************************************/

void MembraneElement2D2N::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Element );
}

/***********************************************************************************/
/***********************************************************************************/

void MembraneElement2D2N::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Element );
}

} // Namespace Kratos
