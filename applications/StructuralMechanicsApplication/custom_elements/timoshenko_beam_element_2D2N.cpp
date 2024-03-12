// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:  Alejandro Cornejo
//
//

// System includes

// External includes

// Project includes
// #include "utilities/geometry_utilities.h"
// #include "utilities/math_utils.h"

// Application includes
#include "custom_elements/timoshenko_beam_element_2D2N.h"
#include "structural_mechanics_application_variables.h"
#include "custom_utilities/structural_mechanics_element_utilities.h"
#include "custom_utilities/constitutive_law_utilities.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos
{

void TimoshenkoBeamElement2D2N::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Initialization should not be done again in a restart!
    if (!rCurrentProcessInfo[IS_RESTARTED]) {
        if (this->UseGeometryIntegrationMethod()) {
            if (GetProperties().Has(INTEGRATION_ORDER) ) {
                const SizeType integration_order = GetProperties()[INTEGRATION_ORDER];
                switch ( integration_order )
                {
                case 1:
                    mThisIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_1;
                    break;
                case 2:
                    mThisIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_2;
                    break;
                case 3:
                    mThisIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_3;
                    break;
                case 4:
                    mThisIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_4;
                    break;
                case 5:
                    mThisIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_5;
                    break;
                default:
                    KRATOS_WARNING("TimoshenkoBeamElement2D2N") << "Integration order "
                        << integration_order << " is not available, using default integration order for the geometry" << std::endl;
                    mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod();
                }
            } else {
                mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod();
            }
        }

        const auto& r_integration_points = this->IntegrationPoints(mThisIntegrationMethod);

        // Constitutive Law initialisation
        if (mConstitutiveLawVector.size() != r_integration_points.size())
            mConstitutiveLawVector.resize(r_integration_points.size());
        InitializeMaterial();
    }
    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void TimoshenkoBeamElement2D2N::InitializeMaterial()
{
    KRATOS_TRY

    if (GetProperties()[CONSTITUTIVE_LAW] != nullptr) {
        const auto& r_geometry   = GetGeometry();
        const auto& r_properties = GetProperties();
        auto N_values            = Vector();
        for (IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number) {
            mConstitutiveLawVector[point_number] = GetProperties()[CONSTITUTIVE_LAW]->Clone();
            mConstitutiveLawVector[point_number]->InitializeMaterial(r_properties, r_geometry, N_values);
        }
    } else
        KRATOS_ERROR << "A constitutive law needs to be specified for the element with ID " << this->Id() << std::endl;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

Element::Pointer TimoshenkoBeamElement2D2N::Clone(
    IndexType NewId,
    NodesArrayType const& rThisNodes
    ) const
{
    KRATOS_TRY

    TimoshenkoBeamElement2D2N::Pointer p_new_elem = Kratos::make_intrusive<TimoshenkoBeamElement2D2N>
        (NewId, GetGeometry().Create(rThisNodes), pGetProperties());
    p_new_elem->SetData(this->GetData());
    p_new_elem->Set(Flags(*this));

    // Currently selected integration methods
    p_new_elem->SetIntegrationMethod(mThisIntegrationMethod);

    // The vector containing the constitutive laws
    p_new_elem->SetConstitutiveLawVector(mConstitutiveLawVector);

    return p_new_elem;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void TimoshenkoBeamElement2D2N::EquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo
    ) const
{
    KRATOS_TRY;

    const auto& r_geom = GetGeometry();
    const SizeType number_of_nodes = r_geom.size();
    const SizeType dofs_per_node = GetDoFsPerNode(); // u, v, theta

    if (rResult.size() != dofs_per_node * number_of_nodes)
        rResult.resize(dofs_per_node * number_of_nodes, false);

    const SizeType pos = r_geom[0].GetDofPosition(DISPLACEMENT_X);

    for (IndexType i = 0; i < number_of_nodes; ++i) {
        const SizeType index = i * 2;
        rResult[index]     = r_geom[i].GetDof(DISPLACEMENT_X, pos    ).EquationId();
        rResult[index + 1] = r_geom[i].GetDof(DISPLACEMENT_Y, pos + 1).EquationId();
        rResult[index + 2] = r_geom[i].GetDof(ROTATION_Z,     pos + 2).EquationId();
    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void TimoshenkoBeamElement2D2N::GetDofList(
    DofsVectorType& rElementalDofList,
    const ProcessInfo& rCurrentProcessInfo
    ) const
{
    KRATOS_TRY;

    const auto& r_geom = GetGeometry();
    const SizeType number_of_nodes = r_geom.size();
    const SizeType dofs_per_node = GetDoFsPerNode(); // u, v, theta
    rElementalDofList.resize(0);
    rElementalDofList.reserve(dofs_per_node * number_of_nodes);

    for (IndexType i = 0; i < number_of_nodes; ++i) {
        rElementalDofList.push_back(r_geom[i].pGetDof(DISPLACEMENT_X));
        rElementalDofList.push_back(r_geom[i].pGetDof(DISPLACEMENT_Y));
        rElementalDofList.push_back(r_geom[i].pGetDof(ROTATION_Z));
    }
    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

const double TimoshenkoBeamElement2D2N::CalculatePhi(ConstitutiveLaw::Parameters &rValues)
{
    const auto &r_material_properties = rValues.GetMaterialProperties();
    const double E   = r_material_properties[YOUNG_MODULUS];
    const double A   = r_material_properties[CROSS_AREA];
    const double I   = r_material_properties[IZ];
    const double k_s = r_material_properties.Has(SHEAR_CORRECTION_XY) ? r_material_properties[SHEAR_CORRECTION_XY] : 5.0 / 6.0; // We assume rectangular shape
    const double A_s = k_s * A;
    const double G   = ConstitutiveLawUtilities<3>::CalculateShearModulus(rValues);
    const double L   = StructuralMechanicsElementUtilities::CalculateReferenceLength2D2N(*this);
    return 12.0 * E * I / (G * A_s * std::pow(L, 2));
}

/***********************************************************************************/
/***********************************************************************************/

void TimoshenkoBeamElement2D2N::GetShapeFunctionsValues(
    VectorType& rN,
    const double Length,
    const double Phi,
    const double xi
    )
{
    if (rN.size() != 4)
        rN.resize(4, false);
    const double one_plus_phi = 1.0 + Phi;
    const double xi_square = xi * xi;
    rN[0] = (xi - 1.0) * (xi + xi_square - 2.0 * one_plus_phi) / (4.0 * one_plus_phi);
    rN[1] = (1.0 - xi_square) * (1.0 - xi + Phi) * Length / (8.0 * one_plus_phi);
    rN[2] = (1.0 + xi) * (xi - xi_square + 2.0 * one_plus_phi) / (4.0 * one_plus_phi);
    rN[3] = (xi_square - 1.0) * (1.0 + xi + Phi) * Length / (8.0 * one_plus_phi);
}

/***********************************************************************************/
/***********************************************************************************/

void TimoshenkoBeamElement2D2N::GetFirstDerivativesShapeFunctionsValues(
    VectorType& rN,
    const double Length,
    const double Phi,
    const double xi
    )
{
    if (rN.size() != 4)
        rN.resize(4, false);
    const double one_plus_phi = 1.0 + Phi;
    const double xi_square = xi * xi;
    rN[0] = (-6.0 + 6.0 * xi_square - 4.0 * Phi) / (4.0 * one_plus_phi * Length);
    rN[1] = (-1.0 + 3.0 * xi_square - 2.0 * xi * one_plus_phi) / (4.0 * one_plus_phi);
    rN[2] = (6.0 - 6.0 * xi_square + 4.0 * Phi) / (4.0 * one_plus_phi * Length);
    rN[3] = (-1.0 + 3.0 * xi_square + 2.0 * xi * one_plus_phi) / (4.0 * one_plus_phi);
}

/***********************************************************************************/
/***********************************************************************************/

void TimoshenkoBeamElement2D2N::GetSecondDerivativesShapeFunctionsValues(
    VectorType& rN,
    const double Length,
    const double Phi,
    const double xi
    )
{
    if (rN.size() != 4)
        rN.resize(4, false);
    const double one_plus_phi = 1.0 + Phi;
    const double xi_square = xi * xi;
    const double L_square = std::pow(Length, 2);
    rN[0] = 6.0 * xi / (one_plus_phi * L_square);
    rN[1] = (-1.0 + 3.0 * xi - Phi) / (one_plus_phi * Length);
    rN[2] = -6.0 * xi / (one_plus_phi * L_square);
    rN[3] = (1.0 + 3.0 * xi + Phi) / (one_plus_phi * Length);
}

/***********************************************************************************/
/***********************************************************************************/

void TimoshenkoBeamElement2D2N::GetThirdDerivativesShapeFunctionsValues(
    VectorType& rN,
    const double Length,
    const double Phi,
    const double xi
    )
{
    if (rN.size() != 4)
        rN.resize(4, false);
    const double one_plus_phi = 1.0 + Phi;
    const double xi_square = xi * xi;
    const double L_square = std::pow(Length, 2);
    rN[0] = 12.0  / (one_plus_phi * Length * L_square);
    rN[1] = 6.0   / (one_plus_phi * L_square);
    rN[2] = -12.0 / (one_plus_phi * Length * L_square);
    rN[3] = 6.0   / (one_plus_phi * L_square);
}

// /***********************************************************************************/
// /***********************************************************************************/

void TimoshenkoBeamElement2D2N::GetNThetaShapeFunctionsValues(
    VectorType& rN,
    const double Length,
    const double Phi,
    const double xi
    )
{
    const double one_plus_phi = 1.0 + Phi;
    if (rN.size() != 4)
        rN.resize(4, false);
    rN[0] = (3.0 * xi * xi - 3.0) / (2.0 * one_plus_phi * Length);
    rN[1] = (xi - 1.0) * (1.0 + 3.0 * xi - 2.0 * Phi) / (4.0 * one_plus_phi);
    rN[2] = (3.0 - 3.0 * xi * xi) / (2 * one_plus_phi * Length);
    rN[3] = (1.0 + xi) * (3.0 * xi - 1.0 + 2.0 * Phi) / (4.0 * one_plus_phi);
}

/***********************************************************************************/
/***********************************************************************************/

void TimoshenkoBeamElement2D2N::GetFirstDerivativesNThetaShapeFunctionsValues(
    VectorType& rN,
    const double Length,
    const double Phi,
    const double xi
    )
{
    const double one_plus_phi = 1.0 + Phi;
    if (rN.size() != 4)
        rN.resize(4, false);
    rN[0] = 3.0 * xi / (one_plus_phi * Length);
    rN[1] = (-0.5 * Phi + 1.5 * xi - 0.5) / (Phi + 1.0);
    rN[2] = (-3.0 * xi) / (Length * Phi + Length);
    rN[3] = (0.5 * Phi + 1.5 * xi + 0.5) / (Phi + 1.0);
    rN *= (2.0 / Length);
}

/***********************************************************************************/
/***********************************************************************************/

void TimoshenkoBeamElement2D2N::GetNu0ShapeFunctionsValues(
    VectorType& rN,
    const double Length,
    const double Phi,
    const double xi
    )
{
    if (rN.size() != 2)
        rN.resize(2, false);
    rN[0] = 0.5 * (1.0 - xi);
    rN[1] = 0.5 * (1.0 + xi);
}

/***********************************************************************************/
/***********************************************************************************/

void TimoshenkoBeamElement2D2N::GetFirstDerivativesNu0ShapeFunctionsValues(
    VectorType& rN,
    const double Length,
    const double Phi,
    const double xi
    )
{
    if (rN.size() != 2)
        rN.resize(2, false);
    rN[0] = -1.0 / Length;
    rN[1] =  1.0 / Length;
}

/***********************************************************************************/
/***********************************************************************************/

void TimoshenkoBeamElement2D2N::GetNodalValuesVector(VectorType& rNodalValues)
{
    if (rNodalValues.size() != 6)
        rNodalValues.resize(6, false);
    const auto &r_geom = GetGeometry();

    rNodalValues[0] = r_geom[0].FastGetSolutionStepValue(DISPLACEMENT_X);
    rNodalValues[1] = r_geom[0].FastGetSolutionStepValue(DISPLACEMENT_Y);
    rNodalValues[2] = r_geom[0].FastGetSolutionStepValue(ROTATION_Z);
    rNodalValues[3] = r_geom[1].FastGetSolutionStepValue(DISPLACEMENT_X);
    rNodalValues[4] = r_geom[1].FastGetSolutionStepValue(DISPLACEMENT_Y);
    rNodalValues[5] = r_geom[1].FastGetSolutionStepValue(ROTATION_Z);
}

/***********************************************************************************/
/***********************************************************************************/

const double TimoshenkoBeamElement2D2N::CalculateAxialStrain(
    const double Length,
    const double Phi,
    const double xi,
    const VectorType& rNodalValues
    )
{
    VectorType N_u0_derivatives(2);
    GetFirstDerivativesNu0ShapeFunctionsValues(N_u0_derivatives, Length, Phi, xi);
    return N_u0_derivatives[0] * rNodalValues[0] + N_u0_derivatives[1] * rNodalValues[3];
}

/***********************************************************************************/
/***********************************************************************************/

const double TimoshenkoBeamElement2D2N::CalculateShearStrain(
    const double Length,
    const double Phi,
    const double xi,
    const VectorType& rNodalValues
    )
{
    VectorType N_derivatives(4), N_theta(4);
    GetFirstDerivativesShapeFunctionsValues(N_derivatives, Length, Phi, xi);
    GetNThetaShapeFunctionsValues(N_theta, Length, Phi, xi);
    const Vector N_s = N_derivatives - N_theta;
    return N_s[0] * rNodalValues[1] + N_s[1] * rNodalValues[2] + N_s[2] * rNodalValues[4] + 
           N_s[3] * rNodalValues[5];
}

/***********************************************************************************/
/***********************************************************************************/

const double TimoshenkoBeamElement2D2N::CalculateBendingCurvature(
    const double Length,
    const double Phi,
    const double xi,
    const VectorType& rNodalValues
    )
{
    VectorType N_theta_derivatives(4);
    GetFirstDerivativesNThetaShapeFunctionsValues(N_theta_derivatives, Length, Phi, xi);
    return N_theta_derivatives[0] * rNodalValues[1] + N_theta_derivatives[1] * rNodalValues[2] +
           N_theta_derivatives[2] * rNodalValues[4] + N_theta_derivatives[3] * rNodalValues[5];
}

/***********************************************************************************/
/***********************************************************************************/

void TimoshenkoBeamElement2D2N::CalculateLocalSystem(
    MatrixType& rLHS,
    VectorType& rRHS,
    const ProcessInfo& rProcessInfo
    )
{
    KRATOS_TRY;
    const auto &r_geometry = GetGeometry();
    const SizeType number_of_nodes = r_geometry.size();
    const SizeType mat_size = GetDoFsPerNode() * number_of_nodes;

    if (rLHS.size1() != mat_size || rLHS.size2() != mat_size) {
        rLHS.resize(mat_size, false);
        rLHS.clear();
    }

    if (rRHS.size() != mat_size) {
        rRHS.resize(mat_size, false);
        rRHS.clear();
    }

    const auto& integration_points = IntegrationPoints(GetIntegrationMethod());

    ConstitutiveLaw::Parameters cl_values(r_geometry, GetProperties(), rProcessInfo);
    auto &r_cl_options = cl_values.GetOptions();
    r_cl_options.Set(ConstitutiveLaw::COMPUTE_STRESS             , true);
    r_cl_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

    const double Phi    = CalculatePhi(cl_values);
    const double length = StructuralMechanicsElementUtilities::CalculateReferenceLength2D2N(*this);
    const double J      = 0.5 * length;

    VectorType strain_vector(3);
    strain_vector.clear();
    cl_values.SetStrainVector(strain_vector);
    VectorType nodal_values(6);
    GetNodalValuesVector(nodal_values);
    VectorType global_size_N(6), N_u_derivatives(2), N_theta_derivatives(4), N_theta(4), N_derivatives(4)

    for (SizeType IP = 0; IP < integration_points.size(); ++IP) {
        global_size_N.clear();
        const double xi     = integration_points[IP].X();
        const double weight = integration_points[IP].Weight();

        strain_vector[0] = CalculateAxialStrain(length, Phi, xi, nodal_values);
        strain_vector[1] = CalculateBendingCurvature(length, Phi, xi, nodal_values);
        strain_vector[2] = CalculateShearStrain(length, Phi, xi, nodal_values);

        mConstitutiveLawVector[IP]->CalculateMaterialResponseCauchy(cl_values);
        const Vector &r_generalized_stresses = cl_values.GetStressVector();
        const double N = r_generalized_stresses[0];
        const double M = r_generalized_stresses[1];
        const double V = r_generalized_stresses[2];

        GetFirstDerivativesNu0ShapeFunctionsValues(N_u_derivatives, length, Phi, xi);
        GetFirstDerivativesNThetaShapeFunctionsValues(N_theta_derivatives, length, Phi, xi);
        GetNThetaShapeFunctionsValues(N_theta, length, Phi, xi);
        GetFirstDerivativesShapeFunctionsValues(N_derivatives, length, Phi, xi);
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void TimoshenkoBeamElement2D2N::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void TimoshenkoBeamElement2D2N::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void TimoshenkoBeamElement2D2N::CalculateOnIntegrationPoints(
    const Variable<double>& rVariable,
    std::vector<double>& rOutput,
    const ProcessInfo& rCurrentProcessInfo
    )
{

}

/***********************************************************************************/
/***********************************************************************************/

void TimoshenkoBeamElement2D2N::CalculateOnIntegrationPoints(
    const Variable<ConstitutiveLaw::Pointer>& rVariable,
    std::vector<ConstitutiveLaw::Pointer>& rValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    if (rVariable == CONSTITUTIVE_LAW) {
        const SizeType integration_points_number = mConstitutiveLawVector.size();
        if (rValues.size() != integration_points_number) {
            rValues.resize(integration_points_number);
        }
        for (IndexType point_number = 0; point_number < integration_points_number; ++point_number) {
            rValues[point_number] = mConstitutiveLawVector[point_number];
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

int  TimoshenkoBeamElement2D2N::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY;

    // TODO
    return 0;

    KRATOS_CATCH( "" );
}

//***********************************************************************
//***********************************************************************

double TimoshenkoBeamElement2D2N::GetIntegrationWeight(
    const GeometryType::IntegrationPointsArrayType& rThisIntegrationPoints,
    const SizeType PointNumber,
    const double detJ
    ) const
{
    return rThisIntegrationPoints[PointNumber].Weight() * detJ;
}

/***********************************************************************************/
/***********************************************************************************/

void TimoshenkoBeamElement2D2N::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);
    int IntMethod = int(this->GetIntegrationMethod());
    rSerializer.save("IntegrationMethod",IntMethod);
    rSerializer.save("ConstitutiveLawVector", mConstitutiveLawVector);
}

/***********************************************************************************/
/***********************************************************************************/

void TimoshenkoBeamElement2D2N::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
    int IntMethod;
    rSerializer.load("IntegrationMethod",IntMethod);
    mThisIntegrationMethod = IntegrationMethod(IntMethod);
    rSerializer.load("ConstitutiveLawVector", mConstitutiveLawVector);
}

/***********************************************************************************/
/***********************************************************************************/

} // Namespace Kratos
