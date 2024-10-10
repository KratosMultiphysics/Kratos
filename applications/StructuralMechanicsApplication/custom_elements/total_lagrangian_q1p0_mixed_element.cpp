// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Alejandro Cornejo
//

// System includes

// External includes

// Project includes
#include "custom_elements/total_lagrangian_q1p0_mixed_element.h"

namespace Kratos
{

/***********************************************************************************/
/***********************************************************************************/

TotalLagrangianQ1P0MixedElement::TotalLagrangianQ1P0MixedElement(
    IndexType NewId,
    GeometryType::Pointer pGeometry
    ) : TotalLagrangian(NewId, pGeometry)
{
    // DO NOT ADD DOFS HERE!!!
    this->SetValue(PRESSURE, 0.0);
}

/***********************************************************************************/
/***********************************************************************************/

TotalLagrangianQ1P0MixedElement::TotalLagrangianQ1P0MixedElement(
    IndexType NewId,
    GeometryType::Pointer pGeometry,
    PropertiesType::Pointer pProperties
    ) : TotalLagrangian(NewId, pGeometry, pProperties)
{
    //DO NOT ADD DOFS HERE!!!
    this->SetValue(PRESSURE, 0.0);
}

/***********************************************************************************/
/***********************************************************************************/

Element::Pointer TotalLagrangianQ1P0MixedElement::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties
    ) const
{
    return Kratos::make_intrusive<TotalLagrangianQ1P0MixedElement>(NewId, GetGeometry().Create(ThisNodes), pProperties);
}

//************************************************************************************
//************************************************************************************

Element::Pointer TotalLagrangianQ1P0MixedElement::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties
    ) const
{
    return Kratos::make_intrusive<TotalLagrangianQ1P0MixedElement>(NewId, pGeom, pProperties);
}

/***********************************************************************************/
/***********************************************************************************/

TotalLagrangianQ1P0MixedElement::~TotalLagrangianQ1P0MixedElement()
{
}

/***********************************************************************************/
/***********************************************************************************/

Element::Pointer TotalLagrangianQ1P0MixedElement::Clone(
    IndexType NewId,
    NodesArrayType const& rThisNodes
    ) const
{
    KRATOS_TRY

    TotalLagrangianQ1P0MixedElement::Pointer p_new_elem = Kratos::make_intrusive<TotalLagrangianQ1P0MixedElement>(NewId, GetGeometry().Create(rThisNodes), pGetProperties());
    p_new_elem->SetData(this->GetData());
    p_new_elem->Set(Flags(*this));

    // Currently selected integration methods
    p_new_elem->SetIntegrationMethod(BaseType::mThisIntegrationMethod);

    // The vector containing the constitutive laws
    p_new_elem->SetConstitutiveLawVector(BaseType::mConstitutiveLawVector);

    return p_new_elem;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void TotalLagrangianQ1P0MixedElement::CalculateAll(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo,
    const bool CalculateStiffnessMatrixFlag,
    const bool CalculateResidualVectorFlag
    )
{
    KRATOS_TRY;

    const auto &r_geometry = GetGeometry();
    const SizeType number_of_nodes = r_geometry.size();
    const SizeType dimension = r_geometry.WorkingSpaceDimension();
    const auto strain_size = GetStrainSize();
    const auto &r_props = GetProperties();

    KinematicVariables this_kinematic_variables(strain_size, dimension, number_of_nodes);
    ConstitutiveVariables this_constitutive_variables(strain_size);

    // Resizing as needed the LHS
    const SizeType mat_size = number_of_nodes * dimension;

    if (CalculateStiffnessMatrixFlag) { // Calculation of the matrix is required
        if (rLeftHandSideMatrix.size1() != mat_size)
            rLeftHandSideMatrix.resize(mat_size, mat_size, false);
        noalias(rLeftHandSideMatrix) = ZeroMatrix(mat_size, mat_size); // resetting LHS
    }

    // Resizing as needed the RHS
    if (CalculateResidualVectorFlag) { // Calculation of the matrix is required
        if (rRightHandSideVector.size() != mat_size)
            rRightHandSideVector.resize(mat_size, false);
        noalias(rRightHandSideVector) = ZeroVector( mat_size ); //resetting RHS
    }

    // Reading integration points
    const GeometryType::IntegrationPointsArrayType& integration_points = r_geometry.IntegrationPoints(this->GetIntegrationMethod());

    ConstitutiveLaw::Parameters Values(r_geometry,r_props,rCurrentProcessInfo);

    // Set constitutive law flags:
    Flags& ConstitutiveLawOptions=Values.GetOptions();
    ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, UseElementProvidedStrain());
    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    if ( CalculateStiffnessMatrixFlag ) {
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
    } else {
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);
    }

    // If strain has to be computed inside of the constitutive law with PK2
    Values.SetStrainVector(this_constitutive_variables.StrainVector); //this is the input  parameter

    // Some declarations
    array_1d<double, 3> body_force;
    double int_to_reference_weight;

    const double bulk_modulus = CalculateBulkModulus(r_props);

    double Kpp = 0.0;
    double Fp = 0.0;
    Vector Kup(mat_size);
    noalias(Kup) = ZeroVector(mat_size);
    Matrix inv_C(dimension, dimension);
    Vector inv_c_voigt(strain_size);
    double det;

    // Computing in all integrations points
    for (IndexType point_number = 0; point_number < integration_points.size(); ++point_number) {

        // Contribution to external forces
        noalias(body_force) = this->GetBodyForce(integration_points, point_number);

        // Compute element kinematics B, F, DN_DX ...
        this->CalculateKinematicVariables(this_kinematic_variables, point_number, this->GetIntegrationMethod());

        // Compute material response
        this->CalculateConstitutiveVariables(this_kinematic_variables, this_constitutive_variables, Values, point_number, integration_points, this->GetStressMeasure(), false);

        const Matrix& r_C = prod(trans(this_kinematic_variables.F), this_kinematic_variables.F);
        MathUtils<double>::InvertMatrix3(r_C, inv_C, det);
        noalias(inv_c_voigt) = MathUtils<double>::StrainTensorToVector(inv_C, strain_size);
        if (dimension == 2) {
            inv_c_voigt[2] /= 2.0;
        } else {
            inv_c_voigt[3] /= 2.0;
            inv_c_voigt[4] /= 2.0;
            inv_c_voigt[5] /= 2.0;
        }

        // Calculating weights for integration on the reference configuration
        int_to_reference_weight = GetIntegrationWeight(integration_points, point_number, this_kinematic_variables.detJ0);

        if (dimension == 2 && r_props.Has(THICKNESS))
            int_to_reference_weight *= r_props[THICKNESS];

        // we compute u-p entities
        noalias(Kup) -= int_to_reference_weight * this_kinematic_variables.detF * prod(trans(this_kinematic_variables.B), inv_c_voigt);
        Kpp          -= int_to_reference_weight / bulk_modulus;
        Fp           -= int_to_reference_weight * ((this_kinematic_variables.detF - 1.0) + (this->GetValue(PRESSURE) / bulk_modulus));

        if (CalculateStiffnessMatrixFlag) { // Calculation of the matrix is required
            // Contributions to stiffness matrix calculated on the reference config
            /* Material stiffness matrix */
            this->CalculateAndAddKm(rLeftHandSideMatrix, this_kinematic_variables.B, this_constitutive_variables.D, int_to_reference_weight);

            /* Geometric stiffness matrix */
            this->CalculateAndAddKg(rLeftHandSideMatrix, this_kinematic_variables.DN_DX, this_constitutive_variables.StressVector, int_to_reference_weight);
        }

        if (CalculateResidualVectorFlag) { // Calculation of the matrix is required
            this->CalculateAndAddResidualVector(rRightHandSideVector, this_kinematic_variables, rCurrentProcessInfo, body_force, this_constitutive_variables.StressVector, int_to_reference_weight);
        }
    } // IP loop

    if (CalculateStiffnessMatrixFlag)
        noalias(rLeftHandSideMatrix) -= outer_prod(Kup, Kup) / Kpp;
    if (CalculateResidualVectorFlag)
        noalias(rRightHandSideVector) += Kup * Fp / Kpp;

    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

void TotalLagrangianQ1P0MixedElement::FinalizeNonLinearIteration(
    const ProcessInfo &rCurrentProcessInfo
    )
{
    KRATOS_TRY

    const auto &r_geometry = GetGeometry();
    const SizeType number_of_nodes = r_geometry.size();
    const SizeType dimension = r_geometry.WorkingSpaceDimension();
    const auto strain_size = GetStrainSize();
    const auto &r_props = GetProperties();

    KinematicVariables this_kinematic_variables(strain_size, dimension, number_of_nodes);
    ConstitutiveVariables this_constitutive_variables(strain_size);

    // Resizing as needed the LHS
    const SizeType mat_size = number_of_nodes * dimension;

    // Reading integration points
    const GeometryType::IntegrationPointsArrayType& integration_points = r_geometry.IntegrationPoints(this->GetIntegrationMethod());

    // Some declarations
    double int_to_reference_weight;
    const double bulk_modulus = CalculateBulkModulus(r_props);

    double Kpp = 0.0;
    double Fp = 0.0;
    Vector Kup(mat_size);
    noalias(Kup) = ZeroVector(mat_size);
    Matrix inv_C(dimension, dimension);
    Vector inv_c_voigt(strain_size);
    double det;

    for (IndexType point_number = 0; point_number < integration_points.size(); ++point_number) {

        // Compute element kinematics B, F, DN_DX ...
        this->CalculateKinematicVariables(this_kinematic_variables, point_number, this->GetIntegrationMethod());

        const Matrix& r_C = prod(trans(this_kinematic_variables.F), this_kinematic_variables.F);
        MathUtils<double>::InvertMatrix3(r_C, inv_C, det);
        noalias(inv_c_voigt) = MathUtils<double>::StrainTensorToVector(inv_C, strain_size);
        if (dimension == 2) {
            inv_c_voigt[2] /= 2.0;
        } else {
            inv_c_voigt[3] /= 2.0;
            inv_c_voigt[4] /= 2.0;
            inv_c_voigt[5] /= 2.0;
        }

        // Calculating weights for integration on the reference configuration
        int_to_reference_weight = GetIntegrationWeight(integration_points, point_number, this_kinematic_variables.detJ0);

        if (dimension == 2 && r_props.Has(THICKNESS))
            int_to_reference_weight *= r_props[THICKNESS];

        // we compute u-p entities
        noalias(Kup) -= int_to_reference_weight * this_kinematic_variables.detF * prod(trans(this_kinematic_variables.B), inv_c_voigt);
        Kpp          -= int_to_reference_weight / bulk_modulus;
        Fp           -= int_to_reference_weight * ((this_kinematic_variables.detF - 1.0) + (this->GetValue(PRESSURE) / bulk_modulus));

    } // IP loop

    Vector displ, displ_old;
    GetValuesVector(displ, 0);
    GetValuesVector(displ_old, 1);
    this->GetValue(PRESSURE) -= (Fp + inner_prod(Kup, displ - displ_old)) / Kpp;
    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

double TotalLagrangianQ1P0MixedElement::CalculateBulkModulus(const Properties& rProperties)
{
    const double E = rProperties[YOUNG_MODULUS];
    const double nu = rProperties[POISSON_RATIO];
    return E / (3.0 * (1.0 - 2.0 * nu));
}

/***********************************************************************************/
/***********************************************************************************/

void TotalLagrangianQ1P0MixedElement::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, TotalLagrangian );
}

/***********************************************************************************/
/***********************************************************************************/

void TotalLagrangianQ1P0MixedElement::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, TotalLagrangian );
}

/***********************************************************************************/
/***********************************************************************************/

void TotalLagrangianQ1P0MixedElement::CalculateOnIntegrationPoints(
    const Variable<double>& rVariable,
    std::vector<double>& rOutput,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    BaseType::CalculateOnIntegrationPoints(rVariable, rOutput, rCurrentProcessInfo);

    const GeometryType::IntegrationPointsArrayType &integration_points = GetGeometry().IntegrationPoints(this->GetIntegrationMethod());

    const std::size_t number_of_integration_points = integration_points.size();

    if (rOutput.size() != number_of_integration_points)
        rOutput.resize(number_of_integration_points, false);

    if (rVariable == PRESSURE) {
        const double pressure  = this->GetValue(PRESSURE);
        for (IndexType i = 0; i < number_of_integration_points; i++) {
            rOutput[i] = pressure;
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

int TotalLagrangianQ1P0MixedElement::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY;
    auto &r_props = GetProperties();
    KRATOS_ERROR_IF_NOT(r_props.Has(YOUNG_MODULUS)) << "The YOUNG_MODULUS is not defined in the properties..." << std::endl;
    KRATOS_ERROR_IF_NOT(r_props.Has(POISSON_RATIO)) << "The POISSON_RATIO is not defined in the properties..." << std::endl;
    BaseType::Check(rCurrentProcessInfo);
    return 0;
    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

} // Namespace Kratos


