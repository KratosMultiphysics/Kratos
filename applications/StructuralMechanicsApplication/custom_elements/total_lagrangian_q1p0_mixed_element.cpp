// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                     license: structural_mechanics_application/license.txt
//
//  Main authors:    Alejandro Cornejo
//

// System includes

// External includes


// Project includes
#include "custom_elements/total_lagrangian_q1p0_mixed_element.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos
{

/***********************************************************************************/
/***********************************************************************************/

TotalLagrangianQ1P0MixedElement::TotalLagrangianQ1P0MixedElement(IndexType NewId, GeometryType::Pointer pGeometry)
    : TotalLagrangian(NewId, pGeometry)
{
    // DO NOT ADD DOFS HERE!!!
}

/***********************************************************************************/
/***********************************************************************************/

TotalLagrangianQ1P0MixedElement::TotalLagrangianQ1P0MixedElement( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
        : TotalLagrangian( NewId, pGeometry, pProperties )
{
    //DO NOT ADD DOFS HERE!!!
}

/***********************************************************************************/
/***********************************************************************************/

Element::Pointer TotalLagrangianQ1P0MixedElement::Create( IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties ) const
{
    return Kratos::make_intrusive<TotalLagrangianQ1P0MixedElement>( NewId, GetGeometry().Create( ThisNodes ), pProperties );
}

//************************************************************************************
//************************************************************************************

Element::Pointer TotalLagrangianQ1P0MixedElement::Create( IndexType NewId,  GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties ) const
{
    return Kratos::make_intrusive<TotalLagrangianQ1P0MixedElement>( NewId, pGeom, pProperties );
}

/***********************************************************************************/
/***********************************************************************************/

TotalLagrangianQ1P0MixedElement::~TotalLagrangianQ1P0MixedElement()
{
}

/***********************************************************************************/
/***********************************************************************************/

Element::Pointer TotalLagrangianQ1P0MixedElement::Clone (
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

    if ( CalculateStiffnessMatrixFlag == true ) { // Calculation of the matrix is required
        if ( rLeftHandSideMatrix.size1() != mat_size )
            rLeftHandSideMatrix.resize( mat_size, mat_size, false );

        noalias( rLeftHandSideMatrix ) = ZeroMatrix( mat_size, mat_size ); //resetting LHS
    }

    // Resizing as needed the RHS
    if ( CalculateResidualVectorFlag == true ) { // Calculation of the matrix is required
        if ( rRightHandSideVector.size() != mat_size )
            rRightHandSideVector.resize( mat_size, false );

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

    const double E = r_props[YOUNG_MODULUS];
    const double nu = r_props[POISSON_RATIO];
    const double bulk_modulus = E / (3.0 * (1.0 - 2.0 * nu));
    const double theta = GetCurrentVolume() / mInitialVolume;
    const double pressure = bulk_modulus * (theta - 1.0);
    Vector I;
    I.resize(6);
    noalias(I) = ZeroVector(6);
    I[0] = 1.0;
    I[1] = 1.0;
    I[2] = 1.0;
    Matrix Emat(6, 6);
    noalias(Emat) = ZeroMatrix(6, 6);
    for (int i = 0; i < 6;i++)
        Emat(i, i) = 1.0;

    // Computing in all integrations points
    for (IndexType point_number = 0; point_number < integration_points.size(); ++point_number) {
        Vector Bv;
        Bv.resize(dimension * number_of_nodes, false);

        Matrix spatial_B(strain_size, dimension * number_of_nodes);
        noalias(spatial_B) = ZeroMatrix(strain_size, dimension * number_of_nodes);

        // Contribution to external forces
        noalias(body_force) = this->GetBodyForce(integration_points, point_number);

        // Compute element kinematics B, F, DN_DX ...
        this->CalculateKinematicVariables(this_kinematic_variables, point_number, this->GetIntegrationMethod());

        // Compute material reponse
        this->CalculateConstitutiveVariables(this_kinematic_variables, this_constitutive_variables, Values, point_number, integration_points, ConstitutiveLaw::StressMeasure_Kirchhoff);

        // Calculating weights for integration on the reference configuration
        int_to_reference_weight = GetIntegrationWeight(integration_points, point_number, this_kinematic_variables.detJ0);

        if (dimension == 2 && r_props.Has(THICKNESS))
            int_to_reference_weight *= r_props[THICKNESS];

        // ----------------------------------------------------------------
        // ----------------------------------------------------------------
        Matrix J, inv_J, DN_Dx;
        CalculateDerivativesOnCurrentConfiguration(J, inv_J, DN_Dx, point_number, this->GetIntegrationMethod());
        // we build Bv
        for (int i = 0; i < number_of_nodes; ++i) {
            const SizeType index = dimension * i;
            if (dimension == 2)
            {
            }
            else
            { // 3D
                Bv[index + 0] = DN_Dx(i, 0);
                Bv[index + 1] = DN_Dx(i, 1);
                Bv[index + 2] = DN_Dx(i, 2);
            }
        } // build end

        // we build spatial B
        for ( IndexType i = 0; i < number_of_nodes; i++ ) {
            const SizeType index = dimension * i;
            spatial_B( 0, index + 0 ) = DN_Dx( i, 0 );
            spatial_B( 1, index + 1 ) = DN_Dx( i, 1 );
            spatial_B( 2, index + 2 ) = DN_Dx( i, 2 );
            spatial_B( 3, index + 0 ) = DN_Dx( i, 1 );
            spatial_B( 3, index + 1 ) = DN_Dx( i, 0 );
            spatial_B( 4, index + 1 ) = DN_Dx( i, 2 );
            spatial_B( 4, index + 2 ) = DN_Dx( i, 1 );
            spatial_B( 5, index + 0 ) = DN_Dx( i, 2 );
            spatial_B( 5, index + 2 ) = DN_Dx( i, 0 );
        }

        const double I1 = this_constitutive_variables.StressVector[0] + this_constitutive_variables.StressVector[1] + this_constitutive_variables.StressVector[2];
        const double p = I1 / 3.0;
        this_constitutive_variables.StressVector[0] -= p;
        this_constitutive_variables.StressVector[1] -= p;
        this_constitutive_variables.StressVector[2] -= p;
        // ----------------------------------------------------------------
        // ----------------------------------------------------------------
        // noalias(this_kinematic_variables.B) = spatial_B;

        if (CalculateStiffnessMatrixFlag) { // Calculation of the matrix is required
            // Contributions to stiffness matrix calculated on the reference config
            /* Material stiffness matrix */
            this->CalculateAndAddKm(rLeftHandSideMatrix, spatial_B, this_constitutive_variables.D + pressure * this_kinematic_variables.detF * (outer_prod(I, I) - 2.0 * Emat), int_to_reference_weight);
            // this->CalculateAndAddKm(rLeftHandSideMatrix, spatial_B, this_constitutive_variables.D, int_to_reference_weight);

            /* Geometric stiffness matrix */
            this->CalculateAndAddKg(rLeftHandSideMatrix, DN_Dx, this_constitutive_variables.StressVector + I * pressure * this_kinematic_variables.detF, int_to_reference_weight);
            // this->CalculateAndAddKg(rLeftHandSideMatrix, this_kinematic_variables.DN_DX, this_constitutive_variables.StressVector, int_to_reference_weight);

            noalias(rLeftHandSideMatrix) += int_to_reference_weight * outer_prod(Bv, Bv) * bulk_modulus * std::pow(this_kinematic_variables.detF, 2.0) / mInitialVolume;
        }

        if ( CalculateResidualVectorFlag ) { // Calculation of the matrix is required
            noalias(this_kinematic_variables.B) = spatial_B;
            this->CalculateAndAddResidualVector(rRightHandSideVector, this_kinematic_variables, rCurrentProcessInfo, body_force, this_constitutive_variables.StressVector, int_to_reference_weight);

            // Pressure contribution
            noalias(rRightHandSideVector) -= int_to_reference_weight * this_kinematic_variables.detF * Bv * pressure;
        }
    }

    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

void TotalLagrangianQ1P0MixedElement::InitializeMaterial()
{
    KRATOS_TRY

    BaseType::InitializeMaterial();

    mInitialVolume = GetCurrentVolume();

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

double TotalLagrangianQ1P0MixedElement::GetCurrentVolume() const
{
    if (GetGeometry().WorkingSpaceDimension() == 3)
        return GetGeometry().Volume();
    else
        return GetGeometry().Area();
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

} // Namespace Kratos


