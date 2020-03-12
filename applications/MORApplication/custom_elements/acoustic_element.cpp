// KRATOS
//
//  License:		 BSD License
//					 license: ../../license.txt
//
//  Main authors:    Ramsubramanian Pazhanisamy
//                   Ricky Aristio
//

// System includes

// External includes


// Project includes
#include "includes/checks.h"

// Application includes
#include "custom_elements/acoustic_element.h"
//#include "custom_utilities/structural_mechanics_element_utilities.h"
#include "mor_application_variables.h"


namespace Kratos
{
AcousticElement::AcousticElement( IndexType NewId, GeometryType::Pointer pGeometry )
    : Element( NewId, pGeometry )
{
    //DO NOT ADD DOFS HERE!!!
}

/***********************************************************************************/
/***********************************************************************************/

AcousticElement::AcousticElement( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
        : Element( NewId, pGeometry, pProperties )
{
    //DO NOT ADD DOFS HERE!!!
}

/***********************************************************************************/
/***********************************************************************************/

Element::Pointer AcousticElement::Create( IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties ) const
{
    return Kratos::make_intrusive<AcousticElement>( NewId, GetGeometry().Create( ThisNodes ), pProperties );
}

/***********************************************************************************/
/***********************************************************************************/

Element::Pointer AcousticElement::Create( IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties ) const
{
    return Kratos::make_intrusive<AcousticElement>( NewId, pGeom, pProperties );
}

/***********************************************************************************/
/***********************************************************************************/

AcousticElement::~AcousticElement()
{
}

/***********************************************************************************/
/***********************************************************************************/

Element::Pointer AcousticElement::Clone (
    IndexType NewId,
    NodesArrayType const& rThisNodes
    ) const
{
    KRATOS_TRY

    AcousticElement::Pointer p_new_elem = Kratos::make_intrusive<AcousticElement>(NewId, GetGeometry().Create(rThisNodes), pGetProperties());
    p_new_elem->SetData(this->GetData());
    p_new_elem->Set(Flags(*this));

    return p_new_elem;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

int  AcousticElement::Check( const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    int ier = Element::Check(rCurrentProcessInfo);

     if(ier != 0) return ier;
  
      // Check that all required variables have been registered
      KRATOS_CHECK_VARIABLE_KEY(ACOUSTIC_PRESSURE)
  
      const SizeType number_of_points = GetGeometry().size();  //added cornejo
      // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
      for ( IndexType i = 0; i < number_of_points; i++ )
      {
          const NodeType &rnode = this->GetGeometry()[i];
          KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(ACOUSTIC_PRESSURE,rnode)
          KRATOS_CHECK_DOF_IN_NODE(ACOUSTIC_PRESSURE,rnode)
      }

    return ier;

    KRATOS_CATCH( "" );
}


/***********************************************************************************/
/***********************************************************************************/

bool AcousticElement::UseElementProvidedStrain() const
{
    return true;
}

/***********************************************************************************/
/***********************************************************************************/

void AcousticElement::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo
    )
{
    //calculation flags
    const bool CalculateStiffnessMatrixFlag = true;
    const bool CalculateResidualVectorFlag = true;

    CalculateAll( rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag );
}

/***********************************************************************************/
/***********************************************************************************/

void AcousticElement::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
                                             ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;
    VectorType RHS;
    CalculateLocalSystem(rLeftHandSideMatrix, RHS, rCurrentProcessInfo);
    
    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void AcousticElement::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo
    )
{
    // Calculation flags
    const bool CalculateStiffnessMatrixFlag = false;
    const bool CalculateResidualVectorFlag = true;
    MatrixType temp = Matrix();

    CalculateAll( temp, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag );
}

/***********************************************************************************/
/***********************************************************************************/
void AcousticElement::CalculateAll(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo,
    const bool CalculateStiffnessMatrixFlag, 
    const bool CalculateResidualVectorFlag
    )
{
    KRATOS_TRY;

    auto& r_geometry = this->GetGeometry();
    const SizeType number_of_nodes = r_geometry.size();
    const SizeType dimension = r_geometry.WorkingSpaceDimension();
    const SizeType pressuregradient_size = number_of_nodes * dimension;

    KinematicVariables this_kinematic_variables(pressuregradient_size, dimension, number_of_nodes);
    ConstitutiveVariables this_constitutive_variables;

    // Resizing as needed the LHS
    const SizeType mat_size = number_of_nodes;

    if ( CalculateStiffnessMatrixFlag ) { // Calculation of the matrix is required
        if ( rLeftHandSideMatrix.size1() != mat_size )
            rLeftHandSideMatrix.resize( mat_size, mat_size, false );
        noalias( rLeftHandSideMatrix ) = ZeroMatrix( mat_size, mat_size ); 
    }

    // Resizing as needed the RHS
    if ( CalculateResidualVectorFlag ) { // Calculation of the matrix is required
        if ( rRightHandSideVector.size() != mat_size )
            rRightHandSideVector.resize( mat_size, false );

        rRightHandSideVector = ZeroVector( mat_size ); //resetting RHS
    }

    // Reading integration points and local gradients
    const GeometryType::IntegrationPointsArrayType& integration_points = r_geometry.IntegrationPoints(this->GetIntegrationMethod());

    ConstitutiveLaw::Parameters Values(r_geometry,GetProperties(),rCurrentProcessInfo);

    // Some declarations
    double body_force;
    double int_to_reference_weight;

    // Computing in all integrations points
    for ( IndexType point_number = 0; point_number < integration_points.size(); ++point_number ) {
        // Contribution to external forces

        body_force = this->GetBodyForce(integration_points, point_number);
        //Get acoustic source vector --> Assumed as zero - Rigid wall conditions
        // Compute element kinematics B, F, DN_DX ...
        CalculateKinematicVariables(this_kinematic_variables, point_number, this->GetIntegrationMethod());

        // Calculating weights for integration on the reference configuration
        int_to_reference_weight = GetIntegrationWeight(integration_points, point_number, this_kinematic_variables.detJ0);

        if ( CalculateStiffnessMatrixFlag ) { // Calculation of the matrix is required
            // Contributions to stiffness matrix calculated on the reference config
            CalculateAndAddKm( rLeftHandSideMatrix, this_kinematic_variables.B,  int_to_reference_weight );
        }

        if ( CalculateResidualVectorFlag ) { // Calculation of the matrix is required
            CalculateAndAddResidualVector(rRightHandSideVector, this_kinematic_variables, rCurrentProcessInfo, body_force,this_kinematic_variables.PressureGradient, int_to_reference_weight);
        }
    }
    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/


void AcousticElement::CalculateKinematicVariables(
    KinematicVariables& rThisKinematicVariables,
    const IndexType PointNumber,
    const GeometryType::IntegrationMethod& rIntegrationMethod
    )
{

    const auto& r_geometry = GetGeometry();

    const GeometryType::IntegrationPointsArrayType& r_integration_points = r_geometry.IntegrationPoints(rIntegrationMethod);
    // Shape functions
    rThisKinematicVariables.N = r_geometry.ShapeFunctionsValues(rThisKinematicVariables.N, r_integration_points[PointNumber].Coordinates());

    rThisKinematicVariables.detJ0 = CalculateDerivativesOnReferenceConfiguration(rThisKinematicVariables.J0, rThisKinematicVariables.InvJ0, rThisKinematicVariables.DN_DX, PointNumber, rIntegrationMethod);

    KRATOS_ERROR_IF(rThisKinematicVariables.detJ0 < 0.0) << "WARNING:: ELEMENT ID: " << this->Id() << " INVERTED. DETJ0: " << rThisKinematicVariables.detJ0 << std::endl;

    // Compute B
    CalculateB( rThisKinematicVariables.B, rThisKinematicVariables.DN_DX, r_integration_points, PointNumber );

    // Compute equivalent F
    GetValuesVector(rThisKinematicVariables.Acoustic_pressure);

}



// /***********************************************************************************/
// /***********************************************************************************/

 void AcousticElement::SetConstitutiveVariables(
        ConstitutiveVariables& rThisConstitutiveVariables,
        ConstitutiveLaw::Parameters& rValues
        )
{
    // Here we set the space on which the results shall be written
    rValues.SetConstitutiveMatrix(rThisConstitutiveVariables.D); //assuming the determinant is computed somewhere else
   
}

// /***********************************************************************************/
// /***********************************************************************************/

void AcousticElement::CalculateConstitutiveVariables(
    ConstitutiveVariables& rThisConstitutiveVariables,
    ConstitutiveLaw::Parameters& rValues
    )
{
    // Set the constitutive variables
    SetConstitutiveVariables(rThisConstitutiveVariables, rValues);

}


/***********************************************************************************/
/***********************************************************************************/


double  AcousticElement::GetBodyForce(
    const GeometryType::IntegrationPointsArrayType& rIntegrationPoints,
    const IndexType PointNumber
    )
{
    double body_force = 0.0;
    return body_force;
}


/***********************************************************************************/
/***********************************************************************************/


void AcousticElement::CalculateAndAddKm(
    MatrixType& rLeftHandSideMatrix,
    const Matrix& B,
    const double IntegrationWeight
    ) const
{
    KRATOS_TRY

    noalias( rLeftHandSideMatrix ) += IntegrationWeight * prod( trans( B ),  B); 

    KRATOS_CATCH( "" )
}


/***********************************************************************************/
/***********************************************************************************/

void AcousticElement::CalculateMassMatrix(
    MatrixType& rMassMatrix,
    ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    const auto& r_prop = GetProperties();

    const auto& r_geom = GetGeometry();
    SizeType number_of_nodes = r_geom.size();
    SizeType mat_size = number_of_nodes;

    // Clear matrix
    if (rMassMatrix.size1() != mat_size || rMassMatrix.size2() != mat_size)
        rMassMatrix.resize( mat_size, mat_size, false );
    rMassMatrix = ZeroMatrix( mat_size, mat_size );

    // Checking density
    KRATOS_ERROR_IF_NOT(r_prop.Has(DENSITY)) << "DENSITY has to be provided for the calculation of the MassMatrix!" << std::endl;

    // LUMPED MASS MATRIX
    VectorType temp_vector(mat_size);
    CalculateLumpedMassVector(temp_vector);
    for (IndexType i = 0; i < mat_size; ++i)
        rMassMatrix(i, i) = temp_vector[i];

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/
void AcousticElement::CalculateLumpedMassVector(VectorType& rMassVector) const
{
    KRATOS_TRY;

    const auto& r_geom = GetGeometry();
    const auto& r_prop = GetProperties();
    const SizeType dimension = r_geom.WorkingSpaceDimension();
    const SizeType number_of_nodes = r_geom.size();
    const SizeType mat_size = number_of_nodes;

    // Clear matrix
    if (rMassVector.size() != mat_size)
        rMassVector.resize( mat_size, false );

    // Getting the density & Young's modulus of the element
    const double density = r_prop[DENSITY];
    const double young_modulus = r_prop[YOUNG_MODULUS];

    const double thickness = (dimension == 2 && r_prop.Has(THICKNESS)) ? r_prop[THICKNESS] : 1.0;

    // LUMPED MASS MATRIX
    const double total_mass = GetGeometry().DomainSize() * density / young_modulus;

    Vector lumping_factors;
    lumping_factors = GetGeometry().LumpingFactors( lumping_factors );

    for ( IndexType i = 0; i < number_of_nodes; ++i ) {
        const double temp = lumping_factors[i] * total_mass;
        rMassVector[i] = temp;
    }

    KRATOS_CATCH("");
}


/***********************************************************************************/
/***********************************************************************************/


void AcousticElement::CalculateAndAddResidualVector(
    VectorType& rRightHandSideVector,
    const KinematicVariables& rThisKinematicVariables,
    const ProcessInfo& rCurrentProcessInfo,
    const double rBodyForce,
    const Vector& rPressureGradientVector,
    const double IntegrationWeight
    ) const
{
    KRATOS_TRY

    // Operation performed: rRightHandSideVector += ExtForce * IntegrationWeight
    this->CalculateAndAddExtForceContribution( rThisKinematicVariables.N, rCurrentProcessInfo, rBodyForce, rRightHandSideVector, IntegrationWeight );

    // Operation performed: rRightHandSideVector -= IntForce * IntegrationWeight
    noalias( rRightHandSideVector ) -= IntegrationWeight * prod( trans( rThisKinematicVariables.B ), rPressureGradientVector);

    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

void AcousticElement::CalculateB(
    Matrix& rB,
    const Matrix& rDN_DX,
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints,
    const IndexType PointNumber
    ) const
{
    KRATOS_TRY;

    const auto& r_geometry = this->GetGeometry();
    const SizeType number_of_nodes = r_geometry.PointsNumber();
    const SizeType dimension = r_geometry.WorkingSpaceDimension();

    rB.clear();


    if(dimension == 2) {
        for ( IndexType i = 0; i < number_of_nodes; ++i ) {
            const IndexType initial_index = i;
            rB(0, initial_index ) = rDN_DX(i, 0);
            rB(1, initial_index ) = rDN_DX(i, 1);
            
        }
    } else if(dimension == 3) {
        for ( IndexType i = 0; i < number_of_nodes; ++i ) {
            const IndexType initial_index = i;
            rB(0, initial_index ) = rDN_DX(i, 0);
            rB(1, initial_index ) = rDN_DX(i, 1);
            rB(2, initial_index ) = rDN_DX(i, 2);
  
        }
    }

    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

double AcousticElement::CalculateDerivativesOnReferenceConfiguration(
    Matrix& rJ0,
    Matrix& rInvJ0,
    Matrix& rDN_DX,
    const IndexType PointNumber,
    IntegrationMethod ThisIntegrationMethod
    ) const
{
    const GeometryType& r_geom = GetGeometry();
    GeometryUtils::JacobianOnInitialConfiguration(
        r_geom,
        r_geom.IntegrationPoints(ThisIntegrationMethod)[PointNumber], rJ0);
    double detJ0;
    MathUtils<double>::InvertMatrix(rJ0, rInvJ0, detJ0);
    const Matrix& rDN_De =
        GetGeometry().ShapeFunctionsLocalGradients(ThisIntegrationMethod)[PointNumber];
    GeometryUtils::ShapeFunctionsGradients(rDN_De, rInvJ0, rDN_DX);
    return detJ0;
}


double AcousticElement::GetIntegrationWeight(
    const GeometryType::IntegrationPointsArrayType& rThisIntegrationPoints,
    const IndexType PointNumber,
    const double detJ
    ) const
{
    return rThisIntegrationPoints[PointNumber].Weight() * detJ;
}



/***********************************************************************************/
/***********************************************************************************/

void AcousticElement::CalculateAndAddExtForceContribution(
    const Vector& rN,
    const ProcessInfo& rCurrentProcessInfo,
    const double rBodyForce,
    VectorType& rRightHandSideVector,
    const double Weight
    ) const
{
    KRATOS_TRY;

    const SizeType number_of_nodes = GetGeometry().PointsNumber();
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();

    for ( IndexType i = 0; i < number_of_nodes; ++i ) {
        const SizeType index = i;
        rRightHandSideVector[index] += Weight * rN[index] * rBodyForce;
    }

    KRATOS_CATCH( "" )
}


/***********************************************************************************/
/***********************************************************************************/


void AcousticElement::EquationIdVector(
    EquationIdVectorType& rResult,
    ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    const SizeType number_of_nodes = GetGeometry().size();

    KRATOS_WATCH(number_of_nodes);

    if (rResult.size() != number_of_nodes)
        rResult.resize(number_of_nodes,false);

    const SizeType pos = this->GetGeometry()[0].GetDofPosition(ACOUSTIC_PRESSURE);

    for (IndexType i = 0; i < number_of_nodes; ++i) {
        rResult[i] = GetGeometry()[i].GetDof(ACOUSTIC_PRESSURE,pos).EquationId();
       }

    KRATOS_CATCH("")
};

/***********************************************************************************/
/***********************************************************************************/

void AcousticElement::GetDofList(
    DofsVectorType& rElementalDofList,
    ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    const SizeType number_of_nodes = GetGeometry().size();
    rElementalDofList.resize(0);
    rElementalDofList.reserve(number_of_nodes);

    for (IndexType i = 0; i < number_of_nodes; ++i) {
        rElementalDofList.push_back(GetGeometry()[i].pGetDof(ACOUSTIC_PRESSURE));
    }

    KRATOS_CATCH("")
};

/***********************************************************************************/
/***********************************************************************************/




} // Namespace Kratos


