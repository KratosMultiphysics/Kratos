// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Athira Vadakkekkara
//

// System includes

// External includes


// Project includes
#include "utilities/math_utils.h"

// Application includes
#include "custom_elements/small_displacement_hex_two_non_local_variables.h"
#include "custom_utilities/structural_mechanics_element_utilities.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos
{
SmallDisplacementHexTwoNonLocalVariables::SmallDisplacementHexTwoNonLocalVariables( IndexType NewId, GeometryType::Pointer pGeometry )
    : BaseSolidElement( NewId, pGeometry )
{
    //DO NOT ADD DOFS HERE!!!
    mConstitutiveMatrix.resize(8,8);
}

/***********************************************************************************/
/***********************************************************************************/

SmallDisplacementHexTwoNonLocalVariables::SmallDisplacementHexTwoNonLocalVariables( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
        : BaseSolidElement( NewId, pGeometry, pProperties )
{
    //DO NOT ADD DOFS HERE!!!
    mConstitutiveMatrix.resize(8,8);
}

/***********************************************************************************/
/***********************************************************************************/

Element::Pointer SmallDisplacementHexTwoNonLocalVariables::Create( IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties ) const
{
    return Kratos::make_intrusive<SmallDisplacementHexTwoNonLocalVariables>( NewId, GetGeometry().Create( ThisNodes ), pProperties );
}

/***********************************************************************************/
/***********************************************************************************/

Element::Pointer SmallDisplacementHexTwoNonLocalVariables::Create( IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties ) const
{
    return Kratos::make_intrusive<SmallDisplacementHexTwoNonLocalVariables>( NewId, pGeom, pProperties );
}

/***********************************************************************************/
/***********************************************************************************/

SmallDisplacementHexTwoNonLocalVariables::~SmallDisplacementHexTwoNonLocalVariables()
{
}

/***********************************************************************************/
/***********************************************************************************/

Element::Pointer SmallDisplacementHexTwoNonLocalVariables::Clone (
    IndexType NewId,
    NodesArrayType const& rThisNodes
    ) const
{
    KRATOS_TRY

    SmallDisplacementHexTwoNonLocalVariables::Pointer p_new_elem = Kratos::make_intrusive<SmallDisplacementHexTwoNonLocalVariables>(NewId, GetGeometry().Create(rThisNodes), pGetProperties());
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

void SmallDisplacementHexTwoNonLocalVariables::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    BaseType::Initialize(rCurrentProcessInfo);
    mpEquivalentStrainsTC = &KratosComponents<Variable<Vector>>::Get(this->GetProperties()[NON_LOCAL_VARIABLE_NAME]);
}

/***********************************************************************************/
/***********************************************************************************/

bool SmallDisplacementHexTwoNonLocalVariables::UseElementProvidedStrain() const
{
    return true;
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementHexTwoNonLocalVariables::EquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo
    ) const
{
    KRATOS_TRY;

    auto& r_geometry = this->GetGeometry();
    const SizeType number_of_nodes = r_geometry.size();
    const SizeType dimension = r_geometry.WorkingSpaceDimension();
    const SizeType mat_size = number_of_nodes * (dimension + 2);

    if (rResult.size() != mat_size){
        rResult.resize(mat_size,false);
    }

    const SizeType pos = r_geometry[0].GetDofPosition(DISPLACEMENT_X);

    if(dimension == 2) {
        for (IndexType i = 0; i < number_of_nodes; ++i) {
            const SizeType index = i * (dimension + 2);
            rResult[index] = r_geometry[i].GetDof(DISPLACEMENT_X,pos).EquationId();
            rResult[index + 1] = r_geometry[i].GetDof(DISPLACEMENT_Y,pos+1).EquationId();
            rResult[index + 2] = r_geometry[i].GetDof(NON_LOCAL_EQUIVALENT_STRAIN_TENSION).EquationId();
            rResult[index + 3] = r_geometry[i].GetDof(NON_LOCAL_EQUIVALENT_STRAIN_COMPRESSION).EquationId();
        }
    } else {
        for (IndexType i = 0; i < number_of_nodes; ++i) {
            const SizeType index = i * (dimension + 2);
            rResult[index] = r_geometry[i].GetDof(DISPLACEMENT_X,pos).EquationId();
            rResult[index + 1] = r_geometry[i].GetDof(DISPLACEMENT_Y,pos+1).EquationId();
            rResult[index + 2] = r_geometry[i].GetDof(DISPLACEMENT_Z,pos+2).EquationId();
            rResult[index + 3] = r_geometry[i].GetDof(NON_LOCAL_EQUIVALENT_STRAIN_TENSION).EquationId();
            rResult[index + 4] = r_geometry[i].GetDof(NON_LOCAL_EQUIVALENT_STRAIN_COMPRESSION).EquationId();
        }
    }

    KRATOS_CATCH("")
};

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementHexTwoNonLocalVariables::GetDofList(
    DofsVectorType& rElementalDofList,
    const ProcessInfo& rCurrentProcessInfo
    ) const
{
    KRATOS_TRY;

    auto& r_geometry = this->GetGeometry();
    const SizeType number_of_nodes = r_geometry.size();
    const SizeType dimension = r_geometry.WorkingSpaceDimension();
    const SizeType mat_size = number_of_nodes * (dimension + 2);
    rElementalDofList.clear();
    rElementalDofList.reserve(mat_size);

    if(dimension == 2) {
        for (IndexType i = 0; i < number_of_nodes; ++i) {
            rElementalDofList.push_back(r_geometry[i].pGetDof(DISPLACEMENT_X));
            rElementalDofList.push_back( r_geometry[i].pGetDof(DISPLACEMENT_Y));
            rElementalDofList.push_back( r_geometry[i].pGetDof(NON_LOCAL_EQUIVALENT_STRAIN_TENSION));
            rElementalDofList.push_back( r_geometry[i].pGetDof(NON_LOCAL_EQUIVALENT_STRAIN_COMPRESSION));
        }
    } else {
        for (IndexType i = 0; i < number_of_nodes; ++i) {
            rElementalDofList.push_back(r_geometry[i].pGetDof(DISPLACEMENT_X));
            rElementalDofList.push_back( r_geometry[i].pGetDof(DISPLACEMENT_Y));
            rElementalDofList.push_back( r_geometry[i].pGetDof(DISPLACEMENT_Z));
            rElementalDofList.push_back( r_geometry[i].pGetDof(NON_LOCAL_EQUIVALENT_STRAIN_TENSION));
            rElementalDofList.push_back( r_geometry[i].pGetDof(NON_LOCAL_EQUIVALENT_STRAIN_COMPRESSION));
        }
    }

    KRATOS_CATCH("")
};

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementHexTwoNonLocalVariables::GetValuesVector(
    Vector& rValues,
    int Step
    ) const
{
    KRATOS_TRY
    auto& r_geometry = this->GetGeometry();
    const SizeType number_of_nodes = r_geometry.size();
    const SizeType dimension = r_geometry.WorkingSpaceDimension();
    const SizeType mat_size = number_of_nodes * (dimension + 2);
    if (rValues.size() != mat_size) {
        rValues.resize(mat_size, false);
    }

    for (IndexType i = 0; i < number_of_nodes; ++i){
        const array_1d<double, 3 >& displacement = r_geometry[i].FastGetSolutionStepValue(DISPLACEMENT, Step);
        const SizeType index = i * (dimension +2);
        for(unsigned int k = 0; k < dimension; ++k)
        {
            rValues[index + k] = displacement[k];
        }
        rValues[index + dimension] = r_geometry[i].FastGetSolutionStepValue(NON_LOCAL_EQUIVALENT_STRAIN_TENSION, Step);
        rValues[index + dimension + 1] = r_geometry[i].FastGetSolutionStepValue(NON_LOCAL_EQUIVALENT_STRAIN_COMPRESSION, Step);
    }
    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementHexTwoNonLocalVariables::GetFirstDerivativesVector(
    Vector& rValues,
    int Step
    ) const
{
    KRATOS_TRY
    auto& r_geometry = this->GetGeometry();
    const SizeType number_of_nodes = r_geometry.size();
    const SizeType dimension = r_geometry.WorkingSpaceDimension();
    const SizeType mat_size = number_of_nodes * (dimension + 2);
    if (rValues.size() != mat_size){
        rValues.resize(mat_size, false);
    }
    for (IndexType i = 0; i < number_of_nodes; ++i) {
        const array_1d<double, 3 >& velocity = r_geometry[i].FastGetSolutionStepValue(VELOCITY, Step);
        const SizeType index = i * (dimension + 2);
        for(unsigned int k = 0; k < dimension; ++k){
            rValues[index + k] = velocity[k];
        }
        rValues[index + dimension] = 0.0;
        rValues[index + dimension + 1] = 0.0;
    }
    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementHexTwoNonLocalVariables::GetSecondDerivativesVector(
    Vector& rValues,
    int Step
    ) const
{
    KRATOS_TRY
    auto& r_geometry = this->GetGeometry();
    const SizeType number_of_nodes = r_geometry.size();
    const SizeType dimension = r_geometry.WorkingSpaceDimension();
    const SizeType mat_size = number_of_nodes * (dimension+2);
    if (rValues.size() != mat_size){
        rValues.resize(mat_size, false);
    }
    for (IndexType i = 0; i < number_of_nodes; ++i) {
        const array_1d<double, 3 >& acceleration = r_geometry[i].FastGetSolutionStepValue(ACCELERATION, Step);
        const SizeType index = i * (dimension + 2);
        for(unsigned int k = 0; k < dimension; ++k){
            rValues[index + k] = acceleration[k];
        }
        rValues[index + dimension] = 0.0;
        rValues[index + dimension + 1] = 0.0;
    }
    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementHexTwoNonLocalVariables::CalculateAll(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo,
    const bool CalculateStiffnessMatrixFlag,
    const bool CalculateResidualVectorFlag
    )
{
    KRATOS_TRY;

    auto& r_geometry = this->GetGeometry();
    const SizeType number_of_nodes = GetGeometry().size();
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();
    const SizeType strain_size = (dimension == 3) ? 6 : 4;

    KinematicVariables this_kinematic_variables(strain_size, dimension, number_of_nodes);
    ConstitutiveVariables this_constitutive_variables(strain_size);
    NonLocalConstitutiveVariables this_non_local_constitutive_variables(strain_size);

    // Resizing as needed the LHS
    const SizeType mat_size = number_of_nodes * (dimension + 2);

    if ( CalculateStiffnessMatrixFlag ) { // Calculation of the matrix is required
        if ( rLeftHandSideMatrix.size1() != mat_size || rLeftHandSideMatrix.size2() != mat_size )
            rLeftHandSideMatrix.resize( mat_size, mat_size, false );

        noalias( rLeftHandSideMatrix ) = ZeroMatrix( mat_size, mat_size ); //resetting LHS
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
    Matrix Kuu = ZeroMatrix(number_of_nodes * dimension , number_of_nodes * dimension);
    Matrix KuNL1 = ZeroMatrix(number_of_nodes * dimension , number_of_nodes);
    Matrix KuNL2 = ZeroMatrix(number_of_nodes * dimension , number_of_nodes);
    Matrix KNL1u = ZeroMatrix(number_of_nodes, number_of_nodes * dimension);
    Matrix KNL2u = ZeroMatrix(number_of_nodes, number_of_nodes * dimension);
    Matrix KNL1NL1 = ZeroMatrix(number_of_nodes , number_of_nodes);
    Matrix KNL2NL2 = ZeroMatrix(number_of_nodes , number_of_nodes);
    Vector Fu = ZeroVector(number_of_nodes * dimension);
    Vector FNL1 = ZeroVector(number_of_nodes);
    Vector FNL2 = ZeroVector(number_of_nodes);

    // Computing in all integrations points
    for ( IndexType point_number = 0; point_number < integration_points.size(); ++point_number ) {
        // Contribution to external forces
        noalias(body_force) = this->GetBodyForce(integration_points, point_number);

        // Compute element kinematics B, F, DN_DX ...
        CalculateKinematicVariables(this_kinematic_variables, point_number, this->GetIntegrationMethod());

        // Compute material reponse
        CalculateAllConstitutiveVariables(this_kinematic_variables, this_constitutive_variables, this_non_local_constitutive_variables, Values, point_number, integration_points, GetStressMeasure());

        // Calculate and set non locl constitutive variables......................................................
        CalculateNonLocalConstitutiveVariables(this_kinematic_variables, this_non_local_constitutive_variables, Values, *mpEquivalentStrainsTC, point_number);

        SetNonLocalConstitutiveVariables(mConstitutiveMatrix, this_constitutive_variables, this_non_local_constitutive_variables);

        // Calculating weights for integration on the reference configuration
        int_to_reference_weight = GetIntegrationWeight(integration_points, point_number, this_kinematic_variables.detJ0);

        if ( dimension == 2 && GetProperties().Has( THICKNESS ))
            int_to_reference_weight *= GetProperties()[THICKNESS];

        if ( CalculateStiffnessMatrixFlag ) { // Calculation of the matrix is required
            // Contributions to stiffness matrix calculated on the reference config
            CalculateAndAddKuu(Kuu, this_kinematic_variables.B, this_constitutive_variables.D, int_to_reference_weight );
            CalculateAndAddKuNL(KuNL1, KuNL2, this_non_local_constitutive_variables, this_kinematic_variables.B, this_kinematic_variables.N, int_to_reference_weight);
            CalculateAndAddKNLu(KNL1u, KNL2u, this_non_local_constitutive_variables, this_kinematic_variables.B, this_kinematic_variables.N, int_to_reference_weight);
            CalculateAndAddKNLNL(KNL1NL1, KNL2NL2, this_non_local_constitutive_variables, this_kinematic_variables.DN_DX, this_kinematic_variables.N, int_to_reference_weight);
        }
        if ( CalculateResidualVectorFlag ) { // Calculation of the matrix is required
            CalculateAndAddResidualForceVector(Fu, this_kinematic_variables, rCurrentProcessInfo, body_force, this_constitutive_variables.StressVector, int_to_reference_weight);
            CalculateAndAddResidualNonLocalVector(FNL1, FNL2, this_kinematic_variables, this_non_local_constitutive_variables, rCurrentProcessInfo, int_to_reference_weight);
        }
    }
    //assemble right hand side and left hand sides
    AssembleRHSAndLHS(rLeftHandSideMatrix, rRightHandSideVector, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag, Kuu, KuNL1, KuNL2, KNL1u, KNL2u, KNL1NL1, KNL2NL2, Fu, FNL1, FNL2);
    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementHexTwoNonLocalVariables::FinalizeSolutionStep( const ProcessInfo& rCurrentProcessInfo )
{
    // We finalize the material reponse if required
    bool required = false;
    for ( IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number ) {
        if (mConstitutiveLawVector[point_number]->RequiresFinalizeMaterialResponse()) {
            required = true;
            break;
        }
    }
    if (required) {
        const GeometryType& r_geometry = GetGeometry();
        const Properties& r_properties = GetProperties();
        const SizeType number_of_nodes = r_geometry.size();
        const SizeType dimension = r_geometry.WorkingSpaceDimension();
        const SizeType strain_size = mConstitutiveLawVector[0]->GetStrainSize();

        KinematicVariables this_kinematic_variables(strain_size, dimension, number_of_nodes);
        ConstitutiveVariables this_constitutive_variables(strain_size);
        NonLocalConstitutiveVariables this_non_local_constitutive_variables(strain_size);

        // Create constitutive law parameters:
        ConstitutiveLaw::Parameters Values(r_geometry, r_properties, rCurrentProcessInfo);

        // Set constitutive law flags:
        Flags& ConstitutiveLawOptions=Values.GetOptions();
        ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, UseElementProvidedStrain());
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);

        Values.SetStrainVector(this_constitutive_variables.StrainVector);
        Values.SetStressVector(this_constitutive_variables.StressVector);
        Values.SetConstitutiveMatrix(mConstitutiveMatrix);
        SetNonLocalConstitutiveVariables(mConstitutiveMatrix, this_constitutive_variables, this_non_local_constitutive_variables);

        // Reading integration points
        const auto& N_values = r_geometry.ShapeFunctionsValues(mThisIntegrationMethod);

        // Reading integration points
        const GeometryType::IntegrationPointsArrayType& integration_points = r_geometry.IntegrationPoints(mThisIntegrationMethod);

        for ( IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number ) {
            if (mConstitutiveLawVector[point_number]->RequiresFinalizeMaterialResponse()) {
                // Compute element kinematics B, F, DN_DX ...
                CalculateKinematicVariables(this_kinematic_variables, point_number, mThisIntegrationMethod);

                // Compute constitutive law variables
                SetConstitutiveVariables(this_kinematic_variables, this_constitutive_variables, this_non_local_constitutive_variables, Values, point_number, integration_points);
                Vector r_stress_vector = this_constitutive_variables.StressVector;
                this->SetValue(STRESSES, r_stress_vector);
                // Call the constitutive law to update material variables
                mConstitutiveLawVector[point_number]->FinalizeMaterialResponse(Values, GetStressMeasure());

                // TODO: Deprecated, remove this
                mConstitutiveLawVector[point_number]->FinalizeSolutionStep( r_properties, r_geometry, row( N_values, point_number ), rCurrentProcessInfo);

            }
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementHexTwoNonLocalVariables::CalculateKinematicVariables(
    KinematicVariables& rThisKinematicVariables,
    const IndexType PointNumber,
    const GeometryType::IntegrationMethod& rIntegrationMethod
    )
{
    KRATOS_TRY;
    const auto& r_geometry = GetGeometry();

    const GeometryType::IntegrationPointsArrayType& r_integration_points = r_geometry.IntegrationPoints(rIntegrationMethod);
    // Shape functions
    rThisKinematicVariables.N = r_geometry.ShapeFunctionsValues(rThisKinematicVariables.N, r_integration_points[PointNumber].Coordinates());

    rThisKinematicVariables.detJ0 = CalculateDerivativesOnReferenceConfiguration(rThisKinematicVariables.J0, rThisKinematicVariables.InvJ0, rThisKinematicVariables.DN_DX, PointNumber, rIntegrationMethod);

    KRATOS_ERROR_IF(rThisKinematicVariables.detJ0 < 0.0) << "WARNING:: ELEMENT ID: " << this->Id() << " INVERTED. DETJ0: " << rThisKinematicVariables.detJ0 << std::endl;

    // Compute B
    CalculateB( rThisKinematicVariables.B, rThisKinematicVariables.DN_DX, r_integration_points, PointNumber );

    // Compute equivalent F
    BaseSolidElement::GetValuesVector(rThisKinematicVariables.Displacements);
    Vector strain_vector = prod(rThisKinematicVariables.B, rThisKinematicVariables.Displacements);
    ComputeEquivalentF(rThisKinematicVariables.F, strain_vector);
    rThisKinematicVariables.detF = MathUtils<double>::Det(rThisKinematicVariables.F);

    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementHexTwoNonLocalVariables::CalculateAllConstitutiveVariables(
    KinematicVariables& rThisKinematicVariables,
    ConstitutiveVariables& rThisConstitutiveVariables,
    NonLocalConstitutiveVariables& rThisNonLocalConstitutiveVariables,
    ConstitutiveLaw::Parameters& rValues,
    const IndexType PointNumber,
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints,
    const ConstitutiveLaw::StressMeasure ThisStressMeasure
    )
{
    KRATOS_TRY;
    // Set the constitutive variables
    SetConstitutiveVariables(rThisKinematicVariables, rThisConstitutiveVariables, rThisNonLocalConstitutiveVariables, rValues, PointNumber, IntegrationPoints);
    // Actually do the computations in the ConstitutiveLaw
    mConstitutiveLawVector[PointNumber]->CalculateMaterialResponse(rValues, ThisStressMeasure); //here the calculations are actually done
    //NOTE: setting NonLocalConstitutiveVariables from the constitutivematrix(7,7)
    SetNonLocalConstitutiveVariables(mConstitutiveMatrix, rThisConstitutiveVariables, rThisNonLocalConstitutiveVariables);

    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementHexTwoNonLocalVariables::SetConstitutiveVariables(
    KinematicVariables& rThisKinematicVariables,
    ConstitutiveVariables& rThisConstitutiveVariables,
    NonLocalConstitutiveVariables& rThisNonLocalConstitutiveVariables,
    ConstitutiveLaw::Parameters& rValues,
    const IndexType PointNumber,
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints
    )
{
    KRATOS_TRY;
    const auto& r_geometry = GetGeometry();
    const SizeType number_of_nodes = r_geometry.size();
    const SizeType dimension = r_geometry.WorkingSpaceDimension();

    //displacement vector
    Vector displacements(number_of_nodes * dimension);
    BaseSolidElement::GetValuesVector(displacements,0);

    // Compute strain
    noalias(rThisConstitutiveVariables.StrainVector) = prod(rThisKinematicVariables.B, displacements);

    // Here we essentially set the input parameters
    rValues.SetShapeFunctionsValues(rThisKinematicVariables.N); // shape functions
    rValues.SetDeterminantF(rThisKinematicVariables.detF); //assuming the determinant is computed somewhere else
    rValues.SetDeformationGradientF(rThisKinematicVariables.F); //F computed somewhere else

    // Here we set the space on which the results shall be written
    rValues.SetConstitutiveMatrix(mConstitutiveMatrix);
    rValues.SetStressVector(rThisConstitutiveVariables.StressVector);

    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementHexTwoNonLocalVariables::CalculateNonLocalConstitutiveVariables(
    KinematicVariables& rThisKinematicVariables,
    NonLocalConstitutiveVariables& rThisNonLocalConstitutiveVariables,
    ConstitutiveLaw::Parameters& rValues,
    const Variable<Vector>& rThisVariable,
    const IndexType PointNumber
    )
{
    KRATOS_TRY;
    const auto& r_geometry = GetGeometry();
    const auto& N = rThisKinematicVariables.N;
    Vector local_variables_gp = ZeroVector(2);
    rThisNonLocalConstitutiveVariables.NonLocal_Variables_GP = ZeroVector(2);
    for (unsigned int i = 0; i < N.size(); ++i) {
        rThisNonLocalConstitutiveVariables.NonLocal_Variables_GP[0] += N[i] * r_geometry[i].FastGetSolutionStepValue(NON_LOCAL_EQUIVALENT_STRAIN_TENSION, 0);
        rThisNonLocalConstitutiveVariables.NonLocal_Variables_GP[1] += N[i] * r_geometry[i].FastGetSolutionStepValue(NON_LOCAL_EQUIVALENT_STRAIN_COMPRESSION, 0);
    }
    mConstitutiveLawVector[PointNumber]->CalculateValue(rValues, rThisVariable, local_variables_gp);
    rThisNonLocalConstitutiveVariables.Local_Variables_GP = local_variables_gp;
    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementHexTwoNonLocalVariables::SetNonLocalConstitutiveVariables(
    const Matrix& ConstitutiveMatrix,
    ConstitutiveVariables& rThisConstitutiveVariables,
    NonLocalConstitutiveVariables& rThisNonLocalConstitutiveVariables
    )const
{
    KRATOS_TRY;

    const SizeType strain_size = mConstitutiveLawVector[0]->GetStrainSize();
    for (std::size_t i = 0; i < strain_size ; ++i) {
        for (std::size_t j = 0; j < strain_size ; ++j) {
            rThisConstitutiveVariables.D(i,j) = ConstitutiveMatrix(i, j) ;
        }
    }
    for (std::size_t i = 0; i < strain_size; ++i) {
        rThisNonLocalConstitutiveVariables.DuNL1[i] = ConstitutiveMatrix(i, strain_size) ;
        rThisNonLocalConstitutiveVariables.DuNL2[i] = ConstitutiveMatrix(i, strain_size+1) ;
    }
    for (std::size_t i = 0; i < strain_size; ++i) {
        rThisNonLocalConstitutiveVariables.DNL1u[i] =  ConstitutiveMatrix(strain_size, i);
        rThisNonLocalConstitutiveVariables.DNL2u[i] =  ConstitutiveMatrix(strain_size+1, i);
    }
    rThisNonLocalConstitutiveVariables.DNL1NL1 =  ConstitutiveMatrix(strain_size, strain_size);
    rThisNonLocalConstitutiveVariables.DNL2NL2 =  ConstitutiveMatrix(strain_size+1, strain_size+1);

    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementHexTwoNonLocalVariables::CalculateB(
    Matrix& rB,
    const Matrix& rDN_DX,
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints,
    const IndexType PointNumber
    ) const
{
    KRATOS_TRY;

    StructuralMechanicsElementUtilities::CalculateB(*this, rDN_DX, rB);

    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementHexTwoNonLocalVariables::ComputeEquivalentF(
    Matrix& rF,
    const Vector& rStrainTensor
    ) const
{
    StructuralMechanicsElementUtilities::ComputeEquivalentF(*this, rStrainTensor, rF);
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementHexTwoNonLocalVariables::CalculateAndAddKuu(
    MatrixType& rKuu,
    const Matrix& B,
    const Matrix& D,
    const double IntegrationWeight
    ) const
{
    KRATOS_TRY
    rKuu  += IntegrationWeight * prod( trans( B ), Matrix(prod(D, B)));
    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementHexTwoNonLocalVariables::CalculateAndAddKuNL(
    Matrix& rStiffnessMatrixKuNL1,
    Matrix& rStiffnessMatrixKuNL2,
    NonLocalConstitutiveVariables& rThisNonLocalConstitutiveVariables,
    const Matrix& B,
    const Vector& N,
    const double IntegrationWeight
    ) const
{
    KRATOS_TRY

    rStiffnessMatrixKuNL1 += IntegrationWeight * prod( trans( B ), outer_prod(rThisNonLocalConstitutiveVariables.DuNL1, N));
    rStiffnessMatrixKuNL2 += IntegrationWeight * prod( trans( B ), outer_prod(rThisNonLocalConstitutiveVariables.DuNL2, N));

    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementHexTwoNonLocalVariables::CalculateAndAddKNLu(
    Matrix& rStiffnessMatrixKNL1u,
    Matrix& rStiffnessMatrixKNL2u,
    NonLocalConstitutiveVariables& rThisNonLocalConstitutiveVariables,
    const Matrix& B,
    const Vector& N,
    const double IntegrationWeight
    ) const
{
    KRATOS_TRY

    rStiffnessMatrixKNL1u  -= IntegrationWeight * prod( outer_prod(N, rThisNonLocalConstitutiveVariables.DNL1u), B);
    rStiffnessMatrixKNL2u  -= IntegrationWeight * prod( outer_prod(N, rThisNonLocalConstitutiveVariables.DNL2u), B);

    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementHexTwoNonLocalVariables::CalculateAndAddKNLNL(
    Matrix& rStiffnessMatrixKNL1NL1,
    Matrix& rStiffnessMatrixKNL2NL2,
    NonLocalConstitutiveVariables& rThisNonLocalConstitutiveVariables,
    const Matrix& DN_DX,
    const Vector& N,
    const double IntegrationWeight
    ) const
{
    KRATOS_TRY
    const Properties& r_properties = GetProperties();
    const double Lc = r_properties[CHARACTERISTIC_INTERNAL_LENGTH];
    rStiffnessMatrixKNL1NL1 += IntegrationWeight *( (1-rThisNonLocalConstitutiveVariables.DNL1NL1) * outer_prod(N, N) + pow(Lc, 2) * prod( DN_DX , trans(DN_DX)));
    rStiffnessMatrixKNL2NL2 += IntegrationWeight *( (1-rThisNonLocalConstitutiveVariables.DNL2NL2) * outer_prod(N, N) + pow(Lc, 2) * prod( DN_DX , trans(DN_DX)));
    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementHexTwoNonLocalVariables::CalculateAndAddResidualForceVector(
    VectorType& rRightHandSideVector,
    const KinematicVariables& rThisKinematicVariables,
    const ProcessInfo& rCurrentProcessInfo,
    const array_1d<double, 3>& rBodyForce,
    const Vector& rStressVector,
    const double IntegrationWeight
    ) const
{
    KRATOS_TRY

    // Operation performed: rRightHandSideVector += ExtForce * IntegrationWeight
    this->CalculateAndAddExtForceContribution( rThisKinematicVariables.N, rCurrentProcessInfo, rBodyForce, rRightHandSideVector, IntegrationWeight );
    // Operation performed: rRightHandSideVector -= IntForce * IntegrationWeight
    noalias( rRightHandSideVector ) -= IntegrationWeight * prod( trans( rThisKinematicVariables.B ), rStressVector );

    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementHexTwoNonLocalVariables::CalculateAndAddResidualNonLocalVector(
    Vector& rResidualNonLocalVector1,
    Vector& rResidualNonLocalVector2,
    const KinematicVariables& rThisKinematicVariables,
    const NonLocalConstitutiveVariables& rThisNonLocalConstitutiveVariables,
    const ProcessInfo& rCurrentProcessInfo,
    const double IntegrationWeight
    )
{
    KRATOS_TRY;
    auto& r_geometry = this->GetGeometry();
    const Properties& r_properties = GetProperties();
    const double Lc = r_properties[CHARACTERISTIC_INTERNAL_LENGTH];
    const Vector rN = rThisKinematicVariables.N;
    const Matrix rdNdX = rThisKinematicVariables.DN_DX;
    const Vector rNonLocal_Variables_GP = rThisNonLocalConstitutiveVariables.NonLocal_Variables_GP;
    const Vector rLocal_Variables_GP = rThisNonLocalConstitutiveVariables.Local_Variables_GP;
    Vector NL_variable_vector1 = ZeroVector(rN.size());
    Vector NL_variable_vector2 = ZeroVector(rN.size());

    for (unsigned int i = 0; i < rN.size(); ++i) {
        NL_variable_vector1[i] = r_geometry[i].FastGetSolutionStepValue(NON_LOCAL_EQUIVALENT_STRAIN_TENSION, 0);
        NL_variable_vector2[i] = r_geometry[i].FastGetSolutionStepValue(NON_LOCAL_EQUIVALENT_STRAIN_COMPRESSION, 0);
    }
    rResidualNonLocalVector1 -= IntegrationWeight * rN * (rNonLocal_Variables_GP[0]-rLocal_Variables_GP[0]);
    rResidualNonLocalVector1 -= IntegrationWeight * pow(Lc, 2) * prod(prod(rdNdX, trans(rdNdX)),NL_variable_vector1);
    rResidualNonLocalVector2 -= IntegrationWeight * rN * (rNonLocal_Variables_GP[1]-rLocal_Variables_GP[1]);
    rResidualNonLocalVector2 -= IntegrationWeight * pow(Lc, 2) * prod(prod(rdNdX, trans(rdNdX)),NL_variable_vector2);
    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementHexTwoNonLocalVariables::AssembleRHSAndLHS(
    MatrixType& LeftHandSideMatrix,
    VectorType& RightHandSideVector,
    const bool StiffnessMatrixFlag,
    const bool ResidualVectorFlag,
    const Matrix& Kuu,
    const Matrix& KuNL1,
    const Matrix& KuNL2,
    const Matrix& KNL1u,
    const Matrix& KNL2u,
    const Matrix& KNL1NL1,
    const Matrix& KNL2NL2,
    const Vector& rFu,
    const Vector& rFNL1,
    const Vector& rFNL2
    )
{
    KRATOS_TRY;
    auto& r_geometry = this->GetGeometry();
    const SizeType nnod = r_geometry.size();
    const SizeType dimension = r_geometry.WorkingSpaceDimension();
    Vector Pos_u = ZeroVector(nnod * dimension);
    Vector Pos_NL1 = ZeroVector(nnod);
    Vector Pos_NL2 = ZeroVector(nnod);

    for(unsigned int a = 0; a < nnod ; a++){
        const SizeType index_u = a * dimension;
        const SizeType index_NL = a * (dimension + 1);
        for(unsigned int b = 0; b < dimension; b++){
            Pos_u(index_u + b) = index_NL + a + b;
        }
        Pos_NL1(a) = index_NL + dimension + a ;
        Pos_NL2(a) = index_NL + dimension + a + 1;
    }
    for(unsigned int i = 0; i < (nnod * dimension); i++){
        if ( ResidualVectorFlag ) {
            RightHandSideVector(Pos_u[i]) = rFu[i];
            if(i < nnod){
                RightHandSideVector(Pos_NL1[i]) = rFNL1[i];
                RightHandSideVector(Pos_NL2[i]) = rFNL2[i];
            }
        }
        if ( StiffnessMatrixFlag ) {
            for(unsigned int j = 0; j < (nnod * dimension); j++){
                LeftHandSideMatrix(Pos_u[i],Pos_u[j]) = Kuu(i,j);
                if(j < nnod){
                    LeftHandSideMatrix(Pos_u[i],Pos_NL1[j]) = KuNL1(i,j);
                    LeftHandSideMatrix(Pos_u[i],Pos_NL2[j]) = KuNL2(i,j);
                }
                if(i < nnod){
                    LeftHandSideMatrix(Pos_NL1[i],Pos_u[j]) = KNL1u(i,j);
                    LeftHandSideMatrix(Pos_NL2[i],Pos_u[j]) = KNL2u(i,j);
                }
                if(i < nnod && j < nnod){
                    LeftHandSideMatrix(Pos_NL1[i],Pos_NL1[j]) = KNL1NL1(i,j);
                    LeftHandSideMatrix(Pos_NL2[i],Pos_NL2[j]) = KNL2NL2(i,j);
                }
            }
        }
    }
    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementHexTwoNonLocalVariables::CalculateOnIntegrationPoints(
    const Variable<Vector>& rVariable,
    std::vector<Vector>& rOutput,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    auto& r_geometry = this->GetGeometry();
    const GeometryType::IntegrationPointsArrayType& integration_points = r_geometry.IntegrationPoints( this->GetIntegrationMethod() );

    const SizeType number_of_integration_points = integration_points.size();
    if ( rOutput.size() != number_of_integration_points ){
        rOutput.resize( number_of_integration_points );
    }
    if ( rVariable == CAUCHY_STRESS_VECTOR || rVariable == PK2_STRESS_VECTOR ){

        const SizeType number_of_nodes = r_geometry.size();
        const SizeType dimension = r_geometry.WorkingSpaceDimension();
        const SizeType strain_size = mConstitutiveLawVector[0]->GetStrainSize();

        KinematicVariables this_kinematic_variables(strain_size, dimension, number_of_nodes);
        ConstitutiveVariables this_constitutive_variables(strain_size);
        NonLocalConstitutiveVariables this_non_local_constitutive_variables(strain_size);

        // Create constitutive law parameters:
        ConstitutiveLaw::Parameters Values(r_geometry,GetProperties(),rCurrentProcessInfo);

        // Set constitutive law flags:
        Flags& ConstitutiveLawOptions=Values.GetOptions();
        ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, UseElementProvidedStrain());
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);

        Values.SetStrainVector(this_constitutive_variables.StrainVector);

        for (IndexType point_number = 0; point_number < number_of_integration_points; ++point_number) {
            // Compute element kinematics B, F, DN_DX ...
            CalculateKinematicVariables(this_kinematic_variables, point_number, mThisIntegrationMethod);

            // Compute material reponse
            CalculateAllConstitutiveVariables(this_kinematic_variables, this_constitutive_variables, this_non_local_constitutive_variables, Values, point_number, integration_points, GetStressMeasure());
            if ( rOutput[point_number].size() != strain_size ){
                rOutput[point_number].resize( strain_size, false );
            }
            rOutput[point_number] = this_constitutive_variables.StressVector;
        }
    } else if( rVariable == GREEN_LAGRANGE_STRAIN_VECTOR  || rVariable == ALMANSI_STRAIN_VECTOR ) {
            // Create and initialize element variables:
            const SizeType number_of_nodes = r_geometry.size();
            const SizeType dimension = r_geometry.WorkingSpaceDimension();
            const SizeType strain_size = mConstitutiveLawVector[0]->GetStrainSize();

            KinematicVariables this_kinematic_variables(strain_size, dimension, number_of_nodes);
            ConstitutiveVariables this_constitutive_variables(strain_size);
            NonLocalConstitutiveVariables this_non_local_constitutive_variables(strain_size);

            // Create constitutive law parameters:
            ConstitutiveLaw::Parameters Values(r_geometry, GetProperties(), rCurrentProcessInfo);

            // Set constitutive law flags:
            Flags &ConstitutiveLawOptions=Values.GetOptions();
            ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, UseElementProvidedStrain());
            ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, false);
            ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);

            Values.SetStrainVector(this_constitutive_variables.StrainVector);

            const ConstitutiveLaw::StressMeasure this_stress_measure = rVariable == GREEN_LAGRANGE_STRAIN_VECTOR ? ConstitutiveLaw::StressMeasure_PK2 : ConstitutiveLaw::StressMeasure_Kirchhoff;

            //reading integration points
            for ( IndexType point_number = 0; point_number < number_of_integration_points; ++point_number ) {
                // Compute element kinematics B, F, DN_DX ...
                CalculateKinematicVariables(this_kinematic_variables, point_number, mThisIntegrationMethod);

                // Compute material reponse
                CalculateAllConstitutiveVariables(this_kinematic_variables, this_constitutive_variables, this_non_local_constitutive_variables, Values, point_number, integration_points, GetStressMeasure());

                if ( rOutput[point_number].size() != strain_size)
                    rOutput[point_number].resize( strain_size, false );

                rOutput[point_number] = this_constitutive_variables.StrainVector;
            }
    }else{
        BaseSolidElement::CalculateOnIntegrationPoints(rVariable, rOutput, rCurrentProcessInfo);
    }
}

/***********************************************************************************/
/***********************************************************************************/

const Parameters SmallDisplacementHexTwoNonLocalVariables::GetSpecifications() const
{
    const Parameters specifications = Parameters(R"({
        "time_integration"           : ["static","implicit","explicit"],
        "framework"                  : "lagrangian",
        "symmetric_lhs"              : true,
        "positive_definite_lhs"      : true,
        "output"                     : {
            "gauss_point"            : ["INTEGRATION_WEIGHT","STRAIN_ENERGY","ERROR_INTEGRATION_POINT","VON_MISES_STRESS","INSITU_STRESS","CAUCHY_STRESS_VECTOR","PK2_STRESS_VECTOR","GREEN_LAGRANGE_STRAIN_VECTOR","ALMANSI_STRAIN_VECTOR","CAUCHY_STRESS_TENSOR","PK2_STRESS_TENSOR","GREEN_LAGRANGE_STRAIN_TENSOR","ALMANSI_STRAIN_TENSOR","CONSTITUTIVE_MATRIX","DEFORMATION_GRADIENT","CONSTITUTIVE_LAW"],
            "nodal_historical"       : ["DISPLACEMENT","VELOCITY","ACCELERATION", "NON_LOCAL_VARIABLE"],
            "nodal_non_historical"   : [],
            "entity"                 : []
        },
        "required_variables"         : ["DISPLACEMENT", "NON_LOCAL_VARIABLE"],
        "required_dofs"              : [],
        "flags_used"                 : [],
        "compatible_geometries"      : ["Triangle2D3", "Triangle2D6", "Quadrilateral2D4", "Quadrilateral2D8", "Quadrilateral2D9","Tetrahedra3D4", "Prism3D6", "Prism3D15", "Hexahedra3D8", "Hexahedra3D20", "Hexahedra3D27", "Tetrahedra3D10"],
        "element_integrates_in_time" : true,
        "compatible_constitutive_laws": {
            "type"        : ["PlaneStrain","ThreeDimensional"],
            "dimension"   : ["2D","3D"],
            "strain_size" : [3,6]
        },
        "required_polynomial_degree_of_geometry" : -1,
        "documentation"   : "This is a pure displacement element"
    })");

    const SizeType dimension = GetGeometry().WorkingSpaceDimension();
    if (dimension == 2) {
        std::vector<std::string> dofs_2d({"DISPLACEMENT_X","DISPLACEMENT_Y", "NON_LOCAL_VARIABLE"});
        specifications["required_dofs"].SetStringArray(dofs_2d);
    } else {
        std::vector<std::string> dofs_3d({"DISPLACEMENT_X","DISPLACEMENT_Y","DISPLACEMENT_Z", "NON_LOCAL_VARIABLE"});
        specifications["required_dofs"].SetStringArray(dofs_3d);
    }

    return specifications;
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementHexTwoNonLocalVariables::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BaseSolidElement );
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementHexTwoNonLocalVariables::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BaseSolidElement );
}

/***********************************************************************************/
/***********************************************************************************/

} // Namespace Kratos


