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
#include "custom_elements/small_displacement_non_local_hex.h"
#include "custom_utilities/structural_mechanics_element_utilities.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos
{
SmallDisplacementNonLocalHex::SmallDisplacementNonLocalHex( IndexType NewId, GeometryType::Pointer pGeometry )
    : BaseSolidElement( NewId, pGeometry )
{
    //DO NOT ADD DOFS HERE!!!
    mConstitutiveMatrix.resize(7,7);
}

/***********************************************************************************/
/***********************************************************************************/

SmallDisplacementNonLocalHex::SmallDisplacementNonLocalHex( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
        : BaseSolidElement( NewId, pGeometry, pProperties )
{
    //DO NOT ADD DOFS HERE!!!
    mConstitutiveMatrix.resize(7,7);
}

/***********************************************************************************/
/***********************************************************************************/

Element::Pointer SmallDisplacementNonLocalHex::Create( IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties ) const
{
    return Kratos::make_intrusive<SmallDisplacementNonLocalHex>( NewId, GetGeometry().Create( ThisNodes ), pProperties );
}

/***********************************************************************************/
/***********************************************************************************/

Element::Pointer SmallDisplacementNonLocalHex::Create( IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties ) const
{
    return Kratos::make_intrusive<SmallDisplacementNonLocalHex>( NewId, pGeom, pProperties );
}

/***********************************************************************************/
/***********************************************************************************/

SmallDisplacementNonLocalHex::~SmallDisplacementNonLocalHex()
{
}

/***********************************************************************************/
/***********************************************************************************/

Element::Pointer SmallDisplacementNonLocalHex::Clone (
    IndexType NewId,
    NodesArrayType const& rThisNodes
    ) const
{
    KRATOS_TRY

    SmallDisplacementNonLocalHex::Pointer p_new_elem = Kratos::make_intrusive<SmallDisplacementNonLocalHex>(NewId, GetGeometry().Create(rThisNodes), pGetProperties());
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

void SmallDisplacementNonLocalHex::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    BaseType::Initialize(rCurrentProcessInfo);
    mpNonLocalVariable = &KratosComponents<Variable<double>>::Get(this->GetProperties()[NON_LOCAL_VARIABLE_NAME]);
}

/***********************************************************************************/
/***********************************************************************************/

bool SmallDisplacementNonLocalHex::UseElementProvidedStrain() const
{
    return true;
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementNonLocalHex::EquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo
    ) const
{
    KRATOS_TRY;

    auto& r_geometry = this->GetGeometry();
    const SizeType number_of_nodes = r_geometry.size();
    const SizeType dimension = r_geometry.WorkingSpaceDimension();
    const SizeType mat_size = number_of_nodes * (dimension+1);

    if (rResult.size() != mat_size){
        rResult.resize(mat_size,false);
    }

    const SizeType pos = r_geometry[0].GetDofPosition(DISPLACEMENT_X);

    if(dimension == 2) {
        for (IndexType i = 0; i < number_of_nodes; ++i) {
            const SizeType index = i * 3;
            rResult[index] = r_geometry[i].GetDof(DISPLACEMENT_X,pos).EquationId();
            rResult[index + 1] = r_geometry[i].GetDof(DISPLACEMENT_Y,pos+1).EquationId();
            rResult[index + 2] = r_geometry[i].GetDof(NON_LOCAL_VARIABLE).EquationId();
        }
    } else {
        for (IndexType i = 0; i < number_of_nodes; ++i) {
            const SizeType index = i * 4;
            rResult[index] = r_geometry[i].GetDof(DISPLACEMENT_X,pos).EquationId();
            rResult[index + 1] = r_geometry[i].GetDof(DISPLACEMENT_Y,pos+1).EquationId();
            rResult[index + 2] = r_geometry[i].GetDof(DISPLACEMENT_Z,pos+2).EquationId();
            rResult[index + 3] = r_geometry[i].GetDof(NON_LOCAL_VARIABLE).EquationId();
        }
    }

    KRATOS_CATCH("")
};

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementNonLocalHex::GetDofList(
    DofsVectorType& rElementalDofList,
    const ProcessInfo& rCurrentProcessInfo
    ) const
{
    KRATOS_TRY;

    auto& r_geometry = this->GetGeometry();
    const SizeType number_of_nodes = r_geometry.size();
    const SizeType dimension = r_geometry.WorkingSpaceDimension();
    const SizeType mat_size = number_of_nodes * (dimension+1);
    rElementalDofList.clear();
    rElementalDofList.reserve(mat_size);

    if(dimension == 2) {
        for (IndexType i = 0; i < number_of_nodes; ++i) {
            rElementalDofList.push_back(r_geometry[i].pGetDof(DISPLACEMENT_X));
            rElementalDofList.push_back( r_geometry[i].pGetDof(DISPLACEMENT_Y));
            rElementalDofList.push_back( r_geometry[i].pGetDof(NON_LOCAL_VARIABLE));
        }
    } else {
        for (IndexType i = 0; i < number_of_nodes; ++i) {
            rElementalDofList.push_back(r_geometry[i].pGetDof(DISPLACEMENT_X));
            rElementalDofList.push_back( r_geometry[i].pGetDof(DISPLACEMENT_Y));
            rElementalDofList.push_back( r_geometry[i].pGetDof(DISPLACEMENT_Z));
            rElementalDofList.push_back( r_geometry[i].pGetDof(NON_LOCAL_VARIABLE));
        }
    }

    KRATOS_CATCH("")
};

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementNonLocalHex::GetValuesVector(
    Vector& rValues,
    int Step
    ) const
{
    KRATOS_TRY
    auto& r_geometry = this->GetGeometry();
    const SizeType number_of_nodes = r_geometry.size();
    const SizeType dimension = r_geometry.WorkingSpaceDimension();
    const SizeType mat_size = number_of_nodes * (dimension+1);
    if (rValues.size() != mat_size) {
        rValues.resize(mat_size, false);
    }

    for (IndexType i = 0; i < number_of_nodes; ++i){
        const array_1d<double, 3 >& displacement = r_geometry[i].FastGetSolutionStepValue(DISPLACEMENT, Step);
        const SizeType index = i * (dimension +1);
        for(unsigned int k = 0; k < dimension; ++k)
        {
            rValues[index + k] = displacement[k];
        }
        rValues[index + dimension] = r_geometry[i].FastGetSolutionStepValue(NON_LOCAL_VARIABLE, Step);
    }
    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementNonLocalHex::GetFirstDerivativesVector(
    Vector& rValues,
    int Step
    ) const
{
    KRATOS_TRY
    auto& r_geometry = this->GetGeometry();
    const SizeType number_of_nodes = r_geometry.size();
    const SizeType dimension = r_geometry.WorkingSpaceDimension();
    const SizeType mat_size = number_of_nodes * (dimension+1);
    if (rValues.size() != mat_size){
        rValues.resize(mat_size, false);
    }
    for (IndexType i = 0; i < number_of_nodes; ++i) {
        const array_1d<double, 3 >& velocity = r_geometry[i].FastGetSolutionStepValue(VELOCITY, Step);
        const SizeType index = i * (dimension + 1);
        for(unsigned int k = 0; k < dimension; ++k){
            rValues[index + k] = velocity[k];
        }
        rValues[index + dimension] = 0.0;
    }
    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementNonLocalHex::GetSecondDerivativesVector(
    Vector& rValues,
    int Step
    ) const
{
    KRATOS_TRY
    auto& r_geometry = this->GetGeometry();
    const SizeType number_of_nodes = r_geometry.size();
    const SizeType dimension = r_geometry.WorkingSpaceDimension();
    const SizeType mat_size = number_of_nodes * (dimension+1);
    if (rValues.size() != mat_size){
        rValues.resize(mat_size, false);
    }
    for (IndexType i = 0; i < number_of_nodes; ++i) {
        const array_1d<double, 3 >& acceleration = r_geometry[i].FastGetSolutionStepValue(ACCELERATION, Step);
        const SizeType index = i * (dimension + 1);
        for(unsigned int k = 0; k < dimension; ++k){
            rValues[index + k] = acceleration[k];
        }
        rValues[index + dimension] = 0.0;
    }
    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementNonLocalHex::CalculateAll(
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
    const SizeType strain_size = (dimension == 3) ? 6 : 4;

    KinematicVariables this_kinematic_variables(strain_size, dimension, number_of_nodes);
    ConstitutiveVariables this_constitutive_variables(strain_size);
    NonLocalConstitutiveVariables this_non_local_constitutive_variables(strain_size);

    // Resizing as needed the LHS
    const SizeType mat_size = number_of_nodes * (dimension+1);

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
    Matrix Kud = ZeroMatrix(number_of_nodes * dimension , number_of_nodes);
    Matrix Kdu = ZeroMatrix(number_of_nodes, number_of_nodes * dimension);
    Matrix Kdd = ZeroMatrix(number_of_nodes , number_of_nodes);;
    Vector Fu = ZeroVector(number_of_nodes * dimension);
    Vector FNL = ZeroVector(number_of_nodes);

    // Computing in all integrations points
    for ( IndexType point_number = 0; point_number < integration_points.size(); ++point_number ) {
        // Contribution to external forces
        noalias(body_force) = this->GetBodyForce(integration_points, point_number);

        // Compute element kinematics B, F, DN_DX ...
        CalculateKinematicVariables(this_kinematic_variables, point_number, this->GetIntegrationMethod());

        // Compute material reponse
        CalculateAllConstitutiveVariables(this_kinematic_variables, this_constitutive_variables, this_non_local_constitutive_variables, Values, point_number, integration_points, GetStressMeasure());

        // Calculate and set non locl constitutive variables......................................................
        CalculateNonLocalConstitutiveVariables(this_kinematic_variables, this_non_local_constitutive_variables, Values, *mpNonLocalVariable, point_number);

        SetNonLocalConstitutiveVariables(mConstitutiveMatrix, this_constitutive_variables, this_non_local_constitutive_variables);

        // Calculating weights for integration on the reference configuration
        int_to_reference_weight = GetIntegrationWeight(integration_points, point_number, this_kinematic_variables.detJ0);

        if ( dimension == 2 && GetProperties().Has( THICKNESS ))
            int_to_reference_weight *= GetProperties()[THICKNESS];

        if ( CalculateStiffnessMatrixFlag ) { // Calculation of the matrix is required
            // Contributions to stiffness matrix calculated on the reference config

            CalculateAndAddKuu(Kuu, this_kinematic_variables.B, this_constitutive_variables.D, int_to_reference_weight );

            CalculateAndAddKmud(Kud, this_kinematic_variables.B, this_non_local_constitutive_variables.DuNL, this_kinematic_variables.N, int_to_reference_weight);

            CalculateAndAddKmdu(Kdu, this_kinematic_variables.B, this_non_local_constitutive_variables.DNLu, this_kinematic_variables.N, int_to_reference_weight);

            CalculateAndAddKmdd(Kdd, this_kinematic_variables.DN_DX, this_non_local_constitutive_variables.DNLNL, this_kinematic_variables.N, int_to_reference_weight);

        }

        if ( CalculateResidualVectorFlag ) { // Calculation of the matrix is required

            CalculateAndAddResidualForceVector(Fu, this_kinematic_variables, rCurrentProcessInfo, body_force, this_constitutive_variables.StressVector, int_to_reference_weight);

            CalculateAndAddResidualNonLocalVector(FNL, this_kinematic_variables, this_non_local_constitutive_variables, rCurrentProcessInfo, int_to_reference_weight);

        }
    }
    //assemble right hand side and left hand sides
    AssembleRHSAndLHS(rLeftHandSideMatrix, rRightHandSideVector, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag, Kuu, Kud, Kdu, Kdd, Fu, FNL);
    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementNonLocalHex::FinalizeSolutionStep( const ProcessInfo& rCurrentProcessInfo )
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

void SmallDisplacementNonLocalHex::CalculateKinematicVariables(
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

void SmallDisplacementNonLocalHex::CalculateAllConstitutiveVariables(
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

void SmallDisplacementNonLocalHex::SetConstitutiveVariables(
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

void SmallDisplacementNonLocalHex::CalculateNonLocalConstitutiveVariables(
    KinematicVariables& rThisKinematicVariables,
    NonLocalConstitutiveVariables& rThisNonLocalConstitutiveVariables,
    ConstitutiveLaw::Parameters& rValues,
    const Variable<double>& rThisVariable,
    const IndexType PointNumber
    )
{
    KRATOS_TRY;
    const auto& r_geometry = GetGeometry();
    const auto& N = rThisKinematicVariables.N;
    double local_variable_gp = 0.0;
    rThisNonLocalConstitutiveVariables.NonLocal_Variable_GP = 0;

    for (unsigned int i = 0; i < N.size(); ++i) {
        rThisNonLocalConstitutiveVariables.NonLocal_Variable_GP += N[i] * r_geometry[i].FastGetSolutionStepValue(NON_LOCAL_VARIABLE, 0);
    }
    mConstitutiveLawVector[PointNumber]->CalculateValue(rValues, rThisVariable, local_variable_gp);
    rThisNonLocalConstitutiveVariables.Local_Variable_GP = local_variable_gp;
    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementNonLocalHex::SetNonLocalConstitutiveVariables(
    const Matrix& ConstitutiveMatrix,
    ConstitutiveVariables& rThisConstitutiveVariables,
    NonLocalConstitutiveVariables& rThisNonLocalConstitutiveVariables
    )const
{
    KRATOS_TRY;
    const std::size_t num_rows_D = rThisConstitutiveVariables.D.size1();
    const std::size_t num_cols_D = rThisConstitutiveVariables.D.size2();
    const std::size_t size_DuNL  = rThisNonLocalConstitutiveVariables.DuNL.size();
    const std::size_t size_DNLu  = rThisNonLocalConstitutiveVariables.DNLu.size();

    for (std::size_t i = 0; i < num_rows_D ; ++i) {
        for (std::size_t j = 0; j < num_cols_D ; ++j) {
            rThisConstitutiveVariables.D(i,j) = ConstitutiveMatrix(i, j) ;
        }
    }
    for (std::size_t i = 0; i < size_DuNL; ++i) {
        rThisNonLocalConstitutiveVariables.DuNL[i] = ConstitutiveMatrix(i, num_cols_D) ;
    }
    for (std::size_t i = 0; i < size_DNLu; ++i) {
        rThisNonLocalConstitutiveVariables.DNLu[i] =  ConstitutiveMatrix(num_rows_D, i);
    }
    rThisNonLocalConstitutiveVariables.DNLNL =  ConstitutiveMatrix(num_rows_D, num_cols_D);

    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementNonLocalHex::CalculateB(
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

void SmallDisplacementNonLocalHex::ComputeEquivalentF(
    Matrix& rF,
    const Vector& rStrainTensor
    ) const
{
    StructuralMechanicsElementUtilities::ComputeEquivalentF(*this, rStrainTensor, rF);
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementNonLocalHex::CalculateAndAddKmud(
    Matrix& rStiffnessMatrixKud,
    const Matrix& B,
    const Vector& DuNL,
    const Vector& N,
    const double IntegrationWeight
    ) const
{
    KRATOS_TRY

    rStiffnessMatrixKud += IntegrationWeight * prod( trans( B ), outer_prod(DuNL, N));

    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementNonLocalHex::CalculateAndAddKmdu(
    Matrix& rStiffnessMatrixKdu,
    const Matrix& B,
    const Vector& DNLu,
    const Vector& N,
    const double IntegrationWeight
    ) const
{
    KRATOS_TRY

    rStiffnessMatrixKdu  -= IntegrationWeight * prod( outer_prod(N, DNLu), B);

    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementNonLocalHex::CalculateAndAddKmdd(
    Matrix& rStiffnessMatrixKdd,
    const Matrix& DN_DX,
    const double& DNLNL,
    const Vector& N,
    const double IntegrationWeight
    ) const
{
    KRATOS_TRY
    const Properties& r_properties = GetProperties();
    const double Lc = r_properties[CHARACTERISTIC_INTERNAL_LENGTH];
    rStiffnessMatrixKdd += IntegrationWeight *( (1-DNLNL) * outer_prod(N, N) + pow(Lc, 2) * prod( DN_DX , trans(DN_DX)));
    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementNonLocalHex::CalculateAndAddResidualForceVector(
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

void SmallDisplacementNonLocalHex::CalculateAndAddResidualNonLocalVector(
    Vector& rResidualNonLocalVector,
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
    const double rNonLocal_Variable_GP = rThisNonLocalConstitutiveVariables.NonLocal_Variable_GP;
    const double rLocal_Variable_GP = rThisNonLocalConstitutiveVariables.Local_Variable_GP;
    Vector NL_variable_vector = ZeroVector(rN.size());

    for (unsigned int i = 0; i < rN.size(); ++i) {
        NL_variable_vector[i] = r_geometry[i].FastGetSolutionStepValue(NON_LOCAL_VARIABLE, 0);
    }
    rResidualNonLocalVector -= IntegrationWeight * rN * (rNonLocal_Variable_GP-rLocal_Variable_GP);
    rResidualNonLocalVector -= IntegrationWeight * pow(Lc, 2) * prod(prod(rdNdX, trans(rdNdX)),NL_variable_vector );
    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementNonLocalHex::AssembleRHSAndLHS(
    MatrixType& LeftHandSideMatrix,
    VectorType& RightHandSideVector,
    const bool StiffnessMatrixFlag,
    const bool ResidualVectorFlag,
    const Matrix& Kuu,
    const Matrix& Kud,
    const Matrix& Kdu,
    const Matrix& Kdd,
    const Vector& rFu,
    const Vector& rFNL
    )
{
    KRATOS_TRY;
    auto& r_geometry = this->GetGeometry();
    const SizeType nnod = r_geometry.size();
    const SizeType dimension = r_geometry.WorkingSpaceDimension();
    Vector Pos_u = ZeroVector(nnod * dimension);
    Vector Pos_NL = ZeroVector(nnod);

    for(unsigned int a = 0; a < nnod ; a++){
        const SizeType index = a * dimension;
        for(unsigned int b = 0; b < dimension; b++){
            Pos_u(index+b) = index + a + b;
        }
        Pos_NL(a) = index + dimension + a ;
    }
    for(unsigned int i = 0; i < (nnod * dimension); i++){
        if ( ResidualVectorFlag ) {
            RightHandSideVector(Pos_u[i]) = rFu[i];
            if(i < nnod){
                RightHandSideVector(Pos_NL[i]) = rFNL[i];
            }
        }
        if ( StiffnessMatrixFlag ) {
            for(unsigned int j = 0; j < (nnod * dimension); j++){
                LeftHandSideMatrix(Pos_u[i],Pos_u[j]) = Kuu(i,j);
                if(j < nnod){
                    LeftHandSideMatrix(Pos_u[i],Pos_NL[j]) = Kud(i,j);
                }
                if(i < nnod){
                    LeftHandSideMatrix(Pos_NL[i],Pos_u[j]) = Kdu(i,j);
                }
                if(i < nnod && j < nnod){
                    LeftHandSideMatrix(Pos_NL[i],Pos_NL[j]) = Kdd(i,j);
                }
            }
        }
    }
    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/
void SmallDisplacementNonLocalHex::CalculateAndAddKuu(
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

void SmallDisplacementNonLocalHex::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BaseSolidElement );
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementNonLocalHex::CalculateOnIntegrationPoints(
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

void SmallDisplacementNonLocalHex::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BaseSolidElement );
}

const Parameters SmallDisplacementNonLocalHex::GetSpecifications() const
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

} // Namespace Kratos


