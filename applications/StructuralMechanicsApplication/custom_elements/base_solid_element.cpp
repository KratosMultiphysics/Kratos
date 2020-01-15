// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Riccardo Rossi
//                   Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h"

// Application includes
#include "custom_elements/base_solid_element.h"
#include "structural_mechanics_application_variables.h"
#include "custom_utilities/structural_mechanics_element_utilities.h"

namespace Kratos
{
void BaseSolidElement::Initialize()
{
    KRATOS_TRY

    if( GetProperties().Has(INTEGRATION_ORDER) ) {
        const SizeType integration_order = GetProperties()[INTEGRATION_ORDER];
        switch ( integration_order )
        {
        case 1:
            mThisIntegrationMethod = GeometryData::GI_GAUSS_1;
            break;
        case 2:
            mThisIntegrationMethod = GeometryData::GI_GAUSS_2;
            break;
        case 3:
            mThisIntegrationMethod = GeometryData::GI_GAUSS_3;
            break;
        case 4:
            mThisIntegrationMethod = GeometryData::GI_GAUSS_4;
            break;
        case 5:
            mThisIntegrationMethod = GeometryData::GI_GAUSS_5;
            break;
        default:
            KRATOS_WARNING("BaseSolidElement") << "Integration order " << integration_order << " is not available, using default integration order for the geometry" << std::endl;
            mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod();
        }
    } else {
        mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod();
    }

    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints(this->GetIntegrationMethod());

    //Constitutive Law initialisation
    if ( mConstitutiveLawVector.size() != integration_points.size() )
        mConstitutiveLawVector.resize( integration_points.size() );

    InitializeMaterial();

    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::InitializeSolutionStep( ProcessInfo& rCurrentProcessInfo )
{
    // We initialize the material reponse if required
    bool required = false;
    for ( IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number ) {
        if (mConstitutiveLawVector[point_number]->RequiresInitializeMaterialResponse()) {
            required = true;
            break;
        }
    }
    if (required) {
        const SizeType number_of_nodes = GetGeometry().size();
        const SizeType dimension = GetGeometry().WorkingSpaceDimension();
        const SizeType strain_size = mConstitutiveLawVector[0]->GetStrainSize();

        KinematicVariables this_kinematic_variables(strain_size, dimension, number_of_nodes);
        ConstitutiveVariables this_constitutive_variables(strain_size);

        // Create constitutive law parameters:
        ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

        // Set constitutive law flags:
        Flags& ConstitutiveLawOptions=Values.GetOptions();
        ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, UseElementProvidedStrain());
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);

        Values.SetStrainVector(this_constitutive_variables.StrainVector);
        Values.SetStressVector(this_constitutive_variables.StressVector);
        Values.SetConstitutiveMatrix(this_constitutive_variables.D);

        // Reading integration points
        const GeometryType& r_geometry = GetGeometry();
        const Properties& r_properties = GetProperties();
        const auto& N_values = r_geometry.ShapeFunctionsValues(mThisIntegrationMethod);

        // Reading integration points
        const GeometryType::IntegrationPointsArrayType& integration_points = r_geometry.IntegrationPoints(mThisIntegrationMethod);

        for ( IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number ) {
            if (mConstitutiveLawVector[point_number]->RequiresInitializeMaterialResponse()) {
                // Compute element kinematics B, F, DN_DX ...
                CalculateKinematicVariables(this_kinematic_variables, point_number, mThisIntegrationMethod);

                // Compute constitutive law variables
                SetConstitutiveVariables(this_kinematic_variables, this_constitutive_variables, Values, point_number, integration_points);

                // Call the constitutive law to update material variables
                mConstitutiveLawVector[point_number]->InitializeMaterialResponse(Values, GetStressMeasure());

                // TODO: Deprecated, remove this
                mConstitutiveLawVector[point_number]->InitializeSolutionStep( r_properties, r_geometry, row( N_values, point_number ), rCurrentProcessInfo);
            }
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::InitializeNonLinearIteration( ProcessInfo& rCurrentProcessInfo )
{
    const GeometryType& r_geometry = GetGeometry();
    const Properties& r_properties = GetProperties();
    const auto& N_values = r_geometry.ShapeFunctionsValues(mThisIntegrationMethod);
    for ( IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number ) {
        // TODO: Deprecated, remove this
        mConstitutiveLawVector[point_number]->InitializeNonLinearIteration( r_properties, r_geometry, row( N_values, point_number ), rCurrentProcessInfo);
    }
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::FinalizeNonLinearIteration( ProcessInfo& rCurrentProcessInfo )
{
    const GeometryType& r_geometry = GetGeometry();
    const Properties& r_properties = GetProperties();
    const auto& N_values = r_geometry.ShapeFunctionsValues(mThisIntegrationMethod);
    for ( IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number ) {
        // TODO: Deprecated, remove this
        mConstitutiveLawVector[point_number]->FinalizeNonLinearIteration( r_properties, r_geometry, row( N_values, point_number ), rCurrentProcessInfo);
    }
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::FinalizeSolutionStep( ProcessInfo& rCurrentProcessInfo )
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
        const SizeType number_of_nodes = GetGeometry().size();
        const SizeType dimension = GetGeometry().WorkingSpaceDimension();
        const SizeType strain_size = mConstitutiveLawVector[0]->GetStrainSize();

        KinematicVariables this_kinematic_variables(strain_size, dimension, number_of_nodes);
        ConstitutiveVariables this_constitutive_variables(strain_size);

        // Create constitutive law parameters:
        ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

        // Set constitutive law flags:
        Flags& ConstitutiveLawOptions=Values.GetOptions();
        ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, UseElementProvidedStrain());
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);

        Values.SetStrainVector(this_constitutive_variables.StrainVector);
        Values.SetStressVector(this_constitutive_variables.StressVector);
        Values.SetConstitutiveMatrix(this_constitutive_variables.D);

        // Reading integration points
        const GeometryType& r_geometry = GetGeometry();
        const Properties& r_properties = GetProperties();
        const auto& N_values = r_geometry.ShapeFunctionsValues(mThisIntegrationMethod);

        // Reading integration points
        const GeometryType::IntegrationPointsArrayType& integration_points = r_geometry.IntegrationPoints(mThisIntegrationMethod);

        for ( IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number ) {
            if (mConstitutiveLawVector[point_number]->RequiresFinalizeMaterialResponse()) {
                // Compute element kinematics B, F, DN_DX ...
                CalculateKinematicVariables(this_kinematic_variables, point_number, mThisIntegrationMethod);

                // Compute constitutive law variables
                SetConstitutiveVariables(this_kinematic_variables, this_constitutive_variables, Values, point_number, integration_points);

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

void BaseSolidElement::InitializeMaterial()
{
    KRATOS_TRY

    if ( GetProperties()[CONSTITUTIVE_LAW] != nullptr ) {
        const GeometryType& r_geometry = GetGeometry();
        const Properties& r_properties = GetProperties();
        const auto& N_values = r_geometry.ShapeFunctionsValues(mThisIntegrationMethod);
        for ( IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number ) {
            mConstitutiveLawVector[point_number] = GetProperties()[CONSTITUTIVE_LAW]->Clone();
            mConstitutiveLawVector[point_number]->InitializeMaterial( r_properties, r_geometry, row(N_values , point_number ));
        }
    } else
        KRATOS_ERROR << "A constitutive law needs to be specified for the element with ID " << this->Id() << std::endl;

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

ConstitutiveLaw::StressMeasure BaseSolidElement::GetStressMeasure() const
{
    return ConstitutiveLaw::StressMeasure_PK2;
}

/***********************************************************************************/
/***********************************************************************************/

bool BaseSolidElement::UseElementProvidedStrain() const
{
    return false;
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::ResetConstitutiveLaw()
{
    KRATOS_TRY

    if ( GetProperties()[CONSTITUTIVE_LAW] != nullptr ) {
        const GeometryType& r_geometry = GetGeometry();
        const Properties& r_properties = GetProperties();
        const auto& N_values = r_geometry.ShapeFunctionsValues(mThisIntegrationMethod);
        for ( IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number )
            mConstitutiveLawVector[point_number]->ResetMaterial( r_properties,  r_geometry, row( N_values, point_number ) );
    }

    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

Element::Pointer BaseSolidElement::Clone (
    IndexType NewId,
    NodesArrayType const& rThisNodes
    ) const
{
    KRATOS_TRY

    KRATOS_WARNING("BaseSolidElement") << " Call BaseSolidElement (base class) Clone " << std::endl;

    BaseSolidElement::Pointer p_new_elem = Kratos::make_intrusive<BaseSolidElement>(NewId, GetGeometry().Create(rThisNodes), pGetProperties());
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

void BaseSolidElement::EquationIdVector(
    EquationIdVectorType& rResult,
    ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    const SizeType number_of_nodes = GetGeometry().size();
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();

    if (rResult.size() != dimension * number_of_nodes)
        rResult.resize(dimension * number_of_nodes,false);

    const SizeType pos = this->GetGeometry()[0].GetDofPosition(DISPLACEMENT_X);

    if(dimension == 2) {
        for (IndexType i = 0; i < number_of_nodes; ++i) {
            const SizeType index = i * 2;
            rResult[index] = GetGeometry()[i].GetDof(DISPLACEMENT_X,pos).EquationId();
            rResult[index + 1] = GetGeometry()[i].GetDof(DISPLACEMENT_Y,pos+1).EquationId();
        }
    } else {
        for (IndexType i = 0; i < number_of_nodes; ++i) {
            const SizeType index = i * 3;
            rResult[index] = GetGeometry()[i].GetDof(DISPLACEMENT_X,pos).EquationId();
            rResult[index + 1] = GetGeometry()[i].GetDof(DISPLACEMENT_Y,pos+1).EquationId();
            rResult[index + 2] = GetGeometry()[i].GetDof(DISPLACEMENT_Z,pos+2).EquationId();
        }
    }

    KRATOS_CATCH("")
};

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::GetDofList(
    DofsVectorType& rElementalDofList,
    ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    const SizeType number_of_nodes = GetGeometry().size();
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();
    rElementalDofList.resize(0);
    rElementalDofList.reserve(dimension*number_of_nodes);

    if(dimension == 2) {
        for (IndexType i = 0; i < number_of_nodes; ++i) {
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_X));
            rElementalDofList.push_back( GetGeometry()[i].pGetDof(DISPLACEMENT_Y));
        }
    } else {
        for (IndexType i = 0; i < number_of_nodes; ++i) {
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_X));
            rElementalDofList.push_back( GetGeometry()[i].pGetDof(DISPLACEMENT_Y));
            rElementalDofList.push_back( GetGeometry()[i].pGetDof(DISPLACEMENT_Z));
        }
    }

    KRATOS_CATCH("")
};

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::GetValuesVector(
    Vector& rValues,
    int Step
    )
{
    const SizeType number_of_nodes = GetGeometry().size();
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();
    const SizeType mat_size = number_of_nodes * dimension;
    if (rValues.size() != mat_size)
        rValues.resize(mat_size, false);
    for (IndexType i = 0; i < number_of_nodes; ++i)
    {
        const array_1d<double, 3 >& displacement = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT, Step);
        const SizeType index = i * dimension;
        for(unsigned int k = 0; k < dimension; ++k)
        {
            rValues[index + k] = displacement[k];
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::GetFirstDerivativesVector(
    Vector& rValues,
    int Step
    )
{
    const SizeType number_of_nodes = GetGeometry().size();
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();
    const SizeType mat_size = number_of_nodes * dimension;
    if (rValues.size() != mat_size)
        rValues.resize(mat_size, false);
    for (IndexType i = 0; i < number_of_nodes; ++i) {
        const array_1d<double, 3 >& velocity = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, Step);
        const SizeType index = i * dimension;
        for(unsigned int k = 0; k < dimension; ++k)
            rValues[index + k] = velocity[k];
    }
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::GetSecondDerivativesVector(
    Vector& rValues,
    int Step
    )
{
    const SizeType number_of_nodes = GetGeometry().size();
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();
    const SizeType mat_size = number_of_nodes * dimension;
    if (rValues.size() != mat_size)
        rValues.resize(mat_size, false);
    for (IndexType i = 0; i < number_of_nodes; ++i) {
        const array_1d<double, 3 >& acceleration = GetGeometry()[i].FastGetSolutionStepValue(ACCELERATION, Step);
        const SizeType index = i * dimension;
        for(unsigned int k = 0; k < dimension; ++k)
            rValues[index + k] = acceleration[k];
    }
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::AddExplicitContribution(
    const VectorType& rRHSVector,
    const Variable<VectorType>& rRHSVariable,
    Variable<double>& rDestinationVariable,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    auto& r_geom = this->GetGeometry();
    const SizeType dimension = r_geom.WorkingSpaceDimension();
    const SizeType number_of_nodes = r_geom.size();
    const SizeType mat_size = number_of_nodes * dimension;

    // Compiting the nodal mass
    if (rDestinationVariable == NODAL_MASS ) {
        VectorType element_mass_vector(mat_size);
        this->CalculateLumpedMassVector(element_mass_vector);

        for (IndexType i = 0; i < number_of_nodes; ++i) {
            const IndexType index = i * dimension;

            #pragma omp atomic
            r_geom[i].GetValue(NODAL_MASS) += element_mass_vector[index];
        }
    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::AddExplicitContribution(
    const VectorType& rRHSVector,
    const Variable<VectorType>& rRHSVariable,
    Variable<array_1d<double, 3>>& rDestinationVariable,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    auto& r_geom = this->GetGeometry();
    const auto& r_prop = this->GetProperties();
    const SizeType dimension = r_geom.WorkingSpaceDimension();
    const SizeType number_of_nodes = r_geom.size();
    const SizeType element_size = dimension * number_of_nodes;

    Vector damping_residual_contribution = ZeroVector(element_size);

    // Calculate damping contribution to residual -->
    if (r_prop.Has(RAYLEIGH_ALPHA) || r_prop.Has(RAYLEIGH_BETA)) {
        Vector current_nodal_velocities = ZeroVector(element_size);
        this->GetFirstDerivativesVector(current_nodal_velocities);

        Matrix damping_matrix(element_size, element_size);
        this->CalculateDampingMatrixWithLumpedMass(damping_matrix, rCurrentProcessInfo);

        // Current residual contribution due to damping
        noalias(damping_residual_contribution) = prod(damping_matrix, current_nodal_velocities);
    }

    // Computing the force residual
    if (rRHSVariable == RESIDUAL_VECTOR && rDestinationVariable == FORCE_RESIDUAL) {
        for (IndexType i = 0; i < number_of_nodes; ++i) {
            const IndexType index = dimension * i;

            array_1d<double, 3>& r_force_residual = r_geom[i].FastGetSolutionStepValue(FORCE_RESIDUAL);

            for (IndexType j = 0; j < dimension; ++j) {
                #pragma omp atomic
                r_force_residual[j] += rRHSVector[index + j] - damping_residual_contribution[index + j];
            }
        }
    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::CalculateLocalSystem(
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

void BaseSolidElement::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
                                             ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;
    VectorType RHS;
    CalculateLocalSystem(rLeftHandSideMatrix, RHS, rCurrentProcessInfo);
    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::CalculateRightHandSide(
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

void BaseSolidElement::CalculateMassMatrix(
    MatrixType& rMassMatrix,
    ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    const auto& r_prop = GetProperties();

    const auto& r_geom = GetGeometry();
    SizeType dimension = r_geom.WorkingSpaceDimension();
    SizeType number_of_nodes = r_geom.size();
    SizeType mat_size = dimension * number_of_nodes;

    // Clear matrix
    if (rMassMatrix.size1() != mat_size || rMassMatrix.size2() != mat_size)
        rMassMatrix.resize( mat_size, mat_size, false );
    rMassMatrix = ZeroMatrix( mat_size, mat_size );

    // Checking density
    KRATOS_ERROR_IF_NOT(r_prop.Has(DENSITY)) << "DENSITY has to be provided for the calculation of the MassMatrix!" << std::endl;

    // Checking if computing lumped mass matrix
    const bool compute_lumped_mass_matrix =  r_prop.Has(COMPUTE_LUMPED_MASS_MATRIX) ? r_prop[COMPUTE_LUMPED_MASS_MATRIX] : false;

    // LUMPED MASS MATRIX
    if (compute_lumped_mass_matrix) {
        VectorType temp_vector(mat_size);
        CalculateLumpedMassVector(temp_vector);
        for (IndexType i = 0; i < mat_size; ++i)
            rMassMatrix(i, i) = temp_vector[i];
    } else { // CONSISTENT MASS
        const double density = StructuralMechanicsElementUtilities::GetDensityForMassMatrixComputation(*this);
        const double thickness = (dimension == 2 && r_prop.Has(THICKNESS)) ? r_prop[THICKNESS] : 1.0;

        Matrix J0(dimension, dimension);

        IntegrationMethod integration_method = IntegrationUtilities::GetIntegrationMethodForExactMassMatrixEvaluation(r_geom);
        const GeometryType::IntegrationPointsArrayType& integration_points = r_geom.IntegrationPoints( integration_method );
        const Matrix& Ncontainer = r_geom.ShapeFunctionsValues(integration_method);

        for ( IndexType point_number = 0; point_number < integration_points.size(); ++point_number ) {
            GeometryUtils::JacobianOnInitialConfiguration(
                r_geom, integration_points[point_number], J0);
            const double detJ0 = MathUtils<double>::DetMat(J0);
            const double integration_weight =
                GetIntegrationWeight(integration_points, point_number, detJ0) * thickness;
            const Vector& rN = row(Ncontainer,point_number);

            for ( IndexType i = 0; i < number_of_nodes; ++i ) {
                const SizeType index_i = i * dimension;

                for ( IndexType j = 0; j < number_of_nodes; ++j ) {
                    const SizeType index_j = j * dimension;
                    const double NiNj_weight = rN[i] * rN[j] * integration_weight * density;

                    for ( IndexType k = 0; k < dimension; ++k )
                        rMassMatrix( index_i + k, index_j + k ) += NiNj_weight;
                }
            }
        }
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::CalculateDampingMatrix(
    MatrixType& rDampingMatrix,
    ProcessInfo& rCurrentProcessInfo
    )
{
    const unsigned int mat_size = GetGeometry().PointsNumber() * GetGeometry().WorkingSpaceDimension();

    StructuralMechanicsElementUtilities::CalculateRayleighDampingMatrix(
        *this,
        rDampingMatrix,
        rCurrentProcessInfo,
        mat_size);
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::CalculateOnIntegrationPoints(
    const Variable<bool>& rVariable,
    std::vector<bool>& rOutput,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    const GeometryType::IntegrationPointsArrayType &integration_points = GetGeometry().IntegrationPoints(this->GetIntegrationMethod());

    const SizeType number_of_integration_points = integration_points.size();
    if (rOutput.size() != number_of_integration_points)
        rOutput.resize(number_of_integration_points);

    if (mConstitutiveLawVector[0]->Has( rVariable)) {
        for ( IndexType point_number = 0; point_number <number_of_integration_points; ++point_number ) {
            bool value;
            mConstitutiveLawVector[point_number]->GetValue( rVariable, value);
            rOutput[point_number] = value;
        }
    } else {
        // Create constitutive law parameters:
        ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

        for ( IndexType ii = 0; ii < mConstitutiveLawVector.size(); ++ii ) {
            bool solution;
            solution = mConstitutiveLawVector[ii]->CalculateValue( Values, rVariable, solution);
            rOutput[ii] = solution;
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::CalculateOnIntegrationPoints(
    const Variable<int>& rVariable,
    std::vector<int>& rOutput,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    const GeometryType::IntegrationPointsArrayType &integration_points = GetGeometry().IntegrationPoints(this->GetIntegrationMethod());

    const SizeType number_of_integration_points = integration_points.size();
    if (rOutput.size() != number_of_integration_points)
        rOutput.resize(number_of_integration_points);

    if (mConstitutiveLawVector[0]->Has( rVariable)) {
        GetValueOnConstitutiveLaw(rVariable, rOutput);
    } else {
        CalculateOnConstitutiveLaw(rVariable, rOutput, rCurrentProcessInfo);
    }
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::CalculateOnIntegrationPoints(
    const Variable<double>& rVariable,
    std::vector<double>& rOutput,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    const GeometryType::IntegrationPointsArrayType &integration_points = GetGeometry().IntegrationPoints(this->GetIntegrationMethod());

    const std::size_t number_of_integration_points = integration_points.size();
    const auto& r_geometry = GetGeometry();

    if ( rOutput.size() != number_of_integration_points )
        rOutput.resize( number_of_integration_points );

    if (mConstitutiveLawVector[0]->Has( rVariable)) {
        GetValueOnConstitutiveLaw(rVariable, rOutput);
    } else {
        if (rVariable == INTEGRATION_WEIGHT) {
            const SizeType number_of_nodes = r_geometry.size();
            const SizeType dimension = r_geometry.WorkingSpaceDimension();
            const SizeType strain_size = mConstitutiveLawVector[0]->GetStrainSize();

            KinematicVariables this_kinematic_variables(strain_size, dimension, number_of_nodes);

            for (IndexType point_number = 0; point_number < number_of_integration_points; ++point_number) {
                this_kinematic_variables.detJ0 = CalculateDerivativesOnReferenceConfiguration(this_kinematic_variables.J0,
                                                                                    this_kinematic_variables.InvJ0,
                                                                                    this_kinematic_variables.DN_DX,
                                                                                    point_number,
                                                                                    this->GetIntegrationMethod());

                double integration_weight = GetIntegrationWeight(integration_points,
                                                                    point_number,
                                                                    this_kinematic_variables.detJ0);

                if (dimension == 2 && this->GetProperties().Has(THICKNESS))
                    integration_weight *= this->GetProperties()[THICKNESS];

                rOutput[point_number] = integration_weight;
            }
        } else if ( rVariable == STRAIN_ENERGY ) {
            const SizeType number_of_nodes = r_geometry.size();
            const SizeType dimension = r_geometry.WorkingSpaceDimension();
            const SizeType strain_size = mConstitutiveLawVector[0]->GetStrainSize();

            KinematicVariables this_kinematic_variables(strain_size, dimension, number_of_nodes);
            ConstitutiveVariables this_constitutive_variables(strain_size);

            // Create constitutive law parameters:
            ConstitutiveLaw::Parameters Values(r_geometry,GetProperties(),rCurrentProcessInfo);

            // Set constitutive law flags:
            Flags& ConstitutiveLawOptions=Values.GetOptions();
            ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, UseElementProvidedStrain());
            ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, false);
            ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);

            // Reading integration points
            const GeometryType::IntegrationPointsArrayType& integration_points = r_geometry.IntegrationPoints(this->GetIntegrationMethod());

            // If strain has to be computed inside of the constitutive law with PK2
            Values.SetStrainVector(this_constitutive_variables.StrainVector); //this is the input  parameter

            for (IndexType point_number = 0; point_number < number_of_integration_points; ++point_number) {
                // Compute element kinematics B, F, DN_DX ...
                CalculateKinematicVariables(this_kinematic_variables, point_number, this->GetIntegrationMethod());

                // Compute constitutive law variables
                SetConstitutiveVariables(this_kinematic_variables, this_constitutive_variables, Values, point_number, integration_points);

                double StrainEnergy = 0.0;

                mConstitutiveLawVector[point_number]->CalculateValue(Values, STRAIN_ENERGY, StrainEnergy);

                rOutput[point_number] = StrainEnergy;
            }
        } else if ( rVariable == ERROR_INTEGRATION_POINT ) {
            const SizeType number_of_nodes = r_geometry.size();
            const SizeType dimension = r_geometry.WorkingSpaceDimension();
            const SizeType strain_size = mConstitutiveLawVector[0]->GetStrainSize();

            KinematicVariables this_kinematic_variables(strain_size, dimension, number_of_nodes);
            ConstitutiveVariables this_constitutive_variables(strain_size);

            // Create constitutive law parameters:
            ConstitutiveLaw::Parameters Values(r_geometry,GetProperties(),rCurrentProcessInfo);

            // Set constitutive law flags:
            Flags& ConstitutiveLawOptions=Values.GetOptions();
            ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, UseElementProvidedStrain());
            ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, false);
            ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

            // Reading integration points
            const GeometryType::IntegrationPointsArrayType& integration_points = r_geometry.IntegrationPoints(  );

            //Calculate Cauchy Stresses from the FE solution
            std::vector<Vector> sigma_FE_solution(number_of_nodes);
            const Variable<Vector>& r_variable_stress = CAUCHY_STRESS_VECTOR;
            CalculateOnIntegrationPoints(r_variable_stress, sigma_FE_solution, rCurrentProcessInfo);

            // calculate the determinatn of the Jacobian in the current configuration
            Vector detJ(number_of_integration_points);
            detJ = r_geometry.DeterminantOfJacobian(detJ);

            // If strain has to be computed inside of the constitutive law with PK2
            Values.SetStrainVector(this_constitutive_variables.StrainVector); //this is the input  parameter

            if (r_geometry[0].Has(RECOVERED_STRESS)) {
                for (IndexType point_number = 0; point_number < number_of_integration_points; point_number++) {
                    // Compute element kinematics B, F, DN_DX ...
                    CalculateKinematicVariables(this_kinematic_variables, point_number, this->GetIntegrationMethod());

                    // Compute material reponse
                    CalculateConstitutiveVariables(this_kinematic_variables, this_constitutive_variables, Values, point_number, integration_points, GetStressMeasure());

                    double integration_weight = GetIntegrationWeight(integration_points, point_number, detJ[point_number]);

                    if (dimension == 2 && this->GetProperties().Has(THICKNESS))
                        integration_weight *= this->GetProperties()[THICKNESS];

                    // Calculate recovered stresses at integration points
                    Vector sigma_recovered = ZeroVector(strain_size);

                    // sigma_recovered = sum(N_i * sigma_recovered_i)
                    for (IndexType node_number=0; node_number<number_of_nodes; node_number++) {
                        const auto& r_sigma_recovered_node = r_geometry[node_number].GetValue(RECOVERED_STRESS);
                        for (IndexType stress_component = 0; stress_component<strain_size; stress_component++) {
                            sigma_recovered[stress_component] += this_kinematic_variables.N[node_number] * r_sigma_recovered_node[stress_component];
                        }
                    }

                    // Calculate error_sigma
                    Vector error_sigma(strain_size);
                    error_sigma = sigma_recovered - sigma_FE_solution[point_number];

                    // For debug
                    KRATOS_TRACE("ERROR_INTEGRATION_POINT")
                    <<"sigma recovered: " << sigma_recovered << std::endl
                    <<"sigma FE: " << sigma_FE_solution[point_number] << std::endl;

                    // Calculate inverse of material matrix
                    Matrix invD(strain_size,strain_size);
                    double detD;
                    MathUtils<double>::InvertMatrix(this_constitutive_variables.D, invD,detD);

                    // Calculate error_energy
                    rOutput[point_number] = integration_weight * inner_prod(error_sigma, prod(invD, error_sigma));
                }
            } else {
                for (IndexType point_number = 0; point_number < number_of_integration_points; point_number++) {
                    rOutput[point_number] = 0.0;
                }
            }
        } else if (rVariable == VON_MISES_STRESS) {
            const SizeType number_of_nodes = r_geometry.size();
            const SizeType dimension = r_geometry.WorkingSpaceDimension();
            const SizeType strain_size = mConstitutiveLawVector[0]->GetStrainSize();

            KinematicVariables this_kinematic_variables(strain_size, dimension, number_of_nodes);
            ConstitutiveVariables this_constitutive_variables(strain_size);

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
                CalculateKinematicVariables(this_kinematic_variables, point_number, this->GetIntegrationMethod());

                // Compute material reponse
                CalculateConstitutiveVariables(this_kinematic_variables, this_constitutive_variables, Values, point_number, integration_points, GetStressMeasure());

                const Matrix stress_tensor = MathUtils<double>::StressVectorToTensor( this_constitutive_variables.StressVector );

                double sigma_equivalent = 0.0;

                if (dimension == 2) {
                    sigma_equivalent = std::pow((stress_tensor(0,0) - stress_tensor(1,1)), 2.0) +
                                                3*(stress_tensor(0,1) * stress_tensor(1,0));
                } else {
                    sigma_equivalent = 0.5*(std::pow((stress_tensor(0,0) - stress_tensor(1,1)), 2.0) +
                                            std::pow((stress_tensor(1,1) - stress_tensor(2,2)), 2.0) +
                                            std::pow((stress_tensor(2,2) - stress_tensor(0,0)), 2.0) +
                                                    6*(stress_tensor(0,1) * stress_tensor(1,0) +
                                                        stress_tensor(1,2) * stress_tensor(2,1) +
                                                        stress_tensor(2,0) * stress_tensor(0,2)));
                }

                if( sigma_equivalent < 0.0 )
                    rOutput[point_number] = 0.0;
                else
                    rOutput[point_number] = std::sqrt(sigma_equivalent);
            }
        } else {
            CalculateOnConstitutiveLaw(rVariable, rOutput, rCurrentProcessInfo);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::CalculateOnIntegrationPoints(
    const Variable<array_1d<double, 3>>& rVariable,
    std::vector<array_1d<double, 3>>& rOutput,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    const GeometryType::IntegrationPointsArrayType &integration_points = GetGeometry().IntegrationPoints(this->GetIntegrationMethod());

    const SizeType number_of_integration_points = integration_points.size();
    if ( rOutput.size() != number_of_integration_points )
        rOutput.resize( number_of_integration_points );

    if (mConstitutiveLawVector[0]->Has( rVariable)) {
        GetValueOnConstitutiveLaw(rVariable, rOutput);
    } else {
        if (rVariable == INTEGRATION_COORDINATES) {
            const SizeType number_of_nodes = GetGeometry().size();
            const SizeType dimension = GetGeometry().WorkingSpaceDimension();
            const SizeType strain_size = mConstitutiveLawVector[0]->GetStrainSize();

            KinematicVariables this_kinematic_variables(strain_size, dimension, number_of_nodes);

            for (IndexType point_number = 0; point_number < number_of_integration_points; ++point_number) {
                Point global_point;
                GetGeometry().GlobalCoordinates(global_point, integration_points[point_number]);

                rOutput[point_number] = global_point.Coordinates();
            }
        } else {
            CalculateOnConstitutiveLaw(rVariable, rOutput, rCurrentProcessInfo);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::CalculateOnIntegrationPoints(
    const Variable<array_1d<double, 6>>& rVariable,
    std::vector<array_1d<double, 6>>& rOutput,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    const GeometryType::IntegrationPointsArrayType &integration_points = GetGeometry().IntegrationPoints(this->GetIntegrationMethod());

    const SizeType number_of_integration_points = integration_points.size();
    if (rOutput.size() != number_of_integration_points)
        rOutput.resize(number_of_integration_points);

    if (mConstitutiveLawVector[0]->Has( rVariable)) {
        GetValueOnConstitutiveLaw(rVariable, rOutput);
    }  else {
        CalculateOnConstitutiveLaw(rVariable, rOutput, rCurrentProcessInfo);
    }
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::CalculateOnIntegrationPoints(
    const Variable<Vector>& rVariable,
    std::vector<Vector>& rOutput,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( this->GetIntegrationMethod() );

    const SizeType number_of_integration_points = integration_points.size();
    if ( rOutput.size() != number_of_integration_points )
        rOutput.resize( number_of_integration_points );

    if (mConstitutiveLawVector[0]->Has( rVariable)) {
        GetValueOnConstitutiveLaw(rVariable, rOutput);
    } else {
        if ( rVariable == INSITU_STRESS ) {
            const SizeType strain_size = mConstitutiveLawVector[0]->GetStrainSize();
            Vector strain_vector( strain_size );

            for ( IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number ) {
                if ( rOutput[point_number].size() != strain_vector.size() )
                    rOutput[point_number].resize( strain_vector.size(), false );

                rOutput[point_number] = mConstitutiveLawVector[point_number]->GetValue( INSITU_STRESS, rOutput[point_number] );
            }
        } else if ( rVariable == CAUCHY_STRESS_VECTOR || rVariable == PK2_STRESS_VECTOR ) {
            // Create and initialize element variables:
            const SizeType number_of_nodes = GetGeometry().size();
            const SizeType dimension = GetGeometry().WorkingSpaceDimension();
            const SizeType strain_size = mConstitutiveLawVector[0]->GetStrainSize();

            KinematicVariables this_kinematic_variables(strain_size, dimension, number_of_nodes);
            ConstitutiveVariables this_constitutive_variables(strain_size);

            // Create constitutive law parameters:
            ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

            // Set constitutive law flags:
            Flags& ConstitutiveLawOptions=Values.GetOptions();
            ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, UseElementProvidedStrain());
            ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
            ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);

            Values.SetStrainVector(this_constitutive_variables.StrainVector);

            // Reading integration points
            for ( IndexType point_number = 0; point_number < number_of_integration_points; ++point_number ) {
                // Compute element kinematics B, F, DN_DX ...
                CalculateKinematicVariables(this_kinematic_variables, point_number, this->GetIntegrationMethod());

                //call the constitutive law to update material variables
                if( rVariable == CAUCHY_STRESS_VECTOR) {
                    // Compute material reponse
                    CalculateConstitutiveVariables(this_kinematic_variables, this_constitutive_variables, Values, point_number, integration_points, ConstitutiveLaw::StressMeasure_Cauchy);
                } else {
                    // Compute material reponse
                    CalculateConstitutiveVariables(this_kinematic_variables, this_constitutive_variables, Values, point_number, integration_points,ConstitutiveLaw::StressMeasure_PK2);
                }

                if ( rOutput[point_number].size() != strain_size )
                    rOutput[point_number].resize( strain_size, false );

                rOutput[point_number] = this_constitutive_variables.StressVector;
            }
        } else if( rVariable == GREEN_LAGRANGE_STRAIN_VECTOR  || rVariable == ALMANSI_STRAIN_VECTOR ) {
            // Create and initialize element variables:
            const SizeType number_of_nodes = GetGeometry().size();
            const SizeType dimension = GetGeometry().WorkingSpaceDimension();
            const SizeType strain_size = mConstitutiveLawVector[0]->GetStrainSize();

            KinematicVariables this_kinematic_variables(strain_size, dimension, number_of_nodes);
            ConstitutiveVariables this_constitutive_variables(strain_size);

            // Create constitutive law parameters:
            ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

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
                CalculateKinematicVariables(this_kinematic_variables, point_number, this->GetIntegrationMethod());

                // Compute material reponse
                CalculateConstitutiveVariables(this_kinematic_variables, this_constitutive_variables, Values, point_number, integration_points, this_stress_measure);

                if ( rOutput[point_number].size() != strain_size)
                    rOutput[point_number].resize( strain_size, false );

                rOutput[point_number] = this_constitutive_variables.StrainVector;
            }
        } else {
            CalculateOnConstitutiveLaw(rVariable, rOutput, rCurrentProcessInfo);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::CalculateOnIntegrationPoints(
    const Variable<Matrix>& rVariable,
    std::vector<Matrix>& rOutput,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( this->GetIntegrationMethod() );
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();

    if ( rOutput.size() != integration_points.size() )
        rOutput.resize( integration_points.size() );

    if (mConstitutiveLawVector[0]->Has( rVariable)) {
        GetValueOnConstitutiveLaw(rVariable, rOutput);
    } else {
        if ( rVariable == CAUCHY_STRESS_TENSOR || rVariable == PK2_STRESS_TENSOR ) {
            std::vector<Vector> stress_vector;

            if( rVariable == CAUCHY_STRESS_TENSOR )
                this->CalculateOnIntegrationPoints( CAUCHY_STRESS_VECTOR, stress_vector, rCurrentProcessInfo );
            else
                this->CalculateOnIntegrationPoints( PK2_STRESS_VECTOR, stress_vector, rCurrentProcessInfo );

            // Loop integration points
            for ( IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number ) {
                if ( rOutput[point_number].size2() != dimension )
                    rOutput[point_number].resize( dimension, dimension, false );

                rOutput[point_number] = MathUtils<double>::StressVectorToTensor(stress_vector[point_number]);
            }
        }
        else if ( rVariable == GREEN_LAGRANGE_STRAIN_TENSOR  || rVariable == ALMANSI_STRAIN_TENSOR) {
            std::vector<Vector> strain_vector;
            if( rVariable == GREEN_LAGRANGE_STRAIN_TENSOR )
                CalculateOnIntegrationPoints( GREEN_LAGRANGE_STRAIN_VECTOR, strain_vector, rCurrentProcessInfo );
            else
                CalculateOnIntegrationPoints( ALMANSI_STRAIN_VECTOR, strain_vector, rCurrentProcessInfo );

            // Loop integration points
            for ( IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number ) {
                if ( rOutput[point_number].size2() != dimension )
                    rOutput[point_number].resize( dimension, dimension, false );

                rOutput[point_number] = MathUtils<double>::StrainVectorToTensor(strain_vector[point_number]);
            }
        } else if ( rVariable == CONSTITUTIVE_MATRIX ) {
            // Create and initialize element variables:
            const SizeType number_of_nodes = GetGeometry().size();
            const SizeType strain_size = mConstitutiveLawVector[0]->GetStrainSize();

            KinematicVariables this_kinematic_variables(strain_size, dimension, number_of_nodes);
            ConstitutiveVariables this_constitutive_variables(strain_size);

            // Create constitutive law parameters:
            ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

            // Set constitutive law flags:
            Flags& ConstitutiveLawOptions=Values.GetOptions();
            ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, UseElementProvidedStrain());
            ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, false);
            ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

            Values.SetStrainVector(this_constitutive_variables.StrainVector);
            Values.SetConstitutiveMatrix(this_constitutive_variables.D); //this is the output parameter

            // Reading integration points
            for ( IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number ) {
                // Compute element kinematics B, F, DN_DX ...
                CalculateKinematicVariables(this_kinematic_variables, point_number, this->GetIntegrationMethod());

                // Compute material reponse
                CalculateConstitutiveVariables(this_kinematic_variables, this_constitutive_variables, Values, point_number, integration_points, GetStressMeasure());

                if( rOutput[point_number].size2() != this_constitutive_variables.D.size2() )
                    rOutput[point_number].resize( this_constitutive_variables.D.size1() , this_constitutive_variables.D.size2() , false );

                rOutput[point_number] = this_constitutive_variables.D;
            }
        } else if ( rVariable == DEFORMATION_GRADIENT ) { // VARIABLE SET FOR TRANSFER PURPOUSES
            // Create and initialize element variables:
            const SizeType number_of_nodes = GetGeometry().size();
            const SizeType strain_size = mConstitutiveLawVector[0]->GetStrainSize();

            KinematicVariables this_kinematic_variables(strain_size, dimension, number_of_nodes);
            ConstitutiveVariables this_constitutive_variables(strain_size);

            // Create constitutive law parameters:
            ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

            // Reading integration points
            for ( IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number ) {
                // Compute element kinematics B, F, DN_DX ...
                CalculateKinematicVariables(this_kinematic_variables, point_number, this->GetIntegrationMethod());

                if( rOutput[point_number].size2() != this_kinematic_variables.F.size2() )
                    rOutput[point_number].resize( this_kinematic_variables.F.size1() , this_kinematic_variables.F.size2() , false );

                rOutput[point_number] = this_kinematic_variables.F;
            }
        }  else {
            CalculateOnConstitutiveLaw(rVariable, rOutput, rCurrentProcessInfo);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::SetValueOnIntegrationPoints(
    const Variable<bool>& rVariable,
    std::vector<bool>& rValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    if (mConstitutiveLawVector[0]->Has( rVariable)) {
        for ( IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number ) {
            mConstitutiveLawVector[point_number]->SetValue( rVariable,rValues[point_number], rCurrentProcessInfo);
        }
    } else {
        KRATOS_WARNING("BaseSolidElement") << "The variable " << rVariable << " is not implemented in the current ConstitutiveLaw" << std::endl;
    }
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::SetValueOnIntegrationPoints(
    const Variable<int>& rVariable,
    std::vector<int>& rValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    if (mConstitutiveLawVector[0]->Has( rVariable)) {
        for ( IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number ) {
            mConstitutiveLawVector[point_number]->SetValue( rVariable,rValues[point_number], rCurrentProcessInfo);
        }
    } else {
        KRATOS_WARNING("BaseSolidElement") << "The variable " << rVariable << " is not implemented in the current ConstitutiveLaw" << std::endl;
    }
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::SetValueOnIntegrationPoints(
    const Variable<double>& rVariable,
    std::vector<double>& rValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    if (mConstitutiveLawVector[0]->Has( rVariable)) {
        for ( IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number ) {
            mConstitutiveLawVector[point_number]->SetValue( rVariable,rValues[point_number], rCurrentProcessInfo);
        }
    } else {
        KRATOS_WARNING("BaseSolidElement") << "The variable " << rVariable << " is not implemented in the current ConstitutiveLaw" << std::endl;
    }
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::SetValueOnIntegrationPoints(
    const Variable<Vector>& rVariable,
    std::vector<Vector>& rValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    if (mConstitutiveLawVector[0]->Has( rVariable)) {
        for ( IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number ) {
            mConstitutiveLawVector[point_number]->SetValue( rVariable,rValues[point_number], rCurrentProcessInfo);
        }
    } else {
        KRATOS_WARNING("BaseSolidElement") << "The variable " << rVariable << " is not implemented in the current ConstitutiveLaw" << std::endl;
    }
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::SetValueOnIntegrationPoints(
    const Variable<ConstitutiveLaw::Pointer>& rVariable,
    std::vector<ConstitutiveLaw::Pointer>& rValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    if (rVariable == CONSTITUTIVE_LAW) {
        const SizeType integration_points_number = mConstitutiveLawVector.size();
        for ( IndexType point_number = 0; point_number < integration_points_number; ++point_number ) {
            mConstitutiveLawVector[point_number] = rValues[point_number];
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::SetValueOnIntegrationPoints(
    const Variable<array_1d<double, 3 > >& rVariable,
    std::vector<array_1d<double, 3 > > rValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    if (mConstitutiveLawVector[0]->Has( rVariable)) {
        for ( IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number ) {
            mConstitutiveLawVector[point_number]->SetValue( rVariable,rValues[point_number], rCurrentProcessInfo);
        }
    } else {
        KRATOS_WARNING("BaseSolidElement") << "The variable " << rVariable << " is not implemented in the current ConstitutiveLaw" << std::endl;
    }
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::SetValueOnIntegrationPoints(
    const Variable<array_1d<double, 6 > >& rVariable,
    std::vector<array_1d<double, 6 > > rValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    if (mConstitutiveLawVector[0]->Has( rVariable)) {
        for ( IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number ) {
            mConstitutiveLawVector[point_number]->SetValue( rVariable,rValues[point_number], rCurrentProcessInfo);
        }
    } else {
        KRATOS_WARNING("BaseSolidElement") << "The variable " << rVariable << " is not implemented in the current ConstitutiveLaw" << std::endl;
    }
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::SetValueOnIntegrationPoints(
    const Variable<Matrix>& rVariable,
    std::vector<Matrix>& rValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    if (mConstitutiveLawVector[0]->Has( rVariable)) {
        for ( IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number ) {
            mConstitutiveLawVector[point_number]->SetValue( rVariable,rValues[point_number], rCurrentProcessInfo);
        }
    } else {
        KRATOS_WARNING("BaseSolidElement") << "The variable " << rVariable << " is not implemented in the current ConstitutiveLaw" << std::endl;
    }
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::GetValueOnIntegrationPoints(
    const Variable<bool>& rVariable,
    std::vector<bool>& rValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::GetValueOnIntegrationPoints(
    const Variable<int>& rVariable,
    std::vector<int>& rValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::GetValueOnIntegrationPoints(
    const Variable<double>& rVariable,
    std::vector<double>& rValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::GetValueOnIntegrationPoints(
    const Variable<array_1d<double, 3>>& rVariable,
    std::vector<array_1d<double, 3>>& rValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::GetValueOnIntegrationPoints(
    const Variable<array_1d<double, 6>>& rVariable,
    std::vector<array_1d<double, 6>>& rValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::GetValueOnIntegrationPoints(
    const Variable<Vector>& rVariable,
    std::vector<Vector>& rValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::GetValueOnIntegrationPoints(
    const Variable<Matrix>& rVariable,
    std::vector<Matrix>& rValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::GetValueOnIntegrationPoints(
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

int  BaseSolidElement::Check( const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY;

    int check = Element::Check(rCurrentProcessInfo);

    // Basic check
    check = StructuralMechanicsElementUtilities::SolidElementCheck(*this, rCurrentProcessInfo, mConstitutiveLawVector);

    return check;

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::CalculateAll(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo,
    const bool CalculateStiffnessMatrixFlag,
    const bool CalculateResidualVectorFlag
    )
{
    KRATOS_ERROR << "You have called to the CalculateAll from the base class for solid elements" << std::endl;
}

//***********************************************************************
//***********************************************************************

double BaseSolidElement::GetIntegrationWeight(
    const GeometryType::IntegrationPointsArrayType& rThisIntegrationPoints,
    const IndexType PointNumber,
    const double detJ
    ) const
{
    return rThisIntegrationPoints[PointNumber].Weight() * detJ;
}

//***********************************************************************
//***********************************************************************

void BaseSolidElement::CalculateShapeGradientOfMassMatrix(MatrixType& rMassMatrix, ShapeParameter Deriv) const
{
    KRATOS_TRY;

    // Properties
    const auto& r_prop = GetProperties();

    // Geometry information
    const auto& r_geom = GetGeometry();
    SizeType dimension = r_geom.WorkingSpaceDimension();
    SizeType number_of_nodes = r_geom.size();
    SizeType mat_size = dimension * number_of_nodes;

    // Clear matrix
    if (rMassMatrix.size1() != mat_size || rMassMatrix.size2() != mat_size)
        rMassMatrix.resize( mat_size, mat_size, false );
    rMassMatrix = ZeroMatrix(mat_size, mat_size);

    // Checking density
    KRATOS_ERROR_IF_NOT(r_prop.Has(DENSITY)) << "DENSITY has to be provided for the calculation of the MassMatrix!" << std::endl;

    // Getting density
    const double density = StructuralMechanicsElementUtilities::GetDensityForMassMatrixComputation(*this);
    const double thickness = (dimension == 2 && r_prop.Has(THICKNESS)) ? r_prop[THICKNESS] : 1.0;

    const IntegrationMethod integration_method =
        IntegrationUtilities::GetIntegrationMethodForExactMassMatrixEvaluation(r_geom);
    const Matrix& Ncontainer = r_geom.ShapeFunctionsValues(integration_method);
    Matrix J0(dimension, dimension), DN_DX0_deriv;
    const auto& integration_points = r_geom.IntegrationPoints(integration_method);
    for (unsigned point_number = 0; point_number < integration_points.size(); ++point_number)
    {
        GeometryUtils::JacobianOnInitialConfiguration(
            r_geom, integration_points[point_number], J0);
        const Matrix& rDN_De = r_geom.ShapeFunctionsLocalGradients(integration_method)[point_number];
        GeometricalSensitivityUtility geometrical_sensitivity(J0, rDN_De);
        double detJ0_deriv;
        geometrical_sensitivity.CalculateSensitivity(Deriv, detJ0_deriv, DN_DX0_deriv);
        const double integration_weight =
            GetIntegrationWeight(integration_points, point_number, detJ0_deriv) * thickness;
        const Vector& rN = row(Ncontainer, point_number);

        for (unsigned i = 0; i < r_geom.size(); ++i)
        {
            const unsigned index_i = i * dimension;

            for (unsigned j = 0; j < r_geom.size(); ++j)
            {
                const unsigned index_j = j * dimension;
                const double NiNj_weight = rN[i] * rN[j] * integration_weight * density;

                for (unsigned k = 0; k < dimension; ++k)
                    rMassMatrix(index_i + k, index_j + k) += NiNj_weight;
            }
        }
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::CalculateKinematicVariables(
    KinematicVariables& rThisKinematicVariables,
    const IndexType PointNumber,
    const GeometryType::IntegrationMethod& rIntegrationMethod
    )
{
    KRATOS_ERROR << "You have called to the CalculateKinematicVariables from the base class for solid elements" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::CalculateConstitutiveVariables(
    KinematicVariables& rThisKinematicVariables,
    ConstitutiveVariables& rThisConstitutiveVariables,
    ConstitutiveLaw::Parameters& rValues,
    const IndexType PointNumber,
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints,
    const ConstitutiveLaw::StressMeasure ThisStressMeasure
    )
{
    // Setting the variables for the CL
    SetConstitutiveVariables(rThisKinematicVariables, rThisConstitutiveVariables, rValues, PointNumber, IntegrationPoints);

    // Actually do the computations in the ConstitutiveLaw
    mConstitutiveLawVector[PointNumber]->CalculateMaterialResponse(rValues, ThisStressMeasure); //here the calculations are actually done
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::SetConstitutiveVariables(
    KinematicVariables& rThisKinematicVariables,
    ConstitutiveVariables& rThisConstitutiveVariables,
    ConstitutiveLaw::Parameters& rValues,
    const IndexType PointNumber,
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints
    )
{
    // Here we essentially set the input parameters
    rValues.SetShapeFunctionsValues(rThisKinematicVariables.N); // shape functions
    rValues.SetDeterminantF(rThisKinematicVariables.detF); // Assuming the determinant is computed somewhere else
    rValues.SetDeformationGradientF(rThisKinematicVariables.F); //F computed somewhere else

    // Here we set the space on which the results shall be written
    rValues.SetConstitutiveMatrix(rThisConstitutiveVariables.D); // Assuming the determinant is computed somewhere else
    rValues.SetStressVector(rThisConstitutiveVariables.StressVector); //F computed somewhere else
}

/***********************************************************************************/
/***********************************************************************************/

Matrix& BaseSolidElement::CalculateDeltaDisplacement(Matrix& DeltaDisplacement) const
{
    KRATOS_TRY

    const SizeType number_of_nodes = GetGeometry().PointsNumber();
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();

    DeltaDisplacement.resize(number_of_nodes , dimension, false);

    for ( IndexType i_node = 0; i_node < number_of_nodes; i_node++ ) {
        const array_1d<double, 3 >& current_displacement  = GetGeometry()[i_node].FastGetSolutionStepValue(DISPLACEMENT);
        const array_1d<double, 3 >& previous_displacement = GetGeometry()[i_node].FastGetSolutionStepValue(DISPLACEMENT,1);

        for ( IndexType j_dim = 0; j_dim < dimension; ++j_dim )
            DeltaDisplacement(i_node, j_dim) = current_displacement[j_dim] - previous_displacement[j_dim];
    }

    return DeltaDisplacement;

    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

double BaseSolidElement::CalculateDerivativesOnReferenceConfiguration(
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

/***********************************************************************************/
/***********************************************************************************/

double BaseSolidElement::CalculateDerivativesOnCurrentConfiguration(
    Matrix& rJ,
    Matrix& rInvJ,
    Matrix& rDN_DX,
    const IndexType PointNumber,
    IntegrationMethod ThisIntegrationMethod
    ) const
{
    double detJ;
    rJ = GetGeometry().Jacobian( rJ, PointNumber, ThisIntegrationMethod );
    const Matrix& DN_De = GetGeometry().ShapeFunctionsLocalGradients(ThisIntegrationMethod)[PointNumber];
    MathUtils<double>::InvertMatrix( rJ, rInvJ, detJ );
    GeometryUtils::ShapeFunctionsGradients(DN_De, rInvJ, rDN_DX);
    return detJ;
}

/***********************************************************************************/
/***********************************************************************************/

array_1d<double, 3> BaseSolidElement::GetBodyForce(
    const GeometryType::IntegrationPointsArrayType& rIntegrationPoints,
    const IndexType PointNumber
    ) const
{
    return StructuralMechanicsElementUtilities::GetBodyForce(*this, rIntegrationPoints, PointNumber);
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::CalculateAndAddKm(
    MatrixType& rLeftHandSideMatrix,
    const Matrix& B,
    const Matrix& D,
    const double IntegrationWeight
    ) const
{
    KRATOS_TRY

    noalias( rLeftHandSideMatrix ) += IntegrationWeight * prod( trans( B ), Matrix(prod(D, B)));

    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::CalculateAndAddKg(
    MatrixType& rLeftHandSideMatrix,
    const Matrix& DN_DX,
    const Vector& StressVector,
    const double IntegrationWeight
    ) const
{
    KRATOS_TRY

    const SizeType dimension = GetGeometry().WorkingSpaceDimension();
    Matrix stress_tensor = MathUtils<double>::StressVectorToTensor( StressVector );
    Matrix reduced_Kg = prod( DN_DX, IntegrationWeight * Matrix( prod( stress_tensor, trans( DN_DX ) ) ) ); //to be optimized
    MathUtils<double>::ExpandAndAddReducedMatrix( rLeftHandSideMatrix, reduced_Kg, dimension );

    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::CalculateAndAddResidualVector(
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

void BaseSolidElement::CalculateAndAddExtForceContribution(
    const Vector& rN,
    const ProcessInfo& rCurrentProcessInfo,
    const array_1d<double, 3>& rBodyForce,
    VectorType& rRightHandSideVector,
    const double Weight
    ) const
{
    KRATOS_TRY;

    const SizeType number_of_nodes = GetGeometry().PointsNumber();
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();

    for ( IndexType i = 0; i < number_of_nodes; ++i ) {
        const SizeType index = dimension * i;

        for ( IndexType j = 0; j < dimension; ++j )
            rRightHandSideVector[index + j] += Weight * rN[i] * rBodyForce[j];
    }

    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::CalculateLumpedMassVector(VectorType& rMassVector) const
{
    KRATOS_TRY;

    const auto& r_geom = GetGeometry();
    const auto& r_prop = GetProperties();
    const SizeType dimension = r_geom.WorkingSpaceDimension();
    const SizeType number_of_nodes = r_geom.size();
    const SizeType mat_size = dimension * number_of_nodes;

    // Clear matrix
    if (rMassVector.size() != mat_size)
        rMassVector.resize( mat_size, false );

    const double density = StructuralMechanicsElementUtilities::GetDensityForMassMatrixComputation(*this);
    const double thickness = (dimension == 2 && r_prop.Has(THICKNESS)) ? r_prop[THICKNESS] : 1.0;

    // LUMPED MASS MATRIX
    const double total_mass = GetGeometry().DomainSize() * density * thickness;

    Vector lumping_factors;
    lumping_factors = GetGeometry().LumpingFactors( lumping_factors );

    for ( IndexType i = 0; i < number_of_nodes; ++i ) {
        const double temp = lumping_factors[i] * total_mass;
        for ( IndexType j = 0; j < dimension; ++j ) {
            IndexType index = i * dimension + j;
            rMassVector[index] = temp;
        }
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::CalculateDampingMatrixWithLumpedMass(
    MatrixType& rDampingMatrix,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    unsigned int number_of_nodes = GetGeometry().size();
    unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    // Resizing as needed the LHS
    unsigned int mat_size = number_of_nodes * dimension;

    if ( rDampingMatrix.size1() != mat_size )
        rDampingMatrix.resize( mat_size, mat_size, false );

    noalias( rDampingMatrix ) = ZeroMatrix( mat_size, mat_size );

    // 1.-Get Damping Coeffitients (RAYLEIGH_ALPHA, RAYLEIGH_BETA)
    double alpha = 0.0;
    if( GetProperties().Has(RAYLEIGH_ALPHA) )
        alpha = GetProperties()[RAYLEIGH_ALPHA];
    else if( rCurrentProcessInfo.Has(RAYLEIGH_ALPHA) )
        alpha = rCurrentProcessInfo[RAYLEIGH_ALPHA];

    double beta  = 0.0;
    if( GetProperties().Has(RAYLEIGH_BETA) )
        beta = GetProperties()[RAYLEIGH_BETA];
    else if( rCurrentProcessInfo.Has(RAYLEIGH_BETA) )
        beta = rCurrentProcessInfo[RAYLEIGH_BETA];

    // Compose the Damping Matrix:
    // Rayleigh Damping Matrix: alpha*M + beta*K

    // 2.-Calculate mass matrix:
    if (alpha > std::numeric_limits<double>::epsilon()) {
        VectorType temp_vector(mat_size);
        CalculateLumpedMassVector(temp_vector);
        for (IndexType i = 0; i < mat_size; ++i)
            rDampingMatrix(i, i) += alpha * temp_vector[i];
    }

    // 3.-Calculate StiffnessMatrix:
    if (beta > std::numeric_limits<double>::epsilon()) {
        MatrixType stiffness_matrix( mat_size, mat_size );
        VectorType residual_vector( mat_size );

        this->CalculateAll(stiffness_matrix, residual_vector, rCurrentProcessInfo, true, false);

        noalias( rDampingMatrix ) += beta  * stiffness_matrix;
    }

    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Element );
    int IntMethod = int(this->GetIntegrationMethod());
    rSerializer.save("IntegrationMethod",IntMethod);
    rSerializer.save("ConstitutiveLawVector", mConstitutiveLawVector);
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Element );
    int IntMethod;
    rSerializer.load("IntegrationMethod",IntMethod);
    mThisIntegrationMethod = IntegrationMethod(IntMethod);
    rSerializer.load("ConstitutiveLawVector", mConstitutiveLawVector);
}
} // Namespace Kratos
