// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Reza Najian Asl
//                   Shahed Rezaei
//

// System includes

// External includes

// Project includes
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h"

// Application includes
#include "custom_elements/base_phase_field_fracture_solid_element.h"
#include "structural_mechanics_application_variables.h"
#include "custom_utilities/structural_mechanics_element_utilities.h"

namespace Kratos
{
void BasePhaseFieldFractureSolidElement::Initialize(const ProcessInfo& rCurrentProcessInfo)
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
            KRATOS_WARNING("BasePhaseFieldFractureSolidElement") << "Integration order " << integration_order << " is not available, using default integration order for the geometry" << std::endl;
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

    //Initialisation of the history variable field H   
    if (mHisVarVector.size() != integration_points.size())
        mHisVarVector.resize( integration_points.size());
    for(int i=0;i<mHisVarVector.size();i++)
        mHisVarVector[i] = 0.0;
    
    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

void BasePhaseFieldFractureSolidElement::InitializeSolutionStep( ProcessInfo& rCurrentProcessInfo )
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

void BasePhaseFieldFractureSolidElement::InitializeNonLinearIteration( ProcessInfo& rCurrentProcessInfo )
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

void BasePhaseFieldFractureSolidElement::FinalizeNonLinearIteration( ProcessInfo& rCurrentProcessInfo )
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

void BasePhaseFieldFractureSolidElement::FinalizeSolutionStep( ProcessInfo& rCurrentProcessInfo )
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

void BasePhaseFieldFractureSolidElement::InitializeMaterial()
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

ConstitutiveLaw::StressMeasure BasePhaseFieldFractureSolidElement::GetStressMeasure() const
{
    return ConstitutiveLaw::StressMeasure_PK2;
}

/***********************************************************************************/
/***********************************************************************************/

bool BasePhaseFieldFractureSolidElement::UseElementProvidedStrain() const
{
    return false;
}

/***********************************************************************************/
/***********************************************************************************/

void BasePhaseFieldFractureSolidElement::ResetConstitutiveLaw()
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

Element::Pointer BasePhaseFieldFractureSolidElement::Clone (
    IndexType NewId,
    NodesArrayType const& rThisNodes
    ) const
{
    KRATOS_TRY

    KRATOS_WARNING("BasePhaseFieldFractureSolidElement") << " Call BasePhaseFieldFractureSolidElement (base class) Clone " << std::endl;

    BasePhaseFieldFractureSolidElement::Pointer p_new_elem = Kratos::make_intrusive<BasePhaseFieldFractureSolidElement>(NewId, GetGeometry().Create(rThisNodes), pGetProperties());
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

void BasePhaseFieldFractureSolidElement::EquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo
    ) const
{
    KRATOS_TRY;

    const SizeType number_of_nodes = GetGeometry().size();
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();

    if (rResult.size() != dimension * number_of_nodes + number_of_nodes)
        rResult.resize(dimension * number_of_nodes + number_of_nodes,false);

    const auto& r_geom = GetGeometry();

    if(dimension == 2) {
        for (IndexType i = 0; i < number_of_nodes; ++i) {
            const SizeType index = i * 3;
            const NodeType& i_node = r_geom[i];
            rResult[index] = i_node.GetDof(DISPLACEMENT_X).EquationId();
            rResult[index + 1] = i_node.GetDof(DISPLACEMENT_Y).EquationId();
            rResult[index + 2] = i_node.GetDof(DAMAGE).EquationId();
        }
    } else {
        for (IndexType i = 0; i < number_of_nodes; ++i) {
            const SizeType index = i * 4;
            const NodeType& i_node = r_geom[i];
            rResult[index] = i_node.GetDof(DISPLACEMENT_X).EquationId();
            rResult[index + 1] = i_node.GetDof(DISPLACEMENT_Y).EquationId();
            rResult[index + 2] = i_node.GetDof(DISPLACEMENT_Z).EquationId();
            rResult[index + 3] = i_node.GetDof(DAMAGE).EquationId();
        }
    }

    KRATOS_CATCH("")
};

/***********************************************************************************/
/***********************************************************************************/

void BasePhaseFieldFractureSolidElement::GetDofList(
    DofsVectorType& rElementalDofList,
    const ProcessInfo& rCurrentProcessInfo
    ) const
{
    KRATOS_TRY;

    const SizeType number_of_nodes = GetGeometry().size();
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();
    rElementalDofList.resize(0);
    rElementalDofList.reserve(dimension*number_of_nodes + number_of_nodes);

    if(dimension == 2) {
        for (IndexType i = 0; i < number_of_nodes; ++i) {
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_X));
            rElementalDofList.push_back( GetGeometry()[i].pGetDof(DISPLACEMENT_Y));
            rElementalDofList.push_back( GetGeometry()[i].pGetDof(DAMAGE));
        }
    } else {
        for (IndexType i = 0; i < number_of_nodes; ++i) {
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_X));
            rElementalDofList.push_back( GetGeometry()[i].pGetDof(DISPLACEMENT_Y));
            rElementalDofList.push_back( GetGeometry()[i].pGetDof(DISPLACEMENT_Z));
            rElementalDofList.push_back( GetGeometry()[i].pGetDof(DAMAGE));
        }
    }

    KRATOS_CATCH("")
};

/***********************************************************************************/
/***********************************************************************************/

void BasePhaseFieldFractureSolidElement::GetValuesVector(
    Vector& rValues,
    int Step
    ) const
{
    const SizeType number_of_nodes = GetGeometry().size();
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();
    const SizeType mat_size = number_of_nodes * dimension + number_of_nodes;
    if (rValues.size() != mat_size)
        rValues.resize(mat_size, false);
    for (IndexType i = 0; i < number_of_nodes; ++i)
    {
        const array_1d<double, 3 >& displacement = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT, Step);
        const SizeType index = i * (dimension + 1);
        for(unsigned int k = 0; k < dimension; ++k)
        {
            rValues[index + k] = displacement[k];
        }

        const double& damagnode = GetGeometry()[i].FastGetSolutionStepValue(DAMAGE, Step);
        rValues[index + dimension] = damagnode;

    }
}

/***********************************************************************************/
/***********************************************************************************/

void BasePhaseFieldFractureSolidElement::GetDisplacementValuesVector(
    Vector& rValues,
    int Step
    ) const
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

void BasePhaseFieldFractureSolidElement::GetPhaseFieldValuesVector(
    Vector& rValues,
    int Step
    ) const
{
    const SizeType number_of_nodes = GetGeometry().size();
    const SizeType mat_size = number_of_nodes;
    if (rValues.size() != mat_size)
        rValues.resize(mat_size, false);
    for (IndexType i = 0; i < number_of_nodes; ++i)
    {   
        const SizeType index = i;
        const double& damagnode = GetGeometry()[i].FastGetSolutionStepValue(DAMAGE, Step);
        rValues[index] = damagnode;
    }
}

/***********************************************************************************/
/***********************************************************************************/

void BasePhaseFieldFractureSolidElement::CalculateLocalSystem(
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

void BasePhaseFieldFractureSolidElement::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
                                             ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;
    VectorType RHS;
    CalculateLocalSystem(rLeftHandSideMatrix, RHS, rCurrentProcessInfo);
    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void BasePhaseFieldFractureSolidElement::CalculateRightHandSide(
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

void BasePhaseFieldFractureSolidElement::CalculateOnIntegrationPoints(
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

void BasePhaseFieldFractureSolidElement::CalculateOnIntegrationPoints(
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

void BasePhaseFieldFractureSolidElement::CalculateOnIntegrationPoints(
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

void BasePhaseFieldFractureSolidElement::CalculateOnIntegrationPoints(
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

void BasePhaseFieldFractureSolidElement::CalculateOnIntegrationPoints(
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

void BasePhaseFieldFractureSolidElement::CalculateOnIntegrationPoints(
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

void BasePhaseFieldFractureSolidElement::CalculateOnIntegrationPoints(
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

void BasePhaseFieldFractureSolidElement::SetValuesOnIntegrationPoints(
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
        KRATOS_WARNING("BasePhaseFieldFractureSolidElement") << "The variable " << rVariable << " is not implemented in the current ConstitutiveLaw" << std::endl;
    }
}

/***********************************************************************************/
/***********************************************************************************/

void BasePhaseFieldFractureSolidElement::SetValuesOnIntegrationPoints(
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
        KRATOS_WARNING("BasePhaseFieldFractureSolidElement") << "The variable " << rVariable << " is not implemented in the current ConstitutiveLaw" << std::endl;
    }
}

/***********************************************************************************/
/***********************************************************************************/

void BasePhaseFieldFractureSolidElement::SetValuesOnIntegrationPoints(
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
        KRATOS_WARNING("BasePhaseFieldFractureSolidElement") << "The variable " << rVariable << " is not implemented in the current ConstitutiveLaw" << std::endl;
    }
}

/***********************************************************************************/
/***********************************************************************************/

void BasePhaseFieldFractureSolidElement::SetValuesOnIntegrationPoints(
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
        KRATOS_WARNING("BasePhaseFieldFractureSolidElement") << "The variable " << rVariable << " is not implemented in the current ConstitutiveLaw" << std::endl;
    }
}

/***********************************************************************************/
/***********************************************************************************/

void BasePhaseFieldFractureSolidElement::SetValuesOnIntegrationPoints(
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

void BasePhaseFieldFractureSolidElement::SetValuesOnIntegrationPoints(
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
        KRATOS_WARNING("BasePhaseFieldFractureSolidElement") << "The variable " << rVariable << " is not implemented in the current ConstitutiveLaw" << std::endl;
    }
}

/***********************************************************************************/
/***********************************************************************************/

void BasePhaseFieldFractureSolidElement::SetValuesOnIntegrationPoints(
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
        KRATOS_WARNING("BasePhaseFieldFractureSolidElement") << "The variable " << rVariable << " is not implemented in the current ConstitutiveLaw" << std::endl;
    }
}

/***********************************************************************************/
/***********************************************************************************/

void BasePhaseFieldFractureSolidElement::SetValuesOnIntegrationPoints(
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
        KRATOS_WARNING("BasePhaseFieldFractureSolidElement") << "The variable " << rVariable << " is not implemented in the current ConstitutiveLaw" << std::endl;
    }
}


/***********************************************************************************/
/***********************************************************************************/

void BasePhaseFieldFractureSolidElement::GetValueOnIntegrationPoints(
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

int  BasePhaseFieldFractureSolidElement::Check( const ProcessInfo& rCurrentProcessInfo ) const
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

void BasePhaseFieldFractureSolidElement::CalculateAll(
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

double BasePhaseFieldFractureSolidElement::GetIntegrationWeight(
    const GeometryType::IntegrationPointsArrayType& rThisIntegrationPoints,
    const IndexType PointNumber,
    const double detJ
    ) const
{
    return rThisIntegrationPoints[PointNumber].Weight() * detJ;
}

/***********************************************************************************/
/***********************************************************************************/

void BasePhaseFieldFractureSolidElement::CalculateKinematicVariables(
    KinematicVariables& rThisKinematicVariables,
    const IndexType PointNumber,
    const GeometryType::IntegrationMethod& rIntegrationMethod
    )
{
    KRATOS_ERROR << "You have called to the CalculateKinematicVariables from the base class for solid elements" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

void BasePhaseFieldFractureSolidElement::CalculateConstitutiveVariables(
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

void BasePhaseFieldFractureSolidElement::SetConstitutiveVariables(
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

double BasePhaseFieldFractureSolidElement::CalculateDerivativesOnReferenceConfiguration(
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

double BasePhaseFieldFractureSolidElement::CalculateDerivativesOnCurrentConfiguration(
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

array_1d<double, 3> BasePhaseFieldFractureSolidElement::GetBodyForce(
    const GeometryType::IntegrationPointsArrayType& rIntegrationPoints,
    const IndexType PointNumber
    ) const
{
    return StructuralMechanicsElementUtilities::GetBodyForce(*this, rIntegrationPoints, PointNumber);
}

/***********************************************************************************/
/***********************************************************************************/

void BasePhaseFieldFractureSolidElement::CalculateAndAddKm(
    MatrixType& rLeftHandSideMatrix,
    const KinematicVariables& rThisKinematicVariables,
    const ConstitutiveVariables& rThisConstitutiveVariables,
    const double IntegrationWeight
    ) const
{
    KRATOS_TRY

    const SizeType dimension = GetGeometry().WorkingSpaceDimension();
    const SizeType number_of_nodes = GetGeometry().PointsNumber();
    const SizeType mat_size = number_of_nodes * (dimension + 1); //damage dimension is included

    //momentum stiffness equation
    Matrix Kuu = IntegrationWeight * prod( trans( rThisKinematicVariables.B ), Matrix(prod(rThisConstitutiveVariables.D, rThisKinematicVariables.B))); //displ res contri.
    
    //damage stiffness equation
    const double gc = (GetProperties().Has( CRITICAL_ENERGY_RELEASE_RATE ) == true) ? this->GetProperties()[CRITICAL_ENERGY_RELEASE_RATE] : 10.0;
    const double lc = (GetProperties().Has( PHASE_FIELD_LENGTH ) == true) ? this->GetProperties()[PHASE_FIELD_LENGTH] : 0.1;
    double Hi = rThisConstitutiveVariables.MaxUndamagedElasticEnergy;
    Matrix Kdd = ZeroMatrix( number_of_nodes, number_of_nodes ); //damage res contri.
    Kdd = IntegrationWeight * gc * ((1.0/lc) * outer_prod(rThisKinematicVariables.N,rThisKinematicVariables.N) + lc * prod(rThisKinematicVariables.DN_DX,trans(rThisKinematicVariables.DN_DX)) );
    noalias(Kdd) += IntegrationWeight * 2 * Hi * outer_prod(rThisKinematicVariables.N,rThisKinematicVariables.N);

    Matrix Kt = ZeroMatrix( mat_size, mat_size );
    for ( IndexType i = 0; i < Kuu.size1(); ++i )
        for ( IndexType j = 0; j < Kuu.size2(); ++j )       
            Kt(i+(i/dimension),j+(j/dimension)) += Kuu(i,j);

    for ( IndexType i = 0; i < Kdd.size1(); ++i )
        for ( IndexType j = 0; j < Kdd.size2(); ++j )
            Kt(i*(dimension+1) + dimension,j*(dimension+1) + dimension) += Kdd(i,j);                        
               

    noalias( rLeftHandSideMatrix ) += Kt;

    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

void BasePhaseFieldFractureSolidElement::CalculateAndAddResidualVector(
    VectorType& rRightHandSideVector,
    const KinematicVariables& rThisKinematicVariables,
    const ConstitutiveVariables& rThisConstitutiveVariables,
    const ProcessInfo& rCurrentProcessInfo,
    const array_1d<double, 3>& rBodyForce,
    const Vector& rStressVector,
    const double IntegrationWeight
    ) const
{
    KRATOS_TRY

    const SizeType number_of_nodes = GetGeometry().PointsNumber();
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();
    const SizeType mat_size = number_of_nodes * (dimension + 1); //damage dimension is included

    // Operation performed: rRightHandSideVector += ExtForce * IntegrationWeight
    Vector forceCont = ZeroVector( number_of_nodes * dimension); //damage res contr.
    this->CalculateAndAddExtForceContribution( rThisKinematicVariables.N, rCurrentProcessInfo, rBodyForce, forceCont, IntegrationWeight );

    // Operation performed: = IntForce * IntegrationWeight
    Vector rHu = forceCont - IntegrationWeight * prod( trans( rThisKinematicVariables.B ), rStressVector ); //displacement res contr.
    //Vector rHd = ZeroVector( number_of_nodes ); //damage res contr.
    //damage params

    const double gc = (GetProperties().Has( CRITICAL_ENERGY_RELEASE_RATE ) == true) ? this->GetProperties()[CRITICAL_ENERGY_RELEASE_RATE] : 10.0;
    const double lc = (GetProperties().Has( PHASE_FIELD_LENGTH ) == true) ? this->GetProperties()[PHASE_FIELD_LENGTH] : 0.1;
    double Hi = rThisConstitutiveVariables.MaxUndamagedElasticEnergy;
    double D = inner_prod(rThisKinematicVariables.N,rThisKinematicVariables.PhaseField);
    double dfDdD = 2 * (D-1);    
    Matrix KDD = IntegrationWeight * gc * ((1.0/lc) * outer_prod(rThisKinematicVariables.N,rThisKinematicVariables.N) + lc * prod(rThisKinematicVariables.DN_DX,trans(rThisKinematicVariables.DN_DX)) );
    Vector rHd = -1 * (prod(KDD, rThisKinematicVariables.PhaseField) + IntegrationWeight * dfDdD * Hi * rThisKinematicVariables.N);

    Vector rHt = ZeroVector( mat_size ); // total RHS composed of rHu and rHt

    for ( IndexType i = 0; i < rHu.size(); ++i )       
        rHt[i+(i/dimension)] += rHu[i];    

    for ( IndexType i = 0; i < rHd.size(); ++i )
        rHt[i*(dimension+1) + dimension] += rHd[i]; 
    
    noalias( rRightHandSideVector ) += rHt;

    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

void BasePhaseFieldFractureSolidElement::CalculateAndAddExtForceContribution(
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
        const SizeType index = (dimension + 1) * i;

        for ( IndexType j = 0; j < dimension; ++j )
            rRightHandSideVector[index + j] += Weight * rN[i] * rBodyForce[j];
    }

    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

void BasePhaseFieldFractureSolidElement::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Element );
    int IntMethod = int(this->GetIntegrationMethod());
    rSerializer.save("IntegrationMethod",IntMethod);
    rSerializer.save("ConstitutiveLawVector", mConstitutiveLawVector);
}

/***********************************************************************************/
/***********************************************************************************/

void BasePhaseFieldFractureSolidElement::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Element );
    int IntMethod;
    rSerializer.load("IntegrationMethod",IntMethod);
    mThisIntegrationMethod = IntegrationMethod(IntMethod);
    rSerializer.load("ConstitutiveLawVector", mConstitutiveLawVector);
}
} // Namespace Kratos
