// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Riccardo Rossi
//                   Vicente Mataix Ferrándiz
//

// System includes

// External includes

// Project includes
#include "includes/checks.h"
#include "includes/define.h"
#include "custom_elements/base_solid_element.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h"
#include "includes/constitutive_law.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos
{
void BaseSolidElement::Initialize()
{
    KRATOS_TRY

    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints(  );

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
    for ( unsigned int point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number ) {
        mConstitutiveLawVector[point_number]->InitializeSolutionStep( GetProperties(),
        GetGeometry(),
        row( GetGeometry().ShapeFunctionsValues(  ), point_number ),
        rCurrentProcessInfo
        );
    }
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::InitializeNonLinearIteration( ProcessInfo& rCurrentProcessInfo )
{
    // TODO: Add somethig if necessary
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::FinalizeNonLinearIteration( ProcessInfo& rCurrentProcessInfo )
{
    // TODO: Add somethig if necessary
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::FinalizeSolutionStep( ProcessInfo& rCurrentProcessInfo )
{
    const unsigned int strain_size = mConstitutiveLawVector[0]->GetStrainSize();
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
    for ( unsigned int point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number ) {
        // Call the constitutive law to update material variables
        mConstitutiveLawVector[point_number]->FinalizeMaterialResponse(Values, GetStressMeasure());

        mConstitutiveLawVector[point_number]->FinalizeSolutionStep( GetProperties(),
        GetGeometry(),
        row( GetGeometry().ShapeFunctionsValues(  ), point_number ),
        rCurrentProcessInfo
        );
    }
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::InitializeMaterial()
{
    KRATOS_TRY

    if ( GetProperties()[CONSTITUTIVE_LAW] != nullptr ) {
        for ( unsigned int point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number ) {
            mConstitutiveLawVector[point_number] = GetProperties()[CONSTITUTIVE_LAW]->Clone();
            mConstitutiveLawVector[point_number]->InitializeMaterial( GetProperties(),
            GetGeometry(),
            row( GetGeometry().ShapeFunctionsValues(  ), point_number )
            );
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

bool BaseSolidElement::UseElementProvidedStrain()
{
    return false;
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::ResetConstitutiveLaw()
{
    KRATOS_TRY

    if ( GetProperties()[CONSTITUTIVE_LAW] != nullptr ) {
        for ( unsigned int point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number )
            mConstitutiveLawVector[point_number]->ResetMaterial( GetProperties(),
            GetGeometry(),
            row( GetGeometry().ShapeFunctionsValues(  ),                                                                                                                             point_number )
            );
    }

    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::EquationIdVector(
    EquationIdVectorType& rResult,
    ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    if (rResult.size() != dimension * number_of_nodes)
        rResult.resize(dimension * number_of_nodes,false);

    const unsigned int pos = this->GetGeometry()[0].GetDofPosition(DISPLACEMENT_X);

    if(dimension == 2) {
        for (unsigned int i = 0; i < number_of_nodes; ++i) {
            const unsigned int index = i * 2;
            rResult[index] = GetGeometry()[i].GetDof(DISPLACEMENT_X,pos).EquationId();
            rResult[index + 1] = GetGeometry()[i].GetDof(DISPLACEMENT_Y,pos+1).EquationId();
        }
    } else {
        for (unsigned int i = 0; i < number_of_nodes; ++i) {
            const unsigned int index = i * 3;
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

    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    rElementalDofList.resize(0);
    rElementalDofList.reserve(dimension*number_of_nodes);

    if(dimension == 2) {
        for (unsigned int i = 0; i < number_of_nodes; ++i) {
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_X));
            rElementalDofList.push_back( GetGeometry()[i].pGetDof(DISPLACEMENT_Y));
        }
    } else {
        for (unsigned int i = 0; i < number_of_nodes; ++i) {
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
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    const unsigned int mat_size = number_of_nodes * dimension;
    if (rValues.size() != mat_size)
        rValues.resize(mat_size, false);
    for (unsigned int i = 0; i < number_of_nodes; ++i)
    {
        const array_1d<double, 3 >& displacement = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT, Step);
        const unsigned int index = i * dimension;
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
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    const unsigned int mat_size = number_of_nodes * dimension;
    if (rValues.size() != mat_size)
        rValues.resize(mat_size, false);
    for (unsigned int i = 0; i < number_of_nodes; ++i) {
        const array_1d<double, 3 >& velocity = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, Step);
        const unsigned int index = i * dimension;
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
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    const unsigned int mat_size = number_of_nodes * dimension;
    if (rValues.size() != mat_size)
        rValues.resize(mat_size, false);
    for (unsigned int i = 0; i < number_of_nodes; ++i) {
        const array_1d<double, 3 >& acceleration = GetGeometry()[i].FastGetSolutionStepValue(ACCELERATION, Step);
        const unsigned int index = i * dimension;
        for(unsigned int k = 0; k < dimension; ++k)
            rValues[index + k] = acceleration[k];
    }
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

void BaseSolidElement::CalculateMassMatrix(
    MatrixType& rMassMatrix,
    ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    const auto& r_geom = GetGeometry();
    const auto& r_prop = GetProperties();
    unsigned int dimension = r_geom.WorkingSpaceDimension();
    unsigned int number_of_nodes = r_geom.size();
    unsigned int mat_size = dimension * number_of_nodes;

    rMassMatrix = ZeroMatrix( mat_size, mat_size );

    Matrix J0(dimension,dimension);

    IntegrationMethod integration_method = IntegrationUtilities::GetIntegrationMethodForExactMassMatrixEvaluation(r_geom);
    const GeometryType::IntegrationPointsArrayType& integration_points = r_geom.IntegrationPoints( integration_method );
    const Matrix& Ncontainer = r_geom.ShapeFunctionsValues(integration_method);

    KRATOS_ERROR_IF_NOT(r_prop.Has( DENSITY ))
        << "DENSITY has to be provided for the calculation of the MassMatrix!" << std::endl;

    const double density = r_prop[DENSITY];
    const double thickness = (dimension == 2 && r_prop.Has(THICKNESS)) ? r_prop[THICKNESS] : 1.0;

    for ( unsigned int point_number = 0; point_number < integration_points.size(); ++point_number ) {
        GeometryUtils::JacobianOnInitialConfiguration(
            r_geom, integration_points[point_number], J0);
        const double detJ0 = MathUtils<double>::DetMat(J0);
        const double integration_weight =
            GetIntegrationWeight(integration_points, point_number, detJ0) * thickness;
        const Vector& rN = row(Ncontainer,point_number);

        for ( unsigned int i = 0; i < number_of_nodes; ++i ) {
            const unsigned int index_i = i * dimension;

            for ( unsigned int j = 0; j < number_of_nodes; ++j ) {
                const unsigned int index_j = j * dimension;
                const double NiNj_weight = rN[i] * rN[j] * integration_weight * density;

                for ( unsigned int k = 0; k < dimension; ++k )
                    rMassMatrix( index_i + k, index_j + k ) += NiNj_weight;
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
    KRATOS_TRY;

    unsigned int number_of_nodes = GetGeometry().size();
    unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    // Resizing as needed the LHS
    unsigned int mat_size = number_of_nodes * dimension;

    if ( rDampingMatrix.size1() != mat_size )
        rDampingMatrix.resize( mat_size, mat_size, false );

    noalias( rDampingMatrix ) = ZeroMatrix( mat_size, mat_size );

    // 1.-Calculate StiffnessMatrix:

    MatrixType StiffnessMatrix  = Matrix();
    VectorType ResidualVector  = Vector();

    this->CalculateAll(StiffnessMatrix, ResidualVector, rCurrentProcessInfo, true, false);

    // 2.-Calculate MassMatrix:

    MatrixType MassMatrix  = Matrix();

    this->CalculateMassMatrix ( MassMatrix, rCurrentProcessInfo );

    // 3.-Get Damping Coeffitients (RAYLEIGH_ALPHA, RAYLEIGH_BETA)
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

    // 4.-Compose the Damping Matrix:

    // Rayleigh Damping Matrix: alpha*M + beta*K
    noalias( rDampingMatrix ) += alpha * MassMatrix;
    noalias( rDampingMatrix ) += beta  * StiffnessMatrix;

    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::CalculateOnIntegrationPoints(
    const Variable<bool>& rVariable,
    std::vector<bool>& rOutput,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    const GeometryType::IntegrationPointsArrayType &integration_points = GetGeometry().IntegrationPoints();

    if (rOutput.size() != GetGeometry().IntegrationPoints(  ).size())
        rOutput.resize(GetGeometry().IntegrationPoints(  ).size());

    for (unsigned int point_number = 0; point_number < integration_points.size(); ++point_number){
        bool flag = false;
        mConstitutiveLawVector[point_number]->GetValue(rVariable, flag );
        rOutput[point_number] = flag;
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
    const GeometryType::IntegrationMethod integration_method =
        GetGeometry().GetDefaultIntegrationMethod();
    const GeometryType::IntegrationPointsArrayType &integration_points = GetGeometry().IntegrationPoints(integration_method);

    if ( rOutput.size() != integration_points.size() )
        rOutput.resize( integration_points.size() );

    if (rVariable == INTEGRATION_WEIGHT) {
        const unsigned int number_of_nodes = GetGeometry().size();
        const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
        const unsigned int strain_size = mConstitutiveLawVector[0]->GetStrainSize();

        KinematicVariables this_kinematic_variables(strain_size, dimension, number_of_nodes);

        for (unsigned int point_number = 0; point_number < integration_points.size(); ++point_number) {
            this_kinematic_variables.detJ0 = CalculateDerivativesOnReferenceConfiguration(this_kinematic_variables.J0,
                                                                                this_kinematic_variables.InvJ0,
                                                                                this_kinematic_variables.DN_DX,
                                                                                point_number,
                                                                                integration_method);

            double integration_weight = GetIntegrationWeight(integration_points,
                                                                point_number,
                                                                this_kinematic_variables.detJ0);

            if (dimension == 2 && this->GetProperties().Has(THICKNESS))
                integration_weight *= this->GetProperties()[THICKNESS];

            rOutput[point_number] = integration_weight;
        }
    } else if ( rVariable == STRAIN_ENERGY ) {
        const unsigned int number_of_nodes = GetGeometry().size();
        const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
        const unsigned int strain_size = mConstitutiveLawVector[0]->GetStrainSize();

        KinematicVariables this_kinematic_variables(strain_size, dimension, number_of_nodes);
        ConstitutiveVariables this_constitutive_variables(strain_size);

        // Create constitutive law parameters:
        ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

        // Reading integration points
        const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints(  );

        // If strain has to be computed inside of the constitutive law with PK2
        Values.SetStrainVector(this_constitutive_variables.StrainVector); //this is the input  parameter

        for (unsigned int point_number = 0; point_number < integration_points.size(); ++point_number) {
            // Compute element kinematics B, F, DN_DX ...
            CalculateKinematicVariables(this_kinematic_variables, point_number, integration_method);

            // Compute material reponse
            CalculateConstitutiveVariables(this_kinematic_variables, this_constitutive_variables, Values, point_number, integration_points, GetStressMeasure());

            double StrainEnergy = 0.0;

            mConstitutiveLawVector[point_number]->CalculateValue(Values, STRAIN_ENERGY, StrainEnergy);

            rOutput[point_number] = StrainEnergy;  // 1/2 * sigma * epsilon
        }
    } else if (rVariable == VON_MISES_STRESS) {
        const unsigned int number_of_nodes = GetGeometry().size();
        const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
        const unsigned int strain_size = mConstitutiveLawVector[0]->GetStrainSize();

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

        for (unsigned int point_number = 0; point_number < integration_points.size(); ++point_number) {
            // Compute element kinematics B, F, DN_DX ...
            CalculateKinematicVariables(this_kinematic_variables, point_number, integration_method);

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
        for ( unsigned int point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number )
            rOutput[point_number] = mConstitutiveLawVector[point_number]->GetValue( rVariable, rOutput[point_number] );
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
    const GeometryType::IntegrationPointsArrayType &integration_points = GetGeometry().IntegrationPoints();

    if ( rOutput.size() != GetGeometry().IntegrationPoints(  ).size() )
        rOutput.resize( GetGeometry().IntegrationPoints(  ).size() );

    if (rVariable == INTEGRATION_COORDINATES) {
        const unsigned int number_of_nodes = GetGeometry().size();
        const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
        const unsigned int strain_size = mConstitutiveLawVector[0]->GetStrainSize();

        KinematicVariables this_kinematic_variables(strain_size, dimension, number_of_nodes);

        for (unsigned int point_number = 0; point_number < integration_points.size(); ++point_number)
        {
            Point global_point;
            GetGeometry().GlobalCoordinates(global_point, integration_points[point_number]);

            rOutput[point_number] = global_point.Coordinates();
        }
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
    const GeometryType::IntegrationMethod integration_method =
        GetGeometry().GetDefaultIntegrationMethod();
    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( integration_method );
    if ( rOutput.size() != integration_points.size() )
        rOutput.resize( integration_points.size() );

    if ( rVariable == INSITU_STRESS ) {
        const unsigned int strain_size = mConstitutiveLawVector[0]->GetStrainSize();
        Vector strain_vector( strain_size );

        for ( unsigned int point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number ) {
            if ( rOutput[point_number].size() != strain_vector.size() )
                rOutput[point_number].resize( strain_vector.size(), false );

            rOutput[point_number] = mConstitutiveLawVector[point_number]->GetValue( INSITU_STRESS, rOutput[point_number] );
        }
    } else if ( rVariable == CAUCHY_STRESS_VECTOR || rVariable == PK2_STRESS_VECTOR ) {
        // Create and initialize element variables:
        const unsigned int number_of_nodes = GetGeometry().size();
        const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
        const unsigned int strain_size = mConstitutiveLawVector[0]->GetStrainSize();

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
        for ( unsigned int point_number = 0; point_number < integration_points.size(); ++point_number ) {
            // Compute element kinematics B, F, DN_DX ...
            CalculateKinematicVariables(this_kinematic_variables, point_number, integration_method);

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
        const unsigned int number_of_nodes = GetGeometry().size();
        const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
        const unsigned int strain_size = mConstitutiveLawVector[0]->GetStrainSize();

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

        //reading integration points
        for ( unsigned int point_number = 0; point_number < integration_points.size(); ++point_number ) {
            // Compute element kinematics B, F, DN_DX ...
            CalculateKinematicVariables(this_kinematic_variables, point_number, integration_method);

            // Compute material reponse
            CalculateConstitutiveVariables(this_kinematic_variables, this_constitutive_variables, Values, point_number, integration_points, GetStressMeasure());

            if ( rOutput[point_number].size() != strain_size)
                rOutput[point_number].resize( strain_size, false );

            rOutput[point_number] = this_constitutive_variables.StrainVector;
        }
    } else {
        for ( unsigned int ii = 0; ii < mConstitutiveLawVector.size(); ++ii )
            rOutput[ii] = mConstitutiveLawVector[ii]->GetValue( rVariable, rOutput[ii] );
    }
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::CalculateOnIntegrationPoints(
    const Variable<Matrix >& rVariable,
    std::vector< Matrix >& rOutput,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    const GeometryType::IntegrationMethod integration_method =
        GetGeometry().GetDefaultIntegrationMethod();
    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( integration_method );
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    if ( rOutput.size() != integration_points.size() )
        rOutput.resize( integration_points.size() );

    if ( rVariable == CAUCHY_STRESS_TENSOR || rVariable == PK2_STRESS_TENSOR ) {
        std::vector<Vector> stress_vector;

        if( rVariable == CAUCHY_STRESS_TENSOR )
            this->CalculateOnIntegrationPoints( CAUCHY_STRESS_VECTOR, stress_vector, rCurrentProcessInfo );
        else
            this->CalculateOnIntegrationPoints( PK2_STRESS_VECTOR, stress_vector, rCurrentProcessInfo );

        // Loop integration points
        for ( unsigned int point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number ) {
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
        for ( unsigned int point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number ) {
            if ( rOutput[point_number].size2() != dimension )
                rOutput[point_number].resize( dimension, dimension, false );

            rOutput[point_number] = MathUtils<double>::StrainVectorToTensor(strain_vector[point_number]);
        }
    } else if ( rVariable == CONSTITUTIVE_MATRIX ) {
        // Create and initialize element variables:
        const unsigned int number_of_nodes = GetGeometry().size();
        const unsigned int strain_size = mConstitutiveLawVector[0]->GetStrainSize();

        KinematicVariables this_kinematic_variables(strain_size, dimension, number_of_nodes);
        ConstitutiveVariables this_constitutive_variables(strain_size);

        // Create constitutive law parameters:
        ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

        // Set constitutive law flags:
        Flags& ConstitutiveLawOptions=Values.GetOptions();
        ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, UseElementProvidedStrain());
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, false);
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

        Values.SetConstitutiveMatrix(this_constitutive_variables.D); //this is the output parameter

        // Reading integration points
        for ( unsigned int point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number ) {
            // Compute element kinematics B, F, DN_DX ...
            CalculateKinematicVariables(this_kinematic_variables, point_number, integration_method);

            // Compute material reponse
            CalculateConstitutiveVariables(this_kinematic_variables, this_constitutive_variables, Values, point_number, integration_points, GetStressMeasure());

            // Call the constitutive law to update material variables
            mConstitutiveLawVector[point_number]->CalculateMaterialResponse(Values, GetStressMeasure());

            if( rOutput[point_number].size2() != this_constitutive_variables.D.size2() )
                rOutput[point_number].resize( this_constitutive_variables.D.size1() , this_constitutive_variables.D.size2() , false );

            rOutput[point_number] = this_constitutive_variables.D;
        }
    } else if ( rVariable == DEFORMATION_GRADIENT ) { // VARIABLE SET FOR TRANSFER PURPOUSES
        // Create and initialize element variables:
        const unsigned int number_of_nodes = GetGeometry().size();
        const unsigned int strain_size = mConstitutiveLawVector[0]->GetStrainSize();

        KinematicVariables this_kinematic_variables(strain_size, dimension, number_of_nodes);
        ConstitutiveVariables this_constitutive_variables(strain_size);

        // Create constitutive law parameters:
        ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

        // Reading integration points
        for ( unsigned int point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number ) {
            // Compute element kinematics B, F, DN_DX ...
            CalculateKinematicVariables(this_kinematic_variables, point_number, integration_method);

            if( rOutput[point_number].size2() != this_kinematic_variables.F.size2() )
                rOutput[point_number].resize( this_kinematic_variables.F.size1() , this_kinematic_variables.F.size2() , false );

            rOutput[point_number] = this_kinematic_variables.F;
        }
    } else {
        for ( unsigned int ii = 0; ii < mConstitutiveLawVector.size(); ++ii )
            rOutput[ii] = mConstitutiveLawVector[ii]->GetValue( rVariable , rOutput[ii] );
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
    for ( unsigned int point_number = 0; point_number < GetGeometry().IntegrationPoints(  ).size(); ++point_number ) {
        mConstitutiveLawVector[point_number]->SetValue( rVariable,
                                                        rValues[point_number],
                                                        rCurrentProcessInfo
                                                        );
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

    for ( unsigned int point_number = 0; point_number < GetGeometry().IntegrationPoints(  ).size(); ++point_number ) {
        mConstitutiveLawVector[point_number]->SetValue( rVariable,
                                                        rValues[point_number],
                                                        rCurrentProcessInfo
                                                        );
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
        const std::size_t integration_points_number = GetGeometry().IntegrationPoints(  ).size();
        for ( std::size_t i_gp = 0; i_gp < integration_points_number; ++i_gp ) {
            mConstitutiveLawVector[i_gp] = rValues[i_gp];
        }
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
    for ( unsigned int point_number = 0; point_number < GetGeometry().IntegrationPoints(  ).size(); ++point_number ) {
        mConstitutiveLawVector[point_number]->SetValue( rVariable,
                                                        rValues[point_number],
                                                        rCurrentProcessInfo
                                                        );
    }

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
    const Variable<Vector>& rVariable,
    std::vector<Vector>& rValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    const unsigned int size = GetGeometry().IntegrationPoints(  ).size();

    if ( rValues.size() != size )
        rValues.resize( size );

    if ( rVariable == STRESSES ) {
        for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); ++i ) {
            if ( rValues[i].size() != 6 )
                rValues[i].resize( 6, false );
            noalias( rValues[i] ) = mConstitutiveLawVector[i]->GetValue( STRESSES, rValues[i] );
        }
    } else if ( rVariable == MATERIAL_PARAMETERS ) {
        for ( unsigned int point_number = 0; point_number < GetGeometry().IntegrationPoints(  ).size(); ++point_number )
            rValues[point_number] = mConstitutiveLawVector[point_number]->GetValue( MATERIAL_PARAMETERS, rValues[point_number] );
    } else if ( rVariable == INTERNAL_VARIABLES ) {
        for ( unsigned int point_number = 0; point_number < GetGeometry().IntegrationPoints(  ).size(); ++point_number )
            rValues[point_number] = mConstitutiveLawVector[point_number]->GetValue( INTERNAL_VARIABLES, rValues[point_number] );
    } else {
        CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
    }
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
        const std::size_t integration_points_number = mConstitutiveLawVector.size();
        if (rValues.size() != mConstitutiveLawVector.size())
            rValues.resize(integration_points_number);
        for (std::size_t i_gp = 0; i_gp < integration_points_number; ++i_gp)
            rValues[i_gp] = mConstitutiveLawVector[i_gp];
    }
}

/***********************************************************************************/
/***********************************************************************************/

int  BaseSolidElement::Check( const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY;

    const unsigned int number_of_nodes = this->GetGeometry().size();
    const unsigned int dimension = this->GetGeometry().WorkingSpaceDimension();

    // Verify that the variables are correctly initialized
    KRATOS_CHECK_VARIABLE_KEY(DISPLACEMENT)
    KRATOS_CHECK_VARIABLE_KEY(VELOCITY)
    KRATOS_CHECK_VARIABLE_KEY(ACCELERATION)
    KRATOS_CHECK_VARIABLE_KEY(DENSITY)
    KRATOS_CHECK_VARIABLE_KEY(VOLUME_ACCELERATION)
    KRATOS_CHECK_VARIABLE_KEY(THICKNESS)

    // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
    for ( unsigned int i = 0; i < number_of_nodes; i++ ) {
        Node<3> &rnode = this->GetGeometry()[i];
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DISPLACEMENT,rnode)
//         KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VELOCITY,rnode)
//         KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(ACCELERATION,rnode)

        KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_X, rnode)
        KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_Y, rnode)
        KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_Z, rnode)
    }

    // Verify that the constitutive law exists
    KRATOS_ERROR_IF_NOT(this->GetProperties().Has( CONSTITUTIVE_LAW )) << "Constitutive law not provided for property " << this->GetProperties().Id() << std::endl;

    // Verify that the constitutive law has the correct dimension
    const unsigned int strain_size = this->GetProperties().GetValue( CONSTITUTIVE_LAW )->GetStrainSize();
    if ( dimension == 2 ) {
        KRATOS_ERROR_IF( strain_size < 3 || strain_size > 4) << "Wrong constitutive law used. This is a 2D element! expected strain size is 3 or 4 (el id = ) " << this->Id() << std::endl;
    } else {
        KRATOS_ERROR_IF_NOT(strain_size == 6) << "Wrong constitutive law used. This is a 3D element! expected strain size is 6 (el id = ) "<<  this->Id() << std::endl;
    }

    // Check constitutive law
    if ( mConstitutiveLawVector.size() > 0 ) {
        return mConstitutiveLawVector[0]->Check( GetProperties(), GetGeometry(), rCurrentProcessInfo );
    }

    return 0;

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::CalculateAll(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo,
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
    const unsigned int point_number,
    const double detJ
    )
{
    return rThisIntegrationPoints[point_number].Weight() * detJ;
}

void BaseSolidElement::CalculateShapeGradientOfMassMatrix(MatrixType& rMassMatrix, ShapeParameter Deriv)
{
    KRATOS_TRY;

    const auto& r_geom = GetGeometry();
    const auto& r_prop = GetProperties();
    unsigned dim = r_geom.WorkingSpaceDimension();
    rMassMatrix = ZeroMatrix(dim * r_geom.size(), dim * r_geom.size());

    KRATOS_ERROR_IF_NOT(r_prop.Has(DENSITY))
        << "DENSITY has to be provided for the calculation of the MassMatrix!"
        << std::endl;

    const double density = r_prop[DENSITY];
    const double thickness =
        (dim == 2 && r_prop.Has(THICKNESS)) ? r_prop[THICKNESS] : 1.0;

    const IntegrationMethod integration_method =
        IntegrationUtilities::GetIntegrationMethodForExactMassMatrixEvaluation(r_geom);
    const Matrix& Ncontainer = r_geom.ShapeFunctionsValues(integration_method);
    Matrix J0(dim, dim), DN_DX0_deriv;
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
            const unsigned index_i = i * dim;

            for (unsigned j = 0; j < r_geom.size(); ++j)
            {
                const unsigned index_j = j * dim;
                const double NiNj_weight = rN[i] * rN[j] * integration_weight * density;

                for (unsigned k = 0; k < dim; ++k)
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
    const unsigned int PointNumber,
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
    const unsigned int PointNumber,
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints,
    const ConstitutiveLaw::StressMeasure ThisStressMeasure
    )
{
    // Here we essentially set the input parameters
    rValues.SetDeterminantF(rThisKinematicVariables.detF); // Assuming the determinant is computed somewhere else
    rValues.SetDeformationGradientF(rThisKinematicVariables.F); //F computed somewhere else

    // Here we set the space on which the results shall be written
    rValues.SetConstitutiveMatrix(rThisConstitutiveVariables.D); // Assuming the determinant is computed somewhere else
    rValues.SetStressVector(rThisConstitutiveVariables.StressVector); //F computed somewhere else

    // Actually do the computations in the ConstitutiveLaw
    mConstitutiveLawVector[PointNumber]->CalculateMaterialResponse(rValues, ThisStressMeasure); //here the calculations are actually done
}

/***********************************************************************************/
/***********************************************************************************/

Matrix& BaseSolidElement::CalculateDeltaDisplacement(Matrix& DeltaDisplacement)
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    DeltaDisplacement.resize(number_of_nodes , dimension);

    for ( unsigned int i_node = 0; i_node < number_of_nodes; i_node++ ) {
        const array_1d<double, 3 >& current_displacement  = GetGeometry()[i_node].FastGetSolutionStepValue(DISPLACEMENT);
        const array_1d<double, 3 >& previous_displacement = GetGeometry()[i_node].FastGetSolutionStepValue(DISPLACEMENT,1);

        for ( unsigned int j_dim = 0; j_dim < dimension; ++j_dim )
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
    const unsigned int PointNumber,
    IntegrationMethod ThisIntegrationMethod
    )
{
    GeometryType& r_geom = GetGeometry();
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
    const unsigned int PointNumber,
    IntegrationMethod ThisIntegrationMethod
    )
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

Vector BaseSolidElement::GetBodyForce(
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints,
    const unsigned int PointNumber
    ) const
{
    Vector body_force = ZeroVector(3);

    double density = 0.0;
    if (GetProperties().Has( DENSITY ) == true)
        density = GetProperties()[DENSITY];

    if (GetProperties().Has( VOLUME_ACCELERATION ) == true)
        body_force += density * GetProperties()[VOLUME_ACCELERATION];

    if( GetGeometry()[0].SolutionStepsDataHas(VOLUME_ACCELERATION) ) {
        Vector N;
        N = GetGeometry().ShapeFunctionsValues(N, IntegrationPoints[PointNumber].Coordinates());
        for (unsigned int i_node = 0; i_node < this->GetGeometry().size(); ++i_node)
            noalias(body_force) += N[i_node] * density * GetGeometry()[i_node].FastGetSolutionStepValue(VOLUME_ACCELERATION);
    }

    return body_force;
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::CalculateAndAddKm(
    MatrixType& rLeftHandSideMatrix,
    const Matrix& B,
    const Matrix& D,
    const double IntegrationWeight
    )
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
    )
{
    KRATOS_TRY

    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
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
    const Vector& rBodyForce,
    const Vector& rStressVector,
    const double IntegrationWeight
    )
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
    const Vector& rBodyForce,
    VectorType& rRightHandSideVector,
    const double Weight
    )
{
    KRATOS_TRY;

    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    for ( unsigned int i = 0; i < number_of_nodes; ++i ) {
        const unsigned int index = dimension * i;

        for ( unsigned int j = 0; j < dimension; ++j )
            rRightHandSideVector[index + j] += Weight * rN[i] * rBodyForce[j];
    }

    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Element );
    rSerializer.save("mConstitutiveLawVector", mConstitutiveLawVector);
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Element );
    rSerializer.load("mConstitutiveLawVector", mConstitutiveLawVector);
}
} // Namespace Kratos


