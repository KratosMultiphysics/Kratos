//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:        BSD License
//                  Kratos default license: kratos/license.txt
//
//  Main authors:    Andi Makarim Katili
//


// System includes
#include <omp.h>
#include <sstream>

// External includes

// Project includes
#include "includes/define.h"
#include "custom_elements/mpm_soft_stiffness.hpp"
#include "utilities/math_utils.h"
#include "includes/constitutive_law.h"
#include "mpm_application_variables.h"
#include "includes/checks.h"
#include "custom_utilities/mpm_energy_calculation_utility.h"
#include "custom_utilities/mpm_explicit_utilities.h"
#include "custom_utilities/mpm_math_utilities.h"

namespace Kratos
{

/**
 * Flags related to the element computation
 */
KRATOS_CREATE_LOCAL_FLAG( MPMSoftStiffness, COMPUTE_RHS_VECTOR,                 0 );
KRATOS_CREATE_LOCAL_FLAG( MPMSoftStiffness, COMPUTE_LHS_MATRIX,                 1 );
KRATOS_CREATE_LOCAL_FLAG( MPMSoftStiffness, COMPUTE_RHS_VECTOR_WITH_COMPONENTS, 2 );
KRATOS_CREATE_LOCAL_FLAG( MPMSoftStiffness, COMPUTE_LHS_MATRIX_WITH_COMPONENTS, 3 );

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

MPMSoftStiffness::MPMSoftStiffness( )
    : Element( )
    , mGridVariables()
{
    //DO NOT CALL IT: only needed for Register and Serialization!!!
}
//******************************CONSTRUCTOR*******************************************
//************************************************************************************
MPMSoftStiffness::MPMSoftStiffness( IndexType NewId, GeometryType::Pointer pGeometry )
    : Element( NewId, pGeometry )
    , mGridVariables()
{
    //DO NOT ADD DOFS HERE!!!
}

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

MPMSoftStiffness::MPMSoftStiffness( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
    : Element( NewId, pGeometry, pProperties )
    , mGridVariables()
{
    mFinalizedStep = true;


}
//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

MPMSoftStiffness::MPMSoftStiffness( MPMSoftStiffness const& rOther)
    :Element(rOther)
    ,mGridVariables(rOther.mGridVariables)
    ,mDeformationGradientF0(rOther.mDeformationGradientF0)
    ,mDeterminantF0(rOther.mDeterminantF0)
    ,mConstitutiveLawVector(rOther.mConstitutiveLawVector)
    ,mFinalizedStep(rOther.mFinalizedStep)
{
}

//******************************ASSIGNMENT OPERATOR***********************************
//************************************************************************************

MPMSoftStiffness&  MPMSoftStiffness::operator=(MPMSoftStiffness const& rOther)
{
    Element::operator=(rOther);

    mGridVariables = rOther.mGridVariables;

    mDeformationGradientF0.clear();
    mDeformationGradientF0 = rOther.mDeformationGradientF0;

    mDeterminantF0 = rOther.mDeterminantF0;
    mConstitutiveLawVector = rOther.mConstitutiveLawVector;

    return *this;
}

//*********************************OPERATIONS*****************************************
//************************************************************************************

Element::Pointer MPMSoftStiffness::Create( IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties ) const
{
    return Element::Pointer( new MPMSoftStiffness( NewId, GetGeometry().Create( ThisNodes ), pProperties ) );
}

Element::Pointer MPMSoftStiffness::Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive< MPMSoftStiffness >(NewId, pGeom, pProperties);
}

//************************************CLONE*******************************************
//************************************************************************************

Element::Pointer MPMSoftStiffness::Clone( IndexType NewId, NodesArrayType const& rThisNodes ) const
{
    MPMSoftStiffness NewElement (NewId, GetGeometry().Create( rThisNodes ), pGetProperties() );

    NewElement.mGridVariables = mGridVariables;

    NewElement.SetConstitutiveLawVector(mConstitutiveLawVector);    

    NewElement.mDeformationGradientF0 = mDeformationGradientF0;

    NewElement.mDeterminantF0 = mDeterminantF0;

    return Element::Pointer( new MPMSoftStiffness(NewElement) );
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************
MPMSoftStiffness::~MPMSoftStiffness()
{
}


//************************************************************************************
//************************************************************************************

void MPMSoftStiffness::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Initialization should not be done again in a restart!
    if (!rCurrentProcessInfo[IS_RESTARTED]) {
        // Initialize parameters
        const unsigned int dimension = this->GetGeometry().WorkingSpaceDimension();
        mDeterminantF0 = 1;
        mDeformationGradientF0 = IdentityMatrix(dimension);



        // Initialize grid volume
        mGridVariables.volume = this->GetGeometry().DomainSize();
        if (GetGeometry().WorkingSpaceDimension() == 2 && this->GetProperties().Has(THICKNESS)) {
            mGridVariables.volume *= this->GetProperties()[THICKNESS];
        }

        // Initialize penalty factor
        mGridVariables.penalty_factor = this->GetProperties()[PENALTY_FACTOR]; // TODO: no error thrown when penalty factor value is not set 

        const auto& gp_size = this->GetGeometry().IntegrationPointsNumber();
        //Constitutive Law initialisation
        if ( mConstitutiveLawVector.size() != gp_size )
            mConstitutiveLawVector.resize(gp_size);
        // Initialize constitutive law and materials
        InitializeMaterial(rCurrentProcessInfo);
    }

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void MPMSoftStiffness::InitializeGeneralVariables (GeneralVariables& rVariables, const ProcessInfo& rCurrentProcessInfo)
{
    const SizeType number_of_nodes = GetGeometry().size();
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();
    const SizeType strain_size = GetProperties().GetValue(CONSTITUTIVE_LAW)->GetStrainSize();
    const bool is_axisymmetric = (rCurrentProcessInfo.Has(IS_AXISYMMETRIC))
        ? rCurrentProcessInfo.GetValue(IS_AXISYMMETRIC)
        : false;
    const SizeType def_grad_dim = (is_axisymmetric)
        ? 3
        : dimension;

    rVariables.detFT = 1;

    rVariables.B.resize(strain_size, number_of_nodes * dimension, false );

    // soft stiffness is always calculated from undeformed configuration
    rVariables.FT.resize(def_grad_dim, def_grad_dim, false );
    rVariables.FT = IdentityMatrix(def_grad_dim);

    rVariables.ConstitutiveMatrix.resize(strain_size, strain_size, false );

    rVariables.StrainVector.resize(strain_size, false );

    rVariables.StressVector.resize(strain_size, false );

    rVariables.N.resize( number_of_nodes, false );
    rVariables.DN_DX.resize( number_of_nodes, dimension, false );

}
//************************************************************************************
//************************************************************************************

void MPMSoftStiffness::SetGeneralVariables(GeneralVariables& rVariables,
                                           ConstitutiveLaw::Parameters& rValues,
                                           const IndexType GaussPointNumber)
{
    rValues.SetDeterminantF(rVariables.detFT); // detFT always = 1
    rValues.SetDeformationGradientF(rVariables.FT); // FT always = identity Matrix
    rValues.SetStrainVector(rVariables.StrainVector);
    rValues.SetStressVector(rVariables.StressVector);
    rValues.SetConstitutiveMatrix(rVariables.ConstitutiveMatrix);
    rValues.SetShapeFunctionsValues(rVariables.N);
    rValues.SetShapeFunctionsDerivatives(rVariables.DN_DX);
}

//************************************************************************************
//************************************************************************************

void MPMSoftStiffness::CalculateElementalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo,
    const bool CalculateStiffnessMatrixFlag,
    const bool CalculateResidualVectorFlag)
{
    KRATOS_TRY
    GeometryType& rGeometry = this->GetGeometry(); // Should return background grid geometry

    // Create and initialize element variables:
    GeneralVariables Variables;
    this->InitializeGeneralVariables(Variables,rCurrentProcessInfo);

    // Create constitutive law parameters:
    ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);
    Flags &ConstitutiveLawOptions=Values.GetOptions();
    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
    Variables.StressMeasure = ConstitutiveLaw::StressMeasure_Cauchy;
    
    // Note: soft stiffness geometry is the background grid geometry, not gauss quadrature
    const GeometryType::IntegrationMethod& rThisIntegrationMethod = rGeometry.GetDefaultIntegrationMethod();
    // KRATOS_INFO("MPMSoftStiffness") << "Integration method used: " << rThisIntegrationMethod << std::endl; // TODO: Remove this

    // Get geometry's default shape functions, it's derivatives and Jacobian
    VectorType DeterminantsOfJacobian;
    Variables.ShapeFunctions = rGeometry.ShapeFunctionsValues();
    rGeometry.ShapeFunctionsIntegrationPointsGradients(Variables.ShapeFunctionsGradients, DeterminantsOfJacobian, rThisIntegrationMethod);
    
    // Computing in all integrations points  
    const GeometryType::IntegrationPointsArrayType& rThisIntegrationPoints = rGeometry.IntegrationPoints(rThisIntegrationMethod);
    for ( IndexType gp_number = 0; gp_number < rThisIntegrationPoints.size(); ++gp_number ) {
        // Set all components to zero
        Variables.N.clear();
        Variables.DN_DX.clear();
        Variables.B.clear();

        // Compute element kinematics B
        this->CalculateKinematics(Variables, gp_number, rCurrentProcessInfo);
        
        // Set general variables to constitutivelaw parameters of current GP
        this->SetGeneralVariables(Variables, Values, gp_number);

        // Compute material response
        mConstitutiveLawVector[gp_number]->CalculateMaterialResponse(Values, Variables.StressMeasure);

        // Calculating weights for integration on the reference configuration
        double integration_weight = CalculateIntegrationWeight(rThisIntegrationPoints[gp_number].Weight(), DeterminantsOfJacobian[gp_number]);
        
        // Contributions to stiffness matrix calculated on the initial reference config
        if ( CalculateStiffnessMatrixFlag ) {
            this->CalculateAndAddLHS(rLeftHandSideMatrix, Variables, integration_weight, rCurrentProcessInfo);
        }
        // 
        if ( CalculateResidualVectorFlag ) {
            this->CalculateAndAddRHS(rRightHandSideVector, Variables, integration_weight, rCurrentProcessInfo);
        }
    }
    rLeftHandSideMatrix *= mGridVariables.stiffness_multiplier;
    




    KRATOS_CATCH( "" )
}
//*********************************COMPUTE KINEMATICS*********************************
//************************************************************************************


void MPMSoftStiffness::CalculateKinematics(GeneralVariables& rVariables, const IndexType GaussPointNumber, const ProcessInfo& rCurrentProcessInfo)

{
    KRATOS_TRY

    // axissymmteric is not considered yet
    // const bool is_axisymmetric = false;
    // const bool is_axisymmetric = (rCurrentProcessInfo.Has(IS_AXISYMMETRIC))
    //     ? rCurrentProcessInfo.GetValue(IS_AXISYMMETRIC)
    //     : false;
    // if (is_axisymmetric) {
    //     rVariables.CurrentRadius = MPMMathUtilities<double>::CalculateRadius(r_N, GetGeometry());
    //     rVariables.ReferenceRadius = MPMMathUtilities<double>::CalculateRadius(r_N, GetGeometry(), Initial);
    // }

    // no need for CalculateDeformationGradient. F is just Identity matrix (undeformed configuration)
    // Calculating the reference jacobian from cartesian coordinates to parent coordinates for the MP element [dx_n/d£]
    Matrix Jacobian;
    GetGeometry().Jacobian(Jacobian, GaussPointNumber);

    // Calculating the inverse of the jacobian and the parameters needed [d£/dx_n]
    Matrix InvJ;
    double detJ;
    MathUtils<double>::InvertMatrix( Jacobian, InvJ, detJ);

    rVariables.N = row(rVariables.ShapeFunctions, GaussPointNumber);
    rVariables.DN_DX = prod(rVariables.ShapeFunctionsGradients[GaussPointNumber], InvJ);
    // Compute the deformation matrix B
    this->CalculateDeformationMatrix(rVariables.B, rVariables.DN_DX);

    KRATOS_CATCH( "" )
}
//************************************************************************************

void MPMSoftStiffness::CalculateDeformationMatrix(Matrix& rB,
        const Matrix& rDN_DX, const bool IsAxisymmetric)
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

    // if (IsAxisymmetric)
    // {
    //     const double radius = MPMMathUtilities<double>::CalculateRadius(rN, GetGeometry());

    //     for (unsigned int i = 0; i < number_of_nodes; i++)
    //     {
    //         const unsigned int index = dimension * i;

    //         rB(0, index + 0) = rDN_DX(i, 0);
    //         rB(1, index + 1) = rDN_DX(i, 1);
    //         rB(2, index + 0) = rN(0, i) / radius;
    //         rB(3, index + 0) = rDN_DX(i, 1);
    //         rB(3, index + 1) = rDN_DX(i, 0);
    //     }
    // }
    if( dimension == 2 )
    {
        for ( unsigned int i = 0; i < number_of_nodes; i++ )
        {
            unsigned int index = 2 * i;
            rB( 0, index + 0 ) = rDN_DX( i, 0 );
            rB( 1, index + 1 ) = rDN_DX( i, 1 );
            rB( 2, index + 0 ) = rDN_DX( i, 1 );
            rB( 2, index + 1 ) = rDN_DX( i, 0 );
        }
    }
    else if( dimension == 3 )
    {
        for ( unsigned int i = 0; i < number_of_nodes; i++ )
        {
            unsigned int index = 3 * i;

            rB( 0, index + 0 ) = rDN_DX( i, 0 );
            rB( 1, index + 1 ) = rDN_DX( i, 1 );
            rB( 2, index + 2 ) = rDN_DX( i, 2 );

            rB( 3, index + 0 ) = rDN_DX( i, 1 );
            rB( 3, index + 1 ) = rDN_DX( i, 0 );

            rB( 4, index + 1 ) = rDN_DX( i, 2 );
            rB( 4, index + 2 ) = rDN_DX( i, 1 );

            rB( 5, index + 0 ) = rDN_DX( i, 2 );
            rB( 5, index + 2 ) = rDN_DX( i, 0 );
        }
    }
    else
    {
        KRATOS_ERROR <<  "Dimension given is wrong!" << std::endl;
    }

    KRATOS_CATCH( "" )
}
//************************************************************************************
//************************************************************************************

void MPMSoftStiffness::CalculateAndAddRHS(
    VectorType& rRightHandSideVector,
    GeneralVariables& rVariables,
    const double& rIntegrationWeight,
    const ProcessInfo& rCurrentProcessInfo)
{
    // Operation performed: rRightHandSideVector -= rLeftHandSideMatrix*CurrentDisp
    this->CalculateAndAddInternalForces(rRightHandSideVector, rVariables, rIntegrationWeight, rCurrentProcessInfo);
}

//************************************************************************************
//************************************************************************************

void MPMSoftStiffness::CalculateAndAddInternalForces(
    VectorType& rRightHandSideVector,
    GeneralVariables& rVariables,
    const double& rIntegrationWeight,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    VectorType internal_forces = rIntegrationWeight * prod( trans( rVariables.B ), rVariables.StressVector );
    noalias( rRightHandSideVector ) -= internal_forces;

    KRATOS_CATCH( "" )
}
//************************************************************************************
//************************************************************************************

void MPMSoftStiffness::CalculateAndAddLHS(
    MatrixType& rLeftHandSideMatrix,
    GeneralVariables& rVariables,
    const double& rIntegrationWeight,
    const ProcessInfo& rCurrentProcessInfo)
{
    // Operation performed: add K_material to the rLefsHandSideMatrix 
    this->CalculateAndAddKuum( rLeftHandSideMatrix, rVariables, rIntegrationWeight );

    // Geometric stiffness is not considered for soft stiffness 
}
//************************************************************************************
//************************************************************************************

void MPMSoftStiffness::CalculateAndAddKuum(
    MatrixType& rLeftHandSideMatrix,
    GeneralVariables& rVariables,
    const double& rIntegrationWeight)
{
    KRATOS_TRY

    noalias( rLeftHandSideMatrix ) += prod( trans( rVariables.B ),  rIntegrationWeight * Matrix( prod( rVariables.ConstitutiveMatrix, rVariables.B ) ) );

    KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************
void MPMSoftStiffness::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    MatrixType left_hand_side_matrix = Matrix(0, 0);

    const SizeType mat_size = GetNumberOfDofs() * GetGeometry().size();
    if (rRightHandSideVector.size() != mat_size) {
        rRightHandSideVector.resize(mat_size, false);
    }
    rRightHandSideVector = ZeroVector(mat_size);

    CalculateElementalSystem(left_hand_side_matrix, rRightHandSideVector,
        rCurrentProcessInfo, false, true);
}

//************************************************************************************
//************************************************************************************


void MPMSoftStiffness::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    VectorType right_hand_side_vector = Vector(0);

    const SizeType mat_size = GetNumberOfDofs() * GetGeometry().size();
    if (rLeftHandSideMatrix.size1() != mat_size && rLeftHandSideMatrix.size2() != mat_size) {
        rLeftHandSideMatrix.resize(mat_size, mat_size, false);
    }
    noalias(rLeftHandSideMatrix) = ZeroMatrix(mat_size, mat_size);

    CalculateElementalSystem(
        rLeftHandSideMatrix, right_hand_side_vector,
        rCurrentProcessInfo, true, false);
}
//************************************************************************************
//************************************************************************************


void MPMSoftStiffness::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    const SizeType mat_size = GetNumberOfDofs() * GetGeometry().size();
    if (rLeftHandSideMatrix.size1() != mat_size && rLeftHandSideMatrix.size2() != mat_size) {
        rLeftHandSideMatrix.resize(mat_size, mat_size, false);
    }
    noalias(rLeftHandSideMatrix) = ZeroMatrix(mat_size, mat_size);

    if (rRightHandSideVector.size() != mat_size) {
        rRightHandSideVector.resize(mat_size, false);
    }
    rRightHandSideVector = ZeroVector(mat_size);

    CalculateElementalSystem(
        rLeftHandSideMatrix, rRightHandSideVector,
        rCurrentProcessInfo, true, true);
}

//*******************************************************************************************
//*******************************************************************************************
void MPMSoftStiffness::InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo )
{
    
    /* NOTE:
    In the InitializeSolutionStep of each time step the nodal initial Elements are evaluated.
    This function is called by the base scheme class.*/
    
    // Get the sum of MP_Volume inside parent geometry
    KRATOS_DEBUG_ERROR_IF(!this->GetGeometry().Has(TOTAL_MP_VOLUME)) << "The TOTAL_MP_VOLUME variable is not defined in the parent geometry." << std::endl;
    mGridVariables.mp_volume_sum = this->GetGeometry().GetValue(TOTAL_MP_VOLUME);
    mGridVariables.volume_ratio = mGridVariables.mp_volume_sum / mGridVariables.volume;
    // Calculate the scalar multiplier for soft stiffness (input scalar * empty_volume_ratio)
    mGridVariables.stiffness_multiplier = mGridVariables.penalty_factor * std::max(0.0, mGridVariables.volume - mGridVariables.mp_volume_sum) / mGridVariables.volume;

    if ((0.0 < mGridVariables.volume_ratio) && (mGridVariables.volume_ratio < GetProperties()[VOLUME_RATIO_THRESHOLD]))
    {
        this->Set(ACTIVE, true);
        KRATOS_INFO("MPMSoftStiffness")<<" Element: "<<this->Id()<< ", is ACTIVE. Volume ratio = " << mGridVariables.volume_ratio <<std::endl;
    }
    else if (mGridVariables.volume_ratio == 0 || mGridVariables.volume_ratio > GetProperties()[VOLUME_RATIO_THRESHOLD])
    {
        this->Set(ACTIVE, false);
        if (mGridVariables.volume_ratio > GetProperties()[VOLUME_RATIO_THRESHOLD])
            {
                KRATOS_INFO("MPMSoftStiffness")<<" Element: "<<this->Id()<< ", is NOT ACTIVE. Volume ratio = " << mGridVariables.volume_ratio <<std::endl;
            }
    }
    else
    {
        KRATOS_ERROR << "MPMSoftStiffness: Negative volume ratio in background grid! \nVolume ratio = "<< mGridVariables.volume_ratio << 
                        "\nTOTAL_MP_VOLUME = " << mGridVariables.mp_volume_sum << "\nGrid Volume = " << mGridVariables.volume << std::endl;
    }
    

    // if (mGridVariables.mp_volume_sum > 0.0){

    //     KRATOS_WATCH(this->Id())
    //     KRATOS_WATCH(mGridVariables.stiffness_multiplier)
    //     KRATOS_WATCH(mGridVariables.penalty_factor)
    //     KRATOS_WATCH(mGridVariables.volume)
    //     KRATOS_WATCH(mGridVariables.mp_volume_sum)
    // }
    mFinalizedStep = false;

}

void MPMSoftStiffness::InitializeMaterial(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    GeneralVariables Variables;

    if ( GetProperties()[CONSTITUTIVE_LAW] != NULL )
    {
        const GeometryType& rGeometry = GetGeometry();
        const Properties& rProperties = GetProperties();

        const Matrix& r_N = rGeometry.ShapeFunctionsValues();
        const auto& gp_size = this->GetGeometry().IntegrationPointsNumber();
        for ( IndexType gp_number = 0; gp_number < gp_size; ++gp_number ) {
            mConstitutiveLawVector[gp_number] = rProperties[CONSTITUTIVE_LAW]->Clone();
            mConstitutiveLawVector[gp_number]->InitializeMaterial( rProperties, rGeometry, row(r_N , gp_number ));
        }
    }
    else
        KRATOS_ERROR <<  "A constitutive law needs to be specified for the element with ID: " << this->Id() << std::endl;

    KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************

void MPMSoftStiffness::ResetConstitutiveLaw()
{
    KRATOS_TRY
    GeneralVariables Variables;

    if ( GetProperties()[CONSTITUTIVE_LAW] != NULL )
    {
        const GeometryType& rGeometry = GetGeometry();
        const Properties& rProperties = GetProperties();
        const Matrix& r_N = rGeometry.ShapeFunctionsValues();
        const auto& gp_size = this->GetGeometry().IntegrationPointsNumber();
        for ( IndexType gp_number = 0; gp_number < gp_size; ++gp_number )
            mConstitutiveLawVector[gp_number]->ResetMaterial( rProperties,  rGeometry, row( r_N, gp_number ) );
    }

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

double MPMSoftStiffness::CalculateIntegrationWeight(const double& rGaussWeight, const double detJ)
{
    const unsigned int dimension = this->GetGeometry().WorkingSpaceDimension();

    double integration_weight = rGaussWeight * detJ;

    // Add thickness if 2D
    KRATOS_DEBUG_ERROR_IF(dimension == 2 && !this->GetProperties().Has(THICKNESS)) << "The geometry is 2D but THICKNESS is not assigned." << std::endl; // TODO: Remove this
    if( dimension == 2 && this->GetProperties().Has( THICKNESS ))
        integration_weight *= this->GetProperties()[THICKNESS];

    return integration_weight;
}

//************************************************************************************
//************************************************************************************

void MPMSoftStiffness::EquationIdVector( EquationIdVectorType& rResult, const ProcessInfo& CurrentProcessInfo ) const
{
    const GeometryType& r_geometry = GetGeometry();
    int number_of_nodes = r_geometry.size();
    int dimension = r_geometry.WorkingSpaceDimension();
    unsigned int matrix_size = number_of_nodes * dimension;

    if ( rResult.size() != matrix_size )
        rResult.resize( matrix_size, false );

    for ( int i = 0; i < number_of_nodes; i++ )
    {
        int index = i * dimension;
        rResult[index] = r_geometry[i].GetDof( DISPLACEMENT_X ).EquationId();
        rResult[index + 1] = r_geometry[i].GetDof( DISPLACEMENT_Y ).EquationId();

        if ( dimension == 3 )
            rResult[index + 2] = r_geometry[i].GetDof( DISPLACEMENT_Z ).EquationId();
    }

}

//************************************************************************************
//************************************************************************************

void MPMSoftStiffness::GetDofList( DofsVectorType& rElementalDofList, const ProcessInfo& CurrentProcessInfo ) const
{
    const GeometryType& r_geometry = GetGeometry();
    rElementalDofList.resize( 0 );

    for ( unsigned int i = 0; i < r_geometry.size(); i++ )
    {
        rElementalDofList.push_back( r_geometry[i].pGetDof( DISPLACEMENT_X ) );
        rElementalDofList.push_back( r_geometry[i].pGetDof( DISPLACEMENT_Y ) );

        if ( r_geometry.WorkingSpaceDimension() == 3 )
        {
            rElementalDofList.push_back( r_geometry[i].pGetDof( DISPLACEMENT_Z ) );
        }
    }

}


//************************************************************************************
//************************************************************************************

void MPMSoftStiffness::GetValuesVector( Vector& values, int Step ) const
{
    const GeometryType& r_geometry = GetGeometry();
    const unsigned int number_of_nodes = r_geometry.size();
    const unsigned int dimension = r_geometry.WorkingSpaceDimension();
    unsigned int matrix_size = number_of_nodes * dimension;

    if ( values.size() != matrix_size ) values.resize( matrix_size, false );

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        unsigned int index = i * dimension;
        values[index] = r_geometry[i].FastGetSolutionStepValue( DISPLACEMENT_X, Step );
        values[index + 1] = r_geometry[i].FastGetSolutionStepValue( DISPLACEMENT_Y, Step );

        if ( dimension == 3 )
            values[index + 2] = r_geometry[i].FastGetSolutionStepValue( DISPLACEMENT_Z, Step );
    }
}

//************************************************************************************
//************************************************************************************

void MPMSoftStiffness::GetFirstDerivativesVector( Vector& values, int Step ) const
{
    const GeometryType& r_geometry = GetGeometry();
    const unsigned int number_of_nodes = r_geometry.size();
    const unsigned int dimension = r_geometry.WorkingSpaceDimension();
    unsigned int matrix_size = number_of_nodes * dimension;

    if ( values.size() != matrix_size ) values.resize( matrix_size, false );

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        unsigned int index = i * dimension;
        values[index] = r_geometry[i].FastGetSolutionStepValue( VELOCITY_X, Step );
        values[index + 1] = r_geometry[i].FastGetSolutionStepValue( VELOCITY_Y, Step );

        if ( dimension == 3 )
            values[index + 2] = r_geometry[i].FastGetSolutionStepValue( VELOCITY_Z, Step );
    }
}

//************************************************************************************
//************************************************************************************

void MPMSoftStiffness::GetSecondDerivativesVector( Vector& values, int Step ) const
{
    const GeometryType& r_geometry = GetGeometry();
    const unsigned int number_of_nodes = r_geometry.size();
    const unsigned int dimension = r_geometry.WorkingSpaceDimension();
    unsigned int matrix_size = number_of_nodes * dimension;

    if ( values.size() != matrix_size ) values.resize( matrix_size, false );

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        unsigned int index = i * dimension;
        values[index] = r_geometry[i].FastGetSolutionStepValue( ACCELERATION_X, Step );
        values[index + 1] = r_geometry[i].FastGetSolutionStepValue( ACCELERATION_Y, Step );

        if ( dimension == 3 )
            values[index + 2] = r_geometry[i].FastGetSolutionStepValue( ACCELERATION_Z, Step );
    }
}

//*************************DECIMAL CORRECTION OF STRAINS******************************
//************************************************************************************

void MPMSoftStiffness::DecimalCorrection(Vector& rVector)
{
    KRATOS_TRY

    for ( unsigned int i = 0; i < rVector.size(); i++ )
    {
        if( rVector[i]*rVector[i]<1e-24 )
        {
            rVector[i]=0;
        }
    }

    KRATOS_CATCH( "" )
}


///@}
///@name Access Get Values
///@{

void MPMSoftStiffness::CalculateOnIntegrationPoints(const Variable<int>& rVariable,
    std::vector<int>& rValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    if (rValues.size() != 1)
        rValues.resize(1);

    if (rVariable == MP_MATERIAL_ID) {
        rValues[0] = GetProperties().Id();
    }
    else
    {
        KRATOS_ERROR << "Variable " << rVariable << " is called in CalculateOnIntegrationPoints, but is not implemented." << std::endl;
    }
}

void MPMSoftStiffness::CalculateOnIntegrationPoints(const Variable<double>& rVariable,
    std::vector<double>& rValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    if (rValues.size() != 1)
        rValues.resize(1);

    
    if (rVariable == TOTAL_MP_VOLUME){
        rValues[0] = mGridVariables.mp_volume_sum;
    }
    else if (rVariable == PENALTY_FACTOR) {
        rValues[0] = mGridVariables.penalty_factor;
    }
    else if (rVariable == MP_MASS) {
        rValues[0] = mGridVariables.mass;
    }
    else
    {
        KRATOS_ERROR << "Variable " << rVariable << " is called in CalculateOnIntegrationPoints, but is not implemented." << std::endl;
    }
}


///@}
///@name Access Set Values
///@{

void MPMSoftStiffness::SetValuesOnIntegrationPoints(const Variable<int>& rVariable,
    const std::vector<int>& rValues,
    const ProcessInfo& rCurrentProcessInfo)
{
}

void MPMSoftStiffness::SetValuesOnIntegrationPoints(const Variable<double>& rVariable,
    const std::vector<double>& rValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_ERROR_IF(rValues.size() > 1)
        << "Only 1 value per integration point allowed! Passed values vector size: "
        << rValues.size() << std::endl;

    if (rVariable == MP_MASS) {
        mGridVariables.mass = rValues[0];
    }
    else if (rVariable == MP_DENSITY) {
        mGridVariables.density = rValues[0];
    }
    else if (rVariable == MP_VOLUME) {
        mGridVariables.volume = rValues[0];
    }
    else
    {
        KRATOS_ERROR << "Variable " << rVariable << " is called in SetValuesOnIntegrationPoints, but is not implemented." << std::endl;
    }
}

void MPMSoftStiffness::SetValuesOnIntegrationPoints(const Variable<Vector>& rVariable,
    const std::vector<Vector>& rValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_ERROR_IF(rValues.size() > 1)
        << "Only 1 value per integration point allowed! Passed values vector size: "
        << rValues.size() << std::endl;

    if (rVariable == MP_CAUCHY_STRESS_VECTOR) {
        mGridVariables.cauchy_stress_vector = rValues[0];
    }
    else if (rVariable == MP_ALMANSI_STRAIN_VECTOR) {
        mGridVariables.almansi_strain_vector = rValues[0];
    }
    else
    {
        KRATOS_ERROR << "Variable " << rVariable << " is called in SetValuesOnIntegrationPoints, but is not implemented." << std::endl;
    }
}

///@}

/**
 * This function provides the place to perform checks on the completeness of the input.
 * It is designed to be called only once (or anyway, not often) typically at the beginning
 * of the calculations, so to verify that nothing is missing from the input
 * or that no common error is found.
 * @param rCurrentProcessInfo
 */
int  MPMSoftStiffness::Check( const ProcessInfo& rCurrentProcessInfo ) const
{
    KRATOS_TRY

    Element::Check(rCurrentProcessInfo);

    const GeometryType& r_geometry = GetGeometry();
    const unsigned int number_of_nodes = r_geometry.size();
    const unsigned int dimension = r_geometry.WorkingSpaceDimension();

    // Verify compatibility with the constitutive law
    ConstitutiveLaw::Features LawFeatures;

    this->GetProperties().GetValue( CONSTITUTIVE_LAW )->GetLawFeatures(LawFeatures);

    const bool is_explicit = (rCurrentProcessInfo.Has(IS_EXPLICIT))
        ? rCurrentProcessInfo.GetValue(IS_EXPLICIT) : false;

    bool correct_strain_measure = false;
    for(unsigned int i=0; i<LawFeatures.mStrainMeasures.size(); i++)
    {
        if(LawFeatures.mStrainMeasures[i] == ConstitutiveLaw::StrainMeasure_Deformation_Gradient) correct_strain_measure = true;
        if (is_explicit && LawFeatures.mStrainMeasures[i] == ConstitutiveLaw::StrainMeasure_Velocity_Gradient) correct_strain_measure = true;

    }
    if (true)
    {

    }

    KRATOS_ERROR_IF(correct_strain_measure == false ) << "Constitutive law is not compatible with the element type: Large Displacements " << std::endl;

    // Verify that the dofs exist
    for ( IndexType i = 0; i < number_of_nodes; i++ ) {
        const NodeType &rnode = this->GetGeometry()[i];
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DISPLACEMENT,rnode)

        KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_X, rnode)
        KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_Y, rnode)
        KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_Z, rnode)
    }

    // Verify that the constitutive law exists
    if( this->GetProperties().Has( CONSTITUTIVE_LAW ) == false)
    {
        KRATOS_ERROR << "Constitutive law not provided for property " << this->GetProperties().Id() << std::endl;
    }
    else
    {
        // Verify that the constitutive law has the correct dimension
        if ( dimension == 2 )
        {
            KRATOS_ERROR_IF_NOT(this->GetProperties().Has( THICKNESS )) << "THICKNESS not provided for element " << this->Id() << std::endl;
        }
        else
        {
            KRATOS_ERROR_IF_NOT(this->GetProperties().GetValue( CONSTITUTIVE_LAW )->GetStrainSize() == 6) << "Wrong constitutive law used. This is a 3D element! expected strain size is 6 (el id = ) " << this->Id() << std::endl;
        }

        // Check constitutive law
        this->GetProperties().GetValue( CONSTITUTIVE_LAW )->Check( this->GetProperties(), r_geometry, rCurrentProcessInfo );
    }

    return 0;

    KRATOS_CATCH( "" );
}

void MPMSoftStiffness::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Element )

    rSerializer.save("ConstitutiveLawVector",mConstitutiveLawVector);
    rSerializer.save("MP",mGridVariables);
}

void MPMSoftStiffness::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Element )
    rSerializer.load("ConstitutiveLawVector",mConstitutiveLawVector);
    rSerializer.load("MP",mGridVariables);
}


} // Namespace Kratos

