//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		BSD License
//					Kratos default license: kratos/license.txt
//
//  Main authors:    Ilaria Iaconeta, Bodhinanda Chandra
//


// System includes
#include <omp.h>
#include <sstream>

// External includes

// Project includes
#include "includes/define.h"
#include "custom_elements/updated_lagrangian_axisymmetry.hpp"
#include "utilities/math_utils.h"
#include "includes/constitutive_law.h"
#include "particle_mechanics_application.h"

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

UpdatedLagrangianAxisymmetry::UpdatedLagrangianAxisymmetry( )
    : UpdatedLagrangian( )
{
    //DO NOT CALL IT: only needed for Register and Serialization!!!
}

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

UpdatedLagrangianAxisymmetry::UpdatedLagrangianAxisymmetry( IndexType NewId, GeometryType::Pointer pGeometry )
        : UpdatedLagrangian( NewId, pGeometry )
{
    //DO NOT ADD DOFS HERE!!!
}

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

UpdatedLagrangianAxisymmetry::UpdatedLagrangianAxisymmetry( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
        : UpdatedLagrangian( NewId, pGeometry, pProperties )
{
mFinalizedStep = true;
}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

UpdatedLagrangianAxisymmetry::UpdatedLagrangianAxisymmetry( UpdatedLagrangianAxisymmetry const& rOther)
    :UpdatedLagrangian(rOther)
{
}

//*********************************OPERATIONS*****************************************
//************************************************************************************

Element::Pointer UpdatedLagrangianAxisymmetry::Create( IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties ) const
{
    return Element::Pointer( new UpdatedLagrangianAxisymmetry( NewId, GetGeometry().Create( ThisNodes ), pProperties ) );
}

//************************************CLONE*******************************************
//************************************************************************************

Element::Pointer UpdatedLagrangianAxisymmetry::Clone( IndexType NewId, NodesArrayType const& rThisNodes ) const
{

    UpdatedLagrangianAxisymmetry NewElement (NewId, GetGeometry().Create( rThisNodes ), pGetProperties() );

    return Element::Pointer( new UpdatedLagrangianAxisymmetry(NewElement) );
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

UpdatedLagrangianAxisymmetry::~UpdatedLagrangianAxisymmetry()
{
}

//************************************************************************************
void UpdatedLagrangianAxisymmetry::Initialize()
{
    KRATOS_TRY
    UpdatedLagrangian::Initialize();

    const array_1d<double,3>& xg = this->GetValue(GAUSS_COORD);
    const double mp_volume = this->GetValue(MP_VOLUME);
    const double mp_mass = mp_volume * 2* GetPI() * xg[0] * GetProperties()[DENSITY];
    this->SetValue(MP_MASS, mp_mass);

    mDeterminantF0 = 1;
    mDeformationGradientF0 = IdentityMatrix(3);

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void UpdatedLagrangianAxisymmetry::InitializeGeneralVariables (GeneralVariables& rVariables, const ProcessInfo& rCurrentProcessInfo)
{
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    const unsigned int strain_size = GetProperties().GetValue( CONSTITUTIVE_LAW )->GetStrainSize();

    rVariables.detF  = 1;

    rVariables.detF0 = 1;

    rVariables.detFT = 1;

    rVariables.detJ = 1;

    rVariables.B.resize( strain_size , number_of_nodes * dimension, false );

    rVariables.F.resize( 3, 3, false );
    rVariables.F = IdentityMatrix(3);

    rVariables.F0.resize( 3, 3, false );
    rVariables.F0 = IdentityMatrix(3);

    rVariables.FT.resize( 3, 3, false );
    rVariables.FT = IdentityMatrix(3);

    rVariables.ConstitutiveMatrix.resize( strain_size, strain_size, false );

    rVariables.StrainVector.resize( strain_size, false );

    rVariables.StressVector.resize( strain_size, false );

    rVariables.DN_DX.resize( number_of_nodes, dimension, false );
    rVariables.DN_De.resize( number_of_nodes, dimension, false );

    array_1d<double,3>& xg = this->GetValue(GAUSS_COORD);

    rVariables.N = this->MPMShapeFunctionPointValues(rVariables.N, xg);

    // Reading shape functions local gradients
    rVariables.DN_De = this->MPMShapeFunctionsLocalGradients( rVariables.DN_De);

    // CurrentDisp is the unknown variable. It represents the nodal delta displacement. When it is predicted is equal to zero.
    rVariables.CurrentDisp = CalculateCurrentDisp(rVariables.CurrentDisp, rCurrentProcessInfo);

    // Calculating the current jacobian from cartesian coordinates to parent coordinates for the MP element [dx_n+1/d£]
    rVariables.j = this->MPMJacobianDelta( rVariables.j, xg, rVariables.CurrentDisp);

    // Calculating the reference jacobian from cartesian coordinates to parent coordinates for the MP element [dx_n/d£]
    rVariables.J = this->MPMJacobian( rVariables.J, xg);

}

//************************************************************************************
//************************************************************************************

void UpdatedLagrangianAxisymmetry::CalculateElementalSystem( LocalSystemComponents& rLocalSystem,
                                ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Create and initialize element variables:
    GeneralVariables Variables;
    this->InitializeGeneralVariables(Variables,rCurrentProcessInfo);

    // Create constitutive law parameters:
    ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

    // Set constitutive law flags:
    Flags &ConstitutiveLawOptions=Values.GetOptions();
    ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);
    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS);
    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);

    // Auxiliary terms
    Vector volume_force;

    // Compute element kinematics B, F, DN_DX ...
    this->CalculateKinematics(Variables,rCurrentProcessInfo);

    // Set general variables to constitutivelaw parameters
    this->SetGeneralVariables(Variables,Values);

    // Calculate Material Response
    /* NOTE:
    The function below will call CalculateMaterialResponseCauchy() by default and then (may)
    call CalculateMaterialResponseKirchhoff() in the constitutive_law.*/
    mConstitutiveLawVector->CalculateMaterialResponse(Values, Variables.StressMeasure);

    /* NOTE:
    The material points will have constant mass as defined at the beginning.
    However, the density and volume (integration weight) are changing every time step.*/
    // Update MP_Density
    double MP_Density = (GetProperties()[DENSITY]) / Variables.detFT;
    this->SetValue(MP_DENSITY, MP_Density);

    // The MP_Volume (integration weight) is evaluated
    double MP_Volume = this->GetValue(MP_MASS)/this->GetValue(MP_DENSITY);
    this->SetValue(MP_VOLUME, MP_Volume);

    if ( rLocalSystem.CalculationFlags.Is(UpdatedLagrangian::COMPUTE_LHS_MATRIX) ) //calculation of the matrix is required
    {
        // Contributions to stiffness matrix calculated on the reference config
        this->CalculateAndAddLHS ( rLocalSystem, Variables, MP_Volume );
    }

    if ( rLocalSystem.CalculationFlags.Is(UpdatedLagrangian::COMPUTE_RHS_VECTOR) ) //calculation of the vector is required
    {
        // Contribution to external forces
        volume_force  = this->CalculateVolumeForce( volume_force, Variables );
        this->CalculateAndAddRHS ( rLocalSystem, Variables, volume_force, MP_Volume );

    }

    KRATOS_CATCH( "" )
}

/**
 * Calculates the radius of axisymmetry
 * @param N: The Gauss Point shape function
 * @param Geom: The geometry studied
 * @return Radius: The radius of axisymmetry
 */

double UpdatedLagrangianAxisymmetry::CalculateRadius(
    Vector& N,
    GeometryType& Geom,
    std::string ThisConfiguration = "Current"
    )
{
    double radius = 0.0;

    for (unsigned int iNode = 0; iNode < Geom.size(); iNode++)
    {
        // Displacement from the reference to the current configuration
        if (ThisConfiguration == "Current")
        {
            const array_1d<double, 3 > DeltaDisplacement = Geom[iNode].FastGetSolutionStepValue(DISPLACEMENT);
            const array_1d<double, 3 > ReferencePosition = Geom[iNode].Coordinates();

            const array_1d<double, 3 > CurrentPosition = ReferencePosition + DeltaDisplacement;
            radius += CurrentPosition[0] * N[iNode];
        }
        else
        {
            const array_1d<double, 3 > ReferencePosition = Geom[iNode].Coordinates();
            radius += ReferencePosition[0] * N[iNode];
        }
    }

    return radius;
}


//*********************************COMPUTE KINEMATICS*********************************
//************************************************************************************

void UpdatedLagrangianAxisymmetry::CalculateKinematics(GeneralVariables& rVariables, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Define the stress measure
    rVariables.StressMeasure = ConstitutiveLaw::StressMeasure_Cauchy;

    // Calculating the inverse of the jacobian and the parameters needed [d£/dx_n]
    Matrix InvJ;
    MathUtils<double>::InvertMatrix( rVariables.J, InvJ, rVariables.detJ);

    rVariables.CurrentDisp = CalculateCurrentDisp(rVariables.CurrentDisp, rCurrentProcessInfo);

    // Compute cartesian derivatives [dN/dx_n]
    rVariables.DN_DX = prod( rVariables.DN_De, InvJ);

    // Compute radius
    const double current_radius = CalculateRadius(rVariables.N, GetGeometry(), "Current");
    const double initial_radius = CalculateRadius(rVariables.N, GetGeometry(), "Initial");

    CalculateDeformationGradient (rVariables.DN_DX, rVariables.F, rVariables.CurrentDisp, current_radius, initial_radius);

    //Calculating the inverse of the jacobian and the parameters needed [d£/(dx_n+1)]
    Matrix Invj;
    MathUtils<double>::InvertMatrix( rVariables.j, Invj, rVariables.detJ ); //overwrites detJ

    // Compute cartesian derivatives [dN/dx_n+1]
    rVariables.DN_DX = prod( rVariables.DN_De, Invj); //overwrites DX now is the current position dx

    rVariables.CurrentRadius   =  current_radius;
    rVariables.ReferenceRadius = initial_radius;

    // Determinant of the previous Deformation Gradient F_n
    rVariables.detF0 = mDeterminantF0;
    rVariables.F0    = mDeformationGradientF0;

    // Compute the deformation matrix B
    this->CalculateDeformationMatrix(rVariables.B, rVariables.F, rVariables.DN_DX, rVariables.N);

    KRATOS_CATCH( "" )
}

//************************************************************************************

void UpdatedLagrangianAxisymmetry::CalculateDeformationMatrix(Matrix& rB,
        Matrix& rF,
        Matrix& rDN_DX,
        Vector& rN)
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
    const double radius = CalculateRadius(rN, GetGeometry(), "Current");

    rB.clear(); // Set all components to zero

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        const unsigned int index = dimension * i;

        rB( 0, index + 0 ) = rDN_DX( i, 0 );
        rB( 1, index + 1 ) = rDN_DX( i, 1 );
        rB( 2, index + 0 ) = rN[i]/radius;
        rB( 3, index + 0 ) = rDN_DX( i, 1 );
        rB( 3, index + 1 ) = rDN_DX( i, 0 );
    }

    KRATOS_CATCH( "" )
}


//*************************COMPUTE DEFORMATION GRADIENT*******************************
//************************************************************************************

void UpdatedLagrangianAxisymmetry::CalculateDeformationGradient(const Matrix& rDN_DX,
        Matrix&  rF,
        Matrix&  rDeltaPosition,
        const double & rCurrentRadius,
        const double & rReferenceRadius)
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

    rF = IdentityMatrix( 3 );

    if( dimension == 2 )
    {
        for ( unsigned int i = 0; i < number_of_nodes; i++ )
        {
            rF ( 0 , 0 ) += rDeltaPosition(i,0)*rDN_DX ( i , 0 );
            rF ( 0 , 1 ) += rDeltaPosition(i,0)*rDN_DX ( i , 1 );
            rF ( 1 , 0 ) += rDeltaPosition(i,1)*rDN_DX ( i , 0 );
            rF ( 1 , 1 ) += rDeltaPosition(i,1)*rDN_DX ( i , 1 );
        }

        rF ( 2 , 2 ) = rCurrentRadius/rReferenceRadius;
    }
    else
    {
        KRATOS_ERROR <<  "Dimension given is wrong: Something is wrong with the given dimension in function: CalculateDeformationGradient" << std::endl;
    }

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void UpdatedLagrangianAxisymmetry::CalculateAndAddKuug(MatrixType& rK,
        GeneralVariables& rVariables,
        double& rIntegrationWeight)
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().size();

    // Axisymmetric geometric matrix
    double alpha_1 = 0;
    double alpha_2 = 0;
    double alpha_3 = 0;

    unsigned int index_i = 0;
    unsigned int index_j = 0;

    const double radius = CalculateRadius(rVariables.N, GetGeometry(), "Current");

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        index_j =0;
        for ( unsigned int j = 0; j < number_of_nodes; j++ )
        {
            alpha_1 = rVariables.DN_DX(j,0) * ( rVariables.DN_DX(i,0) * rVariables.StressVector[0] + rVariables.DN_DX(i,1) * rVariables.StressVector[3] );
            alpha_2 = rVariables.DN_DX(j,1) * ( rVariables.DN_DX(i,0) * rVariables.StressVector[3] + rVariables.DN_DX(i,1) * rVariables.StressVector[1] );
            alpha_3 = rVariables.N[i] * rVariables.N[j] * rVariables.StressVector[2] * (1.0/radius*radius);

            rK(index_i,index_j)     += (alpha_1 + alpha_2 + alpha_3) * rIntegrationWeight ;
            rK(index_i+1,index_j+1) += (alpha_1 + alpha_2) * rIntegrationWeight ;

            index_j+=2;
        }

        index_i+=2;

    }

    KRATOS_CATCH( "" )
}

//************************************CALCULATE VOLUME CHANGE*************************
//************************************************************************************

double& UpdatedLagrangianAxisymmetry::CalculateVolumeChange( double& rVolumeChange, GeneralVariables& rVariables )
{
    KRATOS_TRY

    rVolumeChange = 1.0 / (rVariables.detF * rVariables.detF0);
    return rVolumeChange;

    KRATOS_CATCH( "" )
}

//************************************CALCULATE VOLUME ACCELERATION*******************
//************************************************************************************

Vector& UpdatedLagrangianAxisymmetry::CalculateVolumeForce( Vector& rVolumeForce, GeneralVariables& rVariables )
{
    KRATOS_TRY

    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    rVolumeForce = ZeroVector(dimension);
    rVolumeForce = this->GetValue(MP_VOLUME_ACCELERATION) * this->GetValue(MP_MASS);

    return rVolumeForce;

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void UpdatedLagrangianAxisymmetry::CalculateRightHandSide( VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo )
{
    // Create local system components
    LocalSystemComponents LocalSystem;

    // Calculation flags
    LocalSystem.CalculationFlags.Set(UpdatedLagrangian::COMPUTE_RHS_VECTOR);

    MatrixType LeftHandSideMatrix = Matrix();

    // Initialize sizes for the system components:
    this->InitializeSystemMatrices( LeftHandSideMatrix, rRightHandSideVector, LocalSystem.CalculationFlags );

    // Set Variables to Local system components
    LocalSystem.SetLeftHandSideMatrix(LeftHandSideMatrix);
    LocalSystem.SetRightHandSideVector(rRightHandSideVector);

    // Calculate elemental system
    CalculateElementalSystem( LocalSystem, rCurrentProcessInfo );

}

//************************************************************************************
//************************************************************************************

void UpdatedLagrangianAxisymmetry::CalculateRightHandSide( std::vector< VectorType >& rRightHandSideVectors, const std::vector< Variable< VectorType > >& rRHSVariables, ProcessInfo& rCurrentProcessInfo )
{
    // Create local system components
    LocalSystemComponents LocalSystem;

    // Calculation flags
    LocalSystem.CalculationFlags.Set(UpdatedLagrangian::COMPUTE_RHS_VECTOR);
    LocalSystem.CalculationFlags.Set(UpdatedLagrangian::COMPUTE_RHS_VECTOR_WITH_COMPONENTS);

    MatrixType LeftHandSideMatrix = Matrix();

    // Initialize sizes for the system components:
    if( rRHSVariables.size() != rRightHandSideVectors.size() )
        rRightHandSideVectors.resize(rRHSVariables.size());

    for( unsigned int i=0; i<rRightHandSideVectors.size(); i++ )
        this->InitializeSystemMatrices( LeftHandSideMatrix, rRightHandSideVectors[i], LocalSystem.CalculationFlags );

    // Set Variables to Local system components
    LocalSystem.SetLeftHandSideMatrix(LeftHandSideMatrix);
    LocalSystem.SetRightHandSideVectors(rRightHandSideVectors);

    LocalSystem.SetRightHandSideVariables(rRHSVariables);

    //Calculate elemental system
    CalculateElementalSystem( LocalSystem, rCurrentProcessInfo );
}


//************************************************************************************
//************************************************************************************

void UpdatedLagrangianAxisymmetry::CalculateLeftHandSide( MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo )
{
    // Create local system components
    LocalSystemComponents LocalSystem;

    // Calculation flags
    LocalSystem.CalculationFlags.Set(UpdatedLagrangian::COMPUTE_LHS_MATRIX);

    VectorType RightHandSideVector = Vector();

    // Initialize sizes for the system components:
    this->InitializeSystemMatrices( rLeftHandSideMatrix, RightHandSideVector, LocalSystem.CalculationFlags );

    // Set Variables to Local system components
    LocalSystem.SetLeftHandSideMatrix(rLeftHandSideMatrix);
    LocalSystem.SetRightHandSideVector(RightHandSideVector);

    // Calculate elemental system
    CalculateElementalSystem( LocalSystem, rCurrentProcessInfo );
}

//************************************************************************************
//************************************************************************************

void UpdatedLagrangianAxisymmetry::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo )
{
    // Create local system components
    LocalSystemComponents LocalSystem;

    // Set calculation flags
    LocalSystem.CalculationFlags.Set(UpdatedLagrangian::COMPUTE_LHS_MATRIX);
    LocalSystem.CalculationFlags.Set(UpdatedLagrangian::COMPUTE_RHS_VECTOR);

    // Initialize sizes for the system components:
    this->InitializeSystemMatrices( rLeftHandSideMatrix, rRightHandSideVector, LocalSystem.CalculationFlags );

    // Set Variables to Local system components
    LocalSystem.SetLeftHandSideMatrix(rLeftHandSideMatrix);
    LocalSystem.SetRightHandSideVector(rRightHandSideVector);

    // Calculate elemental system
    CalculateElementalSystem( LocalSystem, rCurrentProcessInfo );
}

//************************************************************************************
//************************************************************************************

void UpdatedLagrangianAxisymmetry::CalculateLocalSystem( std::vector< MatrixType >& rLeftHandSideMatrices,
                                const std::vector< Variable< MatrixType > >& rLHSVariables,
                                std::vector< VectorType >& rRightHandSideVectors,
                                const std::vector< Variable< VectorType > >& rRHSVariables,
                                ProcessInfo& rCurrentProcessInfo )
{
    // Create local system components
    LocalSystemComponents LocalSystem;

    // Calculation flags
    LocalSystem.CalculationFlags.Set(UpdatedLagrangian::COMPUTE_LHS_MATRIX_WITH_COMPONENTS);
    LocalSystem.CalculationFlags.Set(UpdatedLagrangian::COMPUTE_RHS_VECTOR_WITH_COMPONENTS);

    // Initialize sizes for the system components:
    if( rLHSVariables.size() != rLeftHandSideMatrices.size() )
        rLeftHandSideMatrices.resize(rLHSVariables.size());

    if( rRHSVariables.size() != rRightHandSideVectors.size() )
        rRightHandSideVectors.resize(rRHSVariables.size());

    LocalSystem.CalculationFlags.Set(UpdatedLagrangian::COMPUTE_LHS_MATRIX);
    for( unsigned int i=0; i<rLeftHandSideMatrices.size(); i++ )
        this->InitializeSystemMatrices( rLeftHandSideMatrices[i], rRightHandSideVectors[0], LocalSystem.CalculationFlags );

    LocalSystem.CalculationFlags.Set(UpdatedLagrangian::COMPUTE_RHS_VECTOR);
    LocalSystem.CalculationFlags.Set(UpdatedLagrangian::COMPUTE_LHS_MATRIX,false);

    for( unsigned int i=0; i<rRightHandSideVectors.size(); i++ )
        this->InitializeSystemMatrices( rLeftHandSideMatrices[0], rRightHandSideVectors[i], LocalSystem.CalculationFlags );

    LocalSystem.CalculationFlags.Set(UpdatedLagrangian::COMPUTE_LHS_MATRIX,true);

    // Set Variables to Local system components
    LocalSystem.SetLeftHandSideMatrices(rLeftHandSideMatrices);
    LocalSystem.SetRightHandSideVectors(rRightHandSideVectors);

    LocalSystem.SetLeftHandSideVariables(rLHSVariables);
    LocalSystem.SetRightHandSideVariables(rRHSVariables);

    //Calculate elemental system
    CalculateElementalSystem( LocalSystem, rCurrentProcessInfo );

}

//***********************************************************************************
//***********************************************************************************
void UpdatedLagrangianAxisymmetry::Calculate(const Variable<double>& rVariable,
                        double& Output,
                        const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    GeneralVariables Variables;
    this->InitializeGeneralVariables(Variables,rCurrentProcessInfo);

    if (rVariable == DENSITY)
    {
        const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
        const unsigned int number_of_nodes = GetGeometry().PointsNumber();
        array_1d<double,3>& xg = this->GetValue(GAUSS_COORD);
        GeneralVariables Variables;
        Matrix J0 = ZeroMatrix(dimension);

        J0 = this->MPMJacobian(J0, xg);

        // Calculating and storing inverse and the determinant of the jacobian
        MathUtils<double>::InvertMatrix( J0, mInverseJ0, mDeterminantJ0 );

        Variables.N = this->MPMShapeFunctionPointValues(Variables.N, xg);
        double MP_Mass = this->GetValue(MP_MASS);

        for (unsigned int i=0;i<number_of_nodes;i++)
                GetGeometry()[i].GetSolutionStepValue(AUX_R) += Variables.N[i] * (MP_Mass);// - AUX_MP_Mass);

    }

    KRATOS_CATCH( "" )
}


void UpdatedLagrangianAxisymmetry::Calculate(const Variable<array_1d<double, 3 > >& rVariable,
                        array_1d<double, 3 > & Output,
                        const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    array_1d<double,3>& xg = this->GetValue(GAUSS_COORD);
    GeneralVariables Variables;
    Matrix J0 = ZeroMatrix(dimension);

    J0 = this->MPMJacobian(J0, xg);

    // Calculating and storing inverse and the determinant of the jacobian
    MathUtils<double>::InvertMatrix( J0, mInverseJ0, mDeterminantJ0 );

    Variables.N = this->MPMShapeFunctionPointValues(Variables.N, xg);
    double MP_Mass = this->GetValue(MP_MASS);

    if(rVariable == VELOCITY)
    {
        array_1d<double,3>& MP_Velocity = this->GetValue(MP_VELOCITY);
        array_1d<double,3> NodalAuxRVel;

        for (unsigned int i=0;i<number_of_nodes;i++)
        {
            for (unsigned int j = 0; j < dimension; j++)
                NodalAuxRVel[j] = Variables.N[i] * MP_Mass * MP_Velocity[j];

            GetGeometry()[i].GetSolutionStepValue(AUX_R_VEL) += NodalAuxRVel;
        }
    }
    else if(rVariable == ACCELERATION)
    {
        array_1d<double,3>& MP_Acceleration = this->GetValue(MP_ACCELERATION);
        array_1d<double,3> NodalAuxRAcc;

        for (unsigned int i=0;i<number_of_nodes;i++)
        {
            for (unsigned int j = 0; j < dimension; j++)
                NodalAuxRAcc[j] = Variables.N[i] * MP_Mass * MP_Acceleration[j];

            GetGeometry()[i].GetSolutionStepValue(AUX_R_ACC) += NodalAuxRAcc;
        }
    }

    KRATOS_CATCH( "" )
}

//*******************************************************************************************
//*******************************************************************************************

void UpdatedLagrangianAxisymmetry::InitializeSolutionStep( ProcessInfo& rCurrentProcessInfo )
{
    /* NOTE:
    In the InitializeSolutionStep of each time step the nodal initial conditions are evaluated.
    This function is called by the base scheme class.*/

    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    array_1d<double,3>& xg = this->GetValue(GAUSS_COORD);
    GeneralVariables Variables;

    // Calculating and storing inverse and the determinant of the jacobian
    Matrix J0 = ZeroMatrix(dimension);
    J0 = this->MPMJacobian(J0, xg);
    MathUtils<double>::InvertMatrix( J0, mInverseJ0, mDeterminantJ0 );

    // Calculating shape function
    Variables.N = this->MPMShapeFunctionPointValues(Variables.N, xg);

    // Initialize Constitutive Law
    mConstitutiveLawVector->InitializeSolutionStep( GetProperties(),
            GetGeometry(), Variables.N, rCurrentProcessInfo );

    mFinalizedStep = false;

    const array_1d<double,3>& MP_Velocity = this->GetValue(MP_VELOCITY);
    const array_1d<double,3>& MP_Acceleration = this->GetValue(MP_ACCELERATION);
    array_1d<double,3>& AUX_MP_Velocity = this->GetValue(AUX_MP_VELOCITY);
    array_1d<double,3>& AUX_MP_Acceleration = this->GetValue(AUX_MP_ACCELERATION);
    const double MP_Mass = this->GetValue(MP_MASS);
    array_1d<double,3> MP_Momentum;
    array_1d<double,3> MP_Inertia;
    array_1d<double,3> nodal_momentum= ZeroVector(3);
    array_1d<double,3> nodal_inertia = ZeroVector(3);

    for (unsigned int j=0;j<number_of_nodes;j++)
    {
        // These are the values of nodal velocity and nodal acceleration evaluated in the initialize solution step
        for (unsigned int l=0;l<number_of_nodes;l++)
        {
            array_1d<double, 3 > & NodalAcceleration = GetGeometry()[l].FastGetSolutionStepValue(ACCELERATION,1);
            array_1d<double, 3 > & NodalVelocity = GetGeometry()[l].FastGetSolutionStepValue(VELOCITY,1);

            for (unsigned int k = 0; k < dimension; k++)
            {
                AUX_MP_Velocity[k] += Variables.N[j] *Variables.N[l] * NodalVelocity[k];
                AUX_MP_Acceleration[k] += Variables.N[j] *Variables.N[l] * NodalAcceleration[k];
            }
        }
    }

    // Here MP contribution in terms of momentum, inertia and mass are added
    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        for (unsigned int j = 0; j < dimension; j++)
        {
            nodal_momentum[j] = Variables.N[i] * (MP_Velocity[j] - AUX_MP_Velocity[j]) * MP_Mass;
            nodal_inertia[j] = Variables.N[i] * (MP_Acceleration[j] - AUX_MP_Acceleration[j]) * MP_Mass;
        }

        GetGeometry()[i].SetLock();
        GetGeometry()[i].GetSolutionStepValue(NODAL_MOMENTUM, 0) += nodal_momentum;
        GetGeometry()[i].GetSolutionStepValue(NODAL_INERTIA, 0) += nodal_inertia;
        GetGeometry()[i].GetSolutionStepValue(NODAL_MASS, 0) += Variables.N[i] * MP_Mass;
        GetGeometry()[i].UnSetLock();

    }

    AUX_MP_Velocity.clear();
    AUX_MP_Acceleration.clear();

}

//************************************************************************************
//************************************************************************************

void UpdatedLagrangianAxisymmetry::InitializeNonLinearIteration( ProcessInfo& rCurrentProcessInfo )
{
}

//************************************************************************************
//************************************************************************************

void UpdatedLagrangianAxisymmetry::FinalizeNonLinearIteration( ProcessInfo& rCurrentProcessInfo )
{
}

//************************************************************************************************************

void UpdatedLagrangianAxisymmetry::FinalizeStepVariables( GeneralVariables & rVariables, const ProcessInfo& rCurrentProcessInfo)
{
    // Update internal (historical) variables
    mDeterminantF0         = rVariables.detF* rVariables.detF0;
    mDeformationGradientF0 = prod(rVariables.F, rVariables.F0);

    this->SetValue(MP_CAUCHY_STRESS_VECTOR, rVariables.StressVector);
    this->SetValue(MP_ALMANSI_STRAIN_VECTOR, rVariables.StrainVector);
    this->SetValue(MP_JACOBIAN, mDeterminantF0);

    double EquivalentPlasticStrain = mConstitutiveLawVector->GetValue(PLASTIC_STRAIN, EquivalentPlasticStrain );
    this->SetValue(MP_EQUIVALENT_PLASTIC_STRAIN, EquivalentPlasticStrain);

    double Pressure = mConstitutiveLawVector->GetValue(PRESSURE, Pressure );
    this->SetValue(MP_PRESSURE, Pressure);

    MathUtils<double>::InvertMatrix( rVariables.j, mInverseJ, rVariables.detJ );

    this->UpdateGaussPoint(rVariables, rCurrentProcessInfo);
}

//************************************************************************************
//************************************************************************************
/**
 * The position of the Gauss points/Material points is updated
 */

void UpdatedLagrangianAxisymmetry::UpdateGaussPoint( GeneralVariables & rVariables, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    rVariables.CurrentDisp = CalculateCurrentDisp(rVariables.CurrentDisp, rCurrentProcessInfo);
    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    array_1d<double,3> xg = this->GetValue(GAUSS_COORD);
    const array_1d<double,3> MP_PreviousAcceleration = this->GetValue(MP_ACCELERATION);
    const array_1d<double,3> MP_PreviousVelocity = this->GetValue(MP_VELOCITY);

    array_1d<double,3> delta_xg = ZeroVector(3);
    array_1d<double,3> MP_Acceleration = ZeroVector(3);
    array_1d<double,3> MP_Velocity = ZeroVector(3);
    const double delta_time = rCurrentProcessInfo[DELTA_TIME];

    array_1d<double,3> MP_DeltaVelocity = ZeroVector(3);
    rVariables.N = this->MPMShapeFunctionPointValues(rVariables.N, xg);
    int MP_number = this->GetValue(MP_NUMBER);

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        if (rVariables.N[i] > 1.e-16)
        {
            array_1d<double, 3 > & NodalAcceleration = GetGeometry()[i].FastGetSolutionStepValue(ACCELERATION);
            array_1d<double, 3 > & NodalVelocity = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY);
            array_1d<double, 3 > & PreviousNodalVelocity = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY,1);
            double NodalMass = GetGeometry()[i].GetSolutionStepValue(NODAL_MASS, 0);
            array_1d<double,3> nodal_momentum = NodalMass * NodalVelocity;
            array_1d<double,3> nodal_inertia = NodalMass * NodalAcceleration;

            for ( unsigned int j = 0; j < dimension; j++ )
            {
                delta_xg[j] += rVariables.N[i] * rVariables.CurrentDisp(i,j);
                MP_Acceleration[j] += rVariables.N[i] * NodalAcceleration[j];

                /* NOTE: The following interpolation techniques have been tried:
                    MP_Velocity[j]      += rVariables.N[i] * nodal_velocity[j];
                    MP_Acceleration[j]  += nodal_inertia[j]/(rVariables.N[i] * MP_Mass * MP_number);
                    MP_Velocity[j]      += nodal_momentum[j]/(rVariables.N[i] * MP_Mass * MP_number);
                    MP_Velocity[j]      += delta_time * rVariables.N[i] * nodal_acceleration[j];
                */
            }
        }

    }

    /* NOTE:
    Another way to update the MP velocity (see paper Guilkey and Weiss, 2003).
    This assume newmark (or trapezoidal, since n.gamma=0.5) rule of integration*/
    MP_Velocity = MP_PreviousVelocity + 0.5 * delta_time * (MP_Acceleration + MP_PreviousAcceleration);
    this -> SetValue(MP_VELOCITY,MP_Velocity );

    /* NOTE: The following interpolation techniques have been tried:
        MP_Acceleration = 4/(delta_time * delta_time) * delta_xg - 4/delta_time * MP_PreviousVelocity;
        MP_Velocity = 2.0/delta_time * delta_xg - MP_PreviousVelocity;
    */

    // Update the MP Position
    const array_1d<double,3>& new_xg = xg + delta_xg ;
    this -> SetValue(GAUSS_COORD,new_xg);

    // Update the MP Acceleration
    this -> SetValue(MP_ACCELERATION,MP_Acceleration);

    // Update the MP total Displacement
    array_1d<double,3>& MP_Displacement = this->GetValue(MP_DISPLACEMENT);
    MP_Displacement += delta_xg;
    this -> SetValue(MP_DISPLACEMENT,MP_Displacement );

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void UpdatedLagrangianAxisymmetry::InitializeMaterial()
{
    KRATOS_TRY
    array_1d<double,3>& xg = this->GetValue(GAUSS_COORD);
    GeneralVariables Variables;

    if ( GetProperties()[CONSTITUTIVE_LAW] != NULL )
    {
        mConstitutiveLawVector = GetProperties()[CONSTITUTIVE_LAW]->Clone();

        Variables.N = this->MPMShapeFunctionPointValues(Variables.N, xg);

        mConstitutiveLawVector->InitializeMaterial( GetProperties(), GetGeometry(),
                Variables.N );
    }

    else
        KRATOS_ERROR <<  "A constitutive law needs to be specified for the element with ID: " << this->Id() << std::endl;

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void UpdatedLagrangianAxisymmetry::ResetConstitutiveLaw()
{
    KRATOS_TRY
    array_1d<double,3>& xg = this->GetValue(GAUSS_COORD);
    GeneralVariables Variables;

    // Create and initialize element variables:
    if ( GetProperties()[CONSTITUTIVE_LAW] != NULL )
    {
        mConstitutiveLawVector->ResetMaterial( GetProperties(), GetGeometry(), this->MPMShapeFunctionPointValues(Variables.N, xg) );
    }

    KRATOS_CATCH( "" )
}

//*************************COMPUTE CURRENT DISPLACEMENT*******************************
//************************************************************************************

Matrix& UpdatedLagrangianAxisymmetry::CalculateCurrentDisp(Matrix & rCurrentDisp, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    rCurrentDisp = zero_matrix<double>( number_of_nodes , dimension);

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        array_1d<double, 3 > & CurrentDisplacement  = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT);

        for ( unsigned int j = 0; j < dimension; j++ )
        {
            rCurrentDisp(i,j) = CurrentDisplacement[j];
        }
    }

    return rCurrentDisp;

    KRATOS_CATCH( "" )
}


//*************************COMPUTE ALMANSI STRAIN*************************************
//************************************************************************************
void UpdatedLagrangianAxisymmetry::CalculateAlmansiStrain(const Matrix& rF,
    Vector& rStrainVector )
{
    KRATOS_TRY

    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    // Left Cauchy-Green Calculation
    Matrix LeftCauchyGreen = prod( rF, trans( rF ) );

    // Calculating the inverse of the jacobian
    Matrix InverseLeftCauchyGreen ( dimension, dimension );
    double det_b=0;
    MathUtils<double>::InvertMatrix( LeftCauchyGreen, InverseLeftCauchyGreen, det_b);

    if( dimension == 2 )
    {
        // Almansi Strain Calculation
        rStrainVector[0] = 0.5 * (  1.00 - InverseLeftCauchyGreen( 0, 0 ) );
        rStrainVector[1] = 0.5 * (  1.00 - InverseLeftCauchyGreen( 1, 1 ) );
        rStrainVector[2] = - InverseLeftCauchyGreen( 0, 1 ); // xy
    }
    else
    {
        KRATOS_ERROR <<  "Dimension given is wrong: Something is wrong with the given dimension in function: CalculateAlmansiStrain" << std::endl;
    }


    KRATOS_CATCH( "" )
}

//*************************COMPUTE GREEN-LAGRANGE STRAIN*************************************
//************************************************************************************
// Green-Lagrange Strain: E = 0.5 * (U^2 - I) = 0.5 * (C - I)
void UpdatedLagrangianAxisymmetry::CalculateGreenLagrangeStrain(const Matrix& rF,
    Vector& rStrainVector )
{
    KRATOS_TRY

    const unsigned int dimension  = GetGeometry().WorkingSpaceDimension();

    // Right Cauchy-Green Calculation
    Matrix C ( dimension, dimension );
    noalias( C ) = prod( trans( rF ), rF );

    if( dimension == 2 )
    {
        //Green Lagrange Strain Calculation
        if ( rStrainVector.size() != 3 )
            rStrainVector.resize( 3, false );

        rStrainVector[0] = 0.5 * ( C( 0, 0 ) - 1.00 );
        rStrainVector[1] = 0.5 * ( C( 1, 1 ) - 1.00 );
        rStrainVector[2] = C( 0, 1 ); // xy
    }
    else
    {
        KRATOS_ERROR <<  "Dimension given is wrong: Something is wrong with the given dimension in function: CalculateGreenLagrangeStrain" << std::endl;
    }

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

double& UpdatedLagrangianAxisymmetry::CalculateIntegrationWeight(double& rIntegrationWeight)
{
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    if( dimension == 2 )
        rIntegrationWeight *= GetProperties()[THICKNESS];

    return rIntegrationWeight;
}

//************************************************************************************
//************************************************************************************

void UpdatedLagrangianAxisymmetry::EquationIdVector( EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo )
{
    int number_of_nodes = GetGeometry().size();
    int dimension = GetGeometry().WorkingSpaceDimension();
    unsigned int dimension_2 = number_of_nodes * dimension;

    if ( rResult.size() != dimension_2 )
        rResult.resize( dimension_2, false );

    for ( int i = 0; i < number_of_nodes; i++ )
    {
        int index = i * dimension;
        rResult[index] = GetGeometry()[i].GetDof( DISPLACEMENT_X ).EquationId();
        rResult[index + 1] = GetGeometry()[i].GetDof( DISPLACEMENT_Y ).EquationId();

    }

}

//************************************************************************************
//************************************************************************************

void UpdatedLagrangianAxisymmetry::GetDofList( DofsVectorType& rElementalDofList, ProcessInfo& CurrentProcessInfo )
{
    rElementalDofList.resize( 0 );

    for ( unsigned int i = 0; i < GetGeometry().size(); i++ )
    {
        rElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_X ) );
        rElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_Y ) );

        if ( GetGeometry().WorkingSpaceDimension() == 3 )
        {
            rElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_Z ) );
        }
    }
}

//************************************************************************************
//*******************DAMPING MATRIX***************************************************

void UpdatedLagrangianAxisymmetry::CalculateDampingMatrix( MatrixType& rDampingMatrix, ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    //0.-Initialize the DampingMatrix:
    unsigned int number_of_nodes = GetGeometry().size();
    unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    //resizing as needed the LHS
    unsigned int matrix_size = number_of_nodes * dimension;

    if ( rDampingMatrix.size1() != matrix_size )
        rDampingMatrix.resize( matrix_size, matrix_size, false );

    noalias( rDampingMatrix ) = ZeroMatrix( matrix_size );

    //1.-Calculate StiffnessMatrix:
    MatrixType StiffnessMatrix  = Matrix();
    this->CalculateLeftHandSide( StiffnessMatrix, rCurrentProcessInfo );

    //2.-Calculate MassMatrix:
    MatrixType MassMatrix  = Matrix();
    this->CalculateMassMatrix ( MassMatrix, rCurrentProcessInfo );

    //3.-Get Damping Coeffitients (RAYLEIGH_ALPHA, RAYLEIGH_BETA)
    double alpha = 0;
    if( GetProperties().Has(RAYLEIGH_ALPHA) )
    {
        alpha = GetProperties()[RAYLEIGH_ALPHA];
    }
    else if( rCurrentProcessInfo.Has(RAYLEIGH_ALPHA) )
    {
        alpha = rCurrentProcessInfo[RAYLEIGH_ALPHA];
    }

    double beta  = 0;
    if( GetProperties().Has(RAYLEIGH_BETA) )
    {
        beta = GetProperties()[RAYLEIGH_BETA];
    }
    else if( rCurrentProcessInfo.Has(RAYLEIGH_BETA) )
    {
        beta = rCurrentProcessInfo[RAYLEIGH_BETA];
    }

    //4.-Compose the Damping Matrix:
    //Rayleigh Damping Matrix: alpha*M + beta*K
    rDampingMatrix  = alpha * MassMatrix;
    rDampingMatrix += beta  * StiffnessMatrix;

    KRATOS_CATCH( "" )
}

//************************************************************************************
//****************MASS MATRIX*********************************************************

void UpdatedLagrangianAxisymmetry::CalculateMassMatrix( MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    // Call the values of the shape function for the single element
    GeneralVariables Variables;
    this->InitializeGeneralVariables(Variables,rCurrentProcessInfo);

    // Lumped
    unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    unsigned int matrix_size = dimension * number_of_nodes;

    if ( rMassMatrix.size1() != matrix_size )
        rMassMatrix.resize( matrix_size, matrix_size, false );

    rMassMatrix = ZeroMatrix( matrix_size );

    double TotalMass = 0;

    // TOTAL MASS OF ONE MP ELEMENT
    TotalMass = this->GetValue(MP_MASS);

    // LUMPED MATRIX
    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        double temp = Variables.N[i] * TotalMass;

        for ( unsigned int j = 0; j < dimension; j++ )
        {
            unsigned int index = i * dimension + j;
            rMassMatrix( index, index ) = temp;
        }
    }

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

// Function that return Jacobian matrix
Matrix& UpdatedLagrangianAxisymmetry::MPMJacobian( Matrix& rResult, array_1d<double,3>& rPoint)
{

    KRATOS_TRY

    // Derivatives of shape functions
    Matrix shape_functions_gradients;
    shape_functions_gradients =this->MPMShapeFunctionsLocalGradients(
                                    shape_functions_gradients);

    const GeometryType& rGeom = GetGeometry();
    unsigned int number_nodes = rGeom.PointsNumber();
    unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    if (dimension ==2)
    {
        rResult.resize( 2, 2, false );
        rResult = ZeroMatrix(2);

        for ( unsigned int i = 0; i < number_nodes; i++ )
        {
            rResult( 0, 0 ) += ( GetGeometry().GetPoint( i ).X() *  shape_functions_gradients( i, 0 ) );
            rResult( 0, 1 ) += ( GetGeometry().GetPoint( i ).X() *  shape_functions_gradients( i, 1 ) );
            rResult( 1, 0 ) += ( GetGeometry().GetPoint( i ).Y() *  shape_functions_gradients( i, 0 ) );
            rResult( 1, 1 ) += ( GetGeometry().GetPoint( i ).Y() *  shape_functions_gradients( i, 1 ) );
        }

    }
    else
    {
        KRATOS_ERROR <<  "Dimension given is wrong: Something is wrong with the given dimension in function: MPMJacobian" << std::endl;
    }

    return rResult;

    KRATOS_CATCH( "" )
}

/**
     * Jacobian in given point and given a delta position. This method calculate jacobian
     * matrix in given point and a given delta position.
     *
     * @param rPoint point which jacobians has to
* be calculated in it.
*
* @return Matrix of double which is jacobian matrix \f$ J \f$ in given point and a given delta position.
*
* @see DeterminantOfJacobian
* @see InverseOfJacobian
    */
Matrix& UpdatedLagrangianAxisymmetry::MPMJacobianDelta( Matrix& rResult, array_1d<double,3>& rPoint, Matrix & rDeltaPosition )
{
    KRATOS_TRY

    // Derivatives of shape functions
    Matrix shape_functions_gradients;
    shape_functions_gradients = this->MPMShapeFunctionsLocalGradients(
                                    shape_functions_gradients );

    const GeometryType& rGeom = GetGeometry();
    const unsigned int dimension = rGeom.WorkingSpaceDimension();

    if (dimension ==2)
    {
        rResult.resize( 2, 2, false );
        rResult = ZeroMatrix(2);

        for ( unsigned int i = 0; i < rGeom.size(); i++ )
        {
            rResult( 0, 0 ) += ( rGeom.GetPoint( i ).X() + rDeltaPosition(i,0)) * ( shape_functions_gradients( i, 0 ) );
            rResult( 0, 1 ) += ( rGeom.GetPoint( i ).X() + rDeltaPosition(i,0)) * ( shape_functions_gradients( i, 1 ) );
            rResult( 1, 0 ) += ( rGeom.GetPoint( i ).Y() + rDeltaPosition(i,1)) * ( shape_functions_gradients( i, 0 ) );
            rResult( 1, 1 ) += ( rGeom.GetPoint( i ).Y() + rDeltaPosition(i,1)) * ( shape_functions_gradients( i, 1 ) );
        }
    }
    else
    {
        KRATOS_ERROR <<  "Dimension given is wrong: Something is wrong with the given dimension in function: MPMJacobianDelta" << std::endl;
    }

    return rResult;

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

/**
     * Shape function values in given point. This method calculate the shape function
     * vector in given point.
     *
     * @param rPoint point which shape function values have to
* be calculated in it.
*
* @return Vector of double which is shape function vector \f$ N \f$ in given point.
*
    */
Vector& UpdatedLagrangianAxisymmetry::MPMShapeFunctionPointValues( Vector& rResult, array_1d<double,3>& rPoint )
{
    KRATOS_TRY

    const GeometryType& rGeom = GetGeometry();
    const unsigned int dimension = rGeom.WorkingSpaceDimension();

    array_1d<double,3> rPointLocal = ZeroVector(dimension);
    rPointLocal = rGeom.PointLocalCoordinates(rPointLocal, rPoint);

    if (dimension == 2 && rGeom.PointsNumber() == 3)
    {
        rResult.resize(3, false);

        // Shape functions
        rResult( 0 ) = 1 - rPointLocal[0] - rPointLocal[1] ;
        rResult( 1 ) = rPointLocal[0] ;
        rResult( 2 ) = rPointLocal[1];

        if(rResult(0)*rResult(1)*rResult(2) < 0.0 && rResult(0)+rResult(1)+rResult(2) > 1.0)
        {
            KRATOS_ERROR <<"Error in the evaluation of the shape functions: ELEMENT ID "<<this->Id()<<std::endl;
        }
    }
    else if (dimension == 2 && rGeom.PointsNumber() == 4)
    {
        rResult.resize(4, false);
        rResult( 0 ) = 0.25 * (1 - rPointLocal[0]) * (1 - rPointLocal[1]) ;
        rResult( 1 ) = 0.25 * (1 + rPointLocal[0]) * (1 - rPointLocal[1]) ;
        rResult( 2 ) = 0.25 * (1 + rPointLocal[0]) * (1 + rPointLocal[1]) ;
        rResult( 3 ) = 0.25 * (1 - rPointLocal[0]) * (1 + rPointLocal[1]) ;

        if(rResult(0)*rResult(1)*rResult(2)*rResult(3) < 0.0 && rResult(0)+rResult(1)+rResult(2)+rResult(3) > 1.0)
        {
            KRATOS_ERROR <<"Error in the evaluation of the shape functions: ELEMENT ID "<<this->Id()<<std::endl;
        }
    }

    return rResult;

    KRATOS_CATCH( "" )
}

Vector& UpdatedLagrangianAxisymmetry::MPMLocalCoordinates(Vector& rResult, array_1d<double,3>& rPoint)
{
    KRATOS_TRY

    // Only local coordinated of a point in a tetrahedron is computed
    rResult.resize(4,false);

    double x10 = GetGeometry()[1].Coordinates()[0]-GetGeometry()[0].Coordinates()[0];
    double x20 = GetGeometry()[2].Coordinates()[0]-GetGeometry()[0].Coordinates()[0];
    double x30 = GetGeometry()[3].Coordinates()[0]-GetGeometry()[0].Coordinates()[0];
    double y10 = GetGeometry()[1].Coordinates()[1]-GetGeometry()[0].Coordinates()[1];
    double y20 = GetGeometry()[2].Coordinates()[1]-GetGeometry()[0].Coordinates()[1];
    double y30 = GetGeometry()[3].Coordinates()[1]-GetGeometry()[0].Coordinates()[1];
    double z10 = GetGeometry()[1].Coordinates()[2]-GetGeometry()[0].Coordinates()[2];
    double z20 = GetGeometry()[2].Coordinates()[2]-GetGeometry()[0].Coordinates()[2];
    double z30 = GetGeometry()[3].Coordinates()[2]-GetGeometry()[0].Coordinates()[2];

    rResult[3] = (rPoint[0] - GetGeometry()[0].Coordinates()[0])*(y10*z20 - z10*y20) - (rPoint[1] - GetGeometry()[0].Coordinates()[1])*(x10*z20-x20*z10) + (rPoint[2] - GetGeometry()[0].Coordinates()[2])*(y20*x10 - y10*x20)/mDeterminantJ0;

    rResult[2] = (rPoint[0] - GetGeometry()[0].Coordinates()[0])*(y30*z10-y10*z30) + (rPoint[1] - GetGeometry()[0].Coordinates()[1])*(x10*z30-x30*z10) + (rPoint[2] - GetGeometry()[0].Coordinates()[2])*(y10*x30 - y30*x10)/mDeterminantJ0;

    rResult[1] = (rPoint[0] - GetGeometry()[0].Coordinates()[0])*(y20*z30-y30*z20) + (rPoint[1] - GetGeometry()[0].Coordinates()[1])*(x30*z20-x20*z30) + (rPoint[2] - GetGeometry()[0].Coordinates()[2])*(y30*x20 - x30*y20)/mDeterminantJ0;

    rResult[0] = 1 - rResult[1] - rResult[2] -rResult[3];

    return rResult;

    KRATOS_CATCH( "" )
}


// Function which return dN/de
Matrix& UpdatedLagrangianAxisymmetry::MPMShapeFunctionsLocalGradients( Matrix& rResult )
{
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    array_1d<double,3> rPointLocal = ZeroVector(dimension);

    array_1d<double,3>& xg = this->GetValue(GAUSS_COORD);
    rPointLocal = GetGeometry().PointLocalCoordinates(rPointLocal, xg);

    if (dimension == 2 && GetGeometry().PointsNumber() == 3)
    {
        rResult = ZeroMatrix( 3, 2 );
        rResult( 0, 0 ) = -1.0;
        rResult( 0, 1 ) = -1.0;
        rResult( 1, 0 ) =  1.0;
        rResult( 1, 1 ) =  0.0;
        rResult( 2, 0 ) =  0.0;
        rResult( 2, 1 ) =  1.0;
    }
    else if (dimension == 2 && GetGeometry().PointsNumber() == 4)
    {

        rResult = ZeroMatrix( 4, 2 );
        rResult( 0, 0 ) = -0.25 * (1 - rPointLocal[1]);
        rResult( 0, 1 ) = -0.25 * (1 - rPointLocal[0]);
        rResult( 1, 0 ) = 0.25 * (1 - rPointLocal[1]);
        rResult( 1, 1 ) = -0.25 * (1 + rPointLocal[0]);
        rResult( 2, 0 ) = 0.25 * (1 + rPointLocal[1]);
        rResult( 2, 1 ) = 0.25 * (1 + rPointLocal[0]);
        rResult( 3, 0 ) = -0.25 * (1 + rPointLocal[1]);
        rResult( 3, 1 ) = 0.25 * (1 - rPointLocal[0]);
    }

    return rResult;
}

//************************************************************************************
//************************************************************************************

void UpdatedLagrangianAxisymmetry::GetValuesVector( Vector& values, int Step )
{
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    unsigned int matrix_size = number_of_nodes * dimension;

    if ( values.size() != matrix_size ) values.resize( matrix_size, false );

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        unsigned int index = i * dimension;
        values[index] = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_X, Step );
        values[index + 1] = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_Y, Step );
    }
}

//************************************************************************************
//************************************************************************************

void UpdatedLagrangianAxisymmetry::GetFirstDerivativesVector( Vector& values, int Step )
{
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    unsigned int matrix_size = number_of_nodes * dimension;

    if ( values.size() != matrix_size ) values.resize( matrix_size, false );

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        unsigned int index = i * dimension;
        values[index] = GetGeometry()[i].GetSolutionStepValue( VELOCITY_X, Step );
        values[index + 1] = GetGeometry()[i].GetSolutionStepValue( VELOCITY_Y, Step );

    }
}

//************************************************************************************
//************************************************************************************

void UpdatedLagrangianAxisymmetry::GetSecondDerivativesVector( Vector& values, int Step )
{
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    unsigned int matrix_size = number_of_nodes * dimension;

    if ( values.size() != matrix_size ) values.resize( matrix_size, false );

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        unsigned int index = i * dimension;
        values[index] = GetGeometry()[i].GetSolutionStepValue( ACCELERATION_X, Step );
        values[index + 1] = GetGeometry()[i].GetSolutionStepValue( ACCELERATION_Y, Step );

        if ( dimension == 3 )
            values[index + 2] = GetGeometry()[i].GetSolutionStepValue( ACCELERATION_Z, Step );
    }
}

//************************************************************************************
//************************************************************************************
void UpdatedLagrangianAxisymmetry::GetHistoricalVariables( GeneralVariables& rVariables )
{
    //Deformation Gradient F ( set to identity )
    unsigned int size =  rVariables.F.size1();
    rVariables.detF  = 1;
    rVariables.F     = IdentityMatrix(size);

    rVariables.detF0 = mDeterminantF0;
    rVariables.F0    = mDeformationGradientF0;

}

//************************************************************************************
//************************************************************************************

void UpdatedLagrangianAxisymmetry::DecimalCorrection(Vector& rVector)
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

//************************************************************************************
//************************************************************************************

double GetPI()
{
    return std::atan(1.0)*4.0;
}


void UpdatedLagrangianAxisymmetry::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Element )
    rSerializer.save("ConstitutiveLawVector",mConstitutiveLawVector);
    rSerializer.save("DeformationGradientF0",mDeformationGradientF0);
    rSerializer.save("DeterminantF0",mDeterminantF0);
    rSerializer.save("InverseJ0",mInverseJ0);
    rSerializer.save("DeterminantJ0",mDeterminantJ0);

}

void UpdatedLagrangianAxisymmetry::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Element )
    rSerializer.load("ConstitutiveLawVector",mConstitutiveLawVector);
    rSerializer.load("DeformationGradientF0",mDeformationGradientF0);
    rSerializer.load("DeterminantF0",mDeterminantF0);
    rSerializer.load("InverseJ0",mInverseJ0);
    rSerializer.load("DeterminantJ0",mDeterminantJ0);
}

} // Namespace Kratos

