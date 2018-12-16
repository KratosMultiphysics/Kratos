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
#include "custom_elements/updated_lagrangian.hpp"
#include "utilities/math_utils.h"
#include "includes/constitutive_law.h"
#include "particle_mechanics_application_variables.h"

namespace Kratos
{

/**
 * Flags related to the element computation
 */
KRATOS_CREATE_LOCAL_FLAG( UpdatedLagrangian, COMPUTE_RHS_VECTOR,                 0 );
KRATOS_CREATE_LOCAL_FLAG( UpdatedLagrangian, COMPUTE_LHS_MATRIX,                 1 );
KRATOS_CREATE_LOCAL_FLAG( UpdatedLagrangian, COMPUTE_RHS_VECTOR_WITH_COMPONENTS, 2 );
KRATOS_CREATE_LOCAL_FLAG( UpdatedLagrangian, COMPUTE_LHS_MATRIX_WITH_COMPONENTS, 3 );

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

UpdatedLagrangian::UpdatedLagrangian( )
    : Element( )
{
    //DO NOT CALL IT: only needed for Register and Serialization!!!
}
//******************************CONSTRUCTOR*******************************************
//************************************************************************************
UpdatedLagrangian::UpdatedLagrangian( IndexType NewId, GeometryType::Pointer pGeometry )
    : Element( NewId, pGeometry )
{
    //DO NOT ADD DOFS HERE!!!
}

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

UpdatedLagrangian::UpdatedLagrangian( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
    : Element( NewId, pGeometry, pProperties )
{
    mFinalizedStep = true;


}
//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

UpdatedLagrangian::UpdatedLagrangian( UpdatedLagrangian const& rOther)
    :Element(rOther)
    ,mDeformationGradientF0(rOther.mDeformationGradientF0)
    ,mDeterminantF0(rOther.mDeterminantF0)
    ,mInverseJ0(rOther.mInverseJ0)
    ,mInverseJ(rOther.mInverseJ)
    ,mDeterminantJ0(rOther.mDeterminantJ0)
    ,mConstitutiveLawVector(rOther.mConstitutiveLawVector)
    ,mFinalizedStep(rOther.mFinalizedStep)
{
}

//*******************************ASSIGMENT OPERATOR***********************************
//************************************************************************************

UpdatedLagrangian&  UpdatedLagrangian::operator=(UpdatedLagrangian const& rOther)
{
    Element::operator=(rOther);

    mDeformationGradientF0.clear();
    mDeformationGradientF0 = rOther.mDeformationGradientF0;

    mInverseJ0.clear();
    mInverseJ0 = rOther.mInverseJ0;
    mInverseJ.clear();
    mInverseJ = rOther.mInverseJ;

    mDeterminantF0 = rOther.mDeterminantF0;
    mDeterminantJ0 = rOther.mDeterminantJ0;
    mConstitutiveLawVector = rOther.mConstitutiveLawVector;

    return *this;
}

//*********************************OPERATIONS*****************************************
//************************************************************************************

Element::Pointer UpdatedLagrangian::Create( IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties ) const
{
    return Element::Pointer( new UpdatedLagrangian( NewId, GetGeometry().Create( ThisNodes ), pProperties ) );
}

//************************************CLONE*******************************************
//************************************************************************************

Element::Pointer UpdatedLagrangian::Clone( IndexType NewId, NodesArrayType const& rThisNodes ) const
{
    UpdatedLagrangian NewElement (NewId, GetGeometry().Create( rThisNodes ), pGetProperties() );

    NewElement.mConstitutiveLawVector = mConstitutiveLawVector->Clone();

    NewElement.mDeformationGradientF0 = mDeformationGradientF0;

    NewElement.mInverseJ0 = mInverseJ0;
    NewElement.mInverseJ = mInverseJ;

    NewElement.mDeterminantF0 = mDeterminantF0;
    NewElement.mDeterminantJ0 = mDeterminantJ0;

    return Element::Pointer( new UpdatedLagrangian(NewElement) );
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************
UpdatedLagrangian::~UpdatedLagrangian()
{
}


//************************************************************************************
//************************************************************************************

void UpdatedLagrangian::Initialize()
{
    KRATOS_TRY

    // Initial position of the particle
    const array_1d<double,3>& xg = this->GetValue(MP_COORD);

    // Initialize parameters
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    mDeterminantF0 = 1;
    mDeformationGradientF0 = IdentityMatrix(dimension);

    // Compute initial jacobian matrix and inverses
    Matrix J0 = ZeroMatrix(dimension, dimension);
    J0 = this->MPMJacobian(J0, xg);
    MathUtils<double>::InvertMatrix( J0, mInverseJ0, mDeterminantJ0 );

    // Compute current jacobian matrix and inverses
    Matrix j = ZeroMatrix(dimension, dimension);
    j = this->MPMJacobian(j,xg);
    double detj;
    MathUtils<double>::InvertMatrix( j, mInverseJ, detj );

    // Initialize constitutive law and materials
    InitializeMaterial();

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void UpdatedLagrangian::InitializeGeneralVariables (GeneralVariables& rVariables, const ProcessInfo& rCurrentProcessInfo)
{
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    unsigned int voigtsize  = 3;

    if( dimension == 3 )
    {
        voigtsize  = 6;
    }
    rVariables.detF  = 1;

    rVariables.detF0 = 1;

    rVariables.detFT = 1;

    rVariables.detJ = 1;

    rVariables.B.resize( voigtsize, number_of_nodes * dimension, false );

    rVariables.F.resize( dimension, dimension, false );

    rVariables.F0.resize( dimension, dimension, false );

    rVariables.FT.resize( dimension, dimension, false );

    rVariables.ConstitutiveMatrix.resize( voigtsize, voigtsize, false );

    rVariables.StrainVector.resize( voigtsize, false );

    rVariables.StressVector.resize( voigtsize, false );

    rVariables.DN_DX.resize( number_of_nodes, dimension, false );
    rVariables.DN_De.resize( number_of_nodes, dimension, false );

    const array_1d<double,3>& xg = this->GetValue(MP_COORD);

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

void UpdatedLagrangian::SetGeneralVariables(GeneralVariables& rVariables,
        ConstitutiveLaw::Parameters& rValues)
{
    GeometryType& rGeom = GetGeometry();

    // Variables.detF is the determinant of the incremental total deformation gradient
    rVariables.detF  = MathUtils<double>::Det(rVariables.F);

    // Check if detF is negative (element is inverted)
    if(rVariables.detF<0)
    {
        KRATOS_INFO("UpdatedLagrangian")<<" Element: "<<this->Id()<<std::endl;
        KRATOS_INFO("UpdatedLagrangian")<<" Element position: "<<this->GetValue(MP_COORD)<<std::endl;
        const unsigned int number_of_nodes = rGeom.PointsNumber();

        for ( unsigned int i = 0; i < number_of_nodes; i++ )
        {
            const array_1d<double, 3> & current_position      = rGeom[i].Coordinates();
            const array_1d<double, 3> & current_displacement  = rGeom[i].FastGetSolutionStepValue(DISPLACEMENT);
            const array_1d<double, 3> & previous_displacement = rGeom[i].FastGetSolutionStepValue(DISPLACEMENT,1);

            KRATOS_INFO("UpdatedLagrangian")<<" NODE ["<<rGeom[i].Id()<<"]: (Current position: "<<current_position<<") "<<std::endl;
            KRATOS_INFO("UpdatedLagrangian")<<" ---Current Disp: "<<current_displacement<<" (Previour Disp: "<<previous_displacement<<")"<<std::endl;
        }

        for ( unsigned int i = 0; i < number_of_nodes; i++ )
        {
            if( rGeom[i].SolutionStepsDataHas(CONTACT_FORCE) )
            {
                const array_1d<double, 3 > & PreContactForce = rGeom[i].FastGetSolutionStepValue(CONTACT_FORCE,1);
                const array_1d<double, 3 > & ContactForce = rGeom[i].FastGetSolutionStepValue(CONTACT_FORCE);
                KRATOS_INFO("UpdatedLagrangian")<<" ---Contact_Force: (Pre:"<<PreContactForce<<", Current:"<<ContactForce<<") "<<std::endl;
            }
            else
            {
                KRATOS_INFO("UpdatedLagrangian")<<" ---Contact_Force: NULL "<<std::endl;
            }
        }

        KRATOS_ERROR << "MPM UPDATED LAGRANGIAN DISPLACEMENT ELEMENT INVERTED: |F|<0  detF = " << rVariables.detF << std::endl;
    }

    rVariables.detFT = rVariables.detF * rVariables.detF0;
    rVariables.FT    = prod( rVariables.F, rVariables.F0 );

    rValues.SetDeterminantF(rVariables.detFT);
    rValues.SetDeformationGradientF(rVariables.FT);
    rValues.SetStrainVector(rVariables.StrainVector);
    rValues.SetStressVector(rVariables.StressVector);
    rValues.SetConstitutiveMatrix(rVariables.ConstitutiveMatrix);
    rValues.SetShapeFunctionsDerivatives(rVariables.DN_DX);
    rValues.SetShapeFunctionsValues(rVariables.N);

}

//************************************************************************************
//*****************check size of LHS and RHS matrices*********************************

void UpdatedLagrangian::InitializeSystemMatrices(MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        Flags& rCalculationFlags)
{
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

    // Resizing the LHS matrix if needed
    unsigned int matrix_size = number_of_nodes * dimension;   //number of degrees of freedom

    if ( rCalculationFlags.Is(UpdatedLagrangian::COMPUTE_LHS_MATRIX) ) //calculation of the matrix is required
    {
        if ( rLeftHandSideMatrix.size1() != matrix_size )
            rLeftHandSideMatrix.resize( matrix_size, matrix_size, false );

        noalias( rLeftHandSideMatrix ) = ZeroMatrix(matrix_size, matrix_size); //resetting LHS
    }

    // Resizing the RHS vector if needed
    if ( rCalculationFlags.Is(UpdatedLagrangian::COMPUTE_RHS_VECTOR) ) //calculation of the matrix is required
    {
        if ( rRightHandSideVector.size() != matrix_size )
            rRightHandSideVector.resize( matrix_size, false );

        rRightHandSideVector = ZeroVector( matrix_size ); //resetting RHS
    }
}

//************************************************************************************
//************************************************************************************

void UpdatedLagrangian::CalculateElementalSystem( LocalSystemComponents& rLocalSystem,
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
    const double MP_Density = (GetProperties()[DENSITY]) / Variables.detFT;
    this->SetValue(MP_DENSITY, MP_Density);

    // The MP_Volume (integration weight) is evaluated
    const double MP_Volume = this->GetValue(MP_MASS)/this->GetValue(MP_DENSITY);
    this->SetValue(MP_VOLUME, MP_Volume);

    if ( rLocalSystem.CalculationFlags.Is(UpdatedLagrangian::COMPUTE_LHS_MATRIX) ) // if calculation of the matrix is required
    {
        // Contributions to stiffness matrix calculated on the reference configuration
        this->CalculateAndAddLHS ( rLocalSystem, Variables, MP_Volume );
    }

    if ( rLocalSystem.CalculationFlags.Is(UpdatedLagrangian::COMPUTE_RHS_VECTOR) ) // if calculation of the vector is required
    {
        // Contribution to forces (in residual term) are calculated
        volume_force  = this->CalculateVolumeForce( volume_force, Variables );
        this->CalculateAndAddRHS ( rLocalSystem, Variables, volume_force, MP_Volume );
    }

    KRATOS_CATCH( "" )
}
//*********************************COMPUTE KINEMATICS*********************************
//************************************************************************************


void UpdatedLagrangian::CalculateKinematics(GeneralVariables& rVariables, ProcessInfo& rCurrentProcessInfo)

{
    KRATOS_TRY

    // Define the stress measure
    rVariables.StressMeasure = ConstitutiveLaw::StressMeasure_Cauchy;

    // Calculating the inverse of the jacobian and the parameters needed [d£/dx_n]
    Matrix InvJ;
    MathUtils<double>::InvertMatrix( rVariables.J, InvJ, rVariables.detJ);

    // Calculating the inverse of the jacobian and the parameters needed [d£/(dx_n+1)]
    Matrix Invj;
    MathUtils<double>::InvertMatrix( rVariables.j, Invj, rVariables.detJ ); //overwrites detJ

    // Compute cartesian derivatives [dN/dx_n+1]
    rVariables.DN_DX = prod( rVariables.DN_De, Invj); //overwrites DX now is the current position dx

    /* NOTE::
    Deformation Gradient F [(dx_n+1 - dx_n)/dx_n] is to be updated in constitutive law parameter as total deformation gradient.
    The increment of total deformation gradient can be evaluated in 2 ways, which are:
    1. By: noalias( rVariables.F ) = prod( rVariables.j, InvJ);
    2. By means of the gradient of nodal displacement: using this second expression quadratic convergence is not guarantee

    (NOTICE: Here, we are using method no. 2)*/

    // METHOD 1: Update Deformation gradient: F [dx_n+1/dx_n] = [dx_n+1/d£] [d£/dx_n]
    // noalias( rVariables.F ) = prod( rVariables.j, InvJ);

    // METHOD 2: Update Deformation gradient: F_ij = δ_ij + u_i,j
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    Matrix I = IdentityMatrix(dimension);
    Matrix gradient_displacement = ZeroMatrix(dimension, dimension);
    rVariables.CurrentDisp = CalculateCurrentDisp(rVariables.CurrentDisp, rCurrentProcessInfo);
    gradient_displacement = prod(trans(rVariables.CurrentDisp),rVariables.DN_DX);

    noalias( rVariables.F ) = (I + gradient_displacement);

    // Determinant of the previous Deformation Gradient F_n
    rVariables.detF0 = mDeterminantF0;
    rVariables.F0    = mDeformationGradientF0;

    // Compute the deformation matrix B
    this->CalculateDeformationMatrix(rVariables.B, rVariables.F, rVariables.DN_DX);

    KRATOS_CATCH( "" )
}
//************************************************************************************

void UpdatedLagrangian::CalculateDeformationMatrix(Matrix& rB,
        Matrix& rF,
        Matrix& rDN_DX)
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

    rB.clear(); // Set all components to zero

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
        KRATOS_ERROR <<  "Dimension given is wrong: Something is wrong with the given dimension in function: CalculateDeformationMatrix" << std::endl;
    }

    KRATOS_CATCH( "" )
}
//************************************************************************************
//************************************************************************************

void UpdatedLagrangian::CalculateAndAddRHS(LocalSystemComponents& rLocalSystem, GeneralVariables& rVariables, Vector& rVolumeForce, const double& rIntegrationWeight)
{
    // Contribution of the internal and external forces
    if( rLocalSystem.CalculationFlags.Is( UpdatedLagrangian::COMPUTE_RHS_VECTOR_WITH_COMPONENTS ) )
    {
        std::vector<VectorType>& rRightHandSideVectors = rLocalSystem.GetRightHandSideVectors();
        const std::vector< Variable< VectorType > >& rRightHandSideVariables = rLocalSystem.GetRightHandSideVariables();
        for( unsigned int i=0; i<rRightHandSideVariables.size(); i++ )
        {
            bool calculated = false;
            if( rRightHandSideVariables[i] == EXTERNAL_FORCES_VECTOR )
            {
                // Operation performed: rRightHandSideVector += ExtForce*IntToReferenceWeight
                this->CalculateAndAddExternalForces( rRightHandSideVectors[i], rVariables, rVolumeForce, rIntegrationWeight );
                calculated = true;
            }

            if( rRightHandSideVariables[i] == INTERNAL_FORCES_VECTOR )
            {
                // Operation performed: rRightHandSideVector -= IntForce*IntToReferenceWeight
                this->CalculateAndAddInternalForces( rRightHandSideVectors[i], rVariables, rIntegrationWeight );
                calculated = true;
            }

            KRATOS_ERROR_IF(calculated == false) << " ELEMENT can not supply the required local system variable: " << rRightHandSideVariables[i] << std::endl;
        }
    }
    else
    {
        VectorType& rRightHandSideVector = rLocalSystem.GetRightHandSideVector();

        // Operation performed: rRightHandSideVector += ExtForce*IntToReferenceWeight
        this->CalculateAndAddExternalForces( rRightHandSideVector, rVariables, rVolumeForce, rIntegrationWeight );

        // Operation performed: rRightHandSideVector -= IntForce*IntToReferenceWeight
        this->CalculateAndAddInternalForces( rRightHandSideVector, rVariables, rIntegrationWeight );
    }
}

//************************************************************************************
//*********************Calculate the contribution of external force*******************

void UpdatedLagrangian::CalculateAndAddExternalForces(VectorType& rRightHandSideVector,
        GeneralVariables& rVariables,
        Vector& rVolumeForce,
        const double& rIntegrationWeight)
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        int index = dimension * i;

        for ( unsigned int j = 0; j < dimension; j++ )
        {
            rRightHandSideVector[index + j] += rVariables.N[i] * rVolumeForce[j];
        }
    }

    KRATOS_CATCH( "" )
}
//************************************************************************************
//************************************************************************************

void UpdatedLagrangian::CalculateAndAddInternalForces(VectorType& rRightHandSideVector,
        GeneralVariables & rVariables,
        const double& rIntegrationWeight)
{
    KRATOS_TRY

    VectorType internal_forces = rIntegrationWeight * prod( trans( rVariables.B ), rVariables.StressVector );
    noalias( rRightHandSideVector ) -= internal_forces;

    KRATOS_CATCH( "" )
}
//************************************************************************************
//************************************************************************************

void UpdatedLagrangian::CalculateAndAddLHS(LocalSystemComponents& rLocalSystem, GeneralVariables& rVariables, const double& rIntegrationWeight)
{
    // Contributions of the stiffness matrix calculated on the reference configuration
    if( rLocalSystem.CalculationFlags.Is( UpdatedLagrangian::COMPUTE_LHS_MATRIX_WITH_COMPONENTS ) )
    {
        std::vector<MatrixType>& rLeftHandSideMatrices = rLocalSystem.GetLeftHandSideMatrices();
        const std::vector< Variable< MatrixType > >& rLeftHandSideVariables = rLocalSystem.GetLeftHandSideVariables();

        for( unsigned int i=0; i<rLeftHandSideVariables.size(); i++ )
        {
            bool calculated = false;
            if( rLeftHandSideVariables[i] == MATERIAL_STIFFNESS_MATRIX )
            {
                // Operation performed: add K_material to the rLefsHandSideMatrix
                this->CalculateAndAddKuum( rLeftHandSideMatrices[i], rVariables, rIntegrationWeight );
                calculated = true;
            }

            if( rLeftHandSideVariables[i] == GEOMETRIC_STIFFNESS_MATRIX )
            {
                // Operation performed: add K_geometry to the rLefsHandSideMatrix
                this->CalculateAndAddKuug( rLeftHandSideMatrices[i], rVariables, rIntegrationWeight );
                calculated = true;
            }

            KRATOS_ERROR_IF(calculated == false) <<  " ELEMENT can not supply the required local system variable: " << rLeftHandSideVariables[i] << std::endl;
        }
    }
    else
    {
        MatrixType& rLeftHandSideMatrix = rLocalSystem.GetLeftHandSideMatrix();

        // Operation performed: add K_material to the rLefsHandSideMatrix
        this->CalculateAndAddKuum( rLeftHandSideMatrix, rVariables, rIntegrationWeight );

        // Operation performed: add K_geometry to the rLefsHandSideMatrix
        this->CalculateAndAddKuug( rLeftHandSideMatrix, rVariables, rIntegrationWeight );
    }
}
//************************************************************************************
//************************************************************************************

void UpdatedLagrangian::CalculateAndAddKuum(MatrixType& rLeftHandSideMatrix,
        GeneralVariables& rVariables,
        const double& rIntegrationWeight)
{
    KRATOS_TRY

    noalias( rLeftHandSideMatrix ) += prod( trans( rVariables.B ),  rIntegrationWeight * Matrix( prod( rVariables.ConstitutiveMatrix, rVariables.B ) ) );

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void UpdatedLagrangian::CalculateAndAddKuug(MatrixType& rLeftHandSideMatrix,
        GeneralVariables& rVariables,
        const double& rIntegrationWeight)
{
    KRATOS_TRY

    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    Matrix stress_tensor = MathUtils<double>::StressVectorToTensor( rVariables.StressVector );
    Matrix reduced_Kg = prod( rVariables.DN_DX, rIntegrationWeight * Matrix( prod( stress_tensor, trans( rVariables.DN_DX ) ) ) );
    MathUtils<double>::ExpandAndAddReducedMatrix( rLeftHandSideMatrix, reduced_Kg, dimension );

    KRATOS_CATCH( "" )
}

//************************************CALCULATE VOLUME CHANGE*************************
//************************************************************************************

double& UpdatedLagrangian::CalculateVolumeChange( double& rVolumeChange, GeneralVariables& rVariables )
{
    KRATOS_TRY

    rVolumeChange = 1.0 / (rVariables.detF * rVariables.detF0);

    return rVolumeChange;

    KRATOS_CATCH( "" )
}

//************************************CALCULATE VOLUME ACCELERATION*******************
//************************************************************************************

Vector& UpdatedLagrangian::CalculateVolumeForce( Vector& rVolumeForce, GeneralVariables& rVariables )
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
void UpdatedLagrangian::CalculateRightHandSide( VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo )
{
    // Create local system components
    LocalSystemComponents LocalSystem;

    // Set calculation flags
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


void UpdatedLagrangian::CalculateRightHandSide( std::vector< VectorType >& rRightHandSideVectors, const std::vector< Variable< VectorType > >& rRHSVariables, ProcessInfo& rCurrentProcessInfo )
{
    // Create local system components
    LocalSystemComponents LocalSystem;

    // Set calculation flags
    LocalSystem.CalculationFlags.Set(UpdatedLagrangian::COMPUTE_RHS_VECTOR);
    LocalSystem.CalculationFlags.Set(UpdatedLagrangian::COMPUTE_RHS_VECTOR_WITH_COMPONENTS);

    MatrixType LeftHandSideMatrix = Matrix();

    // Initialize sizes for the system components:
    if( rRHSVariables.size() != rRightHandSideVectors.size() )
        rRightHandSideVectors.resize(rRHSVariables.size());

    for( unsigned int i=0; i<rRightHandSideVectors.size(); i++ )
    {
        this->InitializeSystemMatrices( LeftHandSideMatrix, rRightHandSideVectors[i], LocalSystem.CalculationFlags );
    }

    // Set Variables to Local system components
    LocalSystem.SetLeftHandSideMatrix(LeftHandSideMatrix);
    LocalSystem.SetRightHandSideVectors(rRightHandSideVectors);

    LocalSystem.SetRightHandSideVariables(rRHSVariables);

    // Calculate elemental system
    CalculateElementalSystem( LocalSystem, rCurrentProcessInfo );
}


//************************************************************************************
//************************************************************************************


void UpdatedLagrangian::CalculateLeftHandSide( MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo )
{
    // Create local system components
    LocalSystemComponents LocalSystem;

    // Set calculation flags
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


void UpdatedLagrangian::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo )
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

void UpdatedLagrangian::CalculateLocalSystem( std::vector< MatrixType >& rLeftHandSideMatrices,
        const std::vector< Variable< MatrixType > >& rLHSVariables,
        std::vector< VectorType >& rRightHandSideVectors,
        const std::vector< Variable< VectorType > >& rRHSVariables,
        ProcessInfo& rCurrentProcessInfo )
{
    // Create local system components
    LocalSystemComponents LocalSystem;

    // Set calculation flags
    LocalSystem.CalculationFlags.Set(UpdatedLagrangian::COMPUTE_LHS_MATRIX_WITH_COMPONENTS);
    LocalSystem.CalculationFlags.Set(UpdatedLagrangian::COMPUTE_RHS_VECTOR_WITH_COMPONENTS);

    // Initialize sizes for the system components:
    if( rLHSVariables.size() != rLeftHandSideMatrices.size() )
        rLeftHandSideMatrices.resize(rLHSVariables.size());

    if( rRHSVariables.size() != rRightHandSideVectors.size() )
        rRightHandSideVectors.resize(rRHSVariables.size());

    LocalSystem.CalculationFlags.Set(UpdatedLagrangian::COMPUTE_LHS_MATRIX);
    for( unsigned int i=0; i<rLeftHandSideMatrices.size(); i++ )
    {
        // Note: rRightHandSideVectors.size() > 0
        this->InitializeSystemMatrices( rLeftHandSideMatrices[i], rRightHandSideVectors[0], LocalSystem.CalculationFlags );
    }

    LocalSystem.CalculationFlags.Set(UpdatedLagrangian::COMPUTE_RHS_VECTOR);
    LocalSystem.CalculationFlags.Set(UpdatedLagrangian::COMPUTE_LHS_MATRIX,false);

    for( unsigned int i=0; i<rRightHandSideVectors.size(); i++ )
    {
        // Note: rLeftHandSideMatrices.size() > 0
        this->InitializeSystemMatrices( rLeftHandSideMatrices[0], rRightHandSideVectors[i], LocalSystem.CalculationFlags );
    }
    LocalSystem.CalculationFlags.Set(UpdatedLagrangian::COMPUTE_LHS_MATRIX,true);

    // Set Variables to Local system components
    LocalSystem.SetLeftHandSideMatrices(rLeftHandSideMatrices);
    LocalSystem.SetRightHandSideVectors(rRightHandSideVectors);

    LocalSystem.SetLeftHandSideVariables(rLHSVariables);
    LocalSystem.SetRightHandSideVariables(rRHSVariables);

    // Calculate elemental system
    CalculateElementalSystem( LocalSystem, rCurrentProcessInfo );
}


////***********************************************************************************
////***********************************************************************************
void UpdatedLagrangian::Calculate(const Variable<double>& rVariable,
                                  double& Output,
                                  const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    if (rVariable == DENSITY)
    {
        GeometryType& rGeom = GetGeometry();
        const unsigned int dimension = rGeom.WorkingSpaceDimension();
        const unsigned int number_of_nodes = rGeom.PointsNumber();
        const array_1d<double,3>& xg = this->GetValue(MP_COORD);
        GeneralVariables Variables;
        Matrix J0 = ZeroMatrix(dimension, dimension);

        J0 = this->MPMJacobian(J0, xg);

        //calculating and storing inverse and the determinant of the jacobian
        MathUtils<double>::InvertMatrix( J0, mInverseJ0, mDeterminantJ0 );

        Variables.N = this->MPMShapeFunctionPointValues(Variables.N, xg);
        const double & MP_Mass = this->GetValue(MP_MASS);

        for (unsigned int i=0; i<number_of_nodes; i++)
        {
            rGeom[i].SetLock();
            rGeom[i].FastGetSolutionStepValue(AUX_R) += Variables.N[i] * (MP_Mass);
            rGeom[i].UnSetLock();
        }
    }

    KRATOS_CATCH( "" )
}


void UpdatedLagrangian::Calculate(const Variable<array_1d<double, 3 > >& rVariable,
                                  array_1d<double, 3 > & Output,
                                  const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    if(rVariable == VELOCITY)
    {
        GeometryType& rGeom = GetGeometry();
        const unsigned int dimension = rGeom.WorkingSpaceDimension();
        const unsigned int number_of_nodes = rGeom.PointsNumber();
        const array_1d<double,3>& xg = this->GetValue(MP_COORD);
        GeneralVariables Variables;
        Matrix J0 = ZeroMatrix(dimension, dimension);

        J0 = this->MPMJacobian(J0, xg);

        // Calculating and storing inverse and the determinant of the jacobian
        MathUtils<double>::InvertMatrix( J0, mInverseJ0, mDeterminantJ0 );

        Variables.N = this->MPMShapeFunctionPointValues(Variables.N, xg);
        const array_1d<double,3>& MP_Velocity = this->GetValue(MP_VELOCITY);
        const double & MP_Mass = this->GetValue(MP_MASS);

        array_1d<double,3> NodalAuxRVel;
        for (unsigned int i=0; i<number_of_nodes; i++)

        {
            for (unsigned int j = 0; j < dimension; j++)
            {
                NodalAuxRVel[j] = Variables.N[i] * MP_Mass * MP_Velocity[j];
            }
            rGeom[i].SetLock();
            rGeom[i].FastGetSolutionStepValue(AUX_R_VEL) += NodalAuxRVel;
            rGeom[i].UnSetLock();
        }
    }

    if(rVariable == ACCELERATION)
    {
        GeometryType& rGeom = GetGeometry();
        const unsigned int dimension = rGeom.WorkingSpaceDimension();
        const unsigned int number_of_nodes = rGeom.PointsNumber();
        const array_1d<double,3>& xg = this->GetValue(MP_COORD);
        GeneralVariables Variables;
        Matrix J0 = ZeroMatrix(dimension, dimension);

        J0 = this->MPMJacobian(J0, xg);

        //calculating and storing inverse and the determinant of the jacobian
        MathUtils<double>::InvertMatrix( J0, mInverseJ0, mDeterminantJ0 );

        Variables.N = this->MPMShapeFunctionPointValues(Variables.N, xg);
        const array_1d<double,3>& MP_Acceleration = this->GetValue(MP_ACCELERATION);
        const double & MP_Mass = this->GetValue(MP_MASS);

        array_1d<double,3> NodalAuxRAcc;
        for (unsigned int i=0; i<number_of_nodes; i++)

        {
            for (unsigned int j = 0; j < dimension; j++)
            {
                NodalAuxRAcc[j] = Variables.N[i] * MP_Mass * MP_Acceleration[j];
            }
            rGeom[i].SetLock();
            rGeom[i].FastGetSolutionStepValue(AUX_R_ACC) += NodalAuxRAcc;
            rGeom[i].UnSetLock();
        }
    }

    KRATOS_CATCH( "" )
}
//*******************************************************************************************


void UpdatedLagrangian::InitializeSolutionStep( ProcessInfo& rCurrentProcessInfo )
{
    /* NOTE:
    In the InitializeSolutionStep of each time step the nodal initial conditions are evaluated.
    This function is called by the base scheme class.*/
    GeometryType& rGeom = GetGeometry();
    const unsigned int dimension = rGeom.WorkingSpaceDimension();
    const unsigned int number_of_nodes = rGeom.PointsNumber();
    const array_1d<double,3> & xg = this->GetValue(MP_COORD);
    GeneralVariables Variables;

    // Calculating and storing inverse and the determinant of the jacobian
    Matrix J0 = ZeroMatrix(dimension, dimension);
    J0 = this->MPMJacobian(J0, xg);
    MathUtils<double>::InvertMatrix( J0, mInverseJ0, mDeterminantJ0 );

    // Calculating shape function
    Variables.N = this->MPMShapeFunctionPointValues(Variables.N, xg);

    mFinalizedStep = false;

    const array_1d<double,3>& MP_Velocity = this->GetValue(MP_VELOCITY);
    const array_1d<double,3>& MP_Acceleration = this->GetValue(MP_ACCELERATION);
    array_1d<double,3>& AUX_MP_Velocity = this->GetValue(AUX_MP_VELOCITY);
    array_1d<double,3>& AUX_MP_Acceleration = this->GetValue(AUX_MP_ACCELERATION);
    const double & MP_Mass = this->GetValue(MP_MASS);

    array_1d<double,3> nodal_momentum = ZeroVector(3);
    array_1d<double,3> nodal_inertia  = ZeroVector(3);

    for (unsigned int j=0; j<number_of_nodes; j++)
    {
        // These are the values of nodal velocity and nodal acceleration evaluated in the initialize solution step
        const array_1d<double, 3 > & nodal_acceleration = rGeom[j].FastGetSolutionStepValue(ACCELERATION,1);
        const array_1d<double, 3 > & nodal_velocity = rGeom[j].FastGetSolutionStepValue(VELOCITY,1);

        for (unsigned int k = 0; k < dimension; k++)
        {
            AUX_MP_Velocity[k]     += Variables.N[j] * nodal_velocity[k];
            AUX_MP_Acceleration[k] += Variables.N[j] * nodal_acceleration[k];
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

        rGeom[i].SetLock();
        rGeom[i].FastGetSolutionStepValue(NODAL_MOMENTUM, 0) += nodal_momentum;
        rGeom[i].FastGetSolutionStepValue(NODAL_INERTIA, 0)  += nodal_inertia;

        rGeom[i].FastGetSolutionStepValue(NODAL_MASS, 0) += Variables.N[i] * MP_Mass;
        rGeom[i].UnSetLock();

    }

    AUX_MP_Velocity.clear();
    AUX_MP_Acceleration.clear();

}

////************************************************************************************
////************************************************************************************
void UpdatedLagrangian::InitializeNonLinearIteration( ProcessInfo& rCurrentProcessInfo )
{

}

////************************************************************************************
////************************************************************************************

void UpdatedLagrangian::FinalizeNonLinearIteration( ProcessInfo& rCurrentProcessInfo )
{

}

////************************************************************************************
////************************************************************************************

void UpdatedLagrangian::FinalizeSolutionStep( ProcessInfo& rCurrentProcessInfo )
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

    // Compute element kinematics B, F, DN_DX ...
    this->CalculateKinematics(Variables, rCurrentProcessInfo);

    // Set general variables to constitutivelaw parameters
    this->SetGeneralVariables(Variables,Values);

    // Call the constitutive law to update material variables
    mConstitutiveLawVector->FinalizeMaterialResponse(Values, Variables.StressMeasure);

    // Call the element internal variables update
    this->FinalizeStepVariables(Variables, rCurrentProcessInfo);

    mFinalizedStep = true;

    KRATOS_CATCH( "" )
}


////************************************************************************************
////************************************************************************************

void UpdatedLagrangian::FinalizeStepVariables( GeneralVariables & rVariables, const ProcessInfo& rCurrentProcessInfo)
{
    // Update internal (historical) variables
    mDeterminantF0         = rVariables.detF* rVariables.detF0;
    mDeformationGradientF0 = prod(rVariables.F, rVariables.F0);

    this->SetValue(MP_CAUCHY_STRESS_VECTOR, rVariables.StressVector);
    this->SetValue(MP_ALMANSI_STRAIN_VECTOR, rVariables.StrainVector);

    // Delta Plastic Strains
    double DeltaPlasticStrain = mConstitutiveLawVector->GetValue(MP_DELTA_PLASTIC_STRAIN, DeltaPlasticStrain );
    this->SetValue(MP_DELTA_PLASTIC_STRAIN, DeltaPlasticStrain);

    double DeltaPlasticVolumetricStrain = mConstitutiveLawVector->GetValue(MP_DELTA_PLASTIC_VOLUMETRIC_STRAIN, DeltaPlasticVolumetricStrain);
    this->SetValue(MP_DELTA_PLASTIC_VOLUMETRIC_STRAIN, DeltaPlasticVolumetricStrain);

    double DeltaPlasticDeviatoricStrain = mConstitutiveLawVector->GetValue(MP_DELTA_PLASTIC_DEVIATORIC_STRAIN, DeltaPlasticDeviatoricStrain);
    this->SetValue(MP_DELTA_PLASTIC_DEVIATORIC_STRAIN, DeltaPlasticDeviatoricStrain);

    // Total Plastic Strain
    double EquivalentPlasticStrain = mConstitutiveLawVector->GetValue(MP_EQUIVALENT_PLASTIC_STRAIN, EquivalentPlasticStrain );
    this->SetValue(MP_EQUIVALENT_PLASTIC_STRAIN, EquivalentPlasticStrain);

    double AccumulatedPlasticVolumetricStrain = mConstitutiveLawVector->GetValue(MP_ACCUMULATED_PLASTIC_VOLUMETRIC_STRAIN, AccumulatedPlasticVolumetricStrain);
    this->SetValue(MP_ACCUMULATED_PLASTIC_VOLUMETRIC_STRAIN, AccumulatedPlasticVolumetricStrain);

    double AccumulatedPlasticDeviatoricStrain = mConstitutiveLawVector->GetValue(MP_ACCUMULATED_PLASTIC_DEVIATORIC_STRAIN, AccumulatedPlasticDeviatoricStrain);
    this->SetValue(MP_ACCUMULATED_PLASTIC_DEVIATORIC_STRAIN, AccumulatedPlasticDeviatoricStrain);

    MathUtils<double>::InvertMatrix( rVariables.j, mInverseJ, rVariables.detJ );

    this->UpdateGaussPoint(rVariables, rCurrentProcessInfo);

}

//************************************************************************************
//************************************************************************************
/**
 * The position of the Gauss points/Material points is updated
 */

void UpdatedLagrangian::UpdateGaussPoint( GeneralVariables & rVariables, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    rVariables.CurrentDisp = CalculateCurrentDisp(rVariables.CurrentDisp, rCurrentProcessInfo);
    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    const array_1d<double,3> & xg = this->GetValue(MP_COORD);
    const array_1d<double,3> & MP_PreviousAcceleration = this->GetValue(MP_ACCELERATION);
    const array_1d<double,3> & MP_PreviousVelocity = this->GetValue(MP_VELOCITY);

    array_1d<double,3> delta_xg = ZeroVector(3);
    array_1d<double,3> MP_Acceleration = ZeroVector(3);
    array_1d<double,3> MP_Velocity = ZeroVector(3);
    const double delta_time = rCurrentProcessInfo[DELTA_TIME];

    rVariables.N = this->MPMShapeFunctionPointValues(rVariables.N, xg);

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        if (rVariables.N[i] > 1e-16)
        {
            const array_1d<double, 3 > & nodal_acceleration = GetGeometry()[i].FastGetSolutionStepValue(ACCELERATION);

            for ( unsigned int j = 0; j < dimension; j++ )
            {
                delta_xg[j] += rVariables.N[i] * rVariables.CurrentDisp(i,j);
                MP_Acceleration[j] += rVariables.N[i] * nodal_acceleration[j];

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
    this -> SetValue(MP_COORD,new_xg);

    // Update the MP Acceleration
    this -> SetValue(MP_ACCELERATION,MP_Acceleration);

    // Update the MP total displacement
    array_1d<double,3>& MP_Displacement = this->GetValue(MP_DISPLACEMENT);
    MP_Displacement += delta_xg;
    this -> SetValue(MP_DISPLACEMENT,MP_Displacement);

    KRATOS_CATCH( "" )
}


void UpdatedLagrangian::InitializeMaterial()
{
    KRATOS_TRY
    const array_1d<double,3>& xg = this->GetValue(MP_COORD);
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

void UpdatedLagrangian::ResetConstitutiveLaw()
{
    KRATOS_TRY
    const array_1d<double,3>& xg = this->GetValue(MP_COORD);
    GeneralVariables Variables;

    if ( GetProperties()[CONSTITUTIVE_LAW] != NULL )
    {
        mConstitutiveLawVector->ResetMaterial( GetProperties(), GetGeometry(), this->MPMShapeFunctionPointValues(Variables.N, xg) );
    }

    KRATOS_CATCH( "" )
}


//*************************COMPUTE CURRENT DISPLACEMENT*******************************
//************************************************************************************
/*
This function convert the computed nodal displacement into matrix of (number_of_nodes, dimension)
*/
Matrix& UpdatedLagrangian::CalculateCurrentDisp(Matrix & rCurrentDisp, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    GeometryType& rGeom = GetGeometry();
    const unsigned int number_of_nodes = rGeom.PointsNumber();
    const unsigned int dimension = rGeom.WorkingSpaceDimension();

    rCurrentDisp = ZeroMatrix(number_of_nodes, dimension);

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        const array_1d<double, 3 > & current_displacement  = rGeom[i].FastGetSolutionStepValue(DISPLACEMENT);

        for ( unsigned int j = 0; j < dimension; j++ )
        {
            rCurrentDisp(i,j) = current_displacement[j];
        }
    }

    return rCurrentDisp;

    KRATOS_CATCH( "" )
}


//*************************COMPUTE ALMANSI STRAIN*************************************
//************************************************************************************
// Almansi Strain: E = 0.5 (I - U^(-2))
void UpdatedLagrangian::CalculateAlmansiStrain(const Matrix& rF,
        Vector& rStrainVector )
{
    KRATOS_TRY

    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    // Left Cauchy-Green Calculation
    Matrix left_cauchy_green = prod( rF, trans( rF ) );

    // Calculating the inverse of the jacobian
    Matrix inv_left_cauchy_green ( dimension, dimension );
    double det_b=0;
    MathUtils<double>::InvertMatrix( left_cauchy_green, inv_left_cauchy_green, det_b);

    if( dimension == 2 )
    {
        // Almansi Strain Calculation
        rStrainVector[0] = 0.5 * (  1.00 - inv_left_cauchy_green( 0, 0 ) );
        rStrainVector[1] = 0.5 * (  1.00 - inv_left_cauchy_green( 1, 1 ) );
        rStrainVector[2] = - inv_left_cauchy_green( 0, 1 ); // xy
    }
    else if( dimension == 3 )
    {

        // Almansi Strain Calculation
        if ( rStrainVector.size() != 6 ) rStrainVector.resize( 6, false );
        rStrainVector[0] = 0.5 * (  1.00 - inv_left_cauchy_green( 0, 0 ) );
        rStrainVector[1] = 0.5 * (  1.00 - inv_left_cauchy_green( 1, 1 ) );
        rStrainVector[2] = 0.5 * (  1.00 - inv_left_cauchy_green( 2, 2 ) );
        rStrainVector[3] = - inv_left_cauchy_green( 0, 1 ); // xy
        rStrainVector[4] = - inv_left_cauchy_green( 1, 2 ); // yz
        rStrainVector[5] = - inv_left_cauchy_green( 0, 2 ); // xz
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
void UpdatedLagrangian::CalculateGreenLagrangeStrain(const Matrix& rF,
        Vector& rStrainVector )
{
    KRATOS_TRY

    const unsigned int dimension  = GetGeometry().WorkingSpaceDimension();

    // Right Cauchy-Green Calculation
    Matrix C ( dimension, dimension );
    noalias( C ) = prod( trans( rF ), rF );

    if( dimension == 2 )
    {
        // Green Lagrange Strain Calculation
        if ( rStrainVector.size() != 3 ) rStrainVector.resize( 3, false );
        rStrainVector[0] = 0.5 * ( C( 0, 0 ) - 1.00 );
        rStrainVector[1] = 0.5 * ( C( 1, 1 ) - 1.00 );
        rStrainVector[2] = C( 0, 1 ); // xy
    }
    else if( dimension == 3 )
    {
        // Green Lagrange Strain Calculation
        if ( rStrainVector.size() != 6 ) rStrainVector.resize( 6, false );
        rStrainVector[0] = 0.5 * ( C( 0, 0 ) - 1.00 );
        rStrainVector[1] = 0.5 * ( C( 1, 1 ) - 1.00 );
        rStrainVector[2] = 0.5 * ( C( 2, 2 ) - 1.00 );
        rStrainVector[3] = C( 0, 1 ); // xy
        rStrainVector[4] = C( 1, 2 ); // yz
        rStrainVector[5] = C( 0, 2 ); // xz
    }
    else
    {
        KRATOS_ERROR <<  "Dimension given is wrong: Something is wrong with the given dimension in function: CalculateGreenLagrangeStrain" << std::endl;
    }

    KRATOS_CATCH( "" )
}



//************************************************************************************
//************************************************************************************

double& UpdatedLagrangian::CalculateIntegrationWeight(double& rIntegrationWeight)
{
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    if( dimension == 2 )
        rIntegrationWeight *= GetProperties()[THICKNESS];

    return rIntegrationWeight;
}


//************************************************************************************
//************************************************************************************

void UpdatedLagrangian::EquationIdVector( EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo )
{
    GeometryType& rGeom = GetGeometry();
    int number_of_nodes = rGeom.size();
    int dimension = rGeom.WorkingSpaceDimension();
    unsigned int dimension_2 = number_of_nodes * dimension;

    if ( rResult.size() != dimension_2 )
        rResult.resize( dimension_2, false );

    for ( int i = 0; i < number_of_nodes; i++ )
    {
        int index = i * dimension;
        rResult[index] = rGeom[i].GetDof( DISPLACEMENT_X ).EquationId();
        rResult[index + 1] = rGeom[i].GetDof( DISPLACEMENT_Y ).EquationId();

        if ( dimension == 3 )
            rResult[index + 2] = rGeom[i].GetDof( DISPLACEMENT_Z ).EquationId();
    }

}

//************************************************************************************
//************************************************************************************

void UpdatedLagrangian::GetDofList( DofsVectorType& rElementalDofList, ProcessInfo& CurrentProcessInfo )
{
    GeometryType& rGeom = GetGeometry();
    rElementalDofList.resize( 0 );

    for ( unsigned int i = 0; i < rGeom.size(); i++ )
    {
        rElementalDofList.push_back( rGeom[i].pGetDof( DISPLACEMENT_X ) );
        rElementalDofList.push_back( rGeom[i].pGetDof( DISPLACEMENT_Y ) );

        if ( rGeom.WorkingSpaceDimension() == 3 )
        {
            rElementalDofList.push_back( rGeom[i].pGetDof( DISPLACEMENT_Z ) );
        }
    }

}


//************************************************************************************
//*******************DAMPING MATRIX***************************************************

void UpdatedLagrangian::CalculateDampingMatrix( MatrixType& rDampingMatrix, ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    //0.-Initialize the DampingMatrix:
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    //resizing as needed the LHS
    unsigned int matrix_size = number_of_nodes * dimension;

    if ( rDampingMatrix.size1() != matrix_size )
        rDampingMatrix.resize( matrix_size, matrix_size, false );

    noalias( rDampingMatrix ) = ZeroMatrix(matrix_size, matrix_size);

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

void UpdatedLagrangian::CalculateMassMatrix( MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    // Call the values of the shape function for the single element
    GeneralVariables Variables;
    this->InitializeGeneralVariables(Variables,rCurrentProcessInfo);

    // Lumped
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    unsigned int matrix_size = dimension * number_of_nodes;

    if ( rMassMatrix.size1() != matrix_size )
        rMassMatrix.resize( matrix_size, matrix_size, false );

    rMassMatrix = ZeroMatrix(matrix_size, matrix_size);

    // TOTAL MASS OF ONE MP ELEMENT
    const double & TotalMass = this->GetValue(MP_MASS);

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
Matrix& UpdatedLagrangian::MPMJacobian( Matrix& rResult, const array_1d<double,3>& rPoint)
{
    KRATOS_TRY

    // Derivatives of shape functions
    Matrix shape_functions_gradients;
    shape_functions_gradients =this->MPMShapeFunctionsLocalGradients(
                                   shape_functions_gradients);

    const GeometryType& rGeom = GetGeometry();
    const unsigned int number_nodes = rGeom.PointsNumber();
    const unsigned int dimension = rGeom.WorkingSpaceDimension();

    if (dimension ==2)
    {
        rResult.resize( 2, 2, false );
        rResult = ZeroMatrix(2,2);

        for ( unsigned int i = 0; i < number_nodes; i++ )
        {
            rResult( 0, 0 ) += ( rGeom.GetPoint( i ).X() *  shape_functions_gradients( i, 0 ) );
            rResult( 0, 1 ) += ( rGeom.GetPoint( i ).X() *  shape_functions_gradients( i, 1 ) );
            rResult( 1, 0 ) += ( rGeom.GetPoint( i ).Y() *  shape_functions_gradients( i, 0 ) );
            rResult( 1, 1 ) += ( rGeom.GetPoint( i ).Y() *  shape_functions_gradients( i, 1 ) );
        }
    }
    else if(dimension ==3)
    {
        rResult.resize( 3, 3, false );
        rResult = ZeroMatrix(3,3);

        for ( unsigned int i = 0; i < number_nodes; i++ )
        {
            rResult( 0, 0 ) += ( rGeom.GetPoint( i ).X() *  shape_functions_gradients( i, 0 ) );
            rResult( 0, 1 ) += ( rGeom.GetPoint( i ).X() *  shape_functions_gradients( i, 1 ) );
            rResult( 0, 2 ) += ( rGeom.GetPoint( i ).X() *  shape_functions_gradients( i, 2 ) );
            rResult( 1, 0 ) += ( rGeom.GetPoint( i ).Y() *  shape_functions_gradients( i, 0 ) );
            rResult( 1, 1 ) += ( rGeom.GetPoint( i ).Y() *  shape_functions_gradients( i, 1 ) );
            rResult( 1, 2 ) += ( rGeom.GetPoint( i ).Y() *  shape_functions_gradients( i, 2 ) );
            rResult( 2, 0 ) += ( rGeom.GetPoint( i ).Z() *  shape_functions_gradients( i, 0 ) );
            rResult( 2, 1 ) += ( rGeom.GetPoint( i ).Z() *  shape_functions_gradients( i, 1 ) );
            rResult( 2, 2 ) += ( rGeom.GetPoint( i ).Z() *  shape_functions_gradients( i, 2 ) );
        }
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
Matrix& UpdatedLagrangian::MPMJacobianDelta( Matrix& rResult, const array_1d<double,3>& rPoint, const Matrix & rDeltaPosition )
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
        rResult = ZeroMatrix(2,2);

        for ( unsigned int i = 0; i < rGeom.size(); i++ )
        {
            rResult( 0, 0 ) += ( rGeom.GetPoint( i ).X() + rDeltaPosition(i,0)) * ( shape_functions_gradients( i, 0 ) );
            rResult( 0, 1 ) += ( rGeom.GetPoint( i ).X() + rDeltaPosition(i,0)) * ( shape_functions_gradients( i, 1 ) );
            rResult( 1, 0 ) += ( rGeom.GetPoint( i ).Y() + rDeltaPosition(i,1)) * ( shape_functions_gradients( i, 0 ) );
            rResult( 1, 1 ) += ( rGeom.GetPoint( i ).Y() + rDeltaPosition(i,1)) * ( shape_functions_gradients( i, 1 ) );
        }
    }
    else if(dimension ==3)
    {
        rResult.resize( 3, 3, false );
        rResult = ZeroMatrix(3,3);
        for ( unsigned int i = 0; i < rGeom.size(); i++ )
        {
            rResult( 0, 0 ) += ( rGeom.GetPoint( i ).X() + rDeltaPosition(i,0)) * ( shape_functions_gradients( i, 0 ) );
            rResult( 0, 1 ) += ( rGeom.GetPoint( i ).X() + rDeltaPosition(i,0)) * ( shape_functions_gradients( i, 1 ) );
            rResult( 0, 2 ) += ( rGeom.GetPoint( i ).X() + rDeltaPosition(i,0)) * ( shape_functions_gradients( i, 2 ) );
            rResult( 1, 0 ) += ( rGeom.GetPoint( i ).Y() + rDeltaPosition(i,1)) * ( shape_functions_gradients( i, 0 ) );
            rResult( 1, 1 ) += ( rGeom.GetPoint( i ).Y() + rDeltaPosition(i,1)) * ( shape_functions_gradients( i, 1 ) );
            rResult( 1, 2 ) += ( rGeom.GetPoint( i ).Y() + rDeltaPosition(i,1)) * ( shape_functions_gradients( i, 2 ) );
            rResult( 2, 0 ) += ( rGeom.GetPoint( i ).Z() + rDeltaPosition(i,2)) * ( shape_functions_gradients( i, 0 ) );
            rResult( 2, 1 ) += ( rGeom.GetPoint( i ).Z() + rDeltaPosition(i,2)) * ( shape_functions_gradients( i, 1 ) );
            rResult( 2, 2 ) += ( rGeom.GetPoint( i ).Z() + rDeltaPosition(i,2)) * ( shape_functions_gradients( i, 2 ) );
        }
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
Vector& UpdatedLagrangian::MPMShapeFunctionPointValues( Vector& rResult, const array_1d<double,3>& rPoint )
{
    KRATOS_TRY

    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    Vector rPointLocal = ZeroVector(dimension);

    if (dimension == 2)
    {
        rResult.resize(3, false);
        array_1d<double,3> rPointLocal = ZeroVector(3);

        // 1. Obtain the local coordinate of rPoint
        rPointLocal = GetGeometry().PointLocalCoordinates(rPointLocal, rPoint);

        // 2. Get Shape functions: N
        rResult[0] = 1 - rPointLocal[0] - rPointLocal[1] ;
        rResult[1] = rPointLocal[0] ;
        rResult[2] = rPointLocal[1];
    }
    else if (dimension == 3)
    {
        rResult.resize(4, false);
        array_1d<double,3> rPointLocal = ZeroVector(3);

        // 1. Obtain the local coordinate of rPoint
        rPointLocal = GetGeometry().PointLocalCoordinates(rPointLocal, rPoint);

        // 2. Get Shape functions: N
        rResult[0] =  1.0-(rPointLocal[0]+rPointLocal[1]+rPointLocal[2]) ;
        rResult[1] = rPointLocal[0] ;
        rResult[2] = rPointLocal[1];
        rResult[3] = rPointLocal[2];
    }

    return rResult;

    KRATOS_CATCH( "" )
}


Vector& UpdatedLagrangian::MPMLocalCoordinates(Vector& rResult, array_1d<double,3>& rPoint)
{
    KRATOS_TRY

    GeometryType& rGeom = GetGeometry();

    // Only local coordinated of a point in a tetrahedron is computed
    rResult.resize(4,false);

    double x10 = rGeom[1].Coordinates()[0]-rGeom[0].Coordinates()[0];
    double x20 = rGeom[2].Coordinates()[0]-rGeom[0].Coordinates()[0];
    double x30 = rGeom[3].Coordinates()[0]-rGeom[0].Coordinates()[0];
    double y10 = rGeom[1].Coordinates()[1]-rGeom[0].Coordinates()[1];
    double y20 = rGeom[2].Coordinates()[1]-rGeom[0].Coordinates()[1];
    double y30 = rGeom[3].Coordinates()[1]-rGeom[0].Coordinates()[1];
    double z10 = rGeom[1].Coordinates()[2]-rGeom[0].Coordinates()[2];
    double z20 = rGeom[2].Coordinates()[2]-rGeom[0].Coordinates()[2];
    double z30 = rGeom[3].Coordinates()[2]-rGeom[0].Coordinates()[2];

    rResult[3] = (rPoint[0] - rGeom[0].Coordinates()[0])*(y10*z20 - z10*y20) - (rPoint[1] - rGeom[0].Coordinates()[1])*(x10*z20-x20*z10) + (rPoint[2] - rGeom[0].Coordinates()[2])*(y20*x10 - y10*x20)/mDeterminantJ0;

    rResult[2] = (rPoint[0] - rGeom[0].Coordinates()[0])*(y30*z10-y10*z30) + (rPoint[1] - rGeom[0].Coordinates()[1])*(x10*z30-x30*z10) + (rPoint[2] - rGeom[0].Coordinates()[2])*(y10*x30 - y30*x10)/mDeterminantJ0;

    rResult[1] = (rPoint[0] - rGeom[0].Coordinates()[0])*(y20*z30-y30*z20) + (rPoint[1] - rGeom[0].Coordinates()[1])*(x30*z20-x20*z30) + (rPoint[2] - rGeom[0].Coordinates()[2])*(y30*x20 - x30*y20)/mDeterminantJ0;

    rResult[0] = 1 - rResult[1] - rResult[2] -rResult[3];

    return rResult;

    KRATOS_CATCH( "" )
}


// Function which return dN/de
Matrix& UpdatedLagrangian::MPMShapeFunctionsLocalGradients( Matrix& rResult )
{
    unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    if (dimension == 2)
    {
        rResult = ZeroMatrix(3,2);
        rResult( 0, 0 ) = -1.0;
        rResult( 0, 1 ) = -1.0;
        rResult( 1, 0 ) =  1.0;
        rResult( 1, 1 ) =  0.0;
        rResult( 2, 0 ) =  0.0;
        rResult( 2, 1 ) =  1.0;
    }
    else if(dimension == 3)
    {
        rResult = ZeroMatrix(4,3);
        rResult(0,0) = -1.0;
        rResult(0,1) = -1.0;
        rResult(0,2) = -1.0;
        rResult(1,0) =  1.0;
        rResult(1,1) =  0.0;
        rResult(1,2) =  0.0;
        rResult(2,0) =  0.0;
        rResult(2,1) =  1.0;
        rResult(2,2) =  0.0;
        rResult(3,0) =  0.0;
        rResult(3,1) =  0.0;
        rResult(3,2) =  1.0;
    }

    return rResult;
}

//************************************************************************************
//************************************************************************************

void UpdatedLagrangian::GetValuesVector( Vector& values, int Step )
{
    GeometryType& rGeom = GetGeometry();
    const unsigned int number_of_nodes = rGeom.size();
    const unsigned int dimension = rGeom.WorkingSpaceDimension();
    unsigned int matrix_size = number_of_nodes * dimension;

    if ( values.size() != matrix_size ) values.resize( matrix_size, false );

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        unsigned int index = i * dimension;
        values[index] = rGeom[i].FastGetSolutionStepValue( DISPLACEMENT_X, Step );
        values[index + 1] = rGeom[i].FastGetSolutionStepValue( DISPLACEMENT_Y, Step );

        if ( dimension == 3 )
            values[index + 2] = rGeom[i].FastGetSolutionStepValue( DISPLACEMENT_Z, Step );
    }
}

//************************************************************************************
//************************************************************************************

void UpdatedLagrangian::GetFirstDerivativesVector( Vector& values, int Step )
{
    GeometryType& rGeom = GetGeometry();
    const unsigned int number_of_nodes = rGeom.size();
    const unsigned int dimension = rGeom.WorkingSpaceDimension();
    unsigned int matrix_size = number_of_nodes * dimension;

    if ( values.size() != matrix_size ) values.resize( matrix_size, false );

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        unsigned int index = i * dimension;
        values[index] = rGeom[i].FastGetSolutionStepValue( VELOCITY_X, Step );
        values[index + 1] = rGeom[i].FastGetSolutionStepValue( VELOCITY_Y, Step );

        if ( dimension == 3 )
            values[index + 2] = rGeom[i].FastGetSolutionStepValue( VELOCITY_Z, Step );
    }
}

//************************************************************************************
//************************************************************************************

void UpdatedLagrangian::GetSecondDerivativesVector( Vector& values, int Step )
{
    GeometryType& rGeom = GetGeometry();
    const unsigned int number_of_nodes = rGeom.size();
    const unsigned int dimension = rGeom.WorkingSpaceDimension();
    unsigned int matrix_size = number_of_nodes * dimension;

    if ( values.size() != matrix_size ) values.resize( matrix_size, false );

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        unsigned int index = i * dimension;
        values[index] = rGeom[i].FastGetSolutionStepValue( ACCELERATION_X, Step );
        values[index + 1] = rGeom[i].FastGetSolutionStepValue( ACCELERATION_Y, Step );

        if ( dimension == 3 )
            values[index + 2] = rGeom[i].FastGetSolutionStepValue( ACCELERATION_Z, Step );
    }
}
//************************************************************************************
//************************************************************************************

void UpdatedLagrangian::GetHistoricalVariables( GeneralVariables& rVariables )
{
    //Deformation Gradient F ( set to identity )
    unsigned int size =  rVariables.F.size1();
    rVariables.detF  = 1;
    rVariables.F     = IdentityMatrix(size);

    rVariables.detF0 = mDeterminantF0;
    rVariables.F0    = mDeformationGradientF0;
}


//*************************DECIMAL CORRECTION OF STRAINS******************************
//************************************************************************************

void UpdatedLagrangian::DecimalCorrection(Vector& rVector)
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
/**
 * This function provides the place to perform checks on the completeness of the input.
 * It is designed to be called only once (or anyway, not often) typically at the beginning
 * of the calculations, so to verify that nothing is missing from the input
 * or that no common error is found.
 * @param rCurrentProcessInfo
 */
int  UpdatedLagrangian::Check( const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    GeometryType& rGeom = GetGeometry();
    const unsigned int dimension = rGeom.WorkingSpaceDimension();

    // Verify compatibility with the constitutive law
    ConstitutiveLaw::Features LawFeatures;

    this->GetProperties().GetValue( CONSTITUTIVE_LAW )->GetLawFeatures(LawFeatures);

    bool correct_strain_measure = false;
    for(unsigned int i=0; i<LawFeatures.mStrainMeasures.size(); i++)
    {
        if(LawFeatures.mStrainMeasures[i] == ConstitutiveLaw::StrainMeasure_Deformation_Gradient)
            correct_strain_measure = true;
    }

    KRATOS_ERROR_IF(correct_strain_measure == false ) << "Constitutive law is not compatible with the element type: Large Displacements " << std::endl;

    // Verify that the variables are correctly initialized
    KRATOS_ERROR_IF( VELOCITY.Key() == 0 ) << "VELOCITY has Key zero! (check if the application is correctly registered" << std::endl;
    KRATOS_ERROR_IF( DISPLACEMENT.Key() == 0 ) << "DISPLACEMENT has Key zero! (check if the application is correctly registered" << std::endl;
    KRATOS_ERROR_IF( ACCELERATION.Key() == 0 ) << "ACCELERATION has Key zero! (check if the application is correctly registered" << std::endl;
    KRATOS_ERROR_IF( DENSITY.Key() == 0 ) <<  "DENSITY has Key zero! (check if the application is correctly registered" << std::endl;

    // Verify that the dofs exist
    for ( unsigned int i = 0; i < rGeom.size(); i++ )
    {
        KRATOS_ERROR_IF( rGeom[i].SolutionStepsDataHas( DISPLACEMENT ) == false ) << "Missing variable DISPLACEMENT on node " << rGeom[i].Id() << std::endl;
        KRATOS_ERROR_IF( rGeom[i].HasDofFor( DISPLACEMENT_X ) == false || rGeom[i].HasDofFor( DISPLACEMENT_Y ) == false || rGeom[i].HasDofFor( DISPLACEMENT_Z ) == false ) << "Missing one of the dofs for the variable DISPLACEMENT on node " << rGeom[i].Id() << std::endl;
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
            KRATOS_ERROR_IF(THICKNESS.Key() == 0 ) << "THICKNESS has Key zero! (check if the application is correctly registered" << std::endl;
            KRATOS_ERROR_IF( this->GetProperties().Has( THICKNESS ) == false ) << "THICKNESS not provided for element " << this->Id() << std::endl;
        }
        else
        {
            KRATOS_ERROR_IF( this->GetProperties().GetValue( CONSTITUTIVE_LAW )->GetStrainSize() != 6 ) << "Wrong constitutive law used. This is a 3D element! expected strain size is 6 (el id = ) " << this->Id() << std::endl;
        }

        // Check constitutive law
        this->GetProperties().GetValue( CONSTITUTIVE_LAW )->Check( this->GetProperties(), rGeom, rCurrentProcessInfo );

    }

    return 0;

    KRATOS_CATCH( "" );
}

void UpdatedLagrangian::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Element )

    rSerializer.save("ConstitutiveLawVector",mConstitutiveLawVector);
    rSerializer.save("DeformationGradientF0",mDeformationGradientF0);
    rSerializer.save("DeterminantF0",mDeterminantF0);
    rSerializer.save("InverseJ0",mInverseJ0);
    rSerializer.save("DeterminantJ0",mDeterminantJ0);

}

void UpdatedLagrangian::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Element )
    rSerializer.load("ConstitutiveLawVector",mConstitutiveLawVector);
    rSerializer.load("DeformationGradientF0",mDeformationGradientF0);
    rSerializer.load("DeterminantF0",mDeterminantF0);
    rSerializer.load("InverseJ0",mInverseJ0);
    rSerializer.load("DeterminantJ0",mDeterminantJ0);
}


} // Namespace Kratos

