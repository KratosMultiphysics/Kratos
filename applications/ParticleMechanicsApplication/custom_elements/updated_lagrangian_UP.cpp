//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		BSD License
//					Kratos default license: kratos/license.txt
//
//  Main authors:    Ilaria Iaconeta
//


// System includes

// External includes


// Project includes
#include "includes/define.h"
#include "custom_elements/updated_lagrangian_UP.hpp"
#include "utilities/math_utils.h"
#include "includes/constitutive_law.h"
#include "particle_mechanics_application.h"


#include <omp.h>
#include <sstream>

namespace Kratos
{

///**
//* Flags related to the element computation
//*/
//KRATOS_CREATE_LOCAL_FLAG( UpdatedLagrangian, COMPUTE_RHS_VECTOR,                 0 );
//KRATOS_CREATE_LOCAL_FLAG( UpdatedLagrangian, COMPUTE_LHS_MATRIX,                 1 );
//KRATOS_CREATE_LOCAL_FLAG( UpdatedLagrangian, COMPUTE_RHS_VECTOR_WITH_COMPONENTS, 2 );
//KRATOS_CREATE_LOCAL_FLAG( UpdatedLagrangian, COMPUTE_LHS_MATRIX_WITH_COMPONENTS, 3 );

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

UpdatedLagrangianUP::UpdatedLagrangianUP( )
    : UpdatedLagrangian( )
{
    //DO NOT CALL IT: only needed for Register and Serialization!!!
}
//******************************CONSTRUCTOR*******************************************
//************************************************************************************
UpdatedLagrangianUP::UpdatedLagrangianUP( IndexType NewId, GeometryType::Pointer pGeometry )
    : UpdatedLagrangian( NewId, pGeometry )
{
    //DO NOT ADD DOFS HERE!!!
}

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

UpdatedLagrangianUP::UpdatedLagrangianUP( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
    : UpdatedLagrangian( NewId, pGeometry, pProperties )
{
    mFinalizedStep = true;


}
//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

UpdatedLagrangianUP::UpdatedLagrangianUP( UpdatedLagrangianUP const& rOther)
    :UpdatedLagrangian(rOther)
     //,mDeformationGradientF0(rOther.mDeformationGradientF0)
     //,mDeterminantF0(rOther.mDeterminantF0)
{
}

//*******************************ASSIGMENT OPERATOR***********************************
//************************************************************************************

UpdatedLagrangianUP&  UpdatedLagrangianUP::operator=(UpdatedLagrangianUP const& rOther)
{
    //UpdatedLagrangian::operator=(rOther);

    //mDeformationGradientF0.clear();
    //mDeformationGradientF0 = rOther.mDeformationGradientF0;


    mDeterminantF0 = rOther.mDeterminantF0;


    return *this;
}
//*********************************OPERATIONS*****************************************
//************************************************************************************

Element::Pointer UpdatedLagrangianUP::Create( IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties ) const
{
    return Element::Pointer( new UpdatedLagrangianUP( NewId, GetGeometry().Create( ThisNodes ), pProperties ) );
}
//************************************CLONE*******************************************
//************************************************************************************

Element::Pointer UpdatedLagrangianUP::Clone( IndexType NewId, NodesArrayType const& rThisNodes ) const
{

    UpdatedLagrangianUP NewElement (NewId, GetGeometry().Create( rThisNodes ), pGetProperties() );

    //-----------//





    NewElement.mConstitutiveLawVector = mConstitutiveLawVector->Clone();



    //-----------//



    NewElement.mDeformationGradientF0 = mDeformationGradientF0;


    NewElement.mDeterminantF0 = mDeterminantF0;

    return Element::Pointer( new UpdatedLagrangianUP(NewElement) );
}
//*******************************DESTRUCTOR*******************************************
//************************************************************************************
UpdatedLagrangianUP::~UpdatedLagrangianUP()
{
}


//************************************************************************************
//************************************************************************************

void UpdatedLagrangianUP::Initialize()
{
    KRATOS_TRY


    array_1d<double,3>& xg = this->GetValue(GAUSS_COORD);



    const unsigned int dim = GetGeometry().WorkingSpaceDimension();

    //Constitutive Law initialisation

    //if ( mConstitutiveLawVector.size() != 1 )
    //{
    //mConstitutiveLawVector.resize( 1 );

    //}

    mDeterminantF0 = 1;

    mDeformationGradientF0 = identity_matrix<double> (dim);


    //Compute jacobian inverses

    Matrix J0 = ZeroMatrix(dim, dim);

    J0 = this->MPMJacobian(J0, xg);

    //calculating and storing inverse and the determinant of the jacobian
    MathUtils<double>::InvertMatrix( J0, mInverseJ0, mDeterminantJ0 );

    Matrix j = ZeroMatrix(dim,dim);
    j = this->MPMJacobian(j,xg);
    double detj;
    MathUtils<double>::InvertMatrix( j, mInverseJ, detj );

    InitializeMaterial();


    //double MP_KineticEnergy = 0.0;
    //double MP_StrainEnergy = 0.0;

    //for(unsigned int k = 0;k<3;k++)
    //{
    //MP_KineticEnergy += 0.5 * this->GetValue(MP_MASS) * this->GetValue(MP_VELOCITY)[k] * this->GetValue(MP_VELOCITY)[k] ;
    //}
    //for(unsigned int j = 0; j < this->GetValue(MP_CAUCHY_STRESS_VECTOR).size(); j++)
    //{
    //MP_StrainEnergy +=  0.5 * this->GetValue(MP_VOLUME) * this->GetValue(MP_CAUCHY_STRESS_VECTOR)[j] * this->GetValue(MP_ALMANSI_STRAIN_VECTOR)[j];
    //}

    //this->GetValue(MP_KINETIC_ENERGY) = MP_KineticEnergy;
    //this->GetValue(MP_STRAIN_ENERGY) = MP_StrainEnergy;
    this->GetValue(MP_DENSITY) = GetProperties()[DENSITY];

    //std::cout<<"The element is initialized"<<std::endl;


    KRATOS_CATCH( "" )
}
//************************************************************************************
//************************************************************************************

void UpdatedLagrangianUP::InitializeGeneralVariables (GeneralVariables& rVariables, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    UpdatedLagrangian::InitializeGeneralVariables(rVariables,rCurrentProcessInfo);

    //stabilization factor
    double StabilizationFactor = 1.0;
    if( GetProperties().Has(STABILIZATION_FACTOR) )
    {
        StabilizationFactor = GetProperties()[STABILIZATION_FACTOR];
    }
    else if( rCurrentProcessInfo.Has(STABILIZATION_FACTOR) )
    {
        StabilizationFactor = rCurrentProcessInfo[STABILIZATION_FACTOR];
    }
    GetProperties().SetValue(STABILIZATION_FACTOR, StabilizationFactor);

    //ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);
    //Flags &ConstitutiveLawOptions=Values.GetOptions();

    ////std::cout<<"in CalculateElementalSystem 5"<<std::endl;
    //ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);

    //ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS);

    //ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
    //ConstitutiveLawOptions.Set(ConstitutiveLaw::ISOCHORIC_TENSOR_ONLY);



    KRATOS_CATCH( "" )

}
//************************************************************************************
//************************************************************************************
//void UpdatedLagrangianUP::UpdateGaussPoint (GeneralVariables& rVariables, const ProcessInfo& rCurrentProcessInfo)
//{
//KRATOS_TRY

//UpdatedLagrangian::UpdateGaussPoint(rVariables,rCurrentProcessInfo);

//const unsigned int number_of_nodes = GetGeometry().PointsNumber();
//unsigned int dimension = GetGeometry().WorkingSpaceDimension();

//array_1d<double,3> xg = this->GetValue(GAUSS_COORD);
//rVariables.N = this->MPMShapeFunctionPointValues(rVariables.N, xg);
//double MP_Pressure = 0.0;//this->GetValue(MP_PRESSURE);

//for ( unsigned int i = 0; i < number_of_nodes; i++ )
//{
//if (rVariables.N[i] > 1e-16)
//{
//double NodalPressure = GetGeometry()[i].GetSolutionStepValue(PRESSURE, 0);

//MP_Pressure += rVariables.N[i] * NodalPressure;

//}
//}

//KRATOS_CATCH( "" )

//}

/**
    * The position of the Gauss points/Material points is updated
    */

void UpdatedLagrangianUP::UpdateGaussPoint( GeneralVariables & rVariables, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    rVariables.CurrentDisp = CalculateCurrentDisp(rVariables.CurrentDisp, rCurrentProcessInfo);
    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    array_1d<double,3> xg = this->GetValue(GAUSS_COORD);
    array_1d<double,3> MP_PreviousAcceleration = this->GetValue(MP_ACCELERATION);
    array_1d<double,3> MP_PreviousVelocity = this->GetValue(MP_VELOCITY);
    //double MP_Mass = this->GetValue(MP_MASS);
    array_1d<double,3> delta_xg = ZeroVector(3);
    array_1d<double,3> MP_Acceleration = ZeroVector(3);
    array_1d<double,3> MP_Velocity = ZeroVector(3);

    array_1d<double,3> MP_DeltaVelocity = ZeroVector(3);

    double MP_Pressure = 0.0;
    double DeltaTime = rCurrentProcessInfo[DELTA_TIME];


    rVariables.N = this->MPMShapeFunctionPointValues(rVariables.N, xg);
    //int MP_number = this->GetValue(MP_NUMBER);

    //double total_nodal_mass = 0.0;
    //for ( unsigned int i = 0; i < number_of_nodes; i++ )
    //{
    //total_nodal_mass += GetGeometry()[i].GetSolutionStepValue(NODAL_MASS, 0);
    //}
    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        if (rVariables.N[i] > 1e-16)
        {
            array_1d<double, 3 > & NodalAcceleration = GetGeometry()[i].FastGetSolutionStepValue(ACCELERATION);
            array_1d<double, 3 > & NodalVelocity = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY);
            //array_1d<double, 3 > & PreviousNodalVelocity = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY,1);
            double NodalMass = GetGeometry()[i].FastGetSolutionStepValue(NODAL_MASS, 0);
            array_1d<double,3> NodalMomentum = NodalMass * NodalVelocity;
            array_1d<double,3> NodalInertia = NodalMass * NodalAcceleration;

            double NodalPressure = GetGeometry()[i].FastGetSolutionStepValue(PRESSURE, 0);
            MP_Pressure += rVariables.N[i] * NodalPressure;



            for ( unsigned int j = 0; j < dimension; j++ )
            {

                delta_xg[j] += rVariables.N[i] * rVariables.CurrentDisp(i,j);
                MP_Acceleration[j] += rVariables.N[i] * NodalAcceleration[j];
                //MP_Velocity[j] += rVariables.N[i] * NodalVelocity[j];
                //MP_DeltaVelocity[j] += rVariables.N[i] * (NodalVelocity[j]-PreviousNodalVelocity[j]);
                //MP_Acceleration[j] +=NodalInertia[j]/(rVariables.N[i] * MP_Mass * MP_number);//
                //MP_Velocity[j] += NodalMomentum[j]/(rVariables.N[i] * MP_Mass * MP_number);
                //MP_Velocity[j] += DeltaTime * rVariables.N[i] * NodalAcceleration[j];////




            }
        }

    }


    //**************************************************************************************************************************
    //Another way to update the MP velocity (see paper Guilkey and Weiss, 2003)
    MP_Velocity = MP_PreviousVelocity + 0.5 * DeltaTime * (MP_Acceleration + MP_PreviousAcceleration);
    //MP_Acceleration = 2.0/DeltaTime * (MP_Velocity - MP_PreviousVelocity) - MP_PreviousAcceleration;
    //MP_Acceleration = 4/(DeltaTime * DeltaTime) * delta_xg - 4/DeltaTime * MP_PreviousVelocity;
    //MP_Velocity = 2/DeltaTime * delta_xg - MP_PreviousVelocity;

    //MP_Velocity = MP_PreviousVelocity + MP_DeltaVelocity;
    this -> SetValue(MP_PRESSURE,MP_Pressure );
    this -> SetValue(MP_VELOCITY,MP_Velocity );

    const array_1d<double,3>& new_xg = xg + delta_xg ;

    //Update the MP Position
    this -> SetValue(GAUSS_COORD,new_xg);

    //Update the MP Acceleration
    this -> SetValue(MP_ACCELERATION,MP_Acceleration);

    array_1d<double,3>& MP_Displacement = this->GetValue(MP_DISPLACEMENT);

    MP_Displacement += delta_xg;

    //Update the MP Displacement
    this -> SetValue(MP_DISPLACEMENT,MP_Displacement );

    KRATOS_CATCH( "" )
}
//************************************************************************************
//************************************************************************************

//void UpdatedLagrangian::SetGeneralVariables(GeneralVariables& rVariables,
//ConstitutiveLaw::Parameters& rValues)
//{
////Variables.detF is the determinant of the incremental total deformation gradient
//rVariables.detF  = MathUtils<double>::Det(rVariables.F);

//if(rVariables.detF<0){

//std::cout<<" Element: "<<this->Id()<<std::endl;
//std::cout<<" Element position "<<this->GetValue(GAUSS_COORD)<<std::endl;
//unsigned int number_of_nodes = GetGeometry().PointsNumber();

//for ( unsigned int i = 0; i < number_of_nodes; i++ )
//{
//array_1d<double, 3> &CurrentPosition  = GetGeometry()[i].Coordinates();
//array_1d<double, 3> & CurrentDisplacement  = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT);
//array_1d<double, 3> & PreviousDisplacement = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT,1);
////array_1d<double, 3> PreviousPosition  = CurrentPosition - (CurrentDisplacement-PreviousDisplacement);
//std::cout<<" NODE ["<<GetGeometry()[i].Id()<<"]: (Current position: "<<CurrentPosition<<") "<<std::endl;
//std::cout<<" ---Current Disp: "<<CurrentDisplacement<<" (Previour Disp: "<<PreviousDisplacement<<")"<<std::endl;
//}

//for ( unsigned int i = 0; i < number_of_nodes; i++ )
//{
//if( GetGeometry()[i].SolutionStepsDataHas(CONTACT_FORCE) ){
//array_1d<double, 3 > & PreContactForce = GetGeometry()[i].FastGetSolutionStepValue(CONTACT_FORCE,1);
//array_1d<double, 3 > & ContactForce = GetGeometry()[i].FastGetSolutionStepValue(CONTACT_FORCE);
//std::cout<<" ---Contact_Force: (Pre:"<<PreContactForce<<", Cur:"<<ContactForce<<") "<<std::endl;
//}
//else{
//std::cout<<" ---Contact_Force: NULL "<<std::endl;
//}
//}

//KRATOS_THROW_ERROR( std::invalid_argument," MPM UPDATED LAGRANGIAN DISPLACEMENT ELEMENT INVERTED: |F|<0  detF = ", rVariables.detF )
//}

//rVariables.detFT = rVariables.detF * rVariables.detF0;
//rVariables.FT    = prod( rVariables.F, rVariables.F0 );


//rValues.SetDeterminantF(rVariables.detFT);
//rValues.SetDeformationGradientF(rVariables.FT);
//rValues.SetStrainVector(rVariables.StrainVector);
//rValues.SetStressVector(rVariables.StressVector);
//rValues.SetConstitutiveMatrix(rVariables.ConstitutiveMatrix);
//rValues.SetShapeFunctionsDerivatives(rVariables.DN_DX);
//rValues.SetShapeFunctionsValues(rVariables.N);


////std::cout<<"The general variables are set"<<std::endl;

//}
//************************************************************************************
//*****************check size of LHS and RHS matrices*********************************

void UpdatedLagrangianUP::InitializeSystemMatrices(MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        Flags& rCalculationFlags)

{

    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

    //resizing as needed the LHS
    unsigned int MatSize = number_of_nodes * dimension + number_of_nodes;

    if ( rCalculationFlags.Is(UpdatedLagrangian::COMPUTE_LHS_MATRIX) ) //calculation of the matrix is required
    {
        if ( rLeftHandSideMatrix.size1() != MatSize )
            rLeftHandSideMatrix.resize( MatSize, MatSize, false );

        noalias( rLeftHandSideMatrix ) = ZeroMatrix( MatSize, MatSize ); //resetting LHS
    }


    //resizing as needed the RHS
    if ( rCalculationFlags.Is(UpdatedLagrangian::COMPUTE_RHS_VECTOR) ) //calculation of the matrix is required
    {
        if ( rRightHandSideVector.size() != MatSize )
            rRightHandSideVector.resize( MatSize, false );

        rRightHandSideVector = ZeroVector( MatSize ); //resetting RHS
    }
    //std::cout<<"The system matrices are initialized"<<std::endl;
}

//************************************************************************************
//************************************************************************************

//void UpdatedLagrangian::CalculateElementalSystem( LocalSystemComponents& rLocalSystem,
//ProcessInfo& rCurrentProcessInfo)
//{
//KRATOS_TRY

////create and initialize element variables:
//GeneralVariables Variables;

//this->InitializeGeneralVariables(Variables,rCurrentProcessInfo);


////create constitutive law parameters:
//ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);


////set constitutive law flags:
//Flags &ConstitutiveLawOptions=Values.GetOptions();

////std::cout<<"in CalculateElementalSystem 5"<<std::endl;
//ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);

//ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS);

//ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);


////auxiliary terms
//Vector VolumeForce;


////compute element kinematics B, F, DN_DX ...
//this->CalculateKinematics(Variables,rCurrentProcessInfo);

////set general variables to constitutivelaw parameters
//this->SetGeneralVariables(Variables,Values);

//mConstitutiveLawVector->CalculateMaterialResponse(Values, Variables.StressMeasure);

////this->SetValue(MP_CAUCHY_STRESS_VECTOR, Variables.StressVector);
////std::cout<<"AAAAAAAAAAAAAAAAAAAAAAAAAA"<<std::endl;
////std::cout<<"Variables.StressVector in the element "<<Variables.StressVector<<std::endl;

////this->SetValue(MP_ALMANSI_STRAIN_VECTOR, Variables.StrainVector);
////double EquivalentPlasticStrain = mConstitutiveLawVector->GetValue( PLASTIC_STRAIN, EquivalentPlasticStrain );
////this->SetValue(MP_EQUIVALENT_PLASTIC_STRAIN, EquivalentPlasticStrain);
////at the first iteration I recover the previous state of stress and strain
//if(rCurrentProcessInfo[NL_ITERATION_NUMBER] == 1)
//{
//this->SetValue(PREVIOUS_MP_CAUCHY_STRESS_VECTOR, Variables.StressVector);
//this->SetValue(PREVIOUS_MP_ALMANSI_STRAIN_VECTOR, Variables.StrainVector);
//}
////the MP density is updated
//double MP_Density = (GetProperties()[DENSITY]) / Variables.detFT;

////the integration weight is evaluated
//double MP_Volume = this->GetValue(MP_MASS)/this->GetValue(MP_DENSITY);

//this->SetValue(MP_DENSITY, MP_Density);
//this->SetValue(MP_VOLUME, MP_Volume);




//if ( rLocalSystem.CalculationFlags.Is(UpdatedLagrangian::COMPUTE_LHS_MATRIX) ) //calculation of the matrix is required
//{

////contributions to stiffness matrix calculated on the reference config
//this->CalculateAndAddLHS ( rLocalSystem, Variables, MP_Volume );

//}

//if ( rLocalSystem.CalculationFlags.Is(UpdatedLagrangian::COMPUTE_RHS_VECTOR) ) //calculation of the vector is required
//{
////contribution to external forces
//VolumeForce  = this->CalculateVolumeForce( VolumeForce, Variables );

//this->CalculateAndAddRHS ( rLocalSystem, Variables, VolumeForce, MP_Volume );

//}



//KRATOS_CATCH( "" )
//}
//*********************************COMPUTE KINEMATICS*********************************
//************************************************************************************


void UpdatedLagrangianUP::CalculateKinematics(GeneralVariables& rVariables, ProcessInfo& rCurrentProcessInfo)

{
    KRATOS_TRY

    //const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    //Define the stress measure
    rVariables.StressMeasure = ConstitutiveLaw::StressMeasure_Cauchy;

    //Calculating the inverse of the jacobian and the parameters needed [d£/dx_n]
    Matrix InvJ;

    MathUtils<double>::InvertMatrix( rVariables.J, InvJ, rVariables.detJ);



    //Calculating the inverse of the jacobian and the parameters needed [d£/(dx_n+1)]
    Matrix Invj;
    MathUtils<double>::InvertMatrix( rVariables.j, Invj, rVariables.detJ ); //overwrites detJ

    //Compute cartesian derivatives [dN/dx_n+1]
    rVariables.DN_DX = prod( rVariables.DN_De, Invj); //overwrites DX now is the current position dx



    //Deformation Gradient F [(dx_n+1 - dx_n)/dx_n] to be updated in constitutive law parameter as total deformation gradient
    //the increment of total deformation gradient can be evaluated in 2 ways.
    //1 way.
    noalias( rVariables.F ) = prod( rVariables.j, InvJ);

    //2 way by means of the gradient of nodal displacement: using this second expression quadratic convergence is not guarantee

    //Matrix I=identity_matrix<double>( dimension );

    //Matrix GradientDisp = ZeroMatrix(dimension, dimension);
    //rVariables.CurrentDisp = CalculateCurrentDisp(rVariables.CurrentDisp, rCurrentProcessInfo);
    //GradientDisp = prod(trans(rVariables.CurrentDisp),rVariables.DN_DX);

    ////REMEMBER THAT USING JUST ONLY THE FIRST ORDER TERM SOME ISSUES CAN COME UP WHEN FOR PROBLEMS WITH LOTS OF ROTATIONAL MOTION(slender cantilever beam??)
    //noalias( rVariables.F ) = (I + GradientDisp);
  
    //Determinant of the Deformation Gradient F_n

    rVariables.detF0 = mDeterminantF0;
    rVariables.F0    = mDeformationGradientF0;

    //Compute the deformation matrix B
    this->CalculateDeformationMatrix(rVariables.B, rVariables.F, rVariables.DN_DX);


    KRATOS_CATCH( "" )
}
//************************************************************************************

void UpdatedLagrangianUP::CalculateDeformationMatrix(Matrix& rB,
        Matrix& rF,
        Matrix& rDN_DX)
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

    rB.clear(); //set all components to zero

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

        KRATOS_THROW_ERROR( std::invalid_argument, "something is wrong with the dimension", "" )

    }

    KRATOS_CATCH( "" )
}

////************************************************************************************
////************************************************************************************

void UpdatedLagrangianUP::InitializeSolutionStep( ProcessInfo& rCurrentProcessInfo )
{
    // In the Initialize of each time step the nodal initial conditions are evaluated
    //1. first of all I need to evaluate the MP momentum and MP_inertia



    //int MP_bool = this->GetValue(MP_BOOL);

    //std::cout<<" in InitializeSolutionStep2"<<std::endl;
    unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    array_1d<double,3>& xg = this->GetValue(GAUSS_COORD);
    GeneralVariables Variables;
    //this->InitializeGeneralVariables(Variables,rCurrentProcessInfo);



    Matrix J0 = ZeroMatrix(dimension, dimension);

    J0 = this->MPMJacobian(J0, xg);

    //calculating and storing inverse and the determinant of the jacobian
    MathUtils<double>::InvertMatrix( J0, mInverseJ0, mDeterminantJ0 );

    Variables.N = this->MPMShapeFunctionPointValues(Variables.N, xg);

    mConstitutiveLawVector->InitializeSolutionStep( GetProperties(),
            GetGeometry(), Variables.N, rCurrentProcessInfo );

    mFinalizedStep = false;



    array_1d<double,3>& MP_Velocity = this->GetValue(MP_VELOCITY);
    array_1d<double,3>& MP_Acceleration = this->GetValue(MP_ACCELERATION);
    double MP_Pressure = this->GetValue(MP_PRESSURE);
    array_1d<double,3>& AUX_MP_Velocity = this->GetValue(AUX_MP_VELOCITY);
    array_1d<double,3>& AUX_MP_Acceleration = this->GetValue(AUX_MP_ACCELERATION);
    double AUX_MP_Pressure = this->GetValue(AUX_MP_PRESSURE);
    double MP_Mass = this->GetValue(MP_MASS);
    array_1d<double,3> MP_Momentum;
    array_1d<double,3> MP_Inertia;
    //double MP_MPressure;
    array_1d<double,3> NodalMomentum;
    array_1d<double,3> NodalInertia;
    double NodalMPressure;

    for (unsigned int j=0; j<number_of_nodes; j++)
    {
        //these are the values of nodal velocity and nodal acceleration evaluated in the initialize solution step
        array_1d<double, 3 > & NodalAcceleration = GetGeometry()[j].FastGetSolutionStepValue(ACCELERATION,1);
        array_1d<double, 3 > & NodalVelocity = GetGeometry()[j].FastGetSolutionStepValue(VELOCITY,1);
        double & NodalPressure = GetGeometry()[j].FastGetSolutionStepValue(PRESSURE,1);

        AUX_MP_Pressure += Variables.N[j] * NodalPressure;

        //std::cout<<"NodalVelocity "<< GetGeometry()[j].Id()<<std::endl;
        for (unsigned int k = 0; k < dimension; k++)
        {
            AUX_MP_Velocity[k] += Variables.N[j] * NodalVelocity[k];
            AUX_MP_Acceleration[k] += Variables.N[j] * NodalAcceleration[k];
        }
    }

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {

        NodalMPressure =  Variables.N[i] * (MP_Pressure - AUX_MP_Pressure) * MP_Mass;

        for (unsigned int j = 0; j < dimension; j++)
        {
            NodalMomentum[j] = Variables.N[i] * (MP_Velocity[j] - AUX_MP_Velocity[j]) * MP_Mass;
            NodalInertia[j] = Variables.N[i] * (MP_Acceleration[j] - AUX_MP_Acceleration[j]) * MP_Mass;

        }
        GetGeometry()[i].FastGetSolutionStepValue(NODAL_MOMENTUM, 0) += NodalMomentum;
        GetGeometry()[i].FastGetSolutionStepValue(NODAL_INERTIA, 0) += NodalInertia;
        GetGeometry()[i].FastGetSolutionStepValue(NODAL_MPRESSURE, 0) += NodalMPressure;
        GetGeometry()[i].FastGetSolutionStepValue(NODAL_MASS, 0) += Variables.N[i] * MP_Mass;

    }

    AUX_MP_Velocity.clear();
    AUX_MP_Acceleration.clear();
    AUX_MP_Pressure = 0.0;





}
//************************************************************************************
//************************************************************************************

void UpdatedLagrangianUP::CalculateAndAddRHS(LocalSystemComponents& rLocalSystem, GeneralVariables& rVariables, Vector& rVolumeForce, double& rIntegrationWeight)
{

    //contribution of the internal and external forces
    //contribution of the internal and external forces
    VectorType& rRightHandSideVector = rLocalSystem.GetRightHandSideVector();

    rVariables.detF0   *= rVariables.detF;
    double DeterminantF = rVariables.detF;
    rVariables.detF = 1; //in order to simplify updated and spatial lagrangian

    // operation performed: rRightHandSideVector += ExtForce*IntegrationWeight
    CalculateAndAddExternalForces( rRightHandSideVector, rVariables, rVolumeForce, rIntegrationWeight );

    // operation performed: rRightHandSideVector -= IntForce*IntegrationWeight
    CalculateAndAddInternalForces( rRightHandSideVector, rVariables, rIntegrationWeight);

    // operation performed: rRightHandSideVector -= PressureForceBalance*IntegrationWeight
    CalculateAndAddPressureForces( rRightHandSideVector, rVariables, rIntegrationWeight);

    // operation performed: rRightHandSideVector -= Stabilized Pressure Forces
    CalculateAndAddStabilizedPressure( rRightHandSideVector, rVariables, rIntegrationWeight);

    //std::cout<<" rRightHandSideVector " << rRightHandSideVector<<std::endl;

    rVariables.detF     = DeterminantF;
    rVariables.detF0   /= rVariables.detF;

}
//************************************************************************************
//*********************Calculate the contribution of external force*******************

void UpdatedLagrangianUP::CalculateAndAddExternalForces(VectorType& rRightHandSideVector,
        GeneralVariables& rVariables,
        Vector& rVolumeForce,
        double& rIntegrationWeight)

{
    KRATOS_TRY
    unsigned int number_of_nodes = GetGeometry().PointsNumber();

    unsigned int dimension = GetGeometry().WorkingSpaceDimension();



    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        int indexup = dimension * i + i;

        for ( unsigned int j = 0; j < dimension; j++ )
        {
            rRightHandSideVector[indexup + j] += rVariables.N[i] * rVolumeForce[j];

        }

    }
  
    KRATOS_CATCH( "" )
}
//************************************************************************************
//************************************************************************************

void UpdatedLagrangianUP::CalculateAndAddInternalForces(VectorType& rRightHandSideVector,
        GeneralVariables & rVariables,
        double& rIntegrationWeight)
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    VectorType Fh=rRightHandSideVector;

    //in the largedisplacement UP element is defined Vector
    VectorType InternalForces = rIntegrationWeight * prod( trans( rVariables.B ), rVariables.StressVector );

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        unsigned int indexup = dimension * i + i;
        unsigned int indexu  = dimension * i;

        for ( unsigned int j = 0; j < dimension; j++ )
        {
            rRightHandSideVector[indexup + j] -= InternalForces[indexu + j];
        }
    }

    KRATOS_CATCH( "" )
}

//******************************************************************************************************************
//******************************************************************************************************************
double& UpdatedLagrangianUP::CalculatePUCoefficient(double& rCoefficient, GeneralVariables & rVariables)
{
    KRATOS_TRY

    //Mechanical volumetric:

    //Constitutive A:
    //rCoefficient = 0.5*(rVariables.detF0*rVariables.detF0-1)/rVariables.detF0); //(J²-1)/2

    //Constitutive B:
    //rCoefficient = (std::log(rVariables.detF0)/rVariables.detF0);  //(ln(J)) for a free energy function of 0.5 ln^2(J)

    //Constitutive C:free energy function of 0.5 (j-1)^2
    rCoefficient = rVariables.detF0 - 1;

    return rCoefficient;

    KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************

double& UpdatedLagrangianUP::CalculatePUDeltaCoefficient(double &rDeltaCoefficient, GeneralVariables & rVariables)
{

    KRATOS_TRY

    //Mechanical volumetric:

    //Constitutive A:
    //rDeltaCoefficient = (rVariables.detF0*rVariables.detF0 + 1)/(rVariables.detF0*rVariables.detF0); //(J²-1)/2

    //Constitutive B:
    //rDeltaCoefficient = (1.0-std::log(rVariables.detF0))/(rVariables.detF0*rVariables.detF0);   //(ln(J)) for a free energy function of 0.5K ln^2(J)

    //Constitutive C:free energy function of 0.5 (j-1)^2
    rDeltaCoefficient = 1.0;


    return rDeltaCoefficient;


    KRATOS_CATCH( "" )

}

//************************************************************************************
//************************************************************************************
// I changed the term of the pressure: as after a first prediction of the pressure
// variable, this term is evaluated in the return mapping of the DP plastic model.
// After plasticity the pressure is evaluated directly on the integration point
// here called "NewPressure"
//************************************************************************************
//************************************************************************************

void UpdatedLagrangianUP::CalculateAndAddPressureForces(VectorType& rRightHandSideVector,
        GeneralVariables & rVariables,
        double& rIntegrationWeight)
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    unsigned int indexp = dimension;

    VectorType Fh=rRightHandSideVector;

    double BulkModulus= GetProperties()[YOUNG_MODULUS]/(3*(1-2*GetProperties()[POISSON_RATIO]));

    double DeltaCoefficient = 0;
    DeltaCoefficient = this->CalculatePUDeltaCoefficient( DeltaCoefficient, rVariables );

    double Coefficient = 0;
    Coefficient = this->CalculatePUCoefficient( Coefficient, rVariables );

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        for ( unsigned int j = 0; j < number_of_nodes; j++ )
        {

            double& Pressure = GetGeometry()[j].FastGetSolutionStepValue(PRESSURE);
            //if(GetGeometry()[j].Id() == 2276)
            //{
            //std::cout<< "Pressure_in_calculate_pressure_force "<<Pressure<<std::endl;
            //}

            // consistent=1;
            // if(i==j)
            //   consistent=2;

            // if( dimension == 2 ){ //consistent 2D

            //   rRightHandSideVector[indexp] += consistent * (1.0/BulkModulus) * (1.0/12.0) * Pressure * rIntegrationWeight / (rVariables.detF0/rVariables.detF) ; //2D

            // }
            // else{

            //   rRightHandSideVector[indexp] += consistent * (1.0/BulkModulus) * (1.0/20.0) * Pressure * rIntegrationWeight / (rVariables.detF0/rVariables.detF) ; //3D
            // }

            rRightHandSideVector[indexp] += (1.0/(DeltaCoefficient * BulkModulus)) * rVariables.N[i] * rVariables.N[j] * Pressure * rIntegrationWeight / (rVariables.detF0/rVariables.detF) ; //2D-3D
            //rRightHandSideVector[indexp] += (1.0/BulkModulus) * rVariables.N[i] * NewPressure * rIntegrationWeight / (rVariables.detF0/rVariables.detF) ;
            //*********************************************************
            //mixed formulation considering a variation of p: 1st term
            //rRightHandSideVector[indexp] += rVariables.N[i] * rVariables.N[j] * Pressure * rIntegrationWeight / (rVariables.detF0/rVariables.detF) ;
            //*********************************************************
        }

        rRightHandSideVector[indexp] -=  Coefficient/DeltaCoefficient * rVariables.N[i] * rIntegrationWeight / (rVariables.detF0/rVariables.detF);

        //********************************************************
        //mixed formulation considering a variation of p: 2nd term
        //rRightHandSideVector[indexp] -=  rVariables.N[i] * NewPressure * rIntegrationWeight / (rVariables.detF0/rVariables.detF);
        //********************************************************
        indexp += (dimension + 1);
    }

    KRATOS_CATCH( "" )
}
//************************************************************************************
//************************************************************************************

void UpdatedLagrangianUP::CalculateAndAddStabilizedPressure(VectorType& rRightHandSideVector,
        GeneralVariables & rVariables,
        double& rIntegrationWeight)
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    unsigned int indexp = dimension;

    double DeltaCoefficient = 0;
    DeltaCoefficient = this->CalculatePUDeltaCoefficient( DeltaCoefficient, rVariables );
    VectorType Fh=rRightHandSideVector;

    //use of this variable for the complete parameter:
    double AlphaStabilization  = 1.0;
    double StabilizationFactor = GetProperties()[STABILIZATION_FACTOR];
    AlphaStabilization *= StabilizationFactor;

    const double& YoungModulus          = GetProperties()[YOUNG_MODULUS];
    const double& PoissonCoefficient    = GetProperties()[POISSON_RATIO];

    double LameMu =  YoungModulus/(2*(1+PoissonCoefficient));

    //Experimental
    // if(LameMu < rVariables.ConstitutiveMatrix(2,2))
    //   LameMu = rVariables.ConstitutiveMatrix(2,2);

    //double NewPressure = mConstitutiveLawVector->GetValue(MP_PRESSURE, NewPressure );
    double consistent = 1;
    double FactorValue = 8.0; //JMR deffault value
    if( dimension == 3 )
        FactorValue = 10.0; //JMC deffault value


    //NON FUNZIONA
    //consistent = AlphaStabilization*FactorValue / LameMu;
    ////double stabilization_term = 1;

    //for ( unsigned int i = 0; i < number_of_nodes; i++ )
    //{
    ////for ( unsigned int j = 0; j < number_of_nodes; j++ )
    ////{

    ////double& Pressure = GetGeometry()[j].FastGetSolutionStepValue(PRESSURE);


    ////rRightHandSideVector[indexp] += consistent * ( -1/9)* Pressure * rIntegrationWeight / (rVariables.detF0/rVariables.detF); //2D
    ////}
    //rRightHandSideVector[indexp] += consistent * (rVariables.N[i]*rVariables.N[j] -1/9)* NewPressure * rIntegrationWeight / (rVariables.detF0/rVariables.detF); //2D



    //indexp += (dimension + 1);
    //}
    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        for ( unsigned int j = 0; j < number_of_nodes; j++ )
        {

            double& Pressure = GetGeometry()[j].FastGetSolutionStepValue(PRESSURE);

            if( dimension == 2 )  //consistent 2D
            {

                consistent=(-1)*AlphaStabilization*FactorValue/(36.0*LameMu);
                if(i==j)
                    consistent=2*AlphaStabilization*FactorValue/(36.0*LameMu);

                //rRightHandSideVector[indexp] += consistent * (rVariables.N[i]*rVariables.N[j] -1/9)* Pressure * rIntegrationWeight / (rVariables.detF0/rVariables.detF); //2D
                rRightHandSideVector[indexp] += consistent * Pressure * rIntegrationWeight / (DeltaCoefficient * (rVariables.detF0/rVariables.detF)); //2D

            }
            else
            {

                consistent=(-1)*AlphaStabilization*FactorValue/(80.0*LameMu);
                if(i==j)
                    consistent=3*AlphaStabilization*FactorValue/(80.0*LameMu);

                rRightHandSideVector[indexp] += consistent * Pressure * rIntegrationWeight / (rVariables.detF0/rVariables.detF); //3D

            }

            // std::cout<<" Pressure "<<Pressure<<std::endl;
        }


        indexp += (dimension + 1);
    }



    //for ( unsigned int i = 0; i < number_of_nodes; i++ )
    //{

    //rRightHandSideVector[indexp] += consistent * (rVariables.N[i] -1/3) * NewPressure * rIntegrationWeight / (rVariables.detF0/rVariables.detF);

    //indexp += (dimension + 1);

    //}

    // std::cout<<std::endl;
    // std::cout<<" IntegrationWeight "<<rIntegrationWeight<<" detF "<<rVariables.detF0<<std::endl;
    // std::cout<<" FpStab "<<rRightHandSideVector-Fh<<std::endl;

    KRATOS_CATCH( "" )
}
//************************************************************************************
//************************************************************************************

void UpdatedLagrangianUP::CalculateAndAddLHS(LocalSystemComponents& rLocalSystem, GeneralVariables& rVariables, double& rIntegrationWeight)
{

    //contributions of the stiffness matrix calculated on the reference configuration
    //contributions of the stiffness matrix calculated on the reference configuration
    MatrixType& rLeftHandSideMatrix = rLocalSystem.GetLeftHandSideMatrix();
    rVariables.detF0   *= rVariables.detF;
    double DeterminantF = rVariables.detF;
    rVariables.detF = 1; //in order to simplify updated and spatial lagrangian
    // operation performed: add Km to the rLefsHandSideMatrix
    MatrixType Kh;
    //respect to the current configuration n+1
    CalculateAndAddKuum( rLeftHandSideMatrix, rVariables, rIntegrationWeight );
    //std::cout<<" Kuum " << rLeftHandSideMatrix<<std::endl;
    //Kh=rLeftHandSideMatrix;
    //std::cout<<"111111111111111111111"<<std::endl;
    // operation performed: add Kg to the rLefsHandSideMatrix
    CalculateAndAddKuug( rLeftHandSideMatrix, rVariables, rIntegrationWeight );
    //std::cout<<" Kuug " << rLeftHandSideMatrix-Kh<<std::endl;
    //Kh=rLeftHandSideMatrix;
    //std::cout<<"222222222222222222222"<<std::endl;
    // operation performed: add Kup to the rLefsHandSideMatrix
    CalculateAndAddKup( rLeftHandSideMatrix, rVariables, rIntegrationWeight );
    //std::cout<<" Kup " << rLeftHandSideMatrix-Kh<<std::endl;
    //Kh=rLeftHandSideMatrix;
    //std::cout<<"33333333333333333333"<<std::endl;
    // operation performed: add Kpu to the rLefsHandSideMatrix
    CalculateAndAddKpu( rLeftHandSideMatrix, rVariables, rIntegrationWeight );
    //std::cout<<" Kpu " << rLeftHandSideMatrix-Kh<<std::endl;
    //Kh=rLeftHandSideMatrix;
    //std::cout<<"444444444444444444444"<<std::endl;
    // operation performed: add Kpp to the rLefsHandSideMatrix
    CalculateAndAddKpp( rLeftHandSideMatrix, rVariables, rIntegrationWeight );
    //std::cout<<" Kpp " << rLeftHandSideMatrix-Kh<<std::endl;
    //Kh=rLeftHandSideMatrix;

    // operation performed: add Kpp Stab to the rLefsHandSideMatrix
    CalculateAndAddKppStab( rLeftHandSideMatrix, rVariables, rIntegrationWeight );
    //std::cout<<" KppStab " << rLeftHandSideMatrix-Kh<<std::endl;

    //std::cout<<" rLeftHandSideMatrix " << rLeftHandSideMatrix<<std::endl;
    rVariables.detF     = DeterminantF;
    rVariables.detF0   /= rVariables.detF;
    //KRATOS_WATCH( rLeftHandSideMatrix )


}
//************************************************************************************
//************************************************************************************

void UpdatedLagrangianUP::CalculateAndAddKuum(MatrixType& rLeftHandSideMatrix,
        GeneralVariables& rVariables,
        double& rIntegrationWeight
                                             )
{
    KRATOS_TRY
    //std::stringstream ss;

    //unsigned int number_of_nodes = GetGeometry().size();
    //unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    //unsigned int voigtsize  = 3;
    //unsigned int MatSize = number_of_nodes * dimension;

    //Matrix temp = ZeroMatrix(voigtsize, MatSize);


    //temp = prod( rVariables.ConstitutiveMatrix, rVariables.B );

    MatrixType Kh=rLeftHandSideMatrix;

    Matrix Kuum = prod( trans( rVariables.B ),  rIntegrationWeight * Matrix( prod( rVariables.ConstitutiveMatrix, rVariables.B ) ) );
    //std::cout<< " rVariables.B "<< rVariables.B<<std::endl;
    //std::cout<< " rVariables.ConstitutiveMatrix "<< rVariables.ConstitutiveMatrix<<std::endl;
    //assemble into rk the material uu contribution:
    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    unsigned int dimension = GetGeometry().WorkingSpaceDimension();



    unsigned int indexi = 0;
    unsigned int indexj = 0;
    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        for ( unsigned int idim = 0; idim < dimension ; idim ++)
        {
            indexj=0;
            for ( unsigned int j = 0; j < number_of_nodes; j++ )
            {
                for ( unsigned int jdim = 0; jdim < dimension ; jdim ++)
                {
                    rLeftHandSideMatrix(indexi+i,indexj+j)+=Kuum(indexi,indexj);
                    indexj++;
                }
            }
            indexi++;
        }
    }

    KRATOS_CATCH( "" )
}



//************************************************************************************
//************************************************************************************

void UpdatedLagrangianUP::CalculateAndAddKuug(MatrixType& rLeftHandSideMatrix,
        GeneralVariables& rVariables,
        double& rIntegrationWeight)

{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().size();
    unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    int size = number_of_nodes * dimension;
    MatrixType Kh=rLeftHandSideMatrix;  
    Matrix StressTensor = MathUtils<double>::StressVectorToTensor( rVariables.StressVector );

    Matrix ReducedKg = prod( rVariables.DN_DX, rIntegrationWeight * Matrix( prod( StressTensor, trans( rVariables.DN_DX ) ) ) ); //to be optimized

    Matrix Kuug = zero_matrix<double> (size);
    MathUtils<double>::ExpandAndAddReducedMatrix( Kuug, ReducedKg, dimension );

    unsigned int indexi = 0;
    unsigned int indexj = 0;
    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        for ( unsigned int idim = 0; idim < dimension ; idim ++)
        {
            indexj=0;
            for ( unsigned int j = 0; j < number_of_nodes; j++ )
            {
                for ( unsigned int jdim = 0; jdim < dimension ; jdim ++)
                {
                    rLeftHandSideMatrix(indexi+i,indexj+j)+=Kuug(indexi,indexj);
                    indexj++;
                }
            }
            indexi++;
        }
    }
    
    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void UpdatedLagrangianUP::CalculateAndAddKup (MatrixType& rLeftHandSideMatrix,
        GeneralVariables& rVariables,
        double& rIntegrationWeight)
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    //unsigned int voigtsize  = 3;


    //if( dimension == 3 )
    //{
    //voigtsize  = 6;
    //}

    //here I call a public method of CL class which evaluate the additional term due to the
    //derivative of the deviatoric term respect to the pressure

    //Matrix MixedUPTerm = mConstitutiveLawVector->GetValue(MIXED_UP_TERM, MixedUPTerm );
    //std::cout<<"MixedUPTerm "<<MixedUPTerm<<std::endl;
    //Vector MixedUPTermVector = MathUtils<double>::StressTensorToVector( MixedUPTerm, voigtsize );

    MatrixType Kh=rLeftHandSideMatrix;
    //contributions to stiffness matrix calculated on the reference configuration
    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        //double A1 = prod(trans( rVariables.DN_DX( i , k ) ), MixedUPTermVector);
        unsigned int indexp  = dimension;
        unsigned int indexup = dimension * i + i;
        for ( unsigned int j = 0; j < number_of_nodes; j++ )
        {

            for ( unsigned int k = 0; k < dimension; k++ )
            {
                rLeftHandSideMatrix(indexup+k,indexp) +=  rVariables.DN_DX ( i, k ) *  rVariables.N[j] * rIntegrationWeight * rVariables.detF;
                //rLeftHandSideMatrix(indexup+k,indexp) +=  MixedUPTerm(i,k) * rVariables.DN_DX ( i , k ) *  rVariables.N[j] * rIntegrationWeight;

            }
            indexp += (dimension + 1);
        }
    }

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void UpdatedLagrangianUP::CalculateAndAddKpu (MatrixType& rLeftHandSideMatrix,
        GeneralVariables& rVariables,
        double& rIntegrationWeight)

{
    KRATOS_TRY

    //repasar

    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    MatrixType Kh=rLeftHandSideMatrix;

    //contributions to stiffness matrix calculated on the reference configuration
    unsigned int indexp = dimension;

    //***********************************************************************************
    //to make symmetric the system I divide by U'' all the other terms related to the II eq

    //double DeltaCoefficient = 0;
    //DeltaCoefficient = this->CalculatePUDeltaCoefficient( DeltaCoefficient, rVariables );


    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        for ( unsigned int j = 0; j < number_of_nodes; j++ )
        {
            int indexup= dimension*j + j;
            for ( unsigned int k = 0; k < dimension; k++ )
            {
                rLeftHandSideMatrix(indexp,indexup+k) +=  rVariables.N[i] * rVariables.DN_DX ( j, k ) * rIntegrationWeight * rVariables.detF;

                //std::cout<<" value ("<<indexp<<","<<indexup+k<<") "<<(2*detF) * rN[i] * rDN_DX ( j , k ) * rIntegrationWeight<<std::endl;
            }
        }
        indexp += (dimension + 1);
    }

    KRATOS_CATCH( "" )
}
//************************************************************************************
//************************************************************************************

void UpdatedLagrangianUP::CalculateAndAddKpp (MatrixType& rLeftHandSideMatrix,
        GeneralVariables& rVariables,
        double& rIntegrationWeight)
{
    KRATOS_TRY


    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    double BulkModulus= GetProperties()[YOUNG_MODULUS]/(3*(1-2*GetProperties()[POISSON_RATIO]));

    double DeltaCoefficient = 0;
    DeltaCoefficient = this->CalculatePUDeltaCoefficient( DeltaCoefficient, rVariables );

    MatrixType Kh=rLeftHandSideMatrix;

    //contributions to stiffness matrix calculated on the reference configuration
    unsigned int indexpi = dimension;
    //double consistent = 1.0;

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        unsigned int indexpj = dimension;
        for ( unsigned int j = 0; j < number_of_nodes; j++ )
        {

            // consistent=1;
            // if(indexpi==indexpj)
            //   consistent=2;

            // if( dimension == 2 ){ //consistent 2D

            //   rLeftHandSideMatrix(indexpi,indexpj)  -= consistent * ((1.0)/(BulkModulus)) * (1.0/12.0) * rIntegrationWeight / (rVariables.detF0/rVariables.detF); //2D

            // }
            // else{

            //   rLeftHandSideMatrix(indexpi,indexpj)  -= consistent * ((1.0)/(BulkModulus)) * (1.0/20.0) * rIntegrationWeight / (rVariables.detF0/rVariables.detF); //3D

            // }

            rLeftHandSideMatrix(indexpi,indexpj)  -= ((1.0)/(BulkModulus)) * rVariables.N[i] * rVariables.N[j] * rIntegrationWeight /(DeltaCoefficient * (rVariables.detF0/rVariables.detF)); //2D-3D

            //********************************************************
            //mixed formulation considering a variation of p: 2nd term
            //rLeftHandSideMatrix(indexpi,indexpj)  -= rVariables.N[i] * rVariables.N[j] * rIntegrationWeight / (rVariables.detF0/rVariables.detF);

            indexpj += (dimension + 1);
        }

        indexpi += (dimension + 1);
    }


    //for ( unsigned int i = 0; i < number_of_nodes; i++ )
    //{
    //unsigned int indexpj = dimension;
    //for ( unsigned int j = 0; j < number_of_nodes; j++ )
    //{
    //rLeftHandSideMatrix(indexpj,indexpi)  -= ((1.0)/(BulkModulus)) * rVariables.N[i] * rIntegrationWeight / (rVariables.detF0/rVariables.detF);

    //indexpj += (dimension + 1);

    //}

    //indexpi += (dimension + 1);
    //}

    // std::cout<<std::endl;
    //std::cout<<" Kpp "<<rLeftHandSideMatrix-Kh<<std::endl;

    KRATOS_CATCH( "" )
}
//************************************************************************************
//************************************************************************************
// I changed the constant matrix in the stabilized term:
// as in MPM the position of the integration points does not coincide with the
// position of the Gauss points the first matrix is substitute with the product of the
// shape function values of each integration point
//************************************************************************************
//************************************************************************************

void UpdatedLagrangianUP::CalculateAndAddKppStab (MatrixType& rLeftHandSideMatrix,
        GeneralVariables & rVariables,
        double& rIntegrationWeight)
{
    KRATOS_TRY

    //repasar

    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    double DeltaCoefficient = 0;
    DeltaCoefficient = this->CalculatePUDeltaCoefficient( DeltaCoefficient, rVariables );

    MatrixType Kh=rLeftHandSideMatrix;

    //contributions to stiffness matrix calculated on the reference configuration
    unsigned int indexpi = dimension;

    double AlphaStabilization  = 1.0;
    double StabilizationFactor = GetProperties()[STABILIZATION_FACTOR];
    AlphaStabilization *= StabilizationFactor;

    const double& YoungModulus          = GetProperties()[YOUNG_MODULUS];
    const double& PoissonCoefficient    = GetProperties()[POISSON_RATIO];

    double LameMu =  YoungModulus/(2*(1+PoissonCoefficient));

    //Experimental
    // if(LameMu < rVariables.ConstitutiveMatrix(2,2))
    //   LameMu = rVariables.ConstitutiveMatrix(2,2);

    double consistent = 1.0;

    double FactorValue = 8.0; //JMR deffault value
    if( dimension == 3 )
        FactorValue = 10.0; //JMC deffault value

    //consistent = AlphaStabilization*FactorValue/LameMu;

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        unsigned int indexpj = dimension;
        for ( unsigned int j = 0; j < number_of_nodes; j++ )
        {

            if( dimension == 2 )  //consistent 2D
            {

                consistent=(-1)*AlphaStabilization*FactorValue/(36.0*LameMu);
                if(indexpi==indexpj)
                    consistent=2*AlphaStabilization*FactorValue/(36.0*LameMu);

                //rLeftHandSideMatrix(indexpi,indexpj) -= consistent *(rVariables.N[i]*rVariables.N[j] -1/9)*rIntegrationWeight / (rVariables.detF0/rVariables.detF); //2D
                rLeftHandSideMatrix(indexpi,indexpj) -= consistent *rIntegrationWeight / (DeltaCoefficient * (rVariables.detF0/rVariables.detF)); //2D

            }
            else
            {

                consistent=(-1)*AlphaStabilization*FactorValue/(80.0*LameMu);
                if(indexpi==indexpj)
                    consistent=3*AlphaStabilization*FactorValue/(80.0*LameMu);

                rLeftHandSideMatrix(indexpi,indexpj) -= consistent * rIntegrationWeight / (rVariables.detF0/rVariables.detF); //3D

            }


            indexpj += (dimension + 1);
        }

        indexpi += (dimension + 1);
    }



    //for ( unsigned int i = 0; i < number_of_nodes; i++ )
    //{
    //unsigned int indexpj = dimension;
    //for ( unsigned int j = 0; j < number_of_nodes; j++ )
    //{
    //rLeftHandSideMatrix(indexpj,indexpi) -= consistent *(rVariables.N[i]* -1/3)*rIntegrationWeight / (rVariables.detF0/rVariables.detF);


    //indexpj += (dimension + 1);
    //}

    //indexpi += (dimension + 1);
    //}

    // std::cout<<std::endl;
    //std::cout<<" KppStab "<<rLeftHandSideMatrix-Kh<<std::endl;

    KRATOS_CATCH( "" )
}


//************************************CALCULATE VOLUME CHANGE*************************
//************************************************************************************

double& UpdatedLagrangianUP::CalculateVolumeChange( double& rVolumeChange, GeneralVariables& rVariables )
{
    KRATOS_TRY

    rVolumeChange = 1.0 / (rVariables.detF * rVariables.detF0);

    return rVolumeChange;

    KRATOS_CATCH( "" )
}
//************************************CALCULATE VOLUME ACCELERATION*******************
//************************************************************************************

//Vector& UpdatedLagrangian::CalculateVolumeForce( Vector& rVolumeForce, GeneralVariables& rVariables )
//{
//KRATOS_TRY

//const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

//rVolumeForce = ZeroVector(dimension);


//rVolumeForce = this->GetValue(MP_VOLUME_ACCELERATION)* this->GetValue(MP_MASS);


//return rVolumeForce;

//KRATOS_CATCH( "" )
//}
//************************************************************************************
//************************************************************************************
//void UpdatedLagrangian::CalculateRightHandSide( VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo )
//{
////create local system components
//LocalSystemComponents LocalSystem;

////calculation flags
//LocalSystem.CalculationFlags.Set(UpdatedLagrangian::COMPUTE_RHS_VECTOR);

//MatrixType LeftHandSideMatrix = Matrix();

////Initialize sizes for the system components:
//this->InitializeSystemMatrices( LeftHandSideMatrix, rRightHandSideVector, LocalSystem.CalculationFlags );

////Set Variables to Local system components
//LocalSystem.SetLeftHandSideMatrix(LeftHandSideMatrix);
//LocalSystem.SetRightHandSideVector(rRightHandSideVector);

////Calculate elemental system
//CalculateElementalSystem( LocalSystem, rCurrentProcessInfo );
//}

//************************************************************************************
//************************************************************************************


//void UpdatedLagrangian::CalculateRightHandSide( std::vector< VectorType >& rRightHandSideVectors, const std::vector< Variable< VectorType > >& rRHSVariables, ProcessInfo& rCurrentProcessInfo )
//{
////create local system components
//LocalSystemComponents LocalSystem;

////calculation flags
//LocalSystem.CalculationFlags.Set(UpdatedLagrangian::COMPUTE_RHS_VECTOR);
//LocalSystem.CalculationFlags.Set(UpdatedLagrangian::COMPUTE_RHS_VECTOR_WITH_COMPONENTS);

//MatrixType LeftHandSideMatrix = Matrix();

////Initialize sizes for the system components:
//if( rRHSVariables.size() != rRightHandSideVectors.size() )
//rRightHandSideVectors.resize(rRHSVariables.size());

//for( unsigned int i=0; i<rRightHandSideVectors.size(); i++ )
//{
//this->InitializeSystemMatrices( LeftHandSideMatrix, rRightHandSideVectors[i], LocalSystem.CalculationFlags );
//}

////Set Variables to Local system components
//LocalSystem.SetLeftHandSideMatrix(LeftHandSideMatrix);
//LocalSystem.SetRightHandSideVectors(rRightHandSideVectors);

//LocalSystem.SetRightHandSideVariables(rRHSVariables);

////Calculate elemental system
//CalculateElementalSystem( LocalSystem, rCurrentProcessInfo );
//}


//************************************************************************************
//************************************************************************************


//void UpdatedLagrangian::CalculateLeftHandSide( MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo )
//{
////create local system components
//LocalSystemComponents LocalSystem;

////calculation flags
//LocalSystem.CalculationFlags.Set(UpdatedLagrangian::COMPUTE_LHS_MATRIX);

//VectorType RightHandSideVector = Vector();

////Initialize sizes for the system components:
//this->InitializeSystemMatrices( rLeftHandSideMatrix, RightHandSideVector, LocalSystem.CalculationFlags );

////Set Variables to Local system components
//LocalSystem.SetLeftHandSideMatrix(rLeftHandSideMatrix);
//LocalSystem.SetRightHandSideVector(RightHandSideVector);

////Calculate elemental system
//CalculateElementalSystem( LocalSystem, rCurrentProcessInfo );

//}
//************************************************************************************
//************************************************************************************


//void UpdatedLagrangian::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo )
//{

////create local system components
//LocalSystemComponents LocalSystem;

////calculation flags
//LocalSystem.CalculationFlags.Set(UpdatedLagrangian::COMPUTE_LHS_MATRIX);
//LocalSystem.CalculationFlags.Set(UpdatedLagrangian::COMPUTE_RHS_VECTOR);

////Initialize sizes for the system components:
//this->InitializeSystemMatrices( rLeftHandSideMatrix, rRightHandSideVector, LocalSystem.CalculationFlags );

////Set Variables to Local system components
//LocalSystem.SetLeftHandSideMatrix(rLeftHandSideMatrix);
//LocalSystem.SetRightHandSideVector(rRightHandSideVector);

////Calculate elemental system

//CalculateElementalSystem( LocalSystem, rCurrentProcessInfo );
////std::cout<<" in CalculateLocalSystem ends"<<std::endl;


//}


//************************************************************************************
//************************************************************************************

//void UpdatedLagrangian::CalculateLocalSystem( std::vector< MatrixType >& rLeftHandSideMatrices,
//const std::vector< Variable< MatrixType > >& rLHSVariables,
//std::vector< VectorType >& rRightHandSideVectors,
//const std::vector< Variable< VectorType > >& rRHSVariables,
//ProcessInfo& rCurrentProcessInfo )
//{
////create local system components
//LocalSystemComponents LocalSystem;

////calculation flags
//LocalSystem.CalculationFlags.Set(UpdatedLagrangian::COMPUTE_LHS_MATRIX_WITH_COMPONENTS);
//LocalSystem.CalculationFlags.Set(UpdatedLagrangian::COMPUTE_RHS_VECTOR_WITH_COMPONENTS);


////Initialize sizes for the system components:
//if( rLHSVariables.size() != rLeftHandSideMatrices.size() )
//rLeftHandSideMatrices.resize(rLHSVariables.size());

//if( rRHSVariables.size() != rRightHandSideVectors.size() )
//rRightHandSideVectors.resize(rRHSVariables.size());

//LocalSystem.CalculationFlags.Set(UpdatedLagrangian::COMPUTE_LHS_MATRIX);
//for( unsigned int i=0; i<rLeftHandSideMatrices.size(); i++ )
//{
////Note: rRightHandSideVectors.size() > 0
//this->InitializeSystemMatrices( rLeftHandSideMatrices[i], rRightHandSideVectors[0], LocalSystem.CalculationFlags );
//}

//LocalSystem.CalculationFlags.Set(UpdatedLagrangian::COMPUTE_RHS_VECTOR);
//LocalSystem.CalculationFlags.Set(UpdatedLagrangian::COMPUTE_LHS_MATRIX,false);

//for( unsigned int i=0; i<rRightHandSideVectors.size(); i++ )
//{
////Note: rLeftHandSideMatrices.size() > 0
//this->InitializeSystemMatrices( rLeftHandSideMatrices[0], rRightHandSideVectors[i], LocalSystem.CalculationFlags );
//}
//LocalSystem.CalculationFlags.Set(UpdatedLagrangian::COMPUTE_LHS_MATRIX,true);


////Set Variables to Local system components
//LocalSystem.SetLeftHandSideMatrices(rLeftHandSideMatrices);
//LocalSystem.SetRightHandSideVectors(rRightHandSideVectors);

//LocalSystem.SetLeftHandSideVariables(rLHSVariables);
//LocalSystem.SetRightHandSideVariables(rRHSVariables);

////Calculate elemental system
//CalculateElementalSystem( LocalSystem, rCurrentProcessInfo );

//}


////************************************************************************************
////************************************************************************************

//void UpdatedLagrangian::InitializeSolutionStep( ProcessInfo& rCurrentProcessInfo )
//{
//// In the Initialize of each time step the nodal initial conditions are evaluated
////1. first of all I need to evaluate the MP momentum and MP_inertia



//int MP_bool = this->GetValue(MP_BOOL);

////std::cout<<" in InitializeSolutionStep2"<<std::endl;
//unsigned int dimension = GetGeometry().WorkingSpaceDimension();
//const unsigned int number_of_nodes = GetGeometry().PointsNumber();
//array_1d<double,3>& xg = this->GetValue(GAUSS_COORD);
//GeneralVariables Variables;
////this->InitializeGeneralVariables(Variables,rCurrentProcessInfo);



//Matrix J0 = ZeroMatrix(dimension, dimension);

//J0 = this->MPMJacobian(J0, xg);

////calculating and storing inverse and the determinant of the jacobian
//MathUtils<double>::InvertMatrix( J0, mInverseJ0, mDeterminantJ0 );

//Variables.N = this->MPMShapeFunctionPointValues(Variables.N, xg);

//mConstitutiveLawVector->InitializeSolutionStep( GetProperties(),
//GetGeometry(), Variables.N, rCurrentProcessInfo );

//mFinalizedStep = false;



//array_1d<double,3>& MP_Velocity = this->GetValue(MP_VELOCITY);
//array_1d<double,3>& MP_Acceleration = this->GetValue(MP_ACCELERATION);
//array_1d<double,3>& AUX_MP_Velocity = this->GetValue(AUX_MP_VELOCITY);
//array_1d<double,3>& AUX_MP_Acceleration = this->GetValue(AUX_MP_ACCELERATION);
//double MP_Mass = this->GetValue(MP_MASS);
//array_1d<double,3> MP_Momentum;
//array_1d<double,3> MP_Inertia;
//array_1d<double,3> NodalMomentum;
//array_1d<double,3> NodalInertia;

//for (unsigned int j=0;j<number_of_nodes;j++)
//{
////these are the values of nodal velocity and nodal acceleration evaluated in the initialize solution step
//array_1d<double, 3 > & NodalAcceleration = GetGeometry()[j].FastGetSolutionStepValue(ACCELERATION,1);
//array_1d<double, 3 > & NodalVelocity = GetGeometry()[j].FastGetSolutionStepValue(VELOCITY,1);

////std::cout<<"NodalVelocity "<< GetGeometry()[j].Id()<<std::endl;
//for (unsigned int k = 0; k < dimension; k++)
//{
//AUX_MP_Velocity[k] += Variables.N[j] * NodalVelocity[k];
//AUX_MP_Acceleration[k] += Variables.N[j] * NodalAcceleration[k];
//}
//}
////std::cout<<" AUX_MP_Velocity "<<AUX_MP_Velocity<<std::endl;
////std::cout<<" AUX_MP_Acceleration "<<AUX_MP_Acceleration<<std::endl;


////for (unsigned int k = 0; k < 3; k++)
////{
////MP_Momentum[k] = (MP_Velocity[k] - AUX_MP_Velocity[k]) * MP_Mass;
////MP_Inertia[k] = (MP_Acceleration[k] - AUX_MP_Acceleration[k]) * MP_Mass;
////}

//// Here MP contribution in terms of momentum, inertia and mass are added
//for ( unsigned int i = 0; i < number_of_nodes; i++ )
//{
//for (unsigned int j = 0; j < dimension; j++)
//{
//NodalMomentum[j] = Variables.N[i] * (MP_Velocity[j] - AUX_MP_Velocity[j]) * MP_Mass;
//NodalInertia[j] = Variables.N[i] * (MP_Acceleration[j] - AUX_MP_Acceleration[j]) * MP_Mass;

//}
//GetGeometry()[i].GetSolutionStepValue(NODAL_MOMENTUM, 0) += NodalMomentum;
//GetGeometry()[i].GetSolutionStepValue(NODAL_INERTIA, 0) += NodalInertia;
////if(GetGeometry()[i].pGetDof(DISPLACEMENT_X)->IsFixed() == false)
////{
//GetGeometry()[i].GetSolutionStepValue(NODAL_MASS, 0) += Variables.N[i] * MP_Mass;
////}
//}
////for ( unsigned int i = 0; i < number_of_nodes; i++ )

////{
////if(GetGeometry()[i].pGetDof(DISPLACEMENT_X)->IsFixed() == false)
////{
////for (unsigned int j = 0; j < number_of_nodes; j++)
////{
////if(GetGeometry()[j].pGetDof(DISPLACEMENT_X)->IsFixed() == false)
////{
////GetGeometry()[i].GetSolutionStepValue(NODAL_MASS, 0) += Variables.N[i] * Variables.N[j] * MP_Mass;
////}
////}
////}
////}
////MP_Velocity = MP_Velocity - AUX_MP_Velocity;
////MP_Acceleration = MP_Acceleration - AUX_MP_Acceleration;
//AUX_MP_Velocity.clear();
//AUX_MP_Acceleration.clear();



//if(MP_bool == 0)
//{
////std::cout<<" in InitializeSolutionStep1"<<std::endl;
//unsigned int dimension = GetGeometry().WorkingSpaceDimension();
//const unsigned int number_of_nodes = GetGeometry().PointsNumber();
//array_1d<double,3>& xg = this->GetValue(GAUSS_COORD);
//GeneralVariables Variables;
////this->SetValue(PREVIOUS_MP_CAUCHY_STRESS_VECTOR, Variables.StressVector);
////this->SetValue(PREVIOUS_MP_ALMANSI_STRAIN_VECTOR, Variables.StrainVector);
////this->InitializeGeneralVariables(Variables,rCurrentProcessInfo);

//Matrix J0 = ZeroMatrix(dimension, dimension);

//J0 = this->MPMJacobian(J0, xg);

////calculating and storing inverse and the determinant of the jacobian
//MathUtils<double>::InvertMatrix( J0, mInverseJ0, mDeterminantJ0 );

//Variables.N = this->MPMShapeFunctionPointValues(Variables.N, xg);

//mConstitutiveLawVector->InitializeSolutionStep( GetProperties(),
//GetGeometry(), Variables.N, rCurrentProcessInfo );

//mFinalizedStep = false;



//array_1d<double,3>& MP_Velocity = this->GetValue(MP_VELOCITY);
//array_1d<double,3>& MP_Acceleration = this->GetValue(MP_ACCELERATION);

//double MP_Mass = this->GetValue(MP_MASS);
//array_1d<double,3> MP_Momentum;
//array_1d<double,3> MP_Inertia;
//array_1d<double,3> NodalMomentum;
//array_1d<double,3> NodalInertia;




//for (unsigned int k = 0; k < dimension; k++)
//{
//MP_Momentum[k] = MP_Velocity[k] * MP_Mass;
//MP_Inertia[k] = MP_Acceleration[k] * MP_Mass;
//}
////std::cout<<"MP_Momentum "<<MP_Momentum<<std::endl;
////std::cout<<"MP_Inertia "<<MP_Inertia<<std::endl;
//// Here MP contribution in terms of momentum, inertia and mass are added

//for ( unsigned int i = 0; i < number_of_nodes; i++ )
//{
//for (unsigned int j = 0; j < dimension; j++)
//{
//NodalMomentum[j] = Variables.N[i] * MP_Momentum[j];
//NodalInertia[j] = Variables.N[i] * MP_Inertia[j];

//}
//GetGeometry()[i].GetSolutionStepValue(NODAL_MOMENTUM, 0) += NodalMomentum;
//GetGeometry()[i].GetSolutionStepValue(NODAL_INERTIA, 0) += NodalInertia;
////if(GetGeometry()[i].pGetDof(DISPLACEMENT_X)->IsFixed() == false)
////{
//GetGeometry()[i].GetSolutionStepValue(NODAL_MASS, 0) += Variables.N[i] * MP_Mass;
////}
//}
////for ( unsigned int i = 0; i < number_of_nodes; i++ )

////{
////if(GetGeometry()[i].pGetDof(DISPLACEMENT_X)->IsFixed() == false)
////{
////for (unsigned int j = 0; j < number_of_nodes; j++)
////{
////if(GetGeometry()[j].pGetDof(DISPLACEMENT_X)->IsFixed() == false)
////{
////GetGeometry()[i].GetSolutionStepValue(NODAL_MASS, 0) += Variables.N[i] * Variables.N[j] * MP_Mass;
////}
////}
////}
////}



//this->SetValue(MP_BOOL,1);
//}
//else
//{
////std::cout<<" in InitializeSolutionStep2"<<std::endl;
//unsigned int dimension = GetGeometry().WorkingSpaceDimension();
//const unsigned int number_of_nodes = GetGeometry().PointsNumber();
//array_1d<double,3>& xg = this->GetValue(GAUSS_COORD);
//GeneralVariables Variables;
////this->InitializeGeneralVariables(Variables,rCurrentProcessInfo);



//Matrix J0 = ZeroMatrix(dimension, dimension);

//J0 = this->MPMJacobian(J0, xg);

////calculating and storing inverse and the determinant of the jacobian
//MathUtils<double>::InvertMatrix( J0, mInverseJ0, mDeterminantJ0 );

//Variables.N = this->MPMShapeFunctionPointValues(Variables.N, xg);

//mConstitutiveLawVector->InitializeSolutionStep( GetProperties(),
//GetGeometry(), Variables.N, rCurrentProcessInfo );

//mFinalizedStep = false;



//array_1d<double,3>& MP_Velocity = this->GetValue(MP_VELOCITY);
//array_1d<double,3>& MP_Acceleration = this->GetValue(MP_ACCELERATION);
//array_1d<double,3>& AUX_MP_Velocity = this->GetValue(AUX_MP_VELOCITY);
//array_1d<double,3>& AUX_MP_Acceleration = this->GetValue(AUX_MP_ACCELERATION);
//double MP_Mass = this->GetValue(MP_MASS);
//array_1d<double,3> MP_Momentum;
//array_1d<double,3> MP_Inertia;
//array_1d<double,3> NodalMomentum;
//array_1d<double,3> NodalInertia;

//for (unsigned int j=0;j<number_of_nodes;j++)
//{
////these are the values of nodal velocity and nodal acceleration evaluated in the initialize solution step
//array_1d<double, 3 > & NodalAcceleration = GetGeometry()[j].FastGetSolutionStepValue(ACCELERATION,1);
//array_1d<double, 3 > & NodalVelocity = GetGeometry()[j].FastGetSolutionStepValue(VELOCITY,1);

////std::cout<<"NodalVelocity "<< GetGeometry()[j].Id()<<std::endl;
//for (unsigned int k = 0; k < dimension; k++)
//{
//AUX_MP_Velocity[k] += Variables.N[j] * NodalVelocity[k];
//AUX_MP_Acceleration[k] += Variables.N[j] * NodalAcceleration[k];
//}
//}
////std::cout<<" AUX_MP_Velocity "<<AUX_MP_Velocity<<std::endl;
////std::cout<<" AUX_MP_Acceleration "<<AUX_MP_Acceleration<<std::endl;


////for (unsigned int k = 0; k < 3; k++)
////{
////MP_Momentum[k] = (MP_Velocity[k] - AUX_MP_Velocity[k]) * MP_Mass;
////MP_Inertia[k] = (MP_Acceleration[k] - AUX_MP_Acceleration[k]) * MP_Mass;
////}

//// Here MP contribution in terms of momentum, inertia and mass are added
//for ( unsigned int i = 0; i < number_of_nodes; i++ )
//{
//for (unsigned int j = 0; j < dimension; j++)
//{
//NodalMomentum[j] = Variables.N[i] * (MP_Velocity[j] - AUX_MP_Velocity[j]) * MP_Mass;
//NodalInertia[j] = Variables.N[i] * (MP_Acceleration[j] - AUX_MP_Acceleration[j]) * MP_Mass;

//}
//GetGeometry()[i].GetSolutionStepValue(NODAL_MOMENTUM, 0) += NodalMomentum;
//GetGeometry()[i].GetSolutionStepValue(NODAL_INERTIA, 0) += NodalInertia;
////if(GetGeometry()[i].pGetDof(DISPLACEMENT_X)->IsFixed() == false)
////{
//GetGeometry()[i].GetSolutionStepValue(NODAL_MASS, 0) += Variables.N[i] * MP_Mass;
////}
//}
////for ( unsigned int i = 0; i < number_of_nodes; i++ )

////{
////if(GetGeometry()[i].pGetDof(DISPLACEMENT_X)->IsFixed() == false)
////{
////for (unsigned int j = 0; j < number_of_nodes; j++)
////{
////if(GetGeometry()[j].pGetDof(DISPLACEMENT_X)->IsFixed() == false)
////{
////GetGeometry()[i].GetSolutionStepValue(NODAL_MASS, 0) += Variables.N[i] * Variables.N[j] * MP_Mass;
////}
////}
////}
////}
////MP_Velocity = MP_Velocity - AUX_MP_Velocity;
////MP_Acceleration = MP_Acceleration - AUX_MP_Acceleration;
//AUX_MP_Velocity.clear();
//AUX_MP_Acceleration.clear();
//}


//}

////************************************************************************************
////************************************************************************************


////************************************************************************************
////************************************************************************************
//void UpdatedLagrangianUP::InitializeNonLinearIteration( ProcessInfo& rCurrentProcessInfo )
//{



//}

//////************************************************************************************
//////************************************************************************************

//void UpdatedLagrangian::FinalizeNonLinearIteration( ProcessInfo& rCurrentProcessInfo )
//{

//}

////************************************************************************************
////************************************************************************************

//void UpdatedLagrangian::FinalizeSolutionStep( ProcessInfo& rCurrentProcessInfo )
//{
//KRATOS_TRY

//const unsigned int number_of_nodes = GetGeometry().PointsNumber();
//const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
//unsigned int voigtsize  = 3;

//if( dimension == 3 )
//{
//voigtsize  = 6;
//}
////Vector NodalStress = ZeroVector(voigtsize);
//array_1d<double,3> NodalStress;
////create and initialize element variables:
//GeneralVariables Variables;
//this->InitializeGeneralVariables(Variables,rCurrentProcessInfo);

////create constitutive law parameters:
//ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

////set constitutive law flags:
//Flags &ConstitutiveLawOptions=Values.GetOptions();

//ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);
//ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS);
////ConstitutiveLawOptions.Set(ConstitutiveLaw::VOLUMETRIC_TENSOR_ONLY);
////compute element kinematics B, F, DN_DX ...
//this->CalculateKinematics(Variables, rCurrentProcessInfo);

////set general variables to constitutivelaw parameters
//this->SetGeneralVariables(Variables,Values);

////call the constitutive law to update material variables
//mConstitutiveLawVector->FinalizeMaterialResponse(Values, Variables.StressMeasure);

////call the constitutive law to finalize the solution step
//mConstitutiveLawVector->FinalizeSolutionStep( GetProperties(),
//GetGeometry(),
//Variables.N,
//rCurrentProcessInfo );
//for ( unsigned int i = 0; i < number_of_nodes; i++ )
//{
////GetGeometry()[i].GetSolutionStepValue(STRESSES, 0) = ZeroVector(voigtsize);
//if(GetGeometry()[i].GetSolutionStepValue(NODAL_MASS, 0) > 1e-4)
//{
//for ( unsigned int j = 0; j< voigtsize; j++)
//{
//NodalStress[j] = Variables.N[i] * Variables.StressVector[j] * this->GetValue(MP_MASS) / GetGeometry()[i].GetSolutionStepValue(NODAL_MASS, 0);

//}
//GetGeometry()[i].GetSolutionStepValue(NODAL_STRESSES, 0) += NodalStress;
////std::cout<<" GetGeometry()[i].GetSolutionStepValue(STRESSES, 0) "<<  GetGeometry()[i].GetSolutionStepValue(STRESSES, 0)<<std::endl;
//}
//}
////call the element internal variables update
//this->FinalizeStepVariables(Variables, rCurrentProcessInfo);


//mFinalizedStep = true;

//KRATOS_CATCH( "" )
//}

////************************************************************************************************************
//void UpdatedLagrangian::Calculate(const Variable<Vector >& rVariable, Vector& Output, const ProcessInfo& rCurrentProcessInfo)
//{
//KRATOS_TRY


//GeneralVariables Variables;
//this->InitializeGeneralVariables(Variables,rCurrentProcessInfo);
//if (rVariable == MP_CAUCHY_STRESS_VECTOR)
//{
//const unsigned int number_of_nodes = GetGeometry().PointsNumber();
//const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
//unsigned int voigtsize  = 3;


//if( dimension == 3 )
//{
//voigtsize  = 6;
//}
////Vector SmoothMPStress = ZeroVector(voigtsize);
//array_1d<double,3> SmoothMPStress;
//SmoothMPStress.clear();
////Vector NodalStress0 = GetGeometry()[0].GetSolutionStepValue(STRESSES, 0);
////Vector NodalStress1 = GetGeometry()[1].GetSolutionStepValue(STRESSES, 0);
////Vector NodalStress2 = GetGeometry()[2].GetSolutionStepValue(STRESSES, 0);
////Vector NodalStress3 = GetGeometry()[3].GetSolutionStepValue(STRESSES, 0);

////SmoothMPStress = NodalStress0 * Variables.N[0] + NodalStress1 * Variables.N[1] + NodalStress2 * Variables.N[2] + NodalStress3 * Variables.N[3];
//for ( unsigned int i = 0; i < number_of_nodes; i++ )
//{
//array_1d<double,3> NodalStress = GetGeometry()[i].GetSolutionStepValue(NODAL_STRESSES, 0);
//for ( unsigned int j = 0; j< voigtsize; j++)
//{
//SmoothMPStress[j] += NodalStress[j] * Variables.N[i];
//}
//}
//this->SetValue(MP_CAUCHY_STRESS_VECTOR, SmoothMPStress);

//}
//KRATOS_CATCH( "" )
//}

////************************************************************************************
////************************************************************************************

//void UpdatedLagrangian::FinalizeStepVariables( GeneralVariables & rVariables, const ProcessInfo& rCurrentProcessInfo)
//{
////update internal (historical) variables
//mDeterminantF0         = rVariables.detF* rVariables.detF0;
//mDeformationGradientF0 = prod(rVariables.F, rVariables.F0);

//this->SetValue(MP_CAUCHY_STRESS_VECTOR, rVariables.StressVector);


//this->SetValue(MP_ALMANSI_STRAIN_VECTOR, rVariables.StrainVector);
////std::cout<<" before get equivalent plastic strain "<<std::endl;
//double EquivalentPlasticStrain = mConstitutiveLawVector->GetValue(PLASTIC_STRAIN, EquivalentPlasticStrain );
////std::cout<<" EquivalentPlasticStrain in the element "<<EquivalentPlasticStrain<<std::endl;
//this->SetValue(MP_EQUIVALENT_PLASTIC_STRAIN, EquivalentPlasticStrain);



//double DeltaEquivalentPlasticStrain = mConstitutiveLawVector->GetValue(DELTA_PLASTIC_STRAIN, DeltaEquivalentPlasticStrain );
//this->SetValue(MP_EQUIVALENT_DELTA_PLASTIC_STRAIN, DeltaEquivalentPlasticStrain);

//ComparisonUtilities EquivalentStress;
//double MPMStressNorm = EquivalentStress.CalculateStressNorm(rVariables.StressVector);
////std::cout<<" MPMStressNorm "<<MPMStressNorm<<std::endl;
//this->SetValue(MPM_NORM_ISOCHORIC_STRESS, MPMStressNorm);

//MathUtils<double>::InvertMatrix( rVariables.j, mInverseJ, rVariables.detJ );
//this->UpdateGaussPoint(rVariables, rCurrentProcessInfo);

//}

//************************************************************************************
//************************************************************************************
///**
//* The position of the Gauss points/Material points is updated
//*/

//void UpdatedLagrangian::UpdateGaussPoint( GeneralVariables & rVariables, const ProcessInfo& rCurrentProcessInfo)
//{
//KRATOS_TRY

//rVariables.CurrentDisp = CalculateCurrentDisp(rVariables.CurrentDisp, rCurrentProcessInfo);
//const unsigned int number_of_nodes = GetGeometry().PointsNumber();
//unsigned int dimension = GetGeometry().WorkingSpaceDimension();

//array_1d<double,3> xg = this->GetValue(GAUSS_COORD);
//array_1d<double,3> MP_PreviousAcceleration = this->GetValue(MP_ACCELERATION);
//array_1d<double,3> MP_PreviousVelocity = this->GetValue(MP_VELOCITY);
//double MP_Mass = this->GetValue(MP_MASS);
//array_1d<double,3> delta_xg = ZeroVector(3);
//array_1d<double,3> MP_Acceleration = ZeroVector(3);
//array_1d<double,3> MP_Velocity = ZeroVector(3);
//double DeltaTime = rCurrentProcessInfo[DELTA_TIME];


//rVariables.N = this->MPMShapeFunctionPointValues(rVariables.N, xg);
//int MP_number = this->GetValue(MP_NUMBER);

////double total_nodal_mass = 0.0;
////for ( unsigned int i = 0; i < number_of_nodes; i++ )
////{
////total_nodal_mass += GetGeometry()[i].GetSolutionStepValue(NODAL_MASS, 0);
////}
//for ( unsigned int i = 0; i < number_of_nodes; i++ )
//{
//if (rVariables.N[i] > 1e-16)
//{
//array_1d<double, 3 > & NodalAcceleration = GetGeometry()[i].FastGetSolutionStepValue(ACCELERATION);
//array_1d<double, 3 > & NodalVelocity = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY);
//double NodalMass = GetGeometry()[i].GetSolutionStepValue(NODAL_MASS, 0);
//array_1d<double,3> NodalMomentum = NodalMass * NodalVelocity;
//array_1d<double,3> NodalInertia = NodalMass * NodalAcceleration;

//for ( unsigned int j = 0; j < dimension; j++ )
//{

//delta_xg[j] += rVariables.N[i] * rVariables.CurrentDisp(i,j);
//MP_Acceleration[j] += rVariables.N[i] * NodalAcceleration[j];
////PERCHE NON FUNZIONA QUESTA INTERPOLAZIONE?
//MP_Velocity[j] += rVariables.N[i] * NodalVelocity[j];


////MP_Acceleration[j] +=NodalInertia[j]/(rVariables.N[i] * MP_Mass * MP_number);//
////MP_Velocity[j] += NodalMomentum[j]/(rVariables.N[i] * MP_Mass * MP_number);
////MP_Velocity[j] += DeltaTime * rVariables.N[i] * NodalAcceleration[j];////




//}
//}

//}

////for ( unsigned int j = 0; j < dimension; j++ )
////{
////MP_Velocity[j] += delta_xg[j]/DeltaTime;

////}

////**************************************************************************************************************************
////Another way to update the MP velocity (see paper Guilkey and Weiss, 2003) !!!USING THIS EXPRESSION I CONSERVE MORE ENERGY
////THIS EXPRESSION IS OK WHEN NO Dirichlet BOUNDARY CONDITIONS ARE APPLIED.
////WHEN BOUNDARY CONDITIONS ARE APPLIED SOMETHING HAS TO BE CHANGED to conserve the total linear momentum
////MP_Velocity = MP_PreviousVelocity + 0.5 * DeltaTime * (MP_Acceleration + MP_PreviousAcceleration);
////MP_Velocity += MP_PreviousVelocity;
////Update the MP Velocity

////MP_Acceleration = 4/(DeltaTime * DeltaTime) * delta_xg - 4/DeltaTime * MP_PreviousVelocity;
////MP_Velocity = 2/DeltaTime * delta_xg - MP_PreviousVelocity;

//this -> SetValue(MP_VELOCITY,MP_Velocity );

//const array_1d<double,3>& new_xg = xg + delta_xg ;

////Update the MP Position
//this -> SetValue(GAUSS_COORD,new_xg);

////Update the MP Acceleration
//this -> SetValue(MP_ACCELERATION,MP_Acceleration);

//array_1d<double,3>& MP_Displacement = this->GetValue(MP_DISPLACEMENT);

//MP_Displacement += delta_xg;

////Update the MP Displacement
//this -> SetValue(MP_DISPLACEMENT,MP_Displacement );

//KRATOS_CATCH( "" )
//}



//void UpdatedLagrangian::InitializeMaterial()
//{
//KRATOS_TRY
//array_1d<double,3>& xg = this->GetValue(GAUSS_COORD);
//GeneralVariables Variables;
////this->InitializeGeneralVariables(Variables,rCurrentProcessInfo);


//if ( GetProperties()[CONSTITUTIVE_LAW] != NULL )
//{

//mConstitutiveLawVector = GetProperties()[CONSTITUTIVE_LAW]->Clone();


//Variables.N = this->MPMShapeFunctionPointValues(Variables.N, xg);

//mConstitutiveLawVector->InitializeMaterial( GetProperties(), GetGeometry(),
//Variables.N );

////}
//}
//else
//KRATOS_THROW_ERROR( std::logic_error, "a constitutive law needs to be specified for the element with ID ", this->Id() )
////std::cout<< "in initialize material "<<std::endl;
//KRATOS_CATCH( "" )
//}


////************************************************************************************
////************************************************************************************

//void UpdatedLagrangian::ResetConstitutiveLaw()
//{
//KRATOS_TRY
//array_1d<double,3>& xg = this->GetValue(GAUSS_COORD);
//GeneralVariables Variables;
////this->InitializeGeneralVariables(Variables,rCurrentProcessInfo);

////create and initialize element variables:

//if ( GetProperties()[CONSTITUTIVE_LAW] != NULL )
//{

//mConstitutiveLawVector->ResetMaterial( GetProperties(), GetGeometry(), this->MPMShapeFunctionPointValues(Variables.N, xg) );
//}

//KRATOS_CATCH( "" )
//}




////*************************COMPUTE CURRENT DISPLACEMENT*******************************
////************************************************************************************


//Matrix& UpdatedLagrangian::CalculateCurrentDisp(Matrix & rCurrentDisp, const ProcessInfo& rCurrentProcessInfo)
//{
//KRATOS_TRY

//const unsigned int number_of_nodes = GetGeometry().PointsNumber();
//unsigned int dimension = GetGeometry().WorkingSpaceDimension();

//rCurrentDisp = zero_matrix<double>( number_of_nodes , dimension);

//for ( unsigned int i = 0; i < number_of_nodes; i++ )
//{

//array_1d<double, 3 > & CurrentDisplacement  = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT);


//for ( unsigned int j = 0; j < dimension; j++ )
//{

//rCurrentDisp(i,j) = CurrentDisplacement[j];
//}
//}

//return rCurrentDisp;

//KRATOS_CATCH( "" )
//}


////*************************COMPUTE ALMANSI STRAIN*************************************
////************************************************************************************
//void UpdatedLagrangian::CalculateAlmansiStrain(const Matrix& rF,
//Vector& rStrainVector )
//{
//KRATOS_TRY

//const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

////Left Cauchy-Green Calculation
//Matrix LeftCauchyGreen = prod( rF, trans( rF ) );

////Calculating the inverse of the jacobian
//Matrix InverseLeftCauchyGreen ( dimension, dimension );
//double det_b=0;
//MathUtils<double>::InvertMatrix( LeftCauchyGreen, InverseLeftCauchyGreen, det_b);

//if( dimension == 2 )
//{

////Almansi Strain Calculation
//rStrainVector[0] = 0.5 * (  1.00 - InverseLeftCauchyGreen( 0, 0 ) );

//rStrainVector[1] = 0.5 * (  1.00 - InverseLeftCauchyGreen( 1, 1 ) );

//rStrainVector[2] = - InverseLeftCauchyGreen( 0, 1 ); // xy

//}
//else if( dimension == 3 )
//{

////Almansi Strain Calculation
//if ( rStrainVector.size() != 6 ) rStrainVector.resize( 6, false );

//rStrainVector[0] = 0.5 * (  1.00 - InverseLeftCauchyGreen( 0, 0 ) );

//rStrainVector[1] = 0.5 * (  1.00 - InverseLeftCauchyGreen( 1, 1 ) );

//rStrainVector[2] = 0.5 * (  1.00 - InverseLeftCauchyGreen( 2, 2 ) );

//rStrainVector[3] = - InverseLeftCauchyGreen( 0, 1 ); // xy

//rStrainVector[4] = - InverseLeftCauchyGreen( 1, 2 ); // yz

//rStrainVector[5] = - InverseLeftCauchyGreen( 0, 2 ); // xz

//}
//else
//{

//KRATOS_THROW_ERROR( std::invalid_argument, "something is wrong with the dimension", "" );

//}


//KRATOS_CATCH( "" )
//}
////*************************COMPUTE GREEN-LAGRANGE STRAIN*************************************
////************************************************************************************

//void UpdatedLagrangian::CalculateGreenLagrangeStrain(const Matrix& rF,
//Vector& rStrainVector )
//{
//KRATOS_TRY

//const unsigned int dimension  = GetGeometry().WorkingSpaceDimension();

////Right Cauchy-Green Calculation
//Matrix C ( dimension, dimension );
//noalias( C ) = prod( trans( rF ), rF );

//if( dimension == 2 )
//{

////Green Lagrange Strain Calculation
//if ( rStrainVector.size() != 3 ) rStrainVector.resize( 3, false );

//rStrainVector[0] = 0.5 * ( C( 0, 0 ) - 1.00 );

//rStrainVector[1] = 0.5 * ( C( 1, 1 ) - 1.00 );

//rStrainVector[2] = C( 0, 1 ); // xy

//}
//else if( dimension == 3 )
//{

////Green Lagrange Strain Calculation
//if ( rStrainVector.size() != 6 ) rStrainVector.resize( 6, false );

//rStrainVector[0] = 0.5 * ( C( 0, 0 ) - 1.00 );

//rStrainVector[1] = 0.5 * ( C( 1, 1 ) - 1.00 );

//rStrainVector[2] = 0.5 * ( C( 2, 2 ) - 1.00 );

//rStrainVector[3] = C( 0, 1 ); // xy

//rStrainVector[4] = C( 1, 2 ); // yz

//rStrainVector[5] = C( 0, 2 ); // xz

//}
//else
//{

//KRATOS_THROW_ERROR( std::invalid_argument, "something is wrong with the dimension", "" )

//}

//KRATOS_CATCH( "" )
//}



//************************************************************************************
//************************************************************************************

//double& UpdatedLagrangian::CalculateIntegrationWeight(double& rIntegrationWeight)
//{
//const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

//if( dimension == 2 )
//rIntegrationWeight *= GetProperties()[THICKNESS];

//return rIntegrationWeight;
//}


//************************************************************************************
//************************************************************************************

void UpdatedLagrangianUP::EquationIdVector( EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo )
{
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
    unsigned int element_size          = number_of_nodes * dimension + number_of_nodes;

    if ( rResult.size() != element_size )
        rResult.resize( element_size, false );

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        int index = i * dimension + i;
        rResult[index]     = GetGeometry()[i].GetDof( DISPLACEMENT_X ).EquationId();
        rResult[index + 1] = GetGeometry()[i].GetDof( DISPLACEMENT_Y ).EquationId();

        if( dimension == 3)
        {
            rResult[index + 2] = GetGeometry()[i].GetDof( DISPLACEMENT_Z ).EquationId();
            rResult[index + 3] = GetGeometry()[i].GetDof( PRESSURE ).EquationId();
        }
        else
        {
            rResult[index + 2] = GetGeometry()[i].GetDof( PRESSURE ).EquationId();
        }

    }

}

//************************************************************************************
//************************************************************************************

void UpdatedLagrangianUP::GetDofList( DofsVectorType& rElementalDofList, ProcessInfo& CurrentProcessInfo )
{
    rElementalDofList.resize( 0 );

    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    for ( unsigned int i = 0; i < GetGeometry().size(); i++ )
    {
        rElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_X ) );
        rElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_Y ) );

        if( dimension == 3 )
            rElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_Z ) );

        rElementalDofList.push_back( GetGeometry()[i].pGetDof( PRESSURE ));
    }
    //std::cout<< "ElementalDofList.size() "<<rElementalDofList.size()<<std::endl;
    //std::cout<<" in GetDofList of derived class"<<std::endl;
}



//************************************************************************************
//*******************DAMPING MATRIX***************************************************

//void UpdatedLagrangian::CalculateDampingMatrix( MatrixType& rDampingMatrix, ProcessInfo& rCurrentProcessInfo )
//{
//KRATOS_TRY

////0.-Initialize the DampingMatrix:
//unsigned int number_of_nodes = GetGeometry().size();
//unsigned int dimension = GetGeometry().WorkingSpaceDimension();

////resizing as needed the LHS
//unsigned int MatSize = number_of_nodes * dimension;

//if ( rDampingMatrix.size1() != MatSize )
//rDampingMatrix.resize( MatSize, MatSize, false );

//noalias( rDampingMatrix ) = ZeroMatrix( MatSize, MatSize );


////1.-Calculate StiffnessMatrix:

//MatrixType StiffnessMatrix  = Matrix();

//this->CalculateLeftHandSide( StiffnessMatrix, rCurrentProcessInfo );

////2.-Calculate MassMatrix:

//MatrixType MassMatrix  = Matrix();

//this->CalculateMassMatrix ( MassMatrix, rCurrentProcessInfo );


////3.-Get Damping Coeffitients (RAYLEIGH_ALPHA, RAYLEIGH_BETA)
//double alpha = 0;
//if( GetProperties().Has(RAYLEIGH_ALPHA) ){
//alpha = GetProperties()[RAYLEIGH_ALPHA];
//}
//else if( rCurrentProcessInfo.Has(RAYLEIGH_ALPHA) ){
//alpha = rCurrentProcessInfo[RAYLEIGH_ALPHA];
//}

//double beta  = 0;
//if( GetProperties().Has(RAYLEIGH_BETA) ){
//beta = GetProperties()[RAYLEIGH_BETA];
//}
//else if( rCurrentProcessInfo.Has(RAYLEIGH_BETA) ){
//beta = rCurrentProcessInfo[RAYLEIGH_BETA];
//}

////4.-Compose the Damping Matrix:

////Rayleigh Damping Matrix: alpha*M + beta*K
//rDampingMatrix  = alpha * MassMatrix;
//rDampingMatrix += beta  * StiffnessMatrix;
////std::cout<<" rDampingMatrix "<<rDampingMatrix<<std::endl;

//KRATOS_CATCH( "" )
//}
//************************************************************************************
//****************MASS MATRIX*********************************************************

void UpdatedLagrangianUP::CalculateMassMatrix( MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    //I need to call the values of the shape function for the single element
    GeneralVariables Variables;
    this->InitializeGeneralVariables(Variables,rCurrentProcessInfo);

    //lumped
    unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    unsigned int MatSize = number_of_nodes * dimension + number_of_nodes;

    if ( rMassMatrix.size1() != MatSize )
        rMassMatrix.resize( MatSize, MatSize, false );

    rMassMatrix = ZeroMatrix( MatSize, MatSize );

    double TotalMass = 0;

    //TOTAL MASS OF ONE MP ELEMENT

    TotalMass = this->GetValue(MP_MASS);

    //LUMPED MATRIX

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        double temp = Variables.N[i] * TotalMass;
        unsigned int indexup = i * dimension + i;
        for ( unsigned int j = 0; j < dimension; j++ )
        {
            //unsigned int index = i * dimension + j;
            rMassMatrix( indexup+j, indexup+j ) = temp;
        }
    }
    //for ( unsigned int i = 0; i < number_of_nodes; i++ )
    //{
    //unsigned int indexupi = dimension * i + i;

    //for ( unsigned int j = 0; j < number_of_nodes; j++ )
    //{
    //unsigned int indexupj = dimension * j + j;

    //for ( unsigned int k = 0; k < dimension; k++ )
    //{
    //rMassMatrix( indexupi+k , indexupj+k ) += Variables.N[i] * Variables.N[j] * TotalMass;
    //}
    //}
    //}
    //CONSISTENT MATRIX
    //for ( unsigned int i = 0; i < number_of_nodes; i++ )
    //{
    //for ( unsigned int j = 0; j < number_of_nodes; j++ )
    //{

    //rMassMatrix( i*2, j*2 ) += Variables.N[i] * Variables.N[j] * TotalMass;
    //rMassMatrix( i * 2 + 1, j * 2 + 1 ) += Variables.N[i] * Variables.N[j] * TotalMass;

    //}
    //}


    //std::cout<<"rMassMatrix "<<rMassMatrix<<std::endl;

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************



//Matrix& UpdatedLagrangian::MPMJacobian( Matrix& rResult, array_1d<double,3>& rPoint)
//{

//KRATOS_TRY

////derivatives of shape functions
//Matrix shape_functions_gradients;
//shape_functions_gradients =this->MPMShapeFunctionsLocalGradients(
//shape_functions_gradients);
//const GeometryType& rGeom = GetGeometry();

//unsigned int number_nodes = rGeom.PointsNumber();
//unsigned int dimension = GetGeometry().WorkingSpaceDimension();

//if (dimension ==2)
//{
//rResult.resize( 2, 2);
//rResult = ZeroMatrix(2,2);


//for ( unsigned int i = 0; i < number_nodes; i++ )
//{


//rResult( 0, 0 ) += ( GetGeometry().GetPoint( i ).X() *  shape_functions_gradients( i, 0 ) );
//rResult( 0, 1 ) += ( GetGeometry().GetPoint( i ).X() *  shape_functions_gradients( i, 1 ) );
//rResult( 1, 0 ) += ( GetGeometry().GetPoint( i ).Y() *  shape_functions_gradients( i, 0 ) );
//rResult( 1, 1 ) += ( GetGeometry().GetPoint( i ).Y() *  shape_functions_gradients( i, 1 ) );


//}

//}
//else if(dimension ==3)
//{


//rResult.resize( 3,3);
//rResult = ZeroMatrix(3,3);

//for ( unsigned int i = 0; i < number_nodes; i++ )
//{
//rResult( 0, 0 ) += ( GetGeometry().GetPoint( i ).X() *  shape_functions_gradients( i, 0 ) );
//rResult( 0, 1 ) += ( GetGeometry().GetPoint( i ).X() *  shape_functions_gradients( i, 1 ) );
//rResult( 0, 2 ) += ( GetGeometry().GetPoint( i ).X() *  shape_functions_gradients( i, 2 ) );
//rResult( 1, 0 ) += ( GetGeometry().GetPoint( i ).Y() *  shape_functions_gradients( i, 0 ) );
//rResult( 1, 1 ) += ( GetGeometry().GetPoint( i ).Y() *  shape_functions_gradients( i, 1 ) );
//rResult( 1, 2 ) += ( GetGeometry().GetPoint( i ).Y() *  shape_functions_gradients( i, 2 ) );
//rResult( 2, 0 ) += ( GetGeometry().GetPoint( i ).Z() *  shape_functions_gradients( i, 0 ) );
//rResult( 2, 1 ) += ( GetGeometry().GetPoint( i ).Z() *  shape_functions_gradients( i, 1 ) );
//rResult( 2, 2 ) += ( GetGeometry().GetPoint( i ).Z() *  shape_functions_gradients( i, 2 ) );

//}

//}

//return rResult;

//KRATOS_CATCH( "" )
//}
//  ///**
//* Jacobian in given point and given a delta position. This method calculate jacobian
//* matrix in given point and a given delta position.
//*
//* @param rPoint point which jacobians has to
//* be calculated in it.
//*
//* @return Matrix of double which is jacobian matrix \f$ J \f$ in given point and a given delta position.
//*
//* @see DeterminantOfJacobian
//* @see InverseOfJacobian
//*/
//Matrix& UpdatedLagrangian::MPMJacobianDelta( Matrix& rResult, array_1d<double,3>& rPoint, Matrix & rDeltaPosition )
//{
//KRATOS_TRY

//Matrix shape_functions_gradients;

//shape_functions_gradients = this->MPMShapeFunctionsLocalGradients(
//shape_functions_gradients );



//unsigned int dimension = GetGeometry().WorkingSpaceDimension();
////Elements of jacobian matrix (e.g. J(1,1) = dX1/dXi1)

//if (dimension ==2)
//{

//rResult.resize( 2, 2);
//rResult = ZeroMatrix(2,2);

//for ( unsigned int i = 0; i < GetGeometry().size(); i++ )
//{
//rResult( 0, 0 ) += ( GetGeometry().GetPoint( i ).X() + rDeltaPosition(i,0)) * ( shape_functions_gradients( i, 0 ) );
//rResult( 0, 1 ) += ( GetGeometry().GetPoint( i ).X() + rDeltaPosition(i,0)) * ( shape_functions_gradients( i, 1 ) );
//rResult( 1, 0 ) += ( GetGeometry().GetPoint( i ).Y() + rDeltaPosition(i,1)) * ( shape_functions_gradients( i, 0 ) );
//rResult( 1, 1 ) += ( GetGeometry().GetPoint( i ).Y() + rDeltaPosition(i,1)) * ( shape_functions_gradients( i, 1 ) );
//}
//}
//else if(dimension ==3)
//{

//rResult.resize( 3,3);
//rResult = ZeroMatrix(3,3);
//for ( unsigned int i = 0; i < GetGeometry().size(); i++ )
//{
//rResult( 0, 0 ) += ( GetGeometry().GetPoint( i ).X() + rDeltaPosition(i,0)) * ( shape_functions_gradients( i, 0 ) );
//rResult( 0, 1 ) += ( GetGeometry().GetPoint( i ).X() + rDeltaPosition(i,0)) * ( shape_functions_gradients( i, 1 ) );
//rResult( 0, 2 ) += ( GetGeometry().GetPoint( i ).X() + rDeltaPosition(i,0)) * ( shape_functions_gradients( i, 2 ) );
//rResult( 1, 0 ) += ( GetGeometry().GetPoint( i ).Y() + rDeltaPosition(i,1)) * ( shape_functions_gradients( i, 0 ) );
//rResult( 1, 1 ) += ( GetGeometry().GetPoint( i ).Y() + rDeltaPosition(i,1)) * ( shape_functions_gradients( i, 1 ) );
//rResult( 1, 2 ) += ( GetGeometry().GetPoint( i ).Y() + rDeltaPosition(i,1)) * ( shape_functions_gradients( i, 2 ) );
//rResult( 2, 0 ) += ( GetGeometry().GetPoint( i ).Z() + rDeltaPosition(i,2)) * ( shape_functions_gradients( i, 0 ) );
//rResult( 2, 1 ) += ( GetGeometry().GetPoint( i ).Z() + rDeltaPosition(i,2)) * ( shape_functions_gradients( i, 1 ) );
//rResult( 2, 2 ) += ( GetGeometry().GetPoint( i ).Z() + rDeltaPosition(i,2)) * ( shape_functions_gradients( i, 2 ) );
//}
//}



//return rResult;

//KRATOS_CATCH( "" )
//}

////************************************************************************************
////************************************************************************************

//   ///**
//* Shape function values in given point. This method calculate the shape function
//* vector in given point.
//*
//* @param rPoint point which shape function values have to
//* be calculated in it.
//*
//* @return Vector of double which is shape function vector \f$ N \f$ in given point.
//*
//*/
//Vector& UpdatedLagrangian::MPMShapeFunctionPointValues( Vector& rResult, array_1d<double,3>& rPoint )
//{
//KRATOS_TRY

//unsigned int dimension = GetGeometry().WorkingSpaceDimension();
//Vector rPointLocal = ZeroVector(dimension);


//if (dimension == 2)
//{

//rResult.resize(3, false);

////1. I evaluate the local coordinates of a point
//rPointLocal[0] = ((GetGeometry()[2].Coordinates()[1] - GetGeometry()[0].Coordinates()[1])*(rPoint[0] - GetGeometry()[0].Coordinates()[0]) -
//(GetGeometry()[2].Coordinates()[0] - GetGeometry()[0].Coordinates()[0])*(rPoint[1] - GetGeometry()[0].Coordinates()[1]))/mDeterminantJ0;


//rPointLocal[1] = (-(GetGeometry()[1].Coordinates()[1] - GetGeometry()[0].Coordinates()[1])*(rPoint[0] - GetGeometry()[0].Coordinates()[0]) +
//(GetGeometry()[1].Coordinates()[0] - GetGeometry()[0].Coordinates()[0])*(rPoint[1] - GetGeometry()[0].Coordinates()[1]))/mDeterminantJ0;

////2. Shape functions
//rResult( 0 ) = 1 - rPointLocal[0] - rPointLocal[1] ;
//rResult( 1 ) = rPointLocal[0] ;
//rResult( 2 ) = rPointLocal[1];


//}

//else if (dimension == 3)
//{
//rResult.resize(4, false);

//double x10 = GetGeometry()[1].Coordinates()[0]-GetGeometry()[0].Coordinates()[0];
//double x20 = GetGeometry()[2].Coordinates()[0]-GetGeometry()[0].Coordinates()[0];
//double x30 = GetGeometry()[3].Coordinates()[0]-GetGeometry()[0].Coordinates()[0];
//double y10 = GetGeometry()[1].Coordinates()[1]-GetGeometry()[0].Coordinates()[1];
//double y20 = GetGeometry()[2].Coordinates()[1]-GetGeometry()[0].Coordinates()[1];
//double y30 = GetGeometry()[3].Coordinates()[1]-GetGeometry()[0].Coordinates()[1];
//double z10 = GetGeometry()[1].Coordinates()[2]-GetGeometry()[0].Coordinates()[2];
//double z20 = GetGeometry()[2].Coordinates()[2]-GetGeometry()[0].Coordinates()[2];
//double z30 = GetGeometry()[3].Coordinates()[2]-GetGeometry()[0].Coordinates()[2];

//rPointLocal[3] = (rPoint[0] - GetGeometry()[0].Coordinates()[0])*(y10*z20 - z10*y20) - (rPoint[1] - GetGeometry()[0].Coordinates()[1])*(x10*z20-x20*z10) + (rPoint[2] - GetGeometry()[0].Coordinates()[2])*(y20*x10 - y10*x20)/mDeterminantJ0;

//rPointLocal[2] = (rPoint[0] - GetGeometry()[0].Coordinates()[0])*(y30*z10-y10*z30) + (rPoint[1] - GetGeometry()[0].Coordinates()[1])*(x10*z30-x30*z10) + (rPoint[2] - GetGeometry()[0].Coordinates()[2])*(y10*x30 - y30*x10)/mDeterminantJ0;

//rPointLocal[1] = (rPoint[0] - GetGeometry()[0].Coordinates()[0])*(y20*z30-y30*z20) + (rPoint[1] - GetGeometry()[0].Coordinates()[1])*(x30*z20-x20*z30) + (rPoint[2] - GetGeometry()[0].Coordinates()[2])*(y30*x20 - x30*y20)/mDeterminantJ0;

//rPointLocal[0] = 1 - rPointLocal[1] - rPointLocal[2] -rPointLocal[3];


//rResult( 0 ) =  1.0-(rPointLocal[0]+rPointLocal[1]+rPointLocal[2]) ;
//rResult( 1 ) = rPointLocal[0] ;
//rResult( 2 ) = rPointLocal[1];
//rResult( 3 ) = rPointLocal[2];


//}

//return rResult;

//KRATOS_CATCH( "" )
//}



//Vector& UpdatedLagrangian::MPMLocalCoordinates(Vector& rResult, array_1d<double,3>& rPoint)
//{
//KRATOS_TRY

////Only local coordinated of a point in a tetrahedron is computed
//rResult.resize(4,false);

//double x10 = GetGeometry()[1].Coordinates()[0]-GetGeometry()[0].Coordinates()[0];
//double x20 = GetGeometry()[2].Coordinates()[0]-GetGeometry()[0].Coordinates()[0];
//double x30 = GetGeometry()[3].Coordinates()[0]-GetGeometry()[0].Coordinates()[0];
//double y10 = GetGeometry()[1].Coordinates()[1]-GetGeometry()[0].Coordinates()[1];
//double y20 = GetGeometry()[2].Coordinates()[1]-GetGeometry()[0].Coordinates()[1];
//double y30 = GetGeometry()[3].Coordinates()[1]-GetGeometry()[0].Coordinates()[1];
//double z10 = GetGeometry()[1].Coordinates()[2]-GetGeometry()[0].Coordinates()[2];
//double z20 = GetGeometry()[2].Coordinates()[2]-GetGeometry()[0].Coordinates()[2];
//double z30 = GetGeometry()[3].Coordinates()[2]-GetGeometry()[0].Coordinates()[2];

//rResult[3] = (rPoint[0] - GetGeometry()[0].Coordinates()[0])*(y10*z20 - z10*y20) - (rPoint[1] - GetGeometry()[0].Coordinates()[1])*(x10*z20-x20*z10) + (rPoint[2] - GetGeometry()[0].Coordinates()[2])*(y20*x10 - y10*x20)/mDeterminantJ0;

//rResult[2] = (rPoint[0] - GetGeometry()[0].Coordinates()[0])*(y30*z10-y10*z30) + (rPoint[1] - GetGeometry()[0].Coordinates()[1])*(x10*z30-x30*z10) + (rPoint[2] - GetGeometry()[0].Coordinates()[2])*(y10*x30 - y30*x10)/mDeterminantJ0;

//rResult[1] = (rPoint[0] - GetGeometry()[0].Coordinates()[0])*(y20*z30-y30*z20) + (rPoint[1] - GetGeometry()[0].Coordinates()[1])*(x30*z20-x20*z30) + (rPoint[2] - GetGeometry()[0].Coordinates()[2])*(y30*x20 - x30*y20)/mDeterminantJ0;

//rResult[0] = 1 - rResult[1] - rResult[2] -rResult[3];

//return rResult;

//KRATOS_CATCH( "" )
//}




//Matrix& UpdatedLagrangian::MPMShapeFunctionsLocalGradients( Matrix& rResult )
//{
//unsigned int dim = GetGeometry().WorkingSpaceDimension();
//if (dim == 2)
//{
//rResult = ZeroMatrix( 3, 2 );
//rResult( 0, 0 ) = -1.0;
//rResult( 0, 1 ) = -1.0;
//rResult( 1, 0 ) =  1.0;
//rResult( 1, 1 ) =  0.0;
//rResult( 2, 0 ) =  0.0;
//rResult( 2, 1 ) =  1.0;
//}
//else if(dim == 3)
//{
//rResult = ZeroMatrix( 4, 3 );
//rResult(0,0) = -1.0;
//rResult(0,1) = -1.0;
//rResult(0,2) = -1.0;
//rResult(1,0) =  1.0;
//rResult(1,1) =  0.0;
//rResult(1,2) =  0.0;
//rResult(2,0) =  0.0;
//rResult(2,1) =  1.0;
//rResult(2,2) =  0.0;
//rResult(3,0) =  0.0;
//rResult(3,1) =  0.0;
//rResult(3,2) =  1.0;
//}

//return rResult;
//}
//************************************************************************************
//************************************************************************************

//void UpdatedLagrangian::CalculateOnIntegrationPoints( const Variable<double>& rVariable, double& rOutput, ProcessInfo& rCurrentProcessInfo )
//{

//rOutput = mConstitutiveLawVector->GetValue( rVariable, rOutput );
//}

////************************************************************************************
////************************************************************************************

//void UpdatedLagrangian::CalculateOnIntegrationPoints( const Variable<Vector>& rVariable, Vector& rOutput, ProcessInfo& rCurrentProcessInfo )
//{

//KRATOS_TRY





//if ( rVariable == CAUCHY_STRESS_VECTOR || rVariable == PK2_STRESS_VECTOR )
//{
////create and initialize element variables:
//GeneralVariables Variables;
//this->InitializeGeneralVariables(Variables,rCurrentProcessInfo);

////create constitutive law parameters:
//ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

////set constitutive law flags:
//Flags &ConstitutiveLawOptions=Values.GetOptions();

//ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);
//ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS);


////compute element kinematics B, F, DN_DX ...
//this->CalculateKinematics(Variables, rCurrentProcessInfo);

////to take in account previous step writing
//if( mFinalizedStep ){
//this->GetHistoricalVariables(Variables);
//}

////set general variables to constitutivelaw parameters
//this->SetGeneralVariables(Variables,Values);

////call the constitutive law to update material variables
//if( rVariable == CAUCHY_STRESS_VECTOR)
//mConstitutiveLawVector->CalculateMaterialResponseCauchy(Values);
//else
//mConstitutiveLawVector->CalculateMaterialResponsePK2(Values);


//rOutput = Variables.StressVector;


//}


//else if( rVariable == GREEN_LAGRANGE_STRAIN_VECTOR  || rVariable == ALMANSI_STRAIN_VECTOR )
//{
////create and initialize element variables:
//GeneralVariables Variables;
//this->InitializeGeneralVariables(Variables,rCurrentProcessInfo);


////compute element kinematics B, F, DN_DX ...
//this->CalculateKinematics(Variables, rCurrentProcessInfo);

////to take in account previous step writing
//if( mFinalizedStep ){
//this->GetHistoricalVariables(Variables);
//Variables.FT = prod(Variables.F,Variables.F0);
//}

////Compute Green-Lagrange Strain
//if( rVariable == GREEN_LAGRANGE_STRAIN_VECTOR )
//this->CalculateGreenLagrangeStrain( Variables.FT, Variables.StrainVector );
//else
//this->CalculateAlmansiStrain( Variables.FT, Variables.StrainVector );

//if ( rOutput.size() != Variables.StrainVector.size() )
//rOutput.resize( Variables.StrainVector.size(), false );

//rOutput = Variables.StrainVector;



//}
//else
//{

//rOutput = mConstitutiveLawVector->GetValue( rVariable , rOutput );

//}

//KRATOS_CATCH( "" )
//}

////************************************************************************************
////************************************************************************************

//void UpdatedLagrangian::CalculateOnIntegrationPoints( const Variable<Matrix >& rVariable, Matrix& rOutput, ProcessInfo& rCurrentProcessInfo )
//{
//KRATOS_TRY



//const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();



//if ( rVariable == CAUCHY_STRESS_TENSOR || rVariable == PK2_STRESS_TENSOR )
//{
////std::vector<Vector> StressVector;
//Vector StressVector;
//if( rVariable == CAUCHY_STRESS_TENSOR )
//this->CalculateOnIntegrationPoints( CAUCHY_STRESS_VECTOR, StressVector, rCurrentProcessInfo );
//else
//this->CalculateOnIntegrationPoints( PK2_STRESS_VECTOR, StressVector, rCurrentProcessInfo );


//if ( rOutput.size2() != dimension )
//rOutput.resize( dimension, dimension, false );

//rOutput = MathUtils<double>::StressVectorToTensor(StressVector);


//}
//else if ( rVariable == GREEN_LAGRANGE_STRAIN_TENSOR  || rVariable == ALMANSI_STRAIN_TENSOR)
//{


//Vector StrainVector;
//if( rVariable == GREEN_LAGRANGE_STRAIN_TENSOR )
//CalculateOnIntegrationPoints( GREEN_LAGRANGE_STRAIN_VECTOR, StrainVector, rCurrentProcessInfo );
//else
//CalculateOnIntegrationPoints( ALMANSI_STRAIN_VECTOR, StrainVector, rCurrentProcessInfo );



//if ( rOutput.size2() != dimension )
//rOutput.resize( dimension, dimension, false );

//rOutput = MathUtils<double>::StrainVectorToTensor(StrainVector);

//}
//else if ( rVariable == CONSTITUTIVE_MATRIX )
//{
////create and initialize element variables:
//GeneralVariables Variables;
//this->InitializeGeneralVariables(Variables,rCurrentProcessInfo);

////create constitutive law parameters:
//ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

////set constitutive law flags:
//Flags &ConstitutiveLawOptions=Values.GetOptions();

//ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
////ConstitutiveLawOptions.Set(ConstitutiveLaw::LAST_KNOWN_CONFIGURATION); //contact domain formulation UL

////compute element kinematics B, F, DN_DX ...
//this->CalculateKinematics(Variables, rCurrentProcessInfo);

////set general variables to constitutivelaw parameters
//this->SetGeneralVariables(Variables,Values);

////call the constitutive law to update material variables
////mConstitutiveLawVector[PointNumber]->CalculateMaterialResponsePK2(Values); //contact domain formulation UL
//mConstitutiveLawVector->CalculateMaterialResponseCauchy(Values); //contact domain formulation SL

//if( rOutput.size2() != Variables.ConstitutiveMatrix.size2() )
//rOutput.resize( Variables.ConstitutiveMatrix.size1() , Variables.ConstitutiveMatrix.size2() , false );

//rOutput = Variables.ConstitutiveMatrix;

////}


//}
//else if ( rVariable == DEFORMATION_GRADIENT )  // VARIABLE SET FOR TRANSFER PURPOUSES
//{
////create and initialize element variables:
//GeneralVariables Variables;
//this->InitializeGeneralVariables(Variables,rCurrentProcessInfo);


////compute element kinematics B, F, DN_DX ...
//this->CalculateKinematics(Variables, rCurrentProcessInfo);

//if( rOutput.size2() != Variables.F.size2() )
//rOutput.resize( Variables.F.size1() , Variables.F.size2() , false );

//rOutput = Variables.F;


//}
//else
//{

//rOutput = mConstitutiveLawVector->GetValue( rVariable , rOutput );

//}

//KRATOS_CATCH( "" )
//}
//*********************************SET DOUBLE VALUE***********************************
//************************************************************************************

//void UpdatedLagrangian::SetValueOnIntegrationPoints( const Variable<double>& rVariable,
//double& rValues,
//ProcessInfo& rCurrentProcessInfo )
//{
//if (rVariable == DETERMINANT_F){


//mDeterminantF0 = rValues;
//mConstitutiveLawVector->SetValue(rVariable, rValues, rCurrentProcessInfo);


//}
//else{

//UpdatedLagrangian::SetValueOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
//}
//}
////********************************SET VECTOR VALUE************************************
////************************************************************************************

//void UpdatedLagrangian::SetValueOnIntegrationPoints( const Variable<Vector>& rVariable, Vector& rValues, ProcessInfo& rCurrentProcessInfo )
//{

//mConstitutiveLawVector->SetValue( rVariable, rValues, rCurrentProcessInfo );


//}


////*******************************SET MATRIX VALUE*************************************
////************************************************************************************

//void UpdatedLagrangian::SetValueOnIntegrationPoints( const Variable<Matrix>& rVariable, Matrix& rValues, ProcessInfo& rCurrentProcessInfo )
//{

//mConstitutiveLawVector->SetValue( rVariable,
//rValues, rCurrentProcessInfo );


//}
////********************************SET CONSTITUTIVE VALUE******************************
////************************************************************************************

//void UpdatedLagrangian::SetValueOnIntegrationPoints( const Variable<ConstitutiveLaw::Pointer>& rVariable,
//ConstitutiveLaw::Pointer& rValues,
//ProcessInfo& rCurrentProcessInfo )
//{
//if(rVariable == CONSTITUTIVE_LAW)
//{

//mConstitutiveLawVector = rValues->Clone();

//}

//if(rVariable == CONSTITUTIVE_LAW_POINTER)
//{

//mConstitutiveLawVector = rValues;

//}

//}
////***************************GET DOUBLE VALUE*****************************************
////************************************************************************************

//void UpdatedLagrangian::GetValueOnIntegrationPoints( const Variable<double>& rVariable,
//double& rValues,
//ProcessInfo& rCurrentProcessInfo )
//{
//if (rVariable == DETERMINANT_F){


//rValues = mDeterminantF0;

//}
//else{

//rValues = mConstitutiveLawVector->GetValue( rVariable, rValues );
//}
//}

////**************************GET VECTOR VALUE******************************************
////************************************************************************************

//void UpdatedLagrangian::GetValueOnIntegrationPoints( const Variable<Vector>& rVariable, Vector& rValues, ProcessInfo& rCurrentProcessInfo )
//{

//if ( rVariable == PK2_STRESS_VECTOR ||  rVariable == CAUCHY_STRESS_VECTOR )
//{

//CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );

//}
//else if ( rVariable == GREEN_LAGRANGE_STRAIN_TENSOR ||  rVariable == ALMANSI_STRAIN_VECTOR )
//{

//CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );

//}
//else
//{


//rValues = mConstitutiveLawVector->GetValue( rVariable, rValues );


//}


//}

////************************************************************************************
////************************************************************************************

//void UpdatedLagrangian::GetValueOnIntegrationPoints( const Variable<Matrix>& rVariable,
//Matrix& rValues, ProcessInfo& rCurrentProcessInfo )
//{


//if ( rVariable == GREEN_LAGRANGE_STRAIN_TENSOR ||  rVariable == ALMANSI_STRAIN_TENSOR )
//{
//CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
//}

//if ( rVariable == PK2_STRESS_TENSOR ||  rVariable == CAUCHY_STRESS_TENSOR )
//{
//CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
//}
//else
//{

//rValues = mConstitutiveLawVector->GetValue( rVariable, rValues );

//}

//}
////********************************GET CONSTITUTIVE VALUE******************************
////************************************************************************************

//void UpdatedLagrangian::GetValueOnIntegrationPoints( const Variable<ConstitutiveLaw::Pointer>& rVariable,
//ConstitutiveLaw::Pointer& rValues,
//ProcessInfo& rCurrentProcessInfo )
//{

//if(rVariable == CONSTITUTIVE_LAW || rVariable == CONSTITUTIVE_LAW_POINTER)
//{

//rValues = mConstitutiveLawVector;

//}

//}
//************************************************************************************
//************************************************************************************

void UpdatedLagrangianUP::GetValuesVector( Vector& values, int Step )
{
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
    unsigned int       element_size    = number_of_nodes * dimension + number_of_nodes;

    if ( values.size() != element_size ) values.resize( element_size, false );


    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        unsigned int index = i * dimension + i;
        values[index]     = GetGeometry()[i].FastGetSolutionStepValue( DISPLACEMENT_X, Step );
        values[index + 1] = GetGeometry()[i].FastGetSolutionStepValue( DISPLACEMENT_Y, Step );

        if ( dimension == 3 )
        {
            values[index + 2] = GetGeometry()[i].FastGetSolutionStepValue( DISPLACEMENT_Z, Step );
            values[index + 3] = GetGeometry()[i].FastGetSolutionStepValue( PRESSURE, Step );
        }
        else
        {
            values[index + 2] = GetGeometry()[i].FastGetSolutionStepValue( PRESSURE, Step );
        }

    }
}


//************************************************************************************
//************************************************************************************

void UpdatedLagrangianUP::GetFirstDerivativesVector( Vector& values, int Step )
{
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
    unsigned int       element_size    = number_of_nodes * dimension + number_of_nodes;

    if ( values.size() != element_size ) values.resize( element_size, false );

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        unsigned int index = i * dimension + i;
        values[index]     = GetGeometry()[i].FastGetSolutionStepValue( VELOCITY_X, Step );
        values[index + 1] = GetGeometry()[i].FastGetSolutionStepValue( VELOCITY_Y, Step );
        if ( dimension == 3 )
        {
            values[index + 2] = GetGeometry()[i].FastGetSolutionStepValue( VELOCITY_Z, Step );
            values[index + 3] = 0;
        }
        else
        {
            values[index + 2] = 0;
        }
    }
}

//************************************************************************************
//************************************************************************************

void UpdatedLagrangianUP::GetSecondDerivativesVector( Vector& values, int Step )
{
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
    unsigned int       element_size    = number_of_nodes * dimension + number_of_nodes;

    if ( values.size() != element_size ) values.resize( element_size, false );


    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        unsigned int index = i * dimension + i;
        values[index]     = GetGeometry()[i].FastGetSolutionStepValue( ACCELERATION_X, Step );
        values[index + 1] = GetGeometry()[i].FastGetSolutionStepValue( ACCELERATION_Y, Step );

        if ( dimension == 3 )
        {
            values[index + 2] = GetGeometry()[i].FastGetSolutionStepValue( ACCELERATION_Z, Step );
            values[index + 3] = 0;
        }
        else
        {
            values[index + 2] = 0;
        }
    }
}
//************************************************************************************
//************************************************************************************
void UpdatedLagrangianUP::GetHistoricalVariables( GeneralVariables& rVariables )
{
    //Deformation Gradient F ( set to identity )
    unsigned int size =  rVariables.F.size1();
    rVariables.detF  = 1;
    rVariables.F     = IdentityMatrix(size);

    rVariables.detF0 = mDeterminantF0;
    rVariables.F0    = mDeformationGradientF0;

}

void UpdatedLagrangianUP::FinalizeStepVariables( GeneralVariables & rVariables, const ProcessInfo& rCurrentProcessInfo)
{
    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    double voigtsize = 3;
    if ( dimension == 3)
        voigtsize = 6;

    UpdatedLagrangian::FinalizeStepVariables( rVariables, rCurrentProcessInfo);

    //evaluation of the pressure on the material point
    double NodalMeanStress = 0.0;
    for (unsigned int i = 0; i < number_of_nodes; i++)
        NodalMeanStress += GetGeometry()[i].FastGetSolutionStepValue( PRESSURE ) * rVariables.N[i];

    //evaluation of the mean stress on the material point
    double MeanStress = 0.0;
    for (unsigned int i = 0; i < dimension; i++)
        MeanStress += rVariables.StressVector(i);
    MeanStress /= dimension;

    Vector StressVector = ZeroVector(voigtsize);
    StressVector = rVariables.StressVector;
    for (unsigned int i = 0; i < dimension; i++)
        StressVector(i) += (NodalMeanStress - MeanStress);

    this->SetValue(MP_CAUCHY_STRESS_VECTOR, StressVector);

}
//************************************************************************************
//************************************************************************************

//void UpdatedLagrangian::Calculate( const Variable<double>& rVariable, double& Output, const ProcessInfo& rCurrentProcessInfo )
//{

//double lamda = 1.00; // parametro que depende del tipo de problema y del elemento pag 308 libro dinamica de Barbat
//double c1 = 0.00; //sqrt(GetProperties()[YOUNG_MODULUS]/GetProperties()[DENSITY]); velocidad del sonido en el medio
//double c2 = 0.00; // norma de la velocidad actual dentro del elemento
//double c = 0.00;
//double wmax = 0.00;
//Vector Values( GetGeometry().IntegrationPoints( mThisIntegrationMethod ).size() );
//Vector Velocities;

//GetFirstDerivativesVector( Velocities, 0 );

//if ( rVariable == DELTA_TIME )
//{
//for ( unsigned int PointNumber = 0;
//PointNumber < GetGeometry().IntegrationPoints( mThisIntegrationMethod ).size();
//PointNumber++ )
//{
//mConstitutiveLawVector[PointNumber]-> GetValue( DELTA_TIME, c1 );
//Values[PointNumber] = c1;
//}
//}

//c1 = ( *std::max_element( Values.begin(), Values.end() ) );

//c2 = norm_2( Velocities );

//c = ( c1 > c2 ) ? c1 : c2;


//double le = GetGeometry().Length();
////KRATOS_WATCH(le)

///// maxima frecuencia de un elemento
//wmax = ( lamda * c ) / le;
//Output = 2.0 / wmax;
////KRATOS_WATCH(Output)

//}

//************************************************************************************
//************************************************************************************

//void UpdatedLagrangian::Comprobate_State_Vector( Vector& Result )
//{
//for ( unsigned int i = 0.00; i < Result.size(); i++ )
//{
//if ( fabs( Result( i ) ) < 1E-9 )
//{
//Result( i ) = 0.00;
//}
//}
//}

//*************************DECIMAL CORRECTION OF STRAINS******************************
//************************************************************************************

//void UpdatedLagrangian::DecimalCorrection(Vector& rVector)
//{
//KRATOS_TRY

//for ( unsigned int i = 0; i < rVector.size(); i++ )
//{
//if( rVector[i]*rVector[i]<1e-24 )
//{
//rVector[i]=0;
//}

//}

//KRATOS_CATCH( "" )
//}

//************************************************************************************
//************************************************************************************
/**
 * This function provides the place to perform checks on the completeness of the input.
 * It is designed to be called only once (or anyway, not often) typically at the beginning
 * of the calculations, so to verify that nothing is missing from the input
 * or that no common error is found.
 * @param rCurrentProcessInfo
 */
int  UpdatedLagrangianUP::Check( const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY
    int correct = 0;

    correct = UpdatedLagrangian::Check(rCurrentProcessInfo);


    //verify compatibility with the constitutive law
    ConstitutiveLaw::Features LawFeatures;
    this->GetProperties().GetValue( CONSTITUTIVE_LAW )->GetLawFeatures(LawFeatures);

    if(LawFeatures.mOptions.IsNot(ConstitutiveLaw::U_P_LAW))
        KRATOS_THROW_ERROR( std::logic_error, "constitutive law is not compatible with the U-P element type ", " Large Displacements U_P" )

        //verify that the variables are correctly initialized

        if ( PRESSURE.Key() == 0 )
            KRATOS_THROW_ERROR( std::invalid_argument, "PRESSURE has Key zero! (check if the application is correctly registered", "" )

            return correct;



    KRATOS_CATCH( "" );
}




void UpdatedLagrangianUP::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Element )
    //int IntMethod = int(mThisIntegrationMethod);
    //rSerializer.save("IntegrationMethod",IntMethod);
    rSerializer.save("ConstitutiveLawVector",mConstitutiveLawVector);
    rSerializer.save("DeformationGradientF0",mDeformationGradientF0);
    rSerializer.save("DeterminantF0",mDeterminantF0);


}

void UpdatedLagrangianUP::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Element )
    rSerializer.load("ConstitutiveLawVector",mConstitutiveLawVector);
    rSerializer.load("DeformationGradientF0",mDeformationGradientF0);
    rSerializer.load("DeterminantF0",mDeterminantF0);
    //int IntMethod;
    //rSerializer.load("IntegrationMethod",IntMethod);
    //mThisIntegrationMethod = IntegrationMethod(IntMethod);


}





} // Namespace Kratos

