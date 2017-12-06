/*
==============================================================================
KratosStructuralApplication
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi, Janosch Stascheit, Felix Nagel
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
janosch.stascheit@rub.de
nagel@sd.rub.de
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain
- Ruhr-University Bochum, Institute for Structural Mechanics, Germany


Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
 */

//
//   Project Name:        KratosParticleMechanicsApplication $
//   Last modified by:    $Author:                 ilaria $
//   Date:                $Date:                July 2015 $
//   Revision:            $Revision:                  0.0 $
//
//


// System includes

// External includes


// Project includes
#include "includes/define.h"
#include "custom_elements/updated_lagrangian_quadrilateral.hpp"
#include "utilities/math_utils.h"
#include "includes/constitutive_law.h"
#include "particle_mechanics_application.h"


#include <omp.h>
#include <sstream>

namespace Kratos
{

/**
 * Flags related to the element computation
 */
KRATOS_CREATE_LOCAL_FLAG( UpdatedLagrangianQuadrilateral, COMPUTE_RHS_VECTOR,                 0 );
KRATOS_CREATE_LOCAL_FLAG( UpdatedLagrangianQuadrilateral, COMPUTE_LHS_MATRIX,                 1 );
KRATOS_CREATE_LOCAL_FLAG( UpdatedLagrangianQuadrilateral, COMPUTE_RHS_VECTOR_WITH_COMPONENTS, 2 );
KRATOS_CREATE_LOCAL_FLAG( UpdatedLagrangianQuadrilateral, COMPUTE_LHS_MATRIX_WITH_COMPONENTS, 3 );

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

UpdatedLagrangianQuadrilateral::UpdatedLagrangianQuadrilateral( )
    : Element( )
{
    //DO NOT CALL IT: only needed for Register and Serialization!!!
}
//******************************CONSTRUCTOR*******************************************
//************************************************************************************
UpdatedLagrangianQuadrilateral::UpdatedLagrangianQuadrilateral( IndexType NewId, GeometryType::Pointer pGeometry )
    : Element( NewId, pGeometry )
{
    //DO NOT ADD DOFS HERE!!!
}

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

UpdatedLagrangianQuadrilateral::UpdatedLagrangianQuadrilateral( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
    : Element( NewId, pGeometry, pProperties )
{
    mFinalizedStep = true;


}
//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

UpdatedLagrangianQuadrilateral::UpdatedLagrangianQuadrilateral( UpdatedLagrangianQuadrilateral const& rOther)
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

UpdatedLagrangianQuadrilateral&  UpdatedLagrangianQuadrilateral::operator=(UpdatedLagrangianQuadrilateral const& rOther)
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

Element::Pointer UpdatedLagrangianQuadrilateral::Create( IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties ) const
{
    return Element::Pointer( new UpdatedLagrangianQuadrilateral( NewId, GetGeometry().Create( ThisNodes ), pProperties ) );
}
//************************************CLONE*******************************************
//************************************************************************************

Element::Pointer UpdatedLagrangianQuadrilateral::Clone( IndexType NewId, NodesArrayType const& rThisNodes ) const
{

    UpdatedLagrangianQuadrilateral NewElement (NewId, GetGeometry().Create( rThisNodes ), pGetProperties() );

    //-----------//





    NewElement.mConstitutiveLawVector = mConstitutiveLawVector->Clone();



    //-----------//



    NewElement.mDeformationGradientF0 = mDeformationGradientF0;

    NewElement.mInverseJ0 = mInverseJ0;
    NewElement.mInverseJ = mInverseJ;


    NewElement.mDeterminantF0 = mDeterminantF0;
    NewElement.mDeterminantJ0 = mDeterminantJ0;


    return Element::Pointer( new UpdatedLagrangianQuadrilateral(NewElement) );
}
//*******************************DESTRUCTOR*******************************************
//************************************************************************************
UpdatedLagrangianQuadrilateral::~UpdatedLagrangianQuadrilateral()
{
}


//************************************************************************************
//************************************************************************************

void UpdatedLagrangianQuadrilateral::Initialize()
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


    double MP_KineticEnergy = 0.0;
    double MP_StrainEnergy = 0.0;

    for(unsigned int k = 0; k<3; k++)
    {
        MP_KineticEnergy += 0.5 * this->GetValue(MP_MASS) * this->GetValue(MP_VELOCITY)[k] * this->GetValue(MP_VELOCITY)[k] ;
    }
    for(unsigned int j = 0; j < this->GetValue(MP_CAUCHY_STRESS_VECTOR).size(); j++)
    {
        MP_StrainEnergy +=  0.5 * this->GetValue(MP_VOLUME) * this->GetValue(MP_CAUCHY_STRESS_VECTOR)[j] * this->GetValue(MP_ALMANSI_STRAIN_VECTOR)[j];
    }

    this->GetValue(MP_KINETIC_ENERGY) = MP_KineticEnergy;
    this->GetValue(MP_STRAIN_ENERGY) = MP_StrainEnergy;
    this->GetValue(MP_DENSITY) = GetProperties()[DENSITY];


    //std::cout<<"The element is initialized"<<std::endl;


    KRATOS_CATCH( "" )
}
//************************************************************************************
//************************************************************************************

void UpdatedLagrangianQuadrilateral::InitializeGeneralVariables (GeneralVariables& rVariables, const ProcessInfo& rCurrentProcessInfo)
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

    rVariables.B.resize( voigtsize, number_of_nodes * dimension );

    rVariables.F.resize( dimension, dimension );

    rVariables.F0.resize( dimension, dimension );

    rVariables.FT.resize( dimension, dimension );

    rVariables.ConstitutiveMatrix.resize( voigtsize, voigtsize );

    rVariables.Normal.resize( voigtsize, voigtsize);

    rVariables.StrainVector.resize( voigtsize );

    rVariables.StressVector.resize( voigtsize );

    rVariables.IsoStressVector.resize( voigtsize );



    rVariables.DN_DX.resize( number_of_nodes, dimension );
    rVariables.DN_De.resize( number_of_nodes, dimension );

    array_1d<double,3>& xg = this->GetValue(GAUSS_COORD);

    rVariables.N = this->MPMShapeFunctionPointValues(rVariables.N, xg);


    //reading shape functions local gradients
    rVariables.DN_De = this->MPMShapeFunctionsLocalGradients( rVariables.DN_De, xg);


    //**********************************************************************************************************************

    //CurrentDisp is the variable unknown. It represents the nodal delta displacement. When it is predicted is equal to zero.

    rVariables.CurrentDisp = CalculateCurrentDisp(rVariables.CurrentDisp, rCurrentProcessInfo);


    //calculating the current jacobian from cartesian coordinates to parent coordinates for the MP element [dx_n+1/d£]
    rVariables.j = this->MPMJacobianDelta( rVariables.j, xg, rVariables.CurrentDisp);


    //calculating the reference jacobian from cartesian coordinates to parent coordinates for the MP element [dx_n/d£]
    rVariables.J = this->MPMJacobian( rVariables.J, xg);

    //std::cout<<"The general variables are initialized"<<std::endl;
    //*************************************************************************************************************************

}
//************************************************************************************
//************************************************************************************

void UpdatedLagrangianQuadrilateral::SetGeneralVariables(GeneralVariables& rVariables,
        ConstitutiveLaw::Parameters& rValues)
{
    //Variables.detF is the determinant of the incremental total deformation gradient
    rVariables.detF  = MathUtils<double>::Det(rVariables.F);

    if(rVariables.detF<0)
    {

        std::cout<<" Element: "<<this->Id()<<std::endl;
        std::cout<<" Element position "<<this->GetValue(GAUSS_COORD)<<std::endl;
        unsigned int number_of_nodes = GetGeometry().PointsNumber();

        for ( unsigned int i = 0; i < number_of_nodes; i++ )
        {
            array_1d<double, 3> &CurrentPosition  = GetGeometry()[i].Coordinates();
            array_1d<double, 3> & CurrentDisplacement  = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT);
            array_1d<double, 3> & PreviousDisplacement = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT,1);
            //array_1d<double, 3> PreviousPosition  = CurrentPosition - (CurrentDisplacement-PreviousDisplacement);
            std::cout<<" NODE ["<<GetGeometry()[i].Id()<<"]: (Current position: "<<CurrentPosition<<") "<<std::endl;
            std::cout<<" ---Current Disp: "<<CurrentDisplacement<<" (Previour Disp: "<<PreviousDisplacement<<")"<<std::endl;
        }

        for ( unsigned int i = 0; i < number_of_nodes; i++ )
        {
            if( GetGeometry()[i].SolutionStepsDataHas(CONTACT_FORCE) )
            {
                array_1d<double, 3 > & PreContactForce = GetGeometry()[i].FastGetSolutionStepValue(CONTACT_FORCE,1);
                array_1d<double, 3 > & ContactForce = GetGeometry()[i].FastGetSolutionStepValue(CONTACT_FORCE);
                std::cout<<" ---Contact_Force: (Pre:"<<PreContactForce<<", Cur:"<<ContactForce<<") "<<std::endl;
            }
            else
            {
                std::cout<<" ---Contact_Force: NULL "<<std::endl;
            }
        }

        KRATOS_THROW_ERROR( std::invalid_argument," MPM UPDATED LAGRANGIAN DISPLACEMENT ELEMENT INVERTED: |F|<0  detF = ", rVariables.detF )
    }

    rVariables.detFT = rVariables.detF * rVariables.detF0;
    rVariables.FT    = prod( rVariables.F, rVariables.F0 );
    //if (this->Id() == 541 || this->Id() == 534 || this->Id() == 538)
    //{
    //std::cout<<" Total Deformation Gradient "<< rVariables.FT  <<std::endl;
    //std::cout<<" Total Deformation Gradient determinant "<< rVariables.detFT  <<std::endl;
    //}

    rValues.SetDeterminantF(rVariables.detFT);
    rValues.SetDeformationGradientF(rVariables.FT);
    rValues.SetStrainVector(rVariables.StrainVector);
    rValues.SetStressVector(rVariables.StressVector);
    rValues.SetConstitutiveMatrix(rVariables.ConstitutiveMatrix);
    rValues.SetShapeFunctionsDerivatives(rVariables.DN_DX);
    rValues.SetShapeFunctionsValues(rVariables.N);


    //std::cout<<"The general variables are set"<<std::endl;

}
//************************************************************************************
//*****************check size of LHS and RHS matrices*********************************

void UpdatedLagrangianQuadrilateral::InitializeSystemMatrices(MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        Flags& rCalculationFlags)

{

    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

    //resizing as needed the LHS

    unsigned int MatSize = number_of_nodes * dimension;   //number of degrees of freedom

    if ( rCalculationFlags.Is(UpdatedLagrangianQuadrilateral::COMPUTE_LHS_MATRIX) ) //calculation of the matrix is required
    {
        if ( rLeftHandSideMatrix.size1() != MatSize )
            rLeftHandSideMatrix.resize( MatSize, MatSize, false );

        noalias( rLeftHandSideMatrix ) = ZeroMatrix( MatSize, MatSize ); //resetting LHS
    }


    //resizing as needed the RHS
    if ( rCalculationFlags.Is(UpdatedLagrangianQuadrilateral::COMPUTE_RHS_VECTOR) ) //calculation of the matrix is required
    {
        if ( rRightHandSideVector.size() != MatSize )
            rRightHandSideVector.resize( MatSize, false );

        rRightHandSideVector = ZeroVector( MatSize ); //resetting RHS

    }
    //std::cout<<"The system matrices are initialized"<<std::endl;
}

//************************************************************************************
//************************************************************************************

void UpdatedLagrangianQuadrilateral::CalculateElementalSystem( LocalSystemComponents& rLocalSystem,
        ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    //const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    //create and initialize element variables:
    GeneralVariables Variables;

    this->InitializeGeneralVariables(Variables,rCurrentProcessInfo);


    //create constitutive law parameters:
    ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);


    //set constitutive law flags:
    Flags &ConstitutiveLawOptions=Values.GetOptions();

    //std::cout<<"in CalculateElementalSystem 5"<<std::endl;
    ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);

    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS);

    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);


    //auxiliary terms
    Vector VolumeForce;


    //compute element kinematics B, F, DN_DX ...
    //if (this->Id() == 541 || this->Id() == 534 || this->Id() == 538)
    //{
    //std::cout<<" in calculate elemental system "<<std::endl;
    //}
    this->CalculateKinematics(Variables,rCurrentProcessInfo);

    //set general variables to constitutivelaw parameters
    this->SetGeneralVariables(Variables,Values);

    mConstitutiveLawVector->CalculateMaterialResponse(Values, Variables.StressMeasure);

    //double TraceStress = 0;
    //Matrix I=identity_matrix<double>( dimension );
    //Matrix StressTensor = MathUtils<double>::StressVectorToTensor( Variables.StressVector );
    //for( unsigned int i=0; i<StressTensor.size1(); i++)
    //{
    //TraceStress += StressTensor( i , i );
    //}
    //double Pressure = TraceStress/StressTensor.size1();

    //Matrix IsoStressTensor = MathUtils<double>::StressVectorToTensor( Variables.IsoStressVector );
    //IsoStressTensor = StressTensor - Pressure * I;

    //double IsoStressNorm = 0;
    ////IsoStressNorm = sqrt((IsoStressTensor(0,0)*IsoStressTensor(0,0))+(IsoStressTensor(1,1)*IsoStressTensor(1,1))+(IsoStressTensor(2,2)*IsoStressTensor(2,2))+
    ////(IsoStressTensor(0,1)*IsoStressTensor(0,1))+(IsoStressTensor(0,2)*IsoStressTensor(0,2))+(IsoStressTensor(1,2)*IsoStressTensor(1,2))+
    ////(IsoStressTensor(1,0)*IsoStressTensor(1,0))+(IsoStressTensor(2,0)*IsoStressTensor(2,0))+(IsoStressTensor(2,1)*IsoStressTensor(2,1)));

    //IsoStressNorm = sqrt((IsoStressTensor(0,0)*IsoStressTensor(0,0))+(IsoStressTensor(1,1)*IsoStressTensor(1,1))+
    //(IsoStressTensor(0,1)*IsoStressTensor(0,1))+(IsoStressTensor(1,0)*IsoStressTensor(1,0)));
    //Variables.Normal =  IsoStressTensor / IsoStressNorm;
    //if (IsoStressTensor(0,0) == 0 || IsoStressTensor(1,1) == 0 || IsoStressTensor(0,1) == 0 || IsoStressTensor(1,0))
    //{
    //Variables.Normal = ZeroMatrix(dimension);
    //}


    //std::cout<<"IsoStressTensor "<<IsoStressTensor<<std::endl;
    //std::cout<<"IsoStressNorm "<<IsoStressNorm<<std::endl;
    //std::cout<<"Variables.Normal "<<Variables.Normal<<std::endl;

    //this->SetValue(MP_CAUCHY_STRESS_VECTOR, Variables.StressVector);
    //std::cout<<"AAAAAAAAAAAAAAAAAAAAAAAAAA"<<std::endl;
    //std::cout<<"Variables.StressVector in the element "<<Variables.StressVector<<std::endl;

    //this->SetValue(MP_ALMANSI_STRAIN_VECTOR, Variables.StrainVector);
    //double EquivalentPlasticStrain = mConstitutiveLawVector->GetValue( PLASTIC_STRAIN, EquivalentPlasticStrain );
    //this->SetValue(MP_EQUIVALENT_PLASTIC_STRAIN, EquivalentPlasticStrain);
    //at the first iteration I recover the previous state of stress and strain
    if(rCurrentProcessInfo[NL_ITERATION_NUMBER] == 1)
    {
        this->SetValue(PREVIOUS_MP_CAUCHY_STRESS_VECTOR, Variables.StressVector);
        this->SetValue(PREVIOUS_MP_ALMANSI_STRAIN_VECTOR, Variables.StrainVector);
    }
    //the MP density is updated
    double MP_Density = (GetProperties()[DENSITY]) / Variables.detFT;
    //if(this->Id() == 1786 || this->Id() == 1836)
    //{
    //std::cout<<"density "<<this->Id() << GetProperties()[DENSITY]<<std::endl;
    //}
        
    //the integration weight is evaluated
    double MP_Volume = this->GetValue(MP_MASS)/this->GetValue(MP_DENSITY);
    
    this->SetValue(MP_DENSITY, MP_Density);
    this->SetValue(MP_VOLUME, MP_Volume);
    
    
    
    
    if ( rLocalSystem.CalculationFlags.Is(UpdatedLagrangianQuadrilateral::COMPUTE_LHS_MATRIX) ) //calculation of the matrix is required
    {
    
    //contributions to stiffness matrix calculated on the reference config
    this->CalculateAndAddLHS ( rLocalSystem, Variables, MP_Volume );

    }

    if ( rLocalSystem.CalculationFlags.Is(UpdatedLagrangianQuadrilateral::COMPUTE_RHS_VECTOR) ) //calculation of the vector is required
    {
        //contribution to external forces
        VolumeForce  = this->CalculateVolumeForce( VolumeForce, Variables );

        this->CalculateAndAddRHS ( rLocalSystem, Variables, VolumeForce, MP_Volume );

    }

    KRATOS_CATCH( "" )
}
//*********************************COMPUTE KINEMATICS*********************************
//************************************************************************************


void UpdatedLagrangianQuadrilateral::CalculateKinematics(GeneralVariables& rVariables, ProcessInfo& rCurrentProcessInfo)

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
    //if (this->Id() == 541 || this->Id() == 534 || this->Id() == 538)
    //{
    ////std::cout<<" AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA "<<std::endl;
    ////std::cout<<"rVariables.CurrentDisp in calculate kinematic "<<rVariables.CurrentDisp<<std::endl;
    ////std::cout<<" AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA "<<std::endl;
    //std::cout<<" Delta Deformation Gradient "<< rVariables.F <<std::endl;
    //}





    //Determinant of the Deformation Gradient F_n

    rVariables.detF0 = mDeterminantF0;
    rVariables.F0    = mDeformationGradientF0;

    //if(this->Id() == 365)
    //{

    //std::cout<<"rVariables.DN_DX "<<this->Id()<<rVariables.DN_DX<<std::endl;
    //std::cout<<"rVariables.DN_De "<<this->Id()<<rVariables.DN_De<<std::endl;
    //std::cout<<"rVariables.J "<<this->Id()<<rVariables.J<<std::endl;
    //std::cout<<"rVariables.j "<<this->Id()<<rVariables.j<<std::endl;
    //std::cout<<"Invj "<<this->Id()<<Invj<<std::endl;
    //}

    //Compute the deformation matrix B
    this->CalculateDeformationMatrix(rVariables.B, rVariables.F, rVariables.DN_DX);


    KRATOS_CATCH( "" )
}
//************************************************************************************

void UpdatedLagrangianQuadrilateral::CalculateDeformationMatrix(Matrix& rB,
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
        //if(this->Id() == 365)
        //{
        //std::cout<<"rB "<< this->Id()<< rB<<std::endl;
        //}

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
//************************************************************************************
//************************************************************************************

void UpdatedLagrangianQuadrilateral::CalculateAndAddRHS(LocalSystemComponents& rLocalSystem, GeneralVariables& rVariables, Vector& rVolumeForce, double& rIntegrationWeight)
{

    //contribution of the internal and external forces
    if( rLocalSystem.CalculationFlags.Is( UpdatedLagrangianQuadrilateral::COMPUTE_RHS_VECTOR_WITH_COMPONENTS ) )
    {

        std::vector<VectorType>& rRightHandSideVectors = rLocalSystem.GetRightHandSideVectors();
        const std::vector< Variable< VectorType > >& rRightHandSideVariables = rLocalSystem.GetRightHandSideVariables();
        for( unsigned int i=0; i<rRightHandSideVariables.size(); i++ )
        {
            bool calculated = false;
            if( rRightHandSideVariables[i] == EXTERNAL_FORCES_VECTOR )
            {
                // operation performed: rRightHandSideVector += ExtForce*IntToReferenceWeight
                this->CalculateAndAddExternalForces( rRightHandSideVectors[i], rVariables, rVolumeForce, rIntegrationWeight );
                calculated = true;
            }

            if( rRightHandSideVariables[i] == INTERNAL_FORCES_VECTOR )
            {
                // operation performed: rRightHandSideVector -= IntForce*IntToReferenceWeight
                this->CalculateAndAddInternalForces( rRightHandSideVectors[i], rVariables, rIntegrationWeight );
                calculated = true;
            }

            if(calculated == false)
            {
                KRATOS_THROW_ERROR( std::logic_error, " ELEMENT can not supply the required local system variable: ", rRightHandSideVariables[i] )
            }

        }
    }
    else
    {

        VectorType& rRightHandSideVector = rLocalSystem.GetRightHandSideVector();

        // operation performed: rRightHandSideVector += ExtForce*IntToReferenceWeight
        this->CalculateAndAddExternalForces( rRightHandSideVector, rVariables, rVolumeForce, rIntegrationWeight );

        // operation performed: rRightHandSideVector -= IntForce*IntToReferenceWeight
        this->CalculateAndAddInternalForces( rRightHandSideVector, rVariables, rIntegrationWeight );
        //KRATOS_WATCH( rRightHandSideVector )
    }

}
//************************************************************************************
//*********************Calculate the contribution of external force*******************

void UpdatedLagrangianQuadrilateral::CalculateAndAddExternalForces(VectorType& rRightHandSideVector,
        GeneralVariables& rVariables,
        Vector& rVolumeForce,
        double& rIntegrationWeight)

{
    KRATOS_TRY
    unsigned int number_of_nodes = GetGeometry().PointsNumber();

    unsigned int dimension = GetGeometry().WorkingSpaceDimension();



    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        int index = dimension * i;

        for ( unsigned int j = 0; j < dimension; j++ )
        {
            rRightHandSideVector[index + j] += rVariables.N[i] * rVolumeForce[j];

        }

    }
    if(this->Id() == 875)
    {
        std::cout<<"rRightHandSideVector "<<rRightHandSideVector<<std::endl;
    }

    KRATOS_CATCH( "" )
}
//************************************************************************************
//************************************************************************************

void UpdatedLagrangianQuadrilateral::CalculateAndAddInternalForces(VectorType& rRightHandSideVector,
        GeneralVariables & rVariables,
        double& rIntegrationWeight)
{
    KRATOS_TRY

    VectorType InternalForces = rIntegrationWeight * prod( trans( rVariables.B ), rVariables.StressVector );
    noalias( rRightHandSideVector ) -= InternalForces;




    KRATOS_CATCH( "" )
}
//************************************************************************************
//************************************************************************************

void UpdatedLagrangianQuadrilateral::CalculateAndAddLHS(LocalSystemComponents& rLocalSystem, GeneralVariables& rVariables, double& rIntegrationWeight)
{

    //contributions of the stiffness matrix calculated on the reference configuration
    if( rLocalSystem.CalculationFlags.Is( UpdatedLagrangianQuadrilateral::COMPUTE_LHS_MATRIX_WITH_COMPONENTS ) )
    {
        std::vector<MatrixType>& rLeftHandSideMatrices = rLocalSystem.GetLeftHandSideMatrices();
        const std::vector< Variable< MatrixType > >& rLeftHandSideVariables = rLocalSystem.GetLeftHandSideVariables();

        for( unsigned int i=0; i<rLeftHandSideVariables.size(); i++ )
        {
            bool calculated = false;
            if( rLeftHandSideVariables[i] == MATERIAL_STIFFNESS_MATRIX )
            {
                // operation performed: add Km to the rLefsHandSideMatrix
                this->CalculateAndAddKuum( rLeftHandSideMatrices[i], rVariables, rIntegrationWeight );
                calculated = true;
            }

            if( rLeftHandSideVariables[i] == GEOMETRIC_STIFFNESS_MATRIX )
            {
                // operation performed: add Kg to the rLefsHandSideMatrix
                this->CalculateAndAddKuug( rLeftHandSideMatrices[i], rVariables, rIntegrationWeight );
                calculated = true;
            }

            if(calculated == false)
            {
                KRATOS_THROW_ERROR(std::logic_error, " ELEMENT can not supply the required local system variable: ",rLeftHandSideVariables[i])
            }

        }
    }
    else
    {

        MatrixType& rLeftHandSideMatrix = rLocalSystem.GetLeftHandSideMatrix();

        this->CalculateAndAddKuum( rLeftHandSideMatrix, rVariables, rIntegrationWeight );

        this->CalculateAndAddKuug( rLeftHandSideMatrix, rVariables, rIntegrationWeight );

    }


}
//************************************************************************************
//************************************************************************************

void UpdatedLagrangianQuadrilateral::CalculateAndAddKuum(MatrixType& rLeftHandSideMatrix,
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



    noalias( rLeftHandSideMatrix ) += prod( trans( rVariables.B ),  rIntegrationWeight * Matrix( prod( rVariables.ConstitutiveMatrix, rVariables.B ) ) );


    //std::cout << ss.str();

    KRATOS_CATCH( "" )
}



//************************************************************************************
//************************************************************************************

void UpdatedLagrangianQuadrilateral::CalculateAndAddKuug(MatrixType& rLeftHandSideMatrix,
        GeneralVariables& rVariables,
        double& rIntegrationWeight)

{
    KRATOS_TRY

    unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    Matrix StressTensor = MathUtils<double>::StressVectorToTensor( rVariables.StressVector );
    Matrix ReducedKg = prod( rVariables.DN_DX, rIntegrationWeight * Matrix( prod( StressTensor, trans( rVariables.DN_DX ) ) ) ); //to be optimized
    MathUtils<double>::ExpandAndAddReducedMatrix( rLeftHandSideMatrix, ReducedKg, dimension );



    KRATOS_CATCH( "" )
}
//************************************CALCULATE VOLUME CHANGE*************************
//************************************************************************************

double& UpdatedLagrangianQuadrilateral::CalculateVolumeChange( double& rVolumeChange, GeneralVariables& rVariables )
{
    KRATOS_TRY

    rVolumeChange = 1.0 / (rVariables.detF * rVariables.detF0);

    return rVolumeChange;

    KRATOS_CATCH( "" )
}
//************************************CALCULATE VOLUME ACCELERATION*******************
//************************************************************************************

Vector& UpdatedLagrangianQuadrilateral::CalculateVolumeForce( Vector& rVolumeForce, GeneralVariables& rVariables )
{
    KRATOS_TRY

    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    rVolumeForce = ZeroVector(dimension);

    //if (this->Id() == 18274)
    //{
    //std::cout<<"volume acceleration "<<this->GetValue(MP_VOLUME_ACCELERATION)<<std::endl;
    //}
    rVolumeForce = this->GetValue(MP_VOLUME_ACCELERATION)* this->GetValue(MP_MASS);


    return rVolumeForce;

    KRATOS_CATCH( "" )
}
//************************************************************************************
//************************************************************************************
void UpdatedLagrangianQuadrilateral::CalculateRightHandSide( VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo )
{
    //create local system components
    LocalSystemComponents LocalSystem;

    //calculation flags
    LocalSystem.CalculationFlags.Set(UpdatedLagrangianQuadrilateral::COMPUTE_RHS_VECTOR);

    MatrixType LeftHandSideMatrix = Matrix();

    //Initialize sizes for the system components:
    this->InitializeSystemMatrices( LeftHandSideMatrix, rRightHandSideVector, LocalSystem.CalculationFlags );

    //Set Variables to Local system components
    LocalSystem.SetLeftHandSideMatrix(LeftHandSideMatrix);
    LocalSystem.SetRightHandSideVector(rRightHandSideVector);

    //Calculate elemental system
    CalculateElementalSystem( LocalSystem, rCurrentProcessInfo );
}

//************************************************************************************
//************************************************************************************


void UpdatedLagrangianQuadrilateral::CalculateRightHandSide( std::vector< VectorType >& rRightHandSideVectors, const std::vector< Variable< VectorType > >& rRHSVariables, ProcessInfo& rCurrentProcessInfo )
{
    //create local system components
    LocalSystemComponents LocalSystem;

    //calculation flags
    LocalSystem.CalculationFlags.Set(UpdatedLagrangianQuadrilateral::COMPUTE_RHS_VECTOR);
    LocalSystem.CalculationFlags.Set(UpdatedLagrangianQuadrilateral::COMPUTE_RHS_VECTOR_WITH_COMPONENTS);

    MatrixType LeftHandSideMatrix = Matrix();

    //Initialize sizes for the system components:
    if( rRHSVariables.size() != rRightHandSideVectors.size() )
        rRightHandSideVectors.resize(rRHSVariables.size());

    for( unsigned int i=0; i<rRightHandSideVectors.size(); i++ )
    {
        this->InitializeSystemMatrices( LeftHandSideMatrix, rRightHandSideVectors[i], LocalSystem.CalculationFlags );
    }

    //Set Variables to Local system components
    LocalSystem.SetLeftHandSideMatrix(LeftHandSideMatrix);
    LocalSystem.SetRightHandSideVectors(rRightHandSideVectors);

    LocalSystem.SetRightHandSideVariables(rRHSVariables);

    //Calculate elemental system
    CalculateElementalSystem( LocalSystem, rCurrentProcessInfo );
}


//************************************************************************************
//************************************************************************************


void UpdatedLagrangianQuadrilateral::CalculateLeftHandSide( MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo )
{
    //create local system components
    LocalSystemComponents LocalSystem;

    //calculation flags
    LocalSystem.CalculationFlags.Set(UpdatedLagrangianQuadrilateral::COMPUTE_LHS_MATRIX);

    VectorType RightHandSideVector = Vector();

    //Initialize sizes for the system components:
    this->InitializeSystemMatrices( rLeftHandSideMatrix, RightHandSideVector, LocalSystem.CalculationFlags );

    //Set Variables to Local system components
    LocalSystem.SetLeftHandSideMatrix(rLeftHandSideMatrix);
    LocalSystem.SetRightHandSideVector(RightHandSideVector);

    //Calculate elemental system
    CalculateElementalSystem( LocalSystem, rCurrentProcessInfo );

}
//************************************************************************************
//************************************************************************************


void UpdatedLagrangianQuadrilateral::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo )
{

    //create local system components
    LocalSystemComponents LocalSystem;

    //calculation flags
    LocalSystem.CalculationFlags.Set(UpdatedLagrangianQuadrilateral::COMPUTE_LHS_MATRIX);
    LocalSystem.CalculationFlags.Set(UpdatedLagrangianQuadrilateral::COMPUTE_RHS_VECTOR);

    //Initialize sizes for the system components:
    this->InitializeSystemMatrices( rLeftHandSideMatrix, rRightHandSideVector, LocalSystem.CalculationFlags );

    //Set Variables to Local system components
    LocalSystem.SetLeftHandSideMatrix(rLeftHandSideMatrix);
    LocalSystem.SetRightHandSideVector(rRightHandSideVector);

    //Calculate elemental system

    CalculateElementalSystem( LocalSystem, rCurrentProcessInfo );
    //std::cout<<" in CalculateLocalSystem ends"<<std::endl;


}


//************************************************************************************
//************************************************************************************

void UpdatedLagrangianQuadrilateral::CalculateLocalSystem( std::vector< MatrixType >& rLeftHandSideMatrices,
        const std::vector< Variable< MatrixType > >& rLHSVariables,
        std::vector< VectorType >& rRightHandSideVectors,
        const std::vector< Variable< VectorType > >& rRHSVariables,
        ProcessInfo& rCurrentProcessInfo )
{
    //create local system components
    LocalSystemComponents LocalSystem;

    //calculation flags
    LocalSystem.CalculationFlags.Set(UpdatedLagrangianQuadrilateral::COMPUTE_LHS_MATRIX_WITH_COMPONENTS);
    LocalSystem.CalculationFlags.Set(UpdatedLagrangianQuadrilateral::COMPUTE_RHS_VECTOR_WITH_COMPONENTS);


    //Initialize sizes for the system components:
    if( rLHSVariables.size() != rLeftHandSideMatrices.size() )
        rLeftHandSideMatrices.resize(rLHSVariables.size());

    if( rRHSVariables.size() != rRightHandSideVectors.size() )
        rRightHandSideVectors.resize(rRHSVariables.size());

    LocalSystem.CalculationFlags.Set(UpdatedLagrangianQuadrilateral::COMPUTE_LHS_MATRIX);
    for( unsigned int i=0; i<rLeftHandSideMatrices.size(); i++ )
    {
        //Note: rRightHandSideVectors.size() > 0
        this->InitializeSystemMatrices( rLeftHandSideMatrices[i], rRightHandSideVectors[0], LocalSystem.CalculationFlags );
    }

    LocalSystem.CalculationFlags.Set(UpdatedLagrangianQuadrilateral::COMPUTE_RHS_VECTOR);
    LocalSystem.CalculationFlags.Set(UpdatedLagrangianQuadrilateral::COMPUTE_LHS_MATRIX,false);

    for( unsigned int i=0; i<rRightHandSideVectors.size(); i++ )
    {
        //Note: rLeftHandSideMatrices.size() > 0
        this->InitializeSystemMatrices( rLeftHandSideMatrices[0], rRightHandSideVectors[i], LocalSystem.CalculationFlags );
    }
    LocalSystem.CalculationFlags.Set(UpdatedLagrangianQuadrilateral::COMPUTE_LHS_MATRIX,true);


    //Set Variables to Local system components
    LocalSystem.SetLeftHandSideMatrices(rLeftHandSideMatrices);
    LocalSystem.SetRightHandSideVectors(rRightHandSideVectors);

    LocalSystem.SetLeftHandSideVariables(rLHSVariables);
    LocalSystem.SetRightHandSideVariables(rRHSVariables);

    //Calculate elemental system
    CalculateElementalSystem( LocalSystem, rCurrentProcessInfo );

}


////************************************************************************************
////************************************************************************************

void UpdatedLagrangianQuadrilateral::InitializeSolutionStep( ProcessInfo& rCurrentProcessInfo )
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
    array_1d<double,3>& AUX_MP_Velocity = this->GetValue(AUX_MP_VELOCITY);
    array_1d<double,3>& AUX_MP_Acceleration = this->GetValue(AUX_MP_ACCELERATION);
    double MP_Mass = this->GetValue(MP_MASS);
    array_1d<double,3> MP_Momentum;
    array_1d<double,3> MP_Inertia;
    array_1d<double,3> NodalMomentum;
    array_1d<double,3> NodalInertia;

    for (unsigned int j=0; j<number_of_nodes; j++)
    {
        //these are the values of nodal velocity and nodal acceleration evaluated in the initialize solution step
        array_1d<double, 3 > & NodalAcceleration = GetGeometry()[j].FastGetSolutionStepValue(ACCELERATION,1);
        array_1d<double, 3 > & NodalVelocity = GetGeometry()[j].FastGetSolutionStepValue(VELOCITY,1);

        //std::cout<<"NodalVelocity "<< GetGeometry()[j].Id()<<std::endl;
        for (unsigned int k = 0; k < dimension; k++)
        {
            AUX_MP_Velocity[k] += Variables.N[j] * NodalVelocity[k];
            AUX_MP_Acceleration[k] += Variables.N[j] * NodalAcceleration[k];
        }
    }
    //if(this->Id() == 8371)
    //{
    //std::cout<<" AUX_MP_Velocity " <<  AUX_MP_Velocity<<std::endl;
    //std::cout<<" AUX_MP_Acceleration " <<  AUX_MP_Acceleration<<std::endl;
    //}
    //std::cout<<" AUX_MP_Velocity "<<AUX_MP_Velocity<<std::endl;
    //std::cout<<" AUX_MP_Acceleration "<<AUX_MP_Acceleration<<std::endl;


    //for (unsigned int k = 0; k < 3; k++)
    //{
    //MP_Momentum[k] = (MP_Velocity[k] - AUX_MP_Velocity[k]) * MP_Mass;
    //MP_Inertia[k] = (MP_Acceleration[k] - AUX_MP_Acceleration[k]) * MP_Mass;
    //}

    // Here MP contribution in terms of momentum, inertia and mass are added
    //if(this->Id() == 12156)
    //{
    //std::cout<<" AUX_MP_Velocity"<< AUX_MP_Velocity<<std::endl;
    //std::cout<<" AUX_MP_Acceleration"<< AUX_MP_Acceleration<<std::endl;
    //std::cout<<" error velocity"<< MP_Velocity - AUX_MP_Velocity<<std::endl;
    //std::cout<<" error acceleration"<< MP_Acceleration - AUX_MP_Acceleration<<std::endl;
    //}
    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        for (unsigned int j = 0; j < dimension; j++)
        {
            NodalMomentum[j] = Variables.N[i] * (MP_Velocity[j] - AUX_MP_Velocity[j]) * MP_Mass;
            NodalInertia[j] = Variables.N[i] * (MP_Acceleration[j] - AUX_MP_Acceleration[j]) * MP_Mass;

        }
        GetGeometry()[i].GetSolutionStepValue(NODAL_MOMENTUM, 0) += NodalMomentum;
        GetGeometry()[i].GetSolutionStepValue(NODAL_INERTIA, 0) += NodalInertia;
        //if(GetGeometry()[i].pGetDof(DISPLACEMENT_X)->IsFixed() == false)
        //{
        GetGeometry()[i].GetSolutionStepValue(NODAL_MASS, 0) += Variables.N[i] * MP_Mass;
        //}
    }
    //for ( unsigned int i = 0; i < number_of_nodes; i++ )

    //{
    //if(GetGeometry()[i].pGetDof(DISPLACEMENT_X)->IsFixed() == false)
    //{
    //for (unsigned int j = 0; j < number_of_nodes; j++)
    //{
    //if(GetGeometry()[j].pGetDof(DISPLACEMENT_X)->IsFixed() == false)
    //{
    //GetGeometry()[i].GetSolutionStepValue(NODAL_MASS, 0) += Variables.N[i] * Variables.N[j] * MP_Mass;
    //}
    //}
    //}
    //}
    //MP_Velocity = MP_Velocity - AUX_MP_Velocity;
    //MP_Acceleration = MP_Acceleration - AUX_MP_Acceleration;
    //if(this->Id() == 8371)
    //{
    //std::cout<<" MP_Velocity " <<  MP_Velocity<<std::endl;
    //std::cout<<" MP_Acceleration " <<  MP_Acceleration<<std::endl;
    //}
    AUX_MP_Velocity.clear();
    AUX_MP_Acceleration.clear();



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

    //if(this->Id() == 6261 || this->Id() == 1836)
    //{
    //std::cout<<"density "<<this->Id() << GetProperties()[DENSITY]<<std::endl;
    //}

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
    //if(this->Id() == 8371)
    //{
    //std::cout<<" MP_Velocity at iteration 0"<< MP_Velocity<<std::endl;
    //std::cout<<" MP_Acceleration at iteration 0"<< MP_Acceleration<<std::endl;

    //}
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
    ////if(this->Id() == 18274)
    ////{
    ////std::cout<<"GetGeometry()[i].Id()"<<GetGeometry()[i].Id()<<std::endl;
    ////std::cout<<"GetGeometry()[i].GetSolutionStepValue(NODAL_MASS, 0)"<<GetGeometry()[i].GetSolutionStepValue(NODAL_MASS, 0)<<std::endl;
    ////std::cout<<"GetGeometry()[i].GetSolutionStepValue(NODAL_MOMENTUM, 0) "<<GetGeometry()[i].GetSolutionStepValue(NODAL_MOMENTUM, 0) <<std::endl;
    ////std::cout<<"GetGeometry()[i].GetSolutionStepValue(NODAL_INERTIA, 0) "<<GetGeometry()[i].GetSolutionStepValue(NODAL_INERTIA, 0) <<std::endl;
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
    //if(this->Id() == 8371)
    //{
    //std::cout<<" AUX_MP_Velocity " <<  AUX_MP_Velocity<<std::endl;
    //std::cout<<" AUX_MP_Acceleration " <<  AUX_MP_Acceleration<<std::endl;
    //}
    ////std::cout<<" AUX_MP_Velocity "<<AUX_MP_Velocity<<std::endl;
    ////std::cout<<" AUX_MP_Acceleration "<<AUX_MP_Acceleration<<std::endl;


    ////for (unsigned int k = 0; k < 3; k++)
    ////{
    ////MP_Momentum[k] = (MP_Velocity[k] - AUX_MP_Velocity[k]) * MP_Mass;
    ////MP_Inertia[k] = (MP_Acceleration[k] - AUX_MP_Acceleration[k]) * MP_Mass;
    ////}

    //// Here MP contribution in terms of momentum, inertia and mass are added
    ////if(this->Id() == 12156)
    ////{
    ////std::cout<<" AUX_MP_Velocity"<< AUX_MP_Velocity<<std::endl;
    ////std::cout<<" AUX_MP_Acceleration"<< AUX_MP_Acceleration<<std::endl;
    ////std::cout<<" error velocity"<< MP_Velocity - AUX_MP_Velocity<<std::endl;
    ////std::cout<<" error acceleration"<< MP_Acceleration - AUX_MP_Acceleration<<std::endl;
    ////}
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
    ////if(this->Id() == 8371)
    ////{
    ////std::cout<<" MP_Velocity " <<  MP_Velocity<<std::endl;
    ////std::cout<<" MP_Acceleration " <<  MP_Acceleration<<std::endl;
    ////}
    //AUX_MP_Velocity.clear();
    //AUX_MP_Acceleration.clear();
    //}


}

////************************************************************************************
////************************************************************************************

void UpdatedLagrangianQuadrilateral::IterativeExtrapolation( ProcessInfo& rCurrentProcessInfo )
{
    // In the Initialize of each time step the nodal initial conditions are evaluated
    //1. first of all I need to evaluate the MP momentum and MP_inertia

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

    //mConstitutiveLawVector->InitializeSolutionStep( GetProperties(),
    //        GetGeometry(), Variables.N, rCurrentProcessInfo );

    //mFinalizedStep = false;



    array_1d<double,3>& MP_Velocity = this->GetValue(MP_VELOCITY);
    array_1d<double,3>& MP_Acceleration = this->GetValue(MP_ACCELERATION);
    array_1d<double,3>& AUX_MP_Velocity = this->GetValue(AUX_MP_VELOCITY);
    array_1d<double,3>& AUX_MP_Acceleration = this->GetValue(AUX_MP_ACCELERATION);
    double MP_Mass = this->GetValue(MP_MASS);
    array_1d<double,3> MP_Momentum;
    array_1d<double,3> MP_Inertia;
    array_1d<double,3> NodalMomentum;
    array_1d<double,3> NodalInertia;

    for (unsigned int j=0; j<number_of_nodes; j++)
    {
        //these are the values of nodal velocity and nodal acceleration evaluated in the initialize solution step
        array_1d<double, 3 > & NodalAcceleration = GetGeometry()[j].FastGetSolutionStepValue(ACCELERATION,1);
        array_1d<double, 3 > & NodalVelocity = GetGeometry()[j].FastGetSolutionStepValue(VELOCITY,1);
        for (unsigned int k = 0; k < dimension; k++)
        {
            AUX_MP_Velocity[k] += Variables.N[j] * NodalVelocity[k];
            AUX_MP_Acceleration[k] += Variables.N[j] * NodalAcceleration[k];
        }
    }



    for (unsigned int k = 0; k < dimension; k++)
    {
        MP_Momentum[k] = (MP_Velocity[k] - AUX_MP_Velocity[k]) * MP_Mass;
        MP_Inertia[k] = (MP_Acceleration[k] - AUX_MP_Acceleration[k]) * MP_Mass;
    }

    // Here MP contribution in terms of momentum, inertia and mass are added

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        for (unsigned int j = 0; j < dimension; j++)
        {
            NodalMomentum[j] = Variables.N[i] * MP_Momentum[j];
            NodalInertia[j] = Variables.N[i] * MP_Inertia[j];

        }
        GetGeometry()[i].GetSolutionStepValue(NODAL_MOMENTUM, 0) += NodalMomentum;
        GetGeometry()[i].GetSolutionStepValue(NODAL_INERTIA, 0) += NodalInertia;
        //GetGeometry()[i].GetSolutionStepValue(NODAL_MASS, 0) += Variables.N[i] * MP_Mass;


    }
}
////************************************************************************************
////************************************************************************************
void UpdatedLagrangianQuadrilateral::InitializeNonLinearIteration( ProcessInfo& rCurrentProcessInfo )
{

}

////************************************************************************************
////************************************************************************************

void UpdatedLagrangianQuadrilateral::FinalizeNonLinearIteration( ProcessInfo& rCurrentProcessInfo )
{

}

////************************************************************************************
////************************************************************************************

void UpdatedLagrangianQuadrilateral::FinalizeSolutionStep( ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY
    //const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    //const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    //unsigned int voigtsize  = 3;

    //if( dimension == 3 )
    //{
    //voigtsize  = 6;
    //}
    //Vector NodalStress = ZeroVector(voigtsize);
    array_1d<double,3> NodalStress;
    //create and initialize element variables:
    GeneralVariables Variables;
    this->InitializeGeneralVariables(Variables,rCurrentProcessInfo);

    //create constitutive law parameters:
    ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

    //set constitutive law flags:
    Flags &ConstitutiveLawOptions=Values.GetOptions();

    ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);
    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS);
    //ConstitutiveLawOptions.Set(ConstitutiveLaw::ISOCHORIC_TENSOR_ONLY);
    //compute element kinematics B, F, DN_DX ...
    this->CalculateKinematics(Variables, rCurrentProcessInfo);

    //set general variables to constitutivelaw parameters
    this->SetGeneralVariables(Variables,Values);

    //call the constitutive law to update material variables
    mConstitutiveLawVector->FinalizeMaterialResponse(Values, Variables.StressMeasure);

    //call the constitutive law to finalize the solution step
    mConstitutiveLawVector->FinalizeSolutionStep( GetProperties(),
            GetGeometry(),
            Variables.N,
            rCurrentProcessInfo );
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
    //call the element internal variables update
    this->FinalizeStepVariables(Variables, rCurrentProcessInfo);


    mFinalizedStep = true;

    KRATOS_CATCH( "" )
}


////************************************************************************************************************
//void UpdatedLagrangianQuadrilateral::Calculate(const Variable<Vector >& rVariable, Vector& Output, const ProcessInfo& rCurrentProcessInfo)
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
////if (this->Id() == 541 || this->Id() == 534 || this->Id() == 538)
////{
////std::cout<<" GetGeometry()[0].GetSolutionStepValue(STRESSES, 0); "<< " Id "<< GetGeometry()[0].GetSolutionStepValue(STRESSES, 0)<<std::endl;
////}

//for ( unsigned int i = 0; i < number_of_nodes; i++ )
//{
//array_1d<double,3> NodalStress = GetGeometry()[i].GetSolutionStepValue(NODAL_STRESSES, 0);
//for ( unsigned int j = 0; j< voigtsize; j++)
//{
//SmoothMPStress[j] += NodalStress[j] * Variables.N[i];
//}
//}
//this->SetValue(MP_CAUCHY_STRESS_VECTOR, SmoothMPStress);
////if (this->Id() == 541 || this->Id() == 534 || this->Id() == 538)
////{
////std::cout<<" MP_CAUCHY_STRESS_VECTOR "<< " Id "<< this->GetValue(MP_CAUCHY_STRESS_VECTOR)<<std::endl;
////}

//}
//KRATOS_CATCH( "" )
//}

////************************************************************************************
////************************************************************************************

void UpdatedLagrangianQuadrilateral::FinalizeStepVariables( GeneralVariables & rVariables, const ProcessInfo& rCurrentProcessInfo)
{
    //update internal (historical) variables
    mDeterminantF0         = rVariables.detF* rVariables.detF0;
    mDeformationGradientF0 = prod(rVariables.F, rVariables.F0);

    this->SetValue(MP_CAUCHY_STRESS_VECTOR, rVariables.StressVector);
    this->SetValue(MP_ALMANSI_STRAIN_VECTOR, rVariables.StrainVector);
    //std::cout<<" before get equivalent plastic strain "<<std::endl;
    double EquivalentPlasticStrain = mConstitutiveLawVector->GetValue(PLASTIC_STRAIN, EquivalentPlasticStrain );
    ////std::cout<<" EquivalentPlasticStrain in the element "<<EquivalentPlasticStrain<<std::endl;
    this->SetValue(MP_EQUIVALENT_PLASTIC_STRAIN, EquivalentPlasticStrain);



    ////double DeltaEquivalentPlasticStrain = mConstitutiveLawVector->GetValue(DELTA_PLASTIC_STRAIN, DeltaEquivalentPlasticStrain );
    ////this->SetValue(MP_EQUIVALENT_DELTA_PLASTIC_STRAIN, DeltaEquivalentPlasticStrain);

    //ComparisonUtilities EquivalentStress;
    //double MPMStressNorm = EquivalentStress.CalculateStressNorm(rVariables.StressVector);
    ////std::cout<<" MPMStressNorm "<<MPMStressNorm<<std::endl;
    //this->SetValue(MPM_NORM_ISOCHORIC_STRESS, MPMStressNorm);

    MathUtils<double>::InvertMatrix( rVariables.j, mInverseJ, rVariables.detJ );
    this->UpdateGaussPoint(rVariables, rCurrentProcessInfo);

}

//************************************************************************************
//************************************************************************************
/**
 * The position of the Gauss points/Material points is updated
 */

void UpdatedLagrangianQuadrilateral::UpdateGaussPoint( GeneralVariables & rVariables, const ProcessInfo& rCurrentProcessInfo)
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
        if (rVariables.N[i] != 0.0)
        {
            array_1d<double, 3 > & NodalAcceleration = GetGeometry()[i].FastGetSolutionStepValue(ACCELERATION);
            array_1d<double, 3 > & NodalVelocity = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY);
            double NodalMass = GetGeometry()[i].GetSolutionStepValue(NODAL_MASS, 0);
            array_1d<double,3> NodalMomentum = NodalMass * NodalVelocity;
            array_1d<double,3> NodalInertia = NodalMass * NodalAcceleration;

            //if (this->Id() == 469)// || this->Id() == 1513)
            //{
            //std::cout<< "Nodal ID "<< GetGeometry()[i].Id()<<std::endl;
            //std::cout<< "NodalAcceleration "<<NodalAcceleration<<std::endl;
            //std::cout<< "NodalVelocity "<<NodalVelocity<<std::endl;
            //std::cout<< "NodalMass "<<NodalMass<<std::endl;

            //std::cout<< "rVariables.N "<<rVariables.N<<std::endl;
            //}




            for ( unsigned int j = 0; j < dimension; j++ )
            {

                delta_xg[j] += rVariables.N[i] * rVariables.CurrentDisp(i,j);
                MP_Acceleration[j] += rVariables.N[i] * NodalAcceleration[j];
                //MP_Velocity[j] += rVariables.N[i] * NodalVelocity[j];

                //MP_Acceleration[j] +=NodalInertia[j]/(rVariables.N[i] * MP_Mass * MP_number);//
                //MP_Velocity[j] += NodalMomentum[j]/(rVariables.N[i] * MP_Mass * MP_number);
                //MP_Velocity[j] += DeltaTime * rVariables.N[i] * NodalAcceleration[j];////




            }
            //if (this->Id() == 14716)// || this->Id() == 1513)
            //{
            //std::cout<< "Nodal ID "<< GetGeometry()[i].Id()<<std::endl;
            //std::cout<< "rVariables.CurrentDisp(i,0) "<<rVariables.CurrentDisp(i,0)<<std::endl;
            //std::cout<< "rVariables.N[i] "<<rVariables.N[i]<<std::endl;
            //std::cout<< "delta_xg[0] "<<delta_xg[0]<<std::endl;
            //}
        }

    }


    //**************************************************************************************************************************
    //Another way to update the MP velocity (see paper Guilkey and Weiss, 2003) !!!USING THIS EXPRESSION I CONSERVE MORE ENERGY
    MP_Velocity = MP_PreviousVelocity + 0.5 * DeltaTime * (MP_Acceleration + MP_PreviousAcceleration);
    //MP_Velocity += MP_PreviousVelocity;
    //Update the MP Velocity

    //MP_Acceleration = 4/(DeltaTime * DeltaTime) * delta_xg - 4/DeltaTime * MP_PreviousVelocity;
    //MP_Velocity = 2/DeltaTime * delta_xg - MP_PreviousVelocity;
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



    //if (this->Id() == 14716)// || this->Id() == 1513)

    //{
    //std::cout<<" rVariables.N "<<this->Id()<<rVariables.N<<std::endl;
    //std::cout<<" GetGeometry() "<<this->Id()<<GetGeometry()<<std::endl;

    //std::cout<<" delta_xg "<<this->Id()<<delta_xg<<std::endl;

    //std::cout<<" MP_Velocity "<<this->Id()<<this -> GetValue(MP_VELOCITY)<<std::endl;

    //std::cout<<" MP_Acceleration "<<this->Id()<<this -> GetValue(MP_ACCELERATION)<<std::endl;

    //}

    KRATOS_CATCH( "" )
}



void UpdatedLagrangianQuadrilateral::InitializeMaterial()
{
    KRATOS_TRY
    array_1d<double,3>& xg = this->GetValue(GAUSS_COORD);
    GeneralVariables Variables;
    //this->InitializeGeneralVariables(Variables,rCurrentProcessInfo);


    if ( GetProperties()[CONSTITUTIVE_LAW] != NULL )
    {

        mConstitutiveLawVector = GetProperties()[CONSTITUTIVE_LAW]->Clone();


        Variables.N = this->MPMShapeFunctionPointValues(Variables.N, xg);

        mConstitutiveLawVector->InitializeMaterial( GetProperties(), GetGeometry(),
                Variables.N );

        //}
    }
    else
        KRATOS_THROW_ERROR( std::logic_error, "a constitutive law needs to be specified for the element with ID ", this->Id() )
        //std::cout<< "in initialize material "<<std::endl;
        KRATOS_CATCH( "" )
    }


//************************************************************************************
//************************************************************************************

void UpdatedLagrangianQuadrilateral::ResetConstitutiveLaw()
{
    KRATOS_TRY
    array_1d<double,3>& xg = this->GetValue(GAUSS_COORD);
    GeneralVariables Variables;
    //this->InitializeGeneralVariables(Variables,rCurrentProcessInfo);

    //create and initialize element variables:

    if ( GetProperties()[CONSTITUTIVE_LAW] != NULL )
    {

        mConstitutiveLawVector->ResetMaterial( GetProperties(), GetGeometry(), this->MPMShapeFunctionPointValues(Variables.N, xg) );
    }

    KRATOS_CATCH( "" )
}




//*************************COMPUTE CURRENT DISPLACEMENT*******************************
//************************************************************************************


Matrix& UpdatedLagrangianQuadrilateral::CalculateCurrentDisp(Matrix & rCurrentDisp, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    rCurrentDisp = zero_matrix<double>( number_of_nodes, dimension);

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
void UpdatedLagrangianQuadrilateral::CalculateAlmansiStrain(const Matrix& rF,
        Vector& rStrainVector )
{
    KRATOS_TRY

    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    //Left Cauchy-Green Calculation
    Matrix LeftCauchyGreen = prod( rF, trans( rF ) );

    //Calculating the inverse of the jacobian
    Matrix InverseLeftCauchyGreen ( dimension, dimension );
    double det_b=0;
    MathUtils<double>::InvertMatrix( LeftCauchyGreen, InverseLeftCauchyGreen, det_b);

    if( dimension == 2 )
    {

        //Almansi Strain Calculation
        rStrainVector[0] = 0.5 * (  1.00 - InverseLeftCauchyGreen( 0, 0 ) );

        rStrainVector[1] = 0.5 * (  1.00 - InverseLeftCauchyGreen( 1, 1 ) );

        rStrainVector[2] = - InverseLeftCauchyGreen( 0, 1 ); // xy

    }
    else if( dimension == 3 )
    {

        //Almansi Strain Calculation
        if ( rStrainVector.size() != 6 ) rStrainVector.resize( 6, false );

        rStrainVector[0] = 0.5 * (  1.00 - InverseLeftCauchyGreen( 0, 0 ) );

        rStrainVector[1] = 0.5 * (  1.00 - InverseLeftCauchyGreen( 1, 1 ) );

        rStrainVector[2] = 0.5 * (  1.00 - InverseLeftCauchyGreen( 2, 2 ) );

        rStrainVector[3] = - InverseLeftCauchyGreen( 0, 1 ); // xy

        rStrainVector[4] = - InverseLeftCauchyGreen( 1, 2 ); // yz

        rStrainVector[5] = - InverseLeftCauchyGreen( 0, 2 ); // xz

    }
    else
    {

        KRATOS_THROW_ERROR( std::invalid_argument, "something is wrong with the dimension", "" );

    }


    KRATOS_CATCH( "" )
}
//*************************COMPUTE GREEN-LAGRANGE STRAIN*************************************
//************************************************************************************

void UpdatedLagrangianQuadrilateral::CalculateGreenLagrangeStrain(const Matrix& rF,
        Vector& rStrainVector )
{
    KRATOS_TRY

    const unsigned int dimension  = GetGeometry().WorkingSpaceDimension();

    //Right Cauchy-Green Calculation
    Matrix C ( dimension, dimension );
    noalias( C ) = prod( trans( rF ), rF );

    if( dimension == 2 )
    {

        //Green Lagrange Strain Calculation
        if ( rStrainVector.size() != 3 ) rStrainVector.resize( 3, false );

        rStrainVector[0] = 0.5 * ( C( 0, 0 ) - 1.00 );

        rStrainVector[1] = 0.5 * ( C( 1, 1 ) - 1.00 );

        rStrainVector[2] = C( 0, 1 ); // xy

    }
    else if( dimension == 3 )
    {

        //Green Lagrange Strain Calculation
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

        KRATOS_THROW_ERROR( std::invalid_argument, "something is wrong with the dimension", "" )

    }

    KRATOS_CATCH( "" )
}



//************************************************************************************
//************************************************************************************

double& UpdatedLagrangianQuadrilateral::CalculateIntegrationWeight(double& rIntegrationWeight)
{
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    if( dimension == 2 )
        rIntegrationWeight *= GetProperties()[THICKNESS];

    return rIntegrationWeight;
}


//************************************************************************************
//************************************************************************************

void UpdatedLagrangianQuadrilateral::EquationIdVector( EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo )
{
    int number_of_nodes = GetGeometry().size();
    int dim = GetGeometry().WorkingSpaceDimension();
    unsigned int dim2 = number_of_nodes * dim;

    if ( rResult.size() != dim2 )
        rResult.resize( dim2, false );

    for ( int i = 0; i < number_of_nodes; i++ )
    {
        int index = i * dim;
        rResult[index] = GetGeometry()[i].GetDof( DISPLACEMENT_X ).EquationId();
        rResult[index + 1] = GetGeometry()[i].GetDof( DISPLACEMENT_Y ).EquationId();

        if ( dim == 3 )
            rResult[index + 2] = GetGeometry()[i].GetDof( DISPLACEMENT_Z ).EquationId();
    }

}

//************************************************************************************
//************************************************************************************

void UpdatedLagrangianQuadrilateral::GetDofList( DofsVectorType& rElementalDofList, ProcessInfo& CurrentProcessInfo )
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
    //std::cout<< "ElementalDofList.size() "<<rElementalDofList.size()<<std::endl;
}



//************************************************************************************
//*******************DAMPING MATRIX***************************************************

void UpdatedLagrangianQuadrilateral::CalculateDampingMatrix( MatrixType& rDampingMatrix, ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    //0.-Initialize the DampingMatrix:
    unsigned int number_of_nodes = GetGeometry().size();
    unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    //resizing as needed the LHS
    unsigned int MatSize = number_of_nodes * dimension;

    if ( rDampingMatrix.size1() != MatSize )
        rDampingMatrix.resize( MatSize, MatSize, false );

    noalias( rDampingMatrix ) = ZeroMatrix( MatSize, MatSize );


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
    //std::cout<<" rDampingMatrix "<<rDampingMatrix<<std::endl;

    KRATOS_CATCH( "" )
}
//************************************************************************************
//****************MASS MATRIX*********************************************************

void UpdatedLagrangianQuadrilateral::CalculateMassMatrix( MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    //I need to call the values of the shape function for the single element
    GeneralVariables Variables;
    this->InitializeGeneralVariables(Variables,rCurrentProcessInfo);

    //lumped
    unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    unsigned int MatSize = dimension * number_of_nodes;

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

        for ( unsigned int j = 0; j < dimension; j++ )
        {
            unsigned int index = i * dimension + j;
            rMassMatrix( index, index ) = temp;
        }
    }

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



Matrix& UpdatedLagrangianQuadrilateral::MPMJacobian( Matrix& rResult, array_1d<double,3>& rPoint)
{

    KRATOS_TRY

    //derivatives of shape functions
    Matrix shape_functions_gradients;
    shape_functions_gradients =this->MPMShapeFunctionsLocalGradients(
                                   shape_functions_gradients, rPoint);
    const GeometryType& rGeom = GetGeometry();

    unsigned int number_nodes = rGeom.PointsNumber();
    unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    if (dimension ==2)
    {
        rResult.resize( 2, 2);
        rResult = ZeroMatrix(2,2);


        for ( unsigned int i = 0; i < number_nodes; i++ )
        {


            rResult( 0, 0 ) += ( GetGeometry().GetPoint( i ).X() *  shape_functions_gradients( i, 0 ) );
            rResult( 0, 1 ) += ( GetGeometry().GetPoint( i ).X() *  shape_functions_gradients( i, 1 ) );
            rResult( 1, 0 ) += ( GetGeometry().GetPoint( i ).Y() *  shape_functions_gradients( i, 0 ) );
            rResult( 1, 1 ) += ( GetGeometry().GetPoint( i ).Y() *  shape_functions_gradients( i, 1 ) );


        }

    }
    else if(dimension ==3)
    {


        rResult.resize( 3,3);
        rResult = ZeroMatrix(3,3);

        for ( unsigned int i = 0; i < number_nodes; i++ )
        {
            rResult( 0, 0 ) += ( GetGeometry().GetPoint( i ).X() *  shape_functions_gradients( i, 0 ) );
            rResult( 0, 1 ) += ( GetGeometry().GetPoint( i ).X() *  shape_functions_gradients( i, 1 ) );
            rResult( 0, 2 ) += ( GetGeometry().GetPoint( i ).X() *  shape_functions_gradients( i, 2 ) );
            rResult( 1, 0 ) += ( GetGeometry().GetPoint( i ).Y() *  shape_functions_gradients( i, 0 ) );
            rResult( 1, 1 ) += ( GetGeometry().GetPoint( i ).Y() *  shape_functions_gradients( i, 1 ) );
            rResult( 1, 2 ) += ( GetGeometry().GetPoint( i ).Y() *  shape_functions_gradients( i, 2 ) );
            rResult( 2, 0 ) += ( GetGeometry().GetPoint( i ).Z() *  shape_functions_gradients( i, 0 ) );
            rResult( 2, 1 ) += ( GetGeometry().GetPoint( i ).Z() *  shape_functions_gradients( i, 1 ) );
            rResult( 2, 2 ) += ( GetGeometry().GetPoint( i ).Z() *  shape_functions_gradients( i, 2 ) );

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
Matrix& UpdatedLagrangianQuadrilateral::MPMJacobianDelta( Matrix& rResult, array_1d<double,3>& rPoint, Matrix & rDeltaPosition )
{
    KRATOS_TRY

    Matrix shape_functions_gradients;

    shape_functions_gradients = this->MPMShapeFunctionsLocalGradients(
                                    shape_functions_gradients, rPoint );



    unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    //Elements of jacobian matrix (e.g. J(1,1) = dX1/dXi1)

    if (dimension ==2)
    {

        rResult.resize( 2, 2);
        rResult = ZeroMatrix(2,2);

        for ( unsigned int i = 0; i < GetGeometry().size(); i++ )
        {
            rResult( 0, 0 ) += ( GetGeometry().GetPoint( i ).X() + rDeltaPosition(i,0)) * ( shape_functions_gradients( i, 0 ) );
            rResult( 0, 1 ) += ( GetGeometry().GetPoint( i ).X() + rDeltaPosition(i,0)) * ( shape_functions_gradients( i, 1 ) );
            rResult( 1, 0 ) += ( GetGeometry().GetPoint( i ).Y() + rDeltaPosition(i,1)) * ( shape_functions_gradients( i, 0 ) );
            rResult( 1, 1 ) += ( GetGeometry().GetPoint( i ).Y() + rDeltaPosition(i,1)) * ( shape_functions_gradients( i, 1 ) );
        }
    }
    else if(dimension ==3)
    {

        rResult.resize( 3,3);
        rResult = ZeroMatrix(3,3);
        for ( unsigned int i = 0; i < GetGeometry().size(); i++ )
        {
            rResult( 0, 0 ) += ( GetGeometry().GetPoint( i ).X() + rDeltaPosition(i,0)) * ( shape_functions_gradients( i, 0 ) );
            rResult( 0, 1 ) += ( GetGeometry().GetPoint( i ).X() + rDeltaPosition(i,0)) * ( shape_functions_gradients( i, 1 ) );
            rResult( 0, 2 ) += ( GetGeometry().GetPoint( i ).X() + rDeltaPosition(i,0)) * ( shape_functions_gradients( i, 2 ) );
            rResult( 1, 0 ) += ( GetGeometry().GetPoint( i ).Y() + rDeltaPosition(i,1)) * ( shape_functions_gradients( i, 0 ) );
            rResult( 1, 1 ) += ( GetGeometry().GetPoint( i ).Y() + rDeltaPosition(i,1)) * ( shape_functions_gradients( i, 1 ) );
            rResult( 1, 2 ) += ( GetGeometry().GetPoint( i ).Y() + rDeltaPosition(i,1)) * ( shape_functions_gradients( i, 2 ) );
            rResult( 2, 0 ) += ( GetGeometry().GetPoint( i ).Z() + rDeltaPosition(i,2)) * ( shape_functions_gradients( i, 0 ) );
            rResult( 2, 1 ) += ( GetGeometry().GetPoint( i ).Z() + rDeltaPosition(i,2)) * ( shape_functions_gradients( i, 1 ) );
            rResult( 2, 2 ) += ( GetGeometry().GetPoint( i ).Z() + rDeltaPosition(i,2)) * ( shape_functions_gradients( i, 2 ) );
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
Vector& UpdatedLagrangianQuadrilateral::MPMShapeFunctionPointValues( Vector& rResult, array_1d<double,3>& rPoint )
{
    KRATOS_TRY

    unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    array_1d<double,3> rPointLocal = ZeroVector(dimension);
    //local coordinates of the MP position
    //I evaluate the local coordinates of the integration point
    rPointLocal = GetGeometry().PointLocalCoordinates(rPointLocal, rPoint);

    if (dimension == 2)
    {

        rResult.resize(4, false);

        //1. I evaluate the local coordinates of a point
        //rPointLocal[0] = ((GetGeometry()[2].Coordinates()[1] - GetGeometry()[0].Coordinates()[1])*(rPoint[0] - GetGeometry()[0].Coordinates()[0]) -
        //        (GetGeometry()[2].Coordinates()[0] - GetGeometry()[0].Coordinates()[0])*(rPoint[1] - GetGeometry()[0].Coordinates()[1]))/mDeterminantJ0;


        //rPointLocal[1] = (-(GetGeometry()[1].Coordinates()[1] - GetGeometry()[0].Coordinates()[1])*(rPoint[0] - GetGeometry()[0].Coordinates()[0]) +
        //        (GetGeometry()[1].Coordinates()[0] - GetGeometry()[0].Coordinates()[0])*(rPoint[1] - GetGeometry()[0].Coordinates()[1]))/mDeterminantJ0;

        //test (if the first node of the connettivity is the node at the bottom left)
        //rPointLocal[0] = (2 * rPoint[0] - GetGeometry()[2].Coordinates()[0] - GetGeometry()[0].Coordinates()[0])/(GetGeometry()[2].Coordinates()[0] - GetGeometry()[0].Coordinates()[0]);
        //rPointLocal[1] = (2 * rPoint[1] - GetGeometry()[2].Coordinates()[1] - GetGeometry()[0].Coordinates()[1])/(GetGeometry()[2].Coordinates()[1] - GetGeometry()[0].Coordinates()[1]);


        //test (if the first node of the connettivity is the node at the top left)
        //rPointLocal[0] = (2 * rPoint[0] - GetGeometry()[3].Coordinates()[0] - GetGeometry()[1].Coordinates()[0])/(GetGeometry()[3].Coordinates()[0] - GetGeometry()[1].Coordinates()[0]);
        //rPointLocal[1] = (2 * rPoint[1] - GetGeometry()[3].Coordinates()[1] - GetGeometry()[1].Coordinates()[1])/(GetGeometry()[3].Coordinates()[1] - GetGeometry()[1].Coordinates()[1]);

        //test (if the first node of the connettivity is the node at the top right)
        //rPointLocal[0] = (2 * rPoint[0] - GetGeometry()[0].Coordinates()[0] - GetGeometry()[2].Coordinates()[0])/(GetGeometry()[0].Coordinates()[0] - GetGeometry()[2].Coordinates()[0]);
        //rPointLocal[1] = (2 * rPoint[1] - GetGeometry()[0].Coordinates()[1] - GetGeometry()[2].Coordinates()[1])/(GetGeometry()[0].Coordinates()[1] - GetGeometry()[2].Coordinates()[1]);


        //test (if the first node of the connettivity is the node at the bottom right)
        //rPointLocal[0] = (2 * rPoint[0] - GetGeometry()[3].Coordinates()[0] - GetGeometry()[1].Coordinates()[0])/(GetGeometry()[1].Coordinates()[0] - GetGeometry()[3].Coordinates()[0]);
        //rPointLocal[1] = (2 * rPoint[1] - GetGeometry()[3].Coordinates()[1] - GetGeometry()[1].Coordinates()[1])/(GetGeometry()[1].Coordinates()[1] - GetGeometry()[3].Coordinates()[1]);

        //rPointLocal = GetGeometry().PointLocalCoordinates(rPointLocal, rPoint);

        //rResult = GetGeometry().ShapeFunctionsValues(rResult, rPointLocal);

        //if (this->Id() == 14716)// || this->Id() == 1513)
        //{
        //std::cout<<" *********LOCAL COORDINATES*********"<< rPointLocal<<std::endl;
        //}


        //Shape Functions (if the first node of the connettivity is the node at the bottom left)
        rResult( 0 ) = 0.25 * (1 - rPointLocal[0]) * (1 - rPointLocal[1]) ;
        rResult( 1 ) = 0.25 * (1 + rPointLocal[0]) * (1 - rPointLocal[1]) ;
        rResult( 2 ) = 0.25 * (1 + rPointLocal[0]) * (1 + rPointLocal[1]) ;
        rResult( 3 ) = 0.25 * (1 - rPointLocal[0]) * (1 + rPointLocal[1]) ;

        //Shape Function (if the first node of the connettivity is the node at the top left)
        //rResult( 0 ) = 0.25 * (1 - rPointLocal[0]) * (1 + rPointLocal[1]) ;
        //rResult( 1 ) = 0.25 * (1 - rPointLocal[0]) * (1 - rPointLocal[1]) ;
        //rResult( 2 ) = 0.25 * (1 + rPointLocal[0]) * (1 - rPointLocal[1]) ;
        //rResult( 3 ) = 0.25 * (1 + rPointLocal[0]) * (1 + rPointLocal[1]) ;

        //Shape Function (if the first node of the connettivity is the node at the top right)
        //rResult( 0 ) = 0.25 * (1 + rPointLocal[0]) * (1 + rPointLocal[1]) ;
        //rResult( 1 ) = 0.25 * (1 - rPointLocal[0]) * (1 + rPointLocal[1]) ;
        //rResult( 2 ) = 0.25 * (1 - rPointLocal[0]) * (1 - rPointLocal[1]) ;
        //rResult( 3 ) = 0.25 * (1 + rPointLocal[0]) * (1 - rPointLocal[1]) ;

        //Shape Function (if the first node of the connettivity is the node at the bottom right)
        //rResult( 0 ) = 0.25 * (1 + rPointLocal[0]) * (1 - rPointLocal[1]) ;
        //rResult( 1 ) = 0.25 * (1 + rPointLocal[0]) * (1 + rPointLocal[1]) ;
        //rResult( 2 ) = 0.25 * (1 - rPointLocal[0]) * (1 + rPointLocal[1]) ;
        //rResult( 3 ) = 0.25 * (1 - rPointLocal[0]) * (1 - rPointLocal[1]) ;

    }

    else if (dimension == 3)
    {
        rResult.resize(4, false);

        double x10 = GetGeometry()[1].Coordinates()[0]-GetGeometry()[0].Coordinates()[0];
        double x20 = GetGeometry()[2].Coordinates()[0]-GetGeometry()[0].Coordinates()[0];
        double x30 = GetGeometry()[3].Coordinates()[0]-GetGeometry()[0].Coordinates()[0];
        double y10 = GetGeometry()[1].Coordinates()[1]-GetGeometry()[0].Coordinates()[1];
        double y20 = GetGeometry()[2].Coordinates()[1]-GetGeometry()[0].Coordinates()[1];
        double y30 = GetGeometry()[3].Coordinates()[1]-GetGeometry()[0].Coordinates()[1];
        double z10 = GetGeometry()[1].Coordinates()[2]-GetGeometry()[0].Coordinates()[2];
        double z20 = GetGeometry()[2].Coordinates()[2]-GetGeometry()[0].Coordinates()[2];
        double z30 = GetGeometry()[3].Coordinates()[2]-GetGeometry()[0].Coordinates()[2];

        rPointLocal[3] = (rPoint[0] - GetGeometry()[0].Coordinates()[0])*(y10*z20 - z10*y20) - (rPoint[1] - GetGeometry()[0].Coordinates()[1])*(x10*z20-x20*z10) + (rPoint[2] - GetGeometry()[0].Coordinates()[2])*(y20*x10 - y10*x20)/mDeterminantJ0;

        rPointLocal[2] = (rPoint[0] - GetGeometry()[0].Coordinates()[0])*(y30*z10-y10*z30) + (rPoint[1] - GetGeometry()[0].Coordinates()[1])*(x10*z30-x30*z10) + (rPoint[2] - GetGeometry()[0].Coordinates()[2])*(y10*x30 - y30*x10)/mDeterminantJ0;

        rPointLocal[1] = (rPoint[0] - GetGeometry()[0].Coordinates()[0])*(y20*z30-y30*z20) + (rPoint[1] - GetGeometry()[0].Coordinates()[1])*(x30*z20-x20*z30) + (rPoint[2] - GetGeometry()[0].Coordinates()[2])*(y30*x20 - x30*y20)/mDeterminantJ0;

        rPointLocal[0] = 1 - rPointLocal[1] - rPointLocal[2] -rPointLocal[3];


        rResult( 0 ) =  1.0-(rPointLocal[0]+rPointLocal[1]+rPointLocal[2]) ;
        rResult( 1 ) = rPointLocal[0] ;
        rResult( 2 ) = rPointLocal[1];
        rResult( 3 ) = rPointLocal[2];


    }

    return rResult;

    KRATOS_CATCH( "" )
}



Vector& UpdatedLagrangianQuadrilateral::MPMLocalCoordinates(Vector& rResult, array_1d<double,3>& rPoint)
{
    KRATOS_TRY

    //Only local coordinated of a point in a tetrahedron is computed
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




Matrix& UpdatedLagrangianQuadrilateral::MPMShapeFunctionsLocalGradients( Matrix& rResult, array_1d<double,3>& rPoint)
{
    unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    array_1d<double,3> rPointLocal = ZeroVector(dimension);

    //local coordinates of the MP position
    //I evaluate the local coordinates of the integration point
    rPointLocal = GetGeometry().PointLocalCoordinates(rPointLocal, rPoint);
    if (dimension == 2)
    {
        rResult = ZeroMatrix( 4, 2 );


        //test (if the first node of the connettivity is the node at the bottom left)
        //rPointLocal[0] = (2 * rPoint[0] - GetGeometry()[2].Coordinates()[0] - GetGeometry()[0].Coordinates()[0])/(GetGeometry()[2].Coordinates()[0] - GetGeometry()[0].Coordinates()[0]);
        //rPointLocal[1] = (2 * rPoint[1] - GetGeometry()[2].Coordinates()[1] - GetGeometry()[0].Coordinates()[1])/(GetGeometry()[2].Coordinates()[1] - GetGeometry()[0].Coordinates()[1]);


        //test (if the first node of the connettivity is the node at the top left)
        //rPointLocal[0] = (2 * rPoint[0] - GetGeometry()[3].Coordinates()[0] - GetGeometry()[1].Coordinates()[0])/(GetGeometry()[3].Coordinates()[0] - GetGeometry()[1].Coordinates()[0]);
        //rPointLocal[1] = (2 * rPoint[1] - GetGeometry()[3].Coordinates()[1] - GetGeometry()[1].Coordinates()[1])/(GetGeometry()[3].Coordinates()[1] - GetGeometry()[1].Coordinates()[1]);

        //test (if the first node of the connettivity is the node at the top right)
        //rPointLocal[0] = (2 * rPoint[0] - GetGeometry()[0].Coordinates()[0] - GetGeometry()[2].Coordinates()[0])/(GetGeometry()[0].Coordinates()[0] - GetGeometry()[2].Coordinates()[0]);
        //rPointLocal[1] = (2 * rPoint[1] - GetGeometry()[0].Coordinates()[1] - GetGeometry()[2].Coordinates()[1])/(GetGeometry()[0].Coordinates()[1] - GetGeometry()[2].Coordinates()[1]);


        //test (if the first node of the connettivity is the node at the bottom right)
        //rPointLocal[0] = (2 * rPoint[0] - GetGeometry()[3].Coordinates()[0] - GetGeometry()[1].Coordinates()[0])/(GetGeometry()[1].Coordinates()[0] - GetGeometry()[3].Coordinates()[0]);
        //rPointLocal[1] = (2 * rPoint[1] - GetGeometry()[3].Coordinates()[1] - GetGeometry()[1].Coordinates()[1])/(GetGeometry()[1].Coordinates()[1] - GetGeometry()[3].Coordinates()[1]);

        //rPointLocal = GetGeometry().PointLocalCoordinates(rPointLocal, rPoint);

        //rResult = GetGeometry().ShapeFunctionsLocalGradients(rResult, rPointLocal);

        //Gradient Shape Function (if the first node of the connettivity is the node at the bottom left)
        rResult( 0, 0 ) = -0.25 * (1 - rPointLocal[1]);
        rResult( 0, 1 ) = -0.25 * (1 - rPointLocal[0]);
        rResult( 1, 0 ) = 0.25 * (1 - rPointLocal[1]);
        rResult( 1, 1 ) = -0.25 * (1 + rPointLocal[0]);
        rResult( 2, 0 ) = 0.25 * (1 + rPointLocal[1]);
        rResult( 2, 1 ) = 0.25 * (1 + rPointLocal[0]);
        rResult( 3, 0 ) = -0.25 * (1 + rPointLocal[1]);
        rResult( 3, 1 ) = 0.25 * (1 - rPointLocal[0]);

        //Gradient Shape Function (if the first node of the connettivity is the node at the top left)
        //rResult( 0, 0 ) = -0.25 * (1 + rPointLocal[1]);
        //rResult( 0, 1 ) = 0.25 * (1 - rPointLocal[0]);
        //rResult( 1, 0 ) = -0.25 * (1 - rPointLocal[1]);
        //rResult( 1, 1 ) = -0.25 * (1 - rPointLocal[0]);
        //rResult( 2, 0 ) = 0.25 * (1 - rPointLocal[1]);
        //rResult( 2, 1 ) = -0.25 * (1 + rPointLocal[0]);
        //rResult( 3, 0 ) = +0.25 * (1 + rPointLocal[1]);
        //rResult( 3, 1 ) = +0.25 * (1 + rPointLocal[0]);

        //Gradient Shape Function (if the first node of the connettivity is the node at the top right)
        //rResult( 0, 0 ) = 0.25 * (1 + rPointLocal[1]);
        //rResult( 0, 1 ) = 0.25 * (1 + rPointLocal[0]);
        //rResult( 1, 0 ) = -0.25 * (1 + rPointLocal[1]);
        //rResult( 1, 1 ) = 0.25 * (1 - rPointLocal[0]);
        //rResult( 2, 0 ) = -0.25 * (1 - rPointLocal[1]);
        //rResult( 2, 1 ) = -0.25 * (1 - rPointLocal[0]);
        //rResult( 3, 0 ) = +0.25 * (1 - rPointLocal[1]);
        //rResult( 3, 1 ) = -0.25 * (1 + rPointLocal[0]);

        //Gradient Shape Function (if the first node of the connettivity is the node at the bottom right)


        //rResult( 0, 0 ) = 0.25 * (1 - rPointLocal[1]);
        //rResult( 0, 1 ) = - 0.25 * (1 + rPointLocal[0]);
        //rResult( 1, 0 ) = 0.25 * (1 + rPointLocal[1]);
        //rResult( 1, 1 ) = 0.25 * (1 + rPointLocal[0]);
        //rResult( 2, 0 ) = - 0.25 * (1 + rPointLocal[1]);
        //rResult( 2, 1 ) = 0.25 * (1 - rPointLocal[0]);
        //rResult( 3, 0 ) = -0.25 * (1 - rPointLocal[1]);
        //rResult( 3, 1 ) = -0.25 * (1 - rPointLocal[0]);



    }
    //if (this->Id() == 541 || this->Id() == 534 || this->Id() == 538)
    //{
    //std::cout<<" rPointLocal "<< " Id "<< this->Id()<<rPointLocal<<std::endl;
    //}
    else if(dimension == 3)
    {
        rResult = ZeroMatrix( 4, 3 );
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

//void UpdatedLagrangianQuadrilateral::CalculateOnIntegrationPoints( const Variable<double>& rVariable, double& rOutput, ProcessInfo& rCurrentProcessInfo )
//{

//rOutput = mConstitutiveLawVector->GetValue( rVariable, rOutput );
//}

//************************************************************************************
//************************************************************************************

//void UpdatedLagrangianQuadrilateral::CalculateOnIntegrationPoints( const Variable<Vector>& rVariable, Vector& rOutput, ProcessInfo& rCurrentProcessInfo )
//{

//KRATOS_TRY





//if ( rVariable == CAUCHY_STRESS_VECTOR || rVariable == PK2_STRESS_VECTOR )
//{
////create and initialize element variables:
//GeneralVariables Variables;
//this->InitializeGeneralVariables(Variables,rCurrentProcessInfo);

            //ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);
            //ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS);

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

//************************************************************************************
//************************************************************************************

//void UpdatedLagrangianQuadrilateral::CalculateOnIntegrationPoints( const Variable<Matrix >& rVariable, Matrix& rOutput, ProcessInfo& rCurrentProcessInfo )
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

//void UpdatedLagrangianQuadrilateral::SetValueOnIntegrationPoints( const Variable<double>& rVariable,
//double& rValues,
//ProcessInfo& rCurrentProcessInfo )
//{
//if (rVariable == DETERMINANT_F){


//mDeterminantF0 = rValues;
//mConstitutiveLawVector->SetValue(rVariable, rValues, rCurrentProcessInfo);


//}
//else{

//UpdatedLagrangianQuadrilateral::SetValueOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
//}
//}
//********************************SET VECTOR VALUE************************************
//************************************************************************************

//void UpdatedLagrangianQuadrilateral::SetValueOnIntegrationPoints( const Variable<Vector>& rVariable, Vector& rValues, ProcessInfo& rCurrentProcessInfo )
//{

//mConstitutiveLawVector->SetValue( rVariable, rValues, rCurrentProcessInfo );


//}


//*******************************SET MATRIX VALUE*************************************
//************************************************************************************

//void UpdatedLagrangianQuadrilateral::SetValueOnIntegrationPoints( const Variable<Matrix>& rVariable, Matrix& rValues, ProcessInfo& rCurrentProcessInfo )
//{

//mConstitutiveLawVector->SetValue( rVariable,
//rValues, rCurrentProcessInfo );


//}
//********************************SET CONSTITUTIVE VALUE******************************
//************************************************************************************

//void UpdatedLagrangianQuadrilateral::SetValueOnIntegrationPoints( const Variable<ConstitutiveLaw::Pointer>& rVariable,
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
//***************************GET DOUBLE VALUE*****************************************
//************************************************************************************

//void UpdatedLagrangianQuadrilateral::GetValueOnIntegrationPoints( const Variable<double>& rVariable,
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

//**************************GET VECTOR VALUE******************************************
//************************************************************************************

//void UpdatedLagrangianQuadrilateral::GetValueOnIntegrationPoints( const Variable<Vector>& rVariable, Vector& rValues, ProcessInfo& rCurrentProcessInfo )
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

//************************************************************************************
//************************************************************************************

//void UpdatedLagrangianQuadrilateral::GetValueOnIntegrationPoints( const Variable<Matrix>& rVariable,
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
//********************************GET CONSTITUTIVE VALUE******************************
//************************************************************************************

//void UpdatedLagrangianQuadrilateral::GetValueOnIntegrationPoints( const Variable<ConstitutiveLaw::Pointer>& rVariable,
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

void UpdatedLagrangianQuadrilateral::GetValuesVector( Vector& values, int Step )
{
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dim = GetGeometry().WorkingSpaceDimension();
    unsigned int MatSize = number_of_nodes * dim;

    if ( values.size() != MatSize ) values.resize( MatSize, false );

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        unsigned int index = i * dim;
        values[index] = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_X, Step );
        values[index + 1] = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_Y, Step );

        if ( dim == 3 )
            values[index + 2] = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_Z, Step );
    }
}


//************************************************************************************
//************************************************************************************

void UpdatedLagrangianQuadrilateral::GetFirstDerivativesVector( Vector& values, int Step )
{
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dim = GetGeometry().WorkingSpaceDimension();
    unsigned int MatSize = number_of_nodes * dim;

    if ( values.size() != MatSize ) values.resize( MatSize, false );

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        unsigned int index = i * dim;
        values[index] = GetGeometry()[i].GetSolutionStepValue( VELOCITY_X, Step );
        values[index + 1] = GetGeometry()[i].GetSolutionStepValue( VELOCITY_Y, Step );

        if ( dim == 3 )
            values[index + 2] = GetGeometry()[i].GetSolutionStepValue( VELOCITY_Z, Step );
    }
}

//************************************************************************************
//************************************************************************************

void UpdatedLagrangianQuadrilateral::GetSecondDerivativesVector( Vector& values, int Step )
{
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dim = GetGeometry().WorkingSpaceDimension();
    unsigned int MatSize = number_of_nodes * dim;

    if ( values.size() != MatSize ) values.resize( MatSize, false );

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        unsigned int index = i * dim;
        values[index] = GetGeometry()[i].GetSolutionStepValue( ACCELERATION_X, Step );
        values[index + 1] = GetGeometry()[i].GetSolutionStepValue( ACCELERATION_Y, Step );

        if ( dim == 3 )
            values[index + 2] = GetGeometry()[i].GetSolutionStepValue( ACCELERATION_Z, Step );
    }
}
//************************************************************************************
//************************************************************************************
void UpdatedLagrangianQuadrilateral::GetHistoricalVariables( GeneralVariables& rVariables )
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

void UpdatedLagrangianQuadrilateral::DecimalCorrection(Vector& rVector)
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
int  UpdatedLagrangianQuadrilateral::Check( const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY
    //std::cout << " AAAAAAAAAAAA "<<std::endl;
    unsigned int dimension = this->GetGeometry().WorkingSpaceDimension();
    //std::cout << " BBBBBBBBBBBB "<<std::endl;
    //verify compatibility with the constitutive law
    ConstitutiveLaw::Features LawFeatures;
    //std::cout << " CCCCCCCCCCCC "<<std::endl;
    this->GetProperties().GetValue( CONSTITUTIVE_LAW )->GetLawFeatures(LawFeatures);
    //std::cout << " DDDDDDDDDDDDD "<<std::endl;

    bool correct_strain_measure = false;
    //std::cout << " EEEEEEEEEEEE "<<std::endl;
    for(unsigned int i=0; i<LawFeatures.mStrainMeasures.size(); i++)
    {
        if(LawFeatures.mStrainMeasures[i] == ConstitutiveLaw::StrainMeasure_Deformation_Gradient)
            correct_strain_measure = true;
    }
    //std::cout << " FFFFFFFFFFFF "<<std::endl;
    if( correct_strain_measure == false )
        KRATOS_THROW_ERROR( std::logic_error, "constitutive law is not compatible with the element type ", " Large Displacements " )

        //std::cout << " GGGGGGGGGGGGGG "<<std::endl;
        //verify that the variables are correctly initialized

        if ( VELOCITY.Key() == 0 )
            KRATOS_THROW_ERROR( std::invalid_argument, "VELOCITY has Key zero! (check if the application is correctly registered", "" )

            if ( DISPLACEMENT.Key() == 0 )
                KRATOS_THROW_ERROR( std::invalid_argument, "DISPLACEMENT has Key zero! (check if the application is correctly registered", "" )

                if ( ACCELERATION.Key() == 0 )
                    KRATOS_THROW_ERROR( std::invalid_argument, "ACCELERATION has Key zero! (check if the application is correctly registered", "" )

                    if ( DENSITY.Key() == 0 )
                        KRATOS_THROW_ERROR( std::invalid_argument, "DENSITY has Key zero! (check if the application is correctly registered", "" )

                        //if ( CAUCHY_STRESS_VECTOR.Key() == 0 )
                        //KRATOS_THROW_ERROR( std::invalid_argument, "CAUCHY_STRESS_VECTOR has Key zero! (check if the application is correctly registered", "" )

                        // if ( BODY_FORCE.Key() == 0 )
                        //     KRATOS_THROW_ERROR( std::invalid_argument, "BODY_FORCE has Key zero! (check if the application is correctly registered", "" );

                        //std::cout << " the variables have been correctly inizialized "<<std::endl;

                        //verify that the dofs exist

                        for ( unsigned int i = 0; i < this->GetGeometry().size(); i++ )
                        {
                            if ( this->GetGeometry()[i].SolutionStepsDataHas( DISPLACEMENT ) == false )
                                KRATOS_THROW_ERROR( std::invalid_argument, "missing variable DISPLACEMENT on node ", this->GetGeometry()[i].Id() )

                                if ( this->GetGeometry()[i].HasDofFor( DISPLACEMENT_X ) == false || this->GetGeometry()[i].HasDofFor( DISPLACEMENT_Y ) == false || this->GetGeometry()[i].HasDofFor( DISPLACEMENT_Z ) == false )

                                    KRATOS_THROW_ERROR( std::invalid_argument, "missing one of the dofs for the variable DISPLACEMENT on node ", GetGeometry()[i].Id() )
                                }

    //verify that the constitutive law exists
    if ( this->GetProperties().Has( CONSTITUTIVE_LAW ) == false )
    {
        KRATOS_THROW_ERROR( std::logic_error, "constitutive law not provided for property ", this->GetProperties().Id() )
    }

    //Verify that the body force is defined
    // if ( this->GetProperties().Has( BODY_FORCE ) == false )
    // {
    //     KRATOS_THROW_ERROR( std::logic_error, "BODY_FORCE not provided for property ", this->GetProperties().Id() )
    // }

    //verify that the constitutive law has the correct dimension
    if ( dimension == 2 )
    {

        if ( THICKNESS.Key() == 0 )
            KRATOS_THROW_ERROR( std::invalid_argument, "THICKNESS has Key zero! (check if the application is correctly registered", "" )

            if ( this->GetProperties().Has( THICKNESS ) == false )
                KRATOS_THROW_ERROR( std::logic_error, "THICKNESS not provided for element ", this->Id() )
            }
    else
    {
        if ( this->GetProperties().GetValue( CONSTITUTIVE_LAW )->GetStrainSize() != 6 )
            KRATOS_THROW_ERROR( std::logic_error, "wrong constitutive law used. This is a 3D element! expected strain size is 6 (el id = ) ", this->Id() )
        }

    //check constitutive law

    if (mConstitutiveLawVector!= 0)
    {
        return mConstitutiveLawVector->Check( GetProperties(), GetGeometry(), rCurrentProcessInfo );
    }

    return 0;

    KRATOS_CATCH( "" );
}




void UpdatedLagrangianQuadrilateral::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Element )
    //int IntMethod = int(mThisIntegrationMethod);
    //rSerializer.save("IntegrationMethod",IntMethod);
    rSerializer.save("ConstitutiveLawVector",mConstitutiveLawVector);
    rSerializer.save("DeformationGradientF0",mDeformationGradientF0);
    rSerializer.save("DeterminantF0",mDeterminantF0);
    rSerializer.save("InverseJ0",mInverseJ0);
    rSerializer.save("DeterminantJ0",mDeterminantJ0);

}

void UpdatedLagrangianQuadrilateral::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Element )
    rSerializer.load("DeformationGradientF0",mDeformationGradientF0);
    rSerializer.load("DeterminantF0",mDeterminantF0);
    rSerializer.load("InverseJ0",mInverseJ0);
    rSerializer.load("DeterminantJ0",mDeterminantJ0);
    //int IntMethod;
    //rSerializer.load("IntegrationMethod",IntMethod);
    //mThisIntegrationMethod = IntegrationMethod(IntMethod);
    rSerializer.load("ConstitutiveLawVector",mConstitutiveLawVector);

}





} // Namespace Kratos

