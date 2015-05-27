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
/* *********************************************************
*
*   Last Modified by:    $Author: bence
*   Date:                $Date: 2011-09-01 17:14:12 $
*   Revision:            $Revision: 1.2 $
*
* ***********************************************************/

// System includes
#include <iostream>
#include<cmath>

// External includes

// Project includes

#include "includes/define.h"
#include "constitutive_laws/viscoplastic_3d.h"

#include "includes/constitutive_law.h"

#include "utilities/math_utils.h"
#include "includes/variables.h"
#include "includes/process_info.h"
#include "freezing_soil_application.h"
#include "includes/properties.h"
#include "structural_application/custom_utilities/sd_math_utils.h"
// #include "sd_math_utils.h"
// #include "custom_utilities/sd_math_utils.h"

namespace Kratos
{


/**
 * TO BE TESTED!!!
 */
Viscoplastic3D::Viscoplastic3D()
        : ConstitutiveLaw()
{
}

/**
 * TO BE TESTED!!!
 */
Viscoplastic3D::~Viscoplastic3D()
{
}



bool Viscoplastic3D::Has( const Variable<double>& rThisVariable )
{
    return false;
}

bool Viscoplastic3D::Has( const Variable<Vector>& rThisVariable )
{
    if ( rThisVariable == INSITU_STRESS )
        return true;

    if ( rThisVariable == INTERNAL_VARIABLES )
        return true;

    return false;
}

bool Viscoplastic3D::Has( const Variable<Matrix>& rThisVariable )
{
    return false;
}



Vector& Viscoplastic3D::GetValue( const Variable<Vector>& rThisVariable, Vector& rValue )
{
    if ( rThisVariable == INSITU_STRESS )
        //      return mInSituStress;

        if ( rThisVariable == INTERNAL_VARIABLES )
        {
            return( mInternalVariables );   // plastic strains + creeptime SIZE = 7
        }


    KRATOS_THROW_ERROR( std::logic_error, "Vector Variable case not considered", "" );
}



Matrix Viscoplastic3D::GetValue( const Variable<Matrix>& rThisVariable )
{
    KRATOS_THROW_ERROR( std::logic_error, "Vector Variable case not considered", "" );
}



void Viscoplastic3D::SetValue( const Variable<double>& rThisVariable, const double& rValue,
                               const ProcessInfo& rCurrentProcessInfo )
{
}


void Viscoplastic3D::SetValue( const Variable<array_1d<double, 3> >& rThisVariable,
                               const array_1d<double, 3>& rValue,
                               const ProcessInfo& rCurrentProcessInfo )
{
}


void Viscoplastic3D::SetValue( const Variable<Vector>& rThisVariable, const Vector& rValue,
                               const ProcessInfo& rCurrentProcessInfo )
{
    if ( rThisVariable == INSITU_STRESS )
    {
        //        mInSituStress = rValue;
    }
    if ( rThisVariable == INTERNAL_VARIABLES )
    {
        mInternalVariables = rValue;
    }
}


void Viscoplastic3D::SetValue( const Variable<Matrix>& rThisVariable, const Matrix& rValue,
                               const ProcessInfo& rCurrentProcessInfo )
{
}



void Viscoplastic3D::InitializeMaterial( const Properties& props, const GeometryType& geom, const Vector& ShapeFunctionsValues )
{
    //KRATOS_WATCH("viscoplastic3d: line 175");
    ResetMaterial( props, geom, ShapeFunctionsValues );

    mInternalVariables = ZeroVector( 7 );
    mMaterialParameters = props[MATERIAL_PARAMETERS];
    //KRATOS_WATCH("viscoplastic3d: line 180");

    //  CohurrentStress = ZeroVector( 6 );
    // Cohtangent = ZeroMatrix( 6, 6 );
    // mInSituStress = ZeroVector( 6 );

}



void Viscoplastic3D::ResetMaterial( const Properties& props,
                                    const GeometryType& geom,
                                    const Vector& ShapeFunctionsValues )
{
//         CalculateElasticMatrix( Cohtangent, mMaterialParameters[0], mMaterialParameters[1] );
}



void Viscoplastic3D::InitializeSolutionStep( const Properties& props,
        const GeometryType& geom, //this is just to give the array of nodes
        const Vector& ShapeFunctionsValues ,
        const ProcessInfo& CurrentProcessInfo )
{

}



void Viscoplastic3D::FinalizeSolutionStep( const Properties& props,
        const GeometryType& geom, //this is just to give the array of nodes
        const Vector& ShapeFunctionsValues ,
        const ProcessInfo& CurrentProcessInfo )
{
    if ( CurrentProcessInfo[CALCULATE_INSITU_STRESS] )
    {
        //mInSituStress -= CohurrentStress;
        //SetValue( INSITU_STRESS, mInSituStress, CurrentProcessInfo );
    }

    for ( unsigned i = 0; i < 6; i++ )
        mInternalVariables[i] = EpsilonPlNew[i];

    mInternalVariables[6] = tv;
}


void  Viscoplastic3D::CalculateMaterialResponse( const Vector& StrainVector,
        const Matrix& DeformationGradient,
        Vector& StressVector,
        Matrix& AlgorithmicTangent,
        const ProcessInfo& CurrentProcessInfo,
        const Properties& props,
        const GeometryType& geom,
        const Vector& ShapeFunctionsValues,
        bool CalculateStresses,
        int CalculateTangent,
        bool SaveInternalVariables )
{

    //example for obtaining temperature on gauss point
    // using ShapeFunctionsValues and nodal values of the geometry

    // KRATOS_WATCH("viscoplastic3d: line 245");

    temp = 0.0;
    for ( unsigned int i = 0; i < geom.size(); i++ )
    {
        temp += geom[i].GetSolutionStepValue( TEMPERATURE ) * ShapeFunctionsValues[i];
    }

    //retrieving time step
    double time_step = CurrentProcessInfo[DELTA_TIME];

/*  
      std::cout<<"*************************************  INPUTS ***************************************** "<<std::endl;
      std::cout<<"mInternalVariables= "<<mInternalVariables<<", C33 = "<< AlgorithmicTangent(2,2)<<std::endl;
      std::cout<<"StrainVector "<<StrainVector<<std::endl;
      std::cout<<"***************************************************************************** "<<std::endl;
  */    
    

    //=================CONSTITUTIVE CALCULATION=========================
    GetStress( StrainVector, mInternalVariables, time_step, temp, StressVector, AlgorithmicTangent );

    /*
     std::cout<<"****************************************  OUTPUTS  *************************************** "<<std::endl;
     std::cout<<"PlStrains= "<<EpsilonPlNew<<", C33 = "<<AlgorithmicTangent(2,2)<<std::endl;
     std::cout<<", Stress = "<<StressVector<<std::endl;
     std::cout<<"******************************************************************************* "<<std::endl;
    */
}



/* =========================================
 Viscoplastic Returnmapping for Frozen Soils
 ===========================================

 Material Constants: Emod, nu, FR, Coh, beta, CreepC, A1, A2, K1

 Convergence constants: mTOL1, mTOL2, maxITER

 Input data: Old plastic strains, New Total Strains, local time-variable for creep, timestep size, temperature

 Output data: EpsilonPlNew, SigmaNew, Cvp and updated tv

 */


void Viscoplastic3D::GetStress( const Vector& StrainVector,
                                Vector& IntVar,     //Plastic Strains and creep time //Will be updated and returned
                                double deltaT,
                                double temp,     //Matrix& StrainGradient, // Matrix& SigmaOldMatrix
                                Vector& SigmaNew,    //Will be returned
                                Matrix& AlgorithmicTangent )    //Will be returned
{

    //Initialise and change to vector notation
    //========================================

    double Emod = mMaterialParameters[0];		//Elastic constants
    double nu = mMaterialParameters[1];
    double pMelt1 = mMaterialParameters[2];		//Yield function definition
    double pMelt2 = mMaterialParameters[3];
    double qmax = mMaterialParameters[4];
    double beta = mMaterialParameters[5];		//Creep parameters for Orth
    double K1 = mMaterialParameters[6];
    double CreepC = mMaterialParameters[7];
    double A1 = mMaterialParameters[8];
    double A2 = mMaterialParameters[9];
    double temp2 = mMaterialParameters[10];		// Temerature input (for elements without temperature DOF
    double selectMethod = mMaterialParameters[11];	//0: Pure plastic   1: eta = tau*E  2: Cudmani  3:  Dev.Split method
    double b1 = mMaterialParameters[12];		
    
    
    temp= temp2;	// KINEMATIC ELEMENT HAS NO TEMPERATURE SO IT WILL BE READ AS A MATERIAL PARAMETER

    if (temp > 0.0)
    {
        KRATOS_WATCH("ERROR: POSITIVE TEMPERATURE");
    }


    double mTOL1 = 0.001;
    double mTOL2 = 0.01;
    int maxITER = 40;


    Vector EpsilonPlOld( 6 );
    noalias( EpsilonPlOld ) = ZeroVector( 6 );

    Vector EpsilonNew( 6 );
    noalias( EpsilonNew ) = ZeroVector( 6 );

    Matrix EpsilonNewMatrix(6,6);
    noalias(EpsilonNewMatrix) = ZeroMatrix(6,6);

    
    ///TEST NO Plastic strains remaing //	
   
    for ( unsigned int i = 0; i < 6; i++ )
        EpsilonPlOld( i ) = IntVar( i );

    tv = IntVar( 6 );


    EpsilonNew = StrainVector;


    // Compute Elasticity Matrix

    Matrix Cmat( 6, 6 );
    noalias( Cmat ) = ZeroMatrix( 6, 6 );

    //KRATOS_WATCH("viscoplastic3d: line 363");

    double lambda = nu * Emod / (( 1.0 + nu ) * ( 1.0 - 2.0 * nu) );
    double mu = Emod / ( 2.0 * ( 1.0 + nu ) );

    //KRATOS_WATCH("viscoplastic3d: line 368");

    Cmat( 0, 0 ) = 2.0 * mu + lambda;
    Cmat( 0, 1 ) = lambda;
    Cmat( 0, 2 ) = lambda;
    Cmat( 1, 0 ) = lambda;
    Cmat( 1, 1 ) = 2.0 * mu + lambda;
    Cmat( 1, 2 ) = lambda;
    Cmat( 2, 0 ) = lambda;
    Cmat( 2, 1 ) = lambda;
    Cmat( 2, 2 ) = 2.0 * mu + lambda;
    Cmat( 3, 3 ) = mu;
    Cmat( 4, 4 ) = mu;
    Cmat( 5, 5 ) = mu;


    // Define inverse C

    Matrix Cinv(6, 6);

    noalias( Cinv ) = ZeroMatrix( 6, 6 );
    Cinv( 0, 0 ) = 1.0 / Emod;
    Cinv( 0, 1 ) = -nu / Emod;
    Cinv( 0, 2 ) = -nu / Emod;
    Cinv( 1, 0 ) = -nu / Emod;
    Cinv( 1, 1 ) = 1.0 / Emod;
    Cinv( 1, 2 ) = -nu / Emod;
    Cinv( 2, 0 ) = -nu / Emod;
    Cinv( 2, 1 ) = -nu / Emod;
    Cinv( 2, 2 ) = 1.0 / Emod;
    Cinv( 3, 3 ) = 1.0 / mu;
    Cinv( 4, 4 ) = 1.0 / mu;
    Cinv( 5, 5 ) = 1.0 / mu;



    // Compute deviatoric matrix

    Matrix Idev( 6, 6 );

    //KRATOS_WATCH("viscoplastic3d: line 403");

    noalias( Idev ) = ZeroMatrix( 6, 6 );

    Idev( 0, 0 ) = 2.0 / 3.0;
    Idev( 0, 1 ) = -1.0 / 3.0;
    Idev( 0, 2 ) = -1.0 / 3.0;

    Idev( 1, 0 ) = -1.0 / 3.0;
    Idev( 1, 1 ) = 2.0 / 3.0;
    Idev( 1, 2 ) = -1.0 / 3.0;

    Idev( 2, 0 ) = -1.0 / 3.0;
    Idev( 2, 1 ) = -1.0 / 3.0;
    Idev( 2, 2 ) = 2.0 / 3.0;

    Idev( 3, 3 ) = 1.0;
    Idev( 4, 4 ) = 1.0;
    Idev( 5, 5 ) = 1.0;
    
    
    

    // Compute values for the yield function
    //======================================

    /// TEST to exclude them
    
    double pMelt = pMelt1;		// If direct yieldcurve definition
    
    /*
    double yb = 0.04 * exp(4.69 / 16.52 * sqrt(-temp) );
    double yc = 5.37 * 1000000.0 * exp(9.61 / 273 * (-temp) );

    double pMelt = -temp / 0.074 * 1000000.0/meltRatio;
    double pMelt2 = -1.0 / ( ( 1.0/pMelt ) + ( yb/yc) );

    double qmax = yc + yb / 4.0 * pMelt * pMelt / ( yc/yb + pMelt);
    */
    
    
    double pm2 = (pMelt - pMelt2) * (pMelt - pMelt2);



    // Compute trial stress and check flow
    //====================================

    Vector SigmaTr( 6 );
    noalias( SigmaTr ) = ZeroVector( 6 );		//Initialise trial stress vector

    Vector DevTr( 6 );
    noalias( DevTr ) = ZeroVector( 6 );			//Initialise trial deviatoric vector

    for ( unsigned int i = 0; i < 6; i++ )		// Compute trial stress
        for ( unsigned int j = 0; j < 6; j++ )
            SigmaTr( i ) += Cmat( i, j ) * (EpsilonNew( j ) - EpsilonPlOld( j ) );  

    for ( unsigned int i = 0; i < 6; i++ )		//Execute the deviatoric split of sigma trial
        for ( unsigned int j = 0; j < 6; j++ )
            DevTr( i ) += Idev( i, j ) * SigmaTr( j );   

    double meanP = ( SigmaTr( 0 ) + SigmaTr( 1 ) + SigmaTr( 2 ) ) / 3.0;  // Compute the hydrostatic pressure

    double normDev = 0.0;     			       //  Initialise the norm of the deviatoric tensor

    for ( unsigned int i = 0; i < 3; i++ )
        normDev += DevTr( i ) * DevTr( i ) + 2.0 * DevTr( i + 3 ) * DevTr( i + 3 );	// Compute the norm of S

    normDev = sqrt(normDev );			// Finish the computation of the norm

    double qq = sqrt(3.0/2.0)*normDev;		// Get q for the yield curve

    double pp = -meanP;				// Get p for the yield curve
    
    double p1 = pp;				// Save for later use
    
    double q1 = qq;

    Vector Dev1( 6 );
    Dev1 = DevTr;

    double normDev1 = normDev;
    

    double yield = pow( 2.0 * pp - pMelt - pMelt2 ,2.0) / pm2 + qq * qq / qmax / qmax - 1.0;	// Compute the isosurface at the trial stress state
    


    /////TEST : no yield at all
    //yield = -2.0;
    
    
    
    /*std::cout<<"  Sigma=    "<<SigmaTr<<std::endl;
    std::cout<<" 			"<<std::endl;
    std::cout<<"  S=    "<<DevTr<<std::endl;
    std::cout<<" 			"<<std::endl;
    std::cout<<"  p=    "<<pp<<std::endl;
    std::cout<<" 			"<<std::endl;
    std::cout<<"  q=    "<<qq<<std::endl;
    std::cout<<" 			"<<std::endl;*/
    


    if ( yield <= mTOL1 )  //Elastic step, no return mapping and viscous regularisation
    {
        SigmaNew = SigmaTr;
        EpsilonPlNew = EpsilonPlOld;
        AlgorithmicTangent = Cmat;

        // std::cout<<"ELASTIC STEP."<<std::endl;

        return;
    }

    //KRATOS_WATCH("viscoplastic3d: line 460");

    // Start returnmapping
    // ===================

    //Reset and initialise variables

    double deltaGamma = 0.0;

    double delta2Gamma = 0.0;

    Vector DeltaEpsilonPlastic( 6 );
    noalias( DeltaEpsilonPlastic ) = ZeroVector( 6 );

    Vector Resid( 6 );
    noalias( Resid ) = ZeroVector( 6 );

    Vector df1( 6 );
    noalias( df1 ) = ZeroVector( 6 );

    Matrix df2( 6, 6 );
    noalias( df2 ) = ZeroMatrix( 6, 6 );

    Vector Help1( 6 );
    noalias( Help1 ) = ZeroVector( 6 );

    Help1( 0 ) = 1.0 / 3.0;
    Help1( 1 ) = 1.0 / 3.0;
    Help1( 2 ) = 1.0 / 3.0;

    Matrix Ainv( 6, 6 );
    noalias( Ainv ) = ZeroMatrix( 6, 6 );

    Matrix Amat( 6, 6 );
    noalias( Amat ) = ZeroMatrix( 6, 6 );

    Vector dummy( 6 );
    noalias( dummy ) = ZeroVector( 6 );

    Vector dummy2( 6 );
    noalias( dummy2 ) = ZeroVector( 6 );

    EpsilonPlNew = EpsilonPlOld;
    
    bool conv = false;


    for ( unsigned int iter = 0; iter < maxITER; iter++ )
    {

        //Compute SigmaTr and the Yield function
        //======================================


        SigmaTr = ZeroVector( 6 );

        for ( unsigned int i = 0; i < 6; i++ )
            for ( unsigned int j = 0; j < 6; j++ )
                SigmaTr( i ) += Cmat( i, j ) * ( EpsilonNew( j ) - EpsilonPlNew( j ) );

        DevTr = ZeroVector( 6 );

        for ( unsigned int i = 0; i < 6; i++ )
            for ( unsigned int j = 0; j < 6; j++ )
                DevTr( i ) += Idev( i, j ) * SigmaTr( j );

        meanP = ( SigmaTr( 0 ) + SigmaTr( 1 ) + SigmaTr( 2 ) ) / 3.0;

        normDev = 0.0;

        for ( unsigned int i = 0; i < 3; i++ )
            normDev += DevTr( i ) * DevTr( i ) + 2.0 * DevTr( i + 3 ) * DevTr( i + 3 );

        normDev = sqrt( normDev );

        qq = sqrt(3.0/2.0)*normDev;

        pp = -meanP;


        yield = pow( 2.0 * pp - pMelt - pMelt2 , 2.0) / pm2 + qq * qq / qmax / qmax - 1.0;



        //Compute the derivatives of f
        //============================

        df1 = ZeroVector(6);

        for (unsigned int i = 0; i < 3; i++)
        {
            df1(i) = -4.0/3.0 * ( 2.0 * pp - pMelt -pMelt2) / pm2 + 3.0 * DevTr(i) /qmax/qmax;
            df1(i+3) = 6.0 * DevTr(i+3) /qmax/qmax;
        }

        double derA = 8.0/9.0/pm2 + 2.0/qmax/qmax;
        double derB = 8.0/9.0/pm2 - 1.0/qmax/qmax;
        double derC = 6.0/qmax/qmax;

        df2 = ZeroMatrix(6,6);

        for (unsigned int i = 0; i < 3; i++)
        {
            for (unsigned int j = 0; j < 3; j++)
            {
                if (i == j)
                {
                    df2(i,j) = derA;
                    df2(i+3,j+3) = derC;
                }
                else df2(i,j) = derB;
            }
        }



        //Check convergence
        //=================

        Resid = -EpsilonPlNew + EpsilonPlOld + deltaGamma * df1;

        double normResid = 0.0;

        for ( unsigned int i = 0; i < 6; i++ )
            normResid += Resid( i ) * Resid( i );

        normResid = sqrt( normResid );

        if ( (yield <= mTOL1 and normResid <= mTOL2) or (selectMethod == 3.0) )
        {
	    //KRATOS_WATCH("RM Converged: "<<iter);
	    conv = true;
            SigmaNew = SigmaTr;
            EpsilonPlNew = EpsilonPlNew;
            break;
        }


        //Execute the computation
        //=======================

        Ainv = Cinv + deltaGamma * df2;

        SD_MathUtils<double>::InvertMatrix( Ainv, Amat );

        dummy = ZeroVector( 6 );

        for ( unsigned int i = 0; i < 6; i++ )
            for ( unsigned int j = 0; j < 6; j++ )
                dummy( i ) += df1( j ) * Amat( j, i );  // df^T* A

        double dfar = 0.0;

        for ( unsigned int i = 0; i < 6; i++ )
            dfar += dummy( i ) * Resid( i );  // df^T* A* R

        dfaf = 0.0;

        for ( unsigned int i = 0; i < 6; i++ )
            dfaf += dummy( i ) * df1( i );   // df^T* A* df

        delta2Gamma = ( yield - dfar ) / dfaf;

        dummy2 = ZeroVector( 6 );

        for ( unsigned int i = 0; i < 6; i++ )
            for ( unsigned int j = 0; j < 6; j++ )
                dummy2( i ) += Amat( i, j ) * ( Resid( j ) + delta2Gamma * df1( j ) );

        DeltaEpsilonPlastic = ZeroVector( 6 );

        for ( unsigned int i = 0; i < 6; i++ )
            for ( unsigned int j = 0; j < 6; j++ )
                DeltaEpsilonPlastic( i ) += Cinv( i, j ) * dummy2( j );


        //Update the variables
        //====================

        EpsilonPlNew += DeltaEpsilonPlastic;

        deltaGamma += delta2Gamma;

        /*
           std::cout<<"  Iter=    "<<iter<<std::endl;
           std::cout<<"  J2=    "<<normDev<<std::endl;
           std::cout<<"  p =    "<<meanP<<std::endl;
           std::cout<<"   "<<std::endl;
           std::cout<<"  df1=    "<<df1<<std::endl;
           std::cout<<" 							"<<std::endl;
           std::cout<<"  df2=    "<<df2<<std::endl;
           std::cout<<" 							"<<std::endl;
           std::cout<<"  Sigmatr=    "<<SigmaTr<<std::endl;
           std::cout<<" 							"<<std::endl;
           std::cout<<"  yield=    "<<yield<<std::endl;
           std::cout<<" 							"<<std::endl;
           std::cout<<"  Normresid=    "<< normResid<<std::endl;
           std::cout<<" 							"<<std::endl;
           std::cout<<"  dfarrd=    "<< dfar<<std::endl;
           std::cout<<" 							"<<std::endl;
           std::cout<<"  dfaf=    "<< dfaf<<std::endl;
           std::cout<<" 							"<<std::endl;
           std::cout<<"  Normresid=    "<< normResid<<std::endl;
           std::cout<<" 							"<<std::endl;
           std::cout<<"  Normresid=    "<< normResid<<std::endl;
           std::cout<<" 							"<<std::endl;
           std::cout<<"  delta2GammaC=    "<<delta2Gamma<<std::endl;
           std::cout<<" 							"<<std::endl;
           std::cout<<"  deltaEpsilonPL=    "<<DeltaEpsilonPlastic<<std::endl;
           std::cout<<" 							"<<std::endl;
           std::cout<<"  EpsilonPLNew=    "<<EpsilonPlNew<<std::endl;
           std::cout<<" 							"<<std::endl;
        */

        if (iter == maxITER - 1) 
	{
	  KRATOS_WATCH("**DIV**");
	}

    }   /////Return mapping loop's END




    //===========================================
    //Calculate Elastoplastic algorithmic tangent
    //===========================================



    Vector dummy3( 6 );
    noalias( dummy3 ) = ZeroVector( 6 );

    Matrix Dummy4( 6, 6 );
    noalias( Dummy4 ) = ZeroMatrix( 6, 6 );

    Matrix Cep( 6, 6 );
    noalias( Cep ) = ZeroMatrix( 6, 6 );

    for ( unsigned int i = 0; i < 6; i++ )
        for ( unsigned int j = 0; j < 6; j++ )
            dummy3( i ) += Amat( i, j ) * df1( j ); // This is N1

    for ( unsigned int i = 0; i < 6; i++ )
        for ( unsigned int j = 0; j < 6; j++ )
            Dummy4( i, j ) += dummy3( i ) * dummy3( j );   // This is N1 o N1

    Cep = Amat - ( 1.0 / dfaf ) * Dummy4; // Elastoplastic algorithmic tangent




    //========================
    //Compute viscous response
    //========================

    //Initialise
    //==========



    if (selectMethod == 0.0 or !conv)
    {
        SigmaNew = SigmaTr;
        AlgorithmicTangent = Cep;
        return;
    }


    Vector DeltaEpsVp( 6 );
    noalias( DeltaEpsVp ) = ZeroVector( 6 );

    Vector SigmaInf( 6 );
    noalias( SigmaInf ) = ZeroVector( 6 );

    SigmaInf = SigmaNew;       // converged stress from return mapping's SigmaTrial
    
    
    SigmaTr = ZeroVector( 6 );

    for ( unsigned int i = 0; i < 6; i++ )
        for ( unsigned int j = 0; j < 6; j++ )
            SigmaTr( i ) += Cmat( i, j ) * ( EpsilonNew( j ) - EpsilonPlOld( j ) );  // Overstress
	    
    ///
	    
   


    //Compute the creep rate (Orth)
    //=============================

    double sigmaA = A1 * pow( -temp, A2 );

    double cOfT = K1 / ( temp + 273.4 );

    double dotEpsA = 1.0;		// dimension: 1/min

    //test for fixed q1
    //q1 = 8000000.0;


    double dotEpsM = 0.01 * dotEpsA * exp( cOfT * ( q1 / sigmaA - 1.0 ));//   + pow(q1 /3 /p1 ,2.0) - 1.0 );  // from SigmaTR
    
    double dotEpsM2 = 0.01 * dotEpsA * exp( cOfT * ( qq / sigmaA - 1.0 ));//  + pow(qq /3 /pp ,2.0) - 1.0 );  // from SigmaTR
        
    double tc = 0.0;

    if (tv == 0.0) tc = tv + 0.5 * deltaT; 
    else tc = tv;

    //test
    //double tm = CreepC / dotEpsM;
    double tm = CreepC / dotEpsM2;
    
    double dotEps = dotEpsM * exp( -beta ) * exp( beta * tc / tm ) * pow( tc / tm , -beta );
    
    double dotEps2 = dotEpsM2 * exp( -beta ) * exp( beta * tc / tm ) * pow( tc / tm , -beta );
    
    
   
    

//Method 2 		// Only deviatoric stresses induce viscosity
//========


    if ( selectMethod == 2.0 )
    {
        //test
	double eta = abs(q1) / dotEps2;
    
        double tau = eta / mu;	
	
	if ( tau < 0.00001) SigmaNew = SigmaInf;
	else SigmaNew = ( SigmaTr + deltaT / tau * SigmaInf ) / ( 1.0 + deltaT / tau );
	
	Vector SigmaNewDev(6);
	SigmaNewDev = ZeroVector(6);

	//TEST for not only deviatoric evolution

	for ( unsigned int i = 0; i < 6; i++ )
            for ( unsigned int j = 0; j < 6; j++ )
                SigmaNewDev(i) += Idev(i,j) * SigmaNew(j);
	
	
	EpsilonPlNew = ZeroVector(6);
        for ( unsigned int i = 0; i < 6; i++ )
            for ( unsigned int j = 0; j < 6; j++ )
                EpsilonPlNew(i) += Cinv(i,j) * SigmaNew(j);
    
        EpsilonPlNew = EpsilonNew - EpsilonPlNew;
	
	
	
	if ( tau < 0.00001) AlgorithmicTangent = Cmat;
        else AlgorithmicTangent = ( Cmat + deltaT / tau * Cep ) / ( 1.0 + deltaT / tau );
	
	
	
        tv += deltaT;
	
	/*
	if (  conv )
	{
	
	std::cout<<"  qINF    "<<qq<<std::endl;
	std::cout<<" 							"<<std::endl;
	std::cout<<"  SigmaTr3=    "<<SigmaTr(2)<<std::endl;
	std::cout<<" 							"<<std::endl;
	std::cout<<"  SigmaInf=    "<<SigmaInf(2)<<std::endl;
	std::cout<<" 							"<<std::endl;
	std::cout<<"  SigmaNew3=    "<<SigmaNew(2)<<std::endl;
	std::cout<<" 							"<<std::endl;
	std::cout<<"  C33=    "<<Cmat(2,2)<<std::endl;
	std::cout<<"  Cep33=    "<<Cep(2,2)<<std::endl;
	std::cout<<"  Cvp33=    "<<AlgorithmicTangent(2,2)<<std::endl;
	std::cout<<" 							"<<std::endl;
	std::cout<<"  dotepsM=    "<<dotEpsM<<std::endl;
	std::cout<<" 							"<<std::endl;
	//std::cout<<"  dotepsM2=    "<<dotEpsM2<<std::endl;
	std::cout<<" 							"<<std::endl;
	std::cout<<"  doteps=    "<<dotEps<<std::endl;
	std::cout<<" 							"<<std::endl;
	//std::cout<<"  doteps2   "<< dotEps2<<std::endl;
	
	std::cout<<"  q1=    "<<q1<<std::endl;
	std::cout<<" 							"<<std::endl;
	std::cout<<"  dotepsM   "<< dotEpsM<<std::endl;
	std::cout<<" 							"<<std::endl;
	std::cout<<"  tau=    "<<tau<<std::endl;
	std::cout<<"  tm=    "<<tm<<std::endl;
	}
	*/
        return;
    }
    
    
  
//Method 1
//========    // Simple viscoplastic regularisation with tau = eta / E

    if ( selectMethod == 1.0 )
    {


        double eta = abs(meanP) / dotEps;
        double tau = eta / Emod;
	

	///TEST fixed tau ==>. secondary creep
	//tau = 5.0;

        SigmaNew = ( SigmaTr + deltaT / tau * SigmaInf ) / ( 1.0 + deltaT / tau );


	EpsilonPlNew = ZeroVector(6);
        for ( unsigned int i = 0; i < 6; i++ )
            for ( unsigned int j = 0; j < 6; j++ )
                EpsilonPlNew(i) += Cinv(i,j) * SigmaNew(j);

        EpsilonPlNew = EpsilonNew - EpsilonPlNew;
	
        AlgorithmicTangent = ( Cmat + deltaT / tau * Cep ) / ( 1.0 + deltaT / tau );
	
        tv += deltaT;
	
	       
        return;
    }
    
    
    
    
    
    
//Method 3		// Baased on Cudmani, no time integration
//========    

    if (selectMethod == 3.0)
    {
      //double eta = abs(q1) / dotEps;
	
      //double tau = eta / mu;
      
      DeltaEpsilonPlastic = Dev1 / normDev1 * dotEps;
      
      EpsilonPlNew = EpsilonPlOld + DeltaEpsilonPlastic * deltaT;

      
      SigmaNew = ZeroVector(6);
        for ( unsigned int i = 0; i < 6; i++ )
            for ( unsigned int j = 0; j < 6; j++ )
                SigmaNew(i) += Cmat(i,j) * ( EpsilonNew(j) - EpsilonPlNew(j) );

      AlgorithmicTangent =  Cmat; // + deltaT / tau * Cep ) / ( 1.0 + deltaT / tau );
      
      tv += deltaT;

    std::cout<<"  q1=    "<<q1<<std::endl;
	std::cout<<" 							"<<std::endl;
	std::cout<<"  dotepsM   "<< dotEpsM<<std::endl;
	std::cout<<" 							"<<std::endl;
	std::cout<<"  tm=    "<<tm<<std::endl;
   	std::cout<<" 							"<<std::endl;
	std::cout<<"  DeVP=    "<<DeltaEpsilonPlastic(2)<<std::endl;
   	std::cout<<" 							"<<std::endl;
	std::cout<<"  Sigmanew=    "<<SigmaNew(2)<<std::endl;
   	std::cout<<" 							"<<std::endl;

    }






/// TODO to be tested
//Method 4			
//========

// split into deviatoric and volumetric stresses and handle separately with tauE and tauG based on Method 1

    if ( selectMethod == 4.0 )
    {
	double Kmod = (lambda + 2.0/3.0 * mu);
	
        double tauE = 1.0 / dotEps / 3.0 /Kmod;
	
        double tauG = 1.0 / dotEps / mu;

	
        //Split of SigmaInf

        Vector SigmaVolInf( 6 );
        noalias( SigmaVolInf ) = ZeroVector( 6 );
	
	Vector SigmaDevInf( 6 );
	SigmaDevInf = ZeroVector(6);
	
	meanP = ( SigmaInf(0) + SigmaInf(1) + SigmaInf(2) )/ 3.0;

        SigmaVolInf = 3.0 * meanP * Help1;

	for ( unsigned int i = 0; i < 6; i++ )
            for ( unsigned int j = 0; j < 6; j++ )
                SigmaDevInf( i ) += Idev( i, j ) * SigmaInf( j );
	    
	
        //Split of SigmaTr

        Vector SigmaVolTr( 6 );
        noalias( SigmaVolTr ) = ZeroVector( 6 );

        Vector SigmaDevTr( 6 );
        noalias( SigmaDevTr ) = ZeroVector( 6 );

	meanP = ( SigmaTr( 0 ) + SigmaTr( 1 ) + SigmaTr( 2 ) ) / 3.0;

        SigmaVolTr = 3.0 * meanP * Help1;
	
	for ( unsigned int i = 0; i < 6; i++ )
            for ( unsigned int j = 0; j < 6; j++ )
                SigmaDevTr( i ) += Idev( i, j ) * SigmaTr( j );

	 //sum volumetric stress and deviatoric stress 

	SigmaNew = ( SigmaVolTr + deltaT / tauE * SigmaVolInf ) / ( 1.0 + deltaT / tauE ) +  ( SigmaDevTr + deltaT / tauG * SigmaDevInf ) / ( 1.0 + deltaT / tauG ); 

		 
	//Split of DeltaEpsilonVP
	/* Vector DeltaEpsDev( 6 );
	for ( unsigned int i = 0; i < 6; i++ )
	    for ( unsigned int j = 0; j < 6; j++ )
	    {
		DeltaEpsVol( i ) += Cinv( i, j ) * ( SigmaVolTr( j ) - SigmaVolInf( j ) ) / tauE;
		DeltaEpsDev( i ) += Cinv( i, j ) * ( SigmaDevTr( j ) - SigmaDevInf( j ) ) / tauG;
	    }
	DeltaEpsVp = DeltaEpsVol + DeltaEpsDev;*/
       
       
	
	EpsilonPlNew = ZeroVector(6);
        for ( unsigned int i = 0; i < 6; i++ )
            for ( unsigned int j = 0; j < 6; j++ )
                EpsilonPlNew(i) += Cinv(i,j) * SigmaNew(j);

        EpsilonPlNew = EpsilonNew - EpsilonPlNew;

	double tauMean = 2.0 * dotEps / (1.5 * Kmod + mu );
	
        for ( unsigned int i = 0; i < 6; i++ )
	  for ( unsigned int j = 0; j < 6; j++ )
                AlgorithmicTangent(i,j) = ( Cmat(i,j) + deltaT / tauMean * Cep(i,j) ) / ( 1.0 + deltaT / tauMean );
          
        
	  
	  
	 /* 
	std::cout<<"  taue=    "<< tauE<<std::endl;
	std::cout<<" 							"<<std::endl;
	std::cout<<"  taug=    "<< tauG<<std::endl;
	std::cout<<" 							"<<std::endl;
	std::cout<<"  sigma inf=    "<< SigmaTr<<std::endl;
	std::cout<<" 							"<<std::endl;
	std::cout<<"  sigma tr=    "<< SigmaTr<<std::endl;
	std::cout<<" 							"<<std::endl;
	std::cout<<"  sigma volinf=    "<< SigmaVolInf<<std::endl;
	std::cout<<" 							"<<std::endl;
	std::cout<<"  sigma devinf=    "<< SigmaDevInf<<std::endl;
	std::cout<<" 							"<<std::endl;

	  
        
        std::cout<<"  sigma new=    "<< SigmaNew<<std::endl;
	std::cout<<" 							"<<std::endl;
	std::cout<<"  epsilonnewVP=    "<< EpsilonPlNew<<std::endl;
	std::cout<<" 							"<<std::endl;
	std::cout<<"  c11=    "<< AlgorithmicTangent(0,0)<<std::endl;
	std::cout<<" 							"<<std::endl;
	std::cout<<"  c44=    "<<AlgorithmicTangent(3,3)<<std::endl;
	std::cout<<" 							"<<std::endl;
	
	
	std::cout<<"  SigmaTR=    "<<SigmaTr<<std::endl;
	std::cout<<" 							"<<std::endl;
	std::cout<<"  SigmaInf=    "<<SigmaInf<<std::endl;
	std::cout<<" 							"<<std::endl;
	std::cout<<"  SigmaNew=    "<<SigmaNew<<std::endl;
	std::cout<<" 							"<<std::endl;
	std::cout<<"  EpsPlOld=    "<<EpsilonPlOld<<std::endl;
	std::cout<<" 							"<<std::endl;
	std::cout<<"  EpsVpnew   "<< EpsilonPlNew<<std::endl;
	std::cout<<" 							"<<std::endl;

	std::cout<<"  C33=    "<< AlgorithmicTangent(2,2)<<std::endl;
	std::cout<<" 							"<<std::endl;
    *
	std::cout<<"  qq=    "<<qq<<std::endl;
	std::cout<<" 							"<<std::endl;
	std::cout<<"  pp=    "<<pp<<std::endl;
	std::cout<<" 							"<<std::endl;
	std::cout<<"  dotepsM=    "<<dotEpsM<<std::endl;
	std::cout<<" 							"<<std::endl;
	std::cout<<"  doteps   "<< dotEps<<std::endl;
	std::cout<<" 							"<<std::endl;
	std::cout<<"  tc=    "<< tc<<std::endl;
	std::cout<<" 							"<<std::endl;
	std::cout<<"  tm=    "<< tm<<std::endl;
	std::cout<<" 							"<<std::endl;
	std::cout<<"  eta=    "<<eta<<std::endl;
	std::cout<<" 							"<<std::endl;
	std::cout<<"  tau=    "<<tau<<std::endl;
	std::cout<<" 							"<<std::endl;
	  */
	
	
        tv += deltaT;

        return;

    }

}




//**********************************************************************
int Viscoplastic3D::Check( const Properties& props, const GeometryType& geom, const ProcessInfo& CurrentProcessInfo )
{
    KRATOS_TRY
    double nu = props[POISSON_RATIO];

    if ( nu > 0.499 && nu < 0.501 )
    {
        KRATOS_THROW_ERROR( std::logic_error, "invalid poisson ratio in input, close to incompressibility", "" );
        return -1;
    }

    return 0;

    KRATOS_CATCH( "" );
}
} // Namespace Kratos

