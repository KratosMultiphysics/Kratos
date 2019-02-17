
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
*   Last Modified by:    $Author: virginia $
*   Date:                $Date: 2010-02-05 16:10:12 $
*   Revision:            $Revision: 1.4 $
*
* ***********************************************************/

/// BELOW ///////// TANGENT MATRIX

// System includes
#include <iostream>
#include <limits>
// External includes
#include<cmath>

// Project includes

#include "includes/define.h"
#include "constitutive_laws/viscofibers_2d.h" // USAR 2D com "D" maiusculo?

#include "includes/constitutive_law.h"  // nelson
#include "custom_utilities/tensor_utils.h"  // nelson
#include "custom_utilities/sd_math_utils.h"  // nelson

#include "utilities/math_utils.h"
#include "includes/variables.h"
#include "includes/process_info.h"
#include "structural_application.h"
#include "includes/properties.h"

namespace Kratos
{
/// VISCOFIBERS2D
/**Defines a visco-hyperelastic constitutive law with 2 classes of fibers in 2D space,
for quasi-incompressible materials.
//////////////////
*/

namespace Viscofibers2DAuxiliaries
{
}
using namespace Viscofibers2DAuxiliaries;


/**
 *	TO BE TESTED!!!
 */
Viscofibers2D::Viscofibers2D()
    : ConstitutiveLaw()
{
}
/**
 *	TO BE TESTED!!!
 */
Viscofibers2D::~Viscofibers2D()
{
}


bool Viscofibers2D::Has( const Variable<double>& rThisVariable )
{
    return false;
}

bool Viscofibers2D::Has( const Variable<Vector>& rThisVariable )
{
    return false;
}

bool Viscofibers2D::Has( const Variable<Matrix>& rThisVariable )
{
    return false;
}

void Viscofibers2D::SetValue( const Variable<double>& rThisVariable, const double& rValue,
                              const ProcessInfo& rCurrentProcessInfo )
{
}

void Viscofibers2D::SetValue( const Variable<Vector>& rThisVariable, const Vector& rValue,
                              const ProcessInfo& rCurrentProcessInfo )
{
}

void Viscofibers2D::SetValue( const Variable<Matrix>& rThisVariable, const Matrix& rValue,
                              const ProcessInfo& rCurrentProcessInfo )
{
}

void Viscofibers2D::Calculate(const Variable<Matrix >& rVariable, Matrix& rResult,
                              const ProcessInfo& rCurrentProcessInfo)
{
    for (unsigned int ii = 0; ii<3; ii++)
        rResult(0,ii) = mhn1(ii); // What for ??????????????
    KRATOS_WATCH("rResult");
    KRATOS_WATCH(rResult);
    //rResult(0,ii) = mplastic_strain(ii);
}

// NEW GET VALUE !!!!!!!!!
double& Viscofibers2D::GetValue( const Variable<double>& rThisVariable, double& rValue )
{
    KRATOS_THROW_ERROR(std::logic_error, "Vector Variable case not considered" , "");
}
Vector& Viscofibers2D::GetValue( const Variable<Vector>& rThisVariable, Vector& rValue )
{
    if( rThisVariable == INSITU_STRESS )
        return mInSituStress;
    if( rThisVariable == INTERNAL_VARIABLES )
    {
        rValue = ZeroVector(1);
        return( rValue );
    }
    KRATOS_THROW_ERROR(std::logic_error, "Vector Variable case not considered", "");
}

Matrix& Viscofibers2D::GetValue( const Variable<Matrix>& rThisVariable, Matrix& rValue )
{
    KRATOS_THROW_ERROR(std::logic_error,"Vector Variable case not considered", "");
}
// NEW GET VALUE !!!!!!!!!


/**
 *	TO BE TESTED!!!
 */
void Viscofibers2D::InitializeMaterial( const Properties& props,
                                        const GeometryType& geom,
                                        const Vector& ShapeFunctionsValues )
{
    //mCurrentStress = ZeroVector(6);

    mMU = props[MU]; // shear modulus
    mK = props[BULK_MODULUS]; // bulk modulus
    //mE = props[YOUNG_MODULUS];
    //mNU = props[POISSON_RATIO];
    mThickness = props[THICKNESS]; // to check
    mA0= geom.Area();  //geom[A0]; // to check
    mA= geom.Area();  // JUST A TEST, need to check if it will be updated later!!!
    mAlpha=props[ALPHA];  // angle used for retraction
    mRetractionTime=props[RETRACTION_TIME];  // Time to start retraction
    // nelson
    mMaterialParameters.resize((props[MATERIAL_PARAMETERS]).size(), false);
    noalias(mMaterialParameters) = props[MATERIAL_PARAMETERS];
    //KRATOS_WATCH(props[MATERIAL_PARAMETERS])
    //KRATOS_WATCH(mMaterialParameters)
    noalias(md) = ZeroVector(16);
    md[6]=mMaterialParameters[0];
    md[7]=mMaterialParameters[1];
    md[8]=mMaterialParameters[2];
    md[9]=mMaterialParameters[3];
    md[10]=mMaterialParameters[4];
    md[11]=mMaterialParameters[5];
    md[12]=mMaterialParameters[6];
    md[13]=mMaterialParameters[7];
    md[14]=mMaterialParameters[8];
    md[15]=mMaterialParameters[9];
    //noalias(md) = mMaterialParameters;
    //KRATOS_WATCH(md)
    // nelson
    // VISCO FIBERS PARAMETERS K1 and K2
    //md[3]= props[YOUNG_MODULUS];// K1
    //mK1=md[3];
    //md[4]= props[POISSON_RATIO];// K2
    //mK2=md[4];
    mK1=mMaterialParameters[10];
    mK2=mMaterialParameters[11];
    mrvec1_00=mMaterialParameters[12]; // x direction
    mrvec2_01=mMaterialParameters[16]; // y direction
    // 1st fiber
    Vector mrvec1(3);
    mrvec1[0]=mMaterialParameters[12];
    mrvec1[1]=mMaterialParameters[13];
    mrvec1[2]=mMaterialParameters[14];
    //KRATOS_WATCH(mrvec1)
    // 2nd fiber
    Vector mrvec2(3);
    mrvec2[0]=mMaterialParameters[15];
    mrvec2[1]=mMaterialParameters[16];
    mrvec2[2]=mMaterialParameters[17];
    // END VISCO FIBERS PARAMETERS K1 and K2
    //noalias(md) = ZeroVector(16);
    //KRATOS_WATCH("HELLOOOOOOOOOOOOOOOOOOOOOOOOOOOOO");
    // viscous parameters  HOW DO I DECLARE THEIR VALUES????????
    // beta - for each Maxwell device
    //md[6]=props[YOUNG_MODULUS]; // just to test
    //md[6]=0.25; // b1 Christian uses=0.354;
// 		md[7]=0.286; // b2
// 		md[8]=0.298; // b3
// 		md[9]=0.285; // b4
// 		md[10]=0.348; // b5
    // tau (relaxation time) - for each Maxwell device
    //md[11]=props[POISSON_RATIO];// just to test
    //md[11]=0.0001; // ta1 Christian uses=0.001;
// 		md[12]=0.010; // ta2
// 		md[13]=0.010; // ta3
// 		md[14]=1.000; // ta4
// 		md[15]=10.00; // ta5
    /*mCtangent.resize(6,6,false);
    noalias(mCtangent) = ZeroMatrix(6,6);*/
    //mMaterialParameters = props[MATERIAL_PARAMETERS];
    mInSituStress = ZeroVector(6);
    //KRATOS_WATCH("HELLOOOOOOOOOOOOOOOOOOOOOOOOOOOOO 222222222");
    noalias(mhn1)=ZeroVector(36);  // initialize history vector
    noalias(mcurrent_hn1)=ZeroVector(36);  // initialize history vector
    //first six components are S_elastic at time n, following 30 are S_viscous at time n

    //CalculateElasticMatrix(mCtangent, props[YOUNG_MODULUS], props[POISSON_RATIO]);
//                 CalculateElasticMatrix(mCtangent, props[MATERIAL_PARAMETERS][0], props[MATERIAL_PARAMETERS][1]);
}

void Viscofibers2D::InitializeSolutionStep( const Properties& props,
        const GeometryType& geom, //this is just to give the array of nodes
        const Vector& ShapeFunctionsValues ,
        const ProcessInfo& CurrentProcessInfo)
{
}

void Viscofibers2D::FinalizeSolutionStep( const Properties& props,
        const GeometryType& geom, //this is just to give the array of nodes
        const Vector& ShapeFunctionsValues ,
        const ProcessInfo& CurrentProcessInfo)
{
    if( CurrentProcessInfo[CALCULATE_INSITU_STRESS] )
    {
        //noalias(mhn1) =  mcurrent_hn1;  // viscoelasticity !!!!!! original
        //mInSituStress -= mCurrentStress;
        //SetValue( INSITU_STRESS, mInSituStress, CurrentProcessInfo );
    }
}

void Viscofibers2D::ResetMaterial( const Properties& props,
                                   const GeometryType& geom,
                                   const Vector& ShapeFunctionsValues)
{

}

void Viscofibers2D::CalculateMaterialResponse( const Vector& StrainVector,
        const Matrix& DeformationGradient,
        Vector& StressVector,
        Matrix& AlgorithmicTangent,
        const ProcessInfo& CurrentProcessInfo,
        const Properties& props,
        const GeometryType& geom,
        const Vector& ShapeFunctionsValues,
        bool CalculateStresses,
        int CalculateTangent,
        bool SaveInternalVariables)
{
    mA=geom.Area();//props[AREA]; // to check
    //KRATOS_WATCH("HELLOOOOOOOOOOOOOOOOOOOOOOOOOOOOO VISCOUS!!!");
    mcurrentThickness = (mA0*mThickness)/mA; // to check
    //KRATOS_WATCH(mcurrentThickness);
    // RETRACTION ///////////////
    if (CurrentProcessInfo[TIME] >  mRetractionTime)
        //if (CurrentProcessInfo[TIME] > 0.1)   //0.01
    {
        mRetraction=mAlpha*CurrentProcessInfo[TIME]; // [DELTA_TIME]
        //KRATOS_WATCH("retraction");
    }
    else
    {
        mRetraction=0.0;
    }

    // VISCOELASTICITY ///////////////
    UpdateMaterial(StrainVector, props, geom,ShapeFunctionsValues, CurrentProcessInfo); // nelson
    dt=CurrentProcessInfo[DELTA_TIME]; // [DELTA_TIME]
    //KRATOS_WATCH(dt);
    //KRATOS_WATCH(DeformationGradient);
    //CalculateComposite(DeformationGradient,rTau); // NEW VISCO_FIBER
    //CalculateNewComposite(StrainVector, DeformationGradient, rTau); // NEW VISCO_FIBER
    //CalculateElasticMatrix(mCtangent, props[YOUNG_MODULUS], props[POISSON_RATIO]);
    CalculateStress(StrainVector, StressVector);  // this function to pass tau too!!!
    //CalculateStress(StrainVector, StressVector);  // I need to change this function to pass tau too!!!
    //CalculateConstitutiveMatrix(StrainVector, Ctang);  // nelson
    CalculateConstitutiveMatrix(StrainVector, Ctang);  // TEST
    //CalculateConstitutiveMatrix(StrainVector, rResult);  // nelson
// 		mA=geom.Area();//props[AREA]; // to check
//
// 		mcurrentThickness = (mA0*mThickness)/mA; // to check
// 		// RETRACTION ///////////////
// 		if (CurrentProcessInfo[TIME] > 0.0000001)   //0.01
// 		{
// 		mRetraction=mAlpha*CurrentProcessInfo[TIME]; // [DELTA_TIME]
// 		//KRATOS_WATCH("retraction");
// 		}
// 		else
// 		{
// 		mRetraction=0.0;
// 	        }
// 		//KRATOS_WATCH(CurrentProcessInfo[TIME]);
// 		//KRATOS_WATCH(mRetraction);
// 		// RETRACTION ///////////////
// 		//KRATOS_WATCH(mA0);
// 		//KRATOS_WATCH(mA);
// 		//KRATOS_WATCH(mThickness);
// 		//KRATOS_WATCH(mcurrentThickness);

}

/**
 *	TO BE TESTED!!!
 */

// Implementing the Strain Energy for the Hyperelastic material (matrix)

// this function gives W (strain energy), S (2nd Piola-Kirch) and CC (elasticity tensor) or Ctang (tangent matrix)
// Ctang is also known as the Material Tangent Constitutive Tensor

// ORDER:
//#1 void Isotropic2D::CalculateElasticMatrix(Matrix& C, const double E, const double NU)
//#2	void Isotropic2D::CalculateStress(const Vector& StrainVector, Vector& StressVector)
//#3	void Isotropic2D::CalculateConstitutiveMatrix(const Vector& StrainVector, Matrix& rResult)
//#4	void Isotropic2D::CalculateCauchyStresses(Vector& rCauchy_StressVector,
//		const Matrix& rF,
//		const Vector& rPK2_StressVector,
//		const Vector& rGreenLagrangeStrainVector)

//void Viscofibers2D::CalculateConstitutiveMatrix(const Vector& StrainVector, Matrix& rResult)
//*************************************************************************************************************************
void Viscofibers2D::CalculateStress(const Vector& StrainVector, Vector& StressVector)
//void Viscofibers2D::CalculateStress(const Vector& StrainVector, Vector& StressVector)
{
    // this funcion CalculateStress is called by EBST
    if( StressVector.size() != 3 ) // EBST
    {
        StressVector.resize(3,false); // false para no hacer copias
    }
    //noalias(StressVector) = prod(mCtangent,StrainVector);
    //mCurrentStress = StressVector;
    //noalias(StressVector) -= mInSituStress;


    ////TRIAL
    //Viscofibers2D::CalculateStressVector(StrainVector, auxStressVector);
    CalculateStressVector(StrainVector,auxStressVector);
    ////TRIAL
    StressVector[0] = auxStressVector[0];
    StressVector[1] = auxStressVector[1];
    StressVector[2] = auxStressVector[3];
    // new trial to add rTau
// 		StressVector[0] = auxStressVector[0]+rTau[0];
// 		StressVector[1] = auxStressVector[1]+rTau[1];
// 		StressVector[2] = auxStressVector[3]+rTau[2];
    //KRATOS_WATCH(StressVector)

}
// PS NOALIAS SOLO EN LADO ESQUIERDO!
//*************************************************************************************************************************
/**
Calculate Second Piola
*/

void Viscofibers2D::CalculateStressVector(const Vector& StrainVector, Vector& auxStressVector)
//void Viscofibers2D::CalculateStressVector(const Vector& StrainVector, array_1d<double,6>& rResult)
{
    // NEW // EBST
    if( auxStressVector.size() != 6 )
    {
        auxStressVector.resize(6,false); // false para no hacer copias
    }
    // NEW //
    //KRATOS_WATCH("HYPER !!!!!!!!!");
    Matrix I(3,3); // creating I matrix (Identity)
    //I=new IdentityMatrix(3);
    noalias(I)=ZeroMatrix(3,3);
    I(0,0)=1.0;
    I(1,1)=1.0;
    I(2,2)=1.0;

    Matrix C(3,3); // creating matrix C (C=Ft.F)
    noalias(C)=ZeroMatrix(3,3);

    // for EBST and Membrane
    if (mcurrentThickness == 0.0)  // DO I NEED THIS?!?!?!
    {
        mcurrentThickness=(mA0/mA)*mThickness;
    }
    //KRATOS_WATCH(mA0);
    //KRATOS_WATCH(mA);
    //KRATOS_WATCH(mThickness);
    //KRATOS_WATCH(mcurrentThickness);
    double Ez = 0.0;
    // Uncomment below if you want to use Ez
    Ez = (1.0/2.0)*((pow(mcurrentThickness,2.0)-pow(mThickness,2.0))/pow(mThickness,2.0));
    //KRATOS_WATCH(Ez);
    // TO INCLUDE RETRACTION !!!!!!!!!!
    Vector auxStrainVector(6); // 6 components
    auxStrainVector[0] = StrainVector[0]-mRetraction*(0.8*StrainVector[0]);
    auxStrainVector[1] = StrainVector[1]-mRetraction*(0.8*StrainVector[1]);
    auxStrainVector[2] = Ez;
    auxStrainVector[3] = StrainVector[2];
    auxStrainVector[4] = 0.0;
    auxStrainVector[5] = 0.0;
    //KRATOS_WATCH(auxStrainVector);
    //Vector auxStressVector(6); // 6 components

    //// IN TOTAL LAGRANGIAN! (BELOW)
    C(0,0) = 2.0 * auxStrainVector[0] + 1.0;
    C(1,1) = 2.0 * auxStrainVector[1] + 1.0;
    C(2,2) = 2.0 * auxStrainVector[2] + 1.0;
    C(0,1)=auxStrainVector[3];
    C(1,2)=auxStrainVector[4];
    C(0,2)=auxStrainVector[5];
    C(1,0)=C(0,1);
    C(2,1)=C(1,2);
    C(2,0)=C(0,2);
    //KRATOS_WATCH(C)

    Matrix Cinv(3,3);
    noalias(Cinv)=ZeroMatrix(3,3);
    double det=0.0; // det of C  ( det(C)=J2 )
    //MathUtils<double>::InvertMatrix3(Celastic,Cinv,det);
    MathUtils<double>::InvertMatrix3(C,Cinv,det); // Cinv and Det(C) are computed

    double J2=MathUtils<double>::Det3( C ); // determinant of C is computed
    //KRATOS_WATCH(J2);
//    double J=0.0;
//    J=sqrt(J2); // J is computed
    //KRATOS_WATCH(J);


    // computing the derivative of strain energy: dW = S (the deviatoric part)

//    double W=0.0;
    // Vector dW(3);
    Matrix dW(3,3);
    noalias(dW)=ZeroMatrix(3,3);
    double trC=0.0;
    trC=C(0,0)+C(1,1)+C(2,2); // computing Ckk (Trace of C = Sum of the diagonal terms) = I1 (First Invariant)
    //W=0.5*mMU*((pow(J2,0.3333333333333333333333333333))*trC-3.0);
//    W=0.5*mMU*((1.0/(pow(J2,0.3333333333333333333333333333)))*trC-3.0); // ENERGY

    Matrix mdwtemp(3,3);
    noalias(mdwtemp)=ZeroMatrix(3,3);
    noalias(mdwtemp)=Cinv;
    mdwtemp*=trC;
    mdwtemp*=(-1.0/3.0);
    mdwtemp += I;	// parentesis solved
    mdwtemp *= mMU;
    mdwtemp *= 1.0/pow(J2,0.3333333333333333333333333333);
    noalias(dW) = mdwtemp; // S computed!!! :)
    // Second Piola Kirchhoff is 3x3

    // puting second Piola Kirchhodd in Vector form
    Vector SPKir(6); // 6 components
    SD_MathUtils<double>::TensorToVector(dW,SPKir); // Rossi's Vector form
    //KRATOS_WATCH("HYPERRRRRRRRRRRRRR !!!!!!!!!");
    //***************************
    /**
    Calculate the viscoelastic contribution for Second Piola stresses
    */
    // VISCOELASTICITY !!!!!!!!

    // beta makes the correlation between elastic and maxwell bodies
    // tau is the scalar relaxation time
    //Vector d(16);  // need to be defined (6 first components=???, following 5 components for beta, and 5 for tau)
    // aux parameter tid

    // FIBERS !!!!!!!!!
    CalculateComposite(StrainVector, rTau); // NEW VISCO_FIBER
    //CalculateComposite(StrainVector,tau); // it's called in Calculate Material Response !

    double tid = -dt/2.0; //
    int i=0;
    Vector qsi(5);
    noalias(qsi)=ZeroVector(5); // new
    Vector sqr_qsi(5);
    Vector Beta(5);
    noalias(Beta)=ZeroVector(5); // new
    //Vector beta_v(5);
    Vector Se_n(6);
    Vector Sv_n(30);
    Vector Sv(6);
    //KRATOS_WATCH("VISCOUS !!!!!!!!!");
    for  (i=0; i<5; i++)
    {
        //while (md(11+i)!=0.0) // new
        //{
        Beta(i)=md(6+i); // fills beta vector

        // qsi is a dimention-less parameter (usually depends on the time step)
        qsi(i)=exp(tid/md(11+i)); // qsi = xi (Christian's code)
        //KRATOS_WATCH (qsi(i));
        //qsi(i)=md[i]; // qsi = xi (Christian's code)
        sqr_qsi(i)=qsi(i)*qsi(i); // qsi^2

        //}
    }
    //KRATOS_WATCH (md);
    //KRATOS_WATCH (Beta);
    //KRATOS_WATCH (qsi);
    //KRATOS_WATCH (sqr_qsi);
    // 		  pick elastic stress at time n
    // 		  hn1 is a history vector of 36 components for each gauss point
    // 		  the first 6 components are the elastic stresses at time n
    // 		  the other 30 components are the viscous stresses at time n for each Maxwell device
    // 		  (6 comp x 5 devices = 30)
    //KRATOS_WATCH(mhn1);
    for  (i=0; i<6; i++)  // to go through the 6 components
    {
        Se_n(i)=mhn1(i); // elastic stress (HYP) at time n - got from History Parameter

    }
    //KRATOS_WATCH("elastic stress at time n");
    //KRATOS_WATCH(Se_n);
    // update history, store elastic stress at time n

    for  (i=0; i<6; i++)  // to go through the 6 components
    {
        //hn1(i)=tau(i); // elastic stress (HYP) of time n+1
        //mcurrent_hn1(i)=SPKir(i); // elastic stress (HYP) of time n+1 // original
        mhn1(i)=SPKir(i);
    }
    //KRATOS_WATCH("elastic stress at time n+1 !!");
    //KRATOS_WATCH(mhn1);
    // pick overstress at time tn

    int j=0;
    int alp=0;
    // loop over the alpha number of Maxwell elements
    for  (alp=0; alp<5; alp++) // Considering 5 series of Maxwell elements
    {

        for  (i=0; i<6; i++)  // to go through the 6 components
        {
            // computing viscoelastic over stress at time n
            // Sv=viscoelastic contribution of Second Piola stress
            Sv_n(i+j)=mhn1(i+j+6); // to update from the 7th component
        }

        j=j+6;
        //alp=alp+1;
    }
    //KRATOS_WATCH("Viscoelastic stress time n - check 7th component on");
    //KRATOS_WATCH(Sv_n);
    j=0;

    Vector tau2(6); // initialize tau2
    noalias(tau2) = ZeroVector(6);
    //double tau2 = 0.0 // initializen tau2

    noalias(Sv)= ZeroVector(6); // change it to vector!!
    // computing viscoelastic overstress for a particular device
    // loop over the alpha number of Maxwell elements
    for  (alp=0; alp<5; alp++) // Considering 5 series of Maxwell elements
    {
        //noalias(Sv)= ZeroVector(6); // change it to vector!!

        for  (i=0; i<6; i++)  // to go through the 6 components
        {
            // Sv=viscoelastic contribution of Second Piola stress
            //KRATOS_WATCH(sqr_qsi(alp));
            //KRATOS_WATCH(Sv_n(i+j));
            Sv(i)=sqr_qsi(alp)*Sv_n(i+j)+ qsi(alp)*Beta(alp)*(SPKir(i)-Se_n(i)) ; // to update from the 7th component
            //KRATOS_WATCH(Sv(i));
            // TEST
            mhn1(i+j+6)=Sv(i);
            tau2(i)=tau2(i)+Sv(i);
        }
        j=j+6;
        //KRATOS_WATCH(Sv);
        // update history: store overstress at time n
// 			    for  (i=0; i<6; i++)  // to go through the 6 components
// 			      {
// 			      //mcurrent_hn1(i+j+6)=Sv(i); // original
// 			      mhn1(i+j+6)=Sv(i);
// 			      // adding viscous part to the total stress
// 			      tau2(i)=tau2(i)+Sv(i); // tau2 is an aux variable to store the viscous stress
//
// 			      }
//
// 			     j=j+6;
    }
    //KRATOS_WATCH("FINAL Sv");
    //KRATOS_WATCH(Sv);
    //KRATOS_WATCH(mhn1);
// 			 KRATOS_WATCH("TAU OVERSTRESS");
// 			 KRATOS_WATCH(tau2);
// 			 KRATOS_WATCH("SPKir OLD");
// 			 KRATOS_WATCH(SPKir);
    for  (i=0; i<6; i++) // Considering 5 series of Maxwell elements
    {
        // tau is my total stress (StressVector)
        SPKir(i)=tau2(i)+SPKir(i); // here I'm adding the viscous contribution to the total stress
    }
// 			 KRATOS_WATCH("SPKir NEW");
// 			 KRATOS_WATCH(SPKir);
    // END VISCOELASTICITY
    //***************************
    // Volumetric part of the Stresses

    Matrix Svol(3,3);
    noalias(Svol)=ZeroMatrix(3,3);
    double beta=9.0; // previous value: Beta=9
    double auxSvol=0.0;
    auxSvol=(mK/beta)*(1.0-(1.0/pow(J2,4.5))); // previous one!
    //KRATOS_WATCH(Cinv);
    noalias(Svol)=Cinv;
    Svol*=auxSvol;
    //KRATOS_WATCH(Svol);
    Vector vSvol(6); // to put Svol in Vector form
    SD_MathUtils<double>::TensorToVector(Svol,vSvol); // Rossi's Vector form
    // Comment if not considering volumetric part
    //KRATOS_WATCH(SPKir);
    //KRATOS_WATCH(vSvol);
    noalias(auxStressVector)=SPKir+vSvol;
    //KRATOS_WATCH(auxStressVector);
    // NEW !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // Fibers contribution
    auxStressVector[0] = auxStressVector[0]+rTau[0];
    auxStressVector[1] = auxStressVector[1]+rTau[1];
    auxStressVector[3] = auxStressVector[3]+rTau[2];
    // NEW !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //StressVector[0] = auxStressVector[0];
    //StressVector[1] = auxStressVector[1];
    //StressVector[2] = auxStressVector[3]; // xy component
    //KRATOS_WATCH(StressVector);

}
//*******************************************************************************************
/**
Calculate the contribution of two families of fibers
*/
/*
	//void Viscofibers2D::CalculateComposite(const Matrix& rF, Vector& tau)
	void Viscofibers2D::CalculateComposite(const Matrix& DeformationGradient, Vector& rTau)
		{
		  	if( rTau.size() != 3 ) // EBST
			{
			rTau.resize(3,false);
			}

			Vector tau(3); // initialize tau2
			noalias(tau) = ZeroVector(3);

			//double J = MathUtils<double>::Det2( DeformationGradient );
			//double J = MathUtils<double>::Det2( rF );

			// vectors of fibers orientation : anl and ane
			Vector anl(3); // in lagrangian description
			noalias(anl)=ZeroVector(3); // new
			Vector ane(3); // in eulerian description
			noalias(ane)=ZeroVector(3); // new
			int i=0;
			int j=0;
			int ii=0;

			// Just testing some random vectors )perpendicular fiber directions
			//KRATOS_WATCH(mrvec1)
			Matrix rvec(2,3);
// 			rvec(0,0)=mrvec1[0];
// 			rvec(0,1)=mrvec1[1];
// 			rvec(0,2)=mrvec1[2];
// 			KRATOS_WATCH("rvec1 ok !")
// 			rvec(1,0)=mrvec2[0];
// 			rvec(1,1)=mrvec2[1];
// 			rvec(1,2)=mrvec2[2];
			rvec(0,0)=mrvec1_00; // the norm of this vector needs to be 1.0! so I don introduce stretch
			rvec(0,1)=0.0;
			rvec(0,2)=0.0;
			rvec(1,0)=0.0;
			rvec(1,1)=mrvec2_01;  // the norm of this vector needs to be 1.0!
			rvec(1,2)=0.0;

			// trial
			Vector norm(2);
			for  (i=0; i<2; i++)
			  {
			  norm(i)=sqrt(rvec(i,0)*rvec(i,0)+rvec(i,1)*rvec(i,1)+rvec(i,2)*rvec(i,2));
			  }
			//double norm1=sqrt(rvec(0,0)*rvec(0,0)+rvec(0,1)*rvec(0,1)+rvec(0,2)*rvec(0,2));
			KRATOS_WATCH(norm);
			//double norm2=rvec(1,0)+rvec(1,1)+rvec(1,2);
			// trial

			for  (ii=0; ii<2; ii++)
			{
			Vector tau_f(6); // initialize the stress contribution of each fiber
			noalias(tau_f) = ZeroVector(6);



			      for  (j=0; j<3; j++) //
				{
				    // tau is my total stress (StressVector)
				    //anl(j)=rvec(ii,j,n);
				    anl(j)=rvec(ii,j);
				    //anl(j)=rvec(ii,j)*(1.0/norm(ii)); // "ii" is  each fiber direction, "j" is each component, "n" is each element? I think I don't need this n!
				}
				double J = MathUtils<double>::Det2( DeformationGradient );
			        //KRATOS_WATCH (J);
				//KRATOS_WATCH (DeformationGradient);
				double J_onethird = pow(J,(-1.0/3.0));
				//KRATOS_WATCH (J_onethird);
				double J_23 = pow(J,(-2.0/3.0));
				//KRATOS_WATCH (J_23);
				//I don't need below! (the modified deformation gradient)
				boost::numeric::ublas::bounded_matrix<double,2,2> Fbar; // Fbar = J_onethird* F, modified deformation gradient
				//Fbar=J_onethird*DeformationGradient;
				noalias(Fbar)=DeformationGradient;
				//KRATOS_WATCH (Fbar);
				Fbar*=J_onethird;
				//KRATOS_WATCH (Fbar);
				noalias(ane)=ZeroVector(3); // new
				for  (i=0; i<3; i++) //
				{
				  for  (j=0; j<3; j++) //
				  {
				      // tau is my total stress (StressVector)
				      ane(i)=ane(i)+(Fbar(i,j)*anl(j)); // "ii" is  each fiber direction, "j" is each component, "n" is each element? I think I don't need this n!
				  }
				}
				//KRATOS_WATCH (ane);
				//I don't need upper part! (the eulerian form of the fiber)

				// calls a function that computes tau_F (the second piola kirrchoff stress contribution of the fiber)
				//CalculateFiber(ane,tau_f); used for Kirchhoff
				//CalculateFiber(anl, tau_f);
				// TEST
				Vector a(3);
				//noalias(a)=anl; // copy anl in ae ! Is it correct ???? only to compute w4
				noalias(a)=ane; // copy ane in ae !

				double I4 =0.0;
				double I4_1 =0.0;
				//int j=0;
				int i=0;
				for  (i=0; i<3; i++) //
				      {
					// tau is my total stress (StressVector)
					I4=I4+a(i)*a(i); // "ii" is  each fiber direction, "j" is each component, "n" is each element? I think I don't need this n!
				      }
				if ( I4 < 1.0 )
				{
				  //tau=0.0; // which tau is this one?
				  tau_f=ZeroVector(3);
				}
				else
				{

				  // CALL PSI I4 !! 3rd function : derivative of W in respect to I4
				  I4_1 = I4-1.0; // I4 -1
				  KRATOS_WATCH (mK1);
				  double temp= exp(mK2*I4_1*I4_1);
				  // first derivative of the strain energy function in respect to the 4th invariant
				  double w4 = mK1*I4_1*temp;
				  // second derivative of the stored energy function in respect to the 4th invariant
				  //double w44 =  mK1*(1.0+(2.0*mK2*I4_1*I4_1))*temp; // It's used only for the tangent. So I won't use it!

				  // CALL vect_product ae !! 4th function : deviator of "a"
				  Vector a_vs_a(6); // cross product of ae vector, Christian calls it "axa"
				  noalias(a)=anl; // new!
				  int i=0;
				  for  (i=0; i<3; i++)
				      {
					// tau is my total stress (StressVector)
					a_vs_a(i)=a(i)*a(i);
				      }
				  a_vs_a[3]=a[0]*a[1];
				  a_vs_a[4]=a[1]*a[2];
				  a_vs_a[5]=a[0]*a[2];

				  // CALL dev_euler !! 5th function
				  Vector dev_a_vs_a(6); // check the size of this vector !!
				  noalias(dev_a_vs_a)=a_vs_a;
				  i=0;
				  double one3 = 1.0/3.0;
				  double tra= a_vs_a(0)+a_vs_a(1)+a_vs_a(2);
				  tra=one3*tra;
				  for  (i=0; i<3; i++) //
				      {
					// tau is my total stress (StressVector)
					dev_a_vs_a(i)=dev_a_vs_a(i)-tra;
				      }
				  Matrix C(2,2); // creating matrix C (C=Ft.F)
				  noalias(C)=ZeroMatrix(2,2);
				  noalias(C)=prod(trans(DeformationGradient),DeformationGradient);
				  Matrix Cinv(3,3);
				  noalias(Cinv)=ZeroMatrix(3,3);
				  double det=0.0; // det of C  ( det(C)=J2 )
				  MathUtils<double>::InvertMatrix2(C,Cinv,det); // Cinv and Det(C) are computed
				  Vector vCinv(3); // 6 components
				  SD_MathUtils<double>::TensorToVector(Cinv,vCinv);
				  //tau_f=2.0*w4*dev_ae_vs_ae; This is the Kirchoff stress, I need Second Piola Kirchhoff
				  noalias(tau_f) = dev_a_vs_a;
				  noalias(tau_f) -= one3*I4*vCinv;
				  tau_f *= 2.0*w4*J_23;
				  //tau_f=2.0*w4*J_23*(dev_a_vs_a-(one3*I4*Cinv));

				// END TEST


				tau=tau+tau_f;
				}
			noalias(rTau)=tau;

			}

		}*/
//*******************************************************************************************
/**
Calculate the contribution of two families of fibers
*/

//void Viscofibers2D::CalculateComposite(const Matrix& rF, Vector& tau)
void Viscofibers2D::CalculateComposite(const Vector& StrainVector,  Vector& rTau)
{
    if( rTau.size() != 3 ) // EBST
    {
        rTau.resize(3,false);
    }

    Vector tau(3); // initialize tau2
    noalias(tau) = ZeroVector(3);
    Vector tau_f(6); // initialize the stress contribution of each fiber
    noalias(tau_f) = ZeroVector(6);
    //double J = MathUtils<double>::Det2( DeformationGradient );
    //double J = MathUtils<double>::Det2( rF );

    // vectors of fibers orientation : anl and ane
    Vector anl(3); // in lagrangian description
    noalias(anl)=ZeroVector(3); // new
    Vector ane(3); // in eulerian description
    noalias(ane)=ZeroVector(3); // new
    int i=0;
    int j=0;
    int ii=0;

    // Just testing some random vectors )perpendicular fiber directions
    //KRATOS_WATCH(mrvec1)
    Matrix rvec(2,3);
// 			rvec(0,0)=mrvec1[0];
// 			rvec(0,1)=mrvec1[1];
// 			rvec(0,2)=mrvec1[2];
// 			KRATOS_WATCH("rvec1 ok !")
// 			rvec(1,0)=mrvec2[0];
// 			rvec(1,1)=mrvec2[1];
// 			rvec(1,2)=mrvec2[2];
    rvec(0,0)=mrvec1_00; // the norm of this vector needs to be 1.0! so I don introduce stretch
    rvec(0,1)=0.0;
    rvec(0,2)=0.0;
    rvec(1,0)=0.0;
    rvec(1,1)=mrvec2_01;  // the norm of this vector needs to be 1.0!
    rvec(1,2)=0.0;

    // trial
    Vector norm(2);
    for  (i=0; i<2; i++)
    {
        norm(i)=sqrt(rvec(i,0)*rvec(i,0)+rvec(i,1)*rvec(i,1)+rvec(i,2)*rvec(i,2));
    }
    //double norm1=sqrt(rvec(0,0)*rvec(0,0)+rvec(0,1)*rvec(0,1)+rvec(0,2)*rvec(0,2));
    //KRATOS_WATCH(norm);
    //double norm2=rvec(1,0)+rvec(1,1)+rvec(1,2);
    // trial
    //*************************************************
    // filling anl
    //int ii=0;
    for  (ii=0; ii<2; ii++)
    {
        for  (j=0; j<3; j++) //
        {
            //anl(j)=rvec(ii,j);
            anl(j)=rvec(ii,j)*(1.0/norm(ii)); // "ii" is  each fiber direction, "j" is each component, "n" is each element? I think I don't need this n!
        }
        //}
        //*************************************************
        // Below: computing C and J
        Matrix I(3,3); // creating I matrix (Identity)
        //I=new IdentityMatrix(3);
        noalias(I)=ZeroMatrix(3,3);
        I(0,0)=1.0;
        I(1,1)=1.0;
        I(2,2)=1.0;

        Matrix C(3,3); // creating matrix C (C=Ft.F)
        noalias(C)=ZeroMatrix(3,3);

        // for EBST and Membrane
        if (mcurrentThickness == 0.0)  // DO I NEED THIS?!?!?!
        {
            mcurrentThickness=(mA0/mA)*mThickness;
        }
        //KRATOS_WATCH(mA0);
        //KRATOS_WATCH(mA);
        //KRATOS_WATCH(mThickness);
        //KRATOS_WATCH(mcurrentThickness);
        double Ez = 0.0;
        // Uncomment below if you want to use Ez
        Ez = (1.0/2.0)*((pow(mcurrentThickness,2.0)-pow(mThickness,2.0))/pow(mThickness,2.0));
        //KRATOS_WATCH(Ez);
        // TO INCLUDE RETRACTION !!!!!!!!!!
        Vector auxStrainVector(6); // 6 components
        auxStrainVector[0] = StrainVector[0]-mRetraction*(0.8*StrainVector[0]);
        auxStrainVector[1] = StrainVector[1]-mRetraction*(0.8*StrainVector[1]);
        auxStrainVector[2] = Ez;
        auxStrainVector[3] = StrainVector[2];
        auxStrainVector[4] = 0.0;
        auxStrainVector[5] = 0.0;
        //KRATOS_WATCH(auxStrainVector);
        //Vector auxStressVector(6); // 6 components

        //// IN TOTAL LAGRANGIAN! (BELOW)
        C(0,0) = 2.0 * auxStrainVector[0] + 1.0;
        C(1,1) = 2.0 * auxStrainVector[1] + 1.0;
        C(2,2) = 2.0 * auxStrainVector[2] + 1.0;
        C(0,1)=auxStrainVector[3];
        C(1,2)=auxStrainVector[4];
        C(0,2)=auxStrainVector[5];
        C(1,0)=C(0,1);
        C(2,1)=C(1,2);
        C(2,0)=C(0,2);
        //KRATOS_WATCH(C)

        Matrix Cinv(3,3);
        noalias(Cinv)=ZeroMatrix(3,3);
        double det=0.0; // det of C  ( det(C)=J2 )
        //MathUtils<double>::InvertMatrix3(Celastic,Cinv,det);
        MathUtils<double>::InvertMatrix3(C,Cinv,det); // Cinv and Det(C) are computed

        double J2=MathUtils<double>::Det3( C ); // determinant of C is computed
        //KRATOS_WATCH(J2);
        double J=0.0;
        J=sqrt(J2); // J is computed
        //KRATOS_WATCH(J);
        // Up: computing C and J
        //*************************************************
        // computing I4 = C_bar : A
        // where C_bar = J_23 * C
        // and   A = a0 vs a0
        double I4 =0.0; // I4 = A:C' and A= a0 x a0
        // computing C_bar = J_23 * C
        Matrix C_bar(3,3);
        noalias(C_bar)=C;
        double J_23 = pow(J,(-2.0/3.0));
        C_bar *=J_23;
        // computing A= a0 x a0
        Vector a(3);
        noalias(a)=anl;
        Vector b(3);
        noalias(b)=anl;
        // computing cross product
        Vector a_vs_a(6); // cross product of ae vector, Christian calls it "axa"
        //int i=0;
// 			      for  (i=0; i<3; i++)
// 				  {
// 				    // tau is my total stress (StressVector)
// 				    a_vs_a(i)=a(i)*a(i);
// 				  }
// 			      a_vs_a[3]=a[0]*a[1];
// 			      a_vs_a[4]=a[1]*a[2];
// 			      a_vs_a[5]=a[0]*a[2];
// 			      // in matrix form
        Matrix A(3,3);
// 			      Matrix A = MathUtils<double>::StressVectorToTensor(  a_vs_a );
        A(0,0)=a[0]*b[0];
        A(0,1)=a[0]*b[1];
        A(0,2)=a[0]*b[2];
        A(1,0)=a[1]*b[0];
        A(1,1)=a[1]*b[1];
        A(1,2)=a[1]*b[2];
        A(2,0)=a[2]*b[0];
        A(2,1)=a[2]*b[1];
        A(2,2)=a[2]*b[2];
        SD_MathUtils<double>::TensorToVector(A,a_vs_a);
        //TensorProduct3(a, a, A);
        //Matrix A = MathUtils<double>::TensorProduct3(a, a);
        //KRATOS_WATCH(A);
        // finally obtaining I4 = C_bar:A
        I4+=C_bar(0,0)*A(0,0)+C_bar(0,1)*A(1,0)+C_bar(0,2)*A(2,0);
        I4+=C_bar(1,0)*A(0,1)+C_bar(1,1)*A(1,1)+C_bar(1,2)*A(2,1);
        I4+=C_bar(2,0)*A(0,2)+C_bar(2,1)*A(1,2)+C_bar(2,2)*A(2,2);
        //KRATOS_WATCH(I4);
        // computing tau_f
        if ( I4 < 1.0 )
        {
            noalias(tau_f)=ZeroVector(6);
        }
        else
        {

            //int ii=0;
            //for  (ii=0; ii<2; ii++)
            //{
            //Vector tau_f(6); // initialize the stress contribution of each fiber
            noalias(tau_f) = ZeroVector(6);
            //*****************************
            // computing w4
            // CALL PSI I4 !! 3rd function : derivative of W in respect to I4
            double I4_1 = I4-1.0; // I4 -1
            //KRATOS_WATCH (mK1);
            double temp= exp(mK2*I4_1*I4_1);
            // first derivative of the strain energy function in respect to the 4th invariant
            double w4 = mK1*I4_1*temp;
            //*****************************
            // computing tau_f in Second Piola
            Matrix C_bar_inv(3,3);
            noalias(C_bar_inv)=ZeroMatrix(3,3);
            double det_bar = 0.0;
            MathUtils<double>::InvertMatrix3(C_bar, C_bar_inv, det_bar);
            Vector vC_bar_inv(6);
            SD_MathUtils<double>::TensorToVector(C_bar_inv,vC_bar_inv);
            tau_f=a_vs_a;
            tau_f-=(1.0/3.0)*I4*vC_bar_inv;
            tau_f*=2.0*J_23*w4;
            //KRATOS_WATCH(tau_f)
            // second derivative of the stored energy function in respect to the 4th invariant
            //double w44 =  mK1*(1.0+(2.0*mK2*I4_1*I4_1))*temp; // It's used only for the tangent. So I won't use it!
            //*************************************************************
            // CALL vect_product ae !! 4th function : deviator of "a"
            // CALL dev_euler !! 5th function
            // 				Vector dev_a_vs_a(6); // check the size of this vector !!
            // 				noalias(dev_a_vs_a)=a_vs_a;
            // 				i=0;
            // 				double one3 = 1.0/3.0;
            // 				double tra= a_vs_a(0)+a_vs_a(1)+a_vs_a(2);
            // 				tra=one3*tra;
            // 				for  (i=0; i<3; i++) //
            // 				      {
            // 					// tau is my total stress (StressVector)
            // 					dev_a_vs_a(i)=dev_a_vs_a(i)-tra;
            // 				      }
            // 				Matrix C(2,2); // creating matrix C (C=Ft.F)
            // 				noalias(C)=ZeroMatrix(2,2);
            // 				noalias(C)=prod(trans(DeformationGradient),DeformationGradient);
            // 				Matrix Cinv(3,3);
            // 				noalias(Cinv)=ZeroMatrix(3,3);
            // 				double det=0.0; // det of C  ( det(C)=J2 )
            // 				MathUtils<double>::InvertMatrix2(C,Cinv,det); // Cinv and Det(C) are computed
            // 				Vector vCinv(3); // 6 components
            // 				SD_MathUtils<double>::TensorToVector(Cinv,vCinv);
            // 				//tau_f=2.0*w4*dev_ae_vs_ae; This is the Kirchoff stress, I need Second Piola Kirchhoff
            // 				noalias(tau_f) = dev_a_vs_a;
            // 				noalias(tau_f) -= one3*I4*vCinv;
            // 				tau_f *= 2.0*w4*J_23;
            // 				//tau_f=2.0*w4*J_23*(dev_a_vs_a-(one3*I4*Cinv));

            // END TEST
            //tau(0)=tau(0)+tau_f(0);
// 				      tau[0] = tau_f[0]+tau[0];
// 				      tau[1] = tau_f[1]+tau[1];
// 				      tau[2] = tau_f[3]+tau[2];

        }
        //tau=tau+tau_f;
        tau[0] = tau_f[0]+tau[0];
        tau[1] = tau_f[1]+tau[1];
        tau[2] = tau_f[3]+tau[2];
    }

    noalias(rTau)=tau;
    //KRATOS_WATCH(tau)


}
//*******************************************************************************************
// 	void Viscofibers2D::TensorProduct3(const Vector& a, const Vector& b, Matrix& A)
//             {
//                 MatrixType A(3,3);
//                 A(0,0)=a[0]*b[0]; A(0,1)=a[0]*b[1]; A(0,2)=a[0]*b[2];
//                 A(1,0)=a[1]*b[0]; A(1,1)=a[1]*b[1]; A(1,2)=a[1]*b[2];
//                 A(2,0)=a[2]*b[0]; A(2,1)=a[2]*b[1]; A(2,2)=a[2]*b[2];
//                 //return A;
//             }

/*	void Viscofibers2D::CalculateFiber(const Vector& ane, Vector& tau_f)
		{

		  Vector ae(3);
		  noalias(ae)=ane; // copy ane in ae ! Is it correct ????

		  double I4 =0.0;
		  double I4_1 =0.0;
		  int j=0;
		  int i=0;
		  for  (i=0; i<3; i++) //
			{
			   // tau is my total stress (StressVector)
			   I4=I4+ae(i)*ae(i); // "ii" is  each fiber direction, "j" is each component, "n" is each element? I think I don't need this n!
			}
		  if ( I4 < 1.0 )
		  {
		    //tau=0.0; // which tau is this one?
		    tau_f=ZeroVector(3);
		  }
		  else
		  {

		    // CALL PSI I4 !! 3rd function : derivative of W in respect to I4
		    I4_1 = I4-1.0; // I4 -1

		    double temp= exp(mK2*I4_1*I4_1);
		    // first derivative of the strain energy function in respect to the 4th invariant
		    double w4 = mK1*I4_1*temp;
		    // second derivative of the stored energy function in respect to the 4th invariant
		    double w44 =  mK1*(1.0+(2.0*mK2*I4_1*I4_1))*temp; // It's used only for the tangent. So I won't use it!

		    // CALL vect_product ae !! 4th function : deviator of "a"
		    Vector ae_vs_ae(6); // cross product of ae vector, Christian calls it "axa"
		    int i=0;
		    for  (i=0; i<3; i++)
			{
			   // tau is my total stress (StressVector)
			   ae_vs_ae(i)=ae(i)*ae(i);
			}
		    ae_vs_ae[3]=ae[0]*ae[1];
		    ae_vs_ae[4]=ae[1]*ae[2];
		    ae_vs_ae[5]=ae[0]*ae[2];

		    // CALL dev_euler !! 5th function
		    Vector dev_ae_vs_ae(6); // check the size of this vector !!
		    noalias(dev_ae_vs_ae)=ae_vs_ae;
		    i=0;
		    double one3 = 1.0/3.0;
		    double tra= ae_vs_ae(0)+ae_vs_ae(1)+ae_vs_ae(2);
		    tra=one3*tra;
		     for  (i=0; i<3; i++) //
			{
			   // tau is my total stress (StressVector)
			   dev_ae_vs_ae(i)=dev_ae_vs_ae(i)-tra;
			}

		    //tau_f=2.0*w4*dev_ae_vs_ae; This is the Kirchoff stress, I need Second Piola Kirchhoff
		    tau_f=2.0*w4*J_23*a0xa0-one3*I4*C-1; //dev_ae_vs_ae;
		  }

		}*/
//*******************************************************************************************
/**
Calculate Constitutive Matrix, using the pertubation method
*/
void Viscofibers2D::CalculateConstitutiveMatrix(const Vector& StrainVector, Matrix& rResult)

{
    if( rResult.size1() != 3 ) // EBST
    {
        rResult.resize(3,3);
    }
//
    Vector E1(6);
    Vector S1(6);
    Vector E2(6);
    Vector S2(6);
    // rTau
// 			 Vector auxrTau(6); // 6 components
// 			 auxrTau[0] = rTau[0];
// 			 auxrTau[1] = rTau[1];
// 			 auxrTau[2] = 0.0;
// 			 auxrTau[3] = rTau[2];
// 			 auxrTau[4] = 0.0;
// 			 auxrTau[5] = 0.0;
//
    // Vector E1;
// 			array_1d<double,6> S1;
// 			Vector E2;
// 			array_1d<double,6> S2;


    ///////////// TANGENT MATRIX /////////// SERGIO OLLER + NELSON

    // PERTUBATION METHOD
    long double dE=0.00; // delta strain
    //double epsilon_real=std::numeric_limits<double>::epsilon(); //pertubation
    //KRATOS_WATCH(epsilon_real);
    long double epsilon=1E-9; // pertubation
    long double max=1E-14;
    E1.resize(6,false);
    //S1.resize(6,false);
    E2.resize(6,false);
    //S2.resize(6,false);

    // for EBST and Membrane
    mcurrentThickness=(mA0/mA)*mThickness;
    //KRATOS_WATCH(mA0);
    //KRATOS_WATCH(mA);
    //KRATOS_WATCH(mThickness);
    //KRATOS_WATCH(mcurrentThickness);
    double Ez = 0.0;
    // Uncomment below if you want to use Ez (uncomment the part related to y too! - line 587)
    Ez = (1.0/2.0)*((pow(mcurrentThickness,2.0)-pow(mThickness,2.0))/pow(mThickness,2.0));
    //KRATOS_WATCH(Ez);

    Vector auxStrainVector(6); // 6 components
    auxStrainVector[0] = StrainVector[0]-mRetraction*(0.8*StrainVector[0]);  // 0.8*StrainVector[0] represents E0
    auxStrainVector[1] = StrainVector[1]-mRetraction*(0.8*StrainVector[1]);
    auxStrainVector[2] = Ez;
    auxStrainVector[3] = StrainVector[2];
    auxStrainVector[4] = 0.0;
    auxStrainVector[5] = 0.0;
    //KRATOS_WATCH (auxStrainVector);
    Vector auxStressVector(6); // 6 components

    E1=auxStrainVector;
    //KRATOS_WATCH (E1);
    E2=auxStrainVector;
    //KRATOS_WATCH (E2);
    Matrix Ctot(6,6);  //
    noalias(Ctot)=ZeroMatrix(6,6);  // before Ctang
    //KRATOS_WATCH (Ctot);

    //// NEW FORM
    dE=(*std::min_element(auxStrainVector.begin(),auxStrainVector.end()))*epsilon;
    if (dE==0.00)
    {
        dE=epsilon;
    }
    if (dE<max)
    {
        dE=max;
    }
    for (unsigned int i=0; i<6; i++)
    {
        E2(i)+=dE; //pertubed Strain E2
        Viscofibers2D::CalculateStressVector(E1,S1); // computing stress S1
        Viscofibers2D::CalculateStressVector(E2,S2); // computing pertubed stress S2
        noalias(S2)=S2-S1;
        noalias(S2)=S2/(dE);
        for (unsigned int j=0; j<6; j++)
        {
            Ctot(j,i)=S2(j); // Tangent Matrix

        }

        //E2.resize(6,false);
        //S2.resize(6,false);
        noalias(E2) = ZeroVector(6);
        noalias(S2) = ZeroVector(6);
        //noalias(E2) = StrainVector;
    }
    //noalias(rResult)=Ctang; // 2D
    //KRATOS_WATCH (Ctot);
    //noalias(rResult)=Ctang; // 2D

    // TEST !!!
    Matrix auxCtang2(4,4);
    for (unsigned int i=0; i<4; i++)
    {
        for (unsigned int j=0; j<4; j++)
        {
            auxCtang2(i,j)=Ctot(i,j); // aux Tangent Matrix 4X4

        }
    }
    //KRATOS_WATCH (auxCtang2);
    // below putting zeros in last line and column
    Matrix auxCtang3(4,4);
    noalias(auxCtang3)=auxCtang2;
    for (unsigned int i=0; i<4; i++)
    {
        auxCtang3(2,i)=auxCtang2(3,i);
        auxCtang3(3,i)=auxCtang2(2,i);
        auxCtang3(i,2)=auxCtang2(i,3);
        auxCtang3(i,3)=auxCtang2(i,2);
    }
    //auxCtang3(2,2)=auxCtang2(3,3);
    //auxCtang3(2,3)=auxCtang2(3,2);
    //auxCtang3(3,2)=auxCtang2(2,3);
    //auxCtang3(3,3)=auxCtang2(2,2);
    //KRATOS_WATCH (auxCtang3(3,3));
    // end - Matrix 4X4 with last line and column with zeros (ref Sigma_Z)


    Matrix v(3,3);
    for (unsigned int i=0; i<3; i++)
    {
        for (unsigned int j=0; j<3; j++)
        {
            v(i,j)=auxCtang3(i,j);
        }
    }
    //KRATOS_WATCH (v);

    Matrix w(3,1);
    for (unsigned int i=0; i<3; i++)
    {
        w(i,0)=auxCtang3(i,3);
    }
    //KRATOS_WATCH (w);

    Matrix x(1,3);
    for (unsigned int i=0; i<3; i++)
    {
        x(0,i)=auxCtang3(3,i);
    }
    //KRATOS_WATCH (x);

    double y=auxCtang3(3,3);
    //KRATOS_WATCH (y);

    Matrix auxCtang(3,3);
    noalias(auxCtang) += v;
    //KRATOS_WATCH (auxCtang);
    Matrix auxtest(3,3);
    noalias(auxtest) = prod(w,x);
    //KRATOS_WATCH (auxtest);
// 			if (y != 0.0)
// 			{
// 			auxtest *= (1.0)/y;
// 			}
// 			else
// 			{
// 			y = 0.000000000000000000001;
// 			auxtest *= (1.0)/y;
// 			}
    // Uncomment below if considering Ez
    auxtest *= (1.0)/y;
    //KRATOS_WATCH (auxtest);
    noalias(auxCtang) -= auxtest;
    //KRATOS_WATCH (auxCtang);

    noalias(rResult)=auxCtang;
    //KRATOS_WATCH (rResult);
}

//*************************************************************************************************************************
/**
Update Material
*/
void Viscofibers2D::UpdateMaterial( const Vector& StrainVector,
                                    const Properties& props,
                                    const GeometryType& geom,
                                    const Vector& ShapeFunctionsValues,
                                    const ProcessInfo& CurrentProcessInfo)
{
    //noalias(mcurrent_hn1) = mhn1; //original

}

//*************************************************************************************************************************

void Viscofibers2D::CalculateStressAndTangentMatrix(Vector& StressVector,  // check what is this function for!
        const Vector& StrainVector,
        Matrix& algorithmicTangent)
{
    // CalculateConsitutiveMatrix(StrainVector, algorithmicTangent);
}


//*************************************************************************************************************************
/**
Calculate Cauchy Stresses
*/
void Viscofibers2D::CalculateCauchyStresses(
    Vector& rCauchy_StressVector,
    const Matrix& rF,
    const Vector& rPK2_StressVector, // in my case I already have S in a matrix form! How to read it?
    const Vector& rGreenLagrangeStrainVector)
{
    //KRATOS_WATCH (rPK2_StressVector);
    Matrix S = MathUtils<double>::StressVectorToTensor( rPK2_StressVector );
    //KRATOS_WATCH("Second Piola en HYP2D")
    //KRATOS_WATCH (S);
    //KRATOS_WATCH(rF)
    double J = MathUtils<double>::Det2( rF );
    //KRATOS_WATCH("MATRIZ F em HYP2D")
    //KRATOS_WATCH (J);
    boost::numeric::ublas::bounded_matrix<double,2,2> mstemp;
    boost::numeric::ublas::bounded_matrix<double,2,2> msaux;


    noalias(mstemp) = prod(rF,S);
    //KRATOS_WATCH (mstemp);
    noalias(msaux) = prod(mstemp,trans(rF));
    //KRATOS_WATCH (msaux);
    msaux /= J;  // good for Membrane


    if(rCauchy_StressVector.size() != 3)
        rCauchy_StressVector.resize(3, false);

    rCauchy_StressVector[0] = msaux(0,0);
    rCauchy_StressVector[1] = msaux(1,1);
    rCauchy_StressVector[2] = msaux(0,1);
    //KRATOS_WATCH("CAUCHY IN VISCO-HYPERELASTIC")
    //KRATOS_WATCH(rCauchy_StressVector)
}

} // Namespace Kratos

//***********************************************************************************************************
//***********************************************************************************************************
//***********************************************************************************************************
//***********************************************************************************************************
//***********************************************************************************************************
//***********************************************************************************************************
//***********************************************************************************************************

/// UP ///////// TANGENT MATRIX


