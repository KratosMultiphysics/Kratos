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

THE  SOFTWARE IS  PROVIDED  "AS  variablesIS", WITHOUT  WARRANTY  OF ANY  KIND,
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
*   Last Modified by:    $Author: nelson $
*   Date:                $Date: 2008-09-03 
\*   Revision:            $Revision: 1.1 $
*
* ***********************************************************/


// System includes 
#include <iostream>
#include <iomanip>
#include <cmath>

// Project includes 
#include "constitutive_laws/plastic_damage_2d.h"



namespace Kratos
{
	/**
	 *	TO BE TESTED!!!
	 */

         PlasticDamage2D::PlasticDamage2D () 
	: ConstitutiveLaw()
	
	{
	  KRATOS_ERROR(std::logic_error,"Calling the empty constructor.","");
	}

	 PlasticDamage2D::PlasticDamage2D(
         FluencyCriteriaPointer FluencyCriteria,
         FluencyCriteriaPointer FluencyCriteriaTraction,
         SofteningHardeningCriteriaPointer SofteningBehavior, 
         PropertiesPointer Property) 
	: ConstitutiveLaw()
	{
	      mpFluencyCriteria              = FluencyCriteria;
              mpFluencyCriteria_Traction     = FluencyCriteriaTraction;
              mpSofteningBehavior            = SofteningBehavior;
              mpProperties                   = Property;

	}
	
	
	
	/**
	 *	TO BE TESTED!!!
	 */
	PlasticDamage2D::~PlasticDamage2D ()
	{
	}
	
	
	bool PlasticDamage2D::Has( const Variable<double>& rThisVariable )
	{
		return false;
	}
	
	bool PlasticDamage2D::Has( const Variable<Vector>& rThisVariable )
	{
		return false;
	}
	
	bool PlasticDamage2D::Has( const Variable<Matrix>& rThisVariable )
	{
		return false;
	}
	
	double& PlasticDamage2D::GetValue( const Variable<double>& rThisVariable, double& rValue )
	{
	  if( rThisVariable == DAMAGE)
	  {
             rValue = mlocal_fail_factor;  //mplastic_damage ;}
	  }
          else if(rThisVariable == COHESION)
	  
             rValue = mFt[1]; // mcohesion;}
          else if(rThisVariable == DILATANCY_ANGLE)
             rValue = mdilatancy_angle*180.00/PI;
          else if(rThisVariable == FRICTION_INTERNAL_ANGLE)
             rValue = mFt[2]; // mfriction_angle*180.00/PI;}
          return( rValue );
	}   
	

	Vector& PlasticDamage2D::GetValue( const Variable<Vector>& rThisVariable, Vector& rValue )
	{
        return( rValue );
         }
    
    Matrix& PlasticDamage2D::GetValue( const Variable<Matrix>& rThisVariable, Matrix& rValue )
    {
        return( rValue );
	}

    void PlasticDamage2D::SetValue( const Variable<double>& rThisVariable, const double& rValue, 
                                const ProcessInfo& rCurrentProcessInfo )
	{
	}
	
    void PlasticDamage2D::SetValue( const Variable<Vector>& rThisVariable, const Vector& rValue, 
                                const ProcessInfo& rCurrentProcessInfo )
	{
	}
	
    void PlasticDamage2D::SetValue( const Variable<Matrix>& rThisVariable, const Matrix& rValue, 
                                const ProcessInfo& rCurrentProcessInfo )
	{
	}
    
    void PlasticDamage2D::Calculate(const Variable<Matrix >& rVariable, Matrix& rResult, 
                   const ProcessInfo& rCurrentProcessInfo)
    {
      for (unsigned int ii = 0; ii<3; ii++)
           rResult(0,ii) = mplastic_strain(ii); 
    }

    void PlasticDamage2D::Calculate(const Variable<double>& rVariable, 
                                    double& Output, 
                                    const ProcessInfo& rCurrentProcessInfo)
   {
     Output = sqrt(mE/mDE);
   }

//***********************************************************************************************
//***********************************************************************************************

	void PlasticDamage2D::InitializeMaterial( const Properties& props,
	const GeometryType& geom,
	const Vector& ShapeFunctionsValues )
	
	{
	  
	  mE                       = (*mpProperties)[YOUNG_MODULUS];
	  mNU                      = (*mpProperties)[POISSON_RATIO];
	  mDE                      = (*mpProperties)[DENSITY];
	  mlength                  = 1.1283791670955 * sqrt(fabs(geom.Area()));
	  
          mFt[0]                   = (*mpProperties)[FT]; 
          mFt[1]                   = (*mpProperties)[FT];         
          mFt[2]                   = (*mpProperties)[FT]; 

          mplastic_damage          = 0.00;
          mdisipation              = 0.00;
          mefective_plastic_strain = 0.00;
          mlocal_fail_factor       = 0.00; 
	          
          noalias(mplastic_strain)         = ZeroVector(4);
          noalias(mcurrent_plastic_strain) = ZeroVector(4);  

          double Gc           = (*mpProperties)[CRUSHING_ENERGY]/mlength;
          double length_limit = 2.00*mE*Gc/((*mpProperties)[FC]*(*mpProperties)[FC]);       

          if (length_limit<mlength) {std::cout<<"Element length greater than permitted"<<std::endl;}  
          
          mpFluencyCriteria->InitializeMaterial(*mpProperties);


	}

		

//***********************************************************************************************
//***********************************************************************************************

void PlasticDamage2D::InitializeSolutionStep( const Properties& props,
                    const GeometryType& geom,
                    const Vector& ShapeFunctionsValues ,
                    const ProcessInfo& CurrentProcessInfo)
{
     
}


//***********************************************************************************************
//***********************************************************************************************

void PlasticDamage2D::FinalizeSolutionStep( const Properties& props,
		  const GeometryType& geom, 
		  const Vector& ShapeFunctionsValues ,
		  const ProcessInfo& CurrentProcessInfo)
{

        mefective_plastic_strain =  mcurrent_efective_plastic_strain;
        noalias(mplastic_strain) =  mcurrent_plastic_strain;
        mplastic_damage          =  mcurrent_plastic_damage;
        mdisipation              =  mcurrent_disipation;   
        noalias(mFt)             =  mcurrent_Ft;
        
        
        array_1d<double,4> Variables;
        Variables[0] =  mFt[0];
        Variables[1] =  mFt[1];
        Variables[2] =  mFt[2];
        Variables[3] =  mlength;
   
        mpFluencyCriteria->Finalize();   
	mpFluencyCriteria->UpdateVariables(Variables); 
        
	
        
        Matrix Aux_V         = ZeroMatrix(3,3);  
        Matrix PlasticTensor = ZeroMatrix(3,3);  
        Vector Aux_P         = ZeroVector(3); 

      	PlasticTensor(0,0)   = mplastic_strain[0];     PlasticTensor(0,1) = 0.5*mplastic_strain[2]; PlasticTensor(0,2) = 0.00;
	PlasticTensor(1,0)   = 0.5*mplastic_strain[2]; PlasticTensor(1,1) = mplastic_strain[1];     PlasticTensor(1,2) = 0.00;
	PlasticTensor(2,0)   = 0.00;                   PlasticTensor(2,1) = 0.00;                   PlasticTensor(2,2) = mplastic_strain[3];    
	
        SD_MathUtils<double>::EigenVectors(PlasticTensor, Aux_V, Aux_P, 1E-9, 100);
	
	const double Euc   =  2.00 * (*mpProperties)[FRACTURE_ENERGY]/ ( (*mpProperties)[FT] * mlength);
        double Ef          = (*max_element(Aux_P.begin(), Aux_P.end())) + (*mpProperties)[FT]/(*mpProperties)[YOUNG_MODULUS];
        mlocal_fail_factor =  Ef / Euc; 
        if (mlocal_fail_factor > 1.00) {mlocal_fail_factor = 1.00; }  
        
  
}

//***********************************************************************************************
//***********************************************************************************************


void PlasticDamage2D::CalculateElasticMatrix(boost::numeric::ublas::bounded_matrix<double,4,4>& C)
{ 

    /// plane strain and axial symmetric
    double c  =  mE / ((1.00 + mNU)*(1.00-2.00*mNU));
    double c1 =  (1.00 - mNU) * c;
    double c2 =  mNU * c ;    
    double c3 =  0.50 * (1.00 - 2.00*mNU) * c;
    C(0,0) = c1;  C(0,1) = c2;	C(0,2) = 0.0;  C(0,3)  = c2;
    C(1,0) = c2;  C(1,1) = c1;	C(1,2) = 0.0;  C(1,3)  = c2;
    C(2,0) = 0.0; C(2,1) = 0.0;	C(2,2) = c3;   C(2,3)  = 0.00;
    C(3,0) = c2;  C(3,1) = c2;	C(3,2) = 0.00; C(3,3)  = c1;

}    


//***********************************************************************************************
//***********************************************************************************************

void PlasticDamage2D::CalculateElasticStress(array_1d<double,4>& Strain, array_1d<double,4>& Stress)
{
/// plane strain and axial symmetric 
double G          = 0.5*mE / (1.00 + mNU);
double K          = mE / (3.00 * (1.00-2.00*mNU) );
double vol_strain = Strain[0] + Strain[1] + Strain[3];
double pt         = K * vol_strain;     

Stress[0] = 2.00 * G *(Strain[0] - vol_strain/3.00 ) + pt; 
Stress[1] = 2.00 * G *(Strain[1] - vol_strain/3.00 ) + pt;
Stress[2] = G *(Strain[2]);
Stress[3] = 2.00 * G *(Strain[3] - vol_strain/3.00 ) + pt;

}


//***********************************************************************************************
//***********************************************************************************************
void PlasticDamage2D::CalculateMaterialResponse( const Vector& StrainVector,
        const Matrix& DeformationGradient,
        Vector& StressVector,
        Matrix& AlgorithmicTangent,
        const ProcessInfo& CurrentProcessInfo,
        const Properties& props, 
        const GeometryType& geom,
        const Vector& ShapeFunctionsValues,
        bool CalculateStresses,
        int CalculateTangent,
        bool SaveInternalVariables
                                               )
{
  
    UpdateMaterial(StrainVector, props, geom,ShapeFunctionsValues, CurrentProcessInfo);
    if (CalculateStresses==true) { CalculateStress(StrainVector, StressVector);}
    if (CalculateTangent==true){CalculateStressAndTangentMatrix(StressVector,StrainVector, AlgorithmicTangent);}
}


void PlasticDamage2D::CalculateConstitutiveMatrix(const Vector& StrainVector, Matrix& ConstitutiveMatrix)
{

  /// plane strain and axial symmetric  
  ConstitutiveMatrix.resize(3,3, false);
  double  c1 = mE*(1.00-mNU)/((1.00 + mNU)*(1.00-2.00*mNU));
  double  c2 = c1*mNU/(1.00 - mNU);
  double  c3 = 0.5*mE / (1.00 + mNU);

  ConstitutiveMatrix(0,0) = c1;  ConstitutiveMatrix(0,1) = c2;  ConstitutiveMatrix(0,2) = 0.0;
  ConstitutiveMatrix(1,0) = c2;  ConstitutiveMatrix(1,1) = c1;  ConstitutiveMatrix(1,2) = 0.0;
  ConstitutiveMatrix(2,0) = 0.0; ConstitutiveMatrix(2,1) = 0.0; ConstitutiveMatrix(2,2) = c3;
}


//***********************************************************************************************
//***********************************************************************************************

void PlasticDamage2D::CalculateStress(const Vector& StrainVector, 
		  Vector& StressVector)
{
 
  	double ElasticDomain_1               = 0.00;
        double toler                         = 1E-6;      
	array_1d<double,4> ElasticStrain     = ZeroVector(4);                                           
	array_1d<double,4> ElasticStress     = ZeroVector(4);     
        StressVector                         = ZeroVector(3); 
	

        /// calculating elastic strain
        /// plane strain and axial symmetric 
	ElasticStrain[0] =  StrainVector[0];
	ElasticStrain[1] =  StrainVector[1];
	ElasticStrain[2] =  StrainVector[2];
	ElasticStrain[3] =  mcurrent_plastic_strain[3];   // et = ep; ez elastico siempre tiene que ser igual a cero 
 
       
        ///* The elastic strain trial
	noalias(ElasticStrain) = StrainVector - mcurrent_plastic_strain;

	///*calculating elastic stress trial
	CalculateElasticStress(ElasticStrain, ElasticStress);
         

        ///* Comprobado criterio de fluencia   
        mpFluencyCriteria->CalculateEquivalentUniaxialStress(ElasticStress, ElasticDomain_1);     
        
        ///* dominio elastico    
	if( ElasticDomain_1 <= toler ) 
	{
              noalias(mcurrent_Ft)             = mFt;
	      mcurrent_efective_plastic_strain = mefective_plastic_strain;
	      noalias(mcurrent_plastic_strain) = mplastic_strain;
	      mcurrent_plastic_damage          = mplastic_damage;
	      mComputeTangentMatrix            = false;
	}

        else
         {
            mComputeTangentMatrix = false;
	    
	    ///* Spectral Descomposition 
	    Compute_Principal_Stress(ElasticStress);
	    array_1d<double, 3 > PrincipalStress;
            array_1d<double, 3 > Sigma;
	    array_1d<unsigned int, 3 > Order;
            Vector delta_gamma;    
 

            ///* 
            noalias(PrincipalStress) = mPrincipalStress;   

	    ///* Guardando el orden de los cambios
	    IdentifyMaximunAndMinumumPrincipalStres_CalculateOrder(PrincipalStress , Order);

            ///* Regresando ala superficie de fluencia 
            mpFluencyCriteria->ReturnMapping(ElasticStress, ElasticStrain, delta_gamma, Sigma);
              
            ///* Calculando las tensiones elsaticas
	    AssembleUpdateStressAndStrainTensor(Sigma, StrainVector, Order, ElasticStrain, ElasticStress); 
            //KRATOS_WATCH(Sigma)  
	      
	    ///* General evoluation of accumulated hardening  
	    mcurrent_efective_plastic_strain = mefective_plastic_strain + norm_1(delta_gamma); 

	    ///* aculated Von mises plastic
              //double aux_var = sqrt(delta_gamma_a*delta_gamma_a + delta_gamma_b * delta_gamma_b + delta_gamma_c * delta_gamma_c);
	      //mcurrent_efective_plastic_strain = mefective_plastic_strain + (2.00 / sqrt(3.00) ) * aux_var;
 
           ///* updating current variables
            Vector Result(3);
            mpFluencyCriteria->GetValue(Result); noalias(mcurrent_Ft) = Result; 
             

	 }
           StressVector[0] = ElasticStress[0]; // Gxx
	   StressVector[1] = ElasticStress[1]; // Gyy
	   StressVector[2] = ElasticStress[2]; // Txy 

 }     




//***********************************************************************************************
//***********************************************************************************************


void PlasticDamage2D::AssembleUpdateStressAndStrainTensor(
        const array_1d<double,3>& Sigma,
        const array_1d<double,4>& StrainVector_aux,   
        const array_1d<unsigned int,3>& order,
	array_1d<double,4>& ElasticStrain,                                                
	array_1d<double,4>& ElasticStress)     
{
    
Matrix StressTensor      = ZeroMatrix(3,3);
Matrix DesviatoricTensor = ZeroMatrix(3,3);
Matrix VolumnetricTensor = ZeroMatrix(3,3);
Matrix StrainTensor      = ZeroMatrix(3,3); 
const identity_matrix<double> I (3);

//KRATOS_WATCH(order) 

/// Updating the  elastic stress tensor
for(unsigned int i=0; i<3; i++)
{                     
noalias(StressTensor) = StressTensor + Sigma[i] * Matrix(outer_prod(mEigenVectors[order[i]],  mEigenVectors[order[i]]));
} 
                          
/// Warning Solo 2D
ElasticStress[0] = StressTensor(0,0);  //xx
ElasticStress[1] = StressTensor(1,1);  //yy
ElasticStress[2] = StressTensor(0,1);  //xy
ElasticStress[3] = StressTensor(2,2);  //zz


/// Updating the elastic tensor   
double p  = (ElasticStress[0] + ElasticStress[1] + ElasticStress[3])/3.00;  
double G  = 0.5*mE / (1.00 + mNU);
double K  = mE / (3.00 * (1.00-2.00*mNU) );
//noalias(DesviatoricTensor) = StressTensor - p * I;
//noalias(StrainTensor)      = (1.00/(2.00 * G ) ) * DesviatoricTensor + p/(3.00 * K) * I;   

ElasticStrain[0] = (ElasticStress[0] - p ) / (2.00 * G) + p/(3.00 * K);       //StrainTensor(0,0);       //xx
ElasticStrain[1] = (ElasticStress[1] - p ) / (2.00 * G) + p/(3.00 * K);       //StrainTensor(1,1);       //yy
ElasticStrain[2] =  ElasticStress[2]/G;                                       //2.00*StrainTensor(0,1);  //xy
ElasticStrain[3] = (ElasticStress[3] - p )/  (2.00 * G) + p/(3.00 * K);       //StrainTensor(2,2);       //zz  
noalias(mcurrent_plastic_strain) = StrainVector_aux - ElasticStrain; 

//Updated_Internal_Variables(ElasticStress, mcurrent_plastic_strain);
//Tensile_Fracture_Model(ElasticStress,  mcurrent_plastic_strain);


}

//***********************************************************************************************
//***********************************************************************************************

void PlasticDamage2D::IdentifyMaximunAndMinumumPrincipalStres_CalculateOrder(array_1d<double, 3 >& PrincipalStress ,array_1d<unsigned int,3>& order)
{
    order[0] = 0; 
    order[2] = 0;
    PrincipalStress[0] = mPrincipalStress[order[0]];
    PrincipalStress[2] = mPrincipalStress[order[2]];
    for (unsigned int i = 1; i<3; i++)
    {
    if (mPrincipalStress[i]>=PrincipalStress[0]) {order[0] = i; PrincipalStress[0] = mPrincipalStress[order[0]];}
    if (mPrincipalStress[i]<=PrincipalStress[2]) {order[2] = i; PrincipalStress[2] = mPrincipalStress[order[2]];}
    }

    if(order[0]!=0 && order[2]!=0){order[1] = 0;}
    if(order[0]!=1 && order[2]!=1){order[1] = 1;}
    if(order[0]!=2 && order[2]!=2){order[1] = 2;}
    PrincipalStress[1] = mPrincipalStress[order[1]];

    if(PrincipalStress[0] == PrincipalStress[1]) {order[0]= 0; order[1]= 1; order[2]= 2;} 


   /// Evitando redondeos matematicos absurdos
   if ( fabs(PrincipalStress[0]) < 1E-3 ) {PrincipalStress[0] = 0.00; } 
   if ( fabs(PrincipalStress[1]) < 1E-3 ) {PrincipalStress[1] = 0.00; } 
   if ( fabs(PrincipalStress[2]) < 1E-3 ) {PrincipalStress[2] = 0.00; } 

}


//***********************************************************************************************
//***********************************************************************************************


void PlasticDamage2D::CalculateCauchyStresses(
		Vector& rCauchy_StressVector,
		const Matrix& rF,
		const Vector& rPK2_StressVector,
		const Vector& rGreenLagrangeStrainVector)
{
               
		Matrix S = MathUtils<double>::StressVectorToTensor( rPK2_StressVector );

                double J = MathUtils<double>::Det2( rF );
        boost::numeric::ublas::bounded_matrix<double,2,2> temp;
        boost::numeric::ublas::bounded_matrix<double,2,2> aux;
		noalias(temp) = prod(rF,S);
		noalias(aux) = prod(temp,trans(rF));
		aux *= J;

		if(rCauchy_StressVector.size() != 3)
			rCauchy_StressVector.resize(3);  
		
		rCauchy_StressVector[0] = aux(0,0);
		rCauchy_StressVector[1] = aux(1,1);
		rCauchy_StressVector[2] = aux(1,2);
}


//***********************************************************************************************
//***********************************************************************************************

void PlasticDamage2D::CalculateStressAndTangentMatrix(Vector& StressVector,
                    const Vector& StrainVector,
                    Matrix& algorithmicTangent)
{     
CalculateConstitutiveMatrix(StrainVector, algorithmicTangent);
}
                                       

//***********************************************************************************************
//***********************************************************************************************

void PlasticDamage2D::UpdateMaterial( const Vector& StrainVector,
                                     const Properties& props,
                                     const GeometryType& geom,
                                     const Vector& ShapeFunctionsValues,
                                     const ProcessInfo& CurrentProcessInfo)
{
        noalias(mcurrent_plastic_strain) = mplastic_strain;
        mcurrent_plastic_damage          = mplastic_damage;
        mcurrent_disipation              = mdisipation;
        mcurrent_efective_plastic_strain = mefective_plastic_strain;
        noalias(mcurrent_Ft)             = mFt;

        array_1d<double,4> Variables;
        Variables[0] =  mFt[0];
        Variables[1] =  mFt[1];
        Variables[2] =  mFt[2];
        Variables[3] =  mlength;
   
        mpFluencyCriteria->UpdateVariables(Variables);  
  
}
     

//***********************************************************************************************
//***********************************************************************************************

void PlasticDamage2D::Compute_Principal_Stress(const Vector& StressVector)
{

int    iter = 1000;
double zero = 1.0E-12;
Matrix StressTensor    = ZeroMatrix(3,3);
Matrix EigenVectors    = ZeroMatrix(3,3);
mPrincipalStress       = ZeroVector(3);

mpFluencyCriteria->State_Tensor(StressVector,StressTensor);
SD_MathUtils<double>::EigenVectors(StressTensor, EigenVectors,mPrincipalStress, zero, iter);

mEigenVectors[0] = ZeroVector(3);
mEigenVectors[1] = ZeroVector(3);
mEigenVectors[2] = ZeroVector(3);

mEigenVectors[0][0] = EigenVectors(0,0);  
mEigenVectors[0][1] = EigenVectors(0,1);
mEigenVectors[0][2] = EigenVectors(0,2);

mEigenVectors[1][0] = EigenVectors(1,0);  
mEigenVectors[1][1] = EigenVectors(1,1);
mEigenVectors[1][2] = EigenVectors(1,2);

mEigenVectors[2][0] = EigenVectors(2,0);  
mEigenVectors[2][1] = EigenVectors(2,1);
mEigenVectors[2][2] = EigenVectors(2,2);

}

}

