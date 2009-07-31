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

// External includes 
#include<cmath>

// Project includes 
#include "includes/define.h"
#include "constitutive_laws/Isotropic_Damage.h"
#include "includes/constitutive_law.h"
#include "utilities/math_utils.h"
#include "custom_utilities/sd_math_utils.h"
#include "custom_utilities/tensor_utils.h"
#include "includes/variables.h"
#include "includes/process_info.h"
#include "includes/properties.h"


namespace Kratos
{
    namespace Isotropic_Damage_Auxiliaries
    {
        boost::numeric::ublas::bounded_matrix<double,2,2> mstemp;
	#ifdef _OPENMP
	#pragma omp threadprivate(mstemp)
	#endif
        boost::numeric::ublas::bounded_matrix<double,2,2> msaux;
	#ifdef _OPENMP
	#pragma omp threadprivate(msaux)
	#endif
	Matrix ConstitutiveMatrixAux(3,3);
        #ifdef _OPENMP
        #pragma omp threadprivate(ConstitutiveMatrixAux)
        #endif
        Vector StrainVectorPerturbation(3);
	#ifdef _OPENMP
        #pragma omp threadprivate(StrainVectorPerturbation)
        #endif
	Vector StressVectorPerturbation(3);
	#ifdef _OPENMP
        #pragma omp threadprivate(StressVectorPerturbation)
        #endif
	Vector StrainVectorPerturbation_aux(3);
	#ifdef _OPENMP
        #pragma omp threadprivate(StrainVectorPerturbation_aux)
        #endif
	Vector StressVectorPerturbation_aux(3);
	#ifdef _OPENMP
        #pragma omp threadprivate(StressVectorPerturbation_aux)
        #endif



     /* // Ver pagina 57 " Metodologia de Evaluacion de Estructuras de Hormigon Armado"								
	Matrix C(0,0);
	#pragma omp threadprivate(C)
	Matrix D(0,0);
	#pragma omp threadprivate(D)
	Matrix V(0,0);
	#pragma omp threadprivate(V)
        Vector d_Sigma(0);
	#pragma omp threadprivate(d_Sigma)
      */
    } 


    using namespace Isotropic_Damage_Auxiliaries;

	/**
	 *	TO BE TESTED!!!
	 */

         Isotropic_Damage::Isotropic_Damage() 
	: ConstitutiveLaw< Node<3> >()
	
	{
	  KRATOS_ERROR(std::logic_error,"Calling the empty constructor for Isotropic Damage","");
	}

	 Isotropic_Damage::Isotropic_Damage(FluencyCriteriaPointer FluencyCriteria, SofteningHardeningCriteriaPointer SofteningBehavior, PropertiesPointer Property) 
	: ConstitutiveLaw< Node<3> >()
	{
	      mFluencyCriteria   = FluencyCriteria;
              mSofteningBehavior = SofteningBehavior;
              mProperties        = Property;
	      KRATOS_WATCH(*mProperties)
	}
	/**
	 *	TO BE TESTED!!!
	 */
	Isotropic_Damage::~Isotropic_Damage()
	{
	}
	
	
	bool Isotropic_Damage::Has( const Variable<double>& rThisVariable )
	{
		return false;
	}
	
	bool Isotropic_Damage::Has( const Variable<Vector>& rThisVariable )
	{
		return false;
	}
	
	bool Isotropic_Damage::Has( const Variable<Matrix>& rThisVariable )
	{
		return false;
	}
	
	double Isotropic_Damage::GetValue( const Variable<double>& rThisVariable )
	{
	  if( rThisVariable == DAMAGE)
	  {			
	  //KRATOS_WATCH(md)
	  return md;

	  }
	  else return 0.00;
	    //KRATOS_ERROR(std::logic_error, "double Variable case not considered", "");
	}
	

	Vector Isotropic_Damage::GetValue( const Variable<Vector>& rThisVariable )
	{
	    KRATOS_ERROR(std::logic_error, "Vector Variable case not considered", "");
	}
	
	Matrix Isotropic_Damage::GetValue( const Variable<Matrix>& rThisVariable )
	{
	    KRATOS_ERROR(std::logic_error, "Matrix Variable case not considered", "");
	}

    void Isotropic_Damage::SetValue( const Variable<double>& rThisVariable, const double& rValue, 
                                const ProcessInfo& rCurrentProcessInfo )
	{
	}
	
    void Isotropic_Damage::SetValue( const Variable<Vector>& rThisVariable, const Vector& rValue, 
                                const ProcessInfo& rCurrentProcessInfo )
	{
	}
	
    void Isotropic_Damage::SetValue( const Variable<Matrix>& rThisVariable, const Matrix& rValue, 
                                const ProcessInfo& rCurrentProcessInfo )
	{
	}
    
    void Isotropic_Damage::Calculate(const Variable<Matrix >& rVariable, Matrix& rResult, 
                   const ProcessInfo& rCurrentProcessInfo)
    {
    }


//***********************************************************************************************
//***********************************************************************************************

	void Isotropic_Damage::InitializeMaterial( const Properties& props,
	const GeometryType& geom,
	const Vector& ShapeFunctionsValues )
	
	{

	// Nota: Estas Variables no seran necesarias almacenarlas por cada elemento.
	// Basta con llamarlas con sus propiedades.
	mFc    = (*mProperties)[FC];
	mFt    = (*mProperties)[FT];
	mEc    = (*mProperties)[CONCRETE_YOUNG_MODULUS_C];
	mEt    = (*mProperties)[CONCRETE_YOUNG_MODULUS_T];
	mNU    = (*mProperties)[POISSON_RATIO];
	mGE    = (*mProperties)[FRACTURE_ENERGY];
	ml     = sqrt(fabs(geom.Area()));   // longitud del elemento
	mr_old = mFt/sqrt(mEc);
	mr_o   = mr_old; 
	//KRATOS_WATCH(geom.Area())
	
	//mFluencyCriteria->InitializeMaterial(props);
	 mFluencyCriteria->InitializeMaterial(*mProperties);	

	}

		

//***********************************************************************************************
//***********************************************************************************************

void Isotropic_Damage::InitializeSolutionStep( const Properties& props,
                    const GeometryType& geom,
                    const Vector& ShapeFunctionsValues ,
                    const ProcessInfo& CurrentProcessInfo)
{
				
}

//***********************************************************************************************
//***********************************************************************************************


	void Isotropic_Damage::CalculateNoDamageElasticMatrix(Matrix& C, const double E, const double NU)
	{ 

		C.resize(3,3, false);
		double c1 = E / (1.00 - NU*NU);
		double c2 = c1 * NU;
		double c3 = 0.5*E / (1.00 + NU);

		C(0,0) = c1;	    C(0,1) = c2;	    C(0,2) = 0.0;
		C(1,0) = c2;	    C(1,1) = c1;	    C(1,2) = 0.0;
		C(2,0) = 0.0;	   C(2,1) = 0.0;	    C(2,2) = c3;
		
	}
	

//***********************************************************************************************
//***********************************************************************************************


      void Isotropic_Damage::CalculateConstitutiveMatrix(const Vector& StrainVector, Matrix& ConstitutiveMatrix)
      {
      Isotropic_Damage::CalculateNoDamageElasticMatrix(ConstitutiveMatrix,mEc,mNU);
      }


//***********************************************************************************************
//***********************************************************************************************


 void Isotropic_Damage::CalculateDamage(const Matrix& ConstitutiveMatrix, const Vector& StressVector, const Vector& StrainVector )
 { 
       
				
	double Tau 	 = 0.00;
	double A   	 = 0.00;
	double r         = 0.00;
	//double elev      = 0.00; 
         

	
	A  = 1.00/((mGE*mEc)/(ml*mFt*mFt)-0.5);
	if (A < 0.00)
	{    
	A  = 0.00;
	std::cout<< "Warning: A is less than zero"<<std::endl;
	}


        mFluencyCriteria->CalculateEquivalentUniaxialStress(StressVector, StrainVector, ConstitutiveMatrix, Tau);
	r     = std::max(mr_old,Tau);  // rold por ro
        md    = mSofteningBehavior->FunctionSofteningHardeningBehavior(A,mr_o,r);
// 	elev  = A*(1.00-r/mr_o);
// 	md     = 1.00 - (mr_o/r)*exp(elev);
// 	if (md < 0.00) 
// 	{
// 	md = fabs(md);
// 	//std::cout<<"Warning: Damage is less than zero"<<std::endl;
// 	}
	
        this->mr_new  = r;   
 }



//***********************************************************************************************
//***********************************************************************************************


void  Isotropic_Damage::FinalizeSolutionStep( const Properties& props,
                                               const GeometryType& geom, 
                                               const Vector& ShapeFunctionsValues ,
                                             const ProcessInfo& CurrentProcessInfo)

			{
				  
				 //KRATOS_WATCH("---------------") 
				 //KRATOS_WATCH(mr_old);   
				 mr_old = mr_new;
                                 //KRATOS_WATCH(mr_old);
			}


//***********************************************************************************************
//***********************************************************************************************

	void Isotropic_Damage::CalculateNoDamageStress(const Vector& StrainVector, Vector& StressVector)
	{
		

		double c1 = mEc / (1.00 - mNU*mNU);
		double c2 = c1 * mNU;
		double c3 = 0.5*mEc / (1.00 + mNU);
	    
		
		  
		StressVector[0] = c1*StrainVector[0] + c2 * (StrainVector[1])	;
		StressVector[1] = c1*StrainVector[1] + c2 * (StrainVector[0])	;
		StressVector[2] = c3*StrainVector[2];
	}
	

//***********************************************************************************************
//***********************************************************************************************
	

 void Isotropic_Damage::CalculateStress( const Vector& StrainVector, Vector& StressVector)
            
	{
	    ConstitutiveMatrixAux = ZeroMatrix(3,3);
	    Isotropic_Damage::CalculateNoDamageStress(StrainVector, StressVector);
	    Isotropic_Damage::CalculateNoDamageElasticMatrix(ConstitutiveMatrixAux,mEc,mNU);
	    Isotropic_Damage::CalculateDamage(ConstitutiveMatrixAux, StressVector,StrainVector);  
	    noalias(StressVector) = (1.00-md)*StressVector;
            									
        }


//***********************************************************************************************
//***********************************************************************************************	

  void Isotropic_Damage::CalculateCauchyStresses(
		Vector& rCauchy_StressVector,
		const Matrix& rF,
		const Vector& rPK2_StressVector,
		const Vector& rGreenLagrangeStrainVector)
    {
		Matrix S = MathUtils<double>::StressVectorToTensor( rPK2_StressVector );

		double J = MathUtils<double>::Det2( rF );

		noalias(mstemp) = prod(rF,S);
		noalias(msaux)  = prod(mstemp,trans(rF));
		msaux *= J;

		if(rCauchy_StressVector.size() != 3)
			rCauchy_StressVector.resize(3);
		
		rCauchy_StressVector[0] = msaux(0,0);
		rCauchy_StressVector[1] = msaux(1,1);
		rCauchy_StressVector[2] = msaux(1,2);
  }

//***********************************************************************************************
//***********************************************************************************************	
	 void Isotropic_Damage::CalculateStressAndTangentMatrix(Vector& StressVector,
                    const Vector& StrainVector,
                    Matrix& algorithmicTangent)
   {
			 // Using perturbation methods
                         long double delta_strain =  0.00;
			 long double factor       =  1E-14;
                         long double max          =  1E-9;
                         double last_damage       =  md;
			 double last_r            =  mr_new;
                         
			 //Vector StrainVectorPerturbation;
			 //Vector StressVectorPerturbation;
			 //Vector StrainVectorPerturbation_aux;
			 //Vector StressVectorPerturbation_aux;
			 
			 StrainVectorPerturbation.resize(3, false);
			 StressVectorPerturbation.resize(3, false);
			 StrainVectorPerturbation_aux.resize(3, false);
			 StressVectorPerturbation_aux.resize(3, false);
			 
			 // diferencia con caso Damage-3D
			 noalias(StrainVectorPerturbation)     = StrainVector;
			 noalias(StrainVectorPerturbation_aux) = StrainVector;  			  
               
                         algorithmicTangent = ZeroMatrix(3,3);
                         

			 for (unsigned int i=0;i<StrainVectorPerturbation.size();i++)
			    {
				  
				 if  (fabs(StrainVector(i))<1E-10)
				    {
				     delta_strain = (*std::min_element(StrainVector.begin(),StrainVector.end()))*factor;
				     if (delta_strain==0.00)
					 {
					    delta_strain = factor;  //perturbacion
					 }
				    } 
				 else
				    {
				      delta_strain = StrainVector(i)*factor;
				    }
					
				   if (delta_strain < max) {delta_strain=max;}
                                   
                                     
				     StrainVectorPerturbation(i)     += delta_strain;
                                     StrainVectorPerturbation_aux(i) -= delta_strain;
                                     
				     Isotropic_Damage::CalculateStress(StrainVectorPerturbation, StressVectorPerturbation);
				     
                                     Isotropic_Damage::CalculateStress(StrainVectorPerturbation_aux, StressVectorPerturbation_aux);
				     
      


                                     // Antiguo procedimiento
			             //StressVectorPerturbation = StressVectorPerturbation-StressVector;
				     //KRATOS_WATCH(StressVectorPerturbation);
				     //KRATOS_WATCH(StressVectorPerturbation_aux);
				     noalias(StressVectorPerturbation) = StressVectorPerturbation-StressVectorPerturbation_aux;
				     //KRATOS_WATCH(StressVectorPerturbation);
                                     //noalias(StressVectorPerturbation) = StressVectorPerturbation/delta_strain;
				     noalias(StressVectorPerturbation) = StressVectorPerturbation/(2.00*delta_strain);
				     //std::cout<<"Ultimo"<< std::endl; 
                                     //KRATOS_WATCH(StrainVectorPerturbation);
				      //KRATOS_WATCH(StressVectorPerturbation);
				    //KRATOS_WATCH(StrainVector);
				    
                                   
				      
				    for (unsigned int j = 0; j<StrainVectorPerturbation.size(); j++)
					 {
					    algorithmicTangent(j,i) = StressVectorPerturbation(j); 
					 }
				  
				    md     = last_damage;
				    mr_new = last_r;
				    StrainVectorPerturbation = StrainVector; 
				    StrainVectorPerturbation_aux = StrainVector; 
                                    
			    }
                             
			   
		    
			      
			
                        
  }

}



