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
#include "constitutive_laws/brittle_material_2d.h"



namespace Kratos
{

         BrittleMaterial2D::BrittleMaterial2D () 
	: ConstitutiveLaw()
	
	{
	  KRATOS_ERROR(std::logic_error,"Calling the empty constructor.","");
	}

	 BrittleMaterial2D::BrittleMaterial2D(
         FluencyCriteriaPointer FluencyCriteria,
         PropertiesPointer Property) 
	: ConstitutiveLaw()
	{
	      mpFluencyCriteria              = FluencyCriteria;
              mpProperties                   = Property;
	}
	
	
	
	/**
	 *	TO BE TESTED!!!
	 */
	BrittleMaterial2D::~BrittleMaterial2D ()
	{
	}
	
	
	bool BrittleMaterial2D::Has( const Variable<double>& rThisVariable )
	{
		return false;
	}
	
	bool BrittleMaterial2D::Has( const Variable<Vector>& rThisVariable )
	{
		return false;
	}
	
	bool BrittleMaterial2D::Has( const Variable<Matrix>& rThisVariable )
	{
		return false;
	}
	
	double& BrittleMaterial2D::GetValue( const Variable<double>& rThisVariable, double& rValue )
	{
	  if( rThisVariable == DAMAGE)
              mpFluencyCriteria->GetValue(DAMAGE, rValue); 
          else if(rThisVariable == COHESION)
	     mpFluencyCriteria->GetValue(COHESION, rValue); 
          else if(rThisVariable == DILATANCY_ANGLE)
             mpFluencyCriteria->GetValue(DILATANCY_ANGLE, rValue);
	  else if(rThisVariable == INTERNAL_FRICTION_ANGLE)
             mpFluencyCriteria->GetValue(INTERNAL_FRICTION_ANGLE, rValue);
	  else if(rThisVariable==DELTA_TIME)
              rValue = sqrt(mE/mDE);
           
       return rValue; 
	  
	}   
	

	Vector& BrittleMaterial2D::GetValue( const Variable<Vector>& rThisVariable, Vector& rValue )
	{
          return( rValue );
         }
    
    Matrix& BrittleMaterial2D::GetValue( const Variable<Matrix>& rThisVariable, Matrix& rValue )
    {
       
       if (rThisVariable==GREEN_LAGRANGE_PLASTIC_STRAIN_TENSOR)
	  {
	     for(unsigned int i = 0; i< rValue.size2(); i++ )
		     rValue(0,i) = mpFluencyCriteria->mplastic_strain[i]; 
          }
          
        return( rValue );  
       
      
    }

    void BrittleMaterial2D::SetValue( const Variable<double>& rThisVariable, const double& rValue, 
                                const ProcessInfo& rCurrentProcessInfo )
	{
	}
	
    void BrittleMaterial2D::SetValue( const Variable<Vector>& rThisVariable, const Vector& rValue, 
                                const ProcessInfo& rCurrentProcessInfo )
	{
	}
	
    void BrittleMaterial2D::SetValue( const Variable<Matrix>& rThisVariable, const Matrix& rValue, 
                                const ProcessInfo& rCurrentProcessInfo )
	{
	}
    
    void BrittleMaterial2D::Calculate(const Variable<Matrix >& rVariable, Matrix& rResult, 
                   const ProcessInfo& rCurrentProcessInfo)
    {
    }

    void BrittleMaterial2D::Calculate(const Variable<double>& rVariable, 
                                    double& Output, 
                                    const ProcessInfo& rCurrentProcessInfo)
   {
   }

//***********************************************************************************************
//***********************************************************************************************

	void BrittleMaterial2D::InitializeMaterial( const Properties& props,
	const GeometryType& geom,
	const Vector& ShapeFunctionsValues )
	
	{
	  mE                       = (*mpProperties)[YOUNG_MODULUS];
	  mNU                      = (*mpProperties)[POISSON_RATIO];
	  mDE                      = (*mpProperties)[DENSITY];
	  mlength                  = geom.Length();
	  //mfailurefactor           = 0.00;
	  
          double Gc           = (*mpProperties)[CRUSHING_ENERGY]/mlength;
          double length_limit = 2.00*mE*Gc/((*mpProperties)[FC]*(*mpProperties)[FC]);       

          //if (length_limit<mlength) {KRATOS_ERROR(std::logic_error, "Element length greater than permitted" , ""); }          
          mpFluencyCriteria->InitializeMaterial(*mpProperties);
	  mpFluencyCriteria->GetValue(geom.Length());   
	}


//***********************************************************************************************
//***********************************************************************************************
void BrittleMaterial2D::InitializeSolutionStep( const Properties& props,
                    const GeometryType& geom,
                    const Vector& ShapeFunctionsValues ,
                    const ProcessInfo& CurrentProcessInfo)
{
}


//***********************************************************************************************
//***********************************************************************************************

void BrittleMaterial2D::FinalizeSolutionStep( const Properties& props,
		  const GeometryType& geom, 
		  const Vector& ShapeFunctionsValues ,
		  const ProcessInfo& CurrentProcessInfo)
{ 
        mpFluencyCriteria->FinalizeSolutionStep();        
	
	// Calculating the local fail failure
	// double critical_fracture_strain  = 2.00 * (*mpProperties)[FRACTURE_ENERGY]/(mlength * (*mpProperties)[FT]);
	// double inelastic_fracture_strain = 0.00;  
	// mfailurefactor                   = inelastic_fracture_strain/critical_fracture_strain;
	// mfailurefactor                   = mpFluencyCriteria->mpastic_damage_current;
}

//***********************************************************************************************
//***********************************************************************************************


void BrittleMaterial2D::CalculateElasticMatrix(boost::numeric::ublas::bounded_matrix<double,4,4>& C)
{ 

    // plane strain and axial symmetric
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

void BrittleMaterial2D::CalculateElasticStress(array_1d<double,4>& Strain, array_1d<double,4>& Stress)
{
 // plane strain and axial symmetric 
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
void BrittleMaterial2D::CalculateMaterialResponse( const Vector& StrainVector,
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
    if (CalculateTangent==1){CalculateStressAndTangentMatrix(StressVector,StrainVector, AlgorithmicTangent);}
}


void BrittleMaterial2D::CalculateConstitutiveMatrix(const Vector& StrainVector, Matrix& ConstitutiveMatrix)
{

  // plane strain and axial symmetric  
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

void BrittleMaterial2D::CalculateStress(const Vector& StrainVector, 
		  Vector& StressVector)
{
 
	array_1d<double,4> ElasticStrain          = ZeroVector(4);  
	array_1d<double,4> Current_Plastic_Strain = ZeroVector(4);
	array_1d<double,4> ElasticStress          = ZeroVector(4);     
        StressVector                              = ZeroVector(3); 

        // calculating elastic strain
	ElasticStrain[0] =  StrainVector[0];
	ElasticStrain[1] =  StrainVector[1];
	ElasticStrain[2] =  StrainVector[2];
	ElasticStrain[3] =  0.00;
        

        // The elastic strain trial
        noalias(ElasticStrain) -= (mpFluencyCriteria->mplastic_strain);  //Current_Plastic_Strain;

	// calculating elastic stress trial
	CalculateElasticStress(ElasticStrain, ElasticStress); 
	StressVector[0] =  ElasticStress[0];
	StressVector[1] =  ElasticStress[1];
	StressVector[2] =  ElasticStress[2];
	
	
	// return mapping
	mpFluencyCriteria->ReturnMapping(StrainVector, StressVector);

 }     
 
 

void BrittleMaterial2D::CalculateCauchyStresses(
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
		rCauchy_StressVector[2] = aux(0,1);
}


//***********************************************************************************************
//***********************************************************************************************

void BrittleMaterial2D::CalculateStressAndTangentMatrix(Vector& StressVector,
                    const Vector& StrainVector,
                    Matrix& algorithmicTangent)
{     
      CalculateConstitutiveMatrix(StrainVector, algorithmicTangent);
}
                                       

//***********************************************************************************************
//***********************************************************************************************

void BrittleMaterial2D::UpdateMaterial( const Vector& StrainVector,
                                     const Properties& props,
                                     const GeometryType& geom,
                                     const Vector& ShapeFunctionsValues,
                                     const ProcessInfo& CurrentProcessInfo)
{   
        mpFluencyCriteria->UpdateMaterial();  
	
}
     

//***********************************************************************************************
//***********************************************************************************************



}

