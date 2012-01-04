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
#include "constitutive_laws/plane_strain.h"
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
	  else if(rThisVariable == FT)
             mpFluencyCriteria->GetValue(FT, rValue);
	  else if(rThisVariable==DELTA_TIME)
              rValue = sqrt(mE/mDE);
          else if(rThisVariable==PRESSURE)
               mpFluencyCriteria->GetValue(PRESSURE, rValue);
       return rValue; 
	  
	}   
	

	Vector& BrittleMaterial2D::GetValue( const Variable<Vector>& rThisVariable, Vector& rValue )
	{
          return( rValue );
         }
    
    Matrix& BrittleMaterial2D::GetValue( const Variable<Matrix>& rThisVariable, Matrix& rValue )
    {    
        if(rThisVariable==GREEN_LAGRANGE_PLASTIC_STRAIN_TENSOR){
	   rValue.resize(1,4, false); 
	   mpFluencyCriteria->GetValue(GREEN_LAGRANGE_PLASTIC_STRAIN_TENSOR, rValue);
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
          //double Gc           = (*mpProperties)[CRUSHING_ENERGY]/mlength;
          //double length_limit = 2.00*mE*Gc/((*mpProperties)[FC]*(*mpProperties)[FC]);       
          //if (length_limit<mlength) {KRATOS_ERROR(std::logic_error, "Element length greater than permitted" , ""); }          
          
          m_invF_old.resize(2,2, false);
          mF.resize(2,2, false);
	  mStrain_Old.resize(3, false);
	  mStrain_New.resize(3, false);
	  mStrain_Old        = ZeroVector(3);
	  mStrain_New        = ZeroVector(3);
	  noalias(m_invF_old) = IdentityMatrix(2,2);
	  noalias(mF)         = IdentityMatrix(2,2);
          mpFluencyCriteria->InitializeMaterial(*mpProperties);
	  mpFluencyCriteria->GetValue(mlength);   
	  
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
        noalias(mStrain_Old) = mStrain_New; 
        SD_MathUtils<double>::InvertMatrix(mF,m_invF_old);   /// WARNING = descomeatr
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
      double G            = 0.5*mE / (1.00 + mNU);
      double K            = mE / (3.00 * (1.00-2.00*mNU) );
      double vol_strain   = Strain[0] + Strain[1] + Strain[3];
      double pt           = K * vol_strain;     
      double vol_strain3  = vol_strain * 0.333333333333333333; 
      
      Stress[0] = 2.00 * G *(Strain[0] - vol_strain3 ) + pt; 
      Stress[1] = 2.00 * G *(Strain[1] - vol_strain3 ) + pt;
      Stress[2] = G *(Strain[2]);
      Stress[3] = 2.00 * G *(Strain[3] - vol_strain3 ) + pt;

}

void BrittleMaterial2D::ComputeTrialElasticStress(const Matrix& Delta_F, const Matrix& Inv_Delta_F, array_1d<double,4>&  ElasticTrialStress)
{
    const double& Young   = (*mpProperties)[YOUNG_MODULUS];
    const double& Poisson = (*mpProperties)[POISSON_RATIO];
    const double Gmodu    = Young/(2.00 * (1.00 + Poisson) );
    const double Bulk     = Young/(3.00 * (1.00-2.00*Poisson)); 
 
    Vector Desviatoric_Strain(4);
    Vector Elastic_Strain(4);
    ElasticTrialStress.resize(4);
    
    
    Matrix Elastic_Strain_Tensor(3,3);  
    noalias(Elastic_Strain_Tensor) = ZeroMatrix(3,3);
    Matrix Bar_Elastic_Tensor(3,3);     /// strain de push-foward de n a n+1
    noalias(Bar_Elastic_Tensor) = ZeroMatrix(3,3);
    Matrix Delta_Strain(2,2);
    noalias(Delta_Strain) = ZeroMatrix(2,2);
    
    
    mpFluencyCriteria->GetValue(ALMANSI_ELASTIC_STRAIN, Elastic_Strain);
    Elastic_Strain_Tensor(0,0) = Elastic_Strain[0];
    Elastic_Strain_Tensor(1,1) = Elastic_Strain[1];
    Elastic_Strain_Tensor(2,2) = Elastic_Strain[3];
    Elastic_Strain_Tensor(0,1) = 0.50 * Elastic_Strain[2];
    Elastic_Strain_Tensor(1,0) = 0.50 * Elastic_Strain[2];

    
    Matrix Aux(3,3); 
    Aux      = ZeroMatrix(3,3);
    Aux(0,0) =  Inv_Delta_F(0,0);
    Aux(1,1) =  Inv_Delta_F(1,1);
    Aux(2,2) =  1.00;
    Aux(1,0) =  Inv_Delta_F(1,0);
    Aux(0,1) =  Inv_Delta_F(0,1);
    
    
    noalias(Bar_Elastic_Tensor) = prod(Matrix(prod(trans(Aux), Elastic_Strain_Tensor)), Aux);
    const identity_matrix<double> Unit(2);
    /// Delta almansi strain del paso n al paso n+1
    Matrix C(2,2); 
        
     
    /// large strain
    noalias(Delta_Strain) = 0.50 * ( Unit - Matrix(prod(trans(Inv_Delta_F), Inv_Delta_F) ) );
    
    Elastic_Strain     = ZeroVector(4);
    Elastic_Strain[0]  = Bar_Elastic_Tensor(0,0) + Delta_Strain(0,0);  //xx
    Elastic_Strain[1]  = Bar_Elastic_Tensor(1,1) + Delta_Strain(1,1);  //yy
    Elastic_Strain[2]  = 2.00 *(Bar_Elastic_Tensor(0,1) + Delta_Strain(0,1) ); //xy
    Elastic_Strain[3]  = Bar_Elastic_Tensor(2,2);//zz
    

    double volumstrain   = Elastic_Strain[0] + Elastic_Strain[1] + Elastic_Strain[3];
    double Pressurstress = Bulk * volumstrain;

    /// Elastic trial deviatoric strain
    const double  vold3   =  volumstrain/3.00; 
    Desviatoric_Strain[0] =  Elastic_Strain[0] - vold3;     //tau_xx
    Desviatoric_Strain[1] =  Elastic_Strain[1] - vold3;     //tau_yy
    Desviatoric_Strain[3] =  Elastic_Strain[3] - vold3;     //tau_zz
    Desviatoric_Strain[2] =  Elastic_Strain[2];  
    ///Convert engineering shear component into physical component
    double shi_xy         = 0.50 * Desviatoric_Strain[2];            //tau_xy NOT 2txy
    ElasticTrialStress[0] = 2.00*Gmodu*Desviatoric_Strain[0] + Pressurstress;
    ElasticTrialStress[1] = 2.00*Gmodu*Desviatoric_Strain[1] + Pressurstress;
    ElasticTrialStress[2] = 2.00*Gmodu*shi_xy;
    ElasticTrialStress[3] = 2.00*Gmodu*Desviatoric_Strain[3] + Pressurstress;
}




//***********************************************************************************************
//***********************************************************************************************
/// FOR LARGE DEFORMATION USING PLASTIC DEFORMATION
void BrittleMaterial2D::ComputeTrialElasticStress(const Vector& Strain, 
						  const Matrix& Inv_Delta_F,
						  array_1d<double,4>&  ElasticTrialStress)
{
    //const int  dim        = 4; /// size of plastic strain
    const int  dim2       = 3; 
    const double& Young   = (*mpProperties)[YOUNG_MODULUS];
    const double& Poisson = (*mpProperties)[POISSON_RATIO];
    const double Gmodu    = Young/(2.00 * (1.00 + Poisson) );
    const double Bulk     = Young/(3.00 * (1.00-2.00*Poisson)); 
    Vector Desviatoric_Strain(4);
    Vector Elastic_Strain(4);
    Vector Plastic_Strain(4);
    ElasticTrialStress.resize(4);
    ///Elastic predictor: Compute elastic trial state
    /// volumetric Strain and pressure stress
    //Elastic_Strain.resize(dim,false);
    //Elastic_Stress.resize(dim,false);
    //Desviatoric_Strain.resize(dim,false);
    
    Matrix Plastic_Strain_Tensor;  Plastic_Strain_Tensor.resize(dim2,dim2, false);  noalias(Plastic_Strain_Tensor)    = ZeroMatrix(dim2,dim2);
    Matrix Bar_Plastic_Tensor;     Bar_Plastic_Tensor.resize(dim2,dim2, false);     noalias(Bar_Plastic_Tensor)       = ZeroMatrix(dim2,dim2);
    
    mpFluencyCriteria->GetValue(ALMANSI_PLASTIC_STRAIN, Plastic_Strain);
    
    Plastic_Strain_Tensor.resize(3,3, false);
    Plastic_Strain_Tensor      = ZeroMatrix(3,3);
    Plastic_Strain_Tensor(0,0) = Plastic_Strain[0];
    Plastic_Strain_Tensor(1,1) = Plastic_Strain[1];
    Plastic_Strain_Tensor(2,2) = Plastic_Strain[3];
    Plastic_Strain_Tensor(0,1) = 0.50 * Plastic_Strain[2];
    Plastic_Strain_Tensor(1,0) = 0.50 * Plastic_Strain[2];
    Bar_Plastic_Tensor.resize(3,3, false);     noalias(Bar_Plastic_Tensor) = ZeroMatrix(3,3); 
    
    Matrix Aux(3,3); 
    Aux.resize(3,3, false); Aux = ZeroMatrix(3,3);
    Aux(0,0) =  Inv_Delta_F(0,0);
    Aux(1,1) =  Inv_Delta_F(1,1);
    Aux(2,2) =  1.00;
    Aux(1,0) =  Inv_Delta_F(1,0);
    Aux(0,1) =  Inv_Delta_F(0,1);
    
    
    noalias(Bar_Plastic_Tensor) = prod(Matrix(prod(trans(Aux), Plastic_Strain_Tensor)), Aux); 

    Plastic_Strain          = ZeroVector(4);
    Plastic_Strain[0]       = Bar_Plastic_Tensor(0,0);  //xx
    Plastic_Strain[1]       = Bar_Plastic_Tensor(1,1);  //yy
    Plastic_Strain[2]       = 2.00 * Bar_Plastic_Tensor(0,1); //xy
    Plastic_Strain[3]       = Bar_Plastic_Tensor(2,2);//zz
    
    Elastic_Strain[0]       = Strain[0] - Plastic_Strain[0];
    Elastic_Strain[1]       = Strain[1] - Plastic_Strain[1];
    Elastic_Strain[2]       = Strain[2] - Plastic_Strain[2];
    Elastic_Strain[3]       = -Plastic_Strain[3];                

    double volumstrain   = Elastic_Strain[0] + Elastic_Strain[1] + Elastic_Strain[3];
    double Pressurstress = Bulk * volumstrain;

    /// Elastic trial deviatoric strain
    const double  vold3   =  volumstrain/3.00; 
    Desviatoric_Strain[0] =  Elastic_Strain[0] - vold3;     //tau_xx
    Desviatoric_Strain[1] =  Elastic_Strain[1] - vold3;     //tau_yy
    Desviatoric_Strain[3] =  Elastic_Strain[3] - vold3;     //tau_zz
    Desviatoric_Strain[2] =  Elastic_Strain[2];  
    ///Convert engineering shear component into physical component
    double shi_xy         = 0.50 * Desviatoric_Strain[2];            //tau_xy NOT 2txy
    ElasticTrialStress[0] = 2.00*Gmodu*Desviatoric_Strain[0] + Pressurstress;
    ElasticTrialStress[1] = 2.00*Gmodu*Desviatoric_Strain[1] + Pressurstress;
    ElasticTrialStress[2] = 2.00*Gmodu*shi_xy;
    ElasticTrialStress[3] = 2.00*Gmodu*Desviatoric_Strain[3] + Pressurstress;
}

//***********************************************************************************************
//***********************************************************************************************

/// FOR SMALL DEFORMATION
void BrittleMaterial2D::ComputeTrialElasticStress(const Vector& Strain, array_1d<double,4>&  ElasticTrialStress)
{
    //const int  dim        = 4; /// size of plastic strain
    //const int  dim2       = 3; 
    const double& Young   = (*mpProperties)[YOUNG_MODULUS];
    const double& Poisson = (*mpProperties)[POISSON_RATIO];
    const double Gmodu    = Young/(2.00 * (1.00 + Poisson) );
    const double Bulk     = Young/(3.00 * (1.00-2.00*Poisson)); 
    Vector Desviatoric_Strain = ZeroVector(4);
    Vector Elastic_Strain     = ZeroVector(4);
    Vector Plastic_Strain     = ZeroVector(4);
    ElasticTrialStress.resize(4,false);
    
    mpFluencyCriteria->GetValue(ALMANSI_ELASTIC_STRAIN, Elastic_Strain);
    
    Elastic_Strain[0] += Strain[0] - mStrain_Old[0];
    Elastic_Strain[1] += Strain[1] - mStrain_Old[1];
    Elastic_Strain[2] += Strain[2] - mStrain_Old[2];                
    
    double volumstrain   = Elastic_Strain[0] + Elastic_Strain[1] + Elastic_Strain[3];
    double Pressurstress = Bulk * volumstrain;

    /// Elastic trial deviatoric strain
    const double  vold3   =  volumstrain/3.00; 
    Desviatoric_Strain[0] =  Elastic_Strain[0] - vold3;     //tau_xx
    Desviatoric_Strain[1] =  Elastic_Strain[1] - vold3;     //tau_yy
    Desviatoric_Strain[3] =  Elastic_Strain[3] - vold3;     //tau_zz
    Desviatoric_Strain[2] =  Elastic_Strain[2];  
    ///Convert engineering shear component into physical component
    double shi_xy         = 0.50 * Desviatoric_Strain[2];            //tau_xy NOT 2txy
    ElasticTrialStress[0] = 2.00*Gmodu*Desviatoric_Strain[0] + Pressurstress;
    ElasticTrialStress[1] = 2.00*Gmodu*Desviatoric_Strain[1] + Pressurstress;
    ElasticTrialStress[2] = 2.00*Gmodu*shi_xy;
    ElasticTrialStress[3] = 2.00*Gmodu*Desviatoric_Strain[3] + Pressurstress;
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
        bool SaveInternalVariables)
{ 
    noalias(mF)          =  ZeroMatrix(2,2); 
    noalias(mStrain_New) =  ZeroVector(3);
    noalias(mF)          =  DeformationGradient;
    noalias(mStrain_New) =  StrainVector;
    UpdateMaterial(StrainVector, props, geom,ShapeFunctionsValues, CurrentProcessInfo);
    if (CalculateTangent==1){CalculateTangentMatrix(DeformationGradient,StrainVector,StressVector,AlgorithmicTangent);}
    if (CalculateStresses==true) CalculateStress(DeformationGradient, StrainVector, StressVector); 
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

void BrittleMaterial2D::CalculateStress(
                  const Matrix& DeformationGradient,
                  const Vector& StrainVector, 
		  Vector& StressVector)
{
 
        StressVector.resize(3,false);
	noalias(StressVector)     =  ZeroVector(3);
        Vector Almansi            = ZeroVector(3); 
	Matrix Delta_F;     Delta_F.resize(2,2, false);     noalias(Delta_F)     = ZeroMatrix(2,2);
	Matrix Inv_Delta_F; Inv_Delta_F.resize(2,2, false); noalias(Inv_Delta_F) = ZeroMatrix(2,2);	
        CalculateAlmansiStrain(StrainVector, DeformationGradient, Almansi); /// push-foward greeen
	noalias(Delta_F)            = prod(DeformationGradient,m_invF_old);
	SD_MathUtils<double>::InvertMatrix(Delta_F, Inv_Delta_F);
	mpFluencyCriteria->GetValue(Inv_Delta_F);
	

	/// return mapping with kirchhoff Stress
	array_1d<double,4>  ElasticTrialStress; 
	
	///WARNING = For large deformation Using plastic
	//ComputeTrialElasticStress(Almansi, inv_Delta_F, ElasticTrialStress);
	
	///WARNING = For small deformation
	ComputeTrialElasticStress(Almansi, ElasticTrialStress);
	
	///WARNING = For large deformation Using Elastic
 	//ComputeTrialElasticStress(Delta_F, Inv_Delta_F, ElasticTrialStress);

	mpFluencyCriteria->ReturnMapping(Almansi, ElasticTrialStress,  StressVector);
	
	/// Computando Second Piola-Kirchhoff
	Vector SPK             = ZeroVector(3);
	double J               = MathUtils<double>::Det2(DeformationGradient);
	noalias(StressVector)  =(1.00/J)*StressVector; 
	CalculateSPK(StressVector, DeformationGradient, SPK);
	noalias(StressVector) = SPK;
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


   void BrittleMaterial2D::CalculateAlmansiStrain(const Vector& Strain, const Matrix& F,  Vector& Almansi)
   {
     
     int size1 = Strain.size();
     int size2 = F.size1();
     
     Matrix F_inv;
     Matrix Strain_Tensor;
     Matrix Almansi_Tensor;
     
     Almansi.resize(size1);
     F_inv.resize(size2, size2);
     Almansi_Tensor.resize(size2, size2);
     Strain_Tensor.resize(size2, size2);
     
     Strain_Tensor = SD_MathUtils<double>::StrainVectorToTensor(Strain);
     SD_MathUtils<double>::InvertMatrix(F, F_inv);
     noalias(Almansi_Tensor) =  prod( Matrix(prod(trans(F_inv), Strain_Tensor)),F_inv);
     
     noalias(Almansi) = SD_MathUtils<double>::TensorToStrainVector(Almansi_Tensor);
     
   }
   
   void BrittleMaterial2D::CalculateSPK(const Vector& Cauchy_Stress, const Matrix& F, Vector& Stress)
   {
     int size1       = Cauchy_Stress.size();
     int size2       = F.size1();
     const double J  = MathUtils<double>::Det2( F );
     
     Matrix F_inv;
     Matrix Stress_Tensor;
     Matrix Cauchy_Tensor;
     
     Stress.resize(size1);
     F_inv.resize(size2, size2);
     Cauchy_Tensor.resize(size2, size2);
     Stress_Tensor.resize(size2, size2);
     
     Cauchy_Tensor = MathUtils<double>::StressVectorToTensor(Cauchy_Stress);
     SD_MathUtils<double>::InvertMatrix(F, F_inv);
     noalias(Stress_Tensor) =  prod( Matrix(prod(F_inv, Cauchy_Tensor)),trans(F_inv));
     SD_MathUtils<double>::TensorToVector(Stress_Tensor, Stress);
     Stress*= J;  
     
   }

//***********************************************************************************************
//***********************************************************************************************

void BrittleMaterial2D::CalculateTangentMatrix(const Matrix& DeformationGradient, const Vector& StrainVector, Vector& StressVector, Matrix& AlgorithmicTangent)
{     
   
        
	//bool tangent = mpFluencyCriteria->PlasticStep(StressVector);
	//if(tangent==false)
	CalculateConstitutiveMatrix(StrainVector, AlgorithmicTangent);
	/*
	else
	{
         long double delta_strain =  0.00;
         long double factor       =  1E-6;
         long double max          =  1E-8;

	
         Matrix ConstitutiveMatrixAux( 3, 3 );
         Vector StrainVectorPerturbation( 3 );
         Vector StressVectorPerturbation( 3 );
         Vector StrainVectorPerturbation_aux( 3 );
         Vector StressVectorPerturbation_aux( 3 );
         noalias( StrainVectorPerturbation )     = StrainVector;
         noalias( StrainVectorPerturbation_aux ) = StrainVector;
	 AlgorithmicTangent.resize(3,3, false);
         AlgorithmicTangent = ZeroMatrix( 3, 3 );
	
         for ( unsigned int i = 0;i < StrainVectorPerturbation.size();i++ )
          {
            if ( fabs( StrainVector( i ) ) < max )
                delta_strain = ( *std::min_element( StrainVector.begin(), StrainVector.end() ) ) * factor;
                if ( delta_strain == 0.00 )
                    delta_strain = factor;  //perturbacion             
            else
                delta_strain = StrainVector( i ) * factor;
           
            if ( delta_strain < max )  delta_strain = max;

            StrainVectorPerturbation( i )     += delta_strain;
            StrainVectorPerturbation_aux( i ) -= delta_strain;
            CalculateStress(DeformationGradient, StrainVectorPerturbation, StressVectorPerturbation );
            CalculateStress(DeformationGradient, StrainVectorPerturbation_aux, StressVectorPerturbation_aux );
            noalias( StressVectorPerturbation ) = StressVectorPerturbation - StressVectorPerturbation_aux;
            noalias( StressVectorPerturbation ) = StressVectorPerturbation / ( 2.00 * delta_strain );


            for ( unsigned int j = 0; j < StrainVectorPerturbation.size(); j++ )
            {
                AlgorithmicTangent( j, i ) = StressVectorPerturbation( j );
            }

            StrainVectorPerturbation     = StrainVector;
            StrainVectorPerturbation_aux = StrainVector;
          }
          //std::cout<< "Tangent " << std::endl;
          //KRATOS_WATCH(AlgorithmicTangent)
          //KRATOS_WATCH("---------------------------------------")
	}
	*/
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


void BrittleMaterial2D::ResetMaterial(const Properties& props,
                 const GeometryType& geom,
                 const Vector& ShapeFunctionsValues)
     {
       
         //const ConstitutiveLaw::Pointer rLaw = ConstitutiveLaw::Pointer(new PlaneStrain());
         //props.SetValue(CONSTITUTIVE_LAW, rLaw);
         //const ConstitutiveLaw::Pointer rLaw = ConstitutiveLaw::Pointer(new PlaneStrain());
         //boost::shared_ptr<ConstitutiveLaw> p_clone(new PlaneStrain());
	 //props[CONSTITUTIVE_LAW] = rLaw;
	 //this = rLaw->Clone();
	 //this = props[CONSTITUTIVE_LAW]->PlaneStrain::Clone();
         //mpFluencyCriteria->InitializeMaterial(props);
     }

//***********************************************************************************************
//***********************************************************************************************

std::size_t BrittleMaterial2D::GetStrainSize()
{
  return 3;
}

//***********************************************************************************************
//***********************************************************************************************

int BrittleMaterial2D::Check(const Properties& props,
                const GeometryType& geom,
                const ProcessInfo& CurrentProcessInfo)
        {
            KRATOS_TRY
          
            
            if(YOUNG_MODULUS.Key() == 0 || props[YOUNG_MODULUS]<= 0.00)
                KRATOS_ERROR(std::invalid_argument,"YOUNG_MODULUS has Key zero or invalid value ","");

	    const double& nu = props[POISSON_RATIO];
	    const bool check = bool( (nu >0.499 && nu<0.501 ) || (nu < -0.999 && nu > -1.01 ) );
	    if(POISSON_RATIO.Key() == 0 || check==true) // props[POISSON_RATIO] == 1.00 || props[POISSON_RATIO] == -1.00)
                KRATOS_ERROR(std::invalid_argument,"POISSON_RATIO has Key zero invalid value ","");

	    if(FRACTURE_ENERGY.Key() == 0 || props[FRACTURE_ENERGY]<= 0.00)
                KRATOS_ERROR(std::invalid_argument,"FRACTURE_ENERGY has Key zero or invalid value ","");
	    
	    if(CRUSHING_ENERGY.Key() == 0 || props[CRUSHING_ENERGY]<= 0.00)
                KRATOS_ERROR(std::invalid_argument,"CRUSHING_ENERGY has Key zero or invalid value ","");
	  
	    if(DENSITY.Key() == 0 || props[DENSITY]<0.00)
                KRATOS_ERROR(std::invalid_argument,"DENSITY has Key zero or invalid value ","");
            
	    if(INTERNAL_FRICTION_ANGLE.Key() == 0 || props[INTERNAL_FRICTION_ANGLE]<0.00)
                KRATOS_ERROR(std::invalid_argument,"INTERNAL_FRICTION_ANGLE has Key zero or invalid value ","");
	    
	    if(DILATANCY_ANGLE.Key() == 0 || props[DILATANCY_ANGLE]<0.00)
                KRATOS_ERROR(std::invalid_argument,"DILATANCY_ANGLE has Key zero or invalid value ","");
	
	    if(COHESION.Key() == 0 || props[COHESION]<0.00)
                KRATOS_ERROR(std::invalid_argument,"COHESION has Key zero or invalid value ","");
	
	    if(FT.Key() == 0 || props[FT]<0.00)
                KRATOS_ERROR(std::invalid_argument,"FC has Key zero or invalid value ","");
	
	    if(FC.Key() == 0 || props[FC]<0.00)
                KRATOS_ERROR(std::invalid_argument,"FC has Key zero or invalid value ","");
    	    
	    return 0;
	    
            KRATOS_CATCH("");
        }


}

