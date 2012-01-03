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

/* *********************************************************   
*          
*   Last Modified by:    $Author: Nelson $
*   Date:                $Date: 2011-11-9 
\*   Revision:           $Revision: 1.2 $
*
* ***********************************************************/


//#include "utilities/math_utils.h"
//#include "custom_utilities/sd_math_utils.h"
#include "custom_utilities/tensor_utils.h"
#include "fluency_criteria/von_mises_yield_function.h"
#include <cmath>

/***********************************************************************
C STATE UPDATE PROCEDURE FOR THE VON MISES ELASTO-PLASTIC MATERIAL MODEL
C WITH NON-LINEAR (PIECEWISE LINEAR) ISOTROPIC HARDENING:
C IMPLICIT ELASTIC PREDICTOR/RETURN MAPPING ALGORITHM (BOXES 7.3-4).
C PLANE STRAIN AND AXISYMMETRIC IMPLEMENTATIONS.
C***********************************************************************
*/


namespace Kratos
{
    typedef boost::numeric::ublas::vector<Vector> Second_Order_Tensor; // dos opciones: un tensor de segundo orden y/o un vector que almacena un vector		  
    typedef boost::numeric::ublas::vector<Second_Order_Tensor> Third_Order_Tensor;
    
    typedef boost::numeric::ublas::vector<boost::numeric::ublas::vector<Matrix> > Fourth_Order_Tensor;
    
    typedef matrix<Second_Order_Tensor> Matrix_Second_Tensor; // Acumulo un tensor de 2 orden en una matri
    
    /*
    Von_Mises_Yield_Function::Von_Mises_Yield_Function(myState State, myPotencialPlastic PotencialPlastic )
    :FluencyCriteria()
    {
        mState = State; 
        mPotencialPlastic = PotencialPlastic;
    }
    */
        Von_Mises_Yield_Function::Von_Mises_Yield_Function(const myState& State, const SoftHardPointerType& Soft)
    :FluencyCriteria()
    {
        mState          = State; 
        mpSigmaBehavior = Soft;
    }
    
    
    Von_Mises_Yield_Function::~Von_Mises_Yield_Function() {}
    
    //***********************************************************************
    //***********************************************************************
    void Von_Mises_Yield_Function::InitializeMaterial(const Properties& props) 
    { 
        mprops                              = &props;
	maccumulated_plastic_strain_current = 0.00;   
        maccumulated_plastic_strain_old     = 0.00;
        melastic_z_old                      = 0.00;
	mhe                                 = 0.00;
        mSigma_yold                         = (*mprops)[YIELD_STRESS];
	mSigma_ynew                         = (*mprops)[YIELD_STRESS];
        mold_plastic_strain                 = ZeroVector(4);
        mcurrent_plastic_strain             = ZeroVector(4);
        mplastic_damage_old                 = 0.00;
	mplastic_damage_new                 = 0.00;
	mpSigmaBehavior->InitializeMaterial(props); 
    }
    
    //***********************************************************************
    //***********************************************************************
    
    void Von_Mises_Yield_Function::UpdateMaterial()
    { 
        mplastic_damage_new                 = mplastic_damage_old;
        melastic_z_new                      = melastic_z_old;
        mSigma_ynew                         = mSigma_yold;
	maccumulated_plastic_strain_current = maccumulated_plastic_strain_old;   
        noalias(mcurrent_plastic_strain)    = mold_plastic_strain;           
    }
    
    //***********************************************************************
    //***********************************************************************
    
    void Von_Mises_Yield_Function::FinalizeSolutionStep()
    {
        mplastic_damage_old             =   mplastic_damage_new;
        melastic_z_old                  =   melastic_z_new;
        mSigma_yold                     =   mSigma_ynew;
        maccumulated_plastic_strain_old =   maccumulated_plastic_strain_current; 
        noalias(mold_plastic_strain)    =   mcurrent_plastic_strain; 
    }   
    
    //***********************************************************************
    //***********************************************************************
    
    void Von_Mises_Yield_Function::ReturnMapping(const Vector& StrainVector, Vector& StressVector)
    {

        const double toler    = 1E-6; 
        const int  dim        = SizePlasticStrain();
	const double& Young   = (*mprops)[YOUNG_MODULUS];
	const double& Poisson = (*mprops)[POISSON_RATIO];
	const double Gmodu    = Young/(2.00 * (1.00 + Poisson) );
	const double Bulk     = Young/(3.00 * (1.00-2.00*Poisson)); 
	Vector Desviatoric_Strain;	
	Vector Desviatoric_Trial_Stress;
	Vector Elastic_Strain;
	Vector Elastic_Stress;
	
	///Elastic predictor: Compute elastic trial state
	/// volumetric Strain and pressure stress
	Elastic_Strain.resize(dim,false);
	Elastic_Strain[0] = StrainVector[0] - mold_plastic_strain[0];
	Elastic_Strain[1] = StrainVector[1] - mold_plastic_strain[1];
	Elastic_Strain[2] = StrainVector[2] - mold_plastic_strain[2];
	Elastic_Strain[3] = melastic_z_old;                 
	
	double volumstrain   = Elastic_Strain[0] + Elastic_Strain[1] + Elastic_Strain[3];
	double Pressurstress = Bulk * volumstrain;
        
	/// Elastic trial deviatoric strain
	Desviatoric_Strain.resize(dim,false);
	const double  vold3   =  volumstrain/3.00; 
	Desviatoric_Strain[0] =  Elastic_Strain[0] - vold3;     //tau_xx
	Desviatoric_Strain[1] =  Elastic_Strain[1] - vold3;     //tau_yy
	Desviatoric_Strain[3] =  melastic_z_old    - vold3;     //tau_zz
	Desviatoric_Strain[2] =  Elastic_Strain[2];  
	///Convert engineering shear component into physical component
	double shi_xy = 0.50 * Desviatoric_Strain[2];      //tau_xy NOT 2txy
		
	///Compute trial effective stress and uniaxial yield stress
        ///Ambas maneras son análogas
	/*
        double J2 = 4.00 * Gmodu * Gmodu * 
        (                //shi_xy*shi_xy
                        + 0.50*(Desviatoric_Strain[0]*Desviatoric_Strain[0]
                        + Desviatoric_Strain[1]*Desviatoric_Strain[1] + 
                        + Desviatoric_Strain[3]*Desviatoric_Strain[3] ) );
        
	 double Qtrial = std::sqrt(3.00*J2);
	 */
	 const double raiz32 = 1.2247448713916;
	 double Qtrial = 0.00;	 
         Desviatoric_Trial_Stress.resize(dim,false);
	 Desviatoric_Trial_Stress[0] = 2.00*Gmodu*Desviatoric_Strain[0];
	 Desviatoric_Trial_Stress[1] = 2.00*Gmodu*Desviatoric_Strain[1];
	 Desviatoric_Trial_Stress[2] = 2.00*Gmodu*shi_xy;
	 Desviatoric_Trial_Stress[3] = 2.00*Gmodu*Desviatoric_Strain[3];
	 
         Qtrial = Desviatoric_Trial_Stress[0]*Desviatoric_Trial_Stress[0] 
                       + Desviatoric_Trial_Stress[1]*Desviatoric_Trial_Stress[1] 
                       + 2.00 * Desviatoric_Trial_Stress[2]*Desviatoric_Trial_Stress[2]
                       + Desviatoric_Trial_Stress[3]*Desviatoric_Trial_Stress[3];
	 Qtrial       =  raiz32 * std::sqrt(Qtrial);
	 

	double Phi   =  Qtrial - mSigma_ynew;
        if(mSigma_ynew<toler)
	{
	  StressVector[0] = 0.00;
	  StressVector[1] = 0.00;
	  StressVector[2] = 0.00;
	  return;
	}
              
	if(Phi<toler)
	{  
	  mcompute_tangent_matrix = false;
	  StressVector[0] = Desviatoric_Trial_Stress[0] + Pressurstress;
	  StressVector[1] = Desviatoric_Trial_Stress[1] + Pressurstress;
	  StressVector[2] = Desviatoric_Trial_Stress[2];  
	  //KRATOS_WATCH("Elastic")
	  //KRATOS_WATCH(mold_plastic_strain)
	  //KRATOS_WATCH(mcurrent_plastic_strain)
	  //KRATOS_WATCH(Pressurstress)
	  //KRATOS_WATCH(Desviatoric_Trial_Stress)
	  //KRATOS_WATCH(StressVector)
	  //KRATOS_WATCH(StrainVector)            
	}
        else
	{   

          const unsigned int max_iter =  1000;
          unsigned int iter           =  1;
	  double ElasticDomain        =  1.00;
	  double lamda                =  0.00;
	  double d                    =  0.00;
	  double H                    =  0.00; //(*mprops)[ISOTROPIC_HARDENING_MODULUS]; ///Linear Hardening Modulus
	  
	  mcompute_tangent_matrix     = true;
	  Vector Imput_Parameters(2);
	  Imput_Parameters[0] = mhe;
	  Imput_Parameters[1] =  maccumulated_plastic_strain_current;
	  
	  H      =  mpSigmaBehavior->FirstDerivateFunctionBehavior(Imput_Parameters);
          while(ElasticDomain>toler && iter<=max_iter)
	  {
	    d      = -3.00 * Gmodu - H;
	    lamda  = lamda - Phi/d;
	    
	    ///Check for convergence
	    maccumulated_plastic_strain_current  =  maccumulated_plastic_strain_old + lamda;
	    Imput_Parameters[1]                  =  maccumulated_plastic_strain_current;
	    H                                    =  mpSigmaBehavior->FirstDerivateFunctionBehavior(Imput_Parameters);
	    mSigma_ynew  = mpSigmaBehavior->FunctionBehavior(Imput_Parameters); //mSigma_yold  + H * lamda;     
	    if(mSigma_ynew<toler)
	    {
	      mSigma_ynew = 0.00; 
	      break;
	    }
	    
	    Phi                                  =  Qtrial - 3.00 * Gmodu * lamda - mSigma_ynew;
	    ElasticDomain                        =  std::fabs(Phi); 
	    if(iter==max_iter)
	    { KRATOS_WATCH(ElasticDomain)
	      KRATOS_ERROR(std::logic_error, "WARNING = No Convergence Return Mapping Von Mises", " ") ;}
	      
	  }
	   
	   Elastic_Stress.resize(dim,false);
	   /*
	   double factor                     = 1.00 - (3.00  * Gmodu * lamda) / Qtrial;
	   factor                            = (2.00*Gmodu) * factor;
	   StressVector[0]                   = factor * Desviatoric_Strain[0] + Pressurstress;
	   StressVector[1]                   = factor * Desviatoric_Strain[1] + Pressurstress;
	   StressVector[2]                   = factor * shi_xy;
	   /// ingeenering strain 
	   Elastic_Strain[0] = factor * Desviatoric_Strain[0] + vold3;
	   Elastic_Strain[1] = factor * Desviatoric_Strain[1] + vold3;
	   Elastic_Strain[2] = 2.00 * factor * shi_xy;
	   Elastic_Strain[3] = factor * Desviatoric_Strain[3] + vold3;
	   std::cout<< " Two = " << Elastic_Strain<< std::endl;
           */
	   if(mSigma_ynew>toler)
	   {
	   double factor                     = 1.00 - (3.00  * Gmodu * lamda) / Qtrial;
	   noalias(Desviatoric_Trial_Stress) = factor * Desviatoric_Trial_Stress;
           StressVector[0]                   = Desviatoric_Trial_Stress[0] + Pressurstress;
	   StressVector[1]                   = Desviatoric_Trial_Stress[1] + Pressurstress;
	   StressVector[2]                   = Desviatoric_Trial_Stress[2];
	   double sigma_z                    = Desviatoric_Trial_Stress[3] + Pressurstress;
	   Elastic_Stress[0]                 = StressVector[0]; 
	   Elastic_Stress[1]                 = StressVector[1]; 
	   Elastic_Stress[2]                 = StressVector[2]; 
	   Elastic_Stress[3]                 = sigma_z;
	   CalculateElasticStrain(Desviatoric_Trial_Stress, Pressurstress, Elastic_Strain);   
	   CalculatePlasticStrain(StrainVector, StressVector, mcurrent_plastic_strain, sigma_z);
	   CalculatePlasticDamage(Elastic_Stress, Desviatoric_Trial_Stress);
	   }
	   else
	   {
	     StressVector[0] = 0.00;
	     StressVector[1] = 0.00;
	     StressVector[2] = 0.00;
	   }
        }   
        

    }
    
    /// compute te engeenering elastic strain
    void Von_Mises_Yield_Function::CalculateElasticStrain(Vector& DesviatoricStress, double Pressure, Vector& Elastic_Strain)
    {
      const double& Young     = (*mprops)[YOUNG_MODULUS];
      const double& Poisson   = (*mprops)[POISSON_RATIO];
      const double Gmodu      = Young/(2.00 * (1.00 + Poisson) );
      const double Bulk       = Young/(3.00 * (1.00-2.00*Poisson)); 
      const double Gd2        = 1.00/(2.00 *Gmodu);
      double vold3            = Pressure/(3.00*Bulk);
 
      Elastic_Strain[0] =  Gd2 * DesviatoricStress[0] + vold3;
      Elastic_Strain[1] =  Gd2 * DesviatoricStress[1] + vold3;
      Elastic_Strain[2] =  Gd2 * DesviatoricStress[2]; 
      Elastic_Strain[3] =  Gd2 * DesviatoricStress[3] + vold3;  
      melastic_z_new    =  Elastic_Strain[3];
    }
    
    void Von_Mises_Yield_Function::GetValue(const Variable<double>& rVariable, double& Result)
    {
   	if(rVariable == DAMAGE)
	{   
	    Result = mplastic_damage_new;
	}
         if(rVariable == YIELD_STRESS){
	    Result = mSigma_ynew;
	 }
    }
    
    void Von_Mises_Yield_Function::CalculatePlasticDamage(const Vector& Elastic_Stress, const Vector& Desviatoric_Trial_Stress)
{
          /*
          const double toler = 1E-6;
	  array_1d<double, 3> Sigma           = ZeroVector(3);
	  array_1d<double,3>  PrincipalStrain = ZeroVector(3);
	  array_1d<array_1d <double,3 > ,3> EigenVectors;
	  SpectralDecomposition(Elastic_Stress, Sigma, EigenVectors); /// computa las tensiones principales
	  double teta_a     =  Tensor_Utils<double>::Mc_aully(Sigma);
          double teta_b     =  std::fabs(Sigma[0]) + std::fabs(Sigma[1]) + std::fabs(Sigma[2]);
	  double teta       =  0.00;
          array_1d<double, 4>  DeltaPlasticStrain = mcurrent_plastic_strain - mold_plastic_strain;
	  double disipation =  inner_prod(Elastic_Stress, DeltaPlasticStrain);
	  
	  
	  if (fabs(teta_b) > 1E-10)
          {
	   teta = teta_a/teta_b;
	   // computing Kp_punto
	   double gc_p  = (*mprops)[CRUSHING_ENERGY]/mhe; 
	   double gf_p  = (*mprops)[FRACTURE_ENERGY]/mhe;
           double f     =  UniaxialTension(Desviatoric_Trial_Stress);
           double h     = (teta/gf_p + (1.00-teta)/gc_p);
	   double kp_punto  = h * disipation;
	   mplastic_damage_new =  mplastic_damage_old + kp_punto;
	  }
	  */
	  
          const double toler = 1E-6;
	  array_1d<double, 3> hkt = ZeroVector(3);
	  array_1d<double, 3> hkc = ZeroVector(3);
	  array_1d<double, 3> Sigma = ZeroVector(3);
	  array_1d<double,3>  PrincipalStrain = ZeroVector(3);
	  array_1d<array_1d <double,3 > ,3> EigenVectors;
	  SpectralDecomposition(Elastic_Stress, Sigma, EigenVectors); /// computa las tensiones principales
          array_1d<double, 4>  DeltaPlasticStrain = mcurrent_plastic_strain - mold_plastic_strain;
	  SpectralDecompositionStrain(DeltaPlasticStrain, PrincipalStrain, EigenVectors);
	  
	  
	  double teta_at          =  Tensor_Utils<double>::Mc_aully(Sigma);
	  double teta_ac          =  Tensor_Utils<double>::Mc_aully(-Sigma);
          double teta_b           =  std::fabs(Sigma[0]) + std::fabs(Sigma[1]) + std::fabs(Sigma[2]);
           
          // computing Kp_punto
	  double gc_p  = (*mprops)[CRUSHING_ENERGY]/mhe; 
	  double gf_p  = (*mprops)[FRACTURE_ENERGY]/mhe;
          double f     =  UniaxialTension(Desviatoric_Trial_Stress);
	  if(f>toler)
	  {
	    //gc_p         = gc_p * teta_ac/f;
	    //gf_p         = gf_p * teta_at/f;
	    
	  for(unsigned int i = 0; i<3; i++){
	    if(gf_p!=0.00)
	        hkt[i] = Tensor_Utils<double>::Mc_aully(Sigma[i])/gf_p;
	    if(gc_p!=0.00)
	        hkc[i] = Tensor_Utils<double>::Mc_aully(-Sigma[i])/gc_p;
	  }
	  
	  double kp_punto = 0.00; 
	  for(unsigned int i = 0; i<3; i++)
	    kp_punto += ( hkt[i] + hkc[i] ) * PrincipalStrain[i];
	  
	  mplastic_damage_new    =  mplastic_damage_old + kp_punto;
         } 
}

 double Von_Mises_Yield_Function::UniaxialTension(const Vector& Desviatoric_Trial_Stress)
 {
        double Qtrial = 0.00; 
	const double raiz32 = 1.2247448713916;
	if(mState==Plane_Strain){
	Qtrial = Desviatoric_Trial_Stress[0]*Desviatoric_Trial_Stress[0] 
	+ Desviatoric_Trial_Stress[1]*Desviatoric_Trial_Stress[1] 
	+ 2.00 * Desviatoric_Trial_Stress[2]*Desviatoric_Trial_Stress[2]
	+ Desviatoric_Trial_Stress[3]*Desviatoric_Trial_Stress[3];
	Qtrial       =  raiz32 * std::sqrt(Qtrial);
	}
	return Qtrial;
 }
 

 
    void Von_Mises_Yield_Function::GetValue(double& Result)
    {
      mhe = Result;
    }
}


