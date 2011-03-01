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
    
    Von_Mises_Yield_Function::Von_Mises_Yield_Function(myState State, myPotencialPlastic PotencialPlastic )
    :FluencyCriteria()
    {
        mState = State; 
        mPotencialPlastic = PotencialPlastic;
    }
    
    Von_Mises_Yield_Function::~Von_Mises_Yield_Function() {}
    
    //***********************************************************************
    //***********************************************************************
    void Von_Mises_Yield_Function::InitializeMaterial(const Properties& props) 
    { 
        mprops   = &props;
	maccumulated_plastic_strain_current = 0.00;   
        maccumulated_plastic_strain_old     = 0.00;
        mSigma_y                            = (*mprops)[YIELD_STRESS];
	melastic_strain                     = ZeroVector(4);
        mcurrent_elastic_strain             = ZeroVector(4);
        mplastic_strain                     = ZeroVector(4);
        mcurrent_plastic_strain             = ZeroVector(4);
	
    }
    
    void Von_Mises_Yield_Function::ReturnMapping(const Vector& StrainVector, Vector& StressVector)
    {
      
	const double& Young   = (*mprops)[YOUNG_MODULUS];
	const double& Poisson = (*mprops)[POISSON_RATIO];
	const double Gmodu    = Young/(2.00 * (1.00 + Poisson) );
	const double Bulk     = Young/(3.00 * (1.00-2.00*Poisson)); 
	
	///Elastic predictor: Compute elastic trial state
	/// volumetric Strain and pressure stress
	array_1d<double,4> Delta_Strain;
	Delta_Strain[0] = StrainVector[0] - melastic_strain[0];
	Delta_Strain[1] = StrainVector[1] - melastic_strain[1];
	Delta_Strain[2] = StrainVector[2] - melastic_strain[2];
	Delta_Strain[3] = melastic_strain[3];
	
	noalias(mcurrent_elastic_strain) = melastic_strain+ Delta_Strain;
	double volumstrain               = mcurrent_elastic_strain[0] + mcurrent_elastic_strain[1] + mcurrent_elastic_strain[3]; 
	double Pressurstress             = Bulk * volumstrain;
	
	/// desviatoric strain
	array_1d<double,4> Desviatoric_Strain;
	const double  vold3   = volumstrain/3.00; 
	Desviatoric_Strain[0] = mcurrent_elastic_strain[0] - vold3; //tau_xx
	Desviatoric_Strain[1] = mcurrent_elastic_strain[1] - vold3; //tau_yy
	Desviatoric_Strain[2] = 0.50 * mcurrent_elastic_strain[2];  //tau_xy
	Desviatoric_Strain[3] = - vold3;                            //tau_zz
	
	
	/// compute the trial effective stress
	const double Gmodu2  = 2.00 * Gmodu;
	const double Toler   = 1E-6;
        double J2 = Gmodu2 * Gmodu2 * (1.00 *   Desviatoric_Strain[2] * Desviatoric_Strain[2] +
                                       0.50 * ( Desviatoric_Strain[0] * Desviatoric_Strain[0] +
                                                Desviatoric_Strain[1] * Desviatoric_Strain[1] + 
                                                Desviatoric_Strain[3] * Desviatoric_Strain[3]) );


        double Qtrial = std::sqrt(3.00 * J2);
	double Phi    = Qtrial - mSigma_y;
	if( Phi/mSigma_y < Toler)
	{  
	  mcompute_tangent_matrix = false;
	  StressVector[0] = Gmodu2 * Desviatoric_Strain[0] + Pressurstress;
	  StressVector[1] = Gmodu2 * Desviatoric_Strain[1] + Pressurstress;
	  StressVector[2] = Gmodu2 * Desviatoric_Strain[2];  
	}

        else
	{    
          const unsigned int max_iter =  100;
          unsigned int iter           =  1;
	  double ElasticDomain        =  1.00;
	  double lamda                =  0.00;
	  double d                    =  0.00;
	  double H                    = (*mprops)[ISOTROPIC_HARDENING_MODULUS];
	  mcompute_tangent_matrix     = false;
          while(ElasticDomain>Toler && iter<=max_iter)
	  {
	    d      = -3.00 * Gmodu - H;
	    lamda  = lamda - Phi/d;
	    
	    ///Check for convergence
	    Phi                                  =  Qtrial - 3.00 * Gmodu * lamda - mSigma_y;
	    maccumulated_plastic_strain_current  =  maccumulated_plastic_strain_old + lamda;
	  
	    /// Actualiza los valores sigma_y por lamnda
	    Update();
	    
	    ElasticDomain =  std::fabs(Phi/mSigma_y); 
	   
	    
	    if(iter==max_iter){  KRATOS_ERROR(std::logic_error, "WARNING = No Convergence Return Mapping Von Mises", " ") ;}
	  }
	  

	   double factor = 2.00 * Gmodu  * (1.00 - (3.00  * Gmodu * lamda) / Qtrial);
           StressVector[0] = factor * Desviatoric_Strain[0] + Pressurstress;
	   StressVector[1] = factor * Desviatoric_Strain[1] + Pressurstress;
	   StressVector[2] = factor * Desviatoric_Strain[2];
	 
	   
	   factor = (1.00 - (3.00  * Gmodu * lamda) / Qtrial);
	   mcurrent_elastic_strain[0] = factor * Desviatoric_Strain[0] + vold3;
	   mcurrent_elastic_strain[1] = factor * Desviatoric_Strain[1] + vold3;
	   mcurrent_elastic_strain[2] = 2.00 * factor * Desviatoric_Strain[2];
	   mcurrent_elastic_strain[2] = factor * Desviatoric_Strain[2] + vold3;
        }   
      
      
      
    }
    
    
    void Von_Mises_Yield_Function::Update()
    {
      
      mSigma_y  =  (*mprops)[YIELD_STRESS] + (*mprops)[ISOTROPIC_HARDENING_MODULUS] * maccumulated_plastic_strain_current;  
      
    }    
    
    void Von_Mises_Yield_Function::FinalizeSolutionStep()
    {
     
	noalias(melastic_strain)        =  mcurrent_elastic_strain;
	maccumulated_plastic_strain_old =  maccumulated_plastic_strain_current;
      
      
    }     
    
    

}


