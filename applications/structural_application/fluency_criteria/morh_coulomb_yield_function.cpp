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
*   Last Modified by:    $Author: Nelson Lafontaine $
*   Date:                $Date: 26-06-2009 $
*   Revision:            $Revision: 1.0 $
*
* ***********************************************************/

#include "fluency_criteria/morh_coulomb_yield_function.h"

namespace Kratos
  {

           
            Morh_Coulomb_Yield_Function::Morh_Coulomb_Yield_Function(){}
    
            Morh_Coulomb_Yield_Function::Morh_Coulomb_Yield_Function( 
                                                                      const SoftHardPointerType& SofteningBehaviorCohesion,
								      const SoftHardPointerType& SofteningBehaviorFriction,
								      const SoftHardPointerType& SofteningBehaviorDilatancy,
                                                                      const myState& State, 
								      const myPotencialPlastic& PotencialPlastic)
	    :FluencyCriteria()
	    {
              mState                        = State;
	      mpSofteningBehavior_Cohesion  = SofteningBehaviorCohesion;
	      mpSofteningBehavior_Friction  = SofteningBehaviorFriction; 
	      mpSofteningBehavior_Dilatancy = SofteningBehaviorDilatancy;
              mPotencialPlastic             = PotencialPlastic;
	    }

             Morh_Coulomb_Yield_Function::~Morh_Coulomb_Yield_Function() {}

    
//***********************************************************************
//***********************************************************************

void Morh_Coulomb_Yield_Function::InitializeMaterial(const Properties& props) 
{
  
  mprops                                            =  &props;
  mcohesion                                         =  (*mprops)[COHESION]; 
  mdilatancy_angle                                  =  (*mprops)[DILATANCY_ANGLE];
  minternal_friction_angle                          =  (*mprops)[INTERNAL_FRICTION_ANGLE];
  int  size = 4;
  if(mState==Tri_D) size = 6;
  mplastic_strain.resize(size, false); 
  mplastic_strain_old.resize(size, false);
  
  mplastic_step                                     =  false;
  mplastic_strain                                   =  ZeroVector(size);
  mplastic_strain_old                               =  ZeroVector(size); 
  mElastic_strain                                   =  ZeroVector(size);
  mElastic_strain_old                               =  ZeroVector(size);
  mPrincipalPlasticStrain_current                   =  ZeroVector(3);
  mPrincipalPlasticStrain_old                       =  ZeroVector(3);
  mmorh_coulomb_maccumulated_plastic_strain_old     = 0.00; 
  mmorh_coulomb_maccumulated_plastic_strain_current = 0.00; 
  mpastic_damage_old                                = 0.00;
  mpastic_damage_current                            = 0.00;
  mpressure                                         = 0.00;
  melastic_z_old                                    = 0.00;
  melastic_z_new                                    = 0.00;
  msigma_z                                          = 0.00;
  mpSofteningBehavior_Cohesion->InitializeMaterial(props);  
  mpSofteningBehavior_Friction->InitializeMaterial(props);  
  mpSofteningBehavior_Dilatancy->InitializeMaterial(props);
}

//***********************************************************************
//***********************************************************************
 
double Morh_Coulomb_Yield_Function::Toler(const Vector& Stress)
{
  const double tol        = 1.0E-7;
  const double& friction  = mcurrent_minternal_friction_angle; 
  const double& cohe      = mcurrent_cohesion;
  const double cosphi     = std::cos(PI * friction  / 180.00);
  const double sinphi     = std::sin(PI * friction  / 180.00);
  double toler            = tol*((Stress[0] + Stress[2]) * sinphi + 2.00 * cohe* cosphi);
  if(toler>1.00) toler    = 1.00;
  if(toler<0.001)  toler  = 0.001;
  return toler;
}
 
//***********************************************************************
//***********************************************************************


bool Morh_Coulomb_Yield_Function::CheckPlasticAdmisibility(const Vector& Stress)
{ 
  bool result             = false;
  const double tol        = 1.0E-7;
  const double& friction  = mcurrent_minternal_friction_angle; 
  const double& cohe      = mcurrent_cohesion;
  const double sinphi     = std::sin(PI * friction  / 180.00);
  const double cosphi     = std::cos(PI * friction  / 180.00);
  double toler            = tol  * 2.00 * mcurrent_cohesion * cosphi;
  if(toler>0.10) toler    = 0.10;
  if(toler<1E-3) toler    = 1E-3;
  
  // Check plastic admissibility
  double sigma_ef = (Stress[0] - Stress[2]) + (Stress[0] + Stress[2]) * sinphi;
  double phia     = sigma_ef - 2.00 *  cosphi * cohe;
  double res      = phia;
  if(cohe>toler){res =(res/ std::fabs(cohe)); } 
  if(res > toler) result = true;
  return result;
}

//***********************************************************************
//************************************************* ********************** 

void Morh_Coulomb_Yield_Function::ComputeActualStrees(const double& Ppvs, 
						      const array_1d<double,3>& Ppds,
						      const array_1d<double,3>& PrincipalStress,
						      array_1d<double,3>& Sigma)
{
   const double& Young   = (*mprops)[YOUNG_MODULUS];
   const double& Poisson = (*mprops)[POISSON_RATIO];
   const double G        = Young/(2.00 * (1.00 + Poisson) );
   const double K        = Young/(3.00 * (1.00-2.00*Poisson)); 
   array_1d<double,3> I; 
   I[0] = 1.00;
   I[1] = 1.00;
   I[2] = 1.00;
   noalias(Sigma) = PrincipalStress - 2.00 * G * Ppds - K * Ppvs * I;
}    

//***********************************************************************
//*********************************************************************** 

void Morh_Coulomb_Yield_Function::ComputeActualStrain(const array_1d<double,3>& Pps)
{
  Vector PPS_bar(3);
  PPS_bar = ZeroVector(3);
  ComputePlasticStrainBar(mplastic_strain_old ,m_inv_DeltaF, PPS_bar);  
  noalias(mPrincipalPlasticStrain_current) = PPS_bar + Pps;
}


//***********************************************************************
//*********************************************************************** 

void Morh_Coulomb_Yield_Function::ReturnMapping(const Vector& StrainVector, const Vector& TrialStress, Vector& StressVector)
{
    const double& Young   = (*mprops)[YOUNG_MODULUS];
    const double& Poisson = (*mprops)[POISSON_RATIO];
    //const double Gmodu    = Young/(2.00 * (1.00 + Poisson) );
    //const double Bulk     = Young/(3.00 * (1.00-2.00*Poisson)); 
    
    array_1d<double,3> PrincipalStress = ZeroVector(3);
    array_1d<double,3> Sigma = ZeroVector(3);
    array_1d<array_1d < double,3 > ,3> EigenVectors;
    array_1d<unsigned int,3> Order;
    const double d3  = 0.3333333333333333333; 

    const  int dim = TrialStress.size();
    Vector Stress(dim);
    Stress = ZeroVector(dim);
    
    /// computing the trial kirchooff strain
    SpectralDecomposition(TrialStress, PrincipalStress, EigenVectors);
    IdentifyMaximunAndMinumumPrincipalStres_CalculateOrder(PrincipalStress, Order);
    noalias(Sigma)  = PrincipalStress;
    noalias(Stress) = TrialStress;
    
    if(CheckPlasticAdmisibility(PrincipalStress)==true)
    {  
      // plastic Step: Apply return mapping
      // identify possible edge return:either right or left of the main plain  
      const bool edge = ReturnToEdges(PrincipalStress);
      // Apply one-vector return mapping first (return to MAIN PLANE)
      if(ReturnMappingToMainPlane(PrincipalStress, Order,  Sigma)==false){
	 //Apply two-vector return mapping to appropriate EDGES
	 if(TwoVectorReturnToEdges (PrincipalStress, Order, edge, Sigma)==false){
	    //Apply multi-vector return mapping to APEX
	    ReturnMappingToApex(PrincipalStress, Sigma);
	 }
      }   
        
      CalculatePlasticDamage(Sigma);
      AssembleUpdateStressAndStrainTensor(Sigma,  EigenVectors,  Order, StrainVector, Stress);

      
      //CalculatePlasticStrain(StrainVector, StressVector, mplastic_strain, Sigma[2]);

      ///WARNING = Valid for plane Strain 
      Matrix PlasticTensor;  PlasticTensor.resize(3,3,false); 
      PlasticTensor = ZeroMatrix(3,3);
      /// Updating the  elastic plastic strain tensor
      for(unsigned int i=0; i<3; i++)     
          noalias(PlasticTensor) +=  mPrincipalPlasticStrain_current[i] * Matrix(outer_prod(EigenVectors[Order[i]],  EigenVectors[Order[i]]));
      
      mplastic_strain[0] = PlasticTensor(0,0);
      mplastic_strain[1] = PlasticTensor(1,1);
      mplastic_strain[2] = 2.00 * PlasticTensor(1,0);
      mplastic_strain[3] = PlasticTensor(2,2);      
    }   

    CalculateElasticStrain(Stress, mElastic_strain);
    mpressure       =  d3 * (Sigma[0] + Sigma[1] + Sigma[2]);
    //double eevd3    =  d3*mpressure/Bulk;
    //melastic_z_new  = (Sigma[3] - p ) /  (2.00 * Gmodu) + eevd3; 
    
    switch(mState)
      {
      case Plane_Stress:
      {
      KRATOS_ERROR(std::logic_error,  "PLANE STRESS NOT IMPLEMENTED" , "");
      break;
      }
      case Plane_Strain:
      {    
	/// computing the pressure
        StressVector[0] = Stress[0];
        StressVector[1] = Stress[1];
        StressVector[2] = Stress[2];
        msigma_z        = Stress[3];
	melastic_z_new  = mElastic_strain[3]; 
	break;
      }
	
      case Tri_D:
	{
          /// computing the pressure
	  StressVector[0] = Stress[0];
	  StressVector[1] = Stress[1];
	  StressVector[2] = Stress[2];
	  StressVector[3] = Stress[3];
	  StressVector[4] = Stress[4];
	  StressVector[5] = Stress[5];
	  break;
	}
     } 
}

//***********************************************************************
//***********************************************************************

bool Morh_Coulomb_Yield_Function::ReturnMappingToMainPlane(const array_1d<double,3>& PrincipalStress, const array_1d<unsigned int,3>& order, array_1d<double,3>& Sigma)
{
  
UpdateMaterial(); 
const double& E         =   (*mprops)[YOUNG_MODULUS];
const double& NU        =   (*mprops)[POISSON_RATIO];

const double G          = 0.5 * E / (1.00 + NU);
const double K          = E / (3.00 * (1.00-2.00 * NU) );
const double toler      = Toler(PrincipalStress);
const unsigned max      = 10000;

double   dgama    = 0.00;
double   ddgama   = 0.00;
unsigned iter     = 0;
double   denom    = 0.00;

double sinphi             =   std::sin(PI * this->mcurrent_minternal_friction_angle  / 180.00);
double cosphi             =   std::cos(PI * this->mcurrent_minternal_friction_angle  / 180.00);
double sinpsi             =   std::sin(PI * this->mcurrent_dilatancy_angle  / 180.00);
double cospsi             =   std::cos(PI * this->mcurrent_dilatancy_angle  / 180.00); 


// Start Newton-Raphson iterations for DGAMA
double sigma_ef = (PrincipalStress[0] - PrincipalStress[2]) + (PrincipalStress[0] + PrincipalStress[2]) * sinphi;
double phia     = sigma_ef - 2.00 *  cosphi * this->mcurrent_cohesion;
double res      = phia; ///std::fabs(sigma_ef);


double fact              = 0.00; 
double Partial_Cohesion  = 0.00;  
double Partial_Friction  = 0.00;
double Partial_Dilatancy = 0.00;
double Partial_Ep_gama_a = 0.00; 
Vector Imput_Parameters;
Imput_Parameters.resize(4); Imput_Parameters = ZeroVector(4);
Imput_Parameters[0] =  mlength;
Imput_Parameters[1] =  mmorh_coulomb_maccumulated_plastic_strain_current;
Imput_Parameters[2] =  mcurrent_minternal_friction_angle;

Partial_Cohesion    =  mpSofteningBehavior_Cohesion->FirstDerivateFunctionBehavior(Imput_Parameters);  
Partial_Friction    =  mpSofteningBehavior_Friction->FirstDerivateFunctionBehavior(Imput_Parameters);   
Partial_Dilatancy   =  mpSofteningBehavior_Dilatancy->FirstDerivateFunctionBehavior(Imput_Parameters);    

const double d3         = 0.3333333333333333;
const double raiz2d3    = 0.8164958092773; 
double Ppvs             = 0.00;                /// principal plastic volumetric strain  
array_1d<double,3> Ppds = ZeroVector(3);       /// principal plastic desviatoric strain
array_1d<double,3> Pps  = ZeroVector(3);       /// principal plastic  strain
array_1d<double,3> I; 
I[0] = 1.00;
I[1] = 1.00;
I[2] = 1.00;

Vector Imput(3);
Imput    = ZeroVector(3);
//bool mov = false; 
while(fabs(res)>toler && iter++ < max )
   {      
     
     /// IS OK 
    denom   =  (PrincipalStress[0] + PrincipalStress[2]) * cosphi * Partial_Friction * Partial_Ep_gama_a
               - 2.00 * cosphi * Partial_Cohesion * Partial_Ep_gama_a 
               + 2.00 * this->mcurrent_cohesion * sinphi * Partial_Ep_gama_a
               - 4.00 * G - 4.00 * (K + G/3.00)*(dgama * cosphi * sinpsi * Partial_Friction * Partial_Ep_gama_a
               + dgama * cospsi * sinphi * Partial_Dilatancy * Partial_Ep_gama_a)
               - 4.00 * (K +  G/3.00) * sinphi * sinpsi; 
    
    ddgama  =  (-phia)/denom;
    dgama  +=   ddgama;
    
        ///some checks
   if(iter>=max) //|| dgama < 0.00)
    {
     KRATOS_WATCH(iter)
     KRATOS_WATCH(PrincipalStress)
     KRATOS_WATCH(dgama)
     KRATOS_WATCH(toler) 
     KRATOS_WATCH(res) 
     KRATOS_WATCH(Sigma) 
     KRATOS_WATCH(mcohesion) 
     KRATOS_WATCH(mcurrent_cohesion)
     KRATOS_WATCH(mpastic_damage_current)
     KRATOS_WATCH(mpastic_damage_old)
     std::cout<<"RETURN MAPPING TO MAIN PLANE MORH COULOMB  NOT CONVERGED"<< std::endl;
     KRATOS_ERROR(std::logic_error,  "RETURN MAPPING TO MAIN PLANE MORH COULOMB  NOT CONVERGED" , "");
   }
   
     /*
     KRATOS_WATCH(iter)
     KRATOS_WATCH(dgama)
     KRATOS_WATCH(res) 
     KRATOS_WATCH(Sigma) 
     KRATOS_WATCH(mcurrent_cohesion)
     KRATOS_WATCH(mpastic_damage_current)
     KRATOS_WATCH(mcurrent_minternal_friction_angle)
     KRATOS_WATCH(mcurrent_dilatancy_angle)
     KRATOS_WATCH("------------------------")
    */
    
    /// volumetric and desviatoric plastic strain
     Pps[0]        = (1.00 + sinpsi) * dgama;
     Pps[1]        = 0.00;
     Pps[2]        = (sinpsi - 1.00) * dgama;;
     Ppvs          = Pps[0] + Pps[1] + Pps[2]; 
     noalias(Ppds) = Pps - d3 * Ppvs * I;
    
     ///associated accumalated plastic strain
     //fact  = 2.00 * cosphi * dgama; /// debería ser con el angulo de las propiedades
     
    ///von mises acumulated plastic strain
    //fact  = std::fabs(dgama*(std::sqrt((4.00/3.00) * (1.00 + sinpsi*sinpsi)))); 
    fact  = raiz2d3* std::sqrt(inner_prod(Pps,Pps));
    mmorh_coulomb_maccumulated_plastic_strain_current = mmorh_coulomb_maccumulated_plastic_strain_old +  fact;
    ComputeActualStrees(Ppvs, Ppds, PrincipalStress, Sigma);
    ComputeActualStrain(Pps);
    CalculatePlasticDamage(Sigma);
        
    ///variable morh coulomb
    Imput[0] = mpastic_damage_current;
    Imput[1] = mpastic_damage_old;
    Imput[2] = mcohesion; ///old cohesion
    mcurrent_cohesion =  mpSofteningBehavior_Cohesion->EvolucionLaws(Imput, Sigma);
       
    // Compute Internal Variables   
    Imput_Parameters[1]               =  mmorh_coulomb_maccumulated_plastic_strain_current;
    Imput_Parameters[3]               =  mpastic_damage_current;
    //mcurrent_cohesion               =  mpSofteningBehavior_Cohesion->FunctionBehavior(Imput_Parameters); 
    mcurrent_minternal_friction_angle =  mpSofteningBehavior_Friction->FunctionBehavior(Imput_Parameters);   
    Imput_Parameters[2]               =  mcurrent_minternal_friction_angle;
    mcurrent_dilatancy_angle          =  mpSofteningBehavior_Dilatancy->FunctionBehavior(Imput_Parameters);
    
    
     // updating cos and sin
     sinphi             =   std::sin(PI * this->mcurrent_minternal_friction_angle  / 180.00);
     cosphi             =   std::cos(PI * this->mcurrent_minternal_friction_angle  / 180.00);
     sinpsi             =   std::sin(PI * this->mcurrent_dilatancy_angle  / 180.00);
     cospsi             =   std::cos(PI * this->mcurrent_dilatancy_angle  / 180.00); 
    
     
     ///associated accumalated plastic strain 
     //Partial_Ep_gama_a = 2.00 * cosphi; 
     
     /// not  associated accumalated plastic strain
     Partial_Ep_gama_a   = std::sqrt((4.00/3.00)  * (1.00 + sinpsi*sinpsi) );
     Partial_Cohesion    = mpSofteningBehavior_Cohesion->FirstDerivateFunctionBehavior(Imput_Parameters);   
     Partial_Friction    = mpSofteningBehavior_Friction->FirstDerivateFunctionBehavior(Imput_Parameters);   
     Partial_Dilatancy   = mpSofteningBehavior_Dilatancy->FirstDerivateFunctionBehavior(Imput_Parameters); 
    
      /// por si la cohesion tiene un nan
      bool cohe = bool(mcurrent_cohesion==mcurrent_cohesion);
      if(cohe==false){
	 KRATOS_WATCH(mcohesion)
	 KRATOS_WATCH(mcurrent_cohesion)
	 KRATOS_ERROR(std::logic_error,  "ONEEEEEEEEEEEEEEEEEEEEEEEE" , "");
      }
     
     
     if(mcohesion <= toler)
           mcohesion = 0.00;
    
    sigma_ef = (PrincipalStress[0] - PrincipalStress[2]) + (PrincipalStress[0] + PrincipalStress[2]) * sinphi;
    phia     = sigma_ef - 2.00 *  cosphi * this->mcurrent_cohesion - 4.00 * K * dgama * sinphi * sinpsi
               -4.00 * G * dgama * (1.00 + (1.00/3.00) * sinphi * sinpsi ) ;  
     
    // Check convergence
    res = std::fabs(phia);
    //if(sigma_ef!=0) res = res/ std::fabs(sigma_ef);        
   }

  /// updating egain volumetric and desviatoric plastic strain
  Pps[0]        = (1.00 + sinpsi) * dgama;
  Pps[1]        = 0.00;
  Pps[2]        = (sinpsi - 1.00) * dgama;;
  Ppvs          = Pps[0] + Pps[1] + Pps[2]; 
  noalias(Ppds) = Pps - d3 * Ppvs * I;
  
  // Check validity of 1-vector return (check sextant of converged stress)
  //Sigma[0] = PrincipalStress[0] - ( 2.00 * G * (1.00 + sinpsi/3.00) + 2.00 * K * sinpsi) * dgama; 
  //Sigma[1] = PrincipalStress[1] + (4.00 * G /3.00 - 2.00 * K) * sinpsi  * dgama; 
  //Sigma[2] = PrincipalStress[2] + ( 2.00 * G * (1.00 - sinpsi/3.00) - 2.00 * K * sinpsi) * dgama; 
  
  noalias(Sigma) = PrincipalStress - 2.00 * G * Ppds - K * Ppvs * I;  
  
  bool check = CheckValidity(Sigma);
  if(check==true)
  {
  //updating the correct principal pastic strain
  Vector PPS_bar(3);
  PPS_bar = ZeroVector(3);
  ComputePlasticStrainBar(mplastic_strain_old ,m_inv_DeltaF, PPS_bar);  
  noalias(mPrincipalPlasticStrain_current) = PPS_bar + Pps;
  //mPrincipalPlasticStrain_current[0] = /*mPrincipalPlasticStrain_old[0]*/ PPS_bar[0] +  dgama * (1.00  + sinpsi);
  //mPrincipalPlasticStrain_current[1] = /*mPrincipalPlasticStrain_old[1]*/ PPS_bar[1]  + 0.00;
  //mPrincipalPlasticStrain_current[2] = /*mPrincipalPlasticStrain_old[2]*/ PPS_bar[2]  + dgama * (sinpsi - 1.00); 
  }
  
  //KRATOS_WATCH(Sigma)
  return check;
  
}

//***********************************************************************
//***********************************************************************

bool Morh_Coulomb_Yield_Function::TwoVectorReturnToEdges(const array_1d<double,3>& PrincipalStress, const array_1d<unsigned int,3>& order, const bool& edges, array_1d<double,3>& Sigma)
{
 
  UpdateMaterial();
  Sigma = ZeroVector(3);
  const double& E         =   (*mprops)[YOUNG_MODULUS];
  const double& NU        =   (*mprops)[POISSON_RATIO];
  const double G          = 0.5 * E / (1.00 + NU);
  const double K          = E / (3.00 * (1.00-2.00 * NU) );
  const double toler      = Toler(PrincipalStress);
  const unsigned max      = 10000;
  unsigned int iter       = 0;

  
  Matrix d                    = ZeroMatrix(2,2);
  Matrix d_inv                = ZeroMatrix(2,2);
  array_1d<double,2> dgama    = ZeroVector(2);
  array_1d<double,2> ddgama   = ZeroVector(2);
  array_1d<double,2> residual = ZeroVector(2); 
 
  
  double sinphi  =   std::sin(PI * this->mcurrent_minternal_friction_angle  / 180.00);
  double cosphi  =   std::cos(PI * this->mcurrent_minternal_friction_angle  / 180.00);
  double sinpsi  =   std::sin(PI * this->mcurrent_dilatancy_angle  / 180.00);
  double cospsi  =   std::cos(PI * this->mcurrent_dilatancy_angle  / 180.00); 
  

  double sigma_a = 0.00;
  double sigma_b = 0.00;
  double Partial_Cohesion  = 0.00;  
  double Partial_Friction  = 0.00;
  double Partial_Dilatancy = 0.00;
  double Partial_Ep_gama_a = 0.00; 
  double Partial_Ep_gama_b = 0.00; 
  
  const double raiz2d3 = 0.8164958092773; 
  //double aux_1   = 0.00;
  //double aux_2   = 0.00;
  //double aux_3   = 0.00;
  //double aux_4   = 0.00;
  double aux     = 1.00;
  double a       = 0.00;
  double b       = 0.00;

  
  
  double fact1   =  2.00 *  cosphi * this->mcurrent_cohesion;
  double sum     =  0.00; 
  
  array_1d<double, 16> C = ZeroVector(16);
  
  sigma_a     = (PrincipalStress[0] - PrincipalStress[2]) + (PrincipalStress[0] + PrincipalStress[2]) * sinphi; 
  residual[0] = sigma_a - fact1; 
  
  
  if(edges==true)  //rigth edges    
    sigma_b   =  (PrincipalStress[0] - PrincipalStress[1]) + (PrincipalStress[0] + PrincipalStress[1]) * sinphi;
  else  //left edges
    sigma_b   =  (PrincipalStress[1] - PrincipalStress[2]) + (PrincipalStress[1] + PrincipalStress[2]) * sinphi;
    
  
  residual[1]   =  sigma_b - fact1;  
  
  int singular  = 0.00;
  double norma  = norm_2(residual);
  double phipsi = 0.00;
  Vector Imput_Parameters;
  Imput_Parameters.resize(4);
  Imput_Parameters    =  ZeroVector(4); 
  Imput_Parameters[0] =  mlength;
  Imput_Parameters[1] =  mmorh_coulomb_maccumulated_plastic_strain_current;
  Imput_Parameters[2] =  mcurrent_minternal_friction_angle;
  Partial_Cohesion    =  mpSofteningBehavior_Cohesion->FirstDerivateFunctionBehavior(Imput_Parameters);   
  Partial_Friction    =  mpSofteningBehavior_Friction->FirstDerivateFunctionBehavior(Imput_Parameters);   
  Partial_Dilatancy   =  mpSofteningBehavior_Dilatancy->FirstDerivateFunctionBehavior(Imput_Parameters); 

  const double d3         = 0.3333333333333333;
  double Ppvs = 0.00;                          /// principal plastic volumetric strain  
  array_1d<double,3> Ppds = ZeroVector(3);     /// principal plastic desviatoric strain
  array_1d<double,3> Pps  = ZeroVector(3);      /// principal plastic  strain
  array_1d<double,3> I; 
  I[0] = 1.00;
  I[1] = 1.00;
  I[2] = 1.00;
  
  Vector Imput(3);
  Imput = ZeroVector(3);
  
  while(norma > toler && iter++ <max)
{
    if(edges==true)  /// right part
    {
      C[1] = (1.00 - sinpsi) * (1.00 - sinpsi); 
      C[2] = (1.00 + sinpsi) * (1.00 + sinpsi);
      C[3] = (PrincipalStress[0] + PrincipalStress[2] - 2.00 * G * dgama[1]) * cosphi * Partial_Friction; 
      C[4] = (PrincipalStress[0] + PrincipalStress[1] - 2.00 * G * dgama[0]) * cosphi * Partial_Friction; 
      C[5] = 2.00 * G * (1.00 + sinphi + sinpsi);
    }
    else
    {
      
      C[1] = (1.00 + sinpsi) * (1.00 + sinpsi); 
      C[2] = (1.00 - sinpsi) * (1.00 - sinpsi);
      C[3] = (PrincipalStress[0] + PrincipalStress[2] + 2.00 * G * dgama[1]) * cosphi * Partial_Friction; 
      C[4] = (PrincipalStress[0] + PrincipalStress[1] + 2.00 * G * dgama[0]) * cosphi * Partial_Friction; 
      C[5] = 2.00 * G * (1.00 - sinphi - sinpsi); 
    }
  
   
   C[6]  = 2.00 * this->mcurrent_cohesion * sinphi * Partial_Friction - 2.00 * Partial_Cohesion * cosphi;
   C[7]  = -(4.00/3.00) * G * dgama[0] - 4.00 * K * dgama[0] -4.00 * K * dgama[1] + (2.00/3.00) * G * dgama[1]; 
   C[8]  = sinphi*cospsi*Partial_Dilatancy*Partial_Ep_gama_a + sinpsi*cosphi*Partial_Friction*Partial_Ep_gama_a;  
   C[9]  = -(4.00/3.00) * G * dgama[1] - 4.00 * K * dgama[0] -4.00 * K * dgama[1] + (2.00/3.00) * G * dgama[0]; 
   C[10] = sinphi*cospsi*Partial_Dilatancy*Partial_Ep_gama_b + sinpsi*cosphi*Partial_Friction*Partial_Ep_gama_b; 
   C[11] = (4.00/3.00) * G + 4.00 * K;
   C[12] = (2.00/3.00) * G - 4.00 * K;
   C[13] = sinphi * sinpsi;
   C[14] = 2.00 * G * dgama[1] * cospsi * Partial_Dilatancy;
   C[15] = 2.00 * G * dgama[0] * cospsi * Partial_Dilatancy;
  
  
  //Compute residual derivative matrix
  d(0,0) = C[7] * C[8]  - C[11] * C[13] - 4.00 * G +  (C[6] + C[3] -C[14]) * Partial_Ep_gama_a; 
  d(0,1) = C[7] * C[10] + C[12] * C[13] - C[5]     +  (C[6] + C[3] -C[14]) * Partial_Ep_gama_b; 
  d(1,0) = C[9] * C[8]  + C[12] * C[13] - C[5]     +  (C[6] + C[4] -C[15]) * Partial_Ep_gama_a; 
  d(1,1) = C[9] * C[10] - C[11] * C[13] - 4.00 * G +  (C[6] + C[4] -C[15]) * Partial_Ep_gama_b; 
     
  singular       =  SD_MathUtils<double>::InvertMatrix(d, d_inv);
  ddgama         =  -Vector(prod(d_inv, residual));
  
  
  //Compute Newton-Raphson increment and update variables DGAMA and DGAMB
  noalias(dgama) += ddgama;  
  
  
  //some checks 
  if((iter>=max) || (dgama[0]< 0.00 && dgama[1]<0.00) )
{  
   KRATOS_WATCH(dgama) 
   KRATOS_WATCH(PrincipalStress) 
   KRATOS_WATCH(Sigma) 
   KRATOS_WATCH(norma) 
   KRATOS_ERROR(std::logic_error,  "RETURN MAPPING TO MAIN PLANE AND RIGTH O LEFT MORH COULOMB  NOT CONVERGED" , "");
}
  
   // von mises acumulated plastic strain
  if(edges==true)  //rigth edges
   {
    //aux_1   = (dgama[0] + dgama[1]) *  (dgama[0] + dgama[1]);
    //aux_2   = (dgama[0] * dgama[0]) +  (dgama[1] * dgama[1]);
    //aux_3   = (sinpsi+1.00) * (sinpsi+1.00);
    //aux_4   = (sinpsi-1.00) * (sinpsi-1.00);
    //aux     =  std::fabs(raiz2d3 * std::sqrt( aux_1 * aux_3 + aux_2 * aux_4 ));
    
    /// volumetric and desviatoric plastic strain
    Pps[0]        = (dgama[0] + dgama[1])*(1.00 + sinpsi);
    Pps[1]        = dgama[1]*(sinpsi - 1.00);
    Pps[2]        = dgama[0]*(sinpsi - 1.00);
    Ppvs          = Pps[0] + Pps[1] + Pps[2]; 
    noalias(Ppds) = Pps - d3 * Ppvs * I;
    aux = raiz2d3* std::sqrt(inner_prod(Pps,Pps));
    //KRATOS_WATCH(aux)
    mmorh_coulomb_maccumulated_plastic_strain_current = mmorh_coulomb_maccumulated_plastic_strain_old + aux;
   }
    else   //left edges
   {
    //aux_1   = (dgama[0] + dgama[1]) *  (dgama[0] + dgama[1]);
    //aux_2   = (dgama[0] * dgama[0]) +  (dgama[1] * dgama[1]);
    //aux_3   = (sinpsi+1.00) * (sinpsi+1.00);
    //aux_4   = (sinpsi-1.00) * (sinpsi-1.00);
    //aux     =  std::fabs(raiz2d3 * std::sqrt( aux_1 * aux_4 + aux_2 * aux_3 ));

    /// volumetric and desviatoric plastic strain
    Pps[0]        = dgama[0]*(sinpsi + 1.00); 
    Pps[1]        = dgama[1]*(sinpsi + 1.00);
    Pps[2]        = (dgama[0] + dgama[1])*(sinpsi-1.00);
    Ppvs          = Pps[0] + Pps[1] + Pps[2]; 
    noalias(Ppds) = Pps - d3 * Ppvs * I;
    aux           = raiz2d3* std::sqrt(inner_prod(Pps,Pps));
    mmorh_coulomb_maccumulated_plastic_strain_current = mmorh_coulomb_maccumulated_plastic_strain_old + aux;
   }
     
    /// Compute Internal Variables 
    ComputeActualStrees(Ppvs, Ppds, PrincipalStress, Sigma);
    ComputeActualStrain(Pps);
    CalculatePlasticDamage(Sigma);
    
    
    ///variable morh coulomb
    Imput[0] = mpastic_damage_current;
    Imput[1] = mpastic_damage_old;
    Imput[2] = mcohesion; ///old cohesion
    mcurrent_cohesion =  mpSofteningBehavior_Cohesion->EvolucionLaws(Imput, Sigma);
    
    /// Por si la cohesion en un nan
    bool cohe = bool(mcurrent_cohesion==mcurrent_cohesion);
    if(cohe==false){
      KRATOS_WATCH(mcohesion)
      KRATOS_WATCH(mcurrent_cohesion)
      KRATOS_ERROR(std::logic_error,  "TWOOOOOOOOOOOOOOOOOOOOO" , "");
    }
    
    
    // Compute Internal Variables   
    Imput_Parameters[1] =  mmorh_coulomb_maccumulated_plastic_strain_current;
    Imput_Parameters[3] =  mpastic_damage_current;
    //mcurrent_cohesion                 =  mpSofteningBehavior_Cohesion->FunctionBehavior(Imput_Parameters); 
    mcurrent_minternal_friction_angle =  mpSofteningBehavior_Friction->FunctionBehavior(Imput_Parameters);   
    Imput_Parameters[2] =  mcurrent_minternal_friction_angle;
    mcurrent_dilatancy_angle          =  mpSofteningBehavior_Dilatancy->FunctionBehavior(Imput_Parameters);

     // updating cos and sin
     sinphi             =   std::sin(PI * this->mcurrent_minternal_friction_angle  / 180.00);
     cosphi             =   std::cos(PI * this->mcurrent_minternal_friction_angle  / 180.00);
     sinpsi             =   std::sin(PI * this->mcurrent_dilatancy_angle  / 180.00);
     cospsi             =   std::cos(PI * this->mcurrent_dilatancy_angle  / 180.00); 
   

    Partial_Ep_gama_a   = (2.00/3.00) * (C[2] * (dgama[0] + dgama[1] ) + C[1] * dgama[0])/aux;
    Partial_Ep_gama_b   = (2.00/3.00) * (C[2] * (dgama[0] + dgama[1] ) + C[1] * dgama[1])/aux; 
    Partial_Cohesion    = mpSofteningBehavior_Cohesion->FirstDerivateFunctionBehavior(Imput_Parameters);   
    Partial_Friction    = mpSofteningBehavior_Friction->FirstDerivateFunctionBehavior(Imput_Parameters);   
    Partial_Dilatancy   = mpSofteningBehavior_Dilatancy->FirstDerivateFunctionBehavior(Imput_Parameters); 
     
 
  if(mcohesion < 0.0) 
       mcohesion = 0.00;
  

   phipsi = sinphi*sinpsi;
   a      = 4.00 * G * ( 1.00 + (1.00/3.00)*(phipsi) ) + 4.00 * K * phipsi;
  if(edges==true) 
    b     = 2.00 * G * ( 1.00 + sinphi + sinpsi - (1.00/3.00) *  phipsi ) + 4.00 * K * phipsi; 
  else  // left edges
     b =  2.00 * G * ( 1.00 - sinphi - sinpsi - (1.00/3.00) *  phipsi ) + 4.00 * K * phipsi; 


fact1       =  2.00 *  cosphi * this->mcurrent_cohesion;
residual[0] =  sigma_a - a * dgama[0] - b * dgama[1] - fact1;  
residual[1] =  sigma_b - b * dgama[0] - a * dgama[1] - fact1;


//Check convergence
norma = std::fabs(residual[0]) + std::fabs(residual[1]);  
sum   = std::fabs(sigma_a) + std::fabs(sigma_b);
if(sum!=0) {norma = norma / sum;}

}


//const double aux1 = ( 2.00 * G * (1.00 + sinpsi/3.00) + 2.00 * K * sinpsi);
//const double aux2 = ( 4.00 * G /3.00 - 2.00 * K) * sinpsi; 
//const double aux3 = ( 2.00 * G * (1.00 - sinpsi/3.00) - 2.00 * K * sinpsi);

/// updating volumetric and desviatoric plastic strain
if(edges==true)//rigth edges
{
  Pps[0]        = (dgama[0] + dgama[1])*(1.00 + sinpsi);
  Pps[1]        = dgama[1]*(sinpsi - 1.00);
  Pps[2]        = dgama[0]*(sinpsi - 1.00);
  Ppvs          = Pps[0] + Pps[1] + Pps[2]; 
  noalias(Ppds) = Pps - d3 * Ppvs * I;
}
else //left edges
{
  Pps[0]        = dgama[0]*(sinpsi + 1.00); 
  Pps[1]        = dgama[1]*(sinpsi + 1.00);
  Pps[2]        = (dgama[0] + dgama[1])*(sinpsi-1.00);
  Ppvs          = Pps[0] + Pps[1] + Pps[2]; 
  noalias(Ppds) = Pps - d3 * Ppvs * I;
}

noalias(Sigma) = PrincipalStress - 2.00 * G * Ppds - K * Ppvs * I;


/*
if(edges==true)
 { 
  Sigma[0] = PrincipalStress[0] -  aux1 *( dgama[0]   + dgama[1]);
  Sigma[1] = PrincipalStress[1] +  aux2 *  dgama[0]   + aux3 * dgama[1]; 
  Sigma[2] = PrincipalStress[2] +  aux3 *  dgama[0]   + aux2 * dgama[1]; 
 }
// for left edges 
else
 {
  Sigma[0] = PrincipalStress[0] - aux1 *  dgama[0] + aux2 * dgama[1];   
  Sigma[1] = PrincipalStress[1] + aux2 *  dgama[0] - aux1 * dgama[1]; 
  Sigma[2] = PrincipalStress[2] + aux3 * (dgama[0] +        dgama[1]);
 }
*/  
  
  bool check = CheckValidity(Sigma);

  if (check==true)
  { 
    Vector PPS_bar(3);
    PPS_bar = ZeroVector(3);
    ComputePlasticStrainBar(mplastic_strain_old ,m_inv_DeltaF, PPS_bar);
    noalias(mPrincipalPlasticStrain_current) = PPS_bar + Pps;
    
   /* 
   if(edges ==true)  //rigth
   {
    //updating the correct principal pastic strain
    mPrincipalPlasticStrain_current[0] = mPrincipalPlasticStrain_old[0]  PPS_bar[0] +  (dgama[0] + dgama[1])  * (1.00  + sinpsi);
    mPrincipalPlasticStrain_current[1] = mPrincipalPlasticStrain_old[1]  PPS_bar[1] + dgama[1] * (sinpsi - 1.00);
    mPrincipalPlasticStrain_current[2] = mPrincipalPlasticStrain_old[2]  PPS_bar[2]  + dgama[0] * (sinpsi - 1.00); 
   }
   else
   {
    //updating the correct principal pastic strain
    mPrincipalPlasticStrain_current[0] =  mPrincipalPlasticStrain_old[0]  PPS_bar[0]  + dgama[0] * (sinpsi + 1.00); 
    mPrincipalPlasticStrain_current[1] =  mPrincipalPlasticStrain_old[1]  PPS_bar[1]  + dgama[1] * (sinpsi + 1.00);
    mPrincipalPlasticStrain_current[2] =  mPrincipalPlasticStrain_old[2]  PPS_bar[2]  + (dgama[0] + dgama[1])  * (sinpsi - 1.00);
  }
  */
  } 
  return check;

}


///WARNING = For some reason the definition of the desciatoric flow rule in the tesis of owen is erroneuos
/// I take the  The Mohr–Coulomb model. Return mapping to apex given in box 8.7 of book the owen
void Morh_Coulomb_Yield_Function::ReturnMappingToApex(const array_1d<double,3>& PrincipalStress, array_1d<double, 3 >& Sigma)
{
  
UpdateMaterial();
Sigma = ZeroVector(3);
const double& E         =   (*mprops)[YOUNG_MODULUS];
const double& NU        =   (*mprops)[POISSON_RATIO];
const double K          =   E / (3.00 * (1.00-2.00 * NU) );
//const double G          =   0.5 * E / (1.00 + NU);
const double toler      =   Toler(PrincipalStress);
const unsigned max      =   100;
unsigned int iter       =   0;  


 
double Partial_Cohesion  = 0.00;  
double Partial_Friction  = 0.00;
double Partial_Dilatancy = 0.00;
double Partial_Ep_gama_b = 0.00; 

double sinphi               =   std::sin(PI * (mcurrent_minternal_friction_angle)  / 180.00);
double cosphi               =   std::cos(PI * (mcurrent_minternal_friction_angle)  / 180.00);
double cotphi               =   cosphi/sinphi; 
const double& dil_angle     =  (*mprops)[DILATANCY_ANGLE];
const double& int_angle     =  (*mprops)[INTERNAL_FRICTION_ANGLE];
const double sinpsi_o       =   std::sin(PI * (dil_angle) / 180.00);
const double cosphi_o       =   std::cos(PI * (int_angle) / 180.00);
const double alpha          =   cosphi_o/sinpsi_o;

double d           = 0.00;
double r           = 0.00;
double p_trial     = (1.00/3.00) * ( PrincipalStress[0] + PrincipalStress[1] + PrincipalStress[2] );
double p           = 0.00;


//const double dgama_a   = 1.00/(2.00 * G ); 
//const double raiz2d3   = 0.8164958092773;  
double dgama_b         = 0.00;
double ddgama_b        = 0.00;  // volumetric part

array_1d<double, 3> des = ZeroVector(3);
des[0]        =  PrincipalStress[0] - p_trial;
des[1]        =  PrincipalStress[1] - p_trial;
des[2]        =  PrincipalStress[2] - p_trial; 
double fact   =  0.00;
 
array_1d<double,3> Delta_Ep = ZeroVector(3);
 
r  =  p_trial - K *  dgama_b - mcohesion * cotphi; 
Vector Imput_Parameters;
Imput_Parameters.resize(4);
Imput_Parameters    =  ZeroVector(4);
Imput_Parameters[0] =  mlength;
Imput_Parameters[1] =  mmorh_coulomb_maccumulated_plastic_strain_current;
Imput_Parameters[2] =  mcurrent_minternal_friction_angle;
Partial_Cohesion    =  mpSofteningBehavior_Cohesion->FirstDerivateFunctionBehavior(Imput_Parameters);   
Partial_Friction    =  mpSofteningBehavior_Friction->FirstDerivateFunctionBehavior(Imput_Parameters);   
Partial_Dilatancy   =  mpSofteningBehavior_Dilatancy->FirstDerivateFunctionBehavior(Imput_Parameters); 
array_1d<double,3> Pps = ZeroVector(3);       /// principal plastic  strain

Vector Imput(3);
Imput = ZeroVector(3);

while(std::fabs(r) > toler && iter++<max ) 
{
  
   d         = -K - ( Partial_Cohesion * cotphi - ((mcurrent_cohesion)/(sinphi*sinphi)) * (Partial_Friction) ) *Partial_Ep_gama_b;
   ddgama_b  = -r / d;
   dgama_b  += ddgama_b;
   p         = p_trial - K *  dgama_b; 
     
   // Aculutaded von Misses 
   //fact_a =  prod * prod *  dgama_a  *  dgama_a ;
   //fact   =  raiz2d3 * std::sqrt(fact_a + dgama_b*dgama_b/3.00); 
   //Delta_Ep[0] = dgama_a*des[0] + dgama_b/3.00; 
   //Delta_Ep[1] = dgama_a*des[1] + dgama_b/3.00; 
   //Delta_Ep[2] = dgama_a*des[2] + dgama_b/3.00; 
   //fact        = raiz2d3 * std::sqrt(inner_prod(Delta_Ep,Delta_Ep)); 
   
   ///associated hardening 
   fact = alpha *  dgama_b; 
   mmorh_coulomb_maccumulated_plastic_strain_current = mmorh_coulomb_maccumulated_plastic_strain_old +  fact;
   
    /// Compute Internal Variables 
    Sigma[0] = p;
    Sigma[1] = p;
    Sigma[2] = p;
    Pps[0]   = dgama_b/3.00;
    Pps[1]   = dgama_b/3.00;
    Pps[2]   = dgama_b/3.00;
    ComputeActualStrain(Pps);
    CalculatePlasticDamage(Sigma);
    ///variable morh coulomb
    Imput[0] = mpastic_damage_current;
    Imput[1] = mpastic_damage_old;
    Imput[2] = mcohesion; ///old cohesion
    
    mcurrent_cohesion =  mpSofteningBehavior_Cohesion->EvolucionLaws(Imput, Sigma);
    bool cohe = bool(mcurrent_cohesion==mcurrent_cohesion);
    if(cohe==false){
	 KRATOS_WATCH(mcohesion)
	 KRATOS_WATCH(mcurrent_cohesion)
	 KRATOS_ERROR(std::logic_error,  "THREEEEEEEEEEEEEEEEEEEEE" , "");
    }
  
    // Compute Internal Variables   
    Imput_Parameters[1] =  mmorh_coulomb_maccumulated_plastic_strain_current;
    Imput_Parameters[3] =  mpastic_damage_current;
    //mcurrent_cohesion                 =  mpSofteningBehavior_Cohesion->FunctionBehavior(Imput_Parameters); 
    mcurrent_minternal_friction_angle =  mpSofteningBehavior_Friction->FunctionBehavior(Imput_Parameters);   
    Imput_Parameters[2]      =  mcurrent_minternal_friction_angle;
    mcurrent_dilatancy_angle =  mpSofteningBehavior_Dilatancy->FunctionBehavior(Imput_Parameters);
    
   
    // updating cos and sin
    sinphi  =   std::sin(PI * (mcurrent_minternal_friction_angle)  / 180.00);
    cosphi  =   std::cos(PI * (mcurrent_minternal_friction_angle)  / 180.00);
    cotphi  =   cosphi/sinphi; 

    Partial_Ep_gama_b   = alpha; //(2.00/9.00) *  dgama_b / fact;
    Partial_Cohesion    = mpSofteningBehavior_Cohesion->FirstDerivateFunctionBehavior(Imput_Parameters);   
    Partial_Friction    = mpSofteningBehavior_Friction->FirstDerivateFunctionBehavior(Imput_Parameters);  
    Partial_Dilatancy   = mpSofteningBehavior_Dilatancy->FirstDerivateFunctionBehavior(Imput_Parameters); 
     
   p = p_trial - K *  dgama_b;   
   r = p - (mcurrent_cohesion) * cotphi;     

   if(mcurrent_cohesion <=0.00 && std::fabs(r)<toler) 
     {
       mcurrent_cohesion = 0.00;
       Sigma[0] = 1E-9;  
       Sigma[1] = 1E-9;  
       Sigma[2] = 1E-9; 
       break; 
     }
   
}

if(iter>=max && dgama_b<0.00)
{ 
   KRATOS_WATCH(r)
   KRATOS_WATCH(toler)
   KRATOS_WATCH(mcohesion) 
   KRATOS_WATCH(mcurrent_cohesion)
   KRATOS_WATCH(PrincipalStress) 
   std::cout<< "RETURN MAPPING TO APEX  NOT CONVERGED"<< std::endl;
   p = (mcurrent_cohesion) * cotphi;
   KRATOS_ERROR(std::logic_error,  "RETURN MAPPING TO APEX  NOT CONVERGED" , "");
}

   
Sigma[0] = p;
Sigma[1] = p;
Sigma[2] = p;


Vector PPS_bar(3);
PPS_bar = ZeroVector(3);
ComputePlasticStrainBar(mplastic_strain_old ,m_inv_DeltaF, PPS_bar);
//updating the correct principal pastic strain
mPrincipalPlasticStrain_current[0] =  /*mPrincipalPlasticStrain_old[0]*/ PPS_bar[0] + dgama_b/3.00;  //dgama_a * des[0]  + dgama_b/3.00;  
mPrincipalPlasticStrain_current[1] =  /*mPrincipalPlasticStrain_old[1]*/ PPS_bar[1] + dgama_b/3.00;  //dgama_a * des[1]  + dgama_b/3.00;
mPrincipalPlasticStrain_current[2] =  /*mPrincipalPlasticStrain_old[2]*/ PPS_bar[2] + dgama_b/3.00;  //dgama_a * des[2]  + dgama_b/3.00;

}


bool Morh_Coulomb_Yield_Function::CheckValidity(const array_1d<double,3>& Sigma)
{
    bool check   = false;
    array_1d<double,3> Aux_Sigma = ZeroVector(3);
    Aux_Sigma[0] = std::fabs(Sigma[0]);
    Aux_Sigma[1] = std::fabs(Sigma[1]);
    Aux_Sigma[2] = std::fabs(Sigma[2]);
    double delta = (*max_element(Aux_Sigma.begin(), Aux_Sigma.end())) * 1.00E-3; 
    if( (Sigma[0] + delta) >= Sigma[1] && (Sigma[1] + delta ) >= Sigma[2]){ check = true;}
    return check;
}

bool Morh_Coulomb_Yield_Function::ReturnToEdges(const array_1d<double,3>& PrincipalStress)
{
 
    const double& dilatance =  this->mcurrent_dilatancy_angle; //(*mprops)[DILATANCY_ANGLE];   
    const double sinpsi     =   std::sin(PI * dilatance / 180.00);    
    bool  return_rigth      =   false;           // left edges
    
    //array_1d<double, 3 > T  =   ZeroVector(3); 
    //double fact_1           =   3.00 + sinpsi;
    //double fact_2           =   3.00 - sinpsi;
    //T[0]                    =   2.00 / fact_1  - 1.00/fact_2;  
    //T[1]                    =   -1.00 / fact_1 - 1.00/fact_2;
    //T[2]                    =   -1.00 / fact_1 - 2.00/fact_2;
    double scaprd         = (1.00 - sinpsi) * PrincipalStress[0] - 2.00 * PrincipalStress[1] + (1.00 + sinpsi) * PrincipalStress[2]; 
    //double scaprd           = inner_prod(T, PrincipalStress);
    if(scaprd > 0.00) {return_rigth = true;} // return to rigth edges

    return return_rigth;  
}


void Morh_Coulomb_Yield_Function::FinalizeSolutionStep()
{
   mmorh_coulomb_maccumulated_plastic_strain_old =  mmorh_coulomb_maccumulated_plastic_strain_current; 
   mcohesion                                     =  mcurrent_cohesion;
   mdilatancy_angle                              =  mcurrent_dilatancy_angle;
   minternal_friction_angle                      =  mcurrent_minternal_friction_angle;
   noalias(mPrincipalPlasticStrain_old)          =  mPrincipalPlasticStrain_current; 
   noalias(mplastic_strain_old)                  =  mplastic_strain; 
   mpastic_damage_old                            =  mpastic_damage_current;
   melastic_z_old                                =  melastic_z_new;
   noalias(mElastic_strain_old)                  =  mElastic_strain;
}



void Morh_Coulomb_Yield_Function::UpdateMaterial()

{
   mpastic_damage_current                            = mpastic_damage_old;
   mmorh_coulomb_maccumulated_plastic_strain_current = mmorh_coulomb_maccumulated_plastic_strain_old;
   mcurrent_cohesion                                 = mcohesion; 
   mcurrent_dilatancy_angle                          = mdilatancy_angle; 
   mcurrent_minternal_friction_angle                 = minternal_friction_angle;
   noalias(mPrincipalPlasticStrain_current)          = mPrincipalPlasticStrain_old; 
   noalias(mplastic_strain)                          = mplastic_strain_old; 
   melastic_z_new                                    = melastic_z_old;
   noalias(mElastic_strain)                          = mElastic_strain_old;
}


void Morh_Coulomb_Yield_Function::GetValue(const Variable<double>& rVariable, double& Result)
      {
	
	if(rVariable==COHESION){
	  Result = mcurrent_cohesion;
	}
        
	if(rVariable == DILATANCY_ANGLE)
          Result = mcurrent_dilatancy_angle; 
	
	if(rVariable == INTERNAL_FRICTION_ANGLE)
	  Result =  mcurrent_minternal_friction_angle;
	
	if(rVariable == DAMAGE){
	   Result = mpastic_damage_current;
	}
	
	if(rVariable == PRESSURE)
	   Result = mpressure;
      }

void Morh_Coulomb_Yield_Function::GetValue(const Variable<Matrix>& rVariable, Matrix& Result)
 {
      unsigned int size = 0;  
      if(this->mState==Plane_Strain)
        size = 4;
      else
        size = 6;

   if (rVariable==GREEN_LAGRANGE_PLASTIC_STRAIN_TENSOR)
   {
     Result.resize(1, size);
     for(unsigned int i = 0; i< size; i++ )
        Result(0,i) = mplastic_strain(i); 
   }
   
    return;
 }

void Morh_Coulomb_Yield_Function::GetValue(double& Result)
 { 
    mlength = Result; 
    return;
 }

void Morh_Coulomb_Yield_Function::GetValue(Matrix& Result)
{
	m_inv_DeltaF.resize(3,3, false);
	noalias(m_inv_DeltaF) = ZeroMatrix(3,3);
	switch(mState)
         {
	  case Plane_Stress:
	  {
	  KRATOS_ERROR(std::logic_error,  "PLANE STRESS NOT IMPLEMENTED" , "");
	  break;
	  }
          case Plane_Strain:
            {
	      m_inv_DeltaF(0,0)    = Result(0,0);
	      m_inv_DeltaF(0,1)    = Result(0,1);
	      m_inv_DeltaF(1,0)    = Result(1,0);
	      m_inv_DeltaF(1,1)    = Result(1,1);
	      m_inv_DeltaF(2,2)    = 1.00;
	      break;
            }
            
	  case Tri_D:
            {
	      noalias(m_inv_DeltaF) = Result;
	      break;
            }
           
	 }
 }

bool Morh_Coulomb_Yield_Function::PlasticStep(Vector& Stress)
{
    array_1d<double,3> PrincipalStress = ZeroVector(3);
    array_1d<array_1d < double,3 > ,3> EigenVectors;
    array_1d<unsigned int,3> Order;
    
    SpectralDecomposition(Stress, PrincipalStress, EigenVectors);
    IdentifyMaximunAndMinumumPrincipalStres_CalculateOrder(PrincipalStress, Order);
    mplastic_step =  CheckPlasticAdmisibility(PrincipalStress);
    return mplastic_step;
}
 
void Morh_Coulomb_Yield_Function::GetValue(const Variable<Vector>& rVariable, Vector& Result)
    {
      const int size = mplastic_strain.size();
      if(rVariable==ALMANSI_PLASTIC_STRAIN){
	    Result.resize(size);
	    noalias(Result) = mplastic_strain;
          }
      if(rVariable==ALMANSI_ELASTIC_STRAIN){
	    Result.resize(size);
	    noalias(Result) = mElastic_strain;
          }  
          
    }


void Morh_Coulomb_Yield_Function::CalculatePlasticDamage(const array_1d<double,3>& Sigma)
{
          /*
          //double shi             =  0.81649658092773 * (mmorh_coulomb_maccumulated_plastic_strain_current- mmorh_coulomb_maccumulated_plastic_strain_old);
	  array_1d<double, 3> DeltaPlasticStrain = mPrincipalPlasticStrain_current - mPrincipalPlasticStrain_old;
	  double disipation      =  inner_prod(Sigma, DeltaPlasticStrain);  
	  double gf              =  (*mprops)[CRUSHING_ENERGY]/mlength;
	  double kp_punto        =  disipation; 
	  //double d               =  UniaxialTension(Sigma) * shi; /// son igales fijate en el factor 1.22 es raiz 1.5
	  kp_punto               =  kp_punto/gf; 
	  mpastic_damage_current =  mpastic_damage_old + 2.00 * kp_punto / 3.00;  
	  */
	  
          
          const double toler = 1E-10;
	  double teta_a     =  Tensor_Utils<double>::Mc_aully(Sigma);
          double teta_b     =  std::fabs(Sigma[0]) + std::fabs(Sigma[1]) + std::fabs(Sigma[2]);
	  double teta       =  0.00;
          array_1d<double, 3> DeltaPlasticStrain = ZeroVector(3);
          noalias(DeltaPlasticStrain) = mPrincipalPlasticStrain_current - mPrincipalPlasticStrain_old;
	  double disipation  =  inner_prod(Sigma, DeltaPlasticStrain);
	  
	  
	  if (fabs(teta_b) > toler)
          {
	   teta = teta_a/teta_b;
	   // computing Kp_punto
	   double gc_p            = (*mprops)[CRUSHING_ENERGY]/mlength; 
	   double gf_p            = (*mprops)[FRACTURE_ENERGY]/mlength;
           double h               = (teta/gf_p + (1.00-teta)/gc_p);
	   double kp_punto        = h * disipation;    
	   if(disipation > 0.00)
	   {
	     mpastic_damage_current =  mpastic_damage_old + kp_punto;
	     if(mpastic_damage_current > 1.00)
	       mpastic_damage_current = 1.00;
	   }
	  }
 }


 double Morh_Coulomb_Yield_Function::UniaxialTension(const Vector& Stress)
 {
   const double& friction  = (*mprops)[INTERNAL_FRICTION_ANGLE];  
   const double sinphi     = std::sin(PI * friction  / 180.00);
   return (Stress[0] - Stress[2]) + (Stress[0] + Stress[2]) * sinphi;
 }

}




