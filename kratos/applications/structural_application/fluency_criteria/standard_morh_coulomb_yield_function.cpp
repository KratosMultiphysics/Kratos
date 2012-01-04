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

#include "fluency_criteria/standard_morh_coulomb_yield_function.h"

namespace Kratos
  {

            /// WARNING = the file is a test. Only appy for perfect plastic.
           
            Standard_Morh_Coulomb_Yield_Function::Standard_Morh_Coulomb_Yield_Function(){}
    
            Standard_Morh_Coulomb_Yield_Function::Standard_Morh_Coulomb_Yield_Function( 
		const SoftHardPointerType& SofteningBehaviorCohesion,
		const myState& State)
	    :FluencyCriteria()
	    {
              mState                        = State;
	      mpSofteningBehavior_Cohesion  = SofteningBehaviorCohesion;
	    }

             Standard_Morh_Coulomb_Yield_Function::~Standard_Morh_Coulomb_Yield_Function() {}

    
//***********************************************************************
//***********************************************************************

void Standard_Morh_Coulomb_Yield_Function::InitializeMaterial(const Properties& props) 
{
  
  mprops                                            =  &props;
  mcohesion                                         =  (*mprops)[COHESION]; 
  int  size = 4;
  if(mState==Tri_D) size = 6;  
  
  mElastic_strain                                   =  ZeroVector(size);
  mElastic_strain_old                               =  ZeroVector(size);
  mPlastic_strain                                   =  ZeroVector(size);
  mPlastic_strain_old                               =  ZeroVector(size);
  mPrincipalPlasticStrain_current                   =  ZeroVector(3);
  mPrincipalPlasticStrain_old                       =  ZeroVector(3);
  mmorh_coulomb_maccumulated_plastic_strain_old     = 0.00; 
  mmorh_coulomb_maccumulated_plastic_strain_current = 0.00; 
  mpSofteningBehavior_Cohesion->InitializeMaterial(props);  
}

//***********************************************************************
//***********************************************************************
 
bool Standard_Morh_Coulomb_Yield_Function::CheckPlasticAdmisibility(const Vector& Stress)
{
  bool result             = false;
  const double tol        = 1.0E-10; 
  const double d180       = 0.0055555555555555555555;
  const double& friction  = (*mprops)[INTERNAL_FRICTION_ANGLE]; 
  const double& cohe      = mcurrent_cohesion;
  const double sinphi     = std::sin(d180 * PI * friction);
  const double cosphi     = std::cos(d180 * PI * friction);

  // Check plastic admissibility
  double sigma_ef = (Stress[0] - Stress[2]) + (Stress[0] + Stress[2]) * sinphi;
  double phia     = sigma_ef - 2.00 *  cosphi * mcohesion;
  double res      = phia;
  if(cohe>tol) { res = res/std::fabs(cohe); } ///WARNING
  if(res > tol) result = true;
  return result; 
}

//***********************************************************************
//*********************************************************************** 
  
void Standard_Morh_Coulomb_Yield_Function::ReturnMapping(const Vector& StrainVector, const Vector& TrialStress, Vector& StressVector)
{
    
    //const double& Young   = (*mprops)[YOUNG_MODULUS];
    //const double& Poisson = (*mprops)[POISSON_RATIO];
    //const double Gmodu    = Young/(2.00 * (1.00 + Poisson) );
    //const double Bulk     = Young/(3.00 * (1.00-2.00*Poisson)); 
    
    array_1d<double,3> PrincipalStress;
    array_1d<double,3> Sigma;
    array_1d<array_1d < double,3 > ,3> EigenVectors;
    array_1d<unsigned int,3> Order;
    const double d3  = 0.3333333333333333333; 

    Vector Stress(TrialStress.size());
    SpectralDecomposition(TrialStress, PrincipalStress, EigenVectors);
    IdentifyMaximunAndMinumumPrincipalStres_CalculateOrder(PrincipalStress, Order);
    noalias(Sigma)  = PrincipalStress;
    noalias(Stress) = TrialStress;
    if(CheckPlasticAdmisibility(PrincipalStress))
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
      
      AssembleUpdateStressAndStrainTensor(Sigma,  EigenVectors,  Order, StrainVector, Stress);
      CalculatePlasticStrain(StrainVector, StressVector, mPlastic_strain, Sigma[2]);
      
      /*
      Matrix PlasticTensor;  PlasticTensor.resize(3,3,false); 
      PlasticTensor = ZeroMatrix(3,3);
      /// Updating the  elastic plastic strain tensor
      for(unsigned int i=0; i<3; i++)     
          noalias(PlasticTensor) +=  mPrincipalPlasticStrain_current[i] * Matrix(outer_prod(EigenVectors[Order[i]],  EigenVectors[Order[i]]));
      mplastic_strain[0] = PlasticTensor(0,0);
      mplastic_strain[1] = PlasticTensor(1,1);
      mplastic_strain[2] = 2.00 * PlasticTensor(1,0);
      mplastic_strain[3] = PlasticTensor(2,2);
      */
    }   

    CalculateElasticStrain(Stress, mElastic_strain);
    //KRATOS_WATCH(Stress)
    //KRATOS_WATCH(Sigma)
    //KRATOS_WATCH(mElastic_strain_old)
    //KRATOS_WATCH(mElastic_strain)
    //KRATOS_WATCH("--------------------------")
    /// computing the pressure
    mpressure       =  d3 * (Sigma[0] + Sigma[1] + Sigma[2]);
    StressVector[0] = Stress[0];
    StressVector[1] = Stress[1];
    StressVector[2] = Stress[2];
}


//***********************************************************************
//***********************************************************************

bool Standard_Morh_Coulomb_Yield_Function::ReturnMappingToMainPlane(const array_1d<double,3>& PrincipalStress, const array_1d<unsigned int,3>& order, array_1d<double,3>& Sigma)
{
 
const double d180       =   0.005555555555555555555;  
const double& E         =   (*mprops)[YOUNG_MODULUS];
const double& NU        =   (*mprops)[POISSON_RATIO];
const double G          =   0.5 * E / (1.00 + NU);
const double K          =   E / (3.00 * (1.00-2.00 * NU) );
const double toler      =   1E-6;
const double& friction  =  (*mprops)[INTERNAL_FRICTION_ANGLE]; 
const double& dilatancy =  (*mprops)[DILATANCY_ANGLE];
const unsigned max      =  100;


unsigned iter  = 0;
double dgama   = 0.00;
double ddgama  = 0.00;
double denom   = 0.00;
double sinphi  =   std::sin(PI * friction * d180);
double cosphi  =   std::cos(PI * friction * d180);
double sinpsi  =   std::sin(PI * dilatancy* d180);
//double cospsi  =   std::cos(PI * dilatancy* d180); 


// Start Newton-Raphson iterations for DGAMA
double sigma_ef = (PrincipalStress[0] - PrincipalStress[2]) + (PrincipalStress[0] + PrincipalStress[2]) * sinphi;
double phia     = sigma_ef - 2.00 *  cosphi * mcurrent_cohesion;
double res      = phia/std::fabs(sigma_ef);

double fact              = 0.00; 
double Partial_Cohesion  = 0.00;  

Vector Imput_Parameters;
Imput_Parameters.resize(3);
Imput_Parameters[0] =  0.00; 
Imput_Parameters[1] =  0.00; 
Imput_Parameters[2] =  0.00; 


Partial_Cohesion       =  mpSofteningBehavior_Cohesion->FirstDerivateFunctionBehavior(Imput_Parameters);  
const double d3        =  0.3333333333333333;
//const double raiz2d3   =  0.8164958092773; 
double Ppvs = 0.00;           /// principal plastic volumetric strain  
array_1d<double,3> Ppds;      /// principal plastic desviatoric strain
array_1d<double,3> Pps;       /// principal plastic  strain
array_1d<double,3> I; 
I[0] = 1.00;
I[1] = 1.00;
I[2] = 1.00;

while(fabs(res)>toler && iter++ < max )
   {  
     /// IS OK 
    denom   =  -4 * G * (1.00 + d3 * sinpsi * sinphi) - 4.00 * K * sinpsi * sinphi - 4.00 * cosphi * cosphi *  Partial_Cohesion;  
    ddgama  =  (-phia)/denom;
    dgama  +=   ddgama;

    
    /// volumetric and desviatoric plastic strain
     Pps[0]        = 1.00 + sinpsi;
     Pps[1]        = 0.00;
     Pps[2]        = sinpsi - 1.00;
     Pps          *= dgama;
     Ppvs          = Pps[0] + Pps[1] + Pps[2]; 
     noalias(Ppds) = Pps - d3 * Ppvs * I;
    
     ///associated hardening accumalated plastic strain
    fact  = 2.00 * cosphi * dgama; 
    mmorh_coulomb_maccumulated_plastic_strain_current = mmorh_coulomb_maccumulated_plastic_strain_old +  fact;
    
    ///Updating cohesion
    Imput_Parameters[1]  =  mmorh_coulomb_maccumulated_plastic_strain_current;
    mcurrent_cohesion    =  mpSofteningBehavior_Cohesion->FunctionBehavior(Imput_Parameters);
    Partial_Cohesion     =  mpSofteningBehavior_Cohesion->FirstDerivateFunctionBehavior(Imput_Parameters);   
    
     
     if(mcohesion <= toler){
        mcohesion = 0.00;
        noalias(Sigma) = ZeroVector(3);
        return true; }
    
    sigma_ef = (PrincipalStress[0] - PrincipalStress[2]) + (PrincipalStress[0] + PrincipalStress[2]) * sinphi;
    phia     = sigma_ef - 2.00 *  cosphi * mcurrent_cohesion - 4.00 * K * dgama * sinphi * sinpsi -4.00 * G * dgama * (1.00 + d3 * sinphi * sinpsi);  
    
    /// Check convergence
    res = std::fabs(phia);
    if(sigma_ef!=0) res = res / std::fabs(sigma_ef);    
   }

 
  if(iter>=max || dgama<0.00)
  {
    KRATOS_WATCH(res) 
    std::cout<<"RETURN MAPPING TO MAIN PLANE MORH COULOMB  NOT CONVERGED"<< std::endl;
    KRATOS_ERROR(std::logic_error,  "RETURN MAPPING TO MAIN PLANE MORH COULOMB  NOT CONVERGED" , "");
  }
    
  /// Check validity of 1-vector return (check sextant of converged stress)
  //Sigma[0] = PrincipalStress[0] - ( 2.00 * G * (1.00 + sinpsi/3.00) + 2.00 * K * sinpsi) * dgama; 
  //Sigma[1] = PrincipalStress[1] + (4.00 * G /3.00 - 2.00 * K) * sinpsi  * dgama; 
  //Sigma[2] = PrincipalStress[2] + ( 2.00 * G * (1.00 - sinpsi/3.00) - 2.00 * K * sinpsi) * dgama; 
  noalias(Sigma) = PrincipalStress - 2.00 * G * Ppds - K * Ppvs * I;
  bool check = CheckValidity(Sigma);
  
  /*
  if (check==true)
  {
  /// updating the correct principal pastic strain
  Vector PPS_bar(3);
  ComputePlasticStrainBar(mplastic_strain_old ,m_inv_DeltaF, PPS_bar);  
  noalias(mPrincipalPlasticStrain_current) = PPS_bar + Pps;
  //mPrincipalPlasticStrain_current[0] =  PPS_bar[0] +  dgama * (1.00  + sinpsi);
  //mPrincipalPlasticStrain_current[1] =  PPS_bar[1]  + 0.00;
  //mPrincipalPlasticStrain_current[2] =  PPS_bar[2]  + dgama * (sinpsi - 1.00); 
  }
  */
  return check;
}


//***********************************************************************
//***********************************************************************

bool Standard_Morh_Coulomb_Yield_Function::TwoVectorReturnToEdges(const array_1d<double,3>& PrincipalStress, const array_1d<unsigned int,3>& order, const bool& edges, array_1d<double,3>& Sigma)
{
  const double d3         =   0.3333333333333333;
  const double d180       =   0.005555555555555555555;  
  const double& E         =   (*mprops)[YOUNG_MODULUS];
  const double& NU        =   (*mprops)[POISSON_RATIO];
  const double G          =   0.5 * E / (1.00 + NU);
  const double K          =   E / (3.00 * (1.00-2.00 * NU) );
  const double toler      =   1E-6;
  const double& friction  =  (*mprops)[INTERNAL_FRICTION_ANGLE]; 
  const double& dilatancy =  (*mprops)[DILATANCY_ANGLE];
  const unsigned max      =  100;
  

 
  
  Matrix d                    = ZeroMatrix(2,2);
  Matrix d_inv                = ZeroMatrix(2,2);
  array_1d<double,2> dgama    = ZeroVector(2);
  array_1d<double,2> ddgama   = ZeroVector(2);
  array_1d<double,2> residual = ZeroVector(2); 
 
  
  double sinphi  =   std::sin(PI * friction * d180);
  double cosphi  =   std::cos(PI * friction * d180);
  double sinpsi  =   std::sin(PI * dilatancy* d180);
//  double cospsi  =   std::cos(PI * dilatancy* d180); 
  

  double sigma_a = 0.00;
  double sigma_b = 0.00;
  double Partial_Cohesion  = 0.00;  
  
//  const double raiz2d3 = 0.8164958092773; 
  double aux     = 0.00;
  double a       = 0.00;
  double b       = 0.00;
  unsigned iter  = 0;
  
  
  double fact1   =  2.00 *  cosphi * mcurrent_cohesion;
  double sum     =  0.00; 
  
  sigma_a     = (PrincipalStress[0] - PrincipalStress[2]) + (PrincipalStress[0] + PrincipalStress[2]) * sinphi;  
  a           = 4.00 * K * sinphi * sinpsi + 4.00 * G * (1.00 + d3 * sinphi * sinpsi);
  if(edges==true){  ///rigth edges    
    sigma_b   =  (PrincipalStress[0] - PrincipalStress[1]) + (PrincipalStress[0] + PrincipalStress[1]) * sinphi;
    b         = 2.00 * G * (1.00 + sinphi + sinpsi - d3 * sinpsi * sinphi) + 4.00 * K * sinphi * sinpsi;
  }
  else{   ///left edges
    sigma_b   = (PrincipalStress[1] - PrincipalStress[2]) + (PrincipalStress[1] + PrincipalStress[2]) * sinphi;
    b         = 2.00 * G * (1.00 - sinphi - sinpsi - d3 * sinpsi * sinphi) + 4.00 * K * sinphi * sinpsi;
  }
    
  residual[0]   = sigma_a - fact1; 
  residual[1]   = sigma_b - fact1;  
  KRATOS_WATCH(residual)
  
  int singular  = 0.00;
  double norma  = norm_2(residual);
//  double phipsi = 0.00;
  Vector Imput_Parameters;
  Imput_Parameters.resize(3);
  Imput_Parameters[0] =  0.00;
  Imput_Parameters[1] =  0.00;
  Imput_Parameters[2] =  0.00;
  Partial_Cohesion    =  mpSofteningBehavior_Cohesion->FirstDerivateFunctionBehavior(Imput_Parameters);   
 


  double Ppvs      = 0.00;      /// principal plastic volumetric strain  
  array_1d<double,3> Ppds;      /// principal plastic desviatoric strain
  array_1d<double,3> Pps;       /// principal plastic  strain
  array_1d<double,3> I; 
  I[0] = 1.00;
  I[1] = 1.00;
  I[2] = 1.00;
  
  while(norma > toler && iter++ <max)
{

  
  //Compute residual derivative matrix
  d(0,0) = -a - 4.00 * cosphi * cosphi *  Partial_Cohesion;
  d(0,1) = -b - 4.00 * cosphi * cosphi *  Partial_Cohesion;
  d(1,0) = -b - 4.00 * cosphi * cosphi *  Partial_Cohesion;
  d(1,1) = -a - 4.00 * cosphi * cosphi *  Partial_Cohesion;
     
  singular =  SD_MathUtils<double>::InvertMatrix(d, d_inv);
  ddgama   = -Vector(prod(d_inv, residual));
  
  //Compute Newton-Raphson increment and update variables DGAMA and DGAMB
  noalias(dgama) += ddgama;  
  aux = 2.00 * cosphi *(dgama[0] + dgama[1]);
  mmorh_coulomb_maccumulated_plastic_strain_current = mmorh_coulomb_maccumulated_plastic_strain_old + aux;
  
  
  if(edges==true)  //rigth edges
   {
    /// volumetric and desviatoric plastic strain
    Pps[0]        = (dgama[0] + dgama[1])*(1.00 + sinpsi);
    Pps[1]        = dgama[1]*(sinpsi - 1.00);
    Pps[2]        = dgama[0]*(sinpsi - 1.00);
    Ppvs          = Pps[0] + Pps[1] + Pps[2]; 
    noalias(Ppds) = Pps - d3 * Ppvs * I;
   }
    else   //left edges
   {
    /// volumetric and desviatoric plastic strain
    Pps[0]        = dgama[0]*(sinpsi + 1.00); 
    Pps[1]        = dgama[1]*(sinpsi + 1.00);
    Pps[2]        = (dgama[0] + dgama[1])*(sinpsi-1.00);
    Ppvs          = Pps[0] + Pps[1] + Pps[2]; 
    noalias(Ppds) = Pps - d3 * Ppvs * I;
   }
    
    /// Compute Internal Variables 
    Imput_Parameters[1] =  mmorh_coulomb_maccumulated_plastic_strain_current;
    mcurrent_cohesion   =  mpSofteningBehavior_Cohesion->FunctionBehavior(Imput_Parameters);     
 
  if(mcohesion < 0.0) {
    mcohesion = 0.00;
    noalias(Sigma) = ZeroVector(3); 
    return true; }
  

fact1       =  2.00 *  cosphi * mcurrent_cohesion;
residual[0] =  sigma_a - a * dgama[0] - b * dgama[1] - fact1;  
residual[1] =  sigma_b - b * dgama[0] - a * dgama[1] - fact1;


//Check convergence
norma = std::fabs(residual[0]) + std::fabs(residual[1]);  
sum   = std::fabs(sigma_a) + std::fabs(sigma_b);
if(sum!=0) {norma = norma / sum;}
}

if( (iter>=max) || (dgama[0]< 0.00 && dgama[1]<0.00) )
{
   KRATOS_WATCH(dgama) 
   KRATOS_WATCH(iter)
   std::cout<< "RETURN MAPPING TO MAIN PLANE AND RIGTH O LEFT MORH COULOMB  NOT CONVERGED" << std::endl;
   KRATOS_ERROR(std::logic_error,  "RETURN MAPPING TO MAIN PLANE AND RIGTH O LEFT MORH COULOMB  NOT CONVERGED" , "");
}


//const double aux1 = ( 2.00 * G * (1.00 + sinpsi/3.00) + 2.00 * K * sinpsi);
//const double aux2 = ( 4.00 * G /3.00 - 2.00 * K) * sinpsi; 
//const double aux3 = ( 2.00 * G * (1.00 - sinpsi/3.00) - 2.00 * K * sinpsi);
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
  /*
  if (check==true)
  { 
    Vector PPS_bar;
    ComputePlasticStrainBar(mplastic_strain_old ,m_inv_DeltaF, PPS_bar);
    noalias(mPrincipalPlasticStrain_current) = PPS_bar + Pps;
    
   if(edges ==true)  //rigth
   {
    //updating the correct principal pastic strain
    mPrincipalPlasticStrain_current[0] = mPrincipalPlasticStrain_old[0]   PPS_bar[0] +  (dgama[0] + dgama[1])  * (1.00  + sinpsi);
    mPrincipalPlasticStrain_current[1] = mPrincipalPlasticStrain_old[1]   PPS_bar[1] + dgama[1] * (sinpsi - 1.00);
    mPrincipalPlasticStrain_current[2] = mPrincipalPlasticStrain_old[2]   PPS_bar[2]  + dgama[0] * (sinpsi - 1.00); 
   }
   else
   {
    //updating the correct principal pastic strain
    mPrincipalPlasticStrain_current[0] =  mPrincipalPlasticStrain_old[0]  PPS_bar[0]  + dgama[0] * (sinpsi + 1.00); 
    mPrincipalPlasticStrain_current[1] =  mPrincipalPlasticStrain_old[1]  PPS_bar[1]  + dgama[1] * (sinpsi + 1.00);
    mPrincipalPlasticStrain_current[2] =  mPrincipalPlasticStrain_old[2]  PPS_bar[2]  + (dgama[0] + dgama[1])  * (sinpsi - 1.00); 
   }
  }
  */
  return check;
}


///WARNING = For some reason the definition of the desciatoric flow rule in the tesis of owen is erroneuos
/// I take the  The Mohr–Coulomb model. Return mapping to apex given in box 8.7 of book the owen
void Standard_Morh_Coulomb_Yield_Function::ReturnMappingToApex(const array_1d<double,3>& PrincipalStress, array_1d<double, 3 >& Sigma)
{

//const double d3         =   0.3333333333333333;
const double d180       =   0.005555555555555555555;  
const double& E         =   (*mprops)[YOUNG_MODULUS];
const double& NU        =   (*mprops)[POISSON_RATIO];
//const double G          =   0.5 * E / (1.00 + NU);
const double K          =   E / (3.00 * (1.00-2.00 * NU) );
const double toler      =   1E-6;
const double& friction  =  (*mprops)[INTERNAL_FRICTION_ANGLE]; 
const double& dilatancy =  (*mprops)[DILATANCY_ANGLE];
const unsigned max      =  100; 


double Partial_Cohesion =   0.00;  
const double sinphi  =   std::sin(PI * friction * d180);
const double cosphi  =   std::cos(PI * friction * d180);
const double sinpsi  =   std::sin(PI * dilatancy * d180);
const double cotphi  =   cosphi/sinphi; 
const double alpha   =   cosphi/sinpsi;


unsigned int iter  = 0.00;
double d           = 0.00;
double r           = 0.00;
double p_trial     = (1.00/3.00) * ( PrincipalStress[0] + PrincipalStress[1] + PrincipalStress[2] );
double p           = 0.00;
double fact        = 0.00;

//const double raiz2d3   = 0.8164958092773;  
double dgama_b         = 0.00;
double ddgama_b        = 0.00;  // volumetric part


 
r  =  mcohesion * cotphi -  p_trial;
Vector Imput_Parameters;
Imput_Parameters.resize(3);
Imput_Parameters[0] =  0.00;
Imput_Parameters[1] =  0.00;
Imput_Parameters[2] =  0.00;
Partial_Cohesion    =  mpSofteningBehavior_Cohesion->FirstDerivateFunctionBehavior(Imput_Parameters);   



while(std::fabs(r) > toler && iter++<max ) 
{  
   d         = K +  Partial_Cohesion * cotphi * alpha; 
   ddgama_b  = -r / d;
   dgama_b  += ddgama_b;
     
   
   ///associated hardening 
   fact = alpha *  dgama_b; 
   mmorh_coulomb_maccumulated_plastic_strain_current = mmorh_coulomb_maccumulated_plastic_strain_old +  fact;
   
   
   /// Compute Internal Variables   
    Imput_Parameters[1] =  0.00;
    mcurrent_cohesion   =  mpSofteningBehavior_Cohesion->FunctionBehavior(Imput_Parameters);
    Partial_Cohesion    =  mpSofteningBehavior_Cohesion->FirstDerivateFunctionBehavior(Imput_Parameters);   
  
    if(mcurrent_cohesion < 0.0) 
    {
      mcurrent_cohesion = 0.00;
      noalias(Sigma) = ZeroVector(3);   
      break; 
    }
 
   p =  p_trial - K *  dgama_b;   
   r =  mcurrent_cohesion * cotphi - p;
}

if(iter>=max)
{ 
   KRATOS_WATCH(r)
   KRATOS_WATCH(PrincipalStress) 
   std::cout<< "RETURN MAPPING TO APEX  NOT CONVERGED"<< std::endl;
   KRATOS_ERROR(std::logic_error,  "RETURN MAPPING TO APEX  NOT CONVERGED" , "");
}

   
Sigma[0] = p;
Sigma[1] = p;
Sigma[2] = p;

/*
Vector PPS_bar;
ComputePlasticStrainBar(mplastic_strain_old ,m_inv_DeltaF, PPS_bar);
//updating the correct principal pastic strain
mPrincipalPlasticStrain_current[0] =   PPS_bar[0] + dgama_b/3.00;  //dgama_a * des[0]  + dgama_b/3.00;  
mPrincipalPlasticStrain_current[1] =   PPS_bar[1] + dgama_b/3.00;  //dgama_a * des[1]  + dgama_b/3.00;
mPrincipalPlasticStrain_current[2] =   PPS_bar[2] + dgama_b/3.00;  //dgama_a * des[2]  + dgama_b/3.00;
*/
}

    
//***********************************************************************
//***********************************************************************

bool Standard_Morh_Coulomb_Yield_Function::CheckValidity(const array_1d<double,3>& Sigma)
{
  bool check   = false;
  array_1d<double,3> Aux_Sigma;
  Aux_Sigma[0] = fabs(Sigma[0]);
  Aux_Sigma[1] = fabs(Sigma[1]);
  Aux_Sigma[2] = fabs(Sigma[2]); 
  double delta = (*max_element(Aux_Sigma.begin(), Aux_Sigma.end())) * 1.00E-6; 
  if( (Sigma[0] + delta) >= Sigma[1] && (Sigma[1] + delta ) >= Sigma[2]){ check = true;}
  return check;
}
 
 
//***********************************************************************
//***********************************************************************

bool Standard_Morh_Coulomb_Yield_Function::ReturnToEdges(const array_1d<double,3>& PrincipalStress)
{
    const double& dilatance =   (*mprops)[DILATANCY_ANGLE];   
    const double sinpsi     =   std::sin(PI * dilatance / 180.00);    
    bool  return_rigth      =   false;  // left edges
    double scaprd           = (1.00 - sinpsi) * PrincipalStress[0] - 2.00 * PrincipalStress[1] + (1.00 + sinpsi) * PrincipalStress[2]; 
    if(scaprd > 0.00) {return_rigth = true;} // return to rigth edges
    return return_rigth;  
}

    
//***********************************************************************
//***********************************************************************


void Standard_Morh_Coulomb_Yield_Function::FinalizeSolutionStep()
{
   mmorh_coulomb_maccumulated_plastic_strain_old =  mmorh_coulomb_maccumulated_plastic_strain_current; 
   mcohesion                                     =  mcurrent_cohesion;
   noalias(mPrincipalPlasticStrain_old)          =  mPrincipalPlasticStrain_current; 
   noalias(mElastic_strain_old)                  =  mElastic_strain;
   noalias(mPlastic_strain_old)                  =  mPlastic_strain;
}

    
//***********************************************************************
//***********************************************************************


void Standard_Morh_Coulomb_Yield_Function::UpdateMaterial()

{
   mmorh_coulomb_maccumulated_plastic_strain_current = mmorh_coulomb_maccumulated_plastic_strain_old;
   mcurrent_cohesion                                 = mcohesion; 
   noalias(mPrincipalPlasticStrain_current)          = mPrincipalPlasticStrain_old; 
   noalias(mElastic_strain)                          = mElastic_strain_old;
   noalias(mPlastic_strain)                          = mPlastic_strain_old;
}

    
//***********************************************************************
//***********************************************************************
void Standard_Morh_Coulomb_Yield_Function::GetValue(const Variable<double>& rVariable, double& Result)
      {
	
	if(rVariable==COHESION)
	  Result = mcurrent_cohesion;
	
	if(rVariable == DAMAGE)
	   Result = mmorh_coulomb_maccumulated_plastic_strain_current; 
	
	if(rVariable == PRESSURE)
	   Result = mpressure;    
      }


    
//***********************************************************************
//***********************************************************************

void Standard_Morh_Coulomb_Yield_Function::GetValue(const Variable<Matrix>& rVariable, Matrix& Result)
 {   unsigned int size = SizePlasticStrain(); 
     if(rVariable==GREEN_LAGRANGE_PLASTIC_STRAIN_TENSOR){
           Result = ZeroMatrix(1,size);
           for(unsigned int i = 0; i< size; i++){
                 Result(0,i) = mPlastic_strain(i);
	   }
     }
 }

//***********************************************************************
//***********************************************************************

void Standard_Morh_Coulomb_Yield_Function::GetValue(double& Result)
 { 
 }
    
//***********************************************************************
//***********************************************************************

void Standard_Morh_Coulomb_Yield_Function::GetValue(Matrix& Result)
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

         
//***********************************************************************
//***********************************************************************

void Standard_Morh_Coulomb_Yield_Function::GetValue(const Variable<Vector>& rVariable, Vector& Result)
 {
   unsigned int size = SizePlasticStrain(); 
   if(rVariable==ALMANSI_ELASTIC_STRAIN)
   {
     Result.resize(size, false);
     noalias(Result) = mElastic_strain;
   }
   if(rVariable==ALMANSI_PLASTIC_STRAIN)
   {
     Result.resize(size, false);
     noalias(Result) = mPlastic_strain;
   }
   
   
 }

}
