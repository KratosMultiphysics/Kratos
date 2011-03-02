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

#include "fluency_criteria/morh_coulomb_yield_function.h"




namespace Kratos
  {

           
            Morh_Coulomb_Yield_Function::Morh_Coulomb_Yield_Function(){}
    
            Morh_Coulomb_Yield_Function::Morh_Coulomb_Yield_Function( 
                                                                      SoftHardPointerType SofteningBehavior,
                                                                      myState State, 
								      myPotencialPlastic PotencialPlastic)
	    :FluencyCriteria()
	    {
              mState              = State;
	      mpSofteningBehavior = SofteningBehavior;
              mPotencialPlastic   = PotencialPlastic;
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
  mplastic_strain                                   =  ZeroVector(4);
  mplastic_strain_old                               =  ZeroVector(4); 
  mPrincipalPlasticStrain_current                   =  ZeroVector(3);
  mPrincipalPlasticStrain_old                       =  ZeroVector(3);
  mmorh_coulomb_maccumulated_plastic_strain_old     = 0.00; 
  mmorh_coulomb_maccumulated_plastic_strain_current = 0.00; 

}
 
 
bool Morh_Coulomb_Yield_Function::CheckPlasticAdmisibility(const Vector& Stress)
{
  const double tol        = 1.0E-10; 
  const double& friction  = (*mprops)[INTERNAL_FRICTION_ANGLE];  
  const double sinphi     = std::sin(PI * friction  / 180.00);
  const double cosphi     = std::cos(PI * friction  / 180.00);
  
  
  // Check plastic admissibility
  double sigma_ef = (Stress[0] - Stress[2]) + (Stress[0] + Stress[2]) * sinphi;
  double phia     = sigma_ef - 2.00 *  cosphi * mcohesion;
  double res      = phia;
  if(mcohesion!=0.00) { res = res/std::fabs(mcohesion); } 
  if(res > tol) return true;
  else
    return false;
    
}
  
  
void Morh_Coulomb_Yield_Function::ReturnMapping(const Vector& StrainVector, Vector& StressVector)
{
    
    array_1d<double,3> PrincipalStress;
    array_1d<array_1d < double,3 > ,3> EigenVectors;
    array_1d<unsigned int,3> Order;
        
    SpectralDecomposition(StressVector, PrincipalStress, EigenVectors);
    IdentifyMaximunAndMinumumPrincipalStres_CalculateOrder(PrincipalStress, Order);
        
    if(CheckPlasticAdmisibility(PrincipalStress))
    {
      // plastic Step: Apply return mapping
      // identify possible edge return:either right or left of the main plain             
      const bool edge = ReturnToEdges(PrincipalStress);
      
      // Apply one-vector return mapping first (return to MAIN PLANE)
      array_1d<double,3> Sigma;
      if(ReturnMappingToMainPlane(PrincipalStress, Order,  Sigma)==false){
	 //Apply two-vector return mapping to appropriate EDGES
	 if(TwoVectorReturnToEdges (PrincipalStress, Order, edge, Sigma)==false){
	    //Apply multi-vector return mapping to APEX
	    ReturnMappingToApex(PrincipalStress, Sigma);
	 }
	    
      }
      AssembleUpdateStressAndStrainTensor(Sigma,  EigenVectors,  Order, StrainVector, StressVector);

    }
    
     // Updating the plastic strain tensor
     CalculatePlasticStrain(StrainVector, StressVector);
   
}


bool Morh_Coulomb_Yield_Function::ReturnMappingToMainPlane(const array_1d<double,3>& PrincipalStress, const array_1d<unsigned int,3>& order, array_1d<double,3>& Sigma)
{
  
  
const double& E         =   (*mprops)[YOUNG_MODULUS];
const double& NU        =   (*mprops)[POISSON_RATIO];

const double G          = 0.5 * E / (1.00 + NU);
const double K          = E / (3.00 * (1.00-2.00 * NU) );
const double toler      = 1E-6;
const unsigned max      = 1000;

double   dgama    = 0.00;
double   ddgama   = 0.00;
unsigned iter     = 0;
double   denom    = 0.00;

double sinphi             =   std::sin(PI * this->minternal_friction_angle  / 180.00);
double cosphi             =   std::cos(PI * this->minternal_friction_angle  / 180.00);
double sinpsi             =   std::sin(PI * this->mcurrent_dilatancy_angle  / 180.00);
double cospsi             =   std::cos(PI * this->mcurrent_dilatancy_angle  / 180.00); 


// Start Newton-Raphson iterations for DGAMA
double sigma_ef = (PrincipalStress[0] - PrincipalStress[2]) + (PrincipalStress[0] + PrincipalStress[2]) * sinphi;
double phia     = sigma_ef - 2.00 *  cosphi * this->mcurrent_cohesion;
double res      = phia/std::fabs(sigma_ef);

double fact              = 0.00; 
double Partial_Cohesion  = 0.00;  
double Partial_Friction  = 0.00;
double Partial_Dilatancy = 0.00;
double Partial_Ep_gama_a = 0.00; 

while(fabs(res)>toler && iter++ < max )
   {  
    
    Partial_Cohesion  = 0.00;   // double H = 0.00;     
    Partial_Friction  = 0.00; 
    Partial_Dilatancy = 0.00;
    
    Partial_Ep_gama_a = std::sqrt((4.00/3.00)  * (1.00 + sinpsi*sinpsi) ); 
     
    denom   =  (PrincipalStress[0] + PrincipalStress[2]) * cosphi * Partial_Friction * Partial_Ep_gama_a
               - 2.00 * cosphi * Partial_Cohesion * Partial_Ep_gama_a 
               + 2.00 * this->mcurrent_cohesion * sinphi * Partial_Ep_gama_a
               - 4.00 * G - 4.00 * (K + G/3.00)*(dgama * cosphi * sinpsi * Partial_Friction * Partial_Ep_gama_a
               + dgama * cospsi * sinphi * Partial_Dilatancy * Partial_Ep_gama_a)
               - 4.00 * (K + (1.00/3.00)* G ) * sinphi * sinpsi; 
    
       
    ddgama  =  -phia/denom;
    dgama  +=   ddgama;

    
    // von mises acumulated plastic strain
    fact  = dgama*(std::sqrt((4.00/3.00) * (1.00 + sinpsi*sinpsi))); 
    mmorh_coulomb_maccumulated_plastic_strain_current = mmorh_coulomb_maccumulated_plastic_strain_old +  fact;
    
    // Compute Internal Variables   
    /*
    COHESION
    DILATANCY_ANGLE
    FRICTION
    */
    
     // updating cos and sin
     sinphi             =   std::sin(PI * this->minternal_friction_angle  / 180.00);
     cosphi             =   std::cos(PI * this->minternal_friction_angle  / 180.00);
     sinpsi             =   std::sin(PI * this->mcurrent_dilatancy_angle  / 180.00);
     cospsi             =   std::cos(PI * this->mcurrent_dilatancy_angle  / 180.00); 
    
    
    if(mcohesion < 0.0) {
    mcohesion = 0.00;
    Sigma[0] = 1E-14 ;
    Sigma[1] = 1E-14 ; 
    Sigma[2] = 1E-14 ; 
    return true; }
    
    sigma_ef = (PrincipalStress[0] - PrincipalStress[2]) + (PrincipalStress[0] + PrincipalStress[2]) * sinphi;
    phia     = sigma_ef - 2.00 *  cosphi * this->mcurrent_cohesion - 4.00 * K * dgama * sinphi * sinpsi
               -4.00 * G * dgama * (1.00 + (1.00/3.00) * sinphi * sinpsi ) ;  
	        
    // Check convergence
    res = std::fabs(phia);
    if(sigma_ef!=0) res = res / std::fabs(sigma_ef);    
   }


  if(iter>=max)
  {
    KRATOS_ERROR(std::logic_error,  "RETURN MAPPING TO MAIN PLANE MORH COULOMB  NOT CONVERGED" , "");
  }
    
  // Check validity of 1-vector return (check sextant of converged stress)
  Sigma[0] = PrincipalStress[0] - ( 2.00 * G * (1.00 + sinpsi/3.00) + 2.00 * K * sinpsi) * dgama; 
  Sigma[1] = PrincipalStress[1] + (4.00 * G /3.00 - 2.00 * K) * sinpsi  * dgama; 
  Sigma[2] = PrincipalStress[2] + ( 2.00 * G * (1.00 - sinpsi/3.00) - 2.00 * K * sinpsi) * dgama; 
  
  bool check = CheckValidity(Sigma);
  if (check==true)
  {
  //updating the correct principal pastic strain
  mPrincipalPlasticStrain_current[0] = mPrincipalPlasticStrain_old[0] +  dgama * (1.00  + sinpsi);
  mPrincipalPlasticStrain_current[1] = mPrincipalPlasticStrain_old[1] + 0.00;
  mPrincipalPlasticStrain_current[2] = mPrincipalPlasticStrain_old[2] + dgama * (sinpsi - 1.00); 
  }
  
  return check;
  
}




bool Morh_Coulomb_Yield_Function::TwoVectorReturnToEdges(const array_1d<double,3>& PrincipalStress, const array_1d<unsigned int,3>& order, const bool& edges, array_1d<double,3>& Sigma)
{
    
  const double& E         =   (*mprops)[YOUNG_MODULUS];
  const double& NU        =   (*mprops)[POISSON_RATIO];
  const double G          = 0.5 * E / (1.00 + NU);
  const double K          = E / (3.00 * (1.00-2.00 * NU) );
  const double toler      = 1E-6;
  const unsigned max      = 1000;
  unsigned int iter       = 0;

 
  
  Matrix d                    = ZeroMatrix(2,2);
  Matrix d_inv                = ZeroMatrix(2,2);
  array_1d<double,2> dgama    = ZeroVector(2);
  array_1d<double,2> ddgama   = ZeroVector(2);
  array_1d<double,2> residual = ZeroVector(2); 
 
  
  double sinphi  =   std::sin(PI * this->minternal_friction_angle  / 180.00);
  double cosphi  =   std::cos(PI * this->minternal_friction_angle  / 180.00);
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
  double aux_1   = 0.00;
  double aux_2   = 0.00;
  double aux_3   = 0.00;
  double aux_4   = 0.00;
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
  
  while(norma > toler && iter++ <max)
{
   
    
    Partial_Cohesion  = 0.00;  
    Partial_Friction  = 0.00;
    Partial_Dilatancy = 0.00;
    
    Partial_Ep_gama_a = (2.00/3.00) * (C[2] * (dgama[0] + dgama[1] ) + C[1] * dgama[0])/aux;
    Partial_Ep_gama_b = (2.00/3.00) * (C[2] * (dgama[0] + dgama[1] ) + C[1] * dgama[1])/aux;
  
    if(edges==true)  
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
  
   // von mises acumulated plastic strain
  if(edges==true)  //rigth edges
   {
    aux_1   = (dgama[0] + dgama[1]) *  (dgama[0] + dgama[1]);
    aux_2   = (dgama[0] * dgama[0]) +  (dgama[1] * dgama[1]);
    aux_3   = (sinpsi+1.00) * (sinpsi+1.00);
    aux_4   = (sinpsi-1.00) * (sinpsi-1.00);
    aux     =  raiz2d3 * std::sqrt( aux_1 * aux_3 + aux_2 * aux_4 );
    mmorh_coulomb_maccumulated_plastic_strain_current = mmorh_coulomb_maccumulated_plastic_strain_old + aux;
   }
    else   //left edges
   {
    aux_1   = (dgama[0] + dgama[1]) *  (dgama[0] + dgama[1]);
    aux_2   = (dgama[0] * dgama[0]) +  (dgama[1] * dgama[1]);
    aux_3   = (sinpsi+1.00) * (sinpsi+1.00);
    aux_4   = (sinpsi-1.00) * (sinpsi-1.00);
    aux     =  raiz2d3 * std::sqrt( aux_1 * aux_4 + aux_2 * aux_3 );
    mmorh_coulomb_maccumulated_plastic_strain_current = mmorh_coulomb_maccumulated_plastic_strain_old + aux;
   }
      
    // Compute Internal Variables   
    /*
    COHESION
    DILATANCY_ANGLE
    FRICTION
    */ 
     // updating cos and sin
     sinphi             =   std::sin(PI * this->minternal_friction_angle  / 180.00);
     cosphi             =   std::cos(PI * this->minternal_friction_angle  / 180.00);
     sinpsi             =   std::sin(PI * this->mcurrent_dilatancy_angle  / 180.00);
     cospsi             =   std::cos(PI * this->mcurrent_dilatancy_angle  / 180.00); 


 
  if(mcohesion < 0.0) {
    mcohesion = 0.00;
    Sigma[0] = 1E-14 ;
    Sigma[1] = 1E-14 ; 
    Sigma[2] = 1E-14 ; 
    return true; }
  

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

if(iter>=max)
{
   KRATOS_ERROR(std::logic_error,  "RETURN MAPPING TO MAIN PLANE AND RIGTH O LEFT MORH COULOMB  NOT CONVERGED" , "");
}


const double aux1 = ( 2.00 * G * (1.00 + sinpsi/3.00) + 2.00 * K * sinpsi);
const double aux2 = ( 4.00 * G /3.00 - 2.00 * K) * sinpsi; 
const double aux3 = ( 2.00 * G * (1.00 - sinpsi/3.00) - 2.00 * K * sinpsi);

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

  
  bool check = CheckValidity(Sigma);
  if (check==true)
  { 
   if(edges ==true)  //rigth
   {
    //updating the correct principal pastic strain
    mPrincipalPlasticStrain_current[0] = mPrincipalPlasticStrain_old[0] +  (dgama[0] + dgama[1])  * (1.00  + sinpsi);
    mPrincipalPlasticStrain_current[1] = mPrincipalPlasticStrain_old[1] + dgama[1] * (sinpsi - 1.00);
    mPrincipalPlasticStrain_current[2] = mPrincipalPlasticStrain_old[2] + dgama[0] * (sinpsi - 1.00); 
   }
   else
   {
    //updating the correct principal pastic strain
    mPrincipalPlasticStrain_current[0] =  mPrincipalPlasticStrain_old[0] + dgama[0] * (sinpsi + 1.00); 
    mPrincipalPlasticStrain_current[1] =  mPrincipalPlasticStrain_old[1] + dgama[1] * (sinpsi + 1.00);
    mPrincipalPlasticStrain_current[2] =  mPrincipalPlasticStrain_old[2] + (dgama[0] + dgama[1])  * (sinpsi - 1.00);
     
   } 
  }
  return check;

}



void Morh_Coulomb_Yield_Function::ReturnMappingToApex(const array_1d<double,3>& PrincipalStress, array_1d<double, 3 >& Sigma)
{

const double& E         =   (*mprops)[YOUNG_MODULUS];
const double& NU        =   (*mprops)[POISSON_RATIO];
const double K          =   E / (3.00 * (1.00-2.00 * NU) );
const double G          =   0.5 * E / (1.00 + NU);
const double toler      =   1E-6;
const unsigned max      =   1000;
unsigned int iter       =   0;  


double Partial_Cohesion  = 0.00;  
double Partial_Friction  = 0.00;
double Partial_Ep_gama_b = 0.00; 

double sinphi       =   std::sin(PI * (minternal_friction_angle)  / 180.00);
double cosphi       =   std::cos(PI * (minternal_friction_angle)  / 180.00);
double cotphi       =   cosphi/sinphi; 

double d           = 0.00;
double r           = 0.00;
double p_trial     = (1.00/3.00) * ( PrincipalStress[0] + PrincipalStress[1] + PrincipalStress[2] );
double p           = 0.00;


const double dgama_a   = 1.00/(2.00 * G ); 
const double raiz2d3   = 0.8164958092773;  
double dgama_b         = 0.00;
double ddgama_b        = 0.00;  // volumetric part



array_1d<double, 3> des;
des[0]        =  PrincipalStress[0] - p_trial;
des[1]        =  PrincipalStress[1] - p_trial;
des[2]        =  PrincipalStress[2] - p_trial;
double prod   =  inner_prod(des,des);
double fact_a =  prod * prod *  dgama_a  *  dgama_a ;
double fact   =  1.00;
 
r  =  p_trial - K *  dgama_b - mcohesion * cotphi; 
while(std::fabs(r) > toler && iter++<max ) 
{
  
   Partial_Cohesion  = 0.00;    
   Partial_Friction  = 0.00;
   Partial_Ep_gama_b = (2.00/9.00) *  dgama_b / fact;
  
   d         = -K - ( Partial_Cohesion * cotphi - ( (mcurrent_cohesion)/(sinphi*sinphi)) * (Partial_Friction) ) *Partial_Ep_gama_b;
   ddgama_b  = -r / d;
   dgama_b  += ddgama_b;
     
   // Aculutaded von Misses nod efioned for apex  
   fact_a =  prod * prod *  dgama_a  *  dgama_a ;
   fact   = raiz2d3 * std::sqrt(fact_a + dgama_b*dgama_b/3.00);
   mmorh_coulomb_maccumulated_plastic_strain_current = mmorh_coulomb_maccumulated_plastic_strain_old +  fact;
   
    // Compute Internal Variables   
    /*
    COHESION
    DILATANCY_ANGLE
    FRICTION
    */ 
    // updating cos and sin
    sinphi  =   std::sin(PI * (minternal_friction_angle)  / 180.00);
    cosphi  =   std::cos(PI * (minternal_friction_angle)  / 180.00);
    cotphi  =   cosphi/sinphi; 

    
    if(mcurrent_cohesion < 0.0) {mcurrent_cohesion = 0.00;
      Sigma[0] = 1E-9;  
      Sigma[1] = 1E-9;  
      Sigma[2] = 1E-9; 
      break; }
 
   p = p_trial - K *  dgama_b;   
   r = p - (mcurrent_cohesion) * cotphi;      

}

if(iter>=max)
{ 
   std::cout<< "Norma = " << r << std::endl; 
   KRATOS_WATCH(PrincipalStress)
   KRATOS_WATCH("RETURN MAPPING TO APEX  NOT CONVERGED");
   KRATOS_ERROR(std::logic_error,  "RETURN MAPPING TO APEX  NOT CONVERGED" , "");
}

   
Sigma[0] = p;
Sigma[1] = p;
Sigma[2] = p;

//updating the correct principal pastic strain
mPrincipalPlasticStrain_current[0] =  mPrincipalPlasticStrain_old[0] + dgama_a * des[0]  + dgama_b/3.00;  
mPrincipalPlasticStrain_current[1] =  mPrincipalPlasticStrain_old[1] + dgama_a * des[1]  + dgama_b/3.00;
mPrincipalPlasticStrain_current[2] =  mPrincipalPlasticStrain_old[2] + dgama_a * des[2]  + dgama_b/3.00;

}


bool Morh_Coulomb_Yield_Function::CheckValidity(const array_1d<double,3>& Sigma)
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

bool Morh_Coulomb_Yield_Function::ReturnToEdges(const array_1d<double,3>& PrincipalStress)
{
 
    const double& dilatance =  this->mcurrent_dilatancy_angle; //(*mprops)[DILATANCY_ANGLE];   
    const double sinpsi     =   std::sin(PI * dilatance / 180.00);    
    bool  return_rigth      =   false;           // left edges
    array_1d<double, 3 > T  =   ZeroVector(3); 
    double fact_1           =   3.00 + sinpsi;
    double fact_2           =   3.00 - sinpsi;
    T[0]                    =   2.00 / fact_1  - 1.00/fact_2;  
    T[1]                    =   -1.00 / fact_1 - 1.00/fact_2;
    T[2]                    =   -1.00 / fact_1 - 2.00/fact_2;
    
    //double scaprd         = (1.00 - sinpsi) * PrincipalStress[0] - 2.00 * PrincipalStress[1] + (1.00 + sinpsi) * PrincipalStress[2]; 
    double scaprd           = inner_prod(T, PrincipalStress);
    if(scaprd > 0.00) {return_rigth = true;} // return to rigth edges

    return return_rigth;  
}


void Morh_Coulomb_Yield_Function::FinalizeSolutionStep()
{
   mmorh_coulomb_maccumulated_plastic_strain_old = mmorh_coulomb_maccumulated_plastic_strain_current; 
   mcohesion                                     =  mcurrent_cohesion;
   mdilatancy_angle                              =  mcurrent_dilatancy_angle;
   minternal_friction_angle                      =  mcurrent_minternal_friction_angle;
   noalias(mPrincipalPlasticStrain_old)          =  mPrincipalPlasticStrain_current; 
   noalias(mplastic_strain_old)                  =  mplastic_strain; 
}



void Morh_Coulomb_Yield_Function::UpdateMaterial()

{
   mmorh_coulomb_maccumulated_plastic_strain_current = mmorh_coulomb_maccumulated_plastic_strain_old;
   mcurrent_cohesion                                 = mcohesion; 
   mcurrent_dilatancy_angle                          = mdilatancy_angle; 
   mcurrent_minternal_friction_angle                 = minternal_friction_angle;
   noalias(mPrincipalPlasticStrain_current)          = mPrincipalPlasticStrain_old; 
   noalias(mplastic_strain)                          = mplastic_strain_old; 
}


void Morh_Coulomb_Yield_Function::GetValue(const Variable<double>& rVariable, double& Result)
      {
	
	if(rVariable==COHESION)
	  Result = mcurrent_cohesion;
        
	if(rVariable == DILATANCY_ANGLE)
          Result = mcurrent_dilatancy_angle; 
	
	if(rVariable == INTERNAL_FRICTION_ANGLE)
	  Result = mcurrent_minternal_friction_angle;
	          
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

}




