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
    namespace Plastic_Damage_2D_Auxiliaries
    {
        boost::numeric::ublas::bounded_matrix<double,2,2> mstemp;
	#ifdef _OPENMP
	#pragma omp threadprivate(mstemp)
	#endif
        boost::numeric::ublas::bounded_matrix<double,2,2> msaux;
	#ifdef _OPENMP
	#pragma omp threadprivate(msaux)
	#endif
    } 
    
    
    using namespace Plastic_Damage_2D_Auxiliaries;
  

	/**
	 *	TO BE TESTED!!!
	 */

         PlasticDamage2D::PlasticDamage2D () 
	: ConstitutiveLaw< Node<3> >()
	
	{
	  KRATOS_ERROR(std::logic_error,"Calling the empty constructor.","");
	}

	 PlasticDamage2D::PlasticDamage2D(
         FluencyCriteriaPointer FluencyCriteria,
         FluencyCriteriaPointer FluencyCriteriaTraction,
         SofteningHardeningCriteriaPointer SofteningBehavior, 
         PropertiesPointer Property) 
	: ConstitutiveLaw< Node<3> >()
	{
	      mpFluencyCriteria              = FluencyCriteria;
              mpFluencyCriteria_Traction     = FluencyCriteriaTraction;
              mpSofteningBehavior            = SofteningBehavior;
              mpProperties                   = Property;

             // mFluencyCriteria = GetProperties()[CONSTITUTIVE_LAW]->Clone();
	      //KRATOS_WATCH(*mProperties)
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
	
	double PlasticDamage2D::GetValue( const Variable<double>& rThisVariable )
	{
	  if( rThisVariable == DAMAGE)
	      {return mFt[0];}   //mlocal_fail_factor; } //mplastic_damage ;}
	  else if(rThisVariable == COHESION)
              {return mFt[1];} // mcohesion;}
          else if(rThisVariable == DILATANCY_ANGLE)
              {return mdilatancy_angle*180.00/PI;}
          else if(rThisVariable == FRICTION_INTERNAL_ANGLE)
              {return mFt[2];} // mfriction_angle*180.00/PI;}
          else
               {return 0; }
	}   
	

	Vector PlasticDamage2D::GetValue( const Variable<Vector>& rThisVariable )
	{
	    KRATOS_ERROR(std::logic_error, "Vector Variable case not considered", "");     
        }
	
	Matrix PlasticDamage2D::GetValue( const Variable<Matrix>& rThisVariable )
	{
             KRATOS_ERROR(std::logic_error, "Matrix Variable case not considered", "");
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
    //double local_damage =  fabs(1.00-mplastic_damage ); 
    //Output = sqrt(mE/mDE);
    //if(Output==0){Output = 1;}
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
        mmaxfriction_angle       = (*mpProperties)[MAX_FRICTION_INTERNAL_ANGLE]*PI/180.00;
        mmaxdilatancy_angle      = (*mpProperties)[MAX_DILATANCY_ANGLE]*PI/180.00;
        mcohesion                = (*mpProperties)[FC]/(2.00*tan(mmaxfriction_angle/2.00 + PI/4.00));
        mlength                  = 1.1283791670955 * sqrt(fabs(geom.Area())); ///  ( 4 * A / Pi ) ^( 1/ 2 )  

        mFt[0]                   = (*mpProperties)[FT]; 
        mFt[1]                   = (*mpProperties)[FT];         
        mFt[2]                   = (*mpProperties)[FT]; 

        mplastic_damage          = 0.00;
        mdisipation              = 0.00;
        mefective_plastic_strain = 0.00;
        mfriction_angle          = 0.00; //mmaxfriction_angle; 
        mdilatancy_angle         = 0.00; //mmaxdilatancy_angle; 
        mlocal_fail_factor       = 0.00; 
        
        noalias(mplastic_strain)         = ZeroVector(4);
        noalias(mcurrent_plastic_strain) = ZeroVector(4);  

        double Gc           = (*mpProperties)[CRUSHING_ENERGY]/mlength;
        double length_limit = 2.00*mE*Gc/((*mpProperties)[FC]*(*mpProperties)[FC]);       

        if (length_limit<mlength) {std::cout<<"Element length greater than permitted"<<std::endl;}
           
        mpFluencyCriteria->InitializeMaterial(*mpProperties);
        mpFluencyCriteria_Traction->InitializeMaterial(*mpProperties);


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

        //uptating
        //KRATOS_WATCH("FINALIZE_FINALIZE_FINALIZE_FINALIZE_FINALIZE_FINALIZE_FINALIZE")
        array_1d<double,1> Variables;
        mefective_plastic_strain =  mcurrent_efective_plastic_strain;
        noalias(mplastic_strain) =  mcurrent_plastic_strain;
        mplastic_damage          =  mcurrent_plastic_damage;
        mdisipation              =  mcurrent_disipation;  
        mcohesion                =  mcurrent_cohesion;  
        mfriction_angle          =  mcurrent_friction_angle; 
        mdilatancy_angle         =  mcurrent_dilatancy_angle;    
        noalias(mFt)             =  mcurrent_Ft;
        

         // updating internal variables
	Variables[0] = mcurrent_cohesion;
	//Variables[1] = mcurrent_friction_angle;
	//Variables[2] = mcurrent_dilatancy_angle;

	mpFluencyCriteria->UpdateVariables(Variables); 
        //KRATOS_WATCH(mFt)  
        mpFluencyCriteria_Traction->UpdateVariables(mFt); 
        
        
        Matrix PlasticTensor = ZeroMatrix(3,3);  
        Vector Aux_P           = ZeroVector(3); 
        Matrix Aux_V = ZeroMatrix(3,3);   
	PlasticTensor(0,0)   = mplastic_strain[0];     PlasticTensor(0,1) = 0.5*mplastic_strain[2]; PlasticTensor(0,2) = 0.00;
	PlasticTensor(1,0)   = 0.5*mplastic_strain[2]; PlasticTensor(1,1) = mplastic_strain[1];     PlasticTensor(1,2) = 0.00;
	PlasticTensor(2,0)   = 0.00;                   PlasticTensor(2,1) = 0.00;                   PlasticTensor(2,2) = mplastic_strain[3];    

        //if(PlasticTensor(0,0)==0) {PlasticTensor(0,0) = 1E-9;}  
        double Euc        =  2.00 * (*mpProperties)[FRACTURE_ENERGY]/ ( (*mpProperties)[FT] * mlength);   
        SD_MathUtils<double>::EigenVectors(PlasticTensor, Aux_V, Aux_P, 1E-9, 100);
        double Ef = (*max_element(Aux_P.begin(), Aux_P.end()));   
        mlocal_fail_factor =  Ef / Euc; 
        //if (mlocal_fail_factor > 1.00) {mlocal_fail_factor = 1.00; }  

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
/*boost::numeric::ublas::bounded_matrix<double,4,4> C;
CalculateElasticMatrix(C); 
noalias(Stress) = prod(C, Strain);
KRATOS_WATCH(Stress)
*/
/// Owen
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

void PlasticDamage2D::CalculateMaterialResponse(
            const Vector& StrainVector,
            Vector& StressVector,
            Matrix& algorithmicTangent,
            bool calculate_stress_flag,
            bool calculate_tangent_flag,
            bool save_internal_variables
            )
{
if (calculate_stress_flag==true) { CalculateStress(StrainVector, StressVector);}
if (calculate_tangent_flag==true){CalculateStressAndTangentMatrix(StressVector,StrainVector, algorithmicTangent);}
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
        double ElasticDomain_2               = 0.00;
        double toler                         = 1E-6;
	array_1d<double,4> StrainVector_aux  = ZeroVector(4);          
	array_1d<double,4> ElasticStrain     = ZeroVector(4);                                           
	array_1d<double,4> ElasticStress     = ZeroVector(4);     
	array_1d<double,3> Variables         = ZeroVector(3); 
        StressVector                         = ZeroVector(3); 
 
        
	/// calculating elastic strain
        /// plane strain and axial symmetric 
	StrainVector_aux[0] =  StrainVector[0];
	StrainVector_aux[1] =  StrainVector[1];
	StrainVector_aux[2] =  StrainVector[2];
	StrainVector_aux[3] =  mcurrent_plastic_strain[3];   // et = ep; ez elastico siempre tiene que ser igual a cero 

        /// The elastic strain
	noalias(ElasticStrain) = StrainVector_aux - mcurrent_plastic_strain;

	/// calculating elastic stress
	CalculateElasticStress(ElasticStrain, ElasticStress);


        /// Comprobado criterio de fluencia   
        mpFluencyCriteria->CalculateEquivalentUniaxialStress(ElasticStress, ElasticDomain_1); 
        mpFluencyCriteria_Traction->CalculateEquivalentUniaxialStress(ElasticStress, ElasticDomain_2);        


        /// Spectral Descomposition 
        Compute_Principal_Stress(ElasticStress);
        array_1d<double, 3 > PrincipalStress;  
        array_1d<unsigned int, 3 > Order;

	/// Guardando el orden de los cambios
	IdentifyMaximunAndMinumumPrincipalStres_CalculateOrder(PrincipalStress , Order);
	  
        /// Casos 
        /// Traccion Pura             
        if( (PrincipalStress[0]>=0)  &&  (PrincipalStress[1]>=0) && (PrincipalStress[2]>=0) )
	{
	  mCase = Tracc_Pura;
	}  
        /// Compresion Pura 
        else if( (PrincipalStress[0]<=0)  &&  (PrincipalStress[1]<=0) && (PrincipalStress[2]<=0) )
	{
	  mCase =  Tracc_Pura; //Comp_Pura;
	}   
	/// Mixto 
        /// El hecho de que el estado sea mixto no indica que las funciones de ambas suerficies sean activas  
       else { 
             mCase = Tracc_Pura;
            }
        
                     
        switch(mCase)
	      {

              /// None Todo es Elastico 
	      case None:
               {       
                 break;
	       }

	      case Tracc_Pura:
	      {

              /// Caso elastico 
	      if( (ElasticDomain_2 <= toler ) )
		{
		  noalias(mcurrent_Ft)             = mFt;
		  mcurrent_efective_plastic_strain = mefective_plastic_strain;
		  noalias(mcurrent_plastic_strain) = mplastic_strain;
		  mcurrent_plastic_damage          = mplastic_damage;
		  mcurrent_disipation              = mdisipation;
		  mcurrent_cohesion                = mcohesion;
		  mcurrent_friction_angle          = mfriction_angle;
		  mcurrent_dilatancy_angle         = mdilatancy_angle;  
		  mComputeTangentMatrix            = false;
		}


             else
             {
              mComputeTangentMatrix = false; 
              array_1d<double,3> Sigma = ZeroVector(3);           
              array_1d<double,4> Aux_Elastic_Stress = ElasticStress; 

              mactive_surface.resize(0, false); 
	      mactive_surface.reserve(5); 
	      if((mpFluencyCriteria_Traction->mMultisurface_Platicity_Yield[0]) > toler ) {mactive_surface.push_back(0); }
	      if((mpFluencyCriteria_Traction->mMultisurface_Platicity_Yield[1]) > toler ) {mactive_surface.push_back(1); }      
	      if((mpFluencyCriteria_Traction->mMultisurface_Platicity_Yield[2]) > toler ) {mactive_surface.push_back(2); }   

                                 

              double norma      = 1.00;      
              unsigned int iter = 0;    
              double G  = 0.5*mE / (1.00 + mNU);
              double K  = mE / (3.00 * (1.00-2.00*mNU) );
              double ft = (*mpProperties)[FT];
	      double H  =  mlength * ft * ft / ( 2.00 * (*mpProperties)[FRACTURE_ENERGY] );   
	      double delta_gamma_a = 0.00; 
	      double delta_gamma_b = 0.00; 
	      double delta_gamma_c = 0.00; 
	      Matrix d;  
	      Matrix d_inv;
	      Vector delta_gamma;
	      Vector residual;     

	      /// WARNING = Si el hablandamiento es no lineal usar newton Rapshon.  
	      /// Una superficie activa
	      if(mactive_surface.size()==1)
	      {
                 iter  = 0;     
                 norma = 1.00;                
                 delta_gamma     = ZeroVector(1); 
                 residual        = ZeroVector(1); 
                 unsigned int& pos       = mactive_surface[0];                                            

                 while(iter++<=100 && norma>= toler) 
	          {  
                    if(iter>=100){KRATOS_WATCH("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA")}   
                    delta_gamma[0]  += (mpFluencyCriteria_Traction->mMultisurface_Platicity_Yield[pos]) / (4.00 * G /3.00 + K - H ); 
                    //if(delta_gamma[0] < 0.00) {delta_gamma[0] = 0.00; }   
 
                     /// Updatinf mFt
                    if(mcurrent_Ft[pos] <= 0.00) { mcurrent_Ft[pos] = mFt[pos];}
                    else {mcurrent_Ft[pos] = mFt[pos] - H * delta_gamma[0];}                          

 
                    /// comprobando si mft se cumplio   
                                   /// comprobando si mft se cumplio   
		    if(mcurrent_Ft[pos] <= 0.00) 
                       {
		           mcurrent_Ft[pos] = 0.00;
                           delta_gamma_a  = delta_gamma[0];  
                           mpFluencyCriteria_Traction->UpdateVariables(mcurrent_Ft);   
                           break;         
                       } 

                    else
                    {
                    mpFluencyCriteria_Traction->UpdateVariables(mcurrent_Ft);  
                    delta_gamma_a   = delta_gamma[0];                        
                    mpFluencyCriteria_Traction->CalculateEquivalentUniaxialStress(ElasticStress, norma); 
                    residual[0] =  norma - delta_gamma_a * (4.00  * G / 3.00 + K ) ;
                    norma    = fabs(residual[0]);   
                    }
 
                  }   
                     /// Updating Stress  
		    if(mcurrent_Ft[0] == 0.0) 
                    {
		    Sigma[0] = 0.00;
		    Sigma[1] = PrincipalStress[1] - delta_gamma_a*(-2.00 * G / 3.00 + K ); 
		    Sigma[2] = PrincipalStress[2] - delta_gamma_a*(-2.00 * G / 3.00 + K );  
		    } 
                    else
                    { 
	            Sigma[0] = PrincipalStress[0] - delta_gamma_a*(4.00  * G / 3.00 + K );    
		    Sigma[1] = PrincipalStress[1] - delta_gamma_a*(-2.00 * G / 3.00 + K ); 
		    Sigma[2] = PrincipalStress[2] - delta_gamma_a*(-2.00 * G / 3.00 + K );    
                    }      
              } 

	      /// dos superficies activas    
	      if(mactive_surface.size()==2)
	      {
              int singular         =  0;    
	      delta_gamma = ZeroVector(2);
	      residual    = ZeroVector(2);   
              iter  = 0; 
              norma = 1.00;        
              residual[0] = mpFluencyCriteria_Traction->mMultisurface_Platicity_Yield[0];
              residual[1] = mpFluencyCriteria_Traction->mMultisurface_Platicity_Yield[1]; 
              std::cout << std::fixed << std::setprecision(10); 
              KRATOS_WATCH( residual ) 
              KRATOS_WATCH( mcurrent_Ft )
              KRATOS_WATCH(ElasticStress)

              while(iter++<=100 && norma>= toler) 
		  {
   
		      d.resize(2,2);
		      d_inv.resize(2,2);
		      d(0,0) = -( 4.00 * G / 3.00 + K ) + H;  d(0,1) = -(-2.00 * G / 3.00 + K );  
		      d(1,0) = -(-2.00 * G / 3.00 + K );      d(1,1) = -( 4.00 * G / 3.00 + K ) + H;  

		      singular             =  SD_MathUtils<double>::InvertMatrix(d, d_inv);
		      noalias(delta_gamma) =  delta_gamma - Vector(prod(d_inv, residual)); 
                      if(iter>=100){
                        KRATOS_WATCH("BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB")
                        KRATOS_WATCH(residual) 
                        }    



		      //if(delta_gamma[0] < 0.00) {delta_gamma[0] = 0.00; }  
		      //if(delta_gamma[1] < 0.00) {delta_gamma[1] = 0.00; }     

		      delta_gamma_a = delta_gamma[0];
		      delta_gamma_b = delta_gamma[1];
                     
                       /// Updatinf mFt
                       mcurrent_Ft[0] = mFt[0] - H * delta_gamma[0];   
                       mcurrent_Ft[1] = mFt[1] - H * delta_gamma[1]; 
                       
                       /// comprobando si mft se cumplio   
		       if(mcurrent_Ft[0] <= 0.00) {mcurrent_Ft[0] = 0.00; }
                       if(mcurrent_Ft[1] <= 0.00) {mcurrent_Ft[1] = 0.00; } 
                        

                       mpFluencyCriteria_Traction->UpdateVariables(mcurrent_Ft);                      
                       mpFluencyCriteria_Traction->CalculateEquivalentUniaxialStress(ElasticStress, norma); 
                               
                  
                       residual[0] = mpFluencyCriteria_Traction->mMultisurface_Platicity_Yield[0];
                       residual[1] = mpFluencyCriteria_Traction->mMultisurface_Platicity_Yield[1];
                       
                       residual[0]=  residual[0] - delta_gamma_a*( 4.00  * G / 3.00 + K ) - delta_gamma_b*( -2.00  * G / 3.00 + K );   
                       residual[1]=  residual[1] - delta_gamma_a*(-2.00 * G  / 3.00 + K ) - delta_gamma_b*( 4.00  * G / 3.00 + K );   
                       //KRATOS_WATCH(residual);
  
                       //KRATOS_WATCH(residual); 
                       norma = norm_2(residual);

                       }

		      /// Updating Stress
                     /// Updating Stress  
                     Sigma[0] = PrincipalStress[0] - delta_gamma_a*( 4.00  * G / 3.00 + K ) - delta_gamma_b*( -2.00  * G / 3.00 + K );  
                     if(mcurrent_Ft[0] <= 0.00) {Sigma[0] = 1E-14;} 
                     Sigma[1] = PrincipalStress[1] - delta_gamma_a*(-2.00 * G  / 3.00 + K ) - delta_gamma_b*( 4.00  * G / 3.00 + K );   
		     if( mcurrent_Ft[1]<= 0.0)  {Sigma[0] = 1E-14;}  
		     Sigma[2] = PrincipalStress[2] - (delta_gamma_a + delta_gamma_b) * (-2.00 * G / 3.00 + K );    

                    
	      } 

	      /// Muy poco probable en 2D
	      if(mactive_surface.size()==3)
	      {

	      delta_gamma = ZeroVector(3);
	      residual    = ZeroVector(3);   
	      residual[0] = mpFluencyCriteria_Traction->mMultisurface_Platicity_Yield[0];
	      residual[1] = mpFluencyCriteria_Traction->mMultisurface_Platicity_Yield[1];   
	      residual[2] = mpFluencyCriteria_Traction->mMultisurface_Platicity_Yield[2];   
	      d.resize(3,3);
	      d_inv.resize(3,3);
	      d(0,0) = -( 4.00 * G / 3.00 + K ) + H;    d(0,1) = -(-2.00 * G / 3.00 + K );         d(0,2) = -(-2.00 * G / 3.00 + K ); 
	      d(1,0) = -(-2.00 * G / 3.00 + K );        d(1,1) = -( 4.00 * G / 3.00 + K ) + H;     d(1,2) = -(-2.00 * G / 3.00 + K );   
	      d(2,0) = -(-2.00 * G / 3.00 + K );        d(2,1) = -(-2.00 * G / 3.00 + K );         d(2,2) = -( 4.00 * G / 3.00 + K ) + H;   


	      int singular         =  0; 
	      singular             =  SD_MathUtils<double>::InvertMatrix(d, d_inv);
	      noalias(delta_gamma) =  delta_gamma - Vector(prod(d_inv, residual)); 

              if(delta_gamma[0] < 0.00) {delta_gamma[0] = 0.00; }  
              if(delta_gamma[1] < 0.00) {delta_gamma[1] = 0.00; }      
              if(delta_gamma[2] < 0.00) {delta_gamma[2] = 0.00; }  
                  

	      delta_gamma_a = delta_gamma[0];
	      delta_gamma_b = delta_gamma[1];
	      delta_gamma_c = delta_gamma[2];

              /// Updating Stress 
              Sigma[0] = PrincipalStress[0] - delta_gamma_a*( 4.00  * G / 3.00 + K ) - (delta_gamma_b + delta_gamma_c) *( -2.00  * G / 3.00 + K );    
              Sigma[1] = PrincipalStress[1] - (delta_gamma_a + delta_gamma_c) * (-2.00 * G  / 3.00 + K ) - delta_gamma_b*( 4.00  * G / 3.00 + K ); 
	      Sigma[2] = PrincipalStress[2] - (delta_gamma_a + delta_gamma_b) * (-2.00 * G / 3.00 + K )  - delta_gamma_c*( 4.00  * G / 3.00 + K );  

	      } 

               /*
	      /// Spectral Descomposition 
	      Compute_Principal_Stress(ElasticStress);
	      //array_1d<double, 3 > PrincipalStress;  

	      /// Guardando el orden de los cambios
	      IdentifyMaximunAndMinumumPrincipalStres_CalculateOrder(PrincipalStress , Order);
               */
              
	      AssembleUpdateStressAndStrainTensor(Sigma, StrainVector_aux, Order, ElasticStrain, Aux_Elastic_Stress); 
	      noalias(ElasticStress) = Aux_Elastic_Stress;
              

	      /// La cohesion disminuiria en teoria
	      /// General evoluation of accumulated hardening  
	      mcurrent_efective_plastic_strain = mefective_plastic_strain + norm_1(delta_gamma); 

	      /// aculated Von misses plastic
              //double aux_var = sqrt(delta_gamma_a*delta_gamma_a + delta_gamma_b * delta_gamma_b + delta_gamma_c * delta_gamma_c);
	      //mcurrent_efective_plastic_strain = mefective_plastic_strain + (2.00 / sqrt(3.00) ) * aux_var;
 

	      /// Update Cohesion and H
	      ///WARNING = Using linear softening for Cohesion
              /*
	      double cohesion0  = (*mpProperties)[FC]/(2.00*tan(mmaxfriction_angle/2.00 + PI/4.00));
	      double Euc        =  2.00 * (*mpProperties)[FRACTURE_ENERGY]/ ( (*mpProperties)[FT] * mlength) ;
	      mcurrent_cohesion = cohesion0 * ( 1.00 - mcurrent_efective_plastic_strain / Euc);  
	      Vector Variables(1); Variables[0] = mcurrent_cohesion;
	      if(mcurrent_cohesion < 0.0) {mcurrent_cohesion = 0.00;}
	      mpFluencyCriteria->UpdateVariables(Variables);   
              */ 

              //mpFluencyCriteria_Traction->CalculateEquivalentUniaxialStress(ElasticStress, ElasticDomain_2);   
              //KRATOS_WATCH("--------------------------------------------------------------");          

	      }
            break;
           }
              case Comp_Pura:
		      {
                       /// Caso elastico 
	               if( (ElasticDomain_1 <= toler ) )
		        {
			  noalias(mcurrent_Ft)             = mFt;
			  mcurrent_efective_plastic_strain = mefective_plastic_strain;
			  noalias(mcurrent_plastic_strain) = mplastic_strain;
			  mcurrent_plastic_damage          = mplastic_damage;
			  mcurrent_disipation              = mdisipation;
			  mcurrent_cohesion                = mcohesion;
			  mcurrent_friction_angle          = mfriction_angle;
			  mcurrent_dilatancy_angle         = mdilatancy_angle;  
			  mComputeTangentMatrix            = false;
		        }
                      else
                      {

                       mComputeTangentMatrix = false; 
                       array_1d<double,3> Sigma = ZeroVector(3);           
                       array_1d<double,4> Aux_Elastic_Stress = ElasticStress; 
                     
		      /// identufy possible edge return:either right or left of the main plain             
		      const bool edges = ReturnToEdges(Aux_Elastic_Stress);

		      /// One Vector retutn mapping to main plane 
		      ReturnMappingToMainPlane(ElasticStress, PrincipalStress, Order,  Sigma);

		      ///check the validity of main plane return
		      bool check = CheckValidity(Sigma); 
		      if(check==true)
		      {
		      AssembleUpdateStressAndStrainTensor(Sigma, StrainVector_aux, Order, ElasticStrain, Aux_Elastic_Stress); 
		      noalias(ElasticStress) = Aux_Elastic_Stress; 
		      }               
		      /// Return to the edges
		      else
		      {

		      noalias(Aux_Elastic_Stress) = ElasticStress;  
		      //edges = ReturnToEdges(Aux_Elastic_Stress);
		      //KRATOS_WATCH(Sigma)  
		      TwoVectorReturnToEdges(Aux_Elastic_Stress, PrincipalStress, Sigma,  Order, edges);
		      check = CheckValidity(Sigma);
		      /// Solo ira al apex si todas las tensiones son positivas
		      //if(Sigma[0] <=  0.00 or Sigma[1] <= 0.00 or Sigma[2] <=  0.00) { check = true;}                                    
		      //if(check==true) 
		      {

		      AssembleUpdateStressAndStrainTensor(Sigma, StrainVector_aux, Order, ElasticStrain, Aux_Elastic_Stress); 
		      noalias(ElasticStress) = Aux_Elastic_Stress;               
		      }
		      }
                      mpFluencyCriteria->CalculateEquivalentUniaxialStress(ElasticStress, ElasticDomain_1); 
                      }
		      break; 
		    }

              case Mixto:
                         {
                           /// Caso elastico 
	             if( (ElasticDomain_1 <= toler && ElasticDomain_2<= toler ) )
		            {
			      noalias(mcurrent_Ft)             = mFt;
			      mcurrent_efective_plastic_strain = mefective_plastic_strain;
			      noalias(mcurrent_plastic_strain) = mplastic_strain;
			      mcurrent_plastic_damage          = mplastic_damage;
			      mcurrent_disipation              = mdisipation;
			      mcurrent_cohesion                = mcohesion;
			      mcurrent_friction_angle          = mfriction_angle;
			      mcurrent_dilatancy_angle         = mdilatancy_angle;  
			      mComputeTangentMatrix            = false;
		          }
                      else
                          {

                          
                          array_1d<double,3> Sigma = ZeroVector(3);   
                          array_1d<double,4> Aux_Elastic_Stress = ZeroVector(3);   
                          CombinedRankineMohrCoulombSurfaces(ElasticDomain_1, ElasticDomain_2, PrincipalStress, Sigma); 
                          AssembleUpdateStressAndStrainTensor(Sigma, StrainVector_aux, Order, ElasticStrain, Aux_Elastic_Stress); 
		          noalias(ElasticStress) = Aux_Elastic_Stress; 
                          Tensile_Fracture_Model(ElasticStress,  mcurrent_plastic_strain);
                          mpFluencyCriteria_Traction->UpdateVariables(mcurrent_Ft);  

                           }
       
                       break;                  
                       }

                }             

            //KRATOS_WATCH(ElasticStress) 
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

ElasticStrain[0] = (ElasticStress[0] - p ) / (2.00 * G) + p/(3.00 * K);          //StrainTensor(0,0);   //xx
ElasticStrain[1] = (ElasticStress[1] - p ) / (2.00 * G) + p/(3.00 * K);       //StrainTensor(1,1);      //yy
ElasticStrain[2] =  ElasticStress[2]/G;      //2.00*StrainTensor(0,1);                                  //xy
ElasticStrain[3] = (ElasticStress[3] - p )/  (2.00 * G) + p/(3.00 * K);       //StrainTensor(2,2);      //zz  
noalias(mcurrent_plastic_strain) = StrainVector_aux - ElasticStrain; 

Updated_Internal_Variables(ElasticStress, mcurrent_plastic_strain);
Tensile_Fracture_Model(ElasticStress,  mcurrent_plastic_strain);


}

//***********************************************************************************************
//***********************************************************************************************


 bool PlasticDamage2D::CheckValidity(const array_1d<double,3>& Sigma)
{
  
  bool check   = false;
  array_1d<double,3> Aux_Sigma;
  Aux_Sigma[0] = fabs(Sigma[0]);
  Aux_Sigma[1] = fabs(Sigma[1]);
  Aux_Sigma[2] = fabs(Sigma[2]);
  
  double delta = (*max_element(Aux_Sigma.begin(), Aux_Sigma.end())) * 1.00E-6; 
  /* 
  /// WARNING= Tolerancias pueden llevarte al apex
  double toler = fabs(Sigma[1]-Sigma[0]); ///  Dos tensiones se igualan 
  if(toler<0.001) {Aux_Sigma[1] = Aux_Sigma[0];}

  toler = fabs(Sigma[0]-Sigma[2]); /// Dos tensiones se igualan 
  if(toler<0.001) {Aux_Sigma[0] = Aux_Sigma[2];}

  toler = fabs(Sigma[2]-Sigma[1]); /// Dos tensiones se igualan 
  if(toler<0.001) {Aux_Sigma[1] = Aux_Sigma[2];}
  */
     if( (Sigma[0] + delta) >= Sigma[1] && (Sigma[1] + delta ) >= Sigma[2]){ check = true;}
  return check;
}
        
//***********************************************************************************************
//***********************************************************************************************
void PlasticDamage2D::ReturnMappingToMainPlane(const array_1d<double,4>& ElasticStress,
const array_1d<double,3>& PrincipalStress, 
array_1d<unsigned int,3>& order, array_1d<double,3>& Sigma)
{


double deltagamma = 0.00;
double residual   = 0.00;
double G          = 0.5*mE / (1.00 + mNU);
double K          = mE / (3.00 * (1.00-2.00*mNU) );

double sinphi     = sin(mmaxdilatancy_angle);
double sintetha   = sin(mmaxfriction_angle);
double costetha   = cos(mmaxfriction_angle);

double   toler    = 1E-14;
unsigned iter     = 0;
unsigned max      = 1000;
double   d        = 0.00;

double n     = (*mpProperties)[FC] / (*mpProperties)[FT];
double num   = n *  (*mpProperties)[FT] * (*mpProperties)[FT] * mlength ;
double tan2  = tan(mmaxfriction_angle/2.00 + PI/4.00 ); tan2 = tan2 * tan2; 
double denom = 4.00 * tan2 * (*mpProperties)[FRACTURE_ENERGY] ;

double H     =  -num / denom;

mactive_surface.resize(1);
mactive_surface[0] = 0; 

residual = mpFluencyCriteria->mMultisurface_Platicity_Yield[0]; 

/// WARNING  Aculutated von Misses or Normal
double fact = 2.00 * costetha; 
//double fact   = 1.1547005389792 * sqrt(sinphi * sinphi  + 1.00 );   


while(fabs(residual)>toler && iter++ < max )
   {      
    d        = -4.00 * G * ( 1.00 + (sinphi * sintetha)/3.00 ) - 4.00 * K * sinphi * sintetha - 2.00 *costetha * H * fact ;
    deltagamma     += -residual/d;

    /// check for convergence
    /// General evoluation of accumulated hardening  
    mcurrent_efective_plastic_strain = mefective_plastic_strain + 2.00 * costetha * deltagamma;      
     
    /// aculated Von misses plastic
   // mcurrent_efective_plastic_strain = mefective_plastic_strain + (2.00 * sqrt( (sinphi * sinphi + 1.00)/3.00 ) * deltagamma);


    /// Update Cohesion and H
    ///WARNING = Using linear softening for Cohesion
    
    double cohesion0  = (*mpProperties)[FC]/(2.00*tan(mmaxfriction_angle/2.00 + PI/4.00));
    double Euc        =  2.00 * (*mpProperties)[FRACTURE_ENERGY]/ ( (*mpProperties)[FT] * mlength);    
    mcurrent_cohesion = cohesion0 * ( 1.00 - mcurrent_efective_plastic_strain / Euc);  
    Vector Variables(1); Variables[0] = mcurrent_cohesion;
    if(mcurrent_cohesion < 0.0) {
    mcurrent_cohesion = 0.00;
    mpFluencyCriteria->UpdateVariables(Variables); 
    
    Sigma[0] = 1E-9 ; order[0] = 0; 
    Sigma[1] = 1E-9 ; order[1] = 1; 
    Sigma[2] = 1E-9 ; order[2] = 2; 
    break; }
    

    /// Check for convergence
    mpFluencyCriteria->CalculateEquivalentUniaxialStress( ElasticStress, residual);       
    //residual =mpFluencyCriteria->mMultisurface_Platicity_Yield[0];
    residual   -= (4.00 * G * ( 1.00 + (sinphi * sintetha)/3.00 ) + 4.00 * K * sinphi * sintetha) * deltagamma;
    

    /// No puede ser nunca negativo  
    if(deltagamma < 0){ 
    KRATOS_WATCH(deltagamma)
    KRATOS_WATCH("GAMMA_GAMMA_GAMMA_GAMMA_GAMMA_GAMMA_GAMMA_GAMMA_GAMMA_GAMMA" ) }  
   }
   
   /*
    /// Spectral Descomposition 
   Compute_Principal_Stress(ElasticStress);
   array_1d<double, 3 > PrincipalStress;  
   /// Guardando el orden de los cambios
   
   IdentifyMaximunAndMinumumPrincipalStres_CalculateOrder(PrincipalStress , order);
   */
   

   Sigma[0] = PrincipalStress[0] - ( 2.00 * G * (1.00 + sinphi/3.00) + 2.00 * K * sinphi) * deltagamma; 
   Sigma[1] = PrincipalStress[1] + (4.00 * G /3.00 - 2.00 * K) * sinphi  * deltagamma; 
   Sigma[2] = PrincipalStress[2] +  ( 2.00 * G * (1.00 - sinphi/3.00) - 2.00 * K * sinphi) * deltagamma; 
   
}

//***********************************************************************************************
//***********************************************************************************************

void PlasticDamage2D::IdentifyMaximunAndMinumumPrincipalStres_CalculateOrder(array_1d<double, 3 >& PrincipalStress ,array_1d<unsigned int,3>& order)
{
 

   //noalias(PrincipalStress) = mPrincipalStress; 
   /// sorting
   /*sort (PrincipalStress.begin(), PrincipalStress.end()); 
   reverse(PrincipalStress.begin(),PrincipalStress.end());
   fill(order.begin(),order.end(),2);
   for (unsigned int i = 0; i<3; i++){
       for(unsigned int j = 0; j<3; j++)
          {  
	    if( (mPrincipalStress[j] == PrincipalStress[i]) )       //&& (j != order[j]))
	    {  
	      bool repeticion = false;
              for(unsigned int k = 0 ; k < i; k++) { if (j==order[k] ) {repeticion = true;  break;} } 
              if(repeticion==false) {order[i] = j;}
              //order[i] = j;
              //break; 
	    }                    
           }
        }
   */

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


bool PlasticDamage2D::ReturnToEdges(const array_1d<double,4>& ElasticStress)
{

    double sinphi       = sin(mmaxdilatancy_angle); 
    bool   return_rigth = false; /// left edges

    /// Spectral Descomposition 
    Compute_Principal_Stress(ElasticStress);
    array_1d<double, 3 > PrincipalStress = mPrincipalStress; 

    /// sorting
    sort (PrincipalStress.begin(), PrincipalStress.end()); 
    reverse(PrincipalStress.begin(),PrincipalStress.end());
    
    double cond = (1.00 - sinphi) * PrincipalStress[0] - 2.00 * PrincipalStress[1] + (1.00 + sinphi) * PrincipalStress[2];  
    if(cond > 0.00) {return_rigth = true;} /// rigth edges

    return return_rigth;  
   
}

//***********************************************************************************************
//***********************************************************************************************


void PlasticDamage2D::TwoVectorReturnToEdges(const array_1d<double,4>& ElasticStress, const array_1d<double,3>& PrincipalStress, 
array_1d<double,3>& Sigma, array_1d<unsigned int,3>& order, const bool& edges )
{
    
  double sinphi       = sin(mmaxdilatancy_angle);
  //double cosphi       = cos(mmaxdilatancy_angle);
  double sintetha     = sin(mmaxfriction_angle);
  double costetha     = cos(mmaxfriction_angle);
  double G            = 0.5*mE / (1.00 + mNU);
  double K            = mE / (3.00 * (1.00-2.00*mNU) );
  double norma        = 1.00;

  Matrix d            = ZeroMatrix(2,2);
  Matrix d_inv        = ZeroMatrix(2,2);
  array_1d<double,2> delta_gamma = ZeroVector(2);
  array_1d<double,2> residual; 
 
  mactive_surface.resize(2);
  mactive_surface[0] = 0; 
  mactive_surface[1] = 5;

  residual[0] = mpFluencyCriteria->mMultisurface_Platicity_Yield[0]; 
  residual[1] = mpFluencyCriteria->mMultisurface_Platicity_Yield[5];

  ///left edges
  if(edges==false) 
   { residual[1] = mpFluencyCriteria->mMultisurface_Platicity_Yield[1]; 
     mactive_surface[1] = 1; 
   }   
   
  //KRATOS_WATCH(residual[1])
  //KRATOS_WATCH(residual[1])  


  double   fact       = 0.00;
  double   toler      = 1E-9;
  unsigned iter       = 0;
  unsigned max        = 1000;

int singular  = 0.00;
double aux1   = 0.00, aux2 = 0.00, aux3 = 0.00; 
double result = 0; 

  while(norma > toler && iter++ <max)
{
  //WARNING
double n     = (*mpProperties)[FC] / (*mpProperties)[FT];
double num   = n *  (*mpProperties)[FT] * (*mpProperties)[FT] * mlength ;
double tan2  = tan(mmaxfriction_angle/2.00 + PI/4.00 ); tan2 = tan2 * tan2; 
double denom = 4.00 * tan2 * (*mpProperties)[FRACTURE_ENERGY] ;

 double H     =  -num / denom;  

/// WARNING  Aculutated von Misses or Normal
  double fact_1 = 2.00 * costetha; 
  fact = 2.00 * H * costetha * fact_1;
  


  double a = 4.00 * G * ( 1.00 + (sinphi * sintetha)/3.00 ) + 4.00 * K * sinphi * sintetha;
  double b = 2.00 * G * ( 1.00 + sinphi + sintetha -(sinphi *sintetha)/3.00) + 4.00 * K * sinphi * sintetha; 
  
  /// left edges
  if(edges==false) {b = 2.00 * G * ( 1.00 - sinphi - sintetha -(sinphi *sintetha)/3.00) + 4.00 * K * sinphi * sintetha;}
  
  d(0,0) = -a - fact;   d(0,1) = -b - fact;
  d(1,0) = -b - fact;   d(1,1) = -a - fact;  
     

  singular             =  SD_MathUtils<double>::InvertMatrix(d, d_inv);
  noalias(delta_gamma) =  delta_gamma - Vector(prod(d_inv, residual)); 
 
 /// computing acumulated strain
 mcurrent_efective_plastic_strain = mefective_plastic_strain + 2.00 * costetha *(delta_gamma[0] + delta_gamma[1]); 

 /// aculated Von misses plastic
 /*double aux_1 = 2.00 * (sinphi * sinphi + 1.00 ) * (delta_gamma[0]*delta_gamma[0] + delta_gamma[0] * delta_gamma[1] + delta_gamma[1]*delta_gamma[1]);
 double aux_2 = aux_1 + 4.00 * sinphi * delta_gamma[0] * delta_gamma[1];
 if(edges==false) {aux_2 = aux_1 - 4.00 * sinphi * delta_gamma[0] * delta_gamma[1]; }          
 mcurrent_efective_plastic_strain = mefective_plastic_strain +  sqrt( 2.00 * aux_2 / 3.00);
 */
 /// Update Cohesion and H
 ///WARNING
  
  double cohesion0  = (*mpProperties)[FC]/(2.00*tan(mmaxfriction_angle/2.00 + PI/4.00));
  double Euc        =  2.00 * (*mpProperties)[FRACTURE_ENERGY]/ ( (*mpProperties)[FT] * mlength);        
  mcurrent_cohesion = cohesion0 * ( 1.00 - mcurrent_efective_plastic_strain / Euc);   
  Vector Variables(1); Variables[0] = mcurrent_cohesion;
  if(mcurrent_cohesion < 0.0) {mcurrent_cohesion = 0.00;
  mpFluencyCriteria->UpdateVariables(Variables); 
    Sigma[0] = 1E-9 ; order[0] = 0; 
    Sigma[1] = 1E-9 ; order[1] = 1; 
    Sigma[2] = 1E-9 ; order[2] = 2; 
  break; }
  

 /// computing new residual and check for convergence
mpFluencyCriteria->CalculateEquivalentUniaxialStress( ElasticStress, result); 
residual[0] = (mpFluencyCriteria->mMultisurface_Platicity_Yield[mactive_surface[0]]) - a * delta_gamma[0] - b * delta_gamma[1];
residual[1] = (mpFluencyCriteria->mMultisurface_Platicity_Yield[mactive_surface[1]]) - b * delta_gamma[0] - a * delta_gamma[1]; 

norma = norm_2(residual);

/// No puede ser nunca negativo  
if(delta_gamma[0] < 0 or delta_gamma[1] < 0 or iter>=max){ 
KRATOS_WATCH(delta_gamma)
KRATOS_WATCH( (mpFluencyCriteria->mMultisurface_Platicity_Yield[mactive_surface[1]]))
KRATOS_WATCH((mpFluencyCriteria->mMultisurface_Platicity_Yield[1]))
KRATOS_WATCH((mpFluencyCriteria->mMultisurface_Platicity_Yield[5]))
KRATOS_WATCH("RESIDUAL_RESIDUAL_RESIDUAL_RESIDUAL_RESIDUAL_RESIDUAL_RESIDUAL_" ) } 

}

/*
/// Spectral Descomposition 
Compute_Principal_Stress(ElasticStress);
array_1d<double, 3 > PrincipalStress;  

IdentifyMaximunAndMinumumPrincipalStres_CalculateOrder(PrincipalStress, order);
*/

aux1 = ( 2.00 * G * (1.00 + sinphi/3.00) + 2.00 * K * sinphi);
aux2 = (4.00 * G /3.00 - 2.00 * K) * sinphi; 
aux3 = ( 2.00 * G * (1.00 - sinphi/3.00) - 2.00 * K * sinphi);


if(edges==true)
 { 
  //KRATOS_WATCH("rirght_rirght_rirght_rirght_rirght_rirght_rirght_rirght")
  Sigma[0] = PrincipalStress[0] -  aux1 * ( delta_gamma[0] + delta_gamma[1] );
  Sigma[1] = PrincipalStress[1] +  aux2 * delta_gamma[0]   + aux3 * delta_gamma[1]; 
  Sigma[2] = PrincipalStress[2] +  aux3 * delta_gamma[0]   + aux2 * delta_gamma[1]; 
 }
/// for left edges 
else
 {
  //KRATOS_WATCH("LEFT_LEFT_LEFT_LEFT_LEFT_LEFT_LEFT_LEFT_LEFT_LEFT_LEFT_LEFT")
  Sigma[0] = PrincipalStress[0] - aux1 * delta_gamma[0] + aux2 * delta_gamma[1];   
  Sigma[1] = PrincipalStress[1] + aux2 * delta_gamma[0] - aux1 * delta_gamma[1]; 
  Sigma[2] = PrincipalStress[2] + aux3 * (delta_gamma[0] + delta_gamma[1]);
 }

//KRATOS_WATCH(Sigma)
//KRATOS_WATCH(aux1 * ( delta_gamma[0] + delta_gamma[1]) )

}


void PlasticDamage2D::ReturnMappingToApex(array_1d<double,4>& ElasticStress,  array_1d<double, 3 >& Sigma)
{

double vol_strain = 0.00;
double cottehta   = 1.00/tan(mmaxfriction_angle);
double sinphi     = sin(mmaxdilatancy_angle);
//double cosphi     = cos(mmaxdilatancy_angle);
//double sintetha   = sin(mmaxfriction_angle);
double costetha   = cos(mmaxfriction_angle);
double K          = mE / (3.00 * (1.00-2.00*mNU) );

double p  = (ElasticStress[0] + ElasticStress[1] + ElasticStress[3])/3.00; 
  /// WARNING
double H  = 0.00;

double d  = 0.00;
double r  = 1.00;

double   toler      = 1E-9;
unsigned iter       = 0;
unsigned max        = 1000;

r  =  mcurrent_cohesion * cottehta  - p;  
while(r > toler && iter<max ) 
{
  d          = H * cottehta * costetha / sinphi + K;     
  vol_strain = vol_strain -r / d;


/// WARNING = Aculutaded von Misses nod efioned for apex
mcurrent_efective_plastic_strain = mefective_plastic_strain + (costetha/sinphi) * vol_strain;
  
 /// WARNING
double cohesion0  = (*mpProperties)[FC]/(2.00*tan(mmaxfriction_angle/2.00 + PI/4.00));
double Euc        =  2.00 * (*mpProperties)[FRACTURE_ENERGY]/ ( (*mpProperties)[FT] * mlength);        
mcurrent_cohesion = cohesion0 * ( 1.00 - mcurrent_efective_plastic_strain / Euc);   
Vector Variables(1); Variables[0] = mcurrent_cohesion;
if(mcurrent_cohesion < 0.0) {mcurrent_cohesion = 0.00;
mpFluencyCriteria->UpdateVariables(Variables); 
Sigma[0] = 1E-9;  
Sigma[1] = 1E-9;  
Sigma[2] = 1E-9; 
break; }
 
 p = p - K * vol_strain;  
 r = mcurrent_cohesion * cottehta  - p;  
}

ElasticStress[0] = p;
ElasticStress[1] = p;
ElasticStress[2] = 0.00;
ElasticStress[3] = p;

Sigma[0] = p;
Sigma[1] = p;
Sigma[2] = p;


}



void PlasticDamage2D::CombinedRankineMohrCoulombSurfaces(double& ElasticDomain_1,
double& ElasticDomain_2, const array_1d<double, 3>& PrincipalStress, array_1d<double, 3>& Sigma)
{

    Sigma = ZeroVector(3);            
    enum Intersection {Two, Three, Four};
    Intersection Inter; 
    double toler = 1E-6;  
 

    mactive_surface.resize(0, false); 
    mactive_surface_T.resize(0, false); 
    mactive_surface.reserve(6);
    mactive_surface_T.reserve(5);
 
    /// Calculando las superficies activas
    if((mpFluencyCriteria_Traction->mMultisurface_Platicity_Yield[0]) > toler ) {mactive_surface_T.push_back(0); }
    if((mpFluencyCriteria_Traction->mMultisurface_Platicity_Yield[1]) > toler ) {mactive_surface_T.push_back(1); }      
    if((mpFluencyCriteria_Traction->mMultisurface_Platicity_Yield[2]) > toler ) {mactive_surface_T.push_back(2); } 

    if((mpFluencyCriteria->mMultisurface_Platicity_Yield[0]) > toler ) {mactive_surface.push_back(0); }
    if((mpFluencyCriteria->mMultisurface_Platicity_Yield[1]) > toler ) {mactive_surface.push_back(1); }      
    if((mpFluencyCriteria->mMultisurface_Platicity_Yield[5]) > toler ) {mactive_surface.push_back(5); } 
  
    

    unsigned int size_C = mactive_surface.size();  
    unsigned int size_T = mactive_surface_T.size(); 
    if( (size_C==1 && size_T == 1) && (mactive_surface_T[0]==0 && mactive_surface[0]==0) ) {Inter = Two;}
    if( size_C==2 && size_T == 1 ) {Inter = Three;}
    if( size_C==2 && size_T == 2 ) {Inter = Four;}    


                           
    double G             = 0.5*mE / (1.00 + mNU);
    double K             = mE / (3.00 * (1.00-2.00*mNU) );
    double ft            = (*mpProperties)[FT];
    double Hr            =  mlength * ft * ft / ( 2.00 * (*mpProperties)[FRACTURE_ENERGY] );   
    double delta_gamma_a = 0.00; 
    double delta_gamma_b = 0.00; 
    double delta_gamma_c = 0.00; 
    //double delta_gamma_d = 0.00; 
    Matrix d;  
    Matrix d_inv;

    
    /// WARNING = podria cambiar debido a la tension aculada
    double n            = (*mpProperties)[FC] / (*mpProperties)[FT];
    double num          = n *  (*mpProperties)[FT] * (*mpProperties)[FT] * mlength ;
    double tan2         = tan(mmaxfriction_angle/2.00 + PI/4.00 ); tan2 = tan2 * tan2; 
    double denom        = 4.00 * tan2 * (*mpProperties)[FRACTURE_ENERGY] ;
    double Hc           =  -num / denom;  
    double costetha     = cos(mmaxfriction_angle);
    double sinphi       = sin(mmaxdilatancy_angle);
    double sintetha     = sin(mmaxfriction_angle);
   
    double a      = 4.00 * G * ( 1.00 + (sinphi * sintetha)/3.00 ) + 4.00 * K * sinphi * sintetha; 
    double b      = 2.00 * G * ( 1.00 + sintetha/3.00 ) + 2.00 * K * sintetha; 
    double c      = 2.00 * G * ( 1.00 + sinphi  /3.00 ) + 2.00 * K * sinphi; 
    double dd     = 4.00 * G /3.00    + K;
    
    double fact_1 = 0.00;
    double fact_2 = 0.00; 
    //double fact_3 = 0.00;
    double fact   = 2.00 * costetha; 
      
    int singular   =  0; 
     
 
    switch(Inter)
     {
          case Two:
             {
                array_1d<double,3> delta_gamma = ZeroVector(3);
                array_1d<double,3> residual    = ZeroVector(3);  
                d.resize(2,2);
                d_inv.resize(2,2);

                residual[0] = mpFluencyCriteria->mMultisurface_Platicity_Yield[0];    
                residual[1] = mpFluencyCriteria_Traction->mMultisurface_Platicity_Yield[0]; 
   
                /// WARNING acumulated von misses o ambos
                /// ambos
                fact_1  = 2.00 * costetha; 
                fact_2  = 1.00;  

                unsigned int iter = 0;
                              
                double norm       = 1.00; 
                while( iter++<=100 && norm>=toler)
                {            
		    d(0,0) = -a - fact * fact_1 * Hc;     d(0,1) = -b - fact * fact_2 * Hc;         
		    d(1,0) = -c + Hr * (1.00 + sinphi);   d(1,1) = -dd + Hr;        

		  
		    singular             =  SD_MathUtils<double>::InvertMatrix(d, d_inv);
		    noalias(delta_gamma) =  delta_gamma - Vector(prod(d_inv, residual)); 
		    if(delta_gamma[0] < 0.00) {delta_gamma[0] = 0.00; }  
		    if(delta_gamma[1] < 0.00) {delta_gamma[1] = 0.00; }     

		    delta_gamma_a = delta_gamma[0];
		    delta_gamma_b = delta_gamma[1];

		    /// Updating Stress
Sigma[0] = PrincipalStress[0] - ( 2.00 * G * (1.00 + sinphi/3.00) + 2.00 * K * sinphi) * delta_gamma_a - delta_gamma_b * (4.00 * G /3.00 + K ); 
Sigma[1] = PrincipalStress[1] + (4.00 * G /3.00 - 2.00 * K) * sinphi  * delta_gamma_a +  delta_gamma_b * (2.00 * G /3.00 - K ); 
Sigma[2] = PrincipalStress[2] + ( 2.00 * G * (1.00 - sinphi/3.00) - 2.00 * K * sinphi) * delta_gamma_a + delta_gamma_b * (2.00 * G /3.00 - K );   


		    /// acumulated ambos                
		    mcurrent_efective_plastic_strain = mefective_plastic_strain + 2.00 * costetha *  delta_gamma_a +  delta_gamma_b; 
		    /* 
		    /// aculated Von misses plastic
		    double aux_1 = 2.00 * (sinphi * sinphi + 1.00 ) * (delta_gamma[0]*delta_gamma[0] ) + 2.00 * delta_gamma[0] * delta_gamma[1] * (sinphi + 1.00) + delta_gamma[1]*delta_gamma[1] ;
		    mcurrent_efective_plastic_strain = mefective_plastic_strain +  sqrt( 2.00 * aux_1 / 3.00);
		    */ 

		    /// Update Cohesion and H
		    ///WARNING = Using linear softening for Cohesion
		    double cohesion0  = (*mpProperties)[FC]/(2.00*tan(mmaxfriction_angle/2.00 + PI/4.00));
		    double Euc        =  2.00 * (*mpProperties)[FRACTURE_ENERGY]/ ( (*mpProperties)[FT] * mlength);    
		    mcurrent_cohesion = cohesion0 * ( 1.00 - mcurrent_efective_plastic_strain / Euc);  
		    Vector Variables(1); Variables[0] = mcurrent_cohesion;
		    if(mcurrent_cohesion < 0.0) {
		    mcurrent_cohesion = 0.00;
		    mpFluencyCriteria->UpdateVariables(Variables);
		    Sigma[0] = 1E-9; 
		    Sigma[1] = 1E-9; 
		    Sigma[2] = 1E-9;
                    break;
                 }
                 
                 residual[0] = mpFluencyCriteria->mMultisurface_Platicity_Yield[0];    
                 residual[1] = mpFluencyCriteria_Traction->mMultisurface_Platicity_Yield[0];                 
                 norm = norm_2(residual);                
		}   
                break;
             }

            case Three:
             { 
                array_1d<double,3> delta_gamma = ZeroVector(3);
                array_1d<double,3> residual    = ZeroVector(3);   
                d.resize(3,3);
                d_inv.resize(3,3);

                residual[0] = mpFluencyCriteria->mMultisurface_Platicity_Yield[0]; 
                residual[1] = mpFluencyCriteria->mMultisurface_Platicity_Yield[mactive_surface[1]];    
                residual[2] = mpFluencyCriteria_Traction->mMultisurface_Platicity_Yield[0]; 
                
                 double bb = 2.00 * G * ( 1.00 + sinphi + sintetha -(sinphi *sintetha)/3.00) + 4.00 * K * sinphi * sintetha; 
  
                 /// left edges 
                 if(mactive_surface[1]==1) {bb = 2.00 * G * ( 1.00 - sinphi - sintetha -(sinphi *sintetha)/3.00) + 4.00 * K * sinphi * sintetha;}


                /// WARNING acumulated von misses o ambos
                /// ambos
                fact_1  = 2.00 * costetha; 
                fact_2  = 1.00;  

                unsigned int iter = 0;
                double norm       = 1.00;
		double aux1 = 0.00, aux2 = 0, aux3 = 0;
                double cohesion0 = 0,Euc = 0;
 
                while( iter++<=100 && norm>=toler)
                {            
		    d(0,0) = -a - fact * fact_1 * Hc;   d(0,1) = -bb - fact * fact_1 * Hc;  d(0,2) = -b - fact * fact_2 * Hc;  
                    d(1,0) = -bb - fact * fact_1 * Hc;  d(1,1) = -a - fact * fact_1 * Hc;   d(1,2) = -b - fact * fact_2 * Hc;
		    d(2,0) = -c + Hr * (1.00 + sinphi); d(2,1) = -c + Hr * (1.00 + sinphi); d(2,2) = -dd + Hr;        

 		  
		    singular             =  SD_MathUtils<double>::InvertMatrix(d, d_inv);
		    noalias(delta_gamma) =  delta_gamma - Vector(prod(d_inv, residual)); 
		    if(delta_gamma[0] < 0.00) {delta_gamma[0] = 0.00; }  
		    if(delta_gamma[1] < 0.00) {delta_gamma[1] = 0.00; }     
                    if(delta_gamma[2] < 0.00) {delta_gamma[2] = 0.00; }   

		    delta_gamma_a = delta_gamma[0];
		    delta_gamma_b = delta_gamma[1];
                    delta_gamma_c = delta_gamma[2];

		    /// Updating Stress


		   aux1 = ( 2.00 * G * (1.00 + sinphi/3.00) + 2.00 * K * sinphi);
		   aux2 = (4.00 * G /3.00 - 2.00 * K) * sinphi; 
		   aux3 = ( 2.00 * G * (1.00 - sinphi/3.00) - 2.00 * K * sinphi);


		     if(mactive_surface[1]==5)
		    { 
		      //KRATOS_WATCH("rirght_rirght_rirght_rirght_rirght_rirght_rirght_rirght")
		      Sigma[0] = PrincipalStress[0] -  aux1 * ( delta_gamma[0] + delta_gamma[1]) - delta_gamma[2] * (4.00 * G /3.00 + K ); 
		      Sigma[1] = PrincipalStress[1] +  aux2 * delta_gamma[0]   + aux3 * delta_gamma[1] + delta_gamma[2] * (2.00 * G /3.00 - K ); 
		      Sigma[2] = PrincipalStress[2] +  aux3 * delta_gamma[0]   + aux2 * delta_gamma[1] + delta_gamma[2] * (2.00 * G /3.00 - K ); 
		    }
		    /// for left edges 
		    else
		    {
		      //KRATOS_WATCH("LEFT_LEFT_LEFT_LEFT_LEFT_LEFT_LEFT_LEFT_LEFT_LEFT_LEFT_LEFT")
		      Sigma[0] = PrincipalStress[0] - aux1 * delta_gamma[0] + aux2 * delta_gamma[1]- delta_gamma[2] * (4.00 * G /3.00 + K ); 
		      Sigma[1] = PrincipalStress[1] + aux2 * delta_gamma[0] - aux1 * delta_gamma[1] + delta_gamma[2] * (2.00 * G /3.00 - K ); 
		      Sigma[2] = PrincipalStress[2] + aux3 * (delta_gamma[0] + delta_gamma[1]) + delta_gamma[2] * (2.00 * G /3.00 - K ); 
		    }




		    /// acumulated ambos                
		    mcurrent_efective_plastic_strain = mefective_plastic_strain + 2.00 * costetha *  ( delta_gamma_a + delta_gamma_b) +  delta_gamma_b; 
		    /* 
		    /// aculated Von misses plastic
		    double aux_1 = 2.00 * (sinphi * sinphi + 1.00 ) * (delta_gamma[0]*delta_gamma[0] ) + 2.00 * delta_gamma[0] * delta_gamma[1] * (sinphi + 1.00) + delta_gamma[1]*delta_gamma[1] ;
		    mcurrent_efective_plastic_strain = mefective_plastic_strain +  sqrt( 2.00 * aux_1 / 3.00);
		    */ 

		    /// Update Cohesion and H
		    ///WARNING = Using linear softening for Cohesion
		    cohesion0  = (*mpProperties)[FC]/(2.00*tan(mmaxfriction_angle/2.00 + PI/4.00));
		    Euc        =  2.00 * (*mpProperties)[FRACTURE_ENERGY]/ ( (*mpProperties)[FT] * mlength);    
		    mcurrent_cohesion = cohesion0 * ( 1.00 - mcurrent_efective_plastic_strain / Euc);  
		    Vector Variables(1); Variables[0] = mcurrent_cohesion;
		    if(mcurrent_cohesion < 0.0) {
		    mcurrent_cohesion = 0.00;
		    mpFluencyCriteria->UpdateVariables(Variables);
		    Sigma[0] = 1E-9; 
		    Sigma[1] = 1E-9; 
		    Sigma[2] = 1E-9;
                    break;
                 }
                 
                 residual[0] = mpFluencyCriteria->mMultisurface_Platicity_Yield[0]; 
                 residual[1] = mpFluencyCriteria->mMultisurface_Platicity_Yield[mactive_surface[1]];    
                 residual[2] = mpFluencyCriteria_Traction->mMultisurface_Platicity_Yield[0];                 
                 norm = norm_2(residual);
  
                 KRATOS_WATCH(mpFluencyCriteria->mMultisurface_Platicity_Yield[0])
                 KRATOS_WATCH(mpFluencyCriteria->mMultisurface_Platicity_Yield[1])
		 KRATOS_WATCH(mpFluencyCriteria->mMultisurface_Platicity_Yield[2])
		 KRATOS_WATCH(mpFluencyCriteria->mMultisurface_Platicity_Yield[3])
		 KRATOS_WATCH(mpFluencyCriteria->mMultisurface_Platicity_Yield[4])
		 KRATOS_WATCH(mpFluencyCriteria->mMultisurface_Platicity_Yield[5])

                 KRATOS_WATCH(mpFluencyCriteria_Traction->mMultisurface_Platicity_Yield[0])
		 KRATOS_WATCH(mpFluencyCriteria_Traction->mMultisurface_Platicity_Yield[1])
		 KRATOS_WATCH(mpFluencyCriteria_Traction->mMultisurface_Platicity_Yield[2])  
                 KRATOS_WATCH("----------------------------------------------------------")              
		}   
                break;
               
             }
            case Four:
             { 
                KRATOS_WATCH("FOURFOURFOURFOURFOURFOURFOURFOURFOURFOURFOURFOURFOURFOURFOURFOURFOURFOURFOURFOURFOUR")  
                break;
             }
     }        
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

		noalias(mstemp) = prod(rF,S);
		noalias(msaux) = prod(mstemp,trans(rF));
		msaux *= J;

		if(rCauchy_StressVector.size() != 3)
			rCauchy_StressVector.resize(3);  
		
		rCauchy_StressVector[0] = msaux(0,0);
		rCauchy_StressVector[1] = msaux(1,1);
		rCauchy_StressVector[2] = msaux(1,2);
}


//***********************************************************************************************
//***********************************************************************************************

void PlasticDamage2D::CalculateStressAndTangentMatrix(Vector& StressVector,
                    const Vector& StrainVector,
                    Matrix& algorithmicTangent)
{     

if(mComputeTangentMatrix == false)
{
CalculateConstitutiveMatrix(StrainVector, algorithmicTangent);
}
else
{
boost::numeric::ublas::bounded_matrix<double,4,4> C;
boost::numeric::ublas::bounded_matrix<double,4,4> Caux;
 
Matrix A = ZeroMatrix(4,4); 
CalculateElasticMatrix(C);


double var;
array_1d<double, 4 > vect_one;
array_1d<double, 4 > vect_two;
array_1d<double, 4 > Stress;

///WARNING = Solo Plane Strain. No se ha aplicado aun la regla de Koiter
Stress[0] = StressVector[0];
Stress[1] = StressVector[1];
Stress[2] = StressVector[2];
Stress[3] = mNU*(StressVector[0] + StressVector[1]);

/// * Computing flow for Tangent Matrix
vector<Vector> Derivate_Potencial;
vector<Vector> Derivate_Fluency;  
mpFluencyCriteria->CalculateDerivateFluencyCriteriaMultiSurface(Stress, Derivate_Fluency);
mpFluencyCriteria->CalculateDerivatePotencialFlowCriteriaMultiSurface(Stress,  Derivate_Potencial);
Compute_Derivate(Derivate_Fluency, mDerivate_Fluency); 
Compute_Derivate(Derivate_Potencial, mDerivate_Potencial);


/// WARNING = Cohesaion
//Compute_A(Stress, mDerivate_Potencial, A);

noalias(vect_one) = prod(C, mDerivate_Fluency[0]);
noalias(vect_two) = prod(C, mDerivate_Potencial[0]);
var               = inner_prod(mDerivate_Fluency[0], vect_two);
noalias(Caux)     = outer_prod(vect_one, vect_two);

if (mcurrent_cohesion <=0.00 && mcurrent_plastic_damage == 1.00 )
{
   noalias(algorithmicTangent) = C;
}
else
{
noalias(Caux )    = Caux/(var + A(0,0));
noalias(Caux)     = C - Caux;  

for (unsigned int i = 0; i<3; i++)
 {for (unsigned int j = 0; j<3; j++)
    {
     algorithmicTangent(i,j) = Caux(i,j);

    }
}
}
}

}
                                       

//***********************************************************************************************
//***********************************************************************************************

void PlasticDamage2D::UpdateMaterial( const Vector& StrainVector,
                                     const Properties& props,
                                     const GeometryType& geom,
                                     const Vector& ShapeFunctionsValues,
                                     const ProcessInfo& CurrentProcessInfo)
{
        //uptating
        //KRATOS_WATCH("UPDATE_UPDATE_UPDATE_UPDATE_UPDATE_UPDATE_UPDATE_UPDATE")
        Vector Variables(1);
        noalias(mcurrent_plastic_strain) = mplastic_strain;
        mcurrent_plastic_damage          = mplastic_damage;
        mcurrent_disipation              = mdisipation;
        mcurrent_cohesion                = mcohesion;
        mcurrent_friction_angle          = mfriction_angle;
        mcurrent_dilatancy_angle         = mdilatancy_angle;   
        mcurrent_efective_plastic_strain = mefective_plastic_strain;
        noalias(mcurrent_Ft)             = mFt;

        //updating internal variables
	Variables[0] = mcurrent_cohesion;
	//Variables[1] = mcurrent_friction_angle;
	//Variables[2] = mcurrent_dilatancy_angle;
	mpFluencyCriteria->UpdateVariables(Variables); 
        mpFluencyCriteria_Traction->UpdateVariables(mcurrent_Ft);    
        //KRATOS_WATCH(mcurrent_Ft) 

}
     

//***********************************************************************************************
//***********************************************************************************************

/// Stress = Elatsic Stress
/// Strain The plastic Strain
void PlasticDamage2D::Updated_Internal_Variables(const Vector& Stress,  const Vector& Strain)
{
const int    iter  = 50;
const double zero  = 1.0E-9;
double teta_a      = 0.00;
double teta_b      = 0.00;
double teta        = 0.00;


Matrix StressTensor    = ZeroMatrix(3,3);
Matrix EigenVectors    = ZeroMatrix(3,3);
Vector PrincipalStress = ZeroVector(3);

Vector StressVector(3);
StressVector[0] = Stress[0];
StressVector[1] = Stress[1];
StressVector[2] = Stress[2]; 

mpFluencyCriteria->State_Tensor(StressVector,StressTensor);

SD_MathUtils<double>::EigenVectors(StressTensor, EigenVectors,PrincipalStress, zero, iter);

teta_a =  Tensor_Utils<double>::Mc_aully(PrincipalStress);
teta_b =  norm_1(PrincipalStress);  


// Condicion para no tener divison de 0/0 
if (teta_b==0.00)
{teta = 0.00;}
else
{teta = teta_a/teta_b;}

//updating plastic_damage
//computing Kp_punto
double kp_punto = 0.00; 
double gc_p   = (*mpProperties)[CRUSHING_ENERGY]/mlength;
double gf_p   = (*mpProperties)[FRACTURE_ENERGY]/mlength;
//double R_op   = (*mpProperties)[FC]/(*mpProperties)[FT];
double f      = mpFluencyCriteria->mMultisurface_Platicity_Yield[0];



if(f==0.00) 
 {
    mcurrent_plastic_damage = 1.00;
 }
else
 {
    //gc_p = (gc_p/f)*teta_b; 
    //gf_p = (gf_p*R_op/f)*teta_b; 
    array_1d<double, 4> Aux_Strain = Strain - mplastic_strain;    
    double aux_a = inner_prod(Stress, Aux_Strain); 
    mcurrent_disipation = mdisipation + aux_a;  

    kp_punto = ( teta/gf_p + (1.00-teta)/(gc_p) ) * aux_a;
    mcurrent_plastic_damage  = mplastic_damage + kp_punto;                 
    if(mcurrent_plastic_damage>=1.00)
    {
    //KRATOS_WATCH("DAMAGE__DAMAGE__DAMAGE__DAMAGE__DAMAGE__DAMAGE__DAMAGE__DAMAGE__DAMAGE__DAMAGE" ) 
    mcurrent_plastic_damage = 1.00;
    }

 }
//KRATOS_WATCH(kp_punto)
//KRATOS_WATCH(mplastic_damage)
//KRATOS_WATCH(mcurrent_plastic_damage)


/*
//updating cohesion
//mcurrent_cohesion = mcurrent_cohesion;
double hc = 0.00;
double Ct = 0.00;
double Cc = 0.00;

double dCt = 0.00;
double dCc = 0.00;

double cohesion_punto = 0.00;

// WARNING =  Valido solo para morh coulomb
double N  =  2.00*sqrt(R_op);
double cohesion = (*mpProperties)[FC]/(2.00*tan(mmaxfriction_angle/2.00 + PI/4.00));
mpSofteningBehavior->FunctionSofteningHardeningBehavior(mcurrent_plastic_damage, (*mpProperties)[FC] , Cc, dCc);
mpSofteningBehavior->FunctionSofteningHardeningBehavior(mcurrent_plastic_damage, (*mpProperties)[FT] , Ct, dCt); 

Cc = Cc/N;
Ct = Ct * R_op / N;

//KRATOS_WATCH(Cc)
//KRATOS_WATCH(Ct)
if(Cc== 0.00 or Ct == 0.00)
 { hc =  0.00; }
else
{
hc = (teta*dCt/(Ct) + (1.00 - teta)*dCc/(Cc))*cohesion; //mcurrent_cohesion;
}
cohesion_punto      = hc*kp_punto;
mcurrent_cohesion   = mcohesion + cohesion_punto;
if(mcurrent_cohesion<0.00){mcurrent_cohesion = 0.00;}   
//mcurrent_cohesion = mcurrent_cohesion;
//KRATOS_WATCH(mcurrent_cohesion)
*/

//updating friction
mcurrent_friction_angle =  asin(2.00*sqrt(mcurrent_plastic_damage)/(mcurrent_plastic_damage +  1.00)*sin(mmaxfriction_angle));
//KRATOS_WATCH(mcurrent_friction_angle)
//mcurrent_friction_angle = mcurrent_friction_angle;
 
//updating dilatancy
double sin_fric_max      = sin(mmaxfriction_angle);
double sin_dila_max      = sin(mmaxdilatancy_angle);
double sin_fric          = sin(mcurrent_friction_angle);
double sin_phi_cv        = (sin_fric_max-sin_dila_max)/(1.00 - sin_fric_max*sin_dila_max);

mcurrent_dilatancy_angle =  (sin_fric - sin_phi_cv)/(1.00 - sin_phi_cv*sin_fric);
//KRATOS_WATCH(sin_fric - sin_phi_cv);

//mcurrent_dilatancy_angle = mcurrent_dilatancy_angle;

if(mcurrent_dilatancy_angle <0.00) {mcurrent_dilatancy_angle = 0.00;}
else{mcurrent_dilatancy_angle = asin(mcurrent_dilatancy_angle);}
 
mcurrent_dilatancy_angle = std::min(mmaxdilatancy_angle, mcurrent_dilatancy_angle); 

}


void  PlasticDamage2D::Tensile_Fracture_Model(const Vector& Stress, const Vector& plastic_strain )
{
// Anisotropic Rotating crack  Model
const int    iter      = 50;
const double zero      = 1.0E-9;
double& ft = (*mpProperties)[FT];
double& Gf = (*mpProperties)[FRACTURE_ENERGY];

double critical_strain = 2.00 * Gf /(mlength * ft ); // WARNING = Unidades en Kg - cm 
 
array_1d<double,3> H =  ZeroVector(3);

//array_1d<double,3> Rotating;
Vector principal_strain_current = ZeroVector(3);
Vector principal_strain_old     = ZeroVector(3);
Matrix EigenVectors             = ZeroMatrix(3,3);
Matrix StrainTensor             = ZeroMatrix(3,3);

/// plastic strain = current plastic strain
StrainTensor(0,0) = plastic_strain[0];     StrainTensor(0,1) = 0.5*plastic_strain[2];   StrainTensor(0,2) = 0.00;
StrainTensor(1,0) = 0.5*plastic_strain[2]; StrainTensor(1,1) = plastic_strain[1];       StrainTensor(1,2) = 0.00;
StrainTensor(2,0) = 0.00;                  StrainTensor(2,1) = 0.00;                    StrainTensor(2,2) = plastic_strain[3];
SD_MathUtils<double>::EigenVectors(StrainTensor, EigenVectors, principal_strain_current, zero, iter);
sort (principal_strain_current.begin(), principal_strain_current.end()); 
reverse(principal_strain_current.begin(), principal_strain_current.end());


///mplastic_strain = the old strain before do the updtate
StrainTensor(0,0) = mplastic_strain[0];     StrainTensor(0,1) = 0.5*mplastic_strain[2];   StrainTensor(0,2) = 0.00;
StrainTensor(1,0) = 0.5*mplastic_strain[2]; StrainTensor(1,1) = mplastic_strain[1];       StrainTensor(1,2) = 0.00;
StrainTensor(2,0) = 0.00;                   StrainTensor(2,1) = 0.00;                     StrainTensor(2,2) = mplastic_strain[3];
SD_MathUtils<double>::EigenVectors(StrainTensor, EigenVectors, principal_strain_old, zero, iter);
sort (principal_strain_old.begin(), principal_strain_old.end()); 
reverse(principal_strain_old.begin(), principal_strain_old.end());

//KRATOS_WATCH(principal_strain_old)
//KRATOS_WATCH(principal_strain_current)

array_1d<double, 3> Delta_Strain                 = ZeroVector(3);
noalias(Delta_Strain) = principal_strain_current - principal_strain_old; 
for (unsigned int i = 0; i < 3; i++)
{
    if(Delta_Strain[i]<=0.00)  {Delta_Strain[i] = 0.00; }      
    if(principal_strain_current[i]<=critical_strain){ H[i]   = mlength * ft * ft/( 2.00 * Gf  );  } 
}  

          
// Updating mFt in each direction

//KRATOS_WATCH(H)
//KRATOS_WATCH(principal_strain_old)
//KRATOS_WATCH(principal_strain_current)

mcurrent_Ft[0] = mFt[0] -  H[0] * Delta_Strain[0];
mcurrent_Ft[1] = mFt[1] -  H[1] * Delta_Strain[1];
mcurrent_Ft[2] = mFt[2] -  H[2] * Delta_Strain[2];

mpFluencyCriteria_Traction->UpdateVariables(mcurrent_Ft); 

//if(mcurrent_Ft[0] <= 0.00) { mcurrent_Ft[0] = 0.00; }
//if(mcurrent_Ft[1] <= 0.00) { mcurrent_Ft[1] = 0.00; }
//if(mcurrent_Ft[2] <= 0.00) { mcurrent_Ft[2] = 0.00; }
//KRATOS_WATCH(mcurrent_Ft)

}

void PlasticDamage2D::Coupling_Between_Degradetion_In_Compresion_And_Tension(array_1d<double, 3> inelastic_strain)
{

inelastic_strain[0] = 1.00 + sin(mcurrent_dilatancy_angle); 
inelastic_strain[1] = 0.00; 
inelastic_strain[2] = sin(mcurrent_dilatancy_angle) - 1.00;
//noalias(inelastic_strain) = mgamma * inelastic_strain; 
}


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



void PlasticDamage2D::Compute_Derivate(vector<Vector>& Derivate, vector<array_1d<double,4> >& mD)
{
mD.resize(4);
mD[0](0) =   Derivate[0](0);   /// xx
mD[0](1)  =  Derivate[0](1);   /// yy
mD[0](2)  =  Derivate[0](3);   /// xy    
mD[0](3)  =  Derivate[0](2);   /// zz

mD[1](0) =   Derivate[1](0);   
mD[1](1)  =  Derivate[1](1);  
mD[1](2)  =  Derivate[1](3);    
mD[1](3)  =  Derivate[1](2); 

mD[2](0) =   Derivate[2](0);   
mD[2](1)  =  Derivate[2](1);  
mD[2](2)  =  Derivate[2](3);    
mD[2](3)  =  Derivate[2](2); 

mD[3](0) =   Derivate[3](0);   
mD[3](1)  =  Derivate[3](1);  
mD[3](2)  =  Derivate[3](3);    
mD[3](3)  =  Derivate[3](2); 

}
}

