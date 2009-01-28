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
*   Last Modified by:    $Author: janosch $
*   Date:                $Date: 2008-01-24 16:49:04 $
*   Revision:            $Revision: 1.9 $
*
* ***********************************************************/
// System includes 
#include <iostream>

// External includes 
#include<cmath>

// Project includes 

#include "includes/define.h"
#include "constitutive_laws/isotropic_elastic_large_strain.h"

#include "includes/constitutive_law.h"

#include "utilities/math_utils.h"
#include "custom_utilities/sd_math_utils.h"
#include "includes/variables.h"
#include "includes/process_info.h"
#include "structural_application.h"
#include "includes/properties.h"

namespace Kratos
{
    /**
     *	TO BE TESTED!!!
     */
    IsotropicElasticLargeStrain::IsotropicElasticLargeStrain() 
    : ConstitutiveLaw<Node<3> >()
    {
    }
    
    /**
     *	TO BE TESTED!!!
     */

    IsotropicElasticLargeStrain::~IsotropicElasticLargeStrain()
    {
    }

    boost::shared_ptr<ConstitutiveLaw<Node<3> > > IsotropicElasticLargeStrain::Clone() const
    {
        boost::shared_ptr<ConstitutiveLaw<Node<3> > > p_clone(new IsotropicElasticLargeStrain());
        return p_clone;
    }
    
        //**********************************************************************
    void IsotropicElasticLargeStrain::InitializeMaterial( const Properties& props,
					const GeometryType& geom,
					const Vector& ShapeFunctionsValues)
    {
        mE = props[YOUNG_MODULUS];
        mNU = props[POISSON_RATIO];
        
        Matrix Unity(3,3);
        noalias(Unity)= ZeroMatrix(3,3);
        for(int i=0; i<3;i++)
            Unity(i,i)=1.0;
        
        SetVariable(ELASTIC_LEFT_CAUCHY_GREEN_OLD, Unity);

        if(mUpdatedLeftCauchyGreenTensor.size1()!=3 || 
                        mUpdatedLeftCauchyGreenTensor.size2()!=3)
            mUpdatedLeftCauchyGreenTensor.resize(3,3);
        
        noalias(mUpdatedLeftCauchyGreenTensor)= Unity;
        
        if(mCurrentStressTensor.size1()!=3 || 
           mCurrentStressTensor.size2()!=3)
            mCurrentStressTensor.resize(3,3);
        
        noalias(mCurrentStressTensor)= ZeroMatrix(3,3);
        
        if(mInsituStressTensor.size1()!=3 || 
           mInsituStressTensor.size2()!=3)
            mInsituStressTensor.resize(3,3);
        
        noalias(mInsituStressTensor)= ZeroMatrix(3,3);
        
    }
    //**********************************************************************
 	void IsotropicElasticLargeStrain::SetValue( const Variable<Matrix >& rVariable, 
					const Matrix& Value, const ProcessInfo& rCurrentProcessInfo)
    {
        if(rVariable== ELASTIC_LEFT_CAUCHY_GREEN_OLD)
        {
            if(mLeftCauchyGreenTensor_Old.size1()!=3 || mLeftCauchyGreenTensor_Old.size2()!=3)
                mLeftCauchyGreenTensor_Old.resize(3,3);
            noalias(mLeftCauchyGreenTensor_Old)= Value;
        }
    }
    //**********************************************************************
 	void IsotropicElasticLargeStrain::SetValue( const Variable<Vector >& rVariable, 
					const Vector& rValue, const ProcessInfo& rCurrentProcessInfo)
    {
		//INSITU_STRESS is stored as a 6*1 Vector
        if( rVariable == INSITU_STRESS )
        {
            if(mInsituStressTensor.size1()!=3 || mInsituStressTensor.size2())
                mInsituStressTensor.resize(3,3);

			mInsituStressTensor(0,0)=rValue(0);
			mInsituStressTensor(1,1)=rValue(1);
			mInsituStressTensor(2,2)=rValue(2);
			mInsituStressTensor(0,1)=rValue(3);mInsituStressTensor(1,0)=rValue(3);
			mInsituStressTensor(1,2)=rValue(4);mInsituStressTensor(2,1)=rValue(4);
			mInsituStressTensor(2,0)=rValue(5);mInsituStressTensor(0,2)=rValue(5);
        }
    }

    //**********************************************************************
	//Internal use only (private)
    void IsotropicElasticLargeStrain::SetVariable( const Variable<Matrix>& rVariable, Matrix& rValue)
    {
        if(rVariable== ELASTIC_LEFT_CAUCHY_GREEN_OLD)
        {
            if(mLeftCauchyGreenTensor_Old.size1()!=3 || mLeftCauchyGreenTensor_Old.size2()!=3)
                mLeftCauchyGreenTensor_Old.resize(3,3);
            noalias(mLeftCauchyGreenTensor_Old)= rValue;
        }
    }
    //**********************************************************************
	//Internal use only (private)
    void IsotropicElasticLargeStrain::SetVariable( const Variable<Vector>& rVariable, Matrix& rValue)
    {
        if( rVariable == INSITU_STRESS )
        {
            if(mInsituStressTensor.size1()!=3 || mInsituStressTensor.size2())
                mInsituStressTensor.resize(3,3);
            noalias(mInsituStressTensor)= rValue;
        }
    }
    //**********************************************************************
    Matrix IsotropicElasticLargeStrain::GetValue(const Variable<Matrix>& rVariable)
    { 
        if( rVariable == ELASTIC_LEFT_CAUCHY_GREEN_OLD )
        {
            return mLeftCauchyGreenTensor_Old;
        }
		else
			return ZeroMatrix(0,0);
    }

    //**********************************************************************
    Vector IsotropicElasticLargeStrain::GetValue(const Variable<Vector>& rVariable)
    { 
        if( rVariable == INSITU_STRESS )
        {
			Vector rResult(6);
			rResult(0)= mInsituStressTensor(0,0);
			rResult(1)= mInsituStressTensor(1,1);
			rResult(2)= mInsituStressTensor(2,2);
			rResult(3)= mInsituStressTensor(0,1);
			rResult(4)= mInsituStressTensor(1,2);
			rResult(5)= mInsituStressTensor(2,0);

            return rResult;
        }
		else
			return ZeroVector(0);
    }
    //**********************************************************************
    //**********************************************************************
    void IsotropicElasticLargeStrain::FinalizeSolutionStep( const Properties& props,
					const GeometryType& geom, const Vector& ShapeFunctionsValues ,const ProcessInfo& CurrentProcessInfo)
    {

		if( CurrentProcessInfo[CALCULATE_INSITU_STRESS]&& !(CurrentProcessInfo[FIRST_TIME_STEP]))
		{
        	Matrix Unity(3,3);
        	noalias(Unity)= ZeroMatrix(3,3);
        	for(int i=0; i<3;i++)
            		Unity(i,i)=1.0;
			SetVariable(ELASTIC_LEFT_CAUCHY_GREEN_OLD, Unity);
			noalias(mInsituStressTensor) -= mCurrentStressTensor;
		}
		else
        	SetVariable(ELASTIC_LEFT_CAUCHY_GREEN_OLD, mUpdatedLeftCauchyGreenTensor);
    }
     //**********************************************************************       
    void IsotropicElasticLargeStrain::CalculateStressAndTangentMatrix(Matrix& StressTensor, 
            const Matrix& LeftCauchyGreenTensor_trial, 
            std::vector<std::vector<Matrix> >& algorithmicTangent)
    {
        KRATOS_TRY

        unsigned int dim = 3;

        if(StressTensor.size1() != dim || StressTensor.size2() != dim)
            StressTensor.resize(dim, dim);
        noalias(StressTensor)= ZeroMatrix(dim, dim);
//                        CalculateTrialStretches(LeftCauchyGreenTensor,DN_DX );
        //////spectral decomposition of Trial state left elastic Cauchy Green Tensor 
        //Calulate principal stretches and principal directions via Gauss-Seidel-elimination
        Vector stretches(dim);//principal stretches
        noalias(stretches)= ZeroVector(dim);
        Vector henky(dim);//principal henky strains
        noalias(henky)= ZeroVector(dim);
               //Calculate Principal stretches and Henky Strains
        SpectralDecomposition(LeftCauchyGreenTensor_trial, stretches, henky);
               //Calculate principal stresses and material law in principal state
        Vector principalStresses(dim);
        noalias(principalStresses)= ZeroVector(dim);
        Matrix aep(dim,dim);
        noalias(aep)= ZeroMatrix(dim,dim);

                // call Materialroutine for small displacements in principal state
        CalculateStressAndTangentialStiffness_PrincipalState
                (principalStresses, aep, henky); 

        Matrix LeftCauchyGreenTensor_Updated(3,3);
                //inverse spectral decomposition
        InverseSpectralDecomposition(aep, principalStresses, 
                                     stretches,henky, LeftCauchyGreenTensor_trial, algorithmicTangent,
                                     StressTensor, LeftCauchyGreenTensor_Updated);

        if(mLeftCauchyGreenTensor_Old.size1()!=3 || mLeftCauchyGreenTensor_Old.size2()!=3 )
            mLeftCauchyGreenTensor_Old.resize(3,3);

        noalias(mUpdatedLeftCauchyGreenTensor)= LeftCauchyGreenTensor_Updated;

		//Update of Current StressTensor
        noalias(mCurrentStressTensor) = StressTensor;
        //Appliance of In_Situ_Stress
        noalias(StressTensor) -= mInsituStressTensor;

        KRATOS_CATCH("")
    }
    
    void IsotropicElasticLargeStrain::CalculateStressAndTangentMatrix(Matrix& StressTensor, 
            const Matrix& LeftCauchyGreenTensor_trial, 
            BaseType::MaterialTensorType& algorithmicTangent)
            {
                KRATOS_TRY

                        unsigned int dim = 3;

                if(StressTensor.size1() != dim || StressTensor.size2() != dim)
                    StressTensor.resize(dim, dim);
                noalias(StressTensor)= ZeroMatrix(dim, dim);
//                        CalculateTrialStretches(LeftCauchyGreenTensor,DN_DX );
        //////spectral decomposition of Trial state left elastic Cauchy Green Tensor 
        //Calulate principal stretches and principal directions via Gauss-Seidel-elimination
                Vector stretches(dim);//principal stretches
                noalias(stretches)= ZeroVector(dim);
                Vector henky(dim);//principal henky strains
                noalias(henky)= ZeroVector(dim);
               //Calculate Principal stretches and Henky Strains
                SpectralDecomposition(LeftCauchyGreenTensor_trial, stretches, henky);
               //Calculate principal stresses and material law in principal state
                Vector principalStresses(dim);
                noalias(principalStresses)= ZeroVector(dim);
                Matrix aep(dim,dim);
                noalias(aep)= ZeroMatrix(dim,dim);

                // call Materialroutine for small displacements in principal state
                CalculateStressAndTangentialStiffness_PrincipalState
                        (principalStresses, aep, henky); 

                Matrix LeftCauchyGreenTensor_Updated(3,3);
                //inverse spectral decomposition
                InverseSpectralDecomposition(aep, principalStresses, 
                                             stretches,henky, LeftCauchyGreenTensor_trial, algorithmicTangent,
                                             StressTensor, LeftCauchyGreenTensor_Updated);

                if(mLeftCauchyGreenTensor_Old.size1()!=3 || mLeftCauchyGreenTensor_Old.size2()!=3 )
                    mLeftCauchyGreenTensor_Old.resize(3,3);

                noalias(mUpdatedLeftCauchyGreenTensor)= LeftCauchyGreenTensor_Updated;

		//Update of Current StressTensor
                noalias(mCurrentStressTensor) = StressTensor;
        //Appliance of In_Situ_Stress
                noalias(StressTensor) -= mInsituStressTensor;

                KRATOS_CATCH("")
            }
    
    
    
    //**********************************************************************
    void IsotropicElasticLargeStrain::SpectralDecomposition
            (const Matrix& LeftCauchyGreenTensor, 
             Vector& stretches, Vector& henky)
            {
                KRATOS_TRY
                       
            //tolerances for calculation of eigenvalues

           //solving charateristic polynom
                double a= 1.0;
                double b= -(LeftCauchyGreenTensor(0,0)+LeftCauchyGreenTensor(1,1)
                            +LeftCauchyGreenTensor(2,2));
                double c= 0.0;
                for(int i=0; i<3;i++)
                    for(int j=0; j<3; j++)
                        c+= 0.5*(LeftCauchyGreenTensor(i,i)*LeftCauchyGreenTensor(j,j)-
                                LeftCauchyGreenTensor(i,j)*LeftCauchyGreenTensor(j,i));
                double d= 0.0;
                for(int i=0; i<3;i++)
                    for(int j=0; j<3; j++)
                        for(int k=0; k<3; k++)
                            d -= 1.0/3.0*(0.5*
                                    LeftCauchyGreenTensor(i,i)*LeftCauchyGreenTensor(j,j)
                                    * LeftCauchyGreenTensor(k,k)
                                    +LeftCauchyGreenTensor(i,j)*LeftCauchyGreenTensor(j,k)
                                    * LeftCauchyGreenTensor(k,i))
                                    -0.5*(LeftCauchyGreenTensor(i,j)*LeftCauchyGreenTensor(j,i)
                                    * LeftCauchyGreenTensor(k,k));

                Matrix n(3,3);
                noalias(n)= ZeroMatrix(3,3);
                
                if(!(SD_MathUtils<double>::CardanoFormula(a, b, c, d, stretches)))
                    SD_MathUtils<double>::EigenVectors( LeftCauchyGreenTensor, n, stretches);

           
                Matrix newStrain(3,3);
                noalias(newStrain)= ZeroMatrix(3,3);

            //check on repeated principal stretches
                stretches= PerturbateLambda(stretches);

            //calculation of logarithmic strains (Henky Strain) 
                for(int i=0; i<3; i++)
                    henky(i)= log(sqrt(stretches(i)));
            
                KRATOS_CATCH("")

            }
            
    //**********************************************************************            
            
            Vector IsotropicElasticLargeStrain::PerturbateLambda(const Vector& stretches)
            {
            //check on repeated roots and perturbate repeated roots
                double max = stretches(0);//maximal prinpical stretch
                if (stretches(1)> max)
                    max= stretches(1);
                if (stretches(2)> max)
                    max= stretches(2);
            
                double tolerance= max*1e-12;

                Vector result(stretches.size());
                noalias(result)= stretches;

                
                for(int i=0; i<3; i++)
                {  
                    for(int j=0; j<3; j++)
                    {
                        if(i != j)
                        {
                            if(fabs(result(i)-result(j))< tolerance )
                            {
                                result(i) -= 1e-8;
                                result(j) += 1e-8;
                            }
                        }
                    }
                }
            
                return result;
            }
            
    //**********************************************************************  
                       
            void 
                IsotropicElasticLargeStrain::CalculateStressAndTangentialStiffness_PrincipalState
                    (Vector& principalStresses, Matrix& aep,const Vector& logStrains)
                    {
                        KRATOS_TRY

                        unsigned int dim= logStrains.size();
            
                        noalias(principalStresses)= ZeroVector(dim);
                        if(aep.size1()!=dim || aep.size2()!=dim)
                            aep.resize(dim,dim);
                        noalias(aep)=ZeroMatrix(dim,dim);

                        double lambda= mNU*mE/((1+mNU)*(1-2*mNU));
                        double mu= mE/(2*(1+mNU));
                        double kappa= lambda+2.0/3.0*mu;
                
                        Matrix kronecker(3,3);
                        noalias(kronecker)=ZeroMatrix(3,3);
            
                        for(unsigned int i=0; i<3;i++)
                        { 
                            kronecker(i,i)=1;
                        }
                
                        for(unsigned int i=0; i<3; i++)
                        {
                            for(unsigned int j=0; j<3; j++)
                            {
                                aep(i,j)= kappa+2.0*mu*(kronecker(i,j)-1.0/3.0);
                            }
                        }
                
                        for(unsigned int i=0; i< dim; i++)
                        {
                            for(unsigned int j=0; j<dim; j++)
                            {
                                principalStresses(i)+= aep(i,j)*logStrains(j);
                            }
                        }
             
                        KRATOS_CATCH("")
                    }
                    
//**********************************************************************  

         void IsotropicElasticLargeStrain::InverseSpectralDecomposition(const Matrix aep, 
                            const Vector& principalStresses, 
                            const Vector& stretches, const Vector& henky, 
                            const Matrix& LeftCauchyGreenTensor, 
                            std::vector<std::vector<Matrix> >& tanC,  
                            Matrix& StressTensor, Matrix& UpdatedLeftCauchyGreenTensor)
            {
        

            KRATOS_TRY

           unsigned  int dim = 3;

            //m[A]=n[A]*n[A] second order Tensor
            std::vector<Matrix> m(dim);
            
            for(unsigned int i=0; i<dim; i++ )
            {
                   if(m[i].size1()!=dim ||m[i].size2()!=dim)
                          m[i].resize(dim,dim);
                   noalias(m[i]) = ZeroMatrix(dim,dim);
            }

            //Calculate m[A] after Simo
            Matrix kronecker(dim,dim);
            noalias(kronecker)=ZeroMatrix(dim,dim);
            for(unsigned int i=0; i<3; i++)
                 kronecker(i,i)=1.0;
            for(unsigned int A=0; A<3; A++)
            {
                unsigned int B = 0; unsigned int C = 0;
                if(A==0)
                  { B=1; C=2;}
                if(A==1)
                   { B=2; C=0;}
                if(A==2)
                   { B=0; C=1;}
                for(unsigned int i=0; i<3; i++)
                {
                    for(unsigned int j=0; j<3; j++)
                    {
			for(unsigned int k=0; k<3; k++)
			{
				m[A](i,j)+=
                                (LeftCauchyGreenTensor(i,k)-stretches(B)*kronecker(i,k))
                                /(stretches(B)-stretches(A))
                                                    *
                                (LeftCauchyGreenTensor(k,j)-stretches(C)*kronecker(k,j))
                                /(stretches(C)-stretches(A));
			}
                     }
                }
            }
            //Calculate Kirchhoff Stresses
            if(StressTensor.size1()!=dim ||StressTensor.size2()!=dim)
                  StressTensor.resize(dim, dim);
            noalias(StressTensor)= ZeroMatrix(dim, dim);
            
            for(unsigned int A=0; A<dim; A++)
            {
                 for(unsigned int k=0; k<dim; k++)
                 {
                      for(unsigned int l=0; l<dim; l++)
                      {
                          StressTensor(k,l) += principalStresses(A)*m[A](k,l);
                      }
                 }
            }
            //Update LeftCauchyGreenTensor
            noalias(UpdatedLeftCauchyGreenTensor)=ZeroMatrix(3,3);
            for(unsigned int A=0; A<dim; A++)
            {
                 for(unsigned int k=0; k<dim; k++)
                 {
                      for(unsigned int l=0; l<dim; l++)
                      {
                            UpdatedLeftCauchyGreenTensor(k,l) += 
                                    exp(2*henky(A))*m[A](k,l);
                      }
                 }
             }
            
            //set up tangential stiffness Matrix c=sum_A(sum_B( aep_AB*m_A*m_B))
            for(unsigned int i=0; i<dim; i++)
            {
                 for(unsigned int j=0; j<dim; j++)
                 {       
                       if(tanC[i][j].size1()!=dim || tanC[i][j].size2() !=dim)
                                tanC[i][j].resize(dim, dim);
                       noalias(tanC[i][j])=ZeroMatrix(dim, dim);
                 }
            }
            //Calculate material part
            for(unsigned int A=0; A<3; A++)
            {
                  for(unsigned int B=0; B<3; B++)
                  {
                        for(unsigned int i=0; i<3;i++)
                        {
                              for(unsigned int j=0; j<3;j++)
                              {
                                    for(unsigned int k=0; k<3;k++)
                                    {  
                                        for(unsigned int l=0; l<3;l++)
                                        { 
                                           tanC[i][j](k,l)+= aep(A,B)*m[A](i,j)*m[B](k,l);
                                        }
                                    }
                               }
                         }
                   }
            }

            //calculate kinematical part
            std::vector<std::vector<Matrix> > Ib(dim,dim);
            std::vector<std::vector<Matrix> > Unity(dim,dim);
            for(unsigned int i=0; i<3; i++)
            {
                 for(unsigned int j=0; j<3; j++)
                 {
                       if(Ib[i][j].size1()!= dim || Ib[i][j].size2()!= dim)
                             Ib[i][j].resize(dim,dim);
                       noalias(Ib[i][j]) = ZeroMatrix(dim,dim);
                       if(Unity[i][j].size1()!= dim || Unity[i][j].size2()!= dim)
                             Unity[i][j].resize(dim,dim);
                       noalias(Unity[i][j]) = ZeroMatrix(dim,dim);
                 }
            }
            for(unsigned int i=0; i<3; i++)
            {
                for(unsigned int j=0; j<3; j++)
                {
                        for(unsigned int k=0; k<3; k++)
                        {
                              for(unsigned int l=0; l<3; l++)
                              {
                                 Ib[i][j](k,l)=0.5*(
                                    LeftCauchyGreenTensor(i,k)*LeftCauchyGreenTensor(j,l)+
                                    LeftCauchyGreenTensor(i,l)*LeftCauchyGreenTensor(j,k));
                                 Unity[i][j](k,l)= 
					kronecker(i,k)*kronecker(j,l);
// 				0.5*(kronecker(i,k)*kronecker(j,l)+kronecker(i,l)*kronecker(j,k));
                               }
                         }                      
                     }
               }

            //Invariants of the LeftCauchyGreenTensor

             double I_1= (LeftCauchyGreenTensor(0,0)+LeftCauchyGreenTensor(1,1)
                                            +LeftCauchyGreenTensor(2,2));

	     double I_3= MathUtils<double>::Det3(LeftCauchyGreenTensor);

             for(unsigned int A=0; A<3; A++)
             {
                 unsigned int B = 0; unsigned int C = 0;
                 if(A==0)
                   { B=1; C=2;}
                 if(A==1)
                   { B=2; C=0;}
                 if(A==2)
                   { B=0; C=1;}
                 double d_A= 
                   (stretches(A)-stretches(B))*(stretches(A)-stretches(C));
                
                 for(unsigned int i=0; i<3; i++)
                 {
                     for(unsigned int j=0; j<3; j++)
                     { 
                         for(unsigned int k=0; k<3; k++)
                         {
                             for(unsigned int l=0; l<3; l++)
                             {
                                 tanC[i][j](k,l) +=
                                     2*principalStresses(A)*(
                                     (1/d_A)*(Ib[i][j](k,l)-
                                     LeftCauchyGreenTensor(i,j)
                                     *LeftCauchyGreenTensor(k,l)-
                                     I_3/stretches(A)*
                                     (Unity[i][j](k,l)-
                                     (kronecker(i,j)-m[A](i,j))
                                     *(kronecker(k,l)-m[A](k,l)))));
                                 tanC[i][j](k,l) +=
                                     2*principalStresses(A)*(        
                                     stretches(A)/d_A*((
                                     LeftCauchyGreenTensor(i,j)*m[A](k,l)
                                     +m[A](i,j)*LeftCauchyGreenTensor(k,l))+
                                     (I_1-4*stretches(A))*m[A](i,j)*m[A](k,l)));
                              }
                          }
                       }
                  }  
            }

//             InitializeMaterialDummy(tanC);

            //durch zwei aufgrund unbekannter Ursache :-)
            for(unsigned int i=0; i<3; i++)
            {
                for(unsigned int j=0; j<3; j++)
                {
                    for(unsigned int k=0; k<3; k++)
                    {
                         for(unsigned int l=0; l<3; l++)
                         {  
                            tanC[i][j](k,l)= tanC[i][j](k,l)/2.0;
                         }
                    }
                }
           }

           KRATOS_CATCH("")
     }
     
     
     void IsotropicElasticLargeStrain::InverseSpectralDecomposition(const Matrix aep, 
             const Vector& principalStresses, 
             const Vector& stretches, const Vector& henky, 
             const Matrix& LeftCauchyGreenTensor, 
             array_1d<double,81>& tanC,  
             Matrix& StressTensor, Matrix& UpdatedLeftCauchyGreenTensor)
             {
        

                 KRATOS_TRY

                         unsigned  int dim = 3;

            //m[A]=n[A]*n[A] second order Tensor
                 std::vector<Matrix> m(dim);
            
                 for(unsigned int i=0; i<dim; i++ )
                 {
                     if(m[i].size1()!=dim ||m[i].size2()!=dim)
                         m[i].resize(dim,dim);
                     noalias(m[i]) = ZeroMatrix(dim,dim);
                 }

            //Calculate m[A] after Simo
                 Matrix kronecker(dim,dim);
                 noalias(kronecker)=ZeroMatrix(dim,dim);
                 for(unsigned int i=0; i<3; i++)
                     kronecker(i,i)=1.0;
                 for(unsigned int A=0; A<3; A++)
                 {
                     unsigned int B = 0; unsigned int C = 0;
                     if(A==0)
                     { B=1; C=2;}
                     if(A==1)
                     { B=2; C=0;}
                     if(A==2)
                     { B=0; C=1;}
                     for(unsigned int i=0; i<3; i++)
                     {
                         for(unsigned int j=0; j<3; j++)
                         {
                             for(unsigned int k=0; k<3; k++)
                             {
                                 m[A](i,j)+=
                                         (LeftCauchyGreenTensor(i,k)-stretches(B)*kronecker(i,k))
                                /(stretches(B)-stretches(A))
                                         *
                                         (LeftCauchyGreenTensor(k,j)-stretches(C)*kronecker(k,j))
                                /(stretches(C)-stretches(A));
                             }
                         }
                     }
                 }
            //Calculate Kirchhoff Stresses
                 if(StressTensor.size1()!=dim ||StressTensor.size2()!=dim)
                     StressTensor.resize(dim, dim);
                 noalias(StressTensor)= ZeroMatrix(dim, dim);
            
                 for(unsigned int A=0; A<dim; A++)
                 {
                     for(unsigned int k=0; k<dim; k++)
                     {
                         for(unsigned int l=0; l<dim; l++)
                         {
                             StressTensor(k,l) += principalStresses(A)*m[A](k,l);
                         }
                     }
                 }
            //Update LeftCauchyGreenTensor
                 noalias(UpdatedLeftCauchyGreenTensor)=ZeroMatrix(3,3);
                 for(unsigned int A=0; A<dim; A++)
                 {
                     for(unsigned int k=0; k<dim; k++)
                     {
                         for(unsigned int l=0; l<dim; l++)
                         {
                             UpdatedLeftCauchyGreenTensor(k,l) += 
                                     exp(2*henky(A))*m[A](k,l);
                         }
                     }
                 }
            
            //set up tangential stiffness Matrix c=sum_A(sum_B( aep_AB*m_A*m_B))
                 if( tanC.size() != 81 )
                     tanC.resize(81, false);
                 tanC = ZeroVector(81);
                 
            //Calculate material part
                 for(unsigned int A=0; A<3; A++)
                 {
                     for(unsigned int B=0; B<3; B++)
                     {
                         for(unsigned int i=0; i<3;i++)
                         {
                             for(unsigned int j=0; j<3;j++)
                             {
                                 for(unsigned int k=0; k<3;k++)
                                 {  
                                     for(unsigned int l=0; l<3;l++)
                                     { 
                                         tanC[27*i+9*j+3*k+l] += aep(A,B)*m[A](i,j)*m[B](k,l);
                                     }
                                 }
                             }
                         }
                     }
                 }

            //calculate kinematical part
                 std::vector<std::vector<Matrix> > Ib(dim,dim);
                 std::vector<std::vector<Matrix> > Unity(dim,dim);
                 for(unsigned int i=0; i<3; i++)
                 {
                     for(unsigned int j=0; j<3; j++)
                     {
                         if(Ib[i][j].size1()!= dim || Ib[i][j].size2()!= dim)
                             Ib[i][j].resize(dim,dim);
                         noalias(Ib[i][j]) = ZeroMatrix(dim,dim);
                         if(Unity[i][j].size1()!= dim || Unity[i][j].size2()!= dim)
                             Unity[i][j].resize(dim,dim);
                         noalias(Unity[i][j]) = ZeroMatrix(dim,dim);
                     }
                 }
                 for(unsigned int i=0; i<3; i++)
                 {
                     for(unsigned int j=0; j<3; j++)
                     {
                         for(unsigned int k=0; k<3; k++)
                         {
                             for(unsigned int l=0; l<3; l++)
                             {
                                 Ib[i][j](k,l)=0.5*(
                                          LeftCauchyGreenTensor(i,k)*LeftCauchyGreenTensor(j,l)+
                                                  LeftCauchyGreenTensor(i,l)*LeftCauchyGreenTensor(j,k));
                                 Unity[i][j](k,l)= 
                                         kronecker(i,k)*kronecker(j,l);
// 				0.5*(kronecker(i,k)*kronecker(j,l)+kronecker(i,l)*kronecker(j,k));
                             }
                         }                      
                     }
                 }

            //Invariants of the LeftCauchyGreenTensor

                 double I_1= (LeftCauchyGreenTensor(0,0)+LeftCauchyGreenTensor(1,1)
                             +LeftCauchyGreenTensor(2,2));

                 double I_3= MathUtils<double>::Det3(LeftCauchyGreenTensor);

                 for(unsigned int A=0; A<3; A++)
                 {
                     unsigned int B = 0; unsigned int C = 0;
                     if(A==0)
                     { B=1; C=2;}
                     if(A==1)
                     { B=2; C=0;}
                     if(A==2)
                     { B=0; C=1;}
                     double d_A= 
                                 (stretches(A)-stretches(B))*(stretches(A)-stretches(C));
                
                     for(unsigned int i=0; i<3; i++)
                     {
                         for(unsigned int j=0; j<3; j++)
                         { 
                             for(unsigned int k=0; k<3; k++)
                             {
                                 for(unsigned int l=0; l<3; l++)
                                 {
                                     tanC[27*i+9*j+3*k+l] +=
                                             2*principalStresses(A)*(
                                                     (1/d_A)*(Ib[i][j](k,l)-
                                                     LeftCauchyGreenTensor(i,j)
                                                     *LeftCauchyGreenTensor(k,l)-
                                                     I_3/stretches(A)*
                                                     (Unity[i][j](k,l)-
                                                             (kronecker(i,j)-m[A](i,j))
                                                             *(kronecker(k,l)-m[A](k,l)))));
                                     tanC[27*i+9*j+3*k+l] +=
                                             2*principalStresses(A)*(        
                                                     stretches(A)/d_A*((
                                                     LeftCauchyGreenTensor(i,j)*m[A](k,l)
                                                     +m[A](i,j)*LeftCauchyGreenTensor(k,l))+
                                                     (I_1-4*stretches(A))*m[A](i,j)*m[A](k,l)));
                                 }
                             }
                         }
                     }  
                 }

//             InitializeMaterialDummy(tanC);

            //durch zwei aufgrund unbekannter Ursache :-)
                 tanC = 0.5*tanC;

                 KRATOS_CATCH("")
             }
     
                            
     void IsotropicElasticLargeStrain::InitializeMaterialDummy
          (std::vector<std::vector<Matrix> >& C)
     {
     KRATOS_TRY
                    
     unsigned int dimension = 3;
             
     Matrix kronecker(3,3);
     noalias(kronecker)=ZeroMatrix(3,3);
            
     for(unsigned int i=0; i<3;i++)
     { 
            kronecker(i,i)=1;
     }
            
     double E = mE;
     double NU = mNU;
     double lambda= NU*E/((1+NU)*(1-2*NU));
     double mu= E/(2*(1+NU));
            
     for(unsigned int i=0; i<dimension;i++)
     {
              for(unsigned int j=0; j<dimension;j++)
              {
              C[i][j].resize(dimension,dimension);	
              noalias(C[i][j])= ZeroMatrix(dimension,dimension);
              for(unsigned int k=0; k<dimension; k++)
              {
              for(unsigned int l=0; l<dimension; l++)
                   {
                          C[i][j](k,l)=lambda*kronecker(i,j)
                                 *kronecker(k,l)
                                 +mu*(kronecker(i,k)*kronecker(j,l)
                                 +kronecker(i,l)*kronecker(j,k));
                   }
              }
         }
    }

    KRATOS_CATCH("")
    } 
        
} // Namespace Kratos
