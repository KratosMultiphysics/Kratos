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
*   Last Modified by:    $Author: hurga $
*   Date:                $Date: 2009-03-05 12:01:22 $
*   Revision:            $Revision: 1.5 $
*
* ***********************************************************/
// System includes 
#include <iostream>

// External includes 
#include<cmath>

// Project includes 

#include "includes/define.h"
#include "constitutive_laws/drucker_prager_law.h"

#include "includes/constitutive_law.h"

#include "utilities/math_utils.h"
#include "custom_utilities/sd_math_utils.h"
#include "includes/variables.h"
#include "includes/process_info.h"
#include "structural_application.h"
#include "includes/properties.h"
#include "linear_solvers/skyline_lu_factorization_solver.h"

namespace Kratos
{
	//default constructor
    DruckerPragerLaw::DruckerPragerLaw() 
    : ConstitutiveLaw()
    {
        mOldPlasticStrains.resize(3,3,false);
        mCurrentPlasticStrains.resize(3,3,false);
        mInsituStress.resize(3,3,false);
        mCurrentStress.resize(3,3,false);
		SD_MathUtils<double>::DeviatoricUnity(mI_dev);
		mTOL1= 1e-1;
		mTOL2= 1e-8;
		mMaxIter= 10;
    }
    //default desstructor
    DruckerPragerLaw::~DruckerPragerLaw()
    {
    }

    boost::shared_ptr<ConstitutiveLaw> DruckerPragerLaw::Clone() const
    {
        boost::shared_ptr<ConstitutiveLaw> p_clone(new DruckerPragerLaw());
        return p_clone;
    }
    
    void DruckerPragerLaw::CalculateMaterialResponse( const Vector& StrainVector,
                                      const Matrix& DeformationGradient,
                                      Vector& StressVector,
                                      Matrix& AlgorithmicTangent,
                                      const ProcessInfo& CurrentProcessInfo,
                                      const Properties& props, 
                                      const GeometryType& geom,
                                      const Vector& ShapeFunctionsValues,
                                      bool CalculateStresses,
                                      int CalculateTangent,
                                      bool SaveInternalVariables )
    {
        CalculateStressAndTangentMatrix(StressVector, StrainVector, AlgorithmicTangent );
    }
    
    //*********************************************************************
    //*********************************************************************
	// Initialization of the coonstitutive law at the begion of the time step
    //*********************************************************************
    //*********************************************************************
    void DruckerPragerLaw::InitializeMaterial( const Properties& props,
					const GeometryType& geom,
					const Vector& ShapeFunctionsValues)
    {
		Matrix oldPlasticStrains(3,3);
		noalias(oldPlasticStrains)=ZeroMatrix(3,3);
		SetOldPlasticStrains( oldPlasticStrains);

        mE = props[YOUNG_MODULUS];
        mNU = props[POISSON_RATIO];
		mFriction= props[DP_ALPHA1];
		mCohesion= props[DP_K];

		Matrix kronecker(3,3);
     	noalias(kronecker)=ZeroMatrix(3,3);           
     	for(unsigned int i=0; i<3;i++)
     	{ 
            kronecker(i,i)=1;
     	}
		//calculate Lame's parameters
     	double lambda= mNU*mE/((1+mNU)*(1-2*mNU));
     	double mu= mE/(2*(1+mNU));
        
        for(unsigned int i=0; i<3;i++)
            for(unsigned int j=0; j<3;j++)
                for(unsigned int k=0; k<3; k++)
                    for(unsigned int l=0; l<3; l++)
                        mElasticMaterialTensor[27*i+9*j+3*k+l]=lambda*kronecker(i,j)
                                *kronecker(k,l)+mu*(kronecker(i,k)*kronecker(j,l)
                                        +kronecker(i,l)*kronecker(j,k));
        
        noalias(mInsituStress)= ZeroMatrix(3,3);
    }

	void DruckerPragerLaw::InitializeSolutionStep( const Properties& props,
		const GeometryType& geom, //this is just to give the array of nodes
		const Vector& ShapeFunctionsValues ,
		const ProcessInfo& CurrentProcessInfo)
	{
        mE = props[YOUNG_MODULUS];
        mNU = props[POISSON_RATIO];
		Matrix kronecker(3,3);
     	noalias(kronecker)=ZeroMatrix(3,3);           
     	for(unsigned int i=0; i<3;i++)
     	{ 
            kronecker(i,i)=1;
     	}
		//calculate Lame's parameters
     	double lambda= mNU*mE/((1+mNU)*(1-2*mNU));
     	double mu= mE/(2*(1+mNU));
        
        for(unsigned int i=0; i<3;i++)
            for(unsigned int j=0; j<3;j++)
                for(unsigned int k=0; k<3; k++)
                    for(unsigned int l=0; l<3; l++)
                        mElasticMaterialTensor[27*i+9*j+3*k+l]=lambda*kronecker(i,j)
                                *kronecker(k,l)+mu*(kronecker(i,k)*kronecker(j,l)
                                        +kronecker(i,l)*kronecker(j,k));
    }
    //**********************************************************************
    //**********************************************************************
    void DruckerPragerLaw::FinalizeSolutionStep( const Properties& props,
					const GeometryType& geom, const Vector& ShapeFunctionsValues ,const ProcessInfo& CurrentProcessInfo)
    {
		if( CurrentProcessInfo[CALCULATE_INSITU_STRESS]&& !(CurrentProcessInfo[FIRST_TIME_STEP]))
		{
			noalias(mInsituStress) = mCurrentStress;
			return;
		}

		SetOldPlasticStrains( mCurrentPlasticStrains);
    }
    //**********************************************************************
    //**********************************************************************
    double& DruckerPragerLaw::GetValue(const Variable<double>& rVariable, double& rValue)
	{
			if(!(rVariable== DP_EPSILON))
				return rValue;
            rValue = SD_MathUtils<double>::normTensor(mCurrentPlasticStrains);
			return rValue;
	}

    //**********************************************************************
    //**********************************************************************
	//Internal use only (private)
    void DruckerPragerLaw::SetOldPlasticStrains( Matrix& rValue)
    {
		if(mOldPlasticStrains.size1() != 3 || mOldPlasticStrains.size2() != 3 )
			mOldPlasticStrains.resize(3,3);
		noalias(mOldPlasticStrains)= rValue;
    }
    //**********************************************************************
    //**********************************************************************
	// not used in the exercises
    //**********************************************************************
    //**********************************************************************
 	void DruckerPragerLaw::SetValue( const Variable<Matrix >& rVariable, 
					const Matrix& Value, const ProcessInfo& rCurrentProcessInfo)
	{
    }
    //**********************************************************************
    //**********************************************************************
	// not used in the exercises
    //**********************************************************************
    //**********************************************************************
 	void DruckerPragerLaw::SetValue( const Variable<Vector >& rVariable, 
					const Vector& rValue, const ProcessInfo& rCurrentProcessInfo)
    {
       if( rVariable == INSITU_STRESS )
       {
			mInsituStress(0,0)= rValue(0);mInsituStress(0,1)= rValue(3);mInsituStress(0,2)= rValue(5);
			mInsituStress(1,0)= rValue(3);mInsituStress(1,1)= rValue(1);mInsituStress(1,2)= rValue(4);
			mInsituStress(2,0)= rValue(5);mInsituStress(2,1)= rValue(4);mInsituStress(2,2)= rValue(2);
	   }
    }
    //**********************************************************************
    //**********************************************************************
	// not used in the exercises
    //**********************************************************************
    //**********************************************************************
    Matrix& DruckerPragerLaw::GetValue(const Variable<Matrix>& rVariable, Matrix& rValue)
    { 
		return( rValue );
    }

    //**********************************************************************
    //**********************************************************************
	// not used in the exercises
    //**********************************************************************
    //**********************************************************************
    Vector& DruckerPragerLaw::GetValue(const Variable<Vector>& rVariable, Vector& rValue)
    { 
        if( rVariable == INSITU_STRESS )
        {
			rValue.resize( 6, false );
			rValue(0)= mInsituStress(0,0);
			rValue(1)= mInsituStress(1,1);
			rValue(2)= mInsituStress(2,2);
			rValue(3)= mInsituStress(0,1);
			rValue(4)= mInsituStress(1,2);
			rValue(5)= mInsituStress(2,0);
        }
        return( rValue );
    }

    //**********************************************************************
    //**********************************************************************
	// This is the main method called from outside (inside the element at each quadrature point)
    //**********************************************************************
    //**********************************************************************
    void DruckerPragerLaw::CalculateStressAndTangentMatrix(Matrix& StressTensor, 
            const Matrix& StrainTensor, 
            MaterialTensorType& algorithmicTangent)
    {
        KRATOS_TRY

		//Transform to trial PlasticStrainVector
		Matrix PlasticStrainTensor(3,3);
		noalias(PlasticStrainTensor)= mOldPlasticStrains;

		//Calculate trial state
		noalias(StressTensor)= CalculateStressTensor(StressTensor, StrainTensor, PlasticStrainTensor);
		//Calculate deviatoric stresses
		Matrix deviatoric_Stress(6,6);
		noalias(deviatoric_Stress)=ZeroMatrix(6,6);
		double firstInv= StressTensor(0,0)+StressTensor(1,1)+StressTensor(2,2);
		for(unsigned int i=0; i<3; i++)
			for(unsigned int j=0; j<3; j++)
			{
				if(i==j) deviatoric_Stress(i,j)= StressTensor(i,j)-1.0/3.0*firstInv;
				else 	deviatoric_Stress(i,j)= StressTensor(i,j);
			}
		//check yield condition
		double normDev= SD_MathUtils<double>::normTensor(deviatoric_Stress);

		double yield= normDev+mFriction*firstInv-mCohesion;

		Matrix tanC(6,6);

		if(yield <= mTOL1 )//this is a elastic incerement
		{
			StressTensor+= mInsituStress;

			mCurrentStress= StressTensor;

			mCurrentPlasticStrains= PlasticStrainTensor;

			algorithmicTangent= mElasticMaterialTensor;

			return;
		}

		//set up Unity Tensor
		Matrix Unity(3,3);
		noalias(Unity) = ZeroMatrix(3,3);
		for(unsigned int i=0; i<3; i++)
			Unity(i,i) = 1.0;

		//set up df_dsigma
		Matrix df_dsigma(3,3);
		//set up d2f_dsigma2
		std::vector<std::vector<Matrix> > d2f_dsigma2;
		//set up the residuum
		Matrix residuum(3,3);
		//set up consistency parameters
		double delta_gamma= 0.0;
		//set up Xi
		std::vector<std::vector<Matrix> > Xi;
		//set up the increments
		Matrix delta_sigma(3,3);
		Matrix delta_epsilon_plastic(3,3);
		double delta2_gamma;
		//InverseC
		Matrix InvC(6,6);
		noalias(InvC)= InverseC(InvC);
		//Help Values
		Matrix Dummy1(6,6);
		Matrix InverseXi(6,6);
		Matrix Dummy2(6,6);
		Matrix Dummy3(3,3);
		double normXi;
		Vector delta_stress_vector(6);
		Vector delta_strain_vector(6);

		// convergence
		bool convergence= false;

		for(unsigned int iter=0; iter< mMaxIter; iter++)
		{
			//calculate df_dsigma
			noalias(df_dsigma) = ZeroMatrix(3,3);
			for(unsigned int i=0; i<3; i++)
				for(unsigned int j=0; j<3; j++)
					df_dsigma(i,j)= deviatoric_Stress(i,j)/normDev+mFriction*Unity(i,j);

			//calculate d2f_dsigma2
			d2f_dsigma2.resize(3);
			for(unsigned int i=0; i<3; i++)
			{
				d2f_dsigma2[i].resize(3);
           	 	for(unsigned int j=0; j<3;j++)
            	{
              		d2f_dsigma2[i][j].resize(3,3);	

              		noalias(d2f_dsigma2[i][j])= ZeroMatrix(3,3);

              		for(unsigned int k=0; k<3; k++)
              		{
              			for(unsigned int l=0; l<3; l++)
                   		{
							d2f_dsigma2[i][j](k,l)=1/normDev*mI_dev[27*i+9*j+3*k+l]-1/(normDev*normDev*normDev)*deviatoric_Stress(i,j)*deviatoric_Stress(k,l);

						}
					}
				}
			}

			//calculate the residuum
			noalias(residuum)= ZeroMatrix(3,3);
			for(unsigned int i=0; i<3; i++)
				for(unsigned int j=0; j<3; j++)
					residuum(i,j)= -PlasticStrainTensor(i,j)+mOldPlasticStrains(i,j)+delta_gamma*df_dsigma(i,j);

			//calculate Xi
			SD_MathUtils<double>::TensorToMatrix(d2f_dsigma2, Dummy1);
			for(unsigned int i=0; i<6; i++)
				for(unsigned int j=0; j<6; j++)
					InverseXi(i,j)= InvC(i,j)+delta_gamma*Dummy1(i,j);

			SD_MathUtils<double>::InvertMatrix( InverseXi, Dummy2);

			SD_MathUtils<double>::MatrixToTensor(Dummy2, Xi);

			//calculate the increments
			//calculate increment of delta_gamma
			delta2_gamma= 0.0; 
			noalias(Dummy3)= ZeroMatrix(3,3);
			for(unsigned int i=0; i<3; i++)
				for(unsigned int j=0; j<3; j++)
					for(unsigned int k=0; k<3; k++)
						for(unsigned int l=0; l<3; l++)
							Dummy3(i,j)+= Xi[i][j](k,l)*df_dsigma(k,l);

			for(unsigned int k=0; k<3; k++)
				for(unsigned int l=0; l<3; l++)
					delta2_gamma+= residuum(k,l)*Dummy3(k,l);

			delta2_gamma= yield-delta2_gamma;

			normXi=0.0;
			for(unsigned int k=0; k<3; k++)
				for(unsigned int l=0; l<3; l++)
					normXi+= df_dsigma(k,l)*Dummy3(k,l);

			delta2_gamma=delta2_gamma/normXi;
			//calculate increment of delta_sigma
			noalias(delta_sigma)= ZeroMatrix(3,3);

			for(unsigned int i=0; i<3; i++)
				for(unsigned int j=0; j<3; j++)
					for(unsigned int k=0; k<3; k++)
						for(unsigned int l=0; l<3; l++)
							delta_sigma(i,j)+= Xi[i][j](k,l)*(-residuum(k,l)-delta2_gamma*df_dsigma(k,l));

			//calculate increment of delta_epsilon_plastic
			noalias(delta_epsilon_plastic)= ZeroMatrix(3,3);
 			SD_MathUtils<double>::TensorToVector(delta_sigma, delta_stress_vector);
			noalias(delta_strain_vector)= ZeroVector(6);
			for(unsigned int i=0; i<6; i++)
				for(unsigned int k=0; k<6; k++)
					delta_strain_vector(i)+= -InvC(i,k)*delta_stress_vector(k);

 			SD_MathUtils<double>::VectorToTensor(delta_strain_vector, delta_epsilon_plastic);
			//update of variables
			delta_gamma+= delta2_gamma;
			PlasticStrainTensor+= delta_epsilon_plastic;

			noalias(StressTensor)= CalculateStressTensor(StressTensor, StrainTensor, PlasticStrainTensor);

			deviatoric_Stress(3,3);
			noalias(deviatoric_Stress)=ZeroMatrix(3,3);
			firstInv= StressTensor(0,0)+StressTensor(1,1)+StressTensor(2,2);
			for(unsigned int i=0; i<3; i++)
				for(unsigned int j=0; j<3; j++)
				{
					if(i==j) deviatoric_Stress(i,j)= StressTensor(i,j)-1.0/3.0*firstInv;
					else 	deviatoric_Stress(i,j)= StressTensor(i,j);
				}
			//check yield condition
			normDev= SD_MathUtils<double>::normTensor(deviatoric_Stress);

			yield= normDev+mFriction*firstInv-mCohesion;

			if(yield <= mTOL1)//this is an elastic incerement
			{
// 				std::cout<<"return map is converged: "<<iter+1<<" steps needed"<<std::endl;
				convergence= true;
				break;
			}
		}
		if(!convergence)
		{
			std::cout<<"return map is not converged"<<std::endl;
			return;
		}
		//update of plastic strain
		mCurrentPlasticStrains=  PlasticStrainTensor;

		//This is for the appliance of an InsituStress
		StressTensor+= mInsituStress;
		mCurrentStress= StressTensor;
		//Calculation of the algorithmic tangent(same procedure as above)
		//calculate df_dsigma
		noalias(df_dsigma) = ZeroMatrix(3,3);
		for(unsigned int i=0; i<3; i++)
			for(unsigned int j=0; j<3; j++)
				df_dsigma(i,j)= deviatoric_Stress(i,j)/normDev+mFriction*Unity(i,j);

		//calculate d2f_dsigma2
		d2f_dsigma2.resize(3);

		for(unsigned int i=0; i<3; i++)
		{
			d2f_dsigma2[i].resize(3);

           	for(unsigned int j=0; j<3;j++)
           	{
              	d2f_dsigma2[i][j].resize(3,3);	

              	noalias(d2f_dsigma2[i][j])= ZeroMatrix(3,3);

              	for(unsigned int k=0; k<3; k++)
              	{
              		for(unsigned int l=0; l<3; l++)
                   	{
						d2f_dsigma2[i][j](k,l)= 	1/normDev*mI_dev[27*i+9*j+3*k+l]-1/(normDev*normDev*normDev)*deviatoric_Stress(i,j)*deviatoric_Stress(k,l);
					}
				}
			}
		}
		//calculate Xi
		SD_MathUtils<double>::TensorToMatrix(d2f_dsigma2, Dummy1);

		for(unsigned int i=0; i<6; i++)
			for(unsigned int j=0; j<6; j++)
				InverseXi(i,j)= InvC(i,j)+delta_gamma*Dummy1(i,j);

		SD_MathUtils<double>::InvertMatrix( InverseXi, Dummy2);

		SD_MathUtils<double>::MatrixToTensor(Dummy2, Xi);

		//calculate the norm of Xi
		noalias(Dummy3)= ZeroMatrix(3,3);
		for(unsigned int i=0; i<3; i++)
			for(unsigned int j=0; j<3; j++)
				for(unsigned int k=0; k<3; k++)
					for(unsigned int l=0; l<3; l++)
						Dummy3(i,j)+= Xi[i][j](k,l)*df_dsigma(k,l);

		normXi=0.0;
		for(unsigned int k=0; k<3; k++)
			for(unsigned int l=0; l<3; l++)
				normXi+= df_dsigma(k,l)*Dummy3(k,l);

		normXi= sqrt(normXi);

		Matrix N(3,3);
		noalias(N)=ZeroMatrix(3,3);
		for(unsigned int i=0; i<3; i++)
			for(unsigned int j=0; j<3; j++)
				for(unsigned int k=0; k<3; k++)
					for(unsigned int l=0; l<3; l++)
						N(i,j)+= Xi[i][j](k,l)*df_dsigma(k,l)/normXi;

		for(unsigned int i=0; i<3; i++)
			for(unsigned int j=0; j<3; j++)
				for(unsigned int k=0; k<3; k++)
					for(unsigned int l=0; l<3; l++)
						algorithmicTangent[27*i+9*j+3*k+l]= Xi[i][j](k,l)-N(i,j)*N(k,l);

		return;
				
        KRATOS_CATCH("")
    }
    
    /**
     */
    void DruckerPragerLaw::CalculateStressAndTangentMatrix( Vector& StressVector,
            const Vector& StrainVector,
            Matrix& algorithmicTangent)
    {
        unsigned int dim = 3;
        if( algorithmicTangent.size1() != 2*dim || algorithmicTangent.size2() != 2*dim )
            algorithmicTangent.resize(2*dim, 2*dim, false);
        noalias(algorithmicTangent)=ZeroMatrix(2*dim,2*dim);
        //generating strain tensor from vector
        Matrix StrainTensor = SD_MathUtils<double>::StrainVectorToTensor( StrainVector );
        //calling tensorial formulation with member variables
        CalculateStressAndTangentMatrix( mCurrentStress, StrainTensor, mElasticMaterialTensor );
        //copying entries from material tensor to output matrix
        SD_MathUtils<double>::TensorToMatrix( mElasticMaterialTensor, algorithmicTangent );
        //copying entries from stress tensor to output vector
        SD_MathUtils<double>::TensorToVector( mCurrentStress, StressVector );
        return;
    }
            

	Matrix DruckerPragerLaw::CalculateStressTensor(Matrix& StressTensor,const Matrix& StrainTensor, Matrix& PlasticStrainTensor)
	{	
        if(StressTensor.size1() != 3 || StressTensor.size2() != 3)
            StressTensor.resize(3, 3);
        noalias(StressTensor)= ZeroMatrix(3, 3);

		Matrix ElasticStrainTensor(3,3);
		noalias(ElasticStrainTensor)= StrainTensor-PlasticStrainTensor;

		for(unsigned int i=0; i<3; i++)
			for(unsigned int j=0; j<3; j++)
				for(unsigned int k=0; k<3; k++)
					for(unsigned int l=0; l<3; l++)
						StressTensor(i,j)+=mElasticMaterialTensor[27*i+9*j+3*k+l]*ElasticStrainTensor(k,l);

		return StressTensor;
	}

	Matrix DruckerPragerLaw::CalculateElasticTangent(Matrix& tanC)
	{
     	double lambda= mNU*mE/((1+mNU)*(1-2*mNU));
     	double mu= mE/(2*(1+mNU));

		if(tanC.size1() != 6 || tanC.size2() != 6 )
			tanC.resize(6,6);
		
		noalias(tanC)= ZeroMatrix(6,6);

		tanC(0,0)= 2*mu+ lambda; tanC(0,1)= lambda;  tanC(0,2)= lambda;
		tanC(1,0)= lambda; tanC(1,1)= 2*mu+ lambda;  tanC(1,2)= lambda;
		tanC(2,0)= lambda; tanC(2,1)= lambda;  tanC(2,2)= 2*mu+ lambda;
		tanC(3,3)= mu; tanC(4,4)= mu; tanC(5,5)= mu;

		return tanC;
	}

	Matrix DruckerPragerLaw::InverseC(Matrix& InvC)
	{
		if(InvC.size1()!=6 || InvC.size2()!=6)
			InvC.resize(6,6);

		noalias(InvC)= ZeroMatrix(6,6);

     	double lambda= mNU*mE/((1+mNU)*(1-2*mNU));
     	double mu= mE/(2*(1+mNU));

		double a= (4*mu*mu+4*mu*lambda)/(8*mu*mu*mu+12*mu*mu*lambda);
		double b= -(2*mu*lambda)/(8*mu*mu*mu+12*mu*mu*lambda);
	
		InvC(0,0)= a; InvC(0,1)= b;  InvC(0,2)= b;
		InvC(1,0)= b; InvC(1,1)= a;  InvC(1,2)= b;
		InvC(2,0)= b; InvC(2,1)= b;  InvC(2,2)= a;
		InvC(3,3)= 1/(2*mu); InvC(4,4)= 1/(2*mu);  InvC(5,5)= 1/(2*mu);

		return InvC;
	}
	

      int DruckerPragerLaw::Check(const Properties& props, const GeometryType& geom, const ProcessInfo& CurrentProcessInfo)
        {  
            if(DENSITY.Key() == 0 || props[DENSITY]<0.00)
                KRATOS_ERROR(std::invalid_argument,"DENSITY has Key zero or invalid value ","");
	   
	    if(YOUNG_MODULUS.Key() == 0 || props[YOUNG_MODULUS]<= 0.00)
                KRATOS_ERROR(std::invalid_argument,"YOUNG_MODULUS has Key zero or invalid value ","");

	    const double& nu = props[POISSON_RATIO];
	    const bool check = bool( (nu >0.499 && nu<0.501 ) || (nu < -0.999 && nu > -1.01 ) );
	    if(POISSON_RATIO.Key() == 0 || check==true) 
                KRATOS_ERROR(std::invalid_argument,"POISSON_RATIO has Key zero invalid value ","");
    
	    if(DP_ALPHA1.Key() == 0 || props[DP_ALPHA1]<= 0.00)
                KRATOS_ERROR(std::invalid_argument,"DP_ALPHA1 has Key zero or invalid value ","");
	    
	    if(DP_K.Key() == 0 || props[DP_K]<= 0.00)
                KRATOS_ERROR(std::invalid_argument,"DP_K has Key zero or invalid value ","");
	    
	    return 0;
         }

} // Namespace Kratos
