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
*   Date:                $Date: 2008-01-24 16:48:22 $
*   Revision:            $Revision: 1.2 $
*
* ***********************************************************/
// System includes 
#include <iostream>

// External includes 
#include<cmath>

// Project includes 

#include "includes/define.h"
#include "constitutive_laws/hooks_law.h"

#include "includes/constitutive_law.h"

#include "utilities/math_utils.h"
#include "custom_utilities/sd_math_utils.h"
#include "includes/variables.h"
#include "includes/process_info.h"
#include "structural_application.h"
#include "includes/properties.h"

namespace Kratos
{
	//default constructor
    HooksLaw::HooksLaw() 
    : ConstitutiveLaw<Node<3> >()
    {
    }
    //default desstructor
    HooksLaw::~HooksLaw()
    {
    }

    boost::shared_ptr<ConstitutiveLaw<Node<3> > > HooksLaw::Clone() const
    {
        boost::shared_ptr<ConstitutiveLaw<Node<3> > > p_clone(new HooksLaw());
        return p_clone;
    }
    
    //*********************************************************************
    //*********************************************************************
	// Initialization of the coonstitutive law at the begion of the time step
    //*********************************************************************
    //*********************************************************************
    void HooksLaw::InitializeMaterial( const Properties& props,
					const GeometryType& geom,
					const Vector& ShapeFunctionsValues)
    {
		mInsituStress.resize(3,3);
		mCurrentStress.resize(3,3);
        mE = props[YOUNG_MODULUS];
        mNU = props[POISSON_RATIO];
		//set up the material law
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

    //**********************************************************************
    //**********************************************************************
	// not used in the exercises
    //**********************************************************************
    //**********************************************************************
 	void HooksLaw::SetValue( const Variable<Matrix >& rVariable, 
					const Matrix& Value, const ProcessInfo& rCurrentProcessInfo)
	{
    }
    //**********************************************************************
    //**********************************************************************
	// not used in the exercises
    //**********************************************************************
    //**********************************************************************
 	void HooksLaw::SetValue( const Variable<Vector >& rVariable, 
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
    Matrix HooksLaw::GetValue(const Variable<Matrix>& rVariable)
    { 
		return ZeroMatrix(0,0);
    }

    //**********************************************************************
    //**********************************************************************
	// not used in the exercises
    //**********************************************************************
    //**********************************************************************
    Vector HooksLaw::GetValue(const Variable<Vector>& rVariable)
    { 
       if( rVariable == INSITU_STRESS )
        {
			Vector rResult(6);
			rResult(0)= mInsituStress(0,0);
			rResult(1)= mInsituStress(1,1);
			rResult(2)= mInsituStress(2,2);
			rResult(3)= mInsituStress(0,1);
			rResult(4)= mInsituStress(1,2);
			rResult(5)= mInsituStress(2,0);

            return rResult;
        }
		else
			return ZeroVector(0);
    }
    //**********************************************************************
    //**********************************************************************
	// not used in the exercises
    //**********************************************************************
    //**********************************************************************
    double HooksLaw::GetValue(const Variable<double>& rVariable)
    { 
		return 0.0;
    }
    //**********************************************************************
    //**********************************************************************
	void HooksLaw::InitializeSolutionStep( const Properties& props,
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
	//this method is called after each solution step ("time step")
    //**********************************************************************
    //**********************************************************************
    void HooksLaw::FinalizeSolutionStep( const Properties& props,
					const GeometryType& geom, const Vector& ShapeFunctionsValues ,const ProcessInfo& CurrentProcessInfo)
    {
		if( CurrentProcessInfo[CALCULATE_INSITU_STRESS]&& !(CurrentProcessInfo[FIRST_TIME_STEP]))
		{
			noalias(mInsituStress) = mCurrentStress;
			return;
		}
    }
    //**********************************************************************
    //**********************************************************************
	// This is the main method called from outside (inside the element at each quadrature point)
    //**********************************************************************
    //**********************************************************************
    void HooksLaw::CalculateStressAndTangentMatrix(Matrix& StressTensor, 
            const Matrix& StrainTensor, 
            MaterialTensorType& algorithmicTangent)
    {
        KRATOS_TRY

        unsigned int dim = 3;

        if(StressTensor.size1() != dim || StressTensor.size2() != dim)
            StressTensor.resize(dim, dim);
        noalias(StressTensor)= ZeroMatrix(dim, dim);

        algorithmicTangent= mElasticMaterialTensor;

		for(unsigned int i=0; i<3; i++)
			for(unsigned int j=0; j<3; j++)
				for(unsigned int k=0; k<3; k++)
					for(unsigned int l=0; l<3; l++)
						StressTensor(i,j)+=(mElasticMaterialTensor[27*i+9*j+3*k+l])*(StrainTensor(k,l));

		StressTensor+= mInsituStress;

        mCurrentStress= StressTensor;

        KRATOS_CATCH("")
    }
    
    void HooksLaw::CalculateStressAndTangentMatrix( Vector& StressVector,
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
    
} // Namespace Kratos
