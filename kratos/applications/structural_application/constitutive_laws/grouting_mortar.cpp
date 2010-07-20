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
*   Last Modified by:    $Author: nagel $
*   Date:                $Date: 2009-03-20 08:54:08 $
*   Revision:            $Revision: 1.5 $
*
* ***********************************************************/
// System includes 
#include <iostream>

// External includes 
#include<cmath>

// Project includes 

#include "includes/define.h"
#include "constitutive_laws/grouting_mortar.h"

#include "includes/constitutive_law.h"

#include "utilities/math_utils.h"
#include "../structural_application/custom_utilities/sd_math_utils.h"
#include "includes/variables.h"
#include "includes/process_info.h"
#include "ekate_auxiliary_application.h"
#include "includes/properties.h"

namespace Kratos
{
    //*********************************************************************
    //*********************************************************************
    //******* default constructor
    //*********************************************************************
    //*********************************************************************  
    GroutingMortar::GroutingMortar() 
    : ConstitutiveLaw<Node<3> >()
    {
    }
    //*********************************************************************
    //*********************************************************************
    //******* default desstructor
    //*********************************************************************
    //********************************************************************* 
    GroutingMortar::~GroutingMortar()
    {
    }

    boost::shared_ptr<ConstitutiveLaw<Node<3> > > GroutingMortar::Clone() const
    {
        boost::shared_ptr<ConstitutiveLaw<Node<3> > > p_clone(new GroutingMortar());
        return p_clone;
    }
    
    //*********************************************************************
    //*********************************************************************
    //******* Initialization of the coonstitutive law at the start of the simulation
    //*********************************************************************
    //*********************************************************************
    void GroutingMortar::InitializeMaterial( const Properties& props,
                                        const GeometryType& geom,
                                        const Vector& ShapeFunctionsValues)
    {
        mMaterialParameters = props[MATERIAL_PARAMETERS];
        //Initialize some member variables    
        mCurrentStress.resize(6,false);
        noalias(mCurrentStress)= ZeroVector(6);
        mInsituStress.resize(6,false);
        noalias(mInsituStress)= ZeroVector(6);

        mEpsilon_n.resize(6);
        mEpsilon_n= ZeroVector(6);
        mEpsilon_t_n.resize(6);
        mEpsilon_t_n= ZeroVector(6);
        mEpsilon_current.resize(6);
        mEpsilon_current= ZeroVector(6);

        mCurrentTime= 0.0;
        //get the Material parameters defined in the Pre-Processing        
        mE = mMaterialParameters[0];
        mNU = mMaterialParameters[1];
        //calculate Lame's parameters
        double lambda= mNU*mE/((1+mNU)*(1-2*mNU));
        double mu= mE/(2*(1+mNU));

        mC.resize(6,6);
        //calculate Elastic Matrix acc. to Hooks Law after 28days
        mC_28.resize(6,6);
        mC_28(0,0)=2*mu+lambda; mC_28(0,1)=lambda; mC_28(0,2)=lambda; mC_28(0,3)=0.0; mC_28(0,4)=0.0; mC_28(0,5)=0.0;
        mC_28(1,0)=lambda; mC_28(1,1)=2*mu+lambda; mC_28(1,2)=lambda; mC_28(1,3)=0.0; mC_28(1,4)=0.0; mC_28(1,5)=0.0;
        mC_28(2,0)=lambda; mC_28(2,1)=lambda; mC_28(2,2)=2*mu+lambda; mC_28(2,3)=0.0; mC_28(2,4)=0.0; mC_28(2,5)=0.0;
        mC_28(3,0)=0.0; mC_28(3,1)=0.0; mC_28(3,2)=0.0; mC_28(3,3)=mu; mC_28(3,4)=0.0; mC_28(3,5)=0.0;
        mC_28(4,0)=0.0; mC_28(4,1)=0.0; mC_28(4,2)=0.0; mC_28(4,3)=0.0; mC_28(4,4)=mu; mC_28(4,5)=0.0;
        mC_28(5,0)=0.0; mC_28(5,1)=0.0; mC_28(5,2)=0.0; mC_28(5,3)=0.0; mC_28(5,4)=0.0; mC_28(5,5)=mu;
        //calculate and store coeffeicients
        mTe= mMaterialParameters[2]/3600.0; //in hours
        mDTe= mMaterialParameters[3]/3600.0; //in hours
        mRatioE1E28= mMaterialParameters[4];
        mAe=(1.0-(28.0-mDTe/24.0)/(1.0-mDTe/24.0)*mRatioE1E28*mRatioE1E28)/((1.0-(28.0-mDTe/24.0)/(1.0-mDTe/24.0))*mRatioE1E28*mRatioE1E28);
        mBe= (672.0-mDTe)*(1.0-mAe);
        double bE= pow((mAe+mBe/(mTe-mDTe)),-0.5);
        double DbEDt= pow((mAe+mBe/(mTe-mDTe)),-1.5)*0.5*mBe*pow(mTe-mDTe,-2.0);
        mCe= (2.0*bE/mTe-DbEDt);
        mDe= (DbEDt*mTe-bE)/(mTe*mTe);
    }

    //**********************************************************************
    //**********************************************************************
        // not used in the exercises
    //**********************************************************************
    //**********************************************************************
        void GroutingMortar::SetValue( const Variable<Matrix >& rVariable, 
                                        const Matrix& Value, const ProcessInfo& rCurrentProcessInfo)
    {
    }
    //**********************************************************************
    //**********************************************************************
        // not used in the exercises
    //**********************************************************************
    //**********************************************************************
        void GroutingMortar::SetValue( const Variable<Vector >& rVariable, 
                                        const Vector& rValue, const ProcessInfo& rCurrentProcessInfo)
    {
    }

        void GroutingMortar::SetValue( const Variable<double>& rThisVariable, 
                            const double& rValue,
                            const ProcessInfo& rCurrentProcessInfo )
        {
        }
    
        void GroutingMortar::SetValue( const Variable<array_1d<double, 3> >& rThisVariable, 
                            const array_1d<double, 3>& rValue, 
                            const ProcessInfo& rCurrentProcessInfo 
                          )
        {
        }
    
    //**********************************************************************
    //**********************************************************************
        // not used in the exercises
    //**********************************************************************
    //**********************************************************************
    Matrix GroutingMortar::GetValue(const Variable<Matrix>& rVariable)
    { 
        return ZeroMatrix(0,0);
    }

    //**********************************************************************
    //**********************************************************************
        // not used in the exercises
    //**********************************************************************
    //**********************************************************************
    Vector GroutingMortar::GetValue(const Variable<Vector>& rVariable)
    { 
                if( rVariable == INTERNAL_VARIABLES )
                {
                        Vector result(5);
                        result[0] = mCurrentStress[0];
                        result[1] = mCurrentStress[1];
                        result[2] = mCurrentStress[2];
                        result[3] = mCurrentXi/(mE*mDeltaTime);
                        result[4] = mCurrentTime;      
                        return result;                             
                }       
                if( rVariable == MATERIAL_PARAMETERS )
                {
                        return mMaterialParameters;
                }
                if( rVariable == INSITU_STRESS )
                {
                        return mInsituStress;
                }
               return ZeroVector(0);
    }
    //**********************************************************************
    //**********************************************************************
        // not used in the exercises
    //**********************************************************************
    //**********************************************************************
    double GroutingMortar::GetValue(const Variable<double>& rVariable)
    { 
        return 0.0;
    }
    //**********************************************************************
    //**********************************************************************
        void GroutingMortar::InitializeSolutionStep( const Properties& props,
                const GeometryType& geom, //this is just to give the array of nodes
                const Vector& ShapeFunctionsValues ,
                const ProcessInfo& CurrentProcessInfo)
    {
        mDeltaTime= CurrentProcessInfo[DELTA_TIME]/(3600.0);
        mCurrentTime+= mDeltaTime;
        mCurrentXi= CalculateXi();

// std::cout<<"grouting parameters: t="<<mCurrentTime<<"DeltaTime="<<mDeltaTime<<" Xi="<<mCurrentXi<<" ratio="<<mCurrentXi/(mE*mDeltaTime)<<std::endl;

    }
    //**********************************************************************
    //**********************************************************************
        //this method is called after each solution step ("time step")
    //**********************************************************************
    //**********************************************************************
    void GroutingMortar::FinalizeSolutionStep( const Properties& props,
                                        const GeometryType& geom, const Vector& ShapeFunctionsValues ,const ProcessInfo& CurrentProcessInfo)
    {
                if( CurrentProcessInfo[CALCULATE_INSITU_STRESS]&& !(CurrentProcessInfo[FIRST_TIME_STEP]))
                {
                        noalias(mInsituStress) -= mCurrentStress;
                        mCurrentTime-= mDeltaTime;
                }

                noalias(mEpsilon_t_n)=mEpsilon_t_n+(1.0-mCurrentXi/(mE*mDeltaTime))*(mEpsilon_current-mEpsilon_n);
                noalias(mEpsilon_n)=mEpsilon_current;
    }
    //**********************************************************************
    //**********************************************************************
    // This is the main method called from outside (inside the element at each quadrature point)
    //**********************************************************************
    //**********************************************************************
    void GroutingMortar::UpdateMaterial(  const Vector& StrainVector,
                                      const Properties& props,
                                      const GeometryType& geom,
                                      const Vector& ShapeFunctionsValues,
                                      const ProcessInfo& CurrentProcessInfo )
    {
    //not needed for the elastic law
    }
    //**********************************************************************
    //**********************************************************************
    // Computes the stress vector according to the actual plastic strain vector
    //**********************************************************************
    //**********************************************************************     
    void GroutingMortar::CalculateStress(const Vector& StrainVector, Vector& StressVector)
    {
                if(StressVector.size() != 6)
                        StressVector.resize(6);
                noalias(StressVector)= ZeroVector(6);

                for(unsigned int i=0; i<6; i++)
                        for(unsigned int k=0; k<6; k++)
                                StressVector(i)+= mC_28(i,k)*(StrainVector(k)-mEpsilon_t_n(k)-(1.0-mCurrentXi/(mE*mDeltaTime))*(StrainVector(k)-mEpsilon_n(k)));

                noalias(mEpsilon_current)= StrainVector;

                noalias(mC)= mCurrentXi/(mE*mDeltaTime)*mC_28;

                noalias(mCurrentStress)= StressVector;

                return;

    }
    //**********************************************************************
    //**********************************************************************
    // Called by the element, gives back the algorithmic tangent matrix calculated in update()
    //**********************************************************************
    //********************************************************************** 
    void GroutingMortar::CalculateConstitutiveMatrix(const Vector& StrainVector, Matrix& rResult)
    {
        //noalias(rResult) = mCtangent;
        noalias(rResult) = mC;
    }
    //**********************************************************************
    //**********************************************************************
    // Called by the element, gives back the algorithmic tangent matrix and the stress vector calculated in update()
    //**********************************************************************
    //**********************************************************************     
    void GroutingMortar::CalculateStressAndTangentMatrix( Vector& StressVector,
                                          const Vector& StrainVector,
                                          Matrix& algorithmicTangent)
    {
        noalias(StressVector) = mCurrentStress;
        noalias(algorithmicTangent) = mC;
        
        return;
    }

    //**********************************************************************

    double GroutingMortar::CalculateXi()
    {
        double result;
        if(mCurrentTime<= mTe)
        {
//             result= mE*(0.1*(mCurrentTime-(mCurrentTime-mDeltaTime))+
//                 mCe*(pow(mCurrentTime,2)-pow((mCurrentTime-mDeltaTime),2))/2.0+mDe*(pow(mCurrentTime,3)-pow((mCurrentTime-mDeltaTime),3))/3.0);
            result= mE*(mCe*(pow(mCurrentTime,2)-pow((mCurrentTime-mDeltaTime),2))/2.0+mDe*(pow(mCurrentTime,3)-pow((mCurrentTime-mDeltaTime),3))/3.0);

        }
        else if(mCurrentTime<= 672.0)
        {
            double current_phi= 
                pow(mAe+mBe/((mCurrentTime-mDeltaTime)-mDTe),-0.5)*mCurrentTime
                +0.5*mBe/(pow((mCurrentTime-mDeltaTime)-mDTe,2.0)*pow(mAe+mBe/((mCurrentTime-mDeltaTime)-mDTe),1.5))*0.5*pow((mDeltaTime),2.0)
                +(0.75*pow(mBe,2.0)/(pow((mCurrentTime-mDeltaTime)-mDTe,4.0)*pow(mAe+mBe/((mCurrentTime-mDeltaTime)-mDTe),2.5))-mBe/(pow((mCurrentTime-mDeltaTime)-mDTe,3.0)*pow(mAe+mBe/((mCurrentTime-mDeltaTime)-mDTe),1.5))/2.0)*1.0/3.0*pow((mDeltaTime),3.0);

            if(mCurrentTime> mTe && (mCurrentTime-mDeltaTime)< mTe)
            {
                double te_phi= 
                    pow(mAe+mBe/((mCurrentTime-mDeltaTime)-mDTe),-0.5)*mTe
                    +0.5*mBe/(pow((mCurrentTime-mDeltaTime)-mDTe,2.0)*pow(mAe+mBe/((mCurrentTime-mDeltaTime)-mDTe),1.5))*0.5*pow(mTe-(mCurrentTime-mDeltaTime),2.0)
                    +(0.75*pow(mBe,2.0)/(pow((mCurrentTime-mDeltaTime)-mDTe,4.0)*pow(mAe+mBe/((mCurrentTime-mDeltaTime)-mDTe),2.5))-mBe/(pow((mCurrentTime-mDeltaTime)-mDTe,3.0)*pow(mAe+mBe/((mCurrentTime-mDeltaTime)-mDTe),1.5))/2.0)*1.0/3.0*pow(mTe-(mCurrentTime-mDeltaTime),3.0);
                result= mE*(0.1*(mTe-(mCurrentTime-mDeltaTime))+mCe*(pow(mTe,2)-pow((mCurrentTime-mDeltaTime),2))/2.0+mDe*(pow(mTe,3)-pow((mCurrentTime-mDeltaTime),3))/3.0+current_phi-te_phi); 
            }
            else
            {
                double old_phi= 
                    pow(mAe+mBe/((mCurrentTime-mDeltaTime)-mDTe),-0.5)*(mCurrentTime-mDeltaTime);
                result= mE*(current_phi-old_phi);
            }
        }
        else if(mCurrentTime> 672.0)
        {
            double d28_phi= 
            pow(mAe+mBe/((mCurrentTime-mDeltaTime)-mDTe),-0.5)*672.0
            +0.5*mBe/(pow((mCurrentTime-mDeltaTime)-mDTe,2.0)*pow(mAe+mBe/((mCurrentTime-mDeltaTime)-mDTe),1.5))*0.5*pow(672.0-(mCurrentTime-mDeltaTime),2.0)
            +(0.75*pow(mBe,2.0)/(pow((mCurrentTime-mDeltaTime)-mDTe,4.0)*pow(mAe+mBe/((mCurrentTime-mDeltaTime)-mDTe),2.5))-mBe/(pow((mCurrentTime-mDeltaTime)-mDTe,3.0)*pow(mAe+mBe/((mCurrentTime-mDeltaTime)-mDTe),1.5))/2.0)*1.0/3.0*pow(672.0-(mCurrentTime-mDeltaTime),3.0);
            double old_phi= 
            pow(mAe+mBe/((mCurrentTime-mDeltaTime)-mDTe),-0.5)*(mCurrentTime-mDeltaTime);

            if(mCurrentTime> mTe && (mCurrentTime-mDeltaTime)< mTe)
                result= mE*(d28_phi-old_phi+mCurrentTime-672.0); 
            else
                result= mE*mDeltaTime;
        }

        return result;
    }
} // Namespace Kratos
