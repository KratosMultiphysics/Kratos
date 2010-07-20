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
*   Revision:            $Revision: 1.6 $
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
    //*********************************************************************
    //*********************************************************************
    //******* default constructor
    //*********************************************************************
    //*********************************************************************  
    HooksLaw::HooksLaw() 
    : ConstitutiveLaw()
    {
    }
    //*********************************************************************
    //*********************************************************************
    //******* default desstructor
    //*********************************************************************
    //********************************************************************* 
    HooksLaw::~HooksLaw()
    {
    }

    boost::shared_ptr<ConstitutiveLaw> HooksLaw::Clone() const
    {
        boost::shared_ptr<ConstitutiveLaw> p_clone(new HooksLaw());
        return p_clone;
    }
    
    //*********************************************************************
    //*********************************************************************
    //******* Initialization of the coonstitutive law at the start of the simulation
    //*********************************************************************
    //*********************************************************************
    void HooksLaw::InitializeMaterial( const Properties& props,
                                        const GeometryType& geom,
                                        const Vector& ShapeFunctionsValues)
    {
        //Initialize some member variables    
        mCurrentStress.resize(6,false);
        noalias(mCurrentStress)= ZeroVector(6);
        mInsituStress.resize(6,false);
        noalias(mInsituStress)= ZeroVector(6);
        //get the Material parameters defined in the Pre-Processing        
        mMaterialParameters = props[MATERIAL_PARAMETERS];
        mE = mMaterialParameters[0];
        mNU = mMaterialParameters[1];
        
        //calculate Lame's parameters
        double lambda= mNU*mE/((1+mNU)*(1-2*mNU));
        double mu= mE/(2*(1+mNU));
        //calculate Elastic Matrix acc. to Hooks Law
        mC.resize(6,6);

        mC(0,0)=2*mu+lambda; mC(0,1)=lambda; mC(0,2)=lambda; mC(0,3)=0.0; mC(0,4)=0.0; mC(0,5)=0.0;
        mC(1,0)=lambda; mC(1,1)=2*mu+lambda; mC(1,2)=lambda; mC(1,3)=0.0; mC(1,4)=0.0; mC(1,5)=0.0;
        mC(2,0)=lambda; mC(2,1)=lambda; mC(2,2)=2*mu+lambda; mC(2,3)=0.0; mC(2,4)=0.0; mC(2,5)=0.0;
        mC(3,0)=0.0; mC(3,1)=0.0; mC(3,2)=0.0; mC(3,3)=mu; mC(3,4)=0.0; mC(3,5)=0.0;
        mC(4,0)=0.0; mC(4,1)=0.0; mC(4,2)=0.0; mC(4,3)=0.0; mC(4,4)=mu; mC(4,5)=0.0;
        mC(5,0)=0.0; mC(5,1)=0.0; mC(5,2)=0.0; mC(5,3)=0.0; mC(5,4)=0.0; mC(5,5)=mu;
                                    
    }
    
    bool HooksLaw::Has( const Variable<double>& rThisVariable )
    {
        return false;
    }
    
    bool HooksLaw::Has( const Variable<Vector>& rThisVariable )
    {
        if( rThisVariable == INTERNAL_VARIABLES )
            return true;
        return false;
    }
    
    bool HooksLaw::Has( const Variable<Matrix>& rThisVariable )
    {
        return false;
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
    }
    //**********************************************************************
    //**********************************************************************
        // not used in the exercises
    //**********************************************************************
    //**********************************************************************
    //**********************************************************************
    //**********************************************************************
        // not used in the exercises
    //**********************************************************************
    //**********************************************************************
    Vector& HooksLaw::GetValue(const Variable<Vector>& rVariable, Vector& rValue)
    { 
                if( rVariable == INTERNAL_VARIABLES )
                {
                    rValue = mCurrentStress;
                    
                    return( rValue );                               
                }       
               return rValue;
    }
    //**********************************************************************
    //**********************************************************************
        // not used in the exercises
    //**********************************************************************
    //**********************************************************************
    double& HooksLaw::GetValue(const Variable<double>& rVariable, double& rValue)
    { 
        return rValue;
    }
    //**********************************************************************
    //**********************************************************************
        void HooksLaw::InitializeSolutionStep( const Properties& props,
                const GeometryType& geom, //this is just to give the array of nodes
                const Vector& ShapeFunctionsValues ,
                const ProcessInfo& CurrentProcessInfo)
    {
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
                        noalias(mInsituStress) -= mCurrentStress;
                }
    }
    //**********************************************************************
    //**********************************************************************
    // This is the main method called from outside (inside the element at each quadrature point)
    //**********************************************************************
    //**********************************************************************
    void HooksLaw::UpdateMaterial(  const Vector& StrainVector,
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
    void HooksLaw::CalculateStress(const Vector& StrainVector, Vector& StressVector)
    {
        noalias(mCurrentStress)= ZeroVector(6);
        for(unsigned int i=0; i<6; i++)
            for(unsigned int j=0; j<6; j++)
                mCurrentStress(i)+= StrainVector(j)*mC(i,j);
        noalias(StressVector)= mCurrentStress-mInsituStress;
        return;
    }
    //**********************************************************************
    //**********************************************************************
    // Called by the element, gives back the algorithmic tangent matrix calculated in update()
    //**********************************************************************
    //********************************************************************** 
    void HooksLaw::CalculateConstitutiveMatrix(const Vector& StrainVector, Matrix& rResult)
    {
        //noalias(rResult) = mCtangent;
        noalias(rResult) = mC;
    }
    
    void  HooksLaw::CalculateMaterialResponse( const Vector& StrainVector,
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
        CalculateStress(StrainVector, StressVector);
        CalculateConstitutiveMatrix(StrainVector, AlgorithmicTangent);
    }
    //**********************************************************************
    //**********************************************************************
    // Called by the element, gives back the algorithmic tangent matrix and the stress vector calculated in update()
    //**********************************************************************
    //**********************************************************************     
    void HooksLaw::CalculateStressAndTangentMatrix( Vector& StressVector,
                                          const Vector& StrainVector,
                                          Matrix& algorithmicTangent)
    {
        noalias(StressVector) = mCurrentStress;
        noalias(algorithmicTangent) = mC;
        
        return;
    }
} // Namespace Kratos
