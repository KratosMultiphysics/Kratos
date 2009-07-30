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
Nota: This constitutive law use the serial paralell mixing theory.

* ***********************************************************/


// System includes 
#include <iostream>

// External includes 
#include<cmath>

// Project includes 
#include "includes/define.h"
#include "constitutive_laws/compose_material.h"
#include "includes/constitutive_law.h"
#include "utilities/math_utils.h"
#include "custom_utilities/sd_math_utils.h"
#include "custom_utilities/tensor_utils.h"
#include "includes/variables.h"
#include "includes/process_info.h"
#include "includes/properties.h"


namespace Kratos
{
    namespace Compose_Material_Auxiliaries
    {
	  Vector StressVectorAux; 

    }


      using namespace Compose_Material_Auxiliaries;

      /**
      *	TO BE TESTED!!!
      */
     ComposeMaterial::ComposeMaterial() : ConstitutiveLaw< Node<3> >()
      {}

      ComposeMaterial::ComposeMaterial(MaterialsContainer Mat) 
      : ConstitutiveLaw< Node<3> >()
      {
	mMaterials = Mat;
        mNumMat    = mMaterials.size();    
        KRATOS_WATCH(mNumMat )
        
      }

      /**
      *	TO BE TESTED!!!
      */
      ComposeMaterial::~ComposeMaterial()
      {
      }


      bool ComposeMaterial::Has( const Variable<double>& rThisVariable )
      {
      return false;
      }

      bool ComposeMaterial::Has( const Variable<Vector>& rThisVariable )
      {
      return false;
      }

      bool ComposeMaterial::Has( const Variable<Matrix>& rThisVariable )
      {
      return false;
      }


      double ComposeMaterial::GetValue( const Variable<double>& rThisVariable )
      {
      KRATOS_ERROR(std::logic_error, "double Variable case not considered", "");
      }


      Vector ComposeMaterial::GetValue( const Variable<Vector>& rThisVariable )
      {
      KRATOS_ERROR(std::logic_error, "Vector Variable case not considered", "");
      }

      Matrix ComposeMaterial::GetValue( const Variable<Matrix>& rThisVariable )
      {
      KRATOS_ERROR(std::logic_error, "Matrix Variable case not considered", "");
      }

      void ComposeMaterial::SetValue( const Variable<double>& rThisVariable, const double& rValue, 
      const ProcessInfo& rCurrentProcessInfo )
      {
      }

      void ComposeMaterial::SetValue( const Variable<Vector>& rThisVariable, const Vector& rValue, 
      const ProcessInfo& rCurrentProcessInfo )
      {
      }

      void ComposeMaterial::SetValue( const Variable<Matrix>& rThisVariable, const Matrix& rValue, 
      const ProcessInfo& rCurrentProcessInfo )
      {
      }

      void ComposeMaterial::Calculate(const Variable<Matrix >& rVariable, Matrix& rResult, 
      const ProcessInfo& rCurrentProcessInfo)
      {
      }


      //***********************************************************************************************
      //***********************************************************************************************

      void ComposeMaterial::InitializeMaterial( const Properties& props,
      const GeometryType& geom,
      const Vector& ShapeFunctionsValues )
      {
    
	    
	    for (unsigned int i=0; i< mNumMat; i++ )
		{
		     mMaterials[i]->InitializeMaterial(props,geom,ShapeFunctionsValues);
		}
      
      }


      //***********************************************************************************************
      //***********************************************************************************************

      void ComposeMaterial::InitializeSolutionStep( const Properties& props,
      const GeometryType& geom,
      const Vector& ShapeFunctionsValues ,
      const ProcessInfo& CurrentProcessInfo)
      {
	  for (unsigned int i=0; i< mNumMat; i++ )
		{
		  mMaterials[i]->InitializeSolutionStep(props,geom,ShapeFunctionsValues,CurrentProcessInfo); 
		}

      }

      //***********************************************************************************************
      //***********************************************************************************************


      void ComposeMaterial::CalculateConstitutiveMatrix(const Vector& StrainVector, Matrix& ConstitutiveMatrix)
      {

	  for (unsigned int i=0; i< mNumMat; i++ )
		{
		  mMaterials[i]->CalculateConstitutiveMatrix(StrainVector, ConstitutiveMatrix); 
		}

      }


      //***********************************************************************************************
      //***********************************************************************************************


      void  ComposeMaterial::FinalizeSolutionStep( const Properties& props,
      const GeometryType& geom, 
      const Vector& ShapeFunctionsValues ,
      const ProcessInfo& CurrentProcessInfo)

      {
	  for (unsigned int i=0; i< mNumMat; i++ )
		{
		  mMaterials[i]->FinalizeSolutionStep(props,geom,ShapeFunctionsValues,CurrentProcessInfo); 
		}


      }

      //***********************************************************************************************
      //***********************************************************************************************

      void ComposeMaterial::CalculateStress( const Vector& StrainVector, Vector& StressVector)

      {
	StressVectorAux.resize(StrainVector.size());  
	for (unsigned int i=0; i< mNumMat; i++ )
		{
		  mMaterials[i]->CalculateStress(StrainVector, StressVectorAux);
                  StressVector = StressVectorAux + StressVector;  
		}
      }

      //***********************************************************************************************
      //***********************************************************************************************

      void ComposeMaterial::CalculateCauchyStresses(
      Vector& rCauchy_StressVector,
      const Matrix& rF,
      const Vector& rPK2_StressVector,
      const Vector& rGreenLagrangeStrainVector)
      {

	for (unsigned int i=0; i< mNumMat; i++ )
		{
		  mMaterials[i]->CalculateCauchyStresses(rCauchy_StressVector,rF,rPK2_StressVector,rGreenLagrangeStrainVector);
		}
      }

      //***********************************************************************************************
      //***********************************************************************************************	
      void ComposeMaterial::CalculateStressAndTangentMatrix(Vector& StressVector,
      const Vector& StrainVector,
      Matrix& algorithmicTangent)
      {
	for (unsigned int i=0; i< mNumMat; i++)
		{
		  mMaterials[i]->CalculateStressAndTangentMatrix(StressVector,StrainVector,algorithmicTangent); 
		}

      }

      }



