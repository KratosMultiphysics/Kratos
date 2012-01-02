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

// External includes
#include<cmath>

// Project includes
#include "includes/define.h"
#include "constitutive_laws/plasticity_2d.h"
#include "includes/constitutive_law.h"
#include "utilities/math_utils.h"
#include "custom_utilities/sd_math_utils.h"
#include "custom_utilities/tensor_utils.h"
#include "includes/variables.h"
#include "includes/process_info.h"
#include "includes/properties.h"


namespace Kratos
{
    /**
     * TO BE TESTED!!!
     */

    Plasticity2D::Plasticity2D()
            : ConstitutiveLaw()

    {
        KRATOS_ERROR( std::logic_error, "Calling the empty constructor.", "" );
    }

    Plasticity2D::Plasticity2D( FluencyCriteriaPointer FluencyCriteria, PropertiesPointer Property )
            : ConstitutiveLaw()
    {
        mpFluencyCriteria   = FluencyCriteria;
        mpProperties        = Property;
    }

    /**
     * TO BE TESTED!!!
     */
    Plasticity2D::~Plasticity2D()
    {
    }


    bool Plasticity2D::Has( const Variable<double>& rThisVariable )
    {
        return false;
    }

    bool Plasticity2D::Has( const Variable<Vector>& rThisVariable )
    {
        return false;
    }

    bool Plasticity2D::Has( const Variable<Matrix>& rThisVariable )
    {
        return false;
    }

    double& Plasticity2D::GetValue( const Variable<double>& rThisVariable, double& rValue )
    {
        if ( rThisVariable == DAMAGE )
        {
           mpFluencyCriteria->GetValue(rThisVariable, rValue); 
        }
        if ( rThisVariable == YIELD_STRESS)
        {
           mpFluencyCriteria->GetValue(rThisVariable, rValue); 
        }
        return rValue;
    }

    Vector& Plasticity2D::GetValue( const Variable<Vector>& rThisVariable, Vector& rValue )
    {  /*
        if ( rThisVariable == INTERNAL_VARIABLES )
        {
            rValue = ZeroVector( 12 );

            for ( unsigned int i = 0; i < 3; i++ )
                rValue[i] = mplastic_strain[i];
        }
        */
        return( rValue );
    }

    Matrix& Plasticity2D::GetValue( const Variable<Matrix>& rThisVariable, Matrix& rValue )
    {
        return( rValue );
    }

    void Plasticity2D::SetValue( const Variable<double>& rThisVariable, const double& rValue,
                                 const ProcessInfo& rCurrentProcessInfo )
    {
    }

    void Plasticity2D::SetValue( const Variable<Vector>& rThisVariable, const Vector& rValue,
                                 const ProcessInfo& rCurrentProcessInfo )
    {
    }

    void Plasticity2D::SetValue( const Variable<Matrix>& rThisVariable, const Matrix& rValue,
                                 const ProcessInfo& rCurrentProcessInfo )
    {
    }

    void Plasticity2D::Calculate( const Variable<Matrix >& rVariable, Matrix& rResult,
                                  const ProcessInfo& rCurrentProcessInfo )
    {
    }


//***********************************************************************************************
//***********************************************************************************************

    void Plasticity2D::InitializeMaterial( const Properties& props,
                                           const GeometryType& geom,
                                           const Vector& ShapeFunctionsValues )

    {  
          mE   = (*mpProperties)[YOUNG_MODULUS];
	  mNU  = (*mpProperties)[POISSON_RATIO];
	  mDE  = (*mpProperties)[DENSITY];
	  KRATOS_WATCH(mpFluencyCriteria)
          mpFluencyCriteria->InitializeMaterial(*mpProperties); 
	  double he = geom.Length();
	  mpFluencyCriteria->GetValue(he); 
    }


//***********************************************************************************************
//***********************************************************************************************

    void Plasticity2D::InitializeSolutionStep( const Properties& props,
            const GeometryType& geom,
            const Vector& ShapeFunctionsValues ,
            const ProcessInfo& CurrentProcessInfo )
    {
    }

//***********************************************************************************************
//***********************************************************************************************

    void Plasticity2D::CalculateConstitutiveMatrix( const Vector& StrainVector, Matrix& ConstitutiveMatrix )
    {
        Vector StressVector(3);
        CalculateStressAndTangentMatrix( StressVector, StrainVector, ConstitutiveMatrix );
    }



//***********************************************************************************************
//***********************************************************************************************


    void  Plasticity2D::FinalizeSolutionStep( const Properties& props,
            const GeometryType& geom,
            const Vector& ShapeFunctionsValues ,
            const ProcessInfo& CurrentProcessInfo )

    {
      mpFluencyCriteria->FinalizeSolutionStep(); 
    }

//***********************************************************************************************
//***********************************************************************************************

    void Plasticity2D::CalculateElasticStress( const Vector& Strain, array_1d<double, 4>& Stress)
    {
      // plane strain and axial symmetric 
      double G            = 0.5*mE / (1.00 + mNU);
      double K            = mE / (3.00 * (1.00-2.00*mNU) );
      double vol_strain   = Strain[0] + Strain[1] + Strain[3];
      double pt           = K * vol_strain;     
      double vol_strain3  = vol_strain * 0.333333333333333333; 
      
      Stress[0] = 2.00 * G *(Strain[0] - vol_strain3 ) + pt; 
      Stress[1] = 2.00 * G *(Strain[1] - vol_strain3 ) + pt;
      Stress[2] = G *(Strain[2]);
      Stress[3] = 2.00 * G *(Strain[3] - vol_strain3 ) + pt;
    }


//***********************************************************************************************
//***********************************************************************************************

    void Plasticity2D::CalculateStress( const Vector& StrainVector, Vector& StressVector )

    {
        KRATOS_TRY
	mpFluencyCriteria->ReturnMapping(StrainVector, StressVector);
        KRATOS_CATCH( "" )
    }


//***********************************************************************************************
//***********************************************************************************************

    void Plasticity2D::CalculateCauchyStresses(
        Vector& rCauchy_StressVector,
        const Matrix& rF,
        const Vector& rPK2_StressVector,
        const Vector& rGreenLagrangeStrainVector )
    {
        Matrix S = MathUtils<double>::StressVectorToTensor( rPK2_StressVector );
        double J = MathUtils<double>::Det2( rF );
        boost::numeric::ublas::bounded_matrix<double, 2, 2> temp;
        boost::numeric::ublas::bounded_matrix<double, 2, 2> aux;
        noalias( temp ) = prod( rF, S );
        noalias( aux )  = prod( temp, trans( rF ) );
        aux *= J;
        if ( rCauchy_StressVector.size() != 3 )
            rCauchy_StressVector.resize( 3 );

        rCauchy_StressVector[0] = aux( 0, 0 );
        rCauchy_StressVector[1] = aux( 1, 1 );
        rCauchy_StressVector[2] = aux( 0,1 );
    }


//***********************************************************************************************
//***********************************************************************************************
    void Plasticity2D::CalculateStressAndTangentMatrix( Vector& StressVector,
            const Vector& StrainVector,
            Matrix& AlgorithmicTangent )
    {  
          // plane strain and axial symmetric  
          AlgorithmicTangent.resize(3,3, false);
          double  c1 = mE*(1.00-mNU)/((1.00 + mNU)*(1.00-2.00*mNU));
          double  c2 = c1*mNU/(1.00 - mNU);
          double  c3 = 0.5*mE / (1.00 + mNU);
          AlgorithmicTangent(0,0) = c1;  AlgorithmicTangent(0,1) = c2;  AlgorithmicTangent(0,2) = 0.0;
          AlgorithmicTangent(1,0) = c2;  AlgorithmicTangent(1,1) = c1;  AlgorithmicTangent(1,2) = 0.0;
          AlgorithmicTangent(2,0) = 0.0; AlgorithmicTangent(2,1) = 0.0; AlgorithmicTangent(2,2) = c3;
    }




//***********************************************************************************************
//***********************************************************************************************

    void Plasticity2D::UpdateMaterial( const Vector& StrainVector,
                                       const Properties& props,
                                       const GeometryType& geom,
                                       const Vector& ShapeFunctionsValues,
                                       const ProcessInfo& CurrentProcessInfo )
    {
       mpFluencyCriteria->UpdateMaterial(); 
    }


//***********************************************************************************************
//***********************************************************************************************


    void Plasticity2D::Calculate( const Variable<double>& rVariable,
                                  double& Output,
                                  const ProcessInfo& rCurrentProcessInfo )
    {
    }
    
//***********************************************************************************************
//*********************************************************************************************** 

     void Plasticity2D::CalculateMaterialResponse( const Vector& StrainVector,
        const Matrix& DeformationGradient,
        Vector& StressVector,
        Matrix& AlgorithmicTangent,
        const ProcessInfo& CurrentProcessInfo,
        const Properties& props, 
        const GeometryType& geom,
        const Vector& ShapeFunctionsValues,
        bool CalculateStresses,
        int CalculateTangent,
        bool SaveInternalVariables)
        
        {
         UpdateMaterial(StrainVector, props, geom,ShapeFunctionsValues, CurrentProcessInfo);
         if (CalculateStresses==true){CalculateStress(StrainVector, StressVector);}
         if(CalculateTangent==1){CalculateStressAndTangentMatrix(StressVector,StrainVector, AlgorithmicTangent);}
	}



        int Plasticity2D::Check(const Properties& props,
                const GeometryType& geom,
                const ProcessInfo& CurrentProcessInfo)
        {
            KRATOS_TRY

            if(DENSITY.Key() == 0 || props[DENSITY]<0.00)
                KRATOS_ERROR(std::invalid_argument,"DENSITY has Key zero or invalid value ","");
            
            if(YOUNG_MODULUS.Key() == 0 || props[YOUNG_MODULUS]<= 0.00)
                KRATOS_ERROR(std::invalid_argument,"YOUNG_MODULUS has Key zero or invalid value ","");
	    
	    const double& nu = props[POISSON_RATIO];
	    const bool check = bool( (nu >0.499 && nu<0.501 ) || (nu < -0.999 && nu > -1.01 ) );
	    if(POISSON_RATIO.Key() == 0 || check==true) // props[POISSON_RATIO] == 1.00 || props[POISSON_RATIO] == -1.00)
                KRATOS_ERROR(std::invalid_argument,"POISSON_RATIO has Key zero invalid value ","");
         
            if(YIELD_STRESS.Key() == 0 || props[YIELD_STRESS]< 0.00)
                KRATOS_ERROR(std::invalid_argument,"YIELD_STRESS has Key zero or invalid value ","");

	    if(ISOTROPIC_HARDENING_MODULUS.Key())
                KRATOS_ERROR(std::invalid_argument,"ISOTROPIC_HARDENING_MODULUS has Key zero or invalid value ","");
	    
	    return 0;
	    
            KRATOS_CATCH("");
        }
    

}


