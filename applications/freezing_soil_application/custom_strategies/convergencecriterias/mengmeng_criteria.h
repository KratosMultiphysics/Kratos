/*
==============================================================================
KratosR1StructuralApplication 
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
*   Last Modified by:    $Author: mengmeng $
*   Date:                $Date: 2009-02-23 16:03:04 $
*   Revision:            $Revision: 1.3 $
*
* ***********************************************************/


#if !defined(KRATOS_MENGMENG_CRITERIA )
#define  KRATOS_MENGMENG_CRITERIA


/* System includes */


/* External includes */


/* Project includes */
#include "includes/model_part.h"
#include "includes/define.h"
#include "freezing_soil_application.h"

namespace Kratos
{

/**@name Kratos Globals */
/*@{ */


/*@} */
/**@name Type Definitions */       
/*@{ */

/*@} */


/**@name  Enum's */       
/*@{ */


/*@} */
/**@name  Functions */       
/*@{ */



/*@} */
/**@name Kratos Classes */
/*@{ */

	/** Short class definition.
Detail class definition.

\URL[Example of use html]{ extended_documentation/no_ex_of_use.html}

\URL[Example of use pdf]{ extended_documentation/no_ex_of_use.pdf}

\URL[Example of use doc]{ extended_documentation/no_ex_of_use.doc}

\URL[Example of use ps]{ extended_documentation/no_ex_of_use.ps}


\URL[Extended documentation html]{ extended_documentation/no_ext_doc.html}

\URL[Extended documentation pdf]{ extended_documentation/no_ext_doc.pdf}

\URL[Extended documentation doc]{ extended_documentation/no_ext_doc.doc}

\URL[Extended documentation ps]{ extended_documentation/no_ext_doc.ps}


*/
template<class TSparseSpace,
class TDenseSpace > class MengmengCriteria : public ConvergenceCriteria< TSparseSpace, TDenseSpace >
{
public:
/**@name Type Definitions */       
/*@{ */

typedef boost::shared_ptr< MengmengCriteria< TSparseSpace, TDenseSpace > > Pointer;	
typedef ConvergenceCriteria< TSparseSpace, TDenseSpace > BaseType;
typedef typename BaseType::TDataType TDataType;
typedef typename BaseType::DofsArrayType DofsArrayType;
typedef typename BaseType::TSystemMatrixType TSystemMatrixType;
typedef typename BaseType::TSystemVectorType TSystemVectorType;

/*@} */
/**@name Life Cycle 
*/    
/*@{ */

/** Constructor.
*/
MengmengCriteria(
	TDataType NewRatioTolerance,
	TDataType AlwaysConvergedNorm)
: ConvergenceCriteria< TSparseSpace, TDenseSpace >()
{
	mRatioTolerance = NewRatioTolerance;
	mAlwaysConvergedNorm = AlwaysConvergedNorm;

	//mActualizeRHSIsNeeded = false;
}

/** Destructor.
*/
virtual ~MengmengCriteria(){}


/*@} */
/**@name Operators 
*/  
/*@{ */

/*Criterias that need to be called after getting the solution */
bool PostCriteria(
	ModelPart& r_model_part,
	DofsArrayType& rDofSet,
	const TSystemMatrixType& A,
	const TSystemVectorType& Dx,
	const TSystemVectorType& b)
{
	if (Dx.size() != 0) //if we are solving for something
	{
		double mFinalCorrectionNorm_DISP = 0.0;
		double EnergyNorm_DISP = 0.0;
		double referenceNorm_DISP =0.0;
		int counter_DISP = 0;

		double mFinalCorrectionNorm_WATER_PRES = 0.0;
		double EnergyNorm_WATER_PRES= 0.0;
		double referenceNorm_WATER_PRES = 0.0;
		int counter_WATER_PRES= 0;

		double mFinalCorrectionNorm_TEMP = 0.0;
		double EnergyNorm_TEMP = 0.0;
		double referenceNorm_TEMP = 0.0;
		int counter_TEMP= 0;

		bool HasDisp = false;
		bool HasWaterPres = false;
		bool HasTemp = false;

		for(typename DofsArrayType::iterator i_dof = rDofSet.begin() ; i_dof != rDofSet.end() ; ++i_dof)
		{
			if(i_dof->IsFree())
			{	
				if(i_dof->GetVariable().Name()== "DISPLACEMENT_X")
				{	
					HasDisp= true;
						
					mFinalCorrectionNorm_DISP +=
					Dx[i_dof->EquationId()]*Dx[i_dof->EquationId()];
					EnergyNorm_DISP += b[i_dof->EquationId()]*b[i_dof->EquationId()];
					referenceNorm_DISP += 
						i_dof->GetSolutionStepValue(DISPLACEMENT_EINS_X)* i_dof->GetSolutionStepValue(DISPLACEMENT_EINS_X);
					counter_DISP++;
				}
				if( i_dof->GetVariable().Name()== "DISPLACEMENT_Y")
				{	
					HasDisp= true;
						
					mFinalCorrectionNorm_DISP +=
					Dx[i_dof->EquationId()]*Dx[i_dof->EquationId()];
					EnergyNorm_DISP += b[i_dof->EquationId()]*b[i_dof->EquationId()];
					referenceNorm_DISP += 
						i_dof->GetSolutionStepValue(DISPLACEMENT_EINS_Y)* i_dof->GetSolutionStepValue(DISPLACEMENT_EINS_Y);
					counter_DISP++;
				}
				if(i_dof->GetVariable().Name()== "DISPLACEMENT_Z")
				{	
					HasDisp= true;
						
					mFinalCorrectionNorm_DISP +=
					Dx[i_dof->EquationId()]*Dx[i_dof->EquationId()];
					EnergyNorm_DISP += b[i_dof->EquationId()]*b[i_dof->EquationId()];
					referenceNorm_DISP += 
						i_dof->GetSolutionStepValue(DISPLACEMENT_EINS_Z)* i_dof->GetSolutionStepValue(DISPLACEMENT_EINS_Z);
					counter_DISP++;
				}

				if(i_dof->GetVariable().Name()== "WATER_PRESSURE")
				{	
					HasWaterPres= true;
						
					mFinalCorrectionNorm_WATER_PRES += Dx[i_dof->EquationId()]*Dx[i_dof->EquationId()];
					EnergyNorm_WATER_PRES  += b[i_dof->EquationId()]*b[i_dof->EquationId()];
					referenceNorm_WATER_PRES  += i_dof->GetSolutionStepValue(WATER_PRESSURE_EINS)*i_dof->GetSolutionStepValue(WATER_PRESSURE_EINS);
					counter_WATER_PRES ++;
				}

				if(i_dof->GetVariable().Name()== "TEMPERATURE")
				{	
					HasTemp= true;
						
					mFinalCorrectionNorm_TEMP += Dx[i_dof->EquationId()]*Dx[i_dof->EquationId()];
					EnergyNorm_TEMP  += b[i_dof->EquationId()]*b[i_dof->EquationId()];
					referenceNorm_TEMP  += i_dof->GetSolutionStepValue(TEMPERATURE_EINS)*i_dof->GetSolutionStepValue(TEMPERATURE_EINS);
					counter_TEMP ++;
				}
			}
		}

		if(HasDisp)
		{
//		    KRATOS_WATCH(Dx)
//		    KRATOS_WATCH(referenceNorm_DISP)
//		    KRATOS_WATCH(mFinalCorrectionNorm_DISP)
//		    KRATOS_WATCH(counter_DISP)
			referenceNorm_DISP = sqrt(referenceNorm_DISP/counter_DISP);
			mFinalCorrectionNorm_DISP = sqrt(mFinalCorrectionNorm_DISP/counter_DISP);
			EnergyNorm_DISP= sqrt(EnergyNorm_DISP/counter_DISP);
		}

		if(HasWaterPres)
		{
			referenceNorm_WATER_PRES = sqrt(referenceNorm_WATER_PRES/counter_WATER_PRES);
			mFinalCorrectionNorm_WATER_PRES = sqrt(mFinalCorrectionNorm_WATER_PRES/counter_WATER_PRES);
			EnergyNorm_WATER_PRES= sqrt(EnergyNorm_WATER_PRES/counter_WATER_PRES);
		}
		
		if(HasTemp)
		{
			referenceNorm_TEMP = sqrt(referenceNorm_TEMP/counter_TEMP);
			mFinalCorrectionNorm_TEMP = sqrt(mFinalCorrectionNorm_TEMP/counter_TEMP);
			EnergyNorm_TEMP= sqrt(EnergyNorm_TEMP/counter_TEMP);
		}
		
		double ratio_DISP=1.0;
		double ratio_WATER_PRES=1.0;
		double ratio_TEMP=1.0;
		
// 		//---------- START 
// 		if(referenceNorm_DISP >0)
// 			ratio_DISP =mFinalCorrectionNorm_DISP/referenceNorm_DISP;
//                 
//                 if(referenceNorm_WATER_PRES >0)
// 			ratio_WATER_PRES =mFinalCorrectionNorm_WATER_PRES/referenceNorm_WATER_PRES;
//                 
//                 if(referenceNorm_TEMP >0)
// 			ratio_TEMP =mFinalCorrectionNorm_TEMP/referenceNorm_TEMP;		
// 		
// 		std::cout << "********************************************CONVERGENCE CRITERIA FOR FREEZING SOIL PROBLEMS********************************************" <<std::endl;
// 		std::cout.precision(3);
// 		std::cout.setf(std::ios::scientific);
// 		std::cout <<"** expected values: \t\t\t\t\t\t\tchange= " << mAlwaysConvergedNorm <<"\t\t\t\tenergy= "<<mRatioTolerance<< " **"<<std::endl;
// 		if(HasDisp)
// 			std::cout <<"** obtained values Displacement:	\tratio= "<<ratio_DISP <<"\tchange= "<< mFinalCorrectionNorm_DISP<<"\tabsolute= "<< referenceNorm_DISP << "\tenergy= "<< EnergyNorm_DISP <<" **"<<std::endl;
// 		if(HasWaterPres)
// 			std::cout <<"** obtained values WaterPressure:  \tratio= "<< ratio_WATER_PRES <<"\tchange= "<< mFinalCorrectionNorm_WATER_PRES <<"\tabsolute= "<< referenceNorm_WATER_PRES << "\tenergy= "<< EnergyNorm_WATER_PRES <<" **"<<std::endl;    
// 		if(HasTemp)
// 			std::cout <<"** obtained values Temperature:    \tratio= "<< ratio_TEMP <<"\tchange= "<< mFinalCorrectionNorm_TEMP <<"\tabsolute= "<< referenceNorm_TEMP << "\tenergy= "<< EnergyNorm_TEMP <<" **"<<std::endl;    
// 		std::cout << "************************************************************************************************************************************" <<std::endl;
// 	
// 		if ( (EnergyNorm_DISP<= mRatioTolerance || mFinalCorrectionNorm_DISP<= mAlwaysConvergedNorm)  
// 			&& (EnergyNorm_WATER_PRES<= mRatioTolerance || mFinalCorrectionNorm_WATER_PRES<= mAlwaysConvergedNorm)
// 			&& (EnergyNorm_TEMP<= mRatioTolerance || mFinalCorrectionNorm_TEMP<= mAlwaysConvergedNorm)) 
// 		{
// 			std::cout << "Congratulations the time step solution is converged..." << std::endl;
// 			return true;
// 		}
// 		else
// 		{
// 			return false;
// 		}
		//---------- END 

		//---------- START by mm 
		if(referenceNorm_DISP >0)
			ratio_DISP =mFinalCorrectionNorm_DISP/referenceNorm_DISP;
// 		else  	ratio_DISP =mFinalCorrectionNorm_DISP;
		
		if(referenceNorm_WATER_PRES >0)
			ratio_WATER_PRES =mFinalCorrectionNorm_WATER_PRES/referenceNorm_WATER_PRES;
// 		else 	ratio_WATER_PRES =mFinalCorrectionNorm_WATER_PRES;
		
		if(referenceNorm_TEMP >0)
			ratio_TEMP =mFinalCorrectionNorm_TEMP/referenceNorm_TEMP;
// 		else	ratio_TEMP =mFinalCorrectionNorm_TEMP;

		std::cout << "************ CONVERGENCE CRITERIA FOR FREEZING SOIL PROBLEMS ***********" <<std::endl;
		std::cout.precision(3);
		std::cout.setf(std::ios::scientific);
		std::cout <<"expected values:    \t\t\ttol_abs =\t" << mAlwaysConvergedNorm <<";\t\t tol_rel = " << mRatioTolerance <<std::endl;
		if(HasDisp)
			std::cout <<"Displacement:  \t\t\tratio =\t"<< mFinalCorrectionNorm_DISP <<"\t/ "<< referenceNorm_DISP <<"\t\t= "<< ratio_DISP<<std::endl;
		if(HasWaterPres)
			std::cout <<"WaterPressure: \t\t\tratio =\t"<< mFinalCorrectionNorm_WATER_PRES <<"\t/ "<< referenceNorm_WATER_PRES <<"\t\t= "<< ratio_WATER_PRES<<std::endl;    
		if(HasTemp)
			std::cout <<"Temperature:   \t\t\tratio =\t"<< mFinalCorrectionNorm_TEMP <<"\t/ "<< referenceNorm_TEMP <<"\t\t= "<< ratio_TEMP<<std::endl;    
		std::cout << "************************************************************************" <<std::endl;
	
		if ( (!HasDisp || ratio_DISP <= mRatioTolerance || mFinalCorrectionNorm_DISP <= mAlwaysConvergedNorm )  
			&& (!HasWaterPres || ratio_WATER_PRES <= mRatioTolerance || mFinalCorrectionNorm_WATER_PRES <= mAlwaysConvergedNorm )  
			&& (!HasTemp || ratio_TEMP <= mRatioTolerance || mFinalCorrectionNorm_TEMP <= mAlwaysConvergedNorm )  ) 
		{
			std::cout << "Congratulations the time step solution is converged..." << std::endl;
			return true;
		}
		else
		{
			return false;
		}
		//---------- END by mm 
	}
	else //in this case all the displacements are imposed!
	{
		return true;
	}

}

void Initialize(ModelPart& r_model_part) {}

void InitializeSolutionStep(
	ModelPart& r_model_part,
	DofsArrayType& rDofSet,
	const TSystemMatrixType& A,
	const TSystemVectorType& Dx,
	const TSystemVectorType& b){}

void FinalizeSolutionStep(
	ModelPart& r_model_part,
	DofsArrayType& rDofSet,
	const TSystemMatrixType& A,
	const TSystemVectorType& Dx,
	const TSystemVectorType& b){}



/*@} */
/**@name Operations */
/*@{ */


/*@} */  
/**@name Access */
/*@{ */


/*@} */
/**@name Inquiry */
/*@{ */


/*@} */      
/**@name Friends */
/*@{ */


/*@} */

protected:
/**@name Protected static Member Variables */
/*@{ */


/*@} */
/**@name Protected member Variables */
/*@{ */


/*@} */
/**@name Protected Operators*/
/*@{ */


/*@} */
/**@name Protected Operations*/
/*@{ */


/*@} */
/**@name Protected  Access */
/*@{ */


/*@} */     
/**@name Protected Inquiry */
/*@{ */


/*@} */   
/**@name Protected LifeCycle */  
/*@{ */



/*@} */    

private:
/**@name Static Member Variables */
/*@{ */


/*@} */
/**@name Member Variables */
/*@{ */
TDataType mRatioTolerance;
TDataType mAlwaysConvergedNorm;


TDataType mReferenceDispNorm;
/*@} */
/**@name Private Operators*/
/*@{ */

void CalculateReferenceNorm(DofsArrayType& rDofSet)
{
	
}

double CalculateNormDisp(ModelPart& r_model_part)
{ 
	double result=0;
	for(ModelPart::NodeIterator i = r_model_part.NodesBegin() ; 
	i != r_model_part.NodesEnd() ; ++i)
	{
		if( i->HasDofFor(DISPLACEMENT))
		{
			array_1d<double,3>& OldDisp = (i)->GetSolutionStepValue(DISPLACEMENT_NULL);
			array_1d<double,3>& CurrentDisp = (i)->GetSolutionStepValue(DISPLACEMENT);
			
			if( i->pGetDof(DISPLACEMENT_X)->IsFixed() == false )
				result += (CurrentDisp[0]-OldDisp[0])*(CurrentDisp[0]-OldDisp[0]);
			if( i->pGetDof(DISPLACEMENT_Y)->IsFixed() == false )
				result += (CurrentDisp[1]-OldDisp[1])*(CurrentDisp[1]-OldDisp[1]);
			if( i->HasDofFor(DISPLACEMENT_Z))
				if( i->pGetDof(DISPLACEMENT_Z)->IsFixed() == false )
				result += (CurrentDisp[2]-OldDisp[2])*(CurrentDisp[2]-OldDisp[2]);
		}
	}
	return result;
}

double CalculateNormWaterPres(ModelPart& r_model_part)
{ 
	double result=0;
	for(ModelPart::NodeIterator i = r_model_part.NodesBegin() ; 
	i != r_model_part.NodesEnd() ; ++i)
	{
		if( i->HasDofFor(WATER_PRESSURE))
		{
			double& OldWaterPressure = (i)->GetSolutionStepValue(WATER_PRESSURE_NULL);
			double& CurrentWaterPressure = (i)->GetSolutionStepValue(WATER_PRESSURE);
			if( i->pGetDof(WATER_PRESSURE)->IsFixed() == false )
			result += (CurrentWaterPressure-OldWaterPressure)*(CurrentWaterPressure-OldWaterPressure);
		}
	}
	return result;
}

double CalculateNormTemp(ModelPart& r_model_part)
{ 
	double result=0;
	for(ModelPart::NodeIterator i = r_model_part.NodesBegin() ; 
	i != r_model_part.NodesEnd() ; ++i)
	{
		if( i->HasDofFor(TEMPERATURE))
		{
			double& OldTemperature = (i)->GetSolutionStepValue(TEMPERATURE_NULL);
			double& CurrentTemperature = (i)->GetSolutionStepValue(TEMPERATURE);
			if( i->pGetDof(TEMPERATURE)->IsFixed() == false )
			result += (CurrentTemperature-OldTemperature)*(CurrentTemperature-OldTemperature);
		}
	}
	return result;
}	

double CalculateRefNormDisp(ModelPart& r_model_part)
{ 
	double result=0;
	for(ModelPart::NodeIterator i = r_model_part.NodesBegin() ; 
	i != r_model_part.NodesEnd() ; ++i)
	{
		if( i->HasDofFor(DISPLACEMENT))
		{
			array_1d<double,3>& CurrentDisp = (i)->GetSolutionStepValue(DISPLACEMENT);
		
			if( i->pGetDof(DISPLACEMENT_X)->IsFixed() == false )
				result += (CurrentDisp[0])*(CurrentDisp[0]);
			if( i->pGetDof(DISPLACEMENT_Y)->IsFixed() == false )
				result += (CurrentDisp[1])*(CurrentDisp[1]);
			if( i->HasDofFor(DISPLACEMENT_Z))
				if( i->pGetDof(DISPLACEMENT_Z)->IsFixed() == false )
				result += (CurrentDisp[2])*(CurrentDisp[2]);
		}
	}
	return result;
}

double CalculateRefNormWaterPres(ModelPart& r_model_part)
{ 
	double result=0;
	for(ModelPart::NodeIterator i = r_model_part.NodesBegin() ; 
	i != r_model_part.NodesEnd() ; ++i)
	{
		if( i->HasDofFor(WATER_PRESSURE))
		{
			double& CurrentWaterPressure = (i)->GetSolutionStepValue(WATER_PRESSURE);
			if( i->pGetDof(WATER_PRESSURE)->IsFixed() == false )
			result += (CurrentWaterPressure)*(CurrentWaterPressure);
		}
	}
	return result;
}

double CalculateRefNormTemp(ModelPart& r_model_part)
{ 
	double result=0;
	for(ModelPart::NodeIterator i = r_model_part.NodesBegin() ; 
	i != r_model_part.NodesEnd() ; ++i)
	{
		if( i->HasDofFor(TEMPERATURE))
		{
			double& CurrentTemperature = (i)->GetSolutionStepValue(TEMPERATURE);
			if( i->pGetDof(TEMPERATURE)->IsFixed() == false )
			result += (CurrentTemperature)*(CurrentTemperature);
		}
	}
	return result;
}


/*@} */
/**@name Private Operations*/
/*@{ */


/*@} */
/**@name Private  Access */
/*@{ */


/*@} */     
/**@name Private Inquiry */
/*@{ */


/*@} */   
/**@name Un accessible methods */
/*@{ */


/*@} */   

}; /* Class ClassName */

/*@} */

/**@name Type Definitions */       
/*@{ */


/*@} */

}  /* namespace Kratos.*/

#endif /* KRATOS_NEW_DISPLACEMENT_CRITERIA  defined */

