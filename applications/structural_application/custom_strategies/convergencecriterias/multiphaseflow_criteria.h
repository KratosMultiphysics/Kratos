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
*   Date:                $Date: 2007-04-13 15:59:32 $
*   Revision:            $Revision: 1.2 $
*
* ***********************************************************/


#if !defined(KRATOS_MULTIPHASEFLOW_CRITERIA )
#define  KRATOS_MULTIPHASEFLOW_CRITERIA


/* System includes */


/* External includes */


/* Project includes */
#include "includes/model_part.h"
#include "includes/define.h"

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
    class TDenseSpace 
            >
            class MultiPhaseFlowCriteria : public ConvergenceCriteria< TSparseSpace, TDenseSpace >
    {
        public:
            /**@name Type Definitions */       
            /*@{ */

            typedef boost::shared_ptr< MultiPhaseFlowCriteria< TSparseSpace, TDenseSpace > > Pointer;		

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
            MultiPhaseFlowCriteria(
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
            virtual ~MultiPhaseFlowCriteria(){}


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
            const TSystemVectorType& b
                             )
            {
                if (Dx.size() != 0) //if we are solving for something
                {
                    
                    double mFinalCorrectionNorm = 0.0;
                    double EnergyNorm= 0.0;
                    double referenceNorm =0.0;
                    int counter= 0;
                    double mFinalCorrectionNorm_WATER = 0.0;
                    double EnergyNorm_WATER= 0.0;
                    double referenceNorm_WATER =0.0;
                    int counter_WATER= 0;
                    double mFinalCorrectionNorm_AIR = 0.0;
                    double EnergyNorm_AIR= 0.0;
                    double referenceNorm_AIR =0.0;
                    int counter_AIR= 0;
                    
                    bool HasWaterPres = false;
                    bool HasAirPres = false;

                    for(typename DofsArrayType::iterator i_dof = rDofSet.begin() ; i_dof != rDofSet.end() ; ++i_dof)
                        {
                            if(i_dof->IsFree())
                            {
                                if(i_dof->GetVariable().Name()== "DISPLACEMENT_X")
                                {	
                                            mFinalCorrectionNorm +=
                                                Dx[i_dof->EquationId()]*Dx[i_dof->EquationId()];
                                            EnergyNorm += b[i_dof->EquationId()]*b[i_dof->EquationId()];
                                            referenceNorm += 
                                                    i_dof->GetSolutionStepValue(DISPLACEMENT_EINS_X)* i_dof->GetSolutionStepValue(DISPLACEMENT_EINS_X);
                                            counter++;
                                }
                                if( i_dof->GetVariable().Name()== "DISPLACEMENT_Y")
                                {	
                                            mFinalCorrectionNorm +=
                                                Dx[i_dof->EquationId()]*Dx[i_dof->EquationId()];
                                            EnergyNorm += b[i_dof->EquationId()]*b[i_dof->EquationId()];
                                            referenceNorm += 
                                                    i_dof->GetSolutionStepValue(DISPLACEMENT_EINS_Y)* i_dof->GetSolutionStepValue(DISPLACEMENT_EINS_Y);
                                            counter++;
                                }
                                if(i_dof->GetVariable().Name()== "DISPLACEMENT_Z")
                                {	
                                            mFinalCorrectionNorm +=
                                                Dx[i_dof->EquationId()]*Dx[i_dof->EquationId()];
                                            EnergyNorm += b[i_dof->EquationId()]*b[i_dof->EquationId()];
                                            referenceNorm += 
                                                    i_dof->GetSolutionStepValue(DISPLACEMENT_EINS_Z)* i_dof->GetSolutionStepValue(DISPLACEMENT_EINS_Z);
                                            counter++;
                                }
                                if(i_dof->GetVariable().Name()== "WATER_PRESSURE")
                                {	
                                    HasWaterPres= true;
                                            
                                    mFinalCorrectionNorm_WATER +=
                                            Dx[i_dof->EquationId()]*Dx[i_dof->EquationId()];
                                    EnergyNorm_WATER  += b[i_dof->EquationId()]*b[i_dof->EquationId()];
                                    referenceNorm_WATER  += 
                                            i_dof->GetSolutionStepValue(WATER_PRESSURE_EINS)*
						 					i_dof->GetSolutionStepValue(WATER_PRESSURE_EINS);
                                    counter_WATER ++;
                                }
                                if(i_dof->GetVariable().Name()== "AIR_PRESSURE")
                                {	
                                    HasAirPres= true;
                                    
                                    mFinalCorrectionNorm_AIR +=
                                            Dx[i_dof->EquationId()]*Dx[i_dof->EquationId()];
                                    EnergyNorm_AIR  += b[i_dof->EquationId()]*b[i_dof->EquationId()];
                                    referenceNorm_AIR  += 
                                            i_dof->GetSolutionStepValue(AIR_PRESSURE_EINS)* i_dof->GetSolutionStepValue(AIR_PRESSURE_EINS);
                                    counter_AIR ++;
                                }
                            }
                        }

                        referenceNorm = sqrt(referenceNorm/counter);
                        mFinalCorrectionNorm = sqrt(mFinalCorrectionNorm/counter);
                        EnergyNorm= sqrt(EnergyNorm/counter);
                        
                        if(HasWaterPres)
                        {
                            referenceNorm_WATER = sqrt(referenceNorm_WATER/counter_WATER);
                            mFinalCorrectionNorm_WATER = 
                                    sqrt(mFinalCorrectionNorm_WATER/counter_WATER);
                            EnergyNorm_WATER= sqrt(EnergyNorm_WATER/counter_WATER);
                        }
                        
                        if(HasAirPres)
                        {
                            referenceNorm_AIR = sqrt(referenceNorm_AIR/counter_AIR);
                            mFinalCorrectionNorm_AIR = 
                                    sqrt(mFinalCorrectionNorm_AIR/counter_AIR);
                            EnergyNorm_AIR= sqrt(EnergyNorm_AIR/counter_AIR);
                        }

                double ratioDisp=1.0;
                
                double ratioDisp_WATER=1.0;
                
                double ratioDisp_AIR=1.0;
                            
                if(referenceNorm >0)
                    ratioDisp =mFinalCorrectionNorm/referenceNorm;
                
                if(referenceNorm_WATER >0)
                    ratioDisp_WATER =mFinalCorrectionNorm_WATER/referenceNorm_WATER;
                
                if(referenceNorm_AIR >0)
                    ratioDisp_AIR =mFinalCorrectionNorm_AIR/referenceNorm_AIR;

                std::cout << "********************************************CONVERGENCE CRITERIA FOR MULTIPHASE PROBLEMS********************************************" <<std::endl;
				std::cout.precision(3);
				std::cout.setf(std::ios::scientific);
                std::cout <<"** expected values: \t\t\t\t\t\tchange= " << mAlwaysConvergedNorm <<"\t\t\t\tenergy= "<<mRatioTolerance<< " **"<<std::endl;
                std::cout  <<"** obtained values displascement:\tratio= "<<ratioDisp <<"\tchange= "<< mFinalCorrectionNorm<<"\tabsolute= "<< referenceNorm << "\tenergy= "<< EnergyNorm <<" **"<<std::endl;
                if(HasWaterPres)
				{
                    std::cout <<"** obtained values water pressure:\tratio= "<< ratioDisp_WATER <<"\tchange= "<< mFinalCorrectionNorm_WATER <<"\tabsolute= "<< referenceNorm_WATER << " \tenergy= "<< EnergyNorm_WATER <<" **"<<std::endl;    
                	if(HasAirPres)
					{
                    		std::cout <<"** obtained values air pressure:\tratio= "<< ratioDisp_AIR <<"\tchange= "<< mFinalCorrectionNorm_AIR<<"\tabsolute= "<< referenceNorm_AIR << "\tenergy= "<< EnergyNorm_AIR <<" **"<<std::endl;

			        		std::cout <<"** obtained values total:\t\tratio= "<< ratioDisp_AIR+ratioDisp_WATER+ratioDisp <<"\tchange= "<< mFinalCorrectionNorm+mFinalCorrectionNorm_WATER+mFinalCorrectionNorm_AIR<<"\tabsolute= "<< referenceNorm_AIR+referenceNorm_WATER+referenceNorm << "\tenergy= "<< EnergyNorm_AIR+EnergyNorm_WATER+EnergyNorm <<" **"<<std::endl;	
					}
					else
					{
			         	std::cout <<"** obtained values total:\t\tratio= "<< ratioDisp_WATER+ratioDisp <<"\t change= "<< mFinalCorrectionNorm+mFinalCorrectionNorm_WATER<<"\tabsolute= "<< referenceNorm_WATER+referenceNorm << "\tenergy= "<< EnergyNorm_WATER+EnergyNorm <<" **"<<std::endl;	
					}
				}
               std::cout << "************************************************************************************************************************************" <<std::endl;

                if ( 	(EnergyNorm<= mRatioTolerance || mFinalCorrectionNorm<= mAlwaysConvergedNorm) 
					&& 	(EnergyNorm_WATER<= mRatioTolerance || mFinalCorrectionNorm_WATER<= mAlwaysConvergedNorm) 
                    && 	(EnergyNorm_AIR<= mRatioTolerance  || mFinalCorrectionNorm_AIR<= mAlwaysConvergedNorm) )
             	{
                  	std::cout << "Congratulations the time step solution is converged..." << std::endl;
                    return true;
                }
                else
                {
                        return false;
                    }
                }
                else //in this case all the displacements are imposed!
                {
                    return true;
                }

            }

            void Initialize(
                    ModelPart& r_model_part
                           ) 
            {}

            void InitializeSolutionStep(
                    ModelPart& r_model_part,
            DofsArrayType& rDofSet,
            const TSystemMatrixType& A,
            const TSystemVectorType& Dx,
            const TSystemVectorType& b
                                       )
            {
            }

            void FinalizeSolutionStep(
                    ModelPart& r_model_part,
            DofsArrayType& rDofSet,
            const TSystemMatrixType& A,
            const TSystemVectorType& Dx,
            const TSystemVectorType& b
                                     ){}



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
                    array_1d<double,3>& OldDisp = (i)->GetSolutionStepValue(DISPLACEMENT_NULL);
                    array_1d<double,3>& CurrentDisp = (i)->GetSolutionStepValue(DISPLACEMENT);
                    

                    if( (i->pGetDof(DISPLACEMENT_X))->IsFixed() == false )
                        result += (CurrentDisp[0]-OldDisp[0])*(CurrentDisp[0]-OldDisp[0]);
                    if( i->pGetDof(DISPLACEMENT_Y)->IsFixed() == false )
                        result += (CurrentDisp[1]-OldDisp[1])*(CurrentDisp[1]-OldDisp[1]);                        
                    if( i->HasDofFor(DISPLACEMENT_Z))
                        if( i->pGetDof(DISPLACEMENT_Z)->IsFixed() == false )
                            result += (CurrentDisp[2]-OldDisp[2])*(CurrentDisp[2]-OldDisp[2]);
                }
                return result;
            }
                    
            double CalculateNormWater(ModelPart& r_model_part)
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
            
            double CalculateNormAir(ModelPart& r_model_part)
            { 
                double result=0;
                for(ModelPart::NodeIterator i = r_model_part.NodesBegin() ; 
                    i != r_model_part.NodesEnd() ; ++i)
                {
                    if( i->HasDofFor(AIR_PRESSURE))
                    {
                        double& OldAirPressure = (i)->GetSolutionStepValue(AIR_PRESSURE_NULL);
                        double& CurrentAirPressure = (i)->GetSolutionStepValue(AIR_PRESSURE);
                        if( i->pGetDof(AIR_PRESSURE)->IsFixed() == false )
                            result += (CurrentAirPressure-OldAirPressure)*(CurrentAirPressure-OldAirPressure);
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
 
                    array_1d<double,3>& CurrentDisp = (i)->GetSolutionStepValue(DISPLACEMENT);
                    
                    if( (i->pGetDof(DISPLACEMENT_X))->IsFixed() == false )
                        result += (CurrentDisp[0])*(CurrentDisp[0]);
                    if( i->pGetDof(DISPLACEMENT_Y)->IsFixed() == false )
                        result += (CurrentDisp[1])*(CurrentDisp[1]);                        
                    if( i->HasDofFor(DISPLACEMENT_Z))
                        if( i->pGetDof(DISPLACEMENT_Z)->IsFixed() == false )
                            result += (CurrentDisp[2])*(CurrentDisp[2]);
                }
                return result;
            }
                    
            double CalculateRefNormWater(ModelPart& r_model_part)
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
            
            double CalculateRefNormAir(ModelPart& r_model_part)
            { 
                double result=0;
                for(ModelPart::NodeIterator i = r_model_part.NodesBegin() ; 
                    i != r_model_part.NodesEnd() ; ++i)
                {
                    if( i->HasDofFor(AIR_PRESSURE))
                    {

                        double& CurrentAirPressure = (i)->GetSolutionStepValue(AIR_PRESSURE);
                        if( i->pGetDof(AIR_PRESSURE)->IsFixed() == false )
                            result += (CurrentAirPressure)*(CurrentAirPressure);
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

