/*
==============================================================================
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain

Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNER.

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
*   Last Modified by:    $Author: nelson $
*   Date:                $Date: 2009- Mayo-06 $
*   Revision:            $Revision: 1.20 $
*
* ***********************************************************/


#if !defined(KRATOS_LINE_SEARCHES_UTILITY)
#define  KRATOS_LINE_SEARCHES_UTILITY


/* System includes */


/* External includes */
#include "boost/smart_ptr.hpp"



/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/strategies/solving_strategy.h"
#include "solving_strategies/convergencecriterias/convergence_criteria.h"

//default builder and solver
#include "solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver.h"

#include <cmath>

namespace Kratos
{
	template<class TSparseSpace,
	class TDenseSpace,  	
	class TLinearSolver 
	>
	class LineSearchesUtility
			  
	{
	   
		  public:

			typedef ConvergenceCriteria<TSparseSpace,TDenseSpace> TConvergenceCriteriaType;

			typedef SolvingStrategy<TSparseSpace,TDenseSpace,TLinearSolver> BaseType;

			typedef typename BaseType::TBuilderAndSolverType TBuilderAndSolverType;

			typedef typename BaseType::TDataType TDataType;

			typedef TSparseSpace SparseSpaceType;
	  
			typedef typename BaseType::TSchemeType TSchemeType;

			typedef typename BaseType::DofsArrayType DofsArrayType;

			typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

			typedef typename BaseType::TSystemVectorType TSystemVectorType;

			typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

			typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

			typedef typename BaseType::TSystemMatrixPointerType TSystemMatrixPointerType;

			typedef typename BaseType::TSystemVectorPointerType TSystemVectorPointerType;

			/** Counted pointer of ClassName */
			/*KRATOS_CLASS_POINTER_DEFINITION(LineSearchesUtility);
		
			typedef Dof<TDataType> TDofType;
			typedef PointerVectorSet<TDofType, IdentityFunction<TDofType> > DofsArrayType;
			typedef typename PointerVectorSet<TDofType, IndexedObject>::iterator DofIterator;
			typedef typename PointerVectorSet<TDofType, IndexedObject>::const_iterator DofConstantIterator;
			*/


			LineSearchesUtility( 
			//typename TSchemeType::Pointer pScheme,
			//typename TLinearSolver::Pointer pNewLinearSolver,
			unsigned  int& MaxLineSearchIterations,
                        double& tolls,         // energy tolerance factor on LineSearch (0.8 is ok) 
			double& amp,           // maximum amplification factor
			double& etmxa,         // maximum allowed step length
 			double& etmna,       // minimum allowed step length
			bool& MoveMeshFlag,
                        bool& ApplyLineSearches
			)
			 //:rmodel_part(model_part)

			{
			     KRATOS_TRY
			      SetParametersLineSearches(MaxLineSearchIterations,
                                                  tolls,  
						  amp, 
						  etmxa, 
						  etmna, 
						  MoveMeshFlag,       
						  ApplyLineSearches
						  );

                              
			     
			     KRATOS_CATCH("")
			}


		/** Destructor.
		*/
		virtual ~LineSearchesUtility () {}
			

		   protected:

		   // parameters of line searches
		    unsigned int  mMaxLineSearchIterations;
		    double mtolls; 
		    double mamp;
		    double meta;  
		    double metmxa;
		    double metmna; 
		    bool   mMoveMeshFlag; 
                    bool   mApplyLineSearches;
  		
                    // methods of line searches
		    void SetParametersLineSearches(unsigned int& MaxLineSearchIterations,
                                                  double& tolls,  
						  double& amp, 
						  double& etmxa, 
						  double& etmna, 
						  bool& MoveMeshFlag,       
						  bool& ApplyLineSearches)
			{
			    KRATOS_TRY
				mMaxLineSearchIterations = MaxLineSearchIterations;
                                mtolls  = tolls;  
				mamp    = amp; 
				metmxa  = etmxa;  
				metmna  = etmna; 
				mMoveMeshFlag      = MoveMeshFlag;
				mApplyLineSearches = ApplyLineSearches;
				meta= 1.00;
			    KRATOS_CATCH("")
			}


		

		  bool LineSearches( ModelPart& rmodel_part,
				     typename TSchemeType::Pointer& pScheme, 
				     typename TBuilderAndSolverType::Pointer& pBuilderAndSolver,
				     DofsArrayType& rDofSet,
				     TSystemVectorType&  X_old,
				     TSystemVectorType&  Delta_p,
				     TSystemVectorType&  mDx,
				     TSystemVectorType&  mb,
				     TSystemMatrixType&  mA
				    )
		    {

			  KRATOS_TRY
			  double seta               = 0.00;
                          double so                 = 0.00;
                          unsigned int ils          = 0; 
			  unsigned int ico          = 0;
			  unsigned iteration_number = 1;
			  TSystemVectorType prodr(mMaxLineSearchIterations+1); 
			  TSystemVectorType lseta(mMaxLineSearchIterations+2);
			  TSparseSpace::SetToZero(prodr);
			  TSparseSpace::SetToZero(lseta);

			  lseta(1) =  1.00;
			  prodr(0) =  1.00;
                          meta     =  1.00;

			  //KRATOS_WATCH(Delta_p)
			  //KRATOS_WATCH(mDx)
			  //KRATOS_WATCH(mb)			  
			  so = TSparseSpace::Dot(Delta_p,mb); // mb inicial 
			  KRATOS_WATCH(so)


			  while (iteration_number <= mMaxLineSearchIterations)
			    {

			      
			       pScheme->Update(rmodel_part,rDofSet,mA,mDx,mb);	
			       this->BackupDatabase(rDofSet,mDx);
			       //KRATOS_WATCH(mDx) 		       
			       if(mMoveMeshFlag == true) MoveMesh(rmodel_part);	      		      
			       std::cout<<"Line_Search_Iteration_Number:"<<(iteration_number)<<std::endl;
			       TSparseSpace::SetToZero(mb);
			       pBuilderAndSolver->BuildRHS(pScheme,rmodel_part,mb); //mb para calculat eta
			       TSparseSpace::Copy(Delta_p, mDx);	  
			       seta = TSparseSpace::Dot(mDx,mb); //Delta_p =  Const
			       KRATOS_WATCH(seta)
			       //KRATOS_WATCH(mb)
			      // Posible restriccion en el denominador
			       ils = iteration_number -1;
			       prodr(ils+1) = (seta/so);
                               KRATOS_WATCH(prodr(ils+1)) 
			       if (fabs(prodr(ils+1)) < mtolls)
				{ 
				  KRATOS_WATCH(meta)
				  return true;
				  //break;
				 }
			      else 
				{
				 //Search devuelve lseta(ils+2)
				 SetDatabaseToValue(rDofSet, X_old);
				 if(mMoveMeshFlag == true) MoveMesh(rmodel_part);
				 Search(ils,prodr,lseta,ico);
				 
				 if (ico==2)
				   {
				      //ilfail = 2;
                                      meta = 1.00;
				      std::cout<<"********************************"<<std::endl; 
				      std::cout<<"******Line Search Trouble*******"<<std::endl; 
				      std::cout<<"********************************"<<std::endl;
                                      return false;
				       
				  }
				
				 meta = lseta(ils+2);
				 KRATOS_WATCH(meta)
				 TSparseSpace::Assign(mDx,meta, mDx); //mDx =lseta(iteration_number+1)*mDx
				
				 
				 }
				
				/*if (ico==2)
				   {
				      //ilfail = 2;
                                      meta = 1.00;
				      std::cout<<"********************************"<<std::endl; 
				      std::cout<<"******Line Search Trouble*******"<<std::endl; 
				      std::cout<<"********************************"<<std::endl;
                                      return false;
				      //break; 
				  }*/
				
				if (iteration_number>= mMaxLineSearchIterations) 
				    {
				      this->MaxIterationsExceeded();
				      return false;
				      //break; 
				    }
				
				
				iteration_number++;
				}
				return 0;  
				KRATOS_CATCH("")			      
			    }
                            

		     void SetDatabaseToValue(
			DofsArrayType& rDofSet,
			const TSystemVectorType& X_old
			) 
			{ 
			KRATOS_TRY

				for(typename DofsArrayType::iterator i_dof = rDofSet.begin() ; i_dof != rDofSet.end() ; ++i_dof)
				{
					if(i_dof->IsFree())
					{i_dof->GetSolutionStepValue() = X_old[i_dof->EquationId()];}
				}
			KRATOS_CATCH("")
			 }
 
                  
		     // Permite escribir los desplazamientos antiguos en el X_old
		     void BackupDatabase(
			DofsArrayType const& rDofSet,
			TSystemVectorType& X_old
			) 
			{ 
			KRATOS_TRY

				 for(typename DofsArrayType::const_iterator i_dof = rDofSet.begin() ; i_dof != rDofSet.end() ; ++i_dof)
				{
					    if(i_dof->IsFree())
					     { X_old[i_dof->EquationId()] = i_dof->GetSolutionStepValue();}
				}
			KRATOS_CATCH("")
			 }

//realiza Line-Search local para encontrar la longitud del paso en eta(ils+2)
//eta: historia de longitudes de pasos en la iteracion de L-S (eta(1)=0.,eta(2)=1.)
//prodr: productos r
//amp: factor de amplificacion maximo para el paso de carga en caso de extrapolacion
//etmxa y etmna: longitudes de pasos maximos y minimos permitidos
//ico: 1 si en el paso anterior ha sido necesario usar etmxa o etmna, sino es 0
// al salir ico es 1 si en esta iteracion se ha usado etmxa o etmna y 2 si se han 
//usado por segunda vez.

//obtener ineg=No del S-L actual que este mas cerca del origen y que su prodr sea negativo
//tambien busca el S-L maximo en etmaxp
// si no hay existe un eta con prodr<0 ineg=999
	  
		void Search(unsigned int& ils,TSystemVectorType& prodr,TSystemVectorType& lseta, unsigned int& ico)
			  {
			  KRATOS_TRY

			  double etaint      = 0.00;  
			  double etaneg      = 1e5;
			  double etmxt       = 0.00; 
			  double etmaxp      = 0.00;
			  double etaalt      = 0.00;
			  double etaext      = 0.00;
			  unsigned int ipos  = 0;
			  unsigned int ineg  = 999;
			  //KRATOS_WATCH(ils);
			  
			  for (unsigned int i=0;i<=ils+1; i++)
			      {
				if (lseta(i)>etmaxp)
				  {
				    etmaxp=lseta(i);
				  }
				if (prodr(i)<0 && lseta(i)<=etaneg)
				  {
				    etaneg=lseta(i);
				    ineg=i;
				  }
			      }

			  //Decide si interpola o extrapola
			  if (ineg !=999)
			      {
			  //Interpola
			  //entontrar ipos=No del S-L con prodr positivo mas cercano a ineg pero con S-L mas peque�o
			      ipos = 0;
			      for( unsigned int i=0;i<=ils+1;i++)
				{
				  if (prodr(i)>=0 && lseta(i)<=lseta(ineg) && lseta(i)>=lseta(ipos))
				      {
					ipos=i;
				      }
				}
			      //interpola para encontrar S-L etaint
			      etaint=prodr(ineg)*lseta(ipos)-prodr(ipos)*lseta(ineg);
			      if ((prodr(ineg)-prodr(ipos))==0)
				{
				  ico=2;
				  std::cout<<" Warning: Division By Zero. Line Searches Filed"<<std::endl;
				  return;
				}
			      else
				{
				  etaint=etaint/(prodr(ineg)-prodr(ipos));
				}
			      //alternativamente encuentra etaalt para asegurar un cambio razonable
			      etaalt = lseta(ipos) + 0.2*(lseta(ineg)-lseta(ipos));
			      //toma el maximo entre etaint y etaalt
			      etaint=std::max(etaint,etaalt);
			      //la longitud de paso minimo
			      if (etaint<metmna)
				{
				  etaint=metmna;
				  if (ico==1) 
				  {ico=2;
				   std::cout<<"Min Step Length Reached Twice"<<std::endl;}
				  else if (ico==0)
				  {ico=1;}
				}
			      lseta(ils+2)=etaint;
                              meta = lseta(ils+2);
                              return;
			  }

			  else if (ineg==999)
			  {
		  
			    //extrapola
			    //asigna temporalmente la longitud de paso maxima   
			  etmxt = mamp*etmaxp;
			  if(etmxt>metmxa)
			    {
			    etmxt=metmxa;
			    }
			  //extrapola la longitud de paso entre el actual y el anterior
			  etaext=prodr(ils+1)*lseta(ils)-prodr(ils)*lseta(ils+1);
			  if ((prodr(ils+1)-prodr(ils))==0)
			      {
				ico=2;
				std::cout<<" Warning: Division By Zero. Line Searches Filed"<<std::endl;
				return;
			      }
			  else
			      {
			      etaext=etaext/(prodr(ils+1)-prodr(ils));
			      }
		  
			    lseta(ils+2)=etaext;

			    //se acepta eta si esta dentro de los limites
			    if (etaext<=0 || etaext> etmxt)
				{
				  lseta(ils+2) = etmxt;
                                  meta = lseta(ils+2);
                                  return;
				}
			    if (lseta(ils+2)==metmxa && ico==1) 
			      {
				ico=2;
				std::cout<<"Max Step Length Again"<<std::endl;
				return;
			      }
			    if (lseta(ils+2)==metmxa)
				{
				ico=1;
				}
			    return;
			  }
			
			KRATOS_CATCH("")
			}
		
		      void MaxIterationsExceeded()
			{
			std::cout << "*****************************************************************" << std::endl;
			std::cout << "******* ATTENTION: Max Iterations Line Searches Exceeded ********" << std::endl;
			std::cout << "*****************************************************************" << std::endl;
			}


		 void MoveMesh(ModelPart& rmodel_part)
		    {
			KRATOS_TRY

			for(ModelPart::NodeIterator i = rmodel_part.NodesBegin() ; 
				i != rmodel_part.NodesEnd() ; ++i)
			{
				array_1d<double,3>& disp = i->FastGetSolutionStepValue(DISPLACEMENT);
				(i)->X() = (i)->X0() + disp[0];
				(i)->Y() = (i)->Y0() + disp[1];
				(i)->Z() = (i)->Z0() + disp[2];
			}
			KRATOS_CATCH("")
		    }

	  }; // end class lineseaches
} // namespace Kratos
#endif



						      

		    



