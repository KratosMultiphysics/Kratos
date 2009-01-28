/* *********************************************************   
*          
*   Last Modified by:    $Author: pooyan $
*   Date:                $Date: 2006-11-27 16:07:50 $
*   Revision:            $Revision: 1.1.1.1 $
*
* ***********************************************************/


#if !defined(KRATOS_NEW_LAPLACIAN_MESHMOVING_STRATEGY )
#define  KRATOS_NEW_LAPLACIAN_MESHMOVING_STRATEGY


/* System includes */


/* External includes */
#include "boost/smart_ptr.hpp"


/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/strategies/solving_strategy.h"
#include "solving_strategies/strategies/residualbased_linear_strategy.h"
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"
#include "custom_elements/laplacian_meshmoving_element_2d.h"
#include "custom_elements/laplacian_meshmoving_element_3d.h"
#include "ale_application.h"
#include "solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_componentwise.h"


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

	/// Short class definition.
	/**   Detail class definition.

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
	class TDenseSpace, //= DenseSpace<double>,
	class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
	>
	class LaplacianMeshMovingStrategy 
		: public SolvingStrategy<TSparseSpace,TDenseSpace,TLinearSolver>
	{
	public:
		/**@name Type Definitions */       
		/*@{ */

		/** Counted pointer of ClassName */
		//typedef boost::shared_ptr< LaplacianMeshMovingStrategy<TSparseSpace,TDenseSpace,TLinearSolver> > Pointer;
		KRATOS_CLASS_POINTER_DEFINITION( LaplacianMeshMovingStrategy );

		typedef SolvingStrategy<TSparseSpace,TDenseSpace,TLinearSolver> BaseType;

		typedef typename BaseType::TDataType TDataType;

		//typedef typename BaseType::DofSetType DofSetType;

		typedef typename BaseType::DofsArrayType DofsArrayType;

		typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

		typedef typename BaseType::TSystemVectorType TSystemVectorType;

		typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

		typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;



		/*@} */
		/**@name Life Cycle 
		*/    
		/*@{ */

		/** Constructor.
		*/
		LaplacianMeshMovingStrategy(
			ModelPart& model_part, 
			typename TLinearSolver::Pointer pNewLinearSolver,
			int dimension = 3,
			int velocity_order = 1,
			bool reform_dof_at_every_step = true
			)
			: SolvingStrategy<TSparseSpace,TDenseSpace,TLinearSolver>(model_part)
		{
			KRATOS_TRY

				//Generating Mesh Part
				GenerateMeshPart(dimension);

			mdimension = dimension;
			mvel_order = velocity_order;
			mreform_dof_at_every_step = reform_dof_at_every_step;

			typedef Scheme< TSparseSpace,  TDenseSpace > SchemeType;
			typename SchemeType::Pointer pscheme = typename SchemeType::Pointer
				( new ResidualBasedIncrementalUpdateStaticScheme< TSparseSpace,  TDenseSpace >() );

			bool CalculateReactions = false;
			bool ReformDofAtEachIteration = false;
			bool CalculateNormDxFlag = false;
			
 
			//computation of the fractional vel velocity (first step)
//			mstrategy = typename BaseType::Pointer(
//				new ResidualBasedLinearStrategy<TSparseSpace,  TDenseSpace, TLinearSolver >
//				(*mpMeshModelPart,pscheme,pNewLinearSolver,CalculateReactions,ReformDofAtEachIteration,CalculateNormDxFlag)  );
//			mstrategy->SetEchoLevel(2);

			typedef typename BuilderAndSolver<TSparseSpace,TDenseSpace,TLinearSolver>::Pointer BuilderSolverTypePointer;
//			typedef typename Variable< double > VarComponent;

			BuilderSolverTypePointer aux_var_build = BuilderSolverTypePointer(new	ResidualBasedEliminationBuilderAndSolverComponentwise<TSparseSpace,TDenseSpace,TLinearSolver, Variable<double> >(pNewLinearSolver,AUX_MESH_VAR) );
			
			mstrategy = typename BaseType::Pointer( new ResidualBasedLinearStrategy<TSparseSpace,  TDenseSpace, TLinearSolver >				(*mpMeshModelPart,pscheme,pNewLinearSolver,aux_var_build,CalculateReactions,ReformDofAtEachIteration,CalculateNormDxFlag)  );
			mstrategy->SetEchoLevel(2);

			for(ModelPart::NodeIterator i = (*mpMeshModelPart).NodesBegin() ; 
				i != (*mpMeshModelPart).NodesEnd() ; ++i)
			{
				if(!i->IsFixed(DISPLACEMENT_X))
					(i)->GetSolutionStepValue(DISPLACEMENT_X) = 0.00;
				if(!i->IsFixed(DISPLACEMENT_Y))
					(i)->GetSolutionStepValue(DISPLACEMENT_Y) = 0.00;
				if(!i->IsFixed(DISPLACEMENT_Z))
					(i)->GetSolutionStepValue(DISPLACEMENT_Z) = 0.00;
			}

			KRATOS_CATCH("")
		}



		/** Destructor.
		*/
		virtual ~LaplacianMeshMovingStrategy() {}

		/** Destructor.
		*/


		//*********************************************************************************
		//*********************************************************************************
		double Solve()
		{
			KRATOS_TRY

			ReGenerateMeshPart();

			ProcessInfo& rCurrentProcessInfo = (mpMeshModelPart)->GetProcessInfo();

			//updating the time
			rCurrentProcessInfo[TIME] = BaseType::GetModelPart().GetProcessInfo()[TIME];
			rCurrentProcessInfo[DELTA_TIME] = BaseType::GetModelPart().GetProcessInfo()[DELTA_TIME];

			//Fixing Dofs As Needed
			for(ModelPart::NodeIterator i = (*mpMeshModelPart).NodesBegin() ; 
				i != (*mpMeshModelPart).NodesEnd() ; ++i)
			{
				if(i->IsFixed(DISPLACEMENT_X))
					(i)->Fix(AUX_MESH_VAR);
				if(i->IsFixed(DISPLACEMENT_Y))
					(i)->Fix(AUX_MESH_VAR);
				if(i->IsFixed(DISPLACEMENT_Z))
					(i)->Fix(AUX_MESH_VAR);
			}
			//X DIRECTION
			rCurrentProcessInfo[FRACTIONAL_STEP] = 1; //laplacian mesh moving type corresponds to -1
			for(ModelPart::NodeIterator i = (*mpMeshModelPart).NodesBegin() ; 
				i != (*mpMeshModelPart).NodesEnd() ; ++i)
			{
				(i)->FastGetSolutionStepValue(AUX_MESH_VAR) = (i)->GetSolutionStepValue(DISPLACEMENT_X);
			}
			mstrategy->Solve();
			for(ModelPart::NodeIterator i = (*mpMeshModelPart).NodesBegin() ; 
				i != (*mpMeshModelPart).NodesEnd() ; ++i)
			{
				(i)->GetSolutionStepValue(DISPLACEMENT_X) = (i)->FastGetSolutionStepValue(AUX_MESH_VAR); 
				//KRATOS_WATCH((i)->FastGetSolutionStepValue(DISPLACEMENT_X));
			}

			//Y DIRECTION
			if(mdimension > 1)
			{
				rCurrentProcessInfo[FRACTIONAL_STEP] = 2; //laplacian mesh moving type corresponds to -1
				for(ModelPart::NodeIterator i = (*mpMeshModelPart).NodesBegin() ; 
					i != (*mpMeshModelPart).NodesEnd() ; ++i)
				{
					(i)->FastGetSolutionStepValue(AUX_MESH_VAR) = (i)->GetSolutionStepValue(DISPLACEMENT_Y);
				}
				mstrategy->Solve();
				for(ModelPart::NodeIterator i = (*mpMeshModelPart).NodesBegin() ; 
					i != (*mpMeshModelPart).NodesEnd() ; ++i)
				{
					(i)->GetSolutionStepValue(DISPLACEMENT_Y) = (i)->FastGetSolutionStepValue(AUX_MESH_VAR);
				}
			}

			//Z DIRECTION
			if(mdimension > 2)
			{
				rCurrentProcessInfo[FRACTIONAL_STEP] = 3; //laplacian mesh moving type corresponds to -1
				for(ModelPart::NodeIterator i = (*mpMeshModelPart).NodesBegin() ; 
					i != (*mpMeshModelPart).NodesEnd() ; ++i)
				{
					(i)->FastGetSolutionStepValue(AUX_MESH_VAR) = (i)->GetSolutionStepValue(DISPLACEMENT_Z);
				}
				mstrategy->Solve();
				for(ModelPart::NodeIterator i = (*mpMeshModelPart).NodesBegin() ; 
					i != (*mpMeshModelPart).NodesEnd() ; ++i)
				{
					(i)->GetSolutionStepValue(DISPLACEMENT_Z) = (i)->FastGetSolutionStepValue(AUX_MESH_VAR); 
				}
			}

			CalculateMeshVelocities();

			BaseType::MoveMesh();

			//clearing the system if needed
			if(mreform_dof_at_every_step == true)
				mstrategy->Clear();

			return 0.0;

			KRATOS_CATCH("")
		}


		//*********************************************************************************
		//*********************************************************************************
		void CalculateMeshVelocities()
		{
			KRATOS_TRY;

			double DeltaTime = (*mpMeshModelPart).GetProcessInfo()[DELTA_TIME];
			double coeff = 1/DeltaTime;
			if( mvel_order == 1) //mesh velocity calculated as (x(n+1)-x(n))/Dt
			{
				for(ModelPart::NodeIterator i = (*mpMeshModelPart).NodesBegin() ; 
					i != (*mpMeshModelPart).NodesEnd() ; ++i)
				{
					array_1d<double,3>& mesh_v = (i)->FastGetSolutionStepValue(MESH_VELOCITY);
					array_1d<double,3>& disp = (i)->FastGetSolutionStepValue(DISPLACEMENT);
					array_1d<double,3>& dispold = (i)->FastGetSolutionStepValue(DISPLACEMENT,1);
					noalias(mesh_v) =  disp - dispold;
					mesh_v *= coeff;
				}
			}
			else //mesh velocity calculated as (3*x(n+1)-4*x(n)+x(n-1))/(2*Dt)
			{
				double c1 = 1.50*coeff;
				double c2 = -2.0*coeff;
				double c3 = 0.50*coeff;

				for(ModelPart::NodeIterator i = (*mpMeshModelPart).NodesBegin() ; 
					i != (*mpMeshModelPart).NodesEnd() ; ++i)
				{
					array_1d<double,3>& mesh_v = (i)->FastGetSolutionStepValue(MESH_VELOCITY);
					noalias(mesh_v) =  c1 * (i)->FastGetSolutionStepValue(DISPLACEMENT);
					noalias(mesh_v) += c2 * (i)->FastGetSolutionStepValue(DISPLACEMENT,1);
					noalias(mesh_v) += c3 * (i)->FastGetSolutionStepValue(DISPLACEMENT,2);
				}
			}

			KRATOS_CATCH("")
		}

		virtual void SetEchoLevel(int Level) 
		{
			mstrategy->SetEchoLevel(Level);
		}

		void MoveNodes()
		{
			CalculateMeshVelocities();
			BaseType::MoveMesh();
		}

		/*@} */
		/**@name Operators 
		*/  
		/*@{ */

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
		ModelPart::Pointer mpMeshModelPart;

		typename BaseType::Pointer mstrategy;

		int mdimension;
		int mvel_order;
		bool mreform_dof_at_every_step;

		/*@} */
		/**@name Private Operators*/
		/*@{ */

		/*@} */
		/**@name Private Operations*/
		/*@{ */

		void GenerateMeshPart(int dimension)
		{
			mpMeshModelPart = ModelPart::Pointer( new ModelPart("MeshPart",1) );
			//initializing mesh nodes
			mpMeshModelPart->Nodes() = BaseType::GetModelPart().Nodes();
			//mpMeshModelPart->Nodes() = BaseType::GetModelPart().Nodes();

			//creating mesh elements
			ModelPart::ElementsContainerType& MeshElems = mpMeshModelPart->Elements();
			Element::Pointer pElem;

			if(dimension == 2)
				for(ModelPart::ElementsContainerType::iterator it =  BaseType::GetModelPart().ElementsBegin(); 
												it != BaseType::GetModelPart().ElementsEnd(); it++)
				{
					pElem = Element::Pointer(new LaplacianMeshMovingElem2D(
											(*it).Id(), 
											(*it).pGetGeometry(), 
											(*it).pGetProperties() ) );
					MeshElems.push_back(pElem);
				}
			if(dimension == 3)
				for(ModelPart::ElementsContainerType::iterator it =  BaseType::GetModelPart().ElementsBegin(); 
												it != BaseType::GetModelPart().ElementsEnd(); it++)
				{
					pElem = Element::Pointer(new LaplacianMeshMovingElem3D(
											(*it).Id(), 
											(*it).pGetGeometry(), 
											(*it).pGetProperties() ) );
					MeshElems.push_back(pElem);
				}
//KRATOS_WATCH((mpMeshModelPart->Elements()).size());
//KRATOS_WATCH((mpMeshModelPart->Nodes()).size());
//KRATOS_WATCH(mpMeshModelPart->GetMesh());
		}

		void ReGenerateMeshPart()
		{
			std::cout << "regenerating elements for the mesh motion scheme" << std::endl;

			//initializing mesh nodes
			mpMeshModelPart->Nodes().clear();
			mpMeshModelPart->Nodes() = BaseType::GetModelPart().Nodes();

			//creating mesh elements
			ModelPart::ElementsContainerType& MeshElems = mpMeshModelPart->Elements();
			Element::Pointer pElem;

			MeshElems.clear();
			MeshElems.reserve( MeshElems.size() );

			if(mdimension == 2)
				for(ModelPart::ElementsContainerType::iterator it =  BaseType::GetModelPart().ElementsBegin(); 
												it != BaseType::GetModelPart().ElementsEnd(); it++)
				{
					pElem = Element::Pointer(new LaplacianMeshMovingElem2D(
											(*it).Id(), 
											(*it).pGetGeometry(), 
											(*it).pGetProperties() ) );
					MeshElems.push_back(pElem);
				}
			else
				for(ModelPart::ElementsContainerType::iterator it =  BaseType::GetModelPart().ElementsBegin(); 
												it != BaseType::GetModelPart().ElementsEnd(); it++)
				{
					pElem = Element::Pointer(new LaplacianMeshMovingElem3D(
											(*it).Id(), 
											(*it).pGetGeometry(), 
											(*it).pGetProperties() ) );
					MeshElems.push_back(pElem);
				}
				
			KRATOS_WATCH(MeshElems.size());
		}

		/*@} */
		/**@name Private  Access */
		/*@{ */


		/*@} */     
		/**@name Private Inquiry */
		/*@{ */


		/*@} */   
		/**@name Un accessible methods */
		/*@{ */

		/** Copy constructor.
		*/
		LaplacianMeshMovingStrategy(const LaplacianMeshMovingStrategy& Other);


		/*@} */   

	}; /* Class LaplacianMeshMovingStrategy */

	/*@} */

	/**@name Type Definitions */       
	/*@{ */


	/*@} */

}  /* namespace Kratos.*/

#endif /* KRATOS_NEW_LAPLACIAN_MESHMOVING_STRATEGY  defined */

