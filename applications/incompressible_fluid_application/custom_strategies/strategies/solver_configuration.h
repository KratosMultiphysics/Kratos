#if !defined(KRATOS_RESIDUALBASED_SOLVER_CONFIGURATION )
#define  KRATOS_RESIDUALBASED_SOLVER_CONFIGURATION


/* System includes */


/* External includes */
#include "boost/smart_ptr.hpp"


/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/strategies/solving_strategy.h"


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
	class TDenseSpace,
	class TLinearSolver
	>
	class SolverConfiguration
	{
	public:
		/**@name Type Definitions */
		/*@{ */

		/** Counted pointer of ClassName */

		/*@} */
		/**@name Life Cycle
		*/
		/*@{ */

		/** Constructor.
		*/
                SolverConfiguration(ModelPart& model_part, unsigned int domain_size)
                        : mrModelPart(model_part), mDomainSize(domain_size)
                {
                }

		/** Destructor.
		*/

		/*@} */
		/**@name Operators
		*/
		/*@{ */

                unsigned int GetDomainSize()
                {return this->mDomainSize;}

                virtual typename SolvingStrategy<TSparseSpace,TDenseSpace,TLinearSolver>::Pointer pGetStrategy(const std::string& strategy_name )
                {
                    KRATOS_ERROR(std::logic_error,"accessing to the SolverConfiguration base class","");
                }

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
                ModelPart& mrModelPart;
                unsigned int mDomainSize;

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

                /*@{ */
		//this funcion is needed to ensure that all the memory is allocated correctly


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

		/** Copy constructor.
		*/


		/*@} */

	}; /* Class FractionalStepStrategy */

	/*@} */

	/**@name Type Definitions */
	/*@{ */


	/*@} */

}  /* namespace Kratos.*/

#endif /* KRATOS_RESIDUALBASED_SOLVER_CONFIGURATION  defined */
