//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi

#if !defined(KRATOS_RESIDUALBASED_PREDICTOR_CORRECTOR_VELOCITY_BDF_TURBULENT_SCHEME_NO_REACTION )
#define  KRATOS_RESIDUALBASED_PREDICTOR_CORRECTOR_VELOCITY_BDF_TURBULENT_SCHEME_NO_REACTION


/* System includes */


/* External includes */
#include "boost/smart_ptr.hpp"

/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
//#include "solving_strategies/schemes/scheme.h"
#include "residualbased_predictorcorrector_velocity_bdf_scheme_turbulent.h"
#include "includes/variables.h"
#include "includes/cfd_variables.h"
#include "containers/array_1d.h"
#include "utilities/openmp_utils.h"
#include "utilities/coordinate_transformation_utilities.h"


namespace Kratos {


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

	/// BDF2 time scheme for the incompressible flow problem.
	/** This scheme implements update operations and the calculation of the BDF coefficients for variable time step sizes.
	*
	* WARNING: this scheme assumes that the element internally implements the BDF2 scheme and is hence NOT compatible with the
	* elements ASGS2D, ASGS3D, VMS, MonolithicWallConditon
	*
	* the compatible element so far is
	*   @see TwoFluidVMS
	*
	* note also that in the prediction step only the velocity, and NOT the pressure is extrapolated in time.
	*/
	template<class TSparseSpace,
	class TDenseSpace //= DenseSpace<double>
	>
	class ResidualBasedPredictorCorrectorBDFSchemeTurbulentNoReaction : public ResidualBasedPredictorCorrectorBDFSchemeTurbulent<TSparseSpace, TDenseSpace> {
	public:
		/**@name Type Definitions */
		/*@{ */

		KRATOS_CLASS_POINTER_DEFINITION(ResidualBasedPredictorCorrectorBDFSchemeTurbulentNoReaction);

		typedef Scheme<TSparseSpace, TDenseSpace> BaseType;

		typedef typename BaseType::TDataType TDataType;

		typedef typename BaseType::DofsArrayType DofsArrayType;

		typedef typename Element::DofsVectorType DofsVectorType;

		typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

		typedef typename BaseType::TSystemVectorType TSystemVectorType;

		typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

		typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

		typedef Element::GeometryType  GeometryType;


		/*@} */
		/**@name Life Cycle
		*/
		/*@{ */

		/** Constructor without a turbulence model
		*/
		ResidualBasedPredictorCorrectorBDFSchemeTurbulentNoReaction(unsigned int DomainSize): ResidualBasedPredictorCorrectorBDFSchemeTurbulent<TSparseSpace, TDenseSpace>(DomainSize)
		{
		}

		/** Constructor without a turbulence model
		*/
		ResidualBasedPredictorCorrectorBDFSchemeTurbulentNoReaction(unsigned int DomainSize,Variable<double>& rSlipVar)
			:ResidualBasedPredictorCorrectorBDFSchemeTurbulent<TSparseSpace, TDenseSpace>(DomainSize, rSlipVar)
		{
		}

		/** Constructor with a turbulence model
		*/
		ResidualBasedPredictorCorrectorBDFSchemeTurbulentNoReaction(unsigned int DomainSize, Process::Pointer pTurbulenceModel)
			:ResidualBasedPredictorCorrectorBDFSchemeTurbulent<TSparseSpace, TDenseSpace>(DomainSize, pTurbulenceModel) // Second argument is number of matrix rows per node: monolithic elements have velocity and pressure dofs
		{
		}
		//************************************************************************************************
		//************************************************************************************************

		void ComputeReactions(ModelPart &rModelPart, TSystemMatrixType &A, TSystemVectorType &Dx, TSystemVectorType &b) override
		{}
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


		/*@} */
		/**@name Private Operators*/
		/*@{ */

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

	}; /* Class Scheme */

	   /*@} */

	   /**@name Type Definitions */
	   /*@{ */


	   /*@} */

} /* namespace Kratos.*/

#endif