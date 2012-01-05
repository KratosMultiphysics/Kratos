/*
==============================================================================
KratosTrilinosApplication
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

#if !defined(KRATOS_TRILINOS_PREDICTOR_CORRECTOR_VELOCITY_BOSSAK_SCHEME_TURBULENT )
#define  KRATOS_TRILINOS_PREDICTOR_CORRECTOR_VELOCITY_BOSSAK_SCHEME_TURBULENT


/* System includes */


/* External includes */
#include "boost/smart_ptr.hpp"


/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "processes/process.h"
#include "../../../FluidDynamicsApplication/custom_strategies/strategies/residualbased_predictorcorrector_velocity_bossak_scheme_turbulent.h"
#include "includes/variables.h"
#include "containers/array_1d.h"

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

    /// Trilinos version of ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulent.
    template<class TSparseSpace,
    class TDenseSpace //= DenseSpace<double>
    >
    class TrilinosPredictorCorrectorVelocityBossakSchemeTurbulent
        : public ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulent<TSparseSpace, TDenseSpace>
    {
    public:
        /**@name Type Definitions */
        /*@{ */

        //typedef boost::shared_ptr< ResidualBasedPredictorCorrectorBossakScheme<TSparseSpace,TDenseSpace> > Pointer;

        KRATOS_CLASS_POINTER_DEFINITION(TrilinosPredictorCorrectorVelocityBossakSchemeTurbulent );

        typedef ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulent<TSparseSpace, TDenseSpace> BaseType;

        typedef typename BaseType::TDataType TDataType;

        typedef typename BaseType::DofsArrayType DofsArrayType;

        typedef typename Element::DofsVectorType DofsVectorType;

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

        TrilinosPredictorCorrectorVelocityBossakSchemeTurbulent(double NewAlphaBossak, double MoveMeshStrategy, unsigned int DomainSize)
        : ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulent<TSparseSpace,TDenseSpace>(NewAlphaBossak, MoveMeshStrategy, DomainSize)
        {
            std::cout << "using the TRILINOS velocity Bossak Time Integration Scheme (with turbulence model)" << std::endl;
        }

        TrilinosPredictorCorrectorVelocityBossakSchemeTurbulent(double NewAlphaBossak, double MoveMeshStrategy, unsigned int DomainSize, Process::Pointer pTurbulenceModel)
        : ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulent<TSparseSpace,TDenseSpace>(NewAlphaBossak, MoveMeshStrategy, DomainSize, pTurbulenceModel)
        {

            std::cout << "using the TRILINOS velocity Bossak Time Integration Scheme (with turbulence model)" << std::endl;
        }

        /** Destructor.
         */
        virtual ~TrilinosPredictorCorrectorVelocityBossakSchemeTurbulent() {
        }


        /*@} */
        /**@name Operators
         */
        /*@{ */

        virtual int Check(ModelPart& rModelPart)
        {
            KRATOS_TRY

            int ErrorCode = BaseType::Check(rModelPart);
            if (ErrorCode != 0) return ErrorCode;

            // Check buffer size
            if (rModelPart.GetBufferSize() < 2)
                KRATOS_ERROR(std::logic_error, "GearScheme error: Insufficient buffer size for Bossak scheme, should be at least 2, got ",rModelPart.GetBufferSize());

            // Check that all required variables were registered
            if(DELTA_TIME.Key() == 0)
                KRATOS_ERROR(std::invalid_argument,"DELTA_TIME Key is 0. Check if all applications were correctly registered.","");
            if(OSS_SWITCH.Key() == 0)
                KRATOS_ERROR(std::invalid_argument,"OSS_SWITCH Key is 0. Check if all applications were correctly registered.","");

            if(DISPLACEMENT.Key() == 0)
                KRATOS_ERROR(std::invalid_argument,"DISPLACEMENT Key is 0. Check if all applications were correctly registered.","");
            if(VELOCITY.Key() == 0)
                KRATOS_ERROR(std::invalid_argument,"VELOCITY Key is 0. Check if all applications were correctly registered.","");
            if(MESH_VELOCITY.Key() == 0)
                KRATOS_ERROR(std::invalid_argument,"MESH_VELOCITY Key is 0. Check if all applications were correctly registered.","");
            if(ACCELERATION.Key() == 0)
                KRATOS_ERROR(std::invalid_argument,"ACCELERATION Key is 0. Check if all applications were correctly registered.","");

            // Checks on process info
            const ProcessInfo& rCurrentProcessInfo = rModelPart.GetProcessInfo();

            if(rCurrentProcessInfo.Has(DELTA_TIME) == 0)
                KRATOS_ERROR(std::invalid_argument,"ProcessInfo does not contain a value for DELTA_TIME","");

            return 0;
            KRATOS_CATCH("");
        }

        /**
                Performing the update of the solution.
         */

//        virtual void Update(
//                ModelPart& r_model_part,
//                DofsArrayType& rDofSet,
//                TSystemMatrixType& A,
//                TSystemVectorType& Dv,
//                TSystemVectorType& b
//                ) {
//            KRATOS_TRY

//            BasicUpdateOperations(r_model_part, rDofSet, A, Dv, b);

//            this->AdditionalUpdateOperations(r_model_part, rDofSet, A, Dv, b);

//            KRATOS_CATCH("")
//        }

        //***************************************************************************
		void BasicUpdateOperations(
			ModelPart& r_model_part,
			DofsArrayType& rDofSet,
			TSystemMatrixType& A,
			TSystemVectorType& Dx,
			TSystemVectorType& b
			)
		{
			KRATOS_TRY

			int system_size = TSparseSpace::Size(Dx);
			int number_of_dofs = rDofSet.size();
                        std::vector< int > index_array(number_of_dofs);

			//filling the array with the global ids
			int counter = 0;
			for(typename DofsArrayType::iterator i_dof = rDofSet.begin() ; i_dof != rDofSet.end() ; ++i_dof)
			{
				int id = i_dof->EquationId();
				if( id < system_size )
				{
					index_array[counter] = id;
					counter += 1;
				}
			}
			double tot_update_dofs = counter;

//                        double* pValues;
//                        TSparseSpace::GatherValues(Dx, index_array, pValues)

			//defining a map as needed
			Epetra_Map dof_update_map(-1,tot_update_dofs, &(*(index_array.begin())),0,b.Comm() );

			//defining the importer class
			Epetra_Import importer( dof_update_map, Dx.Map() );

			//defining a temporary vector to gather all of the values needed
			Epetra_Vector temp( importer.TargetMap() );

			//importing in the new temp vector the values
			int ierr = temp.Import(Dx,importer,Insert) ;
			if(ierr != 0) KRATOS_ERROR(std::logic_error,"Epetra failure found","");

                        double* temp_values;
			temp.ExtractView( &temp_values );

			b.Comm().Barrier();

			//performing the update
			typename DofsArrayType::iterator dof_begin = rDofSet.begin();
			for(unsigned int iii=0; iii<rDofSet.size(); iii++)
			{
				int global_id = (dof_begin+iii)->EquationId();
				if(global_id < system_size)
				{
					double aaa = temp[dof_update_map.LID(global_id)];
					(dof_begin+iii)->GetSolutionStepValue() += aaa;
				}
			}

//			//performing the update
//			typename DofsArrayType::iterator dof_begin = rDofSet.begin();
//			for(unsigned int iii=0; iii<rDofSet.size(); iii++)
//			{
//                            (dof_begin+iii)->GetSolutionStepValue() += pValues[iii];
//			}


			//removing unnecessary memory
//			delete [] temp_values; //DO NOT DELETE THIS!!
 
			KRATOS_CATCH("")
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
        /*		Matrix mMass;
                        Matrix mDamp;

                        Vector mvel;
                        Vector macc;
                        Vector maccold;

                        DofsVectorType mElementalDofList;
         */

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

#endif /* KRATOS_TRILINOS_PREDICTOR_CORRECTOR_VELOCITY_BOSSAK_SCHEME_TURBULENT  defined */



