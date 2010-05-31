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
 *   Last Modified by:    $Author: Kazem $
 *   Date:                $Date: 2008-07-25 14:48:17 $
 *   Revision:            $Revision: 1.1 $
 *
 * ***********************************************************/


#if !defined(KRATOS_TRILINOS_PREDICTOR_CORRECTOR_VELOCITY_BOSSAK_SCHEME )
#define  KRATOS_TRILINOS_PREDICTOR_CORRECTOR_VELOCITY_BOSSAK_SCHEME


/* System includes */


/* External includes */
#include "boost/smart_ptr.hpp"


/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "../../../incompressible_fluid_application/custom_strategies/strategies/residualbased_predictorcorrector_velocity_bossak_scheme.h"
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

    /** Short class definition.
	
      This class provides the implementation of the basic tasks that are needed by the solution strategy.
      It is intended to be the place for tailoring the solution strategies to problem specific tasks.
	  
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
    class TDenseSpace //= DenseSpace<double>
    >
    class TrilinosPredictorCorrectorVelocityBossakScheme
        : public ResidualBasedPredictorCorrectorVelocityBossakScheme<TSparseSpace, TDenseSpace>
    {
    public:
        /**@name Type Definitions */
        /*@{ */

        //typedef boost::shared_ptr< ResidualBasedPredictorCorrectorBossakScheme<TSparseSpace,TDenseSpace> > Pointer;

        KRATOS_CLASS_POINTER_DEFINITION(TrilinosPredictorCorrectorVelocityBossakScheme );

        typedef ResidualBasedPredictorCorrectorVelocityBossakScheme<TSparseSpace, TDenseSpace> BaseType;

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

        TrilinosPredictorCorrectorVelocityBossakScheme(double NewAlphaBossak, double MoveMeshStrategy)
        : ResidualBasedPredictorCorrectorVelocityBossakScheme<TSparseSpace,TDenseSpace>(NewAlphaBossak, MoveMeshStrategy)
        {

            std::cout << "using the TRILINOS velocity Bossak Time Integration Scheme" << std::endl;
        }

        /** Destructor.
         */
        virtual ~TrilinosPredictorCorrectorVelocityBossakScheme() {
        }


        /*@} */
        /**@name Operators
         */
        /*@{ */

        /**
                Performing the update of the solution.
         */

        //***************************************************************************

        virtual void Update(
                ModelPart& r_model_part,
                DofsArrayType& rDofSet,
                TSystemMatrixType& A,
                TSystemVectorType& Dv,
                TSystemVectorType& b
                ) {
            KRATOS_TRY

            BasicUpdateOperations(r_model_part, rDofSet, A, Dv, b);

            AdditionalUpdateOperations(r_model_part, rDofSet, A, Dv, b);

            KRATOS_CATCH("")
        }

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

#endif /* KRATOS_TRILINOS_PREDICTOR_CORRECTOR_VELOCITY_BOSSAK_SCHEME  defined */



