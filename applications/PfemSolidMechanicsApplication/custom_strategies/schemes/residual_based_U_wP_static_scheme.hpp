//
//   Project Name:        KratosPfemSolidMechanicsApplication $
//   Created by:          $Author:                  LMonforte $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                    July 2015 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined(KRATOS_RESIDUAL_BASED_UwP_STATIC_SCHEME )
#define  KRATOS_RESIDUAL_BASED_UwP_STATIC_SCHEME


/* System includes */


/* External includes */
#include "boost/smart_ptr.hpp"


/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/schemes/scheme.h"
#include "includes/variables.h"

#include "pfem_solid_mechanics_application_variables.h"
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"

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

   template<class TSparseSpace, class TDenseSpace > //= DenseSpace<double>
      class ResidualBasedUwPStaticScheme: public ResidualBasedIncrementalUpdateStaticScheme<TSparseSpace,TDenseSpace>
   {

      public:
         /**@name Type Definitions */
         /*@{ */

         //typedef boost::shared_ptr< ResidualBasedStaticScheme<TSparseSpace,TDenseSpace> > Pointer;
         KRATOS_CLASS_POINTER_DEFINITION( ResidualBasedUwPStaticScheme );

         typedef Scheme<TSparseSpace,TDenseSpace> BaseType;

         typedef typename BaseType::TDataType TDataType;

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
         ResidualBasedUwPStaticScheme()
            : ResidualBasedIncrementalUpdateStaticScheme<TSparseSpace,TDenseSpace>()
         {
         }

         /** Destructor.
          */
         virtual ~ResidualBasedUwPStaticScheme() {}


         /*@} */
         /**@name Operators
          */
         /*@{ */


         /**
           Performing the prediction of the solution.
          */
         virtual void Predict(
               ModelPart& r_model_part,
               DofsArrayType& rDofSet,
               TSystemMatrixType& A,
               TSystemVectorType& Dx,
               TSystemVectorType& b
               )
         {
            KRATOS_TRY


               for (ModelPart::NodeIterator i = r_model_part.NodesBegin();
                     i != r_model_part.NodesEnd(); ++i)
               {
                  //predicting displacement = PreviousDisplacement + PreviousVelocity * DeltaTime;
                  //ATTENTION::: the prediction is performed only on free nodes

                  array_1d<double, 3 > & PreviousDisplacement = (i)->FastGetSolutionStepValue(DISPLACEMENT, 1);
                  array_1d<double, 3 > & CurrentDisplacement  = (i)->FastGetSolutionStepValue(DISPLACEMENT);
                  //array_1d<double, 3 > & ImposedDisplacement  = (i)->FastGetSolutionStepValue(IMPOSED_DISPLACEMENT);
                  //array_1d<double, 3 > & Velocity             = (i)->FastGetSolutionStepValue(VELOCITY);
                  //double time_step = r_model_part.GetProcessInfo()[DELTA_TIME];

                  double               & CurrentWaterPressure   = (i)->FastGetSolutionStepValue(WATER_PRESSURE);
                  const double         & PreviousWaterPressure  = (i)->FastGetSolutionStepValue(WATER_PRESSURE, 1);
                  //const double         & ImposedWaterPressure   = (i)->FastGetSolutionStepValue(IMPOSED_WATER_PRESSURE);

                  if ((i->pGetDof(DISPLACEMENT_X))->IsFixed() == false)
                  {
                     CurrentDisplacement[0] = PreviousDisplacement[0];
                  }
                  else
                  {
                     //CurrentDisplacement[0]  = PreviousDisplacement[0] + ImposedDisplacement[0];//to impose fixed displacements;
                     //CurrentDisplacement[0] +=  time_step * Velocity[0];
                  }

                  if (i->pGetDof(DISPLACEMENT_Y)->IsFixed() == false)
                  {
                     CurrentDisplacement[1] = PreviousDisplacement[1]; 
                  }
                  else
                  {
                     //CurrentDisplacement[1]  = PreviousDisplacement[1] + ImposedDisplacement[1];//to impose fixed displacements;
                     //CurrentDisplacement[1] +=  time_step * Velocity[1];
                  }


                  if (i->HasDofFor(DISPLACEMENT_Z))
                  {
                     if (i->pGetDof(DISPLACEMENT_Z)->IsFixed() == false)
                     {
                        CurrentDisplacement[2] = PreviousDisplacement[2]; 
                     }
                     else
                     {
                        //CurrentDisplacement[2]= PreviousDisplacement[2] + ImposedDisplacement[2];//to impose fixed displacements;
                        //CurrentDisplacement[2] +=  time_step * Velocity[2];
                     }
                  }

                  if (i->HasDofFor(WATER_PRESSURE))
                  {
                     if (i->pGetDof(WATER_PRESSURE)->IsFixed() == false)
                     {
                        CurrentWaterPressure = PreviousWaterPressure ;
                     }
                     else
                     {
                        //CurrentWaterPressure =  ImposedWaterPressure ; // TO BE FIX, have to declare...
                        CurrentWaterPressure =  PreviousWaterPressure ;
                     }
                  }

					// ESTO QUE SIGUE ES UN APAÑO PARA PONER VELOCIDADES EN DIRICHLET
                  /*if (i->HasDofFor(PRESSURE))
                  {
                     double               & CurrentPressure   = (i)->FastGetSolutionStepValue(PRESSURE);
                     const double         & PreviousPressure  = (i)->FastGetSolutionStepValue(PRESSURE, 1);

                     if (i->pGetDof(PRESSURE)->IsFixed() == false)
                     {
                        CurrentPressure = PreviousPressure ;
                     }
                     else
                     {
                        CurrentPressure =  PreviousPressure ; // TO BE FIX, have to declare...
                     }
                  }*/
               }

            KRATOS_CATCH( "" )
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

}  /* namespace Kratos.*/

#endif /* KRATOS_RESIDUAL_BASED_STATIC_SCHEME  defined */

