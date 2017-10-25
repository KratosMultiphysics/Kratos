//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:                 PNavas $
//   Last modified by:    $Co-Author:           LMonforte $
//   Date:                $Date:             October 2017 $
//   Revision:            $Revision:                 -0.1 $
//
//

#if !defined(KRATOS_RESIDUAL_U_W_BASED_BOSSAK_SCHEME )
#define  KRATOS_RESIDUAL_U_W_BASED_BOSSAK_SCHEME

/* System includes */

/* External includes */
#include "boost/smart_ptr.hpp"

/* Project includes */
#include "custom_strategies/schemes/residual_based_bossak_scheme.hpp"

#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/schemes/scheme.h"
#include "includes/variables.h"
#include "containers/array_1d.h"
#include "includes/element.h"
#include "pfem_solid_mechanics_application_variables.h"

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
   /*@} */

   template<class TSparseSpace,  class TDenseSpace >
      class ResidualBasedUWBossakScheme: public ResidualBasedBossakScheme<TSparseSpace,TDenseSpace>
   {

      public:


         /**@name Type Definitions */

         /*@{ */
         KRATOS_CLASS_POINTER_DEFINITION( ResidualBasedUWBossakScheme );

         typedef Scheme<TSparseSpace,TDenseSpace>                      BaseType;

         typedef typename BaseType::TDataType                         TDataType;

         typedef typename BaseType::DofsArrayType                 DofsArrayType;

         typedef typename Element::DofsVectorType                DofsVectorType;

         typedef typename BaseType::TSystemMatrixType         TSystemMatrixType;

         typedef typename BaseType::TSystemVectorType         TSystemVectorType;

         typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

         typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

         typedef ModelPart::ElementsContainerType             ElementsArrayType;

         typedef ModelPart::ConditionsContainerType         ConditionsArrayType;

         typedef typename BaseType::Pointer                     BaseTypePointer;

         /*@} */

         /**
          * Constructor.
          * The bossak method
          */
         ResidualBasedUWBossakScheme(double rAlpham=0,double rDynamic=1)
            :ResidualBasedBossakScheme<TSparseSpace,TDenseSpace>(rAlpham, rDynamic)
         {
         }


         /** Copy Constructor.
          */
         ResidualBasedUWBossakScheme(ResidualBasedUWBossakScheme& rOther)
            :ResidualBasedBossakScheme<TSparseSpace,TDenseSpace>(rOther)
         {
         }


         /** Destructor.
          */
         virtual ~ResidualBasedUWBossakScheme
            () {}

         /*@} */
         /**@name Operators
          */
         /*@{ */


         /**
          * Clone 
          */
         virtual BaseTypePointer Clone()
         {
            return BaseTypePointer( new ResidualBasedUWBossakScheme(*this) );
         }



         //***************************************************************************
         //***************************************************************************

         /**
          * Performing the update of the solution
          * Incremental update within newton iteration. It updates the state variables at the end of the time step: u_{n+1}^{k+1}= u_{n+1}^{k}+ \Delta u
          * @param r_model_part
          * @param rDofSet set of all primary variables
          * @param A	LHS matrix
          * @param Dx incremental update of primary variables
          * @param b RHS Vector
          */
         void Update(
               ModelPart& r_model_part,
               DofsArrayType& rDofSet,
               TSystemMatrixType& A,
               TSystemVectorType& Dx,
               TSystemVectorType& b )
         {
            KRATOS_TRY

            //std::cout << " Update " << std::endl;
            //update of displacement (by DOF)
            for (typename DofsArrayType::iterator i_dof = rDofSet.begin(); i_dof != rDofSet.end(); ++i_dof)
            {
               if (i_dof->IsFree() )
               {
                  i_dof->GetSolutionStepValue() += Dx[i_dof->EquationId()];
               }
             }

            //updating time derivatives (nodally for efficiency)
            array_1d<double, 3 > DeltaDisplacement;
            array_1d<double, 3 > DeltaWaterDisplacement;
            for (ModelPart::NodeIterator i = r_model_part.NodesBegin();
                  i != r_model_part.NodesEnd(); ++i)
            {

               noalias(DeltaDisplacement) = (i)->FastGetSolutionStepValue(DISPLACEMENT) - (i)->FastGetSolutionStepValue(DISPLACEMENT, 1);


               array_1d<double, 3 > & CurrentVelocity      = (i)->FastGetSolutionStepValue(VELOCITY, 0);
               array_1d<double, 3 > & PreviousVelocity     = (i)->FastGetSolutionStepValue(VELOCITY, 1);

               array_1d<double, 3 > & CurrentAcceleration  = (i)->FastGetSolutionStepValue(ACCELERATION, 0);
               array_1d<double, 3 > & PreviousAcceleration = (i)->FastGetSolutionStepValue(ACCELERATION, 1);

               this->UpdateVelocity     (CurrentVelocity, DeltaDisplacement, PreviousVelocity, PreviousAcceleration);

               this->UpdateAcceleration (CurrentAcceleration, DeltaDisplacement, PreviousVelocity, PreviousAcceleration);

               noalias(DeltaWaterDisplacement) = (i)->FastGetSolutionStepValue(WATER_DISPLACEMENT) - (i)->FastGetSolutionStepValue(WATER_DISPLACEMENT, 1);

               array_1d<double, 3 > & CurrentWaterVelocity      = (i)->FastGetSolutionStepValue(WATER_VELOCITY, 0);
               array_1d<double, 3 > & PreviousWaterVelocity     = (i)->FastGetSolutionStepValue(WATER_VELOCITY, 1);

               array_1d<double, 3 > & CurrentWaterAcceleration  = (i)->FastGetSolutionStepValue(WATER_ACCELERATION, 0);
               array_1d<double, 3 > & PreviousWaterAcceleration = (i)->FastGetSolutionStepValue(WATER_ACCELERATION, 1);

               this->UpdateVelocity     (CurrentWaterVelocity, DeltaWaterDisplacement, PreviousWaterVelocity, PreviousWaterAcceleration);

               this->UpdateAcceleration (CurrentWaterAcceleration, DeltaWaterDisplacement, PreviousWaterVelocity, PreviousWaterAcceleration);
            }

            KRATOS_CATCH( "" )
         }


         //***************************************************************************
         //***************************************************************************

         //predicts the solution for the current step:
         // x = xold + vold * Dt

         void Predict(
               ModelPart& r_model_part,
               DofsArrayType& rDofSet,
               TSystemMatrixType& A,
               TSystemVectorType& Dx,
               TSystemVectorType& b
               )
         {

            KRATOS_TRY

            //std::cout << " Prediction " << std::endl;
            array_1d<double, 3 > DeltaDisplacement;


            for (ModelPart::NodeIterator i = r_model_part.NodesBegin();
                  i != r_model_part.NodesEnd(); ++i)
            {
               //predicting displacement = PreviousDisplacement + PreviousVelocity * DeltaTime;
               //ATTENTION::: the prediction is performed only on free nodes

               array_1d<double, 3 > & PreviousVelocity     = (i)->FastGetSolutionStepValue(VELOCITY, 1);
               array_1d<double, 3 > & PreviousDisplacement = (i)->FastGetSolutionStepValue(DISPLACEMENT, 1);
               array_1d<double, 3 > & CurrentDisplacement  = (i)->FastGetSolutionStepValue(DISPLACEMENT);
               //array_1d<double, 3 > & ImposedDisplacement  = (i)->FastGetSolutionStepValue(IMPOSED_DISPLACEMENT);


               if ((i->pGetDof(DISPLACEMENT_X))->IsFixed() == false)
               {
                  CurrentDisplacement[0] = PreviousDisplacement[0];// + DeltaTime * PreviousVelocity[0];
               }
               else
               {
                  //CurrentDisplacement[0]  = PreviousDisplacement[0]; // LMV + ImposedDisplacement[0];//to impose fixed displacements;
                  //PreviousDisplacement[0] = 0;
               }

               if (i->pGetDof(DISPLACEMENT_Y)->IsFixed() == false)
               {
                  CurrentDisplacement[1] = PreviousDisplacement[1]; //+ DeltaTime * PreviousVelocity[1];
               }
               else
               {
                  //CurrentDisplacement[1]  = 0.01;
                     //PreviousDisplacement[1]; // LMV + ImposedDisplacement[1];//to impose fixed displacements;
                  //PreviousDisplacement[1] = 0;
               }


               if (i->HasDofFor(DISPLACEMENT_Z))
               {
                  if (i->pGetDof(DISPLACEMENT_Z)->IsFixed() == false)
                  {
                     CurrentDisplacement[2] = PreviousDisplacement[2]; // + DeltaTime * PreviousVelocity[2];
                  }
                  else
                  {
                     //CurrentDisplacement[2]  = PreviousDisplacement[2]; // LMV + ImposedDisplacement[2];//to impose fixed displacements;
                     //PreviousDisplacement[2] = 0;
                  }
               }


               // std::cout<<" DispPre "<<PreviousDisplacement<<" ID "<<i->Id()<<std::endl;
               // std::cout<<" DispCur "<<CurrentDisplacement<<" ID "<<i->Id()<<std::endl;

               if (i->HasDofFor(PRESSURE))
               {
                  double& PreviousPressure    = (i)->FastGetSolutionStepValue(PRESSURE, 1);
                  double& CurrentPressure     = (i)->FastGetSolutionStepValue(PRESSURE);

                  if ((i->pGetDof(PRESSURE))->IsFixed() == false)
                     CurrentPressure = PreviousPressure;

                  //std::cout<<" PressureCur [1] "<<CurrentPressure<<" PressurePre [1] "<<PreviousPressure<<" ID "<<i->Id()<<std::endl;
               }



               //updating time derivatives ::: please note that displacements and its time derivatives
               //can not be consistently fixed separately

               noalias(DeltaDisplacement) = CurrentDisplacement - PreviousDisplacement;

               array_1d<double, 3 > & PreviousAcceleration  = (i)->FastGetSolutionStepValue(ACCELERATION, 1);
               array_1d<double, 3 > & CurrentVelocity       = (i)->FastGetSolutionStepValue(VELOCITY);
               array_1d<double, 3 > & CurrentAcceleration   = (i)->FastGetSolutionStepValue(ACCELERATION);

               this->UpdateVelocity     (CurrentVelocity, DeltaDisplacement, PreviousVelocity, PreviousAcceleration);

               this->UpdateAcceleration (CurrentAcceleration, DeltaDisplacement, PreviousVelocity, PreviousAcceleration);

            }
            
            // the same but now with the water relative displacement to the soil
            array_1d<double, 3 > DeltaWaterDisplacement;


            for (ModelPart::NodeIterator i = r_model_part.NodesBegin();
                  i != r_model_part.NodesEnd(); ++i)
            {
               //predicting displacement = PreviousDisplacement + PreviousVelocity * DeltaTime;
               //ATTENTION::: the prediction is performed only on free nodes

               array_1d<double, 3 > & PreviousWaterVelocity     = (i)->FastGetSolutionStepValue(WATER_VELOCITY, 1);
               array_1d<double, 3 > & PreviousWaterDisplacement = (i)->FastGetSolutionStepValue(WATER_DISPLACEMENT, 1);
               array_1d<double, 3 > & CurrentWaterDisplacement  = (i)->FastGetSolutionStepValue(WATER_DISPLACEMENT);
               //array_1d<double, 3 > & ImposedWaterDisplacement  = (i)->FastGetSolutionStepValue(IMPOSED_WATER_DISPLACEMENT);


               if ((i->pGetDof(WATER_DISPLACEMENT_X))->IsFixed() == false)
               {
                  CurrentWaterDisplacement[0] = PreviousWaterDisplacement[0];// + DeltaTime * PreviousVelocity[0];
               }
               else
               {
                  //CurrentWaterDisplacement[0]  = PreviousWaterDisplacement[0]; //LMV + ImposedWaterDisplacement[0];//to impose fixed displacements;
                  //PreviousDisplacement[0] = 0;
               }

               if (i->pGetDof(WATER_DISPLACEMENT_Y)->IsFixed() == false)
               {
                  CurrentWaterDisplacement[1] = PreviousWaterDisplacement[1]; //+ DeltaTime * PreviousVelocity[1];
               }
               else
               {
                  //CurrentWaterDisplacement[1]  = PreviousWaterDisplacement[1]; //LMV + ImposedWaterDisplacement[1];//to impose fixed displacements;
                  //PreviousDisplacement[1] = 0;
               }


               if (i->HasDofFor(WATER_DISPLACEMENT_Z))
               {
                  if (i->pGetDof(WATER_DISPLACEMENT_Z)->IsFixed() == false)
                  {
                     CurrentWaterDisplacement[2] = PreviousWaterDisplacement[2]; // + DeltaTime * PreviousVelocity[2];
                  }
                  else
                  {
                     //CurrentWaterDisplacement[2]  = PreviousWaterDisplacement[2]; //LMV + ImposedWaterDisplacement[2];//to impose fixed displacements;
                     //PreviousDisplacement[2] = 0;
                  }
               }


               //updating time derivatives ::: please note that displacements and its time derivatives
               //can not be consistently fixed separately

               noalias(DeltaWaterDisplacement) = CurrentWaterDisplacement - PreviousWaterDisplacement;

               array_1d<double, 3 > & PreviousWaterAcceleration  = (i)->FastGetSolutionStepValue(WATER_ACCELERATION, 1);
               array_1d<double, 3 > & CurrentWaterVelocity       = (i)->FastGetSolutionStepValue(WATER_VELOCITY);
               array_1d<double, 3 > & CurrentWaterAcceleration   = (i)->FastGetSolutionStepValue(WATER_ACCELERATION);

               this->UpdateVelocity     (CurrentWaterVelocity, DeltaWaterDisplacement, PreviousWaterVelocity, PreviousWaterAcceleration);

               this->UpdateAcceleration (CurrentWaterAcceleration, DeltaWaterDisplacement, PreviousWaterVelocity, PreviousWaterAcceleration);

            }

            KRATOS_CATCH("")
         }


         //***************************************************************************
         //***************************************************************************

         /**
          * This function is designed to be called once to perform all the checks needed
          * on the input provided. Checks can be "expensive" as the function is designed
          * to catch user's errors.
          * @param r_model_part
          * @return 0 all ok
          */
         virtual int Check(ModelPart& r_model_part)
         {
            KRATOS_TRY

            return 0;
/*      int err = Scheme<TSparseSpace, TDenseSpace>::Check(r_model_part);
            if(err!=0) return err;

            //check for variables keys
            //verify that the variables are correctly initialized
            if(DISPLACEMENT.Key() == 0)
               KRATOS_THROW_ERROR( std::invalid_argument,"DISPLACEMENT has Key zero! (check if the application is correctly registered", "" )
                  if(VELOCITY.Key() == 0)
                     KRATOS_THROW_ERROR( std::invalid_argument,"VELOCITY has Key zero! (check if the application is correctly registered", "" )
                        if(ACCELERATION.Key() == 0)
                           KRATOS_THROW_ERROR( std::invalid_argument,"ACCELERATION has Key zero! (check if the application is correctly registered", "" )

                              //check that variables are correctly allocated
                              for(ModelPart::NodesContainerType::iterator it=r_model_part.NodesBegin();
                                    it!=r_model_part.NodesEnd(); it++)
                              {
                                 if (it->SolutionStepsDataHas(DISPLACEMENT) == false)
                                    KRATOS_THROW_ERROR( std::logic_error, "DISPLACEMENT variable is not allocated for node ", it->Id() )
                                       if (it->SolutionStepsDataHas(VELOCITY) == false)
                                          KRATOS_THROW_ERROR( std::logic_error, "DISPLACEMENT variable is not allocated for node ", it->Id() )
                                             if (it->SolutionStepsDataHas(ACCELERATION) == false)
                                                KRATOS_THROW_ERROR( std::logic_error, "DISPLACEMENT variable is not allocated for node ", it->Id() )
                              }

            //check that dofs exist
            for(ModelPart::NodesContainerType::iterator it=r_model_part.NodesBegin();
                  it!=r_model_part.NodesEnd(); it++)
            {
               if(it->HasDofFor(DISPLACEMENT_X) == false)
                  KRATOS_THROW_ERROR( std::invalid_argument,"missing DISPLACEMENT_X dof on node ",it->Id() )
                     if(it->HasDofFor(DISPLACEMENT_Y) == false)
                        KRATOS_THROW_ERROR( std::invalid_argument,"missing DISPLACEMENT_Y dof on node ",it->Id() )
                           if(it->HasDofFor(DISPLACEMENT_Z) == false)
                              KRATOS_THROW_ERROR( std::invalid_argument,"missing DISPLACEMENT_Z dof on node ",it->Id() )
            }


            //check for admissible value of the AlphaBossak
            if(mAlpha.m > 0.0 || mAlpha.m < -0.3)
               KRATOS_THROW_ERROR( std::logic_error,"Value not admissible for AlphaBossak. Admissible values should be between 0.0 and -0.3. Current value is ", mAlpha.m )

                  //check for minimum value of the buffer index
                  //verify buffer size
                  if (r_model_part.GetBufferSize() < 2)
                     KRATOS_THROW_ERROR( std::logic_error, "insufficient buffer size. Buffer size should be greater than 2. Current size is", r_model_part.GetBufferSize() )


                        return 0; */
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

   protected:
         /**@name Static Member Variables */
         /*@{ */
         /*@} */
         /**@name Member Variables */
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
         /**@name Unaccessible methods */
         /*@{ */
}; /* Class ResidualBasedUWBossakScheme */
}  /* namespace Kratos.*/

#endif /* KRATOS_RESIDUAL_U_W_BASED_BOSSAK_SCHEME defined */


