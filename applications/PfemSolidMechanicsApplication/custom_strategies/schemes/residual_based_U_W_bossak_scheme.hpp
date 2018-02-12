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

               if ( i->HasDofFor(WATER_PRESSURE) ) {
                  const double& PreviousWaterPressure    = (i)->FastGetSolutionStepValue(WATER_PRESSURE, 1);
                  const double& PreviousWaterPressureVelocity    = (i)->FastGetSolutionStepValue(WATER_PRESSURE_VELOCITY, 1);
                  const double& PreviousWaterPressureAcceleration    = (i)->FastGetSolutionStepValue(WATER_PRESSURE_ACCELERATION, 1);
                  double& CurrentWaterPressure     = (i)->FastGetSolutionStepValue(WATER_PRESSURE);
                  double& CurrentWaterPressureVelocity     = (i)->FastGetSolutionStepValue(WATER_PRESSURE_VELOCITY);
                  double& CurrentWaterPressureAcceleration     = (i)->FastGetSolutionStepValue(WATER_PRESSURE_ACCELERATION);

                  double DeltaWaterPressure = CurrentWaterPressure - PreviousWaterPressure;
                  UpdateVelocityScalar     ( CurrentWaterPressureVelocity, DeltaWaterPressure, PreviousWaterPressureVelocity, PreviousWaterPressureAcceleration);
                  UpdateAccelerationScalar ( CurrentWaterPressureAcceleration, DeltaWaterPressure, PreviousWaterPressureVelocity, PreviousWaterPressureAcceleration);
               }
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
      const double DeltaTime = r_model_part.GetProcessInfo()[DELTA_TIME];

            array_1d<double, 3 > DeltaDisplacement;


            for (ModelPart::NodeIterator i = r_model_part.NodesBegin();
                  i != r_model_part.NodesEnd(); ++i)
            {
               //predicting displacement = PreviousDisplacement + PreviousVelocity * DeltaTime;
               //ATTENTION::: the prediction is performed only on free nodes

               const array_1d<double, 3 > & PreviousAcceleration = (i)->FastGetSolutionStepValue(ACCELERATION, 1);
               const array_1d<double, 3 > & PreviousVelocity     = (i)->FastGetSolutionStepValue(VELOCITY,     1);
               const array_1d<double, 3 > & PreviousDisplacement = (i)->FastGetSolutionStepValue(DISPLACEMENT, 1);
               array_1d<double, 3 > & CurrentAcceleration        = (i)->FastGetSolutionStepValue(ACCELERATION, 0);
               array_1d<double, 3 > & CurrentVelocity            = (i)->FastGetSolutionStepValue(VELOCITY,     0);
               array_1d<double, 3 > & CurrentDisplacement  = (i)->FastGetSolutionStepValue(DISPLACEMENT);


               if (i->pGetDof(ACCELERATION_X)->IsFixed() )
               {
                  CurrentDisplacement[0] = PreviousDisplacement[0] + DeltaTime * PreviousVelocity[0] + std::pow(DeltaTime, 2) * ( 0.5 * (1.0 -  2.0 * this->mNewmark.beta) * PreviousAcceleration[0] + this->mNewmark.beta * CurrentAcceleration[0]);
               }
               else if (i->pGetDof(VELOCITY_X)->IsFixed() )
               {
                  CurrentDisplacement[0] = PreviousDisplacement[0] + 0.5 * DeltaTime * (PreviousVelocity[0] + CurrentVelocity[0]) + 0.5 * std::pow(DeltaTime, 2) * PreviousAcceleration[0];
               }
               else if (i->pGetDof(DISPLACEMENT_X)->IsFixed() == false)
               {
                  //CurrentDisplacement[0] = PreviousDisplacement[0] + DeltaTime * PreviousVelocity[0] + 0.5 * std::pow(DeltaTime, 2) * PreviousAcceleration[0];
               }


               if (i->pGetDof(ACCELERATION_Y)->IsFixed() )
               {
                  CurrentDisplacement[1] = PreviousDisplacement[1] + DeltaTime * PreviousVelocity[1] + std::pow(DeltaTime, 2) * ( 0.5 * (1.0 -  2.0 * this->mNewmark.beta) * PreviousAcceleration[1] + this->mNewmark.beta * CurrentAcceleration[1]);
               }
               else if (i->pGetDof(VELOCITY_Y)->IsFixed() )
               {
                  CurrentDisplacement[1] = PreviousDisplacement[1] + 0.5 * DeltaTime * (PreviousVelocity[1] + CurrentVelocity[1]) + 0.5 * std::pow(DeltaTime, 2) * PreviousAcceleration[1] ;
               }
               else if (i->pGetDof(DISPLACEMENT_Y)->IsFixed() == false )
               {
                  //CurrentDisplacement[1] = PreviousDisplacement[1] + DeltaTime * PreviousVelocity[1] + 0.5 * std::pow(DeltaTime, 2) * PreviousAcceleration[1];
               }

               // For 3D cases
               if (i->HasDofFor(DISPLACEMENT_Z))
               {
                  if (i->pGetDof(ACCELERATION_Z)->IsFixed() )
                  {
                     CurrentDisplacement[2] = PreviousDisplacement[2] + DeltaTime * PreviousVelocity[2] + std::pow(DeltaTime, 2) * ( 0.5 * (1.0 -  2.0 * this->mNewmark.beta) * PreviousAcceleration[2] + this->mNewmark.beta * CurrentAcceleration[2]);
                  }
                  else if (i->pGetDof(VELOCITY_Z)->IsFixed() )
                  {
                     CurrentDisplacement[2] = PreviousDisplacement[2] + 0.5 * DeltaTime * (PreviousVelocity[2] + CurrentVelocity[2]) + 0.5 * std::pow(DeltaTime, 2) * PreviousAcceleration[2] ;
                  }
                  else if (i->pGetDof(DISPLACEMENT_Z)->IsFixed() == false)
                  {
                     //CurrentDisplacement[2] = PreviousDisplacement[2] + DeltaTime * PreviousVelocity[2] + 0.5 * std::pow(DeltaTime, 2) * PreviousAcceleration[2];
                  }
               }


               if (i->HasDofFor(PRESSURE))
               {
                  double& PreviousPressure    = (i)->FastGetSolutionStepValue(PRESSURE, 1);
                  double& CurrentPressure     = (i)->FastGetSolutionStepValue(PRESSURE);

                  if ((i->pGetDof(PRESSURE))->IsFixed() == false)
                     CurrentPressure = PreviousPressure;

                  //std::cout<<" PressureCur [1] "<<CurrentPressure<<" PressurePre [1] "<<PreviousPressure<<" ID "<<i->Id()<<std::endl;
               }

               if (i->HasDofFor(JACOBIAN))
               {
                  double& PreviousJacobian    = (i)->FastGetSolutionStepValue(JACOBIAN, 1);
                  double& CurrentJacobian    = (i)->FastGetSolutionStepValue(JACOBIAN);

                  if ((i->pGetDof(JACOBIAN))->IsFixed() == false)
                     CurrentJacobian = PreviousJacobian;

               }


               //updating time derivatives ::: please note that displacements and its time derivatives
               //can not be consistently fixed separately

               noalias(DeltaDisplacement) = CurrentDisplacement - PreviousDisplacement;

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


               array_1d<double, 3 > & PreviousWaterAcceleration  = (i)->FastGetSolutionStepValue(WATER_ACCELERATION, 1);
               array_1d<double, 3 > & CurrentWaterVelocity       = (i)->FastGetSolutionStepValue(WATER_VELOCITY);
               array_1d<double, 3 > & CurrentWaterAcceleration   = (i)->FastGetSolutionStepValue(WATER_ACCELERATION);

               if (i->pGetDof(WATER_ACCELERATION_X)->IsFixed() )
               {               noalias(DeltaWaterDisplacement) = CurrentWaterDisplacement - PreviousWaterDisplacement;
                  CurrentWaterDisplacement[0] = PreviousWaterDisplacement[0] + DeltaTime * PreviousWaterVelocity[0] + std::pow(DeltaTime, 2) * ( 0.5 * (1.0 -  2.0 * this->mNewmark.beta) * PreviousWaterAcceleration[0] + this->mNewmark.beta * CurrentWaterAcceleration[0]);
               }
               else if (i->pGetDof(WATER_VELOCITY_X)->IsFixed() )
               {
                  CurrentWaterDisplacement[0] = PreviousWaterDisplacement[0] + 0.5 * DeltaTime * (PreviousWaterVelocity[0] + CurrentWaterVelocity[0]) + 0.5 * std::pow(DeltaTime, 2) * PreviousWaterAcceleration[0];
               }
               else if (i->pGetDof(WATER_DISPLACEMENT_X)->IsFixed() == false)
               {
                  //CurrentWaterDisplacement[0] = PreviousWaterDisplacement[0] + DeltaTime * PreviousWaterVelocity[0] + 0.5 * std::pow(DeltaTime, 2) * PreviousWaterAcceleration[0];
               }


               if (i->pGetDof(WATER_ACCELERATION_Y)->IsFixed() )
               {
                  CurrentWaterDisplacement[1] = PreviousWaterDisplacement[1] + DeltaTime * PreviousWaterVelocity[1] + std::pow(DeltaTime, 2) * ( 0.5 * (1.0 -  2.0 * this->mNewmark.beta) * PreviousWaterAcceleration[1] + this->mNewmark.beta * CurrentWaterAcceleration[1]);
               }
               else if (i->pGetDof(WATER_VELOCITY_Y)->IsFixed() )
               {
                  CurrentWaterDisplacement[1] = PreviousWaterDisplacement[1] + 0.5 * DeltaTime * (PreviousWaterVelocity[1] + CurrentWaterVelocity[1]) + 0.5 * std::pow(DeltaTime, 2) * PreviousWaterAcceleration[1] ;
               }
               else if (i->pGetDof(WATER_DISPLACEMENT_Y)->IsFixed() == false)
               {
                  //CurrentWaterDisplacement[1] = PreviousWaterDisplacement[1] + DeltaTime * PreviousWaterVelocity[1] + 0.5 * std::pow(DeltaTime, 2) * PreviousWaterAcceleration[1];
               }

               // For 3D cases
               if (i->HasDofFor(WATER_DISPLACEMENT_Z))
               {
                  if (i->pGetDof(WATER_ACCELERATION_Z)->IsFixed() )
                  {
                     CurrentWaterDisplacement[2] = PreviousWaterDisplacement[2] + DeltaTime * PreviousWaterVelocity[2] + std::pow(DeltaTime, 2) * ( 0.5 * (1.0 -  2.0 * this->mNewmark.beta) * PreviousWaterAcceleration[2] + this->mNewmark.beta * CurrentWaterAcceleration[2]);
                  }
                  else if (i->pGetDof(WATER_VELOCITY_Z)->IsFixed() )
                  {
                     CurrentWaterDisplacement[2] = PreviousWaterDisplacement[2] + 0.5 * DeltaTime * (PreviousWaterVelocity[2] + CurrentWaterVelocity[2]) + 0.5 * std::pow(DeltaTime, 2) * PreviousWaterAcceleration[2] ;
                  }
                  else if (i->pGetDof(WATER_DISPLACEMENT_Z)->IsFixed() == false)
                  {
                     //CurrentWaterDisplacement[2] = PreviousWaterDisplacement[2] + DeltaTime * PreviousWaterVelocity[2] + 0.5 * std::pow(DeltaTime, 2) * PreviousWaterAcceleration[2];
                  }
               }

               if (i->HasDofFor(WATER_PRESSURE))
               {
                  const double& PreviousWaterPressure    = (i)->FastGetSolutionStepValue(WATER_PRESSURE, 1);
                  const double& PreviousWaterPressureVelocity    = (i)->FastGetSolutionStepValue(WATER_PRESSURE_VELOCITY, 1);
                  const double& PreviousWaterPressureAcceleration    = (i)->FastGetSolutionStepValue(WATER_PRESSURE_ACCELERATION, 1);
                  double& CurrentWaterPressure     = (i)->FastGetSolutionStepValue(WATER_PRESSURE);
                  double& CurrentWaterPressureVelocity     = (i)->FastGetSolutionStepValue(WATER_PRESSURE_VELOCITY);
                  double& CurrentWaterPressureAcceleration     = (i)->FastGetSolutionStepValue(WATER_PRESSURE_ACCELERATION);

                  double DeltaWaterPressure = CurrentWaterPressure - PreviousWaterPressure;

                  if ((i->pGetDof(WATER_PRESSURE))->IsFixed() == false)
                     CurrentWaterPressure = PreviousWaterPressure;

                  UpdateVelocityScalar     ( CurrentWaterPressureVelocity, DeltaWaterPressure, PreviousWaterPressureVelocity, PreviousWaterPressureAcceleration);
                  UpdateAccelerationScalar ( CurrentWaterPressureAcceleration, DeltaWaterPressure, PreviousWaterPressureVelocity, PreviousWaterPressureAcceleration);


               }

               //updating time derivatives ::: please note that WaterDisplacements and its time derivatives
               //can not be consistently fixed separately

               noalias(DeltaWaterDisplacement) = CurrentWaterDisplacement - PreviousWaterDisplacement;

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
         //*********************************************************************************
         //Updating first time Derivative. For scalar variables.
         //*********************************************************************************

         inline void UpdateVelocityScalar(double  & CurrentVelocity,
               const double & DeltaDisplacement,
               const double & PreviousVelocity,
               const double & PreviousAcceleration)
         {

            CurrentVelocity =  (this->mNewmark.c1 * DeltaDisplacement - this->mNewmark.c4 * PreviousVelocity
                  - this->mNewmark.c5 * PreviousAcceleration) * this->mNewmark.static_dynamic;

         }


         //*********************************************************************************
         //Updating second time Derivative. For scalar variables.
         //*********************************************************************************

         inline void UpdateAccelerationScalar(double & CurrentAcceleration,
               const double & DeltaDisplacement,
               const double & PreviousVelocity,
               const double & PreviousAcceleration)
         {

            CurrentAcceleration =  (this->mNewmark.c0 * DeltaDisplacement - this->mNewmark.c2 * PreviousVelocity
                  -  this->mNewmark.c3 * PreviousAcceleration) * this->mNewmark.static_dynamic;


         }
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


