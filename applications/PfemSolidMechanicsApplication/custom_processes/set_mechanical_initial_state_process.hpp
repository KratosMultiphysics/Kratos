//
//   Project Name:        KratosPfemSolidMechanicsApplication $
//   Created by:          $Author:                  LMonforte $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                    July 2015 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined ( KRATOS_SET_MECHANICAL_INITIAL_STATE_PROCESS_H_INCLUDED )
#define  KRATOS_SET_MECHANICAL_INITIAL_STATE_PROCESS_H_INCLUDED


/* System includes */


/* External includes */


/* Project includes */
#include "processes/process.h"
#include "includes/model_part.h"
#include "includes/kratos_flags.h"
#include "utilities/math_utils.h"
#include "pfem_solid_mechanics_application_variables.h"

#include "includes/kratos_parameters.h"


namespace Kratos
{

   class KRATOS_API(PFEM_SOLID_MECHANICS_APPLICATION) SetMechanicalInitialStateProcess
      : public Process 
   {
      public:


         /**@name Type Definitions */
         /*@{ */

         // Pointer definition of Process
         KRATOS_CLASS_POINTER_DEFINITION( SetMechanicalInitialStateProcess );


         typedef ModelPart::NodesContainerType               NodesArrayType;
         typedef ModelPart::ConditionsContainerType ConditionsContainerType;
         typedef ModelPart::MeshType                               MeshType;
         /*@} */
         /**@name Life Cycle
          */
         /*@{ */

         // Constructor.

         //SetMechanicalInitialStateProcess(ModelPart& rModelPart);

         //SetMechanicalInitialStateProcess(ModelPart& rModelPart, const bool rGravity, const double rSV = 0.0, const double rSH = 0.0, const double rWaterPressure = 0.0, const bool rYmaxBool = false, const double rYmax = 0.0, const double rWaterLoad = 0.0);

         SetMechanicalInitialStateProcess( ModelPart& rModelPart, Parameters rParameters);


         /** Destructor.
          */

         virtual ~SetMechanicalInitialStateProcess();

         /*@} */
         /**@name Operators
          */
         /*@{ */


         /*@} */
         /**@name Operations */
         /*@{ */



         void operator()()
         {
            Execute();
         }

         void Execute() override;

      protected:

         void SetInitialMechanicalState(ModelPart& rModelPart, int EchoLevel = 0);

         void SetInitialMechanicalStateConstant(ModelPart& rModelPart, double S1, double S2, double WaterPressure, int EchoLevel = 0);

         void SetMechanicalState(ModelPart& rModelPart, int& EchoLevel, const double& rYmax);

         void SetMechanicalStateUwP(ModelPart& rModelPart, int& EchoLevel, const double& rYmax);

         void SetMechanicalStateU(ModelPart& rModelPart, int& EchoLevel, const double& rYmax);

         void SetMechanicalStateUP(ModelPart& rModelPart, int& EchoLevel, const double& rYmax);

         void SetMechanicalStateConstantUP(ModelPart& rModelPart, const double& rS1, const double& rS2, int& EchoLevel);

         void SetMechanicalStateConstant(ModelPart& rModelPart, const double& rS1, const double& rS2, const double& rWaterPressure, int& EchoLevel);

      private:

         // member variables

         bool mGravity;

         std::vector<double> mInitialStress;

         double mInitialWaterPressure; 

         ModelPart& mrModelPart;

         bool mSurfaceLoadBool; 

         double mSurfaceLoad;

         double mWaterLoad;

   }; //end class SetMechanicalInitialStateProcess

   }  // END namespace Kratos

#endif //KRATOS_SET_MECHANICAL_INITIAL_STATE_PROCESS_H_INCLUDED
