//
//   Project Name:        KratosPfemSolidMechanicsApplication $
//   Last modified by:    $Author:                JMCarbonell $
//   Date:                $Date:                    July 2013 $
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



namespace Kratos
{

   // blah blah blah
   class SetMechanicalInitialStateProcess
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

         SetMechanicalInitialStateProcess(ModelPart& rModelPart);

         SetMechanicalInitialStateProcess(ModelPart& rModelPart, const bool rGravity, const double rSV = 0.0, const double rSH = 0.0);


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

         virtual void ExecuteInitialize();

         virtual void ExecuteFinalizeSolutionStep();


      protected:

         void SetInitialMechanicalState(ModelPart& rModelPart, int EchoLevel = 0);

         void SetInitialMechanicalStateConstant(ModelPart& rModelPart, double S1, double S2, int EchoLevel = 0);

         void SetMechanicalState(ModelPart& rModelPart, const unsigned int& rMeshId, int& EchoLevel, const double& rYmax);

         void SetMechanicalStateUwP(ModelPart& rModelPart, const unsigned int& rMeshId, int& EchoLevel, const double& rYmax);

         void SetMechanicalStateU(ModelPart& rModelPart, const unsigned int& rMeshId, int& EchoLevel, const double& rYmax);

         void SetMechanicalStateUP(ModelPart& rModelPart, const unsigned int& rMeshId, int& EchoLevel, const double& rYmax);

         void SetMechanicalStateConstantUP(ModelPart& rModelPart, const unsigned int& rMeshId, const double& rS1, const double& rS2, int& EchoLevel);

         void SetMechanicalStateConstant(ModelPart& rModelPart, const unsigned int& rMeshId, const double& rS1, const double& rS2, int& EchoLevel);

      private:

         // member variables

         bool mGravity;

         std::vector<double> mInitialStress;

         ModelPart mrModelPart;

   }; //end class SetMechanicalInitialStateProcess

   }  // END namespace Kratos

#endif //KRATOS_SET_MECHANICAL_INITIAL_STATE_PROCESS_H_INCLUDED
