//
//   Project Name:        KratosPfemSolidMechanicsApplication $
//   Created by:          $Author:                  LMonforte $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                    July 2018 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined ( KRATOS_NON_LOCAL_PLASTICITY_PROCESS_H_INCLUDED )
#define        KRATOS_NON_LOCAL_PLASTICITY_PROCESS_H_INCLUDED


/* System includes */


/* External includes */


/* Project includes */
#include "processes/process.h"
#include "includes/model_part.h"
#include "includes/kratos_flags.h"
#include "utilities/math_utils.h"
#include "constitutive_models_application_variables.h"

#include "includes/kratos_parameters.h"


namespace Kratos
{

   class KRATOS_API(CONSTITUTIVE_MODELS_APPLICATION) NonLocalPlasticityProcess
      : public Process 
      {

         protected:
            struct GaussPoint
            {
               GaussPoint() {}

               GaussPoint( ConstitutiveLaw::Pointer & pConstLaw,
                     array_1d<double, 3> rCoord) 
               {
                  pConstitutiveLaw = pConstLaw;
                  Coordinates = rCoord;
               }

               void AddNeighbour( const int & rID,
                     double weight)
               {
                  NeighbourGP.push_back( rID);
                  NeighbourWeight.push_back(weight);
               }

               ConstitutiveLaw::Pointer pConstitutiveLaw;
               array_1d<double, 3> Coordinates;
               
               std::vector< int > NeighbourGP;
               std::vector<double>  NeighbourWeight;

            };

         public:


            /**@name Type Definitions */
            /*@{ */

            // Pointer definition of Process
            KRATOS_CLASS_POINTER_DEFINITION( NonLocalPlasticityProcess );


            typedef ModelPart::NodesContainerType               NodesArrayType;
            typedef ModelPart::ConditionsContainerType ConditionsContainerType;
            typedef ModelPart::MeshType                               MeshType;
            /*@} */
            /**@name Life Cycle
             */
            /*@{ */

            // Constructor.


            NonLocalPlasticityProcess( ModelPart& rModelPart, Parameters rParameters);


            /** Destructor.
             */

            virtual ~NonLocalPlasticityProcess();

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




            void PerformGaussPointSearch( std::vector< GaussPoint > & rNeighbourGP, 
                  const double CharacteristicLength);

            double& ComputeWeightFunction( const double& rDistance, const double & rCharacteristicLength, double & rAlpha);

         private:

            // member variables

            ModelPart& mrModelPart;

            double mCharacteristicLength;

            std::vector< Variable<double> > mLocalVariables; 
            std::vector< Variable<double> > mNonLocalVariables; 


      }; //end class NonLocalPlasticityProcess

      }  // END namespace Kratos

#endif //KRATOS_NON_LOCAL_PLASTICITY_PROCESS_H_INCLUDED
