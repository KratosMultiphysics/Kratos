//
//   Project Name:        KratosDamApplication    $
//   Last modified by:    $Author: Lorenzo Gracia $
//   Date:                $Date:     January 2016 $
//   Revision:            $Revision:          0.0 $
//

#if !defined(KRATOS_LINEAR_EVOLUTION_CONDITIONS_TEMPERATURE_PROCESS )
#define  KRATOS_LINEAR_EVOLUTION_CONDITIONS_TEMPERATURE_PROCESS

#include <cmath>

#include "includes/model_part.h"
#include "processes/process.h"

#include "dam_application.h"

namespace Kratos
{

class LinearEvolutionConditionsTemperatureProcess : public Process
{
    
public:

    typedef double result_type; // To be STL conformance.

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Constructor
    LinearEvolutionConditionsTemperatureProcess(ModelPart& r_model_part) : mr_model_part(r_model_part)
    {
    }

    //------------------------------------------------------------------------------------
    
    // Destructor
    virtual ~LinearEvolutionConditionsTemperatureProcess(){}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void Execute()
    {
        this->Uniform_conditions(); 

    }
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    // Member Variables

    ModelPart& mr_model_part;


//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    void Uniform_conditions()
    {
                
        for(ModelPart::NodeIterator i = mr_model_part.pGetMesh(2)->NodesBegin(); i <mr_model_part.pGetMesh(2)->NodesEnd(); i++)
        {
            double& Temperature  = (i)->FastGetSolutionStepValue(TEMPERATURE);
            Temperature += (i)->FastGetSolutionStepValue(IMPOSED_TEMPERATURE);
        }
    }
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

};//Class

} /* namespace Kratos.*/

#endif /* KRATOS_LINEAR_EVOLUTION_CONDITIONS_TEMPERATURE_PROCESS defined */
