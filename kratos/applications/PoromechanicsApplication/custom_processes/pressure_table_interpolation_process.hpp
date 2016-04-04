//
//   Project Name:        KratosPoromechanicsApplication $
//   Last modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:           February 2016 $
//   Revision:            $Revision:                 1.0 $
//

#if !defined(KRATOS_PRESSURE_TABLE_INTERPOLATION_PROCESS )
#define  KRATOS_PRESSURE_TABLE_INTERPOLATION_PROCESS

#include <cmath>

#include "includes/model_part.h"
#include "processes/process.h"

#include "poromechanics_application_variables.h"

namespace Kratos
{

class PressureTableInterpolationProcess : public Process
{
    
public:

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Constructor
    PressureTableInterpolationProcess(ModelPart& r_model_part) : mr_model_part(r_model_part) {}

    //------------------------------------------------------------------------------------
    
    // Destructor
    virtual ~PressureTableInterpolationProcess() {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void Execute()
    {
        const double& Time = mr_model_part.GetProcessInfo()[TIME];
        
        const double& InterpolatedPressure = mr_model_part.pGetTable(4)->GetValue(Time);
        
        for(ModelPart::NodeIterator i = mr_model_part.pGetMesh(4)->NodesBegin(); i <mr_model_part.pGetMesh(4)->NodesEnd(); i++)
        {
            (i)->FastGetSolutionStepValue(WATER_PRESSURE) = InterpolatedPressure;
        }
    }
    
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    // Member Variables

    ModelPart& mr_model_part;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
};//Class

} /* namespace Kratos.*/

#endif /* KRATOS_PRESSURE_TABLE_INTERPOLATION_PROCESS defined */
