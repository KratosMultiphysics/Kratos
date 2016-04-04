//
//   Project Name:        KratosPoromechanicsApplication $
//   Last modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:           February 2016 $
//   Revision:            $Revision:                 1.0 $
//

#if !defined(KRATOS_NORMAL_LOAD_TABLE_INTERPOLATION_PROCESS )
#define  KRATOS_NORMAL_LOAD_TABLE_INTERPOLATION_PROCESS

#include <cmath>

#include "includes/model_part.h"
#include "processes/process.h"

#include "poromechanics_application_variables.h"

namespace Kratos
{

class NormalLoadTableInterpolationProcess : public Process
{
    
public:

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Constructor
    NormalLoadTableInterpolationProcess(ModelPart& r_model_part) : mr_model_part(r_model_part) {}

    //------------------------------------------------------------------------------------
    
    // Destructor
    virtual ~NormalLoadTableInterpolationProcess() {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void Execute()
    {
        const double& Time = mr_model_part.GetProcessInfo()[TIME];
        
        const double& InterpolatedNormalLoad = mr_model_part.pGetTable(11)->GetValue(Time);
        
        for(ModelPart::NodeIterator i = mr_model_part.pGetMesh(11)->NodesBegin(); i <mr_model_part.pGetMesh(11)->NodesEnd(); i++)
        {
            (i)->FastGetSolutionStepValue(NORMAL_CONTACT_STRESS) = InterpolatedNormalLoad;
        }
    }
    
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    // Member Variables

    ModelPart& mr_model_part;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
};//Class

} /* namespace Kratos.*/

#endif /* KRATOS_NORMAL_LOAD_TABLE_INTERPOLATION_PROCESS defined */
