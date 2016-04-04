//
//   Project Name:        KratosPoromechanicsApplication $
//   Last modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:           February 2016 $
//   Revision:            $Revision:                 1.0 $
//

#if !defined(KRATOS_TANGENTIAL_LOAD_TABLE_INTERPOLATION_PROCESS )
#define  KRATOS_TANGENTIAL_LOAD_TABLE_INTERPOLATION_PROCESS

#include <cmath>

#include "includes/model_part.h"
#include "processes/process.h"

#include "poromechanics_application_variables.h"

namespace Kratos
{

class TangentialLoadTableInterpolationProcess : public Process
{
    
public:

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Constructor
    TangentialLoadTableInterpolationProcess(ModelPart& r_model_part) : mr_model_part(r_model_part) {}

    //------------------------------------------------------------------------------------
    
    // Destructor
    virtual ~TangentialLoadTableInterpolationProcess() {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void Execute()
    {
        const double& Time = mr_model_part.GetProcessInfo()[TIME];
        
        const double& InterpolatedTangentialLoad = mr_model_part.pGetTable(12)->GetValue(Time);
        
        for(ModelPart::NodeIterator i = mr_model_part.pGetMesh(12)->NodesBegin(); i <mr_model_part.pGetMesh(12)->NodesEnd(); i++)
        {
            (i)->FastGetSolutionStepValue(TANGENTIAL_CONTACT_STRESS) = InterpolatedTangentialLoad;
        }
    }
    
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    // Member Variables

    ModelPart& mr_model_part;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
};//Class

} /* namespace Kratos.*/

#endif /* KRATOS_TANGENTIAL_LOAD_TABLE_INTERPOLATION_PROCESS defined */
