//
//   Project Name:        KratosPoromechanicsApplication $
//   Last modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:           February 2016 $
//   Revision:            $Revision:                 1.0 $
//

#if !defined(KRATOS_FORCE_TABLE_INTERPOLATION_PROCESS )
#define  KRATOS_FORCE_TABLE_INTERPOLATION_PROCESS

#include <cmath>

#include "includes/model_part.h"
#include "processes/process.h"

#include "poromechanics_application_variables.h"

namespace Kratos
{

class ForceTableInterpolationProcess : public Process
{
    
public:

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Constructor
    ForceTableInterpolationProcess(ModelPart& r_model_part) : mr_model_part(r_model_part) {}

    //------------------------------------------------------------------------------------
    
    // Destructor
    virtual ~ForceTableInterpolationProcess() {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void Execute()
    {
        const double& Time = mr_model_part.GetProcessInfo()[TIME];
        
        const double& InterpolatedPointLoadX = mr_model_part.pGetTable(5)->GetValue(Time);
        
        for(ModelPart::NodeIterator i = mr_model_part.pGetMesh(5)->NodesBegin(); i <mr_model_part.pGetMesh(5)->NodesEnd(); i++)
        {
            (i)->FastGetSolutionStepValue(FORCE_X) = InterpolatedPointLoadX;
        }

        const double& InterpolatedPointLoadY = mr_model_part.pGetTable(6)->GetValue(Time);
        
        for(ModelPart::NodeIterator i = mr_model_part.pGetMesh(6)->NodesBegin(); i <mr_model_part.pGetMesh(6)->NodesEnd(); i++)
        {
            (i)->FastGetSolutionStepValue(FORCE_Y) = InterpolatedPointLoadY;
        }

        const double& InterpolatedPointLoadZ = mr_model_part.pGetTable(7)->GetValue(Time);
        
        for(ModelPart::NodeIterator i = mr_model_part.pGetMesh(7)->NodesBegin(); i <mr_model_part.pGetMesh(7)->NodesEnd(); i++)
        {
            (i)->FastGetSolutionStepValue(FORCE_Z) = InterpolatedPointLoadZ;
        }
    }
    
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    // Member Variables

    ModelPart& mr_model_part;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
};//Class

} /* namespace Kratos.*/

#endif /* KRATOS_FORCE_TABLE_INTERPOLATION_PROCESS defined */
