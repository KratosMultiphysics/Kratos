//
//   Project Name:        KratosPoromechanicsApplication $
//   Last modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:           February 2016 $
//   Revision:            $Revision:                 1.0 $
//

#if !defined(KRATOS_FACE_LOAD_TABLE_INTERPOLATION_PROCESS )
#define  KRATOS_FACE_LOAD_TABLE_INTERPOLATION_PROCESS

#include <cmath>

#include "includes/model_part.h"
#include "processes/process.h"

#include "poromechanics_application_variables.h"

namespace Kratos
{

class FaceLoadTableInterpolationProcess : public Process
{
    
public:

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Constructor
    FaceLoadTableInterpolationProcess(ModelPart& r_model_part) : mr_model_part(r_model_part) {}

    //------------------------------------------------------------------------------------
    
    // Destructor
    virtual ~FaceLoadTableInterpolationProcess() {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void Execute()
    {
        const double& Time = mr_model_part.GetProcessInfo()[TIME];
        
        const double& InterpolatedLineLoadX = mr_model_part.pGetTable(8)->GetValue(Time);
        
        for(ModelPart::NodeIterator i = mr_model_part.pGetMesh(8)->NodesBegin(); i <mr_model_part.pGetMesh(8)->NodesEnd(); i++)
        {
            (i)->FastGetSolutionStepValue(FACE_LOAD_X) = InterpolatedLineLoadX;
        }
        
        const double& InterpolatedLineLoadY = mr_model_part.pGetTable(9)->GetValue(Time);
        
        for(ModelPart::NodeIterator i = mr_model_part.pGetMesh(9)->NodesBegin(); i <mr_model_part.pGetMesh(9)->NodesEnd(); i++)
        {
            (i)->FastGetSolutionStepValue(FACE_LOAD_Y) = InterpolatedLineLoadY;
        }
        
        const double& InterpolatedLineLoadZ = mr_model_part.pGetTable(10)->GetValue(Time);
        
        for(ModelPart::NodeIterator i = mr_model_part.pGetMesh(10)->NodesBegin(); i <mr_model_part.pGetMesh(10)->NodesEnd(); i++)
        {
            (i)->FastGetSolutionStepValue(FACE_LOAD_Z) = InterpolatedLineLoadZ;
        }
    }
    
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    // Member Variables

    ModelPart& mr_model_part;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
};//Class

} /* namespace Kratos.*/

#endif /* KRATOS_FACE_LOAD_TABLE_INTERPOLATION_PROCESS defined */
