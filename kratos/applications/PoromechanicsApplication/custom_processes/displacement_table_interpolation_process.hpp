//
//   Project Name:        KratosPoromechanicsApplication $
//   Last modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:           February 2016 $
//   Revision:            $Revision:                 1.0 $
//

#if !defined(KRATOS_DISPLACEMENT_TABLE_INTERPOLATION_PROCESS )
#define  KRATOS_DISPLACEMENT_TABLE_INTERPOLATION_PROCESS

#include <cmath>

#include "includes/model_part.h"
#include "processes/process.h"

#include "poromechanics_application_variables.h"

namespace Kratos
{

class DisplacementTableInterpolationProcess : public Process
{
    
public:

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Constructor
    DisplacementTableInterpolationProcess(ModelPart& r_model_part) : mr_model_part(r_model_part) {}

    //------------------------------------------------------------------------------------
    
    // Destructor
    virtual ~DisplacementTableInterpolationProcess() {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void Execute()
    {
        const double& Time = mr_model_part.GetProcessInfo()[TIME];
        
        const double& InterpolatedDisplacementX = mr_model_part.pGetTable(1)->GetValue(Time);
        
        for(ModelPart::NodeIterator i = mr_model_part.pGetMesh(1)->NodesBegin(); i <mr_model_part.pGetMesh(1)->NodesEnd(); i++)
        {
            (i)->FastGetSolutionStepValue(DISPLACEMENT_X) = InterpolatedDisplacementX;
        }

        const double& InterpolatedDisplacementY = mr_model_part.pGetTable(2)->GetValue(Time);
        
        for(ModelPart::NodeIterator i = mr_model_part.pGetMesh(2)->NodesBegin(); i <mr_model_part.pGetMesh(2)->NodesEnd(); i++)
        {
            (i)->FastGetSolutionStepValue(DISPLACEMENT_Y) = InterpolatedDisplacementY;
        }

        const double& InterpolatedDisplacementZ = mr_model_part.pGetTable(3)->GetValue(Time);
        
        for(ModelPart::NodeIterator i = mr_model_part.pGetMesh(3)->NodesBegin(); i <mr_model_part.pGetMesh(3)->NodesEnd(); i++)
        {
            (i)->FastGetSolutionStepValue(DISPLACEMENT_Z) = InterpolatedDisplacementZ;
        }
    }
    
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    // Member Variables

    ModelPart& mr_model_part;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
};//Class

} /* namespace Kratos.*/

#endif /* KRATOS_DISPLACEMENT_TABLE_INTERPOLATION_PROCESS defined */
