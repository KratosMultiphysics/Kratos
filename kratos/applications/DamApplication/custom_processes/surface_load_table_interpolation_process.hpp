//
//   Project Name:        					KratosDamApplication $
//   Last modified by:  L Gracia  $Author:    Ignasi de Pouplana $
//   Date: March 2016     		  $Date:           February 2016 $
//   Revision:            		  $Revision:                 1.0 $
//

#if !defined(KRATOS_SURFACE_LOAD_TABLE_INTERPOLATION_PROCESS )
#define  KRATOS_SURFACE_LOAD_TABLE_INTERPOLATION_PROCESS

#include <cmath>

#include "includes/model_part.h"
#include "processes/process.h"

#include "dam_application_variables.h"

namespace Kratos
{

class SurfaceLoadTableInterpolationProcess : public Process
{
    
public:

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Constructor
    SurfaceLoadTableInterpolationProcess(ModelPart& r_model_part, double time_unit_converter) : mr_model_part(r_model_part)
    {
        mtime_unit_converter = time_unit_converter;
    }

    //------------------------------------------------------------------------------------
    
    // Destructor
    virtual ~SurfaceLoadTableInterpolationProcess() {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void Execute()
    {
        double time = mr_model_part.GetProcessInfo()[TIME];
        time = time/mtime_unit_converter;
        
        const double& InterpolatedSurfaceLoadX = mr_model_part.pGetTable(11)->GetValue(time);
        
        for(ModelPart::NodeIterator i = mr_model_part.pGetMesh(11)->NodesBegin(); i <mr_model_part.pGetMesh(11)->NodesEnd(); i++)
        {
            (i)->FastGetSolutionStepValue(SURFACE_LOAD_X) = InterpolatedSurfaceLoadX;
        }

        const double& InterpolatedSurfaceLoadY = mr_model_part.pGetTable(12)->GetValue(time);
        
        for(ModelPart::NodeIterator i = mr_model_part.pGetMesh(12)->NodesBegin(); i <mr_model_part.pGetMesh(12)->NodesEnd(); i++)
        {
            (i)->FastGetSolutionStepValue(SURFACE_LOAD_Y) = InterpolatedSurfaceLoadY;
        }

        const double& InterpolatedSurfaceLoadZ = mr_model_part.pGetTable(13)->GetValue(time);
        
        for(ModelPart::NodeIterator i = mr_model_part.pGetMesh(13)->NodesBegin(); i <mr_model_part.pGetMesh(13)->NodesEnd(); i++)
        {
            (i)->FastGetSolutionStepValue(SURFACE_LOAD_Z) = InterpolatedSurfaceLoadZ;
        }        
    }
    
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    // Member Variables

    ModelPart& mr_model_part;
    double mtime_unit_converter;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
};//Class

} /* namespace Kratos.*/

#endif /* KRATOS_SURFACE_LOAD_TABLE_INTERPOLATION_PROCESS defined */
