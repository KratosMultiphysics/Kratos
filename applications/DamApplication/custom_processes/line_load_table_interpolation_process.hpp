//
//   Project Name:        					KratosDamApplication $
//   Last modified by:  L Gracia  $Author:    Ignasi de Pouplana $
//   Date: March 2016     		  $Date:           February 2016 $
//   Revision:            		  $Revision:                 1.0 $
//

#if !defined(KRATOS_LINE_LOAD_TABLE_INTERPOLATION_PROCESS )
#define  KRATOS_LINE_LOAD_TABLE_INTERPOLATION_PROCESS

#include <cmath>

#include "includes/model_part.h"
#include "processes/process.h"

#include "dam_application_variables.h"

namespace Kratos
{

class LineLoadTableInterpolationProcess : public Process
{
    
public:

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Constructor
    LineLoadTableInterpolationProcess(ModelPart& r_model_part, double time_unit_converter) : mr_model_part(r_model_part)
    {
        mtime_unit_converter = time_unit_converter;
    }

    //------------------------------------------------------------------------------------
    
    // Destructor
    virtual ~LineLoadTableInterpolationProcess() {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void Execute()
    {
        double time = mr_model_part.GetProcessInfo()[TIME];
        time = time/mtime_unit_converter;
        
        const double& InterpolatedLineLoadX = mr_model_part.pGetTable(8)->GetValue(time);
        
        for(ModelPart::NodeIterator i = mr_model_part.pGetMesh(8)->NodesBegin(); i <mr_model_part.pGetMesh(8)->NodesEnd(); i++)
        {
            (i)->FastGetSolutionStepValue(LINE_LOAD_X) = InterpolatedLineLoadX;
        }
        
        const double& InterpolatedLineLoadY = mr_model_part.pGetTable(9)->GetValue(time);
        
        for(ModelPart::NodeIterator i = mr_model_part.pGetMesh(9)->NodesBegin(); i <mr_model_part.pGetMesh(9)->NodesEnd(); i++)
        {
            (i)->FastGetSolutionStepValue(LINE_LOAD_Y) = InterpolatedLineLoadY;
        }
        
        const double& InterpolatedLineLoadZ = mr_model_part.pGetTable(10)->GetValue(time);
        
        for(ModelPart::NodeIterator i = mr_model_part.pGetMesh(10)->NodesBegin(); i <mr_model_part.pGetMesh(10)->NodesEnd(); i++)
        {
            (i)->FastGetSolutionStepValue(LINE_LOAD_Z) = InterpolatedLineLoadZ;
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

#endif /* KRATOS_LINE_LOAD_TABLE_INTERPOLATION_PROCESS defined */
