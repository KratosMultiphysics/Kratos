//
//   Project Name:        					KratosDamApplication $
//   Last modified by:  L Gracia  $Author:    Ignasi de Pouplana $
//   Date: March 2016     		  $Date:           February 2016 $
//   Revision:            		  $Revision:                 1.0 $
//

#if !defined(KRATOS_POINT_LOAD_TABLE_INTERPOLATION_PROCESS )
#define  KRATOS_POINT_LOAD_TABLE_INTERPOLATION_PROCESS

#include <cmath>

#include "includes/model_part.h"
#include "processes/process.h"

#include "dam_application_variables.h"

namespace Kratos
{

class PointLoadTableInterpolationProcess : public Process
{
    
public:

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Constructor
    PointLoadTableInterpolationProcess(ModelPart& r_model_part, double time_unit_converter) : mr_model_part(r_model_part)
    {
        mtime_unit_converter = time_unit_converter;
    }

    //------------------------------------------------------------------------------------
    
    // Destructor
    virtual ~PointLoadTableInterpolationProcess(){}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void Execute()
    {
        double time = mr_model_part.GetProcessInfo()[TIME];
        time = time/mtime_unit_converter;
        
        const double& InterpolatedPointLoadX = mr_model_part.pGetTable(5)->GetValue(time);
        
        for(ModelPart::NodeIterator i = mr_model_part.pGetMesh(5)->NodesBegin(); i <mr_model_part.pGetMesh(5)->NodesEnd(); i++)
        {
            (i)->FastGetSolutionStepValue(POINT_LOAD_X) = InterpolatedPointLoadX;
        }

        const double& InterpolatedPointLoadY = mr_model_part.pGetTable(6)->GetValue(time);
        
        for(ModelPart::NodeIterator i = mr_model_part.pGetMesh(6)->NodesBegin(); i <mr_model_part.pGetMesh(6)->NodesEnd(); i++)
        {
            (i)->FastGetSolutionStepValue(POINT_LOAD_Y) = InterpolatedPointLoadY;
        }

        const double& InterpolatedPointLoadZ = mr_model_part.pGetTable(7)->GetValue(time);
        
        for(ModelPart::NodeIterator i = mr_model_part.pGetMesh(7)->NodesBegin(); i <mr_model_part.pGetMesh(7)->NodesEnd(); i++)
        {
            (i)->FastGetSolutionStepValue(POINT_LOAD_Z) = InterpolatedPointLoadZ;
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

#endif /* KRATOS_POINT_LOAD_TABLE_INTERPOLATION_PROCESS defined */
