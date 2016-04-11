//
//   Project Name:        					KratosDamApplication $
//   Last modified by:  L Gracia  $Author:    Ignasi de Pouplana $
//   Date: March 2016     		  $Date:           February 2016 $
//   Revision:            		  $Revision:                 1.0 $
//

#if !defined(KRATOS_NORMAL_LOAD_TABLE_INTERPOLATION_PROCESS )
#define  KRATOS_NORMAL_LOAD_TABLE_INTERPOLATION_PROCESS

#include <cmath>

#include "includes/model_part.h"
#include "processes/process.h"

#include "dam_application_variables.h"

namespace Kratos
{

class NormalLoadTableInterpolationProcess : public Process
{
    
public:

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Constructor
    NormalLoadTableInterpolationProcess(ModelPart& r_model_part,  double time_unit_converter) : mr_model_part(r_model_part)
    {
         mtime_unit_converter = time_unit_converter;
    }

    //------------------------------------------------------------------------------------
    
    // Destructor
    virtual ~NormalLoadTableInterpolationProcess() {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void Execute()
    {
        double time = mr_model_part.GetProcessInfo()[TIME];
        time = time/mtime_unit_converter;
        
        const double& InterpolatedNormalLoad = mr_model_part.pGetTable(14)->GetValue(time);
        
        for(ModelPart::NodeIterator i = mr_model_part.pGetMesh(14)->NodesBegin(); i <mr_model_part.pGetMesh(14)->NodesEnd(); i++)
        {
            (i)->FastGetSolutionStepValue(NEGATIVE_FACE_PRESSURE) = InterpolatedNormalLoad;
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

#endif /* KRATOS_NORMAL_LOAD_TABLE_INTERPOLATION_PROCESS defined */
