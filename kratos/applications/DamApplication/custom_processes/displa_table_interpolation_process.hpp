//
//   Project Name:        					KratosDamApplication $
//   Last modified by:  L Gracia  $Author:    Ignasi de Pouplana $
//   Date: March 2016     		  $Date:           February 2016 $
//   Revision:            		  $Revision:                 1.0 $
//

#if !defined(KRATOS_DISPLA_TABLE_INTERPOLATION_PROCESS )
#define  KRATOS_DISPLA_TABLE_INTERPOLATION_PROCESS

#include <cmath>

#include "includes/model_part.h"
#include "processes/process.h"

#include "dam_application_variables.h"

namespace Kratos
{

class DisplaTableInterpolationProcess : public Process
{
    
public:

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Constructor
    DisplaTableInterpolationProcess(ModelPart& r_model_part, double time_unit_converter) : mr_model_part(r_model_part) 
    {
        mtime_unit_converter = time_unit_converter;    
    }

    //------------------------------------------------------------------------------------
    
    // Destructor
    virtual ~DisplaTableInterpolationProcess() {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void Execute()
    {
        double time = mr_model_part.GetProcessInfo()[TIME];
        time = time/mtime_unit_converter;
        
        const double& InterpolatedDisplacementX = mr_model_part.pGetTable(1)->GetValue(time);
        
        for(ModelPart::NodeIterator i = mr_model_part.pGetMesh(1)->NodesBegin(); i <mr_model_part.pGetMesh(1)->NodesEnd(); i++)
        {
            (i)->FastGetSolutionStepValue(DISPLACEMENT_X) = InterpolatedDisplacementX;
        }

        const double& InterpolatedDisplacementY = mr_model_part.pGetTable(2)->GetValue(time);
        
        for(ModelPart::NodeIterator i = mr_model_part.pGetMesh(2)->NodesBegin(); i <mr_model_part.pGetMesh(2)->NodesEnd(); i++)
        {
            (i)->FastGetSolutionStepValue(DISPLACEMENT_Y) = InterpolatedDisplacementY;
        }

        const double& InterpolatedDisplacementZ = mr_model_part.pGetTable(3)->GetValue(time);
        
        for(ModelPart::NodeIterator i = mr_model_part.pGetMesh(3)->NodesBegin(); i <mr_model_part.pGetMesh(3)->NodesEnd(); i++)
        {
            (i)->FastGetSolutionStepValue(DISPLACEMENT_Z) = InterpolatedDisplacementZ;
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

#endif /* KRATOS_DISPLA_TABLE_INTERPOLATION_PROCESS defined */
