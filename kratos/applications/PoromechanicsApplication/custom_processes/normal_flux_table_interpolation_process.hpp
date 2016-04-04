//
//   Project Name:        KratosPoromechanicsApplication $
//   Last modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:           February 2016 $
//   Revision:            $Revision:                 1.0 $
//

#if !defined(KRATOS_NORMAL_FLUX_TABLE_INTERPOLATION_PROCESS )
#define  KRATOS_NORMAL_FLUX_TABLE_INTERPOLATION_PROCESS

#include <cmath>

#include "includes/model_part.h"
#include "processes/process.h"

#include "poromechanics_application_variables.h"

namespace Kratos
{

class NormalFluxTableInterpolationProcess : public Process
{
    
public:

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Constructor
    NormalFluxTableInterpolationProcess(ModelPart& r_model_part) : mr_model_part(r_model_part) {}

    //------------------------------------------------------------------------------------
    
    // Destructor
    virtual ~NormalFluxTableInterpolationProcess() {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void Execute()
    {
        const double& Time = mr_model_part.GetProcessInfo()[TIME];
        
        const double& InterpolatedNormalFlux = mr_model_part.pGetTable(13)->GetValue(Time);
        
        for(ModelPart::NodeIterator i = mr_model_part.pGetMesh(13)->NodesBegin(); i <mr_model_part.pGetMesh(13)->NodesEnd(); i++)
        {
            (i)->FastGetSolutionStepValue(NORMAL_FLUID_FLUX) = InterpolatedNormalFlux;
        }
    }
    
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    // Member Variables

    ModelPart& mr_model_part;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
};//Class

} /* namespace Kratos.*/

#endif /* KRATOS_NORMAL_FLUX_TABLE_INTERPOLATION_PROCESS defined */
