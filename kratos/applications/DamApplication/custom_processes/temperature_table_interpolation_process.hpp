//
//   Project Name:        KratosDamApplication    $
//   Last modified by:    $Author: Lorenzo Gracia $
//   Date:                $Date:       March 2016 $
//   Revision:            $Revision:          0.0 $
//

#if !defined(KRATOS_TEMPERATURE_TABLE_INTERPOLATION_PROCESS )
#define  KRATOS_TEMPERATURE_TABLE_INTERPOLATION_PROCESS

#include <cmath>

#include "includes/model_part.h"
#include "processes/process.h"

#include "dam_application_variables.h"

namespace Kratos
{

class TemperatureTableInterpolationProcess : public Process
{
    
public:

    typedef double result_type; // To be STL conformance.

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Constructor
    TemperatureTableInterpolationProcess(ModelPart& r_model_part, double time_unit_converter) : mr_model_part(r_model_part)
    {
        mtime_unit_converter = time_unit_converter;
    }

    //------------------------------------------------------------------------------------
    
    // Destructor
    virtual ~TemperatureTableInterpolationProcess(){}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void Execute()
    {
        double time = mr_model_part.GetProcessInfo()[TIME];
        time = time/mtime_unit_converter;
        
        const result_type& uniform_temp = mr_model_part.pGetTable(4)->GetValue(time);    // Interpolated data with information about uniform temperature
        const result_type& ref_temp = mr_model_part.pGetTable(18)->GetValue(time);       // Interpolated data with information about reference, temperature 
        
        mr_model_part.GetProcessInfo()[REFERENCE_TEMPERATURE] = ref_temp;
    
        this->Uniform_conditions(uniform_temp); 

    }
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    // Member Variables

    ModelPart& mr_model_part;
    double mtime_unit_converter;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void Uniform_conditions(result_type const& uniform_temp)
    {
                
        for(ModelPart::NodeIterator i = mr_model_part.pGetMesh(4)->NodesBegin(); i <mr_model_part.pGetMesh(4)->NodesEnd(); i++)
        {
           (i)->FastGetSolutionStepValue(TEMPERATURE) = uniform_temp;
        }
    }
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

};//Class

} /* namespace Kratos.*/

#endif /* KRATOS_TEMPERATURE_TABLE_INTERPOLATION_PROCESS defined */
