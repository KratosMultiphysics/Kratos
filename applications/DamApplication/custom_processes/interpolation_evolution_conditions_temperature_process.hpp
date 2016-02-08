//
//   Project Name:        KratosDamApplication    $
//   Last modified by:    $Author: Lorenzo Gracia $
//   Date:                $Date:     January 2016 $
//   Revision:            $Revision:          0.0 $
//

#if !defined(KRATOS_INTERPOLATION_EVOLUTION_CONDITIONS_TEMPERATURE_PROCESS )
#define  KRATOS_INTERPOLATION_EVOLUTION_CONDITIONS_TEMPERATURE_PROCESS

#include <cmath>

#include "includes/model_part.h"
#include "processes/process.h"

#include "dam_application.h"

namespace Kratos
{

class InterpolationEvolutionConditionsTemperatureProcess : public Process
{
    
public:

    typedef double result_type; // To be STL conformance.

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Constructor
    InterpolationEvolutionConditionsTemperatureProcess(ModelPart& r_model_part, double time_unit_converter) : mr_model_part(r_model_part)
    {
        mtime_unit_converter = time_unit_converter;
    }

    //------------------------------------------------------------------------------------
    
    // Destructor
    virtual ~InterpolationEvolutionConditionsTemperatureProcess(){}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void Execute()
    {
        double time = mr_model_part.GetProcessInfo()[TIME];
        time = time/mtime_unit_converter;
        
        const result_type& month = mr_model_part.pGetTable(1)->GetValue(time);          // Interpolated data with information about months, used for the computations of Bofang formulation
        const result_type& water_level = mr_model_part.pGetTable(2)->GetValue(time);    // Interpolated data with information about water level in each time
        const result_type& outer_temp = mr_model_part.pGetTable(3)->GetValue(time);     // Interpolated data with information about outer temperature
        const result_type& ref_temp = mr_model_part.pGetTable(4)->GetValue(time);       // Interpolated data with information about reference, temperature 
        
        mr_model_part.GetProcessInfo()[REFERENCE_TEMPERATURE] = ref_temp;
        
        this->Bofang_conditions(month, water_level, outer_temp);
    
        this->Uniform_conditions(outer_temp); 

    }
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    // Member Variables

    ModelPart& mr_model_part;
    double mtime_unit_converter;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void Bofang_conditions(result_type const& month, result_type const& water_level, result_type const& outer_temp)
    {
        // We have to pick up the provided values by the mesh
        std::string direction = (*(mr_model_part.pGetMesh(1)))[GRAVITY_DIRECTION];
        const double& coordinate_base = (*(mr_model_part.pGetMesh(1)))[COORDINATE_BASE_DAM];
        const double& surface_temp = (*(mr_model_part.pGetMesh(1)))[SURFACE_TEMP];
        const double& bottom_temp = (*(mr_model_part.pGetMesh(1)))[BOTTOM_TEMP];
        const double& height_dam = (*(mr_model_part.pGetMesh(1)))[HEIGHT_DAM];
        const double& amplitude = (*(mr_model_part.pGetMesh(1)))[AMPLITUDE];
        const double& freq = (*(mr_model_part.pGetMesh(1)))[FREQUENCY];
        const double& day = (*(mr_model_part.pGetMesh(1)))[DAY_MAXIMUM];

        double aux, aux1;

        for(ModelPart::NodeIterator i = mr_model_part.pGetMesh(1)->NodesBegin(); i <mr_model_part.pGetMesh(1)->NodesEnd(); i++) 
        {
            double& Temperature  = (i)->FastGetSolutionStepValue(TEMPERATURE);
            
            if( direction == "X")
            {
                aux = (coordinate_base + water_level) - (i)->X();
                if(aux >= 0.0)
                {
                    aux1 = ((bottom_temp-(surface_temp*exp(-0.04*height_dam)))/(1-(exp(-0.04*height_dam))));
                    Temperature = (aux1+((surface_temp-aux1)*(exp(-0.04*aux)))+(amplitude*(exp(-0.018*aux))*(cos(freq*(month-(day/30.0)-2.15+(1.30*exp(-0.085*aux)))))));
                }
                else
                    Temperature = outer_temp;
            }
            else if( direction == "Y")
            {
                aux = (coordinate_base + water_level) - (i)->Y();
                if(aux >= 0.0)
                { 
                    aux1 = ((bottom_temp-(surface_temp*exp(-0.04*height_dam)))/(1-(exp(-0.04*height_dam))));
                    Temperature = (aux1+((surface_temp-aux1)*(exp(-0.04*aux)))+(amplitude*(exp(-0.018*aux))*(cos(freq*(month-(day/30.0)-2.15+(1.30*exp(-0.085*aux)))))));
                }
                else
                    Temperature = outer_temp;
            }
            else if( direction == "Z")
            {
                aux = (coordinate_base + water_level) - (i)->Z();
                if(aux >= 0.0)
                { 
                    aux1 = ((bottom_temp-(surface_temp*exp(-0.04*height_dam)))/(1-(exp(-0.04*height_dam))));
                    Temperature = (aux1+((surface_temp-aux1)*(exp(-0.04*aux)))+(amplitude*(exp(-0.018*aux))*(cos(freq*(month-(day/30.0)-2.15+(1.30*exp(-0.085*aux)))))));
                }
                else
                    Temperature = outer_temp;
            }                    
        }
    }
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void Uniform_conditions(result_type const& outer_temp)
    {
                
        for(ModelPart::NodeIterator i = mr_model_part.pGetMesh(2)->NodesBegin(); i <mr_model_part.pGetMesh(2)->NodesEnd(); i++)
        {
            double& Temperature  = (i)->FastGetSolutionStepValue(TEMPERATURE);
            Temperature = outer_temp;
        }
    }
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

};//Class

} /* namespace Kratos.*/

#endif /* KRATOS_INTERPOLATION_EVOLUTION_CONDITIONS_TEMPERATURE_PROCESS defined */
