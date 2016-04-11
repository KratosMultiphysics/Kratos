//
//   Project Name:        KratosDamApplication    $
//   Last modified by:    $Author: Lorenzo Gracia $
//   Date:                $Date:     January 2016 $
//   Revision:            $Revision:          0.0 $
//

#if !defined(KRATOS_EXACT_BOFANG_EVOLUTION_CONDITIONS_TEMPERATURE_PROCESS )
#define  KRATOS_EXACT_BOFANG_EVOLUTION_CONDITIONS_TEMPERATURE_PROCESS

#include <cmath>

#include "includes/model_part.h"
#include "processes/process.h"

#include "dam_application_variables.h"

namespace Kratos
{

class ExactBofangEvolutionConditionsTemperatureProcess : public Process
{
    
public:

    typedef double argument_type; // To be STL conformance.
    typedef double result_type; // To be STL conformance.

    typedef boost::array<result_type, 1>  result_row_type;

    typedef std::pair<argument_type, result_row_type> RecordType;

    typedef std::vector<RecordType> TableContainerType;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Constructor
    ExactBofangEvolutionConditionsTemperatureProcess(ModelPart& r_model_part, double time_unit_converter) : mr_model_part(r_model_part)
    {
        mlast_id = 0;
        mtime_unit_converter = time_unit_converter;
    }

    //------------------------------------------------------------------------------------
    
    // Destructor
    virtual ~ExactBofangEvolutionConditionsTemperatureProcess(){}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void Execute()
    {
        double time = mr_model_part.GetProcessInfo()[TIME];
        time = time/mtime_unit_converter;
        double delta_time = mr_model_part.GetProcessInfo()[DELTA_TIME];
        delta_time = delta_time/mtime_unit_converter;
                
        const TableContainerType& table_month = mr_model_part.pGetTable(15)->Data();          // Table with information about months, used for the computations of Bofang formulation
        const TableContainerType& table_water_level = mr_model_part.pGetTable(16)->Data();    // Table with information about water level in each time
        const TableContainerType& table_outer_temp = mr_model_part.pGetTable(17)->Data();     // Table with information about outer temperature
        const TableContainerType& table_ref_temp = mr_model_part.pGetTable(18)->Data();       // Table with information about reference, temperature 
                
        if( (time + delta_time*1.0e-10) >= table_month[mlast_id+1].first )
        {
            mlast_id = mlast_id + 1;
            
            mr_model_part.GetProcessInfo()[REFERENCE_TEMPERATURE] = table_ref_temp[mlast_id].second[0];
                
            this->Bofang_conditions(table_month, table_water_level, table_outer_temp); 
            
        }
    }
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    // Member Variables

    ModelPart& mr_model_part;
    std::size_t mlast_id;
    double mtime_unit_converter;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void Bofang_conditions(TableContainerType const& table_month, TableContainerType const& table_water_level, TableContainerType const& table_outer_temp)
    {
        // We have to pick up the provided values by the mesh
        std::string direction = (*(mr_model_part.pGetMesh(15)))[GRAVITY_DIRECTION];
        const double& coordinate_base = (*(mr_model_part.pGetMesh(15)))[COORDINATE_BASE_DAM];
        const double& surface_temp = (*(mr_model_part.pGetMesh(15)))[SURFACE_TEMP];
        const double& bottom_temp = (*(mr_model_part.pGetMesh(15)))[BOTTOM_TEMP];
        const double& height_dam = (*(mr_model_part.pGetMesh(15)))[HEIGHT_DAM];
        const double& amplitude = (*(mr_model_part.pGetMesh(15)))[AMPLITUDE];
        const double& freq = 0.5235987756;
        const double& day = (*(mr_model_part.pGetMesh(15)))[DAY_MAXIMUM];

        double aux, aux1;

        // Values from the table
        const double& time_month =  table_month[mlast_id].second[0];                                // Time (month) Bofang
        const double& water_level = table_water_level[mlast_id].second[0];                          // Water Level
        const double& outer_temp =  table_outer_temp[mlast_id].second[0];                           // Outer temperature

        for(ModelPart::NodeIterator i = mr_model_part.pGetMesh(15)->NodesBegin(); i <mr_model_part.pGetMesh(15)->NodesEnd(); i++) 
        {
            double& Temperature  = (i)->FastGetSolutionStepValue(TEMPERATURE);
            
            if( direction == "X")
            {
                aux = (coordinate_base + water_level) - (i)->X();
                if(aux >= 0.0)
                {
                    aux1 = ((bottom_temp-(surface_temp*exp(-0.04*height_dam)))/(1-(exp(-0.04*height_dam))));
                    Temperature = (aux1+((surface_temp-aux1)*(exp(-0.04*aux)))+(amplitude*(exp(-0.018*aux))*(cos(freq*(time_month-(day/30.0)-2.15+(1.30*exp(-0.085*aux)))))));
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
                    Temperature = (aux1+((surface_temp-aux1)*(exp(-0.04*aux)))+(amplitude*(exp(-0.018*aux))*(cos(freq*(time_month-(day/30.0)-2.15+(1.30*exp(-0.085*aux)))))));
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
                    Temperature = (aux1+((surface_temp-aux1)*(exp(-0.04*aux)))+(amplitude*(exp(-0.018*aux))*(cos(freq*(time_month-(day/30.0)-2.15+(1.30*exp(-0.085*aux)))))));
                }
                else
                    Temperature = outer_temp;
            }                    
        }
    }
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
};//Class

} /* namespace Kratos.*/

#endif /* KRATOS_EXACT_BOFANG_EVOLUTION_CONDITIONS_TEMPERATURE_PROCESS defined */
