//
//   Project Name:        KratosDamApplication    $
//   Last modified by:    $Author: Lorenzo Gracia $
//   Date:                $Date:     January 2016 $
//   Revision:            $Revision:          0.0 $
//

#if !defined(KRATOS_INTERPOLATION_EVOLUTION_CONDITIONS_LOAD_PROCESS )
#define  KRATOS_INTERPOLATION_EVOLUTION_CONDITIONS_LOAD_PROCESS

#include <cmath>

#include "includes/model_part.h"
#include "processes/process.h"

#include "dam_application.h"

namespace Kratos
{

class InterpolationEvolutionConditionsLoadProcess : public Process
{
    
public:

    typedef double result_type; // To be STL conformance.

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Constructor
    InterpolationEvolutionConditionsLoadProcess(ModelPart& r_model_part, double time_unit_converter) : mr_model_part(r_model_part)
    {
        mtime_unit_converter = time_unit_converter;
    }

    //------------------------------------------------------------------------------------
    
    // Destructor
    virtual ~InterpolationEvolutionConditionsLoadProcess(){}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void Execute()
    {
        double time = mr_model_part.GetProcessInfo()[TIME];
        time = time/mtime_unit_converter;
        
        const result_type& water_level = mr_model_part.pGetTable(2)->GetValue(time);    // Interpolated data with information about water level in each time

        this->Hydrostatic_pressure(water_level); 
    }
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    // Member Variables

    ModelPart& mr_model_part;
    double mtime_unit_converter;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void Hydrostatic_pressure(result_type const& water_level)
    {
        // We have to pick up the provided values by the mesh
        std::string direction = (*(mr_model_part.pGetMesh(3)))[GRAVITY_DIRECTION];
        const double& coordinate_base = (*(mr_model_part.pGetMesh(3)))[COORDINATE_BASE_DAM];
        const double& spe_weight =  (*(mr_model_part.pGetMesh(3)))[SPECIFIC_WEIGHT];
        
        double ref_coord;
                
        for(ModelPart::NodeIterator i = mr_model_part.pGetMesh(3)->NodesBegin(); i <mr_model_part.pGetMesh(3)->NodesEnd(); i++)
        {
            double& pressure  = (i)->FastGetSolutionStepValue(NORMAL_CONTACT_STRESS);
            
            if(direction=="X")
            {
                ref_coord = coordinate_base + water_level;
                pressure = spe_weight*(ref_coord- (i)->X());
                    if(pressure < 0.0)
                        pressure = 0.0; 
            }
            else if(direction=="Y")
            {
                ref_coord = coordinate_base + water_level;
                pressure = spe_weight*(ref_coord- (i)->Y());
                    if(pressure < 0.0)
                        pressure = 0.0; 
            }
            else if(direction=="Z")
            {
                ref_coord = coordinate_base + water_level;
                pressure = spe_weight*(ref_coord- (i)->Z());
                    if(pressure < 0.0)
                        pressure = 0.0; 
            }
        }
    }
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
};//Class

} /* namespace Kratos.*/

#endif /* KRATOS_INTERPOLATION_EVOLUTION_CONDITIONS_LOAD_PROCESS defined */
