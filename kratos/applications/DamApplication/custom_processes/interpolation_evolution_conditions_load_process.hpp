//
//   Project Name:        KratosDamApplication    $
//   Last modified by:    $Author: Lorenzo Gracia $
//   Date:                $Date:       March 2016 $
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
        
        this->Uplift_pressure(water_level);
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

    void Uplift_pressure(result_type const& water_level)
    {
        // We have to pick up the provided values by the mesh
        std::string direction = (*(mr_model_part.pGetMesh(4)))[GRAVITY_DIRECTION];
        std::string uplift_direction = (*(mr_model_part.pGetMesh(4)))[UPLIFT_DIRECTION];
        const double& coordinate_base = (*(mr_model_part.pGetMesh(4)))[COORDINATE_BASE_DAM];
        const double& coordinate_base_uplift = (*(mr_model_part.pGetMesh(4)))[COORDINATE_BASE_DAM_UPLIFT];
        const double& base_dam = (*(mr_model_part.pGetMesh(4)))[BASE_OF_DAM];
        const double& spe_weight =  (*(mr_model_part.pGetMesh(4)))[SPECIFIC_WEIGHT];
        
        double ref_coord;
        
        for(ModelPart::NodeIterator i = mr_model_part.pGetMesh(4)->NodesBegin(); i <mr_model_part.pGetMesh(4)->NodesEnd(); i++)
        {
            double& uplift_pressure  = (i)->FastGetSolutionStepValue(NORMAL_CONTACT_STRESS);
            
            if(direction=="X")
            {
                ref_coord = coordinate_base + water_level;
                if(uplift_direction=="Y")
                {
                    uplift_pressure = (spe_weight*(ref_coord- (i)->X()))*(1.0-((1.0/base_dam)*(fabs( (i)->Y() - coordinate_base_uplift))));
                        if(uplift_pressure<0.0)
                            uplift_pressure=0.0;
                }            
                else if(uplift_direction=="Z")
                {
                    uplift_pressure = (spe_weight*(ref_coord- (i)->X()))*(1.0-((1.0/base_dam)*(fabs( (i)->Z() - coordinate_base_uplift))));
                        if(uplift_pressure<0.0)
                            uplift_pressure=0.0;
                }            
            }
            else if(direction=="Y")
            {
                ref_coord = coordinate_base + water_level;
                if(uplift_direction=="X")
                {
                    uplift_pressure = (spe_weight*(ref_coord- (i)->Y()))*(1.0-((1.0/base_dam)*(fabs( (i)->X() - coordinate_base_uplift))));
                        if(uplift_pressure<0.0)
                            uplift_pressure=0.0;
                }            
                else if(uplift_direction=="Z")
                {
                    uplift_pressure = (spe_weight*(ref_coord- (i)->Y()))*(1.0-((1.0/base_dam)*(fabs( (i)->Z() - coordinate_base_uplift))));
                        if(uplift_pressure<0.0)
                            uplift_pressure=0.0;
                }
            }
            else if(direction=="Z")
            {
                ref_coord = coordinate_base + water_level;
                if(uplift_direction=="X")
                {
                    uplift_pressure = (spe_weight*(ref_coord- (i)->Z()))*(1.0-((1.0/base_dam)*(fabs( (i)->X() - coordinate_base_uplift))));
                        if(uplift_pressure<0.0)
                            uplift_pressure=0.0;
                }            
                else if(uplift_direction=="Y")
                {
                    uplift_pressure = (spe_weight*(ref_coord- (i)->Z()))*(1.0-((1.0/base_dam)*(fabs( (i)->Y() - coordinate_base_uplift))));
                        if(uplift_pressure<0.0)
                            uplift_pressure=0.0;
                }
            }
        }
    }
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------    

};//Class

} /* namespace Kratos.*/

#endif /* KRATOS_INTERPOLATION_EVOLUTION_CONDITIONS_LOAD_PROCESS defined */
