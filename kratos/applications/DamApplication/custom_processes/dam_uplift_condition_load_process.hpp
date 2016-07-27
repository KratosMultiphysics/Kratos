//
//   Project Name:        KratosDamApplication    $
//   Last modified by:    $Author: Lorenzo Gracia $
//   Date:                $Date:        July 2016 $
//   Revision:            $Revision:          0.0 $
//

#if !defined(KRATOS_DAM_UPLIFT_CONDITION_LOAD_PROCESS )
#define  KRATOS_DAM_UPLIFT_CONDITION_LOAD_PROCESS

#include <cmath>

#include "includes/kratos_flags.h"
#include "includes/kratos_parameters.h"
#include "processes/process.h"

#include "dam_application_variables.h"

namespace Kratos
{

class DamUpliftConditionLoadProcess : public Process
{
    
public:

    KRATOS_CLASS_POINTER_DEFINITION(DamUpliftConditionLoadProcess);

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Constructor
    DamUpliftConditionLoadProcess(ModelPart& model_part,
                                Parameters rParameters
                                ) : Process(Flags()) , mr_model_part(model_part)
    {
        KRATOS_TRY
			 
        //only include validation with c++11 since raw_literals do not exist in c++03
        Parameters default_parameters( R"(
            {
                "model_part_name":"PLEASE_CHOOSE_MODEL_PART_NAME",
                "mesh_id": 0,
                "variable_name": "PLEASE_PRESCRIBE_VARIABLE_NAME",
                "is_fixed"                                              : false,
                "Modify"                                                : true,
                "Gravity_Direction"                                     : "Y",
                "Uplift_Direction"                                      : "X",
                "Reservoir_Bottom_Coordinate_in_Gravity_Direction"      : 0.0,
                "Upstream_Coordinate_at_Base_Dam_in_Uplift_Direction"   : 0.0,
                "Base_of_dam"                                           : 0.0,
                "Spe_weight"                                            : 0.0,
                "Water_level"                                           : 10,
                "Drains"                                                : false,
                "Height_drain"                                          : 0.0,
                "Position_drain"                                        : 0.0,
                "Effectiveness"                                         : 0.0
            }  )" );
            
        
        // Some values need to be mandatorily prescribed since no meaningful default value exist. For this reason try accessing to them
        // So that an error is thrown if they don't exist
        rParameters["Reservoir_Bottom_Coordinate_in_Gravity_Direction"];
        rParameters["Upstream_Coordinate_at_Base_Dam_in_Uplift_Direction"];
        rParameters["variable_name"];
        rParameters["model_part_name"];

        // Now validate agains defaults -- this also ensures no type mismatch
        rParameters.ValidateAndAssignDefaults(default_parameters);
        
        mmesh_id = rParameters["mesh_id"].GetInt();
        mvariable_name = rParameters["variable_name"].GetString();
        mis_fixed = rParameters["is_fixed"].GetBool();
        mgravity_direction = rParameters["Gravity_Direction"].GetString();
        muplift_direction = rParameters["Uplift_Direction"].GetString();
        mreference_coordinate = rParameters["Reservoir_Bottom_Coordinate_in_Gravity_Direction"].GetDouble();
        mreference_coordinate_uplift = rParameters["Upstream_Coordinate_at_Base_Dam_in_Uplift_Direction"].GetDouble();
        mspecific = rParameters["Spe_weight"].GetDouble();
        mbase_dam = rParameters["Base_of_dam"].GetDouble();
        
        // Drains
        mdrain = rParameters["Drains"].GetBool();
        mheight_drain = rParameters["Height_drain"].GetDouble();
        mlength_drain = rParameters["Position_drain"].GetDouble();
        meffectiveness_drain = rParameters["Effectiveness"].GetDouble();
        
        // TODO: PARAMETERS MUST BE GOT FROM THE TABLE
        mwater_level = rParameters["Water_level"].GetDouble();

        KRATOS_CATCH("");
    }

    ///------------------------------------------------------------------------------------
    
    /// Destructor
    virtual ~DamUpliftConditionLoadProcess() {}


//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void Execute()
    {
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void ExecuteInitialize()
    {
        
        KRATOS_TRY;
        
        Variable<double> var = KratosComponents< Variable<double> >::Get(mvariable_name);
        
        const int nnodes = mr_model_part.GetMesh(mmesh_id).Nodes().size();
        
        int direction, up_direction;
        
        if( mgravity_direction == "X")
            direction = 1;
        else if( mgravity_direction == "Y")
            direction = 2;
        else if( mgravity_direction == "Z")
            direction = 3;
            
        if( muplift_direction == "X")
            up_direction = 1;
        else if( muplift_direction == "Y")
            up_direction = 2;
        else if( muplift_direction == "Z")
            up_direction = 3;
              
                
        if(nnodes != 0)
        {
            ModelPart::NodesContainerType::iterator it_begin = mr_model_part.GetMesh(mmesh_id).NodesBegin();
            
            double ref_coord = mreference_coordinate + mwater_level;
            
            if( mdrain == true)
            {
				double coefficient_effectiveness = 1.0 - meffectiveness_drain;
				double aux_drain = coefficient_effectiveness *(mwater_level - mheight_drain)* ((mbase_dam-mlength_drain)/mbase_dam) + mheight_drain;
				
				#pragma omp parallel for
				for(int i = 0; i<nnodes; i++)
				{
					ModelPart::NodesContainerType::iterator it = it_begin + i;

					if(mis_fixed)
					{
						it->Fix(var);
					}
					
					if ( (it->Coordinate(up_direction) +0.000001) <= (mreference_coordinate_uplift + mlength_drain ))
					{						
						muplift_pressure = (mspecific*((ref_coord-aux_drain)- (it->Coordinate(direction))))*(1.0-((1.0/(mlength_drain))*(fabs( (it->Coordinate(up_direction)) - mreference_coordinate_uplift)))) + mspecific*aux_drain;
					}
					else
					{
						muplift_pressure = (mspecific*((mreference_coordinate+aux_drain)- (it->Coordinate(direction))))*(1.0-((1.0/(mbase_dam - mlength_drain))*(fabs( (it->Coordinate(up_direction)) - (mreference_coordinate_uplift+mlength_drain)))));
					}
					
					if(muplift_pressure<0.0)
					{
						it->FastGetSolutionStepValue(var)=0.0;
					}
					else
					{
						it->FastGetSolutionStepValue(var) = muplift_pressure;
					}
				}
					
			}
			else
			{
				#pragma omp parallel for
				for(int i = 0; i<nnodes; i++)
				{
					ModelPart::NodesContainerType::iterator it = it_begin + i;

					if(mis_fixed)
					{
						it->Fix(var);
					}
				
					muplift_pressure = (mspecific*(ref_coord- (it->Coordinate(direction))))*(1.0-((1.0/mbase_dam)*(fabs( (it->Coordinate(up_direction)) - mreference_coordinate_uplift))));
					
					if(muplift_pressure<0.0)
					{
						it->FastGetSolutionStepValue(var)=0.0;
					}
					else
					{
						KRATOS_WATCH(muplift_pressure)
						it->FastGetSolutionStepValue(var) = muplift_pressure;
					}
				}
			}            
        }
        
        KRATOS_CATCH("");
    }
    
    /// Turn back information as a string.
    std::string Info() const
    {
        return "DamUpliftConditionLoadProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "DamUpliftConditionLoadProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const
    {
    }

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    /// Member Variables

    ModelPart& mr_model_part;
    std::size_t mmesh_id;
    std::string mvariable_name;
    std::string mgravity_direction;
    std::string muplift_direction;
    bool mis_fixed;
    double mreference_coordinate;
    double mreference_coordinate_uplift;
    double mspecific;
    double mbase_dam;
    double mwater_level;
    bool mdrain;
    double mheight_drain;
    double mlength_drain;
    double meffectiveness_drain;
    double muplift_pressure;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:

    /// Assignment operator.
    DamUpliftConditionLoadProcess& operator=(DamUpliftConditionLoadProcess const& rOther);

};//Class


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  DamUpliftConditionLoadProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const DamUpliftConditionLoadProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

} /* namespace Kratos.*/

#endif /* KRATOS_DAM_UPLIFT_CONDITION_LOAD_PROCESS defined */
