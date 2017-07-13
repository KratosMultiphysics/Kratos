//
//   Project Name:        KratosDamApplication    $
//   Last modified by:    $Author: Lorenzo Gracia $
//   Date:                $Date:     October 2016 $
//   Revision:            $Revision:          0.0 $
//

#if !defined(KRATOS_DAM_UPLIFT_CIRCULAR_CONDITION_LOAD_PROCESS )
#define  KRATOS_DAM_UPLIFT_CIRCULAR_CONDITION_LOAD_PROCESS

#include <cmath>

#include "includes/kratos_flags.h"
#include "includes/kratos_parameters.h"
#include "processes/process.h"

#include "dam_application_variables.h"

namespace Kratos
{

class DamUpliftCircularConditionLoadProcess : public Process
{
    
public:

    KRATOS_CLASS_POINTER_DEFINITION(DamUpliftCircularConditionLoadProcess);

    typedef Table<double,double> TableType;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Constructor
    DamUpliftCircularConditionLoadProcess(ModelPart& model_part,
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
                "Reservoir_Bottom_Coordinate_in_Gravity_Direction"      : 0.0,
                "Upstream_Coordinate_first_bracket"                     : [0.0,0.0,0.0],
                "Downstream_Coordinate_first_bracket"                   : [0.0,0.0,0.0],
                "Focus"                                                 : [0.0,0.0,0.0],
                "Spe_weight"                                            : 10000,
                "Water_level"                                           : 0.0,
                "Drains"                                                : false,
                "Height_drain"                                          : 0.0,
                "Distance"                                              : 0.0,
                "Effectiveness"                                         : 0.0,
                "table"                                                 : 0 
            }  )" );
            
        
        // Some values need to be mandatorily prescribed since no meaningful default value exist. For this reason try accessing to them
        // So that an error is thrown if they don't exist
        rParameters["Focus"];


        // Now validate agains defaults -- this also ensures no type mismatch
        rParameters.ValidateAndAssignDefaults(default_parameters);
        
        mmesh_id = rParameters["mesh_id"].GetInt();
        mvariable_name = rParameters["variable_name"].GetString();
        mis_fixed = rParameters["is_fixed"].GetBool();
        mgravity_direction = rParameters["Gravity_Direction"].GetString();
        mreference_coordinate = rParameters["Reservoir_Bottom_Coordinate_in_Gravity_Direction"].GetDouble();
        mspecific = rParameters["Spe_weight"].GetDouble();
        
        // Getting the values of the coordinates at first bracket (reference value)
        mupstream.resize(3,false);
        mupstream[0] = rParameters["Upstream_Coordinate_first_bracket"][0].GetDouble();
        mupstream[1] = rParameters["Upstream_Coordinate_first_bracket"][1].GetDouble();
        mupstream[2] = rParameters["Upstream_Coordinate_first_bracket"][2].GetDouble();
        
        mdownstream.resize(3,false);
        mdownstream[0] = rParameters["Downstream_Coordinate_first_bracket"][0].GetDouble();
        mdownstream[1] = rParameters["Downstream_Coordinate_first_bracket"][1].GetDouble();
        mdownstream[2] = rParameters["Downstream_Coordinate_first_bracket"][2].GetDouble();
        
        // Getting the coordinates of the focus (reference value)
        mfocus.resize(3,false);
        mfocus[0] = rParameters["Focus"][0].GetDouble();
        mfocus[1] = rParameters["Focus"][1].GetDouble();
        mfocus[2] = rParameters["Focus"][2].GetDouble();        
        
        // Drains
        mdrain = rParameters["Drains"].GetBool();
        mheight_drain = rParameters["Height_drain"].GetDouble();
        mdistance_drain = rParameters["Distance"].GetDouble();
        meffectiveness_drain = rParameters["Effectiveness"].GetDouble();
        mwater_level = rParameters["Water_level"].GetInt();

        mtime_unit_converter = mr_model_part.GetProcessInfo()[TIME_UNIT_CONVERTER];
        mTableId = rParameters["table"].GetInt();
        
        if(mTableId != 0)
            mpTable = model_part.pGetTable(mTableId);

        KRATOS_CATCH("");
    }

    ///------------------------------------------------------------------------------------
    
    /// Destructor
    virtual ~DamUpliftCircularConditionLoadProcess() {}


//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void Execute()
    {
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void ExecuteInitialize()
    {
        
        KRATOS_TRY;
        
        //Defining necessary variables
        Variable<double> var = KratosComponents< Variable<double> >::Get(mvariable_name);
        const int nnodes = mr_model_part.GetMesh(mmesh_id).Nodes().size();
        array_1d<double,3> auxiliar_vector;
        
        // Gravity direction for computing the hydrostatic pressure
        int direction;
        int radius_comp_1;
        int radius_comp_2;
        
        if( mgravity_direction == "X")
        {
            direction = 1;
            radius_comp_1 = 1;
            radius_comp_2 = 2;
        }
        else if( mgravity_direction == "Y")
        {
            direction = 2;
            radius_comp_1 = 0;
            radius_comp_2 = 2;
        }
        else
        {
            direction = 3;
            radius_comp_1 = 0;
            radius_comp_2 = 1;
        }
        
        // Computation of the angle in radians for each bracket
        //double tetha = (mangle*2*pi)/(mnum_brackets*360);
        
        // Computation of radius for Upstream and Downstream
        double up_radius = norm_2(mfocus - mupstream);
        double down_radius = norm_2(mfocus - mdownstream);
        double width_dam = up_radius - down_radius;
                
        if(nnodes != 0)
        {
            ModelPart::NodesContainerType::iterator it_begin = mr_model_part.GetMesh(mmesh_id).NodesBegin();
            
            double ref_coord = mreference_coordinate + mwater_level;
            
            if( mdrain == true)
            {
				double coefficient_effectiveness = 1.0 - meffectiveness_drain;
				double aux_drain = coefficient_effectiveness *(mwater_level - mheight_drain)* ((width_dam-mdistance_drain)/width_dam) + mheight_drain;

				#pragma omp parallel for
				for(int i = 0; i<nnodes; i++)
				{
					ModelPart::NodesContainerType::iterator it = it_begin + i;
                    
                    auxiliar_vector.resize(3,false);                    
                    auxiliar_vector[0] = mfocus[0] - (it->Coordinate(1));
                    auxiliar_vector[1] = mfocus[1] - (it->Coordinate(2));
                    auxiliar_vector[2] = mfocus[2] - (it->Coordinate(3));
                    
                    //// Computing the new coordinates                                        
                    double current_radius = sqrt(auxiliar_vector[radius_comp_1]*auxiliar_vector[radius_comp_1] + auxiliar_vector[radius_comp_2]*auxiliar_vector[radius_comp_2]);

					if(mis_fixed)
					{
						it->Fix(var);
					}
		
                    //// We compute the first part of the uplift law 
                    muplift_pressure = mspecific*((ref_coord -aux_drain) - (it->Coordinate(direction)))*(1.0 - ((1.0/mdistance_drain)*(fabs(current_radius-up_radius)))) + (mspecific * aux_drain); 
                                       
                    //// If uplift pressure is greater than the limit we compute the second part and we update the value
                    if(muplift_pressure <= mspecific*aux_drain)
                    {
                        muplift_pressure = (mspecific*((mreference_coordinate + aux_drain)-(it->Coordinate(direction))))*(1.0 - ((1.0/(width_dam - mdistance_drain))*(fabs(current_radius-(up_radius-mdistance_drain)))));
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
                    
                    auxiliar_vector.resize(3,false);                    
                    auxiliar_vector[0] = mfocus[0] - (it->Coordinate(1));
                    auxiliar_vector[1] = mfocus[1] - (it->Coordinate(2));
                    auxiliar_vector[2] = mfocus[2] - (it->Coordinate(3));
                    
                    // Computing the current distance to the focus.                    
                    double current_radius = sqrt(auxiliar_vector[radius_comp_1]*auxiliar_vector[radius_comp_1] + auxiliar_vector[radius_comp_2]*auxiliar_vector[radius_comp_2]);
                    
					if(mis_fixed)
					{
						it->Fix(var);
					}
				
					muplift_pressure = mspecific*(ref_coord - (it->Coordinate(direction)))*(1.0 - (1.0/width_dam)*(fabs(current_radius-up_radius)));
                    
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
        }
        
        KRATOS_CATCH("");
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void ExecuteInitializeSolutionStep()
    {
        
        KRATOS_TRY;
        
        //Defining necessary variables
        Variable<double> var = KratosComponents< Variable<double> >::Get(mvariable_name);
        const int nnodes = mr_model_part.GetMesh(mmesh_id).Nodes().size();
        array_1d<double,3> auxiliar_vector;
        
        // Getting the values of table in case that it exist        
        if(mTableId != 0 )
        { 
            double time = mr_model_part.GetProcessInfo()[TIME];
            time = time/mtime_unit_converter;
            mwater_level = mpTable->GetValue(time);
        }


        // Gravity direction for computing the hydrostatic pressure
        int direction;
        int radius_comp_1;
        int radius_comp_2;
        
        if( mgravity_direction == "X")
        {
            direction = 1;
            radius_comp_1 = 1;
            radius_comp_2 = 2;
        }
        else if( mgravity_direction == "Y")
        {
            direction = 2;
            radius_comp_1 = 0;
            radius_comp_2 = 2;
        }
        else
        {
            direction = 3;
            radius_comp_1 = 0;
            radius_comp_2 = 1;
        }
        
        // Computation of the angle in radians for each bracket
        //double tetha = (mangle*2*pi)/(mnum_brackets*360);
        
        // Computation of radius for Upstream and Downstream
        double up_radius = norm_2(mfocus - mupstream);
        double down_radius = norm_2(mfocus - mdownstream);
        double width_dam = up_radius - down_radius;
                
        if(nnodes != 0)
        {
            ModelPart::NodesContainerType::iterator it_begin = mr_model_part.GetMesh(mmesh_id).NodesBegin();
            
            double ref_coord = mreference_coordinate + mwater_level;
            
            if( mdrain == true)
            {
				double coefficient_effectiveness = 1.0 - meffectiveness_drain;
				double aux_drain = coefficient_effectiveness *(mwater_level - mheight_drain)* ((width_dam-mdistance_drain)/width_dam) + mheight_drain;

				#pragma omp parallel for
				for(int i = 0; i<nnodes; i++)
				{
					ModelPart::NodesContainerType::iterator it = it_begin + i;
                    
                    auxiliar_vector.resize(3,false);                    
                    auxiliar_vector[0] = mfocus[0] - (it->Coordinate(1));
                    auxiliar_vector[1] = mfocus[1] - (it->Coordinate(2));
                    auxiliar_vector[2] = mfocus[2] - (it->Coordinate(3));
                    
                    //// Computing the new coordinates                                        
                    double current_radius = sqrt(auxiliar_vector[radius_comp_1]*auxiliar_vector[radius_comp_1] + auxiliar_vector[radius_comp_2]*auxiliar_vector[radius_comp_2]);

					if(mis_fixed)
					{
						it->Fix(var);
					}
		
                    //// We compute the first part of the uplift law 
                    muplift_pressure = mspecific*((ref_coord -aux_drain) - (it->Coordinate(direction)))*(1.0 - ((1.0/mdistance_drain)*(fabs(current_radius-up_radius)))) + (mspecific * aux_drain); 
                                       
                    //// If uplift pressure is greater than the limit we compute the second part and we update the value
                    if(muplift_pressure <= mspecific*aux_drain)
                    {
                        muplift_pressure = (mspecific*((mreference_coordinate + aux_drain)-(it->Coordinate(direction))))*(1.0 - ((1.0/(width_dam - mdistance_drain))*(fabs(current_radius-(up_radius-mdistance_drain)))));
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
                    
                    auxiliar_vector.resize(3,false);                    
                    auxiliar_vector[0] = mfocus[0] - (it->Coordinate(1));
                    auxiliar_vector[1] = mfocus[1] - (it->Coordinate(2));
                    auxiliar_vector[2] = mfocus[2] - (it->Coordinate(3));
                    
                    // Computing the current distance to the focus.                    
                    double current_radius = sqrt(auxiliar_vector[radius_comp_1]*auxiliar_vector[radius_comp_1] + auxiliar_vector[radius_comp_2]*auxiliar_vector[radius_comp_2]);
                    
					if(mis_fixed)
					{
						it->Fix(var);
					}
				
					muplift_pressure = mspecific*(ref_coord - (it->Coordinate(direction)))*(1.0 - (1.0/width_dam)*(fabs(current_radius-up_radius)));
                    
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
        }
        
        KRATOS_CATCH("");
    }
    
    /// Turn back information as a string.
    std::string Info() const
    {
        return "DamUpliftCircularConditionLoadProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "DamUpliftCircularConditionLoadProcess";
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
    bool mis_fixed;
    double mreference_coordinate;
    double mspecific;
    double mbase_dam;
    double mwater_level;
    bool mdrain;
    double mheight_drain;
    double mdistance_drain;
    double meffectiveness_drain;
    double muplift_pressure;
    Vector mupstream;
    Vector mdownstream;
    Vector mfocus;
    double mtime_unit_converter;
    TableType::Pointer mpTable;
    int mTableId;  
    
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:

    /// Assignment operator.
    DamUpliftCircularConditionLoadProcess& operator=(DamUpliftCircularConditionLoadProcess const& rOther);

};//Class


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  DamUpliftCircularConditionLoadProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const DamUpliftCircularConditionLoadProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

} /* namespace Kratos.*/

#endif /* KRATOS_DAM_UPLIFT_CIRCULAR_CONDITION_LOAD_PROCESS defined */
