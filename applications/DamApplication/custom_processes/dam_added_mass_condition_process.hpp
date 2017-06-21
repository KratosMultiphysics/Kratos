//
//   Project Name:        KratosDamApplication    $
//   Last modified by:    $Author: Lorenzo Gracia $
//   Date:                $Date:          May 2017 $
//   Revision:            $Revision:          0.0 $
//

#if !defined(KRATOS_ADDED_MASS_CONDITION_PROCESS )
#define  KRATOS_ADDED_MASS_CONDITION_PROCESS

#include <cmath>

#include "includes/kratos_flags.h"
#include "includes/kratos_parameters.h"
#include "processes/process.h"

#include "dam_application_variables.h"

namespace Kratos
{

class DamAddedMassConditionProcess : public Process
{
    
public:

    KRATOS_CLASS_POINTER_DEFINITION(DamAddedMassConditionProcess);
    
    typedef Table<double,double> TableType;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Constructor
    DamAddedMassConditionProcess(ModelPart& model_part,
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
                "Modify"                                                : true,
                "Gravity_Direction"                                     : "Y",
                "Reservoir_Bottom_Coordinate_in_Gravity_Direction"      : 0.0,
                "Spe_weight"                                            : 0.0,
                "Water_level"                                           : 0.0
            }  )" );
            
        
        // Some values need to be mandatorily prescribed since no meaningful default value exist. For this reason try accessing to them
        // So that an error is thrown if they don't exist
        rParameters["Reservoir_Bottom_Coordinate_in_Gravity_Direction"];
        rParameters["variable_name"];
        rParameters["model_part_name"];

        // Now validate agains defaults -- this also ensures no type mismatch
        rParameters.ValidateAndAssignDefaults(default_parameters);
        
        mmesh_id = rParameters["mesh_id"].GetInt();
        mvariable_name = rParameters["variable_name"].GetString();
        mgravity_direction = rParameters["Gravity_Direction"].GetString();
        mreference_coordinate = rParameters["Reservoir_Bottom_Coordinate_in_Gravity_Direction"].GetDouble();
        mspecific = rParameters["Spe_weight"].GetDouble();
        mwater_level = rParameters["Water_level"].GetDouble();

        KRATOS_CATCH("");
    }

    ///------------------------------------------------------------------------------------
    
    /// Destructor
    virtual ~DamAddedMassConditionProcess() {}


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
        int direction;
        double added_mass;
        
        if( mgravity_direction == "X")
            direction = 1;
        else if( mgravity_direction == "Y")
            direction = 2;
        else
            direction = 3;
        
		double ref_coord = mreference_coordinate + mwater_level;
                  
        if(nnodes != 0)
        {
            ModelPart::NodesContainerType::iterator it_begin = mr_model_part.GetMesh(mmesh_id).NodesBegin();

            #pragma omp parallel for
            for(int i = 0; i<nnodes; i++)
            {
                ModelPart::NodesContainerType::iterator it = it_begin + i;

                double y_water =  ref_coord- (it->Coordinate(direction));
                
                if (y_water<0.0)
                {
                    y_water=0.0;
                }
                
                added_mass = 0.875*mspecific*sqrt(y_water*mwater_level);

                it->FastGetSolutionStepValue(var) = added_mass;

                if(added_mass>0.0)
                {
                    it->FastGetSolutionStepValue(var) = added_mass;
                }
                else
                {
                    it->FastGetSolutionStepValue(var)=0.0;
                }
            }            
        }
        
        KRATOS_CATCH("");
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void ExecuteInitializeSolutionStep()
    {
         KRATOS_TRY;
        
        Variable<double> var = KratosComponents< Variable<double> >::Get(mvariable_name);
        const int nnodes = mr_model_part.GetMesh(mmesh_id).Nodes().size();
        int direction;
        double added_mass;
        
        if( mgravity_direction == "X")
            direction = 1;
        else if( mgravity_direction == "Y")
            direction = 2;
        else
            direction = 3;
        
		double ref_coord = mreference_coordinate + mwater_level;
                  
        if(nnodes != 0)
        {
            ModelPart::NodesContainerType::iterator it_begin = mr_model_part.GetMesh(mmesh_id).NodesBegin();

            #pragma omp parallel for
            for(int i = 0; i<nnodes; i++)
            {
                ModelPart::NodesContainerType::iterator it = it_begin + i;

                double y_water =  ref_coord- (it->Coordinate(direction));
                
                if (y_water<0.0)
                {
                    y_water=0.0;
                }
                
                added_mass = 0.875*mspecific*sqrt(y_water*mwater_level);

                it->FastGetSolutionStepValue(var) = added_mass;

                if(added_mass>0.0)
                {
                    it->FastGetSolutionStepValue(var) = added_mass;
                }
                else
                {
                    it->FastGetSolutionStepValue(var)=0.0;
                }
            }            
        }
        
        KRATOS_CATCH("");
    }
    
    /// Turn back information as a string.
    std::string Info() const
    {
        return "DamAddedMassConditionProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "DamAddedMassConditionProcess";
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
    double mreference_coordinate;
    double mspecific;
    double mwater_level;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:

    /// Assignment operator.
    DamAddedMassConditionProcess& operator=(DamAddedMassConditionProcess const& rOther);

};//Class


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  DamAddedMassConditionProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const DamAddedMassConditionProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

} /* namespace Kratos.*/

#endif /* KRATOS_ADDED_MASS_CONDITION_PROCESS defined */
