//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Lorenzo Gracia
//
//

#if !defined(KRATOS_DAM_RESERVOIR_CONSTANT_TEMPERATURE_PROCESS )
#define  KRATOS_DAM_RESERVOIR_CONSTANT_TEMPERATURE_PROCESS

#include <cmath>

// Project includes
#include "includes/kratos_flags.h"
#include "includes/kratos_parameters.h"
#include "processes/process.h"

// Application include
#include "dam_application_variables.h"

namespace Kratos
{

class DamReservoirConstantTemperatureProcess : public Process
{
    
public:

    KRATOS_CLASS_POINTER_DEFINITION(DamReservoirConstantTemperatureProcess);
    
    typedef Table<double,double> TableType;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Constructor
    DamReservoirConstantTemperatureProcess(ModelPart& model_part,
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
                "is_fixed"                                         : false,
                "Gravity_Direction"                                : "Y",
                "Reservoir_Bottom_Coordinate_in_Gravity_Direction" : 0.0,
                "Water_temp"                                       : 0.0,
                "Water_temp_Table"                                 : 0,                
                "Water_level"                                      : 0.0,
                "Water_level_Table"                                : 0,
                "Outer_temp"                                       : 0.0,
                "Outer_temp_Table"                                 : 0
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
        mis_fixed = rParameters["is_fixed"].GetBool();
        mgravity_direction = rParameters["Gravity_Direction"].GetString();
        mreference_coordinate = rParameters["Reservoir_Bottom_Coordinate_in_Gravity_Direction"].GetDouble();
        mwater_temp = rParameters["Water_temp"].GetDouble();
        mwater_level = rParameters["Water_level"].GetDouble();
        mouter_temp = rParameters["Outer_temp"].GetDouble();

        mtime_unit_converter = mr_model_part.GetProcessInfo()[TIME_UNIT_CONVERTER];
        mTableIdWaterTemp = rParameters["Water_temp_Table"].GetInt();        
        mTableIdWater = rParameters["Water_level_Table"].GetInt();
        mTableIdOuter = rParameters["Outer_temp_Table"].GetInt();

        if(mTableIdWaterTemp != 0)
            mpTableWaterTemp = model_part.pGetTable(mTableIdWaterTemp);

        if(mTableIdWater != 0)
            mpTableWater = model_part.pGetTable(mTableIdWater);
            
        if(mTableIdOuter != 0)
            mpTableOuter = model_part.pGetTable(mTableIdOuter);

        KRATOS_CATCH("");
    }

    ///------------------------------------------------------------------------------------
    
    /// Destructor
    virtual ~DamReservoirConstantTemperatureProcess() {}


//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void Execute()
    {
        
        KRATOS_TRY;
        
        Variable<double> var = KratosComponents< Variable<double> >::Get(mvariable_name);
        const int nnodes = mr_model_part.GetMesh(mmesh_id).Nodes().size();
        int direction;
        
        if( mgravity_direction == "X")
            direction = 1;
        else if( mgravity_direction == "Y")
            direction = 2;
        else
            direction = 3;
              
                
        if(nnodes != 0)
        {
            ModelPart::NodesContainerType::iterator it_begin = mr_model_part.GetMesh(mmesh_id).NodesBegin();
        
            #pragma omp parallel for
            for(int i = 0; i<nnodes; i++)
            {
                ModelPart::NodesContainerType::iterator it = it_begin + i;
                
                if(mis_fixed)
                {
                    it->Fix(var);
                }
                
                double aux = (mreference_coordinate + mwater_level) - it->Coordinate(direction);
                if(aux >= 0.0)
                {                   
                    it->FastGetSolutionStepValue(var) = mwater_temp;
                    
                }
                else
                    it->FastGetSolutionStepValue(var) = mouter_temp;
            }
        }
        
        KRATOS_CATCH("");
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void ExecuteInitializeSolutionStep()
    {
        
        KRATOS_TRY;
        
        Variable<double> var = KratosComponents< Variable<double> >::Get(mvariable_name);
        
        // Getting the values of table in case that it exist

        if(mTableIdWaterTemp != 0)
        { 
            double time = mr_model_part.GetProcessInfo()[TIME];
            time = time/mtime_unit_converter;
            mwater_temp = mpTableWaterTemp->GetValue(time);
        }

        if(mTableIdWater != 0)
        { 
            double time = mr_model_part.GetProcessInfo()[TIME];
            time = time/mtime_unit_converter;
            mwater_level = mpTableWater->GetValue(time);
        }
        
        if(mTableIdOuter != 0)
        { 
            double time = mr_model_part.GetProcessInfo()[TIME];
            time = time/mtime_unit_converter;
            mouter_temp = mpTableOuter->GetValue(time);
        }
        
        const int nnodes = mr_model_part.GetMesh(mmesh_id).Nodes().size();
        int direction;
        
        if( mgravity_direction == "X")
            direction = 1;
        else if( mgravity_direction == "Y")
            direction = 2;
        else
            direction = 3;
              
                
        if(nnodes != 0)
        {
            ModelPart::NodesContainerType::iterator it_begin = mr_model_part.GetMesh(mmesh_id).NodesBegin();
        
            #pragma omp parallel for
            for(int i = 0; i<nnodes; i++)
            {
                ModelPart::NodesContainerType::iterator it = it_begin + i;
                
                if(mis_fixed)
                {
                    it->Fix(var);
                }
                
                double aux = (mreference_coordinate + mwater_level) - it->Coordinate(direction);
                if(aux >= 0.0)
                { 
                    it->FastGetSolutionStepValue(var) = mwater_temp;
                    
                }
                else
                    it->FastGetSolutionStepValue(var) = mouter_temp;
            }
        }
        
        KRATOS_CATCH("");
    }
    
    /// Turn back information as a string.
    std::string Info() const
    {
        return "DamReservoirConstantTemperatureProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "DamReservoirConstantTemperatureProcess";
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
    double mwater_temp;
    double mwater_level;
    double mouter_temp;
    double mtime_unit_converter;
    TableType::Pointer mpTableWaterTemp;
    TableType::Pointer mpTableWater;
    TableType::Pointer mpTableOuter;
    int mTableIdWaterTemp;
    int mTableIdWater;
    int mTableIdOuter;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:

    /// Assignment operator.
    DamReservoirConstantTemperatureProcess& operator=(DamReservoirConstantTemperatureProcess const& rOther);

};//Class


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  DamReservoirConstantTemperatureProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const DamReservoirConstantTemperatureProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

} /* namespace Kratos.*/

#endif /* KRATOS_DAM_RESERVOIR_CONSTANT_TEMPERATURE_PROCESS defined */
