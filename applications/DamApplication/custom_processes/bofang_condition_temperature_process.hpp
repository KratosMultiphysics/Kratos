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

#if !defined(KRATOS_BOFANG_CONDITION_TEMPERATURE_PROCESS )
#define  KRATOS_BOFANG_CONDITION_TEMPERATURE_PROCESS

#include <cmath>

// Project includes
#include "includes/kratos_flags.h"
#include "includes/kratos_parameters.h"
#include "processes/process.h"

// Application includes
#include "dam_application_variables.h"

namespace Kratos
{

class BofangConditionTemperatureProcess : public Process
{
    
public:

    KRATOS_CLASS_POINTER_DEFINITION(BofangConditionTemperatureProcess);
    
    typedef Table<double,double> TableType;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Constructor
    BofangConditionTemperatureProcess(ModelPart& rModelPart,
                                Parameters& rParameters
                                ) : Process(Flags()) , mrModelPart(rModelPart)
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
                "Surface_Temp"                                     : 0.0,
                "Bottom_Temp"                                      : 0.0,
                "Height_Dam"                                       : 0.0,
                "Temperature_Amplitude"                            : 0.0,
                "Day_Ambient_Temp"                                 : 1,
                "Water_level"                                      : 0.0,
                "Water_level_Table"                                : 0,
                "Outer_temp"                                       : 0.0,
                "Outer_temp_Table"                                 : 0,
                "Month"                                            : 1.0,
                "Month_Table"                                      : 0
            }  )" );
        
        // Some values need to be mandatorily prescribed since no meaningful default value exist. For this reason try accessing to them
        // So that an error is thrown if they don't exist
        rParameters["Reservoir_Bottom_Coordinate_in_Gravity_Direction"];
        rParameters["variable_name"];
        rParameters["model_part_name"];

        // Now validate agains defaults -- this also ensures no type mismatch
        rParameters.ValidateAndAssignDefaults(default_parameters);
        
        mMeshId = rParameters["mesh_id"].GetInt();
        mVariableName = rParameters["variable_name"].GetString();
        mIsFixed = rParameters["is_fixed"].GetBool();
        mGravityDirection = rParameters["Gravity_Direction"].GetString();
        mReferenceCoordinate = rParameters["Reservoir_Bottom_Coordinate_in_Gravity_Direction"].GetDouble();
        mSurfaceTemp = rParameters["Surface_Temp"].GetDouble();
        mBottomTemp = rParameters["Bottom_Temp"].GetDouble();
        mHeight = rParameters["Height_Dam"].GetDouble();
        mAmplitude = rParameters["Temperature_Amplitude"].GetDouble();
        mDay = rParameters["Day_Ambient_Temp"].GetInt();
        mWaterLevel = rParameters["Water_level"].GetDouble();
        mOuterTemp = rParameters["Outer_temp"].GetDouble();
        mMonth = rParameters["Month"].GetDouble();
        mFreq = 0.52323;

        mTimeUnitConverter = mrModelPart.GetProcessInfo()[TIME_UNIT_CONVERTER];
        mTableIdWater = rParameters["Water_level_Table"].GetInt();
        mTableIdOuter = rParameters["Outer_temp_Table"].GetInt();
        mTableIdMonth = rParameters["Month_Table"].GetInt();
        
        if(mTableIdWater != 0)
            mpTableWater = mrModelPart.pGetTable(mTableIdWater);
            
        if(mTableIdOuter != 0)
            mpTableOuter = mrModelPart.pGetTable(mTableIdOuter);
        
        if(mTableIdMonth != 0)
            mpTableMonth = mrModelPart.pGetTable(mTableIdMonth);
        
        KRATOS_CATCH("");
    }

    ///------------------------------------------------------------------------------------
    
    /// Destructor
    virtual ~BofangConditionTemperatureProcess() {}


//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void Execute()
    {
        
        KRATOS_TRY;
        
        Variable<double> var = KratosComponents< Variable<double> >::Get(mVariableName);
        const int nnodes = mrModelPart.GetMesh(mMeshId).Nodes().size();
        int direction;
        
        if( mGravityDirection == "X")
            direction = 1;
        else if( mGravityDirection == "Y")
            direction = 2;
        else
            direction = 3;
              
                
        if(nnodes != 0)
        {
            ModelPart::NodesContainerType::iterator it_begin = mrModelPart.GetMesh(mMeshId).NodesBegin();
        
            #pragma omp parallel for
            for(int i = 0; i<nnodes; i++)
            {
                ModelPart::NodesContainerType::iterator it = it_begin + i;
                
                if(mIsFixed)
                {
                    it->Fix(var);
                }
                
                double aux = (mReferenceCoordinate + mWaterLevel) - it->Coordinate(direction);
                if(aux >= 0.0)
                {
                    double aux1 = ((mBottomTemp-(mSurfaceTemp*exp(-0.04*mHeight)))/(1-(exp(-0.04*mHeight))));
                    double Temperature = (aux1+((mSurfaceTemp-aux1)*(exp(-0.04*aux)))+(mAmplitude*(exp(-0.018*aux))*(cos(mFreq*(mMonth-(mDay/30.0)-2.15+(1.30*exp(-0.085*aux)))))));
                    
                    it->FastGetSolutionStepValue(var) = Temperature;
                    
                }
                else
                    it->FastGetSolutionStepValue(var) = mOuterTemp;
            }
        }
        
        KRATOS_CATCH("");
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void ExecuteInitializeSolutionStep()
    {
        
        KRATOS_TRY;
        
        Variable<double> var = KratosComponents< Variable<double> >::Get(mVariableName);
        
        // Getting the values of table in case that it exist        
        if(mTableIdWater != 0)
        { 
            double time = mrModelPart.GetProcessInfo()[TIME];
            time = time/mTimeUnitConverter;
            mWaterLevel = mpTableWater->GetValue(time);
        }
        
        if(mTableIdOuter != 0)
        { 
            double time = mrModelPart.GetProcessInfo()[TIME];
            time = time/mTimeUnitConverter;
            mOuterTemp = mpTableOuter->GetValue(time);
        }
        
        if(mTableIdMonth != 0)
        { 
            double time = mrModelPart.GetProcessInfo()[TIME];
            time = time/mTimeUnitConverter;
            mMonth = mpTableMonth->GetValue(time);
        }

        const int nnodes = mrModelPart.GetMesh(mMeshId).Nodes().size();
        int direction;
        
        if( mGravityDirection == "X")
            direction = 1;
        else if( mGravityDirection == "Y")
            direction = 2;
        else
            direction = 3;
              
                
        if(nnodes != 0)
        {
            ModelPart::NodesContainerType::iterator it_begin = mrModelPart.GetMesh(mMeshId).NodesBegin();
        
            #pragma omp parallel for
            for(int i = 0; i<nnodes; i++)
            {
                ModelPart::NodesContainerType::iterator it = it_begin + i;
                
                if(mIsFixed)
                {
                    it->Fix(var);
                }
                
                double aux = (mReferenceCoordinate + mWaterLevel) - it->Coordinate(direction);
                if(aux >= 0.0)
                {
                    double aux1 = ((mBottomTemp-(mSurfaceTemp*exp(-0.04*mHeight)))/(1-(exp(-0.04*mHeight))));
                    double Temperature = (aux1+((mSurfaceTemp-aux1)*(exp(-0.04*aux)))+(mAmplitude*(exp(-0.018*aux))*(cos(mFreq*(mMonth-(mDay/30.0)-2.15+(1.30*exp(-0.085*aux)))))));
                    
                    it->FastGetSolutionStepValue(var) = Temperature;
                    
                }
                else
                    it->FastGetSolutionStepValue(var) = mOuterTemp;
            }
        }
        
        KRATOS_CATCH("");
    }
    
    /// Turn back information as a string.
    std::string Info() const
    {
        return "BofangConditionTemperatureProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "BofangConditionTemperatureProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const
    {
    }

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    /// Member Variables

    ModelPart& mrModelPart;
    std::size_t mMeshId;
    std::string mVariableName;
    std::string mGravityDirection;
    bool mIsFixed;
    double mReferenceCoordinate;
    double mSurfaceTemp;
    double mBottomTemp;
    double mHeight;
    double mAmplitude;
    int mDay;    
    double mMonth;
    double mWaterLevel;
    double mOuterTemp;
    double mFreq;
    double mTimeUnitConverter;
    TableType::Pointer mpTableWater;
    TableType::Pointer mpTableOuter;
    TableType::Pointer mpTableMonth;
    int mTableIdWater;
    int mTableIdOuter;
    int mTableIdMonth;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:

    /// Assignment operator.
    BofangConditionTemperatureProcess& operator=(BofangConditionTemperatureProcess const& rOther);

};//Class


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  BofangConditionTemperatureProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const BofangConditionTemperatureProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

} /* namespace Kratos.*/

#endif /* KRATOS_BOFANG_CONDITION_TEMPERATURE_PROCESS defined */

