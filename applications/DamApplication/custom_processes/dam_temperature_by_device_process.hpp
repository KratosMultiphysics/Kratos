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

#if !defined(KRATOS_DAM_TEMPERATURE_BY_DEVICE_PROCESS )
#define  KRATOS_DAM_TEMPERATURE_BY_DEVICE_PROCESS

#include <cmath>

// Project includes
#include "includes/kratos_flags.h"
#include "includes/kratos_parameters.h"
#include "processes/process.h"

// Application include
#include "dam_application_variables.h"

namespace Kratos
{

class DamTemperaturebyDeviceProcess : public Process
{
    
public:

    KRATOS_CLASS_POINTER_DEFINITION(DamTemperaturebyDeviceProcess);

    typedef Table<double,double> TableType;  
    
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Constructor
    DamTemperaturebyDeviceProcess(ModelPart& model_part,
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
                "is_fixed"          : false,
                "value"             : 0.0,
                "table"             : 0,
                "position"          : [0.0,0.0,0.0]
            }  )" );
        
        // Some values need to be mandatorily prescribed since no meaningful default value exist. For this reason try accessing to them
        // So that an error is thrown if they don't exist
        rParameters["value"];
        rParameters["variable_name"];
        rParameters["model_part_name"];

        // Now validate agains defaults -- this also ensures no type mismatch
        rParameters.ValidateAndAssignDefaults(default_parameters);

        mmesh_id = rParameters["mesh_id"].GetInt();
        mvariable_name = rParameters["variable_name"].GetString();
        mis_fixed = rParameters["is_fixed"].GetBool();
        mvalue = rParameters["value"].GetDouble();
        // Getting the values of the device coordinates
        mdevice_coordinates.resize(3,false);
        mdevice_coordinates[0] = rParameters["position"][0].GetDouble();
        mdevice_coordinates[1] = rParameters["position"][1].GetDouble();
        mdevice_coordinates[2] = rParameters["position"][2].GetDouble();
        
        mtime_unit_converter = mr_model_part.GetProcessInfo()[TIME_UNIT_CONVERTER];
        mTableId = rParameters["table"].GetInt();
        
        if(mTableId != 0)
            mpTable = model_part.pGetTable(mTableId);

        KRATOS_CATCH("");
    }

    ///------------------------------------------------------------------------------------
    
    /// Destructor
    virtual ~DamTemperaturebyDeviceProcess() {}
  

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void ExecuteInitializeSolutionStep()
    {
        
        KRATOS_TRY;

        const int nelements = mr_model_part.GetMesh(mmesh_id).Elements().size();
        Variable<double> var = KratosComponents< Variable<double> >::Get(mvariable_name);
        bool IsInside = false;
        array_1d<double,3> LocalCoordinates;       
        Element::Pointer pSelectedElement;

        // Getting the values of table in case that it exist        
        if(mTableId != 0 )
        { 
            double time = mr_model_part.GetProcessInfo()[TIME];
            time = time/mtime_unit_converter;
            mvalue = mpTable->GetValue(time);
        }

        KRATOS_WATCH(mvalue)

        if (nelements != 0)
        {
            ModelPart::ElementsContainerType::iterator el_begin = mr_model_part.ElementsBegin();
            int PointsNumber = 0;

            for(int k = 0; k<nelements; k++)
            {
                ModelPart::ElementsContainerType::iterator it = el_begin + k;
                pSelectedElement = (*(it.base()));
                IsInside = pSelectedElement->GetGeometry().IsInside(mdevice_coordinates,LocalCoordinates);
                if(IsInside) break;
            }
            if(IsInside==false)
            {
                for(int k = 0; k<nelements; k++)
                {
                    ModelPart::ElementsContainerType::iterator it = el_begin + k;
                    pSelectedElement = (*(it.base()));
                    IsInside = pSelectedElement->GetGeometry().IsInside(mdevice_coordinates,LocalCoordinates,1.0e-5);
                    if(IsInside) break;
                }
            }
            if(IsInside == false)
            {
                KRATOS_ERROR << "ERROR!!, PLEASE REPEAT THE SEARCH " << std::endl;
            }

            PointsNumber = pSelectedElement->GetGeometry().PointsNumber();

            for(int j = 0; j < PointsNumber; j++)
            {
                if (pSelectedElement->GetGeometry().GetPoint(j).IsFixed(var)==false)
                {
                    pSelectedElement->GetGeometry().GetPoint(j).FastGetSolutionStepValue(var) = mvalue;
                    pSelectedElement->GetGeometry().GetPoint(j).Fix(var);
                    
                }

            }    
        }
    
        KRATOS_CATCH("");
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void ExecuteFinalizeSolutionStep()
    {
    
        KRATOS_TRY;

        const int nelements = mr_model_part.GetMesh(mmesh_id).Elements().size();
        Variable<double> var = KratosComponents< Variable<double> >::Get(mvariable_name);
        bool IsInside = false;
        array_1d<double,3> LocalCoordinates;       
        Element::Pointer pSelectedElement;

        if (nelements != 0)
        {
            ModelPart::ElementsContainerType::iterator el_begin = mr_model_part.ElementsBegin();
            int PointsNumber = 0;

            for(int k = 0; k<nelements; k++)
            {
                ModelPart::ElementsContainerType::iterator it = el_begin + k;
                pSelectedElement = (*(it.base()));
                IsInside = pSelectedElement->GetGeometry().IsInside(mdevice_coordinates,LocalCoordinates);
                if(IsInside) break;
            }
            if(IsInside==false)
            {
                for(int k = 0; k<nelements; k++)
                {
                    ModelPart::ElementsContainerType::iterator it = el_begin + k;
                    pSelectedElement = (*(it.base()));
                    IsInside = pSelectedElement->GetGeometry().IsInside(mdevice_coordinates,LocalCoordinates,1.0e-5);
                    if(IsInside) break;
                }
            }
            if(IsInside == false)
            {
                KRATOS_ERROR << "ERROR!!, PLEASE REPEAT THE SEARCH " << std::endl;
            }

            PointsNumber = pSelectedElement->GetGeometry().PointsNumber();

            for(int j = 0; j < PointsNumber; j++)
            {
                pSelectedElement->GetGeometry().GetPoint(j).Free(var);   
            }    
        }

        KRATOS_CATCH("");
    }


    /// Turn back information as a string.
    std::string Info() const
    {
        return "DamTemperaturebyDeviceProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "DamTemperaturebyDeviceProcess";
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
    bool mis_fixed;
    double mvalue;
    array_1d<double,3> mdevice_coordinates;
    double mtime_unit_converter;
    TableType::Pointer mpTable;
    int mTableId; 
    

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:

    /// Assignment operator.
    DamTemperaturebyDeviceProcess& operator=(DamTemperaturebyDeviceProcess const& rOther);

};//Class


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                DamTemperaturebyDeviceProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const DamTemperaturebyDeviceProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

} /* namespace Kratos.*/

#endif /* KRATOS_DAM_TEMPERATURE_BY_DEVICE_PROCESS defined */

