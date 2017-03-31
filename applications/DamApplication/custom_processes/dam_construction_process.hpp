//
//   Project Name:        KratosDamApplication    $
//   Last modified by:    $Author: Lorenzo Gracia $
//   Date:                $Date:       March 2017 $
//   Revision:            $Revision:          0.0 $
//

#if !defined(KRATOS_DAM_CONSTRUCTION_PROCESS )
#define  KRATOS_DAM_CONSTRUCTION_PROCESS

#include <cmath>

#include "includes/kratos_flags.h"
#include "includes/kratos_parameters.h"
#include "processes/process.h"

#include "dam_application_variables.h"

namespace Kratos
{

class DamConstructionProcess : public Process
{
    
public:

    KRATOS_CLASS_POINTER_DEFINITION(DamConstructionProcess);

    typedef ModelPart::NodesContainerType NodesArrayType;
    typedef ModelPart::ConditionsContainerType ConditionsArrayType;
    
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Constructor
    DamConstructionProcess(ModelPart& model_part,
                                Parameters rParameters
                                ) : Process(Flags()) , mr_model_part(model_part)
    {
        KRATOS_TRY
			 
        //only include validation with c++11 since raw_literals do not exist in c++03
        Parameters default_parameters( R"(
            {
                "model_part_name":"PLEASE_CHOOSE_MODEL_PART_NAME",
                "mesh_id": 0,
                "is_fixed"                                         : false,
                "Gravity_Direction"                                : "Y",
                "Reservoir_Bottom_Coordinate_in_Gravity_Direction" : 0.0,
                "Height_Dam"                                       : 0.0,
                "Number_of_phases                                  : 0.0
            }  )" );
        
        // Some values need to be mandatorily prescribed since no meaningful default value exist. For this reason try accessing to them
        // So that an error is thrown if they don't exist
        rParameters["Number_of_phases"];
        rParameters["model_part_name"];

        // Now validate agains defaults -- this also ensures no type mismatch
        rParameters.ValidateAndAssignDefaults(default_parameters);
        
        mmesh_id = rParameters["mesh_id"].GetInt();
        mis_fixed = rParameters["is_fixed"].GetBool();
        mgravity_direction = rParameters["Gravity_Direction"].GetString();
        mreference_coordinate = rParameters["Reservoir_Bottom_Coordinate_in_Gravity_Direction"].GetDouble();
        mheight = rParameters["Height_Dam"].GetDouble();
        mphases = rParameters["Number_of_phases"].GetDouble();

        mtime_unit_converter = mr_model_part.GetProcessInfo()[TIME_UNIT_CONVERTER];
  
        KRATOS_CATCH("");
    }

    ///------------------------------------------------------------------------------------
    
    /// Destructor
    virtual ~DamConstructionProcess() {}


//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void Execute()
    {
        
        KRATOS_TRY;
        
        const int nelements = mr_model_part.GetMesh(mmesh_id).Elements().size();

        if (nelements != 0)
        {
            ModelPart::ElementsContainerType::iterator el_begin = mr_model_part.ElementsBegin();
            
            #pragma omp parallel for
            for(int k = 0; k<nelements; k++)
            {
                ModelPart::ElementsContainerType::iterator it = el_begin + k;
                it->Set(ACTIVE,false);
            }

        }
        
        KRATOS_CATCH("");
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
 
 void ExecuteInitializeSolutionStep()
        {
        
        KRATOS_TRY;
        
        const int nelements = mr_model_part.GetMesh(mmesh_id).Elements().size();

        int direction;
        
        if( mgravity_direction == "X")
            direction = 1;
        else if( mgravity_direction == "Y")
            direction = 2;
        else
            direction = 3;

        double time = mr_model_part.GetProcessInfo()[TIME];
        time = time/mtime_unit_converter;

        double current_height = mreference_coordinate + (mheight/mphases)*time;

        if (nelements != 0)
        {
            ModelPart::ElementsContainerType::iterator el_begin = mr_model_part.ElementsBegin();
            
            #pragma omp parallel for
            for(int k = 0; k<nelements; k++)
            {
                ModelPart::ElementsContainerType::iterator it = el_begin + k;
                const Geometry< Node<3> >& geom = it->GetGeometry();
                const unsigned int& Dim  = geom.WorkingSpaceDimension();
                Vector central_position = geom.Center();
                central_position.resize(Dim);

                if((central_position(direction) >= mreference_coordinate) && (central_position(direction) <= current_height) )
                {
                    it->Set(ACTIVE, true);
                }
            }
        }
        
        KRATOS_CATCH("");
    }
   

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    /// Turn back information as a string.
    std::string Info() const
    {
        return "DamConstructionProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "DamConstructionProcess";
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
    std::string mgravity_direction;
    bool mis_fixed;
    double mreference_coordinate;
    double mheight;
    double mphases;
    double mtime_unit_converter;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:

    /// Assignment operator.
    DamConstructionProcess& operator=(DamConstructionProcess const& rOther);

};//Class


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  DamConstructionProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const DamConstructionProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

} /* namespace Kratos.*/

#endif /* KRATOS_DAM_CONSTRUCTION_PROCESS defined */