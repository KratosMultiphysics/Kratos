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

#if !defined(KRATOS_DAM_THERMAL_CONSTRUCTION_PROCESS )
#define  KRATOS_DAM_THERMAL_CONSTRUCTION_PROCESS

#include <cmath>

// Project includes
#include "includes/kratos_flags.h"
#include "includes/kratos_parameters.h"
#include "processes/process.h"

// Application include
#include "dam_application_variables.h"

namespace Kratos
{

class DamThermalConstructionProcess : public Process
{
    
public:

    KRATOS_CLASS_POINTER_DEFINITION(DamThermalConstructionProcess);
    
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Constructor
    DamThermalConstructionProcess(ModelPart& rModelPart,
                                Parameters rParameters
                                ) : Process(Flags()) , mrModelPart(rModelPart)
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
                "Number_of_phases"                                 : 0.0
            }  )" );
        
        // Some values need to be mandatorily prescribed since no meaningful default value exist. For this reason try accessing to them
        // So that an error is thrown if they don't exist
        rParameters["Number_of_phases"];
        rParameters["model_part_name"];

        // Now validate agains defaults -- this also ensures no type mismatch
        rParameters.ValidateAndAssignDefaults(default_parameters);
        
        mMeshId = rParameters["mesh_id"].GetInt();
        mIsFixed = rParameters["is_fixed"].GetBool();
        mGravityDirection = rParameters["Gravity_Direction"].GetString();
        mReferenceCoordinate = rParameters["Reservoir_Bottom_Coordinate_in_Gravity_Direction"].GetDouble();
        mHeight = rParameters["Height_Dam"].GetDouble();
        mPhases = rParameters["Number_of_phases"].GetDouble();

        mTimeUnitConverter = mrModelPart.GetProcessInfo()[TIME_UNIT_CONVERTER];
  

        KRATOS_CATCH("");
    }

    ///------------------------------------------------------------------------------------
    
    /// Destructor
    virtual ~DamThermalConstructionProcess() {}


//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void ExecuteInitialize()
    {
        KRATOS_TRY;
        
        const int nelements = mrModelPart.GetMesh(mMeshId).Elements().size();
        const int nnodes = mrModelPart.GetMesh(mMeshId).Nodes().size();
        
        if (nelements != 0)
        {
            ModelPart::ElementsContainerType::iterator el_begin = mrModelPart.ElementsBegin();
            
            #pragma omp parallel for
            for(int k = 0; k<nelements; k++)
            {
                ModelPart::ElementsContainerType::iterator it = el_begin + k;
                it->Set(ACTIVE,false);
            }

            ModelPart::NodesContainerType::iterator it_begin = mrModelPart.NodesBegin();
            
            #pragma omp parallel for
            for(int i = 0; i<nnodes; i++)
            {
                ModelPart::NodesContainerType::iterator it = it_begin + i;
                it->Set(ACTIVE,false);
            }

        }
        
        KRATOS_CATCH("");
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
 
    void ExecuteInitializeSolutionStep()
    {
        KRATOS_TRY;
        
        const int nelements = mrModelPart.GetMesh(mMeshId).Elements().size();
        const int nnodes = mrModelPart.GetMesh(mMeshId).Nodes().size();            
        int direction;
        
        if( mGravityDirection == "X")
            direction = 0;
        else if( mGravityDirection == "Y")
            direction = 1;
        else
            direction = 2;

        double time = mrModelPart.GetProcessInfo()[TIME];
        time = time/mTimeUnitConverter;
        double current_height = mReferenceCoordinate + (mHeight/mPhases)*time;  
        double value;
        int j = 1000;
        std::vector<std::size_t> ConditionNodeIds(3);

        if (nelements != 0)
        {
            //NODES
            ModelPart::NodesContainerType::iterator it_begin = mrModelPart.NodesBegin();
        
            #pragma omp parallel for
            for(int i = 0; i<nnodes; i++)
            {
                ModelPart::NodesContainerType::iterator it = it_begin + i;                    
                if( mGravityDirection == "X")
                    value = it->X();
                else if( mGravityDirection == "Y")
                    value = it->Y();
                else
                    value = it->Z();
                
                if((value >= mReferenceCoordinate) && (value <= (current_height+0.00001) ))
                {
                    it->Set(ACTIVE,true);
                }
            }
            // ELEMENTS
            ModelPart::ElementsContainerType::iterator el_begin = mrModelPart.ElementsBegin();               
            
            #pragma omp parallel for
            for(int k = 0; k<nelements; k++)
            {
                ModelPart::ElementsContainerType::iterator it = el_begin + k;
                const Geometry< Node<3> >& geom = it->GetGeometry();
                array_1d<double,3> central_position = geom.Center();
                if((central_position(direction) >= mReferenceCoordinate) && (central_position(direction) <= current_height) )
                {
                    it->Set(ACTIVE,true);
                }
            }
            #pragma omp parallel for
            for(int k = 0; k<nelements; k++)
            {
                ModelPart::ElementsContainerType::iterator it = el_begin + k;
                // Elements
                if((it)->Is(ACTIVE) == false)
                {
                    for (unsigned int i_face = 0; i_face < (*it).GetGeometry().FacesNumber(); i_face++)
                    {
                        const unsigned int number_of_points = (*it).GetGeometry().Faces()[i_face].PointsNumber();
                        unsigned int count = 0;
                        for (unsigned int i_node = 0; i_node < number_of_points; i_node++)
                        {
                            if ((*it).GetGeometry().Faces()[i_face][i_node].Is(ACTIVE) == true)
                            {
                                count++;
                            }
                        }

                        if (count == number_of_points)
                        {
                            for  (unsigned int m = 0; m < number_of_points; m++)
                            {
                                ConditionNodeIds[m] = (*it).GetGeometry().Faces()[i_face][m].Id();
                                KRATOS_WATCH(ConditionNodeIds[m])
                            }
                            KRATOS_WATCH("OTROOOOOOOOOOOOOOOo")                            
                            //mrModelPart.GetSubModelPart("Thermal_Part_Auto_1").CreateNewCondition("FluxCondition3D3N", j+1, ConditionNodeIds, 0);
                            j++;
                        }
                    }
                }
            }

            KRATOS_WATCH(mrModelPart)
        }
        
        KRATOS_CATCH("");
    }
   

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    /// Turn back information as a string.
    std::string Info() const
    {
        return "DamThermalConstructionProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "DamThermalConstructionProcess";
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
    std::string mGravityDirection;
    bool mIsFixed;
    double mReferenceCoordinate;
    double mHeight;
    double mPhases;
    double mTimeUnitConverter;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:

    /// Assignment operator.
    DamThermalConstructionProcess& operator=(DamThermalConstructionProcess const& rOther);

};//Class


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
    DamThermalConstructionProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const DamThermalConstructionProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

} /* namespace Kratos.*/

#endif /* KRATOS_DAM_THERMAL_CONSTRUCTION_PROCESS defined */