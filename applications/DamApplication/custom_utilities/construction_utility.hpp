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

#if !defined(KRATOS_CONSTRUCTION_UTILITIES )
#define  KRATOS_CONSTRUCTION_UTILITIES

#include <cmath>

// Project includes
#include "includes/kratos_flags.h"
#include "includes/kratos_parameters.h"
#include "processes/process.h"
#include "utilities/openmp_utils.h"
#include "utilities/math_utils.h"
#include "includes/element.h"

// Application includes
#include "dam_application_variables.h"

namespace Kratos
{
    
class ConstructionUtility
{

public:
    
    typedef std::size_t IndexType;
    typedef Table<double,double> TableType;  
    
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Constructor
    ConstructionUtility(ModelPart& rMechanicalModelPart, ModelPart& rThermalModelPart,
                        TableType& rTablePhases, TableType& rTableAmbientTemp, Parameters& rParameters
                        ) : mrMechanicalModelPart(rMechanicalModelPart) , mrThermalModelPart(rThermalModelPart) , mrTablePhases(rTablePhases) , mrTableAmbientTemp(rTableAmbientTemp)
        {
            KRATOS_TRY
    
            // Getting values 
            mMeshId = rParameters["mesh_id"].GetInt();
            mGravityDirection = rParameters["gravity_direction"].GetString();
            mReferenceCoordinate = rParameters["reservoir_bottom_coordinate_in_gravity_direction"].GetDouble();
            mHeight = rParameters["height_dam"].GetDouble();
            mPhases = rParameters["number_of_phases"].GetInt();
            mDensity = rParameters["density"].GetDouble();
            mSpecificHeat = rParameters["specific_heat"].GetDouble();
            mAlpha = rParameters["alpha"].GetDouble();
            mTMax = rParameters["tmax"].GetDouble();
            mH0 = rParameters["h_0"].GetDouble();
            mTimeUnitConverter = mrMechanicalModelPart.GetProcessInfo()[TIME_UNIT_CONVERTER];    
            
            KRATOS_CATCH("");
        }

    ///------------------------------------------------------------------------------------
    
    /// Destructor
    ~ConstructionUtility() {}


//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void Initialize()
{
    KRATOS_TRY;
    
    const unsigned int nelements = mrMechanicalModelPart.GetMesh(mMeshId).Elements().size();   
    const unsigned int nnodes = mrMechanicalModelPart.GetMesh(mMeshId).Nodes().size();

    mMechanicalLastCondition = mrMechanicalModelPart.GetMesh(mMeshId).Conditions().size();
    mThermalLastCondition = mrThermalModelPart.GetMesh(mMeshId).Conditions().size();
   
    if (nelements != 0)
    {
        ModelPart::ElementsContainerType::iterator el_begin = mrMechanicalModelPart.ElementsBegin();
        ModelPart::ElementsContainerType::iterator el_begin_thermal = mrThermalModelPart.ElementsBegin();
        mNumNode = el_begin->GetGeometry().PointsNumber();

        #pragma omp parallel for
        for(unsigned int k = 0; k<nelements; ++k)
        {
            ModelPart::ElementsContainerType::iterator it = el_begin + k;
            ModelPart::ElementsContainerType::iterator it_thermal = el_begin_thermal + k;
            it->Set(ACTIVE,false);
            it_thermal->Set(ACTIVE,false);
        }

        // Same nodes for both computing model part
        ModelPart::NodesContainerType::iterator it_begin = mrMechanicalModelPart.NodesBegin();        
        #pragma omp parallel for
        for(unsigned int i = 0; i<nnodes; ++i)
        {
            ModelPart::NodesContainerType::iterator it = it_begin + i;
            it->Set(ACTIVE,false);
        }

    }

    const unsigned int soil_nelements = mrMechanicalModelPart.GetSubModelPart("Parts_Parts_Auto9").Elements().size();   
    const unsigned int soil_nnodes = mrMechanicalModelPart.GetSubModelPart("Parts_Parts_Auto9").Nodes().size();

    if (soil_nelements != 0)
    {
        ModelPart::ElementsContainerType::iterator el_begin = mrMechanicalModelPart.GetSubModelPart("Parts_Parts_Auto9").ElementsBegin();
        ModelPart::ElementsContainerType::iterator el_begin_thermal = mrThermalModelPart.GetSubModelPart("Thermal_Part_Auto_9").ElementsBegin();
        mNumNode = el_begin->GetGeometry().PointsNumber();

        #pragma omp parallel for
        for(unsigned int k = 0; k<soil_nelements; ++k)
        {
            ModelPart::ElementsContainerType::iterator it = el_begin + k;
            ModelPart::ElementsContainerType::iterator it_thermal = el_begin_thermal + k;
            it->Set(ACTIVE,true);
            it_thermal->Set(ACTIVE,true);
        }

        // Same nodes for both computing model part
        ModelPart::NodesContainerType::iterator it_begin = mrMechanicalModelPart.GetSubModelPart("Parts_Parts_Auto9").NodesBegin();        
        #pragma omp parallel for
        for(unsigned int i = 0; i<soil_nnodes; ++i)
        {
            ModelPart::NodesContainerType::iterator it = it_begin + i;
            it->Set(ACTIVE,true);
        }

    }
    
    KRATOS_CATCH("");
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void InitializeSolutionStep(std::string ThermalSubModelPartName, std::string MechanicalSubModelPartName)
{
    KRATOS_TRY;
    
    const unsigned int nelements = mrThermalModelPart.GetSubModelPart(ThermalSubModelPartName).Elements().size();
    int direction;
    
    if( mGravityDirection == "X")
        direction = 0;
    else if( mGravityDirection == "Y")
        direction = 1;
    else
        direction = 2;

    double time = mrMechanicalModelPart.GetProcessInfo()[TIME];
    int int_time = time/mTimeUnitConverter;
    
    // Getting the value of the table and computing the current height
    unsigned int current_number_of_phases = mrTablePhases.GetValue(int_time-1);
    double current_height = mReferenceCoordinate + (mHeight/mPhases)*current_number_of_phases;

    if (nelements != 0)
    {
        // ELEMENTS
        ModelPart::ElementsContainerType::iterator el_begin = mrMechanicalModelPart.GetSubModelPart(MechanicalSubModelPartName).ElementsBegin();
        ModelPart::ElementsContainerType::iterator el_begin_thermal = mrThermalModelPart.GetSubModelPart(ThermalSubModelPartName).ElementsBegin();
           
        #pragma omp parallel for
        for(unsigned int k = 0; k<nelements; ++k)
        {
            ModelPart::ElementsContainerType::iterator it = el_begin + k;
            ModelPart::ElementsContainerType::iterator it_thermal = el_begin_thermal + k;
            array_1d<double,3> central_position = it->GetGeometry().Center();

            if((central_position(direction) >= mReferenceCoordinate) && (central_position(direction) <= current_height) )
            {
                it->Set(ACTIVE,true);
                it_thermal->Set(ACTIVE,true);

                const unsigned int number_of_points = it_thermal->GetGeometry().PointsNumber();
                for(unsigned int i = 0; i<number_of_points; i++)
                {
                    it->GetGeometry()[i].Set(ACTIVE,true);                    
                }
            }
        }
    }

    // Detecting and creating face heat flux
    this->SearchingFluxes();
           
    KRATOS_CATCH("");
}

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void AfterOutputStep()
{
    KRATOS_TRY;
    
    const unsigned int nelements = mrThermalModelPart.GetMesh(mMeshId).Elements().size();
    const unsigned int Dim = mrMechanicalModelPart.GetProcessInfo()[DOMAIN_SIZE]; 
    std::vector<std::size_t> ConditionNodeIds(Dim);
    if (mNumNode == 8)
        ConditionNodeIds.resize(Dim+1);

    int last_condition_id = mMechanicalLastCondition + mThermalLastCondition;
    
    if (nelements != 0)
    {
        ModelPart::ElementsContainerType::iterator el_begin_thermal = mrThermalModelPart.ElementsBegin();               
        
        if (Dim == 2)
        {
            #pragma omp parallel for
            for(unsigned int k = 0; k<nelements; ++k)
            {
                ModelPart::ElementsContainerType::iterator it_thermal = el_begin_thermal + k;
                // Elements
                if((it_thermal)->Is(ACTIVE) == false)
                {
                    for (unsigned int i_edge = 0; i_edge < (*it_thermal).GetGeometry().EdgesNumber(); ++i_edge)
                    {
                        const unsigned int number_of_points = (*it_thermal).GetGeometry().Edges()[i_edge].PointsNumber();
                        unsigned int count = 0;
                        for (unsigned int i_node = 0; i_node < number_of_points; ++i_node)
                        {
                            if ((*it_thermal).GetGeometry().Edges()[i_edge][i_node].Is(ACTIVE) == true)
                            {
                                count++;
                            }
                        }
                        if (count == number_of_points)
                        {
                            for  (unsigned int m = 0; m < number_of_points; ++m)
                            {
                                ConditionNodeIds[m] = (*it_thermal).GetGeometry().Edges()[i_edge][m].Id();
                            }
                            this->DeactiveFaceHeatFluxStep(ConditionNodeIds);
                            #pragma omp critical
                            {
                                mrThermalModelPart.RemoveConditionFromAllLevels(last_condition_id+1, 0);
                                last_condition_id++;
                            }
                        }
                    }
                }
            }
        } else
        {
            #pragma omp parallel for
            for(unsigned int k = 0; k<nelements; ++k)
            {
                ModelPart::ElementsContainerType::iterator it_thermal = el_begin_thermal + k;
                // Elements
                if((it_thermal)->Is(ACTIVE) == false)
                {
                    for (unsigned int i_face = 0; i_face < (*it_thermal).GetGeometry().FacesNumber(); ++i_face)
                    {
                        const unsigned int number_of_points = (*it_thermal).GetGeometry().Faces()[i_face].PointsNumber();
                        unsigned int count = 0;
                        for (unsigned int i_node = 0; i_node < number_of_points; ++i_node)
                        {
                            if ((*it_thermal).GetGeometry().Faces()[i_face][i_node].Is(ACTIVE) == true)
                            {
                                count++;
                            }
                        }
                        if (count == number_of_points)
                        {
                            for  (unsigned int m = 0; m < number_of_points; ++m)
                            {
                                ConditionNodeIds[m] = (*it_thermal).GetGeometry().Faces()[i_face][m].Id();
                            }
                            this->DeactiveFaceHeatFluxStep(ConditionNodeIds);
                            #pragma omp critical
                            {
                                mrThermalModelPart.RemoveConditionFromAllLevels(last_condition_id+1, 0);
                                last_condition_id++;   
                            }
                        }
                    }
                }
            }
        }
    }

    KRATOS_CATCH("");    
}


///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void SearchingFluxes()
{
    KRATOS_TRY;

    const unsigned int nelements = mrMechanicalModelPart.GetMesh(mMeshId).Elements().size();
    const unsigned int Dim = mrMechanicalModelPart.GetProcessInfo()[DOMAIN_SIZE];
    std::vector<std::size_t> ConditionNodeIds(Dim);
    if (mNumNode == 8)
        ConditionNodeIds.resize(Dim+1);
    
    int last_condition_id = mMechanicalLastCondition + mThermalLastCondition;
        
    if (nelements != 0)
    {
        ModelPart::ElementsContainerType::iterator el_begin_thermal = mrThermalModelPart.ElementsBegin();
        
        if (Dim == 2)
        {
            // Searching for thermal boundary conditions Edges
            #pragma omp parallel for
            for(unsigned int k = 0; k<nelements; ++k)
            {
                ModelPart::ElementsContainerType::iterator it_thermal = el_begin_thermal + k;
                // Elements
                if((it_thermal)->Is(ACTIVE) == false)
                {
                    for (unsigned int i_edge = 0; i_edge < (*it_thermal).GetGeometry().EdgesNumber(); ++i_edge)
                    {
                        const unsigned int number_of_points = (*it_thermal).GetGeometry().Edges()[i_edge].PointsNumber();
                        unsigned int count = 0;
                    
                        for (unsigned int i_node = 0; i_node < number_of_points; ++i_node)
                        {
                            if ((*it_thermal).GetGeometry().Edges()[i_edge][i_node].Is(ACTIVE) == true)
                            {
                                count++;
                            }
                        }
                        if (count == number_of_points)
                        {
                            for  (unsigned int m = 0; m < number_of_points; ++m)
                            {
                                ConditionNodeIds[m] = (*it_thermal).GetGeometry().Edges()[i_edge][m].Id();
                            }                       
                            this->ActiveFaceHeatFluxStep(ConditionNodeIds);
                            #pragma omp critical
                            {
                                mrThermalModelPart.CreateNewCondition("FluxCondition2D2N", last_condition_id+1, ConditionNodeIds, 0);
                                last_condition_id++;
                            }
                        }
                    }
                }
            }

        }
        else
        {
            // Searching for thermal boundary conditions
            #pragma omp parallel for
            for(unsigned int k = 0; k<nelements; ++k)
            {
                ModelPart::ElementsContainerType::iterator it_thermal = el_begin_thermal + k;
                // Elements
                if((it_thermal)->Is(ACTIVE) == false)
                {
                    for (unsigned int i_face = 0; i_face < (*it_thermal).GetGeometry().FacesNumber(); ++i_face)
                    {
                        const unsigned int number_of_points = (*it_thermal).GetGeometry().Faces()[i_face].PointsNumber();
                        unsigned int count = 0;
                    
                        for (unsigned int i_node = 0; i_node < number_of_points; ++i_node)
                        {
                            if ((*it_thermal).GetGeometry().Faces()[i_face][i_node].Is(ACTIVE) == true)
                            {
                                count++;
                            }
                        }
                        if (count == number_of_points)
                        {
                            for  (unsigned int m = 0; m < number_of_points; ++m)
                            {
                                ConditionNodeIds[m] = (*it_thermal).GetGeometry().Faces()[i_face][m].Id();
                            }                       
                            this->ActiveFaceHeatFluxStep(ConditionNodeIds);
                            #pragma omp critical
                            {
                                if (number_of_points==3)
                                {
                                    mrThermalModelPart.CreateNewCondition("FluxCondition3D3N", last_condition_id+1, ConditionNodeIds, 0);
                                    last_condition_id++;
                                }
                                else
                                {
                                    mrThermalModelPart.CreateNewCondition("FluxCondition3D4N", last_condition_id+1, ConditionNodeIds, 0);
                                    last_condition_id++;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    KRATOS_CATCH("");    
}

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void ActiveHeatFlux(std::string ThermalSubModelPartName, int phase, double phase_time)
{
    KRATOS_TRY;
    
    const unsigned int nelements = mrThermalModelPart.GetSubModelPart(ThermalSubModelPartName).Elements().size();

    KRATOS_WATCH(nelements)

    int direction;
    if( mGravityDirection == "X")
        direction = 0;
    else if( mGravityDirection == "Y")
        direction = 1;
    else
        direction = 2;

    if(nelements != 0)
    {
        ModelPart::ElementsContainerType::iterator el_begin_thermal = mrThermalModelPart.GetSubModelPart(ThermalSubModelPartName).ElementsBegin();               
        
        double current_height = mReferenceCoordinate + (mHeight/mPhases)*(phase);
        double previous_height = mReferenceCoordinate + (mHeight/mPhases)*(phase-1);         

        // This formulation is developed using hours as temporal variable
        double real_time = phase_time/3600.0;
        
        // Computing the value of heat flux according the time
        double value = mDensity*mSpecificHeat*mAlpha*mTMax*(exp(-mAlpha*real_time));

        KRATOS_WATCH(value)
        
        #pragma omp parallel for
        for(unsigned int k = 0; k<nelements; ++k)
        {
            ModelPart::ElementsContainerType::iterator it_thermal = el_begin_thermal + k;
            array_1d<double,3> central_position = it_thermal->GetGeometry().Center();
            if((central_position(direction) >= previous_height) && (central_position(direction) <= current_height) )
            {
                const unsigned int number_of_points = it_thermal->GetGeometry().PointsNumber();
                for(unsigned int i = 0; i<number_of_points; ++i)
                {
                    it_thermal->GetGeometry()[i].FastGetSolutionStepValue(HEAT_FLUX) = value;                    
                }
            }
        }             
    }

    KRATOS_CATCH("");
}

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void ActiveFaceHeatFluxStep(std::vector<IndexType> ConditionNodeIds)
{
    KRATOS_TRY;

    const unsigned int size = ConditionNodeIds.size();
    const unsigned int nnodes = mrThermalModelPart.GetMesh(mMeshId).Nodes().size();
    ModelPart::NodesContainerType::iterator it_begin_thermal = mrThermalModelPart.GetMesh(mMeshId).NodesBegin();

    double time = mrThermalModelPart.GetProcessInfo()[TIME];
    time = time/mTimeUnitConverter;
    double ambient_temp = mrTableAmbientTemp.GetValue(time-1);

    if(size != 0)
    {
        for(unsigned int j = 0; j<size; ++j)
        {
            for(unsigned int i = 0; i<nnodes; ++i) 
            {
                ModelPart::NodesContainerType::iterator it_thermal = it_begin_thermal + i;

                if (it_thermal->Id() == ConditionNodeIds[j]) 
                {
                    const double temp_current = it_thermal->FastGetSolutionStepValue(TEMPERATURE);
                    const double heat_flux = mH0*(ambient_temp - temp_current);
                    it_thermal->FastGetSolutionStepValue(FACE_HEAT_FLUX) = heat_flux; 
                }
            }
        }    
    }

    KRATOS_CATCH("");
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void DeactiveFaceHeatFluxStep(std::vector<IndexType> ConditionNodeIds)
{
    KRATOS_TRY;

    const unsigned int size = ConditionNodeIds.size();
    const unsigned int nnodes = mrThermalModelPart.GetMesh(mMeshId).Nodes().size();
    ModelPart::NodesContainerType::iterator it_begin_thermal = mrThermalModelPart.GetMesh(mMeshId).NodesBegin();

    if(size != 0)
    {
        for(unsigned int j = 0; j<size; ++j)
        {
            for(unsigned int i = 0; i<nnodes; ++i) 
            {
                ModelPart::NodesContainerType::iterator it_thermal = it_begin_thermal + i;

                if (it_thermal->Id() == ConditionNodeIds[j]) 
                {   
                    //Setting to 0 the fluxes since in the next step 
                    it_thermal->FastGetSolutionStepValue(FACE_HEAT_FLUX) = 0.0; 
                }
            }
        }    
    }

    KRATOS_CATCH("");
}

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    /// Member Variables

    ModelPart& mrMechanicalModelPart;
    ModelPart& mrThermalModelPart;
    int mNumNode;
    std::size_t mMeshId;
    std::string mGravityDirection;
    double mReferenceCoordinate;
    double mHeight;
    int mPhases;
    double mTimeUnitConverter;    
    double mDensity;
    double mSpecificHeat;
    double mAlpha;
    double mTMax;
    double mH0;
    unsigned int mMechanicalLastCondition;
    unsigned int mThermalLastCondition;
    TableType& mrTablePhases;
    TableType& mrTableAmbientTemp;
    


//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

};//Class

} /* namespace Kratos.*/

#endif /* KRATOS_CONSTRUCTION_UTILITIES defined */