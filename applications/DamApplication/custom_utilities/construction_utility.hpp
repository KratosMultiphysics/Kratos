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
    
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Constructor
    ConstructionUtility(ModelPart& rMechanicalModelPart,
                        ModelPart& rThermalModelPart,
                        Parameters& rParameters
                        ) : mrMechanicalModelPart(rMechanicalModelPart) , mrThermalModelPart(rThermalModelPart)
        {
            KRATOS_TRY
    
            // Getting values 
            mMeshId = rParameters["mesh_id"].GetInt();
            mGravityDirection = rParameters["Gravity_Direction"].GetString();
            mReferenceCoordinate = rParameters["Reservoir_Bottom_Coordinate_in_Gravity_Direction"].GetDouble();
            mHeight = rParameters["Height_Dam"].GetDouble();
            mPhases = rParameters["Number_of_phases"].GetInt();
            mTimeUnitConverter = rMechanicalModelPart.GetProcessInfo()[TIME_UNIT_CONVERTER];
    
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
   
    if (nelements != 0)
    {
        ModelPart::ElementsContainerType::iterator el_begin = mrMechanicalModelPart.ElementsBegin();
        ModelPart::ElementsContainerType::iterator el_begin_thermal = mrThermalModelPart.ElementsBegin();
        
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
    
    KRATOS_CATCH("");
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void InitializeSolutionStep()
{
    KRATOS_TRY;
    
    const unsigned int nelements = mrMechanicalModelPart.GetMesh(mMeshId).Elements().size();
    int direction;
    
    if( mGravityDirection == "X")
        direction = 0;
    else if( mGravityDirection == "Y")
        direction = 1;
    else
        direction = 2;

    double time = mrMechanicalModelPart.GetProcessInfo()[TIME];
    time = time/mTimeUnitConverter;
    int int_time = time;

    KRATOS_WATCH(int_time)

    Vector vector_phase(10);
    vector_phase[0] = 1;
    vector_phase[1] = 1;
    vector_phase[2] = 2;
    vector_phase[3] = 2;
    vector_phase[4] = 3;
    vector_phase[5] = 3;
    vector_phase[6] = 4;
    vector_phase[7] = 4;
    vector_phase[8] = 5;
    vector_phase[9] = 5;
    
    

    double current_height = mReferenceCoordinate + (mHeight/mPhases)*vector_phase[int_time-1];  
    int j = 1000000;

    const unsigned int Dim = mrMechanicalModelPart.GetProcessInfo()[DOMAIN_SIZE];
    std::vector<std::size_t> ConditionNodeIds(Dim);

    if (nelements != 0)
    {
        // ELEMENTS
        ModelPart::ElementsContainerType::iterator el_begin = mrMechanicalModelPart.ElementsBegin();
        ModelPart::ElementsContainerType::iterator el_begin_thermal = mrThermalModelPart.ElementsBegin();
           
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
                            /// ES POSIBLE QUE DEBERIA HABER UN OMP_CRITICAL
                            mrThermalModelPart.CreateNewCondition("FluxCondition2D2N", j+1, ConditionNodeIds, 0);
                            j++;
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
                            /// ES POSIBLE QUE DEBERIA HABER UN OMP_CRITICAL
                            mrThermalModelPart.CreateNewCondition("FluxCondition3D3N", j+1, ConditionNodeIds, 0);
                            j++;
                        }
                    }
                }
            }

        }

    KRATOS_WATCH(mrThermalModelPart)
    }

    // Updating Heat Fluxes and Face Heat FLuxes
    //this->ActiveHeatFlux();
            
    KRATOS_CATCH("");
}

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void AfterOutputStep()
{
    const unsigned int nelements = mrThermalModelPart.GetMesh(mMeshId).Elements().size();
    const unsigned int Dim = mrMechanicalModelPart.GetProcessInfo()[DOMAIN_SIZE]; 
    
    std::vector<std::size_t> ConditionNodeIds(Dim);
    int j = 1000000;  
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
                            mrThermalModelPart.RemoveConditionFromAllLevels(j+1, 0);
                            j++;
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
                            mrThermalModelPart.RemoveConditionFromAllLevels(j+1, 0);
                            j++;
                        }
                    }
                }
            }
        }
    }
}

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void ActiveHeatFlux()
{
    KRATOS_TRY;
    
    const unsigned int nelements = mrThermalModelPart.GetMesh(mMeshId).Elements().size();    

    Vector time_vector(10);
    time_vector[0] = 0.0;
    time_vector[1] = 86400;
    time_vector[2] = 86400*2;
    time_vector[3] = 86400*3;
    time_vector[4] = 86400*4;
    time_vector[5] = 86400*5;
    time_vector[6] = 86400*6;
    time_vector[7] = 86400*7;
    time_vector[8] = 86400*8;
    time_vector[9] = 86400*9;
    
    
    double time = mrMechanicalModelPart.GetProcessInfo()[TIME];        

    int direction;
    if( mGravityDirection == "X")
        direction = 0;
    else if( mGravityDirection == "Y")
        direction = 1;
    else
        direction = 2;


    if(nelements != 0)
    {
        ModelPart::ElementsContainerType::iterator el_begin_thermal = mrThermalModelPart.ElementsBegin();               
        
        for(unsigned int j = 0; j<10; ++j)
        {
            double real_time = time - time_vector[j];
            
            // This formulation is developed using hours as temporal variable
            real_time = real_time/3600.0;

            if (real_time > 0.000001)
            {
                double current_height = mReferenceCoordinate + (mHeight/mPhases)*(j+1);
                double previous_height =mReferenceCoordinate + (mHeight/mPhases)*(j);

                // Computing the value of heat flux according the time
                //double value = mDensity*mSpecificHeat*mAlpha*mTMax*(exp(-mAlpha*time));
                double value = 2400*1*0.025*18*(exp(-0.025*real_time));
                
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

    double mH0=8.0;
    double t_sol_air = 10.0;

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
                    const double heat_flux = mH0*(t_sol_air - temp_current);
                    it_thermal->FastGetSolutionStepValue(FACE_HEAT_FLUX) = -10.0; 
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
    std::size_t mMeshId;
    std::string mGravityDirection;
    double mReferenceCoordinate;
    double mHeight;
    int mPhases;
    double mTimeUnitConverter;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

};//Class

} /* namespace Kratos.*/

#endif /* KRATOS_CONSTRUCTION_UTILITIES defined */