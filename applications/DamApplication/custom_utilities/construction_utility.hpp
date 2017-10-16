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

    Vector vector_phase(150);
    vector_phase[0] = 1;
    vector_phase[1] = 1;
    vector_phase[2] = 1;
    vector_phase[3] = 1;
    vector_phase[4] = 1;
    vector_phase[5] = 2;
    vector_phase[6] = 2;
    vector_phase[7] = 2;
    vector_phase[8] = 2;
    vector_phase[9] = 2;
    vector_phase[10] = 3;
    vector_phase[11] = 3;
    vector_phase[12] = 3;
    vector_phase[13] = 3;
    vector_phase[14] = 3;
    vector_phase[15] = 4;
    vector_phase[16] = 4;
    vector_phase[17] = 4;
    vector_phase[18] = 4;
    vector_phase[19] = 4;
    vector_phase[20] = 5;
    vector_phase[21] = 5;
    vector_phase[22] = 5;
    vector_phase[23] = 5;
    vector_phase[24] = 5;
    vector_phase[25] = 6;
    vector_phase[26] = 6;
    vector_phase[27] = 6;
    vector_phase[28] = 6;
    vector_phase[29] = 6;
    vector_phase[30] = 7;
    vector_phase[31] = 7;
    vector_phase[32] = 7;
    vector_phase[33] = 7;
    vector_phase[34] = 7;
    vector_phase[35] = 8;
    vector_phase[36] = 8;
    vector_phase[37] = 8;
    vector_phase[38] = 8;
    vector_phase[39] = 8;
    vector_phase[40] = 9;
    vector_phase[41] = 9;
    vector_phase[42] = 9;
    vector_phase[43] = 9;
    vector_phase[44] = 9;
    vector_phase[45] = 10;
    vector_phase[46] = 10;
    vector_phase[47] = 10;
    vector_phase[48] = 10;
    vector_phase[49] = 10;
    vector_phase[50] = 10;
    vector_phase[51] = 10;
    vector_phase[52] = 10;
    vector_phase[53] = 10;
    vector_phase[54] = 10;
    vector_phase[55] = 10;
    vector_phase[56] = 10;
    vector_phase[57] = 10;
    vector_phase[58] = 10;
    vector_phase[59] = 10;
    vector_phase[60] = 10;
    vector_phase[61] = 10;
    vector_phase[62] = 10;
    vector_phase[63] = 10;
    vector_phase[64] = 10;
    vector_phase[65] = 10;
    vector_phase[66] = 10;
    vector_phase[67] = 10;
    vector_phase[68] = 10;
    vector_phase[69] = 10;
    vector_phase[70] = 10;
    vector_phase[71] = 10;
    vector_phase[72] = 10;
    vector_phase[73] = 10;
    vector_phase[74] = 10;
    vector_phase[75] = 10;
    vector_phase[76] = 10;
    vector_phase[77] = 10;
    vector_phase[78] = 10;
    vector_phase[79] = 10;
    vector_phase[80] = 10;
    vector_phase[81] = 10;
    vector_phase[82] = 10;
    vector_phase[83] = 10;
    vector_phase[84] = 10;
    vector_phase[85] = 10;
    vector_phase[86] = 10;
    vector_phase[87] = 10;
    vector_phase[88] = 10;
    vector_phase[89] = 10;
    vector_phase[90] = 10;
    vector_phase[91] = 10;
    vector_phase[92] = 10;
    vector_phase[93] = 10;
    vector_phase[94] = 10;
    vector_phase[95] = 10;
    vector_phase[96] = 10;
    vector_phase[97] = 10;
    vector_phase[98] = 10;
    vector_phase[99] = 10;
    vector_phase[100] = 10;
    vector_phase[101] = 10;
    vector_phase[102] = 10;
    vector_phase[103] = 10;
    vector_phase[104] = 10;
    vector_phase[105] = 10;
    vector_phase[106] = 10;
    vector_phase[107] = 10;
    vector_phase[108] = 10;
    vector_phase[109] = 10;
    vector_phase[110] = 10;
    vector_phase[111] = 10;
    vector_phase[112] = 10;
    vector_phase[113] = 10;
    vector_phase[114] = 10;
    vector_phase[115] = 10;
    vector_phase[116] = 10;
    vector_phase[117] = 10;
    vector_phase[118] = 10;
    vector_phase[119] = 10;
    vector_phase[120] = 10;    
    vector_phase[121] = 10;
    vector_phase[122] = 10;
    vector_phase[123] = 10;
    vector_phase[124] = 10;
    vector_phase[125] = 10;
    vector_phase[126] = 10;
    vector_phase[127] = 10;
    vector_phase[128] = 10;
    vector_phase[129] = 10;
    vector_phase[130] = 10;    
    vector_phase[131] = 10;
    vector_phase[132] = 10;
    vector_phase[133] = 10;
    vector_phase[134] = 10;
    vector_phase[135] = 10;
    vector_phase[136] = 10;
    vector_phase[137] = 10;
    vector_phase[138] = 10;
    vector_phase[139] = 10;
    vector_phase[140] = 10;    
    vector_phase[141] = 10;
    vector_phase[142] = 10;
    vector_phase[143] = 10;
    vector_phase[144] = 10;
    vector_phase[145] = 10;
    vector_phase[146] = 10;
    vector_phase[147] = 10;
    vector_phase[148] = 10;
    vector_phase[149] = 10;

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
    this->ActiveHeatFlux();
            
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
    double time = mrMechanicalModelPart.GetProcessInfo()[TIME];    
    const double delta_time = mrMechanicalModelPart.GetProcessInfo()[DELTA_TIME];        
    int int_time = time/mTimeUnitConverter;

    KRATOS_WATCH(int_time)

    Vector vector_phase(150);
    vector_phase[0] = 1;
    vector_phase[1] = 1;
    vector_phase[2] = 1;
    vector_phase[3] = 1;
    vector_phase[4] = 1;
    vector_phase[5] = 2;
    vector_phase[6] = 2;
    vector_phase[7] = 2;
    vector_phase[8] = 2;
    vector_phase[9] = 2;
    vector_phase[10] = 3;
    vector_phase[11] = 3;
    vector_phase[12] = 3;
    vector_phase[13] = 3;
    vector_phase[14] = 3;
    vector_phase[15] = 4;
    vector_phase[16] = 4;
    vector_phase[17] = 4;
    vector_phase[18] = 4;
    vector_phase[19] = 4;
    vector_phase[20] = 5;
    vector_phase[21] = 5;
    vector_phase[22] = 5;
    vector_phase[23] = 5;
    vector_phase[24] = 5;
    vector_phase[25] = 6;
    vector_phase[26] = 6;
    vector_phase[27] = 6;
    vector_phase[28] = 6;
    vector_phase[29] = 6;
    vector_phase[30] = 7;
    vector_phase[31] = 7;
    vector_phase[32] = 7;
    vector_phase[33] = 7;
    vector_phase[34] = 7;
    vector_phase[35] = 8;
    vector_phase[36] = 8;
    vector_phase[37] = 8;
    vector_phase[38] = 8;
    vector_phase[39] = 8;
    vector_phase[40] = 9;
    vector_phase[41] = 9;
    vector_phase[42] = 9;
    vector_phase[43] = 9;
    vector_phase[44] = 9;
    vector_phase[45] = 10;
    vector_phase[46] = 10;
    vector_phase[47] = 10;
    vector_phase[48] = 10;
    vector_phase[49] = 10;
    vector_phase[50] = 10;
    vector_phase[51] = 10;
    vector_phase[52] = 10;
    vector_phase[53] = 10;
    vector_phase[54] = 10;
    vector_phase[55] = 10;
    vector_phase[56] = 10;
    vector_phase[57] = 10;
    vector_phase[58] = 10;
    vector_phase[59] = 10;
    vector_phase[60] = 10;
    vector_phase[61] = 10;
    vector_phase[62] = 10;
    vector_phase[63] = 10;
    vector_phase[64] = 10;
    vector_phase[65] = 10;
    vector_phase[66] = 10;
    vector_phase[67] = 10;
    vector_phase[68] = 10;
    vector_phase[69] = 10;
    vector_phase[70] = 10;
    vector_phase[71] = 10;
    vector_phase[72] = 10;
    vector_phase[73] = 10;
    vector_phase[74] = 10;
    vector_phase[75] = 10;
    vector_phase[76] = 10;
    vector_phase[77] = 10;
    vector_phase[78] = 10;
    vector_phase[79] = 10;
    vector_phase[80] = 10;
    vector_phase[81] = 10;
    vector_phase[82] = 10;
    vector_phase[83] = 10;
    vector_phase[84] = 10;
    vector_phase[85] = 10;
    vector_phase[86] = 10;
    vector_phase[87] = 10;
    vector_phase[88] = 10;
    vector_phase[89] = 10;
    vector_phase[90] = 10;
    vector_phase[91] = 10;
    vector_phase[92] = 10;
    vector_phase[93] = 10;
    vector_phase[94] = 10;
    vector_phase[95] = 10;
    vector_phase[96] = 10;
    vector_phase[97] = 10;
    vector_phase[98] = 10;
    vector_phase[99] = 10;
    vector_phase[100] = 10;
    vector_phase[101] = 10;
    vector_phase[102] = 10;
    vector_phase[103] = 10;
    vector_phase[104] = 10;
    vector_phase[105] = 10;
    vector_phase[106] = 10;
    vector_phase[107] = 10;
    vector_phase[108] = 10;
    vector_phase[109] = 10;
    vector_phase[110] = 10;
    vector_phase[111] = 10;
    vector_phase[112] = 10;
    vector_phase[113] = 10;
    vector_phase[114] = 10;
    vector_phase[115] = 10;
    vector_phase[116] = 10;
    vector_phase[117] = 10;
    vector_phase[118] = 10;
    vector_phase[119] = 10;
    vector_phase[120] = 10;    
    vector_phase[121] = 10;
    vector_phase[122] = 10;
    vector_phase[123] = 10;
    vector_phase[124] = 10;
    vector_phase[125] = 10;
    vector_phase[126] = 10;
    vector_phase[127] = 10;
    vector_phase[128] = 10;
    vector_phase[129] = 10;
    vector_phase[130] = 10;    
    vector_phase[131] = 10;
    vector_phase[132] = 10;
    vector_phase[133] = 10;
    vector_phase[134] = 10;
    vector_phase[135] = 10;
    vector_phase[136] = 10;
    vector_phase[137] = 10;
    vector_phase[138] = 10;
    vector_phase[139] = 10;
    vector_phase[140] = 10;    
    vector_phase[141] = 10;
    vector_phase[142] = 10;
    vector_phase[143] = 10;
    vector_phase[144] = 10;
    vector_phase[145] = 10;
    vector_phase[146] = 10;
    vector_phase[147] = 10;
    vector_phase[148] = 10;
    vector_phase[149] = 10;

    Vector time_vector(10);
    time_vector[0] = 0.0;
    time_vector[1] = 4*delta_time;
    time_vector[2] = 9*delta_time;
    time_vector[3] = 14*delta_time;
    time_vector[4] = 19*delta_time;
    time_vector[5] = 24*delta_time;
    time_vector[6] = 29*delta_time;
    time_vector[7] = 34*delta_time;
    time_vector[8] = 39*delta_time;
    time_vector[9] = 44*delta_time;
    
    int direction;
    if( mGravityDirection == "X")
        direction = 0;
    else if( mGravityDirection == "Y")
        direction = 1;
    else
        direction = 2;

    unsigned int current_number_of_phases = vector_phase[int_time];

    KRATOS_WATCH(current_number_of_phases)

    if(nelements != 0)
    {
        ModelPart::ElementsContainerType::iterator el_begin_thermal = mrThermalModelPart.ElementsBegin();               
        
        for(unsigned int j = 0; j<current_number_of_phases; ++j)
        {
            double current_height = mReferenceCoordinate + (mHeight/mPhases)*(j+1);
            double previous_height =mReferenceCoordinate + (mHeight/mPhases)*(j);
            double real_time = time - time_vector[j];
            KRATOS_WATCH("SALTOOOOOOOOOOOOOOOOOOOOOOOO")
            KRATOS_WATCH(j)
            KRATOS_WATCH(current_height)
            KRATOS_WATCH(previous_height)
            KRATOS_WATCH(real_time)           
            

            // This formulation is developed using hours as temporal variable
            real_time = real_time/3600.0;
        
            // Computing the value of heat flux according the time
            //double value = mDensity*mSpecificHeat*mAlpha*mTMax*(exp(-mAlpha*time));
            double value = 2400*1*0.025*18*(exp(-0.025*real_time));
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

    double mH0=1.0;
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