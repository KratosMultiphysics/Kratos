//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Manuel Messmer
//

#include "mpi.h"

// Project includes
#include "element_communicator_mpi.h"
#include "includes/parallel_environment.h"
#include "particle_mechanics_application_variables.h"
#include "../../custom_conditions/particle_based_conditions/mpm_particle_penalty_dirichlet_condition.h"
#include "../../custom_conditions/particle_based_conditions/mpm_particle_point_load_condition.h"
#include "../mpm_search_element_utility.h"



namespace Kratos {

struct message_struct{
    bool empty = 0; // 0 - not empty, 1 - empty
    bool type = 0;      // 0 - MPMParticlePointLoadCondition2D3N
                        // 1 - MPMParticlePenaltyDirichletCondition
    // General variables
    array_1d<double, 3> xg = ZeroVector(3);
    array_1d<double, 3> mpc_normal = ZeroVector(3);
    array_1d<double, 3> mpc_area = ZeroVector(3);
    // Neumann related variables
    array_1d<double, 3> point_load = ZeroVector(3);
    // Dirichlet related variables
    array_1d<double, 3> mpc_penalty_factor = ZeroVector(3);
    array_1d<double, 3> mpc_displacement = ZeroVector(3);
    array_1d<double, 3> mpc_imposed_displacement = ZeroVector(3);
    array_1d<double, 3> mpc_velocity = ZeroVector(3);
    array_1d<double, 3> mpc_imposed_velocity = ZeroVector(3);
    array_1d<double, 3> mpc_acceleration = ZeroVector(3);
    array_1d<double, 3> mpc_imposed_acceleration = ZeroVector(3);

};

void ElementCommunicatorMPI::MPI_InitialSearch(ModelPart& rMPMModelPart, ModelPart& rBackgroundGridModelPart,  std::vector<typename Condition::Pointer>& rMissingConditions){
    int rank;
    int size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Send message in circle
    int receiver = (rank + 1) % size;
    int sender = (rank - 1);
    if( sender == -1) {sender = size-1;}
    // Construct MPI message_type
    MPI_Datatype message_type;
    MPI_Aint     indices[13];
    MPI_Datatype old_types[13];
    int          blocklens[13];
    // Data types
    old_types[0] = MPI_INT;
    old_types[1] = MPI_INT;
    for(int i = 2; i < 13;i++) {old_types[i] = MPI_DOUBLE;}
    // Number of elements per date type
    blocklens[0] = 1;
    blocklens[1] = 1;
    for(int i = 2; i < 13;i++) {blocklens[i] = 3;}
    indices[0] = 0;
    indices[1] = sizeof(int);
    for(int i = 2; i < 13; i++) {indices[i] = 2*sizeof(int) + (i-2) * 3 *sizeof(double);}

    MPI_Type_create_struct(13,blocklens,indices,old_types,&message_type);
    MPI_Type_commit(&message_type);

    std::vector<typename Element::Pointer> missing_elements = {};
    std::vector<typename Condition::Pointer> missing_conditions = {};

    bool empty = false;
    if( rMissingConditions.size() == 0){
        empty = true;
    }
    // Sent missing condition to next mpi node (in a circle)
    // To keep pattern consistent. If there are no missing conditions, dummy message is sent
    while( rMPMModelPart.GetCommunicator().GetDataCommunicator().SumAll(empty) !=  size){
        if( rMissingConditions.size() == 0){
            empty = true;
            Condition::Pointer dummy_condition_pointer;
            rMissingConditions.push_back(dummy_condition_pointer);
        }
        ElementCommunicatorMPI::SentCondition(rMPMModelPart, rMissingConditions[0], empty, message_type, receiver, 0);
        rMissingConditions.erase(rMissingConditions.begin());

        ElementCommunicatorMPI::RecieveCondition(rMPMModelPart, rBackgroundGridModelPart, missing_conditions, message_type, sender, 0);

        // Todo: Fill the here hard-coded variables
        MPMSearchElementUtility::BinBasedSearchElementsAndConditions<2>(rMPMModelPart,
            rBackgroundGridModelPart, missing_elements, missing_conditions,
                1000, 0.000000001);

        if( missing_conditions.size() != 0){
            rMissingConditions.push_back(missing_conditions[0]);
            missing_conditions.clear();
            empty = false;
        }
    }
}

void ElementCommunicatorMPI::SentCondition(ModelPart& rMPMModelPart, Condition::Pointer cond, bool empty, MPI_Datatype& message_type, int destination, int tag){
        // Get all general condition information
        ProcessInfo process_info = ProcessInfo();
        std::vector<array_1d<double, 3>> xg;
        std::vector<array_1d<double, 3>> mpc_normal;
        std::vector<double> mpc_area;
        std::string submodelpart_name = "";
        message_struct message;
        message.empty = empty;
        if( !empty ){
            cond->CalculateOnIntegrationPoints(MPC_COORD, xg, process_info);
            cond->CalculateOnIntegrationPoints(MPC_NORMAL,mpc_normal, process_info);
            cond->CalculateOnIntegrationPoints(MPC_AREA,mpc_area, process_info);
            message.xg = xg[0];
            message.mpc_normal = mpc_normal[0];
            message.mpc_area[0] = mpc_area[0];
            message.mpc_area[1] = mpc_area[1];
            message.mpc_area[2] = mpc_area[2];

            // Get condition type
            std::string  condition_type_name;
            std::vector<array_1d<double, 3>> point_load;
            std::vector<double> mpc_penalty_factor_vector;
            std::vector<array_1d<double,3>> mpc_displacement;
            std::vector<array_1d<double,3>> mpc_imposed_displacement;
            std::vector<array_1d<double,3>> mpc_velocity;
            std::vector<array_1d<double,3>> mpc_imposed_velocity;
            std::vector<array_1d<double,3>> mpc_acceleration;
            std::vector<array_1d<double,3>> mpc_imposed_acceleration;

            // If condition is MPMParticlePointLoadCondition2D3N
            if( typeid(*cond) == typeid(MPMParticlePointLoadCondition) ){
                condition_type_name = "MPMParticlePointLoadCondition2D3N";
                message.type = 0;
                cond->CalculateOnIntegrationPoints(POINT_LOAD,point_load, process_info);
                message.point_load = point_load[0];
            }
            if( typeid(*cond) == typeid(MPMParticlePenaltyDirichletCondition) ){
                condition_type_name = "MPMParticlePenaltyDirichletCondition2D3N";
                message.type = 1;
                cond->CalculateOnIntegrationPoints(PENALTY_FACTOR,mpc_penalty_factor_vector, process_info);
                message.mpc_penalty_factor[0] = mpc_penalty_factor_vector[0];
                message.mpc_penalty_factor[1] = mpc_penalty_factor_vector[1];
                message.mpc_penalty_factor[2] = mpc_penalty_factor_vector[2];
                cond->CalculateOnIntegrationPoints(MPC_DISPLACEMENT, mpc_displacement, process_info);
                message.mpc_displacement = mpc_displacement[0];
                cond->CalculateOnIntegrationPoints(MPC_IMPOSED_DISPLACEMENT, mpc_imposed_displacement, process_info);
                message.mpc_imposed_displacement = mpc_imposed_displacement[0];
                cond->CalculateOnIntegrationPoints(MPC_VELOCITY, mpc_velocity, process_info);
                message.mpc_velocity = mpc_velocity[0];
                cond->CalculateOnIntegrationPoints(MPC_IMPOSED_VELOCITY, mpc_imposed_velocity, process_info);
                message.mpc_imposed_velocity = mpc_imposed_velocity[0];
                cond->CalculateOnIntegrationPoints(MPC_ACCELERATION, mpc_acceleration, process_info);
                message.mpc_acceleration = mpc_acceleration[0];
                cond->CalculateOnIntegrationPoints(MPC_IMPOSED_ACCELERATION, mpc_imposed_acceleration, process_info);
                message.mpc_imposed_acceleration = mpc_imposed_acceleration[0];

            }

            for( auto& submodelpart : rMPMModelPart.SubModelParts()){
                for( auto& condition : submodelpart.Conditions()){
                    if( condition.Id() == cond->Id()){
                        submodelpart_name = submodelpart.Name();
                    }
                }
            }
            if( submodelpart_name == ""){
                KRATOS_ERROR << "Submodelpart not found!";
            }
        }
        MPI_Send(submodelpart_name.c_str(), submodelpart_name.size(), MPI_CHAR, destination, 1, MPI_COMM_WORLD);
        MPI_Send(&message,1,message_type,destination,0,MPI_COMM_WORLD);
        if( !empty ){
            // Remove condition from model part
            rMPMModelPart.GetSubModelPart(submodelpart_name).RemoveCondition(cond);
            rMPMModelPart.GetSubModelPart(submodelpart_name).GetCommunicator().LocalMesh().RemoveCondition(cond);
            rMPMModelPart.RemoveCondition(cond);
            // Todo: Check if this is required
            cond->GetGeometry().clear();
            cond->Reset(ACTIVE);
            cond->Set(TO_ERASE);
        }
    }

    void ElementCommunicatorMPI::RecieveCondition(ModelPart& rMPMModelPart, ModelPart& rBackgroundGridModelPart, std::vector<typename Condition::Pointer>& rMissingCondition, MPI_Datatype& message_type, int source, int tag){
        // Receive modelpart name
        MPI_Status status;
        MPI_Probe(source,1,MPI_COMM_WORLD,&status);
        int length;
        MPI_Get_count(&status,MPI_CHAR,&length);
        char* number_buf = (char*)malloc(sizeof(char) * length);
        MPI_Recv(number_buf, length,MPI_CHAR,source,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        std::string submodelpart_name(number_buf,length);
        free(number_buf);

        // Receive condition imformation
        message_struct message;
        MPI_Recv(&message,1,message_type,source,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        // If message is not empty
        if( !message.empty){
            // Check for the next free condition Id
            // TODO: Check if this is correct and do better!
            int new_condition_id;
            if(rMPMModelPart.NumberOfConditions() != 0){
                const auto condition_it_begin = rMPMModelPart.Conditions().begin();
                int tmp = condition_it_begin->Id();
                for( int i = 1; i < static_cast<int>(rMPMModelPart.Conditions().size()); ++i ){
                    auto condition_it = condition_it_begin + i;
                    int prev = condition_it->Id();
                    if( prev != tmp + 1 ){
                        break;
                    }
                    tmp = prev;
                }
                new_condition_id = tmp + 1;

            }
            else{
                int number_elements = rMPMModelPart.NumberOfElements();
                int number_nodes = rMPMModelPart.NumberOfNodes();
                int number_conditions = 0;
                if (number_elements > number_nodes && number_elements > number_conditions)
                    new_condition_id = number_elements + 1;
                else if (number_nodes > number_elements && number_nodes > number_conditions)
                    new_condition_id = number_nodes + 1;
                else
                    new_condition_id = number_conditions + 1;

            }

            // Get condition type
            std::string condition_type_name;
            if( message.type == 0){
                condition_type_name = "MPMParticlePointLoadCondition2D3N";
            }
            else if(message.type ==1){
                condition_type_name = "MPMParticlePenaltyDirichletCondition2D3N";
            }

            // Create new condition
            Properties::Pointer properties;
            const Condition& new_condition = KratosComponents<Condition>::Get(condition_type_name);
            Condition::Pointer p_condition = new_condition.Create(new_condition_id, rBackgroundGridModelPart.ElementsBegin()->GetGeometry(), properties);
            // Assign general variables
            ProcessInfo process_info = ProcessInfo();
            p_condition->SetValuesOnIntegrationPoints(MPC_COORD, { message.xg }, process_info);
            std::vector<double> mpc_area_vector = { message.mpc_area[0] };
            p_condition->SetValuesOnIntegrationPoints(MPC_AREA, mpc_area_vector, process_info);
            p_condition->SetValuesOnIntegrationPoints(MPC_NORMAL, { message.mpc_normal }, process_info);
            // Asign neumann related variables
            if( message.type == 0 ){
                p_condition->SetValuesOnIntegrationPoints(POINT_LOAD,{ message.point_load }, process_info);
            }
            // Assign dirichlet related variables
            if( message.type == 1 ){
                p_condition->SetValuesOnIntegrationPoints(MPC_DISPLACEMENT, { message.mpc_displacement }, process_info);
                p_condition->SetValuesOnIntegrationPoints(MPC_IMPOSED_DISPLACEMENT, { message.mpc_imposed_displacement }, process_info);
                p_condition->SetValuesOnIntegrationPoints(MPC_VELOCITY, { message.mpc_velocity }, process_info);
                p_condition->SetValuesOnIntegrationPoints(MPC_IMPOSED_VELOCITY, { message.mpc_imposed_velocity }, process_info);
                p_condition->SetValuesOnIntegrationPoints(MPC_ACCELERATION, { message.mpc_acceleration }, process_info);
                p_condition->SetValuesOnIntegrationPoints(MPC_IMPOSED_ACCELERATION, { message.mpc_imposed_acceleration }, process_info);
                std::vector<double> mpc_penalty_factor = { message.mpc_penalty_factor[0] };
                p_condition->SetValuesOnIntegrationPoints(PENALTY_FACTOR, mpc_penalty_factor, process_info);
            }
            // Add condition to model part
            rMissingCondition.push_back(p_condition);
            if( !rMPMModelPart.HasSubModelPart(submodelpart_name) )
                rMPMModelPart.CreateSubModelPart(submodelpart_name);

            rMPMModelPart.GetSubModelPart(submodelpart_name).AddCondition(p_condition);
            rMPMModelPart.GetSubModelPart(submodelpart_name).GetCommunicator().LocalMesh().AddCondition(p_condition);
        }
    }

} // namespace Kratos