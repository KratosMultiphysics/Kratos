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

// External includes
#include "mpi.h"

// Project includes
#include "element_communicator_mpi.h"
#include "particle_mechanics_application_variables.h"
#include "../mpm_search_element_utility.h"
#include "utilities/variable_utils.h"
#include "utilities/parallel_utilities.h"

// Conditions
#include "../../custom_conditions/particle_based_conditions/mpm_particle_penalty_dirichlet_condition.h"
#include "../../custom_conditions/particle_based_conditions/mpm_particle_point_load_condition.h"


namespace Kratos {

struct ConditionMessage{
    bool empty = 0; // 0 - not empty, 1 - empty
    int type = 0;   // 0 - MPMParticlePointLoadCondition2D3N
                    // 1 - MPMParticlePenaltyDirichletCondition
    int mpi_send_counter = 0;
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

void ElementCommunicatorMPI::CreateMPIDataType(MPI_Datatype& type, bool condition){
    if(condition){
        const unsigned int number_of_elements = 14;
        MPI_Aint     indices[number_of_elements];
        MPI_Datatype old_types[number_of_elements];
        int          blocklens[number_of_elements];
        // Data types
        old_types[0] = MPI_INT;
        old_types[1] = MPI_INT;
        old_types[2] = MPI_INT;
        for(int i = 3; i < number_of_elements;i++) {old_types[i] = MPI_DOUBLE;}
        // Number of elements per date type
        blocklens[0] = 1;
        blocklens[1] = 1;
        blocklens[2] = 1;
        for(int i = 3; i < number_of_elements;i++) {blocklens[i] = 3;}
        indices[0] = 0;
        indices[1] = sizeof(int);
        indices[2] = 2*sizeof(int);
        for(int i = 3; i < number_of_elements; i++) {indices[i] = 3*sizeof(int) + (i-3) * 3 *sizeof(double);}

        MPI_Type_create_struct(number_of_elements,blocklens,indices,old_types,&type);
        MPI_Type_commit(&type);
    }
    else {
        KRATOS_INFO("ElementCommunicator") << "Element 'MPI-Type' is not yet implemented!" << std::endl;
    }
}

void ElementCommunicatorMPI::MPI_Search(ModelPart& rMPMModelPart,
                                        ModelPart& rBackgroundGridModelPart,
                                        std::vector<typename Condition::Pointer>& rMissingConditions,
                                        const std::size_t MaxNumberOfResults, const double Tolerance){
    // Get current rank and size
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
    CreateMPIDataType(message_type, true); // true -> condition, false -> element/particle

    // Construct local arrays
    std::vector<typename Element::Pointer> missing_elements = {};
    std::vector<typename Condition::Pointer> missing_conditions = {};

    // Check if rMissingConditions is empty or not
    bool empty = false;
    if( rMissingConditions.size() == 0){
        empty = true;
    }
    else {
        #pragma omp parallel for //TODO: Use lambda function here
        for( auto& condition : rMissingConditions ){
            condition->SetValue(MPI_SEND_COUNTER,0);
        }
    }

    // Sent missing condition to next mpi node (in a circle)
    // Run BinBasedSearch for received condition.
    // If condition can not be found, sent to next mpi node
    // To keep pattern consistent: If there are no missing conditions, dummy message is sent.
    // TODO: Implement counter to avoid reallocating of rMissingCondtions
    while( rMPMModelPart.GetCommunicator().GetDataCommunicator().SumAll(empty) !=  size){
        if( rMissingConditions.size() == 0){
            empty = true;
            Condition::Pointer dummy_condition_pointer;
            rMissingConditions.push_back(dummy_condition_pointer);
        }
        // Send last condition in rMissingConditions and remove it from vector. If empty==true send empty message.
        ElementCommunicatorMPI::SentCondition(rMPMModelPart, rMissingConditions.back(), empty, message_type, receiver, 0);
        rMissingConditions.pop_back();

        // Receive condition
        bool recieved_message_is_empty;
        unsigned int condition_send_counter = ElementCommunicatorMPI::RecieveCondition(rMPMModelPart,
            rBackgroundGridModelPart, missing_conditions, recieved_message_is_empty, message_type, sender, 0);

        // Only for debugging
        int current_id;
        // Run BinBasedSearch on received condition
        if( !recieved_message_is_empty){
            current_id = missing_conditions[0]->Id();
            MPMSearchElementUtility::BinBasedSearchElementsAndConditions<2>(rMPMModelPart,
                rBackgroundGridModelPart, missing_elements, missing_conditions,
                    MaxNumberOfResults, Tolerance);
        }
        // Only for debugging
        if( missing_conditions.size() == 0 && !recieved_message_is_empty){
            // std::cout << "Rank: " << rank << " found condition " << current_id
            //     << " after " << condition_send_counter << " mpi send operations!!" << std::endl;
        }
        // If condition was not found add to vector
        if( missing_conditions.size() != 0){
            if( condition_send_counter == size - 1){
                std::cout << "hallo" << std::endl;
                KRATOS_INFO_ALL_RANKS("ElementMPICommunicator") << "Warning: MPI Search for condition: " << missing_conditions[0]->Id() << " failed."
                    << "Geometry is cleared" << std::endl;
                // Delete condition
                missing_conditions[0]->GetGeometry().clear();
                missing_conditions[0]->Reset(ACTIVE);
                missing_conditions[0]->Set(TO_ERASE);
                missing_conditions.clear();
            }
            else{
                rMissingConditions.push_back(missing_conditions[0]);
                missing_conditions.clear();
                empty = false;
            }
        }
    }
}

void ElementCommunicatorMPI::SentCondition(ModelPart& rMPMModelPart, Condition::Pointer pCondition, bool empty, MPI_Datatype& message_type, int destination, int tag){
        // Get all general condition information and fill mpi data struct
        ProcessInfo process_info = ProcessInfo();
        std::string submodelpart_name = "";
        ConditionMessage message;
        message.empty = empty;
        if( !empty ){
            // Get general variables
            std::vector<array_1d<double, 3>> xg;
            std::vector<array_1d<double, 3>> mpc_normal;
            std::vector<double> mpc_area;
            pCondition->CalculateOnIntegrationPoints(MPC_COORD, xg, process_info);
            message.xg = xg[0];
            pCondition->CalculateOnIntegrationPoints(MPC_NORMAL,mpc_normal, process_info);
            message.mpc_normal = mpc_normal[0];
            pCondition->CalculateOnIntegrationPoints(MPC_AREA,mpc_area, process_info);
            message.mpc_area[0] = mpc_area[0];
            message.mpc_area[1] = mpc_area[1];
            message.mpc_area[2] = mpc_area[2];
            message.mpi_send_counter = pCondition->GetValue(MPI_SEND_COUNTER) + 1;
            // Get condition specific variables
            std::string  condition_type_name;
            std::vector<array_1d<double, 3>> point_load;
            std::vector<double> mpc_penalty_factor_vector;
            std::vector<array_1d<double,3>> mpc_displacement;
            std::vector<array_1d<double,3>> mpc_imposed_displacement;
            std::vector<array_1d<double,3>> mpc_velocity;
            std::vector<array_1d<double,3>> mpc_imposed_velocity;
            std::vector<array_1d<double,3>> mpc_acceleration;
            std::vector<array_1d<double,3>> mpc_imposed_acceleration;
            // TODO: Implement missing condition types (Only particle conditions are required!!!)
            // If condition is MPMParticlePointLoadCondition2D3N
            if( typeid(*pCondition) == typeid(MPMParticlePointLoadCondition) ){
                condition_type_name = "MPMParticlePointLoadCondition2D3N";
                message.type = 0;
                pCondition->CalculateOnIntegrationPoints(POINT_LOAD,point_load, process_info);
                message.point_load = point_load[0];
            }
            // If condition is MPMParticlePenaltyDirichletCondition2D3N
            if( typeid(*pCondition) == typeid(MPMParticlePenaltyDirichletCondition) ){
                condition_type_name = "MPMParticlePenaltyDirichletCondition2D3N";
                message.type = 1;
                pCondition->CalculateOnIntegrationPoints(PENALTY_FACTOR,mpc_penalty_factor_vector, process_info);
                message.mpc_penalty_factor[0] = mpc_penalty_factor_vector[0];
                message.mpc_penalty_factor[1] = mpc_penalty_factor_vector[1];
                message.mpc_penalty_factor[2] = mpc_penalty_factor_vector[2];
                pCondition->CalculateOnIntegrationPoints(MPC_DISPLACEMENT, mpc_displacement, process_info);
                message.mpc_displacement = mpc_displacement[0];
                pCondition->CalculateOnIntegrationPoints(MPC_IMPOSED_DISPLACEMENT, mpc_imposed_displacement, process_info);
                message.mpc_imposed_displacement = mpc_imposed_displacement[0];
                pCondition->CalculateOnIntegrationPoints(MPC_VELOCITY, mpc_velocity, process_info);
                message.mpc_velocity = mpc_velocity[0];
                pCondition->CalculateOnIntegrationPoints(MPC_IMPOSED_VELOCITY, mpc_imposed_velocity, process_info);
                message.mpc_imposed_velocity = mpc_imposed_velocity[0];
                pCondition->CalculateOnIntegrationPoints(MPC_ACCELERATION, mpc_acceleration, process_info);
                message.mpc_acceleration = mpc_acceleration[0];
                pCondition->CalculateOnIntegrationPoints(MPC_IMPOSED_ACCELERATION, mpc_imposed_acceleration, process_info);
                message.mpc_imposed_acceleration = mpc_imposed_acceleration[0];
            }

            // Find model part name
            // TODO: Is there any better way to find the submodelpart name?
            for( auto& submodelpart : rMPMModelPart.SubModelParts()){
                for( auto& condition : submodelpart.Conditions()){
                    if( condition.Id() == pCondition->Id()){
                        submodelpart_name = submodelpart.Name();
                        break;
                    }
                }
            }
            if( submodelpart_name == ""){
                KRATOS_ERROR << "ElementCommunicatorMPI: " << "Submodelpart of condition ID: " << pCondition->Id()
                    << " not found!";
            }
        }
        // Send name of submodelpart that contains the condition
        MPI_Send(submodelpart_name.c_str(), submodelpart_name.size(), MPI_CHAR, destination, 1, MPI_COMM_WORLD);
        // Send condition information
        MPI_Send(&message, 1, message_type, destination, 0, MPI_COMM_WORLD);

        // If a non empty conditon was sent
        if( !empty ){
            // Remove condition from model part
            rMPMModelPart.GetSubModelPart(submodelpart_name).RemoveCondition(pCondition);
            rMPMModelPart.GetSubModelPart(submodelpart_name).GetCommunicator().LocalMesh().RemoveCondition(pCondition);
            rMPMModelPart.RemoveCondition(pCondition);
            // Delete condition
            pCondition->GetGeometry().clear();
            pCondition->Reset(ACTIVE);
            pCondition->Set(TO_ERASE);
        }
    }

    int ElementCommunicatorMPI::RecieveCondition(ModelPart& rMPMModelPart,
                                                 ModelPart& rBackgroundGridModelPart,
                                                 std::vector<typename Condition::Pointer>& rMissingCondition,
                                                 bool& empty, MPI_Datatype& message_type, int source, int tag){
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
        ConditionMessage message;
        MPI_Recv(&message,1,message_type,source,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        empty = message.empty;
        // If message is not empty
        if( !message.empty){
            // Check for the next free condition Id
            // TODO: Check if this is correct!
            unsigned int new_condition_id;
            if(rMPMModelPart.NumberOfConditions() != 0){
                const auto condition_it_begin = rMPMModelPart.Conditions().begin();
                unsigned int tmp = condition_it_begin->Id();
                for( int i = 1; i < rMPMModelPart.Conditions().size(); ++i ){
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
                const unsigned int number_elements = rMPMModelPart.NumberOfElements();
                const unsigned int number_nodes = rMPMModelPart.NumberOfNodes();
                const unsigned int number_conditions = 0;
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
            p_condition->SetValue(MPI_SEND_COUNTER,message.mpi_send_counter);

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
        return message.mpi_send_counter;
    }

} // namespace Kratos