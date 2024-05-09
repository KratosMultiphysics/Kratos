//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//
//

// System includes
#include <vector>

// External includes

// Project includes
#include "containers/model.h"
#include "includes/parallel_environment.h"
#include "utilities/global_pointer_utilities.h"
#include "utilities/pointer_communicator.h"
#include "utilities/pointer_map_communicator.h"
#include "utilities/retrieve_global_pointers_by_index_functor.h"
#include "utilities/get_value_functor.h"
#include "testing/testing.h"

namespace Kratos::Testing
{

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(PointerCommunicator, KratosMPICoreFastSuite)
{
    DataCommunicator& r_default_comm = ParallelEnvironment::GetDefaultDataCommunicator();
    Model current_model;
    auto& mp = current_model.CreateModelPart("mp");
    mp.AddNodalSolutionStepVariable(PARTITION_INDEX);
    mp.AddNodalSolutionStepVariable(TEMPERATURE);

    const int world_size = r_default_comm.Size();
    const int current_rank = r_default_comm.Rank();

    auto pnode = mp.CreateNewNode(current_rank+1, current_rank,current_rank,current_rank); //the node is equal to the current rank;
    pnode->FastGetSolutionStepValue(PARTITION_INDEX) = current_rank;
    pnode->SetValue(TEMPERATURE, current_rank );

    //we will gather on every node the global pointers of the nodes with index from
    //current_rank(+1) to world_size
    std::vector<int> indices;
    indices.reserve(world_size - current_rank);
    for(int i=current_rank+1; i<=world_size; ++i)
        indices.push_back(i);

    auto gp_list = GlobalPointerUtilities::RetrieveGlobalIndexedPointers(mp.Nodes(), indices, r_default_comm );

    GlobalPointerCommunicator< Node> pointer_comm(r_default_comm, gp_list.ptr_begin(), gp_list.ptr_end());

    auto double_proxy = pointer_comm.Apply(
        [](GlobalPointer< Node >& gp)->double
        {return gp->GetValue(TEMPERATURE);}
    );

    for(unsigned int i=0; i<gp_list.size(); ++i) {
        int expected_id = indices[i];
        auto& gp = gp_list(i);
        KRATOS_EXPECT_EQ(double_proxy.Get(gp), gp.GetRank());
        KRATOS_EXPECT_EQ(double_proxy.Get(gp), expected_id-1);
    }

    //now let's try to retrieve at once TEMPERATURE, and Coordinates of the node
    typedef std::pair<double, array_1d<double,3>> return_type;

    auto pair_proxy = pointer_comm.Apply(
                          [](GlobalPointer< Node >& gp)-> return_type
    {return std::make_pair(gp->GetValue(TEMPERATURE), gp->Coordinates() );}
                      );

    for(unsigned int i=0; i<indices.size(); ++i) {
        auto& gp = gp_list(i);
        return_type result = pair_proxy.Get(gp); //this is now a pair

        KRATOS_EXPECT_EQ(result.first, gp.GetRank());

        for(unsigned int k=0; k<3; ++k)
            KRATOS_EXPECT_EQ(result.second[k], gp.GetRank());
    }
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(PointerCommunicatorPartialPartitions, KratosMPICoreFastSuite)
{
    DataCommunicator& r_default_comm = ParallelEnvironment::GetDefaultDataCommunicator();
    Model current_model;
    auto& mp = current_model.CreateModelPart("mp");
    mp.AddNodalSolutionStepVariable(PARTITION_INDEX);
    mp.AddNodalSolutionStepVariable(TEMPERATURE);

    const int world_size = r_default_comm.Size();
    const int current_rank = r_default_comm.Rank();

    auto pnode = mp.CreateNewNode(current_rank+1, current_rank,current_rank,current_rank); //the node is equal to the current rank;
    pnode->FastGetSolutionStepValue(PARTITION_INDEX) = current_rank;
    pnode->SetValue(TEMPERATURE, current_rank );

    //we will gather on every node the global pointers of the nodes with index from
    //current_rank(+1) to world_size
    std::vector<int> indices, ranks;
    indices.reserve(world_size - current_rank - 1);
    for(int i = current_rank + 1; i < world_size; ++i) {
        indices.push_back(i);
    }
    ranks.reserve(world_size - 1);
    std::string name_data_comm = "SubDataComm_";
    for(int i = 0; i < world_size - 1; ++i) {
        ranks.push_back(i);
        name_data_comm += std::to_string(i) + "_";
    }

    if (current_rank < world_size - 1) {
        auto& r_partial_data_comm = r_default_comm.GetSubDataCommunicator(ranks, name_data_comm);
        auto gp_list = GlobalPointerUtilities::RetrieveGlobalIndexedPointers(mp.Nodes(), indices, r_partial_data_comm );

        GlobalPointerCommunicator< Node> pointer_comm(r_partial_data_comm, gp_list.ptr_begin(), gp_list.ptr_end());

        auto double_proxy = pointer_comm.Apply(
            [](GlobalPointer< Node >& gp)->double
            {return gp->GetValue(TEMPERATURE);}
        );

        for(unsigned int i=0; i<gp_list.size(); ++i) {
            int expected_id = indices[i];
            auto& gp = gp_list(i);
            KRATOS_EXPECT_EQ(double_proxy.Get(gp), gp.GetRank());
            KRATOS_EXPECT_EQ(double_proxy.Get(gp), expected_id-1);
        }

        //now let's try to retrieve at once TEMPERATURE, and Coordinates of the node
        typedef std::pair<double, array_1d<double,3>> return_type;

        auto pair_proxy = pointer_comm.Apply(
                            [](GlobalPointer< Node >& gp)-> return_type
        {return std::make_pair(gp->GetValue(TEMPERATURE), gp->Coordinates() );}
                        );

        for(unsigned int i=0; i<indices.size(); ++i) {
            auto& gp = gp_list(i);
            return_type result = pair_proxy.Get(gp); //this is now a pair

            KRATOS_EXPECT_EQ(result.first, gp.GetRank());

            for(unsigned int k=0; k<3; ++k)
                KRATOS_EXPECT_EQ(result.second[k], gp.GetRank());
        }
    }

    // Extend checks
    auto& r_partial_data_comm_again = r_default_comm.GetSubDataCommunicator(ranks, name_data_comm);
    if (current_rank < world_size - 1) {
        KRATOS_EXPECT_TRUE(r_partial_data_comm_again.IsDefinedOnThisRank());
    } else {
        KRATOS_EXPECT_TRUE(r_partial_data_comm_again.IsNullOnThisRank());
    }

    std::vector<int> ranks_wrong(ranks);
    ranks_wrong.push_back(world_size -1);
    if (current_rank < world_size - 1) {
        KRATOS_EXPECT_EXCEPTION_IS_THROWN(r_default_comm.GetSubDataCommunicator(ranks_wrong, name_data_comm), "Inconsistency between the communicator world size: " + std::to_string(world_size - 1) + " and the number of ranks required: " + std::to_string(world_size));
    } else {
        KRATOS_EXPECT_EXCEPTION_IS_THROWN(r_default_comm.GetSubDataCommunicator(ranks_wrong, name_data_comm), "The rank " + std::to_string(current_rank) + " does not participate in the existing data communicator " + name_data_comm + " despite being in the provided rank list");
    }
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(PointerCommunicatorLocalRetrieveGlobalPointers, KratosMPICoreFastSuite)
{
    DataCommunicator& r_default_comm = ParallelEnvironment::GetDefaultDataCommunicator();
    Model current_model;
    auto& r_mp = current_model.CreateModelPart("mp");
    r_mp.AddNodalSolutionStepVariable(PARTITION_INDEX);
    r_mp.AddNodalSolutionStepVariable(TEMPERATURE);

    const int current_rank = r_default_comm.Rank();

    auto pnode = r_mp.CreateNewNode(current_rank+1, current_rank,current_rank,current_rank); //the node is equal to the current rank;
    pnode->FastGetSolutionStepValue(PARTITION_INDEX) = current_rank;
    pnode->SetValue(TEMPERATURE, current_rank );

    //we will gather on every node the global pointers of the nodes with index from
    // 0 to world_size
    std::vector<int> indices(1, current_rank+1);

    auto gp_list = GlobalPointerUtilities::LocalRetrieveGlobalPointers(r_mp.Nodes(), r_default_comm );
    auto gp_manual_list = GlobalPointerUtilities::RetrieveGlobalIndexedPointers(r_mp.Nodes(), indices, r_default_comm );

    GlobalPointerCommunicator< Node> pointer_comm(r_default_comm, gp_list.ptr_begin(), gp_list.ptr_end());GlobalPointerCommunicator< Node> pointer_comm_manual(r_default_comm, gp_manual_list.ptr_begin(), gp_manual_list.ptr_end());

    auto double_proxy = pointer_comm.Apply(
        [](GlobalPointer< Node >& gp)->double
        {return gp->GetValue(TEMPERATURE);}
    );

    for(unsigned int i=0; i<gp_list.size(); ++i) {
        auto& gp = gp_list(i);
        auto& gp_manual = gp_manual_list(i);
        KRATOS_EXPECT_EQ(double_proxy.Get(gp), gp.GetRank());
        KRATOS_EXPECT_EQ(gp_manual.GetRank(), gp.GetRank());
    }

    //now let's try to retrieve at once TEMPERATURE, and Coordinates of the node
    using return_type = std::pair<double, array_1d<double,3>>;

    auto lambda = [](GlobalPointer< Node >& gp)-> return_type
    {return std::make_pair(gp->GetValue(TEMPERATURE), gp->Coordinates() );};
    auto pair_proxy = pointer_comm.Apply(lambda);
    auto manual_pair_proxy = pointer_comm_manual.Apply(lambda);

    for(unsigned int i=0; i<gp_list.size(); ++i) {
        auto& gp = gp_list(i);
        auto& gp_manual = gp_manual_list(i);
        return_type result = pair_proxy.Get(gp); //this is now a pair
        return_type manual_result = manual_pair_proxy.Get(gp_manual); //this is now a pair

        KRATOS_EXPECT_EQ(result.first, gp.GetRank());
        KRATOS_EXPECT_EQ(result.first, manual_result.first);

        for(unsigned int k=0; k<3; ++k) {
            KRATOS_EXPECT_EQ(result.second[k], gp.GetRank());
            KRATOS_EXPECT_EQ(result.second[k], manual_result.second[k]);
        }
    }
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(PointerCommunicatorGlobalRetrieveGlobalPointers, KratosMPICoreFastSuite)
{
    DataCommunicator& r_default_comm = ParallelEnvironment::GetDefaultDataCommunicator();
    Model current_model;
    auto& r_mp = current_model.CreateModelPart("mp");
    r_mp.AddNodalSolutionStepVariable(PARTITION_INDEX);
    r_mp.AddNodalSolutionStepVariable(TEMPERATURE);

    const int current_rank = r_default_comm.Rank();

    auto pnode = r_mp.CreateNewNode(current_rank+1, current_rank,current_rank,current_rank); //the node is equal to the current rank;
    pnode->FastGetSolutionStepValue(PARTITION_INDEX) = current_rank;
    pnode->SetValue(TEMPERATURE, current_rank );

    //we will gather on every node the global pointers of the nodes with index from
    // 0 to world_size
    const int world_size = r_default_comm.Size();
    std::vector<int> indices;
    for(int i=0; i<world_size; ++i) {
        indices.push_back(i + 1);
    }

    auto gp_list = GlobalPointerUtilities::GlobalRetrieveGlobalPointers(r_mp.Nodes(), r_default_comm );
    auto gp_manual_list = GlobalPointerUtilities::RetrieveGlobalIndexedPointers(r_mp.Nodes(), indices, r_default_comm );

    GlobalPointerCommunicator< Node> pointer_comm(r_default_comm, gp_list.ptr_begin(), gp_list.ptr_end());GlobalPointerCommunicator< Node> pointer_comm_manual(r_default_comm, gp_manual_list.ptr_begin(), gp_manual_list.ptr_end());

    auto double_proxy = pointer_comm.Apply(
        [](GlobalPointer< Node >& gp)->double
        {return gp->GetValue(TEMPERATURE);}
    );

    for(unsigned int i=0; i<gp_list.size(); ++i) {
        auto& gp = gp_list(i);
        auto& gp_manual = gp_manual_list(i);
        KRATOS_EXPECT_EQ(double_proxy.Get(gp), gp.GetRank());
        KRATOS_EXPECT_EQ(gp_manual.GetRank(), gp.GetRank());
    }

    //now let's try to retrieve at once TEMPERATURE, and Coordinates of the node
    using return_type = std::pair<double, array_1d<double,3>>;

    auto lambda = [](GlobalPointer< Node >& gp)-> return_type
    {return std::make_pair(gp->GetValue(TEMPERATURE), gp->Coordinates() );};
    auto pair_proxy = pointer_comm.Apply(lambda);
    auto manual_pair_proxy = pointer_comm_manual.Apply(lambda);

    for(unsigned int i=0; i<gp_list.size(); ++i) {
        auto& gp = gp_list(i);
        auto& gp_manual = gp_manual_list(i);
        return_type result = pair_proxy.Get(gp); //this is now a pair
        return_type manual_result = manual_pair_proxy.Get(gp_manual); //this is now a pair

        KRATOS_EXPECT_EQ(result.first, gp.GetRank());
        KRATOS_EXPECT_EQ(result.first, manual_result.first);

        for(unsigned int k=0; k<3; ++k) {
            KRATOS_EXPECT_EQ(result.second[k], gp.GetRank());
            KRATOS_EXPECT_EQ(result.second[k], manual_result.second[k]);
        }
    }
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(PointerCommunicatorIndexConsistence, KratosMPICoreFastSuite)
{
    DataCommunicator& r_default_comm = ParallelEnvironment::GetDefaultDataCommunicator();
    Model current_model;
    auto& mp = current_model.CreateModelPart("mp");
    mp.AddNodalSolutionStepVariable(PARTITION_INDEX);

    const int world_size = r_default_comm.Size();
    const int current_rank = r_default_comm.Rank();

    // Add 3 Nodes per partition (XYZ) and make X_Z interfaces of the processes on the other side.
    // Nodes X and Z of processes 0 and RANK-1 do not communicate with anyone.
    //
    //      0      1      2    .....    N
    //    -----  -----  -----         -----
    //    X Y Z  X Y Z  X Y Z         X Y Z
    // ID 0 1 2--2 3 4--4 5 6-       -M N K
    // PI 0 0 1--1 1 2--2 2 3-       -N N N
    for(int i = 0; i < 3; i++) {
        int node_id = i + (current_rank * 2);
        auto pnode = mp.CreateNewNode(node_id, current_rank, current_rank, current_rank); //the node is equal to the current rank;
        pnode->SetValue(TEMPERATURE, current_rank );

        int partition = (i != 2) ? current_rank : std::min(current_rank+1,world_size-1);
        pnode->FastGetSolutionStepValue(PARTITION_INDEX) = partition;
    }

    // Build the list
    std::vector<int> indices(3, 0);
    for(int i = 0; i < 3; i++) {
        indices[i] = i + (current_rank * 2);
    }

    auto gp_list = GlobalPointerUtilities::RetrieveGlobalIndexedPointers(mp.Nodes(), indices, r_default_comm );

    GlobalPointerCommunicator< Node> pointer_comm(r_default_comm, gp_list.ptr_begin(), gp_list.ptr_end());

    auto double_proxy = pointer_comm.Apply(
        [](GlobalPointer< Node >& gp)->double {
            return gp->GetSolutionStepValue(PARTITION_INDEX);
        }
    );

    for(unsigned int i=0; i<indices.size(); ++i) {
        auto& gp = gp_list(i);
        KRATOS_EXPECT_EQ(double_proxy.Get(gp), gp.GetRank());
    }
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(PointerCommunicatorConstructByFunctor, KratosMPICoreFastSuite)
{
    DataCommunicator& r_default_comm = ParallelEnvironment::GetDefaultDataCommunicator();
    Model current_model;
    auto& mp = current_model.CreateModelPart("mp");
    mp.AddNodalSolutionStepVariable(PARTITION_INDEX);
    mp.AddNodalSolutionStepVariable(TEMPERATURE);

    const int world_size = r_default_comm.Size();
    const int current_rank = r_default_comm.Rank();

    auto pnode = mp.CreateNewNode(current_rank+1, current_rank,current_rank,current_rank); //the node is equal to the current rank;
    pnode->FastGetSolutionStepValue(PARTITION_INDEX) = current_rank;
    pnode->SetValue(TEMPERATURE, current_rank );

    //we will gather on every node the global pointers of the nodes with index from
    //current_rank(+1) to world_size
    std::vector<int> indices;
    for(int i=current_rank+1; i<=world_size; ++i)
        indices.push_back(i);

    //here i obtain the reference list
    auto gp_list_reference = GlobalPointerUtilities::RetrieveGlobalIndexedPointers(mp.Nodes(), indices, r_default_comm );


    //************************** VERSION WITH MINIMAL BOILERPLATE ***********************************
    //define the pointer communicator - note that a functor is employed here, so to avoid boilerplate
    GlobalPointerCommunicator< Node> pointer_comm(
                                            r_default_comm,
                                            RetrieveGlobalPointersByIndex<ModelPart::NodesContainerType>(mp.Nodes(), indices)  //if we eventually evolve to c++17 the type will be deduced
                                        );

    auto temperature_proxy = pointer_comm.Apply( GetValueFunctor<Variable<double>>(TEMPERATURE) ); //if we eventually evolve to c++17 the template parameter in the functor will be deduced

    for(unsigned int i=0; i<indices.size(); ++i)
    {
        auto& gp = gp_list_reference(i);
        KRATOS_EXPECT_EQ(temperature_proxy.Get(gp), gp.GetRank());
        KRATOS_EXPECT_EQ(temperature_proxy.Get(gp), indices[i]-1);
    }
    //**********************************************************************************************


    //just checking that update works (check is simply repeated after update)
    temperature_proxy.Update();
    for(unsigned int i=0; i<indices.size(); ++i)
    {
        auto& gp = gp_list_reference(i);
        KRATOS_EXPECT_EQ(temperature_proxy.Get(gp), gp.GetRank());
        KRATOS_EXPECT_EQ(temperature_proxy.Get(gp), indices[i]-1);
    }
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(PointerMapCommunicatorAssembly, KratosMPICoreFastSuite)
{
    DataCommunicator& r_default_comm = ParallelEnvironment::GetDefaultDataCommunicator();
    Model current_model;
    auto& mp = current_model.CreateModelPart("mp");
    mp.AddNodalSolutionStepVariable(PARTITION_INDEX);

    const int world_size = r_default_comm.Size();
    const int current_rank = r_default_comm.Rank();

    auto pnode = mp.CreateNewNode(current_rank+1, current_rank,current_rank,current_rank); //the node is equal to the current rank;
    pnode->FastGetSolutionStepValue(PARTITION_INDEX) = current_rank;
    pnode->SetValue(TEMPERATURE, current_rank);
    pnode->SetValue(VELOCITY, array_1d<double, 3>(3, current_rank));

    // creates the gp map
    GlobalPointerMapCommunicator<Node, double> pointer_map_double_comm(r_default_comm);
    GlobalPointerMapCommunicator<Node, array_1d<double, 3>> pointer_map_array3_comm(r_default_comm);

    // creates apply proxy
    auto apply_temperature_assemble_proxy = pointer_map_double_comm.GetApplyProxy(
        [](Node& rNode, const double& NewValue) {
            rNode.GetValue(TEMPERATURE) += NewValue;
        });

    auto apply_velocity_assemble_proxy = pointer_map_array3_comm.GetApplyProxy(
        [](Node& rNode, const array_1d<double, 3>& NewValue) {
            rNode.GetValue(VELOCITY) += NewValue;
        });

    //we will gather on every node the global pointers of the nodes with index from
    //current_rank(+1) to world_size
    std::vector<int> indices;
    for(int i=1; i<=world_size; ++i) {
        indices.push_back(i);
    }
    auto gp_list = GlobalPointerUtilities::RetrieveGlobalIndexedPointers(mp.Nodes(), indices, r_default_comm);

    for (int i = 0; i < current_rank + 1; ++i) {
        auto gp = gp_list(i);
        const double value = (gp.GetRank() + 1.0) * 2;
        apply_temperature_assemble_proxy.Assign(gp, value * 3);
        apply_velocity_assemble_proxy.Assign(gp, array_1d<double, 3>(3, value * 3));
    }

    // Assigns local gps and stores remote gps in a map for future communication
    apply_temperature_assemble_proxy.Assign(gp_list(0), 1);
    apply_temperature_assemble_proxy.Assign(gp_list(0), 2);

    apply_velocity_assemble_proxy.Assign(gp_list(0), array_1d<double, 3>(3, 1));
    apply_velocity_assemble_proxy.Assign(gp_list(0), array_1d<double, 3>(3, 2));

    // communicate remote gps and update them
    apply_temperature_assemble_proxy.SendAndApplyRemotely();
    apply_velocity_assemble_proxy.SendAndApplyRemotely();

    const double check_value =
        (current_rank == 0.0)
            ? world_size * 6 + world_size * 3
            : current_rank + (world_size - current_rank) * (current_rank + 1) * 6;

    KRATOS_EXPECT_EQ(pnode->GetValue(TEMPERATURE), check_value);
    const array_1d<double, 3> check_vector(3, check_value);
    KRATOS_EXPECT_VECTOR_EQ(pnode->GetValue(VELOCITY), check_vector);
}

}