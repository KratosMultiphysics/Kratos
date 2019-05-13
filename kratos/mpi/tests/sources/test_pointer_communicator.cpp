//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//
//
#include <vector>
#include "containers/model.h"
#include "includes/parallel_environment.h"
#include "utilities/global_pointer_utilities.h"
#include "utilities/pointer_communicator.h"

#include "testing/testing.h"

namespace Kratos
{

namespace Testing
{

KRATOS_TEST_CASE_IN_SUITE(TestPointerCommunicator, KratosMPICoreFastSuite)
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

    auto gp_list = GlobalPointerUtilities::RetrieveGlobalIndexedPointers(mp.Nodes(), indices, r_default_comm );

    GlobalPointerCommunicator< Node<3>> pointer_comm(r_default_comm, gp_list.begin(), gp_list.end());

    auto double_proxy = pointer_comm.Apply(
                            [](GlobalPointer< Node<3> >& gp)->double
    {return gp->GetValue(TEMPERATURE);}
                        );

    for(unsigned int i=0; i<indices.size(); ++i)
    {
        int expected_id = indices[i];
        auto& gp = gp_list[i];
        KRATOS_CHECK_EQUAL(double_proxy.Get(gp), gp.GetRank());
        KRATOS_CHECK_EQUAL(double_proxy.Get(gp), expected_id-1);
    }

    //now let's try to retrieve at once TEMPERATURE, and Coordinates of the node
    typedef std::pair<double, array_1d<double,3>> return_type;

    auto pair_proxy = pointer_comm.Apply(
                          [](GlobalPointer< Node<3> >& gp)-> return_type
    {return std::make_pair(gp->GetValue(TEMPERATURE), gp->Coordinates() );}
                      );

    for(unsigned int i=0; i<indices.size(); ++i)
    {
        auto& gp = gp_list[i];
        return_type result = pair_proxy.Get(gp); //this is now a pair

        KRATOS_CHECK_EQUAL(result.first, gp.GetRank());

        for(unsigned int k=0; k<3; ++k)
            KRATOS_CHECK_EQUAL(result.second[k], gp.GetRank());
    }

};

KRATOS_TEST_CASE_IN_SUITE(TestPointerCommunicatorConstructByFunctor, KratosMPICoreFastSuite)
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

    //hide the utility in a functor
    class GatherNodesByIdFunctor
    {
    public:

        GatherNodesByIdFunctor(
            const ModelPart& rModelPart,
            const std::vector<int>& Indices)
            : mrModelPart(rModelPart), mIndices(Indices)
        {}

        std::vector<GlobalPointer<Node<3>>> operator()(const DataCommunicator& rComm) const
        {
            return GlobalPointerUtilities::RetrieveGlobalIndexedPointers(mrModelPart.Nodes(), mIndices, rComm );
        }

    private:
        const ModelPart& mrModelPart;
        const std::vector<int> mIndices;
    };

    GlobalPointerCommunicator< Node<3>> pointer_comm(r_default_comm, GatherNodesByIdFunctor(mp, indices) );

    auto double_proxy = pointer_comm.Apply(
                            [](GlobalPointer< Node<3> >& gp)->double
    {return gp->GetValue(TEMPERATURE);}
                        );

    auto gp_list_reference = GatherNodesByIdFunctor(mp, indices)(r_default_comm);

    for(unsigned int i=0; i<indices.size(); ++i)
    {
        int expected_id = indices[i];
        auto& gp = gp_list_reference[i];
        KRATOS_CHECK_EQUAL(double_proxy.Get(gp), gp.GetRank());
        KRATOS_CHECK_EQUAL(double_proxy.Get(gp), expected_id-1);
    }


    //now let's try to retrieve at once TEMPERATURE, and Coordinates of the node
    //here i define a simple_functor to return the data, even though this is easily achieved by a lambda
    class ResultFunctor
    {
        public:
            std::pair<double, array_1d<double,3>> operator()(GlobalPointer< Node<3> >& gp)
            {
                return std::make_pair(gp->GetValue(TEMPERATURE), gp->Coordinates() );
            }
    };

    auto pair_proxy = pointer_comm.Apply(ResultFunctor());

    for(unsigned int i=0; i<indices.size(); ++i)
    {
        auto& gp = gp_list_reference[i];
        auto result = pair_proxy.Get(gp); //this is now a pair

        KRATOS_CHECK_EQUAL(result.first, gp.GetRank());

        for(unsigned int k=0; k<3; ++k)
            KRATOS_CHECK_EQUAL(result.second[k], gp.GetRank());
    }

    //here we check if Update works correctly
    pair_proxy.Update();
    for(unsigned int i=0; i<indices.size(); ++i)
    {
        auto& gp = gp_list_reference[i];
        auto result = pair_proxy.Get(gp); //this is now a pair

        KRATOS_CHECK_EQUAL(result.first, gp.GetRank());

        for(unsigned int k=0; k<3; ++k)
            KRATOS_CHECK_EQUAL(result.second[k], gp.GetRank());
    }   

};



}
}