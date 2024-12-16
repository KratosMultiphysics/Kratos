//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

// System includes
#include <vector>

// External includes
#include <benchmark/benchmark.h>

// Project includes
#include "containers/model.h"
#include "containers/pointer_vector_set.h"
#include "includes/indexed_object.h"
#include "includes/key_hash.h"
#include "includes/model_part.h"
#include "includes/smart_pointers.h"
#include "utilities/assign_unique_model_part_collection_tag_utility.h"

namespace Kratos {

static void BM_ModelPartCreateNewNode(benchmark::State& state) {
    for (auto _ : state) {
        const IndexType input_size = state.range(0);

        state.PauseTiming();

        auto model = Model();
        auto& r_model_part = model.CreateModelPart("test");

        state.ResumeTiming();

        for (IndexType i = 0; i < input_size; ++i) {
            r_model_part.CreateNewNode(i + 1, 0.0, 0.0, 0.0);
        }
    }
}


static void BM_ModelPartAddNodesToRootFromRange1(benchmark::State& state) {
    for (auto _ : state) {
        const IndexType input_size = state.range(0);
        const bool fill_existing = state.range(1) == 1;

        state.PauseTiming();

        auto model = Model();
        auto& r_model_part = model.CreateModelPart("test");

        std::vector<ModelPart::NodeType::Pointer> nodes_to_be_added, existing_nodes;
        if (!fill_existing) {
            nodes_to_be_added.reserve(input_size);
        } else {
            nodes_to_be_added.reserve(input_size / 2 + 1);
            existing_nodes.reserve(input_size / 2 + 1);
        }

        // doing it in reverse to make sure a proper sort is done always
        for (IndexType i = input_size; i > 0; --i) {
            auto p_node = Kratos::make_intrusive<NodeType>(i, 0.0, 0.0, 0.0);
            if (!fill_existing || i % 2 == 0) {
                nodes_to_be_added.push_back(p_node);
            } else {
                existing_nodes.push_back(p_node);
            }
        }
        r_model_part.AddNodes(existing_nodes.begin(), existing_nodes.end());
        state.ResumeTiming();

        // benchmarking to add nodes for already existing model part with nodes
        r_model_part.AddNodes(nodes_to_be_added.begin(), nodes_to_be_added.end());
    }
}

static void BM_ModelPartAddNodesToRootFromRange2(benchmark::State& state) {
    for (auto _ : state) {
        const IndexType input_size = state.range(0);
        const bool fill_existing = state.range(1) == 1;

        state.PauseTiming();

        auto model = Model();
        auto& r_model_part = model.CreateModelPart("test");

        PointerVectorSet<ModelPart::NodeType, IndexedObject> nodes_to_be_added, existing_nodes;
        if (!fill_existing) {
            nodes_to_be_added.reserve(input_size);
        } else {
            nodes_to_be_added.reserve(input_size / 2 + 1);
            existing_nodes.reserve(input_size / 2 + 1);
        }

        // doing it in reverse to make sure a proper sort is done always
        for (IndexType i = 0; i < input_size; ++i) {
            auto p_node = Kratos::make_intrusive<NodeType>(i + 1, 0.0, 0.0, 0.0);
            if (!fill_existing || i % 2 == 0) {
                nodes_to_be_added.insert(nodes_to_be_added.end(), p_node);
            } else {
                existing_nodes.insert(existing_nodes.end(), p_node);
            }
        }
        r_model_part.AddNodes(existing_nodes.begin(), existing_nodes.end());
        state.ResumeTiming();

        // benchmarking to add nodes for already existing model part with nodes
        r_model_part.AddNodes(nodes_to_be_added.begin(), nodes_to_be_added.end());
    }
}

static void BM_ModelPartAddNodesToSubSubFromRange1(benchmark::State& state) {
    for (auto _ : state) {
        const IndexType input_size = state.range(0);
        const bool fill_existing = state.range(1) == 1;

        state.PauseTiming();

        auto model = Model();
        auto& r_sub_sub_model_part = model.CreateModelPart("test").CreateSubModelPart("sub_test").CreateSubModelPart("sub_sub_test");

        std::vector<ModelPart::NodeType::Pointer> nodes_to_be_added, existing_nodes;
        if (!fill_existing) {
            nodes_to_be_added.reserve(input_size);
        } else {
            nodes_to_be_added.reserve(input_size / 2 + 1);
            existing_nodes.reserve(input_size / 2 + 1);
        }

        // doing it in reverse to make sure a proper sort is done always
        for (IndexType i = 0; i < input_size; ++i) {
            auto p_node = Kratos::make_intrusive<NodeType>(i + 1, 0.0, 0.0, 0.0);
            if (!fill_existing || i % 2 == 0) {
                nodes_to_be_added.insert(nodes_to_be_added.end(), p_node);
            } else {
                existing_nodes.insert(existing_nodes.end(), p_node);
            }
        }
        r_sub_sub_model_part.AddNodes(existing_nodes.begin(), existing_nodes.end());
        state.ResumeTiming();

        r_sub_sub_model_part.AddNodes(nodes_to_be_added.begin(), nodes_to_be_added.end());
    }
}

static void BM_ModelPartAddNodesToSubSubFromRange2(benchmark::State& state) {
    for (auto _ : state) {
        const IndexType input_size = state.range(0);
        const bool fill_existing = state.range(1) == 1;

        state.PauseTiming();

        auto model = Model();
        auto& r_sub_sub_model_part = model.CreateModelPart("test").CreateSubModelPart("sub_test").CreateSubModelPart("sub_sub_test");

        std::vector<ModelPart::NodeType::Pointer> nodes_to_be_added, existing_nodes;
        if (!fill_existing) {
            nodes_to_be_added.reserve(input_size);
        } else {
            nodes_to_be_added.reserve(input_size / 2 + 1);
            existing_nodes.reserve(input_size / 2 + 1);
        }

        // doing it in reverse to make sure a proper sort is done always
        for (IndexType i = 0; i < input_size; ++i) {
            auto p_node = Kratos::make_intrusive<NodeType>(i + 1, 0.0, 0.0, 0.0);
            if (!fill_existing || i % 2 == 0) {
                nodes_to_be_added.insert(nodes_to_be_added.end(), p_node);
            } else {
                existing_nodes.insert(existing_nodes.end(), p_node);
            }
        }
        r_sub_sub_model_part.GetRootModelPart().AddNodes(nodes_to_be_added.begin(), nodes_to_be_added.end());
        r_sub_sub_model_part.GetRootModelPart().AddNodes(existing_nodes.begin(), existing_nodes.end());
        state.ResumeTiming();

        r_sub_sub_model_part.AddNodes(nodes_to_be_added.begin(), nodes_to_be_added.end());
    }
}

static void BM_ModelPartAddNodesToSubSubFromRange3(benchmark::State& state) {
    for (auto _ : state) {
        const IndexType input_size = state.range(0);
        const bool fill_existing = state.range(1) == 1;

        state.PauseTiming();

        auto model = Model();
        auto& r_sub_sub_model_part = model.CreateModelPart("test").CreateSubModelPart("sub_test").CreateSubModelPart("sub_sub_test");

        std::vector<ModelPart::NodeType::Pointer> nodes_to_be_added, existing_nodes;
        if (!fill_existing) {
            nodes_to_be_added.reserve(input_size);
        } else {
            nodes_to_be_added.reserve(input_size / 2 + 1);
            existing_nodes.reserve(input_size / 2 + 1);
        }

        // doing it in reverse to make sure a proper sort is done always
        for (IndexType i = 0; i < input_size; ++i) {
            auto p_node = Kratos::make_intrusive<NodeType>(i + 1, 0.0, 0.0, 0.0);
            if (!fill_existing || i % 2 == 0) {
                nodes_to_be_added.insert(nodes_to_be_added.end(), p_node);
            } else {
                existing_nodes.insert(existing_nodes.end(), p_node);
            }
        }
        r_sub_sub_model_part.GetRootModelPart().AddNodes(nodes_to_be_added.begin(), nodes_to_be_added.end());
        r_sub_sub_model_part.GetRootModelPart().AddNodes(existing_nodes.begin(), existing_nodes.end());
        state.ResumeTiming();

        r_sub_sub_model_part.AddNodes(r_sub_sub_model_part.GetRootModelPart().NodesBegin(), r_sub_sub_model_part.GetRootModelPart().NodesEnd());
    }
}

static void BM_ModelPartAddNodesToSubSubFromId1(benchmark::State& state) {
    for (auto _ : state) {
        const IndexType input_size = state.range(0);
        const bool fill_existing = state.range(1) == 1;

        state.PauseTiming();

        auto model = Model();
        auto& r_root_model_part = model.CreateModelPart("test");
        auto& r_sub_sub_model_part = r_root_model_part.CreateSubModelPart("sub_test").CreateSubModelPart("sub_sub_test");

        std::vector<ModelPart::NodeType::Pointer> nodes;
        std::vector<IndexType> node_ids_to_add, existing_node_ids;
        if (!fill_existing) {
            node_ids_to_add.reserve(input_size);
        } else {
            node_ids_to_add.reserve(input_size / 2 + 1);
            existing_node_ids.reserve(input_size / 2 + 1);
        }

        // doing it in reverse to make sure a proper sort is done always
        for (IndexType i = input_size; i > 0; --i) {
            nodes.push_back(Kratos::make_intrusive<NodeType>(i, 0.0, 0.0, 0.0));
            if (!fill_existing || i % 2 == 0) {
                node_ids_to_add.push_back(i);
            } else {
                existing_node_ids.push_back(i);
            }
        }
        r_root_model_part.AddNodes(nodes.begin(), nodes.end());
        r_sub_sub_model_part.AddNodes(existing_node_ids);
        state.ResumeTiming();

        r_sub_sub_model_part.AddNodes(node_ids_to_add);
    }
}


// Register the function as a benchmark
BENCHMARK(BM_ModelPartCreateNewNode) -> RangeMultiplier(10) -> Range(1e2, 1e+6);
BENCHMARK(BM_ModelPartAddNodesToRootFromRange1) -> ArgsProduct({{100, 1000, 10000, 100000, 1000000}, {0, 1}});
BENCHMARK(BM_ModelPartAddNodesToRootFromRange2) -> ArgsProduct({{100, 1000, 10000, 100000, 1000000}, {0, 1}});
BENCHMARK(BM_ModelPartAddNodesToSubSubFromRange1) -> ArgsProduct({{100, 1000, 10000, 100000, 1000000}, {0, 1}});
BENCHMARK(BM_ModelPartAddNodesToSubSubFromRange2) -> ArgsProduct({{100, 1000, 10000, 100000, 1000000}, {0, 1}});
BENCHMARK(BM_ModelPartAddNodesToSubSubFromRange3) -> ArgsProduct({{100, 1000, 10000, 100000, 1000000}, {0, 1}});
BENCHMARK(BM_ModelPartAddNodesToSubSubFromId1) -> ArgsProduct({{100, 1000, 10000, 100000, 1000000}, {0, 1}});

}  // namespace Kratos

BENCHMARK_MAIN();