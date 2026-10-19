//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ihar Antonau
//

// System includes
#include <cstdio>
#include <fstream>

// External includes

// Project includes
#include "testing/testing.h"
#include "containers/model.h"
#include "includes/stream_serializer.h"
#include "includes/file_serializer.h"
#include "tensor_adaptors/variable_tensor_adaptor.h"
#include "tensor_adaptors/historical_variable_tensor_adaptor.h"
#include "tensor_adaptors/gauss_point_variable_tensor_adaptor.h"
#include "tensor_adaptors/equation_ids_tensor_adaptor.h"
#include "tensor_adaptors/fixity_tensor_adaptor.h"
#include "tensor_adaptors/flags_tensor_adaptor.h"
#include "tensor_adaptors/node_position_tensor_adaptor.h"
#include "tensor_adaptors/geometry_metrics_tensor_adaptor.h"
#include "tensor_adaptors/geometries_tensor_adaptor.h"
#include "tensor_adaptors/connectivity_ids_tensor_adaptor.h"
#include "tensor_adaptors/combined_tensor_adaptor.h"

namespace Kratos::Testing {

namespace {

ModelPart& CreateTestModelPart(Model& rModel, const std::string& rName)
{
    ModelPart& r_model_part = rModel.CreateModelPart(rName);
    r_model_part.AddNodalSolutionStepVariable(TEMPERATURE);
    r_model_part.SetBufferSize(2);
    r_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    r_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
    r_model_part.CreateNewNode(3, 0.0, 1.0, 0.0);
    auto p_properties = r_model_part.CreateNewProperties(1);
    r_model_part.CreateNewElement("Element2D3N", 1, {1, 2, 3}, p_properties);
    return r_model_part;
}

// Saves rOriginal and loads it into rTarget, either with a StreamSerializer or in DataOnly mode.
void SaveAndLoadNDData(
    const NDData<double>& rOriginal,
    NDData<double>& rTarget,
    const bool DataOnly)
{
    if (!DataOnly) {
        StreamSerializer serializer;
        serializer.save("Data", rOriginal);
        serializer.load("Data", rTarget);
    } else {
        const std::string file_name = "test_nd_data_data_only_load";
        struct RemoveFileOnExit { std::string mName; ~RemoveFileOnExit() { std::remove(mName.c_str()); } } remove_file{file_name + ".rest"};
        std::fstream(file_name + ".rest", std::ios::out).close();
        FileSerializer serializer(file_name, Serializer::SERIALIZER_NO_TRACE, /*DataOnly=*/true);
        serializer.save("Data", rOriginal);
        serializer.SetLoadState();
        serializer.load("Data", rTarget);
    }
}

// Root: [level 2, historical TEMPERATURE (step 1)], level 2: [level 3, node position, element PRESSURE],
// level 3: [nodal PRESSURE, historical TEMPERATURE (step 0)]. Every level collects and stores recursively.
CombinedTensorAdaptor<double>::Pointer CreateNestedCombinedTensorAdaptor(ModelPart& rModelPart)
{
    using combined_type = CombinedTensorAdaptor<double>;
    using vector_type = combined_type::TensorAdaptorVectorType;

    auto p_level_3 = Kratos::make_shared<combined_type>(vector_type{
        Kratos::make_shared<VariableTensorAdaptor>(rModelPart.pNodes(), &PRESSURE),
        Kratos::make_shared<HistoricalVariableTensorAdaptor>(rModelPart.pNodes(), &TEMPERATURE, 0)}, true, true, false);
    auto p_level_2 = Kratos::make_shared<combined_type>(vector_type{
        p_level_3,
        Kratos::make_shared<NodePositionTensorAdaptor>(rModelPart.pNodes(), Globals::Configuration::Current),
        Kratos::make_shared<VariableTensorAdaptor>(rModelPart.pElements(), &PRESSURE)}, true, true, false);
    return Kratos::make_shared<combined_type>(vector_type{
        p_level_2,
        Kratos::make_shared<HistoricalVariableTensorAdaptor>(rModelPart.pNodes(), &TEMPERATURE, 1)}, true, true, false);
}

void CheckNestedCombinedTensorAdaptorStructure(const CombinedTensorAdaptor<double>& rRoot)
{
    const auto level_1 = rRoot.GetTensorAdaptors();
    ASSERT_EQ(level_1.size(), 2u);
    const auto p_level_2 = std::dynamic_pointer_cast<CombinedTensorAdaptor<double>>(level_1[0]);
    ASSERT_NE(p_level_2, nullptr);
    KRATOS_EXPECT_NE(std::dynamic_pointer_cast<HistoricalVariableTensorAdaptor>(level_1[1]), nullptr);

    const auto level_2 = p_level_2->GetTensorAdaptors();
    ASSERT_EQ(level_2.size(), 3u);
    const auto p_level_3 = std::dynamic_pointer_cast<CombinedTensorAdaptor<double>>(level_2[0]);
    ASSERT_NE(p_level_3, nullptr);
    KRATOS_EXPECT_NE(std::dynamic_pointer_cast<NodePositionTensorAdaptor>(level_2[1]), nullptr);
    KRATOS_EXPECT_NE(std::dynamic_pointer_cast<VariableTensorAdaptor>(level_2[2]), nullptr);

    const auto level_3 = p_level_3->GetTensorAdaptors();
    ASSERT_EQ(level_3.size(), 2u);
    KRATOS_EXPECT_NE(std::dynamic_pointer_cast<VariableTensorAdaptor>(level_3[0]), nullptr);
    KRATOS_EXPECT_NE(std::dynamic_pointer_cast<HistoricalVariableTensorAdaptor>(level_3[1]), nullptr);
}

} // namespace

KRATOS_TEST_CASE_IN_SUITE(VariableTensorAdaptorSerialization, KratosCoreFastSuite)
{
    Model save_model;
    auto& r_save_model_part = CreateTestModelPart(save_model, "Test");

    VariableTensorAdaptor original(r_save_model_part.pNodes(), &PRESSURE);
    original.ViewData()[0] = 1.0;
    original.ViewData()[1] = 2.0;
    original.ViewData()[2] = 3.0;

    StreamSerializer serializer;
    serializer.Set(Serializer::SHALLOW_GLOBAL_POINTERS_SERIALIZATION);
    serializer.save("Model", save_model);
    serializer.save("TA", original);

    Model load_model;
    VariableTensorAdaptor loaded;
    serializer.load("Model", load_model);
    serializer.load("TA", loaded);

    KRATOS_EXPECT_EQ(loaded.Size(), original.Size());
    for (IndexType i = 0; i < original.Size(); ++i) {
        KRATOS_EXPECT_DOUBLE_EQ(loaded.ViewData()[i], original.ViewData()[i]);
    }
    KRATOS_EXPECT_EQ(loaded.Info(), original.Info());

    // the Model is saved first, so the container must resolve to the loaded model's own Nodes
    // container, not to a detached copy.
    ModelPart& r_load_model_part = load_model.GetModelPart("Test");
    KRATOS_EXPECT_EQ(std::get<ModelPart::NodesContainerType::Pointer>(loaded.GetContainer()).get(), r_load_model_part.pNodes().get());
}

KRATOS_TEST_CASE_IN_SUITE(HistoricalVariableTensorAdaptorSerialization, KratosCoreFastSuite)
{
    Model save_model;
    auto& r_save_model_part = CreateTestModelPart(save_model, "Test");
    for (auto& r_node : r_save_model_part.Nodes()) {
        r_node.FastGetSolutionStepValue(TEMPERATURE) = r_node.Id() * 10.0;
    }

    HistoricalVariableTensorAdaptor original(r_save_model_part.pNodes(), &TEMPERATURE, 0);
    original.CollectData();

    StreamSerializer serializer;
    serializer.save("TA", original);

    HistoricalVariableTensorAdaptor loaded;
    serializer.load("TA", loaded);

    KRATOS_EXPECT_EQ(loaded.Size(), original.Size());
    for (IndexType i = 0; i < original.Size(); ++i) {
        KRATOS_EXPECT_DOUBLE_EQ(loaded.ViewData()[i], original.ViewData()[i]);
    }
    // Info() prints the variable name and step index, so equal Info() means both were restored.
    KRATOS_EXPECT_EQ(loaded.Info(), original.Info());
}

KRATOS_TEST_CASE_IN_SUITE(FlagsTensorAdaptorSerialization, KratosCoreFastSuite)
{
    Model save_model;
    auto& r_save_model_part = CreateTestModelPart(save_model, "Test");

    FlagsTensorAdaptor original(r_save_model_part.pNodes(), ACTIVE);
    original.ViewData()[0] = 1;
    original.ViewData()[1] = -1;
    original.ViewData()[2] = 0;

    StreamSerializer serializer;
    serializer.save("TA", original);

    FlagsTensorAdaptor loaded;
    serializer.load("TA", loaded);

    KRATOS_EXPECT_EQ(loaded.Size(), original.Size());
    for (IndexType i = 0; i < original.Size(); ++i) {
        KRATOS_EXPECT_EQ(loaded.ViewData()[i], original.ViewData()[i]);
    }
    KRATOS_EXPECT_EQ(loaded.Info(), original.Info());
}

KRATOS_TEST_CASE_IN_SUITE(NodePositionTensorAdaptorSerialization, KratosCoreFastSuite)
{
    Model save_model;
    auto& r_save_model_part = CreateTestModelPart(save_model, "Test");

    NodePositionTensorAdaptor original(r_save_model_part.pNodes(), Globals::Configuration::Initial);
    original.CollectData();

    StreamSerializer serializer;
    serializer.save("TA", original);

    NodePositionTensorAdaptor loaded;
    serializer.load("TA", loaded);

    KRATOS_EXPECT_EQ(loaded.Size(), original.Size());
    for (IndexType i = 0; i < original.Size(); ++i) {
        KRATOS_EXPECT_DOUBLE_EQ(loaded.ViewData()[i], original.ViewData()[i]);
    }
    KRATOS_EXPECT_EQ(loaded.Info(), original.Info());
}

KRATOS_TEST_CASE_IN_SUITE(FixityTensorAdaptorSerialization, KratosCoreFastSuite)
{
    Model save_model;
    auto& r_save_model_part = CreateTestModelPart(save_model, "Test");
    for (auto& r_node : r_save_model_part.Nodes()) {
        r_node.AddDof(TEMPERATURE);
    }
    r_save_model_part.GetNode(1).Fix(TEMPERATURE);

    FixityTensorAdaptor original(r_save_model_part.pNodes(), std::vector<const Variable<double>*>{&TEMPERATURE});
    original.CollectData();

    StreamSerializer serializer;
    serializer.save("TA", original);

    FixityTensorAdaptor loaded;
    serializer.load("TA", loaded);

    KRATOS_EXPECT_EQ(loaded.Size(), original.Size());
    for (IndexType i = 0; i < original.Size(); ++i) {
        KRATOS_EXPECT_EQ(loaded.ViewData()[i], original.ViewData()[i]);
    }
}

KRATOS_TEST_CASE_IN_SUITE(EquationIdsTensorAdaptorSerialization, KratosCoreFastSuite)
{
    Model save_model;
    auto& r_save_model_part = CreateTestModelPart(save_model, "Test");

    // build from a shaped base adaptor: the container constructor needs real DOFs, which are not
    // what this test checks (data and mpProcessInfo round trip).
    auto p_nd_data = Kratos::make_shared<NDData<int>>(DenseVector<unsigned int>(1, r_save_model_part.NumberOfElements()));
    TensorAdaptor<int> base(r_save_model_part.pElements(), p_nd_data, false);
    EquationIdsTensorAdaptor original(base, r_save_model_part.pGetProcessInfo(), false);
    original.ViewData()[0] = 42;

    StreamSerializer serializer;
    serializer.save("TA", original);

    EquationIdsTensorAdaptor loaded;
    serializer.load("TA", loaded);

    KRATOS_EXPECT_EQ(loaded.Size(), original.Size());
    KRATOS_EXPECT_EQ(loaded.ViewData()[0], original.ViewData()[0]);
}

KRATOS_TEST_CASE_IN_SUITE(GaussPointVariableTensorAdaptorSerialization, KratosCoreFastSuite)
{
    Model save_model;
    auto& r_save_model_part = CreateTestModelPart(save_model, "Test");

    // build from a shaped base adaptor (see EquationIdsTensorAdaptorSerialization): only the data,
    // mpVariable and mpProcessInfo round trip is tested.
    DenseVector<unsigned int> shape(2);
    shape[0] = r_save_model_part.NumberOfElements();
    shape[1] = 1;
    auto p_nd_data = Kratos::make_shared<NDData<double>>(shape);
    TensorAdaptor<double> base(r_save_model_part.pElements(), p_nd_data, false);
    GaussPointVariableTensorAdaptor original(base, &PRESSURE, r_save_model_part.pGetProcessInfo(), false);
    original.ViewData()[0] = 3.5;

    StreamSerializer serializer;
    serializer.save("TA", original);

    GaussPointVariableTensorAdaptor loaded;
    serializer.load("TA", loaded);

    KRATOS_EXPECT_EQ(loaded.Size(), original.Size());
    KRATOS_EXPECT_DOUBLE_EQ(loaded.ViewData()[0], original.ViewData()[0]);
    KRATOS_EXPECT_EQ(loaded.Info(), original.Info());
}

KRATOS_TEST_CASE_IN_SUITE(GeometryMetricsTensorAdaptorSerialization, KratosCoreFastSuite)
{
    Model save_model;
    auto& r_save_model_part = CreateTestModelPart(save_model, "Test");

    GeometryMetricsTensorAdaptor original(r_save_model_part.pElements(), GeometryMetricsTensorAdaptor::Metric::DomainSize);
    original.CollectData();
    KRATOS_EXPECT_GT(original.ViewData()[0], 0.0); // sanity: real domain size computed, not left at zero

    StreamSerializer serializer;
    serializer.save("TA", original);

    GeometryMetricsTensorAdaptor loaded;
    serializer.load("TA", loaded);

    KRATOS_EXPECT_EQ(loaded.Size(), original.Size());
    KRATOS_EXPECT_DOUBLE_EQ(loaded.ViewData()[0], original.ViewData()[0]);
    KRATOS_EXPECT_EQ(loaded.Info(), original.Info());
}

KRATOS_TEST_CASE_IN_SUITE(GeometriesTensorAdaptorSerialization, KratosCoreFastSuite)
{
    Model save_model;
    auto& r_save_model_part = CreateTestModelPart(save_model, "Test");

    GeometriesTensorAdaptor original(r_save_model_part.pElements(), GeometriesTensorAdaptor::DatumType::IntegrationWeights);
    original.CollectData();
    KRATOS_EXPECT_GT(original.ViewData()[0], 0.0);

    StreamSerializer serializer;
    serializer.save("TA", original);

    GeometriesTensorAdaptor loaded;
    serializer.load("TA", loaded);

    KRATOS_EXPECT_EQ(loaded.Size(), original.Size());
    KRATOS_EXPECT_DOUBLE_EQ(loaded.ViewData()[0], original.ViewData()[0]);
    KRATOS_EXPECT_EQ(loaded.Info(), original.Info());
}

KRATOS_TEST_CASE_IN_SUITE(ConnectivityIdsTensorAdaptorSerialization, KratosCoreFastSuite)
{
    Model save_model;
    auto& r_save_model_part = CreateTestModelPart(save_model, "Test");

    ConnectivityIdsTensorAdaptor original(r_save_model_part.pElements());
    original.CollectData();
    KRATOS_EXPECT_GT(original.ViewData()[0], 0); // sanity: a real node id, not left at zero

    StreamSerializer serializer;
    serializer.save("TA", original);

    ConnectivityIdsTensorAdaptor loaded;
    serializer.load("TA", loaded);

    KRATOS_EXPECT_EQ(loaded.Size(), original.Size());
    for (IndexType i = 0; i < original.Size(); ++i) {
        KRATOS_EXPECT_EQ(loaded.ViewData()[i], original.ViewData()[i]);
    }
}

// children are saved through polymorphic TensorAdaptor<double> pointers. Before the subtypes were
// registered, this threw "no object registered with type id ..."; now each child must come back
// as its own subtype with its data and members.
KRATOS_TEST_CASE_IN_SUITE(CombinedTensorAdaptorWithDerivedChildrenSerialization, KratosCoreFastSuite)
{
    Model save_model;
    auto& r_save_model_part = CreateTestModelPart(save_model, "Test");
    for (auto& r_node : r_save_model_part.Nodes()) {
        r_node.FastGetSolutionStepValue(TEMPERATURE) = r_node.Id() * 10.0;
    }

    auto p_variable_ta = Kratos::make_shared<VariableTensorAdaptor>(r_save_model_part.pNodes(), &PRESSURE);
    p_variable_ta->ViewData()[0] = 1.0;
    p_variable_ta->ViewData()[1] = 2.0;
    p_variable_ta->ViewData()[2] = 3.0;

    auto p_historical_ta = Kratos::make_shared<HistoricalVariableTensorAdaptor>(r_save_model_part.pNodes(), &TEMPERATURE, 0);
    p_historical_ta->CollectData();

    CombinedTensorAdaptor<double> original(
        CombinedTensorAdaptor<double>::TensorAdaptorVectorType{p_variable_ta, p_historical_ta},
        /*PerformCollectDataRecursively=*/false, /*PerformStoreDataRecursively=*/false, /*Copy=*/false);
    original.CollectData();

    StreamSerializer serializer;
    serializer.Set(Serializer::SHALLOW_GLOBAL_POINTERS_SERIALIZATION);
    serializer.save("Model", save_model);
    serializer.save("TA", original);

    Model load_model;
    CombinedTensorAdaptor<double> loaded;
    serializer.load("Model", load_model);
    serializer.load("TA", loaded);

    KRATOS_EXPECT_EQ(loaded.Size(), original.Size());
    for (IndexType i = 0; i < original.Size(); ++i) {
        KRATOS_EXPECT_DOUBLE_EQ(loaded.ViewData()[i], original.ViewData()[i]);
    }

    const auto loaded_children = loaded.GetTensorAdaptors();
    KRATOS_EXPECT_EQ(loaded_children.size(), 2);
    KRATOS_EXPECT_NE(std::dynamic_pointer_cast<VariableTensorAdaptor>(loaded_children[0]), nullptr);
    KRATOS_EXPECT_NE(std::dynamic_pointer_cast<HistoricalVariableTensorAdaptor>(loaded_children[1]), nullptr);

    ModelPart& r_load_model_part = load_model.GetModelPart("Test");
    KRATOS_EXPECT_EQ(std::get<ModelPart::NodesContainerType::Pointer>(loaded_children[0]->GetContainer()).get(), r_load_model_part.pNodes().get());
}

KRATOS_TEST_CASE_IN_SUITE(NestedCombinedTensorAdaptorSerialization, KratosCoreFastSuite)
{
    Model save_model;
    auto& r_save_model_part = CreateTestModelPart(save_model, "Test");
    for (auto& r_node : r_save_model_part.Nodes()) {
        r_node.SetValue(PRESSURE, r_node.Id() * 1.0);
        r_node.FastGetSolutionStepValue(TEMPERATURE, 0) = r_node.Id() * 10.0;
        r_node.FastGetSolutionStepValue(TEMPERATURE, 1) = r_node.Id() * 100.0;
    }
    r_save_model_part.GetElement(1).SetValue(PRESSURE, 5.0);

    const auto p_original = CreateNestedCombinedTensorAdaptor(r_save_model_part);
    p_original->CollectData();
    KRATOS_EXPECT_EQ(p_original->Size(), 3u + 3u + 9u + 1u + 3u);

    StreamSerializer serializer;
    serializer.Set(Serializer::SHALLOW_GLOBAL_POINTERS_SERIALIZATION);
    serializer.save("Model", save_model);
    serializer.save("TA", *p_original);

    Model load_model;
    CombinedTensorAdaptor<double> loaded;
    serializer.load("Model", load_model);
    serializer.load("TA", loaded);

    CheckNestedCombinedTensorAdaptorStructure(loaded);
    // Info() lists the recursive flags of each level and the variables, step indices and containers of the leaves.
    KRATOS_EXPECT_EQ(loaded.Info(), p_original->Info());
    KRATOS_EXPECT_EQ(loaded.Size(), p_original->Size());
    for (IndexType i = 0; i < p_original->Size(); ++i) {
        KRATOS_EXPECT_DOUBLE_EQ(loaded.ViewData()[i], p_original->ViewData()[i]);
    }

    // the recursive StoreData() must reach the loaded model through all levels: an identical structure
    // built on the loaded model collects the stored values, and the original model is untouched.
    for (IndexType i = 0; i < loaded.Size(); ++i) {
        loaded.ViewData()[i] = -1.0 - i;
    }
    loaded.StoreData();

    const auto p_check = CreateNestedCombinedTensorAdaptor(load_model.GetModelPart("Test"));
    p_check->CollectData();
    for (IndexType i = 0; i < loaded.Size(); ++i) {
        KRATOS_EXPECT_DOUBLE_EQ(p_check->ViewData()[i], -1.0 - i);
    }

    const auto p_original_check = CreateNestedCombinedTensorAdaptor(r_save_model_part);
    p_original_check->CollectData();
    for (IndexType i = 0; i < p_original->Size(); ++i) {
        KRATOS_EXPECT_DOUBLE_EQ(p_original_check->ViewData()[i], p_original->ViewData()[i]);
    }
}

// DataOnly mode loads into a prebuilt identical structure, which must keep its own children and containers.
KRATOS_TEST_CASE_IN_SUITE(NestedCombinedTensorAdaptorDataOnlySerialization, KratosCoreFastSuite)
{
    Model model;
    auto& r_save_model_part = CreateTestModelPart(model, "Save");
    auto& r_load_model_part = CreateTestModelPart(model, "Load");
    for (auto& r_node : r_save_model_part.Nodes()) {
        r_node.SetValue(PRESSURE, r_node.Id() * 1.0);
        r_node.FastGetSolutionStepValue(TEMPERATURE, 0) = r_node.Id() * 10.0;
        r_node.FastGetSolutionStepValue(TEMPERATURE, 1) = r_node.Id() * 100.0;
    }
    r_save_model_part.GetElement(1).SetValue(PRESSURE, 5.0);
    // DataOnly also loads the saved containers' entities, which needs the same non-historical variables.
    for (auto& r_node : r_load_model_part.Nodes()) {
        r_node.SetValue(PRESSURE, 0.0);
    }
    r_load_model_part.GetElement(1).SetValue(PRESSURE, 0.0);

    const auto p_original = CreateNestedCombinedTensorAdaptor(r_save_model_part);
    p_original->CollectData();

    const auto p_loaded = CreateNestedCombinedTensorAdaptor(r_load_model_part);
    const auto loaded_children = p_loaded->GetTensorAdaptors();

    const std::string file_name = "test_nested_combined_tensor_adaptor_data_only_load";
    struct RemoveFileOnExit { std::string mName; ~RemoveFileOnExit() { std::remove(mName.c_str()); } } remove_file{file_name + ".rest"};
    std::fstream(file_name + ".rest", std::ios::out).close();
    {
        FileSerializer serializer(file_name, Serializer::SERIALIZER_NO_TRACE, /*DataOnly=*/true);
        serializer.save("TA", *p_original);
        serializer.SetLoadState();
        serializer.load("TA", *p_loaded);
    }

    CheckNestedCombinedTensorAdaptorStructure(*p_loaded);
    KRATOS_EXPECT_EQ(p_loaded->GetTensorAdaptors()[0].get(), loaded_children[0].get());
    KRATOS_EXPECT_EQ(p_loaded->GetTensorAdaptors()[1].get(), loaded_children[1].get());
    KRATOS_EXPECT_EQ(std::get<ModelPart::NodesContainerType::Pointer>(p_loaded->GetTensorAdaptors()[1]->GetContainer()).get(), r_load_model_part.pNodes().get());
    KRATOS_EXPECT_EQ(p_loaded->Size(), p_original->Size());
    for (IndexType i = 0; i < p_original->Size(); ++i) {
        KRATOS_EXPECT_DOUBLE_EQ(p_loaded->ViewData()[i], p_original->ViewData()[i]);
    }
}

// DataOnly mode reads a shared_ptr's data only into a non-null (prebuilt) target. load() used to use
// a fresh null pointer for the container, which read nothing, desynced the stream and dropped the
// preinitialized container. The empty constraint container keeps entity-level DataOnly out of it.
KRATOS_TEST_CASE_IN_SUITE(TensorAdaptorDataOnlyLoadReusesExistingContainer, KratosCoreFastSuite)
{
    Model model;
    ModelPart& r_model_part = model.CreateModelPart("Test");
    auto p_constraints = r_model_part.GetMesh().pMasterSlaveConstraints();

    VariableTensorAdaptor original(p_constraints, &PRESSURE);

    // `loaded` is prebuilt on the same container, as a restart reload does before a DataOnly load.
    VariableTensorAdaptor loaded(p_constraints, &PRESSURE);

    const std::string file_name = "test_tensor_adaptor_data_only_load";
    struct RemoveFileOnExit { std::string mName; ~RemoveFileOnExit() { std::remove(mName.c_str()); } } remove_file{file_name + ".rest"};
    // create the file first: without it FileSerializer falls back to an out-only stream and the
    // load() below reads garbage.
    std::fstream(file_name + ".rest", std::ios::out).close();
    {
        FileSerializer serializer(file_name, Serializer::SERIALIZER_NO_TRACE, /*DataOnly=*/true);
        serializer.save("TA", original);
        // the fstream shares one read/write position (unlike a stringstream), so rewind before reading.
        serializer.SetLoadState();
        serializer.load("TA", loaded);
    }

    KRATOS_EXPECT_EQ(std::get<ModelPart::MasterSlaveConstraintContainerType::Pointer>(loaded.GetContainer()).get(), p_constraints.get());
}

// load(shared_ptr&) loads a non-null target through its existing type. CombinedTensorAdaptor::load()
// used to keep existing children, so the stale VariableTensorAdaptor child below would run the wrong
// load() on the saved HistoricalVariableTensorAdaptor.
KRATOS_TEST_CASE_IN_SUITE(CombinedTensorAdaptorLoadDiscardsStaleDerivedTypeChild, KratosCoreFastSuite)
{
    Model save_model;
    auto& r_save_model_part = CreateTestModelPart(save_model, "Test");
    for (auto& r_node : r_save_model_part.Nodes()) {
        r_node.FastGetSolutionStepValue(TEMPERATURE) = r_node.Id() * 10.0;
    }

    auto p_historical_ta = Kratos::make_shared<HistoricalVariableTensorAdaptor>(r_save_model_part.pNodes(), &TEMPERATURE, 0);
    p_historical_ta->CollectData();

    CombinedTensorAdaptor<double> original(
        CombinedTensorAdaptor<double>::TensorAdaptorVectorType{p_historical_ta},
        /*PerformCollectDataRecursively=*/false, /*PerformStoreDataRecursively=*/false, /*Copy=*/false);
    original.CollectData();

    StreamSerializer serializer;
    serializer.Set(Serializer::SHALLOW_GLOBAL_POINTERS_SERIALIZATION);
    serializer.save("Model", save_model);
    serializer.save("TA", original);

    Model load_model;
    auto p_stale_child = Kratos::make_shared<VariableTensorAdaptor>(r_save_model_part.pNodes(), &PRESSURE);
    CombinedTensorAdaptor<double> loaded(
        CombinedTensorAdaptor<double>::TensorAdaptorVectorType{p_stale_child},
        /*PerformCollectDataRecursively=*/false, /*PerformStoreDataRecursively=*/false, /*Copy=*/false);

    serializer.load("Model", load_model);
    serializer.load("TA", loaded);

    KRATOS_EXPECT_EQ(loaded.Size(), original.Size());
    for (IndexType i = 0; i < original.Size(); ++i) {
        KRATOS_EXPECT_DOUBLE_EQ(loaded.ViewData()[i], original.ViewData()[i]);
    }

    const auto loaded_children = loaded.GetTensorAdaptors();
    KRATOS_EXPECT_EQ(loaded_children.size(), 1);
    KRATOS_EXPECT_NE(std::dynamic_pointer_cast<HistoricalVariableTensorAdaptor>(loaded_children[0]), nullptr);
}

KRATOS_TEST_CASE_IN_SUITE(TensorAdaptorSerializationAsciiTrace, KratosCoreFastSuite)
{
    Model save_model;
    auto& r_save_model_part = CreateTestModelPart(save_model, "Test");

    VariableTensorAdaptor original(r_save_model_part.pNodes(), &PRESSURE);
    original.ViewData()[0] = 1.5;
    original.ViewData()[1] = -2.25;
    original.ViewData()[2] = 3.0;

    StreamSerializer serializer(Serializer::SERIALIZER_TRACE_ALL);
    serializer.Set(Serializer::SHALLOW_GLOBAL_POINTERS_SERIALIZATION);
    serializer.save("Model", save_model);
    serializer.save("TA", original);

    Model load_model;
    VariableTensorAdaptor loaded;
    serializer.load("Model", load_model);
    serializer.load("TA", loaded);

    KRATOS_EXPECT_EQ(loaded.Size(), original.Size());
    for (IndexType i = 0; i < original.Size(); ++i) {
        KRATOS_EXPECT_DOUBLE_EQ(loaded.ViewData()[i], original.ViewData()[i]);
    }
    KRATOS_EXPECT_EQ(std::get<ModelPart::NodesContainerType::Pointer>(loaded.GetContainer()).get(), load_model.GetModelPart("Test").pNodes().get());
}

KRATOS_TEST_CASE_IN_SUITE(BoolTensorAdaptorSerialization, KratosCoreFastSuite)
{
    Model save_model;
    auto& r_save_model_part = CreateTestModelPart(save_model, "Test");

    auto p_nd_data = Kratos::make_shared<NDData<bool>>(DenseVector<unsigned int>(1, r_save_model_part.NumberOfNodes()));
    TensorAdaptor<bool> original(r_save_model_part.pNodes(), p_nd_data, false);
    original.ViewData()[0] = true;
    original.ViewData()[1] = false;
    original.ViewData()[2] = true;

    StreamSerializer serializer;
    serializer.save("TA", original);

    TensorAdaptor<bool> loaded;
    serializer.load("TA", loaded);

    KRATOS_EXPECT_EQ(loaded.Size(), original.Size());
    for (IndexType i = 0; i < original.Size(); ++i) {
        KRATOS_EXPECT_EQ(loaded.ViewData()[i], original.ViewData()[i]);
    }
}

KRATOS_TEST_CASE_IN_SUITE(NDDataUnsignedCharSerializationAsciiTrace, KratosCoreFastSuite)
{
    NDData<unsigned char> original(DenseVector<unsigned int>(1, 3));
    original.ViewData()[0] = 0;
    original.ViewData()[1] = 32; // ' ', skipped by operator>> if written as a character
    original.ViewData()[2] = 255;

    StreamSerializer serializer(Serializer::SERIALIZER_TRACE_ALL);
    serializer.save("Data", original);

    NDData<unsigned char> loaded;
    serializer.load("Data", loaded);

    KRATOS_EXPECT_EQ(loaded.Size(), original.Size());
    for (IndexType i = 0; i < original.Size(); ++i) {
        KRATOS_EXPECT_EQ(loaded.ViewData()[i], original.ViewData()[i]);
    }
}

// A freshly allocated buffer must be owned, even if the saved data was foreign (unmanaged) memory.
KRATOS_TEST_CASE_IN_SUITE(NDDataLoadIntoEmptyOwnsData, KratosCoreFastSuite)
{
    for (const bool data_only : {false, true}) {
        std::vector<double> foreign{1.0, 2.0, 3.0};
        NDData<double> original(foreign.data(), DenseVector<unsigned int>(1, 3), false);
        KRATOS_EXPECT_FALSE(original.pData()->IsManaged());

        NDData<double> loaded;
        SaveAndLoadNDData(original, loaded, data_only);

        KRATOS_EXPECT_EQ(loaded.Size(), 3u);
        KRATOS_EXPECT_TRUE(loaded.pData()->IsManaged());
        KRATOS_EXPECT_NE(loaded.ViewData().data(), foreign.data());
        for (IndexType i = 0; i < 3; ++i) {
            KRATOS_EXPECT_EQ(loaded.ViewData()[i], foreign[i]);
        }
    }
}

// Same-size load overwrites the existing buffer in place (numpy views stay valid), keeping its ownership.
KRATOS_TEST_CASE_IN_SUITE(NDDataLoadInPlaceKeepsForeignMemory, KratosCoreFastSuite)
{
    for (const bool data_only : {false, true}) {
        NDData<double> original(DenseVector<unsigned int>(1, 3));
        for (IndexType i = 0; i < 3; ++i) {
            original.ViewData()[i] = i + 1.0;
        }

        std::vector<double> buffer(3, 0.0);
        NDData<double> target(buffer.data(), DenseVector<unsigned int>(1, 3), false);
        const auto p_wrapper = target.pData();

        SaveAndLoadNDData(original, target, data_only);

        KRATOS_EXPECT_EQ(target.pData().get(), p_wrapper.get());
        KRATOS_EXPECT_EQ(target.ViewData().data(), buffer.data());
        KRATOS_EXPECT_FALSE(target.pData()->IsManaged());
        for (IndexType i = 0; i < 3; ++i) {
            KRATOS_EXPECT_EQ(buffer[i], i + 1.0);
        }
    }
}

// A load with a different size gets a new owning wrapper. The old wrapper (and anything viewing
// its buffer, e.g. a numpy array) is left untouched.
KRATOS_TEST_CASE_IN_SUITE(NDDataLoadSizeMismatchReallocates, KratosCoreFastSuite)
{
    for (const bool data_only : {false, true}) {
        NDData<double> original(DenseVector<unsigned int>(1, 5));
        for (IndexType i = 0; i < 5; ++i) {
            original.ViewData()[i] = i + 1.0;
        }

        std::vector<double> buffer{7.0, 8.0, 9.0};
        NDData<double> foreign_target(buffer.data(), DenseVector<unsigned int>(1, 3), false);

        NDData<double> managed_target(DenseVector<unsigned int>(1, 3), 4.0);
        const auto p_old_wrapper = managed_target.pData();
        const double* p_old_data = p_old_wrapper->Data();

        SaveAndLoadNDData(original, foreign_target, data_only);
        SaveAndLoadNDData(original, managed_target, data_only);

        for (const auto* p_target : {&foreign_target, &managed_target}) {
            KRATOS_EXPECT_EQ(p_target->Shape()[0], 5u);
            KRATOS_EXPECT_EQ(p_target->Size(), 5u);
            KRATOS_EXPECT_TRUE(p_target->pData()->IsManaged());
            for (IndexType i = 0; i < 5; ++i) {
                KRATOS_EXPECT_EQ(p_target->ViewData()[i], i + 1.0);
            }
        }

        KRATOS_EXPECT_NE(foreign_target.ViewData().data(), buffer.data());
        KRATOS_EXPECT_EQ(buffer[0], 7.0);
        KRATOS_EXPECT_EQ(buffer[1], 8.0);
        KRATOS_EXPECT_EQ(buffer[2], 9.0);

        KRATOS_EXPECT_NE(managed_target.pData().get(), p_old_wrapper.get());
        KRATOS_EXPECT_EQ(p_old_wrapper->Data(), p_old_data);
        KRATOS_EXPECT_EQ(p_old_wrapper->Size(), 3u);
        KRATOS_EXPECT_EQ(p_old_data[0], 4.0);
    }
}

// A default-constructed wrapper is only a load target: every accessor and save throw until it is loaded.
KRATOS_TEST_CASE_IN_SUITE(NDDataPointerWrapperDefaultConstructed, KratosCoreFastSuite)
{
    using PointerWrapperType = NDData<double>::PointerWrapper;
    const std::string message = "Uninitialized PointerWrapper (default-constructed, only valid as a Serializer load target).";

    PointerWrapperType wrapper;
    const PointerWrapperType& r_const_wrapper = wrapper;
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(wrapper.Data(), message);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(r_const_wrapper.Data(), message);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(wrapper.Size(), message);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(wrapper.IsManaged(), message);
    StreamSerializer save_serializer;
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(save_serializer.save("Wrapper", wrapper), message);

    // loading into the default wrapper initializes it, with a new managed buffer
    for (const bool data_only : {false, true}) {
        std::vector<double> foreign{1.0, 2.0, 3.0};
        const auto p_original = Kratos::make_intrusive<PointerWrapperType>(foreign.data(), foreign.size(), false);

        auto p_loaded = Kratos::make_intrusive<PointerWrapperType>();
        if (!data_only) {
            StreamSerializer serializer;
            serializer.save("Wrapper", p_original);
            serializer.load("Wrapper", p_loaded);
        } else {
            const std::string file_name = "test_nd_data_pointer_wrapper_data_only_load";
            struct RemoveFileOnExit { std::string mName; ~RemoveFileOnExit() { std::remove(mName.c_str()); } } remove_file{file_name + ".rest"};
            std::fstream(file_name + ".rest", std::ios::out).close();
            FileSerializer serializer(file_name, Serializer::SERIALIZER_NO_TRACE, /*DataOnly=*/true);
            serializer.save("Wrapper", p_original);
            serializer.SetLoadState();
            serializer.load("Wrapper", p_loaded);
        }

        KRATOS_EXPECT_EQ(p_loaded->Size(), 3u);
        KRATOS_EXPECT_TRUE(p_loaded->IsManaged());
        KRATOS_EXPECT_NE(p_loaded->Data(), foreign.data());
        for (IndexType i = 0; i < 3; ++i) {
            KRATOS_EXPECT_EQ(p_loaded->Data()[i], foreign[i]);
        }
    }

    // a default NDData holds a default wrapper, so its data, shape and save throw through the wrapper
    NDData<double> nd_data;
    const NDData<double>& r_const_nd_data = nd_data;
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(nd_data.ViewData(), message);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(r_const_nd_data.ViewData(), message);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(nd_data.Shape(), message);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(nd_data.Size(), message);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(nd_data.Info(), message);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(NDData<double>{nd_data}, message);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(save_serializer.save("Data", nd_data), message);
}

KRATOS_TEST_CASE_IN_SUITE(DofsTensorAdaptorSerialization, KratosCoreFastSuite)
{
    Model save_model;
    auto& r_save_model_part = CreateTestModelPart(save_model, "Test");
    auto p_dofs = Kratos::make_shared<ModelPart::DofsArrayType>();
    for (auto& r_node : r_save_model_part.Nodes()) {
        p_dofs->push_back(r_node.pAddDof(TEMPERATURE));
    }
    p_dofs->Sort();

    auto p_nd_data = Kratos::make_shared<NDData<double>>(DenseVector<unsigned int>(1, p_dofs->size()));
    TensorAdaptor<double> original(p_dofs, p_nd_data, false);

    StreamSerializer serializer;
    serializer.Set(Serializer::SHALLOW_GLOBAL_POINTERS_SERIALIZATION);
    serializer.save("Model", save_model);
    serializer.save("TA", original);

    Model load_model;
    TensorAdaptor<double> loaded;
    serializer.load("Model", load_model);
    serializer.load("TA", loaded);

    // Dofs hold a raw pointer to their node's data, so the loaded container must hold the loaded
    // nodes' own Dof objects (same address), not detached copies.
    auto& r_load_model_part = load_model.GetModelPart("Test");
    const auto& r_loaded_dofs = *std::get<ModelPart::DofsArrayType::Pointer>(loaded.GetContainer());
    KRATOS_EXPECT_EQ(r_loaded_dofs.size(), p_dofs->size());
    for (const auto& r_dof : r_loaded_dofs) {
        auto& r_node = r_load_model_part.GetNode(r_dof.Id());
        KRATOS_EXPECT_EQ(&r_dof, r_node.pGetDof(TEMPERATURE));
    }
}

// A default-constructed TA holds a default NDData, so every data access on it (also through Info(),
// Check() and CollectData()) throws from the PointerWrapper until it is loaded.
KRATOS_TEST_CASE_IN_SUITE(UninitializedTensorAdaptorThrows, KratosCoreFastSuite)
{
    const std::string message = "Uninitialized PointerWrapper";

    TensorAdaptor<double> tensor_adaptor;
    const TensorAdaptor<double>& r_const_tensor_adaptor = tensor_adaptor;
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(tensor_adaptor.Shape(), message);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(tensor_adaptor.DataShape(), message);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(tensor_adaptor.Size(), message);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(tensor_adaptor.ViewData(), message);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(r_const_tensor_adaptor.ViewData(), message);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(tensor_adaptor.Info(), message);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(tensor_adaptor.Check(), message);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(tensor_adaptor.GetContainer(), message);
    StreamSerializer save_serializer;
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(save_serializer.save("TA", tensor_adaptor), message);

    VariableTensorAdaptor variable_tensor_adaptor;
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(variable_tensor_adaptor.Info(), message);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(variable_tensor_adaptor.Check(), message);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(variable_tensor_adaptor.CollectData(), message);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(variable_tensor_adaptor.StoreData(), message);

    ConnectivityIdsTensorAdaptor connectivity_tensor_adaptor;
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(connectivity_tensor_adaptor.Info(), message);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(connectivity_tensor_adaptor.Check(), message);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(connectivity_tensor_adaptor.CollectData(), message);

    FixityTensorAdaptor fixity_tensor_adaptor;
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(fixity_tensor_adaptor.CollectData(), message);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(fixity_tensor_adaptor.StoreData(), message);
}

// The default storage is a valid load target in both modes. A container-less TA is used, because
// DataOnly mode can't create a missing container.
KRATOS_TEST_CASE_IN_SUITE(UninitializedTensorAdaptorLoad, KratosCoreFastSuite)
{
    Model model;
    auto& r_model_part = CreateTestModelPart(model, "Test");
    auto p_child = Kratos::make_shared<VariableTensorAdaptor>(r_model_part.pNodes(), &PRESSURE);
    const CombinedTensorAdaptor<double> combined(CombinedTensorAdaptor<double>::TensorAdaptorVectorType{p_child}, false, false);
    TensorAdaptor<double> original(combined, true);
    KRATOS_EXPECT_FALSE(original.HasContainer());
    for (IndexType i = 0; i < original.Size(); ++i) {
        original.ViewData()[i] = i + 1.0;
    }

    for (const bool data_only : {false, true}) {
        TensorAdaptor<double> loaded;
        if (!data_only) {
            StreamSerializer serializer;
            serializer.save("TA", original);
            serializer.load("TA", loaded);
        } else {
            const std::string file_name = "test_uninitialized_tensor_adaptor_data_only_load";
            struct RemoveFileOnExit { std::string mName; ~RemoveFileOnExit() { std::remove(mName.c_str()); } } remove_file{file_name + ".rest"};
            std::fstream(file_name + ".rest", std::ios::out).close();
            FileSerializer serializer(file_name, Serializer::SERIALIZER_NO_TRACE, /*DataOnly=*/true);
            serializer.save("TA", original);
            serializer.SetLoadState();
            serializer.load("TA", loaded);
        }

        KRATOS_EXPECT_FALSE(loaded.HasContainer());
        KRATOS_EXPECT_VECTOR_EQ(loaded.Shape(), original.Shape());
        for (IndexType i = 0; i < original.Size(); ++i) {
            KRATOS_EXPECT_DOUBLE_EQ(loaded.ViewData()[i], original.ViewData()[i]);
        }
    }
}

KRATOS_TEST_CASE_IN_SUITE(NDDAtaCopyConstructor, KratosCoreFastSuite)
{

    const std::string message = "Uninitialized PointerWrapper (default-constructed, only valid as a Serializer load target)";
    NDData<double> empty;
    NDData<double> valid(DenseVector<unsigned int>(1, 3), 10.0);
    valid = empty;

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(valid.Info(), message);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(valid.Shape(), message);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(valid.Size(), message);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(valid.ViewData(), message);
}

} // namespace Kratos::Testing
