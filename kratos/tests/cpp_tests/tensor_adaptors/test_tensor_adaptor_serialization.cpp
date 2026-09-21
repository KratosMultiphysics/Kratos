//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Claude Sonnet 5
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

KRATOS_TEST_CASE_IN_SUITE(UninitializedTensorAdaptorThrows, KratosCoreFastSuite)
{
    TensorAdaptor<double> tensor_adaptor;
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(tensor_adaptor.Shape(), "Uninitialized TensorAdaptor");
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(tensor_adaptor.Size(), "Uninitialized TensorAdaptor");
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(tensor_adaptor.ViewData(), "Uninitialized TensorAdaptor");
}

} // namespace Kratos::Testing
