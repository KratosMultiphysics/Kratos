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

// One StreamSerializer per test: saves the Model alongside the tensor adaptor(s) under test, so
// SHALLOW_GLOBAL_POINTERS_SERIALIZATION relinks the loaded adaptor's container to the same
// Node/Element objects the loaded Model holds -- not a disconnected copy.
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

    // pointer-identity check: the loaded adaptor's container is the *same* Nodes container the
    // loaded Model holds, not a disconnected reconstruction.
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
    // Info() encodes the variable name and step index -- confirms both round-tripped.
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

    // adopt an already-shaped base TensorAdaptor -- EquationIdsTensorAdaptor's own construction
    // path requires live DOFs, which is unrelated to what this test verifies (that mpProcessInfo
    // round-trips through the Serializer).
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

    // adopt path (see EquationIdsTensorAdaptorSerialization above for why): only mpVariable and
    // mpProcessInfo round-tripping is under test here.
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

// The regression case this whole change is for: CombinedTensorAdaptor::save (already core
// functionality before this change) serializes each child through a polymorphic
// TensorAdaptor<double> pointer. Before this change, only the base/combined prototypes were
// registered, so a combined field holding derived children (as any real control field does, e.g.
// VariableTensorAdaptor/HistoricalVariableTensorAdaptor leaves) threw "no object registered with
// type id ...". This must now round-trip cleanly, each child recovering its own concrete type,
// values, and extra state.
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

    // pointer-identity check on a child, mirroring VariableTensorAdaptorSerialization above.
    ModelPart& r_load_model_part = load_model.GetModelPart("Test");
    KRATOS_EXPECT_EQ(std::get<ModelPart::NodesContainerType::Pointer>(loaded_children[0]->GetContainer()).get(), r_load_model_part.pNodes().get());
}

// Regression for TensorAdaptor::load()'s container-reconstruction switch always declaring a fresh,
// null local pointer before loading into it. Serializer's DataOnly mode only reads a shared_ptr's data
// when the destination is already non-null -- the caller is expected to have prebuilt an identical
// structure first, exactly as ModelPart::load does for its sub-model-parts. A fresh null local made a
// DataOnly load read zero bytes for the container, desynchronizing the stream and discarding whatever
// container the adaptor already held, even when it was correctly preinitialized.
//
// Uses an empty MasterSlaveConstraint container (no entities created) so the test stays focused on
// TensorAdaptor::load()'s own container-pointer handling, without depending on Node/PointerVectorSet's
// own (separate, pre-existing) entity-level DataOnly round-trip behavior.
KRATOS_TEST_CASE_IN_SUITE(TensorAdaptorDataOnlyLoadReusesExistingContainer, KratosCoreFastSuite)
{
    Model model;
    ModelPart& r_model_part = model.CreateModelPart("Test");
    auto p_constraints = r_model_part.GetMesh().pMasterSlaveConstraints();

    VariableTensorAdaptor original(p_constraints, &PRESSURE);

    // `loaded` is preinitialized against the same live (empty) container, as a restart reload
    // would build it before a DataOnly load.
    VariableTensorAdaptor loaded(p_constraints, &PRESSURE);

    const std::string file_name = "test_tensor_adaptor_data_only_load";
    // Pre-create the backing file: FileSerializer opens in|out and only falls back to a
    // (write-only) out-only stream if that fails, which would make the immediately-following
    // load() read garbage.
    std::fstream(file_name + ".rest", std::ios::out).close();
    {
        FileSerializer serializer(file_name, Serializer::SERIALIZER_NO_TRACE, /*DataOnly=*/true);
        serializer.save("TA", original);
        // Unlike StreamSerializer's in-memory stringstream (whose independent get/put positions leave
        // the read position at 0 after a fresh write), FileSerializer's fstream shares one file-level
        // position between reads and writes -- switching from writing to reading requires an explicit
        // seek back to the start.
        serializer.SetLoadState();
        serializer.load("TA", loaded);
    }
    std::remove((file_name + ".rest").c_str());

    KRATOS_EXPECT_EQ(std::get<ModelPart::MasterSlaveConstraintContainerType::Pointer>(loaded.GetContainer()).get(), p_constraints.get());
}

// Regression for CombinedTensorAdaptor::load(): resize(size) preserves any preexisting child pointers,
// and Serializer::load(shared_ptr&) only reconstructs the registered type when the destination pointer
// is null -- otherwise it dispatches virtually through the pointer's *existing* dynamic type. Loading a
// stream saved from one derived type into a stale child of a *different* derived type therefore ran the
// wrong subtype's load() and would have corrupted that child instead of recovering the type actually
// saved. Reproduces in the default (non-DataOnly) mode, via the constructor also used for `original`
// in CombinedTensorAdaptorWithDerivedChildrenSerialization above.
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
    // `loaded` starts with a stale child of the WRONG derived type in the slot the stream actually
    // saved a HistoricalVariableTensorAdaptor to -- before the fix, resize(size) would keep this
    // pointer non-null and Serializer::load(shared_ptr&) would dispatch through its (wrong) dynamic
    // type instead of reconstructing the correct one.
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

} // namespace Kratos::Testing
