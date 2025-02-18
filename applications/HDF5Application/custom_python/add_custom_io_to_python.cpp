//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela
//                   Suneth Warnakulasuriya
//

// System includes

// External includes
#include "pybind11/stl.h"
#include "pybind11/pybind11.h"

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/io.h"
#include "includes/kratos_parameters.h"
#include "includes/communicator.h"
#include "includes/parallel_environment.h"

// Application includes
#include "custom_io/hdf5_file.h"
#include "custom_io/hdf5_model_part_io.h"
#include "custom_io/hdf5_data_value_container_io.h"
#include "custom_io/hdf5_vertex_container_io.h"
#include "custom_io/hdf5_container_component_io.h"

#include "custom_utilities/container_io_utils.h"

#ifdef KRATOS_USING_MPI
#include "custom_io/hdf5_partitioned_model_part_io.h"
#endif

namespace Kratos {
namespace Python {

template<class TContainerType, class TContainerDataIOType>
struct VariableContainerComponentIOWrapper
{
    using ContainerType = TContainerType;

    using ContainerDataIOType = TContainerDataIOType;

    using ContainerIOType = typename HDF5::ContainerComponentIO<
                                        TContainerType,
                                        TContainerDataIOType,
                                        Variable<int>,
                                        Variable<double>,
                                        Variable<array_1d<double, 3>>,
                                        Variable<array_1d<double, 4>>,
                                        Variable<array_1d<double, 6>>,
                                        Variable<array_1d<double, 9>>,
                                        Variable<Kratos::Vector>,
                                        Variable<Kratos::Matrix>>;

};

template<class TContainerType>
struct FlagContainerComponentIOWrapper
{
    using ContainerType = TContainerType;

    using ContainerDataIOType = HDF5::Internals::FlagIO;

    using ContainerIOType = typename HDF5::ContainerComponentIO<
                                        TContainerType,
                                        ContainerDataIOType,
                                        Flags>;

};

template<class TContainerComponentIOWrapper>
void AddContainerComponentIOToPython(
    pybind11::module& m,
    const std::string& rName)
{
    namespace py = pybind11;

    using data_io = typename TContainerComponentIOWrapper::ContainerIOType;

    using container_type = typename TContainerComponentIOWrapper::ContainerType;

    using container_data_io_type = typename TContainerComponentIOWrapper::ContainerDataIOType;

    std::string container_name, arg_name;
    if constexpr(std::is_same_v<container_type, ModelPart::NodesContainerType>) {
        container_name = "Nodal";
        arg_name = "nodes_container";
    } else if constexpr(std::is_same_v<container_type, ModelPart::ConditionsContainerType>) {
        container_name = "Condition";
        arg_name = "conditions_container";
    } else if constexpr(std::is_same_v<container_type, ModelPart::ElementsContainerType>) {
        container_name = "Element";
        arg_name = "elements_container";
    } else {
        static_assert(!std::is_same_v<TContainerComponentIOWrapper, TContainerComponentIOWrapper>, "Unsupported container component io wrapper type.");
    }

    // TODO: Remove legacy suffix once the python side is fixed.
    std::string results_type = "Results", legacy_suffix = "Data";
    if constexpr(std::is_same_v<container_data_io_type, HDF5::Internals::FlagIO>) {
        results_type = "Flags";
        legacy_suffix = "Flag";
    }

    py::class_<data_io, typename data_io::Pointer>(m, rName.c_str())
        // TODO: Remove th below custom init, and use the one without suffix.
        .def(py::init([container_name, legacy_suffix](Parameters Settings, HDF5::File::Pointer pFile) {
            return Kratos::make_shared<data_io>(Settings, pFile, ("/" + container_name + legacy_suffix + "Values/"));
        }))
        // .def(py::init<Parameters, HDF5::File::Pointer>(), py::arg("settings"), py::arg("hdf5_file"))
        .def("Write", [](data_io& rSelf, const ModelPart& rModelPart, const Parameters Attributes) {
                    rSelf.Write(rModelPart, container_data_io_type{}, Attributes);
                },
                py::arg("model_part"),
                py::arg("attributes") = Parameters("""{}"""))
        .def("Read", [](data_io& rSelf, ModelPart& rModelPart) {
                    return rSelf.Read(rModelPart, container_data_io_type{});
                },
                py::arg("model_part"))
        .def("ReadAttributes", &data_io::ReadAttributes)
        .def(("Write" + container_name + results_type).c_str(), [rName, container_name, results_type](data_io& rSelf, const container_type& rContainer) {
                KRATOS_WARNING("DEPRECATION") << "Using deprecated \"Write" << container_name << results_type << "\" method in \"" << rName << "\". Please use \"Write\" method instead.\n";
                rSelf.Write(rContainer, container_data_io_type{}, Parameters("""{}"""));
            }
            , py::arg(arg_name.c_str()))
        .def(("Read" + container_name + results_type).c_str(), [rName, container_name, results_type](data_io& rSelf, container_type& rContainer, Communicator& rCommunicator) {
                KRATOS_WARNING("DEPRECATION") << "Using deprecated \"Read" << container_name << results_type << "\" method in \"" << rName << "\". Please use \"Read\" method instead.\n";
                return rSelf.Read(rContainer, container_data_io_type{}, rCommunicator);
            },
            py::arg(arg_name.c_str()),
            py::arg("communicator"))
        ;
}

class PythonNodalSolutionStepBossakIO : public HDF5::ContainerComponentIO<
                                                                ModelPart::NodesContainerType,
                                                                HDF5::Internals::BossakIO,
                                                                Variable<int>,
                                                                Variable<double>,
                                                                Variable<array_1d<double, 3>>,
                                                                Variable<array_1d<double, 4>>,
                                                                Variable<array_1d<double, 6>>,
                                                                Variable<array_1d<double, 9>>,
                                                                Variable<Kratos::Vector>,
                                                                Variable<Kratos::Matrix>>
{
public:
    ///@name Type Definitions
    ///@{

    using BaseType = HDF5::ContainerComponentIO<
                            ModelPart::NodesContainerType,
                            HDF5::Internals::BossakIO,
                            Variable<int>,
                            Variable<double>,
                            Variable<array_1d<double, 3>>,
                            Variable<array_1d<double, 4>>,
                            Variable<array_1d<double, 6>>,
                            Variable<array_1d<double, 9>>,
                            Variable<Kratos::Vector>,
                            Variable<Kratos::Matrix>>;

    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(PythonNodalSolutionStepBossakIO);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    PythonNodalSolutionStepBossakIO(
        Parameters Settings,
        HDF5::File::Pointer pFile)
        : BaseType(Settings, pFile, "/NodalSolutionStepData/"),
          mAlphaBossak(0.0)
    {
    }

    ///@}
    ///@name Operations
    ///@{

    void SetAlphaBossak(const double AlphaBossak)
    {
        mAlphaBossak = AlphaBossak;
    }

    void Write(const ModelPart& rModelPart)
    {
        BaseType::Write(rModelPart, HDF5::Internals::BossakIO(mAlphaBossak), Parameters("""{}"""));
    }

    void Read(ModelPart& rModelPart)
    {
        BaseType::Read(rModelPart, HDF5::Internals::BossakIO(mAlphaBossak));
    }

    ///@}

private:
    ///@name Private member variables
    ///@{

    double mAlphaBossak;

    ///@}

};

void AddCustomIOToPython(pybind11::module& m)
{
    namespace py = pybind11;

    m.def("WriteDataValueContainer", &HDF5::Internals::WriteDataValueContainer, "");
    m.def("ReadDataValueContainer", &HDF5::Internals::ReadDataValueContainer, "");

    py::class_<HDF5::File, HDF5::File::Pointer >(m,"HDF5File")
        .def(py::init([](Parameters Settings) {
            KRATOS_WARNING("DEPRECATION") << "Using deprecated constructor in \"HDF5File\". Please use constructor with DataCommunicator.\n";
            return Kratos::make_shared<HDF5::File>(ParallelEnvironment::GetDefaultDataCommunicator(), Settings);
        }))
        .def(py::init<const DataCommunicator&, Parameters>())
        .def("HasPath",&HDF5::File::HasPath)
        .def("IsGroup",&HDF5::File::IsGroup)
        .def("IsDataSet",&HDF5::File::IsDataSet)
        .def("CreateGroup",&HDF5::File::CreateGroup)
        .def("AddPath",&HDF5::File::AddPath)
        .def("GetDataDimensions",&HDF5::File::GetDataDimensions)
        .def("HasIntDataType",&HDF5::File::HasIntDataType)
        .def("HasFloatDataType",&HDF5::File::HasFloatDataType)
        .def("Flush",&HDF5::File::Flush)
        .def("Close", &HDF5::File::Close)
        .def("GetFileSize",&HDF5::File::GetFileSize)
        .def("GetFileName",&HDF5::File::GetFileName)
        ;

    py::class_<HDF5::ModelPartIO, HDF5::ModelPartIO::Pointer, IO>(m,"HDF5ModelPartIO")
        .def(py::init<HDF5::File::Pointer, std::string const&>())
        ;

    using nodal_solution_step_data_io = VariableContainerComponentIOWrapper<ModelPart::NodesContainerType, HDF5::Internals::HistoricalIO>::ContainerIOType;
    py::class_<nodal_solution_step_data_io, nodal_solution_step_data_io::Pointer>(m,"HDF5NodalSolutionStepDataIO")
        // TODO: Remove th below custom init, and use the one without suffix.
        .def(py::init([](Parameters Settings, HDF5::File::Pointer pFile) {
            return Kratos::make_shared<nodal_solution_step_data_io>(Settings, pFile, "/NodalSolutionStepData/");
        }))
        // .def(py::init<Parameters, HDF5::File::Pointer>(), py::arg("settings"), py::arg("hdf5_file"))
        .def("Write", [](nodal_solution_step_data_io& rSelf, const ModelPart& rModelPart, const Parameters Attributes, const IndexType StepIndex) {
                rSelf.Write(rModelPart, HDF5::Internals::HistoricalIO(StepIndex), Attributes);
            },
            py::arg("model_part"),
            py::arg("attributes") = Parameters("""{}"""),
            py::arg("step_index") = 0)
        .def("Read", [](nodal_solution_step_data_io& rSelf, ModelPart& rModelPart, const IndexType StepIndex) {
                return rSelf.Read(rModelPart, HDF5::Internals::HistoricalIO(StepIndex));
            },
            py::arg("model_part"),
            py::arg("step_index") = 0)
        .def("ReadAttributes", &nodal_solution_step_data_io::ReadAttributes)
        .def("WriteNodalResults", [](nodal_solution_step_data_io& rSelf, const ModelPart& rModelPart, const IndexType StepIndex) {
                KRATOS_WARNING("DEPRECATION") << "Using deprecated \"WriteNodalResults\" method in \"HDF5NodalSolutionStepDataIO\". Please use \"Write\" method instead.\n";
                rSelf.Write(rModelPart, HDF5::Internals::HistoricalIO(StepIndex), Parameters("""{}"""));
            },
            py::arg("model_part"),
            py::arg("step_index") = 0)
        .def("ReadNodalResults", [](nodal_solution_step_data_io& rSelf, ModelPart& rModelPart, const IndexType StepIndex) {
                KRATOS_WARNING("DEPRECATION") << "Using deprecated \"ReadNodalResults\" method in \"HDF5NodalSolutionStepDataIO\". Please use \"Read\" method instead.\n";
                return rSelf.Read(rModelPart, HDF5::Internals::HistoricalIO(StepIndex));
            },
            py::arg("model_part"),
            py::arg("step_index") = 0)
        ;

    py::class_<PythonNodalSolutionStepBossakIO, PythonNodalSolutionStepBossakIO::Pointer>(m,"HDF5NodalSolutionStepBossakIO")
        .def(py::init<Parameters, HDF5::File::Pointer>(), py::arg("settings"), py::arg("hdf5_file"))
        .def("WriteNodalResults", [](PythonNodalSolutionStepBossakIO& rSelf, const ModelPart& rModelPart) {
                KRATOS_WARNING("DEPRECATION") << "Using deprecated \"WriteNodalResults\" method in \"HDF5NodalSolutionStepBossakIO\". Please use \"Write\" method instead.\n";
                rSelf.Write(rModelPart);
            },
            py::arg("model_part"))
        .def("ReadNodalResults", [](PythonNodalSolutionStepBossakIO& rSelf, ModelPart& rModelPart) {
                KRATOS_WARNING("DEPRECATION") << "Using deprecated \"ReadNodalResults\" method in \"HDF5NodalSolutionStepBossakIO\". Please use \"Read\" method instead.\n";
                return rSelf.Read(rModelPart);
            },
            py::arg("model_part"))
        .def("Write", &PythonNodalSolutionStepBossakIO::Write, py::arg("model_part"))
        .def("Read", &PythonNodalSolutionStepBossakIO::Read, py::arg("model_part"))
        .def("SetAlphaBossak", &PythonNodalSolutionStepBossakIO::SetAlphaBossak, py::arg("alpha_bossak"))
        ;

    AddContainerComponentIOToPython<FlagContainerComponentIOWrapper<ModelPart::ElementsContainerType>>(m, "HDF5ElementFlagValueIO");
    AddContainerComponentIOToPython<VariableContainerComponentIOWrapper<ModelPart::ElementsContainerType, HDF5::Internals::NonHistoricalIO>>(m, "HDF5ElementDataValueIO");
    AddContainerComponentIOToPython<FlagContainerComponentIOWrapper<ModelPart::ConditionsContainerType>>(m, "HDF5ConditionFlagValueIO");
    AddContainerComponentIOToPython<VariableContainerComponentIOWrapper<ModelPart::ConditionsContainerType, HDF5::Internals::NonHistoricalIO>>(m, "HDF5ConditionDataValueIO");
    AddContainerComponentIOToPython<FlagContainerComponentIOWrapper<ModelPart::NodesContainerType>>(m, "HDF5NodalFlagValueIO");
    AddContainerComponentIOToPython<VariableContainerComponentIOWrapper<ModelPart::NodesContainerType, HDF5::Internals::NonHistoricalIO>>(m, "HDF5NodalDataValueIO");

    using element_gauss_io = VariableContainerComponentIOWrapper<ModelPart::ElementsContainerType, HDF5::Internals::GaussPointIO>::ContainerIOType;
    py::class_<element_gauss_io, element_gauss_io::Pointer>(m,"HDF5ElementGaussPointIO")
        // TODO: Remove th below custom init, and use the one without suffix.
        .def(py::init([](Parameters Settings, HDF5::File::Pointer pFile) {
            return Kratos::make_shared<element_gauss_io>(Settings, pFile, "/ElementGaussPointValues/");
        }))
        // .def(py::init<Parameters, HDF5::File::Pointer>(), py::arg("settings"), py::arg("hdf5_file"))
        .def("Write", [](element_gauss_io& rSelf, const ModelPart& rModelPart, const Parameters Attributes) {
                rSelf.Write(rModelPart, HDF5::Internals::GaussPointIO(rModelPart.GetProcessInfo()), Attributes);
            },
            py::arg("model_part"),
            py::arg("attributes") = Parameters("""{}"""))
        .def("WriteElementGaussPointValues", [](element_gauss_io& rSelf, const ModelPart::ElementsContainerType& rContainer, const DataCommunicator& rDataCommunicator, const ProcessInfo& rProcessInfo) {
                KRATOS_WARNING("DEPRECATION") << "Using deprecated \"WriteElementGaussPointValues\" method in \"HDF5ElementGaussPointIO\". Please use \"Write\" method instead.\n";
                rSelf.Write(rContainer, HDF5::Internals::GaussPointIO(rProcessInfo), Parameters("""{}"""));
            },
            py::arg("elements"),
            py::arg("data_communicator"),
            py::arg("process_info"))
        ;

    using condition_gauss_io = VariableContainerComponentIOWrapper<ModelPart::ConditionsContainerType, HDF5::Internals::GaussPointIO>::ContainerIOType;
    py::class_<condition_gauss_io, condition_gauss_io::Pointer>(m,"HDF5ConditionGaussPointIO")
        // TODO: Remove th below custom init, and use the one without suffix.
        .def(py::init([](Parameters Settings, HDF5::File::Pointer pFile) {
            return Kratos::make_shared<condition_gauss_io>(Settings, pFile, "/ConditionGaussPointValues/");
        }))
        // .def(py::init<Parameters, HDF5::File::Pointer>(), py::arg("settings"), py::arg("hdf5_file"))
        .def("Write", [](condition_gauss_io& rSelf, const ModelPart& rModelPart, const Parameters Attributes) {
                rSelf.Write(rModelPart, HDF5::Internals::GaussPointIO(rModelPart.GetProcessInfo()), Attributes);
            },
            py::arg("model_part"),
            py::arg("attributes") = Parameters("""{}"""))
        .def("WriteConditionGaussPointValues", [](condition_gauss_io& rSelf, const ModelPart::ConditionsContainerType& rContainer, const DataCommunicator& rDataCommunicator, const ProcessInfo& rProcessInfo) {
                KRATOS_WARNING("DEPRECATION") << "Using deprecated \"WriteConditionGaussPointValues\" method in \"HDF5ConditionGaussPointIO\". Please use \"Write\" method instead.\n";
                rSelf.Write(rContainer, HDF5::Internals::GaussPointIO(rProcessInfo), Parameters("""{}"""));
            },
            py::arg("conditions"),
            py::arg("data_communicator"),
            py::arg("process_info"))
        ;

    py::class_<HDF5::VertexContainerCoordinateIO, HDF5::VertexContainerCoordinateIO::Pointer>(m, "VertexContainerCoordinateIO")
        .def(py::init<Parameters, HDF5::File::Pointer>(), py::arg("settings"), py::arg("hdf5_file"))
        .def("Write", &HDF5::VertexContainerCoordinateIO::Write, py::arg("vertices"), py::arg("attributes") = Parameters("""{}"""))
        ;

    py::class_<HDF5::VertexContainerVariableIO, HDF5::VertexContainerVariableIO::Pointer>(m, "VertexContainerVariableIO")
        .def(py::init<Parameters, HDF5::File::Pointer>(), py::arg("settings"), py::arg("hdf5_file"))
        .def("Write", &HDF5::VertexContainerVariableIO::Write, py::arg("vertices"), py::arg("attributes") = Parameters("""{}"""))
        ;

#ifdef KRATOS_USING_MPI
    py::class_<HDF5::PartitionedModelPartIO, HDF5::PartitionedModelPartIO::Pointer, HDF5::ModelPartIO>
        (m,"HDF5PartitionedModelPartIO")
        .def(py::init<HDF5::File::Pointer, std::string const&>())
        ;
#endif

}

} // namespace Python.

} // Namespace Kratos
