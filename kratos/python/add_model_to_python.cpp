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
//                   Pooyan Dadvand
//

// System includes

// External includes

// Project includes
#include "includes/define_python.h"
#include "containers/model.h"
#include "python/add_model_to_python.h"
#include "includes/convection_diffusion_settings.h"
#include "includes/radiation_settings.h"

namespace Kratos::Python {

ModelPart& Model_GetModelPart(Model& rModel, const std::string& rFullModelPartName)
{
    return rModel.GetModelPart(rFullModelPartName);
}

template< class TBinderType, typename TVariableType >
void VariableIndexingUtility(TBinderType& binder)
{
    binder.def("Has", [](const Model& rModel, const TVariableType& rV){return rModel.Has(rV);} );
    binder.def("SetValue", [](Model& rModel, const TVariableType& rV, const typename TVariableType::Type rValue){rModel.SetValue(rV, rValue);} );
    binder.def("GetValue", [](Model& rModel, const TVariableType& rV){return rModel.GetValue(rV);} );
}

void  AddModelToPython(pybind11::module& m)
{
    using Array1DVariable3 = Variable<array_1d<double, 3>>;
    using Array1DVariable4 = Variable<array_1d<double, 4>>;
    using Array1DVariable6 = Variable<array_1d<double, 6>>;
    using Array1DVariable9 = Variable<array_1d<double, 9>>;

    namespace py = pybind11;
    using ModelBinder = py::class_<Model, Model::Pointer>;
    auto model_binder = ModelBinder(m, "Model")
        .def(py::init<>())
        .def("Reset", &Model::Reset)
        .def("CreateModelPart", [&](Model &self, const std::string &Name) { return &self.CreateModelPart(Name); }, py::return_value_policy::reference_internal)
        .def("CreateModelPart", [&](Model &self, const std::string &Name, unsigned int BufferSize) { return &self.CreateModelPart(Name, BufferSize); }, py::return_value_policy::reference_internal)
        .def("DeleteModelPart", &Model::DeleteModelPart)
        .def("GetModelPart", &Model_GetModelPart, py::return_value_policy::reference_internal)
        .def("HasModelPart", &Model::HasModelPart)
        .def("GetModelPartNames", &Model::GetModelPartNames)
        .def("__getitem__", &Model_GetModelPart, py::return_value_policy::reference_internal)
        .def("__str__", PrintObject<Model>)
        ;
    VariableIndexingUtility< ModelBinder, Variable<bool> >(model_binder);
    VariableIndexingUtility< ModelBinder, Variable<int> >(model_binder);
    VariableIndexingUtility< ModelBinder, Variable<double> >(model_binder);
    VariableIndexingUtility< ModelBinder, Array1DVariable3 >(model_binder);
    VariableIndexingUtility< ModelBinder, Array1DVariable4 >(model_binder);
    VariableIndexingUtility< ModelBinder, Array1DVariable6 >(model_binder);
    VariableIndexingUtility< ModelBinder, Array1DVariable9 >(model_binder);
    VariableIndexingUtility< ModelBinder, Variable<Vector> >(model_binder);
    VariableIndexingUtility< ModelBinder, Variable<Matrix> >(model_binder);
    VariableIndexingUtility< ModelBinder, Variable<ConvectionDiffusionSettings::Pointer> >(model_binder);
    VariableIndexingUtility< ModelBinder, Variable<RadiationSettings::Pointer> >(model_binder);
    VariableIndexingUtility< ModelBinder, Variable<Quaternion<double>> >(model_binder);
    VariableIndexingUtility< ModelBinder, Variable<std::string> >(model_binder);

}

}  // namespace Kratos::Python.
