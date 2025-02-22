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

// System includes

// External includes

// Project includes
#include "includes/define_python.h"
#include "includes/model_part.h"
#include "python/add_dofs_to_python.h"
#include "utilities/dof_utilities/block_build_dof_array_utility.h"
#include "utilities/dof_utilities/elimination_build_dof_array_utility.h"

namespace Kratos::Python {

void  AddDofsToPython(pybind11::module& m)
{
    namespace py = pybind11;

    py::class_<Dof<double>, std::unique_ptr<Dof<double>, py::nodelete>>(m,"Dof")
        .def("GetVariable", &Dof<double>::GetVariable, py::return_value_policy::reference_internal)
        .def("GetReaction", &Dof<double>::GetReaction, py::return_value_policy::reference_internal)
        .def("Id", &Dof<double>::Id)
        .def("GetSolutionStepValue", [](const Dof<double>& self) {
            return self.GetSolutionStepValue();
        })
        .def("SetSolutionStepValue", [](Dof<double>& self, const double value) {
            return self.GetSolutionStepValue() = value;
        })
        .def_property("EquationId", &Dof<double>::EquationId, &Dof<double>::SetEquationId)
        .def("Fix", &Dof<double>::FixDof)
        .def("Free", &Dof<double>::FreeDof)
        .def("IsFixed", &Dof<double>::IsFixed)
        .def("__lt__", [](Dof<double>& dof1, Dof<double>& dof2){return dof1<dof2; })
        .def("__gt__", [](Dof<double>& dof1, Dof<double>& dof2){return dof1>dof2; })
        .def("__le__", [](Dof<double>& dof1, Dof<double>& dof2){return dof1<=dof2; })
        .def("__ge__", [](Dof<double>& dof1, Dof<double>& dof2){return dof1>=dof2; })
        .def("__eq__", [](Dof<double>& dof1, Dof<double>& dof2){return dof1==dof2; })
        .def("__ne__", [](Dof<double>& dof1, Dof<double>& dof2){return !(dof1==dof2); })
        .def("__str__", PrintObject<Dof<double>>)
        ;

    py::class_<ModelPart::DofsArrayType, ModelPart::DofsArrayType::Pointer>(m, "DofsArrayType")
        .def(py::init<>())
        .def("__len__", [](ModelPart::DofsArrayType &self) {
            return self.size();
        })
        .def("__iter__", [](ModelPart::DofsArrayType &self){
            return py::make_iterator(self.begin(), self.end());
        }, py::keep_alive<0, 1>())
        .def("append", [](ModelPart::DofsArrayType &self, Dof<double> &value) {
            self.push_back(&value);
        })
        .def("unique", [](ModelPart::DofsArrayType &self) {
            self.Unique();
        })
        .def("GetValues", [](ModelPart::DofsArrayType &self) {
            Vector values(self.size());
            int counter = 0;
            for (auto &r_dof : self)
            {
                values[counter++] = r_dof.GetSolutionStepValue();
            }
            return values;
        })
        .def("GetEquationIds", [](ModelPart::DofsArrayType &self) {
            std::vector<std::size_t> values(self.size());
            int counter = 0;
            for (auto &r_dof : self)
            {
                values[counter++] = r_dof.EquationId();
            }
            return values;
        })
        .def("SetValues", [](ModelPart::DofsArrayType &self, const Vector& values) {
            int counter = 0;
            for (auto &r_dof : self)
            {
                r_dof.GetSolutionStepValue() = values[counter++];
            }
        })
        ;

    py::class_<BlockBuildDofArrayUtility>(m, "BlockBuildDofArrayUtility")
        .def_static("SetUpDofArray", [](const ModelPart& rModelPart, BlockBuildDofArrayUtility::DofsArrayType& rDofArray){BlockBuildDofArrayUtility::SetUpDofArray(rModelPart, rDofArray);})
        .def_static("SetUpDofArray", [](const ModelPart& rModelPart, BlockBuildDofArrayUtility::DofsArrayType& rDofArray, const unsigned int EchoLevel){BlockBuildDofArrayUtility::SetUpDofArray(rModelPart, rDofArray, EchoLevel);})
        .def_static("SetUpDofArray", [](const ModelPart& rModelPart, BlockBuildDofArrayUtility::DofsArrayType& rDofArray, const unsigned int EchoLevel, const bool CheckReactionDofs){BlockBuildDofArrayUtility::SetUpDofArray(rModelPart, rDofArray, EchoLevel, CheckReactionDofs);})
    ;

    py::class_<EliminationBuildDofArrayUtility>(m, "EliminationBuildDofArrayUtility")
        .def_static("SetUpDofArray", [](const ModelPart& rModelPart, EliminationBuildDofArrayUtility::DofsArrayType& rDofArray){EliminationBuildDofArrayUtility::SetUpDofArray(rModelPart, rDofArray);})
        .def_static("SetUpDofArray", [](const ModelPart& rModelPart, EliminationBuildDofArrayUtility::DofsArrayType& rDofArray, const unsigned int EchoLevel){EliminationBuildDofArrayUtility::SetUpDofArray(rModelPart, rDofArray, EchoLevel);})
        .def_static("SetUpDofArray", [](const ModelPart& rModelPart, EliminationBuildDofArrayUtility::DofsArrayType& rDofArray, const unsigned int EchoLevel, const bool CheckReactionDofs){EliminationBuildDofArrayUtility::SetUpDofArray(rModelPart, rDofArray, EchoLevel, CheckReactionDofs);})
    ;

}

}  // namespace Kratos::Python.
