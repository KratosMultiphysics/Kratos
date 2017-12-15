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
//                   Pooyan Dadvand
//

// External includes
#include <boost/python.hpp>
#include "boost/python/suite/indexing/map_indexing_suite.hpp"

// System includes
#include <sstream>

// Project includes
#include "includes/define.h"
#include "includes/kernel.h"
#include "python/add_kernel_to_python.h"

namespace Kratos {
namespace Python {
using namespace boost::python;

template <class TVariableType>
bool HasVariable(Kernel& rKernel, const std::string& variable_name) {
    return KratosComponents<TVariableType>::Has(variable_name);
}

template <class TVariableType>
const TVariableType& GetVariable(
    Kernel& rKernel, const std::string& variable_name) {
    if (KratosComponents<TVariableType>::Has(variable_name)) {
        return KratosComponents<TVariableType>::Get(variable_name);
    }

    return TVariableType::StaticObject();
}

template <class TVariableType>
void PrintVariablesName(Kernel& rKernel) {
    KratosComponents<TVariableType> kratos_components;
    kratos_components.PrintData(std::cout);
}
template <class TVariableType>
std::string GetVariableNames(Kernel& rKernel) {
    KratosComponents<TVariableType> kratos_components;
    std::stringstream buffer;
    kratos_components.PrintData(buffer);
    return buffer.str();
}

void AddKernelToPython() {
    //class_<std::map<std::string, const VariableData*> >("VariableDataMap")
    // .def(map_indexing_suite<std::map<std::string, const VariableData> >())
    //  ;
    class_<Kernel, Kernel::Pointer, boost::noncopyable>("Kernel")
        .def("Initialize", &Kernel::Initialize)
        .def("ImportApplication", &Kernel::ImportApplication)
        .def("InitializeApplication", &Kernel::InitializeApplication)
        //.def(""A,&Kernel::Initialize)
        .def("IsImported", &Kernel::IsImported)
        .def("HasBoolVariable", HasVariable<Variable<bool> >)
        .def("GetBoolVariable", GetVariable<Variable<bool> >,
            return_internal_reference<>())
        .def("HasIntVariable", HasVariable<Variable<int> >)
        .def("GetIntVariable", GetVariable<Variable<int> >,
            return_internal_reference<>())
        .def("HasUnsignedIntVariable", HasVariable<Variable<unsigned int> >)
        .def("GetUnsignedIntVariable", GetVariable<Variable<unsigned int> >,
            return_internal_reference<>())
        .def("HasDoubleVariable", HasVariable<Variable<double> >)
        .def("GetDoubleVariable", GetVariable<Variable<double> >,
            return_internal_reference<>())
        .def("HasArrayVariable", HasVariable<Variable<array_1d<double, 3> > >)
        .def("GetArrayVariable", GetVariable<Variable<array_1d<double, 3> > >,
            return_internal_reference<>())
        .def("HasVectorVariable", HasVariable<Variable<Vector> >)
        .def("GetVectorVariable", GetVariable<Variable<Vector> >,
            return_internal_reference<>())
        .def("HasMatrixVariable", HasVariable<Variable<Matrix> >)
        .def("GetMatrixVariable", GetVariable<Variable<Matrix> >,
            return_internal_reference<>())
        .def("HasStringVariable", HasVariable<Variable<std::string> >)
        .def("GetStringVariable", GetVariable<Variable<std::string> >,
            return_internal_reference<>())
        .def("HasVariableComponent",
            HasVariable<VariableComponent<
                VectorComponentAdaptor<array_1d<double, 3> > > >)
        .def("GetVariableComponent",
            GetVariable<VariableComponent<
                VectorComponentAdaptor<array_1d<double, 3> > > >,
            return_internal_reference<>())
        .def("HasFlagsVariable", HasVariable<Variable<Flags> >)
        .def("GetFlagsVariable", GetVariable<Variable<Flags> >,
            return_internal_reference<>())
        .def("HasVariableData", HasVariable<VariableData>)
        .def("PrintAllVariables", PrintVariablesName<VariableData>)
        .def("PrintBoolVariables", PrintVariablesName<Variable<bool> >)
        .def("PrintIntVariables", PrintVariablesName<Variable<int> >)
        .def("PrintUnsignedIntVariables", PrintVariablesName<Variable<int> >)
        .def("PrintDoubleVariables", PrintVariablesName<Variable<double> >)
        .def("PrintArrayVariables",
            PrintVariablesName<Variable<array_1d<double, 3> > >)
        .def("PrintVectorVariables", PrintVariablesName<Variable<Vector> >)
        .def("PrintMatrixVariables", PrintVariablesName<Variable<Matrix> >)
        .def("PrintStringVariables", PrintVariablesName<Variable<std::string> >)
        .def("PrintFlagsVariables", PrintVariablesName<Variable<Flags> >)
        .def("PrintVariableComponentVariables",
            PrintVariablesName<VariableComponent<
                VectorComponentAdaptor<array_1d<double, 3> > > >)
        .def("GetAllVariableNames", GetVariableNames<VariableData>)
        .def("GetBoolVariableNames", GetVariableNames<Variable<bool> >)
        .def("GetIntVariableNames", GetVariableNames<Variable<int> >)
        .def("GetUnsignedIntVariableNames", GetVariableNames<Variable<int> >)
        .def("GetDoubleVariableNames", GetVariableNames<Variable<double> >)
        .def("GetArrayVariableNames",
            GetVariableNames<Variable<array_1d<double, 3> > >)
        .def("GetVectorVariableNames", GetVariableNames<Variable<Vector> >)
        .def("GetMatrixVariableNames", GetVariableNames<Variable<Matrix> >)
        .def("GetStringVariableNames", GetVariableNames<Variable<std::string> >)
        .def("GetFlagsVariableNames", GetVariableNames<Variable<Flags> >)
        .def("GetVariableComponentVariableNames",
            GetVariableNames<VariableComponent<
                VectorComponentAdaptor<array_1d<double, 3> > > >)
        .def(self_ns::str(self));
}

}  // namespace Python.

}  // Namespace Kratos
