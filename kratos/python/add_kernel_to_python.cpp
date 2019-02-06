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

// Project includes
#include "includes/define_python.h"
#include "includes/kernel.h"
#include "python/add_kernel_to_python.h"

// System includes
#include <sstream>

namespace Kratos {
namespace Python {

bool HasFlag(Kernel& rKernel, const std::string& flag_name)
{
    return KratosComponents<Flags>::Has(flag_name);
}

Flags GetFlag(Kernel& rKernel, const std::string& flag_name)
{
    return KratosComponents<Flags>::Get(flag_name);
}

template <class TVariableType>
bool HasVariable(Kernel& rKernel, const std::string& variable_name)
{
    return KratosComponents<TVariableType>::Has(variable_name);
}

template <class TVariableType>
const TVariableType& GetVariable(Kernel& rKernel, const std::string& variable_name)
{
    if (KratosComponents<TVariableType>::Has(variable_name)) {
        return KratosComponents<TVariableType>::Get(variable_name);
    }

    return TVariableType::StaticObject();
}

bool HasConstitutiveLaw(Kernel& rKernel, const std::string& constitutive_law_name)
{
    return KratosComponents<ConstitutiveLaw>::Has(constitutive_law_name);
}

const ConstitutiveLaw& GetConstitutiveLaw(Kernel& rKernel, const std::string& constitutive_law_name)
{
    return KratosComponents<ConstitutiveLaw>::Get(constitutive_law_name);
}

template <class TVariableType>
void PrintVariablesName(Kernel& rKernel)
{
    KratosComponents<TVariableType> kratos_components;
    kratos_components.PrintData(std::cout);
}

template <class TVariableType>
std::string GetVariableNames(Kernel& rKernel)
{
    KratosComponents<TVariableType> kratos_components;
    std::stringstream buffer;
    kratos_components.PrintData(buffer);
    return buffer.str();
}

void RegisterInPythonKernelVariables()
{
    auto comp = KratosComponents<VariableData>::GetComponents();
    auto m = pybind11::module::import("KratosMultiphysics"); //Note that this is added to KratosMultiphysics not to

    for(auto item = comp.begin(); item!=comp.end(); item++) {
        auto& var = (item->second);
        std::string name = item->first;

        m.attr(name.c_str()) = var;
    }
}

void RegisterInPythonApplicationVariables(KratosApplication& Application)
{
    auto comp = KratosComponents<VariableData>::GetComponents();
    auto kernel_module = pybind11::module::import("KratosMultiphysics");
    auto app_module = pybind11::module::import((std::string("KratosMultiphysics.")+Application.Name()).c_str());

    KRATOS_INFO("")
    << "****************************************" << std::endl
    << "Application Name" << Application.Name() << std::endl
    << "****************************************" << std::endl;

    for(auto item = comp.begin(); item!=comp.end(); item++) {
        auto& var = (item->second);
        std::string var_name = item->first;
        KRATOS_INFO("Variable Name") << var_name << std::endl;

        if(! hasattr(kernel_module,var_name.c_str()) ) //variable not present in kernel
            app_module.attr(var_name.c_str()) = var;
    }
}

void AddKernelToPython(pybind11::module& m)
{
    namespace py = pybind11;

    py::class_<Kernel, Kernel::Pointer>(m,"Kernel")
        .def(py::init<>())
        .def("Initialize", [](Kernel& self){ self.Initialize();
        /*RegisterInPythonKernelVariables();*/ }) //&Kernel::Initialize)
        .def("ImportApplication", &Kernel::ImportApplication)
        .def("InitializeApplication",  [](Kernel& self, KratosApplication& App){ self.Initialize();
        /*RegisterInPythonApplicationVariables(App);*/ }) //&Kernel::InitializeApplication)
        //.def(""A,&Kernel::Initialize)
        .def("IsImported", &Kernel::IsImported)
        .def("HasFlag", HasFlag)
        .def("GetFlag", GetFlag)
        .def("HasBoolVariable", HasVariable<Variable<bool> >)
        .def("GetBoolVariable", GetVariable<Variable<bool> >, py::return_value_policy::reference_internal)
        .def("HasIntVariable", HasVariable<Variable<int> >)
        .def("GetIntVariable", GetVariable<Variable<int> >, py::return_value_policy::reference_internal)
        .def("HasUnsignedIntVariable", HasVariable<Variable<unsigned int> >)
        .def("GetUnsignedIntVariable", GetVariable<Variable<unsigned int> >, py::return_value_policy::reference_internal)
        .def("HasDoubleVariable", HasVariable<Variable<double> >)
        .def("GetDoubleVariable", GetVariable<Variable<double> >, py::return_value_policy::reference_internal)
        .def("HasArrayVariable", HasVariable<Variable<array_1d<double, 3> > >)
        .def("HasArray4Variable", HasVariable<Variable<array_1d<double, 4> > >)
        .def("HasArray6Variable", HasVariable<Variable<array_1d<double, 6> > >)
        .def("HasArray9Variable", HasVariable<Variable<array_1d<double, 9> > >)
        .def("GetArrayVariable", GetVariable<Variable<array_1d<double, 3> > >, py::return_value_policy::reference_internal)
        .def("GetArray4Variable", GetVariable<Variable<array_1d<double, 4> > >, py::return_value_policy::reference_internal)
        .def("GetArray6Variable", GetVariable<Variable<array_1d<double, 6> > >, py::return_value_policy::reference_internal)
        .def("GetArray9Variable", GetVariable<Variable<array_1d<double, 9> > >, py::return_value_policy::reference_internal)
        .def("HasVectorVariable", HasVariable<Variable<Vector> >)
        .def("GetVectorVariable", GetVariable<Variable<Vector> >, py::return_value_policy::reference_internal)
        .def("HasMatrixVariable", HasVariable<Variable<Matrix> >)
        .def("GetMatrixVariable", GetVariable<Variable<Matrix> >, py::return_value_policy::reference_internal)
        .def("HasStringVariable", HasVariable<Variable<std::string> >)
        .def("GetStringVariable", GetVariable<Variable<std::string> >, py::return_value_policy::reference_internal)
        .def("HasVariableComponent",HasVariable<VariableComponent< VectorComponentAdaptor<array_1d<double, 3> > > >)
        .def("HasVariableComponent4",HasVariable<VariableComponent< VectorComponentAdaptor<array_1d<double, 4> > > >)
        .def("HasVariableComponent6",HasVariable<VariableComponent< VectorComponentAdaptor<array_1d<double, 6> > > >)
        .def("HasVariableComponent9",HasVariable<VariableComponent< VectorComponentAdaptor<array_1d<double, 9> > > >)
        .def("GetVariableComponent", GetVariable<VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > > >,py::return_value_policy::reference_internal)
        .def("GetVariableComponent4", GetVariable<VariableComponent<VectorComponentAdaptor<array_1d<double, 4> > > >,py::return_value_policy::reference_internal)
        .def("GetVariableComponent6", GetVariable<VariableComponent<VectorComponentAdaptor<array_1d<double, 6> > > >,py::return_value_policy::reference_internal)
        .def("GetVariableComponent9", GetVariable<VariableComponent<VectorComponentAdaptor<array_1d<double, 9> > > >,py::return_value_policy::reference_internal)
        .def("HasFlagsVariable", HasVariable<Variable<Flags> >)
        .def("GetFlagsVariable", GetVariable<Variable<Flags> >, py::return_value_policy::reference_internal)
        .def("HasVariableData", HasVariable<VariableData>)
        .def("PrintAllVariables", PrintVariablesName<VariableData>)
        .def("PrintBoolVariables", PrintVariablesName<Variable<bool> >)
        .def("PrintIntVariables", PrintVariablesName<Variable<int> >)
        .def("PrintUnsignedIntVariables", PrintVariablesName<Variable<int> >)
        .def("PrintDoubleVariables", PrintVariablesName<Variable<double> >)
        .def("PrintArrayVariables",PrintVariablesName<Variable<array_1d<double, 3> > >)
        .def("PrintArray4Variables",PrintVariablesName<Variable<array_1d<double, 4> > >)
        .def("PrintArray6Variables",PrintVariablesName<Variable<array_1d<double, 6> > >)
        .def("PrintArray9Variables",PrintVariablesName<Variable<array_1d<double, 9> > >)
        .def("PrintVectorVariables", PrintVariablesName<Variable<Vector> >)
        .def("PrintMatrixVariables", PrintVariablesName<Variable<Matrix> >)
        .def("PrintStringVariables", PrintVariablesName<Variable<std::string> >)
        .def("PrintFlagsVariables", PrintVariablesName<Variable<Flags> >)
        .def("PrintVariableComponentVariables",PrintVariablesName<VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > > >)
        .def("PrintVariableComponent4Variables",PrintVariablesName<VariableComponent<VectorComponentAdaptor<array_1d<double, 4> > > >)
        .def("PrintVariableComponent6Variables",PrintVariablesName<VariableComponent<VectorComponentAdaptor<array_1d<double, 6> > > >)
        .def("PrintVariableComponent9Variables",PrintVariablesName<VariableComponent<VectorComponentAdaptor<array_1d<double, 9> > > >)
        .def("GetAllVariableNames", GetVariableNames<VariableData>)
        .def("GetBoolVariableNames", GetVariableNames<Variable<bool> >)
        .def("GetIntVariableNames", GetVariableNames<Variable<int> >)
        .def("GetUnsignedIntVariableNames", GetVariableNames<Variable<int> >)
        .def("GetDoubleVariableNames", GetVariableNames<Variable<double> >)
        .def("GetArrayVariableNames",GetVariableNames<Variable<array_1d<double, 3> > >)
        .def("GetArrayVariableNames",GetVariableNames<Variable<array_1d<double, 4> > >)
        .def("GetArrayVariableNames",GetVariableNames<Variable<array_1d<double, 6> > >)
        .def("GetArrayVariableNames",GetVariableNames<Variable<array_1d<double, 9> > >)
        .def("GetVectorVariableNames", GetVariableNames<Variable<Vector> >)
        .def("GetMatrixVariableNames", GetVariableNames<Variable<Matrix> >)
        .def("GetStringVariableNames", GetVariableNames<Variable<std::string> >)
        .def("GetFlagsVariableNames", GetVariableNames<Variable<Flags> >)
        .def("GetVariableComponentVariableNames", GetVariableNames<VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > > >)
        .def("GetVariableComponentVariable4Names", GetVariableNames<VariableComponent<VectorComponentAdaptor<array_1d<double, 4> > > >)
        .def("GetVariableComponentVariable6Names", GetVariableNames<VariableComponent<VectorComponentAdaptor<array_1d<double, 6> > > >)
        .def("GetVariableComponentVariable9Names", GetVariableNames<VariableComponent<VectorComponentAdaptor<array_1d<double, 9> > > >)
        .def("__str__", PrintObject<Kernel>)
        .def("HasConstitutiveLaw", HasConstitutiveLaw)
        .def("GetConstitutiveLaw", GetConstitutiveLaw, py::return_value_policy::reference_internal)
        .def_static("Version", &Kernel::Version)
        .def_static("BuildType", &Kernel::BuildType)
        ;
}

}  // namespace Python.
}  // Namespace Kratos
