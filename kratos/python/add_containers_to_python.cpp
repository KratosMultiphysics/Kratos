
//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//

// System includes

// External includes


// Project includes
#include "includes/define_python.h"
#include "includes/ublas_interface.h"
#include "containers/data_value_container.h"
//#include "containers/hash_data_value_container.h"
#include "containers/variables_list_data_value_container.h"
#include "containers/vector_component_adaptor.h"
#include "containers/flags.h"
//#include "containers/all_variables_data_value_container.h"
#include "includes/kratos_flags.h"
#include "includes/variables.h"
#include "includes/constitutive_law.h"
// #include "python/variable_indexing_python.h"
// #include "python/vector_python_interface.h"
// #include "python/vector_scalar_operator_python.h"
// #include "python/vector_vector_operator_python.h"
// #include "python/bounded_vector_python_interface.h"
#include "python/add_deprecated_variables_to_python.h"
#include "python/add_c2c_variables_to_python.h" //TODO: to be removed eventually
#include "python/add_cfd_variables_to_python.h" //TODO: to be removed eventually
#include "python/add_mesh_moving_variables_to_python.h" //TODO: to be removed eventually
#include "python/add_mapping_variables_to_python.h" //TODO: to be removed eventually
#include "python/add_dem_variables_to_python.h" //TODO: to be removed eventually
#include "python/add_fsi_variables_to_python.h" //TODO: to be removed eventually
#include "python/add_mat_variables_to_python.h" //TODO: to be removed eventually
#include "python/add_legacy_structural_app_vars_to_python.h" //TODO: to be removed eventually

#include "includes/convection_diffusion_settings.h"
#include "includes/radiation_settings.h"
#include "utilities/timer.h"
#include "utilities/quaternion.h"


namespace Kratos
{

namespace Python
{
namespace py = pybind11;

Flags FlagsOr(const Flags& Left, const Flags& Right )
{
    return (Left|Right);
}

Flags FlagsAnd(const Flags& Left, const Flags& Right )
{
    return (Left&Right);
}

void FlagsSet1(Flags& ThisFlag, const Flags& OtherFlag )
{
    ThisFlag.Set(OtherFlag);
}

void FlagsSet2(Flags& ThisFlag, const Flags& OtherFlag, bool Value )
{
    ThisFlag.Set(OtherFlag, Value);
}

template<class TContainerType>
struct Array1DModifier
{
    typedef typename TContainerType::size_type index_type;
    static void Resize( TContainerType& ThisContainer, typename TContainerType::size_type NewSize )
    {
    }
    static void MoveSlice( TContainerType& ThisContainer, index_type Index, index_type From, index_type To )
    {
    }
};

template< class TBinderType, typename TContainerType, typename TVariableType > void VariableIndexingUtility(TBinderType& binder)
    {
        //data value container
        binder.def("__contains__", [](const TContainerType& container, const TVariableType& rV){return container.Has(rV);} );
        binder.def("__setitem__", [](TContainerType& container, const TVariableType& rV, const typename TVariableType::Type rValue){container.SetValue(rV, rValue);} );
        binder.def("__getitem__", [](TContainerType& container, const TVariableType& rV){return container.GetValue(rV);} );
        binder.def("Has", [](const TContainerType& container, const TVariableType& rV){return container.Has(rV);} );
        binder.def("SetValue",  [](TContainerType& container, const TVariableType& rV, const typename TVariableType::Type rValue){container.SetValue(rV, rValue);} );
        binder.def("GetValue", [](TContainerType& container, const TVariableType& rV){return container.GetValue(rV);} );

    }

//     template< typename TVariableType > void RegisterInPythonVariables(pybind11::module& m)
//     {
//         KRATOS_WATCH(KratosComponents<VariableData>::GetComponents().size())
//         for(const auto& item : KratosComponents<VariableData>::GetComponents())
//         {
//             std::cout << "item " << item.first << std::endl;
//             m.attr(item.first.c_str()) = item.second;
//         }
//     }
//
//     void RegisterInPython3DVariablesWithComponents(pybind11::module& m)
//     {
//         for(const auto& item : KratosComponents<Variable<array_1d<double,3>>>::GetComponents())
//         {
//             std::string name = item.first;
//             m.attr(name.c_str()) = item.second;
//
//             std::string xcomponent = name + "_X";
//             std::string ycomponent = name + "_Y";
//             std::string zcomponent = name + "_Z";
//             const auto& xvar = KratosComponents<VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > > >::Get(xcomponent);
//             const auto& yvar = KratosComponents<VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > > >::Get(ycomponent);
//             const auto& zvar = KratosComponents<VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > > >::Get(zcomponent);
//             m.attr(xcomponent.c_str()) = xvar;
//             m.attr(ycomponent.c_str()) = yvar;
//             m.attr(zcomponent.c_str()) = zvar;
//         }
//     }

template <class TVariableType>
TVariableType CreateVariable(const std::string& name)
{
    TVariableType new_var(name);
    // Getting the variable if it exists already
    if(KratosComponents<VariableData>::Has(name)){
        auto& existing_var = KratosComponents<VariableData>::Get(name);
        KRATOS_ERROR_IF(new_var.Key() == existing_var.Key())<<"The variable : "<<name<<" already exists."<<std::endl;
    }

    AddKratosComponent(new_var.Name(), name);
    KratosComponents<VariableData>::Add(new_var.Name(), new_var);
    return new_var;
}


template <class TVariableComponentType, class TVariableComponentAdapterType>
TVariableComponentType CreateVariableComponent(const std::string& name, const std::string& source_name, const int component_index)
{
    TVariableComponentType  new_var(name, source_name, component_index,
                                    TVariableComponentAdapterType(source_name, component_index));

    // Getting the variable if it exists already
    if(KratosComponents<VariableData>::Has(name)){
        auto& existing_var = KratosComponents<TVariableComponentType>::Get(name);
        KRATOS_ERROR_IF(new_var.Key() == existing_var.Key())<<"The variable component : "<<name<<" already exists."<<std::endl;
    }

    AddKratosComponent(new_var.Name(), name);
    KratosComponents<VariableData>::Add(new_var.Name(), new_var);
    return new_var;
}

void  AddContainersToPython(pybind11::module& m)
{
    typedef Variable<array_1d<double, 3> > Array1DVariable3;
    typedef Variable<array_1d<double, 4> > Array1DVariable4;
    typedef Variable<array_1d<double, 6> > Array1DVariable6;
    typedef Variable<array_1d<double, 9> > Array1DVariable9;

    typedef VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > > Array1DComponentVariable;
    typedef VariableComponent<VectorComponentAdaptor<array_1d<double, 4> > > Array1D4ComponentVariable;
    typedef VariableComponent<VectorComponentAdaptor<array_1d<double, 6> > > Array1D6ComponentVariable;
    typedef VariableComponent<VectorComponentAdaptor<array_1d<double, 9> > > Array1D9ComponentVariable;

    //def("TestContainers", TestContainers);

//     BoundedVectorPythonInterface<array_1d<double, 3>, 3>::CreateInterface(m, "Array3" )
//     .def(py::init<vector_expression<array_1d<double, 3> > >() )
//     .def( VectorScalarOperatorPython<array_1d<double, 3>, double, array_1d<double, 3> >() )
//     .def( VectorVectorOperatorPython<array_1d<double, 3>, zero_vector<double>, array_1d<double, 3> >() )
//     .def( VectorVectorOperatorPython<array_1d<double, 3>, unit_vector<double>, array_1d<double, 3> >() )
//     .def( VectorVectorOperatorPython<array_1d<double, 3>, scalar_vector<double>, array_1d<double, 3> >() )
//     .def( VectorVectorOperatorPython<array_1d<double, 3>, mapped_vector<double>, array_1d<double, 3> >() )
//     ;

    py::class_<VariableData>(m, "VariableData" )
    .def("Name", &VariableData::Name, py::return_value_policy::copy)
    .def("__str__", PrintObject<VariableData>)
    ;

    py::class_<Variable<std::string>, VariableData>(m, "StringVariable" )
    .def(py::init<>( [](const std::string& name){return CreateVariable<Variable<std::string>>(name);} ))
    .def("__str__", PrintObject<Variable<std::string>>)
    ;

    py::class_<Variable<bool>, VariableData>(m, "BoolVariable" )
    .def(py::init<>( [](const std::string& name){return CreateVariable<Variable<bool>>(name);} ))
    .def("__str__", PrintObject<Variable<bool>>)
    ;

    py::class_<Variable<int>,VariableData>(m, "IntegerVariable")
    .def(py::init<>( [](const std::string& name){return CreateVariable<Variable<int>>(name);} ))
    .def("__str__", PrintObject<Variable<int>>)
    ;

    py::class_<Variable<DenseVector<int> >,VariableData>(m, "IntegerVectorVariable")
    .def(py::init<>( [](const std::string& name){return CreateVariable<Variable<DenseVector<int>>>(name);} ))
    .def("__str__", PrintObject<Variable<DenseVector<int> >>)
    ;

    py::class_<Variable<double>,VariableData>(m, "DoubleVariable")
    .def(py::init<>( [](const std::string& name){return CreateVariable<Variable<double>>(name);} ))
    .def("__str__", PrintObject<Variable<double>>)
    ;

    py::class_<Variable<Vector >,VariableData>(m, "VectorVariable")
    .def(py::init<>( [](const std::string& name){return CreateVariable<Variable<Vector >>(name);} ))
    .def("__str__", PrintObject<Variable<Vector >>)
    ;

    py::class_<Variable<array_1d<double, 3> >,VariableData>(m, "Array1DVariable3")
    .def(py::init<>( [](const std::string& name)
                 {
                    auto var_x = CreateVariableComponent<VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > >, Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 3> >>
                                    (name+"_X", name, 1);
                    auto var_y = CreateVariableComponent<VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > >, Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 3> >>
                                    (name+"_Y", name, 1);
                    auto var_z = CreateVariableComponent<VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > >, Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 3> >>
                                    (name+"_Z", name, 1);
                    return CreateVariable<Variable<array_1d<double, 3> >>(name);
                 } ))
    .def("__str__", PrintObject<Array1DVariable3>)
    ;

    py::class_<Variable<array_1d<double, 4> >,VariableData>(m, "Array1DVariable4")
    .def(py::init<>( [](const std::string& name){return CreateVariable<Variable<array_1d<double, 4> >>(name);} ))
    .def("__str__", PrintObject<Array1DVariable4>)
    ;

    py::class_<Variable<array_1d<double, 6> >,VariableData>(m, "Array1DVariable6")
    .def(py::init<>( [](const std::string& name){return CreateVariable<Variable<array_1d<double, 6> >>(name);} ))
    .def("__str__", PrintObject<Array1DVariable6>)
    ;

    py::class_<Variable<array_1d<double, 9> >,VariableData>(m, "Array1DVariable9")
    .def(py::init<>( [](const std::string& name){return CreateVariable<Variable<array_1d<double, 9> >>(name);} ))
    .def("__str__", PrintObject<Array1DVariable9>)
    ;

    py::class_<Variable<DenseMatrix<double> >,VariableData>(m, "MatrixVariable")
    .def(py::init<>( [](const std::string& name){return CreateVariable<Variable<DenseMatrix<double> >>(name);} ))
    .def("__str__", PrintObject<Variable<DenseMatrix<double> >>)
    ;

    py::class_<Variable<ConstitutiveLaw::Pointer>,VariableData>(m, "ConstitutuveLawVariable")
    .def("__str__", PrintObject<Variable<ConstitutiveLaw::Pointer>>)
    ;

    py::class_<Variable<ConvectionDiffusionSettings::Pointer > ,VariableData>(m,"ConvectionDiffusionSettingsVariable")
    .def("__str__", PrintObject<Variable<ConvectionDiffusionSettings::Pointer >>)
    ;

    py::class_<Variable<RadiationSettings::Pointer > ,VariableData>(m,"RadiationSettingsVariable")
    .def("__str__", PrintObject<Variable<RadiationSettings::Pointer >>)
    ;
    py::class_<VariableComponent<VectorComponentAdaptor<Vector > >,VariableData>(m, "VectorComponentVariable")
    .def("__str__", PrintObject<VariableComponent<VectorComponentAdaptor<Vector > >>)
    // .def( "GetSourceVariable", &VariableComponent<VectorComponentAdaptor<Vector > >::GetSourceVariable ) // components for vector are not yet fully supported
    ;

    py::class_<VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > >,VariableData>(m, "Array1DComponentVariable")
    .def(py::init<>( [](const std::string& name, const std::string& source_name, const int& component_index)
    {
        return CreateVariableComponent<VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > >, Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 3> >>
            (name, source_name, component_index);
    } ))
    .def("__str__", PrintObject<Array1DComponentVariable>)
    .def( "GetSourceVariable", &Array1DComponentVariable::GetSourceVariable )
    ;

    py::class_<VariableComponent<VectorComponentAdaptor<array_1d<double, 4> > >,VariableData>(m, "Array1D4ComponentVariable")
    .def(py::init<>( [](const std::string& name, const std::string& source_name, const int& component_index)
    {
        return CreateVariableComponent<VariableComponent<VectorComponentAdaptor<array_1d<double, 4> > >, Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 4> >>
            (name, source_name, component_index);
    } ))
    .def("__str__", PrintObject<Array1D4ComponentVariable>)
    .def( "GetSourceVariable", &Array1D4ComponentVariable::GetSourceVariable )
    ;

    py::class_<VariableComponent<VectorComponentAdaptor<array_1d<double, 6> > >,VariableData>(m, "Array1D6ComponentVariable")
    .def(py::init<>( [](const std::string& name, const std::string& source_name, const int& component_index)
    {
        return CreateVariableComponent<VariableComponent<VectorComponentAdaptor<array_1d<double, 6> > >, Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 6> >>
            (name, source_name, component_index);
    } ))
    .def("__str__", PrintObject<Array1D6ComponentVariable>)
    .def( "GetSourceVariable", &Array1D6ComponentVariable::GetSourceVariable )
    ;

    py::class_<VariableComponent<VectorComponentAdaptor<array_1d<double, 9> > >,VariableData>(m, "Array1D9ComponentVariable")
    .def(py::init<>( [](const std::string& name, const std::string& source_name, const int& component_index)
    {
        return CreateVariableComponent<VariableComponent<VectorComponentAdaptor<array_1d<double, 9> > >, Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 9> >>
            (name, source_name, component_index);
    } ))
    .def("__str__", PrintObject<Array1D9ComponentVariable>)
    .def( "GetSourceVariable", &Array1D9ComponentVariable::GetSourceVariable )
    ;

    py::class_<Variable<Quaternion<double> >>(m, "DoubleQuaternionVariable")
    .def("__str__", PrintObject<Variable<Quaternion<double> >>)
    ;

    //***********************************************************************
    //AUTOMATIC REGISTRATION OF VARIABLES_IN_PYTHON
//     RegisterInPythonVariables< Variable<bool> >(m);
//     RegisterInPythonVariables< Variable<int> >(m);
//     RegisterInPythonVariables< Variable<unsigned int> >(m);
//     RegisterInPythonVariables< Variable<double> >(m);
//     RegisterInPythonVariables< Variable<Vector> >(m);
//     RegisterInPythonVariables< Variable<Matrix> >(m);
// //     RegisterInPythonVariables< Variable<ConvectionDiffusionSettings::Pointer> >(m);
// //     RegisterInPythonVariables< Variable<RadiationSettings::Pointer> >(m);
//     RegisterInPythonVariables< VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > > >(m);
//     RegisterInPythonVariables< Variable<Quaternion<double>> >(m);
//     RegisterInPythonVariables< Variable<std::string> >(m);


    //py::class_<AllVariablesDataValueContainer, AllVariablesDataValueContainer::Pointer>( "DataValueContainer" )
    //.def( "__len__", &AllVariablesDataValueContainer::Size )
    //.def( VariableIndexingPython<AllVariablesDataValueContainer, Variable<std::string> >() )
    //.def( VariableIndexingPython<AllVariablesDataValueContainer, Variable<int> >() )
    //.def( VariableIndexingPython<AllVariablesDataValueContainer, Variable<double> >() )
    //.def( VariableIndexingPython<AllVariablesDataValueContainer, Variable<array_1d<double, 3> > >() )
    //.def( VariableIndexingPython<AllVariablesDataValueContainer, Variable<vector<double> > >() )
    //.def( VariableIndexingPython<AllVariablesDataValueContainer, Variable<DenseMatrix<double> > >() )
    //.def( VariableIndexingPython<AllVariablesDataValueContainer, Variable<ConvectionDiffusionSettings::Pointer > >() )
    //.def( VariableIndexingPython<AllVariablesDataValueContainer, Variable<RadiationSettings::Pointer > >() )
    //.def( VariableIndexingPython<AllVariablesDataValueContainer, VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > > >() )
    //.def( self_ns::str( self ) )
    //;

    typedef py::class_<DataValueContainer, DataValueContainer::Pointer> DataValueContainerBinderType;
    DataValueContainerBinderType DataValueBinder(m, "DataValueContainer" );
    DataValueBinder.def( "__len__", &DataValueContainer::Size );
    DataValueBinder.def("__str__", PrintObject<DataValueContainer>);
    VariableIndexingUtility< DataValueContainerBinderType, DataValueContainer, Variable<bool> >(DataValueBinder);
    VariableIndexingUtility< DataValueContainerBinderType, DataValueContainer, Variable<int> >(DataValueBinder);
    VariableIndexingUtility< DataValueContainerBinderType, DataValueContainer, Variable<double> >(DataValueBinder);
    VariableIndexingUtility< DataValueContainerBinderType, DataValueContainer, Array1DVariable3 >(DataValueBinder);
    VariableIndexingUtility< DataValueContainerBinderType, DataValueContainer, Array1DVariable4 >(DataValueBinder);
    VariableIndexingUtility< DataValueContainerBinderType, DataValueContainer, Array1DVariable6 >(DataValueBinder);
    VariableIndexingUtility< DataValueContainerBinderType, DataValueContainer, Array1DVariable9 >(DataValueBinder);
    VariableIndexingUtility< DataValueContainerBinderType, DataValueContainer, Variable<Vector> >(DataValueBinder);
    VariableIndexingUtility< DataValueContainerBinderType, DataValueContainer, Variable<Matrix> >(DataValueBinder);
    VariableIndexingUtility< DataValueContainerBinderType, DataValueContainer, Variable<ConvectionDiffusionSettings::Pointer> >(DataValueBinder);
    VariableIndexingUtility< DataValueContainerBinderType, DataValueContainer, Variable<RadiationSettings::Pointer> >(DataValueBinder);
    VariableIndexingUtility< DataValueContainerBinderType, DataValueContainer, Array1DComponentVariable >(DataValueBinder);
    VariableIndexingUtility< DataValueContainerBinderType, DataValueContainer, Array1D4ComponentVariable >(DataValueBinder);
    VariableIndexingUtility< DataValueContainerBinderType, DataValueContainer, Array1D6ComponentVariable >(DataValueBinder);
    VariableIndexingUtility< DataValueContainerBinderType, DataValueContainer, Array1D9ComponentVariable >(DataValueBinder);
    VariableIndexingUtility< DataValueContainerBinderType, DataValueContainer, Variable<Quaternion<double>> >(DataValueBinder);
    VariableIndexingUtility< DataValueContainerBinderType, DataValueContainer, Variable<std::string> >(DataValueBinder);

    typedef py::class_<VariablesListDataValueContainer, VariablesListDataValueContainer::Pointer> VariableDataValueContainerBinderType;
    VariableDataValueContainerBinderType VariableDataValueBinder(m, "VariablesListDataValueContainer" );
    VariableDataValueBinder.def( "__len__", &VariablesListDataValueContainer::Size );
    VariableDataValueBinder.def("__str__", PrintObject<VariablesListDataValueContainer>);
    VariableIndexingUtility< VariableDataValueContainerBinderType, VariablesListDataValueContainer, Variable<bool> >(VariableDataValueBinder);
    VariableIndexingUtility< VariableDataValueContainerBinderType, VariablesListDataValueContainer, Variable<int> >(VariableDataValueBinder);
    VariableIndexingUtility< VariableDataValueContainerBinderType, VariablesListDataValueContainer, Variable<double> >(VariableDataValueBinder);
    VariableIndexingUtility< VariableDataValueContainerBinderType, VariablesListDataValueContainer, Array1DVariable3 >(VariableDataValueBinder);
    VariableIndexingUtility< VariableDataValueContainerBinderType, VariablesListDataValueContainer, Array1DVariable4 >(VariableDataValueBinder);
    VariableIndexingUtility< VariableDataValueContainerBinderType, VariablesListDataValueContainer, Array1DVariable6 >(VariableDataValueBinder);
    VariableIndexingUtility< VariableDataValueContainerBinderType, VariablesListDataValueContainer, Array1DVariable9 >(VariableDataValueBinder);
    VariableIndexingUtility< VariableDataValueContainerBinderType, VariablesListDataValueContainer, Variable<Vector> >(VariableDataValueBinder);
    VariableIndexingUtility< VariableDataValueContainerBinderType, VariablesListDataValueContainer, Variable<Matrix> >(VariableDataValueBinder);
    VariableIndexingUtility< VariableDataValueContainerBinderType, VariablesListDataValueContainer, Array1DComponentVariable >(VariableDataValueBinder);
    VariableIndexingUtility< VariableDataValueContainerBinderType, VariablesListDataValueContainer, Array1D4ComponentVariable >(VariableDataValueBinder);
    VariableIndexingUtility< VariableDataValueContainerBinderType, VariablesListDataValueContainer, Array1D6ComponentVariable >(VariableDataValueBinder);
    VariableIndexingUtility< VariableDataValueContainerBinderType, VariablesListDataValueContainer, Array1D9ComponentVariable >(VariableDataValueBinder);
    VariableIndexingUtility< VariableDataValueContainerBinderType, VariablesListDataValueContainer, Variable<Quaternion<double>> >(VariableDataValueBinder);
    VariableIndexingUtility< VariableDataValueContainerBinderType, VariablesListDataValueContainer, Variable<std::string> >(VariableDataValueBinder);


    py::class_<Flags, Flags::Pointer>(m,"Flags")
    .def(py::init<>())
    .def(py::init<Flags>())
    .def("Is", &Flags::Is)
    .def("IsNot", &Flags::IsNot)
    .def("Set", FlagsSet1)
    .def("Set", FlagsSet2)
    .def("IsDefined", &Flags::IsDefined)
    .def("IsNotDefined", &Flags::IsNotDefined)
    .def("Reset", &Flags::Reset)
    .def("Flip", &Flags::Flip)
    .def("Clear", &Flags::Clear)
    .def("__or__", FlagsOr)
    .def("__and__", FlagsAnd)
    .def("__str__", PrintObject<Flags>)
    ;

    KRATOS_REGISTER_IN_PYTHON_FLAG(m,STRUCTURE);
    KRATOS_REGISTER_IN_PYTHON_FLAG(m,INTERFACE);
    KRATOS_REGISTER_IN_PYTHON_FLAG(m,FLUID);
    KRATOS_REGISTER_IN_PYTHON_FLAG(m,INLET);
    KRATOS_REGISTER_IN_PYTHON_FLAG(m,OUTLET);
    KRATOS_REGISTER_IN_PYTHON_FLAG(m,VISITED);
    KRATOS_REGISTER_IN_PYTHON_FLAG(m,THERMAL);
    KRATOS_REGISTER_IN_PYTHON_FLAG(m,SELECTED);
    KRATOS_REGISTER_IN_PYTHON_FLAG(m,BOUNDARY);
    KRATOS_REGISTER_IN_PYTHON_FLAG(m,SLIP);
    KRATOS_REGISTER_IN_PYTHON_FLAG(m,CONTACT);
    KRATOS_REGISTER_IN_PYTHON_FLAG(m,TO_SPLIT);
    KRATOS_REGISTER_IN_PYTHON_FLAG(m,TO_ERASE);
    KRATOS_REGISTER_IN_PYTHON_FLAG(m,TO_REFINE);
    KRATOS_REGISTER_IN_PYTHON_FLAG(m,NEW_ENTITY);
    KRATOS_REGISTER_IN_PYTHON_FLAG(m,OLD_ENTITY);
    KRATOS_REGISTER_IN_PYTHON_FLAG(m,ACTIVE);
    KRATOS_REGISTER_IN_PYTHON_FLAG(m,MODIFIED);
    KRATOS_REGISTER_IN_PYTHON_FLAG(m,RIGID);
    KRATOS_REGISTER_IN_PYTHON_FLAG(m,SOLID);
    KRATOS_REGISTER_IN_PYTHON_FLAG(m,MPI_BOUNDARY);
    KRATOS_REGISTER_IN_PYTHON_FLAG(m,INTERACTION);
    KRATOS_REGISTER_IN_PYTHON_FLAG(m,ISOLATED);
    KRATOS_REGISTER_IN_PYTHON_FLAG(m,MASTER);
    KRATOS_REGISTER_IN_PYTHON_FLAG(m,SLAVE);
    KRATOS_REGISTER_IN_PYTHON_FLAG(m,INSIDE);
    KRATOS_REGISTER_IN_PYTHON_FLAG(m,FREE_SURFACE);
    KRATOS_REGISTER_IN_PYTHON_FLAG(m,BLOCKED);
    KRATOS_REGISTER_IN_PYTHON_FLAG(m,MARKER);
    KRATOS_REGISTER_IN_PYTHON_FLAG(m,PERIODIC);


//     AddDeprecatedVariablesToPython();
//     AddC2CVariablesToPython();
//     AddDEMVariablesToPython(); //TODO: move this to the DEM application
//     AddCFDVariablesToPython(); ///@TODO: move variables to CFD application
//     AddALEVariablesToPython(); ///@TODO: move variables to ALE application
//     AddFSIVariablesToPython(); ///@TODO: move variables to FSI application
//     AddMappingVariablesToPython(); ///@TODO: move variables to Mapping application
//     AddMATVariablesToPython(); ///@TODO: move variables to CL application
//     AddLegacyStructuralAppVarsToPython();

    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SPACE_DIMENSION )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, DOMAIN_SIZE )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, IS_RESTARTED )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, COMPUTE_LUMPED_MASS_MATRIX )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, COMPUTE_DYNAMIC_TANGENT )

    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, THERMAL_EXPANSION_COEFFICIENT )

    // These should be moved to applications
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, POWER_LAW_N )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, POWER_LAW_K )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, EQ_STRAIN_RATE )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, YIELD_STRESS )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, MU )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, TAU )


    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, CONSTRAINT_LABELS )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, LOAD_LABELS )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, MARKER_LABELS )

    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, CONSTRAINT_MESHES )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, LOAD_MESHES )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, MARKER_MESHES )

    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, ELEMENTAL_DISTANCES )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, NL_ITERATION_NUMBER )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, FRACTIONAL_STEP )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, STEP )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, TIME )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, START_TIME )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, END_TIME )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, DELTA_TIME )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, PREVIOUS_DELTA_TIME )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, INTERVAL_END_TIME )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, PRINTED_STEP )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, PRINTED_RESTART_STEP )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, RESIDUAL_NORM )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, CONVERGENCE_RATIO )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, TEMPERATURE )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, TEMPERATURE_OLD_IT )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, PRESSURE )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, VISCOSITY_AIR )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, VISCOSITY_WATER )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, ERROR_RATIO )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, TIME_STEPS )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SCALAR_LAGRANGE_MULTIPLIER )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, VECTOR_LAGRANGE_MULTIPLIER )

    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, ANGULAR_ACCELERATION )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, VELOCITY_LAPLACIAN )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, VELOCITY_LAPLACIAN_RATE )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, VELOCITY_COMPONENT_GRADIENT )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, VELOCITY_X_GRADIENT )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, VELOCITY_Y_GRADIENT )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, VELOCITY_Z_GRADIENT )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, ANGULAR_VELOCITY )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, ACCELERATION )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, VELOCITY )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, DISPLACEMENT )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, ROTATION )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, DELTA_ROTATION )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, REACTION_MOMENT )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, REACTION )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, BODY_FORCE )

    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, FORCE_RESIDUAL )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, MOMENT_RESIDUAL )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, INTERNAL_FORCE )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, EXTERNAL_FORCE )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, CONTACT_FORCE )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, CONTACT_NORMAL )

    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, EXTERNAL_FORCES_VECTOR )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, INTERNAL_FORCES_VECTOR )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, CONTACT_FORCES_VECTOR )

    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, LINEAR_MOMENTUM )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, ANGULAR_MOMENTUM )

    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, VOLUME_ACCELERATION )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, SEEPAGE_DRAG )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, NORMAL )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, TANGENT_XI )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, TANGENT_ETA )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, FORCE )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, TORQUE )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, MOMENT )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, FORCE_CM )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, MOMENTUM_CM )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, MOMENTUM )

    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, MASS )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, RHS )

    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, WATER_PRESSURE )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, REACTION_WATER_PRESSURE )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, WATER_PRESSURE_ACCELERATION )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, AIR_PRESSURE )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, REACTION_AIR_PRESSURE )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, RHS_WATER )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, RHS_AIR )

    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, WEIGHT_FATHER_NODES )

    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, INTERNAL_ENERGY )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, STRAIN_ENERGY )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, EXTERNAL_ENERGY )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, TOTAL_ENERGY )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, KINETIC_ENERGY )

    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, LOCAL_INERTIA_TENSOR )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, LOCAL_AXES_MATRIX )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, LOCAL_CONSTITUTIVE_MATRIX )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, CONSTITUTIVE_MATRIX )

    //for structural application TO BE REMOVED
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, CONSTITUTIVE_LAW )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, INTERNAL_VARIABLES )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, MATERIAL_PARAMETERS )

    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, NEGATIVE_FACE_PRESSURE )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, POSITIVE_FACE_PRESSURE )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, POROSITY )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, DIAMETER )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, LIN_DARCY_COEF )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, NONLIN_DARCY_COEF )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, DRAG_FORCE )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, STRUCTURE_VELOCITY )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, K0 )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, NODAL_VOLUME )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, STATIONARY )

    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, FLAG_VARIABLE )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, DISTANCE )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, DISTANCE_GRADIENT )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, INERTIA )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, PERIODIC_PAIR_INDEX )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, PARTITION_INDEX )

    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, LAGRANGE_DISPLACEMENT )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, LAGRANGE_AIR_PRESSURE )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, LAGRANGE_WATER_PRESSURE )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, LAGRANGE_TEMPERATURE )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, INTERNAL_FRICTION_ANGLE )

    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, LAGRANGE_DISPLACEMENT )



    // for MultiScale application
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, INITIAL_STRAIN )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, COEFFICIENT_THERMAL_EXPANSION )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, CHARACTERISTIC_LENGTH_MULTIPLIER )

    //for Incompressible Fluid application
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m,FRACT_VEL)

    //for ALE application
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, DETERMINANT )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, ELEMENTSHAPE )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, MESH_VELOCITY )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, AUX_MESH_VAR )

    //for Adjoint
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, SHAPE_SENSITIVITY )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, NORMAL_SENSITIVITY )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, NUMBER_OF_NEIGHBOUR_ELEMENTS )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, UPDATE_SENSITIVITIES )

    //for electric application

    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, PARTITION_MASK )

    // For MeshingApplication
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, NODAL_ERROR )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, NODAL_ERROR_COMPONENTS )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, ELEMENT_ERROR )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, ELEMENT_H )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, RECOVERED_STRESS )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, ERROR_INTEGRATION_POINT )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, CONTACT_PRESSURE )

    // For explicit time integration
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, RESIDUAL_VECTOR )

    //for PFEM application TO BE REMOVED
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, NODAL_AREA )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, NODAL_H )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, NORMAL_TO_WALL )
    //KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,NEIGHBOUR_NODES)
    //KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,NEIGHBOUR_ELEMENTS)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, FRICTION_COEFFICIENT )

    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, BULK_MODULUS )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SATURATION )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, GRAVITY )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, FACE_LOAD )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, DENSITY_WATER )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, DENSITY_AIR )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, AIR_ENTRY_VALUE )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, FIRST_SATURATION_PARAM )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SECOND_SATURATION_PARAM )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, PERMEABILITY_WATER )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, PERMEABILITY_AIR )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, BULK_AIR )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SCALE )

    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, TEMP_CONV_PROJ )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, CONVECTION_COEFFICIENT)

    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, INSITU_STRESS )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, STRESSES )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, STRAIN )

    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, NODAL_MASS )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, AUX_INDEX )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, EXTERNAL_PRESSURE )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, BDF_COEFFICIENTS )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, ROTATION_CENTER )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, VELOCITY_PERIOD )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, ANGULAR_VELOCITY_PERIOD )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, IDENTIFIER )


    //for xfem application
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, CRACK_OPENING )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, CRACK_TRANSLATION )

    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,ARRHENIUS)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,ARRHENIUSAUX)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,ARRHENIUSAUX_)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,PRESSUREAUX)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,NODAL_MAUX)
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m,NODAL_VAUX)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,NODAL_PAUX)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,FACE_HEAT_FLUX)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,HEAT_FLUX)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,REACTION_FLUX)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,TC)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,CONDUCTIVITY)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,SPECIFIC_HEAT)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,MATERIAL_VARIABLE)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,FUEL)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,YO)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,YF)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,YI)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,Y1)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,Y2)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,YP)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,EMISSIVITY)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,ENTHALPY)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,MIXTURE_FRACTION)

    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,YCH4)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,YO2)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,YCO2)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,YH2O)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,YN2)

    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,INCIDENT_RADIATION_FUNCTION)

    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,ABSORPTION_COEFFICIENT)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,STEFAN_BOLTZMANN_CONSTANT)
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, DIRECTION )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,NODAL_SWITCH)
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m,Y)

    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, LOCAL_AXIS_1 )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, LOCAL_AXIS_2 )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, LOCAL_AXIS_3 )

    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SWITCH_TEMPERATURE )

    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m,EMBEDDED_VELOCITY)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, REFINEMENT_LEVEL )

    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, AIR_SOUND_VELOCITY )

    // for Vulcan application
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, LATENT_HEAT )

    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, ENRICHED_PRESSURES )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, INITIAL_PENALTY )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, DP_EPSILON )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, DP_ALPHA1 )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, DP_K )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, LAMBDA )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, MIN_DT )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, MAX_DT )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, WET_VOLUME)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, CUTTED_AREA)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, NET_INPUT_MATERIAL)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, MIU )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SCALE_FACTOR )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, NORMAL_CONTACT_STRESS )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, TANGENTIAL_CONTACT_STRESS )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, STABILIZATION_FACTOR )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, NEWMARK_BETA )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, NEWMARK_GAMMA )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, BOSSAK_ALPHA )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, EQUILIBRIUM_POINT )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, AIR_SOUND_VELOCITY )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, WATER_SOUND_VELOCITY )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, ACTIVATION_LEVEL )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, FIRST_TIME_STEP )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, QUASI_STATIC_ANALYSIS )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, FRACTIONAL_STEP )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, LOAD_RESTART )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, RIGID_BODY_ID )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, NORMAL_STRESS )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, TANGENTIAL_STRESS )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, PENALTY )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, AMBIENT_TEMPERATURE )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, VEL_ART_VISC )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, PR_ART_VISC )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SOUND_VELOCITY )

    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SEARCH_RADIUS )

    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, INTEGRATION_WEIGHT )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, INTEGRATION_COORDINATES )


    py::class_< ConvectionDiffusionSettings, ConvectionDiffusionSettings::Pointer >	(m,"ConvectionDiffusionSettings")
    .def(py::init<	>() )
    .def("SetDensityVariable",&ConvectionDiffusionSettings::SetDensityVariable)
    .def("SetDiffusionVariable",&ConvectionDiffusionSettings::SetDiffusionVariable)
    .def("SetUnknownVariable",&ConvectionDiffusionSettings::SetUnknownVariable)
    .def("SetVolumeSourceVariable",&ConvectionDiffusionSettings::SetVolumeSourceVariable)
    .def("SetSurfaceSourceVariable",&ConvectionDiffusionSettings::SetSurfaceSourceVariable)
    .def("SetProjectionVariable",&ConvectionDiffusionSettings::SetProjectionVariable)
    .def("SetMeshVelocityVariable",&ConvectionDiffusionSettings::SetMeshVelocityVariable)
    .def("SetConvectionVariable",&ConvectionDiffusionSettings::SetConvectionVariable)
    .def("SetTransferCoefficientVariable",&ConvectionDiffusionSettings::SetTransferCoefficientVariable)
    .def("SetSpecificHeatVariable",&ConvectionDiffusionSettings::SetSpecificHeatVariable)
    .def("SetVelocityVariable",&ConvectionDiffusionSettings::SetVelocityVariable)
    .def("SetReactionVariable",&ConvectionDiffusionSettings::SetReactionVariable)

    .def("GetDensityVariable",&ConvectionDiffusionSettings::GetDensityVariable, py::return_value_policy::reference_internal )
    .def("GetDiffusionVariable",&ConvectionDiffusionSettings::GetDiffusionVariable, py::return_value_policy::reference_internal )
    .def("GetUnknownVariable",&ConvectionDiffusionSettings::GetUnknownVariable, py::return_value_policy::reference_internal )
    .def("GetVolumeSourceVariable",&ConvectionDiffusionSettings::GetVolumeSourceVariable, py::return_value_policy::reference_internal )
    .def("GetSurfaceSourceVariable",&ConvectionDiffusionSettings::GetSurfaceSourceVariable, py::return_value_policy::reference_internal )
    .def("GetProjectionVariable",&ConvectionDiffusionSettings::GetProjectionVariable, py::return_value_policy::reference_internal )
    .def("GetMeshVelocityVariable",&ConvectionDiffusionSettings::GetMeshVelocityVariable, py::return_value_policy::reference_internal )
    .def("GetConvectionVariable",&ConvectionDiffusionSettings::GetConvectionVariable, py::return_value_policy::reference_internal )
    .def("GetTransferCoefficientVariable",&ConvectionDiffusionSettings::GetTransferCoefficientVariable, py::return_value_policy::reference_internal)
    .def("GetSpecificHeatVariable",&ConvectionDiffusionSettings::GetSpecificHeatVariable, py::return_value_policy::reference_internal )
    .def("GetVelocityVariable",&ConvectionDiffusionSettings::GetVelocityVariable, py::return_value_policy::reference_internal )
    .def("GetReactionVariable",&ConvectionDiffusionSettings::GetReactionVariable, py::return_value_policy::reference_internal )

    .def("IsDefinedDensityVariable",&ConvectionDiffusionSettings::IsDefinedDensityVariable)
    .def("IsDefinedDiffusionVariable",&ConvectionDiffusionSettings::IsDefinedDiffusionVariable)
    .def("IsDefinedUnknownVariable",&ConvectionDiffusionSettings::IsDefinedUnknownVariable)
    .def("IsDefinedVolumeSourceVariable",&ConvectionDiffusionSettings::IsDefinedVolumeSourceVariable)
    .def("IsDefinedSurfaceSourceVariable",&ConvectionDiffusionSettings::IsDefinedSurfaceSourceVariable)
    .def("IsDefinedProjectionVariable",&ConvectionDiffusionSettings::IsDefinedProjectionVariable)
    .def("IsDefinedMeshVelocityVariable",&ConvectionDiffusionSettings::IsDefinedMeshVelocityVariable)
    .def("IsDefinedConvectionVariable",&ConvectionDiffusionSettings::IsDefinedConvectionVariable)
    .def("IsDefinedSpecificHeatVariable",&ConvectionDiffusionSettings::IsDefinedSpecificHeatVariable)
    .def("IsDefinedVelocityVariable",&ConvectionDiffusionSettings::IsDefinedVelocityVariable)
    .def("IsDefinedTransferCoefficientVariable",&ConvectionDiffusionSettings::IsDefinedTransferCoefficientVariable)
    .def("IsDefinedReactionVariable",&ConvectionDiffusionSettings::IsDefinedReactionVariable)
    ;

    py::class_< RadiationSettings, RadiationSettings::Pointer>	(m,"RadiationSettings")
    .def(py::init<	>() )
    .def("SetDensityVariable",&RadiationSettings::SetDensityVariable)
    .def("SetDiffusionVariable",&RadiationSettings::SetDiffusionVariable)
    .def("SetUnknownVariable",&RadiationSettings::SetUnknownVariable)
    .def("SetVolumeSourceVariable",&RadiationSettings::SetVolumeSourceVariable)
    .def("SetSurfaceSourceVariable",&RadiationSettings::SetSurfaceSourceVariable)
    .def("SetProjectionVariable",&RadiationSettings::SetProjectionVariable)
    .def("SetMeshVelocityVariable",&RadiationSettings::SetMeshVelocityVariable)
    .def("GetDensityVariable",&RadiationSettings::GetDensityVariable, py::return_value_policy::reference_internal )
    .def("GetDiffusionVariable",&RadiationSettings::GetDiffusionVariable, py::return_value_policy::reference_internal )
    .def("GetUnknownVariable",&RadiationSettings::GetUnknownVariable, py::return_value_policy::reference_internal )
    .def("GetVolumeSourceVariable",&RadiationSettings::GetVolumeSourceVariable, py::return_value_policy::reference_internal )
    .def("GetSurfaceSourceVariable",&RadiationSettings::GetSurfaceSourceVariable, py::return_value_policy::reference_internal )
    //.def("GetSurfaceSourceVariable",&RadiationSettings::GetSurfaceSourceVariable, py::return_value_policy::reference_internal )
    .def("GetProjectionVariable",&RadiationSettings::GetProjectionVariable, py::return_value_policy::reference_internal )
    .def("GetMeshVelocityVariable",&RadiationSettings::GetMeshVelocityVariable, py::return_value_policy::reference_internal )
    ;
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,CONVECTION_DIFFUSION_SETTINGS)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,RADIATION_SETTINGS)

    // KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,CONVECTION_DIFFUSION_SETTINGS)

}
} // namespace Python.
} // Namespace Kratos
