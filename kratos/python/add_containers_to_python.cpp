
//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//

// System includes

// External includes

// Project includes
#include "containers/data_value_container.h"
#include "containers/variables_list_data_value_container.h"
#include "containers/flags.h"
#include "containers/variable.h"
#include "includes/define_python.h"
#include "includes/kratos_flags.h"
#include "includes/constitutive_law.h"
#include "includes/convection_diffusion_settings.h"
#include "includes/radiation_settings.h"
#include "utilities/quaternion.h"

namespace Kratos::Python
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
    // Data container
    binder.def("__contains__", [](const TContainerType& container, const TVariableType& rV){return container.Has(rV);} );
    binder.def("__setitem__", [](TContainerType& container, const TVariableType& rV, const typename TVariableType::Type rValue){container.SetValue(rV, rValue);} );
    binder.def("__getitem__", [](TContainerType& container, const TVariableType& rV){return container.GetValue(rV);} );
    binder.def("Has", [](const TContainerType& container, const TVariableType& rV){return container.Has(rV);} );
    binder.def("SetValue",  [](TContainerType& container, const TVariableType& rV, const typename TVariableType::Type rValue){container.SetValue(rV, rValue);} );
    binder.def("GetValue", [](TContainerType& container, const TVariableType& rV){return container.GetValue(rV);} );
    binder.def("Clear", [](TContainerType& container){container.Clear();} );
}

template< class TBinderType, typename TContainerType, typename TVariableType > void DataValueContainerIndexingUtility(TBinderType& binder)
{
    // Data value container
    VariableIndexingUtility<TBinderType, TContainerType, TVariableType>(binder);
    binder.def("Erase", [](TContainerType& container, const TVariableType& rV){return container.Erase(rV);} );
}

void  AddContainersToPython(pybind11::module& m)
{
    using Array1DVariable3 = Variable<array_1d<double, 3> >;
    using Array1DVariable4 = Variable<array_1d<double, 4> >;
    using Array1DVariable6 = Variable<array_1d<double, 6> >;
    using Array1DVariable9 = Variable<array_1d<double, 9> >;

    py::class_<VariableData>(m, "VariableData" )
    .def("Name", &VariableData::Name, py::return_value_policy::copy)
    .def("Key", &VariableData::Key)
    .def("GetSourceVariable", &VariableData::GetSourceVariable)
    .def("GetComponentIndex", &VariableData::GetComponentIndex)
    .def("IsComponent", &VariableData::IsComponent)
    .def("__str__", PrintObject<VariableData>)
    ;

    py::class_<Variable<std::string>, VariableData>(m, "StringVariable" )
    .def("__str__", PrintObject<Variable<std::string>>)
    ;

    py::class_<Variable<bool>, VariableData>(m, "BoolVariable" )
    .def("__str__", PrintObject<Variable<bool>>)
    ;

    py::class_<Variable<int>,VariableData>(m, "IntegerVariable")
    .def("__str__", PrintObject<Variable<int>>)
    ;

    py::class_<Variable<DenseVector<int> >,VariableData>(m, "IntegerVectorVariable")
    .def("__str__", PrintObject<Variable<DenseVector<int> >>)
    ;

    py::class_<Variable<double>,VariableData>(m, "DoubleVariable")
    .def("__str__", PrintObject<Variable<double>>)
    ;

    py::class_<Variable<Vector >,VariableData>(m, "VectorVariable")
    .def("__str__", PrintObject<Variable<Vector >>)
    ;

    py::class_<Array1DVariable3,VariableData>(m, "Array1DVariable3")
    .def("__str__", PrintObject<Array1DVariable3>)
    ;

    py::class_<Array1DVariable4,VariableData>(m, "Array1DVariable4")
    .def("__str__", PrintObject<Array1DVariable4>)
    ;

    py::class_<Array1DVariable6,VariableData>(m, "Array1DVariable6")
    .def("__str__", PrintObject<Array1DVariable6>)
    ;

    py::class_<Array1DVariable9,VariableData>(m, "Array1DVariable9")
    .def("__str__", PrintObject<Array1DVariable9>)
    ;

    py::class_<Variable<DenseMatrix<double> >,VariableData>(m, "MatrixVariable")
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

    py::class_<Variable<Quaternion<double> >>(m, "DoubleQuaternionVariable")
    .def("__str__", PrintObject<Variable<Quaternion<double> >>)
    ;

    typedef py::class_<DataValueContainer, DataValueContainer::Pointer> DataValueContainerBinderType;
    DataValueContainerBinderType DataValueBinder(m, "DataValueContainer" );
    DataValueBinder.def( "__len__", &DataValueContainer::Size );
    DataValueBinder.def("__str__", PrintObject<DataValueContainer>);
    DataValueContainerIndexingUtility< DataValueContainerBinderType, DataValueContainer, Variable<bool> >(DataValueBinder);
    DataValueContainerIndexingUtility< DataValueContainerBinderType, DataValueContainer, Variable<int> >(DataValueBinder);
    DataValueContainerIndexingUtility< DataValueContainerBinderType, DataValueContainer, Variable<double> >(DataValueBinder);
    DataValueContainerIndexingUtility< DataValueContainerBinderType, DataValueContainer, Array1DVariable3 >(DataValueBinder);
    DataValueContainerIndexingUtility< DataValueContainerBinderType, DataValueContainer, Array1DVariable4 >(DataValueBinder);
    DataValueContainerIndexingUtility< DataValueContainerBinderType, DataValueContainer, Array1DVariable6 >(DataValueBinder);
    DataValueContainerIndexingUtility< DataValueContainerBinderType, DataValueContainer, Array1DVariable9 >(DataValueBinder);
    DataValueContainerIndexingUtility< DataValueContainerBinderType, DataValueContainer, Variable<Vector> >(DataValueBinder);
    DataValueContainerIndexingUtility< DataValueContainerBinderType, DataValueContainer, Variable<Matrix> >(DataValueBinder);
    DataValueContainerIndexingUtility< DataValueContainerBinderType, DataValueContainer, Variable<ConvectionDiffusionSettings::Pointer> >(DataValueBinder);
    DataValueContainerIndexingUtility< DataValueContainerBinderType, DataValueContainer, Variable<RadiationSettings::Pointer> >(DataValueBinder);
    DataValueContainerIndexingUtility< DataValueContainerBinderType, DataValueContainer, Variable<Quaternion<double>> >(DataValueBinder);
    DataValueContainerIndexingUtility< DataValueContainerBinderType, DataValueContainer, Variable<std::string> >(DataValueBinder);

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
    .def("AsFalse", &Flags::AsFalse)
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
    KRATOS_REGISTER_IN_PYTHON_FLAG(m,WALL);

    // Note: using internal macro for these two because they do not have a NOT_ version
    KRATOS_REGISTER_IN_PYTHON_FLAG_IMPLEMENTATION(m,ALL_DEFINED);
    KRATOS_REGISTER_IN_PYTHON_FLAG_IMPLEMENTATION(m,ALL_TRUE);

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
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, ELEMENTAL_EDGE_DISTANCES )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, ELEMENTAL_EDGE_DISTANCES_EXTRAPOLATED )
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
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, RUNGE_KUTTA_STEP )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, RESIDUAL_NORM )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, CONVERGENCE_RATIO )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, BUILD_SCALE_FACTOR )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, CONSTRAINT_SCALE_FACTOR )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, AUXILIAR_CONSTRAINT_SCALE_FACTOR )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, TEMPERATURE )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, TEMPERATURE_OLD_IT )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, PRESSURE )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, VISCOSITY_AIR )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, VISCOSITY_WATER )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, ERROR_RATIO )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, TIME_STEPS )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, PENALTY_COEFFICIENT )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SCALAR_LAGRANGE_MULTIPLIER )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, TIME_INTEGRATION_THETA )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SHAPE_FUNCTIONS_VECTOR)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SHAPE_FUNCTIONS_GRADIENT_MATRIX)
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
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, TEMPERATURE_GRADIENT )

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
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, LOCAL_TANGENT )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, LOCAL_TANGENT_MATRIX )
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

    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, VOLUMETRIC_STRAIN )

    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, LOCAL_INERTIA_TENSOR )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, LOCAL_AXES_MATRIX )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, LOCAL_CONSTITUTIVE_MATRIX )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, CONSTITUTIVE_MATRIX )

    // for geometrical application
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, CHARACTERISTIC_GEOMETRY_LENGTH)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, DETERMINANTS_OF_JACOBIAN_PARENT)

    //for structural application TO BE REMOVED
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, NUMBER_OF_CYCLES )
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
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, AUX_DISTANCE )
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
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, NORMAL_SHAPE_DERIVATIVE )

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
    KRATOS_REGISTER_IN_PYTHON_SYMMETRIC_3D_TENSOR_VARIABLE_WITH_COMPONENTS(m, NODAL_INERTIA_TENSOR )
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
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, LIMITER_COEFFICIENT )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SOUND_VELOCITY )

    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SEARCH_RADIUS )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, ORIENTATION )

    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, INTEGRATION_WEIGHT )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, INTEGRATION_COORDINATES )

    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, PARAMETER_2D_COORDINATES)

    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, VARIATIONAL_REDISTANCE_COEFFICIENT_FIRST)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, VARIATIONAL_REDISTANCE_COEFFICIENT_SECOND)


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
    .def("SetGradientVariable",&ConvectionDiffusionSettings::SetGradientVariable)
    .def("SetTransferCoefficientVariable",&ConvectionDiffusionSettings::SetTransferCoefficientVariable)
    .def("SetSpecificHeatVariable",&ConvectionDiffusionSettings::SetSpecificHeatVariable)
    .def("SetVelocityVariable",&ConvectionDiffusionSettings::SetVelocityVariable)
    .def("SetReactionVariable",&ConvectionDiffusionSettings::SetReactionVariable)
    .def("SetReactionGradientVariable",&ConvectionDiffusionSettings::SetReactionGradientVariable)

    .def("GetDensityVariable",&ConvectionDiffusionSettings::GetDensityVariable, py::return_value_policy::reference_internal )
    .def("GetDiffusionVariable",&ConvectionDiffusionSettings::GetDiffusionVariable, py::return_value_policy::reference_internal )
    .def("GetUnknownVariable",&ConvectionDiffusionSettings::GetUnknownVariable, py::return_value_policy::reference_internal )
    .def("GetVolumeSourceVariable",&ConvectionDiffusionSettings::GetVolumeSourceVariable, py::return_value_policy::reference_internal )
    .def("GetSurfaceSourceVariable",&ConvectionDiffusionSettings::GetSurfaceSourceVariable, py::return_value_policy::reference_internal )
    .def("GetProjectionVariable",&ConvectionDiffusionSettings::GetProjectionVariable, py::return_value_policy::reference_internal )
    .def("GetMeshVelocityVariable",&ConvectionDiffusionSettings::GetMeshVelocityVariable, py::return_value_policy::reference_internal )
    .def("GetConvectionVariable",&ConvectionDiffusionSettings::GetConvectionVariable, py::return_value_policy::reference_internal )
    .def("GetGradientVariable",&ConvectionDiffusionSettings::GetGradientVariable, py::return_value_policy::reference_internal )
    .def("GetTransferCoefficientVariable",&ConvectionDiffusionSettings::GetTransferCoefficientVariable, py::return_value_policy::reference_internal)
    .def("GetSpecificHeatVariable",&ConvectionDiffusionSettings::GetSpecificHeatVariable, py::return_value_policy::reference_internal )
    .def("GetVelocityVariable",&ConvectionDiffusionSettings::GetVelocityVariable, py::return_value_policy::reference_internal )
    .def("GetReactionVariable",&ConvectionDiffusionSettings::GetReactionVariable, py::return_value_policy::reference_internal )
    .def("GetReactionGradientVariable",&ConvectionDiffusionSettings::GetReactionGradientVariable, py::return_value_policy::reference_internal )

    .def("IsDefinedDensityVariable",&ConvectionDiffusionSettings::IsDefinedDensityVariable)
    .def("IsDefinedDiffusionVariable",&ConvectionDiffusionSettings::IsDefinedDiffusionVariable)
    .def("IsDefinedUnknownVariable",&ConvectionDiffusionSettings::IsDefinedUnknownVariable)
    .def("IsDefinedVolumeSourceVariable",&ConvectionDiffusionSettings::IsDefinedVolumeSourceVariable)
    .def("IsDefinedSurfaceSourceVariable",&ConvectionDiffusionSettings::IsDefinedSurfaceSourceVariable)
    .def("IsDefinedProjectionVariable",&ConvectionDiffusionSettings::IsDefinedProjectionVariable)
    .def("IsDefinedMeshVelocityVariable",&ConvectionDiffusionSettings::IsDefinedMeshVelocityVariable)
    .def("IsDefinedConvectionVariable",&ConvectionDiffusionSettings::IsDefinedConvectionVariable)
    .def("IsDefinedGradientVariable",&ConvectionDiffusionSettings::IsDefinedGradientVariable)
    .def("IsDefinedSpecificHeatVariable",&ConvectionDiffusionSettings::IsDefinedSpecificHeatVariable)
    .def("IsDefinedVelocityVariable",&ConvectionDiffusionSettings::IsDefinedVelocityVariable)
    .def("IsDefinedTransferCoefficientVariable",&ConvectionDiffusionSettings::IsDefinedTransferCoefficientVariable)
    .def("IsDefinedReactionVariable",&ConvectionDiffusionSettings::IsDefinedReactionVariable)
    .def("IsDefinedReactionGradientVariable",&ConvectionDiffusionSettings::IsDefinedReactionGradientVariable)
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
}
} // namespace Kratos::Python.
