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
#include <boost/python.hpp>


// Project includes
#include "includes/define.h"
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
#include "python/variable_indexing_python.h"
#include "python/vector_python_interface.h"
#include "python/vector_scalar_operator_python.h"
#include "python/vector_vector_operator_python.h"
#include "python/bounded_vector_python_interface.h"
#include "python/add_deprecated_variables_to_python.h"
#include "python/add_c2c_variables_to_python.h" //TODO: to be removed eventually
#include "python/add_cfd_variables_to_python.h" //TODO: to be removed eventually
#include "python/add_ale_variables_to_python.h" //TODO: to be removed eventually
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
using namespace boost::python;

Flags FlagsOr(const Flags& Left, const Flags& Right )
{
    return (Left|Right);
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

void  AddContainersToPython()
{
    //def("TestContainers", TestContainers);

    BoundedVectorPythonInterface<array_1d<double, 3>, 3>::CreateInterface( "Array3" )
    .def( init<vector_expression<array_1d<double, 3> > >() )
    .def( VectorScalarOperatorPython<array_1d<double, 3>, double, array_1d<double, 3> >() )
    .def( VectorVectorOperatorPython<array_1d<double, 3>, zero_vector<double>, array_1d<double, 3> >() )
    .def( VectorVectorOperatorPython<array_1d<double, 3>, unit_vector<double>, array_1d<double, 3> >() )
    .def( VectorVectorOperatorPython<array_1d<double, 3>, scalar_vector<double>, array_1d<double, 3> >() )
    .def( VectorVectorOperatorPython<array_1d<double, 3>, mapped_vector<double>, array_1d<double, 3> >() )
    ;

    class_<VariableData>( "VariableData", no_init )
    .def("Name", &VariableData::Name, return_value_policy<copy_const_reference>())
    .def( self_ns::str( self ) )
    ;

    class_<Variable<std::string>, bases<VariableData>, boost::noncopyable >( "StringVariable", no_init )
    .def( self_ns::str( self ) )
    ;

    class_<Variable<bool>, bases<VariableData>, boost::noncopyable >( "BoolVariable", no_init )
    .def( self_ns::str( self ) )
    ;

    class_<Variable<int>, bases<VariableData>, boost::noncopyable >( "IntegerVariable", no_init )
    .def( self_ns::str( self ) )
    ;

    class_<Variable<vector<int> >, bases<VariableData>, boost::noncopyable >( "IntegerVectorVariable", no_init )
    .def( self_ns::str( self ) )
    ;

    class_<Variable<double>, bases<VariableData>, boost::noncopyable >( "DoubleVariable", no_init )
    .def( self_ns::str( self ) )
    ;

    class_<Variable<vector<double> >, bases<VariableData>, boost::noncopyable >( "VectorVariable", no_init )
    .def( self_ns::str( self ) )
    ;

    class_<Variable<array_1d<double, 3> >, bases<VariableData>, boost::noncopyable >( "Array1DVariable3", no_init )
    .def( self_ns::str( self ) )
    ;

    class_<Variable<array_1d<double, 6> >, bases<VariableData>, boost::noncopyable >( "Array1DVariable6", no_init )
    .def( self_ns::str( self ) )
    ;

    class_<Variable<matrix<double> >, bases<VariableData>, boost::noncopyable >( "MatrixVariable", no_init )
    .def( self_ns::str( self ) )
    ;

    class_<Variable<ConstitutiveLaw::Pointer>, bases<VariableData>, boost::noncopyable >( "ConstitutuveLawVariable", no_init )
    .def( self_ns::str( self ) )
    ;

    class_<Variable<ConvectionDiffusionSettings::Pointer > , bases<VariableData>, boost::noncopyable >("ConvectionDiffusionSettingsVariable", no_init)
    .def( self_ns::str( self ) )
    ;

    class_<Variable<RadiationSettings::Pointer > , bases<VariableData>, boost::noncopyable >("RadiationSettingsVariable", no_init)
    .def( self_ns::str( self ) )
    ;
    class_<VariableComponent<VectorComponentAdaptor<vector<double> > >, bases<VariableData>, boost::noncopyable >( "VectorComponentVariable", no_init )
    .def( self_ns::str( self ) )
    ;

    class_<VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > >, bases<VariableData>, boost::noncopyable >( "Array1DComponentVariable", no_init )
    .def( self_ns::str( self ) )
    ;

    class_<Variable<Quaternion<double> >, boost::noncopyable >( "DoubleQuaternionVariable", no_init )
    .def( self_ns::str( self ) )
    ;


    //class_<AllVariablesDataValueContainer, AllVariablesDataValueContainer::Pointer>( "DataValueContainer" )
    //.def( "__len__", &AllVariablesDataValueContainer::Size )
    //.def( VariableIndexingPython<AllVariablesDataValueContainer, Variable<std::string> >() )
    //.def( VariableIndexingPython<AllVariablesDataValueContainer, Variable<int> >() )
    //.def( VariableIndexingPython<AllVariablesDataValueContainer, Variable<double> >() )
    //.def( VariableIndexingPython<AllVariablesDataValueContainer, Variable<array_1d<double, 3> > >() )
    //.def( VariableIndexingPython<AllVariablesDataValueContainer, Variable<vector<double> > >() )
    //.def( VariableIndexingPython<AllVariablesDataValueContainer, Variable<matrix<double> > >() )
    //.def( VariableIndexingPython<AllVariablesDataValueContainer, Variable<ConvectionDiffusionSettings::Pointer > >() )
    //.def( VariableIndexingPython<AllVariablesDataValueContainer, Variable<RadiationSettings::Pointer > >() )
    //.def( VariableIndexingPython<AllVariablesDataValueContainer, VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > > >() )
    //.def( self_ns::str( self ) )
    //;

    class_<DataValueContainer, DataValueContainer::Pointer>( "DataValueContainer" )
    .def( "__len__", &DataValueContainer::Size )
    .def( VariableIndexingPython<DataValueContainer, Variable<std::string> >() )
    .def( VariableIndexingPython<DataValueContainer, Variable<bool> >() )
    .def( VariableIndexingPython<DataValueContainer, Variable<int> >() )
    .def( VariableIndexingPython<DataValueContainer, Variable<double> >() )
    .def( VariableIndexingPython<DataValueContainer, Variable<array_1d<double, 3> > >() )
    .def( VariableIndexingPython<DataValueContainer, Variable<vector<double> > >() )
    .def( VariableIndexingPython<DataValueContainer, Variable<matrix<double> > >() )
    .def( VariableIndexingPython<DataValueContainer, Variable<ConvectionDiffusionSettings::Pointer > >() )
    .def( VariableIndexingPython<DataValueContainer, Variable<RadiationSettings::Pointer > >() )
    .def( VariableIndexingPython<DataValueContainer, VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > > >() )
    .def( VariableIndexingPython<DataValueContainer, Variable<Quaternion<double> > >() )
    .def( self_ns::str( self ) )
    ;

    class_<VariablesListDataValueContainer, VariablesListDataValueContainer::Pointer>( "VariablesListDataValueContainer" )
    .def( "__len__", &VariablesListDataValueContainer::Size )
    .def( VariableIndexingPython<VariablesListDataValueContainer, Variable<std::string> >() )
    .def( VariableIndexingPython<VariablesListDataValueContainer, Variable<bool> >() )
    .def( VariableIndexingPython<VariablesListDataValueContainer, Variable<int> >() )
    .def( VariableIndexingPython<VariablesListDataValueContainer, Variable<double> >() )
    .def( VariableIndexingPython<VariablesListDataValueContainer, Variable<array_1d<double, 3> > >() )
    .def( VariableIndexingPython<VariablesListDataValueContainer, Variable<vector<double> > >() )
    .def( VariableIndexingPython<VariablesListDataValueContainer, Variable<matrix<double> > >() )
    .def( VariableIndexingPython<VariablesListDataValueContainer, VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > > >() )
    .def( VariableIndexingPython<VariablesListDataValueContainer, Variable<Quaternion<double> > >() )
    .def( self_ns::str( self ) )
    ;


    class_<Flags, Flags::Pointer>("Flags",init<>())
    .def(init<Flags>())
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
    .def("__and__", FlagsOr) // this is not an error, the and and or are considered both as add. Pooyan.
    .def( self_ns::str( self ) )
    ;



    KRATOS_REGISTER_IN_PYTHON_FLAG(STRUCTURE);
    KRATOS_REGISTER_IN_PYTHON_FLAG(INTERFACE);
    KRATOS_REGISTER_IN_PYTHON_FLAG(FLUID);
    KRATOS_REGISTER_IN_PYTHON_FLAG(INLET);
    KRATOS_REGISTER_IN_PYTHON_FLAG(OUTLET);
    KRATOS_REGISTER_IN_PYTHON_FLAG(VISITED);
    KRATOS_REGISTER_IN_PYTHON_FLAG(THERMAL);
    KRATOS_REGISTER_IN_PYTHON_FLAG(SELECTED);
    KRATOS_REGISTER_IN_PYTHON_FLAG(BOUNDARY);
    KRATOS_REGISTER_IN_PYTHON_FLAG(SLIP);
    KRATOS_REGISTER_IN_PYTHON_FLAG(CONTACT);
    KRATOS_REGISTER_IN_PYTHON_FLAG(TO_SPLIT);
    KRATOS_REGISTER_IN_PYTHON_FLAG(TO_ERASE);
    KRATOS_REGISTER_IN_PYTHON_FLAG(TO_REFINE);
    KRATOS_REGISTER_IN_PYTHON_FLAG(NEW_ENTITY);
    KRATOS_REGISTER_IN_PYTHON_FLAG(OLD_ENTITY);
    KRATOS_REGISTER_IN_PYTHON_FLAG(ACTIVE);
    KRATOS_REGISTER_IN_PYTHON_FLAG(MODIFIED);
    KRATOS_REGISTER_IN_PYTHON_FLAG(RIGID);
    KRATOS_REGISTER_IN_PYTHON_FLAG(SOLID);
    KRATOS_REGISTER_IN_PYTHON_FLAG(MPI_BOUNDARY);
    KRATOS_REGISTER_IN_PYTHON_FLAG(INTERACTION);
    KRATOS_REGISTER_IN_PYTHON_FLAG(ISOLATED);
    KRATOS_REGISTER_IN_PYTHON_FLAG(MASTER);
    KRATOS_REGISTER_IN_PYTHON_FLAG(SLAVE);
    KRATOS_REGISTER_IN_PYTHON_FLAG(INSIDE);
    KRATOS_REGISTER_IN_PYTHON_FLAG(FREE_SURFACE);
    KRATOS_REGISTER_IN_PYTHON_FLAG(BLOCKED);
    KRATOS_REGISTER_IN_PYTHON_FLAG(MARKER);
    KRATOS_REGISTER_IN_PYTHON_FLAG(PERIODIC);


    AddDeprecatedVariablesToPython();
    AddC2CVariablesToPython();
    AddDEMVariablesToPython(); //TODO: move this to the DEM application
    AddCFDVariablesToPython(); ///@TODO: move variables to CFD application
    AddALEVariablesToPython(); ///@TODO: move variables to ALE application
    AddFSIVariablesToPython(); ///@TODO: move variables to FSI application
    AddMappingVariablesToPython(); ///@TODO: move variables to Mapping application
    AddMATVariablesToPython(); ///@TODO: move variables to CL application
    AddLegacyStructuralAppVarsToPython();

    KRATOS_REGISTER_IN_PYTHON_VARIABLE( SPACE_DIMENSION )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( DOMAIN_SIZE )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( IS_RESTARTED )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( COMPUTE_LUMPED_MASS_MATRIX )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( COMPUTE_DYNAMIC_TANGENT )

    KRATOS_REGISTER_IN_PYTHON_VARIABLE( THERMAL_EXPANSION_COEFFICIENT )

    // These should be moved to applications
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( POWER_LAW_N )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( POWER_LAW_K )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( EQ_STRAIN_RATE )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( YIELD_STRESS )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( MU )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( TAU )


    KRATOS_REGISTER_IN_PYTHON_VARIABLE( CONSTRAINT_LABELS )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( LOAD_LABELS )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( MARKER_LABELS )

    KRATOS_REGISTER_IN_PYTHON_VARIABLE( CONSTRAINT_MESHES )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( LOAD_MESHES )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( MARKER_MESHES )

    KRATOS_REGISTER_IN_PYTHON_VARIABLE( ELEMENTAL_DISTANCES )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( NL_ITERATION_NUMBER )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( FRACTIONAL_STEP )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( STEP )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( TIME )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( START_TIME )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( END_TIME )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( DELTA_TIME )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( PREVIOUS_DELTA_TIME )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( INTERVAL_END_TIME )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( PRINTED_STEP )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( PRINTED_RESTART_STEP )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( RESIDUAL_NORM )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( CONVERGENCE_RATIO )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( TEMPERATURE )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( TEMPERATURE_OLD_IT )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( PRESSURE )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( VISCOSITY_AIR )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( VISCOSITY_WATER )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( ERROR_RATIO )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( TIME_STEPS )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( SCALAR_LAGRANGE_MULTIPLIER )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( VECTOR_LAGRANGE_MULTIPLIER )

    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( ANGULAR_ACCELERATION )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( VELOCITY_LAPLACIAN )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( VELOCITY_LAPLACIAN_RATE )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( VELOCITY_COMPONENT_GRADIENT )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( VELOCITY_X_GRADIENT )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( VELOCITY_Y_GRADIENT )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( VELOCITY_Z_GRADIENT )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( ANGULAR_VELOCITY )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( ACCELERATION )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( VELOCITY )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( DISPLACEMENT )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( ROTATION )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( DELTA_ROTATION )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( TORQUE )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( REACTION )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( BODY_FORCE )

    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( FORCE_RESIDUAL )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( MOMENT_RESIDUAL )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( INTERNAL_FORCE )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( EXTERNAL_FORCE )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( CONTACT_FORCE )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( CONTACT_NORMAL )

    KRATOS_REGISTER_IN_PYTHON_VARIABLE( EXTERNAL_FORCES_VECTOR )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( INTERNAL_FORCES_VECTOR )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( CONTACT_FORCES_VECTOR )
      
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( LINEAR_MOMENTUM )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( ANGULAR_MOMENTUM )

    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( VOLUME_ACCELERATION )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( SEEPAGE_DRAG )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( NORMAL )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( TANGENT_XI )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( TANGENT_ETA )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( FORCE )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( MOMENT )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( FORCE_CM )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( MOMENTUM_CM )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( MOMENTUM )

    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( MASS )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( RHS )

    KRATOS_REGISTER_IN_PYTHON_VARIABLE( WATER_PRESSURE )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( REACTION_WATER_PRESSURE )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( WATER_PRESSURE_ACCELERATION )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( AIR_PRESSURE )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( REACTION_AIR_PRESSURE )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( RHS_WATER )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( RHS_AIR )

    KRATOS_REGISTER_IN_PYTHON_VARIABLE( WEIGHT_FATHER_NODES )

    KRATOS_REGISTER_IN_PYTHON_VARIABLE( INTERNAL_ENERGY )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( STRAIN_ENERGY )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( EXTERNAL_ENERGY )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( KINETIC_ENERGY )

    KRATOS_REGISTER_IN_PYTHON_VARIABLE( LOCAL_INERTIA_TENSOR )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( LOCAL_AXES_MATRIX )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( LOCAL_CONSTITUTIVE_MATRIX )

    //for structural application TO BE REMOVED
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( CONSTITUTIVE_LAW )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( INTERNAL_VARIABLES )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( MATERIAL_PARAMETERS )

    KRATOS_REGISTER_IN_PYTHON_VARIABLE( NEGATIVE_FACE_PRESSURE )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( POSITIVE_FACE_PRESSURE )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( POROSITY )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( DIAMETER )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( LIN_DARCY_COEF )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( NONLIN_DARCY_COEF )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( DRAG_FORCE )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( STRUCTURE_VELOCITY )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( K0 )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( NODAL_VOLUME )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( STATIONARY )

    KRATOS_REGISTER_IN_PYTHON_VARIABLE( FLAG_VARIABLE )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( DISTANCE )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( DISTANCE_GRADIENT )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( INERTIA )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( PERIODIC_PAIR_INDEX )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( PARTITION_INDEX )

    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( LAGRANGE_DISPLACEMENT )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( LAGRANGE_AIR_PRESSURE )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( LAGRANGE_WATER_PRESSURE )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( LAGRANGE_TEMPERATURE )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( INTERNAL_FRICTION_ANGLE )

    KRATOS_REGISTER_IN_PYTHON_VARIABLE( LAGRANGE_DISPLACEMENT )



    // for MultiScale application
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( INITIAL_STRAIN )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( COEFFICIENT_THERMAL_EXPANSION )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( CHARACTERISTIC_LENGTH_MULTIPLIER )

    //for Incompressible Fluid application
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(FRACT_VEL)

    //for ALE application
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( DETERMINANT )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( ELEMENTSHAPE )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( MESH_VELOCITY )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( AUX_MESH_VAR )

    //for AdjointFluidApplication
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( ADJOINT_VELOCITY )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( ADJOINT_ACCELERATION )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( AUX_ADJOINT_ACCELERATION )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( ADJOINT_PRESSURE )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( SHAPE_SENSITIVITY )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( NORMAL_SENSITIVITY )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( NUMBER_OF_NEIGHBOUR_ELEMENTS )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( UPDATE_SENSITIVITIES )

    //for electric application

    KRATOS_REGISTER_IN_PYTHON_VARIABLE( PARTITION_MASK )

    // For MeshingApplication
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( NODAL_ERROR )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( NODAL_ERROR_COMPONENTS )

    // For explicit time integration
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( RESIDUAL_VECTOR )
    
    //for PFEM application TO BE REMOVED
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( NODAL_AREA )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( NODAL_H )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( NORMAL_TO_WALL )
    //KRATOS_REGISTER_IN_PYTHON_VARIABLE(NEIGHBOUR_NODES)
    //KRATOS_REGISTER_IN_PYTHON_VARIABLE(NEIGHBOUR_ELEMENTS)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( FRICTION_COEFFICIENT )

    KRATOS_REGISTER_IN_PYTHON_VARIABLE( BULK_MODULUS )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( SATURATION )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( GRAVITY )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( FACE_LOAD )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( DENSITY_WATER )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( DENSITY_AIR )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( AIR_ENTRY_VALUE )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( FIRST_SATURATION_PARAM )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( SECOND_SATURATION_PARAM )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( PERMEABILITY_WATER )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( PERMEABILITY_AIR )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( BULK_AIR )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( SCALE )

    KRATOS_REGISTER_IN_PYTHON_VARIABLE( TEMP_CONV_PROJ )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( CONVECTION_COEFFICIENT)

    KRATOS_REGISTER_IN_PYTHON_VARIABLE( INSITU_STRESS )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( STRESSES )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( STRAIN )

    KRATOS_REGISTER_IN_PYTHON_VARIABLE( NODAL_MASS )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( AUX_INDEX )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( EXTERNAL_PRESSURE )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( BDF_COEFFICIENTS )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( ROTATION_CENTER )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( VELOCITY_PERIOD )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( ANGULAR_VELOCITY_PERIOD )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( IDENTIFIER )


    //for xfem application
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( CRACK_OPENING )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( CRACK_TRANSLATION )

    KRATOS_REGISTER_IN_PYTHON_VARIABLE(ARRHENIUS)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(ARRHENIUSAUX)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(ARRHENIUSAUX_)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(PRESSUREAUX)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(NODAL_MAUX)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(NODAL_PAUX)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(FACE_HEAT_FLUX)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(HEAT_FLUX)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(REACTION_FLUX)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(TC)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(CONDUCTIVITY)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(SPECIFIC_HEAT)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(MATERIAL_VARIABLE)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(FUEL)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(YO)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(YF)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(YI)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(Y1)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(Y2)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(YP)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(EMISSIVITY)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(ENTHALPY)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(MIXTURE_FRACTION)

    KRATOS_REGISTER_IN_PYTHON_VARIABLE(YCH4)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(YO2)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(YCO2)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(YH2O)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(YN2)

    KRATOS_REGISTER_IN_PYTHON_VARIABLE(INCIDENT_RADIATION_FUNCTION)

    KRATOS_REGISTER_IN_PYTHON_VARIABLE(ABSORPTION_COEFFICIENT)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(STEFAN_BOLTZMANN_CONSTANT)
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( DIRECTION )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(NODAL_SWITCH)
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(Y)
    
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( LOCAL_AXIS_1 )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( LOCAL_AXIS_2 )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( LOCAL_AXIS_3 ) 
      
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( SWITCH_TEMPERATURE )

    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(EMBEDDED_VELOCITY)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( REFINEMENT_LEVEL )

    KRATOS_REGISTER_IN_PYTHON_VARIABLE( AIR_SOUND_VELOCITY )

    // for Vulcan application
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( LATENT_HEAT )

    KRATOS_REGISTER_IN_PYTHON_VARIABLE( ENRICHED_PRESSURES )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( INITIAL_PENALTY )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( DP_EPSILON )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( DP_ALPHA1 )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( DP_K )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( LAMBDA )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( MIN_DT )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( MAX_DT )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( WET_VOLUME)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( CUTTED_AREA)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( NET_INPUT_MATERIAL)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( MIU )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( SCALE_FACTOR )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( NORMAL_CONTACT_STRESS )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( TANGENTIAL_CONTACT_STRESS )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( STABILIZATION_FACTOR )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( NEWMARK_BETA )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( NEWMARK_GAMMA )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( BOSSAK_ALPHA )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( EQUILIBRIUM_POINT )      
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( AIR_SOUND_VELOCITY )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( WATER_SOUND_VELOCITY )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( ACTIVATION_LEVEL )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( FIRST_TIME_STEP )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( QUASI_STATIC_ANALYSIS )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( FRACTIONAL_STEP )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( LOAD_RESTART )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( RIGID_BODY_ID )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( NORMAL_STRESS )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( TANGENTIAL_STRESS )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( PENALTY )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( AMBIENT_TEMPERATURE )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( VEL_ART_VISC )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( PR_ART_VISC )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( SOUND_VELOCITY )

    KRATOS_REGISTER_IN_PYTHON_VARIABLE( SEARCH_RADIUS )

    KRATOS_REGISTER_IN_PYTHON_VARIABLE( INTEGRATION_WEIGHT )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( INTEGRATION_COORDINATES )


    class_< ConvectionDiffusionSettings, ConvectionDiffusionSettings::Pointer, boost::noncopyable >	("ConvectionDiffusionSettings", init<	>() )
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

    .def("GetDensityVariable",&ConvectionDiffusionSettings::GetDensityVariable, return_internal_reference<>() )
    .def("GetDiffusionVariable",&ConvectionDiffusionSettings::GetDiffusionVariable, return_internal_reference<>() )
    .def("GetUnknownVariable",&ConvectionDiffusionSettings::GetUnknownVariable, return_internal_reference<>() )
    .def("GetVolumeSourceVariable",&ConvectionDiffusionSettings::GetVolumeSourceVariable, return_internal_reference<>() )
    .def("GetSurfaceSourceVariable",&ConvectionDiffusionSettings::GetSurfaceSourceVariable, return_internal_reference<>() )
    .def("GetProjectionVariable",&ConvectionDiffusionSettings::GetProjectionVariable, return_internal_reference<>() )
    .def("GetMeshVelocityVariable",&ConvectionDiffusionSettings::GetMeshVelocityVariable, return_internal_reference<>() )
    .def("GetConvectionVariable",&ConvectionDiffusionSettings::GetConvectionVariable, return_internal_reference<>() )
    .def("GetTransferCoefficientVariable",&ConvectionDiffusionSettings::GetTransferCoefficientVariable, return_internal_reference<>())
    .def("GetSpecificHeatVariable",&ConvectionDiffusionSettings::GetSpecificHeatVariable, return_internal_reference<>() )
    .def("GetVelocityVariable",&ConvectionDiffusionSettings::GetVelocityVariable, return_internal_reference<>() )
    .def("GetReactionVariable",&ConvectionDiffusionSettings::GetReactionVariable, return_internal_reference<>() )

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

    class_< RadiationSettings, RadiationSettings::Pointer, boost::noncopyable >	("RadiationSettings", init<	>() )
    .def("SetDensityVariable",&RadiationSettings::SetDensityVariable)
    .def("SetDiffusionVariable",&RadiationSettings::SetDiffusionVariable)
    .def("SetUnknownVariable",&RadiationSettings::SetUnknownVariable)
    .def("SetVolumeSourceVariable",&RadiationSettings::SetVolumeSourceVariable)
    .def("SetSurfaceSourceVariable",&RadiationSettings::SetSurfaceSourceVariable)
    .def("SetProjectionVariable",&RadiationSettings::SetProjectionVariable)
    .def("SetMeshVelocityVariable",&RadiationSettings::SetMeshVelocityVariable)
    .def("GetDensityVariable",&RadiationSettings::GetDensityVariable, return_internal_reference<>() )
    .def("GetDiffusionVariable",&RadiationSettings::GetDiffusionVariable, return_internal_reference<>() )
    .def("GetUnknownVariable",&RadiationSettings::GetUnknownVariable, return_internal_reference<>() )
    .def("GetVolumeSourceVariable",&RadiationSettings::GetVolumeSourceVariable, return_internal_reference<>() )
    .def("GetSurfaceSourceVariable",&RadiationSettings::GetSurfaceSourceVariable, return_internal_reference<>() )
    //.def("GetSurfaceSourceVariable",&RadiationSettings::GetSurfaceSourceVariable, return_internal_reference<>() )
    .def("GetProjectionVariable",&RadiationSettings::GetProjectionVariable, return_internal_reference<>() )
    .def("GetMeshVelocityVariable",&RadiationSettings::GetMeshVelocityVariable, return_internal_reference<>() )
    ;
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(CONVECTION_DIFFUSION_SETTINGS)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(RADIATION_SETTINGS)

    // KRATOS_REGISTER_IN_PYTHON_VARIABLE(CONVECTION_DIFFUSION_SETTINGS)

}
} // namespace Python.
} // Namespace Kratos
