/*
==============================================================================
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain

Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNER.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
*/

//
//   Project Name:        Kratos
//   Last modified by:    $Author: mengmeng $
//   Date:                $Date: 2009-02-26 14:28:21 $
//   Revision:            $Revision: 1.13 $
//
//


// System includes

// External includes
#include <boost/python.hpp>


// Project includes
#include "includes/define.h"
#include "includes/ublas_interface.h"
#include "containers/data_value_container.h"
#include "containers/variables_list_data_value_container.h"
#include "containers/fix_data_value_container.h"
#include "containers/vector_component_adaptor.h"
#include "includes/variables.h"
#include "includes/constitutive_law.h"
#include "python/variable_indexing_python.h"
#include "python/vector_python_interface.h"
#include "python/vector_scalar_operator_python.h"
#include "python/vector_vector_operator_python.h"

#include "includes/convection_diffusion_settings.h"
#include "includes/radiation_settings.h"


#ifdef KRATOS_REGISTER_IN_PYTHON_VARIABLE
#undef KRATOS_REGISTER_IN_PYTHON_VARIABLE
#endif
#define KRATOS_REGISTER_IN_PYTHON_VARIABLE(variable) \
 scope().attr(#variable) = boost::ref(variable);

#ifdef KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS
#undef KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS
#endif
#define KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(name) \
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(name) \
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(name##_X) \
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(name##_Y) \
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(name##_Z)

namespace Kratos
{
    namespace Python
    {
        using namespace boost::python;

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
            VectorPythonInterface<array_1d<double, 3>, Array1DModifier<array_1d<double, 3> > >::CreateInterface( "Array3" )
            .def( init<vector_expression<array_1d<double, 3> > >() )
            .def( VectorScalarOperatorPython<array_1d<double, 3>, double, array_1d<double, 3> >() )
            .def( VectorVectorOperatorPython<array_1d<double, 3>, zero_vector<double>, array_1d<double, 3> >() )
            .def( VectorVectorOperatorPython<array_1d<double, 3>, unit_vector<double>, array_1d<double, 3> >() )
            .def( VectorVectorOperatorPython<array_1d<double, 3>, scalar_vector<double>, array_1d<double, 3> >() )
            .def( VectorVectorOperatorPython<array_1d<double, 3>, mapped_vector<double>, array_1d<double, 3> >() )
            ;

            class_<VariableData>( "VariableData", no_init )
            .def( self_ns::str( self ) )
            ;

            class_<Variable<bool>, bases<VariableData>, boost::noncopyable >( "BoolVariable", no_init )
            .def( self_ns::str( self ) )
            ;

            class_<Variable<int>, bases<VariableData>, boost::noncopyable >( "IntegerVariable", no_init )
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
	    

            class_<DataValueContainer, DataValueContainer::Pointer>( "DataValueContainer" )
            .def( "__len__", &DataValueContainer::Size )
            .def( VariableIndexingPython<DataValueContainer, Variable<int> >() )
            .def( VariableIndexingPython<DataValueContainer, Variable<double> >() )
            .def( VariableIndexingPython<DataValueContainer, Variable<array_1d<double, 3> > >() )
            .def( VariableIndexingPython<DataValueContainer, Variable<vector<double> > >() )
            .def( VariableIndexingPython<DataValueContainer, Variable<matrix<double> > >() )
            .def( VariableIndexingPython<DataValueContainer, Variable<ConvectionDiffusionSettings::Pointer > >() )
	    .def( VariableIndexingPython<DataValueContainer, Variable<RadiationSettings::Pointer > >() )
            .def( VariableIndexingPython<DataValueContainer, VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > > >() )
            .def( self_ns::str( self ) )
            ;

            class_<VariablesListDataValueContainer, VariablesListDataValueContainer::Pointer>( "VariablesListDataValueContainer" )
            .def( "__len__", &VariablesListDataValueContainer::Size )
            .def( VariableIndexingPython<VariablesListDataValueContainer, Variable<int> >() )
            .def( VariableIndexingPython<VariablesListDataValueContainer, Variable<double> >() )
            .def( VariableIndexingPython<VariablesListDataValueContainer, Variable<array_1d<double, 3> > >() )
            .def( VariableIndexingPython<VariablesListDataValueContainer, Variable<vector<double> > >() )
            .def( VariableIndexingPython<VariablesListDataValueContainer, Variable<matrix<double> > >() )
            .def( VariableIndexingPython<VariablesListDataValueContainer,    VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > > >() )
            .def( self_ns::str( self ) )
            ;
	    
	    



            KRATOS_REGISTER_IN_PYTHON_VARIABLE( OSS_SWITCH )
            KRATOS_REGISTER_IN_PYTHON_VARIABLE( DYNAMIC_TAU )
            KRATOS_REGISTER_IN_PYTHON_VARIABLE( ERASE_FLAG )
            KRATOS_REGISTER_IN_PYTHON_VARIABLE( NL_ITERATION_NUMBER )
            KRATOS_REGISTER_IN_PYTHON_VARIABLE( FRACTIONAL_STEP )
            KRATOS_REGISTER_IN_PYTHON_VARIABLE( TIME )
            KRATOS_REGISTER_IN_PYTHON_VARIABLE( TIME_STEPS )
            KRATOS_REGISTER_IN_PYTHON_VARIABLE( START_TIME )
            KRATOS_REGISTER_IN_PYTHON_VARIABLE( END_TIME )
            KRATOS_REGISTER_IN_PYTHON_VARIABLE( DELTA_TIME )
            KRATOS_REGISTER_IN_PYTHON_VARIABLE( TEMPERATURE )
            KRATOS_REGISTER_IN_PYTHON_VARIABLE( PRESSURE )
            KRATOS_REGISTER_IN_PYTHON_VARIABLE( DENSITY )
            KRATOS_REGISTER_IN_PYTHON_VARIABLE( VISCOSITY )
            KRATOS_REGISTER_IN_PYTHON_VARIABLE( VISCOSITY_AIR )
            KRATOS_REGISTER_IN_PYTHON_VARIABLE( VISCOSITY_WATER )
            KRATOS_REGISTER_IN_PYTHON_VARIABLE( ERROR_RATIO )
            KRATOS_REGISTER_IN_PYTHON_VARIABLE( TAU )

            KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( ANGULAR_ACCELERATION )
            KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( ANGULAR_VELOCITY )
            KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( ACCELERATION )
            KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( VELOCITY )
            KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( DISPLACEMENT )
            KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( ROTATION )
            KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( MOMENTUM )
            KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( REACTION )
            KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( BODY_FORCE )
            KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( SEEPAGE_DRAG )
            KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( MESH_VELOCITY )
            KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( NORMAL )
            KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( FORCE )
            KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( MOMENT )
            KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( FORCE_CM )
            KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( MOMENTUM_CM )


            KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( MASS )
            KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( RHS )

            KRATOS_REGISTER_IN_PYTHON_VARIABLE( WATER_PRESSURE )
            KRATOS_REGISTER_IN_PYTHON_VARIABLE( REACTION_WATER_PRESSURE )
            KRATOS_REGISTER_IN_PYTHON_VARIABLE( WATER_PRESSURE_DT )
            KRATOS_REGISTER_IN_PYTHON_VARIABLE( WATER_PRESSURE_ACCELERATION )
            KRATOS_REGISTER_IN_PYTHON_VARIABLE( AIR_PRESSURE )
            KRATOS_REGISTER_IN_PYTHON_VARIABLE( REACTION_AIR_PRESSURE )
            KRATOS_REGISTER_IN_PYTHON_VARIABLE( AIR_PRESSURE_DT )
            KRATOS_REGISTER_IN_PYTHON_VARIABLE( AIR_PRESSURE_ACCELERATION )
            KRATOS_REGISTER_IN_PYTHON_VARIABLE( RHS_WATER )
            KRATOS_REGISTER_IN_PYTHON_VARIABLE( RHS_AIR )

            //for structural application TO BE REMOVED
            KRATOS_REGISTER_IN_PYTHON_VARIABLE( INTERNAL_VARIABLES )
            KRATOS_REGISTER_IN_PYTHON_VARIABLE( MATERIAL_PARAMETERS )
            KRATOS_REGISTER_IN_PYTHON_VARIABLE( WRINKLING_APPROACH )
            KRATOS_REGISTER_IN_PYTHON_VARIABLE( GREEN_LAGRANGE_STRAIN_TENSOR )
            KRATOS_REGISTER_IN_PYTHON_VARIABLE( PK2_STRESS_TENSOR )
            KRATOS_REGISTER_IN_PYTHON_VARIABLE( CAUCHY_STRESS_TENSOR )
            KRATOS_REGISTER_IN_PYTHON_VARIABLE( AUXILIARY_MATRIX_1 )
            KRATOS_REGISTER_IN_PYTHON_VARIABLE( YOUNG_MODULUS )
            KRATOS_REGISTER_IN_PYTHON_VARIABLE( POISSON_RATIO )
            KRATOS_REGISTER_IN_PYTHON_VARIABLE( THICKNESS )
            KRATOS_REGISTER_IN_PYTHON_VARIABLE( NEGATIVE_FACE_PRESSURE )
            KRATOS_REGISTER_IN_PYTHON_VARIABLE( POSITIVE_FACE_PRESSURE )
            KRATOS_REGISTER_IN_PYTHON_VARIABLE( IS_POROUS )
            KRATOS_REGISTER_IN_PYTHON_VARIABLE( IS_WATER )
            KRATOS_REGISTER_IN_PYTHON_VARIABLE( POROSITY )
            KRATOS_REGISTER_IN_PYTHON_VARIABLE( DIAMETER )
            KRATOS_REGISTER_IN_PYTHON_VARIABLE( FLAG_VARIABLE )
            KRATOS_REGISTER_IN_PYTHON_VARIABLE( DISTANCE )
            KRATOS_REGISTER_IN_PYTHON_VARIABLE( INERTIA )
            KRATOS_REGISTER_IN_PYTHON_VARIABLE( PARTITION_INDEX )
            KRATOS_REGISTER_IN_PYTHON_VARIABLE( YIELD_STRESS )

            KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( LAGRANGE_DISPLACEMENT )
            KRATOS_REGISTER_IN_PYTHON_VARIABLE( LAGRANGE_AIR_PRESSURE )
            KRATOS_REGISTER_IN_PYTHON_VARIABLE( LAGRANGE_WATER_PRESSURE )
            KRATOS_REGISTER_IN_PYTHON_VARIABLE( LAGRANGE_TEMPERATURE )
            KRATOS_REGISTER_IN_PYTHON_VARIABLE( INTERNAL_FRICTION_ANGLE )

            //for ALE application TO BE REMOVED
            KRATOS_REGISTER_IN_PYTHON_VARIABLE( IS_INTERFACE )
            KRATOS_REGISTER_IN_PYTHON_VARIABLE( IS_VISITED )
            KRATOS_REGISTER_IN_PYTHON_VARIABLE( IS_LAGRANGIAN_INLET )

            //for electric application
            KRATOS_REGISTER_IN_PYTHON_VARIABLE( ELECTRIC_POTENTIAL )
            KRATOS_REGISTER_IN_PYTHON_VARIABLE( ELECTRIC_POTENTIAL_DT )
            KRATOS_REGISTER_IN_PYTHON_VARIABLE( ELECTRIC_POTENTIAL_DT2 )
            KRATOS_REGISTER_IN_PYTHON_VARIABLE( ELECTRIC_POTENTIAL_NULL )
            KRATOS_REGISTER_IN_PYTHON_VARIABLE( ELECTRIC_POTENTIAL_EINS )
            KRATOS_REGISTER_IN_PYTHON_VARIABLE( ELECTRIC_POTENTIAL_NULL_DT )
            KRATOS_REGISTER_IN_PYTHON_VARIABLE( ELECTRIC_POTENTIAL_EINS_DT )
            KRATOS_REGISTER_IN_PYTHON_VARIABLE( ELECTRIC_POTENTIAL_NULL_DT2 )
            KRATOS_REGISTER_IN_PYTHON_VARIABLE( ELECTRIC_POTENTIAL_EINS_DT2 )
            KRATOS_REGISTER_IN_PYTHON_VARIABLE( LAGRANGE_ELECTRIC_POTENTIAL )
            KRATOS_REGISTER_IN_PYTHON_VARIABLE( CONCENTRATION )
            KRATOS_REGISTER_IN_PYTHON_VARIABLE( CONCENTRATION_DT )
            KRATOS_REGISTER_IN_PYTHON_VARIABLE( CONCENTRATION_DT2 )
            KRATOS_REGISTER_IN_PYTHON_VARIABLE( CONCENTRATION_NULL )
            KRATOS_REGISTER_IN_PYTHON_VARIABLE( CONCENTRATION_EINS )
            KRATOS_REGISTER_IN_PYTHON_VARIABLE( CONCENTRATION_NULL_DT )
            KRATOS_REGISTER_IN_PYTHON_VARIABLE( CONCENTRATION_EINS_DT )
            KRATOS_REGISTER_IN_PYTHON_VARIABLE( CONCENTRATION_NULL_DT2 )
            KRATOS_REGISTER_IN_PYTHON_VARIABLE( CONCENTRATION_EINS_DT2 )
            KRATOS_REGISTER_IN_PYTHON_VARIABLE( LAGRANGE_CONCENTRATION )

            //for PFEM application TO BE REMOVED
            KRATOS_REGISTER_IN_PYTHON_VARIABLE( NODAL_AREA )
            KRATOS_REGISTER_IN_PYTHON_VARIABLE( NODAL_H )
            KRATOS_REGISTER_IN_PYTHON_VARIABLE( IS_STRUCTURE )
            KRATOS_REGISTER_IN_PYTHON_VARIABLE( IS_FLUID )
            KRATOS_REGISTER_IN_PYTHON_VARIABLE( IS_BOUNDARY )
            KRATOS_REGISTER_IN_PYTHON_VARIABLE( IS_FREE_SURFACE )
            KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( NORMAL_TO_WALL )
            //KRATOS_REGISTER_IN_PYTHON_VARIABLE(NEIGHBOUR_NODES);
            //KRATOS_REGISTER_IN_PYTHON_VARIABLE(NEIGHBOUR_ELEMENTS);
            KRATOS_REGISTER_IN_PYTHON_VARIABLE( IS_EROSIONABLE )
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
            KRATOS_REGISTER_IN_PYTHON_VARIABLE( INSITU_STRESS )
            KRATOS_REGISTER_IN_PYTHON_VARIABLE( ELASTIC_LEFT_CAUCHY_GREEN_OLD )

            KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( FRACT_VEL )
            KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( PRESS_PROJ )
            KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( CONV_PROJ )
            KRATOS_REGISTER_IN_PYTHON_VARIABLE( NODAL_MASS )
            KRATOS_REGISTER_IN_PYTHON_VARIABLE( AUX_INDEX )
            KRATOS_REGISTER_IN_PYTHON_VARIABLE( EXTERNAL_PRESSURE )
            KRATOS_REGISTER_IN_PYTHON_VARIABLE( PRESSURE_OLD_IT )
            KRATOS_REGISTER_IN_PYTHON_VARIABLE( BDF_COEFFICIENTS )

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
	    KRATOS_REGISTER_IN_PYTHON_VARIABLE(TC)
            KRATOS_REGISTER_IN_PYTHON_VARIABLE(IS_BURN )
	    KRATOS_REGISTER_IN_PYTHON_VARIABLE(CONDUCTIVITY)
	    KRATOS_REGISTER_IN_PYTHON_VARIABLE(SPECIFIC_HEAT)
	    KRATOS_REGISTER_IN_PYTHON_VARIABLE(IS_DRIPPING )
	    KRATOS_REGISTER_IN_PYTHON_VARIABLE( IS_PERMANENT)
	    KRATOS_REGISTER_IN_PYTHON_VARIABLE(MATERIAL_VARIABLE)
	    KRATOS_REGISTER_IN_PYTHON_VARIABLE(IS_WALL )
	    KRATOS_REGISTER_IN_PYTHON_VARIABLE(FUEL)
	    KRATOS_REGISTER_IN_PYTHON_VARIABLE(YO)
	    KRATOS_REGISTER_IN_PYTHON_VARIABLE(YF)
	    KRATOS_REGISTER_IN_PYTHON_VARIABLE(YI)
	    KRATOS_REGISTER_IN_PYTHON_VARIABLE(M)
	    KRATOS_REGISTER_IN_PYTHON_VARIABLE(Y1)
	    KRATOS_REGISTER_IN_PYTHON_VARIABLE(Y2)
	    KRATOS_REGISTER_IN_PYTHON_VARIABLE(YP) 
	    KRATOS_REGISTER_IN_PYTHON_VARIABLE(RADIATIVE_INTENSITY)
	    KRATOS_REGISTER_IN_PYTHON_VARIABLE(RADIATIVE_INTENSITY_1)
            KRATOS_REGISTER_IN_PYTHON_VARIABLE(RADIATIVE_INTENSITY_2)
            KRATOS_REGISTER_IN_PYTHON_VARIABLE(RADIATIVE_INTENSITY_3)
            KRATOS_REGISTER_IN_PYTHON_VARIABLE(RADIATIVE_INTENSITY_4)
	    KRATOS_REGISTER_IN_PYTHON_VARIABLE(RADIATIVE_INTENSITY_5)
            KRATOS_REGISTER_IN_PYTHON_VARIABLE(RADIATIVE_INTENSITY_6)
            KRATOS_REGISTER_IN_PYTHON_VARIABLE(RADIATIVE_INTENSITY_7)
            KRATOS_REGISTER_IN_PYTHON_VARIABLE(RADIATIVE_INTENSITY_8)
            KRATOS_REGISTER_IN_PYTHON_VARIABLE(RADIATIVE_INTENSITY_9)
            KRATOS_REGISTER_IN_PYTHON_VARIABLE(RADIATIVE_INTENSITY_10)
            KRATOS_REGISTER_IN_PYTHON_VARIABLE(RADIATIVE_INTENSITY_11)
            KRATOS_REGISTER_IN_PYTHON_VARIABLE(RADIATIVE_INTENSITY_12)
	    KRATOS_REGISTER_IN_PYTHON_VARIABLE(RADIATIVE_INTENSITY_13)
            KRATOS_REGISTER_IN_PYTHON_VARIABLE(RADIATIVE_INTENSITY_14)
            KRATOS_REGISTER_IN_PYTHON_VARIABLE(RADIATIVE_INTENSITY_15)
            KRATOS_REGISTER_IN_PYTHON_VARIABLE(RADIATIVE_INTENSITY_16)
            KRATOS_REGISTER_IN_PYTHON_VARIABLE(RADIATIVE_INTENSITY_17)
            KRATOS_REGISTER_IN_PYTHON_VARIABLE(RADIATIVE_INTENSITY_18)
            KRATOS_REGISTER_IN_PYTHON_VARIABLE(RADIATIVE_INTENSITY_19)
            KRATOS_REGISTER_IN_PYTHON_VARIABLE(RADIATIVE_INTENSITY_20)
            KRATOS_REGISTER_IN_PYTHON_VARIABLE(RADIATIVE_INTENSITY_21)
            KRATOS_REGISTER_IN_PYTHON_VARIABLE(RADIATIVE_INTENSITY_22)
            KRATOS_REGISTER_IN_PYTHON_VARIABLE(RADIATIVE_INTENSITY_23)
            KRATOS_REGISTER_IN_PYTHON_VARIABLE(RADIATIVE_INTENSITY_24)
	    KRATOS_REGISTER_IN_PYTHON_VARIABLE(EMISSIVITY)
	    KRATOS_REGISTER_IN_PYTHON_VARIABLE(ENTHALPY)
            KRATOS_REGISTER_IN_PYTHON_VARIABLE(MIXTURE_FRACTION)
            KRATOS_REGISTER_IN_PYTHON_VARIABLE(rhoD)
	    KRATOS_REGISTER_IN_PYTHON_VARIABLE(Yfuel)
            KRATOS_REGISTER_IN_PYTHON_VARIABLE(Yox)
            KRATOS_REGISTER_IN_PYTHON_VARIABLE(Ypr)
            KRATOS_REGISTER_IN_PYTHON_VARIABLE(Hfuel)
            KRATOS_REGISTER_IN_PYTHON_VARIABLE(Hpr)
            KRATOS_REGISTER_IN_PYTHON_VARIABLE(Hpr1)
            KRATOS_REGISTER_IN_PYTHON_VARIABLE(Hox)



	    KRATOS_REGISTER_IN_PYTHON_VARIABLE(INCIDENT_RADIATION_FUNCTION)


	    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( MESH_VELOCITY1 )
            KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( MESH_VELOCITY2 )
            KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( MESH_VELOCITY3 )
            KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( MESH_VELOCITY4 )
	    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( MESH_VELOCITY5 )
            KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( MESH_VELOCITY6 )
            KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( MESH_VELOCITY7 )
            KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( MESH_VELOCITY8 )
            KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( MESH_VELOCITY9 )
            KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( MESH_VELOCITY10 )
            KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( MESH_VELOCITY11 )
            KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( MESH_VELOCITY12 )
	    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( MESH_VELOCITY13 )
	    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( MESH_VELOCITY14 )
            KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( MESH_VELOCITY15 )
            KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( MESH_VELOCITY16 )
            KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( MESH_VELOCITY17 )
            KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( MESH_VELOCITY18 )
            KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( MESH_VELOCITY19 )
            KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( MESH_VELOCITY20 )
            KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( MESH_VELOCITY21 )
            KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( MESH_VELOCITY22 )
            KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( MESH_VELOCITY23 )
            KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( MESH_VELOCITY24 )



	    KRATOS_REGISTER_IN_PYTHON_VARIABLE(ABSORPTION_COEFFICIENT)
            KRATOS_REGISTER_IN_PYTHON_VARIABLE(STEFAN_BOLTZMANN_CONSTANT)
            KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( DIRECTION )
            KRATOS_REGISTER_IN_PYTHON_VARIABLE(NODAL_SWITCH)
	    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(Y)
	    KRATOS_REGISTER_IN_PYTHON_VARIABLE( SWITCH_TEMPERATURE )

           KRATOS_REGISTER_IN_PYTHON_VARIABLE( REFINEMENT_LEVEL )
           KRATOS_REGISTER_IN_PYTHON_VARIABLE( IS_INACTIVE )
           KRATOS_REGISTER_IN_PYTHON_VARIABLE( IS_DUPLICATED )
           KRATOS_REGISTER_IN_PYTHON_VARIABLE( SPLIT_NODAL )
           KRATOS_REGISTER_IN_PYTHON_VARIABLE( SPLIT_ELEMENT )
           
           
	    class_< ConvectionDiffusionSettings, ConvectionDiffusionSettings::Pointer, boost::noncopyable >	("ConvectionDiffusionSettings", init<	>() )
 			  .def("SetDensityVariable",&ConvectionDiffusionSettings::SetDensityVariable)
			  .def("SetDiffusionVariable",&ConvectionDiffusionSettings::SetDiffusionVariable)
			  .def("SetUnknownVariable",&ConvectionDiffusionSettings::SetUnknownVariable)
			  .def("SetVolumeSourceVariable",&ConvectionDiffusionSettings::SetVolumeSourceVariable)
			  .def("SetSurfaceSourceVariable",&ConvectionDiffusionSettings::SetSurfaceSourceVariable)
			  .def("SetProjectionVariable",&ConvectionDiffusionSettings::SetConvectionVariable)
			  .def("SetMeshVelocityVariable",&ConvectionDiffusionSettings::SetMeshVelocityVariable)
 			  .def("GetDensityVariable",&ConvectionDiffusionSettings::GetDensityVariable, return_internal_reference<>() )
			  .def("GetDiffusionVariable",&ConvectionDiffusionSettings::GetDiffusionVariable, return_internal_reference<>() )
			  .def("GetUnknownVariable",&ConvectionDiffusionSettings::GetUnknownVariable, return_internal_reference<>() )
			  .def("GetVolumeSourceVariable",&ConvectionDiffusionSettings::GetVolumeSourceVariable, return_internal_reference<>() )
			  .def("GetSurfaceSourceVariable",&ConvectionDiffusionSettings::GetSurfaceSourceVariable, return_internal_reference<>() )
			  .def("GetSurfaceSourceVariable",&ConvectionDiffusionSettings::GetSurfaceSourceVariable, return_internal_reference<>() )
			  .def("GetProjectionVariable",&ConvectionDiffusionSettings::GetConvectionVariable, return_internal_reference<>() )
			  .def("GetMeshVelocityVariable",&ConvectionDiffusionSettings::GetMeshVelocityVariable, return_internal_reference<>() )
			  ;	
 		class_< RadiationSettings, RadiationSettings::Pointer, boost::noncopyable >	("RadiationSettings", init<	>() )
 			  .def("SetDensityVariable",&RadiationSettings::SetDensityVariable)
			  .def("SetDiffusionVariable",&RadiationSettings::SetDiffusionVariable)
			  .def("SetUnknownVariable",&RadiationSettings::SetUnknownVariable)
			  .def("SetVolumeSourceVariable",&RadiationSettings::SetVolumeSourceVariable)
			  .def("SetSurfaceSourceVariable",&RadiationSettings::SetSurfaceSourceVariable)
			  .def("SetProjectionVariable",&RadiationSettings::SetConvectionVariable)
			  .def("SetMeshVelocityVariable",&RadiationSettings::SetMeshVelocityVariable)
 			  .def("GetDensityVariable",&RadiationSettings::GetDensityVariable, return_internal_reference<>() )
			  .def("GetDiffusionVariable",&RadiationSettings::GetDiffusionVariable, return_internal_reference<>() )
			  .def("GetUnknownVariable",&RadiationSettings::GetUnknownVariable, return_internal_reference<>() )
			  .def("GetVolumeSourceVariable",&RadiationSettings::GetVolumeSourceVariable, return_internal_reference<>() )
			  .def("GetSurfaceSourceVariable",&RadiationSettings::GetSurfaceSourceVariable, return_internal_reference<>() )
			  .def("GetSurfaceSourceVariable",&RadiationSettings::GetSurfaceSourceVariable, return_internal_reference<>() )
			  .def("GetProjectionVariable",&RadiationSettings::GetConvectionVariable, return_internal_reference<>() )
			  .def("GetMeshVelocityVariable",&RadiationSettings::GetMeshVelocityVariable, return_internal_reference<>() )
			  ;	
           KRATOS_REGISTER_IN_PYTHON_VARIABLE(CONVECTION_DIFFUSION_SETTINGS)
	   KRATOS_REGISTER_IN_PYTHON_VARIABLE(RADIATION_SETTINGS)

	     // KRATOS_REGISTER_IN_PYTHON_VARIABLE(CONVECTION_DIFFUSION_SETTINGS)

        }
    }  // namespace Python.
} // Namespace Kratos

#undef KRATOS_REGISTER_IN_PYTHON_VARIABLE
#undef KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS
