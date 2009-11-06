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
#include "python/variable_indexing_python.h"
#include "python/vector_python_interface.h"
#include "python/vector_scalar_operator_python.h"
#include "python/vector_vector_operator_python.h"


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
    static void Resize(TContainerType& ThisContainer, typename TContainerType::size_type NewSize)
    {
    }
    static void MoveSlice(TContainerType& ThisContainer, index_type Index, index_type From, index_type To)
    {
    }
  };
  void  AddContainersToPython()
  {

    VectorPythonInterface<array_1d<double, 3>, Array1DModifier<array_1d<double, 3> > >::CreateInterface("Array3")
      .def(init<vector_expression<array_1d<double, 3> > >())
      .def(VectorScalarOperatorPython<array_1d<double, 3>, double, array_1d<double, 3> >())
      .def(VectorVectorOperatorPython<array_1d<double, 3>, zero_vector<double>, array_1d<double, 3> >())
      .def(VectorVectorOperatorPython<array_1d<double, 3>, unit_vector<double>, array_1d<double, 3> >())
      .def(VectorVectorOperatorPython<array_1d<double, 3>, scalar_vector<double>, array_1d<double, 3> >())
      .def(VectorVectorOperatorPython<array_1d<double, 3>, mapped_vector<double>, array_1d<double, 3> >())
       ;

       class_<VariableData>("VariableData", no_init)
	 .def(self_ns::str(self))
	 ;

       class_<Variable<int>, bases<VariableData>, boost::noncopyable >("IntegerVariable", no_init)
	 .def(self_ns::str(self))
	 ;

       class_<Variable<double>, bases<VariableData>, boost::noncopyable >("DoubleVariable", no_init)
	 .def(self_ns::str(self))
	 ;

       class_<Variable<vector<double> >, bases<VariableData>, boost::noncopyable >("VectorVariable", no_init)
	 .def(self_ns::str(self))
	 ;

       class_<Variable<array_1d<double,3> >, bases<VariableData>, boost::noncopyable >("Array1DVariable", no_init)
	 .def(self_ns::str(self))
	 ;

       class_<Variable<matrix<double> >, bases<VariableData>, boost::noncopyable >("MatrixVariable", no_init)
	 .def(self_ns::str(self))
	 ;
       
       class_<VariableComponent<VectorComponentAdaptor<vector<double> > >, bases<VariableData>, boost::noncopyable >("VectorComponentVariable", no_init)
	 .def(self_ns::str(self))
	 ;

	   class_<VariableComponent<VectorComponentAdaptor<array_1d<double,3> > >, bases<VariableData>, boost::noncopyable >("Array1DComponentVariable", no_init)
	 .def(self_ns::str(self))
	 ;

       class_<DataValueContainer, DataValueContainer::Pointer>("DataValueContainer")
		   .def("__len__", &DataValueContainer::Size)
		   .def(VariableIndexingPython<DataValueContainer, Variable<int> >())
		   .def(VariableIndexingPython<DataValueContainer, Variable<double> >())
		   .def(VariableIndexingPython<DataValueContainer, Variable<array_1d<double, 3> > >())
		   .def(VariableIndexingPython<DataValueContainer, Variable<vector<double> > >())
		   .def(VariableIndexingPython<DataValueContainer, Variable<matrix<double> > >())
		   .def(VariableIndexingPython<DataValueContainer, VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > > >())
	       .def(self_ns::str(self))
      ;

  /*     class_<FixDataValueContainer, FixDataValueContainer::Pointer>("FixDataValueContainer")
		   .def("__len__", &FixDataValueContainer::Size)
		   .def(VariableIndexingPython<FixDataValueContainer, Variable<int> >())
		   .def(VariableIndexingPython<FixDataValueContainer, Variable<double> >())
		   .def(VariableIndexingPython<FixDataValueContainer, Variable<array_1d<double, 3> > >())
		   .def(VariableIndexingPython<FixDataValueContainer, Variable<vector<double> > >())
		   .def(VariableIndexingPython<FixDataValueContainer, Variable<matrix<double> > >())
		   .def(VariableIndexingPython<FixDataValueContainer, 			VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > > >())
	       .def(self_ns::str(self))
      ;*/

       class_<VariablesListDataValueContainer, VariablesListDataValueContainer::Pointer>("VariablesListDataValueContainer")
		   .def("__len__", &VariablesListDataValueContainer::Size)
		   .def(VariableIndexingPython<VariablesListDataValueContainer, Variable<int> >())
		   .def(VariableIndexingPython<VariablesListDataValueContainer, Variable<double> >())
		   .def(VariableIndexingPython<VariablesListDataValueContainer, Variable<array_1d<double, 3> > >())
		   .def(VariableIndexingPython<VariablesListDataValueContainer, Variable<vector<double> > >())
		   .def(VariableIndexingPython<VariablesListDataValueContainer, Variable<matrix<double> > >())
		   .def(VariableIndexingPython<VariablesListDataValueContainer, 			VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > > >())
	       .def(self_ns::str(self))
      ;

	KRATOS_REGISTER_IN_PYTHON_VARIABLE( OSS_SWITCH)
        KRATOS_REGISTER_IN_PYTHON_VARIABLE(DYNAMIC_TAU)

	KRATOS_REGISTER_IN_PYTHON_VARIABLE(ERASE_FLAG)
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(NL_ITERATION_NUMBER)
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(FRACTIONAL_STEP)
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(TIME)
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(TIME_STEPS)
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(START_TIME)
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(END_TIME)
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(DELTA_TIME)
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(TEMPERATURE)
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(PRESSURE)
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(DENSITY)
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(VISCOSITY)
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(VISCOSITY_AIR)
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(VISCOSITY_WATER)
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(ERROR_RATIO)
        KRATOS_REGISTER_IN_PYTHON_VARIABLE(TAU)


	KRATOS_REGISTER_IN_PYTHON_VARIABLE( FC )
	KRATOS_REGISTER_IN_PYTHON_VARIABLE( FT )
	KRATOS_REGISTER_IN_PYTHON_VARIABLE( FRACTURE_ENERGY ) 
	KRATOS_REGISTER_IN_PYTHON_VARIABLE( CONCRETE_YOUNG_MODULUS_C)
	KRATOS_REGISTER_IN_PYTHON_VARIABLE( CONCRETE_YOUNG_MODULUS_T)
	KRATOS_REGISTER_IN_PYTHON_VARIABLE( BASE )
	KRATOS_REGISTER_IN_PYTHON_VARIABLE( HEIGHT)
	KRATOS_REGISTER_IN_PYTHON_VARIABLE( YIELD_STRESS)
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(LAMNDA)   
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(DAMAGE )
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(VECTOR_DAMAGE)
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(ORTHOTROPIC_YOUNG_MODULUS_2D) // [E1 E2 G12]
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(ORTHOTROPIC_POISSON_RATIO_2D) // [v12 v21]
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(ORTHOTROPIC_ANGLE)
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(VOLUMEN_FRACTION)
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(FRICTION_INTERNAL_ANGLE)
        KRATOS_REGISTER_IN_PYTHON_VARIABLE(COHESION)
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(ORTHOTROPIC_ELASTIC_LIMIT)
        KRATOS_REGISTER_IN_PYTHON_VARIABLE(ISOTROPIC_ELASTIC_LIMIT)

	
	KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(ANGULAR_ACCELERATION)
	KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(ANGULAR_VELOCITY)
	KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(ACCELERATION)
  	KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(VELOCITY)
  	KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(DISPLACEMENT)
	KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(ROTATION)
	KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(MOMENTUM)
  	KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(REACTION)
  	KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(BODY_FORCE)
	KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(SEEPAGE_DRAG)
	KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(MESH_VELOCITY)
	KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(NORMAL)
	KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(FORCE)
        KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(FORCE_CM)
        KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(MOMENTUM_CM)

	KRATOS_REGISTER_IN_PYTHON_VARIABLE( WATER_PRESSURE )
	KRATOS_REGISTER_IN_PYTHON_VARIABLE( REACTION_WATER_PRESSURE )
	KRATOS_REGISTER_IN_PYTHON_VARIABLE( WATER_PRESSURE_DT )
	KRATOS_REGISTER_IN_PYTHON_VARIABLE( WATER_PRESSURE_ACCELERATION )
			
	KRATOS_REGISTER_IN_PYTHON_VARIABLE( AIR_PRESSURE )
	KRATOS_REGISTER_IN_PYTHON_VARIABLE( REACTION_AIR_PRESSURE )
	KRATOS_REGISTER_IN_PYTHON_VARIABLE( AIR_PRESSURE_DT )
	KRATOS_REGISTER_IN_PYTHON_VARIABLE( AIR_PRESSURE_ACCELERATION )
	//for structural application TO BE REMOVED
        KRATOS_REGISTER_IN_PYTHON_VARIABLE( INTERNAL_VARIABLES )
        KRATOS_REGISTER_IN_PYTHON_VARIABLE( MATERIAL_PARAMETERS )
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(WRINKLING_APPROACH )
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(GREEN_LAGRANGE_STRAIN_TENSOR )
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(PK2_STRESS_TENSOR )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(CAUCHY_STRESS_TENSOR)
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(AUXILIARY_MATRIX_1 )
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(YOUNG_MODULUS )
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(POISSON_RATIO )
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(THICKNESS )
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(NEGATIVE_FACE_PRESSURE )
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(POSITIVE_FACE_PRESSURE )
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(IS_POROUS)	
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(IS_WATER)	
	KRATOS_REGISTER_IN_PYTHON_VARIABLE( POROSITY)


	KRATOS_REGISTER_IN_PYTHON_VARIABLE(FLAG_VARIABLE )
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(DISTANCE )
			
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(INERTIA )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(PARTITION_INDEX )
            
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( LAGRANGE_DISPLACEMENT )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( LAGRANGE_AIR_PRESSURE )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( LAGRANGE_WATER_PRESSURE )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( LAGRANGE_TEMPERATURE )

	//for ALE application TO BE REMOVED
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(IS_INTERFACE )
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(IS_VISITED )

        
    //for electric application
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( ELECTRIC_POTENTIAL )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(ELECTRIC_POTENTIAL_DT )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(ELECTRIC_POTENTIAL_DT2 )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(ELECTRIC_POTENTIAL_NULL )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(ELECTRIC_POTENTIAL_EINS )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(ELECTRIC_POTENTIAL_NULL_DT )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(ELECTRIC_POTENTIAL_EINS_DT )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(ELECTRIC_POTENTIAL_NULL_DT2 )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(ELECTRIC_POTENTIAL_EINS_DT2 )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(LAGRANGE_ELECTRIC_POTENTIAL)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(CONCENTRATION )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(CONCENTRATION_DT )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(CONCENTRATION_DT2 )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(CONCENTRATION_NULL )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(CONCENTRATION_EINS )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(CONCENTRATION_NULL_DT )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(CONCENTRATION_EINS_DT )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(CONCENTRATION_NULL_DT2 )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(CONCENTRATION_EINS_DT2 )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(LAGRANGE_CONCENTRATION)
	
	//for PFEM application TO BE REMOVED
	KRATOS_REGISTER_IN_PYTHON_VARIABLE( NODAL_AREA)
	KRATOS_REGISTER_IN_PYTHON_VARIABLE( NODAL_H)
	KRATOS_REGISTER_IN_PYTHON_VARIABLE( IS_STRUCTURE)
	KRATOS_REGISTER_IN_PYTHON_VARIABLE( IS_FLUID)
	KRATOS_REGISTER_IN_PYTHON_VARIABLE( IS_BOUNDARY)
	KRATOS_REGISTER_IN_PYTHON_VARIABLE( IS_FREE_SURFACE)
	KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(NORMAL_TO_WALL)
	//KRATOS_REGISTER_IN_PYTHON_VARIABLE(NEIGHBOUR_NODES);
	//KRATOS_REGISTER_IN_PYTHON_VARIABLE(NEIGHBOUR_ELEMENTS);
   	KRATOS_REGISTER_IN_PYTHON_VARIABLE( IS_EROSIONABLE)
 	KRATOS_REGISTER_IN_PYTHON_VARIABLE(FRICTION_COEFFICIENT)

	KRATOS_REGISTER_IN_PYTHON_VARIABLE( BULK_MODULUS )
	
	KRATOS_REGISTER_IN_PYTHON_VARIABLE( SATURATION )
	KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( GRAVITY )
	KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( FACE_LOAD )
	KRATOS_REGISTER_IN_PYTHON_VARIABLE( DENSITY_WATER )
	KRATOS_REGISTER_IN_PYTHON_VARIABLE( DENSITY_AIR )
	KRATOS_REGISTER_IN_PYTHON_VARIABLE( POROSITY )
	KRATOS_REGISTER_IN_PYTHON_VARIABLE( AIR_ENTRY_VALUE )
	KRATOS_REGISTER_IN_PYTHON_VARIABLE( FIRST_SATURATION_PARAM )
	KRATOS_REGISTER_IN_PYTHON_VARIABLE( SECOND_SATURATION_PARAM )
	KRATOS_REGISTER_IN_PYTHON_VARIABLE( PERMEABILITY_WATER )
	KRATOS_REGISTER_IN_PYTHON_VARIABLE( PERMEABILITY_AIR )
	KRATOS_REGISTER_IN_PYTHON_VARIABLE( BULK_AIR )
	KRATOS_REGISTER_IN_PYTHON_VARIABLE( SCALE )
	KRATOS_REGISTER_IN_PYTHON_VARIABLE( INSITU_STRESS )
	KRATOS_REGISTER_IN_PYTHON_VARIABLE( ELASTIC_LEFT_CAUCHY_GREEN_OLD )

        KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(FRACT_VEL)
	KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(PRESS_PROJ)
	KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(CONV_PROJ)
         KRATOS_REGISTER_IN_PYTHON_VARIABLE( NODAL_MASS)
	 KRATOS_REGISTER_IN_PYTHON_VARIABLE(  AUX_INDEX)
	 KRATOS_REGISTER_IN_PYTHON_VARIABLE(  EXTERNAL_PRESSURE)
         KRATOS_REGISTER_IN_PYTHON_VARIABLE( PRESSURE_OLD_IT )
	KRATOS_REGISTER_IN_PYTHON_VARIABLE( BDF_COEFFICIENTS )


	// for mengmeng application
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(TEMPERATURE_NULL )
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(TEMPERATURE_EINS )
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(SURFACE_FLOW_HEAT )
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(SURFACE_FLOW_WATER )
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(ICE_VOLUME_FRACTION )
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(ICE_VOLUME_FRACTION_NULL )
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(ICE_VOLUME_FRACTION_EINS )

        //for xfem application
        KRATOS_REGISTER_IN_PYTHON_VARIABLE( CRACK_OPENING )
        KRATOS_REGISTER_IN_PYTHON_VARIABLE( CRACK_TRANSLATION )

	KRATOS_REGISTER_IN_PYTHON_VARIABLE(ARRHENIUS )
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(ARRHENIUSAUX )
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(PRESSUREAUX)
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(NODAL_MAUX)
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(NODAL_PAUX)


   //scope().attr("TIME") = boost::ref(TIME);
  }
  
}  // namespace Python.

} // Namespace Kratos

#undef KRATOS_REGISTER_IN_PYTHON_VARIABLE
#undef KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS

