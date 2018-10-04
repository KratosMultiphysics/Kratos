//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:            February 2016 $
//   Revision:            $Revision:                  0.0 $
//
//

#include "solid_mechanics_application_variables.h"

namespace Kratos
{
  ///@name Type Definitions
  ///@{
  typedef array_1d<double,3> Vector3;
  typedef array_1d<double,6> Vector6;

  typedef array_1d<double, 3>                                                                      VectorType;
  typedef Variable<VectorType>                                                             VariableVectorType;
  typedef Variable<double>                                                                 VariableScalarType;
  typedef VariableComponent<VectorComponentAdaptor<VectorType>>                         VariableComponentType;

  typedef TimeIntegrationMethodsContainer<VariableVectorType, double>      VectorTimeIntegrationContainerType;
  typedef VectorTimeIntegrationContainerType::Pointer               VectorTimeIntegrationContainerPointerType;

  typedef TimeIntegrationMethodsContainer<VariableScalarType, double>      ScalarTimeIntegrationContainerType;
  typedef ScalarTimeIntegrationContainerType::Pointer               ScalarTimeIntegrationContainerPointerType;

  typedef TimeIntegrationMethodsContainer<VariableComponentType, double> ComponentTimeIntegrationContainerType;
  typedef ComponentTimeIntegrationContainerType::Pointer          ComponentTimeIntegrationContainerPointerType;

  ///@}

  ///@name Kratos Globals
  ///@{

  //Create Variables

  // Generalized eigenvalue problem
  KRATOS_CREATE_VARIABLE( int, BUILD_LEVEL )
  KRATOS_CREATE_VARIABLE( Vector, EIGENVALUE_VECTOR)
  KRATOS_CREATE_VARIABLE( Matrix , EIGENVECTOR_MATRIX )

  //integration methods
  KRATOS_CREATE_VARIABLE( VectorTimeIntegrationContainerPointerType, VECTOR_TIME_INTEGRATION_METHODS )
  KRATOS_CREATE_VARIABLE( ScalarTimeIntegrationContainerPointerType, SCALAR_TIME_INTEGRATION_METHODS )
  KRATOS_CREATE_VARIABLE( ComponentTimeIntegrationContainerPointerType, COMPONENT_TIME_INTEGRATION_METHODS )

  //explicit schemes
  KRATOS_CREATE_VARIABLE(bool, COMPUTE_CONSISTENT_MASS_MATRIX)
  KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( MIDDLE_VELOCITY )

  //variables
  KRATOS_CREATE_VARIABLE( double, PRESSURE_VELOCITY )
  KRATOS_CREATE_VARIABLE( double, PRESSURE_ACCELERATION )

  //solution
  KRATOS_CREATE_VARIABLE( bool, DELTA_TIME_CHANGED)
  KRATOS_CREATE_VARIABLE( bool, CONVERGENCE_ACHIEVED)
  KRATOS_CREATE_VARIABLE( int, SEGREGATED_STEP )
  KRATOS_CREATE_VARIABLE( int, WRITE_ID )
  KRATOS_CREATE_VARIABLE( int, TIME_INTEGRATION_ORDER )
  KRATOS_CREATE_VARIABLE( double, RAYLEIGH_ALPHA )
  KRATOS_CREATE_VARIABLE( double, RAYLEIGH_BETA )
  KRATOS_CREATE_VARIABLE( double, MESHING_STEP_TIME )
  KRATOS_CREATE_VARIABLE( double, CONTACT_STEP_TIME )
  KRATOS_CREATE_VARIABLE( double, RESTART_STEP_TIME )

  //geometrical
  KRATOS_CREATE_VARIABLE( Matrix ,GEOMETRIC_STIFFNESS )

  //beam cross section
  //KRATOS_CREATE_VARIABLE( BeamCrossSection::Pointer, BEAM_CROSS_SECTION )
  KRATOS_CREATE_VARIABLE( double, CROSS_SECTION_AREA )
  KRATOS_CREATE_VARIABLE( double, CROSS_SECTION_RADIUS )
  KRATOS_CREATE_VARIABLE( int,    CROSS_SECTION_SIDES )

  //shell cross section
  KRATOS_CREATE_VARIABLE( ShellCrossSection::Pointer, SHELL_CROSS_SECTION )
  KRATOS_CREATE_VARIABLE( int, SHELL_CROSS_SECTION_OUTPUT_PLY_ID )
  KRATOS_CREATE_VARIABLE( double, SHELL_CROSS_SECTION_OUTPUT_PLY_LOCATION )

  //shell generalized variables
  KRATOS_CREATE_VARIABLE( Matrix, SHELL_STRAIN )
  KRATOS_CREATE_VARIABLE( Matrix, SHELL_STRAIN_GLOBAL )
  KRATOS_CREATE_VARIABLE( Matrix, SHELL_CURVATURE )
  KRATOS_CREATE_VARIABLE( Matrix, SHELL_CURVATURE_GLOBAL )
  KRATOS_CREATE_VARIABLE( Matrix, SHELL_FORCE )
  KRATOS_CREATE_VARIABLE( Matrix, SHELL_FORCE_GLOBAL )
  KRATOS_CREATE_VARIABLE( Matrix, SHELL_MOMENT )
  KRATOS_CREATE_VARIABLE( Matrix, SHELL_MOMENT_GLOBAL )

  //nodal load variable (legacy)
  KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( POINT_LOAD )

  //force loads
  KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( FORCE_LOAD )
  KRATOS_CREATE_VARIABLE( Vector, FORCE_LOAD_VECTOR )

  KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( FOLLOWER_FORCE_LOAD )
  KRATOS_CREATE_VARIABLE( Vector, FOLLOWER_FORCE_LOAD_VECTOR )

  //moment loads
  KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( MOMENT_LOAD )
  KRATOS_CREATE_VARIABLE( Vector, MOMENT_LOAD_VECTOR )

  //elastic loads
  KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( ELASTIC_LOAD )
  KRATOS_CREATE_VARIABLE( Vector, ELASTIC_LOAD_VECTOR )

  //force pressure
  KRATOS_CREATE_VARIABLE( Vector, POSITIVE_FACE_PRESSURE_VECTOR )
  KRATOS_CREATE_VARIABLE( Vector, NEGATIVE_FACE_PRESSURE_VECTOR )

  //moment pressures
  KRATOS_CREATE_VARIABLE( double, PLANE_MOMENT_LOAD )
  KRATOS_CREATE_VARIABLE( Vector, PLANE_MOMENT_LOAD_VECTOR )

  //elastic pressures
  KRATOS_CREATE_VARIABLE( double, BALLAST_COEFFICIENT )
  KRATOS_CREATE_VARIABLE( Vector, BALLAST_COEFFICIENT_VECTOR )

  //heat fluxes
  KRATOS_CREATE_VARIABLE( Vector, FACE_HEAT_FLUX_VECTOR )

  //element
  KRATOS_CREATE_VARIABLE( double, VON_MISES_STRESS )

  //nodal dofs
  KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( DISPLACEMENT_REACTION )
  KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( ROTATION_REACTION )
  KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( VELOCITY_REACTION )
  KRATOS_CREATE_VARIABLE( double, PRESSURE_REACTION )
  KRATOS_CREATE_VARIABLE( double, TEMPERATURE_REACTION )

  //explicit beam
  KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( EXTERNAL_MOMENT )

  KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( POSITION_MOMENTUM )
  KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( ROTATION_MOMENTUM )

  KRATOS_CREATE_VARIABLE( Matrix, INERTIA_DYADIC )

  KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( RESIDUAL_LYAPUNOV )
  KRATOS_CREATE_VARIABLE( Matrix, TANGENT_MATRIX )
  KRATOS_CREATE_VARIABLE( Matrix, TANGENT_LYAPUNOV )

  KRATOS_CREATE_VARIABLE( double, ALPHA_TRAPEZOIDAL_RULE )
  KRATOS_CREATE_VARIABLE( bool, POSITION_UPDATE_LABEL )
  KRATOS_CREATE_VARIABLE( bool, ROTATION_UPDATE_LABEL )
  KRATOS_CREATE_VARIABLE( bool, MOMENTUM_UPDATE_LABEL )

  //reading beam section properties
  KRATOS_CREATE_VARIABLE( double, SECTION_HEIGHT )
  KRATOS_CREATE_VARIABLE( double, SECTION_WIDTH  )
  KRATOS_CREATE_VARIABLE( double, INERTIA_X )
  KRATOS_CREATE_VARIABLE( double, INERTIA_Y )
  KRATOS_CREATE_VARIABLE( double, SECTION_SIZE )

  KRATOS_CREATE_VARIABLE( double, YOUNGxAREA )
  KRATOS_CREATE_VARIABLE( double, YOUNGxINERTIA_X )
  KRATOS_CREATE_VARIABLE( double, YOUNGxINERTIA_Y )
  KRATOS_CREATE_VARIABLE( double, SHEARxREDUCED_AREA )
  KRATOS_CREATE_VARIABLE( double, SHEARxPOLAR_INERTIA )

  //boundary definition
  KRATOS_CREATE_VARIABLE( WeakPointerVector< Element >, MASTER_ELEMENTS )

  //thermal properties
  KRATOS_CREATE_VARIABLE( double, HEAT_CAPACITY )
  KRATOS_CREATE_VARIABLE( double, HEAT_CONDUCTIVITY )
  KRATOS_CREATE_VARIABLE( double, HEAT_SOURCE )

  KRATOS_CREATE_VARIABLE( bool,   TEMPERATURE_DEPENDENT )
  KRATOS_CREATE_VARIABLE( double, HEAT_CAPACITY_A )
  KRATOS_CREATE_VARIABLE( double, HEAT_CAPACITY_B )
  KRATOS_CREATE_VARIABLE( double, HEAT_CONDUCTIVITY_A )
  KRATOS_CREATE_VARIABLE( double, HEAT_CONDUCTIVITY_B )

  ///@}

}
