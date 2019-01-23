//
//   Project Name:        KratosSolversApplication $
//   Created by:          $Author:     JMCarbonell $
//   Last modified by:    $Co-Author:              $
//   Date:                $Date:      January 2019 $
//   Revision:            $Revision:           0.0 $
//
//

// System includes

// External includes

// Project includes

#include "solvers_application_variables.h"

namespace Kratos
{
  ///@name Type Definitions
  ///@{
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

  // Create Variables:

  // time settings
  KRATOS_CREATE_VARIABLE( double, MESHING_STEP_TIME )
  KRATOS_CREATE_VARIABLE( double, CONTACT_STEP_TIME )
  KRATOS_CREATE_VARIABLE( double, RESTART_STEP_TIME )

  // time integration methods
  KRATOS_CREATE_VARIABLE( VectorTimeIntegrationContainerPointerType, VECTOR_TIME_INTEGRATION_METHODS )
  KRATOS_CREATE_VARIABLE( ScalarTimeIntegrationContainerPointerType, SCALAR_TIME_INTEGRATION_METHODS )
  KRATOS_CREATE_VARIABLE( ComponentTimeIntegrationContainerPointerType, COMPONENT_TIME_INTEGRATION_METHODS )

  // implicit solution
  KRATOS_CREATE_VARIABLE( bool, CONVERGENCE_ACHIEVED)
  KRATOS_CREATE_VARIABLE( bool, COMPUTE_CONSISTENT_MASS_MATRIX)

  KRATOS_CREATE_VARIABLE( int, SEGREGATED_STEP )
  KRATOS_CREATE_VARIABLE( int, TIME_INTEGRATION_ORDER )

  KRATOS_CREATE_VARIABLE( double, RAYLEIGH_ALPHA )
  KRATOS_CREATE_VARIABLE( double, RAYLEIGH_BETA )

  // explicit solution
  KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( MIDDLE_VELOCITY )

  KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( EXTERNAL_MOMENT )
  KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( POSITION_MOMENTUM )
  KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( ROTATION_MOMENTUM )
  KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( RESIDUAL_LYAPUNOV )

  KRATOS_CREATE_VARIABLE( Matrix, INERTIA_DYADIC )
  KRATOS_CREATE_VARIABLE( Matrix, TANGENT_MATRIX )
  KRATOS_CREATE_VARIABLE( Matrix, TANGENT_LYAPUNOV )

  KRATOS_CREATE_VARIABLE( double, ALPHA_TRAPEZOIDAL_RULE )
  KRATOS_CREATE_VARIABLE( bool, POSITION_UPDATE_LABEL )
  KRATOS_CREATE_VARIABLE( bool, ROTATION_UPDATE_LABEL )
  KRATOS_CREATE_VARIABLE( bool, MOMENTUM_UPDATE_LABEL )

  // eigenvalue solution
  KRATOS_CREATE_VARIABLE( int, BUILD_LEVEL )
  KRATOS_CREATE_VARIABLE( Vector, EIGENVALUE_VECTOR)
  KRATOS_CREATE_VARIABLE( Matrix, EIGENVECTOR_MATRIX )

  ///@}

} // Namespace Kratos
