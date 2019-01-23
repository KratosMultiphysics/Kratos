//------------------------------------------------------------------
//           ___      _                                            .
//   KRATOS / __| ___| |_ _ ___ _ _ ___                            .
//          \__ \/ _ \ \ V / -_) '_|_-<                            .
//          |___/\___/_|\_/\___|_| /__/ APPLICATION                .
//                                                                 .
//   License:(BSD)	  SolversApplication/license.txt           .
//   Main authors:        Josep Maria Carbonell                    .
//                        ..                                       .
//------------------------------------------------------------------
//
//   Project Name:        KratosSolversApplication $
//   Created by:          $Author:     JMCarbonell $
//   Last modified by:    $Co-Author:              $
//   Date:                $Date:      January 2019 $
//   Revision:            $Revision:           0.0 $
//
//

#if !defined(KRATOS_SOLVERS_APPLICATION_VARIABLES_H_INCLUDED )
#define  KRATOS_SOLVERS_APPLICATION_VARIABLES_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/variables.h"
#include "includes/mat_variables.h"
#include "includes/kratos_application.h"
#include "includes/kratos_flags.h"
#include "custom_solvers/time_integration_methods/time_integration_methods_container.hpp"

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

  // Define Variables:

  // time settings
  KRATOS_DEFINE_APPLICATION_VARIABLE( SOLVERS_APPLICATION, double, MESHING_STEP_TIME )
  KRATOS_DEFINE_APPLICATION_VARIABLE( SOLVERS_APPLICATION, double, CONTACT_STEP_TIME )
  KRATOS_DEFINE_APPLICATION_VARIABLE( SOLVERS_APPLICATION, double, RESTART_STEP_TIME )

  // time integration methods
  KRATOS_DEFINE_APPLICATION_VARIABLE( SOLVERS_APPLICATION, VectorTimeIntegrationContainerPointerType, VECTOR_TIME_INTEGRATION_METHODS )
  KRATOS_DEFINE_APPLICATION_VARIABLE( SOLVERS_APPLICATION, ScalarTimeIntegrationContainerPointerType, SCALAR_TIME_INTEGRATION_METHODS )
  KRATOS_DEFINE_APPLICATION_VARIABLE( SOLVERS_APPLICATION, ComponentTimeIntegrationContainerPointerType, COMPONENT_TIME_INTEGRATION_METHODS )

  // implicit solution
  KRATOS_DEFINE_APPLICATION_VARIABLE( SOLVERS_APPLICATION, bool, CONVERGENCE_ACHIEVED )
  KRATOS_DEFINE_APPLICATION_VARIABLE( SOLVERS_APPLICATION, bool, COMPUTE_CONSISTENT_MASS_MATRIX )

  KRATOS_DEFINE_APPLICATION_VARIABLE( SOLVERS_APPLICATION, int, SEGREGATED_STEP )
  KRATOS_DEFINE_APPLICATION_VARIABLE( SOLVERS_APPLICATION, int, TIME_INTEGRATION_ORDER )

  KRATOS_DEFINE_APPLICATION_VARIABLE( SOLVERS_APPLICATION, double, RAYLEIGH_ALPHA )
  KRATOS_DEFINE_APPLICATION_VARIABLE( SOLVERS_APPLICATION, double, RAYLEIGH_BETA )

  // explicit solution
  KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS( SOLVERS_APPLICATION, MIDDLE_VELOCITY )

  KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS( SOLVERS_APPLICATION, EXTERNAL_MOMENT )
  KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS( SOLVERS_APPLICATION, POSITION_MOMENTUM )
  KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS( SOLVERS_APPLICATION, ROTATION_MOMENTUM )
  KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS( SOLVERS_APPLICATION, RESIDUAL_LYAPUNOV )

  KRATOS_DEFINE_APPLICATION_VARIABLE( SOLVERS_APPLICATION, Matrix, INERTIA_DYADIC )
  KRATOS_DEFINE_APPLICATION_VARIABLE( SOLVERS_APPLICATION, Matrix, TANGENT_MATRIX )
  KRATOS_DEFINE_APPLICATION_VARIABLE( SOLVERS_APPLICATION, Matrix, TANGENT_LYAPUNOV )

  KRATOS_DEFINE_APPLICATION_VARIABLE( SOLVERS_APPLICATION, double, ALPHA_TRAPEZOIDAL_RULE )
  KRATOS_DEFINE_APPLICATION_VARIABLE( SOLVERS_APPLICATION, bool, POSITION_UPDATE_LABEL )
  KRATOS_DEFINE_APPLICATION_VARIABLE( SOLVERS_APPLICATION, bool, ROTATION_UPDATE_LABEL )
  KRATOS_DEFINE_APPLICATION_VARIABLE( SOLVERS_APPLICATION, bool, MOMENTUM_UPDATE_LABEL )

  // eigenvalue solution
  KRATOS_DEFINE_APPLICATION_VARIABLE( SOLVERS_APPLICATION, int, BUILD_LEVEL )
  KRATOS_DEFINE_APPLICATION_VARIABLE( SOLVERS_APPLICATION, Vector, EIGENVALUE_VECTOR)
  KRATOS_DEFINE_APPLICATION_VARIABLE( SOLVERS_APPLICATION, Matrix, EIGENVECTOR_MATRIX )


  ///@}

} // Namespace Kratos

#endif	//KRATOS_SOLVERS_APPLICATION_VARIABLES_H_INCLUDED
