//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 license: OptimizationApplication/license.txt
//
//  Main authors:    Reza Najian Asl, https://github.com/RezaNajian
//

#ifndef RESPONSE_H
#define RESPONSE_H

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------
#include <iostream>
#include <string>

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "containers/model.h"
#include "includes/model_part.h"
#include "optimization_application.h"
#include "utilities/integration_utilities.h"
#include "utilities/geometry_utilities.h"
#include "utilities/variable_utils.h"
#include "spaces/ublas_space.h"
#include "linear_solvers/linear_solver.h"
#include "solving_strategies/strategies/residualbased_linear_strategy.h"
#include "includes/define.h"
#include "utilities/math_utils.h"
#include "spatial_containers/spatial_containers.h"
#include "processes/find_conditions_neighbours_process.h"
#include "custom_strategies/strategies/helmholtz_strategy.h"
// ==============================================================================

namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/// Short class definition.
/** Detail class definition.

*/

class Response
{
public:
  ///@name Type Definitions
  ///@{

    // Type definitions for better reading later
    typedef Variable<double> array_1d_component_type;  
    typedef array_1d<double,3> array_3d;
    typedef Element BaseType;
    typedef BaseType::GeometryType GeometryType;
    typedef BaseType::NodesArrayType NodesArrayType;
    typedef BaseType::PropertiesType PropertiesType;
    typedef BaseType::IndexType IndexType;
    typedef BaseType::SizeType SizeType;    
    typedef BaseType::MatrixType MatrixType;
    typedef BaseType::VectorType VectorType;    
    typedef GeometryData::IntegrationMethod IntegrationMethod;
    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
    typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;  
    typedef HelmholtzStrategy<SparseSpaceType, LocalSpaceType,LinearSolverType> StrategyType;

  /// Pointer definition of Response
  KRATOS_CLASS_POINTER_DEFINITION(Response);

  ///@}
  ///@name Life Cycle
  ///@{

  /// Default constructor.
  Response(std::string ResponseName, std::string ResponseType, Model& rModel, Parameters& rResponseSettings): mResponseName(ResponseName), mResponseType(ResponseType), mrModel(rModel), mrResponseSettings(rResponseSettings){}

  /// Destructor.
  virtual ~Response(){}

  ///@}
  ///@name Operators
  ///@{

  ///@}
  ///@name Operations
  ///@{

  // --------------------------------------------------------------------------
  virtual void Initialize(){};

  // --------------------------------------------------------------------------
  virtual double CalculateValue() = 0;

  // --------------------------------------------------------------------------
  virtual void CalculateGradient() {};

  // --------------------------------------------------------------------------

  ///@}
  ///@name Access
  ///@{

  ///@}
  ///@name Inquiry
  ///@{

  ///@}
  ///@name Input and output
  ///@{

  /// Turn back information as a string.
  virtual std::string Info() const
  {
    return "Response base class";
  }

  /// Print information about this object.
  virtual void PrintInfo(std::ostream &rOStream) const
  {
    rOStream << "Response base class";
  }

  /// Print object's data.
  virtual void PrintData(std::ostream &rOStream) const
  {
  }

  ///@}
  ///@name Friends
  ///@{

  ///@}

protected:
  ///@name Protected static Member Variables
  ///@{

  ///@}
  ///@name Protected member Variables
  ///@{
    std::string mResponseName;
    std::string mResponseType;
    Model& mrModel;
    Parameters& mrResponseSettings;

  ///@}
  ///@name Protected Operators
  ///@{

  ///@}
  ///@name Protected Operations
  ///@{

  ///@}
  ///@name Protected  Access
  ///@{

  ///@}
  ///@name Protected Inquiry
  ///@{

  ///@}
  ///@name Protected LifeCycle
  ///@{

  ///@}

private:
  ///@name Static Member Variables
  ///@{

  ///@}
  ///@name Member Variables
  ///@{

  ///@}
  ///@name Private Operators
  ///@{

  ///@}
  ///@name Private Operations
  ///@{

  ///@}
  ///@name Private  Access
  ///@{

  ///@}
  ///@name Private Inquiry
  ///@{

  ///@}
  ///@name Un accessible methods
  ///@{

  /// Assignment operator.
  //      RESPONSE& operator=(RESPONSE const& rOther);

  /// Copy constructor.
  //      RESPONSE(RESPONSE const& rOther);

  ///@}

}; // Class RESPONSE

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

} // namespace Kratos.

#endif // RESPONSE_H
