//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher

#if !defined(KRATOS_MAPPER_H_INCLUDED )
#define  KRATOS_MAPPER_H_INCLUDED

// System includes
#include <string>
#include <iostream>
#include <algorithm>
#include <vector>

// External includes
//@{KRATOS_EXTERNA_INCLUDES}
#include "includes/kratos_flags.h"

// Project includes

#include "includes/define.h"
#include "includes/kratos_flags.h"
#include "includes/element.h"
#include "includes/model_part.h"
#include "geometries/geometry_data.h"

#include "spaces/ublas_space.h"
#include "linear_solvers/linear_solver.h"
// #include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"
// #include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"
// #include "solving_strategies/strategies/residualbased_linear_strategy.h"
// #include "elements/distance_calculation_element_simplex.h"

#include "spatial_containers/bins_dynamic_objects.h"
#include "interface_object.h"
#include "custom_configures/interface_object_configure.h"
#include "interface_object_manager.h"
#include "mapper_communicator.h"

#include <omp.h>
#include "utilities/timer.h"

// For MPI-parallel Mapper
#ifdef KRATOS_USING_MPI
#include "mpi.h"
#include "mapper_mpi_communicator.h"
#endif

namespace Kratos {

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

class Mapper {
public:

  ///@name Type Definitions
  ///@{
  ///@}
  ///@name Pointer Definitions
  /// Pointer definition of Mapper
  KRATOS_CLASS_POINTER_DEFINITION(Mapper);

  typedef matrix<int> GraphType; // GraphColoringProcess

  ///@}
  ///@name Life Cycle
  ///@{

  Mapper(ModelPart& i_model_part_origin, ModelPart& i_model_part_destination) :
	  m_model_part_origin(i_model_part_origin), m_model_part_destination(i_model_part_destination) {

          m_num_conditions_origin = m_model_part_origin.GetCommunicator().LocalMesh().Elements().size();
          m_num_conditions_destination = m_model_part_destination.GetCommunicator().LocalMesh().Elements().size();

          m_num_nodes_origin = m_model_part_origin.GetCommunicator().LocalMesh().Nodes().size();
          m_num_nodes_destination = m_model_part_destination.GetCommunicator().LocalMesh().Nodes().size();

          #ifdef KRATOS_USING_MPI // mpi-parallel compilation
              int comm_size;
              MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
              if (comm_size > 1) { // parallel execution
                  m_mapper_communicator = MapperCommunicator::Pointer ( new MapperMPICommunicator(m_model_part_origin, m_model_part_destination) );
              } else { // serial execution
                  m_mapper_communicator = MapperCommunicator::Pointer ( new MapperCommunicator(m_model_part_origin, m_model_part_destination) );
              }
          #else // serial compilation
              m_mapper_communicator = MapperCommunicator::Pointer ( new MapperCommunicator(m_model_part_origin, m_model_part_destination) );
          #endif
  }

  /// Destructor.
  virtual ~Mapper() {
  }

  ///@}
  ///@name Operators
  ///@{

  ///@}
  ///@name Operations
  ///@{

  virtual void UpdateInterface() = 0;

  /* This function maps from Origin to Destination */
  virtual void Map(const Variable<double>& origin_variable,
                   const Variable<double>& destination_variable,
                   const bool add_value) = 0;

  /* This function maps from Origin to Destination */
  virtual void Map(const Variable< array_1d<double,3> >& origin_variable,
                   const Variable< array_1d<double,3> >& destination_variable,
                   const bool add_value) = 0;

  /* This function maps from Destination to Origin */
  virtual void InverseMap(const Variable<double>& origin_variable,
                          const Variable<double>& destination_variable,
                          const bool add_value) = 0;

  /* This function maps from Destination to Origin */
  virtual void InverseMap(const Variable< array_1d<double,3> >& origin_variable,
                          const Variable< array_1d<double,3> >& destination_variable,
                          const bool add_value) = 0;

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
  virtual std::string Info() const {
      return "Mapper";
  }

  /// Print information about this object.
  virtual void PrintInfo(std::ostream& rOStream) const {
    rOStream << "Mapper";
  }

  /// Print object's data.
  virtual void PrintData(std::ostream& rOStream) const {
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
      ModelPart& m_model_part_origin;
      ModelPart& m_model_part_destination;

      MapperCommunicator::Pointer m_mapper_communicator;

      InterfaceObjectManager::Pointer m_point_comm_manager_origin;
      InterfaceObjectManager::Pointer m_point_comm_manager_destination;

      int m_num_conditions_origin;
      int m_num_conditions_destination;

      int m_num_nodes_origin;
      int m_num_nodes_destination;

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
  Mapper& operator=(Mapper const& rOther);

  /// Copy constructor.
  //Mapper(Mapper const& rOther);

  ///@}

}; // Class Mapper

}  // namespace Kratos.

#endif // KRATOS_MAPPER_H_INCLUDED  defined
