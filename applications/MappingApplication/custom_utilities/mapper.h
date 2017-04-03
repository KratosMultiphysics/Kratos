//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher, Jordi Cotela
//
// See Master-Thesis P.Bucher
// "Development and Implementation of a Parallel
//  Framework for Non-Matching Grid Mapping"

#if !defined(KRATOS_MAPPER_H_INCLUDED )
#define  KRATOS_MAPPER_H_INCLUDED

// System includes
#include <string>
#include <iostream>
#include <vector>

// External includes

// Project includes

#include "includes/define.h"
#include "mapper_communicator.h"
#include "mapper_utilities.h"
#include "mapper_flags.h"

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

  ///@}
  ///@name Life Cycle
  ///@{


  /// Destructor.
  virtual ~Mapper() {
  }

  ///@}
  ///@name Operators
  ///@{

  ///@}
  ///@name Operations
  ///@{

  virtual void UpdateInterface(Kratos::Flags& options, double search_radius) = 0;

  /* This function maps from Origin to Destination */
  virtual void Map(const Variable<double>& origin_variable,
                   const Variable<double>& destination_variable,
                   Kratos::Flags& options) = 0;

  /* This function maps from Origin to Destination */
  virtual void Map(const Variable< array_1d<double,3> >& origin_variable,
                   const Variable< array_1d<double,3> >& destination_variable,
                   Kratos::Flags& options) = 0;

  /* This function maps from Destination to Origin */
  virtual void InverseMap(const Variable<double>& origin_variable,
                          const Variable<double>& destination_variable,
                          Kratos::Flags& options) = 0;

  /* This function maps from Destination to Origin */
  virtual void InverseMap(const Variable< array_1d<double,3> >& origin_variable,
                          const Variable< array_1d<double,3> >& destination_variable,
                          Kratos::Flags& options) = 0;

  MapperCommunicator::Pointer GetMapperCommunicator() {
      return m_p_mapper_communicator;
  }

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

  Parameters m_json_parameters;

  MapperCommunicator::Pointer m_p_mapper_communicator;

  // global, aka of the entire submodel-parts
  int m_num_conditions_origin;
  int m_num_conditions_destination;

  int m_num_nodes_origin;
  int m_num_nodes_destination;

  int m_echo_level = 0;

  ///@}
  ///@name Protected Operators
  ///@{

  ///@}
  ///@name Protected Operations
  ///@{

  // Constructor, can only be called by derived classes (actual mappers)
  Mapper(ModelPart& rModelPartOrigin, ModelPart& rModelPartDestination,
         Parameters JsonParameters) :
  	  m_model_part_origin(rModelPartOrigin),
      m_model_part_destination(rModelPartDestination),
      m_json_parameters(JsonParameters) {

      ComputeNumberOfNodesAndConditions();

      m_echo_level = JsonParameters["echo_level"].GetInt();

      // Create the mapper communicator
      #ifdef KRATOS_USING_MPI // mpi-parallel compilation
          int mpi_initialized;
          MPI_Initialized(&mpi_initialized);
          if (mpi_initialized) { // parallel execution, i.e. mpi imported in python
              int comm_size;
              MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
              if (comm_size > 1) {
                  m_p_mapper_communicator = MapperCommunicator::Pointer (
                      new MapperMPICommunicator(m_model_part_origin,
                                                m_model_part_destination,
                                                m_json_parameters) );
              } else { // mpi importet in python, but execution with one process only
                  InitializeSerialCommunicator();
              }
          } else { // serial execution, i.e. mpi NOT imported in python
              InitializeSerialCommunicator();
          }

      #else // serial compilation
          InitializeSerialCommunicator();
      #endif
  }

  void ComputeNumberOfNodesAndConditions() {
    // Compute the quantities of the local model_parts
    m_num_conditions_origin = m_model_part_origin.GetCommunicator().LocalMesh().NumberOfConditions();
    m_num_conditions_destination = m_model_part_destination.GetCommunicator().LocalMesh().NumberOfConditions();

    m_num_nodes_origin = m_model_part_origin.GetCommunicator().LocalMesh().NumberOfNodes();
    m_num_nodes_destination = m_model_part_destination.GetCommunicator().LocalMesh().NumberOfNodes();

    // Compute the quantities of the global model_parts
    m_model_part_origin.GetCommunicator().SumAll(m_num_conditions_origin);
    m_model_part_destination.GetCommunicator().SumAll(m_num_conditions_destination);

    m_model_part_origin.GetCommunicator().SumAll(m_num_nodes_origin);
    m_model_part_destination.GetCommunicator().SumAll(m_num_nodes_destination);
  }

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

  void InitializeSerialCommunicator() {
      m_p_mapper_communicator = MapperCommunicator::Pointer (
          new MapperCommunicator(m_model_part_origin,
                                 m_model_part_destination,
                                 m_json_parameters) );
  }

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