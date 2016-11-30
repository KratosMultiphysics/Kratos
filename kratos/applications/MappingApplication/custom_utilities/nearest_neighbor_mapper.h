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
//

#if !defined(KRATOS_NEAREST_NEIGHBOR_MAPPER_H_INCLUDED )
#define  KRATOS_NEAREST_NEIGHBOR_MAPPER_H_INCLUDED


// System includes

// External includes
//@{KRATOS_EXTERNA_INCLUDES}
// #include "includes/kratos_flags.h"
// #include <boost/python.hpp>
//
// // Project includes
//
// #include "includes/define.h"
// #include "includes/kratos_flags.h"
// #include "includes/element.h"
// #include "includes/model_part.h"
// #include "geometries/geometry_data.h"

// #include "spaces/ublas_space.h"
// #include "linear_solvers/linear_solver.h"
// #include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"
// #include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"
// #include "solving_strategies/strategies/residualbased_linear_strategy.h"
// #include "elements/distance_calculation_element_simplex.h"
#include "mapper.h"

// #include <unordered_map>

// External includes

namespace Kratos {

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{
	typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
	typedef UblasSpace<double, Matrix, Vector> DenseSpaceType;
	typedef typename SparseSpaceType::MatrixType SparseMatrixType;
	typedef typename SparseSpaceType::VectorType VectorType;

	typedef array_1d<double,3> array_3d;
	typedef Node < 3 > PointType;
	typedef Node < 3 > ::Pointer PointTypePointer;
	typedef std::vector<PointType::Pointer> PointVector;
	typedef std::vector<PointType::Pointer>::iterator PointIterator;
	typedef std::vector<double> DistanceVector;
	typedef std::vector<double>::iterator DistanceIterator;
	typedef ModelPart::ConditionsContainerType ConditionsArrayType;

	typedef Element BaseType;
	typedef BaseType::GeometryType GeometryType;

    typedef std::vector<PointTypePointer> neighborsVector;
    typedef matrix<int> GraphType; // GraphColoringProcess

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

class NearestNeighborMapper : public Mapper
{
public:

  ///@name Type Definitions
  ///@{

  ///@}
  ///@name Pointer Definitions
  /// Pointer definition of NearestNeighborMapper
  KRATOS_CLASS_POINTER_DEFINITION(NearestNeighborMapper);

  ///@}
  ///@name Life Cycle
  ///@{

  NearestNeighborMapper(ModelPart& i_model_part_origin, ModelPart& i_model_part_destination) : Mapper(
                      i_model_part_origin, i_model_part_destination) {

      m_point_comm_manager_origin = Kratos::InterfaceObjectManager::CreateInterfaceNodeManager(m_model_part_origin,
          m_mapper_communicator->MyPID(), m_mapper_communicator->TotalProcesses());
      m_point_comm_manager_destination = Kratos::InterfaceObjectManager::CreateInterfaceNodeManager(m_model_part_destination,
          m_mapper_communicator->MyPID(), m_mapper_communicator->TotalProcesses());

      bool bi_directional_search = true;

      m_mapper_communicator->Initialize(m_point_comm_manager_origin, m_point_comm_manager_destination,
                                        bi_directional_search);

      MPI_Barrier(MPI_COMM_WORLD);
  }

  NearestNeighborMapper(ModelPart& i_model_part_origin, ModelPart& i_model_part_destination,
                        double i_initial_search_radius, int i_max_search_iterations) : Mapper(
                        i_model_part_origin, i_model_part_destination) {

      m_point_comm_manager_origin = Kratos::InterfaceObjectManager::CreateInterfaceNodeManager(m_model_part_origin, m_mapper_communicator->MyPID(), m_mapper_communicator->TotalProcesses());
      m_point_comm_manager_destination = Kratos::InterfaceObjectManager::CreateInterfaceNodeManager(m_model_part_destination, m_mapper_communicator->MyPID(), m_mapper_communicator->TotalProcesses());

      bool bi_directional_search = true;

      m_mapper_communicator->Initialize(m_point_comm_manager_origin, m_point_comm_manager_destination,
                                        bi_directional_search, i_initial_search_radius, i_max_search_iterations);

      MPI_Barrier(MPI_COMM_WORLD);
  }



  /// Destructor.
  virtual ~NearestNeighborMapper() {
  }

  ///@}
  ///@name Operators
  ///@{

  ///@}
  ///@name Operations
  ///@{

  void UpdateInterface() override {
      m_mapper_communicator->ComputeSearchStructure();
  }

  /* This function maps a variable from Origin to Destination */
  void Map(const Variable<double>& origin_variable,
           const Variable<double>& destination_variable,
           const bool add_value) override {
      bool direction = true;
      m_mapper_communicator->TransferData(direction, origin_variable,
                                          destination_variable, add_value);
  }

  /* This function maps a variable from Origin to Destination */
  void Map(const Variable< array_1d<double,3> >& origin_variable,
           const Variable< array_1d<double,3> >& destination_variable,
           const bool add_value) override {
      bool direction = true;
      m_mapper_communicator->TransferData(direction, origin_variable,
                                          destination_variable, add_value);
  }

  /* This function maps a variable from Destination to Origin */
  void InverseMap(const Variable<double>& origin_variable,
                  const Variable<double>& destination_variable,
                  const bool add_value) override { // TvalueType => MPI communicator
      bool direction = false;
      m_mapper_communicator->TransferData(direction, origin_variable,
                                          destination_variable, add_value);
  }

  /* This function maps a variable from Destination to Origin */
  void InverseMap(const Variable< array_1d<double,3> >& origin_variable,
                  const Variable< array_1d<double,3> >& destination_variable,
                  const bool add_value) override {
      bool direction = false;
      m_mapper_communicator->TransferData(direction, origin_variable,
                                          destination_variable, add_value);
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
      return "NearestNeighborMapper";
  }

  /// Print information about this object.
  virtual void PrintInfo(std::ostream& rOStream) const {
    rOStream << "NearestNeighborMapper";
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
  NearestNeighborMapper& operator=(NearestNeighborMapper const& rOther);

  /// Copy constructor.
  //NearestNeighborMapper(NearestNeighborMapper const& rOther);

  ///@}

}; // Class NearestNeighborMapper

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function


inline std::istream & operator >>(std::istream& rIStream,
		NearestNeighborMapper& rThis)
{
    return rIStream;
}

/// output stream function

inline std::ostream & operator <<(std::ostream& rOStream,
                                  const NearestNeighborMapper& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block
}  // namespace Kratos.

#endif // KRATOS_NEAREST_NEIGHBOR_MAPPER_H_INCLUDED  defined
