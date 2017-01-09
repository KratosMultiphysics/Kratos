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

// Project includes

#include "mapper.h"

// #include <unordered_map>

// External includes

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

  NearestNeighborMapper(ModelPart& i_model_part_origin, ModelPart& i_model_part_destination,
                        Parameters& i_json_parameters) : Mapper(
                        i_model_part_origin, i_model_part_destination, i_json_parameters) {
      m_mapper_communicator->InitializeOrigin(MapperUtilities::Node);
      m_mapper_communicator->InitializeDestination(MapperUtilities::Node);
      m_mapper_communicator->Initialize();

      m_inverse_mapper.reset(); // explicitly specified to be safe

      if (i_json_parameters["non_conforming_interface"].GetBool()) {
          KRATOS_ERROR << "MappingApplication; NearestNeighborMapper; invalid "
                       << "option specified for this mapper: "
                       << "\"non_conforming_interface\"" << std::endl;
      }
  }

  /// Destructor.
  virtual ~NearestNeighborMapper() { }

  ///@}
  ///@name Operators
  ///@{

  ///@}
  ///@name Operations
  ///@{

  void UpdateInterface(Kratos::Flags& options, double search_radius) override {
      m_mapper_communicator->UpdateInterface(options, search_radius);
      if (m_inverse_mapper) {
          m_inverse_mapper->UpdateInterface(options, search_radius);
      }

      if (options.Is(MapperFlags::REMESHED)) {
          ComputeNumberOfNodesAndConditions();
      }
  }

  /* This function maps from Origin to Destination */
  void Map(const Variable<double>& origin_variable,
           const Variable<double>& destination_variable,
           Kratos::Flags& options) override {
      double factor = 1.0f;

      if (options.Is(MapperFlags::CONSERVATIVE)) {
          factor = MapperUtilities::ComputeConservativeFactor(
              m_num_nodes_origin,
              m_num_nodes_destination);
      }

      m_mapper_communicator->TransferNodalData(origin_variable,
                                               destination_variable,
                                               options,
                                               factor);
  }

  /* This function maps from Origin to Destination */
  void Map(const Variable< array_1d<double,3> >& origin_variable,
           const Variable< array_1d<double,3> >& destination_variable,
           Kratos::Flags& options) override {
      double factor = 1.0f;

      if (options.Is(MapperFlags::CONSERVATIVE)) {
          factor = MapperUtilities::ComputeConservativeFactor(
              m_num_nodes_origin,
              m_num_nodes_destination);
      }

      m_mapper_communicator->TransferNodalData(origin_variable,
                                               destination_variable,
                                               options,
                                               factor);
    }

  /* This function maps from Destination to Origin */
  void InverseMap(const Variable<double>& origin_variable,
                  const Variable<double>& destination_variable,
                  Kratos::Flags& options) override {
      // Construct the inverse mapper if it hasn't been done before
      // It is constructed with the order of the model_parts changed!
      if (!m_inverse_mapper) {
          m_inverse_mapper = Mapper::Pointer( new NearestNeighborMapper(m_model_part_destination,
                                                                        m_model_part_origin,
                                                                        m_json_parameters) );
      }
      m_inverse_mapper->Map(destination_variable, origin_variable, options);
  }

  /* This function maps from Destination to Origin */
  void InverseMap(const Variable< array_1d<double,3> >& origin_variable,
                  const Variable< array_1d<double,3> >& destination_variable,
                  Kratos::Flags& options) override {
      // Construct the inverse mapper if it hasn't been done before
      // It is constructed with the order of the model_parts changed!
      if (!m_inverse_mapper) {
          m_inverse_mapper = Mapper::Pointer( new NearestNeighborMapper(m_model_part_destination,
                                                                        m_model_part_origin,
                                                                        m_json_parameters) );
      }
      m_inverse_mapper->Map(destination_variable, origin_variable, options);
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

  Mapper::Pointer m_inverse_mapper;

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
