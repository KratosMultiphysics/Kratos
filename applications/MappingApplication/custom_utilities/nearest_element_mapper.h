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

#if !defined(KRATOS_NEAREST_ELEMENT_MAPPER_H_INCLUDED )
#define  KRATOS_NEAREST_ELEMENT_MAPPER_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"

#include "mapper.h"


namespace Kratos
{
  ///@addtogroup ApplicationNameApplication
  ///@{

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
  class NearestElementMapper : public Mapper
    {
    public:
      ///@name Type Definitions
      ///@{

      /// Pointer definition of NearestElementMapper
      KRATOS_CLASS_POINTER_DEFINITION(NearestElementMapper);

      ///@}
      ///@name Life Cycle
      ///@{

      NearestElementMapper(ModelPart& i_model_part_origin, ModelPart& i_model_part_destination,
                           Parameters& rJsonParameters) : Mapper(
                           i_model_part_origin, i_model_part_destination, rJsonParameters) {
          m_p_mapper_communicator->InitializeOrigin(MapperUtilities::Condition_Center);
          m_p_mapper_communicator->InitializeDestination(MapperUtilities::Node);
          m_p_mapper_communicator->Initialize();

          m_p_inverse_mapper.reset(); // explicitly specified to be safe
      }

      /// Destructor.
      virtual ~NearestElementMapper(){ }


      ///@}
      ///@name Operators
      ///@{


      ///@}
      ///@name Operations
      ///@{

      void UpdateInterface(Kratos::Flags& options, double search_radius) override {
          m_p_mapper_communicator->UpdateInterface(options, search_radius);
          if (m_p_inverse_mapper) {
              m_p_inverse_mapper->UpdateInterface(options, search_radius);
          }

          if (options.Is(MapperFlags::REMESHED)) {
              ComputeNumberOfNodesAndConditions();
          }
      }

      /* This function maps from Origin to Destination */
      void Map(const Variable<double>& origin_variable,
               const Variable<double>& destination_variable,
               Kratos::Flags& options) override {
          m_p_mapper_communicator->TransferInterpolatedData(origin_variable,
                                                            destination_variable,
                                                            options);
      }

      /* This function maps from Origin to Destination */
      void Map(const Variable< array_1d<double,3> >& origin_variable,
               const Variable< array_1d<double,3> >& destination_variable,
               Kratos::Flags& options) override {
          m_p_mapper_communicator->TransferInterpolatedData(origin_variable,
                                                            destination_variable,
                                                            options);
        }

      /* This function maps from Destination to Origin */
      void InverseMap(const Variable<double>& origin_variable,
                      const Variable<double>& destination_variable,
                      Kratos::Flags& options) override {
          // Construct the inverse mapper if it hasn't been done before
          // It is constructed with the order of the model_parts changed!
          if (!m_p_inverse_mapper) {
              m_p_inverse_mapper = Mapper::Pointer( new NearestElementMapper(m_model_part_destination,
                                                                             m_model_part_origin,
                                                                             m_json_parameters) );
          }
          m_p_inverse_mapper->Map(destination_variable, origin_variable, options);
      }

      /* This function maps from Destination to Origin */
      void InverseMap(const Variable< array_1d<double,3> >& origin_variable,
                      const Variable< array_1d<double,3> >& destination_variable,
                      Kratos::Flags& options) override {
          // Construct the inverse mapper if it hasn't been done before
          // It is constructed with the order of the model_parts changed!
          if (!m_p_inverse_mapper) {
              m_p_inverse_mapper = Mapper::Pointer( new NearestElementMapper(m_model_part_destination,
                                                                             m_model_part_origin,
                                                                             m_json_parameters) );
          }
          m_p_inverse_mapper->Map(destination_variable, origin_variable, options);
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
      virtual std::string Info() const
      {
	      std::stringstream buffer;
        buffer << "NearestElementMapper" ;
        return buffer.str();
      }

      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const {rOStream << "NearestElementMapper";}

      /// Print object's data.
      virtual void PrintData(std::ostream& rOStream) const {}


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

      Mapper::Pointer m_p_inverse_mapper;

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
      NearestElementMapper& operator=(NearestElementMapper const& rOther);

    //   /// Copy constructor.
    //   NearestElementMapper(NearestElementMapper const& rOther){}


      ///@}

    }; // Class NearestElementMapper

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream,
				    NearestElementMapper& rThis)
    {
        return rIStream;
    }

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream,
				    const NearestElementMapper& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }
  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_NEAREST_ELEMENT_MAPPER_H_INCLUDED  defined
