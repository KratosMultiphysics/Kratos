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

#if !defined(KRATOS_ITERATIVE_MORTAR_MAPPER_H_INCLUDED )
#define  KRATOS_ITERATIVE_MORTAR_MAPPER_H_INCLUDED



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
  class IterativeMortarMapper : public Mapper
    {
    public:
      ///@name Type Definitions
      ///@{

      /// Pointer definition of IterativeMortarMapper
      KRATOS_CLASS_POINTER_DEFINITION(IterativeMortarMapper);

      ///@}
      ///@name Life Cycle
      ///@{

      IterativeMortarMapper(ModelPart& i_model_part_origin, ModelPart& i_model_part_destination,
                            double i_initial_search_radius, int i_max_search_iterations,
                            double i_convergence_tolerance, int i_convergence_iterations) : Mapper(
                            i_model_part_origin, i_model_part_destination) {

          m_point_comm_manager_origin = Kratos::InterfaceObjectManager::CreateInterfaceConditionManager(m_model_part_origin,
              m_mapper_communicator->MyPID(), m_mapper_communicator->TotalProcesses(), GeometryData::GI_GAUSS_2);
          m_point_comm_manager_destination = Kratos::InterfaceObjectManager::CreateInterfaceConditionManager(m_model_part_destination,
              m_mapper_communicator->MyPID(), m_mapper_communicator->TotalProcesses(), GeometryData::GI_GAUSS_2);

          m_mapper_communicator->Initialize(m_point_comm_manager_origin, m_point_comm_manager_destination,
                                            i_initial_search_radius, i_max_search_iterations);

          MPI_Barrier(MPI_COMM_WORLD);
          std::cout << "IterativeMortarMapper Initialized" << std::endl;
      }

      /// Destructor.
      virtual ~IterativeMortarMapper(){}


      ///@}
      ///@name Operators
      ///@{


      ///@}
      ///@name Operations
      ///@{

      void UpdateInterface() override {
          m_mapper_communicator->ComputeSearchStructure();
      }

      /* This function maps from Origin to Destination */
      void Map(const Variable<double>& origin_variable,
               const Variable<double>& destination_variable,
               const bool add_value,
               const bool sign_positive) override {};

      /* This function maps from Origin to Destination */
      void Map(const Variable< array_1d<double,3> >& origin_variable,
               const Variable< array_1d<double,3> >& destination_variable,
               const bool add_value,
               const bool sign_positive) override {};

      // /* This function maps from Destination to Origin */
      // virtual void InverseMap(const Variable<double>& origin_variable,
      //                         const Variable<double>& destination_variable,
      //                         const bool add_value){};
      //
      // /* This function maps from Destination to Origin */
      // virtual void InverseMap(const Variable< array_1d<double,3> >& origin_variable,
      //                         const Variable< array_1d<double,3> >& destination_variable,
      //                         const bool add_value){};


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
        buffer << "IterativeMortarMapper" ;
        return buffer.str();
      }

      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const {rOStream << "IterativeMortarMapper";}

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
      IterativeMortarMapper& operator=(IterativeMortarMapper const& rOther);

    //   /// Copy constructor.
    //   IterativeMortarMapper(IterativeMortarMapper const& rOther){}


      ///@}

    }; // Class IterativeMortarMapper

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream,
				    IterativeMortarMapper& rThis)
    {
        return rIStream;
    }

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream,
				    const IterativeMortarMapper& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }
  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_ITERATIVE_MORTAR_MAPPER_H_INCLUDED  defined
