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

#if !defined(KRATOS_INTERFACE_OBJECT_INCLUDED_H_INCLUDED )
#define  KRATOS_INTERFACE_OBJECT_INCLUDED_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "mapper_utilities.h"


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
  class InterfaceObject : public Point<3>
    {
    public:
      ///@name Type Definitions
      ///@{

      /// Pointer definition of InterfaceObject
      KRATOS_CLASS_POINTER_DEFINITION(InterfaceObject);

      ///@}
      ///@name Life Cycle
      ///@{

      InterfaceObject() : Point<3>(0.0f, 0.0f, 0.0f) {  // Default Constructor
          SetInitialValuesToMembers();
      }

      InterfaceObject(Node<3>& i_node) : Point<3>(i_node) { // constuct from node
          SetInitialValuesToMembers();
      }

      InterfaceObject(array_1d<double, 3> coords) : Point<3>(coords) { // constuct from coordinate-array
          SetInitialValuesToMembers();
      }

      InterfaceObject(double X, double Y, double Z) : Point<3>(X, Y, Z) { // constuct from coordinates
          SetInitialValuesToMembers();
      }

      /// Destructor.
      virtual ~InterfaceObject() { }


      ///@}
      ///@name Operators
      ///@{


      ///@}
      ///@name Operations
      ///@{

      void Reset() {
          SetInitialValuesToMembers();
      }

      bool IsInBoundingBox(double* bounding_box[]){
          // xmax, xmin,  ymax, ymin,  zmax, zmin
          bool is_inside = false;

          if (this->X() < *bounding_box[0] && this->X() > *bounding_box[1]) { // check x-direction
              if (this->Y() < *bounding_box[2] && this->Y() > *bounding_box[3]) { // check y-direction
                  if (this->Z() < *bounding_box[4] && this->Z() > *bounding_box[5]) { // check z-direction
                      is_inside = true;
                  }
              }
          }
          return is_inside;
      }

      void ProcessDistance(double distance, int rank) {
          m_neighbor_found = true;
          if (distance < m_min_distance_neighbor) {
              m_min_distance_neighbor = distance;
              m_neighbor_rank = rank;
          }
      }

      bool NeighborFound() {
          return m_neighbor_found;
      }

      // double MinDistance() {
      //     return m_min_distance_neighbor;
      // }

      bool HasNeighborInPartition(const int partition_index) {
          bool return_value = false;
          if (m_neighbor_found) {
              if (m_neighbor_rank == partition_index)
                  return_value = true;
          }
          return return_value;
      }

      void SetIsBeingSent() {
          m_is_being_sent = true;
      }

      bool GetIsBeingSent() {
          return m_is_being_sent;
      }

      int GetNeighborRank() {
          return m_neighbor_rank;
      }

      virtual int GetObjectId() {
          KRATOS_ERROR << "MappingApplication; InterfaceObject; \"GetObjectId\" "
                       << "of the base class called!" << std::endl;
          return -1;
      }

      virtual void PrintMatchInfo() {
          KRATOS_ERROR << "MappingApplication; InterfaceObject; \"PrintMatchInfo\" "
                       << "of the base class called!" << std::endl;
      }

      virtual bool EvaluateResult(array_1d<double, 3> global_coords, double& min_distance,
                                  double distance, array_1d<double,2>& local_coords,
                                  std::vector<double>& shape_function_values) {
          KRATOS_ERROR << "MappingApplication; InterfaceObject; \"EvaluateResult\" "
                       << "of the base class called!" << std::endl;
          return false;
      }

      // These functions have to be duplicated because virtual templates are not possible in C++
      // Scalars
      virtual double GetObjectValue(const Variable<double>& variable) {
          KRATOS_ERROR << "MappingApplication; InterfaceObject; \"GetObjectValue, "
                       << "double\" of the base class called!" << std::endl;
      }

      virtual void SetObjectValue(const Variable<double>& variable,
                                  const double value,
                                  const Kratos::Flags& options,
                                  const double factor) {
          KRATOS_ERROR << "MappingApplication; InterfaceObject; \"SetObjectValue, "
                       << "double\" of the base class called!" << std::endl;
      }

      virtual double GetObjectValueInterpolated(const Variable<double>& variable,
                                                std::vector<double>& shape_function_values) {
          KRATOS_ERROR << "MappingApplication; InterfaceObject; "
                       << "\"GetObjectValueInterpolated, double\" of the base "
                       << "class called!" << std::endl;
      }

      // Vectors
      virtual array_1d<double,3> GetObjectValue(const Variable< array_1d<double,3> >& variable) {
          KRATOS_ERROR << "MappingApplication; InterfaceObject; \"GetObjectValue, "
                       << "double<3>\" of the base class called!" << std::endl;
      }

      virtual void SetObjectValue(const Variable< array_1d<double,3> >& variable,
                                  const array_1d<double,3>& value,
                                  const Kratos::Flags& options,
                                  const double factor) {
          KRATOS_ERROR << "MappingApplication; InterfaceObject; \"SetObjectValue, "
                       << "double<3>\" of the base class called!" << std::endl;
      }

      virtual array_1d<double,3> GetObjectValueInterpolated(const Variable< array_1d<double,3> >& variable,
                                                            std::vector<double>& shape_function_values) {
          KRATOS_ERROR << "MappingApplication; InterfaceObject; "
                       << "\"GetObjectValueInterpolated, double<3>\" of the base "
                       << "class called!" << std::endl;
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
         buffer << "InterfaceObject" ;
         return buffer.str();
      }

      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const {rOStream << "InterfaceObject";}

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

      double m_min_distance_neighbor;
      bool m_neighbor_found;
      bool m_is_being_sent;
      int m_neighbor_rank;

      ///@}
      ///@name Private Operators
      ///@{


      ///@}
      ///@name Private Operations
      ///@{

      void SetInitialValuesToMembers() {
          m_min_distance_neighbor = std::numeric_limits<double>::max();
          m_neighbor_found = false;
          m_is_being_sent = false;
          m_neighbor_rank = 0;
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
      InterfaceObject& operator=(InterfaceObject const& rOther);

    //   /// Copy constructor.
    //   InterfaceObject(InterfaceObject const& rOther){}


      ///@}

    }; // Class InterfaceObject

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream,
				    InterfaceObject& rThis)
  {
      return rIStream;
  }

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream,
				    const InterfaceObject& rThis)
  {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
  }
  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_INTERFACE_OBJECT_INCLUDED_H_INCLUDED  defined
