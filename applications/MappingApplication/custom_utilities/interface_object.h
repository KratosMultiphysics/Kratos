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
#include "../mapping_application_variables.h"


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
      // TODO hide constructors that are only called from the derived classes
      InterfaceObject() : Point<3>(0.0f, 0.0f, 0.0f) {  // Default Constructor
          SetInitialValuesToMembers();
      }

      InterfaceObject(Node<3>& rNode) : Point<3>(rNode){ // constuct from node
          SetInitialValuesToMembers();
      }

      InterfaceObject(array_1d<double, 3> Coords) : Point<3>(Coords) { // constuct from coordinate-array
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

      bool IsInBoundingBox(double* pBoundingBox[]){
          // xmax, xmin,  ymax, ymin,  zmax, zmin
          bool is_inside = false;

          if (this->X() < *pBoundingBox[0] && this->X() > *pBoundingBox[1]) { // check x-direction
              if (this->Y() < *pBoundingBox[2] && this->Y() > *pBoundingBox[3]) { // check y-direction
                  if (this->Z() < *pBoundingBox[4] && this->Z() > *pBoundingBox[5]) { // check z-direction
                      is_inside = true;
                  }
              }
          }
          return is_inside;
      }

      void ProcessSearchResult(const double Distance, const int PairingStatus, const int Rank) {
          if (mPairingStatus < PairingStatus || (mPairingStatus == PairingStatus 
                                                 && Distance < mMinDistanceNeighbor)) {
              mPairingStatus = PairingStatus;
              mMinDistanceNeighbor = Distance;
              mNeighborRank = Rank;
          }   
      }

      int GetPairingStatus() {
          return mPairingStatus;
      }

      bool NeighborFound() {
          if (mPairingStatus == this->GetIndexNeighborFound())  {
              return true;
          } else {
              return false;
          }
      }

      bool NeighborOrApproximationFound() {
          if (mPairingStatus >= this->GetIndexApproximation())  {
              return true;
          } else {
              return false;
          }
      }

      bool HasNeighborInPartition(const int PartitionIndex) {
          bool return_value = false;
          if (mPairingStatus == this->GetIndexNeighborFound()) {
              if (mNeighborRank == PartitionIndex)
                  return_value = true;
          }
          return return_value;
      }

      int GetIndexNoNeighbor() {
          return 0;
      }

      int GetIndexApproximation() {
          return 1;
      }

      int GetIndexNeighborFound() {
          return 2;
      }

      // void SetIsBeingSent() {
      //     m_is_being_sent = true;
      // }

      // bool GetIsBeingSent() {
      //     return m_is_being_sent;
      // }

      int GetNeighborRank() {
          return mNeighborRank;
      }

      virtual bool EvaluateResult(const array_1d<double, 3> global_coords, double& min_distance,
                                  double distance, array_1d<double,2>& local_coords,
                                  std::vector<double>& shape_function_values) {
          KRATOS_ERROR << "Base class function called!" << std::endl;
          return false;
      }

      virtual bool ComputeApproximation(array_1d<double, 3> global_coords, double& min_distance,
                                        double distance, std::vector<double>& shape_function_values) {
          KRATOS_ERROR << "Base class function called!" << std::endl;
          return false;
      }

      // These functions have to be duplicated because virtual templates are not possible in C++
      // Scalars
      virtual double GetObjectValue(const Variable<double>& variable,
                                    const Kratos::Flags& options) {
          KRATOS_ERROR << "Base class function called!" << std::endl;
      }

      virtual void SetObjectValue(const Variable<double>& variable,
                                  const double value,
                                  const Kratos::Flags& options,
                                  const double factor) {
          KRATOS_ERROR << "Base class function called!" << std::endl;
      }

      virtual double GetObjectValueInterpolated(const Variable<double>& rVariable,
                                                std::vector<double>& rShapeFunctionValues) {
          KRATOS_ERROR << "Base class function called!" << std::endl;
      }

      // Vectors
      virtual array_1d<double,3> GetObjectValue(const Variable< array_1d<double,3> >& variable,
                                                const Kratos::Flags& options) {
          KRATOS_ERROR << "Base class function called!" << std::endl;
      }

      virtual void SetObjectValue(const Variable< array_1d<double,3> >& variable,
                                  const array_1d<double,3>& value,
                                  const Kratos::Flags& options,
                                  const double factor) {
          KRATOS_ERROR << "Base class function called!" << std::endl;
      }

      virtual array_1d<double,3> GetObjectValueInterpolated(const Variable< array_1d<double,3> >& variable,
                                                            std::vector<double>& shape_function_values) {
          KRATOS_ERROR << "Base class function called!" << std::endl;
      }

      // Functions used for Debugging
      virtual int GetObjectId() {
          KRATOS_ERROR << "Base class function called!" << std::endl;
          return -1;
      }

      virtual void PrintNeighbors(const int CommRank) {
          KRATOS_ERROR << "Base class function called!" << std::endl;
      }
      
      virtual void WriteRankAndCoordinatesToVariable(const int CommRank) { 
          KRATOS_ERROR << "Base class function called!" << std::endl;
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

      void PrintMatchInfo(const std::string& rInterfaceObjectType, const int Id, 
                          const int CommRank, const int NeighborCommRank, 
                          array_1d<double, 3>& rNeighborCoordinates)   {

          std::cout << rInterfaceObjectType << " [" 
                    << this->X() << " "
                    << this->Y() << " " 
                    << this->Z() << "], "
                     << "Rank " << CommRank
                    << ", Id " << Id << " || Neighbor ["
                    << rNeighborCoordinates[0] << " " 
                    << rNeighborCoordinates[1] << " "
                    << rNeighborCoordinates[2] << "], Rank "
                    << NeighborCommRank << std::endl;
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

      double mMinDistanceNeighbor;
      int mPairingStatus; // 0 no Neighbor found; 1 approximation (i.e. nearest Node found); 2 match found
      int mNeighborRank;

      ///@}
      ///@name Private Operators
      ///@{


      ///@}
      ///@name Private Operations
      ///@{

      void SetInitialValuesToMembers() {
          mMinDistanceNeighbor = std::numeric_limits<double>::max();
          mPairingStatus = this->GetIndexNoNeighbor();
          mNeighborRank = 0;
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
