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

#if !defined(KRATOS_INTERFACE_CONDITION_INCLUDED_H_INCLUDED )
#define  KRATOS_INTERFACE_CONDITION_INCLUDED_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "interface_object.h"



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
  class InterfaceCondition : public InterfaceObject
    {
    public:
      ///@name Type Definitions
      ///@{

      /// Pointer definition of InterfaceCondition
      KRATOS_CLASS_POINTER_DEFINITION(InterfaceCondition);

      ///@}
      ///@name Life Cycle
      ///@{

      /// Default constructor.
      InterfaceCondition(Condition& i_condition, array_1d<double, 3> coords) :
          InterfaceObject(coords), mpCondition(&i_condition) {

          // in case some calculations have to be done in order to construct the base element (e.g. element center),
          // then first call the empty standard constructor, then do the calculations and populate
          // the base class stuff, although populating the derived class (i.e. "this") should also be ok

          mGeometryFamily = i_condition.GetGeometry().GetGeometryFamily();
          mNumPoints = mpCondition->GetGeometry().PointsNumber();
      }

      /// Destructor.
      virtual ~InterfaceCondition() { }


      ///@}
      ///@name Operators
      ///@{


      ///@}
      ///@name Operations
      ///@{

      bool EvaluateResult(array_1d<double, 3> global_coords, double& min_distance,
                          double distance, array_1d<double,2>& local_coords,
                          std::vector<double>& shape_function_values) override { // I am an object in the bins
          // distance is the distance to the center and not the projection distance, therefore it is unused
          bool is_closer = false;
          bool is_inside;
          double projection_distance;
          array_1d<double,2> projection_local_coords;

          if (mGeometryFamily == GeometryData::Kratos_Linear) { // I am a line condition
              is_inside = MapperUtilities::ProjectPointToLine(mpCondition, global_coords,
                                                              projection_local_coords,
                                                              projection_distance, m_normal);
          } else if (mGeometryFamily == GeometryData::Kratos_Triangle) { // I am a triangular condition
              is_inside = MapperUtilities::ProjectPointToTriangle(mpCondition, global_coords,
                                                                  projection_local_coords,
                                                                  projection_distance, m_normal);
          } else if (mGeometryFamily == GeometryData::Kratos_Quadrilateral) { // I am a quadrilateral condition
              is_inside = MapperUtilities::ProjectPointToQuaddrilateral(mpCondition, global_coords,
                                                                        projection_local_coords,
                                                                        projection_distance, m_normal);
          } else {
              KRATOS_ERROR << "MappingApplication; InterfaceCondition; \"EvaluateResult\" Unsupported Geometry!" << std::endl;
          }

          projection_distance = fabs(projection_distance);

          if (is_inside) {
              if (projection_distance < min_distance){
                  min_distance = projection_distance;
                  local_coords = projection_local_coords;
                  shape_function_values.resize(mNumPoints);
                  for (int i = 0; i < mNumPoints; ++i) {
                      shape_function_values[i] = mpCondition->GetGeometry().ShapeFunctionValue(i, local_coords);
                  }
                  is_closer = true;
              }
          }

          return is_closer;
      }

            bool ComputeApproximation(const array_1d<double, 3> GlobalCoords, double& rMinDistance,
                                double Distance, std::vector<double>& rShapeFunctionValues) override { // I am an object in the bins
          bool is_closer = false;
          double distance_point = std::numeric_limits<double>::max();
          int closest_point_index = -1;
          // Loop over all points of the geometry
          for (int i = 0; i < mNumPoints; ++i) {
              distance_point = sqrt(pow(GlobalCoords[0] - mpCondition->GetGeometry().GetPoint(i).X() , 2) +
                                    pow(GlobalCoords[1] - mpCondition->GetGeometry().GetPoint(i).Y() , 2) +
                                    pow(GlobalCoords[2] - mpCondition->GetGeometry().GetPoint(i).Z() , 2));

              if (distance_point < rMinDistance && distance_point <= mApproximationTolerance) {
                  rMinDistance = distance_point;
                  closest_point_index = i;
                  is_closer = true;
              }
          }

          if (is_closer) {
              rShapeFunctionValues.resize(mNumPoints);
              for (int i = 0; i < mNumPoints; ++i) { 
                  if (i == closest_point_index) {
                      rShapeFunctionValues[i] = 1.0f;
                  } else {
                      rShapeFunctionValues[i] = 0.0f;
                  }
              }
          }
          
          return is_closer;
      }


      double GetObjectValueInterpolated(const Variable<double>& variable,
                                        std::vector<double>& shape_function_values) override {
          double interpolated_value = 0.0f;
          for (int i = 0; i < mNumPoints; ++i) {
              interpolated_value += mpCondition->GetGeometry().GetPoint(i).FastGetSolutionStepValue(variable) * shape_function_values[i];
          }
          return interpolated_value;
      }

      array_1d<double,3> GetObjectValueInterpolated(const Variable< array_1d<double,3> >& variable,
                                                    std::vector<double>& shape_function_values) override {
          array_1d<double,3> interpolated_value;
          interpolated_value[0] = 0.0f;
          interpolated_value[1] = 0.0f;
          interpolated_value[2] = 0.0f;
          for (int i = 0; i < mNumPoints; ++i) {
              interpolated_value[0] += mpCondition->GetGeometry().GetPoint(i).FastGetSolutionStepValue(variable)[0] * shape_function_values[i];
              interpolated_value[1] += mpCondition->GetGeometry().GetPoint(i).FastGetSolutionStepValue(variable)[1] * shape_function_values[i];
              interpolated_value[2] += mpCondition->GetGeometry().GetPoint(i).FastGetSolutionStepValue(variable)[2] * shape_function_values[i];
          }
          return interpolated_value;
      }

      // Functions used for Debugging
      int GetObjectId() override {
          return mpCondition->Id();
      }

      void PrintNeighbors(const int CommRank) override {
          array_1d<double, 3> neighbor_coordinates = mpCondition->GetValue(NEIGHBOR_COORDINATES);
          double neighbor_comm_rank = mpCondition->GetValue(NEIGHBOR_RANK);

          PrintMatchInfo("InterfaceCondition", GetObjectId(), CommRank, 
                         neighbor_comm_rank, neighbor_coordinates);
      }

      void WriteRankAndCoordinatesToVariable(const int CommRank) override {
          // This function writes the coordinates and the rank of the 
          // InterfaceObject to the variables "NEIGHBOR_COORDINATES" 
          // and "NEIGHBOR_RANK", for debugging
          array_1d<double,3> neighbor_coordinates;
          neighbor_coordinates[0] = this->X();
          neighbor_coordinates[1] = this->Y();
          neighbor_coordinates[2] = this->Z();
          mpCondition->SetValue(NEIGHBOR_COORDINATES, neighbor_coordinates);
          mpCondition->SetValue(NEIGHBOR_RANK, CommRank);
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
         buffer << "InterfaceCondition" ;
         return buffer.str();
      }

      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const {rOStream << "InterfaceCondition";}

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

      Condition* mpCondition;
      GeometryData::KratosGeometryFamily mGeometryFamily;
      int mNumPoints;
      double mApproximationTolerance = 0.0f;
      
      array_1d<double,3> m_normal;
      static constexpr double tol = 1e-4;


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
      InterfaceCondition& operator=(InterfaceCondition const& rOther);

    //   /// Copy constructor.
    //   InterfaceCondition(InterfaceCondition const& rOther){}


      ///@}

    }; // Class InterfaceCondition

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream,
				    InterfaceCondition& rThis)
  {
      return rIStream;
  }

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream,
				    const InterfaceCondition& rThis)
  {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
  }
  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_INTERFACE_CONDITION_INCLUDED_H_INCLUDED  defined
