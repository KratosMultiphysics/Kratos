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

// External includes

// Project includes
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
      InterfaceCondition(Condition& rCondition, array_1d<double, 3> Coords, const double ApproximationTolerance) :
          InterfaceObject(Coords), mpCondition(&rCondition), mApproximationTolerance(ApproximationTolerance) {

          // in case some calculations have to be done in order to construct the base element (e.g. element center),
          // then first call the empty standard constructor, then do the calculations and populate
          // the base class stuff, although populating the derived class (i.e. "this") should also be ok

          mGeometryFamily = mpCondition->GetGeometry().GetGeometryFamily();
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

      bool EvaluateResult(const array_1d<double, 3>& GlobalCoords, double& rMinDistance,
                          double Distance, array_1d<double,2>& local_coords,
                          std::vector<double>& rShapeFunctionValues) override { // I am an object in the bins
          // Distance is the distance to the center and not the projection distance, therefore it is unused
          bool is_closer = false;
          bool is_inside = false;
          double projection_distance;
          array_1d<double,3> projection_local_coords;

          if (mGeometryFamily == GeometryData::Kratos_Linear 
              && mNumPoints == 2) { // I am a linear line condition
              is_inside = MapperUtilities::ProjectPointToLine(mpCondition, GlobalCoords,
                                                              projection_local_coords,
                                                              projection_distance);
          } else if (mGeometryFamily == GeometryData::Kratos_Triangle 
                     && mNumPoints == 3) { // I am a linear triangular condition
              is_inside = MapperUtilities::ProjectPointToTriangle(mpCondition, GlobalCoords,
                                                                  projection_local_coords,
                                                                  projection_distance);
          } else if (mGeometryFamily == GeometryData::Kratos_Quadrilateral
                     && mNumPoints == 4) { // I am a linear quadrilateral condition
              is_inside = MapperUtilities::ProjectPointToQuadrilateral(mpCondition, GlobalCoords,
                                                                       projection_local_coords,
                                                                       projection_distance);
          } else {
              std::cout << "MAPPER WARNING, Used Geometry is not implemented, " 
                        << "using an approximation" << std::endl;
              return false;
          }

          projection_distance = fabs(projection_distance);

          if (is_inside) {
              if (projection_distance < rMinDistance){
                  rMinDistance = projection_distance;
                  rShapeFunctionValues.resize(mNumPoints);
                  for (int i = 0; i < mNumPoints; ++i) {
                      rShapeFunctionValues[i] = mpCondition->GetGeometry().ShapeFunctionValue(i, projection_local_coords);
                  }
                  is_closer = true;
              }
          }
          return is_closer;
      }

      bool ComputeApproximation(const array_1d<double, 3>& GlobalCoords, double& rMinDistance,
                                double Distance, std::vector<double>& rShapeFunctionValues) override { // I am an object in the bins
          bool is_closer = false;
          double distance_point = std::numeric_limits<double>::max();
          int closest_point_index = -1;
          // Loop over all points of the geometry
          for (int i = 0; i < mNumPoints; ++i) {
              distance_point = MapperUtilities::ComputeDistance(GlobalCoords, 
                                                                mpCondition->GetGeometry().GetPoint(i).Coordinates());

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

      
      // Scalars
      double GetObjectValue(const Variable<double>& rVariable,
                            const Kratos::Flags& options) override {
          KRATOS_ERROR_IF_NOT(options.Is(MapperFlags::NON_HISTORICAL_DATA))
              << "Only Non-Historical Variables are accessible for Conditions" << std::endl;

          return mpCondition->GetValue(rVariable);
      }

      void SetObjectValue(const Variable<double>& rVariable,
                          const double value,
                          const Kratos::Flags& options,
                          const double factor) override {
          KRATOS_ERROR_IF_NOT(options.Is(MapperFlags::NON_HISTORICAL_DATA))
              << "Only Non-Historical Variables are accessible for Conditions" << std::endl;

          if (options.Is(MapperFlags::ADD_VALUES)) {
              double old_value = mpCondition->GetValue(rVariable);
              mpCondition->SetValue(rVariable, old_value + value * factor);
          } else {
              mpCondition->SetValue(rVariable, value * factor);
          } 
      }

      double GetObjectValueInterpolated(const Variable<double>& rVariable,
                                        std::vector<double>& rShapeFunctionValues) override {
          double interpolated_value = 0.0f;
          double shape_fct_value = 0.0f;
          
          // std::cout << "\n\n" << std::endl;

          for (int i = 0; i < mNumPoints; ++i) {
              shape_fct_value += rShapeFunctionValues[i];
              // std::cout << "rShapeFunctionValues[i] " << rShapeFunctionValues[i] << std::endl;
              // std::cout << "val " << mpCondition->GetGeometry().GetPoint(i).FastGetSolutionStepValue(rVariable) << std::endl;



              interpolated_value += mpCondition->GetGeometry().GetPoint(i).FastGetSolutionStepValue(rVariable) * rShapeFunctionValues[i];
          }
          // std::cout << "[ " << this->X() << " " << this->Y() << " " << this->Z() << " ]" << std::endl;
          // std::cout << "Interpolated Value = " << interpolated_value << " ; shape_fct_value = " << shape_fct_value << "\n" << std::endl;
          return interpolated_value;
      }

      // Vectors
      array_1d<double,3> GetObjectValue(const Variable< array_1d<double,3> >& rVariable,
                                        const Kratos::Flags& options) override {
          KRATOS_ERROR_IF_NOT(options.Is(MapperFlags::NON_HISTORICAL_DATA))
              << "Only Non-Historical Variables are accessible for Conditions" << std::endl;
          
          return mpCondition->GetValue(rVariable);
      }

      void SetObjectValue(const Variable< array_1d<double,3> >& rVariable,
                          const array_1d<double,3>& value,
                          const Kratos::Flags& options,
                          const double factor) override {
          KRATOS_ERROR_IF_NOT(options.Is(MapperFlags::NON_HISTORICAL_DATA))
              << "Only Non-Historical Variables are accessible for Conditions" << std::endl;

          if (options.Is(MapperFlags::ADD_VALUES)) {
              array_1d<double,3> old_value = mpCondition->GetValue(rVariable);
              mpCondition->SetValue(rVariable, old_value + value * factor);
          } else {
              mpCondition->SetValue(rVariable, value * factor);
          }
      }

      array_1d<double,3> GetObjectValueInterpolated(const Variable< array_1d<double,3> >& rVariable,
                                                    std::vector<double>& rShapeFunctionValues) override {
          array_1d<double,3> interpolated_value;
          interpolated_value[0] = 0.0f;
          interpolated_value[1] = 0.0f;
          interpolated_value[2] = 0.0f;
          for (int i = 0; i < mNumPoints; ++i) {
              interpolated_value[0] += mpCondition->GetGeometry().GetPoint(i).FastGetSolutionStepValue(rVariable)[0] * rShapeFunctionValues[i];
              interpolated_value[1] += mpCondition->GetGeometry().GetPoint(i).FastGetSolutionStepValue(rVariable)[1] * rShapeFunctionValues[i];
              interpolated_value[2] += mpCondition->GetGeometry().GetPoint(i).FastGetSolutionStepValue(rVariable)[2] * rShapeFunctionValues[i];
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
          // TODO exchange with "Coordinates()"
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
