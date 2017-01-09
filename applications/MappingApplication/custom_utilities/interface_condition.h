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
          InterfaceObject(coords), m_p_condition(&i_condition) {

          // in case some calculations have to be done in order to construct the base element (e.g. element center),
          // then first call the empty standard constructor, then do the calculations and populate
          // the base class stuff, although populating the derived class (i.e. "this") should also be ok

          m_geometry_family = i_condition.GetGeometry().GetGeometryFamily();
          m_num_points = m_p_condition->GetGeometry().PointsNumber();
      }

      /// Destructor.
      virtual ~InterfaceCondition() { }


      ///@}
      ///@name Operators
      ///@{


      ///@}
      ///@name Operations
      ///@{

      int GetObjectId() override {
          return m_p_condition->Id();
      }

      void PrintMatchInfo() override {
          std::cout << "InteraceCondition; Id = " << GetObjectId()
                    << "; Coordinates = [" << this->X() << " "
                    << this->Y() << " " << this->Z() << "]";
      }

      bool EvaluateResult(array_1d<double, 3> global_coords, double& min_distance,
                          double distance, array_1d<double,2>& local_coords,
                          std::vector<double>& shape_function_values) override { // I am an object in the bins
          // distance is the distance to the center and not the projection distance, therefore it is unused
          bool is_closer = false;
          bool is_inside;
          double projection_distance;
          array_1d<double,2> projection_local_coords;

          if (m_geometry_family == GeometryData::Kratos_Linear) { // I am a line condition
              is_inside = MapperUtilities::ProjectPointToLine(m_p_condition, global_coords,
                                                              projection_local_coords,
                                                              projection_distance, m_normal);
          } else if (m_geometry_family == GeometryData::Kratos_Triangle) { // I am a triangular condition
              is_inside = MapperUtilities::ProjectPointToTriangle(m_p_condition, global_coords,
                                                                  projection_local_coords,
                                                                  projection_distance, m_normal);
          } else if (m_geometry_family == GeometryData::Kratos_Quadrilateral) { // I am a quadrilateral condition
              is_inside = MapperUtilities::ProjectPointToQuaddrilateral(m_p_condition, global_coords,
                                                                        projection_local_coords,
                                                                        projection_distance, m_normal);
          } else {
              KRATOS_ERROR << "MappingApplication; InterfaceCondition; \"EvaluateResult\" Unsupported Geometry!" << std::endl;
          }

          if (is_inside) {
              if (projection_distance < min_distance){
                  min_distance = projection_distance;
                  local_coords = projection_local_coords;
                  shape_function_values.resize(m_num_points);
                  for (int i = 0; i < m_num_points; ++i) {
                      shape_function_values[i] = m_p_condition->GetGeometry().ShapeFunctionValue(i, local_coords);
                  }
                  is_closer = true;
              }
          }
          // array_1d<double,2> local_coords_temppppp;
          // local_coords_temppppp[0] = 0.5;
          // local_coords_temppppp[1] = 0.0;
          // int num_objects = m_p_condition->GetGeometry().PointsNumber();
          // shape_function_values.resize(num_objects);
          // for (int i = 0; i < num_objects; ++i) {
          //     shape_function_values[i] = m_p_condition->GetGeometry().ShapeFunctionValue(i, local_coords);
          //     KRATOS_WATCH(shape_function_values[i])
          // }

          return is_closer;
      }


      double GetObjectValueInterpolated(const Variable<double>& variable,
                                        std::vector<double>& shape_function_values) override {
          double interpolated_value = 0.0f;
          for (int i = 0; i < m_num_points; ++i) {
              interpolated_value += m_p_condition->GetGeometry().GetPoint(i).FastGetSolutionStepValue(variable) * shape_function_values[i];
          }
          return interpolated_value;
      }

      array_1d<double,3> GetObjectValueInterpolated(const Variable< array_1d<double,3> >& variable,
                                                    std::vector<double>& shape_function_values) override {
          array_1d<double,3> interpolated_value;
          interpolated_value[0] = 0.0f;
          interpolated_value[1] = 0.0f;
          interpolated_value[2] = 0.0f;
          for (int i = 0; i < m_num_points; ++i) {
              interpolated_value[0] += m_p_condition->GetGeometry().GetPoint(i).FastGetSolutionStepValue(variable)[0] * shape_function_values[i];
              interpolated_value[1] += m_p_condition->GetGeometry().GetPoint(i).FastGetSolutionStepValue(variable)[1] * shape_function_values[i];
              interpolated_value[2] += m_p_condition->GetGeometry().GetPoint(i).FastGetSolutionStepValue(variable)[2] * shape_function_values[i];
          }
          return interpolated_value;
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

      Condition* m_p_condition;
      GeometryData::KratosGeometryFamily m_geometry_family;
      int m_num_points;
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
