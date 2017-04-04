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

#if !defined(KRATOS_MAPPER_UTILITIES_H_INCLUDED )
#define  KRATOS_MAPPER_UTILITIES_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "utilities/math_utils.h" // Cross Product
#include "utilities/openmp_utils.h" // for GetCurrentTime()
#ifdef KRATOS_USING_MPI
#include "mpi.h" // for GetCurrentTime()
#endif


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
  class MapperUtilities
    {
    public:
      ///@name Type Definitions
      ///@{

      /// Pointer definition of MapperUtilities
      KRATOS_CLASS_POINTER_DEFINITION(MapperUtilities);

      ///@}
      ///@name  Enum's
      ///@{

      enum InterfaceObjectConstructionType
      {
          Node,
          Condition_Center,
          Condition_Gauss_Point,
          // Point or Coordinates or sth => for Contact => create with list of points/coords,
      };

      ///@}
      ///@name Life Cycle
      ///@{

      /// Destructor.
      virtual ~MapperUtilities() { }


      ///@}
      ///@name Operators
      ///@{

      static constexpr bool MAPPER_DEBUG_LEVEL = true;

      ///@}
      ///@name Operations
      ///@{

      static double GetCurrentTime() {
          double current_time;
          #ifdef KRATOS_USING_MPI // mpi-parallel compilation
              int mpi_initialized;
              MPI_Initialized(&mpi_initialized);
              if (mpi_initialized) { // parallel execution, i.e. mpi imported in python
                  current_time = MPI_Wtime();
              } else { // serial execution, i.e. mpi NOT imported in python
                  current_time = OpenMPUtils::GetCurrentTime();
              }
          #else // serial compilation
              current_time = OpenMPUtils::GetCurrentTime();
          #endif


          return current_time;
      }

      static int ComputeNumberOfNodes(ModelPart& model_part) {
          int num_nodes = model_part.GetCommunicator().LocalMesh().NumberOfNodes();
          model_part.GetCommunicator().SumAll(num_nodes); // Compute the sum among the partitions
          return num_nodes;
      }

      static int ComputeNumberOfConditions(ModelPart& model_part) {
          int num_conditions = model_part.GetCommunicator().LocalMesh().NumberOfConditions();
          model_part.GetCommunicator().SumAll(num_conditions); // Compute the sum among the partitions
          return num_conditions;
      }

      static double ComputeSearchRadius(ModelPart& model_part_1, ModelPart& model_part_2) {
          double search_radius = std::max(ComputeSearchRadius(model_part_1),
                                          ComputeSearchRadius(model_part_2));
          return search_radius;
      }

      static double ComputeSearchRadius(ModelPart& model_part) {
          double search_safety_factor = 1.2;
          double max_element_size = 0.0;

          int num_conditions_global = ComputeNumberOfConditions(model_part);

          if (num_conditions_global > 0) {
              // Loop through each edge of a geometrical entity ONCE
              for (auto& condition : model_part.GetCommunicator().LocalMesh().Conditions()) {
                  for (std::size_t i = 0; i < (condition.GetGeometry().size() - 1); ++i) {
                      double node_1_x = condition.GetGeometry()[i].X();
                      double node_1_y = condition.GetGeometry()[i].Y();
                      double node_1_z = condition.GetGeometry()[i].Z();

                      for (std::size_t j = i + 1; j < condition.GetGeometry().size(); ++j) {
                          double node_2_x = condition.GetGeometry()[j].X();
                          double node_2_y = condition.GetGeometry()[j].Y();
                          double node_2_z = condition.GetGeometry()[j].Z();

                          double edge_length = sqrt(pow(node_1_x - node_2_x , 2) +
                                                    pow(node_1_y - node_2_y , 2) +
                                                    pow(node_1_z - node_2_z , 2));

                          max_element_size = std::max(max_element_size, edge_length);
                      }
                  }
              }
          } else {
              // TODO modify and use again
              // std::cout << "MAPPER WARNING, no conditions for search radius "
              //           << "computations found, using nodes (less efficient)"
              //           << std::endl;
              // TODO modify loop such that it loop only once over the nodes
              for (auto& node_1 : model_part.GetCommunicator().LocalMesh().Nodes()) {
                  double node_1_x = node_1.X();
                  double node_1_y = node_1.Y();
                  double node_1_z = node_1.Z();
                  for (auto& node_2 : model_part.GetCommunicator().LocalMesh().Nodes()) {
                      double node_2_x = node_2.X();
                      double node_2_y = node_2.Y();
                      double node_2_z = node_2.Z();

                      double edge_length = sqrt(pow(node_1_x - node_2_x , 2) +
                                                pow(node_1_y - node_2_y , 2) +
                                                pow(node_1_z - node_2_z , 2));

                      max_element_size = std::max(max_element_size, edge_length);
                  }
              }
          }

          model_part.GetCommunicator().MaxAll(max_element_size); // Compute the maximum among the partitions
          return max_element_size * search_safety_factor;
      }

      static double ComputeConservativeFactor(const double NumNodesOrigin,
                                              const double NumNodesDestination) {
          // num_nodes_* are casted to doubles in order to use the double devision
          // if this function would take ints, then the return value would also be an int!
          return NumNodesOrigin / NumNodesDestination;
      }

      static bool ProjectPointToLine(Condition* pCondition,
                                     const array_1d<double, 3> GlobalCoords,
                                     array_1d<double,2>& rLocalCoords,
                                     double& rDistance) {
          // Point<3> point_to_project(GlobalCoords);

          // Point<3> point_projected;
          // Point<3> point_in_plane = pCondition->GetGeometry()[0];

          // array_1d<double,3> vector_points;
          // noalias(vector_points) = point_to_project.Coordinates() - point_in_plane.Coordinates();

          // CalculateLineNormal(pCondition, rNormal);

          // rDistance = inner_prod(vector_points, rNormal);

          // point_projected.Coordinates() = point_to_project.Coordinates() - rNormal * rDistance;

          // rDistance = fabs(rDistance);

          // array_1d<double, 3> point_projected_local_coords;
          // point_projected_local_coords = pCondition->GetGeometry().PointLocalCoordinates(point_projected_local_coords,point_projected);

          // rLocalCoords[0] = point_projected_local_coords[0];
          // rLocalCoords[1] = 0.0;
        

          // if (-1.0-MapperUtilities::tol_local_coords < rLocalCoords[0] && 
          //     rLocalCoords[0] < 1.0+MapperUtilities::tol_local_coords) {
          //     return true;
          // } else {
          //     return false;
          // }
          
          //**********************************************************************
          // xi,yi are Nodal Coordinates, n is the destination condition's unit normal
          // and d is the distance along n from the point to its projection in the condition
          // | DestX-0.5(x1+x2) |   | 0.5(x2-x1)  nx |   | Chi |
          // |                  | = |                | . |     |    
          // | DestY-0.5(y1+y2) |   | 0.5(y2-y1)  ny |   |  d  |

          Matrix transform_matrix(2, 2, false);
          Matrix inv_transform_matrix(2, 2, false);
          double det;

          array_1d<double, 2> RHS;
          array_1d<double, 2> result_vec;
          

          RHS[0] = GlobalCoords[0] - 0.5 * (pCondition->GetGeometry()[0].X() + pCondition->GetGeometry()[1].X());
          RHS[1] = GlobalCoords[1] - 0.5 * (pCondition->GetGeometry()[0].Y() + pCondition->GetGeometry()[1].Y());

          transform_matrix(0, 0) = 0.5 * (pCondition->GetGeometry()[1].X() - pCondition->GetGeometry()[0].X());
          transform_matrix(1, 0) = 0.5 * (pCondition->GetGeometry()[1].Y() - pCondition->GetGeometry()[0].Y());
          
          array_1d<double,3> rNormal;
          CalculateLineNormal(pCondition, rNormal);
          transform_matrix(0, 1) = rNormal[0];
          transform_matrix(1, 1) = rNormal[1];

          MathUtils<double>::InvertMatrix2(transform_matrix,inv_transform_matrix,det);
          noalias(result_vec) = prod(inv_transform_matrix, RHS);
          
          rLocalCoords[0] = result_vec[0];
          rLocalCoords[1] = 0.0f;

          rDistance = result_vec[1];

          // std::cout << "Point to project: [" << GlobalCoords[0] << " " << GlobalCoords[1] << " " << GlobalCoords[2] << "], local points: [" <<
          // pCondition->GetGeometry()[0].X() << " " << pCondition->GetGeometry()[0].Y() << " " << pCondition->GetGeometry()[0].Z() << "] ["<<
          // pCondition->GetGeometry()[1].X() << " " << pCondition->GetGeometry()[1].Y() << " " << pCondition->GetGeometry()[1].Z() << "] local coord: " <<
          // rLocalCoords[0]  << std::endl;

          bool is_inside = false;

          if (rLocalCoords[0] >= -1.0f-MapperUtilities::tol_local_coords) {
              if (rLocalCoords[0] <= 1.0f+MapperUtilities::tol_local_coords) {
                  is_inside = true;
              }
          }
          return is_inside;
      }

      static bool ProjectPointToTriangle(Condition* pCondition,
                                         const array_1d<double, 3> GlobalCoords,
                                         array_1d<double,2>& rLocalCoords,
                                         double& rDistance) {
          // xi,yi,zi are Nodal Coordinates, n is the destination condition's unit normal
          // and d is the distance along n from the point to its projection in the condition
          // | DestX-x1 |   | x2-x1  x3-x1  nx |   | Chi |
          // | DestY-y1 | = | y2-y1  y3-y1  ny | . | Eta |
          // | DestZ-z1 |   | z2-z1  z3-z1  nz |   |  d  |

          Matrix transform_matrix(3, 3, false);
          Matrix inv_transform_matrix(3, 3, false);
          double det;

          array_1d<double, 3> RHS;
          array_1d<double, 3> result_vec;
          

          RHS[0] = GlobalCoords[0] - pCondition->GetGeometry()[0].X();
          RHS[1] = GlobalCoords[1] - pCondition->GetGeometry()[0].Y();
          RHS[2] = GlobalCoords[2] - pCondition->GetGeometry()[0].Z();

          transform_matrix(0, 0) = pCondition->GetGeometry()[1].X() - pCondition->GetGeometry()[0].X();
          transform_matrix(1, 0) = pCondition->GetGeometry()[1].Y() - pCondition->GetGeometry()[0].Y();
          transform_matrix(2, 0) = pCondition->GetGeometry()[1].Z() - pCondition->GetGeometry()[0].Z();

          transform_matrix(0, 1) = pCondition->GetGeometry()[2].X() - pCondition->GetGeometry()[0].X();
          transform_matrix(1, 1) = pCondition->GetGeometry()[2].Y() - pCondition->GetGeometry()[0].Y();
          transform_matrix(2, 1) = pCondition->GetGeometry()[2].Z() - pCondition->GetGeometry()[0].Z();
          
          array_1d<double,3> rNormal;
          CalculateTriangleNormal(pCondition, rNormal);
          transform_matrix(0, 2) = rNormal[0];
          transform_matrix(1, 2) = rNormal[1];
          transform_matrix(2, 2) = rNormal[2];

          MathUtils<double>::InvertMatrix3(transform_matrix,inv_transform_matrix,det);
          noalias(result_vec) = prod(inv_transform_matrix, RHS);

          rLocalCoords[0] = result_vec[0];
          rLocalCoords[1] = result_vec[1];
          
          rDistance = result_vec[2];

          bool is_inside = false;

          if (1.0f-rLocalCoords[0]-rLocalCoords[1] >= 0.0f-MapperUtilities::tol_local_coords) {
              if (rLocalCoords[0] >= 0.0f-MapperUtilities::tol_local_coords) {
                  if (rLocalCoords[1] >= 0.0f-MapperUtilities::tol_local_coords) {
                      is_inside = true;
                  }
              }
          }
          return is_inside;
      }

      static bool ProjectPointToQuadrilateral(Condition* pCondition,
                                              const array_1d<double, 3> GlobalCoords,
                                              array_1d<double,2>& rLocalCoords,
                                              double& rDistance) {
          
          // change localcoords to 3 from 2!
          // Condition::GeometryType& r_condition_geometry = pCondition->GetGeometry();
          // bool is_inside = r_condition_geometry.IsInside(GlobalCoords, rLocalCoords);
          
          // if (is_inside) {
          //     // Calculate Distance
          //     array_1d<double, 3> projection_physical_coords;
          //     projection_physical_coords[0] = 0.0f;
          //     projection_physical_coords[1] = 0.0f;
          //     projection_physical_coords[2] = 0.0f;

          //     double shape_function_value = 0.0f;

          //     for (int i = 0; i < static_cast<int>(r_condition_geometry.PointsNumber()); ++i) {
          //         shape_function_value = r_condition_geometry.ShapeFunctionValue(i, rLocalCoords);
          //         projection_physical_coords[0] += shape_function_value * r_condition_geometry[i].X();
          //         projection_physical_coords[1] += shape_function_value * r_condition_geometry[i].Y();
          //         projection_physical_coords[2] += shape_function_value * r_condition_geometry[i].Z();
          //     }

          //     rDistance = sqrt(pow(GlobalCoords[0] - projection_physical_coords[0] , 2) +
          //                      pow(GlobalCoords[1] - projection_physical_coords[1] , 2) +
          //                      pow(GlobalCoords[2] - projection_physical_coords[2] , 2));

          //     KRATOS_WATCH(rDistance)
          // }

          // return is_inside;

          return false;
      }

      static void CalculateLineNormal(Condition* pCondition,
                                      array_1d<double,3>& rNormal) {
          // TODO use the normal calculation in geometry.h once it is available
          rNormal.clear();
          array_1d<double,3> v1,v2;
          v1[0] = pCondition->GetGeometry()[1].X() - pCondition->GetGeometry()[0].X();
          v1[1] = pCondition->GetGeometry()[1].Y() - pCondition->GetGeometry()[0].Y();
          v1[2] = pCondition->GetGeometry()[1].Z() - pCondition->GetGeometry()[0].Z();
          // Assuming plane X-Y in the 2D-case
          v2[0] = 0.0f;
          v2[1] = 0.0f;
          v2[2] = 1.0f;

          // Compute the condition normal
          MathUtils<double>::CrossProduct(rNormal,v1,v2);

          rNormal /= norm_2(rNormal); // normalize the nomal (i.e. length=1)

          // std::cout << "Normal " << rNormal[0] << " " << rNormal[1] << " " << rNormal[2] << std::endl;
      }

      static void CalculateTriangleNormal(Condition* pCondition,
                                          array_1d<double,3>& rNormal) {
          // TODO use the normal calculation in geometry.h once it is available
          rNormal.clear();
          array_1d<double,3> v1,v2;
          v1[0] = pCondition->GetGeometry()[1].X() - pCondition->GetGeometry()[0].X();
          v1[1] = pCondition->GetGeometry()[1].Y() - pCondition->GetGeometry()[0].Y();
          v1[2] = pCondition->GetGeometry()[1].Z() - pCondition->GetGeometry()[0].Z();

          v2[0] = pCondition->GetGeometry()[2].X() - pCondition->GetGeometry()[0].X();
          v2[1] = pCondition->GetGeometry()[2].Y() - pCondition->GetGeometry()[0].Y();
          v2[2] = pCondition->GetGeometry()[2].Z() - pCondition->GetGeometry()[0].Z();


          // Compute the condition normal
          MathUtils<double>::CrossProduct(rNormal,v1,v2);

          rNormal /= norm_2(rNormal); // normalize the nomal (i.e. length=1)

          // std::cout << "Normal " << rNormal[0] << " " << rNormal[1] << " " << rNormal[2] << std::endl;
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
        buffer << "MapperUtilities" ;
        return buffer.str();
      }

      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const {rOStream << "MapperUtilities";}

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

      static constexpr double tol_local_coords = 1E-15; // TODO what to use here?
      // static constexpr double tol_local_coords = std::numeric_limits<double>::epsilon(); // gives Problems if nodes are in the same location


      ///@}
      ///@name Private Operators
      ///@{


      ///@}
      ///@name Private Operations
      ///@{

      /// Default constructor.
      MapperUtilities() { }

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
      MapperUtilities& operator=(MapperUtilities const& rOther);

    //   /// Copy constructor.
    //   MapperUtilities(MapperUtilities const& rOther){}


      ///@}

    }; // Class MapperUtilities

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream,
				    MapperUtilities& rThis)
    {
        return rIStream;
    }

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream,
				    const MapperUtilities& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }
  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_MAPPER_UTILITIES_H_INCLUDED  defined