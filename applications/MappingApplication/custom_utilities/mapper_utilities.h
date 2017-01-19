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
              std::cout << "MAPPER WARNING, no conditions for search radius "
                        << "computations found, using nodes (less efficient)"
                        << std::endl;
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

      static double ComputeConservativeFactor(const double num_nodes_origin,
                                              const double num_nodes_destination) {
          // num_nodes_* are casted to doubles in order to use the double devision
          // if this function would take ints, then the return value would also be an int!
          return num_nodes_origin / num_nodes_destination;
      }

      static bool ProjectPointToLine(Condition* p_condition,
                                     array_1d<double, 3> global_coords,
                                     array_1d<double,2>& local_coords,
                                     double& distance,
                                     array_1d<double,3>& normal) {
          Point<3> point_to_project(global_coords);

          Point<3> point_projected;
          Point<3> point_in_plane = p_condition->GetGeometry()[0];

          array_1d<double,3> vector_points;
          noalias(vector_points) = point_to_project.Coordinates() - point_in_plane.Coordinates();

          CalculateLineNormal(p_condition, normal);

          distance = inner_prod(vector_points, normal);

          point_projected.Coordinates() = point_to_project.Coordinates() - normal * distance;

          distance = fabs(distance);

          array_1d<double, 3> point_projected_local_coords;
          point_projected_local_coords = p_condition->GetGeometry().PointLocalCoordinates(point_projected_local_coords,point_projected);

          local_coords[0] = point_projected_local_coords[0];
          local_coords[1] = 0.0;

          std::cout << "Point to project: [" << global_coords[0] << " " << global_coords[1] << " " << global_coords[2] << "], local points: [" <<
          p_condition->GetGeometry()[0].X() << " " << p_condition->GetGeometry()[0].Y() << " " << p_condition->GetGeometry()[0].Z() << "] ["<<
          p_condition->GetGeometry()[1].X() << " " << p_condition->GetGeometry()[1].Y() << " " << p_condition->GetGeometry()[1].Z() << "] local coord: " <<
          local_coords[0] << std::endl;


          if (-1.0-MapperUtilities::tol_local_coords < local_coords[0] && local_coords[0] < 1.0+MapperUtilities::tol_local_coords) {
              return true;
          } else {
              return false;
          }

      }

      static bool ProjectPointToTriangle(Condition* p_condition,
                                         array_1d<double, 3> global_coords,
                                         array_1d<double,2>& local_coords,
                                         double& distance,
                                         array_1d<double,3>& normal) {
          return false;
      }

      static bool ProjectPointToQuaddrilateral(Condition* p_condition,
                                               array_1d<double, 3> global_coords,
                                               array_1d<double,2>& local_coords,
                                               double& distance,
                                               array_1d<double,3>& normal) {
          return false;
      }

      static void CalculateLineNormal(Condition* p_condition,
                                      array_1d<double,3>& normal) {
          // TODO use the normal calculation in geometry.h once it is available
          array_1d<double,3> v1,v2;
          v1[0] = p_condition->GetGeometry()[1].X() - p_condition->GetGeometry()[0].X();
          v1[1] = p_condition->GetGeometry()[1].Y() - p_condition->GetGeometry()[0].Y();
          v1[2] = p_condition->GetGeometry()[1].Z() - p_condition->GetGeometry()[0].Z();
          // Assuming plane X-Y in the 2D-case
          v2[0] = 0.0f;
          v2[1] = 0.0f;
          v2[2] = 1.0f;

          // Compute the condition normal
          MathUtils<double>::CrossProduct(normal,v1,v2);

          normal /= norm_2(normal); // normalize the nomal (i.e. length=1)

          // std::cout << "Normal " << normal[0] << " " << normal[1] << " " << normal[2] << std::endl;
      }

      static void CalculateTriangleNormal(Condition* p_condition,
                                          array_1d<double,3>& normal) {
          // TODO use the normal calculation in geometry.h once it is available
          array_1d<double,3> v1,v2;
          v1[0] = p_condition->GetGeometry()[1].X() - p_condition->GetGeometry()[0].X();
          v1[1] = p_condition->GetGeometry()[1].Y() - p_condition->GetGeometry()[0].Y();
          v1[2] = p_condition->GetGeometry()[1].Z() - p_condition->GetGeometry()[0].Z();

          v2[0] = p_condition->GetGeometry()[2].X() - p_condition->GetGeometry()[0].X();
          v2[1] = p_condition->GetGeometry()[2].Y() - p_condition->GetGeometry()[0].Y();
          v2[2] = p_condition->GetGeometry()[2].Z() - p_condition->GetGeometry()[0].Z();


          // Compute the condition normal
          MathUtils<double>::CrossProduct(normal,v1,v2);

          normal /= norm_2(normal); // normalize the nomal (i.e. length=1)

          // std::cout << "Normal " << normal[0] << " " << normal[1] << " " << normal[2] << std::endl;
      }

      static void CalculateQuadrilateralNormal() {

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

      static constexpr double tol_local_coords = 1e-4;


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
