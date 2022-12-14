//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main author:     Sebastian Ares de Parga Regalado
//
// This file is partly copied from
// "DEMApplication/custom_utilities/omp_dem_search.h"
// and modified to create only once the binary tree

#if !defined(KRATOS_ASSIGN_PERIODIC_MPCS_TO_NEIGHBOURS_H_INCLUDED )
#define  KRATOS_ASSIGN_PERIODIC_MPCS_TO_NEIGHBOURS_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// Include kratos definitions
#include "includes/define.h"


// Application includes
#include "utilities/variable_utils.h"
#include "utilities/rbf_shape_functions_utility.h"


// Project includes
#include "spatial_containers/spatial_search.h"
#include "spatial_containers/bins_dynamic.h"
#include "utilities/builtin_timer.h"

// Configures

// External includes

namespace Kratos {
///@addtogroup Kratos Core
///@{

///@name Kratos Classes
///@{

/**
 * @class AssignPeriodicMPCsToNeighboursUtility
 * @ingroup Kratos Core
 * @brief Assing MPCs to Neighbouring Nodes
 * @details This class provides a method to search for neighbouring nodes of one node
 * within a given radius and assign a multipoint constraint (MPC) using radial basis functions.
 * @author Sebastian Ares de Parga Regalado
 */

class KRATOS_API(KRATOS_CORE) AssignPeriodicMPCsToNeighboursUtility
{
    public:
      ///@name Type Definitions
      ///@{

      /// Pointer definition of AssignPeriodicMPCsToNeighboursUtility
      KRATOS_CLASS_POINTER_DEFINITION(AssignPeriodicMPCsToNeighboursUtility);

      //Node Types
      typedef ModelPart::NodeType                           NodeType;
      typedef Node < 3 > ::Pointer NodeTypePointer;

      //Dof types
      typedef NodeType::DofType DofType;

      //Configure Types
      typedef ModelPart::NodesContainerType                 NodesContainerType;

      //Bin
      typedef BinsDynamic<3, NodeType, ModelPart::NodesContainerType::ContainerType> NodeBinsType;

      /// General containers type definitions
      typedef ModelPart::MasterSlaveConstraintContainerType ConstraintContainerType;
      typedef NodesContainerType::ContainerType             ResultNodesContainerType;
      typedef std::vector<ResultNodesContainerType>             VectorResultNodesContainerType;
      typedef std::vector<double>                               RadiusArrayType;
      typedef std::vector<NodeTypePointer> NodeVector;

      typedef std::vector< Dof<double>::Pointer > DofPointerVectorType;

      ///@}
      ///@name Life Cycle
      ///@{

      /// Default constructor.
      AssignPeriodicMPCsToNeighboursUtility(
        ModelPart& rMasterModelPart,
        ModelPart& rAuxiliarModelPart,
        const int NumberOfSlices,
        const double SymmetryAngle,
        IndexType SymmetryNodeId
      ):mrMasterModelPart(rMasterModelPart),
      mrAuxiliarModelPart(rAuxiliarModelPart){
        KRATOS_TRY;
        int aux_node_counter = 1;
        for(int j = 0; j < static_cast<int>(NumberOfSlices); j++){
            double alpha = j*SymmetryAngle;
            Matrix rotation_matrix(3,3);
            rotation_matrix(0,0) = cos(alpha);
            rotation_matrix(0,1) = -sin(alpha);
            rotation_matrix(0,2) = 0.0;
            rotation_matrix(1,0) = sin(alpha);
            rotation_matrix(1,1) = cos(alpha);
            rotation_matrix(1,2) = 0.0;
            rotation_matrix(2,0) = 0.0;
            rotation_matrix(2,1) = 0.0;
            rotation_matrix(2,2) = 1;
            NodesContainerType &r_nodes_master   = mrMasterModelPart.Nodes();

            for(NodeType& node_i : r_nodes_master){
              if (node_i.Id()!=SymmetryNodeId){
                    const auto& master_coordinates = node_i.Coordinates();
                    Vector rotated_coordinates(master_coordinates.size());
                    noalias(rotated_coordinates) = prod(rotation_matrix,master_coordinates);
                    mrAuxiliarModelPart.CreateNewNode(aux_node_counter, rotated_coordinates[0], rotated_coordinates[1], rotated_coordinates[2]);
                    mrAuxiliarModelPart.GetNode(aux_node_counter).SetValue(AUX_INDEX, node_i.Id());
                    mrAuxiliarModelPart.GetNode(aux_node_counter).SetValue(AUX_MESH_VAR, alpha);
                  }
                  aux_node_counter += 1;
            }
        }

        NodesContainerType& rStructureNodes = mrAuxiliarModelPart.Nodes();
        NodesContainerType::ContainerType& nodes_model_part = rStructureNodes.GetContainer();
        mpBins = Kratos::make_unique<NodeBinsType>(nodes_model_part.begin(), nodes_model_part.end());
        mMaxNumberOfNodes = rStructureNodes.size();
        KRATOS_CATCH("");
      }

      /// Destructor.
      ~AssignPeriodicMPCsToNeighboursUtility(){
      }

      ///@}
      ///@name Operations
      ///@{

      /**
       * @brief Perform Node Search
       * @details Searches for nodes within a given 'Radius' of the current nodes
       * @param rStructureNodes Nodes Container.
       * @param rNodes Nodes to obtain their respective neighbouring nodes.
       * @param Radius Search radius.
       * @param rResults Results container.
       **/

      // Search for the neighbouring nodes (in rStructureNodes) of each rNode on a given array of rNodes
      void SearchNodesInRadiusForNodes(
          NodesContainerType const& rNodes,
          const double Radius,
          VectorResultNodesContainerType& rResults)
      {
          KRATOS_TRY;
          NodesContainerType::ContainerType& nodes_array     = const_cast<NodesContainerType::ContainerType&>(rNodes.GetContainer());

          rResults.resize(nodes_array.size());

          #pragma omp parallel
          {
              ResultNodesContainerType  localResults(mMaxNumberOfNodes);
              std::size_t               NumberOfResults = 0;

              #pragma omp for
              for(int i = 0; i < static_cast<int>(nodes_array.size()); i++)
              {
                  ResultNodesContainerType::iterator ResultsPointer    = localResults.begin();

                  NumberOfResults = mpBins->SearchInRadius(*(nodes_array[i]),Radius,ResultsPointer,mMaxNumberOfNodes);

                  rResults[i].insert(rResults[i].begin(),localResults.begin(),localResults.begin()+NumberOfResults);
              }
          }
          KRATOS_CATCH("");
      }

      /**
       * @brief Perform Node Search
       * @details Searches for nodes within a given 'Radius' of the current node
       * @param rStructureNodes Nodes Container.
       * @param rNode Node to obtain its respective neighbouring nodes.
       * @param Radius Search radius.
       * @param rResults Results container.
       **/
      //Search for the neighbouring nodes (in rStructureNodes) of the given rNode
      void SearchNodesInRadiusForNode(
          NodeType::Pointer pNode,
          double const Radius,
          ResultNodesContainerType& rResults)
      {

          KRATOS_TRY;

          ResultNodesContainerType local_results(mMaxNumberOfNodes);
          std::size_t num_of_results = 0;

          auto results_iterator = local_results.begin();

          num_of_results = mpBins->SearchInRadius(*pNode, Radius, results_iterator, mMaxNumberOfNodes);

          rResults.insert(rResults.begin(), local_results.begin(), local_results.begin() + num_of_results);

          KRATOS_CATCH("");
      }

      // Function caller and results initializer for SearchNodesInRadiusForNodes
      VectorResultNodesContainerType SearchCloudOfNodesForNodes(
          NodesContainerType pNodes,
          double const Radius)
      {
          VectorResultNodesContainerType  Results;
          KRATOS_TRY;
          SearchNodesInRadiusForNodes(pNodes, Radius, Results);
          KRATOS_CATCH("");
          return Results;
      }


      /**
       * @brief Collect Dofs and Coordinates
       * @details Collects Dofs and Coordinates for a given node or set of nodes.
       * @param rStructureNodes Nodes Container.
       * @param pNode Obtain respective dofs and coordinates for a given node or set of nodes.
       * @param rVariable Dof variable array or double.
       * @param rCloud_of_dofs Dofs container.
       * @param rCloud_of_nodes_coordinates Coordinates container.
       **/
      void GetDofsAndCoordinatesForNodes(
        NodeType::Pointer pNode,
        const Variable<array_1d<double, 3>>& rVariable,
        DofPointerVectorType& rCloud_of_dofs_x,
        DofPointerVectorType& rCloud_of_dofs_y,
        DofPointerVectorType& rCloud_of_dofs_z,
        array_1d<double,3>& rSlave_coordinates
      ){
        KRATOS_TRY;
        if (rVariable==VELOCITY){
          rCloud_of_dofs_x.push_back(pNode->pGetDof(VELOCITY_X));
          rCloud_of_dofs_y.push_back(pNode->pGetDof(VELOCITY_Y));
          rCloud_of_dofs_z.push_back(pNode->pGetDof(VELOCITY_Z));
        }
        else if (rVariable==DISPLACEMENT){
          // NodesContainerType::ContainerType& nodes_array     = const_cast<NodesContainerType::ContainerType&>(pNodes.GetContainer());
          rCloud_of_dofs_x.push_back(pNode->pGetDof(DISPLACEMENT_X));
          rCloud_of_dofs_y.push_back(pNode->pGetDof(DISPLACEMENT_Y));
          rCloud_of_dofs_z.push_back(pNode->pGetDof(DISPLACEMENT_Z));
        }
        else{
          KRATOS_ERROR << "This function is not yet implemented for variable " << rVariable << std::endl;
        }
        rSlave_coordinates = pNode->Coordinates();
        KRATOS_CATCH("");
      }


      // Get Dofs and coordinates arrays for a given variable array. (For nodes)
      void GetDofsAndCoordinatesForNodes(
        ResultNodesContainerType nodes_array,
        const Variable<array_1d<double, 3>>& rVariable,
        DofPointerVectorType& rCloud_of_dofs_x,
        DofPointerVectorType& rCloud_of_dofs_y,
        DofPointerVectorType& rCloud_of_dofs_z,
        Matrix& rCloud_of_nodes_coordinates
        )
      {
        KRATOS_TRY;
        rCloud_of_dofs_x.resize(nodes_array.size());
        rCloud_of_dofs_y.resize(nodes_array.size());
        rCloud_of_dofs_z.resize(nodes_array.size());
        rCloud_of_nodes_coordinates.resize(nodes_array.size(),3);
        if (rVariable==VELOCITY){
          for(int i = 0; i < static_cast<int>(nodes_array.size()); i++)
          {
            rCloud_of_dofs_x[i] = nodes_array[i]->pGetDof(VELOCITY_X);
            rCloud_of_dofs_y[i] = nodes_array[i]->pGetDof(VELOCITY_Y);
            rCloud_of_dofs_z[i] = nodes_array[i]->pGetDof(VELOCITY_Z);
            noalias(row(rCloud_of_nodes_coordinates,i)) = nodes_array[i]->Coordinates();
          }
        }
        else if (rVariable==DISPLACEMENT){
          for(int i = 0; i < static_cast<int>(nodes_array.size()); i++)
          {
            rCloud_of_dofs_x[i] = nodes_array[i]->pGetDof(DISPLACEMENT_X);
            rCloud_of_dofs_y[i] = nodes_array[i]->pGetDof(DISPLACEMENT_Y);
            rCloud_of_dofs_z[i] = nodes_array[i]->pGetDof(DISPLACEMENT_Z);
            noalias(row(rCloud_of_nodes_coordinates,i)) = nodes_array[i]->Coordinates();
          }
        }
        else{
          KRATOS_ERROR << "This function is not yet implemented for variable " << rVariable << std::endl;
        }
        KRATOS_CATCH("");
      }

      void GetDofsAndCoordinatesForNodes(
        NodeType::Pointer pNode,
        const Variable<double>& rVariable,
        DofPointerVectorType& rCloud_of_dofs,
        array_1d<double,3>& rSlave_coordinates
      ){
        KRATOS_TRY;
        rCloud_of_dofs.push_back(pNode->pGetDof(rVariable));
        rSlave_coordinates = pNode->Coordinates();
        KRATOS_CATCH("");
      }

      // Get Dofs and Coordinates arrays for a given variable double. (For nodes)
      void GetDofsAndCoordinatesForNodes(
        ResultNodesContainerType nodes_array,
        const Variable<double>& rVariable,
        DofPointerVectorType& rCloud_of_dofs,
        Matrix& rCloud_of_nodes_coordinates
        )
      {
        KRATOS_TRY;
        rCloud_of_dofs.resize(nodes_array.size());
        rCloud_of_nodes_coordinates.resize(nodes_array.size(),3);
        for(int i = 0; i < static_cast<int>(nodes_array.size()); i++)
        {
          rCloud_of_dofs[i] = nodes_array[i]->pGetDof(rVariable);
          noalias(row(rCloud_of_nodes_coordinates,i)) = nodes_array[i]->Coordinates();
        }
        KRATOS_CATCH("");
      }

      void GetCoordinatesForNodes(
        ResultNodesContainerType nodes_array,
        Matrix& rCloud_of_nodes_coordinates
        )
      {
        KRATOS_TRY;
        rCloud_of_nodes_coordinates.resize(nodes_array.size(),3);
        for(int i = 0; i < static_cast<int>(nodes_array.size()); i++)
        {
          noalias(row(rCloud_of_nodes_coordinates,i)) = nodes_array[i]->Coordinates();
        }
        KRATOS_CATCH("");
      }

      void GetDofsForNodes(
        ResultNodesContainerType nodes_array,
        const Variable<array_1d<double, 3>>& rVariable,
        DofPointerVectorType& rCloud_of_dofs
      ){
        KRATOS_TRY;
        rCloud_of_dofs.resize(nodes_array.size()*3);
        if (rVariable==VELOCITY){
          for(int i = 0; i < static_cast<int>(nodes_array.size()); i++)
          {
            const int local_index = i*3;
            rCloud_of_dofs[local_index] = nodes_array[i]->pGetDof(VELOCITY_X);
            rCloud_of_dofs[local_index+1] = nodes_array[i]->pGetDof(VELOCITY_Y);
            rCloud_of_dofs[local_index+2] = nodes_array[i]->pGetDof(VELOCITY_Z);
          }
        }
        KRATOS_CATCH("");
      }

      void AssignRotationToNodes(
        NodesContainerType pNodes,
        Matrix RotationMatrix
        )
      {
        //TODO: Do it in parallel
        KRATOS_TRY;
        NodesContainerType::ContainerType& nodes_array     = const_cast<NodesContainerType::ContainerType&>(pNodes.GetContainer());

        // Assign rotation to given nodes
        for(int i = 0; i < static_cast<int>(nodes_array.size()); i++)
        {
          auto coordinate_vector = nodes_array[i]->GetInitialPosition().Coordinates();
          Vector rotated_position(coordinate_vector.size());
          noalias(rotated_position) = prod(RotationMatrix,coordinate_vector);
          nodes_array[i]->X() = rotated_position[0];
          nodes_array[i]->Y() = rotated_position[1];
          nodes_array[i]->Z() = rotated_position[2];
        }
        KRATOS_CATCH("");
      }

      void GetMasterNodesAndAngles(
        ResultNodesContainerType Results,
        NodeVector& ResultsMaster,
        Vector& ListOfAngles
      ){
        ResultsMaster.reserve(Results.size());
        ListOfAngles.resize(Results.size());
        for(int i = 0; i < static_cast<int>(Results.size()); i++)
              {
                double master_id = Results[i]->GetValue(AUX_INDEX);
                double alpha = Results[i]->GetValue(AUX_MESH_VAR);
                ResultsMaster.push_back(mrMasterModelPart.pGetNode(int(master_id)));
                ListOfAngles[i] = alpha;
              }
      }

      void ObtainShapeMatrix(
          Matrix& shape_matrix_x,
          Matrix& shape_matrix_y,
          Matrix& shape_matrix_z,
          Vector& rN_container,
          Vector& ListOfAngles
      )
      {
          std::vector<double> rotation_x_global;
          std::vector<double> rotation_y_global;
          std::vector<double> rotation_z_global;
          std::vector<double> rotation_x(3);
          std::vector<double> rotation_y(3);
          std::vector<double> rotation_z(3);
          for (int i=0; i < static_cast<int>(rN_container.size()); i++){
              //TODO: Any kind of rotation. x-axis, y-axis, z-axis. Make it a function.
              rotation_x[0] = rN_container[i]*cos(ListOfAngles[i]);
              rotation_x[1] = -rN_container[i]*sin(ListOfAngles[i]);
              rotation_x[2] = rN_container[i]*0.0;
              rotation_y[0] = rN_container[i]*sin(ListOfAngles[i]);
              rotation_y[1] = rN_container[i]*cos(ListOfAngles[i]);
              rotation_y[2] = rN_container[i]*0.0;
              rotation_z[0] = rN_container[i]*0.0;
              rotation_z[1] = rN_container[i]*0.0;
              rotation_z[2] = rN_container[i]*1.0;
              // noalias(row(shape_matrix,i)) = aux_vector_x;
              rotation_x_global.insert(rotation_x_global.end(),rotation_x.begin(),rotation_x.end());
              rotation_y_global.insert(rotation_y_global.end(),rotation_y.begin(),rotation_y.end());
              rotation_z_global.insert(rotation_z_global.end(),rotation_z.begin(),rotation_z.end());
          }
          shape_matrix_x.resize(1, rotation_x_global.size());
          shape_matrix_y.resize(1, rotation_y_global.size());
          shape_matrix_z.resize(1, rotation_z_global.size());
          for (int i=0; i<static_cast<int>(rotation_x_global.size());i++){
              shape_matrix_x(0,i) = rotation_x_global[i];
              shape_matrix_y(0,i) = rotation_y_global[i];
              shape_matrix_z(0,i) = rotation_z_global[i];
          }
      }


      /**
       * @brief Assign Multipoint Constraints to a set of Nodes.
       * @details Assign Multipoint Constraints to a set of Nodes given a radius of influence
       * w.r.t MLS or RBF shape functions.
       * @param rStructureNodes Nodes Container.
       * @param pNodes Nodes to set MPCs.
       * @param Radius Search radius.
       * @param rComputinModelPart Model Part to which MPCs are applied.
       * @param rVariable DOFs to assign the MPCs.
       * @param h Shape parameter (to scale the input of the radial kernel)
       **/
      void AssignMPCsToNodes(
          NodesContainerType pNodes,
          double const Radius,
          ModelPart& rComputingModelPart,
          const Variable<array_1d<double, 3>>& rVariable
          )
          //TODO: Do it in parallel and allow rVariable to be a double variable i.e. TEMPERATURE (Make a template)
      {
          KRATOS_TRY;
          // #pragma omp parallel
          {
            BuiltinTimer build_and_assign_mpcs;

            NodesContainerType::ContainerType& nodes_array     = const_cast<NodesContainerType::ContainerType&>(pNodes.GetContainer());

            ModelPart::MasterSlaveConstraintType const& r_clone_constraint = KratosComponents<MasterSlaveConstraint>::Get("LinearMasterSlaveConstraint");
            ConstraintContainerType all_constraints;
            all_constraints.reserve(nodes_array.size()*3);
            int prev_num_mpcs = rComputingModelPart.NumberOfMasterSlaveConstraints();

            // #pragma omp for
              for(int i = 0; i < static_cast<int>(nodes_array.size()); i++)
              {
                //Set slave nodes' slip condition to false
                // nodes_array[i]->Set(SLIP,false);

                // Search cloud of nodes for slave node
                ResultNodesContainerType  Results;
                double localRadius = Radius;
                while (Results.size()<3){
                  Results.clear();
                  SearchNodesInRadiusForNode(nodes_array[i], localRadius, Results);
                  localRadius += Radius;
                }

                NodeVector  ResultsMaster;
                Vector ListOfAngles;
                GetMasterNodesAndAngles(Results, ResultsMaster, ListOfAngles);

                // Get Dofs and Coordinates
                DofPointerVectorType rCloud_of_dofs, rCloud_of_dofs_x,rCloud_of_dofs_y,rCloud_of_dofs_z,rSlave_dof_x,rSlave_dof_y,rSlave_dof_z;
                Matrix rCloud_of_nodes_coordinates_master;
                array_1d<double,3> rSlave_coordinates;
                GetDofsAndCoordinatesForNodes(nodes_array[i], rVariable, rSlave_dof_x,rSlave_dof_y,rSlave_dof_z, rSlave_coordinates);
                GetDofsAndCoordinatesForNodes(ResultsMaster, rVariable, rCloud_of_dofs_x, rCloud_of_dofs_y, rCloud_of_dofs_z, rCloud_of_nodes_coordinates_master);
                GetDofsForNodes(ResultsMaster, rVariable, rCloud_of_dofs);
                if (nodes_array[i]->Id()==1111){
                  KRATOS_WATCH(ListOfAngles)
                }

                Matrix rCloud_of_nodes_coordinates_aux_master;
                GetCoordinatesForNodes(Results, rCloud_of_nodes_coordinates_aux_master);

                // Calculate shape functions
                Vector rN_container;
                RBFShapeFunctionsUtility::CalculateShapeFunctions(rCloud_of_nodes_coordinates_aux_master,rSlave_coordinates,rN_container);

                //Create MPCs
                Matrix shape_matrix_x, shape_matrix_y, shape_matrix_z;
                ObtainShapeMatrix(shape_matrix_x, shape_matrix_y, shape_matrix_z, rN_container, ListOfAngles);
                Matrix shape_matrix(1,rN_container.size());
                noalias(row(shape_matrix,0)) = rN_container;// Shape functions matrix
                const Vector constant_vector = ZeroVector(rN_container.size());
                IndexType it = i*3;
                all_constraints.push_back(r_clone_constraint.Create(prev_num_mpcs + it + 1,rCloud_of_dofs,rSlave_dof_x,shape_matrix_x,constant_vector));
                all_constraints.push_back(r_clone_constraint.Create(prev_num_mpcs + it + 2,rCloud_of_dofs,rSlave_dof_y,shape_matrix_y,constant_vector));
                all_constraints.push_back(r_clone_constraint.Create(prev_num_mpcs + it + 3,rCloud_of_dofs,rSlave_dof_z,shape_matrix_z,constant_vector));
              }
              rComputingModelPart.AddMasterSlaveConstraints(all_constraints.begin(),all_constraints.end());
              KRATOS_INFO("AssignPeriodicMPCsToNeighboursUtility") << "Build and Assign MPCs Time: " << build_and_assign_mpcs.ElapsedSeconds() << std::endl;
          }
          KRATOS_CATCH("");
      }

      void AssignMPCsToNodes(
          NodesContainerType pNodes,
          double const Radius,
          ModelPart& rComputingModelPart,
          const Variable<double>& rVariable
          )
          //TODO: Do it in parallel and allow rVariable to be a double variable i.e. TEMPERATURE (Make a template)
      {
          KRATOS_TRY;
          // #pragma omp parallel
          {
            BuiltinTimer build_and_assign_mpcs;

            NodesContainerType::ContainerType& nodes_array     = const_cast<NodesContainerType::ContainerType&>(pNodes.GetContainer());

            ModelPart::MasterSlaveConstraintType const& r_clone_constraint = KratosComponents<MasterSlaveConstraint>::Get("LinearMasterSlaveConstraint");
            ConstraintContainerType all_constraints;
            all_constraints.reserve(nodes_array.size()*3);
            int prev_num_mpcs = rComputingModelPart.NumberOfMasterSlaveConstraints();

            // #pragma omp for
              for(int i = 0; i < static_cast<int>(nodes_array.size()); i++)
              {
                //Set slave nodes' slip condition to false
                // nodes_array[i]->Set(SLIP,false);

                // Search cloud of nodes for slave node
                ResultNodesContainerType  Results;
                double localRadius = Radius;
                while (Results.size()<3){
                  Results.clear();
                  SearchNodesInRadiusForNode(nodes_array[i], localRadius, Results);
                  localRadius += Radius;
                }

                // Get Dofs and Coordinates
                DofPointerVectorType rCloud_of_dofs,rSlave_dof;
                Matrix rCloud_of_nodes_coordinates;
                array_1d<double,3> rSlave_coordinates;
                GetDofsAndCoordinatesForNodes(nodes_array[i], rVariable, rSlave_dof, rSlave_coordinates);
                GetDofsAndCoordinatesForNodes(Results, rVariable, rCloud_of_dofs, rCloud_of_nodes_coordinates);

                // Calculate shape functions
                Vector rN_container;
                RBFShapeFunctionsUtility::CalculateShapeFunctions(rCloud_of_nodes_coordinates,rSlave_coordinates,rN_container);

                //Create MPCs
                Matrix shape_matrix(1,rN_container.size());
                noalias(row(shape_matrix,0)) = rN_container;// Shape functions matrix
                const Vector constant_vector = ZeroVector(rN_container.size());
                IndexType it = i;
                all_constraints.push_back(r_clone_constraint.Create(prev_num_mpcs + it + 1,rCloud_of_dofs,rSlave_dof,shape_matrix,constant_vector));
              }
              rComputingModelPart.AddMasterSlaveConstraints(all_constraints.begin(),all_constraints.end());
              KRATOS_INFO("AssignPeriodicMPCsToNeighboursUtility") << "Build and Assign MPCs Time: " << build_and_assign_mpcs.ElapsedSeconds() << std::endl;
          }
          KRATOS_CATCH("");
      }



      ///@}
      ///@name Input and output
      ///@{

      /// Turn back information as a string.
      virtual std::string Info() const
      {
          std::stringstream buffer;
          buffer << "AssignPeriodicMPCsToNeighboursUtility" ;

          return buffer.str();
      }

      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const  {rOStream << "AssignPeriodicMPCsToNeighboursUtility";}

      /// Print object's data.
      virtual void PrintData(std::ostream& rOStream) const  {}

      ///@}

    private:
      ///@name Member Variables
      ///@{
      Kratos::unique_ptr<NodeBinsType> mpBins;
      int mMaxNumberOfNodes;
      ModelPart& mrMasterModelPart;
      ModelPart& mrAuxiliarModelPart;
      // ModelPart& mModelPart;


      ///@}
      ///@name Un accessible methods
      ///@{

      /// Assignment operator.
      AssignPeriodicMPCsToNeighboursUtility& operator=(AssignPeriodicMPCsToNeighboursUtility const& rOther)
      {
          return *this;
      }

      /// Copy constructor.
      // AssignPeriodicMPCsToNeighboursUtility(AssignPeriodicMPCsToNeighboursUtility const& rOther)
      // {
      //     *this = rOther;
      // }

      ///@}

    }; // Class AssignPeriodicMPCsToNeighboursUtility



///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_ASSIGN_PERIODIC_MPCS_TO_NEIGHBOURS_H_INCLUDED  defined