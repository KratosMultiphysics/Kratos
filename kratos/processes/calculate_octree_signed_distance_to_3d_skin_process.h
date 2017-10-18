//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//



#if !defined(KRATOS_CALCULATE_DISTANCE_PROCESS_H_INCLUDED )
#define  KRATOS_CALCULATE_DISTANCE_PROCESS_H_INCLUDED



// System includes
#include <string>
#include <iostream> 


// External includes 


// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "includes/model_part.h"
#include "spatial_containers/octree_binary.h"
#include "utilities/spatial_containers_configure.h"
#include "utilities/timer.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h"


namespace Kratos
{


class DistanceSpatialContainersConfigure
{
public:

    class CellNodeData
    {
        double mDistance;
        double mCoordinates[3];
        std::size_t mId;
    public:
        double& Distance(){return mDistance;}
        double& X() {return mCoordinates[0];}
        double& Y() {return mCoordinates[1];}
        double& Z() {return mCoordinates[2];}
        double& Coordinate(int i) {return mCoordinates[i-1];}
        std::size_t& Id(){return mId;}
    };


    ///@name Type Definitions
    ///@{

    enum { Dimension = 3,
           DIMENSION = 3,
           MAX_LEVEL = 12,
           MIN_LEVEL = 2
         };
    typedef Point                                           PointType;  /// always the point 3D
    typedef std::vector<double>::iterator                   DistanceIteratorType;
    typedef ModelPart::ElementsContainerType::ContainerType ContainerType;
    typedef ContainerType::value_type                       PointerType;
    typedef ContainerType::iterator                         IteratorType;
    typedef ModelPart::ElementsContainerType::ContainerType ResultContainerType;
    typedef ResultContainerType::value_type                 ResultPointerType;
    typedef ResultContainerType::iterator                   ResultIteratorType;

    typedef Element::Pointer                                        pointer_type;
    typedef std::vector<CellNodeData*> data_type;



    typedef  std::vector<PointerType>::iterator             PointerTypeIterator;




    /// Pointer definition of DistanceSpatialContainersConfigure
    KRATOS_CLASS_POINTER_DEFINITION(DistanceSpatialContainersConfigure);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    DistanceSpatialContainersConfigure() {}

    /// Destructor.
    virtual ~DistanceSpatialContainersConfigure() {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    static data_type* AllocateData() {
        return new data_type(27, NULL);
    }

    static void CopyData(data_type* source, data_type* destination) {
        destination = source;
    }

    static void DeleteData(data_type* data) {
        delete data;
    }

    static inline void CalculateBoundingBox(const PointerType& rObject, PointType& rLowPoint, PointType& rHighPoint)
    {
        rHighPoint = rObject->GetGeometry().GetPoint(0);
        rLowPoint  = rObject->GetGeometry().GetPoint(0);
        for (unsigned int point = 0; point<rObject->GetGeometry().PointsNumber(); point++)
        {
            for(std::size_t i = 0; i<3; i++)
            {
                rLowPoint[i]  =  (rLowPoint[i]  >  rObject->GetGeometry().GetPoint(point)[i] ) ?  rObject->GetGeometry().GetPoint(point)[i] : rLowPoint[i];
                rHighPoint[i] =  (rHighPoint[i] <  rObject->GetGeometry().GetPoint(point)[i] ) ?  rObject->GetGeometry().GetPoint(point)[i] : rHighPoint[i];
            }
        }
    }

    static inline void GetBoundingBox(const PointerType rObject, double* rLowPoint, double* rHighPoint)
    {

        for(std::size_t i = 0; i<3; i++)
        {
            rLowPoint[i]  =  rObject->GetGeometry().GetPoint(0)[i];
            rHighPoint[i] =  rObject->GetGeometry().GetPoint(0)[i];
        }

        for (unsigned int point = 0; point<rObject->GetGeometry().PointsNumber(); point++)
        {
            for(std::size_t i = 0; i<3; i++)
            {
                rLowPoint[i]  =  (rLowPoint[i]  >  rObject->GetGeometry().GetPoint(point)[i] ) ?  rObject->GetGeometry().GetPoint(point)[i] : rLowPoint[i];
                rHighPoint[i] =  (rHighPoint[i] <  rObject->GetGeometry().GetPoint(point)[i] ) ?  rObject->GetGeometry().GetPoint(point)[i] : rHighPoint[i];
            }
        }
    }

    static inline bool Intersection(const PointerType& rObj_1, const PointerType& rObj_2)
    {
        Element::GeometryType& geom_1 = rObj_1->GetGeometry();
        Element::GeometryType& geom_2 = rObj_2->GetGeometry();
        return  geom_1.HasIntersection(geom_2);

    }

    static inline bool IntersectionBox(const PointerType& rObject,  const PointType& rLowPoint, const PointType& rHighPoint)
    {
        return rObject->GetGeometry().HasIntersection(rLowPoint, rHighPoint);
    }
    
    static inline bool IsIntersected(const PointerType& rObject, const double& tolerance, const double rLowPoint[], const double rHighPoint[])
    {
        Kratos::Element::GeometryType& geom_1 = rObject->GetGeometry();
        
        Kratos::Point rLowPointTolerance;
        Kratos::Point rHighPointTolerance;
        
        for(std::size_t i = 0; i<3; i++)
        {
            rLowPointTolerance[i]  =  rLowPoint[i] * 1+tolerance;
            rHighPointTolerance[i] =  rHighPoint[i] * 1+tolerance;
        }
        
        return  geom_1.HasIntersection(rLowPointTolerance,rHighPointTolerance);
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
        return " Spatial Containers Configure";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const {}

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const {}


    ///@}

protected:

private:

    /// Assignment operator.
    DistanceSpatialContainersConfigure& operator=(DistanceSpatialContainersConfigure const& rOther);

    /// Copy constructor.
    DistanceSpatialContainersConfigure(DistanceSpatialContainersConfigure const& rOther);


}; // Class DistanceSpatialContainersConfigure

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
  class CalculateSignedDistanceTo3DSkinProcess
	: public Process
    {
    public:
      ///@name Type Definitions
      ///@{
      
      /// Pointer definition of CalculateSignedDistanceTo3DSkinProcess
        KRATOS_CLASS_POINTER_DEFINITION(CalculateSignedDistanceTo3DSkinProcess);
        
        typedef DistanceSpatialContainersConfigure ConfigurationType;
        typedef OctreeBinaryCell<ConfigurationType> CellType;
        typedef OctreeBinary<CellType> OctreeType;

      ///@}
      ///@name Life Cycle 
      ///@{ 
      
      /// Constructor.
      CalculateSignedDistanceTo3DSkinProcess(ModelPart& rThisModelPart)
		: mrSkinModelPart(rThisModelPart), mrBodyModelPart(rThisModelPart)
	{
	}

      /// Destructor.
      virtual ~CalculateSignedDistanceTo3DSkinProcess()
	{
	}
      

      ///@}
      ///@name Operators 
      ///@{

      void operator()()
	{
	  Execute();
	}
      
      
      ///@}
      ///@name Operations
      ///@{

      virtual void Execute()
      {
          KRATOS_TRY

	    //std::cout << "Generating the Octree..." << std::endl;
          GenerateOctree();
          //std::cout << "Generating the Octree finished" << std::endl;

          GenerateCellNodalData();

          CalculateDistance();
          CalculateDistance2();

          std::ofstream mesh_file1("octree1.post.msh");
          std::ofstream res_file("octree1.post.res");

          Timer::Start("Writing Gid conform Mesh");

          PrintGiDMesh(mesh_file1);
          PrintGiDResults(res_file);
//           mOctree.PrintGiDMeshNew(mesh_file1);

          Timer::Stop("Writing Gid conform Mesh");

          KRATOS_WATCH(mrBodyModelPart);

          KRATOS_CATCH("");
      }

      void GenerateOctree()
      {
          Timer::Start("Generating Octree");

          for(ModelPart::NodeIterator i_node = mrSkinModelPart.NodesBegin() ; i_node != mrSkinModelPart.NodesEnd() ; i_node++)
          {
              double temp_point[3];
              
              temp_point[0] = i_node->Coordinate(1);
              temp_point[1] = i_node->Coordinate(2);
              temp_point[2] = i_node->Coordinate(3);
              
              mOctree.Insert(temp_point);
          }
          
          mOctree.Constrain2To1();

          for(ModelPart::ElementIterator i_element = mrSkinModelPart.ElementsBegin() ; i_element != mrSkinModelPart.ElementsEnd() ; i_element++)
          {
              mOctree.Insert(*(i_element).base());
          }

          Timer::Stop("Generating Octree");
//          octree.Insert(*(mrSkinModelPart.ElementsBegin().base()));
          KRATOS_WATCH(mOctree);



      }

      void GenerateCellNodalData()
      {
          Timer::Start("Generating Cell Nodal Data");
          std::vector<OctreeType::cell_type*> all_leaves;
          mOctree.GetAllLeavesVector(all_leaves);

#pragma omp parallel for
          for (std::size_t i = 0; i < all_leaves.size(); i++)
          {
              *(all_leaves[i]->pGetDataPointer()) = ConfigurationType::AllocateData();
          }

          std::size_t last_id = mrBodyModelPart.NumberOfNodes() + 1;
          for (std::size_t i = 0; i < all_leaves.size(); i++)
          {
                CellType* cell = all_leaves[i];
                GenerateCellNode(cell,last_id);
          }

          Timer::Stop("Generating Cell Nodal Data");

      }

      void GenerateCellNode(CellType* pCell, std::size_t& LastId)
      {
        for (int i_pos=0; i_pos < 8; i_pos++) // position 8 is for center
        {
            ConfigurationType::CellNodeData* p_node = (*(pCell->pGetData()))[i_pos];
            if(p_node == 0)
            {
                CellType::key_type keys[3];
                pCell->GetKey(i_pos,keys);
                
                double new_point[3];
                
                (*(pCell->pGetData()))[i_pos] = new ConfigurationType::CellNodeData;
                (*(pCell->pGetData()))[i_pos]->Id() = LastId++;
                
                (*(pCell->pGetData()))[i_pos]->X() = pCell->GetCoordinate(keys[0]);
                (*(pCell->pGetData()))[i_pos]->Y() = pCell->GetCoordinate(keys[1]);
                (*(pCell->pGetData()))[i_pos]->Z() = pCell->GetCoordinate(keys[2]);
                
                mOctreeNodes.push_back((*(pCell->pGetData()))[i_pos]);
                
                SetNodeInNeighbours(pCell,i_pos,(*(pCell->pGetData()))[i_pos]);
            }

        }
      }
      
//       void GenerateCellNode(CellType* pCell, std::size_t& LastId)
//       {
//         for (int i_pos=0; i_pos < 8; i_pos++) // position 8 is for center
//         {
//             Node<3>* p_node = (*(pCell->pGetData()))[i_pos];
//             if(p_node == 0)
//             {
//                 CellType::key_type keys[3];
//                 pCell->GetKey(i_pos,keys);
// 
//                 double new_point[3];
// 
//                 new_point[0] = pCell->GetCoordinate(keys[0]);
//                 new_point[1] = pCell->GetCoordinate(keys[1]);
//                 new_point[2] = pCell->GetCoordinate(keys[2]);
//                 
// 
//                 (*(pCell->pGetData()))[i_pos] = (mrBodyModelPart.CreateNewNode(++LastId, new_point[0], new_point[1], new_point[2])).get();
// 
//                 SetNodeInNeighbours(pCell,i_pos,(*(pCell->pGetData()))[i_pos]);
//             }
// 
//         }
//       }

      void SetNodeInNeighbours(CellType* pCell, int Position, ConfigurationType::CellNodeData* pNode)
      {
            CellType::key_type point_key[3];
            pCell->GetKey(Position, point_key);
            
            for (std::size_t i_direction = 0; i_direction < 8; i_direction++) {
                CellType::key_type neighbour_key[3];
                if (pCell->GetNeighbourKey(Position, i_direction, neighbour_key)) {
                    CellType* neighbour_cell = mOctree.pGetCell(neighbour_key);
                    if (!neighbour_cell || (neighbour_cell == pCell))
                        continue;

                    std::size_t position = neighbour_cell->GetLocalPosition(point_key);
                    if((*neighbour_cell->pGetData())[position])
                    {
                        //std::cout << "ERROR!! Bad Position calculated!!!!!!!!!!! position :" << position << std::endl;
                        continue;
                    }
                    
                    (*neighbour_cell->pGetData())[position] = pNode;
                }
            }
      }

      void CalculateDistance()
      {
          Timer::Start("Calculate Distances");
          ModelPart::NodesContainerType::ContainerType& nodes = mrBodyModelPart.NodesArray();
          int nodes_size = nodes.size();
          // first of all we reset the node distance to 1.00 which is the maximum distnace in our normalized space.
#pragma omp parallel for firstprivate(nodes_size)
          for(int i = 0 ; i < nodes_size ; i++)
              nodes[i]->GetSolutionStepValue(DISTANCE) = 1.00;

            std::vector<CellType*> leaves;

          mOctree.GetAllLeavesVector(leaves);
          int leaves_size = leaves.size();

          for(int i = 0 ; i < leaves_size ; i++)
              CalculateNotEmptyLeavesDistance(leaves[i]);

#pragma omp parallel for firstprivate(nodes_size)
          for(int i = 0 ; i < nodes_size ; i++)
          {
              CalculateNodeDistance(*(nodes[i]));
          }
          Timer::Stop("Calculate Distances");

      }

      void CalculateDistance2()
      {
          Timer::Start("Calculate Distances 2");
          ModelPart::NodesContainerType::ContainerType& nodes = mrBodyModelPart.NodesArray();
          int nodes_size = nodes.size();
          // first of all we reste the node distance to 1.00 which is the maximum distnace in our normalized space.
#pragma omp parallel for firstprivate(nodes_size)
          for(int i = 0 ; i < nodes_size ; i++)
              nodes[i]->GetSolutionStepValue(DISTANCE) = 1.00;


            std::vector<CellType*> leaves;

          mOctree.GetAllLeavesVector(leaves);
          int leaves_size = leaves.size();

          for(int i = 0 ; i < leaves_size ; i++)
              CalculateNotEmptyLeavesDistance(leaves[i]);

          for(int i_direction = 0 ; i_direction < 1 ; i_direction++)
          {

//#pragma omp parallel for firstprivate(nodes_size)
            for(int i = 0 ; i < nodes_size ; i++)
            {
                if(nodes[i]->X() < 1.00 && nodes[i]->Y() < 1.00 && nodes[i]->Z() < 1.00)
               // if((*nodes[i])[i_direction] == 0.00)
                    CalculateDistance(*(nodes[i]), i_direction);
            }
          }
          Timer::Stop("Calculate Distances 2");

      }

      void CalculateDistance(Node<3>& rNode, int i_direction)
      {
//           double coords[3] = {rNode.X(), rNode.Y(), rNode.Z()};
//          // KRATOS_WATCH_3(coords);
// 
//              //This function must color the positions in space defined by 'coords'.
//             //coords is of dimension (3) normalized in (0,1)^3 space
// 
//             typedef Element::GeometryType triangle_type;
//             typedef std::vector<std::pair<double, triangle_type*> > intersections_container_type;
// 
//             intersections_container_type intersections;
//             std::vector<Node<3>*> nodes_array;
// 
// 
//             const double epsilon = 1e-12;
// 
//             double distance = 1.0;
// 
//             // Creating the ray
//             double ray[3] = {coords[0], coords[1], coords[2]};
//             ray[i_direction] = 0; // starting from the lower extreme
// 
// //            KRATOS_WATCH_3(ray)
//             GetIntersectionsAndNodes(ray, i_direction, intersections, nodes_array);
// //            KRATOS_WATCH(nodes_array.size())
//             for (int i_node = 0; i_node < nodes_array.size() ; i_node++)
//             {
//                 double coord = nodes_array[i_node]->Coordinate(i_direction+1);
//    //             KRATOS_WATCH(intersections.size());
// 
//                 int ray_color= 1;
//                 std::vector<std::pair<double, Element::GeometryType*> >::iterator i_intersection = intersections.begin();
//                 while (i_intersection != intersections.end()) {
//                     double d = coord - i_intersection->first;
//                     if (d > epsilon) {
// 
//                         ray_color = -ray_color;
//                         distance = d;
//                     } else if (d > -epsilon) {//interface
//                         distance = 0.00;
//                         break;
//                     } else {
//                         if(distance > -d)
//                             distance = -d;
//                         break;
//                     }
// 
//                     i_intersection++;
//                 }
// 
//                 distance *= ray_color;
// 
//                 double& node_distance = nodes_array[i_node]->GetSolutionStepValue(DISTANCE);
//                 if(fabs(distance) < fabs(node_distance))
//                     node_distance = distance;
//                 else if (distance*node_distance < 0.00) // assigning the correct sign
//                     node_distance = -node_distance;
// 
// 
//             }
     }

      void CalculateNotEmptyLeavesDistance(CellType* pCell)
      {
        typedef Element::GeometryType triangle_type;
        typedef OctreeType::cell_type::object_container_type object_container_type;

        object_container_type* objects = (pCell->pGetObjects());

        // There are no intersection in empty cells
        if (objects->empty())
            return;


        for (int i_pos=0; i_pos < 8; i_pos++) // position 8 is for center
        {
            double distance = 1.00; // maximum distance is 1.00

            for(object_container_type::iterator i_object = objects->begin(); i_object != objects->end(); i_object++)
            {
                CellType::key_type keys[3];
                pCell->GetKey(i_pos,keys);

                double cell_point[3];

                //cell_point[0] = pCell->GetCoordinate(keys[0]);
                //cell_point[1] = pCell->GetCoordinate(keys[1]);
                //cell_point[2] = pCell->GetCoordinate(keys[2]);

                double d = GeometryUtils::PointDistanceToTriangle3D((*i_object)->GetGeometry()[0], (*i_object)->GetGeometry()[1], (*i_object)->GetGeometry()[2], Point(cell_point[0], cell_point[1], cell_point[2]));

                if(d < distance)
                    distance = d;
            }

            double& node_distance = (*(pCell->pGetData()))[i_pos]->Distance();
            if(distance < node_distance)
                node_distance = distance;

        }

      }


      void CalculateNodeDistance(Node<3>& rNode)
      {
          double coord[3] = {rNode.X(), rNode.Y(), rNode.Z()};
          double distance = DistancePositionInSpace(coord);
          double& node_distance = rNode.GetSolutionStepValue(DISTANCE);

          const double epsilon = 1.00e-12;
          if(fabs(node_distance) > fabs(distance))
            node_distance = distance;
          else if (distance*node_distance < 0.00) // assigning the correct sign
              node_distance = -node_distance;
      }

      double DistancePositionInSpace(double* coords)
      {
            //This function must color the positions in space defined by 'coords'.
            //coords is of dimension (3) normalized in (0,1)^3 space

            typedef Element::GeometryType triangle_type;
            typedef std::vector<std::pair<double, triangle_type*> > intersections_container_type;

            intersections_container_type intersections;

            const int dimension = 3;
            const double epsilon = 1e-12;

            double distances[3] = {1.0, 1.0, 1.0};

            for (int i_direction = 0; i_direction < dimension; i_direction++)
            {
                // Creating the ray
                double ray[3] = {coords[0], coords[1], coords[2]};
                ray[i_direction] = 0; // starting from the lower extreme

                GetIntersections(ray, i_direction, intersections);

//                if(intersections.size() == 1)
//                    KRATOS_WATCH_3(ray)

   //             KRATOS_WATCH(intersections.size());

                int ray_color= 1;
                std::vector<std::pair<double, Element::GeometryType*> >::iterator i_intersection = intersections.begin();
                while (i_intersection != intersections.end()) {
                    double d = coords[i_direction] - i_intersection->first;
                    if (d > epsilon) {

                        ray_color = -ray_color;
                        distances[i_direction] = d;
//                        if(distances[i_direction] > d) // I think this is redundunt. Pooyan.
//                        {
//                            if(ray_color > 0.00)
//                                    distances[i_direction] = d;
//                            else
//                                distances[i_direction] = -d;
//                        }
                    } else if (d > -epsilon) {//interface
                        distances[i_direction] = 0.00;
                        break;
                    } else {
                        if(distances[i_direction] > -d)
                            distances[i_direction] = -d;
                        break;
                    }

                    i_intersection++;
                }

                distances[i_direction] *= ray_color;
            }

//            if(distances[0]*distances[1] < 0.00 || distances[2]*distances[1] < 0.00)
//                KRATOS_WATCH_3(distances);

#ifdef _DEBUG
            std::cout << "colors : " << colors[0] << ", " << colors[1] << ", " << colors[2] << std::endl;
#endif
            double distance = (fabs(distances[0]) > fabs(distances[1])) ? distances[1] : distances[0];
            distance = (fabs(distance) > fabs(distances[2])) ? distances[2] : distance;

            return distance;

        }


      void GetIntersectionsAndNodes(double* ray, int direction, std::vector<std::pair<double,Element::GeometryType*> >& intersections, std::vector<Node<3>*>& rNodesArray)
      {
    //     //This function passes the ray through the model and gives the hit point to all objects in its way
    //     //ray is of dimension (3) normalized in (0,1)^3 space
    //     // direction can be 0,1,2 which are x,y and z respectively
    // 
    //     const double epsilon = 1.00e-12;
    // 
    //     // first clearing the intersections points vector
    //     intersections.clear();
    // 
    //     OctreeType* octree = &mOctree;
    // 
    //     OctreeType::key_type ray_key[3] = {octree->Key(ray[0]), octree->Key(ray[1]), octree->Key(ray[2])}; //ASK_TOKEN
    //     OctreeType::key_type cell_key[3];
    // 
    //     // getting the entrance cell from lower extreme
    //     ray_key[direction] = 0;
    //     OctreeType::cell_type* cell = octree->pGetCell(ray_key);
    // 
    //     while (cell) {
    //         std::size_t position = cell->GetLocalPosition(ray_key); // Is this the local position!?!?!?!
    //         OctreeType::key_type node_key[3];
    //         cell->GetKey(position, node_key);
    //         if((node_key[0] == ray_key[0]) && (node_key[1] == ray_key[1]) && (node_key[2] == ray_key[2]))
    //         {
    //             if(cell->pGetData())
    //             {
    //                 if(cell->pGetData()->size() > position)
    //                 {
    //                     Node<3>* p_node = (*cell->pGetData())[position];
    //                     if(p_node)
    //                     {
    //                         //KRATOS_WATCH(p_node->Id())
    //                         rNodesArray.push_back(p_node);
    //                     }
    //                 }
    //                 else
    //                     KRATOS_WATCH(cell->pGetData()->size())
    //             }
    //         }
    // 
    // 
    // //        std::cout << ".";
    //       GetCellIntersections(cell, ray, ray_key, direction, intersections);
    // 
    //       // Add the cell's middle node if existed
    // //      cell->GetKey(8, cell_key); // 8 is the central position
    // //      ray_key[direction]=cell_key[direction]; // positioning the ray in the middle of cell in its direction
    // 
    // //      position = cell->GetLocalPosition(ray_key);
    // //      if(position < 27) // principal nodes
    // //      {
    // //          if(cell->pGetData())
    // //          {
    // //              if(cell->pGetData()->size() > position)
    // //              {
    // //                  Node<3>* p_node = (*cell->pGetData())[position];
    // //                  if(p_node)
    // //                  {
    // //                      //KRATOS_WATCH(p_node->Id())
    // //                      rNodesArray.push_back(p_node);
    // //                  }
    // //              }
    // //              else
    // //                  KRATOS_WATCH(cell->pGetData()->size())
    // //          }
    // //      }
    // //      else
    // //      {
    // //          KRATOS_WATCH(position);
    // //          KRATOS_WATCH(*cell);
    // //      }
    // 
    // 
    //       // go to the next cell
    //       if (cell->GetNeighbourKey(1 + direction * 2, cell_key)) {
    //         ray_key[direction] = cell_key[direction];
    //         cell = octree->pGetCell(ray_key);
    //         ray_key[direction] -= 1 ;//the key returned by GetNeighbourKey is inside the cell (minkey +1), to ensure that the corresponding
    //         //cell get in pGetCell is the right one.
    // #ifdef _DEBUG
    //         Octree_Pooyan::key_type min_key[3];
    //         cell->GetMinKey(min_key[0],min_key[1],min_key[2]);
    //         Octree_Pooyan::key_type tmp;
    //         tmp= min_key[direction];
    //         assert(ray_key[direction]==tmp);
    // #endif
    //       } else
    //         cell = NULL;
    //     }
    // 
    // 
    // 
    //  //   KRATOS_WATCH(rNodesArray.size());
    //     // now eliminating the repeated objects
    //     if (!intersections.empty()) {
    //       //sort
    //       std::sort(intersections.begin(), intersections.end());
    //       // unique
    //       std::vector<std::pair<double, Element::GeometryType*> >::iterator i_begin = intersections.begin();
    //       std::vector<std::pair<double, Element::GeometryType*> >::iterator i_intersection = intersections.begin();
    //       while (++i_begin != intersections.end()) {
    //           // considering the very near points as the same points
    //           if (fabs(i_begin->first - i_intersection->first) > epsilon) // if the hit points are far enough they are not the same
    //             *(++i_intersection) = *i_begin;
    //       }
    //       intersections.resize((++i_intersection) - intersections.begin());
    // 
    //     }
      }

      void GetIntersections(double* ray, int direction, std::vector<std::pair<double,Element::GeometryType*> >& intersections)
      {
    //     //This function passes the ray through the model and gives the hit point to all objects in its way
    //     //ray is of dimension (3) normalized in (0,1)^3 space
    //     // direction can be 0,1,2 which are x,y and z respectively
    // 
    //     const double epsilon = 1.00e-12;
    // 
    //     // first clearing the intersections points vector
    //     intersections.clear();
    // 
    //     OctreeType* octree = &mOctree;
    // 
    //     OctreeType::key_type ray_key[3] = {octree->Key(ray[0]), octree->Key(ray[1]), octree->Key(ray[2])};
    //     OctreeType::key_type cell_key[3];
    // 
    //     // getting the entrance cell from lower extreme
    //     OctreeType::cell_type* cell = octree->pGetCell(ray_key);
    // 
    //     while (cell) {
    // //        std::cout << ".";
    //       GetCellIntersections(cell, ray, ray_key, direction, intersections);
    //       // go to the next cell
    //       if (cell->GetNeighbourKey(1 + direction * 2, cell_key)) {
    //         ray_key[direction] = cell_key[direction];
    //         cell = octree->pGetCell(ray_key);
    //         ray_key[direction] -= 1 ;//the key returned by GetNeighbourKey is inside the cell (minkey +1), to ensure that the corresponding
    //         //cell get in pGetCell is the right one.
    // #ifdef _DEBUG
    //         Octree_Pooyan::key_type min_key[3];
    //         cell->GetMinKey(min_key[0],min_key[1],min_key[2]);
    //         Octree_Pooyan::key_type tmp;
    //         tmp= min_key[direction];
    //         assert(ray_key[direction]==tmp);
    // #endif
    //       } else
    //         cell = NULL;
    //     }
    // 
    // 
    //     // now eliminating the repeated objects
    //     if (!intersections.empty()) {
    //       //sort
    //       std::sort(intersections.begin(), intersections.end());
    //       // unique
    //       std::vector<std::pair<double, Element::GeometryType*> >::iterator i_begin = intersections.begin();
    //       std::vector<std::pair<double, Element::GeometryType*> >::iterator i_intersection = intersections.begin();
    //       while (++i_begin != intersections.end()) {
    //           // considering the very near points as the same points
    //           if (fabs(i_begin->first - i_intersection->first) > epsilon) // if the hit points are far enough they are not the same
    //             *(++i_intersection) = *i_begin;
    //       }
    //       intersections.resize((++i_intersection) - intersections.begin());
    // 
    //     }
      }

      int GetCellIntersections(OctreeType::cell_type* cell, double* ray,
        OctreeType::key_type* ray_key, int direction,
        std::vector<std::pair<double, Element::GeometryType*> >& intersections)  {
    //       //This function passes the ray through the cell and gives the hit point to all objects in its way
    //       //ray is of dimension (3) normalized in (0,1)^3 space
    //       // direction can be 0,1,2 which are x,y and z respectively
    // 
    //       typedef Element::GeometryType triangle_type;
    //       typedef OctreeType::cell_type::object_container_type object_container_type;
    // 
    //       object_container_type* objects = (cell->pGetObjects());
    // 
    //       // There are no intersection in empty cells
    //       if (objects->empty())
    //         return 0;
    // 
    // //      std::cout << "X";
    //       // calculating the two extreme of the ray segment inside the cell
    //       double ray_point1[3] = {ray[0], ray[1], ray[2]};
    //       double ray_point2[3] = {ray[0], ray[1], ray[2]};
    //       ray_point1[direction] = cell->GetCoordinate(ray_key[direction]);
    //       ray_point2[direction] = ray_point1[direction] + cell->GetSize();
    // 
    //       for (object_container_type::iterator i_object = objects->begin(); i_object != objects->end(); i_object++) {
    //         double intersection[3]={0.00,0.00,0.00};
    // 
    //         int is_intersected = IntersectionTriangleSegment((*i_object)->GetGeometry(), ray_point1, ray_point2, intersection); // This intersection has to be optimized for axis aligned rays
    // 
    //         if (is_intersected == 1) // There is an intersection but not coplanar
    //           intersections.push_back(std::pair<double, Element::GeometryType*>(intersection[direction], &((*i_object)->GetGeometry())));
    //         //else if(is_intersected == 2) // coplanar case
    //       }
    // 
    //       return 0;
      }

      int IntersectionTriangleSegment(Element::GeometryType& rGeometry, double* RayPoint1, double* RayPoint2, double* IntersectionPoint)
      {
          // This is the adaption of the implemnetation provided in:
          // http://www.softsurfer.com/Archive/algorithm_0105/algorithm_0105.htm#intersect_RayTriangle()

          const double epsilon = 1.00e-12;

        array_1d<double,3>    u, v, n;             // triangle vectors
        array_1d<double,3>    dir, w0, w;          // ray vectors
        double     r, a, b;             // params to calc ray-plane intersect


        // get triangle edge vectors and plane normal
        u = rGeometry[1] - rGeometry[0];
        v = rGeometry[2] - rGeometry[0];
        
        MathUtils<double>::CrossProduct(n, u, v);             // cross product
        
        if (norm_2(n) == 0)            // triangle is degenerate
            return -1;                 // do not deal with this case

        for(int i = 0 ; i < 3 ; i++)
        {
            dir[i] = RayPoint2[i] - RayPoint1[i];             // ray direction vector
            w0[i] = RayPoint1[i] - rGeometry[0][i];
        }

        a = -inner_prod(n,w0);
        b = inner_prod(n,dir);

        if (fabs(b) < epsilon) {     // ray is parallel to triangle plane
            if (a == 0)                // ray lies in triangle plane
                return 2;
            else return 0;             // ray disjoint from plane
        }

        // get intersect point of ray with triangle plane
        r = a / b;
        if (r < 0.0)                   // ray goes away from triangle
            return 0;                  // => no intersect
        // for a segment, also test if (r > 1.0) => no intersect

        for(int i = 0 ; i < 3 ; i++)
          IntersectionPoint[i]  = RayPoint1[i] + r * dir[i];           // intersect point of ray and plane

        // is I inside T?
        double    uu, uv, vv, wu, wv, D;
        uu = inner_prod(u,u);
        uv = inner_prod(u,v);
        vv = inner_prod(v,v);


        for(int i = 0 ; i < 3 ; i++)
            w[i] = IntersectionPoint[i] - rGeometry[0][i];


        wu = inner_prod(w,u);
        wv = inner_prod(w,v);
        D = uv * uv - uu * vv;

        // get and test parametric coords
        double s, t;
        s = (uv * wv - vv * wu) / D;
        if (s < 0.0 - epsilon || s > 1.0 + epsilon)        // I is outside T
            return 0;
        t = (uv * wu - uu * wv) / D;
        if (t < 0.0 - epsilon || (s + t) > 1.0 + epsilon)  // I is outside T
            return 0;

        return 1;                      // I is in T

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
	  return "CalculateSignedDistanceTo3DSkinProcess";
	}
      
      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const
	{
	  rOStream << "CalculateSignedDistanceTo3DSkinProcess";
	}

      /// Print object's data.
      virtual void PrintData(std::ostream& rOStream) const
	{
	}
      
        void PrintGiDMesh(std::ostream & rOStream) const {
            std::vector<CellType*> leaves;

            mOctree.GetAllLeavesVector(leaves);

            std::cout << "writing " << leaves.size() << " leaves" << std::endl;
            rOStream << "MESH \"leaves\" dimension 3 ElemType Hexahedra Nnode 8" << std::endl;
            rOStream << "# color 96 96 96" << std::endl;
            rOStream << "Coordinates" << std::endl;
            rOStream << "# node_number coordinate_x coordinate_y coordinate_z  " << std::endl;

            for(DistanceSpatialContainersConfigure::data_type::const_iterator i_node = mOctreeNodes.begin() ; i_node != mOctreeNodes.end() ; i_node++)
            {
                rOStream << (*i_node)->Id() << "  " << (*i_node)->Coordinate(1) << "  " << (*i_node)->Coordinate(2) << "  " << (*i_node)->Coordinate(3) << std::endl;
            }
            
            std::cout << "Nodes written..." << std::endl;
            rOStream << "end coordinates" << std::endl;
            rOStream << "Elements" << std::endl;
            rOStream << "# element n1 n2 n3 n4 n5 n6 n7 n8" << std::endl;

            for (std::size_t i = 0; i < leaves.size(); i++) 
            {
                if ((leaves[i]->pGetData()))
                { 
                    DistanceSpatialContainersConfigure::data_type& nodes = (*(leaves[i]->pGetData()));

//                     std::cout << "Leave - Level: "  << nodes[0]->Id() << " " << nodes[1]->Id() << " " << nodes[2]->Id() << " etc... " << std::endl;
                    
                    rOStream << i + 1;
                    for(int j = 0 ; j < 8 ; j++)
                        rOStream << "  " << nodes[j]->Id();
                    rOStream << std::endl;
                }
            }
            
            rOStream << "end elements" << std::endl;

        }

        void PrintGiDResults(std::ostream & rOStream) const {
            std::vector<CellType*> leaves;

            mOctree.GetAllLeavesVector(leaves);

            rOStream << "GiD Post Results File 1.0" << std::endl << std::endl;

            rOStream << "Result \"Distance\" \"Kratos\" 1 Scalar OnNodes" << std::endl;

            rOStream << "Values" << std::endl;

          for(ModelPart::NodeIterator i_node = mrBodyModelPart.NodesBegin() ; i_node != mrBodyModelPart.NodesEnd() ; i_node++)
          {
              rOStream << i_node->Id() << "  " << i_node->GetSolutionStepValue(DISTANCE) << std::endl;
          }
            rOStream << "End Values" << std::endl;

        }
            
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
        
      ModelPart& mrSkinModelPart;
      ModelPart& mrBodyModelPart;
      
      DistanceSpatialContainersConfigure::data_type mOctreeNodes;

      OctreeType mOctree;  
        
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
      CalculateSignedDistanceTo3DSkinProcess& operator=(CalculateSignedDistanceTo3DSkinProcess const& rOther);

      /// Copy constructor.
      //CalculateSignedDistanceTo3DSkinProcess(CalculateSignedDistanceTo3DSkinProcess const& rOther);

        
      ///@}    
        
    }; // Class CalculateSignedDistanceTo3DSkinProcess

  ///@} 
  
  ///@name Type Definitions       
  ///@{ 
  
  
  ///@} 
  ///@name Input and output 
  ///@{ 
        
 
  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream, 
				    CalculateSignedDistanceTo3DSkinProcess& rThis);

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream, 
				    const CalculateSignedDistanceTo3DSkinProcess& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }
  ///@} 
  
  
}  // namespace Kratos.

#endif // KRATOS_CALCULATE_DISTANCE_PROCESS_H_INCLUDED  defined 


