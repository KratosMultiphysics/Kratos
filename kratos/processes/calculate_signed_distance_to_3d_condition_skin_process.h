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
//                   Daniel Baumgaertner
//                   Johannes Wolf
//



#if !defined(KRATOS_CALCULATE_DISTANCE_CONDITION_PROCESS_H_INCLUDED )
#define  KRATOS_CALCULATE_DISTANCE_CONDITION_PROCESS_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "includes/model_part.h"
#include "includes/deprecated_variables.h"
#include "spatial_containers/octree_binary.h"
#include "utilities/spatial_containers_configure.h"
#include "utilities/timer.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h"
#include "geometries/triangle_3d_3.h"
#include "utilities/body_normal_calculation_utils.h"


namespace Kratos
{

class DistanceSpatialContainersConditionConfigure
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
        double& operator[](int i) {return mCoordinates[i];}
        std::size_t& Id(){return mId;}
    };




    ///@name Type Definitions
    ///@{

    enum { Dimension = 3,
           DIMENSION = 3,
           MAX_LEVEL = 12,
           MIN_LEVEL = 2    // this cannot be less than 2!!!
         };

    typedef Point                                               PointType;  /// always the point 3D
    typedef std::vector<double>::iterator                       DistanceIteratorType;
    typedef PointerVectorSet<
                GeometricalObject::Pointer, 
                IndexedObject,
                std::less<typename IndexedObject::result_type>,
                std::equal_to<typename IndexedObject::result_type>,
                Kratos::shared_ptr<typename GeometricalObject::Pointer>,
                std::vector< Kratos::shared_ptr<typename GeometricalObject::Pointer> >
                    >  ContainerType;
    typedef ContainerType::value_type                           PointerType;
    typedef ContainerType::iterator                             IteratorType;
    typedef PointerVectorSet<
                GeometricalObject::Pointer, 
                IndexedObject,
                std::less<typename IndexedObject::result_type>,
                std::equal_to<typename IndexedObject::result_type>,
                Kratos::shared_ptr<typename GeometricalObject::Pointer>,
                std::vector< Kratos::shared_ptr<typename GeometricalObject::Pointer> >
                    >  ResultContainerType;
    typedef ResultContainerType::value_type                     ResultPointerType;
    typedef ResultContainerType::iterator                       ResultIteratorType;

    typedef GeometricalObject::Pointer                          pointer_type;
    typedef CellNodeData                                        cell_node_data_type;
    typedef std::vector<CellNodeData*> data_type;

    typedef std::vector<PointerType>::iterator                  PointerTypeIterator;




    /// Pointer definition of DistanceSpatialContainersConditionConfigure
    KRATOS_CLASS_POINTER_DEFINITION(DistanceSpatialContainersConditionConfigure);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    DistanceSpatialContainersConditionConfigure() {}

    /// Destructor.
    virtual ~DistanceSpatialContainersConditionConfigure() {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    static data_type* AllocateData() {
        return new data_type(27, (CellNodeData*)NULL);
    }

    static void CopyData(data_type* source, data_type* destination) {
        *destination = *source;
    }

    static void DeleteData(data_type* data) {
        delete data;
    }
///******************************************************************************************************************
///******************************************************************************************************************

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

///******************************************************************************************************************
///******************************************************************************************************************

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

///******************************************************************************************************************
///******************************************************************************************************************

    static inline bool Intersection(const PointerType& rObj_1, const PointerType& rObj_2)
    {
        Element::GeometryType& geom_1 = rObj_1->GetGeometry();
        Element::GeometryType& geom_2 = rObj_2->GetGeometry();
        return  geom_1.HasIntersection(geom_2);

    }


///******************************************************************************************************************
///******************************************************************************************************************

    static inline bool  IntersectionBox(const PointerType& rObject,  const PointType& rLowPoint, const PointType& rHighPoint)
    {
        return rObject->GetGeometry().HasIntersection(rLowPoint, rHighPoint);
    }

        ///******************************************************************************************************************
        ///******************************************************************************************************************


          static  inline bool  IsIntersected(const Element::Pointer rObject, double Tolerance, const double* rLowPoint, const double* rHighPoint)
         {
             Point low_point(rLowPoint[0] - Tolerance, rLowPoint[1] - Tolerance, rLowPoint[2] - Tolerance);
             Point high_point(rHighPoint[0] + Tolerance, rHighPoint[1] + Tolerance, rHighPoint[2] + Tolerance);

             KRATOS_THROW_ERROR(std::logic_error, "Not Implemented method", "")
             //return HasIntersection(rObject->GetGeometry(), low_point, high_point);
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
    DistanceSpatialContainersConditionConfigure& operator=(DistanceSpatialContainersConditionConfigure const& rOther);

    /// Copy constructor.
    DistanceSpatialContainersConditionConfigure(DistanceSpatialContainersConditionConfigure const& rOther);


}; // Class DistanceSpatialContainersConditionConfigure

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
  class CalculateSignedDistanceTo3DConditionSkinProcess
	: public Process
    {
    public:
      ///@name Type Definitions
      ///@{

      /// Pointer definition of CalculateSignedDistanceTo3DConditionSkinProcess
        KRATOS_CLASS_POINTER_DEFINITION(CalculateSignedDistanceTo3DConditionSkinProcess);

        typedef DistanceSpatialContainersConditionConfigure ConfigurationType;
        typedef OctreeBinaryCell<ConfigurationType> CellType;
        typedef OctreeBinary<CellType> OctreeType;
        typedef ConfigurationType::cell_node_data_type CellNodeDataType;
        typedef Point PointType;  /// always the point 3D
        typedef OctreeType::cell_type::object_container_type object_container_type;
        typedef struct{
            array_1d<double,3>  Coordinates;
            array_1d<double,3>  StructElemNormal;
        }IntersectionNodeStruct;
        typedef struct{
            std::vector<IntersectionNodeStruct> IntNodes;
        }TetEdgeStruct;


      ///@}
      ///@name Life Cycle
      ///@{

      /// Constructor.
      CalculateSignedDistanceTo3DConditionSkinProcess(ModelPart& rThisModelPartStruc, ModelPart& rThisModelPartFluid)
          : mrSkinModelPart(rThisModelPartStruc), mrBodyModelPart(rThisModelPartStruc), mrFluidModelPart(rThisModelPartFluid)
	{
	}

      /// Destructor.
      ~CalculateSignedDistanceTo3DConditionSkinProcess() override
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

      ///******************************************************************************************************************
      ///******************************************************************************************************************

      void Execute() override
      {
          KRATOS_TRY
	/*
            std::cout << "Clearing the list of correspondances between the FIXED MESH ELEMENTS and the EMBEDDED CONDITIONS THAT ARE CROSSING THEM..." << std::endl;


	  for( ModelPart::ElementIterator i_fluidElement = mrFluidModelPart.ElementsBegin();
                                          i_fluidElement != mrFluidModelPart.ElementsEnd();
                                          i_fluidElement++)
          {
              //i_fluidElement->GetValue(NEIGHBOUR_CONDITIONS).resize(0);
		( i_fluidElement->GetValue(NEIGHBOUR_EMBEDDED_FACES)).reserve(6);

		 WeakPointerVector<GeometricalObject >& rE = i_fluidElement->GetValue(NEIGHBOUR_EMBEDDED_FACES);
	            rE.erase(rE.begin(),rE.end() );
          }
	*/


            //std::cout << "Generating the Octree..." << std::endl;
            GenerateOctree();
	  //std::cout << "Generating the Octree finished" << std::endl;

          DistanceFluidStructure();

            //GenerateNodes();
            CalculateDistance2(); // I have to change this. Pooyan.






//            CalculateDistance();
//            CalculateDistance2();

//            double coord[3] = {0.4375, 0.57812, 0.5};
//            double distance = DistancePositionInSpace(coord);
//            KRATOS_WATCH(distance);



            //mrSkinModelPart.GetCommunicator().AssembleCurrentData(DISTANCE);

//            std::ofstream mesh_file1("octree1.post.msh");
//            std::ofstream res_file("octree1.post.res");

//            Timer::Start("Writing Gid conform Mesh");

//            PrintGiDMesh(mesh_file1);
//            PrintGiDResults(res_file);
            //octree.PrintGiDMeshNew(mesh_file2);

//            Timer::Stop("Writing Gid conform Mesh");

//            KRATOS_WATCH(mrBodyModelPart);

            //delete octree. TODO: Carlos

       	    KRATOS_CATCH("");
	}

      ///******************************************************************************************************************
      ///******************************************************************************************************************

      void DistanceFluidStructure()
      {

          // Initialize nodal distances each node in the domain to 1.0
          InitializeDistances();

          // Initialize index table to define line Edges of fluid element
          BoundedMatrix<unsigned int,6,2> TetEdgeIndexTable;
          SetIndexTable(TetEdgeIndexTable);

	for( ModelPart::ElementIterator i_fluidElement = mrFluidModelPart.ElementsBegin();
                                          i_fluidElement != mrFluidModelPart.ElementsEnd();
                                          i_fluidElement++)
          {
		i_fluidElement->GetValue(EMBEDDED_VELOCITY)=ZeroVector(3);
	  }

          // loop over all fluid elements
          for( ModelPart::ElementIterator i_fluidElement = mrFluidModelPart.ElementsBegin();
                                          i_fluidElement != mrFluidModelPart.ElementsEnd();
                                          i_fluidElement++)
          {
              CalcNodalDistancesOfTetNodes( i_fluidElement , TetEdgeIndexTable );
          }
        std::cout<<" finishded distance fluid structure"<<std::endl;

		KRATOS_WATCH("ENDOF LOOP")
      }

      ///******************************************************************************************************************
      ///******************************************************************************************************************

      void InitializeDistances()
      {
          ModelPart::NodesContainerType::ContainerType& nodes = mrFluidModelPart.NodesArray();

          // reset the node distance to 1.0 which is the maximum distance in our normalized space.
          int nodesSize = nodes.size();
          for(int i = 0 ; i < nodesSize ; i++)
                  nodes[i]->GetSolutionStepValue(DISTANCE) = 1.0;

          ModelPart::ElementsContainerType::ContainerType& fluid_elements = mrFluidModelPart.ElementsArray();

          array_1d<double,4> ElementalDistances;
          const double initial_distance = 1.0;
          ElementalDistances[0] = initial_distance;
          ElementalDistances[1] = initial_distance;
          ElementalDistances[2] = initial_distance;
          ElementalDistances[3] = initial_distance;

          // reset the elemental distance to 1.0 which is the maximum distance in our normalized space.
          for(unsigned int i = 0 ; i < fluid_elements.size() ; i++)
          {
              fluid_elements[i]->SetValue(ELEMENTAL_DISTANCES,ElementalDistances);
              fluid_elements[i]->GetValue(SPLIT_ELEMENT) = false;
          }
      }

      ///******************************************************************************************************************
      ///******************************************************************************************************************

      void SetIndexTable( BoundedMatrix<unsigned int,6,2>& TetEdgeIndexTable )
      {
          // Initialize index table to define line Edges of fluid element
          TetEdgeIndexTable(0,0) = 0;
          TetEdgeIndexTable(0,1) = 1;
          TetEdgeIndexTable(1,0) = 0;
          TetEdgeIndexTable(1,1) = 2;
          TetEdgeIndexTable(2,0) = 0;
          TetEdgeIndexTable(2,1) = 3;
          TetEdgeIndexTable(3,0) = 1;
          TetEdgeIndexTable(3,1) = 2;
          TetEdgeIndexTable(4,0) = 1;
          TetEdgeIndexTable(4,1) = 3;
          TetEdgeIndexTable(5,0) = 2;
          TetEdgeIndexTable(5,1) = 3;
      }

      ///******************************************************************************************************************
      ///******************************************************************************************************************

      void CalcNodalDistancesOfTetNodes( ModelPart::ElementsContainerType::iterator& i_fluidElement,
                                         BoundedMatrix<unsigned int,6,2>            TetEdgeIndexTable)
      {
          std::vector<OctreeType::cell_type*> leaves;
          std::vector<TetEdgeStruct>          IntersectedTetEdges;
          unsigned int NumberIntersectionsOnTetCorner = 0;

          // Get leaves of octree intersecting with fluid element
          mOctree.GetIntersectedLeaves(*(i_fluidElement).base(),leaves);

	  int intersection_counter=0;
	  //i_fluidElement->GetValue(EMBEDDED_VELOCITY)=ZeroVector(3);
          // Loop over all 6 line Edges of the tetrahedra
          for(unsigned int i_tetEdge = 0; i_tetEdge < 6; i_tetEdge++)
          {
              IdentifyIntersectionNodes( i_fluidElement , i_tetEdge , leaves , IntersectedTetEdges ,
                                         NumberIntersectionsOnTetCorner , TetEdgeIndexTable, intersection_counter );
          }
	  if (intersection_counter!=0)
		  i_fluidElement->GetValue(EMBEDDED_VELOCITY)/=3.0*intersection_counter;
	//else
	//	i_fluidElement->GetValue(EMBEDDED_VELOCITY)=ZeroVector(3);
	//KRATOS_WATCH("============================================================")
	 // KRATOS_WATCH(i_fluidElement->GetValue(EMBEDDED_VELOCITY))
	//KRATOS_WATCH("???????????????????????????????????????????????????????????????")
	  //if (intersection_counter!=0)
	//	  KRATOS_WATCH(intersection_counter)


          if(IntersectedTetEdges.size() > 0)
              CalcNodalDistanceTo3DSkin( IntersectedTetEdges , i_fluidElement , NumberIntersectionsOnTetCorner );
      }

      ///******************************************************************************************************************
      ///******************************************************************************************************************

      void IdentifyIntersectionNodes( ModelPart::ElementsContainerType::iterator&   i_fluidElement,
                                      unsigned int                                  i_tetEdge,
                                      std::vector<OctreeType::cell_type*>&          leaves,
                                      std::vector<TetEdgeStruct>&                   IntersectedTetEdges,
                                      unsigned int&                                 NumberIntersectionsOnTetCorner,
                                      BoundedMatrix<unsigned int,6,2>              TetEdgeIndexTable,
				      int& intersection_counter)
      {


	  std::vector<unsigned int> IntersectingStructCondID;
          TetEdgeStruct             NewTetEdge;

          // Get nodes of line Edge
          unsigned int EdgeStartIndex = TetEdgeIndexTable(i_tetEdge,0);
          unsigned int EdgeEndIndex   = TetEdgeIndexTable(i_tetEdge,1);

          PointType& P1 = i_fluidElement->GetGeometry()[EdgeStartIndex];
          PointType& P2 = i_fluidElement->GetGeometry()[EdgeEndIndex];

          double EdgeNode1[3] = {P1.X() , P1.Y() , P1.Z()};
          double EdgeNode2[3] = {P2.X() , P2.Y() , P2.Z()};

	  //int count=0;
          // loop over all octree cells which are intersected by the fluid element
          for(unsigned int i_cell = 0 ; i_cell < leaves.size() ; i_cell++)
          {
              // Structural element contained in one cell of the octree
              object_container_type* struct_cond = (leaves[i_cell]->pGetObjects());

              // loop over all structural elements within each octree cell

              for(object_container_type::iterator i_StructCondition = struct_cond->begin(); i_StructCondition != struct_cond->end(); i_StructCondition++)
              {
			//KRATOS_WATCH(struct_cond->size())

                  if( StructuralElementNotYetConsidered( (*i_StructCondition)->Id() , IntersectingStructCondID ) )
                  {

                      // Calculate and associate intersection point to the current fluid element
                      double IntersectionPoint[3] = {0.0 , 0.0 , 0.0};
                      int TetEdgeHasIntersections = IntersectionTriangleSegment( (*i_StructCondition)->GetGeometry() , EdgeNode1 , EdgeNode2 , IntersectionPoint );

                      if( TetEdgeHasIntersections == 1 )
                      {
                          IntersectionNodeStruct NewIntersectionNode;

                          // Assign information to the intersection node
                          NewIntersectionNode.Coordinates[0] = IntersectionPoint[0];
                          NewIntersectionNode.Coordinates[1] = IntersectionPoint[1];
                          NewIntersectionNode.Coordinates[2] = IntersectionPoint[2];

                          if ( IsNewIntersectionNode( NewIntersectionNode , IntersectedTetEdges ) )
                          {

                              if( IsIntersectionNodeOnTetEdge( IntersectionPoint , EdgeNode1 , EdgeNode2 ) )
                              {

                                  // Calculate normal of the structural element at the position of the intersection point
                                  CalculateNormal3D((*i_StructCondition)->GetGeometry(),NewIntersectionNode.StructElemNormal);

                                  // add the new intersection point to the list of intersection points of the fluid element
                                  NewTetEdge.IntNodes.push_back(NewIntersectionNode);

				  //(i_fluidElement->GetValue(NEIGHBOUR_EMBEDDED_FACES)).push_back( GeometricalObject::WeakPointer( *(i_StructCondition.base()) ) );
				  /*
				  array_1d<double,3> emb_vel=(*i_StructCondition)->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);
				  emb_vel+=(*i_StructCondition)->GetGeometry()[1].FastGetSolutionStepValue(VELOCITY);
				  emb_vel+=(*i_StructCondition)->GetGeometry()[2].FastGetSolutionStepValue(VELOCITY);

				  //KRATOS_WATCH(emb_vel)
				  i_fluidElement->GetValue(EMBEDDED_VELOCITY)+=emb_vel;
				  intersection_counter++;
				  */

                                  // check, how many intersection nodes are located on corner points of the tetrahedra
                                  if ( IsIntersectionOnCorner( NewIntersectionNode , EdgeNode1 , EdgeNode2) )
                                      NumberIntersectionsOnTetCorner++;
				  //BY NOW I WANT TO CONSIDER ONLY THE EDGES THAT ARE CUT "NOT AT THE VERTEX"
				  else
					  {
					//	double dummy=0.0;
					  //(i_fluidElement->GetValue(NEIGHBOUR_EMBEDDED_FACES)).push_back( GeometricalObject::WeakPointer( *(i_StructCondition.base()) ) );

					  array_1d<double,3> emb_vel=(*i_StructCondition)->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);
					  emb_vel+=(*i_StructCondition)->GetGeometry()[1].FastGetSolutionStepValue(VELOCITY);
					  emb_vel+=(*i_StructCondition)->GetGeometry()[2].FastGetSolutionStepValue(VELOCITY);

					  //KRATOS_WATCH(emb_vel)

					  i_fluidElement->GetValue(EMBEDDED_VELOCITY)+=emb_vel;
					  intersection_counter++;

				  	  }
					  //(pGeom[i].GetValue(NEIGHBOUR_ELEMENTS)).push_back( Element::WeakPointer( *(ie.base()) ) );
                              }

                          }


                      }

                  }
              }

          }

          // check, if intersection nodes have been found on the tet edge --> if yes, then add these information to the TetEdgeVector
          if( NewTetEdge.IntNodes.size() > 0 )
              IntersectedTetEdges.push_back(NewTetEdge);
      }

      ///******************************************************************************************************************
      ///******************************************************************************************************************

      bool StructuralElementNotYetConsidered( unsigned int                IDCurrentStructCond,
                                              std::vector<unsigned int>&  IntersectingStructCondID )
      {
          // check if the structural element was already considered as intersecting element
          for(unsigned int k = 0 ; k < IntersectingStructCondID.size() ; k++)
          {
              if( IDCurrentStructCond == IntersectingStructCondID[k] )
                  return false;
          }

          // if structural element has not been considered in another octree, which also intersects the fluid element
          // add the new object ID to the vector
          IntersectingStructCondID.push_back( IDCurrentStructCond );
          return true;
      }

      ///******************************************************************************************************************
      ///******************************************************************************************************************

      bool IsIntersectionNodeOnTetEdge( double* IntersectionPoint , double* EdgeNode1 , double* EdgeNode2 )
      {
          // check, if intersection point is located on any edge of the fluid element
          array_1d<double,3> ConnectVectTetNodeIntNode1;
          array_1d<double,3> ConnectVectTetNodeIntNode2;
          array_1d<double,3> EdgeVector;

          ConnectVectTetNodeIntNode1[0] = IntersectionPoint[0] - EdgeNode1[0];
          ConnectVectTetNodeIntNode1[1] = IntersectionPoint[1] - EdgeNode1[1];
          ConnectVectTetNodeIntNode1[2] = IntersectionPoint[2] - EdgeNode1[2];

          ConnectVectTetNodeIntNode2[0] = IntersectionPoint[0] - EdgeNode2[0];
          ConnectVectTetNodeIntNode2[1] = IntersectionPoint[1] - EdgeNode2[1];
          ConnectVectTetNodeIntNode2[2] = IntersectionPoint[2] - EdgeNode2[2];

          double LengthConnectVect1 = norm_2( ConnectVectTetNodeIntNode1 );
          double LengthConnectVect2 = norm_2( ConnectVectTetNodeIntNode2 );

          EdgeVector[0] = EdgeNode2[0] - EdgeNode1[0];
          EdgeVector[1] = EdgeNode2[1] - EdgeNode1[1];
          EdgeVector[2] = EdgeNode2[2] - EdgeNode1[2];

          double MaxEdgeLength = norm_2( EdgeVector );

          // if both connection vectors (corner point --> intersection point)
          // are smaller or equal to the edge length of tetrahedra,
          // then intersection point is located on the edge
          if( (LengthConnectVect1 <= (MaxEdgeLength)) && (LengthConnectVect2 <= (MaxEdgeLength)) )
              return true;
          else
              return false;
      }

      ///******************************************************************************************************************
      ///******************************************************************************************************************

      bool IsNewIntersectionNode(IntersectionNodeStruct&    NewIntersectionNode,
                                 std::vector<TetEdgeStruct> IntersectedTetEdges )
      {
          array_1d<double,3> DiffVector;
          double NormDiffVector;
          unsigned int NumberIntNodes;

          for( unsigned int i_TetEdge = 0 ; i_TetEdge < IntersectedTetEdges.size() ; i_TetEdge++ )
          {
              NumberIntNodes = IntersectedTetEdges[i_TetEdge].IntNodes.size();
              for( unsigned int i_IntNode = 0 ; i_IntNode < NumberIntNodes ; i_IntNode++ )
              {
                  DiffVector[0] = NewIntersectionNode.Coordinates[0] - IntersectedTetEdges[i_TetEdge].IntNodes[i_IntNode].Coordinates[0];
                  DiffVector[1] = NewIntersectionNode.Coordinates[1] - IntersectedTetEdges[i_TetEdge].IntNodes[i_IntNode].Coordinates[1];
                  DiffVector[2] = NewIntersectionNode.Coordinates[2] - IntersectedTetEdges[i_TetEdge].IntNodes[i_IntNode].Coordinates[2];
                  NormDiffVector = norm_2(DiffVector);

                  if( NormDiffVector < epsilon )
                      return false;
              }
          }

          // if the new intersection node is not existing (as intersection with a corner point), then return false
          return true;
      }

      ///******************************************************************************************************************
      ///******************************************************************************************************************

      bool IsIntersectionOnCorner(IntersectionNodeStruct& NewIntersectionNode,
                                  double* EdgeNode1,
                                  double* EdgeNode2 )
      {
          array_1d<double,3> DiffVector;
          double NormDiffVector;

          DiffVector[0] = EdgeNode1[0] - NewIntersectionNode.Coordinates[0];
          DiffVector[1] = EdgeNode1[1] - NewIntersectionNode.Coordinates[1];
          DiffVector[2] = EdgeNode1[2] - NewIntersectionNode.Coordinates[2];
          NormDiffVector = norm_2(DiffVector);

          if( NormDiffVector < epsilon )
              return true;

          DiffVector[0] = EdgeNode2[0] - NewIntersectionNode.Coordinates[0];
          DiffVector[1] = EdgeNode2[1] - NewIntersectionNode.Coordinates[1];
          DiffVector[2] = EdgeNode2[2] - NewIntersectionNode.Coordinates[2];
          NormDiffVector = norm_2(DiffVector);

          if( NormDiffVector < epsilon )
              return true;
          else
              return false;
      }

      ///******************************************************************************************************************
      ///******************************************************************************************************************

      void CalculateNormal3D(Element::GeometryType& rGeometry,
                             array_1d<double,3>&    rResultNormal)
      {
          array_1d<double,3> v1;
          array_1d<double,3> v2 ;

          v1[0] = rGeometry[1].X() - rGeometry[0].X();
          v1[1] = rGeometry[1].Y() - rGeometry[0].Y();
          v1[2] = rGeometry[1].Z() - rGeometry[0].Z();

          v2[0] = rGeometry[2].X() - rGeometry[0].X();
          v2[1] = rGeometry[2].Y() - rGeometry[0].Y();
          v2[2] = rGeometry[2].Z() - rGeometry[0].Z();

          MathUtils<double>::CrossProduct(rResultNormal,v1,v2);
          rResultNormal *= 0.5;
      }

      ///******************************************************************************************************************
      ///******************************************************************************************************************

      void CalcNodalDistanceTo3DSkin(std::vector<TetEdgeStruct>&                 IntersectedTetEdges,
                                     ModelPart::ElementsContainerType::iterator& i_fluid_element,
                                     unsigned int                                NumberIntersectionsOnTetCorner)
      {
          std::vector<IntersectionNodeStruct> NodesOfApproximatedStructure;
          array_1d<double,4> ElementalDistances;

          // Reduce all found intersection nodes located on each tetdrahedra edge to just one intersection node by averaging
          ComputeApproximationNodes(IntersectedTetEdges,NodesOfApproximatedStructure);

          // Intersection with one corner point
          if( NodesOfApproximatedStructure.size() == 1 && NumberIntersectionsOnTetCorner == 1 )
          {
              CalcSignedDistancesToOneIntNode(i_fluid_element,NodesOfApproximatedStructure,ElementalDistances);
              i_fluid_element->GetValue(SPLIT_ELEMENT) = true;
          }

          // Intersection with two corner points / one tetrahedra edge
          if( NodesOfApproximatedStructure.size() == 2 && NumberIntersectionsOnTetCorner == 2 )
          {
              CalcSignedDistancesToTwoIntNodes(i_fluid_element,NodesOfApproximatedStructure,ElementalDistances);
              i_fluid_element->GetValue(SPLIT_ELEMENT) = true;
          }

          // Intersection with three tetrahedra edges
          if( NodesOfApproximatedStructure.size() == 3 )
          {
              CalcSignedDistancesToThreeIntNodes(i_fluid_element,NodesOfApproximatedStructure,IntersectedTetEdges,ElementalDistances);
              i_fluid_element->GetValue(SPLIT_ELEMENT) = true;
          }

          // Intersection with four tetrahedra edges
          if( NodesOfApproximatedStructure.size() == 4 )
          {
              CalcSignedDistancesToFourIntNodes(i_fluid_element,NodesOfApproximatedStructure,IntersectedTetEdges,ElementalDistances);
              i_fluid_element->GetValue(SPLIT_ELEMENT) = true;
          }

          // In case there is NO intersection with fluid element
          if( i_fluid_element->GetValue(SPLIT_ELEMENT) == true )
              AssignDistancesToElements(i_fluid_element,ElementalDistances);
      }

      ///******************************************************************************************************************
      ///******************************************************************************************************************

      void ComputeApproximationNodes(std::vector<TetEdgeStruct>           IntersectedTetEdges,
                                     std::vector<IntersectionNodeStruct>& NodesOfApproximatedStructure)
      {
          unsigned int NumberIntNodes;
          double sum_X;
          double sum_Y;
          double sum_Z;

          // calculate average of all intersection nodes of each tetrahedra edge
          for(unsigned int i_TetEdge = 0 ; i_TetEdge < IntersectedTetEdges.size() ; i_TetEdge++)
          {
              NumberIntNodes = IntersectedTetEdges[i_TetEdge].IntNodes.size();
              sum_X = 0;
              sum_Y = 0;
              sum_Z = 0;

              for( unsigned int i_IntNode = 0 ; i_IntNode < NumberIntNodes ; i_IntNode++ )
              {
                  sum_X += IntersectedTetEdges[i_TetEdge].IntNodes[i_IntNode].Coordinates[0];
                  sum_Y += IntersectedTetEdges[i_TetEdge].IntNodes[i_IntNode].Coordinates[1];
                  sum_Z += IntersectedTetEdges[i_TetEdge].IntNodes[i_IntNode].Coordinates[2];
              }

              IntersectionNodeStruct NewApproximationNode;
              NewApproximationNode.Coordinates[0] = sum_X / NumberIntNodes;
              NewApproximationNode.Coordinates[1] = sum_Y / NumberIntNodes;
              NewApproximationNode.Coordinates[2] = sum_Z / NumberIntNodes;

              if(IntersectedTetEdges.size() <= 2)
                  NewApproximationNode.StructElemNormal = IntersectedTetEdges[i_TetEdge].IntNodes[0].StructElemNormal;

              NodesOfApproximatedStructure.push_back(NewApproximationNode);
          }
      }

      ///******************************************************************************************************************
      ///******************************************************************************************************************

      void CalcSignedDistancesToOneIntNode( ModelPart::ElementsContainerType::iterator& i_fluid_element,
                                            std::vector<IntersectionNodeStruct>         NodesOfApproximatedStructure,
                                            array_1d<double,4>&                         ElementalDistances)
      {
          const array_1d<double,3>& IntersectionNodeCoord = NodesOfApproximatedStructure[0].Coordinates;
          array_1d<double,3> DistVecTetNode;
          array_1d<double,3> TetNode;
          array_1d<double,3> NormalAtIntersectionNode;
          double             NormDistTetNode;
          double             InnerProduct;

          Geometry< Node<3> >& rFluidGeom = i_fluid_element->GetGeometry();

          for(unsigned int i_TetNode = 0 ; i_TetNode < 4 ; i_TetNode++)
          {
              // Get coordinates of the fluid elmenent nodes
              TetNode  = rFluidGeom[i_TetNode].Coordinates();

              // Compute unsigned distance
              DistVecTetNode[0] = TetNode[0] - IntersectionNodeCoord[0];
              DistVecTetNode[1] = TetNode[1] - IntersectionNodeCoord[1];
              DistVecTetNode[2] = TetNode[2] - IntersectionNodeCoord[2];
              NormDistTetNode = norm_2( DistVecTetNode );

              // Get normal at intersection
              NormalAtIntersectionNode = NodesOfApproximatedStructure[0].StructElemNormal;
              InnerProduct = inner_prod(DistVecTetNode,NormalAtIntersectionNode);

              // Assign distances as nodal solution values
              if(InnerProduct>epsilon)
                  ElementalDistances[i_TetNode] = NormDistTetNode;
              else if(InnerProduct>-epsilon)
                  ElementalDistances[i_TetNode] = 0;
              else
                  ElementalDistances[i_TetNode] = -NormDistTetNode;
          }
      }

      ///******************************************************************************************************************
      ///******************************************************************************************************************

      void CalcSignedDistancesToTwoIntNodes( ModelPart::ElementsContainerType::iterator& i_fluid_element,
                                             std::vector<IntersectionNodeStruct>         NodesOfApproximatedStructure,
                                             array_1d<double,4>&                         ElementalDistances)
      {
          const array_1d<double,3>& IntersectionNode1Coord = NodesOfApproximatedStructure[0].Coordinates;
          const array_1d<double,3>& IntersectionNode2Coord = NodesOfApproximatedStructure[1].Coordinates;
          array_1d<double,3> TetNode;
          array_1d<double,3> DistVecTetNode;
          array_1d<double,3> NormalAtIntersectionNode1;
          array_1d<double,3> NormalAtIntersectionNode2;
          array_1d<double,3> ResNormal;
          double             InnerProduct;
          double             NormDistTetNode;

          const Point LinePoint1 = Point(IntersectionNode1Coord[0] , IntersectionNode1Coord[1] , IntersectionNode1Coord[2]);
          const Point LinePoint2 = Point(IntersectionNode2Coord[0] , IntersectionNode2Coord[1] , IntersectionNode2Coord[2]);

          Geometry< Node<3> >& rFluidGeom = i_fluid_element->GetGeometry();

          for(unsigned int i_TetNode = 0 ; i_TetNode < 4 ; i_TetNode++)
          {
              // Get coordinates of the fluid element nodes
              TetNode  = rFluidGeom(i_TetNode)->Coordinates();

              // Compute distance to point
              NormDistTetNode = GeometryUtils::PointDistanceToLineSegment3D(LinePoint1, LinePoint2 , Point(TetNode[0],TetNode[1],TetNode[2]));

              // Compute unsigned distance vector by assuming the mean position vector of the two intersection points
              DistVecTetNode[0] = TetNode[0] - IntersectionNode1Coord[0];
              DistVecTetNode[1] = TetNode[1] - IntersectionNode1Coord[1];
              DistVecTetNode[2] = TetNode[2] - IntersectionNode1Coord[2];

              // Get normal at intersections, average them and check direction of distances
              NormalAtIntersectionNode1 = NodesOfApproximatedStructure[0].StructElemNormal;
              NormalAtIntersectionNode2 = NodesOfApproximatedStructure[1].StructElemNormal;

              // Compute unsigned distance
              ResNormal[0] = 0.5*(NormalAtIntersectionNode1[0] + NormalAtIntersectionNode2[0]);
              ResNormal[1] = 0.5*(NormalAtIntersectionNode1[1] + NormalAtIntersectionNode2[1]);
              ResNormal[2] = 0.5*(NormalAtIntersectionNode1[2] + NormalAtIntersectionNode2[2]);
              InnerProduct = inner_prod(DistVecTetNode,ResNormal);

              // Assign distances as nodal solution values
              if(InnerProduct>epsilon)
                  ElementalDistances[i_TetNode] = NormDistTetNode;
              else if(InnerProduct>-epsilon)
                  ElementalDistances[i_TetNode] = 0;
              else
                  ElementalDistances[i_TetNode] = -NormDistTetNode;
          }
      }

      ///******************************************************************************************************************
      ///******************************************************************************************************************

      void CalcSignedDistancesToThreeIntNodes( ModelPart::ElementsContainerType::iterator& i_fluid_element,
                                               std::vector<IntersectionNodeStruct>         NodesOfApproximatedStructure,
                                               std::vector<TetEdgeStruct>                  IntersectedTetEdges,
                                               array_1d<double,4>&                         ElementalDistances)
      {
          array_1d<unsigned int,3> IndexNodes;

          IndexNodes[0] = 0;
          IndexNodes[1] = 1;
          IndexNodes[2] = 2;

          // Compute distance for each tetrahedra node to the triangle to approximate the structure
          CalcSignedDistancesToApproxTriangle( i_fluid_element , NodesOfApproximatedStructure , IntersectedTetEdges , ElementalDistances , IndexNodes );
      }

      ///******************************************************************************************************************
      ///******************************************************************************************************************

      void CalcSignedDistancesToFourIntNodes( ModelPart::ElementsContainerType::iterator& i_fluid_element,
                                              std::vector<IntersectionNodeStruct>         NodesOfApproximatedStructure,
                                              std::vector<TetEdgeStruct>                  IntersectedTetEdges,
                                              array_1d<double,4>&                         ElementalDistances)
      {
          array_1d<unsigned int,3> IndexNodes_T1; // nodes of first triangle
          array_1d<unsigned int,3> IndexNodes_T2; // nodes of second triangle

          array_1d<double,4> ElementalDistances_T1;
          array_1d<double,4> ElementalDistances_T2;

          double dist_T1;
          double dist_T2;

          // Generate 2 triangles out of the 4 int nodes which form a parallelogram together
          // Therefor first define arbitrarily a set of 3 nodes which form a triangle and search for the 2 nodes
          // of the triangle which form the other triangle together with the fourth node
          IndexNodes_T1[0] = 0;
          IndexNodes_T1[1] = 1;
          IndexNodes_T1[2] = 2;

          FindIndexNodesOfTriangle2(NodesOfApproximatedStructure,IndexNodes_T2);

          // Compute distance for first triangle
          CalcSignedDistancesToApproxTriangle( i_fluid_element , NodesOfApproximatedStructure , IntersectedTetEdges,
                                               ElementalDistances_T1 , IndexNodes_T1 );
          // Compute distance for second triangle
          CalcSignedDistancesToApproxTriangle( i_fluid_element , NodesOfApproximatedStructure , IntersectedTetEdges,
                                               ElementalDistances_T2 , IndexNodes_T2 );

          // Determine about the minimal distance by considering the distances to both triangles
          for(unsigned int i_TetNode = 0 ; i_TetNode < 4 ; i_TetNode++)
          {
              dist_T1 = ElementalDistances_T1[i_TetNode];
              dist_T2 = ElementalDistances_T2[i_TetNode];
              if(fabs(dist_T1) < fabs(dist_T2))
                  ElementalDistances[i_TetNode] = dist_T1;
              else
                  ElementalDistances[i_TetNode] = dist_T2;
          }
      }

      ///******************************************************************************************************************
      ///******************************************************************************************************************

      void FindIndexNodesOfTriangle2(std::vector<IntersectionNodeStruct> NodesOfApproximatedStructure,
                                     array_1d<unsigned int,3>&           IndexNodes_T2)
      {
          double maxDist = 0;
          unsigned int indexExcludedNode = 1000000; // index of the node which is not part of the second triangle
          array_1d<double,3> TrianglePoint;
          array_1d<double,3> RemainingPoint;
          array_1d<double,3> DistVecNode;

          RemainingPoint =  NodesOfApproximatedStructure[3].Coordinates;

          // strategy: these two nodes of the first triangle form a triangle with the 4th node, which are closest to that node
          // --> look for these nodes with the shortest distance to the remaining node
          for(unsigned int i_TriangleNode = 0 ; i_TriangleNode < 3 ; i_TriangleNode++)
          {
              TrianglePoint = NodesOfApproximatedStructure[i_TriangleNode].Coordinates;
              DistVecNode[0] = RemainingPoint[0] - TrianglePoint[0];
              DistVecNode[1] = RemainingPoint[1] - TrianglePoint[1];
              DistVecNode[2] = RemainingPoint[2] - TrianglePoint[2];
              if(norm_2(DistVecNode) > maxDist)
              {
                  maxDist = norm_2(DistVecNode);
                  indexExcludedNode = i_TriangleNode;
              }
          }

          // assign the "not excluded" nodes to the index vector of the second triangle
          unsigned int indexCursor = 0;
          for(unsigned int k = 0 ; k < 3 ; k++)
          {
              if(indexExcludedNode != k)
              {
                  IndexNodes_T2[indexCursor] = k;
                  indexCursor += 1;
              }
          }

          IndexNodes_T2[2] = 3;
      }

      ///******************************************************************************************************************
      ///******************************************************************************************************************

      void CalcSignedDistancesToApproxTriangle( ModelPart::ElementsContainerType::iterator& i_fluid_element,
                                                std::vector<IntersectionNodeStruct>         NodesOfApproximatedStructure,
                                                std::vector<TetEdgeStruct>                  IntersectedTetEdges,
                                                array_1d<double,4>&                         ElementalDistances,
                                                array_1d<unsigned int,3>                    IndexNodes)
      {
          Geometry< Node<3> >& rFluidGeom = i_fluid_element->GetGeometry();

          for(unsigned int i_TetNode = 0 ; i_TetNode < 4 ; i_TetNode++)
          {
              array_1d<double,3> TetNode;
              array_1d<double,3> IntersectionNode1Coord;
              array_1d<double,3> IntersectionNode2Coord;
              array_1d<double,3> IntersectionNode3Coord;
              Point           ApproxTrianglePoint1;
              Point           ApproxTrianglePoint2;
              Point           ApproxTrianglePoint3;
              double             UnsignedDistance;
              double             InnerProduct;
              unsigned int       IndexNode1;
              unsigned int       IndexNode2;
              unsigned int       IndexNode3;

              // Get coordinates of the fluid element nodes
              TetNode = rFluidGeom(i_TetNode)->Coordinates();

              IndexNode1 = IndexNodes[0];
              IndexNode2 = IndexNodes[1];
              IndexNode3 = IndexNodes[2];

              IntersectionNode1Coord = NodesOfApproximatedStructure[IndexNode1].Coordinates;
              IntersectionNode2Coord = NodesOfApproximatedStructure[IndexNode2].Coordinates;
              IntersectionNode3Coord = NodesOfApproximatedStructure[IndexNode3].Coordinates;

              ApproxTrianglePoint1 = Point(IntersectionNode1Coord[0] , IntersectionNode1Coord[1] , IntersectionNode1Coord[2]);
              ApproxTrianglePoint2 = Point(IntersectionNode2Coord[0] , IntersectionNode2Coord[1] , IntersectionNode2Coord[2]);
              ApproxTrianglePoint3 = Point(IntersectionNode3Coord[0] , IntersectionNode3Coord[1] , IntersectionNode3Coord[2]);

              // Compute distance from tet node to current triangle
              UnsignedDistance = GeometryUtils::PointDistanceToTriangle3D(ApproxTrianglePoint1, ApproxTrianglePoint2 , ApproxTrianglePoint3 , Point(TetNode[0],TetNode[1],TetNode[2]));

              bool TetNodeIsInsideStructure = true;
              bool TetNodeIsOnStructure = true;
              array_1d <double,3> DistVec;
              array_1d <double,3> NormalAtIntersectionNode;

              for( unsigned int i_TetEdge = 0 ; i_TetEdge < IntersectedTetEdges.size() ; i_TetEdge++ )
              {
                  for( unsigned int i_IntNode = 0 ; i_IntNode < IntersectedTetEdges[i_TetEdge].IntNodes.size() ; i_IntNode++ )
                  {
                      DistVec[0] = TetNode[0] - IntersectedTetEdges[i_TetEdge].IntNodes[i_IntNode].Coordinates[0];
                      DistVec[1] = TetNode[1] - IntersectedTetEdges[i_TetEdge].IntNodes[i_IntNode].Coordinates[1];
                      DistVec[2] = TetNode[2] - IntersectedTetEdges[i_TetEdge].IntNodes[i_IntNode].Coordinates[2];

                      NormalAtIntersectionNode = IntersectedTetEdges[i_TetEdge].IntNodes[i_IntNode].StructElemNormal;

                      InnerProduct = inner_prod(DistVec,NormalAtIntersectionNode);

                      if(InnerProduct > epsilon)
                      {
                          TetNodeIsInsideStructure = false;
                          TetNodeIsOnStructure = false;
                      }
                      else if (InnerProduct < -epsilon)
                          TetNodeIsOnStructure = false;
                  }
              }

              // Assign distances as nodal solution values ( + = outside of structure, - = inside structure)
              if( TetNodeIsInsideStructure == true )
                  ElementalDistances[i_TetNode] = -UnsignedDistance;
              else if( TetNodeIsOnStructure == true )
                  ElementalDistances[i_TetNode] = 0;
              else
                  ElementalDistances[i_TetNode] = +UnsignedDistance;
          }
      }

      ///******************************************************************************************************************
      ///******************************************************************************************************************

      void AssignDistancesToElements(ModelPart::ElementsContainerType::iterator& i_fluid_element,
                                     array_1d<double,4>                          ElementalDistances)
      {
          Geometry< Node<3> >& rFluidGeom = i_fluid_element->GetGeometry();

          array_1d<double,4> MinElementalDistances;

          for(unsigned int i_TetNode = 0 ; i_TetNode < 4 ; i_TetNode++)
          {
              //Assign distances to the element, if a smaller value could be found
              if( fabs(ElementalDistances[i_TetNode]) < fabs(i_fluid_element->GetValue(ELEMENTAL_DISTANCES)[i_TetNode]) )
                  MinElementalDistances[i_TetNode] = ElementalDistances[i_TetNode];
              else
                  MinElementalDistances[i_TetNode] = i_fluid_element->GetValue(ELEMENTAL_DISTANCES)[i_TetNode];

              //Assign distances to the single nodes (for visualization), if a smaller value could be found
              if( fabs(ElementalDistances[i_TetNode]) < fabs(rFluidGeom[i_TetNode].GetSolutionStepValue(DISTANCE)) )
                  rFluidGeom[i_TetNode].GetSolutionStepValue(DISTANCE) = ElementalDistances[i_TetNode];
          }

          i_fluid_element->SetValue(ELEMENTAL_DISTANCES,MinElementalDistances);
      }

  ///******************************************************************************************************************
  ///******************************************************************************************************************



	/*
      void GenerateOctree()
      {
          Timer::Start("Generating Octree");

	   // Setting the boundingbox for non-normalized coordinates
	   const int dimension = 3;
	   double boundingBox_low[3],boundingBox_high[3];

	   for(int i = 0; i < dimension; i++)
	    {
		boundingBox_low[i]  = mrSkinModelPart.NodesBegin()->Coordinates()[i];
		boundingBox_high[i] = mrSkinModelPart.NodesBegin()->Coordinates()[i];
	    }

	   for(ModelPart::NodeIterator i_node = mrSkinModelPart.NodesBegin();
        	i_node != mrSkinModelPart.NodesEnd();
        	i_node++)
	    {
		for(int i = 0; i < dimension; i++)
		{
		    if(i_node->Coordinates()[i] < boundingBox_low[i])  boundingBox_low[i]  = i_node->Coordinates()[i];
		    if(i_node->Coordinates()[i] > boundingBox_high[i]) boundingBox_high[i] = i_node->Coordinates()[i];
		}
    	    }

	  mOctree.SetBoundingBox(boundingBox_low,boundingBox_high);


          //mOctree.RefineWithUniformSize(0.0625);
          for(ModelPart::NodeIterator i_node = mrSkinModelPart.NodesBegin() ; i_node != mrSkinModelPart.NodesEnd() ; i_node++)
          {
              double temp_point[3];
              const array_1d<double,3>& r_coordinates = i_node->Coordinates();
              temp_point[0] = r_coordinates[0];
              temp_point[1] = r_coordinates[1];
              temp_point[2] = r_coordinates[2];
              mOctree.Insert(temp_point);
          }


          //mOctree.Constrain2To1(); // To be removed. Pooyan.

          for(ModelPart::ElementIterator i_element = mrSkinModelPart.ElementsBegin() ; i_element != mrSkinModelPart.ElementsEnd() ; i_element++)
          {
              mOctree.Insert(*(i_element).base());
          }

          Timer::Stop("Generating Octree");
//          octree.Insert(*(mrSkinModelPart.ElementsBegin().base()));
          KRATOS_WATCH(mOctree);

//          std::cout << "######## WRITING OCTREE MESH #########" << std::endl;
//          std::ofstream myfile;
//          myfile.open ("unserbaum.post.msh");
//          mOctree.PrintGiDMesh(myfile);
//          myfile.close();
      }
*/

	void GenerateOctree()
    {
        Timer::Start("Generating Octree");
        //std::cout << "Generating the Octree..." << std::endl;

        double low[3];
        double high[3];

        for (int i = 0 ; i < 3; i++)
        {
            low[i] = high[i] = mrFluidModelPart.NodesBegin()->Coordinates()[i];
        }

        // loop over all structure nodes
        for(ModelPart::NodeIterator i_node = mrFluidModelPart.NodesBegin();
            i_node != mrFluidModelPart.NodesEnd();
            i_node++)
        {
            const array_1d<double,3>& r_coordinates = i_node->Coordinates();
            for (int i = 0 ; i < 3; i++)
            {
                low[i]  = r_coordinates[i] < low[i]  ? r_coordinates[i] : low[i];
                high[i] = r_coordinates[i] > high[i] ? r_coordinates[i] : high[i];
            }
        }
        for (ModelPart::NodeIterator i_node = mrSkinModelPart.NodesBegin();
             i_node != mrSkinModelPart.NodesEnd();
             i_node++)
        {
            for (int i = 0; i < 3; i++)
            {
                low[i] = i_node->Coordinates()[i] < low[i] ? i_node->Coordinates()[i] : low[i];
                high[i] = i_node->Coordinates()[i] > high[i] ? i_node->Coordinates()[i] : high[i];
            }
        }
// KRATOS_WATCH( low[0] )
// KRATOS_WATCH( low[1] )
// KRATOS_WATCH( low[2] )
// KRATOS_WATCH( "" )
// KRATOS_WATCH( high[0] )
// KRATOS_WATCH( high[1] )
// KRATOS_WATCH( high[2] )
        mOctree.SetBoundingBox(low,high);

        //mOctree.RefineWithUniformSize(0.0625);

        // loop over all structure nodes
        for(ModelPart::NodeIterator i_node = mrSkinModelPart.NodesBegin();
            i_node != mrSkinModelPart.NodesEnd();
            i_node++)
        {
            double temp_point[3];
            temp_point[0] = i_node->X();
            temp_point[1] = i_node->Y();
            temp_point[2] = i_node->Z();
            mOctree.Insert(temp_point);
        }

        //mOctree.Constrain2To1(); // To be removed. Pooyan.

        // loop over all structure elements
        //for(ModelPart::ElementIterator i_element = mrSkinModelPart.ElementsBegin();
         //   i_element != mrSkinModelPart.ElementsEnd();
          //  i_element++)
	for(ModelPart::ConditionIterator i_cond = mrSkinModelPart.ConditionsBegin();
            i_cond != mrSkinModelPart.ConditionsEnd();
            i_cond++)
        {
            mOctree.Insert(*(i_cond).base());
        }

        Timer::Stop("Generating Octree");
        std::cout << "######## Octree generated #########" << std::endl;

//        KRATOS_WATCH(mOctree);

//        std::cout << "######## WRITING OCTREE MESH #########" << std::endl;
//        std::ofstream myfile;
//        myfile.open ("octree.post.msh");
//        mOctree.PrintGiDMesh(myfile);
//        myfile.close();

        //std::cout << "Generating the Octree finished" << std::endl;
    }


      ///******************************************************************************************************************
      ///******************************************************************************************************************

      ///******************************************************************************************************************
      ///******************************************************************************************************************


      void GenerateNodes()
      {
          Timer::Start("Generating Nodes");
          std::vector<OctreeType::cell_type*> all_leaves;
          mOctree.GetAllLeavesVector(all_leaves);

#pragma omp parallel for
          for (int i = 0; i < static_cast<int>(all_leaves.size()); i++)
          {
              *(all_leaves[i]->pGetDataPointer()) = ConfigurationType::AllocateData();
          }


          std::size_t last_id = mrBodyModelPart.NumberOfNodes() + 1;
          //KRATOS_WATCH(all_leaves.size());
          for (std::size_t i = 0; i < all_leaves.size(); i++)
          {
              //KRATOS_WATCH(i)
                CellType* cell = all_leaves[i];
                GenerateCellNode(cell, last_id);
          }

          Timer::Stop("Generating Nodes");

      }

      void GenerateCellNode(CellType* pCell, std::size_t& LastId)
      {
        for (int i_pos=0; i_pos < 8; i_pos++) // position 8 is for center
        {
            ConfigurationType::cell_node_data_type* p_node = (*(pCell->pGetData()))[i_pos];
            if(p_node == 0)
            {
                (*(pCell->pGetData()))[i_pos] = new ConfigurationType::cell_node_data_type;

                (*(pCell->pGetData()))[i_pos]->Id() = LastId++;
                //KRATOS_WATCH(LastId)

                mOctreeNodes.push_back((*(pCell->pGetData()))[i_pos]);

                SetNodeInNeighbours(pCell,i_pos,(*(pCell->pGetData()))[i_pos]);
            }

        }
      }

     void SetNodeInNeighbours(CellType* pCell, int Position, CellNodeDataType* pNode)
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
                        std::cout << "ERROR!! Bad Position calculated!!!!!!!!!!! position :" << position << std::endl;
                        continue;
                    }

                    (*neighbour_cell->pGetData())[position] = pNode;
                }
            }
      }


     void CalculateDistance2()
     {
         Timer::Start("Calculate Distances2");
         ModelPart::NodesContainerType::ContainerType& nodes = mrFluidModelPart.NodesArray();
         int nodes_size = nodes.size();
//         // first of all we reset the node distance to 1.00 which is the maximum distnace in our normalized space.
//#pragma omp parallel for firstprivate(nodes_size)
//         for(int i = 0 ; i < nodes_size ; i++)
//             nodes[i]->GetSolutionStepValue(DISTANCE) = 1.00;

           std::vector<CellType*> leaves;

         mOctree.GetAllLeavesVector(leaves);
         //int leaves_size = leaves.size();

//         for(int i = 0 ; i < leaves_size ; i++)
//             CalculateNotEmptyLeavesDistance(leaves[i]);

#pragma omp parallel for firstprivate(nodes_size)
         for(int i = 0 ; i < nodes_size ; i++)
         {
             CalculateNodeDistance(*(nodes[i]));
         }
         Timer::Stop("Calculate Distances2");

     }
//     void CalculateDistance3()
//     {
//         Timer::Start("Calculate Distances2");
//         ModelPart::NodesContainerType::ContainerType& nodes = mrFluidModelPart.NodesArray();
//         int nodes_size = nodes.size();
////         // first of all we reset the node distance to 1.00 which is the maximum distnace in our normalized space.
//#pragma omp parallel for firstprivate(nodes_size)
//         for(int i = 0 ; i < nodes_size ; i++)
//             nodes[i]->GetSolutionStepValue(DISTANCE) = 1.00;

//           std::vector<CellType*> leaves;

//         mOctree.GetAllLeavesVector(leaves);
//         int leaves_size = leaves.size();

//         for(int i = 0 ; i < leaves_size ; i++)
//             CalculateNotEmptyLeavesDistance(leaves[i]);

//#pragma omp parallel for firstprivate(nodes_size)
//         for(int i = 0 ; i < nodes_size ; i++)
//         {
//             CalculateNodeDistance(*(nodes[i]));
//         }
//         Timer::Stop("Calculate Distances2");

//     }
//     void CalculateDistance4()
//     {
//         Timer::Start("Calculate Distances3");
//         ModelPart::NodesContainerType::ContainerType& nodes = mrFluidModelPart.NodesArray();
//         int nodes_size = nodes.size();
//           std::vector<CellType*> leaves;

//         mOctree.GetAllLeavesVector(leaves);
//         int leaves_size = leaves.size();

//#pragma omp parallel for firstprivate(nodes_size)
//         for(int i = 0 ; i < nodes_size ; i++)
//         {
//             CalculateNodeDistanceFromCell(*(nodes[i]));
//         }
//         Timer::Stop("Calculate Distances3");

//     }


      void CalculateDistance()
      {
          Timer::Start("Calculate Distances");
          ConfigurationType::data_type& nodes = mOctreeNodes;
          int nodes_size = nodes.size();
          // first of all we reste the node distance to 1.00 which is the maximum distnace in our normalized space.
#pragma omp parallel for firstprivate(nodes_size)
          for(int i = 0 ; i < nodes_size ; i++)
              nodes[i]->Distance() = 1.00;


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
          Timer::Stop("Calculate Distances");

      }

      void CalculateDistance(CellNodeDataType& rNode, int i_direction)
      {
          double coords[3] = {rNode.X(), rNode.Y(), rNode.Z()};
         // KRATOS_WATCH_3(coords);

             //This function must color the positions in space defined by 'coords'.
            //coords is of dimension (3) normalized in (0,1)^3 space

            typedef Element::GeometryType triangle_type;
            typedef std::vector<std::pair<double, triangle_type*> > intersections_container_type;

            intersections_container_type intersections;
            ConfigurationType::data_type nodes_array;


            const double epsilon = 1e-12;

            double distance = 1.0;

            // Creating the ray
            double ray[3] = {coords[0], coords[1], coords[2]};
            ray[i_direction] = 0; // starting from the lower extreme

//            KRATOS_WATCH_3(ray)
            GetIntersectionsAndNodes(ray, i_direction, intersections, nodes_array);
//            KRATOS_WATCH(nodes_array.size())
            for (std::size_t i_node = 0; i_node < nodes_array.size() ; i_node++)
            {
                double coord = (*nodes_array[i_node])[i_direction];
   //             KRATOS_WATCH(intersections.size());

                int ray_color= 1;
                std::vector<std::pair<double, Element::GeometryType*> >::iterator i_intersection = intersections.begin();
                while (i_intersection != intersections.end()) {
                    double d = coord - i_intersection->first;
                    if (d > epsilon) {

                        ray_color = -ray_color;
                        distance = d;
                    } else if (d > -epsilon) {//interface
                        distance = 0.00;
                        break;
                    } else {
                        if(distance > -d)
                            distance = -d;
                        break;
                    }

                    i_intersection++;
                }

                distance *= ray_color;

                double& node_distance = nodes_array[i_node]->Distance();
                if(fabs(distance) < fabs(node_distance))
                    node_distance = distance;
                else if (distance*node_distance < 0.00) // assigning the correct sign
                    node_distance = -node_distance;


            }
     }

      void CalculateNotEmptyLeavesDistance(CellType* pCell)
      {
        //typedef Element::GeometryType triangle_type;
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
                mOctree.CalculateCoordinates(keys,cell_point);

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
          double& node_distance =  rNode.GetSolutionStepValue(DISTANCE);

          //const double epsilon = 1.00e-12;
          if(fabs(node_distance) > fabs(distance))
            node_distance = distance;
          else if (distance*node_distance < 0.00) // assigning the correct sign
              node_distance = -node_distance;
      }

//      void CalculateNodeDistanceFromCell(Node<3>& rNode)
//      {
//          OctreeType::key_type node_key[3] = {octree->CalcKeyNormalized(rNode.X()), octree->CalcKeyNormalized(rNode.Y()), octree->CalcKeyNormalized(rNode.Z())};
//          OctreeType::cell_type* pcell = octree->pGetCell(node_key);

//          object_container_type* objects = (pCell->pGetObjects());

//          // We interpolate the cell distances for the node in empty cells
//          if (objects->empty())
//          {

//          }

//          double distance = DistancePositionInSpace(coord);
//          double& node_distance =  rNode.GetSolutionStepValue(DISTANCE);

//          //const double epsilon = 1.00e-12;
//          if(fabs(node_distance) > fabs(distance))
//            node_distance = distance;
//          else if (distance*node_distance < 0.00) // assigning the correct sign
//              node_distance = -node_distance;
//      }

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
				mOctree.NormalizeCoordinates(ray);
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

            double distance = (fabs(distances[0]) > fabs(distances[1])) ? distances[1] : distances[0];
            distance = (fabs(distance) > fabs(distances[2])) ? distances[2] : distance;

            return distance;

        }


  void GetIntersectionsAndNodes(double* ray, int direction, std::vector<std::pair<double,Element::GeometryType*> >& intersections, ConfigurationType::data_type& rNodesArray)
  {
    //This function passes the ray through the model and gives the hit point to all objects in its way
    //ray is of dimension (3) normalized in (0,1)^3 space
    // direction can be 0,1,2 which are x,y and z respectively

    const double epsilon = 1.00e-12;

    // first clearing the intersections points vector
    intersections.clear();

    OctreeType* octree = &mOctree;

    OctreeType::key_type ray_key[3] = {octree->CalcKeyNormalized(ray[0]), octree->CalcKeyNormalized(ray[1]), octree->CalcKeyNormalized(ray[2])};
    OctreeType::key_type cell_key[3];

    // getting the entrance cell from lower extreme
    ray_key[direction] = 0;
    OctreeType::cell_type* cell = octree->pGetCell(ray_key);

    while (cell) {
        std::size_t position = cell->GetLocalPosition(ray_key); // Is this the local position!?!?!?!
        OctreeType::key_type node_key[3];
        cell->GetKey(position, node_key);
        if((node_key[0] == ray_key[0]) && (node_key[1] == ray_key[1]) && (node_key[2] == ray_key[2]))
        {
            if(cell->pGetData())
            {
                if(cell->pGetData()->size() > position)
                {
                    CellNodeDataType* p_node = (*cell->pGetData())[position];
                    if(p_node)
                    {
                        //KRATOS_WATCH(p_node->Id())
                        rNodesArray.push_back(p_node);
                    }
                }
                else
                    KRATOS_WATCH(cell->pGetData()->size())
            }
        }


//        std::cout << ".";
      GetCellIntersections(cell, ray, ray_key, direction, intersections);

      // Add the cell's middle node if existed
//      cell->GetKey(8, cell_key); // 8 is the central position
//      ray_key[direction]=cell_key[direction]; // positioning the ray in the middle of cell in its direction

//      position = cell->GetLocalPosition(ray_key);
//      if(position < 27) // principal nodes
//      {
//          if(cell->pGetData())
//          {
//              if(cell->pGetData()->size() > position)
//              {
//                  Node<3>* p_node = (*cell->pGetData())[position];
//                  if(p_node)
//                  {
//                      //KRATOS_WATCH(p_node->Id())
//                      rNodesArray.push_back(p_node);
//                  }
//              }
//              else
//                  KRATOS_WATCH(cell->pGetData()->size())
//          }
//      }
//      else
//      {
//          KRATOS_WATCH(position);
//          KRATOS_WATCH(*cell);
//      }


      // go to the next cell
      if (cell->GetNeighbourKey(1 + direction * 2, cell_key)) {
        ray_key[direction] = cell_key[direction];
        cell = octree->pGetCell(ray_key);
        ray_key[direction] -= 1 ;//the key returned by GetNeighbourKey is inside the cell (minkey +1), to ensure that the corresponding
        //cell get in pGetCell is the right one.
      } else
        cell = NULL;
    }



 //   KRATOS_WATCH(rNodesArray.size());
    // now eliminating the repeated objects
    if (!intersections.empty()) {
      //sort
      std::sort(intersections.begin(), intersections.end());
      // unique
      std::vector<std::pair<double, Element::GeometryType*> >::iterator i_begin = intersections.begin();
      std::vector<std::pair<double, Element::GeometryType*> >::iterator i_intersection = intersections.begin();
      while (++i_begin != intersections.end()) {
          // considering the very near points as the same points
          if (fabs(i_begin->first - i_intersection->first) > epsilon) // if the hit points are far enough they are not the same
            *(++i_intersection) = *i_begin;
      }
      intersections.resize((++i_intersection) - intersections.begin());

    }
  }

  void GetIntersections(double* ray, int direction, std::vector<std::pair<double,Element::GeometryType*> >& intersections)
  {
    //This function passes the ray through the model and gives the hit point to all objects in its way
    //ray is of dimension (3) normalized in (0,1)^3 space
    // direction can be 0,1,2 which are x,y and z respectively

    const double epsilon = 1.00e-12;

    // first clearing the intersections points vector
    intersections.clear();

    OctreeType* octree = &mOctree;

    OctreeType::key_type ray_key[3] = {octree->CalcKeyNormalized(ray[0]), octree->CalcKeyNormalized(ray[1]), octree->CalcKeyNormalized(ray[2])};
    OctreeType::key_type cell_key[3];

    // getting the entrance cell from lower extreme
    OctreeType::cell_type* cell = octree->pGetCell(ray_key);

    while (cell) {
//        std::cout << ".";
      GetCellIntersections(cell, ray, ray_key, direction, intersections);
      // go to the next cell
      if (cell->GetNeighbourKey(1 + direction * 2, cell_key)) {
        ray_key[direction] = cell_key[direction];
        cell = octree->pGetCell(ray_key);
        ray_key[direction] -= 1 ;//the key returned by GetNeighbourKey is inside the cell (minkey +1), to ensure that the corresponding
        //cell get in pGetCell is the right one.
      } else
        cell = NULL;
    }


    // now eliminating the repeated objects
    if (!intersections.empty()) {
      //sort
      std::sort(intersections.begin(), intersections.end());
      // unique
      std::vector<std::pair<double, Element::GeometryType*> >::iterator i_begin = intersections.begin();
      std::vector<std::pair<double, Element::GeometryType*> >::iterator i_intersection = intersections.begin();
      while (++i_begin != intersections.end()) {
          // considering the very near points as the same points
          if (fabs(i_begin->first - i_intersection->first) > epsilon) // if the hit points are far enough they are not the same
            *(++i_intersection) = *i_begin;
      }
      intersections.resize((++i_intersection) - intersections.begin());

    }
  }

  int GetCellIntersections(OctreeType::cell_type* cell, double* ray,
    OctreeType::key_type* ray_key, int direction,
    std::vector<std::pair<double, Element::GeometryType*> >& intersections)  {
      //This function passes the ray through the cell and gives the hit point to all objects in its way
      //ray is of dimension (3) normalized in (0,1)^3 space
      // direction can be 0,1,2 which are x,y and z respectively

      //typedef Element::GeometryType triangle_type;
      typedef OctreeType::cell_type::object_container_type object_container_type;

      object_container_type* objects = (cell->pGetObjects());

      // There are no intersection in empty cells
      if (objects->empty())
        return 0;

//      std::cout << "X";
      // calculating the two extreme of the ray segment inside the cell
      double ray_point1[3] = {ray[0], ray[1], ray[2]};
      double ray_point2[3] = {ray[0], ray[1], ray[2]};
      double normalized_coordinate;
      mOctree.CalculateCoordinateNormalized(ray_key[direction], normalized_coordinate);
      ray_point1[direction] = normalized_coordinate;
      ray_point2[direction] = ray_point1[direction] + mOctree.CalcSizeNormalized(cell);

      mOctree.ScaleBackToOriginalCoordinate(ray_point1);
      mOctree.ScaleBackToOriginalCoordinate(ray_point2);

      for (object_container_type::iterator i_object = objects->begin(); i_object != objects->end(); i_object++) {
        double intersection[3]={0.00,0.00,0.00};

        int is_intersected = IntersectionTriangleSegment((*i_object)->GetGeometry(), ray_point1, ray_point2, intersection); // This intersection has to be optimized for axis aligned rays

        if (is_intersected == 1) // There is an intersection but not coplanar
          intersections.push_back(std::pair<double, Element::GeometryType*>(intersection[direction], &((*i_object)->GetGeometry())));
        //else if(is_intersected == 2) // coplanar case
      }

      return 0;
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
      std::string Info() const override
	{
	  return "CalculateSignedDistanceTo3DConditionSkinProcess";
	}

      /// Print information about this object.
      void PrintInfo(std::ostream& rOStream) const override
	{
	  rOStream << "CalculateSignedDistanceTo3DConditionSkinProcess";
	}

      /// Print object's data.
      void PrintData(std::ostream& rOStream) const override
	{
	}

        void PrintGiDMesh(std::ostream & rOStream) const {
            std::vector<CellType*> leaves;

            mOctree.GetAllLeavesVector(leaves);

            std::cout << "writing " << leaves.size() << " leaves" << std::endl;
            rOStream << "MESH \"leaves\" dimension 3 ElemType Hexahedra Nnode 8" << std::endl;
            rOStream << "# color 96 96 96" << std::endl;
            rOStream << "Coordinates" << std::endl;
            rOStream << "# node number coordinate_x coordinate_y coordinate_z  " << std::endl;

          for(ConfigurationType::data_type::const_iterator i_node = mOctreeNodes.begin() ; i_node != mOctreeNodes.end() ; i_node++)
          {
              rOStream << (*i_node)->Id() << "  " << (*i_node)->X() << "  " << (*i_node)->Y() << "  " << (*i_node)->Z() << std::endl;
              //mOctree.Insert(temp_point);
          }
            std::cout << "Nodes written..." << std::endl;
            rOStream << "end coordinates" << std::endl;
            rOStream << "Elements" << std::endl;
            rOStream << "# element node_1 node_2 node_3 material_number" << std::endl;

            for (std::size_t i = 0; i < leaves.size(); i++) {
                if ((leaves[i]->pGetData()))
                {
                    ConfigurationType::data_type& nodes = (*(leaves[i]->pGetData()));

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

            for(ConfigurationType::data_type::const_iterator i_node = mOctreeNodes.begin() ; i_node != mOctreeNodes.end() ; i_node++)
          {
                rOStream << (*i_node)->Id() << "  " << (*i_node)->Distance() << std::endl;
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
        ModelPart& mrFluidModelPart;

        ConfigurationType::data_type mOctreeNodes;

        OctreeType mOctree;

        static const double epsilon;

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
      CalculateSignedDistanceTo3DConditionSkinProcess& operator=(CalculateSignedDistanceTo3DConditionSkinProcess const& rOther);

      /// Copy constructor.
      //CalculateSignedDistanceTo3DConditionSkinProcess(CalculateSignedDistanceTo3DConditionSkinProcess const& rOther);


      ///@}

    }; // Class CalculateSignedDistanceTo3DConditionSkinProcess

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream,
				    CalculateSignedDistanceTo3DConditionSkinProcess& rThis);

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream,
				    const CalculateSignedDistanceTo3DConditionSkinProcess& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }
  ///@}

  const double CalculateSignedDistanceTo3DConditionSkinProcess::epsilon = 1e-12;


}  // namespace Kratos.

#endif // KRATOS_CALCULATE_DISTANCE_CONDITION_PROCESS_H_INCLUDED  defined


