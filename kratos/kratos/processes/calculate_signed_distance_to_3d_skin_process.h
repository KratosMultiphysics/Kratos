/*
==============================================================================
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain

Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNER.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
*/
 
//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: pooyan $
//   Date:                $Date: 2007-10-31 17:51:34 $
//   Revision:            $Revision: 1.1 $
//
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
#include "geometries/triangle_3d_3.h"
#include "utilities/body_normal_calculation_utils.h"


namespace Kratos
{


class DistanceSpatialContainersConfigure
{

    class CellNodeData
    {
        double mDistance;
    public:
        double& Distance(){return mDistance;}
    };

public:


    ///@name Type Definitions
    ///@{

    enum { Dimension = 3,
           DIMENSION = 3,
           MAX_LEVEL = 10,
           MIN_LEVEL = 2
         };
    typedef Point<3, double>                                PointType;  /// always the point 3D
    typedef std::vector<double>::iterator                   DistanceIteratorType;
    typedef ModelPart::ElementsContainerType::ContainerType ContainerType;
    typedef ContainerType::value_type                       PointerType;
    typedef ContainerType::iterator                         IteratorType;
    typedef ModelPart::ElementsContainerType::ContainerType ResultContainerType;
    typedef ResultContainerType::value_type                 ResultPointerType;
    typedef ResultContainerType::iterator                   ResultIteratorType;

    typedef Element::Pointer                                        pointer_type;
    typedef CellNodeData                cell_node_data_type;
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
        typedef DistanceSpatialContainersConfigure::cell_node_data_type CellNodeDataType;
        typedef Point<3, double> PointType;  /// always the point 3D
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
      CalculateSignedDistanceTo3DSkinProcess(ModelPart& rThisModelPartStruc, ModelPart& rThisModelPartFluid)
          : mrSkinModelPart(rThisModelPartStruc), mrBodyModelPart(rThisModelPartStruc), mrFluidModelPart(rThisModelPartFluid)
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

      ///******************************************************************************************************************
      ///******************************************************************************************************************

      virtual void Execute()
      {
          KRATOS_TRY

          GenerateOctree();

          DistanceFluidStructure();

          KRATOS_CATCH("");
      }

      ///******************************************************************************************************************
      ///******************************************************************************************************************

      void DistanceFluidStructure()
      {
          // Initialize nodal distances each node in the domain to 1.0
          InitializeDistances();

          // Initialize index table to define line Edges of fluid element
          bounded_matrix<unsigned int,6,2> TetEdgeIndexTable;
          SetIndexTable(TetEdgeIndexTable);

          // loop over all fluid elements
          for( ModelPart::ElementIterator i_fluidElement = mrFluidModelPart.ElementsBegin();
                                          i_fluidElement != mrFluidModelPart.ElementsEnd();
                                          i_fluidElement++)
          {            
              CalcNodalDistancesOfTetNodes( i_fluidElement , TetEdgeIndexTable );
          }
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

          ModelPart::ElementsContainerType::ContainerType& elements = mrFluidModelPart.ElementsArray();

          array_1d<double,4> ElementalDistances;
          const double initial_distance = 1.0;
          ElementalDistances[0] = initial_distance;
          ElementalDistances[1] = initial_distance;
          ElementalDistances[2] = initial_distance;
          ElementalDistances[3] = initial_distance;

          // reset the elemental distance to 1.0 which is the maximum distance in our normalized space.
          for(unsigned int i = 0 ; i < elements.size() ; i++)
              elements[i]->SetValue(ELEMENTAL_DISTANCES,ElementalDistances);
      }      

      ///******************************************************************************************************************
      ///******************************************************************************************************************

      void SetIndexTable( bounded_matrix<unsigned int,6,2>& TetEdgeIndexTable )
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
                                         bounded_matrix<unsigned int,6,2>             TetEdgeIndexTable)
      {
          std::vector<OctreeType::cell_type*> leaves;
          std::vector<TetEdgeStruct>          IntersectedTetEdges;
          unsigned int NumberIntersectionsOnTetCorner = 0;

          // Get leaves of octree intersecting with fluid element
          mOctree.GetIntersectedLeaves(*(i_fluidElement).base(),leaves);

          // Loop over all 6 line Edges of the tetrahedra
          for(unsigned int i_tetEdge = 0; i_tetEdge < 6; i_tetEdge++)
          {
              IdentifyIntersectionNodes( i_fluidElement , i_tetEdge , leaves , IntersectedTetEdges ,
                                         NumberIntersectionsOnTetCorner , TetEdgeIndexTable );  
          }

          CalcNodalDistanceTo3DSkin( IntersectedTetEdges , i_fluidElement , NumberIntersectionsOnTetCorner );
      }

      ///******************************************************************************************************************
      ///******************************************************************************************************************

      void IdentifyIntersectionNodes( ModelPart::ElementsContainerType::iterator&   i_fluidElement,
                                      unsigned int                                  i_tetEdge,
                                      std::vector<OctreeType::cell_type*>&          leaves,
                                      std::vector<TetEdgeStruct>&                   IntersectedTetEdges,
                                      unsigned int                                  NumberIntersectionsOnTetCorner,
                                      bounded_matrix<unsigned int,6,2>               TetEdgeIndexTable)
      {
          std::vector<unsigned int> IntersectingStructElemID;
          TetEdgeStruct             NewTetEdge;

          // Get nodes of line Edge
          unsigned int EdgeStartIndex = TetEdgeIndexTable(i_tetEdge,0);
          unsigned int EdgeEndIndex   = TetEdgeIndexTable(i_tetEdge,1);

          PointType& P1 = i_fluidElement->GetGeometry()[EdgeStartIndex];
          PointType& P2 = i_fluidElement->GetGeometry()[EdgeEndIndex];

          double EdgeNode1[3] = {P1.X() , P1.Y() , P1.Z()};
          double EdgeNode2[3] = {P2.X() , P2.Y() , P2.Z()};

          // loop over all octree cells which are intersected by the fluid element
          for(unsigned int i_cell = 0 ; i_cell < leaves.size() ; i_cell++)
          {
              // Structural element contained in one cell of the octree
              object_container_type* struct_elem = (leaves[i_cell]->pGetObjects());

              // loop over all structural elements within each octree cell
              for(object_container_type::iterator i_StructElement = struct_elem->begin(); i_StructElement != struct_elem->end(); i_StructElement++)
              {

                  if( StructuralElementNotYetConsidered( (*i_StructElement)->Id() , IntersectingStructElemID ) )
                  {

                      // Calculate and associate intersection point to the current fluid element
                      double IntersectionPoint[3] = {0.0 , 0.0 , 0.0};
                      int TetEdgeHasIntersections = IntersectionTriangleSegment( (*i_StructElement)->GetGeometry() , EdgeNode1 , EdgeNode2 , IntersectionPoint );

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
                                  CalculateNormal3D((*i_StructElement)->GetGeometry(),NewIntersectionNode.StructElemNormal);

                                  // add the new intersection point to the list of intersection points of the fluid element
                                  NewTetEdge.IntNodes.push_back(NewIntersectionNode);

                                  // check, how many intersection nodes are located on corner points of the tetrahedra
                                  if ( IsIntersectionOnCorner( NewIntersectionNode, i_fluidElement->GetGeometry() ) )
                                      NumberIntersectionsOnTetCorner++;

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

      bool StructuralElementNotYetConsidered( unsigned int                IDCurrentStructElem,
                                              std::vector<unsigned int>&  IntersectingStructElemID )
      {
          // check if the structural element was already considered as intersecting element
          for(unsigned int k = 0 ; k < IntersectingStructElemID.size() ; k++)
          {
              if( IDCurrentStructElem == IntersectingStructElemID[k] )
                  return false;
          }

          // if structural element has not been considered in another octree, which also intersects the fluid element
          // add the new object ID to the vector
          IntersectingStructElemID.push_back( IDCurrentStructElem );
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

                  if( fabs(NormDiffVector) < epsilon )
                      return false;
              }
          }

          // if the new intersection node is not existing (as intersection with a corner point), then return false
          return true;
      }

      ///******************************************************************************************************************
      ///******************************************************************************************************************

      bool IsIntersectionOnCorner(IntersectionNodeStruct& NewIntersectionNode,
                                  Geometry< Node<3> >&    rFluidGeom )
      {
          array_1d<double,3> TetNode;
          array_1d<double,3> DiffVector;
          double NormDiffVector;
          //const double tolerance = 1e-8;

          for(unsigned int i_TetNode = 0 ; i_TetNode < 4 ; i_TetNode++)
          {
              TetNode  = rFluidGeom[i_TetNode].Coordinates();
              DiffVector[0] = TetNode[0] - NewIntersectionNode.Coordinates[0];
              DiffVector[1] = TetNode[1] - NewIntersectionNode.Coordinates[1];
              DiffVector[2] = TetNode[2] - NewIntersectionNode.Coordinates[2];
              NormDiffVector = norm_2(DiffVector);

              if( fabs(NormDiffVector) < epsilon )
                  return true;
          }

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

          Geometry< Node<3> >& rFluidGeom = i_fluid_element->GetGeometry();

          // Set Flag for intersected fluid elements
          i_fluid_element->GetValue(SPLIT_ELEMENT) = true;

          // Reduce all found intersection nodes located on each tetdrahedra edge to just one intersection node by averaging
          ComputeApproximationNodes(IntersectedTetEdges,NodesOfApproximatedStructure);

          // Intersection with one corner point
          if( (NodesOfApproximatedStructure.size() == 1) && (NumberIntersectionsOnTetCorner == 1) )
              CalcSignedDistancesToOneIntNode(rFluidGeom,NodesOfApproximatedStructure,ElementalDistances);

          // Intersection with two corner points / one tetrahedra edge
          else if( (NodesOfApproximatedStructure.size() == 2) && (NumberIntersectionsOnTetCorner == 2) )
              CalcSignedDistancesToTwoIntNodes(rFluidGeom,NodesOfApproximatedStructure,ElementalDistances);

          // Intersection with three tetrahedra edges
          else if( NodesOfApproximatedStructure.size() == 3 )
              CalcSignedDistancesToThreeIntNodes(i_fluid_element,NodesOfApproximatedStructure,IntersectedTetEdges,ElementalDistances);

          // Intersection with four tetrahedra edges
          else if( NodesOfApproximatedStructure.size() == 4 )
              CalcSignedDistancesToFourIntNodes(i_fluid_element,NodesOfApproximatedStructure,IntersectedTetEdges,ElementalDistances);

          // In case there is NO intersection with fluid element
          else
              i_fluid_element->GetValue(SPLIT_ELEMENT) = false;

          // Assign elemental distances
          if( i_fluid_element->GetValue(SPLIT_ELEMENT) == true )
              AssignDistancesToElements(i_fluid_element,ElementalDistances);

//          if(NodesOfApproximatedStructure.size()>2)
//          {
//              KRATOS_WATCH(NodesOfApproximatedStructure.size());
//              KRATOS_WATCH(ElementalDistances);
//              KRATOS_WATCH(i_fluid_element->GetValue(SPLIT_ELEMENT));
//              KRATOS_WATCH(i_fluid_element->GetValue(ELEMENTAL_DISTANCES));
//          }



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

              IntersectionNodeStruct NewAppxoimationNode;
              NewAppxoimationNode.Coordinates[0] = sum_X / NumberIntNodes;
              NewAppxoimationNode.Coordinates[1] = sum_Y / NumberIntNodes;
              NewAppxoimationNode.Coordinates[2] = sum_Z / NumberIntNodes;

              NodesOfApproximatedStructure.push_back(NewAppxoimationNode);
          }
      }

      ///******************************************************************************************************************
      ///******************************************************************************************************************

      void CalcSignedDistancesToOneIntNode( Geometry< Node<3> >&                rFluidGeom,
                                            std::vector<IntersectionNodeStruct> NodesOfApproximatedStructure,
                                            array_1d<double,4>&                 ElementalDistances)
      {
          const array_1d<double,3>& IntersectionNodeCoord = NodesOfApproximatedStructure[0].Coordinates;
          array_1d<double,3> DistVecTetNode;
          array_1d<double,3> TetNode;
          array_1d<double,3> NormalAtIntersectionNode;
          double             NormDistTetNode;
          double             InnerProduct;

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

      void CalcSignedDistancesToTwoIntNodes( Geometry< Node<3> >&                rFluidGeom,
                                             std::vector<IntersectionNodeStruct> NodesOfApproximatedStructure,
                                             array_1d<double,4>&                 ElementalDistances)
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

          const Point<3> LinePoint1 = Point<3>(IntersectionNode1Coord[0] , IntersectionNode1Coord[1] , IntersectionNode1Coord[2]);
          const Point<3> LinePoint2 = Point<3>(IntersectionNode2Coord[0] , IntersectionNode2Coord[1] , IntersectionNode2Coord[2]);

          for(unsigned int i_TetNode = 0 ; i_TetNode < 4 ; i_TetNode++)
          {
              // Get coordinates of the fluid element nodes
              TetNode  = rFluidGeom(i_TetNode)->Coordinates();

              // Compute distance to point
              NormDistTetNode = GeometryUtils::PointDistanceToLineSegment3D(LinePoint1, LinePoint2 , Point<3>(TetNode[0],TetNode[1],TetNode[2]));

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
          array_1d<unsigned int,3> IndexNodes;

          // Having 4 Intersection nodes, one can build 4 triangles. The 4 possible combinations are defined by indices in a matrix.
          bounded_matrix<unsigned int,4,3> IndexNodesTriangles;

          IndexNodesTriangles(0,0) = 0;
          IndexNodesTriangles(0,1) = 1;
          IndexNodesTriangles(0,2) = 2;

          IndexNodesTriangles(1,0) = 1;
          IndexNodesTriangles(1,1) = 2;
          IndexNodesTriangles(1,2) = 3;

          IndexNodesTriangles(2,0) = 0;
          IndexNodesTriangles(2,1) = 1;
          IndexNodesTriangles(2,2) = 3;

          IndexNodesTriangles(3,0) = 0;
          IndexNodesTriangles(3,1) = 2;
          IndexNodesTriangles(3,2) = 3;

          for(unsigned int i_Triangle = 0 ; i_Triangle < 4 ; i_Triangle++)
          {
              IndexNodes[0] = IndexNodesTriangles(i_Triangle,0);
              IndexNodes[1] = IndexNodesTriangles(i_Triangle,1);
              IndexNodes[2] = IndexNodesTriangles(i_Triangle,2);

              // Compute distance for each tetrahedra node to the triangle to approximate the structure
              CalcSignedDistancesToApproxTriangle( i_fluid_element , NodesOfApproximatedStructure , IntersectedTetEdges , ElementalDistances , IndexNodes );
          }
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
              Point<3>           ApproxTrianglePoint1;
              Point<3>           ApproxTrianglePoint2;
              Point<3>           ApproxTrianglePoint3;
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

              ApproxTrianglePoint1 = Point<3>(IntersectionNode1Coord[0] , IntersectionNode1Coord[1] , IntersectionNode1Coord[2]);
              ApproxTrianglePoint2 = Point<3>(IntersectionNode2Coord[0] , IntersectionNode2Coord[1] , IntersectionNode2Coord[2]);
              ApproxTrianglePoint3 = Point<3>(IntersectionNode3Coord[0] , IntersectionNode3Coord[1] , IntersectionNode3Coord[2]);

              // Compute distance from tet node to current triangle
              UnsignedDistance = GeometryUtils::PointDistanceToTriangle3D(ApproxTrianglePoint1, ApproxTrianglePoint2 , ApproxTrianglePoint3 , Point<3>(TetNode[0],TetNode[1],TetNode[2]));

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

      void GenerateOctree()
      {
          for(ModelPart::NodeIterator i_node = mrSkinModelPart.NodesBegin() ; i_node != mrSkinModelPart.NodesEnd() ; i_node++)
          {
              double temp_point[3];
              temp_point[0] = i_node->Coordinate(1);
              temp_point[1] = i_node->Coordinate(2);
              temp_point[2] = i_node->Coordinate(3);
              mOctree.Insert(temp_point);
          }

          mOctree.Constrain2To1(); // To be removed. Pooyan.

          //mOctree.RefineWithUniformSize(1.0/(std::pow(2.0,8)));

          for(ModelPart::ElementIterator i_element = mrSkinModelPart.ElementsBegin() ; i_element != mrSkinModelPart.ElementsEnd() ; i_element++)
          {
              mOctree.Insert(*(i_element).base());
          }
          KRATOS_WATCH(mOctree);

          std::cout << "######## WRITING OCTREE MESH #########" << std::endl;
          std::ofstream myfile;
          myfile.open ("unserbaum.post.msh");
          mOctree.PrintGiDMesh(myfile);
          myfile.close();
      }

      ///******************************************************************************************************************
      ///******************************************************************************************************************

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

      ///******************************************************************************************************************
      ///******************************************************************************************************************

      ///@}
      ///@name Access
      ///@{
      
      
      ///@}
      ///@name Inquiry
      ///@{
      
      
      ///@}
      ///@name Input and output
      ///@{

            
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

        DistanceSpatialContainersConfigure::data_type mOctreeNodes;

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

  const double CalculateSignedDistanceTo3DSkinProcess::epsilon = 1e-12;
  
  
}  // namespace Kratos.

#endif // KRATOS_CALCULATE_DISTANCE_PROCESS_H_INCLUDED  defined 


