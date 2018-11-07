//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Navaneeth K Narayanan
//
//

#if !defined(KRATOS_CALCULATE_SIGNED_DISTANCE_CONDITION_2D_PROCESS_H_INCLUDED)
#define KRATOS_CALCULATE_SIGNED_DISTANCE_CONDITION_2D_PROCESS_H_INCLUDED

// System includes
#include <string>
#include <iostream>
#include <ctime>

// External includes

// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "includes/model_part.h"
#include "includes/deprecated_variables.h"

#include "custom_utilities/quadtree_binary.h"
#include "utilities/spatial_containers_configure.h"
#include "utilities/timer.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h"
#include "geometries/triangle_3d_3.h"
#include "geometries/quadrilateral_3d_4.h"
#include "utilities/body_normal_calculation_utils.h"
#include "includes/kratos_flags.h"
#include "utilities/binbased_fast_point_locator.h"
#include "utilities/binbased_nodes_in_element_locator.h"
#include "elements/distance_calculation_element_simplex.h"
//#include "utilities/parallel_levelset_distance_calculator.h"

#ifdef _OPENMP
#include "omp.h"
#endif

using namespace boost::numeric::ublas;

namespace Kratos
{

class DistanceSpatialContainersConditionConfigure2d
{
  public:
    class CellNodeData
    {
        double mDistance;
        double mCoordinates[3];
        std::size_t mId;

      public:
        double &Distance() { return mDistance; }
        double &X() { return mCoordinates[0]; }
        double &Y() { return mCoordinates[1]; }
        double &Z() { return mCoordinates[2]; }
        double &Coordinate(int i) { return mCoordinates[i - 1]; }
        std::size_t &Id() { return mId; }
    };

    ///@name Type Definitions
    ///@{

    enum
    {
        Dimension = 3,
        DIMENSION = 2,
        MAX_LEVEL = 12,
        MIN_LEVEL = 2 // this cannot be less than 2!!!
    };

    typedef Point PointType; /// always the point 3D
    typedef std::vector<double>::iterator DistanceIteratorType;
    typedef PointerVectorSet<
                GeometricalObject::Pointer, 
                IndexedObject,
                std::less<typename IndexedObject::result_type>,
                std::equal_to<typename IndexedObject::result_type>,
                Kratos::shared_ptr<typename GeometricalObject::Pointer>,
                std::vector< Kratos::shared_ptr<typename GeometricalObject::Pointer> >
                    >  ContainerType;
    typedef ContainerType::value_type PointerType;
    typedef ContainerType::iterator IteratorType;
    typedef PointerVectorSet<
                GeometricalObject::Pointer, 
                IndexedObject,
                std::less<typename IndexedObject::result_type>,
                std::equal_to<typename IndexedObject::result_type>,
                Kratos::shared_ptr<typename GeometricalObject::Pointer>,
                std::vector< Kratos::shared_ptr<typename GeometricalObject::Pointer> >
                    >  ResultContainerType;
    typedef ResultContainerType::value_type ResultPointerType;
    typedef ResultContainerType::iterator ResultIteratorType;

    typedef GeometricalObject::Pointer pointer_type; //nav
    typedef CellNodeData cell_node_data_type;
    typedef std::vector<CellNodeData *> data_type;

    typedef std::vector<PointerType>::iterator PointerTypeIterator;

    /// Pointer definition of DistanceSpatialContainersConditionConfigure2d
    KRATOS_CLASS_POINTER_DEFINITION(DistanceSpatialContainersConditionConfigure2d);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    DistanceSpatialContainersConditionConfigure2d() {}

    /// Destructor.
    virtual ~DistanceSpatialContainersConditionConfigure2d() {}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    static data_type *AllocateData()
    {
        return new data_type(8, (CellNodeData *)NULL);
    }

    static void CopyData(data_type *source, data_type *destination)
    {
        *destination = *source;
    }

    static void DeleteData(data_type *data)
    {
        delete data;
    }

    static inline void CalculateBoundingBox(const PointerType &rObject, PointType &rLowPoint, PointType &rHighPoint)
    {
        rHighPoint = rObject->GetGeometry().GetPoint(0);
        rLowPoint = rObject->GetGeometry().GetPoint(0);

        for (std::size_t point = 0; point < rObject->GetGeometry().PointsNumber(); point++)
        {
            for (std::size_t i = 0; i < 3; i++)
            {
                rLowPoint[i] = (rLowPoint[i] > rObject->GetGeometry().GetPoint(point)[i]) ? rObject->GetGeometry().GetPoint(point)[i] : rLowPoint[i];
                rHighPoint[i] = (rHighPoint[i] < rObject->GetGeometry().GetPoint(point)[i]) ? rObject->GetGeometry().GetPoint(point)[i] : rHighPoint[i];
            }
        }
    }

    static inline void GetBoundingBox(const PointerType rObject, double *rLowPoint, double *rHighPoint)
    {

        for (std::size_t i = 0; i < 3; i++)
        {
            rLowPoint[i] = rObject->GetGeometry().GetPoint(0)[i];
            rHighPoint[i] = rObject->GetGeometry().GetPoint(0)[i];
        }

        for (std::size_t point = 0; point < rObject->GetGeometry().PointsNumber(); point++)
        {
            for (std::size_t i = 0; i < 3; i++)
            {
                rLowPoint[i] = (rLowPoint[i] > rObject->GetGeometry().GetPoint(point)[i]) ? rObject->GetGeometry().GetPoint(point)[i] : rLowPoint[i];
                rHighPoint[i] = (rHighPoint[i] < rObject->GetGeometry().GetPoint(point)[i]) ? rObject->GetGeometry().GetPoint(point)[i] : rHighPoint[i];
            }
        }
    }

    static inline bool Intersection(const PointerType &rObj_1, const PointerType &rObj_2)
    {
        Element::GeometryType &geom_1 = rObj_1->GetGeometry();
        Element::GeometryType &geom_2 = rObj_2->GetGeometry();
        return geom_1.HasIntersection(geom_2);
    }

    static inline bool IntersectionBox(const PointerType &rObject, const PointType &rLowPoint, const PointType &rHighPoint)
    {
        return rObject->GetGeometry().HasIntersection(rLowPoint, rHighPoint);
    }

    static inline bool IsIntersected(const Element::Pointer rObject, double Tolerance, const double *rLowPoint, const double *rHighPoint)
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
    virtual void PrintInfo(std::ostream &rOStream) const {}

    /// Print object's data.
    virtual void PrintData(std::ostream &rOStream) const {}

    ///@}

  protected:
  private:
    /// Assignment operator.
    DistanceSpatialContainersConditionConfigure2d &operator=(DistanceSpatialContainersConditionConfigure2d const &rOther);

    /// Copy constructor.
    DistanceSpatialContainersConditionConfigure2d(DistanceSpatialContainersConditionConfigure2d const &rOther);

}; // Class DistanceSpatialContainersConditionConfigure2d

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
class CalculateSignedDistanceTo2DConditionSkinProcess
    : public Process
{
  public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of CalculateSignedDistanceTo2DConditionSkinProcess
    KRATOS_CLASS_POINTER_DEFINITION(CalculateSignedDistanceTo2DConditionSkinProcess);

    typedef DistanceSpatialContainersConditionConfigure2d ConfigurationType;
    typedef QuadtreeBinaryCell<ConfigurationType> CellType;
    typedef QuadtreeBinary<CellType> QuadtreeType;
    typedef ConfigurationType::cell_node_data_type CellNodeDataType;
    typedef Point PointType; /// always the point 3D
    typedef QuadtreeType::cell_type::object_container_type object_container_type;
    typedef struct
    {
        array_1d<double, 3> Coordinates;
        array_1d<double, 3> StructElemNormal;
        std::size_t EdgeNode1;
        std::size_t EdgeNode2;
    } IntersectionNodeStruct;
    typedef struct
    {
        std::vector<IntersectionNodeStruct> IntNodes;
    } TriangleEdgeStruct;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    CalculateSignedDistanceTo2DConditionSkinProcess(ModelPart &rThisModelPartStruc, ModelPart &rThisModelPartFluid)
        : mrSkinModelPart(rThisModelPartStruc), mrFluidModelPart(rThisModelPartFluid)
    {
        //this->pDistanceCalculator = typename ParallelDistanceCalculator<2>::Pointer(new ParallelDistanceCalculator<2>());
        dist_limit = 1e-10;
    }

    /// Destructor.
    virtual ~CalculateSignedDistanceTo2DConditionSkinProcess()
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

    virtual void Execute() override
    {
        KRATOS_TRY;

        //GenerateFluidModelPartbasedOnBoundingBox();
        GenerateQuadtree();

        DistanceFluidStructure();
        //std::size_t max_level = 100;
        //double max_distance = 200;

        //          ------------------------------------------------------------------
        //          GenerateNodes();
        CalculateDistance2(); // I have to change this. Pooyan.
        //mrSkinModelPart.GetCommunicator().AssembleCurrentData(DISTANCE);
        //          std::ofstream mesh_file1("quadtree1.post.msh");
        //          std::ofstream res_file("quadtree1.post.res");
        //          Timer::Start("Writing Gid conform Mesh");
        //          PrintGiDMesh(mesh_file1);
        //          PrintGiDResults(res_file);
        //          quadtree.PrintGiDMeshNew(mesh_file2);
        //          Timer::Stop("Writing Gid conform Mesh");
        //          delete quadtree. TODO: Carlos
        //          ------------------------------------------------------------------
        //this->pDistanceCalculator->CalculateDistances(mrFluidModelPart,DISTANCE,NODAL_AREA,max_level,max_distance);
        AvoidZeroNodalDistances();
        KRATOS_CATCH("");
    }

    ///******************************************************************************************************************

    ///******************************************************************************************************************

    // generates
    /*void GenerateFluidModelPartbasedOnBoundingBox()
    {

        array_1d<double,3> low;
        array_1d<double,3> high;
        array_1d<double,3> size;


        // loop over all skin nodes
        for(ModelPart::NodeIterator i_node = mrSkinModelPart.NodesBegin();
            i_node != mrSkinModelPart.NodesEnd();
            i_node++)
        {
            for (int i = 0 ; i < 3; i++)
            {
                low[i]  = i_node->Coordinate(i+1) < low[i]  ? i_node->Coordinate(i+1) : low[i];
                high[i] = i_node->Coordinate(i+1) > high[i] ? i_node->Coordinate(i+1) : high[i];
            }
        }*/

    /*for (int i = 0 ; i < 3; i++)
            {
                //considering 10% more length
                size[i] = high[i]-low[i];
                low[i]  = low[i]- 0.1*size[i];
                high[i] = high[i]*0.1*size[i];
            }
   */

    /*std::cout<<"Lowest Dimension"<<low[0]<<","<<low[1]<<","<<low[2]<<","<<std::endl;

        std::cout<<"Highest Dimension"<<high[0]<<","<<high[1]<<","<<high[2]<<","<<std::endl;

        for(ModelPart::ElementIterator elem = mrBackgroundModelPart.ElementsBegin(); elem != mrBackgroundModelPart.ElementsEnd(); elem ++)

        {

            Geometry<Node<3> >& geom = elem->GetGeometry();
            array_1d<double,3> coords;
            std::size_t numOfPointsInside = 0;
			for (int j = 0 ; j < geom.size(); j++){
				coords = elem->GetGeometry()[j].Coordinates();

				{
				  bool isInside = IsInside(coords,high,low);
                  numOfPointsInside++;

				}
			}

            if (numOfPointsInside > 0)
            {

            Element::Pointer pElem = Element::Pointer(new Element(*elem));
			mrFluidModelPart.Elements().push_back(pElem);

            }



        }

        std::cout<<"Fluid Model Part that is considered"<<mrFluidModelPart.Elements().size()<<std::endl;


    }    */

    //Checks if the node is inside the bounding box implemented only for 2d case
    /*  bool IsInside (array_1d<double,3> coords, array_1d<double,3> high,array_1d<double,3> low)

        {

            array_1d<double, 3 > ab;
            array_1d<double, 3 > ad;
            array_1d<double, 3 > am;

            ab[0] = high[0]-low[0];
            ad[1] = high[1]-low[1];
            ab[2] = 0.0;
            ad[2] = 0.0;
            for (int i = 0; i <3 ; i++ )
            am[i] = coords[i]-low[i];

            if (inner_prod(ab,am) < 0 && inner_prod(ab,am)>inner_prod(ab,ab))
            return false;

            if (inner_prod(ad,am) < 0 && inner_prod(ad,am)>inner_prod(ad,ad))
            return false;

            return true;


        }

*/

    void DistanceFluidStructure()
    {
        //std::cout << "Start calculating Elemental distances..." << std::endl;

        // Initialize nodal distances in the domain
        InitializeDistances();

        // Initialize index table that defines line Edges of fluid Element
        bounded_matrix<std::size_t, 3, 2> TriangleEdgeIndexTable;
        SetIndexTable(TriangleEdgeIndexTable);

// loop over all fluid Elements
// this loop is parallelized using openmp
#ifdef _OPENMP
        int number_of_threads = omp_get_max_threads();
#else
        int number_of_threads = 1;
#endif

        ModelPart::ElementsContainerType &pElements = mrFluidModelPart.Elements();

        vector<std::size_t> Element_partition;
        CreatePartition(number_of_threads, pElements.size(), Element_partition);

#pragma omp parallel for
        for (int k = 0; k < number_of_threads; k++)
        {
            ModelPart::ElementsContainerType::iterator it_begin = pElements.ptr_begin() + Element_partition[k];
            ModelPart::ElementsContainerType::iterator it_end = pElements.ptr_begin() + Element_partition[k + 1];

            // assemble all Elements
            for (ModelPart::ElementIterator it = it_begin; it != it_end; ++it)
            {
                CalcNodalDistancesOfTriangleNodes(it, TriangleEdgeIndexTable);
            }
        }

        // Finally, each triangle Element has 3 distance values. But each node belongs to
        // several Elements, such that it is assigned several distance values
        // --> now synchronize these values by finding the minimal distance and assign to each node a minimal nodal distance
        AssignMinimalNodalDistance(); // revisit -nav
    }


    void InitializeDistances()
    {
        const double initial_distance = 10000;

        ModelPart::NodesContainerType::ContainerType &nodes = mrFluidModelPart.NodesArray();

        // reset the node distance to 1.0 which is the maximum distance in our normalized space.
        int nodesSize = nodes.size();

#pragma omp parallel for firstprivate(nodesSize)
        for (int i = 0; i < nodesSize; i++)
            nodes[i]->GetSolutionStepValue(DISTANCE) = initial_distance;

        ModelPart::ElementsContainerType::ContainerType &fluid_Elements = mrFluidModelPart.ElementsArray();

        array_1d<double, 3> ElementalDistances;
        ElementalDistances[0] = initial_distance;
        ElementalDistances[1] = initial_distance;
        ElementalDistances[2] = initial_distance;

        // reset the Elemental distance to 1.0 which is the maximum distance in our normalized space.
        // also initialize the embedded velocity of the fluid Element
        int ElementsSize = fluid_Elements.size();

        //#pragma omp parallel for firstprivate(ElementsSize) no omp because it allocates the memory!
        for (int i = 0; i < ElementsSize; i++)
        {
            fluid_Elements[i]->GetValue(ELEMENTAL_DISTANCES) = ElementalDistances;
            //fluid_Elements[i]->GetValue(SPLIT_ELEMENT) = false;
            //fluid_Elements[i]->GetValue(EMBEDDED_VELOCITY)=ZeroVector(3);
        }
    }

    ///******************************************************************************************************************
    ///******************************************************************************************************************

    void SetIndexTable(bounded_matrix<std::size_t, 3, 2> &TriangleEdgeIndexTable)
    {
        // Initialize index table to define line Edges of fluid Element
        TriangleEdgeIndexTable(0, 0) = 0;
        TriangleEdgeIndexTable(0, 1) = 1;
        TriangleEdgeIndexTable(1, 0) = 0;
        TriangleEdgeIndexTable(1, 1) = 2;
        TriangleEdgeIndexTable(2, 0) = 1;
        TriangleEdgeIndexTable(2, 1) = 2;
    }

    ///******************************************************************************************************************
    ///******************************************************************************************************************

    void CalcNodalDistancesOfTriangleNodes(ModelPart::ElementsContainerType::iterator &i_fluidElement,
                                           bounded_matrix<std::size_t, 3, 2> TriangleEdgeIndexTable)
    {
        std::vector<QuadtreeType::cell_type *> leaves;
        std::vector<TriangleEdgeStruct> IntersectedTriangleEdges;
        std::size_t NumberIntersectionsOnTriangleCorner = 0;

        // Get leaves of quadtree intersecting with fluid Element
        mpQuadtree->GetIntersectedLeaves(*(i_fluidElement).base(), leaves);

        //std::cout<<"Fluid element"<<i_fluidElement->Id()<<std::endl;

        int intersection_counter = 0;

        // Loop over all 6 line Edges of the triangle
        for (std::size_t i_triangleEdge = 0;
             i_triangleEdge < 3;
             i_triangleEdge++)
        {
            IdentifyIntersectionNodes(i_fluidElement, i_triangleEdge, leaves, IntersectedTriangleEdges, NumberIntersectionsOnTriangleCorner, TriangleEdgeIndexTable, intersection_counter);
        }

        /*if (intersection_counter!=0)
        {
            i_fluidElement->GetValue(EMBEDDED_VELOCITY) /= intersection_counter;
        }*/

        if (IntersectedTriangleEdges.size() > 0)

            CalcDistanceTo2DSkin(IntersectedTriangleEdges, i_fluidElement, NumberIntersectionsOnTriangleCorner);
    }

    ///******************************************************************************************************************
    ///******************************************************************************************************************

    void IdentifyIntersectionNodes(ModelPart::ElementsContainerType::iterator &i_fluidElement,
                                   std::size_t i_triangleEdge,
                                   std::vector<QuadtreeType::cell_type *> &leaves,
                                   std::vector<TriangleEdgeStruct> &IntersectedTriangleEdges,
                                   std::size_t &NumberIntersectionsOnTriangleCorner,
                                   bounded_matrix<std::size_t, 3, 2> TriangleEdgeIndexTable,
                                   int &intersection_counter)

    {
        std::vector<std::size_t> IntersectingStructCondID;

        TriangleEdgeStruct NewTriangleEdge;
        std::size_t NumberIntersectionsOnTriangleCornerCurrentEdge = 0;

        // Get nodes of line Edge
        std::size_t EdgeStartIndex = TriangleEdgeIndexTable(i_triangleEdge, 0);
        std::size_t EdgeEndIndex = TriangleEdgeIndexTable(i_triangleEdge, 1);

        PointType &P1 = i_fluidElement->GetGeometry()[EdgeStartIndex];
        PointType &P2 = i_fluidElement->GetGeometry()[EdgeEndIndex];

        double EdgeNode1[3] = {P1.X(), P1.Y(), P1.Z()};
        double EdgeNode2[3] = {P2.X(), P2.Y(), P2.Z()};

        // loop over all quadtree cells which are intersected by the fluid Element
        for (std::size_t i_cell = 0; i_cell < leaves.size(); i_cell++)
        {
            // Structural Element contained in one cell of the quadtree

            object_container_type *struct_cond = (leaves[i_cell]->pGetObjects());

            // loop over all structural Conditions within each quadtree cell
            for (object_container_type::iterator i_StructCondition = struct_cond->begin(); i_StructCondition != struct_cond->end(); i_StructCondition++)
            {

                if (StructuralElementNotYetConsidered((*i_StructCondition)->Id(), IntersectingStructCondID))
                {

                    // Calculate and associate intersection point to the current fluid Element
                    double IntersectionPoint[3] = {0.0, 0.0, 0.0};

                    int TriangleEdgeHasIntersections = IntersectionLineSegment((*i_StructCondition)->GetGeometry(), EdgeNode1, EdgeNode2, IntersectionPoint);

                    if (TriangleEdgeHasIntersections == 1)
                    {

                        IntersectionNodeStruct NewIntersectionNode;

                        // Assign information to the intersection node
                        NewIntersectionNode.Coordinates[0] = IntersectionPoint[0];
                        NewIntersectionNode.Coordinates[1] = IntersectionPoint[1];
                        NewIntersectionNode.Coordinates[2] = IntersectionPoint[2];

                        if (IsIntersectionNodeOnTriangleEdge(IntersectionPoint, EdgeNode1, EdgeNode2))
                        {

                            if (IsNewIntersectionNode(NewIntersectionNode, IntersectedTriangleEdges))
                            {

                                CalculateNormal2D((*i_StructCondition)->GetGeometry()[0],
                                                  (*i_StructCondition)->GetGeometry()[1],
                                                  NewIntersectionNode.StructElemNormal);

                                // check, how many intersection nodes are located on corner points of the triangle
                                if (IsIntersectionOnCorner(NewIntersectionNode, EdgeNode1, EdgeNode2))
                                {
                                    NumberIntersectionsOnTriangleCornerCurrentEdge++;

                                    //std::cout << " ####Intersection on corner####" << std::endl;

                                    // only allow one intersection node on a triangle edge
                                    if (NumberIntersectionsOnTriangleCornerCurrentEdge < 2)
                                    {
                                        // add the new intersection point to the list of intersection points of the fluid Element
                                        NewIntersectionNode.EdgeNode1 = EdgeStartIndex;
                                        NewIntersectionNode.EdgeNode2 = EdgeEndIndex;
                                        NewTriangleEdge.IntNodes.push_back(NewIntersectionNode);

                                        // if triangle edge belonging to this intersection point is not already marked as "IntersectedTriangleEdge" --> put it into the respective container
                                        // when a second intersection node is found, then it is not necessary to push_back again
                                        if (NewTriangleEdge.IntNodes.size() == 1)
                                            IntersectedTriangleEdges.push_back(NewTriangleEdge);
                                    }

                                    // this corner intersection node is only considered once for each triangle edge
                                    if (NumberIntersectionsOnTriangleCornerCurrentEdge == 1)
                                    {
                                        NumberIntersectionsOnTriangleCorner++;
                                    }
                                }
                                else
                                {
                                    // add the new intersection point to the list of intersection points of the fluid Element
                                    NewIntersectionNode.EdgeNode1 = EdgeStartIndex;
                                    NewIntersectionNode.EdgeNode2 = EdgeEndIndex;
                                    NewTriangleEdge.IntNodes.push_back(NewIntersectionNode);

                                    // velocity mapping structure --> fluid
                                    /*                                    array_1d<double,3> emb_vel = (*i_StructCondition)->GetGeometry()[0].GetSolutionStepValue(VELOCITY);
                                    emb_vel += (*i_StructCondition)->GetGeometry()[1].GetSolutionStepValue(VELOCITY);
                                    emb_vel += (*i_StructCondition)->GetGeometry()[2].GetSolutionStepValue(VELOCITY);*/
                                    /*
                                    i_fluidElement->GetValue(EMBEDDED_VELOCITY) += emb_vel/3;*/
                                    intersection_counter++;
                                }
                            }
                        }
                    }
                }
            }
        }

        // Finally put the found intersection nodes into the container
        if (NewTriangleEdge.IntNodes.size() > 0)
        {
            if (NumberIntersectionsOnTriangleCornerCurrentEdge == 0)
                IntersectedTriangleEdges.push_back(NewTriangleEdge);
        }
    }

    ///******************************************************************************************************************
    ///******************************************************************************************************************

    bool StructuralElementNotYetConsidered(std::size_t IDCurrentStructCond,
                                           std::vector<std::size_t> &IntersectingStructCondID)
    {
        // check if the structural Element was already considered as intersecting Element
        for (std::size_t k = 0; k < IntersectingStructCondID.size(); k++)
        {
            if (IDCurrentStructCond == IntersectingStructCondID[k])
                return false;
        }

        // if structural Element has not been considered in another quadtree, which also intersects the fluid Element
        // add the new object ID to the vector
        IntersectingStructCondID.push_back(IDCurrentStructCond);
        return true;
    }

    ///******************************************************************************************************************
    ///******************************************************************************************************************

    bool IsIntersectionNodeOnTriangleEdge(double *IntersectionPoint,
                                          double *EdgeNode1,
                                          double *EdgeNode2)
    {
        // check, if intersection point is located on any edge of the fluid Element
        array_1d<double, 3> ConnectVectTriangleNodeIntNode1;
        array_1d<double, 3> ConnectVectTriangleNodeIntNode2;
        array_1d<double, 3> EdgeVector;

        ConnectVectTriangleNodeIntNode1[0] = IntersectionPoint[0] - EdgeNode1[0];
        ConnectVectTriangleNodeIntNode1[1] = IntersectionPoint[1] - EdgeNode1[1];
        ConnectVectTriangleNodeIntNode1[2] = IntersectionPoint[2] - EdgeNode1[2];

        ConnectVectTriangleNodeIntNode2[0] = IntersectionPoint[0] - EdgeNode2[0];
        ConnectVectTriangleNodeIntNode2[1] = IntersectionPoint[1] - EdgeNode2[1];
        ConnectVectTriangleNodeIntNode2[2] = IntersectionPoint[2] - EdgeNode2[2];

        double LengthConnectVect1 = norm_2(ConnectVectTriangleNodeIntNode1);
        double LengthConnectVect2 = norm_2(ConnectVectTriangleNodeIntNode2);

        EdgeVector[0] = EdgeNode2[0] - EdgeNode1[0];
        EdgeVector[1] = EdgeNode2[1] - EdgeNode1[1];
        EdgeVector[2] = EdgeNode2[2] - EdgeNode1[2];

        double MaxEdgeLength = norm_2(EdgeVector);

        // if both connection vectors (corner point --> intersection point)
        // are smaller or equal to the edge length of triangle,
        // then intersection point is located on the edge
        if ((LengthConnectVect1 <= (MaxEdgeLength)) && (LengthConnectVect2 <= (MaxEdgeLength)))
            return true;
        else
            return false;
    }

    ///******************************************************************************************************************
    ///******************************************************************************************************************

    bool IsNewIntersectionNode(IntersectionNodeStruct &NewIntersectionNode,
                               std::vector<TriangleEdgeStruct> &IntersectedTriangleEdges)
    {
        array_1d<double, 3> DiffVector;
        double NormDiffVector = 0;
        std::size_t NumberIntNodes = 0;

        for (std::size_t i_TriangleEdge = 0; i_TriangleEdge < IntersectedTriangleEdges.size(); i_TriangleEdge++)
        {
            NumberIntNodes = IntersectedTriangleEdges[i_TriangleEdge].IntNodes.size();
            for (std::size_t i_IntNode = 0; i_IntNode < NumberIntNodes; i_IntNode++)
            {
                DiffVector[0] = NewIntersectionNode.Coordinates[0] - IntersectedTriangleEdges[i_TriangleEdge].IntNodes[i_IntNode].Coordinates[0];
                DiffVector[1] = NewIntersectionNode.Coordinates[1] - IntersectedTriangleEdges[i_TriangleEdge].IntNodes[i_IntNode].Coordinates[1];
                DiffVector[2] = NewIntersectionNode.Coordinates[2] - IntersectedTriangleEdges[i_TriangleEdge].IntNodes[i_IntNode].Coordinates[2];

                NormDiffVector = norm_2(DiffVector);

                if (NormDiffVector < epsilon)
                    return false;
            }
        }

        // if the new intersection node is not existing (as intersection with a corner point), then return false
        return true;
    }

    ///******************************************************************************************************************
    ///******************************************************************************************************************

    bool IsIntersectionOnCorner(IntersectionNodeStruct &NewIntersectionNode,
                                double *EdgeNode1,
                                double *EdgeNode2)
    {
        array_1d<double, 3> DiffVector;
        double NormDiffVector;

        DiffVector[0] = EdgeNode1[0] - NewIntersectionNode.Coordinates[0];
        DiffVector[1] = EdgeNode1[1] - NewIntersectionNode.Coordinates[1];
        DiffVector[2] = EdgeNode1[2] - NewIntersectionNode.Coordinates[2];
        NormDiffVector = norm_2(DiffVector);

        if (NormDiffVector < epsilon)
            return true;

        DiffVector[0] = EdgeNode2[0] - NewIntersectionNode.Coordinates[0];
        DiffVector[1] = EdgeNode2[1] - NewIntersectionNode.Coordinates[1];
        DiffVector[2] = EdgeNode2[2] - NewIntersectionNode.Coordinates[2];
        NormDiffVector = norm_2(DiffVector);

        if (NormDiffVector < epsilon)
            return true;
        else
            return false;
    }

    ///******************************************************************************************************************
    ///******************************************************************************************************************

    void CalculateNormal2D(Point &Point1,
                           Point &Point2,
                           array_1d<double, 3> &rResultNormal)
    {
        rResultNormal[0] = (Point2[1] - Point1[1]);
        rResultNormal[1] = -(Point2[0] - Point1[0]);
        rResultNormal[2] = 0.00;
        //std::cout<<" Normal "<<rResultNormal[0]<<","<<rResultNormal[1]<<std::endl;
    }

    ///******************************************************************************************************************
    ///******************************************************************************************************************

    void CalcDistanceTo2DSkin(std::vector<TriangleEdgeStruct> &IntersectedTriangleEdges,
                              ModelPart::ElementsContainerType::iterator &i_fluid_Element,
                              std::size_t NumberIntersectionsOnTriangleCorner)
    {
        std::vector<IntersectionNodeStruct> NodesOfApproximatedStructure;
        array_1d<double, 3> ElementalDistances;

        FillIntNodesContainer(IntersectedTriangleEdges, NodesOfApproximatedStructure);

        // Intersection with one corner point
        if (NodesOfApproximatedStructure.size() == 1 && NumberIntersectionsOnTriangleCorner == 1)
        {
            CalcSignedDistancesToOneIntNode(i_fluid_Element, NodesOfApproximatedStructure, ElementalDistances);
            i_fluid_Element->GetValue(SPLIT_ELEMENT) = true;
            //std::cout<<"CalcSignedDistancesToOneIntNode"<<"ED "<<ElementalDistances[0]<<","<<ElementalDistances[1]<<","<<ElementalDistances[2]<<std::endl;
        }

        // Intersection with two corner points / one triangle edge
        /*if (NodesOfApproximatedStructure.size() == 2 && NumberIntersectionsOnTriangleCorner == 2)
        {
            CalcSignedDistancesToTwoEdgeIntNodes(i_fluid_Element, NodesOfApproximatedStructure, ElementalDistances); //changed
            i_fluid_Element->GetValue(SPLIT_ELEMENT) = true;
            //_DEBUG
        if (i_fluid_Element->Id() == 3504 || i_fluid_Element->Id() == 544)
            std::cout<<"CalcSignedDistancesToTwoIntNodes"<<"ED "<<ElementalDistances[0]<<","<<ElementalDistances[1]<<","<<ElementalDistances[2]<<std::endl;
        }*/

        // Intersection with two triangle edges
        if (NodesOfApproximatedStructure.size() == 2 && NumberIntersectionsOnTriangleCorner <= 2)
        {
            CalcSignedDistancesToTwoEdgeIntNodes(i_fluid_Element, NodesOfApproximatedStructure, ElementalDistances);
            i_fluid_Element->GetValue(SPLIT_ELEMENT) = true;

            //std::cout<<"CalcSignedDistancesToTwoEdgeIntNodes"<<"ED "<<ElementalDistances[0]<<","<<ElementalDistances[1]<<","<<ElementalDistances[2]<<std::endl;
            /*if(i_fluid_Element->Id()==2541)

            {
                std::cout<<"CalcSignedDistancesToTwoEdgeIntNodes"<<"ED "<<ElementalDistances[0]<<","<<ElementalDistances[1]<<","<<ElementalDistances[2]<<std::endl;

                std::cout<<"Normal1 "<<NodesOfApproximatedStructure[0].StructElemNormal<<"Normal2 "<<NodesOfApproximatedStructure[1].StructElemNormal<<std::endl;
                std::exit(-1);


            }*/
        }

        // Intersection with more than two triangle edges
        if (NodesOfApproximatedStructure.size() > 2)
        {
            CalcSignedDistancesToMoreThanTwoIntNodes(i_fluid_Element, NodesOfApproximatedStructure, ElementalDistances, IntersectedTriangleEdges);
            i_fluid_Element->GetValue(SPLIT_ELEMENT) = true;
        }

        // Postprocessing treatment of Elemental distances
        if (i_fluid_Element->GetValue(SPLIT_ELEMENT) == true)
            AvoidZeroDistances(i_fluid_Element, ElementalDistances);

        // In case there is intersection with fluid Element: assign distances to the Element
        if (i_fluid_Element->GetValue(SPLIT_ELEMENT) == true)
        {
            //std::cout<<"Assigning distances to element"<<std::endl;

            i_fluid_Element->GetValue(ELEMENTAL_DISTANCES) = ElementalDistances;

            //std::cout<<"Assignment done"<<std::endl;
        }
    }

    ///******************************************************************************************************************
    ///******************************************************************************************************************

    void FillIntNodesContainer(std::vector<TriangleEdgeStruct> &IntersectedTriangleEdges,
                               std::vector<IntersectionNodeStruct> &NodesOfApproximatedStructure)
    {
        const std::size_t NumberCutEdges = IntersectedTriangleEdges.size();

        for (std::size_t i_TriangleEdge = 0; i_TriangleEdge < NumberCutEdges; i_TriangleEdge++)
        {
            std::size_t NumberIntNodes = IntersectedTriangleEdges[i_TriangleEdge].IntNodes.size();

            for (std::size_t i_IntNode = 0; i_IntNode < NumberIntNodes; i_IntNode++)
            {
                NodesOfApproximatedStructure.push_back(IntersectedTriangleEdges[i_TriangleEdge].IntNodes[i_IntNode]);
            }
        }
    }

    ///******************************************************************************************************************
    ///******************************************************************************************************************

    void CalcSignedDistancesToOneIntNode(ModelPart::ElementsContainerType::iterator &i_fluid_Element,
                                         std::vector<IntersectionNodeStruct> NodesOfApproximatedStructure,
                                         array_1d<double, 3> &ElementalDistances)
    {
        Geometry<Node<3>> &rFluidGeom = i_fluid_Element->GetGeometry();

        Point P1;
        P1.Coordinates() = NodesOfApproximatedStructure[0].Coordinates;

        array_1d<double, 3> &Normal = NodesOfApproximatedStructure[0].StructElemNormal;

        // Compute distance values for all triangle-nodes
        for (std::size_t i_TriangleNode = 0; i_TriangleNode < 3; i_TriangleNode++)
        {
            ElementalDistances[i_TriangleNode] = PointDistanceToLine(P1, Normal, rFluidGeom[i_TriangleNode]);
        }
    }

    ///******************************************************************************************************************
    ///******************************************************************************************************************

    void CalcSignedDistancesToTwoIntNodes(ModelPart::ElementsContainerType::iterator &i_fluid_Element,
                                          std::vector<IntersectionNodeStruct> NodesOfApproximatedStructure,
                                          array_1d<double, 3> &ElementalDistances)
    {
        Geometry<Node<3>> &rFluidGeom = i_fluid_Element->GetGeometry();

        Point P1;
        P1.Coordinates() = NodesOfApproximatedStructure[0].Coordinates;

        // Get normal at intersections, average them and check direction of distances
        array_1d<double, 3> NormalAtIntersectionNode1 = NodesOfApproximatedStructure[0].StructElemNormal;
        array_1d<double, 3> NormalAtIntersectionNode2 = NodesOfApproximatedStructure[1].StructElemNormal;

        // Compute normal of surface plane
        array_1d<double, 3> Normal;
        Normal[0] = 0.5 * (NormalAtIntersectionNode1[0] + NormalAtIntersectionNode2[0]);
        Normal[1] = 0.5 * (NormalAtIntersectionNode1[1] + NormalAtIntersectionNode2[1]);
        Normal[2] = 0.5 * (NormalAtIntersectionNode1[2] + NormalAtIntersectionNode2[2]);

        // Check whether orientation of normal is in direction of the normal of the intersecting structure
        // Note: The normal of the approx. surface can be max. 90deg to every surrounding normal of the structure at the intersection nodes
        const array_1d<double, 3> NormalAtOneIntersectionNode = NodesOfApproximatedStructure[0].StructElemNormal;

        bool NormalWrongOriented = false;
        if (inner_prod(NormalAtOneIntersectionNode, Normal) < 0)
            NormalWrongOriented = true;

        // switch direction of normal
        if (NormalWrongOriented)
            Normal *= -1;

        // Compute distance values for all triangle-nodes
        for (std::size_t i_TriangleNode = 0; i_TriangleNode < 3; i_TriangleNode++)
        {
            ElementalDistances[i_TriangleNode] = PointDistanceToLine(P1, Normal, rFluidGeom[i_TriangleNode]);
        }
    }

    ///******************************************************************************************************************
    ///******************************************************************************************************************

    void CalcSignedDistancesToTwoEdgeIntNodes(ModelPart::ElementsContainerType::iterator &i_fluid_Element,
                                              std::vector<IntersectionNodeStruct> &NodesOfApproximatedStructure,
                                              array_1d<double, 3> &ElementalDistances)
    {
        Geometry<Node<3>> &rFluidGeom = i_fluid_Element->GetGeometry();

        Point P1;
        Point P2;

        P1.Coordinates() = NodesOfApproximatedStructure[0].Coordinates;
        P2.Coordinates() = NodesOfApproximatedStructure[1].Coordinates;

        array_1d<double, 3> Normal;
        CalculateNormal2D(P1, P2, Normal);
        //nav
        //std::cout<<" Calculated  Normal " << Normal<<std::endl;
        // Check whether orientation of normal is in direction of the normal of the intersecting structure
        // Note: The normal of the approx. surface can be max. 90deg to every surrounding normal of the structure at the intersection nodes
        const array_1d<double, 3> NormalAtOneIntersectionNode = NodesOfApproximatedStructure[0].StructElemNormal;
        //nav
        //std::cout<<" Intersection Normal " << NodesOfApproximatedStructure[0].StructElemNormal<<std::endl;
        //std::cout<<" P1 " <<P1.Coordinates()<<std::endl;

        bool NormalWrongOriented = false;
        if (inner_prod(NormalAtOneIntersectionNode, Normal) < 0)
            NormalWrongOriented = true;

        // switch direction of normal
        if (NormalWrongOriented)
            Normal *= -1;

        // Compute distance values for all triangle-nodes
        for (std::size_t i_TriangleNode = 0; i_TriangleNode < 3; i_TriangleNode++)
        {
            ElementalDistances[i_TriangleNode] = PointDistanceToLine(P1, Normal, rFluidGeom[i_TriangleNode]);
        }
    }

    ///******************************************************************************************************************
    ///******************************************************************************************************************

    void CalcSignedDistancesToMoreThanTwoIntNodes(ModelPart::ElementsContainerType::iterator &i_fluid_Element,
                                                  std::vector<IntersectionNodeStruct> NodesOfApproximatedStructure,
                                                  array_1d<double, 3> &ElementalDistances,
                                                  std::vector<TriangleEdgeStruct> &IntersectedTriangleEdges)
    {
        std::size_t numberCutEdges = NodesOfApproximatedStructure.size();

        // Compute average of the intersection nodes which is a node on the line we look for
        Point P_mean;
        for (std::size_t k = 0; k < numberCutEdges; k++)
            for (std::size_t i = 0; i < 2; i++)
                P_mean.Coordinates()[i] += NodesOfApproximatedStructure[k].Coordinates[i];

        P_mean.Coordinates()[2] = 0.0;
        /*
        std::cout<<"Element Id "<<i_fluid_Element->Id()<<std::endl;
        std::cout<<"Number of cut edges "<<NodesOfApproximatedStructure.size()<<std::endl;*/
        //nav
        /*for(std::size_t k=0; k<numberCutEdges; k++)
            {

                std::cout<<"Coordinate "<< k <<" "<< NodesOfApproximatedStructure[k].Coordinates[0]<<","<<NodesOfApproximatedStructure[k].Coordinates[1]<<","<<NodesOfApproximatedStructure[k].Coordinates[2]<<std::endl;

                std::cout<<"Normal "<< k <<" "<< NodesOfApproximatedStructure[k].StructElemNormal[0]<<","<<NodesOfApproximatedStructure[k].StructElemNormal[1]<<","<<NodesOfApproximatedStructure[k].StructElemNormal[2]<<std::endl;
            }*/

        for (std::size_t i = 0; i < 2; i++)
            P_mean.Coordinates()[i] /= numberCutEdges;

        // Compute normal for the best-fitted plane
        array_1d<double, 3> N_mean;

        Matrix coordinates(numberCutEdges, 2);
        for (std::size_t i = 0; i < numberCutEdges; i++)
            for (std::size_t j = 0; j < 2; j++)
                coordinates(i, j) = NodesOfApproximatedStructure[i].Coordinates[j] - P_mean[j];

        Matrix A = prod(trans(coordinates), coordinates);
        Matrix V(2, 2);
        Vector lambda(2);

        // Calculate the eigenvectors V and the corresponding eigenvalues lambda
        EigenVectors(A, V, lambda);

        // Look for the minimal eigenvalue all lambdas
        std::size_t min_pos = 0;
        double min_lambda = lambda[min_pos];

        for (std::size_t i = 1; i < 2; i++)
            if (min_lambda > lambda[i])
            {
                min_lambda = lambda[i];
                min_pos = i;
            }

        // the normal equals to the eigenvector which corresponds to the minimal eigenvalue
        for (std::size_t i = 0; i < 2; i++)
            N_mean[i] = V(min_pos, i);

        N_mean[2] = 0.0;
        N_mean /= norm_2(N_mean);

        //nav
        /*std::cout<<"Min lamda "<<min_lambda<<std::endl;
        std::cout<<"Min eigenvector "<<V(min_pos,0)<<","<<V(min_pos,1)<<std::endl;
        std::cout<<"N_mean "<< N_mean[0]<<","<<N_mean[1]<<","<<N_mean[2]<<std::endl;*/
        // Check whether orientation of normal is in direction of the normal of the intersecting structure
        // Note: The normal of the approx. surface can be max. 90deg to every surrounding normal of the structure at the intersection nodes
        array_1d<double, 3> NormalAtOneIntersectionNode;
        NormalAtOneIntersectionNode = NodesOfApproximatedStructure[0].StructElemNormal;

        bool NormalWrongOriented = false;
        if (inner_prod(NormalAtOneIntersectionNode, N_mean) < 0)
            NormalWrongOriented = true;

        // switch direction of normal
        if (NormalWrongOriented)
            N_mean *= -1;

        // Determine about the minimal distance by considering the distances to both triangles
        for (std::size_t i_TriangleNode = 0; i_TriangleNode < 3; i_TriangleNode++)
        {
            ElementalDistances[i_TriangleNode] = PointDistanceToLine(P_mean, N_mean, i_fluid_Element->GetGeometry()[i_TriangleNode]);
        }

        // #################################################
        std::size_t numberDoubleCutEdges = 0;
        std::size_t indexDoubleCutEdge = 0;

        // figure out the edges which are cut more than once
        for (std::size_t i_TriangleEdge = 0; i_TriangleEdge < IntersectedTriangleEdges.size(); i_TriangleEdge++)
        {
            std::size_t NumberIntNodes = IntersectedTriangleEdges[i_TriangleEdge].IntNodes.size();
            if (NumberIntNodes == 2)
            {
                numberDoubleCutEdges++;
                indexDoubleCutEdge = i_TriangleEdge;
            }
        }

        if ((numberDoubleCutEdges >= 1))

        {

            array_1d<double, 3> normal_1 = IntersectedTriangleEdges[indexDoubleCutEdge].IntNodes[0].StructElemNormal;
            array_1d<double, 3> normal_2 = IntersectedTriangleEdges[indexDoubleCutEdge].IntNodes[1].StructElemNormal;

            // normalize normals
            normal_1 /= norm_2(normal_1);
            normal_2 /= norm_2(normal_2);

            const double pi = 3.1415926;

            // compute angle between normals
            double angle_n1n2 = acos(inner_prod(normal_1, normal_2));
            // rad --> degree
            angle_n1n2 *= 180 / pi;

            // if angle between -60º and 120º, take the mean
            if ((angle_n1n2 > -60) && (angle_n1n2 < 120))
            {
                // take the mean of the normals
                N_mean = 0.5 * (normal_1 + normal_2);
            }
            else
            {
                N_mean = 0.5 * (normal_1 - normal_2);
            }

            // Based on N_mean and P_mean compute the distances to that plane
            for (std::size_t i_TriangleNode = 0; i_TriangleNode < 3; i_TriangleNode++)
            {
                ElementalDistances[i_TriangleNode] = PointDistanceToLine(P_mean, N_mean, i_fluid_Element->GetGeometry()[i_TriangleNode]);
            }
        }
    }

    ///******************************************************************************************************************
    ///******************************************************************************************************************

    /**
       * This function calculates the distance of a 3D point to a line  spanned by a line segment
       * @param linie base point
       * @param lineNormal
       * @param ToPoint The point which distance is required
       * @return The distance between the point and the line spanned by the line segment
       */
    double PointDistanceToLine(Point &lineBasePoint,
                               array_1d<double, 3> &lineNormal,
                               Point &ToPoint)
    {
        // calculate vector pointing from a node in the line (e.g. line segment point 1) to the considered node ToPoint
        array_1d<double, 3> lineToPointVec = ToPoint - lineBasePoint;

        // projection of node on the plane
        const double sn = inner_prod(lineToPointVec, lineNormal);
        //std::cout<<"nav dot product"<<sn<<std::endl;
        const double sd = inner_prod(lineNormal, lineNormal);
        double DistanceToLine = sn / sqrt(sd);

        if (fabs(DistanceToLine) < epsilon)
            DistanceToLine = 0;

        return DistanceToLine;
    }

    ///******************************************************************************************************************
    ///******************************************************************************************************************

    void AssignMinimalNodalDistance()
    {
        // loop over all fluid Elements
        for (ModelPart::ElementIterator i_fluid_Element = mrFluidModelPart.ElementsBegin();
             i_fluid_Element != mrFluidModelPart.ElementsEnd();
             i_fluid_Element++)
        {
            Geometry<Node<3>> &geom = i_fluid_Element->GetGeometry();
            const Vector &ElementalDistances = i_fluid_Element->GetValue(ELEMENTAL_DISTANCES);

            // Assign distances to the single nodes, if a smaller value is found
            for (std::size_t i_TriangleNode = 0; i_TriangleNode < 3; i_TriangleNode++)
            {
                double currentNodeDist = geom[i_TriangleNode].GetSolutionStepValue(DISTANCE);
                double nodeInElemDist = ElementalDistances[i_TriangleNode];

                if (fabs(nodeInElemDist) < fabs(currentNodeDist))
                    geom[i_TriangleNode].GetSolutionStepValue(DISTANCE) = nodeInElemDist; // overwrite nodal distance (which is global)
                                                                                          //std::cout<<geom[i_TriangleNode].GetSolutionStepValue(DISTANCE)<<std::endl;
            }                                                                             // loop i_TriangleNode
        }                                                                                 // loop i_fluidElement
    }

    ///******************************************************************************************************************
    ///******************************************************************************************************************

    /**
       * If structure directly passes through the corner point of a trianglerahedra (leading to zero distances
       * in the respective node), then a small distance value (different from zero) will be stored for
       * that point. This is necessary since the embedded solver cannot handle zero distances.
       * @param Element            current Element which was cut by the structure (flag SPLIT_ELEMENT is set to one)
       * @param ElementalDistances Elemental distances calculated by the intersection pattern
       */
    void AvoidZeroDistances(ModelPart::ElementsContainerType::iterator &Element,
                            array_1d<double, 3> &ElementalDistances)
    {
        // Assign a distance limit
        //double dist_limit = 1e-5;*/
        bool distChangedToLimit = false; //variable to indicate that a distance value < tolerance is set to a limit distance = tolerance
                                         /*
        for (std::size_t i_node = 0; i_node < 3; i_node++)
        {
            if (fabs(ElementalDistances[i_node]) < dist_limit)
            {
                ElementalDistances[i_node] = dist_limit;
                distChangedToLimit = true;
            }
        }*/

        // Check, if this approach changes the split-flag (might be, that Element is not cut anymore if node with zero distance gets a positive limit distance value
        std::size_t numberNodesPositiveDistance = 0;

        for (std::size_t i_node = 0; i_node < 3; i_node++)
        {
            double &di = ElementalDistances[i_node];
            if (fabs(di) < dist_limit)
            {
                distChangedToLimit = true;
                /*if (di >= 0)
                {
                    di = dist_limit;

                }

                else
                    di = -dist_limit;
            }*/

                di = dist_limit;
            }
        }

        for (std::size_t i_node = 0; i_node < 3; i_node++)
        {
            if ((ElementalDistances[i_node]) > 0)
                numberNodesPositiveDistance++;
        }

        // Element is not set
        if (numberNodesPositiveDistance == 3 && distChangedToLimit == true)
            Element->GetValue(SPLIT_ELEMENT) = false;
    }

    ///******************************************************************************************************************
    ///******************************************************************************************************************

    void AvoidZeroNodalDistances()
    {

        ModelPart::NodesContainerType::ContainerType &nodes = mrFluidModelPart.NodesArray();

        // reset the node distance to 1.0 which is the maximum distance in our normalized space.
        int nodesSize = nodes.size();

#pragma omp parallel for firstprivate(nodesSize)
        for (int i = 0; i < nodesSize; i++)
        {
            double &distance = nodes[i]->GetSolutionStepValue(DISTANCE);

            if (fabs(distance) < dist_limit)
                distance = dist_limit;
        }
    }
    void GenerateSkinModelPart(ModelPart &mrNewSkinModelPart)
    {
        std::size_t id_node = 1;
        std::size_t id_condition = 1;

        mrNewSkinModelPart.Nodes().reserve(mrFluidModelPart.Nodes().size());
        mrNewSkinModelPart.Conditions().reserve(mrFluidModelPart.Elements().size());

        for (ModelPart::ElementIterator i_fluid_element = mrFluidModelPart.ElementsBegin();
             i_fluid_element != mrFluidModelPart.ElementsEnd();
             i_fluid_element++)
        {
            bool is_split = i_fluid_element->GetValue(SPLIT_ELEMENT);
            if (is_split == true)
            {
                const Vector &distances = i_fluid_element->GetValue(ELEMENTAL_DISTANCES);
                Geometry<Node<3>> &geom = i_fluid_element->GetGeometry();

                // generate the points on the edges at the zero of the distance function
                std::vector<Point> edge_points;
                edge_points.reserve(2);

                // loop over all 6 edges of the triangle
                for (std::size_t i = 0; i < 2; i++)
                {
                    for (std::size_t j = i + 1; j < 3; j++) // go through the edges 01, 02, 03, 12
                    {
                        double di = distances[i];
                        double dj = distances[j];

                        if (di * dj < 0) //edge is cut
                        {
                            // generate point on edge by linear interpolation
                            double Ni = fabs(dj) / (fabs(di) + fabs(dj));
                            double Nj = 1.0 - Ni;
                            Point edge_point(Ni * geom[i] + Nj * geom[j]);
                            edge_points.push_back(edge_point);
                        }
                    }
                }

                // three intersection nodes
                if (edge_points.size() == 2)
                {
                    // ######## ADDING NEW NODE #########
                    Node<3>::Pointer pnode1 = mrNewSkinModelPart.CreateNewNode(id_node++, edge_points[0].X(), edge_points[0].Y(), edge_points[0].Z());
                    Node<3>::Pointer pnode2 = mrNewSkinModelPart.CreateNewNode(id_node++, edge_points[1].X(), edge_points[1].Y(), edge_points[1].Z());

                    // ######## ADDING NEW CONDITION #########
                    //form a triangle
                    Line2D2<Node<3>> line(pnode1, pnode2);

                    Condition const &rReferenceCondition = KratosComponents<Condition>::Get("Condition2D");
                    Properties::Pointer properties = mrNewSkinModelPart.rProperties()(0);
                    Condition::Pointer p_condition = rReferenceCondition.Create(id_condition++, line, properties);

                    mrNewSkinModelPart.Conditions().push_back(p_condition);
                }

                else
                {
                    std::cout << "Error :: Found more than two intersection in the triagular element" << std::endl;
                    std::exit(-1);
                }

                // four intersection nodes
                /*    if (edge_points.size() == 4)
                {
                    //form a quadrilatera with the 4 cut nodes
                    array_1d<double, 3> x21 = edge_points[1] - edge_points[0];
                    array_1d<double, 3> x31 = edge_points[2] - edge_points[0];
                    array_1d<double, 3> x41 = edge_points[3] - edge_points[0];

                    //define a vector oriented as x21
                    array_1d<double, 3> v1 = x21 / norm_2(x21);

                    boost::numeric::ublas::bounded_matrix<double, 4, 3> DN_DX;
                    array_1d<double, 4> msN;
                    double Area;
                    GeometryUtils::CalculateGeometryData(geom, DN_DX, msN, Area);

                    array_1d<double, 3> n = prod(trans(DN_DX), distances);
                    n /= norm_2(n);

                    array_1d<double, 3> v2 = MathUtils<double>::CrossProduct(n, v1);

                    array_1d<double, 3> angles;
                    angles[0] = 0.0;                                             //angle between x21 and v1
                    angles[1] = atan2(inner_prod(x31, v2), inner_prod(x31, v1)); //angle between x31 and v1
                    angles[2] = atan2(inner_prod(x41, v2), inner_prod(x41, v1)); //angle between x31 and v1

                    double max_angle = 0.0;
                    double min_angle = 0.0;
                    std::size_t min_pos = 1;
                    std::size_t max_pos = 1;
                    for (std::size_t i = 1; i < 3; i++)
                    {
                        if (angles[i] < min_angle)
                        {
                            min_pos = i + 1; //this is the local index of the edge point which forms the minimal angle
                            min_angle = angles[i];
                        }
                        else if (angles[i] > max_angle)
                        {
                            max_pos = i + 1; //this is the local index of the edge point which forms the maximal angle
                            max_angle = angles[i];
                        }
                    }

                    //find the pos of the center node
                    std::size_t center_pos = 0;
                    for (std::size_t i = 1; i < 4; i++)
                    {
                        if ((i != min_pos) && (i != max_pos))
                        {
                            center_pos = i;
                        }
                    }

                    // ######## ADDING NEW NODE #########
                    Node<3>::Pointer pnode1 = mrNewSkinModelPart.CreateNewNode(id_node++, edge_points[0].X(), edge_points[0].Y(), edge_points[0].Z());
                    Node<3>::Pointer pnode2 = mrNewSkinModelPart.CreateNewNode(id_node++, edge_points[min_pos].X(), edge_points[min_pos].Y(), edge_points[min_pos].Z());
                    Node<3>::Pointer pnode3 = mrNewSkinModelPart.CreateNewNode(id_node++, edge_points[center_pos].X(), edge_points[center_pos].Y(), edge_points[center_pos].Z());
                    Node<3>::Pointer pnode4 = mrNewSkinModelPart.CreateNewNode(id_node++, edge_points[max_pos].X(), edge_points[max_pos].Y(), edge_points[max_pos].Z());

                    // ######## ADDING NEW CONDITION #########
                    //form two triangles
                    Triangle3D3<Node<3>> triangle1(pnode1, pnode2, pnode3);
                    Triangle3D3<Node<3>> triangle2(pnode1, pnode3, pnode4);

                    Condition const &rReferenceCondition = KratosComponents<Condition>::Get("Condition3D");

                    Properties::Pointer properties = mrNewSkinModelPart.rProperties()(0);

                    Condition::Pointer p_condition1 = rReferenceCondition.Create(id_condition++, triangle1, properties);
                    Condition::Pointer p_condition2 = rReferenceCondition.Create(id_condition++, triangle2, properties);

                    mrNewSkinModelPart.Conditions().push_back(p_condition1);
                    mrNewSkinModelPart.Conditions().push_back(p_condition2);
                }*/
            }
        }
    }

    ///******************************************************************************************************************
    ///******************************************************************************************************************

    void GenerateQuadtree()
    {
        Timer::Start("Generating Quadtree");
        //std::cout << "Generating the Quadtree..." << std::endl;
        boost::shared_ptr<QuadtreeType> temp_quadtree = boost::shared_ptr<QuadtreeType>(new QuadtreeType());
        //QuadtreeType::Pointer temp_quadtree = QuadtreeType::Pointer(new QuadtreeType() );
        mpQuadtree.swap(temp_quadtree);

        std::cout << "Inside Quadtree" << std::endl;

        double low[2];
        double high[2];

        std::cout << "Fluid Model part" << mrFluidModelPart << std::endl;

        for (int i = 0; i < 2; i++)
        {
            low[i] = high[i] = mrFluidModelPart.NodesBegin()->Coordinates()[i];
        }

        // loop over all nodes in the bounding box
        for (ModelPart::NodeIterator i_node = mrFluidModelPart.NodesBegin();
             i_node != mrFluidModelPart.NodesEnd();
             i_node++)
        {
            for (int i = 0; i < 2; i++)
            {
                low[i] = i_node->Coordinates()[i] < low[i] ? i_node->Coordinates()[i] : low[i];
                high[i] = i_node->Coordinates()[i] > high[i] ? i_node->Coordinates()[i] : high[i];
            }
        }

        // loop over all skin nodes
        for (ModelPart::NodeIterator i_node = mrSkinModelPart.NodesBegin();
             i_node != mrSkinModelPart.NodesEnd();
             i_node++)
        {
            for (int i = 0; i < 2; i++)
            {
                low[i] = i_node->Coordinates()[i] < low[i] ? i_node->Coordinates()[i] : low[i];
                high[i] = i_node->Coordinates()[i] > high[i] ? i_node->Coordinates()[i] : high[i];
            }
        }

        ///rishtih
        /*
        low[0] -=  0.3*(high[0]-low[0] );
        low[1] -=  0.3*(high[1]-low[1] );

        high[0] +=  0.3*(high[0]-low[0] );
        high[1] +=  0.3*(high[1]-low[1] );
 */

        std::cout << "Skin added" << std::endl;
        mpQuadtree->SetBoundingBox(low, high);
        std::cout << "Bounding Box Lowest Dimension" << low[0] << "," << low[1] << std::endl;

        std::cout << "Bounding Box Highest Dimension" << high[0] << "," << high[1] << std::endl;

        //mpQuadtree->RefineWithUniformSize(0.0625);

        // loop over all structure nodes
        for (ModelPart::NodeIterator i_node = mrSkinModelPart.NodesBegin();
             i_node != mrSkinModelPart.NodesEnd();
             i_node++)
        {
            double temp_point[2];
            temp_point[0] = i_node->X();
            temp_point[1] = i_node->Y();

            mpQuadtree->Insert(temp_point);
        }

        // loop over all structure elements
        for (ModelPart::ConditionIterator i_cond = mrSkinModelPart.ConditionsBegin(); //nav
             i_cond != mrSkinModelPart.ConditionsEnd();                               //nav
             i_cond++)
        {
            mpQuadtree->Insert(*(i_cond).base());
        }

        Timer::Stop("Generating Quadtree");
        std::cout << "Quadtree generation finished" << std::endl;
        KRATOS_WATCH(mpQuadtree);

        /*GenerateNodes();

        std::cout << "######## WRITING QUADTREE MESH #########" << std::endl;
        std::ofstream myfile;
        myfile.open ("quadtree.post.msh");
        mpQuadtree->PrintGiDMeshNew(myfile);
        myfile.close();*/

        std::cout << "Generating the Quadtree finished" << std::endl;
    }

    ///******************************************************************************************************************
    ///******************************************************************************************************************

    /*void GenerateNodes()
    {
        Timer::Start("Generating Nodes");
        std::vector<QuadtreeType::cell_type*> all_leaves;
        mpQuadtree->GetAllLeavesVector(all_leaves);

        int leaves_size = all_leaves.size();

#pragma omp parallel for
        for (int i = 0; i < leaves_size; i++)
        {
            *(all_leaves[i]->pGetDataPointer()) = ConfigurationType::AllocateData();
        }


        std::size_t last_id = mrSkinModelPart.NumberOfNodes() +mrFluidModelPart.numberOfNodes()+ 1;

        for (std::size_t i = 0; i < all_leaves.size(); i++)
        {
                    CellType* cell = all_leaves[i];
            GenerateCellNode(cell, last_id);
        }

        Timer::Stop("Generating Nodes");

    }*/

    ///******************************************************************************************************************
    ///******************************************************************************************************************

    void GenerateCellNode(CellType *pCell, std::size_t &LastId)
    {

        for (int i_pos = 0; i_pos < 4; i_pos++) // position 4 is for center
        {

            DistanceSpatialContainersConditionConfigure2d::cell_node_data_type *p_node = (*(pCell->pGetData()))[i_pos];

            if (p_node == 0)
            {
                (*(pCell->pGetData()))[i_pos] = new DistanceSpatialContainersConditionConfigure2d::cell_node_data_type;

                (*(pCell->pGetData()))[i_pos]->Id() = LastId++;

                mQuadtreeNodes.push_back((*(pCell->pGetData()))[i_pos]);

                SetNodeInNeighbours(pCell, i_pos, (*(pCell->pGetData()))[i_pos]);
            }
        }
    }

    ///******************************************************************************************************************
    ///******************************************************************************************************************

    void SetNodeInNeighbours(CellType *pCell, int Position, CellNodeDataType *pNode)
    {
        CellType::key_type point_key[2];
        pCell->GetKey(Position, point_key);

        for (std::size_t i_direction = 0; i_direction < 4; i_direction++)
        {
            CellType::key_type neighbour_key[2];
            if (pCell->GetNeighbourKey(Position, i_direction, neighbour_key))
            {
                CellType *neighbour_cell = mpQuadtree->pGetCell(neighbour_key);
                if (!neighbour_cell || (neighbour_cell == pCell))
                    continue;

                std::size_t position = neighbour_cell->GetLocalPosition(point_key);
                if ((*neighbour_cell->pGetData())[position])
                {
                    std::cout << "ERROR!! Bad Position calculated!!!!!!!!!!! position :" << position << std::endl;
                    continue;
                }

                (*neighbour_cell->pGetData())[position] = pNode;
            }
        }
    }

    ///******************************************************************************************************************
    ///******************************************************************************************************************

    void CalculateDistance2()
    {

        Timer::Start("Calculate Distances2");
        ModelPart::NodesContainerType::ContainerType &nodes = mrFluidModelPart.NodesArray();
        int nodes_size = nodes.size();
        //         // first of all we reset the node distance to 1.00 which is the maximum distnace in our normalized space.
        //#pragma omp parallel for firstprivate(nodes_size)
        //         for(int i = 0 ; i < nodes_size ; i++)
        //             nodes[i]->GetSolutionStepValue(DISTANCE) = 1.00;

        //std::vector<CellType*> leaves;

        //mpQuadtree->GetAllLeavesVector(leaves);
        //int leaves_size = leaves.size();

        //         for(int i = 0 ; i < leaves_size ; i++)
        //             CalculateNotEmptyLeavesDistance(leaves[i]);

#pragma omp parallel for firstprivate(nodes_size)
        for (int i = 0; i < nodes_size; i++)
        {

            /*std::cout<<"nav Node Id"<<nodes[i]->Id()<<std::endl;*/

            //debug

            CalculateNodeDistance(*(nodes[i]));
        }
        Timer::Stop("Calculate Distances2");
    }

    ///******************************************************************************************************************
    ///******************************************************************************************************************

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

    //         mpQuadtree->GetAllLeavesVector(leaves);
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

    //         mpQuadtree->GetAllLeavesVector(leaves);
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
        DistanceSpatialContainersConditionConfigure2d::data_type &nodes = mQuadtreeNodes;
        int nodes_size = nodes.size();
// first of all we reste the node distance to 1.00 which is the maximum distnace in our normalized space.
#pragma omp parallel for firstprivate(nodes_size)
        for (int i = 0; i < nodes_size; i++)
            nodes[i]->Distance() = 1.00;

        std::vector<CellType *> leaves;

        mpQuadtree->GetAllLeavesVector(leaves);
        int leaves_size = leaves.size();

        for (int i = 0; i < leaves_size; i++)
            CalculateNotEmptyLeavesDistance(leaves[i]);

        for (int i_direction = 0; i_direction < 1; i_direction++)
        {

            //#pragma omp parallel for firstprivate(nodes_size)
            for (int i = 0; i < nodes_size; i++)
            {
                if (nodes[i]->X() < 1.00 && nodes[i]->Y() < 1.00 && nodes[i]->Z() < 1.00)
                    // if((*nodes[i])[i_direction] == 0.00)
                    CalculateDistance(*(nodes[i]), i_direction);
            }
        }
        Timer::Stop("Calculate Distances");
    }

    void CalculateDistance(CellNodeDataType &rNode, int i_direction)
    {

        double coords[3] = {rNode.X(), rNode.Y(), rNode.Z()};
        // KRATOS_WATCH_3(coords);

        //This function must color the positions in space defined by 'coords'.
        //coords is of dimension (3) normalized in (0,1)^3 space

        typedef Element::GeometryType triangle_type;
        typedef std::vector<std::pair<double, triangle_type *>> intersections_container_type;

        intersections_container_type intersections;
        DistanceSpatialContainersConditionConfigure2d::data_type nodes_array;

        const double epsilon = 1e-12;

        double distance = 1.0;

        // Creating the ray
        double ray[3] = {coords[0], coords[1], coords[2]};

        mpQuadtree->NormalizeCoordinates(ray);
        ray[i_direction] = 0; // starting from the lower extreme

        //            KRATOS_WATCH_3(ray)
        GetIntersectionsAndNodes(ray, i_direction, intersections, nodes_array);
        //            KRATOS_WATCH(nodes_array.size())
        for (std::size_t i_node = 0; i_node < nodes_array.size(); i_node++)
        {
            double coord = nodes_array[i_node]->Coordinate(i_direction + 1);
            //             KRATOS_WATCH(intersections.size());

            int ray_color = 1;
            std::vector<std::pair<double, Element::GeometryType *>>::iterator i_intersection = intersections.begin();
            while (i_intersection != intersections.end())
            {
                double d = coord - i_intersection->first;
                if (d > epsilon)
                {

                    ray_color = -ray_color;
                    distance = d;
                }
                else if (d > -epsilon)
                { //interface
                    distance = 0.00;
                    break;
                }
                else
                {
                    if (distance > -d)
                        distance = -d;
                    break;
                }

                i_intersection++;
            }

            distance *= ray_color;

            double &node_distance = nodes_array[i_node]->Distance();
            if (fabs(distance) < fabs(node_distance))
                node_distance = distance;
            else if (distance * node_distance < 0.00) // assigning the correct sign
                node_distance = -node_distance;
        }
    }

    void CalculateNotEmptyLeavesDistance(CellType *pCell)
    {
        //typedef Element::GeometryType triangle_type;
        typedef QuadtreeType::cell_type::object_container_type object_container_type;

        object_container_type *objects = (pCell->pGetObjects());

        // There are no intersection in empty cells
        if (objects->empty())
            return;

        for (int i_pos = 0; i_pos < 8; i_pos++) // position 8 is for center
        {
            double distance = 1.00; // maximum distance is 1.00

            for (object_container_type::iterator i_object = objects->begin(); i_object != objects->end(); i_object++)
            {
                CellType::key_type keys[3];
                pCell->GetKey(i_pos, keys);

                double cell_point[3];
                mpQuadtree->CalculateCoordinates(keys, cell_point);

                //                cell_point[0] = pCell->GetCoordinate(keys[0]);
                //                cell_point[1] = pCell->GetCoordinate(keys[1]);
                //                cell_point[2] = pCell->GetCoordinate(keys[2]);

                double d = GeometryUtils::PointDistanceToTriangle3D((*i_object)->GetGeometry()[0], (*i_object)->GetGeometry()[1], (*i_object)->GetGeometry()[2], Point(cell_point[0], cell_point[1], cell_point[2]));

                if (d < distance)
                    distance = d;
            }

            double &node_distance = (*(pCell->pGetData()))[i_pos]->Distance();
            if (distance < node_distance)
                node_distance = distance;
        }
    }

    void CalculateNodeDistance(Node<3> &rNode)
    {
        double coord[3] = {rNode.X(), rNode.Y(), rNode.Z()};
        double distance = DistancePositionInSpace(coord);
        double &node_distance = rNode.GetSolutionStepValue(DISTANCE);

        //debug
        /*if(rNode.Id() == 1481)
        {
            std::cout<<" previous distance "<<node_distance<<std::endl;
            std::cout<<" From quadtree"<< distance<<std::endl;
            std::exit(-1);

        }*/

        if (fabs(node_distance) > fabs(distance) || fabs(node_distance < dist_limit))
            node_distance = distance;
        else if (distance * node_distance < 0.00) // assigning the correct sign
            node_distance = -node_distance;
    }

    //      void CalculateNodeDistanceFromCell(Node<3>& rNode)
    //      {
    //          QuadtreeType::key_type node_key[3] = {quadtree->CalcKeyNormalized(rNode.X()), quadtree->CalcKeyNormalized(rNode.Y()), quadtree->CalcKeyNormalized(rNode.Z())};
    //          QuadtreeType::cell_type* pcell = quadtree->pGetCell(node_key);

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

    double DistancePositionInSpace(double *coords)
    {
        //This function must color the positions in space defined by 'coords'.
        //coords is of dimension (3) normalized in (0,1)^3 space

        typedef Element::GeometryType line_type;
        typedef std::vector<std::pair<double, line_type *>> intersections_container_type;

        intersections_container_type intersections;

        const int dimension = 2;
        const double epsilon = 1e-12;

        double distances[2] = {10000, 10000};

        for (int i_direction = 0; i_direction < dimension; i_direction++)
        {
            // Creating the ray
            double ray[2] = {coords[0], coords[1]};

            /*std::cout<<"navray "<<ray[0]<<","<<ray[1]<<std::endl;*/

            mpQuadtree->NormalizeCoordinates(ray);

            ray[i_direction] = 0; // starting from the lower extreme

            GetIntersections(ray, i_direction, intersections);

            //                if(intersections.size() == 1)
            //                    KRATOS_WATCH_3(ray)

            //             KRATOS_WATCH(intersections.size());

            int ray_color = 1;
            std::vector<std::pair<double, Element::GeometryType *>>::iterator i_intersection = intersections.begin();
            while (i_intersection != intersections.end())
            {
                double d = coords[i_direction] - i_intersection->first;

                if (d > epsilon)
                {

                    ray_color = -ray_color;
                    distances[i_direction] = d;
                    //                        if(distances[i_direction] > d) // I think this is redundunt. Pooyan.
                    //                        {
                    //                            if(ray_color > 0.00)
                    //                                    distances[i_direction] = d;
                    //                            else
                    //                                distances[i_direction] = -d;
                    //                        }
                }
                else if (d > -epsilon)
                { //interface
                    distances[i_direction] = 0.00;
                    break;
                }
                else
                {
                    if (distances[i_direction] > -d)
                        distances[i_direction] = -d;
                    break;
                }

                i_intersection++;
            }

            distances[i_direction] *= ray_color;
        }

        //            if(distances[0]*distances[1] < 0.00 || distances[2]*distances[1] < 0.00)
        //                KRATOS_WATCH_3(distances);

        //#ifdef _DEBUG
        //        std::cout << "colors : " << colors[0] << ", " << colors[1] << ", " << colors[2] << std::endl;
        //#endif
        double distance = (fabs(distances[0]) > fabs(distances[1])) ? distances[1] : distances[0];

         if (distances[0] * distances[1] < 0)
            distance = fabs(distance);
        //distance = (fabs(distance) > fabs(distances[2])) ? distances[2] : distance;

        return distance;
    }

    void GetIntersectionsAndNodes(double *ray, int direction, std::vector<std::pair<double, Element::GeometryType *>> &intersections, DistanceSpatialContainersConditionConfigure2d::data_type &rNodesArray)
    {
        //This function passes the ray through the model and gives the hit point to all objects in its way
        //ray is of dimension (3) normalized in (0,1)^3 space
        // direction can be 0,1,2 which are x,y and z respectively

        const double epsilon = 1.00e-12;

        // first clearing the intersections points vector
        intersections.clear();

        //QuadtreeType* quadtree = &mQuadtree;
        QuadtreeType *quadtree = mpQuadtree.get();

        QuadtreeType::key_type ray_key[3] = {quadtree->CalcKeyNormalized(ray[0]), quadtree->CalcKeyNormalized(ray[1]), quadtree->CalcKeyNormalized(ray[2])};
        QuadtreeType::key_type cell_key[3];

        // getting the entrance cell from lower extreme
        ray_key[direction] = 0;
        QuadtreeType::cell_type *cell = quadtree->pGetCell(ray_key);

        while (cell)
        {
            std::size_t position = cell->GetLocalPosition(ray_key); // Is this the local position!?!?!?!
            QuadtreeType::key_type node_key[3];
            cell->GetKey(position, node_key);
            if ((node_key[0] == ray_key[0]) && (node_key[1] == ray_key[1]) && (node_key[2] == ray_key[2]))
            {
                if (cell->pGetData())
                {
                    if (cell->pGetData()->size() > position)
                    {
                        CellNodeDataType *p_node = (*cell->pGetData())[position];
                        if (p_node)
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
            if (cell->GetNeighbourKey(1 + direction * 2, cell_key))
            {
                ray_key[direction] = cell_key[direction];
                cell = quadtree->pGetCell(ray_key);
                ray_key[direction] -= 1; //the key returned by GetNeighbourKey is inside the cell (minkey +1), to ensure that the corresponding
                //cell get in pGetCell is the right one.
                //#ifdef _DEBUG
                //                Quadtree_Pooyan::key_type min_key[3];
                //                cell->GetMinKey(min_key[0],min_key[1],min_key[2]);
                //                Quadtree_Pooyan::key_type tmp;
                //                tmp= min_key[direction];
                //                assert(ray_key[direction]==tmp);
                //#endif
            }
            else
                cell = NULL;
        }

        //   KRATOS_WATCH(rNodesArray.size());
        // now eliminating the repeated objects
        if (!intersections.empty())
        {
            //sort
            std::sort(intersections.begin(), intersections.end());
            // unique
            std::vector<std::pair<double, Element::GeometryType *>>::iterator i_begin = intersections.begin();
            std::vector<std::pair<double, Element::GeometryType *>>::iterator i_intersection = intersections.begin();
            while (++i_begin != intersections.end())
            {
                // considering the very near points as the same points
                if (fabs(i_begin->first - i_intersection->first) > epsilon) // if the hit points are far enough they are not the same
                    *(++i_intersection) = *i_begin;
            }
            intersections.resize((++i_intersection) - intersections.begin());
        }
    }

    void GetIntersections(double *ray, int direction, std::vector<std::pair<double, Element::GeometryType *>> &intersections)
    {
        //This function passes the ray through the model and gives the hit point to all objects in its way
        //ray is of dimension (2) normalized in (0,1)^2 space
        // direction can be 0,1which are x,y respectively

        const double epsilon = 1.00e-12;

        // first clearing the intersections points vector
        intersections.clear();

        //QuadtreeType* quadtree = &mQuadtree;
        QuadtreeType *quadtree = mpQuadtree.get();

        QuadtreeType::key_type ray_key[2] = {quadtree->CalcKeyNormalized(ray[0]), quadtree->CalcKeyNormalized(ray[1])};
        QuadtreeType::key_type cell_key[2];

        // getting the entrance cell from lower extreme
        QuadtreeType::cell_type *cell = quadtree->pGetCell(ray_key);

        while (cell)
        {
            //        std::cout << ".";
            GetCellIntersections(cell, ray, ray_key, direction, intersections);
            // go to the next cell
            if (cell->GetNeighbourKey(1 + direction * 2, cell_key))
            {
                ray_key[direction] = cell_key[direction];
                cell = quadtree->pGetCell(ray_key);
                ray_key[direction] -= 1; //the key returned by GetNeighbourKey is inside the cell (minkey +1), to ensure that the corresponding
                //cell get in pGetCell is the right one.
                //#ifdef _DEBUG
                //                Quadtree_Pooyan::key_type min_key[3];
                //                cell->GetMinKey(min_key[0],min_key[1],min_key[2]);
                //                Quadtree_Pooyan::key_type tmp;
                //                tmp= min_key[direction];
                //                assert(ray_key[direction]==tmp);
                //#endif
            }
            else
                cell = NULL;
        }

        // now eliminating the repeated objects
        if (!intersections.empty())
        {
            //sort
            std::sort(intersections.begin(), intersections.end());
            // unique
            std::vector<std::pair<double, Element::GeometryType *>>::iterator i_begin = intersections.begin();
            std::vector<std::pair<double, Element::GeometryType *>>::iterator i_intersection = intersections.begin();
            while (++i_begin != intersections.end())
            {
                // considering the very near points as the same points
                if (fabs(i_begin->first - i_intersection->first) > epsilon) // if the hit points are far enough they are not the same
                    {

                            *(++i_intersection) = *i_begin;
                    }

                //Not yet tested , this can occur in some rare cases (depending on mesh):  once it happens we need to test this code
                // else
                //{
                   // std::cout<<" test this area to solve the distance calculation isse"<<std::endl;

                    /* double Line1 = (fabs(i_begin->first - i_begin->second->Points()[0][direction]) > fabs(i_begin->first - i_begin->second->Points()[1][direction])) ? i_begin->first - i_begin->second->Points()[0][direction] : i_begin->first - i_begin->second->Points()[1][direction];
                    double Line2 = (fabs(i_intersection->first - i_intersection->second->Points()[0][direction]) > fabs(i_intersection->first - i_intersection->second->Points()[1][direction])) ? i_intersection->first - i_intersection->second->Points()[0][direction] : i_intersection->first - i_intersection->second->Points()[1][direction];

                    std::cout<<"##########################################################################"<<Line1<<std::endl;
                    std::cout<<"##########################################################################"<<Line2<<std::endl;
                    if( (Line1>0 && Line2>0) || (Line1<0 && Line2<0) )
                    {
                        std::cout<<"Encountered a corner case beacause of this particular mesh, this part of the code taes care of it. But the code isnt tested yet"<<std::endl;
                        std::cout<<i_begin->first<<"######"<<std::endl;
                        std::cout<<i_begin->second->Points()[0][direction]<<std::endl;
                        std::cout<<i_begin->second->Points()[1][direction]<<std::endl;
                        *i_intersection=*(i_begin+1);
                        *i_begin = *(i_begin+1);
                        i_begin--;
                        std::exit(-1);

                    } */

               // }
            }
            intersections.resize((++i_intersection) - intersections.begin());

        }
    }

    int GetCellIntersections(QuadtreeType::cell_type *cell, double *ray,
                             QuadtreeType::key_type *ray_key, int direction,
                             std::vector<std::pair<double, Element::GeometryType *>> &intersections)
    {
        //This function passes the ray through the cell and gives the hit point to all objects in its way
        //ray is of dimension (2) normalized in (0,1)^2 space
        // direction can be 0,1 which are x,y  respectively

        //typedef Element::GeometryType triangle_type;
        typedef QuadtreeType::cell_type::object_container_type object_container_type;

        object_container_type *objects = (cell->pGetObjects());

        // There are no intersection in empty cells
        if (objects->empty())
            return 0;

        //      std::cout << "X";
        // calculating the two extreme of the ray segment inside the cell
        double ray_point1[3] = {ray[0], ray[1]};
        double ray_point2[3] = {ray[0], ray[1]};
        double normalized_coordinate;
        mpQuadtree->CalculateCoordinateNormalized(ray_key[direction], normalized_coordinate);
        ray_point1[direction] = normalized_coordinate;
        ray_point2[direction] = ray_point1[direction] + mpQuadtree->CalcSizeNormalized(cell);

        mpQuadtree->ScaleBackToOriginalCoordinate(ray_point1);
        mpQuadtree->ScaleBackToOriginalCoordinate(ray_point2);

        for (object_container_type::iterator i_object = objects->begin(); i_object != objects->end(); i_object++)
        {
            double intersection[3] = {0.00, 0.00, 0.00};

            int is_intersected = IntersectionLineSegment((*i_object)->GetGeometry(), ray_point1, ray_point2, intersection); // This intersection has to be optimized for axis aligned rays

            if (is_intersected == 1) // There is an intersection but not coplanar
            {

                /*                std::cout<<"navIntersection happened"<<std::endl;
                std::cout<<"navPoint 1 "<<(*i_object)->GetGeometry()[0][0]<<","<<(*i_object)->GetGeometry()[0][1]<<std::endl;
                std::cout<<"navPoint 2 "<<(*i_object)->GetGeometry()[1][0]<<","<<(*i_object)->GetGeometry()[1][1]<<std::endl;
                std::cout<<"navIntersection  "<<intersection[direction]<<" Direction "<<direction<<std::endl;*/
                //KRATOS_WATCH(intersection[direction]);

                //KRATOS_WATCH(((*i_object)->GetGeometry()));



                intersections.push_back(std::pair<double, Element::GeometryType *>(intersection[direction], &((*i_object)->GetGeometry())));
            }

            //else if(is_intersected == 2) // coplanar case
        }

        return 0;
    }

    int IntersectionLineSegment(Element::GeometryType &rGeometry, double *RayPoint1, double *RayPoint2, double *IntersectionPoint)
    {
        // This is the adaption of the algorithm provided in:
        // http://mathworld.wolfram.com/Line-LineIntersection.html

        /*               D
                        /
                       /
                      / I
                A----o-------B
                    /
                   /
                  /
                 C
*/

        array_1d<double, 3> vec_a, vec_ab;               // struc line segment vector vec_a + s vec_ab
        array_1d<double, 3> vec_c, vec_cd;               // fluid line edge vector vec_c + t vec_cd
        double s, t;                                     // params to calc ray-ray intersection
        array_1d<double, 3> vec_ac, vec_g, vec_h, vec_k; // intermediate vector for intersection calculation

        // get struc edge vectors

        vec_a = rGeometry[0];
        vec_ab = rGeometry[1] - rGeometry[0];

        if (norm_2(vec_ab) == 0) //line segment is degenerate
            return -1;           // do not deal with this case

        // get fluid edge vectors
        for (int i = 0; i < 3; i++)
        {
            vec_c[i] = RayPoint1[i];
            vec_cd[i] = RayPoint2[i] - RayPoint1[i]; // ray direction vector CD
        }

        vec_ac = vec_c - vec_a;

        MathUtils<double>::CrossProduct(vec_g, vec_cd, vec_ac);
        MathUtils<double>::CrossProduct(vec_h, vec_ab, vec_ac);
        MathUtils<double>::CrossProduct(vec_k, vec_cd, vec_ab);
        double mag_cd, gk, hk, kk;
        mag_cd = norm_2(vec_cd);
        double mag_ab = norm_2(vec_ab);

        gk = inner_prod(vec_g, vec_k);
        hk = inner_prod(vec_h, vec_k);
        kk = inner_prod(vec_k, vec_k);

        //if (kk < epsilon)
        if (kk/(mag_ab*mag_cd) < epsilon)
        {                    // ray is parallel to the struc line segment
            /* if (mag_cd == 0) // ray is coincident with struc edge
                return 2;
            else */
                return 0; // ray disjoint from struc edge
        }
        // get intersect point of ray with the struc edge
        s = gk / kk;
        t = hk / kk;
        //if (r < 0.0)                   // ray goes away from triangle
        //return 0;                  // => no intersect
        // for a segment, also test if (r > 1.0) => no intersect

        for (int i = 0; i < 3; i++)
        {

            IntersectionPoint[i] = vec_a[i] + s * vec_ab[i]; // intersect point of ray and struc edge
        }

        // is I inside line segment?

        if (s < 0.0 - epsilon || s > 1.0 + epsilon) // I is outside struc edge
            return 0;

        if (t < 0.0 - epsilon || t > 1.0 + epsilon) // I is outside fluid edge
            return 0;

        /*std::cout<<"str edge "<<rGeometry[0][0]<<","<<rGeometry[0][1]<<")  ("<<rGeometry[1][0]<<","<<rGeometry[1][1]<<std::endl;
        std::cout<<"Fluid edge "<<RayPoint1[0]<<","<<RayPoint1[1]<<")  ("<<RayPoint2[0]<<","<<RayPoint2[1]<<std::endl;
        std::cout<<" s "<<s<<" t "<<t<<" gk "<<gk<<" hk "<<hk<<" kk "<<kk<<std::endl;
        std::cout<<" g "<<vec_g[0]<<","<<vec_g[1]<<","<<vec_g[2]<<std::endl;
        std::cout<<" h "<<vec_h[0]<<","<<vec_h[1]<<","<<vec_h[2]<<std::endl;
        std::cout<<" k "<<vec_k[0]<<","<<vec_k[1]<<","<<vec_k[2]<<std::endl;
        std::cout<<" Intersection"<<IntersectionPoint[0]<<","<<IntersectionPoint[1]<<std::endl;
        std::exit(-1);*/

        return 1; // I is in struc edge
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
    virtual std::string Info() const override
    {
        return "CalculateSignedDistanceTo2DConditionSkinProcess";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream &rOStream) const override
    {
        rOStream << "CalculateSignedDistanceTo2DConditionSkinProcess";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream &rOStream) const override
    {
    }

    void PrintGiDMesh(std::ostream &rOStream) const
    {
        std::vector<CellType *> leaves;

        mpQuadtree->GetAllLeavesVector(leaves);

        std::cout << "writing " << leaves.size() << " leaves" << std::endl;
        rOStream << "MESH \"leaves\" dimension 3 ElemType Hexahedra Nnode 8" << std::endl;
        rOStream << "# color 96 96 96" << std::endl;
        rOStream << "Coordinates" << std::endl;
        rOStream << "# node number coordinate_x coordinate_y coordinate_z  " << std::endl;

        for (DistanceSpatialContainersConditionConfigure2d::data_type::const_iterator i_node = mQuadtreeNodes.begin(); i_node != mQuadtreeNodes.end(); i_node++)
        {
            rOStream << (*i_node)->Id() << "  " << (*i_node)->Coordinate(1) << "  " << (*i_node)->Coordinate(2) << "  " << (*i_node)->Coordinate(3) << std::endl;
            //mpQuadtree->Insert(temp_point);
        }
        std::cout << "Nodes written..." << std::endl;
        rOStream << "end coordinates" << std::endl;
        rOStream << "Elements" << std::endl;
        rOStream << "# Element node_1 node_2 node_3 material_number" << std::endl;

        for (std::size_t i = 0; i < leaves.size(); i++)
        {
            if ((leaves[i]->pGetData()))
            {
                DistanceSpatialContainersConditionConfigure2d::data_type &nodes = (*(leaves[i]->pGetData()));

                rOStream << i + 1;
                for (int j = 0; j < 8; j++)
                    rOStream << "  " << nodes[j]->Id();
                rOStream << std::endl;
            }
        }
        rOStream << "end Elements" << std::endl;
    }

    void PrintGiDResults(std::ostream &rOStream) const
    {
        std::vector<CellType *> leaves;

        mpQuadtree->GetAllLeavesVector(leaves);

        rOStream << "GiD Post Results File 1.0" << std::endl
                 << std::endl;

        rOStream << "Result \"Distance\" \"Kratos\" 1 Scalar OnNodes" << std::endl;

        rOStream << "Values" << std::endl;

        for (DistanceSpatialContainersConditionConfigure2d::data_type::const_iterator i_node = mQuadtreeNodes.begin(); i_node != mQuadtreeNodes.end(); i_node++)
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
    ModelPart &mrSkinModelPart;

    //ModelPart& mrBackgroundModelPart;
    ModelPart &mrFluidModelPart;
    DistanceSpatialContainersConditionConfigure2d::data_type mQuadtreeNodes;

    boost::shared_ptr<QuadtreeType> mpQuadtree;

    static const double epsilon;
    double dist_limit;
    //typename ParallelDistanceCalculator<2>::Pointer pDistanceCalculator;

    /**
         * @}
         */
    /**
         * calculates the eigenvectors and eigenvalues of given symmetric matrix A.
         * The eigenvectors and eigenvalues are calculated using the iterative
         * Gauss-Seidel-method
         * @param A the given symmetric matrix the eigenvectors are to be calculated.
         * :WARNING: Matrix A will be overwritten and has to be symmetric
         * @param V the result matrix (will be overwritten with the eigenvectors)
         * @param zero_tolerance the largest value considered to be zero
         */

    static inline void EigenVectors(const Matrix &A, Matrix &vectors, Vector &lambda, double zero_tolerance = 1e-9, int max_iterations = 10)
    {
        Matrix Help = A;

        for (int i = 0; i < 2; i++)
            for (int j = 0; j < 2; j++)
                Help(i, j) = Help(i, j);

        vectors.resize(Help.size1(), Help.size2(), false);

        lambda.resize(Help.size1(), false);

        Matrix HelpDummy(Help.size1(), Help.size2());

        bool is_converged = false;

        Matrix unity = ZeroMatrix(Help.size1(), Help.size2());

        for (std::size_t i = 0; i < Help.size1(); i++)
            unity(i, i) = 1.0;

        Matrix V = unity;

        Matrix VDummy(Help.size1(), Help.size2());

        Matrix Rotation(Help.size1(), Help.size2());

        for (int iterations = 0; iterations < max_iterations; iterations++)
        {

            is_converged = true;

            double a = 0.0;

            std::size_t index1 = 0;

            std::size_t index2 = 1;

            for (std::size_t i = 0; i < Help.size1(); i++)
            {
                for (std::size_t j = (i + 1); j < Help.size2(); j++)
                {
                    if ((fabs(Help(i, j)) > a) && (fabs(Help(i, j)) > zero_tolerance))
                    {
                        a = fabs(Help(i, j));

                        index1 = i;
                        index2 = j;

                        is_converged = false;
                    }
                }
            }

            //                 KRATOS_WATCH(Help);

            if (is_converged)
                break;

            //Calculation of Rotationangle

            double gamma = (Help(index2, index2) - Help(index1, index1)) / (2 * Help(index1, index2));

            double u = 1.0;

            if (fabs(gamma) > zero_tolerance && fabs(gamma) < (1 / zero_tolerance))
            {
                u = gamma / fabs(gamma) * 1.0 / (fabs(gamma) + sqrt(1.0 + gamma * gamma));
            }
            else
            {
                if (fabs(gamma) >= (1.0 / zero_tolerance))
                    u = 0.5 / gamma;
            }

            double c = 1.0 / (sqrt(1.0 + u * u));

            double s = c * u;

            double trianglea = s / (1.0 + c);

            //Ratotion of the Matrix
            HelpDummy = Help;

            HelpDummy(index2, index2) = Help(index2, index2) + u * Help(index1, index2);
            HelpDummy(index1, index1) = Help(index1, index1) - u * Help(index1, index2);
            HelpDummy(index1, index2) = 0.0;
            HelpDummy(index2, index1) = 0.0;

            for (std::size_t i = 0; i < Help.size1(); i++)
            {
                if ((i != index1) && (i != index2))
                {
                    HelpDummy(index2, i) = Help(index2, i) + s * (Help(index1, i) - trianglea * Help(index2, i));
                    HelpDummy(i, index2) = Help(index2, i) + s * (Help(index1, i) - trianglea * Help(index2, i));

                    HelpDummy(index1, i) = Help(index1, i) - s * (Help(index2, i) + trianglea * Help(index1, i));
                    HelpDummy(i, index1) = Help(index1, i) - s * (Help(index2, i) + trianglea * Help(index1, i));
                }
            }

            Help = HelpDummy;

            //Calculation of the eigenvectors V
            Rotation = unity;
            Rotation(index2, index1) = -s;
            Rotation(index1, index2) = s;
            Rotation(index1, index1) = c;
            Rotation(index2, index2) = c;

            //                 Help=ZeroMatrix(A.size1(),A.size1());

            VDummy = ZeroMatrix(Help.size1(), Help.size2());

            for (std::size_t i = 0; i < Help.size1(); i++)
            {
                for (std::size_t j = 0; j < Help.size1(); j++)
                {
                    for (std::size_t k = 0; k < Help.size1(); k++)
                    {
                        VDummy(i, j) += V(i, k) * Rotation(k, j);
                    }
                }
            }
            V = VDummy;
        }

        if (!(is_converged))
        {
            std::cout << "########################################################" << std::endl;
            std::cout << "Max_Iterations exceed in Jacobi-Seidel-Iteration (eigenvectors)" << std::endl;
            std::cout << "########################################################" << std::endl;
        }

        for (std::size_t i = 0; i < Help.size1(); i++)
        {
            for (std::size_t j = 0; j < Help.size1(); j++)
            {
                vectors(i, j) = V(j, i);
            }
        }

        for (std::size_t i = 0; i < Help.size1(); i++)
            lambda(i) = Help(i, i);

        return;
    }

    inline void CreatePartition(std::size_t number_of_threads, const int number_of_rows, vector<std::size_t> &partitions)
    {
        partitions.resize(number_of_threads + 1);
        int partition_size = number_of_rows / number_of_threads;
        partitions[0] = 0;
        partitions[number_of_threads] = number_of_rows;
        for (std::size_t i = 1; i < number_of_threads; i++)
            partitions[i] = partitions[i - 1] + partition_size;
    }

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
    CalculateSignedDistanceTo2DConditionSkinProcess &operator=(CalculateSignedDistanceTo2DConditionSkinProcess const &rOther);

    /// Copy constructor.
    //CalculateSignedDistanceTo2DConditionSkinProcess(CalculateSignedDistanceTo2DConditionSkinProcess const& rOther);

    ///@}

}; // Class CalculateSignedDistanceTo2DConditionSkinProcess

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream &operator>>(std::istream &rIStream,
                                CalculateSignedDistanceTo2DConditionSkinProcess &rThis);

/// output stream function
inline std::ostream &operator<<(std::ostream &rOStream,
                                const CalculateSignedDistanceTo2DConditionSkinProcess &rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

const double CalculateSignedDistanceTo2DConditionSkinProcess::epsilon = 1e-12;

} // namespace Kratos.

#endif // KRATOS_CALCULATE_DISTANCE_PROCESS_H_INCLUDED  defined
