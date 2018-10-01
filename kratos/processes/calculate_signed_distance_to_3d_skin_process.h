//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Daniel Baumgaertner
//                   Johannes Wolf
//


#if !defined(KRATOS_CALCULATE_DISTANCE_PROCESS_H_INCLUDED )
#define  KRATOS_CALCULATE_DISTANCE_PROCESS_H_INCLUDED



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

#include "spatial_containers/octree_binary.h"
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
#include "processes/calculate_distance_to_skin_process.h"

#ifdef _OPENMP
#include "omp.h"
#endif

using namespace boost::numeric::ublas;


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

    typedef Point                                           PointType;  /// always the point 3D
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

    typedef std::vector<PointerType>::iterator             PointerTypeIterator;




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
        return new data_type(27, (CellNodeData*)NULL);
    }

    static void CopyData(data_type* source, data_type* destination) {
        *destination = *source;
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


    static inline bool  IntersectionBox(const PointerType& rObject,  const PointType& rLowPoint, const PointType& rHighPoint)
    {
        return rObject->GetGeometry().HasIntersection(rLowPoint, rHighPoint);
    }


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
    typedef ConfigurationType::cell_node_data_type CellNodeDataType;
    typedef Point PointType;  /// always the point 3D
    typedef OctreeType::cell_type::object_container_type object_container_type;
    typedef struct{
        array_1d<double,3>  Coordinates;
        array_1d<double,3>  StructElemNormal;
        unsigned int EdgeNode1;
        unsigned int EdgeNode2;
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
    ~CalculateSignedDistanceTo3DSkinProcess() override
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
        KRATOS_TRY;

        GenerateOctree();

        //DistanceFluidStructure();

		CalculateDistanceToSkinProcess<3> distance_process(mrFluidModelPart, mrBodyModelPart);
        distance_process.Execute();

        //          ------------------------------------------------------------------
        //          GenerateNodes();
        CalculateDistance2(); // I have to change this. Pooyan.
        //mrSkinModelPart.GetCommunicator().AssembleCurrentData(DISTANCE);
        //          std::ofstream mesh_file1("octree1.post.msh");
        //          std::ofstream res_file("octree1.post.res");
        //          Timer::Start("Writing Gid conform Mesh");
        //          PrintGiDMesh(mesh_file1);
        //          PrintGiDResults(res_file);
        //          octree.PrintGiDMeshNew(mesh_file2);
        //          Timer::Stop("Writing Gid conform Mesh");
        //          delete octree. TODO: Carlos
        //          ------------------------------------------------------------------

        KRATOS_CATCH("");
    }

    ///******************************************************************************************************************
    ///******************************************************************************************************************
    /**
       * This function maps the nodal pressure values computed in the CFD analysis to the respective
       * structural nodes, i.e. for each structural node inside a fluid tetrahedra positive and negative
       * face pressure is computed by mapping between the nodal values of the tetrahedra. Afterwards
       * the resulting delta is applied as new nodal pressure.
       */
    void MappingPressureToStructure(BinBasedFastPointLocator<3>& node_locator)
    {
        //loop over nodes and find the tetra in which it falls, than do interpolation
        Vector N;
        const int max_results = 10000;
        BinBasedFastPointLocator<3>::ResultContainerType results(max_results);
        const int n_structure_nodes = mrSkinModelPart.Nodes().size();

#pragma omp parallel for firstprivate(results,N)
	//MY NEW LOOP: reset the viisted flaf
	for (int i = 0; i < n_structure_nodes; i++)
        {
            ModelPart::NodesContainerType::iterator iparticle = mrSkinModelPart.NodesBegin() + i;
            Node < 3 > ::Pointer p_structure_node = *(iparticle.base());
	    p_structure_node->Set(VISITED, false);
	}
        for (int i = 0; i < n_structure_nodes; i++)
        {
            ModelPart::NodesContainerType::iterator iparticle = mrSkinModelPart.NodesBegin() + i;
            Node < 3 > ::Pointer p_structure_node = *(iparticle.base());
            BinBasedFastPointLocator<3>::ResultIteratorType result_begin = results.begin();
            Element::Pointer pElement;

            bool is_found = node_locator.FindPointOnMesh(p_structure_node->Coordinates(), N, pElement, result_begin, max_results);

            if (is_found == true)
            {
                array_1d<double,4> nodalPressures;
                const Vector& ElementalDistances = pElement->GetValue(ELEMENTAL_DISTANCES);

                Geometry<Node<3> >& geom = pElement->GetGeometry();

                for(unsigned int j=0; j<geom.size(); j++)
                {
                    nodalPressures[j] = geom[j].FastGetSolutionStepValue(PRESSURE);
                }

                if(pElement->GetValue(SPLIT_ELEMENT)==true)
                {
                    array_1d<double,4> Npos,Nneg;

                    // Do mapping
                    ComputeDiscontinuousInterpolation((*p_structure_node),pElement->GetGeometry(),ElementalDistances,Npos,Nneg);

                    // Compute face pressure
                    double p_positive_structure = inner_prod(nodalPressures,Npos);
                    double p_negative_structure = inner_prod(nodalPressures,Nneg);

                    // Assign ModelPart::ElementIteratorface pressure to structure node
                    p_structure_node->FastGetSolutionStepValue(POSITIVE_FACE_PRESSURE) = p_positive_structure;
                    p_structure_node->FastGetSolutionStepValue(NEGATIVE_FACE_PRESSURE) = p_negative_structure;
		    p_structure_node->Set(VISITED);
                }
                else
                {
                    double p = inner_prod(nodalPressures,N);
                    p_structure_node->FastGetSolutionStepValue(POSITIVE_FACE_PRESSURE) = p;
                    p_structure_node->FastGetSolutionStepValue(NEGATIVE_FACE_PRESSURE) = p;
	            p_structure_node->Set(VISITED);
                }
            }
        }
	//AND NOW WE "TREAT" the bad nodes, the ones that belong to the structural faces that by some chance did not cross the fluid elements
	//to such nodes we simply extrapolate the pressure from the neighbors
	int n_bad_nodes=0;
	for (int i = 0; i < n_structure_nodes; i++)
        {
            ModelPart::NodesContainerType::iterator iparticle = mrSkinModelPart.NodesBegin() + i;
            Node < 3 > ::Pointer p_structure_node = *(iparticle.base());
	    if (p_structure_node->IsNot(VISITED))
		n_bad_nodes++;
	}
	//KRATOS_WATCH("THERE WERE THIS MANY BAD NODES ORIGINALLY")
	//KRATOS_WATCH(n_bad_nodes)
	while (n_bad_nodes >= 1.0) {
                int n_bad_nodes_backup = n_bad_nodes;

                for (int i = 0; i < n_structure_nodes; i++) {
                    ModelPart::NodesContainerType::iterator iparticle = mrSkinModelPart.NodesBegin() + i;
                    Node < 3 > ::Pointer p_structure_node = *(iparticle.base());

                    //here we store the number of neigbor nodes that were given the pressure in the previous loop (i.e. were found)
                    if (p_structure_node->IsNot(VISITED)) {
                        int n_good_neighbors = 0;
                        double pos_pres = 0.0;
                        double neg_pres = 0.0;
                        WeakPointerVector< Node < 3 > >& neighours = p_structure_node->GetValue(NEIGHBOUR_NODES);
                        
                        for (WeakPointerVector< Node < 3 > >::iterator j = neighours.begin(); j != neighours.end(); j++) {
                            if (j->Is(VISITED)) {
                                n_good_neighbors++;
                                pos_pres += j->FastGetSolutionStepValue(POSITIVE_FACE_PRESSURE);
                                neg_pres += j->FastGetSolutionStepValue(NEGATIVE_FACE_PRESSURE);
                                //KRATOS_WATCH("Good neighbor found")
                            }
                        }
                        if (n_good_neighbors != 0) {
                            pos_pres /= n_good_neighbors;
                            neg_pres /= n_good_neighbors;
                            p_structure_node->FastGetSolutionStepValue(POSITIVE_FACE_PRESSURE) = pos_pres;
                            p_structure_node->FastGetSolutionStepValue(NEGATIVE_FACE_PRESSURE) = neg_pres;
                            p_structure_node->Set(VISITED);
                            n_bad_nodes--;
                        }
                        //KRATOS_WATCH(pos_pres)
                        //KRATOS_WATCH(neg_pres)          
                    }
                }
                
                if(n_bad_nodes == n_bad_nodes_backup) break; //WE BREAK THE WHILE HERE, OTHERWISE THE CODE HANGS (it was not able to remove any other node)

                /*int n_bad_nodes=0;
                for (int i = 0; i < n_structure_nodes; i++)
                {
                   ModelPart::NodesContainerType::iterator iparticle = mrSkinModelPart.NodesBegin() + i;
                   Node < 3 > ::Pointer p_structure_node = *(iparticle.base());
                   if (p_structure_node->IsNot(VISITED))
                       n_bad_nodes++;		  
                }
                 */
                //KRATOS_WATCH(n_bad_nodes)

            }
		//THE BELOW ONE IS A "CHEAT".. THERE IS A PROBLEM OF INCORRECT PROJECTION BETWEEN THE MESHES AT SOME POINTS
		//FOR NODES WITH PRESSURE VERY DIFFERENT FROM THAT OF THE NEIGHBORS, I JUST TAKE THE NEIGHBOR PRESSURE AVERAGED
		for (int i = 0; i < n_structure_nodes; i++)
		{
		    ModelPart::NodesContainerType::iterator iparticle = mrSkinModelPart.NodesBegin() + i;
		    Node < 3 > ::Pointer p_structure_node = *(iparticle.base());

		    double pos_pressure=p_structure_node->FastGetSolutionStepValue(POSITIVE_FACE_PRESSURE);
		    double neg_pressure=p_structure_node->FastGetSolutionStepValue(NEGATIVE_FACE_PRESSURE);
                    
                    WeakPointerVector< Node < 3 > >& neighours = p_structure_node->GetValue(NEIGHBOUR_NODES);
			
	  	    if (neighours.size()>=1.0)
			{			    
			    double av_pos_pres=0.0;
			    double av_neg_pres=0.0;
			    for( WeakPointerVector< Node<3> >::iterator j = neighours.begin();
				        j != neighours.end(); j++)
			    {
				
					av_pos_pres+=j->FastGetSolutionStepValue(POSITIVE_FACE_PRESSURE);
					av_neg_pres+=j->FastGetSolutionStepValue(NEGATIVE_FACE_PRESSURE);
										
			    }
			    av_pos_pres/=neighours.size();
			    av_neg_pres/=neighours.size();
			    
			    //IF the average pressure of the neighbors is 10 times lower than of the given node, something is bad and we reset its value
			    if (fabs(pos_pressure)>3.0*fabs(av_pos_pres))
					{
					p_structure_node->FastGetSolutionStepValue(POSITIVE_FACE_PRESSURE) = av_pos_pres;
					//KRATOS_WATCH("BAD NODE")
					}
			    if (fabs(neg_pressure)>3.0*fabs(av_neg_pres))
					{
					p_structure_node->FastGetSolutionStepValue(NEGATIVE_FACE_PRESSURE) = av_neg_pres;
					//KRATOS_WATCH("BAD NODE")
					}				    	
			 
			}
		 }
		


    }

    ///******************************************************************************************************************
    ///******************************************************************************************************************

    void ComputeDiscontinuousInterpolation( const Node<3>& pNode,
                                            Geometry< Node<3> >& geom,
                                            const array_1d<double,4>& distances,
                                            array_1d<double,4>& Npos,
                                            array_1d<double,4>& Nneg)
    {
        //count positives
        int n_positives = 0;
        for(unsigned int i=0; i<distances.size(); i++)
            if(distances[i]>0) n_positives++;

        //generate the points on the edges at the zero of the distance function
        //generate "father nodes", defined as the end nodes of the edge on which the local point is located
        std::vector< Point > edge_points;
        edge_points.reserve(4);
        array_1d<unsigned int, 4> positive_fathers,  negative_fathers;	//there are at most 4 cut edges
        unsigned int k=0;
        unsigned int l=0;

        for(unsigned int i=0; i<3; i++)
        {
            for(unsigned int j=i+1; j<4; j++) // go through the edges 01, 02, 03, 12, 13, 23
            {
                double di = distances[i];
                double dj = distances[j];

                if(di*dj < 0) //edge is cut
                {
                    //generate point on edge by linear interpolation
                    double Ni = fabs(dj) / ( fabs(di) + fabs(dj) );
                    double Nj = 1.0 - Ni;
                    Point edge_point(Ni * geom[i] + Nj * geom[j]);
                    edge_points.push_back(edge_point);

                    //store the id of the positive and negative fathers
                    if(di > 0.0)
                    {
                        positive_fathers[k++] = i;
                        negative_fathers[l++] = j;
                    }
                    else
                    {
                        positive_fathers[k++] = j;
                        negative_fathers[l++] = i;
                    }
                }
            }
        }

        if(edge_points.size() == 3)
        {
            //compute local shape functions (tell how to interpolate from the edge nodes)
            Vector Nlocal(3);

            //form a triangle with the edge nodes
            Triangle3D3< Point > triangle(Point::Pointer(new Point(edge_points[0])), 
					     Point::Pointer(new Point(edge_points[1])), 
					     Point::Pointer(new Point(edge_points[2]))
					     );

            array_1d<double,3> local_coords;
            local_coords = triangle.PointLocalCoordinates(local_coords, pNode);

            for(unsigned int i=0; i<3;i++)
                Nlocal[i] = triangle.ShapeFunctionValue(i, local_coords );

            noalias(Npos) = ZeroVector(4);
            noalias(Nneg) = ZeroVector(4);
            for(unsigned int i=0; i<3; i++)
            {
                Npos[ positive_fathers[i] ] += Nlocal[i];
                Nneg[ negative_fathers[i] ] += Nlocal[i];
            }
        }

        if(edge_points.size() == 4)
        {
            //compute local shape functions (tell how to interpolate from the edge nodes)
            Vector Nlocal(4);

            //form a quadrilatera with the 4 cut nodes
            array_1d<double,3> x21 = edge_points[1] - edge_points[0];
            array_1d<double,3> x31 = edge_points[2] - edge_points[0];
            array_1d<double,3> x41 = edge_points[3] - edge_points[0];

            //define a vector oriented as x21
            array_1d<double,3> v1 = x21 / norm_2(x21);

            BoundedMatrix<double,4,3> DN_DX;
            array_1d<double,4> msN;
            double Area;
            GeometryUtils::CalculateGeometryData( geom, DN_DX, msN, Area );

            array_1d<double,3> n = prod(trans(DN_DX),distances);
            n /= norm_2(n);

            array_1d<double,3> v2;
            MathUtils<double>::CrossProduct(v2,v1,n); // v2 = v1 x n

            array_1d<double,3> angles;
            angles[0] = 0.0; //angle between x21 and v1
            angles[1] = atan2( inner_prod(x31,v2), inner_prod(x31,v1) ); //angle between x31 and v1
            angles[2] = atan2( inner_prod(x41,v2), inner_prod(x41,v1) ); //angle between x31 and v1

            double max_angle = 0.0;
            double min_angle = 0.0;
            unsigned int min_pos = 1;
            unsigned int max_pos = 1;
            for(unsigned int i=1; i<3; i++)
            {
                if(angles[i] < min_angle)
                {
                    min_pos = i+1; //this is the local index of the edge point which forms the minimal angle
                    min_angle = angles[i];
                }
                else if(angles[i] > max_angle)
                {
                    max_pos = i+1; //this is the local index of the edge point which forms the maximal angle
                    max_angle = angles[i];
                }
            }

            //find the pos of the center node
            unsigned int center_pos = 0;
            for(unsigned int i=1; i<4; i++)
            {
                if((i!= min_pos) && (i!=max_pos))
                { center_pos = i; }
            }

            //form a quadrilateral with the edge nodes
            Quadrilateral3D4< Point > quad = Quadrilateral3D4< Point >(
			Point::Pointer(new Point(edge_points[0])),
			Point::Pointer(new Point(edge_points[min_pos])),
			Point::Pointer(new Point(edge_points[center_pos])), 
			Point::Pointer(new Point(edge_points[max_pos]))
			);

            array_1d<double,3> local_coords;
            local_coords = quad.PointLocalCoordinates(local_coords, pNode);

            array_1d<unsigned int, 4> indices;
            indices[0] = 0;
            indices[1] = min_pos;
            indices[2] = center_pos;
            indices[3] = max_pos;

            for(unsigned int i=0; i<4;i++)
                Nlocal[ i ]  = quad.ShapeFunctionValue(i, local_coords );

            noalias(Npos) = ZeroVector(4);
            noalias(Nneg) = ZeroVector(4);
            for(unsigned int i=0; i<4; i++)
            {
                Npos[ positive_fathers[i] ] += Nlocal[indices[i]];
                Nneg[ negative_fathers[i] ] += Nlocal[indices[i]];
            }
        }
    }

    ///******************************************************************************************************************
    ///******************************************************************************************************************

    void AveragePressureToNode(BinBasedFastPointLocator<3>& node_locator,
                               Node<3>& node)
    {
        //loop over nodes and find the tetra in which it falls, than do interpolation
        Vector N;
        const int max_results = 10000;
        BinBasedFastPointLocator<3>::ResultContainerType results(max_results);
        BinBasedFastPointLocator<3>::ResultIteratorType result_begin = results.begin();
        Element::Pointer pElement;

        bool is_found = node_locator.FindPointOnMesh(node.Coordinates(), N, pElement, result_begin, max_results);

        if (is_found == true)
        {
            array_1d<double,4> nodalPressures;
            const Vector& ElementalDistances = pElement->GetValue(ELEMENTAL_DISTANCES);
            Geometry<Node<3> >& geom = pElement->GetGeometry();

            for(unsigned int i=0; i<4; i++)
                nodalPressures[i] = geom[i].GetSolutionStepValue(PRESSURE);

            if(pElement->GetValue(SPLIT_ELEMENT)==true)
            {
                // Compute average of all positive and all negative values
                double positiveAverage = 0;
                double negativeAverage = 0;
                unsigned int nPos = 0;
                unsigned int nNeg = 0;

                for(unsigned int i=0 ; i<4 ; i++)
                {
                    if(ElementalDistances[i]>=0)
                    {
                        positiveAverage += nodalPressures[i];
                        nPos++;
                    }
                    else
                    {
                        negativeAverage += nodalPressures[i];
                        nNeg++;
                    }
                }

                positiveAverage /= nPos;
                negativeAverage /= nNeg;

                // Assign Pressures
                node.GetSolutionStepValue(POSITIVE_FACE_PRESSURE,0) = positiveAverage;
                node.GetSolutionStepValue(NEGATIVE_FACE_PRESSURE,0) = negativeAverage;
            }
            else
            {
                // Compute average of all positive and all negative values
                double Average = 0;

                // for output of
                for(unsigned int i = 0 ; i<4 ; i++)
                    Average += nodalPressures[i];

                Average /= 4;

                // Assign Pressures
                node.GetSolutionStepValue(POSITIVE_FACE_PRESSURE,0) = Average;
                node.GetSolutionStepValue(NEGATIVE_FACE_PRESSURE,0) = Average;
            }
        }
    }

    ///******************************************************************************************************************
    ///******************************************************************************************************************

    void DistanceFluidStructure()
    {
      //std::cout << "Start calculating Elemental distances..." << std::endl;

        // Initialize Elemental distances in the domain
        Initialize();

        // Initialize index table that defines line Edges of fluid Element
        BoundedMatrix<unsigned int,6,2> TetEdgeIndexTable;
        SetIndexTable(TetEdgeIndexTable);

        // loop over all fluid Elements
        // this loop is parallelized using openmp
#ifdef _OPENMP
        int number_of_threads = omp_get_max_threads();
#else
        int number_of_threads = 1;
#endif

        ModelPart::ElementsContainerType& pElements = mrFluidModelPart.Elements();

        DenseVector<unsigned int> Element_partition;
        CreatePartition(number_of_threads, pElements.size(), Element_partition);

#pragma omp parallel for
        for (int k = 0; k < number_of_threads; k++)
        {
            ModelPart::ElementsContainerType::iterator it_begin = pElements.ptr_begin() + Element_partition[k];
            ModelPart::ElementsContainerType::iterator it_end = pElements.ptr_begin() + Element_partition[k+1];

            // assemble all Elements
            for (ModelPart::ElementIterator it = it_begin; it != it_end; ++it)
            {
                CalcElementDistances( it , TetEdgeIndexTable );
            }
        }        

        // Finally, each tetrahedral Element has 4 distance values. But each node belongs to
        // several Elements, such that it is assigned several distance values
        // --> now synchronize these values by finding the minimal distance and assign to each node a minimal nodal distance
        AssignMinimalNodalDistance();

        //std::cout << "Finished calculating Elemental distances..." << std::endl;
    }

    ///******************************************************************************************************************
    ///******************************************************************************************************************

    void Initialize()
    {
        const double initial_distance = 1.0;

        ModelPart::NodesContainerType::ContainerType& nodes = mrFluidModelPart.NodesArray();

        // reset the node distance to 1.0 which is the maximum distance in our normalized space.
        int nodesSize = nodes.size();

#pragma omp parallel for firstprivate(nodesSize)
        for(int i = 0 ; i < nodesSize ; i++)
            nodes[i]->GetSolutionStepValue(DISTANCE) = initial_distance;

        ModelPart::ElementsContainerType::ContainerType& fluid_Elements = mrFluidModelPart.ElementsArray();

        array_1d<double,4> ElementalDistances;
        ElementalDistances[0] = initial_distance;
        ElementalDistances[1] = initial_distance;
        ElementalDistances[2] = initial_distance;
        ElementalDistances[3] = initial_distance;

        // reset the Elemental distance to 1.0 which is the maximum distance in our normalized space.
        // also initialize the embedded velocity of the fluid Element
        int ElementsSize = fluid_Elements.size();

#pragma omp parallel for firstprivate(ElementsSize)
        for(int i = 0 ; i < ElementsSize ; i++)
        {
            fluid_Elements[i]->GetValue(ELEMENTAL_DISTANCES) = ElementalDistances;
            fluid_Elements[i]->GetValue(SPLIT_ELEMENT) = false;
            fluid_Elements[i]->GetValue(EMBEDDED_VELOCITY)=ZeroVector(3);
        }

    }

    ///******************************************************************************************************************
    ///******************************************************************************************************************

    void SetIndexTable( BoundedMatrix<unsigned int,6,2>& TetEdgeIndexTable )
    {
        // Initialize index table to define line Edges of fluid Element
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

    void CalcElementDistances( ModelPart::ElementsContainerType::iterator& i_fluidElement,
                               BoundedMatrix<unsigned int,6,2>            TetEdgeIndexTable )
    {
        std::vector<OctreeType::cell_type*> leaves;
        std::vector<TetEdgeStruct>          IntersectedTetEdges;
        unsigned int NumberIntersectionsOnTetCorner = 0;

        // Get leaves of octree intersecting with fluid Element
        mpOctree->GetIntersectedLeaves(*(i_fluidElement).base(),leaves);

        int intersection_counter = 0;

        // Loop over all 6 line Edges of the tetrahedra
        for(unsigned int i_tetEdge = 0;
            i_tetEdge < 6;
            i_tetEdge++)
        {
            IdentifyIntersectionNodes( i_fluidElement, i_tetEdge, leaves, IntersectedTetEdges, NumberIntersectionsOnTetCorner, TetEdgeIndexTable, intersection_counter );
        }

        if (intersection_counter!=0)
        {
            i_fluidElement->GetValue(EMBEDDED_VELOCITY) /= intersection_counter;
        }

        if(IntersectedTetEdges.size() > 0)
            CalcDistanceTo3DSkin( IntersectedTetEdges , i_fluidElement , NumberIntersectionsOnTetCorner );
    }

    ///******************************************************************************************************************
    ///******************************************************************************************************************

    void IdentifyIntersectionNodes( ModelPart::ElementsContainerType::iterator&   i_fluidElement,
                                    unsigned int                                  i_tetEdge,
                                    std::vector<OctreeType::cell_type*>&          leaves,
                                    std::vector<TetEdgeStruct>&                   IntersectedTetEdges,
                                    unsigned int&                                 NumberIntersectionsOnTetCorner,
                                    BoundedMatrix<unsigned int,6,2>              TetEdgeIndexTable,
                                    int&                                          intersection_counter )
    {
        std::vector<unsigned int> IntersectingStructElemID;
        TetEdgeStruct             NewTetEdge;
        unsigned int              NumberIntersectionsOnTetCornerCurrentEdge = 0;

        // Get nodes of line Edge
        unsigned int EdgeStartIndex = TetEdgeIndexTable(i_tetEdge,0);
        unsigned int EdgeEndIndex   = TetEdgeIndexTable(i_tetEdge,1);

        PointType& P1 = i_fluidElement->GetGeometry()[EdgeStartIndex];
        PointType& P2 = i_fluidElement->GetGeometry()[EdgeEndIndex];

        double EdgeNode1[3] = {P1.X() , P1.Y() , P1.Z()};
        double EdgeNode2[3] = {P2.X() , P2.Y() , P2.Z()};
	
        // loop over all octree cells which are intersected by the fluid Element
        for(unsigned int i_cell = 0 ; i_cell < leaves.size() ; i_cell++)
        {
            // Structural Element contained in one cell of the octree
            object_container_type* struct_elem = (leaves[i_cell]->pGetObjects());

            // loop over all structural Elements within each octree cell
            for(object_container_type::iterator i_StructElement = struct_elem->begin(); i_StructElement != struct_elem->end(); i_StructElement++)
            {

                if( StructuralElementNotYetConsidered( (*i_StructElement)->Id() , IntersectingStructElemID ) )
                {

                    // Calculate and associate intersection point to the current fluid Element
                    double IntersectionPoint[3] = {0.0 , 0.0 , 0.0};
                    int TetEdgeHasIntersections = IntersectionTriangleSegment( (*i_StructElement)->GetGeometry() , EdgeNode1 , EdgeNode2 , IntersectionPoint );

                    if( TetEdgeHasIntersections == 1 )
                    {
                        IntersectionNodeStruct NewIntersectionNode;

                        // Assign information to the intersection node
                        NewIntersectionNode.Coordinates[0] = IntersectionPoint[0];
                        NewIntersectionNode.Coordinates[1] = IntersectionPoint[1];
                        NewIntersectionNode.Coordinates[2] = IntersectionPoint[2];
			 
                        if( IsIntersectionNodeOnTetEdge( IntersectionPoint , EdgeNode1 , EdgeNode2 ) )
                        {
                            if ( IsNewIntersectionNode( NewIntersectionNode , IntersectedTetEdges ) )
                            {
                                // Calculate normal of the structural Element at the position of the intersection point
                                CalculateNormal3D((*i_StructElement)->GetGeometry()[0],
                                                  (*i_StructElement)->GetGeometry()[1],
                                                  (*i_StructElement)->GetGeometry()[2],
                                                  NewIntersectionNode.StructElemNormal);

                                // check, how many intersection nodes are located on corner points of the tetrahedra
                                if ( IsIntersectionOnCorner( NewIntersectionNode , EdgeNode1 , EdgeNode2) )
                                {
                                    NumberIntersectionsOnTetCornerCurrentEdge++;

                                    // only allow one intersection node on a tet edge
                                    if(NumberIntersectionsOnTetCornerCurrentEdge < 2)
                                    {
                                        // add the new intersection point to the list of intersection points of the fluid Element
                                        NewIntersectionNode.EdgeNode1 = EdgeStartIndex;
                                        NewIntersectionNode.EdgeNode2 = EdgeEndIndex;
                                        NewTetEdge.IntNodes.push_back(NewIntersectionNode);

                                        // if tet edge belonging to this intersection point is not already marked as "IntersectedTetEdge" --> put it into the respective container
                                        // when a second intersection node is found, then it is not necessary to push_back again
                                        if( NewTetEdge.IntNodes.size() == 1 )
                                            IntersectedTetEdges.push_back(NewTetEdge);
                                    }

                                    // this corner intersection node is only considered once for each tet edge
                                    if(NumberIntersectionsOnTetCornerCurrentEdge==1)
                                    {
                                        NumberIntersectionsOnTetCorner++;
                                    }
                                }
                                else
                                {
                                    // add the new intersection point to the list of intersection points of the fluid Element
                                    NewIntersectionNode.EdgeNode1 = EdgeStartIndex;
                                    NewIntersectionNode.EdgeNode2 = EdgeEndIndex;
                                    NewTetEdge.IntNodes.push_back(NewIntersectionNode);

                                    // velocity mapping structure --> fluid
                                    array_1d<double,3> emb_vel = (*i_StructElement)->GetGeometry()[0].GetSolutionStepValue(VELOCITY);
                                    emb_vel += (*i_StructElement)->GetGeometry()[1].GetSolutionStepValue(VELOCITY);
                                    emb_vel += (*i_StructElement)->GetGeometry()[2].GetSolutionStepValue(VELOCITY);

                                    i_fluidElement->GetValue(EMBEDDED_VELOCITY) += emb_vel/3;
                                    intersection_counter++;
                                }
                            }
                        }
                    }
                }
            }
        }

        // Finally put the found intersection nodes into the container
        if( NewTetEdge.IntNodes.size() > 0 )
	{
	  if(NumberIntersectionsOnTetCornerCurrentEdge == 0)
            IntersectedTetEdges.push_back(NewTetEdge);
	}
    }

    ///******************************************************************************************************************
    ///******************************************************************************************************************

    bool StructuralElementNotYetConsidered( unsigned int                IDCurrentStructElem,
                                            std::vector<unsigned int>&  IntersectingStructElemID )
    {
        // check if the structural Element was already considered as intersecting Element
        for(unsigned int k = 0 ; k < IntersectingStructElemID.size() ; k++)
        {
            if( IDCurrentStructElem == IntersectingStructElemID[k] )
                return false;
        }

        // if structural Element has not been considered in another octree, which also intersects the fluid Element
        // add the new object ID to the vector
        IntersectingStructElemID.push_back( IDCurrentStructElem );
        return true;
    }

    ///******************************************************************************************************************
    ///******************************************************************************************************************

    bool IsIntersectionNodeOnTetEdge( double* IntersectionPoint,
                                      double* EdgeNode1,
                                      double* EdgeNode2 )
    {
        // check, if intersection point is located on any edge of the fluid Element
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

    bool IsNewIntersectionNode( IntersectionNodeStruct&     NewIntersectionNode,
                                std::vector<TetEdgeStruct>& IntersectedTetEdges )
    {
        array_1d<double,3> DiffVector;
        double NormDiffVector = 0;
        unsigned int NumberIntNodes = 0;

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

    bool IsIntersectionOnCorner( IntersectionNodeStruct& NewIntersectionNode,
                                 double*                 EdgeNode1,
                                 double*                 EdgeNode2 )
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

    void CalculateNormal3D( Point&       Point1,
                            Point&       Point2,
                            Point&       Point3,
                            array_1d<double,3>&   rResultNormal )
    {
        array_1d<double,3> v1 = Point2 - Point1;
        array_1d<double,3> v2 = Point3 - Point1;

        MathUtils<double>::CrossProduct(rResultNormal,v1,v2);
        rResultNormal *= 0.5;
    }

    ///******************************************************************************************************************
    ///******************************************************************************************************************

    void CalcDistanceTo3DSkin( std::vector<TetEdgeStruct>&                 IntersectedTetEdges,
                               ModelPart::ElementsContainerType::iterator& i_fluid_Element,
                               unsigned int                                NumberIntersectionsOnTetCorner )
    {
        std::vector<IntersectionNodeStruct> NodesOfApproximatedStructure;
        array_1d<double,4> ElementalDistances;

        FillIntNodesContainer(IntersectedTetEdges,NodesOfApproximatedStructure);
	
        // Intersection with one corner point
        if( NodesOfApproximatedStructure.size() == 1 && NumberIntersectionsOnTetCorner == 1 )
        {
            CalcSignedDistancesToOneIntNode(i_fluid_Element,NodesOfApproximatedStructure,ElementalDistances);
            i_fluid_Element->GetValue(SPLIT_ELEMENT) = true;
        }

        // Intersection with two corner points / one tetrahedra edge
        if( NodesOfApproximatedStructure.size() == 2 && NumberIntersectionsOnTetCorner == 2 )
        {
            CalcSignedDistancesToTwoIntNodes(i_fluid_Element,NodesOfApproximatedStructure,ElementalDistances);
            i_fluid_Element->GetValue(SPLIT_ELEMENT) = true;
        }

        // Intersection with three tetrahedra edges
        if( NodesOfApproximatedStructure.size() == 3 )
        {
            CalcSignedDistancesToThreeIntNodes(i_fluid_Element,NodesOfApproximatedStructure,ElementalDistances);
            i_fluid_Element->GetValue(SPLIT_ELEMENT) = true;
        }

        // Intersection with more than three tetrahedra edges
        if( NodesOfApproximatedStructure.size() > 3 )
        {
            CalcSignedDistancesToMoreThanThreeIntNodes(i_fluid_Element,NodesOfApproximatedStructure,ElementalDistances,IntersectedTetEdges);
            i_fluid_Element->GetValue(SPLIT_ELEMENT) = true;
        }

        // Postprocessing treatment of Elemental distances
        if( i_fluid_Element->GetValue(SPLIT_ELEMENT) == true )
            AvoidZeroDistances(i_fluid_Element, ElementalDistances);

        // In case there is intersection with fluid Element: assign distances to the Element
        if( i_fluid_Element->GetValue(SPLIT_ELEMENT) == true )
            i_fluid_Element->GetValue(ELEMENTAL_DISTANCES) = ElementalDistances;
    }

    ///******************************************************************************************************************
    ///******************************************************************************************************************

    void FillIntNodesContainer( std::vector<TetEdgeStruct>&          IntersectedTetEdges,
                                std::vector<IntersectionNodeStruct>& NodesOfApproximatedStructure )
    {
        const unsigned int NumberCutEdges = IntersectedTetEdges.size();

        for(unsigned int i_TetEdge = 0 ; i_TetEdge < NumberCutEdges ; i_TetEdge++)
        {
            unsigned int NumberIntNodes = IntersectedTetEdges[i_TetEdge].IntNodes.size();

            for( unsigned int i_IntNode = 0 ; i_IntNode < NumberIntNodes ; i_IntNode++ )
            {
                NodesOfApproximatedStructure.push_back(IntersectedTetEdges[i_TetEdge].IntNodes[i_IntNode]);
            }
        }
    }

    ///******************************************************************************************************************
    ///******************************************************************************************************************

    void CalcSignedDistancesToOneIntNode( ModelPart::ElementsContainerType::iterator& i_fluid_Element,
                                          std::vector<IntersectionNodeStruct>         NodesOfApproximatedStructure,
                                          array_1d<double,4>&                         ElementalDistances )
    {
        Geometry< Node<3> >& rFluidGeom = i_fluid_Element->GetGeometry();

        Point  P1;
        P1.Coordinates() = NodesOfApproximatedStructure[0].Coordinates;

        array_1d<double,3>&  Normal = NodesOfApproximatedStructure[0].StructElemNormal;

        // Compute distance values for all tet-nodes
        for(unsigned int i_TetNode = 0 ; i_TetNode < 4 ; i_TetNode++)
        {
            ElementalDistances[i_TetNode] = PointDistanceToPlane(P1, Normal, rFluidGeom[i_TetNode]);
        }
    }

    ///******************************************************************************************************************
    ///******************************************************************************************************************

    void CalcSignedDistancesToTwoIntNodes( ModelPart::ElementsContainerType::iterator& i_fluid_Element,
                                           std::vector<IntersectionNodeStruct>         NodesOfApproximatedStructure,
                                           array_1d<double,4>&                         ElementalDistances )
    {
        Geometry< Node<3> >& rFluidGeom = i_fluid_Element->GetGeometry();

        Point  P1;
        P1.Coordinates() = NodesOfApproximatedStructure[0].Coordinates;

        // Get normal at intersections, average them and check direction of distances
        array_1d<double,3> NormalAtIntersectionNode1 = NodesOfApproximatedStructure[0].StructElemNormal;
        array_1d<double,3> NormalAtIntersectionNode2 = NodesOfApproximatedStructure[1].StructElemNormal;

        // Compute normal of surface plane
        array_1d<double,3> Normal;
        Normal[0] = 0.5*(NormalAtIntersectionNode1[0] + NormalAtIntersectionNode2[0]);
        Normal[1] = 0.5*(NormalAtIntersectionNode1[1] + NormalAtIntersectionNode2[1]);
        Normal[2] = 0.5*(NormalAtIntersectionNode1[2] + NormalAtIntersectionNode2[2]);

        // Check whether orientation of normal is in direction of the normal of the intersecting structure
        // Note: The normal of the approx. surface can be max. 90deg to every surrounding normal of the structure at the intersection nodes
        const array_1d<double,3> NormalAtOneIntersectionNode = NodesOfApproximatedStructure[0].StructElemNormal;

        bool NormalWrongOriented = false;
        if(inner_prod(NormalAtOneIntersectionNode,Normal)<0)
            NormalWrongOriented = true;

        // switch direction of normal
        if(NormalWrongOriented)
            Normal *=-1;	
	
        // Compute distance values for all tet-nodes
        for(unsigned int i_TetNode = 0 ; i_TetNode < 4 ; i_TetNode++)
        {
            ElementalDistances[i_TetNode] = PointDistanceToPlane(P1, Normal, rFluidGeom[i_TetNode]);
        }
    }

    ///******************************************************************************************************************
    ///******************************************************************************************************************

    void CalcSignedDistancesToThreeIntNodes( ModelPart::ElementsContainerType::iterator& i_fluid_Element,
                                             std::vector<IntersectionNodeStruct>&        NodesOfApproximatedStructure,
                                             array_1d<double,4>&                         ElementalDistances )
    {
        Geometry< Node<3> >& rFluidGeom = i_fluid_Element->GetGeometry();

        Point P1;
        Point P2;
        Point P3;

        P1.Coordinates() = NodesOfApproximatedStructure[0].Coordinates;
        P2.Coordinates() = NodesOfApproximatedStructure[1].Coordinates;
        P3.Coordinates() = NodesOfApproximatedStructure[2].Coordinates;

        array_1d<double,3> Normal;
        CalculateNormal3D(P1,P2,P3,Normal);

        // Check whether orientation of normal is in direction of the normal of the intersecting structure
        // Note: The normal of the approx. surface can be max. 90deg to every surrounding normal of the structure at the intersection nodes
        const array_1d<double,3> NormalAtOneIntersectionNode = NodesOfApproximatedStructure[0].StructElemNormal;

        bool NormalWrongOriented = false;
        if(inner_prod(NormalAtOneIntersectionNode,Normal)<0)
            NormalWrongOriented = true;

        // switch direction of normal
        if(NormalWrongOriented)
            Normal *=-1;

        // Compute distance values for all tet-nodes
        for(unsigned int i_TetNode = 0 ; i_TetNode < 4 ; i_TetNode++)
        {
            ElementalDistances[i_TetNode] = PointDistanceToPlane(P1, Normal, rFluidGeom[i_TetNode] );
        }
    }

    ///******************************************************************************************************************
    ///******************************************************************************************************************

    void CalcSignedDistancesToMoreThanThreeIntNodes(  ModelPart::ElementsContainerType::iterator& i_fluid_Element,
                                                      std::vector<IntersectionNodeStruct>         NodesOfApproximatedStructure,
                                                      array_1d<double,4>&                         ElementalDistances,
                                                      std::vector<TetEdgeStruct>&                 IntersectedTetEdges )
    {
        unsigned int numberCutEdges = NodesOfApproximatedStructure.size();

        // Compute average of the intersection nodes which is a node on the plane we look for
        Point P_mean;
        for(unsigned int k=0; k<numberCutEdges; k++)
            for(unsigned int i=0; i<3; i++)
                P_mean.Coordinates()[i] += NodesOfApproximatedStructure[k].Coordinates[i];

        for(unsigned int i=0; i<3; i++)
            P_mean.Coordinates()[i] /= numberCutEdges;

        // Compute normal for the best-fitted plane
        array_1d<double,3> N_mean;

        Matrix coordinates(numberCutEdges,3);
        for(unsigned int i=0; i<numberCutEdges; i++)
            for(unsigned int j=0; j<3; j++)
                coordinates(i,j) = NodesOfApproximatedStructure[i].Coordinates[j] - P_mean[j];

        Matrix A = prod(trans(coordinates),coordinates);
        Matrix V(3,3);
        Vector lambda(3);

        // Calculate the eigenvectors V and the corresponding eigenvalues lambda
        EigenVectors(A, V, lambda);

        // Look for the minimal eigenvalue all lambdas
        unsigned int min_pos = 0;
        double min_lambda = lambda[min_pos];
        for(unsigned int i=1;i<3; i++)
            if(min_lambda > lambda[i])
            {
                min_lambda = lambda[i];
                min_pos = i;
            }

        // the normal equals to the eigenvector which corresponds to the minimal eigenvalue
        for(unsigned int i=0;i<3; i++) N_mean[i] = V(min_pos,i);
        N_mean /= norm_2(N_mean);

        // Check whether orientation of normal is in direction of the normal of the intersecting structure
        // Note: The normal of the approx. surface can be max. 90deg to every surrounding normal of the structure at the intersection nodes
        array_1d<double,3> NormalAtOneIntersectionNode;
        NormalAtOneIntersectionNode = NodesOfApproximatedStructure[0].StructElemNormal;

        bool NormalWrongOriented = false;
        if(inner_prod(NormalAtOneIntersectionNode,N_mean)<0)
            NormalWrongOriented = true;

        // switch direction of normal
        if(NormalWrongOriented)
            N_mean *=-1;

        // Determine about the minimal distance by considering the distances to both triangles
        for(unsigned int i_TetNode = 0 ; i_TetNode < 4 ; i_TetNode++)
        {
            ElementalDistances[i_TetNode] = PointDistanceToPlane(P_mean, N_mean, i_fluid_Element->GetGeometry()[i_TetNode] );
        }

        // #################################################
        unsigned int numberDoubleCutEdges = 0;
        unsigned int indexDoubleCutEdge = 0;

        // figure out the edges which are cut more than once
        for(unsigned int i_TetEdge = 0 ; i_TetEdge < IntersectedTetEdges.size() ; i_TetEdge++)
        {
            unsigned int NumberIntNodes = IntersectedTetEdges[i_TetEdge].IntNodes.size();
            if(NumberIntNodes == 2)
            {
                numberDoubleCutEdges++;
                indexDoubleCutEdge = i_TetEdge;
            }
        }

        if((numberDoubleCutEdges >= 1))
        {
            array_1d<double,3> normal_1 = IntersectedTetEdges[indexDoubleCutEdge].IntNodes[0].StructElemNormal;
            array_1d<double,3> normal_2 = IntersectedTetEdges[indexDoubleCutEdge].IntNodes[1].StructElemNormal;

            // normalize normals
            normal_1 /= norm_2(normal_1);
            normal_2 /= norm_2(normal_2);

            const double pi = 3.1415926;

            // compute angle between normals
            double angle_n1n2 = acos( inner_prod(normal_1,normal_2) );
            // rad --> degree
            angle_n1n2 *= 180 / pi;

            // if angle between -60º and 120º, take the mean
            if( (angle_n1n2 > -60) && (angle_n1n2 < 120) )
            {
                // take the mean of the normals
                N_mean = 0.5 * (normal_1 + normal_2);
            }
            else
            {
                N_mean = 0.5 * (normal_1 - normal_2);
            }

            // Based on N_mean and P_mean compute the distances to that plane
            for(unsigned int i_TetNode = 0 ; i_TetNode < 4 ; i_TetNode++)
            {
                ElementalDistances[i_TetNode] = PointDistanceToPlane(P_mean, N_mean, i_fluid_Element->GetGeometry()[i_TetNode] );
            }
        }
    }

    ///******************************************************************************************************************
    ///******************************************************************************************************************

    /**
       * This function calculates the distance of a 3D point to a plane spanned by a 3D triangle
       * @param Plane base point
       * @param planeNormal
       * @param ToPoint The point which distance is required
       * @return The distance between the point and the plane spanned by the 3D triangle
       */
    double PointDistanceToPlane( Point&            planeBasePoint,
                                 array_1d<double, 3>& planeNormal,
                                 Point&            ToPoint)
    {
        // calculate vector pointing from a node in the plane (e.g. triangle point 1) to the considered node ToPoint
        array_1d<double,3> planeToPointVec = ToPoint - planeBasePoint;

        // projection of node on the plane
        const double sn = inner_prod(planeToPointVec,planeNormal);
        const double sd = inner_prod(planeNormal,planeNormal);
        double DistanceToPlane = sn / sqrt(sd);

        if( fabs(DistanceToPlane) < epsilon )
            DistanceToPlane = 0;

        return DistanceToPlane;
    }

    ///******************************************************************************************************************
    ///******************************************************************************************************************

    void AssignMinimalNodalDistance()
    {
        // loop over all fluid Elements
        for( ModelPart::ElementIterator i_fluid_Element = mrFluidModelPart.ElementsBegin();
             i_fluid_Element != mrFluidModelPart.ElementsEnd();
             i_fluid_Element++)
        {
            Geometry< Node<3> >& geom = i_fluid_Element->GetGeometry();
            const Vector& ElementalDistances = i_fluid_Element->GetValue(ELEMENTAL_DISTANCES);

            // Assign distances to the single nodes, if a smaller value is found
            for( unsigned int i_TetNode = 0; i_TetNode < 4; i_TetNode++ )
            {
                double currentNodeDist = geom[i_TetNode].GetSolutionStepValue(DISTANCE);
                double nodeInElemDist  = ElementalDistances[i_TetNode];

                if( fabs( nodeInElemDist ) < fabs( currentNodeDist ) )
                    geom[i_TetNode].GetSolutionStepValue(DISTANCE) = nodeInElemDist; // overwrite nodal distance (which is global)
            } // loop i_TetNode
        } // loop i_fluidElement
    }

    ///******************************************************************************************************************
    ///******************************************************************************************************************

    /**
       * If structure directly passes through the corner point of a tetrahedra (leading to zero distances
       * in the respective node), then a small distance value (different from zero) will be stored for
       * that point. This is necessary since the embedded solver cannot handle zero distances.
       * @param Element            current Element which was cut by the structure (flag SPLIT_ELEMENT is set to one)
       * @param ElementalDistances Elemental distances calculated by the intersection pattern
       */
    void AvoidZeroDistances( ModelPart::ElementsContainerType::iterator& Element,
                             array_1d<double,4>&                         ElementalDistances)
    {
        // Assign a distance limit
        double dist_limit = 1e-5;
//         bool distChangedToLimit = false; //variable to indicate that a distance value < tolerance is set to a limit distance = tolerance
// 
//         for(unsigned int i_node = 0; i_node < 4; i_node++)
//         {
//             if(fabs(ElementalDistances[i_node]) < dist_limit)
//             {
//                 ElementalDistances[i_node] = dist_limit;
//                 distChangedToLimit = true;
//             }
//         }
// 
//         // Check, if this approach changes the split-flag (might be, that Element is not cut anymore if node with zero distance gets a positive limit distance value
//         unsigned int numberNodesPositiveDistance = 0;
//         for(unsigned int i_node = 0; i_node < 4; i_node++)
//         {
//             if((ElementalDistances[i_node]) > 0)
//                 numberNodesPositiveDistance++;
//         }
        
        for(unsigned int i_node = 0; i_node < 4; i_node++)
        {
            double & di = ElementalDistances[i_node];
            if(fabs(di) < dist_limit)
            {
                if(di >= 0) di = dist_limit;
                else di = -dist_limit;
            }
        }

        // Element is not set
//         if(numberNodesPositiveDistance == 4 && distChangedToLimit == true)
//             Element->GetValue(SPLIT_ELEMENT) = false;
    }

    ///******************************************************************************************************************
    ///******************************************************************************************************************

    void GenerateSkinModelPart( ModelPart& mrNewSkinModelPart )
    {
        unsigned int id_node = mrFluidModelPart.NumberOfNodes() + 1;
        unsigned int id_condition = mrFluidModelPart.NumberOfConditions() + 1;

        mrNewSkinModelPart.Nodes().reserve(mrFluidModelPart.Nodes().size());
        mrNewSkinModelPart.Conditions().reserve(mrFluidModelPart.Elements().size());

        for(ModelPart::ElementIterator i_fluid_element = mrFluidModelPart.ElementsBegin();
            i_fluid_element != mrFluidModelPart.ElementsEnd();
            i_fluid_element++)
        {
            bool is_split = i_fluid_element->Is(TO_SPLIT);
            if(is_split == true)
            {
                const Vector& distances = i_fluid_element->GetValue(ELEMENTAL_DISTANCES);
                Geometry< Node<3> >& geom = i_fluid_element->GetGeometry();

                // generate the points on the edges at the zero of the distance function
                std::vector< Point > edge_points;
                edge_points.reserve(4);

                // loop over all 6 edges of the tetrahedra
                for(unsigned int i=0; i<3; i++)
                {
                    for(unsigned int j=i+1; j<4; j++) // go through the edges 01, 02, 03, 12, 13, 23
                    {
                        double di = distances[i];
                        double dj = distances[j];

                        if(di*dj < 0) //edge is cut
                        {
                            // generate point on edge by linear interpolation
                            double Ni = fabs(dj) / ( fabs(di) + fabs(dj) );
                            double Nj = 1.0 - Ni;
                            Point edge_point(Ni * geom[i] + Nj * geom[j]);
                            edge_points.push_back(edge_point);
                        }
                    }
                }

                // three intersection nodes
                if(edge_points.size() == 3)
                {
                    // ######## ADDING NEW NODE #########
                    Node < 3 >::Pointer pnode1 = mrNewSkinModelPart.CreateNewNode(id_node++,edge_points[0].X(),edge_points[0].Y(),edge_points[0].Z());
                    Node < 3 >::Pointer pnode2 = mrNewSkinModelPart.CreateNewNode(id_node++,edge_points[1].X(),edge_points[1].Y(),edge_points[1].Z());
                    Node < 3 >::Pointer pnode3 = mrNewSkinModelPart.CreateNewNode(id_node++,edge_points[2].X(),edge_points[2].Y(),edge_points[2].Z());

                    // ######## ADDING NEW CONDITION #########
                    //form a triangle
                    Triangle3D3< Node<3> > triangle(pnode1, pnode2, pnode3);

                    Condition const& rReferenceCondition = KratosComponents<Condition>::Get("Condition3D");
                    Properties::Pointer properties = mrNewSkinModelPart.rProperties()(0);
                    Condition::Pointer p_condition = rReferenceCondition.Create(id_condition++, triangle, properties);

                    mrNewSkinModelPart.Conditions().push_back(p_condition);
                }

                // four intersection nodes
                if(edge_points.size() == 4)
                {
                    //form a quadrilatera with the 4 cut nodes
                    array_1d<double,3> x21 = edge_points[1] - edge_points[0];
                    array_1d<double,3> x31 = edge_points[2] - edge_points[0];
                    array_1d<double,3> x41 = edge_points[3] - edge_points[0];

                    //define a vector oriented as x21
                    array_1d<double,3> v1 = x21 / norm_2(x21);

                    BoundedMatrix<double,4,3> DN_DX;
                    array_1d<double,4> msN;
                    double Area;
                    GeometryUtils::CalculateGeometryData( geom, DN_DX, msN, Area );

                    array_1d<double,3> n = prod(trans(DN_DX),distances);
                    n /= norm_2(n);

                    array_1d<double,3> v2;
                    MathUtils<double>::CrossProduct(v2,v1,n); // v2 = v1 x n

                    array_1d<double,3> angles;
                    angles[0] = 0.0; //angle between x21 and v1
                    angles[1] = atan2( inner_prod(x31,v2), inner_prod(x31,v1) ); //angle between x31 and v1
                    angles[2] = atan2( inner_prod(x41,v2), inner_prod(x41,v1) ); //angle between x31 and v1

                    double max_angle = 0.0;
                    double min_angle = 0.0;
                    unsigned int min_pos = 1;
                    unsigned int max_pos = 1;
                    for(unsigned int i=1; i<3; i++)
                    {
                        if(angles[i] < min_angle)
                        {
                            min_pos = i+1; //this is the local index of the edge point which forms the minimal angle
                            min_angle = angles[i];
                        }
                        else if(angles[i] > max_angle)
                        {
                            max_pos = i+1; //this is the local index of the edge point which forms the maximal angle
                            max_angle = angles[i];
                        }
                    }

                    //find the pos of the center node
                    unsigned int center_pos = 0;
                    for(unsigned int i=1; i<4; i++)
                    {
                        if((i!= min_pos) && (i!=max_pos))
                        { center_pos = i; }
                    }

                    // ######## ADDING NEW NODE #########
                    Node < 3 >::Pointer pnode1 = mrNewSkinModelPart.CreateNewNode(id_node++,edge_points[0].X(),edge_points[0].Y(),edge_points[0].Z());
                    Node < 3 >::Pointer pnode2 = mrNewSkinModelPart.CreateNewNode(id_node++,edge_points[min_pos].X(),edge_points[min_pos].Y(),edge_points[min_pos].Z());
                    Node < 3 >::Pointer pnode3 = mrNewSkinModelPart.CreateNewNode(id_node++,edge_points[center_pos].X(),edge_points[center_pos].Y(),edge_points[center_pos].Z());
                    Node < 3 >::Pointer pnode4 = mrNewSkinModelPart.CreateNewNode(id_node++,edge_points[max_pos].X(),edge_points[max_pos].Y(),edge_points[max_pos].Z());

                    // ######## ADDING NEW CONDITION #########
                    //form two triangles
                    Triangle3D3< Node<3> > triangle1(pnode1, pnode2, pnode3);
                    Triangle3D3< Node<3> > triangle2(pnode1, pnode3, pnode4);

                    Condition const& rReferenceCondition = KratosComponents<Condition>::Get("Condition3D");
                 
                    Properties::Pointer properties = mrNewSkinModelPart.rProperties()(0);

                    Condition::Pointer p_condition1 = rReferenceCondition.Create(id_condition++, triangle1, properties);
                    Condition::Pointer p_condition2 = rReferenceCondition.Create(id_condition++, triangle2, properties);

                    mrNewSkinModelPart.Conditions().push_back(p_condition1);
                    mrNewSkinModelPart.Conditions().push_back(p_condition2);
                    
                }

            }
        }

    }

    ///******************************************************************************************************************
    ///******************************************************************************************************************

    void GenerateOctree()
    {
        Timer::Start("Generating Octree");
        //std::cout << "Generating the Octree..." << std::endl;
        auto temp_octree =  Kratos::make_shared<OctreeType>();
        //OctreeType::Pointer temp_octree = OctreeType::Pointer(new OctreeType() );
        mpOctree.swap(temp_octree);
        
        double low[3];
        double high[3];
        
        for (int i = 0 ; i < 3; i++)
        {
            low[i] = high[i] = mrFluidModelPart.NodesBegin()->Coordinates()[i];
        }
        
        // loop over all nodes in the bounding box
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
        
        // loop over all skin nodes
        for(ModelPart::NodeIterator i_node = mrSkinModelPart.NodesBegin();
            i_node != mrSkinModelPart.NodesEnd();
            i_node++)
        {
            const array_1d<double,3>& r_coordinates = i_node->Coordinates();
            for (int i = 0 ; i < 3; i++)
            {
                low[i]  = r_coordinates[i] < low[i]  ? r_coordinates[i] : low[i];
                high[i] = r_coordinates[i] > high[i] ? r_coordinates[i] : high[i];
            }
        }
                
        mpOctree->SetBoundingBox(low,high);

        //mpOctree->RefineWithUniformSize(0.0625);

        // loop over all structure nodes
        for(ModelPart::NodeIterator i_node = mrSkinModelPart.NodesBegin();
            i_node != mrSkinModelPart.NodesEnd();
            i_node++)
        {
            double temp_point[3];
            temp_point[0] = i_node->X();
            temp_point[1] = i_node->Y();
            temp_point[2] = i_node->Z();
            mpOctree->Insert(temp_point);
        }

        //mpOctree->Constrain2To1(); // To be removed. Pooyan.

        // loop over all structure elements
        for(ModelPart::ElementIterator i_element = mrSkinModelPart.ElementsBegin();
            i_element != mrSkinModelPart.ElementsEnd();
            i_element++)
        {
            mpOctree->Insert(*(i_element).base());
        }

        Timer::Stop("Generating Octree");

//        KRATOS_WATCH(mpOctree);

//        std::cout << "######## WRITING OCTREE MESH #########" << std::endl;
//        std::ofstream myfile;
//        myfile.open ("octree.post.msh");
//        mpOctree.PrintGiDMesh(myfile);
//        myfile.close();

        //std::cout << "Generating the Octree finished" << std::endl;
    }

    ///******************************************************************************************************************
    ///******************************************************************************************************************

    void GenerateNodes()
    {
        Timer::Start("Generating Nodes");
        std::vector<OctreeType::cell_type*> all_leaves;
        mpOctree->GetAllLeavesVector(all_leaves);

        int leaves_size = all_leaves.size();

#pragma omp parallel for
        for (int i = 0; i < leaves_size; i++)
        {
            *(all_leaves[i]->pGetDataPointer()) = ConfigurationType::AllocateData();
        }


        std::size_t last_id = mrBodyModelPart.NumberOfNodes() + 1;

        for (std::size_t i = 0; i < all_leaves.size(); i++)
        {
                    CellType* cell = all_leaves[i];
            GenerateCellNode(cell, last_id);
        }

        Timer::Stop("Generating Nodes");

    }

    ///******************************************************************************************************************
    ///******************************************************************************************************************

    void GenerateCellNode(CellType* pCell, std::size_t& LastId)
    {
        for (int i_pos=0; i_pos < 8; i_pos++) // position 8 is for center
        {
            DistanceSpatialContainersConfigure::cell_node_data_type* p_node = (*(pCell->pGetData()))[i_pos];
            if(p_node == 0)
            {
                (*(pCell->pGetData()))[i_pos] = new DistanceSpatialContainersConfigure::cell_node_data_type;

                (*(pCell->pGetData()))[i_pos]->Id() = LastId++;

                        mOctreeNodes.push_back((*(pCell->pGetData()))[i_pos]);

                SetNodeInNeighbours(pCell,i_pos,(*(pCell->pGetData()))[i_pos]);
            }

        }
    }

    ///******************************************************************************************************************
    ///******************************************************************************************************************

    void SetNodeInNeighbours(CellType* pCell, int Position, CellNodeDataType* pNode)
    {
        CellType::key_type point_key[3];
        pCell->GetKey(Position, point_key);

        for (std::size_t i_direction = 0; i_direction < 8; i_direction++) {
            CellType::key_type neighbour_key[3];
            if (pCell->GetNeighbourKey(Position, i_direction, neighbour_key)) {
                CellType* neighbour_cell = mpOctree->pGetCell(neighbour_key);
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

    ///******************************************************************************************************************
    ///******************************************************************************************************************

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

        mpOctree->GetAllLeavesVector(leaves);
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

    //         mpOctree->GetAllLeavesVector(leaves);
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

    //         mpOctree->GetAllLeavesVector(leaves);
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
        DistanceSpatialContainersConfigure::data_type& nodes = mOctreeNodes;
        int nodes_size = nodes.size();
        // first of all we reste the node distance to 1.00 which is the maximum distnace in our normalized space.
#pragma omp parallel for firstprivate(nodes_size)
        for(int i = 0 ; i < nodes_size ; i++)
            nodes[i]->Distance() = 1.00;


        std::vector<CellType*> leaves;

        mpOctree->GetAllLeavesVector(leaves);
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
        DistanceSpatialContainersConfigure::data_type nodes_array;


        const double epsilon = 1e-12;

        double distance = 1.0;

        // Creating the ray
        double ray[3] = {coords[0], coords[1], coords[2]};
        
        mpOctree->NormalizeCoordinates(ray);
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
                mpOctree->CalculateCoordinates(keys,cell_point);

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
        //if(fabs(node_distance) > fabs(distance))
        //    node_distance = distance;
        /*else*/ if (distance*node_distance < 0.00) // assigning the correct sign
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
            
            mpOctree->NormalizeCoordinates(ray);
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

//#ifdef _DEBUG
//        std::cout << "colors : " << colors[0] << ", " << colors[1] << ", " << colors[2] << std::endl;
//#endif
        double distance = (fabs(distances[0]) > fabs(distances[1])) ? distances[1] : distances[0];
        distance = (fabs(distance) > fabs(distances[2])) ? distances[2] : distance;

        return distance;

    }

    void GetIntersectionsAndNodes(double* ray, int direction, std::vector<std::pair<double,Element::GeometryType*> >& intersections, DistanceSpatialContainersConfigure::data_type& rNodesArray)
    {
        //This function passes the ray through the model and gives the hit point to all objects in its way
        //ray is of dimension (3) normalized in (0,1)^3 space
        // direction can be 0,1,2 which are x,y and z respectively

        const double epsilon = 1.00e-12;

        // first clearing the intersections points vector
        intersections.clear();

        //OctreeType* octree = &mOctree;
        OctreeType* octree = mpOctree.get();

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
//#ifdef _DEBUG
//                Octree_Pooyan::key_type min_key[3];
//                cell->GetMinKey(min_key[0],min_key[1],min_key[2]);
//                Octree_Pooyan::key_type tmp;
//                tmp= min_key[direction];
//                assert(ray_key[direction]==tmp);
//#endif
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

        //OctreeType* octree = &mOctree;
        OctreeType* octree = mpOctree.get();

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
//#ifdef _DEBUG
//                Octree_Pooyan::key_type min_key[3];
//                cell->GetMinKey(min_key[0],min_key[1],min_key[2]);
//                Octree_Pooyan::key_type tmp;
//                tmp= min_key[direction];
//                assert(ray_key[direction]==tmp);
//#endif
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
        mpOctree->CalculateCoordinateNormalized(ray_key[direction], normalized_coordinate);
        ray_point1[direction] = normalized_coordinate;
        ray_point2[direction] = ray_point1[direction] + mpOctree->CalcSizeNormalized(cell);

        mpOctree->ScaleBackToOriginalCoordinate(ray_point1);
        mpOctree->ScaleBackToOriginalCoordinate(ray_point2);

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

		double triangle_origin_distance = -inner_prod(n, rGeometry[0]);
		Point ray_point_1, ray_point_2;
		
		for(int i = 0 ; i < 3 ; i++)
        {
            dir[i] = RayPoint2[i] - RayPoint1[i];             // ray direction vector
            w0[i] = RayPoint1[i] - rGeometry[0][i];
			ray_point_1[i] = RayPoint1[i];
			ray_point_2[i] = RayPoint2[i];
		}

		double sign_distance_1 = inner_prod(n, ray_point_1) + triangle_origin_distance;
		double sign_distance_2 = inner_prod(n, ray_point_2) + triangle_origin_distance;

		if (sign_distance_1*sign_distance_2 > epsilon) // segment line point on the same side of plane
			return 0;
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
        return "CalculateSignedDistanceTo3DSkinProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "CalculateSignedDistanceTo3DSkinProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

    void PrintGiDMesh(std::ostream & rOStream) const {
        std::vector<CellType*> leaves;

        mpOctree->GetAllLeavesVector(leaves);

        std::cout << "writing " << leaves.size() << " leaves" << std::endl;
        rOStream << "MESH \"leaves\" dimension 3 ElemType Hexahedra Nnode 8" << std::endl;
        rOStream << "# color 96 96 96" << std::endl;
        rOStream << "Coordinates" << std::endl;
        rOStream << "# node number coordinate_x coordinate_y coordinate_z  " << std::endl;

        for(DistanceSpatialContainersConfigure::data_type::const_iterator i_node = mOctreeNodes.begin() ; i_node != mOctreeNodes.end() ; i_node++)
        {
            rOStream << (*i_node)->Id() << "  " << (*i_node)->X() << "  " << (*i_node)->Y() << "  " << (*i_node)->Z() << std::endl;
            //mpOctree->Insert(temp_point);
        }
        std::cout << "Nodes written..." << std::endl;
        rOStream << "end coordinates" << std::endl;
        rOStream << "Elements" << std::endl;
        rOStream << "# Element node_1 node_2 node_3 material_number" << std::endl;

        for (std::size_t i = 0; i < leaves.size(); i++) {
            if ((leaves[i]->pGetData()))
            {
                DistanceSpatialContainersConfigure::data_type& nodes = (*(leaves[i]->pGetData()));

                rOStream << i + 1;
                for(int j = 0 ; j < 8 ; j++)
                    rOStream << "  " << nodes[j]->Id();
                rOStream << std::endl;
            }
        }
        rOStream << "end Elements" << std::endl;

    }

    void PrintGiDResults(std::ostream & rOStream) const {
        std::vector<CellType*> leaves;

        mpOctree->GetAllLeavesVector(leaves);

        rOStream << "GiD Post Results File 1.0" << std::endl << std::endl;

        rOStream << "Result \"Distance\" \"Kratos\" 1 Scalar OnNodes" << std::endl;

        rOStream << "Values" << std::endl;

        for(DistanceSpatialContainersConfigure::data_type::const_iterator i_node = mOctreeNodes.begin() ; i_node != mOctreeNodes.end() ; i_node++)
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

    DistanceSpatialContainersConfigure::data_type mOctreeNodes;

    Kratos::shared_ptr<OctreeType> mpOctree;

    static const double epsilon;

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

    static inline void EigenVectors(const Matrix& A, Matrix& vectors, Vector& lambda, double zero_tolerance =1e-9, int max_iterations = 10)
    {
        Matrix Help= A;

        for(int i=0; i<3; i++)
            for(int j=0; j<3; j++)
                Help(i,j)= Help(i,j);


        vectors.resize(Help.size1(),Help.size2(),false);

        lambda.resize(Help.size1(),false);

        Matrix HelpDummy(Help.size1(),Help.size2());

        bool is_converged = false;

        Matrix unity=ZeroMatrix(Help.size1(),Help.size2());

        for(unsigned int i=0; i< Help.size1(); i++)
            unity(i,i)= 1.0;

        Matrix V= unity;

        Matrix VDummy(Help.size1(),Help.size2());

        Matrix Rotation(Help.size1(),Help.size2());


        for(int iterations=0; iterations<max_iterations; iterations++)
        {

            is_converged= true;

            double a= 0.0;

            unsigned int index1= 0;

            unsigned int index2= 1;

            for(unsigned int i=0; i< Help.size1(); i++)
            {
                for(unsigned int j=(i+1); j< Help.size2(); j++)
                {
                    if((fabs(Help(i,j)) > a ) && (fabs(Help(i,j)) > zero_tolerance))
                    {
                        a= fabs(Help(i,j));

                        index1= i;
                        index2= j;

                        is_converged= false;
                    }
                }
            }

            //                 KRATOS_WATCH(Help);

            if(is_converged)
                break;

            //Calculation of Rotationangle

            double gamma= (Help(index2,index2)-Help(index1,index1))/(2*Help(index1,index2));

            double u=1.0;

            if(fabs(gamma) > zero_tolerance && fabs(gamma)< (1/zero_tolerance))
            {
                u= gamma/fabs(gamma)*1.0/(fabs(gamma)+sqrt(1.0+gamma*gamma));
            }
            else
            {
                if  (fabs(gamma)>= (1.0/zero_tolerance))
                    u= 0.5/gamma;
            }

            double c= 1.0/(sqrt(1.0+u*u));

            double s= c*u;

            double teta= s/(1.0+c);

            //Ratotion of the Matrix
            HelpDummy= Help;

            HelpDummy(index2,index2)= Help(index2,index2)+u*Help(index1,index2);
            HelpDummy(index1,index1)= Help(index1,index1)-u*Help(index1,index2);
            HelpDummy(index1,index2)= 0.0;
            HelpDummy(index2,index1)= 0.0;

            for(unsigned int i=0; i<Help.size1(); i++)
            {
                if((i!= index1) && (i!= index2))
                {
                    HelpDummy(index2,i)=Help(index2,i)+s*(Help(index1,i)- teta*Help(index2,i));
                    HelpDummy(i,index2)=Help(index2,i)+s*(Help(index1,i)- teta*Help(index2,i));

                    HelpDummy(index1,i)=Help(index1,i)-s*(Help(index2,i)+ teta*Help(index1,i));
                    HelpDummy(i,index1)=Help(index1,i)-s*(Help(index2,i)+ teta*Help(index1,i));
                }
            }


            Help= HelpDummy;

            //Calculation of the eigenvectors V
            Rotation =unity;
            Rotation(index2,index1)=-s;
            Rotation(index1,index2)=s;
            Rotation(index1,index1)=c;
            Rotation(index2,index2)=c;

            //                 Help=ZeroMatrix(A.size1(),A.size1());

            VDummy = ZeroMatrix(Help.size1(), Help.size2());

            for(unsigned int i=0; i< Help.size1(); i++)
            {
                for(unsigned int j=0; j< Help.size1(); j++)
                {
                    for(unsigned int k=0; k< Help.size1(); k++)
                    {
                        VDummy(i,j) += V(i,k)*Rotation(k,j);
                    }
                }
            }
            V= VDummy;
        }

        if(!(is_converged))
        {
            std::cout<<"########################################################"<<std::endl;
            std::cout<<"Max_Iterations exceed in Jacobi-Seidel-Iteration (eigenvectors)"<<std::endl;
            std::cout<<"########################################################"<<std::endl;
        }

        for(unsigned int i=0; i< Help.size1(); i++)
        {
            for(unsigned int j=0; j< Help.size1(); j++)
            {
                vectors(i,j)= V(j,i);
            }
        }

        for(unsigned int i=0; i<Help.size1(); i++)
            lambda(i)= Help(i,i);

        return;
    }


    inline void CreatePartition(unsigned int number_of_threads, const int number_of_rows, DenseVector<unsigned int>& partitions)
    {
        partitions.resize(number_of_threads + 1);
        int partition_size = number_of_rows / number_of_threads;
        partitions[0] = 0;
        partitions[number_of_threads] = number_of_rows;
        for (unsigned int i = 1; i < number_of_threads; i++)
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

const double CalculateSignedDistanceTo3DSkinProcess::epsilon = 1e-18;


}  // namespace Kratos.

#endif // KRATOS_CALCULATE_DISTANCE_PROCESS_H_INCLUDED  defined 


