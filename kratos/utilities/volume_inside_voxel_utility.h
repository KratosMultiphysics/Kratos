//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Ariadna Cortés
//
//

#if !defined(KRATOS_VOLUME_INSIDE_VOXEL)
#define  KRATOS_VOLUME_INSIDE_VOXEL

// System includes

// External includes

// Project includes
#include "includes/geometrical_object.h"
#include "includes/node.h"
#include "geometries/geometry.h"
#include "geometries/hexahedra_3d_8.h"
#include "intersection_utilities.h"

namespace Kratos
{ 
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

/**
 * @class VolumeInsideVoxel
 * @ingroup KratosCore
 * @brief Utilities to compute the real volume inside a voxel
 * @details This class provides static methods to compute (using different approximations) the portion of a 
 * voxel that is actually filled with volume.
 * @author Ariadna Cortés
 */
class VolumeInsideVoxelUtility
{
public:

    ///@name Type Definitions
    ///@{

    typedef Node<3> NodeType;
    typedef Node<3>::Pointer NodePtrType;
    typedef Geometry<NodeType> GeometryType;
    typedef GeometryType::Pointer GeometryPtrType;
    typedef GeometryType::GeometriesArrayType GeometryArrayType;
    typedef GeometryType::PointsArrayType PointsArrayType;


    /// Pointer definition of VoxelInsideVolume
    KRATOS_CLASS_POINTER_DEFINITION( VolumeInsideVoxelUtility );

    ///@}
    ///@name Life Cycle
    ///@{

    ///@}
    ///@name Operators
    ///@{

    /**
     * @brief Default constructor
     */
    VolumeInsideVoxelUtility(){}

    /// Destructor
    virtual ~VolumeInsideVoxelUtility(){}

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Aproximates the actual volume inside the voxel 
     * @param rVoxel references to the voxel whose actual volume will be approximated
     * @return Approximated volume 
     * @note This approximation assigns a fraction of volume (1/8) to each node of the
     * voxel, and counts it as volume if the node is inside the object (NodeDistance > 0)
     */  
    template<class TGeometryType>
    static double NodesApproximation(
        const TGeometryType& rVoxel        
    ) {
        double volume = 0;
        PointsArrayType nodes = rVoxel.Points();
        for (int i = 0; i < nodes.size(); i++) {
            //std::cout << "Node " << i << " is at " << nodes[i] << std::endl;
            if (nodes[i].GetSolutionStepValue(DISTANCE) > 0) {
                volume+=(1.0/nodes.size()); //heyho
            } 
        }
        return volume;
    }

    /*This method is completly useless since it does the same calculation as the previous 
    one but in a different way. Helps to illustrate use of edges
    */
    template<class TGeometryType>
    static double EdgesApproximation(
        const TGeometryType& rVoxel        
    ) {
        double volume = 0;
        GeometryArrayType edges = rVoxel.GenerateEdges();
        PointsArrayType nodes = rVoxel.Points();
        for (int i = 0; i < edges.size(); i++) {
            PointsArrayType ends = edges[i].Points();
            if(ends[0].GetSolutionStepValue(DISTANCE) > 0 && ends[1].GetSolutionStepValue(DISTANCE) > 0) {
                volume+=1.0/edges.size();
            } else if(
                ends[0].GetSolutionStepValue(DISTANCE) > 0 && ends[1].GetSolutionStepValue(DISTANCE) < 0 || 
                ends[0].GetSolutionStepValue(DISTANCE) < 0 && ends[1].GetSolutionStepValue(DISTANCE) > 0 ) {
                volume+=1.0/(edges.size()*2);
            }
        }
        return volume;
    }
    
    /**
     * @brief Aproximates the actual volume inside the voxel 
     * @param rVoxel references to the voxel whose actual volume will be approximated
     * @param rTriangles references to the triangles which intersect the voxel at some edge.
     * @return Approximated volume 
     * @note This approximation finds the portion of each edge that is part of the volume (using
     * intersection point with triangles of the mesh). Even if this class is templated for both 
     * parameters, it will only work with intersecting TRIANGLES, since the utility used to compute
     * the intersection does not allow templating.
     */  
    template<class TGeometryType, class TGeometryArrayType>
    static double EdgesPortionApproximation(
        const TGeometryType& rVoxel,  
        const TGeometryArrayType& rTriangles     
    ) {
        double volume = 0;
        GeometryArrayType edges = rVoxel.GenerateEdges();
        std::vector<double> Distances;

        for (int i = 0; i < edges.size(); i++) {
            Distances.push_back(0);
            PointsArrayType ends = edges[i].Points();
            //std::cout << "Edge " << i << " has nodes " << ends[0] << " " << ends[1] << std::endl;

            for (auto triangle : rTriangles) {
                array_1d<double,3> intersection;
                int result = IntersectionUtilities::ComputeTriangleLineIntersection(triangle,ends[0],ends[1],intersection);
                
                if(result == 1) {
                    double Dist = Distance(ends[0], intersection);
                    Distances.push_back(Dist);
                }  
            } 
            Distances.push_back(Distance(ends[0],ends[1]));
            std::sort(Distances.begin(),Distances.end());       //WOULD A SET BE MORE EFFICIENT?     
            double edgePortion = VolumeInsideVoxelUtility::EdgeFilledPortion(Distances, ends);
            volume += edgePortion/edges.size();  
            Distances.clear();              
        }
        return volume;
    }

     /**
     * @brief Aproximates the actual volume inside the voxel 
     * @param rVoxel references to the voxel whose actual volume will be approximated
     * @param rTriangles references to the triangles which intersect the voxel at some edge.
     * @return Approximated volume 
     * @note This approximation finds the portion of each node that represents volume (using
     * intersection points with triangles of the mesh). Even if this class is templated for both 
     * parameters, it will only work with intersecting TRIANGLES, since the utility used to compute
     * the intersection does not allow templating.
     */  
    static double NodesGeometricalApproximation(
        const GeometryType& rVoxel,  
        const GeometryArrayType& rTriangles     
    ) {
        return 0;    
    }

    /**
     * @brief Aproximates the actual area inside a quadrilateral 
     * @param rFace references to the quadrilateral3D4 whose actual area will be approximated
     * @param rTriangles references to the triangles which intersect the quadrilateral at some edge.
     * @return Approximated area 
     */  
    static double NodesGeometrical2D(
        const GeometryType& rFace,  
        const GeometryArrayType& rTriangles     
    ) {
        double area = 0;
        GeometryArrayType edges = rFace.GenerateEdges();
        PointsArrayType nodes = rFace.Points(); 
        std::vector<std::pair<double,double>> MinDistanceToNode(edges.size(),{0.5,0.5}); 
        //each pair represents an edge and contains as first() the minimum distance between
        //ends[1] and an intersection and as second() the minimum distance between ends[1] and an intersection 

        std::vector<double> Length(edges.size()); 
        for(int i = 0; i < edges.size(); i++) {
            PointsArrayType ends = edges[i].Points();
            double l = Distance(ends[0], ends[1]);
            Length[i] = l;
        }

        for (int i = 0; i < rTriangles.size(); i++) {
            //We will iterate through the edges using a while loop, so that if a triangles intersects 2 edges (unlikely 
            //but possible), only one will be taken into account to create the matrixes.
            int result = 0; 
            array_1d<double,3> intersection;
            int j = 0;
            while(!result && j < edges.size()) { 
                PointsArrayType ends = edges[j].Points();
                result = IntersectionUtilities::ComputeTriangleLineIntersection(rTriangles[i],ends[0],ends[1],intersection);

                if (result) {
                    double dist = Distance(ends[0], intersection);
                    if ( dist < (MinDistanceToNode[j].first*Length[j])) {
                        MinDistanceToNode[j].first = dist/Length[j];
                    } 

                    double dist2 = Distance(ends[1], intersection);
                    if (dist2 < (MinDistanceToNode[j].second*Length[j])) {
                        MinDistanceToNode[j].second = dist2/Length[j];
                    } 
                }
                j++;
            }
        }

        std::vector<std::vector<double>> neighbours{{3,1},{0,2},{1,3},{0,2}};  
        for(int i = 0; i < nodes.size(); i++ ) {
            double factor = GetFactor(nodes, neighbours,i);
            double PartialArea;
            if (nodes[i].GetSolutionStepValue(DISTANCE) > 0) {
                PartialArea = factor*MinDistanceToNode[(i+3)%4].second*MinDistanceToNode[i].first;
            } else  {
                PartialArea = 1.0/nodes.size() - factor*MinDistanceToNode[(i+3)%4].second*MinDistanceToNode[i].first;
            }
            area += PartialArea;
        }
        return area;    
    }

    /**
     * @brief Returns the distance between two 3D points.
     * @param rPoint0 reference to the first point
     * @param rPoint1 reference to the second point
     * @return Distance 
     */  
    static double Distance(const NodeType& Point0, const NodeType& Point1) {
        const double lx = Point0.X() - Point1.X();
        const double ly = Point0.Y() - Point1.Y();
        const double lz = Point0.Z() - Point1.Z();

        const double length = lx * lx + ly * ly + lz * lz;

        return std::sqrt( length );
    }

    /**
     * @brief Returns the distance between two 3D points.
     * @param rPoint0 reference to the first point
     * @param rPoint1 reference an array of 3 coordinates representing the second point
     * @return Distance 
     */  
    static double Distance(const NodeType& Point0, const array_1d<double,3>& Point1) {
        const double lx = Point0.X() - Point1[0];
        const double ly = Point0.Y() - Point1[1];
        const double lz = Point0.Z() - Point1[2];

        const double length = lx * lx + ly * ly + lz * lz;

        return std::sqrt( length );
    }
    
    /**
     * @brief Aproximates the portion of the edge that represents volume
     * @param rDistances references to a sorted vector containing the distances of each intersecting point with the edge
     * @param rEnds references to the nodes at both sides of the edge
     * @return Approximated volume 
     */  
    static const double EdgeFilledPortion(std::vector<double>& Distances, const PointsArrayType& rEnds) {
        double Length = Distances[Distances.size() - 1];
        double portion = 0; 
        bool inside;       

        //casuistics: Both nodes inside, both nodes outside, one node inside and one outside
        if (rEnds[0].GetSolutionStepValue(DISTANCE) > 0 && rEnds[1].GetSolutionStepValue(DISTANCE) > 0) {
            bool inside = true;
            if(Distances.size() % 2 == 1) Distances.pop_back();
            if (Distances.size() == 2) return 1;
        } else if (rEnds[0].GetSolutionStepValue(DISTANCE) > 0 && rEnds[1].GetSolutionStepValue(DISTANCE) < 0) {
            if(Distances.size() % 2 == 0) Distances.pop_back();
            inside = true;
        } else if (rEnds[0].GetSolutionStepValue(DISTANCE) < 0 && rEnds[1].GetSolutionStepValue(DISTANCE) > 0) {
            if(Distances.size() % 2 == 0) Distances.pop_back();
            inside = false;
        } else {    //rEnds[0].GetSolutionStepValue(DISTANCE) < 0 && rEnds[1].GetSolutionStepValue(DISTANCE) < 0
            if (Distances.size() % 2 == 1) Distances.pop_back();
            if (Distances.size() == 2) return 0;
            inside = false;
        }

        for(int i = 1; i < Distances.size(); i++) {
            if (inside) {
                portion += abs(Distances[i]-Distances[i-1])/Length;
                inside = false;
            } else inside = true;             
        }
        return portion;
    }

private:

    ///@name Private static Member Variables
    ///@{

    ///@}
    ///@name Private member Variables
    ///@{

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    static double GetFactor(const PointsArrayType& nodes, const std::vector<std::vector<double>>& neighbours,const int node) {

        if( nodes[node].GetSolutionStepValue(DISTANCE) > 0 && nodes[neighbours[node][0]].GetSolutionStepValue(DISTANCE) < 0 &&
            nodes[neighbours[node][1]].GetSolutionStepValue(DISTANCE) < 0 || nodes[node].GetSolutionStepValue(DISTANCE) < 0 && 
            nodes[neighbours[node][0]].GetSolutionStepValue(DISTANCE) > 0 && nodes[neighbours[node][1]].GetSolutionStepValue(DISTANCE) > 0) {
                return 0.5;
            }
        return 1;
    }

}; /* Class VoxelInsideVolumeUtility */

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

}  /* namespace Kratos.*/

#endif /* KRATOS_VOXEL_INSIDE_VOLUME  defined */