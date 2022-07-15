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
        for (int i = 0; i < 8; i++) {
            if (nodes[i].GetSolutionStepValue(DISTANCE) > 0) {
                volume+=0.125; //heyho
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
        for (int i = 0; i < 12; i++) {
            PointsArrayType ends = edges[i].Points();
            if(ends[0].GetSolutionStepValue(DISTANCE) > 0 && ends[1].GetSolutionStepValue(DISTANCE) > 0) {
                volume+=1.0/12;
            } else if(
                ends[0].GetSolutionStepValue(DISTANCE) > 0 && ends[1].GetSolutionStepValue(DISTANCE) < 0 || 
                ends[0].GetSolutionStepValue(DISTANCE) < 0 && ends[1].GetSolutionStepValue(DISTANCE) > 0 ) {
                volume+=1.0/24;
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
     * intersection point with triangles of the mesh)
     */  
    template<class TGeometryType, class TGeometryArrayType>
    static double EdgesPortionApproximation(
        const TGeometryType& rVoxel,  
        const TGeometryArrayType& rTriangles     
    ) {
        double volume = 0;
        GeometryArrayType edges = rVoxel.GenerateEdges();
        for (int i = 0; i < 12; i++) {
            std::vector<double> Distances;
            Distances.push_back(0);
            PointsArrayType ends = edges[i].Points();

            for (auto triangle : rTriangles) {
                array_1d<double,3> intersection;
                int result = IntersectionUtilities::ComputeTriangleLineIntersection(triangle,ends[0],ends[1],intersection);
                
                if(result == 1) {
                    double Dist = Distance(ends[0], intersection);
                    Distances.push_back(Dist);
                    //std::cout << std::endl << "Intersection with edge " << i << ": " << ends[0] << " " << ends[1] << std::endl;  
                }  
            } 
            Distances.push_back(Distance(ends[0],ends[1]));
            std::sort(Distances.begin(),Distances.end());       //WOULD A MAP BE MORE EFFICIENT?     
            double edgePortion = VolumeInsideVoxelUtility::EdgeFilledPortion(Distances, ends);
            volume += edgePortion/12;                
        }
        return volume;
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

    /**
     * @brief Aproximates the portion of the edge that represents volume
     * @param rDistances references to a vector containing the distances of each intersecting point with the edge
     * @param rEnds references to the nodes at both sides of the edge
     * @return Approximated volume 
     * @note This approximation assigns a fraction of volume (1/8) to each node of the
     * voxel, and counts it as volume if the node is inside the object (NodeDistance > 0)
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
        //std::cout << "edge: " << rEnds[0] << " " << rEnds[1] << "contributed with " << portion << std::endl;  
        return portion;
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

}; /* Class VoxelInsideVolumeUtility */

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

}  /* namespace Kratos.*/

#endif /* KRATOS_VOXEL_INSIDE_VOLUME  defined */