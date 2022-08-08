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
    one but in a different way. Helps to illustrate use of Edges
    */
    template<class TGeometryType>
    static double EdgesApproximation(
        const TGeometryType& rVoxel        
    ) {
        double volume = 0;
        GeometryArrayType Edges = rVoxel.GenerateEdges();
        PointsArrayType nodes = rVoxel.Points();
        for (int i = 0; i < Edges.size(); i++) {
            PointsArrayType ends = Edges[i].Points();
            if(ends[0].GetSolutionStepValue(DISTANCE) > 0 && ends[1].GetSolutionStepValue(DISTANCE) > 0) {
                volume+=1.0/Edges.size();
            } else if(
                ends[0].GetSolutionStepValue(DISTANCE) > 0 && ends[1].GetSolutionStepValue(DISTANCE) < 0 || 
                ends[0].GetSolutionStepValue(DISTANCE) < 0 && ends[1].GetSolutionStepValue(DISTANCE) > 0 ) {
                volume+=1.0/(Edges.size()*2);
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
        GeometryArrayType Edges = rVoxel.GenerateEdges();
        std::vector<double> Distances;

        for (int i = 0; i < Edges.size(); i++) {
            Distances.push_back(0);
            PointsArrayType ends = Edges[i].Points();
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
            volume += edgePortion/Edges.size();  
            Distances.clear();              
        }
        return volume;
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
        double Area = 0;
        GeometryArrayType Edges = rFace.GenerateEdges();
        PointsArrayType Nodes = rFace.Points(); 
        std::vector<std::pair<double,double>> MinDistanceToNode(Edges.size(),{0.5,0.5}); 
        //each pair represents an edge and contains as first() the minimum distance between
        //ends[1] and an intersection and as second() the minimum distance between ends[1] and an intersection 

        std::vector<double> Length(Edges.size()); 
        for(int i = 0; i < Edges.size(); i++) {
            PointsArrayType ends = Edges[i].Points();
            double l = Distance(ends[0], ends[1]);
            Length[i] = l;
        }

        for (int i = 0; i < rTriangles.size(); i++) {
            //We will iterate through the edges using a while loop, so that if a triangles intersects 2 edges (unlikely 
            //but possible), only one will be taken into account to create the matrixes.
            int Result = 0; 
            array_1d<double,3> Intersection;
            int j = 0;
            while(!Result && j < Edges.size()) { 
                PointsArrayType ends = Edges[j].Points();
                Result = IntersectionUtilities::ComputeTriangleLineIntersection(rTriangles[i],ends[0],ends[1],Intersection);

                if (Result) {
                    double Dist = Distance(ends[0], Intersection);
                    if ( Dist < (MinDistanceToNode[j].first*Length[j])) {
                        MinDistanceToNode[j].first = Dist/Length[j];
                    } 

                    double Dist2 = Distance(ends[1], Intersection);
                    if (Dist2 < (MinDistanceToNode[j].second*Length[j])) {
                        MinDistanceToNode[j].second = Dist2/Length[j];
                    } 
                }
                j++;
            }
        }

        std::vector<std::vector<double>> Neighbours{{3,1},{0,2},{1,3},{0,2}};  
        for(int i = 0; i < Nodes.size(); i++ ) {
            double Factor = GetFactor(Nodes, Neighbours,i);
            double PartialArea;
            if (Nodes[i].GetSolutionStepValue(DISTANCE) > 0) {
                PartialArea = Factor*MinDistanceToNode[(i+3)%4].second*MinDistanceToNode[i].first;
            } else  {
                PartialArea = 1.0/Nodes.size() - Factor*MinDistanceToNode[(i+3)%4].second*MinDistanceToNode[i].first;
            }
            Area += PartialArea;
        }
        return Area;    
    }

    /**
     * @brief Aproximates the actual area inside a quadrilateral 
     * @param rFace references to the quadrilateral3D4 whose actual area will be approximated
     * @param rTriangles references to the triangles which intersect the quadrilateral at some edge.
     * @return Approximated area 
     */  
    static double NodesGeometricalCases2D(
        const GeometryType& rFace,  
        const GeometryArrayType& rTriangles     
    ) {
        double Area = 0;
        GeometryArrayType Edges = rFace.GenerateEdges();
        PointsArrayType Nodes = rFace.Points(); 
        double FaceArea = TetraVolume(Nodes);
        std::vector<std::pair<double,double>> MinDistanceToNode(Edges.size(),{1,1}); 
        
        int NodesInside = 0;
        for (int i = 0; i < Nodes.size(); i++)  if (Nodes[i].GetSolutionStepValue(DISTANCE) > 0) NodesInside++;

        if(NodesInside == 3) {
            for (int i = 0; i < Nodes.size(); i++)  Nodes[i].GetSolutionStepValue(DISTANCE) = (-1)*Nodes[i].GetSolutionStepValue(DISTANCE);
            double Area = NodesGeometricalCases2D(rFace,rTriangles);
            for (int i = 0; i < Nodes.size(); i++)  Nodes[i].GetSolutionStepValue(DISTANCE) = (-1)*Nodes[i].GetSolutionStepValue(DISTANCE);
            return 1-Area;
        }

        std::vector<double> Length(Edges.size()); 
        for(int i = 0; i < Edges.size(); i++) {
            PointsArrayType ends = Edges[i].Points();
            double l = Distance(ends[0], ends[1]);
            Length[i] = l;
        }

        for (int i = 0; i < rTriangles.size(); i++) {
            int Result = 0; 
            array_1d<double,3> Intersection;
            int j = 0;
            while(!Result && j < Edges.size()) { 
                PointsArrayType ends = Edges[j].Points();
                Result = IntersectionUtilities::ComputeTriangleLineIntersection(rTriangles[i],ends[0],ends[1],Intersection);

                if (Result) {
                    double Dist = Distance(ends[0], Intersection);
                    if ( Dist < (MinDistanceToNode[j].first*Length[j])) {
                        MinDistanceToNode[j].first = Dist/Length[j];
                    } 

                    double Dist2 = Distance(ends[1], Intersection);
                    if (Dist2 < (MinDistanceToNode[j].second*Length[j])) {
                        MinDistanceToNode[j].second = Dist2/Length[j];
                    } 
                }
                j++;
            }
        }

        std::vector<std::vector<double>> Neighbours{{3,1},{0,2},{1,3},{2,0}};  
        for(int i = 0; i < Nodes.size(); i++ ) {
            array_1d<double,3> v_left{Nodes[i].X() -Nodes[Neighbours[i][0]].X(), Nodes[i].Y() -Nodes[Neighbours[i][0]].Y(), Nodes[i].Z() -Nodes[Neighbours[i][0]].Z()};
            array_1d<double,3> v_right{Nodes[i].X() -Nodes[Neighbours[i][1]].X(), Nodes[i].Y() -Nodes[Neighbours[i][1]].Y(), Nodes[i].Z() -Nodes[Neighbours[i][1]].Z()};

            double Case = GetCase(Nodes, Neighbours,i);
            double PartialArea;
            double factor = 1;
            double left = 0;
            double right = 0;
            if (Nodes[i].GetSolutionStepValue(DISTANCE) > 0) {
                if (Case == 0) {
                    factor = 0.5;
                    left = MinDistanceToNode[(i+3)%4].second;
                    right = MinDistanceToNode[i].first;
                }
                else if (Case == 1)  {
                    left = MinDistanceToNode[(i+3)%4].second;
                    right = std::min(0.5,MinDistanceToNode[i].first);
                }
                else if (Case == 2) {
                    left = std::min(0.5,MinDistanceToNode[(i+3)%4].second);
                    right = MinDistanceToNode[i].first;
                }
                else if (Case == 3) {
                    left = std::min(0.5,MinDistanceToNode[(i+3)%4].second);
                    right = std::min(0.5,MinDistanceToNode[i].first);
                }
                NodePtrType Int_left(new Node<3>(1, Nodes[i].X() + left*v_left[0], Nodes[i].Y() + left*v_left[1], Nodes[i].Z() + left*v_left[2]));
                NodePtrType Int_right(new Node<3>(2, Nodes[i].X() + right*v_right[0], Nodes[i].Y() + right*v_right[1], Nodes[i].Z() + right*v_right[2]));
                NodePtrType c(new Node<3>(2, Nodes[i].X() + left*v_left[0] + right*v_right[0], 
                                                    Nodes[i].Y() + left*v_left[1] + right*v_right[1], 
                                                    Nodes[i].Z() + left*v_left[2] + right*v_right[2]));
                PointsArrayType points;
                points.push_back(Int_left);
                points.push_back(&Nodes[i]);
                points.push_back(Int_right);
                points.push_back(c);
                PartialArea = factor*TetraVolume(points)/FaceArea;
            } else  {
                if (Case == 1) {
                    left = 0.5;
                    right = std::min(0.5,MinDistanceToNode[i].first);
                }
                if (Case == 2) {
                    left = std::min(0.5,MinDistanceToNode[(i+3)%4].second);
                    right = 0.5;
                }
                if (Case == 3) {
                    left = std::min(0.5,MinDistanceToNode[(i+3)%4].second);
                    right = std::min(0.5,MinDistanceToNode[i].first);
                }
                NodePtrType Int_left(new Node<3>(1, Nodes[i].X() + left*v_left[0], Nodes[i].Y() + left*v_left[1], Nodes[i].Z() + left*v_left[2]));
                NodePtrType Int_right(new Node<3>(2, Nodes[i].X() + right*v_right[0], Nodes[i].Y() + right*v_right[1], Nodes[i].Z() + right*v_right[2]));
                NodePtrType c(new Node<3>(2, Nodes[i].X() + left*v_left[0] + right*v_right[0], 
                                                    Nodes[i].Y() + left*v_left[1] + right*v_right[1], 
                                                    Nodes[i].Z() + left*v_left[2] + right*v_right[2]));
                PointsArrayType points;
                points.push_back(Int_left);
                points.push_back(&Nodes[i]);
                points.push_back(Int_right);
                points.push_back(c);
                if (Case != 0) PartialArea =  1.0/Nodes.size() -factor*TetraVolume(points)/FaceArea;
                else PartialArea = 0;
            }
            Area += PartialArea;
        }
        return Area;    
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
     * @brief Returns the volume enclosed by a set of points
     * @param rPoints the array of points
     * @return Volume/area inside this points
     */  
    static double TetraVolume(const PointsArrayType& rPoints) {
        GeometryPtrType pGeom =Kratos::make_shared<Quadrilateral3D4<NodeType>>(rPoints);
        return pGeom->Volume();
    }
    static double TriangleVolume(const PointsArrayType& rPoints) {
        GeometryPtrType pGeom =Kratos::make_shared<Triangle3D3<NodeType>>(rPoints);
        return pGeom->Volume();
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
    
    static int GetCase(const PointsArrayType& nodes, const std::vector<std::vector<double>>& neighbours,const int node) {
        if (nodes[node].GetSolutionStepValue(DISTANCE) > 0 ) {
            if (nodes[neighbours[node][0]].GetSolutionStepValue(DISTANCE) < 0 &&
            nodes[neighbours[node][1]].GetSolutionStepValue(DISTANCE) < 0) {
                return 0;
            }
            else if (nodes[neighbours[node][0]].GetSolutionStepValue(DISTANCE) < 0 &&
            nodes[neighbours[node][1]].GetSolutionStepValue(DISTANCE) > 0) {
                return 1;
            }
            else if (nodes[neighbours[node][0]].GetSolutionStepValue(DISTANCE) > 0 &&
            nodes[neighbours[node][1]].GetSolutionStepValue(DISTANCE) < 0) {
                return 2;
            } 
            else  {
                return 3;
            }
        }
        else {
            if (nodes[neighbours[node][0]].GetSolutionStepValue(DISTANCE) > 0 && 
            nodes[neighbours[node][1]].GetSolutionStepValue(DISTANCE) > 0) {
                return 0;
            }
            else if (nodes[neighbours[node][0]].GetSolutionStepValue(DISTANCE) > 0 && 
            nodes[neighbours[node][1]].GetSolutionStepValue(DISTANCE) < 0) {
                return 1;
            }
            else if (nodes[neighbours[node][0]].GetSolutionStepValue(DISTANCE) < 0 && 
            nodes[neighbours[node][1]].GetSolutionStepValue(DISTANCE) > 0) {
                return 2;
            }
            else return 3;
        }
    }

}; /* Class VoxelInsideVolumeUtility */

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

}  /* namespace Kratos.*/

#endif /* KRATOS_VOXEL_INSIDE_VOLUME  defined */