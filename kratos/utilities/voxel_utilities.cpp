//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Ariadna Cortes
//

// System includes

// External includes

// Project includes
#include "voxel_utilities.h"

namespace Kratos 
{
constexpr std::size_t VoxelUtilities::mNeighbours[4][2];

/**********************************************************************************/
/**********************************************************************************/

double VoxelUtilities::NodesApproximation(const GeometryType& rVoxel) {
    double volume = 0;
    const auto& nodes = rVoxel.Points(); //PointsArrayType
    const std::size_t n = nodes.size();
    for (std::size_t i = 0; i < n; i++) {
        if (nodes[i].GetSolutionStepValue(DISTANCE) > 0) {
            volume+=(1.0/nodes.size()); 
        } 
    }
    return volume;
}

/**********************************************************************************/
/**********************************************************************************/

double VoxelUtilities::EdgesPortionApproximation(
    const GeometryType& rVoxel,  
    const GeometryArrayType& rTriangles     
) {
    double volume = 0;
    GeometryArrayType edges = rVoxel.GenerateEdges();
    std::vector<double> distances;
    const std::size_t n = edges.size();
    for (std::size_t i = 0; i < n; i++) {
        distances.push_back(0);
        const auto& ends = edges[i].Points();

        for (auto& r_triangle : rTriangles) {
            array_1d<double,3> intersection;
            int result = IntersectionUtilities::ComputeTriangleLineIntersection(r_triangle,ends[0],ends[1],intersection);
            
            if(result == 1) {
                const double dist = norm_2(ends[0].Coordinates() - intersection);
                distances.push_back(dist);
            }  
        } 
        distances.push_back(norm_2(ends[0].Coordinates() - ends[1].Coordinates()));
        std::sort(distances.begin(),distances.end());        
        const double edge_portion = VoxelUtilities::EdgeFilledPortion(distances, ends);
        volume += edge_portion/edges.size();  
        distances.clear();              
    }
    return volume;
}

/**********************************************************************************/
/**********************************************************************************/

double VoxelUtilities::FaceArea(
    const GeometryType& rFace,  
    const GeometryArrayType& rTriangles     
) {
    double area = 0;
    GeometryArrayType edges = rFace.GenerateEdges();
    auto nodes = rFace.Points(); 
    const double face_area = PointsArea(nodes);
    std::vector<std::pair<double,double>> min_distance_to_node(edges.size(),{1,1}); 
    
    int nodes_inside = 0;
    const std::size_t nodes_size = nodes.size();
    for (std::size_t i = 0; i < nodes_size; i++) {
        if (nodes[i].GetSolutionStepValue(DISTANCE) > 0) nodes_inside++;
    }

    if(nodes_inside == 3) {
        for (std::size_t i = 0; i < nodes_size; i++)  {
            nodes[i].GetSolutionStepValue(DISTANCE) = -nodes[i].GetSolutionStepValue(DISTANCE);
        }
        const double area = FaceArea(rFace,rTriangles);
        for (std::size_t i = 0; i < nodes_size; i++) {
            nodes[i].GetSolutionStepValue(DISTANCE) = -nodes[i].GetSolutionStepValue(DISTANCE);
        }
        return 1-area;
    }
    
    const std::size_t edges_size = edges.size();
    std::vector<double> length(edges_size); 
    for(std::size_t i = 0; i < edges_size; i++) {
        const auto& ends = edges[i].Points();
        const double l = norm_2(ends[0].Coordinates() - ends[1].Coordinates());
        length[i] = l;
    }

    for (std::size_t i = 0; i < rTriangles.size(); i++) {
        int result = 0; 
        array_1d<double,3> intersection;
        std::size_t j = 0;
        while(!result && j < edges_size) { 
            const auto& ends = edges[j].Points();
            result = IntersectionUtilities::ComputeTriangleLineIntersection(rTriangles[i],ends[0],ends[1],intersection);

            if (result) {
                const double dist =  norm_2(ends[0].Coordinates() - intersection);
                if ( dist < (min_distance_to_node[j].first*length[j])) {
                    min_distance_to_node[j].first = dist/length[j];
                } 

                const double dist2 =  norm_2(ends[1].Coordinates() - intersection);
                if (dist2 < (min_distance_to_node[j].second*length[j])) {
                    min_distance_to_node[j].second = dist2/length[j];
                } 
            }
            j++;
        }
    }

    for(std::size_t i = 0; i < nodes_size; i++ ) {
        array_1d<double,3> v_left{nodes[i].X() -nodes[VoxelUtilities::mNeighbours[i][0]].X(), nodes[i].Y() -nodes[VoxelUtilities::mNeighbours[i][0]].Y(), nodes[i].Z() -nodes[VoxelUtilities::mNeighbours[i][0]].Z()};
        array_1d<double,3> v_right{nodes[i].X() -nodes[VoxelUtilities::mNeighbours[i][1]].X(), nodes[i].Y() -nodes[VoxelUtilities::mNeighbours[i][1]].Y(), nodes[i].Z() -nodes[VoxelUtilities::mNeighbours[i][1]].Z()};

        const double Case = GetCase(nodes,i);
        double partial_area;
        double factor = 1;
        double left = 0;
        double right = 0;
        if (nodes[i].GetSolutionStepValue(DISTANCE) > 0) {
            if (Case == 0) {
                factor = 0.5;
                left = min_distance_to_node[(i+3)%4].second;
                right = min_distance_to_node[i].first;
            } else if (Case == 1)  {
                left = min_distance_to_node[(i+3)%4].second;
                right = std::min(0.5,min_distance_to_node[i].first);
            } else if (Case == 2) {
                left = std::min(0.5,min_distance_to_node[(i+3)%4].second);
                right = min_distance_to_node[i].first;
            } else if (Case == 3) {
                left = std::min(0.5,min_distance_to_node[(i+3)%4].second);
                right = std::min(0.5,min_distance_to_node[i].first);
            }
            NodePtrType int_left(new Node<3>(1, nodes[i].X() + left*v_left[0], nodes[i].Y() + left*v_left[1], nodes[i].Z() + left*v_left[2]));
            NodePtrType int_right(new Node<3>(2, nodes[i].X() + right*v_right[0], nodes[i].Y() + right*v_right[1], nodes[i].Z() + right*v_right[2]));
            NodePtrType c(new Node<3>(2, nodes[i].X() + left*v_left[0] + right*v_right[0], 
                                                nodes[i].Y() + left*v_left[1] + right*v_right[1], 
                                                nodes[i].Z() + left*v_left[2] + right*v_right[2]));
            PointsArrayType points(4);
            points(0) = int_left;
            points(1) = &nodes[i];
            points(2) = int_right;
            points(3) = c;
            partial_area = factor*PointsArea(points)/face_area;
        } else {
            NodePtrType max_left(new Node<3>(1, nodes[i].X() + 0.5*v_left[0], nodes[i].Y() + 0.5*v_left[1], nodes[i].Z() + 0.5*v_left[2]));
            NodePtrType max_right(new Node<3>(2, nodes[i].X() + 0.5*v_right[0], nodes[i].Y() + 0.5*v_right[1], nodes[i].Z() + 0.5*v_right[2]));
            NodePtrType max_c(new Node<3>(2, nodes[i].X() + 0.5*v_left[0] + 0.5*v_right[0], 
                                                nodes[i].Y() + 0.5*v_left[1] + 0.5*v_right[1], 
                                                nodes[i].Z() + 0.5*v_left[2] + 0.5*v_right[2]));
            PointsArrayType points(4);
            points(3) = max_left;
            points(0) = &nodes[i];
            points(1) = max_right;
            points(2) = max_c;
            const double max_volume = PointsArea(points)/face_area;

            if (Case == 1) {
                left = 0.5;
                right = std::min(0.5,min_distance_to_node[i].first);
            } else if (Case == 2) {
                left = std::min(0.5,min_distance_to_node[(i+3)%4].second);
                right = 0.5;
            } else if (Case == 3) {
                left = std::min(0.5,min_distance_to_node[(i+3)%4].second);
                right = std::min(0.5,min_distance_to_node[i].first);
            }
            NodePtrType int_left(new Node<3>(1, nodes[i].X() + left*v_left[0], nodes[i].Y() + left*v_left[1], nodes[i].Z() + left*v_left[2]));
            NodePtrType int_right(new Node<3>(2, nodes[i].X() + right*v_right[0], nodes[i].Y() + right*v_right[1], nodes[i].Z() + right*v_right[2]));
            NodePtrType c(new Node<3>(2, nodes[i].X() + left*v_left[0] + right*v_right[0], 
                                                nodes[i].Y() + left*v_left[1] + right*v_right[1], 
                                                nodes[i].Z() + left*v_left[2] + right*v_right[2]));
            points(3) = int_left;
            points(0) = &nodes[i];
            points(1) = int_right;
            points(2) = c;

            if (Case != 0) {
                partial_area =  max_volume -factor*PointsArea(points)/face_area;
            } else {
                partial_area = 0;
            }
        }
        area += partial_area;
    }
    return area;    
}

/**********************************************************************************/
/**********************************************************************************/

double VoxelUtilities::PointsArea(const PointsArrayType& rPoints) {
    const NodeType& p0 = rPoints[0];
    const NodeType& p1 = rPoints[1];
    const NodeType& p2 = rPoints[2];
    const NodeType& p3 = rPoints[3];

    const array_1d<double,3> u1 = (p1.Coordinates() - p0.Coordinates());
    const array_1d<double,3> v1 = (p3.Coordinates() - p0.Coordinates());

    const array_1d<double,3> u2 = (p3.Coordinates() - p2.Coordinates());
    const array_1d<double,3> v2 = (p1.Coordinates() - p2.Coordinates());

    array_1d<double,3> area1;
    array_1d<double,3> area2;

    MathUtils<double>::CrossProduct(area1,u1,v1);
    MathUtils<double>::CrossProduct(area2,u2,v2);

    return 0.5*(norm_2(area1) + norm_2(area2));
}

/**********************************************************************************/
/**********************************************************************************/

double VoxelUtilities::EdgeFilledPortion(std::vector<double>& Distances, const PointsArrayType& rEnds) {
    const double length = Distances[Distances.size() - 1];
    double portion = 0; 
    bool inside = true;       

    const bool left_inside = rEnds[0].GetSolutionStepValue(DISTANCE) > 0;
    const bool right_inside = rEnds[1].GetSolutionStepValue(DISTANCE) > 0; 
    
    //casuistics: Both nodes inside, both nodes outside, one node inside and one outside
    if (left_inside && right_inside ) {
        if(Distances.size() % 2 == 1) Distances.pop_back();
        if (Distances.size() == 2) return 1;
    } else if (left_inside && !right_inside) {
        if(Distances.size() % 2 == 0) Distances.pop_back();
    } else if (!left_inside && right_inside) {
        if(Distances.size() % 2 == 0) Distances.pop_back();
        inside = false;
    } else {    //rEnds[0].GetSolutionStepValue(DISTANCE) < 0 && rEnds[1].GetSolutionStepValue(DISTANCE) < 0
        if (Distances.size() % 2 == 1) Distances.pop_back();
        if (Distances.size() == 2) return 0;
        inside = false;
    }

    for(std::size_t i = 1; i < Distances.size(); i++) {
        if (inside) {
            portion += std::abs(Distances[i]-Distances[i-1])/length;
            inside = false;
        } else inside = true;             
    }
    return portion;
}

/**********************************************************************************/
/**********************************************************************************/

double VoxelUtilities::GetFactor(
    const PointsArrayType& rNodes, 
    const int NodeIndex) 
{
    const bool me_inside = rNodes[NodeIndex].GetSolutionStepValue(DISTANCE) > 0;
    const bool left_inside = rNodes[VoxelUtilities::mNeighbours[NodeIndex][0]].GetSolutionStepValue(DISTANCE) > 0;
    const bool right_inside =  rNodes[VoxelUtilities::mNeighbours[NodeIndex][1]].GetSolutionStepValue(DISTANCE) < 0;

    if ((me_inside && !left_inside && !right_inside) || (!me_inside && left_inside && right_inside)) {
            return 0.5;
    }
    return 1.0;
}

/**********************************************************************************/
/**********************************************************************************/

int VoxelUtilities::GetCase(
    const PointsArrayType& rNodes,
    const int NodeIndex) 
{
    const int me_inside = (int) rNodes[NodeIndex].GetSolutionStepValue(DISTANCE) > 0;
    const int left_inside = (int) (rNodes[VoxelUtilities::mNeighbours[NodeIndex][0]].GetSolutionStepValue(DISTANCE) > 0);
    const int right_inside = (int) (rNodes[VoxelUtilities::mNeighbours[NodeIndex][1]].GetSolutionStepValue(DISTANCE) > 0);

    if (me_inside) {
        return left_inside*2 + right_inside;
    }
    else {
        return (!left_inside)*2 + !right_inside;
    }
}

}  /* namespace Kratos.*/