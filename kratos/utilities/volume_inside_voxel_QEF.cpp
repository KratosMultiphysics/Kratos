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

//Project includes
#include "volume_inside_voxel_qef.h"

namespace Kratos { 

    double VolumeInsideVoxelQEF::SimpleNodesQEFApproximation(
        const GeometryType& rVoxel,  
        const GeometryArrayType& rTriangles     
    ) {
        double volume = 0;
        PointsArrayType nodes = rVoxel.Points();

        array_1d<double,3> qef = QuadraticErrorFunction::QuadraticErrorFunctionPoint(rVoxel,rTriangles); 
        //this is unefficient since we will repeat the same calculations to find the intersections afterwards 
        
        DenseMatrix<unsigned int> nodes_in_faces(4,6);
        rVoxel.NodesInFaces(nodes_in_faces);

        for(int i = 0; i < nodes_in_faces.size2(); i++) {
            double portion = GetPortion(nodes,nodes_in_faces,i);
            double dist = NormalizedDistanceToQEF(nodes, nodes_in_faces, qef, i);
            
            double partial_volume = portion*abs(dist)/3.0;   //Volume of a piramid
            volume += partial_volume;
        }
        
        return volume;
    }

    /***********************************************************************************
     **********************************************************************************/

    double VolumeInsideVoxelQEF::VolumeQEFApproximation(
        const GeometryType& rVoxel,  
        const GeometryArrayType& rTriangles     
    ) {
        double volume = 0;
        GeometryArrayType faces = rVoxel.GenerateFaces();

        array_1d<double,3> qef = QuadraticErrorFunction::QuadraticErrorFunctionPoint(rVoxel,rTriangles); 
        //this is unefficient since we will repeat the same calculations to find the intersections afterwards 

        for(int i = 0; i < faces.size(); i++) {
            double portion = FaceArea(faces[i],rTriangles);
            double dist = NormalizedDistanceToQEF(faces[i], qef, i);
            
            double partial_volume = portion*abs(dist)/3.0;   //Volume of a piramid
            volume += partial_volume;
        }
        
        if (volume > 1) return NodesApproximation(rVoxel); //if our approximation fails, use a simpler one with nearly no additional cost
        KRATOS_ERROR_IF(volume < 0) << "Volume of a mesh element less than 0" << std::endl;
        return volume;
    }

    /***********************************************************************************
     **********************************************************************************/

     double VolumeInsideVoxelQEF::GetPortion(
        PointsArrayType& rNodes,
        const DenseMatrix<unsigned int>& rNodesInFaces,
        int Face) 
    {
        double portion = 0;
        for(int i = 0; i < rNodesInFaces.size1(); i++) {
            if( rNodes[rNodesInFaces(i,Face)].GetSolutionStepValue(DISTANCE) > 0) portion += 1.0/4;
        }
        return portion;
    }

    /***********************************************************************************
     **********************************************************************************/

    double VolumeInsideVoxelQEF::NormalizedDistanceToQEF(
        PointsArrayType& rNodes,
        const DenseMatrix<unsigned int>& rNodesInFaces, 
        const array_1d<double,3>& rPoint, 
        int Face) 
    {
        array_1d<double, 3> edge1 = rNodes[rNodesInFaces(1,Face)] - rNodes[rNodesInFaces(0,Face)];
        array_1d<double, 3> edge2 = rNodes[rNodesInFaces(2,Face)] - rNodes[rNodesInFaces(0,Face)];
        array_1d<double, 3> v_normal;
        MathUtils<double>::UnitCrossProduct(v_normal, edge1, edge2);

        const double mConstant =  inner_prod(v_normal, rNodes[rNodesInFaces(0,Face)]);
        double side = norm_2(rNodes[rNodesInFaces(1,Face)].Coordinates() - rNodes[rNodesInFaces(0,Face)].Coordinates());
        double distance = inner_prod(v_normal,rPoint) - mConstant;

        return distance/side;
    }

    /***********************************************************************************
     **********************************************************************************/
    
    double VolumeInsideVoxelQEF::NormalizedDistanceToQEF(
        GeometryType& rFace, 
        const array_1d<double,3>& rPoint,
        int face) 
    {
        PointsArrayType nodes = rFace.Points();
        array_1d<double, 3> edge1 = nodes[1]- nodes[0];
        array_1d<double, 3> edge2 = nodes[2] - nodes[0];
        
        array_1d<double, 3> v_normal;
        MathUtils<double>::UnitCrossProduct(v_normal, edge1, edge2);
        
        const double v_constant =  inner_prod(v_normal, nodes[0]);
        double side = norm_2(nodes[1].Coordinates() - nodes[0].Coordinates());
        double distance = inner_prod(v_normal,rPoint) - v_constant;

        return distance/side;
    }

} //Namespace Kratos