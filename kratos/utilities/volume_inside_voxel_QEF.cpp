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

//Project includes
#include "volume_inside_voxel_qef.h"

namespace Kratos { 

    double VolumeInsideVoxelQEF::SimpleNodesQEFApproximation(
        const GeometryType& rVoxel,  
        const GeometryArrayType& rTriangles     
    ) {
        double volume = 0;
        GeometryArrayType edges = rVoxel.GenerateEdges();
        PointsArrayType Nodes = rVoxel.Points();

        array_1d<double,3> qef = QEF::QEFPoint(rVoxel,rTriangles); 
        //this is unefficient since we will repeat the same calculations to find the intersections afterwards 
        
        DenseMatrix<unsigned int> NodesInFaces(4,6);
        rVoxel.NodesInFaces(NodesInFaces);

        for(int i = 0; i < NodesInFaces.size2(); i++) {
            double Portion = GetPortion(Nodes,NodesInFaces,i);
            double dist = NormalizedDistanceToQEF(Nodes, NodesInFaces, qef, i);
            
            double PartialVolume = Portion*abs(dist)/3.0;   //Volume of a piramid
            volume += PartialVolume;
        }
        //if (volume == 0) return EdgesPortionApproximation(rVoxel,rTriangles);
        
        return volume;
    }

    /***********************************************************************************
     **********************************************************************************/

    double VolumeInsideVoxelQEF::FacesPortionQEFApproximation(
        const GeometryType& rVoxel,  
        const GeometryArrayType& rTriangles     
    ) {
        double volume = 0;
        GeometryArrayType Faces = rVoxel.GenerateFaces();

        array_1d<double,3> qef = QEF::QEFPoint(rVoxel,rTriangles);
        //this is unefficient since we will repeat the same calculations to find the intersections afterwards 

        for(int i = 0; i < Faces.size(); i++) {
            double Portion = EdgesPortionApproximation(Faces[i],rTriangles);
            double dist = NormalizedDistanceToQEF(Faces[i], qef, i);
            
            double PartialVolume = Portion*abs(dist)/3.0;   //Volume of a piramid
            volume += PartialVolume;
        }
        
        return volume;
    }

    /***********************************************************************************
     **********************************************************************************/

    double VolumeInsideVoxelQEF::VoxelVolumeQEFApproximation(
        const GeometryType& rVoxel,  
        const GeometryArrayType& rTriangles     
    ) {
        double Volume = 0;
        GeometryArrayType Faces = rVoxel.GenerateFaces();

        array_1d<double,3> QEF = QEF::QEFPoint(rVoxel,rTriangles); 
        //this is unefficient since we will repeat the same calculations to find the intersections afterwards 

        for(int i = 0; i < Faces.size(); i++) {
            double Portion = VoxelVolume2D(Faces[i],rTriangles);
            double Dist = NormalizedDistanceToQEF(Faces[i], QEF, i);
            
            double PartialVolume = Portion*abs(Dist)/3.0;   //Volume of a piramid
            Volume += PartialVolume;
            //KRATOS_WATCH(PartialVolume);    
        }
        
        if (Volume > 1) return NodesApproximation(rVoxel); //if our approximation fails, use a simpler one with nearly no additional cost
        KRATOS_ERROR_IF(Volume < 0) << "Volume of a mesh element less than 0" << std::endl;
        return Volume;
    }

    /***********************************************************************************
     **********************************************************************************/

    double VolumeInsideVoxelQEF::HexaVolumeQEFApproximation(
        const GeometryType& rVoxel,  
        const GeometryArrayType& rTriangles     
    ) {
        double Volume = 0;
        GeometryArrayType Faces = rVoxel.GenerateFaces();

        array_1d<double,3> QEF = QEF::QEFPoint(rVoxel,rTriangles); 
        //this is unefficient since we will repeat the same calculations to find the intersections afterwards 

        for(int i = 0; i < Faces.size(); i++) {
            double Portion = HexaVolume2D(Faces[i],rTriangles);
            double Dist = NormalizedDistanceToQEF(Faces[i], QEF, i);
            
            double PartialVolume = Portion*abs(Dist)/3.0;   //Volume of a piramid
            Volume += PartialVolume;
            //KRATOS_WATCH(PartialVolume);    
        }
        
        if (Volume > 1) return NodesApproximation(rVoxel); //if our approximation fails, use a simpler one with nearly no additional cost
        KRATOS_ERROR_IF(Volume < 0) << "Volume of a mesh element less than 0" << std::endl;
        return Volume;
    }

    /***********************************************************************************
     **********************************************************************************/

     double VolumeInsideVoxelQEF::GetPortion(PointsArrayType& Nodes,const DenseMatrix<unsigned int>& NodesInFaces, int& face) 
    {
        double Portion = 0;
        for(int i = 0; i < NodesInFaces.size1(); i++) {
            if( Nodes[NodesInFaces(i,face)].GetSolutionStepValue(DISTANCE) > 0) Portion += 1.0/4;
        }
        return Portion;
    }

    /***********************************************************************************
     **********************************************************************************/

    double VolumeInsideVoxelQEF::NormalizedDistanceToQEF(
        PointsArrayType& Nodes,
        const DenseMatrix<unsigned int>& NodesInFaces, 
        const array_1d<double,3>& Point, int& face) 
    {
        array_1d<double, 3> edge1 = Nodes[NodesInFaces(1,face)] - Nodes[NodesInFaces(0,face)];
        array_1d<double, 3> edge2 = Nodes[NodesInFaces(2,face)] - Nodes[NodesInFaces(0,face)];
        array_1d<double, 3> mNormal;
        MathUtils<double>::UnitCrossProduct(mNormal, edge1, edge2);

        const double mConstant =  inner_prod(mNormal, Nodes[NodesInFaces(0,face)]);
        double Side = Distance(Nodes[NodesInFaces(1,face)],Nodes[NodesInFaces(0,face)]);
        double Distance = inner_prod(mNormal,Point) - mConstant;

        return Distance/Side;
    }

    /***********************************************************************************
     **********************************************************************************/
    
    double VolumeInsideVoxelQEF::NormalizedDistanceToQEF(
        GeometryType& rFace, 
        const array_1d<double,3>& Point, int& face) 
    {
        PointsArrayType Nodes = rFace.Points();
        array_1d<double, 3> edge1 = Nodes[1]- Nodes[0];
        array_1d<double, 3> edge2 = Nodes[2] - Nodes[0];
        array_1d<double, 3> mNormal;
        MathUtils<double>::UnitCrossProduct(mNormal, edge1, edge2);

        const double mConstant =  inner_prod(mNormal, Nodes[0]);
        double Side = Distance(Nodes[1],Nodes[0]);
        double Distance = inner_prod(mNormal,Point) - mConstant;

        return Distance/Side;
    }

} //Namespace Kratos