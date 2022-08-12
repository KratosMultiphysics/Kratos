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

// Project includes
#include "qef_utility.h"

namespace Kratos {

    array_1d<double,3> QEF::QEFPoint (
        const GeometryType& rVoxel,  
        const GeometryArrayType& rTriangles     
    ) {
        array_1d<double,3> xPoint;
        array_1d<double,3> center = CalculateCenter(rVoxel);
        VectorType vCenter = center; 
        MatrixType mCenter(3,1);
        column(mCenter,0) = center;
        GeometryArrayType edges = rVoxel.GenerateEdges();

        //Initialize the corresponding matrixes
        MatrixType AtA(3,3,0);  //3x3 matrix initialized to 0
        MatrixType AtB(3,1,0);  //3x1 matrix
        MatrixType BtB(1,1,0);  //1x1 matrix

        for (int i = 0; i < rTriangles.size(); i++) {
            array_1d<double,3> normal = CalculateNormal(rTriangles[i]);
            VectorType vNormal = normal; 
            MatrixType mNormal(3,1);
            column(mNormal,0) = normal;

            //We will iterate through the edges using a while loop, so that if a triangles intersects 2 edges (unlikely 
            //but possible), only one will be taken into account to create the matrixes.
            int result = 0; 
            array_1d<double,3> intersection;
            VectorType vIntersection;
            int j = 0;
            while(!result && j < edges.size()) { 
                PointsArrayType ends = edges[j++].Points();
                result = IntersectionUtilities::ComputeTriangleLineIntersection(rTriangles[i],ends[0],ends[1],intersection);
                vIntersection = intersection;
            }
            if (result) {
                //Fill the matrixes with the corresponding information from the intersection and normal
                MatrixType mNormalTrans = trans(mNormal);
                MatrixType help = prod(mNormal,mNormalTrans);
                AtA = AtA + help;
                double aux = MathUtils<double>::Dot(normal,intersection);
                MatrixType mAux(1,1,aux);
                AtB = AtB + prod(mNormal,mAux);
                BtB = BtB + prod(mAux,mAux); 
            }
        }

        //Find the eigenvalues and eigenvectors to AtA
        MatrixType mEigenvectors;
        MatrixType mEigenvalues;
        const bool converged = MathUtils<double>::GaussSeidelEigenSystem(AtA, mEigenvectors, mEigenvalues);

        MatrixType D(3,3,0);
        for (int i : {0,1,2}) mEigenvalues(i,i) < 1e-12 ? D(i,i) = 0 : D(i,i) = check(1.0/mEigenvalues(i,i), 1e-12);

        MatrixType AtAInverse;  
        MathUtils<double>::BDBtProductOperation(AtAInverse, D, mEigenvectors);
       
        MatrixType AtAc = prod(AtA,mCenter);
        MatrixType solution = prod(AtAInverse, AtB - AtAc) + mCenter;        
        xPoint = column(solution,0);

        return xPoint;
    }


    array_1d<double,3> QEF::CalculateCenter(const GeometryType& rVoxel) {
        PointsArrayType nodes = rVoxel.Points();
        double x, y, z;
        x = y = z = 0;
        for(int i = 0; i < nodes.size(); i++) {
            x += nodes[i].X();
            y += nodes[i].Y();
            z += nodes[i].Z();
        }
        array_1d<double,3> center({x/nodes.size(),y/nodes.size(),z/nodes.size()});
        return center;
    }

    array_1d<double,3> QEF::CalculateNormal(const GeometryType& triangle) {
        PointsArrayType nodes = triangle.Points();
        array_1d<double,3> u = nodes[1] - nodes[0];
        array_1d<double,3> v = nodes[2] - nodes[0];
        array_1d<double,3> normal;
        MathUtils<double>::CrossProduct<array_1d<double,3>,array_1d<double,3>,array_1d<double,3>>(normal,u,v);
        return normal;
    }
}