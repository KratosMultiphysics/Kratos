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
#include "qef_utility.h"

namespace Kratos {

    array_1d<double,3> QuadraticErrorFunction::QuadraticErrorFunctionPoint (
        const GeometryType& rVoxel,  
        const GeometryArrayType& rTriangles     
    ) {
        array_1d<double,3> x_point;
        array_1d<double,3> center = rVoxel.Center().Coordinates();
        MatrixType mat_center(3,1);
        column(mat_center,0) = center;
        GeometryArrayType edges = rVoxel.GenerateEdges();

        //Initialize the corresponding matrixes
        BoundedMatrix<double,3,3>ata = ZeroMatrix(3,3);  
        BoundedMatrix<double,3,1>atb = ZeroMatrix(3,1);

        for (std::size_t i = 0; i < rTriangles.size(); i++) {
            array_1d<double,3> normal = CalculateNormal(rTriangles[i]);
            BoundedMatrix<double,3,1> mat_normal;
            column(mat_normal,0) = normal;

            //We will iterate through the edges using a while loop, so that if a triangles intersects 2 edges (unlikely 
            //but possible), only one will be taken into account to create the matrixes.
            int result = 0; 
            array_1d<double,3> intersection;
            std::size_t j = 0;
            while(!result && j < edges.size()) { 
                PointsArrayType ends = edges[j++].Points();
                result = IntersectionUtilities::ComputeTriangleLineIntersection(rTriangles[i],ends[0],ends[1],intersection);
            }
            if (result) {
                //Fill the matrixes with the corresponding information from the intersection and normal
                MatrixType mat_normal_trans = trans(mat_normal);
                ata += prod(mat_normal,mat_normal_trans);
                double aux = MathUtils<double>::Dot(normal,intersection);
                MatrixType mat_aux(1,1,aux);
                atb += prod(mat_normal,mat_aux);
            }
        }

        //Find the eigenvalues and eigenvectors to ata
        MatrixType mat_eigenvectors;
        MatrixType mat_eigenvalues;
        const bool converged = MathUtils<double>::GaussSeidelEigenSystem(ata, mat_eigenvectors, mat_eigenvalues);

        KRATOS_WARNING_IF("QuadraticErrorFunctionPoint", !converged) << "Method for matrix eigenvalues didn't converge" << std::endl;

        MatrixType d(3,3,0);
        for (int i : {0,1,2}) {
            if (mat_eigenvalues(i,i) < 1e-12) {
                d(i,i) = 0; 
            } else {
                d(i,i) = 1.0/mat_eigenvalues(i,i) > 1.0e-12 ? 1.0/mat_eigenvalues(i,i) : 0.0;
            }
        }

        MatrixType ata_inverse;  
        MathUtils<double>::BDBtProductOperation(ata_inverse, d, mat_eigenvectors);
       
        MatrixType ata_c = prod(ata,mat_center);
        MatrixType solution = prod(ata_inverse, atb - ata_c) + mat_center;        
        x_point = column(solution,0);

        return x_point;
    }

    array_1d<double,3> QuadraticErrorFunction::QuadraticErrorFunctionPoint (
        const BoundingBox<Point>& rBox,  
        const std::vector<GeometricalObject*>& rTriangles     
    ) { 
        array_1d<double,3> x_point;
        array_1d<double,3> center =  rBox.GetMinPoint() + (rBox.GetMaxPoint() - rBox.GetMinPoint())/2;
        MatrixType mat_center(3,1);
        column(mat_center,0) = center;
        
        //GeometryArrayType edges = rVoxel.GenerateEdges();

        //Initialize the corresponding matrixes
        BoundedMatrix<double,3,3>ata = ZeroMatrix(3,3);  
        BoundedMatrix<double,3,1>atb = ZeroMatrix(3,1);

        for (std::size_t i = 0; i < rTriangles.size(); i++) {
            array_1d<double,3> normal = CalculateNormal(*rTriangles[i]->pGetGeometry());
            BoundedMatrix<double,3,1> mat_normal;
            column(mat_normal,0) = normal; 

            //We will iterate through the edges using a while loop, so that if a triangles intersects 2 edges (unlikely 
            //but possible), only one will be taken into account to create the matrixes.
            int result = 0; 
            array_1d<double,3> intersection;
            std::size_t j = 0; 
            while(!result && j < 12) { 
                Point end0 = FirstEnd(j, rBox);
                Point end1 = SecondEnd(j, rBox);
                //KRATOS_WATCH(end0);
                //KRATOS_WATCH(end1);
                result = IntersectionUtilities::ComputeTriangleLineIntersection(*rTriangles[i]->pGetGeometry(),end0,end1,intersection);
                j++;
            }
            if (result) {
                //Fill the matrixes with the corresponding information from the intersection and normal
                MatrixType mat_normal_trans = trans(mat_normal);
                ata += prod(mat_normal,mat_normal_trans);
                double aux = MathUtils<double>::Dot(normal,intersection);
                MatrixType mat_aux(1,1,aux);
                atb += prod(mat_normal,mat_aux);
            }
        }

        //Find the eigenvalues and eigenvectors to ata
        MatrixType mat_eigenvectors;
        MatrixType mat_eigenvalues;
        const bool converged = MathUtils<double>::GaussSeidelEigenSystem(ata, mat_eigenvectors, mat_eigenvalues);

        KRATOS_WARNING_IF("QuadraticErrorFunctionPoint", !converged) << "Method for matrix eigenvalues didn't converge" << std::endl;

        MatrixType d(3,3,0);
        for (int i : {0,1,2}) {
            if (mat_eigenvalues(i,i) < 1e-12) {
                d(i,i) = 0; 
            } else {
                d(i,i) = 1.0/mat_eigenvalues(i,i) > 1.0e-12 ? 1.0/mat_eigenvalues(i,i) : 0.0;
            }
        }

        MatrixType ata_inverse;  
        MathUtils<double>::BDBtProductOperation(ata_inverse, d, mat_eigenvectors);
       
        MatrixType ata_c = prod(ata,mat_center);
        MatrixType solution = prod(ata_inverse, atb - ata_c) + mat_center;        
        x_point = column(solution,0);
        return x_point;
    }

    array_1d<double,3> QuadraticErrorFunction::CalculateNormal(const GeometryType& rTriangle) {
        PointsArrayType nodes = rTriangle.Points();
        array_1d<double,3> u = nodes[1] - nodes[0];
        array_1d<double,3> v = nodes[2] - nodes[0];
        array_1d<double,3> normal;
        MathUtils<double>::CrossProduct<array_1d<double,3>,array_1d<double,3>,array_1d<double,3>>(normal,u,v);
        return normal;
    }

    Point QuadraticErrorFunction::FirstEnd(int i, const BoundingBox<Point>& rBox) {
        double xmin = rBox.GetMinPoint()[0];
        double ymin = rBox.GetMinPoint()[1];
        double zmin = rBox.GetMinPoint()[2];
        double xmax = rBox.GetMaxPoint()[0];
        double ymax = rBox.GetMaxPoint()[1];
        double zmax = rBox.GetMaxPoint()[2];

        Point res;
        if (i == 1 || i == 2 || i == 5 || i == 6 || i == 9 || i == 10) res.X() = xmax;
        else res.X() = xmin; 
        if (i == 2 || i == 6 || i == 10 || i == 11) res.Y() = ymax;
        else res.Y() = ymin;
        if (i == 4 || i == 5 || i == 6 || i == 7) res.Z() = zmax;
        else res.Z() = zmin;

        return res;
    }

    Point QuadraticErrorFunction::SecondEnd(int i, const BoundingBox<Point>& rBox) {
        double xmin = rBox.GetMinPoint()[0];
        double ymin = rBox.GetMinPoint()[1];
        double zmin = rBox.GetMinPoint()[2];
        double xmax = rBox.GetMaxPoint()[0];
        double ymax = rBox.GetMaxPoint()[1];
        double zmax = rBox.GetMaxPoint()[2];

        Point res;
        if (i == 0 || i == 1 || i == 4 || i == 5 || i == 9 || i == 10) res.X() = xmax;
        else res.X() = xmin; 
        if (i == 0 || i == 4 || i == 8 || i == 9 ) res.Y() = ymin;
        else res.Y() = ymax;
        if (i == 0 || i == 1 || i == 2 || i == 3) res.Z() = zmin;
        else res.Z() = zmax;

        return res;
    }
}