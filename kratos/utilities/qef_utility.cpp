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
        array_1d<double,3> x_point;
        array_1d<double,3> center = CalculateCenter(rVoxel);
        VectorType v_center = center; 
        MatrixType m_center(3,1);
        column(m_center,0) = center;
        GeometryArrayType edges = rVoxel.GenerateEdges();

        //Initialize the corresponding matrixes
        MatrixType ata(3,3,0);  //3x3 matrix initialized to 0
        MatrixType atb(3,1,0);  //3x1 matrix
        MatrixType btb(1,1,0);  //1x1 matrix

        for (std::size_t i = 0; i < rTriangles.size(); i++) {
            array_1d<double,3> normal = CalculateNormal(rTriangles[i]);
            VectorType v_normal = normal; 
            MatrixType m_normal(3,1);
            column(m_normal,0) = normal;

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
                MatrixType m_normal_trans = trans(m_normal);
                MatrixType help = prod(m_normal,m_normal_trans);
                ata = ata + help;
                double aux = MathUtils<double>::Dot(normal,intersection);
                MatrixType m_aux(1,1,aux);
                atb = atb + prod(m_normal,m_aux);
                btb = btb + prod(m_aux,m_aux); 
            }
        }

        //Find the eigenvalues and eigenvectors to AtA
        MatrixType m_eigenvectors;
        MatrixType m_eigenvalues;
        const bool converged = MathUtils<double>::GaussSeidelEigenSystem(ata, m_eigenvectors, m_eigenvalues);

        MatrixType d(3,3,0);
        for (int i : {0,1,2}) {
            m_eigenvalues(i,i) < 1e-12 ? d(i,i) = 0 : d(i,i) = Check(1.0/m_eigenvalues(i,i), 1e-12);
        }

        MatrixType ata_inverse;  
        MathUtils<double>::BDBtProductOperation(ata_inverse, d, m_eigenvectors);
       
        MatrixType ata_c = prod(ata,m_center);
        MatrixType solution = prod(ata_inverse, atb - ata_c) + m_center;        
        x_point = column(solution,0);

        return x_point;
    }


    array_1d<double,3> QEF::CalculateCenter(const GeometryType& rVoxel) {
        PointsArrayType nodes = rVoxel.Points();
        double x, y, z;
        x = y = z = 0;
        for(std::size_t i = 0; i < nodes.size(); i++) {
            x += nodes[i].X();
            y += nodes[i].Y();
            z += nodes[i].Z();
        }
        array_1d<double,3> center({x/nodes.size(),y/nodes.size(),z/nodes.size()});
        return center;
    }

    array_1d<double,3> QEF::CalculateNormal(const GeometryType& rTriangle) {
        PointsArrayType nodes = rTriangle.Points();
        array_1d<double,3> u = nodes[1] - nodes[0];
        array_1d<double,3> v = nodes[2] - nodes[0];
        array_1d<double,3> normal;
        MathUtils<double>::CrossProduct<array_1d<double,3>,array_1d<double,3>,array_1d<double,3>>(normal,u,v);
        return normal;
    }

    double QEF::Check(
        const double d, 
        const double epsilon)
    {
        if (d > epsilon)  {
            return d;
        } else {
            return 0;
        } 
    }
}