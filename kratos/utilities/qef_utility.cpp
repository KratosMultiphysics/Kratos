//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ariadna Cortes
//

// System includes

// External includes

// Project includes
#include "utilities/qef_utility.h"
#include "utilities/intersection_utilities.h"
#include "includes/ublas_interface.h"
#include "includes/geometrical_object.h"

namespace Kratos
{

array_1d<double,3> QuadraticErrorFunction::QuadraticErrorFunctionPoint (
    const GeometryType& rVoxel,
    const GeometryArrayType& rTriangles
    )
{
    const array_1d<double,3> center = rVoxel.Center().Coordinates();
    BoundedMatrix<double,3,1> mat_center(3,1);
    column(mat_center,0) = center;
    GeometryArrayType edges = rVoxel.GenerateEdges();

    //Initialize the corresponding matrices
    BoundedMatrix<double,3,3> ata = ZeroMatrix(3,3);
    BoundedMatrix<double,3,1> atb = ZeroMatrix(3,1);

    // Define auxiliary variables
    array_1d<double,3> normal;
    array_1d<double,3> intersection;
    BoundedMatrix<double,3,1> mat_normal;
    BoundedMatrix<double,1,3> mat_normal_trans;
    BoundedMatrix<double,1, 1> mat_aux;
    GeometryType::CoordinatesArrayType aux_center = ZeroVector(3);
    for (std::size_t i = 0; i < rTriangles.size(); i++) {
        noalias(normal) = rTriangles[i].Normal(aux_center);
        column(mat_normal,0) = normal;

        //We will iterate through the edges using a while loop, so that if a triangles intersects 2 edges (unlikely 
        //but possible), only one will be taken into account to create the matrixes.
        int result = 0;
        std::size_t j = 0;
        while(!result && j < edges.size()) {
            PointsArrayType ends = edges[j++].Points();
            result = IntersectionUtilities::ComputeTriangleLineIntersection(rTriangles[i], ends[0], ends[1], intersection);
        }
        if (result) {
            // Fill the matrixes with the corresponding information from the intersection and normal
            noalias(mat_normal_trans) = trans(mat_normal);
            ata += prod(mat_normal, mat_normal_trans);
            const double aux = MathUtils<double>::Dot(normal, intersection);
            mat_aux(0, 0) = aux;
            atb += prod(mat_normal, mat_aux);
        }
    }

    return ComputeQuadraticErrorFunctionPoint(ata, atb, mat_center);
}

/***********************************************************************************/
/***********************************************************************************/

array_1d<double,3> QuadraticErrorFunction::QuadraticErrorFunctionPoint (
    const BoundingBox<Point>& rBox,
    const std::vector<GeometricalObject*>& rTriangles
    )
{
    const array_1d<double,3> center = rBox.GetMinPoint() + (rBox.GetMaxPoint() - rBox.GetMinPoint())/2;
    BoundedMatrix<double,3,1> mat_center;
    column(mat_center,0) = center;

    //Initialize the corresponding matrices
    BoundedMatrix<double,3,3> ata = ZeroMatrix(3,3);
    BoundedMatrix<double,3,1> atb = ZeroMatrix(3,1);

    // Define auxiliary variables
    array_1d<double,3> normal;
    array_1d<double,3> intersection;
    BoundedMatrix<double,3,1> mat_normal;
    BoundedMatrix<double,1,3> mat_normal_trans;
    BoundedMatrix<double,1, 1> mat_aux;
    GeometryType::CoordinatesArrayType aux_center = ZeroVector(3);
    Point end0, end1;
    for (std::size_t i = 0; i < rTriangles.size(); i++) {
        const auto& r_geometry = rTriangles[i]->GetGeometry();
        noalias(normal) = r_geometry.Normal(aux_center);
        column(mat_normal,0) = normal;

        //We will iterate through the edges using a while loop, so that if a triangles intersects 2 edges (unlikely 
        //but possible), only one will be taken into account to create the matrixes.
        int result = 0;
        std::size_t j = 0;
        while(!result && j < 12) {
            FirstEnd(end0, j, rBox);
            SecondEnd(end1, j, rBox);
            result = IntersectionUtilities::ComputeTriangleLineIntersection(r_geometry, end0, end1, intersection);
            j++;
        }
        if (result) {
            // Fill the matrixes with the corresponding information from the intersection and normal
            noalias(mat_normal_trans) = trans(mat_normal);
            ata += prod(mat_normal, mat_normal_trans);
            const double aux = MathUtils<double>::Dot(normal, intersection);
            mat_aux(0, 0) = aux;
            atb += prod(mat_normal, mat_aux);
        }
    }

    return ComputeQuadraticErrorFunctionPoint(ata, atb, mat_center);
}

/***********************************************************************************/
/***********************************************************************************/

void QuadraticErrorFunction::FirstEnd(
    Point& rPoint,
    const int i,
    const BoundingBox<Point>& rBox
    )
{
    if (i == 1 || i == 2 || i == 5 || i == 6 || i == 9 || i == 10) {
        rPoint.X() = rBox.GetMaxPoint()[0];
    } else {
        rPoint.X() = rBox.GetMinPoint()[0];
    }
    if (i == 2 || i == 6 || i == 10 || i == 11) {
        rPoint.Y() = rBox.GetMaxPoint()[1];
    } else {
        rPoint.Y() = rBox.GetMinPoint()[1];
    }
    if (i == 4 || i == 5 || i == 6 || i == 7) {
        rPoint.Z() = rBox.GetMaxPoint()[2];
    } else {
        rPoint.Z() = rBox.GetMinPoint()[2];
    }
}

/***********************************************************************************/
/***********************************************************************************/

void QuadraticErrorFunction::SecondEnd(
    Point& rPoint,
    const int i,
    const BoundingBox<Point>& rBox
    )
{
    if (i == 0 || i == 1 || i == 4 || i == 5 || i == 9 || i == 10) {
        rPoint.X() = rBox.GetMaxPoint()[0];
    } else {
        rPoint.X() = rBox.GetMinPoint()[0];
    }
    if (i == 0 || i == 4 || i == 8 || i == 9 ) {
        rPoint.Y() = rBox.GetMinPoint()[1];
    } else {
        rPoint.Y() = rBox.GetMaxPoint()[1];
    }
    if (i == 0 || i == 1 || i == 2 || i == 3) {
        rPoint.Z() = rBox.GetMinPoint()[2];
    } else {
        rPoint.Z() = rBox.GetMaxPoint()[2];
    }
}

/***********************************************************************************/
/***********************************************************************************/

array_1d<double,3> QuadraticErrorFunction::ComputeQuadraticErrorFunctionPoint(
    const BoundedMatrix<double,3,3>& rATA,
    const BoundedMatrix<double,3,1>& rATB,
    const BoundedMatrix<double,3,1>& rMatCenter
    )
{
    // Compute the eigen-decomposition of the ATA matrix.
    // mat_eigenvectors will contain the eigenvectors,
    // and mat_eigenvalues is a diagonal matrix with the corresponding eigenvalues.
    BoundedMatrix<double,3,3> mat_eigenvectors, mat_eigenvalues;
    const bool converged = MathUtils<double>::GaussSeidelEigenSystem(rATA, mat_eigenvectors, mat_eigenvalues);
    KRATOS_WARNING_IF("QuadraticErrorFunctionPoint", !converged) << "Method for matrix eigenvalues didn't converge" << std::endl;

    // Define a tolerance to avoid division by very small or huge numbers.
    const double tolerance = 1e-12;

    // Construct the diagonal matrix 'd' for computing the pseudo-inverse.
    // For each eigenvalue, if it is too small (or its inverse is too small), set its inverse to zero.
    BoundedMatrix<double,3,3> d = ZeroMatrix(3,3);
    for (std::size_t i = 0; i < 3; ++i) {
        const double eigenvalue = mat_eigenvalues(i,i);
        const double inv_eigenvalue = 1.0 / eigenvalue;
        // If the eigenvalue is below tolerance or its inverse is negligible, set to zero.
        d(i,i) = (eigenvalue < tolerance || inv_eigenvalue < tolerance) ? 0.0 : inv_eigenvalue;
    }

    // Compute the pseudo-inverse of the ATA matrix using the eigen-decomposition.
    BoundedMatrix<double,3,3> ata_inverse;
    MathUtils<double>::BDBtProductOperation(ata_inverse, d, mat_eigenvectors);

    // Compute the product of ATA and the center vector.
    const BoundedMatrix<double,3,3> ata_c = prod(rATA, rMatCenter);
    // Compute the solution vector: center + ATA_inverse * (ATB - ATA * center)
    const BoundedMatrix<double,3,1> solution = prod(ata_inverse, rATB - ata_c) + rMatCenter;

    // Extract the resulting point from the solution vector.
    return column(solution, 0);
}

}