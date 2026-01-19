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
//                   Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "utilities/qef_utility.h"
#include "includes/geometrical_object.h"
#include "utilities/intersection_utilities.h"
#include "utilities/math_utils.h"
#include "utilities/svd_utils.h"

namespace Kratos
{

array_1d<double, 3> QuadraticErrorFunction::QuadraticErrorFunctionPoint (
    const GeometryType& rVoxel,
    const GeometryArrayType& rTriangles
    )
{
    // Initialize the corresponding matrices
    AuxiliaryClasses aux;

    // Compute the center of the box
    const array_1d<double, 3> center = rVoxel.Center().Coordinates();
    column(aux.MatCenter, 0) = center;

    // Generate the edges of the voxel
    GeometryArrayType edges = rVoxel.GenerateEdges();

    // Iterate over the triangles
    for (std::size_t i = 0; i < rTriangles.size(); i++) {
        noalias(aux.Normal) = rTriangles[i].Normal(aux.AuxCenter);
        column(aux.MatNormal,0) = aux.Normal;

        // We will iterate through the edges using a while loop, so that if a triangles intersects 2 edges (unlikely
        // but possible), only one will be taken into account to create the matrices.
        // TODO: Intersection should be in every side, not just one
        int result = 0;
        std::size_t j = 0;
        while(!result && j < edges.size()) {
            PointsArrayType ends = edges[j++].Points();
            result = IntersectionUtilities::ComputeTriangleLineIntersection(rTriangles[i], ends[0], ends[1], aux.Intersection);
        }

        // If the triangle intersects the edge, update the matrices
        if (result) {
            UpdateATAATB(aux);
        }
    }

    return ComputeQuadraticErrorFunctionPoint(aux.ATA, aux.ATB, aux.MatCenter);
}

/***********************************************************************************/
/***********************************************************************************/

array_1d<double, 3> QuadraticErrorFunction::QuadraticErrorFunctionPoint (
    const BoundingBox<Point>& rBox,
    const std::vector<GeometricalObject*>& rTriangles
    )
{
    // Initialize the corresponding matrices
    AuxiliaryClasses aux;

    // Compute the center of the box
    const array_1d<double, 3> center = rBox.GetMinPoint() + (rBox.GetMaxPoint() - rBox.GetMinPoint())/2;
    column(aux.MatCenter, 0) = center;

    // Iterate over the triangles
    Point end0, end1;
    for (std::size_t i = 0; i < rTriangles.size(); i++) {
        const auto& r_geometry = rTriangles[i]->GetGeometry();
        noalias(aux.Normal) = r_geometry.Normal(aux.AuxCenter);
        column(aux.MatNormal,0) = aux.Normal;

        // We will iterate through the edges using a while loop, so that if a triangles intersects 2 edges (unlikely
        // but possible), only one will be taken into account to create the matrices.
        // TODO: Intersection should be in every side, not just one
        int result = 0;
        std::size_t j = 0;
        while(!result && j < 12) {
            FirstEnd(end0, j, rBox);
            SecondEnd(end1, j, rBox);
            result = IntersectionUtilities::ComputeTriangleLineIntersection(r_geometry, end0, end1, aux.Intersection);
            ++j;
        }

        // If the triangle intersects the edge, update the matrices
        if (result) {
            UpdateATAATB(aux);
        }
    }

    return ComputeQuadraticErrorFunctionPoint(aux.ATA, aux.ATB, aux.MatCenter);
}

/***********************************************************************************/
/***********************************************************************************/

void QuadraticErrorFunction::UpdateATAATB(AuxiliaryClasses& rAuxiliaryClasses)
{
    // Fill the matrices with the corresponding information from the intersection and normal
    noalias(rAuxiliaryClasses.MatNormalTrans) = trans(rAuxiliaryClasses.MatNormal);
    rAuxiliaryClasses.ATA += prod(rAuxiliaryClasses.MatNormal, rAuxiliaryClasses.MatNormalTrans);
    const double dot_product = MathUtils<double>::Dot(rAuxiliaryClasses.Normal, rAuxiliaryClasses.Intersection);
    rAuxiliaryClasses.MatAux(0, 0) = dot_product;
    rAuxiliaryClasses.ATB += prod(rAuxiliaryClasses.MatNormal, rAuxiliaryClasses.MatAux);
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

array_1d<double, 3> QuadraticErrorFunction::ComputeQuadraticErrorFunctionPoint(
    const BoundedMatrix<double,3,3>& rATA,
    const BoundedMatrix<double,3,1>& rATB,
    const BoundedMatrix<double,3,1>& rMatCenter
    )
{
    // LEGACY: Compute eigenvalues and eigenvectors of the ATA matrix
    // // Compute the eigen-decomposition of the ATA matrix.
    // // mat_eigenvectors will contain the eigenvectors,
    // // and mat_eigenvalues is a diagonal matrix with the corresponding eigenvalues.
    // BoundedMatrix<double,3,3> mat_eigenvectors, mat_eigenvalues;
    // const bool converged = MathUtils<double>::GaussSeidelEigenSystem(rATA, mat_eigenvectors, mat_eigenvalues);
    // KRATOS_WARNING_IF("QuadraticErrorFunctionPoint", !converged) << "Method for matrix eigenvalues didn't converge" << std::endl;

    // Compute SVD decomposition of the ATA matrix
    // The SVD decomposition is used to compute the pseudo-inverse of the ATA matrix.
    // The pseudo-inverse is used to solve the linear system of equations. (https://www.johndcook.com/blog/2018/05/05/svd/)
    // U_matrix and V_matrix are orthogonal matrices, and S_matrix is a diagonal matrix with the singular values.
    BoundedMatrix<double, 3, 3> U_matrix, S_matrix, V_matrix;
    SVDUtils<double>::JacobiSingularValueDecomposition(rATA, U_matrix, S_matrix, V_matrix);

    // Define a tolerance to avoid division by very small or huge numbers.
    const double tolerance = 1e-12;

    // Construct the diagonal matrix 'd' for computing the pseudo-inverse.
    // For each singularvalue, if it is too small (or its inverse is too small), set its inverse to zero.
    BoundedMatrix<double,3,3> d = ZeroMatrix(3, 3);
    for (std::size_t i = 0; i < 3; ++i) {
        const double singularvalue = S_matrix(i, i); // mat_eigenvalues(i, i); NOTE: LEGACY
        const double inv_singularvalue = 1.0 / singularvalue;
        // If the singularvalue is below tolerance or its inverse is negligible, set to zero.
        d(i,i) = (singularvalue < tolerance || inv_singularvalue < tolerance) ? 0.0 : inv_singularvalue;
    }

    // Compute the pseudo-inverse of the ATA matrix using the SVD-decomposition.
    BoundedMatrix<double,3,3> ata_pseudo_inverse;
    // MathUtils<double>::BDBtProductOperation(ata_pseudo_inverse, d, mat_eigenvectors); NOTE: LEGACY
    noalias(ata_pseudo_inverse) = prod(trans(V_matrix), BoundedMatrix<double,3,3>(prod(d, trans(U_matrix))));

    // Compute the product of ATA and the center vector.
    const BoundedMatrix<double, 3, 3> ata_c = prod(rATA, rMatCenter);
    // Compute the solution vector: center + ATA_inverse * (ATB - ATA * center)
    const BoundedMatrix<double, 3, 1> solution = prod(ata_pseudo_inverse, rATB - ata_c) + rMatCenter;

    // Extract the resulting point from the solution vector.
    return column(solution, 0);
}

}