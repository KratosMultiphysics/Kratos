//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//
//

// System includes

// External includes

// Project includes
#include "utilities/math_utils.h"

// Include base h
#include "geometry_metric_calculator.h"

namespace Kratos
{

template<>
void KRATOS_API(KRATOS_CORE) GeometryMetricCalculator::CalculateMetricTensorData<2,3>(
    const GeometryType& rGeometry,
    BoundedMatrix<double, 2, 2>& rMetricTensor,
    double& rReferenceElementSize,
    double& rMetricInfimum,
    double& rMetricSupremum)
{
    const array_1d<double, 3>& r_p_1 = rGeometry[0].Coordinates();
    const array_1d<double, 3>& r_p_2 = rGeometry[1].Coordinates();
    const array_1d<double, 3>& r_p_3 = rGeometry[2].Coordinates();

    // Solve the metric problem trans(e)*M*e = 1
    // This means, find the coefficients of the matrix M such that all the edges have unit length
    Vector sol; //FIXME: WE SHOULD MAKE THE MATH UTILS COMPATIBLE WITH THIS
    array_1d<double,3> aux_vect;
    BoundedMatrix<double,3,3> aux_mat;
    aux_mat(0,0) = std::pow(r_p_1[0]-r_p_2[0], 2); aux_mat(0,1) = 2.0*(r_p_1[0]-r_p_2[0])*(r_p_1[1]-r_p_2[1]); aux_mat(0,2) = std::pow(r_p_1[1]-r_p_2[1], 2);
    aux_mat(1,0) = std::pow(r_p_1[0]-r_p_3[0], 2); aux_mat(1,1) = 2.0*(r_p_1[0]-r_p_3[0])*(r_p_1[1]-r_p_3[1]); aux_mat(1,2) = std::pow(r_p_1[1]-r_p_3[1], 2);
    aux_mat(2,0) = std::pow(r_p_2[0]-r_p_3[0], 2); aux_mat(2,1) = 2.0*(r_p_2[0]-r_p_3[0])*(r_p_2[1]-r_p_3[1]); aux_mat(2,2) = std::pow(r_p_2[1]-r_p_3[1], 2);
    aux_vect[0] = 1.0;
    aux_vect[1] = 1.0;
    aux_vect[2] = 1.0;        
    MathUtils<double>::Solve(aux_mat, sol, aux_vect);

    // Set the metric tensor
    rMetricTensor(0,0) = sol[0]; rMetricTensor(0,1) = sol[1];
    rMetricTensor(1,0) = sol[1]; rMetricTensor(1,1) = sol[2];

    // Calculate the eigenvalues of the metric tensor to obtain the ellipsis of inertia axes lengths
    BoundedMatrix<double,2,2> eigenvects, eigenvals;
    MathUtils<double>::GaussSeidelEigenSystem(rMetricTensor, eigenvects, eigenvals);

    // Calculate the reference element size as the average of the ellipsis of intertia axes lengths
    const double h_1 = std::sqrt(1.0 / eigenvals(0,0));
    const double h_2 = std::sqrt(1.0 / eigenvals(1,1));
    rReferenceElementSize = 0.5 * (h_1 + h_2);

    // Get te infimum and supremum of the metric
    rMetricInfimum = std::min(eigenvals(0,0), eigenvals(1,1));
    rMetricSupremum = std::max(eigenvals(0,0), eigenvals(1,1));
    
    // Make the metric dimensionless
    rMetricTensor *= std::pow(rReferenceElementSize,2);
}

template<>
void KRATOS_API(KRATOS_CORE) GeometryMetricCalculator::CalculateMetricTensorData<3,4>(
    const GeometryType& rGeometry,
    BoundedMatrix<double, 3, 3>& rMetricTensor,
    double& rReferenceElementSize,
    double& rMetricInfimum,
    double& rMetricSupremum)
{
        // Solve the metric problem trans(e)*M*e = 1
        // This means, find the coefficients of the matrix M such that all the edges have unit length
        //TODO: CHECK THIS, I THINK IT'S WRONG, WE NEED TO FOLLOW THE KRATOS EDGES CONVENTION...
        Vector sol;
        array_1d<double, 6> aux_vect;
        BoundedMatrix<double, 6, 6> aux_mat;
        unsigned int row = 0;
        for (unsigned int i = 0; i < 3; ++i) {
            const auto& i_coord = rGeometry[i].Coordinates();
            for (unsigned int j = i + 1; j < 4; ++j) {
                const auto& j_coord = rGeometry[j].Coordinates();
                aux_mat(row, 0) = std::pow(i_coord[0]-j_coord[0], 2);
                aux_mat(row, 1) = 2.0*(i_coord[0]-j_coord[0])*(i_coord[1]-j_coord[1]);
                aux_mat(row, 2) = 2.0*(i_coord[0]-j_coord[0])*(i_coord[2]-j_coord[2]);
                aux_mat(row, 3) = std::pow(i_coord[1]-j_coord[1], 2);
                aux_mat(row, 4) = 2.0*(i_coord[1]-j_coord[1])*(i_coord[2]-j_coord[2]);
                aux_mat(row, 5) = std::pow(i_coord[2]-j_coord[2], 2);
                aux_vect(row) = 1.0;
                row++;
            }
        }
        MathUtils<double>::Solve(aux_mat, sol, aux_vect);

        // Set the metric tensor
        rMetricTensor(0,0) = sol[0]; rMetricTensor(0,1) = sol[1]; rMetricTensor(0,2) = sol[2];
        rMetricTensor(1,0) = sol[1]; rMetricTensor(1,1) = sol[3]; rMetricTensor(1,2) = sol[4];
        rMetricTensor(2,0) = sol[2]; rMetricTensor(2,1) = sol[4]; rMetricTensor(2,2) = sol[5];

        // Calculate the eigenvalues of the metric tensor to obtain the ellipsis of inertia axes lengths
        BoundedMatrix<double,3,3> eigenvects, eigenvals;
        MathUtils<double>::GaussSeidelEigenSystem(rMetricTensor, eigenvects, eigenvals);

        // Calculate the reference element size as the average of the ellipsis of intertia axes lengths
        const double h_1 = std::sqrt(1.0 / eigenvals(0,0));
        const double h_2 = std::sqrt(1.0 / eigenvals(1,1));
        const double h_3 = std::sqrt(1.0 / eigenvals(2,2));
        rReferenceElementSize = (h_1 + h_2 + h_3) / 3.0;

        // Get te infimum and supremum of the metric
        rMetricInfimum = std::min(eigenvals(0,0), std::min(eigenvals(1,1), eigenvals(2,2)));
        rMetricSupremum = std::max(eigenvals(0,0), std::max(eigenvals(1,1), eigenvals(2,2)));
}

} // namespace Kratos.
