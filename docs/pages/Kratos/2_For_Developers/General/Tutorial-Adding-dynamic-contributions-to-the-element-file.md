---
title: Adding dynamic contributions to the element file
keywords: 
tags: [Tutorial-Adding-dynamic-contributions-to-the-element-file.md]
sidebar: kratos_for_developers
summary: 
---

# Overview

Here we add the corresponding dynamic contributions to `MyLaplacianElement`.

# "Mass" matrix

There is no second term derivative for the temperature, so the corresponding "mass" matrix will be a zero matrix. 

We need to add to `my_laplacian_element.h`:

```cpp
/**
 * @brief This is called during the assembling process in order to calculate the elemental mass matrix
 * @param rMassMatrix The elemental mass matrix
 * @param rCurrentProcessInfo The current process info instance
 */
void CalculateMassMatrix(
    MatrixType& rMassMatrix,
    ProcessInfo& rCurrentProcessInfo 
    ) override;
```

Then we will add the following to our `my_laplacian_element.cpp`, this is a trivial code with the code already computed on the base class `element.h` it will be enough, we are doing this for consistency.

```cpp
void MyLaplacianElement::CalculateMassMatrix(
    MatrixType& rMassMatrix,
    ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    SizeType dimension = r_geom.WorkingSpaceDimension();
    SizeType number_of_nodes = r_geom.size();
    SizeType mat_size = dimension * number_of_nodes;

    rMassMatrix = ZeroMatrix( mat_size, mat_size );

    KRATOS_CATCH("");
}
```

# "Damping" matrix

This "damping matrix" will be computed similar to a regular mass matrix, but using the specific heat instead of the density.

We need to add to `my_laplacian_element.h`:

```cpp
/**
 * @brief This is called during the assembling process in order to calculate the elemental damping matrix
 * @param rDampingMatrix The elemental damping matrix
 * @param rCurrentProcessInfo The current process info instance
 */
void CalculateDampingMatrix(
    MatrixType& rDampingMatrix,
    ProcessInfo& rCurrentProcessInfo 
    ) override;
```

Now we can add the following to `my_laplacian_element.cpp`:

```cpp
void MyLaplacianElement::CalculateDampingMatrix(
    MatrixType& rDampingMatrix,
    ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    const auto& r_geom = GetGeometry();
    const auto& r_prop = GetProperties();
    SizeType dimension = r_geom.WorkingSpaceDimension();
    SizeType number_of_nodes = r_geom.size();
    SizeType mat_size = dimension * number_of_nodes;

    rDampingMatrix = ZeroMatrix( mat_size, mat_size );

    KRATOS_ERROR_IF_NOT(r_prop.Has( SPECIFIC_HEAT )) << "SPECIFIC_HEAT has to be provided for the calculation of the DampingMatrix!" << std::endl;

    const double specific_heat = r_prop[SPECIFIC_HEAT];
    const double thickness = (dimension == 2 && r_prop.Has(THICKNESS)) ? r_prop[THICKNESS] : 1.0;

    const bool compute_lumped_damped_matrix =  rCurrentProcessInfo.Has(COMPUTE_LUMPED_MASS_MATRIX) ? rCurrentProcessInfo[COMPUTE_LUMPED_MASS_MATRIX] : false;

    // LUMPED DAMPING MATRIX
    if (compute_lumped_damped_matrix == true) {
        const double total_specific_heat = GetGeometry().Volume() * specific_heat * thickness;

        Vector lumping_factors;
        lumping_factors = GetGeometry().LumpingFactors( lumping_factors );

        for ( IndexType i = 0; i < number_of_nodes; ++i ) {
            const double temp = lumping_factors[i] * total_specific_heat;
            for ( IndexType j = 0; j < dimension; ++j ) {
                IndexType index = i * dimension + j;
                rDampingMatrix( index, index ) = temp;
            }
        }
    } else { // CONSISTENT DAMPING
        Matrix J0(dimension, dimension);

        IntegrationMethod integration_method = IntegrationUtilities::GetIntegrationMethodForExactMassMatrixEvaluation(r_geom);
        const GeometryType::IntegrationPointsArrayType& integration_points = r_geom.IntegrationPoints( integration_method );
        const Matrix& Ncontainer = r_geom.ShapeFunctionsValues(integration_method);

        for ( IndexType point_number = 0; point_number < integration_points.size(); ++point_number ) {
            GeometryUtils::JacobianOnInitialConfiguration(r_geom, integration_points[point_number], J0);
            const double detJ0 = MathUtils<double>::DetMat(J0);
            const double integration_weight = GetIntegrationWeight(integration_points, point_number, detJ0) * thickness;
            const Vector& rN = row(Ncontainer,point_number);

            for ( IndexType i = 0; i < number_of_nodes; ++i ) {
                const SizeType index_i = i * dimension;

                for ( IndexType j = 0; j < number_of_nodes; ++j ) {
                    const SizeType index_j = j * dimension;
                    const double NiNj_weight = rN[i] * rN[j] * integration_weight * specific_heat;

                    for ( IndexType k = 0; k < dimension; ++k )
                        rDampingMatrix( index_i + k, index_j + k ) += NiNj_weight;
                }
            }
        }
    }

    KRATOS_CATCH("");
}
```
