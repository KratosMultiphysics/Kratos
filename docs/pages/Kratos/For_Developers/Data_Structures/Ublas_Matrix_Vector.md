---
title: Using alias for Ublas matrix and vector
keywords: 
tags: [Using-alias-for-Ublas-matrix-and-vector.md]
sidebar: kratos_for_developers
summary: 
---

The first step to be taken in order to migrate from ublas dense matrix and vectors to AMatrix consists in use of Kratos alias for different matrix and vector classes of ublas doing the following replaces which are recommended to be done in following order:

1. Replacing the `boost::numeric::ublas::vector` with `DenseVector`
2. Replacing the `boost::numeric::ublas::matrix` with `DenseMatrix`
3. Replacing the `boost::numeric::ublas::bounded_vector` with `BoundedVector`
4. Replacing the `boost::numeric::ublas::bounded_matrix` with `BoundedMatrix`
5. Replacing the `bounded_vector` with `BoundedVector`
6. Replacing the `bounded_matrix` with `BoundedMatrix`
7. Relplacing the `matrix<` with `DenseMatrix<` **NOTE:** the `<` character should be included in this step to avoid unnecessary changes
8. Replacing the `boost::numeric::ublas::prod` with `prod`
9. Replacing the `boost::numeric::ublas::trans` with `trans`


