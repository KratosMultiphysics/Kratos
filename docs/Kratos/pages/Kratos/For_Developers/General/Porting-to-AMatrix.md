---
title: Porting to AMatrix
keywords: 
tags: [Porting-to-AMatrix.md]
sidebar: kratos_for_developers
summary: 
---

This page provides a list of common changes for migration from ublas to AMatrix

## Compiling

In  order to compile AMatrix with Kratos add the following to your `configure.sh`:

```sh
-DAMATRIX_DIR="/your_kratos_location/Kratos/external_libraries/a_matrix/" \
```

## Using Kratos vector and matrix aliases
It is recommended to use Kratos vector and matrix aliases instead of ones define by ublas:

* Using `Vector` for a dynamic size vector of double
* Using `Matrix` for a dynamic size matrix of double
* Replacing the `boost::numeric::ublas::vector` with `DenseVector`
* Replacing the `boost::numeric::ublas::matrix` with `DenseMatrix`
* Replacing the `boost::numeric::ublas::bounded_vector` with `BoundedVector`
* Replacing the `boost::numeric::ublas::bounded_matrix` with `BoundedMatrix`
* Replacing the `bounded_vector` with `BoundedVector`
* Replacing the `bounded_matrix` with `BoundedMatrix`
* Relplacing the `matrix<` with `DenseMatrix<` **NOTE:** the `<` character should be included in this step to avoid unnecessary changes
* Replacing the `boost::numeric::ublas::prod` with `prod`
* Replacing the `boost::numeric::ublas::trans` with `trans`

## Other common changes
There are differences in interface and also behavior of these to libraries. The most important ones are the resize which does not preserve the values and there is no operator() for vector. Here is the list of common changes:

* Passing false to vector resize: `vector.resize(new_sizse)` to `vector.resize(new_size, false)`
* Passing false to matrix resize: `matrix.resize(size1,size2)` to `matrix.resize(size1,size2, false)`
* There is no constructor with value for matrix and vector. Please use `ScalarMatrix` for that: `Matrix m(size1,size2, value)` to `Matrix m = ScalarMatrix(size1,size2, value)`
* `IdentityMatrix` is only square and accept one size: `IdentityMatrix(size,size)` to `IdentityMatrix(size)`
* The `array_1d` does not have constructor with size: `array_1d<double,3> a(3)` to `array_1d<double,3> a`
* Changing the ublas `vector<double>`  to `Vector`
* Changing the ublas `vector<other_type>` to `DenseVector<other_type>`
* Amatrix classes does not have operator `()`; use `[]` instead: `v(i)` to `v[i]`
* `data()` method returns a raw pointer to the first member and not an iterator
* when dealing with spaces please take care of using correctly the `SparseSpace::DenseVector` and `LocalSpace::DenseVector` as they are different types now.
* There is no `project` function in amatrix. Instead there is a SubMatrix expression:  `AMatrix::SubMatrix<MatrixType>a_sub_matrix(original_matrix, sub_index1, sub_index2, sub_size1, sub_size2);` The arguments are different so it might be used with `KRATOS_USE_AMATRIX` macro in order to maintain the compatibility to ublas
*  Replacing the `matrix_row(rVariables.Strain.Eigen.Vectors,i)` with `row((rVariables.Strain.Eigen.Vectors,i)` 

In order to maintain the backward compatibility until all applications are migrated a A new `KRATOS_USE_AMATRIX` macro is added which is used when the AMatrix interface is different than Ublas. The idea is to remove `KRATOS_USE_AMATRIX` macro after the migration, so please use it wisely and as less as possible.

## Tips and tricks

### Using Matrix and Vector from within an Internals namespace

Sometimes you may want to use `Matrix` and `Vector` types within an `Internals` namespace for your class, but the names `Kratos::Internals::Matrix` and `Kratos::Internals::Vectors` are already taken when you compile Kratos with AMatrix. The name clash may manifest with an error like:

```
error: invalid use of template-name ‘Kratos::Internals::Matrix’ without an argument list
     typedef Matrix MatrixType;
             ^~~~~~
/home/jcotela/Kratos/kratos/utilities/helper_classes_for_constraint_builder.h:431:13: note: class template argument deduction is only available with -std=c++17 or -std=gnu++17
In file included from /home/jcotela/Kratos/kratos/includes/ublas_interface.h:35,
                 from /home/jcotela/Kratos/kratos/includes/serializer.h:31,
                 from /home/jcotela/Kratos/kratos/includes/model_part.h:33,
                 from /home/jcotela/Kratos/kratos/python/add_strategies_to_python.cpp:25:
/home/jcotela/Kratos/kratos/includes/amatrix_interface.h:39:7: note: ‘template<class TDataType, long unsigned int TSize1, long unsigned int TSize2> class Kratos::Internals::Matrix’ declared here
 class Matrix : public AMatrix::MatrixExpression<Matrix<TDataType, TSize1, TSize2>,
       ^~~~~~
```

If this happens, declare your local MatrixType using the complete type name:
```cpp
namespace Kratos
{
    namespace Internals
    {
        //...
        typedef Kratos::Matrix MatrixType;
        typedef Kratos::Vector VectorType;
        //...
    }
}
```

### Runtime errors when importing application strategies, builder and solvers or processes

Keep in mind that the SparseSpace type used as template argument to define strategies and builder and solvers still uses ublas types. This means that you will probably need to update the `add_whatever_to_python.cpp` files in your application.

If you have a type definition that looks like
```cpp
typedef UblasSpace<double, CompressedMatrix, Vector > SparseSpaceType;
```
keep in mind that the `Vector` type has been redefined in AMatrix, so you need to point back to the ublas type:
```cpp
typedef UblasSpace<double, CompressedMatrix, boost::numeric::ublas::vector<double>> SparseSpaceType;
```
Note that this may require adding an include for `"boost/numeric/ublas/matrix.hpp"`.