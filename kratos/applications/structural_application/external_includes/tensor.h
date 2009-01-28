/**
 * Definition of a tensor.
 *
 * JBC. March 2003.
 */

/* --- Begin of licensing stuff ---
 * Tensor - a library to manage tensors in C++.
 * Copyright (C) 2003  Jordi Burguet Castell
 *
 * This program is free software; you can redistribute it
 * and/or modify it under the terms of the GNU General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later
 * version.
 *
 * The explicit terms of the license can be found at
 * http://www.gnu.org/licenses/gpl.html
 *  --- end of licensing stuff ---
 */


#ifndef TENSOR_H
#define TENSOR_H

#include <iostream>
#include <vector>

#include <cstdarg>


// Forward declarations
template <class T> class Tensor;

template <class T> std::ostream& operator << (std::ostream& s,
					      const Tensor<T>& t);

template <class T> Tensor<T> operator + (const Tensor<T>& t1,
					 const Tensor<T>& t2);
template <class T> Tensor<T> operator - (const Tensor<T>& t1,
					 const Tensor<T>& t2);

template <class T> Tensor<T> operator - (const Tensor<T>& t);

template <class T> Tensor<T> operator * (const Tensor<T>& t1,
					 const Tensor<T>& t2);

template <class T> Tensor<T> contract(const Tensor<T>& t_ij, int i, int j);

template <class T> Tensor<T> swap(const Tensor<T>& t_ij, int i, int j);


/**
 * Class Tensor.
 *
 * Represents a tensor of arbitrary rank and dimension of its
 * associated vectorial space.
 */
template <class T>
class Tensor {
private:
    T* m_value;          //!< array with all the values

    int m_rank;          //!< rank of the tensor
    long m_dimension;    //!< dimension of the vectorial space

public:
    Tensor(int rank, long dimension);
    Tensor(const Tensor<T>& t);
    ~Tensor();
    
    /* Read */
    int rank() const;           //!< rank of the tensor
    long dimension() const;     //!< dimension of the tensor

    /* Read and write */
    T& operator () (long i0, ...);  //!< A(i0, i1, i2, ..., i_rank-1)
    T& operator () ();              //!< A()

    /* Friend functions */
    friend std::ostream& operator << <> (std::ostream& s, const Tensor<T>& t);

    friend Tensor<T> operator + <> (const Tensor<T>& t1, const Tensor<T>& t2);
    friend Tensor<T> operator - <> (const Tensor<T>& t1, const Tensor<T>& t2);

    friend Tensor<T> operator - <> (const Tensor<T>& t);

    friend Tensor<T> operator * <> (const Tensor<T>& t1, const Tensor<T>& t2);
    friend Tensor<T> contract <> (const Tensor<T>& t_ij, int i, int j);
    friend Tensor<T> swap <> (const Tensor<T>& t_ij, int i, int j);
};



/*
 * Utilities.
 */

inline long quick_pow(long d, int r) {
    long p = 1; for (int i = 0; i < r; i++) p *= d; return p;
}


std::vector<long> decimal2base(int n_digits, long base, long n);

long base2decimal(long base, const std::vector<long>& digits);

void vec_insert(std::vector<long>& v, int position, long value);

void vec_swap(std::vector<long>& v, int i, int j);

#endif  // TENSOR_H
