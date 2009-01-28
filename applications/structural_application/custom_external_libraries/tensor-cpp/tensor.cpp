/**
 * Implementation of the functions from class Tensor.
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


#include "tensor.h"




/* ---------- Member Functions --------- */
using namespace std;

/**
 * Constructor.
 *
 * It reads the rank of the tensor and the dimension of the vectorial
 * space in which it is defined.
 *
 * For example, to create a rank 3 tensor in R^4:
 *    Tensor<double> t(3, 4);
 */
template <class T>
Tensor<T>::Tensor(int rank, long dimension)
{
    m_rank = rank;
    m_dimension = dimension;

    long n = quick_pow(m_dimension, m_rank);  // number of elements

    m_value = new T[n];

    if (m_value == NULL) {
	cerr << "ERROR: not enough space to allocate tensor.\n";
    }
}


/**
 * Copy constructor.
 */
template <class T>
Tensor<T>::Tensor(const Tensor<T>& t)
{
    m_rank      = t.m_rank;
    m_dimension = t.m_dimension;

    long n = quick_pow(m_dimension, m_rank);  // number of elements

    m_value = new T[n];
    memcpy(m_value, t.m_value, n*sizeof(T));
}


/**
 * Destructor.
 */
template <class T>
Tensor<T>::~Tensor()
{
    delete[] m_value;
}


/**
 * Returns the rank of the tensor (number of indices that the tensor has).
 */
template <class T>
int Tensor<T>::rank() const
{
    return m_rank;
}


/**
 * Dimension of the vectorial space in which the tensor is defined.
 */
template <class T>
long Tensor<T>::dimension() const
{
    return m_dimension;
}


/**
 * Get an element, either to set its value or to read it.
 *
 * It uses "C notation" (instead of Fortran notation)
 * that is, numbering from 0 to dimension-1.
 *
 * For example:
 *   Tensor<double> t(3, 8);
 *   t(2,1,6) = 74.346;
 *   cout << t(2,1,6) << endl;
 */
template <class T>
T& Tensor<T>::operator () (long i0, ...)
{
    va_list argptr;
    va_start(argptr, i0);

    long shift = quick_pow(m_dimension, m_rank-1);

    long position = i0 * shift;

    for (int i = 0; i < rank()-1; i++) {
	shift /= m_dimension;
	position += shift * va_arg(argptr, int);
    }

    return m_value[position];
}

/**
 * Same thing as above, for a rank 0 tensor (scalar).
 */
template <class T>
T& Tensor<T>::operator () ()
{
    return m_value[0];
}



/* -------------- */

/**
 * Tensor to a stream, for output purposes.
 *
 * A 3rd rank tensor of dimension 2 would look like:
 *   1  2
 *   3  4
 *
 *   5  6
 *   7  8
 */
template <class T>
ostream& operator << (ostream& s, const Tensor<T>& t)
{
    long n = quick_pow(t.dimension(), t.rank());  // number of elements

    for (long pos = 0; pos < n; pos++) {
	if (pos != 0 && pos % t.dimension() == 0) {
	    long p = pos;
	    do {                     // put number of spaces according to
		s << endl;           // the index we are jumping through
		p /= t.dimension();
	    } while (p % t.dimension() == 0);
	}

	s << t.m_value[pos] << "\t";
    }

    return s;
}



/* ------ Tensor Machine ------ */


/**
 * Addition of two tensors.
 *
 * This returns C_ij = A_ij + B_ij
 */
template <class T>
Tensor<T> operator + (const Tensor<T>& t1, const Tensor<T>& t2)
{
    if (t1.rank() != t2.rank()) {
	cerr << "ERROR: trying to add tensors of different rank\n";
	return Tensor<T>(1, 0);
    }
    if (t1.dimension() != t2.dimension()) {
	cerr << "ERROR: trying to add tensors of different dimension\n";
	return Tensor<T>(1, 0);
    }

    int rank       = t1.rank();
    long dimension = t1.dimension();
    
    Tensor<T> s(rank, dimension);

    long n = quick_pow(dimension, rank);  // number of elements

    for (long pos = 0; pos < n; pos++)
	s.m_value[pos] = t1.m_value[pos] + t2.m_value[pos];

    return s;
}


/**
 * Substraction of two tensors.
 *
 * This returns C_ij = A_ij - B_ij
 */
template <class T>
Tensor<T> operator - (const Tensor<T>& t1, const Tensor<T>& t2)
{
    if (t1.rank() != t2.rank()) {
	cerr << "ERROR: trying to substract tensors of different rank\n";
	return Tensor<T>(1, 0);
    }
    if (t1.dimension() != t2.dimension()) {
	cerr << "ERROR: trying to substract tensors of different dimension\n";
	return Tensor<T>(1, 0);
    }

    int rank       = t1.rank();
    long dimension = t1.dimension();
    
    Tensor<T> s(rank, dimension);

    long n = quick_pow(dimension, rank);  // number of elements

    for (long pos = 0; pos < n; pos++)
	s.m_value[pos] = t1.m_value[pos] - t2.m_value[pos];

    return s;
}


/**
 * Negative of a tensor.
 *
 * This returns C_ij = - A_ij
 */
template <class T>
Tensor<T> operator - (const Tensor<T>& t)
{
    Tensor<T> mt(t.rank(), t.dimension());

    long n = quick_pow(mt.dimension(), mt.rank());  // number of elements

    for (long pos = 0; pos < n; pos++)
	mt.m_value[pos] = - t.m_value[pos];

    return mt;
}


/**
 * Tensorial product of two tensors.
 *
 * This returns
 *   C_i..jk..l = A_i..j * B_k..l
 */
template <class T>
Tensor<T> operator * (const Tensor<T>& t1, const Tensor<T>& t2)
{
    if (t1.dimension() != t2.dimension()) {
	cerr << "ERROR: trying to multiply tensors of different dimension\n";
	return Tensor<T>(1, 0);
    }

    long dimension = t1.dimension();

    Tensor<T> prod(t1.rank() + t2.rank(), dimension);

    // number of elements of input tensors
    long n1 = quick_pow(dimension, t1.rank());
    long n2 = quick_pow(dimension, t2.rank());

    for (long i = 0; i < n1; i++)
	for (long j = 0; j < n2; j++) {
	    long pos = n2*i + j;

	    prod.m_value[pos] = t1.m_value[i] * t2.m_value[j];
	}

    return prod;
}


/**
 * Contracts the indices i and j of a tensor.
 *
 * This function, unlike the others above, is not trivial. The
 * algorithm is based on the correspondence of the indices in the
 * tensor and the serialized correspondig position.
 */
template <class T>
Tensor<T> contract(const Tensor<T>& t_ij, int i, int j)
{
    if (i >= t_ij.rank() || j >= t_ij.rank() || i == j) {
	cerr << "ERROR: Trying to contract tensor of rank "
	     << t_ij.rank() << " with bad indices: "
	     << i << ", " << j << endl;

	return t_ij;
    }

    long dimension = t_ij.dimension();  // to write less later

    // Create a new tensor with the appropiate rank
    Tensor<T> t_ii(t_ij.rank() - 2, dimension);
    

    // Number of elements of the new tensor
    long n = quick_pow(dimension, t_ii.rank());


    // Start counting the indices in the opposite direction,
    // for clarity in the algorithm (but potentially confusing!)
    // Now the last two indices, for instance, will be 1,0.
    i = (t_ij.rank() - 1) - i;
    j = (t_ij.rank() - 1) - j;

    if (i > j) {         // swap indices if necessary, so j > i
	int k = i;
	i = j;
	j = k;
    }
    
    long step = quick_pow(dimension, i) + quick_pow(dimension, j);

    for (long pos = 0; pos < n; pos++) {
	vector<long> pos_in_base = decimal2base(t_ii.rank()+1, dimension, pos);

	vec_insert(pos_in_base, i, 0);	
	vec_insert(pos_in_base, j, 0);

	long init = base2decimal(dimension, pos_in_base);

	T sum = 0;

	for (long l = 0; l < dimension; l++) {
	    long old_pos = init + step * l;
	    sum += t_ij.m_value[old_pos];
	}

	t_ii.m_value[pos] = sum;
    }

    return t_ii;
}


/**
 * Swaps the indices i and j of a tensor.
 *
 * T'_..i..j.. = T_..j..i..
 */
template <class T>
Tensor<T> swap(const Tensor<T>& t_ij, int i, int j)
{
    if (i >= t_ij.rank() || j >= t_ij.rank() || i == j) {
	cerr << "ERROR: Trying to swap tensor of rank "
	     << t_ij.rank() << " with bad indices: "
	     << i << ", " << j << endl;

	return t_ij;
    }

    long dimension = t_ij.dimension();  // to write less later
    int rank = t_ij.rank();

    // Create a new tensor with the appropiate rank
    Tensor<T> t_ji(rank, dimension);
    

    // Number of elements of the new tensor
    long n = quick_pow(dimension, rank);

    // Start counting the indices in the opposite direction,
    // for clarity in the algorithm (but potentially confusing!)
    // Now the last two indices, for instance, will be 1,0.
    i = (rank - 1) - i;
    j = (rank - 1) - j;
    
    for (long pos = 0; pos < n; pos++) {
	vector<long> pos_in_base = decimal2base(rank+1, dimension, pos);
	
	vec_swap(pos_in_base, i, j);

	long newpos = base2decimal(dimension, pos_in_base);

	t_ji.m_value[newpos] = t_ij.m_value[pos];
    }

    return t_ji;
}   



/* ------ General Utilities ------ */

/**
 * Changes base of number "n" to "base".
 *
 * It returns a vector "digits" with digits[i] being the ith digit
 * of number n in base "base".
 */
vector<long> decimal2base(int n_digits, long base, long n)
{
    vector<long> digits;

    for (int i = 0; i < n_digits; i++) {
	digits.push_back(n - base * (n / base));
	n /= base;
    }

    return digits;
}


/**
 * Changes number in base "base" to a long decimal.
 *
 * This is the complementary of function decimal2base().
 */
long base2decimal(long base, const vector<long>& digits)
{
    long n = 0;
    
    long scale = 1;
    for (int i = 0; i < (int) digits.size(); i++) {
	n += digits[i] * scale;
	scale *= base;
    }
    
    return n;
}


/**
 * Inserts an element (value) into a vector (v) at a given position.
 */
void vec_insert(vector<long>& v, int position, long value)
{
    if (position >= (int) v.size()) {
	cerr << "ERROR: trying to insert at position " << position
	     << " on a vector of size " << v.size() << endl;

	return;
    }

    vector<long>::iterator current_pos = v.begin();

    for (int i = 0; i < position && current_pos < v.end(); i++)
	current_pos++;

    v.insert(current_pos, value);
}


/**
 * Swaps indices i and j of vector.
 */
void vec_swap(vector<long>& v, int i, int j)
{
    long temp = v[i];
    v[i] = v[j];
    v[j] = temp;
}



/* -------------- */


/*
 * Explicit instantiation for the case of doubles.
 */
template class Tensor<double>;

template Tensor<double> operator + (const Tensor<double>& t1,
				    const Tensor<double>& t2);
template Tensor<double> operator - (const Tensor<double>& t1,
				    const Tensor<double>& t2);

template Tensor<double> operator - (const Tensor<double>& t);

template Tensor<double> operator * (const Tensor<double>& t1,
				    const Tensor<double>& t2);

template Tensor<double> contract(const Tensor<double>& t_ij, int i, int j);
template Tensor<double> swap(const Tensor<double>& t_ij, int i, int j);

template ostream& operator << (ostream& s, const Tensor<double>& t);


/*
 * And for the same price, explicit instantiation for complex<double>.
 */
#include <complex>

typedef complex<double> dcomplex;

template class Tensor<dcomplex>;

template Tensor<dcomplex> operator + (const Tensor<dcomplex>& t1,
				      const Tensor<dcomplex>& t2);
template Tensor<dcomplex> operator - (const Tensor<dcomplex>& t1,
				      const Tensor<dcomplex>& t2);

template Tensor<dcomplex> operator - (const Tensor<dcomplex>& t);

template Tensor<dcomplex> operator * (const Tensor<dcomplex>& t1,
				      const Tensor<dcomplex>& t2);

template Tensor<dcomplex> contract(const Tensor<dcomplex>& t_ij, int i, int j);
template Tensor<dcomplex> swap(const Tensor<dcomplex>& t_ij, int i, int j);

template ostream& operator << (ostream& s, const Tensor<dcomplex>& t);
