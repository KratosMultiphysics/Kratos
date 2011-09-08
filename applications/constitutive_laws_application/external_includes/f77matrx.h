#include <assert.h>

/*
  class FMATRIX
  =============
  A minimal class used when passing multi-dimensional array
  arguments from C++ to FORTRAN 77 (received as FORTRAN arrays),
  and subsequently returned back to C++ as properly aranged
  C++ arrays.

  Problem : FORTRAN organises data in a "column-first" order,
            while C++ organises data in a "row-first" order.

  Solution:
       (1)  The FMATRIX class can take a C++ array as a constructor
        parameter. A FORTRAN compatible copy of the array is
        then made. The destructor will then copy the result back
        to the original c++ array.

     (2)  The FMATRIX class provides "subscript operators" allowing
        the programmer to read and write from the array, using
        FORTRAN-like syntax and indexing semantics.

  Author: Carsten A. Arnholm, 04-MAR-1996
*/

template <class T>
class FMATRIX {
public:
   FMATRIX(size_t dim1, size_t dim2);
   FMATRIX(T* cpparr, size_t dim1, size_t dim2);
   operator T*();
   T& operator()(size_t index1, size_t index2);
  ~FMATRIX();
public:
   const size_t ndim;  // number of array dimensions
   size_t dim[7];      // size of each dimension
   T*  cpprep;         // original c++ array
   T*  f77rep;         // array used by FORTRAN
};

template <class T>
FMATRIX<T>::FMATRIX(size_t dim1, size_t dim2)
: cpprep(NULL),f77rep(new T[dim1*dim2]),ndim(2)
{
   dim[0]=dim1;
   dim[1]=dim2;
   dim[2]=0;
   dim[3]=0;
   dim[4]=0;
   dim[5]=0;
   dim[6]=0;
}

template <class T>
FMATRIX<T>::FMATRIX(T* cpparr, size_t dim1, size_t dim2)
: cpprep(cpparr),f77rep(new T[dim1*dim2]),ndim(2)
{
   dim[0]=dim1;
   dim[1]=dim2;
   dim[2]=0;
   dim[3]=0;
   dim[4]=0;
   dim[5]=0;
   dim[6]=0;

   // make a FORTRAN-compatible copy of the array
   size_t index_cpp=0;
   size_t index_f77;
   for(size_t i=0;i<dim[0];i++) {
      for(size_t j=0;j<dim[1];j++) {
       index_f77 = j*dim[0] + i;
       f77rep[index_f77] = cpprep[index_cpp++];
    }
   }
}

template <class T>
FMATRIX<T>::operator T*()
{
   // Pass the FORTRAN representation when calling a function
   return f77rep;
}

template <class T>
T&  FMATRIX<T>::operator()(size_t index1, size_t index2)
{
   assert(ndim==2);  // only 2d arrays supported (so far)

   // indexing according to F77 conventions
   size_t index_f77 = index2*dim[0] + index1;

   // return a reference to the array element
   return *(f77rep+index_f77);
}

template <class T>
FMATRIX<T>::~FMATRIX()
{
   if(cpprep) {
      assert(ndim==2);  // only 2d arrays supported (so far)

      // copy back from FORTRAN to C++ array
      size_t index_cpp;
      size_t index_f77=0;
      for(size_t j=0;j<dim[1];j++) {
         for(size_t i=0;i<dim[0];i++) {
            index_cpp = i*dim[1] + j;
            cpprep[index_cpp] = f77rep[index_f77++];
         }
      }
   }

   // delete the FORTRAN copy of the arry
   delete[] f77rep;
}
